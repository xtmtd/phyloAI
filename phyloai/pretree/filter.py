"""Marker-level filtering: TAPER, TreeShrink, metric rules, and clustering."""

from __future__ import annotations

import csv
import io
import json
import shlex
import shutil
import subprocess
import sys
import tempfile
import time
import math
import re
from collections.abc import Callable
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Any

from Bio import SeqIO
from rich.table import Table
from rich.console import Console

from phyloai.core.checkpoint import (
    Checkpoint,
    load_checkpoint,
    save_checkpoint_atomic,
    validate_resume_params,
)
from phyloai.core.env import ToolEnv
from phyloai.core.file_matching import (
    logical_msa_locus_name,
    pair_msa_and_tree_maps,
    scan_msa_dir,
    scan_tree_dir,
)
from phyloai.core.runner import Runner
from phyloai.core.sequence_normalization import (
    resolve_seq_type,
)
from phyloai.core.sequence_output_validation import validate_fasta_output
from phyloai.core.schema import write_result_json
from phyloai.pretree.checkpoint_helpers import (
    build_initial_checkpoint,
    mark_task,
    plan_resume,
)

_CHECKPOINT_FLUSH_INTERVAL = 2.0
_console = Console()


# --- Shared output helpers ---

def _write_csv_table(rows: list[dict], path: Path, columns: list[str], delimiter: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=columns, delimiter=delimiter, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def _table_delimiter(table_format: str) -> str:
    return "\t" if table_format == "tsv" else ","


def _table_suffix(table_format: str) -> str:
    return ".tsv" if table_format == "tsv" else ".csv"


def _common_output_conflict(output_dir: Path, overwrite: bool) -> None:
    if output_dir.exists() and any(output_dir.iterdir()):
        if not overwrite:
            raise ValueError(
                f"Output directory '{output_dir}' already exists and is non-empty. "
                "Use --overwrite to replace it."
            )
        shutil.rmtree(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)



def render_filter_summary_table(summary: dict) -> Table:
    table = Table(title="Filter Results")
    table.add_column("Metric", style="cyan")
    table.add_column("Value", style="green")
    for key, value in summary.items():
        table.add_row(str(key), str(value))
    return table


def _compute_retained_msa_stats(msa_paths: list[Path]) -> dict:
    if not msa_paths:
        return {"n_msa": 0, "total_length": 0, "mean_length": 0,
                "min_length": 0, "max_length": 0, "mean_taxa": 0}
    lengths = []
    taxa_counts = []
    for p in msa_paths:
        try:
            records = list(SeqIO.parse(str(p), "fasta"))
            if records:
                lengths.append(len(records[0].seq))
                taxa_counts.append(len(records))
        except Exception:
            continue
    if not lengths:
        return {"n_msa": 0, "total_length": 0, "mean_length": 0,
                "min_length": 0, "max_length": 0, "mean_taxa": 0}
    return {
        "n_msa": len(lengths),
        "total_length": sum(lengths),
        "mean_length": round(sum(lengths) / len(lengths), 2),
        "min_length": min(lengths),
        "max_length": max(lengths),
        "mean_taxa": round(sum(taxa_counts) / len(taxa_counts), 2),
    }


# --- TAPER ---

_TAPER_CUTOFF_DEFAULT = 3
_TAPER_MANAGED_FLAGS = {"-m", "-a", "-c", "-l"}
_TAPER_NT_CMD_EXTRA = ["-m", "N", "-a", "N"]
_STANDARD_AA = set("ARNDCQEGHILKMFPSTWYV")


def _build_taper_cmd(
    input_file: Path, output_file: Path, seq_type: str, cutoff: int,
    julia_exe: str, taper_script: str, tool_args: str | None,
) -> list[str]:
    cmd = [julia_exe, taper_script, "-c", str(cutoff)]
    if seq_type == "NT":
        cmd.extend(_TAPER_NT_CMD_EXTRA)
    if tool_args:
        extra = shlex.split(tool_args)
        for flag in _TAPER_MANAGED_FLAGS:
            if flag in extra:
                raise ValueError(f"Flag {flag!r} is managed by PhyloAI; remove from --tool-args.")
        cmd.extend(extra)
    cmd.append(str(input_file))
    return cmd


def _build_taper_command(
    msa_dir: Path, output_dir: Path, seq_type: str, cutoff: int,
    nt_dir: Path | None = None,
    taper_path: Path | None = None, julia_path: Path | None = None,
    threads: int = 4, tool_args: str | None = None,
    resume: bool = False, overwrite: bool = False,
    dry_run: bool = False, quiet: bool = False,
    table_format: str = "csv",
    show_masked_sites: bool = False,
) -> str:
    cmd = ["phyloai", "pretree", "filter", "taper",
           "--msa-dir", str(msa_dir), "--output-dir", str(output_dir),
           "--seq-type", seq_type, "--cutoff", str(cutoff),
           "--threads", str(threads)]
    if nt_dir:
        cmd.extend(["--nt-dir", str(nt_dir)])
    if taper_path:
        cmd.extend(["--taper-path", str(taper_path)])
    if julia_path:
        cmd.extend(["--julia-path", str(julia_path)])
    if tool_args:
        cmd.extend(["--tool-args", tool_args])
    if resume:
        cmd.append("--resume")
    if overwrite:
        cmd.append("--overwrite")
    if dry_run:
        cmd.append("--dry-run")
    if quiet:
        cmd.append("--quiet")
    if table_format != "csv":
        cmd.extend(["--table-format", table_format])
    if show_masked_sites:
        cmd.append("--show-masked-sites")
    return shlex.join(cmd)


def _run_taper_one(
    input_file: Path, output_file: Path, seq_type: str, cutoff: int,
    julia_exe: str, taper_script: str, tool_args: str | None,
) -> dict:
    t0 = time.monotonic()
    cmd = _build_taper_cmd(input_file, output_file, seq_type, cutoff, julia_exe, taper_script, tool_args)
    in_recs_all = list(SeqIO.parse(str(input_file), "fasta"))
    n_taxa_before = len(in_recs_all)
    length_before = len(in_recs_all[0].seq) if in_recs_all else 0
    output_file.parent.mkdir(parents=True, exist_ok=True)
    with open(output_file, "w") as fh:
        proc = subprocess.run(cmd, stdout=fh, stderr=subprocess.PIPE, text=True, timeout=86400)
    wall_time = round(time.monotonic() - t0, 2)
    n_taxa_after = 0
    length_after = 0
    out_recs_list: list = []
    if proc.returncode == 0 and output_file.exists():
        out_recs_list = list(SeqIO.parse(str(output_file), "fasta"))
        SeqIO.write(out_recs_list, str(output_file), "fasta")
        n_taxa_after = len(out_recs_list)
        length_after = len(out_recs_list[0].seq) if out_recs_list else 0
    new_mask_count = 0
    masked_taxa: list[dict] = []
    warnings_list: list[str] = []
    if proc.returncode == 0 and output_file.exists() and seq_type == "AA":
        in_recs = {rec.id: str(rec.seq) for rec in in_recs_all}
        out_recs_dict = {rec.id: str(rec.seq) for rec in (out_recs_list if out_recs_list else SeqIO.parse(str(output_file), "fasta"))}
        for taxon in in_recs:
            if taxon in out_recs_dict:
                taxon_mask_count = 0
                for i, (in_ch, out_ch) in enumerate(zip(in_recs[taxon], out_recs_dict[taxon])):
                    if in_ch != "X" and out_ch == "X":
                        taxon_mask_count += 1
                new_mask_count += taxon_mask_count
                if taxon_mask_count > 0:
                    masked_taxa.append({"taxon": taxon, "masked_sites": taxon_mask_count})
    if proc.returncode != 0:
        warnings_list.append(f"TAPER exited with code {proc.returncode}")
    return {
        "locus": logical_msa_locus_name(input_file),
        "status": "success" if proc.returncode == 0 else "failed",
        "returncode": proc.returncode,
        "cmd": " ".join(cmd),
        "stderr": proc.stderr if proc.stderr else "",
        "new_masked_sites": new_mask_count,
        "masked_taxa_count": len(masked_taxa),
        "masked_taxa": masked_taxa,
        "output": str(output_file),
        "n_taxa_before": n_taxa_before,
        "n_taxa_after": n_taxa_after,
        "length_before": length_before,
        "length_after": length_after,
        "wall_time": wall_time,
        "warnings": warnings_list,
    }


def _project_taper_masks_to_cds(
    aa_original_path: Path, aa_masked_path: Path,
    nt_input_path: Path, nt_output_path: Path,
) -> dict:
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord

    aa_orig = {rec.id: str(rec.seq) for rec in SeqIO.parse(str(aa_original_path), "fasta")}
    aa_masked = {rec.id: str(rec.seq) for rec in SeqIO.parse(str(aa_masked_path), "fasta")}
    nt_recs = {rec.id: str(rec.seq) for rec in SeqIO.parse(str(nt_input_path), "fasta")}

    if aa_orig.keys() != nt_recs.keys() or aa_masked.keys() != nt_recs.keys():
        raise ValueError(f"AA/NT taxa mismatch for {aa_masked_path.stem}")
    length = len(next(iter(aa_orig.values())))
    for nt_seq in nt_recs.values():
        if len(nt_seq) != length * 3:
            raise ValueError("NT alignment length != AA length * 3")

    projected = 0
    warnings_list: list[str] = []
    new_nt_records = []

    for taxon in aa_orig:
        orig_seq = aa_orig[taxon]
        masked_seq = aa_masked[taxon]
        nt_chars = list(nt_recs[taxon])
        for i, (orig_ch, mask_ch) in enumerate(zip(orig_seq, masked_seq)):
            codon_start = i * 3
            if orig_ch == "X" and mask_ch == "X":
                pass
            elif orig_ch in _STANDARD_AA and mask_ch == "X":
                original_codon = "".join(nt_chars[codon_start:codon_start + 3])
                if original_codon not in ("---", "NNN"):
                    nt_chars[codon_start:codon_start + 3] = ["N", "N", "N"]
                    projected += 1
            elif orig_ch == "X" and mask_ch != "X":
                warnings_list.append(
                    f"Original X at pos {i} for {taxon} changed to {mask_ch!r}"
                )
            elif orig_ch == "-" and mask_ch == "X":
                warnings_list.append(
                    f"Gap at pos {i} for {taxon} became X -- no CDS change"
                )
        new_nt_records.append(SeqRecord(Seq("".join(nt_chars)), id=taxon, description=""))
    nt_output_path.parent.mkdir(parents=True, exist_ok=True)
    SeqIO.write(new_nt_records, str(nt_output_path), "fasta")
    return {"projected_codons": projected, "warnings": warnings_list}


def _verify_taper_outputs(aa_path: Path, nt_path: Path | None) -> bool:
    aa_ok = validate_fasta_output(aa_path, require_aligned=True).ok
    if nt_path is not None:
        return aa_ok and validate_fasta_output(nt_path, require_aligned=True).ok
    return aa_ok


def run_taper(
    msa_dir: Path, output_dir: Path, *,
    seq_type: str = "auto", nt_dir: Path | None = None,
    cutoff: int = _TAPER_CUTOFF_DEFAULT,
    taper_path: Path | None = None, julia_path: Path | None = None,
    threads: int = 4, tool_args: str | None = None,
    resume: bool = False, overwrite: bool = False,
    dry_run: bool = False, quiet: bool = False,
    table_format: str = "csv",
    show_masked_sites: bool = False,
    progress_callback: Callable[[Path], None] | None = None,
) -> dict[str, Any]:
    start = time.monotonic()
    if overwrite and resume:
        raise ValueError("--overwrite and --resume are mutually exclusive.")

    output_dir = output_dir.resolve()
    env = ToolEnv()
    julia_exe = str(julia_path) if julia_path else str(env.require("julia"))
    taper_script = str(taper_path) if taper_path else str(env.require("correction_multi.jl"))

    julia_version = "unknown"
    try:
        proc = subprocess.run([julia_exe, "-v"], capture_output=True, text=True, timeout=30)
        if proc.returncode == 0 and proc.stdout.strip():
            raw = proc.stdout.strip().splitlines()[0].strip()
            m = re.match(r"julia version (\S+)", raw)
            julia_version = m.group(1) if m else raw
    except Exception:
        pass

    delimiter = _table_delimiter(table_format)
    suffix = _table_suffix(table_format)

    msa_map = scan_msa_dir(msa_dir)
    if not msa_map:
        raise ValueError(f"No valid MSA files in {msa_dir}")
    if seq_type == "auto":
        first = list(msa_map.values())[0]
        sample = list(SeqIO.parse(str(first), "fasta"))
        seq_type = resolve_seq_type([str(r.seq) for r in sample])[0]

    nt_map: dict[str, Path] = {}
    if nt_dir is not None:
        nt_map = scan_msa_dir(nt_dir)
        if not nt_map:
            raise ValueError(f"No valid NT MSA files in {nt_dir}")
    is_aa_cds = nt_dir is not None

    if is_aa_cds and seq_type != "AA":
        raise ValueError(
            "AA+CDS mode (--nt-dir) requires amino-acid MSA input.  "
            f"Detected --seq-type {seq_type}.  If your input is AA, pass "
            "--seq-type AA explicitly.  NT alignments are not compatible "
            "with AA+CDS projection."
        )

    params = {
        "msa_dir": str(msa_dir), "output_dir": str(output_dir),
        "nt_dir": str(nt_dir) if nt_dir else None,
        "seq_type": seq_type, "cutoff": cutoff,
        "taper_path": taper_script, "julia_path": julia_exe,
        "threads": threads, "tool_args": tool_args, "table_format": table_format,
        "show_masked_sites": show_masked_sites,
        "resume": resume, "overwrite": overwrite, "dry_run": dry_run, "quiet": quiet,
    }
    command = _build_taper_command(
        msa_dir, output_dir, seq_type, cutoff,
        nt_dir=nt_dir, taper_path=taper_path, julia_path=julia_path,
        threads=threads, tool_args=tool_args,
        resume=resume, overwrite=overwrite, dry_run=dry_run, quiet=quiet,
        table_format=table_format, show_masked_sites=show_masked_sites,
    )

    ckpt_path = output_dir / "checkpoint.json"
    locus_list = sorted(msa_map.keys())
    input_files = [msa_map[loc] for loc in locus_list]

    def _output_for(inp: Path) -> Path:
        prefix = "seqs/faa" if is_aa_cds else "seqs"
        return output_dir / prefix / inp.name

    def _nt_output_for(inp: Path) -> Path | None:
        if not is_aa_cds:
            return None
        locus = logical_msa_locus_name(inp)
        if locus in nt_map:
            return output_dir / "seqs" / "fna" / nt_map[locus].name
        return None

    checkpoint: Checkpoint | None = None
    resume_success_results: list[dict] = []

    if dry_run:
        cmds = [" ".join(_build_taper_cmd(inp, _output_for(inp), seq_type, cutoff, julia_exe, taper_script, tool_args)) for inp in input_files]
        return {
            "status": "success", "command": command, "wall_time": time.monotonic() - start,
            "tool_versions": {}, "params": params,
            "key_results": {"n_input": len(input_files)}, "error": None,
            "data": {"dry_run_cmds": cmds, "summary": {"n_input_files": len(input_files)}},
        }

    if resume:
        checkpoint = load_checkpoint(ckpt_path)
        validate_resume_params(checkpoint, params, step="pretree.filter.taper")
        to_run_ids, _ = plan_resume(checkpoint, _verify_taper_outputs)
        for task in checkpoint.tasks:
            if task.task_id in set(to_run_ids):
                continue
            if task.status == "success":
                resume_success_results.append({
                    "locus": task.task_id, "status": "success",
                    "output": task.outputs.get("aa", ""),
                    "nt_output": task.outputs.get("nt"),
                    "new_masked_sites": 0,
                })
        input_files = [Path(task.input) for task in checkpoint.tasks if task.task_id in set(to_run_ids)]
    else:
        _common_output_conflict(output_dir, overwrite)
        checkpoint = build_initial_checkpoint(
            step="pretree.filter.taper", command=command, params=params,
            inputs=input_files, output_for=_output_for, nt_output_for=_nt_output_for,
        )
        save_checkpoint_atomic(checkpoint, ckpt_path)

    to_run_ids = [logical_msa_locus_name(f) for f in input_files]
    if checkpoint and to_run_ids:
        for tid in to_run_ids:
            mark_task(checkpoint, tid, status="running")
        save_checkpoint_atomic(checkpoint, ckpt_path)

    if not dry_run:
        logs_dir = output_dir / "logs"
        logs_dir.mkdir(parents=True, exist_ok=True)

    file_results: list[dict] = list(resume_success_results)
    last_flush = time.monotonic()
    interrupted = False

    try:
        if to_run_ids:
            with ProcessPoolExecutor(max_workers=threads) as pool:
                futures = {}
                for inp in input_files:
                    fut = pool.submit(
                        _run_taper_one, input_file=inp, output_file=_output_for(inp),
                        seq_type=seq_type, cutoff=cutoff, julia_exe=julia_exe,
                        taper_script=taper_script, tool_args=tool_args,
                    )
                    futures[fut] = logical_msa_locus_name(inp)
                for fut in as_completed(futures):
                    task_id = futures[fut]
                    try:
                        result = fut.result()
                    except Exception as exc:
                        result = {"locus": task_id, "status": "failed", "reason": str(exc)[:200], "output": ""}
                    if is_aa_cds and result.get("status") == "success":
                        locus = result["locus"]
                        if locus in nt_map:
                            nt_out = _nt_output_for(msa_map[locus])
                            try:
                                proj = _project_taper_masks_to_cds(msa_map[locus], Path(result["output"]), nt_map[locus], nt_out)
                                result["nt_output"] = str(nt_out)
                                result["projected_codons"] = proj["projected_codons"]
                            except Exception as exc:
                                result["status"] = "failed"
                                result["nt_error"] = str(exc)[:200]
                    file_results.append(result)
                    if not dry_run and result.get("stderr"):
                        (logs_dir / f"{task_id}.log").write_text(result["stderr"])
                    if checkpoint:
                        mark_task(checkpoint, task_id, status=result.get("status", "failed"), reason=result.get("reason"))
                        now = time.monotonic()
                        if now - last_flush >= _CHECKPOINT_FLUSH_INTERVAL:
                            save_checkpoint_atomic(checkpoint, ckpt_path)
                            last_flush = now
                    if progress_callback:
                        progress_callback(Path(result.get("output", "")))
    except KeyboardInterrupt:
        interrupted = True
        if checkpoint:
            checkpoint.status = "interrupted"
            save_checkpoint_atomic(checkpoint, ckpt_path, fsync=True)
        raise

    if checkpoint:
        checkpoint.status = "success" if not interrupted else "interrupted"
        checkpoint.completed_at = None if interrupted else checkpoint.touch()
        save_checkpoint_atomic(checkpoint, ckpt_path, fsync=True)

    retained = [r for r in file_results if r.get("status") == "success"]
    dropped = [r for r in file_results if r.get("status") != "success"]

    if not dry_run:
        retained_csv = output_dir / f"retained_loci{suffix}"
        dropped_csv = output_dir / f"dropped_loci{suffix}"
        filter_decisions_csv = output_dir / f"filter_decisions{suffix}"
        _write_csv_table([{"locus": r["locus"]} for r in retained], retained_csv, ["locus"], delimiter)
        _write_csv_table([{"locus": r["locus"], "reason": r.get("reason", "")} for r in dropped], dropped_csv, ["locus", "reason"], delimiter)
        decision_columns = ["locus", "status", "new_masked_sites", "masked_taxa_count"]
        decisions = [
            {
                "locus": r.get("locus", ""),
                "status": r.get("status", ""),
                "new_masked_sites": r.get("new_masked_sites", 0),
                "masked_taxa_count": r.get("masked_taxa_count", 0),
                **({"masked_taxa_detail": "; ".join(f"{t['taxon']}:{t['masked_sites']}" for t in r.get("masked_taxa", []))} if show_masked_sites else {}),
            }
            for r in file_results
        ]
        if show_masked_sites:
            decision_columns.append("masked_taxa_detail")
        _write_csv_table(decisions, filter_decisions_csv, decision_columns, delimiter)
        taper_output_files = {
            "retained_loci": {"path": str(retained_csv), "description": "Loci that passed TAPER masking and were retained"},
            "dropped_loci": {"path": str(dropped_csv), "description": "Loci excluded for failing TAPER masking criteria"},
            "filter_decisions": {"path": str(filter_decisions_csv), "description": "Per-locus TAPER masking decisions with site-level detail"},
        }
    else:
        taper_output_files = {}

    wall_time = time.monotonic() - start
    total_masked_sites = sum(r.get("new_masked_sites", 0) for r in file_results)
    total_masked_taxa = sum(r.get("masked_taxa_count", 0) for r in file_results)
    masked_loci_count = sum(1 for r in retained if r.get("new_masked_sites", 0) > 0)

    files_list: list[dict] = []
    for r in file_results:
        locus = r["locus"]
        entry = {
            "locus": locus,
            "status": "retained" if r.get("status") == "success" else "dropped",
            "cmd": shlex.split(r.get("cmd", "")),
            "log_file": f"logs/{locus}.log",
            "n_taxa_before": r.get("n_taxa_before", 0),
            "n_taxa_after": r.get("n_taxa_after", 0),
            "length_before": r.get("length_before", 0),
            "length_after": r.get("length_after", 0),
            "masked_sites": r.get("new_masked_sites", 0),
            "wall_time": r.get("wall_time", 0),
            "warnings": r.get("warnings", []),
        }
        files_list.append(entry)

    payload = {
        "status": "success" if retained else "error",
        "command": command, "wall_time": round(wall_time, 2),
        "tool_versions": {"julia": julia_version, "correction_multi.jl": "1.0.0"},
        "params": params,
        "key_results": {
            "n_input": len(file_results), "n_retained": len(retained), "n_dropped": len(dropped),
            "total_masked_aa_sites": total_masked_sites,
            "total_masked_taxa": total_masked_taxa,
            "masked_loci": masked_loci_count,
        },
        "error": None if retained else "All loci failed TAPER.",
        "data": {
            "files": files_list,
            "output_files": taper_output_files,
            "summary": {
                "n_input": len(file_results),
                "n_retained": len(retained),
                "n_dropped": len(dropped),
                "total_masked_aa_sites": total_masked_sites,
                "total_masked_taxa": total_masked_taxa,
                "masked_loci": masked_loci_count,
            },
        },
    }
    if not dry_run:
        write_result_json(payload, output_dir)
    return payload


# --- TreeShrink ---

_TREESHRINK_MANAGED_FLAGS = {"-i", "-t", "-a", "-q", "-m", "-o", "-O"}


def _build_treeshrink_command(
    tree_dir: Path, output_dir: Path,
    msa_dir: Path | None = None, threshold: float = 0.05,
    treeshrink_mode: str = "auto", treeshrink_path: Path | None = None,
    tool_args: str | None = None, keep_work_dir: bool = False,
    overwrite: bool = False, dry_run: bool = False,
    quiet: bool = False, table_format: str = "csv",
) -> str:
    cmd = ["phyloai", "pretree", "filter", "treeshrink",
           "--tree-dir", str(tree_dir), "--output-dir", str(output_dir),
           "--threshold", str(threshold)]
    if msa_dir:
        cmd.extend(["--msa-dir", str(msa_dir)])
    if treeshrink_mode != "auto":
        cmd.extend(["--treeshrink-mode", treeshrink_mode])
    if treeshrink_path:
        cmd.extend(["--treeshrink-path", str(treeshrink_path)])
    if tool_args:
        if " " in tool_args:
            cmd.append(f"--tool-args '{tool_args}'")
        else:
            cmd.extend(["--tool-args", tool_args])
    if keep_work_dir:
        cmd.append("--keep-work-dir")
    if overwrite:
        cmd.append("--overwrite")
    if dry_run:
        cmd.append("--dry-run")
    if quiet:
        cmd.append("--quiet")
    if table_format != "csv":
        cmd.extend(["--table-format", table_format])
    return shlex.join(cmd)


def run_treeshrink(
    tree_dir: Path, output_dir: Path, *,
    msa_dir: Path | None = None, threshold: float = 0.05,
    treeshrink_mode: str = "auto", treeshrink_path: Path | None = None,
    tool_args: str | None = None, keep_work_dir: bool = False,
    overwrite: bool = False, dry_run: bool = False,
    quiet: bool = False, table_format: str = "csv",
) -> dict[str, Any]:
    start = time.monotonic()
    output_dir = output_dir.resolve()
    tool_paths = {"run_treeshrink.py": treeshrink_path} if treeshrink_path else {}
    env = ToolEnv(tool_paths=tool_paths)
    treeshrink_exe = str(env.require("run_treeshrink.py"))
    info = env._detect_tool("run_treeshrink.py", version_flag="--version")
    treeshrink_version = info.version or "unknown"

    delimiter = _table_delimiter(table_format)
    suffix = _table_suffix(table_format)
    tree_map = scan_tree_dir(tree_dir)
    if not tree_map:
        raise ValueError(f"No valid tree files in {tree_dir}")
    msa_map: dict[str, Path] = scan_msa_dir(msa_dir) if msa_dir else {}
    pairing = pair_msa_and_tree_maps(msa_map, list(tree_map.values()))

    params = {"tree_dir": str(tree_dir), "output_dir": str(output_dir),
              "msa_dir": str(msa_dir) if msa_dir else None,
              "threshold": threshold, "treeshrink_mode": treeshrink_mode,
              "treeshrink_path": str(treeshrink_path) if treeshrink_path else None,
              "tool_args": tool_args, "keep_work_dir": keep_work_dir,
              "overwrite": overwrite, "dry_run": dry_run, "quiet": quiet,
              "table_format": table_format}
    command = _build_treeshrink_command(
        tree_dir, output_dir,
        msa_dir=msa_dir, threshold=threshold, treeshrink_mode=treeshrink_mode,
        treeshrink_path=treeshrink_path, tool_args=tool_args,
        keep_work_dir=keep_work_dir, overwrite=overwrite, dry_run=dry_run,
        quiet=quiet, table_format=table_format,
    )

    if dry_run:
        work_dir_display = output_dir / "work" if keep_work_dir else Path("/tmp/treeshrink_tmp")
        cmd_display = [treeshrink_exe, "-i", str(work_dir_display / "input"), "-t", "input.tree", "-q", str(threshold)]
        if msa_dir:
            cmd_display.extend(["-a", "input.fasta"])
        if treeshrink_mode != "auto":
            cmd_display.extend(["-m", treeshrink_mode])
        return {"status": "success", "command": command, "wall_time": 0, "tool_versions": {"run_treeshrink.py": treeshrink_version},
                "params": params, "key_results": {"n_input": len(pairing.paired)}, "error": None,
                "data": {"dry_run_cmd": " ".join(cmd_display), "summary": {"n_input_files": len(pairing.paired)}}}

    _common_output_conflict(output_dir, overwrite)
    output_dir.mkdir(parents=True, exist_ok=True)

    work_dir = output_dir / "work" if keep_work_dir else Path(tempfile.mkdtemp(prefix="treeshrink_"))
    input_dir = work_dir / "input"

    for locus, (msa_path, tree_path) in pairing.paired.items():
        if tree_path is None:
            continue
        gene_dir = input_dir / locus
        gene_dir.mkdir(parents=True, exist_ok=True)
        shutil.copy2(tree_path, gene_dir / "input.tree")
        if msa_path is not None:
            shutil.copy2(msa_path, gene_dir / "input.fasta")

    cmd = [treeshrink_exe, "-i", str(input_dir), "-t", "input.tree", "-q", str(threshold)]
    if msa_dir:
        cmd.extend(["-a", "input.fasta"])
    if treeshrink_mode != "auto":
        cmd.extend(["-m", treeshrink_mode])
    if tool_args:
        extra = shlex.split(tool_args)
        for flag in _TREESHRINK_MANAGED_FLAGS:
            if flag in extra:
                raise ValueError(f"Flag {flag!r} is managed by PhyloAI; remove from --tool-args.")
        cmd.extend(extra)

    runner = Runner()
    ts_result = runner.run(cmd, tool_name="run_treeshrink.py")

    trees_out = output_dir / "trees"
    seqs_out = output_dir / "seqs"
    trees_out.mkdir(parents=True, exist_ok=True)

    file_results, retained, dropped, modified_loci, removed_taxa = [], [], [], [], []
    for locus, (msa_path, tree_path) in pairing.paired.items():
        if tree_path is None:
            dropped.append({"locus": locus, "reason": "no tree input"})
            continue
        src_tree = input_dir / locus / "output.tree"
        if src_tree.exists():
            dst_tree = trees_out / f"{locus}.tre"
            shutil.copy2(src_tree, dst_tree)
            entry = {"locus": locus, "status": "success", "output_tree": str(dst_tree)}
            try:
                from Bio import Phylo
                in_tree = Phylo.read(str(tree_path), "newick")
                out_tree = Phylo.read(str(src_tree), "newick")
                in_taxa = {c.name for c in in_tree.get_terminals()}
                out_taxa = {c.name for c in out_tree.get_terminals()}
                removed = in_taxa - out_taxa
                if removed:
                    modified_loci.append({"locus": locus, "removed_count": len(removed)})
                    for t in sorted(removed):
                        removed_taxa.append({"locus": locus, "taxon": t})
            except Exception:
                pass
            if msa_path:
                src_fa = input_dir / locus / "output.fasta"
                if src_fa.exists():
                    seqs_out.mkdir(parents=True, exist_ok=True)
                    dst_fa = seqs_out / f"{locus}.fa"
                    shutil.copy2(src_fa, dst_fa)
                    records = list(SeqIO.parse(str(dst_fa), "fasta"))
                    SeqIO.write(records, str(dst_fa), "fasta")
                    entry["output_msa"] = str(dst_fa)
            retained.append(entry)
            file_results.append(entry)
        else:
            dropped.append({"locus": locus, "reason": "output missing"})
            file_results.append({"locus": locus, "status": "failed", "reason": "output missing"})

    if not dry_run:
        ts_retained_csv = output_dir / f"retained_loci{suffix}"
        ts_dropped_csv = output_dir / f"dropped_loci{suffix}"
        ts_modified_csv = output_dir / f"modified_loci{suffix}"
        ts_removed_taxa_csv = output_dir / f"removed_taxa{suffix}"
        ts_filter_decisions_csv = output_dir / f"filter_decisions{suffix}"
        _write_csv_table([{"locus": r["locus"]} for r in retained], ts_retained_csv, ["locus"], delimiter)
        _write_csv_table(dropped, ts_dropped_csv, ["locus", "reason"], delimiter)
        _write_csv_table(modified_loci, ts_modified_csv, ["locus", "removed_count"], delimiter)
        _write_csv_table(removed_taxa, ts_removed_taxa_csv, ["locus", "taxon"], delimiter)
        decisions = [{"locus": r.get("locus", ""), "status": r.get("status", "failed"), "removed_count": sum(1 for t in removed_taxa if t["locus"] == r.get("locus", ""))} for r in file_results]
        _write_csv_table(decisions, ts_filter_decisions_csv, ["locus", "status", "removed_count"], delimiter)
        ts_output_files = {
            "retained_loci": {"path": str(ts_retained_csv), "description": "Loci retained after TreeShrink taxon pruning"},
            "dropped_loci": {"path": str(ts_dropped_csv), "description": "Loci fully excluded by TreeShrink"},
            "modified_loci": {"path": str(ts_modified_csv), "description": "Loci where some taxa were pruned but the locus was retained"},
            "removed_taxa": {"path": str(ts_removed_taxa_csv), "description": "Taxa removed by TreeShrink across all loci"},
            "filter_decisions": {"path": str(ts_filter_decisions_csv), "description": "Per-locus TreeShrink pruning decisions"},
        }
    else:
        ts_output_files = {}

    if not keep_work_dir:
        shutil.rmtree(work_dir, ignore_errors=True)

    wall_time = time.monotonic() - start
    modified_locus_names = {m["locus"] for m in modified_loci}
    results: list[dict] = []
    for r in file_results:
        locus = r.get("locus", "")
        status = r.get("status", "failed")
        entry: dict = {"locus": locus, "status": status if status != "failed" else "dropped"}
        if locus in modified_locus_names:
            entry["status"] = "modified"
            entry["removed_taxa"] = [t["taxon"] for t in removed_taxa if t["locus"] == locus]
        if r.get("output_tree"):
            entry["output_tree"] = r["output_tree"]
        if r.get("output_msa"):
            entry["output_msa"] = r["output_msa"]
        if r.get("reason"):
            entry["reason"] = r["reason"]
        results.append(entry)

    merged_stderr = "\n".join(p for p in (ts_result.stdout.strip(), ts_result.stderr.strip()) if p)
    payload = {
        "status": "success" if retained else "error",
        "command": command, "wall_time": round(wall_time, 2),
        "tool_versions": {"run_treeshrink.py": treeshrink_version}, "params": params,
        "key_results": {
            "n_input": len(pairing.paired), "n_retained": len(retained) - len(modified_loci),
            "n_modified": len(modified_loci), "n_dropped": len(dropped),
            "n_removed_taxa_total": len(removed_taxa),
        },
        "error": None if retained else "All loci failed.",
        "data": {
            "cmd": cmd,
            "tool_stderr": merged_stderr,
            "output_files": ts_output_files,
            "summary": {
                "n_input": len(pairing.paired),
                "n_retained": len(retained) - len(modified_loci),
                "n_modified": len(modified_loci),
                "n_dropped": len(dropped),
                "n_removed_taxa_total": len(removed_taxa),
            },
            "results": results,
        },
    }
    write_result_json(payload, output_dir)
    return payload


# --- Metrics rule filtering ---

import re as _re  # noqa: E402

_OP_PATTERN = _re.compile(r"^([^><=!]+?)(>=|<=|!=|==|>|<)(.+)$")
_NUMERIC_OPS = {">=": lambda a, b: a >= b, "<=": lambda a, b: a <= b, ">": lambda a, b: a > b, "<": lambda a, b: a < b, "==": lambda a, b: a == b, "!=": lambda a, b: a != b}


class FilterCondition:
    def __init__(self, col: str, op: str, value: float | str):
        self.col = col
        self.op = op
        self.value = value

    def __str__(self) -> str:
        value_repr = self.value if isinstance(self.value, float) else repr(self.value)
        return f"{self.col}{self.op}{value_repr}"

    def evaluate(self, row: dict) -> bool:
        raw = row.get(self.col, "")
        if raw in (None, "", "NA"):
            return False
        if isinstance(self.value, str):
            return _NUMERIC_OPS[self.op](str(raw).strip(), self.value)
        try:
            return _NUMERIC_OPS[self.op](float(raw), self.value)
        except (ValueError, TypeError):
            return False


def parse_keep_conditions(keep: str, known_columns: set[str]) -> list[FilterCondition]:
    conditions = []
    for part in keep.split(","):
        part = part.strip()
        if not part:
            continue
        m = _OP_PATTERN.match(part)
        if not m:
            raise ValueError(f"Malformed condition: {part!r}. Expected form: col>=val")
        col, op, val_str = m.group(1).strip(), m.group(2), m.group(3).strip()
        if col not in known_columns:
            raise ValueError(f"Unknown column {col!r} in --keep. Known: {sorted(known_columns)}")
        try:
            val: float | str = float(val_str)
        except ValueError:
            val = val_str.strip("\"'")
            if op not in ("==", "!="):
                raise ValueError(
                    f"String value {val!r} in condition {part!r} can only use "
                    f"== or != operators.  Numeric operators ({op}) require a "
                    "numeric value."
                )
        conditions.append(FilterCondition(col, op, val))
    return conditions


def _apply_metric_filters(
    rows: list[dict], conditions: list[FilterCondition], loci_column: str = "loci"
) -> tuple[list[dict], list[dict], dict[str, int]]:
    retained, dropped = [], []
    failure_counts = {str(c): 0 for c in conditions}
    for row in rows:
        failed_conditions = [c for c in conditions if not c.evaluate(row)]
        if not failed_conditions:
            retained.append(row)
        else:
            failures = [str(c) for c in failed_conditions]
            for failure in failures:
                failure_counts[failure] += 1
            dropped.append({**row, "_filter_reason": "FAIL: " + ", ".join(failures)})
    return retained, dropped, failure_counts


def _detect_input_delimiter(path: Path, input_format: str) -> str:
    if input_format == "csv":
        return ","
    if input_format == "tsv":
        return "\t"
    with open(path, newline="") as fh:
        sample = fh.read(4096)
    tabs = sample.count("\t")
    commas = sample.count(",")
    if tabs == 0 and commas == 0:
        raise ValueError(f"Cannot detect delimiter in {path}: no tabs or commas found.")
    if tabs > 0 and commas > 0:
        if max(tabs, commas) < 2 * min(tabs, commas):
            raise ValueError(
                f"Ambiguous delimiter in {path}: {tabs} tabs, {commas} commas. "
                "Use --input-format csv|tsv to specify explicitly."
            )
    return "\t" if tabs > commas else ","


def _build_metrics_filter_command(
    table_path: Path, output_dir: Path, keep: str,
    input_format: str = "auto", loci_column: str = "loci",
    msa_dir: Path | None = None, tree_dir: Path | None = None,
    copy: bool = False, overwrite: bool = False,
    dry_run: bool = False, quiet: bool = False,
    table_format: str = "csv",
) -> str:
    cmd = ["phyloai", "pretree", "filter", "metrics",
           "--table", str(table_path), "--output-dir", str(output_dir),
           "--keep", keep]
    if input_format != "auto":
        cmd.extend(["--input-format", input_format])
    if loci_column != "loci":
        cmd.extend(["--loci-column", loci_column])
    if msa_dir:
        cmd.extend(["--msa-dir", str(msa_dir)])
    if tree_dir:
        cmd.extend(["--tree-dir", str(tree_dir)])
    if copy:
        cmd.append("--copy")
    if overwrite:
        cmd.append("--overwrite")
    if dry_run:
        cmd.append("--dry-run")
    if quiet:
        cmd.append("--quiet")
    if table_format != "csv":
        cmd.extend(["--table-format", table_format])
    return shlex.join(cmd)


def run_metrics_filter(
    table_path: Path, output_dir: Path, *, keep: str,
    input_format: str = "auto", loci_column: str = "loci",
    msa_dir: Path | None = None, tree_dir: Path | None = None,
    copy: bool = False, overwrite: bool = False,
    dry_run: bool = False, quiet: bool = False,
    table_format: str = "csv",
) -> dict[str, Any]:
    start = time.monotonic()
    output_dir = output_dir.resolve()
    if copy and not msa_dir and not tree_dir:
        raise ValueError("--copy requires at least one of --msa-dir or --tree-dir.")
    delimiter_in = _detect_input_delimiter(table_path, input_format)
    delimiter_out = _table_delimiter(table_format)
    suffix = _table_suffix(table_format)
    rows = []
    with open(table_path, newline="") as fh:
        for row in csv.DictReader(fh, delimiter=delimiter_in):
            rows.append(row)
    if not rows:
        raise ValueError(f"No data rows in {table_path}")
    columns = list(rows[0].keys())
    conditions = parse_keep_conditions(keep, set(columns))
    retained, dropped, failure_counts = _apply_metric_filters(rows, conditions, loci_column)
    params = {"table": str(table_path), "output_dir": str(output_dir), "keep": keep, "input_format": input_format, "loci_column": loci_column, "msa_dir": str(msa_dir) if msa_dir else None, "tree_dir": str(tree_dir) if tree_dir else None, "copy": copy, "overwrite": overwrite, "dry_run": dry_run, "quiet": quiet, "table_format": table_format}
    command = _build_metrics_filter_command(
        table_path, output_dir, keep,
        input_format=input_format, loci_column=loci_column,
        msa_dir=msa_dir, tree_dir=tree_dir,
        copy=copy, overwrite=overwrite, dry_run=dry_run, quiet=quiet,
        table_format=table_format,
    )
    if dry_run:
        return {"status": "success", "command": command, "wall_time": 0, "tool_versions": {}, "params": params, "key_results": {"n_total": len(rows), "n_retained": len(retained), "n_dropped": len(dropped)}, "error": None, "data": {"condition_failure_counts": failure_counts}}
    _common_output_conflict(output_dir, overwrite)
    output_dir.mkdir(parents=True, exist_ok=True)
    mf_retained_csv = output_dir / f"retained_loci{suffix}"
    mf_dropped_csv = output_dir / f"dropped_loci{suffix}"
    mf_filter_decisions_csv = output_dir / f"filter_decisions{suffix}"
    _write_csv_table([{loci_column: r[loci_column]} for r in retained], mf_retained_csv, [loci_column], delimiter_out)
    _write_csv_table([{loci_column: d[loci_column], "reason": d.get("_filter_reason", "")} for d in dropped], mf_dropped_csv, [loci_column, "reason"], delimiter_out)
    decisions = [{loci_column: r[loci_column], "status": "retained", "reason": ""} for r in retained] + [{loci_column: d[loci_column], "status": "dropped", "reason": d.get("_filter_reason", "")} for d in dropped]
    _write_csv_table(decisions, mf_filter_decisions_csv, [loci_column, "status", "reason"], delimiter_out)
    copied_msa, copied_tree = 0, 0
    msa_map = scan_msa_dir(msa_dir) if msa_dir else {}
    tree_map = scan_tree_dir(tree_dir) if tree_dir else {}
    retained_set = {r[loci_column] for r in retained}

    if copy:
        if msa_map:
            (output_dir / "seqs").mkdir(parents=True, exist_ok=True)
            for locus in retained_set:
                if locus in msa_map:
                    shutil.copy2(msa_map[locus], output_dir / "seqs" / msa_map[locus].name)
                    copied_msa += 1
        if tree_map:
            (output_dir / "trees").mkdir(parents=True, exist_ok=True)
            for locus in retained_set:
                if locus in tree_map:
                    shutil.copy2(tree_map[locus], output_dir / "trees" / tree_map[locus].name)
                    copied_tree += 1

    msa_stats = _compute_retained_msa_stats(
        [msa_map[loc] for loc in retained_set if loc in msa_map]
    ) if msa_map else {}
    wall_time = time.monotonic() - start

    files_list: list[dict] = []
    for r in retained:
        locus = r.get(loci_column, "")
        files_list.append({"locus": locus, "status": "retained", "warnings": []})
    for d in dropped:
        locus = d.get(loci_column, "")
        files_list.append({
            "locus": locus,
            "status": "dropped",
            "warnings": [d.get("_filter_reason", "")],
        })

    payload = {
        "status": "success", "command": command, "wall_time": round(wall_time, 2),
        "tool_versions": {}, "params": params,
        "key_results": {
            "n_total": len(rows), "n_retained": len(retained),
            "n_dropped": len(dropped), "condition_failure_counts": failure_counts,
        },
        "error": None,
        "data": {
            "files": files_list,
            "output_files": {
                "retained_loci": {"path": str(mf_retained_csv), "description": "Loci that matched the metric-rule filter criteria"},
                "dropped_loci": {"path": str(mf_dropped_csv), "description": "Loci excluded for failing one or more metric-rule conditions"},
                "filter_decisions": {"path": str(mf_filter_decisions_csv), "description": "Per-locus metric values evaluated against the filtering rules"},
            },
            "summary": {
                "n_total": len(rows),
                "n_retained": len(retained),
                "n_dropped": len(dropped),
                "copied_msa": copied_msa,
                "copied_tree": copied_tree,
                "retained_msa_stats": msa_stats,
                "condition_failure_counts": failure_counts,
            },
        },
    }
    write_result_json(payload, output_dir)
    return payload


# --- Cluster-based filtering ---

import numpy as np  # noqa: E402


def _select_features(rows: list[dict], columns: list[str], metrics: str | None, exclude_regex: list[str], loci_column: str) -> tuple[list[str], list[dict]]:
    """Return (included_features, excluded_entries) where each excluded entry
    has ``column``, ``included`` (False), and ``reason``."""
    exclude_patterns = [_re.compile(p) for p in (exclude_regex or [])]
    all_numeric: list[str] = []
    excluded: list[dict] = []
    for col in columns:
        if col == loci_column or col == "DataType":
            excluded.append({"column": col, "included": False, "reason": "locus_id_or_DataType"})
            continue
        if any(pat.search(col) for pat in exclude_patterns):
            excluded.append({"column": col, "included": False, "reason": "exclude_regex"})
            continue
        vals = set()
        for row in rows:
            v = row.get(col, "")
            if v not in (None, "", "NA"):
                try:
                    vals.add(float(v))
                except (ValueError, TypeError):
                    pass
        if len(vals) <= 1:
            excluded.append({"column": col, "included": False, "reason": "constant_or_non_numeric"})
            continue
        all_numeric.append(col)
    if metrics and metrics != "all":
        requested = [m.strip() for m in metrics.split(",")]
        included = [c for c in requested if c in all_numeric]
        for c in all_numeric:
            if c not in included:
                excluded.append({"column": c, "included": False, "reason": "not_in_metrics_list"})
    else:
        included = all_numeric
    for c in included:
        excluded.append({"column": c, "included": True, "reason": ""})
    return included, excluded


def _extract_feature_matrix(rows: list[dict], features: list[str], loci_column: str) -> tuple[np.ndarray, list[str], list[dict]]:
    data, labels, valid_rows = [], [], []
    for row in rows:
        vals = []
        all_valid = True
        for f in features:
            v = row.get(f, "")
            if v in (None, "", "NA"):
                all_valid = False
                break
            try:
                vals.append(float(v))
            except (ValueError, TypeError):
                all_valid = False
                break
        if all_valid:
            data.append(vals)
            labels.append(row.get(loci_column, ""))
            valid_rows.append(row)
    return np.array(data, dtype=float), labels, valid_rows


def _scale_features(matrix: np.ndarray) -> np.ndarray:
    from sklearn.preprocessing import StandardScaler
    return StandardScaler().fit_transform(matrix)


def _reduce_pca(scaled: np.ndarray) -> np.ndarray:
    from sklearn.decomposition import PCA
    return PCA(n_components=min(3, scaled.shape[1])).fit_transform(scaled)


def _hierarchical_clustering(
    reduced: np.ndarray, n_clusters: int | None, max_clusters: int | None,
    linkage: str, distance: str, n_loci: int,
) -> tuple[int, np.ndarray, list[dict]]:
    from sklearn.cluster import AgglomerativeClustering
    from sklearn.metrics import silhouette_score, calinski_harabasz_score, davies_bouldin_score

    if max_clusters is None:
        max_clusters = min(30, max(6, int(np.ceil(np.sqrt(n_loci) / 3))))

    if n_clusters is not None:
        cl = AgglomerativeClustering(n_clusters=n_clusters, metric=distance, linkage=linkage)
        labels = cl.fit_predict(reduced)
        return n_clusters, labels, []

    k_min, k_max = 2, min(max_clusters, n_loci - 1)
    if k_max < k_min:
        return 1, np.zeros(n_loci, dtype=int), []

    results: list[dict] = []
    for k in range(k_min, k_max + 1):
        cl = AgglomerativeClustering(n_clusters=k, metric=distance, linkage=linkage)
        lbs = cl.fit_predict(reduced)
        sil = silhouette_score(reduced, lbs, metric=distance)
        ch = calinski_harabasz_score(reduced, lbs)
        db = davies_bouldin_score(reduced, lbs)
        results.append({"k": k, "silhouette": sil, "calinski_harabasz": ch, "davies_bouldin": db})

    ks = [r["k"] for r in results]
    sil_ordered = sorted(ks, key=lambda k: -next(r["silhouette"] for r in results if r["k"] == k))
    ch_ordered = sorted(ks, key=lambda k: -next(r["calinski_harabasz"] for r in results if r["k"] == k))
    db_ordered = sorted(ks, key=lambda k: next(r["davies_bouldin"] for r in results if r["k"] == k))

    rank_sums: dict[int, int] = {}
    for k in ks:
        rank_sums[k] = sil_ordered.index(k) + ch_ordered.index(k) + db_ordered.index(k)

    best_k = min(ks, key=lambda k: (rank_sums[k],
                                     -next(r["silhouette"] for r in results if r["k"] == k), k))

    cl = AgglomerativeClustering(n_clusters=best_k, metric=distance, linkage=linkage)
    labels = cl.fit_predict(reduced)
    return best_k, labels, results


def _select_best_umap_replicate(
    scaled: np.ndarray,
    n_replicates: int,
    base_random_state: int,
    n_neighbors: int,
    min_dist: float,
    n_clusters: int | None,
    max_clusters: int | None,
    linkage: str,
    distance: str,
    threads: int = 1,
) -> tuple[np.ndarray, int, int, list[dict], list[dict]]:
    from umap import UMAP

    replicate_rows: list[dict] = []
    for replicate_index in range(n_replicates):
        if n_replicates == 1:
            kwargs = dict(n_components=3, n_neighbors=n_neighbors, min_dist=min_dist, random_state=base_random_state, n_jobs=1)
        else:
            kwargs = dict(n_components=3, n_neighbors=n_neighbors, min_dist=min_dist, n_jobs=threads)
        reducer = UMAP(**kwargs)
        reduced = reducer.fit_transform(scaled)
        selected_k, labels, selection_rows = _hierarchical_clustering(
            reduced, n_clusters, max_clusters, linkage, distance, len(scaled)
        )
        if selection_rows:
            selected_row = next(row for row in selection_rows if row["k"] == selected_k)
            silhouette = selected_row["silhouette"]
            ch = selected_row["calinski_harabasz"]
            db = selected_row["davies_bouldin"]
        else:
            silhouette = float("nan")
            ch = float("nan")
            db = float("nan")
        replicate_rows.append({
            "replicate": replicate_index,
            "random_state": base_random_state + replicate_index if n_replicates == 1 else None,
            "selected_k": selected_k,
            "silhouette": silhouette,
            "calinski_harabasz": ch,
            "davies_bouldin": db,
        })

    valid_metric_rows = [row for row in replicate_rows if row["silhouette"] == row["silhouette"]]
    if valid_metric_rows:
        sil_order = sorted(valid_metric_rows, key=lambda row: (-row["silhouette"], row["replicate"]))
        ch_order = sorted(valid_metric_rows, key=lambda row: (-row["calinski_harabasz"], row["replicate"]))
        db_order = sorted(valid_metric_rows, key=lambda row: (row["davies_bouldin"], row["replicate"]))
        for row in valid_metric_rows:
            row["rank_sum"] = sil_order.index(row) + ch_order.index(row) + db_order.index(row)
        best_row = min(valid_metric_rows, key=lambda row: (row["rank_sum"], -row["silhouette"], row["davies_bouldin"], row["replicate"]))
        best_index = best_row["replicate"]
    else:
        for row in replicate_rows:
            row["rank_sum"] = "NA"
        best_index = 0

    if n_replicates == 1:
        best_reduced = reduced
    else:
        reducer = UMAP(n_components=3, n_neighbors=n_neighbors, min_dist=min_dist, n_jobs=threads)
        best_reduced = reducer.fit_transform(scaled)
    best_k, best_labels, best_selection_rows = _hierarchical_clustering(
        best_reduced, n_clusters, max_clusters, linkage, distance, len(scaled)
    )
    return best_reduced, best_k, best_index, replicate_rows, best_selection_rows


# --- Cluster diagnostics ---


def _generate_cluster_plots(
    reduced: np.ndarray, labels: np.ndarray, output_dir: Path, n_clusters: int,
    plots_dir: Path | None = None,
) -> list[str]:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plots_dir = plots_dir or (output_dir / "plots")
    plots_dir.mkdir(parents=True, exist_ok=True)
    paths = []
    cmap = plt.get_cmap("tab20")
    # 2D scatter
    fig, ax = plt.subplots(figsize=(10, 8))
    for c in range(n_clusters):
        mask = labels == c
        ax.scatter(reduced[mask, 0], reduced[mask, 1], color=cmap(c % 20), label=f"Cluster {c}", alpha=0.7, s=30)
    ax.set_xlabel("Dim 1")
    ax.set_ylabel("Dim 2")
    ax.set_title(f"Cluster Scatter (2D) -- k={n_clusters}")
    ax.legend(bbox_to_anchor=(1.05, 1), loc="upper left", fontsize=8)
    fig.tight_layout()
    pdf_2d = plots_dir / "cluster_2d.pdf"
    fig.savefig(pdf_2d, dpi=150, bbox_inches="tight")
    plt.close(fig)
    paths.append(str(pdf_2d))
    # 3D scatter
    if reduced.shape[1] >= 3:
        fig = plt.figure(figsize=(10, 8))
        ax = fig.add_subplot(111, projection="3d")
        for c in range(n_clusters):
            mask = labels == c
            ax.scatter(reduced[mask, 0], reduced[mask, 1], reduced[mask, 2],
                       color=cmap(c % 20), label=f"Cluster {c}", alpha=0.7, s=20)
        ax.set_xlabel("Dim 1")
        ax.set_ylabel("Dim 2")
        ax.set_zlabel("Dim 3")
        ax.set_title(f"Cluster Scatter (3D) -- k={n_clusters}")
        ax.legend(bbox_to_anchor=(1.05, 1), loc="upper left", fontsize=8)
        fig.tight_layout()
        pdf_3d = plots_dir / "cluster_3d.pdf"
        fig.savefig(pdf_3d, dpi=150, bbox_inches="tight")
        plt.close(fig)
        paths.append(str(pdf_3d))
    return paths


def _generate_cluster_metric_means(
    valid_rows: list[dict], labels: np.ndarray, features: list[str],
    output_dir: Path, table_format: str,
    plot_label_angle: float = 45.0,
    plots_dir: Path | None = None,
) -> str:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from sklearn.preprocessing import StandardScaler

    n_clusters = len(set(labels))
    cluster_means = np.zeros((n_clusters, len(features)))
    for c in range(n_clusters):
        indices = [i for i, lb in enumerate(labels) if lb == c]
        for j, f in enumerate(features):
            vals = [float(valid_rows[i].get(f, 0)) for i in indices if valid_rows[i].get(f, "") not in ("", "NA")]
            cluster_means[c, j] = np.mean(vals) if vals else 0.0

    means_rows = []
    for c in range(n_clusters):
        entry = {"cluster": c, "n_loci": int((labels == c).sum())}
        for j, f in enumerate(features):
            entry[f] = round(float(cluster_means[c, j]), 6)
        means_rows.append(entry)
    delimiter_out = _table_delimiter(table_format)
    suffix = _table_suffix(table_format)
    _write_csv_table(means_rows, output_dir / f"cluster_metric_means{suffix}",
                     ["cluster", "n_loci"] + features, delimiter_out)

    scaler = StandardScaler()
    heat_data = scaler.fit_transform(cluster_means)
    ncols, nrows_cells = len(features), n_clusters
    fig_w = max(10, ncols * 0.7)
    fig_h = max(3, nrows_cells * 0.6)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    im = ax.imshow(heat_data, aspect="auto", cmap="RdBu_r", interpolation="nearest")
    ax.set_xticks(range(ncols))
    ax.set_xticklabels(features, rotation=plot_label_angle, ha="right", fontsize=8)
    ax.set_yticks(range(nrows_cells))
    ax.set_yticklabels([f"Cluster {c}" for c in range(n_clusters)])
    ax.set_title("Standardized Per-Cluster Metric Means\n(blue=below mean, red=above mean)", fontsize=10)
    cbar = plt.colorbar(im, ax=ax, shrink=0.8)
    cbar.set_label("z-score", fontsize=8)
    # Annotate cell values
    for i in range(nrows_cells):
        for j in range(ncols):
            ax.text(j, i, f"{heat_data[i, j]:.2f}", ha="center", va="center",
                    fontsize=6, color="black" if abs(heat_data[i, j]) < 1.5 else "white")
    fig.tight_layout()
    pdf = (plots_dir or output_dir / "plots") / "cluster_metric_heatmap.pdf"
    pdf.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(pdf, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return str(pdf)


def _auto_drop_outlier_clusters(labels: np.ndarray, rows: list[dict], n_loci: int, outlier_metric: str, outlier_direction: str, max_drop_fraction: float) -> tuple[set[int], list[dict]]:
    max_drop = int(np.floor(n_loci * max_drop_fraction))
    if max_drop == 0:
        return set(), []
    n_clusters = len(set(labels))
    cluster_means = {}
    cluster_sizes = {}
    for c in range(n_clusters):
        indices = [i for i, lb in enumerate(labels) if lb == c]
        cluster_sizes[c] = len(indices)
        vals = []
        for i in indices:
            v = rows[i].get(outlier_metric, "")
            try:
                vals.append(float(v))
            except (ValueError, TypeError):
                pass
        cluster_means[c] = np.mean(vals) if vals else 0.0
    reverse = outlier_direction == "high"
    sorted_clusters = sorted(cluster_means.items(), key=lambda x: x[1], reverse=reverse)
    drop = set()
    dropped_count = 0
    for c, _ in sorted_clusters:
        if dropped_count + cluster_sizes[c] <= max_drop:
            drop.add(c)
            dropped_count += cluster_sizes[c]
        else:
            break
    return drop, []


def _generate_cluster_metric_boxplots(
    valid_rows: list[dict], labels: np.ndarray, features: list[str],
    output_dir: Path, n_clusters: int,
    plot_metrics_cols: int = 2,
    plots_dir: Path | None = None,
) -> list[str]:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    n_features = len(features)
    ncols = plot_metrics_cols
    nrows = int(math.ceil(n_features / ncols))

    box_dir = plots_dir or output_dir / "plots"
    box_dir.mkdir(parents=True, exist_ok=True)
    cmap = plt.get_cmap("tab10")
    fig, axes = plt.subplots(nrows, ncols, figsize=(ncols * 4, nrows * 2.8))
    if nrows * ncols == 1:
        axes = np.array([[axes]])
    elif nrows == 1:
        axes = np.array([axes])
    elif ncols == 1:
        axes = np.array([[ax] for ax in axes])
    for idx, feature in enumerate(features):
        r, c = divmod(idx, ncols)
        ax = axes[r, c]
        grouped = []
        cluster_colors = []
        for cl in range(n_clusters):
            values = []
            for i, lb in enumerate(labels):
                if lb != cl:
                    continue
                raw = valid_rows[i].get(feature, "")
                if raw in ("", "NA"):
                    continue
                try:
                    values.append(float(raw))
                except (TypeError, ValueError):
                    continue
            grouped.append(values)
            cluster_colors.append(cmap(cl % 10))
        bp = ax.boxplot(grouped, tick_labels=[f"C{c}" for c in range(n_clusters)], patch_artist=True)
        for patch, color in zip(bp["boxes"], cluster_colors):
            patch.set_facecolor(color)
            patch.set_alpha(0.5)
        ax.set_title(feature, fontsize=9)
        ax.tick_params(axis="x", labelsize=7)
    for idx in range(n_features, nrows * ncols):
        r, c = divmod(idx, ncols)
        axes[r, c].set_visible(False)
    fig.tight_layout(pad=2.0)
    pdf_path = box_dir / "cluster_metric_boxplots.pdf"
    fig.savefig(pdf_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return [str(pdf_path)]


def _write_outlier_diagnostics(
    valid_rows: list[dict], labels: np.ndarray, drop_clusters: set[int],
    features: list[str], output_dir: Path, table_format: str,
    outlier_boxplot_cols: int = 4,
    plots_dir: Path | None = None,
) -> list[str]:
    from scipy.stats import mannwhitneyu
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    delimiter_out = _table_delimiter(table_format)
    suffix = _table_suffix(table_format)
    outlier_flags = [label in drop_clusters for label in labels]

    comparison_rows = []
    wilcoxon_rows = []
    for feature in features:
        normal_vals = []
        outlier_vals = []
        for row, is_outlier in zip(valid_rows, outlier_flags):
            raw = row.get(feature, "")
            if raw in ("", "NA"):
                continue
            try:
                value = float(raw)
            except (TypeError, ValueError):
                continue
            if is_outlier:
                outlier_vals.append(value)
            else:
                normal_vals.append(value)
        comparison_rows.append({
            "metric": feature,
            "normal_mean": round(float(np.mean(normal_vals)), 6) if normal_vals else "NA",
            "normal_median": round(float(np.median(normal_vals)), 6) if normal_vals else "NA",
            "normal_std": round(float(np.std(normal_vals)), 6) if normal_vals else "NA",
            "normal_n": len(normal_vals),
            "outlier_mean": round(float(np.mean(outlier_vals)), 6) if outlier_vals else "NA",
            "outlier_median": round(float(np.median(outlier_vals)), 6) if outlier_vals else "NA",
            "outlier_std": round(float(np.std(outlier_vals)), 6) if outlier_vals else "NA",
            "outlier_n": len(outlier_vals),
        })
        if normal_vals and outlier_vals:
            stat = mannwhitneyu(normal_vals, outlier_vals, alternative="two-sided")
            direction = "outlier_higher" if np.mean(outlier_vals) > np.mean(normal_vals) else "outlier_lower"
            wilcoxon_rows.append({"metric": feature, "u_statistic": round(float(stat.statistic), 6), "p_value": round(float(stat.pvalue), 6), "direction": direction})
        else:
            wilcoxon_rows.append({"metric": feature, "u_statistic": "NA", "p_value": "NA", "direction": "insufficient_data"})

    _write_csv_table(comparison_rows, output_dir / f"outlier_comparison{suffix}", ["metric", "normal_mean", "normal_median", "normal_std", "normal_n", "outlier_mean", "outlier_median", "outlier_std", "outlier_n"], delimiter_out)
    _write_csv_table(wilcoxon_rows, output_dir / f"outlier_wilcoxon{suffix}", ["metric", "u_statistic", "p_value", "direction"], delimiter_out)

    # One figure with all features in a grid layout
    n_features = len(features)
    ncols = outlier_boxplot_cols
    nrows = int(math.ceil(n_features / ncols))
    box_dir = plots_dir or output_dir / "plots"
    box_dir.mkdir(parents=True, exist_ok=True)
    normal_color = "#a6cee3"
    outlier_color = "#fb9a99"
    fig, axes = plt.subplots(nrows, ncols, figsize=(ncols * 3.5, nrows * 2.8))
    if nrows * ncols == 1:
        axes = np.array([[axes]])
    elif nrows == 1:
        axes = np.array([axes])
    elif ncols == 1:
        axes = np.array([[ax] for ax in axes])
    for idx, feature in enumerate(features):
        r, c = divmod(idx, ncols)
        ax = axes[r, c]
        normal_vals = []
        outlier_vals = []
        for row, is_outlier in zip(valid_rows, outlier_flags):
            raw = row.get(feature, "")
            if raw in ("", "NA"):
                continue
            try:
                value = float(raw)
            except (TypeError, ValueError):
                continue
            (outlier_vals if is_outlier else normal_vals).append(value)
        bp = ax.boxplot([normal_vals, outlier_vals], tick_labels=["normal", "outlier"], patch_artist=True)
        bp["boxes"][0].set_facecolor(normal_color)
        bp["boxes"][0].set_alpha(0.6)
        bp["boxes"][1].set_facecolor(outlier_color)
        bp["boxes"][1].set_alpha(0.6)
        ax.set_title(feature, fontsize=9)
        ax.tick_params(axis="x", labelsize=7)
        # significance annotation — placed above the box
        wrow = next((row for row in wilcoxon_rows if row["metric"] == feature), None)
        if wrow and isinstance(wrow.get("p_value"), (int, float)):
            p = wrow["p_value"]
            if p < 0.001:
                sig = "***"
            elif p < 0.01:
                sig = "**"
            elif p < 0.05:
                sig = "*"
            else:
                sig = f"p={p:.3f}"
            all_vals = normal_vals + outlier_vals
            y_top = np.percentile(all_vals, 95) if all_vals else 0
            y_range = max(all_vals) - min(all_vals) if all_vals else 1
            y_pos = y_top + y_range * 0.08
            ax.annotate(sig, xy=(1.5, y_pos), ha="center", fontsize=9, fontweight="bold",
                        va="bottom")
    # Hide unused subplots
    for idx in range(n_features, nrows * ncols):
        r, c = divmod(idx, ncols)
        axes[r, c].set_visible(False)
    fig.tight_layout(pad=2.0)
    pdf_path = box_dir / "outlier_comparison_boxplots.pdf"
    fig.savefig(pdf_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return [str(pdf_path)]


# --- Cluster main entry point ---


def _build_cluster_command(
    table_path: Path, output_dir: Path,
    input_format: str = "auto", metrics: str | None = None,
    exclude_regex: list[str] | None = None, reduction: str = "pca",
    n_clusters: int | None = None, max_clusters: int | None = None,
    cluster_linkage: str = "ward", cluster_distance: str = "euclidean",
    drop_outlier_clusters: str = "none", outlier_metric: str = "average_BS",
    outlier_direction: str = "low", max_drop_fraction: float = 0.2,
    plot_metrics_cols: int = 2, plot_label_angle: float = 45.0,
    outlier_boxplot_cols: int = 4,
    umap_n_neighbors: int = 15, umap_min_dist: float = 0.001,
    umap_replicates: int = 1, umap_random_state: int = 42,
    threads: int = 1,
    msa_dir: Path | None = None, tree_dir: Path | None = None,
    copy: bool = False, overwrite: bool = False,
    dry_run: bool = False, quiet: bool = False,
    table_format: str = "csv",
) -> str:
    cmd = ["phyloai", "pretree", "filter", "cluster",
           "--table", str(table_path), "--output-dir", str(output_dir),
           "--reduction", reduction]
    if input_format != "auto":
        cmd.extend(["--input-format", input_format])
    if metrics:
        cmd.extend(["--metrics", metrics])
    if exclude_regex:
        for pattern in exclude_regex:
            cmd.extend(["--exclude-regex", pattern])
    if n_clusters is not None:
        cmd.extend(["--n-clusters", str(n_clusters)])
    if max_clusters is not None:
        cmd.extend(["--max-clusters", str(max_clusters)])
    if cluster_linkage != "ward":
        cmd.extend(["--cluster-linkage", cluster_linkage])
    if cluster_distance != "euclidean":
        cmd.extend(["--cluster-distance", cluster_distance])
    if drop_outlier_clusters != "none":
        cmd.extend(["--drop-outlier-clusters", drop_outlier_clusters])
    if outlier_metric != "average_BS":
        cmd.extend(["--outlier-metric", outlier_metric])
    if outlier_direction != "low":
        cmd.extend(["--outlier-direction", outlier_direction])
    if max_drop_fraction != 0.2:
        cmd.extend(["--max-drop-fraction", str(max_drop_fraction)])
    if plot_metrics_cols != 2:
        cmd.extend(["--plot-metrics-cols", str(plot_metrics_cols)])
    if plot_label_angle != 45.0:
        cmd.extend(["--plot-label-angle", str(plot_label_angle)])
    if outlier_boxplot_cols != 4:
        cmd.extend(["--outlier-boxplot-cols", str(outlier_boxplot_cols)])
    if umap_n_neighbors != 15:
        cmd.extend(["--umap-n-neighbors", str(umap_n_neighbors)])
    if umap_min_dist != 0.001:
        cmd.extend(["--umap-min-dist", str(umap_min_dist)])
    if umap_replicates != 1:
        cmd.extend(["--umap-replicates", str(umap_replicates)])
    if umap_random_state != 42:
        cmd.extend(["--umap-random-state", str(umap_random_state)])
    if threads != 1:
        cmd.extend(["--threads", str(threads)])
    if msa_dir:
        cmd.extend(["--msa-dir", str(msa_dir)])
    if tree_dir:
        cmd.extend(["--tree-dir", str(tree_dir)])
    if copy:
        cmd.append("--copy")
    if overwrite:
        cmd.append("--overwrite")
    if dry_run:
        cmd.append("--dry-run")
    if quiet:
        cmd.append("--quiet")
    if table_format != "csv":
        cmd.extend(["--table-format", table_format])
    return shlex.join(cmd)


def run_cluster_filter(
    table_path: Path, output_dir: Path, *,
    input_format: str = "auto", metrics: str | None = None,
    exclude_regex: list[str] | None = None, reduction: str = "pca",
    n_clusters: int | None = None, max_clusters: int | None = None,
    cluster_linkage: str = "ward", cluster_distance: str = "euclidean",
    drop_outlier_clusters: str = "none", outlier_metric: str = "average_BS",
    outlier_direction: str = "low", max_drop_fraction: float = 0.2,
    plot_metrics_cols: int = 2,
    plot_label_angle: float = 45.0,
    outlier_boxplot_cols: int = 4,
    umap_n_neighbors: int = 15, umap_min_dist: float = 0.001,
    umap_replicates: int = 1, umap_random_state: int = 42,
    threads: int = 1,
    msa_dir: Path | None = None, tree_dir: Path | None = None,
    copy: bool = False, overwrite: bool = False,
    dry_run: bool = False, quiet: bool = False,
    table_format: str = "csv",
) -> dict[str, Any]:
    start = time.monotonic()
    output_dir = output_dir.resolve()
    if reduction == "umap":
        try:
            import umap  # noqa: F401
        except ImportError:
            raise ImportError("umap-learn required for --reduction umap. pip install umap-learn")
    if cluster_linkage == "ward" and cluster_distance != "euclidean":
        raise ValueError("Ward linkage requires Euclidean distance.")
    delimiter_in = _detect_input_delimiter(table_path, input_format)
    delimiter_out = _table_delimiter(table_format)
    suffix = _table_suffix(table_format)
    loci_column = "loci"
    if not quiet:
        _console.print("[bold]Selecting features from metrics table ...[/bold]")
    rows = []
    with open(table_path, newline="") as fh:
        for row in csv.DictReader(fh, delimiter=delimiter_in):
            rows.append(row)
    if not rows:
        raise ValueError(f"No data rows in {table_path}")
    columns = list(rows[0].keys())
    features, feature_entries = _select_features(rows, columns, metrics, exclude_regex or [], loci_column)
    if len(features) < 2:
        raise ValueError(f"Need >=2 features; found {len(features)}.")
    if not quiet:
        _console.print(f"  {len(features)} features selected, {len(feature_entries) - len(features)} excluded")
    params = {"table": str(table_path), "output_dir": str(output_dir), "input_format": input_format, "metrics": metrics, "exclude_regex": list(exclude_regex) if exclude_regex else None, "reduction": reduction, "n_clusters": n_clusters, "max_clusters": max_clusters, "cluster_linkage": cluster_linkage, "cluster_distance": cluster_distance, "drop_outlier_clusters": drop_outlier_clusters, "outlier_metric": outlier_metric, "outlier_direction": outlier_direction, "max_drop_fraction": max_drop_fraction, "plot_metrics_cols": plot_metrics_cols, "plot_label_angle": plot_label_angle, "outlier_boxplot_cols": outlier_boxplot_cols, "umap_n_neighbors": umap_n_neighbors, "umap_min_dist": umap_min_dist, "umap_replicates": umap_replicates, "umap_random_state": umap_random_state, "threads": threads, "msa_dir": str(msa_dir) if msa_dir else None, "tree_dir": str(tree_dir) if tree_dir else None, "copy": copy, "overwrite": overwrite, "dry_run": dry_run, "quiet": quiet, "table_format": table_format}
    command = _build_cluster_command(
        table_path, output_dir,
        input_format=input_format, metrics=metrics, exclude_regex=exclude_regex,
        reduction=reduction, n_clusters=n_clusters, max_clusters=max_clusters,
        cluster_linkage=cluster_linkage, cluster_distance=cluster_distance,
        drop_outlier_clusters=drop_outlier_clusters, outlier_metric=outlier_metric,
        outlier_direction=outlier_direction, max_drop_fraction=max_drop_fraction,
        plot_metrics_cols=plot_metrics_cols, plot_label_angle=plot_label_angle,
        outlier_boxplot_cols=outlier_boxplot_cols,
        umap_n_neighbors=umap_n_neighbors, umap_min_dist=umap_min_dist,
        umap_replicates=umap_replicates, umap_random_state=umap_random_state,
        threads=threads, msa_dir=msa_dir, tree_dir=tree_dir,
        copy=copy, overwrite=overwrite, dry_run=dry_run, quiet=quiet,
        table_format=table_format,
    )
    if dry_run:
        k_range = [n_clusters, n_clusters] if n_clusters is not None else [2, min(max_clusters or min(30, max(6, int(np.ceil(np.sqrt(len(rows)) / 3)))), max(2, len(rows) - 1))]
        return {"status": "success", "command": command, "wall_time": 0, "tool_versions": {}, "params": params, "key_results": {"n_loci": len(rows), "n_features": len(features)}, "error": None, "data": {"features": features, "reduction": reduction, "k_range": k_range, "drop_outlier_clusters": drop_outlier_clusters, "copy": copy}}
    _common_output_conflict(output_dir, overwrite)
    output_dir.mkdir(parents=True, exist_ok=True)
    input_dir = output_dir / "01-input"
    reduction_dir = output_dir / "02-reduction"
    clustering_dir = output_dir / "03-clustering"
    diagnostics_dir = output_dir / "04-diagnostics"
    outlier_dir = output_dir / "05-outlier-drop"
    if not quiet:
        _console.print(f"[bold]Extracting feature matrix ({len(rows)} loci, {len(features)} features) ...[/bold]")
    matrix, valid_loci, valid_rows = _extract_feature_matrix(rows, features, loci_column)
    scaled = _scale_features(matrix)
    if not quiet:
        _console.print(f"  Valid loci: {len(valid_loci)}")
        _console.print(f"  01-input/features_used{suffix} — feature audit trail")
        _console.print(f"[bold]Reducing dimensions ({reduction.upper()}) ...[/bold]")
    if reduction == "pca":
        reduced = _reduce_pca(scaled)
        selected_replicate = None
        umap_replicate_rows: list[dict] = []
        if not quiet:
            _console.print("[bold]Clustering ...[/bold]")
        selected_k, labels, selection_rows = _hierarchical_clustering(reduced, n_clusters, max_clusters, cluster_linkage, cluster_distance, len(valid_loci))
    else:
        reduced, selected_k, selected_replicate, umap_replicate_rows, selection_rows = _select_best_umap_replicate(
            scaled,
            umap_replicates,
            umap_random_state,
            umap_n_neighbors,
            umap_min_dist,
            n_clusters,
            max_clusters,
            cluster_linkage,
            cluster_distance,
            threads=threads,
        )
    if not quiet:
        _console.print(
            f"  02-reduction/reduction{suffix} — reduced coordinates with cluster labels"
        )
        _console.print(
            f"  02-reduction/cluster_selection{suffix} — k-selection scoring metrics"
        )
        if reduction == "umap" and umap_replicate_rows:
            _console.print(
                f"  02-reduction/umap_replicates{suffix} — per-replicate UMAP metrics"
            )
        _console.print(
            f"[bold]Clustering (selected k={selected_k}) ...[/bold]"
        )
    _, labels, _ = _hierarchical_clustering(reduced, selected_k, max_clusters, cluster_linkage, cluster_distance, len(valid_loci))

    coord_names = ["PC1", "PC2", "PC3"] if reduction == "pca" else ["UMAP1", "UMAP2", "UMAP3"]

    if not quiet:
        _console.print(f"[bold]Writing cluster assignments ({selected_k} clusters) ...[/bold]")
    _write_csv_table(
        feature_entries,
        input_dir / f"features_used{suffix}", ["column", "included", "reason"], delimiter_out,
    )
    red_rows = []
    for i, locus in enumerate(valid_loci):
        row = {loci_column: locus}
        for j, cname in enumerate(coord_names):
            if j < reduced.shape[1]:
                row[cname] = round(float(reduced[i, j]), 6)
        row["cluster"] = int(labels[i])
        red_rows.append(row)
    _write_csv_table(red_rows, reduction_dir / f"reduction{suffix}",
                     [loci_column] + coord_names[:reduced.shape[1]] + ["cluster"], delimiter_out)
    if selection_rows:
        _write_csv_table(selection_rows, reduction_dir / f"cluster_selection{suffix}",
                         ["k", "silhouette", "calinski_harabasz", "davies_bouldin"], delimiter_out)
    if umap_replicate_rows:
        _write_csv_table(umap_replicate_rows, reduction_dir / f"umap_replicates{suffix}", ["replicate", "random_state", "selected_k", "silhouette", "calinski_harabasz", "davies_bouldin", "rank_sum"], delimiter_out)
    cluster_assign_rows = [{loci_column: valid_loci[i], "cluster": int(labels[i])} for i in range(len(valid_loci))]
    _write_csv_table(cluster_assign_rows, clustering_dir / f"clusters{suffix}", [loci_column, "cluster"], delimiter_out)
    if not quiet:
        _console.print(f"  03-clustering/clusters{suffix} — per-locus cluster assignments")
    _write_csv_table([{"cluster": int(c), "n_loci": int((labels == c).sum())} for c in range(selected_k)],
                     clustering_dir / f"cluster_summary{suffix}", ["cluster", "n_loci"], delimiter_out)
    if not quiet:
        _console.print(f"  03-clustering/cluster_summary{suffix} — cluster sizes")
    cluster_loci_dir = clustering_dir / "cluster_loci"
    cluster_loci_dir.mkdir(parents=True, exist_ok=True)
    for c in range(selected_k):
        mask = labels == c
        loci_in = [valid_loci[i] for i in range(len(valid_loci)) if mask[i]]
        _write_csv_table([{loci_column: locus} for locus in loci_in], cluster_loci_dir / f"cluster_{c}{suffix}", [loci_column], delimiter_out)
    if not quiet:
        _console.print("  03-clustering/cluster_loci/cluster_*.csv — per-cluster locus lists")

    if not quiet:
        _console.print("[bold]Generating diagnostic plots ...[/bold]")
    means_path = _generate_cluster_metric_means(valid_rows, labels, features, diagnostics_dir, table_format, plot_label_angle, diagnostics_dir)
    if not quiet:
        _console.print("  04-diagnostics/cluster_metric_means.csv + plots/cluster_metric_heatmap.pdf — per-cluster metric means and heatmap")
    boxplot_paths = _generate_cluster_metric_boxplots(valid_rows, labels, features, diagnostics_dir, selected_k, plot_metrics_cols, diagnostics_dir)
    if not quiet:
        _console.print("  04-diagnostics/plots/cluster_metric_boxplots.pdf — per-metric distributions by cluster")
    plot_paths = _generate_cluster_plots(reduced, labels, diagnostics_dir, selected_k, diagnostics_dir / "plots")
    if not quiet:
        _console.print("  04-diagnostics/plots/ — cluster_2d.pdf + cluster_3d.pdf scatter plots")

    drop_clusters: set[int] = set()
    outer_plot_paths: list[str] = []
    retained_set: list[str] = []
    if drop_outlier_clusters == "auto":
        if not quiet:
            _console.print(f"[bold]Running outlier cluster detection (by {outlier_metric}, direction={outlier_direction}) ...[/bold]")
        drop_clusters, _ = _auto_drop_outlier_clusters(labels, valid_rows, len(valid_rows), outlier_metric, outlier_direction, max_drop_fraction)
        if drop_clusters:
            retained_set = [valid_loci[i] for i in range(len(valid_loci)) if labels[i] not in drop_clusters]
            dropped_set = [valid_loci[i] for i in range(len(valid_loci)) if labels[i] in drop_clusters]
            _write_csv_table([{loci_column: locus} for locus in retained_set], outlier_dir / f"retained_loci{suffix}", [loci_column], delimiter_out)
            _write_csv_table([{loci_column: locus, "reason": "outlier_cluster"} for locus in dropped_set], outlier_dir / f"dropped_loci{suffix}", [loci_column, "reason"], delimiter_out)
            decisions = [{loci_column: valid_loci[i], "status": "dropped" if labels[i] in drop_clusters else "retained", "cluster": int(labels[i])} for i in range(len(valid_loci))]
            _write_csv_table(decisions, outlier_dir / f"filter_decisions{suffix}", [loci_column, "status", "cluster"], delimiter_out)
            outer_plot_paths = _write_outlier_diagnostics(valid_rows, labels, drop_clusters, features, outlier_dir, table_format, outlier_boxplot_cols, outlier_dir)
            if not quiet:
                _console.print(f"  05-outlier-drop/retained_loci{suffix} — {len(retained_set)} loci kept")
                _console.print(f"  05-outlier-drop/dropped_loci{suffix} — {len(dropped_set)} loci removed (clusters {sorted(drop_clusters)})")
                _console.print(f"  05-outlier-drop/filter_decisions{suffix} — per-locus keep/drop status")
                _console.print(f"  05-outlier-drop/outlier_comparison{suffix} + outlier_wilcoxon{suffix} — normal vs outlier stats")
                _console.print("  05-outlier-drop/plots/outlier_comparison_boxplots.pdf — full comparison figure")
            if copy:
                retained_locus_names = set(retained_set)
                if msa_dir:
                    msa_map = scan_msa_dir(msa_dir)
                    (outlier_dir / "seqs").mkdir(parents=True, exist_ok=True)
                    for locus in retained_locus_names:
                        if locus in msa_map:
                            shutil.copy2(msa_map[locus], outlier_dir / "seqs" / msa_map[locus].name)
                if tree_dir:
                    tree_map = scan_tree_dir(tree_dir)
                    (outlier_dir / "trees").mkdir(parents=True, exist_ok=True)
                    for locus in retained_locus_names:
                        if locus in tree_map:
                            shutil.copy2(tree_map[locus], outlier_dir / "trees" / tree_map[locus].name)
        elif copy:
            if not quiet:
                _console.print("[WARN] No outlier clusters dropped (all within max_drop_fraction). Copy skipped.")
            print("[WARN] No outlier clusters dropped (all within max_drop_fraction). Copy skipped.",
                  file=sys.stderr)

    wall_time = time.monotonic() - start
    n_dropped = int(sum((labels == c).sum() for c in drop_clusters)) if drop_clusters else 0
    n_retained = len(valid_loci) - n_dropped
    msa_stats = {}
    if msa_dir and retained_set:
        msa_paths = [msa_dir / f"{locus}.fa" for locus in retained_set]
        msa_stats = _compute_retained_msa_stats([p for p in msa_paths if p.exists()])

    files_list: list[dict] = []
    for i in range(len(valid_loci)):
        locus = valid_loci[i]
        is_dropped = labels[i] in drop_clusters if drop_clusters else False
        entry: dict = {"locus": locus, "status": "dropped" if is_dropped else "retained", "warnings": []}
        if is_dropped:
            entry["warnings"].append("outlier_cluster")
        files_list.append(entry)

    cluster_output_files = {
        "features_used": {"path": str(input_dir / f"features_used{suffix}"), "description": "Features selected for dimensionality reduction and clustering"},
        "reduction": {"path": str(reduction_dir / f"reduction{suffix}"), "description": "Dimensionality reduction coordinates and parameters"},
        "clusters": {"path": str(clustering_dir / f"clusters{suffix}"), "description": "Cluster assignment for each locus"},
        "cluster_summary": {"path": str(clustering_dir / f"cluster_summary{suffix}"), "description": "Per-cluster size and statistics summary"},
        "cluster_metric_means": {"path": str(diagnostics_dir / f"cluster_metric_means{suffix}"), "description": "Per-cluster mean value of each phylogenetic metric, used to characterise cluster properties"},
        "cluster_metric_heatmap": {"path": str(means_path), "description": "Heatmap of z-score normalised mean metric values: rows are clusters, columns are metrics, cell intensity shows how each cluster deviates from the global mean for that metric"},
    }
    if selection_rows:
        cluster_output_files["cluster_selection"] = {"path": str(reduction_dir / f"cluster_selection{suffix}"), "description": "Cluster count evaluation: scores from three indices (silhouette, Davies-Bouldin, Calinski-Harabasz) per candidate k and UMAP parameter combination, ranked to select the optimal number of clusters"}
    if umap_replicate_rows:
        cluster_output_files["umap_replicates"] = {"path": str(reduction_dir / f"umap_replicates{suffix}"), "description": "UMAP projection coordinates from replicate runs with different random seeds, used to assess projection stability and confirm the chosen cluster count"}
    for i, pdf in enumerate(plot_paths):
        cluster_output_files[f"cluster_scatter_{i}"] = {"path": pdf, "description": f"UMAP {'3D perspective view' if i > 0 else '2D scatter plot'}: loci coloured by cluster assignment, showing cluster separation in {'three' if i > 0 else 'two'} dimensions"}
    for i, pdf in enumerate(boxplot_paths):
        cluster_output_files[f"cluster_metric_boxplots_{i}"] = {"path": pdf, "description": "Per-metric boxplots comparing distributions across clusters; each subplot shows one metric, bars coloured by cluster"}
    for c in range(selected_k):
        cluster_output_files[f"cluster_{c}"] = {"path": str(cluster_loci_dir / f"cluster_{c}{suffix}"), "description": f"Loci assigned to cluster {c}"}
    if drop_outlier_clusters == "auto" and drop_clusters:
        cluster_output_files.update({
            "outlier_retained_loci": {"path": str(outlier_dir / f"retained_loci{suffix}"), "description": "Loci retained after outlier cluster removal"},
            "outlier_dropped_loci": {"path": str(outlier_dir / f"dropped_loci{suffix}"), "description": "Loci dropped as outlier clusters"},
            "outlier_filter_decisions": {"path": str(outlier_dir / f"filter_decisions{suffix}"), "description": "Per-locus outlier filtering decisions"},
            "outlier_comparison": {"path": str(outlier_dir / f"outlier_comparison{suffix}"), "description": "Per-metric descriptive statistics comparing retained and outlier clusters (mean, median, std for each group)"},
            "outlier_wilcoxon": {"path": str(outlier_dir / f"outlier_wilcoxon{suffix}"), "description": "Wilcoxon rank-sum test results per metric: U statistic, p-value, and direction of difference between retained and outlier clusters"},
        })
        for i, pdf in enumerate(outer_plot_paths):
            cluster_output_files[f"outlier_boxplots_{i}"] = {"path": pdf, "description": "Per-metric boxplots comparing retained vs outlier clusters; each subplot shows one metric with Wilcoxon test significance stars"}

    payload = {
        "status": "success", "command": command, "wall_time": round(wall_time, 2),
        "tool_versions": {}, "params": params,
        "key_results": {
            "n_loci": len(rows), "n_valid_loci": len(valid_loci),
            "n_features": len(features), "n_clusters": int(selected_k),
            "reduction": reduction,
            "selected_umap_replicate": selected_replicate,
            "n_retained": n_retained, "n_dropped": n_dropped,
        },
        "error": None,
        "data": {
            "files": files_list,
            "output_files": cluster_output_files,
            "summary": {
                "n_loci": len(rows),
                "n_valid_loci": len(valid_loci),
                "n_features": len(features),
                "n_clusters": int(selected_k),
                "reduction": reduction,
                "n_retained": n_retained,
                "n_dropped": n_dropped,
                "features": features,
                "cluster_sizes": {int(c): int((labels == c).sum()) for c in range(selected_k)},
                "drop_clusters": [int(c) for c in sorted(drop_clusters)],
                "retained_msa_stats": msa_stats,
                "umap_replicates": umap_replicate_rows,
            },
        },
    }
    write_result_json(payload, output_dir)
    return payload


# --- Symmetry test (symtest) ---

_EXPECTED_SYMTEST_COLUMNS = {"Name", "SymSig", "SymNon", "SymPval",
                              "MarSig", "MarNon", "MarPval",
                              "IntSig", "IntNon", "IntPval"}


def _parse_symtest_csv(fileobj) -> list[dict[str, Any]]:
    """Parse IQ-TREE ``.symtest.csv`` output into a list of per-partition dicts.

    Skips comment lines (starting with ``#``).  P-value columns are
    parsed as float; ``NA`` or unparseable values become None.
    """
    lines = [line for line in fileobj if not line.startswith("#")]
    if not lines:
        return []

    reader = csv.DictReader(io.StringIO("".join(lines)))
    if not reader.fieldnames:
        raise ValueError("Empty CSV header in symtest output")
    missing = _EXPECTED_SYMTEST_COLUMNS - set(reader.fieldnames)
    if missing:
        raise ValueError(
            f"Symtest CSV missing expected columns: {', '.join(sorted(missing))}"
        )

    results: list[dict[str, Any]] = []
    for row in reader:
        entry: dict[str, Any] = {}
        for key, value in row.items():
            if key in ("SymPval", "MarPval", "IntPval"):
                try:
                    entry[key] = float(value)
                except (ValueError, TypeError):
                    entry[key] = None
            elif key in ("SymSig", "SymNon", "MarSig", "MarNon", "IntSig", "IntNon"):
                try:
                    entry[key] = int(value)
                except (ValueError, TypeError):
                    entry[key] = 0
            else:
                entry[key] = value
        results.append(entry)
    return results


def _filter_by_symtest_pval(
    results: list[dict[str, Any]],
    symtest_type: str,
    threshold: float,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]]]:
    """Filter parsed symtest results by the selected p-value column.

    symtest_type is one of ``"Sym"``, ``"MAR"``, ``"INT"``.
    """
    pval_col = {"Sym": "SymPval", "MAR": "MarPval", "INT": "IntPval"}[symtest_type]

    retained: list[dict[str, Any]] = []
    dropped: list[dict[str, Any]] = []
    decisions: list[dict[str, Any]] = []

    for entry in results:
        locus = entry.get("Name", "")
        p_value = entry.get(pval_col)

        decision = {
            "locus": locus,
            "status": "",  # set below
            "p_value": p_value,
            "symtest_type": symtest_type,
            "sym_pval": entry.get("SymPval"),
            "mar_pval": entry.get("MarPval"),
            "int_pval": entry.get("IntPval"),
            "sym_sig": entry.get("SymSig", 0),
            "sym_non": entry.get("SymNon", 0),
            "mar_sig": entry.get("MarSig", 0),
            "mar_non": entry.get("MarNon", 0),
            "int_sig": entry.get("IntSig", 0),
            "int_non": entry.get("IntNon", 0),
        }

        if p_value is None:
            decision["status"] = "dropped"
            dropped.append({"locus": locus, "reason": "p_value is null"})
            decisions.append(decision)
        elif p_value >= threshold:
            decision["status"] = "retained"
            retained.append({"locus": locus, "p_value": p_value})
            decisions.append(decision)
        else:
            decision["status"] = "dropped"
            dropped.append({"locus": locus, "reason": f"{pval_col}={p_value} < {threshold}"})
            decisions.append(decision)

    return retained, dropped, decisions


def _build_symtest_supermatrix(
    msa_map: dict[str, Path],
) -> tuple[str, list[tuple[str, int, int]], str]:
    """Build a supermatrix string and partition list from a dict of MSAs.

    Returns ``(matrix_fasta_str, genes, prefix_type)`` where *genes* is
    ``[(name, start1, end1), ...]`` with 1-based positions.  Uses
    ``_read_msa`` from concat.py for format-agnostic reading.
    """
    from phyloai.pretree.concat import _read_msa

    if not msa_map:
        raise ValueError("No valid MSA files found")

    all_taxa: set[str] = set()
    msa_records: dict[str, tuple[list[str], list[str], int]] = {}

    for locus, path in sorted(msa_map.items()):
        taxa, seqs, length = _read_msa(path)
        if not taxa:
            continue
        all_taxa.update(taxa)
        msa_records[locus] = (taxa, seqs, length)

    if not msa_records:
        raise ValueError("No valid MSA files found")

    # Auto-detect seq_type from first 3 loci
    sample_seqs: list[str] = []
    for locus in list(msa_records.keys())[:3]:
        _, seqs, _ = msa_records[locus]
        sample_seqs.extend(seqs[:10])
    from phyloai.core.sequence_normalization import detect_seq_type
    seq_type = detect_seq_type(sample_seqs)

    if seq_type == "other":
        raise ValueError(
            "Could not determine sequence type from MSA files. "
            "Detected type: 'other'. Ensure input files contain "
            "valid AA or NT sequences."
        )

    prefix_type = "DNA" if seq_type in ("NT", "CODON") else "LG"

    # Build supermatrix
    matrix_parts: dict[str, list[str]] = {taxon: [] for taxon in all_taxa}
    genes: list[tuple[str, int, int]] = []
    pos = 1

    for locus, (taxa, seqs, length) in sorted(msa_records.items()):
        genes.append((locus, pos, pos + length - 1))
        pos += length
        taxon_to_seq = dict(zip(taxa, seqs))
        for taxon in all_taxa:
            seq = taxon_to_seq.get(taxon, "?" * length)
            matrix_parts[taxon].append(seq)

    taxon_order = sorted(all_taxa)
    lines: list[str] = []
    for taxon in taxon_order:
        seq = "".join(matrix_parts[taxon])
        wrapped = "\n".join(seq[i:i + 60] for i in range(0, len(seq), 60))
        lines.append(f">{taxon}\n{wrapped}")

    matrix_str = "\n".join(lines) + "\n"
    return matrix_str, genes, prefix_type


def _build_symtest_command(
    msa_dir: Path, output_dir: Path,
    symtest_type: str | None = None,
    symtest_pval: float = 0.05,
    symtest_keep_zero: bool = False,
    iqtree_path: Path | None = None,
    threads: int = 4,
    tree_dir: Path | None = None,
    table_format: str = "csv",
    dry_run: bool = False,
    overwrite: bool = False,
    quiet: bool = False,
) -> str:
    cmd = ["phyloai", "pretree", "filter", "symtest",
           "--msa-dir", str(msa_dir), "--output-dir", str(output_dir),
           "--symtest-pval", str(symtest_pval)]
    if symtest_type:
        cmd.extend(["--symtest-type", symtest_type])
    if symtest_keep_zero:
        cmd.append("--symtest-keep-zero")
    if iqtree_path:
        cmd.extend(["--iqtree-path", str(iqtree_path)])
    if threads != 4:
        cmd.extend(["--threads", str(threads)])
    if tree_dir:
        cmd.extend(["--tree-dir", str(tree_dir)])
    if table_format != "csv":
        cmd.extend(["--table-format", table_format])
    if dry_run:
        cmd.append("--dry-run")
    if overwrite:
        cmd.append("--overwrite")
    if quiet:
        cmd.append("--quiet")
    return shlex.join(cmd)


def run_symtest(
    msa_dir: Path, output_dir: Path, *,
    symtest_type: str | None = None,
    symtest_pval: float = 0.05,
    symtest_keep_zero: bool = False,
    iqtree_path: Path | None = None,
    threads: int = 4,
    tree_dir: Path | None = None,
    msa_map: dict[str, Path] | None = None,
    table_format: str = "csv",
    dry_run: bool = False,
    overwrite: bool = False,
    quiet: bool = False,
) -> dict[str, Any]:
    """Run IQ-TREE symmetry test on all MSAs and filter by p-value.

    *msa_map* is an optional pre-scanned ``{locus: path}`` dict; when not
    provided it is built via :func:`scan_msa_dir`.
    """
    start = time.monotonic()
    output_dir = output_dir.resolve()
    tool_paths = {"iqtree3": iqtree_path} if iqtree_path else {}
    env = ToolEnv(tool_paths=tool_paths)
    iqtree_exe = str(env.require("iqtree3"))
    info = env._detect_tool("iqtree3", version_flag="--version")
    iqtree_version = info.version or "unknown"

    msa_map = scan_msa_dir(msa_dir) if msa_map is None else msa_map
    if not msa_map:
        raise ValueError(f"No valid MSA files found in {msa_dir}")

    # Resolve symtest_type: None -> "Sym"
    resolved_type = symtest_type if symtest_type else "Sym"

    params = {
        "msa_dir": str(msa_dir), "output_dir": str(output_dir),
        "symtest_type": symtest_type,
        "symtest_pval": symtest_pval, "symtest_keep_zero": symtest_keep_zero,
        "iqtree_path": str(iqtree_path) if iqtree_path else None,
        "threads": threads, "tree_dir": str(tree_dir) if tree_dir else None,
        "overwrite": overwrite, "dry_run": dry_run, "quiet": quiet,
        "table_format": table_format,
    }
    command = _build_symtest_command(
        msa_dir, output_dir,
        symtest_type=symtest_type, symtest_pval=symtest_pval,
        symtest_keep_zero=symtest_keep_zero, iqtree_path=iqtree_path,
        threads=threads, tree_dir=tree_dir, table_format=table_format,
        dry_run=dry_run, overwrite=overwrite, quiet=quiet,
    )

    if dry_run:
        sym_extra = f" --symtest-type {symtest_type}" if symtest_type else ""
        return {
            "status": "success", "command": command, "wall_time": 0,
            "tool_versions": {"iqtree3": iqtree_version}, "params": params,
            "key_results": {"n_input": len(msa_map)},
            "error": None,
            "data": {"dry_run_cmd": f"{iqtree_exe} -s <matrix> -p <partitions> "
                     f"--symtest-only{sym_extra} -T {threads}"},
        }

    _common_output_conflict(output_dir, overwrite)

    # Build supermatrix + partition files in temp dir
    matrix_str, genes, prefix_type = _build_symtest_supermatrix(msa_map)

    work_dir = Path(tempfile.mkdtemp(prefix="symtest_"))
    try:
        matrix_path = work_dir / "symtest_matrix.fa"
        partitions_path = work_dir / "symtest_partitions.txt"
        matrix_path.write_text(matrix_str)
        from phyloai.pretree.concat import _write_partitions
        _write_partitions(partitions_path, genes, prefix_type)

        # Build IQ-TREE command (--symtest-pval NOT passed; used Python-side only)
        cmd = [
            iqtree_exe,
            "-s", str(matrix_path),
            "-p", str(partitions_path),
            "--symtest-only",
        ]
        if symtest_type:
            cmd.extend(["--symtest-type", symtest_type])
        if symtest_keep_zero:
            cmd.append("--symtest-keep-zero")
        if threads > 1:
            cmd.extend(["-T", str(threads)])

        # Run IQ-TREE
        runner = Runner()
        result = runner.run(cmd, tool_name="iqtree3", cwd=work_dir)

        if result.returncode != 0:
            raise RuntimeError(
                f"iqtree3 exited with code {result.returncode}.\n"
                f"STDERR:\n{result.stderr}"
            )

        # Parse symtest output
        symtest_csv = work_dir / "symtest_partitions.txt.symtest.csv"
        if not symtest_csv.exists():
            raise RuntimeError(
                f"Expected symtest output not found: {symtest_csv}\n"
                f"STDERR:\n{result.stderr}"
            )

        with open(symtest_csv) as fh:
            symtest_results = _parse_symtest_csv(fh)

        if not symtest_results:
            raise RuntimeError("Symtest CSV is empty -- no partitions parsed.")

        # Cross-validate CSV names against MSA map
        csv_names = {r["Name"] for r in symtest_results}
        msa_names = set(msa_map.keys())
        missing_in_csv = msa_names - csv_names
        extra_in_csv = csv_names - msa_names
        if missing_in_csv:
            raise RuntimeError(
                f"Loci in MSA directory but missing from symtest CSV: "
                f"{', '.join(sorted(missing_in_csv))}. "
                f"IQ-TREE may have dropped these partitions."
            )
        extra_in_csv_sorted = sorted(extra_in_csv) if extra_in_csv else []
        if extra_in_csv_sorted:
            import warnings
            warnings.warn(
                f"Partition names in symtest CSV not found in MSA map: "
                f"{', '.join(extra_in_csv_sorted)}. These will be skipped."
            )
            symtest_results = [r for r in symtest_results if r["Name"] in msa_names]

        # Filter
        retained, dropped, decisions = _filter_by_symtest_pval(
            symtest_results, resolved_type, symtest_pval,
        )

        # Copy retained MSAs
        seqs_out = output_dir / "seqs"
        seqs_out.mkdir(parents=True, exist_ok=True)
        for r in retained:
            locus = r["locus"]
            if locus in msa_map:
                shutil.copy2(msa_map[locus], seqs_out / f"{locus}.fa")

        # Copy retained trees (if --tree-dir)
        retained_tree_count = 0
        missed_tree_count = 0
        if tree_dir:
            tree_map = scan_tree_dir(tree_dir)
            trees_out = output_dir / "trees"
            trees_out.mkdir(parents=True, exist_ok=True)
            retained_loci = {r["locus"] for r in retained}
            for locus in sorted(retained_loci):
                if locus in tree_map:
                    shutil.copy2(tree_map[locus], trees_out / tree_map[locus].name)
                    retained_tree_count += 1
                else:
                    missed_tree_count += 1

        # Write decision tables
        delimiter = _table_delimiter(table_format)
        suffix = _table_suffix(table_format)

        sym_retained_csv = output_dir / f"retained_loci{suffix}"
        sym_dropped_csv = output_dir / f"dropped_loci{suffix}"
        sym_filter_decisions_csv = output_dir / f"filter_decisions{suffix}"

        _write_csv_table(
            retained, sym_retained_csv,
            ["locus"], delimiter,
        )
        _write_csv_table(
            dropped, sym_dropped_csv,
            ["locus", "reason"], delimiter,
        )
        _write_csv_table(
            decisions, sym_filter_decisions_csv,
            ["locus", "status", "p_value", "symtest_type",
             "sym_pval", "mar_pval", "int_pval",
             "sym_sig", "sym_non", "mar_sig", "mar_non", "int_sig", "int_non"],
            delimiter,
        )

        # MSA stats
        retained_paths = [seqs_out / f"{r['locus']}.fa" for r in retained]
        msa_stats = _compute_retained_msa_stats(retained_paths)

        wall_time = time.monotonic() - start
        merged_stderr = "\n".join(p for p in (result.stdout.strip(), result.stderr.strip()) if p)
        results: list[dict] = []
        for d in decisions:
            locus = d.get("locus", "")
            entry: dict = {"locus": locus, "status": d.get("status", "dropped")}
            if d.get("status") == "dropped":
                reason = next((dr["reason"] for dr in dropped if dr["locus"] == locus), "")
                if reason:
                    entry["reason"] = reason
            results.append(entry)

        payload = {
            "status": "success",
            "command": command,
            "wall_time": round(wall_time, 2),
            "tool_versions": {"iqtree3": iqtree_version},
            "params": params,
            "key_results": {
                "n_input": len(symtest_results),
                "n_retained": len(retained),
                "n_dropped": len(dropped),
                "p_value_threshold": symtest_pval,
                "symtest_type": resolved_type,
                "retained_trees_copied": retained_tree_count,
            },
            "error": None,
            "data": {
                "cmd": cmd,
                "tool_stderr": merged_stderr,
                "output_files": {
                    "retained_loci": {"path": str(sym_retained_csv), "description": "Loci that passed the symmetry test and were retained"},
                    "dropped_loci": {"path": str(sym_dropped_csv), "description": "Loci excluded for failing the symmetry test at the given p-value threshold"},
                    "filter_decisions": {"path": str(sym_filter_decisions_csv), "description": "Per-locus symmetry test results including p-values and decision status"},
                },
                "summary": {
                    "n_input": len(symtest_results),
                    "n_retained": len(retained),
                    "n_dropped": len(dropped),
                    "p_value_threshold": symtest_pval,
                    "symtest_type": resolved_type,
                    "retained_msa_stats": msa_stats,
                    "retained_tree_count": retained_tree_count,
                    "missed_tree_count": missed_tree_count,
                    "skipped_names": extra_in_csv_sorted,
                },
                "results": results,
            },
        }
        write_result_json(payload, output_dir)

    finally:
        shutil.rmtree(work_dir, ignore_errors=True)

    return payload
