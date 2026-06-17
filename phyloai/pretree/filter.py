"""Marker-level filtering: TAPER, TreeShrink, metric rules, and clustering."""

from __future__ import annotations

import csv
import datetime
import json
import shlex
import shutil
import subprocess
import sys
import tempfile
import time
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


def _write_result_json(payload: dict, output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    with open(output_dir / "result.json", "w") as fh:
        json.dump(payload, fh, indent=2)


def _common_output_conflict(output_dir: Path, overwrite: bool) -> None:
    if output_dir.exists() and any(output_dir.iterdir()):
        if not overwrite:
            raise ValueError(
                f"Output directory '{output_dir}' already exists and is non-empty. "
                "Use --overwrite to replace it."
            )
        shutil.rmtree(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)


def _write_filter_log(output_dir: Path, command: str, wall_time: float,
                      tool_versions: dict, success: bool) -> None:
    log_path = output_dir / "filter.log"
    with open(log_path, "a") as fh:
        fh.write(f"# {command}\n")
        fh.write(f"# Started: {datetime.datetime.now().isoformat()}\n")
        fh.write(f"# Tool versions: {json.dumps(tool_versions)}\n")
        fh.write(f"# Wall time: {wall_time:.1f}s\n")
        fh.write(f"# Exit code: {'0' if success else '1'}\n")
        fh.write("---\n")


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


def _run_taper_one(
    input_file: Path, output_file: Path, seq_type: str, cutoff: int,
    julia_exe: str, taper_script: str, tool_args: str | None,
) -> dict:
    cmd = _build_taper_cmd(input_file, output_file, seq_type, cutoff, julia_exe, taper_script, tool_args)
    output_file.parent.mkdir(parents=True, exist_ok=True)
    with open(output_file, "w") as fh:
        proc = subprocess.run(cmd, stdout=fh, stderr=subprocess.PIPE, text=True, timeout=86400)
    new_mask_count = 0
    masked_taxa: list[dict] = []
    if proc.returncode == 0 and output_file.exists() and seq_type == "AA":
        in_recs = {rec.id: str(rec.seq) for rec in SeqIO.parse(str(input_file), "fasta")}
        out_recs = {rec.id: str(rec.seq) for rec in SeqIO.parse(str(output_file), "fasta")}
        for taxon in in_recs:
            if taxon in out_recs:
                taxon_mask_count = 0
                for i, (in_ch, out_ch) in enumerate(zip(in_recs[taxon], out_recs[taxon])):
                    if in_ch != "X" and out_ch == "X":
                        taxon_mask_count += 1
                new_mask_count += taxon_mask_count
                if taxon_mask_count > 0:
                    masked_taxa.append({"taxon": taxon, "masked_sites": taxon_mask_count})
    return {
        "locus": logical_msa_locus_name(input_file),
        "status": "success" if proc.returncode == 0 else "failed",
        "returncode": proc.returncode,
        "cmd": " ".join(cmd),
        "stderr": proc.stderr[:500] if proc.stderr else "",
        "new_masked_sites": new_mask_count,
        "masked_taxa_count": len(masked_taxa),
        "masked_taxa": masked_taxa,
        "output": str(output_file),
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

    env = ToolEnv()
    julia_exe = str(julia_path) if julia_path else str(env.require("julia"))
    taper_script = str(taper_path) if taper_path else str(env.require("correction_multi.jl"))

    julia_version = "unknown"
    try:
        proc = subprocess.run([julia_exe, "-v"], capture_output=True, text=True, timeout=30)
        if proc.returncode == 0 and proc.stdout.strip():
            julia_version = proc.stdout.strip().splitlines()[0].strip()
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
        "msa_dir": str(msa_dir), "nt_dir": str(nt_dir) if nt_dir else None,
        "seq_type": seq_type, "cutoff": cutoff,
        "taper_path": taper_script, "julia_path": julia_exe,
        "threads": threads, "tool_args": tool_args, "table_format": table_format,
        "show_masked_sites": show_masked_sites,
    }
    command = f"phyloai pretree filter taper --msa-dir {msa_dir} --seq-type {seq_type} --cutoff {cutoff}"

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
        _write_csv_table([{"locus": r["locus"]} for r in retained], output_dir / f"retained_loci{suffix}", ["locus"], delimiter)
        _write_csv_table([{"locus": r["locus"], "reason": r.get("reason", "")} for r in dropped], output_dir / f"dropped_loci{suffix}", ["locus", "reason"], delimiter)
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
        _write_csv_table(decisions, output_dir / f"filter_decisions{suffix}", decision_columns, delimiter)

    wall_time = time.monotonic() - start
    total_masked_sites = sum(r.get("new_masked_sites", 0) for r in file_results)
    total_masked_taxa = sum(r.get("masked_taxa_count", 0) for r in file_results)
    masked_loci_count = sum(1 for r in retained if r.get("new_masked_sites", 0) > 0)
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
            "retained_loci": [r["locus"] for r in retained],
            "dropped_loci": [r["locus"] for r in dropped],
            "file_results": file_results,
            "retained_msa_stats": _compute_retained_msa_stats(
                [Path(r["output"]) for r in retained if r.get("output")]),
        },
    }
    if not dry_run:
        _write_result_json(payload, output_dir)
        _write_filter_log(output_dir, command, wall_time, payload["tool_versions"], payload["status"] == "success")
    return payload


# --- TreeShrink ---

_TREESHRINK_MANAGED_FLAGS = {"-i", "-t", "-a", "-q", "-m", "-o", "-O"}


def run_treeshrink(
    tree_dir: Path, output_dir: Path, *,
    msa_dir: Path | None = None, threshold: float = 0.05,
    treeshrink_mode: str = "auto", treeshrink_path: Path | None = None,
    tool_args: str | None = None, keep_work_dir: bool = False,
    overwrite: bool = False, dry_run: bool = False,
    quiet: bool = False, table_format: str = "csv",
) -> dict[str, Any]:
    start = time.monotonic()
    env = ToolEnv()
    treeshrink_exe = str(treeshrink_path) if treeshrink_path else str(env.require("run_treeshrink.py"))

    delimiter = _table_delimiter(table_format)
    suffix = _table_suffix(table_format)
    tree_map = scan_tree_dir(tree_dir)
    if not tree_map:
        raise ValueError(f"No valid tree files in {tree_dir}")
    msa_map: dict[str, Path] = scan_msa_dir(msa_dir) if msa_dir else {}
    pairing = pair_msa_and_tree_maps(msa_map, list(tree_map.values()))

    params = {"tree_dir": str(tree_dir), "msa_dir": str(msa_dir) if msa_dir else None,
              "threshold": threshold, "treeshrink_mode": treeshrink_mode, "table_format": table_format}
    command = f"phyloai pretree filter treeshrink --tree-dir {tree_dir} --threshold {threshold}"

    if dry_run:
        work_dir_display = output_dir / "work" if keep_work_dir else Path("/tmp/treeshrink_tmp")
        cmd_display = [treeshrink_exe, "-i", str(work_dir_display / "input"), "-t", "input.tree", "-q", str(threshold)]
        if msa_dir:
            cmd_display.extend(["-a", "input.fasta"])
        if treeshrink_mode != "auto":
            cmd_display.extend(["-m", treeshrink_mode])
        return {"status": "success", "command": command, "wall_time": 0, "tool_versions": {},
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
    runner.run(cmd, tool_name="run_treeshrink.py")

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
                    shutil.copy2(src_fa, seqs_out / f"{locus}.fa")
                    entry["output_msa"] = str(seqs_out / f"{locus}.fa")
            retained.append(entry)
            file_results.append(entry)
        else:
            dropped.append({"locus": locus, "reason": "output missing"})
            file_results.append({"locus": locus, "status": "failed", "reason": "output missing"})

    if not dry_run:
        _write_csv_table([{"locus": r["locus"]} for r in retained], output_dir / f"retained_loci{suffix}", ["locus"], delimiter)
        _write_csv_table(dropped, output_dir / f"dropped_loci{suffix}", ["locus", "reason"], delimiter)
        _write_csv_table(modified_loci, output_dir / f"modified_loci{suffix}", ["locus", "removed_count"], delimiter)
        _write_csv_table(removed_taxa, output_dir / f"removed_taxa{suffix}", ["locus", "taxon"], delimiter)
        decisions = [{"locus": r.get("locus", ""), "status": r.get("status", "failed"), "removed_count": sum(1 for t in removed_taxa if t["locus"] == r.get("locus", ""))} for r in file_results]
        _write_csv_table(decisions, output_dir / f"filter_decisions{suffix}", ["locus", "status", "removed_count"], delimiter)

    if not keep_work_dir:
        shutil.rmtree(work_dir, ignore_errors=True)

    msa_stats = _compute_retained_msa_stats(list(seqs_out.glob("*.fa"))) if msa_dir else {}
    wall_time = time.monotonic() - start
    payload = {"status": "success" if retained else "error", "command": command, "wall_time": round(wall_time, 2),
               "tool_versions": {"run_treeshrink.py": "unknown"}, "params": params,
               "key_results": {"n_input": len(pairing.paired), "n_retained": len(retained), "n_modified": len(modified_loci), "n_dropped": len(dropped), "n_removed_taxa_total": len(removed_taxa)},
               "error": None if retained else "All loci failed.", "data": {"retained_loci": [r["locus"] for r in retained], "modified_loci": modified_loci, "dropped_loci": dropped, "removed_taxa": removed_taxa, "retained_msa_stats": msa_stats}}
    _write_result_json(payload, output_dir)
    _write_filter_log(output_dir, command, wall_time, payload["tool_versions"], payload["status"] == "success")
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


def run_metrics_filter(
    table_path: Path, output_dir: Path, *, keep: str,
    input_format: str = "auto", loci_column: str = "loci",
    msa_dir: Path | None = None, tree_dir: Path | None = None,
    copy: bool = False, overwrite: bool = False,
    dry_run: bool = False, quiet: bool = False,
    table_format: str = "csv",
) -> dict[str, Any]:
    start = time.monotonic()
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
    params = {"table": str(table_path), "keep": keep, "input_format": input_format, "loci_column": loci_column, "copy": copy, "table_format": table_format}
    command = f"phyloai pretree filter metrics --table {table_path} --keep {keep!r}"
    if dry_run:
        return {"status": "success", "command": command, "wall_time": 0, "tool_versions": {}, "params": params, "key_results": {"n_total": len(rows), "n_retained": len(retained), "n_dropped": len(dropped)}, "error": None, "data": {"condition_failure_counts": failure_counts}}
    _common_output_conflict(output_dir, overwrite)
    output_dir.mkdir(parents=True, exist_ok=True)
    _write_csv_table([{loci_column: r[loci_column]} for r in retained], output_dir / f"retained_loci{suffix}", [loci_column], delimiter_out)
    _write_csv_table([{loci_column: d[loci_column], "reason": d.get("_filter_reason", "")} for d in dropped], output_dir / f"dropped_loci{suffix}", [loci_column, "reason"], delimiter_out)
    decisions = [{loci_column: r[loci_column], "status": "retained", "reason": ""} for r in retained] + [{loci_column: d[loci_column], "status": "dropped", "reason": d.get("_filter_reason", "")} for d in dropped]
    _write_csv_table(decisions, output_dir / f"filter_decisions{suffix}", [loci_column, "status", "reason"], delimiter_out)
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
    payload = {"status": "success", "command": command, "wall_time": round(wall_time, 2), "tool_versions": {}, "params": params, "key_results": {"n_total": len(rows), "n_retained": len(retained), "n_dropped": len(dropped), "condition_failure_counts": failure_counts}, "error": None, "data": {"copied_msa": copied_msa, "copied_tree": copied_tree, "retained_msa_stats": msa_stats, "condition_failure_counts": failure_counts}}
    _write_result_json(payload, output_dir)
    _write_filter_log(output_dir, command, wall_time, {}, True)
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


def _generate_cluster_plots(reduced: np.ndarray, labels: np.ndarray, output_dir: Path, n_clusters: int) -> list[str]:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plots_dir = output_dir / "cluster_plots"
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
    fig, ax = plt.subplots(figsize=(max(8, len(features) * 0.5), max(3, n_clusters * 0.5)))
    im = ax.imshow(heat_data, aspect="auto", cmap="RdBu_r", interpolation="nearest")
    ax.set_xticks(range(len(features)))
    ax.set_xticklabels(features, rotation=plot_label_angle, ha="right", fontsize=8)
    ax.set_yticks(range(n_clusters))
    ax.set_yticklabels([f"Cluster {c}" for c in range(n_clusters)])
    ax.set_title("Standardized Cluster Metric Means")
    plt.colorbar(im, ax=ax)
    fig.tight_layout()
    pdf = output_dir / "cluster_metric_heatmap.pdf"
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
    plot_metrics_rows: str = "auto", plot_metrics_cols: int = 2,
) -> list[str]:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    if plot_metrics_rows != "auto":
        nrows = int(plot_metrics_rows)
    else:
        nrows = 12 if n_clusters <= 6 else 6 if n_clusters <= 12 else 4 if n_clusters <= 20 else 2
    ncols = plot_metrics_cols
    per_page = nrows * ncols

    box_dir = output_dir / "cluster_metric_boxplots"
    box_dir.mkdir(parents=True, exist_ok=True)
    cmap = plt.get_cmap("tab10")
    pdf_paths: list[str] = []
    for page_start in range(0, len(features), per_page):
        page_features = features[page_start:page_start + per_page]
        n_page_rows = min(nrows, int(np.ceil(len(page_features) / ncols)))
        if len(page_features) <= ncols:
            n_page_rows = 1
            ncols_actual = len(page_features)
        else:
            ncols_actual = ncols
        fig, axes = plt.subplots(n_page_rows, ncols_actual, figsize=(ncols_actual * 4, n_page_rows * 2.8))
        if n_page_rows * ncols_actual == 1:
            axes = [[axes]]
        elif n_page_rows == 1:
            axes = [axes]
        elif ncols_actual == 1:
            axes = [[ax] for ax in axes]
        for idx, feature in enumerate(page_features):
            r, c = divmod(idx, ncols_actual)
            if r < len(axes) and c < len(axes[r]):
                ax = axes[r][c]
            else:
                continue
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
            bp = ax.boxplot(grouped, labels=[f"C{c}" for c in range(n_clusters)], patch_artist=True)
            for patch, color in zip(bp["boxes"], cluster_colors):
                patch.set_facecolor(color)
                patch.set_alpha(0.5)
            ax.set_title(feature, fontsize=9)
            ax.tick_params(axis="x", labelsize=7)
        for idx in range(len(page_features), n_page_rows * ncols_actual):
            r, c = divmod(idx, ncols_actual)
            if r < len(axes) and c < len(axes[r]):
                axes[r][c].set_visible(False)
        fig.tight_layout()
        pdf_path = box_dir / f"cluster_metric_boxplots_{(page_start // per_page) + 1:03d}.pdf"
        fig.savefig(pdf_path, dpi=150, bbox_inches="tight")
        plt.close(fig)
        pdf_paths.append(str(pdf_path))
    return pdf_paths


def _write_outlier_diagnostics(
    valid_rows: list[dict], labels: np.ndarray, drop_clusters: set[int],
    features: list[str], output_dir: Path, table_format: str,
    outlier_boxplot_rows: str = "auto", outlier_boxplot_cols: int = 4,
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

    if outlier_boxplot_rows != "auto":
        nrows = int(outlier_boxplot_rows)
    else:
        nrows = 12 if len(drop_clusters) <= 6 else 6
    ncols = outlier_boxplot_cols
    per_page = nrows * ncols

    box_dir = output_dir / "outlier_comparison_boxplots"
    box_dir.mkdir(parents=True, exist_ok=True)
    pdf_paths: list[str] = []
    normal_color = "#a6cee3"
    outlier_color = "#fb9a99"
    for page_start in range(0, len(features), per_page):
        page_features = features[page_start:page_start + per_page]
        n_page_rows = min(nrows, int(np.ceil(len(page_features) / ncols)))
        if len(page_features) <= ncols:
            n_page_rows = 1
            ncols_actual = len(page_features)
        else:
            ncols_actual = ncols
        fig, axes = plt.subplots(n_page_rows, ncols_actual, figsize=(ncols_actual * 3.5, n_page_rows * 2.8))
        if n_page_rows * ncols_actual == 1:
            axes = [[axes]]
        elif n_page_rows == 1:
            axes = [axes]
        elif ncols_actual == 1:
            axes = [[ax] for ax in axes]
        for idx, feature in enumerate(page_features):
            r, c = divmod(idx, ncols_actual)
            if r >= len(axes) or c >= len(axes[r]):
                continue
            ax = axes[r][c]
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
            bp = ax.boxplot([normal_vals, outlier_vals], labels=["normal", "outlier"], patch_artist=True)
            bp["boxes"][0].set_facecolor(normal_color)
            bp["boxes"][0].set_alpha(0.6)
            bp["boxes"][1].set_facecolor(outlier_color)
            bp["boxes"][1].set_alpha(0.6)
            ax.set_title(feature, fontsize=9)
            ax.tick_params(axis="x", labelsize=7)
            # significance annotation
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
                y_max = max(max(normal_vals) if normal_vals else 0, max(outlier_vals) if outlier_vals else 0)
                ax.annotate(sig, xy=(1.5, y_max * 1.02), ha="center", fontsize=10, fontweight="bold")
        for idx in range(len(page_features), n_page_rows * ncols_actual):
            r, c = divmod(idx, ncols_actual)
            if r < len(axes) and c < len(axes[r]):
                axes[r][c].set_visible(False)
        fig.tight_layout()
        pdf_path = box_dir / f"outlier_comparison_boxplots_{(page_start // per_page) + 1:03d}.pdf"
        fig.savefig(pdf_path, dpi=150, bbox_inches="tight")
        plt.close(fig)
        pdf_paths.append(str(pdf_path))
    return pdf_paths


# --- Cluster main entry point ---


def run_cluster_filter(
    table_path: Path, output_dir: Path, *,
    input_format: str = "auto", metrics: str | None = None,
    exclude_regex: list[str] | None = None, reduction: str = "pca",
    n_clusters: int | None = None, max_clusters: int | None = None,
    cluster_linkage: str = "ward", cluster_distance: str = "euclidean",
    drop_outlier_clusters: str = "none", outlier_metric: str = "average_BS",
    outlier_direction: str = "low", max_drop_fraction: float = 0.2,
    plot_metrics_rows: str = "auto", plot_metrics_cols: int = 2,
    plot_label_angle: float = 45.0,
    outlier_boxplot_rows: str = "auto", outlier_boxplot_cols: int = 4,
    umap_n_neighbors: int = 15, umap_min_dist: float = 0.001,
    umap_replicates: int = 1, umap_random_state: int = 42,
    threads: int = 1,
    msa_dir: Path | None = None, tree_dir: Path | None = None,
    copy: bool = False, overwrite: bool = False,
    dry_run: bool = False, quiet: bool = False,
    table_format: str = "csv",
) -> dict[str, Any]:
    start = time.monotonic()
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
    params = {"table": str(table_path), "metrics": metrics, "reduction": reduction, "n_clusters": n_clusters, "max_clusters": max_clusters, "cluster_linkage": cluster_linkage, "cluster_distance": cluster_distance, "drop_outlier_clusters": drop_outlier_clusters, "outlier_metric": outlier_metric, "outlier_direction": outlier_direction, "max_drop_fraction": max_drop_fraction, "plot_metrics_rows": plot_metrics_rows, "plot_metrics_cols": plot_metrics_cols, "plot_label_angle": plot_label_angle, "outlier_boxplot_rows": outlier_boxplot_rows, "outlier_boxplot_cols": outlier_boxplot_cols, "umap_n_neighbors": umap_n_neighbors, "umap_min_dist": umap_min_dist, "umap_replicates": umap_replicates, "umap_random_state": umap_random_state, "threads": threads, "table_format": table_format}
    command = f"phyloai pretree filter cluster --table {table_path} --reduction {reduction}"
    if dry_run:
        k_range = [n_clusters, n_clusters] if n_clusters is not None else [2, min(max_clusters or min(30, max(6, int(np.ceil(np.sqrt(len(rows)) / 3)))), max(2, len(rows) - 1))]
        return {"status": "success", "command": command, "wall_time": 0, "tool_versions": {}, "params": params, "key_results": {"n_loci": len(rows), "n_features": len(features)}, "error": None, "data": {"features": features, "reduction": reduction, "k_range": k_range, "drop_outlier_clusters": drop_outlier_clusters, "copy": copy}}
    _common_output_conflict(output_dir, overwrite)
    output_dir.mkdir(parents=True, exist_ok=True)
    if not quiet:
        _console.print(f"[bold]Extracting feature matrix ({len(rows)} loci, {len(features)} features) ...[/bold]")
    matrix, valid_loci, valid_rows = _extract_feature_matrix(rows, features, loci_column)
    scaled = _scale_features(matrix)
    if not quiet:
        _console.print(f"  Valid loci: {len(valid_loci)}")
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
            _console.print(f"[bold]Clustering (selected k={selected_k}) ...[/bold]")
        _, labels, _ = _hierarchical_clustering(reduced, selected_k, max_clusters, cluster_linkage, cluster_distance, len(valid_loci))

    coord_names = ["PC1", "PC2", "PC3"] if reduction == "pca" else ["UMAP1", "UMAP2", "UMAP3"]

    if not quiet:
        _console.print(f"[bold]Writing cluster assignments ({selected_k} clusters) ...[/bold]")
    _write_csv_table(
        feature_entries,
        output_dir / f"features_used{suffix}", ["column", "included", "reason"], delimiter_out,
    )
    red_rows = []
    for i, locus in enumerate(valid_loci):
        row = {loci_column: locus}
        for j, cname in enumerate(coord_names):
            if j < reduced.shape[1]:
                row[cname] = round(float(reduced[i, j]), 6)
        row["cluster"] = int(labels[i])
        red_rows.append(row)
    _write_csv_table(red_rows, output_dir / f"reduction{suffix}",
                     [loci_column] + coord_names[:reduced.shape[1]] + ["cluster"], delimiter_out)
    if not quiet:
        _console.print(f"  reduction{suffix} — reduced coordinates with cluster labels")
    if selection_rows:
        _write_csv_table(selection_rows, output_dir / f"cluster_selection{suffix}",
                         ["k", "silhouette", "calinski_harabasz", "davies_bouldin"], delimiter_out)
        if not quiet:
            _console.print(f"  cluster_selection{suffix} — k-selection scoring metrics")
    if umap_replicate_rows:
        _write_csv_table(umap_replicate_rows, output_dir / f"umap_replicates{suffix}", ["replicate", "random_state", "selected_k", "silhouette", "calinski_harabasz", "davies_bouldin", "rank_sum"], delimiter_out)
        if not quiet:
            _console.print(f"  umap_replicates{suffix} — per-replicate UMAP metrics")
    cluster_assign_rows = [{loci_column: valid_loci[i], "cluster": int(labels[i])} for i in range(len(valid_loci))]
    _write_csv_table(cluster_assign_rows, output_dir / f"clusters{suffix}", [loci_column, "cluster"], delimiter_out)
    if not quiet:
        _console.print(f"  clusters{suffix} — per-locus cluster assignments")
    _write_csv_table([{"cluster": int(c), "n_loci": int((labels == c).sum())} for c in range(selected_k)],
                     output_dir / f"cluster_summary{suffix}", ["cluster", "n_loci"], delimiter_out)
    if not quiet:
        _console.print(f"  cluster_summary{suffix} — cluster sizes")
    cluster_loci_dir = output_dir / "cluster_loci"
    cluster_loci_dir.mkdir(parents=True, exist_ok=True)
    for c in range(selected_k):
        mask = labels == c
        loci_in = [valid_loci[i] for i in range(len(valid_loci)) if mask[i]]
        _write_csv_table([{loci_column: locus} for locus in loci_in], cluster_loci_dir / f"cluster_{c}{suffix}", [loci_column], delimiter_out)

    if not quiet:
        _console.print("[bold]Generating diagnostic plots ...[/bold]")
    means_path = _generate_cluster_metric_means(valid_rows, labels, features, output_dir, table_format, plot_label_angle)
    if not quiet:
        _console.print(f"  cluster_metric_means{suffix} + cluster_metric_heatmap.pdf — per-cluster mean values (standardized heatmap shows relative patterns across metrics)")
    boxplot_paths = _generate_cluster_metric_boxplots(valid_rows, labels, features, output_dir, selected_k, plot_metrics_rows, plot_metrics_cols)
    if not quiet:
        _console.print(f"  cluster_metric_boxplots/ — per-metric distributions by cluster ({len(boxplot_paths)} page(s))")
    plot_paths = _generate_cluster_plots(reduced, labels, output_dir, selected_k)
    if not quiet:
        _console.print("  cluster_plots/ — 2D scatter + 3D scatter colored by cluster")

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
            _write_csv_table([{loci_column: locus} for locus in retained_set], output_dir / f"retained_loci{suffix}", [loci_column], delimiter_out)
            _write_csv_table([{loci_column: locus, "reason": "outlier_cluster"} for locus in dropped_set], output_dir / f"dropped_loci{suffix}", [loci_column, "reason"], delimiter_out)
            decisions = [{loci_column: valid_loci[i], "status": "dropped" if labels[i] in drop_clusters else "retained", "cluster": int(labels[i])} for i in range(len(valid_loci))]
            _write_csv_table(decisions, output_dir / f"filter_decisions{suffix}", [loci_column, "status", "cluster"], delimiter_out)
            outer_plot_paths = _write_outlier_diagnostics(valid_rows, labels, drop_clusters, features, output_dir, table_format, outlier_boxplot_rows, outlier_boxplot_cols)
            if not quiet:
                _console.print(f"  retained_loci{suffix} — {len(retained_set)} loci kept")
                _console.print(f"  dropped_loci{suffix} — {len(dropped_set)} loci removed (clusters {sorted(drop_clusters)})")
                _console.print(f"  filter_decisions{suffix} — per-locus keep/drop status")
                _console.print(f"  outlier_comparison{suffix} + outlier_wilcoxon{suffix} — normal vs outlier stats")
                _console.print(f"  outlier_comparison_boxplots/ — {len(outer_plot_paths)} page(s)")
            if copy:
                retained_locus_names = set(retained_set)
                if msa_dir:
                    msa_map = scan_msa_dir(msa_dir)
                    (output_dir / "seqs").mkdir(parents=True, exist_ok=True)
                    for locus in retained_locus_names:
                        if locus in msa_map:
                            shutil.copy2(msa_map[locus], output_dir / "seqs" / msa_map[locus].name)
                if tree_dir:
                    tree_map = scan_tree_dir(tree_dir)
                    (output_dir / "trees").mkdir(parents=True, exist_ok=True)
                    for locus in retained_locus_names:
                        if locus in tree_map:
                            shutil.copy2(tree_map[locus], output_dir / "trees" / tree_map[locus].name)
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
            "features": features,
            "cluster_sizes": {int(c): int((labels == c).sum()) for c in range(selected_k)},
            "drop_clusters": [int(c) for c in sorted(drop_clusters)],
            "retained_loci": retained_set,
            "retained_msa_stats": msa_stats,
            "plot_paths": plot_paths + boxplot_paths + [means_path] + outer_plot_paths,
            "umap_replicates": umap_replicate_rows,
        },
    }
    _write_result_json(payload, output_dir)
    _write_filter_log(output_dir, command, wall_time, {}, True)
    return payload
