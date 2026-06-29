"""Batch sequence alignment using MAFFT or MAGUS."""

from __future__ import annotations

import platform
import shlex
import shutil
import subprocess
import tempfile
import time
from collections.abc import Callable
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Any

from Bio import SeqIO

from phyloai.core.checkpoint import Checkpoint
from phyloai.core.schema import COMMON_ALIGNMENT_EXTENSIONS
from phyloai.core.sequence_output_validation import validate_fasta_output
from phyloai.core.sequence_normalization import detect_seq_type


MAFFT_METHODS = {"fftns1", "fftns2", "auto", "linsi", "einsi", "ginsi"}
MAGUS_METHODS = {"magus"}
ALL_METHODS = MAFFT_METHODS | MAGUS_METHODS

INPUT_EXTENSIONS = {
    ext for ext in COMMON_ALIGNMENT_EXTENSIONS
    if ext in {".fa", ".fas", ".fasta", ".faa", ".fna"}
}

STOP_CODONS = {"TAA", "TAG", "TGA"}

# Maximum interval between checkpoint flushes during the alignment loop.
# Flushing the full checkpoint after every gene is O(N^2) in disk I/O; throttling
# bounds the worst-case lost progress on interrupt to roughly this many seconds.
CHECKPOINT_FLUSH_INTERVAL = 2.0


def _detect_seq_type_from_files(files: list[Path], max_files: int = 3) -> str:
    sequences: list[str] = []
    for f in files[:max_files]:
        try:
            for rec in SeqIO.parse(str(f), "fasta"):
                sequences.append(str(rec.seq))
                if len(sequences) >= 10:
                    break
        except Exception:
            continue
        if len(sequences) >= 10:
            break
    return detect_seq_type(sequences) if sequences else "AA"


def _scan_input(
    seq_dir: Path,
) -> tuple[list[Path], list[dict[str, str]]]:
    found: list[Path] = []
    skipped: list[dict[str, str]] = []

    for entry in sorted(seq_dir.iterdir(), key=lambda p: p.name):
        if entry.is_dir():
            skipped.append({"path": str(entry), "reason": "directory"})
            continue
        if not entry.is_file():
            skipped.append({"path": str(entry), "reason": "not a file"})
            continue
        if entry.stat().st_size == 0:
            skipped.append({"path": str(entry), "reason": "empty file"})
            continue
        if entry.suffix.lower() not in INPUT_EXTENSIONS:
            skipped.append({"path": str(entry), "reason": "unrecognized extension"})
            continue
        found.append(entry)

    return found, skipped


def _build_mafft_cmd(
    input_file: Path,
    output_file: Path,
    method: str,
    executable: str = "mafft",
) -> list[str]:
    base = [executable]

    if method == "fftns1":
        base += ["--retree", "1", "--thread", "1"]
    elif method == "fftns2":
        base += ["--retree", "2", "--thread", "1"]
    elif method == "auto":
        base += ["--auto", "--thread", "1"]
    elif method == "linsi":
        base += ["--maxiterate", "1000", "--localpair", "--thread", "1"]
    elif method == "einsi":
        base += ["--maxiterate", "1000", "--genafpair", "--thread", "1"]
    elif method == "ginsi":
        base += ["--maxiterate", "1000", "--globalpair", "--thread", "1"]
    else:
        raise ValueError(f"Unknown MAFFT method: {method!r}")

    base.append(str(input_file))
    return base


def _build_magus_cmd(
    input_file: Path,
    output_file: Path,
    work_dir: Path,
    seq_type: str,
    tool_args: str | None,
    executable: str = "magus",
) -> list[str]:
    datatype = "protein" if seq_type == "AA" else "dna"

    internal: dict[str, str] = {
        "-i": str(input_file),
        "-o": str(output_file),
        "-d": str(work_dir),
        "--datatype": datatype,
        "-np": "1",
    }
    extra: list[str] = shlex.split(tool_args) if tool_args else []
    managed = {"-i", "-o", "-d", "--datatype", "-np"}
    blocked = managed.intersection(extra)
    if blocked:
        raise ValueError(f"--tool-args cannot include PhyloAI-managed MAGUS argument(s): {', '.join(sorted(blocked))}")

    cmd = [executable]
    for key, val in internal.items():
        cmd += [key, val]
    return cmd + extra


def _validate_cds(
    sequences: dict[str, str],
    n_aa_taxa: int | None = None,
    aa_taxa: set[str] | None = None,
) -> list[str]:
    warnings: list[str] = []

    if n_aa_taxa is None and aa_taxa is not None:
        n_aa_taxa = len(aa_taxa)

    if n_aa_taxa is not None and len(sequences) != n_aa_taxa:
        warnings.append(
            f"taxon count mismatch: CDS file has {len(sequences)} sequences "
            f"but AA alignment has {n_aa_taxa}"
        )

    if aa_taxa is not None:
        cds_taxa = set(sequences)
        if cds_taxa != aa_taxa:
            missing = sorted(aa_taxa - cds_taxa)
            extra = sorted(cds_taxa - aa_taxa)
            warnings.append(
                "taxon ID mismatch: "
                f"missing in CDS={missing}; extra in CDS={extra}"
            )

    for name, seq in sequences.items():
        seq_upper = seq.upper().replace("-", "")

        if len(seq_upper) % 3 != 0:
            warnings.append(
                f"{name}: CDS length {len(seq_upper)} is not a multiple of 3"
            )
            continue

        codons = [seq_upper[i:i+3] for i in range(0, len(seq_upper) - 3, 3)]
        for pos, codon in enumerate(codons):
            if codon in STOP_CODONS:
                warnings.append(
                    f"{name}: internal stop codon '{codon}' at codon position {pos + 1}"
                )
                break

    return warnings


def _validate_msa_output(path: Path) -> tuple[int, int, list[str]]:
    result = validate_fasta_output(path, require_aligned=True)
    return result.n_records, result.length, result.warnings


def verify_align_outputs(aa_path: Path, nt_path: Path | None) -> bool:
    aa_result = validate_fasta_output(aa_path, require_aligned=True)
    if not aa_result.ok:
        return False
    if nt_path is None:
        return True
    if not nt_path.exists() or nt_path.stat().st_size == 0:
        return False
    nt_result = validate_fasta_output(nt_path, require_aligned=True)
    return nt_result.ok


def _resolved_align_params(
    *,
    seq_dir: Path,
    output_dir: Path,
    method: str,
    resolved_seq_type: str,
    backtrans: bool,
    nt_dir: Path | None,
    threads: int,
    tool_args: str | None,
    mafft_path: str,
    magus_path: str,
    trimal_path: str,
    quiet: bool,
) -> dict[str, Any]:
    return {
        "seq_dir": str(seq_dir),
        "output_dir": str(output_dir),
        "method": method,
        "seq_type": resolved_seq_type,
        "backtrans": backtrans,
        "nt_dir": str(nt_dir) if nt_dir is not None else None,
        "threads": int(threads),
        "tool_args": tool_args,
        "mafft_path": mafft_path,
        "magus_path": magus_path,
        "trimal_path": trimal_path,
        "quiet": quiet,
    }


def _read_alignment_metrics(path: Path) -> tuple[int, int]:
    result = validate_fasta_output(path, require_aligned=True)
    return result.n_records, result.length


def reconstruct_align_result(
    *,
    checkpoint: Checkpoint,
    params: dict[str, Any],
    tool_versions: dict[str, str],
    wall_time: float,
    skipped_inputs: list[dict[str, str]],
    scan_warnings: list[str],
    file_results: list[dict[str, Any]] | None = None,
) -> dict[str, Any]:
    file_results_out: list[dict[str, Any]] = []
    failed: list[dict[str, str]] = []
    n_backtrans = 0

    _method = params.get("method", "")
    _mafft_exe = params.get("mafft_path", "mafft")
    _magus_exe = params.get("magus_path", "magus")
    _seq_type = params.get("seq_type", "AA")
    _tool_args = params.get("tool_args")

    _wall_map: dict[str, float] = {}
    if file_results:
        _wall_map = {r["input"]: r.get("wall_time", 0.0) for r in file_results}

    for task in checkpoint.tasks:
        if task.status != "success":
            failed.append({"path": task.input, "reason": task.reason or task.status})
            continue

        aa_path = Path(task.outputs["aa"]) if task.outputs.get("aa") else None
        nt_path = Path(task.outputs["nt"]) if task.outputs.get("nt") else None
        if aa_path is None or not aa_path.exists():
            failed.append({"path": task.input, "reason": "missing AA output"})
            continue

        n_taxa, alignment_length = _read_alignment_metrics(aa_path)
        if nt_path is not None and nt_path.exists():
            n_backtrans += 1

        input_path = Path(task.input)
        if _method in MAFFT_METHODS:
            cmd_list = _build_mafft_cmd(input_path, aa_path, _method, executable=_mafft_exe)
        else:
            work_dir = aa_path.parent / f".{input_path.stem}.magus"
            cmd_list = _build_magus_cmd(input_path, aa_path, work_dir, _seq_type, _tool_args, executable=_magus_exe)

        file_results_out.append(
            {
                "input": task.input,
                "output_aa": str(aa_path),
                "output_nt": str(nt_path) if nt_path else None,
                "n_taxa": n_taxa,
                "alignment_length": alignment_length,
                "wall_time": _wall_map.get(task.input, -1.0),
                "warnings": [],
                "cmd": cmd_list,
                "log_file": f"logs/{input_path.stem}.log",
            }
        )
        log_path = aa_path.parent.parent / "logs" / f"{input_path.stem}.log"
        if not log_path.exists():
            log_path.parent.mkdir(parents=True, exist_ok=True)
            log_path.write_text("# resumed from checkpoint — original stderr unavailable\n")

    aligned_lengths = [r["alignment_length"] for r in file_results_out if r["alignment_length"]]
    aligned_taxa = [r["n_taxa"] for r in file_results_out if r["n_taxa"]]
    mean_len = round(sum(aligned_lengths) / len(aligned_lengths), 1) if aligned_lengths else 0.0
    mean_taxa = round(sum(aligned_taxa) / len(aligned_taxa), 1) if aligned_taxa else 0.0
    skipped = list(skipped_inputs) + failed

    return {
        "status": "success" if file_results_out else "error",
        "command": checkpoint.command,
        "wall_time": wall_time,
        "tool_versions": tool_versions,
        "params": params,
        "key_results": {
            "n_aligned": len(file_results_out),
            "n_skipped": len(skipped),
            "method": params.get("method"),
            "backtrans": params.get("backtrans", False),
            "mean_alignment_length": mean_len,
            "mean_n_taxa": mean_taxa,
        },
        "error": None if file_results_out else "No genes were aligned.",
        "data": {
            "summary": {
                "n_input_files": len(checkpoint.tasks) + len(skipped_inputs),
                "n_aligned": len(file_results_out),
                "n_backtrans": n_backtrans,
                "n_skipped": len(skipped),
            },
            "files": file_results_out,
            "skipped": skipped,
            "warnings": list(scan_warnings),
        },
    }


def _align_one(
    gene_path: Path,
    output_dir: Path,
    method: str,
    seq_type: str,
    tool_args: str | None,
    dry_run: bool,
    mafft_executable: str = "mafft",
    magus_executable: str = "magus",
) -> dict[str, Any]:
    out_aa = output_dir / f"{gene_path.stem}.fa"

    if method in MAFFT_METHODS:
        cmd = _build_mafft_cmd(gene_path, out_aa, method, executable=mafft_executable)
    else:
        work_dir = (
            output_dir / f".{gene_path.stem}.magus-dry-run"
            if dry_run
            else Path(tempfile.mkdtemp(prefix="phyloai_magus_"))
        )
        cmd = _build_magus_cmd(gene_path, out_aa, work_dir, seq_type, tool_args, executable=magus_executable)

    if dry_run:
        return {"status": "dry_run", "input": str(gene_path), "cmd": cmd, "tool_cmd": " ".join(cmd), "warnings": [], "wall_time": 0.0}

    start = time.monotonic()
    try:
        if method in MAFFT_METHODS:
            proc = subprocess.run(
                cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
            )
            wall_time = time.monotonic() - start
            if proc.returncode != 0:
                return {
                    "status": "skipped",
                    "input": str(gene_path),
                    "reason": f"mafft exited with code {proc.returncode}: {proc.stderr[:200]}",
                    "tool_cmd": " ".join(cmd),
                    "tool_stdout": proc.stdout,
                    "tool_stderr": proc.stderr,
                    "wall_time": wall_time,
                }
            out_aa.parent.mkdir(parents=True, exist_ok=True)
            out_aa.write_text(proc.stdout)
            stdout, stderr = "", proc.stderr
        else:
            out_aa.parent.mkdir(parents=True, exist_ok=True)
            try:
                proc = subprocess.run(
                    cmd,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    text=True,
                )
                wall_time = time.monotonic() - start
            finally:
                _cleanup_magus_workdir(cmd)
            if proc.returncode != 0:
                return {
                    "status": "skipped",
                    "input": str(gene_path),
                    "reason": f"magus exited with code {proc.returncode}: {proc.stderr[:200]}",
                    "tool_cmd": " ".join(cmd),
                    "tool_stdout": proc.stdout,
                    "tool_stderr": proc.stderr,
                    "wall_time": wall_time,
                }
            stdout, stderr = "", proc.stderr

    except Exception as exc:
        return {
            "status": "skipped",
            "input": str(gene_path),
            "reason": str(exc),
            "tool_cmd": " ".join(cmd),
            "tool_stdout": "",
            "tool_stderr": "",
            "wall_time": time.monotonic() - start,
        }

    n_taxa, alignment_length, warnings = _validate_msa_output(out_aa)
    if warnings:
        return {
            "status": "skipped",
            "input": str(gene_path),
            "reason": "; ".join(warnings),
            "output_aa": str(out_aa),
            "tool_cmd": " ".join(cmd),
            "tool_stdout": stdout,
            "tool_stderr": stderr,
            "wall_time": wall_time,
        }

    return {
        "status": "success",
        "input": str(gene_path),
        "output_aa": str(out_aa),
        "output_nt": None,
        "n_taxa": n_taxa,
        "alignment_length": alignment_length,
        "wall_time": wall_time,
        "tool_cmd": " ".join(cmd),
        "tool_stdout": stdout,
        "tool_stderr": stderr,
        "warnings": warnings,
    }


def _align_one_worker(args: tuple[Path, Path, str, str, str | None, bool, str, str]) -> dict[str, Any]:
    gene_path, output_dir, method, seq_type, tool_args, dry_run, mafft_exe, magus_exe = args
    return _align_one(
        gene_path,
        output_dir,
        method=method,
        seq_type=seq_type,
        tool_args=tool_args,
        dry_run=dry_run,
        mafft_executable=mafft_exe,
        magus_executable=magus_exe,
    )


def _cleanup_magus_workdir(cmd: list[str]) -> None:
    try:
        idx = cmd.index("-d")
        work = Path(cmd[idx + 1])
        if work.exists():
            shutil.rmtree(work, ignore_errors=True)
    except (ValueError, IndexError):
        pass


def _run_backtrans_for_gene(
    res: dict[str, Any],
    *,
    nt_dir: Path,
    nt_out_dir: Path,
    trimal_exe: str,
) -> int:
    """Run backtrans for a single aligned gene, mutating ``res`` in place.

    Returns 1 if an NT alignment was produced, else 0. Records warnings on
    ``res['warnings']`` and sets ``res['output_nt']`` on success.
    """
    gene_stem = Path(res["input"]).stem
    nt_candidates = [
        p for p in nt_dir.iterdir()
        if p.is_file() and p.stem == gene_stem
    ]
    if not nt_candidates:
        res["warnings"].append(
            f"no matching CDS file found in --nt-dir for gene '{gene_stem}'"
        )
        return 0

    nt_path = nt_candidates[0]
    aa_aln_path = Path(res["output_aa"])

    try:
        aa_records = list(SeqIO.parse(str(aa_aln_path), "fasta"))
        aa_taxa = {r.id for r in aa_records}
        nt_records = list(SeqIO.parse(str(nt_path), "fasta"))
        nt_seqs = {r.id: str(r.seq) for r in nt_records}
    except Exception as exc:
        res["warnings"].append(f"could not read CDS file: {exc}")
        return 0

    cds_warnings = _validate_cds(nt_seqs, n_aa_taxa=res["n_taxa"], aa_taxa=aa_taxa)
    if cds_warnings:
        res["warnings"].extend(cds_warnings)
        res["warnings"].append("backtrans skipped due to CDS validation errors")
        return 0

    out_nt = nt_out_dir / f"{gene_stem}.fa"
    bt_result = _backtrans_one(aa_aln_path, nt_path, out_nt, dry_run=False, executable=trimal_exe)
    res["_bt_tool_result"] = {**bt_result, "input": str(nt_path)}
    if bt_result["status"] == "success":
        res["output_nt"] = bt_result["output_nt"]
        return 1
    res["warnings"].append(bt_result["reason"])
    return 0


def _backtrans_one(
    aa_aln_path: Path,
    nt_path: Path,
    output_nt_path: Path,
    dry_run: bool,
    executable: str = "trimal",
) -> dict[str, Any]:
    cmd = [
        executable,
        "-in", str(aa_aln_path),
        "-backtrans", str(nt_path),
        "-ignorestopcodon",
        "-out", str(output_nt_path),
        "-fasta",
    ]

    if dry_run:
        return {"status": "dry_run", "cmd": cmd}

    output_nt_path.parent.mkdir(parents=True, exist_ok=True)
    start = time.monotonic()
    try:
        proc = subprocess.run(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
        wall_time = time.monotonic() - start
    except Exception as exc:
        return {
            "status": "skipped",
            "reason": f"trimal error: {exc}",
            "tool_stdout": "",
            "tool_stderr": "",
            "wall_time": time.monotonic() - start,
        }

    if proc.returncode != 0:
        return {
            "status": "skipped",
            "reason": f"trimal -backtrans exited with code {proc.returncode}: {proc.stderr[:300]}",
            "tool_stdout": proc.stdout,
            "tool_stderr": proc.stderr,
            "wall_time": wall_time,
        }

    return {
        "status": "success",
        "output_nt": str(output_nt_path),
        "tool_cmd": " ".join(cmd),
        "tool_stdout": proc.stdout,
        "tool_stderr": proc.stderr,
        "wall_time": wall_time,
    }


def _resolve_tool_paths(
    method: str,
    backtrans: bool,
    mafft_path: Path | None,
    magus_path: Path | None,
    trimal_path: Path | None,
    dry_run: bool,
) -> tuple[str, str, str]:
    from phyloai.core.env import ToolEnv

    env = ToolEnv()
    if method == "magus":
        magus_exe = str(_validate_executable_path(magus_path, "magus")) if magus_path else (str(env.get("magus") or "magus") if dry_run else str(env.require("magus")))
        mafft_exe = str(mafft_path) if mafft_path else (str(env.get("mafft") or "mafft") if dry_run else str(env.require("mafft")))
    else:
        mafft_exe = str(_validate_executable_path(mafft_path, "mafft")) if mafft_path else (str(env.get("mafft") or "mafft") if dry_run else str(env.require("mafft")))
        magus_exe = str(magus_path) if magus_path else (str(env.get("magus") or "magus") if dry_run else str(env.require("magus")))

    if backtrans:
        trimal_exe = str(_validate_executable_path(trimal_path, "trimal")) if trimal_path else (str(env.get("trimal") or "trimal") if dry_run else str(env.require("trimal")))
    else:
        trimal_exe = str(trimal_path) if trimal_path else (str(env.get("trimal") or "trimal") if dry_run else str(env.require("trimal")))

    return mafft_exe, magus_exe, trimal_exe


def _validate_executable_path(path: Path, tool_name: str) -> Path:
    if not path.exists() or path.is_dir():
        raise FileNotFoundError(f"Required tool '{tool_name}' not found at explicit path: {path}")
    return path


def _detect_tool_versions(
    method: str,
    backtrans: bool,
    mafft_path: Path | None,
    magus_path: Path | None,
    trimal_path: Path | None,
) -> dict[str, str]:
    from phyloai.core.env import TOOL_REGISTRY, ToolEnv, ToolStatus

    tool_paths = {}
    if mafft_path:
        tool_paths["mafft"] = mafft_path
    if magus_path:
        tool_paths["magus"] = magus_path
    if trimal_path:
        tool_paths["trimal"] = trimal_path

    env = ToolEnv(tool_paths=tool_paths)
    versions: dict[str, str] = {}
    names = ["magus" if method == "magus" else "mafft"]
    if backtrans:
        names.append("trimal")

    for name in names:
        meta = TOOL_REGISTRY[name]
        info = env._detect_tool(
            name,
            version_flag=meta.get("version_flag", ""),
            version_args=meta.get("version_args"),
            bundled=False,
            path_aliases=meta.get("path_aliases"),
        )
        if info.status == ToolStatus.OK and info.version:
            versions[name] = info.version
    return versions


def run_align(
    seq_dir: Path,
    output_dir: Path,
    method: str,
    seq_type: str,
    backtrans: bool = False,
    nt_dir: Path | None = None,
    threads: int = 4,
    tool_args: str | None = None,
    mafft_path: Path | None = None,
    magus_path: Path | None = None,
    trimal_path: Path | None = None,
    overwrite: bool = False,
    resume: bool = False,
    dry_run: bool = False,
    quiet: bool = False,
    progress_callback: Callable[[Path], None] | None = None,
) -> dict[str, Any]:
    run_start = time.monotonic()

    from phyloai.core.checkpoint import (
        load_checkpoint,
        save_checkpoint_atomic,
        validate_resume_params,
    )
    from phyloai.pretree.checkpoint_helpers import (
        build_initial_checkpoint,
        mark_task,
        plan_resume,
    )

    if backtrans and nt_dir is None:
        raise ValueError("--nt-dir is required when --backtrans is set.")
    if backtrans and seq_type == "NT":
        raise ValueError("--backtrans requires --seq-type AA (backtrans produces NT from AA alignment).")
    if resume and overwrite:
        raise ValueError("--overwrite and --resume are mutually exclusive.")

    if method == "magus" and platform.system() != "Linux":
        raise ValueError(
            "--method magus requires Linux (bundled MAGUS binaries are Linux-only). "
            "On macOS, use a MAFFT method or run MAGUS manually with custom MAFFT/MCL paths."
        )

    global_warnings: list[str] = []
    if tool_args and method in MAFFT_METHODS:
        global_warnings.append(
            f"--tool-args is ignored for MAFFT method '{method}'; "
            "it is only used with --method magus."
        )
        tool_args = None

    found, scan_skipped = _scan_input(seq_dir)

    if not found and not dry_run:
        raise ValueError("No genes were aligned: no valid input files found.")

    if seq_type == "auto":
        if found:
            seq_type = _detect_seq_type_from_files(found)
        else:
            seq_type = "AA"
        global_warnings.append(f"seq_type auto-detected as '{seq_type}' from input files.")

    if backtrans and seq_type == "NT":
        raise ValueError("--backtrans requires --seq-type AA (backtrans produces NT from AA alignment).")

    mafft_exe, magus_exe, trimal_exe = _resolve_tool_paths(
        method=method,
        backtrans=backtrans,
        mafft_path=mafft_path,
        magus_path=magus_path,
        trimal_path=trimal_path,
        dry_run=dry_run,
    )

    _cmd_parts = [
        "phyloai", "pretree", "align",
        "--seq-dir", str(seq_dir),
        "--output-dir", str(output_dir),
        "--method", method,
        "--seq-type", seq_type,
        "--threads", str(threads),
    ]
    if backtrans:
        _cmd_parts.append("--backtrans")
    if nt_dir is not None:
        _cmd_parts += ["--nt-dir", str(nt_dir)]
    if tool_args is not None:
        _cmd_parts += ["--tool-args", tool_args]
    if mafft_path is not None:
        _cmd_parts += ["--mafft-path", str(mafft_path)]
    if magus_path is not None:
        _cmd_parts += ["--magus-path", str(magus_path)]
    if trimal_path is not None:
        _cmd_parts += ["--trimal-path", str(trimal_path)]
    if overwrite:
        _cmd_parts.append("--overwrite")
    if resume:
        _cmd_parts.append("--resume")
    if dry_run:
        _cmd_parts.append("--dry-run")
    if quiet:
        _cmd_parts.append("--quiet")
    full_command = " ".join(_cmd_parts)

    resolved = _resolved_align_params(
        seq_dir=seq_dir,
        output_dir=output_dir,
        method=method,
        resolved_seq_type=seq_type,
        backtrans=backtrans,
        nt_dir=nt_dir,
        threads=threads,
        tool_args=tool_args,
        mafft_path=mafft_exe,
        magus_path=magus_exe,
        trimal_path=trimal_exe,
        quiet=quiet,
    )

    if backtrans:
        aa_out_dir = output_dir / "seqs" / "faa"
        nt_out_dir = output_dir / "seqs" / "fna"
    else:
        aa_out_dir = output_dir / "seqs"
        nt_out_dir = None
    logs_dir = output_dir / "logs"

    skipped: list[dict[str, str]] = list(scan_skipped)
    checkpoint = None
    ckpt_path = output_dir / "checkpoint.json"
    to_run_ids: list[str] | None = None

    if resume:
        try:
            checkpoint = load_checkpoint(ckpt_path)
        except FileNotFoundError as exc:
            raise ValueError(str(exc)) from exc
        validate_resume_params(checkpoint, resolved, step="pretree.align")
        if checkpoint.status == "success":
            return reconstruct_align_result(
                checkpoint=checkpoint,
                params={**checkpoint.params, "overwrite": overwrite, "resume": resume, "dry_run": dry_run},
                tool_versions=_detect_tool_versions(
                    method=method,
                    backtrans=backtrans,
                    mafft_path=mafft_path,
                    magus_path=magus_path,
                    trimal_path=trimal_path,
                ),
                wall_time=0.0,
                skipped_inputs=list(scan_skipped),
                scan_warnings=list(global_warnings),
            )
        to_run_ids, _skipped_ids = plan_resume(checkpoint, verify_align_outputs)
        if not to_run_ids:
            checkpoint.status = "success"
            checkpoint.completed_at = checkpoint.touch()
            save_checkpoint_atomic(checkpoint, ckpt_path)
            return reconstruct_align_result(
                checkpoint=checkpoint,
                params={**checkpoint.params, "overwrite": overwrite, "resume": resume, "dry_run": dry_run},
                tool_versions=_detect_tool_versions(
                    method=method,
                    backtrans=backtrans,
                    mafft_path=mafft_path,
                    magus_path=magus_path,
                    trimal_path=trimal_path,
                ),
                wall_time=0.0,
                skipped_inputs=list(scan_skipped),
                scan_warnings=list(global_warnings),
            )
        found = [Path(task.input) for task in checkpoint.tasks if task.task_id in set(to_run_ids)]

    else:
        if not dry_run:
            if output_dir.exists() and any(output_dir.iterdir()):
                if not overwrite:
                    raise ValueError(
                        f"Output directory '{output_dir}' already exists and is non-empty. "
                        "Use --overwrite to replace it."
                    )
                shutil.rmtree(output_dir)
            output_dir.mkdir(parents=True, exist_ok=True)

    if not dry_run:
        aa_out_dir.mkdir(parents=True, exist_ok=True)
        if nt_out_dir is not None:
            nt_out_dir.mkdir(parents=True, exist_ok=True)
        logs_dir.mkdir(parents=True, exist_ok=True)

    if not resume and not dry_run:
        checkpoint = build_initial_checkpoint(
            step="pretree.align",
            command=full_command,
            params=resolved,
            inputs=found,
            output_for=lambda p: aa_out_dir / f"{p.stem}.fa",
            nt_output_for=(lambda p: None) if nt_out_dir is None else (lambda p: nt_out_dir / f"{p.stem}.fa"),
        )
        save_checkpoint_atomic(checkpoint, ckpt_path)
        to_run_ids = [path.stem for path in found]

    _ckpt_write = (
        checkpoint is not None and ckpt_path is not None and to_run_ids and not dry_run
    )
    _to_run_set = set(to_run_ids) if to_run_ids else set()

    if _ckpt_write:
        for task_id in to_run_ids:
            mark_task(checkpoint, task_id, status="running", reason=None)
        save_checkpoint_atomic(checkpoint, ckpt_path)

    file_results: list[dict[str, Any]] = []

    worker_args = [
        (g, aa_out_dir, method, seq_type, tool_args, dry_run, mafft_exe, magus_exe)
        for g in found
    ]
    n_backtrans = 0
    do_backtrans = bool(backtrans and nt_dir and not dry_run)

    # Throttled checkpoint flushing: writing the full checkpoint after every task
    # is O(N^2) and dominated by fsync on large files. Instead flush at most once
    # every CHECKPOINT_FLUSH_INTERVAL seconds (no fsync), force a flush on
    # interrupt, and do one durable fsync flush at the very end.
    _last_flush = time.monotonic()

    def _maybe_flush(*, force: bool = False) -> None:
        nonlocal _last_flush
        if not _ckpt_write:
            return
        now = time.monotonic()
        if force or (now - _last_flush) >= CHECKPOINT_FLUSH_INTERVAL:
            save_checkpoint_atomic(checkpoint, ckpt_path)
            _last_flush = now

    interrupted = False
    try:
        with ProcessPoolExecutor(max_workers=threads) as pool:
            futures = {pool.submit(_align_one_worker, arg): arg[0] for arg in worker_args}
            for future in as_completed(futures):
                gene_path = futures[future]
                try:
                    result = future.result()
                except Exception as exc:
                    result = {
                        "status": "skipped",
                        "input": str(gene_path),
                        "reason": str(exc),
                    }

                if result["status"] == "skipped":
                    skipped.append({"path": result["input"], "reason": result.get("reason", "unknown")})
                    if _ckpt_write:
                        task_id = Path(result["input"]).stem
                        if task_id in _to_run_set:
                            mark_task(checkpoint, task_id, status="failed", reason=result.get("reason"))
                    if progress_callback:
                        progress_callback(gene_path)
                    _maybe_flush()
                    continue

                # AA alignment succeeded. Run this gene's backtrans inline so each
                # gene is a single complete, resumable unit (AA + NT together).
                file_results.append(result)
                if do_backtrans:
                    n_backtrans += _run_backtrans_for_gene(
                        result, nt_dir=nt_dir, nt_out_dir=nt_out_dir, trimal_exe=trimal_exe
                    )
                    result.pop("_bt_tool_result", None)

                if not dry_run:
                    log_path = logs_dir / f"{gene_path.stem}.log"
                    _out = result.get("tool_stdout", "").strip()
                    _err = result.get("tool_stderr", "").strip()
                    if _out and _err:
                        log_path.write_text(f"{_out}\n{_err}")
                    else:
                        log_path.write_text(_out or _err)

                if _ckpt_write:
                    task_id = Path(result["input"]).stem
                    if task_id in _to_run_set:
                        mark_task(checkpoint, task_id, status="success", reason=None)
                if progress_callback:
                    progress_callback(gene_path)
                _maybe_flush()
    except KeyboardInterrupt:
        interrupted = True

    if _ckpt_write:
        if interrupted:
            checkpoint.status = "interrupted"
        save_checkpoint_atomic(checkpoint, ckpt_path, fsync=True)

    if interrupted:
        raise KeyboardInterrupt

    if not dry_run and not file_results:
        if checkpoint is not None:
            checkpoint.status = "error"
            checkpoint.touch()
            save_checkpoint_atomic(checkpoint, ckpt_path)
        raise ValueError("No genes were aligned: all input files failed or were skipped.")

    if checkpoint is not None and not dry_run:
        checkpoint.status = "success"
        checkpoint.completed_at = checkpoint.touch()
        save_checkpoint_atomic(checkpoint, ckpt_path, fsync=True)
        return reconstruct_align_result(
            checkpoint=checkpoint,
            params={**resolved, "overwrite": overwrite, "resume": resume, "dry_run": dry_run},
            tool_versions=_detect_tool_versions(
                method=method,
                backtrans=backtrans,
                mafft_path=mafft_path,
                magus_path=magus_path,
                trimal_path=trimal_path,
            ),
            wall_time=time.monotonic() - run_start,
            skipped_inputs=skipped,
            scan_warnings=list(global_warnings),
            file_results=file_results,
        )

    n_aligned = len(file_results)
    aligned_lengths = [r["alignment_length"] for r in file_results if r.get("alignment_length")]
    aligned_taxa = [r["n_taxa"] for r in file_results if r.get("n_taxa")]
    mean_len = round(sum(aligned_lengths) / len(aligned_lengths), 1) if aligned_lengths else 0.0
    mean_taxa = round(sum(aligned_taxa) / len(aligned_taxa), 1) if aligned_taxa else 0.0

    all_warnings = list(global_warnings)
    for r in file_results:
        all_warnings.extend(r.get("warnings", []))

    payload: dict[str, Any] = {
        "status": "success",
        "command": full_command,
        "wall_time": time.monotonic() - run_start,
        "tool_versions": _detect_tool_versions(
            method=method,
            backtrans=backtrans,
            mafft_path=mafft_path,
            magus_path=magus_path,
            trimal_path=trimal_path,
        ),
        "params": {
            "seq_dir": str(seq_dir),
            "method": method,
            "seq_type": seq_type,
            "backtrans": backtrans,
            "nt_dir": str(nt_dir) if nt_dir else None,
            "output_dir": str(output_dir),
            "threads": threads,
            "tool_args": tool_args,
            "mafft_path": str(mafft_path) if mafft_path else None,
            "magus_path": str(magus_path) if magus_path else None,
            "trimal_path": str(trimal_path) if trimal_path else None,
            "overwrite": overwrite,
            "resume": resume,
            "dry_run": dry_run,
            "quiet": quiet,
        },
        "key_results": {
            "n_aligned": n_aligned,
            "n_skipped": len(skipped),
            "method": method,
            "backtrans": backtrans,
            "mean_alignment_length": mean_len,
            "mean_n_taxa": mean_taxa,
        },
        "error": None,
        "data": {
            "summary": {
                "n_input_files": len(found) + len(scan_skipped),
                "n_aligned": n_aligned,
                "n_backtrans": n_backtrans,
                "n_skipped": len(skipped),
            },
            "files": [
                {
                    "input": r["input"],
                    "output_aa": r.get("output_aa"),
                    "output_nt": r.get("output_nt"),
                    "cmd": r.get("cmd") if isinstance(r.get("cmd"), list) else shlex.split(r.get("tool_cmd", "")),
                    "log_file": f"logs/{Path(r['input']).stem}.log",
                    "n_taxa": r.get("n_taxa", 0),
                    "alignment_length": r.get("alignment_length", 0),
                    "wall_time": r.get("wall_time", 0.0),
                    "warnings": r.get("warnings", []),
                }
                for r in file_results
            ],
            "skipped": skipped,
            "warnings": all_warnings,
        },
    }
    return payload


def render_align_summary_table(summary: dict[str, Any]) -> "Table":
    from rich.table import Table as _Table
    table = _Table(title="pretree align summary")
    table.add_column("Metric")
    table.add_column("Value")
    for key in ["n_input_files", "n_aligned", "n_backtrans", "n_skipped"]:
        table.add_row(key, str(summary.get(key, "")))
    return table
