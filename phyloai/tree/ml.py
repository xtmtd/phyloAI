"""Maximum-likelihood tree inference with FastTree."""

from __future__ import annotations

import json
import os
import shlex
import shutil
import subprocess
import time
import time as _time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Any, Callable

from Bio import SeqIO

from phyloai.core.checkpoint import load_checkpoint, save_checkpoint_atomic, validate_resume_params
from phyloai.core.env import ToolEnv
from phyloai.core.schema import COMMON_ALIGNMENT_EXTENSIONS
from phyloai.core.sequence_normalization import detect_seq_type
from phyloai.tree.checkpoint_helpers import build_initial_checkpoint, mark_task, plan_resume

FASTTREE_COMPATIBLE_EXTENSIONS = frozenset({
    ".fa", ".fas", ".fasta", ".faa", ".fna",
    ".phy", ".phylip",
})

FASTTREE_MANAGED_FLAGS = frozenset({
    "-nt", "-expert", "-help",
})

CHECKPOINT_FLUSH_INTERVAL = 2.0


def _scan_input(msa_dir: Path) -> tuple[list[Path], list[dict[str, str]]]:
    if not msa_dir.exists():
        return [], []

    found: list[Path] = []
    skipped: list[dict[str, str]] = []

    for entry in sorted(msa_dir.iterdir()):
        if entry.is_dir():
            skipped.append({"path": str(entry), "reason": "directory"})
            continue
        if not entry.is_file():
            skipped.append({"path": str(entry), "reason": "not a regular file"})
            continue
        if entry.stat().st_size == 0:
            skipped.append({"path": str(entry), "reason": "empty file"})
            continue

        ext = entry.suffix.lower()
        if ext in FASTTREE_COMPATIBLE_EXTENSIONS:
            found.append(entry)
        elif ext in {".nex", ".nxs", ".nexus"}:
            skipped.append({
                "path": str(entry),
                "reason": "NEXUS format not supported by FastTree; use pretree convert first",
            })
        elif ext in set(COMMON_ALIGNMENT_EXTENSIONS):
            skipped.append({"path": str(entry), "reason": f"unrecognized extension: {ext}"})
        else:
            skipped.append({"path": str(entry), "reason": f"unrecognized extension: {ext}"})

    return found, skipped


def _build_fasttree_cmd(
    input_path: Path,
    output_path: Path,
    *,
    executable: str = "FastTree",
    seq_type: str = "AA",
    model: str = "lg",
    mode: str = "normal",
    boot: int = 1000,
    cat: int = 20,
    gamma: bool = True,
    tool_args: str | None = None,
) -> list[str]:
    cmd = [executable]

    tool_tokens = set(shlex.split(tool_args)) if tool_args else set()
    has_model_override = bool({"-lg", "-wag", "-gtr"} & tool_tokens)
    has_mode_override = bool({"-fastest", "-slow"} & tool_tokens)
    has_gamma_override = "-gamma" in tool_tokens
    has_cat_override = "-cat" in tool_tokens
    has_boot_override = bool({"-boot", "-nosupport"} & tool_tokens)

    if seq_type == "NT":
        cmd.append("-nt")
        if model == "gtr" and not has_model_override:
            cmd.append("-gtr")
    else:
        if not has_model_override:
            if model == "lg":
                cmd.append("-lg")
            elif model == "wag":
                cmd.append("-wag")

    if not has_mode_override:
        if mode == "fastest":
            cmd.append("-fastest")
        elif mode == "slow":
            cmd.append("-slow")

    if gamma and not has_gamma_override:
        cmd.append("-gamma")

    if not has_cat_override:
        cmd.extend(["-cat", str(cat)])

    if not has_boot_override:
        if boot > 0:
            cmd.extend(["-boot", str(boot)])
        else:
            cmd.append("-nosupport")

    if tool_args:
        _check_managed_flag_conflict(tool_args)
        cmd.extend(shlex.split(tool_args))

    cmd.append(str(input_path))
    return cmd


def _check_managed_flag_conflict(tool_args: str) -> None:
    tokens = shlex.split(tool_args)
    managed_set = FASTTREE_MANAGED_FLAGS
    for token in tokens:
        if token in managed_set:
            raise ValueError(f"Blocked managed flag in --tool-args: {token}")
        if "/" in token or ">" in token:
            raise ValueError(f"Blocked I/O override in --tool-args: {token}")


def _run_one_fasttree(
    gene_path: Path,
    *,
    seq_type: str,
    model: str,
    mode: str,
    boot: int,
    cat: int,
    gamma: bool,
    tool_args: str | None,
    log_dir: Path,
    fasttree_executable: str = "FastTree",
    output_dir: Path | None = None,
    dry_run: bool = False,
) -> dict[str, Any]:
    if output_dir is None:
        output_dir = gene_path.parent

    result: dict[str, Any] = {
        "input": str(gene_path),
        "output_tree": None,
        "log_file": None,
        "n_taxa": 0,
    }

    if not gene_path.exists():
        return {
            **result,
            "status": "failed",
            "reason": f"input file not found: {gene_path}",
            "wall_time": 0,
            "warnings": [],
        }

    stem = gene_path.stem
    out_tree = output_dir / f"{stem}.tre"
    out_log = log_dir / f"{stem}.log"

    try:
        n_taxa = sum(1 for _ in SeqIO.parse(str(gene_path), "fasta"))
        result["n_taxa"] = n_taxa
    except Exception:
        pass

    cmd = _build_fasttree_cmd(
        gene_path, out_tree,
        executable=fasttree_executable,
        seq_type=seq_type, model=model, mode=mode,
        boot=boot, cat=cat, gamma=gamma,
        tool_args=tool_args,
    )

    result.update({
        "output_tree": str(out_tree),
        "log_file": str(out_log),
        "cmd": cmd,
    })

    if dry_run:
        return {**result, "status": "dry_run", "wall_time": 0, "warnings": []}

    warnings: list[str] = []
    start = _time.monotonic()
    try:
        out_tree.parent.mkdir(parents=True, exist_ok=True)
        out_log.parent.mkdir(parents=True, exist_ok=True)

        proc = subprocess.run(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
        wall_time = _time.monotonic() - start

        out_tree.write_text(proc.stdout)
        out_log.write_text(proc.stderr)
        result["log_file"] = str(out_log)

        if proc.returncode != 0:
            return {
                **result,
                "status": "failed",
                "reason": f"FastTree exited with code {proc.returncode}: {proc.stderr[:200]}",
                "wall_time": wall_time,
                "warnings": warnings,
            }

        from Bio import Phylo
        try:
            Phylo.read(str(out_tree), "newick")
        except Exception as e:
            return {
                **result,
                "status": "failed",
                "reason": f"FastTree produced unparseable Newick output: {e}",
                "wall_time": wall_time,
                "warnings": warnings,
            }

        return {
            **result,
            "status": "success",
            "wall_time": wall_time,
            "tool_stderr": proc.stderr,
            "warnings": warnings,
        }

    except Exception as exc:
        return {
            **result,
            "status": "failed",
            "reason": str(exc),
            "wall_time": _time.monotonic() - start,
            "warnings": warnings,
        }


def _validate_seq_types(
    files: list[Path],
    *,
    declared_type: str | None,
) -> tuple[str | None, list[dict[str, Any]]]:
    if not files:
        return (declared_type or "AA"), []

    all_types: dict[str, str] = {}
    offending: list[dict[str, Any]] = []

    for f in files:
        try:
            ext = f.suffix.lower()
            if ext in {".phy", ".phylip"}:
                seqs = [str(r.seq) for r in SeqIO.parse(str(f), "phylip-relaxed")]
            else:
                seqs = [str(r.seq) for r in SeqIO.parse(str(f), "fasta")]
            if not seqs:
                offending.append({"file": str(f), "reason": "no sequences found"})
                continue
            dt = detect_seq_type(seqs)
            all_types[str(f)] = dt
        except Exception:
            offending.append({"file": str(f), "reason": "failed to parse input file"})
            continue

    if declared_type:
        for f_str, dt in all_types.items():
            if dt != declared_type:
                offending.append({"file": f_str, "expected": declared_type, "detected": dt})
        return declared_type, offending

    type_counts: dict[str, int] = {}
    for dt in all_types.values():
        type_counts[dt] = type_counts.get(dt, 0) + 1

    if len(type_counts) == 1:
        resolved = next(iter(type_counts))
        return resolved, []

    majority = max(type_counts, key=type_counts.get)
    for f_str, dt in all_types.items():
        if dt != majority:
            offending.append({"file": f_str, "expected": majority, "detected": dt})

    return None, offending


def _format_offender(o: dict[str, Any]) -> str:
    if "reason" in o:
        return f"{o['file']}: {o['reason']}"
    return f"{o['file']}: {o['detected']} (expected {o['expected']})"


def _resolved_fasttree_params(
    msa_dir: Path | None,
    matrix: Path | None,
    seq_type: str,
    model: str,
    mode: str,
    boot: int,
    cat: int,
    gamma: bool,
    output_dir: Path,
    threads: int,
    fasttree_path: str | None,
    tool_args: str | None,
) -> dict[str, Any]:
    return {
        "msa_dir": str(msa_dir) if msa_dir else None,
        "matrix": str(matrix) if matrix else None,
        "seq_type": seq_type,
        "model": model,
        "mode": mode,
        "boot": boot,
        "cat": cat,
        "gamma": gamma,
        "output_dir": str(output_dir),
        "threads": threads,
        "fasttree_path": fasttree_path,
        "tool_args": tool_args,
    }


def run_fasttree(
    *,
    msa_dir: Path | None = None,
    matrix: Path | None = None,
    output_dir: Path,
    seq_type: str = "auto",
    model: str | None = None,
    mode: str = "normal",
    boot: int = 1000,
    cat: int = 20,
    gamma: bool = True,
    threads: int = 4,
    fasttree_path: str | None = None,
    tool_args: str | None = None,
    overwrite: bool = False,
    resume: bool = False,
    dry_run: bool = False,
    quiet: bool = False,
    progress_callback: Callable[[Path], None] | None = None,
) -> dict[str, Any]:
    run_start = time.monotonic()
    output_dir = output_dir.resolve()

    if (msa_dir is None and matrix is None) or (msa_dir is not None and matrix is not None):
        raise ValueError("Either --msa-dir or --matrix must be provided (not both).")

    fasttree_exe = _resolve_fasttree(fasttree_path, dry_run)

    batch_mode = msa_dir is not None
    n_resume_skipped = 0

    trees_dir = output_dir / "trees"
    logs_dir = output_dir / "logs"

    resolved_seq_type = seq_type
    if not batch_mode:
        assert matrix is not None
        matrix_ext = matrix.suffix.lower()
        if matrix_ext not in FASTTREE_COMPATIBLE_EXTENSIONS:
            raise ValueError(
                f"--matrix has unsupported extension: {matrix.suffix}. "
                f"Supported: {sorted(FASTTREE_COMPATIBLE_EXTENSIONS)}"
            )
        matrix_format = "phylip-relaxed" if matrix_ext in {".phy", ".phylip"} else "fasta"
        try:
            recs = list(SeqIO.parse(str(matrix), matrix_format))
        except Exception:
            recs = []
        if not recs:
            raise ValueError(
                f"Cannot parse --matrix as {matrix_format}: {matrix}"
            )
        if seq_type == "auto":
            resolved_seq_type = detect_seq_type([str(r.seq) for r in recs])
        else:
            sample = [str(r.seq) for r in recs[:10]]
            if sample:
                resolved_sample = detect_seq_type(sample)
                if resolved_sample != seq_type:
                    raise ValueError(
                        f"--seq-type {seq_type} but detected {resolved_sample} in {matrix}"
                    )

    found: list[Path] = []
    skipped_input: list[dict[str, str]] = []
    if batch_mode:
        assert msa_dir is not None
        found, skipped_input = _scan_input(msa_dir)
        if not found:
            raise ValueError("No valid input files found in --msa-dir")
        stems = [p.stem for p in found]
        dupes = {s for s in stems if stems.count(s) > 1}
        if dupes:
            colliding = sorted(p.name for p in found if p.stem in dupes)
            raise ValueError(
                "Duplicate output stems detected in --msa-dir. Files with different extensions "
                "but the same base name collide:\n  " + "\n  ".join(colliding) +
                "\nRename files so each stem is unique within the directory."
            )
        declared = None if seq_type == "auto" else seq_type
        resolved_seq_type, offending = _validate_seq_types(found, declared_type=declared)
        if resolved_seq_type is None:
            offending_strs = [_format_offender(o) for o in offending[:10]]
            raise ValueError("Mixed sequence types in --msa-dir:\n" + "\n".join(offending_strs))
        if offending:
            offending_strs = [_format_offender(o) for o in offending[:10]]
            raise ValueError(
                f"Files with wrong --seq-type ({declared}) in --msa-dir:\n" + "\n".join(offending_strs)
            )

    if model is None:
        model = "gtr" if resolved_seq_type == "NT" else "lg"

    aa_models = {"jtt", "lg", "wag"}
    nt_models = {"jc", "gtr"}
    if resolved_seq_type == "AA":
        if model not in aa_models:
            raise ValueError(f"Invalid model for AA: {model}. Choose from {aa_models}")
    elif resolved_seq_type == "NT":
        if model not in nt_models:
            raise ValueError(f"Invalid model for NT: {model}. Choose from {nt_models}")

    checkpoint: Any = None
    ckpt_path = output_dir / "checkpoint.json"

    if not dry_run:
        if overwrite and resume:
            raise ValueError("--overwrite and --resume are mutually exclusive")
        if resume:
            if not ckpt_path.exists():
                raise ValueError(f"--resume requires {ckpt_path}, not found")
            checkpoint = load_checkpoint(ckpt_path)
            resolved_params = _resolved_fasttree_params(
                msa_dir=msa_dir, matrix=matrix,
                seq_type=resolved_seq_type, model=model,
                mode=mode, boot=boot, cat=cat, gamma=gamma,
                output_dir=output_dir, threads=threads,
                fasttree_path=fasttree_path, tool_args=tool_args,
            )
            validate_resume_params(checkpoint, resolved_params, step="tree.ml.fasttree")
            if checkpoint.status == "success":
                return _reconstruct_result(output_dir, run_start)

            to_run_ids, skipped_ids = plan_resume(checkpoint)
            n_resume_skipped = len(skipped_ids)
            if not to_run_ids:
                checkpoint.status = "success"
                save_checkpoint_atomic(checkpoint, ckpt_path)
                return _reconstruct_result(output_dir, run_start)
            found = [Path(task.input) for task in checkpoint.tasks if task.task_id in to_run_ids]
        else:
            if overwrite and output_dir.exists():
                shutil.rmtree(output_dir)
            if output_dir.exists() and any(output_dir.iterdir()):
                raise ValueError(
                    f"Output directory {output_dir} already exists and is non-empty. "
                    "Use --overwrite to replace."
                )
            if batch_mode:
                trees_dir.mkdir(parents=True, exist_ok=True)
                logs_dir.mkdir(parents=True, exist_ok=True)
            else:
                output_dir.mkdir(parents=True, exist_ok=True)

    if not batch_mode:
        assert matrix is not None
        result = _run_one_fasttree(
            gene_path=matrix,
            seq_type=resolved_seq_type,
            model=model,
            mode=mode,
            boot=boot,
            cat=cat,
            gamma=gamma,
            tool_args=tool_args,
            log_dir=output_dir,
            fasttree_executable=fasttree_exe,
            output_dir=output_dir,
            dry_run=dry_run,
        )
        return _assemble_result(
            run_start=run_start, fasttree_exe=fasttree_exe,
            batch_mode=False, results=[result],
            resolved_seq_type=resolved_seq_type, model=model, mode=mode, boot=boot,
            cat=cat, gamma=gamma, output_dir=output_dir,
            msa_dir=msa_dir, matrix=matrix,
            fasttree_path=fasttree_path, tool_args=tool_args,
            overwrite=overwrite, threads=threads,
            skipped_input=[],
            resume=resume,
            dry_run=dry_run,
            quiet=quiet,
        )

    if not resume and not dry_run:
        resolved_params = _resolved_fasttree_params(
            msa_dir=msa_dir, matrix=matrix,
            seq_type=resolved_seq_type, model=model,
            mode=mode, boot=boot, cat=cat, gamma=gamma,
            output_dir=output_dir, threads=threads,
            fasttree_path=fasttree_path, tool_args=tool_args,
        )
        cmd_parts = ["phyloai", "tree", "ml", "fasttree", "--msa-dir", str(msa_dir)]
        cmd_parts.extend([
            "--seq-type", resolved_seq_type, "--model", model,
            "--mode", mode, "--boot", str(boot), "--cat", str(cat),
        ])
        if not gamma:
            cmd_parts.append("--no-gamma")
        cmd_parts.extend(["-o", str(output_dir)])
        if threads != 4:
            cmd_parts.extend(["-t", str(threads)])
        if fasttree_path:
            cmd_parts.extend(["--fasttree-path", fasttree_path])
        if tool_args:
            if " " in tool_args:
                cmd_parts.append(f"--tool-args '{tool_args}'")
            else:
                cmd_parts.extend(["--tool-args", tool_args])
        if overwrite:
            cmd_parts.append("--overwrite")
        checkpoint = build_initial_checkpoint(
            step="tree.ml.fasttree",
            command=" ".join(cmd_parts),
            params=resolved_params,
            inputs=found,
            trees_dir=trees_dir,
            logs_dir=logs_dir,
        )
        save_checkpoint_atomic(checkpoint, ckpt_path)

    _ckpt_write = checkpoint is not None and not dry_run
    _last_flush = time.monotonic()

    def _maybe_flush(*, force: bool = False) -> None:
        nonlocal _last_flush
        if not _ckpt_write:
            return
        now = time.monotonic()
        if force or (now - _last_flush) >= CHECKPOINT_FLUSH_INTERVAL:
            save_checkpoint_atomic(checkpoint, ckpt_path)
            _last_flush = now

    file_results: list[dict[str, Any]] = []
    failed_results: list[dict[str, Any]] = []

    worker_args = [
        (p, resolved_seq_type, model, mode, boot, cat, gamma,
         tool_args, logs_dir, fasttree_exe, trees_dir, dry_run)
        for p in found
    ]

    interrupted = False
    try:
        if dry_run:
            for arg in worker_args:
                result = _run_one_fasttree(
                    gene_path=arg[0], seq_type=arg[1], model=arg[2], mode=arg[3],
                    boot=arg[4], cat=arg[5], gamma=arg[6], tool_args=arg[7],
                    log_dir=arg[8], fasttree_executable=arg[9], output_dir=arg[10],
                    dry_run=arg[11],
                )
                file_results.append(result)
                if progress_callback:
                    progress_callback(arg[0])
        else:
            with ProcessPoolExecutor(max_workers=threads) as pool:
                futures = {
                    pool.submit(_run_one_fasttree,
                        gene_path=arg[0], seq_type=arg[1], model=arg[2], mode=arg[3],
                        boot=arg[4], cat=arg[5], gamma=arg[6], tool_args=arg[7],
                        log_dir=arg[8], fasttree_executable=arg[9], output_dir=arg[10],
                        dry_run=arg[11],
                    ): arg[0]
                    for arg in worker_args
                }
                for future in as_completed(futures):
                    gene_path = futures[future]
                    result = future.result()
                    task_id = gene_path.stem

                    if result["status"] == "success":
                        file_results.append(result)
                        mark_task(checkpoint, task_id, status="success")
                    elif result["status"] == "failed":
                        failed_results.append(result)
                        mark_task(checkpoint, task_id, status="failed",
                                  reason=result.get("reason"))
                    else:
                        skipped_input.append({
                            "path": result.get("input", ""),
                            "reason": result.get("reason", "unknown"),
                        })

                    if progress_callback:
                        progress_callback(gene_path)
                    _maybe_flush()

    except KeyboardInterrupt:
        interrupted = True

    from datetime import datetime as _dt_cls, timezone as _tz
    if _ckpt_write:
        if interrupted:
            checkpoint.status = "interrupted"
        elif failed_results:
            checkpoint.status = "running"
        else:
            checkpoint.status = "success"
            checkpoint.completed_at = _dt_cls.now(_tz.utc).isoformat(timespec="seconds")
        save_checkpoint_atomic(checkpoint, ckpt_path, fsync=True)
    if interrupted:
        raise KeyboardInterrupt

    return _assemble_result(
        run_start=run_start, fasttree_exe=fasttree_exe,
        batch_mode=True, results=file_results,
        failed_results=failed_results,
        resolved_seq_type=resolved_seq_type, model=model, mode=mode, boot=boot,
        cat=cat, gamma=gamma, output_dir=output_dir,
        msa_dir=msa_dir, matrix=matrix,
        fasttree_path=fasttree_path, tool_args=tool_args,
        overwrite=overwrite, threads=threads,
        skipped_input=skipped_input,
        n_resume_skipped=n_resume_skipped,
        resume=resume,
        dry_run=dry_run,
        quiet=quiet,
    )


def _resolve_fasttree(fasttree_path: str | None, dry_run: bool) -> str:
    if fasttree_path:
        p = Path(fasttree_path)
        if not p.exists():
            raise ValueError(f"--fasttree-path does not exist: {fasttree_path}")
        if not os.access(p, os.X_OK):
            raise ValueError(f"--fasttree-path is not executable: {fasttree_path}")
        return fasttree_path
    if dry_run:
        return "FastTree"
    try:
        env = ToolEnv()
        return str(env.require("FastTree"))
    except FileNotFoundError:
        raise FileNotFoundError("FastTree not found. Install it or use --fasttree-path.")


def _assemble_result(
    *,
    run_start: float,
    fasttree_exe: str,
    batch_mode: bool,
    results: list[dict[str, Any]],
    failed_results: list[dict[str, Any]] | None = None,
    resolved_seq_type: str,
    model: str,
    mode: str,
    boot: int,
    cat: int,
    gamma: bool,
    output_dir: Path,
    msa_dir: Path | None,
    matrix: Path | None,
    fasttree_path: str | None,
    tool_args: str | None,
    overwrite: bool,
    threads: int,
    skipped_input: list[dict[str, str]],
    n_resume_skipped: int = 0,
    resume: bool = False,
    dry_run: bool = False,
    quiet: bool = False,
) -> dict[str, Any]:
    if failed_results is None:
        failed_results = []

    all_ok = [r for r in results if r["status"] in {"success", "dry_run"}]
    n_trees = len(all_ok) + n_resume_skipped
    n_failed = len(failed_results)
    n_skipped = len(skipped_input)

    is_error = n_trees == 0 and (n_failed > 0 or n_skipped > 0)
    if is_error:
        error_msg = "All FastTree runs failed"
    else:
        error_msg = None

    mean_n_taxa = 0.0
    mean_wall_time = 0.0
    if n_trees > 0:
        total_n_taxa = sum(r.get("n_taxa", 0) for r in all_ok)
        mean_n_taxa = total_n_taxa / n_trees if n_trees else 0.0
        total_wall = sum(r.get("wall_time", 0.0) for r in all_ok)
        mean_wall_time = total_wall / n_trees if n_trees else 0.0

    try:
        versions = _detect_fasttree_version(fasttree_exe)
    except Exception:
        versions = {"FastTree": "unknown"}

    cmd_parts = ["phyloai", "tree", "ml", "fasttree"]
    if batch_mode:
        cmd_parts.extend(["--msa-dir", str(msa_dir)])
    else:
        cmd_parts.extend(["--matrix", str(matrix)])
    cmd_parts.extend([
        "--seq-type", resolved_seq_type, "--model", model,
        "--mode", mode, "--boot", str(boot), "--cat", str(cat),
    ])
    if not gamma:
        cmd_parts.append("--no-gamma")
    cmd_parts.extend(["-o", str(output_dir)])
    if threads != 4:
        cmd_parts.extend(["-t", str(threads)])
    if fasttree_path:
        cmd_parts.extend(["--fasttree-path", fasttree_path])
    if tool_args:
        if " " in tool_args:
            cmd_parts.append(f"--tool-args '{tool_args}'")
        else:
            cmd_parts.extend(["--tool-args", tool_args])
    if overwrite:
        cmd_parts.append("--overwrite")
    if dry_run:
        cmd_parts.append("--dry-run")
    if resume:
        cmd_parts.append("--resume")
    cmd_str = " ".join(cmd_parts)

    if batch_mode:
        ok_files = []
        for r in all_ok:
            entry = {k: v for k, v in r.items() if k != "tool_stderr"}
            log_path = Path(entry.get("log_file", ""))
            if "log_file" in entry and log_path.is_absolute():
                try:
                    entry["log_file"] = str(log_path.relative_to(output_dir))
                except ValueError:
                    pass
            ok_files.append(entry)
        ok_failed = [
            {k: v for k, v in r.items() if k != "tool_stderr"}
            for r in failed_results
        ]
        for entry in ok_failed:
            log_path = Path(entry.get("log_file", ""))
            if "log_file" in entry and log_path.is_absolute():
                try:
                    entry["log_file"] = str(log_path.relative_to(output_dir))
                except ValueError:
                    pass
        data_block: dict[str, Any] = {
            "summary": {
                "n_input_files": len(results) + n_failed + n_skipped + n_resume_skipped,
                "n_trees": n_trees,
                "n_failed": n_failed,
                "n_skipped": n_skipped,
                "n_resume_skipped": n_resume_skipped,
                "mean_n_taxa": mean_n_taxa,
                "mean_wall_time": mean_wall_time,
                "mode": "--msa-dir" if batch_mode else "--matrix",
            },
            "files": ok_files,
            "failed": ok_failed,
            "skipped": skipped_input,
            "warnings": [],
        }
    else:
        first = results[0] if results else {}
        data_block = {
            "cmd": first.get("cmd", []),
            "tool_stderr": "",
            "output": first.get("output_tree", ""),
            "warnings": first.get("warnings", []),
        }

    payload: dict[str, Any] = {
        "status": "error" if is_error else "success",
        "command": cmd_str,
        "wall_time": time.monotonic() - run_start,
        "tool_versions": versions,
        "params": {
            "msa_dir": str(msa_dir) if msa_dir else None,
            "matrix": str(matrix) if matrix else None,
            "seq_type": resolved_seq_type,
            "model": model,
            "mode": mode,
            "boot": boot,
            "cat": cat,
            "gamma": gamma,
            "output_dir": str(output_dir),
            "threads": threads,
            "overwrite": overwrite,
            "resume": resume,
            "dry_run": dry_run,
            "quiet": quiet,
            "fasttree_path": fasttree_path,
            "tool_args": tool_args,
        },
        "key_results": {
            "n_input": len(results) + n_failed + n_skipped + n_resume_skipped,
            "n_trees": n_trees,
            "n_failed": n_failed,
            "n_skipped": n_skipped,
            "seq_type": resolved_seq_type,
            "model": model,
            "mode": mode,
            "boot": boot,
        },
        "error": error_msg,
        "data": data_block,
    }

    return payload


def _detect_fasttree_version(executable: str) -> dict[str, str]:
    import subprocess as _sp
    import re as _re

    exe_name = Path(executable).name
    try:
        proc = _sp.run([executable], capture_output=True, text=True, timeout=10)
        combined = proc.stdout + proc.stderr
    except Exception:
        return {exe_name: "unknown"}

    m = _re.search(r"(?:version|FastTree)\s*([\d.]+)", combined, _re.IGNORECASE)
    if m:
        return {exe_name: m.group(1)}

    m = _re.search(r"([\d]+\.[\d]+(?:\.[\d]+)?)", combined)
    if m:
        return {exe_name: m.group(1)}

    return {exe_name: "unknown"}


def _reconstruct_result(output_dir: Path, run_start: float) -> dict[str, Any]:
    result_path = output_dir / "result.json"
    if result_path.exists():
        return json.loads(result_path.read_text())
    ckpt_path = output_dir / "checkpoint.json"
    if ckpt_path.exists():
        ckpt = json.loads(ckpt_path.read_text())
        cmd = ckpt.get("command", "")
    else:
        cmd = ""
    if not cmd:
        return {
            "status": "error",
            "command": "phyloai tree ml fasttree",
            "wall_time": time.monotonic() - run_start,
            "tool_versions": {},
            "params": {},
            "key_results": {},
            "error": "Cannot reconstruct result: result.json not found and checkpoint.json missing or has no command",
            "data": {"summary": {}, "files": [], "failed": [], "skipped": [], "warnings": []},
        }
    return {
        "status": "success",
        "command": cmd,
        "wall_time": time.monotonic() - run_start,
        "tool_versions": {},
        "params": {},
        "key_results": {},
        "error": None,
        "data": {"summary": {}, "files": [], "failed": [], "skipped": [], "warnings": []},
    }
