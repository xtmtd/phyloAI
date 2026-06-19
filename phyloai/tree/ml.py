"""Maximum-likelihood tree inference with FastTree."""

from __future__ import annotations

import shlex
import time
from pathlib import Path
from typing import Any, Callable

from phyloai.core.schema import COMMON_ALIGNMENT_EXTENSIONS

FASTTREE_COMPATIBLE_EXTENSIONS = frozenset({
    ".fa", ".fas", ".fasta", ".faa", ".fna",
    ".phy", ".phylip",
})

FASTTREE_MANAGED_FLAGS = frozenset({
    "-nt", "-gtr", "-lg", "-wag",
    "-cat", "-gamma", "-boot", "-nosupport",
    "-fastest", "-slow",
    "-expert", "-help",
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

    if seq_type == "NT":
        cmd.append("-nt")
        if model == "gtr":
            cmd.append("-gtr")
    else:
        if model == "lg":
            cmd.append("-lg")
        elif model == "wag":
            cmd.append("-wag")

    if mode == "fastest":
        cmd.append("-fastest")
    elif mode == "slow":
        cmd.append("-slow")

    if gamma:
        cmd.append("-gamma")

    cmd.extend(["-cat", str(cat)])

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


import subprocess
import time as _time


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

        if proc.returncode != 0:
            return {
                **result,
                "status": "failed",
                "reason": f"FastTree exited with code {proc.returncode}: {proc.stderr[:200]}",
                "tool_stderr": proc.stderr,
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
                "tool_stderr": proc.stderr,
                "wall_time": wall_time,
                "warnings": warnings,
            }

        return {
            **result,
            "status": "success",
            "wall_time": wall_time,
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
