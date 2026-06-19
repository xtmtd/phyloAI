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
