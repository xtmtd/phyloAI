"""Shared IQ-TREE helpers used by tree and posttree modules."""
from __future__ import annotations

import os
import re as _re
import subprocess
from pathlib import Path

from phyloai.core.env import ToolEnv

IQTREE_COMPATIBLE_EXTENSIONS = frozenset({
    ".fa", ".fas", ".fasta", ".faa", ".fna",
    ".phy", ".phylip",
    ".nex", ".nxs", ".nexus",
    ".aln",
})


def _resolve_iqtree_path(iqtree_path: str | None, dry_run: bool) -> str:
    """Resolve IQ-TREE executable: custom path > PATH > bundled.

    Returns an absolute path string.  Raises ValueError for missing
    or non-executable custom paths; raises FileNotFoundError when no
    IQ-TREE can be found on the system.
    """
    if iqtree_path:
        p = Path(iqtree_path).resolve()
        if not p.exists():
            raise ValueError(f"--iqtree-path does not exist: {iqtree_path}")
        if not os.access(str(p), os.X_OK):
            raise ValueError(f"--iqtree-path is not executable: {iqtree_path}")
        return str(p)
    if dry_run:
        return "iqtree3"
    try:
        env = ToolEnv()
        return str(env.require("iqtree3"))
    except FileNotFoundError:
        raise FileNotFoundError(
            "iqtree3 not found. Install from https://github.com/iqtree/iqtree3/releases "
            "or use --iqtree-path."
        )


def _detect_iqtree_version(executable: str) -> dict[str, str]:
    """Detect IQ-TREE version from --version output."""
    exe_name = Path(executable).name
    try:
        proc = subprocess.run(
            [executable, "--version"],
            capture_output=True, text=True, timeout=10,
        )
        combined = proc.stdout + proc.stderr
    except Exception:
        return {exe_name: "unknown"}

    m = _re.search(r"version\s*([\d.]+)", combined, _re.IGNORECASE)
    if m:
        return {exe_name: m.group(1)}
    m = _re.search(r"([\d]+\.[\d]+(?:\.[\d]+)?)", combined)
    if m:
        return {exe_name: m.group(1)}
    return {exe_name: "unknown"}
