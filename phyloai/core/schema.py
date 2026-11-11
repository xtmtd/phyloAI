# phyloai/core/schema.py
"""Shared data structures for PhyloAI."""

from __future__ import annotations

import json
import os
import secrets
import stat
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any


COMMON_ALIGNMENT_EXTENSIONS = (
    ".fa",
    ".fas",
    ".fasta",
    ".faa",
    ".fna",
    ".phy",
    ".phylip",
    ".nex",
    ".nxs",
    ".nexus",
)


@dataclass
class MSACollection:
    """A directory of multiple sequence alignment files."""
    directory: Path
    seq_type: str = "AA"          # "AA" or "NT"
    file_extension: str = ".fa"   # ".fa", ".faa", ".fna", ".phy", ".nex"
    count: int = 0                # number of alignment files found

    def __post_init__(self):
        self.directory = Path(self.directory)
        if self.directory.exists():
            if self.file_extension == ".fa":
                self.count = sum(
                    len(list(self.directory.glob(f"*{suffix}")))
                    for suffix in COMMON_ALIGNMENT_EXTENSIONS
                )
            else:
                self.count = len(list(
                    self.directory.glob(f"*{self.file_extension}")
                ))


@dataclass
class TreeSet:
    """A directory of phylogenetic tree files."""
    directory: Path
    format: str = "newick"        # "newick" or "nexus"
    file_extension: str = ".treefile"
    count: int = 0

    def __post_init__(self):
        self.directory = Path(self.directory)
        if self.directory.exists():
            self.count = len(list(
                self.directory.glob(f"*{self.file_extension}")
            ))


@dataclass
class ToolResult:
    """Result of a single external tool invocation."""
    tool: str
    command: str
    returncode: int
    stdout: str
    stderr: str
    wall_time: float             # seconds

    @property
    def success(self) -> bool:
        return self.returncode == 0

    def to_dict(self) -> dict[str, Any]:
        return {
            "tool": self.tool,
            "command": self.command,
            "returncode": self.returncode,
            "stdout": self.stdout,
            "stderr": self.stderr,
            "wall_time": self.wall_time,
            "success": self.success,
        }


@dataclass
class RunRecord:
    """Full record of a PhyloAI analysis run."""
    run_dir: Path
    phyloai_version: str = ""
    steps: list[dict[str, Any]] = field(default_factory=list)

    def __post_init__(self):
        self.run_dir = Path(self.run_dir)
        if not self.phyloai_version:
            from phyloai import __version__
            self.phyloai_version = __version__

    def add_step(self, step_name: str, params: dict, result: ToolResult) -> None:
        self.steps.append({
            "step": step_name,
            "params": params,
            "result": result.to_dict(),
        })

    def to_dict(self) -> dict[str, Any]:
        return {
            "run_dir": str(self.run_dir),
            "phyloai_version": self.phyloai_version,
            "steps": self.steps,
        }


def _create_temp_target(path: Path) -> tuple[Path, int]:
    """Create a unique temporary file in the target directory."""
    for _ in range(100):
        tmp = path.with_name(f".{path.name}.{os.getpid()}.{secrets.token_hex(6)}.tmp")
        try:
            fd = os.open(tmp, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o666)
        except FileExistsError:
            continue
        return tmp, fd
    raise OSError(f"Could not create a temporary file next to {path}")


def save_json_atomic(path: Path, payload: dict[str, Any], *, fsync: bool = False) -> None:
    """Atomically publish JSON, preserving an existing target and its mode.

    A reader never sees a half-written target: the payload is written to a
    temporary file in the same directory and then moved into place with
    ``os.replace``.  On any failure the old target is left untouched.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    mode = stat.S_IMODE(path.stat().st_mode) if path.exists() else None
    tmp, fd = _create_temp_target(path)
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as fh:
            json.dump(payload, fh, indent=2)
            fh.flush()
            if fsync:
                os.fsync(fh.fileno())
        if mode is not None:
            os.chmod(tmp, mode)
        os.replace(tmp, path)
    except BaseException:
        try:
            os.unlink(tmp)
        except OSError:
            pass
        raise


def write_result_json(payload: dict[str, Any], output_dir: Path, *, fsync: bool = False) -> None:
    """Write result.json payload to output directory (shared helper)."""
    output_dir.mkdir(parents=True, exist_ok=True)
    save_json_atomic(output_dir / "result.json", payload, fsync=fsync)
