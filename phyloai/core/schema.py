# phyloai/core/schema.py
"""Shared data structures for PhyloAI."""

from __future__ import annotations
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
