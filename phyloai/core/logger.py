"""Per-step log file writer for PhyloAI runs."""

from __future__ import annotations
import datetime
from pathlib import Path

from phyloai.core.schema import ToolResult


class StepLogger:
    def __init__(self, run_dir: Path):
        self.log_dir = Path(run_dir) / "logs"
        self.log_dir.mkdir(parents=True, exist_ok=True)

    def write(self, step_name: str, result: ToolResult) -> Path:
        log_path = self.log_dir / f"{step_name}.log"
        timestamp = datetime.datetime.now().isoformat(timespec="seconds")
        entry = (
            f"{'='*60}\n"
            f"timestamp: {timestamp}\n"
            f"tool:      {result.tool}\n"
            f"command:   {result.command}\n"
            f"returncode: {result.returncode}\n"
            f"wall_time: {result.wall_time:.2f}s\n"
            f"--- stdout ---\n{result.stdout}\n"
            f"--- stderr ---\n{result.stderr}\n"
        )
        with open(log_path, "a") as fh:
            fh.write(entry)
        return log_path
