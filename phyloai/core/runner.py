"""Unified external tool call interface."""

from __future__ import annotations
import subprocess
import time
from pathlib import Path
from typing import Optional

from phyloai.core.schema import ToolResult


class Runner:
    def __init__(self, timeout: int = 86400):
        self.timeout = timeout

    def run(
        self,
        cmd: list[str],
        tool_name: str,
        cwd: Optional[Path] = None,
        env: Optional[dict] = None,
    ) -> ToolResult:
        command_str = " ".join(str(c) for c in cmd)
        start = time.monotonic()
        try:
            proc = subprocess.run(
                [str(c) for c in cmd],
                capture_output=True,
                text=True,
                timeout=self.timeout,
                cwd=str(cwd) if cwd else None,
                env=env,
            )
        except FileNotFoundError:
            raise FileNotFoundError(
                f"Executable not found: '{cmd[0]}'. "
                f"Check 'phyloai doctor' for installation status."
            )
        except subprocess.TimeoutExpired:
            raise TimeoutError(
                f"Tool '{tool_name}' exceeded timeout of {self.timeout}s. "
                f"Command: {command_str}"
            )
        wall_time = time.monotonic() - start
        return ToolResult(
            tool=tool_name,
            command=command_str,
            returncode=proc.returncode,
            stdout=proc.stdout,
            stderr=proc.stderr,
            wall_time=wall_time,
        )
