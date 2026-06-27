"""MCP job lifecycle helpers."""

from __future__ import annotations

import json
import os
import shlex
import subprocess
import time
from pathlib import Path
from typing import Any

import click


def write_job_json(output_dir: Path, pid: int, command: str, *, early_exit_stderr: str = "") -> dict[str, Any]:
    """Write `job.json` and return its payload."""
    output_dir.mkdir(parents=True, exist_ok=True)
    payload: dict[str, Any] = {
        "pid": pid,
        "started_at": time.strftime("%Y-%m-%dT%H:%M:%S"),
        "command": command,
    }
    if early_exit_stderr:
        payload["early_exit_stderr"] = early_exit_stderr
    with open(output_dir / "job.json", "w") as fh:
        json.dump(payload, fh, indent=2)
    return payload


def read_job_json(output_dir: Path) -> dict[str, Any] | None:
    """Read `job.json`; return None when missing or invalid."""
    path = output_dir / "job.json"
    if not path.exists():
        return None
    try:
        with open(path) as fh:
            return json.load(fh)
    except (OSError, json.JSONDecodeError):
        return None


def build_cli_argv(descriptor: dict[str, Any], params: dict[str, Any]) -> list[str]:
    """Build a CLI argv list using Click parameter metadata when available."""
    argv = [descriptor.get("executable", "phyloai"), *descriptor["command_path"]]
    command = descriptor.get("click_command")
    if command is None:
        for key, value in params.items():
            _append_param(argv, f"--{key.replace('_', '-')}", value, is_flag=isinstance(value, bool))
        return argv

    for param in command.params:
        if getattr(param, "hidden", False):
            continue
        value = params.get(param.name, param.default if param.default is not ... else None)
        flag = _longest_option(param) or f"--{param.name.replace('_', '-')}"
        _append_param(argv, flag, value, is_flag=isinstance(param, click.Option) and param.is_flag)
    return argv


def _longest_option(param: Any) -> str | None:
    opts = [opt for opt in getattr(param, "opts", []) if opt.startswith("--")]
    if not opts:
        opts = list(getattr(param, "opts", []))
    return max(opts, key=len) if opts else None


def _append_param(argv: list[str], flag: str, value: Any, *, is_flag: bool) -> None:
    if value is None or value is False:
        return
    if is_flag:
        if value is True:
            argv.append(flag)
        return
    argv.extend([flag, str(value)])


def launch_cli(
    descriptor: dict[str, Any],
    params: dict[str, Any],
    output_dir: Path,
    *,
    env: dict[str, str] | None = None,
) -> tuple[Path, int]:
    """Launch a detached CLI command and write `job.json`."""
    output_dir = output_dir.resolve()
    if not output_dir.parent.exists():
        raise ValueError(f"Parent directory does not exist: {output_dir.parent}")
    params = dict(params)
    params.setdefault("output_dir", str(output_dir))
    argv = build_cli_argv(descriptor, params)

    proc_env = dict(os.environ)
    if env:
        proc_env.update(env)
    try:
        proc = subprocess.Popen(
            argv,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.PIPE,
            stdin=subprocess.DEVNULL,
            env=proc_env,
            start_new_session=True,
        )
    except (FileNotFoundError, OSError) as exc:
        raise ValueError(f"Failed to start subprocess: {exc}") from exc

    try:
        _, stderr = proc.communicate(timeout=0.2)
    except subprocess.TimeoutExpired:
        write_job_json(output_dir, proc.pid, shlex.join(argv))
        return output_dir, proc.pid

    early = stderr.decode("utf-8", errors="replace")[:1000] if stderr else ""
    write_job_json(output_dir, proc.pid, shlex.join(argv), early_exit_stderr=early)
    raise ValueError(f"Process exited immediately with code {proc.returncode}. stderr: {early}")
