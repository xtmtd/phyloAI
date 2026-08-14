"""MCP job lifecycle helpers."""

from __future__ import annotations

import hashlib
import json
import os
import shlex
import subprocess
import threading
import time
from pathlib import Path
from typing import Any

import click


def _job_path(output_dir: Path) -> Path:
    key = hashlib.sha256(str(output_dir.resolve()).encode()).hexdigest()
    return output_dir.parent / ".phyloai-jobs" / f"{key}.json"


def write_job_json(output_dir: Path, pid: int, command: str, *, early_exit_stderr: str = "") -> dict[str, Any]:
    """Write MCP lifecycle metadata outside the command output directory."""
    path = _job_path(output_dir)
    path.parent.mkdir(parents=True, exist_ok=True)
    payload: dict[str, Any] = {
        "pid": pid,
        "started_at": time.strftime("%Y-%m-%dT%H:%M:%S"),
        "command": command,
    }
    if early_exit_stderr:
        payload["early_exit_stderr"] = early_exit_stderr
    with open(path, "w") as fh:
        json.dump(payload, fh, indent=2)
    return payload


def read_job_json(output_dir: Path) -> dict[str, Any] | None:
    """Read ``job.json``; return None when missing or invalid."""
    path = _job_path(output_dir)
    if not path.exists():
        path = output_dir / "job.json"
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
    """Launch a detached CLI command and track only processes still running."""
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
        # Keep consuming stderr so a verbose child cannot block on the pipe.
        threading.Thread(target=proc.stderr.read, daemon=True).start()
        # Wait until CLI passes its conflict check; metadata stays outside output_dir.
        _thread = threading.Thread(target=_write_job_when_ready, args=(output_dir, proc.pid, shlex.join(argv)), daemon=True)
        _thread.start()
        return output_dir, proc.pid

    early = stderr.decode("utf-8", errors="replace")[:1000] if stderr else ""
    if proc.returncode == 0:
        return output_dir, proc.pid
    raise ValueError(f"Process exited immediately with code {proc.returncode}. stderr: {early}")


def _write_job_when_ready(output_dir: Path, pid: int, command: str) -> None:
    """Write lifecycle metadata once *output_dir* exists without polluting it."""
    deadline = time.time() + 30
    while time.time() < deadline:
        if output_dir.is_dir():
            write_job_json(output_dir, pid, command)
            return
        time.sleep(0.05)
