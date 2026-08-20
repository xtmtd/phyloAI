from __future__ import annotations

import json
import sys
import tempfile
import time
from pathlib import Path

import click

from phyloai.mcp.job import build_cli_argv, launch_cli, read_job_json, write_job_json


def test_write_and_read_job_json_keeps_metadata_outside_output_dir() -> None:
    with tempfile.TemporaryDirectory() as tmp:
        output_dir = Path(tmp) / "output"
        payload = write_job_json(output_dir, pid=12345, command="phyloai pretree align")

        assert not output_dir.exists()
        assert payload["pid"] == 12345
        assert payload["command"] == "phyloai pretree align"
        assert "started_at" in payload
        assert read_job_json(output_dir) == payload


def test_build_cli_argv_uses_click_params() -> None:
    cmd = click.Command(
        "align",
        params=[
            click.Option(["--seq-dir"], required=True),
            click.Option(["--method"], default="linsi"),
            click.Option(["-t", "--threads"], type=int, default=4),
            click.Option(["--tool-args"], default=None),
            click.Option(["--overwrite"], is_flag=True, default=False),
        ],
    )
    argv = build_cli_argv(
        {"command_path": ["pretree", "align"], "click_command": cmd},
        {"seq_dir": "./raw", "method": "auto", "threads": 8, "tool_args": None, "overwrite": True},
    )

    assert argv == ["phyloai", "pretree", "align", "--seq-dir", "./raw", "--method", "auto", "--threads", "8", "--overwrite"]


def test_launch_cli_does_not_report_fast_success_as_failure() -> None:
    cmd = click.Command(
        "success",
        params=[click.Option(["--output-dir"], type=click.Path(path_type=Path), required=True)],
    )
    code = "import sys\n"
    descriptor = {"command_path": ["-c", code], "click_command": cmd, "executable": sys.executable}
    with tempfile.TemporaryDirectory() as tmp:
        output_dir = Path(tmp) / "job"
        result_dir, pid = launch_cli(descriptor, {"output_dir": str(output_dir)}, output_dir)

        assert result_dir == output_dir.resolve()
        assert pid > 0
        assert not (output_dir / "job.json").exists()
        assert read_job_json(result_dir) is None


def test_launch_cli_drains_stderr_for_running_process() -> None:
    cmd = click.Command(
        "stderr",
        params=[click.Option(["--output-dir"], type=click.Path(path_type=Path), required=True)],
    )
    code = (
        "import pathlib, sys, time\n"
        "sys.stderr.write('x' * 131072)\n"
        "sys.stderr.flush()\n"
        "out = pathlib.Path(sys.argv[-1])\n"
        "out.mkdir(parents=True, exist_ok=True)\n"
        "time.sleep(0.5)\n"
    )
    descriptor = {"command_path": ["-c", code], "click_command": cmd, "executable": sys.executable}
    with tempfile.TemporaryDirectory() as tmp:
        output_dir = Path(tmp) / "job"
        result_dir, pid = launch_cli(descriptor, {"output_dir": str(output_dir)}, output_dir)

        deadline = time.time() + 10
        job = None
        while time.time() < deadline:
            job = read_job_json(result_dir)
            if job is not None:
                break
            time.sleep(0.05)

    assert result_dir == output_dir.resolve()
    assert job is not None
    assert job["pid"] == pid


def test_launch_cli_writes_job_json_for_running_process() -> None:
    cmd = click.Command(
        "sleep",
        params=[click.Option(["--output-dir"], type=click.Path(path_type=Path), required=True)],
    )
    code = (
        "import pathlib, sys, time\n"
        "out = pathlib.Path(sys.argv[-1])\n"
        "out.mkdir(parents=True, exist_ok=True)\n"
        "time.sleep(2)\n"
    )
    descriptor = {"command_path": ["-c", code], "click_command": cmd, "executable": sys.executable}
    with tempfile.TemporaryDirectory() as tmp:
        output_dir = Path(tmp) / "job"
        result_dir, pid = launch_cli(descriptor, {"output_dir": str(output_dir)}, output_dir)

        deadline = time.time() + 10
        job = None
        while time.time() < deadline:
            job = read_job_json(result_dir)
            if job is not None:
                break
            time.sleep(0.05)

    assert result_dir == output_dir.resolve()
    assert job is not None
    assert job["pid"] == pid
