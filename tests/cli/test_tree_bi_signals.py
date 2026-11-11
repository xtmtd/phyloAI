"""Scheduler-stop behaviour for `tree bi pb` (POSIX signals only).

A scheduler cancels a job with a signal, not with Ctrl+C, and result.json used to
be written only after the shutdown returned. These tests drive the real CLI in a
subprocess because signals cannot be delivered to the in-process runner; the
purely in-process checks live next to the code they cover (`tests/tree/test_bi.py`
for `_install_stop_handlers`, `tests/cli/test_tree_bi.py` for the CLI backup).
"""

from __future__ import annotations

import json
import os
import signal
import subprocess
import sys
import time
from pathlib import Path

import pytest

from tests.helpers import validate_result_json

pytestmark = pytest.mark.skipif(os.name == "nt", reason="POSIX signals are required")

REPO_ROOT = Path(__file__).resolve().parents[2]
CLI_CODE = "from phyloai.cli.main import cli; cli()"

MPIRUN = 'shift 2\nexec "$@"\n'

PB_MPI = r"""for last do :; done
name="$last"
echo $$ > "${name}.pid"
printf 'iter time loglik\n' > "${name}.trace"
i=0
while [ $i -lt 40 ]; do
    i=$((i+1))
    printf '%s 0 -10\n' "$i" >> "${name}.trace"
done
while [ ! -f "${name}.run" ]; do sleep 1; done
printf 'left\n' > "${name}.left"
printf '(a,b);\n' > "${name}.treelist"
printf 'state\n' > "${name}.chain"
"""

# Ignores the soft stop and never leaves its cycle on its own.
PB_MPI_STUCK = PB_MPI.replace(
    'while [ ! -f "${name}.run" ]; do sleep 1; done',
    'while [ ! -f "${name}.run" ]; do sleep 60; done',
)

# Leaves the cycle a few seconds after the soft stop.
PB_MPI_SLOW_STOP = PB_MPI.replace(
    'while [ ! -f "${name}.run" ]; do sleep 1; done',
    'while [ ! -f "${name}.run" ]; do sleep 1; done\nsleep 3',
)

BPCOMP = r"""while [ $# -gt 0 ]; do
    if [ "$1" = "-o" ]; then shift; out=$1; fi
    shift
done
printf 'maxdiff 0.081\nmeandiff 0.006\n' > "${out}.bpdiff"
printf '(a,b);\n' > "${out}.con.tre"
printf 'split\n' > "${out}.bplist"
"""

BPCOMP_SLOW = r"""while [ $# -gt 0 ]; do
    if [ "$1" = "-o" ]; then shift; out=$1; fi
    shift
done
printf 'started\n' > "${out}.started"
sleep 30
"""

TRACECOMP = r"""printf 'name effsize rel_diff\nloglik 312 0.094\n'
"""


def _write_tool(path: Path, body: str) -> None:
    path.write_text("#!/bin/sh\n" + body)
    path.chmod(0o755)


def _make_tools(tmp_path: Path, stuck: bool = False, slow_stop: bool = False,
                slow_bpcomp: bool = False) -> Path:
    tools = tmp_path / "tools"
    tools.mkdir()
    _write_tool(tools / "mpirun", MPIRUN)
    if stuck:
        chain_stub = PB_MPI_STUCK
    elif slow_stop:
        chain_stub = PB_MPI_SLOW_STOP
    else:
        chain_stub = PB_MPI
    _write_tool(tools / "pb_mpi", chain_stub)
    _write_tool(tools / "bpcomp", BPCOMP_SLOW if slow_bpcomp else BPCOMP)
    _write_tool(tools / "tracecomp", TRACECOMP)
    return tools


def _child_env(tools: Path) -> dict[str, str]:
    """Run with the tool directory first on PATH so `mpirun` is the stub, not a real MPI."""
    env = dict(os.environ)
    existing = env.get("PYTHONPATH")
    env["PYTHONPATH"] = str(REPO_ROOT) + (os.pathsep + existing if existing else "")
    env["PATH"] = str(tools) + os.pathsep + env.get("PATH", "")
    return env


def _start_cli(args: list[str], cwd: Path, tools: Path, *, ignore_hup: bool = False) -> subprocess.Popen[str]:
    command = [sys.executable, "-c", CLI_CODE, *args]
    if ignore_hup:
        wrapper = tools / "ignore_hup.sh"
        _write_tool(wrapper, "trap '' HUP\nexec \"$@\"\n")
        command = [str(wrapper), *command]
    return subprocess.Popen(
        command,
        cwd=cwd,
        env=_child_env(tools),
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
    )


def _pb_args(tmp_path: Path, tools: Path, out: Path) -> list[str]:
    matrix = tmp_path / "m.phy"
    matrix.write_text("2 3\na AAA\nb AAA\n")
    return [
        "tree", "bi", "pb",
        "--matrix", str(matrix),
        "--output-dir", str(out),
        "--chains", "2",
        "--threads", "2",
        "--monitor-freq", "100000",
        "--poll-interval", "1",
        "--quiet",
        "--pb-path", str(tools),
    ]


def _wait_for(predicate, timeout: float = 30.0) -> bool:
    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        if predicate():
            return True
        time.sleep(0.1)
    return False


def _wait_started(out: Path, names: tuple[str, str] = ("chain1", "chain2"), timeout: float = 30.0) -> list[int]:
    paths = [out / "chains" / f"{name}.pid" for name in names]
    pids: list[int] | None = None
    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        if all(p.exists() for p in paths):
            try:
                pids = [int(p.read_text().strip()) for p in paths]
            except ValueError:  # the stub is still writing the file
                pids = None
            if pids:
                break
        time.sleep(0.1)
    assert pids, "chains did not start"
    assert _wait_for(lambda: all((out / "chains" / f"{n}.trace").exists() for n in names))
    return pids


def _stop(proc: subprocess.Popen[str], sig: int, timeout: float = 60.0) -> int:
    try:
        proc.send_signal(sig)
        return proc.wait(timeout=timeout)
    finally:
        if proc.poll() is None:
            proc.kill()
            proc.wait()


def _assert_soft_stopped(out: Path, proc: subprocess.Popen[str], rc: int) -> dict:
    assert rc == 0, proc.stdout.read() if proc.stdout else ""
    payload = json.loads((out / "result.json").read_text())
    validate_result_json(payload)
    assert payload["data"]["interrupted"] is True
    assert payload["data"]["shutdown_state"] == "stopped"
    assert "stop_signal" in payload["data"]
    assert not (out / "stop_state.json").exists()
    for name in ("chain1", "chain2"):
        assert (out / "chains" / f"{name}.run").read_text() == "0\n"
    return payload


def test_normal_completion_publishes_completed_result_and_removes_sidecar(tmp_path: Path):
    tools = _make_tools(tmp_path)
    out = tmp_path / "out"
    args = _pb_args(tmp_path, tools, out) + ["--nsamples", "40"]
    proc = _start_cli(args, tmp_path, tools)
    try:
        rc = proc.wait(timeout=60)
    finally:
        if proc.poll() is None:
            proc.kill()
            proc.wait()

    assert rc == 0
    payload = json.loads((out / "result.json").read_text())
    validate_result_json(payload)
    assert payload["status"] == "success"
    assert payload["data"]["interrupted"] is False
    assert payload["data"]["shutdown_state"] == "completed"
    assert not (out / "stop_state.json").exists()


@pytest.mark.parametrize("sig", [signal.SIGTERM, signal.SIGHUP], ids=["sigterm", "sighup"])
def test_stop_signal_stops_the_run_and_writes_result(tmp_path: Path, sig: int):
    """The scheduler's TERM/HUP take the Ctrl+C path and write result.json."""
    tools = _make_tools(tmp_path)
    out = tmp_path / "out"
    proc = _start_cli(_pb_args(tmp_path, tools, out), tmp_path, tools)
    try:
        _wait_started(out)
        rc = _stop(proc, sig)
    finally:
        if proc.poll() is None:
            proc.kill()
            proc.wait()

    _assert_soft_stopped(out, proc, rc)


def test_stop_state_survives_a_hard_kill(tmp_path: Path):
    """The sidecar is written before the wait, so SIGKILL keeps the record."""
    tools = _make_tools(tmp_path, stuck=True)
    out = tmp_path / "out"
    proc = _start_cli(_pb_args(tmp_path, tools, out), tmp_path, tools)
    try:
        _wait_started(out)
        proc.send_signal(signal.SIGTERM)
        assert _wait_for(lambda: (out / "stop_state.json").exists(), timeout=5), "no stop_state.json before the wait"
        assert proc.poll() is None, "the shutdown finished before the kill could land"
        proc.kill()
    finally:
        proc.wait()

    payload = json.loads((out / "stop_state.json").read_text())
    assert payload["state"] == "stop_requested"
    assert payload["stop_signal"] == signal.SIGTERM

    from phyloai.mcp.tools.utils import check_status

    status = check_status(str(out))
    assert status["status"] == "unknown"
    assert status["stop_state"]["state"] == "stop_requested"


def test_inherited_sighup_ignore_is_not_overridden(tmp_path: Path):
    """nohup-style SIG_IGN must survive controller installation."""
    tools = _make_tools(tmp_path)
    out = tmp_path / "out"
    proc = _start_cli(_pb_args(tmp_path, tools, out), tmp_path, tools, ignore_hup=True)
    try:
        _wait_started(out)
        proc.send_signal(signal.SIGHUP)
        time.sleep(1)
        assert proc.poll() is None, "HUP was handled despite inherited SIG_IGN"
        assert not (out / "stop_state.json").exists()
        rc = _stop(proc, signal.SIGTERM)
    finally:
        if proc.poll() is None:
            proc.kill()
            proc.wait()

    _assert_soft_stopped(out, proc, rc)


def test_second_signal_stops_immediately(tmp_path: Path):
    """First signal protects the cycle; the second ends the wait."""
    tools = _make_tools(tmp_path, stuck=True)
    out = tmp_path / "out"
    proc = _start_cli(_pb_args(tmp_path, tools, out), tmp_path, tools)
    try:
        _wait_started(out)
        proc.send_signal(signal.SIGTERM)
        assert _wait_for(lambda: (out / "stop_state.json").exists(), timeout=5)
        assert proc.poll() is None, "the first signal must not end the shutdown"
        proc.send_signal(signal.SIGTERM)
        rc = proc.wait(timeout=40)
    finally:
        if proc.poll() is None:
            proc.kill()
            proc.wait()

    assert rc == 143
    payload = json.loads((out / "result.json").read_text())
    validate_result_json(payload)
    assert payload["status"] == "error"
    assert payload["data"]["shutdown_state"] == "forced"
    assert payload["data"]["force_signal"] == signal.SIGTERM
    assert any("Second stop request" in w for w in payload["data"]["warnings"])


def test_ctrl_c_keeps_waiting_for_the_chain(tmp_path: Path):
    """Ctrl+C must wait for the running cycle instead of being cut short."""
    tools = _make_tools(tmp_path, slow_stop=True)
    out = tmp_path / "out"
    proc = _start_cli(_pb_args(tmp_path, tools, out), tmp_path, tools)
    try:
        _wait_started(out)
        proc.send_signal(signal.SIGINT)
        assert _wait_for(lambda: (out / "stop_state.json").exists(), timeout=10)
        assert proc.poll() is None, "Ctrl+C stopped waiting for the chains"
        left = [out / "chains" / f"{name}.left" for name in ("chain1", "chain2")]
        assert _wait_for(lambda: all(p.exists() for p in left), timeout=30), "chains were cut short"
        rc = proc.wait(timeout=40)
    finally:
        if proc.poll() is None:
            proc.kill()
            proc.wait()

    assert rc == 0
    payload = json.loads((out / "result.json").read_text())
    validate_result_json(payload)
    assert payload["data"]["interrupted"] is True
    assert payload["data"]["shutdown_state"] == "stopped"
    assert payload["data"]["stop_signal"] == signal.SIGINT


def test_second_signal_stops_a_running_convergence_check(tmp_path: Path):
    """The chains leave quickly; the second signal must end the convergence check."""
    tools = _make_tools(tmp_path, slow_bpcomp=True)
    out = tmp_path / "out"
    proc = _start_cli(_pb_args(tmp_path, tools, out), tmp_path, tools)
    try:
        _wait_started(out)
        proc.send_signal(signal.SIGTERM)  # chains leave at once, bpcomp starts
        marker = out / "convergence" / "bpcomp_all.started"
        assert _wait_for(lambda: marker.exists(), timeout=20), "the convergence check did not start"
        assert proc.poll() is None, "the convergence check finished before the second signal"
        proc.send_signal(signal.SIGTERM)
        rc = proc.wait(timeout=40)
    finally:
        if proc.poll() is None:
            proc.kill()
            proc.wait()

    assert rc == 143
    payload = json.loads((out / "result.json").read_text())
    validate_result_json(payload)
    assert payload["status"] == "error"
    assert payload["data"]["shutdown_state"] == "forced"
    assert payload["data"]["force_signal"] == signal.SIGTERM
    assert any("final convergence check skipped" in w for w in payload["data"]["warnings"])
