"""Bayesian phylogenetic inference with PhyloBayes-MPI."""

from __future__ import annotations

import itertools
import json
import math
import os
import re
import shutil
import signal
import subprocess
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from phyloai.core.env import ToolEnv
from phyloai.core.schema import save_json_atomic

MODEL_FLAGS = {
    "gtr": "-gtr",
    "poisson": "-poisson",
    "lg": "-lg",
    "wag": "-wag",
    "jtt": "-jtt",
    "mtrev": "-mtrev",
    "mtzoa": "-mtzoa",
    "mtart": "-mtart",
}

RESUME_ALL = "__ALL__"


# ---------------------------------------------------------------------------
# Pure helpers (Task 2)
# ---------------------------------------------------------------------------


def _resolve_chain_names(chains: int, chain_prefix: str, chain_names: str | None) -> list[str]:
    if chain_names:
        names = [item.strip() for item in chain_names.split(",") if item.strip()]
        if not names:
            raise ValueError("--chain-names must contain at least one non-empty name")
        if len(set(names)) != len(names):
            raise ValueError("--chain-names contains duplicate names")
        return names
    if chains < 1:
        raise ValueError("--chains must be at least 1")
    return [f"{chain_prefix}{i}" for i in range(1, chains + 1)]


def _build_model_flags(model: str, mixture: str, gamma_cats: int, start_tree: Path | None, fix_tree: Path | None) -> list[str]:
    if model not in MODEL_FLAGS:
        raise ValueError(f"Invalid --model: {model}")
    if gamma_cats < 1:
        raise ValueError("--gamma-cats must be at least 1")
    flags: list[str] = []
    if mixture == "auto":
        flags.append("-cat")
    else:
        try:
            ncat = int(mixture)
        except ValueError as exc:
            raise ValueError("--mixture must be 'auto' or a positive integer") from exc
        if ncat < 1:
            raise ValueError("--mixture integer must be at least 1")
        flags.extend(["-ncat", str(ncat)])
    flags.append(MODEL_FLAGS[model])
    flags.extend(["-dgam", str(gamma_cats)])
    if start_tree is not None and fix_tree is not None:
        raise ValueError("--start-tree and --fix-tree are mutually exclusive")
    if start_tree is not None:
        flags.extend(["-t", str(start_tree.resolve())])
    if fix_tree is not None:
        flags.extend(["-T", str(fix_tree.resolve())])
    return flags


def _count_trace_samples(trace_path: Path) -> int:
    if not trace_path.exists():
        return 0
    raw = trace_path.read_bytes()
    if not raw:
        return 0
    lines = raw.splitlines(keepends=True)
    complete = [line.decode(errors="ignore").strip() for line in lines if line.endswith((b"\n", b"\r\n"))]
    rows = [line for line in complete if line]
    if len(rows) <= 1:
        return 0
    return max(0, len(rows) - 1)


# ---------------------------------------------------------------------------
# Tool resolution & command builders (Task 3)
# ---------------------------------------------------------------------------


def _detect_tools(pb_path: Path | None, dry_run: bool) -> dict[str, str]:
    if pb_path is not None:
        tool_paths = {
            "pb_mpi": pb_path / "pb_mpi",
            "bpcomp": pb_path / "bpcomp",
            "tracecomp": pb_path / "tracecomp",
        }
        readpb = pb_path / "readpb_mpi"
        if readpb.exists():
            tool_paths["readpb_mpi"] = readpb
        env = ToolEnv(tool_paths=tool_paths)
    else:
        env = ToolEnv()
    if dry_run:
        return {"pb_mpi": "pb_mpi", "bpcomp": "bpcomp", "tracecomp": "tracecomp", "mpirun": "mpirun"}
    return {name: str(env.require(name)) for name in ["pb_mpi", "bpcomp", "tracecomp", "mpirun"]}


def _build_chain_cmd(mpirun: str, pb_mpi: str, threads: int, matrix: Path, model_flags: list[str], sample_freq: int, nsamples: int, chain_name: str) -> list[str]:
    if threads < 2:
        raise ValueError("--threads must be at least 2")
    if sample_freq < 1:
        raise ValueError("--sample-freq must be at least 1")
    if nsamples != -1 and nsamples < 1:
        raise ValueError("--nsamples must be -1 or a positive integer")
    cmd = [mpirun, "-np", str(threads), pb_mpi, "-d", str(matrix.resolve()), *model_flags, "-x", str(sample_freq), str(nsamples)]
    cmd.append(chain_name)
    return cmd


def _build_resume_cmd(mpirun: str, pb_mpi: str, threads: int, chain_name: str) -> list[str]:
    if threads < 2:
        raise ValueError("--threads must be at least 2")
    return [mpirun, "-np", str(threads), pb_mpi, chain_name]


# ---------------------------------------------------------------------------
# Version detection (Task 3 extension / Task 7)
# ---------------------------------------------------------------------------


def _detect_pb_version(pb_dir: str) -> str | None:
    resolved = Path(pb_dir).resolve() if Path(pb_dir).exists() else Path(pb_dir)
    tool_dir = resolved.parent
    search_dirs = [tool_dir, tool_dir.parent]
    patterns = ["pb_mpi*Manual*.pdf", "pb_mpi*README*", "VERSION", "CHANGELOG", "*.pdf"]
    for sd in list(dict.fromkeys(search_dirs)):
        if not sd.is_dir():
            continue
        for pat in patterns:
            for path in sorted(sd.glob(pat), key=lambda p: str(p).lower()):
                if not path.is_file():
                    continue
                for re_pat in [r"[Vv]ersion\s*(\d+\.\d+(?:\.\d+)?)", r"(\d+\.\d+)"]:
                    m = re.search(re_pat, path.name)
                    if m:
                        return m.group(1)
    return None


def _detect_mpirun_version(mpirun_path: str) -> str | None:
    try:
        proc = subprocess.run([mpirun_path, "--version"], capture_output=True, text=True, timeout=5)
        if proc.returncode == 0:
            m = re.search(r"(\d+\.\d+[\.\d]*)", proc.stdout[:200])
            if m:
                return m.group(1)
    except Exception:
        pass
    return None


def _detect_tool_versions(tools: dict[str, str]) -> dict[str, str | None]:
    pb_ver = _detect_pb_version(tools["pb_mpi"])
    mpirun_ver = _detect_mpirun_version(tools["mpirun"])
    return {
        "pb_mpi": pb_ver,
        "bpcomp": pb_ver,
        "tracecomp": pb_ver,
        "mpirun": mpirun_ver,
    }


def _build_x_flag(burnin: int, sample_freq: int, until: str) -> list[str]:
    parts: list[str] = ["-x", str(burnin)]
    append_sample = sample_freq != 1 or until != "all"
    if append_sample:
        parts.append(str(sample_freq))
    if until != "all":
        parts.append(str(until))
    return parts


# ---------------------------------------------------------------------------
# Run state (Task 4)
# ---------------------------------------------------------------------------


def _state_payload(chain_names: list[str], matrix: Path, model_flags: list[str], sample_freq: int, nsamples: int, threads: int, model: str = "", mixture: str = "", gamma_cats: int = 0, start_tree: str | None = None, fix_tree: str | None = None) -> dict[str, Any]:
    return {
        "chain_names": chain_names,
        "matrix": str(matrix.resolve()),
        "model_flags": model_flags,
        "sample_freq": sample_freq,
        "nsamples": nsamples,
        "threads": threads,
        "model": model,
        "mixture": mixture,
        "gamma_cats": gamma_cats,
        "start_tree": start_tree,
        "fix_tree": fix_tree,
    }


def _write_run_state(output_dir: Path, payload: dict[str, Any]) -> None:
    (output_dir / "run_state.json").write_text(json.dumps(payload, indent=2))


def _read_run_state(output_dir: Path) -> dict[str, Any]:
    state_path = output_dir / "run_state.json"
    if not state_path.exists():
        raise ValueError(f"Missing run_state.json in {output_dir}")
    return json.loads(state_path.read_text())


def _update_run_state_for_new_chains(output_dir: Path, new_names: list[str], current_payload: dict[str, Any]) -> dict[str, Any]:
    existing = _read_run_state(output_dir)
    for key in ["matrix", "model_flags", "sample_freq", "nsamples", "threads", "model", "mixture", "gamma_cats"]:
        if existing.get(key) != current_payload.get(key):
            raise ValueError("Model parameters conflict with existing run_state.json. Use --resume to continue existing chains or choose a different --output-dir.")
    current_names = list(existing.get("chain_names", []))
    overlap = sorted(set(current_names) & set(new_names))
    if overlap:
        raise ValueError(f"Chain name(s) already exist in run_state.json: {', '.join(overlap)}")
    existing["chain_names"] = current_names + new_names
    _write_run_state(output_dir, existing)
    return existing


def _resolve_resume_names(resume: str | None, state: dict[str, Any]) -> list[str]:
    if resume is None:
        return []
    available = list(state.get("chain_names", []))
    if resume == RESUME_ALL:
        return available
    requested = [item.strip() for item in resume.split(",") if item.strip()]
    missing = sorted(set(requested) - set(available))
    if missing:
        raise ValueError(f"Resume chain(s) not found in run_state.json: {', '.join(missing)}")
    return requested


# ---------------------------------------------------------------------------
# FASTA-to-PHYLIP conversion (Task 7)
# ---------------------------------------------------------------------------


def _prepare_matrix(matrix: Path, output_dir: Path, dry_run: bool) -> Path:
    from phyloai.core.formats import FormatConverter, AlignmentFormat

    converter = FormatConverter()
    fmt = converter.detect(matrix)
    if fmt == AlignmentFormat.PHYLIP:
        return matrix
    if fmt == AlignmentFormat.FASTA:
        if dry_run:
            return matrix
        phylip_path = output_dir / "matrix.phy"
        converter.convert(matrix, phylip_path, target=AlignmentFormat.PHYLIP)
        return phylip_path
    raise ValueError(f"Unsupported matrix format: {fmt}. Expected FASTA or PHYLIP.")


# ---------------------------------------------------------------------------
# Convergence parsers & plot (Task 5)
# ---------------------------------------------------------------------------


def _parse_bpcomp_bpdiff(path: Path) -> dict[str, float | None]:
    text = path.read_text() if path.exists() else ""
    maxdiff = None
    meandiff = None
    for line in text.splitlines():
        low = line.lower()
        nums = re.findall(r"[-+]?\d*\.\d+|[-+]?\d+", line)
        if "maxdiff" in low and nums:
            maxdiff = float(nums[-1])
        if "meandiff" in low and nums:
            meandiff = float(nums[-1])
    return {"maxdiff": maxdiff, "meandiff": meandiff}


def _parse_tracecomp_contdiff(path: Path) -> dict[str, float | None]:
    if not path.exists():
        return {"min_effsize": None, "max_rel_diff": None}
    min_effsize = None
    max_rel_diff = None
    for line in path.read_text().splitlines():
        if not line.strip() or line.lower().startswith("name"):
            continue
        parts = line.split()
        if len(parts) < 3:
            continue
        try:
            effsize = float(parts[-2])
            rel_diff = float(parts[-1])
        except ValueError:
            continue
        min_effsize = effsize if min_effsize is None else min(min_effsize, effsize)
        max_rel_diff = rel_diff if max_rel_diff is None else max(max_rel_diff, rel_diff)
    return {"min_effsize": min_effsize, "max_rel_diff": max_rel_diff}


def _status_from_metrics(bp_maxdiff: float | None, min_effsize: float | None, max_rel_diff: float | None) -> str:
    if bp_maxdiff is None or min_effsize is None or max_rel_diff is None:
        return "not converged"
    if bp_maxdiff < 0.1 and min_effsize > 300 and max_rel_diff < 0.1:
        return "good"
    if bp_maxdiff < 0.3 and min_effsize > 50 and max_rel_diff < 0.3:
        return "ok"
    return "not converged"


def _bpcomp_status(maxdiff: float | None) -> str:
    if maxdiff is None:
        return "no"
    if maxdiff < 0.1:
        return "good"
    if maxdiff < 0.3:
        return "ok"
    return "no"


def _tracecomp_status(min_effsize: float | None, max_rel_diff: float | None) -> str:
    if min_effsize is None or max_rel_diff is None:
        return "no"
    if min_effsize > 300 and max_rel_diff < 0.1:
        return "good"
    if min_effsize > 50 and max_rel_diff < 0.3:
        return "ok"
    return "no"


def _generate_trace_plots(trace_paths: list[Path], output_pdf: Path, burnin: int) -> bool:
    try:
        import matplotlib
        matplotlib.use("Agg")
        from matplotlib.backends.backend_pdf import PdfPages
        import matplotlib.pyplot as plt
    except Exception:
        return False
    rows_by_chain: dict[str, tuple[list[str], list[list[str]]]] = {}
    for path in trace_paths:
        if not path.exists():
            continue
        lines = [line.strip().split() for line in path.read_text().splitlines() if line.strip()]
        if len(lines) <= 1:
            continue
        rows_by_chain[path.stem] = (lines[0], lines[1:])
    if not rows_by_chain:
        return False
    first_header = next(iter(rows_by_chain.values()))[0]
    iter_idx = first_header.index("iter") if "iter" in first_header else 0
    columns = [col for col in first_header if col not in {"iter", "time"}]
    # Find iteration value at burnin-th sample from first chain
    burnin_iter: float | None = None
    first_key = next(iter(rows_by_chain.keys()))
    _, first_rows = rows_by_chain[first_key]
    if 0 < burnin <= len(first_rows) and iter_idx < len(first_rows[burnin - 1]):
        try:
            burnin_iter = float(first_rows[burnin - 1][iter_idx])
        except (ValueError, IndexError):
            burnin_iter = None
    output_pdf.parent.mkdir(parents=True, exist_ok=True)
    with PdfPages(output_pdf) as pdf:
        for column in columns:
            fig, ax = plt.subplots(figsize=(8, 4))
            for chain, (header, rows) in rows_by_chain.items():
                if column not in header:
                    continue
                idx = header.index(column)
                y = []
                x = []
                for row in rows:
                    if idx >= len(row) or iter_idx >= len(row):
                        continue
                    try:
                        y.append(float(row[idx]))
                        x.append(float(row[iter_idx]))
                    except ValueError:
                        continue
                ax.plot(x, y, label=chain)
            if burnin_iter is not None:
                ax.axvline(burnin_iter, linestyle="--", color="black", linewidth=0.8)
            ax.set_title(column)
            ax.set_xlabel("iteration")
            ax.set_ylabel(column)
            ax.legend()
            pdf.savefig(fig)
            plt.close(fig)
    return True


# ---------------------------------------------------------------------------
# Convergence runner (Task 6)
# ---------------------------------------------------------------------------


class _ConvergenceInterrupted(Exception):
    """Raised when a stop request interrupts a convergence tool."""

    def __init__(self, *, forced: bool) -> None:
        super().__init__("convergence interrupted")
        self.forced = forced


def _run_tool_interruptible(
    cmd: list[str],
    cwd: Path,
    stop: "StopController | None",
    *,
    abort_if_requested: bool,
) -> subprocess.CompletedProcess[str]:
    if stop is None:
        return subprocess.run(cmd, cwd=cwd, capture_output=True, text=True)
    proc = subprocess.Popen(cmd, cwd=cwd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    try:
        while True:
            try:
                stdout, stderr = proc.communicate(timeout=0.5)
                return subprocess.CompletedProcess(cmd, proc.returncode, stdout, stderr)
            except subprocess.TimeoutExpired:
                if stop.forced or (abort_if_requested and stop.requested):
                    _terminate_end(proc)
                    try:
                        proc.communicate(timeout=2)
                    except (OSError, subprocess.TimeoutExpired):
                        pass
                    raise _ConvergenceInterrupted(forced=stop.forced)
    except BaseException:
        if proc.poll() is None:
            _terminate_end(proc)
        raise


def _run_convergence_check(output_dir: Path, chain_names: list[str], tools: dict[str, str], burnin: int, stop: "StopController | None" = None, *, abort_if_requested: bool = True) -> dict[str, Any]:
    conv_dir = output_dir / "convergence"
    conv_dir.mkdir(parents=True, exist_ok=True)
    comparisons: dict[str, list[str]] = {"all": chain_names}
    for a, b in itertools.combinations(chain_names, 2):
        comparisons[f"{a}_{b}"] = [a, b]
    warnings: list[str] = []
    result: dict[str, Any] = {"all_chains": {}, "pairwise": {}}
    for label, names in comparisons.items():
        base = f"bpcomp_{'all' if label == 'all' else label}"
        bp_base = conv_dir / base
        bp_cmd = [tools["bpcomp"], "-x", str(burnin), "-o", base, *[str(Path("../chains") / name) for name in names]]
        bp_proc = _run_tool_interruptible(bp_cmd, conv_dir, stop, abort_if_requested=abort_if_requested)
        if bp_proc.returncode != 0:
            warnings.append(f"bpcomp {label} exited with code {bp_proc.returncode}: {bp_proc.stderr[:200]}")
        trace_base = f"tracecomp_{'all' if label == 'all' else label}"
        trace_out = conv_dir / f"{trace_base}.contdiff"
        trace_cmd = [tools["tracecomp"], "-x", str(burnin), *[str(Path("../chains") / f"{name}.trace") for name in names]]
        trace_proc = _run_tool_interruptible(trace_cmd, conv_dir, stop, abort_if_requested=abort_if_requested)
        if trace_proc.returncode != 0:
            warnings.append(f"tracecomp {label} exited with code {trace_proc.returncode}: {trace_proc.stderr[:200]}")
        trace_out.write_text((trace_proc.stdout or "") + (("\n" + trace_proc.stderr) if trace_proc.stderr else ""))
        bp = _parse_bpcomp_bpdiff(bp_base.with_suffix(".bpdiff"))
        tr = _parse_tracecomp_contdiff(trace_out)
        metrics = {
            "bpcomp_maxdiff": bp["maxdiff"],
            "bpcomp_meandiff": bp["meandiff"],
            "tracecomp_min_effsize": tr["min_effsize"],
            "tracecomp_max_reldiff": tr["max_rel_diff"],
            "status": _status_from_metrics(bp["maxdiff"], tr["min_effsize"], tr["max_rel_diff"]),
        }
        if label == "all":
            result["all_chains"] = metrics
        else:
            result["pairwise"][label] = metrics
    try:
        _generate_trace_plots([output_dir / "chains" / f"{name}.trace" for name in chain_names], conv_dir / "trace_plots.pdf", burnin)
    except Exception:
        pass
    result["warnings"] = warnings
    return result


# ---------------------------------------------------------------------------
# Result assembly (Task 7)
# ---------------------------------------------------------------------------


def _assemble_result(params: dict[str, Any], command: str, wall_time: float, tool_versions: dict[str, str | None], chain_cmds: dict[str, list[str]], chain_lengths: dict[str, int], final_convergence: dict[str, Any] | None, tool_outputs: dict[str, str], interrupted: bool, status: str = "success", error: str | None = None, warnings: list[str] | None = None, output_dir: Path | None = None, shutdown_state: str | None = "completed", stop_signal: int | None = None, force_signal: int | None = None) -> dict[str, Any]:
    all_warnings = list(warnings or [])
    if final_convergence:
        all_warnings.extend(final_convergence.get("warnings", []))
    consensus_tree: str | None = None
    output_files: dict[str, dict[str, str]] = {}
    if final_convergence is not None and output_dir is not None:
        candidate = output_dir / "convergence" / "bpcomp_all.con.tre"
        if candidate.exists():
            consensus_tree = "convergence/bpcomp_all.con.tre"
        trace_pdf = output_dir / "convergence" / "trace_plots.pdf"
        if trace_pdf.exists():
            output_files["trace_plots"] = {"path": str(trace_pdf), "description": "MCMC trace plots showing parameter sampling over iterations for all chains"}
        bpcomp_patterns = sorted((output_dir / "convergence").glob("bpcomp_*"))
        for p in bpcomp_patterns:
            key = p.name.replace(".", "_")
            output_files[key] = {"path": str(p), "description": f"PhyloBayes bpcomp output: {p.name}"}
        tracecomp_patterns = sorted((output_dir / "convergence").glob("tracecomp_*"))
        for p in tracecomp_patterns:
            key = p.name.replace(".", "_")
            output_files[key] = {"path": str(p), "description": f"PhyloBayes tracecomp output: {p.name}"}
    data: dict[str, Any] = {
        "chain_cmds": chain_cmds,
        "tool_stderr": tool_outputs,
        "tool_logs": {name: f"chains/{name}.log" for name in chain_cmds},
        "output_files": output_files,
        "interrupted": interrupted,
        "warnings": all_warnings,
    }
    if shutdown_state is not None:
        data["shutdown_state"] = shutdown_state
    if stop_signal is not None:
        data["stop_signal"] = stop_signal
    if force_signal is not None:
        data["force_signal"] = force_signal
    return {
        "status": status,
        "command": command,
        "wall_time": wall_time,
        "tool_versions": tool_versions,
        "params": params,
        "key_results": {
            "chain_names": list(chain_cmds.keys()),
            "chain_lengths": chain_lengths,
            "final_convergence": {k: v for k, v in (final_convergence or {}).items() if k != "warnings"},
            "consensus_tree": consensus_tree,
        },
        "error": error,
        "data": data,
    }


# ---------------------------------------------------------------------------
# Process launch & monitoring (Task 8)
# ---------------------------------------------------------------------------


def _write_stop_state(
    output_dir: Path,
    *,
    state: str,
    pid: int,
    stop_signal: int,
    force_signal: int | None = None,
) -> None:
    """Write/update the minimal stop sidecar atomically."""
    if state not in ("stop_requested", "forced"):
        raise ValueError(f"Unknown stop state: {state}")
    path = output_dir / "stop_state.json"
    now = datetime.now(timezone.utc).isoformat()
    requested_at = now
    if path.exists():
        try:
            existing = json.loads(path.read_text())
            requested_at = existing.get("requested_at", now)
        except (OSError, json.JSONDecodeError):
            pass
    payload: dict[str, Any] = {
        "state": state,
        "pid": pid,
        "requested_at": requested_at,
        "stop_signal": stop_signal,
    }
    if force_signal is not None:
        payload["forced_at"] = now
        payload["force_signal"] = force_signal
    save_json_atomic(path, payload, fsync=True)


def _clear_stop_state(output_dir: Path) -> None:
    (output_dir / "stop_state.json").unlink(missing_ok=True)


def finalize_stop_payload(
    stop: "StopController",
    payload: dict[str, Any],
    *,
    output_dir: Path,
) -> dict[str, Any]:
    """Process a late stop request after ``run_bi_pb()`` returned.

    This is the library-side handoff used by the CLI final publisher; it never
    restarts analysis or chains. A late first request after a completed result
    is treated as already handled (no chains remain); a late second request
    converts the payload to the forced-stop error form.
    """
    data = payload.setdefault("data", {})
    warnings = data.setdefault("warnings", [])
    if stop.forced:
        stop_signal = stop.first_signal if stop.first_signal is not None else getattr(signal, "SIGINT", 2)
        force_signal = stop.force_signal if stop.force_signal is not None else stop_signal
        payload["status"] = "error"
        if not payload.get("error"):
            payload["error"] = "Run forcibly stopped by a second stop signal"
        data["interrupted"] = True
        data["shutdown_state"] = "forced"
        data["stop_signal"] = stop_signal
        data["force_signal"] = force_signal
        if "Second stop request" not in " ".join(warnings):
            warnings.append("Second stop request: terminating the chains now")
        try:
            _write_stop_state(
                output_dir,
                state="forced",
                pid=os.getpid(),
                stop_signal=stop_signal,
                force_signal=force_signal,
            )
        except OSError as exc:
            warnings.append(f"Could not update stop_state.json: {exc}")
        stop.mark_handled(2)
        return payload

    if stop.requested:
        # No chains remain in this handoff; a first request cannot do more
        # than be recorded as handled before commit.
        stop.mark_handled(1)
    return payload


class StopController:
    """Command-local stop state and signal lifecycle for ``tree bi pb``.

    The handler only records signals; it does no I/O and never raises. The
    library drives soft/forced shutdown by polling the state below.
    """

    def __init__(self, *, install_handlers: bool = True) -> None:
        self._install_handlers = install_handlers
        self.first_signal: int | None = None
        self.force_signal: int | None = None
        self.stop_level = 0
        self.handled_level = 0
        self.phase = "running"
        self._latched: list[int] = []
        self._previous: dict[int, Any] = {}
        self._installed = False

    @property
    def requested(self) -> bool:
        return self.stop_level >= 1

    @property
    def forced(self) -> bool:
        return self.stop_level >= 2

    @property
    def pending_request(self) -> bool:
        return self.stop_level > self.handled_level

    @property
    def latched_signals(self) -> tuple[int, ...]:
        return tuple(self._latched)

    def _signals(self) -> list[int]:
        signals: list[int] = []
        for name in ("SIGINT", "SIGTERM", "SIGHUP"):
            signum = getattr(signal, name, None)
            if signum is not None:
                signals.append(signum)
        return signals

    def install(self) -> None:
        if not self._install_handlers or self._installed:
            return
        for signum in self._signals():
            try:
                current = signal.getsignal(signum)
                if current == signal.SIG_IGN:
                    continue
                self._previous[signum] = current
                signal.signal(signum, self._handle)
            except (ValueError, OSError):
                self._previous.pop(signum, None)
        self._installed = True

    def restore(self) -> None:
        for signum, handler in list(self._previous.items()):
            try:
                signal.signal(signum, handler)
            except (ValueError, OSError):
                pass
        self._previous.clear()
        self._installed = False

    def _handle(self, signum: int, _frame: Any) -> None:
        self.note_signal(signum)

    def note_signal(self, signum: int) -> None:
        """Record a stop request without doing any I/O."""
        if self.phase == "committing":
            self._latched.append(signum)
            return
        if self.stop_level == 0:
            self.first_signal = signum
            self.stop_level = 1
            self.phase = "stop_requested"
        elif self.stop_level == 1:
            self.force_signal = signum
            self.stop_level = 2
            self.phase = "forced"
        else:
            self._latched.append(signum)

    def note_keyboard_interrupt(self) -> None:
        """Direct-library-call fallback: treat KeyboardInterrupt as SIGINT."""
        self.note_signal(getattr(signal, "SIGINT", 2))

    def mark_handled(self, level: int) -> None:
        """Mark the stop action as submission-ready up to *level*."""
        self.handled_level = max(self.handled_level, level)

    def begin_commit(self) -> bool:
        """Enter the commit attempt and report whether a request is pending."""
        self.phase = "committing"
        if self.pending_request:
            if self.stop_level >= 2:
                self.phase = "forced"
            elif self.stop_level == 1:
                self.phase = "stop_requested"
            else:
                self.phase = "running"
            return False
        return True

    def __enter__(self) -> "StopController":
        self.install()
        return self

    def __exit__(self, _exc_type: Any, _exc: Any, _tb: Any) -> None:
        self.restore()


def _soft_stop_chains(output_dir: Path, chain_names: list[str]) -> None:
    for name in chain_names:
        (output_dir / "chains" / f"{name}.run").write_text("0\n")


def _terminate_end(proc: subprocess.Popen[str]) -> None:
    """Last resort for a chain that ignored the soft stop: TERM, then KILL."""
    try:
        proc.terminate()
        proc.wait(timeout=2)
    except (OSError, subprocess.TimeoutExpired):
        try:
            proc.kill()
            try:
                proc.wait(timeout=2)
            except subprocess.TimeoutExpired:
                pass
        except OSError:
            pass


def _format_convergence_text(conv_result: dict[str, Any]) -> str:
    """Render convergence stats as ASCII-only text matching the design spec layout.
    Handles None metrics (from failed convergence tools) gracefully."""
    lines: list[str] = []
    all_c = conv_result.get("all_chains", {})
    if all_c:
        bpmax = all_c.get("bpcomp_maxdiff")
        bpmean = all_c.get("bpcomp_meandiff")
        eff = all_c.get("tracecomp_min_effsize")
        rel = all_c.get("tracecomp_max_reldiff")
        bp_st = _bpcomp_status(bpmax)
        tc_st = _tracecomp_status(eff, rel)
        bp = f"  bpcomp    maxdiff  {bpmax:.3f}" if bpmax is not None else "  bpcomp    maxdiff  --"
        if bpmean is not None:
            bp += f"   meandiff  {bpmean:.3f}"
        tr = f"  tracecomp  min effsize  {eff:.0f}" if eff is not None else "  tracecomp  min effsize  --"
        if rel is not None:
            tr += f"   max rel_diff  {rel:.3f}"
        lines.append(f"\n  All chains")
        lines.append(bp + f"   [{bp_st}]")
        lines.append(tr + f"   [{tc_st}]")
    pw = conv_result.get("pairwise", {})
    if pw:
        lines.append("\n  Pairwise")
        lines.append("    pair              maxdiff  min effsize  max rel_diff  bpcomp  tracecomp")
        for label, m in pw.items():
            md = m.get("bpcomp_maxdiff")
            es = m.get("tracecomp_min_effsize")
            rd = m.get("tracecomp_max_reldiff")
            md_s = f"{md:.3f}" if md is not None else "--"
            es_s = f"{es:.0f}" if es is not None else "--"
            rd_s = f"{rd:.3f}" if rd is not None else "--"
            bp_st2 = _bpcomp_status(md)
            tc_st2 = _tracecomp_status(es, rd)
            pair_name = label.replace("_", " x ")
            lines.append(f"    {pair_name:<18} {md_s:^7}  {es_s:^11}  {rd_s:^12}  {bp_st2:^6}  {tc_st2:^9}")
    return "\n".join(lines)


def _run_bi_processes(output_dir: Path, chain_names: list[str], chain_cmds: dict[str, list[str]], tools: dict[str, str], tool_versions: dict[str, str | None], params: dict[str, Any], command: str, run_start: float, nsamples: int, monitor_freq: int, burnin_frac: float, poll_interval: int, quiet: bool, stop: StopController | None = None) -> dict[str, Any]:
    stop = stop or StopController(install_handlers=False)
    stop.install()

    procs: dict[str, subprocess.Popen[str]] = {}
    outputs: dict[str, str] = {name: "" for name in chain_names}
    conv_warnings: list[str] = []
    live_display: Any = None
    sidecar_written = False
    soft_stop_done = False
    forced_cleanup_done = False

    def _close_log(proc: subprocess.Popen[str]) -> None:
        log = getattr(proc, "_phyloai_log", None)
        if log is not None and not log.closed:
            log.close()

    def _stop_display() -> None:
        if live_display is not None:
            try:
                live_display.stop()
            except Exception:
                pass

    def _terminate_all() -> None:
        for proc in list(procs.values()):
            if proc.poll() is None:
                _terminate_end(proc)
            _close_log(proc)

    def _write_requested_sidecar() -> None:
        nonlocal sidecar_written
        if sidecar_written:
            return
        try:
            _write_stop_state(
                output_dir,
                state="stop_requested",
                pid=os.getpid(),
                stop_signal=stop.first_signal if stop.first_signal is not None else getattr(signal, "SIGINT", 2),
            )
        except OSError as exc:
            conv_warnings.append(f"Could not write stop_state.json: {exc}")
        sidecar_written = True

    def _warn_once(message: str) -> None:
        if message not in conv_warnings:
            conv_warnings.append(message)

    def _write_forced_sidecar() -> None:
        try:
            _write_stop_state(
                output_dir,
                state="forced",
                pid=os.getpid(),
                stop_signal=stop.first_signal if stop.first_signal is not None else getattr(signal, "SIGINT", 2),
                force_signal=stop.force_signal if stop.force_signal is not None else getattr(signal, "SIGINT", 2),
            )
        except OSError as exc:
            conv_warnings.append(f"Could not update stop_state.json: {exc}")

    def _force_cleanup() -> None:
        nonlocal forced_cleanup_done
        if forced_cleanup_done:
            return
        _warn_once("Second stop request: terminating the chains now")
        _write_forced_sidecar()
        _terminate_all()
        forced_cleanup_done = True
        stop.mark_handled(2)

    def execute() -> dict[str, Any] | None:
        nonlocal live_display
        if stop.requested:
            return None
        from rich.console import Console
        from rich.progress import BarColumn, Progress, TaskProgressColumn, TextColumn, TimeElapsedColumn
        from rich.live import Live

        console = Console()
        progress = Progress(
            TextColumn(" {task.description}"),
            BarColumn(),
            TaskProgressColumn(),
            TextColumn(" {task.fields[samples]} samples"),
            TimeElapsedColumn(),
            console=console,
        )
        for name in chain_names:
            init = 0
            trace_path = output_dir / "chains" / f"{name}.trace"
            if trace_path.exists():
                init = _count_trace_samples(trace_path)
            progress.add_task(name, samples=init,
                              total=None if nsamples == -1 else nsamples,
                              completed=init)
        live_display = Live(progress, console=console, refresh_per_second=1 / 60)
        if not quiet:
            live_display.start()

        for name, cmd in chain_cmds.items():
            if stop.requested:
                return None
            log = open(output_dir / "chains" / f"{name}.log", "w")
            proc = subprocess.Popen(cmd, cwd=output_dir / "chains", stdout=log, stderr=subprocess.STDOUT, text=True)
            proc._phyloai_log = log  # type: ignore[attr-defined]
            procs[name] = proc
            if not quiet:
                print(f"  Started {name} (pid {proc.pid}): {' '.join(cmd)}", flush=True)

        last_check = 0
        trace_lengths: dict[str, int] = {}
        for name in chain_names:
            trace_path = output_dir / "chains" / f"{name}.trace"
            trace_lengths[name] = _count_trace_samples(trace_path) if trace_path.exists() else 0
        last_trace_read = time.monotonic()
        final_convergence: dict[str, Any] | None = None
        while procs:
            if stop.requested:
                return None
            now = time.monotonic()
            if now - last_trace_read >= poll_interval:
                trace_lengths = {name: _count_trace_samples(output_dir / "chains" / f"{name}.trace") for name in chain_names}
                last_trace_read = now
                for tid, task in enumerate(progress.tasks):
                    name = chain_names[tid]
                    progress.update(task.id, samples=trace_lengths.get(name, 0),
                                    completed=trace_lengths.get(name, 0),
                                    total=nsamples if nsamples != -1 else None)
                if nsamples != -1:
                    for name in list(procs):
                        if trace_lengths.get(name, 0) >= nsamples:
                            (output_dir / "chains" / f"{name}.run").write_text("0\n")
                min_len = min(trace_lengths.values()) if trace_lengths else 0
                if min_len - last_check >= monitor_freq:
                    burnin = math.floor(min_len * burnin_frac)
                    if burnin >= 10:
                        try:
                            conv = _run_convergence_check(output_dir, chain_names, tools, burnin, stop=stop, abort_if_requested=True)
                        except _ConvergenceInterrupted:
                            return None
                        final_convergence = conv
                        if not quiet:
                            conv_text = _format_convergence_text(conv)
                            (output_dir / "convergence" / "convergence_render.txt").write_text(conv_text)
                            _stop_display()
                            print(f"\n--- Convergence Check @ {min_len} samples (burnin {int(burnin_frac*100)}% = {burnin}) ---{conv_text}")
                            print("-" * 60)
                            pairwise_statuses = [m.get("status") for m in conv.get("pairwise", {}).values()]
                            all_status = conv.get("all_chains", {}).get("status")
                            n_good = pairwise_statuses.count("good")
                            n_ok = pairwise_statuses.count("ok")
                            n_not = pairwise_statuses.count("not converged")
                            if all_status == "good" and n_not == 0 and n_ok == 0:
                                print("  *** All convergence criteria met (all pairs good). You may stop chains with Ctrl+C when ready. ***")
                            elif n_not == 0:
                                print(f"  Convergence acceptable across all chain pairs ({n_good} good, {n_ok} ok). Consider stopping when ready.")
                            elif n_good + n_ok >= 1:
                                print(f"  Some chain pairs agree ({n_good} good, {n_ok} ok, {n_not} not converged).")
                            print("\n\n\n\n")
                            live_display.start()
                    else:
                        msg = f"Skipping convergence check: chains too short (burnin {burnin} < 10)"
                        conv_warnings.append(msg)
                        if not quiet:
                            _stop_display()
                            print(f"  Warning: {msg}")
                            print("\n\n\n\n")
                            live_display.start()
                    last_check = min_len
            for name, proc in list(procs.items()):
                if proc.poll() is not None:
                    _close_log(proc)
                    if proc.returncode != 0:
                        remaining = [n for n in procs if n != name]
                        if remaining:
                            if not quiet:
                                _stop_display()
                                print(f"  Chain {name} failed -- stopping remaining chains...")
                            _soft_stop_chains(output_dir, remaining)
                            if stop.requested:
                                _write_requested_sidecar()
                            for rem_name in remaining:
                                rem_proc = procs[rem_name]
                                while rem_proc.poll() is None:
                                    if stop.forced:
                                        break
                                    if stop.requested:
                                        _write_requested_sidecar()
                                    try:
                                        rem_proc.wait(timeout=1)
                                    except subprocess.TimeoutExpired:
                                        continue
                                if stop.forced:
                                    break
                                _close_log(rem_proc)
                            if stop.forced:
                                _force_cleanup()
                            elif stop.requested:
                                _write_requested_sidecar()
                                _warn_once("Stop requested while handling a failed chain")
                                stop.mark_handled(1)
                        if not quiet:
                            _stop_display()
                        trace_lengths = {n: _count_trace_samples(output_dir / "chains" / f"{n}.trace") for n in chain_names}
                        for n in chain_names:
                            log_path = output_dir / "chains" / f"{n}.log"
                            outputs[n] = log_path.read_text() if log_path.exists() else ""
                        shutdown_state = "forced" if stop.forced else ("stopped" if stop.requested else None)
                        return _assemble_result(
                            params, command, time.monotonic() - run_start, tool_versions,
                            chain_cmds, trace_lengths, final_convergence, outputs,
                            False, "error", f"pb_mpi chain {name} exited with code {proc.returncode}",
                            warnings=conv_warnings, output_dir=output_dir,
                            shutdown_state=shutdown_state,
                            stop_signal=stop.first_signal if stop.requested else None,
                            force_signal=stop.force_signal if stop.forced else None,
                        )
                    del procs[name]
            if stop.requested:
                return None
            time.sleep(1)
        _stop_display()
        lengths = {name: _count_trace_samples(output_dir / "chains" / f"{name}.trace") for name in chain_names}
        if lengths:
            burnin = math.floor(min(lengths.values()) * burnin_frac)
            if burnin >= 10:
                try:
                    final_convergence = _run_convergence_check(output_dir, chain_names, tools, burnin, stop=stop, abort_if_requested=True)
                except _ConvergenceInterrupted:
                    return None
            else:
                conv_warnings.append(f"Skipping convergence check: chains too short (burnin {burnin} < 10)")
        for name in chain_names:
            log_path = output_dir / "chains" / f"{name}.log"
            outputs[name] = log_path.read_text() if log_path.exists() else ""
        if stop.requested:
            return None
        return _assemble_result(params, command, time.monotonic() - run_start, tool_versions, chain_cmds, lengths, final_convergence, outputs, False, warnings=conv_warnings, output_dir=output_dir, shutdown_state="completed")

    def stop_flow() -> dict[str, Any]:
        nonlocal sidecar_written, soft_stop_done, forced_cleanup_done
        if not sidecar_written:
            _stop_display()
            if not quiet:
                print("\n  Caught interrupt -- sending soft-stop to all chains...")
            _write_requested_sidecar()
        if not soft_stop_done:
            _soft_stop_chains(output_dir, chain_names)
            soft_stop_done = True
            if not quiet:
                for name in chain_names:
                    print(f"    Wrote 0 -> chains/{name}.run")
                print("    Waiting for chains to finish current cycle (signal again to stop now)...")

        if not stop.forced:
            for name, proc in list(procs.items()):
                while proc.poll() is None:
                    if stop.forced:
                        break
                    try:
                        proc.wait(timeout=1)
                    except subprocess.TimeoutExpired:
                        continue
                if stop.forced:
                    break
                _close_log(proc)
                if not quiet:
                    samples = _count_trace_samples(output_dir / "chains" / f"{name}.trace")
                    print(f"    {name} stopped at {samples} samples.")
            if not stop.forced:
                stop.mark_handled(1)

        if stop.forced and not forced_cleanup_done:
            _force_cleanup()

        lengths = {name: _count_trace_samples(output_dir / "chains" / f"{name}.trace") for name in chain_names}
        final_convergence: dict[str, Any] | None = None
        if not stop.forced and lengths:
            burnin = math.floor(min(lengths.values()) * burnin_frac)
            if burnin >= 10:
                if not quiet:
                    print("    Running final convergence check...")
                try:
                    final_convergence = _run_convergence_check(output_dir, chain_names, tools, burnin, stop=stop, abort_if_requested=False)
                except _ConvergenceInterrupted:
                    conv_warnings.append("Second stop request: final convergence check skipped")
                    stop.mark_handled(2)
                except KeyboardInterrupt:
                    stop.note_keyboard_interrupt()
                    raise
            else:
                conv_warnings.append(f"Skipping convergence check: chains too short (burnin {burnin} < 10)")

        if stop.forced and not forced_cleanup_done:
            _force_cleanup()

        if not stop.forced:
            stop.mark_handled(1)
            for name in chain_names:
                log_path = output_dir / "chains" / f"{name}.log"
                outputs[name] = log_path.read_text() if log_path.exists() else ""
        if stop.forced:
            final_convergence = None
            if not forced_cleanup_done:
                _force_cleanup()

        shutdown_state = "forced" if stop.forced else "stopped"
        status = "error" if stop.forced else "success"
        error = "Run forcibly stopped by a second stop signal" if stop.forced else None
        if not quiet:
            print(f"    Writing result.json  (status: {status})")
        return _assemble_result(
            params, command, time.monotonic() - run_start, tool_versions,
            chain_cmds, lengths, final_convergence, outputs,
            True, status, error, warnings=conv_warnings, output_dir=output_dir,
            shutdown_state=shutdown_state,
            stop_signal=stop.first_signal,
            force_signal=stop.force_signal,
        )

    while True:
        try:
            if stop.requested:
                return stop_flow()
            result = execute()
            if result is not None:
                return result
        except KeyboardInterrupt:
            stop.note_keyboard_interrupt()
        except _ConvergenceInterrupted:
            return stop_flow()


# ---------------------------------------------------------------------------
# Main entry point (Task 7)
# ---------------------------------------------------------------------------


def run_bi_pb(
    matrix: Path | None,
    output_dir: Path = Path("runs/tree/bi"),
    overwrite: bool = False,
    model: str = "gtr",
    mixture: str = "auto",
    gamma_cats: int = 4,
    start_tree: Path | None = None,
    fix_tree: Path | None = None,
    chains: int = 3,
    chain_prefix: str = "chain",
    chain_names: str | None = None,
    threads: int = 4,
    sample_freq: int = 1,
    nsamples: int | None = None,
    resume: str | None = None,
    monitor_freq: int = 100,
    burnin_frac: float = 0.5,
    poll_interval: int = 60,
    pb_path: Path | None = None,
    dry_run: bool = False,
    quiet: bool = False,
    *,
    stop: StopController | None = None,
) -> dict[str, Any]:
    run_start = time.monotonic()
    output_dir = output_dir.resolve()
    if overwrite and resume is not None:
        raise ValueError("--overwrite and --resume are mutually exclusive")
    if resume is None and matrix is None:
        raise ValueError("--matrix is required unless --resume is used")
    if matrix is not None and resume is None and not matrix.exists():
        raise ValueError(f"--matrix does not exist: {matrix}")
    if burnin_frac < 0.0 or burnin_frac >= 1.0:
        raise ValueError("--burnin-frac must be 0.0 <= x < 1.0")
    model_flags = _build_model_flags(model, mixture, gamma_cats, start_tree, fix_tree)
    tools = _detect_tools(pb_path, dry_run)
    tool_versions = _detect_tool_versions(tools) if not dry_run else {"pb_mpi": None, "bpcomp": None, "tracecomp": None, "mpirun": None}
    nsamples_updated = False
    if resume is None:
        if nsamples is None:
            nsamples = -1
        names = _resolve_chain_names(chains, chain_prefix, chain_names)
        assert matrix is not None
        state = _state_payload(names, matrix, model_flags, sample_freq, nsamples, threads,
                               model=model, mixture=mixture, gamma_cats=gamma_cats,
                               start_tree=str(start_tree) if start_tree else None,
                               fix_tree=str(fix_tree) if fix_tree else None)
    else:
        state = _read_run_state(output_dir)
        names = _resolve_resume_names(resume, state)
        matrix = Path(state["matrix"])
        model_flags = list(state["model_flags"])
        sample_freq = int(state["sample_freq"])
        stored_nsamples = int(state["nsamples"])
        threads = int(state["threads"])
        model = state.get("model", model)
        mixture = state.get("mixture", mixture)
        gamma_cats = state.get("gamma_cats", gamma_cats)
        start_tree = Path(state["start_tree"]) if state.get("start_tree") else None
        fix_tree = Path(state["fix_tree"]) if state.get("fix_tree") else None
        if nsamples is not None:
            if nsamples != -1 and nsamples < 1:
                raise ValueError("--nsamples must be -1 or a positive integer")
            state["nsamples"] = nsamples
            nsamples_updated = True
        else:
            nsamples = stored_nsamples
            nsamples_updated = False
    if not dry_run:
        # --- conflict check BEFORE any filesystem writes ---
        if not overwrite and output_dir.exists():
            if resume is not None:
                pass
            elif chain_names and (output_dir / "run_state.json").exists():
                pass
            elif any(output_dir.iterdir()):
                raise ValueError(f"Output directory {output_dir} already exists and is non-empty. Use --overwrite to replace.")
        if overwrite and output_dir.exists():
            shutil.rmtree(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        (output_dir / "chains").mkdir(exist_ok=True)
        (output_dir / "convergence").mkdir(exist_ok=True)
        if resume is None:
            if (output_dir / "run_state.json").exists() and chain_names:
                _update_run_state_for_new_chains(output_dir, names, state)
            else:
                _write_run_state(output_dir, state)
        elif nsamples_updated:
            _write_run_state(output_dir, state)
    # --- auto-convert FASTA to PHYLIP for pb_mpi (skip on resume: matrix unused) ---
    if matrix is not None and resume is None:
        matrix_for_chains = _prepare_matrix(matrix, output_dir, dry_run)
    else:
        matrix_for_chains = None

    params = {
        "matrix": str(matrix) if matrix else None,
        "output_dir": str(output_dir),
        "overwrite": overwrite,
        "model": model,
        "mixture": mixture,
        "gamma_cats": gamma_cats,
        "start_tree": str(start_tree) if start_tree else None,
        "fix_tree": str(fix_tree) if fix_tree else None,
        "chains": chains,
        "chain_prefix": chain_prefix,
        "chain_names": chain_names,
        "threads": threads,
        "sample_freq": sample_freq,
        "nsamples": nsamples,
        "resume": resume,
        "monitor_freq": monitor_freq,
        "burnin_frac": burnin_frac,
        "poll_interval": poll_interval,
        "pb_path": str(pb_path) if pb_path else None,
        "dry_run": dry_run,
        "quiet": quiet,
    }
    command_parts = ["phyloai", "tree", "bi", "pb"]
    if resume is None and matrix is not None:
        command_parts.extend(["--matrix", str(matrix)])
    command_parts.extend(["--output-dir", str(output_dir)])
    if overwrite:
        command_parts.append("--overwrite")
    command_parts.extend(["--model", model, "--mixture", mixture, "--gamma-cats", str(gamma_cats)])
    if start_tree is not None:
        command_parts.extend(["--start-tree", str(start_tree)])
    if fix_tree is not None:
        command_parts.extend(["--fix-tree", str(fix_tree)])
    if resume is None:
        if chain_names:
            command_parts.extend(["--chain-names", chain_names])
        else:
            command_parts.extend(["--chains", str(chains), "--chain-prefix", chain_prefix])
    command_parts.extend(["--threads", str(threads), "--sample-freq", str(sample_freq), "--nsamples", str(nsamples)])
    command_parts.extend(["--monitor-freq", str(monitor_freq), "--burnin-frac", str(burnin_frac), "--poll-interval", str(poll_interval)])
    if resume is not None:
        command_parts.append("--resume")
        if resume != RESUME_ALL:
            command_parts.append(resume)
    if pb_path is not None:
        command_parts.extend(["--pb-path", str(pb_path)])
    if dry_run:
        command_parts.append("--dry-run")
    if quiet:
        command_parts.append("--quiet")
    command = " ".join(command_parts)

    chain_cmds: dict[str, list[str]] = {}
    skipped_names: list[str] = []
    for name in names:
        if resume is not None:
            current_len = _count_trace_samples(output_dir / "chains" / f"{name}.trace")
            if nsamples != -1 and current_len >= nsamples:
                skipped_names.append(name)
                continue
        if resume is None:
            chain_cmds[name] = _build_chain_cmd(tools["mpirun"], tools["pb_mpi"], threads, matrix_for_chains, model_flags, sample_freq, nsamples, name)
        else:
            chain_cmds[name] = _build_resume_cmd(tools["mpirun"], tools["pb_mpi"], threads, name)
    if not dry_run:
        # A new accepted run supersedes any stale stop record from a previous
        # interrupted attempt. Rejected runs, input-conversion failures and
        # chain-command validation failures never reach this point.
        _clear_stop_state(output_dir)
    if skipped_names:
        names = [n for n in names if n not in skipped_names]
        if not names:
            if not quiet:
                print("All requested chains already reached the NSAMPLES target; nothing to resume.")
            lengths = {n: _count_trace_samples(output_dir / "chains" / f"{n}.trace") for n in state.get("chain_names", [])}
            return _assemble_result(params, command, time.monotonic() - run_start, tool_versions, chain_cmds, lengths, None, {n: "" for n in state.get("chain_names", [])}, False, warnings=[f"Skipped {len(skipped_names)} chain(s) already at target: {', '.join(skipped_names)}"], output_dir=output_dir)
        if not quiet:
            msg = f"  Skipped {len(skipped_names)} chain(s) already at target: {', '.join(skipped_names)}"
            print(msg)
    if dry_run:
        return _assemble_result(params, command, time.monotonic() - run_start, {"pb_mpi": None, "bpcomp": None, "tracecomp": None, "mpirun": None}, chain_cmds, {name: 0 for name in names}, None, {name: "" for name in names}, False, output_dir=output_dir)
    return _run_bi_processes(output_dir, names, chain_cmds, tools, tool_versions, params, command, run_start, nsamples, monitor_freq, burnin_frac, poll_interval, quiet, stop=stop)
