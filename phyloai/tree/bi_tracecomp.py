"""Parameter convergence analysis with tracecomp."""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
import time
from pathlib import Path
from typing import Any

from phyloai.core.env import ToolEnv


def _discover_trace_names(chain_dir: Path) -> list[str]:
    names = []
    for entry in sorted(chain_dir.iterdir()):
        if entry.is_file() and entry.suffix == ".trace":
            names.append(entry.stem)
    if not names:
        raise FileNotFoundError(f"No .trace files found in {chain_dir}")
    return names


def _validate_trace_names(chain_dir: Path, names: list[str]) -> None:
    missing = []
    for name in names:
        if not (chain_dir / f"{name}.trace").exists():
            missing.append(name)
    if missing:
        raise FileNotFoundError(
            f"Trace file(s) not found in {chain_dir}: {', '.join(missing)}"
        )


def _annotate_tracecomp_output(stdout: str) -> tuple[str, float | None, float | None]:
    from phyloai.tree.bi import _tracecomp_status

    lines = stdout.strip().splitlines()
    if not lines:
        return ("", None, None)

    header = lines[0]
    annotated_lines = [header + "\tstatus"]

    min_effsize: float | None = None
    max_rel_diff: float | None = None

    for line in lines[1:]:
        stripped = line.strip()
        if not stripped:
            continue
        parts = stripped.split()
        if len(parts) < 3:
            annotated_lines.append(stripped)
            continue
        try:
            effsize = float(parts[-2])
            rel_diff = float(parts[-1])
        except (ValueError, IndexError):
            annotated_lines.append(stripped)
            continue

        if min_effsize is None:
            min_effsize = effsize
        else:
            min_effsize = min(min_effsize, effsize)

        if max_rel_diff is None:
            max_rel_diff = rel_diff
        else:
            max_rel_diff = max(max_rel_diff, rel_diff)

        status = _tracecomp_status(effsize, rel_diff)
        annotated_lines.append(f"{stripped}\t[{status}]")

    return ("\n".join(annotated_lines), min_effsize, max_rel_diff)


def run_bi_tracecomp(
    chain_dir: Path = Path("runs/tree/bi/chains"),
    chain_names: str = "all",
    output_dir: Path = Path("runs/tree/bi/tracecomp"),
    overwrite: bool = False,
    burnin: int = 0,
    pb_path: Path | None = None,
    dry_run: bool = False,
    quiet: bool = False,
) -> dict[str, Any]:
    start = time.monotonic()

    if burnin < 0:
        raise ValueError("--burnin must be >= 0")

    if chain_names == "all":
        resolved = _discover_trace_names(chain_dir)
    else:
        resolved = [n.strip() for n in chain_names.split(",") if n.strip()]
        if not resolved:
            raise ValueError("--chain-names must contain at least one non-empty name")
        _validate_trace_names(chain_dir, resolved)

    if overwrite and output_dir.exists():
        shutil.rmtree(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    if dry_run:
        cmd = ["tracecomp", "-x", str(burnin)]
        for name in resolved:
            cmd.append(os.path.relpath(chain_dir / f"{name}.trace", output_dir))
        return {
            "status": "success",
            "command": "phyloai tree bi tracecomp",
            "wall_time": time.monotonic() - start,
            "tool_versions": {},
            "params": {
                "chain_dir": str(chain_dir),
                "chain_names": chain_names,
                "output_dir": str(output_dir),
                "overwrite": overwrite,
                "burnin": burnin,
                "pb_path": str(pb_path) if pb_path else None,
                "dry_run": dry_run,
                "quiet": quiet,
            },
            "key_results": {
                "chains_used": resolved,
                "tracecomp_min_effsize": None,
                "tracecomp_max_reldiff": None,
                "tracecomp_status": None,
            },
            "error": None,
            "data": {
                "cmd": cmd,
                "output_files": {},
                "tool_stderr": "",
            },
        }

    if pb_path is not None:
        env = ToolEnv(tool_paths={"tracecomp": pb_path / "tracecomp"})
    else:
        env = ToolEnv()
    tracecomp_exe = str(env.require("tracecomp"))

    from phyloai.tree.bi import _detect_pb_version
    tracecomp_ver = _detect_pb_version(tracecomp_exe)

    cmd = [tracecomp_exe, "-x", str(burnin)]
    for name in resolved:
        cmd.append(os.path.relpath(chain_dir / f"{name}.trace", output_dir))

    if not quiet:
        print(f"PhyloAI: executing {' '.join(cmd)}")

    proc = subprocess.run(
        cmd,
        cwd=output_dir,
        stdout=subprocess.PIPE,
        stderr=None,
        text=True,
    )
    if proc.returncode != 0:
        return {
            "status": "error",
            "command": "phyloai tree bi tracecomp",
            "wall_time": time.monotonic() - start,
            "tool_versions": {},
            "params": {
                "chain_dir": str(chain_dir),
                "chain_names": chain_names,
                "output_dir": str(output_dir),
                "overwrite": overwrite,
                "burnin": burnin,
                "pb_path": str(pb_path) if pb_path else None,
                "dry_run": dry_run,
                "quiet": quiet,
            },
            "key_results": {"chains_used": resolved},
            "error": f"tracecomp exited with code {proc.returncode}",
            "data": {"cmd": cmd, "output_files": {}, "tool_stderr": proc.stdout or ""},
        }

    stdout = proc.stdout or ""
    (output_dir / "tracecomp.contdiff").write_text(stdout)

    annotated, min_effsize, max_rel_diff = _annotate_tracecomp_output(stdout)

    from phyloai.tree.bi import _tracecomp_status
    tc_status = _tracecomp_status(min_effsize, max_rel_diff)

    if not quiet:
        print(annotated)
        me_str = f"{min_effsize:.0f}" if min_effsize is not None else "--"
        mr_str = f"{max_rel_diff:.7f}" if max_rel_diff is not None else "--"
        print(f"PhyloAI: min effsize {me_str}  max rel_diff {mr_str}  [{tc_status}]")

    output_files = {}
    contdiff_path = output_dir / "tracecomp.contdiff"
    if contdiff_path.exists():
        output_files["contdiff"] = {"path": str(contdiff_path), "description": "tracecomp output"}

    return {
        "status": "success",
        "command": "phyloai tree bi tracecomp",
        "wall_time": time.monotonic() - start,
        "tool_versions": {"tracecomp": tracecomp_ver},
        "params": {
            "chain_dir": str(chain_dir),
            "chain_names": chain_names,
            "output_dir": str(output_dir),
            "overwrite": overwrite,
            "burnin": burnin,
            "pb_path": str(pb_path) if pb_path else None,
            "dry_run": dry_run,
            "quiet": quiet,
        },
        "key_results": {
            "chains_used": resolved,
            "tracecomp_min_effsize": min_effsize,
            "tracecomp_max_reldiff": max_rel_diff,
            "tracecomp_status": tc_status,
        },
        "error": None,
        "data": {
            "cmd": cmd,
            "output_files": output_files,
            "tool_stderr": stdout,
        },
    }
