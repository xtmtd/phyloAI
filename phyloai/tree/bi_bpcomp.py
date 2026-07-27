"""Topology convergence analysis with bpcomp."""

from __future__ import annotations

import os
import shutil
import subprocess
import time
from pathlib import Path
from typing import Any

from phyloai.core.env import ToolEnv


def _discover_chain_names(chain_dir: Path) -> list[str]:
    names = []
    for entry in sorted(chain_dir.iterdir()):
        if entry.is_file() and entry.suffix == ".chain":
            names.append(entry.stem)
    if not names:
        raise FileNotFoundError(f"No .chain files found in {chain_dir}")
    return names


def _validate_chain_names(chain_dir: Path, names: list[str]) -> None:
    missing = []
    for name in names:
        if not (chain_dir / f"{name}.chain").exists():
            missing.append(name)
    if missing:
        raise FileNotFoundError(
            f"Chain file(s) not found in {chain_dir}: {', '.join(missing)}"
        )


def run_bi_bpcomp(
    chain_dir: Path,
    chain_names: str = "all",
    output_dir: Path = Path("runs/tree/bi/bpcomp"),
    overwrite: bool = False,
    burnin: int = 0,
    sample_freq: int = 1,
    until: str = "all",
    cutoff: float = 0.5,
    pb_path: Path | None = None,
    dry_run: bool = False,
    quiet: bool = False,
) -> dict[str, Any]:
    start = time.monotonic()

    if burnin < 0:
        raise ValueError("--burnin must be >= 0")
    if sample_freq < 1:
        raise ValueError("--sample-freq must be >= 1")
    if not (0.0 < cutoff < 1.0):
        raise ValueError("--cutoff must be strictly between 0 and 1")
    if until != "all":
        try:
            _until = int(until)
            if _until <= 0:
                raise ValueError("--until must be 'all' or a positive integer")
        except ValueError:
            raise ValueError("--until must be 'all' or a positive integer")

    if chain_names == "all":
        resolved = _discover_chain_names(chain_dir)
    else:
        resolved = [n.strip() for n in chain_names.split(",") if n.strip()]
        if not resolved:
            raise ValueError("--chain-names must contain at least one non-empty name")
        _validate_chain_names(chain_dir, resolved)

    if not dry_run:
        if not overwrite and output_dir.exists() and any(output_dir.iterdir()):
            raise ValueError(
                f"Output directory {output_dir} already exists and is non-empty. "
                "Use --overwrite to replace."
            )
        if overwrite and output_dir.exists():
            shutil.rmtree(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

    # --- build full command string ---
    _cmd_tokens = ["phyloai", "tree", "bi", "bpcomp",
                   "--chain-dir", str(chain_dir),
                   "--chain-names", chain_names,
                   "--output-dir", str(output_dir),
                   "--burnin", str(burnin),
                   "--sample-freq", str(sample_freq),
                   "--until", str(until),
                   "--cutoff", str(cutoff)]
    if pb_path is not None:
        _cmd_tokens.extend(["--pb-path", str(pb_path)])
    if overwrite:
        _cmd_tokens.append("--overwrite")
    if dry_run:
        _cmd_tokens.append("--dry-run")
    if quiet:
        _cmd_tokens.append("--quiet")
    command_str = " ".join(_cmd_tokens)

    if dry_run:
        from phyloai.tree.bi import _build_x_flag
        cmd = ["bpcomp", *_build_x_flag(burnin, sample_freq, until),
               "-c", str(cutoff), "-o", "bpcomp"]
        for name in resolved:
            cmd.append(os.path.relpath(chain_dir / name, output_dir))
        return {
            "status": "success",
            "command": command_str,
            "wall_time": time.monotonic() - start,
            "tool_versions": {},
            "params": {
                "chain_dir": str(chain_dir),
                "chain_names": chain_names,
                "output_dir": str(output_dir),
                "overwrite": overwrite,
                "burnin": burnin,
                "sample_freq": sample_freq,
                "until": until,
                "cutoff": cutoff,
                "pb_path": str(pb_path) if pb_path else None,
                "dry_run": dry_run,
                "quiet": quiet,
            },
            "key_results": {
                "chains_used": resolved,
                "bpcomp_maxdiff": None,
                "bpcomp_meandiff": None,
                "bpcomp_status": None,
                "consensus_tree": None,
            },
            "error": None,
            "data": {
                "cmd": cmd,
                "output_files": {},
                "tool_stderr": "",
            },
        }

    if pb_path is not None:
        env = ToolEnv(tool_paths={"bpcomp": pb_path / "bpcomp"})
    else:
        env = ToolEnv()
    bpcomp_exe = str(env.require("bpcomp"))

    from phyloai.tree.bi import _build_x_flag, _detect_pb_version
    bpcomp_ver = _detect_pb_version(bpcomp_exe)

    cmd = [bpcomp_exe, *_build_x_flag(burnin, sample_freq, until),
           "-c", str(cutoff), "-o", "bpcomp"]
    for name in resolved:
        cmd.append(os.path.relpath(chain_dir / name, output_dir))

    if not quiet:
        print(f"PhyloAI: executing {' '.join(cmd)}")

    proc = subprocess.run(cmd, cwd=output_dir, capture_output=False)
    if proc.returncode != 0:
        return {
            "status": "error",
            "command": command_str,
            "wall_time": time.monotonic() - start,
            "tool_versions": {"bpcomp": bpcomp_ver},
            "params": {
                "chain_dir": str(chain_dir),
                "chain_names": chain_names,
                "output_dir": str(output_dir),
                "overwrite": overwrite,
                "burnin": burnin,
                "sample_freq": sample_freq,
                "until": until,
                "cutoff": cutoff,
                "pb_path": str(pb_path) if pb_path else None,
                "dry_run": dry_run,
                "quiet": quiet,
            },
            "key_results": {"chains_used": resolved},
            "error": f"bpcomp exited with code {proc.returncode}",
            "data": {"cmd": cmd, "output_files": {}, "tool_stderr": ""},
        }

    from phyloai.tree.bi import _parse_bpcomp_bpdiff, _bpcomp_status

    bpdiff_path = output_dir / "bpcomp.bpdiff"
    bplist_path = output_dir / "bpcomp.bplist"
    contree_path = output_dir / "bpcomp.con.tre"

    if not (bpdiff_path.exists() and bplist_path.exists() and contree_path.exists()):
        missing = []
        for name, p in [("bpcomp.bpdiff", bpdiff_path), ("bpcomp.bplist", bplist_path), ("bpcomp.con.tre", contree_path)]:
            if not p.exists():
                missing.append(name)
        return {
            "status": "error",
            "command": command_str,
            "wall_time": time.monotonic() - start,
            "tool_versions": {"bpcomp": bpcomp_ver},
            "params": {
                "chain_dir": str(chain_dir),
                "chain_names": chain_names,
                "output_dir": str(output_dir),
                "overwrite": overwrite,
                "burnin": burnin,
                "sample_freq": sample_freq,
                "until": until,
                "cutoff": cutoff,
                "pb_path": str(pb_path) if pb_path else None,
                "dry_run": dry_run,
                "quiet": quiet,
            },
            "key_results": {"chains_used": resolved},
            "error": f"bpcomp exited 0 but expected output files not found: {', '.join(missing)}",
            "data": {"cmd": cmd, "output_files": {}, "tool_stderr": ""},
        }

    bpdiff = _parse_bpcomp_bpdiff(bpdiff_path)
    maxdiff = bpdiff["maxdiff"]
    meandiff = bpdiff["meandiff"]
    status = _bpcomp_status(maxdiff)

    consensus_tree = str(contree_path) if contree_path.exists() else None
    if not quiet:
        from phyloai.tree.bi import _bpcomp_status as _st
        md_str = f"{maxdiff:.3f}" if maxdiff is not None else "--"
        mn_str = f"{meandiff:.3f}" if meandiff is not None else "--"
        st = _st(maxdiff)
        print(f"PhyloAI: maxdiff {md_str}  meandiff {mn_str}  [{st}]  -> bpcomp/bpcomp.con.tre")

    output_files = {}
    for name, p in [("bpdiff", bpdiff_path), ("bplist", bplist_path), ("consensus_tree", contree_path)]:
        if p.exists():
            output_files[name] = {"path": str(p), "description": f"bpcomp output: {name}"}

    return {
        "status": "success",
        "command": command_str,
        "wall_time": time.monotonic() - start,
        "tool_versions": {"bpcomp": bpcomp_ver},
        "params": {
            "chain_dir": str(chain_dir),
            "chain_names": chain_names,
            "output_dir": str(output_dir),
            "overwrite": overwrite,
            "burnin": burnin,
            "sample_freq": sample_freq,
            "until": until,
            "cutoff": cutoff,
            "pb_path": str(pb_path) if pb_path else None,
            "dry_run": dry_run,
            "quiet": quiet,
        },
        "key_results": {
            "chains_used": resolved,
            "bpcomp_maxdiff": maxdiff,
            "bpcomp_meandiff": meandiff,
            "bpcomp_status": status,
            "consensus_tree": consensus_tree,
        },
        "error": None,
        "data": {
            "cmd": cmd,
            "output_files": output_files,
            "tool_stderr": "",
        },
    }
