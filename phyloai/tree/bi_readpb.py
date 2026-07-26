"""Posterior analysis with readpb_mpi."""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
import time
from pathlib import Path
from typing import Any

import numpy as np

from phyloai.core.env import ToolEnv

_VALID_MODES: set[str] = {
    "rr", "ss", "r", "sitelogl", "ppred",
    "div", "sitecomp", "siteconvprob", "comp", "allppred",
}

_MODE_FLAG_MAP: dict[str, str] = {
    "rr": "-rr", "ss": "-ss", "r": "-r", "sitelogl": "-sitelogl",
    "ppred": "-ppred", "div": "-div", "sitecomp": "-sitecomp",
    "siteconvprob": "-siteconvprob", "comp": "-comp", "allppred": "-allppred",
}

_PREDICTIVE_MODES: set[str] = {"div", "sitecomp", "siteconvprob", "comp"}

PAML_ORDER = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]
PHYLOBAYES_AA = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]
IQTREE_AA = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]


def _validate_modes(modes: list[str]) -> None:
    seen: set[str] = set()
    for mode in modes:
        if mode not in _VALID_MODES:
            raise ValueError(
                f"Unrecognised --mode value: {mode}. "
                f"Valid modes: {', '.join(sorted(_VALID_MODES))}"
            )
        if mode in seen:
            raise ValueError(f"Duplicate --mode value: {mode}")
        seen.add(mode)

    has_allppred = "allppred" in seen
    conflict = seen & _PREDICTIVE_MODES
    if has_allppred and conflict:
        raise ValueError(
            "--mode allppred is mutually exclusive with "
            f"{', '.join(sorted(conflict))}"
        )


def _convert_meanrr_to_exchangeabilities(meanrr_path: Path) -> Path:
    text = meanrr_path.read_text()
    lines = [l for l in text.splitlines() if l.strip()]
    if not lines:
        raise ValueError(f"Empty meanrr file: {meanrr_path}")

    order_of_aa = lines[0].split()
    if len(order_of_aa) != 20:
        raise ValueError(f"Expected 20 AA symbols in meanrr header, got {len(order_of_aa)}")

    aa_to_idx = {aa: i for i, aa in enumerate(order_of_aa)}
    exch = np.zeros((20, 20), dtype=np.float64)

    for line in lines[1:]:
        parts = line.split()
        if len(parts) < 3:
            continue
        source, target = parts[0], parts[1]
        try:
            val = float(parts[2])
        except ValueError:
            continue
        si = aa_to_idx[source]
        ti = aa_to_idx[target]
        exch[si, ti] = val
        exch[ti, si] = val

    paml_idx = [aa_to_idx[aa] for aa in PAML_ORDER]
    paml_exch = exch[np.ix_(paml_idx, paml_idx)]

    output_path = meanrr_path.with_suffix(".exchangeabilities")
    with open(output_path, "w") as fh:
        for i in range(20):
            row_vals = [f"{paml_exch[i, j]:08.6f}" for j in range(i)]
            fh.write(" ".join(row_vals))
            if row_vals:
                fh.write(" ")
            fh.write("\n")
        fh.write("\n")
        fh.write(" ".join(["0.050000"] * 20) + " \n")

    return output_path


def _convert_siteprofiles_to_sitefreq(siteprofiles_path: Path) -> Path:
    text = siteprofiles_path.read_text()
    lines = text.splitlines()
    if len(lines) < 3:
        raise ValueError(f"Too few lines in siteprofiles: {siteprofiles_path}")

    pb_order = PHYLOBAYES_AA
    iq_order = IQTREE_AA
    pb_to_iq = [pb_order.index(aa) for aa in iq_order]

    output_path = siteprofiles_path.with_suffix(".sitefreq")
    with open(output_path, "w") as fh:
        for line in lines[2:]:
            stripped = line.strip()
            if not stripped:
                continue
            parts = stripped.split()
            if len(parts) < 21:
                continue
            site_idx = parts[0]
            try:
                freqs = np.array([float(x) for x in parts[1:21]], dtype=np.float64)
            except ValueError:
                continue
            reindexed = freqs[pb_to_iq]
            reindexed = np.maximum(reindexed, 1e-8)
            reindexed = reindexed / reindexed.sum()
            freq_str = " ".join(f"{v:.8f}" for v in reindexed)
            fh.write(f"{site_idx} {freq_str}\n")

    return output_path


def run_bi_readpb(
    chain: Path,
    mode: str,
    output_dir: Path = Path("runs/tree/bi/readpb"),
    overwrite: bool = False,
    burnin: int = 0,
    sample_freq: int = 1,
    until: str = "all",
    threads: int = 4,
    pb_path: Path | None = None,
    dry_run: bool = False,
    quiet: bool = False,
) -> dict[str, Any]:
    start = time.monotonic()

    chain_path = Path(str(chain) + ".chain")
    if not chain_path.exists():
        raise ValueError(f"Chain file not found: {chain_path}")

    if burnin < 0:
        raise ValueError("--burnin must be >= 0")
    if sample_freq < 1:
        raise ValueError("--sample-freq must be >= 1")
    if threads < 2:
        raise ValueError("--threads must be at least 2")
    if until != "all":
        try:
            _until = int(until)
            if _until <= 0:
                raise ValueError("--until must be 'all' or a positive integer")
        except ValueError:
            raise ValueError("--until must be 'all' or a positive integer")

    modes = [m.strip() for m in mode.split(",") if m.strip()]
    if not modes:
        raise ValueError("--mode must contain at least one non-empty mode")
    _validate_modes(modes)

    if overwrite and output_dir.exists():
        shutil.rmtree(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    if dry_run:
        from phyloai.tree.bi import _build_x_flag
        cmds: dict[str, list[str]] = {}
        x_flag = _build_x_flag(burnin, sample_freq, until)
        work_dir = chain.parent
        chain_stem = chain.name
        for m in modes:
            cmds[m] = ["mpirun", "-np", str(threads), "readpb_mpi",
                        *x_flag, _MODE_FLAG_MAP[m], chain_stem]
        return {
            "status": "success",
            "command": "phyloai tree bi readpb",
            "wall_time": time.monotonic() - start,
            "tool_versions": {},
            "params": {
                "chain": str(chain),
                "mode": mode,
                "output_dir": str(output_dir),
                "overwrite": overwrite,
                "burnin": burnin,
                "sample_freq": sample_freq,
                "until": until,
                "threads": threads,
                "pb_path": str(pb_path) if pb_path else None,
                "dry_run": dry_run,
                "quiet": quiet,
            },
            "key_results": {"modes_run": modes, "output_files": {}},
            "error": None,
            "data": {"cmds": cmds, "post_processing": {}, "output_files": {}, "tool_stderr": "", "warnings": []},
        }

    if pb_path is not None:
        tool_paths = {"readpb_mpi": pb_path / "readpb_mpi", "mpirun": pb_path / "mpirun"}
        env = ToolEnv(tool_paths=tool_paths)
    else:
        env = ToolEnv()
    readpb_exe = str(env.require("readpb_mpi"))
    mpirun_exe = str(env.require("mpirun"))

    from phyloai.tree.bi import _build_x_flag, _detect_pb_version, _detect_mpirun_version
    readpb_ver = _detect_pb_version(readpb_exe)
    mpirun_ver = _detect_mpirun_version(mpirun_exe)

    x_flag = _build_x_flag(burnin, sample_freq, until)
    work_dir = chain.parent
    chain_stem = chain.name

    cmds: dict[str, list[str]] = {}
    all_stdout_parts: list[str] = []
    post_processing: dict[str, Any] = {}
    output_files: dict[str, str] = {}
    _pre_paths: dict[str, int] = {}
    _stat_warnings: list[str] = []
    for entry in work_dir.iterdir():
        try:
            _pre_paths[entry.name] = entry.stat().st_mtime_ns
        except OSError:
            _stat_warnings.append(f"cannot stat {entry.name}: skipped from pre-run snapshot")
            continue

    def _is_from_this_run(path: Path) -> bool:
        if path.name not in _pre_paths:
            return True
        try:
            return path.stat().st_mtime_ns != _pre_paths[path.name]
        except OSError:
            return False

    def _is_new_nonempty_file(path: Path) -> bool:
        try:
            return path.is_file() and _is_from_this_run(path) and path.stat().st_size > 0
        except OSError as exc:
            _stat_warnings.append(f"Could not inspect {path.name}: {exc}")
            return False

    for m in modes:
        cmd = [mpirun_exe, "-np", str(threads), readpb_exe, *x_flag, _MODE_FLAG_MAP[m], chain_stem]
        cmds[m] = cmd

        proc = subprocess.run(
            cmd,
            cwd=work_dir,
            stdout=subprocess.PIPE,
            stderr=None,
            text=True,
        )

        mode_stdout = proc.stdout or ""
        all_stdout_parts.append(f"--- {m} ---\n{mode_stdout}")
        if mode_stdout.strip():
            (output_dir / f"{m}.stdout").write_text(mode_stdout)
            if not quiet:
                sys.stdout.write(mode_stdout)

        if proc.returncode != 0:
            _data_of = {k: {"path": v, "description": f"readpb_mpi output: {k}"} for k, v in output_files.items()}
            return {
                "status": "error",
                "command": "phyloai tree bi readpb",
                "wall_time": time.monotonic() - start,
                "tool_versions": {},
                "params": {
                    "chain": str(chain),
                    "mode": mode,
                    "output_dir": str(output_dir),
                    "overwrite": overwrite,
                    "burnin": burnin,
                    "sample_freq": sample_freq,
                    "until": until,
                    "threads": threads,
                    "pb_path": str(pb_path) if pb_path else None,
                    "dry_run": dry_run,
                    "quiet": quiet,
                },
                "key_results": {"modes_run": modes, "output_files": output_files},
                "error": f"readpb_mpi --mode {m} exited with code {proc.returncode}",
                "data": {"cmds": cmds, "output_files": _data_of, "post_processing": post_processing, "tool_stderr": "\n".join(all_stdout_parts), "warnings": _stat_warnings},
            }

    pp_status: dict[str, Any] = {}
    _MODE_OUTPUT_PREFIX: dict[str, str] = {
        "rr": ".meanrr", "ss": ".siteprofiles", "r": ".meanrate",
        "sitelogl": ".sitelogl",
    }
    _MODE_OUTPUT_GLOB: set[str] = {
        "ppred", "div", "sitecomp", "siteconvprob", "comp",
    }
    _ALLPPRED_GLOB = ["div", "sitecomp", "siteconvprob", "comp"]
    for m in modes:
        if m == "allppred":
            missing_subs: list[str] = []
            for sub in _ALLPPRED_GLOB:
                found = [p for p in work_dir.glob(f"{chain_stem}.{sub}.*")
                         if _is_new_nonempty_file(p)]
                if not found:
                    missing_subs.append(sub)
                for p in found:
                    output_files[p.name.replace(".", "_")] = str(p)
            if missing_subs:
                pp_status[m] = {
                    "status": "error",
                    "error": f"tool exited 0 but no new non-empty output matching {chain_stem}.{{{','.join(missing_subs)}}}.* found",
                }
            else:
                pp_status[m] = {"status": "success"}
        elif m in _MODE_OUTPUT_GLOB:
            found = [p for p in work_dir.glob(f"{chain_stem}.{m}.*")
                     if _is_new_nonempty_file(p)]
            if not found:
                pp_status[m] = {
                    "status": "error",
                    "error": f"tool exited 0 but no new non-empty output matching {chain_stem}.{m}.* found",
                }
            else:
                pp_status[m] = {"status": "success"}
                for p in found:
                    output_files[p.name.replace(".", "_")] = str(p)
        else:
            suffix = _MODE_OUTPUT_PREFIX.get(m)
            if suffix is not None:
                expected = work_dir / f"{chain_stem}{suffix}"
                if not _is_new_nonempty_file(expected):
                    pp_status[m] = {
                        "status": "error",
                        "error": f"tool exited 0 but expected new non-empty output not found: {expected.name}",
                    }
                else:
                    output_files[expected.name.replace(".", "_")] = str(expected)

    if "rr" in modes:
        meanrr_path = work_dir / f"{chain_stem}.meanrr"
        if _is_new_nonempty_file(meanrr_path):
            try:
                exch_path = _convert_meanrr_to_exchangeabilities(meanrr_path)
                pp_status["rr"] = {
                    "input": f"{chain_stem}.meanrr",
                    "output": f"{chain_stem}.exchangeabilities",
                    "status": "success",
                }
                output_files["exchangeabilities"] = str(exch_path)
            except Exception as exc:
                pp_status["rr"] = {
                    "input": f"{chain_stem}.meanrr",
                    "output": f"{chain_stem}.exchangeabilities",
                    "status": "error",
                    "error": str(exc),
                }
        else:
            pp_status.setdefault("rr", {"status": "error", "error": f"new output not found: {chain_stem}.meanrr"})
        if meanrr_path.exists() and _is_from_this_run(meanrr_path):
            output_files["meanrr"] = str(meanrr_path)

    if "ss" in modes:
        siteprof_path = work_dir / f"{chain_stem}.siteprofiles"
        if _is_new_nonempty_file(siteprof_path):
            try:
                sitefreq_path = _convert_siteprofiles_to_sitefreq(siteprof_path)
                pp_status["ss"] = {
                    "input": f"{chain_stem}.siteprofiles",
                    "output": f"{chain_stem}.sitefreq",
                    "status": "success",
                }
                output_files["sitefreq"] = str(sitefreq_path)
            except Exception as exc:
                pp_status["ss"] = {
                    "input": f"{chain_stem}.siteprofiles",
                    "output": f"{chain_stem}.sitefreq",
                    "status": "error",
                    "error": str(exc),
                }
        else:
            pp_status.setdefault("ss", {"status": "error", "error": f"new output not found: {chain_stem}.siteprofiles"})
        if siteprof_path.exists() and _is_from_this_run(siteprof_path):
            output_files["siteprofiles"] = str(siteprof_path)

    pp_errors = [
        f"{m}: {info['error']}"
        for m, info in pp_status.items()
        if isinstance(info, dict) and info.get("status") == "error" and info.get("error")
    ]
    status = "error" if pp_errors else "success"
    error_msg = "; ".join(pp_errors) if pp_errors else None

    return {
        "status": status,
        "command": "phyloai tree bi readpb",
        "wall_time": time.monotonic() - start,
        "tool_versions": {"readpb_mpi": readpb_ver, "mpirun": mpirun_ver},
        "params": {
            "chain": str(chain),
            "mode": mode,
            "output_dir": str(output_dir),
            "overwrite": overwrite,
            "burnin": burnin,
            "sample_freq": sample_freq,
            "until": until,
            "threads": threads,
            "pb_path": str(pb_path) if pb_path else None,
            "dry_run": dry_run,
            "quiet": quiet,
        },
        "key_results": {
            "modes_run": modes,
            "output_files": output_files,
            "post_processing": pp_status,
        },
        "error": error_msg,
        "data": {
            "cmds": cmds,
            "post_processing": pp_status,
            "output_files": {k: {"path": v, "description": f"readpb_mpi output: {k}"} for k, v in output_files.items()},
            "tool_stderr": "\n".join(all_stdout_parts),
            "warnings": _stat_warnings,
        },
    }
