"""Posterior analysis with readpb_mpi."""

from __future__ import annotations

import os
import re
import shutil
import subprocess
import sys
import time
from math import isfinite
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
        for line_number, line in enumerate(lines[2:], start=3):
            stripped = line.strip()
            if not stripped:
                continue
            parts = stripped.split()
            if len(parts) != 21:
                raise ValueError(
                    f"Invalid siteprofiles row {line_number}: expected site ID plus 20 frequencies"
                )
            site_idx = parts[0]
            try:
                freqs = np.array([float(x) for x in parts[1:21]], dtype=np.float64)
            except ValueError as exc:
                raise ValueError(f"Invalid siteprofiles row {line_number}: invalid frequency") from exc
            if not np.isfinite(freqs).all():
                raise ValueError(f"Invalid siteprofiles row {line_number}: invalid frequency")
            reindexed = freqs[pb_to_iq]
            reindexed = np.maximum(reindexed, 1e-8)
            reindexed = reindexed / reindexed.sum()
            freq_str = " ".join(f"{v:.8f}" for v in reindexed)
            fh.write(f"{site_idx} {freq_str}\n")

    return output_path


def _write_pmsf_partition(
    sitefreq_path: Path,
    exchangeabilities_path: Path,
    meansiterates_path: Path,
    trace_path: Path,
    log_path: Path,
    burnin: int,
    sample_freq: int,
    until: str,
) -> tuple[Path, float, int]:
    frequencies: dict[int, list[str]] = {}
    for line_number, line in enumerate(sitefreq_path.read_text().splitlines(), start=1):
        parts = line.split()
        if not parts:
            continue
        if len(parts) != 21:
            raise ValueError(f"Invalid sitefreq row {line_number}: expected site ID plus 20 frequencies")
        try:
            site = int(parts[0])
        except ValueError as exc:
            raise ValueError(f"Invalid sitefreq site ID on row {line_number}: {parts[0]}") from exc
        if site <= 0:
            raise ValueError(f"Invalid sitefreq site ID on row {line_number}: {site}")
        if site in frequencies:
            raise ValueError(f"Duplicate sitefreq site ID: {site}")
        try:
            values = [float(value) for value in parts[1:]]
        except ValueError as exc:
            raise ValueError(f"Invalid sitefreq frequency on row {line_number}") from exc
        if not all(isfinite(value) for value in values):
            raise ValueError(f"Invalid sitefreq frequency on row {line_number}")
        frequencies[site] = parts[1:]
    if not frequencies:
        raise ValueError("sitefreq contains no sites")

    zero_based_rates: dict[int, str] = {}
    for line_number, line in enumerate(meansiterates_path.read_text().splitlines(), start=1):
        parts = line.split()
        if len(parts) != 2:
            raise ValueError(f"Invalid meansiterates row {line_number}")
        try:
            site = int(parts[0])
        except ValueError as exc:
            raise ValueError(f"Invalid meansiterates site index on row {line_number}") from exc
        if site < 0:
            raise ValueError(f"Invalid meansiterates site index on row {line_number}: {site}")
        if site in zero_based_rates:
            raise ValueError(f"Duplicate meansiterates site index: {site}")
        try:
            numeric_rate = float(parts[1])
        except ValueError as exc:
            raise ValueError(f"Invalid meansiterates rate on row {line_number}") from exc
        if not isfinite(numeric_rate) or numeric_rate <= 0:
            raise ValueError(f"Invalid meansiterates rate on row {line_number}")
        zero_based_rates[site] = parts[1]
    rates = {site + 1: rate for site, rate in zero_based_rates.items()}

    try:
        trace_lines = trace_path.read_text().splitlines()
        log_text = log_path.read_text()
    except OSError as exc:
        raise ValueError(f"Could not read PMSF metadata: {exc.filename}") from exc
    if not trace_lines:
        raise ValueError(f"Empty trace file: {trace_path}")
    trace_header = trace_lines[0].split()
    if "iter" not in trace_header or "alpha" not in trace_header:
        raise ValueError("Trace file must contain iter and alpha columns")
    iter_column = trace_header.index("iter")
    alpha_column = trace_header.index("alpha")
    until_value = None if until == "all" else int(until)
    alphas: list[float] = []
    for line_number, line in enumerate(trace_lines[1:], start=2):
        parts = line.split()
        if len(parts) != len(trace_header):
            raise ValueError(f"Invalid trace row {line_number}")
        try:
            iteration = int(parts[iter_column])
            alpha = float(parts[alpha_column])
        except ValueError as exc:
            raise ValueError(f"Invalid trace alpha on row {line_number}") from exc
        if not isfinite(alpha):
            raise ValueError(f"Invalid trace alpha on row {line_number}")
        if iteration < burnin or (iteration - burnin) % sample_freq:
            continue
        if until_value is not None and iteration > until_value:
            continue
        alphas.append(alpha)
    if not alphas:
        raise ValueError("No trace alpha samples remain after burnin/sample-freq/until filtering")
    mean_alpha = sum(alphas) / len(alphas)

    match = re.search(r"discrete gamma distribution of rates across sites \((\d+) categories\)", log_text)
    if not match:
        raise ValueError("Chain log does not declare discrete gamma categories")
    gamma_categories = int(match.group(1))

    if set(frequencies) != set(rates):
        raise ValueError("meansiterates and sitefreq site IDs do not match")

    output_path = sitefreq_path.parent / "partition.PMSF.nex"
    temporary_path = output_path.with_suffix(".nex.tmp")
    with open(temporary_path, "w") as fh:
        fh.write("#nexus\nbegin sets;\n")
        for site in sorted(frequencies):
            fh.write(f"    charset site_{site} = {site};\n")
        fh.write("    charpartition PMSF =\n")
        for index, site in enumerate(sorted(frequencies)):
            rate = rates[site]
            ending = ";" if index == len(frequencies) - 1 else ","
            freq_text = "/".join(frequencies[site])
            fh.write(
                f"        {exchangeabilities_path.name}+F{{{freq_text}}}+G{gamma_categories}{{{mean_alpha:.8f}}}:"
                f"site_{site}{{{rate}}}{ending}\n"
            )
        fh.write("end;\n")
    temporary_path.replace(output_path)
    return output_path, mean_alpha, gamma_categories


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
    _cmd_tokens = ["phyloai", "tree", "bi", "readpb",
                   "--chain", str(chain),
                   "--mode", mode,
                   "--output-dir", str(output_dir),
                   "--burnin", str(burnin),
                   "--sample-freq", str(sample_freq),
                   "--until", str(until),
                   "--threads", str(threads)]
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
        cmds: dict[str, list[str]] = {}
        x_flag = _build_x_flag(burnin, sample_freq, until)
        work_dir = chain.parent
        chain_stem = chain.name
        for m in modes:
            cmds[m] = ["mpirun", "-np", str(threads), "readpb_mpi",
                        *x_flag, _MODE_FLAG_MAP[m], chain_stem]
        return {
            "status": "success",
            "command": command_str,
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
    raw_output_files: dict[str, str] = {}
    partition_written = False
    make_pmsf_partition = {"ss", "rr", "r"}.issubset(modes)
    _MODE_OUTPUT_PREFIX: dict[str, str] = {
        "rr": ".meanrr", "ss": ".siteprofiles", "r": ".meansiterates",
        "sitelogl": ".sitelogl",
    }
    _MODE_OUTPUT_GLOB: set[str] = {
        "ppred", "div", "sitecomp", "siteconvprob", "comp",
    }

    _detection_threshold = time.time()

    def _file_exists(path: Path) -> bool:
        try:
            return path.is_file() and path.stat().st_size > 0 and path.stat().st_mtime >= _detection_threshold
        except OSError:
            return False

    def _move_mode_outputs(mode_name: str) -> str | None:
        found: list[Path] = []
        if mode_name == "allppred":
            found = [
                path for path in work_dir.glob(f"{chain_stem}.ppred")
                if _file_exists(path)
            ]
        elif mode_name == "ppred":
            found = [
                path for path in work_dir.glob(f"{chain_stem}*ppred*")
                if _file_exists(path)
            ]
        elif mode_name in _MODE_OUTPUT_GLOB:
            candidates = [work_dir / f"{chain_stem}.{mode_name}"]
            candidates.extend(work_dir.glob(f"{chain_stem}.{mode_name}.*"))
            found = [path for path in candidates if _file_exists(path)]
        else:
            suffix = _MODE_OUTPUT_PREFIX[mode_name]
            path = work_dir / f"{chain_stem}{suffix}"
            found = [path] if _file_exists(path) else []
            if mode_name == "sitelogl":
                cpo_path = work_dir / f"{chain_stem}.cpo"
                if _file_exists(cpo_path):
                    found.append(cpo_path)

        if not found:
            return f"tool exited 0 but no non-empty output found for --mode {mode_name}"

        try:
            for source in found:
                destination = output_dir / "ppred" if mode_name == "ppred" else output_dir
                destination.mkdir(parents=True, exist_ok=True)
                target = destination / source.name
                shutil.move(str(source), str(target))
                raw_output_files[source.name.replace(".", "_")] = str(target)
        except OSError as exc:
            return f"could not move --mode {mode_name} output: {exc}"

        if mode_name == "rr":
            meanrr_path = destination / f"{chain_stem}.meanrr"
            try:
                exchangeabilities = _convert_meanrr_to_exchangeabilities(meanrr_path)
                post_processing["rr"] = {
                    "input": meanrr_path.name,
                    "output": exchangeabilities.name,
                    "status": "success",
                }
                raw_output_files["exchangeabilities"] = str(exchangeabilities)
            except Exception as exc:
                return f"rr conversion failed: {exc}"
        elif mode_name == "ss":
            siteprofiles_path = destination / f"{chain_stem}.siteprofiles"
            try:
                sitefreq = _convert_siteprofiles_to_sitefreq(siteprofiles_path)
                post_processing["ss"] = {
                    "input": siteprofiles_path.name,
                    "output": sitefreq.name,
                    "status": "success",
                }
                raw_output_files["sitefreq"] = str(sitefreq)
            except Exception as exc:
                return f"ss conversion failed: {exc}"
        else:
            post_processing[mode_name] = {"status": "success"}

        return None

    def _write_partition_if_ready() -> str | None:
        nonlocal partition_written
        if not make_pmsf_partition or partition_written:
            return None
        sitefreq_path = output_dir / f"{chain_stem}.sitefreq"
        exchangeabilities_path = output_dir / f"{chain_stem}.exchangeabilities"
        meansiterates_path = output_dir / f"{chain_stem}.meansiterates"
        if not sitefreq_path.exists() or not exchangeabilities_path.exists() or not meansiterates_path.exists():
            return None
        try:
            partition, mean_alpha, gamma_categories = _write_pmsf_partition(
                sitefreq_path, exchangeabilities_path, meansiterates_path,
                work_dir / f"{chain_stem}.trace", work_dir / f"{chain_stem}.log",
                burnin, sample_freq, until,
            )
        except ValueError as exc:
            post_processing["pmsf_partition"] = {"status": "error", "error": str(exc)}
            return str(exc)
        partition_written = True
        raw_output_files["pmsf_partition"] = str(partition)
        post_processing["pmsf_partition"] = {
            "inputs": [meansiterates_path.name, f"{chain_stem}.trace", f"{chain_stem}.log"],
            "output": partition.name, "mean_alpha": mean_alpha,
            "gamma_categories": gamma_categories, "status": "success",
        }
        return None

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
        if mode_stdout:
            all_stdout_parts.append(f"--- {m} ---\n{mode_stdout}")
        if mode_stdout.strip():
            (output_dir / f"{m}.stdout").write_text(mode_stdout)
            if not quiet:
                sys.stdout.write(mode_stdout)

        if proc.returncode != 0:
            _data_of = {k: {"path": v, "description": f"readpb_mpi output: {k}"} for k, v in raw_output_files.items()}
            return {
                "status": "error",
                "command": command_str,
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
                "key_results": {"modes_run": modes, "output_files": raw_output_files},
                "error": f"readpb_mpi --mode {m} exited with code {proc.returncode}",
                "data": {"cmds": cmds, "output_files": _data_of, "post_processing": post_processing, "tool_stderr": "\n".join(all_stdout_parts), "warnings": []},
            }

        move_error = _move_mode_outputs(m)
        if move_error:
            post_processing[m] = {"status": "error", "error": move_error}
            return {
                "status": "error",
                "command": command_str,
                "wall_time": time.monotonic() - start,
                "tool_versions": {"readpb_mpi": readpb_ver, "mpirun": mpirun_ver},
                "params": {
                    "chain": str(chain), "mode": mode, "output_dir": str(output_dir),
                    "overwrite": overwrite, "burnin": burnin, "sample_freq": sample_freq,
                    "until": until, "threads": threads,
                    "pb_path": str(pb_path) if pb_path else None,
                    "dry_run": dry_run, "quiet": quiet,
                },
                "key_results": {"modes_run": modes, "output_files": raw_output_files},
                "error": move_error,
                "data": {
                    "cmds": cmds, "output_files": {
                        key: {"path": path, "description": f"readpb_mpi output: {key}"}
                        for key, path in raw_output_files.items()
                    },
                    "post_processing": post_processing,
                    "tool_stderr": "\n".join(all_stdout_parts), "warnings": [],
                },
            }

        partition_error = _write_partition_if_ready()
        if partition_error:
            return {
                "status": "error", "command": command_str,
                "wall_time": time.monotonic() - start,
                "tool_versions": {"readpb_mpi": readpb_ver, "mpirun": mpirun_ver},
                "params": {"chain": str(chain), "mode": mode, "output_dir": str(output_dir),
                           "overwrite": overwrite, "burnin": burnin, "sample_freq": sample_freq,
                           "until": until, "threads": threads,
                           "pb_path": str(pb_path) if pb_path else None, "dry_run": dry_run, "quiet": quiet},
                "key_results": {"modes_run": modes, "output_files": raw_output_files},
                "error": partition_error,
                "data": {"cmds": cmds, "post_processing": post_processing,
                         "output_files": {k: {"path": v, "description": f"readpb_mpi output: {k}"} for k, v in raw_output_files.items()},
                         "tool_stderr": "\n".join(all_stdout_parts), "warnings": []},
            }

    if not quiet and raw_output_files:
        print("PhyloAI: output files:")
        for key in sorted(raw_output_files):
            print(f"  {key}: {raw_output_files[key]}")

    return {
        "status": "success",
        "command": command_str,
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
            "output_files": raw_output_files,
            "post_processing": post_processing,
        },
        "error": None,
        "data": {
            "cmds": cmds,
            "post_processing": post_processing,
            "output_files": {k: {"path": v, "description": f"readpb_mpi output: {k}"} for k, v in raw_output_files.items()},
            "tool_stderr": "\n".join(all_stdout_parts),
            "warnings": [],
        },
    }
