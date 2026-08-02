"""Model comparison via PhyloBayes LOO-CV / wAIC."""
from __future__ import annotations

import csv
import json
import math
import shlex
import shutil
import time as _time
from pathlib import Path
from typing import Any

# Student's t 95% critical values for df 1..30 (from standard tables)
# For df > 30: linear interpolation toward 1.96 (z_0.975).
# Formula: t30 + (1.96 - t30) * (1 - 30/df) — monotonically decreases from t30=2.042 to 1.96.
_STUDENT_T95 = {
    1: 12.706, 2: 4.303, 3: 3.182, 4: 2.776, 5: 2.571,
    6: 2.447, 7: 2.365, 8: 2.306, 9: 2.262, 10: 2.228,
    11: 2.201, 12: 2.179, 13: 2.160, 14: 2.145, 15: 2.131,
    16: 2.120, 17: 2.110, 18: 2.101, 19: 2.093, 20: 2.086,
    21: 2.080, 22: 2.074, 23: 2.069, 24: 2.064, 25: 2.060,
    26: 2.056, 27: 2.052, 28: 2.048, 29: 2.045, 30: 2.042,
}

_MIN_ESS = 10.0

_LOGL_IDX = 1
_VAR_IDX = 2
_LOGCPO_IDX = 3
_ESS_LOGCPO_IDX = 4
_LOGPOSTMEANL_IDX = 5
_ESS_POSTMEANL_IDX = 6


def _mean(x: list[float]) -> float:
    if not x:
        return 0.0
    return sum(x) / len(x)


def _unbiased_var(x: list[float]) -> float:
    n = len(x)
    if n < 2:
        return 0.0
    m = sum(x) / n
    return sum((v - m) ** 2 for v in x) / (n - 1)


def _student_t95(df: int) -> float:
    if df <= 1:
        return _STUDENT_T95[1]
    if df < 30:
        return _STUDENT_T95[df]
    if df == 30:
        return _STUDENT_T95[30]
    return _STUDENT_T95[30] + (1.96 - _STUDENT_T95[30]) * (1 - 30 / df)


def _parse_sitelogl(path: Path) -> list[list[float]]:
    """Parse one .sitelogl file (skip header) into rows of 7 floats.

    Rejects duplicate site identifiers within a single file (duplicate sites
    would double-count and distort LOO-CV/wAIC scores).
    """
    rows: list[list[float]] = []
    seen_sites: set[float] = set()
    with open(path) as fh:
        fh.readline()  # header: site logl var logcpo ess logpostmeanl ess
        for line in fh:
            line = line.strip()
            if not line:
                continue
            parts = line.split()
            if len(parts) != 7:
                raise ValueError(
                    f"Malformed sitelogl line in {path}: expected 7 columns, got {len(parts)}"
                )
            try:
                row = [float(p) for p in parts]
            except ValueError:
                raise ValueError(f"Malformed numeric value in sitelogl line in {path}: {line}")
            if row[0] in seen_sites:
                raise ValueError(
                    f"Duplicate site identifier {row[0]:g} in sitelogl file: {path}"
                )
            seen_sites.add(row[0])
            rows.append(row)
    if not rows:
        raise ValueError(f"sitelogl file has no data rows: {path}")
    return rows


def _compute_loocv_waic(runs: list[list[list[float]]]) -> dict[str, Any]:
    """Integrated readwaic() logic (Lartillot 2023) for a single model group.

    `runs` is a list of chains; each chain is a list of site rows, each row
    is ``[site, logl, var, logcpo, ess, logpostmeanl, ess]``.
    """
    nrun = len(runs)
    m = len(runs[0])

    # ---- LOO-CV ----
    loocvs = [_mean([row[_LOGCPO_IDX] for row in run]) for run in runs]
    loocv = _mean(loocvs)
    var_loocv = _unbiased_var(loocvs)
    across_vars_loocv = [_unbiased_var([run[i][_LOGCPO_IDX] for run in runs]) for i in range(m)]
    bias_loocv = 0.5 * _mean(across_vars_loocv)
    debiased_loocv = loocv - bias_loocv
    ess_loocv = [_mean([run[i][_ESS_LOGCPO_IDX] for run in runs]) for i in range(m)]
    mean_ess_loocv = _mean(ess_loocv)
    pct_ess_lt10_loocv = sum(1 for s in ess_loocv if s < _MIN_ESS) / m
    site_logcpo_mean = [_mean([run[i][_LOGCPO_IDX] for run in runs]) for i in range(m)]
    frac_ess_lt10_loocv = (
        abs(_mean([(ess_loocv[i] < _MIN_ESS) * site_logcpo_mean[i] for i in range(m)]) / loocv)
        if loocv != 0 else 0.0
    )

    # ---- wAIC (Watanabe) ----
    logpostmeanls = [_mean([row[_LOGPOSTMEANL_IDX] for row in run]) for run in runs]
    postvar_logls = [_mean([row[_VAR_IDX] for row in run]) for run in runs]
    waics = [logpostmeanls[i] - postvar_logls[i] for i in range(nrun)]
    waic = _mean(waics)
    var_waic = _unbiased_var(waics)
    across_vars_waic = [_unbiased_var([run[i][_LOGPOSTMEANL_IDX] for run in runs]) for i in range(m)]
    bias_waic = -0.5 * _mean(across_vars_waic)
    debiased_waic = waic - bias_waic
    ess_waic = [_mean([run[i][_ESS_POSTMEANL_IDX] for run in runs]) for i in range(m)]
    mean_ess_waic = _mean(ess_waic)
    pct_ess_lt10_waic = sum(1 for s in ess_waic if s < _MIN_ESS) / m
    sitewaics = [_mean([run[i][_LOGPOSTMEANL_IDX] - run[i][_VAR_IDX] for run in runs]) for i in range(m)]
    frac_ess_lt10_waic = (
        abs(_mean([(ess_waic[i] < _MIN_ESS) * sitewaics[i] for i in range(m)]) / waic)
        if waic != 0 else 0.0
    )

    t = _student_t95(nrun - 1)
    stdev_loocv = math.sqrt(var_loocv)
    stdev_waic = math.sqrt(var_waic)

    return {
        "loocv": {
            "score": debiased_loocv,
            "bias": bias_loocv,
            "stdev": stdev_loocv,
            "ci95_min": debiased_loocv - t * stdev_loocv,
            "ci95_max": debiased_loocv + t * stdev_loocv,
            "ess": mean_ess_loocv,
            "pct_ess_lt10": pct_ess_lt10_loocv,
            "frac_ess_lt10": frac_ess_lt10_loocv,
            "quality": _classify_quality(pct_ess_lt10_loocv, frac_ess_lt10_loocv),
        },
        "waic": {
            "score": debiased_waic,
            "bias": bias_waic,
            "stdev": stdev_waic,
            "ci95_min": debiased_waic - t * stdev_waic,
            "ci95_max": debiased_waic + t * stdev_waic,
            "ess": mean_ess_waic,
            "pct_ess_lt10": pct_ess_lt10_waic,
            "frac_ess_lt10": frac_ess_lt10_waic,
            "quality": _classify_quality(pct_ess_lt10_waic, frac_ess_lt10_waic),
        },
    }


def _classify_quality(pct_ess_lt10: float, frac_ess_lt10: float) -> str:
    """good / ok / no based on max(%(ess<10), f(ess<10)) thresholds."""
    worst = max(pct_ess_lt10, frac_ess_lt10)
    if worst < 0.1:
        return "good"
    if worst < 0.3:
        return "ok"
    return "no"


def _groups_from_dirs(sitelogl_dir: str) -> list[list[Path]]:
    dirs = [d.strip() for d in sitelogl_dir.split(",") if d.strip()]
    if not dirs:
        raise ValueError("--sitelogl-dir must contain at least one directory.")
    groups: list[list[Path]] = []
    for d in dirs:
        p = Path(d)
        if not p.is_dir():
            raise ValueError(f"--sitelogl-dir entry is not a directory: {p}")
        files = sorted(
            f for f in p.iterdir()
            if f.is_file() and f.suffix.lower() == ".sitelogl"
        )
        if len(files) < 2:
            raise ValueError(
                f"Directory '{p}' must contain at least 2 *.sitelogl files, found {len(files)}."
            )
        groups.append(files)
    return groups


def _groups_from_files(sitelogl: list[str]) -> list[list[Path]]:
    groups: list[list[Path]] = []
    for gi, seg in enumerate(sitelogl, 1):
        files = [Path(p.strip()) for p in seg.split(",") if p.strip()]
        if not files:
            raise ValueError(f"--sitelogl group #{gi} is empty.")
        for f in files:
            if not f.is_file():
                raise ValueError(f"--sitelogl group #{gi} file does not exist: {f}")
        if len(files) < 2:
            raise ValueError(
                f"--sitelogl group #{gi} must contain at least 2 .sitelogl files, found {len(files)}."
            )
        groups.append(files)
    return groups


def _is_safe_name(name: str) -> bool:
    """True if `name` is a single path component safe for use as a subdirectory.

    Rejects path separators, absolute paths, and ``.``/``..`` so that a
    ``--model-names`` label can never escape the output directory.
    """
    if not name or name in (".", ".."):
        return False
    return Path(name).name == name


def _dedupe_basename(basename: str, used: set[str]) -> str:
    if basename not in used:
        used.add(basename)
        return basename
    stem = Path(basename).stem
    suffix = Path(basename).suffix
    counter = 1
    candidate = f"{stem}_{counter}{suffix}"
    while candidate in used:
        counter += 1
        candidate = f"{stem}_{counter}{suffix}"
    used.add(candidate)
    return candidate


def _write_model_fit_csv(path: Path, model_results: list[dict[str, Any]]) -> None:
    if len(model_results) == 1:
        _write_single_model_csv(path, model_results[0])
    else:
        _write_multi_model_csv(path, model_results)


def _write_single_model_csv(path: Path, mr: dict[str, Any]) -> None:
    fieldnames = [
        "Metric", "Score", "Bias", "StDev", "CI95_min", "CI95_max",
        "ESS", "Pct_ESS_lt10", "Frac_ESS_lt10", "Quality",
    ]
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        for metric, label in (("loocv", "LOO-CV"), ("waic", "wAIC")):
            m = mr[metric]
            writer.writerow({
                "Metric": label,
                "Score": m["score"],
                "Bias": m["bias"],
                "StDev": m["stdev"],
                "CI95_min": m["ci95_min"],
                "CI95_max": m["ci95_max"],
                "ESS": m["ess"],
                "Pct_ESS_lt10": m["pct_ess_lt10"],
                "Frac_ESS_lt10": m["frac_ess_lt10"],
                "Quality": m["quality"],
            })


def _write_multi_model_csv(path: Path, model_results: list[dict[str, Any]]) -> None:
    fieldnames = [
        "Model", "LOO-CV", "LOO-CV_Bias", "LOO-CV_StDev", "LOO-CV_CI95min",
        "LOO-CV_CI95max", "LOO-CV_ESS", "LOO-CV_Pct_ESS_lt10",
        "LOO-CV_Frac_ESS_lt10", "LOO-CV_Quality",
        "wAIC", "wAIC_Bias", "wAIC_StDev", "wAIC_CI95min", "wAIC_CI95max",
        "wAIC_ESS", "wAIC_Pct_ESS_lt10", "wAIC_Frac_ESS_lt10", "wAIC_Quality",
        "Delta_LOOCV", "Delta_wAIC",
    ]
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        for mr in model_results:
            loocv, waic = mr["loocv"], mr["waic"]
            writer.writerow({
                "Model": mr["model"],
                "LOO-CV": loocv["score"],
                "LOO-CV_Bias": loocv["bias"],
                "LOO-CV_StDev": loocv["stdev"],
                "LOO-CV_CI95min": loocv["ci95_min"],
                "LOO-CV_CI95max": loocv["ci95_max"],
                "LOO-CV_ESS": loocv["ess"],
                "LOO-CV_Pct_ESS_lt10": loocv["pct_ess_lt10"],
                "LOO-CV_Frac_ESS_lt10": loocv["frac_ess_lt10"],
                "LOO-CV_Quality": loocv["quality"],
                "wAIC": waic["score"],
                "wAIC_Bias": waic["bias"],
                "wAIC_StDev": waic["stdev"],
                "wAIC_CI95min": waic["ci95_min"],
                "wAIC_CI95max": waic["ci95_max"],
                "wAIC_ESS": waic["ess"],
                "wAIC_Pct_ESS_lt10": waic["pct_ess_lt10"],
                "wAIC_Frac_ESS_lt10": waic["frac_ess_lt10"],
                "wAIC_Quality": waic["quality"],
                "Delta_LOOCV": mr["delta_loocv"],
                "Delta_wAIC": mr["delta_waic"],
            })


def _write_result_json(output_dir: Path, result: dict[str, Any]) -> None:
    (output_dir / "result.json").write_text(json.dumps(result, indent=2))


def run_modelcompare_pb(
    *,
    sitelogl_dir: str | None = None,
    sitelogl: list[str] | None = None,
    model_names: str | None = None,
    output_dir: Path | None = None,
    overwrite: bool = False,
    quiet: bool = False,
) -> dict[str, Any]:
    run_start = _time.time()
    output_dir = (output_dir or Path("runs/posttree/modelcompare/pb")).resolve()

    params: dict[str, Any] = {
        "sitelogl_dir": sitelogl_dir,
        "sitelogl": list(sitelogl) if sitelogl else None,
        "model_names": model_names,
        "output_dir": str(output_dir),
        "overwrite": overwrite,
        "quiet": quiet,
    }

    def error_result(message: str, category: str) -> dict[str, Any]:
        return {
            "status": "error", "command": "", "wall_time": 0.0,
            "tool_versions": {}, "params": params, "key_results": {},
            "error": message, "error_category": category,
            "data": {"cmd": [], "tool_stderr": "", "output_files": {}},
        }

    if bool(sitelogl_dir) == bool(sitelogl):
        return error_result(
            "Exactly one of --sitelogl-dir or --sitelogl must be provided.",
            "input",
        )

    try:
        if sitelogl_dir:
            groups = _groups_from_dirs(sitelogl_dir)
        else:
            groups = _groups_from_files(list(sitelogl))
    except ValueError as exc:
        return error_result(str(exc), "input")

    if model_names:
        names = [n.strip() for n in model_names.split(",") if n.strip()]
        if len(names) != len(groups):
            return error_result(
                f"--model-names provides {len(names)} label(s) but {len(groups)} "
                "model group(s) were given.",
                "input",
            )
        bad = [n for n in names if not _is_safe_name(n)]
        if bad:
            return error_result(
                f"--model-names contains unsafe label(s): {', '.join(bad)}. "
                "Labels must be single path components (no '/', '..', etc.) as "
                "they are used as output subdirectory names.",
                "input",
            )
        if len(set(names)) != len(names):
            return error_result("--model-names labels must be unique.", "input")
        labels = names
    else:
        labels = [f"model_{i + 1}" for i in range(len(groups))]

    try:
        parsed = [[_parse_sitelogl(p) for p in group] for group in groups]
    except (OSError, ValueError) as exc:
        return error_result(str(exc), "input")

    site_counts: list[int] = []
    for gi, group in enumerate(parsed):
        m = len(group[0])
        ref_sites = [row[0] for row in group[0]]
        for run in group[1:]:
            if len(run) != m:
                return error_result(
                    f"Model group '{labels[gi]}' has inconsistent site counts "
                    "between its .sitelogl files.",
                    "input",
                )
            if [row[0] for row in run] != ref_sites:
                return error_result(
                    f"Model group '{labels[gi]}' has .sitelogl files whose site "
                    "identifiers differ in order or content within the group.",
                    "input",
                )
        site_counts.append(m)

    if len(groups) >= 2:
        if any(m != site_counts[0] for m in site_counts[1:]):
            return error_result(
                "Cross-model site count mismatch: all model groups must have the "
                "same number of sites. Scores from different alignments are not "
                "comparable.",
                "input",
            )
        ref_sites = [row[0] for row in parsed[0][0]]
        for gi in range(1, len(parsed)):
            sites = [row[0] for row in parsed[gi][0]]
            if sites != ref_sites:
                return error_result(
                    f"Model group '{labels[gi]}' has site identifiers that differ "
                    "in order or content from the first model group. Scores from "
                    "different alignments or site orders are not comparable.",
                    "input",
                )

    n_sites = site_counts[0]

    if overwrite and output_dir.exists():
        shutil.rmtree(output_dir)
    elif output_dir.exists() and any(output_dir.iterdir()):
        return error_result(
            f"Output directory '{output_dir}' already exists and is non-empty. "
            "Use --overwrite to replace it.",
            "input",
        )
    output_dir.mkdir(parents=True, exist_ok=True)

    sitelogl_root = output_dir / "sitelogl"
    for gi, group in enumerate(groups):
        dest = sitelogl_root / labels[gi]
        dest.mkdir(parents=True, exist_ok=True)
        used: set[str] = set()
        for path in group:
            shutil.copy2(path, dest / _dedupe_basename(path.name, used))

    model_results: list[dict[str, Any]] = []
    for gi, runs in enumerate(parsed):
        stats = _compute_loocv_waic(runs)
        model_results.append({"model": labels[gi], "n_runs": len(runs), **stats})

    if len(model_results) >= 2:
        best_loocv = max(model_results, key=lambda mr: mr["loocv"]["score"])
        best_waic = max(model_results, key=lambda mr: mr["waic"]["score"])
        for mr in model_results:
            mr["delta_loocv"] = mr["loocv"]["score"] - best_loocv["loocv"]["score"]
            mr["delta_waic"] = mr["waic"]["score"] - best_waic["waic"]["score"]
        key_results: dict[str, Any] = {
            "n_models": len(model_results),
            "n_sites": n_sites,
            "best_model_loocv": best_loocv["model"],
            "best_model_waic": best_waic["model"],
            "best_loocv_score": best_loocv["loocv"]["score"],
            "best_loocv_quality": best_loocv["loocv"]["quality"],
            "best_waic_score": best_waic["waic"]["score"],
            "best_waic_quality": best_waic["waic"]["quality"],
        }
    else:
        only = model_results[0]
        key_results = {
            "n_models": 1,
            "n_sites": n_sites,
            "n_runs": only["n_runs"],
            "best_model_loocv": labels[0],
            "best_model_waic": labels[0],
            "best_loocv_score": only["loocv"]["score"],
            "best_loocv_quality": only["loocv"]["quality"],
            "best_loocv_ess": only["loocv"]["ess"],
            "best_loocv_pct_ess_lt10": only["loocv"]["pct_ess_lt10"],
            "best_loocv_frac_ess_lt10": only["loocv"]["frac_ess_lt10"],
            "best_waic_score": only["waic"]["score"],
            "best_waic_quality": only["waic"]["quality"],
        }

    cli_parts = ["phyloai", "posttree", "modelcompare", "pb"]
    if sitelogl_dir:
        cli_parts += ["--sitelogl-dir", sitelogl_dir]
    else:
        for seg in sitelogl:
            cli_parts += ["--sitelogl", seg]
    if model_names:
        cli_parts += ["--model-names", model_names]
    cli_parts += ["-o", str(output_dir)]
    if overwrite:
        cli_parts.append("--overwrite")
    if quiet:
        cli_parts.append("-q")
    command = shlex.join(cli_parts)

    model_fit_path = output_dir / "model_fit.csv"
    _write_model_fit_csv(model_fit_path, model_results)

    data_models: list[dict[str, Any]] = []
    for mr in model_results:
        entry: dict[str, Any] = {"model": mr["model"], "n_runs": mr["n_runs"]}
        for metric in ("loocv", "waic"):
            entry[metric] = {k: v for k, v in mr[metric].items()}
        if "delta_loocv" in mr:
            entry["delta_loocv"] = mr["delta_loocv"]
            entry["delta_waic"] = mr["delta_waic"]
        data_models.append(entry)

    result = {
        "status": "success",
        "command": command,
        "wall_time": _time.time() - run_start,
        "tool_versions": {},
        "params": params,
        "key_results": key_results,
        "error": None,
        "error_category": None,
        "data": {
            "output_files": {
                "model_fit_csv": {
                    "path": str(model_fit_path),
                    "description": "LOO-CV and wAIC scores per model with Δ values",
                }
            },
            "models": data_models,
        },
    }
    _write_result_json(output_dir, result)
    return result
