"""Model comparison via IQ-TREE ModelFinder (BIC/AIC/AICc)."""
from __future__ import annotations

import csv
import json
import shlex
import shutil
import subprocess
import time as _time
from pathlib import Path
from typing import Any

from phyloai.core.formats import FormatConverter
from phyloai.core.iqtree import (
    _detect_iqtree_version,
    _resolve_iqtree_path,
    IQTREE_COMPATIBLE_EXTENSIONS,
)
from phyloai.core.sequence_normalization import detect_seq_type

_BLOCKED_FLAGS = frozenset({"-s", "--prefix"})

AA_STANDARD_MODELS = frozenset({
    "LG", "Poisson", "cpREV", "mtREV", "Dayhoff", "mtMAM", "JTT", "WAG",
    "mtART", "mtZOA", "VT", "rtREV", "DCMut", "PMB", "HIVb", "HIVw",
    "JTTDCMut", "FLU", "Blosum62", "GTR20", "mtMet", "mtVer", "mtInv",
    "FLAVI", "Q.LG", "Q.pfam", "Q.pfam_gb", "Q.bird", "Q.mammal",
    "Q.insect", "Q.plant", "Q.yeast",
})

NT_STANDARD_MODELS = frozenset({
    "GTR", "HKY", "JC", "F81", "K2P", "K3P", "K81uf", "TN", "TNef",
    "TIM", "TIMef", "TVM", "TVMef", "SYM",
})

AA_HETEROGENEOUS_MODELS = frozenset({
    "C10", "C20", "C30", "C40", "C50", "C60",
    "EX2", "EX3", "EHO", "UL2", "UL3", "EX_EHO", "LG4M", "LG4X",
})

VALID_MRATE_TOKENS = frozenset({"E", "G", "R"})
VALID_HET_MRATE_TOKENS = frozenset({"E", "G", "R"})


def _expand_heterogeneous_models(models: list[str], het_mrate: str) -> list[str]:
    """Expand each mixture model M into variants selected by het-mrate tokens.

    Each token selects one variant family, mirroring how IQ-TREE ``-mrate``
    controls homogeneous model variants:
      E → M, M+F                  (empirical state frequencies, no rate categories)
      G → M+G4, M+F+G4
      R → M+R4, M+F+R4
    Only the requested families are produced. Default ``E,G`` yields all four.
    """
    rates = [r.strip().upper() for r in het_mrate.split(",") if r.strip()]
    expanded: list[str] = []
    for model in models:
        m = model.strip()
        for rate in rates:
            if rate == "E":
                expanded.append(m)
                expanded.append(f"{m}+F")
            elif rate == "G":
                expanded.append(f"{m}+G4")
                expanded.append(f"{m}+F+G4")
            elif rate == "R":
                expanded.append(f"{m}+R4")
                expanded.append(f"{m}+F+R4")
    return expanded


def _detect_seq_type(matrix: Path) -> str:
    converter = FormatConverter()
    alignment = converter.read(matrix)
    return detect_seq_type([str(record.seq) for record in alignment])


def _validate_inputs(
    *,
    matrix: Path,
    homogeneous_models: list[str],
    mrate: str,
    heterogeneous_models: list[str] | None,
    het_mrate: str,
    seq_type: str,
    threads: str,
    prefix: str = "modelcompare",
    tool_args: str | None = None,
    overwrite: bool = False,
    resume: bool = False,
) -> list[str]:
    """Validate inputs; `seq_type` must already be resolved to AA or NT."""
    errors: list[str] = []

    if not matrix.exists():
        errors.append(f"--matrix does not exist: {matrix}")
        return errors
    if not matrix.is_file():
        errors.append(f"--matrix is not a regular file: {matrix}")
        return errors
    ext = matrix.suffix.lower()
    if ext not in IQTREE_COMPATIBLE_EXTENSIONS:
        errors.append(
            f"Unsupported matrix extension: {ext}. "
            f"Supported: {', '.join(sorted(IQTREE_COMPATIBLE_EXTENSIONS))}"
        )
        return errors

    if prefix and (Path(prefix).name != prefix or prefix in (".", "..")):
        errors.append(
            f"--prefix must be a single filename (no path separators, '..', etc.); "
            f"got {prefix!r}."
        )

    if not homogeneous_models:
        errors.append(
            "--homogeneous-model is required (comma-separated list of standard models)."
        )

    mrate_tokens = [t.strip().upper() for t in mrate.split(",") if t.strip()]
    if not mrate_tokens:
        errors.append("--mrate must contain at least one of E, G, R.")
    else:
        bad = [t for t in mrate_tokens if t not in VALID_MRATE_TOKENS]
        if bad:
            errors.append(
                f"--mrate contains invalid token(s): {', '.join(bad)}. Valid: E, G, R."
            )

    if heterogeneous_models:
        het_mrate_tokens = [t.strip().upper() for t in het_mrate.split(",") if t.strip()]
        if not het_mrate_tokens:
            errors.append("--het-mrate must contain at least one of E, G, R.")
        else:
            bad = [t for t in het_mrate_tokens if t not in VALID_HET_MRATE_TOKENS]
            if bad:
                errors.append(
                    f"--het-mrate contains invalid token(s): {', '.join(bad)}. "
                    f"Valid: E (no rate variants), G, R, or G,R."
                )

    st = seq_type.upper()
    if st not in ("AA", "NT"):
        errors.append(f"Unresolved seq-type '{seq_type}'; expected AA or NT.")
        return errors

    if threads != "auto":
        try:
            if int(threads) < 1:
                raise ValueError
        except ValueError:
            errors.append(f"--threads must be a positive integer or 'auto', got {threads!r}.")

    if homogeneous_models:
        model_set = AA_STANDARD_MODELS if st == "AA" else NT_STANDARD_MODELS
        bad = [m for m in homogeneous_models if m not in model_set]
        if bad:
            errors.append(
                f"--homogeneous-model contains model(s) not valid for {st} data: {', '.join(bad)}."
            )

    if heterogeneous_models:
        if st != "AA":
            errors.append(
                f"--heterogeneous-model requires amino-acid (AA) data; got seq-type {st}. "
                "Heterogeneous mixture models are AA-only."
            )
        else:
            bad = [m for m in heterogeneous_models if m not in AA_HETEROGENEOUS_MODELS]
            if bad:
                errors.append(
                    f"--heterogeneous-model contains unknown mixture model(s): {', '.join(bad)}."
                )

    if overwrite and resume:
        errors.append("--overwrite and --resume are mutually exclusive.")

    if tool_args:
        tokens = shlex.split(tool_args)
        for tok in tokens:
            flag = tok.split("=", 1)[0]
            if flag in _BLOCKED_FLAGS:
                errors.append(
                    f"--tool-args contains blocked flag '{tok}' (managed by PhyloAI)."
                )

    return errors


def _tool_arg_flags(tool_args: str | None) -> set[str]:
    """Flag names present in --tool-args (handles both '--flag value' and '--flag=value')."""
    if not tool_args:
        return set()
    return {tok.split("=", 1)[0] for tok in shlex.split(tool_args)}


def _build_cmd(
    *,
    executable: str,
    matrix: Path,
    homogeneous_models: list[str],
    mrate: str,
    expanded_het: list[str] | None,
    prefix: str,
    threads: str,
    tool_args: str | None,
) -> list[str]:
    overridden = _tool_arg_flags(tool_args)
    cmd = [executable, "-s", str(matrix)]
    if "-m" not in overridden:
        cmd += ["-m", "MF"]
    if "-mset" not in overridden:
        cmd += ["-mset", ",".join(homogeneous_models)]
    if "-mrate" not in overridden:
        cmd += ["-mrate", mrate]
    if "-cmin" not in overridden:
        cmd += ["-cmin", "4"]
    if "-cmax" not in overridden:
        cmd += ["-cmax", "4"]
    if expanded_het and "-madd" not in overridden:
        cmd += ["-madd", ",".join(expanded_het)]
    if "--prefix" not in overridden:
        cmd += ["--prefix", prefix]
    if "-T" not in overridden:
        cmd += ["-T", threads]
    if tool_args:
        cmd += shlex.split(tool_args)
    return cmd


def _pop_weight(tokens: list[str]) -> tuple[str, float]:
    """Pop one [+/-]<value> weight pair from the end of a token list.

    IQ-TREE writes the 95%-confidence-set marker (+/-) as a separate token
    before the weight value; also tolerate the sign being attached.
    """
    value = tokens.pop()
    sign = ""
    if tokens and tokens[-1] in ("+", "-"):
        sign = tokens.pop()
    elif value[:1] in ("+", "-"):
        sign = value[0]
        value = value[1:]
    return sign, float(value)


def _parse_model_row(line: str) -> dict[str, Any] | None:
    tokens = line.split()
    if len(tokens) < 8:
        return None
    try:
        w_bic_sign, w_bic = _pop_weight(tokens)
        bic = float(tokens.pop())
        w_aicc_sign, w_aicc = _pop_weight(tokens)
        aicc = float(tokens.pop())
        w_aic_sign, w_aic = _pop_weight(tokens)
        aic = float(tokens.pop())
        logl = float(tokens.pop())
        model = " ".join(tokens).strip()
    except ValueError:
        return None
    if not model:
        return None
    return {
        "model": model,
        "logl": logl,
        "aic": aic,
        "w_aic": w_aic,
        "in_aic_95": w_aic_sign == "+",
        "aicc": aicc,
        "w_aicc": w_aicc,
        "in_aicc_95": w_aicc_sign == "+",
        "bic": bic,
        "w_bic": w_bic,
        "in_bic_95": w_bic_sign == "+",
    }


def _parse_modelfinder_results(iqtree_file: Path) -> list[dict[str, Any]]:
    """Parse the 'List of models sorted by BIC scores:' section of an .iqtree file."""
    lines = iqtree_file.read_text(errors="ignore").splitlines()
    start = -1
    for i, line in enumerate(lines):
        if line.strip().startswith("List of models sorted by BIC scores:"):
            start = i + 2  # skip blank line + column header
            break
    if start < 0:
        raise ValueError(
            f"No 'List of models sorted by BIC scores:' section found in {iqtree_file}"
        )
    models: list[dict[str, Any]] = []
    for line in lines[start:]:
        if not line.strip():
            break
        row = _parse_model_row(line)
        if row is not None:
            models.append(row)
    if not models:
        raise ValueError(f"No model rows parsed from {iqtree_file}")
    return models


def _write_model_fit_csv(path: Path, models: list[dict[str, Any]]) -> None:
    fieldnames = [
        "Rank", "Model", "LogL", "AIC", "w_AIC", "In_AIC_95",
        "AICc", "w_AICc", "In_AICc_95", "BIC", "w_BIC", "In_BIC_95",
    ]
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        for i, m in enumerate(models, 1):
            writer.writerow({
                "Rank": i,
                "Model": m["model"],
                "LogL": m["logl"],
                "AIC": m["aic"],
                "w_AIC": m["w_aic"],
                "In_AIC_95": m["in_aic_95"],
                "AICc": m["aicc"],
                "w_AICc": m["w_aicc"],
                "In_AICc_95": m["in_aicc_95"],
                "BIC": m["bic"],
                "w_BIC": m["w_bic"],
                "In_BIC_95": m["in_bic_95"],
            })


def _write_result_json(output_dir: Path, result: dict[str, Any]) -> None:
    (output_dir / "result.json").write_text(json.dumps(result, indent=2))


def run_modelcompare_iqtree(
    *,
    matrix: Path,
    homogeneous_model: str,
    mrate: str = "E,G",
    heterogeneous_model: str | None = None,
    het_mrate: str = "E,G",
    seq_type: str = "auto",
    prefix: str = "modelcompare",
    output_dir: Path | None = None,
    threads: str = "auto",
    iqtree_path: str | None = None,
    tool_args: str | None = None,
    overwrite: bool = False,
    resume: bool = False,
    dry_run: bool = False,
    quiet: bool = False,
) -> dict[str, Any]:
    run_start = _time.time()
    output_dir = (output_dir or Path("runs/posttree/modelcompare/iqtree")).resolve()
    matrix = matrix.resolve()

    homogeneous_models = [m.strip() for m in homogeneous_model.split(",") if m.strip()]
    heterogeneous_models = (
        [m.strip() for m in heterogeneous_model.split(",") if m.strip()]
        if heterogeneous_model else None
    )
    mrate_norm = ",".join(t.strip().upper() for t in mrate.split(",") if t.strip())
    het_mrate_norm = ",".join(t.strip().upper() for t in het_mrate.split(",") if t.strip())

    params: dict[str, Any] = {
        "matrix": str(matrix),
        "homogeneous_model": ",".join(homogeneous_models),
        "mrate": mrate_norm,
        "heterogeneous_model": ",".join(heterogeneous_models) if heterogeneous_models else None,
        "het_mrate": het_mrate_norm if heterogeneous_models else None,
        "seq_type": seq_type,
        "detected_seq_type": None,
        "prefix": prefix,
        "output_dir": str(output_dir),
        "threads": threads,
        "iqtree_path": iqtree_path,
        "tool_args": tool_args,
        "overwrite": overwrite,
        "resume": resume,
        "dry_run": dry_run,
        "quiet": quiet,
    }

    def error_result(message: str, category: str) -> dict[str, Any]:
        return {
            "status": "error", "command": "", "wall_time": 0.0,
            "tool_versions": {}, "params": params, "key_results": {},
            "error": message, "error_category": category,
            "data": {"cmd": [], "tool_stderr": "", "tool_log": None, "output_files": {}},
        }

    detected_seq_type: str
    if seq_type == "auto":
        try:
            detected_seq_type = _detect_seq_type(matrix)
        except Exception as exc:
            return error_result(f"Cannot read --matrix as FASTA/PHYLIP/NEXUS: {exc}", "input")
    else:
        detected_seq_type = seq_type.upper()
        try:
            actual = _detect_seq_type(matrix)
        except Exception as exc:
            return error_result(f"Cannot read --matrix as FASTA/PHYLIP/NEXUS: {exc}", "input")
        if actual != detected_seq_type:
            return error_result(
                f"--seq-type {detected_seq_type} conflicts with the matrix's "
                f"detected molecule type ({actual}).",
                "input",
            )
    params["detected_seq_type"] = detected_seq_type

    errors = _validate_inputs(
        matrix=matrix,
        homogeneous_models=homogeneous_models,
        mrate=mrate_norm,
        heterogeneous_models=heterogeneous_models,
        het_mrate=het_mrate_norm,
        seq_type=detected_seq_type,
        threads=threads,
        prefix=prefix,
        tool_args=tool_args,
        overwrite=overwrite,
        resume=resume,
    )
    if errors:
        return error_result("; ".join(errors), "input")

    expanded_het = (
        _expand_heterogeneous_models(heterogeneous_models, het_mrate_norm)
        if heterogeneous_models else None
    )

    if not dry_run:
        if overwrite and output_dir.exists():
            shutil.rmtree(output_dir)
        elif not resume and output_dir.exists() and any(output_dir.iterdir()):
            return error_result(
                f"Output directory '{output_dir}' already exists and is non-empty. "
                "Use --overwrite to replace it.",
                "input",
            )
        output_dir.mkdir(parents=True, exist_ok=True)

    try:
        executable = _resolve_iqtree_path(iqtree_path, dry_run)
    except (ValueError, FileNotFoundError) as exc:
        return error_result(str(exc), "env")
    tool_versions = {"iqtree3": "dry-run"} if dry_run else _detect_iqtree_version(executable)

    cmd = _build_cmd(
        executable=executable,
        matrix=matrix,
        homogeneous_models=homogeneous_models,
        mrate=mrate_norm,
        expanded_het=expanded_het,
        prefix=prefix,
        threads=threads,
        tool_args=tool_args,
    )

    cli_parts = [
        "phyloai", "posttree", "modelcompare", "iqtree",
        "--matrix", str(matrix),
        "--homogeneous-model", ",".join(homogeneous_models),
        "--mrate", mrate_norm,
        "--seq-type", seq_type,
        "--prefix", prefix,
        "--threads", str(threads),
        "-o", str(output_dir),
    ]
    if heterogeneous_models:
        cli_parts += ["--heterogeneous-model", ",".join(heterogeneous_models), "--het-mrate", het_mrate_norm]
    if iqtree_path:
        cli_parts += ["--iqtree-path", iqtree_path]
    if tool_args:
        cli_parts += ["--tool-args", tool_args]
    if overwrite:
        cli_parts.append("--overwrite")
    if resume:
        cli_parts.append("--resume")
    if dry_run:
        cli_parts.append("--dry-run")
    if quiet:
        cli_parts.append("-q")
    command = shlex.join(cli_parts)

    if dry_run:
        return {
            "status": "success", "command": command, "wall_time": 0.0,
            "tool_versions": tool_versions, "params": params, "key_results": {},
            "error": None,
            "data": {"cmd": cmd, "tool_stderr": "", "tool_log": None, "output_files": {}},
        }

    iqtree_dir = output_dir / "iqtree"
    iqtree_dir.mkdir(parents=True, exist_ok=True)
    proc = subprocess.run(
        cmd, stdout=None, stderr=subprocess.PIPE, text=True, cwd=str(iqtree_dir),
    )
    tool_stderr = proc.stderr.strip() if proc.stderr else ""
    iqtree_report = iqtree_dir / f"{prefix}.iqtree"
    output_files: dict[str, Any] = {
        "model_fit_csv": {"path": str(output_dir / "model_fit.csv"), "description": "Model comparison table sorted by BIC"},
        "iqtree_report": {"path": str(iqtree_report), "description": "IQ-TREE native report"},
        "iqtree_log": {"path": str(iqtree_dir / f"{prefix}.log"), "description": "IQ-TREE console log"},
    }

    if proc.returncode:
        result = error_result(f"IQ-TREE exited with code {proc.returncode}", "tool")
        result.update({
            "command": command, "wall_time": _time.time() - run_start,
            "tool_versions": tool_versions,
        })
        result["data"] = {
            "cmd": cmd, "tool_stderr": tool_stderr,
            "tool_log": str(iqtree_dir / f"{prefix}.log"), "output_files": output_files,
        }
        _write_result_json(output_dir, result)
        return result

    try:
        models = _parse_modelfinder_results(iqtree_report)
    except (OSError, ValueError) as exc:
        result = error_result(str(exc), "output")
        result.update({
            "command": command, "wall_time": _time.time() - run_start,
            "tool_versions": tool_versions,
        })
        result["data"] = {
            "cmd": cmd, "tool_stderr": tool_stderr,
            "tool_log": str(iqtree_dir / f"{prefix}.log"), "output_files": output_files,
        }
        _write_result_json(output_dir, result)
        return result

    best_bic = min(models, key=lambda m: m["bic"])
    best_aic = min(models, key=lambda m: m["aic"])
    best_aicc = min(models, key=lambda m: m["aicc"])

    model_fit_path = output_dir / "model_fit.csv"
    _write_model_fit_csv(model_fit_path, models)

    key_results: dict[str, Any] = {
        "best_model_bic": best_bic["model"],
        "best_model_aic": best_aic["model"],
        "best_model_aicc": best_aicc["model"],
        "best_model_bic_value": best_bic["bic"],
        "best_model_bic_weight": best_bic["w_bic"],
        "n_models_tested": len(models),
        "madd_expanded": ",".join(expanded_het) if expanded_het else None,
        "n_madd_expanded": len(expanded_het) if expanded_het else 0,
    }

    result = {
        "status": "success",
        "command": command,
        "wall_time": _time.time() - run_start,
        "tool_versions": tool_versions,
        "params": params,
        "key_results": key_results,
        "error": None,
        "data": {
            "cmd": cmd,
            "tool_stderr": tool_stderr,
            "tool_log": str(iqtree_dir / f"{prefix}.log"),
            "output_files": output_files,
            "models": [{"rank": i, **m} for i, m in enumerate(models, 1)],
        },
    }
    _write_result_json(output_dir, result)
    return result
