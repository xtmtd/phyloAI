"""AliSim alignment simulation via IQ-TREE3.

`phyloai posttree simulate alisim iqtree` runs IQ-TREE3 ``--alisim`` in
single-parameter mode (one tree/model/length) or batch mode (rows sampled
from a ``params.tsv`` table under one of three strategies: ``complete``,
``mixed``, or ``pdf``).  Batch mode records full per-simulation provenance
in ``params_sampled.tsv`` and supports checkpoint/resume.
"""
from __future__ import annotations

import csv
import json
import math
import random
import shlex
import shutil
import time as _time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Any, Mapping

from Bio import SeqIO

from phyloai.core.checkpoint import (
    Checkpoint,
    CheckpointTask,
    load_checkpoint,
    save_checkpoint_atomic,
    validate_resume_params,
)
from phyloai.core.iqtree import _detect_iqtree_version, _resolve_iqtree_path
from phyloai.core.runner import Runner
from phyloai.core.sequence_output_validation import validate_fasta_output

REQUIRED_COLUMNS = (
    "id", "seqtype", "length", "subs_model", "subs_rate", "freq",
    "prop_inv", "rate_heterogeneity", "rate_categories", "rate_param",
    "tree_path",
)

SAMPLED_COLUMNS = (
    "simulation_id", "source_id", "seqtype", "length", "subs_model",
    "subs_rate", "freq", "prop_inv", "rate_heterogeneity",
    "rate_categories", "rate_param", "tree_path", "seed",
)

_BLOCKED_FLAGS = frozenset({
    "--alisim", "-t", "-m", "-p", "-q", "-Q", "--seqtype", "--length",
    "--out-format", "-af", "--num-alignments", "-T", "--seed", "--prefix",
})

VALID_OVERRIDE_KEYS = frozenset({"length", "prop_inv"})
VALID_PDF_PARAMS = frozenset({"length", "prop_inv", "rate_param"})
OUT_FORMAT_EXT = {"fasta": ".fa", "phy": ".phy"}


# ===================================================================
# Model string reconstruction
# ===================================================================

def build_model_string(row: Mapping[str, str]) -> str:
    """Reconstruct an IQ-TREE ``-m`` string from one params-table row."""
    model = row["subs_model"]
    if row["subs_rate"]:
        model += "{" + row["subs_rate"].replace("/", ",") + "}"
    if row["freq"]:
        model += "+F{" + row["freq"].replace("/", ",") + "}"
    if row["prop_inv"]:
        model += "+I{" + row["prop_inv"] + "}"
    if row["rate_heterogeneity"]:
        model += (
            f'+{row["rate_heterogeneity"]}{row["rate_categories"]}'
            + "{" + row["rate_param"].replace("/", ",") + "}"
        )
    return model


# ===================================================================
# Table loading and validation
# ===================================================================

def load_params_table(path: Path) -> list[dict[str, str]]:
    """Read a ``params.tsv`` table, validating columns and row content."""
    if not path.exists():
        raise ValueError(f"--model-params does not exist: {path}")
    if not path.is_file():
        raise ValueError(f"--model-params is not a regular file: {path}")
    if path.stat().st_size == 0:
        raise ValueError(f"--model-params is empty: {path}")

    rows: list[dict[str, str]] = []
    with open(path, encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"--model-params is not a valid TSV: {path}")
        missing = [col for col in REQUIRED_COLUMNS if col not in reader.fieldnames]
        if missing:
            raise ValueError(
                f"--model-params is missing required columns: {', '.join(missing)}"
            )
        for row in reader:
            rows.append({key: (row.get(key) or "").strip() for key in REQUIRED_COLUMNS})

    if not rows:
        raise ValueError(f"--model-params contains no data rows: {path}")

    for idx, row in enumerate(rows):
        for col in ("id", "seqtype", "length", "subs_model"):
            if not row[col]:
                raise ValueError(
                    f"--model-params row {idx + 1}: required column {col!r} is empty"
                )
        if row["seqtype"] not in {"AA", "DNA"}:
            raise ValueError(
                f"--model-params row {idx + 1}: seqtype must be AA or DNA, got {row['seqtype']!r}"
            )
        try:
            if int(row["length"]) < 1:
                raise ValueError
        except ValueError:
            raise ValueError(
                f"--model-params row {idx + 1}: length must be a positive integer, "
                f"got {row['length']!r}"
            )
    return rows


# ===================================================================
# Sampling
# ===================================================================

def _fd_bins(values: list[float]) -> tuple[list[float], list[int]]:
    """Freedman-Diaconis bins with safe fallback for constant/singleton input."""
    if not values:
        return [], []
    if len(values) == 1:
        v = values[0]
        return [v - 0.5, v + 0.5], [1]
    lo, hi = min(values), max(values)
    if hi == lo:
        return [lo - 0.5, lo + 0.5], [len(values)]
    sorted_values = sorted(values)
    q25 = sorted_values[(len(values) - 1) // 4]
    q75 = sorted_values[(3 * (len(values) - 1)) // 4]
    iqr = q75 - q25
    width = (2.0 * iqr / (len(values) ** (1.0 / 3.0))) if iqr > 0 else (hi - lo) / math.sqrt(len(values))
    if not math.isfinite(width) or width <= 0:
        width = (hi - lo) / math.sqrt(len(values))
    n_bins = max(1, math.ceil((hi - lo) / width))
    edges = [lo + i * width for i in range(n_bins + 1)]
    counts = [0] * n_bins
    for value in values:
        index = min(int((value - lo) / width), n_bins - 1)
        counts[index] += 1
    return edges, counts


def _sample_density(
    values: list[float], rng: random.Random, noise_scale: float,
) -> float:
    """Draw one value from the empirical distribution via FD-binned resampling.

    The draw is clamped to the empirical value span so bin midpoints never
    fall outside the observed range of the parameter.
    """
    lo, hi = min(values), max(values)
    edges, counts = _fd_bins(values)
    total = sum(counts)
    if total <= 0:
        return rng.choice(values)
    index = rng.choices(range(len(counts)), weights=counts, k=1)[0]
    left, right = edges[index], edges[index + 1]
    if noise_scale <= 0.0 or right <= left:
        value = (left + right) / 2.0
    else:
        value = left + rng.uniform(0.0, right - left) * noise_scale
    return max(lo, min(value, hi))


def _sample_rate_group(
    rows: list[dict[str, str]],
    rng: random.Random,
) -> tuple[str, str, str]:
    """Sample (rate_heterogeneity, rate_categories, rate_param) as one unit.

    The three fields always come from the same source row, keeping Gamma
    alpha vs FreeRate pairs intact.  Distinct (type, categories, param)
    combinations are the sampling units so that model configurations are
    not diluted by row repetition in the table.
    """
    by_key: dict[tuple[str, str, str], dict[str, str]] = {}
    for row in rows:
        key = (row["rate_heterogeneity"], row["rate_categories"], row["rate_param"])
        if key not in by_key:
            by_key[key] = row
    if not by_key:
        return "", "", ""
    source = rng.choice(list(by_key.values()))
    return source["rate_heterogeneity"], source["rate_categories"], source["rate_param"]


def _sample_prop_inv(
    rows: list[dict[str, str]],
    rng: random.Random,
    *,
    density: bool,
    noise_scale: float,
    empirical_values: list[float],
) -> str:
    """Two-step +I sampling: presence decision, then value.

    The presence/absence ratio is the empirical fraction of rows with
    non-empty ``prop_inv`` values (all rows, not just those with trees,
    since ``prop_inv`` is sampled independently of ``tree_path``).
    """
    present_ratio = sum(1 for row in rows if row["prop_inv"]) / len(rows)
    if present_ratio <= 0.0:
        return ""
    if present_ratio >= 1.0:
        present = True
    else:
        present = rng.random() < present_ratio
    if not present:
        return ""
    if density and empirical_values:
        return f"{_sample_density(empirical_values, rng, noise_scale):.6f}"
    non_empty = [row["prop_inv"] for row in rows if row["prop_inv"]]
    return rng.choice(non_empty)


def sample_batch_rows(
    rows: list[dict[str, str]],
    *,
    strategy: str,
    n: int,
    rng: random.Random,
    pdf_params: tuple[str, ...],
    noise_scale: float,
    overrides: dict[str, str],
) -> list[dict[str, str]]:
    """Sample ``n`` complete rows from a validated params table.

    ``complete``: rows kept intact (tree_path must be non-empty).
    ``mixed``: model core group, rate group, and independent parameters are
    sampled separately; ``prop_inv`` and the rate group preserve empirical
    presence/type ratios.
    ``pdf``: built on ``mixed`` with histogram density resampling for the
    parameters in ``pdf_params``.
    ``overrides`` fix specific keys and are never sampled.
    """
    if strategy == "complete":
        usable = [row for row in rows if row["tree_path"]]
        if not usable:
            raise ValueError(
                "no rows with a non-empty tree_path available for 'complete' sampling"
            )
        sampled = [dict(rng.choice(usable)) for _ in range(n)]
    else:
        core_groups: dict[tuple[str, str, str, str], dict[str, str]] = {}
        for row in rows:
            key = (row["seqtype"], row["subs_model"], row["subs_rate"], row["freq"])
            if key not in core_groups:
                core_groups[key] = row
        tree_rows = [row for row in rows if row["tree_path"]]
        empirical_lengths = [float(row["length"]) for row in rows if row["length"]]
        empirical_prop_inv = [float(row["prop_inv"]) for row in rows if row["prop_inv"]]
        gamma_alphas = [
            float(row["rate_param"]) for row in rows
            if row["rate_heterogeneity"] == "G" and row["rate_param"]
        ]

        sampled = []
        for _ in range(n):
            core = rng.choice(list(core_groups.values()))
            rate_het, rate_cats, rate_param = _sample_rate_group(rows, rng)
            row: dict[str, str] = {
                "id": core["id"],
                "seqtype": core["seqtype"],
                "subs_model": core["subs_model"],
                "subs_rate": core["subs_rate"],
                "freq": core["freq"],
                "rate_heterogeneity": rate_het,
                "rate_categories": rate_cats,
                "rate_param": rate_param,
                "tree_path": rng.choice(tree_rows)["tree_path"],
                "length": rng.choice(rows)["length"],
                "prop_inv": _sample_prop_inv(
                    rows, rng, density=("prop_inv" in pdf_params),
                    noise_scale=noise_scale,
                    empirical_values=empirical_prop_inv,
                ),
            }

            if strategy == "pdf":
                if "length" in pdf_params and empirical_lengths:
                    value = _sample_density(empirical_lengths, rng, noise_scale)
                    row["length"] = str(max(1, int(round(value))))
                if "prop_inv" in pdf_params and empirical_prop_inv and row["prop_inv"]:
                    value = _sample_density(empirical_prop_inv, rng, noise_scale)
                    row["prop_inv"] = f"{max(0.0, min(value, 1.0)):.6f}"
                if "rate_param" in pdf_params and gamma_alphas and rate_het == "G":
                    value = _sample_density(gamma_alphas, rng, noise_scale)
                    row["rate_param"] = f"{value:.6f}"

            for key, value in overrides.items():
                row[key] = value
            sampled.append(row)

    return sampled


# ===================================================================
# Command construction
# ===================================================================

def _check_managed_flag_conflict(tool_args: str) -> None:
    tokens = shlex.split(tool_args)
    for token in tokens:
        if token in _BLOCKED_FLAGS:
            raise ValueError(f"Blocked managed flag in --tool-args: {token}")


def _build_alisim_cmd(
    *,
    executable: str,
    msa_prefix: str,
    seq_type: str,
    ref_tree: str,
    model: str | None,
    model_partitions: str | None,
    length: int | None,
    out_format: str,
    num_alignments: int,
    iqtree_threads: int,
    seed: int,
    tool_args: str | None,
) -> list[str]:
    cmd = [
        executable, "--alisim", msa_prefix,
        "--seqtype", seq_type,
        "-t", ref_tree,
    ]
    if model is not None:
        cmd.extend(["-m", model])
    if model_partitions is not None:
        cmd.extend(["-p", model_partitions])
    if length is not None:
        cmd.extend(["--length", str(length)])
    cmd.extend(["--out-format", out_format])
    cmd.extend(["--num-alignments", str(num_alignments)])
    cmd.extend(["-T", str(iqtree_threads)])
    cmd.extend(["--seed", str(seed)])
    if tool_args:
        cmd.extend(shlex.split(tool_args))
    return cmd


# ===================================================================
# Output validation
# ===================================================================

def _validate_output(path: Path, out_format: str) -> tuple[bool, list[str]]:
    if out_format == "fasta":
        result = validate_fasta_output(path, require_aligned=True)
        return result.ok, list(result.warnings)
    if not path.exists() or path.stat().st_size == 0:
        return False, ["generated PHYLIP output is empty"]
    try:
        records = list(SeqIO.parse(str(path), "phylip"))
    except Exception as exc:
        return False, [f"could not parse generated PHYLIP output: {exc}"]
    if not records:
        return False, ["generated PHYLIP output contains no records"]
    lengths = {len(rec.seq) for rec in records}
    if any(length == 0 for length in lengths):
        return False, ["generated PHYLIP output contains empty sequences"]
    if len(lengths) > 1:
        return False, [f"generated PHYLIP MSA has unequal sequence lengths: {sorted(lengths)}"]
    return True, []


def _wrap_fasta(path: Path) -> None:
    records = list(SeqIO.parse(str(path), "fasta"))
    with open(path, "w", encoding="utf-8") as fh:
        for record in records:
            fh.write(f">{record.id}\n")
            seq = str(record.seq)
            for i in range(0, len(seq), 60):
                fh.write(seq[i:i + 60] + "\n")


# ===================================================================
# Single-parameter mode execution
# ===================================================================

def _run_single_mode(
    *,
    ref_tree: Path,
    model: str | None,
    model_partitions: Path | None,
    seq_type: str,
    length: int | None,
    msa_prefix: str,
    num_alignments: int,
    out_format: str,
    iqtree_threads: int,
    seed: int,
    iqtree_exe: str,
    tool_args: str | None,
    output_dir: Path,
    dry_run: bool,
) -> dict[str, Any]:
    ext = OUT_FORMAT_EXT[out_format]
    if model_partitions is not None:
        resolved_partitions = str(model_partitions.resolve())
        model_arg: str | None = None
    else:
        resolved_partitions = None
        model_arg = model

    cmd = _build_alisim_cmd(
        executable=iqtree_exe,
        msa_prefix=msa_prefix,
        seq_type=seq_type,
        ref_tree=str(ref_tree.resolve()),
        model=model_arg,
        model_partitions=resolved_partitions,
        length=length,
        out_format=out_format,
        num_alignments=num_alignments,
        iqtree_threads=iqtree_threads,
        seed=seed,
        tool_args=tool_args,
    )

    tool_stderr = ""
    output_files: dict[str, Any] = {}
    n_generated = 0

    if dry_run:
        data: dict[str, Any] = {
            "cmd": cmd,
            "tool_stderr": "",
            "output_files": {},
        }
        return data, n_generated

    work_dir = output_dir / "_work"
    work_dir.mkdir(parents=True, exist_ok=True)

    result = Runner().run(cmd, "iqtree3", cwd=work_dir)
    tool_stderr = result.stderr

    msa_dir = output_dir / "MSAs"
    logs_dir = output_dir / "logs"
    msa_dir.mkdir(exist_ok=True)
    logs_dir.mkdir(exist_ok=True)

    if num_alignments == 1:
        names = [f"{msa_prefix}{ext}"]
    else:
        names = [f"{msa_prefix}_{i}{ext}" for i in range(1, num_alignments + 1)]

    for name in names:
        source = work_dir / name
        if source.exists():
            target = msa_dir / name
            shutil.move(str(source), str(target))
            if out_format == "fasta":
                _wrap_fasta(target)
            ok, warnings = _validate_output(target, out_format)
            if ok:
                n_generated += 1
            else:
                tool_stderr = (tool_stderr + "\n" + "; ".join(warnings)).strip()

    for suffix in (".iqtree", ".log"):
        source = work_dir / f"{msa_prefix}{suffix}"
        if source.exists():
            shutil.move(str(source), str(logs_dir / f"{msa_prefix}{suffix}"))

    if result.returncode != 0:
        raise RuntimeError(
            f"IQ-TREE AliSim failed with exit code {result.returncode}.\n"
            f"Command: {' '.join(cmd)}\n{result.stderr}"
        )

    output_files = {
        "msas": {"path": str(msa_dir), "description": "Simulated MSA files"},
        "iqtree_report": {
            "path": str(logs_dir / f"{msa_prefix}.iqtree"),
            "description": "IQ-TREE AliSim report",
        },
        "iqtree_log": {
            "path": str(logs_dir / f"{msa_prefix}.log"),
            "description": "IQ-TREE console log",
        },
    }
    shutil.rmtree(work_dir, ignore_errors=True)

    data = {
        "cmd": cmd,
        "tool_stderr": tool_stderr,
        "output_files": output_files,
    }
    return data, n_generated


# ===================================================================
# Batch-mode worker
# ===================================================================

def _run_simulation_worker(args: tuple[str, str, str, str, int, int, str, str | None]) -> dict[str, Any]:
    """Run one batch simulation in its own working directory."""
    (
        simulation_id, seq_type, ref_tree, model, length, seed,
        iqtree_exe, tool_args,
    ) = args
    work_dir = Path(f"./_work_{simulation_id}")
    work_dir.mkdir(parents=True, exist_ok=True)

    cmd = _build_alisim_cmd(
        executable=iqtree_exe,
        msa_prefix=simulation_id,
        seq_type=seq_type,
        ref_tree=ref_tree,
        model=model,
        model_partitions=None,
        length=length,
        out_format="fasta",
        num_alignments=1,
        iqtree_threads=1,
        seed=seed,
        tool_args=tool_args,
    )

    start = _time.time()
    try:
        result = Runner().run(cmd, "iqtree3", cwd=work_dir)
        wall_time = _time.time() - start
        status = "success" if result.returncode == 0 else "failed"
        log_text = f"{result.stdout}\n{result.stderr}".strip()
        reason = None if status == "success" else result.stderr.strip() or "IQ-TREE failed"
        return {
            "simulation_id": simulation_id,
            "status": status,
            "wall_time": round(wall_time, 3),
            "cmd": cmd,
            "log_file": f"logs/{simulation_id}.log",
            "log_text": log_text,
            "output_file": f"MSAs/{simulation_id}.fa",
            "output_path": work_dir / f"{simulation_id}.fa",
            "reason": reason,
        }
    finally:
        shutil.rmtree(work_dir, ignore_errors=True)


def _read_params_sampled(output_dir: Path) -> list[dict[str, str]]:
    path = output_dir / "params_sampled.tsv"
    if not path.exists():
        return []
    with open(path, encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        return [dict(row) for row in reader]


def _write_params_sampled(output_dir: Path, rows: list[dict[str, str]]) -> None:
    with open(output_dir / "params_sampled.tsv", "w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=SAMPLED_COLUMNS, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


# ===================================================================
# PDF density plots
# ===================================================================

def _plot_density(empirical: list[float], simulated: list[float], title: str, path: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    combined = empirical + simulated
    if not combined or max(combined) == min(combined):
        fig, ax = plt.subplots(figsize=(6, 4))
        ax.text(0.5, 0.5, "constant parameter distribution", ha="center", va="center")
        ax.axis("off")
    else:
        fig, ax = plt.subplots(figsize=(6, 4))
        bins = _fd_bins(combined)[0]
        if len(bins) < 2:
            bins = 20
        ax.hist(empirical, bins=bins, density=True, alpha=0.6, color="blue", label="empirical")
        ax.hist(simulated, bins=bins, density=True, alpha=0.6, color="orange", label="simulated")
        ax.set_xlabel(title)
        ax.set_ylabel("density")
        ax.legend()
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def _generate_density_plots(
    rows: list[dict[str, str]],
    sampled: list[dict[str, str]],
    pdf_params: tuple[str, ...],
    output_dir: Path,
) -> dict[str, Any]:
    plots: dict[str, Any] = {}
    param_columns: dict[str, tuple[str, str]] = {
        "length": ("alignment length", "length_density"),
        "prop_inv": ("proportion of invariable sites", "prop_inv_density"),
        "rate_param": ("Gamma alpha parameter", "rate_param_density"),
    }
    for param in pdf_params:
        label, filename = param_columns[param]
        empirical = []
        for row in rows:
            value = row[param]
            if value and param != "rate_param":
                try:
                    empirical.append(float(value))
                except ValueError:
                    pass
        if param == "rate_param":
            empirical = [float(row["rate_param"]) for row in rows
                         if row["rate_heterogeneity"] == "G" and row["rate_param"]]
        simulated = []
        for row in sampled:
            value = row[param]
            if value and param != "rate_param":
                try:
                    simulated.append(float(row[param]))
                except ValueError:
                    pass
        if param == "rate_param":
            simulated = [float(row["rate_param"]) for row in sampled
                         if row["rate_heterogeneity"] == "G" and row["rate_param"]]
        path = output_dir / "plots" / f"{filename}.pdf"
        path.parent.mkdir(parents=True, exist_ok=True)
        _plot_density(empirical, simulated, label, path)
        plots[filename] = {
            "path": str(path),
            "description": f"Empirical vs simulated density: {label}",
        }
    return plots


# ===================================================================
# Batch-mode orchestration
# ===================================================================

def _run_batch_mode(
    *,
    rows: list[dict[str, str]],
    strategy: str,
    num_simulations: int,
    overrides: dict[str, str],
    noise_scale: float,
    pdf_params: tuple[str, ...],
    msa_prefix: str,
    iqtree_threads: int,
    threads: int,
    seed: int | None,
    iqtree_exe: str,
    tool_args: str | None,
    output_dir: Path,
    command: str,
    params: dict[str, Any],
    overwrite: bool,
    resume: bool,
    dry_run: bool,
    quiet: bool,
) -> dict[str, Any]:
    run_start = _time.time()

    if dry_run:
        master_seed = seed if seed is not None else random.randint(1, 2**31 - 1)
        rng = random.Random(master_seed)
        sampled = sample_batch_rows(
            rows, strategy=strategy, n=num_simulations, rng=rng,
            pdf_params=pdf_params, noise_scale=noise_scale, overrides=overrides,
        )
        for index, row in enumerate(sampled):
            row["simulation_id"] = f"{msa_prefix}{index + 1:03d}"
            row["seed"] = master_seed + index
        payload = _assemble_batch_result(
            run_start=run_start, tool_versions={"iqtree3": "dry-run"},
            params=params, rows=rows, sampled=sampled, files=[],
            n_completed=0, n_failed=0, strategy=strategy,
            dry_run=True, plots={},
        )
        return payload

    master_seed = seed if seed is not None else random.randint(1, 2**31 - 1)
    rng = random.Random(master_seed)

    ckpt_path = output_dir / "checkpoint.json"
    step = "posttree.simulate.alisim.iqtree"

    def _verifier(task: CheckpointTask) -> bool:
        msa = Path(task.outputs.get("output_file") or "")
        if msa.exists():
            ok, _ = _validate_output(msa, "fasta")
            return ok
        return False

    completed_rows: list[dict[str, str]] = []
    files: list[dict[str, Any]] = []

    if resume:
        checkpoint = load_checkpoint(ckpt_path)
        validate_resume_params(checkpoint, params, step=step)
        completed_rows = _read_params_sampled(output_dir)
        files = []
        to_run: list[str] = []
        for task in checkpoint.tasks:
            if task.status == "success" and _verifier(task):
                files.append({
                    "simulation_id": task.task_id,
                    "status": "success",
                    "wall_time": task.outputs.get("wall_time") or 0.0,
                    "cmd": task.input,
                    "log_file": f"logs/{task.task_id}.log",
                    "output_file": f"MSAs/{task.task_id}.fa",
                })
            else:
                to_run.append(task.task_id)
        if not to_run:
            checkpoint.status = "success"
            checkpoint.completed_at = checkpoint.touch()
            save_checkpoint_atomic(checkpoint, ckpt_path)
            return _assemble_batch_result(
                run_start=run_start,
                tool_versions=_detect_iqtree_version(iqtree_exe),
                params=params, rows=rows, sampled=completed_rows,
                files=files, n_completed=len(files),
                n_failed=0, strategy=strategy, dry_run=False,
                plots=_generate_density_plots(
                    rows, completed_rows, pdf_params, output_dir,
                ),
            )
        num_simulations = len(to_run)
    else:
        if output_dir.exists() and any(output_dir.iterdir()):
            if not overwrite:
                raise ValueError(
                    f"Output directory '{output_dir}' already exists and is non-empty. "
                    "Use --overwrite to replace it."
                )
            shutil.rmtree(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        (output_dir / "MSAs").mkdir()
        (output_dir / "logs").mkdir()

        sampled = sample_batch_rows(
            rows, strategy=strategy, n=num_simulations, rng=rng,
            pdf_params=pdf_params, noise_scale=noise_scale,
            overrides=overrides,
        )
        completed_rows = []
        for index, row in enumerate(sampled):
            simulation_id = f"{msa_prefix}{index + 1:03d}"
            completed_rows.append({
                "simulation_id": simulation_id,
                "source_id": row["id"],
                "seqtype": row["seqtype"],
                "length": row["length"],
                "subs_model": row["subs_model"],
                "subs_rate": row["subs_rate"],
                "freq": row["freq"],
                "prop_inv": row["prop_inv"],
                "rate_heterogeneity": row["rate_heterogeneity"],
                "rate_categories": row["rate_categories"],
                "rate_param": row["rate_param"],
                "tree_path": row["tree_path"],
                "seed": str(master_seed + index),
            })
        _write_params_sampled(output_dir, completed_rows)

        tasks = [
            CheckpointTask(
                task_id=row["simulation_id"],
                status="pending",
                input="",
                outputs={"output_file": str(output_dir / "MSAs" / f"{row['simulation_id']}.fa")},
            )
            for row in completed_rows
        ]
        checkpoint = Checkpoint(
            schema_version=1,
            step=step,
            command=command,
            status="running",
            params_hash="",
            params=params,
            started_at=_time.strftime("%Y-%m-%dT%H:%M:%SZ", _time.gmtime()),
            updated_at="",
            completed_at=None,
            tasks=tasks,
        )
        save_checkpoint_atomic(checkpoint, ckpt_path)
        to_run = [row["simulation_id"] for row in completed_rows]

    rows_by_id = {row["simulation_id"]: row for row in completed_rows}
    msa_dir = output_dir / "MSAs"
    logs_dir = output_dir / "logs"

    from rich.progress import BarColumn, Progress, TextColumn

    def _maybe_flush() -> None:
        save_checkpoint_atomic(checkpoint, ckpt_path)

    with Progress(
        TextColumn("[bold blue]Simulating..."), BarColumn(),
        TextColumn("{task.completed}/{task.total}"), console=None,
        disable=quiet,
    ) as progress:
        task_progress = progress.add_task("simulate", total=len(to_run))
        with ProcessPoolExecutor(max_workers=threads) as pool:
            worker_args = []
            for simulation_id in to_run:
                row = rows_by_id[simulation_id]
                worker_args.append((
                    simulation_id, row["seqtype"], row["tree_path"],
                    build_model_string(row), int(row["length"]),
                    int(row["seed"]), iqtree_exe, tool_args,
                ))
            futures = {pool.submit(_run_simulation_worker, arg): arg[0] for arg in worker_args}
            for future in as_completed(futures):
                simulation_id = futures[future]
                try:
                    worker_result = future.result()
                except Exception as exc:
                    worker_result = {
                        "simulation_id": simulation_id,
                        "status": "failed",
                        "wall_time": 0.0,
                        "cmd": [],
                        "log_file": f"logs/{simulation_id}.log",
                        "log_text": f"worker error: {exc}",
                        "output_file": f"MSAs/{simulation_id}.fa",
                        "output_path": None,
                        "reason": str(exc),
                    }

                if worker_result["status"] == "success":
                    output_path = worker_result.get("output_path")
                    ok = False
                    if output_path and output_path.exists():
                        target = msa_dir / f"{simulation_id}.fa"
                        shutil.move(str(output_path), str(target))
                        _wrap_fasta(target)
                        ok, _warnings = _validate_output(target, "fasta")
                    if ok:
                        files.append({
                            "simulation_id": simulation_id,
                            "status": "success",
                            "wall_time": worker_result["wall_time"],
                            "cmd": worker_result["cmd"],
                            "log_file": worker_result["log_file"],
                            "output_file": worker_result["output_file"],
                        })
                        mark = "success"
                    else:
                        worker_result["status"] = "failed"
                        worker_result["reason"] = worker_result.get("reason") or "generated MSA failed validation"
                        mark = "failed"
                else:
                    mark = "failed"

                log_text = worker_result.get("log_text") or ""
                if log_text:
                    (logs_dir / f"{simulation_id}.log").write_text(log_text)

                for task in checkpoint.tasks:
                    if task.task_id == simulation_id:
                        task.status = mark
                        task.reason = worker_result.get("reason")
                        task.input = " ".join(worker_result.get("cmd") or [])
                        task.outputs["wall_time"] = str(worker_result["wall_time"])
                        task.outputs["log_file"] = f"logs/{simulation_id}.log"
                        task.updated_at = _time.strftime("%Y-%m-%dT%H:%M:%SZ", _time.gmtime())
                _maybe_flush()
                progress.update(task_progress, advance=1)

    checkpoint.status = "success"
    checkpoint.completed_at = checkpoint.touch()
    save_checkpoint_atomic(checkpoint, ckpt_path)

    plots = _generate_density_plots(rows, completed_rows, pdf_params, output_dir)
    n_failed = len(to_run) - len([f for f in files if f["status"] == "success"])
    return _assemble_batch_result(
        run_start=run_start,
        tool_versions=_detect_iqtree_version(iqtree_exe),
        params=params, rows=rows, sampled=completed_rows, files=files,
        n_completed=len(files), n_failed=n_failed, strategy=strategy,
        dry_run=False, plots=plots,
    )


def _assemble_batch_result(
    *,
    run_start: float,
    tool_versions: dict[str, str],
    params: dict[str, Any],
    rows: list[dict[str, str]],
    sampled: list[dict[str, str]],
    files: list[dict[str, Any]],
    n_completed: int,
    n_failed: int,
    strategy: str,
    dry_run: bool,
    plots: dict[str, Any],
) -> dict[str, Any]:
    output_files: dict[str, Any] = {
        "msas_dir": {
            "path": str(Path(params["output_dir"]) / "MSAs"),
            "description": "Simulated MSA files (one per simulation)",
        },
        "params_sampled": {
            "path": str(Path(params["output_dir"]) / "params_sampled.tsv"),
            "description": "Actual parameters used per simulation (TSV)",
        },
    }
    output_files.update(plots)
    return {
        "status": "success",
        "command": params["_command"],
        "wall_time": round(_time.time() - run_start, 3),
        "tool_versions": tool_versions,
        "params": {k: v for k, v in params.items() if not k.startswith("_")},
        "key_results": {
            "n_simulations_requested": params["num_simulations"],
            "n_simulations_completed": n_completed,
            "n_simulations_failed": n_failed,
            "strategy": strategy,
            "source_loci": len(rows),
        },
        "error": None,
        "error_category": None,
        "data": {
            "output_files": output_files,
            "sampled_rows": sampled if dry_run else [],
            "files": files,
        },
    }


# ===================================================================
# Entry point
# ===================================================================

def _parse_overrides(override: str | None) -> dict[str, str]:
    if not override:
        return {}
    result: dict[str, str] = {}
    for segment in override.split(","):
        segment = segment.strip()
        if not segment:
            continue
        if "=" not in segment:
            raise ValueError(
                f"invalid --override segment {segment!r}; expected key=value"
            )
        key, _, value = segment.partition("=")
        key = key.strip()
        if key not in VALID_OVERRIDE_KEYS:
            raise ValueError(
                f"invalid --override key {key!r}; valid keys: {', '.join(sorted(VALID_OVERRIDE_KEYS))}"
            )
        result[key] = value.strip()
    return result


def run_alisim_iqtree(
    *,
    ref_tree: Path | None = None,
    model: str | None = None,
    model_partitions: Path | None = None,
    seq_type: str | None = None,
    length: int | None = None,
    model_params: Path | None = None,
    strategy: str | None = None,
    num_simulations: int | None = None,
    override: str | None = None,
    noise_scale: float = 1.0,
    pdf_params: str = "length,prop_inv,rate_param",
    msa_prefix: str = "sim",
    out_format: str = "fasta",
    num_alignments: int = 1,
    iqtree_threads: int = 1,
    threads: int = 4,
    seed: int | None = None,
    iqtree_path: str | None = None,
    tool_args: str | None = None,
    output_dir: Path = Path("runs/posttree/simulate/alisim/iqtree"),
    overwrite: bool = False,
    resume: bool = False,
    dry_run: bool = False,
    quiet: bool = False,
) -> dict[str, Any]:
    """Run one AliSim invocation or a resumable batch of one-alignment invocations.

    Returns a result.json payload dict.  Raises ValueError on validation errors.
    """
    run_start = _time.time()

    if overwrite and resume:
        raise ValueError("--overwrite and --resume are mutually exclusive.")

    batch_mode = model_params is not None
    single_required = {"ref_tree": ref_tree, "model": model,
                       "model_partitions": model_partitions, "seq_type": seq_type,
                       "length": length}
    if batch_mode and any(value is not None for value in single_required.values()):
        raise ValueError(
            "--model-params is mutually exclusive with --ref-tree, --model, "
            "--model-partitions, --seq-type, and --length (batch mode vs single mode)."
        )

    output_dir = output_dir.resolve()

    if batch_mode:
        if strategy is None:
            raise ValueError("--strategy is required in batch mode (--model-params).")
        if strategy not in {"complete", "mixed", "pdf"}:
            raise ValueError(f"invalid --strategy {strategy!r}; expected complete, mixed, or pdf")
        if num_simulations is None:
            raise ValueError("--num-simulations is required in batch mode (--model-params).")
        if num_simulations < 1:
            raise ValueError(f"--num-simulations must be >= 1, got {num_simulations}")
        if resume and dry_run:
            raise ValueError("--resume and --dry-run are mutually exclusive.")
        if num_alignments != 1:
            raise ValueError("--num-alignments is single-parameter mode only.")
        if strategy != "pdf":
            if noise_scale != 1.0 or pdf_params != "length,prop_inv,rate_param":
                raise ValueError("--noise-scale and --pdf-params require --strategy pdf.")
        overrides = _parse_overrides(override)
    else:
        if ref_tree is None:
            raise ValueError("--ref-tree is required in single mode.")
        if not ref_tree.exists():
            raise ValueError(f"--ref-tree does not exist: {ref_tree}")
        if model is None and model_partitions is None:
            raise ValueError("--model or --model-partitions (exactly one) is required.")
        if model is not None and model_partitions is not None:
            raise ValueError("--model and --model-partitions are mutually exclusive.")
        if seq_type is None:
            raise ValueError("--seq-type is required in single mode.")
        if seq_type not in {"AA", "DNA"}:
            raise ValueError(f"--seq-type must be AA or DNA, got {seq_type!r}")
        if model_partitions is None:
            if length is None:
                raise ValueError("--length is required unless --model-partitions is used.")
            if length < 1:
                raise ValueError(f"--length must be >= 1, got {length}")
        if resume:
            raise ValueError("--resume is batch mode only.")
        if num_alignments < 1:
            raise ValueError(f"--num-alignments must be >= 1, got {num_alignments}")
        strategy = "single"

    if out_format not in OUT_FORMAT_EXT:
        raise ValueError(f"--out-format must be fasta or phy, got {out_format!r}")
    if iqtree_threads < 1:
        raise ValueError(f"--iqtree-threads must be >= 1, got {iqtree_threads}")
    if threads < 1:
        raise ValueError(f"--threads must be >= 1, got {threads}")
    if not (0.0 <= noise_scale <= 1.0):
        raise ValueError(f"--noise-scale must be in [0.0, 1.0], got {noise_scale}")

    pdf_param_tuple = tuple(p.strip() for p in pdf_params.split(",") if p.strip())
    for param in pdf_param_tuple:
        if param not in VALID_PDF_PARAMS:
            raise ValueError(
                f"invalid --pdf-params value {param!r}; "
                f"valid values: {', '.join(sorted(VALID_PDF_PARAMS))}"
            )

    if tool_args:
        _check_managed_flag_conflict(tool_args)

    resolved_seed = seed if seed is not None else random.randint(1, 2**31 - 1)

    try:
        iqtree_exe = _resolve_iqtree_path(iqtree_path, dry_run)
    except (ValueError, FileNotFoundError) as exc:
        raise ValueError(str(exc)) from exc

    if batch_mode:
        rows = load_params_table(model_params)

    _command_parts = ["phyloai", "posttree", "simulate", "alisim", "iqtree"]
    if model_params:
        _command_parts += ["--model-params", str(model_params)]
    if strategy and batch_mode:
        _command_parts += ["--strategy", strategy]
    if num_simulations and batch_mode:
        _command_parts += ["--num-simulations", str(num_simulations)]
    if override:
        _command_parts += ["--override", override]
    if batch_mode:
        _command_parts += ["--seed", str(resolved_seed), "--output-dir", str(output_dir)]
    else:
        _command_parts += ["--ref-tree", str(ref_tree)]
        if model:
            _command_parts += ["--model", model]
        if model_partitions:
            _command_parts += ["--model-partitions", str(model_partitions)]
        _command_parts += ["--seq-type", seq_type]
        if length is not None:
            _command_parts += ["--length", str(length)]
        _command_parts += ["--out-format", out_format, "--num-alignments", str(num_alignments)]
        if seed is not None:
            _command_parts += ["--seed", str(seed)]
        _command_parts += ["--output-dir", str(output_dir)]
    full_command = " ".join(_command_parts)

    params: dict[str, Any] = {
        "model_params": str(model_params.resolve()) if model_params else None,
        "strategy": strategy,
        "num_simulations": num_simulations,
        "override": override,
        "noise_scale": noise_scale,
        "pdf_params": pdf_params,
        "msa_prefix": msa_prefix,
        "out_format": out_format,
        "iqtree_threads": iqtree_threads,
        "threads": threads,
        "seed": resolved_seed,
        "output_dir": str(output_dir),
        "iqtree_path": iqtree_path,
        "tool_args": tool_args,
        "overwrite": overwrite,
        "resume": resume,
        "dry_run": dry_run,
        "quiet": quiet,
    }
    if batch_mode:
        params.update({
            "ref_tree": None,
            "model": None,
            "model_partitions": None,
            "seq_type": None,
            "length": None,
            "num_alignments": num_alignments,
        })
    else:
        params.update({
            "ref_tree": str(ref_tree.resolve()),
            "model": model,
            "model_partitions": str(model_partitions.resolve()) if model_partitions else None,
            "seq_type": seq_type,
            "length": length,
            "num_alignments": num_alignments,
        })
    params["_command"] = full_command

    if batch_mode:
        payload = _run_batch_mode(
            rows=rows,
            strategy=strategy,
            num_simulations=num_simulations,
            overrides=overrides,
            noise_scale=noise_scale,
            pdf_params=pdf_param_tuple,
            msa_prefix=msa_prefix,
            iqtree_threads=iqtree_threads,
            threads=threads,
            seed=resolved_seed,
            iqtree_exe=iqtree_exe,
            tool_args=tool_args,
            output_dir=output_dir,
            command=full_command,
            params=params,
            overwrite=overwrite,
            resume=resume,
            dry_run=dry_run,
            quiet=quiet,
        )
        if not dry_run:
            with open(output_dir / "result.json", "w", encoding="utf-8") as fh:
                json.dump(payload, fh, indent=2)
        return payload

    if not dry_run:
        if output_dir.exists() and any(output_dir.iterdir()):
            raise ValueError(
                f"Output directory '{output_dir}' already exists and is non-empty. "
                "Use --overwrite to replace it."
            )
        output_dir.mkdir(parents=True, exist_ok=True)

    data, n_generated = _run_single_mode(
        ref_tree=ref_tree,
        model=model,
        model_partitions=model_partitions,
        seq_type=seq_type,
        length=length,
        msa_prefix=msa_prefix,
        num_alignments=num_alignments,
        out_format=out_format,
        iqtree_threads=iqtree_threads,
        seed=resolved_seed,
        iqtree_exe=iqtree_exe,
        tool_args=tool_args,
        output_dir=output_dir,
        dry_run=dry_run,
    )

    payload: dict[str, Any] = {
        "status": "success",
        "command": full_command,
        "wall_time": round(_time.time() - run_start, 3),
        "tool_versions": (
            {"iqtree3": "dry-run"} if dry_run else _detect_iqtree_version(iqtree_exe)
        ),
        "params": {k: v for k, v in params.items() if not k.startswith("_")},
        "key_results": {
            "n_msas_generated": n_generated,
            "seq_type": seq_type,
            "length": length,
            "model": model or str(model_partitions),
        },
        "error": None,
        "error_category": None,
        "data": data,
    }

    if not dry_run:
        with open(output_dir / "result.json", "w", encoding="utf-8") as fh:
            json.dump(payload, fh, indent=2)

    return payload
