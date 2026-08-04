"""Model adequacy checks from an observed MSA against simulated replicates.

Pure-Python implementation of the PhyloBayes AllPostPred summary statistics
(PPA-DIV, PPA-CONV, PPA-VAR, PPA-COMP) described in the adequacy design spec.
"""

from __future__ import annotations

import csv
import json
import math
import shlex
import shutil
import statistics
import time as _time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Any, Callable

from Bio.Align import MultipleSeqAlignment

from phyloai.core.checkpoint import (
    Checkpoint,
    CheckpointTask,
    canonical_params_hash,
    load_checkpoint,
    save_checkpoint_atomic,
    validate_resume_params,
)
from phyloai.core.formats import FormatConverter
from phyloai.core.sequence_normalization import detect_seq_type

AA_STATES = "ACDEFGHIKLMNPQRSTVWY"
NT_STATES = "ACGT"
SCALAR_NAMES = ("div", "siteconvprob", "sitecomp", "comp_max", "comp_mean")


class PreflightError(ValueError):
    """A pre-flight refusal raised before the output directory is claimed."""


_FORMAT_CONVERTER = FormatConverter()


def _read_alignment(path: Path) -> MultipleSeqAlignment:
    try:
        alignment = _FORMAT_CONVERTER.read(path)
    except Exception as exc:
        raise ValueError(f"unable to parse alignment file {path}: {exc}") from exc
    if not len(alignment) or not alignment.get_alignment_length():
        raise ValueError(f"alignment file {path} is empty")
    return alignment


def _alignment_ids(alignment: MultipleSeqAlignment, label: str) -> list[str]:
    names = [record.id for record in alignment]
    duplicates = sorted({name for name in names if names.count(name) > 1})
    if duplicates:
        raise ValueError(f"{label} alignment has duplicate taxon IDs: {', '.join(duplicates)}")
    if len({len(record.seq) for record in alignment}) != 1:
        raise ValueError(f"{label} alignment sequences have unequal lengths")
    return names


def _process_simulation(path: Path, original_ids: list[str], original_length: int, seq_type: str) -> dict[str, Any]:
    alignment = _read_alignment(path)
    ids = _alignment_ids(alignment, "simulated")
    if set(ids) != set(original_ids):
        raise ValueError("taxon name mismatch between original and simulated MSAs")
    if alignment.get_alignment_length() != original_length:
        raise ValueError("length mismatch between original and simulated MSAs")
    records = {record.id: record for record in alignment}
    ordered = MultipleSeqAlignment([records[name] for name in original_ids])
    return _compute_statistics(ordered, seq_type)


def _resolved_seq_type(original: MultipleSeqAlignment, requested: str) -> str:
    normalized = requested.upper()
    if normalized == "AUTO":
        return detect_seq_type([str(record.seq) for record in original])
    if normalized not in {"AA", "NT"}:
        raise ValueError(f"invalid seq_type: {requested!r}")
    return normalized


def _fingerprint(path: Path) -> str:
    stat = path.stat()
    return f"{path.resolve()}|{stat.st_size}|{stat.st_mtime_ns}"


def _task_outputs(stats: dict[str, Any]) -> dict[str, str]:
    return {
        **{name: repr(float(stats[name])) for name in SCALAR_NAMES},
        "taxon_dist_j": json.dumps({name: repr(float(value)) for name, value in stats["taxon_dist_j"].items()}, sort_keys=True),
    }


def _task_stats(task: CheckpointTask) -> dict[str, Any]:
    return {
        **{name: float(task.outputs[name]) for name in SCALAR_NAMES},
        "taxon_dist_j": {name: float(value) for name, value in json.loads(task.outputs["taxon_dist_j"]).items()},
    }


def _write_csv(path: Path, fieldnames: list[str], rows: list[dict[str, Any]], delimiter: str = ",") -> None:
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="raise", delimiter=delimiter)
        writer.writeheader()
        writer.writerows(rows)


def _list_alignment_files(directory: Path) -> list[Path]:
    """Sorted regular non-empty files in a directory (non-recursive)."""
    if not directory.exists():
        raise PreflightError(f"--simulated-dir does not exist: {directory}")
    if not directory.is_dir():
        raise PreflightError(f"--simulated-dir is not a directory: {directory}")
    files = sorted(
        path for path in directory.iterdir()
        if path.is_file() and path.stat().st_size > 0
    )
    if not files:
        raise PreflightError(f"--simulated-dir contains no non-empty files: {directory}")
    return files


def _reconcile_tasks(
    checkpoint: Checkpoint,
    disk_files: list[Path],
    warnings: list[str],
) -> list[CheckpointTask]:
    """Reconcile checkpoint tasks against the on-disk simulated files.

    Newly discovered files become pending tasks; stale success tasks whose
    fingerprint changed are reset to pending; checkpoint tasks whose file no
    longer exists are dropped with a warning.
    """
    disk_ids = {str(path.resolve()): path for path in disk_files}
    by_id: dict[str, CheckpointTask] = {task.task_id: task for task in checkpoint.tasks}
    for task_id in list(by_id):
        if task_id not in disk_ids:
            warnings.append(f"dropping checkpoint task {task_id}: simulated file no longer present")
            del by_id[task_id]
    for path in disk_files:
        task_id = str(path.resolve())
        fingerprint = _fingerprint(path)
        task = by_id.get(task_id)
        if task is None:
            by_id[task_id] = CheckpointTask(task_id=task_id, status="pending", input=fingerprint)
        elif task.status == "success" and task.input == fingerprint:
            continue
        else:
            task.status = "pending"
            task.input = fingerprint
            task.reason = None
    return [by_id[task_id] for task_id in sorted(by_id)]


def run_simulate_adequacy(
    original_msa: Path,
    simulated_dir: Path,
    seq_type: str = "auto",
    threads: int = 4,
    table_format: str = "csv",
    output_dir: Path = Path("runs/posttree/simulate/adequacy"),
    overwrite: bool = False,
    resume: bool = False,
    dry_run: bool = False,
    quiet: bool = False,
    progress_callback: Callable[[int, int], None] | None = None,
) -> dict[str, Any]:
    """Assess model adequacy from an observed MSA and simulated replicates.

    Returns the standard result.json payload. Raises ``PreflightError`` for
    output-lifecycle refusals (before the output directory is claimed) and
    ``ValueError`` for ordinary validation failures after the directory is
    claimed.
    """
    run_start = _time.time()

    if overwrite and resume:
        raise PreflightError("--overwrite and --resume are mutually exclusive.")
    if table_format not in {"csv", "tsv"}:
        raise PreflightError(f"invalid table_format: {table_format!r}")
    if seq_type.upper() not in {"AA", "NT", "AUTO"}:
        raise PreflightError(f"invalid seq_type: {seq_type!r}")

    if not original_msa.exists():
        raise PreflightError(f"--original-msa does not exist: {original_msa}")
    if not original_msa.is_file():
        raise PreflightError(f"--original-msa is not a file: {original_msa}")

    resume_params = {
        "original_msa": str(original_msa.resolve()),
        "simulated_dir": str(simulated_dir.resolve()),
        "seq_type": seq_type,
    }
    output_dir = output_dir.resolve()
    ckpt_path = output_dir / "checkpoint.json"
    step = "posttree.simulate.adequacy"

    command_parts = [
        "phyloai", "posttree", "simulate", "adequacy",
        "--original-msa", str(original_msa.resolve()),
        "--simulated-dir", str(simulated_dir.resolve()),
        "--seq-type", seq_type,
        "--threads", str(threads),
        "-o", str(output_dir.resolve()),
    ]
    if table_format != "csv":
        command_parts.extend(["--table-format", table_format])
    if overwrite:
        command_parts.append("--overwrite")
    if resume:
        command_parts.append("--resume")
    if dry_run:
        command_parts.append("--dry-run")
    if quiet:
        command_parts.append("--quiet")
    full_command = shlex.join(command_parts)

    original_fingerprint = _fingerprint(original_msa)
    disk_files = _list_alignment_files(simulated_dir)

    checkpoint: Checkpoint | None = None
    if dry_run:
        pass
    elif resume:
        try:
            checkpoint = load_checkpoint(ckpt_path)
        except (FileNotFoundError, ValueError) as exc:
            raise PreflightError(str(exc)) from exc
        try:
            validate_resume_params(checkpoint, resume_params, step=step)
        except ValueError as exc:
            raise PreflightError(str(exc)) from exc
        stored_fingerprint = checkpoint.original_msa_fingerprint
        if stored_fingerprint is None:
            raise PreflightError(
                "checkpoint has no original_msa fingerprint; use --overwrite to start fresh"
            )
        if stored_fingerprint != original_fingerprint:
            raise PreflightError(
                "original MSA has changed since the run started; use --overwrite to start fresh"
            )
    elif output_dir.exists() and any(output_dir.iterdir()):
        if not overwrite:
            raise PreflightError(
                f"Output directory '{output_dir}' already exists and is non-empty. "
                "Use --overwrite to replace it."
            )

    original = _read_alignment(original_msa)
    original_ids = _alignment_ids(original, "original")
    resolved_seq_type = _resolved_seq_type(original, seq_type)
    original_length = original.get_alignment_length()
    original_stats = _compute_statistics(original, resolved_seq_type)
    payload_params = {
        "original_msa": str(original_msa.resolve()),
        "simulated_dir": str(simulated_dir.resolve()),
        "seq_type": seq_type,
        "detected_seq_type": resolved_seq_type,
        "threads": threads,
        "table_format": table_format,
        "output_dir": str(output_dir),
        "overwrite": overwrite,
        "resume": resume,
        "dry_run": dry_run,
        "quiet": quiet,
    }

    if dry_run:
        payload: dict[str, Any] = {
            "status": "success",
            "command": full_command,
            "wall_time": round(_time.time() - run_start, 3),
            "tool_versions": {},
            "params": payload_params,
            "key_results": {
                "seq_type": resolved_seq_type,
                "n_taxa": len(original_ids),
                "n_sites": original_length,
                "n_simulated_files": len(disk_files),
            },
            "error": None,
            "error_category": None,
            "data": {"cmd": [], "tool_stderr": "", "warnings": [], "output_files": {}},
        }
        if not quiet:
            _print_summary(payload, dry_run=True)
        return payload

    warnings: list[str] = []

    if resume:
        assert checkpoint is not None
    else:
        if overwrite and output_dir.exists():
            shutil.rmtree(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        checkpoint = Checkpoint(
            schema_version=1,
            step=step,
            command=full_command,
            status="running",
            params_hash=canonical_params_hash(resume_params),
            params=resume_params,
            started_at=_time.strftime("%Y-%m-%dT%H:%M:%SZ", _time.gmtime()),
            updated_at=_time.strftime("%Y-%m-%dT%H:%M:%SZ", _time.gmtime()),
            completed_at=None,
            tasks=[
                CheckpointTask(
                    task_id=str(path.resolve()),
                    status="pending",
                    input=_fingerprint(path),
                )
                for path in disk_files
            ],
            original_msa_fingerprint=original_fingerprint,
        )
        save_checkpoint_atomic(checkpoint, ckpt_path)

    checkpoint.tasks = _reconcile_tasks(checkpoint, disk_files, warnings)
    save_checkpoint_atomic(checkpoint, ckpt_path)

    pending = [task for task in checkpoint.tasks if task.status == "pending"]
    if pending:
        with ProcessPoolExecutor(max_workers=threads) as pool:
            futures = {
                pool.submit(
                    _process_simulation,
                    Path(task.task_id),
                    original_ids,
                    original_length,
                    resolved_seq_type,
                ): task
                for task in pending
            }
            completed = 0
            total = len(pending)
            for future in as_completed(futures):
                task = futures[future]
                try:
                    stats = future.result()
                except Exception as exc:
                    task.status = "failed"
                    task.reason = str(exc)
                else:
                    task.outputs = _task_outputs(stats)
                    task.status = "success"
                    task.reason = None
                task.updated_at = _time.strftime("%Y-%m-%dT%H:%M:%SZ", _time.gmtime())
                completed += 1
                save_checkpoint_atomic(checkpoint, ckpt_path)
                if progress_callback is not None:
                    progress_callback(completed, total)

    success_tasks = [task for task in checkpoint.tasks if task.status == "success"]
    if len(success_tasks) < 10:
        raise ValueError(
            "at least 10 valid simulated MSAs are required; "
            f"only {len(success_tasks)} of {len(checkpoint.tasks)} succeeded"
        )

    stats_list = [_task_stats(task) for task in success_tasks]

    summary_rows: list[dict[str, Any]] = []
    statistics_kr: dict[str, Any] = {}
    for name in SCALAR_NAMES:
        direction = "div" if name == "div" else "high"
        obs_value = float(original_stats[name])
        summary = _summarize_distribution(
            [float(stats[name]) for stats in stats_list],
            obs_value,
            direction,
        )
        summary_row = {
            "statistic": name,
            "obs": obs_value,
            "mean_sim": summary["mean_sim"],
            "sd_sim": summary["sd_sim"],
            "ci_lower": summary["ci_lower"],
            "ci_upper": summary["ci_upper"],
            "z_score": summary["z_score"],
            "pp": "" if summary["pp"] is None else summary["pp"],
            "n_simulations": len(success_tasks),
        }
        summary_rows.append(summary_row)
        summary = {"obs": obs_value, **summary}
        if name in {"comp_max", "comp_mean"}:
            statistics_kr.setdefault("comp", {})[name.replace("comp_", "")] = summary
        else:
            statistics_kr[name] = summary

    taxon_rows: list[dict[str, Any]] = []
    for taxon in original_ids:
        values = [float(stats["taxon_dist_j"][taxon]) for stats in stats_list]
        obs_value = float(original_stats["taxon_dist_j"][taxon])
        summary = _summarize_distribution(values, obs_value, "high")
        taxon_rows.append({
            "taxon": taxon,
            "obs": obs_value,
            "mean_pred": summary["mean_sim"],
            "sd_pred": summary["sd_sim"],
            "ci_lower": summary["ci_lower"],
            "ci_upper": summary["ci_upper"],
            "z_score": summary["z_score"],
            "pp": "" if summary["pp"] is None else summary["pp"],
        })

    per_sim_rows: list[dict[str, Any]] = []
    for task in success_tasks:
        stats = _task_stats(task)
        per_sim_rows.append({
            "file": Path(task.task_id).name,
            **{name: float(stats[name]) for name in SCALAR_NAMES},
        })

    delimiter = "\t" if table_format == "tsv" else ","
    suffix = ".tsv" if table_format == "tsv" else ".csv"
    summary_path = output_dir / f"adequacy_summary{suffix}"
    taxon_path = output_dir / f"adequacy_taxon_comp{suffix}"
    per_sim_path = output_dir / f"per_simulation_stats{suffix}"

    if not dry_run:
        _write_csv(summary_path, ["statistic", "obs", "mean_sim", "sd_sim", "ci_lower", "ci_upper", "z_score", "pp", "n_simulations"], summary_rows, delimiter=delimiter)
        _write_csv(taxon_path, ["taxon", "obs", "mean_pred", "sd_pred", "ci_lower", "ci_upper", "z_score", "pp"], taxon_rows, delimiter=delimiter)
        _write_csv(per_sim_path, ["file", *SCALAR_NAMES], per_sim_rows, delimiter=delimiter)

        checkpoint.status = "success"
        checkpoint.completed_at = checkpoint.touch()
        save_checkpoint_atomic(checkpoint, ckpt_path)

    failed_tasks = [task for task in checkpoint.tasks if task.status == "failed"]
    n_failed = len(failed_tasks)
    warnings.extend(
        f"skipping simulated file {Path(task.task_id).name}: {task.reason}"
        for task in failed_tasks
    )

    payload = {
        "status": "success",
        "command": full_command,
        "wall_time": round(_time.time() - run_start, 3),
        "tool_versions": {},
        "params": payload_params,
        "key_results": {
            "n_simulations": len(success_tasks),
            "n_failed": n_failed,
            "seq_type": resolved_seq_type,
            "n_taxa": len(original_ids),
            "n_sites": original_length,
            "statistics": statistics_kr,
        },
        "error": None,
        "error_category": None,
        "data": {
            "cmd": [],
            "tool_stderr": "",
            "warnings": warnings,
            "output_files": {
                "adequacy_summary": {
                    "path": str(summary_path),
                    "description": "Model adequacy summary: obs, mean_sim, sd_sim, CI, z-score, pp for 5 statistics",
                },
                "adequacy_taxon_comp": {
                    "path": str(taxon_path),
                    "description": "Per-taxon PPA-COMP statistics: obs, mean_pred, sd_pred, CI, z-score, pp",
                },
                "per_simulation_stats": {
                    "path": str(per_sim_path),
                    "description": "Raw statistic values for all simulated replicates (null distribution)",
                },
            },
        },
    }

    with open(output_dir / "result.json", "w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2)

    if not quiet:
        _print_summary(payload, dry_run=False)

    return payload


def _print_summary(payload: dict[str, Any], dry_run: bool) -> None:
    click_echo = __import__("click").echo
    kr = payload["key_results"]
    click_echo("Assessing model adequacy...")
    click_echo(f"  Sequence type: {kr['seq_type']}")
    click_echo(f"  Taxa: {kr['n_taxa']}")
    click_echo(f"  Sites: {kr['n_sites']}")
    if dry_run:
        click_echo(f"  Simulated files found: {kr.get('n_simulated_files', 0)}")
        click_echo(f"  Output directory: {payload['params']['output_dir']}")
        click_echo("Dry run: no files written.")
        return
    click_echo(f"  Valid simulations: {kr['n_simulations']}")
    click_echo(f"  Failed files: {kr['n_failed']}")
    click_echo("  Outputs:")
    for label, info in payload["data"]["output_files"].items():
        click_echo(f"    {label}: {info['path']}")
    click_echo("Result written to " + str(Path(payload["params"]["output_dir"]) / "result.json"))


def _compute_statistics(alignment: MultipleSeqAlignment, seq_type: str) -> dict[str, Any]:
    states = AA_STATES if seq_type == "AA" else NT_STATES
    state_index = {state: index for index, state in enumerate(states)}
    n_states = len(states)
    names = [record.id for record in alignment]
    sequences = [str(record.seq).upper() for record in alignment]
    n_sites = alignment.get_alignment_length()

    site_freqs: list[list[float]] = []
    diversity_total = 0.0
    squared_freq_total = 0.0
    for site in range(n_sites):
        counts = [0] * n_states
        for sequence in sequences:
            index = state_index.get(sequence[site])
            if index is not None:
                counts[index] += 1
        observed = sum(counts)
        if not observed:
            continue
        freqs = [count / observed for count in counts]
        site_freqs.append(freqs)
        diversity_total += sum(count > 0 for count in counts)
        squared_freq_total += sum(freq * freq for freq in freqs)

    if not site_freqs:
        raise ValueError("alignment has no informative sites")
    n_informative = len(site_freqs)
    means = [sum(freq[k] for freq in site_freqs) / n_informative for k in range(n_states)]
    sitecomp = sum(
        sum(freq[k] * freq[k] for freq in site_freqs) / n_informative - means[k] * means[k]
        for k in range(n_states)
    ) / n_states

    taxon_freqs: dict[str, list[float]] = {}
    for name, sequence in zip(names, sequences):
        counts = [0] * n_states
        for state in sequence:
            index = state_index.get(state)
            if index is not None:
                counts[index] += 1
        total = sum(counts)
        if not total:
            raise ValueError(f"taxon {name!r} has no valid characters")
        taxon_freqs[name] = [count / total for count in counts]
    global_freq = [sum(freq[k] for freq in taxon_freqs.values()) / len(taxon_freqs) for k in range(n_states)]
    taxon_dist = {
        name: sum((freq[k] - global_freq[k]) ** 2 for k in range(n_states))
        for name, freq in taxon_freqs.items()
    }
    return {
        "div": diversity_total / n_informative,
        "siteconvprob": squared_freq_total / n_informative,
        "sitecomp": sitecomp,
        "comp_max": max(taxon_dist.values()),
        "comp_mean": sum(taxon_dist.values()) / len(taxon_dist),
        "taxon_dist_j": taxon_dist,
        "n_informative_sites": n_informative,
    }


def _summarize_distribution(values: list[float], obs: float, direction: str) -> dict[str, float | int | None]:
    if len(values) < 10:
        raise ValueError("at least 10 valid simulated MSAs are required")
    mean_sim = sum(values) / len(values)
    sd_sim = math.sqrt(sum(value * value for value in values) / len(values) - mean_sim * mean_sim)
    quantiles = statistics.quantiles(values, n=40, method="inclusive")
    ci_lower, ci_upper = quantiles[0], quantiles[-1]
    if sd_sim == 0:
        return {"mean_sim": mean_sim, "sd_sim": 0.0, "ci_lower": mean_sim, "ci_upper": mean_sim, "z_score": 0.0, "pp": None}
    if direction == "div":
        z_score = (mean_sim - obs) / sd_sim
        pp = sum(value <= obs for value in values) / len(values)
    else:
        z_score = (obs - mean_sim) / sd_sim
        pp = sum(value > obs for value in values) / len(values)
    return {"mean_sim": mean_sim, "sd_sim": sd_sim, "ci_lower": ci_lower, "ci_upper": ci_upper, "z_score": z_score, "pp": pp}
