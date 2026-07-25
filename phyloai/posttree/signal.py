"""Phylogenetic signal distribution analysis: lnl, consistent, fclm subcommands."""
from __future__ import annotations

import csv
import json
import os
import re
import shlex
import shutil
import subprocess
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Any

import numpy as np
from scipy.stats import mannwhitneyu

from phyloai.core.formats import FormatConverter
from phyloai.core.iqtree import _detect_iqtree_version, _resolve_iqtree_path
from phyloai.core.file_matching import logical_tree_locus_candidates, scan_tree_dir

_PARTITION_RE = re.compile(r"^\s*([^,]+)\s*,\s*(.+?)\s*=\s*(\d+)\s*-\s*(\d+)\s*$")
_NEXUS_CHARSET_RE = re.compile(r"^\s*charset\s+(\S+)\s*=\s*(\d+)\s*-\s*(\d+)\s*;?\s*$", re.IGNORECASE)
_FLOAT_TOL = 1e-9
_LNL_BLOCKED_FLAGS = frozenset({"-s", "-z", "-wslr", "--prefix", "-p", "-Q"})
_FCLM_BLOCKED_FLAGS = frozenset({"-s", "-lmap", "-lmclust", "-n", "-p", "-Q", "--prefix"})
_CONSISTENT_BLOCKED_FLAGS = frozenset({"-s", "-z", "-wslr", "-p", "-Q", "--prefix"})


def _parse_partition_ranges(path: Path) -> list[dict[str, Any]]:
    """Parse RAxML-like or NEXUS (charset-only) partitions into 1-based ranges.

    RAxML style: ``LG, geneA = 1-235``.  NEXUS style: ``charset geneA = 1-235;``.
    NEXUS ``charpartition`` lines (model assignments) are ignored — this function
    only extracts locus-name / start-end boundaries.
    """
    records: list[dict[str, Any]] = []
    with open(path) as handle:
        first_line = handle.readline().strip()
        is_nexus = first_line.upper().startswith("#NEXUS")
        handle.seek(0)
        for lineno, raw in enumerate(handle, 1):
            line = raw.strip()
            if not line:
                continue
            if is_nexus:
                lower = line.lower()
                if line.startswith("[") or lower.startswith("#nexus") or lower.startswith("begin") or lower.startswith("end;") or lower.startswith("charpartition"):
                    continue
                match = _NEXUS_CHARSET_RE.match(line)
                if match is None:
                    continue
                locus, start, end = match.groups()
            else:
                match = _PARTITION_RE.match(line)
                if match is None:
                    raise ValueError(f"Unparseable partition line {lineno}: {raw.rstrip()}")
                _, locus, start, end = match.groups()
            start_int, end_int = int(start), int(end)
            if start_int < 1 or end_int < start_int:
                raise ValueError(f"Invalid range line {lineno}: {start_int}-{end_int}")
            records.append({"locus": locus.strip(), "start": start_int, "end": end_int})
    if not records:
        raise ValueError(f"Partition file is empty: {path}")
    return records


def _parse_sitelh(path: Path) -> tuple[list[str], list[list[float]]]:
    """Parse an IQ-TREE .sitelh file into labels and per-tree site scores."""
    with open(path) as handle:
        lines = [line for line in handle if line.strip()]
    header = lines[0].split()
    n_trees, n_sites = int(header[0]), int(header[1])
    labels: list[str] = []
    scores: list[list[float]] = []
    for line in lines[1:]:
        parts = line.split()
        if len(parts) != n_sites + 1:
            raise ValueError(f"Expected {n_sites + 1} columns in sitelh line, got {len(parts)}")
        labels.append(parts[0])
        scores.append([float(value) for value in parts[1:]])
    if len(labels) != n_trees:
        raise ValueError(f"Expected {n_trees} tree rows, got {len(labels)}")
    return labels, scores


def _sum_gene_lnl(site_scores: list[list[float]], start: int, end: int) -> list[float]:
    """Sum site lnL values for a 1-based inclusive locus range per tree."""
    return [sum(scores[start - 1:end]) for scores in site_scores]


def _delta_score(scores: list[float]) -> float:
    """Return signed two-tree delta or mean absolute pairwise score difference."""
    if len(scores) == 2:
        return scores[0] - scores[1]
    differences = [abs(first - second) for i, first in enumerate(scores) for second in scores[i + 1:]]
    return sum(differences) / len(differences) if differences else 0.0


def _support_label(scores: list[float], labels: list[str]) -> str:
    """Return the winning tree label, or ``ambiguous`` for a tie."""
    best = float(max(scores))
    if sum(abs(score - best) < _FLOAT_TOL for score in scores) > 1:
        return "ambiguous"
    return labels[scores.index(best)]


def _outlier_loci(gene_delta: list[float]) -> list[bool]:
    """Return a mask for absolute deltas outside Tukey 1.5-IQR whiskers."""
    values = np.abs(np.asarray(gene_delta))
    q1, q3 = np.percentile(values, [25, 75])
    iqr = q3 - q1
    low, high = q1 - 1.5 * iqr, q3 + 1.5 * iqr
    return [bool(value < low or value > high) for value in values]


def _compare_groups(
    group_a: list[str],
    group_b: list[str],
    metrics_csv: Path,
    label_a: str,
    label_b: str,
    output_dir: Path,
) -> tuple[Path, Path, dict[str, Any]]:
    """Write metric means, Mann-Whitney p-values, and comparison boxplots.

    Returns (csv_path, pdf_path, sig_info) where sig_info contains
    ``n_sig_metrics`` and ``sig_metric_names`` for p < 0.05 metrics.
    """
    import math
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    with open(metrics_csv, newline="") as handle:
        reader = csv.DictReader(handle)
        rows = list(reader)
        fieldnames = reader.fieldnames or []
    locus_column = "loci" if "loci" in fieldnames else fieldnames[0]
    requested = set(group_a) | set(group_b)
    missing = requested - {row[locus_column] for row in rows}
    if missing:
        raise ValueError(f"Loci missing from --metrics file: {', '.join(sorted(missing))}")

    numeric_columns = [
        column for column in fieldnames
        if column != locus_column and all(_is_numeric(row.get(column)) for row in rows if row[locus_column] in requested)
    ]
    rows_a = [row for row in rows if row[locus_column] in set(group_a)]
    rows_b = [row for row in rows if row[locus_column] in set(group_b)]
    csv_path = output_dir / f"{label_a}_comparison.csv"
    pdf_path = output_dir / f"{label_a}_comparison.pdf"
    output_dir.mkdir(parents=True, exist_ok=True)
    columns = ["metric", f"{label_a}_mean", f"{label_a}_n", f"{label_b}_mean", f"{label_b}_n", "wilcoxon_p"]
    empty_sig: dict[str, Any] = {"n_sig_metrics": 0, "sig_metric_names": []}
    if not group_a:
        comparison_rows = [
            {
                "metric": column,
                f"{label_a}_mean": "NA",
                f"{label_a}_n": 0,
                f"{label_b}_mean": round(float(np.mean([float(row[column]) for row in rows_b])), 6) if rows_b else "NA",
                f"{label_b}_n": len(rows_b),
                "wilcoxon_p": "NA",
            }
            for column in numeric_columns
        ]
        with open(csv_path, "w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=columns)
            writer.writeheader()
            writer.writerows(comparison_rows)
        figure, axis = plt.subplots(figsize=(4, 3))
        axis.text(0.5, 0.5, f"No {label_a} loci", ha="center", va="center", transform=axis.transAxes)
        figure.savefig(pdf_path, dpi=150, bbox_inches="tight")
        plt.close(figure)
        return csv_path, pdf_path, empty_sig

    comparison_rows = []
    wilcoxon_rows = []
    for column in numeric_columns:
        values_a = [float(row[column]) for row in rows_a]
        values_b = [float(row[column]) for row in rows_b]
        comparison_rows.append({
            "metric": column,
            f"{label_a}_mean": round(float(np.mean(values_a)), 6) if values_a else "NA",
            f"{label_a}_n": len(values_a),
            f"{label_b}_mean": round(float(np.mean(values_b)), 6) if values_b else "NA",
            f"{label_b}_n": len(values_b),
            "wilcoxon_p": round(float(mannwhitneyu(values_a, values_b).pvalue), 6) if values_a and values_b else "NA",
        })
        if values_a and values_b:
            stat = mannwhitneyu(values_a, values_b, alternative="two-sided")
            direction = f"{label_a}_higher" if np.mean(values_a) > np.mean(values_b) else f"{label_b}_higher"
            wilcoxon_rows.append({
                "metric": column,
                "u_statistic": round(float(stat.statistic), 6),
                "p_value": round(float(stat.pvalue), 6),
                "direction": direction,
            })
        else:
            wilcoxon_rows.append({
                "metric": column,
                "u_statistic": "NA",
                "p_value": "NA",
                "direction": "insufficient_data",
            })
    with open(csv_path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        writer.writerows(comparison_rows)

    n_features = len(numeric_columns)
    if n_features == 0:
        figure, axis = plt.subplots(figsize=(4, 3))
        axis.text(0.5, 0.5, "No numeric metrics", ha="center", va="center")
        figure.savefig(pdf_path, dpi=150, bbox_inches="tight")
        plt.close(figure)
        return csv_path, pdf_path, empty_sig

    ncols = min(4, n_features)
    nrows = int(math.ceil(n_features / ncols))
    color_a = "#a6cee3"
    color_b = "#fb9a99"
    fig, axes = plt.subplots(nrows, ncols, figsize=(ncols * 3.5, nrows * 2.8))
    if nrows * ncols == 1:
        axes = np.atleast_2d(axes).reshape(1, 1)
    elif nrows == 1:
        axes = np.array([axes])
    elif ncols == 1:
        axes = np.array([[ax] for ax in axes])

    for idx, column in enumerate(numeric_columns):
        r, c = divmod(idx, ncols)
        ax = axes[r, c] if nrows > 1 or ncols > 1 else axes[0, 0]
        if nrows == 1 or ncols == 1:
            ax = axes[0][idx] if nrows == 1 else axes[idx][0]
        else:
            ax = axes[r, c]
        vals_a = [float(row[column]) for row in rows_a]
        vals_b = [float(row[column]) for row in rows_b]
        bp = ax.boxplot([vals_a, vals_b], tick_labels=[label_a, label_b], patch_artist=True)
        bp["boxes"][0].set_facecolor(color_a)
        bp["boxes"][0].set_alpha(0.6)
        bp["boxes"][1].set_facecolor(color_b)
        bp["boxes"][1].set_alpha(0.6)
        ax.set_title(column, fontsize=9)
        ax.tick_params(axis="x", labelsize=7)
        wrow = next((row for row in wilcoxon_rows if row["metric"] == column), None)
        if wrow and isinstance(wrow.get("p_value"), (int, float)):
            p = wrow["p_value"]
            if p < 0.001:
                sig = "***"
            elif p < 0.01:
                sig = "**"
            elif p < 0.05:
                sig = "*"
            else:
                sig = f"p={p:.3f}"
            all_vals = vals_a + vals_b
            y_top = np.percentile(all_vals, 95) if all_vals else 0
            y_range = max(all_vals) - min(all_vals) if all_vals else 1
            y_pos = y_top + y_range * 0.08
            ax.annotate(sig, xy=(1.5, y_pos), ha="center", fontsize=9, fontweight="bold",
                        va="bottom")

    for idx in range(n_features, nrows * ncols):
        r, c = divmod(idx, ncols)
        axes[r, c].set_visible(False)

    fig.tight_layout(pad=2.0)
    fig.savefig(pdf_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    sig_metrics = [row["metric"] for row in wilcoxon_rows if isinstance(row.get("p_value"), (int, float)) and row["p_value"] < 0.05]
    return csv_path, pdf_path, {"n_sig_metrics": len(sig_metrics), "sig_metric_names": sig_metrics}


def _compare_multiple_groups(
    groups: dict[str, list[str]],
    metrics_csv: Path,
    output_dir: Path,
    prefix: str = "support_comparison",
) -> tuple[Path, Path, dict[str, Any]]:
    """Compare N support groups across metrics — one merged CSV + one merged PDF.

    Returns (csv_path, pdf_path, sig_info) where sig_info maps pair name
    (e.g. ``T1_vs_T2``) to list of significant metric names.
    """
    import math
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    with open(metrics_csv, newline="") as handle:
        reader = csv.DictReader(handle)
        rows = list(reader)
        fieldnames = reader.fieldnames or []
    locus_column = "loci" if "loci" in fieldnames else fieldnames[0]
    all_loci = [loc for locs in groups.values() for loc in locs]
    requested = set(all_loci)
    missing = requested - {row[locus_column] for row in rows}
    if missing:
        raise ValueError(f"Loci missing from --metrics file: {', '.join(sorted(missing))}")

    numeric_columns = [
        column for column in fieldnames
        if column != locus_column and all(_is_numeric(row.get(column)) for row in rows if row[locus_column] in requested)
    ]
    group_order = list(groups.keys())
    group_values: dict[str, dict[str, list[float]]] = {label: {} for label in group_order}
    for label in group_order:
        set_label = set(groups[label])
        for column in numeric_columns:
            group_values[label][column] = [float(r[column]) for r in rows if r[locus_column] in set_label]

    csv_path = output_dir / f"{prefix}.csv"
    pdf_path = output_dir / f"{prefix}.pdf"
    output_dir.mkdir(parents=True, exist_ok=True)

    sig_entries: dict[str, list[str]] = {}
    csv_rows: list[dict[str, Any]] = []
    for column in numeric_columns:
        row: dict[str, Any] = {"metric": column}
        vals_by_label: dict[str, list[float]] = {}
        for label in group_order:
            vals = group_values[label][column]
            vals_by_label[label] = vals
            row[f"{label}_mean"] = round(float(np.mean(vals)), 6) if vals else "NA"
            row[f"{label}_n"] = len(vals)
        for i in range(len(group_order)):
            for j in range(i + 1, len(group_order)):
                la, lb = group_order[i], group_order[j]
                va, vb = vals_by_label[la], vals_by_label[lb]
                pair_key = f"{la}_vs_{lb}_wilcoxon_p"
                if va and vb:
                    p = float(mannwhitneyu(va, vb, alternative="two-sided").pvalue)
                    row[pair_key] = round(p, 6)
                    if p < 0.05:
                        sig_entries.setdefault(f"{la}_vs_{lb}", []).append(column)
                else:
                    row[pair_key] = "NA"
        csv_rows.append(row)

    fieldnames_out = ["metric"]
    for label in group_order:
        fieldnames_out.extend([f"{label}_mean", f"{label}_n"])
    for i in range(len(group_order)):
        for j in range(i + 1, len(group_order)):
            fieldnames_out.append(f"{group_order[i]}_vs_{group_order[j]}_wilcoxon_p")
    with open(csv_path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames_out)
        writer.writeheader()
        writer.writerows(csv_rows)

    n_metrics = len(numeric_columns)
    if n_metrics == 0:
        fig, ax = plt.subplots(figsize=(4, 3))
        ax.text(0.5, 0.5, "No numeric metrics", ha="center", va="center")
        fig.savefig(pdf_path, dpi=150, bbox_inches="tight")
        plt.close(fig)
        return csv_path, pdf_path, sig_entries

    n_groups = len(group_order)
    palette = ["#a6cee3", "#fb9a99", "#b2df8a", "#fdbf6f", "#cab2d6", "#ffff99",
               "#1f78b4", "#e31a1c", "#33a02c", "#ff7f00"][:n_groups]
    ncols = min(4, n_metrics)
    nrows = int(math.ceil(n_metrics / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(ncols * 3.5, nrows * 2.8))
    if nrows * ncols == 1:
        axes = np.atleast_2d(axes).reshape(1, 1)
    elif nrows == 1:
        axes = np.array([axes])
    elif ncols == 1:
        axes = np.array([[ax] for ax in axes])
    pos = list(range(1, n_groups + 1))

    for idx, column in enumerate(numeric_columns):
        r, c = divmod(idx, ncols)
        ax = axes[r, c] if nrows > 1 or ncols > 1 else axes[0, 0]
        data = [group_values[label][column] for label in group_order]
        bp = ax.boxplot(data, positions=pos, patch_artist=True, widths=0.5)
        for patch, color in zip(bp["boxes"], palette):
            patch.set_facecolor(color)
            patch.set_alpha(0.6)
        ax.set_xticks(pos)
        ax.set_xticklabels(group_order, fontsize=8)
        ax.set_title(column, fontsize=9)
        ax.tick_params(axis="x", labelsize=7)
        all_vals = [v for d in data for v in d]
        y_max = max(all_vals) if all_vals else 1
        y_min = min(all_vals) if all_vals else 0
        y_range = y_max - y_min
        y_top = y_max + y_range * 0.05
        bracket_y = y_top
        step = y_range * 0.08
        for i in range(n_groups):
            for j in range(i + 1, n_groups):
                la, lb = group_order[i], group_order[j]
                pair_key = f"{la}_vs_{lb}"
                pair_metrics = sig_entries.get(pair_key, [])
                if column in pair_metrics:
                    entry = next((r for r in csv_rows if r["metric"] == column), None)
                    if entry:
                        p = entry.get(f"{la}_vs_{lb}_wilcoxon_p", "NA")
                        if isinstance(p, (int, float)):
                            if p < 0.001:
                                sig = "***"
                            elif p < 0.01:
                                sig = "**"
                            elif p < 0.05:
                                sig = "*"
                            else:
                                sig = f"p={p:.3f}"
                        else:
                            continue
                    else:
                        continue
                    ax.plot([pos[i], pos[j]], [bracket_y, bracket_y], "k-", lw=0.8)
                    ax.annotate(sig, xy=((pos[i] + pos[j]) / 2, bracket_y),
                                ha="center", va="bottom", fontsize=8, fontweight="bold")
                    bracket_y += step
        ax.margins(y=0.15)

    for idx in range(n_metrics, nrows * ncols):
        r, c = divmod(idx, ncols)
        axes[r, c].set_visible(False)

    fig.tight_layout(pad=2.0)
    fig.savefig(pdf_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return csv_path, pdf_path, sig_entries


def _is_numeric(value: Any) -> bool:
    try:
        float(value)
    except (TypeError, ValueError):
        return False
    return True


def _validate_lnl_inputs(
    *,
    matrix: Path,
    candidate_trees: list[Path],
    model_expr: str | None,
    partitions: Path | None,
    partition_mode: str | None = None,
    locus_ranges: Path | None = None,
    tool_args: str | None = None,
    guide_tree: Path | None = None,
    metrics: Path | None = None,
    threads: str = "auto",
) -> list[str]:
    errors: list[str] = []
    if not matrix.is_file():
        errors.append(f"--matrix does not exist: {matrix}")
    if not candidate_trees:
        errors.append("--candidate-trees must not be empty")
    for index, tree in enumerate(candidate_trees, 1):
        if not tree.is_file():
            errors.append(f"--candidate-trees #{index} does not exist: {tree}")
        elif tree.stat().st_size == 0:
            errors.append(f"--candidate-trees #{index} is empty: {tree}")
    if partitions and locus_ranges:
        errors.append("--partitions and --locus-ranges are mutually exclusive")
    if partition_mode and not partitions:
        errors.append("--partition-mode is only valid when --partitions is provided")
    for option, path in (
        ("--partitions", partitions),
        ("--locus-ranges", locus_ranges),
        ("--metrics", metrics),
        ("--guide-tree", guide_tree),
    ):
        if path and not path.is_file():
            errors.append(f"{option} must be a readable regular file: {path}")
        elif path:
            try:
                path.read_text()
            except OSError as exc:
                errors.append(f"{option} is not readable: {path} ({exc})")
    tokens = set(shlex.split(tool_args)) if tool_args else set()
    for flag in _LNL_BLOCKED_FLAGS & tokens:
        errors.append(f"Blocked flag in --tool-args: {flag}")
    if not (model_expr or partitions or "-m" in tokens or "-p" in tokens or "-Q" in tokens):
        errors.append("Must specify --model-expr, --partitions, or -m/-p/-Q in --tool-args")
    if threads != "auto":
        try:
            n = int(threads)
            if n < 1:
                errors.append(f"--threads must be a positive integer or 'auto', got {threads!r}")
        except ValueError:
            errors.append(f"--threads must be a positive integer or 'auto', got {threads!r}")
    return errors


def _build_lnl_cmd(
    *,
    executable: str,
    matrix: Path,
    candidate_trees: Path,
    prefix: str,
    model_expr: str | None,
    partitions: str | None,
    partition_mode: str | None = None,
    guide_tree: str | None = None,
    threads: str = "auto",
    tool_args: str | None = None,
) -> list[str]:
    cmd = [executable, "-s", str(matrix), "-z", str(candidate_trees)]
    tokens = set(shlex.split(tool_args)) if tool_args else set()
    if "--prefix" not in tokens:
        cmd.extend(["--prefix", prefix])
    if model_expr and "-m" not in tokens:
        cmd.extend(["-m", model_expr])
    if partitions:
        p_flag = f"-{partition_mode or 'p'}"
        if p_flag not in tokens:
            cmd.extend([p_flag, partitions])
    if guide_tree and "-ft" not in tokens:
        cmd.extend(["-ft", guide_tree])
    cmd.append("-wslr")
    if "-T" not in tokens:
        cmd.extend(["-T", str(threads)])
    if tool_args:
        cmd.extend(shlex.split(tool_args))
    return cmd


def _write_csv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str]) -> None:
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def _write_result_json(output_dir: Path, result: dict[str, Any]) -> None:
    (output_dir / "result.json").write_text(json.dumps(result, indent=2))


def _load_result_json(output_dir: Path) -> dict[str, Any]:
    return json.loads((output_dir / "result.json").read_text())


def _plot_support_bar(
    support_values: list[str], tree_labels: list[str], output_path: Path,
    xlabel: str = "Supported topology", ylabel: str = "Count",
) -> Path:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    categories = tree_labels + ["ambiguous"]
    counts = {category: 0 for category in categories}
    for support in support_values:
        counts[support] = counts.get(support, 0) + 1
    figure, axis = plt.subplots(figsize=(max(4, len(categories) * 1.5), 5))
    axis.bar(counts.keys(), counts.values())
    axis.set_xlabel(xlabel)
    axis.set_ylabel(ylabel)
    figure.tight_layout()
    figure.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close(figure)
    return output_path


def run_signal_lnl(
    *,
    matrix: Path,
    candidate_trees: list[Path],
    model_expr: str | None = None,
    partitions: Path | None = None,
    partition_mode: str | None = None,
    locus_ranges: Path | None = None,
    guide_tree: Path | None = None,
    threads: str = "auto",
    iqtree_path: str | None = None,
    tool_args: str | None = None,
    metrics: Path | None = None,
    prefix: str = "lnl",
    resume: bool = False,
    output_dir: Path | None = None,
    overwrite: bool = False,
    dry_run: bool = False,
    quiet: bool = False,
) -> dict[str, Any]:
    run_start = time.time()
    output_dir = (output_dir or Path("runs/posttree/signal/lnl")).resolve()
    matrix = matrix.resolve()
    candidate_trees = [tree.resolve() for tree in candidate_trees]
    errors = _validate_lnl_inputs(
        matrix=matrix, candidate_trees=candidate_trees, model_expr=model_expr,
        partitions=partitions, partition_mode=partition_mode,
        locus_ranges=locus_ranges, tool_args=tool_args,
        guide_tree=guide_tree, metrics=metrics,
        threads=threads,
    )
    params: dict[str, Any] = {
        "matrix": str(matrix),
        "candidate_trees_raw": ",".join(map(str, candidate_trees)),
        "model_expr": model_expr,
        "partitions": str(partitions.resolve()) if partitions else None,
        "partition_mode": partition_mode if partitions else None,
        "locus_ranges": str(locus_ranges.resolve()) if locus_ranges else None,
        "guide_tree": str(guide_tree.resolve()) if guide_tree else None,
        "threads": threads, "iqtree_path": iqtree_path, "tool_args": tool_args,
        "metrics": str(metrics.resolve()) if metrics else None,
        "prefix": prefix, "resume": resume,
        "output_dir": str(output_dir), "overwrite": overwrite,
        "dry_run": dry_run, "quiet": quiet,
    }

    def error_result(message: str, category: str) -> dict[str, Any]:
        return {"status": "error", "command": "", "wall_time": 0.0,
                "tool_versions": {}, "params": params, "key_results": {},
                "error": message, "error_category": category,
                "data": {"cmd": [], "tool_stderr": "", "tool_log": None, "output_files": {}}}

    if errors:
        return error_result("; ".join(errors), "input")
    if overwrite and resume:
        return error_result("--overwrite and --resume are mutually exclusive", "input")
    if not dry_run:
        if overwrite and output_dir.exists():
            shutil.rmtree(output_dir)
        elif not resume and output_dir.exists() and any(output_dir.iterdir()):
            return error_result(
                f"Output directory '{output_dir}' already exists and is non-empty. Use --overwrite to replace it.",
                "input",
            )
        output_dir.mkdir(parents=True, exist_ok=True)

    candidate_trees_path = candidate_trees[0] if len(candidate_trees) == 1 else output_dir / "candidate.trees"
    candidates_merged = len(candidate_trees) > 1
    if candidates_merged and not dry_run:
        candidate_trees_path.write_text("".join(f"{tree.read_text().strip()}\n" for tree in candidate_trees))
    try:
        executable = _resolve_iqtree_path(iqtree_path, dry_run)
    except (ValueError, FileNotFoundError) as exc:
        return error_result(str(exc), "env")
    tool_versions = {"iqtree3": "dry-run"} if dry_run else _detect_iqtree_version(executable)
    cmd = _build_lnl_cmd(
        executable=executable, matrix=matrix, candidate_trees=candidate_trees_path,
        prefix=prefix, model_expr=model_expr,
        partitions=str(partitions.resolve()) if partitions else None,
        partition_mode=partition_mode,
        guide_tree=str(guide_tree.resolve()) if guide_tree else None,
        threads=threads, tool_args=tool_args,
    )
    cli_parts = ["phyloai", "posttree", "signal", "lnl", "--matrix", str(matrix),
                 "--candidate-trees", params["candidate_trees_raw"]]
    if model_expr:
        cli_parts.extend(["--model-expr", model_expr])
    if partitions:
        cli_parts.extend(["--partitions", str(partitions)])
        cli_parts.extend(["--partition-mode", partition_mode or "p"])
    if locus_ranges:
        cli_parts.extend(["--locus-ranges", str(locus_ranges)])
    if guide_tree:
        cli_parts.extend(["--guide-tree", str(guide_tree)])
    if tool_args:
        cli_parts.extend(["--tool-args", tool_args])
    if metrics:
        cli_parts.extend(["--metrics", str(metrics)])
    if prefix != "lnl":
        cli_parts.extend(["--prefix", prefix])
    cli_parts.extend(["--threads", str(threads), "-o", str(output_dir)])
    if iqtree_path:
        cli_parts.extend(["--iqtree-path", iqtree_path])
    if resume:
        cli_parts.append("--resume")
    if overwrite:
        cli_parts.append("--overwrite")
    if dry_run:
        cli_parts.append("--dry-run")
    if quiet:
        cli_parts.append("-q")
    command = shlex.join(cli_parts)
    if dry_run:
        return {"status": "success", "command": command, "wall_time": 0.0,
                "tool_versions": tool_versions, "params": params, "key_results": {},
                "error": None, "data": {"cmd": cmd, "tool_stderr": "", "tool_log": None, "output_files": {}}}

    iqtree_dir = output_dir / "iqtree"
    iqtree_dir.mkdir(parents=True, exist_ok=True)
    proc = subprocess.run(cmd, stdout=None, stderr=subprocess.PIPE, text=True, cwd=str(iqtree_dir))
    tool_stderr = proc.stderr.strip() if proc.stderr else ""
    output_files: dict[str, Any] = {
        "iqtree_report": {"path": str(iqtree_dir / f"{prefix}.iqtree"), "description": "IQ-TREE native report"},
        "iqtree_sitelh": {"path": str(iqtree_dir / f"{prefix}.sitelh"), "description": "IQ-TREE raw site log-likelihoods"},
        "iqtree_log": {"path": str(iqtree_dir / f"{prefix}.log"), "description": "IQ-TREE console log"},
    }
    if candidates_merged:
        output_files["candidate_trees"] = {"path": str(candidate_trees_path), "description": "Merged candidate trees (IQ-TREE -z input)"}
    if proc.returncode:
        result = error_result(f"IQ-TREE exited with code {proc.returncode}", "tool")
        result.update({"command": command, "wall_time": time.time() - run_start, "tool_versions": tool_versions})
        result["data"] = {"cmd": cmd, "tool_stderr": tool_stderr, "tool_log": str(iqtree_dir / f"{prefix}.log"), "output_files": output_files}
        _write_result_json(output_dir, result)
        return result

    try:
        tree_labels, site_scores = _parse_sitelh(iqtree_dir / f"{prefix}.sitelh")
        site_rows = []
        for site, scores in enumerate(zip(*site_scores), 1):
            row = {"site": site, **{f"lnL_{label}": round(score, 6) for label, score in zip(tree_labels, scores)},
                   "ΔSLS": round(_delta_score(list(scores)), 6), "support": _support_label(list(scores), tree_labels)}
            site_rows.append(row)
        site_rows.sort(key=lambda row: row["ΔSLS"], reverse=True)
        site_csv = output_dir / "site_lnl.csv"
        _write_csv(site_csv, site_rows, ["site"] + [f"lnL_{label}" for label in tree_labels] + ["ΔSLS", "support"])
        output_files["site_lnl"] = {"path": str(site_csv), "description": "Site-wise lnL scores per tree, ΔSLS, support; sorted by ΔSLS descending"}
        site_plot = _plot_support_bar([row["support"] for row in site_rows], tree_labels, output_dir / "site_support.pdf")
        output_files["site_support_plot"] = {"path": str(site_plot), "description": "Site support distribution bar chart"}
        site_support_counts = {label: sum(1 for r in site_rows if r["support"] == label) for label in tree_labels}
        site_support_counts["ambiguous"] = sum(1 for r in site_rows if r["support"] == "ambiguous")
        support_summary_rows = [{"tree": label, "n_sites": site_support_counts.get(label, 0)} for label in tree_labels] + [{"tree": "ambiguous", "n_sites": site_support_counts["ambiguous"]}]
        support_summary_csv = output_dir / "support_summary_sites.csv"
        _write_csv(support_summary_csv, support_summary_rows, ["tree", "n_sites"])
        output_files["support_summary_sites"] = {"path": str(support_summary_csv), "description": "Number of sites supporting each topology"}
        key_results: dict[str, Any] = {"n_trees": len(tree_labels), "n_sites": len(site_rows),
                                        "site_support_counts": {label: count for label, count in site_support_counts.items()}}

        boundary_path = partitions or locus_ranges
        if boundary_path:
            records = _parse_partition_ranges(boundary_path.resolve())
            for record in records:
                if record["end"] > len(site_rows):
                    raise ValueError(f"Locus range {record['locus']} ends at {record['end']}, beyond sitelh site count {len(site_rows)}")
            gene_rows = []
            for record in records:
                scores = _sum_gene_lnl(site_scores, record["start"], record["end"])
                row = {"locus": record["locus"], **{f"lnL_{label}": round(score, 6) for label, score in zip(tree_labels, scores)},
                       "ΔGLS": round(_delta_score(scores), 6), "support": _support_label(scores, tree_labels)}
                if len(tree_labels) == 2:
                    row["support_sig"] = abs(row["ΔGLS"]) >= 2.0
                gene_rows.append(row)
            gene_rows.sort(key=lambda row: row["ΔGLS"], reverse=True)
            fields = ["locus"] + [f"lnL_{label}" for label in tree_labels] + ["ΔGLS", "support"]
            if len(tree_labels) == 2:
                fields.append("support_sig")
            gene_csv = output_dir / "gene_lnl.csv"
            _write_csv(gene_csv, gene_rows, fields)
            output_files["gene_lnl"] = {"path": str(gene_csv), "description": "Gene-wise lnL scores per tree, ΔGLS, support; sorted by ΔGLS descending"}
            gene_plot = _plot_support_bar([row["support"] for row in gene_rows], tree_labels, output_dir / "gene_support.pdf", ylabel="Number of genes")
            output_files["gene_support_plot"] = {"path": str(gene_plot), "description": "Gene support distribution bar chart"}
            gene_support_counts = {label: sum(1 for r in gene_rows if r["support"] == label) for label in tree_labels}
            gene_support_counts["ambiguous"] = sum(1 for r in gene_rows if r["support"] == "ambiguous")
            gene_summary_rows = [{"tree": label, "n_genes": gene_support_counts.get(label, 0)} for label in tree_labels] + [{"tree": "ambiguous", "n_genes": gene_support_counts["ambiguous"]}]
            gene_summary_csv = output_dir / "support_summary_genes.csv"
            _write_csv(gene_summary_csv, gene_summary_rows, ["tree", "n_genes"])
            output_files["support_summary_genes"] = {"path": str(gene_summary_csv), "description": "Number of genes supporting each topology"}
            outlier_loci = [row["locus"] for row, flag in zip(gene_rows, _outlier_loci([row["ΔGLS"] for row in gene_rows])) if flag]
            outlier_path = output_dir / "outlier_genes.txt"
            outlier_path.write_text("\n".join(outlier_loci) + ("\n" if outlier_loci else ""))
            output_files["outlier_genes"] = {"path": str(outlier_path), "description": "Loci with |ΔGLS| outside boxplot whiskers (Shen 2017 eq. 3/4)"}
            key_results.update(n_loci=len(gene_rows), n_outlier_genes=len(outlier_loci),
                               gene_support_counts={label: count for label, count in gene_support_counts.items()})
            if metrics:
                others = [row["locus"] for row in gene_rows if row["locus"] not in outlier_loci]
                comparison_csv, comparison_plot, sig_info = _compare_groups(outlier_loci, others, metrics.resolve(), "outlier", "non_outlier", output_dir)
                output_files["outlier_comparison"] = {"path": str(comparison_csv), "description": "Outlier vs non-outlier per-metric means and Wilcoxon p-values"}
                output_files["outlier_comparison_plot"] = {"path": str(comparison_plot), "description": "Outlier vs non-outlier metric distribution boxplots"}
                key_results.update(n_sig_metrics_outlier=sig_info["n_sig_metrics"],
                                   sig_metric_names_outlier=sig_info["sig_metric_names"])
                support_groups: dict[str, list[str]] = {}
                for row in gene_rows:
                    if row["support"] != "ambiguous":
                        support_groups.setdefault(row["support"], []).append(row["locus"])
                if len(support_groups) >= 2:
                    support_labels = [lbl for lbl in tree_labels if lbl in support_groups]
                    ordered_groups = {lbl: support_groups[lbl] for lbl in support_labels}
                    comp_csv, comp_plot, pair_sigs = _compare_multiple_groups(
                        ordered_groups, metrics.resolve(), output_dir, prefix="support_comparison",
                    )
                    output_files["support_comparison"] = {
                        "path": str(comp_csv),
                        "description": "Per-metric means and pairwise Wilcoxon p-values across tree support groups",
                    }
                    output_files["support_comparison_plot"] = {
                        "path": str(comp_plot),
                        "description": "Multi-group metric distribution boxplots with pairwise significance",
                    }
                    if pair_sigs:
                        key_results["support_comparison_sig_metrics"] = pair_sigs
    except (OSError, ValueError, IndexError) as exc:
        result = error_result(str(exc), "output")
        result.update({"command": command, "wall_time": time.time() - run_start, "tool_versions": tool_versions})
        result["data"] = {"cmd": cmd, "tool_stderr": tool_stderr, "tool_log": str(iqtree_dir / f"{prefix}.log"), "output_files": output_files}
        _write_result_json(output_dir, result)
        return result

    result = {"status": "success", "command": command, "wall_time": time.time() - run_start,
              "tool_versions": tool_versions, "params": params, "key_results": key_results, "error": None,
               "data": {"cmd": cmd, "tool_stderr": "", "tool_log": str(iqtree_dir / f"{prefix}.log"),
                       "summary": key_results, "output_files": output_files}}
    _write_result_json(output_dir, result)
    return result


# ---------------------------------------------------------------------------
#  signal consistent — GLS + GQS consistency analysis
# ---------------------------------------------------------------------------


def _prune_reference_tree(ref_tree_str: str, taxa_to_remove: set[str]) -> str:
    """Prune taxa from a newick tree string using Bio.Phylo. Returns pruned newick."""
    from io import StringIO

    from Bio import Phylo

    tree = Phylo.read(StringIO(ref_tree_str), "newick")
    for taxon in taxa_to_remove:
        try:
            tree.prune(taxon)
        except Exception:
            pass
    out = StringIO()
    Phylo.write(tree, out, "newick")
    return out.getvalue().strip()


def _count_tree_taxa(newick_str: str) -> int:
    from io import StringIO

    from Bio import Phylo

    tree = Phylo.read(StringIO(newick_str), "newick")
    return len(tree.get_terminals())


def _get_tree_taxa(newick_str: str) -> set[str]:
    from io import StringIO

    from Bio import Phylo

    tree = Phylo.read(StringIO(newick_str), "newick")
    return {c.name for c in tree.get_terminals() if c.name}


def _run_wastral_gqs(
    gene_tree_path: Path,
    ref_tree_str: str,
    wastral_exe: str,
    work_dir: Path,
    locus: str,
) -> float:
    """Run wastral quartet score for one gene tree vs one reference tree.

    Returns score float. Raises RuntimeError if wastral exits non-zero
    or Score: line is absent (external-tool failure -> caller writes error result).
    """
    ref_path = work_dir / f"ref_{locus}.nwk"
    ref_path.write_text(ref_tree_str)
    proc = subprocess.run(
        [wastral_exe, "-i", str(gene_tree_path), "-C", "-c", str(ref_path), "--mode", "4"],
        capture_output=True, text=True,
    )
    combined = (proc.stdout or "") + "\n" + (proc.stderr or "")
    if proc.returncode != 0:
        raise RuntimeError(
            f"wastral exited {proc.returncode} for locus {locus!r}.\n{combined[:500]}"
        )
    for line in combined.splitlines():
        if line.strip().startswith("Score:"):
            try:
                return float(line.split(":")[1].strip())
            except ValueError:
                pass
    raise RuntimeError(
        f"wastral produced no 'Score:' line for locus {locus!r}.\n{combined[:500]}"
    )


def _compute_gqs_for_locus(
    gene_tree_path: Path,
    logical_locus: str,
    t1_str: str,
    t2_str: str,
    ref_taxa: set[str],
    wastral_exe: str,
    work_dir: Path,
    tree_labels: list[str],
) -> dict[str, Any]:
    """Compute GQS for a single gene tree. Returns gqs row dict.

    Raises ValueError for input errors (extra taxa -> hard error).
    Raises RuntimeError for wastral tool failures -> caller writes error result.
    """
    gene_str = gene_tree_path.read_text()
    gene_taxa = _get_tree_taxa(gene_str)

    extra = gene_taxa - ref_taxa
    if extra:
        raise ValueError(
            f"Gene tree {gene_tree_path.name} contains taxa not in reference trees: "
            f"{', '.join(sorted(extra))}. Matrix and gene trees must share the same taxon set."
        )

    missing = ref_taxa - gene_taxa
    t1_pruned = _prune_reference_tree(t1_str, missing)
    t2_pruned = _prune_reference_tree(t2_str, missing)

    if _count_tree_taxa(t1_pruned) < 4:
        return {
            "locus": logical_locus,
            "GQS_T1": None,
            "GQS_T2": None,
            "\u0394GQS": None,
            "support": "ambiguous",
            "status": "skipped",
            "reason": "pruned_tree_too_small",
        }

    locus_work = work_dir / logical_locus
    locus_work.mkdir(parents=True, exist_ok=True)
    gqs_t1 = _run_wastral_gqs(gene_tree_path, t1_pruned, wastral_exe, locus_work, logical_locus)
    gqs_t2 = _run_wastral_gqs(gene_tree_path, t2_pruned, wastral_exe, locus_work, logical_locus)

    delta = gqs_t1 - gqs_t2
    if abs(delta) < _FLOAT_TOL:
        support = "ambiguous"
    elif delta > 0:
        support = tree_labels[0]
    else:
        support = tree_labels[1]

    return {
        "locus": logical_locus,
        "GQS_T1": round(gqs_t1, 6),
        "GQS_T2": round(gqs_t2, 6),
        "\u0394GQS": round(delta, 6),
        "support": support,
        "status": "success",
        "reason": None,
    }


def _validate_consistent_inputs(
    *,
    matrix: Path,
    candidate_trees: list[Path],
    tree_dir: Path,
    model_expr: str | None,
    partitions: Path | None,
    partition_mode: str | None = None,
    locus_ranges: Path | None = None,
    tool_args: str | None = None,
    guide_tree: Path | None = None,
    threads: str = "auto",
    metrics: Path | None = None,
) -> tuple[list[str], list[str]]:
    """Validate consistent inputs and normalise candidate trees.

    Returns (errors, normalised_trees) where normalised_trees is always
    exactly 2 newick strings (T1, T2) on success, empty on error.
    """
    errors: list[str] = []
    norm_trees: list[str] = []
    if not matrix.exists() or not matrix.is_file():
        errors.append(f"--matrix does not exist or is not a regular file: {matrix}")
    if partitions and (not partitions.exists() or not partitions.is_file()):
        errors.append(f"--partitions does not exist or is not a regular file: {partitions}")
    if locus_ranges and (not locus_ranges.exists() or not locus_ranges.is_file()):
        errors.append(f"--locus-ranges does not exist or is not a regular file: {locus_ranges}")
    if guide_tree and (not guide_tree.exists() or not guide_tree.is_file()):
        errors.append(f"--guide-tree does not exist or is not a regular file: {guide_tree}")
    if metrics and (not metrics.exists() or not metrics.is_file()):
        errors.append(f"--metrics does not exist or is not a regular file: {metrics}")
    # Normalise candidate trees before counting: a single tree-list file
    # is expanded in-place.
    if len(candidate_trees) == 1:
        try:
            content = candidate_trees[0].read_text().strip()
        except OSError:
            errors.append(f"--candidate-trees #{1} is not readable: {candidate_trees[0]}")
        else:
            norm_trees = [t.strip() for t in content.splitlines() if t.strip()]
    elif len(candidate_trees) == 2:
        for i, ct in enumerate(candidate_trees):
            try:
                norm_trees.append(ct.read_text().strip())
            except OSError:
                errors.append(f"--candidate-trees #{i + 1} is not readable: {ct}")
    if len(norm_trees) != 2:
        if not errors:
            errors.append(
                f"--candidate-trees must contain exactly 2 trees, "
                f"got {len(norm_trees) if norm_trees else len(candidate_trees)}"
            )
    for i, ct in enumerate(candidate_trees):
        if not ct.exists() or not ct.is_file():
            errors.append(f"--candidate-trees #{i + 1} does not exist: {ct}")
        elif ct.stat().st_size == 0:
            errors.append(f"--candidate-trees #{i + 1} is empty: {ct}")
    if not tree_dir.exists() or not tree_dir.is_dir():
        errors.append(f"--tree-dir does not exist: {tree_dir}")
    if partitions and locus_ranges:
        errors.append("--partitions and --locus-ranges are mutually exclusive")
    if partition_mode and not partitions:
        errors.append("--partition-mode is only valid when --partitions is provided")
    if not partitions and not locus_ranges:
        errors.append("Must provide --partitions or --locus-ranges (GLS requires locus boundaries)")
    has_tool_model = False
    if tool_args:
        toks = set(shlex.split(tool_args))
        has_tool_model = "-m" in toks or "-p" in toks or "-Q" in toks
        for flag in _CONSISTENT_BLOCKED_FLAGS:
            if flag in toks:
                errors.append(f"Blocked flag in --tool-args: {flag}")
    if not model_expr and not partitions and not has_tool_model:
        errors.append("Must specify --model-expr, --partitions, or -m/-p/-Q in --tool-args")
    if threads != "auto":
        try:
            n = int(threads)
            if n < 1:
                errors.append(f"--threads must be a positive integer or 'auto', got {threads!r}")
        except ValueError:
            errors.append(f"--threads must be a positive integer or 'auto', got {threads!r}")
    return errors, norm_trees


def run_signal_consistent(
    *,
    matrix: Path,
    candidate_trees: list[Path],
    tree_dir: Path,
    model_expr: str | None = None,
    partitions: Path | None = None,
    partition_mode: str | None = None,
    locus_ranges: Path | None = None,
    guide_tree: Path | None = None,
    threads: str = "auto",
    iqtree_path: str | None = None,
    wastral_path: str | None = None,
    tool_args: str | None = None,
    metrics: Path | None = None,
    prefix: str = "consistent",
    resume: bool = False,
    output_dir: Path | None = None,
    overwrite: bool = False,
    dry_run: bool = False,
    quiet: bool = False,
) -> dict[str, Any]:
    run_start = time.time()
    if output_dir is None:
        output_dir = Path("runs/posttree/signal/consistent")
    output_dir = output_dir.resolve()
    matrix = matrix.resolve()
    candidate_trees = [ct.resolve() for ct in candidate_trees]
    tree_dir = tree_dir.resolve()

    errors, norm_trees = _validate_consistent_inputs(
        matrix=matrix,
        candidate_trees=candidate_trees,
        tree_dir=tree_dir,
        model_expr=model_expr,
        partitions=partitions,
        partition_mode=partition_mode,
        locus_ranges=locus_ranges,
        tool_args=tool_args,
        guide_tree=guide_tree,
        threads=threads,
        metrics=metrics,
    )

    params: dict[str, Any] = {
        "matrix": str(matrix),
        "candidate_trees_raw": ",".join(str(ct) for ct in candidate_trees),
        "tree_dir": str(tree_dir),
        "model_expr": model_expr,
        "partitions": str(partitions.resolve()) if partitions else None,
        "partition_mode": partition_mode if partitions else None,
        "locus_ranges": str(locus_ranges.resolve()) if locus_ranges else None,
        "guide_tree": str(guide_tree.resolve()) if guide_tree else None,
        "threads": threads,
        "iqtree_path": iqtree_path,
        "wastral_path": wastral_path,
        "tool_args": tool_args,
        "metrics": str(metrics.resolve()) if metrics else None,
        "prefix": prefix,
        "resume": resume,
        "output_dir": str(output_dir),
        "overwrite": overwrite,
        "dry_run": dry_run,
        "quiet": quiet,
    }

    if errors:
        return {
            "status": "error",
            "command": "",
            "wall_time": 0.0,
            "tool_versions": {},
            "params": params,
            "key_results": {},
            "error": "; ".join(errors),
            "error_category": "input",
            "data": {"cmd": [], "tool_stderr": "", "tool_log": None, "output_files": {}},
        }
    if overwrite and resume:
        return {
            "status": "error",
            "command": "",
            "wall_time": 0.0,
            "tool_versions": {},
            "params": params,
            "key_results": {},
            "error": "--overwrite and --resume are mutually exclusive",
            "error_category": "input",
            "data": {"cmd": [], "tool_stderr": "", "tool_log": None, "output_files": {}},
        }

    if not dry_run:
        if overwrite and output_dir.exists():
            shutil.rmtree(output_dir)
        elif not resume and output_dir.exists() and any(output_dir.iterdir()):
            return {
                "status": "error",
                "command": "",
                "wall_time": 0.0,
                "tool_versions": {},
                "params": params,
                "key_results": {},
                "error": f"Output directory '{output_dir}' already exists and is non-empty. Use --overwrite to replace it.",
                "error_category": "input",
                "data": {"cmd": [], "tool_stderr": "", "tool_log": None, "output_files": {}},
            }
        output_dir.mkdir(parents=True, exist_ok=True)

    if len(candidate_trees) == 1:
        candidate_trees_path = candidate_trees[0]
    else:
        candidate_trees_path = output_dir / "candidate.trees"
        if not dry_run:
            with open(candidate_trees_path, "w") as fh:
                for ct in candidate_trees:
                    fh.write(ct.read_text().strip() + "\n")

    iqtree_flag = f"-{partition_mode}" if partitions and partition_mode else "-m"
    resolved_partitions = str(partitions.resolve()) if partitions else None
    resolved_guide_tree = str(guide_tree.resolve()) if guide_tree else None

    try:
        iqtree_exe = _resolve_iqtree_path(iqtree_path, dry_run)
    except (ValueError, FileNotFoundError) as e:
        return {
            "status": "error",
            "command": "",
            "wall_time": 0.0,
            "tool_versions": {},
            "params": params,
            "key_results": {},
            "error": str(e),
            "error_category": "env",
            "data": {"cmd": [], "tool_stderr": "", "tool_log": None, "output_files": {}},
        }

    from phyloai.tree.msc import _resolve_wastral_path, _detect_wastral_version

    try:
        wastral_exe = _resolve_wastral_path(wastral_path, dry_run)
    except (ValueError, FileNotFoundError) as e:
        return {
            "status": "error",
            "command": "",
            "wall_time": 0.0,
            "tool_versions": {},
            "params": params,
            "key_results": {},
            "error": str(e),
            "error_category": "env",
            "data": {"cmd": [], "tool_stderr": "", "tool_log": None, "output_files": {}},
        }

    if not dry_run:
        tool_versions = {**_detect_iqtree_version(iqtree_exe), **_detect_wastral_version(wastral_exe)}
    else:
        tool_versions = {"iqtree3": "dry-run", "wastral": "dry-run"}

    cmd: list[str] = [iqtree_exe, "-s", str(matrix), "-z", str(candidate_trees_path)]
    tool_toks = set(shlex.split(tool_args)) if tool_args else set()
    if "--prefix" not in tool_toks:
        cmd.extend(["--prefix", prefix])
    if model_expr and "-m" not in tool_toks:
        cmd.extend(["-m", model_expr])
    if partitions:
        p_flag = f"-{partition_mode or 'p'}"
        if p_flag not in tool_toks:
            cmd.extend([p_flag, resolved_partitions])
    cmd.append("-wslr")
    if resolved_guide_tree and "-ft" not in tool_toks:
        cmd.extend(["-ft", resolved_guide_tree])
    if "-T" not in tool_toks:
        cmd.extend(["-T", str(threads)])
    if tool_args:
        cmd.extend(shlex.split(tool_args))

    cli_parts = [
        "phyloai", "posttree", "signal", "consistent",
        "--matrix", str(matrix),
        "--candidate-trees", params["candidate_trees_raw"],
        "--tree-dir", str(tree_dir),
    ]
    if model_expr:
        cli_parts.extend(["--model-expr", model_expr])
    if partitions:
        cli_parts.extend(["--partitions", str(partitions), "--partition-mode", partition_mode or "p"])
    if locus_ranges:
        cli_parts.extend(["--locus-ranges", str(locus_ranges)])
    if guide_tree:
        cli_parts.extend(["--guide-tree", str(guide_tree)])
    if tool_args:
        cli_parts.extend(["--tool-args", tool_args])
    if metrics:
        cli_parts.extend(["--metrics", str(metrics)])
    if prefix != "consistent":
        cli_parts.extend(["--prefix", prefix])
    cli_parts.extend(["--threads", str(threads), "-o", str(output_dir)])
    if iqtree_path:
        cli_parts.extend(["--iqtree-path", iqtree_path])
    if wastral_path:
        cli_parts.extend(["--wastral-path", wastral_path])
    if resume:
        cli_parts.append("--resume")
    if overwrite:
        cli_parts.append("--overwrite")
    if dry_run:
        cli_parts.append("--dry-run")
    if quiet:
        cli_parts.append("-q")
    full_command = shlex.join(cli_parts)

    t1_str, t2_str = norm_trees

    try:
        t1_taxa = _get_tree_taxa(t1_str)
        t2_taxa = _get_tree_taxa(t2_str)
    except Exception as exc:
        result = {
            "status": "error",
            "command": full_command,
            "wall_time": time.time() - run_start,
            "tool_versions": tool_versions,
            "params": params,
            "key_results": {},
            "error": f"Failed to parse candidate tree Newick: {exc}",
            "error_category": "input",
            "data": {"cmd": cmd, "tool_stderr": "", "tool_log": None, "output_files": {}},
        }
        if not dry_run:
            _write_result_json(output_dir, result)
        return result
    if t1_taxa != t2_taxa:
        diff = t1_taxa.symmetric_difference(t2_taxa)
        result = {
            "status": "error",
            "command": full_command,
            "wall_time": time.time() - run_start,
            "tool_versions": tool_versions,
            "params": params,
            "key_results": {},
            "error": f"Candidate trees T1 and T2 have different taxon sets: {', '.join(sorted(diff))}",
            "error_category": "input",
            "data": {"cmd": cmd, "tool_stderr": "", "tool_log": None, "output_files": {}},
        }
        if not dry_run:
            _write_result_json(output_dir, result)
        return result

    # Validate partition ranges and locus<->gene-tree match BEFORE IQ-TREE.
    boundary_path = (partitions or locus_ranges).resolve()
    try:
        partition_recs = _parse_partition_ranges(boundary_path)
    except ValueError as exc:
        result = {
            "status": "error",
            "command": full_command,
            "wall_time": time.time() - run_start,
            "tool_versions": tool_versions,
            "params": params,
            "key_results": {},
            "error": f"Failed to parse partition file: {exc}",
            "error_category": "input",
            "data": {"cmd": cmd, "tool_stderr": "", "tool_log": None, "output_files": {}},
        }
        if not dry_run:
            _write_result_json(output_dir, result)
        return result
    partition_loci = {rec["locus"]: rec for rec in partition_recs}

    try:
        gene_tree_map = scan_tree_dir(tree_dir)
    except ValueError as exc:
        result = {
            "status": "error",
            "command": full_command,
            "wall_time": time.time() - run_start,
            "tool_versions": tool_versions,
            "params": params,
            "key_results": {},
            "error": f"--tree-dir contains ambiguous or duplicate filenames: {exc}",
            "error_category": "input",
            "data": {"cmd": cmd, "tool_stderr": "", "tool_log": None, "output_files": {}},
        }
        _write_result_json(output_dir, result) if not dry_run else None
        return result

    tree_loci = set(gene_tree_map.keys())
    # ponytail: expand with two-suffix reductions for files like gene.fa.treefile → gene
    for tree_path in list(gene_tree_map.values()):
        _, cand2 = logical_tree_locus_candidates(tree_path)
        if cand2 and cand2 not in gene_tree_map:
            tree_loci.add(cand2)
            gene_tree_map[cand2] = tree_path
    partition_loci_set = set(partition_loci.keys())
    missing_trees = partition_loci_set - tree_loci
    if missing_trees:
        result = {
            "status": "error",
            "command": full_command,
            "wall_time": time.time() - run_start,
            "tool_versions": tool_versions,
            "params": params,
            "key_results": {},
            "error": f"Loci in partition file with no gene tree: {', '.join(sorted(missing_trees))}",
            "error_category": "input",
            "data": {"cmd": cmd, "tool_stderr": "", "tool_log": None, "output_files": {}},
        }
        if not dry_run:
            _write_result_json(output_dir, result)
        return result

    if dry_run:
        return {
            "status": "success",
            "command": full_command,
            "wall_time": 0.0,
            "tool_versions": tool_versions,
            "params": params,
            "key_results": {},
            "error": None,
            "data": {"cmd": cmd, "tool_stderr": "", "tool_log": None, "output_files": {}},
        }

    iqtree_dir = output_dir / "iqtree"
    iqtree_dir.mkdir(parents=True, exist_ok=True)
    proc = subprocess.run(
        cmd, cwd=str(iqtree_dir),
        stdout=None,
        stderr=subprocess.PIPE,
        text=True,
    )
    tool_stderr = proc.stderr.strip() if proc.stderr else ""
    tool_log_path = iqtree_dir / f"{prefix}.log"

    if proc.returncode != 0:
        result = {
            "status": "error",
            "command": full_command,
            "wall_time": time.time() - run_start,
            "tool_versions": tool_versions,
            "params": params,
            "key_results": {},
            "error": f"IQ-TREE exited {proc.returncode}",
            "error_category": "tool",
            "data": {"cmd": cmd, "tool_stderr": tool_stderr, "tool_log": str(tool_log_path), "output_files": {}},
        }
        _write_result_json(output_dir, result)
        return result

    sitelh_path = iqtree_dir / f"{prefix}.sitelh"
    output_files: dict[str, Any] = {
        "iqtree_report": {"path": str(iqtree_dir / f"{prefix}.iqtree"), "description": "IQ-TREE native report"},
        "iqtree_sitelh": {"path": str(sitelh_path), "description": "IQ-TREE raw site log-likelihoods"},
        "iqtree_log": {"path": str(tool_log_path), "description": "IQ-TREE console log"},
    }

    try:
        tree_labels, site_scores = _parse_sitelh(sitelh_path)

        # Validate partition ranges do not exceed .sitelh site count.
        n_sites = len(site_scores[0]) if site_scores else 0
        for rec in partition_recs:
            if rec["end"] > n_sites:
                raise ValueError(
                    f"Partition range {rec['start']}-{rec['end']} for locus "
                    f"{rec['locus']!r} exceeds site count ({n_sites}) in .sitelh"
                )

        gls_rows = []
        for rec in partition_recs:
            gene_scores = _sum_gene_lnl(site_scores, rec["start"], rec["end"])
            delta = _delta_score(gene_scores)
            support = _support_label(gene_scores, tree_labels)
            gls_rows.append({
                "locus": rec["locus"],
                "lnL_T1": round(gene_scores[0], 6),
                "lnL_T2": round(gene_scores[1], 6),
                "\u0394GLS": round(delta, 6),
                "support": support,
                "support_sig": abs(delta) >= 2.0,
            })

        gls_csv_path = output_dir / "gls.csv"
        _write_csv(gls_csv_path, gls_rows,
                   ["locus", "lnL_T1", "lnL_T2", "\u0394GLS", "support", "support_sig"])
    except (OSError, ValueError, IndexError) as exc:
        result = {
            "status": "error",
            "command": full_command,
            "wall_time": time.time() - run_start,
            "tool_versions": tool_versions,
            "params": params,
            "key_results": {},
            "error": str(exc),
            "error_category": "output",
            "data": {"cmd": cmd, "tool_stderr": tool_stderr,
                     "tool_log": str(tool_log_path), "output_files": output_files},
        }
        _write_result_json(output_dir, result)
        return result

    ref_taxa = t1_taxa
    n_workers = os.cpu_count() if threads == "auto" else int(threads)
    work_dir = output_dir / "_gqs_work"
    work_dir.mkdir(parents=True, exist_ok=True)

    try:
        gqs_rows: list[dict[str, Any]] = []
        with ProcessPoolExecutor(max_workers=n_workers) as executor:
            futures = {
                executor.submit(
                    _compute_gqs_for_locus,
                    gene_tree_map[locus],
                    locus,
                    t1_str,
                    t2_str,
                    ref_taxa,
                    wastral_exe,
                    work_dir,
                    tree_labels,
                ): locus
                for locus in partition_loci_set
            }
            locus_to_gqs: dict[str, dict] = {}
            for future in as_completed(futures):
                locus = futures[future]
                try:
                    locus_to_gqs[locus] = future.result()
                except ValueError as exc:
                    result = {
                        "status": "error",
                        "command": full_command,
                        "wall_time": time.time() - run_start,
                        "tool_versions": tool_versions,
                        "params": params,
                        "key_results": {},
                        "error": str(exc),
                        "error_category": "input",
                        "data": {"cmd": cmd, "tool_stderr": "", "tool_log": None, "output_files": {}},
                    }
                    _write_result_json(output_dir, result)
                    return result
                except RuntimeError as exc:
                    result = {
                        "status": "error",
                        "command": full_command,
                        "wall_time": time.time() - run_start,
                        "tool_versions": tool_versions,
                        "params": params,
                        "key_results": {},
                        "error": str(exc),
                        "error_category": "tool",
                        "data": {"cmd": cmd, "tool_stderr": str(exc), "tool_log": None, "output_files": {}},
                    }
                    _write_result_json(output_dir, result)
                    return result

        for rec in partition_recs:
            gqs_rows.append(locus_to_gqs[rec["locus"]])

        gqs_csv_path = output_dir / "gqs.csv"
        _write_csv(gqs_csv_path, gqs_rows,
                   ["locus", "GQS_T1", "GQS_T2", "\u0394GQS", "support", "status", "reason"])

        gls_map = {r["locus"]: r["support"] for r in gls_rows}
        gqs_map = {r["locus"]: r for r in gqs_rows}
        consistent_loci, inconsistent_loci = [], []
        for locus in partition_loci_set:
            gls_sup = gls_map.get(locus, "ambiguous")
            gqs_row = gqs_map.get(locus, {})
            gqs_sup = gqs_row.get("support", "ambiguous")
            gqs_status = gqs_row.get("status", "skipped")
            if gls_sup == gqs_sup and gls_sup != "ambiguous" and gqs_status == "success":
                consistent_loci.append(locus)
            else:
                inconsistent_loci.append(locus)

        (output_dir / "consistent_genes.txt").write_text("\n".join(sorted(consistent_loci)) + "\n")
        (output_dir / "inconsistent_genes.txt").write_text("\n".join(sorted(inconsistent_loci)) + "\n")

        gls_pdf = _plot_support_bar(
            [r["support"] for r in gls_rows],
            tree_labels,
            output_dir / "gls_support.pdf",
            ylabel="Number of genes",
        )
        gqs_pdf = _plot_support_bar(
            [r["support"] for r in gqs_rows if r["status"] == "success"],
            tree_labels,
            output_dir / "gqs_support.pdf",
            ylabel="Number of genes",
        )

        n_gqs_skipped = sum(1 for r in gqs_rows if r["status"] == "skipped")
        output_files.update({
            "gls": {"path": str(gls_csv_path), "description": "Gene-wise lnL scores, \u0394GLS, support"},
            "gqs": {"path": str(gqs_csv_path), "description": "Gene quartet scores, \u0394GQS, support, status"},
            "consistent_genes": {
                "path": str(output_dir / "consistent_genes.txt"),
                "description": "Loci where GLS and GQS support agree",
            },
            "inconsistent_genes": {
                "path": str(output_dir / "inconsistent_genes.txt"),
                "description": "Loci where GLS and GQS support disagree or ambiguous",
            },
            "gls_support_plot": {"path": str(gls_pdf), "description": "GLS support distribution bar chart"},
            "gqs_support_plot": {"path": str(gqs_pdf), "description": "GQS support distribution bar chart"},
        })
        if len(candidate_trees) > 1:
            output_files["candidate_trees"] = {"path": str(candidate_trees_path), "description": "Merged candidate trees (IQ-TREE -z input)"}

        key_results: dict[str, Any] = {
            "n_loci": len(partition_recs),
            "n_consistent": len(consistent_loci),
            "n_inconsistent": len(inconsistent_loci),
            "n_gqs_skipped": n_gqs_skipped,
            "gls_support_counts": {label: sum(1 for r in gls_rows if r["support"] == label) for label in tree_labels}
                              | {"ambiguous": sum(1 for r in gls_rows if r["support"] == "ambiguous")},
            "gqs_support_counts": {label: sum(1 for r in gqs_rows if r["support"] == label and r["status"] == "success") for label in tree_labels}
                              | {"ambiguous": sum(1 for r in gqs_rows if r["support"] == "ambiguous" or r["status"] == "skipped")},
        }

        if metrics:
            try:
                csv_p, pdf_p, sig_info = _compare_groups(
                    consistent_loci,
                    inconsistent_loci,
                    metrics.resolve(),
                    "consistent",
                    "inconsistent",
                    output_dir,
                )
            except Exception as exc:
                result = {
                    "status": "error",
                    "command": full_command,
                    "wall_time": time.time() - run_start,
                    "tool_versions": tool_versions,
                    "params": params,
                    "key_results": {},
                    "error": f"--metrics comparison failed: {exc}",
                    "error_category": "output",
                    "data": {"cmd": cmd, "tool_stderr": tool_stderr,
                             "tool_log": str(tool_log_path), "output_files": output_files},
                }
                _write_result_json(output_dir, result)
                return result
            output_files["consistent_comparison"] = {
                "path": str(csv_p),
                "description": "Consistent vs inconsistent per-metric means and Wilcoxon p-values",
            }
            output_files["consistent_comparison_plot"] = {
                "path": str(pdf_p),
                "description": "Consistent vs inconsistent metric distribution boxplots",
            }
            key_results.update(n_sig_metrics_consistent=sig_info["n_sig_metrics"],
                               sig_metric_names_consistent=sig_info["sig_metric_names"])

    finally:
        shutil.rmtree(work_dir, ignore_errors=True)

    wall_time = time.time() - run_start
    result = {
        "status": "success",
        "command": full_command,
        "wall_time": wall_time,
        "tool_versions": tool_versions,
        "params": params,
        "key_results": key_results,
        "error": None,
        "data": {
            "cmd": cmd,
            "tool_stderr": "",
            "tool_log": str(tool_log_path),
            "summary": {
                **key_results,
                "wastral_n_gene_trees": len(partition_recs),
                "wastral_threads_used": n_workers,
            },
            "output_files": output_files,
        },
    }
    _write_result_json(output_dir, result)
    return result


def _read_matrix_taxa(matrix: Path) -> set[str]:
    alignment = FormatConverter().read(matrix)
    return {record.id for record in alignment}


def _parse_fclm_report(report_path: Path) -> dict[str, Any]:
    """Extract FcLM region proportions from an IQ-TREE .iqtree report.

    Returns dict with keys: n_quartets, regions (list of {region, percent, count}).
    Returns empty dict if the report cannot be parsed.
    """
    import re as _re

    try:
        text = report_path.read_text()
    except OSError:
        return {}

    n_quartets = None
    regions = []
    in_lmap = False
    region_re = _re.compile(
        r"^\s*(?P<region>\d+[a-z]?)\..*?[:\s]+\s*(?P<percent>[\d.]+)%\s*\((?P<count>\d+)\)"
    )
    for line in text.splitlines():
        stripped = line.strip()
        if stripped.startswith("Likelihood mapping analysis"):
            in_lmap = True
            continue
        if in_lmap:
            if stripped.startswith("Number of quartets"):
                try:
                    n_quartets = int(stripped.split(":")[-1].strip().split()[0])
                except (ValueError, IndexError):
                    pass
                continue
            m = region_re.match(stripped)
            if m:
                regions.append({
                    "region": m.group("region"),
                    "percent": float(m.group("percent")),
                    "count": int(m.group("count")),
                })
            elif regions and not stripped:
                break
    return {"n_quartets": n_quartets, "regions": regions}


def _csv_to_taxset_map(taxset_csv: Path) -> dict[str, list[str]]:
    """Parse taxset CSV into {taxset: [taxa, ...]} without writing any files."""
    taxset_map: dict[str, list[str]] = {}
    with open(taxset_csv, newline="") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            taxon = row["taxon"].strip()
            taxset = row["taxset"].strip()
            taxset_map.setdefault(taxset, []).append(taxon)
    return taxset_map


def _csv_to_nexus(taxset_csv: Path, output_path: Path) -> dict[str, list[str]]:
    taxset_map: dict[str, list[str]] = {}
    with open(taxset_csv, newline="") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            taxon = row["taxon"].strip()
            taxset = row["taxset"].strip()
            taxset_map.setdefault(taxset, []).append(taxon)
    lines = ["#NEXUS", "begin sets;"]
    for ts, taxa in taxset_map.items():
        lines.append(f"  taxset {ts} = {' '.join(taxa)};")
    lines.append("end;")
    output_path.write_text("\n".join(lines) + "\n")
    return taxset_map


def _validate_fclm_inputs(
    *,
    matrix: Path,
    taxset_csv: Path,
    model_expr: str | None = None,
    partitions: Path | None = None,
    partition_mode: str | None = None,
    tool_args: str | None = None,
    guide_tree: Path | None = None,
    lmap: str | None = None,
    threads: str = "auto",
) -> list[str]:
    errors: list[str] = []
    if not matrix.exists() or not matrix.is_file():
        errors.append(f"--matrix does not exist or is not a regular file: {matrix}")
    if not taxset_csv.exists():
        errors.append(f"--taxset-csv does not exist: {taxset_csv}")
        return errors
    if not taxset_csv.is_file():
        errors.append(f"--taxset-csv is not a regular file: {taxset_csv}")
    if guide_tree and (not guide_tree.exists() or not guide_tree.is_file()):
        errors.append(f"--guide-tree does not exist or is not a regular file: {guide_tree}")
    if partitions and (not partitions.exists() or not partitions.is_file()):
        errors.append(f"--partitions does not exist or is not a regular file: {partitions}")
    if partition_mode and not partitions:
        errors.append("--partition-mode is only valid when --partitions is provided")
    if lmap is not None and lmap != "ALL":
        try:
            lmap_val = int(lmap)
            if lmap_val < 1:
                errors.append(f"--lmap must be ALL or a positive integer, got {lmap_val}")
        except ValueError:
            errors.append(f"--lmap must be ALL or a positive integer, got {lmap!r}")

    try:
        matrix_taxa = _read_matrix_taxa(matrix)
    except (OSError, ValueError) as exc:
        errors.append(f"Failed to parse --matrix: {exc}")
        return errors
    csv_taxa: dict[str, str] = {}
    taxsets: dict[str, list[str]] = {}
    try:
        with open(taxset_csv, newline="") as fh:
            reader = csv.DictReader(fh)
            for row in reader:
                taxon = row["taxon"].strip()
                taxset = row["taxset"].strip()
                if taxon in csv_taxa:
                    errors.append(f"Taxon {taxon!r} appears in multiple taxsets (not mutually exclusive)")
                csv_taxa[taxon] = taxset
                taxsets.setdefault(taxset, []).append(taxon)
    except Exception as exc:
        errors.append(f"Cannot read --taxset-csv: {exc}")
        return errors

    extra_taxa = set(csv_taxa) - matrix_taxa
    if extra_taxa:
        errors.append(f"Taxa in --taxset-csv not found in --matrix: {', '.join(sorted(extra_taxa))}")
    missing_taxa = matrix_taxa - set(csv_taxa)
    if missing_taxa:
        errors.append(f"Taxa in --matrix not assigned in --taxset-csv: {', '.join(sorted(missing_taxa))}")
    if len(taxsets) < 4:
        errors.append(f"FcLM requires at least 4 taxsets, got {len(taxsets)}")

    has_model = model_expr or partitions
    has_tool_model = False
    if tool_args:
        toks = set(shlex.split(tool_args))
        has_tool_model = "-m" in toks or "-p" in toks or "-Q" in toks
        for flag in _FCLM_BLOCKED_FLAGS:
            if flag in toks:
                errors.append(f"Blocked flag in --tool-args: {flag}")
    if not has_model and not has_tool_model:
        errors.append("Must specify --model-expr, --partitions, or -m/-p/-Q in --tool-args")
    if threads != "auto":
        try:
            n = int(threads)
            if n < 1:
                errors.append(f"--threads must be a positive integer or 'auto', got {threads!r}")
        except ValueError:
            errors.append(f"--threads must be a positive integer or 'auto', got {threads!r}")
    return errors


def run_signal_fclm(
    *,
    matrix: Path,
    taxset_csv: Path,
    model_expr: str | None = None,
    partitions: Path | None = None,
    partition_mode: str | None = None,
    lmap: str | None = None,
    guide_tree: Path | None = None,
    threads: str = "auto",
    iqtree_path: str | None = None,
    tool_args: str | None = None,
    prefix: str = "fclm",
    resume: bool = False,
    output_dir: Path | None = None,
    overwrite: bool = False,
    dry_run: bool = False,
    quiet: bool = False,
) -> dict[str, Any]:
    run_start = time.time()
    if output_dir is None:
        output_dir = Path("runs/posttree/signal/fclm")
    output_dir = output_dir.resolve()
    matrix = matrix.resolve()
    taxset_csv = taxset_csv.resolve()

    errors = _validate_fclm_inputs(
        matrix=matrix, taxset_csv=taxset_csv,
        model_expr=model_expr, partitions=partitions,
        partition_mode=partition_mode,
        tool_args=tool_args, guide_tree=guide_tree, lmap=lmap,
        threads=threads,
    )

    params: dict[str, Any] = {
        "matrix": str(matrix),
        "taxset_csv": str(taxset_csv),
        "model_expr": model_expr if model_expr else None,
        "partitions": str(partitions.resolve()) if partitions else None,
        "partition_mode": partition_mode if partitions else None,
        "lmap": lmap,
        "guide_tree": str(guide_tree.resolve()) if guide_tree else None,
        "threads": threads,
        "iqtree_path": iqtree_path,
        "tool_args": tool_args,
        "prefix": prefix,
        "resume": resume,
        "output_dir": str(output_dir),
        "overwrite": overwrite,
        "dry_run": dry_run,
        "quiet": quiet,
    }

    if errors:
        return {"status": "error", "command": "", "wall_time": 0.0,
                "tool_versions": {}, "params": params, "key_results": {},
                "error": "; ".join(errors), "error_category": "input",
                "data": {"cmd": [], "tool_stderr": "", "tool_log": None, "output_files": {}}}

    if overwrite and resume:
        return {"status": "error", "command": "", "wall_time": 0.0,
                "tool_versions": {}, "params": params, "key_results": {},
                "error": "--overwrite and --resume are mutually exclusive",
                "error_category": "input",
                "data": {"cmd": [], "tool_stderr": "", "tool_log": None, "output_files": {}}}

    if not dry_run:
        if overwrite and output_dir.exists():
            shutil.rmtree(output_dir)
        elif not resume and output_dir.exists() and any(output_dir.iterdir()):
            return {"status": "error", "command": "", "wall_time": 0.0,
                    "tool_versions": {}, "params": params, "key_results": {},
                    "error": f"Output directory '{output_dir}' already exists and is non-empty. Use --overwrite to replace it.",
                    "error_category": "input",
                    "data": {"cmd": [], "tool_stderr": "", "tool_log": None, "output_files": {}}}
        output_dir.mkdir(parents=True, exist_ok=True)

    nexus_path = output_dir / "cluster.nexus"
    if not dry_run:
        try:
            taxset_map = _csv_to_nexus(taxset_csv, nexus_path)
        except (OSError, ValueError, csv.Error) as exc:
            return {"status": "error", "command": "", "wall_time": 0.0,
                    "tool_versions": {}, "params": params, "key_results": {},
                    "error": f"Failed to read --taxset-csv: {exc}",
                    "error_category": "input",
                    "data": {"cmd": [], "tool_stderr": "", "tool_log": None, "output_files": {}}}
    else:
        try:
            taxset_map = _csv_to_taxset_map(taxset_csv)
        except (OSError, ValueError, csv.Error) as exc:
            return {"status": "error", "command": "", "wall_time": 0.0,
                    "tool_versions": {}, "params": params, "key_results": {},
                    "error": f"Failed to read --taxset-csv: {exc}",
                    "error_category": "input",
                    "data": {"cmd": [], "tool_stderr": "", "tool_log": None, "output_files": {}}}
    n_taxsets = len(taxset_map)

    try:
        matrix_taxa = _read_matrix_taxa(matrix)
        n_taxa = len(matrix_taxa)
    except (OSError, ValueError) as exc:
        return {"status": "error", "command": "", "wall_time": 0.0,
                "tool_versions": {}, "params": params, "key_results": {},
                "error": f"Failed to parse --matrix: {exc}",
                "error_category": "input",
                "data": {"cmd": [], "tool_stderr": "", "tool_log": None, "output_files": {}}}

    if lmap is None:
        lmap_val = str(50 * n_taxa)
    else:
        lmap_val = lmap

    try:
        iqtree_exe = _resolve_iqtree_path(iqtree_path, dry_run)
    except (ValueError, FileNotFoundError) as e:
        return {"status": "error", "command": "", "wall_time": 0.0,
                "tool_versions": {}, "params": params, "key_results": {},
                "error": str(e), "error_category": "env",
                "data": {"cmd": [], "tool_stderr": "", "tool_log": None, "output_files": {}}}

    tool_versions = _detect_iqtree_version(iqtree_exe) if not dry_run else {"iqtree3": "dry-run"}

    tool_toks = set(shlex.split(tool_args)) if tool_args else set()
    cmd = [iqtree_exe, "-s", str(matrix)]
    if model_expr and "-m" not in tool_toks:
        cmd.extend(["-m", model_expr])
    if partitions:
        p_flag = f"-{partition_mode or 'p'}"
        if p_flag not in tool_toks:
            cmd.extend([p_flag, str(partitions.resolve())])
    if guide_tree and "-ft" not in tool_toks:
        cmd.extend(["-ft", str(guide_tree.resolve())])
    cmd.extend(["-lmap", lmap_val, "-lmclust", str(nexus_path), "-n", "0"])
    if "--prefix" not in tool_toks:
        cmd.extend(["--prefix", prefix])
    if "-T" not in tool_toks:
        cmd.extend(["-T", str(threads)])
    if tool_args:
        cmd.extend(shlex.split(tool_args))

    cli_parts = ["phyloai", "posttree", "signal", "fclm",
                 "--matrix", str(matrix), "--taxset-csv", str(taxset_csv)]
    if model_expr:
        cli_parts.extend(["--model-expr", model_expr])
    if partitions:
        cli_parts.extend(["--partitions", str(partitions), "--partition-mode", partition_mode or "p"])
    if lmap:
        cli_parts.extend(["--lmap", lmap])
    if guide_tree:
        cli_parts.extend(["--guide-tree", str(guide_tree)])
    if tool_args:
        cli_parts.extend(["--tool-args", tool_args])
    if prefix != "fclm":
        cli_parts.extend(["--prefix", prefix])
    cli_parts.extend(["--threads", str(threads), "-o", str(output_dir)])
    if iqtree_path:
        cli_parts.extend(["--iqtree-path", iqtree_path])
    if resume:
        cli_parts.append("--resume")
    if overwrite:
        cli_parts.append("--overwrite")
    if dry_run:
        cli_parts.append("--dry-run")
    if quiet:
        cli_parts.append("-q")
    full_command = shlex.join(cli_parts)

    if dry_run:
        return {"status": "success", "command": full_command,
                "wall_time": 0.0, "tool_versions": tool_versions,
                "params": params,
                "key_results": {"n_taxsets": n_taxsets, "n_quartets": lmap_val,
                                "n_taxa": len(matrix_taxa)},
                "error": None,
                "data": {"cmd": cmd, "tool_stderr": "", "tool_log": None,
                         "output_files": {"cluster_nexus": {"path": str(nexus_path),
                                          "description": "NEXUS cluster file for IQ-TREE -lmclust"}}}}

    iqtree_dir = output_dir / "iqtree"
    iqtree_dir.mkdir(parents=True, exist_ok=True)
    proc = subprocess.run(
        cmd, cwd=str(iqtree_dir),
        stdout=None,
        stderr=subprocess.PIPE,
        text=True,
    )
    tool_stderr = proc.stderr.strip() if proc.stderr else ""
    tool_log_path = iqtree_dir / f"{prefix}.log"

    output_files: dict[str, Any] = {
        "cluster_nexus": {"path": str(nexus_path), "description": "NEXUS cluster file for IQ-TREE -lmclust"},
        "lmap_figure": {"path": str(iqtree_dir / f"{prefix}.lmap.svg"),
                        "description": "IQ-TREE likelihood mapping figure (SVG)"},
        "lmap_figure_eps": {"path": str(iqtree_dir / f"{prefix}.lmap.eps"),
                            "description": "IQ-TREE likelihood mapping figure (EPS)"},
        "iqtree_report": {"path": str(iqtree_dir / f"{prefix}.iqtree"),
                          "description": "IQ-TREE native report (contains all lmap statistics)"},
        "iqtree_log": {"path": str(iqtree_dir / f"{prefix}.log"), "description": "IQ-TREE console log"},
    }

    fclm_data: dict[str, Any] = {}
    if proc.returncode == 0:
        lmap_path = iqtree_dir / f"{prefix}.lmap.eps"
        iqtree_report_path = iqtree_dir / f"{prefix}.iqtree"
        if not lmap_path.exists() or not iqtree_report_path.exists():
            status = "error"
            missing = []
            if not lmap_path.exists():
                missing.append(str(lmap_path))
            if not iqtree_report_path.exists():
                missing.append(str(iqtree_report_path))
            error_msg = f"IQ-TREE exited 0 but missing expected output files: {', '.join(missing)}"
        else:
            status = "success"
            error_msg = None
            fclm_data = _parse_fclm_report(iqtree_report_path)
    else:
        status = "error"
        error_msg = f"IQ-TREE exited {proc.returncode}"

    result = {
        "status": status, "command": full_command,
        "wall_time": time.time() - run_start,
        "tool_versions": tool_versions, "params": params,
        "key_results": {"n_taxsets": n_taxsets, "n_quartets": lmap_val,
                        "n_taxa": len(matrix_taxa),
                        "lmap_proportions": fclm_data.get("regions", [])},
        "error": error_msg,
        "error_category": None if status == "success" else "tool",
        "data": {"cmd": cmd, "tool_stderr": tool_stderr if status == "error" else "",
                 "tool_log": str(tool_log_path), "output_files": output_files},
    }
    _write_result_json(output_dir, result)
    return result
