"""Taxon composition heterogeneity diagnostics (taxcomp).

Pure-Python screen for compositional heterogeneity across taxa in one
nucleotide or amino-acid alignment. Provides a Pearson common-composition
chi-square test (overall + per-taxon row contributions with nominal and
Holm-adjusted p-values), a sparse-cell diagnostic, and the observed PPA-COMP
descriptive statistics shared with ``posttree simulate adequacy``.

All p-values are nominal exploratory values; the command makes no taxon
removal, recoding, model, or topology decision.
"""

from __future__ import annotations

import csv
import shlex
import shutil
import time as _time
from pathlib import Path
from typing import Any

from Bio.Align import MultipleSeqAlignment
from scipy.stats import chi2

from phyloai.core.formats import FormatConverter
from phyloai.core.schema import write_result_json
from phyloai.core.sequence_normalization import detect_seq_type
from phyloai.posttree.simulate_adequacy import (
    AA_STATES,
    NT_STATES,
    _compute_taxon_composition,
)

OVERALL_FIELDS = [
    "n_taxa", "n_states", "x2", "df", "p_nominal",
    "sparse_count_check", "expected_cells_total", "expected_cells_below_1",
    "expected_cells_below_5", "expected_cells_below_5_fraction",
    "comp_max", "comp_mean",
]
TAXON_FIELDS = [
    "taxon", "x2_contribution", "df", "p_nominal", "p_holm",
    "squared_composition_distance",
]

_FORMAT_CONVERTER = FormatConverter()


def _resolve_seq_type(matrix: Path, requested: str) -> str:
    """Resolve ``auto`` to the alignment's detected type; validate AA/NT."""
    normalized = requested.upper()
    if normalized == "AUTO":
        alignment = _FORMAT_CONVERTER.read(matrix)
        return detect_seq_type([str(record.seq) for record in alignment])
    if normalized not in {"AA", "NT"}:
        raise ValueError(f"invalid seq_type: {requested!r}")
    return normalized


def _read_alignment(matrix: Path) -> MultipleSeqAlignment:
    try:
        alignment = _FORMAT_CONVERTER.read(matrix)
    except Exception as exc:
        raise ValueError(f"unable to parse alignment file {matrix}: {exc}") from exc
    if not len(alignment) or not alignment.get_alignment_length():
        raise ValueError(f"alignment file {matrix} is empty")
    return alignment


def compute_taxcomp_statistics(
    alignment: MultipleSeqAlignment,
    seq_type: str,
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    """Compute Pearson, Holm, sparse-cell, and PPA-COMP summaries.

    Returns ``(overall, rows)`` where ``overall`` carries the OVERALL_FIELDS
    columns and ``rows`` one TAXON_FIELDS row per taxon in input order.
    """
    names = [record.id for record in alignment]
    duplicates = sorted({name for name in names if names.count(name) > 1})
    if duplicates:
        raise ValueError(f"duplicate taxon IDs: {', '.join(duplicates)}")
    if len(names) < 2:
        raise ValueError("at least two taxa are required")
    if len({len(record.seq) for record in alignment}) != 1:
        raise ValueError("alignment sequences have unequal lengths")

    states = AA_STATES if seq_type == "AA" else NT_STATES
    state_index = {state: index for index, state in enumerate(states)}
    n_full = len(states)

    counts = [[0] * n_full for _ in names]
    for row, record in enumerate(alignment):
        for state in str(record.seq).upper():
            index = state_index.get(state)
            if index is not None:
                counts[row][index] += 1

    row_totals = [sum(row) for row in counts]
    for name, total in zip(names, row_totals):
        if not total:
            raise ValueError(f"taxon {name!r} has no valid characters")

    column_totals = [
        sum(counts[row][col] for row in range(len(names)))
        for col in range(n_full)
    ]
    retained = [col for col in range(n_full) if column_totals[col] > 0]
    if len(retained) < 2:
        raise ValueError("at least two globally observed standard states are required")

    n_taxa = len(names)
    n_states = len(retained)
    grand_total = sum(row_totals)

    row_x2: list[float] = []
    expected_cells: list[float] = []
    for row in range(n_taxa):
        contribution = 0.0
        for col in retained:
            expected = row_totals[row] * column_totals[col] / grand_total
            expected_cells.append(expected)
            observed = counts[row][col]
            contribution += (observed - expected) ** 2 / expected
        row_x2.append(contribution)

    x2_overall = sum(row_x2)
    df_overall = (n_taxa - 1) * (n_states - 1)
    df_taxon = n_states - 1
    p_overall = float(chi2.sf(x2_overall, df_overall))
    p_nominal_rows = [float(chi2.sf(value, df_taxon)) for value in row_x2]
    p_holm = holm_adjust(p_nominal_rows)

    sparse = sparse_count_check(expected_cells)

    composition = _compute_taxon_composition(alignment, seq_type)
    taxon_dist = composition["taxon_dist_j"]

    rows = [
        {
            "taxon": name,
            "x2_contribution": row_x2[index],
            "df": df_taxon,
            "p_nominal": p_nominal_rows[index],
            "p_holm": p_holm[index],
            "squared_composition_distance": taxon_dist[name],
        }
        for index, name in enumerate(names)
    ]
    overall = {
        "n_taxa": n_taxa,
        "n_states": n_states,
        "x2": x2_overall,
        "df": df_overall,
        "p_nominal": p_overall,
        **sparse,
        "comp_max": composition["comp_max"],
        "comp_mean": composition["comp_mean"],
    }
    return overall, rows


def holm_adjust(p_values: list[float]) -> list[float]:
    """Holm step-down adjustment preserving input order (1-based ranks)."""
    if not p_values:
        return []
    ordered = sorted(range(len(p_values)), key=lambda i: p_values[i])
    adjusted_sorted: list[float] = []
    running_max = 0.0
    m = len(p_values)
    for rank, original_index in enumerate(ordered, start=1):
        adjusted = (m - rank + 1) * p_values[original_index]
        running_max = max(running_max, adjusted)
        adjusted_sorted.append(min(1.0, running_max))
    restored = [0.0] * len(p_values)
    for sorted_index, original_index in enumerate(ordered):
        restored[original_index] = adjusted_sorted[sorted_index]
    return restored


def sparse_count_check(expected_cells: list[float]) -> dict[str, int | float | str]:
    """Evaluate the conventional sparse-cell rule with strict boundaries."""
    total = len(expected_cells)
    below_1 = sum(1 for value in expected_cells if value < 1)
    below_5 = sum(1 for value in expected_cells if value < 5)
    fraction = below_5 / total if total else 0.0
    triggered = below_1 > 0 or fraction > 0.2
    return {
        "sparse_count_check": "triggered" if triggered else "not_triggered",
        "expected_cells_total": total,
        "expected_cells_below_1": below_1,
        "expected_cells_below_5": below_5,
        "expected_cells_below_5_fraction": fraction,
    }


def build_taxcomp_command(
    matrix: Path,
    seq_type: str,
    table_format: str,
    output_dir: Path,
    overwrite: bool,
    dry_run: bool,
    quiet: bool,
) -> str:
    parts = [
        "phyloai", "posttree", "syserror", "taxcomp",
        "--matrix", str(matrix.resolve()),
        "--seq-type", seq_type,
        "--table-format", table_format,
        "--output-dir", str(output_dir.resolve()),
    ]
    if overwrite:
        parts.append("--overwrite")
    if dry_run:
        parts.append("--dry-run")
    if quiet:
        parts.append("--quiet")
    return shlex.join(parts)


def run_taxcomp(
    matrix: Path,
    seq_type: str = "auto",
    table_format: str = "csv",
    output_dir: Path = Path("runs/posttree/syserror/taxcomp"),
    overwrite: bool = False,
    dry_run: bool = False,
    quiet: bool = False,
) -> dict[str, Any]:
    """Run the taxcomp diagnostic; returns the standard result.json payload.

    Raises ``ValueError`` for validation and output-conflict failures; the CLI
    handler decides whether an error ``result.json`` may be written.
    """
    start = _time.monotonic()

    if table_format not in {"csv", "tsv"}:
        raise ValueError(f"invalid table_format: {table_format!r}")
    if seq_type.upper() not in {"AA", "NT", "AUTO"}:
        raise ValueError(f"invalid seq_type: {seq_type!r}")
    if not matrix.exists():
        raise ValueError(f"--matrix does not exist: {matrix}")
    if not matrix.is_file():
        raise ValueError(f"--matrix is not a file: {matrix}")

    resolved_seq_type = _resolve_seq_type(matrix, seq_type)
    alignment = _read_alignment(matrix)
    overall, rows = compute_taxcomp_statistics(alignment, resolved_seq_type)

    output_dir = output_dir.resolve()
    params = {
        "matrix": str(matrix.resolve()),
        "seq_type": seq_type,
        "detected_seq_type": resolved_seq_type,
        "table_format": table_format,
        "output_dir": str(output_dir),
        "overwrite": overwrite,
        "dry_run": dry_run,
        "quiet": quiet,
    }
    alphabet = AA_STATES if resolved_seq_type == "AA" else NT_STATES
    largest_x2 = max(rows, key=lambda row: row["x2_contribution"])
    largest_dist = max(rows, key=lambda row: row["squared_composition_distance"])
    key_results = {
        **overall,
        "n_sites": alignment.get_alignment_length(),
        "seq_type": resolved_seq_type,
        "largest_x2_contribution_taxon": largest_x2["taxon"],
        "largest_squared_composition_distance_taxon": largest_dist["taxon"],
    }
    command = build_taxcomp_command(
        matrix, resolved_seq_type, table_format, output_dir, overwrite, dry_run, quiet,
    )
    warnings: list[str] = []
    if overall["sparse_count_check"] == "triggered":
        warnings.append(
            "sparse expected counts make the nominal chi-square p-values "
            "especially unreliable"
        )
    payload: dict[str, Any] = {
        "status": "success",
        "command": command,
        "wall_time": 0.0,
        "tool_versions": {},
        "params": params,
        "key_results": key_results,
        "error": None,
        "error_category": None,
        "data": {
            "cmd": [],
            "tool_stderr": "",
            "warnings": warnings,
            "output_files": {},
            "character_policy": {
                "alphabet": alphabet,
                "seq_type": resolved_seq_type,
                "missing_characters": (
                    "all non-standard characters (gaps, unknowns, ambiguity "
                    "codes, stops) were excluded as missing"
                ),
            },
        },
    }
    if dry_run:
        payload["wall_time"] = round(_time.monotonic() - start, 3)
        return payload

    if output_dir.exists() and not output_dir.is_dir():
        raise ValueError(f"output path is not a directory: {output_dir}")
    if output_dir.is_dir() and any(output_dir.iterdir()):
        if not overwrite:
            raise ValueError(f"output directory is not empty: {output_dir}")
        shutil.rmtree(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    delimiter = "\t" if table_format == "tsv" else ","
    suffix = ".tsv" if table_format == "tsv" else ".csv"
    overall_path = output_dir / f"overall_summary{suffix}"
    taxon_path = output_dir / f"taxon_summary{suffix}"

    with overall_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=OVERALL_FIELDS, extrasaction="raise", delimiter=delimiter)
        writer.writeheader()
        writer.writerow(overall)
    with taxon_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=TAXON_FIELDS, extrasaction="raise", delimiter=delimiter)
        writer.writeheader()
        writer.writerows(rows)

    payload["data"]["output_files"] = {
        "overall_summary": {
            "path": str(overall_path.resolve()),
            "description": "Overall Pearson common-composition X2, sparse-cell diagnostic, and PPA-COMP comp_max/comp_mean",
        },
        "taxon_summary": {
            "path": str(taxon_path.resolve()),
            "description": "Per-taxon X2 contribution, nominal and Holm-adjusted p-values, and squared composition distance",
        },
    }
    payload["wall_time"] = round(_time.monotonic() - start, 3)
    write_result_json(payload, output_dir)

    if not quiet:
        click_echo = __import__("click").echo
        click_echo("Taxon composition heterogeneity screen:")
        click_echo(f"  Sequence type: {resolved_seq_type}")
        click_echo(f"  Taxa: {overall['n_taxa']}  States: {overall['n_states']}  Sites: {alignment.get_alignment_length()}")
        click_echo(f"  Overall X2: {overall['x2']}  df: {overall['df']}  p (nominal): {overall['p_nominal']}")
        click_echo(f"  sparse-cell rule {overall['sparse_count_check']}")
        click_echo(f"  comp_max: {overall['comp_max']}  comp_mean: {overall['comp_mean']}")
        for label, info in payload["data"]["output_files"].items():
            click_echo(f"  {label}: {info['path']}")
        click_echo(f"  Result: {output_dir / 'result.json'}")
    return payload
