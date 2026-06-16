"""Sequence and alignment statistics for pretree workflows."""

from __future__ import annotations

import csv
import json
import statistics
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from typing import Any

from Bio import SeqIO
from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord
from rich.panel import Panel
from rich.table import Table

from phyloai.core.formats import AlignmentFormat, FormatConverter
from phyloai.core.schema import COMMON_ALIGNMENT_EXTENSIONS
from phyloai.core.sequence_normalization import (
    gap_chars,
    resolve_seq_type,
    standard_chars,
)


PER_GENE_COLUMNS = [
    "gene",
    "n_taxa",
    "n_taxa_ratio",
    "length_type",
    "alignment_length",
    "seq_length_min",
    "seq_length_max",
    "seq_length_mean",
    "seq_length_median",
    "seq_length_stdev",
    "gap_ratio",
    "ambiguous_ratio",
    "gap_ambiguous_ratio",
    "missing_taxa",
    "missing_taxa_ratio",
]


def _is_blank(value: Any) -> bool:
    return value in {"", None}


# Columns that only apply to one alignment regime.
_ALIGNED_ONLY = frozenset({"alignment_length"})
_UNALIGNED_ONLY = frozenset({
    "seq_length_min", "seq_length_max",
    "seq_length_mean", "seq_length_median", "seq_length_stdev",
})
_SITE_PATTERN_KEYS = (
    "distinct_patterns", "constant_sites",
    "parsimony_informative", "singleton_sites",
)


def per_gene_columns_for_rows(
    rows: list[dict[str, Any]],
    is_aligned: bool = True,
) -> list[str]:
    """Return the ``PER_GENE_COLUMNS`` subset relevant for *is_aligned*.

    When *is_aligned* is ``True`` (default) the output includes
    ``alignment_length`` and site-pattern columns (when present) but
    excludes per-sequence length statistics.  When ``False`` the output
    includes the per-sequence length columns but excludes alignment-
    specific fields.
    """
    if not rows:
        return list(PER_GENE_COLUMNS)

    excluded = _UNALIGNED_ONLY if is_aligned else _ALIGNED_ONLY
    columns = [
        c
        for c in PER_GENE_COLUMNS
        if c not in excluded and any(not _is_blank(r.get(c, "")) for r in rows)
    ]

    if is_aligned:
        for sp_key in _SITE_PATTERN_KEYS:
            if any(r.get(sp_key) not in ("", None) for r in rows):
                columns.append(sp_key)

    return columns


def _median(values: list[float]) -> float:
    return statistics.median(values) if values else 0.0


def _mean(values: list[float]) -> float:
    return statistics.mean(values) if values else 0.0


def _stdev(values: list[float]) -> float:
    return statistics.stdev(values) if len(values) > 1 else 0.0


def _round(value: float) -> float:
    return round(value, 6)


def _read_records(path: Path, input_format: str | AlignmentFormat | None = None) -> tuple[list[SeqRecord], AlignmentFormat]:
    converter = FormatConverter()
    declared_format = _normalize_format(input_format)
    fmt = converter.detect(path, declared_format=declared_format)
    if fmt == AlignmentFormat.FASTA:
        records = list(SeqIO.parse(str(path), fmt.value))
    else:
        alignment: MultipleSeqAlignment = converter.read(path, source_format=fmt)
        records = list(alignment)
    if not records:
        raise ValueError(f"No sequences found in '{path}'.")
    return records, fmt


def _normalize_format(input_format: str | AlignmentFormat | None) -> AlignmentFormat | None:
    if input_format is None or isinstance(input_format, AlignmentFormat):
        return input_format
    normalized = input_format.lower()
    for fmt in AlignmentFormat:
        if normalized == fmt.value or normalized == fmt.name.lower():
            return fmt
    raise ValueError(f"Unsupported input format '{input_format}'.")


def _sequence_strings(records: list[SeqRecord]) -> list[str]:
    return [str(record.seq).upper() for record in records]


def check_stop_codons(sequences: list[str], filename: str) -> list[str]:
    if any("*" in sequence for sequence in sequences):
        return [
            f"[WARN] Stop codon (*) found in {filename}. This may indicate upstream processing errors."
        ]
    return []


def per_taxon_stats(record: SeqRecord, seq_type: str) -> dict[str, Any]:
    sequence = str(record.seq).upper()
    raw_length = len(sequence)
    if raw_length == 0:
        raise ValueError(f"Sequence '{record.id}' is empty.")
    gap_count = sum(sequence.count(char) for char in gap_chars(seq_type))
    standard_count = sum(sequence.count(char) for char in standard_chars(seq_type))
    ambiguous_count = raw_length - gap_count - standard_count
    ungapped_length = raw_length - gap_count
    return {
        "name": record.id,
        "raw_length": raw_length,
        "ungapped_length": ungapped_length,
        "gap_ratio": _round(gap_count / raw_length),
        "ambiguous_ratio": _round(ambiguous_count / raw_length),
        "standard_ratio": _round(standard_count / raw_length),
    }


def compute_site_patterns(sequences: list[str], seq_type: str) -> dict[str, Any]:
    alignment_length = len(sequences[0]) if sequences else 0
    parsimony_informative = 0
    singleton_sites = 0
    standard_codes = {ord(char) for char in standard_chars(seq_type)}
    gap_code = ord("-")
    distinct_pattern_set: set[bytes] = set()
    for column in zip(*(sequence.encode("ascii") for sequence in sequences)):
        distinct_pattern_set.add(bytes(gap_code if ch not in standard_codes else ch for ch in column))

        standard_counts: dict[int, int] = {}
        for char in column:
            if char in standard_codes:
                standard_counts[char] = standard_counts.get(char, 0) + 1
        if sum(standard_counts.values()) < 2:
            continue
        if len(standard_counts) == 1:
            continue
        repeated = sum(1 for value in standard_counts.values() if value >= 2)
        if repeated >= 2:
            parsimony_informative += 1
        else:
            singleton_sites += 1
    variable_sites = parsimony_informative + singleton_sites
    constant_sites = alignment_length - variable_sites
    distinct_patterns = len(distinct_pattern_set)
    return {
        "alignment_length": alignment_length,
        "distinct_patterns": _count_and_ratio(distinct_patterns, alignment_length),
        "constant_sites": _count_and_ratio(constant_sites, alignment_length),
        "parsimony_informative": _count_and_ratio(parsimony_informative, alignment_length),
        "singleton_sites": _count_and_ratio(singleton_sites, alignment_length),
    }


def _count_and_ratio(count: int, total: int) -> dict[str, Any]:
    ratio = 0.0 if total == 0 else _round(count / total)
    return {"count": count, "ratio": ratio}


def _summarize_per_taxon(per_taxon: list[dict[str, Any]]) -> dict[str, float]:
    return {
        "gap_ratio": _round(_mean([entry["gap_ratio"] for entry in per_taxon])),
        "ambiguous_ratio": _round(_mean([entry["ambiguous_ratio"] for entry in per_taxon])),
        "gap_ambiguous_ratio": _round(_mean([entry["gap_ratio"] + entry["ambiguous_ratio"] for entry in per_taxon])),
        "standard_ratio": _round(_mean([entry["standard_ratio"] for entry in per_taxon])),
    }


def file_stats_unaligned(
    path: Path,
    seq_type: str | None = None,
    input_format: str | AlignmentFormat | None = None,
) -> dict[str, Any]:
    records, fmt = _read_records(path, input_format)
    sequences = _sequence_strings(records)
    detected_seq_type, seq_type_warnings = (seq_type, []) if seq_type else resolve_seq_type(sequences)
    warnings = seq_type_warnings + check_stop_codons(sequences, path.name)
    per_taxon = [per_taxon_stats(record, detected_seq_type) for record in records]
    ungapped_lengths = [entry["ungapped_length"] for entry in per_taxon]
    summary = _summarize_per_taxon(per_taxon)
    return {
        "filename": str(path),
        "gene": path.stem,
        "format": fmt.value,
        "seq_type": detected_seq_type,
        "is_aligned": False,
        "n_taxa": len(records),
        "taxon_names": [record.id for record in records],
        "character_summary": summary,
        "per_taxon": per_taxon,
        "seq_length": {
            "min": min(ungapped_lengths),
            "max": max(ungapped_lengths),
            "mean": _round(_mean(ungapped_lengths)),
            "median": _round(_median(ungapped_lengths)),
            "stdev": _round(_stdev(ungapped_lengths)),
        },
        "total_length": sum(ungapped_lengths),
        "length": _round(_median(ungapped_lengths)),
        "length_type": "seq_length",
        "alignment_length": "",
        "seq_length_min": min(ungapped_lengths),
        "seq_length_max": max(ungapped_lengths),
        "seq_length_mean": _round(_mean(ungapped_lengths)),
        "seq_length_median": _round(_median(ungapped_lengths)),
        "seq_length_stdev": _round(_stdev(ungapped_lengths)),
        "gap_ratio": summary["gap_ratio"],
        "ambiguous_ratio": summary["ambiguous_ratio"],
        "gap_ambiguous_ratio": summary["gap_ambiguous_ratio"],
        "standard_ratio": summary["standard_ratio"],
        "warnings": warnings,
    }


def file_stats_aligned(
    path: Path,
    seq_type: str | None = None,
    input_format: str | AlignmentFormat | None = None,
) -> dict[str, Any]:
    records, fmt = _read_records(path, input_format)
    sequences = _sequence_strings(records)
    lengths = {len(sequence) for sequence in sequences}
    if len(lengths) != 1:
        raise ValueError(f"Alignment '{path}' contains unequal sequence lengths.")
    detected_seq_type, seq_type_warnings = (seq_type, []) if seq_type else resolve_seq_type(sequences)
    warnings = seq_type_warnings + check_stop_codons(sequences, path.name)
    per_taxon = [per_taxon_stats(record, detected_seq_type) for record in records]
    summary = _summarize_per_taxon(per_taxon)
    patterns = compute_site_patterns(sequences, detected_seq_type)
    return {
        "filename": str(path),
        "gene": path.stem,
        "format": fmt.value,
        "seq_type": detected_seq_type,
        "is_aligned": True,
        "n_taxa": len(records),
        "taxon_names": [record.id for record in records],
        "character_summary": summary,
        "per_taxon": per_taxon,
        "alignment_length": patterns["alignment_length"],
        "length": patterns["alignment_length"],
        "length_type": "alignment_length",
        "seq_length_min": "",
        "seq_length_max": "",
        "seq_length_mean": "",
        "seq_length_median": "",
        "seq_length_stdev": "",
        "gap_ratio": summary["gap_ratio"],
        "ambiguous_ratio": summary["ambiguous_ratio"],
        "gap_ambiguous_ratio": summary["gap_ambiguous_ratio"],
        "standard_ratio": summary["standard_ratio"],
        "warnings": warnings,
        **{key: value for key, value in patterns.items() if key != "alignment_length"},
    }


def stats_single_file(
    path: Path,
    seq_type: str | None = None,
    input_format: str | AlignmentFormat | None = None,
) -> dict[str, Any]:
    records, fmt = _read_records(path, input_format)
    sequences = _sequence_strings(records)
    lengths = {len(sequence) for sequence in sequences}
    is_aligned = len(records) > 1 and len(lengths) == 1
    if is_aligned:
        return file_stats_aligned(path, seq_type=seq_type, input_format=fmt)
    return file_stats_unaligned(path, seq_type=seq_type, input_format=fmt)


def _worker(args: tuple[Path, str | None, str | AlignmentFormat | None]) -> dict[str, Any]:
    path, seq_type, input_format = args
    try:
        return stats_single_file(path, seq_type=seq_type, input_format=input_format)
    except Exception as exc:  # pragma: no cover - exercised via directory error test later
        return {"gene": path.stem, "filename": str(path), "error": str(exc), "warnings": []}


def collect_seq_files(directory: Path) -> list[Path]:
    return sorted(
        [path for path in directory.iterdir() if path.is_file() and path.suffix.lower() in COMMON_ALIGNMENT_EXTENSIONS],
        key=lambda path: path.name,
    )


def stats_directory(
    directory: Path,
    seq_type: str | None,
    input_format: str | AlignmentFormat | None,
    threads: int,
    progress_callback: Any | None = None,
) -> tuple[list[dict[str, Any]], list[str]]:
    files = collect_seq_files(directory)
    args = [(path, seq_type, input_format) for path in files]
    if threads == 1:
        results = []
        for path, arg in zip(files, args):
            results.append(_worker(arg))
            if progress_callback is not None:
                progress_callback(path)
    else:
        with ProcessPoolExecutor(max_workers=threads) as executor:
            results = []
            for path, result in zip(files, executor.map(_worker, args)):
                results.append(result)
                if progress_callback is not None:
                    progress_callback(path)
    warnings: list[str] = []
    for result in results:
        warnings.extend(result.get("warnings", []))
    return results, warnings


def aggregate_summary(per_file_results: list[dict[str, Any]]) -> dict[str, Any]:
    successful = [result for result in per_file_results if "error" not in result]
    formats = {result["format"] for result in successful}
    aligned_states = {result["is_aligned"] for result in successful}
    seq_types = {result["seq_type"] for result in successful}
    total_taxa = sorted({name for result in successful for name in result["taxon_names"]})
    n_taxa_values = [result["n_taxa"] for result in successful]
    length_values = [result["length"] for result in successful]
    gap_values = [result["gap_ratio"] for result in successful]
    ambiguous_values = [result["ambiguous_ratio"] for result in successful]
    gap_ambiguous_values = [result["gap_ambiguous_ratio"] for result in successful]
    max_taxa = max(n_taxa_values) if n_taxa_values else 0
    per_gene_with_coverage: list[dict[str, Any]] = []
    for result in per_file_results:
        enriched = dict(result)
        if "error" not in result and max_taxa > 0:
            missing_taxa = max_taxa - result["n_taxa"]
            enriched["missing_taxa"] = missing_taxa
            enriched["missing_taxa_ratio"] = _round(missing_taxa / max_taxa)
            enriched["n_taxa_ratio"] = _round(result["n_taxa"] / max_taxa)
        per_gene_with_coverage.append(enriched)
        result.update({
            "missing_taxa": enriched.get("missing_taxa", 0),
            "missing_taxa_ratio": enriched.get("missing_taxa_ratio", 0.0),
            "n_taxa_ratio": enriched.get("n_taxa_ratio", 0.0),
        })
    missing_taxa_ratios = [result["missing_taxa_ratio"] for result in successful] if successful else []
    return {
        "n_genes": len(per_file_results),
        "n_genes_ok": len(successful),
        "n_errors": len(per_file_results) - len(successful),
        "format": _mixed_value(formats),
        "is_aligned": _mixed_value(aligned_states),
        "seq_type": _mixed_value(seq_types),
        "total_taxa": len(total_taxa),
        "taxa_per_gene": _range_summary(n_taxa_values),
        "length": _range_summary(length_values, include_total=True),
        "gap_ratio": _mean_median_summary(gap_values),
        "ambiguous_ratio": _mean_median_summary(ambiguous_values),
        "gap_ambiguous_ratio": _mean_median_summary(gap_ambiguous_values),
        "missing_taxa_ratio": _mean_median_summary(missing_taxa_ratios),
        "warnings": [],
    }


def _mixed_value(values: set[Any]) -> Any:
    if not values:
        return None
    if len(values) == 1:
        return next(iter(values))
    return "mixed"


def _range_summary(values: list[float], include_total: bool = False) -> dict[str, float]:
    if not values:
        summary: dict[str, float] = {"min": 0, "max": 0, "mean": 0.0, "median": 0.0}
    else:
        summary = {
            "min": min(values),
            "max": max(values),
            "mean": _round(_mean(values)),
            "median": _round(_median(values)),
        }
    if include_total:
        summary["total"] = sum(values) if values else 0
    return summary


def _mean_median_summary(values: list[float]) -> dict[str, float]:
    if not values:
        return {"mean": 0.0, "median": 0.0}
    return {"mean": _round(_mean(values)), "median": _round(_median(values))}


def render_summary_table(summary: dict[str, Any]) -> Table:
    table = Table(title="pretree stats summary")
    table.add_column("Metric")
    table.add_column("Value")
    table.add_row("n_genes", str(summary["n_genes"]))
    table.add_row("n_genes_ok", str(summary["n_genes_ok"]))
    table.add_row("format", str(summary["format"]))
    table.add_row("is_aligned", str(summary["is_aligned"]))
    table.add_row("seq_type", str(summary["seq_type"]))
    table.add_row("total_taxa", str(summary["total_taxa"]))
    table.add_row("length_mean", str(summary["length"]["mean"]))
    table.add_row("gap_ratio_mean", str(summary["gap_ratio"]["mean"]))
    table.add_row("ambiguous_ratio_mean", str(summary["ambiguous_ratio"]["mean"]))
    table.add_row("gap_ambiguous_ratio_mean", str(summary["gap_ambiguous_ratio"]["mean"]))
    return table


def render_per_gene_table(per_file: list[dict[str, Any]]) -> Table:
    table = Table(title="Per-gene statistics")
    columns = per_gene_columns_for_rows(per_file)
    for column in columns:
        table.add_column(column)
    for result in sorted(per_file, key=lambda item: item["gene"]):
        table.add_row(*(str(result.get(column, "")) for column in columns))
    return table


def render_single_file_panels(stats: dict[str, Any]) -> list[Panel]:
    overview = Table(show_header=False)
    overview.add_column("Metric")
    overview.add_column("Value")
    for key in ["filename", "format", "seq_type", "is_aligned", "n_taxa"]:
        overview.add_row(key, str(stats[key]))

    character = Table(show_header=False)
    character.add_column("Metric")
    character.add_column("Value")
    for key, value in stats["character_summary"].items():
        character.add_row(key, str(value))

    per_taxon = Table(title="Per-taxon")
    for column in ["name", "raw_length", "ungapped_length", "gap_ratio", "ambiguous_ratio", "gap_ambiguous_ratio"]:
        per_taxon.add_column(column)
    for entry in stats["per_taxon"]:
        per_taxon.add_row(
            str(entry["name"]),
            str(entry["raw_length"]),
            str(entry["ungapped_length"]),
            str(entry["gap_ratio"]),
            str(entry["ambiguous_ratio"]),
            str(_round(entry["gap_ratio"] + entry["ambiguous_ratio"])),
        )

    panels = [
        Panel(overview, title="Overview"),
        Panel(character, title="Character Summary"),
    ]
    if stats["is_aligned"]:
        site_table = Table(title="Site Patterns")
        site_table.add_column("Metric")
        site_table.add_column("Count")
        site_table.add_column("Ratio")
        site_table.add_row("MSA length", str(stats["alignment_length"]), "1.0")
        for key in [
            "distinct_patterns",
            "constant_sites",
            "parsimony_informative",
            "singleton_sites",
        ]:
            site_table.add_row(key, str(stats[key]["count"]), str(stats[key]["ratio"]))
        panels.append(Panel(site_table, title="Site Patterns"))
    panels.append(Panel(per_taxon, title="Per-taxon"))
    return panels


def write_output(data: dict[str, Any], path: Path, mode: str, per_gene: bool = False, force_json: bool = False, output_format: str = "json") -> None:
    """Write stats payload to *path* in the format dictated by *output_format*.

    ``output_format="json"`` (default) writes pretty-printed JSON regardless of
    the file extension.  ``output_format="text"`` writes a human-readable
    key=value text report.  The legacy *force_json* flag is retained for
    backwards compatibility and forces JSON when ``True``.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    if force_json or output_format == "json":
        path.write_text(json.dumps(data, indent=2, sort_keys=True))
        return
    # text format
    if mode == "directory":
        lines = ["[summary]"]
        lines.extend(f"{key}={value}" for key, value in data["data"]["summary"].items())
    else:
        lines = [f"{key}={value}" for key, value in data["data"].items() if key != "per_taxon"]
        lines.append("")
        lines.append("[per_taxon]")
        lines.append("name,raw_length,ungapped_length,gap_ratio,ambiguous_ratio,gap_ambiguous_ratio")
        for entry in data["data"].get("per_taxon", []):
            lines.append(
                ",".join(
                    str(value)
                    for value in [
                        entry.get("name", ""),
                        entry.get("raw_length", ""),
                        entry.get("ungapped_length", ""),
                        entry.get("gap_ratio", ""),
                        entry.get("ambiguous_ratio", ""),
                        _round(entry.get("gap_ratio", 0.0) + entry.get("ambiguous_ratio", 0.0)),
                    ]
                )
            )
    path.write_text("\n".join(lines) + "\n")


def per_gene_output_path(summary_path: Path, output_format: str = "csv") -> Path:
    table_suffix = f".{output_format}"
    return summary_path.with_name(f"{summary_path.stem}.per-gene{table_suffix}")


def write_per_gene_output(data: dict[str, Any], path: Path) -> None:
    suffix = path.suffix.lower()
    if suffix not in {".csv", ".tsv"}:
        raise ValueError(f"Unsupported per-gene output extension '{path.suffix}'.")
    path.parent.mkdir(parents=True, exist_ok=True)
    delimiter = "," if suffix == ".csv" else "\t"
    columns = per_gene_columns_for_rows(data["data"]["per_gene"])
    rows = [{column: result.get(column, "") for column in columns} for result in data["data"]["per_gene"]]
    if not rows:
        path.write_text("")
        return
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter=delimiter)
        writer.writeheader()
        writer.writerows(rows)
