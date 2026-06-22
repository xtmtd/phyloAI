"""Format conversion and sequence normalization for pretree workflows."""

from __future__ import annotations

import shutil
import time
import shlex
from collections.abc import Callable
from pathlib import Path
from typing import Any

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from rich.table import Table

from phyloai.core.formats import AlignmentFormat, FormatConverter
from phyloai.core.sequence_output_validation import validate_fasta_output
from phyloai.core.sequence_normalization import detect_seq_type, expand_dots_from_first_sequence, normalize_sequences


TARGET_SUFFIX = {
    "fasta": ".fa",
    "phylip-relaxed": ".phy",
    "phylip-paml": ".paml.phy",
    "nexus": ".nex",
}


def _build_convert_command(
    input_path: Path,
    output_dir: Path,
    target_format: str,
    input_format: str | None,
    seq_type: str | None,
    aa_special: str,
    threads: int,
    overwrite: bool,
    quiet: bool = False,
) -> str:
    parts = [
        "phyloai", "pretree", "convert",
        "--input", str(input_path),
        "--output-dir", str(output_dir),
        "--to", target_format,
        "--input-format", input_format or "auto",
        "--seq-type", seq_type or "auto",
        "--aa-special", aa_special,
        "--threads", str(threads),
    ]
    if overwrite:
        parts.append("--overwrite")
    if quiet:
        parts.append("--quiet")
    return " ".join(parts)


def convert_input(
    input_path: Path,
    output_dir: Path,
    target_format: str = "fasta",
    input_format: str | None = None,
    seq_type: str | None = None,
    aa_special: str = "x",
    threads: int = 4,
    overwrite: bool = False,
    quiet: bool = False,
    progress_callback: Callable[[Path], None] | None = None,
) -> dict[str, Any]:
    start = time.monotonic()
    if output_dir.exists() and any(output_dir.iterdir()):
        if not overwrite:
            raise ValueError(f"Output directory '{output_dir}' already exists and is non-empty. Use --overwrite to replace it.")
        shutil.rmtree(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    seqs_dir = output_dir / "seqs"
    seqs_dir.mkdir(parents=True, exist_ok=True)
    entries = [input_path] if input_path.is_file() else sorted(input_path.iterdir(), key=lambda path: path.name)
    files: list[dict[str, Any]] = []
    skipped: list[dict[str, str]] = []
    for entry in entries:
        result = _convert_one(entry, seqs_dir, target_format, input_format, seq_type, aa_special)
        if "skipped" in result:
            skipped.append(result["skipped"])
        else:
            files.append(result)
        if progress_callback is not None:
            progress_callback(entry)
    if not files:
        raise ValueError("All input entries failed or were skipped during conversion.")
    total_replacements = sum(sum(item["replacements"].values()) for item in files)
    payload = {
        "status": "success",
        "command": _build_convert_command(input_path, output_dir, target_format, input_format, seq_type, aa_special, threads, overwrite, quiet=quiet),
        "wall_time": round(time.monotonic() - start, 3),
        "tool_versions": {},
        "params": {
            "input_path": str(input_path),
            "output_dir": str(output_dir),
            "target_format": target_format,
            "input_format": input_format or "auto",
            "seq_type": seq_type,
            "aa_special": aa_special,
            "threads": threads,
            "overwrite": overwrite,
            "quiet": quiet,
        },
        "key_results": {},
        "error": None,
        "data": {
            "summary": {
                "n_input_entries": len(entries),
                "n_converted": len(files),
                "n_skipped": len(skipped),
                "target_format": target_format,
                "seq_type_summary": _mixed({item["seq_type"] for item in files}),
                "total_replacements": total_replacements,
                "total_taxon_name_changes": sum(item["taxon_name_changes"] for item in files),
            },
            "files": files,
            "skipped": skipped,
            "warnings": [warning for item in files for warning in item["warnings"]],
        },
    }
    return payload


def _convert_one(entry: Path, output_dir: Path, target_format: str, input_format: str | None, seq_type: str | None, aa_special: str) -> dict[str, Any]:
    if entry.is_dir():
        return {"skipped": {"path": str(entry), "reason": "directory"}}
    if not entry.is_file():
        return {"skipped": {"path": str(entry), "reason": "not a file"}}
    if entry.stat().st_size == 0:
        return {"skipped": {"path": str(entry), "reason": "empty file"}}
    converter = FormatConverter()
    try:
        fmt = converter.detect(entry, declared_format=_alignment_format(input_format))
        records = list(SeqIO.parse(str(entry), "fasta")) if fmt == AlignmentFormat.FASTA else list(converter.read(entry, source_format=fmt))
        if not records:
            return {"skipped": {"path": str(entry), "reason": "no sequences found"}}
        raw_sequences = [str(record.seq) for record in records]
        detected_seq_type = seq_type or detect_seq_type(raw_sequences)
        missing_char = "N" if detected_seq_type == "NT" else "X"
        expanded_sequences, dot_counts = expand_dots_from_first_sequence(raw_sequences, missing_char=missing_char)
        normalized = normalize_sequences(expanded_sequences, detected_seq_type, aa_special=aa_special)
        replacements = {**normalized.replacements}
        for key, value in dot_counts.items():
            replacements[key] = replacements.get(key, 0) + value
        normalized_records = [SeqRecord(Seq(sequence), id=_safe_record_id(record.description, record.id), description="") for record, sequence in zip(records, normalized.sequences)]
        out = output_dir / f"{entry.stem}{TARGET_SUFFIX[target_format]}"
        writer_name_changes = _write_records(normalized_records, out, target_format, detected_seq_type)
        output_warnings = _validate_convert_output(out, target_format)
        if output_warnings:
            return {"skipped": {"path": str(entry), "reason": "; ".join(output_warnings)}}
        generic_name_changes = [
            {"input": original.id, "output": record.id}
            for original, record in zip(records, normalized_records)
            if original.id != record.id
        ]
        taxon_name_changes = _merge_name_changes(generic_name_changes, writer_name_changes)
        return {
            "input": str(entry),
            "output": str(out),
            "input_format": fmt.value,
            "target_format": target_format,
            "seq_type": detected_seq_type,
            "replacements": replacements,
            "taxon_name_changes": len(taxon_name_changes),
            "taxon_name_change_details": taxon_name_changes,
            "warnings": normalized.warnings,
        }
    except Exception as exc:
        return {"skipped": {"path": str(entry), "reason": str(exc)}}


def _write_records(records: list[SeqRecord], out: Path, target_format: str, seq_type: str) -> list[dict[str, str]]:
    out.parent.mkdir(parents=True, exist_ok=True)
    if target_format == "fasta":
        SeqIO.write(records, str(out), "fasta")
        return []
    from Bio.Align import MultipleSeqAlignment
    alignment = MultipleSeqAlignment(records)
    converter = FormatConverter()
    if target_format == "nexus":
        molecule_type = "DNA" if seq_type == "NT" else "protein"
        for record in alignment:
            record.annotations["molecule_type"] = molecule_type
    return converter.write_alignment(alignment, out, target=_alignment_format(target_format), molecule_type="DNA" if seq_type == "NT" else "protein")


def _validate_convert_output(out: Path, target_format: str) -> list[str]:
    if target_format != "fasta":
        return []
    return validate_fasta_output(out, require_aligned=False).warnings


def _merge_name_changes(generic_changes: list[dict[str, str]], writer_changes: list[dict[str, str]]) -> list[dict[str, str]]:
    if not writer_changes:
        return generic_changes
    by_intermediate = {change["input"]: change for change in generic_changes}
    merged: list[dict[str, str]] = []
    writer_inputs = set()
    for writer_change in writer_changes:
        writer_inputs.add(writer_change["input"])
        generic_change = by_intermediate.get(writer_change["input"])
        merged.append({
            "input": generic_change["input"] if generic_change else writer_change["input"],
            "output": writer_change["output"],
        })
    merged.extend(change for change in generic_changes if change["output"] not in writer_inputs)
    return merged


def _alignment_format(value: str | None) -> AlignmentFormat | None:
    if value in {None, "auto"}:
        return None
    for fmt in AlignmentFormat:
        if value == fmt.value:
            return fmt
    raise ValueError(f"Unsupported alignment format '{value}'.")


def _safe_record_id(*candidates: str) -> str:
    for name in candidates:
        joined = "_".join(name.strip().split())
        if joined:
            return joined
    return "taxon"


def _mixed(values: set[str]) -> str | None:
    if not values:
        return None
    if len(values) == 1:
        return next(iter(values))
    return "mixed"


def render_convert_summary_table(summary: dict[str, Any]) -> Table:
    table = Table(title="pretree convert summary")
    table.add_column("Metric")
    table.add_column("Value")
    for key in [
        "n_input_entries",
        "n_converted",
        "n_skipped",
        "target_format",
        "seq_type_summary",
        "total_replacements",
        "total_taxon_name_changes",
    ]:
        table.add_row(key, str(summary[key]))
    return table
