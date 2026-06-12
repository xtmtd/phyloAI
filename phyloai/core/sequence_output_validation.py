"""Validation helpers for generated sequence/alignment files."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from Bio import SeqIO


@dataclass(frozen=True)
class SequenceOutputValidation:
    ok: bool
    n_records: int
    length: int
    warnings: list[str]


def validate_fasta_output(path: Path, require_aligned: bool) -> SequenceOutputValidation:
    if not path.exists() or path.stat().st_size == 0:
        return SequenceOutputValidation(False, 0, 0, ["generated sequence output is empty"])

    try:
        records = list(SeqIO.parse(str(path), "fasta"))
    except Exception as exc:
        return SequenceOutputValidation(False, 0, 0, [f"could not parse generated sequence output as FASTA: {exc}"])

    if not records:
        return SequenceOutputValidation(False, 0, 0, ["generated sequence output contains no FASTA records"])

    lengths = [len(rec.seq) for rec in records]
    warnings: list[str] = []
    if any(length == 0 for length in lengths):
        warnings.append("generated sequence output contains empty sequences")
    if require_aligned and len(set(lengths)) > 1:
        warnings.append(f"generated MSA has unequal sequence lengths: {sorted(set(lengths))}")

    return SequenceOutputValidation(not warnings, len(records), lengths[0], warnings)
