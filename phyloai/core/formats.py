"""Sequence alignment format detection and conversion."""

from __future__ import annotations
from enum import Enum
from pathlib import Path
from typing import Optional

from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord


class AlignmentFormat(str, Enum):
    FASTA = "fasta"
    PHYLIP = "phylip-relaxed"
    PHYLIP_PAML = "phylip-paml"
    NEXUS = "nexus"


_BIOPYTHON_FORMATS = {
    AlignmentFormat.FASTA: "fasta",
    AlignmentFormat.PHYLIP: "phylip-relaxed",
    AlignmentFormat.NEXUS: "nexus",
}


_EXT_MAP: dict[str, AlignmentFormat] = {
    ".fa":    AlignmentFormat.FASTA,
    ".fas":   AlignmentFormat.FASTA,
    ".fasta": AlignmentFormat.FASTA,
    ".faa":   AlignmentFormat.FASTA,
    ".fna":   AlignmentFormat.FASTA,
    ".phy":   AlignmentFormat.PHYLIP,
    ".phylip": AlignmentFormat.PHYLIP,
    ".nex":   AlignmentFormat.NEXUS,
    ".nxs":   AlignmentFormat.NEXUS,
    ".nexus": AlignmentFormat.NEXUS,
}


_COMPOUND_EXT_MAP: dict[str, AlignmentFormat] = {
    ".paml.phy": AlignmentFormat.PHYLIP_PAML,
}


def detect_alignment_format(path: Path, declared_format: Optional[AlignmentFormat] = None) -> AlignmentFormat:
    if declared_format is not None:
        return declared_format
    name = path.name.lower()
    for suffix, fmt in _COMPOUND_EXT_MAP.items():
        if name.endswith(suffix):
            return fmt
    suffix = path.suffix.lower()
    if suffix in _EXT_MAP:
        return _EXT_MAP[suffix]
    if path.exists():
        first = path.read_text(errors="ignore")[:200]
        if first.strip().startswith(">"):
            return AlignmentFormat.FASTA
        if first.strip().upper().startswith("#NEXUS"):
            return AlignmentFormat.NEXUS
        lines = [line for line in first.splitlines() if line.strip()]
        if lines and all(part.isdigit() for part in lines[0].strip().split()[:2]):
            return AlignmentFormat.PHYLIP
    raise ValueError(
        f"Cannot detect alignment format for '{path}'. "
        f"Supported extensions: {list(_EXT_MAP.keys()) + list(_COMPOUND_EXT_MAP.keys())}"
    )


def write_phylip_paml(alignment: MultipleSeqAlignment, dst: Path) -> list[dict[str, str]]:
    dst.parent.mkdir(parents=True, exist_ok=True)
    name_changes = _normalize_paml_names([record.id for record in alignment])
    length = alignment.get_alignment_length()
    with open(dst, "w") as fh:
        fh.write(f"{len(alignment)} {length}\n")
        for record, change in zip(alignment, name_changes):
            fh.write(f"{change['output']}  {str(record.seq)}\n")
    return [change for change in name_changes if change["input"] != change["output"]]


def _normalize_paml_names(names: list[str]) -> list[dict[str, str]]:
    used: set[str] = set()
    changes: list[dict[str, str]] = []
    for raw_name in names:
        safe = _safe_paml_name(raw_name)
        base = safe[:30] or "taxon"
        candidate = base
        suffix_number = 2
        while candidate in used:
            suffix = f"_{suffix_number}"
            candidate = f"{base[:30 - len(suffix)]}{suffix}"
            suffix_number += 1
        used.add(candidate)
        changes.append({"input": raw_name, "output": candidate})
    return changes


def _safe_paml_name(name: str) -> str:
    bad = {'"', ',', ':', '#', '(', ')', '$', '='}
    safe = "".join("_" if char.isspace() or char in bad else char for char in name.strip())
    while "__" in safe:
        safe = safe.replace("__", "_")
    return safe.strip("_")


class FormatConverter:
    def detect(
        self,
        path: Path,
        declared_format: Optional[AlignmentFormat] = None,
    ) -> AlignmentFormat:
        return detect_alignment_format(path, declared_format=declared_format)

    def convert(
        self,
        src: Path,
        dst: Path,
        target: AlignmentFormat,
        source_format: Optional[AlignmentFormat] = None,
        molecule_type: str = "protein",
    ) -> Path:
        src_fmt = source_format or self.detect(src)
        alignment = self.read(src, source_format=src_fmt)
        self.write_alignment(alignment, dst, target=target, molecule_type=molecule_type)
        return dst

    def read(
        self,
        path: Path,
        source_format: Optional[AlignmentFormat] = None,
    ) -> MultipleSeqAlignment:
        fmt = source_format or self.detect(path)
        if fmt == AlignmentFormat.PHYLIP_PAML:
            fmt = AlignmentFormat.PHYLIP
        return AlignIO.read(str(path), _BIOPYTHON_FORMATS[fmt])

    def write_alignment(
        self,
        alignment: MultipleSeqAlignment,
        dst: Path,
        target: AlignmentFormat,
        molecule_type: str = "protein",
    ) -> list[dict[str, str]]:
        dst.parent.mkdir(parents=True, exist_ok=True)
        if target == AlignmentFormat.PHYLIP_PAML:
            return write_phylip_paml(alignment, dst)
        if target == AlignmentFormat.NEXUS:
            for record in alignment:
                if "molecule_type" not in record.annotations:
                    record.annotations["molecule_type"] = molecule_type
        with open(dst, "w") as fh:
            AlignIO.write(alignment, fh, _BIOPYTHON_FORMATS[target])
        return []
