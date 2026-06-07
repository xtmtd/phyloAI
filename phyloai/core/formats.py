"""Sequence alignment format detection and conversion."""

from __future__ import annotations
from enum import Enum
from pathlib import Path
from typing import Optional

from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment


class AlignmentFormat(str, Enum):
    FASTA = "fasta"
    PHYLIP = "phylip-relaxed"
    PHYLIP_PAML = "phylip"
    NEXUS = "nexus"


_EXT_MAP: dict[str, AlignmentFormat] = {
    ".fa":    AlignmentFormat.FASTA,
    ".fasta": AlignmentFormat.FASTA,
    ".faa":   AlignmentFormat.FASTA,
    ".fna":   AlignmentFormat.FASTA,
    ".phy":   AlignmentFormat.PHYLIP,
    ".nex":   AlignmentFormat.NEXUS,
    ".nxs":   AlignmentFormat.NEXUS,
    ".nexus": AlignmentFormat.NEXUS,
}


class FormatConverter:
    def detect(self, path: Path) -> AlignmentFormat:
        suffix = path.suffix.lower()
        if suffix in _EXT_MAP:
            return _EXT_MAP[suffix]
        if path.exists():
            first = path.read_text(errors="ignore")[:200]
            if first.strip().startswith(">"):
                return AlignmentFormat.FASTA
            if first.strip().upper().startswith("#NEXUS"):
                return AlignmentFormat.NEXUS
            lines = [l for l in first.splitlines() if l.strip()]
            if lines and all(p.isdigit() for p in lines[0].strip().split()[:2]):
                return AlignmentFormat.PHYLIP
        raise ValueError(
            f"Cannot detect alignment format for '{path}'. "
            f"Supported extensions: {list(_EXT_MAP.keys())}"
        )

    def convert(
        self,
        src: Path,
        dst: Path,
        target: AlignmentFormat,
        source_format: Optional[AlignmentFormat] = None,
        molecule_type: str = "protein",
    ) -> Path:
        src_fmt = source_format or self.detect(src)
        alignment = AlignIO.read(str(src), src_fmt.value)
        dst.parent.mkdir(parents=True, exist_ok=True)
        # Nexus writer requires molecule_type annotation on each record
        if target == AlignmentFormat.NEXUS:
            for record in alignment:
                if "molecule_type" not in record.annotations:
                    record.annotations["molecule_type"] = molecule_type
        with open(dst, "w") as fh:
            AlignIO.write(alignment, fh, target.value)
        return dst

    def read(self, path: Path) -> MultipleSeqAlignment:
        fmt = self.detect(path)
        return AlignIO.read(str(path), fmt.value)
