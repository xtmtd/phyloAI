"""Concatenate multiple MSAs into a supermatrix for phylogenetic inference."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from phyloai.core.formats import AlignmentFormat, FormatConverter
from phyloai.core.schema import COMMON_ALIGNMENT_EXTENSIONS
from phyloai.core.sequence_normalization import detect_seq_type, normalize_sequences

AA_RECODING_TABLES: dict[str, dict[str, str]] = {
    "Dayhoff-6": {
        "A": "0", "G": "0", "P": "0", "S": "0", "T": "0",
        "D": "1", "E": "1", "N": "1", "Q": "1",
        "H": "2", "K": "2", "R": "2",
        "I": "3", "L": "3", "M": "3", "V": "3",
        "F": "4", "W": "4", "Y": "4",
        "C": "5",
        "X": "X",
    },
    "Dayhoff-9": {
        "D": "0", "E": "0", "H": "0", "N": "0", "Q": "0",
        "I": "1", "L": "1", "M": "1", "V": "1",
        "F": "2", "Y": "2",
        "A": "3", "S": "3", "T": "3",
        "K": "4", "R": "4",
        "G": "5",
        "P": "6",
        "C": "7",
        "W": "8",
        "X": "X",
    },
    "Dayhoff-12": {
        "D": "0", "E": "0", "Q": "0",
        "M": "1", "L": "1", "I": "1", "V": "1",
        "F": "2", "Y": "2",
        "K": "3", "H": "3", "R": "3",
        "G": "4",
        "A": "5",
        "P": "6",
        "S": "7",
        "T": "8",
        "N": "9",
        "W": "A",
        "C": "B",
        "X": "X",
    },
    "Dayhoff-15": {
        "D": "0", "E": "0", "Q": "0",
        "M": "1", "L": "1",
        "I": "2", "V": "2",
        "F": "3", "Y": "3",
        "G": "4",
        "A": "5",
        "P": "6",
        "S": "7",
        "T": "8",
        "N": "9",
        "K": "A",
        "H": "B",
        "R": "C",
        "W": "D",
        "C": "E",
        "X": "X",
    },
    "Dayhoff-18": {
        "F": "0", "Y": "0",
        "M": "1", "L": "1",
        "I": "2",
        "V": "3",
        "G": "4",
        "A": "5",
        "P": "6",
        "S": "7",
        "T": "8",
        "D": "9",
        "E": "A",
        "Q": "B",
        "N": "C",
        "K": "D",
        "H": "E",
        "R": "F",
        "W": "G",
        "C": "H",
        "X": "X",
    },
    "SandR-6": {
        "A": "0", "P": "0", "S": "0", "T": "0",
        "D": "1", "E": "1", "N": "1", "G": "1",
        "Q": "2", "K": "2", "R": "2",
        "M": "3", "I": "3", "V": "3", "L": "3",
        "W": "4", "C": "4",
        "F": "5", "Y": "5", "H": "5",
        "X": "X",
    },
    "KGB-6": {
        "A": "0", "G": "0", "P": "0", "S": "0",
        "D": "1", "E": "1", "N": "1", "Q": "1", "H": "1", "K": "1", "R": "1", "T": "1",
        "M": "2", "I": "2", "L": "2",
        "W": "3",
        "F": "4", "Y": "4",
        "C": "5", "V": "5",
        "X": "X",
    },
}

NT_RECODING_TABLES: dict[str, dict[str, str]] = {
    "RY-nucleotide": {
        "A": "R", "G": "R",
        "C": "Y", "T": "Y", "U": "Y",
        "N": "?",
        "X": "?",
        "R": "R", "Y": "Y",
        "S": "?", "W": "?", "K": "?", "M": "?",
        "B": "?", "D": "?", "H": "?", "V": "?",
        "-": "-",
        "?": "?",
        ".": ".",
    },
}

_GAP_CHARS = frozenset("-?.*")


def _apply_recoding(
    matrix: dict[str, str], scheme: str
) -> tuple[dict[str, str], list[str]]:
    if scheme in AA_RECODING_TABLES:
        table = AA_RECODING_TABLES[scheme]
    elif scheme in NT_RECODING_TABLES:
        table = NT_RECODING_TABLES[scheme]
    else:
        raise ValueError(f"Unknown recoding scheme: {scheme!r}")

    warnings_set: set[str] = set()
    result: dict[str, str] = {}
    for taxon, seq in matrix.items():
        chars: list[str] = []
        for ch in seq:
            if ch in _GAP_CHARS:
                chars.append(ch)
            elif ch in table:
                chars.append(table[ch])
            else:
                chars.append(ch)
                warnings_set.add(
                    f"Character '{ch}' not in recoding table '{scheme}', passed through unchanged"
                )
        result[taxon] = "".join(chars)
    return result, sorted(warnings_set)


def _translate_codon(seq: str) -> str:
    n_complete = (len(seq) // 3) * 3
    result: list[str] = []
    for i in range(0, n_complete, 3):
        codon = seq[i:i + 3]
        if "-" in codon:
            result.append("-")
        else:
            result.append(str(Seq(codon).translate()))
    return "".join(result)


def _exclude_codon3(seq: str) -> str:
    return "".join(ch for i, ch in enumerate(seq) if i % 3 != 2)


def _scan_msa_files(msa_dir: Path) -> list[Path]:
    found = []
    for ext in COMMON_ALIGNMENT_EXTENSIONS:
        found.extend(msa_dir.glob(f"*{ext}"))
    return sorted(set(found))


def _read_msa(path: Path) -> tuple[list[str], list[str], int]:
    converter = FormatConverter()
    fmt = converter.detect(path)
    alignment = converter.read(path, source_format=fmt)
    taxa = [record.id for record in alignment]
    seqs = [str(record.seq).upper() for record in alignment]
    length = alignment.get_alignment_length()
    return taxa, seqs, length


def _read_msa_headers(path: Path) -> list[str]:
    suffix = path.suffix.lower()
    if suffix in (".fa", ".faa", ".ffn", ".frn", ".fasta", ".fas"):
        ids = []
        with open(path) as fh:
            for line in fh:
                line = line.strip()
                if line.startswith(">"):
                    ids.append(line[1:].split(None, 1)[0])
        return ids
    else:
        converter = FormatConverter()
        fmt = converter.detect(path)
        alignment = converter.read(path, source_format=fmt)
        return [record.id for record in alignment]


def _filter_by_occupancy(
    msa_paths: list[Path],
    msa_taxa: dict[str, set[str]],
    total_taxa: set[str],
    threshold: float,
) -> tuple[list[Path], list[dict[str, Any]]]:
    n_total = len(total_taxa)
    kept: list[Path] = []
    dropped: list[dict[str, Any]] = []
    for path in msa_paths:
        taxa = msa_taxa[str(path)]
        n_taxa = len(taxa)
        ratio = n_taxa / n_total if n_total > 0 else 0.0
        if ratio >= threshold:
            kept.append(path)
        else:
            dropped.append({
                "filename": path.name,
                "n_taxa": n_taxa,
                "occupancy_ratio": round(ratio, 4),
                "total_taxa": n_total,
            })
    return kept, dropped
