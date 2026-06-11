"""Shared sequence character normalization for PhyloAI commands."""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass, field


AA_STANDARD = set("ARNDCQEGHILKMFPSTWYV")
AA_SPECIAL = set("BZJXUO")
NT_STANDARD = set("ACGT")
NT_AMBIGUOUS = set("RYSWKMBDHVN")
NT_AUTO_CHARS = set("ACGTURYSWKMBDHVN")
GAP_CHARS = set("-?")
NT_MISSING = set("-?N")


@dataclass
class NormalizationResult:
    sequences: list[str]
    seq_type: str
    replacements: dict[str, int] = field(default_factory=dict)
    warnings: list[str] = field(default_factory=list)


def normalize_sequences(
    sequences: list[str],
    seq_type: str,
    aa_special: str = "x",
) -> NormalizationResult:
    counters: Counter[str] = Counter()
    normalized = [_normalize_one(sequence, seq_type, aa_special, counters) for sequence in sequences]
    return NormalizationResult(sequences=normalized, seq_type=seq_type, replacements=dict(counters), warnings=[])


def _normalize_one(sequence: str, seq_type: str, aa_special: str, counters: Counter[str]) -> str:
    output: list[str] = []
    for raw_char in sequence:
        if raw_char.isspace():
            counters["whitespace_removed"] += 1
            continue
        char = raw_char.upper()
        if seq_type == "NT":
            output.append(_normalize_nt_char(char, counters))
        else:
            output.append(_normalize_aa_char(char, aa_special, counters))
    return "".join(output)


def _normalize_nt_char(char: str, counters: Counter[str]) -> str:
    if char in NT_STANDARD or char in NT_AMBIGUOUS or char == "-":
        return char
    if char == "U":
        counters["u_to_t"] += 1
        return "T"
    if char == "?":
        counters["question_to_missing"] += 1
        return "N"
    counters["invalid_to_missing"] += 1
    return "N"


def _normalize_aa_char(char: str, aa_special: str, counters: Counter[str]) -> str:
    if char in AA_STANDARD or char == "-":
        return char
    if char in AA_SPECIAL:
        if aa_special == "keep":
            return char
        counters["aa_special_to_x"] += 1
        return "X"
    if char == "?":
        counters["question_to_missing"] += 1
        return "X"
    if char == "*":
        counters["stop_to_x"] += 1
        return "X"
    counters["invalid_to_missing"] += 1
    return "X"


def expand_dots_from_first_sequence(sequences: list[str], missing_char: str) -> tuple[list[str], dict[str, int]]:
    counters: Counter[str] = Counter()
    if not sequences:
        return [], {}
    reference = sequences[0]
    same_length = all(len(sequence) == len(reference) for sequence in sequences)
    expanded: list[str] = []
    for sequence in sequences:
        chars: list[str] = []
        for index, char in enumerate(sequence):
            if char != ".":
                chars.append(char)
                continue
            if same_length and index < len(reference) and reference[index] != ".":
                chars.append(reference[index])
                counters["dot_expanded"] += 1
            else:
                chars.append(missing_char)
                counters["dot_to_missing"] += 1
        expanded.append("".join(chars))
    return expanded, dict(counters)


def detect_seq_type(sequences: list[str]) -> str:
    seq_type, _warnings = resolve_seq_type(sequences)
    return seq_type


def resolve_seq_type(sequences: list[str]) -> tuple[str, list[str]]:
    observed: set[str] = set()
    for sequence in sequences:
        observed.update(sequence.upper())
    observed.difference_update(GAP_CHARS)
    observed.discard(".")
    if observed & set("EFILPQWYZ"):
        return "AA", []
    if observed and observed <= NT_AUTO_CHARS:
        return "NT", []
    return "AA", ["[WARN] Cannot determine seq_type, defaulting to AA"]


def normalize_pattern_char(char: str) -> str:
    upper = char.upper()
    return "-" if upper == "?" else upper


def standard_chars(seq_type: str) -> set[str]:
    return NT_STANDARD if seq_type == "NT" else AA_STANDARD


def gap_chars(seq_type: str) -> set[str]:
    return NT_MISSING if seq_type == "NT" else GAP_CHARS


def classify_char(char: str, seq_type: str) -> str:
    upper = char.upper()
    if upper in standard_chars(seq_type):
        return "standard"
    if upper in gap_chars(seq_type):
        return "gap"
    return "ambiguous"
