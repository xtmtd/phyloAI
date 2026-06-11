from __future__ import annotations

from phyloai.core.sequence_normalization import (
    classify_char,
    detect_seq_type,
    expand_dots_from_first_sequence,
    normalize_pattern_char,
    normalize_sequences,
)


def test_normalize_nt_preserves_iupac_and_converts_u_question_and_invalid() -> None:
    result = normalize_sequences(["acguryswkmbdhvn?.!"], seq_type="NT")

    assert result.sequences == ["ACGTRYSWKMBDHVNNNN"]
    assert result.seq_type == "NT"
    assert result.replacements["u_to_t"] == 1
    assert result.replacements["question_to_missing"] == 1
    assert result.replacements["invalid_to_missing"] == 2


def test_normalize_aa_converts_special_question_stop_and_invalid_to_x_by_default() -> None:
    result = normalize_sequences(["arndbzjxuop?*1"], seq_type="AA")

    assert result.sequences == ["ARNDXXXXXXPXXX"]
    assert result.replacements["aa_special_to_x"] == 6
    assert result.replacements["question_to_missing"] == 1
    assert result.replacements["stop_to_x"] == 1
    assert result.replacements["invalid_to_missing"] == 1


def test_normalize_aa_special_keep_preserves_bzjxuO() -> None:
    result = normalize_sequences(["BZJXUO"], seq_type="AA", aa_special="keep")

    assert result.sequences == ["BZJXUO"]
    assert result.replacements == {}


def test_expand_dots_from_first_sequence_when_lengths_match() -> None:
    expanded, counts = expand_dots_from_first_sequence(["ACGT", "A..T", "...."], missing_char="N")

    assert expanded == ["ACGT", "ACGT", "ACGT"]
    assert counts == {"dot_expanded": 6}


def test_expand_dots_uses_missing_char_when_lengths_do_not_match() -> None:
    expanded, counts = expand_dots_from_first_sequence(["ACGT", "A.."], missing_char="N")

    assert expanded == ["ACGT", "ANN"]
    assert counts == {"dot_to_missing": 2}


def test_shared_detect_seq_type_nt_with_iupac() -> None:
    assert detect_seq_type(["ACGT", "ACGR", "NNNN"]) == "NT"


def test_shared_detect_seq_type_defaults_to_aa_when_x_seen() -> None:
    assert detect_seq_type(["ACGTX"]) == "AA"


def test_shared_classify_char_matches_stats_terms() -> None:
    assert classify_char("A", "AA") == "standard"
    assert classify_char("B", "AA") == "ambiguous"
    assert classify_char("-", "AA") == "gap"
    assert classify_char("N", "NT") == "gap"
    assert classify_char("R", "NT") == "ambiguous"


def test_shared_normalize_pattern_char_treats_question_mark_as_gap() -> None:
    assert normalize_pattern_char("?") == "-"
    assert normalize_pattern_char("a") == "A"
