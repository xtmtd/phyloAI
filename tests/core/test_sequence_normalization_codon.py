from __future__ import annotations

from phyloai.core.sequence_normalization import validate_codon_msa


def test_valid_codon_msa_no_issues():
    seqs = {
        "seq1": "ATGGCTTCT",
        "seq2": "ATGGCTTCT",
    }
    result = validate_codon_msa(seqs)
    assert result.skip is False
    assert result.warnings == []
    assert result.sequences["seq1"] == "ATGGCTTCT"


def test_alignment_length_not_divisible_by_3_triggers_skip():
    seqs = {
        "seq1": "ATGGCT",
        "seq2": "ATGGCT",
    }
    result = validate_codon_msa(seqs)
    assert result.skip is False

    seqs_bad = {
        "seq1": "ATGGCTA",
        "seq2": "ATGGCTA",
    }
    result_bad = validate_codon_msa(seqs_bad)
    assert result_bad.skip is True
    assert any("multiple of 3" in warning for warning in result_bad.warnings)


def test_internal_stop_codon_warns_but_continues():
    seqs = {
        "seq1": "ATGTAAGCT",
        "seq2": "ATGGCTTCT",
    }
    result = validate_codon_msa(seqs)
    assert result.skip is False
    assert any("internal stop" in warning for warning in result.warnings)


def test_terminal_stop_codon_is_removed():
    seqs = {
        "seq1": "ATGGCTTAA",
        "seq2": "ATGGCTTAA",
    }
    result = validate_codon_msa(seqs)
    assert result.skip is False
    assert result.sequences["seq1"] == "ATGGCT"
    assert result.sequences["seq2"] == "ATGGCT"


def test_terminal_stop_in_some_taxa_preserves_equal_alignment_length():
    seqs = {
        "seq1": "ATGGCTTAA",
        "seq2": "ATGGCTTCT",
    }
    result = validate_codon_msa(seqs)
    assert result.skip is False
    assert len({len(seq) for seq in result.sequences.values()}) == 1
    assert "TAA" not in result.sequences["seq1"]


def test_terminal_stop_with_trailing_gaps_preserved():
    seqs = {
        "seq1": "ATGGCT---TAA",
        "seq2": "ATGGCTTCT---",
    }
    result = validate_codon_msa(seqs)
    assert result.skip is False
    assert len(result.sequences["seq1"].replace("-", "")) == 6
    assert result.sequences["seq2"].replace("-", "") == "ATGGCTTCT"
    assert len({len(seq) for seq in result.sequences.values()}) == 1


def test_gap_columns_in_msa_handled():
    seqs = {
        "seq1": "ATG---GCT",
        "seq2": "ATGGCTTCT",
    }
    result = validate_codon_msa(seqs)
    assert result.skip is False


def test_all_gapped_sequence_does_not_raise():
    seqs = {
        "seq1": "---------",
        "seq2": "ATGGCTTCT",
    }
    result = validate_codon_msa(seqs)
    assert result.skip is False


def test_empty_seqs_dict_returns_skip():
    result = validate_codon_msa({})
    assert result.skip is True
