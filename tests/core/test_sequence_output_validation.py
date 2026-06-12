from __future__ import annotations

from pathlib import Path


def test_validate_fasta_output_rejects_empty_file(tmp_path: Path) -> None:
    from phyloai.core.sequence_output_validation import validate_fasta_output

    out = tmp_path / "empty.fa"
    out.write_text("")

    result = validate_fasta_output(out, require_aligned=True)

    assert not result.ok
    assert result.n_records == 0
    assert any("empty" in w.lower() for w in result.warnings)


def test_validate_fasta_output_rejects_unequal_lengths_for_msa(tmp_path: Path) -> None:
    from phyloai.core.sequence_output_validation import validate_fasta_output

    out = tmp_path / "bad.fa"
    out.write_text(">a\nMKT\n>b\nMKTA\n")

    result = validate_fasta_output(out, require_aligned=True)

    assert not result.ok
    assert any("unequal" in w.lower() for w in result.warnings)


def test_validate_fasta_output_allows_unequal_lengths_for_sequences(tmp_path: Path) -> None:
    from phyloai.core.sequence_output_validation import validate_fasta_output

    out = tmp_path / "seqs.fa"
    out.write_text(">a\nMKT\n>b\nMKTA\n")

    result = validate_fasta_output(out, require_aligned=False)

    assert result.ok
    assert result.n_records == 2
    assert result.length == 3
