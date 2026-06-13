from __future__ import annotations

from pathlib import Path

import pytest


def test_dayhoff6_recoding_maps_all_standard_amino_acids() -> None:
    from phyloai.pretree.concat import AA_RECODING_TABLES

    table = AA_RECODING_TABLES["Dayhoff-6"]
    assert table["A"] == "0"
    assert table["G"] == "0"
    assert table["P"] == "0"
    assert table["S"] == "0"
    assert table["T"] == "0"
    assert table["D"] == "1"
    assert table["E"] == "1"
    assert table["N"] == "1"
    assert table["Q"] == "1"
    assert table["H"] == "2"
    assert table["K"] == "2"
    assert table["R"] == "2"
    assert table["I"] == "3"
    assert table["L"] == "3"
    assert table["M"] == "3"
    assert table["V"] == "3"
    assert table["F"] == "4"
    assert table["W"] == "4"
    assert table["Y"] == "4"
    assert table["C"] == "5"


def test_ry_nucleotide_recoding_includes_ambiguous() -> None:
    from phyloai.pretree.concat import NT_RECODING_TABLES

    table = NT_RECODING_TABLES["RY-nucleotide"]
    assert table["A"] == "R"
    assert table["G"] == "R"
    assert table["C"] == "Y"
    assert table["T"] == "Y"
    assert table["U"] == "Y"
    assert table["N"] == "?"
    assert table["X"] == "?"
    assert table["-"] == "-"
    assert table["?"] == "?"
    assert table["."] == "."


def test_apply_recoding_dayhoff6() -> None:
    from phyloai.pretree.concat import _apply_recoding

    matrix = {"tax1": "ACDEFG", "tax2": "HIKLMN"}
    recoded, warnings = _apply_recoding(matrix, "Dayhoff-6")
    assert recoded["tax1"] == "051140"
    assert recoded["tax2"] == "232331"
    assert warnings == []


def test_apply_recoding_preserves_gaps_and_special_chars() -> None:
    from phyloai.pretree.concat import _apply_recoding

    matrix = {"tax1": "A-H.?*X"}
    recoded, warnings = _apply_recoding(matrix, "Dayhoff-6")
    assert recoded["tax1"] == "0-2.?*X"
    assert warnings == []


def test_apply_recoding_unknown_scheme_raises() -> None:
    from phyloai.pretree.concat import _apply_recoding

    with pytest.raises(ValueError, match="Unknown recoding scheme"):
        _apply_recoding({"tax1": "ACG"}, "FakeScheme")


def test_apply_recoding_ry_nucleotide() -> None:
    from phyloai.pretree.concat import _apply_recoding

    matrix = {"tax1": "ACGTN--?."}
    recoded, warnings = _apply_recoding(matrix, "RY-nucleotide")
    assert recoded["tax1"] == "RYRY?--?."


def test_translate_codon_standard_genetic_code() -> None:
    from phyloai.pretree.concat import _translate_codon

    assert _translate_codon("ATGCGTAAA") == "MRK"
    assert _translate_codon("TTTGGGCCC") == "FGP"


def test_translate_codon_with_gaps_preserves_codon_structure() -> None:
    from phyloai.pretree.concat import _translate_codon

    assert _translate_codon("ATG---AAA") == "M-K"
    assert _translate_codon("---ATGCGT") == "-MR"
    assert _translate_codon("ATG-AA-TAA") == "M--"


def test_translate_codon_trims_incomplete_codon_at_end() -> None:
    from phyloai.pretree.concat import _translate_codon

    assert _translate_codon("ATGCG") == "M"


def test_translate_codon_all_gaps_then_bases() -> None:
    from phyloai.pretree.concat import _translate_codon

    assert _translate_codon("---ATGAAA") == "-MK"


def test_exclude_codon3_drops_every_third_position() -> None:
    from phyloai.pretree.concat import _exclude_codon3

    assert _exclude_codon3("ATGCGTAAATTT") == "ATCGAATT"


def test_exclude_codon3_preserves_gaps_at_kept_positions() -> None:
    from phyloai.pretree.concat import _exclude_codon3

    assert _exclude_codon3("A-GC-T---") == "A-C---"


def test_scan_msa_files_finds_supported_extensions(tmp_path: Path) -> None:
    from phyloai.pretree.concat import _scan_msa_files

    (tmp_path / "gene1.fa").write_text(">a\nACGT\n")
    (tmp_path / "gene2.fasta").write_text(">a\nACGT\n")
    (tmp_path / "gene3.phy").write_text("2 4\na ACGT\nb ACGT\n")
    (tmp_path / "notes.txt").write_text("not an alignment")
    (tmp_path / "subdir").mkdir()

    found = _scan_msa_files(tmp_path)
    names = sorted(p.name for p in found)
    assert names == ["gene1.fa", "gene2.fasta", "gene3.phy"]


def test_scan_msa_files_empty_dir_returns_empty(tmp_path: Path) -> None:
    from phyloai.pretree.concat import _scan_msa_files

    assert _scan_msa_files(tmp_path) == []


def test_read_msa_fasta(tmp_path: Path) -> None:
    from phyloai.pretree.concat import _read_msa

    msa_path = tmp_path / "gene.fa"
    msa_path.write_text(">tax1\nACGT\n>tax2\nACGT\n")

    taxa, seqs, length = _read_msa(msa_path)
    assert taxa == ["tax1", "tax2"]
    assert seqs == ["ACGT", "ACGT"]
    assert length == 4


def test_read_msa_phylip_paml(tmp_path: Path) -> None:
    from phyloai.pretree.concat import _read_msa

    msa_path = tmp_path / "gene.phy"
    msa_path.write_text("2 4\ntax1  ACGT\ntax2  ACGT\n")

    taxa, seqs, length = _read_msa(msa_path)
    assert taxa == ["tax1", "tax2"]
    assert seqs == ["ACGT", "ACGT"]
    assert length == 4


def test_read_msa_headers_fasta(tmp_path: Path) -> None:
    from phyloai.pretree.concat import _read_msa_headers

    msa_path = tmp_path / "gene.fa"
    msa_path.write_text(">tax1\nACGT\n>tax2\nGGCC\n")

    taxa = _read_msa_headers(msa_path)
    assert taxa == ["tax1", "tax2"]


def test_read_msa_headers_phylip(tmp_path: Path) -> None:
    from phyloai.pretree.concat import _read_msa_headers

    msa_path = tmp_path / "gene.phy"
    msa_path.write_text("2 4\ntax1  ACGT\ntax2  ACGT\n")

    taxa = _read_msa_headers(msa_path)
    assert taxa == ["tax1", "tax2"]


def test_read_msa_headers_is_fast_and_does_not_parse_sequences(tmp_path: Path) -> None:
    from phyloai.pretree.concat import _read_msa_headers

    msa_path = tmp_path / "big.fa"
    msa_path.write_text(">tax1\n" + ("A" * 10000 + "\n") * 100)

    taxa = _read_msa_headers(msa_path)
    assert taxa == ["tax1"]
