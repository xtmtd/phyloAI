from __future__ import annotations

from pathlib import Path

import pytest


def test_convert_single_file_defaults_to_fasta_and_normalizes_nt(tmp_path: Path) -> None:
    from phyloai.pretree.convert import convert_input

    src = tmp_path / "gene.fna"
    src.write_text(">tax one\nacgu?ry!\n")
    out_dir = tmp_path / "runs" / "run001" / "pretree" / "convert"

    payload = convert_input(src, out_dir, target_format="fasta", seq_type="NT", threads=1, overwrite=False)

    out = out_dir / "seqs" / "gene.fa"
    assert out.exists()
    assert out.read_text() == ">tax_one\nACGTNRYN\n"
    assert payload["data"]["summary"]["n_converted"] == 1
    assert payload["data"]["summary"]["total_replacements"] == 3
    assert payload["data"]["files"][0]["output"] == str(out)


def test_convert_directory_skips_non_sequence_empty_and_subdirectory(tmp_path: Path) -> None:
    from phyloai.pretree.convert import convert_input

    src_dir = tmp_path / "raw"
    src_dir.mkdir()
    (src_dir / "ok.fa").write_text(">a\nACGT\n")
    (src_dir / "notes.txt").write_text("not sequence")
    (src_dir / "empty.fa").write_text("")
    (src_dir / "nested").mkdir()
    out_dir = tmp_path / "out"

    payload = convert_input(src_dir, out_dir, target_format="fasta", seq_type="NT", threads=1, overwrite=False)

    assert payload["data"]["summary"]["n_converted"] == 1
    assert payload["data"]["summary"]["n_skipped"] == 3
    reasons = {item["reason"] for item in payload["data"]["skipped"]}
    assert "empty file" in reasons
    assert "directory" in reasons


def test_convert_expands_paml_dots_before_normalization(tmp_path: Path) -> None:
    from phyloai.pretree.convert import convert_input

    src = tmp_path / "dots.fa"
    src.write_text(">ref\nACGT\n>second\nA..T\n")
    out_dir = tmp_path / "out"

    payload = convert_input(src, out_dir, target_format="fasta", seq_type="NT", threads=1, overwrite=False)

    assert (out_dir / "seqs" / "dots.fa").read_text() == ">ref\nACGT\n>second\nACGT\n"
    assert payload["data"]["files"][0]["replacements"]["dot_expanded"] == 2


def test_convert_output_dir_conflict_requires_overwrite(tmp_path: Path) -> None:
    from phyloai.pretree.convert import convert_input

    src = tmp_path / "gene.fa"
    src.write_text(">a\nACGT\n")
    out_dir = tmp_path / "out"
    out_dir.mkdir()
    (out_dir / "old.fa").write_text("old")

    try:
        convert_input(src, out_dir, target_format="fasta", seq_type="NT", threads=1, overwrite=False)
    except ValueError as exc:
        assert "already exists and is non-empty" in str(exc)
    else:
        raise AssertionError("Expected output directory conflict")

    payload = convert_input(src, out_dir, target_format="fasta", seq_type="NT", threads=1, overwrite=True)
    assert payload["data"]["summary"]["n_converted"] == 1
    assert not (out_dir / "old.fa").exists()


def test_convert_phylip_relaxed(tmp_path: Path) -> None:
    from phyloai.pretree.convert import convert_input

    src = tmp_path / "gene.fa"
    src.write_text(">taxon1\nACGT\n>taxon2\nACGA\n")
    out_dir = tmp_path / "out"

    payload = convert_input(src, out_dir, target_format="phylip-relaxed", seq_type="NT", threads=1, overwrite=False)

    out = out_dir / "seqs" / "gene.phy"
    assert out.exists()
    content = out.read_text()
    lines = content.strip().split("\n")
    assert lines[0] == "2 4"
    assert "taxon1" in lines[1]
    assert "ACGT" in lines[1]
    assert "taxon2" in lines[2]
    assert "ACGA" in lines[2]
    assert payload["data"]["summary"]["n_converted"] == 1


def test_convert_reads_phylip_paml_input_and_expands_dots(tmp_path: Path) -> None:
    from phyloai.pretree.convert import convert_input

    src = tmp_path / "gene.paml.phy"
    src.write_text("2 4\ntaxon1  ACGT\ntaxon2  A..T\n")
    out_dir = tmp_path / "out"

    payload = convert_input(src, out_dir, target_format="fasta", seq_type="NT", threads=1, overwrite=False)

    out = out_dir / "seqs" / "gene.paml.fa"
    assert out.read_text() == ">taxon1\nACGT\n>taxon2\nACGT\n"
    assert payload["data"]["summary"]["n_converted"] == 1
    assert payload["data"]["files"][0]["input_format"] == "phylip-paml"
    assert payload["data"]["files"][0]["replacements"]["dot_expanded"] == 2


def test_convert_reports_paml_specific_taxon_name_changes(tmp_path: Path) -> None:
    from phyloai.pretree.convert import convert_input

    src = tmp_path / "gene.fa"
    src.write_text(
        ">Taxon name with spaces and very long suffix\nACGT\n"
        ">Taxon:bad#chars\nACGA\n"
    )
    out_dir = tmp_path / "out"

    payload = convert_input(src, out_dir, target_format="phylip-paml", seq_type="NT", threads=1, overwrite=False)

    file_entry = payload["data"]["files"][0]
    paml_lines = (out_dir / "seqs" / "gene.paml.phy").read_text().splitlines()
    assert paml_lines[0] == "2  4  S"
    assert paml_lines[1][30:32] == "  "
    assert file_entry["taxon_name_changes"] == 2
    assert payload["data"]["summary"]["total_taxon_name_changes"] == 2
    assert len(file_entry["taxon_name_change_details"]) == 2


def test_convert_aa_special_keep_preserves_special_codes(tmp_path: Path) -> None:
    from phyloai.pretree.convert import convert_input

    src = tmp_path / "protein.fa"
    src.write_text(">taxon1\nBZJXUO?*1\n")
    out_dir = tmp_path / "out"

    payload = convert_input(src, out_dir, target_format="fasta", seq_type="AA", aa_special="keep", threads=1, overwrite=False)

    assert (out_dir / "seqs" / "protein.fa").read_text() == ">taxon1\nBZJXUOXXX\n"
    replacements = payload["data"]["files"][0]["replacements"]
    assert "aa_special_to_x" not in replacements
    assert replacements["question_to_missing"] == 1
    assert replacements["stop_to_x"] == 1
    assert replacements["invalid_to_missing"] == 1


def test_convert_skips_invalid_generated_fasta(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    from phyloai.pretree import convert

    src = tmp_path / "gene.fa"
    src.write_text(">taxon1\nACGT\n")
    out_dir = tmp_path / "out"

    def write_empty(_records, out: Path, _target_format: str, _seq_type: str):
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text("")
        return []

    monkeypatch.setattr(convert, "_write_records", write_empty)

    with pytest.raises(ValueError, match="All input entries failed"):
        convert.convert_input(src, out_dir, target_format="fasta", seq_type="NT", threads=1, overwrite=False)
