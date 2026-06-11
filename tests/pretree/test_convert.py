from __future__ import annotations

from pathlib import Path


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
