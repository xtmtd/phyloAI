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
    assert recoded["tax1"] == "0-2.?*?"
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
    (tmp_path / "not_a_file.fa").mkdir()

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


def test_filter_by_occupancy_keeps_msas_at_or_above_threshold(tmp_path: Path) -> None:
    from phyloai.pretree.concat import _filter_by_occupancy

    msa_paths = [tmp_path / "gene1.fa", tmp_path / "gene2.fa", tmp_path / "gene3.fa"]
    msa_taxa = {
        str(msa_paths[0]): {"A", "B", "C", "D"},
        str(msa_paths[1]): {"A", "B", "C", "D"},
        str(msa_paths[2]): {"A", "B"},
    }
    total_taxa = {"A", "B", "C", "D"}

    kept, dropped = _filter_by_occupancy(msa_paths, msa_taxa, total_taxa, 0.5)
    assert len(kept) == 3
    assert len(dropped) == 0


def test_filter_by_occupancy_drops_msas_below_threshold(tmp_path: Path) -> None:
    from phyloai.pretree.concat import _filter_by_occupancy

    msa_paths = [tmp_path / "gene1.fa", tmp_path / "gene2.fa", tmp_path / "gene3.fa"]
    msa_taxa = {
        str(msa_paths[0]): {"A", "B", "C", "D"},
        str(msa_paths[1]): {"A", "B", "C"},
        str(msa_paths[2]): {"A"},
    }
    total_taxa = {"A", "B", "C", "D"}

    kept, dropped = _filter_by_occupancy(msa_paths, msa_taxa, total_taxa, 0.5)
    assert len(kept) == 2
    assert msa_paths[0] in kept
    assert msa_paths[1] in kept
    assert len(dropped) == 1
    assert dropped[0]["filename"] == "gene3.fa"
    assert dropped[0]["n_taxa"] == 1
    assert dropped[0]["occupancy_ratio"] == 0.25


def test_filter_by_occupancy_zero_keeps_all(tmp_path: Path) -> None:
    from phyloai.pretree.concat import _filter_by_occupancy

    msa_paths = [tmp_path / "gene1.fa"]
    msa_taxa = {str(msa_paths[0]): {"A"}}
    total_taxa = {"A", "B", "C", "D"}

    kept, dropped = _filter_by_occupancy(msa_paths, msa_taxa, total_taxa, 0.0)
    assert len(kept) == 1
    assert len(dropped) == 0


def test_filter_by_occupancy_one_keeps_only_full(tmp_path: Path) -> None:
    from phyloai.pretree.concat import _filter_by_occupancy

    msa_paths = [tmp_path / "gene1.fa", tmp_path / "gene2.fa"]
    msa_taxa = {
        str(msa_paths[0]): {"A", "B", "C"},
        str(msa_paths[1]): {"A", "B"},
    }
    total_taxa = {"A", "B", "C"}

    kept, dropped = _filter_by_occupancy(msa_paths, msa_taxa, total_taxa, 1.0)
    assert len(kept) == 1
    assert msa_paths[0] in kept
    assert len(dropped) == 1


def test_concat_alignments_builds_supermatrix(tmp_path: Path) -> None:
    from phyloai.pretree.concat import _concat_alignments

    msa_paths = [tmp_path / "gene1.fa", tmp_path / "gene2.fa"]
    msa_data = {
        str(msa_paths[0]): (["A", "B"], ["ACGT", "ACGT"], 4),
        str(msa_paths[1]): (["A", "C"], ["GGCC", "GGCC"], 4),
    }
    total_taxa = {"A", "B", "C"}

    matrix, taxon_order = _concat_alignments(msa_paths, msa_data, total_taxa)
    assert set(matrix.keys()) == {"A", "B", "C"}
    assert matrix["A"] == "ACGTGGCC"
    assert matrix["B"] == "ACGT????"
    assert matrix["C"] == "????GGCC"
    assert taxon_order == ["A", "B", "C"]


def test_reorder_outgroup_moves_taxon_to_first() -> None:
    from phyloai.pretree.concat import _reorder_outgroup

    matrix = {"A": "ACGT", "B": "ACGT", "C": "ACGT"}
    reordered = _reorder_outgroup(matrix, "B")
    assert list(reordered.keys())[0] == "B"


def test_reorder_outgroup_not_found_raises() -> None:
    from phyloai.pretree.concat import _reorder_outgroup

    matrix = {"A": "ACGT", "B": "ACGT"}
    with pytest.raises(ValueError, match="Outgroup taxon 'X' not found"):
        _reorder_outgroup(matrix, "X")


def test_reorder_outgroup_none_returns_unchanged() -> None:
    from phyloai.pretree.concat import _reorder_outgroup

    matrix = {"A": "ACGT", "B": "ACGT"}
    result = _reorder_outgroup(matrix, None)
    assert list(result.keys()) == ["A", "B"]


def test_write_matrix_fasta(tmp_path: Path) -> None:
    from phyloai.pretree.concat import _write_matrix

    matrix = {"tax1": "ACGT", "tax2": "ACGT"}
    out_path = tmp_path / "out.fa"
    _write_matrix(matrix, out_path, "fasta", "NT")

    content = out_path.read_text()
    assert ">tax1" in content
    assert "ACGT" in content


def test_write_matrix_phylip_relaxed(tmp_path: Path) -> None:
    from phyloai.pretree.concat import _write_matrix

    matrix = {"tax1": "ACGT", "tax2": "ACGT"}
    out_path = tmp_path / "out.phy"
    _write_matrix(matrix, out_path, "phylip-relaxed", "NT")

    content = out_path.read_text()
    assert "2 4" in content
    assert "tax1" in content


def test_compute_concat_stats_uses_per_taxon_stats() -> None:
    from phyloai.pretree.concat import _compute_concat_stats

    matrix = {"tax1": "ACGT", "tax2": "ACGT", "tax3": "ACGT"}
    stats = _compute_concat_stats(matrix, "NT")

    assert stats["n_taxa"] == 3
    assert stats["alignment_length"] == 4
    assert stats["character_summary"]["gap_ratio"] == 0.0
    assert stats["character_summary"]["ambiguous_ratio"] == 0.0
    assert stats["site_patterns"]["constant_sites"]["count"] == 4


def test_render_concat_panels_shows_all_overview_fields() -> None:
    from phyloai.pretree.concat import _render_concat_panels
    from rich.panel import Panel

    overview = {
        "prefix": "matrix",
        "to_format": "fasta",
        "n_taxa": 10,
        "n_msa_input": 50,
        "n_msa_used": 45,
        "n_msa_dropped": 5,
        "taxon_occupancy_threshold": 0.5,
        "recoding": "Dayhoff-6",
        "outgroup": "Sp_A",
        "variants_produced": ["matrix.fa", "matrix.recoded.fa"],
    }
    variant_stats = [
        {
            "variant": "original",
            "seq_type": "AA",
            "total_length": 100,
            "character_summary": {
                "gap_ratio": 0.1, "ambiguous_ratio": 0.02,
                "gap_ambiguous_ratio": 0.12, "standard_ratio": 0.88,
            },
            "site_patterns": {
                "alignment_length": 100,
                "distinct_patterns": {"count": 50, "ratio": 0.5},
                "constant_sites": {"count": 30, "ratio": 0.3},
                "parsimony_informative": {"count": 15, "ratio": 0.15},
                "singleton_sites": {"count": 5, "ratio": 0.05},
            },
        },
        {
            "variant": "recoded",
            "seq_type": "AA",
            "total_length": 100,
            "character_summary": {
                "gap_ratio": 0.1, "ambiguous_ratio": 0.02,
                "gap_ambiguous_ratio": 0.12, "standard_ratio": 0.88,
            },
            "site_patterns": {
                "alignment_length": 100,
                "distinct_patterns": {"count": 40, "ratio": 0.4},
                "constant_sites": {"count": 50, "ratio": 0.5},
                "parsimony_informative": {"count": 8, "ratio": 0.08},
                "singleton_sites": {"count": 2, "ratio": 0.02},
            },
        },
    ]
    panels = _render_concat_panels(overview, variant_stats)
    assert len(panels) == 3
    assert all(isinstance(p, Panel) for p in panels)


def test_render_concat_panels_hides_recoding_when_none() -> None:
    from phyloai.pretree.concat import _render_concat_panels

    overview = {
        "prefix": "matrix",
        "to_format": "fasta",
        "n_taxa": 5,
        "n_msa_input": 10,
        "n_msa_used": 10,
        "n_msa_dropped": 0,
        "taxon_occupancy_threshold": 0.5,
        "recoding": None,
        "outgroup": None,
        "variants_produced": ["matrix.fa"],
    }
    variant_stats = [
        {
            "variant": "original",
            "seq_type": "NT",
            "total_length": 200,
            "character_summary": {
                "gap_ratio": 0.0, "ambiguous_ratio": 0.0,
                "gap_ambiguous_ratio": 0.0, "standard_ratio": 1.0,
            },
            "site_patterns": {
                "alignment_length": 200,
                "distinct_patterns": {"count": 10, "ratio": 0.05},
                "constant_sites": {"count": 180, "ratio": 0.9},
                "parsimony_informative": {"count": 5, "ratio": 0.025},
                "singleton_sites": {"count": 5, "ratio": 0.025},
            },
        },
    ]
    panels = _render_concat_panels(overview, variant_stats)
    assert len(panels) == 3


def test_run_concat_basic(tmp_path: Path) -> None:
    from phyloai.pretree.concat import run_concat

    msa_dir = tmp_path / "msas"
    msa_dir.mkdir()
    (msa_dir / "gene1.fa").write_text(">A\nACGT\n>B\nACGT\n>C\nACGT\n")
    (msa_dir / "gene2.fa").write_text(">A\nGGCC\n>B\nGGCC\n>C\nGGCC\n")

    output_dir = tmp_path / "out"
    payload = run_concat(
        msa_dir=msa_dir,
        output_dir=output_dir,
        prefix="matrix",
        seq_type="NT",
        taxa_occupancy=0.5,
        recoding=None,
        outgroup=None,
        to_format="fasta",
        translate_codon=False,
        exclude_codon3=False,
        dry_run=False,
        overwrite=False,
    )

    assert payload["status"] == "success"
    assert payload["key_results"]["n_taxa"] == 3
    assert payload["key_results"]["n_msa_used"] == 2
    assert (output_dir / "matrix.fa").exists()
    content = (output_dir / "matrix.fa").read_text()
    assert ">A" in content
    assert "ACGTGGCC" in content
    assert (output_dir / "concat.log").exists()


def test_run_concat_with_recoding_and_warnings(tmp_path: Path) -> None:
    from phyloai.pretree.concat import run_concat

    msa_dir = tmp_path / "msas"
    msa_dir.mkdir()
    (msa_dir / "gene1.fa").write_text(">A\nACGT\n>B\nACGT\n")

    output_dir = tmp_path / "out"
    payload = run_concat(
        msa_dir=msa_dir, output_dir=output_dir, prefix="matrix",
        seq_type="NT", taxa_occupancy=0.0, recoding="RY-nucleotide",
        outgroup=None, to_format="fasta",
        translate_codon=False, exclude_codon3=False,
        dry_run=False, overwrite=False,
    )

    assert payload["status"] == "success"
    recoded_path = output_dir / "matrix.recoded.fa"
    assert recoded_path.exists()
    content = recoded_path.read_text()
    assert "RYRY" in content
    assert "recoding_warnings" in payload["data"]


def test_run_concat_occupancy_filtering(tmp_path: Path) -> None:
    from phyloai.pretree.concat import run_concat

    msa_dir = tmp_path / "msas"
    msa_dir.mkdir()
    (msa_dir / "gene1.fa").write_text(">A\nACGT\n>B\nACGT\n>C\nACGT\n")
    (msa_dir / "gene2.fa").write_text(">A\nGGCC\n")

    output_dir = tmp_path / "out"
    payload = run_concat(
        msa_dir=msa_dir, output_dir=output_dir, prefix="matrix",
        seq_type="NT", taxa_occupancy=0.5, recoding=None,
        outgroup=None, to_format="fasta",
        translate_codon=False, exclude_codon3=False,
        dry_run=False, overwrite=False,
    )

    assert payload["key_results"]["n_msa_used"] == 1
    assert payload["key_results"]["n_msa_dropped"] == 1
    assert (output_dir / "dropped_alignments.csv").exists()


def test_run_concat_outgroup_reordering(tmp_path: Path) -> None:
    from phyloai.pretree.concat import run_concat

    msa_dir = tmp_path / "msas"
    msa_dir.mkdir()
    (msa_dir / "gene1.fa").write_text(">A\nACGT\n>B\nACGT\n>C\nACGT\n")

    output_dir = tmp_path / "out"
    run_concat(
        msa_dir=msa_dir, output_dir=output_dir, prefix="matrix",
        seq_type="NT", taxa_occupancy=0.0, recoding=None,
        outgroup="C", to_format="fasta",
        translate_codon=False, exclude_codon3=False,
        dry_run=False, overwrite=False,
    )

    content = (output_dir / "matrix.fa").read_text()
    lines = content.strip().split("\n")
    assert lines[0] == ">C"


def test_run_concat_dry_run_writes_no_files(tmp_path: Path) -> None:
    from phyloai.pretree.concat import run_concat

    msa_dir = tmp_path / "msas"
    msa_dir.mkdir()
    (msa_dir / "gene1.fa").write_text(">A\nACGT\n>B\nACGT\n")

    output_dir = tmp_path / "out"
    payload = run_concat(
        msa_dir=msa_dir, output_dir=output_dir, prefix="matrix",
        seq_type="NT", taxa_occupancy=0.0, recoding=None,
        outgroup=None, to_format="fasta",
        translate_codon=False, exclude_codon3=False,
        dry_run=True, overwrite=False,
    )

    assert payload["status"] == "success"
    assert payload["key_results"]["n_taxa"] == 2
    assert not (output_dir / "matrix.fa").exists()
    assert not (output_dir / "result.json").exists()
    assert not (output_dir / "concat.log").exists()


def test_run_concat_dry_run_overwrite_does_not_delete_existing_output(tmp_path: Path) -> None:
    from phyloai.pretree.concat import run_concat

    msa_dir = tmp_path / "msas"
    msa_dir.mkdir()
    (msa_dir / "gene1.fa").write_text(">A\nACGT\n>B\nACGT\n")

    output_dir = tmp_path / "out"
    output_dir.mkdir()
    sentinel = output_dir / "old.txt"
    sentinel.write_text("keep")

    payload = run_concat(
        msa_dir=msa_dir, output_dir=output_dir, prefix="matrix",
        seq_type="NT", taxa_occupancy=0.0, recoding=None,
        outgroup=None, to_format="fasta",
        translate_codon=False, exclude_codon3=False,
        dry_run=True, overwrite=True,
    )

    assert payload["status"] == "success"
    assert sentinel.read_text() == "keep"


def test_run_concat_dry_run_reports_planned_variants(tmp_path: Path) -> None:
    from phyloai.pretree.concat import run_concat

    msa_dir = tmp_path / "msas"
    msa_dir.mkdir()
    (msa_dir / "gene1.fa").write_text(">A\nATGAAA\n>B\nATGAAA\n")

    output_dir = tmp_path / "out"
    payload = run_concat(
        msa_dir=msa_dir, output_dir=output_dir, prefix="matrix",
        seq_type="CODON", taxa_occupancy=0.0, recoding="RY-nucleotide",
        outgroup=None, to_format="fasta",
        translate_codon=True, exclude_codon3=True,
        dry_run=True, overwrite=False,
    )

    assert payload["key_results"]["variants_produced"] == [
        str(output_dir / "matrix.fa"),
        str(output_dir / "matrix.recoded.fa"),
        str(output_dir / "matrix.translated.fa"),
        str(output_dir / "matrix.cds12.fa"),
    ]
    assert [v["variant"] for v in payload["data"]["variants"]] == [
        "original", "recoded", "translated", "cds12",
    ]
    assert not output_dir.exists()


def test_run_concat_recoding_validation_rejects_aa_scheme_on_nt(tmp_path: Path) -> None:
    from phyloai.pretree.concat import run_concat

    msa_dir = tmp_path / "msas"
    msa_dir.mkdir()
    (msa_dir / "gene1.fa").write_text(">A\nACGT\n")

    output_dir = tmp_path / "out"
    with pytest.raises(ValueError, match="requires AA seq_type"):
        run_concat(
            msa_dir=msa_dir, output_dir=output_dir, prefix="matrix",
            seq_type="NT", taxa_occupancy=0.0, recoding="Dayhoff-6",
            outgroup=None, to_format="fasta",
            translate_codon=False, exclude_codon3=False,
            dry_run=False, overwrite=False,
        )


def test_run_concat_output_dir_conflict(tmp_path: Path) -> None:
    from phyloai.pretree.concat import run_concat

    msa_dir = tmp_path / "msas"
    msa_dir.mkdir()
    (msa_dir / "gene1.fa").write_text(">A\nACGT\n")

    output_dir = tmp_path / "out"
    output_dir.mkdir()
    (output_dir / "old.txt").write_text("old")

    with pytest.raises(ValueError, match="non-empty"):
        run_concat(
            msa_dir=msa_dir, output_dir=output_dir, prefix="matrix",
            seq_type="NT", taxa_occupancy=0.0, recoding=None,
            outgroup=None, to_format="fasta",
            translate_codon=False, exclude_codon3=False,
            dry_run=False, overwrite=False,
        )

    payload = run_concat(
        msa_dir=msa_dir, output_dir=output_dir, prefix="matrix",
        seq_type="NT", taxa_occupancy=0.0, recoding=None,
        outgroup=None, to_format="fasta",
        translate_codon=False, exclude_codon3=False,
        dry_run=False, overwrite=True,
    )
    assert payload["status"] == "success"


def test_run_concat_rejects_unknown_recoding_scheme(tmp_path: Path) -> None:
    from phyloai.pretree.concat import run_concat

    msa_dir = tmp_path / "msas"
    msa_dir.mkdir()
    (msa_dir / "gene1.fa").write_text(">A\nACGT\n")

    output_dir = tmp_path / "out"
    with pytest.raises(ValueError, match="Unknown recoding scheme"):
        run_concat(
            msa_dir=msa_dir, output_dir=output_dir, prefix="matrix",
            seq_type="NT", taxa_occupancy=0.0, recoding="NotARealScheme",
            outgroup=None, to_format="fasta",
            translate_codon=False, exclude_codon3=False,
            dry_run=False, overwrite=False,
        )


def test_run_concat_validation_error_writes_result_json(tmp_path: Path) -> None:
    from phyloai.pretree.concat import run_concat
    import json

    msa_dir = tmp_path / "msas"
    msa_dir.mkdir()
    (msa_dir / "gene1.fa").write_text(">A\nACGT\n")

    output_dir = tmp_path / "out"
    with pytest.raises(ValueError):
        run_concat(
            msa_dir=msa_dir, output_dir=output_dir, prefix="matrix",
            seq_type="NT", taxa_occupancy=0.0, recoding="NotARealScheme",
            outgroup=None, to_format="fasta",
            translate_codon=False, exclude_codon3=False,
            dry_run=False, overwrite=False,
        )

    assert output_dir.exists()
    result_path = output_dir / "result.json"
    assert result_path.exists()
    payload = json.loads(result_path.read_text())
    assert payload["status"] == "error"
    assert "NotARealScheme" in payload["error"]


def test_run_concat_validation_error_leaves_no_partial_matrix(tmp_path: Path) -> None:
    from phyloai.pretree.concat import run_concat

    msa_dir = tmp_path / "msas"
    msa_dir.mkdir()
    (msa_dir / "gene1.fa").write_text(">A\nACGT\n>B\nACGT\n")

    output_dir = tmp_path / "out"
    with pytest.raises(ValueError):
        run_concat(
            msa_dir=msa_dir, output_dir=output_dir, prefix="matrix",
            seq_type="NT", taxa_occupancy=0.0, recoding="NotARealScheme",
            outgroup=None, to_format="fasta",
            translate_codon=False, exclude_codon3=False,
            dry_run=False, overwrite=False,
        )

    assert not (output_dir / "matrix.fa").exists()
    assert not (output_dir / "matrix.partitions").exists()


def test_write_partitions_dna_three_genes(tmp_path: Path) -> None:
    from phyloai.pretree.concat import _write_partitions

    out_path = tmp_path / "matrix.partitions"
    genes = [("COI", 1, 654), ("16S", 655, 1203), ("CYTB", 1204, 1980)]
    _write_partitions(out_path, genes, "DNA")

    content = out_path.read_text()
    assert "DNA, COI = 1-654\n" in content
    assert "DNA, 16S = 655-1203\n" in content
    assert "DNA, CYTB = 1204-1980\n" in content


def test_write_partitions_lg_single_gene(tmp_path: Path) -> None:
    from phyloai.pretree.concat import _write_partitions

    out_path = tmp_path / "matrix.partitions"
    genes = [("GENE1", 1, 500)]
    _write_partitions(out_path, genes, "LG")

    content = out_path.read_text()
    assert content == "LG, GENE1 = 1-500\n"


def test_write_partitions_auto_prefix(tmp_path: Path) -> None:
    from phyloai.pretree.concat import _write_partitions

    out_path = tmp_path / "matrix.recoded.partitions"
    genes = [("geneA", 1, 100), ("geneB", 101, 200)]
    _write_partitions(out_path, genes, "AUTO")

    content = out_path.read_text()
    assert "AUTO, geneA = 1-100\n" in content
    assert "AUTO, geneB = 101-200\n" in content


def test_run_concat_writes_partitions_for_all_variants(tmp_path: Path) -> None:
    from phyloai.pretree.concat import run_concat

    msa_dir = tmp_path / "msas"
    msa_dir.mkdir()
    (msa_dir / "gene1.fa").write_text(">A\nACGT\n>B\nACGT\n>C\nACGT\n")
    (msa_dir / "gene2.fa").write_text(">A\nGGCC\n>B\nGGCC\n>C\nGGCC\n")

    output_dir = tmp_path / "out"
    run_concat(
        msa_dir=msa_dir, output_dir=output_dir, prefix="matrix",
        seq_type="NT", taxa_occupancy=0.5, recoding="RY-nucleotide",
        outgroup=None, to_format="fasta",
        translate_codon=False, exclude_codon3=False,
        dry_run=False, overwrite=False,
    )

    assert (output_dir / "matrix.partitions").exists()
    assert (output_dir / "matrix.recoded.partitions").exists()

    orig = (output_dir / "matrix.partitions").read_text()
    assert "DNA, gene1 = 1-4\n" in orig
    assert "DNA, gene2 = 5-8\n" in orig

    recoded = (output_dir / "matrix.recoded.partitions").read_text()
    assert "AUTO, gene1 = 1-4\n" in recoded
    assert "AUTO, gene2 = 5-8\n" in recoded


def test_run_concat_partitions_dry_run_no_files(tmp_path: Path) -> None:
    from phyloai.pretree.concat import run_concat

    msa_dir = tmp_path / "msas"
    msa_dir.mkdir()
    (msa_dir / "gene1.fa").write_text(">A\nACGT\n>B\nACGT\n")

    output_dir = tmp_path / "out"
    run_concat(
        msa_dir=msa_dir, output_dir=output_dir, prefix="matrix",
        seq_type="NT", taxa_occupancy=0.0, recoding=None,
        outgroup=None, to_format="fasta",
        translate_codon=False, exclude_codon3=False,
        dry_run=True, overwrite=False,
    )

    assert not (output_dir / "matrix.partitions").exists()


def test_run_concat_partitions_with_codon_variants(tmp_path: Path) -> None:
    from phyloai.pretree.concat import run_concat

    msa_dir = tmp_path / "msas"
    msa_dir.mkdir()
    (msa_dir / "gene1.fa").write_text(">A\nATGCGT\n>B\nATGCGT\n")

    output_dir = tmp_path / "out"
    run_concat(
        msa_dir=msa_dir, output_dir=output_dir, prefix="matrix",
        seq_type="CODON", taxa_occupancy=0.0, recoding=None,
        outgroup=None, to_format="fasta",
        translate_codon=True, exclude_codon3=True,
        dry_run=False, overwrite=False,
    )

    orig_p = (output_dir / "matrix.partitions").read_text()
    assert "DNA, gene1 = 1-6\n" in orig_p

    trans_p = (output_dir / "matrix.translated.partitions").read_text()
    assert "LG, gene1 = 1-2\n" in trans_p

    cds12_p = (output_dir / "matrix.cds12.partitions").read_text()
    assert "DNA, gene1 = 1-4\n" in cds12_p
