"""Tests for phyloai.core.file_matching."""

from pathlib import Path

import pytest

from phyloai.core.file_matching import (
    logical_msa_locus_name,
    logical_tree_locus_candidates,
    pair_msa_and_tree_maps,
    scan_msa_dir,
    scan_tree_dir,
)


class TestLogicalMsaLocusName:
    def test_basic(self):
        assert logical_msa_locus_name(Path("gene1.fa")) == "gene1"

    def test_dotted(self):
        assert logical_msa_locus_name(Path("gene.v1.ALI")) == "gene.v1"

    def test_no_suffix(self):
        assert logical_msa_locus_name(Path("gene1")) == "gene1"

    def test_uppercase(self):
        assert logical_msa_locus_name(Path("gene1.FASTA")) == "gene1"


class TestLogicalTreeLocusCandidates:
    def test_one_suffix(self):
        a, b = logical_tree_locus_candidates(Path("gene.treefile"))
        assert a == "gene"
        assert b is None

    def test_two_suffix(self):
        a, b = logical_tree_locus_candidates(Path("gene.fa.treefile"))
        assert a == "gene.fa"
        assert b == "gene"

    def test_dotted_two_suffix(self):
        a, b = logical_tree_locus_candidates(Path("gene.v1.fa.treefile"))
        assert a == "gene.v1.fa"
        assert b == "gene.v1"


class TestPairMsaAndTreeMaps:
    def test_direct_match(self):
        msa_map = {"gene1": Path("gene1.fa"), "gene2": Path("gene2.fa")}
        tree_paths = [Path("gene1.treefile"), Path("gene2.tre")]
        result = pair_msa_and_tree_maps(msa_map, tree_paths)
        assert result.paired["gene1"][1] == tree_paths[0]
        assert result.paired["gene2"][1] == tree_paths[1]

    def test_ambiguous_raises(self):
        msa_map = {"gene.fa": Path("g1.fa"), "gene": Path("g2.fa")}
        tree_paths = [Path("gene.fa.treefile")]
        with pytest.raises(ValueError, match="ambiguous"):
            pair_msa_and_tree_maps(msa_map, tree_paths)

    def test_unpaired_tree_warns(self):
        msa_map = {"gene1": Path("gene1.fa")}
        tree_paths = [Path("gene2.treefile")]
        result = pair_msa_and_tree_maps(msa_map, tree_paths)
        assert result.paired.get("gene2") is not None
        assert result.paired["gene2"][0] is None
        assert len(result.warnings) == 1


class TestScanMsaDir:
    def test_returns_locus_map(self, tmp_path):
        (tmp_path / "gene1.fa").write_text(">a\nACGT\n")
        (tmp_path / "gene2.FASTA").write_text(">b\nACGT\n")
        (tmp_path / "gene3.tre").write_text("(a,b);")
        (tmp_path / "notes.txt").write_text("hello")
        (tmp_path / "subdir").mkdir()

        result = scan_msa_dir(tmp_path)
        assert "gene1" in result
        assert "gene2" in result
        assert "gene3" in result
        assert "notes" in result
        assert len(result) == 4

    def test_empty_on_nonexistent(self):
        assert scan_msa_dir(Path("/nonexistent")) == {}

    def test_skips_empty_files(self, tmp_path):
        (tmp_path / "empty.fa").write_text("")
        (tmp_path / "gene1.fa").write_text(">a\nACGT\n")
        result = scan_msa_dir(tmp_path)
        assert "empty" not in result
        assert "gene1" in result


class TestScanTreeDir:
    def test_returns_locus_map(self, tmp_path):
        (tmp_path / "gene1.treefile").write_text("(a,b);")
        (tmp_path / "gene2.tre").write_text("(c,d);")
        (tmp_path / "gene3.fa").write_text(">a\nACGT\n")

        result = scan_tree_dir(tmp_path)
        assert "gene1" in result
        assert "gene2" in result
        assert "gene3" in result
        assert len(result) == 3

    def test_handles_two_suffix(self, tmp_path):
        (tmp_path / "gene.v1.treefile").write_text("(a,b);")
        result = scan_tree_dir(tmp_path)
        assert len(result) == 1

    def test_handles_two_suffix_with_conflict(self, tmp_path):
        (tmp_path / "gene.fa.treefile").write_text("(a,b);")
        (tmp_path / "gene.treefile").write_text("(c,d);")
        result = scan_tree_dir(tmp_path)
        assert len(result) == 2
        assert "gene.fa" in result
        assert "gene" in result

    def test_duplicate_candidates_raises(self, tmp_path):
        (tmp_path / "gene.treefile").write_text("(a,b);")
        (tmp_path / "gene.tre").write_text("(c,d);")
        with pytest.raises(ValueError, match="Duplicate or ambiguous tree file"):
            scan_tree_dir(tmp_path)

    def test_empty_on_nonexistent(self):
        assert scan_tree_dir(Path("/nonexistent")) == {}
