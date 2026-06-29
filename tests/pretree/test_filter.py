"""Tests for phyloai.pretree.filter."""

from pathlib import Path
from unittest.mock import patch

import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


class TestBuildTaperCmd:
    def test_aa_default(self):
        from phyloai.pretree.filter import _build_taper_cmd
        cmd = _build_taper_cmd(Path("a.fa"), Path("out.fa"), "AA", 3, "julia", "/t/c.jl", None)
        assert cmd[:3] == ["julia", "/t/c.jl", "-c"]
        assert "3" in cmd

    def test_nt_has_mode_flags(self):
        from phyloai.pretree.filter import _build_taper_cmd
        cmd = _build_taper_cmd(Path("a.fa"), Path("out.fa"), "NT", 3, "julia", "/t/c.jl", None)
        assert "-m" in cmd and "-a" in cmd

    def test_blocks_managed_flags(self):
        from phyloai.pretree.filter import _build_taper_cmd
        with pytest.raises(ValueError, match="managed"):
            _build_taper_cmd(Path("a.fa"), Path("out.fa"), "AA", 3, "julia", "/t/c.jl", "-c 5")

    def test_detect_taper_version_uses_shared_tool_probe(self, tmp_path):
        from phyloai.pretree.filter import _detect_taper_version
        taper = tmp_path / "correction_multi.jl"
        taper.write_text("# taper placeholder")

        with patch("subprocess.run") as mock_run:
            mock_run.return_value = type(
                "Result", (),
                {"stdout": "Version 1.0.0\n(type: Float64, default: 3.0)\n", "stderr": "", "returncode": 0}
            )()
            version = _detect_taper_version(str(taper), "/opt/julia/bin/julia")

        assert version == "1.0.0"
        mock_run.assert_called_once_with(
            ["/opt/julia/bin/julia", str(taper), "-h"],
            capture_output=True,
            text=True,
            timeout=5,
        )


class TestTaperCDSProjection:
    def test_projection(self, tmp_path):
        from phyloai.pretree.filter import _project_taper_masks_to_cds
        aa_original = tmp_path / "aa_original.fa"
        aa_masked = tmp_path / "aa_masked.fa"
        SeqIO.write([SeqRecord(Seq("AA-"), id="t1", description=""), SeqRecord(Seq("AX-"), id="t2", description="")], str(aa_original), "fasta")
        SeqIO.write([SeqRecord(Seq("XA-"), id="t1", description=""), SeqRecord(Seq("AX-"), id="t2", description="")], str(aa_masked), "fasta")
        nt_path = tmp_path / "nt.fa"
        SeqIO.write([SeqRecord(Seq("GCANNN---"), id="t1", description=""), SeqRecord(Seq("GCANNN---"), id="t2", description="")], str(nt_path), "fasta")
        out = tmp_path / "out.fna"
        result = _project_taper_masks_to_cds(aa_original, aa_masked, nt_path, out)
        assert result["projected_codons"] == 1
        assert out.exists()


class TestRetainedMsaStats:
    def test_empty(self):
        from phyloai.pretree.filter import _compute_retained_msa_stats
        s = _compute_retained_msa_stats([])
        assert s["n_msa"] == 0

    def test_one_fasta(self, tmp_path):
        from phyloai.pretree.filter import _compute_retained_msa_stats
        p = tmp_path / "g.fa"
        SeqIO.write([SeqRecord(Seq("ACGT"), id="a", description=""), SeqRecord(Seq("ACGT"), id="b", description="")], str(p), "fasta")
        s = _compute_retained_msa_stats([p])
        assert s["n_msa"] == 1
        assert s["mean_taxa"] == 2.0
        assert s["total_length"] == 4


class TestTreeshrink:
    def test_empty_tree_dir_raises(self, tmp_path):
        from phyloai.pretree.filter import run_treeshrink
        (tmp_path / "trees").mkdir()
        with pytest.raises(ValueError, match="No valid tree files"):
            run_treeshrink(tree_dir=tmp_path / "trees", output_dir=tmp_path / "out", dry_run=True)


class TestFilterCondition:
    def test_numeric_gte(self):
        from phyloai.pretree.filter import FilterCondition
        c = FilterCondition("dvmc", ">=", 0.3)
        assert c.evaluate({"dvmc": "0.5"})
        assert not c.evaluate({"dvmc": "0.1"})
        assert not c.evaluate({"dvmc": ""})

    def test_string_eq(self):
        from phyloai.pretree.filter import FilterCondition
        c = FilterCondition("DataType", "==", "AA")
        assert c.evaluate({"DataType": "AA"})
        assert not c.evaluate({"DataType": "NT"})


class TestParseKeepConditions:
    def test_simple(self):
        from phyloai.pretree.filter import parse_keep_conditions
        conds = parse_keep_conditions("dvmc>=0,dvmc<=0.3,average_BS>=0.8", {"dvmc", "average_BS", "num_sites"})
        assert len(conds) == 3

    def test_unknown_column_raises(self):
        from phyloai.pretree.filter import parse_keep_conditions
        with pytest.raises(ValueError, match="Unknown column"):
            parse_keep_conditions("badcol>=0", {"dvmc"})

    def test_malformed_raises(self):
        from phyloai.pretree.filter import parse_keep_conditions
        with pytest.raises(ValueError, match="Malformed"):
            parse_keep_conditions("not_valid", {"dvmc"})

    def test_rejects_non_eq_on_string(self):
        from phyloai.pretree.filter import parse_keep_conditions
        with pytest.raises(ValueError, match="can only use"):
            parse_keep_conditions("DataType>=AA", {"DataType", "dvmc"})


class TestClusterFeatureSelection:
    def test_selects_numeric_excludes_loci(self):
        from phyloai.pretree.filter import _select_features
        rows = [
            {"loci": "g1", "DataType": "AA", "dvmc": "0.1", "bs": "0.9", "name": "abc"},
            {"loci": "g2", "DataType": "AA", "dvmc": "0.2", "bs": "0.8", "name": "xyz"},
        ]
        cols = list(rows[0].keys())
        feats, entries = _select_features(rows, cols, None, [], "loci")
        assert "loci" not in feats
        assert "DataType" not in feats
        assert "name" not in feats
        assert "dvmc" in feats
        assert "bs" in feats
        assert len(entries) > 0
        assert any(e["column"] == "loci" and not e["included"] for e in entries)
        assert any(e["column"] == "name" and not e["included"] for e in entries)
        assert any(e["column"] == "dvmc" and e["included"] for e in entries)

    def test_cluster_outputs_are_grouped_by_stage(self, tmp_path):
        from phyloai.pretree.filter import run_cluster_filter

        table = tmp_path / "metrics.csv"
        table.write_text(
            "loci,dvmc,average_BS,gc_content\n"
            "g1,0.10,0.95,0.40\n"
            "g2,0.20,0.90,0.45\n"
            "g3,1.10,0.40,0.70\n"
            "g4,1.20,0.35,0.75\n"
        )
        out = tmp_path / "out"

        payload = run_cluster_filter(table_path=table, output_dir=out, n_clusters=2, quiet=True)

        assert payload["status"] == "success"
        assert (out / "01-input" / "features_used.csv").exists()
        assert (out / "02-reduction" / "reduction.csv").exists()
        assert (out / "03-clustering" / "clusters.csv").exists()
        assert (out / "03-clustering" / "cluster_loci").is_dir()
        assert (out / "04-diagnostics" / "cluster_metric_means.csv").exists()
        assert (out / "04-diagnostics" / "plots" / "cluster_2d.pdf").exists()
        assert not (out / "features_used.csv").exists()
        assert not (out / "clusters.csv").exists()
