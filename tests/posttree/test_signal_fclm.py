"""Tests for phyloai.posttree.signal fclm subcommand."""
from __future__ import annotations

from pathlib import Path

import pytest


class TestFclmValidation:
    def test_reads_phylip_taxa_for_taxset_validation(self, tmp_path: Path) -> None:
        from phyloai.posttree.signal import run_signal_fclm

        matrix = tmp_path / "m.phy"
        matrix.write_text("4 3\nA MKT\nB MKA\nC MKL\nD MKP\n")
        csv = tmp_path / "taxsets.csv"
        csv.write_text("taxon,taxset\nA,G1\nB,G2\nC,G3\nD,G4\n")
        result = run_signal_fclm(
            matrix=matrix, taxset_csv=csv, model_expr="LG+F+R4", lmap=None,
            guide_tree=None, threads="auto", iqtree_path=None, tool_args=None,
            output_dir=tmp_path / "out", overwrite=False, dry_run=True, quiet=True,
        )
        assert result["status"] == "success"
        assert result["key_results"]["n_taxsets"] == 4

    def test_taxset_csv_extra_taxa_errors(self, tmp_path: Path) -> None:
        from phyloai.posttree.signal import run_signal_fclm

        matrix = tmp_path / "m.fa"
        matrix.write_text(">A\nMKT\n>B\nMKA\n>C\nMKL\n>D\nMKP\n")
        csv = tmp_path / "taxsets.csv"
        csv.write_text("taxon,taxset\nA,G1\nB,G2\nC,G3\nD,G4\nX,G1\n")
        result = run_signal_fclm(
            matrix=matrix, taxset_csv=csv,
            model_expr="LG+F+R4", lmap=None,
            guide_tree=None, threads="auto",
            iqtree_path=None, tool_args=None,
            output_dir=tmp_path / "out",
            overwrite=False, dry_run=True, quiet=True,
        )
        assert result["status"] == "error"
        assert "X" in result["error"] or "taxa" in result["error"].lower()

    def test_taxset_csv_fewer_than_4_taxsets_errors(self, tmp_path: Path) -> None:
        from phyloai.posttree.signal import run_signal_fclm

        matrix = tmp_path / "m.fa"
        matrix.write_text(">A\nMKT\n>B\nMKA\n>C\nMKL\n>D\nMKP\n")
        csv = tmp_path / "taxsets.csv"
        csv.write_text("taxon,taxset\nA,G1\nB,G1\nC,G2\nD,G3\n")
        result = run_signal_fclm(
            matrix=matrix, taxset_csv=csv,
            model_expr="LG+F+R4", lmap=None,
            guide_tree=None, threads="auto",
            iqtree_path=None, tool_args=None,
            output_dir=tmp_path / "out",
            overwrite=False, dry_run=True, quiet=True,
        )
        assert result["status"] == "error"
        assert "4" in result["error"]

    def test_malformed_matrix_returns_input_error(self, tmp_path: Path) -> None:
        from phyloai.posttree.signal import run_signal_fclm

        matrix = tmp_path / "m.fa"
        matrix.write_text("not a valid fasta")
        csv = tmp_path / "taxsets.csv"
        csv.write_text("taxon,taxset\nA,G1\nB,G2\nC,G3\nD,G4\n")
        result = run_signal_fclm(
            matrix=matrix, taxset_csv=csv,
            model_expr="LG+F+R4", lmap=None,
            guide_tree=None, threads="auto",
            iqtree_path=None, tool_args=None,
            output_dir=tmp_path / "out",
            overwrite=False, dry_run=True, quiet=True,
        )
        assert result["status"] == "error"
        assert result["error_category"] == "input"
        assert "matrix" in result["error"].lower()

    def test_dry_run_shows_command_without_writing_files(self, tmp_path: Path) -> None:
        from phyloai.posttree.signal import run_signal_fclm

        matrix = tmp_path / "m.fa"
        matrix.write_text(">A\nMKT\n>B\nMKA\n>C\nMKL\n>D\nMKP\n>E\nMKQ\n")
        csv = tmp_path / "taxsets.csv"
        csv.write_text("taxon,taxset\nA,G1\nB,G2\nC,G3\nD,G4\nE,G1\n")
        result = run_signal_fclm(
            matrix=matrix, taxset_csv=csv,
            model_expr="LG+F+R4", lmap=None,
            guide_tree=None, threads="auto",
            iqtree_path=None, tool_args=None,
            output_dir=tmp_path / "out",
            overwrite=False, dry_run=True, quiet=True,
        )
        assert result["status"] == "success"
        assert "-lmap" in result["data"]["cmd"]
        assert "-lmclust" in result["data"]["cmd"]
        assert result["key_results"]["n_taxsets"] == 4
        nexus = tmp_path / "out" / "cluster.nexus"
        assert not nexus.exists()
        assert not (tmp_path / "out").exists()
