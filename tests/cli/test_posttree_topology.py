"""CLI tests for phyloai posttree topology."""
from __future__ import annotations

import json
import shutil
from pathlib import Path

import pytest
from click.testing import CliRunner

from phyloai.cli.main import cli


@pytest.fixture
def runner() -> CliRunner:
    return CliRunner()


@pytest.fixture
def matrix_file(tmp_path: Path) -> Path:
    p = tmp_path / "matrix.fa"
    p.write_text(">a\nMKTLLL\n>b\nMKTLLL\n")
    return p


@pytest.fixture
def candidate_trees_file(tmp_path: Path) -> Path:
    p = tmp_path / "candidates.trees"
    p.write_text("(a,b);\n")
    return p


# ------------------------------------------------------------------
# Help content
# ------------------------------------------------------------------

class TestCLIHelp:
    def test_help_contains_key_sections(self, runner: CliRunner) -> None:
        result = runner.invoke(cli, ["posttree", "topology", "--help"])
        assert result.exit_code == 0
        # Help now uses docstring-style text (no \\b markers).
        # Check that the essential content is present.
        assert "Does NOT infer new trees" in result.output
        assert "PhyloAI does NOT expose ModelFinder" in result.output
        assert "Default tests:" in result.output
        assert "Examples:" in result.output
        assert "p-values" in result.output
        assert "Recommended: AU, WSH, WKH" in result.output

    def test_help_contains_all_six_examples(self, runner: CliRunner) -> None:
        result = runner.invoke(cli, ["posttree", "topology", "--help"])
        assert result.exit_code == 0
        assert "LG+F+R4" in result.output
        assert "C20+F+R4" in result.output
        assert "LG+C20+F+R4" in result.output
        assert "--partitions raw.best_model.nex" in result.output
        assert "h1.nwk" in result.output
        assert "h3.nwk" in result.output
        assert "custom.exchangeabilities" in result.output

    def test_help_shows_all_options(self, runner: CliRunner) -> None:
        result = runner.invoke(cli, ["posttree", "topology", "--help"])
        assert result.exit_code == 0
        for opt in ("--matrix", "--candidate-trees", "--model-expr", "--partitions",
                     "--guide-tree", "--replicates", "--prefix", "--output-dir",
                     "--threads", "--iqtree-path", "--tool-args",
                     "--overwrite", "--resume", "--dry-run"):
            assert opt in result.output


# ------------------------------------------------------------------
# Input validation (all exit code 1 per spec)
# ------------------------------------------------------------------

class TestCLIValidation:
    def test_missing_matrix_exits_1(self, runner: CliRunner) -> None:
        result = runner.invoke(cli, ["posttree", "topology"])
        assert result.exit_code == 1
        assert "matrix" in result.output.lower()

    def test_missing_candidate_trees_exits_1(
        self, runner: CliRunner, matrix_file: Path,
    ) -> None:
        result = runner.invoke(cli, [
            "posttree", "topology", "--matrix", str(matrix_file),
        ])
        assert result.exit_code == 1
        assert "candidate" in result.output.lower()

    def test_missing_model_source_exits_1(
        self, runner: CliRunner, matrix_file: Path, candidate_trees_file: Path,
    ) -> None:
        result = runner.invoke(cli, [
            "posttree", "topology",
            "--matrix", str(matrix_file),
            "--candidate-trees", str(candidate_trees_file),
        ])
        assert result.exit_code == 1
        assert "model" in result.output.lower()

    def test_replicates_below_minimum(
        self, runner: CliRunner, matrix_file: Path, candidate_trees_file: Path,
    ) -> None:
        result = runner.invoke(cli, [
            "posttree", "topology",
            "--matrix", str(matrix_file),
            "--candidate-trees", str(candidate_trees_file),
            "--model-expr", "LG+F+R4",
            "--replicates", "999",
        ])
        assert result.exit_code == 1

    def test_overwrite_and_resume_mutually_exclusive(
        self, runner: CliRunner, matrix_file: Path, candidate_trees_file: Path,
    ) -> None:
        result = runner.invoke(cli, [
            "posttree", "topology",
            "--matrix", str(matrix_file),
            "--candidate-trees", str(candidate_trees_file),
            "--model-expr", "LG+F+R4",
            "--overwrite", "--resume",
        ])
        assert result.exit_code == 1
        assert "mutually exclusive" in result.output.lower()

    def test_both_model_expr_and_partitions(
        self, runner: CliRunner, matrix_file: Path, candidate_trees_file: Path,
        tmp_path: Path,
    ) -> None:
        (tmp_path / "m.best_model.nex").write_text("#nexus\n")
        result = runner.invoke(cli, [
            "posttree", "topology",
            "--matrix", str(matrix_file),
            "--candidate-trees", str(candidate_trees_file),
            "--model-expr", "LG+F+R4",
            "--partitions", str(tmp_path / "m.best_model.nex"),
        ])
        assert result.exit_code == 1
        assert "mutually exclusive" in result.output.lower()

    def test_blocked_s_in_tool_args(
        self, runner: CliRunner, matrix_file: Path, candidate_trees_file: Path,
    ) -> None:
        result = runner.invoke(cli, [
            "posttree", "topology",
            "--matrix", str(matrix_file),
            "--candidate-trees", str(candidate_trees_file),
            "--model-expr", "LG+F+R4",
            "--tool-args", "-s other.fa",
        ])
        assert result.exit_code == 1
        assert "-s" in result.output

    def test_blocked_z_in_tool_args(
        self, runner: CliRunner, matrix_file: Path, candidate_trees_file: Path,
    ) -> None:
        result = runner.invoke(cli, [
            "posttree", "topology",
            "--matrix", str(matrix_file),
            "--candidate-trees", str(candidate_trees_file),
            "--model-expr", "LG+F+R4",
            "--tool-args", "-z other.trees",
        ])
        assert result.exit_code == 1
        assert "-z" in result.output

    def test_tool_args_accepted(
        self, runner: CliRunner, matrix_file: Path, candidate_trees_file: Path,
        tmp_path: Path,
    ) -> None:
        out = tmp_path / "out"
        result = runner.invoke(cli, [
            "posttree", "topology",
            "--matrix", str(matrix_file),
            "--candidate-trees", str(candidate_trees_file),
            "--model-expr", "LG+F+R4",
            "--tool-args", "--prefix custom -T 30 -fs custom.sitefreq",
            "--output-dir", str(out),
            "--dry-run",
        ])
        assert result.exit_code == 0
        assert "iqtree3" in result.output


# ------------------------------------------------------------------
# Dry-run
# ------------------------------------------------------------------

class TestCLIDryRun:
    def test_dry_run_single_tree_file(
        self, runner: CliRunner, matrix_file: Path, candidate_trees_file: Path,
        tmp_path: Path,
    ) -> None:
        out = tmp_path / "out"
        result = runner.invoke(cli, [
            "posttree", "topology",
            "--matrix", str(matrix_file),
            "--candidate-trees", str(candidate_trees_file),
            "--model-expr", "LG+F+R4",
            "--output-dir", str(out),
            "--dry-run",
        ])
        assert result.exit_code == 0
        assert "Would run:" in result.output
        assert "iqtree3" in result.output

    def test_dry_run_multiple_tree_files(
        self, runner: CliRunner, matrix_file: Path, tmp_path: Path,
    ) -> None:
        (tmp_path / "h1.nwk").write_text("(a,b);\n")
        (tmp_path / "h2.nwk").write_text("(a,c);\n")
        out = tmp_path / "out"
        result = runner.invoke(cli, [
            "posttree", "topology",
            "--matrix", str(matrix_file),
            "--candidate-trees", f"{tmp_path / 'h1.nwk'},{tmp_path / 'h2.nwk'}",
            "--model-expr", "LG+F+R4",
            "--output-dir", str(out),
            "--dry-run",
        ])
        assert result.exit_code == 0
        assert "Would run:" in result.output


# ------------------------------------------------------------------
# Integration (real IQ-TREE)
# ------------------------------------------------------------------

class TestCLIIntegration:
    @pytest.mark.skipif(
        not shutil.which("iqtree3"),
        reason="iqtree3 not found in PATH",
    )
    def test_successful_run_writes_result_json(
        self, runner: CliRunner, tmp_path: Path,
    ) -> None:
        matrix = tmp_path / "matrix.fa"
        matrix.write_text(
            ">t1\nMKTLLLTLWVV\n>t2\nMKTLLLTLWVI\n>t3\nMKTLLLSLWVI\n>t4\nMKTLLLTLWVA\n"
        )
        (tmp_path / "trees").write_text(
            "(t1,t2,(t3,t4));\n(t1,t3,(t2,t4));\n"
        )
        out = tmp_path / "out"

        result = runner.invoke(cli, [
            "posttree", "topology",
            "--matrix", str(matrix),
            "--candidate-trees", str(tmp_path / "trees"),
            "--model-expr", "LG",
            "--replicates", "1000",
            "--output-dir", str(out),
            "--threads", "1",
        ])
        assert result.exit_code == 0
        assert (out / "result.json").exists()

        with open(out / "result.json") as fh:
            payload = json.load(fh)
        assert payload["status"] == "success"
        assert "iqtree3" in payload["tool_versions"]
        assert len(payload["data"]["tests"]) == 2
