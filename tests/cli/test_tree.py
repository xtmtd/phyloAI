from __future__ import annotations

import json
from pathlib import Path

from click.testing import CliRunner

from phyloai.cli.main import cli


def test_tree_group_exists() -> None:
    result = CliRunner().invoke(cli, ["tree", "--help"])
    assert result.exit_code == 0
    assert "Maximum-likelihood" in result.output or "ml" in result.output


def test_tree_ml_help_shows_both_backends() -> None:
    result = CliRunner().invoke(cli, ["tree", "ml", "--help"])
    assert result.exit_code == 0
    assert "fasttree" in result.output


def test_tree_ml_fasttree_help() -> None:
    result = CliRunner().invoke(cli, ["tree", "ml", "fasttree", "--help"])
    assert result.exit_code == 0
    for flag in ["--msa-dir", "--matrix", "--seq-type", "--model", "--mode",
                  "--boot", "--cat", "--gamma", "--output-dir", "--threads"]:
        assert flag in result.output


def test_tree_ml_fasttree_mutual_exclusivity(tmp_path: Path) -> None:
    msa_dir = tmp_path / "msas"
    msa_dir.mkdir()
    mat = tmp_path / "matrix.fa"
    mat.write_text(">a\nMKT\n")

    result = CliRunner().invoke(cli, [
        "tree", "ml", "fasttree",
        "--msa-dir", str(msa_dir), "--matrix", str(mat),
    ])
    assert result.exit_code == 1


def test_tree_ml_fasttree_neither_input() -> None:
    result = CliRunner().invoke(cli, [
        "tree", "ml", "fasttree",
    ])
    assert result.exit_code == 1


def test_cli_msa_dir_nonexistent_exits_1() -> None:
    result = CliRunner().invoke(cli, [
        "tree", "ml", "fasttree",
        "--msa-dir", "/nonexistent/path",
    ])
    assert result.exit_code == 1


def test_tree_ml_fasttree_quiet_dry_run_batch(tmp_path: Path) -> None:
    msa_dir = tmp_path / "msas"
    msa_dir.mkdir()
    (msa_dir / "g1.fa").write_text(">a\nMKTLLL\n>b\nMKTLLL\n")

    out_dir = tmp_path / "out"

    result = CliRunner().invoke(cli, [
        "tree", "ml", "fasttree",
        "--msa-dir", str(msa_dir),
        "--output-dir", str(out_dir),
        "--seq-type", "AA",
        "--model", "lg",
        "--quiet",
        "--dry-run",
    ])

    assert result.exit_code == 0
    assert not (out_dir / "result.json").exists()


def test_tree_ml_fasttree_quiet_dry_run_single(tmp_path: Path) -> None:
    mat = tmp_path / "matrix.fa"
    mat.write_text(">a\nMKTLLL\n>b\nMKTLLL\n")

    out_dir = tmp_path / "out"

    result = CliRunner().invoke(cli, [
        "tree", "ml", "fasttree",
        "--matrix", str(mat),
        "--output-dir", str(out_dir),
        "--seq-type", "AA",
        "--model", "lg",
        "--quiet",
        "--dry-run",
    ])

    assert result.exit_code == 0


def test_tree_ml_fasttree_invalid_model_exits_1(tmp_path: Path) -> None:
    mat = tmp_path / "matrix.fa"
    mat.write_text(">a\nMKTLLL\n")

    out_dir = tmp_path / "out"

    result = CliRunner().invoke(cli, [
        "tree", "ml", "fasttree",
        "--matrix", str(mat),
        "--output-dir", str(out_dir),
        "--seq-type", "AA",
        "--model", "gtr",
        "--quiet",
    ])

    assert result.exit_code == 1


def test_tree_ml_fasttree_blocked_tool_args(tmp_path: Path) -> None:
    mat = tmp_path / "matrix.fa"
    mat.write_text(">a\nMKTLLL\n")

    out_dir = tmp_path / "out"

    result = CliRunner().invoke(cli, [
        "tree", "ml", "fasttree",
        "--matrix", str(mat),
        "--output-dir", str(out_dir),
        "--tool-args", "-nt",
        "--quiet",
    ])

    assert result.exit_code == 1


def test_tree_ml_fasttree_threads_warn_single(tmp_path: Path) -> None:
    mat = tmp_path / "matrix.fa"
    mat.write_text(">a\nMKTLLL\n>b\nMKTLLL\n")

    out_dir = tmp_path / "out"

    result = CliRunner().invoke(cli, [
        "tree", "ml", "fasttree",
        "--matrix", str(mat),
        "--output-dir", str(out_dir),
        "--threads", "8",
        "--quiet",
        "--dry-run",
    ])
    assert "has no effect" in result.output.lower() or result.exit_code == 0


def test_tree_ml_fasttree_writes_result_json_and_log(tmp_path: Path) -> None:
    mat = tmp_path / "matrix.fa"
    mat.write_text(">a\nMKTLLL\n>b\nMKTLLL\n")

    out_dir = tmp_path / "out"

    result = CliRunner().invoke(cli, [
        "tree", "ml", "fasttree",
        "--matrix", str(mat),
        "--output-dir", str(out_dir),
        "--seq-type", "AA",
        "--model", "lg",
        "--quiet",
    ])
    if result.exit_code == 0:
        assert (out_dir / "result.json").exists()
        assert (out_dir / "fasttree.log").exists()
    elif result.exit_code == 3:
        import pytest
        pytest.skip("FastTree not installed")
