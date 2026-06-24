from __future__ import annotations

from pathlib import Path

from click.testing import CliRunner

from phyloai.cli.main import cli


def test_tree_bi_help_shows_flags():
    result = CliRunner().invoke(cli, ["tree", "bi", "--help"])
    assert result.exit_code == 0
    for flag in ["--matrix", "--model", "--mixture", "--gamma-cats", "--chains", "--threads", "--resume", "--pb-path", "--poll-interval"]:
        assert flag in result.output


def test_tree_group_shows_bi():
    result = CliRunner().invoke(cli, ["tree", "--help"])
    assert result.exit_code == 0
    assert "bi" in result.output


def test_tree_bi_dry_run(tmp_path: Path):
    matrix = tmp_path / "m.phy"
    matrix.write_text("2 3\na AAA\nb AAA\n")
    result = CliRunner().invoke(cli, ["tree", "bi", "--matrix", str(matrix), "--output-dir", str(tmp_path / "out"), "--dry-run"])
    assert result.exit_code == 0
    assert "pb_mpi" in result.output


def test_tree_bi_resume_bare_parses(tmp_path: Path):
    out = tmp_path / "out"
    out.mkdir()
    (out / "run_state.json").write_text('{"chain_names": ["chain1"], "matrix": "m.phy", "model_flags": ["-cat", "-gtr", "-dgam", "4"], "sample_freq": 1, "nsamples": -1, "threads": 4}')
    result = CliRunner().invoke(cli, ["tree", "bi", "--output-dir", str(out), "--resume", "--dry-run", "--quiet"])
    assert result.exit_code == 0


def test_tree_bi_dry_run_fails_missing_matrix():
    result = CliRunner().invoke(cli, ["tree", "bi", "--matrix", "/nonexistent.phy", "--dry-run"])
    assert result.exit_code == 1
    assert "does not exist" in result.output


def test_tree_bi_dry_run_phylogeny_with_fasta(tmp_path: Path):
    matrix = tmp_path / "test.fa"
    matrix.write_text(">a\nAAA\n>b\nAAA\n")
    result = CliRunner().invoke(cli, ["tree", "bi", "--matrix", str(matrix), "--output-dir", str(tmp_path / "out"), "--dry-run"])
    assert result.exit_code == 0
    assert "pb_mpi" in result.output
