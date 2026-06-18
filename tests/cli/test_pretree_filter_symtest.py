"""CLI tests for phyloai pretree filter symtest."""

from pathlib import Path

from click.testing import CliRunner

from phyloai.cli.main import cli


def _make_msa_dir(tmp_path: Path) -> Path:
    """Create a minimal MSA directory for testing."""
    msa_dir = tmp_path / "msa"
    msa_dir.mkdir()
    (msa_dir / "gene1.fa").write_text(">t1\nACGT\n>t2\nACGT\n")
    (msa_dir / "gene2.fa").write_text(">t1\nTTTT\n>t2\nGGGG\n")
    return msa_dir


def test_symtest_help():
    runner = CliRunner()
    result = runner.invoke(cli, ["pretree", "filter", "symtest", "--help"])
    assert result.exit_code == 0
    assert "--msa-dir" in result.output
    assert "--symtest-pval" in result.output


def test_symtest_requires_msa_dir():
    runner = CliRunner()
    result = runner.invoke(cli, ["pretree", "filter", "symtest"])
    assert result.exit_code != 0


def test_symtest_dry_run(tmp_path):
    msa_dir = _make_msa_dir(tmp_path)
    output_dir = tmp_path / "out"
    runner = CliRunner()
    result = runner.invoke(cli, [
        "pretree", "filter", "symtest",
        "--msa-dir", str(msa_dir),
        "--output-dir", str(output_dir),
        "--dry-run",
    ])
    assert result.exit_code == 0
    assert "Dry run" in result.output
    assert not (output_dir / "result.json").exists()


def test_symtest_invalid_pval(tmp_path):
    msa_dir = _make_msa_dir(tmp_path)
    runner = CliRunner()
    result = runner.invoke(cli, [
        "pretree", "filter", "symtest",
        "--msa-dir", str(msa_dir),
        "--symtest-pval", "2.0",
    ])
    assert result.exit_code != 0


def test_symtest_invalid_pval_zero(tmp_path):
    msa_dir = _make_msa_dir(tmp_path)
    runner = CliRunner()
    result = runner.invoke(cli, [
        "pretree", "filter", "symtest",
        "--msa-dir", str(msa_dir),
        "--symtest-pval", "0",
    ])
    assert result.exit_code != 0


def test_symtest_output_dir_conflict(tmp_path):
    msa_dir = _make_msa_dir(tmp_path)
    output_dir = tmp_path / "out"
    output_dir.mkdir()
    (output_dir / "existing.txt").write_text("data")
    runner = CliRunner()
    result = runner.invoke(cli, [
        "pretree", "filter", "symtest",
        "--msa-dir", str(msa_dir),
        "--output-dir", str(output_dir),
    ])
    assert result.exit_code != 0
    assert "already exists" in result.output


def test_symtest_overwrite(tmp_path):
    msa_dir = _make_msa_dir(tmp_path)
    output_dir = tmp_path / "out"
    output_dir.mkdir()
    (output_dir / "old.txt").write_text("old")
    runner = CliRunner()
    result = runner.invoke(cli, [
        "pretree", "filter", "symtest",
        "--msa-dir", str(msa_dir),
        "--output-dir", str(output_dir),
        "--overwrite",
    ])
    assert "already exists" not in result.output


def test_symtest_missing_msa_dir(tmp_path):
    runner = CliRunner()
    result = runner.invoke(cli, [
        "pretree", "filter", "symtest",
        "--msa-dir", str(tmp_path / "nonexistent"),
    ])
    assert result.exit_code != 0


def test_symtest_threads_negative(tmp_path):
    msa_dir = _make_msa_dir(tmp_path)
    runner = CliRunner()
    result = runner.invoke(cli, [
        "pretree", "filter", "symtest",
        "--msa-dir", str(msa_dir),
        "--threads", "-1",
    ])
    assert result.exit_code != 0
