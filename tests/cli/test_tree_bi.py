from __future__ import annotations

from pathlib import Path

from click.testing import CliRunner

from phyloai.cli.main import cli


def test_tree_bi_pb_help_shows_flags():
    result = CliRunner().invoke(cli, ["tree", "bi", "pb", "--help"])
    assert result.exit_code == 0
    for flag in ["--matrix", "--model", "--mixture", "--gamma-cats", "--chains", "--threads", "--resume", "--pb-path", "--poll-interval"]:
        assert flag in result.output


def test_tree_group_shows_bi():
    result = CliRunner().invoke(cli, ["tree", "--help"])
    assert result.exit_code == 0
    assert "bi" in result.output


def test_tree_bi_group_help():
    result = CliRunner().invoke(cli, ["tree", "bi", "--help"])
    assert result.exit_code == 0
    for sub in ["pb", "bpcomp", "tracecomp", "readpb"]:
        assert sub in result.output


def test_tree_bi_bpcomp_help():
    result = CliRunner().invoke(cli, ["tree", "bi", "bpcomp", "--help"])
    assert result.exit_code == 0
    for flag in ["--chain-dir", "--burnin", "--cutoff"]:
        assert flag in result.output


def test_tree_bi_tracecomp_help():
    result = CliRunner().invoke(cli, ["tree", "bi", "tracecomp", "--help"])
    assert result.exit_code == 0
    for flag in ["--chain-dir", "--burnin"]:
        assert flag in result.output


def test_tree_bi_readpb_help():
    result = CliRunner().invoke(cli, ["tree", "bi", "readpb", "--help"])
    assert result.exit_code == 0
    for flag in ["--chain", "--mode", "--threads"]:
        assert flag in result.output


def test_tree_bi_pb_dry_run(tmp_path: Path):
    matrix = tmp_path / "m.phy"
    matrix.write_text("2 3\na AAA\nb AAA\n")
    result = CliRunner().invoke(cli, ["tree", "bi", "pb", "--matrix", str(matrix), "--output-dir", str(tmp_path / "out"), "--dry-run"])
    assert result.exit_code == 0
    assert "pb_mpi" in result.output


def test_tree_bi_pb_resume_bare_parses(tmp_path: Path):
    out = tmp_path / "out"
    out.mkdir()
    (out / "run_state.json").write_text('{"chain_names": ["chain1"], "matrix": "m.phy", "model_flags": ["-cat", "-gtr", "-dgam", "4"], "sample_freq": 1, "nsamples": -1, "threads": 4}')
    result = CliRunner().invoke(cli, ["tree", "bi", "pb", "--output-dir", str(out), "--resume", "--dry-run", "--quiet"])
    assert result.exit_code == 0


def test_tree_bi_pb_dry_run_fails_missing_matrix():
    result = CliRunner().invoke(cli, ["tree", "bi", "pb", "--matrix", "/nonexistent.phy", "--dry-run"])
    assert result.exit_code == 1
    assert "does not exist" in result.output


def test_tree_bi_pb_dry_run_phylogeny_with_fasta(tmp_path: Path):
    matrix = tmp_path / "test.fa"
    matrix.write_text(">a\nAAA\n>b\nAAA\n")
    result = CliRunner().invoke(cli, ["tree", "bi", "pb", "--matrix", str(matrix), "--output-dir", str(tmp_path / "out"), "--dry-run"])
    assert result.exit_code == 0
    assert "pb_mpi" in result.output


def test_tree_bi_bpcomp_dry_run(tmp_path: Path):
    chain_dir = tmp_path / "chains"
    chain_dir.mkdir()
    (chain_dir / "chain1.chain").write_text("")
    (chain_dir / "chain2.chain").write_text("")
    result = CliRunner().invoke(cli, ["tree", "bi", "bpcomp", "--chain-dir", str(chain_dir), "--output-dir", str(tmp_path / "bpcomp"), "--burnin", "1000", "--dry-run"])
    assert result.exit_code == 0


def test_tree_bi_tracecomp_dry_run(tmp_path: Path):
    chain_dir = tmp_path / "chains"
    chain_dir.mkdir()
    (chain_dir / "chain1.trace").write_text("")
    (chain_dir / "chain2.trace").write_text("")
    result = CliRunner().invoke(cli, ["tree", "bi", "tracecomp", "--chain-dir", str(chain_dir), "--output-dir", str(tmp_path / "tracecomp"), "--dry-run"])
    assert result.exit_code == 0


def test_tree_bi_readpb_dry_run(tmp_path: Path):
    chain = tmp_path / "chain1"
    (tmp_path / "chain1.chain").write_text("fake state")
    result = CliRunner().invoke(cli, ["tree", "bi", "readpb", "--chain", str(chain), "--mode", "rr", "--output-dir", str(tmp_path / "readpb"), "--dry-run"])
    assert result.exit_code == 0
