"""Smoke tests for `phyloai posttree dating` CLI."""
from __future__ import annotations
from pathlib import Path
from click.testing import CliRunner
from phyloai.cli.main import cli
import json


def test_dating_group_help():
    runner = CliRunner()
    result = runner.invoke(cli, ["posttree", "dating", "--help"])
    assert result.exit_code == 0
    assert "hessian" in result.output
    assert "mcmc" in result.output


def test_hessian_help():
    runner = CliRunner()
    result = runner.invoke(cli, ["posttree", "dating", "hessian", "--help"])
    assert result.exit_code == 0
    assert "--matrix" in result.output
    assert "--rooted-tree" in result.output
    assert "--seq-type" in result.output
    assert "--partitions" in result.output


def test_mcmc_help():
    runner = CliRunner()
    result = runner.invoke(cli, ["posttree", "dating", "mcmc", "--help"])
    assert result.exit_code == 0
    assert "--hessian-dir" in result.output
    assert "--ctl" in result.output
    assert "--clock" in result.output
    assert "--burnin" in result.output
    assert "--nsamples" in result.output


def test_mcmc_ctl_conflict_rejected(tmp_path):
    hdir = tmp_path / "hessian"
    hdir.mkdir()
    for f in ("iqtree.dummy.phy", "iqtree.rooted.nwk", "iqtree.mcmctree.hessian"):
        (hdir / f).write_text("x")
    ctl = tmp_path / "mine.ctl"
    ctl.write_text("seed = -1\n")
    runner = CliRunner()
    result = runner.invoke(cli, [
        "posttree", "dating", "mcmc",
        "--hessian-dir", str(hdir),
        "--ctl", str(ctl),
        "--burnin", "50000",
    ])
    assert result.exit_code != 0
    assert "--ctl" in result.output or "mutually exclusive" in result.output.lower()


def test_mcmc_dry_run_with_ctl_uses_user_file(tmp_path):
    hdir = tmp_path / "hessian"
    hdir.mkdir()
    for f in ("iqtree.dummy.phy", "iqtree.rooted.nwk", "iqtree.mcmctree.hessian"):
        (hdir / f).write_text("x")
    ctl = tmp_path / "mine.ctl"
    ctl.write_text("      seed = -1\n      model = 7\n")
    runner = CliRunner()
    result = runner.invoke(cli, [
        "posttree", "dating", "mcmc",
        "--hessian-dir", str(hdir),
        "--ctl", str(ctl),
        "--dry-run",
        "-o", str(tmp_path / "out"),
    ])
    assert result.exit_code == 0
    assert "model = 7" in result.output


def test_hessian_missing_matrix(tmp_path):
    runner = CliRunner()
    result = runner.invoke(cli, [
        "posttree", "dating", "hessian",
        "--matrix", str(tmp_path / "nope.fa"),
        "--rooted-tree", str(tmp_path / "t.nwk"),
        "-o", str(tmp_path / "out"),
    ])
    assert result.exit_code != 0


def test_mcmc_missing_hessian_dir(tmp_path):
    runner = CliRunner()
    result = runner.invoke(cli, [
        "posttree", "dating", "mcmc",
        "--hessian-dir", str(tmp_path / "nope"),
        "-o", str(tmp_path / "out"),
    ])
    assert result.exit_code != 0


def test_hessian_dry_run(tmp_path):
    matrix = tmp_path / "m.fa"
    matrix.write_text(">sp1\nMKTVFLGEI\n>sp2\nMLTVFLGEI\n")
    tree = tmp_path / "t.nwk"
    tree.write_text("(sp1,(sp2,sp3))'<4.2';\n")
    runner = CliRunner()
    result = runner.invoke(cli, [
        "posttree", "dating", "hessian",
        "--matrix", str(matrix),
        "--rooted-tree", str(tree),
        "--dry-run",
        "-o", str(tmp_path / "out"),
    ])
    assert result.exit_code == 0
    assert "iqtree" in result.output.lower()


def test_hessian_auto_seq_type_accepts_phylip_dry_run(tmp_path):
    matrix = tmp_path / "m.phy"
    matrix.write_text("2 9\nsp1  MKTVFLGEI\nsp2  MLTVFLGEI\n")
    tree = tmp_path / "t.nwk"
    tree.write_text("(sp1,sp2)'<4.2';\n")
    result = CliRunner().invoke(cli, [
        "posttree", "dating", "hessian",
        "--matrix", str(matrix),
        "--rooted-tree", str(tree),
        "--dry-run",
        "-o", str(tmp_path / "out"),
    ])
    assert result.exit_code == 0
    assert "LG+F+G4" in result.output


def test_hessian_auto_seq_type_accepts_nexus_dry_run(tmp_path):
    matrix = tmp_path / "m.nex"
    matrix.write_text(
        "#NEXUS\n"
        "begin data;\n"
        "dimensions ntax=2 nchar=8;\n"
        "format datatype=dna gap=- missing=?;\n"
        "matrix\n"
        "sp1 ACGTACGT\n"
        "sp2 ACGTACGT\n"
        ";\n"
        "end;\n"
    )
    tree = tmp_path / "t.nwk"
    tree.write_text("(sp1,sp2)'<4.2';\n")
    result = CliRunner().invoke(cli, [
        "posttree", "dating", "hessian",
        "--matrix", str(matrix),
        "--rooted-tree", str(tree),
        "--dry-run",
        "-o", str(tmp_path / "out"),
    ])
    assert result.exit_code == 0
    assert "GTR+G4" in result.output


def test_hessian_iqtree_missing_writes_env_result_json(tmp_path, monkeypatch):
    matrix = tmp_path / "m.fa"
    matrix.write_text(">sp1\nMKTV\n>sp2\nMLTV\n")
    tree = tmp_path / "t.nwk"
    tree.write_text("(sp1,sp2)'<4.2';\n")
    out = tmp_path / "out"

    def _missing_iqtree(*args, **kwargs):
        raise FileNotFoundError("iqtree3 not found")

    monkeypatch.setattr(
        "phyloai.posttree.dating_hessian._resolve_iqtree_path",
        _missing_iqtree,
    )
    result = CliRunner().invoke(cli, [
        "posttree", "dating", "hessian",
        "--matrix", str(matrix),
        "--rooted-tree", str(tree),
        "-o", str(out),
    ])
    assert result.exit_code == 3
    payload = json.loads((out / "result.json").read_text())
    assert payload["status"] == "error"
    assert payload["error_category"] == "env"
    assert "iqtree3 not found" in payload["error"]
