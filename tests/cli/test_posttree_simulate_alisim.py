"""CLI tests for phyloai posttree simulate alisim."""
from __future__ import annotations

import json
from pathlib import Path

from click.testing import CliRunner

from phyloai.cli.main import cli


def _write_report(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        "Input data: amino-acid\nTo simulate an alignment of the same\n"
        'iqtree3 --alisim sim -m "LG+G4{.6}" --length 100\n'
    )


def test_simulate_alisim_commands_are_registered() -> None:
    result = CliRunner().invoke(cli, ["posttree", "simulate", "alisim", "--help"])
    assert result.exit_code == 0
    assert {"params", "iqtree", "transfergaps"} <= set(result.output.split())


def test_simulate_future_commands_registered() -> None:
    result = CliRunner().invoke(cli, ["posttree", "simulate", "--help"])
    assert result.exit_code == 0
    assert {"alisim", "adequacy", "phybase"} <= set(result.output.split())


def test_adequacy_and_phybase_report_not_implemented() -> None:
    for name in ("adequacy", "phybase"):
        result = CliRunner().invoke(cli, ["posttree", "simulate", name])
        assert result.exit_code == 0
        assert "not yet implemented" in result.output


def test_alisim_params_dry_run_writes_nothing(tmp_path: Path) -> None:
    iqtree_dir = tmp_path / "reports"
    tree_dir = tmp_path / "trees"
    tree_dir.mkdir()
    _write_report(iqtree_dir / "g1.iqtree")
    (tree_dir / "g1.treefile").write_text("(A,B);\n")
    output_dir = tmp_path / "out"

    result = CliRunner().invoke(cli, [
        "posttree", "simulate", "alisim", "params",
        "--iqtree-dir", str(iqtree_dir), "--tree-dir", str(tree_dir),
        "--output-dir", str(output_dir), "--dry-run", "--quiet",
    ])
    assert result.exit_code == 0
    assert not output_dir.exists()


def test_alisim_params_no_report_files_fails(tmp_path: Path) -> None:
    iqtree_dir = tmp_path / "reports"
    iqtree_dir.mkdir()
    tree_dir = tmp_path / "trees"
    tree_dir.mkdir()
    output_dir = tmp_path / "out"

    result = CliRunner().invoke(cli, [
        "posttree", "simulate", "alisim", "params",
        "--iqtree-dir", str(iqtree_dir), "--tree-dir", str(tree_dir),
        "--output-dir", str(output_dir),
    ])
    assert result.exit_code == 1
    assert "no .iqtree files" in result.output
    payload = json.loads((output_dir / "result.json").read_text())
    assert payload["status"] == "error"


def test_alisim_iqtree_single_dry_run(tmp_path: Path) -> None:
    ref_tree = tmp_path / "ref.nwk"
    ref_tree.write_text("(A,B);\n")
    result = CliRunner().invoke(cli, [
        "posttree", "simulate", "alisim", "iqtree",
        "--ref-tree", str(ref_tree), "--model", "LG+G4", "--seq-type", "AA",
        "--length", "100", "--dry-run", "--quiet",
    ])
    assert result.exit_code == 0
    assert "--alisim" in result.output or result.output == ""


def test_alisim_iqtree_batch_dry_run(tmp_path: Path) -> None:
    table = tmp_path / "params.tsv"
    table.write_text(
        "id\tseqtype\tlength\tsubs_model\tsubs_rate\tfreq\tprop_inv\t"
        "rate_heterogeneity\trate_categories\trate_param\ttree_path\n"
        "g1\tAA\t100\tLG\t\t\t\tG\t4\t0.6\t/t/g1.tre\n"
    )
    result = CliRunner().invoke(cli, [
        "posttree", "simulate", "alisim", "iqtree",
        "--model-params", str(table), "--strategy", "complete",
        "--num-simulations", "2", "--seed", "1", "--dry-run", "--quiet",
    ])
    assert result.exit_code == 0


def test_alisim_iqtree_rejects_overwrite_and_resume(tmp_path: Path) -> None:
    table = tmp_path / "params.tsv"
    table.write_text(
        "id\tseqtype\tlength\tsubs_model\tsubs_rate\tfreq\tprop_inv\t"
        "rate_heterogeneity\trate_categories\trate_param\ttree_path\n"
        "g1\tAA\t100\tLG\t\t\t\tG\t4\t0.6\t/t/g1.tre\n"
    )
    result = CliRunner().invoke(cli, [
        "posttree", "simulate", "alisim", "iqtree",
        "--model-params", str(table), "--strategy", "complete",
        "--num-simulations", "2", "--overwrite", "--resume", "--dry-run",
    ])
    assert result.exit_code == 1
    assert "mutually exclusive" in result.output


def test_alisim_transfergaps_dry_run(tmp_path: Path) -> None:
    original = tmp_path / "original.fa"
    original.write_text(">A\nAC-GT-\n>B\nACG-TA\n")
    simulated = tmp_path / "sim001.fa"
    simulated.write_text(">A\nACGTGT\n>B\nACGTAC\n")
    output_dir = tmp_path / "out"
    result = CliRunner().invoke(cli, [
        "posttree", "simulate", "alisim", "transfergaps",
        "--original-msa", str(original), "--simulated-msa", str(simulated),
        "--output-dir", str(output_dir), "--dry-run",
    ])
    assert result.exit_code == 0
    assert "Sequences: 2" in result.output
    assert not output_dir.exists()


def test_alisim_transfergaps_error_writes_result_json(tmp_path: Path) -> None:
    original = tmp_path / "original.fa"
    original.write_text(">A\nAC-GT-\n")
    missing = tmp_path / "nope.fa"
    output_dir = tmp_path / "out"
    result = CliRunner().invoke(cli, [
        "posttree", "simulate", "alisim", "transfergaps",
        "--original-msa", str(original), "--simulated-msa", str(missing),
        "--output-dir", str(output_dir),
    ])
    assert result.exit_code == 1
    assert "does not exist" in result.output
    payload = json.loads((output_dir / "result.json").read_text())
    assert payload["status"] == "error"
    assert payload["error_category"] == "input"
