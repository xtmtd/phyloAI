from __future__ import annotations

import json

from click.testing import CliRunner

from phyloai.cli.main import cli
from phyloai.mcp.schema_gen import walk_click_tree
from phyloai.mcp.tools.stubs import STUB_TOOL_NAMES


def _inputs(tmp_path):
    frequencies = " ".join(["0.05"] * 20)
    freq = tmp_path / "chain.sitefreq"
    freq.write_text(f"1 {frequencies}\n")
    first = tmp_path / "first.csv"
    first.write_text("site,lnL_Tree1,lnL_Tree2\n1,-4,-3\n")
    second = tmp_path / "second.csv"
    second.write_text("site,lnL_Tree1,lnL_Tree2\n1,-5,-4\n")
    return freq, first, second


def test_cca_help_and_generated_mcp_leaf() -> None:
    result = CliRunner().invoke(cli, ["posttree", "syserror", "cca", "--help"])

    assert result.exit_code == 0
    for option in ("--site-freq", "--site-lnl1", "--site-lnl2", "--model1-name", "--fig-width", "--dry-run"):
        assert option in result.output
    assert "site, lnL_Tree1, and lnL_Tree2" in result.output
    assert "posttree_syserror_cca" in {item["tool_name"] for item in walk_click_tree(cli)}
    assert "posttree_syserror_cca" not in STUB_TOOL_NAMES


def test_cca_cli_writes_outputs_and_error_result(tmp_path) -> None:
    freq, first, second = _inputs(tmp_path)
    output = tmp_path / "out"
    result = CliRunner().invoke(cli, [
        "posttree", "syserror", "cca", "--site-freq", str(freq), "--site-lnl1", str(first),
        "--site-lnl2", str(second), "--model1-name", "LG", "--model2-name", "C20", "-o", str(output),
    ])

    assert result.exit_code == 0, result.output
    assert (output / "cca.csv").exists() and (output / "cca.pdf").exists()
    assert f"Output: {output / 'cca.csv'}" in result.output
    assert f"Output: {output / 'cca.pdf'}" in result.output
    assert f"Output: {output / 'result.json'}" in result.output
    assert set(json.loads((output / "result.json").read_text())["data"]["output_files"]) == {"cca_table", "cca_figure"}

    failed = CliRunner().invoke(cli, [
        "posttree", "syserror", "cca", "--site-freq", str(tmp_path / "missing"),
        "--site-lnl1", str(first), "--site-lnl2", str(second), "-o", str(tmp_path / "bad"),
    ])
    assert failed.exit_code == 1
    assert (tmp_path / "bad" / "result.json").exists()




def test_cca_cli_invalid_figure_parameter_writes_error_result(tmp_path) -> None:
    freq, first, second = _inputs(tmp_path)
    output = tmp_path / "bad-figure"

    result = CliRunner().invoke(cli, [
        "posttree", "syserror", "cca", "--site-freq", str(freq), "--site-lnl1", str(first),
        "--site-lnl2", str(second), "--fig-width", "0", "-o", str(output),
    ])

    assert result.exit_code == 1
    assert json.loads((output / "result.json").read_text())["status"] == "error"


    _, first, second = _inputs(tmp_path)
    output = tmp_path / "out"
    output.mkdir()
    keep = output / "keep.txt"
    keep.write_text("existing data")

    result = CliRunner().invoke(cli, [
        "posttree", "syserror", "cca", "--site-freq", str(tmp_path / "missing"),
        "--site-lnl1", str(first), "--site-lnl2", str(second),
        "--overwrite", "-o", str(output),
    ])

    assert result.exit_code == 1
    assert keep.read_text() == "existing data"
    error = json.loads((output / "result.json").read_text())
    assert error["status"] == "error"
    assert error["error_category"] == "input"


def test_cca_cli_dry_run_writes_nothing(tmp_path) -> None:
    freq, first, second = _inputs(tmp_path)
    output = tmp_path / "dry"

    result = CliRunner().invoke(cli, [
        "posttree", "syserror", "cca", "--site-freq", str(freq), "--site-lnl1", str(first),
        "--site-lnl2", str(second), "-o", str(output), "--dry-run",
    ])

    assert result.exit_code == 0, result.output
    assert not output.exists()
