from __future__ import annotations

import json
from pathlib import Path

from click.testing import CliRunner

from phyloai.cli.main import cli
from phyloai.mcp.schema_gen import walk_click_tree
from phyloai.mcp.tools.stubs import STUB_TOOL_NAMES


def _matrix(tmp_path) -> None:
    matrix = tmp_path / "matrix.fa"
    matrix.write_text(">A\nAAAA\n>B\nCCCC\n>C\nAACC\n")
    return matrix


def test_taxcomp_help_and_generated_mcp_leaf() -> None:
    result = CliRunner().invoke(cli, ["posttree", "syserror", "taxcomp", "--help"])

    assert result.exit_code == 0
    for option in ("--matrix", "--seq-type", "--table-format", "--overwrite", "--dry-run", "--quiet"):
        assert option in result.output
    assert "exploratory" in result.output.lower()
    assert "sparse-cell" in result.output.lower()
    # parameters are few, so taxcomp help must not group options under section headers
    for header in ("Required Inputs:", "Analysis Options:", "Common Options:"):
        assert header not in result.output
    # detailed-help contract (design 8.2)
    assert "Clustal is unsupported" in result.output
    assert "auto resolves AA/NT" in result.output
    assert "NT only ACGT" in result.output
    assert "phylogenetically dependent" in result.output
    assert "posttree_syserror_taxcomp" in {item["tool_name"] for item in walk_click_tree(cli)}
    assert "posttree_syserror_taxcomp" not in STUB_TOOL_NAMES


def test_taxcomp_cli_writes_tables_and_result_json(tmp_path) -> None:
    matrix = _matrix(tmp_path)
    output = tmp_path / "out"

    result = CliRunner().invoke(cli, [
        "posttree", "syserror", "taxcomp", "--matrix", str(matrix),
        "--seq-type", "NT", "-o", str(output),
    ])

    assert result.exit_code == 0, result.output
    assert (output / "overall_summary.csv").exists()
    assert (output / "taxon_summary.csv").exists()
    assert (output / "result.json").exists()
    payload = json.loads((output / "result.json").read_text())
    assert payload["status"] == "success"
    assert payload["error_category"] is None
    assert payload["key_results"]["seq_type"] == "NT"
    assert payload["key_results"]["x2"] == 8.0


def test_taxcomp_cli_dry_run_writes_nothing(tmp_path) -> None:
    matrix = _matrix(tmp_path)
    output = tmp_path / "dry"

    result = CliRunner().invoke(cli, [
        "posttree", "syserror", "taxcomp", "--matrix", str(matrix),
        "--seq-type", "NT", "-o", str(output), "--dry-run",
    ])

    assert result.exit_code == 0, result.output
    assert not output.exists()
    payload = json.loads(result.output)
    assert payload["data"]["output_files"] == {}


def test_taxcomp_cli_invalid_matrix_writes_error_result(tmp_path) -> None:
    output = tmp_path / "bad"

    result = CliRunner().invoke(cli, [
        "posttree", "syserror", "taxcomp", "--matrix", str(tmp_path / "missing.fa"),
        "-o", str(output),
    ])

    assert result.exit_code == 1
    payload = json.loads((output / "result.json").read_text())
    assert payload["status"] == "error"
    assert payload["error_category"] == "input"


def test_taxcomp_command_docs_cover_interpretation_boundaries() -> None:
    root = Path(__file__).resolve().parents[2]
    en = (root / "docs/commands/posttree-syserror-taxcomp.md").read_text()
    zh = (root / "docs/commands/posttree-syserror-taxcomp.zh.md").read_text()

    for text in (en, zh):
        assert "phyloai posttree syserror taxcomp" in text
        for token in ("p_holm", "sparse_count_check", "comp_max", "comp_mean", "Dayhoff-6"):
            assert token in text
    assert "automatic" in en
    assert "自动" in zh
