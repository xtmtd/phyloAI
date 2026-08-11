from __future__ import annotations

import json

from click.testing import CliRunner

from phyloai.cli.main import cli
from phyloai.mcp.schema_gen import walk_click_tree
from phyloai.mcp.tools.stubs import STUB_TOOL_NAMES
from phyloai.report.collector import STEP_ORDER, parse_step_id
from phyloai.report.templates import generate_all_methods


def test_brlen_help_lists_modes_map_and_patristic_guard() -> None:
    result = CliRunner().invoke(cli, ["posttree", "syserror", "brlen", "--help"])
    assert result.exit_code == 0
    for text in ("--tree-dir", "tip-to-tip", "node-to-tip", "--map", "--max-rows", "O(n²)"):
        assert text in result.output


def test_brlen_click_leaves_replace_stub() -> None:
    names = {item["tool_name"] for item in walk_click_tree(cli)}
    assert {"posttree_syserror_brlen", "posttree_syserror_brlen_label_nodes"} <= names
    assert "posttree_syserror_brlen" not in STUB_TOOL_NAMES


def test_brlen_methods_text_is_quantitative() -> None:
    text = generate_all_methods(
        "posttree.syserror.brlen",
        {"map": "nodes.map"},
        {"n_trees": 20, "modes": ["terminal"], "summary": {"terminal": {"mean": 0.1, "sd": 0.02}}},
        {},
    )
    assert "20" in text and "map" in text.lower() and "0.1000" in text


def test_label_nodes_has_distinct_report_step_without_methods_text() -> None:
    command = "phyloai posttree syserror brlen label-nodes --tree species.nwk"
    assert parse_step_id(command) == "posttree.syserror.brlen.label-nodes"
    assert "posttree.syserror.brlen.label-nodes" in STEP_ORDER
    assert generate_all_methods("posttree.syserror.brlen.label-nodes", {}, {}, {}) == ""


def test_brlen_cli_runs_end_to_end(tmp_path) -> None:
    tree_file = tmp_path / "one.tre"
    tree_file.write_text("((A:1,B:2):3,C:4);")
    result = CliRunner().invoke(cli, [
        "posttree", "syserror", "brlen", "--tree", str(tree_file),
        "--mode", "terminal", "-o", str(tmp_path / "out"),
    ])
    assert result.exit_code == 0
    assert (tmp_path / "out" / "tables" / "terminal.csv").exists()
    assert (tmp_path / "out" / "result.json").exists()


def test_brlen_cli_rejects_missing_mode(tmp_path) -> None:
    tree_file = tmp_path / "one.tre"
    tree_file.write_text("((A:1,B:2):3,C:4);")
    result = CliRunner().invoke(cli, [
        "posttree", "syserror", "brlen", "--tree", str(tree_file),
        "-o", str(tmp_path / "out"),
    ])
    assert result.exit_code != 0


def test_brlen_cli_endpoint_failure_writes_error_result_json(tmp_path) -> None:
    tree_file = tmp_path / "one.tre"
    tree_file.write_text("((A:1,B:2):3,C:4);")
    result = CliRunner().invoke(cli, [
        "posttree", "syserror", "brlen", "--tree", str(tree_file),
        "--mode", "node-to-tip", "--node1", "Missing", "-o", str(tmp_path / "out"),
    ])
    assert result.exit_code == 1
    error_json = json.loads((tmp_path / "out" / "result.json").read_text())
    assert error_json["status"] == "error"
    assert "Missing" in error_json["error"]


def test_brlen_cli_dry_run_failure_writes_nothing(tmp_path) -> None:
    tree_file = tmp_path / "one.tre"
    tree_file.write_text("((A:1,B:2):3,C:4);")
    result = CliRunner().invoke(cli, [
        "posttree", "syserror", "brlen", "--tree", str(tree_file),
        "--mode", "node-to-tip", "--node1", "Missing",
        "-o", str(tmp_path / "out"), "--dry-run",
    ])
    assert result.exit_code == 1
    assert not (tmp_path / "out").exists()


def test_brlen_cli_refused_nonempty_dir_is_not_polluted(tmp_path) -> None:
    tree_file = tmp_path / "one.tre"
    tree_file.write_text("((A:1,B:2):3,C:4);")
    out = tmp_path / "out"
    out.mkdir()
    (out / "keep.txt").write_text("existing data")
    result = CliRunner().invoke(cli, [
        "posttree", "syserror", "brlen", "--tree", str(tree_file),
        "--mode", "terminal", "-o", str(out),
    ])
    assert result.exit_code == 1
    assert (out / "keep.txt").exists()
    assert not (out / "result.json").exists()
