from __future__ import annotations

from pathlib import Path

import click

from phyloai.cli.main import cli
from phyloai.mcp.schema_gen import build_mcp_tool, click_param_to_json_schema, walk_click_tree


def test_click_param_to_json_schema_handles_common_click_types() -> None:
    assert click_param_to_json_schema(click.Option(["--name"], default="x", help="Name")) == {
        "description": "Name",
        "type": "string",
        "default": "x",
    }

    choice = click_param_to_json_schema(click.Option(["--method"], type=click.Choice(["linsi", "auto"]), default="linsi"))
    assert choice["type"] == "string"
    assert choice["enum"] == ["linsi", "auto"]
    assert choice["default"] == "linsi"

    path = click_param_to_json_schema(click.Option(["--matrix"], type=click.Path(path_type=Path), required=True))
    assert path["type"] == "string"
    assert path["format"] == "path"

    assert click_param_to_json_schema(click.Option(["--threads"], type=int, default=4))["type"] == "integer"
    assert click_param_to_json_schema(click.Option(["--overwrite"], is_flag=True, default=False))["type"] == "boolean"


def test_walk_click_tree_finds_leaf_commands_only() -> None:
    tool_names = {d["tool_name"] for d in walk_click_tree(cli)}

    assert "pretree_convert" in tool_names
    assert "pretree_align" in tool_names
    assert "tree_ml_iqtree" in tool_names
    assert "doctor" in tool_names
    assert "report" in tool_names
    assert "posttree_dating_hessian" in tool_names
    assert "posttree_dating_mcmc" in tool_names
    assert "pretree" not in tool_names
    assert "tree" not in tool_names
    assert "posttree_dating" not in tool_names


def test_walk_click_tree_excludes_non_analysis_commands() -> None:
    tool_names = {d["tool_name"] for d in walk_click_tree(cli)}

    assert "mcp-server" not in tool_names
    assert "mcp_server" not in tool_names
    assert "completion_bash" not in tool_names
    assert "completion_zsh" not in tool_names
    assert "completion_fish" not in tool_names


def test_build_mcp_tool_doctor_exposes_runtime_schema() -> None:
    doctor_descriptor = next(d for d in walk_click_tree(cli) if d["tool_name"] == "doctor")
    tool_def = build_mcp_tool(doctor_descriptor)

    assert tool_def["name"] == "doctor"
    assert tool_def["inputSchema"]["type"] == "object"
    assert "output_format" in tool_def["inputSchema"]["properties"]


def test_walk_click_tree_finds_pretree_concat_jackknife() -> None:
    descriptor = next(d for d in walk_click_tree(cli) if d["tool_name"] == "pretree_concat_jackknife")
    tool_def = build_mcp_tool(descriptor)

    assert tool_def["name"] == "pretree_concat_jackknife"
    props = tool_def["inputSchema"]["properties"]
    assert props["replicates"]["default"] == 100
    assert props["target_length"]["default"] == 50000
    assert props["seed"]["default"] == 42
    assert props["table_format"]["enum"] == ["csv", "tsv"]


def test_tree_ml_iqtree_schema_exposes_site_freq_file() -> None:
    descriptor = next(d for d in walk_click_tree(cli) if d["tool_name"] == "tree_ml_iqtree")
    props = build_mcp_tool(descriptor)["inputSchema"]["properties"]

    assert props["site_freq_file"]["type"] == "string"
    assert props["site_freq_file"]["format"] == "path"


def test_walk_click_tree_finds_signal_subcommands() -> None:
    tool_names = {d["tool_name"] for d in walk_click_tree(cli)}

    assert "posttree_signal_lnl" in tool_names
    assert "posttree_signal_consistent" in tool_names
    assert "posttree_signal_fclm" in tool_names
    assert "posttree_signal" not in tool_names


def test_signal_schemas_expose_required_params() -> None:
    for name in ("posttree_signal_lnl", "posttree_signal_consistent", "posttree_signal_fclm"):
        descriptor = next(d for d in walk_click_tree(cli) if d["tool_name"] == name)
        tool_def = build_mcp_tool(descriptor)
        props = tool_def["inputSchema"]["properties"]
        required = tool_def["inputSchema"].get("required", [])

        assert "matrix" in props
        assert "matrix" in required
        assert "output_dir" in props
        assert "overwrite" in props
        assert "dry_run" in props
        assert "quiet" in props

    lnl = next(d for d in walk_click_tree(cli) if d["tool_name"] == "posttree_signal_lnl")
    lnl_tool = build_mcp_tool(lnl)
    assert "candidate_trees_raw" in lnl_tool["inputSchema"]["required"]

    fclm = next(d for d in walk_click_tree(cli) if d["tool_name"] == "posttree_signal_fclm")
    fclm_tool = build_mcp_tool(fclm)
    assert "taxset_csv" in fclm_tool["inputSchema"]["required"]


def test_adequacy_mcp_tool_is_generated_from_click() -> None:
    descriptor = next(
        item for item in walk_click_tree(cli)
        if item["tool_name"] == "posttree_simulate_adequacy"
    )
    tool = build_mcp_tool(descriptor)

    assert tool["inputSchema"]["required"] == ["original_msa", "simulated_dir"]
    assert "threads" in tool["inputSchema"]["properties"]


def test_brlen_main_and_label_nodes_schemas_are_generated() -> None:
    tool_names = {d["tool_name"] for d in walk_click_tree(cli)}
    assert "posttree_syserror_brlen" in tool_names
    assert "posttree_syserror_brlen_label_nodes" in tool_names

    brlen = next(d for d in walk_click_tree(cli) if d["tool_name"] == "posttree_syserror_brlen")
    props = build_mcp_tool(brlen)["inputSchema"]["properties"]
    for key in ("tree", "tree_dir", "mode", "map", "node1", "node2", "tip1", "tip2",
                "table_format", "threads", "max_rows", "overwrite", "dry_run", "quiet"):
        assert key in props
    assert props["table_format"]["enum"] == ["csv", "tsv"]
    assert props["max_rows"]["default"] == 5000000

    label = next(d for d in walk_click_tree(cli) if d["tool_name"] == "posttree_syserror_brlen_label_nodes")
    label_props = build_mcp_tool(label)["inputSchema"]
    assert label_props["required"] == ["tree"]
    assert set(label_props["properties"]) == {"tree", "output_dir", "overwrite", "quiet"}


def test_taxcomp_mcp_tool_is_generated_from_click() -> None:
    descriptor = next(
        item for item in walk_click_tree(cli)
        if item["tool_name"] == "posttree_syserror_taxcomp"
    )
    tool = build_mcp_tool(descriptor)

    assert tool["inputSchema"]["required"] == ["matrix"]
    props = tool["inputSchema"]["properties"]
    assert props["seq_type"]["enum"] == ["AA", "NT", "auto"]
    assert props["table_format"]["enum"] == ["csv", "tsv"]
