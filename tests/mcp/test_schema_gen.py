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


def test_build_mcp_tool_doctor_exposes_runtime_schema() -> None:
    doctor_descriptor = next(d for d in walk_click_tree(cli) if d["tool_name"] == "doctor")
    tool_def = build_mcp_tool(doctor_descriptor)

    assert tool_def["name"] == "doctor"
    assert tool_def["inputSchema"]["type"] == "object"
    assert "output_format" in tool_def["inputSchema"]["properties"]
