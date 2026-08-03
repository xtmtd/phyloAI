from __future__ import annotations

import json

import pytest

from phyloai.cli.main import cli
from phyloai.mcp.schema_gen import walk_click_tree
from phyloai.mcp.tools import cli_tools
from phyloai.mcp.tools.cli_tools import get_tool_definitions, make_tool_handlers
from phyloai.mcp.tools.stubs import STUB_TOOL_NAMES


def test_make_tool_handlers_includes_dynamic_and_non_overlapping_stubs() -> None:
    handlers = make_tool_handlers()
    dynamic = {d["tool_name"] for d in walk_click_tree(cli)}

    assert "doctor" in handlers
    assert "read_result" not in handlers
    assert not (dynamic & STUB_TOOL_NAMES)
    assert STUB_TOOL_NAMES <= set(handlers)


@pytest.mark.anyio
async def test_stub_handler_returns_json() -> None:
    handlers = make_tool_handlers()
    result = await handlers[next(iter(STUB_TOOL_NAMES))]()

    assert json.loads(result)["status"] == "not_implemented"


def test_tree_bi_tools_replace_legacy_tool() -> None:
    definitions = get_tool_definitions()

    assert {
        "tree_bi_pb",
        "tree_bi_bpcomp",
        "tree_bi_tracecomp",
        "tree_bi_readpb",
    } <= set(definitions)
    assert "tree_bi" not in definitions
    assert "iqtree_rate" not in definitions["tree_bi_readpb"]["inputSchema"]["properties"]


@pytest.mark.anyio
@pytest.mark.parametrize(
    "tool_name",
    ["tree_bi_pb", "tree_bi_bpcomp", "tree_bi_tracecomp", "tree_bi_readpb"],
)
async def test_tree_bi_handlers_accept_default_output_dir(monkeypatch: pytest.MonkeyPatch, tool_name: str) -> None:
    monkeypatch.setattr(cli_tools, "launch_cli", lambda descriptor, params, output_dir: (output_dir, 123))

    result = json.loads(await make_tool_handlers()[tool_name]())

    assert result["status"] == "launched"


@pytest.mark.anyio
@pytest.mark.parametrize(
    "tool_name",
    [
        "posttree_simulate_alisim_params",
        "posttree_simulate_alisim_iqtree",
        "posttree_simulate_alisim_transfergaps",
    ],
)
async def test_alisim_handlers_accept_default_output_dir(
    monkeypatch: pytest.MonkeyPatch, tool_name: str,
) -> None:
    monkeypatch.setattr(cli_tools, "launch_cli", lambda descriptor, params, output_dir: (output_dir, 123))

    result = json.loads(await make_tool_handlers()[tool_name]())

    assert result["status"] == "launched"
