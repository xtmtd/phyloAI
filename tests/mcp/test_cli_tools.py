from __future__ import annotations

import json

import pytest

from phyloai.cli.main import cli
from phyloai.mcp.schema_gen import walk_click_tree
from phyloai.mcp.tools.cli_tools import make_tool_handlers
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
