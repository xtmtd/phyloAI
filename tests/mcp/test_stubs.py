from __future__ import annotations

from phyloai.cli.main import cli
from phyloai.mcp.schema_gen import walk_click_tree
from phyloai.mcp.tools.stubs import STUB_TOOL_NAMES, STUB_TOOLS, handle_stub


def test_stub_tools_are_valid_and_do_not_overlap_real_cli() -> None:
    dynamic = {d["tool_name"] for d in walk_click_tree(cli)}

    assert STUB_TOOL_NAMES
    assert not (STUB_TOOL_NAMES & dynamic)
    assert {t["name"] for t in STUB_TOOLS} == STUB_TOOL_NAMES
    assert handle_stub(next(iter(STUB_TOOL_NAMES)))["status"] == "not_implemented"
    assert handle_stub("pretree_align") is None
