from __future__ import annotations

from phyloai.cli.main import cli
from phyloai.mcp.schema_gen import walk_click_tree
from phyloai.mcp.tools.stubs import STUB_TOOL_NAMES, STUB_TOOLS, handle_stub


def test_stub_registry_is_empty_and_sites_is_not_a_stub() -> None:
    dynamic = {d["tool_name"] for d in walk_click_tree(cli)}

    assert STUB_TOOL_NAMES == frozenset()
    assert STUB_TOOLS == []
    assert "posttree_syserror_sites" not in STUB_TOOL_NAMES
    assert "posttree_syserror_sites" not in dynamic
    assert "posttree_syserror_taxcomp" in dynamic
    assert handle_stub("posttree_syserror_sites") is None
    assert handle_stub("pretree_align") is None


def test_brlen_is_no_longer_a_stub() -> None:
    assert "posttree_syserror_brlen" not in STUB_TOOL_NAMES
    assert "posttree_syserror_brlen" in {d["tool_name"] for d in walk_click_tree(cli)}
