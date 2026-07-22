"""Stub MCP tools for public commands not implemented in this package."""

from __future__ import annotations

STUB_TOOL_NAMES: frozenset[str] = frozenset(
    {
        "posttree_simulate",
        "posttree_syserror_brlen",
        "posttree_syserror_cca",
        "posttree_syserror_sites",
    }
)

_DESCRIPTIONS = {
    "posttree_simulate": "AliSim simulation and gene-jackknife resampling (not yet available).",
    "posttree_syserror_brlen": "Systematic error diagnosis: branch-length screen (not yet available).",
    "posttree_syserror_cca": "Systematic error diagnosis: composition analysis (not yet available).",
    "posttree_syserror_sites": "Systematic error diagnosis: site-wise analysis (not yet available).",
}

STUB_TOOLS: list[dict] = [
    {"name": name, "description": _DESCRIPTIONS[name], "inputSchema": {"type": "object", "properties": {}, "required": []}}
    for name in sorted(STUB_TOOL_NAMES)
]


def handle_stub(tool_name: str) -> dict | None:
    """Return a not-implemented response for known stub tools."""
    if tool_name not in STUB_TOOL_NAMES:
        return None
    return {"status": "not_implemented", "message": f"The '{tool_name}' command is not yet available in the installed version."}
