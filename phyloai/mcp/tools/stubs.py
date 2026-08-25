"""Stub MCP tools for public commands not implemented in this package.

The registry is intentionally empty: every current public CLI command is
implemented and discovered dynamically from the Click tree. ``handle_stub``
is kept so ``cli_tools.py`` can route unknown tool names safely.
"""

from __future__ import annotations

STUB_TOOL_NAMES: frozenset[str] = frozenset()

_DESCRIPTIONS: dict[str, str] = {}

STUB_TOOLS: list[dict] = []


def handle_stub(tool_name: str) -> dict | None:
    """Return a not-implemented response for known stub tools."""
    if tool_name not in STUB_TOOL_NAMES:
        return None
    return {"status": "not_implemented", "message": f"The '{tool_name}' command is not yet available in the installed version."}
