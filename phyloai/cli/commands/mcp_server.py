"""phyloai mcp-server command."""

from __future__ import annotations

import click


@click.command(name="mcp-server")
def mcp_server() -> None:
    """Start the PhyloAI MCP server over stdio."""
    from phyloai.mcp.server import entry

    entry()
