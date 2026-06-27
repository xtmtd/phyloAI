"""PhyloAI MCP stdio server."""

from __future__ import annotations

import json
import sys
from typing import Any

from mcp.server import Server
from mcp.server.stdio import stdio_server
from mcp.types import TextContent, Tool

from phyloai.mcp.tools.cli_tools import get_tool_definitions, make_tool_handlers
from phyloai.mcp.tools.utils import check_status, get_command_schema, read_report, read_result


def create_server() -> Server:
    """Create a stdio MCP server with PhyloAI tools registered."""
    definitions = get_tool_definitions()
    handlers = make_tool_handlers()
    utility_definitions = _utility_definitions()
    server = Server("phyloai")

    @server.list_tools()
    async def list_tools() -> list[Tool]:
        all_defs = {**definitions, **utility_definitions}
        return [
            Tool(name=d["name"], description=d.get("description", ""), inputSchema=d["inputSchema"])
            for d in all_defs.values()
        ]

    @server.call_tool()
    async def call_tool(name: str, arguments: dict[str, Any] | None) -> list[TextContent]:
        args = arguments or {}
        if name == "check_status":
            result = check_status(args["output_dir"])
        elif name == "read_result":
            result = read_result(args["output_dir"])
        elif name == "read_report":
            result = read_report(args["run_dir"])
        elif name == "get_command_schema":
            result = get_command_schema(args["command_name"])
        elif name in handlers:
            return [TextContent(type="text", text=await handlers[name](**args))]
        else:
            result = {"status": "error", "message": f"Unknown tool: {name}"}
        return [TextContent(type="text", text=json.dumps(result))]

    return server


def _utility_definitions() -> dict[str, dict]:
    return {
        "check_status": {"name": "check_status", "description": "Check a PhyloAI job status by output_dir.", "inputSchema": _schema({"output_dir": "string"}, ["output_dir"])},
        "read_result": {"name": "read_result", "description": "Read result.json from an output directory.", "inputSchema": _schema({"output_dir": "string"}, ["output_dir"])},
        "read_report": {"name": "read_report", "description": "Read report/report.json from a run directory.", "inputSchema": _schema({"run_dir": "string"}, ["run_dir"])},
        "get_command_schema": {"name": "get_command_schema", "description": "Return runtime schema for a PhyloAI tool.", "inputSchema": _schema({"command_name": "string"}, ["command_name"])},
    }


def _schema(properties: dict[str, str], required: list[str]) -> dict:
    return {"type": "object", "properties": {k: {"type": v} for k, v in properties.items()}, "required": required}


async def main() -> None:
    server = create_server()
    async with stdio_server() as streams:
        await server.run(streams[0], streams[1], server.create_initialization_options())


def entry() -> None:
    import asyncio

    try:
        asyncio.run(main())
    except KeyboardInterrupt:
        pass
    except Exception as exc:
        print(f"MCP server error: {exc}", file=sys.stderr)
        raise SystemExit(1) from exc
