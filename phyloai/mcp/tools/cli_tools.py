"""MCP wrappers around PhyloAI CLI commands."""

from __future__ import annotations

import json
import subprocess
from pathlib import Path
from typing import Any, Awaitable, Callable

from phyloai.mcp.job import launch_cli
from phyloai.mcp.schema_gen import build_mcp_tool, walk_click_tree
from phyloai.mcp.tools.stubs import STUB_TOOLS, handle_stub

Handler = Callable[..., Awaitable[str]]

_CLICK_DEFAULT_OUTPUT_DIR_TOOLS = frozenset(
    {"tree_bi_pb", "tree_bi_bpcomp", "tree_bi_tracecomp", "tree_bi_readpb"}
)


def make_tool_handlers() -> dict[str, Handler]:
    """Build async handlers for CLI-backed tools and unavailable stubs."""
    from phyloai.cli.main import cli

    descriptors = walk_click_tree(cli)
    handlers = {d["tool_name"]: _make_launch_handler(d) for d in descriptors}
    for stub in STUB_TOOLS:
        if stub["name"] not in handlers:
            handlers[stub["name"]] = _make_stub_handler(stub["name"])
    return handlers


def get_tool_definitions() -> dict[str, dict]:
    """Return tool definitions keyed by name."""
    from phyloai.cli.main import cli

    definitions = {d["tool_name"]: build_mcp_tool(d) for d in walk_click_tree(cli)}
    for stub in STUB_TOOLS:
        definitions.setdefault(stub["name"], stub)
    return definitions


def _make_stub_handler(tool_name: str) -> Handler:
    async def handler(**_: Any) -> str:
        return json.dumps(handle_stub(tool_name))

    return handler


def _make_launch_handler(descriptor: dict[str, Any]) -> Handler:
    tool_name = descriptor["tool_name"]

    async def handler(**kwargs: Any) -> str:
        if tool_name == "doctor":
            return _run_sync(["phyloai", "doctor", "--output-format", "json"], timeout=30)
        if tool_name == "report":
            run_dir = kwargs.get("run_dir") or kwargs.get("run-dir")
            if run_dir is None:
                return json.dumps({"status": "error", "message": "run_dir is required"})
            result = _run_sync(["phyloai", "report", "--run-dir", str(run_dir), "--overwrite"], timeout=300)
            report_path = Path(run_dir).resolve() / "report" / "report.json"
            if report_path.exists():
                with open(report_path) as fh:
                    return json.dumps(json.load(fh))
            return result

        output_dir = kwargs.get("output_dir") or kwargs.get("output-dir")
        if output_dir is None and tool_name in _CLICK_DEFAULT_OUTPUT_DIR_TOOLS:
            output_dir = next(param.default for param in descriptor["click_command"].params if param.name == "output_dir")
        if output_dir is None:
            return json.dumps({"status": "error", "message": f"output_dir is required for {tool_name}"})
        try:
            result_dir, pid = launch_cli(descriptor, kwargs, Path(output_dir))
        except ValueError as exc:
            return json.dumps({"status": "error", "message": str(exc)})
        return json.dumps({"status": "launched", "output_dir": str(result_dir), "pid": pid, "message": f"Track progress with check_status('{result_dir}')."})

    return handler


def _run_sync(argv: list[str], *, timeout: int) -> str:
    result = subprocess.run(argv, capture_output=True, text=True, timeout=timeout)
    if result.returncode != 0:
        return json.dumps({"status": "error", "message": result.stderr[:1000], "stdout": result.stdout})
    return result.stdout or json.dumps({"status": "success"})
