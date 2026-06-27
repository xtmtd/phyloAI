"""Read-only MCP utility tools."""

from __future__ import annotations

import json
import os
from pathlib import Path

from phyloai.mcp.job import read_job_json
from phyloai.mcp.schema_gen import build_mcp_tool, walk_click_tree

_schema_cache: dict[str, dict] | None = None


def check_status(output_dir: str) -> dict:
    """Inspect an output directory and return job state."""
    od = Path(output_dir).resolve()
    result = _read_json(od / "result.json")
    if result is not None:
        return {"status": result.get("status", "unknown"), "output_dir": str(od), "result": result}

    checkpoint = _read_json(od / "checkpoint.json")
    job = read_job_json(od)
    if job is None:
        return {"status": "not_started", "output_dir": str(od)}

    response = {"output_dir": str(od)}
    if checkpoint is not None:
        response["checkpoint"] = checkpoint
    pid = job.get("pid")
    if pid is not None:
        try:
            os.kill(int(pid), 0)
            response["status"] = "running"
            return response
        except (OSError, ValueError):
            pass
    response.update({"status": "unknown", "message": "Process exited but result.json not found. Check logs/ for tool stderr."})
    return response


def read_result(output_dir: str) -> dict:
    """Read `<output_dir>/result.json`."""
    path = Path(output_dir).resolve() / "result.json"
    payload = _read_json(path)
    if payload is None:
        return {"status": "error", "message": f"result.json not found at {path}"}
    return payload


def read_report(run_dir: str) -> dict:
    """Read `<run_dir>/report/report.json`."""
    path = Path(run_dir).resolve() / "report" / "report.json"
    payload = _read_json(path)
    if payload is None:
        return {"status": "error", "message": f"report.json not found at {path}. Run 'phyloai report --run-dir {Path(run_dir).resolve()}' to generate it."}
    return payload


def get_command_schema(command_name: str) -> dict:
    """Return runtime MCP schema for a Click command tool name."""
    cache = _load_schema_cache()
    if command_name in cache:
        return cache[command_name]
    return {"error": f"Unknown command: {command_name}", "available": sorted(cache)}


def _read_json(path: Path) -> dict | None:
    if not path.exists():
        return None
    try:
        with open(path) as fh:
            return json.load(fh)
    except (OSError, json.JSONDecodeError):
        return None


def _load_schema_cache() -> dict[str, dict]:
    global _schema_cache
    if _schema_cache is None:
        from phyloai.cli.main import cli

        _schema_cache = {d["tool_name"]: build_mcp_tool(d) for d in walk_click_tree(cli)}
    return _schema_cache
