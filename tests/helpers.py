"""Shared test helpers for result.json structural validation.

Per JSON Output Standard (docs/superpowers/specs/2026-06-21-phyloai-json-output-standard.md)
Section 8: Testing Assertions.
"""

from __future__ import annotations

import shlex
from typing import Any


def validate_result_json(payload: dict[str, Any]) -> None:
    """Validate a result.json payload against the structural spec.

    Raises AssertionError on non-compliance.
    """
    assert isinstance(payload, dict), "result.json must be a dict"

    assert payload.get("status") in ("success", "error"), \
        f"status must be 'success' or 'error', got {payload.get('status')!r}"

    assert "wall_time" in payload, "wall_time field required"
    assert isinstance(payload["wall_time"], (int, float)), \
        "wall_time must be numeric"

    assert "tool_versions" in payload, "tool_versions field required"
    assert isinstance(payload["tool_versions"], dict), \
        "tool_versions must be a dict"

    assert "params" in payload, "params field required"
    assert isinstance(payload["params"], dict), "params must be a dict"

    assert "key_results" in payload, "key_results field required"
    assert isinstance(payload["key_results"], dict), "key_results must be a dict"

    assert "error" in payload, "error field required"
    if payload["status"] == "success":
        assert payload["error"] is None, "error must be null on success"
    elif payload["status"] == "error":
        assert isinstance(payload["error"], str), "error must be a string on error status"

    assert "data" in payload, "data field required"
    assert isinstance(payload["data"], dict), "data must be a dict"

    assert "command" in payload, "command field required"
    cmd = payload["command"]
    assert isinstance(cmd, str), "command must be a string"
    assert cmd.startswith("phyloai "), f"command must start with 'phyloai ', got: {cmd[:50]!r}"
    assert len(cmd.split()) >= 3, f"command must have >= 3 tokens, got: {cmd!r}"

    _validate_command_flags(cmd)
    _validate_data_section(payload["data"])


def validate_params_completeness(
    payload: dict[str, Any],
    expected_keys: set[str],
) -> set[str]:
    """Assert params dict contains all expected signature keys.

    Returns the set of absent keys (empty = all present).
    """
    params = payload.get("params", {})
    absent = expected_keys - set(params.keys())
    assert not absent, f"params missing keys: {sorted(absent)}"
    return absent


def _validate_command_flags(cmd: str) -> None:
    """Check command string for common mistakes.

    - No `...` placeholders
    - Flags match Click decorator pattern (--lower-case-with-dashes)
    - Known bad patterns (e.g. --to-format when flag is --to)
    """
    assert "..." not in cmd, f"command contains '...' placeholder: {cmd!r}"

    try:
        tokens = shlex.split(cmd)
    except ValueError:
        return

    _KNOWN_BAD: dict[str, str] = {
        "--to-format": "--to",
        "--codon-view": "(nonexistent flag)",
    }
    for token in tokens:
        if token.startswith("--") and token in _KNOWN_BAD:
            raise AssertionError(
                f"command uses incorrect flag {token!r}; correct flag is {_KNOWN_BAD[token]!r}. "
                f"command: {cmd!r}"
            )

    for token in tokens:
        if " " in token and token.startswith("-"):
            raise AssertionError(
                f"command flag is joined with its value as a single token {token!r}. "
                f"flag and value must be separate tokens. command: {cmd!r}"
            )


def _validate_data_section(data: dict[str, Any]) -> None:
    if "cmd" in data:
        assert isinstance(data["cmd"], list), \
            f"data.cmd must be a list, got {type(data['cmd'])}"

    if "tool_stderr" in data:
        stderr = data["tool_stderr"]
        assert isinstance(stderr, str), \
            f"data.tool_stderr must be a str, got {type(stderr)}"
        for field in ("wall_time", "exit_code", "command"):
            assert field not in stderr, \
                f"data.tool_stderr must not contain '{field}' field"

    if "tool_log" in data:
        assert isinstance(data["tool_log"], str), \
            f"data.tool_log must be a str, got {type(data['tool_log'])}"

    for f in data.get("files", []):
        assert isinstance(f, dict), "files[] entries must be dicts"
        if "cmd" in f:
            assert isinstance(f["cmd"], list), \
                f"files[].cmd must be a list, got {type(f['cmd'])}"
            assert f.get("wall_time", 0) > 0 or f.get("status") == "dry_run", \
                f"files[].wall_time must be > 0 for non-dry-run tasks, got {f.get('wall_time')}"
        if "log_file" in f:
            log = f["log_file"]
            assert isinstance(log, str), \
                f"files[].log_file must be a string, got {type(log)}"
            assert log.startswith("logs/"), \
                f"files[].log_file must start with 'logs/', got: {log!r}"
