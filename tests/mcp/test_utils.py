from __future__ import annotations

import json
import os
import tempfile
from pathlib import Path

from phyloai.mcp.tools.utils import check_status, get_command_schema, read_report, read_result


def test_check_status_states_and_checkpoint() -> None:
    with tempfile.TemporaryDirectory() as tmp:
        output_dir = Path(tmp)
        assert check_status(str(output_dir))["status"] == "not_started"

        (output_dir / "job.json").write_text(json.dumps({"pid": 99999999, "command": "x"}))
        (output_dir / "checkpoint.json").write_text(json.dumps({"completed": 1, "total": 2}))
        unknown = check_status(str(output_dir))
        assert unknown["status"] == "unknown"
        assert unknown["checkpoint"] == {"completed": 1, "total": 2}

        (output_dir / "result.json").write_text(json.dumps({"status": "success", "key_results": {"n": 1}}))
        assert check_status(str(output_dir))["status"] == "success"


def test_check_status_stop_sidecar_takes_precedence() -> None:
    with tempfile.TemporaryDirectory() as tmp:
        output_dir = Path(tmp)
        (output_dir / "job.json").write_text(json.dumps({"pid": 99999999, "command": "x"}))
        (output_dir / "result.json").write_text(json.dumps({"status": "success", "old": True}))
        (output_dir / "stop_state.json").write_text(json.dumps({
            "state": "stop_requested",
            "pid": 123,
            "stop_signal": 15,
        }))

        payload = check_status(str(output_dir))
        assert payload["status"] == "unknown"
        assert payload["stop_state"]["state"] == "stop_requested"
        assert payload["result"]["old"] is True

        (output_dir / "stop_state.json").unlink()
        assert check_status(str(output_dir))["status"] == "success"


def test_check_status_corrupt_stop_sidecar_does_not_fall_back_to_success() -> None:
    with tempfile.TemporaryDirectory() as tmp:
        output_dir = Path(tmp)
        (output_dir / "result.json").write_text(json.dumps({"status": "success", "old": True}))
        (output_dir / "stop_state.json").write_text("{broken")

        payload = check_status(str(output_dir))
        assert payload["status"] == "unknown"
        assert "stop_state_error" in payload
        assert payload["result"]["old"] is True


def test_check_status_running_includes_stop_sidecar() -> None:
    with tempfile.TemporaryDirectory() as tmp:
        output_dir = Path(tmp)
        (output_dir / "job.json").write_text(json.dumps({"pid": os.getpid(), "command": "x"}))
        (output_dir / "stop_state.json").write_text(json.dumps({"state": "forced", "stop_signal": 15, "force_signal": 2}))

        payload = check_status(str(output_dir))
        assert payload["status"] == "running"
        assert payload["stop_state"]["state"] == "forced"


def test_read_result_and_report() -> None:
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        assert read_result(str(root))["status"] == "error"
        (root / "result.json").write_text(json.dumps({"status": "success"}))
        assert read_result(str(root))["status"] == "success"

        assert read_report(str(root))["status"] == "error"
        (root / "report").mkdir()
        (root / "report" / "report.json").write_text(json.dumps({"steps": []}))
        assert read_report(str(root))["steps"] == []


def test_get_command_schema() -> None:
    schema = get_command_schema("doctor")
    assert schema["name"] == "doctor"
    assert "output_format" in schema["inputSchema"]["properties"]
    assert "error" in get_command_schema("missing")
