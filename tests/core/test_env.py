import pytest
from unittest.mock import patch
from pathlib import Path
from phyloai.core.env import ToolEnv, ToolStatus, ToolInfo


def test_tool_status_values():
    assert ToolStatus.OK == "ok"
    assert ToolStatus.WARN == "warn"
    assert ToolStatus.MISSING == "missing"


def test_tool_info_defaults():
    info = ToolInfo(name="iqtree3")
    assert info.status == ToolStatus.MISSING
    assert info.path is None
    assert info.version is None


def test_detect_present_tool():
    env = ToolEnv()
    with patch("shutil.which", return_value="/usr/bin/echo"):
        with patch.object(env, "_get_version", return_value="1.0"):
            result = env._detect_tool("echo", version_flag="--version")
    assert result.status == ToolStatus.OK
    assert result.path == Path("/usr/bin/echo")


def test_detect_missing_tool():
    env = ToolEnv()
    with patch("shutil.which", return_value=None):
        result = env._detect_tool("nonexistent_tool_xyz")
    assert result.status == ToolStatus.MISSING
    assert result.path is None


def test_check_all_returns_dict():
    env = ToolEnv()
    results = env.check_all()
    assert isinstance(results, dict)
    for key in ["iqtree3", "mafft", "trimal"]:
        assert key in results


def test_check_all_omits_removed_simphy_entry():
    env = ToolEnv()
    results = env.check_all()
    assert "simphy" not in results


def test_explicit_tool_path_override_is_preferred(tmp_path):
    custom_tool = tmp_path / "iqtree3"
    custom_tool.write_text("#!/bin/sh\n")

    env = ToolEnv(tool_paths={"iqtree3": custom_tool})

    with patch("shutil.which", return_value="/usr/bin/iqtree3"):
        with patch.object(env, "_get_version", return_value="3.0.1"):
            result = env._detect_tool("iqtree3", version_flag="--version")

    assert result.status == ToolStatus.OK
    assert result.path == custom_tool


def test_explicit_missing_tool_path_reports_missing(tmp_path):
    env = ToolEnv(tool_paths={"iqtree3": tmp_path / "missing-iqtree3"})

    with patch("shutil.which", return_value="/usr/bin/iqtree3"):
        result = env._detect_tool("iqtree3", version_flag="--version")

    assert result.status == ToolStatus.MISSING
    assert result.path is None


def test_get_tool_path_raises_when_missing():
    env = ToolEnv()
    with patch("shutil.which", return_value=None):
        env._tools["fake_tool"] = ToolInfo(name="fake_tool")
        with pytest.raises(FileNotFoundError, match="fake_tool"):
            env.require("fake_tool")


def test_check_all_includes_runtime_and_taper_entries():
    env = ToolEnv()
    results = env.check_all()

    for key in ["correction_multi.jl", "java", "julia"]:
        assert key in results


def test_get_version_uses_alternative_args_when_needed(tmp_path):
    tool = tmp_path / "astral-hybrid"
    tool.write_text("#!/bin/sh\n")

    env = ToolEnv()

    with patch("subprocess.run") as mock_run:
        mock_run.side_effect = [
            type("Result", (), {"stdout": "", "stderr": "", "returncode": 0})(),
            type("Result", (), {"stdout": "ASTRAL-HYBRID version 5.7.8\n", "stderr": "", "returncode": 0})(),
        ]
        version = env._get_version(tool, [["--version"], ["-h"]])

    assert version == "5.7.8"


def test_get_version_runs_bmge_jar_via_java(tmp_path):
    tool = tmp_path / "BMGE.jar"
    tool.write_text("jar placeholder")

    env = ToolEnv()

    with patch("subprocess.run") as mock_run:
        mock_run.return_value = type(
            "Result", (),
            {"stdout": "", "stderr": "BMGE (version 1.12) arguments :\n", "returncode": 0} 
        )()
        version = env._get_version(tool, [["-?"]])

    assert version == "1.12"
    mock_run.assert_called_once_with(
        ["java", "-jar", str(tool), "-?"],
        capture_output=True,
        text=True,
        timeout=5,
    )
