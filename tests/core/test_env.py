import pytest
from unittest.mock import patch
from pathlib import Path
from phyloai.core.env import ToolEnv, ToolStatus, ToolInfo


def test_tool_status_values():
    assert ToolStatus.OK == "ok"
    assert ToolStatus.WARN == "warn"
    assert ToolStatus.MISSING == "missing"


def test_tool_info_defaults():
    info = ToolInfo(name="iqtree2")
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
    for key in ["iqtree2", "mafft", "trimal"]:
        assert key in results


def test_get_tool_path_raises_when_missing():
    env = ToolEnv()
    with patch("shutil.which", return_value=None):
        env._tools["fake_tool"] = ToolInfo(name="fake_tool")
        with pytest.raises(FileNotFoundError, match="fake_tool"):
            env.require("fake_tool")
