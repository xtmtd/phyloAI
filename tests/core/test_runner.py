import pytest
import time
from pathlib import Path
from phyloai.core.runner import Runner
from phyloai.core.schema import ToolResult


def test_runner_run_success():
    runner = Runner()
    result = runner.run(["echo", "hello"], tool_name="echo")
    assert result.success is True
    assert "hello" in result.stdout
    assert result.tool == "echo"
    assert result.wall_time > 0


def test_runner_run_failure():
    runner = Runner()
    result = runner.run(["false"], tool_name="false")
    assert result.success is False
    assert result.returncode != 0


def test_runner_run_returns_tool_result():
    runner = Runner()
    result = runner.run(["echo", "test"], tool_name="echo")
    assert isinstance(result, ToolResult)


def test_runner_raises_on_missing_executable():
    runner = Runner()
    with pytest.raises(FileNotFoundError):
        runner.run(["nonexistent_binary_xyz_abc"], tool_name="fake")


def test_runner_timeout():
    runner = Runner(timeout=1)
    with pytest.raises(TimeoutError):
        runner.run(["sleep", "10"], tool_name="sleep")


def test_runner_captures_stderr():
    runner = Runner()
    result = runner.run(
        ["ls", "/nonexistent_path_xyz_abc_123"],
        tool_name="ls"
    )
    assert result.success is False
    assert result.stderr != "" or result.returncode != 0
