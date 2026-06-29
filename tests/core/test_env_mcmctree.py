"""Test mcmctree version extraction from stdout."""
from __future__ import annotations
from unittest.mock import MagicMock, patch
from pathlib import Path

from phyloai.core.env import ToolEnv, TOOL_REGISTRY


def _make_result(stdout: str):
    r = MagicMock()
    r.stdout = stdout
    r.stderr = ""
    return r


def test_mcmctree_version_pattern_extracts_paml_version():
    env = ToolEnv()
    fake_output = "MCMCTREE in paml version 4.10.10, 27 Jan 2026\n"
    with patch("subprocess.run", return_value=_make_result(fake_output)):
        ver = env._get_version(
            Path("/usr/bin/mcmctree"), [],
            version_pattern=r"paml version (\d+(?:\.\d+)+)",
        )
    assert ver == "4.10.10"


def test_mcmctree_registry_uses_version_args_list():
    meta = TOOL_REGISTRY["mcmctree"]
    assert "version_args" in meta
    assert meta["version_args"] == []
    assert "version_flag" not in meta


def test_mcmctree_version_falls_back_to_install_dir(tmp_path):
    tool_dir = tmp_path / "paml-4.10.10" / "bin"
    tool_dir.mkdir(parents=True)
    tool = tool_dir / "mcmctree"
    tool.write_text("#!/bin/sh\n")

    env = ToolEnv()

    with patch("subprocess.run", return_value=_make_result("")):
        ver = env._get_version(
            tool,
            [],
            version_pattern=r"paml version (\d+(?:\.\d+)+)",
        )

    assert ver == "4.10.10"
