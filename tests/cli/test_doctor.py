import json
import pytest
from click.testing import CliRunner
from unittest.mock import patch
from phyloai.cli.main import cli
from phyloai.core.env import ToolStatus, ToolInfo
from pathlib import Path


def _mock_tools():
    return {
        "iqtree2": ToolInfo("iqtree2", ToolStatus.OK,
                            Path("/usr/bin/iqtree2"), "2.3.1"),
        "mafft":   ToolInfo("mafft",   ToolStatus.OK,
                            Path("/usr/bin/mafft"), "7.520"),
        "pb_mpi":  ToolInfo("pb_mpi",  ToolStatus.MISSING,
                            note="install: http://www.phylobayes.org"),
    }


def test_doctor_exits_zero():
    runner = CliRunner()
    with patch("phyloai.cli.doctor.ToolEnv") as MockEnv:
        MockEnv.return_value.check_all.return_value = _mock_tools()
        result = runner.invoke(cli, ["doctor"])
    assert result.exit_code == 0


def test_doctor_shows_ok_tools():
    runner = CliRunner()
    with patch("phyloai.cli.doctor.ToolEnv") as MockEnv:
        MockEnv.return_value.check_all.return_value = _mock_tools()
        result = runner.invoke(cli, ["doctor"])
    assert "iqtree2" in result.output
    assert "mafft" in result.output


def test_doctor_shows_missing_tools():
    runner = CliRunner()
    with patch("phyloai.cli.doctor.ToolEnv") as MockEnv:
        MockEnv.return_value.check_all.return_value = _mock_tools()
        result = runner.invoke(cli, ["doctor"])
    assert "pb_mpi" in result.output


def test_doctor_json_output():
    runner = CliRunner()
    with patch("phyloai.cli.doctor.ToolEnv") as MockEnv:
        MockEnv.return_value.check_all.return_value = _mock_tools()
        result = runner.invoke(cli, ["doctor", "--output-format", "json"])
    assert result.exit_code == 0
    data = json.loads(result.output)
    assert "iqtree2" in data
    assert data["iqtree2"]["status"] == "ok"
    assert data["pb_mpi"]["status"] == "missing"
