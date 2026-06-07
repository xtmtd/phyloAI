import json
import pytest
from click.testing import CliRunner
from unittest.mock import patch
from phyloai.cli.main import cli
from phyloai.core.env import ToolStatus, ToolInfo
from pathlib import Path


def _mock_tools():
    return {
        "iqtree3": ToolInfo("iqtree3", ToolStatus.OK,
                            Path("/usr/bin/iqtree3"), "3.0.1"),
        "astral-hybrid": ToolInfo("astral-hybrid", ToolStatus.OK,
                            Path("/usr/bin/astral-hybrid"), "5.7.8"),
        "mafft":   ToolInfo("mafft",   ToolStatus.OK,
                            Path("/usr/bin/mafft"), "7.520"),
        "phykit":  ToolInfo("phykit",  ToolStatus.OK,
                            Path("/usr/bin/phykit"), "2.1.2"),
        "java":    ToolInfo("java",    ToolStatus.OK,
                            Path("/usr/bin/java"), "21.0.2"),
        "julia":   ToolInfo("julia",   ToolStatus.MISSING,
                            note="install: https://julialang.org/downloads/"),
        "correction_multi.jl": ToolInfo("correction_multi.jl", ToolStatus.MISSING,
                            note="install: https://github.com/chaoszhang/TAPER"),
        "pb_mpi":  ToolInfo("pb_mpi",  ToolStatus.MISSING,
                            note="install: https://github.com/bayesiancook/pbmpi"),
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
    assert "iqtree3" in result.output
    assert "astral-hybrid" in result.output
    assert "phykit" in result.output
    assert "mafft" in result.output


def test_doctor_shows_missing_tools():
    runner = CliRunner()
    with patch("phyloai.cli.doctor.ToolEnv") as MockEnv:
        MockEnv.return_value.check_all.return_value = _mock_tools()
        result = runner.invoke(cli, ["doctor"])
    assert "pb_mpi" in result.output
    assert "julia" in result.output


def test_doctor_json_output():
    runner = CliRunner()
    with patch("phyloai.cli.doctor.ToolEnv") as MockEnv:
        MockEnv.return_value.check_all.return_value = _mock_tools()
        result = runner.invoke(cli, ["doctor", "--output-format", "json"])
    assert result.exit_code == 0
    data = json.loads(result.output)
    assert "iqtree3" in data
    assert data["iqtree3"]["status"] == "ok"
    assert data["pb_mpi"]["status"] == "missing"


def test_doctor_json_version_is_plain_number():
    runner = CliRunner()
    with patch("phyloai.cli.doctor.ToolEnv") as MockEnv:
        MockEnv.return_value.check_all.return_value = _mock_tools()
        result = runner.invoke(cli, ["doctor", "--output-format", "json"])

    data = json.loads(result.output)
    assert data["iqtree3"]["version"] == "3.0.1"
    assert data["astral-hybrid"]["version"] == "5.7.8"
    assert data["phykit"]["version"] == "2.1.2"


def test_doctor_help_mentions_default_text_output():
    runner = CliRunner()
    result = runner.invoke(cli, ["doctor", "-h"])

    assert result.exit_code == 0
    assert "Default: text." in result.output


def test_doctor_json_includes_runtime_checks():
    runner = CliRunner()
    with patch("phyloai.cli.doctor.ToolEnv") as MockEnv:
        MockEnv.return_value.check_all.return_value = _mock_tools()
        result = runner.invoke(cli, ["doctor", "--output-format", "json"])

    data = json.loads(result.output)
    assert data["java"]["status"] == "ok"
    assert data["julia"]["status"] == "missing"
    assert data["correction_multi.jl"]["status"] == "missing"
    assert data["correction_multi.jl"]["note"] == "install: https://github.com/chaoszhang/TAPER"
