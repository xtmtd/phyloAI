import json
import pytest
from click.testing import CliRunner
from unittest.mock import patch
from phyloai.cli.main import cli
from phyloai.core.env import ToolStatus, ToolInfo, ToolEnv
from pathlib import Path


def _mock_tools():
    return {
        "iqtree3": ToolInfo("iqtree3", ToolStatus.OK,
                            Path("/usr/bin/iqtree3"), "3.0.1"),
        "wastral": ToolInfo("wastral", ToolStatus.OK,
                            Path("/usr/bin/wastral"), "5.7.8"),
        "mafft":   ToolInfo("mafft",   ToolStatus.OK,
                            Path("/usr/bin/mafft"), "7.520"),
        "bmge":    ToolInfo("bmge",    ToolStatus.OK,
                            Path("/usr/local/bin/BMGE.jar"), "1.12"),
        "java":    ToolInfo("java",    ToolStatus.OK,
                            Path("/usr/bin/java"), "21.0.2"),
        "julia":   ToolInfo("julia",   ToolStatus.MISSING,
                            note="install: https://julialang.org/downloads/"),
        "correction_multi.jl": ToolInfo("correction_multi.jl", ToolStatus.OK,
                            Path("/usr/local/bin/correction_multi.jl"), "1.0.0"),
        "pb_mpi":  ToolInfo("pb_mpi",  ToolStatus.MISSING,
                            note="install: https://github.com/bayesiancook/pbmpi"),
        "bpcomp":  ToolInfo("bpcomp",  ToolStatus.MISSING,
                            note="install: https://github.com/bayesiancook/pbmpi"),
        "tracecomp":  ToolInfo("tracecomp",  ToolStatus.MISSING,
                            note="install: https://github.com/bayesiancook/pbmpi"),
        "readpb_mpi":  ToolInfo("readpb_mpi",  ToolStatus.MISSING,
                            note="install: https://github.com/bayesiancook/pbmpi"),
        "mpirun":  ToolInfo("mpirun",  ToolStatus.MISSING,
                            note="install: https://www.open-mpi.org"),
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
    assert "wastral" in result.output
    assert "mafft" in result.output


def test_doctor_shows_bmge_jar_display_name():
    runner = CliRunner()
    with patch("phyloai.cli.doctor.ToolEnv") as MockEnv:
        MockEnv.return_value.check_all.return_value = _mock_tools()
        result = runner.invoke(cli, ["doctor"])

    assert "BMGE.jar" in result.output


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
    assert data["wastral"]["version"] == "5.7.8"


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
    assert data["correction_multi.jl"]["status"] == "ok"
    assert data["correction_multi.jl"]["note"] == ""
    assert data["bmge"]["note"] == ""


def test_bmge_jar_resolves_from_path_alias(tmp_path):
    bundled_root = tmp_path / "bundled" / "BMGE-1.12"
    bundled_root.mkdir(parents=True)
    (bundled_root / "BMGE.jar").write_text("jar placeholder")

    env = ToolEnv()
    env._bundled_dir = tmp_path / "bundled"

    with patch("shutil.which", side_effect=lambda name: "/usr/local/bin/BMGE.jar" if name == "BMGE.jar" else None):
        info = env.check_all()["bmge"]

    assert info.status == ToolStatus.OK
    assert info.path == Path("/usr/local/bin/BMGE.jar")


def test_explicit_bmge_path_override_is_preferred_over_path(tmp_path):
    bmge_jar = tmp_path / "BMGE.jar"
    bmge_jar.write_text("jar placeholder")

    env = ToolEnv(tool_paths={"bmge": bmge_jar})

    with patch.object(env, "_get_version", return_value="1.12"):
        with patch("shutil.which", return_value="/usr/local/bin/BMGE.jar"):
            info = env.check_all()["bmge"]

    assert info.status == ToolStatus.OK
    assert info.path == bmge_jar
    assert info.version == "1.12"


def test_doctor_json_includes_phylobayes_mpi_tools():
    runner = CliRunner()
    with patch("phyloai.cli.doctor.ToolEnv") as MockEnv:
        MockEnv.return_value.check_all.return_value = _mock_tools()
        result = runner.invoke(cli, ["doctor", "--output-format", "json"])

    data = json.loads(result.output)
    for name in ["pb_mpi", "bpcomp", "tracecomp", "readpb_mpi", "mpirun"]:
        assert name in data
        assert data[name]["status"] in ("ok", "missing", "warn")
