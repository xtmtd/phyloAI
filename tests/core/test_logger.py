import pytest
from pathlib import Path
from phyloai.core.logger import StepLogger
from phyloai.core.schema import ToolResult


def test_logger_creates_log_file(tmp_path):
    logger = StepLogger(run_dir=tmp_path)
    logger.write("align", ToolResult(
        tool="mafft",
        command="mafft --auto input.fa > output.fa",
        returncode=0,
        stdout="Done.",
        stderr="",
        wall_time=3.2,
    ))
    log_file = tmp_path / "logs" / "align.log"
    assert log_file.exists()


def test_logger_log_contains_required_fields(tmp_path):
    logger = StepLogger(run_dir=tmp_path)
    result = ToolResult(
        tool="iqtree2",
        command="iqtree2 -s matrix.fa",
        returncode=0,
        stdout="Analysis done",
        stderr="",
        wall_time=45.1,
    )
    logger.write("iqtree", result)
    content = (tmp_path / "logs" / "iqtree.log").read_text()
    assert "iqtree2" in content
    assert "iqtree2 -s matrix.fa" in content
    assert "45.1" in content
    assert "returncode: 0" in content
    assert "Analysis done" in content


def test_logger_appends_on_retry(tmp_path):
    logger = StepLogger(run_dir=tmp_path)
    result = ToolResult("mafft", "mafft in.fa", 0, "ok", "", 1.0)
    logger.write("align", result)
    logger.write("align", result)
    content = (tmp_path / "logs" / "align.log").read_text()
    assert content.count("mafft in.fa") == 2


def test_logger_logs_dir_is_created(tmp_path):
    run_dir = tmp_path / "runs" / "run001"
    logger = StepLogger(run_dir=run_dir)
    result = ToolResult("echo", "echo hi", 0, "hi", "", 0.01)
    logger.write("test_step", result)
    assert (run_dir / "logs").is_dir()
