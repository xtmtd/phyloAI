# tests/core/test_schema.py
import pytest
from pathlib import Path
from phyloai.core.schema import MSACollection, TreeSet, RunRecord, ToolResult


def test_msa_collection_requires_directory():
    with pytest.raises(TypeError):
        MSACollection()  # missing required field


def test_msa_collection_defaults():
    c = MSACollection(directory=Path("./alignments"))
    assert c.seq_type == "AA"
    assert c.file_extension == ".fa"
    assert c.count == 0


def test_tool_result_success():
    r = ToolResult(
        tool="iqtree2",
        command="iqtree2 -s matrix.fa",
        returncode=0,
        stdout="Analysis done",
        stderr="",
        wall_time=12.5,
    )
    assert r.success is True


def test_tool_result_failure():
    r = ToolResult(
        tool="iqtree2",
        command="iqtree2 -s missing.fa",
        returncode=1,
        stdout="",
        stderr="ERROR: file not found",
        wall_time=0.1,
    )
    assert r.success is False


def test_run_record_to_dict():
    record = RunRecord(run_dir=Path("./runs/run001"))
    d = record.to_dict()
    assert "run_dir" in d
    assert "steps" in d
    assert isinstance(d["steps"], list)


def test_tree_set_defaults():
    ts = TreeSet(directory=Path("./trees"))
    assert ts.format == "newick"
    assert ts.count == 0
