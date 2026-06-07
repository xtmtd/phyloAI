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


def test_msa_collection_auto_count(tmp_path):
    # create 3 .fa files in a temp dir
    for i in range(3):
        (tmp_path / f"gene_{i}.fa").write_text(f">Taxon\nACGT\n")
    c = MSACollection(directory=tmp_path)
    assert c.count == 3


def test_tree_set_auto_count(tmp_path):
    for i in range(2):
        (tmp_path / f"gene_{i}.treefile").write_text("(A,B);")
    ts = TreeSet(directory=tmp_path)
    assert ts.count == 2


def test_run_record_add_step():
    from pathlib import Path
    record = RunRecord(run_dir=Path("./runs/run001"))
    result = ToolResult(
        tool="mafft", command="mafft in.fa", returncode=0,
        stdout="done", stderr="", wall_time=1.5,
    )
    record.add_step("align", {"method": "linsi"}, result)
    assert len(record.steps) == 1
    assert record.steps[0]["step"] == "align"
    assert record.steps[0]["result"]["success"] is True


def test_run_record_explicit_version():
    from pathlib import Path
    record = RunRecord(run_dir=Path("./runs/run001"), phyloai_version="9.9.9")
    assert record.phyloai_version == "9.9.9"
