"""Tests for phyloai run command."""

from __future__ import annotations

import json
from pathlib import Path
from unittest.mock import patch, MagicMock

import click
import pytest
from click.testing import CliRunner

from phyloai.cli.main import cli
from phyloai.cli.commands._run_pipeline import (
    _build_run_params,
    _build_run_checkpoint,
    _load_run_checkpoint,
    _validate_run_resume,
)
from phyloai.core.checkpoint import canonical_params_hash


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_seq_dir(tmp_path: Path) -> Path:
    d = tmp_path / "markers"
    d.mkdir()
    (d / "gene1.fa").write_text(">sp1\nMKT\n>sp2\nMKA\n")
    (d / "gene2.fa").write_text(">sp1\nGHT\n>sp2\nGHA\n")
    return d


def _mock_step_result(n_files: int = 2) -> dict:
    return {
        "status": "success",
        "command": "phyloai pretree ...",
        "wall_time": 1.0,
        "tool_versions": {},
        "params": {},
        "key_results": {},
        "error": None,
        "data": {"files": [{"input": f"g{i}.fa"} for i in range(n_files)]},
    }


# ---------------------------------------------------------------------------
# Task 1: CLI help
# ---------------------------------------------------------------------------

def test_run_help() -> None:
    runner = CliRunner()
    result = runner.invoke(cli, ["run", "--help"])
    assert result.exit_code == 0
    assert "--seq-dir" in result.output
    assert "--mode" in result.output
    assert "--speed" in result.output
    assert "--resume" in result.output
    assert "supermatrix" in result.output
    assert "supertree" in result.output
    assert "normal" in result.output
    assert "fast" in result.output


# ---------------------------------------------------------------------------
# Task 2: Checkpoint helpers
# ---------------------------------------------------------------------------

def test_build_run_params_keys() -> None:
    params = _build_run_params(
        seq_dir=Path("./markers"),
        mode="supermatrix",
        speed="normal",
        threads=4,
        output_dir=Path("runs/run"),
    )
    assert "mode" in params
    assert "speed" in params
    assert "threads" in params
    assert params["mode"] == "supermatrix"
    assert params["speed"] == "normal"


def test_build_run_checkpoint_schema(tmp_path: Path) -> None:
    params = _build_run_params(Path("m"), "supermatrix", "normal", 4, Path("r"))
    ckpt = _build_run_checkpoint("phyloai run ...", params, mode="supermatrix", speed="normal")
    assert ckpt["schema_version"] == 1
    assert ckpt["step"] == "run"
    assert ckpt["status"] == "running"
    assert "steps" in ckpt
    assert isinstance(ckpt["steps"], list)


def test_validate_resume_mismatch_raises(tmp_path: Path) -> None:
    params = _build_run_params(Path("m"), "supermatrix", "normal", 4, Path("r"))
    ckpt = _build_run_checkpoint("cmd", params, mode="supermatrix", speed="normal")
    ckpt["params_hash"] = "sha256:deadbeef"
    with pytest.raises(click.ClickException, match="Parameter mismatch"):
        _validate_run_resume(ckpt, canonical_params_hash(params))


def test_load_run_checkpoint_missing_raises(tmp_path: Path) -> None:
    with pytest.raises(click.ClickException, match="run_checkpoint.json"):
        _load_run_checkpoint(tmp_path / "run_checkpoint.json")


# ---------------------------------------------------------------------------
# Task 3: Output directory setup and guards
# ---------------------------------------------------------------------------

def test_run_resume_and_overwrite_mutually_exclusive(tmp_path: Path) -> None:
    runner = CliRunner()
    seq_dir = _make_seq_dir(tmp_path)
    result = runner.invoke(cli, [
        "run", "--seq-dir", str(seq_dir),
        "--resume", "--overwrite",
    ])
    assert result.exit_code == 1


def test_run_resume_without_checkpoint_exits_1(tmp_path: Path) -> None:
    runner = CliRunner()
    seq_dir = _make_seq_dir(tmp_path)
    out_dir = tmp_path / "run"
    out_dir.mkdir()
    result = runner.invoke(cli, [
        "run", "--seq-dir", str(seq_dir),
        "--output-dir", str(out_dir),
        "--resume",
    ])
    assert result.exit_code == 1


def test_run_nonempty_output_dir_exits_1(tmp_path: Path) -> None:
    runner = CliRunner()
    seq_dir = _make_seq_dir(tmp_path)
    out_dir = tmp_path / "run"
    out_dir.mkdir()
    (out_dir / "somefile").write_text("x")
    result = runner.invoke(cli, [
        "run", "--seq-dir", str(seq_dir),
        "--output-dir", str(out_dir),
    ])
    assert result.exit_code == 1


# ---------------------------------------------------------------------------
# Task 4: Dry-run
# ---------------------------------------------------------------------------

def test_run_dry_run_supermatrix_normal(tmp_path: Path) -> None:
    runner = CliRunner()
    seq_dir = _make_seq_dir(tmp_path)
    result = runner.invoke(cli, [
        "run", "--seq-dir", str(seq_dir),
        "--mode", "supermatrix", "--speed", "normal",
        "--dry-run",
    ])
    assert result.exit_code == 0
    assert "convert" in result.output.lower()
    assert "align" in result.output.lower()
    assert "trim" in result.output.lower()
    assert "filter" in result.output.lower()
    assert "concat" in result.output.lower()
    assert "iq-tree" in result.output.lower()


def test_run_dry_run_supermatrix_fast_no_filter(tmp_path: Path) -> None:
    runner = CliRunner()
    seq_dir = _make_seq_dir(tmp_path)
    result = runner.invoke(cli, [
        "run", "--seq-dir", str(seq_dir),
        "--mode", "supermatrix", "--speed", "fast",
        "--dry-run",
    ])
    assert result.exit_code == 0
    assert "fasttree" in result.output.lower()


def test_run_dry_run_supertree_normal(tmp_path: Path) -> None:
    runner = CliRunner()
    seq_dir = _make_seq_dir(tmp_path)
    result = runner.invoke(cli, [
        "run", "--seq-dir", str(seq_dir),
        "--mode", "supertree", "--speed", "normal",
        "--dry-run",
    ])
    assert result.exit_code == 0
    output_lower = result.output.lower()
    assert "gene" in output_lower or "genetree" in output_lower
    assert "iq-tree" in output_lower  # normal mode uses IQ-TREE3 for gene trees
    assert "wastral" in output_lower


# ---------------------------------------------------------------------------
# Task 5-6: Full pipeline mocked
# ---------------------------------------------------------------------------

def test_run_calls_convert_and_align_and_trim(tmp_path: Path) -> None:
    runner = CliRunner()
    seq_dir = _make_seq_dir(tmp_path)
    out_dir = tmp_path / "run"

    r = _mock_step_result()

    with patch("phyloai.pretree.convert.convert_input", return_value=r) as mock_conv, \
         patch("phyloai.pretree.align.run_align", return_value=r) as mock_align, \
         patch("phyloai.pretree.trim.run_trim", return_value=r) as mock_trim, \
         patch("phyloai.pretree.filter.run_taper", side_effect=NotImplementedError) as _mock_taper, \
         patch("phyloai.pretree.concat.run_concat", side_effect=NotImplementedError) as _mock_concat, \
         patch("phyloai.tree.ml_iqtree.run_iqtree", side_effect=NotImplementedError) as _mock_iq:
        result = runner.invoke(cli, [
            "run", "--seq-dir", str(seq_dir),
            "--mode", "supermatrix", "--speed", "normal",
            "--output-dir", str(out_dir),
        ])

    assert mock_conv.called
    assert mock_align.called
    assert mock_trim.called
    assert (out_dir / "1-convert").exists()
    assert (out_dir / "2-align").exists()
    assert (out_dir / "3-trim").exists()


def test_run_supermatrix_normal_full_pipeline_mocked(tmp_path: Path) -> None:
    runner = CliRunner()
    seq_dir = _make_seq_dir(tmp_path)
    out_dir = tmp_path / "run"

    r = _mock_step_result()

    def _mock_convert(*args, **kwargs) -> dict:
        out = kwargs.get("output_dir", args[1] if len(args) > 1 else None)
        if out:
            (Path(out) / "seqs").mkdir(parents=True, exist_ok=True)
            (Path(out) / "seqs" / "gene1.fa").write_text(">sp1\nA")
        return r

    def _mock_iqtree(*args, **kwargs) -> dict:
        out = kwargs.get("output_dir", Path("."))
        tf = Path(out) / "iqtree.treefile"
        tf.parent.mkdir(parents=True, exist_ok=True)
        tf.write_text("(sp1,sp2);")
        return {**r, "data": {"output": str(tf)}}

    with patch("phyloai.pretree.convert.convert_input", side_effect=_mock_convert), \
         patch("phyloai.pretree.align.run_align", return_value=r), \
         patch("phyloai.pretree.trim.run_trim", return_value=r), \
         patch("phyloai.pretree.filter.run_taper", return_value=r) as mock_taper, \
         patch("phyloai.pretree.concat.run_concat",
               return_value={**r, "data": {"matrix_file": str(out_dir / "5-concat" / "matrix.fa")},
                             "key_results": {"total_length": 10, "n_taxa": 2}}) as mock_concat, \
         patch("phyloai.tree.ml_iqtree.run_iqtree", side_effect=_mock_iqtree) as mock_iqtree:
        result = runner.invoke(cli, [
            "run", "--seq-dir", str(seq_dir),
            "--mode", "supermatrix", "--speed", "normal",
            "--output-dir", str(out_dir),
        ])

    assert result.exit_code == 0, result.output
    assert mock_taper.called
    assert mock_concat.called
    assert mock_iqtree.called
    result_json = out_dir / "result.json"
    assert result_json.exists()
    data = json.loads(result_json.read_text())
    assert data["status"] == "success"
    assert "final_tree" in data["key_results"]


def test_run_supertree_fast_full_pipeline_mocked(tmp_path: Path) -> None:
    runner = CliRunner()
    seq_dir = _make_seq_dir(tmp_path)
    out_dir = tmp_path / "run"

    r = _mock_step_result()

    def _mock_convert(*args, **kwargs) -> dict:
        out = kwargs.get("output_dir", args[1] if len(args) > 1 else None)
        if out:
            (Path(out) / "seqs").mkdir(parents=True, exist_ok=True)
            (Path(out) / "seqs" / "gene1.fa").write_text(">sp1\nA")
        return r

    def _mock_wastral(*args, **kwargs) -> dict:
        out = kwargs.get("output_dir", Path("."))
        tf = Path(out) / "wastral.tre"
        tf.parent.mkdir(parents=True, exist_ok=True)
        tf.write_text("(sp1,sp2);")
        return {**r, "data": {"output_tree": str(tf)}}

    with patch("phyloai.pretree.convert.convert_input", side_effect=_mock_convert), \
         patch("phyloai.pretree.align.run_align", return_value=r), \
         patch("phyloai.pretree.trim.run_trim", return_value=r), \
         patch("phyloai.pretree.filter.run_taper",
               side_effect=AssertionError("should not be called")) as mock_taper, \
         patch("phyloai.tree.ml.run_fasttree",
               return_value={**r, "data": {"trees_dir": str(out_dir / "5-genetrees" / "trees")}}) as mock_ft, \
         patch("phyloai.tree.msc.run_wastral", side_effect=_mock_wastral) as mock_wastral:
        result = runner.invoke(cli, [
            "run", "--seq-dir", str(seq_dir),
            "--mode", "supertree", "--speed", "fast",
            "--output-dir", str(out_dir),
        ])

    assert result.exit_code == 0, result.output
    assert not mock_taper.called
    assert mock_ft.called
    assert mock_wastral.called
    result_json = out_dir / "result.json"
    assert result_json.exists()
    data = json.loads(result_json.read_text())
    assert data["status"] == "success"
    assert not (out_dir / "5-concat").exists()
    assert not (out_dir / "4-filter").exists()


# ---------------------------------------------------------------------------
# Task 7: Error handling
# ---------------------------------------------------------------------------

def test_run_step_failure_writes_error_result_json(tmp_path: Path) -> None:
    runner = CliRunner()
    seq_dir = _make_seq_dir(tmp_path)
    out_dir = tmp_path / "run"

    r = _mock_step_result()

    with patch("phyloai.pretree.convert.convert_input", return_value=r), \
         patch("phyloai.pretree.align.run_align",
               side_effect=RuntimeError("MAFFT not found")):
        result = runner.invoke(cli, [
            "run", "--seq-dir", str(seq_dir),
            "--output-dir", str(out_dir),
        ])

    assert result.exit_code == 2, f"expected exit 2, got {result.exit_code}"
    result_json_path = out_dir / "result.json"
    assert result_json_path.exists()
    data = json.loads(result_json_path.read_text())
    assert data["status"] == "error"


# ---------------------------------------------------------------------------
# Task 8: Resume
# ---------------------------------------------------------------------------

def test_run_resume_skips_completed_steps(tmp_path: Path) -> None:
    runner = CliRunner()
    seq_dir = _make_seq_dir(tmp_path)
    out_dir = tmp_path / "run"
    out_dir.mkdir()

    # Pre-create convert and align dirs with success result.json
    for subdir in ["1-convert", "2-align"]:
        d = out_dir / subdir
        d.mkdir()
        (d / "seqs").mkdir()
        (d / "result.json").write_text(json.dumps({"status": "success", "data": {}}))
    (out_dir / "1-convert" / "seqs" / "gene1.fa").write_text(">sp1\nMKT\n")
    (out_dir / "2-align" / "seqs" / "gene1.fa").write_text(">sp1\nMKT\n")

    # Write a run_checkpoint.json with convert+align success, trim pending
    from phyloai.cli.commands._run_pipeline import (
        _build_run_params, _build_run_checkpoint, _save_run_checkpoint,
    )
    params = _build_run_params(seq_dir, "supermatrix", "normal", 4, out_dir)
    ckpt = _build_run_checkpoint("phyloai run ...", params, mode="supermatrix", speed="normal")
    for s in ckpt["steps"]:
        s["output_dir"] = str(out_dir / {
            "convert": "1-convert", "align": "2-align", "trim": "3-trim",
            "filter_taper": "4-filter", "concat": "5-concat", "tree": "6-tree",
        }[s["name"]])
    ckpt["steps"][0]["status"] = "success"
    ckpt["steps"][1]["status"] = "success"
    _save_run_checkpoint(ckpt, out_dir / "run_checkpoint.json")

    r = _mock_step_result()
    mock_convert = MagicMock(return_value=r)
    mock_align = MagicMock(return_value=r)

    def _mock_iqtree(*args, **kwargs) -> dict:
        out = kwargs.get("output_dir", Path("."))
        tf = Path(out) / "iqtree.treefile"
        tf.parent.mkdir(parents=True, exist_ok=True)
        tf.write_text("(sp1,sp2);")
        return {**r, "data": {"output": str(tf)}}

    with patch("phyloai.pretree.convert.convert_input", mock_convert), \
         patch("phyloai.pretree.align.run_align", mock_align), \
         patch("phyloai.pretree.trim.run_trim", return_value=r), \
         patch("phyloai.pretree.filter.run_taper", return_value=r), \
         patch("phyloai.pretree.concat.run_concat", return_value=r), \
         patch("phyloai.tree.ml_iqtree.run_iqtree", side_effect=_mock_iqtree):
        result = runner.invoke(cli, [
            "run", "--seq-dir", str(seq_dir),
            "--output-dir", str(out_dir),
            "--resume",
        ])

    assert result.exit_code == 0, result.output
    assert not mock_convert.called
    assert not mock_align.called


# ---------------------------------------------------------------------------
# Additional tests from code review
# ---------------------------------------------------------------------------

def test_resume_treats_running_like_interrupted(tmp_path: Path) -> None:
    """--resume with step status 'running' should resume, not overwrite."""
    runner = CliRunner()
    seq_dir = _make_seq_dir(tmp_path)
    out_dir = tmp_path / "run"
    out_dir.mkdir()

    # Pre-create convert dir with success, align dir empty (running will get reset)
    for subdir in ["1-convert"]:
        d = out_dir / subdir
        d.mkdir()
        (d / "seqs").mkdir()
        (d / "result.json").write_text(json.dumps({"status": "success", "data": {}}))
    (out_dir / "1-convert" / "seqs" / "gene1.fa").write_text(">sp1\nA")

    from phyloai.cli.commands._run_pipeline import (
        _build_run_params, _build_run_checkpoint, _save_run_checkpoint,
    )
    params = _build_run_params(seq_dir, "supermatrix", "normal", 4, out_dir)
    ckpt = _build_run_checkpoint("phyloai run ...", params, mode="supermatrix", speed="normal")
    for s in ckpt["steps"]:
        s["output_dir"] = str(out_dir / {
            "convert": "1-convert", "align": "2-align", "trim": "3-trim",
            "filter_taper": "4-filter", "concat": "5-concat", "tree": "6-tree",
        }[s["name"]])
    ckpt["steps"][0]["status"] = "success"
    ckpt["steps"][1]["status"] = "running"   # align was interrupted
    _save_run_checkpoint(ckpt, out_dir / "run_checkpoint.json")

    r = _mock_step_result()

    def _mock_iqtree(*args, **kwargs) -> dict:
        out = kwargs.get("output_dir", Path("."))
        tf = Path(out) / "iqtree.treefile"
        tf.parent.mkdir(parents=True, exist_ok=True)
        tf.write_text("(sp1,sp2);")
        return {**r, "data": {"output": str(tf)}}

    mocked_align = MagicMock(return_value=r)

    with patch("phyloai.pretree.convert.convert_input", MagicMock(return_value=r)), \
         patch("phyloai.pretree.align.run_align", mocked_align) as mock_align, \
         patch("phyloai.pretree.trim.run_trim", return_value=r), \
         patch("phyloai.pretree.filter.run_taper", return_value=r), \
         patch("phyloai.pretree.concat.run_concat", return_value={**r, "key_results": {"total_length": 10, "n_taxa": 2}}), \
         patch("phyloai.tree.ml_iqtree.run_iqtree", side_effect=_mock_iqtree):
        result = runner.invoke(cli, [
            "run", "--seq-dir", str(seq_dir),
            "--output-dir", str(out_dir),
            "--resume",
        ])

    assert result.exit_code == 0, result.output
    assert mock_align.called
    call_kwargs = mock_align.call_args.kwargs
    assert call_kwargs["resume"] is True, "align should be called with resume=True"
    assert call_kwargs["overwrite"] is False, "align with resume should NOT have overwrite=True"


def test_tool_missing_exits_3(tmp_path: Path) -> None:
    """Missing tool should exit 3 (environment error), not exit 2."""
    runner = CliRunner()
    seq_dir = _make_seq_dir(tmp_path)
    out_dir = tmp_path / "run"

    with patch("phyloai.pretree.convert.convert_input",
               side_effect=FileNotFoundError("Tool 'mafft' not found")):
        result = runner.invoke(cli, [
            "run", "--seq-dir", str(seq_dir),
            "--output-dir", str(out_dir),
        ])

    assert result.exit_code == 3, f"expected exit 3, got {result.exit_code}"


def test_final_tree_missing_marks_tree_step_failed(tmp_path: Path) -> None:
    """When final tree is missing, the tree step in checkpoint should be 'failed'."""
    runner = CliRunner()
    seq_dir = _make_seq_dir(tmp_path)
    out_dir = tmp_path / "run"

    r = _mock_step_result()

    def _mock_convert(*args, **kwargs) -> dict:
        out = kwargs.get("output_dir", Path("."))
        (Path(out) / "seqs").mkdir(parents=True, exist_ok=True)
        (Path(out) / "seqs" / "gene1.fa").write_text(">sp1\nA")
        return r

    # iqtree succeeds but does NOT create the tree file
    def _mock_iqtree_no_file(*args, **kwargs) -> dict:
        return {**r, "data": {"output": str(Path(kwargs.get("output_dir", ".")) / "missing.treefile")}}

    with patch("phyloai.pretree.convert.convert_input", side_effect=_mock_convert), \
         patch("phyloai.pretree.align.run_align", return_value=r), \
         patch("phyloai.pretree.trim.run_trim", return_value=r), \
         patch("phyloai.pretree.filter.run_taper", return_value=r), \
         patch("phyloai.pretree.concat.run_concat", return_value={**r, "key_results": {"total_length": 10, "n_taxa": 2}}), \
         patch("phyloai.tree.ml_iqtree.run_iqtree", side_effect=_mock_iqtree_no_file):
        result = runner.invoke(cli, [
            "run", "--seq-dir", str(seq_dir),
            "--output-dir", str(out_dir),
        ])

    assert result.exit_code != 0

    # Check checkpoint: tree step should be "failed"
    ckpt_path = out_dir / "run_checkpoint.json"
    assert ckpt_path.exists()
    ckpt = json.loads(ckpt_path.read_text())
    tree_step = [s for s in ckpt["steps"] if s["name"] == "tree"][0]
    assert tree_step["status"] == "failed", (
        f"tree step should be 'failed' but is '{tree_step['status']}'"
    )


def test_supermatrix_result_has_matrix_fields(tmp_path: Path) -> None:
    """Supermatrix result.json should include matrix_length and matrix_taxa."""
    runner = CliRunner()
    seq_dir = _make_seq_dir(tmp_path)
    out_dir = tmp_path / "run"

    r = _mock_step_result()

    def _mock_convert(*args, **kwargs) -> dict:
        out = kwargs.get("output_dir", Path("."))
        (Path(out) / "seqs").mkdir(parents=True, exist_ok=True)
        (Path(out) / "seqs" / "gene1.fa").write_text(">sp1\nA")
        return r

    def _mock_iqtree(*args, **kwargs) -> dict:
        out = kwargs.get("output_dir", Path("."))
        tf = Path(out) / "iqtree.treefile"
        tf.parent.mkdir(parents=True, exist_ok=True)
        tf.write_text("(sp1,sp2);")
        return {**r, "data": {"output": str(tf)}}

    with patch("phyloai.pretree.convert.convert_input", side_effect=_mock_convert), \
         patch("phyloai.pretree.align.run_align", return_value=r), \
         patch("phyloai.pretree.trim.run_trim", return_value=r), \
         patch("phyloai.pretree.filter.run_taper", return_value=r), \
         patch("phyloai.pretree.concat.run_concat",
               return_value={**r, "key_results": {"total_length": 43820, "n_taxa": 52}}), \
         patch("phyloai.tree.ml_iqtree.run_iqtree", side_effect=_mock_iqtree):
        result = runner.invoke(cli, [
            "run", "--seq-dir", str(seq_dir),
            "--mode", "supermatrix", "--speed", "normal",
            "--output-dir", str(out_dir),
        ])

    assert result.exit_code == 0, result.output
    data = json.loads((out_dir / "result.json").read_text())
    kr = data["key_results"]
    assert kr["matrix_length"] == 43820
    assert kr["matrix_taxa"] == 52
