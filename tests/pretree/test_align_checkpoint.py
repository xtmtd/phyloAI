from __future__ import annotations

from pathlib import Path

from phyloai.core.checkpoint import Checkpoint, CheckpointTask


def test_verify_align_outputs_accepts_valid_aa(tmp_path: Path) -> None:
    from phyloai.pretree.align import verify_align_outputs

    aa = tmp_path / "gene1.fa"
    aa.write_text(">a\nMKT\n>b\nMKA\n")
    nt = tmp_path / "gene1.nt.fa"
    nt.write_text(">a\nATGAAGACT\n>b\nATGAAGGCT\n")

    assert verify_align_outputs(aa, nt) is True


def test_verify_align_outputs_accepts_missing_nt_when_none(tmp_path: Path) -> None:
    from phyloai.pretree.align import verify_align_outputs

    aa = tmp_path / "gene1.fa"
    aa.write_text(">a\nMKT\n>b\nMKA\n")

    assert verify_align_outputs(aa, None) is True


def test_verify_align_outputs_rejects_missing_nt(tmp_path: Path) -> None:
    from phyloai.pretree.align import verify_align_outputs

    aa = tmp_path / "gene1.fa"
    aa.write_text(">a\nMKT\n>b\nMKA\n")

    assert verify_align_outputs(aa, tmp_path / "missing.nt.fa") is False


def test_resolved_align_params_includes_required_keys() -> None:
    from phyloai.pretree.align import _resolved_align_params

    params = _resolved_align_params(
        seq_dir=Path("raw"),
        output_dir=Path("out"),
        method="linsi",
        resolved_seq_type="AA",
        backtrans=False,
        nt_dir=None,
        threads=8,
        tool_args=None,
        mafft_executable="/usr/bin/mafft",
        magus_executable="magus",
        trimal_executable="trimal",
        quiet=True,
    )

    assert params == {
        "seq_dir": "raw",
        "output_dir": "out",
        "method": "linsi",
        "seq_type": "AA",
        "backtrans": False,
        "nt_dir": None,
        "threads": 8,
        "tool_args": None,
        "mafft_executable": "/usr/bin/mafft",
        "magus_executable": "magus",
        "trimal_executable": "trimal",
        "quiet": True,
    }


def test_resolved_align_params_excludes_mode_flags() -> None:
    from phyloai.pretree.align import _resolved_align_params

    params = _resolved_align_params(
        seq_dir=Path("raw"),
        output_dir=Path("out"),
        method="linsi",
        resolved_seq_type="AA",
        backtrans=False,
        nt_dir=None,
        threads=1,
        tool_args=None,
        mafft_executable="mafft",
        magus_executable="magus",
        trimal_executable="trimal",
        quiet=False,
    )

    assert "overwrite" not in params
    assert "resume" not in params
    assert "dry_run" not in params


def test_build_initial_checkpoint(tmp_path: Path) -> None:
    from phyloai.pretree.align import _resolved_align_params
    from phyloai.pretree.checkpoint_helpers import build_initial_checkpoint

    inputs = [tmp_path / "gene1.fa", tmp_path / "gene2.fa"]
    for path in inputs:
        path.write_text(">a\nMKT\n")

    params = _resolved_align_params(
        seq_dir=tmp_path,
        output_dir=tmp_path / "out",
        method="linsi",
        resolved_seq_type="AA",
        backtrans=False,
        nt_dir=None,
        threads=2,
        tool_args=None,
        mafft_executable="mafft",
        magus_executable="magus",
        trimal_executable="trimal",
        quiet=False,
    )

    checkpoint = build_initial_checkpoint(
        step="pretree.align",
        command="phyloai pretree align",
        params=params,
        inputs=inputs,
        output_for=lambda p: tmp_path / "out" / f"{p.stem}.fa",
        nt_output_for=lambda p: None,
    )

    assert checkpoint.step == "pretree.align"
    assert checkpoint.status == "running"
    assert [task.task_id for task in checkpoint.tasks] == ["gene1", "gene2"]


def test_mark_task_updates_status() -> None:
    from phyloai.pretree.checkpoint_helpers import mark_task

    checkpoint = Checkpoint(
        schema_version=1,
        step="pretree.align",
        command="phyloai pretree align",
        status="running",
        params_hash="sha256:abc",
        params={},
        started_at="2026-06-12T10:00:00+00:00",
        updated_at="2026-06-12T10:00:00+00:00",
        completed_at=None,
        tasks=[CheckpointTask("g1", "pending", "raw/g1.fa", {"aa": "out/g1.fa", "nt": None})],
    )

    mark_task(checkpoint, "g1", status="running", reason=None)

    assert checkpoint.tasks[0].status == "running"
    assert checkpoint.tasks[0].attempts == 1
    assert checkpoint.tasks[0].updated_at is not None


def test_plan_resume_marks_invalid_success_for_rerun(tmp_path: Path) -> None:
    from phyloai.pretree.align import verify_align_outputs
    from phyloai.pretree.checkpoint_helpers import plan_resume

    checkpoint = Checkpoint(
        schema_version=1,
        step="pretree.align",
        command="phyloai pretree align",
        status="running",
        params_hash="sha256:abc",
        params={},
        started_at="2026-06-12T10:00:00+00:00",
        updated_at="2026-06-12T10:00:00+00:00",
        completed_at=None,
        tasks=[
            CheckpointTask("g1", "success", "raw/g1.fa", {"aa": str(tmp_path / 'missing.fa'), "nt": None}),
            CheckpointTask("g2", "failed", "raw/g2.fa", {"aa": str(tmp_path / 'g2.fa'), "nt": None}),
            CheckpointTask("g3", "pending", "raw/g3.fa", {"aa": str(tmp_path / 'g3.fa'), "nt": None}),
        ],
    )

    to_run, skipped = plan_resume(checkpoint, verify_align_outputs)

    assert sorted(to_run) == ["g1", "g2", "g3"]
    assert skipped == []


def test_reconstruct_align_result_aggregates_states(tmp_path: Path) -> None:
    from phyloai.pretree.align import reconstruct_align_result

    aa = tmp_path / "g1.fa"
    aa.write_text(">a\nMKT\n>b\nMKA\n")
    nt = tmp_path / "g1.nt.fa"
    nt.write_text(">a\nATGAAGACT\n>b\nATGAAGGCT\n")

    checkpoint = Checkpoint(
        schema_version=1,
        step="pretree.align",
        command="phyloai pretree align",
        status="success",
        params_hash="sha256:abc",
        params={"method": "linsi", "backtrans": True},
        started_at="2026-06-12T10:00:00+00:00",
        updated_at="2026-06-12T10:10:00+00:00",
        completed_at="2026-06-12T10:10:00+00:00",
        tasks=[
            CheckpointTask("g1", "success", "raw/g1.fa", {"aa": str(aa), "nt": str(nt)}, attempts=1),
        ],
    )

    payload = reconstruct_align_result(
        checkpoint=checkpoint,
        params=checkpoint.params,
        tool_versions={"mafft": "7.526"},
        wall_time=3.2,
        skipped_inputs=[{"path": "raw/bad.fa", "reason": "empty file"}],
        scan_warnings=[],
    )

    assert payload["status"] == "success"
    assert payload["key_results"]["n_aligned"] == 1
    assert payload["key_results"]["n_skipped"] == 1
    assert payload["data"]["files"][0]["alignment_length"] == 3
    assert payload["data"]["files"][0]["n_taxa"] == 2
