from __future__ import annotations

import json
from pathlib import Path

import pytest


def test_canonical_params_hash_is_stable() -> None:
    from phyloai.core.checkpoint import canonical_params_hash

    a = {"b": 1, "a": 2, "c": [1, 2, 3]}
    b = {"a": 2, "b": 1, "c": [1, 2, 3]}
    assert canonical_params_hash(a) == canonical_params_hash(b)


def test_canonical_params_hash_changes_with_values() -> None:
    from phyloai.core.checkpoint import canonical_params_hash

    assert canonical_params_hash({"a": 1}) != canonical_params_hash({"a": 2})


def test_save_checkpoint_atomic_writes_valid_json(tmp_path: Path) -> None:
    from phyloai.core.checkpoint import Checkpoint, CheckpointTask, save_checkpoint_atomic

    checkpoint = Checkpoint(
        schema_version=1,
        step="pretree.align",
        command="phyloai pretree align --seq-dir /data --output-dir /out --method linsi --seq-type AA --threads 4",
        status="running",
        params_hash="sha256:abc",
        params={"method": "linsi"},
        started_at="2026-06-12T10:00:00+00:00",
        updated_at="2026-06-12T10:00:00+00:00",
        completed_at=None,
        tasks=[
            CheckpointTask(
                task_id="gene1",
                status="pending",
                input="raw/gene1.fa",
                outputs={"aa": "out/gene1.fa", "nt": None},
            )
        ],
    )

    path = tmp_path / "checkpoint.json"
    save_checkpoint_atomic(checkpoint, path)

    assert path.exists()
    assert not path.with_suffix(".json.tmp").exists()
    loaded = json.loads(path.read_text())
    assert loaded["step"] == "pretree.align"
    assert loaded["tasks"][0]["task_id"] == "gene1"


def test_load_checkpoint_rejects_malformed_json(tmp_path: Path) -> None:
    from phyloai.core.checkpoint import load_checkpoint

    path = tmp_path / "checkpoint.json"
    path.write_text("{not json")

    with pytest.raises(ValueError, match="malformed JSON"):
        load_checkpoint(path)


def test_validate_resume_params_rejects_mismatch() -> None:
    from phyloai.core.checkpoint import (
        Checkpoint,
        canonical_params_hash,
        validate_resume_params,
    )

    params = {"method": "linsi", "threads": 4, "quiet": False}
    checkpoint = Checkpoint(
        schema_version=1,
        step="pretree.align",
        command="phyloai pretree align --seq-dir /data --output-dir /out --method linsi --seq-type AA --threads 4",
        status="running",
        params_hash=canonical_params_hash(params),
        params=params,
        started_at="2026-06-12T10:00:00+00:00",
        updated_at="2026-06-12T10:00:00+00:00",
        completed_at=None,
        tasks=[],
    )

    with pytest.raises(ValueError, match="Resume parameter mismatch"):
        validate_resume_params(
            checkpoint,
            {"method": "fftns1", "threads": 4, "quiet": False},
            step="pretree.align",
        )


def test_load_checkpoint_rejects_unsupported_schema(tmp_path: Path) -> None:
    from phyloai.core.checkpoint import load_checkpoint

    path = tmp_path / "checkpoint.json"
    path.write_text(
        json.dumps(
            {
                "schema_version": 999,
                "step": "pretree.align",
                "command": "phyloai pretree align",
                "status": "running",
                "params_hash": "sha256:abc",
                "params": {},
                "started_at": "2026-06-12T10:00:00+00:00",
                "updated_at": "2026-06-12T10:00:00+00:00",
                "completed_at": None,
                "tasks": [],
            }
        )
    )

    with pytest.raises(ValueError, match="Unsupported checkpoint schema_version"):
        load_checkpoint(path)


def test_summarize_resume_tasks_counts_invalid_successes(tmp_path: Path) -> None:
    from phyloai.core.checkpoint import Checkpoint, CheckpointTask, summarize_resume_tasks

    ok = tmp_path / "ok.fa"
    ok.write_text(">a\nMKT\n>b\nMKA\n")

    checkpoint = Checkpoint(
        schema_version=1,
        step="pretree.align",
        command="phyloai pretree align --seq-dir /data --output-dir /out --method linsi --seq-type AA --threads 4",
        status="running",
        params_hash="sha256:abc",
        params={},
        started_at="2026-06-12T10:00:00+00:00",
        updated_at="2026-06-12T10:00:00+00:00",
        completed_at=None,
        tasks=[
            CheckpointTask("g1", "success", "raw/g1.fa", {"aa": str(ok), "nt": None}),
            CheckpointTask("g2", "success", "raw/g2.fa", {"aa": str(tmp_path / 'missing.fa'), "nt": None}),
            CheckpointTask("g3", "failed", "raw/g3.fa", {"aa": str(tmp_path / 'g3.fa'), "nt": None}),
        ],
    )

    summary = summarize_resume_tasks(
        checkpoint,
        lambda task: Path(task.outputs["aa"]).exists(),
    )

    assert summary == {"skip": 1, "rerun": 1, "invalid": 1}


def test_validate_resume_params_excludes_control_flags() -> None:
    from phyloai.core.checkpoint import (
        Checkpoint,
        canonical_params_hash,
        validate_resume_params,
    )

    params_full = {
        "msa_dir": "/data",
        "tool": "trimal",
        "threads": 4,
        "resume": True,
        "overwrite": False,
        "dry_run": False,
        "quiet": False,
    }
    stored_params = dict(params_full)
    stored_params["resume"] = False

    cp = Checkpoint(
        schema_version=1,
        step="pretree.trim",
        command="",
        status="running",
        params_hash=canonical_params_hash(stored_params),
        params=stored_params,
        started_at="",
        updated_at="",
        completed_at=None,
        tasks=[],
    )

    validate_resume_params(cp, params_full, step="pretree.trim")
