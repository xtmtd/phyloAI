"""pretree-specific checkpoint helpers."""

from __future__ import annotations

import datetime as _dt
from collections.abc import Callable
from pathlib import Path
from typing import Any

from phyloai.core.checkpoint import (
    CHECKPOINT_SCHEMA_VERSION,
    Checkpoint,
    CheckpointTask,
    canonical_params_hash,
)


def _utc_now_iso() -> str:
    return _dt.datetime.now(_dt.timezone.utc).isoformat(timespec="seconds")


def build_initial_checkpoint(
    *,
    step: str,
    command: str,
    params: dict[str, Any],
    inputs: list[Path],
    output_for: Callable[[Path], Path],
    nt_output_for: Callable[[Path], Path | None],
) -> Checkpoint:
    now = _utc_now_iso()
    tasks = [
        CheckpointTask(
            task_id=inp.stem,
            status="pending",
            input=str(inp),
            outputs={
                "aa": str(output_for(inp)),
                "nt": str(nt_out) if (nt_out := nt_output_for(inp)) is not None else None,
            },
        )
        for inp in inputs
    ]
    return Checkpoint(
        schema_version=CHECKPOINT_SCHEMA_VERSION,
        step=step,
        command=command,
        status="running",
        params_hash=canonical_params_hash(params),
        params=params,
        started_at=now,
        updated_at=now,
        completed_at=None,
        tasks=tasks,
    )


def mark_task(
    checkpoint: Checkpoint,
    task_id: str,
    *,
    status: str,
    reason: str | None = None,
) -> CheckpointTask:
    for task in checkpoint.tasks:
        if task.task_id == task_id:
            task.status = status
            task.reason = reason
            task.attempts += 1
            task.updated_at = _utc_now_iso()
            checkpoint.touch()
            return task
    raise KeyError(f"Task {task_id!r} not found in checkpoint")


def resume_verifier(
    verify_outputs: Callable[[Path, Path | None], bool],
) -> Callable[[CheckpointTask], bool]:
    def _verifier(task: CheckpointTask) -> bool:
        aa = Path(task.outputs["aa"]) if task.outputs.get("aa") else None
        nt = Path(task.outputs["nt"]) if task.outputs.get("nt") else None
        if aa is None:
            return False
        return verify_outputs(aa, nt)

    return _verifier


def plan_resume(
    checkpoint: Checkpoint,
    verify_outputs: Callable[[Path, Path | None], bool],
) -> tuple[list[str], list[str]]:
    to_run: list[str] = []
    skipped: list[str] = []
    for task in checkpoint.tasks:
        if task.status in {"pending", "running", "failed"}:
            to_run.append(task.task_id)
        elif task.status == "success":
            aa = Path(task.outputs["aa"]) if task.outputs.get("aa") else None
            nt = Path(task.outputs["nt"]) if task.outputs.get("nt") else None
            if aa is not None and verify_outputs(aa, nt):
                skipped.append(task.task_id)
            else:
                to_run.append(task.task_id)
        else:
            skipped.append(task.task_id)
    return to_run, skipped
