"""Tree module checkpoint helpers."""

from __future__ import annotations

import datetime as _dt
from pathlib import Path
from typing import Any, Callable

from Bio import Phylo

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
    trees_dir: Path,
    logs_dir: Path,
) -> Checkpoint:
    now = _utc_now_iso()
    tasks = [
        CheckpointTask(
            task_id=inp.stem,
            status="pending",
            input=str(inp),
            outputs={
                "tree": str(trees_dir / f"{inp.stem}.tre"),
                "log": str(logs_dir / f"{inp.stem}.log"),
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


def resume_verifier() -> Callable[[Path], bool]:
    def _verify(tree_path: Path) -> bool:
        if not tree_path.exists() or tree_path.stat().st_size == 0:
            return False
        try:
            Phylo.read(str(tree_path), "newick")
            return True
        except Exception:
            return False

    return _verify


def plan_resume(checkpoint: Checkpoint) -> tuple[list[str], list[str]]:
    to_run: list[str] = []
    skipped: list[str] = []
    verifier = resume_verifier()

    for task in checkpoint.tasks:
        if task.status in {"pending", "running", "failed"}:
            to_run.append(task.task_id)
        elif task.status == "success":
            tree_path = Path(task.outputs["tree"]) if task.outputs.get("tree") else None
            if tree_path is not None and verifier(tree_path):
                skipped.append(task.task_id)
            else:
                to_run.append(task.task_id)
        else:
            skipped.append(task.task_id)

    return to_run, skipped


def resume_verifier_iqtree(validate_tree: bool = True) -> Callable[[Path], bool]:
    """Return verifier for IQ-TREE output. validate_tree=True checks .treefile with Bio.Phylo.
    validate_tree=False (MF mode) checks .iqtree exists and is non-empty."""
    if validate_tree:
        def _verify(tree_path: Path) -> bool:
            if not tree_path.exists() or tree_path.stat().st_size == 0:
                return False
            try:
                trees = Phylo.parse(str(tree_path), "newick")
                next(trees)
                for _ in trees:
                    pass
                return True
            except Exception:
                return False
        return _verify
    else:
        def _verify(iqtree_path: Path) -> bool:
            if not iqtree_path.exists() or iqtree_path.stat().st_size == 0:
                return False
            return True
        return _verify


def build_initial_iqtree_checkpoint(
    *,
    step: str,
    command: str,
    params: dict[str, Any],
    inputs: list[Path],
    trees_dir: Path,
    logs_dir: Path,
) -> Checkpoint:
    """Like build_initial_checkpoint but records both tree and iqtree outputs per task."""
    now = _utc_now_iso()
    tasks = [
        CheckpointTask(
            task_id=inp.stem,
            status="pending",
            input=str(inp),
            outputs={
                "tree": str(trees_dir / f"{inp.stem}.treefile"),
                "iqtree": str(logs_dir / f"{inp.stem}.iqtree"),
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


def plan_resume_iqtree(
    checkpoint: Checkpoint,
    is_mf_only: bool = False,
) -> tuple[list[str], list[str]]:
    """Like plan_resume but for IQ-TREE. Uses .treefile for tree-producing,
    .iqtree for MF-only workflows."""
    to_run: list[str] = []
    skipped: list[str] = []
    verifier = resume_verifier_iqtree(validate_tree=not is_mf_only)

    for task in checkpoint.tasks:
        if task.status in {"pending", "running", "failed"}:
            to_run.append(task.task_id)
        elif task.status == "success":
            output_key = "iqtree" if is_mf_only else "tree"
            output_path_str = task.outputs.get(output_key)
            if output_path_str and verifier(Path(output_path_str)):
                skipped.append(task.task_id)
            else:
                to_run.append(task.task_id)
        else:
            skipped.append(task.task_id)

    return to_run, skipped
