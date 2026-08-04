"""Shared checkpoint helpers for resumable commands."""

from __future__ import annotations

import datetime as _dt
import hashlib
import json
import os
from collections.abc import Callable
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any


CHECKPOINT_SCHEMA_VERSION = 1


def _utc_now_iso() -> str:
    return _dt.datetime.now(_dt.timezone.utc).isoformat(timespec="seconds")


def canonical_params_hash(params: dict[str, Any]) -> str:
    payload = json.dumps(
        params,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=False,
    )
    return f"sha256:{hashlib.sha256(payload.encode('utf-8')).hexdigest()}"


@dataclass(eq=True)
class CheckpointTask:
    task_id: str
    status: str
    input: str
    outputs: dict[str, str | None] = field(default_factory=dict)
    attempts: int = 0
    reason: str | None = None
    updated_at: str | None = None

    def to_dict(self) -> dict[str, Any]:
        return {
            "task_id": self.task_id,
            "status": self.status,
            "input": self.input,
            "outputs": dict(self.outputs),
            "attempts": self.attempts,
            "reason": self.reason,
            "updated_at": self.updated_at,
        }

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> "CheckpointTask":
        return cls(
            task_id=str(data["task_id"]),
            status=str(data["status"]),
            input=str(data["input"]),
            outputs={str(k): (None if v is None else str(v)) for k, v in data.get("outputs", {}).items()},
            attempts=int(data.get("attempts", 0)),
            reason=None if data.get("reason") is None else str(data.get("reason")),
            updated_at=None if data.get("updated_at") is None else str(data.get("updated_at")),
        )


@dataclass(eq=True)
class Checkpoint:
    schema_version: int
    step: str
    command: str
    status: str
    params_hash: str
    params: dict[str, Any]
    started_at: str
    updated_at: str
    completed_at: str | None
    tasks: list[CheckpointTask]
    original_msa_fingerprint: str | None = None

    def touch(self) -> str:
        self.updated_at = _utc_now_iso()
        return self.updated_at

    def to_dict(self) -> dict[str, Any]:
        return {
            "schema_version": self.schema_version,
            "step": self.step,
            "command": self.command,
            "status": self.status,
            "params_hash": self.params_hash,
            "params": self.params,
            "started_at": self.started_at,
            "updated_at": self.updated_at,
            "completed_at": self.completed_at,
            "tasks": [task.to_dict() for task in self.tasks],
            "original_msa_fingerprint": self.original_msa_fingerprint,
        }

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> "Checkpoint":
        return cls(
            schema_version=int(data["schema_version"]),
            step=str(data["step"]),
            command=str(data["command"]),
            status=str(data["status"]),
            params_hash=str(data["params_hash"]),
            params=dict(data.get("params", {})),
            started_at=str(data["started_at"]),
            updated_at=str(data["updated_at"]),
            completed_at=None if data.get("completed_at") is None else str(data.get("completed_at")),
            tasks=[CheckpointTask.from_dict(task) for task in data.get("tasks", [])],
            original_msa_fingerprint=(
                None if data.get("original_msa_fingerprint") is None
                else str(data.get("original_msa_fingerprint"))
            ),
        )


def save_checkpoint_atomic(
    checkpoint: Checkpoint, path: Path, *, fsync: bool = False
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = path.with_suffix(path.suffix + ".tmp")
    payload = json.dumps(checkpoint.to_dict(), indent=2, ensure_ascii=False)
    with open(tmp_path, "w", encoding="utf-8") as fh:
        fh.write(payload)
        fh.flush()
        if fsync:
            os.fsync(fh.fileno())
    os.replace(tmp_path, path)


def load_checkpoint(path: Path) -> Checkpoint:
    if not path.exists():
        raise FileNotFoundError(
            f"No checkpoint found at {path}. Use --overwrite to start fresh."
        )
    try:
        with open(path, "r", encoding="utf-8") as fh:
            data = json.load(fh)
    except json.JSONDecodeError as exc:
        raise ValueError(
            f"Checkpoint file {path} is malformed JSON: {exc}. Use --overwrite to start fresh."
        ) from exc

    version = int(data.get("schema_version", -1))
    if version != CHECKPOINT_SCHEMA_VERSION:
        raise ValueError(
            f"Unsupported checkpoint schema_version {version}; expected {CHECKPOINT_SCHEMA_VERSION}. "
            "Use --overwrite to start fresh."
        )

    return Checkpoint.from_dict(data)


_RESUME_EXCLUDED_KEYS = frozenset({"resume", "overwrite", "dry_run", "quiet",
                                   "_command"})


def _clean_params_for_resume(params: dict[str, Any]) -> dict[str, Any]:
    return {k: v for k, v in params.items() if k not in _RESUME_EXCLUDED_KEYS}


def validate_resume_params(
    checkpoint: Checkpoint,
    params: dict[str, Any],
    *,
    step: str | None = None,
) -> None:
    if step is not None and checkpoint.step != step:
        raise ValueError(
            f"Checkpoint step is {checkpoint.step!r}, current command step is {step!r}. "
            "Use --overwrite to start fresh."
        )
    current_hash = canonical_params_hash(_clean_params_for_resume(params))
    stored_hash = canonical_params_hash(_clean_params_for_resume(checkpoint.params))
    if current_hash != stored_hash:
        raise ValueError(
            "Resume parameter mismatch: current invocation does not match the checkpoint. "
            "To change parameters, restart with --overwrite."
        )


def summarize_resume_tasks(
    checkpoint: Checkpoint,
    verifier: Callable[[CheckpointTask], bool],
) -> dict[str, int]:
    summary = {"skip": 0, "rerun": 0, "invalid": 0}
    for task in checkpoint.tasks:
        if task.status in {"pending", "running", "failed"}:
            summary["rerun"] += 1
        elif task.status == "success":
            if verifier(task):
                summary["skip"] += 1
            else:
                summary["invalid"] += 1
        else:
            summary["skip"] += 1
    return summary
