# Checkpoint and Resume Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a shared checkpoint and `--resume` mechanism to PhyloAI, and adopt it in `phyloai pretree align` so interrupted long runs can be continued without redoing completed work.

**Architecture:** Add a small `phyloai/core/checkpoint.py` module that owns the on-disk `checkpoint.json` format, atomic writes, parameter hashing, and task-state summarization. `phyloai pretree align` gains a `--resume` flag, treats one alignment task per input gene, and rebuilds its final `result.json` from checkpoint task states plus the validated output files. Strict resume requires exact parameter equality, skips verified successes, and reruns failed or interrupted tasks.

**Tech Stack:** Python 3.10+, dataclasses, `json`, `hashlib`, `tempfile`, `os.replace`, `pytest`, Click CliRunner, Biopython `SeqIO`, `phyloai.core.sequence_output_validation`.

---

## File Structure

- Create: `phyloai/core/checkpoint.py` — `Checkpoint`, `CheckpointTask`, hash helper, atomic save/load, resume validation
- Modify: `phyloai/core/__init__.py` — export `Checkpoint`, `CheckpointTask`, and resume helpers
- Create: `tests/core/test_checkpoint.py` — unit tests for the shared checkpoint module
- Modify: `phyloai/pretree/align.py` — checkpoint creation, per-task updates, `--resume` flow, result reconstruction
- Create: `tests/pretree/test_align_checkpoint.py` — resume-specific align tests
- Modify: `phyloai/cli/commands/pretree.py` — add `--resume` flag for `align`, mutual-exclusion check, dry-run resume summary
- Create: `docs/commands/pretree-align.md` — user-facing resume behavior documentation (if not already present)
- Modify: `docs/commands/pretree-align.md` (or create if missing) — `--resume` and checkpoint section

---

## Required Design Decisions Applied In This Plan

- Resume is **explicit**: `--resume` must be supplied; no auto-resume on non-empty directories.
- Resume is **strict**: every resolved parameter (analysis and run-control, including `--quiet`) must match `params_hash` exactly. Mismatch exits 1.
- Failed, interrupted, and unstarted tasks are **retried by default**; success tasks are **skipped only after output verification**.
- The `core` checkpoint module is **command-agnostic**. `pretree align` provides its own verifier and result reconstruction.
- `result.json` remains the final structured command result; `checkpoint.json` only stores the minimum state needed to continue.
- `--overwrite` and `--resume` are **mutually exclusive**.
- The shared module supports a stable **schema version 1**. Loading a higher version exits with a clear error.

---

## Task 1: Core Checkpoint Dataclasses and Hash Helper

**Files:**
- Create: `phyloai/core/checkpoint.py`
- Create: `tests/core/test_checkpoint.py`

- [ ] **Step 1: Write failing tests for dataclasses and parameter hash**

Create `tests/core/test_checkpoint.py`:

```python
from __future__ import annotations

from pathlib import Path

import pytest


def test_canonical_params_hash_is_stable() -> None:
    from phyloai.core.checkpoint import canonical_params_hash

    a = {"b": 1, "a": 2, "c": [1, 2, 3]}
    b = {"a": 2, "b": 1, "c": (1, 2, 3)}
    assert canonical_params_hash(a) == canonical_params_hash(b)


def test_canonical_params_hash_changes_with_values() -> None:
    from phyloai.core.checkpoint import canonical_params_hash

    assert canonical_params_hash({"a": 1}) != canonical_params_hash({"a": 2})


def test_checkpoint_task_to_dict_drops_none_reason() -> None:
    from phyloai.core.checkpoint import CheckpointTask

    task = CheckpointTask(
        task_id="gene1",
        status="success",
        input="raw/gene1.fa",
        outputs={"aa": "out/gene1.fa", "nt": None},
    )
    d = task.to_dict()
    assert d["task_id"] == "gene1"
    assert d["status"] == "success"
    assert d["outputs"] == {"aa": "out/gene1.fa", "nt": None}
    assert d["reason"] is None
    assert d["attempts"] == 0


def test_checkpoint_to_dict_contains_required_keys() -> None:
    import datetime as _dt
    from phyloai.core.checkpoint import Checkpoint, CheckpointTask

    cp = Checkpoint(
        schema_version=1,
        step="pretree.align",
        command="phyloai pretree align",
        status="running",
        params_hash="sha256:abc",
        params={"method": "linsi"},
        started_at="2026-06-12T10:00:00",
        updated_at="2026-06-12T10:00:00",
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
    d = cp.to_dict()
    for key in [
        "schema_version",
        "step",
        "command",
        "status",
        "params_hash",
        "params",
        "started_at",
        "updated_at",
        "completed_at",
        "tasks",
    ]:
        assert key in d
    assert d["schema_version"] == 1
    assert d["tasks"][0]["task_id"] == "gene1"


def test_checkpoint_from_dict_round_trip() -> None:
    from phyloai.core.checkpoint import Checkpoint, CheckpointTask

    original = Checkpoint(
        schema_version=1,
        step="pretree.align",
        command="phyloai pretree align",
        status="success",
        params_hash="sha256:abc",
        params={"method": "linsi"},
        started_at="2026-06-12T10:00:00",
        updated_at="2026-06-12T10:30:00",
        completed_at="2026-06-12T10:30:00",
        tasks=[
            CheckpointTask(
                task_id="gene1",
                status="success",
                input="raw/gene1.fa",
                outputs={"aa": "out/gene1.fa", "nt": None},
                attempts=1,
            )
        ],
    )
    reloaded = Checkpoint.from_dict(original.to_dict())
    assert reloaded == original
```

- [ ] **Step 2: Run tests to verify they fail**

```bash
pytest tests/core/test_checkpoint.py -v 2>&1 | head -20
```

Expected: FAIL with `ModuleNotFoundError: No module named 'phyloai.core.checkpoint'`

- [ ] **Step 3: Implement dataclasses and hash helper**

Create `phyloai/core/checkpoint.py`:

```python
"""Shared checkpoint and resume helpers for long-running PhyloAI commands."""

from __future__ import annotations

import datetime as _dt
import hashlib
import json
from collections.abc import Callable
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

CHECKPOINT_SCHEMA_VERSION = 1


def _utc_now_iso() -> str:
    return _dt.datetime.now(_dt.timezone.utc).isoformat(timespec="seconds")


def canonical_params_hash(params: dict[str, Any]) -> str:
    """Compute a stable SHA-256 hash of a resolved parameter dictionary."""
    payload = json.dumps(
        params,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=False,
        default=str,
    )
    digest = hashlib.sha256(payload.encode("utf-8")).hexdigest()
    return f"sha256:{digest}"


@dataclass
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
            reason=data.get("reason"),
            updated_at=data.get("updated_at"),
        )


@dataclass
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

    def touch(self) -> None:
        self.updated_at = _utc_now_iso()

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
            completed_at=data.get("completed_at"),
            tasks=[CheckpointTask.from_dict(t) for t in data.get("tasks", [])],
        )


def save_checkpoint_atomic(checkpoint: Checkpoint, path: Path) -> None:
    """Atomically write the checkpoint to ``path`` as JSON."""
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = json.dumps(checkpoint.to_dict(), indent=2, ensure_ascii=False)
    tmp = path.with_suffix(path.suffix + ".tmp")
    with open(tmp, "w", encoding="utf-8") as fh:
        fh.write(payload)
        fh.flush()
        import os

        os.fsync(fh.fileno())
    import os

    os.replace(tmp, path)


def load_checkpoint(path: Path) -> Checkpoint:
    """Load a checkpoint from ``path`` and validate the schema version."""
    if not path.exists():
        raise FileNotFoundError(
            f"No checkpoint found at {path}. Use --overwrite to start fresh."
        )
    try:
        with open(path, "r", encoding="utf-8") as fh:
            data = json.load(fh)
    except json.JSONDecodeError as exc:
        raise ValueError(
            f"Checkpoint file {path} is malformed JSON: {exc}. "
            "Use --overwrite to start fresh."
        ) from exc
    version = int(data.get("schema_version", -1))
    if version != CHECKPOINT_SCHEMA_VERSION:
        raise ValueError(
            f"Unsupported checkpoint schema_version {version}; "
            f"expected {CHECKPOINT_SCHEMA_VERSION}. Use --overwrite to start fresh."
        )
    return Checkpoint.from_dict(data)


def validate_resume_params(
    checkpoint: Checkpoint,
    params: dict[str, Any],
    *,
    step: str | None = None,
) -> None:
    """Verify that the current invocation parameters match the checkpoint."""
    if step is not None and checkpoint.step != step:
        raise ValueError(
            f"Checkpoint step is {checkpoint.step!r}, current command step is "
            f"{step!r}. Use --overwrite to start fresh."
        )
    if canonical_params_hash(params) != checkpoint.params_hash:
        raise ValueError(
            "Resume parameter mismatch: current invocation does not match the "
            "checkpoint. To change parameters, restart with --overwrite."
        )


def summarize_resume_tasks(
    checkpoint: Checkpoint,
    verifier: Callable[[CheckpointTask], bool],
) -> dict[str, int]:
    """Classify checkpoint tasks into skip/rerun/invalid buckets."""
    skip = 0
    rerun = 0
    invalid = 0
    for task in checkpoint.tasks:
        if task.status in ("pending", "failed", "running"):
            rerun += 1
        elif task.status == "success":
            if verifier(task):
                skip += 1
            else:
                invalid += 1
        else:
            # skipped tasks keep their status unless the command reclassifies them
            skip += 1
    return {"skip": skip, "rerun": rerun, "invalid": invalid}
```

- [ ] **Step 4: Run tests to verify they pass**

```bash
pytest tests/core/test_checkpoint.py -v
```

Expected: all PASS

- [ ] **Step 5: Commit**

```bash
git add phyloai/core/checkpoint.py tests/core/test_checkpoint.py
git commit -m "feat(core): add shared checkpoint module and parameter hash helper"
```

---

## Task 2: Core Package Exports for Checkpoint

**Files:**
- Modify: `phyloai/core/__init__.py`
- Create: `tests/core/test_core_exports.py`

- [ ] **Step 1: Write failing test for core exports**

Create `tests/core/test_core_exports.py`:

```python
def test_core_exposes_checkpoint_helpers() -> None:
    from phyloai.core import (
        Checkpoint,
        CheckpointTask,
        canonical_params_hash,
        load_checkpoint,
        save_checkpoint_atomic,
        validate_resume_params,
    )

    assert Checkpoint is not None
    assert CheckpointTask is not None
    assert callable(canonical_params_hash)
    assert callable(load_checkpoint)
    assert callable(save_checkpoint_atomic)
    assert callable(validate_resume_params)
```

- [ ] **Step 2: Run test to verify it fails**

```bash
pytest tests/core/test_core_exports.py -v
```

Expected: FAIL with `ImportError: cannot import name 'Checkpoint' from 'phyloai.core'`

- [ ] **Step 3: Update core package exports**

Replace `phyloai/core/__init__.py` with:

```python
"""Core infrastructure for PhyloAI."""

from phyloai.core.checkpoint import (
    CHECKPOINT_SCHEMA_VERSION,
    Checkpoint,
    CheckpointTask,
    canonical_params_hash,
    load_checkpoint,
    save_checkpoint_atomic,
    validate_resume_params,
)
from phyloai.core.schema import MSACollection, RunRecord, ToolResult, TreeSet
from phyloai.core.env import ToolEnv
from phyloai.core.runner import Runner
from phyloai.core.formats import FormatConverter
from phyloai.core.logger import StepLogger

__all__ = [
    "MSACollection", "TreeSet", "RunRecord", "ToolResult",
    "ToolEnv", "Runner", "FormatConverter", "StepLogger",
    "Checkpoint", "CheckpointTask",
    "CHECKPOINT_SCHEMA_VERSION",
    "canonical_params_hash",
    "load_checkpoint",
    "save_checkpoint_atomic",
    "validate_resume_params",
]
```

- [ ] **Step 4: Run test to verify it passes**

```bash
pytest tests/core/test_core_exports.py -v
```

Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add phyloai/core/__init__.py tests/core/test_core_exports.py
git commit -m "feat(core): export checkpoint helpers from core package"
```

---

## Task 3: Pretree Align Output Verifier

**Files:**
- Modify: `phyloai/pretree/align.py`
- Create: `tests/pretree/test_align_checkpoint.py`

- [ ] **Step 1: Write failing tests for the align output verifier**

Create `tests/pretree/test_align_checkpoint.py`:

```python
from __future__ import annotations

from pathlib import Path


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


def test_verify_align_outputs_rejects_missing_aa(tmp_path: Path) -> None:
    from phyloai.pretree.align import verify_align_outputs

    nt = tmp_path / "gene1.nt.fa"
    nt.write_text(">a\nATGAAGACT\n")
    assert verify_align_outputs(tmp_path / "missing.fa", nt) is False


def test_verify_align_outputs_rejects_empty_aa(tmp_path: Path) -> None:
    from phyloai.pretree.align import verify_align_outputs

    aa = tmp_path / "gene1.fa"
    aa.write_text("")
    assert verify_align_outputs(aa, None) is False


def test_verify_align_outputs_rejects_unaligned_aa(tmp_path: Path) -> None:
    from phyloai.pretree.align import verify_align_outputs

    aa = tmp_path / "gene1.fa"
    aa.write_text(">a\nMKT\n>b\nMKA\nM\n")
    assert verify_align_outputs(aa, None) is False


def test_verify_align_outputs_rejects_missing_nt(tmp_path: Path) -> None:
    from phyloai.pretree.align import verify_align_outputs

    aa = tmp_path / "gene1.fa"
    aa.write_text(">a\nMKT\n>b\nMKA\n")
    nt = tmp_path / "missing.nt.fa"
    assert verify_align_outputs(aa, nt) is False
```

- [ ] **Step 2: Run tests to verify they fail**

```bash
pytest tests/pretree/test_align_checkpoint.py -v
```

Expected: FAIL with `ImportError: cannot import name 'verify_align_outputs'`

- [ ] **Step 3: Implement the align output verifier**

Append the following function to `phyloai/pretree/align.py` (imports are already present at the top of the file):

```python
def verify_align_outputs(aa_path: Path, nt_path: Path | None) -> bool:
    """Validate a gene's AA and (optional) NT outputs for resume."""
    aa = validate_fasta_output(aa_path, require_aligned=True)
    if not aa.ok:
        return False
    if nt_path is None:
        return True
    if not nt_path.exists() or nt_path.stat().st_size == 0:
        return False
    nt = validate_fasta_output(nt_path, require_aligned=True)
    return nt.ok
```

- [ ] **Step 4: Run tests to verify they pass**

```bash
pytest tests/pretree/test_align_checkpoint.py -v
```

Expected: all PASS

- [ ] **Step 5: Commit**

```bash
git add phyloai/pretree/align.py tests/pretree/test_align_checkpoint.py
git commit -m "feat(align): add output verifier for resume verification"
```

---

## Task 4: Pretree Align Resolved Parameters Helper

**Files:**
- Modify: `phyloai/pretree/align.py`
- Create: `tests/pretree/test_align_checkpoint.py`

- [ ] **Step 1: Write failing tests for `_resolved_align_params`**

Append to `tests/pretree/test_align_checkpoint.py`:

```python
def test_resolved_align_params_includes_required_keys() -> None:
    from phyloai.pretree.align import _resolved_align_params

    params = _resolved_align_params(
        seq_dir=Path("raw"),
        output_dir=Path("runs/run001/pretree/align"),
        method="linsi",
        resolved_seq_type="AA",
        backtrans=False,
        nt_dir=None,
        threads=8,
        extra_args=None,
        mafft_executable="/usr/bin/mafft",
        magus_executable="magus",
        trimal_executable="trimal",
        quiet=True,
    )
    expected = {
        "seq_dir", "output_dir", "method", "seq_type", "backtrans",
        "nt_dir", "threads", "extra_args", "mafft_executable",
        "magus_executable", "trimal_executable", "quiet",
    }
    assert set(params) == expected
    assert params["method"] == "linsi"
    assert params["seq_type"] == "AA"
    assert params["backtrans"] is False
    assert params["threads"] == 8
    assert params["quiet"] is True


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
        extra_args=None,
        mafft_executable="mafft",
        magus_executable="magus",
        trimal_executable="trimal",
        quiet=False,
    )
    assert "overwrite" not in params
    assert "resume" not in params
    assert "dry_run" not in params
    assert "progress_callback" not in params
```

- [ ] **Step 2: Run tests to verify they fail**

```bash
pytest tests/pretree/test_align_checkpoint.py::test_resolved_align_params_includes_required_keys \
       tests/pretree/test_align_checkpoint.py::test_resolved_align_params_excludes_mode_flags -v
```

Expected: FAIL with `ImportError: cannot import name '_resolved_align_params'`

- [ ] **Step 3: Implement `_resolved_align_params`**

Append to `phyloai/pretree/align.py`:

```python
def _resolved_align_params(
    *,
    seq_dir: Path,
    output_dir: Path,
    method: str,
    resolved_seq_type: str,
    backtrans: bool,
    nt_dir: Path | None,
    threads: int,
    extra_args: str | None,
    mafft_executable: str,
    magus_executable: str,
    trimal_executable: str,
    quiet: bool,
) -> dict[str, Any]:
    """Build the resolved parameter dictionary that determines resume eligibility."""
    return {
        "seq_dir": str(seq_dir),
        "output_dir": str(output_dir),
        "method": method,
        "seq_type": resolved_seq_type,
        "backtrans": backtrans,
        "nt_dir": str(nt_dir) if nt_dir is not None else None,
        "threads": int(threads),
        "extra_args": extra_args,
        "mafft_executable": mafft_executable,
        "magus_executable": magus_executable,
        "trimal_executable": trimal_executable,
        "quiet": bool(quiet),
    }
```

This helper must receive already-resolved values, not raw CLI inputs. In particular:

- `resolved_seq_type` must be the post-auto-detection value (`AA` or `NT`), never the literal string `"auto"`
- `mafft_executable`, `magus_executable`, and `trimal_executable` must be the actual executable strings returned by `_resolve_tool_paths`
- mode flags such as `resume`, `overwrite`, and `dry_run` must not be included in the hashable params dictionary

- [ ] **Step 4: Run tests to verify they pass**

```bash
pytest tests/pretree/test_align_checkpoint.py -v
```

Expected: all PASS

- [ ] **Step 5: Commit**

```bash
git add phyloai/pretree/align.py tests/pretree/test_align_checkpoint.py
git commit -m "feat(align): add resolved parameter helper for resume hashing"
```

---

## Task 5: Pretree Align Checkpoint Creation and Per-Task Updates

**Files:**
- Modify: `phyloai/pretree/align.py`
- Modify: `tests/pretree/test_align_checkpoint.py`

- [ ] **Step 1: Write failing tests for checkpoint creation and update**

Append to `tests/pretree/test_align_checkpoint.py`:

```python
def test_build_initial_checkpoint(tmp_path: Path) -> None:
    from phyloai.pretree.checkpoint_helpers import build_initial_checkpoint
    from phyloai.pretree.align import _resolved_align_params, verify_align_outputs

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
        extra_args=None,
        mafft_executable="mafft",
        magus_executable="magus",
        trimal_executable="trimal",
        quiet=False,
    )

    cp = build_initial_checkpoint(
        step="pretree.align",
        command="phyloai pretree align",
        params=params,
        inputs=inputs,
        output_for=lambda p: tmp_path / "out" / f"{p.stem}.fa",
        nt_output_for=lambda p: None,
    )

    assert cp.step == "pretree.align"
    assert cp.status == "running"
    assert len(cp.tasks) == 2
    assert {t.task_id for t in cp.tasks} == {"gene1", "gene2"}
    assert all(t.status == "pending" for t in cp.tasks)
    assert all(t.outputs["nt"] is None for t in cp.tasks)


def test_mark_task_updates_status(tmp_path: Path) -> None:
    from phyloai.core.checkpoint import Checkpoint, CheckpointTask
    from phyloai.pretree.checkpoint_helpers import mark_task

    cp = Checkpoint(
        schema_version=1,
        step="pretree.align",
        command="phyloai pretree align",
        status="running",
        params_hash="sha256:abc",
        params={},
        started_at="2026-06-12T10:00:00",
        updated_at="2026-06-12T10:00:00",
        completed_at=None,
        tasks=[CheckpointTask(task_id="g1", status="pending", input="raw/g1.fa",
                              outputs={"aa": "out/g1.fa", "nt": None})],
    )

    mark_task(cp, "g1", status="success", reason=None)
    task = next(t for t in cp.tasks if t.task_id == "g1")
    assert task.status == "success"
    assert task.attempts == 1
    assert task.updated_at is not None
    assert cp.updated_at != "2026-06-12T10:00:00"


def test_resume_verifier_uses_align_outputs(tmp_path: Path) -> None:
    from phyloai.core.checkpoint import Checkpoint, CheckpointTask
    from phyloai.pretree.align import verify_align_outputs
    from phyloai.pretree.checkpoint_helpers import resume_verifier

    aa = tmp_path / "g1.fa"
    aa.write_text(">a\nMKT\n>b\nMKA\n")
    cp = Checkpoint(
        schema_version=1,
        step="pretree.align",
        command="phyloai pretree align",
        status="running",
        params_hash="sha256:abc",
        params={},
        started_at="2026-06-12T10:00:00",
        updated_at="2026-06-12T10:00:00",
        completed_at=None,
        tasks=[CheckpointTask(task_id="g1", status="success", input="raw/g1.fa",
                              outputs={"aa": str(aa), "nt": None})],
    )

    verifier = resume_verifier(verify_align_outputs)
    assert verifier(cp.tasks[0]) is True

    aa.unlink()
    assert verifier(cp.tasks[0]) is False
```

- [ ] **Step 2: Run tests to verify they fail**

```bash
pytest tests/pretree/test_align_checkpoint.py::test_build_initial_checkpoint \
       tests/pretree/test_align_checkpoint.py::test_mark_task_updates_status \
       tests/pretree/test_align_checkpoint.py::test_resume_verifier_uses_align_outputs -v
```

Expected: FAIL with `ModuleNotFoundError: No module named 'phyloai.pretree.checkpoint_helpers'`

- [ ] **Step 3: Create align-specific checkpoint helpers**

Create `phyloai/pretree/checkpoint_helpers.py`:

```python
"""Checkpoint helpers specific to the pretree phase."""

from __future__ import annotations

import datetime as _dt
from collections.abc import Callable
from pathlib import Path
from typing import Any

from phyloai.core.checkpoint import (
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
    """Create the initial checkpoint for a fresh batch run."""
    now = _utc_now_iso()
    tasks: list[CheckpointTask] = []
    for inp in inputs:
        aa_out = output_for(inp)
        nt_out = nt_output_for(inp)
        tasks.append(
            CheckpointTask(
                task_id=inp.stem,
                status="pending",
                input=str(inp),
                outputs={"aa": str(aa_out), "nt": str(nt_out) if nt_out else None},
            )
        )
    return Checkpoint(
        schema_version=1,
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
    """Update a single task in-place and refresh checkpoint timestamps."""
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
    """Build a verifier closure that delegates to a command-specific function."""

    def _verifier(task: CheckpointTask) -> bool:
        aa = Path(task.outputs["aa"]) if task.outputs.get("aa") else None
        nt = Path(task.outputs["nt"]) if task.outputs.get("nt") else None
        if aa is None:
            return False
        return verify_outputs(aa, nt)

    return _verifier
```

- [ ] **Step 4: Run tests to verify they pass**

```bash
pytest tests/pretree/test_align_checkpoint.py -v
```

Expected: all PASS

- [ ] **Step 5: Commit**

```bash
git add phyloai/pretree/checkpoint_helpers.py tests/pretree/test_align_checkpoint.py
git commit -m "feat(pretree): add align checkpoint helpers and resume verifier"
```

---

## Task 6: Pretree Align Resume and Result Reconstruction

**Files:**
- Modify: `phyloai/pretree/align.py`
- Modify: `tests/pretree/test_align_checkpoint.py`

- [ ] **Step 1: Write failing tests for the resume helper functions**

Append to `tests/pretree/test_align_checkpoint.py`:

```python
def test_plan_resume_marks_invalid_success_for_rerun(tmp_path: Path) -> None:
    from phyloai.core.checkpoint import Checkpoint, CheckpointTask
    from phyloai.pretree.align import verify_align_outputs
    from phyloai.pretree.checkpoint_helpers import plan_resume

    cp = Checkpoint(
        schema_version=1,
        step="pretree.align",
        command="phyloai pretree align",
        status="running",
        params_hash="sha256:abc",
        params={},
        started_at="2026-06-12T10:00:00",
        updated_at="2026-06-12T10:00:00",
        completed_at=None,
        tasks=[
            CheckpointTask(task_id="g1", status="success", input="raw/g1.fa",
                           outputs={"aa": str(tmp_path / "missing.fa"), "nt": None}),
            CheckpointTask(task_id="g2", status="failed", input="raw/g2.fa",
                           outputs={"aa": str(tmp_path / "g2.fa"), "nt": None}),
            CheckpointTask(task_id="g3", status="pending", input="raw/g3.fa",
                           outputs={"aa": str(tmp_path / "g3.fa"), "nt": None}),
        ],
    )

    to_run, skipped = plan_resume(cp, verify_align_outputs)
    assert sorted(to_run) == ["g1", "g2", "g3"]
    assert skipped == []


def test_plan_resume_keeps_valid_success_skipped(tmp_path: Path) -> None:
    from phyloai.core.checkpoint import Checkpoint, CheckpointTask
    from phyloai.pretree.align import verify_align_outputs
    from phyloai.pretree.checkpoint_helpers import plan_resume

    aa = tmp_path / "g1.fa"
    aa.write_text(">a\nMKT\n>b\nMKA\n")
    cp = Checkpoint(
        schema_version=1,
        step="pretree.align",
        command="phyloai pretree align",
        status="running",
        params_hash="sha256:abc",
        params={},
        started_at="2026-06-12T10:00:00",
        updated_at="2026-06-12T10:00:00",
        completed_at=None,
        tasks=[CheckpointTask(task_id="g1", status="success", input="raw/g1.fa",
                               outputs={"aa": str(aa), "nt": None})],
    )

    to_run, skipped = plan_resume(cp, verify_align_outputs)
    assert to_run == []
    assert skipped == ["g1"]


def test_reconstruct_align_result_aggregates_states(tmp_path: Path) -> None:
    from phyloai.core.checkpoint import Checkpoint, CheckpointTask
    from phyloai.pretree.align import reconstruct_align_result

    aa = tmp_path / "g1.fa"
    aa.write_text(">a\nMKT\n>b\nMKA\n")
    nt = tmp_path / "g1.nt.fa"
    nt.write_text(">a\nATGAAGACT\n>b\nATGAAGGCT\n")
    cp = Checkpoint(
        schema_version=1,
        step="pretree.align",
        command="phyloai pretree align",
        status="success",
        params_hash="sha256:abc",
        params={},
        started_at="2026-06-12T10:00:00",
        updated_at="2026-06-12T10:30:00",
        completed_at="2026-06-12T10:30:00",
        tasks=[CheckpointTask(task_id="g1", status="success", input="raw/g1.fa",
                               outputs={"aa": str(aa), "nt": str(nt)},
                               attempts=1)],
    )

    payload = reconstruct_align_result(
        checkpoint=cp,
        params={"method": "linsi"},
        tool_versions={"mafft": "7.526"},
        wall_time=12.5,
        skipped_inputs=[{"path": "raw/bad.fa", "reason": "empty file"}],
        scan_warnings=[],
    )
    assert payload["status"] == "success"
    assert payload["key_results"]["n_aligned"] == 1
    assert payload["key_results"]["n_skipped"] == 1
    assert payload["data"]["files"][0]["alignment_length"] == 3
    assert payload["data"]["files"][0]["n_taxa"] == 2
```

- [ ] **Step 2: Run tests to verify they fail**

```bash
pytest tests/pretree/test_align_checkpoint.py::test_plan_resume_marks_invalid_success_for_rerun \
       tests/pretree/test_align_checkpoint.py::test_plan_resume_keeps_valid_success_skipped \
       tests/pretree/test_align_checkpoint.py::test_reconstruct_align_result_aggregates_states -v
```

Expected: FAIL with `ImportError: cannot import name 'plan_resume' from 'phyloai.pretree.checkpoint_helpers'`

- [ ] **Step 3: Implement `plan_resume` and `reconstruct_align_result`**

Append to `phyloai/pretree/checkpoint_helpers.py`:

```python
def plan_resume(
    checkpoint: Checkpoint,
    verify_outputs: Callable[[Path, Path | None], bool],
) -> tuple[list[str], list[str]]:
    """Decide which task IDs must be rerun vs. skipped after output verification."""
    to_run: list[str] = []
    skipped: list[str] = []
    for task in checkpoint.tasks:
        if task.status in ("pending", "running", "failed"):
            to_run.append(task.task_id)
            continue
        if task.status == "success":
            aa = Path(task.outputs["aa"]) if task.outputs.get("aa") else None
            nt = Path(task.outputs["nt"]) if task.outputs.get("nt") else None
            if aa is not None and verify_outputs(aa, nt):
                skipped.append(task.task_id)
            else:
                to_run.append(task.task_id)
            continue
        # skipped tasks: leave as skipped; resume will not rerun them
        skipped.append(task.task_id)
    return to_run, skipped
```

Append the result reconstruction helper to `phyloai/pretree/align.py`:

```python
def _read_alignment_metrics(path: Path) -> tuple[int, int]:
    from phyloai.core.sequence_output_validation import validate_fasta_output

    result = validate_fasta_output(path, require_aligned=True)
    return result.n_records, result.length


def reconstruct_align_result(
    *,
    checkpoint: Checkpoint,
    params: dict[str, Any],
    tool_versions: dict[str, str],
    wall_time: float,
    skipped_inputs: list[dict[str, str]],
    scan_warnings: list[str],
) -> dict[str, Any]:
    """Reconstruct the final ``result.json`` payload from checkpoint + outputs."""
    file_results: list[dict[str, Any]] = []
    failed: list[dict[str, str]] = []
    backtrans_count = 0

    for task in checkpoint.tasks:
        if task.status != "success":
            failed.append({
                "path": task.input,
                "reason": task.reason or task.status,
            })
            continue
        aa_path = Path(task.outputs["aa"]) if task.outputs.get("aa") else None
        nt_path = Path(task.outputs["nt"]) if task.outputs.get("nt") else None
        if aa_path is None or not aa_path.exists():
            failed.append({"path": task.input, "reason": "missing AA output"})
            continue
        n_taxa, alignment_length = _read_alignment_metrics(aa_path)
        if nt_path is not None and nt_path.exists():
            backtrans_count += 1
        file_results.append({
            "input": task.input,
            "output_aa": str(aa_path),
            "output_nt": str(nt_path) if nt_path else None,
            "n_taxa": n_taxa,
            "alignment_length": alignment_length,
            "wall_time": 0.0,
            "warnings": [],
        })

    aligned_lengths = [r["alignment_length"] for r in file_results if r["alignment_length"]]
    aligned_taxa = [r["n_taxa"] for r in file_results if r["n_taxa"]]
    mean_len = round(sum(aligned_lengths) / len(aligned_lengths), 1) if aligned_lengths else 0.0
    mean_taxa = round(sum(aligned_taxa) / len(aligned_taxa), 1) if aligned_taxa else 0.0

    skipped = list(skipped_inputs)
    skipped.extend(failed)

    all_warnings: list[str] = list(scan_warnings)
    for fr in file_results:
        all_warnings.extend(fr.get("warnings", []))

    return {
        "status": "success" if file_results else "error",
        "command": checkpoint.command,
        "wall_time": wall_time,
        "tool_versions": tool_versions,
        "params": params,
        "key_results": {
            "n_aligned": len(file_results),
            "n_skipped": len(skipped),
            "method": params.get("method"),
            "backtrans": params.get("backtrans", False),
            "mean_alignment_length": mean_len,
            "mean_n_taxa": mean_taxa,
        },
        "error": None if file_results else "No genes were aligned.",
        "data": {
            "summary": {
                "n_input_files": len(checkpoint.tasks) + len(skipped_inputs),
                "n_aligned": len(file_results),
                "n_backtrans": backtrans_count,
                "n_skipped": len(skipped),
            },
            "files": file_results,
            "skipped": skipped,
            "warnings": all_warnings,
        },
    }
```

- [ ] **Step 4: Run tests to verify they pass**

```bash
pytest tests/pretree/test_align_checkpoint.py -v
```

Expected: all PASS

- [ ] **Step 5: Commit**

```bash
git add phyloai/pretree/align.py phyloai/pretree/checkpoint_helpers.py tests/pretree/test_align_checkpoint.py
git commit -m "feat(align): add resume planning and result reconstruction helpers"
```

---

## Task 7: Wire `--resume` Into `run_align`

**Files:**
- Modify: `phyloai/pretree/align.py`
- Create: `tests/pretree/test_run_align_resume.py`

- [ ] **Step 1: Write failing tests for `run_align` resume behavior**

Create `tests/pretree/test_run_align_resume.py`:

```python
from __future__ import annotations

import shutil
from pathlib import Path

import pytest


def test_run_align_resume_requires_checkpoint(tmp_path: Path) -> None:
    from phyloai.pretree.align import run_align

    seq_dir = tmp_path / "seqs"
    seq_dir.mkdir()
    (seq_dir / "gene1.fa").write_text(">a\nMKT\n>b\nMKA\n")
    out_dir = tmp_path / "out"
    out_dir.mkdir()

    with pytest.raises(ValueError, match="No checkpoint"):
        run_align(
            seq_dir=seq_dir,
            output_dir=out_dir,
            method="linsi",
            seq_type="AA",
            resume=True,
        )


def test_run_align_resume_rejects_overwrite(tmp_path: Path) -> None:
    from phyloai.pretree.align import run_align

    seq_dir = tmp_path / "seqs"
    seq_dir.mkdir()
    (seq_dir / "gene1.fa").write_text(">a\nMKT\n>b\nMKA\n")
    out_dir = tmp_path / "out"
    out_dir.mkdir()
    cp = out_dir / "checkpoint.json"
    cp.write_text("{}")

    with pytest.raises(ValueError, match="--overwrite and --resume"):
        run_align(
            seq_dir=seq_dir,
            output_dir=out_dir,
            method="linsi",
            seq_type="AA",
            resume=True,
            overwrite=True,
        )


def test_run_align_resume_detects_param_mismatch(tmp_path: Path) -> None:
    from phyloai.core.checkpoint import canonical_params_hash
    from phyloai.pretree.align import run_align
    from phyloai.pretree.checkpoint_helpers import build_initial_checkpoint
    from phyloai.pretree.align import _resolved_align_params
    import json

    seq_dir = tmp_path / "seqs"
    seq_dir.mkdir()
    (seq_dir / "gene1.fa").write_text(">a\nMKT\n>b\nMKA\n")
    out_dir = tmp_path / "out"
    out_dir.mkdir()
    aa_dir = out_dir / "seqs"
    aa_dir.mkdir(parents=True)

    params = _resolved_align_params(
        seq_dir=seq_dir,
        output_dir=out_dir,
        method="linsi",
        resolved_seq_type="AA",
        backtrans=False,
        nt_dir=None,
        threads=2,
        extra_args=None,
        mafft_executable="mafft",
        magus_executable="magus",
        trimal_executable="trimal",
        quiet=False,
    )
    cp = build_initial_checkpoint(
        step="pretree.align",
        command="phyloai pretree align",
        params=params,
        inputs=[seq_dir / "gene1.fa"],
        output_for=lambda p: aa_dir / f"{p.stem}.fa",
        nt_output_for=lambda p: None,
    )
    cp.params_hash = canonical_params_hash({**params, "method": "fftns1"})
    (out_dir / "checkpoint.json").write_text(json.dumps(cp.to_dict()))

    with pytest.raises(ValueError, match="Resume parameter mismatch"):
        run_align(
            seq_dir=seq_dir,
            output_dir=out_dir,
            method="linsi",
            seq_type="AA",
            threads=2,
            resume=True,
        )


def test_run_align_resume_skips_successful_tasks(tmp_path: Path) -> None:
    if not shutil.which("mafft"):
        pytest.skip("mafft not found")
    from phyloai.core.checkpoint import canonical_params_hash
    from phyloai.pretree.align import run_align, _resolved_align_params
    from phyloai.pretree.checkpoint_helpers import build_initial_checkpoint, mark_task
    import json

    seq_dir = tmp_path / "seqs"
    seq_dir.mkdir()
    (seq_dir / "gene1.fa").write_text(">a\nMKTLLLTLVVVTIVC\n>b\nMKTLLLTLAAVTIVC\n>c\nMKTLLLTLVVVTIVC\n")
    (seq_dir / "gene2.fa").write_text(">a\nGHTLLLTLVVVTIVC\n>b\nGHTLLLTLAAVTIVC\n>c\nGHTLLLTLVVVTIVC\n")
    out_dir = tmp_path / "out"
    out_dir.mkdir()
    aa_dir = out_dir / "seqs"
    aa_dir.mkdir(parents=True)

    params = _resolved_align_params(
        seq_dir=seq_dir,
        output_dir=out_dir,
        method="linsi",
        resolved_seq_type="AA",
        backtrans=False,
        nt_dir=None,
        threads=2,
        extra_args=None,
        mafft_executable="mafft",
        magus_executable="magus",
        trimal_executable="trimal",
        quiet=False,
    )
    cp = build_initial_checkpoint(
        step="pretree.align",
        command="phyloai pretree align",
        params=params,
        inputs=[seq_dir / "gene1.fa", seq_dir / "gene2.fa"],
        output_for=lambda p: aa_dir / f"{p.stem}.fa",
        nt_output_for=lambda p: None,
    )

    # pre-align gene1 to simulate a completed task
    from phyloai.pretree.align import _align_one
    first = _align_one(
        seq_dir / "gene1.fa",
        aa_dir,
        method="linsi",
        seq_type="AA",
        extra_args=None,
        dry_run=False,
    )
    assert first["status"] == "success"
    mark_task(cp, "gene1", status="success", reason=None)
    mark_task(cp, "gene2", status="failed", reason="previous run error")
    (out_dir / "checkpoint.json").write_text(json.dumps(cp.to_dict()))

    payload = run_align(
        seq_dir=seq_dir,
        output_dir=out_dir,
        method="linsi",
        seq_type="AA",
        threads=2,
        resume=True,
    )

    assert payload["status"] == "success"
    assert payload["key_results"]["n_aligned"] == 2
    assert (aa_dir / "gene1.fa").exists()
    assert (aa_dir / "gene2.fa").exists()
```

- [ ] **Step 2: Run tests to verify they fail**

```bash
pytest tests/pretree/test_run_align_resume.py -v
```

Expected: FAIL with `TypeError: run_align() got an unexpected keyword argument 'resume'`

- [ ] **Step 3: Add `resume` parameter and resume flow to `run_align`**

Modify the signature of `run_align` in `phyloai/pretree/align.py`:

```python
def run_align(
    seq_dir: Path,
    output_dir: Path,
    method: str,
    seq_type: str,
    backtrans: bool = False,
    nt_dir: Path | None = None,
    threads: int = 4,
    extra_args: str | None = None,
    mafft_path: Path | None = None,
    magus_path: Path | None = None,
    trimal_path: Path | None = None,
    overwrite: bool = False,
    resume: bool = False,
    dry_run: bool = False,
    quiet: bool = False,
    progress_callback: Callable[[Path], None] | None = None,
) -> dict[str, Any]:
```

Move the new resume logic after input scanning, auto-detection, and tool-path resolution so parameter matching uses fully resolved values.

1. Keep the early `resume`/`overwrite` mutual-exclusion check near the top.
2. Run `_scan_input(seq_dir)` before branching so both fresh runs and resume runs have the normal scan results.
3. If `seq_type == "auto"`, resolve it before computing checkpoint params.
4. Call `_resolve_tool_paths(...)` before computing checkpoint params so the stored hash uses actual executable strings.
5. Build `resolved = _resolved_align_params(...)` with `resolved_seq_type`, resolved executable strings, and `quiet`.

The shared input scan should be captured once and reused in both branches:

```python
    found, scan_skipped = _scan_input(seq_dir)
```

Insert the new resume block in `run_align` after those resolution steps:

```python
    if resume and overwrite:
        raise ValueError("--overwrite and --resume are mutually exclusive.")

    if resume:
        from phyloai.core.checkpoint import (
            load_checkpoint,
            save_checkpoint_atomic,
            validate_resume_params,
        )
        from phyloai.pretree.checkpoint_helpers import (
            mark_task,
            plan_resume,
        )

        resolved = _resolved_align_params(
            seq_dir=seq_dir,
            output_dir=output_dir,
            method=method,
            resolved_seq_type=seq_type,
            backtrans=backtrans,
            nt_dir=nt_dir,
            threads=threads,
            extra_args=extra_args,
            mafft_executable=mafft_exe,
            magus_executable=magus_exe,
            trimal_executable=trimal_exe,
            quiet=quiet,
        )
        ckpt_path = output_dir / "checkpoint.json"
        if not ckpt_path.exists():
            raise ValueError(
                f"No checkpoint found at {ckpt_path}. "
                "Use --overwrite to start fresh."
            )
        checkpoint = load_checkpoint(ckpt_path)
        validate_resume_params(checkpoint, resolved, step="pretree.align")
        if checkpoint.status == "success":
            payload = reconstruct_align_result(
                checkpoint=checkpoint,
                params=checkpoint.params,
                tool_versions=_detect_tool_versions(
                    method=method,
                    backtrans=backtrans,
                    mafft_path=mafft_path,
                    magus_path=magus_path,
                    trimal_path=trimal_path,
                ),
                wall_time=0.0,
                skipped_inputs=list(scan_skipped),
                scan_warnings=list(scan_skipped),
            )
            return payload
        to_run, skipped_ids = plan_resume(checkpoint, verify_align_outputs)
        if not to_run:
            checkpoint.status = "success"
            checkpoint.completed_at = checkpoint.touch() and checkpoint.updated_at
            save_checkpoint_atomic(checkpoint, ckpt_path)
            payload = reconstruct_align_result(
                checkpoint=checkpoint,
                params=checkpoint.params,
                tool_versions=_detect_tool_versions(
                    method=method,
                    backtrans=backtrans,
                    mafft_path=mafft_path,
                    magus_path=magus_path,
                    trimal_path=trimal_path,
                ),
                wall_time=0.0,
                skipped_inputs=list(scan_skipped),
                scan_warnings=list(scan_skipped),
            )
            return payload

        found = [Path(t.input) for t in checkpoint.tasks if t.task_id in set(to_run)]
        skipped: list[dict[str, str]] = []
    else:
        if not dry_run:
            if output_dir.exists() and any(output_dir.iterdir()):
                if not overwrite:
                    raise ValueError(
                        f"Output directory '{output_dir}' already exists and is non-empty. "
                        "Use --overwrite to replace it."
                    )
                shutil.rmtree(output_dir)
            output_dir.mkdir(parents=True, exist_ok=True)
        skipped: list[dict[str, str]] = list(scan_skipped)
        checkpoint = None
        ckpt_path = None
        to_run = None
        skipped_ids = []
```

In the fresh-run branch, once `found` is known and `checkpoint is not None`, create the initial checkpoint from the resolved params and save it before launching workers.

Before submitting work to the process pool, mark all tasks in `to_run` as `running` in the parent process and save the checkpoint once. This is required so an interrupted run preserves the `running` state required by the spec.

```python
    if checkpoint is not None and ckpt_path is not None and to_run and not dry_run:
        for task_id in to_run:
            mark_task(checkpoint, task_id, status="running", reason=None)
        save_checkpoint_atomic(checkpoint, ckpt_path)
```

After the worker pool is filled and before the backtranslation loop, when `checkpoint is not None`, add:

```python
    if checkpoint is not None and ckpt_path is not None and to_run:
        to_run_set = set(to_run)
        for res in file_results:
            stem = Path(res["input"]).stem
            if stem not in to_run_set:
                continue
            if res["status"] == "success":
                mark_task(checkpoint, stem, status="success", reason=None)
            else:
                mark_task(checkpoint, stem, status="failed", reason=res.get("reason"))
        save_checkpoint_atomic(checkpoint, ckpt_path)
```

At the end of `run_align`, before the final payload assembly, add:

```python
    if checkpoint is not None and ckpt_path is not None and not dry_run:
        checkpoint.status = "success"
        checkpoint.completed_at = checkpoint.touch() and checkpoint.updated_at
        save_checkpoint_atomic(checkpoint, ckpt_path)
        payload = reconstruct_align_result(
            checkpoint=checkpoint,
            params=resolved,
            tool_versions=_detect_tool_versions(
                method=method,
                backtrans=backtrans,
                mafft_path=mafft_path,
                magus_path=magus_path,
                trimal_path=trimal_path,
            ),
            wall_time=time.monotonic() - run_start,
            skipped_inputs=skipped,
            scan_warnings=list(scan_skipped),
        )
        return payload
```

The existing final payload block (`payload: dict[str, Any] = { ... }`) should now be the non-resume branch only.

Also add a targeted test in `tests/pretree/test_run_align_resume.py` that a checkpoint with `status == "success"` but mismatched `params_hash` still raises `ValueError("Resume parameter mismatch")` rather than returning early.

- [ ] **Step 4: Run tests to verify they pass**

```bash
pytest tests/pretree/test_run_align_resume.py -v
```

Expected: all PASS (tests requiring `mafft` skip if absent)

- [ ] **Step 5: Run the existing align test suite to ensure no regressions**

```bash
pytest tests/pretree/test_align.py -v
```

Expected: all PASS

- [ ] **Step 6: Commit**

```bash
git add phyloai/pretree/align.py tests/pretree/test_run_align_resume.py
git commit -m "feat(align): add --resume flow to run_align with checkpoint wiring"
```

---

## Task 8: CLI `--resume` Flag and Mutual Exclusion

**Files:**
- Modify: `phyloai/cli/commands/pretree.py`
- Create: `tests/cli/test_pretree_align_resume.py`

- [ ] **Step 1: Write failing CLI tests for `--resume`**

Create `tests/cli/test_pretree_align_resume.py`:

```python
from __future__ import annotations

from pathlib import Path
import json

from click.testing import CliRunner

from phyloai.cli.commands.pretree import pretree


def test_align_cli_rejects_overwrite_with_resume(tmp_path: Path) -> None:
    runner = CliRunner()
    seq_dir = tmp_path / "seqs"
    seq_dir.mkdir()
    (seq_dir / "gene1.fa").write_text(">a\nMKT\n>b\nMKA\n")
    out_dir = tmp_path / "out"
    out_dir.mkdir()
    (out_dir / "old.txt").write_text("old")

    result = runner.invoke(
        pretree,
        [
            "align",
            "--seq-dir", str(seq_dir),
            "--method", "linsi",
            "--seq-type", "AA",
            "--output-dir", str(out_dir),
            "--resume",
            "--overwrite",
        ],
    )
    assert result.exit_code != 0
    assert "--overwrite and --resume" in result.output


def test_align_cli_resume_requires_checkpoint(tmp_path: Path) -> None:
    runner = CliRunner()
    seq_dir = tmp_path / "seqs"
    seq_dir.mkdir()
    (seq_dir / "gene1.fa").write_text(">a\nMKT\n>b\nMKA\n")
    out_dir = tmp_path / "out"
    out_dir.mkdir()

    result = runner.invoke(
        pretree,
        [
            "align",
            "--seq-dir", str(seq_dir),
            "--method", "linsi",
            "--seq-type", "AA",
            "--output-dir", str(out_dir),
            "--resume",
        ],
    )
    assert result.exit_code != 0
    assert "No checkpoint" in result.output


def test_align_cli_resume_dry_run_validates_params(tmp_path: Path) -> None:
    from phyloai.pretree.align import _resolved_align_params
    from phyloai.pretree.checkpoint_helpers import build_initial_checkpoint

    runner = CliRunner()
    seq_dir = tmp_path / "seqs"
    seq_dir.mkdir()
    (seq_dir / "gene1.fa").write_text(">a\nMKT\n>b\nMKA\n")
    out_dir = tmp_path / "out"
    out_dir.mkdir()
    aa_dir = out_dir / "seqs"
    aa_dir.mkdir()

    params = _resolved_align_params(
        seq_dir=seq_dir,
        output_dir=out_dir,
        method="linsi",
        resolved_seq_type="AA",
        backtrans=False,
        nt_dir=None,
        threads=4,
        extra_args=None,
        mafft_executable="mafft",
        magus_executable="magus",
        trimal_executable="trimal",
        quiet=False,
    )
    checkpoint = build_initial_checkpoint(
        step="pretree.align",
        command="phyloai pretree align",
        params=params,
        inputs=[seq_dir / "gene1.fa"],
        output_for=lambda p: aa_dir / f"{p.stem}.fa",
        nt_output_for=lambda p: None,
    )
    (out_dir / "checkpoint.json").write_text(json.dumps(checkpoint.to_dict()))

    result = runner.invoke(
        pretree,
        [
            "align",
            "--seq-dir", str(seq_dir),
            "--method", "fftns1",
            "--seq-type", "AA",
            "--output-dir", str(out_dir),
            "--resume",
            "--dry-run",
        ],
    )
    assert result.exit_code == 1
    assert "Resume parameter mismatch" in result.output
```

- [ ] **Step 2: Run tests to verify they fail**

```bash
pytest tests/cli/test_pretree_align_resume.py -v
```

Expected: FAIL because `--resume` is not a recognized option

- [ ] **Step 3: Add `--resume` to the align CLI command**

In `phyloai/cli/commands/pretree.py`, insert the new option just before `--overwrite` in `align_command`:

```python
@click.option("--resume", is_flag=True, default=False,
              help=(
                  "Resume from checkpoint.json in the output directory. "
                  "Requires the same parameters as the original run. "
                  "Mutually exclusive with --overwrite."
              ))
```

Update the `align_command` signature to add `resume: bool`:

```python
def align_command(
    seq_dir: Path,
    method: str,
    seq_type: str,
    backtrans: bool,
    nt_dir: Path | None,
    output_dir: Path,
    threads: int,
    extra_args: str | None,
    mafft_path: Path | None,
    magus_path: Path | None,
    trimal_path: Path | None,
    resume: bool,
    overwrite: bool,
    dry_run: bool,
    quiet: bool,
) -> None:
```

Pass both `resume=resume` and `quiet=quiet` into the `run_align` call inside `_invoke`. If `resume` is set, skip the dry-run early return path so resume dry-run summaries still print.

In the dry-run block, after the existing `dry_run` print, add a resume summary section that performs the same validation ordering as `run_align`: load checkpoint, validate step and params, then summarize skip/rerun/invalid counts.

```python
    if resume and dry_run:
        from phyloai.core.checkpoint import (
            load_checkpoint,
            summarize_resume_tasks,
            validate_resume_params,
        )
        from phyloai.pretree.align import (
            _detect_seq_type_from_files,
            _resolved_align_params,
            _resolve_tool_paths,
            _scan_input,
            verify_align_outputs,
        )
        from phyloai.pretree.checkpoint_helpers import resume_verifier

        ckpt_path = output_dir / "checkpoint.json"
        found, _ = _scan_input(seq_dir)
        resolved_seq_type = seq_type
        if resolved_seq_type == "auto":
            resolved_seq_type = _detect_seq_type_from_files(found) if found else "AA"
        mafft_exe, magus_exe, trimal_exe = _resolve_tool_paths(
            method=method,
            backtrans=backtrans,
            mafft_path=mafft_path,
            magus_path=magus_path,
            trimal_path=trimal_path,
            dry_run=True,
        )
        resolved = _resolved_align_params(
            seq_dir=seq_dir,
            output_dir=output_dir,
            method=method,
            resolved_seq_type=resolved_seq_type,
            backtrans=backtrans,
            nt_dir=nt_dir,
            threads=threads,
            extra_args=extra_args,
            mafft_executable=mafft_exe,
            magus_executable=magus_exe,
            trimal_executable=trimal_exe,
            quiet=quiet,
        )
        checkpoint = load_checkpoint(ckpt_path)
        validate_resume_params(checkpoint, resolved, step="pretree.align")
        summary = summarize_resume_tasks(
            checkpoint,
            resume_verifier(verify_align_outputs),
        )
        if checkpoint.status == "success":
            click.echo("Resume dry-run: run already complete; 0 tasks rerun.", err=True)
            return
        click.echo(
            f"Resume dry-run: skip {summary['skip']} tasks, "
            f"rerun {summary['rerun']} tasks, "
            f"invalidate {summary['invalid']} recorded successes.",
            err=True,
        )
        return
```

- [ ] **Step 4: Run tests to verify they pass**

```bash
pytest tests/cli/test_pretree_align_resume.py -v
```

Expected: all PASS

- [ ] **Step 5: Run the full existing align CLI test suite**

```bash
pytest tests/cli -v
```

Expected: all PASS

- [ ] **Step 6: Commit**

```bash
git add phyloai/cli/commands/pretree.py tests/cli/test_pretree_align_resume.py
git commit -m "feat(cli): expose --resume on pretree align command"
```

---

## Task 9: Documentation Update for `--resume`

**Files:**
- Modify: `docs/commands/pretree-align.md` (create if missing)
- Modify: `README.md` (only if command index needs an entry)

- [ ] **Step 1: Confirm command doc exists; if not, scaffold it**

```bash
test -f docs/commands/pretree-align.md && echo "exists" || echo "missing"
```

If the file is missing, create it with the standard sections described in the parent design Section 4.5 (Purpose, Usage, Inputs, Outputs, Examples, Warnings and Errors, Notes).

- [ ] **Step 2: Add the Resume behavior section**

Append a `## Resume behavior` section to `docs/commands/pretree-align.md`:

```markdown
## Resume behavior

`pretree align` supports `--resume` to recover from interruption, power loss, or
external tool failure without redoing completed work.

- A run interrupted after some alignments completed can be continued with
  `phyloai pretree align ... --resume`.
- The output directory must already contain `checkpoint.json`. The current
  invocation's resolved parameters (analysis parameters and run-control
  parameters such as `--threads`) must match the checkpoint exactly. A mismatch
  exits with code 1.
- Tasks with status `success` and valid output files are skipped. Tasks with
  status `failed`, `pending`, `running`, or `success` whose outputs are
  missing or invalid are rerun.
- `--resume` and `--overwrite` are mutually exclusive. Use `--overwrite` to
  start a fresh run.
- Resume appends to `align.log` and rewrites `result.json` on completion.
```

Add a small example to the Examples section:

```markdown
# Resume an interrupted run
phyloai pretree align --seq-dir ./raw_aa --method linsi --seq-type AA \
  --output-dir ./runs/run001/pretree/align --resume
```

- [ ] **Step 3: Update README command index entry if it lists `pretree align`**

Open `README.md`. Locate the align entry in the command index. If it lists
`pretree align` as a link to `docs/commands/pretree-align.md`, no change is
required. If the entry does not mention resume, append a one-line note:
"`--resume` supported for long runs."

- [ ] **Step 4: Commit**

```bash
git add docs/commands/pretree-align.md README.md
git commit -m "docs(align): document --resume and checkpoint behavior"
```

---

## Task 10: Final Smoke Test and Verification

**Files:**
- Read-only verification

- [ ] **Step 1: Run the full test suite**

```bash
pytest -v
```

Expected: all PASS; tool-dependent tests skip cleanly if `mafft` or `trimal`
are absent.

- [ ] **Step 2: Run a manual resume scenario**

Create a minimal input directory of two short genes, run `pretree align`,
simulate an interruption by deleting one of the output FASTA files, then run
again with `--resume`:

```bash
mkdir -p /tmp/phyloai-resume/raw
printf '>a\nMKTLLLTLVVVTIVC\n>b\nMKTLLLTLAAVTIVC\n' > /tmp/phyloai-resume/raw/g1.fa
printf '>a\nGHTLLLTLVVVTIVC\n>b\nGHTLLLTLAAVTIVC\n' > /tmp/phyloai-resume/raw/g2.fa
rm -rf /tmp/phyloai-resume/out
phyloai pretree align --seq-dir /tmp/phyloai-resume/raw --method linsi \
  --seq-type AA --output-dir /tmp/phyloai-resume/out
rm /tmp/phyloai-resume/out/seqs/g2.fa
phyloai pretree align --seq-dir /tmp/phyloai-resume/raw --method linsi \
  --seq-type AA --output-dir /tmp/phyloai-resume/out --resume
```

Expected:
- After the second run, both `g1.fa` and `g2.fa` exist under `seqs/`.
- `result.json` lists `n_aligned == 2`.
- `checkpoint.json` has `status: success`.

- [ ] **Step 3: Run a parameter-mismatch smoke test**

```bash
phyloai pretree align --seq-dir /tmp/phyloai-resume/raw --method fftns1 \
  --seq-type AA --output-dir /tmp/phyloai-resume/out --resume
```

Expected: exit code 1 with `Resume parameter mismatch` in stderr.

- [ ] **Step 4: Commit the verification artifacts (if any new ones were produced)**

This task usually does not produce new tracked files. If a scratch file was
created outside `/tmp`, remove it and confirm `git status` shows no unintended
verification artifacts. Do not require a globally clean working tree; this
repository may already contain unrelated or intentionally modified files.

```bash
git status
```

Expected: no unexpected scratch artifacts from Task 10 remain in the worktree.

---

## Spec Coverage Notes

The plan implements every requirement in `docs/superpowers/specs/2026-06-12-checkpoint-resume-design.md`:

- Spec §3 CLI semantics: covered by Task 8 and Task 7.
- Spec §4 checkpoint file model: covered by Task 1 and Task 5.
- Spec §5 atomic writes: covered by `save_checkpoint_atomic` in Task 1 and used in Task 6 and Task 7.
- Spec §6 parameter matching: covered by `canonical_params_hash` and `validate_resume_params` (Tasks 1 and 7).
- Spec §7 `pretree align` resume flow: covered by Tasks 3-7.
- Spec §8 future command mapping: documented in the new spec; not implemented in this plan by design.
- Spec §9 core API design: Task 1.
- Spec §10 error handling: Task 7 raises user-input errors with exit 1 paths and Task 8 surfaces them via the CLI.
- Spec §11 documentation: Task 9.
- Spec §12 testing: Tasks 1-8 include unit and CLI tests; Task 10 performs a manual smoke test.
