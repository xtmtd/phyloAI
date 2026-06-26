# `phyloai run` Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Implement `phyloai run` — a one-click pipeline command that orchestrates convert → align → trim → [filter taper] → [concat | gene trees] → [iqtree | wastral] for supermatrix and supertree modes, with `--speed normal|fast`, `--resume`, Rich progress display, and a `run_checkpoint.json` step-level checkpoint.

**Architecture:** `phyloai run` lives in `phyloai/cli/commands/run.py` (CLI) and calls library functions directly (no subprocess). The orchestrator manages a `run_checkpoint.json` at the run output root using the existing `phyloai/core/checkpoint.py` helpers. Each step's detailed checkpoint (align, trim, gene trees, iqtree) remains in that step's own subdirectory and is managed by the step's own library function.

**Tech Stack:** Python 3.11+, Click, Rich, existing `phyloai` library layer (`convert_input`, `run_align`, `run_trim`, `run_taper`, `run_concat`, `run_fasttree`, `run_iqtree`, `run_wastral`), `phyloai/core/checkpoint.py`.

## Global Constraints

- All parameter names use kebab-case in CLI; snake_case in Python.
- `--output-dir` default: `./runs/run` (fixed, no timestamp).
- `--resume` and `--overwrite` are mutually exclusive; both together → exit 1.
- Internal calls use the Python library layer, not subprocess.
- Step subdirectories use numeric prefixes: `1-convert`, `2-align`, `3-trim`, `4-filter`, `5-concat` or `5-genetrees`, `6-tree`.
- `4-filter/` directory is NOT created when `--speed fast` is used.
- `result.json` is written at `<output-dir>/result.json` only on successful completion (or on failure with `status: "error"`).
- `run_checkpoint.json` uses `phyloai/core/checkpoint.py` `save_checkpoint_atomic` with `fsync=True` on terminal writes only.
- Rich progress: each step prints a `[N/total]` header; batch steps show a progress bar.
- Exit codes: 0 success, 1 user/input error, 2 tool failure, 3 environment error.
- All `result.json` fields follow `2026-06-21-phyloai-json-output-standard.md`.
- Tests go in `tests/cli/test_run.py` using Click's `CliRunner`.

---

## File Map

| File | Action | Responsibility |
|------|--------|---------------|
| `phyloai/cli/commands/run.py` | **Create** | CLI command definition, orchestrator logic, Rich progress, `run_checkpoint.json` management |
| `phyloai/cli/main.py` | **Modify** | Register `run` command in the root `cli` group |
| `tests/cli/test_run.py` | **Create** | CLI tests for `run` command |
| `docs/commands/run.md` | **Create** | User-facing command reference |

No new library modules needed — all step logic is already in the existing library layer.

---

## Task 1: Scaffolding — `run.py` and CLI registration

**Files:**
- Create: `phyloai/cli/commands/run.py`
- Modify: `phyloai/cli/main.py`
- Test: `tests/cli/test_run.py`

**Interfaces:**
- Produces: `run` Click command importable from `phyloai.cli.commands.run`

- [ ] **Step 1: Write the failing test**

```python
# tests/cli/test_run.py
from __future__ import annotations

from click.testing import CliRunner
from phyloai.cli.main import cli


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
```

- [ ] **Step 2: Run test to verify it fails**

```bash
cd /Users/zf/data/coding/phyloAI
uv run pytest tests/cli/test_run.py::test_run_help -v
```

Expected: FAIL — `run` command not found.

- [ ] **Step 3: Create `phyloai/cli/commands/run.py` with the Click command skeleton**

```python
"""One-click phylogenomics pipeline."""

from __future__ import annotations

from pathlib import Path

import click
from rich.console import Console

console = Console()


@click.command(
    "run",
    help=(
        "One-click phylogenomics pipeline from raw sequences to a species tree.\n\n"
        "Runs the full pipeline automatically using sensible defaults. "
        "For fine-grained control over individual steps, use the constituent "
        "subcommands (phyloai pretree align, phyloai tree ml iqtree, etc.).\n\n"
        "Modes:\n\n"
        "  supermatrix  convert → align → trim → [filter] → concat → iqtree\n\n"
        "  supertree    convert → align → trim → [filter] → gene trees → wastral\n\n"
        "The [filter] step (TAPER error-site masking) is included in --speed normal "
        "and skipped in --speed fast.\n\n"
        "Speed modes:\n\n"
        "  normal  MAFFT linsi, trimAl -automated1, TAPER filter, IQ-TREE3 / FastTree\n\n"
        "  fast    MAFFT auto, trimAl -automated1, no filter, FastTree\n\n"
        "Examples:\n\n"
        "  phyloai run --seq-dir ./markers --mode supermatrix\n\n"
        "  phyloai run --seq-dir ./markers --mode supertree --speed fast --threads 16\n\n"
        "  phyloai run --seq-dir ./markers --mode supermatrix --resume"
    ),
)
@click.option(
    "--seq-dir",
    required=True,
    type=click.Path(exists=True, file_okay=False, path_type=Path),
    help="Input sequence directory (any format; always converted first).",
)
@click.option(
    "--mode",
    type=click.Choice(["supermatrix", "supertree"]),
    default="supermatrix",
    show_default=True,
    help="Pipeline mode: supermatrix (concat → iqtree) or supertree (gene trees → wastral).",
)
@click.option(
    "--speed",
    type=click.Choice(["normal", "fast"]),
    default="normal",
    show_default=True,
    help=(
        "Speed/accuracy trade-off. normal: MAFFT linsi + TAPER + IQ-TREE3. "
        "fast: MAFFT auto, no TAPER, FastTree."
    ),
)
@click.option(
    "--output-dir",
    "-o",
    type=click.Path(path_type=Path),
    default=Path("runs/run"),
    show_default=True,
    help="Root output directory for all pipeline steps.",
)
@click.option(
    "--threads",
    "-t",
    type=int,
    default=4,
    show_default=True,
    help="Thread count passed to all steps.",
)
@click.option("--resume", is_flag=True, default=False, help="Resume from run_checkpoint.json.")
@click.option("--overwrite", is_flag=True, default=False, help="Delete and recreate output directory.")
@click.option("--dry-run", is_flag=True, default=False, help="Show steps without running.")
@click.option("--quiet", "-q", is_flag=True, default=False, help="Suppress non-error output.")
def run(
    seq_dir: Path,
    mode: str,
    speed: str,
    output_dir: Path,
    threads: int,
    resume: bool,
    overwrite: bool,
    dry_run: bool,
    quiet: bool,
) -> None:
    from phyloai.cli.commands._run_pipeline import execute_pipeline
    execute_pipeline(
        seq_dir=seq_dir,
        mode=mode,
        speed=speed,
        output_dir=output_dir,
        threads=threads,
        resume=resume,
        overwrite=overwrite,
        dry_run=dry_run,
        quiet=quiet,
    )
```

- [ ] **Step 4: Register `run` in `phyloai/cli/main.py`**

In `phyloai/cli/main.py`, add the import and `cli.add_command(run)`:

```python
# Add after the existing imports:
from phyloai.cli.commands.run import run

# Add in _RootGroup.list_commands — insert "run" before "pretree":
def list_commands(self, ctx: click.Context) -> list[str]:
    return ["completion", "doctor", "run", "pretree", "tree", "posttree"]

# Add after cli.add_command(posttree):
cli.add_command(run)
```

- [ ] **Step 5: Create stub `phyloai/cli/commands/_run_pipeline.py`**

```python
"""Pipeline orchestration for phyloai run."""

from __future__ import annotations

from pathlib import Path
import sys


def execute_pipeline(
    *,
    seq_dir: Path,
    mode: str,
    speed: str,
    output_dir: Path,
    threads: int,
    resume: bool,
    overwrite: bool,
    dry_run: bool,
    quiet: bool,
) -> None:
    raise NotImplementedError("Pipeline not yet implemented")
```

- [ ] **Step 6: Run test to verify it passes**

```bash
uv run pytest tests/cli/test_run.py::test_run_help -v
```

Expected: PASS — `--help` shows all options including `--speed`, `--mode`, `--resume`.

- [ ] **Step 7: Commit**

```bash
git add phyloai/cli/commands/run.py phyloai/cli/commands/_run_pipeline.py phyloai/cli/main.py tests/cli/test_run.py
git commit -m "feat(run): scaffold phyloai run CLI command with --help"
```

---

## Task 2: Run-level checkpoint helpers

**Files:**
- Modify: `phyloai/cli/commands/_run_pipeline.py`
- Test: `tests/cli/test_run.py`

**Interfaces:**
- Consumes: `phyloai.core.checkpoint.canonical_params_hash`, `save_checkpoint_atomic`, `load_checkpoint`, `CHECKPOINT_SCHEMA_VERSION`
- Produces:
  - `_build_run_params(seq_dir, mode, speed, threads, output_dir) -> dict`
  - `_build_run_checkpoint(command_str, params) -> dict`  returns a plain dict (not a `Checkpoint` dataclass) matching the `run_checkpoint.json` schema
  - `_load_run_checkpoint(checkpoint_path) -> dict`
  - `_validate_run_resume(checkpoint, params_hash) -> None` raises `click.ClickException` on mismatch

- [ ] **Step 1: Write the failing tests**

```python
# tests/cli/test_run.py  (append)
import json
from pathlib import Path
from phyloai.cli.commands._run_pipeline import (
    _build_run_params,
    _build_run_checkpoint,
    _load_run_checkpoint,
    _validate_run_resume,
)
from phyloai.core.checkpoint import canonical_params_hash
import click
import pytest


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
    # Change the hash to simulate mismatch
    ckpt["params_hash"] = "sha256:deadbeef"
    with pytest.raises(click.ClickException, match="Parameter mismatch"):
        _validate_run_resume(ckpt, canonical_params_hash(params))


def test_load_run_checkpoint_missing_raises(tmp_path: Path) -> None:
    with pytest.raises(click.ClickException, match="run_checkpoint.json"):
        _load_run_checkpoint(tmp_path / "run_checkpoint.json")
```

- [ ] **Step 2: Run tests to verify they fail**

```bash
uv run pytest tests/cli/test_run.py -k "checkpoint or params or resume_mismatch or load_run" -v
```

Expected: FAIL — functions not defined.

- [ ] **Step 3: Implement checkpoint helpers in `_run_pipeline.py`**

```python
"""Pipeline orchestration for phyloai run."""

from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Any

import click

from phyloai.core.checkpoint import (
    CHECKPOINT_SCHEMA_VERSION,
    canonical_params_hash,
    save_checkpoint_atomic,
)


# Steps per mode × speed. "skipped" steps are included so the list is always
# canonical; the orchestrator filters them when computing the step counter.
_STEP_DEFINITIONS: dict[tuple[str, str], list[dict[str, Any]]] = {
    ("supermatrix", "normal"): [
        {"name": "convert",      "subdir": "1-convert"},
        {"name": "align",        "subdir": "2-align"},
        {"name": "trim",         "subdir": "3-trim"},
        {"name": "filter_taper", "subdir": "4-filter"},
        {"name": "concat",       "subdir": "5-concat"},
        {"name": "tree",         "subdir": "6-tree"},
    ],
    ("supermatrix", "fast"): [
        {"name": "convert",      "subdir": "1-convert"},
        {"name": "align",        "subdir": "2-align"},
        {"name": "trim",         "subdir": "3-trim"},
        {"name": "concat",       "subdir": "5-concat"},
        {"name": "tree",         "subdir": "6-tree"},
    ],
    ("supertree", "normal"): [
        {"name": "convert",      "subdir": "1-convert"},
        {"name": "align",        "subdir": "2-align"},
        {"name": "trim",         "subdir": "3-trim"},
        {"name": "filter_taper", "subdir": "4-filter"},
        {"name": "genetrees",    "subdir": "5-genetrees"},
        {"name": "tree",         "subdir": "6-tree"},
    ],
    ("supertree", "fast"): [
        {"name": "convert",      "subdir": "1-convert"},
        {"name": "align",        "subdir": "2-align"},
        {"name": "trim",         "subdir": "3-trim"},
        {"name": "genetrees",    "subdir": "5-genetrees"},
        {"name": "tree",         "subdir": "6-tree"},
    ],
}


def _build_run_params(
    seq_dir: Path,
    mode: str,
    speed: str,
    threads: int,
    output_dir: Path,
) -> dict[str, Any]:
    return {
        "seq_dir": str(seq_dir.resolve()),
        "mode": mode,
        "speed": speed,
        "threads": threads,
        "output_dir": str(output_dir.resolve()),
    }


def _build_run_checkpoint(
    command_str: str,
    params: dict[str, Any],
    *,
    mode: str,
    speed: str,
) -> dict[str, Any]:
    import datetime as _dt

    now = _dt.datetime.now(_dt.timezone.utc).isoformat(timespec="seconds")
    step_defs = _STEP_DEFINITIONS[(mode, speed)]
    steps = [
        {
            "name": defn["name"],
            "status": "pending",
            "output_dir": None,
        }
        for defn in step_defs
    ]
    return {
        "schema_version": CHECKPOINT_SCHEMA_VERSION,
        "step": "run",
        "command": command_str,
        "status": "running",
        "params_hash": canonical_params_hash(params),
        "params": params,
        "started_at": now,
        "updated_at": now,
        "completed_at": None,
        "steps": steps,
    }


def _load_run_checkpoint(checkpoint_path: Path) -> dict[str, Any]:
    if not checkpoint_path.exists():
        raise click.ClickException(
            f"run_checkpoint.json not found at {checkpoint_path}. "
            "Use --overwrite to start a clean run."
        )
    try:
        data = json.loads(checkpoint_path.read_text(encoding="utf-8"))
    except Exception as exc:
        raise click.ClickException(f"Malformed run_checkpoint.json: {exc}") from exc
    if data.get("schema_version") != CHECKPOINT_SCHEMA_VERSION:
        raise click.ClickException(
            f"Unsupported run_checkpoint.json schema version: {data.get('schema_version')}. "
            "Use --overwrite to start a clean run."
        )
    if data.get("step") != "run":
        raise click.ClickException(
            "run_checkpoint.json belongs to a different command. "
            "Use --overwrite to start a clean run."
        )
    return data


def _validate_run_resume(checkpoint: dict[str, Any], current_hash: str) -> None:
    if checkpoint["params_hash"] != current_hash:
        raise click.ClickException(
            "Parameter mismatch: current parameters differ from the original run. "
            "Use --overwrite to start a clean run with the new parameters."
        )


def _save_run_checkpoint(checkpoint: dict[str, Any], path: Path, *, fsync: bool = False) -> None:
    import datetime as _dt
    checkpoint["updated_at"] = _dt.datetime.now(_dt.timezone.utc).isoformat(timespec="seconds")
    save_checkpoint_atomic(path, checkpoint, fsync=fsync)


def execute_pipeline(
    *,
    seq_dir: Path,
    mode: str,
    speed: str,
    output_dir: Path,
    threads: int,
    resume: bool,
    overwrite: bool,
    dry_run: bool,
    quiet: bool,
) -> None:
    raise NotImplementedError("Pipeline not yet implemented")
```

- [ ] **Step 4: Run tests to verify they pass**

```bash
uv run pytest tests/cli/test_run.py -k "checkpoint or params or resume_mismatch or load_run" -v
```

Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add phyloai/cli/commands/_run_pipeline.py tests/cli/test_run.py
git commit -m "feat(run): add run-level checkpoint helpers and step definitions"
```

---

## Task 3: Output directory setup and `--resume`/`--overwrite` guard

**Files:**
- Modify: `phyloai/cli/commands/_run_pipeline.py`
- Test: `tests/cli/test_run.py`

**Interfaces:**
- Consumes: `_build_run_params`, `_build_run_checkpoint`, `_load_run_checkpoint`, `_validate_run_resume`, `_save_run_checkpoint`
- Produces: `execute_pipeline` handles directory setup and resume/overwrite conflict before dispatching steps

- [ ] **Step 1: Write the failing tests**

```python
# tests/cli/test_run.py  (append)
import shutil


def _make_seq_dir(tmp_path: Path) -> Path:
    d = tmp_path / "markers"
    d.mkdir()
    (d / "gene1.fa").write_text(">sp1\nMKT\n>sp2\nMKA\n")
    (d / "gene2.fa").write_text(">sp1\nGHT\n>sp2\nGHA\n")
    return d


def test_run_resume_and_overwrite_mutually_exclusive(tmp_path: Path) -> None:
    runner = CliRunner()
    seq_dir = _make_seq_dir(tmp_path)
    result = runner.invoke(cli, [
        "run", "--seq-dir", str(seq_dir),
        "--resume", "--overwrite",
    ])
    assert result.exit_code == 1
    assert "mutually exclusive" in result.output.lower() or "mutually exclusive" in (result.exception or "")


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
    assert "run_checkpoint.json" in result.output


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
    assert "non-empty" in result.output or "--overwrite" in result.output
```

- [ ] **Step 2: Run tests to verify they fail**

```bash
uv run pytest tests/cli/test_run.py -k "mutually_exclusive or without_checkpoint or nonempty_output" -v
```

Expected: FAIL.

- [ ] **Step 3: Implement directory setup and guards in `execute_pipeline`**

Replace the `execute_pipeline` stub in `_run_pipeline.py`:

```python
def execute_pipeline(
    *,
    seq_dir: Path,
    mode: str,
    speed: str,
    output_dir: Path,
    threads: int,
    resume: bool,
    overwrite: bool,
    dry_run: bool,
    quiet: bool,
) -> None:
    import shutil
    from rich.console import Console
    console = Console(quiet=quiet)

    # --- Guard: mutually exclusive flags ---
    if resume and overwrite:
        raise click.ClickException("--resume and --overwrite are mutually exclusive.")

    checkpoint_path = output_dir / "run_checkpoint.json"

    # --- Build resolved params (for hash) ---
    params = _build_run_params(seq_dir, mode, speed, threads, output_dir)
    params_hash = canonical_params_hash(params)

    # --- Build command string for result.json ---
    import sys
    command_str = " ".join(sys.argv)

    # --- Load or initialise checkpoint ---
    if resume:
        checkpoint = _load_run_checkpoint(checkpoint_path)
        _validate_run_resume(checkpoint, params_hash)
    else:
        # Default: fail if non-empty
        if output_dir.exists() and any(output_dir.iterdir()):
            if overwrite:
                if not dry_run:
                    shutil.rmtree(output_dir)
            else:
                raise click.ClickException(
                    f"Output directory '{output_dir}' is non-empty. "
                    "Use --overwrite to replace it or --resume to continue."
                )
        if not dry_run:
            output_dir.mkdir(parents=True, exist_ok=True)
        checkpoint = _build_run_checkpoint(command_str, params, mode=mode, speed=speed)
        step_defs = _STEP_DEFINITIONS[(mode, speed)]
        for i, defn in enumerate(step_defs):
            checkpoint["steps"][i]["output_dir"] = str(output_dir / defn["subdir"])
        if not dry_run:
            _save_run_checkpoint(checkpoint, checkpoint_path)

    # --- Placeholder: dispatch steps (Task 4+) ---
    raise NotImplementedError("Step dispatch not yet implemented")
```

- [ ] **Step 4: Run tests to verify they pass**

```bash
uv run pytest tests/cli/test_run.py -k "mutually_exclusive or without_checkpoint or nonempty_output" -v
```

Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add phyloai/cli/commands/_run_pipeline.py tests/cli/test_run.py
git commit -m "feat(run): add output directory setup, --resume/--overwrite guards"
```

---

## Task 4: Dry-run mode and step listing

**Files:**
- Modify: `phyloai/cli/commands/_run_pipeline.py`
- Test: `tests/cli/test_run.py`

**Interfaces:**
- Consumes: `_STEP_DEFINITIONS`, `execute_pipeline`
- Produces: `--dry-run` prints ordered step list with tool names and exits 0

- [ ] **Step 1: Write the failing tests**

```python
# tests/cli/test_run.py  (append)
def test_run_dry_run_supermatrix_normal(tmp_path: Path) -> None:
    runner = CliRunner()
    seq_dir = _make_seq_dir(tmp_path)
    result = runner.invoke(cli, [
        "run", "--seq-dir", str(seq_dir),
        "--mode", "supermatrix", "--speed", "normal",
        "--dry-run",
    ])
    assert result.exit_code == 0
    assert "convert" in result.output
    assert "align" in result.output
    assert "linsi" in result.output
    assert "filter" in result.output   # TAPER present in normal
    assert "concat" in result.output
    assert "iqtree" in result.output.lower()


def test_run_dry_run_supermatrix_fast_no_filter(tmp_path: Path) -> None:
    runner = CliRunner()
    seq_dir = _make_seq_dir(tmp_path)
    result = runner.invoke(cli, [
        "run", "--seq-dir", str(seq_dir),
        "--mode", "supermatrix", "--speed", "fast",
        "--dry-run",
    ])
    assert result.exit_code == 0
    assert "filter" not in result.output.lower() or "skipped" in result.output.lower()
    assert "fasttree" in result.output.lower()
    # Step counter shows 5 not 6
    assert "[1/5]" in result.output


def test_run_dry_run_supertree_normal(tmp_path: Path) -> None:
    runner = CliRunner()
    seq_dir = _make_seq_dir(tmp_path)
    result = runner.invoke(cli, [
        "run", "--seq-dir", str(seq_dir),
        "--mode", "supertree", "--speed", "normal",
        "--dry-run",
    ])
    assert result.exit_code == 0
    assert "genetrees" in result.output.lower() or "gene trees" in result.output.lower()
    assert "wastral" in result.output.lower()
```

- [ ] **Step 2: Run tests to verify they fail**

```bash
uv run pytest tests/cli/test_run.py -k "dry_run" -v
```

Expected: FAIL — NotImplementedError.

- [ ] **Step 3: Implement dry-run output in `execute_pipeline`**

Add a `_STEP_TOOL_LABELS` dict and dry-run printing block inside `execute_pipeline`, replacing the `NotImplementedError` placeholder at the end:

```python
# At module level in _run_pipeline.py:
_STEP_TOOL_LABELS: dict[tuple[str, str, str], str] = {
    # (mode, speed, step_name) -> display label
    ("supermatrix", "normal", "convert"):      "pretree convert",
    ("supermatrix", "normal", "align"):        "pretree align (MAFFT linsi)",
    ("supermatrix", "normal", "trim"):         "pretree trim (trimAl -automated1)",
    ("supermatrix", "normal", "filter_taper"): "pretree filter taper (TAPER)",
    ("supermatrix", "normal", "concat"):       "pretree concat",
    ("supermatrix", "normal", "tree"):         "tree ml iqtree (unpartitioned)",
    ("supermatrix", "fast", "convert"):        "pretree convert",
    ("supermatrix", "fast", "align"):          "pretree align (MAFFT auto)",
    ("supermatrix", "fast", "trim"):           "pretree trim (trimAl -automated1)",
    ("supermatrix", "fast", "concat"):         "pretree concat",
    ("supermatrix", "fast", "tree"):           "tree ml fasttree --matrix",
    ("supertree", "normal", "convert"):        "pretree convert",
    ("supertree", "normal", "align"):          "pretree align (MAFFT linsi)",
    ("supertree", "normal", "trim"):           "pretree trim (trimAl -automated1)",
    ("supertree", "normal", "filter_taper"):   "pretree filter taper (TAPER)",
    ("supertree", "normal", "genetrees"):      "tree ml fasttree --msa-dir (normal)",
    ("supertree", "normal", "tree"):           "tree msc (wASTRAL mode 1)",
    ("supertree", "fast", "convert"):          "pretree convert",
    ("supertree", "fast", "align"):            "pretree align (MAFFT auto)",
    ("supertree", "fast", "trim"):             "pretree trim (trimAl -automated1)",
    ("supertree", "fast", "genetrees"):        "tree ml fasttree --msa-dir (fast)",
    ("supertree", "fast", "tree"):             "tree msc (wASTRAL mode 1)",
}

# Inside execute_pipeline, replace the final NotImplementedError:
    step_defs = _STEP_DEFINITIONS[(mode, speed)]
    total = len(step_defs)

    if dry_run:
        console.print(f"\n[bold]Dry run — {mode} / {speed}[/bold]  ({total} steps)\n")
        for i, defn in enumerate(step_defs, 1):
            label = _STEP_TOOL_LABELS.get((mode, speed, defn["name"]), defn["name"])
            console.print(f"  [{i}/{total}] {defn['name']:15s}  {label}")
        console.print()
        return
```

- [ ] **Step 4: Run tests to verify they pass**

```bash
uv run pytest tests/cli/test_run.py -k "dry_run" -v
```

Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add phyloai/cli/commands/_run_pipeline.py tests/cli/test_run.py
git commit -m "feat(run): implement --dry-run step listing with tool labels"
```

---

## Task 5: Step dispatch — convert, align, trim

**Files:**
- Modify: `phyloai/cli/commands/_run_pipeline.py`
- Test: `tests/cli/test_run.py`

**Interfaces:**
- Consumes:
  - `phyloai.pretree.convert.convert_input(input_path, output_dir, target_format, overwrite, quiet, progress_callback) -> dict`
  - `phyloai.pretree.align.run_align(seq_dir, output_dir, method, seq_type, threads, overwrite, resume, dry_run, quiet, progress_callback) -> dict`
  - `phyloai.pretree.trim.run_trim(msa_dir, output_dir, tool, trimal_method, threads, overwrite, resume, dry_run, quiet, progress_callback) -> dict`
- Produces: `_dispatch_step()` helper; convert/align/trim steps run and their `result.json` is written by the library functions

- [ ] **Step 1: Write the failing integration test (dry-run skips tool execution, so use mocks for real execution)**

```python
# tests/cli/test_run.py  (append)
from unittest.mock import patch, MagicMock


def _mock_step_result(n_files: int = 2) -> dict:
    return {
        "status": "success",
        "command": "phyloai pretree convert ...",
        "wall_time": 1.0,
        "tool_versions": {},
        "params": {},
        "key_results": {},
        "error": None,
        "data": {"files": [{"input": f"g{i}.fa"} for i in range(n_files)]},
    }


def test_run_calls_convert_and_align_and_trim(tmp_path: Path) -> None:
    """Verify execute_pipeline calls convert, align, trim in order."""
    runner = CliRunner()
    seq_dir = _make_seq_dir(tmp_path)
    out_dir = tmp_path / "run"

    with patch("phyloai.pretree.convert.convert_input", return_value=_mock_step_result()) as mock_conv, \
         patch("phyloai.pretree.align.run_align", return_value=_mock_step_result()) as mock_align, \
         patch("phyloai.pretree.trim.run_trim", return_value=_mock_step_result()) as mock_trim, \
         patch("phyloai.pretree.filter.run_taper", side_effect=NotImplementedError) as _mock_taper, \
         patch("phyloai.pretree.concat.run_concat", side_effect=NotImplementedError) as _mock_concat, \
         patch("phyloai.tree.ml_iqtree.run_iqtree", side_effect=NotImplementedError) as _mock_iq:
        result = runner.invoke(cli, [
            "run", "--seq-dir", str(seq_dir),
            "--mode", "supermatrix", "--speed", "normal",
            "--output-dir", str(out_dir),
        ])

    # Should fail at filter step (NotImplementedError not yet implemented),
    # but convert/align/trim must have been called before that.
    assert mock_conv.called
    assert mock_align.called
    assert mock_trim.called
    # Step dirs created
    assert (out_dir / "1-convert").exists()
    assert (out_dir / "2-align").exists()
    assert (out_dir / "3-trim").exists()
```

- [ ] **Step 2: Run test to verify it fails**

```bash
uv run pytest tests/cli/test_run.py::test_run_calls_convert_and_align_and_trim -v
```

Expected: FAIL — NotImplementedError from step dispatch.

- [ ] **Step 3: Implement `_dispatch_step` and the convert/align/trim steps**

Add to `_run_pipeline.py`, and replace the final `NotImplementedError` at the end of `execute_pipeline` with the step loop:

```python
# --- Step dispatch helper ---
def _dispatch_step(
    *,
    checkpoint: dict[str, Any],
    checkpoint_path: Path,
    step_index: int,
    step_name: str,
    step_label: str,
    step_number: int,
    total_steps: int,
    console: Any,
    quiet: bool,
    runner: Any,  # callable that runs the step, returns result dict
) -> dict[str, Any]:
    """Run one step, updating checkpoint before and after."""
    step = checkpoint["steps"][step_index]

    # Check if already done (resume)
    if step["status"] == "success":
        result_path = Path(step["output_dir"]) / "result.json"
        if result_path.exists():
            try:
                data = json.loads(result_path.read_text())
                if data.get("status") == "success":
                    if not quiet:
                        console.print(f"  [{step_number}/{total_steps}] {step_label}  [dim](already done, skipping)[/dim]")
                    return data
            except Exception:
                pass

    if not quiet:
        console.print(f"\n[bold][{step_number}/{total_steps}][/bold] {step_label} ...")

    # Mark running
    step["status"] = "running"
    _save_run_checkpoint(checkpoint, checkpoint_path)

    try:
        result = runner()
    except click.ClickException:
        step["status"] = "failed"
        _save_run_checkpoint(checkpoint, checkpoint_path, fsync=True)
        raise
    except Exception as exc:
        step["status"] = "failed"
        _save_run_checkpoint(checkpoint, checkpoint_path, fsync=True)
        raise click.ClickException(f"Step '{step_name}' failed: {exc}") from exc

    step["status"] = "success"
    _save_run_checkpoint(checkpoint, checkpoint_path)
    return result


# --- Inside execute_pipeline, replace the final section (after dry_run block) ---
    # Set output_dir on steps if not already set (fresh run)
    if not resume:
        for i, defn in enumerate(step_defs):
            checkpoint["steps"][i]["output_dir"] = str(output_dir / defn["subdir"])
        _save_run_checkpoint(checkpoint, checkpoint_path)

    # Build step index map for easy lookup
    step_map = {s["name"]: (i, s) for i, s in enumerate(checkpoint["steps"])}

    import time as _time
    run_start = _time.monotonic()

    # Collect tool versions across steps
    all_tool_versions: dict[str, str] = {}

    from phyloai.pretree.convert import convert_input
    from phyloai.pretree.align import run_align
    from phyloai.pretree.trim import run_trim

    # --- Step 1: Convert ---
    step_idx, step_info = step_map["convert"]
    step_out = Path(step_info["output_dir"])
    step_out.mkdir(parents=True, exist_ok=True)

    def _run_convert() -> dict[str, Any]:
        return convert_input(
            input_path=seq_dir,
            output_dir=step_out,
            target_format="fasta",
            overwrite=True,  # always overwrite within run (step dir managed by run)
            quiet=quiet,
        )

    convert_result = _dispatch_step(
        checkpoint=checkpoint,
        checkpoint_path=checkpoint_path,
        step_index=step_idx,
        step_name="convert",
        step_label=_STEP_TOOL_LABELS[(mode, speed, "convert")],
        step_number=1,
        total_steps=total,
        console=console,
        quiet=quiet,
        runner=_run_convert,
    )
    all_tool_versions.update(convert_result.get("tool_versions") or {})
    converted_seqs_dir = step_out / "seqs"

    # --- Step 2: Align ---
    step_idx, step_info = step_map["align"]
    step_out = Path(step_info["output_dir"])
    step_out.mkdir(parents=True, exist_ok=True)
    align_method = "linsi" if speed == "normal" else "auto"

    def _run_align_step() -> dict[str, Any]:
        return run_align(
            seq_dir=converted_seqs_dir,
            output_dir=step_out,
            method=align_method,
            seq_type="auto",
            threads=threads,
            overwrite=True,
            resume=(step_info["status"] == "interrupted"),
            quiet=quiet,
        )

    align_result = _dispatch_step(
        checkpoint=checkpoint,
        checkpoint_path=checkpoint_path,
        step_index=step_idx,
        step_name="align",
        step_label=_STEP_TOOL_LABELS[(mode, speed, "align")],
        step_number=2,
        total_steps=total,
        console=console,
        quiet=quiet,
        runner=_run_align_step,
    )
    all_tool_versions.update(align_result.get("tool_versions") or {})
    aligned_seqs_dir = step_out / "seqs"

    # --- Step 3: Trim ---
    step_idx, step_info = step_map["trim"]
    step_out = Path(step_info["output_dir"])
    step_out.mkdir(parents=True, exist_ok=True)

    def _run_trim_step() -> dict[str, Any]:
        return run_trim(
            msa_dir=aligned_seqs_dir,
            output_dir=step_out,
            tool="trimal",
            trimal_method="automated1",
            threads=threads,
            overwrite=True,
            resume=(step_info["status"] == "interrupted"),
            quiet=quiet,
        )

    trim_result = _dispatch_step(
        checkpoint=checkpoint,
        checkpoint_path=checkpoint_path,
        step_index=step_idx,
        step_name="trim",
        step_label=_STEP_TOOL_LABELS[(mode, speed, "trim")],
        step_number=3,
        total_steps=total,
        console=console,
        quiet=quiet,
        runner=_run_trim_step,
    )
    all_tool_versions.update(trim_result.get("tool_versions") or {})
    trimmed_seqs_dir = step_out / "seqs"

    # Placeholder: remaining steps
    raise NotImplementedError("Remaining steps not yet implemented")
```

- [ ] **Step 4: Run test to verify it passes**

```bash
uv run pytest tests/cli/test_run.py::test_run_calls_convert_and_align_and_trim -v
```

Expected: PASS (crashes at filter step with NotImplementedError as expected by the test).

- [ ] **Step 5: Commit**

```bash
git add phyloai/cli/commands/_run_pipeline.py tests/cli/test_run.py
git commit -m "feat(run): implement convert/align/trim step dispatch"
```

---

## Task 6: Step dispatch — filter taper, concat / gene trees, tree

**Files:**
- Modify: `phyloai/cli/commands/_run_pipeline.py`
- Test: `tests/cli/test_run.py`

**Interfaces:**
- Consumes:
  - `phyloai.pretree.filter.run_taper(msa_dir, output_dir, threads, overwrite, quiet) -> dict`
  - `phyloai.pretree.concat.run_concat(msa_dir, output_dir, threads, overwrite, quiet) -> dict`
  - `phyloai.tree.ml.run_fasttree(msa_dir, output_dir, mode, threads, overwrite, resume, quiet) -> dict`
  - `phyloai.tree.ml.run_fasttree(matrix, output_dir, mode, threads, overwrite, quiet) -> dict`
  - `phyloai.tree.ml_iqtree.run_iqtree(matrix, output_dir, threads, overwrite, resume, quiet) -> dict`
  - `phyloai.tree.msc.run_wastral(tree_dir, output_dir, mode, threads, overwrite, quiet) -> dict`
- Produces: full pipeline runs end-to-end; `result.json` written at `<output_dir>/result.json`

- [ ] **Step 1: Write the failing end-to-end mock test**

```python
# tests/cli/test_run.py  (append)
def test_run_supermatrix_normal_full_pipeline_mocked(tmp_path: Path) -> None:
    """All six steps called in order; result.json written at output root."""
    runner = CliRunner()
    seq_dir = _make_seq_dir(tmp_path)
    out_dir = tmp_path / "run"

    r = _mock_step_result()

    with patch("phyloai.pretree.convert.convert_input", return_value=r), \
         patch("phyloai.pretree.align.run_align", return_value=r), \
         patch("phyloai.pretree.trim.run_trim", return_value=r), \
         patch("phyloai.pretree.filter.run_taper", return_value=r) as mock_taper, \
         patch("phyloai.pretree.concat.run_concat", return_value={**r, "data": {"matrix_file": str(out_dir / "5-concat" / "matrix.fa"), "n_taxa": 2, "n_sites": 10}}) as mock_concat, \
         patch("phyloai.tree.ml_iqtree.run_iqtree", return_value={**r, "data": {"tree_file": str(out_dir / "6-tree" / "iqtree.treefile")}}) as mock_iqtree:
        result = runner.invoke(cli, [
            "run", "--seq-dir", str(seq_dir),
            "--mode", "supermatrix", "--speed", "normal",
            "--output-dir", str(out_dir),
        ])

    assert result.exit_code == 0, result.output
    assert mock_taper.called
    assert mock_concat.called
    assert mock_iqtree.called
    # result.json written
    result_json = out_dir / "result.json"
    assert result_json.exists()
    data = json.loads(result_json.read_text())
    assert data["status"] == "success"
    assert data["key_results"]["n_input_genes"] >= 1
    assert "final_tree" in data["key_results"]


def test_run_supertree_fast_full_pipeline_mocked(tmp_path: Path) -> None:
    """Supertree fast: no filter, fasttree --msa-dir, wastral. 5 steps."""
    runner = CliRunner()
    seq_dir = _make_seq_dir(tmp_path)
    out_dir = tmp_path / "run"

    r = _mock_step_result()

    with patch("phyloai.pretree.convert.convert_input", return_value=r), \
         patch("phyloai.pretree.align.run_align", return_value=r), \
         patch("phyloai.pretree.trim.run_trim", return_value=r), \
         patch("phyloai.pretree.filter.run_taper", side_effect=AssertionError("should not be called")) as mock_taper, \
         patch("phyloai.tree.ml.run_fasttree", return_value={**r, "data": {"trees_dir": str(out_dir / "5-genetrees" / "trees")}}) as mock_ft, \
         patch("phyloai.tree.msc.run_wastral", return_value={**r, "data": {"tree_file": str(out_dir / "6-tree" / "wastral.tre")}}) as mock_wastral:
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
    # 5-genetrees present, not 5-concat
    assert (out_dir / "5-genetrees").exists()
    assert not (out_dir / "5-concat").exists()
    assert not (out_dir / "4-filter").exists()
```

- [ ] **Step 2: Run tests to verify they fail**

```bash
uv run pytest tests/cli/test_run.py -k "full_pipeline_mocked" -v
```

Expected: FAIL — NotImplementedError.

- [ ] **Step 3: Implement remaining steps in `execute_pipeline`**

Replace the final `raise NotImplementedError(...)` at the end of `execute_pipeline` in `_run_pipeline.py` with the following block (continues after `trimmed_seqs_dir`):

```python
    from phyloai.pretree.filter import run_taper
    from phyloai.pretree.concat import run_concat
    from phyloai.tree.ml import run_fasttree
    from phyloai.tree.ml_iqtree import run_iqtree
    from phyloai.tree.msc import run_wastral

    # Input to concat or gene-trees: after filter (normal) or trim (fast)
    filtered_seqs_dir = trimmed_seqs_dir  # updated below if filter runs

    step_number = 4  # steps 1-3 already done above

    # --- Step 4: Filter (taper) — normal mode only ---
    if "filter_taper" in step_map:
        step_idx, step_info = step_map["filter_taper"]
        step_out = Path(step_info["output_dir"])
        step_out.mkdir(parents=True, exist_ok=True)

        def _run_taper() -> dict[str, Any]:
            return run_taper(
                msa_dir=trimmed_seqs_dir,
                output_dir=step_out,
                threads=threads,
                overwrite=True,
                quiet=quiet,
            )

        taper_result = _dispatch_step(
            checkpoint=checkpoint,
            checkpoint_path=checkpoint_path,
            step_index=step_idx,
            step_name="filter_taper",
            step_label=_STEP_TOOL_LABELS[(mode, speed, "filter_taper")],
            step_number=step_number,
            total_steps=total,
            console=console,
            quiet=quiet,
            runner=_run_taper,
        )
        all_tool_versions.update(taper_result.get("tool_versions") or {})
        filtered_seqs_dir = step_out / "seqs"
        step_number += 1

    n_genes_after_filter: int = len(list(filtered_seqs_dir.glob("*.fa"))) if filtered_seqs_dir.exists() else 0
    final_tree_path: str = ""

    # --- Branch: supermatrix vs supertree ---
    if mode == "supermatrix":
        # --- Step 5: Concat ---
        step_idx, step_info = step_map["concat"]
        step_out = Path(step_info["output_dir"])
        step_out.mkdir(parents=True, exist_ok=True)

        def _run_concat() -> dict[str, Any]:
            return run_concat(
                msa_dir=filtered_seqs_dir,
                output_dir=step_out,
                overwrite=True,
                quiet=quiet,
            )

        concat_result = _dispatch_step(
            checkpoint=checkpoint,
            checkpoint_path=checkpoint_path,
            step_index=step_idx,
            step_name="concat",
            step_label=_STEP_TOOL_LABELS[(mode, speed, "concat")],
            step_number=step_number,
            total_steps=total,
            console=console,
            quiet=quiet,
            runner=_run_concat,
        )
        all_tool_versions.update(concat_result.get("tool_versions") or {})
        matrix_file = Path(step_out / "matrix.fa")
        step_number += 1

        # --- Step 6: Tree (iqtree normal / fasttree fast) ---
        step_idx, step_info = step_map["tree"]
        step_out = Path(step_info["output_dir"])
        step_out.mkdir(parents=True, exist_ok=True)

        if speed == "normal":
            def _run_tree() -> dict[str, Any]:
                return run_iqtree(
                    matrix=matrix_file,
                    output_dir=step_out,
                    threads=threads,
                    overwrite=True,
                    resume=(step_info["status"] == "interrupted"),
                    quiet=quiet,
                )
        else:
            def _run_tree() -> dict[str, Any]:
                return run_fasttree(
                    matrix=matrix_file,
                    output_dir=step_out,
                    mode="normal",
                    threads=threads,
                    overwrite=True,
                    quiet=quiet,
                )

        tree_result = _dispatch_step(
            checkpoint=checkpoint,
            checkpoint_path=checkpoint_path,
            step_index=step_idx,
            step_name="tree",
            step_label=_STEP_TOOL_LABELS[(mode, speed, "tree")],
            step_number=step_number,
            total_steps=total,
            console=console,
            quiet=quiet,
            runner=_run_tree,
        )
        all_tool_versions.update(tree_result.get("tool_versions") or {})
        # Resolve final tree path
        final_tree_path = str(
            tree_result.get("data", {}).get("tree_file")
            or tree_result.get("data", {}).get("treefile")
            or ""
        )

    else:  # supertree
        # --- Step 5: Gene trees (fasttree --msa-dir) ---
        step_idx, step_info = step_map["genetrees"]
        step_out = Path(step_info["output_dir"])
        step_out.mkdir(parents=True, exist_ok=True)
        ft_mode = "normal" if speed == "normal" else "fast"

        def _run_genetrees() -> dict[str, Any]:
            return run_fasttree(
                msa_dir=filtered_seqs_dir,
                output_dir=step_out,
                mode=ft_mode,
                threads=threads,
                overwrite=True,
                resume=(step_info["status"] == "interrupted"),
                quiet=quiet,
            )

        gt_result = _dispatch_step(
            checkpoint=checkpoint,
            checkpoint_path=checkpoint_path,
            step_index=step_idx,
            step_name="genetrees",
            step_label=_STEP_TOOL_LABELS[(mode, speed, "genetrees")],
            step_number=step_number,
            total_steps=total,
            console=console,
            quiet=quiet,
            runner=_run_genetrees,
        )
        all_tool_versions.update(gt_result.get("tool_versions") or {})
        trees_dir = step_out / "trees"
        step_number += 1

        # --- Step 6: Species tree (wASTRAL) ---
        step_idx, step_info = step_map["tree"]
        step_out = Path(step_info["output_dir"])
        step_out.mkdir(parents=True, exist_ok=True)

        def _run_wastral() -> dict[str, Any]:
            return run_wastral(
                tree_dir=trees_dir,
                output_dir=step_out,
                mode=1,
                threads=threads,
                overwrite=True,
                quiet=quiet,
            )

        wastral_result = _dispatch_step(
            checkpoint=checkpoint,
            checkpoint_path=checkpoint_path,
            step_index=step_idx,
            step_name="tree",
            step_label=_STEP_TOOL_LABELS[(mode, speed, "tree")],
            step_number=step_number,
            total_steps=total,
            console=console,
            quiet=quiet,
            runner=_run_wastral,
        )
        all_tool_versions.update(wastral_result.get("tool_versions") or {})
        final_tree_path = str(
            wastral_result.get("data", {}).get("tree_file") or ""
        )

    # --- Mark run complete ---
    import time as _time
    wall_time = round(_time.monotonic() - run_start, 3)

    checkpoint["status"] = "success"
    import datetime as _dt
    checkpoint["completed_at"] = _dt.datetime.now(_dt.timezone.utc).isoformat(timespec="seconds")
    _save_run_checkpoint(checkpoint, checkpoint_path, fsync=True)

    # --- Write result.json ---
    n_input = len(list((output_dir / "1-convert" / "seqs").glob("*.fa"))) if (output_dir / "1-convert" / "seqs").exists() else 0

    result_payload: dict[str, Any] = {
        "status": "success",
        "command": command_str,
        "wall_time": wall_time,
        "tool_versions": all_tool_versions,
        "params": params,
        "key_results": {
            "n_input_genes": n_input,
            "n_genes_after_filter": n_genes_after_filter,
            "final_tree": final_tree_path,
        },
        "error": None,
        "data": {
            "mode": mode,
            "speed": speed,
            "steps": [
                {
                    "name": s["name"],
                    "status": s["status"],
                    "output_dir": s["output_dir"],
                    "result_json": str(Path(s["output_dir"]) / "result.json") if s["output_dir"] else None,
                }
                for s in checkpoint["steps"]
            ],
        },
    }

    result_file = output_dir / "result.json"
    result_file.write_text(
        json.dumps(result_payload, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )

    if not quiet:
        console.print(f"\n[bold green]✓ Pipeline complete[/bold green]  [dim][total: {wall_time:.0f}s][/dim]")
        if final_tree_path:
            console.print(f"  Species tree:  {final_tree_path}")
        console.print(f"  Results:       {result_file}")
```

- [ ] **Step 4: Run tests to verify they pass**

```bash
uv run pytest tests/cli/test_run.py -k "full_pipeline_mocked" -v
```

Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add phyloai/cli/commands/_run_pipeline.py tests/cli/test_run.py
git commit -m "feat(run): implement full step dispatch for all modes and speeds"
```

---

## Task 7: Error handling — step failure writes `result.json` with `status: error`

**Files:**
- Modify: `phyloai/cli/commands/_run_pipeline.py`
- Test: `tests/cli/test_run.py`

**Interfaces:**
- Consumes: `execute_pipeline` (existing)
- Produces: on step failure, `result.json` written with `status: "error"`, exit code 2

- [ ] **Step 1: Write the failing tests**

```python
# tests/cli/test_run.py  (append)
def test_run_step_failure_writes_error_result_json(tmp_path: Path) -> None:
    runner = CliRunner()
    seq_dir = _make_seq_dir(tmp_path)
    out_dir = tmp_path / "run"

    r = _mock_step_result()

    with patch("phyloai.pretree.convert.convert_input", return_value=r), \
         patch("phyloai.pretree.align.run_align", side_effect=RuntimeError("MAFFT not found")):
        result = runner.invoke(cli, [
            "run", "--seq-dir", str(seq_dir),
            "--output-dir", str(out_dir),
        ])

    assert result.exit_code != 0
    result_json_path = out_dir / "result.json"
    assert result_json_path.exists()
    data = json.loads(result_json_path.read_text())
    assert data["status"] == "error"
    assert "align" in (data["error"] or "").lower() or "mafft" in (data["error"] or "").lower()
```

- [ ] **Step 2: Run test to verify it fails**

```bash
uv run pytest tests/cli/test_run.py::test_run_step_failure_writes_error_result_json -v
```

Expected: FAIL — no `result.json` written on failure.

- [ ] **Step 3: Add error `result.json` write and error handling to `execute_pipeline`**

Wrap the entire step-dispatch section in a try/except inside `execute_pipeline`. Add this after the import section, wrapping from the convert step through the final `result.json` write:

```python
    # --- Pipeline execution (wrapped for error result.json) ---
    import time as _time
    run_start = _time.monotonic()
    all_tool_versions: dict[str, str] = {}

    def _write_error_result(error_msg: str) -> None:
        wall_time = round(_time.monotonic() - run_start, 3)
        payload: dict[str, Any] = {
            "status": "error",
            "command": command_str,
            "wall_time": wall_time,
            "tool_versions": all_tool_versions,
            "params": params,
            "key_results": {},
            "error": error_msg,
            "data": {
                "mode": mode,
                "speed": speed,
                "steps": [
                    {
                        "name": s["name"],
                        "status": s["status"],
                        "output_dir": s["output_dir"],
                        "result_json": str(Path(s["output_dir"]) / "result.json") if s["output_dir"] else None,
                    }
                    for s in checkpoint["steps"]
                ],
            },
        }
        if not dry_run:
            result_file = output_dir / "result.json"
            result_file.write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")

    try:
        # ... all step dispatch code ...
    except click.ClickException as exc:
        _write_error_result(exc.format_message())
        raise
    except Exception as exc:
        _write_error_result(str(exc))
        raise click.ClickException(str(exc))
```

Note: refactor `execute_pipeline` to move the `run_start = _time.monotonic()` and `all_tool_versions` declarations before the try block, and move `_write_error_result` definition before the try block. The existing success path's `wall_time = ...` stays inside the try block.

- [ ] **Step 4: Run test to verify it passes**

```bash
uv run pytest tests/cli/test_run.py::test_run_step_failure_writes_error_result_json -v
```

Expected: PASS.

- [ ] **Step 5: Run all run tests**

```bash
uv run pytest tests/cli/test_run.py -v
```

Expected: all PASS.

- [ ] **Step 6: Commit**

```bash
git add phyloai/cli/commands/_run_pipeline.py tests/cli/test_run.py
git commit -m "feat(run): write error result.json on step failure"
```

---

## Task 8: Resume flow — skip completed steps

**Files:**
- Modify: `phyloai/cli/commands/_run_pipeline.py`
- Test: `tests/cli/test_run.py`

**Interfaces:**
- Consumes: `_dispatch_step` (already checks `status == "success"` and verifies `result.json`)
- Produces: `--resume` skips steps whose `run_checkpoint.json` entry is `success` with valid `result.json`

- [ ] **Step 1: Write the failing tests**

```python
# tests/cli/test_run.py  (append)
def test_run_resume_skips_completed_steps(tmp_path: Path) -> None:
    """If run_checkpoint shows convert+align done, only remaining steps run."""
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

    with patch("phyloai.pretree.convert.convert_input", mock_convert), \
         patch("phyloai.pretree.align.run_align", mock_align), \
         patch("phyloai.pretree.trim.run_trim", return_value=r), \
         patch("phyloai.pretree.filter.run_taper", return_value=r), \
         patch("phyloai.pretree.concat.run_concat", return_value=r), \
         patch("phyloai.tree.ml_iqtree.run_iqtree", return_value=r):
        result = runner.invoke(cli, [
            "run", "--seq-dir", str(seq_dir),
            "--output-dir", str(out_dir),
            "--resume",
        ])

    assert result.exit_code == 0, result.output
    # convert and align should NOT have been called again
    assert not mock_convert.called
    assert not mock_align.called
```

- [ ] **Step 2: Run test to verify it fails**

```bash
uv run pytest tests/cli/test_run.py::test_run_resume_skips_completed_steps -v
```

Expected: FAIL — resume doesn't fully work yet.

- [ ] **Step 3: Fix step output_dir assignment on resume**

In `execute_pipeline`, after `checkpoint = _load_run_checkpoint(checkpoint_path)` and `_validate_run_resume(...)`, ensure step output dirs are set from the checkpoint (they already are in the JSON) — the key is that on resume, we do NOT re-set `output_dir` on steps. Verify the resume path doesn't overwrite step output dirs. The existing logic (only setting dirs on `not resume`) already handles this — ensure the resume branch skips the step-dir assignment block.

The `_dispatch_step` helper already checks `step["status"] == "success"` and validates `result.json`. This test passes once the resume branch correctly feeds the loaded checkpoint's step list with existing `output_dir` values.

- [ ] **Step 4: Run test to verify it passes**

```bash
uv run pytest tests/cli/test_run.py::test_run_resume_skips_completed_steps -v
```

Expected: PASS.

- [ ] **Step 5: Run all run tests**

```bash
uv run pytest tests/cli/test_run.py -v
```

Expected: all PASS.

- [ ] **Step 6: Commit**

```bash
git add phyloai/cli/commands/_run_pipeline.py tests/cli/test_run.py
git commit -m "feat(run): verify --resume skips already-completed steps"
```

---

## Task 9: Main design update and `docs/commands/run.md`

**Files:**
- Modify: `docs/superpowers/specs/2026-06-07-phyloai-design.md`
- Create: `docs/commands/run.md`

**Interfaces:** Documentation only.

- [ ] **Step 1: Update `2026-06-07-phyloai-design.md`**

In Section 4.1, update the `phyloai run` examples to replace `--mode coalescent` with `--mode supertree`:

```markdown
# One-click pipeline
phyloai run --seq-dir ./markers --output-dir ./runs/run --mode supermatrix
phyloai run --seq-dir ./markers --output-dir ./runs/run --mode supertree
```

In Section 4.2, replace the pipeline table with:

```markdown
| Mode | Steps | Notes |
|------|-------|-------|
| `--mode supermatrix` | convert → align → trim → filter taper → concat → iqtree (unpartitioned) | `--speed normal` (default): MAFFT linsi + TAPER + IQ-TREE3; `--speed fast`: MAFFT auto, no TAPER, FastTree |
| `--mode supertree` | convert → align → trim → filter taper → gene trees → wastral | `--speed normal`: MAFFT linsi + TAPER + FastTree gene trees + wASTRAL; `--speed fast`: no TAPER, FastTree fast mode |
```

In Section 9.2, add `--speed` to the shared parameter registry:

```markdown
| `--speed` | | `normal\|fast` | `normal` | `run` only |
```

- [ ] **Step 2: Create `docs/commands/run.md`**

```markdown
# `phyloai run`

## Purpose

One-click phylogenomics pipeline from raw sequence files to a species tree. Orchestrates all preprocessing and inference steps using sensible defaults. For fine-grained control over any individual step, use the constituent subcommands directly.

## Usage

```bash
phyloai run --seq-dir ./markers [OPTIONS]
```

## Inputs

- `--seq-dir`: Directory of raw sequence files in any format (FASTA, Nexus, PHYLIP, Phylip-PAML). All files are converted to normalized FASTA as the first step.

## Pipeline Modes

### `--mode supermatrix` (default)

```
convert → align → trim → [filter taper] → concat → species tree
```

| Speed | Align | Filter | Tree |
|-------|-------|--------|------|
| `normal` | MAFFT linsi | TAPER | IQ-TREE3 (unpartitioned) |
| `fast` | MAFFT auto | skipped | FastTree |

### `--mode supertree`

```
convert → align → trim → [filter taper] → gene trees → species tree
```

| Speed | Align | Filter | Gene Trees | Species Tree |
|-------|-------|--------|------------|--------------|
| `normal` | MAFFT linsi | TAPER | FastTree (normal) | wASTRAL (mode 1) |
| `fast` | MAFFT auto | skipped | FastTree (fast) | wASTRAL (mode 1) |

## Options

| Option | Default | Description |
|--------|---------|-------------|
| `--seq-dir PATH` | required | Input sequence directory |
| `--mode supermatrix\|supertree` | `supermatrix` | Pipeline mode |
| `--speed normal\|fast` | `normal` | Speed/accuracy trade-off |
| `-o, --output-dir PATH` | `./runs/run` | Output directory root |
| `-t, --threads INT` | 4 | Thread count for all steps |
| `--resume` | off | Resume from checkpoint |
| `--overwrite` | off | Overwrite output directory |
| `--dry-run` | off | Show steps without running |
| `-q, --quiet` | off | Suppress non-error output |

## Outputs

```
runs/run/
├── run_checkpoint.json
├── result.json
├── 1-convert/
├── 2-align/
├── 3-trim/
├── 4-filter/          (--speed normal only)
├── 5-concat/          (supermatrix) or 5-genetrees/ (supertree)
└── 6-tree/
```

Each subdirectory contains its own `result.json` with detailed step results.

## Resume Behaviour

`--resume` loads `run_checkpoint.json`, verifies that parameters match exactly, and skips steps whose output is already validated. Steps that were interrupted continue from their own subcommand checkpoint where supported.

`--resume` and `--overwrite` are mutually exclusive.

## Examples

```bash
# Default: supermatrix, normal speed
phyloai run --seq-dir ./markers

# Supertree with fast speed and 16 threads
phyloai run --seq-dir ./markers --mode supertree --speed fast --threads 16

# Resume a previously interrupted run
phyloai run --seq-dir ./markers --output-dir ./runs/run --resume

# Preview steps without running
phyloai run --seq-dir ./markers --mode supertree --dry-run
```

## Warnings / Errors

- Exit 1: `--seq-dir` does not exist, `--resume` without checkpoint, `--resume` + `--overwrite` together, non-empty output directory without `--overwrite`.
- Exit 2: A pipeline step (external tool) failed. Check the failing step's `result.json` and logs.
- Exit 3: A required tool is not installed. Run `phyloai doctor` to check your environment.

## Notes

- `phyloai run` uses each step's default parameters. For non-default settings (e.g. partitioned IQ-TREE, custom TAPER cutoff), run steps individually via `phyloai pretree align`, `phyloai tree ml iqtree`, etc.
- IQ-TREE3 in supermatrix normal mode runs automatic model selection (ModelFinder) without a partition file. This is a first-pass result; partitioned analyses require `phyloai tree ml iqtree` directly.
- The `4-filter/` directory is not created in `--speed fast` mode.
```

- [ ] **Step 3: Commit**

```bash
git add docs/superpowers/specs/2026-06-07-phyloai-design.md docs/commands/run.md
git commit -m "docs: update main design for run command; add docs/commands/run.md"
```

---

## Task 10: Final verification

- [ ] **Step 1: Run the full test suite**

```bash
uv run pytest tests/cli/test_run.py -v
```

Expected: all PASS.

- [ ] **Step 2: Verify `--help` output**

```bash
uv run phyloai run --help
```

Verify: `--seq-dir`, `--mode`, `--speed`, `--resume`, `--overwrite`, `--dry-run`, `--threads`, `--quiet` all present; mode step sequences described; examples present.

- [ ] **Step 3: Verify `run` appears in top-level help**

```bash
uv run phyloai --help
```

Verify: `run` listed between `doctor` and `pretree`.

- [ ] **Step 4: Smoke test with `--dry-run`**

```bash
mkdir -p /tmp/test_markers
echo ">sp1\nMKTLL\n>sp2\nMKTAA" > /tmp/test_markers/gene1.fa
uv run phyloai run --seq-dir /tmp/test_markers --dry-run
uv run phyloai run --seq-dir /tmp/test_markers --mode supertree --speed fast --dry-run
```

Verify: step list printed, no tool execution, exit 0.

- [ ] **Step 5: Final commit**

```bash
git add -A
git commit -m "feat(run): phyloai run pipeline complete (supermatrix/supertree, normal/fast, resume)"
```
