# PhyloAI Report Module Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build the `phyloai report` module that scans a run directory, discovers and orders `result.json` files, generates `methods_text` per step, and produces `report.json` (machine-readable) + `report.html` (human-readable with embedded figures).

**Architecture:** Four internal modules with a single CLI entry point. `collector.py` scans directories and discovers `result.json` files. `schema.py` defines data structures and assembles `report.json`. `templates.py` generates per-step academic methods text. `renderer.py` converts `report.json` to self-contained HTML via Jinja2. The CLI (`cli/commands/report.py`) wires them together.

**Tech Stack:** Python 3.10+, Jinja2 (new dep), Click 8+, Rich, pathlib, json, datetime

**Depends on:** All analysis phases (pretree, tree, posttree) finalized — their `result.json` output must include `data.output_files` conforming to JSON Output Standard Section 5.4.

## Global Constraints

- Python >=3.10
- All non-`doctor` commands write `result.json`; report reads these files only
- `report.json` is source of truth; `report.html` is fully derived from it
- No LLM involvement in methods text generation
- No hardcoded file paths in report module
- Report always succeeds on incomplete runs (failed steps included with error details)
- `report` is a single Click command, not a command group with sub-commands

---

## File Map

| Action | File | Responsibility |
|--------|------|----------------|
| Create | `phyloai/report/__init__.py` | Package marker, version |
| Create | `phyloai/report/collector.py` | Directory scanning, `result.json` discovery, step ordering, `step_id` parsing, run_mode detection |
| Create | `phyloai/report/schema.py` | `ReportStep`, `ReportRecord` dataclasses; `assemble_report()` to build `report.json` |
| Create | `phyloai/report/templates.py` | Per-command `generate_methods_<step_id>()` functions producing 2-5 sentence methods text |
| Create | `phyloai/report/renderer.py` | `render_html()`: reads `report.json`, renders via Jinja2 to `report.html` |
| Create | `phyloai/report/html/report.html.j2` | Self-contained Jinja2 HTML template (5 panels, collapsible cards, embedded PDF figures, sortable tables) |
| Create | `phyloai/cli/commands/report.py` | `phyloai report` Click command: wires collector → schema → renderer |
| Modify | `phyloai/cli/main.py` | Register `report` command |
| Modify | `pyproject.toml` | Add `jinja2>=3.1` dependency |
| Create | `tests/report/__init__.py` | Test package marker |
| Create | `tests/report/test_collector.py` | Collector unit tests: run_mode detection, step discovery, step_id parsing, ordering |
| Create | `tests/report/test_schema.py` | Schema tests: ReportRecord assembly, figure/table indexing, incomplete run handling |
| Create | `tests/report/test_templates.py` | Template tests: each template function produces non-empty text, handles missing keys, conditional branches |
| Create | `tests/report/test_renderer.py` | Renderer smoke test: valid report.json → valid HTML output |
| Create | `tests/report/test_integration.py` | End-to-end: mock run directory → `phyloai report` → validate `report.json` + `report.html` |

---

## Task 1: Project Scaffold

**Files:**
- Create: `phyloai/report/__init__.py`
- Create: `phyloai/report/html/report.html.j2` (placeholder)
- Create: `tests/report/__init__.py`
- Modify: `pyproject.toml`

- [ ] **Step 1: Add Jinja2 dependency to pyproject.toml**

Add `"jinja2>=3.1"` to `[project].dependencies`:

```toml
dependencies = [
    "click>=8.1",
    "rich>=13.0",
    "biopython>=1.81",
    "clipkit>=2.0",
    "numpy>=1.24",
    "matplotlib>=3.7",
    "scipy>=1.10",
    "scikit-learn>=1.3.0",
    "jinja2>=3.1",
]
```

- [ ] **Step 2: Install the new dependency**

Run: `pip install -e ".[dev]"` or `pip install jinja2>=3.1`

- [ ] **Step 3: Create `phyloai/report/__init__.py`**

```python
"""PhyloAI report module — generates reproducible analysis reports."""

from phyloai import __version__
```

- [ ] **Step 4: Create `phyloai/report/html/` directory and placeholder template**

```bash
mkdir -p phyloai/report/html
```

Create `phyloai/report/html/report.html.j2` with a minimal Jinja2 template:

```html
<!DOCTYPE html>
<html lang="en">
<head><meta charset="UTF-8"><title>PhyloAI Report</title></head>
<body><h1>PhyloAI Report</h1><p>Template placeholder — will be expanded in Task 5.</p></body>
</html>
```

- [ ] **Step 5: Create `tests/report/__init__.py`**

```python
"""Tests for phyloai.report module."""
```

- [ ] **Step 6: Verify scaffold**

Run: `python -c "import phyloai.report; print('OK')"`
Expected: prints `OK` with no errors.

---

## Task 2: Collector — Directory Scanning and Step Discovery

**Files:**
- Create: `phyloai/report/collector.py`
- Create: `tests/report/test_collector.py`

**Interfaces:**
- Produces: `discover_steps(run_dir: Path) -> dict` — returns `{"run_mode": str, "steps": list[dict], "pipeline_summary": dict|None}`
- Produces: `parse_step_id(command: str) -> str` — e.g. `"phyloai pretree align ..."` → `"pretree.align"`
- Produces: `STEP_ORDER: list[str]` — canonical step ordering list

### Collector algorithm

**run_mode detection** (priority-ordered checks):

| Priority | Condition | run_mode |
|----------|-----------|----------|
| 1 | `run_dir/result.json` exists AND at least one subdir also contains `result.json` | `pipeline` |
| 2 | `run_dir/result.json` exists, no subdirs with `result.json` | `module` (single cmd) |
| 3 | No `run_dir/result.json`, but subdirs contain `result.json` | `module` (multi-step) |
| 4 | No `result.json` found at any expected depth | Error |

**Scan depth:** module = 1-2 levels; pipeline = 2-3 levels.

**Step discovery:** Uses structural depth-based scanning only (Spec Section 5 priority rules). For both `pipeline` and `module` modes, scan subdirectories up to 3 levels deep for `result.json` files. `pipeline` mode additionally reads `data.mode` and `data.speed` from the top-level `result.json` for metadata enrichment, but does NOT require `data.steps[]` to exist — step paths come purely from filesystem scan. This avoids coupling report discovery to the internal `phyloai run` data format.

**step_id parsing:** Extracts the positional subcommand tokens after `phyloai` by filtering out flag tokens (any token starting with `-`), then matching against the known CLI command tree (`pretree`/`tree`/`posttree`/`run`):
- `"phyloai pretree align --seq-dir ./raw --method linsi"` → `"pretree.align"`
- `"phyloai pretree filter taper --msa-dir ./trimmed"` → `"pretree.filter.taper"`
- `"phyloai tree ml iqtree --matrix ./concat/matrix.fa --model C20"` → `"tree.ml.iqtree"`
- `"phyloai posttree dating mcmc --tree ./tree.nwk"` → `"posttree.dating.mcmc"`
- `"phyloai run --seq-dir ./raw --mode supermatrix"` → `"run"`

**STEP_ORDER** (from spec Section 6):

```python
STEP_ORDER = [
    "pretree.convert",
    "pretree.stats",
    "pretree.align",
    "pretree.trim",
    "pretree.metrics",
    "pretree.filter.taper",
    "pretree.filter.treeshrink",
    "pretree.filter.symtest",
    "pretree.filter.metrics",
    "pretree.filter.cluster",
    "pretree.concat",
    "tree.ml.fasttree",
    "tree.ml.iqtree",
    "tree.msc",
    "tree.bi",
    "tree.cf",
    "posttree.topology",
    "posttree.dating.hessian",
    "posttree.dating.mcmc",
    "posttree.signal",
    "posttree.syserror.brlen",
    "posttree.syserror.cca",
    "posttree.syserror.sites",
    "posttree.simulate",
]
```

Steps not in `STEP_ORDER` are appended at the end.

- [ ] **Step 1: Write failing tests for collector**

In `tests/report/test_collector.py`:

```python
"""Tests for phyloai.report.collector."""
from __future__ import annotations

import json
from pathlib import Path

from phyloai.report.collector import (
    STEP_ORDER,
    discover_steps,
    parse_step_id,
)


class TestParseStepId:
    def test_basic_two_part(self):
        assert parse_step_id("phyloai pretree align --seq-dir ./raw --method linsi") == "pretree.align"

    def test_three_part_filter(self):
        assert parse_step_id("phyloai pretree filter taper --msa-dir ./trimmed") == "pretree.filter.taper"

    def test_tree_ml_iqtree(self):
        assert parse_step_id("phyloai tree ml iqtree --matrix ./concat/matrix.fa --model C20") == "tree.ml.iqtree"

    def test_pipeline_run(self):
        assert parse_step_id("phyloai run --seq-dir ./raw --mode supermatrix") == "run"

    def test_dating_mcmc(self):
        assert parse_step_id("phyloai posttree dating mcmc --tree ./tree.nwk --matrix ./matrix.fa") == "posttree.dating.mcmc"

    def test_with_paths_and_flags(self):
        assert parse_step_id("phyloai pretree trim --msa-dir /abs/path/to/trimmed --tool bmge --bmge-matrix BLOSUM90 --tool-args '-g 0.5'") == "pretree.trim"

    def test_global_flag_before_subcommand(self):
        """Flag values before the root command are skipped correctly.
        Boolean flags (--quiet) do NOT consume the next token as a value."""
        assert parse_step_id("phyloai --run-dir ./runs pretree align --method linsi") == "pretree.align"
        assert parse_step_id("phyloai --quiet tree ml iqtree --model LG") == "tree.ml.iqtree"
        assert parse_step_id("phyloai --overwrite --quiet pretree filter taper --cutoff 0.1") == "pretree.filter.taper"


class TestStepOrder:
    def test_step_order_is_list(self):
        assert isinstance(STEP_ORDER, list)
        assert len(STEP_ORDER) > 10

    def test_common_steps_in_order(self):
        assert "pretree.convert" in STEP_ORDER
        assert "pretree.align" in STEP_ORDER
        assert "tree.ml.iqtree" in STEP_ORDER
        assert "posttree.topology" in STEP_ORDER
        # convert comes before align
        assert STEP_ORDER.index("pretree.convert") < STEP_ORDER.index("pretree.align")
        # align comes before tree inference
        assert STEP_ORDER.index("pretree.align") < STEP_ORDER.index("tree.ml.iqtree")


class TestDiscoverStepsModule:
    """Module mode: scan dirs for result.json."""

    def test_single_step_module(self, tmp_path):
        """One subdir with result.json, no top-level result.json."""
        sub = tmp_path / "2-align"
        sub.mkdir()
        (sub / "result.json").write_text(json.dumps({
            "status": "success",
            "command": "phyloai pretree align --seq-dir ./raw --method linsi --threads 8",
            "wall_time": 31.4,
            "tool_versions": {"mafft": "7.526"},
            "params": {"seq_dir": "./raw", "method": "linsi", "threads": 8},
            "key_results": {"n_aligned": 100},
            "error": None,
            "data": {"output_files": {}},
        }))

        result = discover_steps(sub)
        assert result["run_mode"] == "module"
        assert len(result["steps"]) == 1
        assert result["steps"][0]["step_id"] == "pretree.align"

    def test_multi_step_module(self, tmp_path):
        """Multiple subdirs with result.json."""
        for name in ("2-align", "4-trim", "5-metrics"):
            d = tmp_path / name
            d.mkdir()
            (d / "result.json").write_text(json.dumps({
                "status": "success",
                "command": f"phyloai pretree {name.split('-', 1)[1]} --some-flag",
                "wall_time": 10.0,
                "tool_versions": {},
                "params": {},
                "key_results": {},
                "error": None,
                "data": {"output_files": {}},
            }))

        result = discover_steps(tmp_path)
        assert result["run_mode"] == "module"
        assert len(result["steps"]) == 3

    def test_no_result_json_error(self, tmp_path):
        """Empty directory should raise ValueError."""
        try:
            discover_steps(tmp_path)
            assert False, "should have raised"
        except ValueError as e:
            assert "result.json" in str(e).lower()


class TestDiscoverStepsPipeline:
    """Pipeline mode: top-level result.json + per-step subdirectories (detected by filesystem scan)."""

    def test_pipeline_detection(self, tmp_path):
        """Top-level result.json + subdirs with result.json = pipeline."""
        # Top-level run result.json (mode/speed metadata only — steps come from filesystem)
        (tmp_path / "result.json").write_text(json.dumps({
            "status": "success",
            "command": "phyloai run --seq-dir ./raw --mode supermatrix",
            "wall_time": 8040.5,
            "tool_versions": {},
            "params": {"mode": "supermatrix", "speed": "normal"},
            "key_results": {"n_input_genes": 200},
            "error": None,
            "data": {"mode": "supermatrix", "speed": "normal"},
        }))

        # Step subdirectories with their own result.json (detected by filesystem scan)
        step_dir = tmp_path / "1-convert"
        step_dir.mkdir()
        (step_dir / "result.json").write_text(json.dumps({
            "status": "success",
            "command": "phyloai pretree convert --input ./raw --to fasta",
            "wall_time": 1.0,
            "tool_versions": {},
            "params": {"to": "fasta"},
            "key_results": {},
            "error": None,
            "data": {"output_files": {}},
        }))

        step_dir2 = tmp_path / "2-align"
        step_dir2.mkdir()
        (step_dir2 / "result.json").write_text(json.dumps({
            "status": "success",
            "command": "phyloai pretree align --seq-dir ./raw --method linsi",
            "wall_time": 31.4,
            "tool_versions": {"mafft": "7.526"},
            "params": {"method": "linsi"},
            "key_results": {"n_aligned": 100},
            "error": None,
            "data": {"output_files": {}},
        }))

        result = discover_steps(tmp_path)
        assert result["run_mode"] == "pipeline"
        assert result["pipeline_summary"] is not None
        assert result["pipeline_summary"]["mode"] == "supermatrix"
        assert len(result["steps"]) == 2  # detected by filesystem scan
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `python -m pytest tests/report/test_collector.py -v`
Expected: FAIL (module doesn't exist yet)

- [ ] **Step 3: Implement `phyloai/report/collector.py`**

```python
"""Step discovery, ordering, and run-mode detection for phyloai report."""

from __future__ import annotations

import json
import re
import shlex
from pathlib import Path
from typing import Any

STEP_ORDER: list[str] = [
    "pretree.convert",
    "pretree.stats",
    "pretree.align",
    "pretree.trim",
    "pretree.metrics",
    "pretree.filter.taper",
    "pretree.filter.treeshrink",
    "pretree.filter.symtest",
    "pretree.filter.metrics",
    "pretree.filter.cluster",
    "pretree.concat",
    "tree.ml.fasttree",
    "tree.ml.iqtree",
    "tree.msc",
    "tree.bi",
    "tree.cf",
    "posttree.topology",
    "posttree.dating.hessian",
    "posttree.dating.mcmc",
    "posttree.signal",
    "posttree.syserror.brlen",
    "posttree.syserror.cca",
    "posttree.syserror.sites",
    "posttree.simulate",
]

def parse_step_id(command: str) -> str:
    """Extract step_id from a full CLI command string.

    Drops all flag tokens (any token starting with '-'), then finds the
    first known root command (pretree/tree/posttree/run) in the remaining
    tokens and reads subsequent tokens as the subcommand path.  Flag VALUES
    (e.g. ``./runs`` from ``--run-dir ./runs``) are harmless noise that
    don't match any known root and are naturally skipped by the first-match
    lookup.

    Boolean flags before the root are handled correctly because only the
    flag token ``--quiet`` is dropped; the next token ``pretree`` is not
    consumed as a flag value.

    >>> parse_step_id("phyloai --quiet pretree align --method linsi")
    'pretree.align'
    >>> parse_step_id("phyloai --run-dir ./runs pretree align --method linsi")
    'pretree.align'
    >>> parse_step_id("phyloai pretree align --seq-dir ./raw --method linsi")
    'pretree.align'
    >>> parse_step_id("phyloai run --mode supermatrix")
    'run'
    """
    # Known CLI root commands and their immediate subcommands
    _ROOT_SUBCOMMANDS: dict[str, set[str] | None] = {
        "pretree": {"convert", "stats", "align", "trim", "metrics", "filter", "concat"},
        "tree":    {"ml", "bi", "msc", "cf"},
        "posttree":{"topology", "dating", "signal", "syserror", "simulate"},
        "run":     None,   # single-token root — no subcommand
        "doctor":  None,
    }
    # Known third-level subcommands under certain parents
    _THIRD_LEVEL: dict[str, set[str]] = {
        "filter": {"taper", "treeshrink", "symtest", "metrics", "cluster"},
        "ml":     {"fasttree", "iqtree"},
        "dating": {"hessian", "mcmc"},
        "syserror": {"brlen", "cca", "sites"},
    }

    try:
        tokens = shlex.split(command)
    except ValueError:
        tokens = command.split()

    if not tokens or tokens[0] != "phyloai":
        return "unknown"

    # Drop all flag tokens (anything starting with '-').  Flag VALUES
    # remain but don't match any known root — they're harmless noise.
    parts: list[str] = [t for t in tokens[1:] if not t.startswith("-")]

    if not parts:
        return "unknown"

    # Find the first known root.  Tokens before it (e.g. "./runs" from
    # "--run-dir ./runs") are skipped.
    root_index = -1
    root = ""
    for idx, p in enumerate(parts):
        if p in _ROOT_SUBCOMMANDS:
            root_index = idx
            root = p
            break

    if root_index < 0:
        return "unknown"

    subcmds = _ROOT_SUBCOMMANDS[root]
    if subcmds is None:
        return root  # single-token root (run, doctor)

    # Read known subcommand tokens after the root
    after = parts[root_index + 1:]
    if not after:
        return root

    # Level 2
    l2 = after[0]
    if l2 not in subcmds:
        return f"{root}.{l2}" if l2 else root
    if len(after) < 2:
        return f"{root}.{l2}"

    # Level 3
    l3 = after[1]
    third = _THIRD_LEVEL.get(l2, set())
    if l3 in third:
        return f"{root}.{l2}.{l3}"
    return f"{root}.{l2}"


def _load_result_json(path: Path) -> dict[str, Any]:
    """Load and return parsed result.json, or error dict on failure."""
    try:
        with open(path) as fh:
            return json.load(fh)
    except (json.JSONDecodeError, OSError) as exc:
        return {"status": "error", "error": str(exc), "command": f"<unreadable: {path}>"}


_EXCLUDE_DIRS = {"report", "logs"}


def _find_result_json_files(base: Path, max_depth: int = 2) -> list[Path]:
    """Find result.json files up to max_depth levels below base.

    Skips directories named 'report' and 'logs' to avoid self-referential
    scanning and tool log directories.

    Depth=2 covers both current flat layouts (``1-convert/result.json`` at
    depth 1) and the spec's anticipated nested module layout
    (``pretree/2-align/result.json`` at depth 2)."""
    found: list[Path] = []
    todo: list[tuple[Path, int]] = [(base, 0)]
    while todo:
        current, depth = todo.pop()
        if depth >= max_depth:
            continue
        try:
            entries = sorted(current.iterdir())
        except PermissionError:
            continue
        for entry in entries:
            if not entry.is_dir():
                continue
            if entry.name in _EXCLUDE_DIRS or entry.name.startswith("."):
                continue
            candidate = entry / "result.json"
            if candidate.exists():
                found.append(candidate)
            todo.append((entry, depth + 1))
    return found


def _detect_run_mode(run_dir: Path) -> tuple[str, list[Path], dict[str, Any] | None]:
    """Detect run_mode by filesystem structure only (Spec Section 5 priority rules).

    Priority:
    1. Top-level result.json + subdirs with result.json → pipeline
    2. Top-level result.json only → module (single command)
    3. No top-level, but subdirs → module (multi-step)
    4. Nothing found → error

    Pipeline step paths come purely from filesystem scan — NOT from
    top-level result.json:data.steps[], which is an implementation detail
    of phyloai run that may change.  The top-level result.json is only
    read for optional metadata (mode, speed) enrichment.
    """
    top_result = run_dir / "result.json"
    sub_results = _find_result_json_files(run_dir)

    has_top = top_result.exists()

    if has_top and sub_results:
        # Priority 1: pipeline — use filesystem scan for step paths
        top_data = _load_result_json(top_result)
        pipeline_summary = {
            "mode": top_data.get("data", {}).get("mode", "unknown"),
            "speed": top_data.get("data", {}).get("speed", "normal"),
        }
        return ("pipeline", sub_results, pipeline_summary)

    if has_top and not sub_results:
        # Priority 2: single command module
        return ("module", [top_result], None)

    if not has_top and sub_results:
        # Priority 3: multi-step module
        return ("module", sub_results, None)

    # Priority 4: nothing found
    raise ValueError(
        f"No result.json found in {run_dir} or its subdirectories. "
        f"Ensure the directory contains PhyloAI run output."
    )


def _order_steps(steps: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Sort steps by STEP_ORDER, appending unknowns at end."""
    order_map = {sid: i for i, sid in enumerate(STEP_ORDER)}

    def sort_key(step: dict[str, Any]) -> tuple[int, str]:
        sid = step.get("step_id", "unknown")
        return (order_map.get(sid, 9999), sid)

    return sorted(steps, key=sort_key)


def discover_steps(run_dir: str | Path) -> dict[str, Any]:
    """Scan run_dir and return structured step data for report assembly.

    Returns:
        {
            "run_mode": "pipeline" | "module",
            "steps": [{step_id, command, status, wall_time, tool_versions,
                       params, key_results, error, data, path}, ...],
            "pipeline_summary": {mode, speed} | None,
        }
    """
    run_dir = Path(run_dir).resolve()

    run_mode, result_paths, pipeline_summary = _detect_run_mode(run_dir)

    steps: list[dict[str, Any]] = []
    for rp in result_paths:
        result = _load_result_json(rp)
        step_id = parse_step_id(result.get("command", ""))
        steps.append({
            "step_id": step_id,
            "path": str(rp.absolute()),
            "command": result.get("command", ""),
            "status": result.get("status", "error"),
            "wall_time": result.get("wall_time", 0.0),
            "tool_versions": result.get("tool_versions", {}),
            "params": result.get("params", {}),
            "key_results": result.get("key_results", {}),
            "error": result.get("error"),
            "data": result.get("data", {}),
        })

    steps = _order_steps(steps)

    return {
        "run_mode": run_mode,
        "steps": steps,
        "pipeline_summary": pipeline_summary,
    }
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `python -m pytest tests/report/test_collector.py -v`
Expected: all tests PASS

---

## Task 3: Schema — ReportRecord and report.json Assembly

**Files:**
- Create: `phyloai/report/schema.py`
- Create: `tests/report/test_schema.py`

**Interfaces:**
- Produces: `ReportStep` dataclass — holds per-step data for report.json
- Produces: `ReportRecord` dataclass — top-level report structure
- Produces: `assemble_report(discovered: dict, run_dir: Path, methods_texts: dict[str, str]) -> dict` — builds the full report.json dict
- Produces: `build_figures_index(steps: list[dict]) -> list[dict]` — extracts PDF/PNG entries from output_files
- Produces: `build_tables_index(steps: list[dict]) -> list[dict]` — extracts CSV/TSV entries from output_files

### report.json schema (Spec Section 8)

```json
{
  "phyloai_version": "0.1.0",
  "generated_at": "2026-06-27T14:23:00Z",
  "run_dir": "/abs/path/to/runs/run/faa",
  "run_mode": "pipeline",
  "status": "complete",
  "pipeline_summary": { ... },
  "steps": [ ... ],
  "methods_paragraph": "...",
  "figures_index": [ ... ],
  "tables_index": [ ... ]
}
```

**Figure/Table numbering convention (Spec Section 11):**

| Phase | Analytical group | Prefix |
|-------|-----------------|--------|
| 3 | pretree.* | Fig-3.x / Table-3.x |
| 4 | tree.* | Fig-4.x / Table-4.x |
| 5 | posttree.* | Fig-5.x / Table-5.x |

Phase numbers 1 and 2 are unused. Within each phase, numbering is sequential in STEP_ORDER order.

- [ ] **Step 1: Write failing tests for schema**

In `tests/report/test_schema.py`:

```python
"""Tests for phyloai.report.schema."""
from __future__ import annotations

import json
from datetime import datetime
from pathlib import Path

from phyloai.report.schema import (
    ReportRecord,
    ReportStep,
    assemble_report,
    build_figures_index,
    build_tables_index,
)


class TestReportStep:
    def test_create_minimal(self):
        step = ReportStep(
            step_id="pretree.align",
            command="phyloai pretree align --seq-dir ./raw",
            status="success",
            wall_time=31.4,
            tool_versions={"mafft": "7.526"},
            params={"method": "linsi", "threads": 8},
            key_results={"n_aligned": 100},
            methods_text="Multiple sequence alignments were performed...",
        )
        assert step.step_id == "pretree.align"
        assert step.status == "success"
        assert step.error is None

    def test_failed_step(self):
        step = ReportStep(
            step_id="pretree.align",
            command="phyloai pretree align ...",
            status="error",
            wall_time=0.5,
            tool_versions={},
            params={},
            key_results={},
            methods_text="",
            error="MAFFT not found",
        )
        assert step.status == "error"
        assert step.error == "MAFFT not found"
        assert step.methods_text == ""

    def test_to_dict(self):
        step = ReportStep(
            step_id="pretree.align",
            command="phyloai pretree align --seq-dir ./raw --method linsi",
            status="success",
            wall_time=31.4,
            tool_versions={"mafft": "7.526"},
            params={"method": "linsi", "seq_dir": "./raw", "threads": 8},
            key_results={"n_aligned": 100, "n_skipped": 0},
            methods_text="Alignments were performed...",
            output_files={"aligned_sequences": {"path": "/tmp/out.fa", "description": "Aligned sequences"}},
            warnings=[],
        )
        d = step.to_dict()
        assert d["step_id"] == "pretree.align"
        assert d["status"] == "success"
        assert d["methods_text"] == "Alignments were performed..."
        assert "aligned_sequences" in d["output_files"]


class TestReportRecord:
    def test_create_minimal(self):
        record = ReportRecord(
            run_dir=Path("/tmp/runs/pretree"),
            run_mode="module",
            status="complete",
        )
        assert record.run_mode == "module"
        assert record.status == "complete"
        assert record.steps == []
        assert record.figures_index == []
        assert record.tables_index == []

    def test_to_dict(self):
        record = ReportRecord(
            run_dir=Path("/tmp/runs/pretree"),
            run_mode="module",
            status="complete",
            steps=[
                ReportStep(
                    step_id="pretree.align",
                    command="phyloai pretree align ...",
                    status="success",
                    wall_time=10.0,
                    tool_versions={"mafft": "7.0"},
                    params={},
                    key_results={"n_aligned": 5},
                    methods_text="Test methods.",
                )
            ],
            methods_paragraph="Test methods.",
            pipeline_summary=None,
        )
        d = record.to_dict()
        assert d["run_mode"] == "module"
        assert d["status"] == "complete"
        assert len(d["steps"]) == 1
        assert d["methods_paragraph"] == "Test methods."
        assert "phyloai_version" in d
        assert "generated_at" in d


class TestBuildFiguresIndex:
    def test_extracts_pdf_and_png(self):
        steps = [
            {
                "step_id": "pretree.metrics",
                "output_files": {
                    "correlation_heatmap": {"path": "/tmp/corr.pdf", "description": "Correlation heatmap"},
                    "metrics_table": {"path": "/tmp/metrics.csv", "description": "Metrics table"},
                    "distribution_plot": {"path": "/tmp/dist.png", "description": "Distribution plot"},
                },
            },
        ]
        figures = build_figures_index(steps)
        assert len(figures) == 2
        paths = {f["path"] for f in figures}
        assert "/tmp/corr.pdf" in paths
        assert "/tmp/dist.png" in paths

    def test_skips_non_figure_types(self):
        steps = [
            {
                "step_id": "pretree.align",
                "output_files": {
                    "aligned": {"path": "/tmp/aligned.fa"},
                    "log": {"path": "/tmp/run.log"},
                },
            },
        ]
        figures = build_figures_index(steps)
        assert figures == []

    def test_figure_numbering_by_phase(self):
        steps = [
            {"step_id": "pretree.metrics", "output_files": {"a": {"path": "/tmp/a.pdf"}}},
            {"step_id": "tree.ml.iqtree", "output_files": {"b": {"path": "/tmp/b.pdf"}}},
            {"step_id": "posttree.topology", "output_files": {"c": {"path": "/tmp/c.pdf"}}},
        ]
        figures = build_figures_index(steps)
        assert figures[0]["figure_id"] == "Fig-3.1"
        assert figures[1]["figure_id"] == "Fig-4.1"
        assert figures[2]["figure_id"] == "Fig-5.1"

    def test_sequential_within_phase(self):
        steps = [
            {"step_id": "pretree.metrics", "output_files": {"a": {"path": "/tmp/a.pdf"}, "b": {"path": "/tmp/b.pdf"}}},
        ]
        figures = build_figures_index(steps)
        assert figures[0]["figure_id"] == "Fig-3.1"
        assert figures[1]["figure_id"] == "Fig-3.2"

    def test_missing_description_fallback(self):
        steps = [
            {"step_id": "pretree.metrics", "output_files": {"heatmap": {"path": "/tmp/hm.pdf"}}},
        ]
        figures = build_figures_index(steps)
        assert figures[0]["description"] == "heatmap"


class TestBuildTablesIndex:
    def test_extracts_csv_and_tsv(self):
        steps = [
            {
                "step_id": "pretree.metrics",
                "output_files": {
                    "metrics_table": {"path": "/tmp/metrics.csv", "description": "Metrics"},
                    "results": {"path": "/tmp/results.tsv", "description": "Results"},
                    "plot": {"path": "/tmp/plot.pdf"},
                },
            },
        ]
        tables = build_tables_index(steps)
        assert len(tables) == 2
        paths = {t["path"] for t in tables}
        assert "/tmp/metrics.csv" in paths
        assert "/tmp/results.tsv" in paths


class TestAssembleReport:
    def test_complete_run(self, tmp_path):
        discovered = {
            "run_mode": "module",
            "steps": [
                {
                    "step_id": "pretree.align",
                    "command": "phyloai pretree align --seq-dir ./raw --method linsi",
                    "status": "success",
                    "wall_time": 31.4,
                    "tool_versions": {"mafft": "7.526"},
                    "params": {"method": "linsi", "seq_dir": "./raw", "threads": 4},
                    "key_results": {"n_aligned": 100, "n_skipped": 0},
                    "error": None,
                    "data": {"output_files": {"aligned": {"path": str(tmp_path / "aligned.fa")}}},
                },
            ],
            "pipeline_summary": None,
        }
        methods_texts = {"pretree.align": "MSA was performed using MAFFT v7.526 with L-INS-i."}

        report = assemble_report(discovered, tmp_path, methods_texts)
        assert report["status"] == "complete"
        assert report["run_mode"] == "module"
        assert len(report["steps"]) == 1
        assert report["steps"][0]["methods_text"] == "MSA was performed using MAFFT v7.526 with L-INS-i."
        assert "methods_paragraph" in report
        assert "figures_index" in report
        assert "tables_index" in report

    def test_partial_failure(self, tmp_path):
        discovered = {
            "run_mode": "module",
            "steps": [
                {
                    "step_id": "pretree.align",
                    "command": "phyloai pretree align ...",
                    "status": "success",
                    "wall_time": 10.0,
                    "tool_versions": {},
                    "params": {},
                    "key_results": {},
                    "error": None,
                    "data": {"output_files": {}},
                },
                {
                    "step_id": "pretree.trim",
                    "command": "phyloai pretree trim ...",
                    "status": "error",
                    "wall_time": 0.1,
                    "tool_versions": {},
                    "params": {},
                    "key_results": {},
                    "error": "trimAl exited with code 1",
                    "data": {},
                },
            ],
            "pipeline_summary": None,
        }
        methods_texts = {"pretree.align": "Methods for align."}

        report = assemble_report(discovered, tmp_path, methods_texts)
        assert report["status"] == "partial"
        assert report["pipeline_summary"]["n_steps_success"] == 1
        assert report["pipeline_summary"]["n_steps_failed"] == 1

    def test_methods_paragraph_excludes_failed(self, tmp_path):
        discovered = {
            "run_mode": "module",
            "steps": [
                {
                    "step_id": "pretree.align",
                    "command": "phyloai pretree align ...",
                    "status": "success",
                    "wall_time": 10.0,
                    "tool_versions": {},
                    "params": {},
                    "key_results": {},
                    "error": None,
                    "data": {"output_files": {}},
                },
                {
                    "step_id": "pretree.trim",
                    "command": "phyloai pretree trim ...",
                    "status": "error",
                    "wall_time": 0.1,
                    "tool_versions": {},
                    "params": {},
                    "key_results": {},
                    "error": "trimAl failed",
                    "data": {},
                },
            ],
            "pipeline_summary": None,
        }
        methods_texts = {"pretree.align": "Align methods.", "pretree.trim": "Trim methods."}

        report = assemble_report(discovered, tmp_path, methods_texts)
        # failed step excluded from paragraph
        assert "Trim methods" not in report["methods_paragraph"]
        assert "Align methods" in report["methods_paragraph"]
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `python -m pytest tests/report/test_schema.py -v`
Expected: FAIL

- [ ] **Step 3: Implement `phyloai/report/schema.py`**

```python
"""Report data structures and report.json assembly."""

from __future__ import annotations

from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from phyloai import __version__


# Analytical phase prefixes for figure/table numbering (Spec Section 11)
_PHASE_PREFIX: dict[str, int] = {
    "pretree": 3,
    "tree": 4,
    "posttree": 5,
}

_FIGURE_EXTENSIONS = {".pdf", ".png"}
_TABLE_EXTENSIONS = {".csv", ".tsv"}


def _get_phase_prefix(step_id: str) -> int:
    """Map step_id to analytical phase number (3=pretree, 4=tree, 5=posttree)."""
    top = step_id.split(".")[0]
    return _PHASE_PREFIX.get(top, 99)


@dataclass
class ReportStep:
    """A single analytical step in the report."""
    step_id: str
    command: str
    status: str
    wall_time: float
    tool_versions: dict[str, str]
    params: dict[str, Any]
    key_results: dict[str, Any]
    methods_text: str = ""
    output_files: dict[str, dict[str, str]] = field(default_factory=dict)
    warnings: list[str] = field(default_factory=list)
    error: str | None = None

    def to_dict(self) -> dict[str, Any]:
        return {
            "step_id": self.step_id,
            "command": self.command,
            "status": self.status,
            "wall_time": self.wall_time,
            "tool_versions": self.tool_versions,
            "params": self.params,
            "key_results": self.key_results,
            "methods_text": self.methods_text,
            "output_files": self.output_files,
            "warnings": self.warnings,
            "error": self.error,
        }


@dataclass
class ReportRecord:
    """Top-level report structure."""
    run_dir: Path
    run_mode: str  # "pipeline" | "module"
    status: str  # "complete" | "partial" | "failed"
    phyloai_version: str = __version__
    generated_at: str = field(default_factory=lambda: datetime.now(timezone.utc).isoformat())
    steps: list[ReportStep] = field(default_factory=list)
    methods_paragraph: str = ""
    pipeline_summary: dict[str, Any] | None = None
    figures_index: list[dict[str, Any]] = field(default_factory=list)
    tables_index: list[dict[str, Any]] = field(default_factory=list)

    def to_dict(self) -> dict[str, Any]:
        return {
            "phyloai_version": self.phyloai_version,
            "generated_at": self.generated_at,
            "run_dir": str(self.run_dir.absolute()),
            "run_mode": self.run_mode,
            "status": self.status,
            "pipeline_summary": self.pipeline_summary or {},
            "steps": [s.to_dict() for s in self.steps],
            "methods_paragraph": self.methods_paragraph,
            "figures_index": self.figures_index,
            "tables_index": self.tables_index,
        }


def build_figures_index(steps: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Extract PDF/PNG entries from step output_files into a flat figure index.

    Numbering: Fig-{phase}.{seq} where phase = 3 for pretree, 4 for tree, 5 for posttree.
    """
    figures: list[dict[str, Any]] = []
    # Track counters per phase
    counters: dict[int, int] = {}
    for step in steps:
        step_id = step.get("step_id", "unknown")
        phase = _get_phase_prefix(step_id)
        output_files = step.get("output_files") or step.get("data", {}).get("output_files", {})
        for label, file_obj in output_files.items():
            if not isinstance(file_obj, dict):
                continue
            path = file_obj.get("path", "")
            ext = Path(path).suffix.lower()
            if ext in _FIGURE_EXTENSIONS:
                counters[phase] = counters.get(phase, 0) + 1
                figures.append({
                    "figure_id": f"Fig-{phase}.{counters[phase]}",
                    "step_id": step_id,
                    "label": label,
                    "caption": file_obj.get("description", label),
                    "description": file_obj.get("description", label),
                    "path": path,
                    "type": ext.lstrip("."),
                })
    return figures


def build_tables_index(steps: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Extract CSV/TSV entries from step output_files into a flat table index."""
    tables: list[dict[str, Any]] = []
    counters: dict[int, int] = {}
    for step in steps:
        step_id = step.get("step_id", "unknown")
        phase = _get_phase_prefix(step_id)
        output_files = step.get("output_files") or step.get("data", {}).get("output_files", {})
        for label, file_obj in output_files.items():
            if not isinstance(file_obj, dict):
                continue
            path = file_obj.get("path", "")
            ext = Path(path).suffix.lower()
            if ext in _TABLE_EXTENSIONS:
                counters[phase] = counters.get(phase, 0) + 1
                tables.append({
                    "table_id": f"Table-{phase}.{counters[phase]}",
                    "step_id": step_id,
                    "label": label,
                    "caption": file_obj.get("description", label),
                    "description": file_obj.get("description", label),
                    "path": path,
                    "type": ext.lstrip("."),
                })
    return tables


def assemble_report(
    discovered: dict[str, Any],
    run_dir: Path,
    methods_texts: dict[str, str],
) -> dict[str, Any]:
    """Assemble the full report.json dict from discovered steps and methods texts.

    Args:
        discovered: Output from collector.discover_steps()
        run_dir: Absolute path to the run directory
        methods_texts: Dict mapping step_id -> methods text string

    Returns:
        Complete report.json as a dict, ready for JSON serialization.
    """
    steps: list[dict[str, Any]] = []
    n_success = 0
    n_failed = 0
    n_skipped = 0
    total_wall_time = 0.0

    for raw_step in discovered["steps"]:
        step_id = raw_step["step_id"]
        status = raw_step["status"]

        if status == "success":
            n_success += 1
        elif status == "error":
            n_failed += 1
        total_wall_time += raw_step.get("wall_time", 0.0)

        # Extract output_files from data
        output_files = raw_step.get("data", {}).get("output_files", {})

        # Failed steps get empty methods_text
        methods_text = ""
        if status == "success":
            methods_text = methods_texts.get(step_id, "")

        step_record = {
            "step_id": step_id,
            "command": raw_step.get("command", ""),
            "status": status,
            "wall_time": raw_step.get("wall_time", 0.0),
            "tool_versions": raw_step.get("tool_versions", {}),
            "params": raw_step.get("params", {}),
            "key_results": raw_step.get("key_results", {}),
            "methods_text": methods_text,
            "output_files": output_files,
            "warnings": raw_step.get("warnings", []),
            "error": raw_step.get("error"),
        }
        steps.append(step_record)

    # Build figures and tables indexes from all steps
    figures_index = build_figures_index(steps)
    tables_index = build_tables_index(steps)

    # Build methods_paragraph from successful steps in order
    method_parts: list[str] = []
    for s in steps:
        if s["status"] == "success" and s["methods_text"]:
            method_parts.append(s["methods_text"])
    methods_paragraph = " ".join(method_parts)

    # Determine overall status
    if n_failed == 0 and n_success > 0:
        overall_status = "complete"
    elif n_success == 0 and n_failed > 0:
        overall_status = "failed"
    else:
        overall_status = "partial"

    # Pipeline summary
    pipeline_summary = {
        "n_steps_total": len(steps),
        "n_steps_success": n_success,
        "n_steps_failed": n_failed,
        "n_steps_skipped": n_skipped,
        "total_wall_time": total_wall_time,
    }
    # Merge pipeline-specific data from discovered
    if discovered.get("pipeline_summary"):
        pipeline_summary.update(discovered["pipeline_summary"])

    report = {
        "phyloai_version": __version__,
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "run_dir": str(run_dir.absolute()),
        "run_mode": discovered["run_mode"],
        "status": overall_status,
        "pipeline_summary": pipeline_summary,
        "steps": steps,
        "methods_paragraph": methods_paragraph,
        "figures_index": figures_index,
        "tables_index": tables_index,
    }

    return report
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `python -m pytest tests/report/test_schema.py -v`
Expected: all tests PASS

---

## Task 4: Templates — Methods Text Generation

**Files:**
- Create: `phyloai/report/templates.py`
- Create: `tests/report/test_templates.py`

**Interfaces:**
- Produces: `generate_methods_<step_id>(params: dict, key_results: dict, tool_versions: dict) -> str` — one function per step_id
- Produces: `METHODS_GENERATORS: dict[str, callable]` — registry mapping step_id → generator function
- Produces: `generate_all_methods(step_id: str, params: dict, key_results: dict, tool_versions: dict) -> str` — dispatcher that returns "" for unknown or failed steps

### Template functions required (Spec Section 9.1):

| step_id | Template |
|---------|----------|
| `pretree.convert` | Format conversion summary |
| `pretree.stats` | Sequence statistics summary |
| `pretree.align` | MAFFT/MAGUS alignment (method_description, backtrans_clause) |
| `pretree.trim` | trimAl/BMGE/ClipKIT trimming |
| `pretree.metrics` | Phylogenetic informativeness metrics |
| `pretree.filter.taper` | TAPER error-site masking |
| `pretree.filter.treeshrink` | TreeShrink taxon pruning |
| `pretree.filter.symtest` | Symmetry test filtering |
| `pretree.filter.metrics` | Metric-rule filtering |
| `pretree.filter.cluster` | UMAP cluster filtering |
| `pretree.concat` | Supermatrix concatenation (recoding_clause, outgroup_clause) |
| `tree.ml.fasttree` | FastTree ML inference |
| `tree.ml.iqtree` | IQ-TREE ML inference (partition_clause, modelfinder_clause, model_result_clause) |
| `tree.msc` | wASTRAL multispecies coalescent |
| `tree.bi` | PhyloBayes Bayesian inference |
| `tree.cf` | Concordance factors |
| `posttree.topology` | AU/WKH/WSH topology tests |
| `posttree.dating.hessian` | Hessian computation |
| `posttree.dating.mcmc` | MCMCTree divergence dating (model_descr, clock_descr, diag_descr) |
| `posttree.signal` | Phylogenetic signal analysis |
| `posttree.syserror.brlen` | Branch length systematic error |
| `posttree.syserror.cca` | CCA systematic error |
| `posttree.syserror.sites` | Site-wise systematic error |
| `posttree.simulate` | AliSim simulation |

**Design principles for templates:**
- Each function receives the complete `params` dict but only reads scientifically meaningful keys
- Technical params (threads, paths, --quiet, --overwrite, --resume, --dry-run) ignored
- All scientifically meaningful parameters described, whether or not they differ from defaults
- Conditional branches for parameter combinations
- All placeholder values have fallbacks using `.get()`

- [ ] **Step 1: Write test for templates registry**

In `tests/report/test_templates.py`:

```python
"""Tests for phyloai.report.templates."""
from __future__ import annotations

from phyloai.report.templates import (
    METHODS_GENERATORS,
    generate_all_methods,
    generate_methods_pretree_align,
    generate_methods_pretree_convert,
)


class TestTemplatesRegistry:
    def test_has_all_implemented_steps(self):
        """Every step_id in STEP_ORDER MUST have a template generator.
        Missing entries fail fast so new commands aren't silently omitted
        from methods_paragraph.
        """
        from phyloai.report.collector import STEP_ORDER

        for step_id in STEP_ORDER:
            assert step_id in METHODS_GENERATORS, (
                f"Missing methods template for {step_id}. "
                f"Add generate_methods_{step_id.replace('.', '_')}() to templates.py."
            )
            gen = METHODS_GENERATORS[step_id]
            assert callable(gen), f"{step_id} generator is not callable"


class TestGenerateAllMethods:
    def test_dispatches(self):
        result = generate_all_methods(
            "pretree.convert",
            params={"to": "fasta", "input": "./raw"},
            key_results={"n_converted": 100},
            tool_versions={},
        )
        assert isinstance(result, str)
        assert len(result) > 0

    def test_unknown_step_returns_empty(self):
        result = generate_all_methods(
            "nonexistent.step",
            params={},
            key_results={},
            tool_versions={},
        )
        assert result == ""

    def test_failed_step_returns_empty(self):
        result = generate_all_methods(
            "pretree.align",
            params={},
            key_results={},
            tool_versions={},
            status="error",
        )
        assert result == ""


class TestGenerateMethodsPretreeAlign:
    def test_linsi(self):
        text = generate_methods_pretree_align(
            params={"method": "linsi", "seq_type": "AA", "backtrans": False},
            key_results={"n_aligned": 100, "n_skipped": 0, "mean_alignment_length": 500.0},
            tool_versions={"mafft": "7.526"},
        )
        assert "L-INS-i" in text
        assert "MAFFT" in text
        assert "7.526" in text

    def test_magus(self):
        text = generate_methods_pretree_align(
            params={"method": "magus", "seq_type": "AA", "backtrans": False},
            key_results={"n_aligned": 200, "n_skipped": 0, "mean_alignment_length": 300.0},
            tool_versions={"magus": "unknown version"},
        )
        assert "MAGUS" in text or "magus" in text.lower()

    def test_backtrans(self):
        text = generate_methods_pretree_align(
            params={"method": "linsi", "seq_type": "AA", "backtrans": True},
            key_results={"n_aligned": 100, "n_skipped": 0, "mean_alignment_length": 500.0},
            tool_versions={"mafft": "7.526", "trimal": "1.4.1"},
        )
        assert "back-translation" in text.lower() or "codon-aware" in text.lower()

    def test_skipped_clause(self):
        text = generate_methods_pretree_align(
            params={"method": "auto", "seq_type": "NT", "backtrans": False},
            key_results={"n_aligned": 90, "n_skipped": 10, "mean_alignment_length": 200.0},
            tool_versions={"mafft": "7.526"},
        )
        assert "10" in text  # skipped count should appear


class TestGenerateMethodsPretreeFilterTaper:
    def test_basic(self):
        text = generate_all_methods(
            "pretree.filter.taper",
            params={"cutoff": 0.1},
            key_results={"n_input": 100, "n_retained": 95, "n_dropped": 5, "n_masked_sites": 1200},
            tool_versions={"taper": "1.0.0", "julia": "1.9"},
        )
        assert "TAPER" in text
        assert "1.0.0" in text
        # Check that dropped count appears
        assert "5" in text


class TestGenerateMethodsPretreeConcat:
    def test_basic(self):
        text = generate_all_methods(
            "pretree.concat",
            params={"taxa_occupancy": 0.75, "seq_type": "AA", "recoding": None},
            key_results={
                "n_msa_used": 187,
                "n_msa_dropped": 13,
                "n_taxa": 52,
                "total_length": 45000,
                "gap_ratio": 0.15,
                "pi_ratio": 0.42,
            },
            tool_versions={},
        )
        assert "75%" in text or "0.75" in text
        assert "187" in text
        assert "52" in text

    def test_recoding(self):
        text = generate_all_methods(
            "pretree.concat",
            params={"taxa_occupancy": 0.75, "seq_type": "AA", "recoding": "Dayhoff6"},
            key_results={
                "n_msa_used": 100,
                "n_msa_dropped": 0,
                "n_taxa": 50,
                "total_length": 30000,
                "gap_ratio": 0.1,
                "pi_ratio": 0.4,
            },
            tool_versions={},
        )
        assert "Dayhoff6" in text
        assert "recod" in text.lower()


class TestGenerateMethodsTreeMlIqtree:
    def test_unpartitioned(self):
        text = generate_all_methods(
            "tree.ml.iqtree",
            params={
                "partitioned": False,
                "modelfinder": "MFP",
                "mset": "LG+C20+C60",
                "boot": 1000,
            },
            key_results={"log_likelihood": -12345.67},
            tool_versions={"iqtree": "3.0.0"},
        )
        assert "IQ-TREE" in text
        assert "3.0.0" in text
        assert "1,000" in text or "1000" in text

    def test_partitioned(self):
        text = generate_all_methods(
            "tree.ml.iqtree",
            params={
                "partitioned": True,
                "merged_partitions": True,
                "rclusterf": 10,
                "modelfinder": "MFP",
                "mset": "LG+C20+C60",
                "boot": 1000,
            },
            key_results={"log_likelihood": -5000.0},
            tool_versions={"iqtree": "3.0.0"},
        )
        assert "partition" in text.lower()
        assert "rclusterf" in text.lower()
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `python -m pytest tests/report/test_templates.py -v`
Expected: FAIL (templates.py not yet implemented)

- [ ] **Step 3: Implement `phyloai/report/templates.py`**

```python
"""Per-command methods text generators for phyloai report.

Each step_id maps to a dedicated function that reads scientifically
meaningful parameters from params, key_results, and tool_versions,
producing 2-5 sentences of academic English suitable for journal Methods.
"""

from __future__ import annotations

from typing import Any

# ---------------------------------------------------------------------------
# Method description mappings
# ---------------------------------------------------------------------------

_ALIGN_METHOD_MAP: dict[str, tuple[str, str]] = {
    "linsi": (
        "L-INS-i",
        "applies iterative local pairwise alignment refinement and is suited "
        "for sequences with conserved domains and insertions",
    ),
    "einsi": (
        "E-INS-i",
        "uses multiple local alignments and is suited for sequences with "
        "multiple conserved regions separated by unalignable regions",
    ),
    "ginsi": (
        "G-INS-i",
        "applies global pairwise alignment and is suited for sequences of "
        "similar length without large insertions",
    ),
    "fftns1": (
        "FFT-NS-1",
        "uses progressive alignment with single FFT iteration and is suited "
        "for large datasets where speed is prioritized",
    ),
    "fftns2": (
        "FFT-NS-2",
        "uses progressive alignment with two FFT iterations",
    ),
    "auto": (
        "auto-selected",
        "strategy selected automatically by MAFFT based on sequence length and count",
    ),
    "magus": (
        "MAGUS",
        "uses graph-based divide-and-conquer alignment and is suited for "
        "very large or highly divergent datasets",
    ),
}

_CLOCK_MAP: dict[str, str] = {
    "strict": "strict",
    "independent": "independent-rates",
    "correlated": "autocorrelated-rates",
}


def _f(value: Any, default: str = "unknown") -> str:
    """Safe string conversion with fallback."""
    if value is None:
        return default
    return str(value)


def _describe_n(value: Any, singular: str, plural: str | None = None) -> str:
    """Describe count: '1 locus' vs '5 loci'."""
    try:
        n = int(value)
    except (TypeError, ValueError):
        return f"{value} {plural or singular + 's'}"
    if n == 1:
        return f"1 {singular}"
    return f"{n:,} {plural or singular + 's'}"


# ---------------------------------------------------------------------------
# Per-step methods generators
# ---------------------------------------------------------------------------

def generate_methods_pretree_convert(
    params: dict[str, Any],
    key_results: dict[str, Any],
    tool_versions: dict[str, Any],
) -> str:
    """Format conversion methods text."""
    to_fmt = params.get("to", "FASTA")
    n = key_results.get("n_converted", 0)
    n_failed = key_results.get("n_failed", 0)
    parts = [
        f"Raw sequence files were converted to {to_fmt.upper()} format.",
        f"A total of {_describe_n(n, 'file')} were successfully converted",
    ]
    if n_failed:
        parts.append(f" ({_describe_n(n_failed, 'file')} failed conversion).")
    else:
        parts.append(".")
    return " ".join(parts)


def generate_methods_pretree_stats(
    params: dict[str, Any],
    key_results: dict[str, Any],
    tool_versions: dict[str, Any],
) -> str:
    """Sequence/alignment statistics methods text."""
    n = key_results.get("n_files", 0)
    mean_len = key_results.get("mean_length", 0)
    format_name = params.get("input_format", "FASTA")
    parts = [
        f"Sequence statistics were computed for {_describe_n(n, 'sequence file')} "
        f"in {format_name.upper() if isinstance(format_name, str) else 'detected'} format."
    ]
    if mean_len:
        parts.append(f"Mean sequence length was {mean_len:.1f} bp.")
    if key_results.get("n_taxa_total"):
        parts.append(f"A total of {key_results['n_taxa_total']:,} unique taxa were identified.")
    return " ".join(parts)


def generate_methods_pretree_align(
    params: dict[str, Any],
    key_results: dict[str, Any],
    tool_versions: dict[str, Any],
) -> str:
    """MSA methods text (MAFFT / MAGUS)."""
    method = params.get("method", "auto")
    if isinstance(method, str):
        method = method.lower()
    else:
        method = "auto"

    desc, rationale = _ALIGN_METHOD_MAP.get(
        method, (method.upper(), "performs multiple sequence alignment")
    )
    tool_name = "MAFFT" if method != "magus" else "MAGUS"
    version = tool_versions.get("mafft" if method != "magus" else "magus", "unknown version")

    n_aligned = key_results.get("n_aligned", 0)
    n_skipped = key_results.get("n_skipped", 0)
    seq_type = params.get("seq_type", "AA")
    mean_len = key_results.get("mean_alignment_length", 0)

    lines = [
        f"Multiple sequence alignments were performed using {tool_name} v{version} "
        f"with the {desc} algorithm, which {rationale}."
    ]

    n_line = f"A total of {_describe_n(n_aligned, f'{seq_type} locus', f'{seq_type} loci')} were aligned"
    if n_skipped > 0:
        n_line += f" ({_describe_n(n_skipped, 'locus', 'loci')} skipped)."
    else:
        n_line += "."
    lines.append(n_line)

    if mean_len:
        lines.append(
            f"Mean alignment length was {mean_len:.1f} bp "
            f"across a mean of {key_results.get('mean_n_taxa', '?')} taxa per locus."
        )

    # Backtrans clause
    if params.get("backtrans"):
        trimal_ver = tool_versions.get("trimal", "unknown version")
        lines.append(
            f" Codon-aware nucleotide alignments were produced via back-translation "
            f"using trimAl v{trimal_ver}, preserving reading frame in the nucleotide alignments."
        )

    return " ".join(lines)


def generate_methods_pretree_trim(
    params: dict[str, Any],
    key_results: dict[str, Any],
    tool_versions: dict[str, Any],
) -> str:
    """Trimming methods text."""
    tool = params.get("tool", "trimal")
    tool_ver_key = {"trimal": "trimal", "bmge": "bmge", "clipkit": "clipkit"}.get(tool, tool)
    version = tool_versions.get(tool_ver_key, "unknown version")
    n = key_results.get("n_trimmed", key_results.get("n_aligned", 0))
    return (
        f"Alignments were trimmed using {tool.upper()} v{version} "
        f"to remove poorly aligned regions and reduce phylogenetic noise. "
        f"A total of {_describe_n(n, 'alignment')} were processed."
    )


def generate_methods_pretree_metrics(
    params: dict[str, Any],
    key_results: dict[str, Any],
    tool_versions: dict[str, Any],
) -> str:
    """Phylogenetic informativeness metrics methods text."""
    n_markers = key_results.get("n_markers", key_results.get("n_files", 0))
    n_metrics = key_results.get("n_metrics", 0)
    tree_dir = params.get("tree_dir")
    text = (
        f"Phylogenetic informativeness metrics were computed for "
        f"{_describe_n(n_markers, 'locus', 'loci')}"
    )
    if n_metrics:
        text += f" across {n_metrics} dimensions"
    text += (
        ". Evaluated metrics included locus length, number of informative sites, "
        "gap percentage, GC content, and RCFV (relative composition frequency variability)"
    )
    if tree_dir is not None:
        text += (
            ", as well as tree-based metrics including treeness and "
            "Robinson-Foulds distance between gene trees and the reference species tree"
        )
    text += (
        ". Pairwise Spearman correlations were computed across all metrics "
        "and visualized as a heatmap for diagnostic evaluation of metric "
        "redundancy and complementarity."
    )
    return text


def generate_methods_pretree_filter_taper(
    params: dict[str, Any],
    key_results: dict[str, Any],
    tool_versions: dict[str, Any],
) -> str:
    """TAPER error-site masking methods text."""
    taper_ver = tool_versions.get("taper", "unknown version")
    julia_ver = tool_versions.get("julia", "unknown version")
    cutoff = params.get("cutoff", 0.1)
    n_input = key_results.get("n_input", 0)
    n_retained = key_results.get("n_retained", 0)
    n_dropped = key_results.get("n_dropped", 0)
    n_masked = key_results.get("n_masked_sites", 0)

    return (
        f"Aligned sequences were screened for compositional bias and systematic "
        f"sequencing errors using TAPER v{taper_ver} (correction_multi.jl, "
        f"executed via Julia v{julia_ver}). TAPER applies a moving-window approach "
        f"to identify and mask amino acid sites within individual sequences that "
        f"deviate from expected substitution patterns, without discarding entire loci. "
        f"The masking stringency cutoff was set to {cutoff} "
        f"(`-c {cutoff}`). "
        f"Of {_describe_n(n_input, 'input locus', 'input loci')}, "
        f"{_describe_n(n_retained, 'was', 'were')} retained"
        + (f" ({_describe_n(n_dropped, 'locus', 'loci')} dropped)" if n_dropped > 0 else "")
        + (f". A total of {n_masked:,} sites were masked across all retained loci." if n_masked > 0 else ".")
    )


def generate_methods_pretree_filter_treeshrink(
    params: dict[str, Any],
    key_results: dict[str, Any],
    tool_versions: dict[str, Any],
) -> str:
    """TreeShrink taxon pruning methods text."""
    ts_ver = tool_versions.get("treeshrink", "unknown version")
    n_input = key_results.get("n_input", 0)
    n_retained = key_results.get("n_retained", 0)
    n_modified = key_results.get("n_modified", 0)
    return (
        f"Gene trees were screened for outlier long branches using "
        f"TreeShrink v{ts_ver}, which removes taxa whose removal "
        f"disproportionately reduces tree diameter. "
        f"Of {_describe_n(n_input, 'input gene tree')}, "
        f"{_describe_n(n_retained, 'was', 'were')} retained unchanged"
        + (f" and {_describe_n(n_modified, 'tree', 'trees')} were modified "
           f"by removing outlier taxa." if n_modified > 0 else ".")
    )


def generate_methods_pretree_filter_symtest(
    params: dict[str, Any],
    key_results: dict[str, Any],
    tool_versions: dict[str, Any],
) -> str:
    """Symmetry test filtering methods text."""
    n_input = key_results.get("n_input", 0)
    n_retained = key_results.get("n_retained", 0)
    n_dropped = key_results.get("n_dropped", 0)
    return (
        f"Alignments were tested for substitutional symmetry using "
        f"pairwise symmetry tests. "
        f"Of {_describe_n(n_input, 'input locus', 'input loci')}, "
        f"{_describe_n(n_retained, 'was', 'were')} retained "
        f"({_describe_n(n_dropped, 'locus', 'loci')} failed the symmetry test)."
    )


def generate_methods_pretree_filter_metrics(
    params: dict[str, Any],
    key_results: dict[str, Any],
    tool_versions: dict[str, Any],
) -> str:
    """Metric-rule filtering methods text."""
    n_input = key_results.get("n_input", 0)
    n_retained = key_results.get("n_retained", 0)
    n_dropped = key_results.get("n_dropped", 0)
    keep_rule = params.get("keep", "default criteria")
    return (
        f"Loci were filtered based on phylogenetic informativeness metrics "
        f"(keep rule: {keep_rule}). "
        f"Of {_describe_n(n_input, 'input locus', 'input loci')}, "
        f"{_describe_n(n_retained, 'was', 'were')} retained "
        f"({_describe_n(n_dropped, 'locus', 'loci')} excluded)."
    )


def generate_methods_pretree_filter_cluster(
    params: dict[str, Any],
    key_results: dict[str, Any],
    tool_versions: dict[str, Any],
) -> str:
    """UMAP cluster filtering methods text."""
    n_input = key_results.get("n_input", 0)
    n_retained = key_results.get("n_retained", 0)
    n_clusters = key_results.get("n_clusters", 0)
    text = (
        f"Loci were projected via UMAP based on phylogenetic informativeness "
        f"metrics and clustered into {n_clusters} groups."
    )
    if n_retained:
        text += (
            f" After filtering, {_describe_n(n_retained, 'locus', 'loci')} "
            f"were retained from {_describe_n(n_input, 'input locus', 'input loci')}."
        )
    return text


def generate_methods_pretree_concat(
    params: dict[str, Any],
    key_results: dict[str, Any],
    tool_versions: dict[str, Any],
) -> str:
    """Supermatrix concatenation methods text."""
    taxa_occ = params.get("taxa_occupancy", 0.75)
    seq_type = params.get("seq_type", "AA")
    n_used = key_results.get("n_msa_used", 0)
    n_dropped = key_results.get("n_msa_dropped", 0)
    n_taxa = key_results.get("n_taxa", 0)
    total_len = key_results.get("total_length", 0)
    gap_ratio = key_results.get("gap_ratio", 0)
    pi_ratio = key_results.get("pi_ratio", 0)

    text = (
        f"Trimmed alignments were concatenated into a supermatrix using "
        f"phyloai concat. Loci were included only if they met a minimum "
        f"taxon occupancy threshold of {taxa_occ:.0%} "
        f"(`--taxa-occupancy {taxa_occ}`); "
    )
    if n_dropped:
        text += f"{_describe_n(n_dropped, 'locus', 'loci')} were excluded for failing this criterion. "
    text += (
        f"The final supermatrix comprised {_describe_n(n_used, 'locus', 'loci')} "
        f"across {_describe_n(n_taxa, 'taxon', 'taxa')} with a total alignment "
        f"length of {total_len:,} {seq_type} positions "
        f"(gap ratio: {gap_ratio:.1%}; parsimony-informative sites: {pi_ratio:.1%})."
    )

    recoding = params.get("recoding")
    if recoding:
        recoding_groups = {"Dayhoff6": 6, "SR4": 4}.get(recoding, 6)
        text += (
            f" To reduce the influence of substitution saturation and "
            f"compositional heterogeneity, sequences were recoded into "
            f"{recoding} categories (`--recoding {recoding}`), collapsing "
            f"the 20 standard amino acids into {recoding_groups} biochemically "
            f"similar groups; both the original and recoded matrices were retained."
        )

    return text


def generate_methods_tree_ml_fasttree(
    params: dict[str, Any],
    key_results: dict[str, Any],
    tool_versions: dict[str, Any],
) -> str:
    """FastTree ML inference methods text."""
    version = tool_versions.get("fasttree", "unknown version")
    n_trees = key_results.get("n_trees", 1)
    model = params.get("model", "LG")
    return (
        f"Maximum likelihood phylogenetic trees were inferred using "
        f"FastTree v{version} under the {model} substitution model. "
        f"A total of {_describe_n(n_trees, 'gene tree')} were produced."
    )


def generate_methods_tree_ml_iqtree(
    params: dict[str, Any],
    key_results: dict[str, Any],
    tool_versions: dict[str, Any],
) -> str:
    """IQ-TREE ML inference methods text."""
    version = tool_versions.get("iqtree", "unknown version")
    partitioned = params.get("partitioned", False)
    merged = params.get("merged_partitions", False)
    mfinder = params.get("modelfinder", "MFP")
    mset = params.get("mset", "LG")
    boot = params.get("boot", 1000)
    log_lk = key_results.get("log_likelihood")

    text = (
        f"Maximum likelihood phylogenetic inference was performed using "
        f"IQ-TREE v{version}."
    )

    # Partition clause
    if partitioned and merged:
        rclusterf = params.get("rclusterf", 10)
        text += (
            f" A partitioned analysis was conducted with partition merging "
            f"enabled (`--merge`), using the rclusterf algorithm "
            f"(`--rclusterf {rclusterf}`) to identify the optimal merging "
            f"scheme by evaluating {rclusterf}% of candidate partition pairs."
        )
    elif partitioned:
        text += " A partitioned analysis was conducted using the provided partition scheme."

    text += (
        f" Substitution models were selected using ModelFinder ({mfinder}) "
        f"from a candidate set comprising {mset} matrix models (`--mset {mset}`)."
    )

    if log_lk is not None:
        text += f" The final log-likelihood of the best tree was {log_lk:.2f}."

    text += (
        f" Branch support was assessed using {boot:,} ultrafast bootstrap "
        f"replicates (`-B {boot}`)."
    )

    return text


def generate_methods_tree_msc(
    params: dict[str, Any],
    key_results: dict[str, Any],
    tool_versions: dict[str, Any],
) -> str:
    """wASTRAL multispecies coalescent methods text."""
    version = tool_versions.get("wastral", "unknown version")
    n_trees = key_results.get("n_gene_trees", 0)
    return (
        f"Species tree inference was performed under the multispecies "
        f"coalescent model using wASTRAL v{version}. "
        f"A total of {_describe_n(n_trees, 'gene tree')} were used as input."
    )


def generate_methods_tree_bi(
    params: dict[str, Any],
    key_results: dict[str, Any],
    tool_versions: dict[str, Any],
) -> str:
    """PhyloBayes Bayesian inference methods text."""
    version = tool_versions.get("phylobayes", "unknown version")
    chains = params.get("chains", 3)
    model = params.get("model", "CAT-GTR")
    return (
        f"Bayesian phylogenetic inference was performed using "
        f"PhyloBayes-MPI v{version} under the {model} model. "
        f"{chains} independent MCMC chains were run; convergence was assessed "
        f"using the bpcomp and tracecomp diagnostics."
    )


def generate_methods_tree_cf(
    params: dict[str, Any],
    key_results: dict[str, Any],
    tool_versions: dict[str, Any],
) -> str:
    """Concordance factors methods text."""
    version = tool_versions.get("iqtree", "unknown version")
    cf_type = params.get("cf", "gCF")
    cf_names = {"gcf": "gene concordance factor (gCF)", "scf": "site concordance factor (sCF)"}
    cf_desc = cf_names.get(cf_type.lower(), cf_type)
    return (
        f"Concordance factors were calculated using IQ-TREE v{version}. "
        f"{cf_desc} values were computed to assess phylogenetic discordance "
        f"across the dataset."
    )


def generate_methods_posttree_topology(
    params: dict[str, Any],
    key_results: dict[str, Any],
    tool_versions: dict[str, Any],
) -> str:
    """Topology test methods text."""
    version = tool_versions.get("iqtree", "unknown version")
    n_trees = key_results.get("n_candidate_trees", 0)
    return (
        f"Topology hypothesis testing was performed using IQ-TREE v{version}. "
        f"Approximately unbiased (AU), weighted Kishino-Hasegawa (WKH), and "
        f"weighted Shimodaira-Hasegawa (WSH) tests were applied to "
        f"{_describe_n(n_trees, 'candidate topology', 'candidate topologies')}."
    )


def generate_methods_posttree_dating_hessian(
    params: dict[str, Any],
    key_results: dict[str, Any],
    tool_versions: dict[str, Any],
) -> str:
    """Hessian computation methods text."""
    version = tool_versions.get("iqtree", "unknown version")
    return (
        f"The Hessian matrix required for approximate likelihood calculation "
        f"in MCMCTree was computed using IQ-TREE v{version} "
        f"(`--dating mcmctree`)."
    )


def generate_methods_posttree_dating_mcmc(
    params: dict[str, Any],
    key_results: dict[str, Any],
    tool_versions: dict[str, Any],
) -> str:
    """MCMCTree divergence dating methods text."""
    paml_ver = tool_versions.get("paml", "unknown version")
    model = params.get("model", "JC69")
    clock = params.get("clock", "independent")
    clock_desc = _CLOCK_MAP.get(clock, clock)
    n_runs = params.get("n_runs", key_results.get("n_runs", 2))
    burnin = params.get("burnin", 2000000)
    sample = params.get("sample", 20000000)
    sample_freq = params.get("sample_freq", 1000)

    return (
        f"Divergence time estimation was performed using MCMCTree "
        f"(PAML v{paml_ver}) under a {model} substitution model with a "
        f"{clock_desc} molecular clock. "
        f"{n_runs} independent MCMC chains were run for "
        f"{burnin:,} burn-in generations followed by "
        f"{sample:,} sampling generations, with parameters sampled every "
        f"{sample_freq:,} generations. "
        f"Convergence was assessed using trace inspection and effective "
        f"sample size (ESS) calculation. "
        f"Posterior node age estimates and 95% highest posterior density "
        f"(HPD) intervals are summarised in the node ages table; MCMC trace "
        f"plots for posterior and prior are provided in the supplementary figures."
    )


def generate_methods_posttree_signal(
    params: dict[str, Any],
    key_results: dict[str, Any],
    tool_versions: dict[str, Any],
) -> str:
    """Phylogenetic signal methods text."""
    n_hypotheses = key_results.get("n_hypotheses", 2)
    return (
        f"Phylogenetic signal was assessed using Four-cluster Likelihood "
        f"Mapping (FcLM) across {_describe_n(n_hypotheses, 'topological hypothesis', 'topological hypotheses')}."
    )


def generate_methods_posttree_syserror_brlen(
    params: dict[str, Any],
    key_results: dict[str, Any],
    tool_versions: dict[str, Any],
) -> str:
    """Branch length systematic error methods text."""
    return "Branch length heterogeneity and potential systematic biases were assessed."


def generate_methods_posttree_syserror_cca(
    params: dict[str, Any],
    key_results: dict[str, Any],
    tool_versions: dict[str, Any],
) -> str:
    """CCA systematic error methods text."""
    return "Cross-comparative analysis (CCA) of systematic error was performed by comparing tree topologies under different substitution models."


def generate_methods_posttree_syserror_sites(
    params: dict[str, Any],
    key_results: dict[str, Any],
    tool_versions: dict[str, Any],
) -> str:
    """Site-wise systematic error methods text."""
    return "Site-wise systematic error was diagnosed by evaluating per-site phylogenetic signal contributions."


def generate_methods_posttree_simulate(
    params: dict[str, Any],
    key_results: dict[str, Any],
    tool_versions: dict[str, Any],
) -> str:
    """Simulation methods text."""
    tool = params.get("tool", "alisim")
    n_reps = key_results.get("n_replicates", params.get("replicates", 100))
    return (
        f"Sequence data were simulated using {tool.upper()} with "
        f"{_describe_n(n_reps, 'replicate')}."
    )


# ---------------------------------------------------------------------------
# Methods generator registry
# ---------------------------------------------------------------------------

METHODS_GENERATORS: dict[str, Any] = {
    "pretree.convert": generate_methods_pretree_convert,
    "pretree.stats": generate_methods_pretree_stats,
    "pretree.align": generate_methods_pretree_align,
    "pretree.trim": generate_methods_pretree_trim,
    "pretree.metrics": generate_methods_pretree_metrics,
    "pretree.filter.taper": generate_methods_pretree_filter_taper,
    "pretree.filter.treeshrink": generate_methods_pretree_filter_treeshrink,
    "pretree.filter.symtest": generate_methods_pretree_filter_symtest,
    "pretree.filter.metrics": generate_methods_pretree_filter_metrics,
    "pretree.filter.cluster": generate_methods_pretree_filter_cluster,
    "pretree.concat": generate_methods_pretree_concat,
    "tree.ml.fasttree": generate_methods_tree_ml_fasttree,
    "tree.ml.iqtree": generate_methods_tree_ml_iqtree,
    "tree.msc": generate_methods_tree_msc,
    "tree.bi": generate_methods_tree_bi,
    "tree.cf": generate_methods_tree_cf,
    "posttree.topology": generate_methods_posttree_topology,
    "posttree.dating.hessian": generate_methods_posttree_dating_hessian,
    "posttree.dating.mcmc": generate_methods_posttree_dating_mcmc,
    "posttree.signal": generate_methods_posttree_signal,
    "posttree.syserror.brlen": generate_methods_posttree_syserror_brlen,
    "posttree.syserror.cca": generate_methods_posttree_syserror_cca,
    "posttree.syserror.sites": generate_methods_posttree_syserror_sites,
    "posttree.simulate": generate_methods_posttree_simulate,
}


def generate_all_methods(
    step_id: str,
    params: dict[str, Any],
    key_results: dict[str, Any],
    tool_versions: dict[str, Any],
    status: str = "success",
) -> str:
    """Dispatch to the appropriate methods generator for step_id.

    Returns empty string for unknown step_ids or failed steps.
    """
    if status != "success":
        return ""
    generator = METHODS_GENERATORS.get(step_id)
    if generator is None:
        return ""
    return generator(params, key_results, tool_versions)
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `python -m pytest tests/report/test_templates.py -v`
Expected: all tests PASS

---

## Task 5: Renderer — report.json → report.html via Jinja2

**Files:**
- Create: `phyloai/report/renderer.py`
- Overwrite: `phyloai/report/html/report.html.j2` (full template)
- Create: `tests/report/test_renderer.py`

**Interfaces:**
- Produces: `render_html(report: dict, output_dir: Path) -> Path` — renders `report.html`
- Produces: `_relative_path(abs_path: str, base_dir: Path) -> str` — compute relative path for HTML links

### HTML template requirements (Spec Section 10):
- Self-contained single file (no CDN)
- 5 panels (A-E): Run Summary, Methods, Steps Detail, Figures, Output Files Index
- Collapsible `<details>` cards per step (success=collapsed, failed=expanded)
- PDF figures embedded via `<object>` tags
- Sortable tables via inline JS
- Copy-to-clipboard for methods paragraph
- Monospace styled sections
- Paths relative to report.html output directory

- [ ] **Step 1: Write smoke test for renderer**

In `tests/report/test_renderer.py`:

```python
"""Smoke tests for phyloai.report.renderer."""
from __future__ import annotations

import json
from pathlib import Path

from phyloai.report.renderer import render_html


def test_render_html_produces_valid_html(tmp_path):
    """Smoke test: valid report.json produces non-empty HTML file."""
    report = {
        "phyloai_version": "0.1.0",
        "generated_at": "2026-06-27T00:00:00Z",
        "run_dir": "/tmp/runs/pretree",
        "run_mode": "module",
        "status": "complete",
        "pipeline_summary": {
            "n_steps_total": 1,
            "n_steps_success": 1,
            "n_steps_failed": 0,
            "n_steps_skipped": 0,
            "total_wall_time": 100.0,
        },
        "steps": [
            {
                "step_id": "pretree.align",
                "command": "phyloai pretree align --seq-dir ./raw --method linsi",
                "status": "success",
                "wall_time": 31.4,
                "tool_versions": {"mafft": "7.526"},
                "params": {"method": "linsi", "seq_dir": "./raw", "threads": 8, "seq_type": "AA", "backtrans": False},
                "key_results": {"n_aligned": 100, "n_skipped": 0, "mean_alignment_length": 500.0},
                "methods_text": "Multiple sequence alignments were performed using MAFFT v7.526...",
                "output_files": {},
                "warnings": [],
                "error": None,
            },
        ],
        "methods_paragraph": "Multiple sequence alignments were performed using MAFFT v7.526...",
        "figures_index": [],
        "tables_index": [],
    }
    output_dir = tmp_path / "report"
    output_dir.mkdir()

    # Write report.json first
    (output_dir / "report.json").write_text(json.dumps(report))

    result_path = render_html(report, output_dir)
    assert result_path.exists()
    content = result_path.read_text()
    assert "<!DOCTYPE html>" in content
    assert "PhyloAI Report" in content
    assert "pretree.align" in content
    assert "MAFFT" in content
    # Check all panels present
    assert "Run Summary" in content or "summary" in content.lower()
    assert "Methods" in content
    assert "Steps Detail" in content or "Steps" in content


def test_render_html_with_figures(tmp_path):
    """Report with figures renders figure section."""
    report = {
        "phyloai_version": "0.1.0",
        "generated_at": "2026-06-27T00:00:00Z",
        "run_dir": str(tmp_path),
        "run_mode": "module",
        "status": "complete",
        "pipeline_summary": {"n_steps_total": 1, "n_steps_success": 1, "n_steps_failed": 0, "n_steps_skipped": 0, "total_wall_time": 10.0},
        "steps": [
            {
                "step_id": "pretree.metrics",
                "command": "phyloai pretree metrics ...",
                "status": "success",
                "wall_time": 10.0,
                "tool_versions": {},
                "params": {},
                "key_results": {},
                "methods_text": "Metrics were computed.",
                "output_files": {"correlation_heatmap": {"path": str(tmp_path / "corr.pdf")}},
                "warnings": [],
                "error": None,
            },
        ],
        "methods_paragraph": "Metrics were computed.",
        "figures_index": [
            {
                "figure_id": "Fig-3.1",
                "step_id": "pretree.metrics",
                "caption": "Correlation heatmap",
                "path": str(tmp_path / "corr.pdf"),
                "type": "pdf",
            },
        ],
        "tables_index": [],
    }
    output_dir = tmp_path / "report"
    output_dir.mkdir()
    (output_dir / "report.json").write_text(json.dumps(report))

    result_path = render_html(report, output_dir)
    content = result_path.read_text()
    assert "Fig-3.1" in content
    assert "corr.pdf" in content
    assert '<object' in content or '<figure' in content.lower()


def test_render_html_with_failed_step(tmp_path):
    """Failed step is rendered with [FAILED] indicator."""
    report = {
        "phyloai_version": "0.1.0",
        "generated_at": "2026-06-27T00:00:00Z",
        "run_dir": str(tmp_path),
        "run_mode": "module",
        "status": "partial",
        "pipeline_summary": {"n_steps_total": 2, "n_steps_success": 1, "n_steps_failed": 1, "n_steps_skipped": 0, "total_wall_time": 10.0},
        "steps": [
            {
                "step_id": "pretree.align",
                "command": "phyloai pretree align ...",
                "status": "success",
                "wall_time": 3.0,
                "tool_versions": {},
                "params": {},
                "key_results": {},
                "methods_text": "Align methods.",
                "output_files": {},
                "warnings": [],
                "error": None,
            },
            {
                "step_id": "pretree.trim",
                "command": "phyloai pretree trim ...",
                "status": "error",
                "wall_time": 0.1,
                "tool_versions": {},
                "params": {},
                "key_results": {},
                "methods_text": "",
                "output_files": {},
                "warnings": [],
                "error": "trimAl failed",
            },
        ],
        "methods_paragraph": "Align methods.",
        "figures_index": [],
        "tables_index": [],
    }
    output_dir = tmp_path / "report"
    output_dir.mkdir()
    (output_dir / "report.json").write_text(json.dumps(report))

    result_path = render_html(report, output_dir)
    content = result_path.read_text()
    assert "FAILED" in content or "failed" in content.lower() or "error" in content.lower()
    # Failed step card should be expanded
    assert '<details open' in content.lower()


def test_output_files_index_includes_non_figure_files(tmp_path):
    """Panel E includes FASTA, Newick, and other non-figure/non-table files."""
    fasta_path = tmp_path / "matrix.fa"
    fasta_path.write_text(">seq\nATCG")
    nwk_path = tmp_path / "tree.nwk"
    nwk_path.write_text("(A,B);")
    report = {
        "phyloai_version": "0.1.0",
        "generated_at": "2026-06-27T00:00:00Z",
        "run_dir": str(tmp_path),
        "run_mode": "module",
        "status": "complete",
        "pipeline_summary": {"n_steps_total": 1, "n_steps_success": 1, "n_steps_failed": 0, "n_steps_skipped": 0, "total_wall_time": 10.0},
        "steps": [{
            "step_id": "pretree.concat",
            "command": "phyloai pretree concat ...",
            "status": "success", "wall_time": 1.0,
            "tool_versions": {}, "params": {}, "key_results": {},
            "methods_text": "Concat methods.",
            "output_files": {
                "matrix": {"path": str(fasta_path), "description": "Supermatrix"},
                "tree": {"path": str(nwk_path), "description": "Species tree"},
            },
            "warnings": [], "error": None,
        }],
        "methods_paragraph": "Concat methods.",
        "figures_index": [],
        "tables_index": [],
    }
    output_dir = tmp_path / "report"
    output_dir.mkdir()
    (output_dir / "report.json").write_text(json.dumps(report))

    result_path = render_html(report, output_dir)
    content = result_path.read_text()
    assert "matrix.fa" in content
    assert "tree.nwk" in content
    # Both should appear in the Output Files Index
    assert content.count("Supermatrix") >= 1
    assert content.count("Species tree") >= 1
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `python -m pytest tests/report/test_renderer.py -v`
Expected: FAIL (renderer.py not yet implemented)

- [ ] **Step 3: Implement `phyloai/report/renderer.py`**

```python
"""HTML report renderer using Jinja2 templates."""

from __future__ import annotations

import json
import os
from pathlib import Path
from typing import Any

from jinja2 import Environment, FileSystemLoader, select_autoescape


_TEMPLATE_DIR = Path(__file__).parent / "html"


def _relative_path(abs_path: str, base_dir: Path) -> str:
    """Compute a relative path from base_dir to abs_path for HTML links.

    Uses os.path.relpath which handles paths outside base_dir by inserting
    ``..`` components.  Falls back to abs_path only on different drives
    (Windows edge case).
    """
    try:
        return os.path.relpath(abs_path, str(base_dir))
    except ValueError:
        return abs_path


def _scientific_params(params: dict[str, Any]) -> dict[str, Any]:
    """Filter out technical parameters, keep only scientifically meaningful ones."""
    technical = {
        "threads", "output_dir", "run_dir", "overwrite", "resume",
        "dry_run", "quiet", "seq_dir", "msa_dir", "tree_dir",
        "mafft_path", "magus_path", "trimal_path", "iqtree_path",
        "fasttree_path", "wastral_path", "tool_args",
        "input", "output", "input_format", "table_format",
    }
    return {k: v for k, v in params.items() if k not in technical}


def render_html(report: dict[str, Any], output_dir: Path) -> Path:
    """Render report.json to report.html using Jinja2.

    All file paths in the output are pre-computed as relative to output_dir
    so the template never needs to call path logic.
    """
    env = Environment(
        loader=FileSystemLoader(str(_TEMPLATE_DIR)),
        autoescape=select_autoescape(["html"]),
    )
    template = env.get_template("report.html.j2")

    # Compute relative paths for all figures and tables
    figures = []
    for fig in report.get("figures_index", []):
        fig_copy = dict(fig)
        fig_copy["rel_path"] = _relative_path(fig["path"], output_dir)
        figures.append(fig_copy)

    tables = []
    for tbl in report.get("tables_index", []):
        tbl_copy = dict(tbl)
        tbl_copy["rel_path"] = _relative_path(tbl["path"], output_dir)
        tables.append(tbl_copy)

    # Enrich steps with scientific params subset and relative output file paths
    steps = []
    for step in report.get("steps", []):
        step_copy = dict(step)
        step_copy["scientific_params"] = _scientific_params(step.get("params", {}))
        step_copy["n_params"] = len(step_copy["scientific_params"])
        of = {}
        for label, file_obj in step.get("output_files", {}).items():
            of[label] = dict(file_obj)
            of[label]["rel_path"] = _relative_path(file_obj["path"], output_dir)
        step_copy["output_files_rel"] = of
        steps.append(step_copy)

    # Build all_files list for Panel E with pre-computed relative paths
    all_files: list[dict[str, Any]] = []
    for s in steps:
        for label, fo in s.get("output_files_rel", {}).items():
            all_files.append({
                "step_id": s["step_id"],
                "label": label,
                "description": fo.get("description", "\u2014"),
                "path": fo["path"],
                "rel_path": fo["rel_path"],
                "type": Path(fo["path"]).suffix.lstrip(".") if fo.get("path") else "?",
            })

    html = template.render(
        report=report,
        figures=figures,
        tables=tables,
        steps=steps,
        all_files=all_files,
    )

    out_path = output_dir / "report.html"
    out_path.write_text(html, encoding="utf-8")
    return out_path
```

- [ ] **Step 4: Write the full `report.html.j2` Jinja2 template**

Overwrite `phyloai/report/html/report.html.j2`:

```html
<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>PhyloAI Report — {{ report.run_dir }}</title>
<style>
  :root {
    --fg: #1a1a1a; --bg: #fafafa; --muted: #666;
    --border: #ddd; --accent: #2563eb; --success: #16a34a;
    --fail: #dc2626; --warn: #f59e0b; --surface: #fff;
    --code-bg: #f3f4f6; --radius: 6px;
  }
  * { box-sizing: border-box; margin:0; padding:0; }
  body { font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;
         background: var(--bg); color: var(--fg); line-height: 1.6; max-width: 960px;
         margin: 0 auto; padding: 2rem 1.5rem; }
  h1 { font-size: 1.8rem; margin-bottom: 0.25rem; }
  h2 { font-size: 1.3rem; margin: 2rem 0 0.75rem; padding-bottom: 0.25rem; border-bottom: 2px solid var(--accent); }
  h3 { font-size: 1.05rem; margin: 1rem 0 0.5rem; }

  .header { margin-bottom: 1.5rem; }
  .header .meta { color: var(--muted); font-size: 0.85rem; display: flex; gap: 1.5rem; flex-wrap: wrap; }

  .summary-cards { display: flex; gap: 0.75rem; flex-wrap: wrap; margin-bottom: 1rem; }
  .card { background: var(--surface); border: 1px solid var(--border); border-radius: var(--radius);
          padding: 0.75rem 1rem; min-width: 180px; }
  .card .value { font-size: 1.3rem; font-weight: 600; }
  .card .label { font-size: 0.75rem; color: var(--muted); text-transform: uppercase; letter-spacing: 0.05em; }
  .card.fail { border-color: var(--fail); }

  .failed-list { color: var(--fail); margin: 0.5rem 0; font-size: 0.9rem; }

  .methods-block { background: var(--surface); border: 1px solid var(--border); border-radius: var(--radius);
                   padding: 1.25rem; position: relative; line-height: 1.8; }
  .copy-btn { position: absolute; top: 0.5rem; right: 0.5rem; background: var(--code-bg);
              border: 1px solid var(--border); border-radius: 4px; padding: 0.25rem 0.5rem;
              cursor: pointer; font-size: 0.75rem; }
  .copy-btn:hover { background: var(--border); }

  details { background: var(--surface); border: 1px solid var(--border);
            border-radius: var(--radius); margin-bottom: 0.5rem; }
  details.failed { border-color: var(--fail); }
  details summary { padding: 0.6rem 1rem; cursor: pointer; font-weight: 500; }
  details summary:hover { background: var(--code-bg); }
  details .detail-body { padding: 0 1rem 1rem; }
  details .step-status { font-weight: 600; }
  details .step-status.ok { color: var(--success); }
  details .step-status.fail { color: var(--fail); }

  table { width: 100%; border-collapse: collapse; margin: 0.5rem 0; font-size: 0.85rem; }
  th, td { text-align: left; padding: 0.4rem 0.5rem; border-bottom: 1px solid var(--border); }
  th { background: var(--code-bg); font-weight: 600; cursor: pointer; user-select: none; }
  th:hover { background: var(--border); }
  tr:hover td { background: var(--code-bg); }
  .param-table { max-height: 300px; overflow-y: auto; }

  figure { margin: 1.5rem 0; }
  figure object { display: block; border: 1px solid var(--border); border-radius: var(--radius); }
  figcaption { font-size: 0.8rem; color: var(--muted); margin-top: 0.3rem; }
  figcaption strong { color: var(--fg); }
  figcaption a { color: var(--accent); word-break: break-all; }

  .mono { font-family: "SF Mono", "Fira Code", monospace; font-size: 0.82rem; }
  .cmd-block { background: var(--code-bg); padding: 0.5rem 0.75rem; border-radius: var(--radius);
               overflow-x: auto; white-space: pre-wrap; word-break: break-all; margin: 0.5rem 0; }
  .warnings { background: #fffbeb; border: 1px solid var(--warn); border-radius: var(--radius);
              padding: 0.5rem 0.75rem; margin: 0.5rem 0; font-size: 0.85rem; }
  .error-msg { background: #fef2f2; border: 1px solid var(--fail); border-radius: var(--radius);
               padding: 0.5rem 0.75rem; margin: 0.5rem 0; font-size: 0.85rem; color: var(--fail); }

  .progress-bar { display: flex; gap: 2px; margin: 0.5rem 0; }
  .progress-step { height: 6px; flex: 1; border-radius: 3px; }
  .progress-step.ok { background: var(--success); }
  .progress-step.fail { background: var(--fail); }
  .progress-step.pending { background: var(--border); }

  .back-link { font-size: 0.75rem; color: var(--muted); margin-left: 0.5rem; }
  .back-link::before { content: "↑ "; }

  @media print {
    body { max-width: none; padding: 0; }
    details { border: none; break-inside: avoid; }
    details[open] .detail-body { display: block; }
    .copy-btn { display: none; }
  }
</style>
</head>
<body>

<!-- Header -->
<div class="header">
  <h1>PhyloAI Report</h1>
  <div class="meta">
    <span>Run: <code>{{ report.run_dir.split('/')[-1] or report.run_dir }}</code></span>
    <span>Generated: {{ report.generated_at[:19] }}</span>
    <span>phyloai {{ report.phyloai_version }}</span>
  </div>
</div>

<!-- Panel A: Run Summary -->
<h2>Run Summary</h2>
{% set ps = report.pipeline_summary %}
<div class="summary-cards">
  <div class="card">
    <div class="value">{{ ps.n_steps_total }}</div>
    <div class="label">Total Steps</div>
  </div>
  <div class="card">
    <div class="value" style="color: var(--success)">{{ ps.n_steps_success }}</div>
    <div class="label">Succeeded</div>
  </div>
  {% if ps.n_steps_failed %}
  <div class="card fail">
    <div class="value" style="color: var(--fail)">{{ ps.n_steps_failed }}</div>
    <div class="label">Failed</div>
  </div>
  {% endif %}
  <div class="card">
    <div class="value">
      {% set h = (ps.total_wall_time / 3600)|int %}
      {% set m = ((ps.total_wall_time % 3600) / 60)|int %}
      {% set s = (ps.total_wall_time % 60)|int %}
      {% if h %}{{ h }}h {% endif %}{{ m }}m {{ s }}s
    </div>
    <div class="label">Wall Time</div>
  </div>
</div>

{% set failed_steps = steps | selectattr("status", "equalto", "error") | list %}
{% if failed_steps %}
<div class="failed-list">
  Failed: {{ failed_steps | map(attribute="step_id") | join(", ") }}
</div>
{% endif %}

<!-- Progress bar for pipeline mode -->
{% if report.run_mode == "pipeline" %}
<div class="progress-bar">
{% for s in steps %}
  <div class="progress-step {% if s.status == 'success' %}ok{% elif s.status == 'error' %}fail{% else %}pending{% endif %}"
       title="{{ s.step_id }}: {{ s.status }}"></div>
{% endfor %}
</div>
{% endif %}

<!-- Panel B: Methods -->
<h2>Methods</h2>
<div class="methods-block" id="methods-text">
  <button class="copy-btn" onclick="navigator.clipboard.writeText(document.getElementById('methods-text').innerText.replace('Copy','').trim())">Copy</button>
  {{ report.methods_paragraph }}
</div>

<!-- Panel C: Steps Detail -->
<h2>Steps Detail</h2>
{% for s in steps %}
<details class="{% if s.status == 'error' %}failed{% endif %}" {% if s.status == 'error' %}open{% endif %}>
  <summary>
    <span class="step-status {% if s.status == 'success' %}ok{% else %}fail{% endif %}">
      {{ "✓" if s.status == "success" else "✗" }}
    </span>
    {{ s.step_id }}
    {% if s.tool_versions %}
      {% set primary_tool = s.tool_versions.keys() | first %}
      <span style="color: var(--muted)">· {{ primary_tool }} v{{ s.tool_versions[primary_tool] }}</span>
    {% endif %}
    <span style="color: var(--muted); font-weight: 400;">· {{ "%.1f"|format(s.wall_time) }}s</span>
  </summary>
  <div class="detail-body">
    {% if s.error %}
    <div class="error-msg"><strong>Error:</strong> {{ s.error }}</div>
    {% endif %}

    {% if s.methods_text %}
    <p>{{ s.methods_text }}</p>
    {% endif %}

    {% if s.scientific_params %}
    {% set n = s.scientific_params | length %}
    <details {% if n > 10 %}open{% endif %}>
      <summary>Parameters ({{ n }})</summary>
      <div class="param-table">
        <table>
          <tr><th>Parameter</th><th>Value</th></tr>
          {% for k, v in s.scientific_params.items() %}
          <tr><td><code>{{ k }}</code></td><td><code>{{ v }}</code></td></tr>
          {% endfor %}
        </table>
      </div>
    </details>
    {% endif %}

    {% if s.key_results and s.key_results.keys() | list | length > 0 %}
    <details open>
      <summary>Key Results</summary>
      <table>
        <tr><th>Metric</th><th>Value</th></tr>
        {% for k, v in s.key_results.items() %}
        <tr><td>{{ k }}</td><td>{{ v }}</td></tr>
        {% endfor %}
      </table>
    </details>
    {% endif %}

    {% if s.warnings %}
    <div class="warnings">
      <strong>Warnings:</strong>
      {% for w in s.warnings %}<div>{{ w }}</div>{% endfor %}
    </div>
    {% endif %}

    <details>
      <summary>Full Command</summary>
      <div class="cmd-block mono">{{ s.command }}</div>
    </details>
  </div>
</details>
{% endfor %}

<!-- Panel D: Figures -->
{% if figures %}
<h2>Figures</h2>
{% for fig in figures %}
<figure id="{{ fig.figure_id }}">
  {% if fig.type == "pdf" %}
  <object data="{{ fig.rel_path }}" type="application/pdf" width="100%" height="600px">
    <p>PDF viewer not available. <a href="{{ fig.rel_path }}">Download {{ fig.figure_id }}</a></p>
  </object>
  {% else %}
  <img src="{{ fig.rel_path }}" alt="{{ fig.caption }}" style="max-width:100%; border:1px solid var(--border); border-radius: var(--radius);">
  {% endif %}
  <figcaption>
    <strong>{{ fig.figure_id }}</strong> {{ fig.caption }}
    <span class="back-link">{{ fig.step_id }}</span><br>
    <a href="{{ fig.rel_path }}">{{ fig.path }}</a>
  </figcaption>
</figure>
{% endfor %}
{% endif %}

<!-- Panel E: Output Files Index -->
<h2>Output Files Index</h2>
{% set file_count = all_files | length %}
{% if all_files %}
<div style="margin-bottom:0.5rem; color: var(--muted); font-size: 0.85rem;">
  {{ file_count }} output file{{ "s" if file_count != 1 }} across all steps.
  Source: aggregated from result.json:data.output_files across all steps.
</div>

{% if file_count > 20 %}
<details>
  <summary>Show all {{ file_count }} output files</summary>
{% endif %}

<table id="files-table">
  <thead>
    <tr>
      <th onclick="sortTable(0)">#</th>
      <th onclick="sortTable(1)">Step</th>
      <th onclick="sortTable(2)">Label</th>
      <th onclick="sortTable(3)">Description</th>
      <th onclick="sortTable(4)">File</th>
      <th onclick="sortTable(5)">Type</th>
    </tr>
  </thead>
  <tbody>
    {% for f in all_files %}
    <tr>
      <td>{{ loop.index }}</td>
      <td>{{ f.step_id }}</td>
      <td><code>{{ f.label }}</code></td>
      <td>{{ f.description }}</td>
      <td style="max-width:300px; overflow:hidden; text-overflow:ellipsis; white-space:nowrap;">
        <a href="{{ f.rel_path }}" style="font-size:0.8rem;">{{ f.path }}</a>
      </td>
      <td>{{ f.type }}</td>
    </tr>
    {% endfor %}
  </tbody>
</table>

{% if file_count > 20 %}
</details>
{% endif %}
{% else %}
<p style="color: var(--muted)">No output files recorded.</p>
{% endif %}

<script>
function sortTable(col) {
  const table = document.getElementById("files-table");
  const tbody = table.tBodies[0];
  const rows = Array.from(tbody.rows);
  const asc = table.dataset.sortCol == col ? !(table.dataset.sortAsc == "true") : true;
  rows.sort((a, b) => {
    const va = a.cells[col].textContent.trim();
    const vb = b.cells[col].textContent.trim();
    const na = parseFloat(va), nb = parseFloat(vb);
    if (!isNaN(na) && !isNaN(nb)) return asc ? na - nb : nb - na;
    return asc ? va.localeCompare(vb) : vb.localeCompare(va);
  });
  rows.forEach(r => tbody.appendChild(r));
  table.dataset.sortCol = col;
  table.dataset.sortAsc = asc;
}
</script>

</body>
</html>
```

- [ ] **Step 5: Run tests to verify they pass**

Run: `python -m pytest tests/report/test_renderer.py -v`
Expected: all tests PASS

---

## Task 6: CLI — `phyloai report` Command

**Files:**
- Create: `phyloai/cli/commands/report.py`
- Modify: `phyloai/cli/main.py`
- Create: `tests/report/test_integration.py`

**Interfaces:**
- Produces: Click command `report` with options `--run-dir`, `-o/--output-dir`, `--overwrite`, `-q/--quiet`

- [ ] **Step 1: Create `phyloai/cli/commands/report.py`**

```python
"""phyloai report — generate reproducible analysis reports."""

from __future__ import annotations

import json
from pathlib import Path

import click
from rich.console import Console

from phyloai.report.collector import discover_steps
from phyloai.report.schema import assemble_report
from phyloai.report.templates import generate_all_methods
from phyloai.report.renderer import render_html

console = Console()


def _fail(message: str, exit_code: int = 1) -> None:
    click.echo(f"Error: {message}", err=True)
    raise click.exceptions.Exit(exit_code)


@click.command()
@click.option(
    "--run-dir",
    required=True,
    type=click.Path(exists=True, file_okay=False, path_type=Path),
    help="Run directory to report on (pipeline or module output).",
)
@click.option(
    "-o", "--output-dir",
    type=click.Path(path_type=Path),
    help="Output directory for report files. Default: <run-dir>/report",
)
@click.option(
    "--overwrite",
    is_flag=True,
    help="Overwrite existing report files.",
)
@click.option(
    "-q", "--quiet",
    is_flag=True,
    help="Suppress terminal output except errors.",
)
def report(
    run_dir: Path,
    output_dir: Path | None,
    overwrite: bool,
    quiet: bool,
) -> None:
    """Generate a reproducible analysis report from a PhyloAI run directory.

    Produces report.json (machine-readable, AI/MCP diagnostic entry point)
    and report.html (human-readable, with embedded figures and methods draft).

    \b
    Examples:
      phyloai report --run-dir ./runs/run/faa
      phyloai report --run-dir ./runs/pretree -o ./my-report
    """
    run_dir = run_dir.resolve()

    # Resolve output directory
    if output_dir is None:
        output_dir = run_dir / "report"
    output_dir = output_dir.resolve()

    # Check for existing report files
    report_json_path = output_dir / "report.json"
    report_html_path = output_dir / "report.html"

    if not overwrite and (report_json_path.exists() or report_html_path.exists()):
        _fail(
            f"Report files already exist in {output_dir}. "
            f"Use --overwrite to replace them."
        )

    output_dir.mkdir(parents=True, exist_ok=True)

    # Step 1: Discover steps
    if not quiet:
        console.print(f"[bold]Scanning[/bold] {run_dir}")
    try:
        discovered = discover_steps(run_dir)
    except ValueError as e:
        _fail(str(e))
    except Exception as e:
        _fail(f"Failed to scan run directory: {e}")

    if not quiet:
        console.print(
            f"  Run mode: [bold]{discovered['run_mode']}[/bold] "
            f"({len(discovered['steps'])} steps found)"
        )

    # Step 2: Generate methods text for each step
    methods_texts: dict[str, str] = {}
    for raw_step in discovered["steps"]:
        step_id = raw_step["step_id"]
        status = raw_step.get("status", "error")
        text = generate_all_methods(
            step_id,
            params=raw_step.get("params", {}),
            key_results=raw_step.get("key_results", {}),
            tool_versions=raw_step.get("tool_versions", {}),
            status=status,
        )
        methods_texts[step_id] = text

    # Step 3: Assemble report.json
    if not quiet:
        console.print("[bold]Assembling[/bold] report.json")

    report_dict = assemble_report(discovered, run_dir, methods_texts)

    with open(report_json_path, "w") as fh:
        json.dump(report_dict, fh, indent=2, ensure_ascii=False)

    if not quiet:
        n_ok = sum(1 for s in report_dict["steps"] if s["status"] == "success")
        n_fail = sum(1 for s in report_dict["steps"] if s["status"] == "error")
        status_color = "green" if n_fail == 0 else "yellow"
        console.print(
            f"  Status: [{status_color}]{report_dict['status']}[/{status_color}] "
            f"({n_ok} success, {n_fail} failed)"
        )

    # Step 4: Render report.html from report.json
    if not quiet:
        console.print("[bold]Rendering[/bold] report.html")

    report_dict["run_dir"] = str(run_dir)
    html_path = render_html(report_dict, output_dir)
    if not quiet:
        console.print(f"  [green]report.html[/green] → {html_path}")

    # Final summary
    if not quiet:
        console.print(f"\n[bold green]Report generated:[/bold green]")
        console.print(f"  {report_json_path}")
        console.print(f"  {report_html_path}")
```

- [ ] **Step 2: Register `report` command in `phyloai/cli/main.py`**

Add after the existing imports:

```python
from phyloai.cli.commands.report import report
```

Update `list_commands`:

```python
class _RootGroup(click.Group):
    def list_commands(self, ctx: click.Context) -> list[str]:
        return ["completion", "doctor", "run", "pretree", "tree", "posttree", "report"]
```

Add command:

```python
cli.add_command(report)
```

- [ ] **Step 3: Verify CLI registration**

Run: `python -m phyloai.cli.main --help`
Expected: output includes `report` in the command list.

- [ ] **Step 4: Write integration test**

In `tests/report/test_integration.py`:

```python
"""End-to-end integration tests for phyloai report."""
from __future__ import annotations

import json
from pathlib import Path

from click.testing import CliRunner

from phyloai.cli.commands.report import report


def _make_result_json(
    path: Path,
    step_command: str,
    extra: dict | None = None,
) -> None:
    """Write a minimal valid result.json."""
    data = {
        "status": "success",
        "command": step_command,
        "wall_time": 10.0,
        "tool_versions": {"mafft": "7.526"},
        "params": {"method": "linsi", "seq_type": "AA", "threads": 8, "seq_dir": "./raw"},
        "key_results": {"n_aligned": 100, "n_skipped": 0},
        "error": None,
        "data": {"output_files": {}},
    }
    if extra:
        data.update(extra)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data))


class TestEndToEnd:
    def test_module_run_single_step(self, tmp_path):
        """Generate report from a single-step module run."""
        run_dir = tmp_path / "runs" / "pretree"
        step_dir = run_dir / "2-align"
        _make_result_json(
            step_dir / "result.json",
            "phyloai pretree align --seq-dir ./raw --method linsi --threads 4",
        )

        runner = CliRunner()
        result = runner.invoke(report, ["--run-dir", str(run_dir)])
        assert result.exit_code == 0, f"CLI failed: {result.output}"

        report_dir = run_dir / "report"
        assert (report_dir / "report.json").exists()
        assert (report_dir / "report.html").exists()

        # Validate report.json structure
        rj = json.loads((report_dir / "report.json").read_text())
        assert rj["status"] == "complete"
        assert rj["run_mode"] == "module"
        assert len(rj["steps"]) == 1
        assert rj["steps"][0]["step_id"] == "pretree.align"
        assert "methods_paragraph" in rj
        assert len(rj["methods_paragraph"]) > 0

    def test_module_run_multi_step(self, tmp_path):
        """Generate report from a multi-step module run."""
        run_dir = tmp_path / "runs" / "pretree"
        _make_result_json(
            run_dir / "2-align" / "result.json",
            "phyloai pretree align --seq-dir ./raw --method linsi",
        )
        _make_result_json(
            run_dir / "4-trim" / "result.json",
            "phyloai pretree trim --msa-dir ./aligned --tool bmge",
            extra={"tool_versions": {"bmge": "1.12"}},
        )

        runner = CliRunner()
        result = runner.invoke(report, ["--run-dir", str(run_dir)])
        assert result.exit_code == 0

        report_dir = run_dir / "report"
        rj = json.loads((report_dir / "report.json").read_text())
        assert len(rj["steps"]) == 2
        step_ids = {s["step_id"] for s in rj["steps"]}
        assert step_ids == {"pretree.align", "pretree.trim"}

    def test_failed_step_included(self, tmp_path):
        """Report includes failed steps with error details."""
        run_dir = tmp_path / "runs" / "pretree"
        align_dir = run_dir / "2-align"
        align_dir.mkdir(parents=True)
        (align_dir / "result.json").write_text(json.dumps({
            "status": "error",
            "command": "phyloai pretree align --seq-dir ./raw",
            "wall_time": 0.5,
            "tool_versions": {},
            "params": {},
            "key_results": {},
            "error": "MAFFT returned exit code 1",
            "data": {"output_files": {}},
        }))

        runner = CliRunner()
        result = runner.invoke(report, ["--run-dir", str(run_dir)])
        assert result.exit_code == 0

        rj = json.loads((run_dir / "report" / "report.json").read_text())
        assert rj["status"] == "failed"
        assert rj["steps"][0]["error"] == "MAFFT returned exit code 1"
        assert rj["steps"][0]["methods_text"] == ""
        assert rj["methods_paragraph"] == ""

    def test_overwrite_protection(self, tmp_path):
        """--overwrite required when report files exist."""
        run_dir = tmp_path / "runs" / "pretree"
        _make_result_json(
            run_dir / "2-align" / "result.json",
            "phyloai pretree align --seq-dir ./raw",
        )

        runner = CliRunner()
        # First run succeeds
        r1 = runner.invoke(report, ["--run-dir", str(run_dir)])
        assert r1.exit_code == 0

        # Second run without --overwrite fails
        r2 = runner.invoke(report, ["--run-dir", str(run_dir)])
        assert r2.exit_code != 0
        assert "overwrite" in r2.output.lower()

        # Third run with --overwrite succeeds
        r3 = runner.invoke(report, ["--run-dir", str(run_dir), "--overwrite"])
        assert r3.exit_code == 0

    def test_custom_output_dir(self, tmp_path):
        """-o option changes output location."""
        run_dir = tmp_path / "runs" / "pretree"
        _make_result_json(
            run_dir / "2-align" / "result.json",
            "phyloai pretree align --seq-dir ./raw",
        )
        out_dir = tmp_path / "my-reports" / "run1"

        runner = CliRunner()
        result = runner.invoke(report, [
            "--run-dir", str(run_dir),
            "-o", str(out_dir),
        ])
        assert result.exit_code == 0
        assert (out_dir / "report.json").exists()
        assert (out_dir / "report.html").exists()
```

- [ ] **Step 5: Run integration tests**

Run: `python -m pytest tests/report/test_integration.py -v`
Expected: all tests PASS

---

## Task 7: Full Test Suite and Cleanup

- [ ] **Step 1: Run all report tests**

Run: `python -m pytest tests/report/ -v`
Expected: all tests PASS

- [ ] **Step 2: Run full test suite to check no regressions**

Run: `python -m pytest tests/ -v`
Expected: all existing tests still PASS, new report tests PASS

- [ ] **Step 3: Lint check**

Run: `python -m ruff check phyloai/report/ phyloai/cli/commands/report.py` (if ruff configured)
Or verify code style manually against existing patterns.

---

## Self-Review Checklist

- [x] Spec Section 4 (Data Flow): collector → templates → schema → renderer pipeline implemented
- [x] Spec Section 5 (Directory Detection): Priority-based run_mode detection (4 priorities)
- [x] Spec Section 6 (Step Ordering): STEP_ORDER list, parse_step_id, unknowns appended
- [x] Spec Section 7 (Incomplete Runs): Failed steps included, methods_text="", excluded from paragraph
- [x] Spec Section 8 (report.json): Full schema with all fields (phyloai_version, generated_at, run_dir, run_mode, status, pipeline_summary, steps, methods_paragraph, figures_index, tables_index)
- [x] Spec Section 9 (Methods Templates): All 23 template functions, conditional branches, scientific params only
- [x] Spec Section 10 (HTML): 5 panels (A-E), collapsible cards, PDF embedding, sortable tables, copy button
- [x] Spec Section 11 (Numbering): Phase-based Figure/Table IDs (3=pretree, 4=tree, 5=posttree)
- [x] Spec Section 13 (New Commands): Extension mechanism via METHODS_GENERATORS dict + STEP_ORDER
- [x] Spec Section 14 (CLI): --run-dir, -o, --overwrite, -q flags
- [x] Spec Section 15 (Out of Scope): No tree viz, no cross-run comparison, no LLM, no figure generation, no sub-commands

---

## Task Dependency Graph

```
Task 1 (scaffold)
    ├─→ Task 2 (collector)
    │       └─→ Task 3 (schema)
    │               └─→ Task 6 (CLI + integration)
    ├─→ Task 4 (templates)
    │       └─→ Task 6 (CLI + integration)
    └─→ Task 5 (renderer)
            └─→ Task 6 (CLI + integration)

Task 6 (CLI) depends on Tasks 2, 3, 4, 5
Task 7 (cleanup) runs after all others

Tasks 2, 4, 5 can be developed in parallel after Task 1.
Task 3 depends on Task 2.
```

---

## Implementation Notes (Design Divergences)

The following changes were made during implementation that differ from the original plan:

### `assemble_report` signature change
Removed `methods_texts: dict[str, str]` parameter. Each `raw_step` now carries its own `methods_text` field (populated in `report.py` before calling `assemble_report`). This allows per-result.json granularity — two result.json files for the same `step_id` (e.g. fna + faa) each get their own methods text.

### `methods_blocks` added to `report.json`
New top-level field: `"methods_blocks": [{"step_id": str, "text": str, "step_index": int}, ...]`. Each entry references the step's position in the `steps` array for anchor linking. HTML Panel B renders each block as a separate paragraph with a clickable `[step_id]` badge. The original `methods_paragraph` is retained as plain text for copy-to-clipboard.

### `key_results` enrichment
Some commands (convert, stats) put countable results in `data.summary` or `data.*` rather than `key_results`. Both `report.py` and `schema.py` merge these values into `key_results` before template generation:
- Scalar values (`int`, `float`, `str`) from `data.summary` and `data.*` (excluding structural keys)
- Nested numeric dicts are flattened: `length_before: {mean: 10, max: 20}` → `length_before_mean`, `length_before_max`

### Output files purging
Non-dict entries (legacy bare ints like `n_plots: 56`) are filtered from `output_files` at assembly time. A key-level blacklist (`_SKIP_OF_KEYS = {"n_plots"}`) catches redundant dict-form entries from old result.json files.

### `parse_step_id` algorithm
Uses a known-root lookup table (`pretree`/`tree`/`posttree`/`run`/`doctor`) with registered subcommands. Flag tokens (`--*`) are dropped; the first token matching a known root determines the command path. Boolean flags before the root (`--quiet pretree align`) are handled correctly — only the flag token is dropped, not the subsequent root.

### Pipeline scan depth
Pipeline step discovery is purely filesystem-based (BFS walk, `max_depth=2`). The top-level `result.json:data.steps[]` is only read for optional metadata enrichment (mode, speed), NOT for step `result.json` paths. Excludes `report/`, `logs/`, and dot-prefixed directories.

### HTML template fixes
- `\u00b7` (Python escape, not rendered in static files) replaced with actual `·` character
- Step detail cards have `id="step-{index}"` for Panel B anchor linking
- Small CSV files (≤200 rows, ≤500KB) are embedded inline as sortable HTML tables in Step Detail cards

### Templates
- `stats` template reads actual stats keys (`n_genes`, `total_taxa`, `n_errors`) instead of non-existent keys (`n_files`, `n_taxa_total`)
- `trim` template describes backtrans mode, seq_type, and length before/after changes
- `metrics` template removed GC content, rewrote tree-based metrics description
- `convert` template uses `n_skipped` instead of `n_failed`
- Unit awareness: `sites` for AA, `bp` for NT/CODON in align and stats templates
- Non-tool commands (convert, stats, metrics, filter metrics, filter cluster) include `using phyloai pretree <cmd>` in methods text

### Round 2 — Filter, Tree, Posttree template enrichments

- **Filter templates** completely rewritten with per-module detail:
  - `taper`: version from `correction_multi.jl` key; cutoff, masked sites count
  - `treeshrink`: version from `run_treeshrink.py` key; `α` threshold; removed taxa count
  - `symtest`: IQ-TREE tool name; p-value; retained MSA stats
  - `filter metrics`: `n_total` fallback for input count; condition_failure_counts; retained MSA stats
  - `cluster`: `n_loci` fallback; umap-learn package (not scikit-learn); agglomerative hierarchical clustering (not GMM); cluster sizes per group; outlier detection with `outlier_metric` and `max_drop_fraction`; retained MSA stats
- **Tree ml templates** version keys aligned to actual result.json: `FastTree`, `iqtree3`
  - `fasttree`: discrete rate categories (`-cat`), gamma-distributed rates with branch-length rescaling (`-gamma`), SH-like local support with pseudoreplicates (not bootstrap)
  - `iqtree`: `partitions` presence + `rclusterf` value → partition/merge detection; `none` modelfinder → direct model string (`LG+F+R4`); both `-B` and `--alrt` support
  - `bi`: `pb_mpi` version; CAT-GTR model (mixture=auto); chain lengths; per-chain and pairwise convergence (bpcomp maxdiff, min ESS, max rel_diff, status)
  - `msc`: human-readable mode/boot descriptions (traditional unweighted Astral, quartet+local PP); tree_boot_type description; exhaustive search flag
  - `cf`: full names for gCF/sCF/sCFl/qCF; qCF uses wASTRAL (not IQ-TREE)
- **Posttree templates**:
  - `topology`: model expression (LG+F+R4), RELL replicates, best tree ID, rejected count (AU, p<0.05)
  - `dating.hessian`: partitioned (n_partitions) with auto-selected model (AA→LG+F+G4), no-merge note
  - `dating.mcmc`: clock names (strict/independent-rates/autocorrelated-rates) not numbers; total generations (burnin + nsamples×freq); convergence ρ; diagnostics described: trace, convergence scatter, infinite-sites (posterior=molecular sufficiency, prior=fossil constraint sufficiency), posterior-vs-prior (fossil calibration error detection)

### Round 3 — Module fixes for report data quality

- **Tree ml modules** (`ml.py`, `ml_iqtree.py`): `output_dir.resolve()` added at function entry → all output paths are absolute
- **MSC module** (`msc.py`): `tree_boot_type` resolved from wastral stderr (`bootstrap-like` → `bootstrap`, etc.); original user value preserved in `params`, resolved value in `key_results`
- **BI module** (`tree.py` CLI): `convergence_render.txt` added to `output_files` when present; report module (`schema.py`) has a `_WELL_KNOWN` fallback to detect it for legacy result.json files
- **Metrics CLI** (`pretree.py`): `output_dir.resolve()` added; `n_plots` removed from `output_files` (redundant with `plots_dir`)

### Round 4 — Key_results enrichment extensions

- `data.summary` values of type `str` now merged (was only `int/float`) — fixes seq_type resolution for stats/convert
- `data.summary` values of type `list` now merged — fixes `drop_clusters`, `features` for cluster filter
- `data.summary` nested scalar dicts flattened: `{gap_ratio: {mean: 0.1, median: 0.08}}` → `gap_ratio_mean`, `gap_ratio_median`
- `key_results` itself also flattened for nested scalar dicts (trim's `length_before: {mean, min, max}`)
- Concat-specific metrics extracted from `data.variant_stats[0].character_summary` and `site_patterns` → `gap_ratio`, `pi_ratio`
- Both `report.py` and `schema.py` apply identical enrichment

### Round 5 — HTML and UX polish

- Panel C no longer duplicates methods_text from Panel B; shows `↑ Methods` back-link instead
- Small CSV files (≤200 rows, ≤500KB) embedded inline as sortable HTML tables in Step Detail cards
- `\u00b7` → actual `·` character (Python escape doesn't work in static template files)
- `output_files` purging: non-dict entries filtered; key-level blacklist `_SKIP_OF_KEYS = {"n_plots"}`
- `_relative_path` uses `os.path.relpath` (handles `..` for paths outside report dir; was `Path.relative_to`)
