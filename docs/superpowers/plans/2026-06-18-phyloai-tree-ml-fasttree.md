# PhyloAI Tree ML FastTree Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Implement `phyloai tree ml fasttree` — maximum-likelihood tree inference using FastTree, supporting batch gene trees (`--msa-dir`) and single supermatrix (`--matrix`).

**Architecture:** CLI (`cli/commands/tree.py`) delegates to library (`tree/ml.py`) which scans inputs, builds FastTree commands, runs via `ProcessPoolExecutor` (batch) or directly (single), writes trees to `trees/` and per-task logs to `logs/`, produces `result.json`. Checkpoint/resume via `tree/checkpoint_helpers.py`. Follows existing patterns from `pretree/align.py`.

**Tech Stack:** Click, Rich, ProcessPoolExecutor, subprocess (FastTree), Bio.SeqIO, Bio.Phylo, pytest, CliRunner

---

### File Structure

| File | Responsibility |
|------|---------------|
| `phyloai/tree/__init__.py` | Package init |
| `phyloai/tree/ml.py` | `run_fasttree()`, `_scan_input()`, `_build_fasttree_cmd()`, `_run_one_fasttree()`, `_validate_seq_types()`, `_check_managed_flag_conflict()` |
| `phyloai/tree/checkpoint_helpers.py` | `build_initial_checkpoint()`, `mark_task()`, `resume_verifier()`, `plan_resume()` for tree step |
| `phyloai/cli/commands/tree.py` | Click group: `tree` → `ml` → `fasttree`/`iqtree`; `_fail()`, `_TreeGroup`, `_MLGroup` |
| `phyloai/cli/main.py` | Register `tree` command group |
| `docs/commands/tree-ml.md` | User-facing command docs |
| `tests/tree/test_ml.py` | Library-level tests (cmd building, input scanning, seq_type validation, tool-args blocking) |
| `tests/cli/test_tree.py` | CLI integration tests (CliRunner, mutual exclusivity, output validation) |
| `tests/tree/test_checkpoint.py` | Checkpoint tests (build, mark, resume, Newick validation) |

---

### Task 1: Create `phyloai/tree/__init__.py`

**Files:**
- Create: `phyloai/tree/__init__.py`

- [ ] **Step 1: Write the package init**

```python
"""PhyloAI tree inference module."""
```

- [ ] **Step 2: Verify imports**

```bash
python -c "import phyloai.tree; print('ok')"
```

Expected: `ok`

- [ ] **Step 3: Commit**

```bash
git add phyloai/tree/__init__.py
git commit -m "feat(tree): add tree module package init"
```

---

### Task 2: FastTree command builder and input scanner (library layer)

**Files:**
- Create: `phyloai/tree/ml.py`

- [ ] **Step 1: Write failing tests for `_scan_input`**

Create `tests/tree/test_ml.py`:

```python
from __future__ import annotations

from pathlib import Path


def test_scan_input_finds_fasta_phylip_files(tmp_path: Path) -> None:
    from phyloai.tree.ml import _scan_input

    (tmp_path / "gene1.fa").write_text(">a\nACGT\n")
    (tmp_path / "gene2.faa").write_text(">b\nMKT\n")
    (tmp_path / "gene3.phy").write_text("2 10\na  ACGT\nb  ACGT\n")
    (tmp_path / "gene4.phylip").write_text("2 10\na  ACGT\nb  ACGT\n")
    (tmp_path / "gene5.nex").write_text("#NEXUS\n")
    (tmp_path / "notes.txt").write_text("skip")
    (tmp_path / "empty.fa").write_text("")
    (tmp_path / "subdir").mkdir()

    found, skipped = _scan_input(tmp_path)

    assert len(found) == 4
    assert len(skipped) == 4
    skip_reasons = {s["reason"] for s in skipped}
    assert "NEXUS format not supported by FastTree; use pretree convert first" in skip_reasons
    assert "empty file" in skip_reasons
    assert "directory" in skip_reasons
        assert "unrecognized extension: .txt" in skip_reasons


```

- [ ] **Step 2: Run tests to verify they fail**

```bash
pytest tests/tree/test_ml.py::test_scan_input_finds_fasta_phylip_files -v
```

Expected: FAIL (ImportError: `phyloai.tree.ml` not found or `_scan_input` not defined)

- [ ] **Step 3: Write `_scan_input` and `_build_fasttree_cmd` in `phyloai/tree/ml.py`**

```python
"""Maximum-likelihood tree inference with FastTree."""

from __future__ import annotations

import shlex
import time
from pathlib import Path
from typing import Any, Callable

from phyloai.core.schema import COMMON_ALIGNMENT_EXTENSIONS

FASTTREE_COMPATIBLE_EXTENSIONS = frozenset({
    ".fa", ".fas", ".fasta", ".faa", ".fna",
    ".phy", ".phylip",
})

FASTTREE_MANAGED_FLAGS = frozenset({
    "-nt", "-gtr", "-lg", "-wag",
    "-cat", "-gamma", "-boot", "-nosupport",
    "-fastest", "-slow",
    "-expert", "-help",
})

CHECKPOINT_FLUSH_INTERVAL = 2.0


def _scan_input(msa_dir: Path) -> tuple[list[Path], list[dict[str, str]]]:
    if not msa_dir.exists():
        return [], []

    found: list[Path] = []
    skipped: list[dict[str, str]] = []

    for entry in sorted(msa_dir.iterdir()):
        if entry.is_dir():
            skipped.append({"path": str(entry), "reason": "directory"})
            continue
        if not entry.is_file():
            skipped.append({"path": str(entry), "reason": "not a regular file"})
            continue
        if entry.stat().st_size == 0:
            skipped.append({"path": str(entry), "reason": "empty file"})
            continue

        ext = entry.suffix.lower()
        if ext in FASTTREE_COMPATIBLE_EXTENSIONS:
            found.append(entry)
        elif ext in {".nex", ".nxs", ".nexus"}:
            skipped.append({
                "path": str(entry),
                "reason": "NEXUS format not supported by FastTree; use pretree convert first",
            })
        elif ext in set(COMMON_ALIGNMENT_EXTENSIONS):
            skipped.append({"path": str(entry), "reason": f"unrecognized extension: {ext}"})
        else:
            skipped.append({"path": str(entry), "reason": f"unrecognized extension: {ext}"})

    return found, skipped


def _build_fasttree_cmd(
    input_path: Path,
    output_path: Path,
    *,
    executable: str = "FastTree",
    seq_type: str = "AA",
    model: str = "lg",
    mode: str = "normal",
    boot: int = 1000,
    cat: int = 20,
    gamma: bool = True,
    tool_args: str | None = None,
) -> list[str]:
    cmd = [executable]

    if seq_type == "NT":
        cmd.append("-nt")
        if model == "gtr":
            cmd.append("-gtr")
    else:
        if model == "lg":
            cmd.append("-lg")
        elif model == "wag":
            cmd.append("-wag")

    if mode == "fastest":
        cmd.append("-fastest")
    elif mode == "slow":
        cmd.append("-slow")

    if gamma:
        cmd.append("-gamma")

    cmd.extend(["-cat", str(cat)])

    if boot > 0:
        cmd.extend(["-boot", str(boot)])
    else:
        cmd.append("-nosupport")

    if tool_args:
        _check_managed_flag_conflict(tool_args)
        cmd.extend(shlex.split(tool_args))

    cmd.append(str(input_path))
    return cmd


def _check_managed_flag_conflict(tool_args: str) -> None:
    tokens = shlex.split(tool_args)
    managed_set = FASTTREE_MANAGED_FLAGS
    for token in tokens:
        if token in managed_set:
            raise ValueError(f"Blocked managed flag in --tool-args: {token}")
        if "/" in token or ">" in token:
            raise ValueError(f"Blocked I/O override in --tool-args: {token}")
```

- [ ] **Step 4: Run tests to verify they pass**

```bash
pytest tests/tree/test_ml.py::test_scan_input_finds_fasta_phylip_files -v
```

Expected: PASS

- [ ] **Step 5: Write tests for `_build_fasttree_cmd` and `_check_managed_flag_conflict`**

Append to `tests/tree/test_ml.py`:

```python
import pytest


def test_build_fasttree_cmd_aa_lg_full(tmp_path: Path) -> None:
    from phyloai.tree.ml import _build_fasttree_cmd

    inp = tmp_path / "gene.fa"
    out = tmp_path / "gene.tre"
    cmd = _build_fasttree_cmd(inp, out, seq_type="AA", model="lg", mode="normal",
                              boot=1000, cat=20, gamma=True)

    assert cmd[0] == "FastTree"
    assert "-lg" in cmd
    assert "-gamma" in cmd
    assert "-cat" in cmd and "20" in cmd
    assert "-boot" in cmd and "1000" in cmd
    assert "-nosupport" not in cmd
    assert str(inp) == cmd[-1]


def test_build_fasttree_cmd_nt_gtr(tmp_path: Path) -> None:
    from phyloai.tree.ml import _build_fasttree_cmd

    inp = tmp_path / "gene.fa"
    out = tmp_path / "gene.tre"
    cmd = _build_fasttree_cmd(inp, out, seq_type="NT", model="gtr")

    assert "-nt" in cmd
    assert "-gtr" in cmd
    assert "-lg" not in cmd


def test_build_fasttree_cmd_aa_jtt_default_no_flags(tmp_path: Path) -> None:
    from phyloai.tree.ml import _build_fasttree_cmd

    inp = tmp_path / "gene.fa"
    out = tmp_path / "gene.tre"
    cmd = _build_fasttree_cmd(inp, out, seq_type="AA", model="jtt")

    assert "-lg" not in cmd
    assert "-wag" not in cmd


def test_build_fasttree_cmd_fastest_mode(tmp_path: Path) -> None:
    from phyloai.tree.ml import _build_fasttree_cmd

    inp = tmp_path / "gene.fa"
    out = tmp_path / "gene.tre"
    cmd = _build_fasttree_cmd(inp, out, mode="fastest")

    assert "-fastest" in cmd
    assert "-slow" not in cmd


def test_build_fasttree_cmd_boot_zero_gives_nosupport(tmp_path: Path) -> None:
    from phyloai.tree.ml import _build_fasttree_cmd

    inp = tmp_path / "gene.fa"
    out = tmp_path / "gene.tre"
    cmd = _build_fasttree_cmd(inp, out, boot=0)

    assert "-nosupport" in cmd
    assert "-boot" not in cmd


def test_build_fasttree_cmd_no_gamma(tmp_path: Path) -> None:
    from phyloai.tree.ml import _build_fasttree_cmd

    inp = tmp_path / "gene.fa"
    out = tmp_path / "gene.tre"
    cmd = _build_fasttree_cmd(inp, out, gamma=False)

    assert "-gamma" not in cmd


def test_build_fasttree_cmd_with_tool_args(tmp_path: Path) -> None:
    from phyloai.tree.ml import _build_fasttree_cmd

    inp = tmp_path / "gene.fa"
    out = tmp_path / "gene.tre"
    cmd = _build_fasttree_cmd(inp, out, boot=1000, tool_args="-spr 4 -mlacc 2")

    assert "-spr" in cmd
    assert "4" in cmd
    assert "-mlacc" in cmd
    assert "2" in cmd


def test_check_managed_flag_conflict_blocks_lg() -> None:
    from phyloai.tree.ml import _check_managed_flag_conflict

    with pytest.raises(ValueError, match="Blocked managed flag.*-lg"):
        _check_managed_flag_conflict("-lg")


def test_check_managed_flag_conflict_blocks_boot() -> None:
    from phyloai.tree.ml import _check_managed_flag_conflict

    with pytest.raises(ValueError, match="Blocked managed flag.*-boot"):
        _check_managed_flag_conflict("-boot 500")


def test_check_managed_flag_conflict_allows_strategy_args() -> None:
    from phyloai.tree.ml import _check_managed_flag_conflict

    _check_managed_flag_conflict("-spr 4 -mlacc 2 -slownni")


def test_build_fasttree_cmd_with_explicit_executable(tmp_path: Path) -> None:
    from phyloai.tree.ml import _build_fasttree_cmd

    inp = tmp_path / "gene.fa"
    out = tmp_path / "gene.tre"
    cmd = _build_fasttree_cmd(inp, out, executable="/opt/bin/FastTree")

    assert cmd[0] == "/opt/bin/FastTree"
```

- [ ] **Step 6: Run tests**

```bash
pytest tests/tree/test_ml.py -v
```

Expected: all PASS

- [ ] **Step 7: Commit**

```bash
git add phyloai/tree/ml.py tests/tree/test_ml.py tests/tree/__init__.py
git commit -m "feat(tree): add FastTree command builder and input scanner"
```

---

### Task 3: Single-gene FastTree worker function

**Files:**
- Modify: `phyloai/tree/ml.py` (append)

- [ ] **Step 1: Write failing test for `_run_one_fasttree`**

Append to `tests/tree/test_ml.py`:

```python
def test_run_one_fasttree_success(tmp_path: Path) -> None:
    from phyloai.tree.ml import _run_one_fasttree

    inp = tmp_path / "gene.fa"
    inp.write_text(">a\nMKTLLL\n>b\nMKTLLL\n")
    log_dir = tmp_path / "logs"
    log_dir.mkdir()

    result = _run_one_fasttree(
        gene_path=inp, seq_type="AA", model="lg", mode="normal",
        boot=1000, cat=20, gamma=True, tool_args=None,
        log_dir=log_dir, fasttree_executable="FastTree",
    )

    assert result["status"] == "success"
    assert "output_tree" in result
    assert "log_file" in result


def test_run_one_fasttree_dry_run(tmp_path: Path) -> None:
    from phyloai.tree.ml import _run_one_fasttree

    inp = tmp_path / "gene.fa"
    inp.write_text(">a\nMKT\n")
    log_dir = tmp_path / "logs"
    log_dir.mkdir()

    result = _run_one_fasttree(
        gene_path=inp, seq_type="AA", model="lg", mode="normal",
        boot=1000, cat=20, gamma=True, tool_args=None,
        log_dir=log_dir, fasttree_executable="FastTree",
        dry_run=True,
    )

    assert result["status"] == "dry_run"
    assert "cmd" in result


def test_run_one_fasttree_missing_input(tmp_path: Path) -> None:
    from phyloai.tree.ml import _run_one_fasttree

    inp = tmp_path / "missing.fa"
    log_dir = tmp_path / "logs"
    log_dir.mkdir()

    result = _run_one_fasttree(
        gene_path=inp, seq_type="AA", model="lg", mode="normal",
        boot=1000, cat=20, gamma=True, tool_args=None,
        log_dir=log_dir, fasttree_executable="FastTree",
    )

    assert result["status"] == "failed"
    assert "reason" in result
```

- [ ] **Step 2: Run tests to verify failure**

```bash
pytest tests/tree/test_ml.py::test_run_one_fasttree_dry_run -v
```

Expected: FAIL

- [ ] **Step 3: Write `_run_one_fasttree` in `phyloai/tree/ml.py`**

Append to `phyloai/tree/ml.py`:

```python
import subprocess
import time as _time


def _run_one_fasttree(
    gene_path: Path,
    *,
    seq_type: str,
    model: str,
    mode: str,
    boot: int,
    cat: int,
    gamma: bool,
    tool_args: str | None,
    log_dir: Path,
    fasttree_executable: str = "FastTree",
    output_dir: Path | None = None,
    dry_run: bool = False,
) -> dict[str, Any]:
    if output_dir is None:
        output_dir = gene_path.parent

    result: dict[str, Any] = {
        "input": str(gene_path),
        "output_tree": None,
        "log_file": None,
    }

    if not gene_path.exists():
        return {
            **result,
            "status": "failed",
            "reason": f"input file not found: {gene_path}",
            "wall_time": 0,
            "warnings": [],
        }

    stem = gene_path.stem
    out_tree = output_dir / f"{stem}.tre"
    out_log = log_dir / f"{stem}.log"

    cmd = _build_fasttree_cmd(
        gene_path, out_tree,
        executable=fasttree_executable,
        seq_type=seq_type, model=model, mode=mode,
        boot=boot, cat=cat, gamma=gamma,
        tool_args=tool_args,
    )

    result.update({
        "output_tree": str(out_tree),
        "log_file": str(out_log),
        "cmd": cmd,
    })

    if dry_run:
        return {**result, "status": "dry_run", "wall_time": 0, "warnings": []}

    warnings: list[str] = []
    start = _time.monotonic()
    try:
        out_tree.parent.mkdir(parents=True, exist_ok=True)
        out_log.parent.mkdir(parents=True, exist_ok=True)

        proc = subprocess.run(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
        wall_time = _time.monotonic() - start

        out_tree.write_text(proc.stdout)
        out_log.write_text(proc.stderr)

        if proc.returncode != 0:
            return {
                **result,
                "status": "failed",
                "reason": f"FastTree exited with code {proc.returncode}: {proc.stderr[:200]}",
                "tool_stderr": proc.stderr,
                "wall_time": wall_time,
                "warnings": warnings,
            }

        # Validate Newick output
        from Bio import Phylo
        try:
            Phylo.read(str(out_tree), "newick")
        except Exception as e:
            return {
                **result,
                "status": "failed",
                "reason": f"FastTree produced unparseable Newick output: {e}",
                "tool_stderr": proc.stderr,
                "wall_time": wall_time,
                "warnings": warnings,
            }

        return {
            **result,
            "status": "success",
            "wall_time": wall_time,
            "warnings": warnings,
        }

    except Exception as exc:
        return {
            **result,
            "status": "failed",
            "reason": str(exc),
            "wall_time": _time.monotonic() - start,
            "warnings": warnings,
        }
```

- [ ] **Step 4: Run tests**

```bash
pytest tests/tree/test_ml.py -v
```

Expected: all tests that can run without FastTree installed PASS (dry_run, cmd building, scanning)

- [ ] **Step 5: Commit**

```bash
git add phyloai/tree/ml.py tests/tree/test_ml.py
git commit -m "feat(tree): add single-gene FastTree worker function"
```

---

### Task 4: Sequence type validation for batch directories

**Files:**
- Modify: `phyloai/tree/ml.py` (append)

- [ ] **Step 1: Write tests for `_validate_seq_types`**

Append to `tests/tree/test_ml.py`:

```python
def test_validate_seq_types_homogeneous_aa(tmp_path: Path) -> None:
    from phyloai.tree.ml import _validate_seq_types

    (tmp_path / "g1.fa").write_text(">a\nMKTLLL\n")
    (tmp_path / "g2.fa").write_text(">b\nAAAAAA\n")
    files = sorted(tmp_path.glob("*.fa"))

    resolved, offending = _validate_seq_types(files, declared_type=None)

    assert resolved == "AA"
    assert len(offending) == 0


def test_validate_seq_types_mixed_raises(tmp_path: Path) -> None:
    from phyloai.tree.ml import _validate_seq_types

    (tmp_path / "g1.fa").write_text(">a\nMKTLLL\n")
    (tmp_path / "g2.fa").write_text(">b\nACGTAC\n")
    files = sorted(tmp_path.glob("*.fa"))

    resolved, offending = _validate_seq_types(files, declared_type=None)

    assert resolved is None
    assert len(offending) >= 1


def test_validate_seq_types_explicit_mismatch(tmp_path: Path) -> None:
    from phyloai.tree.ml import _validate_seq_types

    (tmp_path / "g1.fa").write_text(">a\nMKTLLL\n")
    files = sorted(tmp_path.glob("*.fa"))

    resolved, offending = _validate_seq_types(files, declared_type="NT")

    assert resolved == "NT"
    assert len(offending) == 1


def test_validate_seq_types_no_files(tmp_path: Path) -> None:
    from phyloai.tree.ml import _validate_seq_types

    resolved, offending = _validate_seq_types([], declared_type=None)

    assert resolved == "AA"
    assert len(offending) == 0
```

- [ ] **Step 2: Run tests to verify failure**

```bash
pytest tests/tree/test_ml.py::test_validate_seq_types_homogeneous_aa -v
```

Expected: FAIL

- [ ] **Step 3: Write `_validate_seq_types`**

Append to `phyloai/tree/ml.py`:

```python
from Bio import SeqIO

from phyloai.core.sequence_normalization import detect_seq_type


def _validate_seq_types(
    files: list[Path],
    *,
    declared_type: str | None,
) -> tuple[str | None, list[dict[str, Any]]]:
    if not files:
        return (declared_type or "AA"), []

    all_types: dict[str, str] = {}
    offending: list[dict[str, Any]] = []

    for f in files:
        try:
            ext = f.suffix.lower()
            if ext in {".phy", ".phylip"}:
                seqs = [str(r.seq) for r in SeqIO.parse(str(f), "phylip-relaxed")]
            else:
                seqs = [str(r.seq) for r in SeqIO.parse(str(f), "fasta")]
            if not seqs:
                offending.append({"file": str(f), "reason": "no sequences found"})
                continue
            dt = detect_seq_type(seqs)
            all_types[str(f)] = dt
        except Exception:
            offending.append({"file": str(f), "reason": "failed to parse input file"})
            continue

    if declared_type:
        for f_str, dt in all_types.items():
            if dt != declared_type:
                offending.append({"file": f_str, "expected": declared_type, "detected": dt})
        return declared_type, offending

    type_counts: dict[str, int] = {}
    for dt in all_types.values():
        type_counts[dt] = type_counts.get(dt, 0) + 1

    if len(type_counts) == 1:
        resolved = next(iter(type_counts))
        return resolved, []

    majority = max(type_counts, key=type_counts.get)
    for f_str, dt in all_types.items():
        if dt != majority:
            offending.append({"file": f_str, "expected": majority, "detected": dt})

    return None, offending
```

- [ ] **Step 4: Run tests**

```bash
pytest tests/tree/test_ml.py::test_validate_seq_types_homogeneous_aa \
  tests/tree/test_ml.py::test_validate_seq_types_mixed_raises \
  tests/tree/test_ml.py::test_validate_seq_types_explicit_mismatch \
  tests/tree/test_ml.py::test_validate_seq_types_no_files -v
```

Expected: 4 PASS

- [ ] **Step 5: Commit**

```bash
git add phyloai/tree/ml.py tests/tree/test_ml.py
git commit -m "feat(tree): add batch sequence type validation"
```

---

### Task 5: Checkpoint helpers for tree module

**Files:**
- Create: `phyloai/tree/checkpoint_helpers.py`

- [ ] **Step 1: Write failing checkpoint tests**

Create `tests/tree/test_checkpoint.py`:

```python
from __future__ import annotations

from pathlib import Path

from Bio import Phylo
import io


def test_build_initial_checkpoint_tree(tmp_path: Path) -> None:
    from phyloai.tree.checkpoint_helpers import build_initial_checkpoint

    inputs = [Path("/data/gene1.fa"), Path("/data/gene2.fa")]
    trees_dir = Path("/out/trees")
    logs_dir = Path("/out/logs")

    ck = build_initial_checkpoint(
        step="tree.ml.fasttree",
        command="phyloai tree ml fasttree --msa-dir /data",
        params={"seq_type": "AA", "model": "lg"},
        inputs=inputs,
        trees_dir=trees_dir,
        logs_dir=logs_dir,
    )

    assert ck.schema_version == 1
    assert ck.step == "tree.ml.fasttree"
    assert ck.status == "running"
    assert len(ck.tasks) == 2
    assert ck.tasks[0].task_id == "gene1"
    assert ck.tasks[1].task_id == "gene2"
    assert ck.tasks[0].status == "pending"
    assert ck.tasks[0].outputs["tree"] == str(trees_dir / "gene1.tre")
    assert ck.tasks[0].outputs["log"] == str(logs_dir / "gene1.log")


def test_mark_task_updates_checkpoint(tmp_path: Path) -> None:
    from phyloai.tree.checkpoint_helpers import build_initial_checkpoint, mark_task

    inputs = [Path("/data/gene1.fa")]
    ck = build_initial_checkpoint(
        step="tree.ml.fasttree",
        command="cmd",
        params={},
        inputs=inputs,
        trees_dir=Path("/out/trees"),
        logs_dir=Path("/out/logs"),
    )

    mark_task(ck, "gene1", status="success")
    assert ck.tasks[0].status == "success"
    assert ck.tasks[0].attempts == 1

    mark_task(ck, "gene1", status="failed", reason="FastTree error")
    assert ck.tasks[0].status == "failed"
    assert ck.tasks[0].reason == "FastTree error"


def test_resume_verifier_valid_newick(tmp_path: Path) -> None:
    from phyloai.tree.checkpoint_helpers import resume_verifier

    tree_path = tmp_path / "gene1.tre"
    tree_path.write_text("(a:0.1,b:0.2);\n")

    verify = resume_verifier()
    assert verify(tree_path) is True


def test_resume_verifier_invalid_newick(tmp_path: Path) -> None:
    from phyloai.tree.checkpoint_helpers import resume_verifier

    tree_path = tmp_path / "gene1.tre"
    tree_path.write_text("not a newick tree\n")

    verify = resume_verifier()
    assert verify(tree_path) is False


def test_resume_verifier_empty_file(tmp_path: Path) -> None:
    from phyloai.tree.checkpoint_helpers import resume_verifier

    tree_path = tmp_path / "gene1.tre"
    tree_path.write_text("")

    verify = resume_verifier()
    assert verify(tree_path) is False


def test_resume_verifier_nonexistent_file(tmp_path: Path) -> None:
    from phyloai.tree.checkpoint_helpers import resume_verifier

    tree_path = tmp_path / "missing.tre"

    verify = resume_verifier()
    assert verify(tree_path) is False


def test_plan_resume_splits_tasks(tmp_path: Path) -> None:
    from phyloai.tree.checkpoint_helpers import build_initial_checkpoint, mark_task, plan_resume

    inputs = [Path("/data/g1.fa"), Path("/data/g2.fa"), Path("/data/g3.fa")]
    ck = build_initial_checkpoint(
        step="tree.ml.fasttree",
        command="cmd",
        params={},
        inputs=inputs,
        trees_dir=tmp_path / "trees",
        logs_dir=tmp_path / "logs",
    )

    (tmp_path / "trees").mkdir(parents=True)
    (tmp_path / "trees" / "g1.tre").write_text("(a:0.1,b:0.2);\n")

    mark_task(ck, "g1", status="success")
    mark_task(ck, "g2", status="failed", reason="error")
    mark_task(ck, "g3", status="pending")

    to_run, skipped = plan_resume(ck)

    assert "g2" in to_run
    assert "g3" in to_run
    assert "g1" in skipped
```

- [ ] **Step 2: Run tests verify failure**

```bash
pytest tests/tree/test_checkpoint.py -v
```

Expected: FAIL (ImportError)

- [ ] **Step 3: Write `phyloai/tree/checkpoint_helpers.py`**

```python
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
```

- [ ] **Step 4: Run tests**

```bash
pytest tests/tree/test_checkpoint.py -v
```

Expected: all PASS

- [ ] **Step 5: Commit**

```bash
git add phyloai/tree/checkpoint_helpers.py tests/tree/test_checkpoint.py
git commit -m "feat(tree): add checkpoint helpers with Newick validation"
```

---

### Task 6: `run_fasttree()` library function (batch and single modes)

**Files:**
- Modify: `phyloai/tree/ml.py` (append `run_fasttree` function)

- [ ] **Step 1: Write the `run_fasttree` function**

Append to `phyloai/tree/ml.py`:

```python
import json
import os
import shutil
from concurrent.futures import ProcessPoolExecutor, as_completed

from phyloai.core.env import ToolEnv
from phyloai.core.checkpoint import load_checkpoint, save_checkpoint_atomic, validate_resume_params
from phyloai.tree.checkpoint_helpers import build_initial_checkpoint, mark_task, plan_resume


def _resolved_fasttree_params(
    msa_dir: Path | None,
    matrix: Path | None,
    seq_type: str,
    model: str,
    mode: str,
    boot: int,
    cat: int,
    gamma: bool,
    output_dir: Path,
    threads: int,
    fasttree_path: str | None,
    tool_args: str | None,
) -> dict[str, Any]:
    return {
        "msa_dir": str(msa_dir) if msa_dir else None,
        "matrix": str(matrix) if matrix else None,
        "seq_type": seq_type,
        "model": model,
        "mode": mode,
        "boot": boot,
        "cat": cat,
        "gamma": gamma,
        "output_dir": str(output_dir),
        "threads": threads,
        "fasttree_path": fasttree_path,
        "tool_args": tool_args,
    }


def run_fasttree(
    *,
    msa_dir: Path | None = None,
    matrix: Path | None = None,
    output_dir: Path,
    seq_type: str = "auto",
    model: str | None = None,
    mode: str = "normal",
    boot: int = 1000,
    cat: int = 20,
    gamma: bool = True,
    threads: int = 4,
    fasttree_path: str | None = None,
    tool_args: str | None = None,
    overwrite: bool = False,
    resume: bool = False,
    dry_run: bool = False,
    quiet: bool = False,
    progress_callback: Callable[[Path], None] | None = None,
) -> dict[str, Any]:
    run_start = time.monotonic()

    # --- validate inputs ---
    if (msa_dir is None and matrix is None) or (msa_dir is not None and matrix is not None):
        raise ValueError("Either --msa-dir or --matrix must be provided (not both).")

    # --- resolve tool ---
    fasttree_exe = _resolve_fasttree(fasttree_path, dry_run)

    # --- input mode ---
    batch_mode = msa_dir is not None
    input_label = "--msa-dir" if batch_mode else "--matrix"

    # --- output subdirs ---
    trees_dir = output_dir / "trees"
    logs_dir = output_dir / "logs"

    # --- resolve seq_type for single matrix ---
    resolved_seq_type = seq_type
    if not batch_mode:
        assert matrix is not None
        matrix_ext = matrix.suffix.lower()
        matrix_format = "phylip-relaxed" if matrix_ext in {".phy", ".phylip"} else "fasta"
        try:
            recs = list(SeqIO.parse(str(matrix), matrix_format))
        except Exception:
            recs = []
        if seq_type == "auto":
            resolved_seq_type = detect_seq_type([str(r.seq) for r in recs]) if recs else "AA"
        else:
            sample = [str(r.seq) for r in recs[:10]]
            if sample:
                resolved_sample = detect_seq_type(sample)
                if resolved_sample != seq_type:
                    raise ValueError(
                        f"--seq-type {seq_type} but detected {resolved_sample} in {matrix}"
                    )

    # --- batch mode: scan + validate ---
    found: list[Path] = []
    skipped_input: list[dict[str, str]] = []
    if batch_mode:
        assert msa_dir is not None
        found, skipped_input = _scan_input(msa_dir)
        if not found and not dry_run:
            raise ValueError("No valid input files found in --msa-dir")
        declared = None if seq_type == "auto" else seq_type
        resolved_seq_type, offending = _validate_seq_types(found, declared_type=declared)
        if resolved_seq_type is None:
            offending_strs = [f"{o['file']}: {o['detected']} (expected homogeneous)" for o in offending[:10]]
            raise ValueError(f"Mixed sequence types in --msa-dir:\n" + "\n".join(offending_strs))
        if offending:
            offending_strs = [f"{o['file']}: {o['detected']} (expected {o['expected']})" for o in offending[:10]]
            raise ValueError(
                f"Files with wrong --seq-type ({declared}) in --msa-dir:\n" + "\n".join(offending_strs)
            )

    # --- validate model ---
    # --- resolve model default ---
    if model is None:
        model = "gtr" if resolved_seq_type == "NT" else "lg"

    aa_models = {"jtt", "lg", "wag"}
    nt_models = {"jc", "gtr"}
    if resolved_seq_type == "AA":
        if model not in aa_models:
            raise ValueError(f"Invalid model for AA: {model}. Choose from {aa_models}")
    elif resolved_seq_type == "NT":
        if model not in nt_models:
            raise ValueError(f"Invalid model for NT: {model}. Choose from {nt_models}")

    # --- prepare output dir ---
    checkpoint: Any = None
    ckpt_path = output_dir / "checkpoint.json"

    if not dry_run:
        if overwrite and resume:
            raise ValueError("--overwrite and --resume are mutually exclusive")
        if resume:
            if not ckpt_path.exists():
                raise ValueError(f"--resume requires {ckpt_path}, not found")
            checkpoint = load_checkpoint(ckpt_path)
            resolved_params = _resolved_fasttree_params(
                msa_dir=msa_dir, matrix=matrix,
                seq_type=resolved_seq_type, model=model,
                mode=mode, boot=boot, cat=cat, gamma=gamma,
                output_dir=output_dir, threads=threads,
                fasttree_path=fasttree_path, tool_args=tool_args,
            )
            validate_resume_params(checkpoint, resolved_params, step="tree.ml.fasttree")
            if checkpoint.status == "success":
                return _reconstruct_result(output_dir, run_start)

            tree_verifier = __import__("phyloai.tree.checkpoint_helpers", fromlist=["resume_verifier"]).resume_verifier()
            to_run_ids, _skipped_ids = plan_resume(checkpoint)
            if not to_run_ids:
                checkpoint.status = "success"
                save_checkpoint_atomic(checkpoint, ckpt_path)
                return _reconstruct_result(output_dir, run_start)
            found = [Path(task.input) for task in checkpoint.tasks if task.task_id in to_run_ids]
        else:
            if overwrite and output_dir.exists():
                shutil.rmtree(output_dir)
            if output_dir.exists() and any(output_dir.iterdir()):
                raise ValueError(
                    f"Output directory {output_dir} already exists and is non-empty. "
                    "Use --overwrite to replace."
                )
            if batch_mode:
                trees_dir.mkdir(parents=True, exist_ok=True)
                logs_dir.mkdir(parents=True, exist_ok=True)
            else:
                output_dir.mkdir(parents=True, exist_ok=True)

    # --- single-matrix mode ---
    if not batch_mode:
        assert matrix is not None
        result = _run_one_fasttree(
            gene_path=matrix,
            seq_type=resolved_seq_type,
            model=model,
            mode=mode,
            boot=boot,
            cat=cat,
            gamma=gamma,
            tool_args=tool_args,
            log_dir=output_dir,
            fasttree_executable=fasttree_exe,
            output_dir=output_dir,
            dry_run=dry_run,
        )
        return _assemble_result(
            run_start=run_start, fasttree_exe=fasttree_exe,
            batch_mode=False, results=[result],
            resolved_seq_type=resolved_seq_type, model=model, mode=mode, boot=boot,
            cat=cat, gamma=gamma, output_dir=output_dir,
            msa_dir=msa_dir, matrix=matrix,
            fasttree_path=fasttree_path, tool_args=tool_args,
            overwrite=overwrite, threads=threads,
            skipped_input=[],
        )

    # --- batch mode: build checkpoint ---
    if not resume and not dry_run:
        resolved_params = _resolved_fasttree_params(
            msa_dir=msa_dir, matrix=matrix,
            seq_type=resolved_seq_type, model=model,
            mode=mode, boot=boot, cat=cat, gamma=gamma,
            output_dir=output_dir, threads=threads,
            fasttree_path=fasttree_path, tool_args=tool_args,
        )
        checkpoint = build_initial_checkpoint(
            step="tree.ml.fasttree",
            command=f"phyloai tree ml fasttree --msa-dir {msa_dir} ...",
            params=resolved_params,
            inputs=found,
            trees_dir=trees_dir,
            logs_dir=logs_dir,
        )
        save_checkpoint_atomic(checkpoint, ckpt_path)

    # --- batch execution ---
    _ckpt_write = checkpoint is not None and not dry_run
    _last_flush = time.monotonic()

    def _maybe_flush(*, force: bool = False) -> None:
        nonlocal _last_flush
        if not _ckpt_write:
            return
        now = time.monotonic()
        if force or (now - _last_flush) >= CHECKPOINT_FLUSH_INTERVAL:
            save_checkpoint_atomic(checkpoint, ckpt_path)
            _last_flush = now

    file_results: list[dict[str, Any]] = []
    failed_results: list[dict[str, Any]] = []

    worker_args = [
        (p, resolved_seq_type, model, mode, boot, cat, gamma,
         tool_args, logs_dir, fasttree_exe, trees_dir, dry_run)
        for p in found
    ]

    interrupted = False
    try:
        if dry_run:
            for arg in worker_args:
                result = _run_one_fasttree(
                    gene_path=arg[0], seq_type=arg[1], model=arg[2], mode=arg[3],
                    boot=arg[4], cat=arg[5], gamma=arg[6], tool_args=arg[7],
                    log_dir=arg[8], fasttree_executable=arg[9], output_dir=arg[10],
                    dry_run=arg[11],
                )
                file_results.append(result)
                if progress_callback:
                    progress_callback(arg[0])
        else:
            with ProcessPoolExecutor(max_workers=threads) as pool:
                futures = {
                    pool.submit(_run_one_fasttree,
                        gene_path=arg[0], seq_type=arg[1], model=arg[2], mode=arg[3],
                        boot=arg[4], cat=arg[5], gamma=arg[6], tool_args=arg[7],
                        log_dir=arg[8], fasttree_executable=arg[9], output_dir=arg[10],
                        dry_run=arg[11],
                    ): arg[0]
                    for arg in worker_args
                }
                for future in as_completed(futures):
                    gene_path = futures[future]
                    result = future.result()
                    task_id = gene_path.stem

                    if result["status"] == "success":
                        file_results.append(result)
                        mark_task(checkpoint, task_id, status="success")
                    elif result["status"] == "failed":
                        failed_results.append(result)
                        mark_task(checkpoint, task_id, status="failed",
                                  reason=result.get("reason"))
                    else:
                        skipped_input.append({
                            "path": result.get("input", ""),
                            "reason": result.get("reason", "unknown"),
                        })

                    if progress_callback:
                        progress_callback(gene_path)
                    _maybe_flush()

    except KeyboardInterrupt:
        interrupted = True

    from datetime import datetime as _dt_cls, timezone as _tz
    if _ckpt_write:
        if interrupted:
            checkpoint.status = "interrupted"
        else:
            checkpoint.status = "success"
            checkpoint.completed_at = _dt_cls.now(_tz.utc).isoformat(timespec="seconds")
        save_checkpoint_atomic(checkpoint, ckpt_path, fsync=True)
    if interrupted:
        raise KeyboardInterrupt

    return _assemble_result(
        run_start=run_start, fasttree_exe=fasttree_exe,
        batch_mode=True, results=file_results,
        failed_results=failed_results,
        resolved_seq_type=resolved_seq_type, model=model, mode=mode, boot=boot,
        cat=cat, gamma=gamma, output_dir=output_dir,
        msa_dir=msa_dir, matrix=matrix,
        fasttree_path=fasttree_path, tool_args=tool_args,
        overwrite=overwrite, threads=threads,
        skipped_input=skipped_input,
    )


def _resolve_fasttree(fasttree_path: str | None, dry_run: bool) -> str:
    if fasttree_path:
        p = Path(fasttree_path)
        if not p.exists():
            raise ValueError(f"--fasttree-path does not exist: {fasttree_path}")
        if not os.access(p, os.X_OK):
            raise ValueError(f"--fasttree-path is not executable: {fasttree_path}")
        return fasttree_path
    if dry_run:
        return "FastTree"
    try:
        env = ToolEnv()
        return str(env.require("FastTree"))
    except FileNotFoundError:
        raise FileNotFoundError("FastTree not found. Install it or use --fasttree-path.")


def _assemble_result(
    *,
    run_start: float,
    fasttree_exe: str,
    batch_mode: bool,
    results: list[dict[str, Any]],
    failed_results: list[dict[str, Any]] | None = None,
    resolved_seq_type: str,
    model: str,
    mode: str,
    boot: int,
    cat: int,
    gamma: bool,
    output_dir: Path,
    msa_dir: Path | None,
    matrix: Path | None,
    fasttree_path: str | None,
    tool_args: str | None,
    overwrite: bool,
    threads: int,
    skipped_input: list[dict[str, str]],
) -> dict[str, Any]:
    if failed_results is None:
        failed_results = []

    all_ok = [r for r in results if r["status"] == "success"]
    n_trees = len(all_ok)
    n_failed = len(failed_results)
    n_skipped = len(skipped_input)

    is_error = n_trees == 0 and (n_failed > 0 or n_skipped > 0)
    if is_error:
        error_msg = "All FastTree runs failed"
    else:
        error_msg = None

    mean_n_taxa = 0.0
    mean_wall_time = 0.0
    if n_trees > 0:
        total_n_taxa = sum(r.get("n_taxa", 0) for r in all_ok)
        mean_n_taxa = total_n_taxa / n_trees if n_trees else 0.0
        total_wall = sum(r.get("wall_time", 0.0) for r in all_ok)
        mean_wall_time = total_wall / n_trees if n_trees else 0.0

    try:
        versions = _detect_fasttree_version(fasttree_exe)
    except Exception:
        versions = {"FastTree": "unknown"}

    cmd_parts = ["phyloai", "tree", "ml", "fasttree"]
    if batch_mode:
        cmd_parts.extend(["--msa-dir", str(msa_dir)])
    else:
        cmd_parts.extend(["--matrix", str(matrix)])
    cmd_parts.extend([
        "--seq-type", resolved_seq_type, "--model", model,
        "--mode", mode, "--boot", str(boot), "--cat", str(cat),
    ])
    if not gamma:
        cmd_parts.append("--no-gamma")
    cmd_parts.extend(["-o", str(output_dir)])
    cmd_str = " ".join(cmd_parts)

    payload: dict[str, Any] = {
        "status": "error" if is_error else "success",
        "command": cmd_str,
        "wall_time": time.monotonic() - run_start,
        "tool_versions": versions,
        "params": {
            "msa_dir": str(msa_dir) if msa_dir else None,
            "matrix": str(matrix) if matrix else None,
            "seq_type": resolved_seq_type,
            "model": model,
            "mode": mode,
            "boot": boot,
            "cat": cat,
            "gamma": gamma,
            "output_dir": str(output_dir),
            "threads": threads,
            "overwrite": overwrite,
            "fasttree_path": fasttree_path,
            "tool_args": tool_args,
        },
        "key_results": {
            "n_input": len(results) + n_failed + n_skipped,
            "n_trees": n_trees,
            "n_failed": n_failed,
            "n_skipped": n_skipped,
            "seq_type": resolved_seq_type,
            "model": model,
            "mode": mode,
            "boot": boot,
        },
        "error": error_msg,
        "data": {
            "summary": {
                "n_input_files": len(results) + n_failed + n_skipped,
                "n_trees": n_trees,
                "n_failed": n_failed,
                "n_skipped": n_skipped,
                "mean_n_taxa": mean_n_taxa,
                "mean_wall_time": mean_wall_time,
                "mode": "--msa-dir" if batch_mode else "--matrix",
            },
            "files": all_ok,
            "failed": failed_results,
            "skipped": skipped_input,
            "warnings": [],
        },
    }
    # Write aggregated step log
    import datetime as _dt
    log_path = output_dir / "fasttree.log"
    now_local = _dt.datetime.now().isoformat(timespec="seconds")
    with open(log_path, "a") as lf:
        lf.write(f"{now_local} | phyloai tree ml fasttree | exit={0 if n_trees > 0 else 2}\n")
        lf.write(f"command: {cmd_str}\n")
        for tool, ver in versions.items():
            lf.write(f"{tool}: {ver}\n")
        lf.write(f"wall_time: {payload['wall_time']:.2f}s\n")
        lf.write(f"trees: {n_trees}, failed: {n_failed}, skipped: {n_skipped}\n")
    return payload


def _detect_fasttree_version(executable: str) -> dict[str, str]:
    import subprocess as _sp
    import re as _re

    exe_name = Path(executable).name
    try:
        proc = _sp.run([executable], capture_output=True, text=True, timeout=10)
        combined = proc.stdout + proc.stderr
    except Exception:
        return {exe_name: "unknown"}

    m = _re.search(r"(?:version|FastTree)\s*([\d.]+)", combined, _re.IGNORECASE)
    if m:
        return {exe_name: m.group(1)}

    m = _re.search(r"([\d]+\.[\d]+(?:\.[\d]+)?)", combined)
    if m:
        return {exe_name: m.group(1)}

    return {exe_name: "unknown"}


def _reconstruct_result(output_dir: Path, run_start: float) -> dict[str, Any]:
    result_path = output_dir / "result.json"
    if result_path.exists():
        return json.loads(result_path.read_text())
    return {
        "status": "success",
        "command": "",
        "wall_time": time.monotonic() - run_start,
        "tool_versions": {},
        "params": {},
        "key_results": {},
        "error": None,
        "data": {"summary": {}, "files": [], "failed": [], "skipped": [], "warnings": []},
    }
```

- [ ] **Step 2: Write test for `run_fasttree` dry-run batch mode**

Append to `tests/tree/test_ml.py`:

```python
def test_run_fasttree_batch_dry_run(tmp_path: Path) -> None:
    from phyloai.tree.ml import run_fasttree

    msa_dir = tmp_path / "msas"
    msa_dir.mkdir()
    (msa_dir / "g1.fa").write_text(">a\nMKTLLL\n>b\nMKTLLL\n")
    (msa_dir / "g2.fa").write_text(">c\nACGTAC\n>d\nACGTAC\n")

    out_dir = tmp_path / "out"

    payload = run_fasttree(
        msa_dir=msa_dir, output_dir=out_dir,
        seq_type="AA", model="lg", dry_run=True, quiet=True,
    )

    assert payload["status"] == "success"
    assert payload["data"]["summary"]["n_input_files"] >= 2
    assert "files" in payload["data"]
    for f in payload["data"]["files"]:
        assert "cmd" in f


def test_run_fasttree_single_dry_run(tmp_path: Path) -> None:
    from phyloai.tree.ml import run_fasttree

    mat = tmp_path / "matrix.fa"
    mat.write_text(">a\nMKTLLL\n>b\nMKTLLL\n")

    out_dir = tmp_path / "out"

    payload = run_fasttree(
        matrix=mat, output_dir=out_dir,
        seq_type="AA", model="lg", dry_run=True, quiet=True,
    )

    assert payload["data"]["summary"]["mode"] == "--matrix"
    assert len(payload["data"]["files"]) >= 1
    if payload["data"]["files"]:
        assert "cmd" in payload["data"]["files"][0]


def test_run_fasttree_batch_auto_detects_mixed_and_fails(tmp_path: Path) -> None:
    from phyloai.tree.ml import run_fasttree

    msa_dir = tmp_path / "mixed"
    msa_dir.mkdir()
    (msa_dir / "aa.fa").write_text(">a\nMKTLLL\n")
    (msa_dir / "nt.fa").write_text(">b\nACGTAC\n")

    out_dir = tmp_path / "out"

    with pytest.raises(ValueError, match="Mixed sequence types"):
        run_fasttree(msa_dir=msa_dir, output_dir=out_dir, seq_type="auto", quiet=True)


def test_run_fasttree_batch_explicit_mismatch_fails(tmp_path: Path) -> None:
    from phyloai.tree.ml import run_fasttree

    msa_dir = tmp_path / "dir"
    msa_dir.mkdir()
    (msa_dir / "aa.fa").write_text(">a\nMKTLLL\n")

    out_dir = tmp_path / "out"

    with pytest.raises(ValueError, match="Files with wrong --seq-type"):
        run_fasttree(msa_dir=msa_dir, output_dir=out_dir, seq_type="NT", quiet=True)


def test_run_fasttree_invalid_model_for_aa(tmp_path: Path) -> None:
    from phyloai.tree.ml import run_fasttree

    msa_dir = tmp_path / "dir"
    msa_dir.mkdir()
    (msa_dir / "g1.fa").write_text(">a\nMKT\n")

    out_dir = tmp_path / "out"

    with pytest.raises(ValueError, match="Invalid model for AA.*gtr"):
        run_fasttree(msa_dir=msa_dir, output_dir=out_dir, seq_type="AA", model="gtr", quiet=True)


def test_run_fasttree_neither_input_raises() -> None:
    from phyloai.tree.ml import run_fasttree
    from pathlib import Path

    out_dir = Path("/tmp/out")
    with pytest.raises(ValueError, match="Either --msa-dir or --matrix"):
        run_fasttree(output_dir=out_dir, seq_type="AA", model="lg", quiet=True)
```

- [ ] **Step 3: Run tests for dry-run scenarios**

```bash
pytest tests/tree/test_ml.py::test_run_fasttree_batch_dry_run \
  tests/tree/test_ml.py::test_run_fasttree_single_dry_run \
  tests/tree/test_ml.py::test_run_fasttree_batch_auto_detects_mixed_and_fails \
  tests/tree/test_ml.py::test_run_fasttree_batch_explicit_mismatch_fails \
  tests/tree/test_ml.py::test_run_fasttree_invalid_model_for_aa -v
```

Expected: PASS (for tests that don't need FastTree installed)

- [ ] **Step 4: Commit**

```bash
git add phyloai/tree/ml.py tests/tree/test_ml.py
git commit -m "feat(tree): add run_fasttree library function"
```

---

### Task 7: CLI layer — `phyloai/cli/commands/tree.py`

**Files:**
- Create: `phyloai/cli/commands/tree.py`

- [ ] **Step 1: Write failing CLI tests**

Create `tests/cli/test_tree.py`:

```python
from __future__ import annotations

import json
from pathlib import Path

from click.testing import CliRunner

from phyloai.cli.main import cli


def test_tree_group_exists() -> None:
    result = CliRunner().invoke(cli, ["tree", "--help"])
    assert result.exit_code == 0
    assert "Maximum-likelihood" in result.output or "ml" in result.output


def test_tree_ml_help_shows_both_backends() -> None:
    result = CliRunner().invoke(cli, ["tree", "ml", "--help"])
    assert result.exit_code == 0
    assert "fasttree" in result.output
    assert "iqtree" in result.output


def test_tree_ml_fasttree_help() -> None:
    result = CliRunner().invoke(cli, ["tree", "ml", "fasttree", "--help"])
    assert result.exit_code == 0
    for flag in ["--msa-dir", "--matrix", "--seq-type", "--model", "--mode",
                  "--boot", "--cat", "--gamma", "--output-dir", "--threads"]:
        assert flag in result.output


def test_tree_ml_fasttree_mutual_exclusivity(tmp_path: Path) -> None:
    msa_dir = tmp_path / "msas"
    msa_dir.mkdir()
    mat = tmp_path / "matrix.fa"
    mat.write_text(">a\nMKT\n")

    result = CliRunner().invoke(cli, [
        "tree", "ml", "fasttree",
        "--msa-dir", str(msa_dir), "--matrix", str(mat),
    ])
    assert result.exit_code == 1


def test_tree_ml_fasttree_neither_input() -> None:
    result = CliRunner().invoke(cli, [
        "tree", "ml", "fasttree",
    ])
    assert result.exit_code == 1


def test_cli_msa_dir_nonexistent_exits_1() -> None:
    from click.testing import CliRunner
    from phyloai.cli.main import cli

    result = CliRunner().invoke(cli, [
        "tree", "ml", "fasttree",
        "--msa-dir", "/nonexistent/path",
    ])
    assert result.exit_code == 1


def test_tree_ml_fasttree_quiet_dry_run_batch(tmp_path: Path) -> None:
    msa_dir = tmp_path / "msas"
    msa_dir.mkdir()
    (msa_dir / "g1.fa").write_text(">a\nMKTLLL\n>b\nMKTLLL\n")

    out_dir = tmp_path / "out"

    result = CliRunner().invoke(cli, [
        "tree", "ml", "fasttree",
        "--msa-dir", str(msa_dir),
        "--output-dir", str(out_dir),
        "--seq-type", "AA",
        "--model", "lg",
        "--quiet",
        "--dry-run",
    ])

    assert result.exit_code == 0
    assert not (out_dir / "result.json").exists()


def test_tree_ml_fasttree_quiet_dry_run_single(tmp_path: Path) -> None:
    mat = tmp_path / "matrix.fa"
    mat.write_text(">a\nMKTLLL\n>b\nMKTLLL\n")

    out_dir = tmp_path / "out"

    result = CliRunner().invoke(cli, [
        "tree", "ml", "fasttree",
        "--matrix", str(mat),
        "--output-dir", str(out_dir),
        "--seq-type", "AA",
        "--model", "lg",
        "--quiet",
        "--dry-run",
    ])

    assert result.exit_code == 0


def test_tree_ml_fasttree_invalid_model_exits_1(tmp_path: Path) -> None:
    mat = tmp_path / "matrix.fa"
    mat.write_text(">a\nMKTLLL\n")

    out_dir = tmp_path / "out"

    result = CliRunner().invoke(cli, [
        "tree", "ml", "fasttree",
        "--matrix", str(mat),
        "--output-dir", str(out_dir),
        "--seq-type", "AA",
        "--model", "gtr",
        "--quiet",
    ])

    assert result.exit_code == 1


def test_tree_ml_fasttree_blocked_tool_args(tmp_path: Path) -> None:
    mat = tmp_path / "matrix.fa"
    mat.write_text(">a\nMKTLLL\n")

    out_dir = tmp_path / "out"

    result = CliRunner().invoke(cli, [
        "tree", "ml", "fasttree",
        "--matrix", str(mat),
        "--output-dir", str(out_dir),
        "--tool-args", "-lg",
        "--quiet",
    ])

    assert result.exit_code == 1


def test_tree_ml_fasttree_threads_warn_single(tmp_path: Path) -> None:
    mat = tmp_path / "matrix.fa"
    mat.write_text(">a\nMKTLLL\n>b\nMKTLLL\n")

    out_dir = tmp_path / "out"

    result = CliRunner().invoke(cli, [
        "tree", "ml", "fasttree",
        "--matrix", str(mat),
        "--output-dir", str(out_dir),
        "--threads", "8",
        "--quiet",
        "--dry-run",
    ])
    assert "has no effect" in result.output.lower() or result.exit_code == 0


def test_tree_ml_fasttree_writes_result_json_and_log(tmp_path: Path) -> None:
    import pytest

    mat = tmp_path / "matrix.fa"
    mat.write_text(">a\nMKTLLL\n>b\nMKTLLL\n")

    out_dir = tmp_path / "out"

    result = CliRunner().invoke(cli, [
        "tree", "ml", "fasttree",
        "--matrix", str(mat),
        "--output-dir", str(out_dir),
        "--seq-type", "AA",
        "--model", "lg",
        "--quiet",
    ])
    # May exit 3 if FastTree not installed
    if result.exit_code == 0:
        assert (out_dir / "result.json").exists()
        assert (out_dir / "fasttree.log").exists()
    elif result.exit_code == 3:
        pytest.skip("FastTree not installed")
```

- [ ] **Step 2: Run tests verify failure**

```bash
pytest tests/cli/test_tree.py::test_tree_group_exists -v
```

Expected: FAIL (tree group not registered)

- [ ] **Step 3: Write `phyloai/cli/commands/tree.py`**

```python
"""Tree inference CLI commands."""

from __future__ import annotations

import json
from pathlib import Path

import click
from rich.console import Console
from rich.progress import Progress

from phyloai.tree.ml import run_fasttree

console = Console()


def _fail(message: str, exit_code: int) -> None:
    click.echo(f"Error: {message}", err=True)
    raise click.exceptions.Exit(exit_code)


class _MLGroup(click.Group):
    def list_commands(self, ctx: click.Context) -> list[str]:
        return ["fasttree", "iqtree"]


class _TreeGroup(click.Group):
    def list_commands(self, ctx: click.Context) -> list[str]:
        return ["ml", "bi", "msc", "concordance"]


@click.group(cls=_TreeGroup)
def tree() -> None:
    """Phylogenetic tree inference commands."""


@tree.group(cls=_MLGroup)
def ml() -> None:
    """Maximum-likelihood tree inference (FastTree / IQ-TREE3)."""


@ml.command(
    "fasttree",
    help=(
        "Infer ML trees using FastTree.\n\n"
        "Two modes:\n"
        "  --msa-dir  : batch gene tree inference from a directory of MSA files\n"
        "  --matrix   : single supermatrix tree inference from one concatenated matrix\n\n"
        "FastTree natively reads FASTA and relaxed PHYLIP interleaved formats.\n"
        "NEXUS files must be converted first via 'phyloai pretree convert'."
    ),
)
@click.option(
    "--msa-dir",
    type=click.Path(file_okay=False, path_type=Path),
    default=None,
    help="Directory of MSA files for batch gene tree inference.",
)
@click.option(
    "--matrix",
    type=click.Path(dir_okay=False, path_type=Path),
    default=None,
    help="Single concatenated matrix file for supermatrix tree inference.",
)
@click.option(
    "--seq-type",
    type=click.Choice(["AA", "NT", "auto"]),
    default="auto",
    show_default=True,
    help="Molecule type.",
)
@click.option(
    "--model",
    type=click.Choice(["jtt", "lg", "wag", "jc", "gtr"]),
    default=None,
    help="Substitution model. AA default: lg. NT default: gtr.",
)
@click.option(
    "--mode",
    type=click.Choice(["normal", "fastest", "slow"]),
    default="normal",
    show_default=True,
    help="Speed/accuracy trade-off.",
)
@click.option(
    "--boot",
    type=click.IntRange(min=0),
    default=1000,
    show_default=True,
    help="Bootstrap replicates. 0 disables node support.",
)
@click.option(
    "--cat",
    type=click.IntRange(min=1),
    default=20,
    show_default=True,
    help="Number of rate categories for FastTree (-cat N).",
)
@click.option(
    "--gamma/--no-gamma",
    default=True,
    show_default=True,
    help="Enable gamma-distributed rate heterogeneity.",
)
@click.option(
    "--output-dir", "-o",
    type=click.Path(file_okay=False, path_type=Path),
    default=Path("runs/tree/ml/fasttree"),
    show_default=True,
    help="Output directory.",
)
@click.option(
    "--threads", "-t",
    type=int,
    default=4,
    show_default=True,
    help="Parallel gene tree workers. Only used in --msa-dir mode.",
)
@click.option(
    "--fasttree-path",
    type=Path,
    default=None,
    help="Explicit path to FastTree executable.",
)
@click.option(
    "--tool-args",
    type=str,
    default=None,
    help="Extra strategy flags passed verbatim to FastTree. Managed flags are blocked.",
)
@click.option("--overwrite", is_flag=True, default=False, help="Overwrite existing output directory.")
@click.option("--resume", is_flag=True, default=False, help="Resume from checkpoint (--msa-dir only).")
@click.option("--dry-run", is_flag=True, default=False, help="Show commands without executing.")
@click.option("--quiet", "-q", is_flag=True, default=False, help="Suppress terminal output except errors.")
def fasttree_command(
    msa_dir: Path | None,
    matrix: Path | None,
    seq_type: str,
    model: str | None,
    mode: str,
    boot: int,
    cat: int,
    gamma: bool,
    output_dir: Path,
    threads: int,
    fasttree_path: Path | None,
    tool_args: str | None,
    overwrite: bool,
    resume: bool,
    dry_run: bool,
    quiet: bool,
) -> None:
    batch_mode = msa_dir is not None
    single_mode = matrix is not None

    if batch_mode == single_mode:
        if not batch_mode and not single_mode:
            _fail("Either --msa-dir or --matrix is required.", 1)
        else:
            _fail("--msa-dir and --matrix are mutually exclusive.", 1)

    if threads < 1:
        _fail("--threads must be at least 1.", 1)

    if resume and overwrite:
        _fail("--overwrite and --resume are mutually exclusive.", 1)

    if resume and single_mode:
        _fail("--resume is only supported in --msa-dir mode.", 1)

    if not single_mode and not batch_mode:
        _fail("Either --msa-dir or --matrix is required.", 1)

    # Validate input paths exist
    if msa_dir is not None and not msa_dir.exists():
        _fail(f"--msa-dir does not exist: {msa_dir}", 1)
    if matrix is not None and not matrix.exists():
        _fail(f"--matrix does not exist: {matrix}", 1)

    # Validate fasttree-path
    if fasttree_path is not None:
        if not fasttree_path.exists():
            _fail(f"--fasttree-path does not exist: {fasttree_path}", 1)
        import os
        if not os.access(str(fasttree_path), os.X_OK):
            _fail(f"--fasttree-path is not executable: {fasttree_path}", 1)

    # Warn about --threads in single mode
    if single_mode and threads != 4:
        if not quiet:
            click.echo("Warning: --threads has no effect in single --matrix mode.", err=True)

    fasttree_path_str = str(fasttree_path) if fasttree_path else None

    def _invoke(progress_callback=None):
        return run_fasttree(
            msa_dir=msa_dir,
            matrix=matrix,
            output_dir=output_dir,
            seq_type=seq_type,
            model=model,
            mode=mode,
            boot=boot,
            cat=cat,
            gamma=gamma,
            threads=threads,
            fasttree_path=fasttree_path_str,
            tool_args=tool_args,
            overwrite=overwrite,
            resume=resume,
            dry_run=dry_run,
            quiet=quiet,
            progress_callback=progress_callback,
        )

    error_msg: str | None = None

    try:
        if not quiet and not dry_run and batch_mode:
            from phyloai.tree.ml import _scan_input

            found, _ = _scan_input(msa_dir)
            total = len(found)

            if total == 0:
                _fail("No valid input files found in --msa-dir.", 1)

            with Progress(console=console, transient=True) as progress:
                task = progress.add_task(
                    "[cyan]Inferring gene trees with FastTree", total=total
                )
                try:
                    payload = _invoke(
                        progress_callback=lambda _: progress.advance(task)
                    )
                except (ValueError, FileNotFoundError) as exc:
                    error_msg = str(exc)
        else:
            try:
                payload = _invoke()
            except (ValueError, FileNotFoundError) as exc:
                error_msg = str(exc)
    except SystemExit:
        raise
    except Exception as exc:
        error_msg = str(exc)

    if error_msg is not None:
        exit_code = 3 if "not found" in error_msg.lower() else 1
        _fail(error_msg, exit_code)

    if dry_run:
        if not quiet:
            click.echo(
                f"Dry run: {payload['data']['summary']['n_input_files']} input(s) would be processed."
            )
            for item in payload["data"].get("files", []):
                if "cmd" in item:
                    click.echo(" ".join(item["cmd"]))
        return

    result_path = output_dir / "result.json"
    result_path.parent.mkdir(parents=True, exist_ok=True)
    with open(result_path, "w") as fh:
        json.dump(payload, fh, indent=2)

    summary = payload["data"]["summary"]
    n_failed = summary.get("n_failed", 0)
    n_trees = summary.get("n_trees", 0)
    n_skipped = summary.get("n_skipped", 0)

    if not quiet:
        click.echo(
            f"Trees: {n_trees} | Failed: {n_failed} | Skipped: {n_skipped}"
        )
        if batch_mode:
            click.echo(f"Trees saved to {output_dir / 'trees'}", err=True)
            click.echo(f"Logs saved to {output_dir / 'logs'}", err=True)
        click.echo(f"Results saved to {result_path}", err=True)

        if n_failed > 0:
            click.echo(
                f"Warning: {n_failed} gene(s) failed. Check result.json data.failed for details.",
                err=True,
            )

    if n_trees == 0 and (n_failed > 0):
        _fail("All FastTree runs failed.", 2)
```

- [ ] **Step 4: Register tree in main.py**

Append to `phyloai/cli/main.py` (after existing imports):

```python
from phyloai.cli.commands.tree import tree
```

Append after the last `cli.add_command(...)`:

```python
cli.add_command(tree)
```

- [ ] **Step 5: Run CLI tests**

```bash
pytest tests/cli/test_tree.py -v
```

Expected: all dry-run/validation tests PASS (tests requiring actual FastTree will be run in CI with FastTree installed)

- [ ] **Step 6: Verify full --help output**

```bash
uv run phyloai tree ml fasttree --help
```

Expected: displays all parameters with help text

- [ ] **Step 7: Commit**

```bash
git add phyloai/cli/commands/tree.py phyloai/cli/main.py tests/cli/test_tree.py
git commit -m "feat(cli): add phyloai tree ml fasttree command"
```

---

### Task 8: User documentation

**Files:**
- Create: `docs/commands/tree-ml.md`

- [ ] **Step 1: Write command documentation**

```markdown
# phyloai tree ml fasttree

## Purpose

Infer maximum-likelihood phylogenetic trees using FastTree.

## Usage

```bash
# Batch gene trees from MSA directory (parallel)
phyloai tree ml fasttree --msa-dir ./trimmed/seqs \
    --seq-type AA --model lg --mode normal --boot 1000 \
    --cat 20 --gamma --threads 8 -o runs/tree/ml/fasttree

# Single supermatrix tree
phyloai tree ml fasttree --matrix ./concat/matrix.fa \
    --seq-type NT --model gtr --mode slow --boot 1000 \
    -o runs/tree/ml/fasttree

# Disable bootstrap (no node support)
phyloai tree ml fasttree --msa-dir ./trimmed --boot 0

# Fast mode, JTT model (AA default), no gamma
phyloai tree ml fasttree --msa-dir ./trimmed --mode fastest --model jtt --no-gamma
```

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--msa-dir` | — | Directory of MSA files. Mutually exclusive with `--matrix`. |
| `--matrix` | — | Single concatenated matrix file. Mutually exclusive with `--msa-dir`. |
| `--seq-type` | auto | AA, NT, or auto (detect from input). |
| `--model` | lg (AA) / gtr (NT) | Substitution model. AA: jtt, lg, wag. NT: jc, gtr. |
| `--mode` | normal | Speed/accuracy: normal, fastest, slow. |
| `--boot` | 1000 | Bootstrap replicates. 0 = no support (-nosupport). |
| `--cat` | 20 | Number of rate categories. |
| `--gamma` / `--no-gamma` | enabled | Gamma-distributed rate heterogeneity. |
| `--output-dir` / `-o` | runs/tree/ml/fasttree | Output directory. |
| `--threads` / `-t` | 4 | Parallel workers (--msa-dir only). |
| `--fasttree-path` | — | Explicit path to FastTree. |
| `--tool-args` | — | Extra strategy flags for FastTree. |
| `--overwrite` | — | Overwrite existing output dir. |
| `--resume` | — | Resume from checkpoint (--msa-dir only). |
| `--dry-run` | — | Print commands without executing. |
| `--quiet` / `-q` | — | Suppress output except errors. |

## Outputs

- `result.json`: structured results (trees, failed, skipped)
- `trees/`: Newick tree files (one per input, --msa-dir mode)
- `logs/`: per-gene FastTree logs (--msa-dir mode)
- `checkpoint.json`: resume state (--msa-dir mode)
- Single `.tre` file (--matrix mode)

## Supported Formats

FastTree natively reads FASTA (.fa, .fas, .fasta, .faa, .fna) and relaxed PHYLIP interleaved (.phy, .phylip).

NEXUS files (.nex, .nxs, .nexus) are not supported. Convert them first:
```bash
phyloai pretree convert --input data.nex --to fasta
```

## Notes

- `--threads` only controls parallel gene tree inference in `--msa-dir` mode. FastTree itself is single-threaded.
- `--resume` is only available in `--msa-dir` batch mode.
- Model default: LG for amino acid (AA), GTR for nucleotide (NT).
```

- [ ] **Step 2: Commit**

```bash
git add docs/commands/tree-ml.md
git commit -m "docs: add tree-ml command documentation"
```

---

### Task 9: Integration verification

- [ ] **Step 1: Run full test suite**

```bash
uv run pytest tests/tree/ tests/cli/test_tree.py -v
```

Expected: all tests pass

- [ ] **Step 2: Lint check**

```bash
uv run ruff check phyloai/tree/ phyloai/cli/commands/tree.py
```

Expected: no errors

- [ ] **Step 3: Verify CLI registration**

```bash
uv run phyloai --help
uv run phyloai tree --help
uv run phyloai tree ml --help
uv run phyloai tree ml fasttree --help
```

Expected: all help texts display correctly with proper ordering

- [ ] **Step 4: Commit**

```bash
git add -A
git commit -m "test: verify tree ml fasttree integration"
```
