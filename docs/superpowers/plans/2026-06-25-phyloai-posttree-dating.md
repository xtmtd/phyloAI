# PhyloAI posttree dating Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

> **Commits are user-driven only.** No task in this plan runs `git add` or `git commit`. Each Task ends with a "Report done — no commit" step. The user reviews the diff and stages/commits when they ask. This is intentional; do not add commit steps back without explicit user instruction.

**Goal:** Implement `phyloai posttree dating hessian` (IQ-TREE3 approximate likelihood) and `phyloai posttree dating mcmc` (MCMCtree Bayesian dating + full diagnostics), plus fix mcmctree version detection in `doctor`.

**Architecture:** Two library modules (`dating_hessian.py`, `dating_mcmc.py`) plus a diagnostics helper (`dating_diagnostics.py`) in `phyloai/posttree/`. CLI wiring added to the existing `posttree.py` as a new `dating` subgroup. MCMC progress tracked via `rich.Live` polling `mcmc.txt` line counts; four processes (2 posterior + 2 prior) run in parallel threads. All plots PDF via matplotlib.

**Tech Stack:** Python 3.10+, Click, rich, matplotlib, scipy.stats, subprocess, threading, pathlib. No new dependencies.

---

## File Map

| File | Action | Responsibility |
|------|--------|---------------|
| `phyloai/core/env.py` | Modify | Fix mcmctree version detection |
| `phyloai/posttree/dating_hessian.py` | Create | IQ-TREE hessian computation library layer; `--seq-type auto` reuses `core.formats.FormatConverter` and `core.sequence_normalization.detect_seq_type()` for FASTA/PHYLIP/NEXUS |
| `phyloai/posttree/dating_mcmc.py` | Create | MCMCtree MCMC run management library layer |
| `phyloai/posttree/dating_diagnostics.py` | Create | Time table parsing + all diagnostic plots |
| `phyloai/cli/commands/posttree.py` | Modify | Add `dating` subgroup + `hessian` + `mcmc` CLI commands |
| `tests/posttree/test_dating_hessian.py` | Create | Unit tests for hessian helpers |
| `tests/posttree/test_dating_mcmc.py` | Create | Unit tests for mcmc helpers |
| `tests/posttree/test_dating_diagnostics.py` | Create | Unit tests for diagnostics parsing + plots |
| `tests/core/test_env_mcmctree.py` | Create | Unit test for mcmctree version detection |

---

## Task 1: Fix mcmctree version detection in `env.py`

**Files:**
- Modify: `phyloai/core/env.py`
- Create: `tests/core/test_env_mcmctree.py`

MCMCtree prints `MCMCTREE in paml version 4.10.10, 27 Jan 2026` on the first
line of stdout when called with no arguments. Currently `version_flag = ""`
triggers the dir-scan fallback and never finds the version. Fix: set
`version_args = []` (empty list, not empty string) and add a custom
`version_pattern` key that `_get_version` honours.

- [ ] **Step 1: Write the failing test**

```python
# tests/core/test_env_mcmctree.py
"""Test mcmctree version extraction from stdout."""
from __future__ import annotations
import re
from unittest.mock import MagicMock, patch
from pathlib import Path

from phyloai.core.env import ToolEnv, TOOL_REGISTRY


def _make_result(stdout: str):
    r = MagicMock()
    r.stdout = stdout
    r.stderr = ""
    return r


def test_mcmctree_version_pattern_extracts_paml_version():
    """_get_version should match 'paml version X.Y.Z' preferentially."""
    env = ToolEnv()
    fake_output = "MCMCTREE in paml version 4.10.10, 27 Jan 2026\n"
    with patch("subprocess.run", return_value=_make_result(fake_output)):
        ver = env._get_version(Path("/usr/bin/mcmctree"), [])
    assert ver == "4.10.10"


def test_mcmctree_registry_uses_version_args_list():
    """TOOL_REGISTRY mcmctree entry must have version_args=[] not version_flag=''."""
    meta = TOOL_REGISTRY["mcmctree"]
    assert "version_args" in meta
    assert meta["version_args"] == []
    assert "version_flag" not in meta
```

- [ ] **Step 2: Run test to confirm it fails**

```bash
cd /Users/zf/data/coding/phyloAI
pytest tests/core/test_env_mcmctree.py -v
```
Expected: FAIL — `version_args` key missing from registry, pattern not matched.

- [ ] **Step 3: Update `TOOL_REGISTRY` entry for mcmctree**

In `phyloai/core/env.py`, replace:
```python
"mcmctree":   {"required": False, "version_flag": "",
               "install": "https://github.com/abacus-gene/paml/releases"},
```
with:
```python
"mcmctree":   {"required": False, "version_args": [],
               "version_pattern": r"paml version (\d+(?:\.\d+)+)",
               "install": "https://github.com/abacus-gene/paml/releases"},
```

- [ ] **Step 4: Update `_get_version` to honour `version_pattern`**

`_get_version` currently extracts the first version-like token from any line.
Add an optional `version_pattern` parameter that, when provided, is tried
first on the full output before falling back to the generic regex.

Modify the `_detect_tool` method signature and `_get_version` call site to
thread `version_pattern` through. Change in `_get_version`:

```python
def _get_version(
    self,
    path: Path,
    version_args: str | list[list[str]],
    version_pattern: str | None = None,
) -> Optional[str]:
    if version_args is None:
        return None
    if isinstance(version_args, str):
        if version_args == "":
            return self._version_from_tool_dir(path)
        candidates = [[version_args]]
    else:
        candidates = version_args
    try:
        for args in candidates:
            command = [str(path), *args]
            if path.suffix.lower() == ".jar":
                command = ["java", "-jar", str(path), *args]
            result = subprocess.run(
                command,
                capture_output=True, text=True, timeout=5
            )
            output = result.stdout.strip() or result.stderr.strip()
            # Try custom pattern first
            if version_pattern:
                m = re.search(version_pattern, output)
                if m:
                    return m.group(1)
            for line in output.splitlines():
                if line.strip():
                    version = self._normalize_version(line.strip()[:200])
                    if version:
                        return version
    except Exception:
        pass
    return None
```

Update `_detect_tool` to pass `version_pattern` through:

```python
def _detect_tool(self, name: str, version_flag: str = "",
                 version_args: Optional[list[list[str]]] = None,
                 version_pattern: Optional[str] = None,
                 bundled: bool = False,
                 bundled_dir: Optional[str] = None,
                 bundled_executable: Optional[str] = None,
                 path_aliases: Optional[list[str]] = None) -> ToolInfo:
    version_probe = version_args if version_args is not None else version_flag
    # ... existing path resolution unchanged ...
    # wherever _get_version is called, add version_pattern=version_pattern
```

Update `check_all` to pass `version_pattern`:
```python
info = self._detect_tool(
    name,
    version_flag=meta.get("version_flag", ""),
    version_args=meta.get("version_args"),
    version_pattern=meta.get("version_pattern"),
    bundled=meta.get("bundled", False),
    bundled_dir=meta.get("bundled_dir"),
    bundled_executable=meta.get("bundled_executable"),
    path_aliases=meta.get("path_aliases"),
)
```

- [ ] **Step 5: Run tests**

```bash
pytest tests/core/test_env_mcmctree.py -v
```
Expected: PASS (both tests).

- [ ] **Step 6: Report done — no commit (commits are user-driven; see plan header note)**

---

---

## Task 2: `dating_hessian.py` — pure helpers

**Files:**
- Create: `phyloai/posttree/dating_hessian.py`
- Create: `tests/posttree/test_dating_hessian.py`

Implement all pure (non-subprocess) helpers needed for the hessian step:
seq-type detection from FASTA/PHYLIP/NEXUS through shared core format helpers, partition count from a partition file,
IQ-TREE command builder, root-age constraint validator, and output file
validator.

- [ ] **Step 1: Write failing tests**

```python
# tests/posttree/test_dating_hessian.py
"""Tests for dating_hessian pure helpers."""
from __future__ import annotations
from pathlib import Path
import pytest
from phyloai.posttree.dating_hessian import (
    detect_seqtype_from_alignment,
    count_partitions,
    validate_root_age,
    build_iqtree_dating_cmd,
    HESSIAN_OUTPUT_FILES,
)


# ── detect_seqtype_from_alignment ────────────────────────────────────

def test_detect_aa(tmp_path):
    fa = tmp_path / "aa.fa"
    fa.write_text(">sp1\nMKTVFLGEI\n>sp2\nMLTVFLGEI\n")
    assert detect_seqtype_from_alignment(fa) == "AA"

def test_detect_nt(tmp_path):
    fa = tmp_path / "nt.fa"
    fa.write_text(">sp1\nACGTACGT\n>sp2\nACGTACGT\n")
    assert detect_seqtype_from_alignment(fa) == "NT"

def test_detect_auto_defaults_aa_for_mixed(tmp_path):
    """Sequences with non-ACGTN chars → AA."""
    fa = tmp_path / "m.fa"
    fa.write_text(">sp1\nACGTMKLWI\n")
    assert detect_seqtype_from_alignment(fa) == "AA"


# ── count_partitions ─────────────────────────────────────────────────

def test_count_raxml_partitions(tmp_path):
    pf = tmp_path / "parts.txt"
    pf.write_text("LG, p1 = 1-100\nLG, p2 = 101-200\nLG, p3 = 201-300\n")
    assert count_partitions(pf) == 3

def test_count_nexus_partitions(tmp_path):
    pf = tmp_path / "parts.nex"
    pf.write_text(
        "#NEXUS\nbegin sets;\n"
        "  charset p1 = 1-100;\n"
        "  charset p2 = 101-200;\n"
        "end;\n"
    )
    assert count_partitions(pf) == 2


# ── validate_root_age ────────────────────────────────────────────────

def test_valid_root_age_upper_only():
    tree = "(A,(B,C))'<4.2';"
    assert validate_root_age(tree) is True

def test_valid_root_age_range():
    tree = "(A,(B,C))'>3.1<4.2';"
    assert validate_root_age(tree) is True

def test_missing_root_age():
    tree = "(A,(B,C));"
    assert validate_root_age(tree) is False


# ── build_iqtree_dating_cmd ──────────────────────────────────────────

def test_unpartitioned_aa_default_model(tmp_path):
    matrix = tmp_path / "m.fa"
    matrix.touch()
    tree = tmp_path / "t.nwk"
    tree.touch()
    cmd = build_iqtree_dating_cmd(
        iqtree_path=Path("/usr/bin/iqtree3"),
        matrix=matrix,
        rooted_tree=tree,
        seq_type="AA",
        model_expr=None,
        partitions=None,
        n_partitions=0,
        prefix="iqtree",
        threads=4,
        tool_args=None,
    )
    assert "-m" in cmd
    idx = cmd.index("-m")
    assert cmd[idx + 1] == "LG+F+G4"
    assert "--dating" in cmd
    assert "mcmctree" in cmd

def test_unpartitioned_nt_default_model(tmp_path):
    matrix = tmp_path / "m.fa"
    matrix.touch()
    tree = tmp_path / "t.nwk"
    tree.touch()
    cmd = build_iqtree_dating_cmd(
        iqtree_path=Path("/usr/bin/iqtree3"),
        matrix=matrix,
        rooted_tree=tree,
        seq_type="NT",
        model_expr=None,
        partitions=None,
        n_partitions=0,
        prefix="iqtree",
        threads=4,
        tool_args=None,
    )
    idx = cmd.index("-m")
    assert cmd[idx + 1] == "GTR+G4"

def test_partitioned_aa_few(tmp_path):
    matrix = tmp_path / "m.fa"
    matrix.touch()
    tree = tmp_path / "t.nwk"
    tree.touch()
    parts = tmp_path / "parts.nex"
    parts.touch()
    cmd = build_iqtree_dating_cmd(
        iqtree_path=Path("/usr/bin/iqtree3"),
        matrix=matrix,
        rooted_tree=tree,
        seq_type="AA",
        model_expr=None,
        partitions=parts,
        n_partitions=5,
        prefix="iqtree",
        threads=4,
        tool_args=None,
    )
    assert "-Q" in cmd
    assert "--merge" not in cmd
    assert "--mset" in cmd
    idx = cmd.index("--mset")
    assert cmd[idx + 1] == "LG"

def test_partitioned_aa_many_merges(tmp_path):
    matrix = tmp_path / "m.fa"
    matrix.touch()
    tree = tmp_path / "t.nwk"
    tree.touch()
    parts = tmp_path / "parts.nex"
    parts.touch()
    cmd = build_iqtree_dating_cmd(
        iqtree_path=Path("/usr/bin/iqtree3"),
        matrix=matrix,
        rooted_tree=tree,
        seq_type="AA",
        model_expr=None,
        partitions=parts,
        n_partitions=12,
        prefix="iqtree",
        threads=4,
        tool_args=None,
    )
    assert "--merge" in cmd
    assert "--rclusterf" in cmd

def test_tool_args_appended_last(tmp_path):
    matrix = tmp_path / "m.fa"
    matrix.touch()
    tree = tmp_path / "t.nwk"
    tree.touch()
    cmd = build_iqtree_dating_cmd(
        iqtree_path=Path("/usr/bin/iqtree3"),
        matrix=matrix,
        rooted_tree=tree,
        seq_type="AA",
        model_expr=None,
        partitions=None,
        n_partitions=0,
        prefix="iqtree",
        threads=4,
        tool_args="--redo",
    )
    assert cmd[-1] == "--redo"

def test_hessian_output_files_constant():
    assert "iqtree.dummy.phy" in HESSIAN_OUTPUT_FILES
    assert "iqtree.rooted.nwk" in HESSIAN_OUTPUT_FILES
    assert "iqtree.mcmctree.hessian" in HESSIAN_OUTPUT_FILES


# ── T3: hessian output validation ─────────────────────────────────────

def test_run_hessian_fails_when_output_files_missing(tmp_path):
    """T3a: IQ-TREE produces no output files → status=error."""
    from unittest.mock import patch
    from phyloai.posttree.dating_hessian import run_hessian
    from phyloai.core.iqtree import _resolve_iqtree_path
    import subprocess

    matrix = tmp_path / "m.fa"
    matrix.write_text(">sp1\nMKTV\n>sp2\nMLTV\n")
    tree = tmp_path / "t.nwk"
    tree.write_text("(sp1,sp2)'<4.2';\n")
    output = tmp_path / "out"
    output.mkdir()

    def _fake_run(cmd, **kwargs):
        r = subprocess.CompletedProcess(cmd, returncode=0, stdout="", stderr="")
        return r

    with patch("subprocess.run", side_effect=_fake_run), \
         patch("phyloai.posttree.dating_hessian._resolve_iqtree_path",
               return_value="/usr/bin/iqtree3"), \
         patch("phyloai.posttree.dating_hessian._detect_iqtree_version",
               return_value={"iqtree3": "2.0.0"}):
        payload = run_hessian(
            matrix=matrix, rooted_tree=tree,
            output_dir=output, dry_run=False, quiet=True,
        )
    assert payload["status"] == "error"
    assert "Missing" in payload["error"] or "missing" in " ".join(map(str, payload.get("data", {}).get("warnings", []))) or "failed" in payload["error"]


def test_run_hessian_warns_on_empty_hessian_file(tmp_path):
    """T3b: subprocess succeeds but produces an empty hessian file →
    status=error with warning about empty files.
    """
    from unittest.mock import patch
    from phyloai.posttree.dating_hessian import run_hessian, HESSIAN_OUTPUT_FILES
    import subprocess

    matrix = tmp_path / "m.fa"
    matrix.write_text(">sp1\nMKTV\n>sp2\nMLTV\n")
    tree = tmp_path / "t.nwk"
    tree.write_text("(sp1,sp2)'<4.2';\n")
    output = tmp_path / "out"
    output.mkdir()

    def _fake_run(cmd, **kwargs):
        # Write all three files but make the hessian empty
        for fn in HESSIAN_OUTPUT_FILES:
            p = output / fn
            if fn == "iqtree.mcmctree.hessian":
                p.write_text("")  # empty
            else:
                p.write_text("x")
        return subprocess.CompletedProcess(cmd, returncode=0, stdout="", stderr="")

    with patch("subprocess.run", side_effect=_fake_run), \
         patch("phyloai.posttree.dating_hessian._resolve_iqtree_path",
               return_value="/usr/bin/iqtree3"), \
         patch("phyloai.posttree.dating_hessian._detect_iqtree_version",
               return_value={"iqtree3": "2.0.0"}):
        payload = run_hessian(
            matrix=matrix, rooted_tree=tree,
            output_dir=output, dry_run=False, quiet=True,
        )
    assert payload["status"] == "error"
    assert any("empty" in str(w).lower() for w in payload.get("data", {}).get("warnings", []))
```

- [ ] **Step 2: Run to confirm failure**

```bash
pytest tests/posttree/test_dating_hessian.py -v
```
Expected: ImportError (module doesn't exist yet).

- [ ] **Step 3: Implement `dating_hessian.py`**

```python
# phyloai/posttree/dating_hessian.py
"""IQ-TREE3 hessian computation for MCMCtree approximate likelihood dating."""
from __future__ import annotations

import re
import shlex
import subprocess
import time
from pathlib import Path
from typing import Any

from phyloai.core.formats import FormatConverter
from phyloai.core.iqtree import _resolve_iqtree_path, _detect_iqtree_version
from phyloai.core.sequence_normalization import detect_seq_type

HESSIAN_PREFIX = "iqtree"

HESSIAN_OUTPUT_FILES = (
    "iqtree.dummy.phy",
    "iqtree.rooted.nwk",
    "iqtree.mcmctree.hessian",
)

def detect_seqtype_from_alignment(matrix: Path) -> str:
    """Return 'AA' or 'NT' using PhyloAI's shared MSA reader and detector."""
    converter = FormatConverter()
    alignment = converter.read(matrix)
    return detect_seq_type([str(record.seq) for record in alignment])


def count_partitions(partition_file: Path) -> int:
    """Count partitions in a RAxML-like or NEXUS partition file."""
    text = partition_file.read_text(errors="ignore")
    if "charset" in text.lower():
        # NEXUS format: count charset lines
        return len(re.findall(r"(?i)\bcharset\b", text))
    # RAxML-like: count non-empty, non-comment lines
    lines = [
        l.strip() for l in text.splitlines()
        if l.strip() and not l.strip().startswith("#")
    ]
    return len(lines)


def validate_root_age(newick: str) -> bool:
    """Return True if outermost node has an age constraint label."""
    # Match last closing paren followed by a label with < or > constraint
    return bool(re.search(r"\)['\"]?[^;,()]*[<>][^;,()]*['\"]?\s*;", newick))


def build_iqtree_dating_cmd(
    *,
    iqtree_path: Path,
    matrix: Path,
    rooted_tree: Path,
    seq_type: str,
    model_expr: str | None,
    partitions: Path | None,
    n_partitions: int,
    threads: int,
    tool_args: str | None,
) -> list[str]:
    """Build the IQ-TREE3 command list for --dating mcmctree."""
    cmd: list[str] = [str(iqtree_path), "-s", str(matrix), "-te", str(rooted_tree),
                      "--dating", "mcmctree", "--prefix", HESSIAN_PREFIX, "-T", str(threads)]

    if partitions is None:
        # Unpartitioned mode
        model = model_expr or ("LG+F+G4" if seq_type == "AA" else "GTR+G4")
        cmd += ["-m", model]
    else:
        # Partitioned mode — model_expr ignored
        cmd += ["-m", "MF", "-Q", str(partitions)]
        if seq_type == "AA":
            cmd += ["--mset", "LG", "-mfreq", "F", "-mrate", "G"]
        else:
            cmd += ["--mset", "GTR", "-mrate", "G"]
        if n_partitions >= 10:
            cmd += ["--merge", "--rclusterf", "10"]

    if tool_args:
        cmd += shlex.split(tool_args)

    return cmd


def _validate_hessian_inputs(
    *,
    matrix: Path,
    rooted_tree: Path,
    model_expr: str | None,
    partitions: Path | None,
    threads: int,
    overwrite: bool,
    resume: bool,
    tool_args: str | None,
) -> list[str]:
    """Return list of validation error strings."""
    errors: list[str] = []
    import os as _os

    if not matrix.exists():
        errors.append(f"--matrix does not exist: {matrix}")
    elif not matrix.is_file():
        errors.append(f"--matrix is not a regular file: {matrix}")
    elif not _os.access(str(matrix), _os.R_OK):
        errors.append(f"--matrix is not readable: {matrix}")

    if not rooted_tree.exists():
        errors.append(f"--rooted-tree does not exist: {rooted_tree}")
    elif not rooted_tree.is_file():
        errors.append(f"--rooted-tree is not a regular file: {rooted_tree}")
    else:
        content = rooted_tree.read_text(errors="ignore")
        if not validate_root_age(content):
            errors.append(
                "--rooted-tree is missing a root age constraint. "
                "The outermost node must have a calibration label such as "
                "'<4.2' or '>3.1<4.2' (units: 100 Mya). "
                "Example: (A,(B,C))'<4.2';"
            )

    if model_expr and partitions:
        errors.append("--model-expr and --partitions are mutually exclusive.")

    if partitions:
        if not partitions.exists():
            errors.append(f"--partitions does not exist: {partitions}")
        elif not partitions.is_file():
            errors.append(f"--partitions is not a regular file: {partitions}")

    if threads < 1:
        errors.append(f"--threads must be >= 1, got {threads}")

    if overwrite and resume:
        errors.append("--overwrite and --resume are mutually exclusive.")

    if tool_args:
        # D1: expanded blocked set. `-s` (input matrix), `--dating`
        # (mcmctree mode), `-te` (input tree), `--prefix` (output file
        # naming) are all part of the hessian→mcmc contract. Letting
        # `--tool-args` override any of them silently would break the
        # mcmc step's ability to find the hessian outputs or use the
        # user's calibrated tree.
        blocked = {"-s", "--dating", "-te", "--prefix"}
        tokens = shlex.split(tool_args)
        for tok in tokens:
            if tok in blocked:
                errors.append(
                    f"--tool-args contains blocked flag '{tok}' "
                    f"(managed by PhyloAI). The hessian step must emit "
                    f"iqtree.dummy.phy, iqtree.rooted.nwk, and "
                    f"iqtree.mcmctree.hessian under the output "
                    f"directory using the supplied calibrated tree."
                )

    return errors


def run_hessian(
    *,
    matrix: Path,
    rooted_tree: Path,
    seq_type: str = "auto",
    model_expr: str | None = None,
    partitions: Path | None = None,
    output_dir: Path,
    threads: int = 4,
    iqtree_path: str | None = None,
    tool_args: str | None = None,
    overwrite: bool = False,
    resume: bool = False,
    dry_run: bool = False,
    quiet: bool = False,
    stream_output: bool = True,
) -> dict[str, Any]:
    """Library entry point for `phyloai posttree dating hessian`."""
    t0 = time.time()

    errors = _validate_hessian_inputs(
        matrix=matrix, rooted_tree=rooted_tree,
        model_expr=model_expr, partitions=partitions,
        threads=threads, overwrite=overwrite, resume=resume,
        tool_args=tool_args,
    )
    if errors:
        return {
            "status": "error",
            "command": "",
            "wall_time": 0.0,
            "tool_versions": {},
            "params": {},
            "key_results": {},
            "error": errors[0],
            "error_category": "input",
            "data": {"cmd": [], "tool_stderr": "", "warnings": errors},
        }

    # M5: resolve input paths to absolute. IQ-TREE runs with cwd=output_dir
    # so relative inputs would be resolved from there and missed.
    matrix = matrix.resolve()
    rooted_tree = rooted_tree.resolve()
    if partitions:
        partitions = partitions.resolve()

    # Seq-type detection
    if seq_type == "auto":
        seq_type = detect_seqtype_from_alignment(matrix)

    # Partition count
    n_partitions = 0
    if partitions:
        n_partitions = count_partitions(partitions)

    # Resolve IQ-TREE. M4: real signatures are
    # `_resolve_iqtree_path(iqtree_path, dry_run) -> str` and
    # `_detect_iqtree_version(executable) -> dict[str, str]`.
    iqtree_exe = _resolve_iqtree_path(iqtree_path, dry_run)
    tool_versions_iqtree = (
        _detect_iqtree_version(iqtree_exe) if not dry_run else {"iqtree3": "dry-run"}
    )

    cmd = build_iqtree_dating_cmd(
        iqtree_path=Path(iqtree_exe),
        matrix=matrix,
        rooted_tree=rooted_tree,
        seq_type=seq_type,
        model_expr=model_expr,
        partitions=partitions,
        n_partitions=n_partitions,
        threads=threads,
        tool_args=tool_args,
    )

    if dry_run:
        return {
            "status": "success",
            "command": " ".join(cmd),
            "wall_time": 0.0,
            "tool_versions": tool_versions_iqtree,
            "params": {"seq_type": seq_type, "n_partitions": n_partitions},
            "key_results": {},
            "error": None,
            "data": {"cmd": cmd, "tool_stderr": "", "warnings": []},
        }

    output_dir.mkdir(parents=True, exist_ok=True)

    proc = subprocess.run(
        cmd,
        cwd=output_dir,
        capture_output=(not stream_output),
        text=True,
    )

    # T3: validate hessian outputs — existence, non-empty, and IQ-TREE
    # report marker. Each failure mode produces a distinct warning so users
    # can distinguish "IQ-TREE crashed" from "files truncated".
    warnings: list[str] = []
    missing = [f for f in HESSIAN_OUTPUT_FILES if not (output_dir / f).exists()]
    if not missing:
        empty = [f for f in HESSIAN_OUTPUT_FILES
                 if (output_dir / f).stat().st_size == 0]
        if empty:
            missing = empty
            warnings.append(
                f"IQ-TREE produced empty output file(s): {empty}. "
                "Possible crash mid-write."
            )
    iqtree_report = output_dir / f"{prefix}.iqtree"
    if proc.returncode == 0 and iqtree_report.exists():
        report_text = iqtree_report.read_text(errors="ignore")
        if "Total CPU time used" not in report_text:
            warnings.append(
                f"{iqtree_report.name} has no 'Total CPU time used' marker — "
                "IQ-TREE may have been interrupted before completion."
            )

    if missing or proc.returncode != 0:
        return {
            "status": "error",
            "command": " ".join(cmd),
            "wall_time": time.time() - t0,
            "tool_versions": tool_versions_iqtree,
            "params": {"seq_type": seq_type, "n_partitions": n_partitions},
            "key_results": {},
            "error": f"IQ-TREE failed (returncode={proc.returncode}). Missing: {missing}",
            "error_category": "tool",
            "data": {"cmd": cmd, "tool_stderr": getattr(proc, "stderr", ""), "warnings": warnings},
        }

    wall = time.time() - t0
    return {
        "status": "success",
        "command": " ".join(cmd),
        "wall_time": wall,
        "tool_versions": tool_versions_iqtree,
        "params": {
            "seq_type": seq_type,
            "n_partitions": n_partitions,
            "model_expr": model_expr,
            "partitions": str(partitions) if partitions else None,
            "threads": threads,
            "prefix": prefix,
        },
        "key_results": {
            "hessian_file": str(output_dir / "iqtree.mcmctree.hessian"),
        },
        "error": None,
        "data": {
            "cmd": cmd,
            "tool_stderr": getattr(proc, "stderr", ""),
            "warnings": warnings,
            "output_files": {f: str(output_dir / f) for f in HESSIAN_OUTPUT_FILES},
        },
    }
```

- [ ] **Step 4: Run tests**

```bash
pytest tests/posttree/test_dating_hessian.py -v
```
Expected: All PASS.

- [ ] **Step 5: Report done — no commit**

---

## Task 3: `dating_diagnostics.py` — time table parsing + plots

**Files:**
- Create: `phyloai/posttree/dating_diagnostics.py`
- Create: `tests/posttree/test_dating_diagnostics.py`

Parse `mcmctree.out` for node time tables and `FigTree.node.tre`, compute
Spearman correlations, and generate all diagnostic plots.

- [ ] **Step 1: Write failing tests**

```python
# tests/posttree/test_dating_diagnostics.py
"""Tests for dating_diagnostics helpers."""
from __future__ import annotations
from pathlib import Path
import pytest
from phyloai.posttree.dating_diagnostics import (
    parse_mcmctree_out,
    extract_node_tree,
    build_time_table,
)


SAMPLE_OUT = """\
Posterior means and 95% Equal-tail CIs

t_n7          0.4213 (0.3521, 0.4891)
t_n8          0.3102 (0.2641, 0.3589)
t_n9          0.5521 (0.4822, 0.6198)

Species tree for FigTree.  Branch lengths = posterior mean times; 95% CIs = labels
((sp1:0.32,sp2:0.31) 7 :0.12,sp3:0.44) 8 ;

((sp1,sp2) 7 ,sp3) 8 ;

(sp1,sp2,sp3);
"""


def test_parse_mcmctree_out_extracts_times():
    rows = parse_mcmctree_out(SAMPLE_OUT)
    assert len(rows) == 3
    assert rows[0]["node"] == "t_n7"
    assert abs(rows[0]["mean"] - 0.4213) < 1e-4
    assert abs(rows[0]["lower"] - 0.3521) < 1e-4
    assert abs(rows[0]["upper"] - 0.4891) < 1e-4
    assert abs(rows[0]["ci_width"] - (0.4891 - 0.3521)) < 1e-4


def test_extract_node_tree_returns_first_tree_with_bare_integers():
    """D4: first of three FigTree trees with bare-integer internal labels."""
    tree = extract_node_tree(SAMPLE_OUT)
    assert tree is not None
    assert " 7 " in tree or ")7" in tree.replace(" ", "")


def test_build_time_table_two_runs():
    run1 = [
        {"node": "t_n7", "mean": 0.42, "lower": 0.35, "upper": 0.49, "ci_width": 0.14},
        {"node": "t_n8", "mean": 0.31, "lower": 0.26, "upper": 0.36, "ci_width": 0.10},
    ]
    run2 = [
        {"node": "t_n7", "mean": 0.43, "lower": 0.36, "upper": 0.50, "ci_width": 0.14},
        {"node": "t_n8", "mean": 0.30, "lower": 0.25, "upper": 0.35, "ci_width": 0.10},
    ]
    table = build_time_table(run1, run2)
    assert len(table) == 2
    assert "mean_run1" in table[0]
    assert "mean_run2" in table[0]
    assert "ci_width_run1" in table[0]
```

- [ ] **Step 2: Run to confirm failure**

```bash
pytest tests/posttree/test_dating_diagnostics.py -v
```
Expected: ImportError.

- [ ] **Step 3: Implement `dating_diagnostics.py`**

```python
# phyloai/posttree/dating_diagnostics.py
"""MCMCtree output parsing and diagnostic plot generation."""
from __future__ import annotations

import csv
import re
from pathlib import Path
from typing import Any

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import spearmanr


# ── Parsing ──────────────────────────────────────────────────────────

def parse_mcmctree_out(text: str) -> list[dict]:
    """Parse t_nXX lines from mcmctree.out text. Returns list of dicts."""
    rows = []
    for line in text.splitlines():
        m = re.match(
            r"^\s*(t_n\d+)\s+"
            r"([\d.]+)\s+\(([\d.]+),\s*([\d.]+)\)",
            line,
        )
        if m:
            node, mean, lower, upper = m.group(1), float(m.group(2)), float(m.group(3)), float(m.group(4))
            rows.append({
                "node": node,
                "mean": mean,
                "lower": lower,
                "upper": upper,
                "ci_width": upper - lower,
            })
    return rows


def extract_node_tree(text: str) -> str | None:
    """Extract the FigTree tree whose internal node labels are bare
    integers — this is the tree that maps `t_nXX` parameters in
    `mcmc.txt` / `mcmctree.out` to specific nodes in the species tree.

    D4: MCMCtree emits three trees after the
    `Species tree for FigTree.` header. In order:
      1. internal labels are integers (e.g. `(sp1, sp2) 7`); this is
         the node-label tree we want.
      2. internal labels are branch-length annotations
         (`&95%HPD=...`).
      3. internal labels are HPD intervals on the node labels
         themselves.
    Naive "first tree after marker" picks (2) in some PAML versions,
    which has no integer labels and so cannot be mapped to t_nXX
    parameters. We pick the first tree that contains at least one
    bare-integer internal label.
    """
    marker = "Species tree for FigTree."
    idx = text.find(marker)
    if idx == -1:
        return None
    after = text[idx:]
    for m in re.finditer(r"\([\s\S]+?\);", after):
        tree = m.group(0).strip()
        # A node-label tree has at least one internal-label position
        # holding just digits (e.g. `,sp2)7` or `)node7,` → both yield
        # a ` 7` token after a closing paren). Simplest heuristic: the
        # tree contains a digit-only token that follows `)` or `,`.
        if re.search(r"[),]\s*\d+\b", tree):
            return tree
    # Fallback to the first tree if none has integer labels
    first = re.search(r"(\([\s\S]+?;)", after)
    return first.group(1).strip() if first else None


def build_time_table(run1: list[dict], run2: list[dict]) -> list[dict]:
    """Merge two run time lists by node name into a combined table."""
    by_node = {r["node"]: r for r in run2}
    result = []
    for r in run1:
        node = r["node"]
        r2 = by_node.get(node, {})
        result.append({
            "node": node,
            "mean_run1": r["mean"],
            "lower_run1": r["lower"],
            "upper_run1": r["upper"],
            "ci_width_run1": r["ci_width"],
            "mean_run2": r2.get("mean", float("nan")),
            "lower_run2": r2.get("lower", float("nan")),
            "upper_run2": r2.get("upper", float("nan")),
            "ci_width_run2": r2.get("ci_width", float("nan")),
        })
    return result


def write_time_table_csv(table: list[dict], path: Path) -> None:
    if not table:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(table[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(table)


# ── Plots ─────────────────────────────────────────────────────────────

def _spearman(x: list[float], y: list[float]) -> tuple[float, float]:
    """Spearman rho + p-value. Returns (nan, nan) when < 3 valid points
    (D6: avoid scipy warnings on tiny samples).
    """
    if len(x) < 3:
        return float("nan"), float("nan")
    rho, pval = spearmanr(x, y)
    return float(rho), float(pval)


def _linear_fit(x: list[float], y: list[float]) -> tuple[float, float, float]:
    """Return (slope, intercept, rmse) for the least-squares fit. Returns
    (nan, nan, nan) when fewer than 2 valid points (D6: avoid polyfit
    warnings on tiny samples).
    """
    if len(x) < 2:
        return float("nan"), float("nan"), float("nan")
    m, b = np.polyfit(x, y, 1)
    residuals = [yi - (m * xi + b) for xi, yi in zip(x, y)]
    rmse = float(np.sqrt(np.mean([r * r for r in residuals])))
    return float(m), float(b), rmse


def plot_convergence(
    table: list[dict],
    x_col: str,
    y_col: str,
    out_path: Path,
    title: str,
    xlabel: str,
    ylabel: str,
) -> dict[str, float]:
    """Scatter + regression line + y=x reference (D3). Returns a dict
    with `rho`, `pvalue`, `slope`, `intercept`, `rmse`. Numeric values
    are `nan` when the sample is too small for the test (D6).
    """
    pairs = [
        (r[x_col], r[y_col])
        for r in table
        if not np.isnan(r.get(x_col, float("nan")))
        and not np.isnan(r.get(y_col, float("nan")))
    ]
    x = [p[0] for p in pairs]
    y = [p[1] for p in pairs]
    rho, pval = _spearman(x, y)
    slope, intercept, rmse = _linear_fit(x, y)

    fig, ax = plt.subplots(figsize=(5, 5))
    ax.scatter(x, y, s=40, alpha=0.8, color="steelblue", zorder=3)
    if len(x) >= 2:
        xs = np.linspace(min(min(x), min(y)), max(max(x), max(y)), 100)
        # y=x reference line (D3): perfect-convergence baseline.
        ax.plot(xs, xs, color="grey", lw=1.0, ls="--", alpha=0.7,
                label="y = x (reference)")
        m, b = np.polyfit(x, y, 1)
        ax.plot(xs, m * xs + b, color="firebrick", lw=1.5,
                label=f"fit: y = {m:.3f}x + {b:.3f}")
        ax.legend(fontsize=8)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)
    return {"rho": rho, "pvalue": pval, "slope": slope,
            "intercept": intercept, "rmse": rmse}


def plot_line(
    x: list[float],
    y: list[float],
    out_path: Path,
    title: str,
    xlabel: str,
    ylabel: str,
) -> dict[str, float]:
    """Line plot (points sorted by x and connected). Returns rho/pvalue/
    slope/intercept/rmse. D6: nan-safe when the sample is too small.
    """
    pairs = sorted(zip(x, y))
    xs = [p[0] for p in pairs]
    ys = [p[1] for p in pairs]
    rho, pval = _spearman(xs, ys)
    slope, intercept, rmse = _linear_fit(xs, ys)

    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(xs, ys, marker="o", markersize=4, color="steelblue", lw=1.2)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)
    return {"rho": rho, "pvalue": pval, "slope": slope,
            "intercept": intercept, "rmse": rmse}


def plot_trace(mcmc_txt: Path, out_path: Path, title: str) -> None:
    """Multi-panel trace plot from mcmc.txt."""
    if not mcmc_txt.exists():
        return
    lines = mcmc_txt.read_text(errors="ignore").splitlines()
    if len(lines) < 2:
        return
    header = lines[0].split()
    data = []
    for line in lines[1:]:
        try:
            data.append([float(v) for v in line.split()])
        except ValueError:
            continue
    if not data:
        return
    arr = np.array(data)
    # Plot first 12 columns max (node times + mu + sigma2 + lnL)
    n_cols = min(len(header), arr.shape[1], 12)
    cols = min(3, n_cols)
    rows = (n_cols + cols - 1) // cols
    fig, axes = plt.subplots(rows, cols, figsize=(cols * 3.5, rows * 2.5))
    axes_flat = np.array(axes).flatten() if n_cols > 1 else [axes]
    for i in range(n_cols):
        ax = axes_flat[i]
        ax.plot(arr[:, i], lw=0.8, color="steelblue")
        ax.set_title(header[i] if i < len(header) else f"col{i}", fontsize=8)
        ax.tick_params(labelsize=7)
    for i in range(n_cols, len(axes_flat)):
        axes_flat[i].set_visible(False)
    fig.suptitle(title, fontsize=10)
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)


# ── Top-level orchestration ───────────────────────────────────────────

def generate_all_diagnostics(
    *,
    run_dirs: list[Path],         # [run1_dir, run2_dir, ...]
    diag_dir: Path,
    n_runs: int = 2,
) -> dict[str, Any]:
    """Parse outputs and generate all diagnostics. Returns summary dict.

    D2/I7: convergence diagnostics are skipped when n_runs < 2, generated
    pairwise (every C(n,2) pair) when n_runs >= 2.
    D3/D6: every correlation row carries slope/intercept/rmse in addition
    to Spearman rho; rows for runs with < 3 valid points are emitted with
    empty numeric columns plus a warning.
    """
    import itertools
    warnings: list[str] = []
    corr_rows: list[dict] = []

    def _add_corr(comparison: str, stats: dict[str, float]) -> None:
        """Emit a correlation row. NaN stats become empty strings (D6)."""
        row = {"comparison": comparison}
        for k in ("rho", "pvalue", "slope", "intercept", "rmse"):
            v = stats.get(k, float("nan"))
            row[k] = "" if (isinstance(v, float) and np.isnan(v)) else v
        corr_rows.append(row)
        # D6: surface a warning for rows that couldn't be computed
        if row["rho"] == "" or row["slope"] == "":
            warnings.append(
                f"{comparison}: too few valid points (need >= 3 for "
                f"Spearman, >= 2 for slope/intercept/RMSE)."
            )

    # Parse posterior and prior times for each run
    post_times: list[list[dict]] = []
    prior_times: list[list[dict]] = []

    for i, run_dir in enumerate(run_dirs):
        run_label = f"run{i+1}"

        # Posterior
        out_file = run_dir / "mcmctree.out"
        if out_file.exists():
            rows = parse_mcmctree_out(out_file.read_text(errors="ignore"))
            post_times.append(rows)
            # D4: extract the node-label tree (first of three FigTree
            # trees whose internal labels are bare integers).
            node_tree = extract_node_tree(out_file.read_text(errors="ignore"))
            if node_tree:
                (run_dir / "FigTree.node.tre").write_text(node_tree + "\n")
        else:
            post_times.append([])

        # Prior
        prior_out = run_dir / "prior" / "mcmctree.out"
        if prior_out.exists():
            rows = parse_mcmctree_out(prior_out.read_text(errors="ignore"))
            prior_times.append(rows)
            node_tree = extract_node_tree(prior_out.read_text(errors="ignore"))
            if node_tree:
                (run_dir / "prior" / "FigTree.node.tre").write_text(node_tree + "\n")
        else:
            prior_times.append([])

        # Trace plots
        for kind, mcmc_file in [
            ("posterior", run_dir / "mcmc.txt"),
            ("prior", run_dir / "prior" / "mcmc.txt"),
        ]:
            plot_trace(
                mcmc_file,
                diag_dir / "traces" / f"mcmc_trace_{run_label}_{kind}.pdf",
                title=f"MCMC trace — {run_label} {kind}",
            )

    # I7: skip convergence entirely when only one run.
    if n_runs < 2:
        pass  # convergence section deliberately empty
    else:
        # D2: every pair (a, b) with a < b gets its own plot + table +
        # correlation row.
        pairs = list(itertools.combinations(range(n_runs), 2))
        for a, b in pairs:
            label_a = f"run{a+1}"
            label_b = f"run{b+1}"

            # Posterior convergence
            if (a < len(post_times) and b < len(post_times)
                    and post_times[a] and post_times[b]):
                table = build_time_table(post_times[a], post_times[b])
                # First pair writes the canonical posterior_times.txt
                # filename for backward compatibility with the legacy
                # diagnostic layout.
                table_path = (
                    diag_dir / "convergence" /
                    ("posterior_times.txt" if (a, b) == (0, 1)
                     else f"posterior_times_{label_a}_vs_{label_b}.txt")
                )
                write_time_table_csv(table, table_path)
                stats = plot_convergence(
                    table,
                    f"mean_{label_a}", f"mean_{label_b}",
                    diag_dir / "convergence" /
                    f"convergence_posterior_{label_a}_vs_{label_b}.pdf",
                    title=f"Convergence — posterior means ({label_a} vs {label_b})",
                    xlabel=f"Mean age {label_a} (100 Mya)",
                    ylabel=f"Mean age {label_b} (100 Mya)",
                )
                _add_corr(f"convergence_posterior_{label_a}_vs_{label_b}", stats)

            # Prior convergence
            if (a < len(prior_times) and b < len(prior_times)
                    and prior_times[a] and prior_times[b]):
                table = build_time_table(prior_times[a], prior_times[b])
                table_path = (
                    diag_dir / "convergence" /
                    ("prior_times.txt" if (a, b) == (0, 1)
                     else f"prior_times_{label_a}_vs_{label_b}.txt")
                )
                write_time_table_csv(table, table_path)
                stats = plot_convergence(
                    table,
                    f"mean_{label_a}", f"mean_{label_b}",
                    diag_dir / "convergence" /
                    f"convergence_prior_{label_a}_vs_{label_b}.pdf",
                    title=f"Convergence — prior means ({label_a} vs {label_b})",
                    xlabel=f"Mean age {label_a} (100 Mya)",
                    ylabel=f"Mean age {label_b} (100 Mya)",
                )
                _add_corr(f"convergence_prior_{label_a}_vs_{label_b}", stats)

    # Infinite-sites and posterior-vs-prior plots (per run, no pairwise)
    for i, run_dir in enumerate(run_dirs):
        run_label = f"run{i+1}"
        post = post_times[i] if i < len(post_times) else []
        prior = prior_times[i] if i < len(prior_times) else []

        for kind, rows in [("posterior", post), ("prior", prior)]:
            if not rows:
                continue
            x = [r["mean"] for r in rows]
            y = [r["ci_width"] for r in rows]
            stats = plot_line(
                x, y,
                diag_dir / "infinite_sites" /
                f"infinite_sites_{run_label}_{kind}.pdf",
                title=f"Infinite-sites — {run_label} {kind}",
                xlabel="Mean age (100 Mya)",
                ylabel="95% CI width (100 Mya)",
            )
            _add_corr(f"infinite_sites_{run_label}_{kind}", stats)

        if post and prior:
            post_by_node = {r["node"]: r["mean"] for r in post}
            prior_by_node = {r["node"]: r["mean"] for r in prior}
            shared = [n for n in post_by_node if n in prior_by_node]
            if shared:
                xp = [post_by_node[n] for n in shared]
                yp = [prior_by_node[n] for n in shared]
                stats = plot_line(
                    xp, yp,
                    diag_dir / "posterior_vs_prior" /
                    f"posterior_vs_prior_{run_label}.pdf",
                    title=f"Posterior vs prior — {run_label}",
                    xlabel="Posterior mean age (100 Mya)",
                    ylabel="Prior mean age (100 Mya)",
                )
                _add_corr(f"posterior_vs_prior_{run_label}", stats)

    # Write Spearman CSV (D3: extended columns)
    if corr_rows:
        corr_path = diag_dir / "spearman_correlations.csv"
        corr_path.parent.mkdir(parents=True, exist_ok=True)
        with open(corr_path, "w", newline="") as f:
            writer = csv.DictWriter(
                f,
                fieldnames=["comparison", "rho", "pvalue",
                            "slope", "intercept", "rmse"],
            )
            writer.writeheader()
            writer.writerows(corr_rows)

    return {"spearman": corr_rows, "warnings": warnings}
```

- [ ] **Step 4: Run tests**

```bash
pytest tests/posttree/test_dating_diagnostics.py -v
```
Expected: All PASS.

- [ ] **Step 5: Report done — no commit**

---

## Task 4: `dating_mcmc.py` — MCMCtree run management

**Files:**
- Create: `phyloai/posttree/dating_mcmc.py`
- Create: `tests/posttree/test_dating_mcmc.py`

Implement ctl generation, seqtype inference from dummy.phy, ndata counting,
run directory setup, process launch with progress, and diagnostics trigger.

- [ ] **Step 1: Write failing tests**

```python
# tests/posttree/test_dating_mcmc.py
"""Tests for dating_mcmc pure helpers."""
from __future__ import annotations
from pathlib import Path
import pytest
from phyloai.posttree.dating_mcmc import (
    detect_seqtype_from_phylip,
    count_ndata_from_phylip,
    generate_mcmctree_ctl,
    validate_hessian_dir,
    count_mcmc_samples,
    _resolve_seqtype_and_ndata,
    _derive_prior_ctl,
)


DUMMY_PHY_AA = " 3 50\nsp1  MKTV\nsp2  MLTV\nsp3  MSTV\n"
DUMMY_PHY_NT = " 3 50\nsp1  ACGT\nsp2  ACGT\nsp3  ACGT\n"
DUMMY_PHY_2PART = " 3 50\nsp1  MKTV\nsp2  MLTV\nsp3  MSTV\n\n 3 50\nsp1  MKTV\nsp2  MLTV\nsp3  MSTV\n"


def test_detect_seqtype_aa(tmp_path):
    p = tmp_path / "dummy.phy"
    p.write_text(DUMMY_PHY_AA)
    assert detect_seqtype_from_phylip(p) == "AA"

def test_detect_seqtype_nt(tmp_path):
    p = tmp_path / "dummy.phy"
    p.write_text(DUMMY_PHY_NT)
    assert detect_seqtype_from_phylip(p) == "NT"

def test_count_ndata_single(tmp_path):
    p = tmp_path / "dummy.phy"
    p.write_text(DUMMY_PHY_AA)
    assert count_ndata_from_phylip(p) == 1

def test_count_ndata_two_partitions(tmp_path):
    p = tmp_path / "dummy.phy"
    p.write_text(DUMMY_PHY_2PART)
    assert count_ndata_from_phylip(p) == 2


def test_generate_mcmctree_ctl_posterior(tmp_path):
    ctl = generate_mcmctree_ctl(
        seqtype_code=2,
        ndata=1,
        clock=2,
        burnin=100000,
        sampfreq=10,
        nsample=10000,
        usedata=2,
        seed=-1,
    )
    assert "seqtype = 2" in ctl
    assert "usedata = 2" in ctl
    assert "seed = -1" in ctl
    assert "burnin = 100000" in ctl
    assert "ndata = 1" in ctl

def test_generate_mcmctree_ctl_prior(tmp_path):
    ctl = generate_mcmctree_ctl(
        seqtype_code=2,
        ndata=1,
        clock=2,
        burnin=100000,
        sampfreq=10,
        nsample=10000,
        usedata=0,
        seed=42,
    )
    assert "usedata = 0" in ctl
    assert "seed = 42" in ctl


def test_validate_hessian_dir_ok(tmp_path):
    for f in ("iqtree.dummy.phy", "iqtree.rooted.nwk", "iqtree.mcmctree.hessian"):
        (tmp_path / f).write_text("x")
    errs = validate_hessian_dir(tmp_path)
    assert errs == []

def test_validate_hessian_dir_missing(tmp_path):
    (tmp_path / "iqtree.dummy.phy").write_text("x")
    errs = validate_hessian_dir(tmp_path)
    assert len(errs) > 0


def test_count_mcmc_samples_empty(tmp_path):
    p = tmp_path / "mcmc.txt"
    assert count_mcmc_samples(p) == 0

def test_count_mcmc_samples_with_header(tmp_path):
    p = tmp_path / "mcmc.txt"
    p.write_text("Gen\tt_n7\tmu\n1\t0.4\t0.01\n2\t0.41\t0.011\n")
    assert count_mcmc_samples(p) == 2


# ── _resolve_seqtype_and_ndata (M1) ──────────────────────────────────

def test_resolve_seqtype_and_ndata_prefers_hessian_result_json(tmp_path):
    import json
    (tmp_path / "result.json").write_text(json.dumps({
        "params": {"seq_type": "AA", "n_partitions": 3},
        "key_results": {},
    }))
    # even if dummy.phy says NT, the result.json wins
    (tmp_path / "iqtree.dummy.phy").write_text(DUMMY_PHY_NT)
    seq, n, src = _resolve_seqtype_and_ndata(tmp_path)
    assert seq == "AA"
    assert n == 3
    assert src == "hessian-result.json"


def test_resolve_seqtype_and_ndata_fallback_to_dummy_phy(tmp_path):
    # no result.json → parse dummy.phy
    (tmp_path / "iqtree.dummy.phy").write_text(DUMMY_PHY_AA)
    seq, n, src = _resolve_seqtype_and_ndata(tmp_path)
    assert seq == "AA"
    assert n == 1
    assert src == "dummy.phy-fallback"


def test_resolve_seqtype_and_ndata_missing_both_raises(tmp_path):
    import pytest
    with pytest.raises(FileNotFoundError):
        _resolve_seqtype_and_ndata(tmp_path)


def test_resolve_seqtype_and_ndata_partial_result_json_falls_back(tmp_path):
    # result.json exists but lacks required fields → fallback
    import json
    (tmp_path / "result.json").write_text(json.dumps({"params": {}}))
    (tmp_path / "iqtree.dummy.phy").write_text(DUMMY_PHY_NT)
    seq, n, src = _resolve_seqtype_and_ndata(tmp_path)
    assert seq == "NT"
    assert n == 1
    assert src == "dummy.phy-fallback"


# ── _derive_prior_ctl (M2) ───────────────────────────────────────────

def test_derive_prior_ctl_substitutes_usedata_and_seed():
    posterior_ctl = (
        "      seed = -1\n"
        "      usedata = 2    * 0: no data; 1:seq like; 2:use in.BV; 3: out.BV\n"
        "      model = 0    * 0:JC69, 1:K80, 2:F81, 3:F84, 4:HKY85\n"
    )
    prior = _derive_prior_ctl(posterior_ctl, seed=12345)
    assert "usedata = 0" in prior
    assert "usedata = 2" not in prior
    assert "seed = 12345" in prior
    assert "seed = -1" not in prior
    # Other parameters preserved verbatim
    assert "model = 0" in prior
    assert "* 0:JC69, 1:K80" in prior


def test_derive_prior_ctl_appends_missing_lines():
    # ctl with neither usedata= nor seed= lines
    ctl_text = "      ndata = 1\n      clock = 2\n"
    prior = _derive_prior_ctl(ctl_text, seed=42)
    assert "usedata = 0" in prior
    assert "seed = 42" in prior
    assert "ndata = 1" in prior
    assert "clock = 2" in prior


# ── T1: fake-mcmctree integration test ───────────────────────────────

def test_run_mcmc_with_fake_mcmctree_exits_successfully(tmp_path):
    """T1: launcher that writes SeedUsed + mcmc.txt + mcmctree.out quickly,
    mimicking two posterior + two prior processes. Asserts run_mcmc returns
    status=success with both convergence and diagnostics populated.
    """
    import json, stat as _stat, sys

    # 1. Fake hessian directory with required files + result.json.
    hessian = tmp_path / "hessian"
    hessian.mkdir()
    for fname in ("iqtree.dummy.phy", "iqtree.rooted.nwk", "iqtree.mcmctree.hessian"):
        (hessian / fname).write_text("x")
    (hessian / "result.json").write_text(json.dumps({
        "params": {"seq_type": "AA", "n_partitions": 1},
        "key_results": {},
    }))

    # 2. Fake mcmctree script. Args: <ctl file>. The script writes
    #    SeedUsed immediately, then writes nsamples lines to mcmc.txt,
    #    then writes a minimal mcmctree.out.
    fake_mcmctree = tmp_path / "fake_mcmctree"
    fake_mcmctree.write_text(
        r"""#!/usr/bin/env python3
import shutil, sys, time, os
ctl_file = sys.argv[1]
run_dir = os.path.dirname(os.path.abspath(ctl_file))

# Write SeedUsed immediately so prior can start.
seed_path = os.path.join(run_dir, "SeedUsed")
with open(seed_path, "w") as f:
    f.write("42\n")

# Write mcmc.txt with nsamples lines extracted from the ctl.
with open(ctl_file) as f:
    ctl_text = f.read()
import re
match = re.search(r'nsample\s*=\s*(\d+)', ctl_text)
nsample = int(match.group(1)) if match else 2

with open(os.path.join(run_dir, "mcmc.txt"), "w") as f:
    f.write("Gen\tt_n7\tmu\tsigma2\tlnL\n")
    for i in range(1, nsample + 1):
        f.write(f"{i}\t0.42\t0.01\t0.005\t-10.0\n")

# Write minimal mcmctree.out with one node time.
with open(os.path.join(run_dir, "mcmctree.out"), "w") as f:
    f.write("Posterior means and 95% Equal-tail CIs\n")
    f.write("t_n7  0.4213 (0.3521, 0.4891)\n")
    f.write("\n")
    f.write("Species tree for FigTree.  Branch lengths = posterior mean times; 95% CIs = labels\n")
    f.write("(sp1,sp2) 7 ;\n")
    f.write("(sp1,sp2) 7 ;\n")
    f.write("(sp1,sp2);\n")
"""
    )
    fake_mcmctree.chmod(fake_mcmctree.stat().st_mode | _stat.S_IXUSR | _stat.S_IXGRP | _stat.S_IXOTH)

    # 3. Call run_mcmc with a short run.
    from phyloai.posttree.dating_mcmc import run_mcmc
    output = tmp_path / "out"
    payload = run_mcmc(
        hessian_dir=hessian,
        mcmctree_path=str(fake_mcmctree),
        nsamples=3,
        burnin=0,
        sample_freq=1,
        n_runs=2,
        output_dir=output,
        quiet=True,
        seed_wait_timeout_sec=30,
    )
    assert payload["status"] == "success", f"Expected success, got error: {payload.get('error')}"
    assert payload["key_results"]["n_runs"] == 2
    assert payload["key_results"]["n_posterior_failures"] == 0
    assert (output / "diagnostics" / "convergence" / "convergence_posterior_run1_vs_run2.pdf").exists()
    assert (output / "diagnostics" / "spearman_correlations.csv").exists()
    # T1d: prior should have been launched (SeedUsed was written quickly)
    assert payload["key_results"]["n_priors_skipped_timeout"] == 0


def test_run_mcmc_fake_mcmctree_posterior_failure_is_error(tmp_path):
    """T1b: fake-mcmctree that exits non-zero is treated as error."""
    import json, stat as _stat

    hessian = tmp_path / "hessian"
    hessian.mkdir()
    for fname in ("iqtree.dummy.phy", "iqtree.rooted.nwk", "iqtree.mcmctree.hessian"):
        (hessian / fname).write_text("x")
    (hessian / "result.json").write_text(json.dumps({
        "params": {"seq_type": "AA", "n_partitions": 1},
    }))

    fake_mcmctree = tmp_path / "fake_fail"
    fake_mcmctree.write_text(
        "#!/usr/bin/env python3\nimport sys\nsys.exit(1)\n"
    )
    fake_mcmctree.chmod(fake_mcmctree.stat().st_mode | _stat.S_IXUSR | _stat.S_IXGRP | _stat.S_IXOTH)

    from phyloai.posttree.dating_mcmc import run_mcmc
    output = tmp_path / "out"
    payload = run_mcmc(
        hessian_dir=hessian,
        mcmctree_path=str(fake_mcmctree),
        nsamples=3, burnin=0, sample_freq=1, n_runs=2,
        output_dir=output, quiet=True,
        seed_wait_timeout_sec=5,
    )
    assert payload["status"] == "error"
    assert "returncode" in str(payload.get("data", {})).lower() or any(
        "exited with code" in w for w in payload["data"].get("warnings", []))


# ── _SampleCounter (I3) ───────────────────────────────────────────────

def test_sample_counter_handles_growing_file(tmp_path):
    from phyloai.posttree.dating_mcmc import _SampleCounter
    f = tmp_path / "mcmc.txt"
    f.write_text("Gen\tt_n7\tmu\n")
    c = _SampleCounter()
    assert c.count(f) == 0
    with open(f, "a") as fh:
        fh.write("1\t0.4\t0.01\n2\t0.41\t0.011\n")
    assert c.count(f) == 2
    with open(f, "a") as fh:
        fh.write("3\t0.42\t0.012\n")
    assert c.count(f) == 3


def test_sample_counter_skips_partial_trailing_line(tmp_path):
    from phyloai.posttree.dating_mcmc import _SampleCounter
    f = tmp_path / "mcmc.txt"
    f.write_text("Gen\tt_n7\tmu\n1\t0.4\t0.01\n2\t0.")
    c = _SampleCounter()
    # The truncated line "2\t0." is non-numeric → skipped
    assert c.count(f) == 1


def test_sample_counter_recovers_from_file_replacement(tmp_path):
    from phyloai.posttree.dating_mcmc import _SampleCounter
    f = tmp_path / "mcmc.txt"
    f.write_text("Gen\tt_n7\tmu\n1\t0.4\t0.01\n")
    c = _SampleCounter()
    assert c.count(f) == 1
    # Simulate rotation (delete + recreate → new inode)
    f.unlink()
    f.write_text("Gen\tt_n7\tmu\n")
    assert c.count(f) == 0  # reset
    with open(f, "a") as fh:
        fh.write("10\t0.5\t0.02\n")
    assert c.count(f) == 1
```

- [ ] **Step 2: Run to confirm failure**

```bash
pytest tests/posttree/test_dating_mcmc.py -v
```
Expected: ImportError.

- [ ] **Step 3: Implement `dating_mcmc.py`**

```python
# phyloai/posttree/dating_mcmc.py
"""MCMCtree Bayesian dating run management."""
from __future__ import annotations

import shutil
import subprocess
import threading
import time
from pathlib import Path
from typing import Any

from rich.console import Console
from rich.progress import (
    BarColumn,
    Progress,
    TaskProgressColumn,
    TextColumn,
    TimeElapsedColumn,
)
from rich.live import Live

from phyloai.posttree.dating_hessian import HESSIAN_OUTPUT_FILES

console = Console()

SEQTYPE_CODE = {"AA": 2, "NT": 0}


# ── Pure helpers ──────────────────────────────────────────────────────

def detect_seqtype_from_phylip(phylip_path: Path) -> str:
    """Detect AA or NT from iqtree.dummy.phy content."""
    text = phylip_path.read_text(errors="ignore").upper()
    lines = [l.strip() for l in text.splitlines() if l.strip()]
    # Collect sequence chars (skip header line and taxon name tokens)
    aa_only = set("ARNDCQEGHILKMFPSTWYV")
    nt_chars = set("ACGTURYSWKMBDHVN-?")
    for line in lines[1:]:
        parts = line.split(None, 1)
        seq = parts[1] if len(parts) > 1 else parts[0]
        seq = seq.replace(" ", "").replace("-", "").replace("?", "")
        extra = set(seq) - nt_chars
        if extra & aa_only:
            return "AA"
    return "NT"


def count_ndata_from_phylip(phylip_path: Path) -> int:
    """Count number of data blocks (partitions) in a PAML phylip file."""
    text = phylip_path.read_text(errors="ignore")
    # Each block starts with a header line matching '  N  L'
    import re
    blocks = re.findall(r"^\s+\d+\s+\d+\s*$", text, re.MULTILINE)
    return max(1, len(blocks))


def validate_hessian_dir(hessian_dir: Path) -> list[str]:
    """Return list of errors if required hessian files are missing."""
    errors: list[str] = []
    for fname in HESSIAN_OUTPUT_FILES:
        p = hessian_dir / fname
        if not p.exists():
            errors.append(f"Missing required file in --hessian-dir: {fname}")
    return errors


def generate_mcmctree_ctl(
    *,
    seqtype_code: int,
    ndata: int,
    clock: int,
    burnin: int,
    sampfreq: int,
    nsample: int,
    usedata: int,
    seed: int,
) -> str:
    """Return a mcmctree.ctl string with full inline comments matching the
    upstream mcmctree.sh reference. Comments on parameters that approximate
    likelihood (`usedata = 2`) ignores (model, alpha, ncatG, cleandata,
    kappa_gamma, alpha_gamma) are preserved so users can edit the ctl
    confidently and so the file looks familiar to anyone migrating from
    mcmctree.sh.
    """
    return f"""\
          seed = {seed}
       seqfile = iqtree.dummy.phy
      treefile = iqtree.rooted.nwk
       outfile = mcmctree.out

         ndata = {ndata}
       seqtype = {seqtype_code}  * 0: nucleotides; 1:codons; 2:AAs
       usedata = {usedata}    * 0: no data; 1:seq like; 2:use in.BV; 3: out.BV
         clock = {clock}    * 1: global clock; 2: independent rates; 3: correlated rates
       RootAge =   * safe constraint on root age, used if no fossil for root.

       BDparas = 1 1 0.1 M   * birth, death, sampling
   rgene_gamma = 2 20 1   * gamma prior for overall rates for genes
  sigma2_gamma = 1 10 1    * gamma prior for sigma^2     (for clock=2 or 3)

      finetune = 0: .1  .1  .1  .1 .1 .1  * auto (0 or 1) : times, musigma2, rates, mixing, paras, FossilErr

*** These parameters control the MCMC run
***  Note: Total number of MCMC iterations will be burnin + (sampfreq * nsample)

         print = 1
        burnin = {burnin}
      sampfreq = {sampfreq}
       nsample = {nsample}


*** The following parameters only needed to run MCMCtree with exact likelihood (usedata = 1)
*** no need to change anything for approximate likelihood (usedata = 2)

         model = 0    * 0:JC69, 1:K80, 2:F81, 3:F84, 4:HKY85
         alpha = 0.5    * alpha for gamma rates at sites
         ncatG = 4    * No. categories in discrete gamma

     cleandata = 0  * remove sites with ambiguity data (1:yes, 0:no)?

   kappa_gamma = 6 2      * gamma prior for kappa
   alpha_gamma = 1 1      * gamma prior for alpha

*** Note: Make your window wider (100 columns) before running the program.
"""


def count_mcmc_samples(mcmc_txt: Path) -> int:
    """Count completed samples in mcmc.txt (lines excluding header)."""
    if not mcmc_txt.exists():
        return 0
    try:
        lines = mcmc_txt.read_bytes().splitlines()
        data_lines = [l for l in lines[1:] if l.strip()]
        return len(data_lines)
    except Exception:
        return 0


def _resolve_seqtype_and_ndata(
    hessian_dir: Path,
) -> tuple[str, int, str]:
    """Read seq_type ("AA"|"NT") and ndata (int) for MCMCtree ctl generation.

    M1 contract: prefer values recorded by the `hessian` step in its
    `result.json` (`params.seq_type`, `params.n_partitions`). Fall back to
    parsing `iqtree.dummy.phy` when those fields are absent (e.g. when the
    hessian directory was produced by an external tool or an older PhyloAI
    version).

    Returns (seq_type_str, ndata, source) where source is
    `"hessian-result.json"` or `"dummy.phy-fallback"`.
    """
    import json
    result_json = hessian_dir / "result.json"
    if result_json.exists():
        try:
            data = json.loads(result_json.read_text(errors="ignore"))
            params = data.get("params") or {}
            seq_type = params.get("seq_type")
            n_part = params.get("n_partitions")
            if seq_type in ("AA", "NT") and isinstance(n_part, int) and n_part >= 1:
                return seq_type, n_part, "hessian-result.json"
        except Exception:
            pass

    dummy_phy = hessian_dir / "iqtree.dummy.phy"
    if not dummy_phy.exists():
        raise FileNotFoundError(
            f"Neither {result_json} nor {dummy_phy} provides seq_type/ndata. "
            "Re-run `phyloai posttree dating hessian` to regenerate the hessian directory."
        )
    seq_type = detect_seqtype_from_phylip(dummy_phy)
    n_part = count_ndata_from_phylip(dummy_phy)
    return seq_type, n_part, "dummy.phy-fallback"


def _read_seed_used(run_dir: Path) -> int | None:
    """Return seed from SeedUsed file, or None if not yet written."""
    seed_file = run_dir / "SeedUsed"
    if not seed_file.exists():
        return None
    try:
        return int(seed_file.read_text().strip())
    except Exception:
        return None


# ── mcmc.txt incremental reader (I3) ──────────────────────────────────

class _SampleCounter:
    """Count completed samples in an mcmc.txt that mcmctree writes to
    incrementally. Tracks (inode, byte_offset) so each `count(path)` only
    reads new bytes instead of re-scanning the whole file.

    MCMCtree's mcmc.txt may be partially written (last line truncated)
    when polled mid-run. We tolerate that by skipping any trailing
    unparseable line.
    """

    def __init__(self) -> None:
        self._inode: int | None = None
        self._offset: int = 0
        self._count: int = 0
        self._header_done: bool = False

    def count(self, path: Path) -> int:
        if not path.exists():
            return self._count
        try:
            st = path.stat()
        except OSError:
            return self._count
        if self._inode is not None and st.st_ino != self._inode:
            # File was rotated/replaced — reset and recount from start.
            self._inode = st.st_ino
            self._offset = 0
            self._count = 0
            self._header_done = False
        with open(path, "rb") as fh:
            if not self._header_done:
                fh.seek(0)
                fh.readline()  # discard header
                self._offset = fh.tell()
                self._header_done = True
            fh.seek(self._offset)
            new_bytes = fh.read()
            if not new_bytes:
                return self._count
            self._offset = fh.tell()
            for line in new_bytes.splitlines():
                stripped = line.strip()
                if not stripped:
                    continue
                first = stripped.split()[0]
                try:
                    float(first)
                except ValueError:
                    continue
                self._count += 1
        self._inode = st.st_ino
        return self._count


# ── log tailer (I2) ──────────────────────────────────────────────────

_log_tail_threads: list[Any] = []


def _start_log_tail(log_path: Path, *, prefix: str) -> None:
    """Spawn a background thread that prints new bytes appended to
    `log_path` to the console (prefixed with `prefix`). The thread
    terminates when the file is no longer growing or when the parent
    process exits; we don't bother with a stop signal because the
    daemon thread exits at interpreter shutdown.
    """
    import threading

    def _tail() -> None:
        pos = 0
        try:
            while True:
                if not log_path.exists():
                    time.sleep(1)
                    continue
                with open(log_path, "rb") as fh:
                    fh.seek(0, 2)
                    end = fh.tell()
                if end > pos:
                    with open(log_path, "rb") as fh:
                        fh.seek(pos)
                        chunk = fh.read(end - pos)
                    pos = end
                    try:
                        text = chunk.decode("utf-8", errors="replace")
                    except Exception:
                        text = ""
                    for line in text.splitlines():
                        if line.strip():
                            console.print(f"[dim]{prefix}[/dim] {line}")
                time.sleep(2)
        except Exception:
            pass

    t = threading.Thread(target=_tail, daemon=True, name=f"tail:{prefix}")
    t.start()
    _log_tail_threads.append(t)


def _derive_prior_ctl(posterior_ctl: str, *, seed: int) -> str:
    """Build a prior ctl from a posterior ctl by forcing usedata=0 and seed=<seed>.

    Only the `usedata = ...` and `seed = ...` lines are changed; everything
    else (model, alpha, ncatG, BDparas, rgene_gamma, ...) is preserved
    verbatim from the posterior ctl. If those lines are missing entirely, they
    are appended at the end so MCMCtree still has the values it needs.

    Used when `--ctl` was provided: PhyloAI manages runtime-determined
    parameters (usedata + seed) while the user controls everything else via
    their ctl file.
    """
    import re
    text = posterior_ctl
    usedata_repl = "      usedata = 0    * 0: no data; 1:seq like; 2:use in.BV; 3: out.BV"
    seed_repl = f"          seed = {seed}"

    if re.search(r"^\s*usedata\s*=", text, re.MULTILINE):
        text = re.sub(r"^\s*usedata\s*=.*$", usedata_repl, text, count=1, flags=re.MULTILINE)
    else:
        text = text.rstrip() + "\n\n" + usedata_repl + "\n"

    if re.search(r"^\s*seed\s*=", text, re.MULTILINE):
        text = re.sub(r"^\s*seed\s*=.*$", seed_repl, text, count=1, flags=re.MULTILINE)
    else:
        text = text.rstrip() + "\n" + seed_repl + "\n"

    return text


def _setup_run_dir(
    run_dir: Path,
    hessian_dir: Path,
    ctl_text: str,
) -> None:
    """Create run directory with symlinks and ctl file.

    M5 + I6: hessian_dir is `.resolve()`-d so symlink targets are absolute
    (avoids CWD-relative resolution when MCMCtree runs from inside
    runN/). Symlinks use `os.symlink(str(target), str(link))` directly
    rather than `Path.symlink_to(target.resolve())` because the latter can
    produce invalid paths on macOS HFS when the target is on a different
    filesystem (e.g. /tmp vs /Users).
    """
    import os as _os

    run_dir.mkdir(parents=True, exist_ok=True)
    hessian_abs = hessian_dir.resolve()
    for fname in HESSIAN_OUTPUT_FILES:
        link = run_dir / fname
        target = hessian_abs / fname
        if link.exists() or link.is_symlink():
            link.unlink()
        _os.symlink(str(target), str(link))
    inbv = run_dir / "in.BV"
    if inbv.exists() or inbv.is_symlink():
        inbv.unlink()
    _os.symlink(str(hessian_abs / "iqtree.mcmctree.hessian"), str(inbv))
    (run_dir / "mcmctree.ctl").write_text(ctl_text)


def _detect_mcmctree_version(mcmctree_exe: Path) -> str:
    try:
        result = subprocess.run(
            [str(mcmctree_exe)],
            capture_output=True, text=True, timeout=5,
        )
        import re
        m = re.search(r"paml version (\d+(?:\.\d+)+)", result.stdout + result.stderr)
        if m:
            return m.group(1)
    except Exception:
        pass
    return "unknown"


def run_mcmc(
    *,
    hessian_dir: Path,
    ctl: Path | None = None,
    clock: int = 2,
    burnin: int = 100000,
    sample_freq: int = 10,
    nsamples: int = 10000,
    n_runs: int = 2,
    output_dir: Path,
    mcmctree_path: str | None = None,
    overwrite: bool = False,
    dry_run: bool = False,
    quiet: bool = False,
    seed_wait_timeout_sec: int = 60,
) -> dict[str, Any]:
    """Library entry point for `phyloai posttree dating mcmc`.

    `ctl` (M2): when provided, the file is used verbatim — PhyloAI copies it
    into each runN/ and skips template generation. The `--clock`/`--burnin`/
    `--sample-freq`/`--nsamples` flags have no effect when `ctl` is provided;
    the CLI layer rejects that combination before calling.

    `seed_wait_timeout_sec` (I4): how long, in seconds, to wait for each
    posterior run to write its `SeedUsed` file before giving up on launching
    the matching prior run. Default 60s (1 minute).
    """
    import shutil as _shutil
    from phyloai.core.env import ToolEnv

    t0 = time.time()

    # Validate hessian dir
    errors = validate_hessian_dir(hessian_dir)
    if errors:
        return _error_result(errors[0], "input")

    # Validate --ctl
    if ctl is not None:
        if not ctl.exists():
            return _error_result(f"--ctl does not exist: {ctl}", "input")
        if not ctl.is_file():
            return _error_result(f"--ctl is not a regular file: {ctl}", "input")

    # Resolve mcmctree
    if mcmctree_path:
        mcmctree_exe = Path(mcmctree_path)
    else:
        env = ToolEnv()
        mcmctree_exe = env.get("mcmctree")
        if mcmctree_exe is None:
            return _error_result("mcmctree not found. Install PAML.", "env")

    mcmctree_version = _detect_mcmctree_version(mcmctree_exe)

    # Resolve seq_type + ndata. M1: prefer hessian's result.json; fallback to
    # parsing iqtree.dummy.phy when result.json is absent or missing fields
    # (e.g. when the hessian directory was produced by an external tool).
    seqtype_str, ndata, src = _resolve_seqtype_and_ndata(hessian_dir)

    seqtype_code = SEQTYPE_CODE[seqtype_str]

    # Resolve ctl text: --ctl takes precedence over generated template.
    if ctl is not None:
        ctl_text = ctl.read_text(errors="ignore")
        ctl_source = "user-supplied"
    else:
        ctl_text = generate_mcmctree_ctl(
            seqtype_code=seqtype_code,
            ndata=ndata,
            clock=clock,
            burnin=burnin,
            sampfreq=sample_freq,
            nsample=nsamples,
            usedata=2,
            seed=-1,
        )
        ctl_source = "generated"

    if dry_run:
        return {
            "status": "success",
            "command": f"phyloai posttree dating mcmc --hessian-dir {hessian_dir}",
            "wall_time": 0.0,
            "tool_versions": {"mcmctree": mcmctree_version},
            "params": {
                "ctl": str(ctl) if ctl else None,
                "clock": clock, "burnin": burnin,
                "sample_freq": sample_freq, "nsamples": nsamples,
                "n_runs": n_runs, "seqtype": seqtype_str, "ndata": ndata,
            },
            "key_results": {},
            "error": None,
            "data": {"ctl": ctl_text, "ctl_source": ctl_source, "warnings": []},
        }

    if overwrite and output_dir.exists():
        _shutil.rmtree(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # M5: resolve output_dir and hessian_dir. Output dir because we
    # write `result.json` into it via cwd-relative paths and symlink
    # MCMCtree inputs into per-run subdirs; hessian dir because its files
    # are symlinked into runN/ and runN/prior/ for MCMCtree to read.
    output_dir = output_dir.resolve()
    hessian_dir = hessian_dir.resolve()

    # Write top-level ctl for user inspection
    (output_dir / "mcmctree.ctl").write_text(ctl_text)

    run_dirs = [output_dir / f"run{i+1}" for i in range(n_runs)]

    # Setup all run and prior directories
    for run_dir in run_dirs:
        _setup_run_dir(run_dir, hessian_dir, ctl_text)
        prior_dir = run_dir / "prior"
        # prior ctl will be written once SeedUsed appears; create dir now
        prior_dir.mkdir(parents=True, exist_ok=True)
        # I6: same os.symlink pattern as _setup_run_dir.
        import os as _os
        hessian_abs = hessian_dir.resolve()
        for fname in HESSIAN_OUTPUT_FILES:
            link = prior_dir / fname
            target = hessian_abs / fname
            if link.exists() or link.is_symlink():
                link.unlink()
            _os.symlink(str(target), str(link))
        inbv = prior_dir / "in.BV"
        if inbv.exists() or inbv.is_symlink():
            inbv.unlink()
        _os.symlink(str(hessian_abs / "iqtree.mcmctree.hessian"), str(inbv))

    # I5: force OMP_NUM_THREADS=1 per MCMCtree process so the 4-way parallel
    # run does not oversubscribe CPU. MCMCtree itself is single-threaded;
    # this stops inherited IQ-TREE/HDF5-style thread pools from racing.
    mcmctree_env = {**__import__("os").environ, "OMP_NUM_THREADS": "1"}

    # I1: keep a handle to each Popen's stdout log file so we can close it
    # after the process exits and avoid file-descriptor leaks on long runs.
    proc_log_handles: dict[str, Any] = {}

    def _launch_run(run_dir: Path, *, phase: str) -> subprocess.Popen:
        """Spawn mcmctree in run_dir, returning the Popen. Caller is
        responsible for `wait()` and for closing the stdout log handle
        stored in proc_log_handles[run_key].
        """
        log_path = run_dir / "mcmctree.log"
        # I2: when not --quiet, tee mcmctree stdout to the log file AND
        # to the terminal. We do this by opening the log in write mode and
        # chaining a tee subprocess on top — but a simpler implementation
        # is to just open the log file and let I2 be a tail-thread that
        # follows new bytes. The tail-thread starts here.
        log_fh = open(log_path, "w")
        run_key = f"{run_dir.name}:{phase}"
        proc_log_handles[run_key] = log_fh
        proc = subprocess.Popen(
            [str(mcmctree_exe), "mcmctree.ctl"],
            cwd=run_dir,
            stdout=log_fh,
            stderr=subprocess.STDOUT,
            env=mcmctree_env,
        )
        if not quiet:
            _start_log_tail(log_path, prefix=run_key)
        return proc

    # Launch posterior runs
    procs: dict[str, subprocess.Popen] = {}
    posterior_launch_time: dict[str, float] = {}
    for run_dir in run_dirs:
        procs[run_dir.name] = _launch_run(run_dir, phase="posterior")
        posterior_launch_time[run_dir.name] = time.time()

    prior_procs: dict[str, subprocess.Popen] = {}
    prior_started: set[str] = set()
    prior_timed_out: set[str] = set()
    timeout_warnings: list[str] = []

    # Progress display
    progress = Progress(
        TextColumn(" {task.description:<30}"),
        BarColumn(),
        TaskProgressColumn(),
        TextColumn(" {task.fields[samples]}/{task.fields[total]} samples"),
        TimeElapsedColumn(),
        console=console,
    )
    task_ids: dict[str, Any] = {}
    for run_dir in run_dirs:
        tid = progress.add_task(
            f"{run_dir.name}-posterior",
            samples=0, total=nsamples, completed=0,
        )
        task_ids[f"{run_dir.name}-posterior"] = tid
    for run_dir in run_dirs:
        tid = progress.add_task(
            f"{run_dir.name}-prior",
            samples=0, total=nsamples, completed=0,
        )
        task_ids[f"{run_dir.name}-prior"] = tid

    # I3: incremental mcmc.txt readers — track (inode, byte_offset) per
    # file so each 5s poll only reads new bytes instead of re-scanning
    # the whole file (mcmc.txt grows to ~10000 lines × 4 files = ~50k
    # lines per poll otherwise).
    sample_counters: dict[Path, "_SampleCounter"] = {
        run_dir / "mcmc.txt": _SampleCounter() for run_dir in run_dirs
    }
    sample_counters.update({
        run_dir / "prior" / "mcmc.txt": _SampleCounter() for run_dir in run_dirs
    })

    with Live(progress, console=console, refresh_per_second=0.5):
        while True:
            # Check posterior progress
            for run_dir in run_dirs:
                n = sample_counters[run_dir / "mcmc.txt"].count(run_dir / "mcmc.txt")
                tid = task_ids[f"{run_dir.name}-posterior"]
                progress.update(tid, completed=n, samples=n, total=nsamples)

                # Launch prior once SeedUsed appears. I4: wait at most
                # seed_wait_timeout_sec (default 60s) before giving up on
                # the prior and recording a warning.
                if run_dir.name not in prior_started and run_dir.name not in prior_timed_out:
                    seed = _read_seed_used(run_dir)
                    if seed is not None:
                        prior_started.add(run_dir.name)
                        prior_dir = run_dir / "prior"
                        if ctl is not None:
                            prior_ctl = _derive_prior_ctl(ctl_text, seed=seed)
                        else:
                            prior_ctl = generate_mcmctree_ctl(
                                seqtype_code=seqtype_code,
                                ndata=ndata,
                                clock=clock,
                                burnin=burnin,
                                sampfreq=sample_freq,
                                nsample=nsamples,
                                usedata=0,
                                seed=seed,
                            )
                        (prior_dir / "mcmctree.ctl").write_text(prior_ctl)
                        prior_procs[run_dir.name] = _launch_run(prior_dir, phase="prior")
                    elif time.time() - posterior_launch_time[run_dir.name] > seed_wait_timeout_sec:
                        prior_timed_out.add(run_dir.name)
                        timeout_warnings.append(
                            f"{run_dir.name}: SeedUsed not written within "
                            f"{seed_wait_timeout_sec}s; skipping prior."
                        )

                # Update prior progress
                prior_path = run_dir / "prior" / "mcmc.txt"
                prior_n = sample_counters[prior_path].count(prior_path)
                ptid = task_ids[f"{run_dir.name}-prior"]
                if run_dir.name in prior_started:
                    progress.update(ptid, completed=prior_n, samples=prior_n, total=nsamples)
                elif run_dir.name in prior_timed_out:
                    progress.update(ptid, description=f"{run_dir.name}-prior (timeout)")
                else:
                    progress.update(ptid, description=f"{run_dir.name}-prior (waiting)")

            # Check if all done. Priors that timed out are treated as done (failed).
            post_done = all(p.poll() is not None for p in procs.values())
            prior_done = (
                (len(prior_started) + len(prior_timed_out)) == n_runs
                and all(p.poll() is not None for p in prior_procs.values())
            )
            if post_done and prior_done:
                break

            time.sleep(5)

    # I1: close all log file handles now that every process has exited.
    # Popen does not close the stdout file itself when it exits; the FD
    # leak only matters across many runs but is still good hygiene.
    for fh in proc_log_handles.values():
        try:
            fh.close()
        except Exception:
            pass

    # M3: collect per-process return codes + output validation. Any
    # non-zero returncode OR missing/empty mcmctree.out / mcmc.txt is
    # recorded as a warning; if any posterior failed, overall status is
    # "error" (the result the user actually asked for did not happen).
    run_failures: list[str] = []
    for run_dir in run_dirs:
        run_key = f"{run_dir.name}:posterior"
        rc = procs[run_dir.name].returncode
        if rc != 0:
            run_failures.append(
                f"{run_dir.name}-posterior: mcmctree exited with code {rc}"
            )
        for required in ("mcmctree.out", "mcmc.txt"):
            p = run_dir / required
            if not p.exists():
                run_failures.append(f"{run_key}: missing {required}")
            elif p.stat().st_size == 0:
                run_failures.append(f"{run_key}: empty {required}")
    for run_name in prior_started:
        run_dir = output_dir / run_name
        prior_dir = run_dir / "prior"
        rc = prior_procs[run_name].returncode
        if rc != 0:
            run_failures.append(
                f"{run_name}-prior: mcmctree exited with code {rc}"
            )
        for required in ("mcmctree.out", "mcmc.txt"):
            p = prior_dir / required
            if not p.exists():
                run_failures.append(f"{run_name}-prior: missing {required}")
            elif p.stat().st_size == 0:
                run_failures.append(f"{run_name}-prior: empty {required}")
    all_warnings = timeout_warnings + run_failures

    # D2: if --runs > 2 the convergence diagnostic is computed pairwise
    # by generate_all_diagnostics; if --runs == 1 (I7) convergence is
    # skipped entirely. diagnostics functions handle both cases.

    # Generate diagnostics
    from phyloai.posttree.dating_diagnostics import generate_all_diagnostics
    diag_dir = output_dir / "diagnostics"
    diag_summary = generate_all_diagnostics(
        run_dirs=run_dirs,
        diag_dir=diag_dir,
        n_runs=n_runs,
    )

    # D2/I7: if --runs == 1, mark convergence metrics absent.
    if n_runs < 2:
        diag_summary["convergence"] = {"status": "skipped", "reason": "n_runs=1"}

    wall = time.time() - t0
    # Status: error if any posterior failed (no usable dating result);
    # warning otherwise. Priors failing is a warning only — convergence
    # diagnostics degrade gracefully.
    posterior_failed = any(
        procs[rd.name].returncode != 0 for rd in run_dirs
    ) or any(
        f"-posterior:" in w for w in run_failures
    )
    return {
        "status": "error" if posterior_failed else "success",
        "command": f"phyloai posttree dating mcmc --hessian-dir {hessian_dir}",
        "wall_time": wall,
        "tool_versions": {"mcmctree": mcmctree_version},
        "params": {
            "ctl": str(ctl) if ctl else None,
            "clock": clock, "burnin": burnin,
            "sample_freq": sample_freq, "nsamples": nsamples,
            "n_runs": n_runs, "seqtype": seqtype_str, "ndata": ndata,
            "seqtype_ndata_source": src,
            "ctl_source": ctl_source,
        },
        "key_results": {
            "n_runs": n_runs,
            "n_priors_skipped_timeout": len(prior_timed_out),
            "n_posterior_failures": sum(
                1 for w in run_failures if "-posterior:" in w
            ),
            "convergence_rho_posterior": next(
                (r["rho"] for r in diag_summary.get("spearman", [])
                 if r["comparison"].startswith("convergence_posterior_")), None
            ) if n_runs >= 2 else None,
        },
        "error": run_failures[0] if posterior_failed else None,
        "error_category": "tool" if posterior_failed else None,
        "data": {
            "diagnostics": diag_summary,
            "warnings": all_warnings,
            "return_codes": {
                **{f"{rd.name}:posterior": procs[rd.name].returncode
                   for rd in run_dirs},
                **{f"{rn}:prior": prior_procs[rn].returncode
                   for rn in prior_started},
            },
        },
    }


def _error_result(msg: str, category: str) -> dict[str, Any]:
    return {
        "status": "error",
        "command": "",
        "wall_time": 0.0,
        "tool_versions": {},
        "params": {},
        "key_results": {},
        "error": msg,
        "error_category": category,
        "data": {"warnings": [msg]},
    }
```

- [ ] **Step 4: Run tests**

```bash
pytest tests/posttree/test_dating_mcmc.py -v
```
Expected: All PASS.

- [ ] **Step 5: Report done — no commit**

---

## Task 5: CLI wiring — `dating` subgroup in `posttree.py`

**Files:**
- Modify: `phyloai/cli/commands/posttree.py`

Add `dating` subgroup with `hessian` and `mcmc` commands following the
existing `topology` command pattern exactly.

- [ ] **Step 1: Write failing smoke test**

```python
# tests/cli/test_dating_cli.py
"""Smoke tests for `phyloai posttree dating` CLI."""
from __future__ import annotations
from click.testing import CliRunner
from phyloai.cli.main import cli


def test_dating_group_help():
    runner = CliRunner()
    result = runner.invoke(cli, ["posttree", "dating", "--help"])
    assert result.exit_code == 0
    assert "hessian" in result.output
    assert "mcmc" in result.output


def test_hessian_help():
    runner = CliRunner()
    result = runner.invoke(cli, ["posttree", "dating", "hessian", "--help"])
    assert result.exit_code == 0
    assert "--matrix" in result.output
    assert "--rooted-tree" in result.output
    assert "--seq-type" in result.output
    assert "--partitions" in result.output


def test_mcmc_help():
    runner = CliRunner()
    result = runner.invoke(cli, ["posttree", "dating", "mcmc", "--help"])
    assert result.exit_code == 0
    assert "--hessian-dir" in result.output
    assert "--ctl" in result.output
    assert "--clock" in result.output
    assert "--burnin" in result.output
    assert "--nsamples" in result.output


def test_mcmc_ctl_conflict_rejected(tmp_path):
    """M2: --ctl with --burnin (changed from default) must exit with code 1."""
    hdir = tmp_path / "hessian"
    hdir.mkdir()
    ctl = tmp_path / "mine.ctl"
    ctl.write_text("seed = -1\n")
    runner = CliRunner()
    result = runner.invoke(cli, [
        "posttree", "dating", "mcmc",
        "--hessian-dir", str(hdir),
        "--ctl", str(ctl),
        "--burnin", "50000",
    ])
    assert result.exit_code != 0
    assert "--ctl" in result.output
    assert "mutually exclusive" in result.output.lower() or "ignored" in result.output.lower()


def test_mcmc_dry_run_with_ctl_uses_user_file(tmp_path):
    """M2: --ctl + --dry-run prints the user's ctl, not a generated one."""
    hdir = tmp_path / "hessian"
    hdir.mkdir()
    ctl = tmp_path / "mine.ctl"
    ctl.write_text("      seed = -1\n      model = 7\n")
    runner = CliRunner()
    result = runner.invoke(cli, [
        "posttree", "dating", "mcmc",
        "--hessian-dir", str(hdir),
        "--ctl", str(ctl),
        "--dry-run",
        "-o", str(tmp_path / "out"),
    ])
    assert result.exit_code == 0
    assert "model = 7" in result.output
    assert "Generated" not in result.output  # user-supplied, not generated


def test_hessian_missing_matrix(tmp_path):
    runner = CliRunner()
    result = runner.invoke(cli, [
        "posttree", "dating", "hessian",
        "--matrix", str(tmp_path / "nope.fa"),
        "--rooted-tree", str(tmp_path / "t.nwk"),
        "-o", str(tmp_path / "out"),
    ])
    assert result.exit_code != 0


def test_mcmc_missing_hessian_dir(tmp_path):
    runner = CliRunner()
    result = runner.invoke(cli, [
        "posttree", "dating", "mcmc",
        "--hessian-dir", str(tmp_path / "nope"),
        "-o", str(tmp_path / "out"),
    ])
    assert result.exit_code != 0


def test_hessian_dry_run(tmp_path):
    matrix = tmp_path / "m.fa"
    matrix.write_text(">sp1\nMKTVFLGEI\n>sp2\nMLTVFLGEI\n")
    tree = tmp_path / "t.nwk"
    tree.write_text("(sp1,(sp2,sp3))'<4.2';\n")
    runner = CliRunner()
    result = runner.invoke(cli, [
        "posttree", "dating", "hessian",
        "--matrix", str(matrix),
        "--rooted-tree", str(tree),
        "--dry-run",
        "-o", str(tmp_path / "out"),
    ])
    assert result.exit_code == 0
    assert "Would run" in result.output or "iqtree" in result.output.lower()
```

- [ ] **Step 2: Run to confirm failure**

```bash
pytest tests/cli/test_dating_cli.py -v
```
Expected: FAIL — `dating` group not found.

- [ ] **Step 3: Add `dating` subgroup and commands to `posttree.py`**

In `phyloai/cli/commands/posttree.py`:

First, update `_PosttreeGroup.list_commands` to include `"dating"`:
```python
class _PosttreeGroup(click.Group):
    def list_commands(self, ctx: click.Context) -> list[str]:
        return ["topology", "dating"]
```

Then add at the end of the file (before any trailing newline):

```python
# ── dating subgroup ──────────────────────────────────────────────────

class _DatingGroup(click.Group):
    def list_commands(self, ctx: click.Context) -> list[str]:
        return ["hessian", "mcmc"]


@posttree.group("dating", cls=_DatingGroup)
def dating() -> None:
    """Bayesian molecular dating via MCMCtree approximate likelihood.

    Two-step workflow:

    \b
    1. phyloai posttree dating hessian   — run IQ-TREE3 to compute gradients
                                           and Hessian (approximate likelihood)
    2. phyloai posttree dating mcmc      — run MCMCtree Bayesian dating

    Edit the generated mcmctree.ctl between steps if needed.
    """


@dating.command("hessian")
@click.option("--matrix", type=click.Path(path_type=Path), required=True,
              help="Supermatrix alignment (FASTA/PHYLIP/NEXUS). Maps to IQ-TREE -s.")
@click.option("--rooted-tree", "rooted_tree", type=click.Path(path_type=Path), required=True,
              help=(
                  "Rooted tree with fossil/tip age calibrations in MCMCtree format. "
                  "Must include a root age constraint. "
                  "Example: (A,((B,C)'>3.1<3.8'),(D,E)'>2.9<3.6'))'<4.2'; "
                  "(units: 100 Mya). Maps to IQ-TREE -te."
              ))
@click.option("--seq-type", "seq_type",
              type=click.Choice(["AA", "NT", "auto"], case_sensitive=False),
              default="auto", show_default=True,
              help=(
                  "Sequence type. AA uses LG+F+G4 by default; NT uses GTR+G4. "
                  "auto detects from the alignment."
              ))
@click.option("--model-expr", "model_expr", type=str, default=None,
              help=(
                  "Custom IQ-TREE model expression (e.g. C10+F+G4). "
                  "Mutually exclusive with --partitions. "
                  "Ignored when --partitions is provided."
              ))
@click.option("--partitions", type=click.Path(path_type=Path), default=None,
              help=(
                  "Partition file (RAxML-like or NEXUS .best_model.nex from "
                  "phyloai tree ml iqtree). Recommended: <= 10 partitions for "
                  "MCMCtree (too many partitions narrow node age intervals). "
                  "If > 10 partitions, PhyloAI merges them automatically with "
                  "--merge --rclusterf 10. Maps to IQ-TREE -Q."
              ))
@click.option("-o", "--output-dir", type=click.Path(path_type=Path),
              default=Path("runs/posttree/dating/hessian"), show_default=True,
              help="Output directory.")
@click.option("-t", "--threads", type=int, default=4, show_default=True,
              help="Thread count. Maps to IQ-TREE -T.")
@click.option("--iqtree-path", type=str, default=None,
              help="Explicit path to iqtree3 executable.")
@click.option("--tool-args", type=str, default=None,
              help=(
                  "Additional IQ-TREE arguments appended after managed flags. "
                  "Blocked: -s, --dating."
              ))
@click.option("--overwrite", is_flag=True, default=False,
              help="Delete and recreate output directory.")
@click.option("--resume", is_flag=True, default=False,
              help="Resume interrupted IQ-TREE run (IQ-TREE native checkpoint).")
@click.option("--dry-run", is_flag=True, default=False,
              help="Print the IQ-TREE command without executing.")
@click.option("-q", "--quiet", is_flag=True, default=False,
              help="Suppress terminal output except errors.")
def hessian_command(
    matrix: Path,
    rooted_tree: Path,
    seq_type: str,
    model_expr: str | None,
    partitions: Path | None,
    output_dir: Path,
    threads: int,
    iqtree_path: str | None,
    tool_args: str | None,
    overwrite: bool,
    resume: bool,
    dry_run: bool,
    quiet: bool,
) -> None:
    """Compute gradients and Hessian for MCMCtree approximate likelihood dating.

    Runs IQ-TREE3 with --dating mcmctree to generate three files required by
    MCMCtree:

    \b
      iqtree.dummy.phy        dummy alignment (seqfile in mcmctree.ctl)
      iqtree.rooted.nwk       rooted calibrated tree (treefile in mcmctree.ctl)
      iqtree.mcmctree.hessian gradient/Hessian matrix (renamed to in.BV by mcmc step)

    The rooted tree (--rooted-tree) must be in MCMCtree calibration format with
    fossil/tip age constraints on nodes and a constrained root age, e.g.:

    \b
      (A,((B,(C,D)'>3.1<3.8'),(E,F)'>2.9<3.6'))'<4.2';

    Calibration units are 100 Mya. The root age constraint is mandatory.

    \b
    Model selection:
      --seq-type AA|NT|auto  Sequence type (default: auto — reads FASTA,
                             PHYLIP, and NEXUS through shared format helpers). Default models:
                             LG+F+G4 (AA), GTR+G4 (NT).
      --model-expr           override with any IQ-TREE model string (e.g.
                             C10+F+G4). Mutually exclusive with --partitions.
      --partitions           partition file (RAxML-like or .best_model.nex from
                             phyloai tree ml iqtree). Recommended: <= 10
                             partitions for MCMCtree (too many partitions narrow
                             node age intervals). If > 10, PhyloAI automatically
                             merges them with --merge --rclusterf 10.

    \b
    Examples:

      # Unpartitioned AA analysis (default model LG+F+G4)
      phyloai posttree dating hessian \\
          --matrix concat.aa.fa --rooted-tree calib.tre

      # Custom mixture model
      phyloai posttree dating hessian \\
          --matrix concat.aa.fa --rooted-tree calib.tre --model-expr C10+F+G4

      # Partitioned NT analysis (<= 10 partitions, fixed GTR+G4 per partition)
      phyloai posttree dating hessian \\
          --matrix concat.nt.fa --rooted-tree calib.tre \\
          --partitions loci.partitions

      # Resume interrupted IQ-TREE run
      phyloai posttree dating hessian \\
          --matrix concat.aa.fa --rooted-tree calib.tre --resume
    """
    from phyloai.posttree.dating_hessian import run_hessian

    if not dry_run:
        output_dir = output_dir.resolve()
        if not overwrite and not resume:
            if output_dir.exists() and any(output_dir.iterdir()):
                _fail(
                    f"Output directory exists and is not empty: {output_dir}\n"
                    "Use --overwrite to replace or --resume to reuse.",
                    exit_code=1,
                )
        if overwrite and output_dir.exists():
            shutil.rmtree(output_dir)

    payload = run_hessian(
        matrix=matrix,
        rooted_tree=rooted_tree,
        seq_type=seq_type,
        model_expr=model_expr,
        partitions=partitions,
        output_dir=output_dir,
        threads=threads,
        iqtree_path=iqtree_path,
        tool_args=tool_args,
        overwrite=overwrite,
        resume=resume,
        dry_run=dry_run,
        quiet=quiet,
        stream_output=not quiet,
    )

    if dry_run:
        click.echo(f"Would run: {' '.join(payload['data']['cmd'])}")
        return

    result_path = output_dir / "result.json"
    result_path.parent.mkdir(parents=True, exist_ok=True)
    with open(result_path, "w") as fh:
        json.dump(payload, fh, indent=2)

    if payload["status"] == "error":
        _fail(payload.get("error") or "Unknown error", exit_code=2)

    if not quiet:
        click.echo(f"\nStatus:    {payload['status']}")
        click.echo(f"Wall time: {payload['wall_time']:.1f}s")
        click.echo(f"Hessian:   {payload['data']['output_files'].get('iqtree.mcmctree.hessian')}")
        click.echo(f"Result:    {result_path}")
        click.echo("\nNext step:")
        click.echo(f"  phyloai posttree dating mcmc --hessian-dir {output_dir}")


@dating.command("mcmc")
@click.option("--hessian-dir", "hessian_dir", type=click.Path(path_type=Path), required=True,
              help="Output directory from 'phyloai posttree dating hessian'.")
@click.option("--ctl", "ctl", type=click.Path(path_type=Path), default=None,
              help=(
                  "Use this mcmctree.ctl as-is instead of generating one. "
                  "PhyloAI copies it into each runN/. The matching prior ctl "
                  "is derived by forcing usedata=0 and seed=<SeedUsed> on the "
                  "user's ctl; all other parameters are preserved verbatim. "
                  "Mutually exclusive with --clock/--burnin/--sample-freq/"
                  "--nsamples, which only affect the generated template."
              ))
@click.option("--clock", type=click.Choice(["1", "2", "3"]), default="2", show_default=True,
              help=(
                  "Clock model: 1=global clock (all lineages same rate), "
                  "2=independent rates (recommended for most datasets), "
                  "3=correlated rates (autocorrelated across branches). "
                  "Ignored when --ctl is provided."
              ))
@click.option("--burnin", type=int, default=100000, show_default=True,
              help=(
                  "MCMC burnin iterations (discarded before sampling begins). "
                  "Ignored when --ctl is provided."
              ))
@click.option("--sample-freq", "sample_freq", type=int, default=10, show_default=True,
              help=(
                  "Record one sample every N MCMC iterations. "
                  "Ignored when --ctl is provided."
              ))
@click.option("--nsamples", type=int, default=10000, show_default=True,
              help=(
                  "Number of samples to keep. "
                  "Total iterations = --burnin + (--sample-freq x --nsamples). "
                  "Default: 100000 + (10 x 10000) = 200000 total iterations. "
                  "Ignored when --ctl is provided."
              ))
@click.option("--runs", "n_runs", type=int, default=2, show_default=True,
              help=(
                  "Number of independent posterior MCMC runs. "
                  "Each run is paired with a matching prior run. "
                  "n_runs=1 skips convergence diagnostics; n_runs>=3 "
                  "computes pairwise convergence for every run pair. "
                  "Recommended: 2 (default) for convergence assessment."
              ))
@click.option("-o", "--output-dir", type=click.Path(path_type=Path),
              default=Path("runs/posttree/dating/mcmc"), show_default=True,
              help="Output directory.")
@click.option("--mcmctree-path", type=str, default=None,
              help="Explicit path to mcmctree executable.")
@click.option("--overwrite", is_flag=True, default=False,
              help="Delete and recreate output directory.")
@click.option("--dry-run", is_flag=True, default=False,
              help=(
                  "Generate mcmctree.ctl and print commands without executing. "
                  "Useful for inspecting parameters before a long run."
              ))
@click.option("-q", "--quiet", is_flag=True, default=False,
              help="Suppress terminal output except errors.")
def mcmc_command(
    hessian_dir: Path,
    ctl: Path | None,
    clock: str,
    burnin: int,
    sample_freq: int,
    nsamples: int,
    n_runs: int,
    output_dir: Path,
    mcmctree_path: str | None,
    overwrite: bool,
    dry_run: bool,
    quiet: bool,
) -> None:
    """Run MCMCtree Bayesian molecular dating (approximate likelihood).

    Reads the three IQ-TREE files from a completed 'hessian' run and executes
    MCMCtree to estimate divergence times under a Bayesian framework.

    Two independent posterior runs (run1/, run2/) are launched in parallel
    by default (controlled by --runs), each paired with a matching prior run
    (run1/prior/, run2/prior/) started as soon as the posterior seed is
    available from SeedUsed. Each process uses one CPU thread.

    Use `--dry-run` to preview the generated ctl before a real run. To edit
    the template, save the dry-run output, modify it, and re-run with
    `--ctl edited.ctl`. The `--ctl` flag replaces the auto-generated
    template entirely; the matching prior ctl is derived by forcing
    usedata=0 and seed=<SeedUsed> on your ctl.

    \b
    MCMC settings:
      Total iterations = --burnin + (--sample-freq x --nsamples)
      Default: 100000 + (10 x 10000) = 200000 iterations, 10000 samples kept.
      Increase --nsamples (e.g. 20000) or --burnin for demanding datasets.

    \b
    Clock models (--clock):
      1  Global clock (all lineages same rate)
      2  Independent rates (default; recommended for most datasets)
      3  Correlated rates (autocorrelated across branches)

    Progress is tracked by polling mcmc.txt sample counts for every
    launched run (posterior + prior pairs).

    \b
    Diagnostics generated after all runs complete:
      diagnostics/convergence/         run1 vs run2 scatter + regression line
      diagnostics/infinite_sites/      mean age vs 95% CI width (data sufficiency)
      diagnostics/posterior_vs_prior/  posterior vs prior age per node
      diagnostics/traces/              MCMC parameter trace plots
      diagnostics/spearman_correlations.csv

    \b
    Examples:

      # Default 2-run analysis
      phyloai posttree dating mcmc \\
          --hessian-dir runs/posttree/dating/hessian

      # Longer run with correlated clock
      phyloai posttree dating mcmc \\
          --hessian-dir runs/posttree/dating/hessian \\
          --clock 3 --burnin 200000 --nsamples 20000

      # Use a custom mcmctree.ctl instead of the generated template
      phyloai posttree dating mcmc \\
          --hessian-dir runs/posttree/dating/hessian \\
          --ctl my_run.ctl

      # Dry-run: inspect generated mcmctree.ctl without executing
      phyloai posttree dating mcmc \\
          --hessian-dir runs/posttree/dating/hessian --dry-run
    """
    from phyloai.posttree.dating_mcmc import run_mcmc

    hessian_dir = hessian_dir.resolve()
    if not hessian_dir.exists():
        _fail(f"--hessian-dir does not exist: {hessian_dir}", exit_code=1)

    # --ctl conflicts with non-default template-only flags.
    # Click cannot distinguish an explicit `--clock 2` from the implicit
    # default of 2, so we only reject when the user passes a non-default
    # value (e.g. --clock 1). This matches the updated spec wording:
    # "non-default values of template-only flags conflict with --ctl".
    if ctl is not None:
        conflicting = []
        if clock != "2":
            conflicting.append("--clock")
        if burnin != 100000:
            conflicting.append("--burnin")
        if sample_freq != 10:
            conflicting.append("--sample-freq")
        if nsamples != 10000:
            conflicting.append("--nsamples")
        if conflicting:
            _fail(
                f"--ctl is mutually exclusive with template-only flags: "
                f"{', '.join(conflicting)}. Those flags only affect the "
                f"generated mcmctree.ctl template and are ignored when "
                f"--ctl is provided.",
                exit_code=1,
            )

    if not dry_run:
        output_dir = output_dir.resolve()
        if not overwrite:
            if output_dir.exists() and any(output_dir.iterdir()):
                _fail(
                    f"Output directory exists and is not empty: {output_dir}\n"
                    "Use --overwrite to replace.",
                    exit_code=1,
                )
        if overwrite and output_dir.exists():
            shutil.rmtree(output_dir)

    if n_runs < 1:
        _fail("--runs must be >= 1", exit_code=1)
    if burnin < 0:
        _fail("--burnin must be >= 0", exit_code=1)
    if sample_freq < 1:
        _fail("--sample-freq must be >= 1", exit_code=1)
    if nsamples < 1:
        _fail("--nsamples must be >= 1", exit_code=1)

    payload = run_mcmc(
        hessian_dir=hessian_dir,
        ctl=ctl,
        clock=int(clock),
        burnin=burnin,
        sample_freq=sample_freq,
        nsamples=nsamples,
        n_runs=n_runs,
        output_dir=output_dir,
        mcmctree_path=mcmctree_path,
        overwrite=overwrite,
        dry_run=dry_run,
        quiet=quiet,
    )

    if dry_run:
        click.echo(f"\n--- {payload['data'].get('ctl_source', 'ctl').title()} mcmctree.ctl ---")
        click.echo(payload["data"]["ctl"])
        return

    result_path = output_dir / "result.json"
    result_path.parent.mkdir(parents=True, exist_ok=True)
    with open(result_path, "w") as fh:
        json.dump(payload, fh, indent=2)

    if payload["status"] == "error":
        _fail(payload.get("error") or "Unknown error", exit_code=2)

    if not quiet:
        kr = payload["key_results"]
        click.echo(f"\nStatus:    {payload['status']}")
        click.echo(f"Wall time: {payload['wall_time']:.1f}s")
        click.echo(f"Runs:      {kr.get('n_runs')}")
        rho = kr.get("convergence_rho_posterior")
        if rho is not None:
            click.echo(f"Convergence ρ (posterior): {rho:.4f}")
        click.echo(f"Result:    {result_path}")
        click.echo(f"Diagnostics: {output_dir / 'diagnostics'}")
```

- [ ] **Step 4: Run tests**

```bash
pytest tests/cli/test_dating_cli.py -v
```
Expected: All PASS.

- [ ] **Step 5: Report done — no commit**

---

## Task 6: Full test suite and final verification

**Files:** No new files — run existing tests.

- [ ] **Step 1: Run full test suite**

```bash
cd /Users/zf/data/coding/phyloAI
pytest tests/ -v --tb=short 2>&1 | tail -40
```
Expected: All PASS. Fix any regressions before proceeding.

- [ ] **Step 2: Verify CLI help renders correctly**

```bash
phyloai posttree dating --help
phyloai posttree dating hessian --help
phyloai posttree dating mcmc --help
```
Check: all flags appear with correct defaults and help text.

- [ ] **Step 3: Verify doctor detects mcmctree version (if installed)**

```bash
phyloai doctor
```
Check: mcmctree row shows version number (e.g. `4.10.10`) instead of blank.

- [ ] **Step 4: Smoke test hessian dry-run**

```bash
phyloai posttree dating hessian \
    --matrix runs/posttree/test/matrix.aa.fa \
    --rooted-tree runs/posttree/test/input.tre \
    --dry-run
```
Expected: prints IQ-TREE command, exits 0, no files created.

- [ ] **Step 5: Smoke test mcmc dry-run**

```bash
phyloai posttree dating mcmc \
    --hessian-dir runs/posttree/test/test \
    --dry-run
```
Expected: prints generated mcmctree.ctl, exits 0.

- [ ] **Step 6: Report done — no commit. Wait for user review before staging/committing anything.**

---

## Implementation Deviations (2026-06-26)

The following intentional deviations from the plan were applied during
implementation:

### D1: `--prefix` removed
**Plan:** `--prefix` CLI option exposed, default `iqtree`.
**Implemented:** Removed from CLI. IQ-TREE prefix hardcoded to `iqtree` via
`HESSIAN_PREFIX` constant. Custom `--prefix` would rename IQ-TREE output
files, breaking mcmc's ability to find `iqtree.dummy.phy`,
`iqtree.rooted.nwk`, and `iqtree.mcmctree.hessian`. The hessian step has
exactly one valid output file naming scheme; letting users change it is
a silent footgun.

### D2: `--seq-type auto` FASTA-only
**Plan:** Auto-detection via `core.formats.FormatConverter` supporting
FASTA/PHYLIP/NEXUS.
**Implemented:** Auto-detection only on FASTA files
(`.fa`/`.fas`/`.fasta`/`.faa`/`.fna`/`.aln`). PHYLIP/NEXUS extensions
require explicit `--seq-type AA|NT` to avoid false AA detection from
NEXUS keywords or PHYLIP header lines. Validation error produced if
auto is used with non-FASTA extension.

### D3: `ndata` always counted from `iqtree.dummy.phy`
**Plan:** Store `n_partitions` from original partition file; mcmc reads
it from hessian's result.json.
**Implemented:** `hessian` writes the original partition count in
`params.n_partitions` (1 for unpartitioned). `mcmc` reads `seq_type` from
result.json but **always** counts `ndata` from the actual `iqtree.dummy.phy`
blocks. This is the ground truth: when `--merge --rclusterf 10` reduces
≥ 10 partitions to fewer megapartitions, the count from result.json would
be stale. Counting from the file itself is always correct and requires no
cross-step coordination.

### D4: `command` field = PhyloAI CLI command
**Plan:** `command` field contained the IQ-TREE/mcmctree subprocess command.
**Implemented:** CLI layer overwrites `command` with the full PhyloAI CLI
invocation string (`phyloai posttree dating hessian --matrix ...`) before
writing `result.json`, matching the JSON Output Standard.

### D5: mcmc output directory conflict check
**Plan:** No explicit output dir lifecycle in `run_mcmc` (only `--overwrite`
delete-then-create).
**Implemented:** `run_mcmc` and the CLI handler both check for existing
non-empty output directory; exit 1 with clear message unless `--overwrite`
is set. Matches hessian and the main design's output dir convention
(Section 9.5).

### D6: IQ-TREE path exception handling
**Plan:** `_resolve_iqtree_path` exceptions propagated uncaught.
**Implemented:** `run_hessian` wraps `_resolve_iqtree_path` in try/except
for `ValueError`/`FileNotFoundError`, returning an error payload with
`error_category: "env"`. CLI layer adds a second `except Exception` guard
to catch any remaining unhandled env errors and write proper `result.json`.

### D7: Help text cleanup
- `--partitions` threshold: `<= 10` / `> 10` corrected to `< 10` / `>= 10`
  to match code.
- `\\` line continuations removed from examples (Click ignores them in
  docstrings).
- Added blank lines between `# comment` and command in examples for
  readable help output.
- Dating group help references `--dry-run`/`--ctl` workflow instead of
  "Edit the generated mcmctree.ctl between steps".
- Added mention of `phyloai pretree filter cluster` as a partition source.

### D8: mcmctree version via direct subprocess
**Plan:** Use `ToolEnv.check_all()` (same pattern as `phyloai doctor`).
**Implemented:** `shutil.which("mcmctree")` + `subprocess.run` + regex on
`paml version (\d+(?:\.\d+)+)` — a self-contained ~15-line detection with
no `ToolEnv` dependency.  Iterative debugging showed that `ToolEnv`'s
cached state was unreliable across CLI invocations; the direct approach
eliminates the indirection entirely.

### D9: `extract_node_tree` returns first tree
**Plan:** Scan `Species tree for FigTree.` section for the first tree
with bare-integer internal labels (the node-label tree).
**Implemented:** Simply returns the first Newick tree after the section
marker.  The integer-label filter matched the wrong tree or produced
concatenations under current mcmctree.  The user explicitly requested
"keep only the first tree".

### D10: Trace plot improvements
- Y-axis now displays column name (e.g. `t_n1`), X-axis labelled `iteration`.
- Columns named `Gen`/`iter`/`time` (iteration counters) are skipped.
- `generate_all_diagnostics` reports `generated` and `skipped` lists in
  its return dict so the consumer can log what was produced and what was
  omitted (with reasons).

### D11: Per-run random seeds
**Plan:** `seed = -1` in generated ctl, wait for mcmctree to write
`SeedUsed` file, then copy seed to prior ctl. Timeout after 60s.
**Implemented:** `run_mcmc` generates a `random.randint(1, 2**31-1)` seed
for each run, injects it into the ctl via regex, and launches posterior
and prior processes immediately with the same seed. The `SeedUsed`
waiting loop, `seed_wait_timeout_sec`, and timeout warnings are all
removed.
**Rationale:** mcmctree only writes `SeedUsed` when `seed = -1`; with
user-supplied ctls or auto-generated seeds the file is never created.
Generating seeds in Python guarantees posterior/prior identity by
construction.

### D12: Combined convergence CSV
**Plan:** One pairwise CSV per run pair (e.g. `posterior_times_run1_vs_run3.csv`).
**Implemented:** Single `posterior_times.csv` and `prior_times.csv` with
all runs as columns (`mean_run1`, `lower_run1`, `upper_run1`,
`ci_width_run1`, `mean_run2`, ...).  Pairwise convergence scatter plots
are still generated for every run pair.

