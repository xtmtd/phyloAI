# PhyloAI posttree dating Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Implement `phyloai posttree dating hessian` (IQ-TREE3 approximate likelihood) and `phyloai posttree dating mcmc` (MCMCtree Bayesian dating + full diagnostics), plus fix mcmctree version detection in `doctor`.

**Architecture:** Two library modules (`dating_hessian.py`, `dating_mcmc.py`) plus a diagnostics helper (`dating_diagnostics.py`) in `phyloai/posttree/`. CLI wiring added to the existing `posttree.py` as a new `dating` subgroup. MCMC progress tracked via `rich.Live` polling `mcmc.txt` line counts; four processes (2 posterior + 2 prior) run in parallel threads. All plots PDF via matplotlib.

**Tech Stack:** Python 3.10+, Click, rich, matplotlib, scipy.stats, subprocess, threading, pathlib. No new dependencies.

---

## File Map

| File | Action | Responsibility |
|------|--------|---------------|
| `phyloai/core/env.py` | Modify | Fix mcmctree version detection |
| `phyloai/posttree/dating_hessian.py` | Create | IQ-TREE hessian computation library layer |
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

- [ ] **Step 6: Commit**

```bash
git add phyloai/core/env.py tests/core/test_env_mcmctree.py
git commit -m "fix(doctor): detect mcmctree version from paml stdout"
```

---

## Task 2: `dating_hessian.py` — pure helpers

**Files:**
- Create: `phyloai/posttree/dating_hessian.py`
- Create: `tests/posttree/test_dating_hessian.py`

Implement all pure (non-subprocess) helpers needed for the hessian step:
seq-type detection from a PHYLIP file, partition count from a partition file,
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
    detect_seqtype_from_fasta,
    count_partitions,
    validate_root_age,
    build_iqtree_dating_cmd,
    HESSIAN_OUTPUT_FILES,
)


# ── detect_seqtype_from_fasta ────────────────────────────────────────

def test_detect_aa(tmp_path):
    fa = tmp_path / "aa.fa"
    fa.write_text(">sp1\nMKTVFLGEI\n>sp2\nMLTVFLGEI\n")
    assert detect_seqtype_from_fasta(fa) == "AA"

def test_detect_nt(tmp_path):
    fa = tmp_path / "nt.fa"
    fa.write_text(">sp1\nACGTACGT\n>sp2\nACGTACGT\n")
    assert detect_seqtype_from_fasta(fa) == "NT"

def test_detect_auto_defaults_aa_for_mixed(tmp_path):
    """Sequences with non-ACGTN chars → AA."""
    fa = tmp_path / "m.fa"
    fa.write_text(">sp1\nACGTMKLWI\n")
    assert detect_seqtype_from_fasta(fa) == "AA"


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
    assert cmd[idx + 1] == "HKY+G4"

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

from phyloai.core.iqtree import _resolve_iqtree_path, _detect_iqtree_version

HESSIAN_OUTPUT_FILES = (
    "iqtree.dummy.phy",
    "iqtree.rooted.nwk",
    "iqtree.mcmctree.hessian",
)

_AA_ONLY_CHARS = set("ARNDCQEGHILKMFPSTWYV")


def detect_seqtype_from_fasta(matrix: Path) -> str:
    """Return 'AA' or 'NT' by scanning sequence characters."""
    content = matrix.read_text(errors="ignore").upper()
    seqs = []
    for line in content.splitlines():
        line = line.strip()
        if line and not line.startswith(">"):
            seqs.append(line)
    sample = "".join(seqs)[:2000].replace("-", "").replace("?", "").replace("N", "")
    nt_chars = set("ACGTURYSWKMBDHV")
    aa_only = set(sample) - nt_chars
    if aa_only & _AA_ONLY_CHARS:
        return "AA"
    return "NT"


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
    prefix: str,
    threads: int,
    tool_args: str | None,
) -> list[str]:
    """Build the IQ-TREE3 command list for --dating mcmctree."""
    cmd: list[str] = [str(iqtree_path), "-s", str(matrix), "-te", str(rooted_tree),
                      "--dating", "mcmctree", "--prefix", prefix, "-T", str(threads)]

    if partitions is None:
        # Unpartitioned mode
        model = model_expr or ("LG+F+G4" if seq_type == "AA" else "HKY+G4")
        cmd += ["-m", model]
    else:
        # Partitioned mode — model_expr ignored
        cmd += ["-m", "MF", "-Q", str(partitions)]
        if seq_type == "AA":
            cmd += ["--mset", "LG", "-mfreq", "F", "-mrate", "G"]
        else:
            cmd += ["--mset", "HKY", "-mrate", "G"]
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
        blocked = {"-s", "--dating"}
        tokens = shlex.split(tool_args)
        for tok in tokens:
            if tok in blocked:
                errors.append(
                    f"--tool-args contains blocked flag '{tok}' "
                    f"(managed by PhyloAI)."
                )

    return errors


def run_hessian(
    *,
    matrix: Path,
    rooted_tree: Path,
    seq_type: str = "auto",
    model_expr: str | None = None,
    partitions: Path | None = None,
    prefix: str = "iqtree",
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

    # Seq-type detection
    if seq_type == "auto":
        seq_type = detect_seqtype_from_fasta(matrix)

    # Partition count
    n_partitions = 0
    if partitions:
        n_partitions = count_partitions(partitions)

    # Resolve IQ-TREE
    iqtree_exe = _resolve_iqtree_path(iqtree_path)
    iqtree_version = _detect_iqtree_version(iqtree_exe)

    cmd = build_iqtree_dating_cmd(
        iqtree_path=iqtree_exe,
        matrix=matrix,
        rooted_tree=rooted_tree,
        seq_type=seq_type,
        model_expr=model_expr,
        partitions=partitions,
        n_partitions=n_partitions,
        prefix=prefix,
        threads=threads,
        tool_args=tool_args,
    )

    if dry_run:
        return {
            "status": "success",
            "command": " ".join(cmd),
            "wall_time": 0.0,
            "tool_versions": {"iqtree3": iqtree_version},
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

    # Validate outputs
    missing = [f for f in HESSIAN_OUTPUT_FILES if not (output_dir / f).exists()]
    if missing or proc.returncode != 0:
        return {
            "status": "error",
            "command": " ".join(cmd),
            "wall_time": time.time() - t0,
            "tool_versions": {"iqtree3": iqtree_version},
            "params": {"seq_type": seq_type, "n_partitions": n_partitions},
            "key_results": {},
            "error": f"IQ-TREE failed (returncode={proc.returncode}). Missing: {missing}",
            "error_category": "tool",
            "data": {"cmd": cmd, "tool_stderr": getattr(proc, "stderr", ""), "warnings": []},
        }

    wall = time.time() - t0
    return {
        "status": "success",
        "command": " ".join(cmd),
        "wall_time": wall,
        "tool_versions": {"iqtree3": iqtree_version},
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
            "warnings": [],
            "output_files": {f: str(output_dir / f) for f in HESSIAN_OUTPUT_FILES},
        },
    }
```

- [ ] **Step 4: Run tests**

```bash
pytest tests/posttree/test_dating_hessian.py -v
```
Expected: All PASS.

- [ ] **Step 5: Commit**

```bash
git add phyloai/posttree/dating_hessian.py tests/posttree/test_dating_hessian.py
git commit -m "feat(dating): dating_hessian library layer with IQ-TREE command builder"
```

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
(
(sp1:0.32,sp2:0.31)node7:0.12,sp3:0.44)node8;

(
(sp1,sp2)node7,sp3)node8;

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


def test_extract_node_tree_returns_first_tree():
    tree = extract_node_tree(SAMPLE_OUT)
    assert tree is not None
    assert "node7" in tree
    assert "node8" in tree


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
    """Extract first Newick tree after the FigTree node-label section."""
    marker = "Species tree for FigTree."
    idx = text.find(marker)
    if idx == -1:
        return None
    after = text[idx:]
    # Find first complete newick (ends with ;)
    m = re.search(r"(\([\s\S]+?;)", after)
    if m:
        return m.group(1).strip()
    return None


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
    if len(x) < 3:
        return float("nan"), float("nan")
    rho, pval = spearmanr(x, y)
    return float(rho), float(pval)


def plot_convergence(
    table: list[dict],
    x_col: str,
    y_col: str,
    out_path: Path,
    title: str,
    xlabel: str,
    ylabel: str,
) -> tuple[float, float]:
    """Scatter + regression line. Returns (rho, pvalue)."""
    x = [r[x_col] for r in table if not np.isnan(r.get(x_col, float("nan")))]
    y = [r[y_col] for r in table if not np.isnan(r.get(y_col, float("nan")))]
    rho, pval = _spearman(x, y)

    fig, ax = plt.subplots(figsize=(5, 5))
    ax.scatter(x, y, s=40, alpha=0.8, color="steelblue", zorder=3)
    if len(x) >= 2:
        m, b = np.polyfit(x, y, 1)
        xs = np.linspace(min(x), max(x), 100)
        ax.plot(xs, m * xs + b, color="firebrick", lw=1.5, label=f"ρ={rho:.3f}, p={pval:.3g}")
        ax.legend(fontsize=8)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)
    return rho, pval


def plot_line(
    x: list[float],
    y: list[float],
    out_path: Path,
    title: str,
    xlabel: str,
    ylabel: str,
) -> tuple[float, float]:
    """Line plot (points sorted by x and connected). Returns (rho, pvalue)."""
    pairs = sorted(zip(x, y))
    xs = [p[0] for p in pairs]
    ys = [p[1] for p in pairs]
    rho, pval = _spearman(xs, ys)

    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(xs, ys, marker="o", markersize=4, color="steelblue", lw=1.2)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)
    return rho, pval


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
    """Parse outputs and generate all diagnostics. Returns summary dict."""
    corr_rows: list[dict] = []

    def _add_corr(comparison: str, rho: float, pval: float) -> None:
        corr_rows.append({"comparison": comparison, "rho": rho, "pvalue": pval})

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
            # Extract node tree
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

    # Convergence plots (run1 vs run2)
    if len(post_times) >= 2 and post_times[0] and post_times[1]:
        table = build_time_table(post_times[0], post_times[1])
        write_time_table_csv(table, diag_dir / "convergence" / "posterior_times.txt")
        rho, pval = plot_convergence(
            table, "mean_run1", "mean_run2",
            diag_dir / "convergence" / "convergence_posterior.pdf",
            title="Convergence — posterior means",
            xlabel="Mean age run1 (100 Mya)",
            ylabel="Mean age run2 (100 Mya)",
        )
        _add_corr("convergence_posterior", rho, pval)

    if len(prior_times) >= 2 and prior_times[0] and prior_times[1]:
        table = build_time_table(prior_times[0], prior_times[1])
        write_time_table_csv(table, diag_dir / "convergence" / "prior_times.txt")
        rho, pval = plot_convergence(
            table, "mean_run1", "mean_run2",
            diag_dir / "convergence" / "convergence_prior.pdf",
            title="Convergence — prior means",
            xlabel="Mean age run1 (100 Mya)",
            ylabel="Mean age run2 (100 Mya)",
        )
        _add_corr("convergence_prior", rho, pval)

    # Infinite-sites and posterior-vs-prior plots
    for i, run_dir in enumerate(run_dirs):
        run_label = f"run{i+1}"
        post = post_times[i] if i < len(post_times) else []
        prior = prior_times[i] if i < len(prior_times) else []

        for kind, rows in [("posterior", post), ("prior", prior)]:
            if not rows:
                continue
            x = [r["mean"] for r in rows]
            y = [r["ci_width"] for r in rows]
            rho, pval = plot_line(
                x, y,
                diag_dir / "infinite_sites" / f"infinite_sites_{run_label}_{kind}.pdf",
                title=f"Infinite-sites — {run_label} {kind}",
                xlabel="Mean age (100 Mya)",
                ylabel="95% CI width (100 Mya)",
            )
            _add_corr(f"infinite_sites_{run_label}_{kind}", rho, pval)

        if post and prior:
            post_by_node = {r["node"]: r["mean"] for r in post}
            prior_by_node = {r["node"]: r["mean"] for r in prior}
            shared = [n for n in post_by_node if n in prior_by_node]
            if shared:
                xp = [post_by_node[n] for n in shared]
                yp = [prior_by_node[n] for n in shared]
                rho, pval = plot_line(
                    xp, yp,
                    diag_dir / "posterior_vs_prior" / f"posterior_vs_prior_{run_label}.pdf",
                    title=f"Posterior vs prior — {run_label}",
                    xlabel="Posterior mean age (100 Mya)",
                    ylabel="Prior mean age (100 Mya)",
                )
                _add_corr(f"posterior_vs_prior_{run_label}", rho, pval)

    # Write Spearman CSV
    if corr_rows:
        corr_path = diag_dir / "spearman_correlations.csv"
        corr_path.parent.mkdir(parents=True, exist_ok=True)
        with open(corr_path, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=["comparison", "rho", "pvalue"])
            writer.writeheader()
            writer.writerows(corr_rows)

    return {"spearman": corr_rows}
```

- [ ] **Step 4: Run tests**

```bash
pytest tests/posttree/test_dating_diagnostics.py -v
```
Expected: All PASS.

- [ ] **Step 5: Commit**

```bash
git add phyloai/posttree/dating_diagnostics.py tests/posttree/test_dating_diagnostics.py
git commit -m "feat(dating): dating_diagnostics — parsing and plot generation"
```

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
    """Return a mcmctree.ctl string."""
    return f"""\
          seed = {seed}
       seqfile = iqtree.dummy.phy
      treefile = iqtree.rooted.nwk
       outfile = mcmctree.out

         ndata = {ndata}
       seqtype = {seqtype_code}   * 0: nucleotides; 1:codons; 2:AAs
       usedata = {usedata}        * 0: no data; 1:seq like; 2:use in.BV; 3: out.BV
         clock = {clock}          * 1: global clock; 2: independent rates; 3: correlated rates
       RootAge =                  * safe constraint on root age

     cleandata = 0

       BDparas = 1 1 0.1 M
   rgene_gamma = 2 20 1
  sigma2_gamma = 1 10 1

      finetune = 0: .1  .1  .1  .1 .1 .1

*** These parameters control the MCMC run
***  Note: Total number of MCMC iterations will be burnin + (sampfreq * nsample)

         print = 1
        burnin = {burnin}
      sampfreq = {sampfreq}
       nsample = {nsample}


*** The following parameters only needed to run MCMCtree with exact likelihood (usedata = 1)
*** no need to change anything for approximate likelihood (usedata = 2)

         model = 0
         alpha = 0.5
         ncatG = 4

   kappa_gamma = 6 2
   alpha_gamma = 1 1
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


def _read_seed_used(run_dir: Path) -> int | None:
    """Return seed from SeedUsed file, or None if not yet written."""
    seed_file = run_dir / "SeedUsed"
    if not seed_file.exists():
        return None
    try:
        return int(seed_file.read_text().strip())
    except Exception:
        return None


def _setup_run_dir(
    run_dir: Path,
    hessian_dir: Path,
    ctl_text: str,
) -> None:
    """Create run directory with symlinks and ctl file."""
    run_dir.mkdir(parents=True, exist_ok=True)
    # Symlinks to hessian files
    for fname in HESSIAN_OUTPUT_FILES:
        link = run_dir / fname
        target = hessian_dir / fname
        if link.exists() or link.is_symlink():
            link.unlink()
        link.symlink_to(target.resolve())
    # in.BV symlink
    inbv = run_dir / "in.BV"
    if inbv.exists() or inbv.is_symlink():
        inbv.unlink()
    inbv.symlink_to((hessian_dir / "iqtree.mcmctree.hessian").resolve())
    # ctl file
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
) -> dict[str, Any]:
    """Library entry point for `phyloai posttree dating mcmc`."""
    import shutil as _shutil
    from phyloai.core.env import ToolEnv

    t0 = time.time()

    # Validate hessian dir
    errors = validate_hessian_dir(hessian_dir)
    if errors:
        return _error_result(errors[0], "input")

    # Resolve mcmctree
    if mcmctree_path:
        mcmctree_exe = Path(mcmctree_path)
    else:
        env = ToolEnv()
        mcmctree_exe = env.get("mcmctree")
        if mcmctree_exe is None:
            return _error_result("mcmctree not found. Install PAML.", "env")

    mcmctree_version = _detect_mcmctree_version(mcmctree_exe)

    # Infer seqtype and ndata from dummy.phy
    dummy_phy = hessian_dir / "iqtree.dummy.phy"
    seqtype_str = detect_seqtype_from_phylip(dummy_phy)
    seqtype_code = SEQTYPE_CODE[seqtype_str]
    ndata = count_ndata_from_phylip(dummy_phy)

    # Generate top-level ctl (usedata=2 template)
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

    if dry_run:
        return {
            "status": "success",
            "command": f"phyloai posttree dating mcmc --hessian-dir {hessian_dir}",
            "wall_time": 0.0,
            "tool_versions": {"mcmctree": mcmctree_version},
            "params": {
                "clock": clock, "burnin": burnin,
                "sample_freq": sample_freq, "nsamples": nsamples,
                "n_runs": n_runs, "seqtype": seqtype_str, "ndata": ndata,
            },
            "key_results": {},
            "error": None,
            "data": {"ctl": ctl_text, "warnings": []},
        }

    if overwrite and output_dir.exists():
        _shutil.rmtree(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Write top-level ctl for user inspection
    (output_dir / "mcmctree.ctl").write_text(ctl_text)

    run_dirs = [output_dir / f"run{i+1}" for i in range(n_runs)]

    # Setup all run and prior directories
    for run_dir in run_dirs:
        _setup_run_dir(run_dir, hessian_dir, ctl_text)
        prior_dir = run_dir / "prior"
        # prior ctl will be written once SeedUsed appears; create dir now
        prior_dir.mkdir(parents=True, exist_ok=True)
        for fname in HESSIAN_OUTPUT_FILES:
            link = prior_dir / fname
            target = hessian_dir / fname
            if link.exists() or link.is_symlink():
                link.unlink()
            link.symlink_to(target.resolve())
        inbv = prior_dir / "in.BV"
        if inbv.exists() or inbv.is_symlink():
            inbv.unlink()
        inbv.symlink_to((hessian_dir / "iqtree.mcmctree.hessian").resolve())

    # Launch posterior runs
    procs: dict[str, subprocess.Popen] = {}
    for run_dir in run_dirs:
        log = open(run_dir / "mcmctree.log", "w")
        proc = subprocess.Popen(
            [str(mcmctree_exe), "mcmctree.ctl"],
            cwd=run_dir,
            stdout=log,
            stderr=subprocess.STDOUT,
        )
        procs[run_dir.name] = proc

    prior_procs: dict[str, subprocess.Popen] = {}
    prior_started: set[str] = set()

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

    with Live(progress, console=console, refresh_per_second=0.5):
        while True:
            # Check posterior progress
            for run_dir in run_dirs:
                n = count_mcmc_samples(run_dir / "mcmc.txt")
                tid = task_ids[f"{run_dir.name}-posterior"]
                progress.update(tid, completed=n, samples=n, total=nsamples)

                # Launch prior once SeedUsed appears
                if run_dir.name not in prior_started:
                    seed = _read_seed_used(run_dir)
                    if seed is not None:
                        prior_started.add(run_dir.name)
                        prior_dir = run_dir / "prior"
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
                        log = open(prior_dir / "mcmctree.log", "w")
                        prior_proc = subprocess.Popen(
                            [str(mcmctree_exe), "mcmctree.ctl"],
                            cwd=prior_dir,
                            stdout=log,
                            stderr=subprocess.STDOUT,
                        )
                        prior_procs[run_dir.name] = prior_proc

                # Update prior progress
                prior_n = count_mcmc_samples(run_dir / "prior" / "mcmc.txt")
                ptid = task_ids[f"{run_dir.name}-prior"]
                if run_dir.name in prior_started:
                    progress.update(ptid, completed=prior_n, samples=prior_n, total=nsamples)
                else:
                    progress.update(ptid, description=f"{run_dir.name}-prior (waiting)")

            # Check if all done
            post_done = all(p.poll() is not None for p in procs.values())
            prior_done = (
                len(prior_procs) == n_runs
                and all(p.poll() is not None for p in prior_procs.values())
            )
            if post_done and prior_done:
                break

            time.sleep(5)

    # Generate diagnostics
    from phyloai.posttree.dating_diagnostics import generate_all_diagnostics
    diag_dir = output_dir / "diagnostics"
    diag_summary = generate_all_diagnostics(
        run_dirs=run_dirs,
        diag_dir=diag_dir,
        n_runs=n_runs,
    )

    wall = time.time() - t0
    return {
        "status": "success",
        "command": f"phyloai posttree dating mcmc --hessian-dir {hessian_dir}",
        "wall_time": wall,
        "tool_versions": {"mcmctree": mcmctree_version},
        "params": {
            "clock": clock, "burnin": burnin,
            "sample_freq": sample_freq, "nsamples": nsamples,
            "n_runs": n_runs, "seqtype": seqtype_str, "ndata": ndata,
        },
        "key_results": {
            "n_runs": n_runs,
            "convergence_rho_posterior": next(
                (r["rho"] for r in diag_summary.get("spearman", [])
                 if r["comparison"] == "convergence_posterior"), None
            ),
        },
        "error": None,
        "data": {"diagnostics": diag_summary, "warnings": []},
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

- [ ] **Step 5: Commit**

```bash
git add phyloai/posttree/dating_mcmc.py tests/posttree/test_dating_mcmc.py
git commit -m "feat(dating): dating_mcmc library layer with parallel run management"
```

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
    assert "--clock" in result.output
    assert "--burnin" in result.output
    assert "--nsamples" in result.output


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
                  "Sequence type. AA uses LG+F+G4 by default; NT uses HKY+G4. "
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
@click.option("--prefix", type=str, default="iqtree", show_default=True,
              help="IQ-TREE output prefix. Maps to IQ-TREE --prefix.")
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
    prefix: str,
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
      --seq-type AA|NT|auto  detects sequence type from the alignment (default:
                             auto). Default models: LG+F+G4 (AA), HKY+G4 (NT).
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

      # Partitioned NT analysis (<= 10 partitions, fixed HKY+G4 per partition)
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
        prefix=prefix,
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
@click.option("--clock", type=click.Choice(["1", "2", "3"]), default="2", show_default=True,
              help=(
                  "Clock model: 1=global clock (all lineages same rate), "
                  "2=independent rates (recommended for most datasets), "
                  "3=correlated rates (autocorrelated across branches)."
              ))
@click.option("--burnin", type=int, default=100000, show_default=True,
              help="MCMC burnin iterations (discarded before sampling begins).")
@click.option("--sample-freq", "sample_freq", type=int, default=10, show_default=True,
              help="Record one sample every N MCMC iterations.")
@click.option("--nsamples", type=int, default=10000, show_default=True,
              help=(
                  "Number of samples to keep. "
                  "Total iterations = --burnin + (--sample-freq x --nsamples). "
                  "Default: 100000 + (10 x 10000) = 200000 total iterations."
              ))
@click.option("--runs", "n_runs", type=int, default=2, show_default=True,
              help=(
                  "Number of independent posterior MCMC runs. "
                  "Each run is paired with a matching prior run. "
                  "Use >= 2 to assess convergence."
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

    Two independent posterior runs (run1/, run2/) are launched in parallel,
    each paired with a matching prior run (run1/prior/, run2/prior/) started
    as soon as the posterior seed is available from SeedUsed. All four runs
    use one CPU thread each (4 threads total).

    A mcmctree.ctl control file is generated in the output directory before
    any run starts. You may inspect and edit it freely — the file is copied
    into each run directory at launch time.

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

    Progress is tracked by polling mcmc.txt sample counts for all four runs.

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

      # Dry-run: inspect generated mcmctree.ctl without executing
      phyloai posttree dating mcmc \\
          --hessian-dir runs/posttree/dating/hessian --dry-run
    """
    from phyloai.posttree.dating_mcmc import run_mcmc

    hessian_dir = hessian_dir.resolve()
    if not hessian_dir.exists():
        _fail(f"--hessian-dir does not exist: {hessian_dir}", exit_code=1)

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
        click.echo("\n--- Generated mcmctree.ctl ---")
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

- [ ] **Step 5: Commit**

```bash
git add phyloai/cli/commands/posttree.py tests/cli/test_dating_cli.py
git commit -m "feat(dating): CLI wiring for posttree dating hessian + mcmc commands"
```

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

- [ ] **Step 6: Final commit**

```bash
git add -A
git commit -m "feat(dating): complete posttree dating implementation"
```
