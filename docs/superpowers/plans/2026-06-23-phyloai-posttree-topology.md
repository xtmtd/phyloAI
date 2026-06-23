# Posttree Topology Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Implement `phyloai posttree topology` — IQ-TREE tree topology tests (AU/KH/SH/WKH/WSH/c-ELW) comparing a set of candidate trees against a supermatrix alignment.

**Architecture:** Single-mode IQ-TREE invocation. Shared IQ-TREE helpers (`_resolve_iqtree_path`, `_detect_iqtree_version`) extracted to `phyloai/core/iqtree.py` as the canonical location; both `ml_iqtree.py` and `cf.py` are refactored to import from there. Library function `run_topology()` enforces spec-mandated input validation, builds the IQ-TREE command, executes it, parses the USER TREES table via header-offset-based tolerant parsing, and assembles `result.json`. CLI layer handles output-directory lifecycle (conflict, overwrite, resume) and writes `result.json`.

**Tech Stack:** Python 3.12+, Click, Rich, BioPython, pytest, IQ-TREE3

**Spec:** `docs/superpowers/specs/2026-06-22-phyloai-posttree-topology-design.md`

**Commit policy:** All code is written across tasks first. A single final commit is made at Task 9 after full test suite passes. Individual task "commit" steps are marked `[skip commit]` — they are checkpoints, not required commit points.

---

## File Structure

| File | Purpose |
|------|---------|
| `phyloai/core/iqtree.py` (CREATE) | Shared IQ-TREE helpers: path resolution, version detection |
| `phyloai/tree/ml_iqtree.py` (MODIFY) | Replace local `_resolve_iqtree_path` / `_detect_iqtree_version` with import from `core/iqtree` |
| `phyloai/tree/cf.py` (MODIFY) | Replace local `_resolve_iqtree_path` / `_detect_iqtree_version` with import from `core/iqtree` |
| `phyloai/posttree/__init__.py` (CREATE) | Posttree package marker |
| `phyloai/posttree/topology.py` (CREATE) | Core library: validation, command builder, runner, USER TREES parser, result assembly |
| `phyloai/cli/commands/posttree.py` (CREATE) | Click group + topology command with rich 8-section help |
| `phyloai/cli/main.py` (MODIFY) | Register `posttree` group |
| `tests/posttree/__init__.py` (CREATE) | Test package marker |
| `tests/posttree/test_topology.py` (CREATE) | Library-level unit + integration tests |
| `tests/cli/test_posttree_topology.py` (CREATE) | CLI validation, help content, dry-run, integration tests |
| `docs/commands/posttree-topology.md` (CREATE) | User-facing command documentation |

---

### Task 1: Create shared IQ-TREE helpers and refactor existing modules

**Files:**
- Create: `phyloai/core/iqtree.py`
- Modify: `phyloai/tree/ml_iqtree.py` (lines 582-619)
- Modify: `phyloai/tree/cf.py` (lines 386-420)

- [ ] **Step 1: Create `phyloai/core/iqtree.py`**

Combine the best of each existing version:
- `_resolve_iqtree_path`: adopt cf.py's `.resolve()` pattern (always returns an absolute path; safer when subprocess changes cwd)
- `_detect_iqtree_version`: adopt ml_iqtree.py's dynamic `exe_name` key and two-round regex fallback (more robust)

```python
"""Shared IQ-TREE helpers used by tree and posttree modules."""
from __future__ import annotations

import os
import re as _re
import subprocess
from pathlib import Path

from phyloai.core.env import ToolEnv

IQTREE_COMPATIBLE_EXTENSIONS = frozenset({
    ".fa", ".fas", ".fasta", ".faa", ".fna",
    ".phy", ".phylip",
    ".nex", ".nxs", ".nexus",
    ".aln",
})


def _resolve_iqtree_path(iqtree_path: str | None, dry_run: bool) -> str:
    """Resolve IQ-TREE executable: custom path > PATH > bundled.

    Returns an absolute path string.  Raises ValueError for missing
    or non-executable custom paths; raises FileNotFoundError when no
    IQ-TREE can be found on the system.
    """
    if iqtree_path:
        p = Path(iqtree_path).resolve()
        if not p.exists():
            raise ValueError(f"--iqtree-path does not exist: {iqtree_path}")
        if not os.access(str(p), os.X_OK):
            raise ValueError(f"--iqtree-path is not executable: {iqtree_path}")
        return str(p)
    if dry_run:
        return "iqtree3"
    try:
        env = ToolEnv()
        return str(env.require("iqtree3"))
    except FileNotFoundError:
        raise FileNotFoundError(
            "iqtree3 not found. Install from https://github.com/iqtree/iqtree3/releases "
            "or use --iqtree-path."
        )


def _detect_iqtree_version(executable: str) -> dict[str, str]:
    """Detect IQ-TREE version from --version output."""
    exe_name = Path(executable).name
    try:
        proc = subprocess.run(
            [executable, "--version"],
            capture_output=True, text=True, timeout=10,
        )
        combined = proc.stdout + proc.stderr
    except Exception:
        return {exe_name: "unknown"}

    m = _re.search(r"version\s*([\d.]+)", combined, _re.IGNORECASE)
    if m:
        return {exe_name: m.group(1)}
    m = _re.search(r"([\d]+\.[\d]+(?:\.[\d]+)?)", combined)
    if m:
        return {exe_name: m.group(1)}
    return {exe_name: "unknown"}
```

- [ ] **Step 2: Run a syntax check on the new module**

```bash
python -c "from phyloai.core.iqtree import _resolve_iqtree_path, _detect_iqtree_version, IQTREE_COMPATIBLE_EXTENSIONS; print('OK')"
```

Expected: `OK`

- [ ] **Step 3: Refactor `phyloai/tree/ml_iqtree.py`**

Remove the local definitions of `_resolve_iqtree_path` (lines 582-599) and `_detect_iqtree_version` (lines 602-619). Add an import at the top of the file (near the existing `from phyloai.core.env import ToolEnv`):

```python
from phyloai.core.iqtree import (
    _resolve_iqtree_path,
    _detect_iqtree_version,
)
```

- [ ] **Step 4: Refactor `phyloai/tree/cf.py`**

Remove the local definitions of `_resolve_iqtree_path` (lines 386-404) and `_detect_iqtree_version` (lines 407-420). Add the same import near existing `ToolEnv` import:

```python
from phyloai.core.iqtree import (
    _resolve_iqtree_path,
    _detect_iqtree_version,
)
```

- [ ] **Step 5: Run existing IQ-TREE and CF test suites to verify zero regression**

```bash
pytest tests/tree/test_ml_iqtree.py -x -q --tb=short
pytest tests/tree/test_cf.py -x -q --tb=short
```

Expected: both PASS.

- [ ] **Step 6: [skip commit]** — final commit will happen at Task 9

---

### Task 2: Create posttree package skeleton

**Files:**
- Create: `phyloai/posttree/__init__.py`
- Create: `tests/posttree/__init__.py`

- [ ] **Step 1: Create `phyloai/posttree/__init__.py`**

```python
"""PhyloAI post-tree analysis module."""
```

- [ ] **Step 2: Create `tests/posttree/__init__.py`**

```python
"""Tests for posttree analysis modules."""
```

- [ ] **Step 3: [skip commit]**

---

### Task 3: Implement `run_topology()` library function with input validation

**Files:**
- Create: `phyloai/posttree/topology.py`

- [ ] **Step 1: Write the complete file**

```python
"""Tree topology tests (AU / KH / SH / WKH / WSH / c-ELW) via IQ-TREE."""
from __future__ import annotations

import shlex
import subprocess
import time as _time
from pathlib import Path
from typing import Any

from phyloai.core.iqtree import (
    _detect_iqtree_version,
    _resolve_iqtree_path,
    IQTREE_COMPATIBLE_EXTENSIONS,
)


# ===================================================================
# Library-layer input validation
# ===================================================================

def _validate_inputs(
    *,
    matrix: Path,
    candidate_trees: list[Path],
    replicates: int,
    threads: int,
    overwrite: bool,
    resume: bool,
    model_expr: str | None,
    partitions: str | None,
    tool_args: str | None,
    guide_tree: str | None,
) -> list[str]:
    """Validate all spec-mandated inputs.  Returns a list of error messages.

    Covers: matrix existence, candidate-tree existence / non-empty / readable,
    numeric bounds, parameter mutual-exclusion, model-source completeness.
    Output-directory lifecycle (conflict / delete / resume-dir) is left to
    the CLI layer.
    """
    errors: list[str] = []

    import os as _os

    if not matrix.exists():
        errors.append(f"--matrix does not exist: {matrix}")
    elif not matrix.is_file():
        errors.append(f"--matrix is not a regular file: {matrix}")
    elif not _os.access(str(matrix), _os.R_OK):
        errors.append(f"--matrix is not readable: {matrix}")

    if not candidate_trees:
        errors.append("At least one --candidate-trees value is required")
    else:
        for i, ct in enumerate(candidate_trees):
            if not ct.exists():
                errors.append(f"--candidate-trees #{i + 1} does not exist: {ct}")
            elif not ct.is_file():
                errors.append(f"--candidate-trees #{i + 1} is not a regular file: {ct}")
            elif ct.stat().st_size == 0:
                errors.append(f"--candidate-trees #{i + 1} is empty: {ct}")
            elif not _os.access(str(ct), _os.R_OK):
                errors.append(f"--candidate-trees #{i + 1} is not readable: {ct}")

    if replicates < 1000:
        errors.append(f"--replicates must be >= 1000, got {replicates}")

    if threads < 1:
        errors.append(f"--threads must be >= 1, got {threads}")

    if overwrite and resume:
        errors.append("--overwrite and --resume are mutually exclusive")

    has_explicit_model = model_expr is not None or partitions is not None
    has_tool_args_model = False
    if tool_args:
        tokens = shlex.split(tool_args)
        has_tool_args_model = "-m" in tokens or "-p" in tokens
    if not has_explicit_model and not has_tool_args_model:
        errors.append(
            "Neither --model-expr, --partitions, nor -m/-p in --tool-args provided. "
            "Must specify one model source."
        )
    if model_expr and partitions:
        errors.append("--model-expr and --partitions are mutually exclusive")

    # Cross-source model conflict: if a high-level model-source flag is
    # given, --tool-args must not contain the OTHER model-source flag,
    # otherwise the effective command would have both -m and -p.
    if tool_args:
        tokens = shlex.split(tool_args)
        if model_expr and "-p" in tokens:
            errors.append(
                "--model-expr is set but --tool-args contains -p. "
                "Remove --model-expr if you want -p from --tool-args to take effect."
            )
        if partitions and "-m" in tokens:
            errors.append(
                "--partitions is set but --tool-args contains -m. "
                "Remove --partitions if you want -m from --tool-args to take effect."
            )

    if guide_tree:
        gt = Path(guide_tree)
        if not gt.exists():
            errors.append(f"--guide-tree does not exist: {guide_tree}")
        elif not gt.is_file():
            errors.append(f"--guide-tree is not a regular file: {guide_tree}")
        elif not _os.access(str(gt), _os.R_OK):
            errors.append(f"--guide-tree is not readable: {guide_tree}")

    return errors


# ===================================================================
# Candidate tree merging
# ===================================================================

def _merge_candidate_trees(
    tree_files: list[Path],
    output_dir: Path,
) -> Path:
    """Merge multiple NEWICK tree files into one candidate.trees file.

    Each input file should contain one NEWICK tree. Trees are concatenated
    in order, one per line.
    """
    merged_path = output_dir / "candidate.trees"
    with open(merged_path, "w") as out:
        for tf in tree_files:
            text = tf.read_text().strip()
            if not text:
                raise ValueError(f"Empty candidate tree file: {tf}")
            out.write(text)
            if not text.endswith("\n"):
                out.write("\n")
    return merged_path


# ===================================================================
# Command builder
# ===================================================================

_TOPOLOGY_BLOCKED_FLAGS = frozenset({"-s", "-z"})
_TOPOLOGY_BLOCKED_IO_CHARS = frozenset({"<", ">", "|"})


def _check_managed_flag_conflict(
    tool_args: str,
    *,
    blocked_flags: frozenset[str],
) -> None:
    """Reject BLOCKED flags and I/O redirects in --tool-args."""
    tokens = shlex.split(tool_args)
    for token in tokens:
        if token in blocked_flags:
            raise ValueError(f"Blocked managed flag in --tool-args: {token}")
        if any(c in token for c in _TOPOLOGY_BLOCKED_IO_CHARS):
            raise ValueError(f"Blocked I/O override in --tool-args: {token}")


def _is_flag_overridden(flag: str, tool_tokens: set[str]) -> bool:
    return flag in tool_tokens


def _build_topology_cmd(
    *,
    executable: str,
    matrix: Path,
    candidate_trees: Path,
    prefix: str,
    model_expr: str | None,
    partitions: str | None,
    guide_tree: str | None,
    replicates: int,
    threads: int,
    tool_args: str | None,
) -> list[str]:
    """Build the IQ-TREE topology test command line.

    Assembly order per design spec section 6.3:
    1. executable, -s <matrix>, -z <candidate-trees>
    2. --prefix (suppress-if-present)
    3. model source: -m or -p (suppress-if-present)
    4. -ft <guide-tree> (suppress-if-present)
    5. topology-test defaults: -n 0 -zb <replicates> -zw -au (each suppressible)
    6. -T <threads> (suppress-if-present)
    7. raw --tool-args appended last
    """
    cmd = [executable]
    tool_tokens = set(shlex.split(tool_args)) if tool_args else set()

    if tool_args:
        _check_managed_flag_conflict(
            tool_args, blocked_flags=_TOPOLOGY_BLOCKED_FLAGS,
        )

    # 1. Required inputs (never overridable)
    cmd.extend(["-s", str(matrix)])
    cmd.extend(["-z", str(candidate_trees)])

    # 2. Prefix
    if not _is_flag_overridden("--prefix", tool_tokens):
        cmd.extend(["--prefix", prefix])

    # 3. Model source
    if model_expr and not _is_flag_overridden("-m", tool_tokens):
        cmd.extend(["-m", model_expr])
    elif partitions and not _is_flag_overridden("-p", tool_tokens):
        cmd.extend(["-p", partitions])

    # 4. Guide tree
    if guide_tree and not _is_flag_overridden("-ft", tool_tokens):
        cmd.extend(["-ft", guide_tree])

    # 5. Topology test defaults
    if not _is_flag_overridden("-n", tool_tokens):
        cmd.extend(["-n", "0"])
    if not _is_flag_overridden("-zb", tool_tokens):
        cmd.extend(["-zb", str(replicates)])
    if not _is_flag_overridden("-zw", tool_tokens):
        cmd.append("-zw")
    if not _is_flag_overridden("-au", tool_tokens):
        cmd.append("-au")

    # 6. Threads
    if not _is_flag_overridden("-T", tool_tokens):
        cmd.extend(["-T", str(threads)])

    # 7. Raw tool-args
    if tool_args:
        cmd.extend(shlex.split(tool_args))

    return cmd


# ===================================================================
# USER TREES table parser (header-offset tolerant)
# ===================================================================

# Ordered list of columns that IQ-TREE typically emits in the USER TREES
# section.  Each entry is (canonical_name, header_token).  The parser uses
# the header line to discover offsets and then slices data rows by column
# position, not by shlex splitting.  This is robust across IQ-TREE versions
# with varying column spacing and additional weighted columns.
_USER_TREES_COLUMNS = [
    ("tree_id", "Tree"),
    ("log_likelihood", "logL"),
    ("delta_likelihood", "deltaL"),
    ("bp_rell", "bp-RELL"),
    ("p_kh", "p-KH"),
    ("p_sh", "p-SH"),
    ("p_wkh", "p-WKH"),
    ("p_wsh", "p-WSH"),
    ("c_elw", "c-ELW"),
    ("p_au", "p-AU"),
]


def _parse_user_trees_table(iqtree_path: Path) -> tuple[list[dict[str, Any]], list[str]]:
    """Parse the USER TREES section of an .iqtree report.

    Returns (tests, warnings).  Each test dict maps column names to parsed
    values; absent columns are stored as None.  Each test dict also contains
    a ``raw_line`` key with the original text of the data row.
    """
    if not iqtree_path.exists():
        return [], [f".iqtree file not found: {iqtree_path}"]

    text = iqtree_path.read_text()
    idx = text.find("USER TREES")
    if idx < 0:
        return [], [".iqtree file contains no USER TREES section"]

    lines = text[idx:].splitlines()

    # Find header row — the first line containing "Tree", "logL", "deltaL"
    header_idx = None
    for i, line in enumerate(lines):
        if "Tree" in line and "logL" in line and "deltaL" in line:
            header_idx = i
            break
    if header_idx is None:
        return [], ["Could not find USER TREES table header in .iqtree"]

    header_line = lines[header_idx]

    # Build a mapping: canonical_name → (start_pos, end_pos) based on
    # where each header token appears in the header line.
    col_spans: dict[str, tuple[int, int]] = {}
    for canonical, token in _USER_TREES_COLUMNS:
        pos = header_line.find(token)
        if pos >= 0:
            end = pos + len(token)
            col_spans[canonical] = (pos, end)

    if not col_spans:
        return [], ["No known column tokens found in USER TREES header"]

    # Sort columns by their start position for ordered slicing
    sorted_cols = sorted(col_spans.items(), key=lambda kv: kv[1][0])

    # For each data column, compute a value-extraction span:
    # from its end position to the start position of the next column
    # (or end of line for the last column).
    extractors: list[tuple[str, int]] = []
    for i, (name, (start, end)) in enumerate(sorted_cols):
        if i + 1 < len(sorted_cols):
            next_start = sorted_cols[i + 1][1][0]
        else:
            next_start = None
        extractors.append((name, end, next_start))

    # Parse data rows
    tests: list[dict[str, Any]] = []
    warnings: list[str] = []

    for line in lines[header_idx + 1:]:
        stripped = line.strip()
        if not stripped:
            break  # blank line ends the table
        if all(c in "-= " for c in stripped):
            continue  # separator line

        row: dict[str, Any] = {"raw_line": line}
        # Pre-fill all known columns (and sign variants) with None so that
        # absent columns never cause KeyError — the caller always sees None.
        for canonical, _token in _USER_TREES_COLUMNS:
            row[canonical] = None
            if canonical not in ("tree_id", "log_likelihood", "delta_likelihood"):
                row[canonical + "_sign"] = None

        for name, slice_start, slice_end in extractors:
            if slice_end is not None:
                token = line[slice_start:slice_end].strip()
            else:
                token = line[slice_start:].strip()

            # Split token into value and optional sign for p-value columns
            if name in ("bp_rell", "p_kh", "p_sh", "p_wkh", "p_wsh",
                        "c_elw", "p_au"):
                parts = token.split()
                if parts:
                    try:
                        row[name] = float(parts[0])
                    except ValueError:
                        row[name] = None
                    row[name + "_sign"] = parts[1] if len(parts) > 1 and parts[1] in ("+", "-") else None
                else:
                    row[name] = None
                    row[name + "_sign"] = None
            elif name == "tree_id":
                try:
                    row[name] = int(token)
                except ValueError:
                    row[name] = None
            else:
                try:
                    row[name] = float(token)
                except ValueError:
                    row[name] = None

        tests.append(row)

    return tests, warnings


# ===================================================================
# Output-file discovery
# ===================================================================

def _discover_optimized_trees(output_dir: Path, prefix: str) -> str | None:
    """Find IQ-TREE-optimized candidate tree file by prefix and known suffix patterns."""
    for suffix in (".treels.trees", ".trees", ".treels"):
        path = output_dir / f"{prefix}{suffix}"
        if path.exists():
            return str(path.name)
    return None


# ===================================================================
# result.json assembly
# ===================================================================

def _assemble_topology_result(
    *,
    run_start: float,
    iqtree_exe: str,
    tool_versions: dict[str, str],
    params: dict[str, Any],
    matrix: Path,
    candidate_trees_raw: list[Path],
    candidate_trees_mode: str,
    candidate_trees_effective: Path,
    effective_cmd: list[str],
    tool_stderr: str,
    output_dir: Path,
    prefix: str,
    returncode: int,
    dry_run: bool,
    warnings: list[str],
    error_category: str | None = None,
) -> dict[str, Any]:
    """Assemble the result.json payload.

    error_category: None for success, 'input' for validation errors,
    'env' for missing executable, 'tool' for IQ-TREE execution failure.
    """
    wall_time = round(_time.time() - run_start, 2)

    # Reconstruct CLI command string from resolved params
    cmd_parts = ["phyloai", "posttree", "topology"]
    cmd_parts.extend(["--matrix", str(matrix)])
    for ct in candidate_trees_raw:
        cmd_parts.extend(["--candidate-trees", str(ct)])
    if params.get("input_format") and params["input_format"] != "auto":
        cmd_parts.extend(["--input-format", params["input_format"]])
    if params.get("model_expr"):
        cmd_parts.extend(["--model-expr", params["model_expr"]])
    if params.get("partitions"):
        cmd_parts.extend(["--partitions", params["partitions"]])
    if params.get("guide_tree"):
        cmd_parts.extend(["--guide-tree", params["guide_tree"]])
    if params.get("replicates") != 10000:
        cmd_parts.extend(["--replicates", str(params["replicates"])])
    if params.get("prefix") and params["prefix"] != matrix.stem:
        cmd_parts.extend(["--prefix", params["prefix"]])
    cmd_parts.extend(["-o", str(output_dir)])
    if params.get("threads") != 4:
        cmd_parts.extend(["-t", str(params["threads"])])
    if params.get("iqtree_path"):
        cmd_parts.extend(["--iqtree-path", params["iqtree_path"]])
    if params.get("tool_args"):
        cmd_parts.extend(["--tool-args", params["tool_args"]])
    if params.get("overwrite"):
        cmd_parts.append("--overwrite")
    if params.get("resume"):
        cmd_parts.append("--resume")
    if params.get("dry_run"):
        cmd_parts.append("--dry-run")
    if params.get("quiet"):
        cmd_parts.append("-q")
    command = " ".join(cmd_parts)

    iqtree_log = output_dir / f"{prefix}.log"
    iqtree_report = output_dir / f"{prefix}.iqtree"

    # Parse USER TREES table
    tests: list[dict[str, Any]] = []
    parse_warnings: list[str] = []
    parsed_ok = True
    if iqtree_report.exists() and not dry_run:
        tests, parse_warnings = _parse_user_trees_table(iqtree_report)
        if parse_warnings:
            warnings.extend(parse_warnings)
        if not tests and returncode == 0:
            warnings.append("USER TREES parsing produced no rows; IQ-TREE completed successfully")

    # Discover optimized trees
    optimized_trees = _discover_optimized_trees(output_dir, prefix)

    # Key results
    n_candidate = len(tests) if tests else None
    best_tree_id = None
    n_rejected_au = None
    if tests:
        candidates_with_delta = [
            t for t in tests if t.get("delta_likelihood") is not None
        ]
        if candidates_with_delta:
            best = min(candidates_with_delta, key=lambda t: t["delta_likelihood"])
            best_tree_id = best.get("tree_id")
        au_vals = [
            t.get("p_au") for t in tests
            if t.get("p_au") is not None
        ]
        if au_vals:
            n_rejected_au = sum(1 for v in au_vals if v < 0.05)

    # Model source tracking — reflect which source actually supplied the
    # -m/-p flag in the final IQ-TREE command, accounting for suppress-if-present.
    _model_expr = params.get("model_expr")
    _partitions = params.get("partitions")
    _tool_args = params.get("tool_args")
    tool_tokens_for_source = set(shlex.split(_tool_args)) if _tool_args else set()
    if _model_expr and _is_flag_overridden("-m", tool_tokens_for_source):
        model_source = "tool-args"
    elif _partitions and _is_flag_overridden("-p", tool_tokens_for_source):
        model_source = "tool-args"
    elif _model_expr:
        model_source = "model-expr"
    elif _partitions:
        model_source = "partitions"
    else:
        model_source = "tool-args"

    key_results: dict[str, Any] = {
        "n_candidate_trees": n_candidate,
        "best_tree_id": best_tree_id,
        "n_rejected_au_0_05": n_rejected_au,
        "replicates": params["replicates"],
        "model_source": model_source,
    }

    data: dict[str, Any] = {
        "cmd": effective_cmd,
        "tool_stderr": tool_stderr,
        "log_iqtree": str(iqtree_report.name) if iqtree_report.exists() else None,
        "tool_log": str(iqtree_log.name) if iqtree_log.exists() else None,
        "optimized_trees": optimized_trees,
        "merged_candidate_trees": (
            str(candidate_trees_effective.name)
            if candidate_trees_mode == "individual-files"
            else None
        ),
        "tests": tests,
        "warnings": warnings,
    }

    return {
        "status": "success" if returncode == 0 else "error",
        "command": command,
        "wall_time": wall_time,
        "tool_versions": tool_versions,
        "params": params,
        "key_results": key_results,
        "error": None if returncode == 0 else f"IQ-TREE exited with code {returncode}",
        "error_category": error_category,
        "data": data,
    }


# ===================================================================
# Main entry point
# ===================================================================

def run_topology(
    *,
    matrix: Path,
    candidate_trees: list[Path],
    input_format: str = "auto",
    model_expr: str | None = None,
    partitions: str | None = None,
    guide_tree: str | None = None,
    replicates: int = 10000,
    prefix: str | None = None,
    output_dir: Path | None = None,
    threads: int = 4,
    iqtree_path: str | None = None,
    tool_args: str | None = None,
    overwrite: bool = False,
    resume: bool = False,
    dry_run: bool = False,
    quiet: bool = False,
) -> dict[str, Any]:
    """Run IQ-TREE tree topology tests.

    Returns a result.json payload dict.
    """
    run_start = _time.time()

    # --- Library-layer validation ---
    errors = _validate_inputs(
        matrix=matrix,
        candidate_trees=candidate_trees,
        replicates=replicates,
        threads=threads,
        overwrite=overwrite,
        resume=resume,
        model_expr=model_expr,
        partitions=partitions,
        tool_args=tool_args,
        guide_tree=guide_tree,
    )
    if errors:
        return {
            "status": "error",
            "command": "",
            "wall_time": 0.0,
            "tool_versions": {},
            "params": {},
            "key_results": {},
            "error": "; ".join(errors),
            "error_category": "input",
            "data": {"cmd": [], "tool_stderr": "", "tests": [], "warnings": errors},
        }

    # --- Resolve paths ---
    if output_dir is None:
        output_dir = Path("runs/posttree/topology")
    output_dir = output_dir.resolve()
    if not dry_run:
        output_dir.mkdir(parents=True, exist_ok=True)

    matrix = matrix.resolve()

    candidate_trees_raw = [ct.resolve() for ct in candidate_trees]
    if len(candidate_trees) == 1:
        candidate_trees_mode = "tree-list"
        candidate_trees_effective = candidate_trees_raw[0]
    else:
        candidate_trees_mode = "individual-files"
        if not dry_run:
            candidate_trees_effective = _merge_candidate_trees(
                candidate_trees_raw, output_dir,
            )
        else:
            candidate_trees_effective = output_dir / "candidate.trees"

    resolved_prefix = prefix if prefix else matrix.stem

    # --- Resolve IQ-TREE ---
    try:
        iqtree_exe = _resolve_iqtree_path(iqtree_path, dry_run)
    except (ValueError, FileNotFoundError) as e:
        return {
            "status": "error",
            "command": "",
            "wall_time": 0.0,
            "tool_versions": {},
            "params": {},
            "key_results": {},
            "error": str(e),
            "error_category": "env",
            "data": {"cmd": [], "tool_stderr": "", "tests": [], "warnings": [str(e)]},
        }
    tool_versions = (
        _detect_iqtree_version(iqtree_exe)
        if not dry_run
        else {"iqtree3": "dry-run"}
    )

    # --- Params dict (complete, all parameters) ---
    params: dict[str, Any] = {
        "matrix": str(matrix),
        "candidate_trees": [str(ct) for ct in candidate_trees_raw],
        "candidate_trees_mode": candidate_trees_mode,
        "candidate_trees_effective": str(candidate_trees_effective),
        "input_format": input_format,
        "model_expr": model_expr,
        "partitions": partitions,
        "guide_tree": guide_tree,
        "replicates": replicates,
        "prefix": resolved_prefix,
        "output_dir": str(output_dir),
        "threads": threads,
        "iqtree_path": iqtree_path,
        "tool_args": tool_args,
        "overwrite": overwrite,
        "resume": resume,
        "dry_run": dry_run,
        "quiet": quiet,
    }

    # --- Build IQ-TREE command ---
    effective_cmd = _build_topology_cmd(
        executable=iqtree_exe,
        matrix=matrix,
        candidate_trees=candidate_trees_effective,
        prefix=resolved_prefix,
        model_expr=model_expr,
        partitions=partitions,
        guide_tree=guide_tree,
        replicates=replicates,
        threads=threads,
        tool_args=tool_args,
    )

    if dry_run:
        return _assemble_topology_result(
            run_start=run_start,
            iqtree_exe=iqtree_exe,
            tool_versions=tool_versions,
            params=params,
            matrix=matrix,
            candidate_trees_raw=candidate_trees_raw,
            candidate_trees_mode=candidate_trees_mode,
            candidate_trees_effective=candidate_trees_effective,
            effective_cmd=effective_cmd,
            tool_stderr="",
            output_dir=output_dir,
            prefix=resolved_prefix,
            returncode=0,
            dry_run=True,
            warnings=[],
        )

    # --- Execute IQ-TREE ---
    try:
        proc = subprocess.run(
            effective_cmd,
            capture_output=True,
            text=True,
            cwd=str(output_dir),
        )
    except FileNotFoundError:
        return _assemble_topology_result(
            run_start=run_start,
            iqtree_exe=iqtree_exe,
            tool_versions=tool_versions,
            params=params,
            matrix=matrix,
            candidate_trees_raw=candidate_trees_raw,
            candidate_trees_mode=candidate_trees_mode,
            candidate_trees_effective=candidate_trees_effective,
            effective_cmd=effective_cmd,
            tool_stderr=f"Executable not found: {iqtree_exe}",
            output_dir=output_dir,
            prefix=resolved_prefix,
            returncode=3,
            dry_run=False,
            warnings=[],
            error_category="env",
        )

    merged_output = ""
    if proc.stdout:
        merged_output += proc.stdout
    if proc.stderr:
        if merged_output and not merged_output.endswith("\n"):
            merged_output += "\n"
        merged_output += proc.stderr

    return _assemble_topology_result(
        run_start=run_start,
        iqtree_exe=iqtree_exe,
        tool_versions=tool_versions,
        params=params,
        matrix=matrix,
        candidate_trees_raw=candidate_trees_raw,
        candidate_trees_mode=candidate_trees_mode,
        candidate_trees_effective=candidate_trees_effective,
        effective_cmd=effective_cmd,
        tool_stderr=merged_output,
        output_dir=output_dir,
        prefix=resolved_prefix,
        returncode=proc.returncode,
        dry_run=False,
        warnings=[],
        error_category="tool" if proc.returncode != 0 else None,
    )
```

- [ ] **Step 2: Syntax check**

```bash
python -c "from phyloai.posttree.topology import run_topology; print('OK')"
```

Expected: `OK`

- [ ] **Step 3: [skip commit]**

---

### Task 4: Write library-level unit tests

**Files:**
- Create: `tests/posttree/test_topology.py`

- [ ] **Step 1: Write tests**

```python
"""Tests for phyloai.posttree.topology."""
from __future__ import annotations

import shutil
from pathlib import Path

import pytest

from phyloai.posttree.topology import (
    _build_topology_cmd,
    _merge_candidate_trees,
    _parse_user_trees_table,
    _validate_inputs,
    run_topology,
)


# ------------------------------------------------------------------
# _validate_inputs
# ------------------------------------------------------------------

class TestValidateInputs:
    def test_matrix_does_not_exist(self, tmp_path: Path) -> None:
        errs = _validate_inputs(
            matrix=tmp_path / "nope.fa",
            candidate_trees=[tmp_path / "t.nwk"],
            replicates=10000, threads=4,
            overwrite=False, resume=False,
            model_expr="LG+F+R4", partitions=None,
            tool_args=None, guide_tree=None,
        )
        assert any("matrix" in e.lower() for e in errs)

    def test_candidate_tree_empty(self, tmp_path: Path) -> None:
        matrix = tmp_path / "m.fa"; matrix.write_text(">a\nMKT\n")
        ct = tmp_path / "empty.nwk"; ct.write_text("")
        errs = _validate_inputs(
            matrix=matrix, candidate_trees=[ct],
            replicates=10000, threads=4,
            overwrite=False, resume=False,
            model_expr="LG+F+R4", partitions=None,
            tool_args=None, guide_tree=None,
        )
        assert any("empty" in e for e in errs)

    def test_replicates_too_low(self, tmp_path: Path) -> None:
        matrix = tmp_path / "m.fa"; matrix.write_text(">a\nMKT\n")
        ct = tmp_path / "t.nwk"; ct.write_text("(a,b);\n")
        errs = _validate_inputs(
            matrix=matrix, candidate_trees=[ct],
            replicates=999, threads=4,
            overwrite=False, resume=False,
            model_expr="LG+F+R4", partitions=None,
            tool_args=None, guide_tree=None,
        )
        assert any("replicates" in e.lower() for e in errs)

    def test_overwrite_and_resume_mutually_exclusive(self, tmp_path: Path) -> None:
        matrix = tmp_path / "m.fa"; matrix.write_text(">a\nMKT\n")
        ct = tmp_path / "t.nwk"; ct.write_text("(a,b);\n")
        errs = _validate_inputs(
            matrix=matrix, candidate_trees=[ct],
            replicates=10000, threads=4,
            overwrite=True, resume=True,
            model_expr="LG+F+R4", partitions=None,
            tool_args=None, guide_tree=None,
        )
        assert any("mutually exclusive" in e.lower() for e in errs)

    def test_no_model_source(self, tmp_path: Path) -> None:
        matrix = tmp_path / "m.fa"; matrix.write_text(">a\nMKT\n")
        ct = tmp_path / "t.nwk"; ct.write_text("(a,b);\n")
        errs = _validate_inputs(
            matrix=matrix, candidate_trees=[ct],
            replicates=10000, threads=4,
            overwrite=False, resume=False,
            model_expr=None, partitions=None,
            tool_args=None, guide_tree=None,
        )
        assert any("model" in e.lower() for e in errs)

    def test_both_model_expr_and_partitions(self, tmp_path: Path) -> None:
        matrix = tmp_path / "m.fa"; matrix.write_text(">a\nMKT\n")
        ct = tmp_path / "t.nwk"; ct.write_text("(a,b);\n")
        errs = _validate_inputs(
            matrix=matrix, candidate_trees=[ct],
            replicates=10000, threads=4,
            overwrite=False, resume=False,
            model_expr="LG", partitions="m.nex",
            tool_args=None, guide_tree=None,
        )
        assert any("mutually exclusive" in e.lower() for e in errs)

    def test_all_valid(self, tmp_path: Path) -> None:
        matrix = tmp_path / "m.fa"; matrix.write_text(">a\nMKT\n")
        ct = tmp_path / "t.nwk"; ct.write_text("(a,b);\n")
        errs = _validate_inputs(
            matrix=matrix, candidate_trees=[ct],
            replicates=10000, threads=4,
            overwrite=False, resume=False,
            model_expr="LG+F+R4", partitions=None,
            tool_args=None, guide_tree=None,
        )
        assert errs == []


# ------------------------------------------------------------------
# _merge_candidate_trees
# ------------------------------------------------------------------

class TestMergeCandidateTrees:
    def test_merge_two_files(self, tmp_path: Path) -> None:
        (tmp_path / "h1.nwk").write_text("(A,B);\n")
        (tmp_path / "h2.nwk").write_text("(A,C);\n")
        merged = _merge_candidate_trees(
            [tmp_path / "h1.nwk", tmp_path / "h2.nwk"], tmp_path,
        )
        assert merged.name == "candidate.trees"
        content = merged.read_text()
        assert "(A,B);" in content
        assert "(A,C);" in content

    def test_empty_file_raises(self, tmp_path: Path) -> None:
        (tmp_path / "empty.nwk").write_text("")
        with pytest.raises(ValueError, match="Empty candidate tree file"):
            _merge_candidate_trees([tmp_path / "empty.nwk"], tmp_path)


# ------------------------------------------------------------------
# _build_topology_cmd
# ------------------------------------------------------------------

class TestBuildTopologyCmd:
    def test_basic_command(self, tmp_path: Path) -> None:
        cmd = _build_topology_cmd(
            executable="iqtree3",
            matrix=tmp_path / "matrix.fa",
            candidate_trees=tmp_path / "trees",
            prefix="matrix",
            model_expr="LG+F+R4",
            partitions=None,
            guide_tree=None,
            replicates=10000, threads=20,
            tool_args=None,
        )
        assert cmd[0] == "iqtree3"
        assert "-s" in cmd; assert "-z" in cmd
        assert "-m" in cmd and "LG+F+R4" in cmd
        assert "-n" in cmd and "0" in cmd
        assert "-zb" in cmd and "10000" in cmd
        assert "-zw" in cmd; assert "-au" in cmd
        assert "-T" in cmd and "20" in cmd

    def test_partitions_mode(self, tmp_path: Path) -> None:
        cmd = _build_topology_cmd(
            executable="iqtree3",
            matrix=tmp_path / "matrix.fa",
            candidate_trees=tmp_path / "trees",
            prefix="m",
            model_expr=None, partitions="m.best_model.nex",
            guide_tree=None,
            replicates=10000, threads=4,
            tool_args=None,
        )
        assert "-p" in cmd; assert "-m" not in cmd

    def test_guide_tree(self, tmp_path: Path) -> None:
        cmd = _build_topology_cmd(
            executable="iqtree3",
            matrix=tmp_path / "matrix.fa",
            candidate_trees=tmp_path / "trees",
            prefix="m",
            model_expr="LG+C20+F+R4", partitions=None,
            guide_tree="guide.nwk",
            replicates=10000, threads=4,
            tool_args=None,
        )
        assert "-ft" in cmd

    def test_suppress_threads_via_tool_args(self, tmp_path: Path) -> None:
        cmd = _build_topology_cmd(
            executable="iqtree3",
            matrix=tmp_path / "matrix.fa",
            candidate_trees=tmp_path / "trees",
            prefix="m",
            model_expr="LG+F+R4", partitions=None, guide_tree=None,
            replicates=10000, threads=20,
            tool_args="-T 30",
        )
        t_indices = [i for i, t in enumerate(cmd) if t == "-T"]
        assert len(t_indices) == 1
        assert cmd[t_indices[0] + 1] == "30"

    def test_suppress_zb_via_tool_args(self, tmp_path: Path) -> None:
        cmd = _build_topology_cmd(
            executable="iqtree3",
            matrix=tmp_path / "matrix.fa",
            candidate_trees=tmp_path / "trees",
            prefix="m",
            model_expr="LG+F+R4", partitions=None, guide_tree=None,
            replicates=10000, threads=4,
            tool_args="-zb 5000",
        )
        zb_indices = [i for i, t in enumerate(cmd) if t == "-zb"]
        assert len(zb_indices) == 1
        assert cmd[zb_indices[0] + 1] == "5000"

    def test_blocked_s_flag_raises(self, tmp_path: Path) -> None:
        with pytest.raises(ValueError, match="-s"):
            _build_topology_cmd(
                executable="iqtree3",
                matrix=tmp_path / "matrix.fa",
                candidate_trees=tmp_path / "trees",
                prefix="m",
                model_expr="LG+F+R4", partitions=None, guide_tree=None,
                replicates=10000, threads=4,
                tool_args="-s other.fa",
            )

    def test_blocked_z_flag_raises(self, tmp_path: Path) -> None:
        with pytest.raises(ValueError, match="-z"):
            _build_topology_cmd(
                executable="iqtree3",
                matrix=tmp_path / "matrix.fa",
                candidate_trees=tmp_path / "trees",
                prefix="m",
                model_expr="LG+F+R4", partitions=None, guide_tree=None,
                replicates=10000, threads=4,
                tool_args="-z other.trees",
            )

    def test_replicates_value(self, tmp_path: Path) -> None:
        cmd = _build_topology_cmd(
            executable="iqtree3",
            matrix=tmp_path / "matrix.fa",
            candidate_trees=tmp_path / "trees",
            prefix="m",
            model_expr="LG+F+R4", partitions=None, guide_tree=None,
            replicates=2000, threads=4,
            tool_args=None,
        )
        zb_idx = cmd.index("-zb")
        assert cmd[zb_idx + 1] == "2000"


# ------------------------------------------------------------------
# _parse_user_trees_table
# ------------------------------------------------------------------

class TestParseUserTreesTable:
    def test_parse_standard_iqtree_output(self, tmp_path: Path) -> None:
        content = """\
USER TREES:

See http://www.iqtree.org/doc/Topology-Tests

Tree      logL    deltaL  bp-RELL    p-KH     p-SH    p-WKH    p-WSH    c-ELW       p-AU
------------------------------------------------------------------------------------------
 1  -21152.617    0.000  0.7110 +  0.7400 +  1.0000 +  0.7380 +  1.0000 +  0.6954 +  0.7939 +
 2  -21158.123    5.506  0.2299 +  0.2590 +  0.1260 +  0.2690 +  0.1330 +  0.2275 +  0.2336 +
 3  -21162.987   10.370  0.0392 +  0.0080 -  0.0070 -  0.0080 -  0.0060 -  0.0404 +  0.0140 -
 4  -21168.456   15.839  0.0199 -  0.0010 -  0.0000 -  0.0010 -  0.0000 -  0.0367 -  0.0030 -

------------------------------------------------------------------------------------------
"""
        p = tmp_path / "test.iqtree"
        p.write_text(content)
        tests, warnings = _parse_user_trees_table(p)
        assert len(tests) == 4
        t1 = tests[0]
        assert t1["tree_id"] == 1
        assert t1["log_likelihood"] == -21152.617
        assert t1["bp_rell"] == 0.7110
        assert t1["bp_rell_sign"] == "+"
        assert t1["p_au"] == 0.7939
        # raw_line preserved
        assert "raw_line" in t1
        assert "1  -21152.617" in t1["raw_line"]

    def test_missing_file(self, tmp_path: Path) -> None:
        tests, warnings = _parse_user_trees_table(tmp_path / "nope.iqtree")
        assert tests == []
        assert len(warnings) > 0

    def test_no_user_trees_section(self, tmp_path: Path) -> None:
        (tmp_path / "test.iqtree").write_text("OTHER SECTION\nNo trees\n")
        tests, warnings = _parse_user_trees_table(tmp_path / "test.iqtree")
        assert tests == []
        assert len(warnings) > 0

    def test_missing_column_stores_none(self, tmp_path: Path) -> None:
        # IQ-TREE may omit some columns; absent values must be None
        content = """\
USER TREES:
Tree      logL    deltaL
 1  -100.0  0.0
 2  -105.0  5.0
"""
        (tmp_path / "test.iqtree").write_text(content)
        tests, _ = _parse_user_trees_table(tmp_path / "test.iqtree")
        assert len(tests) == 2
        assert tests[0]["tree_id"] == 1
        assert tests[0]["p_au"] is None  # column absent entirely


# ------------------------------------------------------------------
# run_topology (dry-run + integration)
# ------------------------------------------------------------------

class TestRunTopology:
    def test_dry_run_single_tree_file(self, tmp_path: Path) -> None:
        from tests.helpers import validate_result_json

        (tmp_path / "matrix.fa").write_text(">a\nMKTLLL\n>b\nMKTLLL\n")
        (tmp_path / "trees").write_text("(a,b);\n")

        result = run_topology(
            matrix=tmp_path / "matrix.fa",
            candidate_trees=[tmp_path / "trees"],
            model_expr="LG+F+R4",
            output_dir=tmp_path / "out",
            dry_run=True,
        )
        validate_result_json(result)
        assert result["status"] == "success"
        assert result["params"]["candidate_trees_mode"] == "tree-list"
        assert isinstance(result["data"]["cmd"], list)

    def test_dry_run_multiple_tree_files(self, tmp_path: Path) -> None:
        (tmp_path / "matrix.fa").write_text(">a\nMKTLLL\n>b\nMKTLLL\n")
        (tmp_path / "h1.nwk").write_text("(a,b);\n")
        (tmp_path / "h2.nwk").write_text("(a,c);\n")

        result = run_topology(
            matrix=tmp_path / "matrix.fa",
            candidate_trees=[tmp_path / "h1.nwk", tmp_path / "h2.nwk"],
            model_expr="LG+F+R4",
            output_dir=tmp_path / "out",
            dry_run=True,
        )
        assert result["params"]["candidate_trees_mode"] == "individual-files"
        assert len(result["params"]["candidate_trees"]) == 2

    def test_validation_error_returns_error_payload(self, tmp_path: Path) -> None:
        result = run_topology(
            matrix=tmp_path / "nope.fa",
            candidate_trees=[],
            model_expr="LG+F+R4",
            dry_run=True,
        )
        assert result["status"] == "error"
        assert result["error"] is not None

    @pytest.mark.skipif(
        not shutil.which("iqtree3"),
        reason="iqtree3 not found in PATH",
    )
    def test_run_topology_real_iqtree(self, tmp_path: Path) -> None:
        from tests.helpers import validate_result_json

        matrix = tmp_path / "matrix.fa"
        matrix.write_text(
            ">t1\nMKTLLLTLWVV\n>t2\nMKTLLLTLWVI\n>t3\nMKTLLLSLWVI\n>t4\nMKTLLLTLWVA\n"
        )
        (tmp_path / "t1.nwk").write_text("(t1,t2,(t3,t4));\n")
        (tmp_path / "t2.nwk").write_text("(t1,t3,(t2,t4));\n")
        (tmp_path / "t3.nwk").write_text("(t1,t4,(t2,t3));\n")
        (tmp_path / "t4.nwk").write_text("(t2,t3,(t1,t4));\n")

        out = tmp_path / "out"
        result = run_topology(
            matrix=matrix,
            candidate_trees=[
                tmp_path / "t1.nwk", tmp_path / "t2.nwk",
                tmp_path / "t3.nwk", tmp_path / "t4.nwk",
            ],
            model_expr="LG",
            replicates=1000, output_dir=out, threads=1,
        )

        validate_result_json(result)
        assert result["status"] == "success"
        assert "iqtree3" in result["tool_versions"]
        assert len(result["data"]["tests"]) == 4
        assert result["data"]["merged_candidate_trees"] == "candidate.trees"
        assert (out / "candidate.trees").exists()
        assert (out / "matrix.iqtree").exists()
        assert (out / "matrix.log").exists()
```

- [ ] **Step 2: Run library tests**

```bash
pytest tests/posttree/test_topology.py -v --tb=short
```

Expected: unit tests PASS; integration test PASS or SKIP.

- [ ] **Step 3: [skip commit]**

---

### Task 5: Create CLI command with 8-section help and manual input validation

**Files:**
- Create: `phyloai/cli/commands/posttree.py`

- [ ] **Step 1: Write the CLI file**

```python
"""Post-tree analysis CLI commands."""
from __future__ import annotations

import json
import shlex
import shutil
from pathlib import Path

import click

from phyloai.core.iqtree import IQTREE_COMPATIBLE_EXTENSIONS


class _PosttreeGroup(click.Group):
    def list_commands(self, ctx: click.Context) -> list[str]:
        return ["topology"]


@click.group(cls=_PosttreeGroup)
def posttree() -> None:
    """Post-tree analysis commands."""


def _fail(message: str, exit_code: int = 1) -> None:
    click.echo(f"Error: {message}", err=True)
    raise click.exceptions.Exit(exit_code)


_TOPOLOGY_HELP = """
Tree topology tests (AU / KH / SH / WKH / WSH / c-ELW).

Compares a set of candidate trees against a supermatrix alignment using
IQ-TREE's built-in topology test framework.

\\b
PURPOSE
  This command tests whether alternative topologies are significantly worse
  than the best-scoring candidate. It does NOT infer new trees — use
  `phyloai tree ml iqtree` for ML tree inference and model selection.

\\b
INPUT
  --matrix             Single supermatrix alignment (FASTA/PHYLIP/NEXUS/CLUSTAL).
  --candidate-trees    Accepts either one tree-list file (one NEWICK tree per
                       line) or multiple individual NEWICK tree files passed
                       in order. Multiple files are merged by PhyloAI into
                       candidate.trees before invoking IQ-TREE.

\\b
MODEL SOURCE  (exactly one required)
  --model-expr          Complete IQ-TREE -m expression (e.g. LG+F+R4, C20+F+R4).
  --partitions PATH     Previously optimized partition model (e.g. .best_model.nex).
                        Maps to IQ-TREE -p.
  --guide-tree PATH     Guide tree for PMSF models. Maps to IQ-TREE -ft.

  PhyloAI does NOT expose ModelFinder here because topology tests are run
  after model inference — you should already have a preferred model.

\\b
DEFAULT TESTS
  PhyloAI generates standard topology-test flags:

      -n 0 -zb <replicates> -zw -au

  This produces: bp-RELL, KH, SH, weighted KH, weighted SH, c-ELW, and AU.
  Each test is individually suppressible via --tool-args.

\\b
ADVANCED IQ-TREE ARGS
  --tool-args TEXT    Additional IQ-TREE strategy parameters. Blocked flags:
                      -s, -z (managed by --matrix and --candidate-trees).
  --iqtree-path PATH  Explicit path to iqtree3 executable.

  PhyloAI-built flags are suppressed when the same flag appears in
  --tool-args (suppress-if-present).  Overrideable: -m, -p, -ft, -n, -zb,
  -zw, -au, -T, --prefix.

\\b
EXAMPLES

  # Homogeneous unpartitioned model
  phyloai posttree topology --matrix raw.fa --candidate-trees trees \\
      --model-expr LG+F+R4 --replicates 10000 -t 20

  # Heterogeneous model
  phyloai posttree topology --matrix raw.fa --candidate-trees trees \\
      --model-expr C20+F+R4 -t 20

  # PMSF model with guide tree
  phyloai posttree topology --matrix raw.fa --candidate-trees trees \\
      --model-expr LG+C20+F+R4 --guide-tree guide.tree -t 4

  # Previously optimized partition model
  phyloai posttree topology --matrix raw.fa --candidate-trees trees \\
      --partitions raw.best_model.nex -t 20

  # Multiple individual tree files (merged by PhyloAI)
  phyloai posttree topology --matrix raw.fa \\
      --candidate-trees h1.nwk --candidate-trees h2.nwk \\
      --candidate-trees h3.nwk --model-expr LG+F+R4 -t 20

  # Custom exchangeabilities + site frequencies via --tool-args
  phyloai posttree topology --matrix raw.fa --candidate-trees trees \\
      --model-expr custom.exchangeabilities+R4 \\
      --tool-args "-fs custom.sitefreq" -t 30

\\b
INPUT FORMAT AND SEQUENCE TYPE
  --input-format only affects PhyloAI's own matrix preflight validation;
  it is NOT passed to IQ-TREE.  IQ-TREE's --seqtype flag can be passed
  via --tool-args when needed (e.g. --tool-args "--seqtype AA").

\\b
INTERPRETATION
  KH / SH / WKH / WSH / AU are p-values.  Trees with p < 0.05 are rejected
  by that test.  bp-RELL and c-ELW are weights (not p-values).  The AU test
  is generally considered the most reliable.
"""


@posttree.command("topology", help=_TOPOLOGY_HELP)
@click.option(
    "--matrix", type=click.Path(path_type=Path), default=None,
    help="Single supermatrix alignment (FASTA/PHYLIP/NEXUS/CLUSTAL).  Maps to IQ-TREE -s.",
)
@click.option(
    "--candidate-trees", "candidate_trees_raw",
    multiple=True, type=click.Path(path_type=Path),
    help=(
        "Candidate tree input. Accepts either one tree-list file (one NEWICK tree "
        "per line) or multiple individual NEWICK tree files (merged in order by PhyloAI)."
    ),
)
@click.option(
    "--input-format",
    type=click.Choice(["auto", "fasta", "phylip-relaxed", "nexus", "clustal"]),
    default="auto",
    help="PhyloAI-side matrix format hint for preflight validation. Not passed to IQ-TREE.",
)
@click.option(
    "--model-expr", type=str, default=None,
    help="Complete IQ-TREE -m model expression (e.g. LG+F+R4, C20+F+R4).",
)
@click.option(
    "--partitions", type=click.Path(path_type=Path), default=None,
    help="Previously optimized partition model. Maps to IQ-TREE -p.",
)
@click.option(
    "--guide-tree", type=click.Path(path_type=Path), default=None,
    help="Guide tree for PMSF-style model expressions. Maps to IQ-TREE -ft.",
)
@click.option(
    "--replicates", type=int, default=10000,
    help="RELL replicates (min 1000, default 10000). Maps to IQ-TREE -zb.",
)
@click.option(
    "--prefix", type=str, default=None,
    help="IQ-TREE output prefix (default: matrix file stem).",
)
@click.option(
    "-o", "--output-dir", type=click.Path(path_type=Path),
    default=Path("runs/posttree/topology"),
    help="Output directory.",
)
@click.option(
    "-t", "--threads", type=int, default=4,
    help="Thread count. Maps to IQ-TREE -T unless overridden by --tool-args.",
)
@click.option(
    "--iqtree-path", type=str, default=None,
    help="Explicit path to iqtree3 executable.",
)
@click.option(
    "--tool-args", type=str, default=None,
    help="Additional IQ-TREE strategy parameters. Blocked flags: -s, -z.",
)
@click.option("--overwrite", is_flag=True, default=False,
              help="Delete and recreate output directory.")
@click.option("--resume", is_flag=True, default=False,
              help="Reuse existing output directory with IQ-TREE native resume.")
@click.option("--dry-run", is_flag=True, default=False,
              help="Print the IQ-TREE command without executing.")
@click.option("-q", "--quiet", is_flag=True, default=False,
              help="Suppress terminal output except errors.")
def topology_command(
    matrix: Path | None,
    candidate_trees_raw: tuple[Path, ...],
    input_format: str,
    model_expr: str | None,
    partitions: Path | None,
    guide_tree: Path | None,
    replicates: int,
    prefix: str | None,
    output_dir: Path,
    threads: int,
    iqtree_path: str | None,
    tool_args: str | None,
    overwrite: bool,
    resume: bool,
    dry_run: bool,
    quiet: bool,
) -> None:
    """Tree topology tests (AU / KH / SH / WKH / WSH / c-ELW)."""
    from phyloai.posttree.topology import run_topology

    # ---- Manual validation (exit code 1 for all user input errors) ----

    if matrix is None:
        _fail("--matrix is required", exit_code=1)
    if not candidate_trees_raw:
        _fail("At least one --candidate-trees is required", exit_code=1)

    matrix_path: Path = matrix

    # Matrix extension
    ext = matrix_path.suffix.lower()
    if ext not in IQTREE_COMPATIBLE_EXTENSIONS:
        _fail(
            f"Unsupported matrix extension: {ext}. "
            f"Supported: {', '.join(sorted(IQTREE_COMPATIBLE_EXTENSIONS))}",
            exit_code=1,
        )

    # Candidate trees (existence / non-empty / readability)
    candidate_trees_list = list(candidate_trees_raw)
    for i, ct in enumerate(candidate_trees_list):
        if not ct.is_file():
            _fail(f"--candidate-trees #{i + 1} is not a regular file: {ct}", exit_code=1)
        if ct.stat().st_size == 0:
            _fail(f"--candidate-trees #{i + 1} is empty: {ct}", exit_code=1)

    # Model source
    has_explicit = model_expr is not None or partitions is not None
    has_tool_args_model = False
    if tool_args:
        tokens = shlex.split(tool_args)
        has_tool_args_model = "-m" in tokens or "-p" in tokens
    if not has_explicit and not has_tool_args_model:
        _fail(
            "Neither --model-expr, --partitions, nor -m/-p in --tool-args provided. "
            "Must specify one model source.",
            exit_code=1,
        )
    if model_expr and partitions:
        _fail("--model-expr and --partitions are mutually exclusive.", exit_code=1)

    # Cross-source model conflict: if a high-level model source is given,
    # --tool-args must not contain the other model-source flag.
    if tool_args:
        tokens = shlex.split(tool_args)
        if model_expr and "-p" in tokens:
            _fail(
                "--model-expr is set but --tool-args contains -p. "
                "Remove --model-expr if you want -p from --tool-args to take effect.",
                exit_code=1,
            )
        if partitions and "-m" in tokens:
            _fail(
                "--partitions is set but --tool-args contains -m. "
                "Remove --partitions if you want -m from --tool-args to take effect.",
                exit_code=1,
            )

    # partitions / guide-tree existence and readability (no longer validated by Click)
    if partitions:
        if not partitions.exists():
            _fail(f"--partitions does not exist: {partitions}", exit_code=1)
        if not partitions.is_file():
            _fail(f"--partitions is not a regular file: {partitions}", exit_code=1)
        if not __import__("os").access(str(partitions), __import__("os").R_OK):
            _fail(f"--partitions is not readable: {partitions}", exit_code=1)
    if guide_tree:
        if not guide_tree.exists():
            _fail(f"--guide-tree does not exist: {guide_tree}", exit_code=1)
        if not guide_tree.is_file():
            _fail(f"--guide-tree is not a regular file: {guide_tree}", exit_code=1)
        if not __import__("os").access(str(guide_tree), __import__("os").R_OK):
            _fail(f"--guide-tree is not readable: {guide_tree}", exit_code=1)

    # overwrite / resume
    if overwrite and resume:
        _fail("--overwrite and --resume are mutually exclusive.", exit_code=1)

    # Numeric bounds
    if replicates < 1000:
        _fail(f"--replicates must be >= 1000, got {replicates}", exit_code=1)
    if threads < 1:
        _fail(f"--threads must be >= 1, got {threads}", exit_code=1)

    # --tool-args blocked flags
    if tool_args:
        from phyloai.posttree.topology import _check_managed_flag_conflict
        try:
            _check_managed_flag_conflict(
                tool_args, blocked_flags=frozenset({"-s", "-z"}),
            )
        except ValueError as e:
            _fail(str(e), exit_code=1)

    # ---- Output directory lifecycle (CLI layer) ----

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

    # ---- Execute ----

    guide_tree_str = str(guide_tree) if guide_tree else None

    payload = run_topology(
        matrix=matrix_path,
        candidate_trees=candidate_trees_list,
        input_format=input_format,
        model_expr=model_expr,
        partitions=str(partitions) if partitions else None,
        guide_tree=guide_tree_str,
        replicates=replicates,
        prefix=prefix,
        output_dir=output_dir,
        threads=threads,
        iqtree_path=iqtree_path,
        tool_args=tool_args,
        overwrite=overwrite,
        resume=resume,
        dry_run=dry_run,
        quiet=quiet,
    )

    # ---- Write / display ----

    if not dry_run:
        result_path = output_dir / "result.json"
        result_path.parent.mkdir(parents=True, exist_ok=True)
        with open(result_path, "w") as fh:
            json.dump(payload, fh, indent=2)
        if not quiet:
            click.echo(f"Result written to {result_path}")

    if dry_run:
        cmd_str = " ".join(payload["data"]["cmd"])
        click.echo(f"Would run: {cmd_str}")
    elif payload["status"] == "error":
        err_msg = payload.get("error") or "Unknown error"
        cat = payload.get("error_category")
        if cat == "input":
            code = 1
        elif cat == "env":
            code = 3
        else:
            code = 2
        _fail(err_msg, exit_code=code)

    if not quiet:
        kr = payload["key_results"]
        click.echo(f"Status: {payload['status']}")
        click.echo(f"Wall time: {payload['wall_time']:.1f}s")
        if kr.get("n_candidate_trees"):
            click.echo(f"Candidate trees tested: {kr['n_candidate_trees']}")
        if kr.get("best_tree_id") is not None:
            click.echo(f"Best tree: #{kr['best_tree_id']}")
        if kr.get("n_rejected_au_0_05") is not None:
            click.echo(f"Rejected (AU < 0.05): {kr['n_rejected_au_0_05']}")
```

- [ ] **Step 2: Verify CLI help renders properly**

```bash
python -c "
from click.testing import CliRunner
from phyloai.cli.commands.posttree import posttree
r = CliRunner().invoke(posttree, ['topology', '--help'])
print(r.output)
print('Exit:', r.exit_code)
"
```

Expected: exit 0, output contains 8 section headers: PURPOSE, INPUT, MODEL SOURCE, DEFAULT TESTS, ADVANCED IQ-TREE ARGS, EXAMPLES, INPUT FORMAT AND SEQUENCE TYPE, INTERPRETATION. 6 example blocks present.

- [ ] **Step 3: [skip commit]**

---

### Task 6: Register posttree group in CLI main

**Files:**
- Modify: `phyloai/cli/main.py`

- [ ] **Step 1: Add import and registration**

In `phyloai/cli/main.py`, add after the existing imports:

```python
from phyloai.cli.commands.posttree import posttree
```

Add after the existing `cli.add_command(tree)`:

```python
cli.add_command(posttree)
```

- [ ] **Step 2: Verify CLI tree**

```bash
python -m phyloai.cli.main --help
```

Expected: shows `posttree` in command list.

```bash
python -m phyloai.cli.main posttree --help
```

Expected: shows `topology` subcommand.

```bash
python -m phyloai.cli.main posttree topology --help
```

Expected: full 8-section help with 6 examples.

- [ ] **Step 3: [skip commit]**

---

### Task 7: Write CLI validation and integration tests

**Files:**
- Create: `tests/cli/test_posttree_topology.py`

- [ ] **Step 1: Write CLI tests**

```python
"""CLI tests for phyloai posttree topology."""
from __future__ import annotations

import json
import shutil
from pathlib import Path

import pytest
from click.testing import CliRunner

from phyloai.cli.main import cli


@pytest.fixture
def runner() -> CliRunner:
    return CliRunner()


@pytest.fixture
def matrix_file(tmp_path: Path) -> Path:
    p = tmp_path / "matrix.fa"
    p.write_text(">a\nMKTLLL\n>b\nMKTLLL\n")
    return p


@pytest.fixture
def candidate_trees_file(tmp_path: Path) -> Path:
    p = tmp_path / "candidates.trees"
    p.write_text("(a,b);\n")
    return p


# ------------------------------------------------------------------
# Help content
# ------------------------------------------------------------------

class TestCLIHelp:
    def test_help_contains_eight_sections(self, runner: CliRunner) -> None:
        result = runner.invoke(cli, ["posttree", "topology", "--help"])
        assert result.exit_code == 0
        required_sections = [
            "PURPOSE", "INPUT", "MODEL SOURCE", "DEFAULT TESTS",
            "ADVANCED IQ-TREE ARGS", "EXAMPLES",
            "INPUT FORMAT AND SEQUENCE TYPE", "INTERPRETATION",
        ]
        for section in required_sections:
            assert section in result.output, f"Missing help section: {section}"

    def test_help_contains_all_six_examples(self, runner: CliRunner) -> None:
        result = runner.invoke(cli, ["posttree", "topology", "--help"])
        assert result.exit_code == 0
        assert "LG+F+R4" in result.output
        assert "C20+F+R4" in result.output
        assert "LG+C20+F+R4" in result.output
        assert "--partitions raw.best_model.nex" in result.output
        assert "h1.nwk --candidate-trees h2.nwk" in result.output
        assert "custom.exchangeabilities" in result.output

    def test_help_shows_all_options(self, runner: CliRunner) -> None:
        result = runner.invoke(cli, ["posttree", "topology", "--help"])
        assert result.exit_code == 0
        for opt in ("--matrix", "--candidate-trees", "--model-expr", "--partitions",
                     "--guide-tree", "--replicates", "--prefix", "--output-dir",
                     "--threads", "--iqtree-path", "--tool-args",
                     "--overwrite", "--resume", "--dry-run"):
            assert opt in result.output


# ------------------------------------------------------------------
# Input validation (all exit code 1 per spec)
# ------------------------------------------------------------------

class TestCLIValidation:
    def test_missing_matrix_exits_1(self, runner: CliRunner) -> None:
        result = runner.invoke(cli, ["posttree", "topology"])
        assert result.exit_code == 1
        assert "matrix" in result.output.lower()

    def test_missing_candidate_trees_exits_1(
        self, runner: CliRunner, matrix_file: Path,
    ) -> None:
        result = runner.invoke(cli, [
            "posttree", "topology", "--matrix", str(matrix_file),
        ])
        assert result.exit_code == 1
        assert "candidate" in result.output.lower()

    def test_missing_model_source_exits_1(
        self, runner: CliRunner, matrix_file: Path, candidate_trees_file: Path,
    ) -> None:
        result = runner.invoke(cli, [
            "posttree", "topology",
            "--matrix", str(matrix_file),
            "--candidate-trees", str(candidate_trees_file),
        ])
        assert result.exit_code == 1
        assert "model" in result.output.lower()

    def test_replicates_below_minimum(
        self, runner: CliRunner, matrix_file: Path, candidate_trees_file: Path,
    ) -> None:
        result = runner.invoke(cli, [
            "posttree", "topology",
            "--matrix", str(matrix_file),
            "--candidate-trees", str(candidate_trees_file),
            "--model-expr", "LG+F+R4",
            "--replicates", "999",
        ])
        assert result.exit_code == 1

    def test_overwrite_and_resume_mutually_exclusive(
        self, runner: CliRunner, matrix_file: Path, candidate_trees_file: Path,
    ) -> None:
        result = runner.invoke(cli, [
            "posttree", "topology",
            "--matrix", str(matrix_file),
            "--candidate-trees", str(candidate_trees_file),
            "--model-expr", "LG+F+R4",
            "--overwrite", "--resume",
        ])
        assert result.exit_code == 1
        assert "mutually exclusive" in result.output.lower()

    def test_both_model_expr_and_partitions(
        self, runner: CliRunner, matrix_file: Path, candidate_trees_file: Path,
        tmp_path: Path,
    ) -> None:
        (tmp_path / "m.best_model.nex").write_text("#nexus\n")
        result = runner.invoke(cli, [
            "posttree", "topology",
            "--matrix", str(matrix_file),
            "--candidate-trees", str(candidate_trees_file),
            "--model-expr", "LG+F+R4",
            "--partitions", str(tmp_path / "m.best_model.nex"),
        ])
        assert result.exit_code == 1
        assert "mutually exclusive" in result.output.lower()

    def test_blocked_s_in_tool_args(
        self, runner: CliRunner, matrix_file: Path, candidate_trees_file: Path,
    ) -> None:
        result = runner.invoke(cli, [
            "posttree", "topology",
            "--matrix", str(matrix_file),
            "--candidate-trees", str(candidate_trees_file),
            "--model-expr", "LG+F+R4",
            "--tool-args", "-s other.fa",
        ])
        assert result.exit_code == 1
        assert "-s" in result.output

    def test_blocked_z_in_tool_args(
        self, runner: CliRunner, matrix_file: Path, candidate_trees_file: Path,
    ) -> None:
        result = runner.invoke(cli, [
            "posttree", "topology",
            "--matrix", str(matrix_file),
            "--candidate-trees", str(candidate_trees_file),
            "--model-expr", "LG+F+R4",
            "--tool-args", "-z other.trees",
        ])
        assert result.exit_code == 1
        assert "-z" in result.output

    def test_tool_args_accepted(
        self, runner: CliRunner, matrix_file: Path, candidate_trees_file: Path,
        tmp_path: Path,
    ) -> None:
        out = tmp_path / "out"
        result = runner.invoke(cli, [
            "posttree", "topology",
            "--matrix", str(matrix_file),
            "--candidate-trees", str(candidate_trees_file),
            "--model-expr", "LG+F+R4",
            "--tool-args", "--prefix custom -T 30 -fs custom.sitefreq",
            "--output-dir", str(out),
            "--dry-run",
        ])
        assert result.exit_code == 0
        assert "iqtree3" in result.output


# ------------------------------------------------------------------
# Dry-run
# ------------------------------------------------------------------

class TestCLIDryRun:
    def test_dry_run_single_tree_file(
        self, runner: CliRunner, matrix_file: Path, candidate_trees_file: Path,
        tmp_path: Path,
    ) -> None:
        out = tmp_path / "out"
        result = runner.invoke(cli, [
            "posttree", "topology",
            "--matrix", str(matrix_file),
            "--candidate-trees", str(candidate_trees_file),
            "--model-expr", "LG+F+R4",
            "--output-dir", str(out),
            "--dry-run",
        ])
        assert result.exit_code == 0
        assert "Would run:" in result.output
        assert "iqtree3" in result.output

    def test_dry_run_multiple_tree_files(
        self, runner: CliRunner, matrix_file: Path, tmp_path: Path,
    ) -> None:
        (tmp_path / "h1.nwk").write_text("(a,b);\n")
        (tmp_path / "h2.nwk").write_text("(a,c);\n")
        out = tmp_path / "out"
        result = runner.invoke(cli, [
            "posttree", "topology",
            "--matrix", str(matrix_file),
            "--candidate-trees", str(tmp_path / "h1.nwk"),
            "--candidate-trees", str(tmp_path / "h2.nwk"),
            "--model-expr", "LG+F+R4",
            "--output-dir", str(out),
            "--dry-run",
        ])
        assert result.exit_code == 0
        assert "Would run:" in result.output


# ------------------------------------------------------------------
# Integration (real IQ-TREE)
# ------------------------------------------------------------------

class TestCLIIntegration:
    @pytest.mark.skipif(
        not shutil.which("iqtree3"),
        reason="iqtree3 not found in PATH",
    )
    def test_successful_run_writes_result_json(
        self, runner: CliRunner, tmp_path: Path,
    ) -> None:
        matrix = tmp_path / "matrix.fa"
        matrix.write_text(
            ">t1\nMKTLLLTLWVV\n>t2\nMKTLLLTLWVI\n>t3\nMKTLLLSLWVI\n>t4\nMKTLLLTLWVA\n"
        )
        (tmp_path / "trees").write_text(
            "(t1,t2,(t3,t4));\n(t1,t3,(t2,t4));\n"
        )
        out = tmp_path / "out"

        result = runner.invoke(cli, [
            "posttree", "topology",
            "--matrix", str(matrix),
            "--candidate-trees", str(tmp_path / "trees"),
            "--model-expr", "LG",
            "--replicates", "1000",
            "--output-dir", str(out),
            "--threads", "1",
        ])
        assert result.exit_code == 0
        assert (out / "result.json").exists()

        with open(out / "result.json") as fh:
            payload = json.load(fh)
        assert payload["status"] == "success"
        assert "iqtree3" in payload["tool_versions"]
        assert len(payload["data"]["tests"]) == 2
```

- [ ] **Step 2: Run CLI tests**

```bash
pytest tests/cli/test_posttree_topology.py -v --tb=short
```

Expected: help content tests PASS, validation tests PASS (missing required options → exit 1), integration PASS or SKIP.

- [ ] **Step 3: [skip commit]**

---

### Task 8: Write user-facing command documentation

**Files:**
- Create: `docs/commands/posttree-topology.md`

- [ ] **Step 1: Write documentation**

```markdown
# phyloai posttree topology

## Purpose

Performs IQ-TREE tree topology tests (AU / KH / SH / WKH / WSH / c-ELW)
comparing a set of candidate trees against a supermatrix alignment. This
command tests whether alternative topologies are significantly worse than
the best-scoring candidate — it does **not** infer new trees.

## Usage

```bash
# Homogeneous model
phyloai posttree topology \
  --matrix matrix.fa \
  --candidate-trees candidates.trees \
  --model-expr LG+F+R4

# PMSF model with guide tree
phyloai posttree topology \
  --matrix matrix.fa \
  --candidate-trees candidates.trees \
  --model-expr LG+C20+F+R4 \
  --guide-tree guide.nwk

# Previously optimized partition model
phyloai posttree topology \
  --matrix matrix.fa \
  --candidate-trees candidates.trees \
  --partitions matrix.best_model.nex

# Multiple individual tree files (merged by PhyloAI)
phyloai posttree topology \
  --matrix matrix.fa \
  --candidate-trees h1.nwk \
  --candidate-trees h2.nwk \
  --candidate-trees h3.nwk \
  --model-expr LG+F+R4

# Custom exchangeabilities + site frequencies via --tool-args
phyloai posttree topology \
  --matrix matrix.fa \
  --candidate-trees trees \
  --model-expr custom.exchangeabilities+R4 \
  --tool-args "-fs custom.sitefreq" -t 30

# Heterogeneous model
phyloai posttree topology \
  --matrix matrix.fa \
  --candidate-trees trees \
  --model-expr C20+F+R4
```

## Inputs

| Input | Description |
|-------|-------------|
| `--matrix` | Single supermatrix alignment (FASTA, PHYLIP, NEXUS, or CLUSTAL). Maps to IQ-TREE `-s`. |
| `--candidate-trees` | One tree-list file (one NEWICK tree per line) or multiple individual NEWICK files. Multiple files are merged in order by PhyloAI into `candidate.trees`. Maps to IQ-TREE `-z`. |
| `--input-format` | Optional PhyloAI-side matrix format hint (`auto` by default). Not passed to IQ-TREE. |

## Model Source

Provide exactly one model source. PhyloAI does **not** re-run ModelFinder — use
`phyloai tree ml iqtree` for model selection.

| Option | Description |
|--------|-------------|
| `--model-expr` | Complete IQ-TREE `-m` expression. Examples: `LG+F+R4`, `C20+F+R4`, `LG+C20+F+R4`, `custom.exchangeabilities+R4`. |
| `--partitions` | Previously optimized partition file (e.g., `.best_model.nex` from IQ-TREE). Maps to IQ-TREE `-p`. |
| `--tool-args "-m ..."` or `-p` | Model source provided through raw IQ-TREE flags. If a high-level `--model-expr` or `--partitions` is also set,
the same-named flag in `--tool-args` takes precedence (suppress-if-present); providing the opposite flag
(e.g. `--model-expr` with `-p` in `--tool-args`) is rejected as a cross-source conflict. |

`--guide-tree` is used with PMSF models (e.g., `LG+C20+F+R4`). Maps to IQ-TREE `-ft`.

## Default Tests

PhyloAI generates the standard topology-test flags:

```
-n 0 -zb <replicates> -zw -au
```

| Test | Description |
|------|-------------|
| bp-RELL | Bootstrap proportion (RELL) |
| KH | Kishino-Hasegawa test |
| SH | Shimodaira-Hasegawa test |
| WKH | Weighted KH test |
| WSH | Weighted SH test |
| c-ELW | Expected likelihood weight |
| AU | Approximately unbiased test |

KH, SH, WKH, WSH, and AU are **p-values**. Trees with p < 0.05 are rejected
by that test. bp-RELL and c-ELW are **weights**, not p-values. The AU test is
generally considered the most reliable.

## Advanced IQ-TREE Args

| Flag | Description |
|------|-------------|
| `--tool-args` | Additional IQ-TREE strategy parameters. **Blocked flags:** `-s` (matrix), `-z` (candidate trees). Shell I/O redirects (`<`, `>`, `|`) rejected. |
| `--iqtree-path` | Explicit path to `iqtree3` executable. |
| `--prefix` | IQ-TREE output prefix (default: matrix file stem). |

PhyloAI-built flags are suppressed when the same flag appears in `--tool-args`
(suppress-if-present). Overrideable flags: `-m`, `-p`, `-ft`, `-n`, `-zb`,
`-zw`, `-au`, `-T`, `--prefix`.

## Input Format and Sequence Type

`--input-format` (`auto|fasta|phylip-relaxed|nexus|clustal`) only affects
PhyloAI's own matrix preflight validation; it is **not** passed to IQ-TREE.
Explicit IQ-TREE `--seqtype` belongs in `--tool-args` when needed
(e.g., `--tool-args "--seqtype AA"`).

## Outputs

IQ-TREE native files:
- `<prefix>.iqtree` — full IQ-TREE report with topology test table
- `<prefix>.log` — IQ-TREE log
- `<prefix>.treels.trees` — IQ-TREE optimized candidate trees (suffix may vary)

PhyloAI files:
- `result.json` — structured result with parsed test table
- `candidate.trees` — merged tree file (only when multiple `--candidate-trees` were provided)

## Exit Codes

| Code | Meaning |
|------|---------|
| 0 | Success |
| 1 | User input error (missing files, invalid parameters, output conflict) |
| 2 | IQ-TREE execution failed |
| 3 | IQ-TREE executable not found |

## Notes

- This command is single-matrix only (no batch mode).
- `--replicates` defaults to 10000. Very large values can make RELL resampling slow.
- IQ-TREE native resume (via `.ckp.gz`) is supported with `--resume`.
```

- [ ] **Step 2: [skip commit]**

---

### Task 9: Run full test suite and make final commit

- [ ] **Step 1: Run all topology tests**

```bash
pytest tests/posttree/ tests/cli/test_posttree_topology.py -v --tb=short
```

Expected: all PASS or SKIP.

- [ ] **Step 2: Run full test suite for regressions**

```bash
pytest tests/ -x -q --tb=short
```

Expected: no regressions. All previously-passing tests continue to pass.

- [ ] **Step 3: Single final commit**

```bash
git status
git add phyloai/core/iqtree.py
git add phyloai/tree/ml_iqtree.py phyloai/tree/cf.py
git add phyloai/posttree/__init__.py phyloai/posttree/topology.py
git add phyloai/cli/commands/posttree.py phyloai/cli/main.py
git add tests/posttree/__init__.py tests/posttree/test_topology.py
git add tests/cli/test_posttree_topology.py
git add docs/commands/posttree-topology.md
git commit -m "feat: implement phyloai posttree topology

- Extract shared IQ-TREE helpers to phyloai/core/iqtree.py
- Refactor ml_iqtree.py and cf.py to import from shared module
- Implement run_topology() with library-layer input validation
- Add header-offset tolerant USER TREES table parser
- Add CLI command with 8-section help and manual validation
- Comprehensive unit and integration tests
- User-facing command documentation"
```
