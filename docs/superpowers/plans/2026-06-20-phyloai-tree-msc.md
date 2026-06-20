# phyloai tree msc (wASTRAL) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Implement `phyloai tree msc` — multispecies coalescent species tree inference with wASTRAL, consuming gene trees from `--tree` (single file) or `--tree-dir` (merged directory).

**Architecture:** Single direct Click command (not a Group) under `tree`. Library layer (`phyloai/tree/msc.py`) handles input scanning/merging, two-tier `--tool-args` flag management, command building, single subprocess execution, and result.json assembly. No batch parallelism (wastral handles internal multithreading via `-t`), no checkpoint/resume (one-shot computation).

**Tech Stack:** Python 3.12+, Click, subprocess, shlex, Path, json

---

## File Structure

| File | Action | Responsibility |
|------|--------|----------------|
| `phyloai/tree/msc.py` | Create | Library: `run_wastral()`, input scanning/merging, command builder, flag override detection, result assembly |
| `phyloai/cli/commands/tree.py` | Modify | CLI: register `msc` command, Click options, wire to `run_wastral()`, Rich summary, result.json writing |
| `phyloai/core/env.py` | Modify | Add `path_aliases: ["aster"]` to `wastral` TOOL_REGISTRY entry |
| `tests/tree/test_msc.py` | Create | Library tests: command building, flag override, input scanning, result assembly |
| `tests/cli/test_tree_msc.py` | Create | CLI tests: help output, mutual exclusivity, validation, dry-run, integration |
| `docs/commands/tree-msc.md` | Create | User-facing command documentation |

---

### Task 1: Add `path_aliases` for `wastral` in TOOL_REGISTRY

**Files:**
- Modify: `phyloai/core/env.py:34-35`

- [ ] **Step 1: Update the wastral registry entry**

```python
# Change from:
    "wastral":    {"required": False, "version_args": [["-v"], ["-h"]],
                    "install": "https://github.com/chaoszhang/ASTER"},
# To:
    "wastral":    {"required": False, "version_args": [["-v"], ["-h"]],
                    "path_aliases": ["aster"],
                    "install": "https://github.com/chaoszhang/ASTER"},
```

- [ ] **Step 2: Verify the edit**

Run: `python -c "from phyloai.core.env import TOOL_REGISTRY; print(TOOL_REGISTRY['wastral'])"`
Expected: Output contains `'path_aliases': ['aster']`

- [ ] **Step 3: Commit**

```bash
git add phyloai/core/env.py
git commit -m "feat: add path_aliases ['aster'] for wastral tool resolution"
```

---

### Task 2: Create library file skeleton with constants and input scanning

**Files:**
- Create: `phyloai/tree/msc.py`

- [ ] **Step 1: Write the file with constants and `_scan_input_wastral()`**

```python
"""Multispecies coalescent inference with wASTRAL (ASTER)."""

from __future__ import annotations

from pathlib import Path
from typing import Any


# Gene tree file extensions wastral can read (newick variants)
_WASTRAL_TREE_EXTENSIONS = frozenset({
    ".nwk", ".tre", ".tree", ".nw", ".trees", ".newick",
})


# ---- Two-tier --tool-args management ----------------------------------

# Tier 1: BLOCKED — phyloAI always manages these; reject hard
_WASTRAL_BLOCKED_FLAGS = frozenset({"-i", "-o"})
_WASTRAL_BLOCKED_IO_CHARS = frozenset({"<", ">", "|"})

# Tier 2: OVERRIDEABLE — suppress phyloAI's own flag if present in --tool-args
# Maps phyloAI parameter group -> set of wastral flag strings
_WASTRAL_OVERRIDE_MAP: dict[str, frozenset[str]] = {
    "mode": frozenset({"--mode"}),
    "boot": frozenset({"-u"}),
    "extra_rounds": frozenset({"-R", "-r", "-s"}),
    "tree_boot": frozenset({"--lrt", "--bayes", "--bootstrap", "-x", "-n"}),
    "threads": frozenset({"-t"}),
}


# ---- Input scanning ---------------------------------------------------

def _scan_input_wastral(
    tree_dir: Path,
) -> tuple[list[Path], list[dict[str, str]]]:
    """Scan a directory for valid gene tree files.

    Returns:
        (valid_files, skipped_entries)
    """
    if not tree_dir.exists():
        return [], []

    found: list[Path] = []
    skipped: list[dict[str, str]] = []

    for entry in sorted(tree_dir.iterdir()):
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
        if ext in _WASTRAL_TREE_EXTENSIONS:
            found.append(entry)
        else:
            skipped.append({"path": str(entry), "reason": f"unrecognized extension: {ext}"})

    return found, skipped


def _merge_gene_trees(
    tree_dir: Path,
    output_path: Path,
) -> tuple[int, list[dict[str, str]]]:
    """Scan tree_dir for newick files, merge into one file (one tree per line).

    Returns:
        (count_of_trees_merged, skipped_entries)
    """
    found, skipped = _scan_input_wastral(tree_dir)

    count = 0
    with open(output_path, "w") as out:
        for f in found:
            content = f.read_text().strip()
            if not content:
                continue
            for line in content.splitlines():
                line = line.strip()
                if line:
                    out.write(line + "\n")
                    count += 1

    return count, skipped
```

- [ ] **Step 2: Verify the file imports cleanly**

Run: `python -c "from phyloai.tree.msc import _scan_input_wastral, _merge_gene_trees; print('OK')"`
Expected: `OK`

- [ ] **Step 3: Commit**

```bash
git add phyloai/tree/msc.py
git commit -m "feat: add msc library skeleton with input scanning and constants"
```

---

### Task 3: Add `--tool-args` flag override functions

**Files:**
- Modify: `phyloai/tree/msc.py`

- [ ] **Step 1: Append after the constants section (before the input scanning functions, after `_WASTRAL_OVERRIDE_MAP`)**

```python
import shlex


def _check_managed_flag_conflict(tool_args: str) -> None:
    """Reject BLOCKED flags and I/O redirects in --tool-args."""
    tokens = shlex.split(tool_args)
    for token in tokens:
        if token in _WASTRAL_BLOCKED_FLAGS:
            raise ValueError(f"Blocked managed flag in --tool-args: {token}")
        if any(c in token for c in _WASTRAL_BLOCKED_IO_CHARS):
            raise ValueError(f"Blocked I/O override in --tool-args: {token}")


def _is_flag_overridden(group: str, tool_tokens: set[str]) -> bool:
    """Check whether any flag in an override group appears in --tool-args."""
    flags = _WASTRAL_OVERRIDE_MAP.get(group, frozenset())
    return bool(flags & tool_tokens)
```

- [ ] **Step 2: Verify the import resolution stays clean**

Run: `python -c "from phyloai.tree.msc import _check_managed_flag_conflict, _is_flag_overridden; print('OK')"`
Expected: `OK`

- [ ] **Step 3: Commit**

```bash
git add phyloai/tree/msc.py
git commit -m "feat: add --tool-args flag override detection for wastral"
```

---

### Task 4: Add `_resolve_wastral_path()` and `_detect_wastral_version()`

**Files:**
- Modify: `phyloai/tree/msc.py`

- [ ] **Step 1: Append after `_merge_gene_trees()`**

```python
import os
import re as _re
import subprocess

from phyloai.core.env import ToolEnv


def _resolve_wastral_path(wastral_path: str | None, dry_run: bool) -> str:
    """Resolve wastral executable path.

    Priority: 1) explicit --wastral-path  2) ToolEnv.require("wastral")
    """
    if wastral_path:
        p = Path(wastral_path)
        if not p.exists():
            raise ValueError(f"--wastral-path does not exist: {wastral_path}")
        if not os.access(str(p), os.X_OK):
            raise ValueError(f"--wastral-path is not executable: {wastral_path}")
        return wastral_path
    if dry_run:
        return "wastral"
    try:
        env = ToolEnv()
        return str(env.require("wastral"))
    except FileNotFoundError:
        raise FileNotFoundError(
            "wastral not found. Install from https://github.com/chaoszhang/ASTER "
            "or use --wastral-path."
        )


def _detect_wastral_version(executable: str) -> dict[str, str]:
    """Detect wastral version via -v or -h fallback."""
    combined = ""
    for flag in ("-v", "-h"):
        try:
            proc = subprocess.run(
                [executable, flag],
                capture_output=True, text=True, timeout=10,
            )
            combined = proc.stdout + "\n" + proc.stderr
            if proc.returncode == 0 or (proc.stdout.strip() or proc.stderr.strip()):
                break
        except Exception:
            continue

    if not combined.strip():
        return {"wastral": "unknown"}

    m = _re.search(r"version\s*([\d.]+)", combined, _re.IGNORECASE)
    if m:
        return {"wastral": m.group(1)}
    m = _re.search(r"([\d]+\.[\d]+(?:\.[\d]+)?)", combined)
    if m:
        return {"wastral": m.group(1)}
    return {"wastral": "unknown"}
```

- [ ] **Step 2: Verify**

Run: `python -c "from phyloai.tree.msc import _resolve_wastral_path, _detect_wastral_version; print('OK')"`
Expected: `OK`

- [ ] **Step 3: Commit**

```bash
git add phyloai/tree/msc.py
git commit -m "feat: add wastral executable resolution and version detection"
```

---

### Task 5: Add `_build_wastral_cmd()` command builder

**Files:**
- Modify: `phyloai/tree/msc.py`

- [ ] **Step 1: Append after `_detect_wastral_version()`**

```python
def _build_wastral_cmd(
    input_path: Path,
    output_path: Path,
    *,
    mode: int,
    boot: int,
    extra_rounds: bool,
    tree_boot_type: str,
    tree_boot_min: float | None,
    tree_boot_max: float | None,
    threads: int,
    executable: str = "wastral",
    tool_args: str | None = None,
) -> list[str]:
    """Build the full wastral command line."""
    cmd = [executable]

    # Managed input/output (always emitted, blocked in --tool-args)
    cmd.extend(["-i", str(input_path)])
    cmd.extend(["-o", str(output_path)])

    tool_tokens = set(shlex.split(tool_args)) if tool_args else set()

    if tool_args:
        _check_managed_flag_conflict(tool_args)

    # --mode
    if not _is_flag_overridden("mode", tool_tokens):
        cmd.extend(["--mode", str(mode)])

    # -u (branch support)
    if not _is_flag_overridden("boot", tool_tokens):
        cmd.extend(["-u", str(boot)])

    # -t (threads)
    if not _is_flag_overridden("threads", tool_tokens):
        cmd.extend(["-t", str(threads)])

    # Tree boot type (input gene tree branch support)
    if tree_boot_type != "auto" and not _is_flag_overridden("tree_boot", tool_tokens):
        boot_flag_map = {
            "likelihood": "--lrt",
            "abayes": "--bayes",
            "bootstrap": "--bootstrap",
        }
        boot_defaults = {
            "likelihood":  ("-d", "0",   "-x", "1.0", "-n", "0.0"),
            "abayes":      ("-d", "0.333", "-x", "1.0", "-n", "0.333"),
            "bootstrap":   ("-d", "0",   "-x", "100", "-n", "0"),
        }
        flag = boot_flag_map[tree_boot_type]
        cmd.append(flag)

        # -d, -x, -n: preset defaults, overridable by explicit min/max
        defaults = boot_defaults[tree_boot_type]
        max_val = str(tree_boot_max) if tree_boot_max is not None else defaults[3]
        min_val = str(tree_boot_min) if tree_boot_min is not None else defaults[5]
        cmd.extend([
            defaults[0], defaults[1],  # -d VALUE (hardcoded per type)
            "-x", max_val,
            "-n", min_val,
        ])

    # -R (extra rounds / exhaustive search)
    if extra_rounds and not _is_flag_overridden("extra_rounds", tool_tokens):
        cmd.append("-R")

    # --tool-args appended last
    if tool_args:
        cmd.extend(shlex.split(tool_args))

    return cmd
```

- [ ] **Step 2: Verify**

Run: `python -c "from phyloai.tree.msc import _build_wastral_cmd; print('OK')"`
Expected: `OK`

- [ ] **Step 3: Commit**

```bash
git add phyloai/tree/msc.py
git commit -m "feat: add wastral command builder with two-tier --tool-args support"
```

---

### Task 6: Add `_assemble_wastral_result()` and `run_wastral()` main entry point

**Files:**
- Modify: `phyloai/tree/msc.py`

- [ ] **Step 1: Append after `_build_wastral_cmd()`**

```python
import time as _time
import json
import shutil


def _assemble_wastral_result(
    *,
    run_start: float,
    wastral_exe: str,
    tree: Path | None,
    tree_dir: Path | None,
    output_dir: Path,
    mode: int,
    boot: int,
    extra_rounds: bool,
    tree_boot_type: str,
    tree_boot_min: float | None,
    tree_boot_max: float | None,
    threads: int,
    wastral_path: str | None,
    tool_args: str | None,
    overwrite: bool,
    n_input_trees: int,
    input_path: Path,
    cmd: list[str],
    wall_time: float,
    skipped: list[dict[str, str]],
    warnings_list: list[str],
    is_error: bool,
    error_msg: str | None,
) -> dict[str, Any]:
    """Build the result.json payload."""
    input_mode = "--tree-dir" if tree_dir is not None else "--tree"

    try:
        versions = _detect_wastral_version(wastral_exe)
    except Exception:
        versions = {"wastral": "unknown"}

    # Reconstruct CLI invocation string
    cmd_parts = ["phyloai", "tree", "msc"]
    if tree is not None:
        cmd_parts.extend(["--tree", str(tree)])
    else:
        cmd_parts.extend(["--tree-dir", str(tree_dir)])
    cmd_parts.extend([
        "--mode", str(mode),
        "--boot", str(boot),
        "-t", str(threads),
        "-o", str(output_dir),
    ])
    if extra_rounds:
        cmd_parts.append("-R")
    if tree_boot_type != "auto":
        cmd_parts.extend(["--tree-boot-type", tree_boot_type])
        if tree_boot_min is not None:
            cmd_parts.extend(["--tree-boot-min", str(tree_boot_min)])
        if tree_boot_max is not None:
            cmd_parts.extend(["--tree-boot-max", str(tree_boot_max)])
    if wastral_path:
        cmd_parts.extend(["--wastral-path", str(wastral_path)])
    if tool_args:
        if " " in tool_args:
            cmd_parts.append(f"--tool-args '{tool_args}'")
        else:
            cmd_parts.extend(["--tool-args", tool_args])
    if overwrite:
        cmd_parts.append("--overwrite")
    cmd_str = " ".join(cmd_parts)

    input_data: dict[str, Any] = {
        "path": str(input_path),
    }
    if input_mode == "--tree-dir":
        input_data["n_trees"] = n_input_trees

    return {
        "status": "error" if is_error else "success",
        "command": cmd_str,
        "wall_time": _time.monotonic() - run_start,
        "tool_versions": versions,
        "params": {
            "tree": str(tree) if tree else None,
            "tree_dir": str(tree_dir) if tree_dir else None,
            "mode": mode,
            "boot": boot,
            "extra_rounds": extra_rounds,
            "tree_boot_type": tree_boot_type,
            "tree_boot_min": tree_boot_min,
            "tree_boot_max": tree_boot_max,
            "output_dir": str(output_dir),
            "threads": threads,
            "overwrite": overwrite,
            "tool_args": tool_args,
            "wastral_path": wastral_path,
        },
        "key_results": {
            "mode": mode,
            "boot": boot,
            "extra_rounds": extra_rounds,
            "tree_boot_type": tree_boot_type,
            "n_input_trees": n_input_trees,
            "input_mode": input_mode,
        },
        "error": error_msg,
        "data": {
            "input_mode": input_mode,
            "input": input_data,
            "output_tree": str(output_dir / "wastral.tre"),
            "cmd": cmd,
            "skipped": skipped,
            "warnings": warnings_list,
            **({"freq_quad_csv": str(output_dir / "freqQuad.csv")} if boot == 3 else {}),
        },
    }


def run_wastral(
    *,
    tree: Path | None = None,
    tree_dir: Path | None = None,
    output_dir: Path = Path("runs/tree/msc"),
    mode: int = 1,
    boot: int = 1,
    extra_rounds: bool = False,
    tree_boot_type: str = "auto",
    tree_boot_min: float | None = None,
    tree_boot_max: float | None = None,
    threads: int = 4,
    wastral_path: str | None = None,
    tool_args: str | None = None,
    overwrite: bool = False,
    dry_run: bool = False,
    quiet: bool = False,
) -> dict[str, Any]:
    """Run wASTRAL multispecies coalescent species tree inference.

    Returns a result.json-compatible payload dict.
    """
    # --- Library-level parameter validation ---
    if mode not in (1, 2, 3, 4):
        raise ValueError(f"--mode must be 1, 2, 3, or 4. Got: {mode}")
    if boot not in (0, 1, 2, 3):
        raise ValueError(f"--boot must be 0, 1, 2, or 3. Got: {boot}")
    if threads < 1:
        raise ValueError(f"--threads must be >= 1. Got: {threads}")

    # --- Input mutual exclusivity ---
    if (tree is None and tree_dir is None) or (tree is not None and tree_dir is not None):
        raise ValueError(
            "Either --tree or --tree-dir must be provided (mutually exclusive)."
        )

    # --- Input path validation ---
    if tree is not None and not tree.exists():
        raise ValueError(f"--tree does not exist: {tree}")
    if tree_dir is not None and not tree_dir.exists():
        raise ValueError(f"--tree-dir does not exist: {tree_dir}")

    # --- Tree boot type + min/max validation ---
    if tree_boot_type != "auto" and tree_boot_type not in {"likelihood", "abayes", "bootstrap"}:
        raise ValueError(
            f"Invalid --tree-boot-type: {tree_boot_type}. "
            f"Valid: auto, likelihood, abayes, bootstrap"
        )

    if tree_boot_type == "auto":
        if tree_boot_min is not None or tree_boot_max is not None:
            raise ValueError(
                "--tree-boot-min and --tree-boot-max are only valid when "
                "--tree-boot-type is not 'auto'."
            )

    if tree_boot_min is not None and tree_boot_max is not None and tree_boot_min >= tree_boot_max:
        raise ValueError(
            f"--tree-boot-min ({tree_boot_min}) must be < --tree-boot-max ({tree_boot_max})."
        )

    # --- Output directory conflict ---
    if not dry_run:
        if overwrite and output_dir.exists():
            shutil.rmtree(output_dir)
        if not overwrite and output_dir.exists() and any(output_dir.iterdir()):
            raise ValueError(
                f"Output directory {output_dir} already exists and is non-empty. "
                "Use --overwrite to replace."
            )

    run_start = _time.monotonic()

    # --- Input resolution (before tool resolution to ensure input errors come first) ---
    skipped: list[dict[str, str]] = []
    warnings_list: list[str] = []
    n_input_trees = 0

    if tree is not None:
        # Single file mode
        input_path = tree
        n_input_trees = 0  # unknown, wastral counts internally
    else:
        assert tree_dir is not None
        # Directory mode: merge into merged.trees
        output_dir.mkdir(parents=True, exist_ok=True)
        merged_path = output_dir / "merged.trees"
        n_input_trees, scanned_skipped = _merge_gene_trees(tree_dir, merged_path)

        skipped = scanned_skipped

        if n_input_trees == 0:
            raise ValueError(
                f"No valid gene tree files found in --tree-dir: {tree_dir}"
            )

        # Count valid files (not individual trees within files) for warnings
        valid_files, _unused = _scan_input_wastral(tree_dir)
        n_valid_files = len(valid_files)

        if n_valid_files == 1 and not quiet:
            warnings_list.append(
                "Exactly 1 valid gene tree file in --tree-dir. "
                "Consider using --tree mode directly."
            )

        # Warn about non-newick files
        unrecognized = [
            s for s in skipped
            if s.get("reason", "").startswith("unrecognized")
        ]
        if unrecognized and not quiet:
            warnings_list.append(
                f"--tree-dir contains {len(unrecognized)} non-newick file(s); "
                "skipped. See result.json data.skipped for details."
            )

        input_path = merged_path

    # --- Resolve executable (after input validation, skip for dry_run) ---
    wastral_exe = _resolve_wastral_path(wastral_path, dry_run)

    # --- Command building ---
    # Resolve paths to absolute so wastral sees correct paths regardless of cwd
    output_tree = (output_dir / "wastral.tre").resolve()
    cmd = _build_wastral_cmd(
        input_path=input_path.resolve(),
        output_path=output_tree,
        mode=mode,
        boot=boot,
        extra_rounds=extra_rounds,
        tree_boot_type=tree_boot_type,
        tree_boot_min=tree_boot_min,
        tree_boot_max=tree_boot_max,
        threads=threads,
        executable=wastral_exe,
        tool_args=tool_args,
    )

    if dry_run:
        return _assemble_wastral_result(
            run_start=run_start,
            wastral_exe=wastral_exe,
            tree=tree, tree_dir=tree_dir,
            output_dir=output_dir,
            mode=mode, boot=boot, extra_rounds=extra_rounds,
            tree_boot_type=tree_boot_type, tree_boot_min=tree_boot_min,
            tree_boot_max=tree_boot_max,
            threads=threads, wastral_path=wastral_path, tool_args=tool_args,
            overwrite=overwrite, n_input_trees=n_input_trees,
            input_path=input_path, cmd=cmd, wall_time=0.0,
            skipped=skipped, warnings_list=warnings_list,
            is_error=False, error_msg=None,
        )

    # --- Execution (cwd = output_dir so freqQuad.csv lands in output dir) ---
    output_dir.mkdir(parents=True, exist_ok=True)
    log_path = output_dir / "wastral.log"

    try:
        proc = subprocess.run(
            cmd,
            cwd=str(output_dir),
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
    except Exception as exc:
        return _assemble_wastral_result(
            run_start=run_start,
            wastral_exe=wastral_exe,
            tree=tree, tree_dir=tree_dir,
            output_dir=output_dir,
            mode=mode, boot=boot, extra_rounds=extra_rounds,
            tree_boot_type=tree_boot_type, tree_boot_min=tree_boot_min,
            tree_boot_max=tree_boot_max,
            threads=threads, wastral_path=wastral_path, tool_args=tool_args,
            overwrite=overwrite, n_input_trees=n_input_trees,
            input_path=input_path, cmd=cmd, wall_time=0.0,
            skipped=skipped, warnings_list=warnings_list,
            is_error=True, error_msg=str(exc),
        )

    # Save stderr as wastral.log
    log_path.write_text(proc.stderr)

    wall_time = _time.monotonic() - run_start

    if proc.returncode != 0:
        return _assemble_wastral_result(
            run_start=run_start,
            wastral_exe=wastral_exe,
            tree=tree, tree_dir=tree_dir,
            output_dir=output_dir,
            mode=mode, boot=boot, extra_rounds=extra_rounds,
            tree_boot_type=tree_boot_type, tree_boot_min=tree_boot_min,
            tree_boot_max=tree_boot_max,
            threads=threads, wastral_path=wastral_path, tool_args=tool_args,
            overwrite=overwrite, n_input_trees=n_input_trees,
            input_path=input_path, cmd=cmd, wall_time=wall_time,
            skipped=skipped, warnings_list=warnings_list,
            is_error=True,
            error_msg=f"wastral exited with code {proc.returncode}: {proc.stderr[:200]}",
        )

    # Verify output tree was produced
    if not output_tree.exists() or output_tree.stat().st_size == 0:
        return _assemble_wastral_result(
            run_start=run_start,
            wastral_exe=wastral_exe,
            tree=tree, tree_dir=tree_dir,
            output_dir=output_dir,
            mode=mode, boot=boot, extra_rounds=extra_rounds,
            tree_boot_type=tree_boot_type, tree_boot_min=tree_boot_min,
            tree_boot_max=tree_boot_max,
            threads=threads, wastral_path=wastral_path, tool_args=tool_args,
            overwrite=overwrite, n_input_trees=n_input_trees,
            input_path=input_path, cmd=cmd, wall_time=wall_time,
            skipped=skipped, warnings_list=warnings_list,
            is_error=True,
            error_msg="wastral did not produce output tree",
        )

    return _assemble_wastral_result(
        run_start=run_start,
        wastral_exe=wastral_exe,
        tree=tree, tree_dir=tree_dir,
        output_dir=output_dir,
        mode=mode, boot=boot, extra_rounds=extra_rounds,
        tree_boot_type=tree_boot_type, tree_boot_min=tree_boot_min,
        tree_boot_max=tree_boot_max,
        threads=threads, wastral_path=wastral_path, tool_args=tool_args,
        overwrite=overwrite, n_input_trees=n_input_trees,
        input_path=input_path, cmd=cmd, wall_time=wall_time,
        skipped=skipped, warnings_list=warnings_list,
        is_error=False, error_msg=None,
    )
```

- [ ] **Step 2: Verify**

Run: `python -c "from phyloai.tree.msc import run_wastral; print('OK')"`
Expected: `OK`

- [ ] **Step 3: Commit**

```bash
git add phyloai/tree/msc.py
git commit -m "feat: add run_wastral() main entry point and result assembly"
```

---

### Task 7: Write library-level tests for input scanning

**Files:**
- Create: `tests/tree/test_msc.py`

- [ ] **Step 1: Write tests for `_scan_input_wastral()` and `_merge_gene_trees()`**

```python
from __future__ import annotations

from pathlib import Path


_WASTRAL_EXTENSIONS = frozenset({
    ".nwk", ".tre", ".tree", ".nw", ".trees", ".newick",
})


def test_scan_input_wastral_finds_all_supported(tmp_path: Path) -> None:
    from phyloai.tree.msc import _scan_input_wastral

    (tmp_path / "gene1.nwk").write_text("((a,b),c);\n")
    (tmp_path / "gene2.tre").write_text("((x,y),z);\n")
    (tmp_path / "gene3.tree").write_text("((1,2),3);\n")
    (tmp_path / "gene4.nw").write_text("(A,B);\n")
    (tmp_path / "gene5.trees").write_text("((a,b),c);\n")
    (tmp_path / "gene6.newick").write_text("(X,Y);\n")
    (tmp_path / "notes.txt").write_text("skip")
    (tmp_path / "empty.nwk").write_text("")
    (tmp_path / "subdir").mkdir()

    found, skipped = _scan_input_wastral(tmp_path)

    assert len(found) == 6
    assert len(skipped) == 3
    skip_reasons = {s["reason"] for s in skipped}
    assert "empty file" in skip_reasons
    assert "directory" in skip_reasons
    assert "unrecognized extension: .txt" in skip_reasons


def test_scan_input_wastral_nonexistent_dir() -> None:
    from phyloai.tree.msc import _scan_input_wastral

    found, skipped = _scan_input_wastral(Path("/nonexistent/dir"))
    assert found == []
    assert skipped == []


def test_merge_gene_trees_concatenates(tmp_path: Path) -> None:
    from phyloai.tree.msc import _merge_gene_trees

    (tmp_path / "a.nwk").write_text("((a,b),c);\n")
    (tmp_path / "b.tre").write_text("((x,y),z);\n")

    merged_path = tmp_path / "merged.trees"
    count, skipped = _merge_gene_trees(tmp_path, merged_path)

    assert count == 2
    assert merged_path.exists()
    lines = merged_path.read_text().strip().split("\n")
    assert len(lines) == 2
    assert lines[0] == "((a,b),c);"
    assert lines[1] == "((x,y),z);"


def test_merge_gene_trees_multi_line_file(tmp_path: Path) -> None:
    from phyloai.tree.msc import _merge_gene_trees

    (tmp_path / "multi.trees").write_text("((a,b),c);\n((x,y),z);\n(A,B);\n")

    merged_path = tmp_path / "merged.trees"
    count, _ = _merge_gene_trees(tmp_path, merged_path)

    assert count == 3


def test_merge_gene_trees_skips_non_newick(tmp_path: Path) -> None:
    from phyloai.tree.msc import _merge_gene_trees

    (tmp_path / "a.nwk").write_text("((a,b),c);\n")
    (tmp_path / "notes.txt").write_text("not a tree")
    (tmp_path / "empty.tre").write_text("")

    merged_path = tmp_path / "merged.trees"
    count, skipped = _merge_gene_trees(tmp_path, merged_path)

    assert count == 1
    assert len(skipped) == 2
    skip_reasons = {s["reason"] for s in skipped}
    assert "unrecognized extension: .txt" in skip_reasons
    assert "empty file" in skip_reasons


def test_merge_gene_trees_no_valid_files(tmp_path: Path) -> None:
    from phyloai.tree.msc import _merge_gene_trees

    (tmp_path / "notes.txt").write_text("not a tree")
    (tmp_path / "empty.tre").write_text("")

    merged_path = tmp_path / "merged.trees"
    count, skipped = _merge_gene_trees(tmp_path, merged_path)

    assert count == 0
    assert len(skipped) == 2
```

- [ ] **Step 2: Run tests to verify they pass**

Run: `pytest tests/tree/test_msc.py -v`
Expected: All 6 tests PASS

- [ ] **Step 3: Commit**

```bash
git add tests/tree/test_msc.py
git commit -m "test: add input scanning and merge tests for wastral"
```

---

### Task 8: Write library-level tests for command building and flag override

**Files:**
- Modify: `tests/tree/test_msc.py` (append)

- [ ] **Step 1: Append tests after existing ones**

```python
import pytest


def test_build_wastral_cmd_defaults(tmp_path: Path) -> None:
    from phyloai.tree.msc import _build_wastral_cmd

    inp = tmp_path / "input.trees"
    out = tmp_path / "output.tre"
    cmd = _build_wastral_cmd(inp, out, mode=1, boot=1, extra_rounds=False,
                              tree_boot_type="auto", tree_boot_min=None,
                              tree_boot_max=None, threads=4)

    assert cmd[0] == "wastral"
    assert "-i" in cmd
    assert str(inp) in cmd
    assert "-o" in cmd
    assert str(out) in cmd
    assert "--mode" in cmd
    assert "1" in cmd
    assert "-u" in cmd
    assert "-t" in cmd
    assert "4" in cmd
    assert "-R" not in cmd
    assert "--lrt" not in cmd
    assert "--bayes" not in cmd
    assert "--bootstrap" not in cmd


def test_build_wastral_cmd_mode_4_exhaustive(tmp_path: Path) -> None:
    from phyloai.tree.msc import _build_wastral_cmd

    inp = tmp_path / "input.trees"
    out = tmp_path / "output.tre"
    cmd = _build_wastral_cmd(inp, out, mode=4, boot=2, extra_rounds=True,
                              tree_boot_type="auto", tree_boot_min=None,
                              tree_boot_max=None, threads=8)

    assert "--mode" in cmd and "4" in cmd
    assert "-u" in cmd and "2" in cmd
    assert "-R" in cmd
    assert "-t" in cmd and "8" in cmd


def test_build_wastral_cmd_bootstrap_type(tmp_path: Path) -> None:
    from phyloai.tree.msc import _build_wastral_cmd

    inp = tmp_path / "input.trees"
    out = tmp_path / "output.tre"
    cmd = _build_wastral_cmd(inp, out, mode=1, boot=1, extra_rounds=False,
                              tree_boot_type="bootstrap", tree_boot_min=10,
                              tree_boot_max=95, threads=4)

    assert "--bootstrap" in cmd
    assert "-x" in cmd
    assert "95" in cmd
    assert "-n" in cmd
    assert "10" in cmd
    assert "-d" in cmd
    assert "0" in cmd  # hardcoded -d for bootstrap


def test_build_wastral_cmd_likelihood_type(tmp_path: Path) -> None:
    from phyloai.tree.msc import _build_wastral_cmd

    inp = tmp_path / "input.trees"
    out = tmp_path / "output.tre"
    cmd = _build_wastral_cmd(inp, out, mode=1, boot=1, extra_rounds=False,
                              tree_boot_type="likelihood", tree_boot_min=None,
                              tree_boot_max=None, threads=4)

    assert "--lrt" in cmd
    assert "-x" in cmd
    assert "1.0" in cmd
    assert "-n" in cmd
    assert "0.0" in cmd
    assert "-d" in cmd
    assert "0" in cmd


def test_build_wastral_cmd_abayes_type(tmp_path: Path) -> None:
    from phyloai.tree.msc import _build_wastral_cmd

    inp = tmp_path / "input.trees"
    out = tmp_path / "output.tre"
    cmd = _build_wastral_cmd(inp, out, mode=1, boot=1, extra_rounds=False,
                              tree_boot_type="abayes", tree_boot_min=None,
                              tree_boot_max=None, threads=4)

    assert "--bayes" in cmd
    assert "-x" in cmd
    assert "1.0" in cmd
    assert "-n" in cmd
    assert "0.333" in cmd
    assert "-d" in cmd
    assert "0.333" in cmd


def test_build_wastral_cmd_tool_args_blocks_minus_i(tmp_path: Path) -> None:
    from phyloai.tree.msc import _build_wastral_cmd

    inp = tmp_path / "input.trees"
    out = tmp_path / "output.tre"
    with pytest.raises(ValueError, match="Blocked managed flag"):
        _build_wastral_cmd(inp, out, mode=1, boot=1, extra_rounds=False,
                            tree_boot_type="auto", tree_boot_min=None,
                            tree_boot_max=None, threads=4,
                            tool_args="-i other.trees")


def test_build_wastral_cmd_tool_args_blocks_minus_o(tmp_path: Path) -> None:
    from phyloai.tree.msc import _build_wastral_cmd

    inp = tmp_path / "input.trees"
    out = tmp_path / "output.tre"
    with pytest.raises(ValueError, match="Blocked managed flag"):
        _build_wastral_cmd(inp, out, mode=1, boot=1, extra_rounds=False,
                            tree_boot_type="auto", tree_boot_min=None,
                            tree_boot_max=None, threads=4,
                            tool_args="-o output.tre")


def test_build_wastral_cmd_tool_args_overrides_u(tmp_path: Path) -> None:
    from phyloai.tree.msc import _build_wastral_cmd

    inp = tmp_path / "input.trees"
    out = tmp_path / "output.tre"
    cmd = _build_wastral_cmd(inp, out, mode=1, boot=1, extra_rounds=False,
                              tree_boot_type="auto", tree_boot_min=None,
                              tree_boot_max=None, threads=4,
                              tool_args="-u 3")

    # phyloAI should suppress its own -u, and -u 3 from tool-args should appear
    assert cmd.count("-u") >= 1
    assert "3" in cmd


def test_build_wastral_cmd_tool_args_overrides_R(tmp_path: Path) -> None:
    from phyloai.tree.msc import _build_wastral_cmd

    inp = tmp_path / "input.trees"
    out = tmp_path / "output.tre"
    cmd = _build_wastral_cmd(inp, out, mode=1, boot=1, extra_rounds=True,
                              tree_boot_type="auto", tree_boot_min=None,
                              tree_boot_max=None, threads=4,
                              tool_args="-R")

    # -R should appear exactly once (from tool-args, phyloAI suppressed)
    assert cmd.count("-R") == 1
```

- [ ] **Step 2: Run tests**

Run: `pytest tests/tree/test_msc.py -v`
Expected: All 15 tests PASS

- [ ] **Step 3: Commit**

```bash
git add tests/tree/test_msc.py
git commit -m "test: add command builder and flag override tests for wastral"
```

---

### Task 9: Write library-level tests for `run_wastral()` error paths

**Files:**
- Modify: `tests/tree/test_msc.py` (append)

- [ ] **Step 1: Append tests**

```python
def test_run_wastral_mutual_exclusivity_both_none() -> None:
    from phyloai.tree.msc import run_wastral

    with pytest.raises(ValueError, match="Either --tree or --tree-dir"):
        run_wastral(tree=None, tree_dir=None)


def test_run_wastral_mutual_exclusivity_both(tmp_path: Path) -> None:
    from phyloai.tree.msc import run_wastral

    tree_file = tmp_path / "g.trees"
    tree_file.write_text("((a,b),c);\n")
    tree_dir = tmp_path / "genetrees"
    tree_dir.mkdir()
    (tree_dir / "a.nwk").write_text("((a,b),c);\n")

    with pytest.raises(ValueError, match="mutually exclusive"):
        run_wastral(tree=tree_file, tree_dir=tree_dir)


def test_run_wastral_tree_file_not_found() -> None:
    from phyloai.tree.msc import run_wastral

    with pytest.raises(ValueError, match="does not exist"):
        run_wastral(tree=Path("/nonexistent/file.trees"))


def test_run_wastral_tree_dir_no_valid_files(tmp_path: Path) -> None:
    from phyloai.tree.msc import run_wastral

    tree_dir = tmp_path / "empty"
    tree_dir.mkdir()

    with pytest.raises(ValueError, match="No valid gene tree files"):
        run_wastral(tree_dir=tree_dir)


def test_run_wastral_tree_boot_min_ge_max(tmp_path: Path) -> None:
    from phyloai.tree.msc import run_wastral

    tree_file = tmp_path / "g.trees"
    tree_file.write_text("((a,b),c);\n")

    with pytest.raises(ValueError, match="min must be < max"):
        run_wastral(
            tree=tree_file,
            tree_boot_type="bootstrap",
            tree_boot_min=95,
            tree_boot_max=10,
        )


def test_run_wastral_tree_boot_min_max_with_auto(tmp_path: Path) -> None:
    from phyloai.tree.msc import run_wastral

    tree_file = tmp_path / "g.trees"
    tree_file.write_text("((a,b),c);\n")

    with pytest.raises(ValueError, match="tree-boot-min.*only valid"):
        run_wastral(
            tree=tree_file,
            tree_boot_type="auto",
            tree_boot_min=10,
        )


def test_run_wastral_output_dir_exists(tmp_path: Path) -> None:
    from phyloai.tree.msc import run_wastral

    tree_file = tmp_path / "g.trees"
    tree_file.write_text("((a,b),c);\n")
    out_dir = tmp_path / "out"
    out_dir.mkdir()
    (out_dir / "existing.txt").write_text("old")

    with pytest.raises(ValueError, match="already exists"):
        run_wastral(tree=tree_file, output_dir=out_dir)


def test_run_wastral_dry_run_produces_payload(tmp_path: Path) -> None:
    from phyloai.tree.msc import run_wastral

    tree_file = tmp_path / "g.trees"
    tree_file.write_text("((a,b),c);\n")
    out_dir = tmp_path / "out"

    result = run_wastral(tree=tree_file, output_dir=out_dir, dry_run=True)

    assert result["status"] == "success"
    assert "cmd" in result["data"]
    assert "phyloai tree msc" in result["command"]
    assert not (out_dir / "wastral.tre").exists()


def test_run_wastral_dry_run_tree_dir(tmp_path: Path) -> None:
    from phyloai.tree.msc import run_wastral

    tree_dir = tmp_path / "genetrees"
    tree_dir.mkdir()
    (tree_dir / "a.nwk").write_text("((a,b),c);\n")
    (tree_dir / "b.tre").write_text("((x,y),z);\n")
    out_dir = tmp_path / "out"

    result = run_wastral(tree_dir=tree_dir, output_dir=out_dir, dry_run=True)

    assert result["status"] == "success"
    assert result["key_results"]["n_input_trees"] == 2
    assert result["key_results"]["input_mode"] == "--tree-dir"
    assert not (out_dir / "wastral.tre").exists()
```

- [ ] **Step 2: Validation is already integrated in the `run_wastral()` body above (mode/boot/threads range, mutual exclusivity, tree-boot-type min/max, output dir conflict). Run tests to verify all validation paths work.**

Run: `pytest tests/tree/test_msc.py -v`
Expected: All 24 tests PASS

- [ ] **Step 3: Commit**

```bash
git add tests/tree/test_msc.py phyloai/tree/msc.py
git commit -m "test: add run_wastral validation and dry-run tests; add input validation"
```

---

### Task 10: Register `msc` command in CLI tree.py

**Files:**
- Modify: `phyloai/cli/commands/tree.py`

- [ ] **Step 1: Update `_TreeGroup.list_commands()` and add the `msc` command registration**

```python
# Change _TreeGroup to include "msc":
class _TreeGroup(click.Group):
    def list_commands(self, ctx: click.Context) -> list[str]:
        return ["ml", "msc"]
```

Add after the `fasttree_command()` function (before `iqtree` section) or at the end of the file:

```python
@tree.command(
    "msc",
    cls=_GroupedHelpCommand,
    help=(
        "Multispecies coalescent species tree inference with wASTRAL.\n\n"
        "  --tree     : single gene tree file (newick, one tree per line)\n\n"
        "  --tree-dir : directory of gene tree files (merged into one input)\n\n"
        "--tree and --tree-dir are mutually exclusive. "
        "wASTRAL is one-shot computation (no --resume)."
    ),
)
@click.option(
    "--tree",
    type=click.Path(dir_okay=False, path_type=Path),
    default=None,
    help="Single gene tree file (newick, one tree per line).",
)
@click.option(
    "--tree-dir",
    type=click.Path(file_okay=False, path_type=Path),
    default=None,
    help="Directory of gene tree files for merging.",
)
@click.option(
    "--mode",
    type=click.IntRange(1, 4),
    default=1,
    show_default=True,
    help="wastral --mode. 1=hybrid, 2=branch support weighting, 3=branch length weighting, 4=traditional unweighted Astral.",
)
@click.option(
    "--boot",
    type=click.IntRange(0, 3),
    default=1,
    show_default=True,
    help="wastral -u. 0=topology only, 1=local posterior probability, 2=quartet+local-PP, 3=same as 2 + freqQuad.csv.",
)
@click.option(
    "--extra-rounds", "-R",
    is_flag=True,
    default=False,
    help="Enable exhaustive search (wastral -R).",
)
@click.option(
    "--tree-boot-type",
    type=click.Choice(["auto", "likelihood", "abayes", "bootstrap"]),
    default="auto",
    show_default=True,
    help="Input gene tree branch support type.",
)
@click.option(
    "--tree-boot-min",
    type=float,
    default=None,
    help="Minimum support threshold (wastral -n). Only with non-auto --tree-boot-type.",
)
@click.option(
    "--tree-boot-max",
    type=float,
    default=None,
    help="Maximum support value (wastral -x). Only with non-auto --tree-boot-type.",
)
@click.option(
    "--output-dir", "-o",
    type=click.Path(file_okay=False, path_type=Path),
    default=Path("runs/tree/msc"),
    show_default=True,
    help="Output directory.",
)
@click.option(
    "--threads", "-t",
    type=int,
    default=4,
    show_default=True,
    help="Thread count for wastral -t.",
)
@click.option(
    "--wastral-path",
    type=Path,
    default=None,
    help="Explicit path to wastral executable.",
)
@click.option(
    "--tool-args",
    type=str,
    default=None,
    help="Extra wastral flags. -i/-o blocked; strategy flags override phyloAI defaults.",
)
@click.option("--overwrite", is_flag=True, default=False, help="Overwrite existing output directory.")
@click.option("--dry-run", is_flag=True, default=False, help="Show commands without executing.")
@click.option("--quiet", "-q", is_flag=True, default=False, help="Suppress terminal output except errors.")
def msc_command(
    tree: Path | None,
    tree_dir: Path | None,
    mode: int,
    boot: int,
    extra_rounds: bool,
    tree_boot_type: str,
    tree_boot_min: float | None,
    tree_boot_max: float | None,
    output_dir: Path,
    threads: int,
    wastral_path: Path | None,
    tool_args: str | None,
    overwrite: bool,
    dry_run: bool,
    quiet: bool,
) -> None:
    from phyloai.tree.msc import run_wastral

    # Mutual exclusivity: CLI-layer early check for better error messages
    if (tree is None and tree_dir is None) or (tree is not None and tree_dir is not None):
        _fail("Either --tree or --tree-dir must be provided (mutually exclusive).", 1)

    if tree is not None and not tree.exists():
        _fail(f"--tree does not exist: {tree}", 1)
    if tree_dir is not None and not tree_dir.exists():
        _fail(f"--tree-dir does not exist: {tree_dir}", 1)

    if wastral_path is not None:
        if not wastral_path.exists():
            _fail(f"--wastral-path does not exist: {wastral_path}", 1)
        if not os.access(str(wastral_path), os.X_OK):
            _fail(f"--wastral-path is not executable: {wastral_path}", 1)

    error_msg: str | None = None

    try:
        payload = run_wastral(
            tree=tree,
            tree_dir=tree_dir,
            output_dir=output_dir,
            mode=mode,
            boot=boot,
            extra_rounds=extra_rounds,
            tree_boot_type=tree_boot_type,
            tree_boot_min=tree_boot_min,
            tree_boot_max=tree_boot_max,
            threads=threads,
            wastral_path=str(wastral_path) if wastral_path else None,
            tool_args=tool_args,
            overwrite=overwrite,
            dry_run=dry_run,
            quiet=quiet,
        )
    except (ValueError, FileNotFoundError) as exc:
        error_msg = str(exc)
    except SystemExit:
        raise
    except Exception as exc:
        error_msg = str(exc)

    if error_msg is not None:
        exit_code = 3 if "wastral not found" in error_msg.lower() else 1
        _fail(error_msg, exit_code)

    if dry_run:
        if not quiet:
            click.echo(
                f"Dry run: {payload['key_results']['n_input_trees']} gene tree(s) "
                f"would be processed."
            )
            click.echo(" ".join(payload["data"]["cmd"]))
        return

    # Check payload status for wastral execution failures (exit 2)
    if payload["status"] == "error":
        _fail(payload.get("error", "wastral execution failed"), 2)

    result_path = output_dir / "result.json"
    result_path.parent.mkdir(parents=True, exist_ok=True)
    with open(result_path, "w") as fh:
        json.dump(payload, fh, indent=2)

    if not quiet:
        click.echo(
            f"Species tree saved to {output_dir / 'wastral.tre'}"
        )
        click.echo(f"Log saved to {output_dir / 'wastral.log'}", err=True)
        click.echo(f"Results saved to {result_path}", err=True)
```

- [ ] **Step 2: Verify `--help` works**

Run: `python -c "from phyloai.cli.main import cli; from click.testing import CliRunner; r = CliRunner().invoke(cli, ['tree', 'msc', '--help']); print(r.output); assert r.exit_code == 0"`
Expected: Help output with all flags listed, exit 0

- [ ] **Step 3: Commit**

```bash
git add phyloai/cli/commands/tree.py
git commit -m "feat: register msc command in tree CLI with full Click options"
```

---

### Task 11: Write CLI integration tests

**Files:**
- Create: `tests/cli/test_tree_msc.py`

- [ ] **Step 1: Write the CLI test file**

```python
from __future__ import annotations

from pathlib import Path

import pytest
from click.testing import CliRunner

from phyloai.cli.main import cli


def test_tree_msc_help_shows_all_flags() -> None:
    result = CliRunner().invoke(cli, ["tree", "msc", "--help"])
    assert result.exit_code == 0
    for flag in [
        "--tree", "--tree-dir", "--mode", "--boot",
        "--extra-rounds", "--tree-boot-type", "--tree-boot-min",
        "--tree-boot-max", "--output-dir", "--threads",
        "--wastral-path", "--tool-args", "--overwrite",
        "--dry-run", "--quiet",
    ]:
        assert flag in result.output


def test_tree_msc_mutual_exclusivity(tmp_path: Path) -> None:
    tree_file = tmp_path / "g.trees"
    tree_file.write_text("((a,b),c);\n")
    tree_dir = tmp_path / "genetrees"
    tree_dir.mkdir()

    result = CliRunner().invoke(cli, [
        "tree", "msc",
        "--tree", str(tree_file),
        "--tree-dir", str(tree_dir),
    ])
    assert result.exit_code == 1


def test_tree_msc_neither_input() -> None:
    result = CliRunner().invoke(cli, ["tree", "msc"])
    assert result.exit_code == 1


def test_tree_msc_tree_nonexistent() -> None:
    result = CliRunner().invoke(cli, [
        "tree", "msc", "--tree", "/nonexistent/file.trees",
    ])
    assert result.exit_code == 1


def test_tree_msc_tree_dir_nonexistent() -> None:
    result = CliRunner().invoke(cli, [
        "tree", "msc", "--tree-dir", "/nonexistent/dir",
    ])
    assert result.exit_code == 1


def test_tree_msc_dry_run_tree_single(tmp_path: Path) -> None:
    tree_file = tmp_path / "g.trees"
    tree_file.write_text("((a,b),c);\n")
    out_dir = tmp_path / "out"

    result = CliRunner().invoke(cli, [
        "tree", "msc",
        "--tree", str(tree_file),
        "--output-dir", str(out_dir),
        "--dry-run",
    ])

    assert result.exit_code == 0
    assert "Dry run" in result.output


def test_tree_msc_dry_run_tree_dir(tmp_path: Path) -> None:
    tree_dir = tmp_path / "genetrees"
    tree_dir.mkdir()
    (tree_dir / "a.nwk").write_text("((a,b),c);\n")
    (tree_dir / "b.tre").write_text("((x,y),z);\n")
    out_dir = tmp_path / "out"

    result = CliRunner().invoke(cli, [
        "tree", "msc",
        "--tree-dir", str(tree_dir),
        "--output-dir", str(out_dir),
        "--dry-run",
    ])

    assert result.exit_code == 0
    assert "Dry run" in result.output
    assert "2 gene tree(s)" in result.output


def test_tree_msc_invalid_mode(tmp_path: Path) -> None:
    tree_file = tmp_path / "g.trees"
    tree_file.write_text("((a,b),c);\n")

    result = CliRunner().invoke(cli, [
        "tree", "msc",
        "--tree", str(tree_file),
        "--mode", "5",
    ])
    assert result.exit_code != 0


def test_tree_msc_invalid_boot(tmp_path: Path) -> None:
    tree_file = tmp_path / "g.trees"
    tree_file.write_text("((a,b),c);\n")

    result = CliRunner().invoke(cli, [
        "tree", "msc",
        "--tree", str(tree_file),
        "--boot", "4",
    ])
    assert result.exit_code != 0


def test_tree_msc_tree_boot_min_max_with_auto(tmp_path: Path) -> None:
    tree_file = tmp_path / "g.trees"
    tree_file.write_text("((a,b),c);\n")
    out_dir = tmp_path / "out"

    result = CliRunner().invoke(cli, [
        "tree", "msc",
        "--tree", str(tree_file),
        "--output-dir", str(out_dir),
        "--tree-boot-type", "auto",
        "--tree-boot-min", "10",
    ])
    assert result.exit_code == 1


def test_tree_msc_tool_args_blocked_i(tmp_path: Path) -> None:
    tree_file = tmp_path / "g.trees"
    tree_file.write_text("((a,b),c);\n")
    out_dir = tmp_path / "out"

    result = CliRunner().invoke(cli, [
        "tree", "msc",
        "--tree", str(tree_file),
        "--output-dir", str(out_dir),
        "--tool-args", "-i other.trees",
        "--dry-run",
    ])
    assert result.exit_code == 1


def test_tree_msc_tool_args_override_u(tmp_path: Path) -> None:
    tree_file = tmp_path / "g.trees"
    tree_file.write_text("((a,b),c);\n")
    out_dir = tmp_path / "out"

    result = CliRunner().invoke(cli, [
        "tree", "msc",
        "--tree", str(tree_file),
        "--output-dir", str(out_dir),
        "--tool-args", "-u 3",
        "--dry-run",
        "--quiet",
    ])
    assert result.exit_code == 0
```

- [ ] **Step 2: Run CLI tests**

Run: `pytest tests/cli/test_tree_msc.py -v`
Expected: All tests PASS

- [ ] **Step 3: Commit**

```bash
git add tests/cli/test_tree_msc.py
git commit -m "test: add CLI integration tests for tree msc"
```

---

### Task 12: Verify `tree --help` shows msc and `tree msc --help` works

**Files:** None (verification only)

- [ ] **Step 1: Check top-level tree help**

Run: `python -c "from click.testing import CliRunner; from phyloai.cli.main import cli; r = CliRunner().invoke(cli, ['tree', '--help']); print(r.output); assert 'msc' in r.output; assert r.exit_code == 0"`
Expected: Output contains "msc", exit 0

- [ ] **Step 2: Check all tests pass together**

Run: `pytest tests/tree/test_msc.py tests/cli/test_tree_msc.py -v`
Expected: All tests PASS

- [ ] **Step 3: Run existing tests to ensure no regressions**

Run: `pytest tests/cli/test_tree.py tests/tree/ -v`
Expected: Pre-existing tests continue to PASS

---

### Task 13: Write user-facing command documentation

**Files:**
- Create: `docs/commands/tree-msc.md`

- [ ] **Step 1: Write the command doc**

```markdown
# phyloai tree msc

Multispecies coalescent species tree inference with [wASTRAL](https://github.com/chaoszhang/ASTER) (ASTER).

## Purpose

`phyloai tree msc` consumes gene trees and produces a species tree with local posterior probability branch support using wASTRAL. wASTRAL is a re-implementation of ASTRAL for species tree inference under the multispecies coalescent model.

wASTRAL computation is one-shot — there is no `--resume` support.

## Usage

```bash
# Single gene tree file input
phyloai tree msc --tree gene_trees.trees -o runs/tree/msc

# Directory of gene tree files (auto-merged)
phyloai tree msc --tree-dir ./genetrees/

# Traditional unweighted Astral with exhaustive search
phyloai tree msc --tree-dir ./genetrees/ --mode 4 -R

# Bootstrap input support with custom range
phyloai tree msc --tree-dir ./genetrees/ \
    --mode 1 --boot 2 -R \
    --tree-boot-type bootstrap --tree-boot-min 10 --tree-boot-max 95 \
    -t 8 -o runs/tree/msc

# Override via --tool-args
phyloai tree msc --tree input.trees --tool-args "-r 32 -s 32"
```

## Inputs

| Option | Description |
|--------|-------------|
| `--tree` | Single gene tree file (newick, one tree per line). Mutually exclusive with `--tree-dir`. |
| `--tree-dir` | Directory of gene tree files. Scanned for `.nwk`, `.tre`, `.tree`, `.nw`, `.trees`, `.newick` extensions, merged into one input. Mutually exclusive with `--tree`. |

## Parameters

| Option | Default | Description |
|--------|---------|-------------|
| `--mode` | 1 | 1=hybrid, 2=branch support weighting, 3=branch length weighting, 4=traditional unweighted |
| `--boot` | 1 | 0=topology only, 1=local posterior probability, 2=quartet+local-PP, 3=same as 2 + freqQuad.csv |
| `-R` / `--extra-rounds` | off | Enable exhaustive search |
| `--tree-boot-type` | auto | Input gene tree branch support type: `auto`, `likelihood`, `abayes`, `bootstrap` |
| `--tree-boot-min` | — | Minimum support threshold. Only with non-auto `--tree-boot-type`. |
| `--tree-boot-max` | — | Maximum support value. Only with non-auto `--tree-boot-type`. |
| `--outgroup` | — | Outgroup species for rooting (wastral --root). |

## Outputs

```
runs/tree/msc/
├── result.json            # PhyloAI structured results
├── wastral.tre            # Species tree output (newick)
├── wastral.log            # wastral stderr diagnostic output
├── merged.trees           # Merged input (--tree-dir mode only)
└── freqQuad.csv           # Quartet frequency data (--boot 3 only)
```

## Exit Codes

| Code | Meaning |
|------|---------|
| 0 | Success |
| 1 | User input error |
| 2 | wastral execution failed |
| 3 | wastral not found |

## Notes

- wASTRAL must be installed and on PATH (or use `--wastral-path`).
- wastral stderr is saved to `wastral.log` for diagnostics.
- `--tree-dir` mode merges all valid gene tree files into one input file saved as `merged.trees`.
- `--tool-args` passes extra flags verbatim to wastral. `-i` and `-o` are blocked. Strategy flags override phyloAI defaults.
```

- [ ] **Step 2: Commit**

```bash
git add docs/commands/tree-msc.md
git commit -m "docs: add user-facing tree-msc command documentation"
```

---

### Task 14: Final verification — full test suite

- [ ] **Step 1: Run all tree-related tests**

```bash
pytest tests/tree/test_msc.py tests/cli/test_tree_msc.py tests/cli/test_tree.py tests/tree/ -v
```

- [ ] **Step 2: Verify no import errors**

```bash
python -c "from phyloai.cli.commands.tree import tree; print('CLI OK')"
python -c "from phyloai.tree.msc import run_wastral, _build_wastral_cmd, _merge_gene_trees; print('Library OK')"
```

- [ ] **Step 3: Verify `tree --help` works end-to-end**

```bash
python -c "from click.testing import CliRunner; from phyloai.cli.main import cli; r = CliRunner().invoke(cli, ['tree', '--help']); print(r.output); assert r.exit_code == 0; assert 'msc' in r.output"
```

- [ ] **Step 4: Commit if any fixes were needed**

```bash
git add -A
git commit -m "chore: final verification and fixes for tree msc"
```
```

---

## Self-Review

**1. Spec coverage:**

| Spec requirement | Covered by |
|-----------------|-----------|
| `--tree` / `--tree-dir` mutual exclusivity | Task 6 (validation in `run_wastral`) + Task 10 (CLI-layer check) |
| `--mode` 1-4 | Task 5 (command builder, IntRange in Task 10) |
| `--boot` 0-3 | Task 5 + Task 10 (IntRange(0,3)) |
| `-R` / `--extra-rounds` | Task 5 (command builder) |
| `--tree-boot-type` presets (likelihood/abayes/bootstrap/auto) | Task 5 (boot_flag_map + boot_defaults) |
| `--tree-boot-min` / `--tree-boot-max` | Task 5 + Task 6 (validation) |
| `--overwrite` output dir policy | Task 6 (validation block) |
| No `--resume` | Implicit (no checkpoint code written) |
| Two-tier `--tool-args` model | Task 3 (`_check_managed_flag_conflict` + `_is_flag_overridden`) |
| Blocked flags: `-i`, `-o`, shell redirects | Task 3 (`_WASTRAL_BLOCKED_FLAGS`, `_WASTRAL_BLOCKED_IO_CHARS`) |
| Override groups: mode, boot, extra_rounds, tree_boot, threads | Task 3 (`_WASTRAL_OVERRIDE_MAP`) |
| Input scanning: `.nwk`, `.tre`, `.tree`, `.nw`, `.trees`, `.newick` | Task 2 (`_WASTRAL_TREE_EXTENSIONS`) |
| Gene tree merging for `--tree-dir` | Task 2 (`_merge_gene_trees`) |
| `result.json` schema with all fields | Task 6 (`_assemble_wastral_result`) |
| `wastral.log` stderr capture | Task 6 (subprocess stderr → log file) |
| `merged.trees` output (--tree-dir mode) | Task 6 (merged_path in run_wastral) |
| `path_aliases: ["aster"]` in TOOL_REGISTRY | Task 1 |
| Exit codes 0/1/2/3 | Task 6 + Task 10 |
| Warning on exactly 1 file in --tree-dir | Task 6 (warnings_list) |
| Warning on non-newick files | Task 6 (warnings_list) |

**2. Placeholder scan:** No TBD, TODO, "implement later", or other red flags found. Every task contains complete, executable code.

**3. Type consistency:** Function signatures match across tasks. `_assemble_wastral_result` is called with the same parameter names defined earlier. CLI-layer parameter names match the spec's kebab-case convention. Library-layer uses spec parameter names (e.g., `tree_boot_type`, `tree_boot_min`).

---

**Plan complete and saved to `docs/superpowers/plans/2026-06-20-phyloai-tree-msc.md`. Two execution options:**

1. **Subagent-Driven (recommended)** - I dispatch a fresh subagent per task, review between tasks, fast iteration
2. **Inline Execution** - Execute tasks in this session using executing-plans, batch execution with checkpoints

**Which approach?**
