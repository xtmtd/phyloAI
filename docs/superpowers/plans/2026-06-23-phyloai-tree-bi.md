# phyloai tree bi Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Implement `phyloai tree bi` for PhyloBayes-MPI Bayesian tree inference with multi-chain execution, resume, convergence monitoring, trace plots, and result.json output.

**Architecture:** Add a direct `tree bi` Click command that delegates all behavior to `phyloai/tree/bi.py`. Keep command construction, run-state handling, trace parsing, convergence parsing, and result assembly as small testable functions. Use fake executable scripts in tests so CI does not require real PhyloBayes-MPI or MPI.

**Tech Stack:** Python 3.12+, Click, subprocess, pathlib, json, itertools, matplotlib optional, pytest, Click CliRunner.

---

## File Structure

| File | Action | Responsibility |
|------|--------|----------------|
| `phyloai/core/env.py` | Modify | Add `bpcomp`, `tracecomp`, `readpb_mpi`, and `mpirun` registry entries. |
| `phyloai/cli/doctor.py` | Modify | Keep JSON output flat; text output should naturally include new tools from registry. |
| `phyloai/cli/commands/tree.py` | Modify | Register direct `tree bi` command, grouped help, options, validation, result writing. |
| `phyloai/tree/bi.py` | Create | Library implementation: validation, tool resolution, command building, run_state, monitoring helpers, convergence parsing, result assembly. |
| `tests/tree/test_bi.py` | Create | Unit tests for pure library helpers and fake-tool integration. |
| `tests/cli/test_tree_bi.py` | Create | CLI tests for help, validation, dry-run, resume option parsing. |
| `tests/cli/test_doctor.py` | Modify | Include new PhyloBayes MPI tools in mocked doctor results and JSON checks. |
| `docs/commands/tree-bi.md` | Create | User-facing command docs with examples, resume behavior, stopping advice, output files. |
| `docs/superpowers/specs/2026-06-07-phyloai-design.md` | Modify | Replace stale `tree bi phylobayes` references and output path. |
| `docs/superpowers/specs/2026-06-17-phyloai-tree-design.md` | Modify | Replace stale `bi` group hierarchy with direct command. |

Implementation note: the real observed `bpcomp -o <basename>` summary file is `<basename>.bpdiff`; parse `.bpdiff` and make fake `bpcomp` scripts write `.bpdiff`.

---

### Task 1: Add PhyloBayes-MPI Tools To Registry And Doctor Tests

**Files:**
- Modify: `phyloai/core/env.py:37` (add pb_mpi group table; move existing pb_mpi into `phylobayes_mpi` group)
- Modify: `phyloai/cli/doctor.py` (text mode: group bpcomp/tracecomp/readpb_mpi/mpirun under "PhyloBayes MPI" section; move existing pb_mpi there)
- Modify: `tests/cli/test_doctor.py:10`

- [ ] **Step 1: Extend `TOOL_REGISTRY` — add a `phylobayes_mpi` group and register bpcomp/tracecomp/readpb_mpi/mpirun; move existing pb_mpi into this group**

```python
    "pb_mpi":     {"required": False, "version_flag": "",
                   "install": "https://github.com/bayesiancook/pbmpi"},
    "bpcomp":     {"required": False, "version_flag": "",
                   "install": "https://github.com/bayesiancook/pbmpi"},
    "tracecomp":  {"required": False, "version_flag": "",
                   "install": "https://github.com/bayesiancook/pbmpi"},
    "readpb_mpi": {"required": False, "version_flag": "",
                   "install": "https://github.com/bayesiancook/pbmpi"},
    "mpirun":     {"required": False, "version_flag": "--version",
                   "install": "https://www.open-mpi.org  (or: brew install open-mpi / apt install openmpi-bin)"},
```

- [ ] **Step 2: Extend `_mock_tools()` in `tests/cli/test_doctor.py`**

```python
        "pb_mpi": ToolInfo("pb_mpi", ToolStatus.OK, Path("/opt/pbmpi/bin/pb_mpi"), "1.9"),
        "bpcomp": ToolInfo("bpcomp", ToolStatus.OK, Path("/opt/pbmpi/bin/bpcomp"), "1.9"),
        "tracecomp": ToolInfo("tracecomp", ToolStatus.OK, Path("/opt/pbmpi/bin/tracecomp"), "1.9"),
        "readpb_mpi": ToolInfo("readpb_mpi", ToolStatus.OK, Path("/opt/pbmpi/bin/readpb_mpi"), "1.9"),
        "mpirun": ToolInfo("mpirun", ToolStatus.OK, Path("/usr/local/bin/mpirun"), "4.1.2"),
```

- [ ] **Step 3: Add doctor assertions**

```python
def test_doctor_json_includes_phylobayes_mpi_tools():
    runner = CliRunner()
    with patch("phyloai.cli.doctor.ToolEnv") as MockEnv:
        MockEnv.return_value.check_all.return_value = _mock_tools()
        result = runner.invoke(cli, ["doctor", "--output-format", "json"])

    data = json.loads(result.output)
    for name in ["pb_mpi", "bpcomp", "tracecomp", "readpb_mpi", "mpirun"]:
        assert name in data
        assert data[name]["status"] == "ok"
```

- [ ] **Step 3b: Ensure text-mode doctor groups these tools under a "PhyloBayes MPI" section**

The doctor text output must show a labelled section:

```
PhyloBayes MPI
  pb_mpi      found  /opt/pbmpi/bin/pb_mpi        version: 1.9
  bpcomp      found  /opt/pbmpi/bin/bpcomp         version: 1.9
  tracecomp   found  /opt/pbmpi/bin/tracecomp      version: 1.9
  readpb_mpi  found  /opt/pbmpi/bin/readpb_mpi     version: 1.9
  mpirun      found  /usr/local/bin/mpirun          version: 4.1.2
```

If the existing doctor uses a data-driven group table (e.g. `TOOL_GROUPS`), add a `phylobayes_mpi` group entry and remove `pb_mpi` from its prior group. If the text output is generated by iterating a flat list of groups, add the new group to that list. Update the corresponding doctor text-mode test.

- [ ] **Step 4: Run focused tests**

Run: `pytest tests/cli/test_doctor.py -q`
Expected: all doctor tests pass.

---

### Task 2: Create `phyloai/tree/bi.py` Core Constants And Pure Helpers

**Files:**
- Create: `phyloai/tree/bi.py`
- Create: `tests/tree/test_bi.py`

- [ ] **Step 1: Create library skeleton**

```python
"""Bayesian phylogenetic inference with PhyloBayes-MPI."""

from __future__ import annotations

import itertools
import json
import math
import os
import re
import shutil
import subprocess
import time
from pathlib import Path
from typing import Any

from phyloai.core.env import ToolEnv


MODEL_FLAGS = {
    "gtr": "-gtr",
    "poisson": "-poisson",
    "lg": "-lg",
    "wag": "-wag",
    "jtt": "-jtt",
    "mtrev": "-mtrev",
    "mtzoa": "-mtzoa",
    "mtart": "-mtart",
}

RESUME_ALL = "__ALL__"


def _resolve_chain_names(chains: int, chain_prefix: str, chain_names: str | None) -> list[str]:
    if chain_names:
        names = [item.strip() for item in chain_names.split(",") if item.strip()]
        if not names:
            raise ValueError("--chain-names must contain at least one non-empty name")
        if len(set(names)) != len(names):
            raise ValueError("--chain-names contains duplicate names")
        return names
    if chains < 1:
        raise ValueError("--chains must be at least 1")
    return [f"{chain_prefix}{i}" for i in range(1, chains + 1)]


def _build_model_flags(model: str, mixture: str, gamma_cats: int, start_tree: Path | None, fix_tree: Path | None) -> list[str]:
    if model not in MODEL_FLAGS:
        raise ValueError(f"Invalid --model: {model}")
    if gamma_cats < 1:
        raise ValueError("--gamma-cats must be at least 1")
    flags: list[str] = []
    if mixture == "auto":
        flags.append("-cat")
    else:
        try:
            ncat = int(mixture)
        except ValueError as exc:
            raise ValueError("--mixture must be 'auto' or a positive integer") from exc
        if ncat < 1:
            raise ValueError("--mixture integer must be at least 1")
        flags.extend(["-ncat", str(ncat)])
    flags.append(MODEL_FLAGS[model])
    flags.extend(["-dgam", str(gamma_cats)])
    if start_tree is not None and fix_tree is not None:
        raise ValueError("--start-tree and --fix-tree are mutually exclusive")
    if start_tree is not None:
        flags.extend(["-t", str(start_tree.resolve())])
    if fix_tree is not None:
        flags.extend(["-T", str(fix_tree.resolve())])
    return flags


def _count_trace_samples(trace_path: Path) -> int:
    if not trace_path.exists():
        return 0
    raw = trace_path.read_bytes()
    if not raw:
        return 0
    lines = raw.splitlines(keepends=True)
    complete = [line.decode(errors="ignore").strip() for line in lines if line.endswith((b"\n", b"\r\n"))]
    rows = [line for line in complete if line]
    if len(rows) <= 1:
        return 0
    return max(0, len(rows) - 1)
```

- [ ] **Step 2: Add pure helper tests**

```python
from __future__ import annotations

from pathlib import Path

import pytest

from phyloai.tree.bi import _build_model_flags, _count_trace_samples, _resolve_chain_names


def test_resolve_chain_names_auto():
    assert _resolve_chain_names(3, "chain", None) == ["chain1", "chain2", "chain3"]


def test_resolve_chain_names_explicit():
    assert _resolve_chain_names(3, "chain", "a,b") == ["a", "b"]


def test_build_model_flags_default():
    assert _build_model_flags("gtr", "auto", 4, None, None) == ["-cat", "-gtr", "-dgam", "4"]


def test_build_model_flags_homogeneous_lg():
    assert _build_model_flags("lg", "1", 4, None, None) == ["-ncat", "1", "-lg", "-dgam", "4"]


def test_build_model_flags_rejects_both_trees(tmp_path: Path):
    tree = tmp_path / "t.nwk"
    tree.write_text("(a,b);\n")
    with pytest.raises(ValueError, match="mutually exclusive"):
        _build_model_flags("gtr", "auto", 4, tree, tree)


def test_count_trace_samples_ignores_partial_line(tmp_path: Path):
    trace = tmp_path / "chain1.trace"
    trace.write_bytes(b"iter\ttime\tloglik\n1\t0\t-10\n2\t1\t-9")
    assert _count_trace_samples(trace) == 1
```

- [ ] **Step 3: Run tests**

Run: `pytest tests/tree/test_bi.py -q`
Expected: all tests pass.

---

### Task 3: Add Tool Resolution And Command Builders

**Files:**
- Modify: `phyloai/tree/bi.py`
- Modify: `tests/tree/test_bi.py`

- [ ] **Step 1: Add tool and command helpers**

```python
def _detect_tools(pb_path: Path | None, dry_run: bool) -> dict[str, str]:
    if pb_path is not None:
        tool_paths = {
            "pb_mpi": pb_path / "pb_mpi",
            "bpcomp": pb_path / "bpcomp",
            "tracecomp": pb_path / "tracecomp",
        }
        readpb = pb_path / "readpb_mpi"
        if readpb.exists():
            tool_paths["readpb_mpi"] = readpb
        env = ToolEnv(tool_paths=tool_paths)
    else:
        env = ToolEnv()
    if dry_run:
        return {"pb_mpi": "pb_mpi", "bpcomp": "bpcomp", "tracecomp": "tracecomp", "mpirun": "mpirun"}
    return {name: str(env.require(name)) for name in ["pb_mpi", "bpcomp", "tracecomp", "mpirun"]}


def _build_chain_cmd(mpirun: str, pb_mpi: str, threads: int, matrix: Path, model_flags: list[str], sample_freq: int, nsamples: int, chain_name: str) -> list[str]:
    if threads < 2:
        raise ValueError("--threads must be at least 2")
    if sample_freq < 1:
        raise ValueError("--sample-freq must be at least 1")
    if nsamples != -1 and nsamples < 1:
        raise ValueError("--nsamples must be -1 or a positive integer")
    cmd = [mpirun, "-np", str(threads), pb_mpi, "-d", str(matrix.resolve()), *model_flags, "-x", str(sample_freq), str(nsamples)]
    cmd.append(chain_name)
    return cmd


def _build_resume_cmd(mpirun: str, pb_mpi: str, threads: int, chain_name: str) -> list[str]:
    if threads < 2:
        raise ValueError("--threads must be at least 2")
    return [mpirun, "-np", str(threads), pb_mpi, chain_name]
```

- [ ] **Step 2: Add command builder tests**

```python
from phyloai.tree.bi import _build_chain_cmd, _build_resume_cmd


def test_build_chain_cmd_forever(tmp_path: Path):
    matrix = tmp_path / "m.phy"
    matrix.write_text("2 3\na AAA\nb AAA\n")
    cmd = _build_chain_cmd("mpirun", "pb_mpi", 4, matrix, ["-cat", "-gtr", "-dgam", "4"], 1, -1, "chain1")
    assert cmd[-3:] == ["1", "-1", "chain1"]


def test_build_chain_cmd_with_target(tmp_path: Path):
    matrix = tmp_path / "m.phy"
    matrix.write_text("2 3\na AAA\nb AAA\n")
    cmd = _build_chain_cmd("mpirun", "pb_mpi", 4, matrix, ["-ncat", "1", "-lg", "-dgam", "4"], 1, 10000, "chain1")
    assert cmd[-3:] == ["1", "10000", "chain1"]


def test_build_resume_cmd():
    assert _build_resume_cmd("mpirun", "pb_mpi", 4, "chain1") == ["mpirun", "-np", "4", "pb_mpi", "chain1"]
```

- [ ] **Step 3: Run tests**

Run: `pytest tests/tree/test_bi.py -q`
Expected: all tests pass.

---

### Task 4: Add Run State Handling

**Files:**
- Modify: `phyloai/tree/bi.py`
- Modify: `tests/tree/test_bi.py`

- [ ] **Step 1: Add run_state helpers**

```python
def _state_payload(chain_names: list[str], matrix: Path, model_flags: list[str], sample_freq: int, nsamples: int, threads: int) -> dict[str, Any]:
    return {
        "chain_names": chain_names,
        "matrix": str(matrix.resolve()),
        "model_flags": model_flags,
        "sample_freq": sample_freq,
        "nsamples": nsamples,
        "threads": threads,
    }


def _write_run_state(output_dir: Path, payload: dict[str, Any]) -> None:
    (output_dir / "run_state.json").write_text(json.dumps(payload, indent=2))


def _read_run_state(output_dir: Path) -> dict[str, Any]:
    state_path = output_dir / "run_state.json"
    if not state_path.exists():
        raise ValueError(f"Missing run_state.json in {output_dir}")
    return json.loads(state_path.read_text())


def _update_run_state_for_new_chains(output_dir: Path, new_names: list[str], current_payload: dict[str, Any]) -> dict[str, Any]:
    existing = _read_run_state(output_dir)
    for key in ["matrix", "model_flags", "sample_freq", "nsamples", "threads"]:
        if existing.get(key) != current_payload.get(key):
            raise ValueError("Model parameters conflict with existing run_state.json. Use --resume to continue existing chains or choose a different --output-dir.")
    current_names = list(existing.get("chain_names", []))
    overlap = sorted(set(current_names) & set(new_names))
    if overlap:
        raise ValueError(f"Chain name(s) already exist in run_state.json: {', '.join(overlap)}")
    existing["chain_names"] = current_names + new_names
    _write_run_state(output_dir, existing)
    return existing


def _resolve_resume_names(resume: str | None, state: dict[str, Any]) -> list[str]:
    if resume is None:
        return []
    available = list(state.get("chain_names", []))
    if resume == RESUME_ALL:
        return available
    requested = [item.strip() for item in resume.split(",") if item.strip()]
    missing = sorted(set(requested) - set(available))
    if missing:
        raise ValueError(f"Resume chain(s) not found in run_state.json: {', '.join(missing)}")
    return requested
```

- [ ] **Step 2: Add run_state tests**

```python
from phyloai.tree.bi import _read_run_state, _resolve_resume_names, _state_payload, _update_run_state_for_new_chains, _write_run_state, RESUME_ALL


def test_run_state_roundtrip(tmp_path: Path):
    matrix = tmp_path / "m.phy"
    matrix.write_text("2 3\na AAA\nb AAA\n")
    payload = _state_payload(["chain1"], matrix, ["-cat", "-gtr", "-dgam", "4"], 1, -1, 4)
    _write_run_state(tmp_path, payload)
    assert _read_run_state(tmp_path)["chain_names"] == ["chain1"]


def test_update_run_state_adds_new_chain(tmp_path: Path):
    matrix = tmp_path / "m.phy"
    matrix.write_text("2 3\na AAA\nb AAA\n")
    payload = _state_payload(["chain1"], matrix, ["-cat", "-gtr", "-dgam", "4"], 1, -1, 4)
    _write_run_state(tmp_path, payload)
    updated = _update_run_state_for_new_chains(tmp_path, ["chain2"], payload)
    assert updated["chain_names"] == ["chain1", "chain2"]


def test_resolve_resume_all():
    state = {"chain_names": ["chain1", "chain2"]}
    assert _resolve_resume_names(RESUME_ALL, state) == ["chain1", "chain2"]
```

- [ ] **Step 3: Run tests**

Run: `pytest tests/tree/test_bi.py -q`
Expected: all tests pass.

---

### Task 5: Add Convergence Parsers And Plot Helper

**Files:**
- Modify: `phyloai/tree/bi.py`
- Modify: `tests/tree/test_bi.py`

- [ ] **Step 1: Add parser helpers**

```python
def _parse_bpcomp_bpdiff(path: Path) -> dict[str, float | None]:
    text = path.read_text() if path.exists() else ""
    maxdiff = None
    meandiff = None
    for line in text.splitlines():
        low = line.lower()
        nums = re.findall(r"[-+]?\d*\.\d+|[-+]?\d+", line)
        if "maxdiff" in low and nums:
            maxdiff = float(nums[-1])
        if "meandiff" in low and nums:
            meandiff = float(nums[-1])
    return {"maxdiff": maxdiff, "meandiff": meandiff}


def _parse_tracecomp_contdiff(path: Path) -> dict[str, float | None]:
    if not path.exists():
        return {"min_effsize": None, "max_rel_diff": None}
    min_effsize = None
    max_rel_diff = None
    for line in path.read_text().splitlines():
        if not line.strip() or line.lower().startswith("name"):
            continue
        parts = line.split()
        if len(parts) < 3:
            continue
        try:
            effsize = float(parts[-2])
            rel_diff = float(parts[-1])
        except ValueError:
            continue
        min_effsize = effsize if min_effsize is None else min(min_effsize, effsize)
        max_rel_diff = rel_diff if max_rel_diff is None else max(max_rel_diff, rel_diff)
    return {"min_effsize": min_effsize, "max_rel_diff": max_rel_diff}


def _status_from_metrics(bp_maxdiff: float | None, min_effsize: float | None, max_rel_diff: float | None) -> str:
    if bp_maxdiff is None or min_effsize is None or max_rel_diff is None:
        return "not converged"
    if bp_maxdiff < 0.1 and min_effsize > 300 and max_rel_diff < 0.1:
        return "good"
    if bp_maxdiff < 0.3 and min_effsize > 50 and max_rel_diff < 0.3:
        return "ok"
    return "not converged"


def _bpcomp_status(maxdiff: float | None) -> str:
    if maxdiff is None:
        return "no"
    if maxdiff < 0.1:
        return "good"
    if maxdiff < 0.3:
        return "ok"
    return "no"


def _tracecomp_status(min_effsize: float | None, max_rel_diff: float | None) -> str:
    if min_effsize is None or max_rel_diff is None:
        return "no"
    if min_effsize > 300 and max_rel_diff < 0.1:
        return "good"
    if min_effsize > 50 and max_rel_diff < 0.3:
        return "ok"
    return "no"
```

- [ ] **Step 2: Add parser tests**

```python
from phyloai.tree.bi import _parse_bpcomp_bpdiff, _parse_tracecomp_contdiff, _status_from_metrics


def test_parse_bpcomp_bpdiff(tmp_path: Path):
    path = tmp_path / "bpcomp_all.bpdiff"
    path.write_text("maxdiff 0.081\nmeandiff 0.006\n")
    assert _parse_bpcomp_bpdiff(path) == {"maxdiff": 0.081, "meandiff": 0.006}


def test_parse_tracecomp_contdiff(tmp_path: Path):
    path = tmp_path / "tracecomp_all.contdiff"
    path.write_text("name effsize rel_diff\nloglik 312 0.094\nlength 400 0.050\n")
    assert _parse_tracecomp_contdiff(path) == {"min_effsize": 312.0, "max_rel_diff": 0.094}


def test_status_from_metrics_good():
    assert _status_from_metrics(0.081, 312, 0.094) == "good"
```

- [ ] **Step 3: Add optional plot helper**

```python
def _generate_trace_plots(trace_paths: list[Path], output_pdf: Path, burnin: int) -> bool:
    try:
        from matplotlib.backends.backend_pdf import PdfPages
        import matplotlib.pyplot as plt
    except Exception:
        return False
    rows_by_chain: dict[str, tuple[list[str], list[list[str]]]] = {}
    for path in trace_paths:
        if not path.exists():
            continue
        lines = [line.strip().split() for line in path.read_text().splitlines() if line.strip()]
        if len(lines) <= 1:
            continue
        rows_by_chain[path.stem] = (lines[0], lines[1:])
    if not rows_by_chain:
        return False
    first_header = next(iter(rows_by_chain.values()))[0]
    iter_idx = first_header.index("iter") if "iter" in first_header else 0
    columns = [col for col in first_header if col not in {"iter", "time"}]
    output_pdf.parent.mkdir(parents=True, exist_ok=True)
    with PdfPages(output_pdf) as pdf:
        for column in columns:
            fig, ax = plt.subplots(figsize=(8, 4))
            for chain, (header, rows) in rows_by_chain.items():
                if column not in header:
                    continue
                idx = header.index(column)
                y = []
                x = []
                for row in rows:
                    if idx >= len(row) or iter_idx >= len(row):
                        continue
                    try:
                        y.append(float(row[idx]))
                        x.append(float(row[iter_idx]))
                    except ValueError:
                        continue
                ax.plot(x, y, label=chain)
            ax.axvline(burnin, linestyle="--", color="black", linewidth=0.8)
            ax.set_title(column)
            ax.set_xlabel("iteration")
            ax.set_ylabel(column)
            ax.legend()
            pdf.savefig(fig)
            plt.close(fig)
    return True
```

- [ ] **Step 4: Run tests**

Run: `pytest tests/tree/test_bi.py -q`
Expected: all tests pass.

---

### Task 6: Add Convergence Runner

**Files:**
- Modify: `phyloai/tree/bi.py`
- Modify: `tests/tree/test_bi.py`

- [ ] **Step 1: Add convergence runner**

```python
def _run_convergence_check(output_dir: Path, chain_names: list[str], tools: dict[str, str], burnin: int) -> dict[str, Any]:
    conv_dir = output_dir / "convergence"
    conv_dir.mkdir(parents=True, exist_ok=True)
    comparisons: dict[str, list[str]] = {"all": chain_names}
    for a, b in itertools.combinations(chain_names, 2):
        comparisons[f"{a}_{b}"] = [a, b]
    warnings: list[str] = []
    result: dict[str, Any] = {"all_chains": {}, "pairwise": {}}
    for label, names in comparisons.items():
        base = f"bpcomp_{'all' if label == 'all' else label}"
        bp_base = conv_dir / base
        bp_cmd = [tools["bpcomp"], "-x", str(burnin), "-o", base, *[str(Path("../chains") / name) for name in names]]
        bp_proc = subprocess.run(bp_cmd, cwd=conv_dir, capture_output=True, text=True)
        if bp_proc.returncode != 0:
            warnings.append(f"bpcomp {label} exited with code {bp_proc.returncode}: {bp_proc.stderr[:200]}")
        trace_base = f"tracecomp_{'all' if label == 'all' else label}"
        trace_out = conv_dir / f"{trace_base}.contdiff"
        trace_cmd = [tools["tracecomp"], "-x", str(burnin), *[str(Path("../chains") / f"{name}.trace") for name in names]]
        trace_proc = subprocess.run(trace_cmd, cwd=conv_dir, capture_output=True, text=True)
        if trace_proc.returncode != 0:
            warnings.append(f"tracecomp {label} exited with code {trace_proc.returncode}: {trace_proc.stderr[:200]}")
        trace_out.write_text((trace_proc.stdout or "") + (("\n" + trace_proc.stderr) if trace_proc.stderr else ""))
        bp = _parse_bpcomp_bpdiff(bp_base.with_suffix(".bpdiff"))
        tr = _parse_tracecomp_contdiff(trace_out)
        metrics = {
            "bpcomp_maxdiff": bp["maxdiff"],
            "bpcomp_meandiff": bp["meandiff"],
            "tracecomp_min_effsize": tr["min_effsize"],
            "tracecomp_max_reldiff": tr["max_rel_diff"],
            "status": _status_from_metrics(bp["maxdiff"], tr["min_effsize"], tr["max_rel_diff"]),
        }
        if label == "all":
            result["all_chains"] = metrics
        else:
            result["pairwise"][label] = metrics
    try:
        _generate_trace_plots([output_dir / "chains" / f"{name}.trace" for name in chain_names], conv_dir / "trace_plots.pdf", burnin)
    except Exception:
        pass
    result["warnings"] = warnings
    return result
```

- [ ] **Step 2: Add fake convergence tool test**

```python
def test_run_convergence_check_with_fake_tools(tmp_path: Path):
    from phyloai.tree.bi import _run_convergence_check
    chains = tmp_path / "chains"
    chains.mkdir()
    for name in ["chain1", "chain2"]:
        (chains / f"{name}.trace").write_text("iter time loglik\n1 0 -10\n2 1 -9\n")
    bpcomp = tmp_path / "bpcomp"
    bpcomp.write_text("#!/bin/sh\nwhile [ $# -gt 0 ]; do if [ \"$1\" = \"-o\" ]; then shift; out=$1; fi; shift; done\nprintf 'maxdiff 0.081\\nmeandiff 0.006\\n' > ${out}.bpdiff\nprintf '(a,b);\\n' > ${out}.con.tre\nprintf 'split\\n' > ${out}.bplist\n")
    tracecomp = tmp_path / "tracecomp"
    tracecomp.write_text("#!/bin/sh\nprintf 'name effsize rel_diff\\nloglik 312 0.094\\n'\n")
    bpcomp.chmod(0o755)
    tracecomp.chmod(0o755)
    result = _run_convergence_check(tmp_path, ["chain1", "chain2"], {"bpcomp": str(bpcomp), "tracecomp": str(tracecomp)}, 1, True)
    assert result["all_chains"]["status"] == "good"
    assert (tmp_path / "convergence" / "bpcomp_all.bpdiff").exists()
```

- [ ] **Step 3: Run tests**

Run: `pytest tests/tree/test_bi.py::test_run_convergence_check_with_fake_tools -q`
Expected: test passes.

---

### Task 7: Add `run_bi()` Dry-Run And Result Assembly

**Files:**
- Modify: `phyloai/tree/bi.py`
- Modify: `tests/tree/test_bi.py`

- [ ] **Step 1: Add FASTA conversion, version detection, result assembly, and directory-safe `run_bi()`**

```python
def _detect_pb_version(pb_dir: str) -> str | None:
    """File-based version heuristic for PhyloBayes tools.
    Searches the tool's parent directory for version-bearing filenames
    (e.g. pb_mpiManual1.9.pdf, VERSION, CHANGELOG, *.pdf) using regex
    [Vv]ersion\s*(\d+\.\d+[\.\d]*) or (\d+\.\d+) against filenames.
    Returns first match or None."""
    root = Path(pb_dir).parent
    candidates = sorted(root.glob("*"), key=lambda p: str(p).lower())
    for path in candidates:
        if not path.is_file():
            continue
        for pattern in [r"[Vv]ersion\s*(\d+\.\d+(?:\.\d+)?)", r"(\d+\.\d+)"]:
            m = re.search(pattern, path.name)
            if m:
                return m.group(1)
    return None


def _detect_mpirun_version(mpirun_path: str) -> str | None:
    try:
        proc = subprocess.run([mpirun_path, "--version"], capture_output=True, text=True, timeout=5)
        if proc.returncode == 0:
            m = re.search(r"(\d+\.\d+[\.\d]*)", proc.stdout[:200])
            if m:
                return m.group(1)
    except Exception:
        pass
    return None


def _detect_tool_versions(tools: dict[str, str]) -> dict[str, str | None]:
    """Detect versions for pb_mpi, bpcomp, tracecomp, and mpirun.
    PhyloBayes tools use file-based heuristic; mpirun uses --version."""
    pb_ver = _detect_pb_version(tools["pb_mpi"])
    mpirun_ver = _detect_mpirun_version(tools["mpirun"])
    return {
        "pb_mpi": pb_ver,
        "bpcomp": pb_ver,       # same distro as pb_mpi
        "tracecomp": pb_ver,
        "mpirun": mpirun_ver,
    }


def _prepare_matrix(matrix: Path, output_dir: Path, dry_run: bool) -> Path:
    """Detect input format; convert FASTA to PHYLIP if needed.
    pb_mpi only supports PHYLIP natively."""
    from phyloai.core.formats import FormatConverter, AlignmentFormat

    converter = FormatConverter()
    fmt = converter.detect(matrix)
    if fmt == AlignmentFormat.PHYLIP:
        return matrix
    if fmt == AlignmentFormat.FASTA:
        if dry_run:
            return matrix  # dry-run skips conversion
        phylip_path = output_dir / "matrix.phy"
        converter.convert(matrix, phylip_path, target=AlignmentFormat.PHYLIP)
        return phylip_path
    raise ValueError(f"Unsupported matrix format: {fmt}. Expected FASTA or PHYLIP.")


def _assemble_result(params: dict[str, Any], command: str, wall_time: float, tool_versions: dict[str, str | None], chain_cmds: dict[str, list[str]], chain_lengths: dict[str, int], final_convergence: dict[str, Any] | None, tool_outputs: dict[str, str], interrupted: bool, status: str = "success", error: str | None = None, warnings: list[str] | None = None) -> dict[str, Any]:
    # Merge convergence warnings if present
    all_warnings = list(warnings or [])
    if final_convergence:
        all_warnings.extend(final_convergence.get("warnings", []))
    return {
        "status": status,
        "command": command,
        "wall_time": wall_time,
        "tool_versions": tool_versions,
        "params": params,
        "key_results": {
            "chain_names": list(chain_cmds.keys()),
            "chain_lengths": chain_lengths,
            "final_convergence": {k: v for k, v in (final_convergence or {}).items() if k != "warnings"},
            "consensus_tree": "convergence/bpcomp_all.con.tre",
        },
        "error": error,
        "data": {
            "chain_cmds": chain_cmds,
            "tool_stderr": tool_outputs,
            "tool_logs": {name: f"chains/{name}.log" for name in chain_cmds},
            "interrupted": interrupted,
            "warnings": all_warnings,
        },
    }


def run_bi(
    matrix: Path | None,
    output_dir: Path = Path("runs/tree/bi"),
    overwrite: bool = False,
    model: str = "gtr",
    mixture: str = "auto",
    gamma_cats: int = 4,
    start_tree: Path | None = None,
    fix_tree: Path | None = None,
    chains: int = 3,
    chain_prefix: str = "chain",
    chain_names: str | None = None,
    threads: int = 4,
    sample_freq: int = 1,
    nsamples: int = -1,
    resume: str | None = None,
    monitor_freq: int = 100,
    burnin_frac: float = 0.5,
    poll_interval: int = 60,
    pb_path: Path | None = None,
    dry_run: bool = False,
    quiet: bool = False,
) -> dict[str, Any]:
    run_start = time.monotonic()
    if resume is None and matrix is None:
        raise ValueError("--matrix is required unless --resume is used")
    if matrix is not None and not matrix.exists():
        raise ValueError(f"--matrix does not exist: {matrix}")
    if overwrite and resume is not None:
        raise ValueError("--overwrite and --resume are mutually exclusive")
    model_flags = _build_model_flags(model, mixture, gamma_cats, start_tree, fix_tree)
    tools = _detect_tools(pb_path, dry_run)
    tool_versions = _detect_tool_versions(tools) if not dry_run else {"pb_mpi": None, "bpcomp": None, "tracecomp": None, "mpirun": None}
    if resume is None:
        names = _resolve_chain_names(chains, chain_prefix, chain_names)
        assert matrix is not None
        state = _state_payload(names, matrix, model_flags, sample_freq, nsamples, threads)
    else:
        state = _read_run_state(output_dir)
        names = _resolve_resume_names(resume, state)
        matrix = Path(state["matrix"])
        model_flags = list(state["model_flags"])
        sample_freq = int(state["sample_freq"])
        nsamples = int(state["nsamples"])
        threads = int(state["threads"])
    if not dry_run:
        # --- conflict check BEFORE any filesystem writes ---
        if not overwrite and output_dir.exists():
            if resume is not None:
                pass  # --resume always allowed in existing dir
            elif chain_names and (output_dir / "run_state.json").exists():
                pass  # adding chains to existing run
            elif any(output_dir.iterdir()):
                raise ValueError(f"Output directory {output_dir} already exists and is non-empty. Use --overwrite to replace.")
        # --- now safe to create directory structure ---
        if overwrite and output_dir.exists():
            shutil.rmtree(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        (output_dir / "chains").mkdir(exist_ok=True)
        (output_dir / "convergence").mkdir(exist_ok=True)
        if resume is None:
            if (output_dir / "run_state.json").exists() and chain_names:
                _update_run_state_for_new_chains(output_dir, names, state)
            else:
                _write_run_state(output_dir, state)
    # --- auto-convert FASTA to PHYLIP for pb_mpi ---
    if matrix is not None:
        matrix_for_chains = _prepare_matrix(matrix, output_dir, dry_run)
    else:
        matrix_for_chains = matrix
    chain_cmds: dict[str, list[str]] = {}
    for name in names:
        if resume is None:
            chain_cmds[name] = _build_chain_cmd(tools["mpirun"], tools["pb_mpi"], threads, matrix_for_chains, model_flags, sample_freq, nsamples, name)
        else:
            chain_cmds[name] = _build_resume_cmd(tools["mpirun"], tools["pb_mpi"], threads, name)
    params = {
        "matrix": str(matrix) if matrix else None,
        "output_dir": str(output_dir),
        "overwrite": overwrite,
        "model": model,
        "mixture": mixture,
        "gamma_cats": gamma_cats,
        "start_tree": str(start_tree) if start_tree else None,
        "fix_tree": str(fix_tree) if fix_tree else None,
        "chains": chains,
        "chain_prefix": chain_prefix,
        "chain_names": chain_names,
        "threads": threads,
        "sample_freq": sample_freq,
        "nsamples": nsamples,
        "resume": resume,
        "monitor_freq": monitor_freq,
        "burnin_frac": burnin_frac,
        "poll_interval": poll_interval,
        "pb_path": str(pb_path) if pb_path else None,
        "dry_run": dry_run,
        "quiet": quiet,
    }
    command = "phyloai tree bi " + " ".join(["--matrix", str(matrix), "--output-dir", str(output_dir), "--model", model, "--mixture", mixture, "--gamma-cats", str(gamma_cats), "--chains", str(chains), "--chain-prefix", chain_prefix, "--threads", str(threads), "--sample-freq", str(sample_freq), "--nsamples", str(nsamples), "--monitor-freq", str(monitor_freq), "--burnin-frac", str(burnin_frac), "--poll-interval", str(poll_interval)])
    if dry_run:
        return _assemble_result(params, command, time.monotonic() - run_start, {"pb_mpi": None, "bpcomp": None, "tracecomp": None, "mpirun": None}, chain_cmds, {name: 0 for name in names}, None, {name: "" for name in names}, False)
    return _run_bi_processes(output_dir, names, chain_cmds, tools, tool_versions, params, command, run_start, nsamples, monitor_freq, burnin_frac, poll_interval, quiet)
```

- [ ] **Step 2: Add dry-run test**

```python
from tests.helpers import validate_params_completeness, validate_result_json
from phyloai.tree.bi import run_bi


def test_run_bi_dry_run_result_json_shape(tmp_path: Path):
    matrix = tmp_path / "m.phy"
    matrix.write_text("2 3\na AAA\nb AAA\n")
    payload = run_bi(matrix=matrix, output_dir=tmp_path / "out", dry_run=True)
    assert payload["status"] == "success"
    assert payload["data"]["chain_cmds"]["chain1"][-1] == "chain1"
    validate_result_json({**payload, "data": {**payload["data"], "tool_stderr": ""}})
    validate_params_completeness(payload, {"matrix", "output_dir", "overwrite", "model", "mixture", "gamma_cats", "start_tree", "fix_tree", "chains", "chain_prefix", "chain_names", "threads", "sample_freq", "nsamples", "resume", "monitor_freq", "burnin_frac", "poll_interval", "pb_path", "dry_run", "quiet"})
```

- [ ] **Step 3: Run tests**

Run: `pytest tests/tree/test_bi.py::test_run_bi_dry_run_result_json_shape -q`
Expected: test passes.

---

### Task 8: Add Process Launch And Monitoring Loop

**Files:**
- Modify: `phyloai/tree/bi.py`
- Modify: `tests/tree/test_bi.py`

- [ ] **Step 1: Add soft-stop, Rich Live progress display, and Ctrl+C-safe process runner**

```python
def _soft_stop_chains(output_dir: Path, chain_names: list[str]) -> None:
    for name in chain_names:
        (output_dir / "chains" / f"{name}.run").write_text("0\n")


def _format_convergence_text(conv_result: dict[str, Any]) -> str:
    """Render convergence stats as ASCII-only text. 6-column pairwise table with
    per-column bpcomp/tracecomp status."""
    lines: list[str] = []
    all_c = conv_result.get("all_chains", {})
    if all_c:
        bpmax = all_c.get("bpcomp_maxdiff")
        bpmean = all_c.get("bpcomp_meandiff")
        eff = all_c.get("tracecomp_min_effsize")
        rel = all_c.get("tracecomp_max_reldiff")
        bp_st = _bpcomp_status(bpmax)
        tc_st = _tracecomp_status(eff, rel)
        bp = f"  bpcomp    maxdiff  {bpmax:.3f}" if bpmax is not None else "  bpcomp    maxdiff  --"
        if bpmean is not None:
            bp += f"   meandiff  {bpmean:.3f}"
        tr = f"  tracecomp  min effsize  {eff:.0f}" if eff is not None else "  tracecomp  min effsize  --"
        if rel is not None:
            tr += f"   max rel_diff  {rel:.3f}"
        lines.append(f"\n  All chains")
        lines.append(bp + f"   [{bp_st}]")
        lines.append(tr + f"   [{tc_st}]")
    pw = conv_result.get("pairwise", {})
    if pw:
        lines.append("\n  Pairwise")
        lines.append("    pair              maxdiff  min effsize  max rel_diff  bpcomp  tracecomp")
        for label, m in pw.items():
            md = m.get("bpcomp_maxdiff")
            es = m.get("tracecomp_min_effsize")
            rd = m.get("tracecomp_max_reldiff")
            md_s = f"{md:.3f}" if md is not None else "--"
            es_s = f"{es:.0f}" if es is not None else "--"
            rd_s = f"{rd:.3f}" if rd is not None else "--"
            bp_st2 = _bpcomp_status(md)
            tc_st2 = _tracecomp_status(es, rd)
            pair_name = label.replace("_", " x ")
            lines.append(f"    {pair_name:<18} {md_s:^7}  {es_s:^11}  {rd_s:^12}  {bp_st2:^6}  {tc_st2:^9}")
    return "\n".join(lines)


def _run_bi_processes(output_dir: Path, chain_names: list[str], chain_cmds: dict[str, list[str]], tools: dict[str, str], tool_versions: dict[str, str | None], params: dict[str, Any], command: str, run_start: float, nsamples: int, monitor_freq: int, burnin_frac: float, poll_interval: int, quiet: bool) -> dict[str, Any]:
    procs: dict[str, subprocess.Popen[str]] = {}
    outputs: dict[str, str] = {name: "" for name in chain_names}
    try:
        # Build Rich progress bars
        from rich.console import Console
        from rich.progress import BarColumn, Progress, TaskProgressColumn, TextColumn, TimeElapsedColumn
        from rich.live import Live

        console = Console()
        progress = Progress(
            TextColumn(" {task.description}"),
            BarColumn(),
            TaskProgressColumn(),
            TextColumn(" {task.fields[samples]} samples"),
            TimeElapsedColumn(),
            console=console,
        )
        for name in chain_names:
            init = 0
            trace_path = output_dir / "chains" / f"{name}.trace"
            if trace_path.exists():
                init = _count_trace_samples(trace_path)
            progress.add_task(name, samples=init,
                              total=None if nsamples == -1 else nsamples,
                              completed=init if nsamples != -1 else None)
        live_display = Live(progress, console=console, refresh_per_second=1 / 60)
        if not quiet:
            live_display.start()

        for name, cmd in chain_cmds.items():
            log = open(output_dir / "chains" / f"{name}.log", "w")
            proc = subprocess.Popen(cmd, cwd=output_dir / "chains", stdout=log, stderr=subprocess.STDOUT, text=True)
            proc._phyloai_log = log  # type: ignore[attr-defined]
            procs[name] = proc

        last_check = 0
        trace_lengths: dict[str, int] = {}
        for name in chain_names:
            trace_path = output_dir / "chains" / f"{name}.trace"
            trace_lengths[name] = _count_trace_samples(trace_path) if trace_path.exists() else 0
        last_trace_read = time.monotonic()  # delay first poll by --poll-interval
        final_convergence: dict[str, Any] | None = None
        while procs:
            now = time.monotonic()
            # --- Poll .trace files every --poll-interval seconds ---
            if now - last_trace_read >= poll_interval:
                trace_lengths = {name: _count_trace_samples(output_dir / "chains" / f"{name}.trace") for name in chain_names}
                last_trace_read = now
                # Update progress bars
                for tid, task in enumerate(progress.tasks):
                    name = chain_names[tid]
                    progress.update(task.id, samples=trace_lengths.get(name, 0),
                                    completed=trace_lengths.get(name, 0) if nsamples != -1 else None,
                                    total=nsamples if nsamples != -1 else None)
                # Stop chains that reached nsamples target
                if nsamples != -1:
                    for name in list(procs):
                        if trace_lengths.get(name, 0) >= nsamples:
                            (output_dir / "chains" / f"{name}.run").write_text("0\n")
                # Convergence check
                min_len = min(trace_lengths.values()) if trace_lengths else 0
                if min_len - last_check >= monitor_freq:
                    burnin = math.floor(min_len * burnin_frac)
                    if burnin >= 10:
                        conv = _run_convergence_check(output_dir, chain_names, tools, burnin)
                        final_convergence = conv
                        if not quiet:
                            conv_text = _format_convergence_text(conv)
                            live_display.stop()
                            console.print(f"\n--- Convergence Check @ {min_len} samples (burnin {int(burnin_frac*100)}% = {burnin}) ---{conv_text}")
                            console.print("-" * 60)
                            pairwise_statuses = [m.get("status") for m in conv.get("pairwise", {}).values()]
                            all_status = conv.get("all_chains", {}).get("status")
                            n_good = pairwise_statuses.count("good")
                            n_ok = pairwise_statuses.count("ok")
                            n_not = pairwise_statuses.count("not converged")
                            if all_status == "good" and n_not == 0 and n_ok == 0:
                                console.print("  *** All convergence criteria met (all pairs good). You may stop chains with Ctrl+C when ready. ***")
                            elif n_not == 0:
                                console.print(f"  Convergence acceptable across all chain pairs ({n_good} good, {n_ok} ok). Consider stopping when ready.")
                            elif n_good + n_ok >= 1:
                                console.print(f"  Some chain pairs agree ({n_good} good, {n_ok} ok, {n_not} not converged).")
                            live_display.start()
                    last_check = min_len
            # --- Fast poll: process exit checks every 1s for responsiveness ---
            for name, proc in list(procs.items()):
                if proc.poll() is not None:
                    getattr(proc, "_phyloai_log").close()
                    if proc.returncode != 0:
                        remaining = [n for n in procs if n != name]
                        if remaining:
                            if not quiet:
                                live_display.console.print(f"  Chain {name} failed -- stopping remaining chains...")
                            _soft_stop_chains(output_dir, remaining)
                            for rem_name in remaining:
                                procs[rem_name].wait()
                                getattr(procs[rem_name], "_phyloai_log").close()
                        if not quiet:
                            live_display.stop()
                        # Read trace one last time for accurate chain_lengths in error result
                        trace_lengths = {n: _count_trace_samples(output_dir / "chains" / f"{n}.trace") for n in chain_names}
                        return _assemble_result(params, command, time.monotonic() - run_start, tool_versions, chain_cmds, trace_lengths, final_convergence, outputs, False, "error", f"pb_mpi chain {name} exited with code {proc.returncode}")
                    del procs[name]
            time.sleep(1)
        # All chains exited normally
        if not quiet:
            live_display.stop()
        lengths = {name: _count_trace_samples(output_dir / "chains" / f"{name}.trace") for name in chain_names}
        if lengths:
            burnin = math.floor(min(lengths.values()) * burnin_frac)
            if burnin >= 10:
                final_convergence = _run_convergence_check(output_dir, chain_names, tools, burnin)
        for name in chain_names:
            log_path = output_dir / "chains" / f"{name}.log"
            outputs[name] = log_path.read_text() if log_path.exists() else ""
        return _assemble_result(params, command, time.monotonic() - run_start, tool_versions, chain_cmds, lengths, final_convergence, outputs, False)
    except KeyboardInterrupt:
        # --- Ctrl+C soft-stop flow ---
        if not quiet:
            live_display.stop()
            console.print("\n  Caught interrupt -- sending soft-stop to all chains...")
        _soft_stop_chains(output_dir, chain_names)
        if not quiet:
            for name in chain_names:
                console.print(f"    Wrote 0 -> chains/{name}.run")
            console.print("    Waiting for chains to finish current cycle...")
        # Wait for all subprocesses to exit
        for name, proc in procs.items():
            proc.wait()
            getattr(proc, "_phyloai_log").close()
            if not quiet:
                console.print(f"    {name} stopped at {_count_trace_samples(output_dir / 'chains' / f'{name}.trace')} samples.")
        if not quiet:
            live_display.stop()
            console.print("    Running final convergence check...")
        lengths = {name: _count_trace_samples(output_dir / "chains" / f"{name}.trace") for name in chain_names}
        final_convergence = None
        if lengths:
            burnin = math.floor(min(lengths.values()) * burnin_frac)
            if burnin >= 10:
                final_convergence = _run_convergence_check(output_dir, chain_names, tools, burnin)
        for name in chain_names:
            log_path = output_dir / "chains" / f"{name}.log"
            outputs[name] = log_path.read_text() if log_path.exists() else ""
        if not quiet:
            console.print("    Writing result.json  (status: success)")
        return _assemble_result(params, command, time.monotonic() - run_start, tool_versions, chain_cmds, lengths, final_convergence, outputs, True)
```

- [ ] **Step 2: Add fake pb_mpi integration test**

```python
def test_run_bi_fake_tools_executes_chains(tmp_path: Path):
    from phyloai.tree.bi import run_bi
    tool_dir = tmp_path / "tools"
    tool_dir.mkdir()
    mpirun = tool_dir / "mpirun"
    mpirun.write_text("#!/bin/sh\nshift 2\nexec \"$@\"\n")
    pb = tool_dir / "pb_mpi"
    pb.write_text("#!/bin/sh\nfor last do :; done\nname=\"$last\"\nprintf 'iter time loglik\\n1 0 -10\\n2 1 -9\\n' > ${name}.trace\nprintf '(a,b);\\n' > ${name}.treelist\nprintf 'state\\n' > ${name}.chain\nprintf 'stdout from pb\\n'\n")
    bp = tool_dir / "bpcomp"
    bp.write_text("#!/bin/sh\nwhile [ $# -gt 0 ]; do if [ \"$1\" = \"-o\" ]; then shift; out=$1; fi; shift; done\nprintf 'maxdiff 0.081\\nmeandiff 0.006\\n' > ${out}.bpdiff\nprintf '(a,b);\\n' > ${out}.con.tre\nprintf 'split\\n' > ${out}.bplist\n")
    tr = tool_dir / "tracecomp"
    tr.write_text("#!/bin/sh\nprintf 'name effsize rel_diff\\nloglik 312 0.094\\n'\n")
    for path in [mpirun, pb, bp, tr]:
        path.chmod(0o755)
    matrix = tmp_path / "m.phy"
    matrix.write_text("2 3\na AAA\nb AAA\n")
    payload = run_bi(matrix=matrix, output_dir=tmp_path / "out", chains=2, nsamples=2, monitor_freq=1, burnin_frac=0.5, pb_path=tool_dir, quiet=True)
    assert payload["status"] == "success"
    assert payload["key_results"]["chain_lengths"] == {"chain1": 2, "chain2": 2}
```

- [ ] **Step 3: Run fake integration test**

Run: `pytest tests/tree/test_bi.py::test_run_bi_fake_tools_executes_chains -q`
Expected: test passes.

---

### Task 9: Register CLI Command

**Files:**
- Modify: `phyloai/cli/commands/tree.py`
- Create: `tests/cli/test_tree_bi.py`

- [ ] **Step 1: Add `bi` to `_TreeGroup.list_commands()`**

```python
class _TreeGroup(click.Group):
    def list_commands(self, ctx: click.Context) -> list[str]:
        return ["ml", "bi", "msc", "cf"]
```

- [ ] **Step 2: Add command before `msc_command`**

```python
@tree.command("bi", cls=_GroupedHelpCommand, help="Bayesian phylogenetic inference with PhyloBayes-MPI (pb_mpi).")
@click.option("--matrix", "-m", type=click.Path(dir_okay=False, path_type=Path), default=None, help="Input alignment, PHYLIP or FASTA.")
@click.option("--output-dir", "-o", type=click.Path(file_okay=False, path_type=Path), default=Path("runs/tree/bi"), show_default=True, help="Output directory.")
@click.option("--overwrite", is_flag=True, default=False, help="Overwrite existing output directory.")
@click.option("--model", type=click.Choice(["gtr", "poisson", "lg", "wag", "jtt", "mtrev", "mtzoa", "mtart"]), default="gtr", show_default=True, help="Relative exchangeabilities.")
@click.option("--mixture", type=str, default="auto", show_default=True, help="Profile mixture: auto or positive integer.")
@click.option("--gamma-cats", type=click.IntRange(1, None), default=4, show_default=True, help="Discrete gamma categories.")
@click.option("--start-tree", type=click.Path(dir_okay=False, path_type=Path), default=None, help="Starting tree in Newick format.")
@click.option("--fix-tree", type=click.Path(dir_okay=False, path_type=Path), default=None, help="Fixed topology tree in Newick format.")
@click.option("--chains", type=click.IntRange(1, None), default=3, show_default=True, help="Number of independent chains.")
@click.option("--chain-prefix", type=str, default="chain", show_default=True, help="Prefix for auto chain names.")
@click.option("--chain-names", type=str, default=None, help="Explicit comma-separated chain names.")
@click.option("--threads", "-t", type=click.IntRange(2, None), default=4, show_default=True, help="MPI processes per chain.")
@click.option("--sample-freq", type=click.IntRange(1, None), default=1, show_default=True, help="pb_mpi -x every value.")
@click.option("--nsamples", type=int, default=-1, show_default=True, help="Stop after N saved points; -1 runs forever.")
@click.option("--monitor-freq", type=click.IntRange(1, None), default=100, show_default=True, help="Convergence check interval in samples.")
@click.option("--burnin-frac", type=click.FloatRange(0.0, 0.95), default=0.5, show_default=True, help="Dynamic burn-in fraction.")
@click.option("--poll-interval", type=click.IntRange(1, None), default=60, show_default=True, help="Seconds between trace file polls.")
@click.option("--no-plot", is_flag=True, default=False, help="Disable trace_plots.pdf generation.")
@click.option("--resume", default=None, is_flag=False, flag_value="__ALL__", help="Resume all chains, or comma-separated selected chains.")
@click.option("--pb-path", type=click.Path(file_okay=False, path_type=Path), default=None, help="Directory containing pb_mpi, bpcomp, tracecomp.")
@click.option("--dry-run", is_flag=True, default=False, help="Print commands without executing.")
@click.option("--quiet", "-q", is_flag=True, default=False, help="Suppress terminal output except errors.")
def bi_command(**kwargs) -> None:
    from phyloai.tree.bi import run_bi

    try:
        payload = run_bi(**kwargs)
    except FileNotFoundError as exc:
        _fail(str(exc), 3)
    except ValueError as exc:
        _fail(str(exc), 1)
    if kwargs.get("dry_run"):
        if not kwargs.get("quiet"):
            for cmd in payload["data"]["chain_cmds"].values():
                click.echo(" ".join(cmd))
        return
    result_path = kwargs["output_dir"] / "result.json"
    result_path.parent.mkdir(parents=True, exist_ok=True)
    if result_path.exists():
        from datetime import datetime
        ts = datetime.now().strftime("%Y%m%d_%H%M%S")
        backup_path = result_path.with_name(f"result_{ts}.json")
        import shutil
        shutil.copy2(str(result_path), str(backup_path))
    with open(result_path, "w") as fh:
        json.dump(payload, fh, indent=2)
    if payload["status"] == "error":
        _fail(payload.get("error", "pb_mpi execution failed"), 2)
    if not kwargs.get("quiet"):
        click.echo(f"Results saved to {result_path}", err=True)
```

- [ ] **Step 3: Add CLI tests**

```python
from __future__ import annotations

from pathlib import Path

from click.testing import CliRunner

from phyloai.cli.main import cli


def test_tree_bi_help_shows_flags():
    result = CliRunner().invoke(cli, ["tree", "bi", "--help"])
    assert result.exit_code == 0
    for flag in ["--matrix", "--model", "--mixture", "--gamma-cats", "--chains", "--threads", "--resume", "--pb-path", "--poll-interval"]:
        assert flag in result.output


def test_tree_group_shows_bi():
    result = CliRunner().invoke(cli, ["tree", "--help"])
    assert result.exit_code == 0
    assert "bi" in result.output


def test_tree_bi_dry_run(tmp_path: Path):
    matrix = tmp_path / "m.phy"
    matrix.write_text("2 3\na AAA\nb AAA\n")
    result = CliRunner().invoke(cli, ["tree", "bi", "--matrix", str(matrix), "--output-dir", str(tmp_path / "out"), "--dry-run"])
    assert result.exit_code == 0
    assert "pb_mpi" in result.output


def test_tree_bi_resume_bare_parses(tmp_path: Path):
    out = tmp_path / "out"
    out.mkdir()
    (out / "run_state.json").write_text('{"chain_names": ["chain1"], "matrix": "m.phy", "model_flags": ["-cat", "-gtr", "-dgam", "4"], "sample_freq": 1, "nsamples": -1, "threads": 4}')
    result = CliRunner().invoke(cli, ["tree", "bi", "--output-dir", str(out), "--resume", "--dry-run", "--quiet"])
    assert result.exit_code == 0
```

- [ ] **Step 4: Run CLI tests**

Run: `pytest tests/cli/test_tree_bi.py -q`
Expected: all tests pass.

---

### Task 10: Add User Documentation And Update Stale Design References

**Files:**
- Create: `docs/commands/tree-bi.md`
- Modify: `docs/superpowers/specs/2026-06-07-phyloai-design.md`
- Modify: `docs/superpowers/specs/2026-06-17-phyloai-tree-design.md`

- [ ] **Step 1: Create `docs/commands/tree-bi.md`**

````markdown
# phyloai tree bi

Bayesian phylogenetic inference with PhyloBayes-MPI.

## Usage

```bash
phyloai tree bi --matrix concat/matrix.phy
```

## Common Examples

```bash
phyloai tree bi --matrix concat/matrix.phy
phyloai tree bi --matrix concat/matrix.phy --model lg --mixture 1 --nsamples 10000
phyloai tree bi --output-dir runs/tree/bi --resume
phyloai tree bi --output-dir runs/tree/bi --resume chain1,chain3
```

## Monitoring

`--poll-interval` (default 60 s) controls how often trace files are read to update progress. `--monitor-freq` controls how many new samples must accumulate before the next `bpcomp`+`tracecomp` convergence check.

## Safe Stopping

Use Ctrl+C. PhyloAI writes `0` to each `chains/<chain>.run` file and waits for pb_mpi to finish its current cycle. Direct interruption of pb_mpi can leave incomplete samples.

## Outputs

- `chains/<chain>.trace`: MCMC trace.
- `chains/<chain>.treelist`: sampled trees.
- `chains/<chain>.chain`: saved chain state.
- `chains/<chain>.log`: merged stdout and stderr.
- `convergence/bpcomp_all.bpdiff`: bpcomp summary parsed by PhyloAI.
- `convergence/tracecomp_all.contdiff`: tracecomp summary.
- `convergence/trace_plots.pdf`: trace plots when matplotlib is available.
- `run_state.json`: resume metadata.
- `result.json`: structured PhyloAI result.
````

- [ ] **Step 2: Update stale design references**

Replace `phyloai tree bi phylobayes --matrix ./concat/matrix.phy --chains 3 --threads 8` with `phyloai tree bi --matrix ./concat/matrix.phy --chains 3 --threads 8`.

Replace `runs/tree/bi/phylobayes/` with `runs/tree/bi/`.

In the tree design hierarchy, replace:

```text
├── bi (click.Group)
│   └── phylobayes
```

with:

```text
├── bi                        # Bayesian inference (PhyloBayes-MPI, direct command)
```

- [ ] **Step 3: Verify docs references**

Run: `rg "tree bi phylobayes|bi/phylobayes" docs/superpowers/specs docs/commands`
Expected: no matches.

---

### Task 11: Full Verification

**Files:**
- No edits expected.

- [ ] **Step 1: Run focused BI tests**

Run: `pytest tests/tree/test_bi.py tests/cli/test_tree_bi.py -q`
Expected: all tests pass.

- [ ] **Step 2: Run affected CLI tests**

Run: `pytest tests/cli/test_tree.py tests/cli/test_doctor.py tests/cli/test_tree_bi.py -q`
Expected: all tests pass.

- [ ] **Step 3: Run all tree tests**

Run: `pytest tests/tree tests/cli/test_tree*.py -q`
Expected: all tests pass.

- [ ] **Step 4: Manual dry-run smoke test**

Run: `python -m phyloai.cli.main tree bi --matrix /tmp/nonexistent.phy --dry-run`
Expected: exits 1 with `--matrix does not exist`.

  If no files changed, do not commit.

---

## Self-Review

- Spec coverage: CLI, model flags, chain naming, sampling, resume, run_state, pb-path, convergence commands, `.bpdiff` parsing, tracecomp parsing, plotting, output structure, result.json, doctor integration, stale references, **FASTA-to-PHYLIP auto-conversion**, **tool version detection**, **Rich Live terminal display**, **Ctrl+C soft-stop with final convergence check**, **convergence tool return-code checking**, and **POSIX-portable fake test scripts** are covered.
- Placeholder scan: no unfinished placeholder markers or unspecified implementation steps are present.
- Type consistency: `run_bi()` params match the design signature except `matrix` is `Path | None` to support bare `--resume`; result params still record `matrix` as a nullable field.
- Testing strategy: pure helpers are tested without external tools; process behavior is tested with fake `mpirun`, `pb_mpi`, `bpcomp`, and `tracecomp` scripts.
- Directory safety: conflict check runs BEFORE any filesystem writes, preventing false-positive non-empty-directory errors.
- FASTA handling: `_prepare_matrix()` detects input format and converts FASTA→PHYLIP transparently via `core/formats.py`. `pb_mpi` only sees PHYLIP.
- Version detection: `_detect_pb_version()` uses file-based heuristic; `_detect_mpirun_version()` parses `--version`; results stored in `tool_versions` dict and written to `result.json`.
- Display: `_run_bi_processes` uses `rich.live.Live` + `rich.progress.Progress` with 60s refresh; convergence stats appended as ASCII-only text below progress bars.
- Soft-stop: Ctrl+C writes `.run` files, waits for all subprocesses, runs final convergence check, writes `result.json` with `status: "success"` and `interrupted: true`. No re-raise.
