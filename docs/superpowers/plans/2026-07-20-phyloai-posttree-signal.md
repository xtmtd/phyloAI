# Posttree Signal Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Implement `phyloai posttree signal lnl`, `consistent`, and `fclm` subcommands for phylogenetic signal distribution analysis.

**Architecture:** Three subcommands share a common `posttree/signal.py` library module; CLI wiring lives in `cli/commands/posttree.py`; shared helpers (_parse_sitelh, _sum_gene_lnl, _parse_partition_ranges, _compare_groups) are extracted into signal.py for reuse across lnl and consistent. Each subcommand is independently testable.

**Tech Stack:** Python 3.11+, click, biopython (Bio.Phylo), matplotlib, scipy (mannwhitneyu), concurrent.futures, IQ-TREE3, wASTRAL.

## Global Constraints

- All result.json must conform to JSON output standard (`2026-06-21-phyloai-json-output-standard.md`): Single Pattern with `status`, `command`, `wall_time`, `tool_versions`, `params`, `key_results`, `error`, `data.cmd`, `data.tool_stderr`, `data.tool_log`, `data.output_files`
- `params` dict must include EVERY parameter from `run_*` function signature; null for unused mutually-exclusive params
- `command` field must be full re-executable CLI string with all resolved defaults
- Reuse `phyloai.core.iqtree._resolve_iqtree_path`, `_detect_iqtree_version`, `IQTREE_COMPATIBLE_EXTENSIONS`
- Reuse `phyloai.core.file_matching.scan_tree_dir`, `logical_tree_locus_candidates` for logical locus name resolution (suffix-agnostic, 1-2 dot segment removal)
- Reuse `phyloai.core.formats.FormatConverter().read()` to get matrix taxon IDs; do not implement FASTA-only parsing because FcLM supports IQ-TREE-compatible FASTA, PHYLIP, NEXUS, and PAML inputs
- Reuse `phyloai.pretree.concat._parse_partitions` regex pattern for partition file parsing (RAxML-like: `MODEL, LOCUS = START-END`)
- Reuse `phyloai.pretree.filter._write_outlier_diagnostics` pattern for group comparison plots
- `ambiguous` tie-break threshold: `abs(score_a - score_b) < 1e-9`
- `--threads` passed as-is to IQ-TREE `-T`; for wastral parallelism: int → that value, `"auto"` → `os.cpu_count()`
- No new pip dependencies
- Tests live in `tests/posttree/test_signal_lnl.py`, `test_signal_consistent.py`, `test_signal_fclm.py`
- Spec: `docs/superpowers/specs/2026-07-20-phyloai-posttree-signal-design.md`

---

## File Map

| Action | File | Responsibility |
|--------|------|----------------|
| Create | `phyloai/posttree/signal.py` | All library logic: helpers, run_signal_lnl, run_signal_consistent, run_signal_fclm |
| Modify | `phyloai/cli/commands/posttree.py` | Add signal Click group + 3 subcommand wrappers |
| Modify | `phyloai/report/collector.py` | Add signal.lnl/consistent/fclm to STEP_ORDER and _THIRD_LEVEL |
| Modify | `phyloai/report/templates.py` | Replace posttree.signal stub with 3 generators; update registry |
| Create | `tests/posttree/test_signal_lnl.py` | Tests for lnl subcommand |
| Create | `tests/posttree/test_signal_consistent.py` | Tests for consistent subcommand |
| Create | `tests/posttree/test_signal_fclm.py` | Tests for fclm subcommand |
| Modify | `phyloai/mcp/tools/stubs.py` | Remove `posttree_signal` stub after Click subcommands auto-generate MCP tools |
| Modify | `tests/mcp/test_cli_tools.py` | Verify three generated signal MCP tools replace the stub |
| Modify | `skills/phyloai-workflow/SKILL.md` | Add signal subcommands usage, inputs/outputs, parameter rules |
| Modify | `README.md` | Replace signal stub with 3 example lines |
| Create | `docs/commands/posttree-signal.md` | English command documentation |
| Create | `docs/commands/posttree-signal.zh.md` | Chinese command documentation |

---

## Task 1: Shared helpers in `signal.py`

These helpers are used by both `lnl` and `consistent`. Build and test them first.

**Files:**
- Create: `phyloai/posttree/signal.py`
- Create: `tests/posttree/test_signal_lnl.py` (helpers section only)

**Interfaces produced:**
- `_parse_partition_ranges(path: Path) -> list[dict]` — returns `[{locus, start, end}, ...]` (1-based, inclusive)
- `_parse_sitelh(path: Path) -> tuple[list[str], list[list[float]]]` — returns `(tree_labels, site_scores)` where `site_scores[i]` is list of per-site lnL for tree i
- `_sum_gene_lnl(site_scores: list[list[float]], start: int, end: int) -> list[float]` — returns per-tree gene lnL sums (1-based start/end inclusive)
- `_delta_score(scores: list[float]) -> float` — ΔSLS/ΔGLS: if len==2 returns scores[0]-scores[1]; else mean of all pairwise |a-b|
- `_support_label(scores: list[float], tol: float = 1e-9) -> str` — returns "T1"/"T2"/... or "ambiguous"
- `_outlier_loci(gene_delta: list[float]) -> list[bool]` — True where |ΔGLS| > Q3+1.5*IQR or < Q1-1.5*IQR
- `_compare_groups(group_a: list[str], group_b: list[str], metrics_csv: Path, label_a: str, label_b: str, output_dir: Path) -> tuple[Path, Path]` — writes `<label_a>_comparison.csv` and `<label_a>_comparison.pdf`; returns (csv_path, pdf_path)

- [ ] **Step 1: Create `phyloai/posttree/signal.py` with module docstring and imports**

```python
"""Phylogenetic signal distribution analysis: lnl, consistent, fclm subcommands."""
from __future__ import annotations

import csv
import os
import re
import shlex
import subprocess
import time as _time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Any

import numpy as np
from scipy.stats import mannwhitneyu

from phyloai.core.iqtree import (
    IQTREE_COMPATIBLE_EXTENSIONS,
    _detect_iqtree_version,
    _resolve_iqtree_path,
)
from phyloai.core.file_matching import scan_tree_dir
from phyloai.core.formats import FormatConverter

_PARTITION_RE = re.compile(r"^\s*([^,]+)\s*,\s*(.+?)\s*=\s*(\d+)\s*-\s*(\d+)\s*$")
_FLOAT_TOL = 1e-9
```

- [ ] **Step 2: Implement `_parse_partition_ranges`**

```python
def _parse_partition_ranges(path: Path) -> list[dict[str, Any]]:
    """Parse RAxML-like partition file → [{locus, start, end}, ...] (1-based inclusive)."""
    records: list[dict[str, Any]] = []
    with open(path) as fh:
        for lineno, raw in enumerate(fh, 1):
            line = raw.strip()
            if not line:
                continue
            m = _PARTITION_RE.match(line)
            if m is None:
                raise ValueError(f"Unparseable partition line {lineno}: {raw.rstrip()}")
            _, locus, s, e = m.groups()
            start, end = int(s), int(e)
            if start < 1 or end < start:
                raise ValueError(f"Invalid range line {lineno}: {start}-{end}")
            records.append({"locus": locus.strip(), "start": start, "end": end})
    if not records:
        raise ValueError(f"Partition file is empty: {path}")
    return records
```

- [ ] **Step 3: Implement `_parse_sitelh`**

```python
def _parse_sitelh(path: Path) -> tuple[list[str], list[list[float]]]:
    """Parse IQ-TREE .sitelh file.

    Returns (tree_labels, site_scores) where site_scores[tree_idx] = [lnL_s1, lnL_s2, ...].
    Line 1: '<n_trees> <n_sites>'
    Lines 2+: 'TreeN  lnL1 lnL2 ...'
    """
    with open(path) as fh:
        lines = [l for l in fh if l.strip()]
    header = lines[0].split()
    n_trees, n_sites = int(header[0]), int(header[1])
    tree_labels: list[str] = []
    site_scores: list[list[float]] = []
    for line in lines[1:]:
        parts = line.split()
        if len(parts) != n_sites + 1:
            raise ValueError(f"Expected {n_sites+1} columns in sitelh line, got {len(parts)}")
        tree_labels.append(parts[0])
        site_scores.append([float(x) for x in parts[1:]])
    if len(tree_labels) != n_trees:
        raise ValueError(f"Expected {n_trees} tree rows, got {len(tree_labels)}")
    return tree_labels, site_scores
```

- [ ] **Step 4: Implement scoring helpers**

```python
def _sum_gene_lnl(site_scores: list[list[float]], start: int, end: int) -> list[float]:
    """Sum site lnL values for a locus (1-based inclusive start/end) per tree."""
    return [sum(scores[start - 1:end]) for scores in site_scores]


def _delta_score(scores: list[float]) -> float:
    """ΔSLS/ΔGLS: signed diff for 2 trees; mean pairwise |diff| for >2."""
    if len(scores) == 2:
        return scores[0] - scores[1]
    total, count = 0.0, 0
    for i in range(len(scores)):
        for j in range(i + 1, len(scores)):
            total += abs(scores[i] - scores[j])
            count += 1
    return total / count if count else 0.0


def _support_label(scores: list[float], labels: list[str]) -> str:
    """Return label of best-scoring tree, or 'ambiguous' if tied."""
    best_idx = int(np.argmax(scores))
    best_val = scores[best_idx]
    ties = [i for i, s in enumerate(scores) if abs(s - best_val) < _FLOAT_TOL]
    if len(ties) > 1:
        return "ambiguous"
    return labels[best_idx]


def _outlier_loci(gene_deltas: list[float]) -> list[bool]:
    """Return bool mask: True where |ΔGLS| is outside Tukey 1.5*IQR whiskers."""
    arr = np.array([abs(d) for d in gene_deltas])
    q1, q3 = float(np.percentile(arr, 25)), float(np.percentile(arr, 75))
    iqr = q3 - q1
    lo, hi = q1 - 1.5 * iqr, q3 + 1.5 * iqr
    return [bool(v < lo or v > hi) for v in arr]
```

- [ ] **Step 5: Implement `_compare_groups`**

```python
def _compare_groups(
    group_a_loci: list[str],
    group_b_loci: list[str],
    metrics_csv: Path,
    label_a: str,
    label_b: str,
    output_dir: Path,
    csv_filename: str,
    pdf_filename: str,
) -> tuple[Path, Path]:
    """Compare two locus groups across metrics from a pretree metrics CSV.

    Validates that all loci in both groups are present in the metrics CSV.
    When group_a is empty, writes an NA-filled CSV and a blank PDF (no error).
    Writes <csv_filename> (means + Wilcoxon p-value per metric) and
    <pdf_filename> (boxplots). Returns (csv_path, pdf_path).
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    # Read metrics CSV
    with open(metrics_csv, newline="") as fh:
        reader = csv.DictReader(fh)
        rows = list(reader)
        fieldnames = reader.fieldnames or []

    locus_col = "loci" if "loci" in fieldnames else fieldnames[0]
    metrics_loci = {r[locus_col] for r in rows}

    # Validate all analysis loci present in metrics file
    all_analysis_loci = set(group_a_loci) | set(group_b_loci)
    missing = all_analysis_loci - metrics_loci
    if missing:
        raise ValueError(
            f"Loci missing from --metrics file: {', '.join(sorted(missing))}"
        )

    numeric_cols = [
        c for c in fieldnames
        if c != locus_col and all(
            r.get(c, "") not in ("", "NA") and _is_numeric(r.get(c, ""))
            for r in rows if r.get(locus_col) in all_analysis_loci
        )
    ]

    set_a = set(group_a_loci)
    set_b = set(group_b_loci)
    rows_a = [r for r in rows if r.get(locus_col) in set_a]
    rows_b = [r for r in rows if r.get(locus_col) in set_b]

    # Handle empty group_a: write NA CSV and blank PDF
    csv_path = output_dir / csv_filename
    pdf_path = output_dir / pdf_filename
    if not group_a_loci:
        with open(csv_path, "w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=["metric", f"{label_a}_mean", f"{label_a}_n",
                                                      f"{label_b}_mean", f"{label_b}_n", "wilcoxon_p"])
            writer.writeheader()
            for col in numeric_cols:
                vals_b = [float(r[col]) for r in rows_b if _is_numeric(r.get(col, ""))]
                writer.writerow({"metric": col, f"{label_a}_mean": "NA", f"{label_a}_n": 0,
                                  f"{label_b}_mean": round(float(np.mean(vals_b)), 6) if vals_b else "NA",
                                  f"{label_b}_n": len(vals_b), "wilcoxon_p": "NA"})
        fig, ax = plt.subplots(figsize=(4, 3))
        ax.text(0.5, 0.5, f"No {label_a} loci", ha="center", va="center", transform=ax.transAxes)
        fig.savefig(pdf_path, dpi=150, bbox_inches="tight")
        plt.close(fig)
        return csv_path, pdf_path

    comparison_rows = []
    for col in numeric_cols:
        vals_a = [float(r[col]) for r in rows_a if _is_numeric(r.get(col, ""))]
        vals_b = [float(r[col]) for r in rows_b if _is_numeric(r.get(col, ""))]
        row: dict[str, Any] = {
            "metric": col,
            f"{label_a}_mean": round(float(np.mean(vals_a)), 6) if vals_a else "NA",
            f"{label_a}_n": len(vals_a),
            f"{label_b}_mean": round(float(np.mean(vals_b)), 6) if vals_b else "NA",
            f"{label_b}_n": len(vals_b),
            "wilcoxon_p": "NA",
        }
        if vals_a and vals_b:
            stat = mannwhitneyu(vals_a, vals_b, alternative="two-sided")
            row["wilcoxon_p"] = round(float(stat.pvalue), 6)
        comparison_rows.append(row)

    # Write CSV
    csv_path = output_dir / csv_filename
    fieldnames_out = ["metric", f"{label_a}_mean", f"{label_a}_n", f"{label_b}_mean", f"{label_b}_n", "wilcoxon_p"]
    with open(csv_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames_out)
        writer.writeheader()
        writer.writerows(comparison_rows)

    # Write PDF boxplots
    pdf_path = output_dir / pdf_filename
    n_cols = 4
    n_rows = max(1, -(-len(numeric_cols) // n_cols))
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(n_cols * 3.5, n_rows * 2.8))
    axes_flat = np.array(axes).flatten()
    for idx, col in enumerate(numeric_cols):
        ax = axes_flat[idx]
        vals_a = [float(r[col]) for r in rows_a if _is_numeric(r.get(col, ""))]
        vals_b = [float(r[col]) for r in rows_b if _is_numeric(r.get(col, ""))]
        bp = ax.boxplot([vals_a, vals_b], tick_labels=[label_a, label_b], patch_artist=True)
        bp["boxes"][0].set_facecolor("#a6cee3"); bp["boxes"][0].set_alpha(0.6)
        bp["boxes"][1].set_facecolor("#fb9a99"); bp["boxes"][1].set_alpha(0.6)
        ax.set_title(col, fontsize=9)
        ax.tick_params(axis="x", labelsize=7)
    for idx in range(len(numeric_cols), len(axes_flat)):
        axes_flat[idx].set_visible(False)
    fig.tight_layout(pad=2.0)
    fig.savefig(pdf_path, dpi=150, bbox_inches="tight")
    plt.close(fig)

    return csv_path, pdf_path


def _is_numeric(val: Any) -> bool:
    try:
        float(val)
        return True
    except (TypeError, ValueError):
        return False
```

- [ ] **Step 6: Write tests for helpers**

```python
# tests/posttree/test_signal_lnl.py
"""Tests for phyloai.posttree.signal helpers."""
from __future__ import annotations
from pathlib import Path
import pytest
from phyloai.posttree.signal import (
    _parse_partition_ranges,
    _parse_sitelh,
    _sum_gene_lnl,
    _delta_score,
    _support_label,
    _outlier_loci,
)

class TestParsePartitionRanges:
    def test_basic(self, tmp_path: Path) -> None:
        f = tmp_path / "p.txt"
        f.write_text("LG, geneA = 1-235\nLG, geneB = 236-461\n")
        recs = _parse_partition_ranges(f)
        assert len(recs) == 2
        assert recs[0] == {"locus": "geneA", "start": 1, "end": 235}
        assert recs[1] == {"locus": "geneB", "start": 236, "end": 461}

    def test_empty_raises(self, tmp_path: Path) -> None:
        f = tmp_path / "empty.txt"
        f.write_text("")
        with pytest.raises(ValueError, match="empty"):
            _parse_partition_ranges(f)

    def test_bad_line_raises(self, tmp_path: Path) -> None:
        f = tmp_path / "bad.txt"
        f.write_text("not a partition line\n")
        with pytest.raises(ValueError, match="Unparseable"):
            _parse_partition_ranges(f)


class TestParseSitelh:
    def test_two_trees(self, tmp_path: Path) -> None:
        f = tmp_path / "site.sitelh"
        f.write_text("2 3\nTree1 -1.0 -2.0 -3.0\nTree2 -1.5 -2.5 -3.5\n")
        labels, scores = _parse_sitelh(f)
        assert labels == ["Tree1", "Tree2"]
        assert scores[0] == pytest.approx([-1.0, -2.0, -3.0])
        assert scores[1] == pytest.approx([-1.5, -2.5, -3.5])

    def test_wrong_column_count_raises(self, tmp_path: Path) -> None:
        f = tmp_path / "bad.sitelh"
        f.write_text("2 3\nTree1 -1.0 -2.0\nTree2 -1.5 -2.5 -3.5\n")
        with pytest.raises(ValueError, match="columns"):
            _parse_sitelh(f)


class TestDeltaScore:
    def test_two_trees_signed(self) -> None:
        assert _delta_score([-5.0, -8.0]) == pytest.approx(3.0)

    def test_three_trees_mean_pairwise(self) -> None:
        # pairs: |a-b|=1, |a-c|=2, |b-c|=1 → mean=4/3
        result = _delta_score([-1.0, -2.0, -3.0])
        assert result == pytest.approx(4/3)


class TestSupportLabel:
    def test_clear_winner(self) -> None:
        assert _support_label([-5.0, -8.0], ["T1", "T2"]) == "T1"

    def test_tie_returns_ambiguous(self) -> None:
        assert _support_label([-5.0, -5.0], ["T1", "T2"]) == "ambiguous"

    def test_tie_within_tolerance(self) -> None:
        assert _support_label([-5.0, -5.0 + 5e-10], ["T1", "T2"]) == "ambiguous"


class TestOutlierLoci:
    def test_no_outliers_all_same(self) -> None:
        deltas = [1.0] * 10
        assert not any(_outlier_loci(deltas))

    def test_extreme_value_is_outlier(self) -> None:
        deltas = [1.0] * 9 + [100.0]
        flags = _outlier_loci(deltas)
        assert flags[-1] is True
        assert not any(flags[:-1])
```

- [ ] **Step 7: Run helper tests**

```bash
cd /Users/zf/data/coding/phyloAI-syserror
python -m pytest tests/posttree/test_signal_lnl.py -v -k "TestParse or TestDelta or TestSupport or TestOutlier"
```

Expected: all pass.

---

## Task 2: `signal lnl` — library function `run_signal_lnl`

**Files:**
- Modify: `phyloai/posttree/signal.py`
- Modify: `tests/posttree/test_signal_lnl.py`

**Interfaces consumed:** helpers from Task 1, `_resolve_iqtree_path`, `_detect_iqtree_version`

**Interfaces produced:**
- `run_signal_lnl(*, matrix, candidate_trees, model_expr, partitions, locus_ranges, guide_tree, threads, iqtree_path, tool_args, metrics, output_dir, overwrite, dry_run, quiet) -> dict`

- [ ] **Step 1: Write failing test for dry-run**

```python
# Add to tests/posttree/test_signal_lnl.py
class TestRunSignalLnlDryRun:
    def test_dry_run_returns_cmd(self, tmp_path: Path) -> None:
        from phyloai.posttree.signal import run_signal_lnl
        matrix = tmp_path / "m.fa"
        matrix.write_text(">A\nMKT\n>B\nMKA\n")
        t1 = tmp_path / "T1.nwk"
        t1.write_text("(A,B);\n")
        t2 = tmp_path / "T2.nwk"
        t2.write_text("(B,A);\n")
        result = run_signal_lnl(
            matrix=matrix,
            candidate_trees=[t1, t2],
            model_expr="LG+F+R4",
            partitions=None, locus_ranges=None,
            guide_tree=None, threads="auto",
            iqtree_path=None, tool_args=None,
            metrics=None,
            output_dir=tmp_path / "out",
            overwrite=False, dry_run=True, quiet=True,
        )
        assert result["status"] == "success"
        assert "iqtree3" in result["data"]["cmd"][0]
        assert "-wslr" in result["data"]["cmd"]
        assert result["params"]["model_expr"] == "LG+F+R4"
        assert result["params"]["partitions"] is None
        assert result["params"]["locus_ranges"] is None
```

- [ ] **Step 2: Run to verify failure**

```bash
python -m pytest tests/posttree/test_signal_lnl.py::TestRunSignalLnlDryRun -v
```

Expected: ImportError or AttributeError (function not yet defined).

- [ ] **Step 3: Implement `run_signal_lnl`**

Add to `phyloai/posttree/signal.py`:

```python
_LNL_BLOCKED_FLAGS = frozenset({"-s", "-z", "-wslr"})


def _validate_lnl_inputs(
    *,
    matrix: Path,
    candidate_trees: list[Path],
    model_expr: str | None,
    partitions: Path | None,
    locus_ranges: Path | None,
    tool_args: str | None,
    guide_tree: Path | None,
) -> list[str]:
    errors: list[str] = []
    import os as _os
    if not matrix.exists() or not matrix.is_file():
        errors.append(f"--matrix does not exist: {matrix}")
    for i, ct in enumerate(candidate_trees):
        if not ct.exists() or not ct.is_file():
            errors.append(f"--candidate-trees #{i+1} does not exist: {ct}")
        elif ct.stat().st_size == 0:
            errors.append(f"--candidate-trees #{i+1} is empty: {ct}")
    if partitions and locus_ranges:
        errors.append("--partitions and --locus-ranges are mutually exclusive")
    if model_expr and partitions:
        errors.append("--model-expr and --partitions are mutually exclusive")
    has_model = model_expr or partitions
    has_tool_model = False
    if tool_args:
        toks = set(shlex.split(tool_args))
        has_tool_model = "-m" in toks or "-p" in toks
        for flag in _LNL_BLOCKED_FLAGS:
            if flag in toks:
                errors.append(f"Blocked flag in --tool-args: {flag}")
    if not has_model and not has_tool_model:
        errors.append("Must specify --model-expr, --partitions, or -m/-p in --tool-args")
    if guide_tree and not Path(guide_tree).exists():
        errors.append(f"--guide-tree does not exist: {guide_tree}")
    return errors


def _build_lnl_cmd(
    *,
    executable: str,
    matrix: Path,
    candidate_trees: Path,
    prefix: str,
    model_expr: str | None,
    partitions: str | None,
    guide_tree: str | None,
    threads: str,
    tool_args: str | None,
) -> list[str]:
    cmd = [executable, "-s", str(matrix), "-z", str(candidate_trees)]
    tool_toks = set(shlex.split(tool_args)) if tool_args else set()

    if "--prefix" not in tool_toks:
        cmd.extend(["--prefix", prefix])
    if model_expr and "-m" not in tool_toks:
        cmd.extend(["-m", model_expr])
    elif partitions and "-p" not in tool_toks:
        cmd.extend(["-p", partitions])
    if guide_tree and "-ft" not in tool_toks:
        cmd.extend(["-ft", guide_tree])
    cmd.append("-wslr")
    if "-T" not in tool_toks:
        cmd.extend(["-T", str(threads)])
    if tool_args:
        cmd.extend(shlex.split(tool_args))
    return cmd


def run_signal_lnl(
    *,
    matrix: Path,
    candidate_trees: list[Path],
    model_expr: str | None = None,
    partitions: Path | None = None,
    partition_mode: str = "p",
    locus_ranges: Path | None = None,
    guide_tree: Path | None = None,
    threads: str = "auto",
    iqtree_path: str | None = None,
    tool_args: str | None = None,
    metrics: Path | None = None,
    output_dir: Path | None = None,
    overwrite: bool = False,
    dry_run: bool = False,
    quiet: bool = False,
) -> dict[str, Any]:
    run_start = _time.time()
    if output_dir is None:
        output_dir = Path("runs/posttree/signal/lnl")
    output_dir = output_dir.resolve()
    matrix = matrix.resolve()
    candidate_trees = [ct.resolve() for ct in candidate_trees]

    errors = _validate_lnl_inputs(
        matrix=matrix, candidate_trees=candidate_trees,
        model_expr=model_expr, partitions=partitions,
        locus_ranges=locus_ranges, tool_args=tool_args,
        guide_tree=guide_tree,
    )

    params: dict[str, Any] = {
        "matrix": str(matrix),
        "candidate_trees_raw": ",".join(str(ct) for ct in candidate_trees),
        "model_expr": model_expr,
        "partitions": str(partitions.resolve()) if partitions else None,
        "locus_ranges": str(locus_ranges.resolve()) if locus_ranges else None,
        "guide_tree": str(guide_tree.resolve()) if guide_tree else None,
        "threads": threads,
        "iqtree_path": iqtree_path,
        "tool_args": tool_args,
        "metrics": str(metrics.resolve()) if metrics else None,
        "output_dir": str(output_dir),
        "overwrite": overwrite,
        "dry_run": dry_run,
        "quiet": quiet,
    }

    if errors:
        return {"status": "error", "command": "", "wall_time": 0.0,
                "tool_versions": {}, "params": params, "key_results": {},
                "error": "; ".join(errors), "error_category": "input",
                "data": {"cmd": [], "tool_stderr": "", "tool_log": None, "output_files": {}}}

    if not dry_run:
        if overwrite and output_dir.exists():
            import shutil; shutil.rmtree(output_dir)
        elif output_dir.exists() and any(output_dir.iterdir()):
            return {"status": "error", "command": "", "wall_time": 0.0,
                    "tool_versions": {}, "params": params, "key_results": {},
                    "error": f"Output directory '{output_dir}' already exists and is non-empty. Use --overwrite to replace it.",
                    "error_category": "input",
                    "data": {"cmd": [], "tool_stderr": "", "tool_log": None, "output_files": {}}}
        output_dir.mkdir(parents=True, exist_ok=True)

    # Merge candidate trees if multiple files
    if len(candidate_trees) == 1:
        candidate_trees_path = candidate_trees[0]
    else:
        candidate_trees_path = output_dir / "candidate.trees"
        if not dry_run:
            with open(candidate_trees_path, "w") as fh:
                for ct in candidate_trees:
                    text = ct.read_text().strip()
                    fh.write(text + "\n")

    prefix = "lnl"
    resolved_partitions = str(partitions.resolve()) if partitions else None
    resolved_guide_tree = str(guide_tree.resolve()) if guide_tree else None

    try:
        iqtree_exe = _resolve_iqtree_path(iqtree_path, dry_run)
    except (ValueError, FileNotFoundError) as e:
        return {"status": "error", "command": "", "wall_time": 0.0,
                "tool_versions": {}, "params": params, "key_results": {},
                "error": str(e), "error_category": "env",
                "data": {"cmd": [], "tool_stderr": "", "tool_log": None, "output_files": {}}}

    tool_versions = _detect_iqtree_version(iqtree_exe) if not dry_run else {"iqtree3": "dry-run"}

    cmd = _build_lnl_cmd(
        executable=iqtree_exe, matrix=matrix,
        candidate_trees=candidate_trees_path,
        prefix=prefix, model_expr=model_expr,
        partitions=resolved_partitions,
        guide_tree=resolved_guide_tree,
        threads=threads, tool_args=tool_args,
    )

    # Build full CLI command string for result.json
    cli_parts = ["phyloai", "posttree", "signal", "lnl",
                 "--matrix", str(matrix),
                 "--candidate-trees", params["candidate_trees_raw"]]
    if model_expr: cli_parts.extend(["--model-expr", model_expr])
    if partitions: cli_parts.extend(["--partitions", str(partitions)])
    if locus_ranges: cli_parts.extend(["--locus-ranges", str(locus_ranges)])
    if guide_tree: cli_parts.extend(["--guide-tree", str(guide_tree)])
    if tool_args: cli_parts.extend(["--tool-args", tool_args])
    if metrics: cli_parts.extend(["--metrics", str(metrics)])
    cli_parts.extend(["--threads", str(threads), "-o", str(output_dir)])
    if iqtree_path: cli_parts.extend(["--iqtree-path", iqtree_path])
    if overwrite: cli_parts.append("--overwrite")
    if dry_run: cli_parts.append("--dry-run")
    if quiet: cli_parts.append("-q")
    full_command = shlex.join(cli_parts)

    if dry_run:
        return {"status": "success", "command": full_command,
                "wall_time": 0.0, "tool_versions": tool_versions,
                "params": params, "key_results": {}, "error": None,
                "data": {"cmd": cmd, "tool_stderr": "", "tool_log": None,
                         "output_files": {}}}

    # Run IQ-TREE
    proc = subprocess.run(cmd, capture_output=True, text=True, cwd=str(output_dir))
    tool_stderr = (proc.stdout or "") + ("\n" if proc.stdout and proc.stderr else "") + (proc.stderr or "")
    tool_log_path = output_dir / f"{prefix}.log"

    output_files: dict[str, Any] = {
        "iqtree_report": {"path": str(output_dir / f"{prefix}.iqtree"),
                          "description": "IQ-TREE native report"},
        "iqtree_sitelh": {"path": str(output_dir / f"{prefix}.sitelh"),
                          "description": "IQ-TREE raw site log-likelihoods"},
    }

    if proc.returncode != 0:
        wall_time = _time.time() - run_start
        _write_result_json(output_dir, {
            "status": "error", "command": full_command, "wall_time": wall_time,
            "tool_versions": tool_versions, "params": params, "key_results": {},
            "error": f"IQ-TREE exited with code {proc.returncode}",
            "error_category": "tool",
            "data": {"cmd": cmd, "tool_stderr": tool_stderr,
                     "tool_log": str(tool_log_path), "output_files": output_files},
        })
        return _load_result_json(output_dir)

    # Parse .sitelh and generate outputs
    sitelh_path = output_dir / f"{prefix}.sitelh"
    tree_labels, site_scores = _parse_sitelh(sitelh_path)
    n_trees = len(tree_labels)
    n_sites = len(site_scores[0]) if site_scores else 0

    # Build site_lnl.csv
    site_csv_path = output_dir / "site_lnl.csv"
    site_rows = []
    for si in range(n_sites):
        scores = [site_scores[ti][si] for ti in range(n_trees)]
        delta = _delta_score(scores)
        support = _support_label(scores, tree_labels)
        row: dict[str, Any] = {"site": si + 1}
        for ti, lbl in enumerate(tree_labels):
            row[f"lnL_{lbl}"] = round(scores[ti], 6)
        row["ΔSLS"] = round(delta, 6)
        row["support"] = support
        site_rows.append(row)
    site_rows.sort(key=lambda r: r["ΔSLS"], reverse=True)
    _write_csv(site_csv_path, site_rows,
               ["site"] + [f"lnL_{l}" for l in tree_labels] + ["ΔSLS", "support"])
    output_files["site_lnl"] = {"path": str(site_csv_path),
                                 "description": "Site-wise lnL scores per tree, ΔSLS, support; sorted by ΔSLS descending"}

    # Site support bar chart
    site_pdf = _plot_support_bar(
        [r["support"] for r in site_rows], tree_labels,
        output_dir / "site_support.pdf",
        xlabel="Supported topology", ylabel="Number of sites",
    )
    output_files["site_support_plot"] = {"path": str(site_pdf),
                                          "description": "Site support distribution bar chart"}

    # Gene-wise (if locus boundaries available)
    boundary_path = partitions or locus_ranges
    n_loci, n_outlier_genes = None, None
    if boundary_path:
        partition_recs = _parse_partition_ranges(boundary_path.resolve())
        gene_rows = []
        for rec in partition_recs:
            gene_scores = _sum_gene_lnl(site_scores, rec["start"], rec["end"])
            delta = _delta_score(gene_scores)
            support = _support_label(gene_scores, tree_labels)
            row = {"locus": rec["locus"]}
            for ti, lbl in enumerate(tree_labels):
                row[f"lnL_{lbl}"] = round(gene_scores[ti], 6)
            row["ΔGLS"] = round(delta, 6)
            row["support"] = support
            if n_trees == 2:
                row["support_sig"] = abs(delta) >= 2.0
            gene_rows.append(row)
        gene_rows.sort(key=lambda r: r["ΔGLS"], reverse=True)

        gene_cols = ["locus"] + [f"lnL_{l}" for l in tree_labels] + ["ΔGLS", "support"]
        if n_trees == 2:
            gene_cols.append("support_sig")
        gene_csv_path = output_dir / "gene_lnl.csv"
        _write_csv(gene_csv_path, gene_rows, gene_cols)
        output_files["gene_lnl"] = {"path": str(gene_csv_path),
                                     "description": "Gene-wise lnL scores per tree, ΔGLS, support; sorted by ΔGLS descending"}

        gene_pdf = _plot_support_bar(
            [r["support"] for r in gene_rows], tree_labels,
            output_dir / "gene_support.pdf",
            xlabel="Supported topology", ylabel="Number of genes",
        )
        output_files["gene_support_plot"] = {"path": str(gene_pdf),
                                              "description": "Gene support distribution bar chart"}

        # Outliers
        deltas = [r["ΔGLS"] for r in gene_rows]
        outlier_flags = _outlier_loci(deltas)
        outlier_loci = [gene_rows[i]["locus"] for i, f in enumerate(outlier_flags) if f]
        outlier_txt = output_dir / "outlier_genes.txt"
        outlier_txt.write_text("\n".join(outlier_loci) + ("\n" if outlier_loci else ""))
        output_files["outlier_genes"] = {"path": str(outlier_txt),
                                          "description": "Loci with |ΔGLS| outside boxplot whiskers (Shen 2017 eq. 3/4)"}
        n_loci = len(gene_rows)
        n_outlier_genes = len(outlier_loci)

        # Metrics comparison (optional)
        if metrics and outlier_loci:
            all_loci = [r["locus"] for r in gene_rows]
            non_outlier = [l for l in all_loci if l not in set(outlier_loci)]
            csv_p, pdf_p = _compare_groups(
                outlier_loci, non_outlier, metrics.resolve(),
                "outlier", "non_outlier",
                output_dir, "outlier_comparison.csv", "outlier_comparison.pdf",
            )
            output_files["outlier_comparison"] = {"path": str(csv_p),
                "description": "Outlier vs non-outlier per-metric means and Wilcoxon p-values"}
            output_files["outlier_comparison_plot"] = {"path": str(pdf_p),
                "description": "Outlier vs non-outlier metric distribution boxplots"}

    key_results: dict[str, Any] = {"n_trees": n_trees, "n_sites": n_sites}
    if n_loci is not None:
        key_results["n_loci"] = n_loci
        key_results["n_outlier_genes"] = n_outlier_genes

    wall_time = _time.time() - run_start
    result = {
        "status": "success", "command": full_command, "wall_time": wall_time,
        "tool_versions": tool_versions, "params": params,
        "key_results": key_results, "error": None,
        "data": {
            "cmd": cmd, "tool_stderr": tool_stderr,
            "tool_log": str(tool_log_path),
            "summary": key_results,
            "output_files": output_files,
        },
    }
    _write_result_json(output_dir, result)
    return result
```

- [ ] **Step 4: Implement small I/O helpers (add to signal.py)**

```python
def _write_csv(path: Path, rows: list[dict], fieldnames: list[str]) -> None:
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def _write_result_json(output_dir: Path, result: dict) -> None:
    import json
    (output_dir / "result.json").write_text(json.dumps(result, indent=2))


def _load_result_json(output_dir: Path) -> dict:
    import json
    return json.loads((output_dir / "result.json").read_text())


def _plot_support_bar(
    support_values: list[str],
    tree_labels: list[str],
    output_path: Path,
    xlabel: str = "Supported topology",
    ylabel: str = "Count",
) -> Path:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    categories = tree_labels + ["ambiguous"]
    counts = {c: 0 for c in categories}
    for s in support_values:
        counts[s] = counts.get(s, 0) + 1

    fig, ax = plt.subplots(figsize=(max(4, len(categories) * 1.5), 5))
    ax.bar(list(counts.keys()), list(counts.values()))
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    fig.tight_layout()
    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return output_path
```

- [ ] **Step 5: Run dry-run test**

```bash
python -m pytest tests/posttree/test_signal_lnl.py::TestRunSignalLnlDryRun -v
```

Expected: PASS.

- [ ] **Step 6: Write integration test using real .sitelh fixture**

```python
class TestRunSignalLnlWithFixture:
    """Uses runs/signal/ test data."""
    SIGNAL_DIR = Path("runs/signal")

    def test_site_lnl_csv_columns(self, tmp_path: Path) -> None:
        from phyloai.posttree.signal import run_signal_lnl
        if not self.SIGNAL_DIR.exists():
            pytest.skip("Signal test data not present")
        matrix = self.SIGNAL_DIR / "matrix.aa.fa"
        trees = self.SIGNAL_DIR / "trees"
        result = run_signal_lnl(
            matrix=matrix, candidate_trees=[trees],
            model_expr="LG+F+R4", partitions=None,
            locus_ranges=self.SIGNAL_DIR / "matrix.aa.partitions",
            guide_tree=None, threads="auto",
            iqtree_path=None, tool_args=None,
            metrics=None,
            output_dir=tmp_path / "lnl_out",
            overwrite=False, dry_run=False, quiet=True,
        )
        if result["status"] == "error" and result.get("error_category") == "env":
            pytest.skip("iqtree3 not available")
        assert result["status"] == "success"
        assert result["key_results"]["n_trees"] == 3
        assert result["key_results"]["n_sites"] == 5604
        assert result["key_results"]["n_loci"] == 20
        import csv as csv_mod
        with open(tmp_path / "lnl_out" / "site_lnl.csv") as fh:
            reader = csv_mod.DictReader(fh)
            cols = reader.fieldnames or []
        assert "ΔSLS" in cols
        assert "support" in cols
        assert "lnL_Tree1" in cols
```

- [ ] **Step 7: Run all lnl tests**

```bash
python -m pytest tests/posttree/test_signal_lnl.py -v
```

Expected: helpers PASS; integration PASS or skip (no iqtree3).

---

## Task 3: `signal consistent` — library function `run_signal_consistent`

**Files:**
- Modify: `phyloai/posttree/signal.py`
- Create: `tests/posttree/test_signal_consistent.py`

**Interfaces consumed:** helpers from Task 1, `scan_tree_dir` from `core.file_matching`, `_resolve_iqtree_path`, wastral detection from `tree.msc._resolve_wastral_path`

**Interfaces produced:**
- `run_signal_consistent(*, matrix, candidate_trees, tree_dir, model_expr, partitions, partition_mode, locus_ranges, guide_tree, threads, iqtree_path, wastral_path, tool_args, metrics, output_dir, overwrite, dry_run, quiet) -> dict`

- [ ] **Step 1: Write failing dry-run test**

```python
# tests/posttree/test_signal_consistent.py
"""Tests for phyloai.posttree.signal consistent subcommand."""
from __future__ import annotations
from pathlib import Path
import pytest


class TestRunSignalConsistentValidation:
    def test_more_than_two_trees_errors(self, tmp_path: Path) -> None:
        from phyloai.posttree.signal import run_signal_consistent
        matrix = tmp_path / "m.fa"
        matrix.write_text(">A\nMKT\n>B\nMKA\n")
        trees = [tmp_path / f"T{i}.nwk" for i in range(3)]
        for t in trees:
            t.write_text("(A,B);\n")
        tree_dir = tmp_path / "gtrees"
        tree_dir.mkdir()
        result = run_signal_consistent(
            matrix=matrix, candidate_trees=trees,
            tree_dir=tree_dir, model_expr="LG+F+R4",
            partitions=None, partition_mode="p",
            locus_ranges=tmp_path / "fake.partitions",
            guide_tree=None, threads="auto",
            iqtree_path=None, wastral_path=None,
            tool_args=None, metrics=None,
            output_dir=tmp_path / "out",
            overwrite=False, dry_run=True, quiet=True,
        )
        assert result["status"] == "error"
        assert "2" in result["error"]

    def test_no_locus_boundaries_errors(self, tmp_path: Path) -> None:
        from phyloai.posttree.signal import run_signal_consistent
        matrix = tmp_path / "m.fa"
        matrix.write_text(">A\nMKT\n>B\nMKA\n")
        trees = [tmp_path / f"T{i}.nwk" for i in range(2)]
        for t in trees:
            t.write_text("(A,B);\n")
        tree_dir = tmp_path / "gtrees"
        tree_dir.mkdir()
        result = run_signal_consistent(
            matrix=matrix, candidate_trees=trees,
            tree_dir=tree_dir, model_expr="LG+F+R4",
            partitions=None, partition_mode="p",
            locus_ranges=None,
            guide_tree=None, threads="auto",
            iqtree_path=None, wastral_path=None,
            tool_args=None, metrics=None,
            output_dir=tmp_path / "out",
            overwrite=False, dry_run=True, quiet=True,
        )
        assert result["status"] == "error"
        assert "locus" in result["error"].lower() or "partition" in result["error"].lower()
```

- [ ] **Step 2: Run to verify failure**

```bash
python -m pytest tests/posttree/test_signal_consistent.py -v
```

Expected: ImportError (function not defined).

- [ ] **Step 3: Implement GQS helper and `run_signal_consistent`**

Add to `phyloai/posttree/signal.py`:

```python
_CONSISTENT_BLOCKED_FLAGS = frozenset({"-s", "-z", "-wslr"})


def _prune_reference_tree(ref_tree_str: str, taxa_to_remove: set[str]) -> str:
    """Prune taxa from a newick tree string using Bio.Phylo. Returns pruned newick."""
    from Bio import Phylo
    from io import StringIO
    tree = Phylo.read(StringIO(ref_tree_str), "newick")
    for taxon in taxa_to_remove:
        try:
            tree.prune(taxon)
        except Exception:
            pass  # taxon not in tree is fine
    out = StringIO()
    Phylo.write(tree, out, "newick")
    return out.getvalue().strip()


def _count_tree_taxa(newick_str: str) -> int:
    from Bio import Phylo
    from io import StringIO
    tree = Phylo.read(StringIO(newick_str), "newick")
    return len(tree.get_terminals())


def _get_tree_taxa(newick_str: str) -> set[str]:
    from Bio import Phylo
    from io import StringIO
    tree = Phylo.read(StringIO(newick_str), "newick")
    return {c.name for c in tree.get_terminals() if c.name}


def _run_wastral_gqs(
    gene_tree_path: Path,
    ref_tree_str: str,
    wastral_exe: str,
    work_dir: Path,
    locus: str,
) -> float:
    """Run wastral quartet score for one gene tree vs one reference tree.

    Returns score float. Raises RuntimeError if wastral exits non-zero
    or Score: line is absent (external-tool failure → caller writes error result).
    """
    ref_path = work_dir / f"ref_{locus}.nwk"
    ref_path.write_text(ref_tree_str)
    proc = subprocess.run(
        [wastral_exe, "-i", str(gene_tree_path), "-C", "-c", str(ref_path), "--mode", "4"],
        capture_output=True, text=True,
    )
    combined = (proc.stdout or "") + "\n" + (proc.stderr or "")
    if proc.returncode != 0:
        raise RuntimeError(
            f"wastral exited {proc.returncode} for locus {locus!r}.\n{combined[:500]}"
        )
    for line in combined.splitlines():
        if line.strip().startswith("Score:"):
            try:
                return float(line.split(":")[1].strip())
            except ValueError:
                pass
    raise RuntimeError(
        f"wastral produced no 'Score:' line for locus {locus!r}.\n{combined[:500]}"
    )


def _compute_gqs_for_locus(
    gene_tree_path: Path,
    logical_locus: str,
    t1_str: str,
    t2_str: str,
    ref_taxa: set[str],
    wastral_exe: str,
    work_dir: Path,
) -> dict[str, Any]:
    """Compute GQS for a single gene tree. Returns gqs row dict.

    Raises ValueError for input errors (extra taxa → hard error).
    Raises RuntimeError for wastral tool failures → caller writes error result.
    """
    # Read gene tree taxa
    gene_str = gene_tree_path.read_text()
    gene_taxa = _get_tree_taxa(gene_str)

    # Check for extra taxa (gene tree has taxa not in reference) → hard error
    extra = gene_taxa - ref_taxa
    if extra:
        raise ValueError(
            f"Gene tree {gene_tree_path.name} contains taxa not in reference trees: "
            f"{', '.join(sorted(extra))}. Matrix and gene trees must share the same taxon set."
        )

    missing = ref_taxa - gene_taxa
    t1_pruned = _prune_reference_tree(t1_str, missing)
    t2_pruned = _prune_reference_tree(t2_str, missing)

    if _count_tree_taxa(t1_pruned) < 4:
        return {
            "locus": logical_locus, "GQS_T1": None, "GQS_T2": None,
            "ΔGQS": None, "support": "ambiguous",
            "status": "skipped", "reason": "pruned_tree_too_small",
        }

    locus_work = work_dir / logical_locus
    locus_work.mkdir(parents=True, exist_ok=True)
    # RuntimeError propagates to caller → error result written at run_signal_consistent level
    gqs_t1 = _run_wastral_gqs(gene_tree_path, t1_pruned, wastral_exe, locus_work, logical_locus)
    gqs_t2 = _run_wastral_gqs(gene_tree_path, t2_pruned, wastral_exe, locus_work, logical_locus)

    delta = gqs_t1 - gqs_t2
    if abs(delta) < _FLOAT_TOL:
        support = "ambiguous"
    elif delta > 0:
        support = "T1"
    else:
        support = "T2"

    return {
        "locus": logical_locus, "GQS_T1": round(gqs_t1, 6), "GQS_T2": round(gqs_t2, 6),
        "ΔGQS": round(delta, 6), "support": support,
        "status": "success", "reason": None,
    }


def _validate_consistent_inputs(
    *, matrix: Path, candidate_trees: list[Path],
    tree_dir: Path, model_expr: str | None,
    partitions: Path | None, locus_ranges: Path | None,
    tool_args: str | None, guide_tree: Path | None,
) -> list[str]:
    errors: list[str] = []
    if not matrix.exists():
        errors.append(f"--matrix does not exist: {matrix}")
    if len(candidate_trees) != 2:
        errors.append(f"--candidate-trees must be exactly 2 trees, got {len(candidate_trees)}")
    for i, ct in enumerate(candidate_trees):
        if not ct.exists() or not ct.is_file():
            errors.append(f"--candidate-trees #{i+1} does not exist: {ct}")
        elif ct.stat().st_size == 0:
            errors.append(f"--candidate-trees #{i+1} is empty: {ct}")
    if not tree_dir.exists() or not tree_dir.is_dir():
        errors.append(f"--tree-dir does not exist: {tree_dir}")
    if partitions and locus_ranges:
        errors.append("--partitions and --locus-ranges are mutually exclusive")
    if model_expr and partitions:
        errors.append("--model-expr and --partitions are mutually exclusive")
    if not partitions and not locus_ranges:
        errors.append("Must provide --partitions or --locus-ranges (GLS requires locus boundaries)")
    has_tool_model = False
    if tool_args:
        toks = set(shlex.split(tool_args))
        has_tool_model = "-m" in toks or "-p" in toks
        for flag in _CONSISTENT_BLOCKED_FLAGS:
            if flag in toks:
                errors.append(f"Blocked flag in --tool-args: {flag}")
    if not model_expr and not partitions and not has_tool_model:
        errors.append("Must specify --model-expr, --partitions, or -m/-p in --tool-args")
    return errors


def run_signal_consistent(
    *,
    matrix: Path,
    candidate_trees: list[Path],
    tree_dir: Path,
    model_expr: str | None = None,
    partitions: Path | None = None,
    partition_mode: str = "p",
    locus_ranges: Path | None = None,
    guide_tree: Path | None = None,
    threads: str = "auto",
    iqtree_path: str | None = None,
    wastral_path: str | None = None,
    tool_args: str | None = None,
    metrics: Path | None = None,
    output_dir: Path | None = None,
    overwrite: bool = False,
    dry_run: bool = False,
    quiet: bool = False,
) -> dict[str, Any]:
    run_start = _time.time()
    if output_dir is None:
        output_dir = Path("runs/posttree/signal/consistent")
    output_dir = output_dir.resolve()
    matrix = matrix.resolve()
    candidate_trees = [ct.resolve() for ct in candidate_trees]
    tree_dir = tree_dir.resolve()

    errors = _validate_consistent_inputs(
        matrix=matrix, candidate_trees=candidate_trees,
        tree_dir=tree_dir, model_expr=model_expr,
        partitions=partitions, locus_ranges=locus_ranges,
        tool_args=tool_args, guide_tree=guide_tree,
    )

    params: dict[str, Any] = {
        "matrix": str(matrix),
        "candidate_trees_raw": ",".join(str(ct) for ct in candidate_trees),
        "tree_dir": str(tree_dir),
        "model_expr": model_expr,
        "partitions": str(partitions.resolve()) if partitions else None,
        "partition_mode": partition_mode if partitions else None,
        "locus_ranges": str(locus_ranges.resolve()) if locus_ranges else None,
        "guide_tree": str(guide_tree.resolve()) if guide_tree else None,
        "threads": threads,
        "iqtree_path": iqtree_path,
        "wastral_path": wastral_path,
        "tool_args": tool_args,
        "metrics": str(metrics.resolve()) if metrics else None,
        "output_dir": str(output_dir),
        "overwrite": overwrite,
        "dry_run": dry_run,
        "quiet": quiet,
    }

    if errors:
        return {"status": "error", "command": "", "wall_time": 0.0,
                "tool_versions": {}, "params": params, "key_results": {},
                "error": "; ".join(errors), "error_category": "input",
                "data": {"cmd": [], "tool_stderr": "", "tool_log": None, "output_files": {}}}

    if not dry_run:
        if overwrite and output_dir.exists():
            import shutil; shutil.rmtree(output_dir)
        elif output_dir.exists() and any(output_dir.iterdir()):
            return {"status": "error", "command": "", "wall_time": 0.0,
                    "tool_versions": {}, "params": params, "key_results": {},
                    "error": f"Output directory '{output_dir}' already exists and is non-empty. Use --overwrite to replace it.",
                    "error_category": "input",
                    "data": {"cmd": [], "tool_stderr": "", "tool_log": None, "output_files": {}}}
        output_dir.mkdir(parents=True, exist_ok=True)

    # Merge candidate trees
    if len(candidate_trees) == 1:
        candidate_trees_path = candidate_trees[0]
    else:
        candidate_trees_path = output_dir / "candidate.trees"
        if not dry_run:
            with open(candidate_trees_path, "w") as fh:
                for ct in candidate_trees:
                    fh.write(ct.read_text().strip() + "\n")

    prefix = "consistent"
    iqtree_flag = f"-{partition_mode}" if partitions else "-m"
    resolved_partitions = str(partitions.resolve()) if partitions else None
    resolved_guide_tree = str(guide_tree.resolve()) if guide_tree else None

    try:
        iqtree_exe = _resolve_iqtree_path(iqtree_path, dry_run)
    except (ValueError, FileNotFoundError) as e:
        return {"status": "error", "command": "", "wall_time": 0.0,
                "tool_versions": {}, "params": params, "key_results": {},
                "error": str(e), "error_category": "env",
                "data": {"cmd": [], "tool_stderr": "", "tool_log": None, "output_files": {}}}

    # Resolve wastral
    from phyloai.tree.msc import _resolve_wastral_path
    try:
        wastral_exe = _resolve_wastral_path(wastral_path, dry_run)
    except (ValueError, FileNotFoundError) as e:
        return {"status": "error", "command": "", "wall_time": 0.0,
                "tool_versions": {}, "params": params, "key_results": {},
                "error": str(e), "error_category": "env",
                "data": {"cmd": [], "tool_stderr": "", "tool_log": None, "output_files": {}}}

    if not dry_run:
        from phyloai.tree.msc import _detect_wastral_version
        tool_versions = {**_detect_iqtree_version(iqtree_exe), **_detect_wastral_version(wastral_exe)}
    else:
        tool_versions = {"iqtree3": "dry-run", "wastral": "dry-run"}

    # Build IQ-TREE command (same as lnl, -wslr)
    cmd: list[str] = [iqtree_exe, "-s", str(matrix), "-z", str(candidate_trees_path)]
    tool_toks = set(shlex.split(tool_args)) if tool_args else set()
    if "--prefix" not in tool_toks:
        cmd.extend(["--prefix", prefix])
    if model_expr and "-m" not in tool_toks:
        cmd.extend(["-m", model_expr])
    elif partitions:
        p_flag = f"-{partition_mode}"
        if p_flag not in tool_toks:
            cmd.extend([p_flag, resolved_partitions])
    cmd.append("-wslr")
    if resolved_guide_tree and "-ft" not in tool_toks:
        cmd.extend(["-ft", resolved_guide_tree])
    if "-T" not in tool_toks:
        cmd.extend(["-T", str(threads)])
    if tool_args:
        cmd.extend(shlex.split(tool_args))

    # CLI command string
    cli_parts = ["phyloai", "posttree", "signal", "consistent",
                 "--matrix", str(matrix),
                 "--candidate-trees", params["candidate_trees_raw"],
                 "--tree-dir", str(tree_dir)]
    if model_expr: cli_parts.extend(["--model-expr", model_expr])
    if partitions: cli_parts.extend(["--partitions", str(partitions), "--partition-mode", partition_mode])
    if locus_ranges: cli_parts.extend(["--locus-ranges", str(locus_ranges)])
    if guide_tree: cli_parts.extend(["--guide-tree", str(guide_tree)])
    if tool_args: cli_parts.extend(["--tool-args", tool_args])
    if metrics: cli_parts.extend(["--metrics", str(metrics)])
    cli_parts.extend(["--threads", str(threads), "-o", str(output_dir)])
    if iqtree_path: cli_parts.extend(["--iqtree-path", iqtree_path])
    if wastral_path: cli_parts.extend(["--wastral-path", wastral_path])
    if overwrite: cli_parts.append("--overwrite")
    if dry_run: cli_parts.append("--dry-run")
    if quiet: cli_parts.append("-q")
    full_command = shlex.join(cli_parts)

    # Fail-fast: validate T1/T2 taxon sets BEFORE running IQ-TREE
    t1_str = candidate_trees[0].read_text().strip()
    t2_str = candidate_trees[1].read_text().strip()
    t1_taxa = _get_tree_taxa(t1_str)
    t2_taxa = _get_tree_taxa(t2_str)
    if t1_taxa != t2_taxa:
        diff = t1_taxa.symmetric_difference(t2_taxa)
        result = {"status": "error", "command": full_command,
                  "wall_time": _time.time() - run_start,
                  "tool_versions": tool_versions, "params": params, "key_results": {},
                  "error": f"Candidate trees T1 and T2 have different taxon sets: {', '.join(sorted(diff))}",
                  "error_category": "input",
                  "data": {"cmd": cmd, "tool_stderr": "", "tool_log": None, "output_files": {}}}
        _write_result_json(output_dir, result)
        return result

    if dry_run:
        return {"status": "success", "command": full_command,
                "wall_time": 0.0, "tool_versions": tool_versions,
                "params": params, "key_results": {}, "error": None,
                "data": {"cmd": cmd, "tool_stderr": "", "tool_log": None, "output_files": {}}}

    # Run IQ-TREE
    proc = subprocess.run(cmd, capture_output=True, text=True, cwd=str(output_dir))
    tool_stderr = (proc.stdout or "") + ("\n" if proc.stdout and proc.stderr else "") + (proc.stderr or "")
    tool_log_path = output_dir / f"{prefix}.log"

    if proc.returncode != 0:
        result = {"status": "error", "command": full_command,
                  "wall_time": _time.time() - run_start,
                  "tool_versions": tool_versions, "params": params, "key_results": {},
                  "error": f"IQ-TREE exited {proc.returncode}", "error_category": "tool",
                  "data": {"cmd": cmd, "tool_stderr": tool_stderr,
                           "tool_log": str(tool_log_path), "output_files": {}}}
        _write_result_json(output_dir, result)
        return result

    # Parse .sitelh
    sitelh_path = output_dir / f"{prefix}.sitelh"
    tree_labels, site_scores = _parse_sitelh(sitelh_path)

    # GLS computation
    boundary_path = (partitions or locus_ranges).resolve()
    partition_recs = _parse_partition_ranges(boundary_path)
    partition_loci = {rec["locus"]: rec for rec in partition_recs}

    # Validate locus <-> gene tree matching
    try:
        gene_tree_map = scan_tree_dir(tree_dir)  # {logical_locus: Path}
    except ValueError as exc:
        result = {"status": "error", "command": full_command,
                  "wall_time": _time.time() - run_start,
                  "tool_versions": tool_versions, "params": params, "key_results": {},
                  "error": f"--tree-dir contains ambiguous or duplicate filenames: {exc}",
                  "error_category": "input",
                  "data": {"cmd": cmd, "tool_stderr": "", "tool_log": None, "output_files": {}}}
        _write_result_json(output_dir, result)
        return result
    tree_loci = set(gene_tree_map.keys())
    partition_loci_set = set(partition_loci.keys())
    missing_trees = partition_loci_set - tree_loci
    extra_trees = tree_loci - partition_loci_set
    if missing_trees or extra_trees:
        msg_parts = []
        if missing_trees:
            msg_parts.append(f"Loci in partition file with no gene tree: {', '.join(sorted(missing_trees))}")
        if extra_trees:
            msg_parts.append(f"Gene tree files with no matching partition locus: {', '.join(sorted(extra_trees))}")
        result = {"status": "error", "command": full_command,
                  "wall_time": _time.time() - run_start,
                  "tool_versions": tool_versions, "params": params, "key_results": {},
                  "error": "; ".join(msg_parts), "error_category": "input",
                  "data": {"cmd": cmd, "tool_stderr": "", "tool_log": None, "output_files": {}}}
        _write_result_json(output_dir, result)
        return result

    # GLS rows
    gls_rows = []
    for rec in partition_recs:
        gene_scores = _sum_gene_lnl(site_scores, rec["start"], rec["end"])
        delta = _delta_score(gene_scores)  # always 2 trees
        support = _support_label(gene_scores, tree_labels)
        gls_rows.append({
            "locus": rec["locus"],
            "lnL_T1": round(gene_scores[0], 6),
            "lnL_T2": round(gene_scores[1], 6),
            "ΔGLS": round(delta, 6),
            "support": support,
            "support_sig": abs(delta) >= 2.0,
        })

    gls_csv_path = output_dir / "gls.csv"
    _write_csv(gls_csv_path, gls_rows,
               ["locus", "lnL_T1", "lnL_T2", "ΔGLS", "support", "support_sig"])

    # T1/T2 taxon validation already done before IQ-TREE; ref_taxa set there
    ref_taxa = t1_taxa
    n_workers = os.cpu_count() if threads == "auto" else int(threads)
    work_dir = output_dir / "_gqs_work"
    work_dir.mkdir(parents=True, exist_ok=True)

    gqs_rows: list[dict[str, Any]] = []
    with ProcessPoolExecutor(max_workers=n_workers) as executor:
        futures = {
            executor.submit(
                _compute_gqs_for_locus,
                gene_tree_map[locus], locus, t1_str, t2_str, ref_taxa, wastral_exe, work_dir,
            ): locus
            for locus in partition_loci_set
        }
        locus_to_gqs: dict[str, dict] = {}
        for future in as_completed(futures):
            locus = futures[future]
            try:
                locus_to_gqs[locus] = future.result()
            except ValueError as exc:
                # input error (extra taxa) → exit code 1
                result = {"status": "error", "command": full_command,
                          "wall_time": _time.time() - run_start,
                          "tool_versions": tool_versions, "params": params, "key_results": {},
                          "error": str(exc), "error_category": "input",
                          "data": {"cmd": cmd, "tool_stderr": "", "tool_log": None, "output_files": {}}}
                _write_result_json(output_dir, result)
                return result
            except RuntimeError as exc:
                # tool failure (wastral non-zero exit or missing Score:) → exit code 2
                result = {"status": "error", "command": full_command,
                          "wall_time": _time.time() - run_start,
                          "tool_versions": tool_versions, "params": params, "key_results": {},
                          "error": str(exc), "error_category": "tool",
                          "data": {"cmd": cmd, "tool_stderr": str(exc), "tool_log": None, "output_files": {}}}
                _write_result_json(output_dir, result)
                return result

    # Maintain partition order
    for rec in partition_recs:
        gqs_rows.append(locus_to_gqs[rec["locus"]])

    gqs_csv_path = output_dir / "gqs.csv"
    _write_csv(gqs_csv_path, gqs_rows,
               ["locus", "GQS_T1", "GQS_T2", "ΔGQS", "support", "status", "reason"])

    # Consistency determination
    gls_map = {r["locus"]: r["support"] for r in gls_rows}
    gqs_map = {r["locus"]: r for r in gqs_rows}
    consistent_loci, inconsistent_loci = [], []
    for locus in partition_loci_set:
        gls_sup = gls_map.get(locus, "ambiguous")
        gqs_row = gqs_map.get(locus, {})
        gqs_sup = gqs_row.get("support", "ambiguous")
        gqs_status = gqs_row.get("status", "skipped")
        if gls_sup == gqs_sup and gls_sup != "ambiguous" and gqs_status == "success":
            consistent_loci.append(locus)
        else:
            inconsistent_loci.append(locus)

    (output_dir / "consistent_genes.txt").write_text("\n".join(sorted(consistent_loci)) + "\n")
    (output_dir / "inconsistent_genes.txt").write_text("\n".join(sorted(inconsistent_loci)) + "\n")

    # Support plots
    gls_pdf = _plot_support_bar([r["support"] for r in gls_rows], tree_labels,
                                 output_dir / "gls_support.pdf",
                                 ylabel="Number of genes")
    gqs_pdf = _plot_support_bar(
        [r["support"] for r in gqs_rows if r["status"] == "success"],
        tree_labels, output_dir / "gqs_support.pdf", ylabel="Number of genes",
    )

    n_gqs_skipped = sum(1 for r in gqs_rows if r["status"] == "skipped")
    output_files: dict[str, Any] = {
        "gls": {"path": str(gls_csv_path), "description": "Gene-wise lnL scores, ΔGLS, support"},
        "gqs": {"path": str(gqs_csv_path), "description": "Gene quartet scores, ΔGQS, support, status"},
        "consistent_genes": {"path": str(output_dir / "consistent_genes.txt"),
                              "description": "Loci where GLS and GQS support agree"},
        "inconsistent_genes": {"path": str(output_dir / "inconsistent_genes.txt"),
                                "description": "Loci where GLS and GQS support disagree or ambiguous"},
        "gls_support_plot": {"path": str(gls_pdf), "description": "GLS support distribution bar chart"},
        "gqs_support_plot": {"path": str(gqs_pdf), "description": "GQS support distribution bar chart"},
        "iqtree_report": {"path": str(output_dir / f"{prefix}.iqtree"), "description": "IQ-TREE native report"},
        "iqtree_sitelh": {"path": str(sitelh_path), "description": "IQ-TREE raw site log-likelihoods"},
    }

    # Metrics comparison (optional)
    if metrics:
        csv_p, pdf_p = _compare_groups(
            consistent_loci, inconsistent_loci, metrics.resolve(),
            "consistent", "inconsistent",
            output_dir, "consistent_comparison.csv", "consistent_comparison.pdf",
        )
        output_files["consistent_comparison"] = {"path": str(csv_p),
            "description": "Consistent vs inconsistent per-metric means and Wilcoxon p-values"}
        output_files["consistent_comparison_plot"] = {"path": str(pdf_p),
            "description": "Consistent vs inconsistent metric distribution boxplots"}

    key_results = {
        "n_loci": len(partition_recs),
        "n_consistent": len(consistent_loci),
        "n_inconsistent": len(inconsistent_loci),
        "n_gqs_skipped": n_gqs_skipped,
    }
    wall_time = _time.time() - run_start
    result = {
        "status": "success", "command": full_command, "wall_time": wall_time,
        "tool_versions": tool_versions, "params": params,
        "key_results": key_results, "error": None,
        "data": {
            "cmd": cmd, "tool_stderr": tool_stderr,
            "tool_log": str(tool_log_path),
            "summary": {**key_results, "wastral_n_gene_trees": len(partition_recs),
                        "wastral_threads_used": n_workers},
            "output_files": output_files,
        },
    }
    _write_result_json(output_dir, result)
    return result
```

- [ ] **Step 4: Run validation tests**

```bash
python -m pytest tests/posttree/test_signal_consistent.py -v
```

Expected: PASS.

- [ ] **Step 5: Write integration test**

```python
class TestRunSignalConsistentFixture:
    SIGNAL_DIR = Path("runs/signal")

    def test_consistent_with_fixture(self, tmp_path: Path) -> None:
        from phyloai.posttree.signal import run_signal_consistent
        if not self.SIGNAL_DIR.exists():
            pytest.skip("Signal test data not present")
        matrix = self.SIGNAL_DIR / "matrix.aa.fa"
        t1 = self.SIGNAL_DIR / "T1.tre"
        t2 = self.SIGNAL_DIR / "T2.tre"
        tree_dir = self.SIGNAL_DIR / "gene_trees1066"
        partitions = self.SIGNAL_DIR / "matrix.aa.partitions"
        result = run_signal_consistent(
            matrix=matrix, candidate_trees=[t1, t2],
            tree_dir=tree_dir, model_expr="LG+F+R4",
            partitions=None, partition_mode="p",
            locus_ranges=partitions,
            guide_tree=None, threads="auto",
            iqtree_path=None, wastral_path=None,
            tool_args=None, metrics=None,
            output_dir=tmp_path / "consistent_out",
            overwrite=False, dry_run=False, quiet=True,
        )
        if result["status"] == "error" and result.get("error_category") == "env":
            pytest.skip("iqtree3 or wastral not available")
        assert result["status"] == "success"
        assert result["key_results"]["n_loci"] == 20
        assert (tmp_path / "consistent_out" / "consistent_genes.txt").exists()
        assert (tmp_path / "consistent_out" / "inconsistent_genes.txt").exists()
```

- [ ] **Step 6: Run all consistent tests**

```bash
python -m pytest tests/posttree/test_signal_consistent.py -v
```

---

## Task 4: `signal fclm` — library function `run_signal_fclm`

**Files:**
- Modify: `phyloai/posttree/signal.py`
- Create: `tests/posttree/test_signal_fclm.py`

**Interfaces produced:**
- `run_signal_fclm(*, matrix, taxset_csv, model_expr, partitions, partition_mode, lmap, guide_tree, threads, iqtree_path, tool_args, output_dir, overwrite, dry_run, quiet) -> dict`

- [ ] **Step 1: Write failing tests**

```python
# tests/posttree/test_signal_fclm.py
"""Tests for phyloai.posttree.signal fclm subcommand."""
from __future__ import annotations
from pathlib import Path
import pytest


class TestFclmValidation:
    def test_reads_phylip_taxa_for_taxset_validation(self, tmp_path: Path) -> None:
        from phyloai.posttree.signal import run_signal_fclm
        matrix = tmp_path / "m.phy"
        matrix.write_text("4 3\nA MKT\nB MKA\nC MKL\nD MKP\n")
        csv = tmp_path / "taxsets.csv"
        csv.write_text("taxon,taxset\nA,G1\nB,G2\nC,G3\nD,G4\n")
        result = run_signal_fclm(
            matrix=matrix, taxset_csv=csv, model_expr="LG+F+R4", lmap=None,
            guide_tree=None, threads="auto", iqtree_path=None, tool_args=None,
            output_dir=tmp_path / "out", overwrite=False, dry_run=True, quiet=True,
        )
        assert result["status"] == "success"
        assert result["key_results"]["n_taxsets"] == 4

    def test_taxset_csv_extra_taxa_errors(self, tmp_path: Path) -> None:
        from phyloai.posttree.signal import run_signal_fclm
        matrix = tmp_path / "m.fa"
        matrix.write_text(">A\nMKT\n>B\nMKA\n>C\nMKL\n>D\nMKP\n")
        csv = tmp_path / "taxsets.csv"
        csv.write_text("taxon,taxset\nA,G1\nB,G2\nC,G3\nD,G4\nX,G1\n")  # X not in matrix
        result = run_signal_fclm(
            matrix=matrix, taxset_csv=csv,
            model_expr="LG+F+R4", lmap=None,
            guide_tree=None, threads="auto",
            iqtree_path=None, tool_args=None,
            output_dir=tmp_path / "out",
            overwrite=False, dry_run=True, quiet=True,
        )
        assert result["status"] == "error"
        assert "X" in result["error"] or "taxa" in result["error"].lower()

    def test_taxset_csv_fewer_than_4_taxsets_errors(self, tmp_path: Path) -> None:
        from phyloai.posttree.signal import run_signal_fclm
        matrix = tmp_path / "m.fa"
        matrix.write_text(">A\nMKT\n>B\nMKA\n>C\nMKL\n>D\nMKP\n")
        csv = tmp_path / "taxsets.csv"
        csv.write_text("taxon,taxset\nA,G1\nB,G1\nC,G2\nD,G3\n")  # only 3 taxsets
        result = run_signal_fclm(
            matrix=matrix, taxset_csv=csv,
            model_expr="LG+F+R4", lmap=None,
            guide_tree=None, threads="auto",
            iqtree_path=None, tool_args=None,
            output_dir=tmp_path / "out",
            overwrite=False, dry_run=True, quiet=True,
        )
        assert result["status"] == "error"
        assert "4" in result["error"]

    def test_dry_run_produces_nexus(self, tmp_path: Path) -> None:
        from phyloai.posttree.signal import run_signal_fclm
        matrix = tmp_path / "m.fa"
        matrix.write_text(">A\nMKT\n>B\nMKA\n>C\nMKL\n>D\nMKP\n>E\nMKQ\n")
        csv = tmp_path / "taxsets.csv"
        csv.write_text("taxon,taxset\nA,G1\nB,G2\nC,G3\nD,G4\nE,G1\n")
        result = run_signal_fclm(
            matrix=matrix, taxset_csv=csv,
            model_expr="LG+F+R4", lmap=None,
            guide_tree=None, threads="auto",
            iqtree_path=None, tool_args=None,
            output_dir=tmp_path / "out",
            overwrite=False, dry_run=True, quiet=True,
        )
        assert result["status"] == "success"
        assert "-lmap" in result["data"]["cmd"]
        assert "-lmclust" in result["data"]["cmd"]
        nexus = tmp_path / "out" / "cluster.nexus"
        assert nexus.exists()
        content = nexus.read_text()
        assert "taxset G1" in content
```

- [ ] **Step 2: Run to verify failure**

```bash
python -m pytest tests/posttree/test_signal_fclm.py -v
```

Expected: ImportError.

- [ ] **Step 3: Implement `run_signal_fclm`**

Add to `phyloai/posttree/signal.py`:

```python
_FCLM_BLOCKED_FLAGS = frozenset({"-s", "-lmap", "-lmclust", "-n"})


def _read_matrix_taxa(matrix: Path) -> set[str]:
    """Read taxon IDs from any supported alignment format.

    Uses the project format reader rather than parsing FASTA headers, so FcLM
    accepts the same FASTA, PHYLIP, PHYLIP-PAML, and NEXUS inputs as IQ-TREE.
    """
    alignment = FormatConverter().read(matrix)
    return {record.id for record in alignment}


def _csv_to_nexus(taxset_csv: Path, output_path: Path) -> dict[str, list[str]]:
    """Convert taxon,taxset CSV to NEXUS cluster file. Returns {taxset: [taxa]}."""
    taxset_map: dict[str, list[str]] = {}
    with open(taxset_csv, newline="") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            taxon = row["taxon"].strip()
            taxset = row["taxset"].strip()
            taxset_map.setdefault(taxset, []).append(taxon)
    lines = ["#NEXUS", "begin sets;"]
    for ts, taxa in taxset_map.items():
        lines.append(f"  taxset {ts} = {' '.join(taxa)};")
    lines.append("end;")
    output_path.write_text("\n".join(lines) + "\n")
    return taxset_map


def _validate_fclm_inputs(
    *, matrix: Path, taxset_csv: Path,
    model_expr: str | None, partitions: Path | None = None,
    partition_mode: str | None = None,
    tool_args: str | None = None,
    guide_tree: Path | None = None,
    lmap: str | None = None,
) -> list[str]:
    errors: list[str] = []
    if not matrix.exists():
        errors.append(f"--matrix does not exist: {matrix}")
    if model_expr and partitions:
        errors.append("--model-expr and --partitions are mutually exclusive")
    if not taxset_csv.exists():
        errors.append(f"--taxset-csv does not exist: {taxset_csv}")
        return errors  # can't validate further

    if guide_tree and (not guide_tree.exists() or not guide_tree.is_file()):
        errors.append(f"--guide-tree does not exist or is not a regular file: {guide_tree}")
    if partitions and (not partitions.exists() or not partitions.is_file()):
        errors.append(f"--partitions does not exist or is not a regular file: {partitions}")
    if lmap is not None and lmap != "ALL":

    # Read taxa
    matrix_taxa = _read_matrix_taxa(matrix) if matrix.exists() else set()
    csv_taxa: dict[str, str] = {}
    taxsets: dict[str, list[str]] = {}
    try:
        with open(taxset_csv, newline="") as fh:
            reader = csv.DictReader(fh)
            for row in reader:
                taxon = row["taxon"].strip()
                taxset = row["taxset"].strip()
                if taxon in csv_taxa:
                    errors.append(f"Taxon {taxon!r} appears in multiple taxsets (not mutually exclusive)")
                csv_taxa[taxon] = taxset
                taxsets.setdefault(taxset, []).append(taxon)
    except Exception as exc:
        errors.append(f"Cannot read --taxset-csv: {exc}")
        return errors

    extra_taxa = set(csv_taxa) - matrix_taxa
    if extra_taxa:
        errors.append(f"Taxa in --taxset-csv not found in --matrix: {', '.join(sorted(extra_taxa))}")
    missing_taxa = matrix_taxa - set(csv_taxa)
    if missing_taxa:
        errors.append(f"Taxa in --matrix not assigned in --taxset-csv: {', '.join(sorted(missing_taxa))}")
    if len(taxsets) < 4:
        errors.append(f"FcLM requires at least 4 taxsets, got {len(taxsets)}")

    has_model = model_expr or partitions
    has_tool_model = False
    if tool_args:
        toks = set(shlex.split(tool_args))
        has_tool_model = "-m" in toks or "-p" in toks or "-Q" in toks
        for flag in _FCLM_BLOCKED_FLAGS:
            if flag in toks:
                errors.append(f"Blocked flag in --tool-args: {flag}")
    if not has_model and not has_tool_model:
        errors.append("Must specify --model-expr, --partitions, or -m/-p/-Q in --tool-args")
    return errors


def run_signal_fclm(
    *,
    matrix: Path,
    taxset_csv: Path,
    model_expr: str | None = None,
    partitions: Path | None = None,
    partition_mode: str = "p",
    lmap: str | None = None,
    guide_tree: Path | None = None,
    threads: str = "auto",
    iqtree_path: str | None = None,
    tool_args: str | None = None,
    output_dir: Path | None = None,
    overwrite: bool = False,
    dry_run: bool = False,
    quiet: bool = False,
) -> dict[str, Any]:
    run_start = _time.time()
    if output_dir is None:
        output_dir = Path("runs/posttree/signal/fclm")
    output_dir = output_dir.resolve()
    matrix = matrix.resolve()
    taxset_csv = taxset_csv.resolve()

    errors = _validate_fclm_inputs(
        matrix=matrix, taxset_csv=taxset_csv,
        model_expr=model_expr, partitions=partitions,
        partition_mode=partition_mode,
        tool_args=tool_args, guide_tree=guide_tree, lmap=lmap,
    )

    params: dict[str, Any] = {
        "matrix": str(matrix),
        "taxset_csv": str(taxset_csv),
        "model_expr": model_expr if model_expr else None,
        "partitions": str(partitions.resolve()) if partitions else None,
        "partition_mode": partition_mode if partitions else None,
        "lmap": lmap,
        "guide_tree": str(guide_tree.resolve()) if guide_tree else None,
        "threads": threads,
        "iqtree_path": iqtree_path,
        "tool_args": tool_args,
        "output_dir": str(output_dir),
        "overwrite": overwrite,
        "dry_run": dry_run,
        "quiet": quiet,
    }

    if errors:
        return {"status": "error", "command": "", "wall_time": 0.0,
                "tool_versions": {}, "params": params, "key_results": {},
                "error": "; ".join(errors), "error_category": "input",
                "data": {"cmd": [], "tool_stderr": "", "tool_log": None, "output_files": {}}}

    if not dry_run:
        if overwrite and output_dir.exists():
            import shutil; shutil.rmtree(output_dir)
        elif output_dir.exists() and any(output_dir.iterdir()):
            return {"status": "error", "command": "", "wall_time": 0.0,
                    "tool_versions": {}, "params": params, "key_results": {},
                    "error": f"Output directory '{output_dir}' already exists and is non-empty. Use --overwrite to replace it.",
                    "error_category": "input",
                    "data": {"cmd": [], "tool_stderr": "", "tool_log": None, "output_files": {}}}
    output_dir.mkdir(parents=True, exist_ok=True)

    # Convert CSV to NEXUS
    nexus_path = output_dir / "cluster.nexus"
    taxset_map = _csv_to_nexus(taxset_csv, nexus_path)
    n_taxsets = len(taxset_map)

    # Resolve lmap value
    if lmap is None:
        matrix_taxa = _read_matrix_taxa(matrix)
        lmap_val = str(50 * len(matrix_taxa))
    else:
        lmap_val = lmap  # "ALL" or integer string

    prefix = "fclm"
    try:
        iqtree_exe = _resolve_iqtree_path(iqtree_path, dry_run)
    except (ValueError, FileNotFoundError) as e:
        return {"status": "error", "command": "", "wall_time": 0.0,
                "tool_versions": {}, "params": params, "key_results": {},
                "error": str(e), "error_category": "env",
                "data": {"cmd": [], "tool_stderr": "", "tool_log": None, "output_files": {}}}

    tool_versions = _detect_iqtree_version(iqtree_exe) if not dry_run else {"iqtree3": "dry-run"}

    tool_toks = set(shlex.split(tool_args)) if tool_args else set()
    cmd = [iqtree_exe, "-s", str(matrix)]
    if model_expr and "-m" not in tool_toks:
        cmd.extend(["-m", model_expr])
    elif partitions:
        cmd.extend([f"-{partition_mode}", str(partitions.resolve())])
    if guide_tree and "-ft" not in tool_toks:
        cmd.extend(["-ft", str(guide_tree.resolve())])
    cmd.extend(["-lmap", lmap_val, "-lmclust", str(nexus_path), "-n", "0"])
    if "--prefix" not in tool_toks:
        cmd.extend(["--prefix", prefix])
    if "-T" not in tool_toks:
        cmd.extend(["-T", str(threads)])
    if tool_args:
        cmd.extend(shlex.split(tool_args))

    # CLI command string
    cli_parts = ["phyloai", "posttree", "signal", "fclm",
                 "--matrix", str(matrix), "--taxset-csv", str(taxset_csv)]
    if model_expr: cli_parts.extend(["--model-expr", model_expr])
    if lmap: cli_parts.extend(["--lmap", lmap])
    if guide_tree: cli_parts.extend(["--guide-tree", str(guide_tree)])
    if tool_args: cli_parts.extend(["--tool-args", tool_args])
    cli_parts.extend(["--threads", str(threads), "-o", str(output_dir)])
    if iqtree_path: cli_parts.extend(["--iqtree-path", iqtree_path])
    if overwrite: cli_parts.append("--overwrite")
    if dry_run: cli_parts.append("--dry-run")
    if quiet: cli_parts.append("-q")
    full_command = shlex.join(cli_parts)

    if dry_run:
        return {"status": "success", "command": full_command,
                "wall_time": 0.0, "tool_versions": tool_versions,
                "params": params,
                "key_results": {"n_taxsets": n_taxsets, "n_quartets": lmap_val},
                "error": None,
                "data": {"cmd": cmd, "tool_stderr": "", "tool_log": None,
                         "output_files": {"cluster_nexus": {"path": str(nexus_path),
                                          "description": "NEXUS cluster file for IQ-TREE -lmclust"}}}}

    proc = subprocess.run(cmd, capture_output=True, text=True, cwd=str(output_dir))
    tool_stderr = (proc.stdout or "") + ("\n" if proc.stdout and proc.stderr else "") + (proc.stderr or "")
    tool_log_path = output_dir / f"{prefix}.log"

    output_files: dict[str, Any] = {
        "cluster_nexus": {"path": str(nexus_path), "description": "NEXUS cluster file for IQ-TREE -lmclust"},
        "lmap_figure": {"path": str(output_dir / f"{prefix}.lmap.eps"),
                        "description": "IQ-TREE likelihood mapping figure"},
        "iqtree_report": {"path": str(output_dir / f"{prefix}.iqtree"),
                          "description": "IQ-TREE native report (contains all lmap statistics)"},
    }

    status = "success" if proc.returncode == 0 else "error"
    error_msg = None if proc.returncode == 0 else f"IQ-TREE exited {proc.returncode}"

    result = {
        "status": status, "command": full_command,
        "wall_time": _time.time() - run_start,
        "tool_versions": tool_versions, "params": params,
        "key_results": {"n_taxsets": n_taxsets, "n_quartets": lmap_val},
        "error": error_msg,
        "error_category": None if status == "success" else "tool",
        "data": {"cmd": cmd, "tool_stderr": tool_stderr,
                 "tool_log": str(tool_log_path), "output_files": output_files},
    }
    _write_result_json(output_dir, result)
    return result
```

- [ ] **Step 4: Run fclm tests**

```bash
python -m pytest tests/posttree/test_signal_fclm.py -v
```

Expected: all PASS.

---

## Task 5: CLI wiring — Click commands in `posttree.py`

**Files:**
- Modify: `phyloai/cli/commands/posttree.py`

- [ ] **Step 1: Add signal Click group and three subcommands**

Add after the existing `dating` group at the bottom of `phyloai/cli/commands/posttree.py`:

```python
# ===================================================================
# Signal group
# ===================================================================

class _SignalGroup(click.Group):
    def list_commands(self, ctx: click.Context) -> list[str]:
        return ["lnl", "fclm", "consistent"]


@posttree.group("signal", cls=_SignalGroup)
def signal() -> None:
    """Phylogenetic signal distribution analysis."""


@signal.command("lnl")
@click.option("--matrix", required=True, type=click.Path(exists=True, path_type=Path), help="Supermatrix alignment.")
@click.option("--candidate-trees", "candidate_trees_raw", required=True, type=str, help="Tree-list file or comma-separated NEWICK files.")
@click.option("--model-expr", type=str, default=None, help="IQ-TREE model expression. Mutually exclusive with --partitions.")
@click.option("--partitions", type=click.Path(path_type=Path), default=None, help="Partition file; passed to IQ-TREE -p/-Q. Also used for locus boundary extraction. Mutually exclusive with --model-expr and --locus-ranges.")
@click.option("--partition-mode", type=click.Choice(["p", "Q"]), default="p", show_default=True, help="p=-p (edge-linked); Q=-Q (edge-unlinked). Only valid with --partitions.")
@click.option("--locus-ranges", type=click.Path(path_type=Path), default=None, help="Partition file for locus boundary extraction only (not passed to IQ-TREE). Mutually exclusive with --partitions.")
@click.option("--guide-tree", type=click.Path(path_type=Path), default=None, help="Guide tree for PMSF models.")
@click.option("--metrics", type=click.Path(path_type=Path), default=None, help="Metrics CSV from phyloai pretree metrics for outlier and tree-support-group comparisons.")
@click.option("--threads", "-t", default="auto", show_default=True, help="IQ-TREE -T value (integer or auto).")
@click.option("--iqtree-path", type=str, default=None, help="Explicit path to iqtree3 executable.")
@click.option("--tool-args", type=str, default=None, help="Extra IQ-TREE flags. Blocked: -s, -z, -wslr, -p, -Q.")
@click.option("--output-dir", "-o", type=click.Path(path_type=Path), default=Path("runs/posttree/signal/lnl"), show_default=True)
@click.option("--overwrite", is_flag=True, default=False)
@click.option("--dry-run", is_flag=True, default=False)
@click.option("--quiet", "-q", is_flag=True, default=False)
def lnl_command(
    matrix: Path, candidate_trees_raw: str,
    model_expr: str | None, partitions: Path | None,
    partition_mode: str,
    locus_ranges: Path | None, guide_tree: Path | None,
    metrics: Path | None, threads: str, iqtree_path: str | None,
    tool_args: str | None, output_dir: Path,
    overwrite: bool, dry_run: bool, quiet: bool,
) -> None:
    """Site-wise and gene-wise log-likelihood score distribution.

    Examples:

      # Homogeneous model, site-wise only

      phyloai posttree signal lnl --matrix matrix.fa --candidate-trees trees --model-expr LG+F+R4

      # With gene-wise output via locus ranges

      phyloai posttree signal lnl --matrix matrix.fa --candidate-trees trees --model-expr LG+F+R4 --locus-ranges partitions.txt

      # With outlier and tree-support-group metrics comparisons

      phyloai posttree signal lnl --matrix matrix.fa --candidate-trees trees --model-expr LG+F+R4 --locus-ranges partitions.txt --metrics metrics.csv
    """
    from phyloai.posttree.signal import run_signal_lnl

    # Parse candidate trees (comma-separated or single)
    candidate_trees = _parse_candidate_trees(candidate_trees_raw)

    if partitions is None:
        partition_mode = None
    result = run_signal_lnl(
        matrix=matrix, candidate_trees=candidate_trees,
        model_expr=model_expr, partitions=partitions,
        partition_mode=partition_mode,
        locus_ranges=locus_ranges, guide_tree=guide_tree,
        threads=threads, iqtree_path=iqtree_path,
        tool_args=tool_args, metrics=metrics,
        output_dir=output_dir, overwrite=overwrite,
        dry_run=dry_run, quiet=quiet,
    )

    import json
    result_path = output_dir.resolve() / "result.json"
    result_path.parent.mkdir(parents=True, exist_ok=True)
    with open(result_path, "w") as fh:
        json.dump(result, fh, indent=2)

    if result["status"] == "error":
        _fail(result.get("error") or "Unknown error",
              exit_code=1 if result.get("error_category") == "input" else
              3 if result.get("error_category") == "env" else 2)

    if not quiet:
        kr = result["key_results"]
        click.echo(f"\nStatus: {result['status']}")
        click.echo(f"Wall time: {result['wall_time']:.1f}s")
        click.echo(f"Trees: {kr.get('n_trees')}  Sites: {kr.get('n_sites')}")
        if kr.get("n_loci"):
            click.echo(f"Loci: {kr['n_loci']}  Outliers: {kr.get('n_outlier_genes')}")
        click.echo(f"Result: {result_path}")


@signal.command("consistent")
@click.option("--matrix", required=True, type=click.Path(exists=True, path_type=Path))
@click.option("--candidate-trees", "candidate_trees_raw", required=True, type=str, help="Exactly 2 candidate trees.")
@click.option("--tree-dir", required=True, type=click.Path(exists=True, path_type=Path), help="Directory of gene tree files for GQS.")
@click.option("--model-expr", type=str, default=None)
@click.option("--partitions", type=click.Path(path_type=Path), default=None, help="Partition file passed to IQ-TREE -p/-Q. Also extracts locus boundaries. Mutually exclusive with --model-expr and --locus-ranges.")
@click.option("--partition-mode", type=click.Choice(["p", "Q"]), default="p", show_default=True, help="p=-p (edge-linked proportional); Q=-Q (edge-unlinked, independent branch lengths per partition). Only valid with --partitions.")
@click.option("--locus-ranges", type=click.Path(path_type=Path), default=None, help="Partition file for locus boundary extraction only. Not passed to IQ-TREE. Mutually exclusive with --partitions.")
@click.option("--guide-tree", type=click.Path(path_type=Path), default=None, help="Guide tree for PMSF models.")
@click.option("--metrics", type=click.Path(path_type=Path), default=None)
@click.option("--threads", "-t", default="auto", show_default=True)
@click.option("--iqtree-path", type=str, default=None)
@click.option("--wastral-path", type=str, default=None, help="Explicit path to wastral executable.")
@click.option("--tool-args", type=str, default=None)
@click.option("--output-dir", "-o", type=click.Path(path_type=Path), default=Path("runs/posttree/signal/consistent"), show_default=True)
@click.option("--overwrite", is_flag=True, default=False)
@click.option("--dry-run", is_flag=True, default=False)
@click.option("--quiet", "-q", is_flag=True, default=False)
def consistent_command(
    matrix: Path, candidate_trees_raw: str, tree_dir: Path,
    model_expr: str | None, partitions: Path | None, partition_mode: str,
    locus_ranges: Path | None, guide_tree: Path | None,
    metrics: Path | None, threads: str,
    iqtree_path: str | None, wastral_path: str | None,
    tool_args: str | None, output_dir: Path,
    overwrite: bool, dry_run: bool, quiet: bool,
) -> None:
    """Consistent gene identification via GLS + GQS (Shen et al. 2021).

    Requires exactly 2 candidate trees. Identifies genes where both
    likelihood-based (GLS) and quartet-based (GQS) signal agree.

    Examples:

      phyloai posttree signal consistent --matrix matrix.fa --candidate-trees T1.tre,T2.tre --tree-dir gene_trees/ --model-expr LG+F+R4 --locus-ranges partitions.txt
    """
    from phyloai.posttree.signal import run_signal_consistent

    candidate_trees = _parse_candidate_trees(candidate_trees_raw)
    result = run_signal_consistent(
        matrix=matrix, candidate_trees=candidate_trees,
        tree_dir=tree_dir, model_expr=model_expr,
        partitions=partitions, partition_mode=partition_mode,
        locus_ranges=locus_ranges, guide_tree=guide_tree,
        threads=threads, iqtree_path=iqtree_path,
        wastral_path=wastral_path, tool_args=tool_args,
        metrics=metrics, output_dir=output_dir,
        overwrite=overwrite, dry_run=dry_run, quiet=quiet,
    )

    import json
    result_path = output_dir.resolve() / "result.json"
    result_path.parent.mkdir(parents=True, exist_ok=True)
    with open(result_path, "w") as fh:
        json.dump(result, fh, indent=2)

    if result["status"] == "error":
        _fail(result.get("error") or "Unknown error",
              exit_code=1 if result.get("error_category") == "input" else
              3 if result.get("error_category") == "env" else 2)

    if not quiet:
        kr = result["key_results"]
        click.echo(f"\nStatus: {result['status']}")
        click.echo(f"Wall time: {result['wall_time']:.1f}s")
        click.echo(f"Loci: {kr.get('n_loci')}  Consistent: {kr.get('n_consistent')}  Inconsistent: {kr.get('n_inconsistent')}")
        click.echo(f"Result: {result_path}")


@signal.command("fclm")
@click.option("--matrix", required=True, type=click.Path(exists=True, path_type=Path),
              help="Single supermatrix alignment (FASTA/PHYLIP/NEXUS).")
@click.option("--taxset-csv", required=True, type=click.Path(exists=True, path_type=Path),
              help="Two-column CSV (taxon,taxset) defining cluster membership. Minimum 4 taxsets.")
@click.option("--model-expr", type=str, default=None,
              help="IQ-TREE model expression. Mutually exclusive with --partitions.")
@click.option("--partitions", type=click.Path(path_type=Path), default=None,
              help="Partition file (e.g. .best_model.nex). Passed to IQ-TREE as -p or -Q.")
@click.option("--partition-mode", type=click.Choice(["p", "Q"]), default="p", show_default=True,
              help="p=-p (edge-linked); Q=-Q (edge-unlinked).")
@click.option("--lmap", type=str, default=None, help="Quartet count: ALL, integer, or omit for 50*n_taxa.")
@click.option("--guide-tree", type=click.Path(path_type=Path), default=None, help="Guide tree for PMSF models.")
@click.option("--threads", "-t", default="auto", show_default=True, help="IQ-TREE -T value (integer or auto).")
@click.option("--iqtree-path", type=str, default=None, help="Explicit path to iqtree3 executable.")
@click.option("--tool-args", type=str, default=None, help="Extra IQ-TREE flags. Blocked: -s, -lmap, -lmclust, -n, -p, -Q.")
@click.option("--output-dir", "-o", type=click.Path(path_type=Path), default=Path("runs/posttree/signal/fclm"), show_default=True)
@click.option("--overwrite", is_flag=True, default=False)
@click.option("--dry-run", is_flag=True, default=False)
@click.option("--quiet", "-q", is_flag=True, default=False)
def fclm_command(
    matrix: Path, taxset_csv: Path, model_expr: str | None,
    partitions: Path | None, partition_mode: str,
    lmap: str | None, guide_tree: Path | None, threads: str,
    iqtree_path: str | None, tool_args: str | None, output_dir: Path,
    overwrite: bool, dry_run: bool, quiet: bool,
) -> None:
    """Four-cluster Likelihood Mapping (FcLM).

    Assesses phylogenetic signal supporting alternative hypotheses among
    four taxon clusters using IQ-TREE3 -lmap -lmclust.

    Examples:

      # Homogeneous model

      phyloai posttree signal fclm --matrix matrix.fa --taxset-csv taxsets.csv --model-expr LG+C60+F+R4

      # Partition model

      phyloai posttree signal fclm --matrix matrix.fa --taxset-csv taxsets.csv --partitions matrix.best_model.nex
    """
    from phyloai.posttree.signal import run_signal_fclm

    result = run_signal_fclm(
        matrix=matrix, taxset_csv=taxset_csv,
        model_expr=model_expr, partitions=partitions,
        partition_mode=partition_mode,
        lmap=lmap, guide_tree=guide_tree,
        threads=threads, iqtree_path=iqtree_path,
        tool_args=tool_args, output_dir=output_dir,
        overwrite=overwrite, dry_run=dry_run, quiet=quiet,
    )

    import json
    result_path = output_dir.resolve() / "result.json"
    result_path.parent.mkdir(parents=True, exist_ok=True)
    with open(result_path, "w") as fh:
        json.dump(result, fh, indent=2)

    if result["status"] == "error":
        _fail(result.get("error") or "Unknown error",
              exit_code=1 if result.get("error_category") == "input" else
              3 if result.get("error_category") == "env" else 2)

    if not quiet:
        kr = result["key_results"]
        click.echo(f"\nStatus: {result['status']}")
        click.echo(f"Wall time: {result['wall_time']:.1f}s")
        click.echo(f"Taxsets: {kr.get('n_taxsets')}")
        click.echo(f"Report: {output_dir.resolve() / 'fclm.iqtree'}")
        click.echo(f"Result: {result_path}")
```

Also add the helper for parsing candidate trees (reuse topology pattern):

```python
def _parse_candidate_trees(candidate_trees_raw: str) -> list[Path]:
    """Parse comma-separated tree paths or single path."""
    if "," in candidate_trees_raw:
        return [Path(p.strip()) for p in candidate_trees_raw.split(",")]
    return [Path(candidate_trees_raw.strip())]
```

- [ ] **Step 2: Verify CLI is importable**

```bash
cd /Users/zf/data/coding/phyloAI-syserror
python -c "from phyloai.cli.commands.posttree import posttree; print('OK')"
```

Expected: `OK`

- [ ] **Step 3: Smoke test CLI help**

```bash
python -m phyloai posttree signal --help
python -m phyloai posttree signal lnl --help
python -m phyloai posttree signal consistent --help
python -m phyloai posttree signal fclm --help
```

Expected: help text displayed without errors.

- [ ] **Step 4: Dry-run CLI smoke test**

```bash
python -m phyloai posttree signal lnl \
  --matrix runs/signal/matrix.aa.fa \
  --candidate-trees runs/signal/trees \
  --model-expr LG+F+R4 \
  --dry-run --quiet \
  -o /tmp/signal_lnl_test
```

Expected: prints iqtree3 command with `-wslr`, exits 0.

---

## Task 6: Report integration

**Files:**
- Modify: `phyloai/report/collector.py`
- Modify: `phyloai/report/templates.py`

- [ ] **Step 1: Update `collector.py` STEP_ORDER and `_THIRD_LEVEL`**

In `phyloai/report/collector.py`, replace:

```python
    "posttree.signal",
```

with:

```python
    "posttree.signal.lnl",
    "posttree.signal.consistent",
    "posttree.signal.fclm",
```

And add `"signal"` to `_THIRD_LEVEL`:

```python
    _THIRD_LEVEL: dict[str, set[str]] = {
        "filter": {"taper", "treeshrink", "symtest", "metrics", "cluster"},
        "concat": {"jackknife"},
        "ml":     {"fasttree", "iqtree"},
        "dating": {"hessian", "mcmc"},
        "syserror": {"brlen", "cca", "sites"},
        "signal": {"lnl", "consistent", "fclm"},   # ADD THIS LINE
    }
```

- [ ] **Step 2: Update `templates.py` — replace stub with three generators**

In `phyloai/report/templates.py`, replace the existing `generate_methods_posttree_signal` function:

```python
def generate_methods_posttree_signal_lnl(
    params: dict[str, Any],
    key_results: dict[str, Any],
    tool_versions: dict[str, Any],
) -> str:
    n_trees = key_results.get("n_trees", 2)
    n_sites = key_results.get("n_sites", "")
    model = params.get("model_expr") or "a partition model"
    iqtree_ver = tool_versions.get("iqtree3", "IQ-TREE3")
    parts = [
        f"Site-wise and gene-wise log-likelihood scores were computed using {iqtree_ver} "
        f"({model}) across {_describe_n(n_trees, 'candidate topology', 'candidate topologies')} "
        f"and {_describe_n(n_sites, 'alignment site')} "
        f"following Shen et al. (2017)."
    ]
    if key_results.get("n_loci"):
        parts.append(
            f"Gene-wise ΔLnL (ΔGLS) was computed for {key_results['n_loci']} loci; "
            f"{key_results.get('n_outlier_genes', 0)} outlier genes were identified "
            f"using Tukey's 1.5×IQR criterion."
        )
    return " ".join(parts)


def generate_methods_posttree_signal_consistent(
    params: dict[str, Any],
    key_results: dict[str, Any],
    tool_versions: dict[str, Any],
) -> str:
    n_loci = key_results.get("n_loci", "")
    n_consistent = key_results.get("n_consistent", "")
    model = params.get("model_expr") or "a partition model"
    iqtree_ver = tool_versions.get("iqtree3", "IQ-TREE3")
    wastral_ver = tool_versions.get("wastral", "wASTRAL")
    return (
        f"Consistent genes were identified across {_describe_n(n_loci, 'locus', 'loci')} "
        f"following Shen et al. (2021) using gene-wise log-likelihood scores (GLS) "
        f"computed with {iqtree_ver} ({model}) and quartet concordance scores (GQS) "
        f"computed with {wastral_ver}. "
        f"{_describe_n(n_consistent, 'locus', 'loci')} showed concordant support between GLS and GQS."
    )


def generate_methods_posttree_signal_fclm(
    params: dict[str, Any],
    key_results: dict[str, Any],
    tool_versions: dict[str, Any],
) -> str:
    n_taxsets = key_results.get("n_taxsets", 4)
    model = params.get("model_expr") or "a substitution model"
    iqtree_ver = tool_versions.get("iqtree3", "IQ-TREE3")
    lmap = params.get("lmap") or "50 × n_taxa"
    return (
        f"Four-cluster Likelihood Mapping (FcLM) was performed using {iqtree_ver} "
        f"({model}) across {_describe_n(n_taxsets, 'taxon cluster', 'taxon clusters')} "
        f"with {lmap} quartets sampled."
    )
```

Then update the `METHODS_GENERATORS` registry — replace:

```python
    "posttree.signal": generate_methods_posttree_signal,
```

with:

```python
    "posttree.signal.lnl": generate_methods_posttree_signal_lnl,
    "posttree.signal.consistent": generate_methods_posttree_signal_consistent,
    "posttree.signal.fclm": generate_methods_posttree_signal_fclm,
```

- [ ] **Step 3: Verify report integration**

```bash
python -c "
from phyloai.report.collector import STEP_ORDER, parse_step_id
print('posttree.signal.lnl' in STEP_ORDER)
cmd = 'phyloai posttree signal lnl --matrix m.fa --candidate-trees t --model-expr LG+F+R4'
sid = parse_step_id(cmd)
print(sid)
"
```

Expected: `True` and `posttree.signal.lnl`

- [ ] **Step 4: Verify templates**

```bash
python -c "
from phyloai.report.templates import METHODS_GENERATORS
assert 'posttree.signal.lnl' in METHODS_GENERATORS
assert 'posttree.signal.consistent' in METHODS_GENERATORS
assert 'posttree.signal.fclm' in METHODS_GENERATORS
assert 'posttree.signal' not in METHODS_GENERATORS
print('OK')
"
```

Expected: `OK`

---

## Task 7: Full test suite + result.json compliance check

**Files:**
- Modify: `tests/posttree/test_signal_lnl.py`
- Modify: `tests/posttree/test_signal_consistent.py`
- Modify: `tests/posttree/test_signal_fclm.py`

- [ ] **Step 1: Add result.json compliance tests**

```python
# Add to each test file:

class TestResultJsonCompliance:
    """Verify result.json conforms to JSON output standard."""

    def test_lnl_dry_run_params_complete(self, tmp_path: Path) -> None:
        from phyloai.posttree.signal import run_signal_lnl
        matrix = tmp_path / "m.fa"
        matrix.write_text(">A\nMKT\n>B\nMKA\n")
        trees = tmp_path / "trees"
        trees.write_text("(A,B);\n(B,A);\n")
        result = run_signal_lnl(
            matrix=matrix, candidate_trees=[trees],
            model_expr="LG+F+R4", partitions=None, locus_ranges=None,
            guide_tree=None, threads="auto", iqtree_path=None,
            tool_args=None, metrics=None,
            output_dir=tmp_path / "out", overwrite=False,
            dry_run=True, quiet=True,
        )
        required_params = {
            "matrix", "candidate_trees_raw", "model_expr", "partitions",
            "locus_ranges", "guide_tree", "threads", "iqtree_path",
            "tool_args", "metrics", "output_dir", "overwrite", "dry_run", "quiet",
        }
        assert required_params <= set(result["params"].keys()), \
            f"Missing params: {required_params - set(result['params'].keys())}"
        assert result["command"].startswith("phyloai posttree signal lnl")
        assert isinstance(result["data"]["cmd"], list)
        assert isinstance(result["data"]["tool_stderr"], str)
        assert "wall_time" in result
        assert "tool_versions" in result
        assert "key_results" in result
        assert result["error"] is None
```

- [ ] **Step 2: Run full test suite**

```bash
python -m pytest tests/posttree/test_signal_lnl.py tests/posttree/test_signal_consistent.py tests/posttree/test_signal_fclm.py -v
```

Expected: all unit/validation tests PASS; integration tests PASS or skip if tools unavailable.

- [ ] **Step 3: Run full posttree test suite to confirm no regressions**

```bash
python -m pytest tests/posttree/ -v
```

Expected: all existing tests still PASS.

---

## Self-Review

Spec coverage check:

| Spec requirement | Task |
|-----------------|------|
| `signal lnl` ΔSLS 2-tree signed / >2-tree pairwise | Task 1 `_delta_score` |
| `signal lnl` site_lnl.csv + gene_lnl.csv | Task 2 `run_signal_lnl` |
| `signal lnl` support_sig column (2-tree only) | Task 2 |
| `signal lnl` outlier identification (IQR whisker) | Task 1 `_outlier_loci` + Task 2 |
| `signal lnl` --metrics comparison plots | Task 1 `_compare_groups` + Task 2 |
| `signal consistent` exactly 2 trees validation | Task 3 `_validate_consistent_inputs` |
| `signal consistent` locus↔gene_tree matching (global policy) | Task 3 `scan_tree_dir` + matching check |
| `signal consistent` extra-taxa hard error | Task 3 `_compute_gqs_for_locus` |
| `signal consistent` post-prune < 4 taxa skip | Task 3 `_compute_gqs_for_locus` |
| `signal consistent` GQS Score: extraction | Task 3 `_run_wastral_gqs` |
| `signal consistent` float tolerance 1e-9 | Task 1 `_FLOAT_TOL` |
| `signal consistent` gqs.csv status/reason columns | Task 3 |
| `signal fclm` taxset CSV validation | Task 4 `_validate_fclm_inputs` |
| `signal fclm` CSV → NEXUS conversion | Task 4 `_csv_to_nexus` |
| `signal fclm` lmap default 50*n_taxa | Task 4 `run_signal_fclm` |
| CLI --partition-mode p\|Q | Task 5 `consistent_command` |
| CLI --locus-ranges separate from --partitions | Task 5 `lnl_command` |
| result.json Single Pattern compliance | Task 7 compliance tests |
| collector STEP_ORDER + _THIRD_LEVEL | Task 6 |
| templates 3 generators + registry | Task 6 |
| output_dir non-empty conflict → error | Task 2/3/4 `run_signal_*` |
| wastral non-zero exit → error result exit code 2 | Task 3 `_run_wastral_gqs` + RuntimeError handler |
| T1/T2 taxon set mismatch → hard error | Task 3 `run_signal_consistent` |
| `--metrics` validation + empty group_a NA output | Task 1 `_compare_groups` |
| MCP 3 tools replacing stub | Task 8 |
| README signal stub replacement | Task 8 |
| `docs/commands/posttree-signal.md` + `.zh.md` | Task 8 |
| workflow skill update | Task 8 |
| FcLM taxa extraction supports non-FASTA matrix formats | Task 4 `_read_matrix_taxa` via `FormatConverter` + PHYLIP regression test |
| FcLM `n_quartets` in `key_results` | Task 4 `run_signal_fclm` (dry_run + success paths) |
| `scan_tree_dir` ValueError → structured error result | Task 3 `run_signal_consistent` try/except |
| T1/T2 taxon mismatch checked before IQ-TREE | Task 3 `run_signal_consistent` fail-fast |
| `_parse_candidate_trees` unused `output_dir` removed | Task 5 |
| MCP signal stub replaced by generated lnl/consistent/fclm tools | Task 8 MCP stub removal + regression test |

---

## Task 8: Associated updates (MCP, README, docs, skill)

**Files:**
- Modify: `skills/phyloai-workflow/SKILL.md`
- Modify: `README.md`
- Create: `docs/commands/posttree-signal.md`
- Create: `docs/commands/posttree-signal.zh.md`

- [ ] **Step 1: Update `skills/phyloai-workflow/SKILL.md`**

Find the section that describes `posttree topology` usage and add a parallel section for signal. Add after the topology entry:

```markdown
### posttree signal lnl
**Purpose:** Site-wise and gene-wise log-likelihood score distribution.
**Required:** `--matrix`, `--candidate-trees`, one of `--model-expr`/`--partitions`/`-m` in `--tool-args`
**Gene-wise:** requires `--partitions` (also passed to IQ-TREE -p) or `--locus-ranges` (boundary-only, compatible with `--model-expr`)
**`--partitions` and `--locus-ranges` are mutually exclusive**
**Output:** `site_lnl.csv`, `site_support.pdf`; if locus ranges: `gene_lnl.csv`, `gene_support.pdf`, `outlier_genes.txt`; if `--metrics`: `outlier_comparison.csv/pdf` and, with at least two non-ambiguous support groups, `support_comparison.csv/pdf`

### posttree signal consistent
**Purpose:** Consistent gene identification (Shen et al. 2021) via GLS + GQS.
**Required:** `--matrix`, `--candidate-trees` (exactly 2), `--tree-dir`, one of `--partitions`/`--locus-ranges` for GLS
**`--partition-mode p|Q`:** only valid with `--partitions`; `p`=edge-linked proportional, `Q`=edge-unlinked
**Output:** `gls.csv`, `gqs.csv`, `consistent_genes.txt`, `inconsistent_genes.txt`, support PDFs; if `--metrics`: `consistent_comparison.csv/pdf`

### posttree signal fclm
**Purpose:** Four-cluster Likelihood Mapping.
**Required:** `--matrix`, `--taxset-csv` (taxon,taxset; min 4 taxsets; mutually exclusive; all matrix taxa assigned), one of `--model-expr`/`-m` in `--tool-args`
**Output:** `cluster.nexus`, `<prefix>.lmap.eps`, `<prefix>.iqtree` (contains all lmap statistics)
```

- [ ] **Step 2: Update `README.md`**

Find the existing `phyloai posttree signal` stub line and replace with:

```bash
phyloai posttree signal lnl        --matrix ./matrix.fa --candidate-trees trees --model-expr LG+F+R4
phyloai posttree signal consistent --matrix ./matrix.fa --candidate-trees T1.tre,T2.tre --tree-dir ./gene_trees --model-expr LG+F+R4 --locus-ranges partitions.txt
phyloai posttree signal fclm       --matrix ./matrix.fa --taxset-csv taxsets.csv --model-expr LG+C60+F+R4
```

- [ ] **Step 3: Create `docs/commands/posttree-signal.md`**

```markdown
# phyloai posttree signal

Phylogenetic signal distribution analysis. Three independent subcommands sharing the same model parameter interface as `posttree topology`.

## Subcommands

- [`lnl`](#lnl) — Site-wise and gene-wise log-likelihood score distribution
- [`consistent`](#consistent) — Consistent gene identification (Shen et al. 2021)
- [`fclm`](#fclm) — Four-cluster Likelihood Mapping

---

## lnl

### Purpose
Computes site-wise and optionally gene-wise log-likelihood scores across candidate topologies using IQ-TREE3 `-wslr` (Shen et al. 2017).

### Usage
```bash
phyloai posttree signal lnl \
  --matrix matrix.aa.fa \
  --candidate-trees trees \
  --model-expr LG+F+R4 \
  [--locus-ranges partitions.txt] \
  [--threads auto] \
  [--output-dir runs/posttree/signal/lnl]
```

### Inputs
| Flag | Required | Description |
|------|----------|-------------|
| `--matrix` | Yes | Supermatrix alignment |
| `--candidate-trees` | Yes | Tree-list file or comma-separated NEWICK files |
| `--model-expr` | One of | IQ-TREE model expression; mutually exclusive with `--partitions` |
| `--partitions` | One of | Partition file → IQ-TREE `-p`/`-Q` (per `--partition-mode`); also extracts locus boundaries; mutually exclusive with `--locus-ranges` |
| `--partition-mode` | No | `p` = `-p` (edge-linked), `Q` = `-Q` (edge-unlinked); default `p`; only valid with `--partitions` |
| `--locus-ranges` | No | Partition file for locus boundary extraction only (not passed to IQ-TREE) |
| `--metrics` | No | Metrics CSV from `pretree metrics` for outlier-vs-nonoutlier and tree-support-group comparisons |

### Outputs
- `site_lnl.csv` — site-wise lnL, ΔSLS, support; ΔSLS = lnL_T1−lnL_T2 (2 trees) or mean pairwise |diff| (>2 trees)
- `site_support.pdf` — bar chart of site support distribution
- `gene_lnl.csv`, `gene_support.pdf`, `outlier_genes.txt` — if locus ranges provided
- `outlier_comparison.csv/pdf` — if `--metrics` provided
- `support_comparison.csv/pdf` — if `--metrics` provided and at least two non-ambiguous support groups exist

### Notes
- `support_sig` column (|ΔGLS| ≥ 2) only appears in `gene_lnl.csv` for 2-tree comparisons
- Outlier genes defined by Tukey 1.5×IQR on |ΔGLS| (Shen et al. 2017 eq. 3/4)

### Examples
```bash
# Site-wise only
phyloai posttree signal lnl --matrix m.fa --candidate-trees trees --model-expr LG+F+R4

# With gene-wise output
phyloai posttree signal lnl --matrix m.fa --candidate-trees trees --model-expr LG+F+R4 --locus-ranges partitions.txt

# With outlier and tree-support-group comparisons
phyloai posttree signal lnl --matrix m.fa --candidate-trees trees --model-expr LG+F+R4 --locus-ranges partitions.txt --metrics metrics.csv
```

---

## consistent

### Purpose
Identifies consistent genes where both likelihood-based (GLS) and quartet-based (GQS) signal agree, following Shen et al. (2021).

### Usage
```bash
phyloai posttree signal consistent \
  --matrix matrix.aa.fa \
  --candidate-trees T1.tre,T2.tre \
  --tree-dir gene_trees/ \
  --model-expr LG+F+R4 \
  --locus-ranges partitions.txt
```

### Inputs
| Flag | Required | Description |
|------|----------|-------------|
| `--matrix` | Yes | Supermatrix alignment |
| `--candidate-trees` | Yes | Exactly 2 candidate trees |
| `--tree-dir` | Yes | Directory of gene tree files (logical locus name = filename with 1-2 dot segments removed) |
| `--partitions` or `--locus-ranges` | Yes | Locus boundary source (one required) |
| `--partition-mode p\|Q` | No | Only with `--partitions`: `p`=edge-linked (default), `Q`=edge-unlinked |

### Outputs
- `gls.csv` — GLS per locus
- `gqs.csv` — GQS per locus (with `status`/`reason` columns for skipped loci)
- `consistent_genes.txt`, `inconsistent_genes.txt`
- `gls_support.pdf`, `gqs_support.pdf`
- `consistent_comparison.csv/pdf` — if `--metrics` provided

### Warnings/Errors
- Error if `--candidate-trees` ≠ 2 trees
- Error if gene tree contains taxa not in reference trees
- Error if locus set in partition file ≠ gene tree file set
- Skipped (not error) if pruned reference tree has < 4 taxa

---

## fclm

### Purpose
Four-cluster Likelihood Mapping using IQ-TREE3 `-lmap -lmclust`.

### Usage
```bash
phyloai posttree signal fclm \
  --matrix matrix.aa.fa \
  --taxset-csv taxsets.csv \
  --model-expr LG+C60+F+R4
```

### Inputs
| Flag | Required | Description |
|------|----------|-------------|
| `--matrix` | Yes | Supermatrix alignment |
| `--taxset-csv` | Yes | Two-column CSV (`taxon,taxset`); min 4 taxsets; each taxon in exactly one taxset |
| `--model-expr` | One of | Model expression; mutually exclusive with `--partitions` |
| `--partitions` | One of | Partition file (e.g. `.best_model.nex`); passed to IQ-TREE as `-p`/`-Q` (per `--partition-mode`) |
| `--partition-mode` | No | `p` = `-p` (edge-linked), `Q` = `-Q` (edge-unlinked); default `p`; only valid with `--partitions` |
| `--lmap` | No | Quartet count: `ALL`, integer, or omit for `50 × n_taxa` |

### Outputs
- `cluster.nexus` — generated NEXUS cluster file
- `<prefix>.lmap.eps` — IQ-TREE likelihood mapping figure
- `<prefix>.iqtree` — full IQ-TREE report containing all lmap statistics

### Notes
- All lmap statistics (per-taxon and overall quartet resolution) are in `<prefix>.iqtree`
- `--model-expr` and `--partitions` are mutually exclusive; `--partition-mode` (default `p`) controls the IQ-TREE flag
```

- [ ] **Step 4: Create `docs/commands/posttree-signal.zh.md`**

Create Chinese translation of the above document following the same structure as existing `posttree-topology.zh.md`.

- [ ] **Step 5: Verify docs/commands directory**

```bash
ls /Users/zf/data/coding/phyloAI-syserror/docs/commands/ | grep signal
```

Expected: `posttree-signal.md` and `posttree-signal.zh.md`

- [ ] **Step 6: Replace the MCP stub with auto-generated signal tools**

MCP tool definitions are generated from Click commands by
`phyloai.mcp.schema_gen.walk_click_tree`, so do not write bespoke MCP handlers.
After Task 5 adds the three Click commands, edit `phyloai/mcp/tools/stubs.py`:

```python
STUB_TOOL_NAMES: frozenset[str] = frozenset(
    {
        # Remove "posttree_signal"; posttree_signal_lnl, posttree_signal_consistent,
        # and posttree_signal_fclm now come from the Click tree.
        "posttree_simulate",
        "posttree_syserror_brlen",
        "posttree_syserror_cca",
        "posttree_syserror_sites",
    }
)

_DESCRIPTIONS = {
    # Remove the "posttree_signal" entry.
    "posttree_simulate": "AliSim simulation and gene-jackknife resampling (not yet available).",
    "posttree_syserror_brlen": "Systematic error diagnosis: branch-length screen (not yet available).",
    "posttree_syserror_cca": "Systematic error diagnosis: composition analysis (not yet available).",
    "posttree_syserror_sites": "Systematic error diagnosis: site-wise analysis (not yet available).",
}
```

- [ ] **Step 7: Add MCP tool-generation regression test**

Add to the existing `tests/mcp/test_cli_tools.py` (file already exists at `tests/mcp/test_cli_tools.py`):

```python
def test_signal_commands_replace_signal_stub() -> None:
    from phyloai.mcp.tools.cli_tools import get_tool_definitions

    definitions = get_tool_definitions()
    assert "posttree_signal_lnl" in definitions
    assert "posttree_signal_consistent" in definitions
    assert "posttree_signal_fclm" in definitions
    assert "posttree_signal" not in definitions

    assert "matrix" in definitions["posttree_signal_lnl"]["inputSchema"]["properties"]
    assert "tree_dir" in definitions["posttree_signal_consistent"]["inputSchema"]["properties"]
    assert "taxset_csv" in definitions["posttree_signal_fclm"]["inputSchema"]["properties"]
```

- [ ] **Step 8: Run MCP regression test**

```bash
python -m pytest tests/mcp/test_cli_tools.py -v
```

Expected: PASS; MCP definitions include exactly the three signal subcommands and no legacy `posttree_signal` stub.
