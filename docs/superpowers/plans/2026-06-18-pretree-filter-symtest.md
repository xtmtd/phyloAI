# phyloai pretree filter symtest Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a `symtest` subcommand to `phyloai pretree filter` that runs IQ-TREE3's `--symtest-only` to test phylogenetic symmetry assumptions, filters loci by p-value, and copies retained MSAs (and optionally trees).

**Architecture:** Single IQ-TREE invocation on a temporary supermatrix + partition file (built by reusing concat internals). Parses `.symtest.csv` to select per-locus p-values from `SymPval` (default), `MarPval`, or `IntPval`. Retains loci with `p >= --symtest-pval`. No `--resume` (single deterministic invocation).

**Tech Stack:** Python stdlib (csv, json, shutil, tempfile, pathlib), BioPython (SeqIO), Click, Rich, existing phyloai core/runner/pretree modules.

---

## File Structure

| File | Action | Responsibility |
|------|--------|----------------|
| `phyloai/pretree/filter.py` | Modify | Add `run_symtest()` core function |
| `phyloai/cli/commands/pretree.py` | Modify | Add `filter_symtest_command` CLI handler + register in `_FilterGroup` |
| `tests/pretree/test_filter.py` | Modify | Add unit tests for `.symtest.csv` parsing, filtering logic |
| `tests/cli/test_pretree_filter.py` | Create | Add CLI tests for symtest subcommand |
| `docs/commands/pretree-filter.md` | Modify | Add symtest section |
| `docs/superpowers/specs/2026-06-15-phyloai-pretree-filter-design.md` | Modify | Update subcommand list, output dirs |
| `docs/superpowers/specs/2026-06-07-phyloai-design.md` | Modify | Add symtest to CLI examples |
| `README.md` | Modify | Add symtest to filter subcommand list |

---

### Task 1: Core library function `run_symtest` in filter.py

**Files:**
- Modify: `phyloai/pretree/filter.py` (append ~180 lines after cluster section, before end of file)

- [ ] **Step 1: Write the failing tests**

Create `tests/pretree/test_symtest.py` (or add to `tests/pretree/test_filter.py`):

```python
"""Tests for symtest filtering in filter.py."""

import csv
import io
import tempfile
from pathlib import Path

import pytest

from phyloai.pretree.filter import _parse_symtest_csv, _build_symtest_supermatrix


# --- _parse_symtest_csv ---

def test_parse_symtest_csv_all_columns():
    csv_content = (
        "Name,SymSig,SymNon,SymPval,MarSig,MarNon,MarPval,IntSig,IntNon,IntPval\n"
        "gene1,44,92,0.475,50,86,0.722,4,132,0.239\n"
        "gene2,43,93,0.142,49,87,0.205,5,131,0.170\n"
        "gene3,53,83,0.005,58,78,0.002,6,130,0.343\n"
    )
    fp = io.StringIO(csv_content)
    results = _parse_symtest_csv(fp)
    assert len(results) == 3
    assert results[0]["Name"] == "gene1"
    assert results[0]["SymPval"] == 0.475
    assert results[2]["SymPval"] == 0.005
    assert results[1]["MarPval"] == 0.205
    assert results[1]["IntPval"] == 0.170


def test_parse_symtest_csv_skips_comment_lines():
    csv_content = (
        "# Matched-pair tests of symmetry\n"
        "# comment\n"
        "Name,SymSig,SymNon,SymPval,MarSig,MarNon,MarPval,IntSig,IntNon,IntPval\n"
        "gene1,44,92,0.475,50,86,0.722,4,132,0.239\n"
    )
    fp = io.StringIO(csv_content)
    results = _parse_symtest_csv(fp)
    assert len(results) == 1
    assert results[0]["Name"] == "gene1"


def test_parse_symtest_csv_empty():
    csv_content = (
        "# comment only\n"
        "Name,SymSig,SymNon,SymPval,MarSig,MarNon,MarPval,IntSig,IntNon,IntPval\n"
    )
    fp = io.StringIO(csv_content)
    results = _parse_symtest_csv(fp)
    assert results == []


def test_parse_symtest_csv_missing_header_raises():
    csv_content = "bad,header,here\nval1,val2,val3\n"
    fp = io.StringIO(csv_content)
    with pytest.raises(ValueError, match="missing expected columns"):
        _parse_symtest_csv(fp)


def test_parse_symtest_csv_non_numeric_pval():
    csv_content = (
        "Name,SymSig,SymNon,SymPval,MarSig,MarNon,MarPval,IntSig,IntNon,IntPval\n"
        "gene1,44,92,NA,50,86,NA,4,132,NA\n"
    )
    fp = io.StringIO(csv_content)
    results = _parse_symtest_csv(fp)
    assert results[0]["SymPval"] is None


# --- _filter_by_symtest_pval ---

def test_filter_symtest_retain_above_threshold():
    from phyloai.pretree.filter import _filter_by_symtest_pval
    results = [
        {"Name": "gene1", "SymPval": 0.475, "SymSig": 44, "SymNon": 92,
         "MarSig": 50, "MarNon": 86, "MarPval": 0.722,
         "IntSig": 4, "IntNon": 132, "IntPval": 0.239},
        {"Name": "gene2", "SymPval": 0.005, "SymSig": 53, "SymNon": 83,
         "MarSig": 58, "MarNon": 78, "MarPval": 0.002,
         "IntSig": 6, "IntNon": 130, "IntPval": 0.343},
    ]
    retained, dropped, decisions = _filter_by_symtest_pval(results, "Sym", 0.05)
    assert len(retained) == 1
    assert retained[0]["locus"] == "gene1"
    assert len(dropped) == 1
    assert dropped[0]["locus"] == "gene2"
    assert decisions[0]["status"] == "retained"
    assert decisions[1]["status"] == "dropped"


def test_filter_symtest_uses_mar_column():
    from phyloai.pretree.filter import _filter_by_symtest_pval
    results = [
        {"Name": "gene1", "SymPval": 0.001, "SymSig": 44, "SymNon": 92,
         "MarSig": 50, "MarNon": 86, "MarPval": 0.722,
         "IntSig": 4, "IntNon": 132, "IntPval": 0.239},
    ]
    retained, dropped, decisions = _filter_by_symtest_pval(results, "MAR", 0.05)
    assert len(retained) == 1  # MarPval=0.722 >= 0.05
    assert retained[0]["p_value"] == 0.722


def test_filter_symtest_none_pval_dropped():
    from phyloai.pretree.filter import _filter_by_symtest_pval
    results = [
        {"Name": "gene1", "SymPval": None, "SymSig": 44, "SymNon": 92,
         "MarSig": 50, "MarNon": 86, "MarPval": None,
         "IntSig": 4, "IntNon": 132, "IntPval": None},
    ]
    retained, dropped, decisions = _filter_by_symtest_pval(results, "Sym", 0.05)
    assert len(retained) == 0
    assert dropped[0]["reason"] == "p_value is null"


# --- _build_symtest_supermatrix ---

def test_build_symtest_supermatrix(tmp_path):
    msa1 = tmp_path / "gene1.fa"
    msa1.write_text(">taxa1\nACGT\n>taxa2\nACGT\n")
    msa2 = tmp_path / "gene2.fa"
    msa2.write_text(">taxa1\nTTTT\n>taxa2\nGGGG\n>taxa3\nCCCC\n")
    msa_map = {"gene1": msa1, "gene2": msa2}

    matrix_str, genes, prefix_type = _build_symtest_supermatrix(msa_map)

    assert prefix_type in ("DNA", "LG")
    assert len(genes) == 2
    assert genes[0][0] == "gene1"
    assert genes[0][1] == 1
    assert genes[0][2] == 4
    assert genes[1][0] == "gene2"
    assert genes[1][1] == 5
    assert genes[1][2] == 8
    assert ">taxa1" in matrix_str
    assert ">taxa2" in matrix_str
    assert ">taxa3" in matrix_str


def test_build_symtest_supermatrix_empty_dir():
    with pytest.raises(ValueError, match="No valid MSA files"):
        _build_symtest_supermatrix({})
```

- [ ] **Step 2: Run tests to verify they fail**

```bash
cd /Users/zf/data/coding/phyloAI && python -m pytest tests/pretree/test_symtest.py -v 2>&1 | head -30
```

Expected: All tests FAIL with `ImportError` (functions not defined).

- [ ] **Step 3: Update filter.py imports and add helper functions**

At the top of `filter.py`, add `import io` alongside existing `import csv` (line 5). Also update the `file_matching` import (line ~31) to include `scan_tree_dir`:

```python
# Before (line ~5):
import csv

# After:
import csv
import io

# Before (line ~31):
from phyloai.core.file_matching import (
    logical_msa_locus_name,
    pair_msa_and_tree_maps,
    scan_msa_dir,
    scan_tree_dir,    # <-- ADD THIS
)
```

Append the following code at the end of `filter.py`:

```python
# --- Symmetry test (symtest) ---

import io  # added to top of file alongside existing `import csv`

_EXPECTED_SYMTEST_COLUMNS = {"Name", "SymSig", "SymNon", "SymPval",
                              "MarSig", "MarNon", "MarPval",
                              "IntSig", "IntNon", "IntPval"}


def _parse_symtest_csv(fileobj) -> list[dict[str, Any]]:
    """Parse IQ-TREE ``.symtest.csv`` output into a list of per-partition dicts.

    Skips comment lines (starting with ``#``).  P-value columns are
    parsed as float; ``NA`` or unparseable values become None.
    """
    lines = [line for line in fileobj if not line.startswith("#")]
    if not lines:
        return []

    reader = csv.DictReader(io.StringIO("".join(lines)))
    if not reader.fieldnames:
        raise ValueError("Empty CSV header in symtest output")
    missing = _EXPECTED_SYMTEST_COLUMNS - set(reader.fieldnames)
    if missing:
        raise ValueError(
            f"Symtest CSV missing expected columns: {', '.join(sorted(missing))}"
        )

    results: list[dict[str, Any]] = []
    for row in reader:
        entry: dict[str, Any] = {}
        for key, value in row.items():
            if key in ("SymPval", "MarPval", "IntPval"):
                try:
                    entry[key] = float(value)
                except (ValueError, TypeError):
                    entry[key] = None
            elif key in ("SymSig", "SymNon", "MarSig", "MarNon", "IntSig", "IntNon"):
                try:
                    entry[key] = int(value)
                except (ValueError, TypeError):
                    entry[key] = 0
            else:
                entry[key] = value
        results.append(entry)
    return results


def _filter_by_symtest_pval(
    results: list[dict[str, Any]],
    symtest_type: str,
    threshold: float,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]]]:
    """Filter parsed symtest results by the selected p-value column.

    symtest_type is one of ``"Sym"``, ``"MAR"``, ``"INT"``.
    """
    pval_col = {"Sym": "SymPval", "MAR": "MarPval", "INT": "IntPval"}[symtest_type]

    retained: list[dict[str, Any]] = []
    dropped: list[dict[str, Any]] = []
    decisions: list[dict[str, Any]] = []

    for entry in results:
        locus = entry.get("Name", "")
        p_value = entry.get(pval_col)

        decision = {
            "locus": locus,
            "p_value": p_value,
            "symtest_type": symtest_type,
            "sym_sig": entry.get("SymSig", 0),
            "sym_non": entry.get("SymNon", 0),
            "mar_sig": entry.get("MarSig", 0),
            "mar_non": entry.get("MarNon", 0),
            "int_sig": entry.get("IntSig", 0),
            "int_non": entry.get("IntNon", 0),
        }

        if p_value is None:
            decision["status"] = "dropped"
            dropped.append({"locus": locus, "reason": "p_value is null"})
            decisions.append(decision)
        elif p_value >= threshold:
            decision["status"] = "retained"
            retained.append({"locus": locus})
            decisions.append(decision)
        else:
            decision["status"] = "dropped"
            dropped.append({"locus": locus, "reason": f"{pval_col}={p_value} < {threshold}"})
            decisions.append(decision)

    return retained, dropped, decisions


def _build_symtest_supermatrix(
    msa_map: dict[str, Path],
) -> tuple[str, list[tuple[str, int, int]], str]:
    """Build a supermatrix string and partition list from a dict of MSAs.

    Returns ``(matrix_fasta_str, genes, prefix_type)`` where *genes* is
    ``[(name, start1, end1), ...]`` with 1-based positions.  Uses
    ``_read_msa`` from concat.py for format-agnostic reading.
    """
    from phyloai.pretree.concat import _read_msa

    if not msa_map:
        raise ValueError("No valid MSA files found")

    all_taxa: set[str] = set()
    msa_records: dict[str, tuple[list[str], list[str], int]] = {}

    for locus, path in sorted(msa_map.items()):
        taxa, seqs, length = _read_msa(path)
        if not taxa:
            continue
        all_taxa.update(taxa)
        msa_records[locus] = (taxa, seqs, length)

    if not msa_records:
        raise ValueError("No valid MSA files found")

    # Auto-detect seq_type from first 3 loci
    sample_seqs: list[str] = []
    for locus in list(msa_records.keys())[:3]:
        _, seqs, _ = msa_records[locus]
        sample_seqs.extend(seqs[:10])
    from phyloai.core.sequence_normalization import detect_seq_type
    seq_type = detect_seq_type(sample_seqs)

    if seq_type == "other":
        raise ValueError(
            f"Could not determine sequence type from MSA files. "
            f"Detected type: 'other'. Ensure input files contain "
            f"valid AA or NT sequences."
        )

    prefix_type = "DNA" if seq_type in ("NT", "CODON") else "LG"

    # Build supermatrix
    matrix_parts: dict[str, list[str]] = {taxon: [] for taxon in all_taxa}
    genes: list[tuple[str, int, int]] = []
    pos = 1

    for locus, (taxa, seqs, length) in sorted(msa_records.items()):
        genes.append((locus, pos, pos + length - 1))
        pos += length
        taxon_to_seq = dict(zip(taxa, seqs))
        for taxon in all_taxa:
            seq = taxon_to_seq.get(taxon, "?" * length)
            matrix_parts[taxon].append(seq)

    taxon_order = sorted(all_taxa)
    lines: list[str] = []
    for taxon in taxon_order:
        seq = "".join(matrix_parts[taxon])
        # 60-char line wrapping per design §9.11
        wrapped = "\n".join(seq[i:i + 60] for i in range(0, len(seq), 60))
        lines.append(f">{taxon}\n{wrapped}")

    matrix_str = "\n".join(lines) + "\n"
    return matrix_str, genes, prefix_type
```

- [ ] **Step 4: Run tests to verify helper functions pass**

```bash
cd /Users/zf/data/coding/phyloAI && python -m pytest tests/pretree/test_symtest.py -v
```

Expected: All tests PASS.

- [ ] **Step 5: Add `run_symtest` main function to filter.py**

Append after the helper functions added in Step 3:

```python
def run_symtest(
    msa_dir: Path, output_dir: Path, *,
    symtest_type: str | None = None,
    symtest_pval: float = 0.05,
    symtest_keep_zero: bool = False,
    iqtree_path: Path | None = None,
    threads: int = 4,
    tree_dir: Path | None = None,
    msa_map: dict[str, Path] | None = None,
    table_format: str = "csv",
    dry_run: bool = False,
    overwrite: bool = False,
    quiet: bool = False,
) -> dict[str, Any]:
    """Run IQ-TREE symmetry test on all MSAs and filter by p-value.

    *msa_map* is an optional pre-scanned ``{locus: path}`` dict; when not
    provided it is built via :func:`scan_msa_dir`.
    """
    start = time.monotonic()
    env = ToolEnv()
    iqtree_exe = str(iqtree_path) if iqtree_path else str(env.require("iqtree3"))

    msa_map = scan_msa_dir(msa_dir) if msa_map is None else msa_map
    if not msa_map:
        raise ValueError(f"No valid MSA files found in {msa_dir}")

    # Resolve symtest_type: None -> "Sym"
    resolved_type = symtest_type if symtest_type else "Sym"

    command = f"phyloai pretree filter symtest --msa-dir {msa_dir} --symtest-pval {symtest_pval}"
    if symtest_type:
        command += f" --symtest-type {symtest_type}"

    params = {
        "msa_dir": str(msa_dir), "symtest_type": symtest_type,
        "symtest_pval": symtest_pval, "symtest_keep_zero": symtest_keep_zero,
        "threads": threads, "tree_dir": str(tree_dir) if tree_dir else None,
        "table_format": table_format,
    }

    if dry_run:
        return {
            "status": "success", "command": command, "wall_time": 0,
            "tool_versions": {"iqtree3": "unknown"}, "params": params,
            "key_results": {"n_input": len(msa_map)},
            "error": None,
            "data": {"dry_run_cmd": f"{iqtree_exe} -s <matrix> -p <partitions> "
                     f"--symtest-only {'--symtest-type ' + symtest_type if symtest_type else ''} "
                     f"-T {threads}"},
        }

    _common_output_conflict(output_dir, overwrite)

    # Build supermatrix + partition files in temp dir
    matrix_str, genes, prefix_type = _build_symtest_supermatrix(msa_map)

    work_dir = Path(tempfile.mkdtemp(prefix="symtest_"))
    try:
        matrix_path = work_dir / "symtest_matrix.fa"
        partitions_path = work_dir / "symtest_partitions.txt"
        matrix_path.write_text(matrix_str)
        from phyloai.pretree.concat import _write_partitions
        _write_partitions(partitions_path, genes, prefix_type)

        # Build IQ-TREE command (--symtest-pval NOT passed; used Python-side only)
        cmd = [
            iqtree_exe,
            "-s", str(matrix_path),
            "-p", str(partitions_path),
            "--symtest-only",
        ]
        if symtest_type:
            cmd.extend(["--symtest-type", symtest_type])
        if symtest_keep_zero:
            cmd.append("--symtest-keep-zero")
        if threads > 1:
            cmd.extend(["-T", str(threads)])

        # Run IQ-TREE
        runner = Runner()
        result = runner.run(cmd, tool_name="iqtree3", cwd=work_dir)

        if result.returncode != 0:
            raise RuntimeError(
                f"iqtree3 exited with code {result.returncode}.\n"
                f"STDERR:\n{result.stderr}"
            )

        # Parse symtest output
        symtest_csv = work_dir / "symtest_partitions.txt.symtest.csv"
        if not symtest_csv.exists():
            raise RuntimeError(
                f"Expected symtest output not found: {symtest_csv}\n"
                f"STDERR:\n{result.stderr}"
            )

        with open(symtest_csv) as fh:
            symtest_results = _parse_symtest_csv(fh)

        if not symtest_results:
            raise RuntimeError("Symtest CSV is empty -- no partitions parsed.")

        # Cross-validate CSV names against MSA map
        csv_names = {r["Name"] for r in symtest_results}
        msa_names = set(msa_map.keys())
        missing_in_csv = msa_names - csv_names
        extra_in_csv = csv_names - msa_names
        if missing_in_csv:
            raise RuntimeError(
                f"Loci in MSA directory but missing from symtest CSV: "
                f"{', '.join(sorted(missing_in_csv))}. "
                f"IQ-TREE may have dropped these partitions."
            )
        if extra_in_csv:
            # Warn but don't fail -- IQ-TREE may have renamed partitions
            import warnings
            warnings.warn(
                f"Partition names in symtest CSV not found in MSA map: "
                f"{', '.join(sorted(extra_in_csv))}. These will be skipped."
            )
            symtest_results = [r for r in symtest_results if r["Name"] in msa_names]

        # Filter
        retained, dropped, decisions = _filter_by_symtest_pval(
            symtest_results, resolved_type, symtest_pval,
        )

        # Copy retained MSAs
        seqs_out = output_dir / "seqs"
        seqs_out.mkdir(parents=True, exist_ok=True)
        for r in retained:
            locus = r["locus"]
            if locus in msa_map:
                shutil.copy2(msa_map[locus], seqs_out / f"{locus}.fa")

        # Copy retained trees (if --tree-dir)
        retained_tree_count = 0
        missed_tree_count = 0
        if tree_dir:
            tree_map = scan_tree_dir(tree_dir)
            trees_out = output_dir / "trees"
            trees_out.mkdir(parents=True, exist_ok=True)
            retained_loci = {r["locus"] for r in retained}
            for locus in sorted(retained_loci):
                if locus in tree_map:
                    shutil.copy2(tree_map[locus], trees_out / tree_map[locus].name)
                    retained_tree_count += 1
                else:
                    missed_tree_count += 1

        # Write decision tables
        delimiter = _table_delimiter(table_format)
        suffix = _table_suffix(table_format)

        _write_csv_table(
            retained, output_dir / f"retained_loci{suffix}",
            ["locus"], delimiter,
        )
        _write_csv_table(
            dropped, output_dir / f"dropped_loci{suffix}",
            ["locus", "reason"], delimiter,
        )
        _write_csv_table(
            decisions, output_dir / f"filter_decisions{suffix}",
            ["locus", "status", "p_value", "symtest_type",
             "sym_sig", "sym_non", "mar_sig", "mar_non", "int_sig", "int_non"],
            delimiter,
        )

        # MSA stats
        retained_paths = [seqs_out / f"{r['locus']}.fa" for r in retained]
        msa_stats = _compute_retained_msa_stats(retained_paths)

        wall_time = time.monotonic() - start
        payload = {
            "status": "success",
            "command": command,
            "wall_time": round(wall_time, 2),
            "tool_versions": {"iqtree3": "unknown"},
            "params": params,
            "key_results": {
                "n_input": len(symtest_results),
                "n_retained": len(retained),
                "n_dropped": len(dropped),
                "p_value_threshold": symtest_pval,
                "symtest_type": resolved_type,
                "retained_trees_copied": retained_tree_count,
            },
            "error": None,
            "data": {
                "retained_msa_stats": msa_stats,
                "retained_loci": [r["locus"] for r in retained],
                "dropped_loci": dropped,
                "decisions": decisions,
                "retained_tree_count": retained_tree_count,
                "missed_tree_count": missed_tree_count,
                "skipped": [],
            },
        }
        _write_result_json(payload, output_dir)
        _write_filter_log(output_dir, command, wall_time,
                          payload["tool_versions"], True)

    finally:
        shutil.rmtree(work_dir, ignore_errors=True)

    return payload
```

- [ ] **Step 6: Run tests again to confirm they pass**

```bash
cd /Users/zf/data/coding/phyloAI && python -m pytest tests/pretree/test_symtest.py -v
```

Expected: All PASS.

- [ ] **Step 7: Commit**

```bash
git add phyloai/pretree/filter.py tests/pretree/test_symtest.py
git commit -m "feat(filter): add run_symtest core function with symtest CSV parsing and p-value filtering"
```

---

---

### Task 2: CLI subcommand `filter symtest` in pretree.py

**Files:**
- Modify: `phyloai/cli/commands/pretree.py` (update imports, _FilterGroup, add command + registration)

- [ ] **Step 1: Update the filter import line**

Change line 1366 from:
```python
from phyloai.pretree.filter import render_filter_summary_table, run_taper, run_treeshrink, run_metrics_filter, run_cluster_filter  # noqa: E402
```
to:
```python
from phyloai.pretree.filter import render_filter_summary_table, run_taper, run_treeshrink, run_metrics_filter, run_symtest, run_cluster_filter  # noqa: E402
```

- [ ] **Step 2: Update `_FilterGroup.list_commands`**

Change line 1371 from:
```python
        return ["taper", "treeshrink", "metrics", "cluster"]
```
to:
```python
        return ["taper", "treeshrink", "metrics", "symtest", "cluster"]
```

- [ ] **Step 3: Add symtest command (insert after `filter_metrics_command`, before `filter_cluster_command`)**

Insert the following ~100 lines after the `filter_metrics_command` function (before line ~1900, the `filter_cluster_command`):

```python
# ---- filter symtest ----

_SYMTEST_HELP = (
    "Test phylogenetic symmetry assumptions per locus using IQ-TREE3's "
    "--symtest-only, then filter loci by p-value.\n\n"
    "Workflow:\n"
    "  1. Scan --msa-dir for per-locus MSA files\n"
    "  2. Build a temporary supermatrix + partition file\n"
    "  3. Run 'iqtree -s <matrix> -p <partitions> --symtest-only'\n"
    "  4. Parse .symtest.csv for per-partition p-values\n"
    "  5. Retain loci with p >= --symtest-pval, drop those below\n"
    "  6. Copy retained MSAs to seqs/, optionally trees to trees/\n\n"
    "The p-value column used depends on --symtest-type:\n"
    "  (default) -> SymPval (combined stationarity + homogeneity)\n"
    "  MAR       -> MarPval (marginal / stationarity)\n"
    "  INT       -> IntPval (internal / homogeneity)\n\n"
    "References: Naser-Khdour et al. (2019) doi:10.1093/gbe/evz193"
)


def _validate_symtest_pval(ctx, param, value):
    if value <= 0 or value > 1:
        raise click.BadParameter("must be > 0 and <= 1")
    return value


@filter_group.command("symtest", help=_SYMTEST_HELP)
@click.option(
    "--msa-dir", type=click.Path(exists=True, file_okay=False, path_type=Path),
    required=True,
    help="Directory containing per-locus MSA files (any suffix).",
)
@click.option(
    "--symtest-type", type=click.Choice(["MAR", "INT"]),
    default=None,
    help="Which symmetry test to use for filtering.  When omitted (default), "
    "the combined Sym test is used (SymPval column).  MAR uses marginal "
    "(stationarity) test.  INT uses internal (homogeneity) test.",
)
@click.option(
    "--symtest-pval", type=float, default=0.05, show_default=True,
    callback=_validate_symtest_pval,
    help="P-value threshold.  Loci with p >= threshold are retained.",
)
@click.option(
    "--symtest-keep-zero", is_flag=True, default=False,
    help="Pass --symtest-keep-zero to IQ-TREE (keep NAs in the tests).",
)
@click.option(
    "--iqtree-path", type=click.Path(exists=True, dir_okay=False, path_type=Path),
    default=None,
    help="Explicit path to iqtree binary.  When omitted, resolved via "
    "PATH ('phyloai doctor' for detection status).",
)
@click.option(
    "--threads", "-t", type=int, default=4, show_default=True,
    help="Number of threads for IQ-TREE (-T).",
)
@click.option(
    "--tree-dir", type=click.Path(exists=True, file_okay=False, path_type=Path),
    default=None,
    help="Optional directory of gene tree files.  Trees matching retained "
    "loci (by logical locus name) are copied to trees/.",
)
@click.option(
    "--output-dir", "-o", type=click.Path(file_okay=False, path_type=Path),
    default=Path("runs/pretree/filter/symtest"), show_default=True,
    help="Directory for retained MSAs, optional trees, decision tables, "
    "result.json, and filter.log.",
)
@click.option(
    "--table-format", type=click.Choice(["csv", "tsv"]),
    default="csv", show_default=True,
    help="Delimiter and file suffix for auxiliary tables.",
)
@click.option(
    "--overwrite", is_flag=True, default=False,
    help="Delete and recreate --output-dir if it already exists.",
)
@click.option(
    "--dry-run", is_flag=True, default=False,
    help="Validate inputs and show the planned IQ-TREE command without writing files.",
)
@click.option(
    "--quiet", "-q", is_flag=True, default=False,
    help="Suppress all terminal output except errors.",
)
def filter_symtest_command(msa_dir, symtest_type, symtest_pval, symtest_keep_zero,
                           iqtree_path, threads, tree_dir, output_dir, table_format,
                           overwrite, dry_run, quiet):
    if threads < 1:
        _fail("--threads must be at least 1.", 1)

    if not quiet and not dry_run:
        from phyloai.core.file_matching import scan_msa_dir
        msa_map = scan_msa_dir(msa_dir)
        total = max(len(msa_map), 1)
        with Progress(console=console, transient=True) as progress:
            progress.add_task("IQ-TREE symmetry test running...", total=None)
            try:
                payload = run_symtest(
                    msa_dir=msa_dir, output_dir=output_dir,
                    symtest_type=symtest_type, symtest_pval=symtest_pval,
                    symtest_keep_zero=symtest_keep_zero,
                    iqtree_path=iqtree_path, threads=threads,
                    tree_dir=tree_dir, msa_map=msa_map,
                    table_format=table_format,
                    dry_run=dry_run, overwrite=overwrite, quiet=quiet,
                )
            except (ValueError, FileNotFoundError, RuntimeError) as exc:
                msg = str(exc)
                exit_code = 3 if "not found" in msg.lower() else (
                    2 if "exited with code" in msg.lower() else 1)
                _fail(msg, exit_code)
    else:
        try:
            payload = run_symtest(
                msa_dir=msa_dir, output_dir=output_dir,
                symtest_type=symtest_type, symtest_pval=symtest_pval,
                symtest_keep_zero=symtest_keep_zero,
                iqtree_path=iqtree_path, threads=threads,
                tree_dir=tree_dir, table_format=table_format,
                dry_run=dry_run, overwrite=overwrite, quiet=quiet,
            )
        except (ValueError, FileNotFoundError, RuntimeError) as exc:
            msg = str(exc)
            exit_code = 3 if "not found" in msg.lower() else (
                2 if "exited with code" in msg.lower() else 1)
            _fail(msg, exit_code)

    if dry_run:
        click.echo(f"Dry run: {payload['key_results']['n_input']} loci would be processed.")
        click.echo(payload["data"]["dry_run_cmd"])
        return

    if not quiet:
        console.print(render_filter_summary_table({
            "Input": payload["key_results"]["n_input"],
            "Retained": payload["key_results"]["n_retained"],
            "Dropped": payload["key_results"]["n_dropped"],
            "P-value threshold": payload["key_results"]["p_value_threshold"],
            "Symtest type": payload["key_results"]["symtest_type"],
        }))
        msa_stats = payload["data"].get("retained_msa_stats", {})
        if msa_stats and msa_stats.get("n_msa", 0) > 0:
            console.print(render_filter_summary_table({
                "Retained MSAs": msa_stats["n_msa"],
                "Total length": msa_stats["total_length"],
                "Mean length": msa_stats["mean_length"],
                "Min length": msa_stats["min_length"],
                "Max length": msa_stats["max_length"],
                "Mean taxa": msa_stats["mean_taxa"],
            }))
        if payload["key_results"].get("retained_trees_copied", 0) > 0:
            mt = payload["data"].get("missed_tree_count", 0)
            console.print(render_filter_summary_table({
                "Trees copied": payload["key_results"]["retained_trees_copied"],
                "Trees missed": mt,
            }))
        click.echo(f"Retained MSAs saved to {output_dir / 'seqs'}", err=True)
        if payload["key_results"].get("retained_trees_copied", 0) > 0:
            click.echo(f"Retained trees saved to {output_dir / 'trees'}", err=True)
        click.echo(f"Results saved to {output_dir / 'result.json'}", err=True)
```

- [ ] **Step 4: Add `_SYMTEST_HELP` to help text and update filter_group help text**

Update the filter_group help text at line ~1377 to mention symtest:
```python
@click.group(
    "filter",
    cls=_FilterGroup,
    help="TAPER site masking, TreeShrink taxa pruning, "
    "metric-rule loci filtering, symmetry test filtering, "
    "cluster-based exploration.",
)
```

- [ ] **Step 5: Verify CLI loads correctly**

```bash
cd /Users/zf/data/coding/phyloAI && python -m phyloai pretree filter symtest --help
```

Expected: Prints help text for symtest subcommand.

- [ ] **Step 6: Commit**

```bash
git add phyloai/cli/commands/pretree.py
git commit -m "feat(cli): add filter symtest subcommand with IQ-TREE symmetry test integration"
```

---

### Task 3: CLI integration tests

**Files:**
- Create: `tests/cli/test_pretree_filter_symtest.py`

- [ ] **Step 1: Write CLI tests**

```python
"""CLI tests for phyloai pretree filter symtest."""

from pathlib import Path

from click.testing import CliRunner

from phyloai.cli.main import main


def _make_msa_dir(tmp_path: Path) -> Path:
    """Create a minimal MSA directory for testing."""
    msa_dir = tmp_path / "msa"
    msa_dir.mkdir()
    (msa_dir / "gene1.fa").write_text(">t1\nACGT\n>t2\nACGT\n")
    (msa_dir / "gene2.fa").write_text(">t1\nTTTT\n>t2\nGGGG\n")
    return msa_dir


def test_symtest_help():
    runner = CliRunner()
    result = runner.invoke(main, ["pretree", "filter", "symtest", "--help"])
    assert result.exit_code == 0
    assert "--msa-dir" in result.output
    assert "--symtest-pval" in result.output


def test_symtest_requires_msa_dir():
    runner = CliRunner()
    result = runner.invoke(main, ["pretree", "filter", "symtest"])
    assert result.exit_code != 0


def test_symtest_dry_run(tmp_path):
    msa_dir = _make_msa_dir(tmp_path)
    output_dir = tmp_path / "out"
    runner = CliRunner()
    result = runner.invoke(main, [
        "pretree", "filter", "symtest",
        "--msa-dir", str(msa_dir),
        "--output-dir", str(output_dir),
        "--dry-run",
    ])
    assert result.exit_code == 0
    assert "Dry run" in result.output
    assert not (output_dir / "result.json").exists()


def test_symtest_invalid_pval(tmp_path):
    msa_dir = _make_msa_dir(tmp_path)
    runner = CliRunner()
    result = runner.invoke(main, [
        "pretree", "filter", "symtest",
        "--msa-dir", str(msa_dir),
        "--symtest-pval", "2.0",
    ])
    assert result.exit_code != 0


def test_symtest_invalid_pval_zero(tmp_path):
    msa_dir = _make_msa_dir(tmp_path)
    runner = CliRunner()
    result = runner.invoke(main, [
        "pretree", "filter", "symtest",
        "--msa-dir", str(msa_dir),
        "--symtest-pval", "0",
    ])
    assert result.exit_code != 0


def test_symtest_output_dir_conflict(tmp_path):
    msa_dir = _make_msa_dir(tmp_path)
    output_dir = tmp_path / "out"
    output_dir.mkdir()
    (output_dir / "existing.txt").write_text("data")
    runner = CliRunner()
    result = runner.invoke(main, [
        "pretree", "filter", "symtest",
        "--msa-dir", str(msa_dir),
        "--output-dir", str(output_dir),
    ])
    assert result.exit_code != 0
    assert "already exists" in result.output


def test_symtest_overwrite(tmp_path):
    msa_dir = _make_msa_dir(tmp_path)
    output_dir = tmp_path / "out"
    output_dir.mkdir()
    (output_dir / "old.txt").write_text("old")
    # This will fail at IQ-TREE invocation (tool not found), but not at
    # output-dir conflict -- it will overwrite and proceed.
    runner = CliRunner()
    result = runner.invoke(main, [
        "pretree", "filter", "symtest",
        "--msa-dir", str(msa_dir),
        "--output-dir", str(output_dir),
        "--overwrite",
    ])
    # Expects failure because no real iqtree, but not "already exists"
    assert "already exists" not in result.output


def test_symtest_missing_msa_dir(tmp_path):
    runner = CliRunner()
    result = runner.invoke(main, [
        "pretree", "filter", "symtest",
        "--msa-dir", str(tmp_path / "nonexistent"),
    ])
    assert result.exit_code != 0


def test_symtest_threads_negative(tmp_path):
    msa_dir = _make_msa_dir(tmp_path)
    runner = CliRunner()
    result = runner.invoke(main, [
        "pretree", "filter", "symtest",
        "--msa-dir", str(msa_dir),
        "--threads", "-1",
    ])
    assert result.exit_code != 0
```

- [ ] **Step 2: Run CLI tests**

```bash
cd /Users/zf/data/coding/phyloAI && python -m pytest tests/cli/test_pretree_filter_symtest.py -v
```

Expected: All PASS (dry-run and input validation tests don't need IQ-TREE; the overwrite test may exit with tool-not-found which is fine as it proves no conflict error).

- [ ] **Step 3: Commit**

```bash
git add tests/cli/test_pretree_filter_symtest.py
git commit -m "test(cli): add symtest CLI validation and dry-run tests"
```

---

### Task 4: Symtest integration test with mock iqtree

**Files:**
- Modify: `tests/pretree/test_symtest.py` (add integration test)

- [ ] **Step 1: Add mock-iqtree integration test**

Append to `tests/pretree/test_symtest.py`:

```python
# --- run_symtest integration with mock IQ-TREE ---

def test_run_symtest_with_mock_iqtree(tmp_path, monkeypatch):
    """Full integration: build supermatrix, invoke mock iqtree, parse output, filter."""
    from phyloai.pretree.filter import run_symtest
    import shutil

    # Create MSA files
    msa_dir = tmp_path / "msa"
    msa_dir.mkdir()
    (msa_dir / "gene1.fa").write_text(">t1\nACGT\n>t2\nACGT\n")
    (msa_dir / "gene2.fa").write_text(">t1\nTTTT\n>t2\nGGGG\n")
    (msa_dir / "gene3.fa").write_text(">t1\nCCCC\n>t2\nAAAA\n")

    # Create mock iqtree that writes a .symtest.csv on --symtest-only
    mock_iqtree = tmp_path / "mock_iqtree"
    mock_script = (
        "#!/usr/bin/env bash\n"
        "# Find the partitions file argument (-p <path>)\n"
        "partfile=\"\"\n"
        "while [ $# -gt 0 ]; do\n"
        '  if [ "$1" = "-p" ]; then shift; partfile="$1"; fi\n'
        "  shift\n"
        "done\n"
        'symcsv="${partfile}.symtest.csv"\n'
        "cat > \"$symcsv\" <<'CSVEOF'\n"
        "# comment\n"
        "Name,SymSig,SymNon,SymPval,MarSig,MarNon,MarPval,IntSig,IntNon,IntPval\n"
        "gene1,44,92,0.475,50,86,0.722,4,132,0.239\n"
        "gene2,43,93,0.004,49,87,0.003,5,131,0.170\n"
        "gene3,53,83,0.620,58,78,0.550,6,130,0.343\n"
        "CSVEOF\n"
        "exit 0\n"
    )
    mock_iqtree.write_text(mock_script)
    mock_iqtree.chmod(0o755)
    monkeypatch.setattr("shutil.which", lambda x, path=None: str(mock_iqtree))

    output_dir = tmp_path / "out"
    payload = run_symtest(
        msa_dir=msa_dir, output_dir=output_dir,
        symtest_type=None, symtest_pval=0.05,
        iqtree_path=mock_iqtree, threads=1,
        table_format="csv",
    )

    assert payload["status"] == "success"
    assert payload["key_results"]["n_input"] == 3
    assert payload["key_results"]["n_retained"] == 2  # gene1 + gene3
    assert payload["key_results"]["n_dropped"] == 1   # gene2

    # Check output files
    assert (output_dir / "result.json").exists()
    assert (output_dir / "filter.log").exists()
    assert (output_dir / "retained_loci.csv").exists()
    assert (output_dir / "dropped_loci.csv").exists()
    assert (output_dir / "filter_decisions.csv").exists()

    # Check seqs/
    seqs = output_dir / "seqs"
    assert (seqs / "gene1.fa").exists()
    assert (seqs / "gene3.fa").exists()
    assert not (seqs / "gene2.fa").exists()

    # Check filter_decisions content
    import csv
    with open(output_dir / "filter_decisions.csv") as fh:
        rows = list(csv.DictReader(fh))
    assert len(rows) == 3
    statuses = {r["locus"]: r["status"] for r in rows}
    assert statuses["gene1"] == "retained"
    assert statuses["gene2"] == "dropped"
    assert statuses["gene3"] == "retained"


def test_run_symtest_with_tree_dir(tmp_path, monkeypatch):
    """Tree copy: tree dir provided, only retained-locus trees copied."""
    from phyloai.pretree.filter import run_symtest

    msa_dir = tmp_path / "msa"
    msa_dir.mkdir()
    (msa_dir / "gene1.fa").write_text(">t1\nACGT\n>t2\nACGT\n")
    (msa_dir / "gene2.fa").write_text(">t1\nTTTT\n>t2\nGGGG\n")

    tree_dir = tmp_path / "trees"
    tree_dir.mkdir()
    (tree_dir / "gene1.tre").write_text("(t1,t2);")
    (tree_dir / "gene2.tre").write_text("(t1,t2);")

    mock_iqtree = tmp_path / "mock_iqtree"
    mock_script = (
        "#!/usr/bin/env bash\n"
        "partfile=\"\"\n"
        "while [ $# -gt 0 ]; do\n"
        '  if [ "$1" = "-p" ]; then shift; partfile="$1"; fi\n'
        "  shift\n"
        "done\n"
        'symcsv="${partfile}.symtest.csv"\n'
        "cat > \"$symcsv\" <<'CSVEOF'\n"
        "Name,SymSig,SymNon,SymPval,MarSig,MarNon,MarPval,IntSig,IntNon,IntPval\n"
        "gene1,44,92,0.475,50,86,0.722,4,132,0.239\n"
        "gene2,43,93,0.004,49,87,0.003,5,131,0.170\n"
        "CSVEOF\n"
        "exit 0\n"
    )
    mock_iqtree.write_text(mock_script)
    mock_iqtree.chmod(0o755)
    monkeypatch.setattr("shutil.which", lambda x, path=None: str(mock_iqtree))

    output_dir = tmp_path / "out"
    payload = run_symtest(
        msa_dir=msa_dir, output_dir=output_dir,
        symtest_type=None, symtest_pval=0.05,
        iqtree_path=mock_iqtree, threads=1,
        tree_dir=tree_dir, table_format="csv",
    )

    assert payload["status"] == "success"
    assert payload["key_results"]["retained_trees_copied"] == 1
    trees = output_dir / "trees"
    assert (trees / "gene1.tre").exists()
    assert not (trees / "gene2.tre").exists()


def test_run_symtest_iqtree_nonzero_exit(tmp_path, monkeypatch):
    """IQ-TREE non-zero exit code raises RuntimeError."""
    from phyloai.pretree.filter import run_symtest

    msa_dir = tmp_path / "msa"
    msa_dir.mkdir()
    (msa_dir / "gene1.fa").write_text(">t1\nACGT\n>t2\nACGT\n")

    mock_iqtree = tmp_path / "mock_iqtree"
    mock_iqtree.write_text("#!/usr/bin/env bash\necho 'SIMULATED ERROR' >&2\nexit 1\n")
    mock_iqtree.chmod(0o755)
    monkeypatch.setattr("shutil.which", lambda x, path=None: str(mock_iqtree))

    output_dir = tmp_path / "out"
    with pytest.raises(RuntimeError, match="exited with code 1"):
        run_symtest(
            msa_dir=msa_dir, output_dir=output_dir,
            iqtree_path=mock_iqtree, threads=1, table_format="csv",
        )
```

- [ ] **Step 2: Run integration tests**

```bash
cd /Users/zf/data/coding/phyloAI && python -m pytest tests/pretree/test_symtest.py -v
```

Expected: All tests PASS (mock iqtree writes predictable CSV).

- [ ] **Step 3: Commit**

```bash
git add tests/pretree/test_symtest.py
git commit -m "test(symtest): add integration tests with mock iqtree for supermatrix build and filtering"
```

---

### Task 5: Documentation updates

**Files:**
- Modify: `docs/commands/pretree-filter.md`
- Modify: `docs/superpowers/specs/2026-06-15-phyloai-pretree-filter-design.md`
- Modify: `docs/superpowers/specs/2026-06-07-phyloai-design.md`
- Modify: `README.md`

- [ ] **Step 1: Update pretree-filter.md -- add symtest section between metrics and cluster**

Insert the following content after line 433 (end of `filter metrics` notes section, before `## filter cluster`):

```markdown
---

## `filter symtest` -- Symmetry Test Filtering

### Purpose

Run IQ-TREE3's tests of symmetry (Naser-Khdour et al., 2019) via `--symtest-only` to detect loci that violate phylogenetic assumptions of stationarity, homogeneity, or reversibility, then filter out loci with p-value below a configurable threshold.

This is locus-level filtering by statistical test: loci that fail symmetry are discarded. This does not compute metrics (use `pretree metrics`) or mask/prune individual sites/taxa (use `filter taper` or `filter treeshrink`).

### Usage

```bash
phyloai pretree filter symtest \
  --msa-dir <msa_dir> \
  [--symtest-type MAR|INT] \
  [--symtest-pval 0.05] \
  [--symtest-keep-zero] \
  [--iqtree-path <path>] \
  [--threads 4] \
  [--tree-dir <tree_dir>] \
  [--output-dir runs/pretree/filter/symtest] \
  [--table-format csv|tsv] \
  [--dry-run] [--quiet] [--overwrite]
```

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--msa-dir` | required | Directory of per-locus MSA files. Every regular non-empty file is scanned. |
| `--symtest-type` | — | `MAR` (marginal/stationarity test), `INT` (internal/homogeneity test). When omitted, the combined Sym test is used (SymPval column). |
| `--symtest-pval` | 0.05 | P-value threshold. Loci with p >= threshold are retained; those below are dropped. Must be > 0 and <= 1. |
| `--symtest-keep-zero` | off | Pass `--symtest-keep-zero` to IQ-TREE (retain NAs in the tests). |
| `--iqtree-path` | — | Explicit path to iqtree binary. Resolved via PATH. |
| `--threads` / `-t` | 4 | IQ-TREE threads (`-T`). |
| `--tree-dir` | — | Optional gene tree directory. Trees matching retained loci are copied to `trees/`. |
| `--output-dir` / `-o` | `runs/pretree/filter/symtest` | Output directory. |
| `--table-format` | `csv` | Format for auxiliary tables. |
| `--overwrite` | off | Delete and recreate output directory. |
| `--dry-run` | off | Show resolved IQ-TREE command and locus count. |
| `--quiet` / `-q` | off | Suppress terminal output. |

### Inputs

`--msa-dir` is the only required input. The command:

1. Builds a temporary supermatrix + RAxML-style partition file from all MSAs
2. Runs `iqtree -s <matrix> -p <partitions> --symtest-only`
3. Parses `<partitions>.symtest.csv` for per-partition p-values
4. Applies the p-value threshold on the selected test column

The p-value column used depends on `--symtest-type`:
- (default) -> `SymPval` (combined stationarity + homogeneity)
- `MAR` -> `MarPval` (marginal / stationarity test)
- `INT` -> `IntPval` (internal / homogeneity test)

Temporary files are cleaned up after the run. Unlike TAPER, `--resume` is not supported (single IQ-TREE invocation, deterministic parsing).

### Outputs

```
runs/pretree/filter/symtest/
├── seqs/                              (retained MSAs)
├── trees/                             (only when --tree-dir provided)
├── retained_loci.csv|tsv
├── dropped_loci.csv|tsv               (locus, reason)
├── filter_decisions.csv|tsv           (locus, status, p_value, symtest_type,
│                                       sym_sig, sym_non, mar_sig, mar_non,
│                                       int_sig, int_non)
├── filter.log
└── result.json
```

Terminal output: Filter Results table (input/retained/dropped/p-value threshold/symtest type) + Retained MSA Statistics table + optionally Trees Copied table.

`result.json.key_results`: `n_input`, `n_retained`, `n_dropped`, `p_value_threshold`, `symtest_type`, `retained_trees_copied`.

### Examples

```bash
# Default symmetry test (Sym), p < 0.05 dropped
phyloai pretree filter symtest --msa-dir ./trimmed

# Marginal symmetry (stationarity) test
phyloai pretree filter symtest --msa-dir ./trimmed --symtest-type MAR

# Stricter threshold
phyloai pretree filter symtest --msa-dir ./trimmed --symtest-pval 0.01

# With tree directory: retain matching gene trees
phyloai pretree filter symtest \
  --msa-dir ./trimmed --tree-dir ./genetrees

# Internal homogeneity test
phyloai pretree filter symtest --msa-dir ./trimmed --symtest-type INT

# Dry-run to inspect command
phyloai pretree filter symtest --msa-dir ./trimmed --dry-run
```

### Warnings and Errors

| Condition | Behaviour |
|-----------|-----------|
| `--msa-dir` is empty or no valid MSA files | Exit 1 |
| IQ-TREE not found | Exit 3 |
| IQ-TREE exits non-zero | Exit 2 |
| `.symtest.csv` missing or unparseable | Exit 2 |
| Partition name / locus name mismatch | Exit 1 |
| All loci dropped (all p < threshold) | Success with warning |
| `--tree-dir` provided but no matching trees | Success with warning |
| Non-empty output directory without `--overwrite` | Exit 1 |
| `--symtest-pval` <= 0 or > 1 | Exit 1 |
| `--threads` < 1 | Exit 1 |

### Notes

Symmetry testing should be run after alignment and trimming but before supermatrix concatenation, since violations of stationarity or homogeneity can bias phylogenetic inference. The `--symtest-type` default (combined Sym test) is the most general and widely applicable.

References: Naser-Khdour et al. (2019) "Assessing the Goodness of Fit of Phylogenetic Models..." doi:10.1093/gbe/evz193.
```

- [ ] **Step 2: Update filter design spec header (lines 1-18)**

In `docs/superpowers/specs/2026-06-15-phyloai-pretree-filter-design.md`, update Section 1 (Purpose) line 11:
```
`phyloai pretree filter` runs after `pretree metrics` and before `pretree concat`. It provides four filtering workflows:
```
Change "four" to "five":
```
`phyloai pretree filter` runs after `pretree metrics` and before `pretree concat`. It provides five filtering workflows:
```

After line 15 (end of item 4), add:
```
5. **Symtest filtering**: run IQ-TREE3's `--symtest-only` to test phylogenetic symmetry assumptions (stationarity, homogeneity, reversibility) per locus, then filter by p-value threshold.
```

- [ ] **Step 3: Update filter design spec command structure (line 62)**

After `phyloai pretree filter metrics` line (line 64), add:
```bash
phyloai pretree filter symtest
```
so the list reads: taper, treeshrink, metrics, symtest, cluster.

- [ ] **Step 4: Update filter design spec output dirs (line 86)**

After the `metrics` line, add:
```
- `filter symtest`: `runs/pretree/filter/symtest`
```

- [ ] **Step 5: Update main design spec (Section 4.1)**

In `docs/superpowers/specs/2026-06-07-phyloai-design.md` line 111 (after `phyloai pretree filter metrics`), add:
```bash
phyloai pretree filter symtest  --msa-dir ./trimmed [--tree-dir ./genetrees]
```

- [ ] **Step 6: Update README.md filter row (line 90)**

Replace:
```
| `phyloai pretree filter`  | Marker-level filtering: TAPER error-site masking, TreeShrink taxon pruning, metric-rule filtering, cluster-based exploration. | [docs/commands/pretree-filter.md](docs/commands/pretree-filter.md) |
```
with:
```
| `phyloai pretree filter`  | Marker-level filtering: TAPER error-site masking, TreeShrink taxon pruning, metric-rule filtering, symmetry test filtering, cluster-based exploration. | [docs/commands/pretree-filter.md](docs/commands/pretree-filter.md) |
```

- [ ] **Step 7: Commit**

```bash
git add docs/commands/pretree-filter.md docs/superpowers/specs/2026-06-15-phyloai-pretree-filter-design.md docs/superpowers/specs/2026-06-07-phyloai-design.md README.md
git commit -m "docs: add symtest to filter documentation, design specs, and README"
```

---

### Task 6: Final verification

- [ ] **Step 1: Run all filter tests**

```bash
cd /Users/zf/data/coding/phyloAI && python -m pytest tests/pretree/test_symtest.py tests/pretree/test_filter.py tests/cli/test_pretree_filter_symtest.py -v
```

Expected: All PASS.

- [ ] **Step 2: Verify lint**

```bash
cd /Users/zf/data/coding/phyloAI && ruff check phyloai/pretree/filter.py phyloai/cli/commands/pretree.py tests/pretree/test_symtest.py tests/cli/test_pretree_filter_symtest.py
```

Expected: No errors.

- [ ] **Step 3: Run full test suite**

```bash
cd /Users/zf/data/coding/phyloAI && python -m pytest -x --timeout=60 2>&1 | tail -20
```

Expected: All pass (some may skip due to missing external tools).

- [ ] **Step 4: Commit any fixes**

```bash
git add -A && git commit -m "chore(symtest): fix lint and test issues from final verification"
```
