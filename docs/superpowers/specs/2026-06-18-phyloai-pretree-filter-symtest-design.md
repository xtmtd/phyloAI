# phyloai pretree filter symtest Design Specification

**Date:** 2026-06-18  
**Status:** Draft for review  
**Reference:** `docs/superpowers/specs/2026-06-07-phyloai-design.md`, `docs/superpowers/specs/2026-06-15-phyloai-pretree-filter-design.md`, `docs/superpowers/specs/2026-06-13-phyloai-pretree-concat-design.md`

---

## 1. Purpose

`phyloai pretree filter symtest` is a 5th subcommand under `pretree filter`, placed between `metrics` and `cluster`. It implements IQ-TREE3's tests of symmetry (Naser-Khdour et al., 2019) via `--symtest-only` to detect loci that violate stationarity or homogeneity assumptions. Loci with p-value below a configurable threshold are filtered out.

The command:
1. Scans `--msa-dir` for input MSAs
2. Generates a temporary supermatrix + IQ-TREE partition file (reusing `concat.py` internal helpers)
3. Runs a single `iqtree -s -p --symtest-only` invocation
4. Parses `<partition>.symtest.csv` for per-partition p-values
5. Retains loci with `p >= --symtest-pval`, drops loci with `p < --symtest-pval`
6. Copies retained MSAs to `seqs/` and optionally retained gene trees to `trees/`

What it does **not** do:
- compute marker metrics; use `pretree metrics`
- infer gene trees; use later tree commands
- concatenate retained MSAs; use `pretree concat`
- remove individual taxa or mask sites; use `taper` or `treeshrink`

---

## 2. CLI

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

### 2.1 Parameter Details

| Parameter | Type | Default | Notes |
|-----------|------|---------|-------|
| `--msa-dir` | Path | **required** | Directory containing per-locus MSA files |
| `--symtest-type` | `MAR\|INT` | none (use default `Sym`) | IQ-TREE `--symtest-type`; when omitted IQ-TREE runs all three tests and the p-value column used for filtering is `SymPval` |
| `--symtest-pval` | float | `0.05` | P-value threshold for PhyloAI-side filtering only (NOT passed to IQ-TREE) |
| `--symtest-keep-zero` | flag | off | IQ-TREE `--symtest-keep-zero`; keep NAs in the tests |
| `--iqtree-path` | Path | auto-detect | Custom path to `iqtree` binary; falls back to `core/env.py` detection. IQ-TREE3 >= 2.3.0 required (introduced `--symtest-only`). |
| `--threads`, `-t` | int | `4` | IQ-TREE `-T` |
| `--tree-dir` | Path | none | If provided, copy retained gene trees to `trees/` under output dir; matched to retained MSAs by locus name |
| `--output-dir`, `-o` | Path | `runs/pretree/filter/symtest` | Output directory |
| `--table-format` | `csv\|tsv` | `csv` | Format for auxiliary tables; file suffix follows format |
| `--dry-run` | flag | off | Validate and show planned work without writing files |
| `--quiet`, `-q` | flag | off | Suppress terminal output except errors |
| `--overwrite` | flag | off | Replace a non-empty output directory |

`--seq-type` is not needed: IQ-TREE auto-detects sequence type from the input. `--resume` is not applicable: the command runs one IQ-TREE invocation (deterministic parsing), not a batch of per-locus tasks.

### 2.2 IQ-TREE Command

P-values are parsed from the selected column in the symmetry test output, not via `--symtest-remove-bad` inside IQ-TREE. This keeps filtering logic in PhyloAI for transparency and consistent output tables.

Constructed command:

```bash
iqtree -s <temp_matrix.fa> -p <temp_partitions.txt> --symtest-only \
  [--symtest-type MAR|INT] [--symtest-keep-zero] \
  -T <threads>
```

`--symtest-pval` is NOT passed to IQ-TREE. It is used only for PhyloAI-side filtering after parsing the CSV output, keeping filtering logic transparent and visible in the decision tables.

| `--symtest-type`  | Column used | Test meaning                      |
| ----------------- | ----------- | --------------------------------- |
| not specified     | `SymPval`   | Combined stationarity + homogeneity |
| `MAR`             | `MarPval`   | Marginal (stationarity)            |
| `INT`             | `IntPval`   | Internal (homogeneity)             |

The output file is `<partition_file>.symtest.csv` with columns: `Name`, `SymSig`, `SymNon`, `SymPval`, `MarSig`, `MarNon`, `MarPval`, `IntSig`, `IntNon`, `IntPval`.

Managed IQ-TREE flags (blocked in `--tool-args`): `-s`, `-p`, `--symtest-only`, `--symtest-type`, `--symtest-keep-zero`, `-T`, `--redo`.

The p-value column used for filtering depends on `--symtest-type`:

---

## 3. Architecture

```
phyloai/core/file_matching.py          # scan_msa_dir
phyloai/core/env.py                    # iqtree detection (TOOL_REGISTRY key: iqtree3)
phyloai/core/runner.py                 # Runner.run() for single IQ-TREE invocation
phyloai/pretree/filter.py              # run_symtest + shared helpers
phyloai/pretree/concat.py              # internal reuse: _read_msa, _concat_alignments,
                                       #   _write_partitions
phyloai/cli/commands/pretree.py        # filter_symtest_command + registration
docs/commands/pretree-filter.md        # update with symtest section
```

### 3.1 Concat Integration

Two internal functions from `concat.py` are reused:

| Function | Purpose |
|----------|---------|
| `_read_msa(filepath)` | Format-agnostic MSA reading (FASTA/Nexus/Phylip via FormatConverter); returns `(taxa, seqs, length)` |
| `_write_partitions(out_path, genes, prefix_type)` | Write RAxML-style partition file |

Supermatrix assembly and partition position tracking are built inline in `run_symtest` (approximately 30 lines) because `_concat_alignments` does not return partition position information. The inline logic follows the same pattern as `run_concat` (concat.py lines 554-584).

### 3.2 Execution Flow

```
1. scan_msa_dir(--msa-dir) → {locus_name: Path}
2. Read all MSAs via _read_msa from concat.py (format-agnostic: FASTA/Nexus/Phylip), discover full taxon set
3. Build supermatrix + partition positions inline (reusing _write_partitions from concat.py)
4. Write supermatrix to temp file (FASTA, 60-char line wrap per §9.11)
5. _write_partitions → temp partition file (RAxML format)
6. Runner.run(iqtree cmd) → ToolResult
7. Parse <partition_file>.symtest.csv
8. Cross-validate: all CSV Name entries must exist in msa_map; unmatched entries are warned
9. Select p-value column per --symtest-type
10. Filter: retain loci with p >= --symtest-pval
11. Copy retained MSAs → seqs/
12. If --tree-dir: match retained loci to trees, copy → trees/
13. Clean up temp dir
14. Write result.json (single pattern: data.cmd, data.tool_stderr, data.results), decision tables
```

Temporary files are written to a `tempfile.mkdtemp()` directory and cleaned up after step 12. IQ-TREE working files (`.log`, `.iqtree`, etc.) are also in the temp dir and cleaned up.

### 3.3 Tree Copy Behavior

When `--tree-dir` is provided:
- Only trees matching retained loci (by logical locus name) are copied
- Locus matching follows the global file-matching policy (§9.7): suffix-agnostic, 1-2 dot reductions
- Ambiguous tree matches for a given locus exit with code 1
- Tree files without a matching retained locus are silently skipped
- This mirrors `metrics --copy` behavior

---

## 4. Outputs

### 4.1 Output Layout

```
runs/pretree/filter/symtest/
├── result.json                     # data.cmd, data.tool_stderr (single pattern)
├── seqs/                           # retained MSAs
├── trees/                          # retained trees (only when --tree-dir)
├── retained_loci.csv|tsv
├── dropped_loci.csv|tsv
└── filter_decisions.csv|tsv
```

### 4.2 filter_decisions.csv

| Column | Description |
|--------|-------------|
| `locus` | Logical locus name |
| `status` | `retained` or `dropped` |
| `p_value` | P-value from the selected test column |
| `symtest_type` | Which column was used: `Sym`, `MAR`, or `INT` |
| `sym_sig` | Number of significant pairs (symmetry test) |
| `sym_non` | Number of non-significant pairs (symmetry test) |
| `mar_sig` | Number of significant pairs (marginal test) |
| `mar_non` | Number of non-significant pairs (marginal test) |
| `int_sig` | Number of significant pairs (internal test) |
| `int_non` | Number of non-significant pairs (internal test) |

### 4.3 Terminal Summary

Rich table output includes:
- Total input loci
- Retained count, dropped count
- Symtest type used
- P-value threshold
- Retained MSA statistics: retained MSA count, total concat length, mean/min/max marker length, mean taxa count
- When `--tree-dir`: retained tree count, missing tree match count

### 4.4 result.json

Follows single-pattern schema (JSON Output Standard §5.2). Key results:
```json
{
  "status": "success",
  "command": "phyloai pretree filter symtest ...",
  "wall_time": 45.2,
  "tool_versions": {"iqtree3": "2.3.x"},
  "params": {"msa_dir": "...", "symtest_pval": 0.05, "symtest_type": "Sym", "threads": 4},
  "key_results": {
    "n_input": 100,
    "n_retained": 85,
    "n_dropped": 15,
    "p_value_threshold": 0.05,
    "symtest_type": "Sym",
    "retained_trees_copied": 0
  },
  "error": null,
  "data": {
    "cmd": ["iqtree3", "-s", "symtest_matrix.fa", "-p", "symtest_partitions.txt", "--symtest-only", "-T", "4"],
    "tool_stderr": "IQ-TREE symtest output ...",
    "summary": {
      "n_input": 100,
      "n_retained": 85,
      "n_dropped": 15,
      "p_value_threshold": 0.05,
      "symtest_type": "Sym",
      "retained_msa_stats": {"n_msa": 85, "total_length": 42500, "mean_length": 500},
      "retained_tree_count": 0,
      "missed_tree_count": 0,
      "skipped_names": []
    },
    "results": [
      {"locus": "gene1", "status": "retained"},
      {"locus": "gene2", "status": "dropped", "reason": "SymPval 0.12 > threshold 0.05"}
    ]
  }
}
```

---

## 5. Error Handling

| Case | Exit code | Behavior |
|------|-----------|----------|
| `--msa-dir` missing or empty | 1 | Error: no valid MSA input |
| IQ-TREE not found | 3 | Environment error with install guidance |
| IQ-TREE < 2.3.0 (no --symtest-only) | 2 | Tool error from IQ-TREE stderr |
| IQ-TREE non-zero exit | 2 | External tool failure; stderr in log |
| `.symtest.csv` not generated | 2 | Tool error; stderr in log |
| `.symtest.csv` unparseable | 2 | Parse error with file excerpt |
| Partition name / locus name mismatch | 1 | Ambiguity error with mismatched names |
| All loci dropped (all p < threshold) | 0 | Success with empty `seqs/`; warning in terminal |
| `--tree-dir` provided but no tree matches | 0 | Success; warning about 0 trees copied |
| Ambiguous tree match for a locus | 1 | Error with conflicting tree candidates |
| `--tool-args` blocks managed flag | 1 | Error with blocked flag name |

---

## 6. Shared Component Changes

### 6.1 filter.py

Add `run_symtest` function following the pattern of existing `run_taper`, `run_treeshrink`, etc.:
- Signature: `run_symtest(msa_dir, symtest_type, symtest_pval, symtest_keep_zero, iqtree_path, threads, tree_dir, output_dir, table_format, dry_run, quiet, tool_args=None) -> None`
- Writes `result.json`, decision tables
- Renders Rich terminal summary

### 6.2 concat.py

If `_read_msa`, `_concat_alignments`, or `_write_partitions` have implicit dependencies on CLI context or global state, refactor them to accept explicit parameters. The underlying logic should not change; only the function signatures may need adjustment for reuse.

### 6.3 pretree.py CLI

Register `symtest` as the 5th filter subcommand overall (placed after `metrics`, before `cluster` in the `list_commands` order). Follow the existing `_FilterGroup` pattern with Click options.

### 6.4 Documentation Updates

- `docs/commands/pretree-filter.md`: add symtest section with examples
- `docs/superpowers/specs/2026-06-15-phyloai-pretree-filter-design.md` (main filter spec): update Section 1 (Purpose) and Section 3 (Command Structure) to reflect 5 subcommands; update Section 2 architecture to mention iqtree as external tool; add symtest to Section 3.1 per-subcommand default output dirs
- `docs/superpowers/specs/2026-06-07-phyloai-design.md`: update Section 4.1 command examples to include `phyloai pretree filter symtest`
- `README.md`: update filter subcommand list

---

## 7. Testing Strategy

### 7.1 Unit Tests

- `.symtest.csv` parsing (valid and malformed)
- P-value column selection per `--symtest-type`
- Locus name matching between partition map and CSV
- Filtering logic (p >= threshold retained, p < threshold dropped)
- All-loci-dropped edge case
- All-loci-retained edge case

### 7.2 CLI Tests

- `--msa-dir` required validation
- `--symtest-pval` must be > 0
- `--dry-run` writes no files
- `--tree-dir` copies only retained trees
- Output directory conflict (no `--overwrite` → exit 1)
- `--overwrite` replaces existing output dir

### 7.3 External Tool Tests

- Mock `iqtree` executable that:
  - Accepts `-s <matrix> -p <partitions> --symtest-only`
  - Writes a predictable `<partitions>.symtest.csv`
  - Returns exit code 0
- Mock `iqtree` that returns non-zero exit code
- Mock `iqtree` that writes no `.symtest.csv`
- Do not require IQ-TREE in CI

---

## 8. Design Decisions

| Decision | Rationale |
|----------|-----------|
| Single IQ-TREE invocation, not per-locus | `--symtest-only` with `-s -p` processes all partitions in one run; per-locus invocations would duplicate model initialization |
| No `--seq-type` | IQ-TREE auto-detects sequence type; parameter adds no value |
| No `--resume` | Single deterministic invocation; no per-locus checkpoint needed |
| No `--copy` flag; MSAs always copied | Symtest is a filtering step that should always produce a usable output directory; unlike `metrics` exploration where you might just want the table |
| `--tree-dir` copies only retained trees | Matches `metrics --copy` semantics; dropped loci' trees are excluded |
| Parse `.symtest.csv` in PhyloAI, not use `--symtest-remove-bad` | PhyloAI controls the filtering logic and writes consistent decision tables. `--symtest-pval` NOT passed to IQ-TREE; used Python-side only. |
| Temp dir with automatic cleanup | Keeps output directory clean; IQ-TREE working files are large and not useful for post-hoc inspection |
| `_read_msa` / `_write_partitions` from concat.py | Reuse existing, format-agnostic MSA reading (FASTA/Nexus/Phylip) and partition writing logic |
| Supermatrix assembly inline in `run_symtest` | `_concat_alignments` does not return partition positions; inline logic (~30 lines) follows the same pattern as concat.py:554-584 |
| Cross-validate CSV names against msa_map | IQ-TREE may sanitize partition names containing special characters; cross-validation catches mismatches before silent data loss |
| Symtest placed between metrics and cluster in subcommand order | Logical progression: rule-based filtering (metrics) → statistical test filtering (symtest) → exploratory cluster filtering (cluster) |
