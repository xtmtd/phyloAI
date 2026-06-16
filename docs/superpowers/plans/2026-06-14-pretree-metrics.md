# pretree metrics Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Implement `phyloai pretree metrics`, `metrics plot`, and `metrics correlate` — molecular marker MSA+tree attribute computation, distribution visualization, and correlation analysis.

**Architecture:** Core library in `phyloai/pretree/metrics.py` (compute, plot, correlate, helpers); CLI registration appended to `phyloai/cli/commands/pretree.py`; FastTree added to `core/env.py`. Follows `pretree trim` and `pretree stats` patterns throughout.

**Tech Stack:** Python 3.10+, BioPython, Click 8+, Rich, ProcessPoolExecutor, numpy, scipy, matplotlib

---

## File Map

| Action | File | Responsibility |
|--------|------|----------------|
| Create | `phyloai/pretree/metrics.py` | All metrics logic: MSA metrics, tree metrics, saturation, plot generation, correlation analysis |
| Modify | `phyloai/cli/commands/pretree.py` | Register `metrics` Click group with `plot` and `correlate` subcommands |
| Modify | `phyloai/core/env.py` | Add `fasttree` to `TOOL_REGISTRY` |
| Create | `tests/pretree/test_metrics.py` | Unit + integration tests for metrics module |

---

## Task 1: Core Data Structures and Helpers

**Files:** Create `phyloai/pretree/metrics.py` (helpers section), Create `tests/pretree/test_metrics.py` (test stubs)

Implement shared utilities for MSA loading, file pairing, and taxon name checks.

- [ ] **Step 1.1: FASTA header scanning (two-pass taxa pool)**

  ```python
  def _scan_msa_headers(msa_dir: Path) -> tuple[dict[str, set[str]], set[str]]:
      """Pass 1: scan all MSA files, collect taxon names per marker + total pool."""
      per_marker: dict[str, set[str]] = {}
      total_pool: set[str] = set()
      for path in ...
          records = list(SeqIO.parse(path, "fasta"))
          taxa = {record.id for record in records}
          per_marker[path.stem] = taxa
          total_pool.update(taxa)
      return per_marker, total_pool
  ```

- [ ] **Step 1.2: File pairing logic**

  ```python
  def _pair_files(
      msa_dir: Path | None,
      tree_dir: Path | None,
  ) -> tuple[dict[str, tuple[Path | None, Path | None]], list[str]]:
      """Match MSA and tree files by logical locus name via the shared helper."""
  ```

  - Scan `--msa-dir` for `.fa`, `.fasta`, `.fas`, `.fna`, `.faa`, `.aln`
  - Scan `--tree-dir` for `.tre`, `.tree`, `.nwk`, `.newick`, `.treefile`, `.bestTree`, `.contree`
  - MSA logical loci derive from the filename text before the final dot (for example, `gene.v1.ALI` → `gene.v1`)
  - Tree files try two logical-locus candidates from the shared helper: remove one suffix, and when possible also remove two suffixes (for example, `gene.fa.treefile` → `gene.fa` and `gene`)
  - If both tree candidates match different MSA logical loci, raise an ambiguity error instead of guessing
  - Match by logical locus name (case-sensitive) using the shared helper
  - Unmatched MSA or tree files are collected as warnings
  - At least one of `--msa-dir` or `--tree-dir` required

- [ ] **Step 1.3: MSA–tree taxon consistency check**

  ```python
  def _check_taxon_consistency(msa_path: Path, tree_path: Path) -> dict | None:
      """Verify set(msa_taxa) == set(tree_taxa). Return mismatch info or None."""
  ```

  - Read MSA headers via `SeqIO.parse`
  - Read tree tips via `Phylo.read` + `tree.get_terminals()`
  - Compare sets, return `{"marker": stem, "msa_only": [...], "tree_only": [...]}` or None

- [ ] **Step 1.4: Sequence type resolution**

  Reuse `core/sequence_normalization.resolve_seq_type` for `--seq-type auto`. Map `NT` → `DNA` (BioPython convention), pass through `AA`.

- [ ] **Step 1.5: Standardize sequence helper**

  ```python
  def _standardize_sequence(sequence: str, seq_type: str) -> str:
      """Replace non-standard chars with gap (-)."""
  ```

  - AA: valid = `ARNDCQEGHILKMFPSTWYV`, others → `-`
  - DNA: valid = `ACGT`, U→T, others → `-`
  - Mirror `extract_msa_tree_features.py:standardize_sequence`

### Verification for Task 1

```bash
# Run the test stubs - expect ImportError/skip until implementation complete
python -c "from phyloai.pretree.metrics import _scan_msa_headers, _pair_files, _check_taxon_consistency; print('imports OK')"
```

---

## Task 2: MSA Metric Computation

**Files:** Modify `phyloai/pretree/metrics.py` (MSA section), Modify `tests/pretree/test_metrics.py` (MSA tests)

Implement all 22 MSA metrics as defined in the design spec (Section 4.1).

- [ ] **Step 2.1: Core counts and frequencies**

  ```python
  def _compute_msa_metrics(
      msa_path: Path,
      seq_type: str,
      total_taxa_pool: int,
      skip_freq: bool,
  ) -> dict[str, Any]:
  ```

  **Metrics computed:**
  1. `num_taxa` — `len(msa)`
  2. `taxa_occupancy` — `num_taxa / total_taxa_pool`
  3. `num_sites` — `msa.get_alignment_length()`
  4. `num_patterns` — normalize non-standard→`-`, hash columns, `len(set(...))`
  5. `proportion_patterns` — `num_patterns / num_sites`
  6. `num_parsimony_sites` — chars with count≥2, at least 2 such chars
  7. `proportion_parsimony` — `num_parsimony_sites / num_sites`
  8. `num_singletons` — variable sites not PI (stats.py logic)
  9. `proportion_singletons` — `num_singletons / num_sites`
  10. `num_sites/num_taxa`
  11. `num_patterns/num_taxa`
  12. `num_parsimony_sites/num_taxa`
  13. `num_singletons/num_taxa`
  14. `proportion_gaps` — `total_gaps / (nsites × ntaxa)`
  15. `proportion_invariant` — sites with ≤1 unique valid char
  16. `entropy` — mean per-site Shannon (valid chars only)
  17. `bollback` — `pattern_entropy - nsites × log(nsites)`
  18. `pattern_entropy` — `Σ N_i × log(N_i)`
  19. `rcfv` — classic RCFV (frequency-based, per `extract_msa_tree_features.py`)
  20. `nrcfv` — `RCFV / (p^(-0.5) × n^(0.01) × c) × 100`
  21. `average_pairwise_identity` — mean pairwise identity
  22. `GC_content` — NT only, `(G+C)/total_valid`; `""` for AA

  **Key implementation detail — `rcfv`:**
  ```python
  rcfv = 0.0
  for state in valid_chars:
      state_freqs = []
      for seq in msa:
          seq_str = str(seq.seq)
          total = sum(1 for c in seq_str if c in valid_chars)
          if total > 0:
              state_freqs.append(seq_str.count(state) / total)
      if state_freqs:
          mean_freq = sum(state_freqs) / len(state_freqs)
          rcfv += sum(abs(f - mean_freq) for f in state_freqs) / ntaxa
  ```

  **Key implementation detail — `nrcfv`:**
  ```python
  c = 4 if seq_type == "DNA" else 20
  nrcfv = rcfv / (num_sites ** (-0.5) * ntaxa ** 0.01 * c) * 100
  ```

  **Key implementation detail — `num_patterns` (stats.py logic):**
  ```python
  gap_byte = ord('-')
  standard_codes = {ord(ch) for ch in standard_chars}
  pattern_set = set()
  for col in zip(*(s.encode() for s in sequences)):
      normalized = bytes(gap_byte if ch not in standard_codes else ch for ch in col)
      pattern_set.add(normalized)
  num_patterns = len(pattern_set)
  ```

  **Key implementation detail — `num_singletons`:**
  Follow `stats.py:compute_site_patterns`:
  ```python
  singleton_count = 0
  for column in ...:
      standard_counts = Counter(ch for ch in column if ord(ch) in standard_codes)
      if sum(standard_counts.values()) < 2:
          continue  # constant
      if len(standard_counts) == 1:
          continue  # constant
      repeated = sum(1 for v in standard_counts.values() if v >= 2)
      if repeated >= 2:
          parsimony_informative += 1
      else:
          singleton_count += 1
  ```

- [ ] **Step 2.2: Frequency statistics (optional)**

  When `--skip-freq-statistics` is False:
  - NT: `freqA`, `freqC`, `freqG`, `freqT`
  - AA: `freqA`, `freqR`, `freqN`, `freqD`, `freqC`, `freqQ`, `freqE`, `freqG`, `freqH`, `freqI`, `freqL`, `freqK`, `freqM`, `freqF`, `freqP`, `freqS`, `freqT`, `freqW`, `freqY`, `freqV`
  - All frequencies sum to 1.0 (gaps excluded)

- [ ] **Step 2.3: MSA tests**

  - Create a small test MSA FASTA file (temporary)
  - Verify `num_taxa`, `num_sites`
  - Verify `proportion_gaps` equals `stats.py:gap_ratio`
  - Verify `num_patterns` matches stats.py approach (non-standard→gap)
  - Verify `rcfv` matches manual calculation
  - Verify `nrcfv` formula correctness
  - Verify `GC_content` for NT, empty for AA
  - Verify singleton detection logic

### Verification for Task 2

```bash
python -m pytest tests/pretree/test_metrics.py -k "test_msa" -v
```

---

## Task 3: Tree Metric Computation

**Files:** Modify `phyloai/pretree/metrics.py` (tree section), Modify `tests/pretree/test_metrics.py` (tree tests)

Implement all 15 tree metrics as defined in the design spec (Section 4.2).

- [ ] **Step 3.1: Basic tree metrics (items 1–12)**

  ```python
  def _compute_tree_metrics(
      tree_path: Path,
      outgroup_list: Path | None,
      ref_tree_path: Path | None,
  ) -> dict[str, Any]:
  ```

  1. `average_BS` — mean of `clade.confidence` for internal nodes. Note: `Bio.Phylo` parses Newick comments following node labels as `.confidence` values (e.g., `(A,B)0.95:0.1` → `confidence=0.95`). Nodes without confidence values are excluded from the mean.
  2. `sd_BS` — stdev of bootstrap values
  3. `total_tree_length` — sum of all `clade.branch_length`
  4. `average_internal_branch_length` — mean of non-terminal brlen
  5. `sd_internal_branch_length` — stdev
  6. `average_terminal_branch_length` — mean of terminal brlen
  7. `sd_terminal_branch_length` — stdev
  8. `tree_diameter` — max patristic distance
  9. `average_patristic_distance` — mean across all pairs
  10. `sd_patristic_distance` — stdev
  11. `evo_rate` — `total_tree_length / n_taxa`
  12. `treeness` — `sum(internal_brlen) / total_tree_length`

- [ ] **Step 3.2: DVMC (item 13)**

  ```python
  def _compute_dvmc(tree: Tree, outgroup_list: Path | None) -> float:
  ```

  - If `outgroup_list` provided: read outgroup names, prune from tree
  - Compute: `distances = [tree.distance(term) for term in tree.get_terminals()]`
  - `dvmc = np.std(distances, ddof=1)`  # ≡ sqrt formula from both implementations
  - Handle edge case: `num_spp < 2` → `""` (NA)

- [ ] **Step 3.3: RF distance (item 15, with `--ref-tree`)**

  ```python
  def _compute_rf_distance(tree: Tree, ref_tree: Tree) -> float:
  ```

  - Prune both trees to shared tips
  - Root both with first shared tip
  - Compute plain RF by comparing clade sets
  - `normalized_rf = plain_rf / (2 * (tip_count - 3))`

- [ ] **Step 3.4: Saturation (item 14)**

  ```python
  def _compute_saturation(msa_path: Path, tree_path: Path, exclude_gaps: bool = True) -> float:
  ```

  Implementation following phykit `saturation.py` algorithm:

  ```python
  alignment = AlignIO.read(msa_path, "fasta")
  tree = Phylo.read(tree_path, "newick")
  tips = [term.name for term in tree.get_terminals()]

  # Build seq_arrays dict for fast lookup
  seq_arrays = {record.name: np.array(list(str(record.seq).upper())) for record in alignment}

  patristic = []
  uncorrected = []
  for i in range(len(tips)):
      for j in range(i + 1, len(tips)):
          pd = tree.distance(tips[i], tips[j])
          patristic.append(pd)

          seq1 = seq_arrays[tips[i]]
          seq2 = seq_arrays[tips[j]]
          if exclude_gaps:
              valid = (seq1 != '-') & (seq2 != '-')
              if valid.any():
                  identity = np.sum(seq1[valid] == seq2[valid])
                  ud = 1.0 - identity / valid.sum()
              else:
                  ud = float('nan')
          else:
              identity = np.sum(seq1 == seq2)
              ud = 1.0 - identity / len(seq1)
          uncorrected.append(ud)

  x = np.array(patristic, dtype=float)
  y = np.array(uncorrected, dtype=float)
  mask = np.isfinite(x) & np.isfinite(y)
  x, y = x[mask], y[mask]
  denom = float(np.dot(x, x))
  slope = float(np.dot(x, y) / denom) if denom != 0.0 else 0.0
  return slope
  ```

- [ ] **Step 3.5: Tree tests**

  - Create test Newick tree strings
  - Verify `total_tree_length` sums correctly
  - Verify `treeness` computation
  - Verify DVMC with mocked `tree.distance()` values
  - Verify saturation with simple MSA+tree (small dataset)
  - Verify RF_distance edge cases (identical trees → 0, disjoint → high)

### Verification for Task 3

```bash
python -m pytest tests/pretree/test_metrics.py -k "test_tree" -v
```

---

## Task 4: Pseudo-Tree Metrics (FastTree)

**Files:** Modify `phyloai/pretree/metrics.py` (FastTree section), Modify `phyloai/core/env.py` (FastTree registry)

- [ ] **Step 4.1: FastTree invocation**

  ```python
  def _compute_pseudo_tree_metrics(msa_path: Path, seq_type: str, fasttree_path: str) -> dict:
  ```

  - AA: `FastTree -lg -noml -boot 500 < msa_file`
  - NT: `FastTree -nt -gtr -noml -boot 500 < msa_file`
  - Pipe MSA content via stdin
  - Parse output Newick tree with `Bio.Phylo`
  - Extract same 10 metrics as `calculate_tree_features` but with `_FT` suffix
  - Handle FastTree failure → return `{... all "NA"}` for each metric

- [ ] **Step 4.2: FastTree in doctor**

  Modify `phyloai/core/env.py`:
  ```python
  TOOL_REGISTRY["fasttree"] = {
      "required": False,
      "version_flag": "",
      "install": "Download from http://www.microbesonline.org/fasttree/",
  }
  ```

- [ ] **Step 4.3: FastTree tests**

  - Mock `subprocess.Popen` to return a test Newick tree
  - Verify `_FT` suffix in all metric names
  - Verify FastTree failure → "NA" values, not crash

### Verification for Task 4

```bash
python -m pytest tests/pretree/test_metrics.py -k "test_fasttree" -v
```

---

## Task 5: Orchestration — `run_metrics`

**Files:** Modify `phyloai/pretree/metrics.py` (orchestrator)

- [ ] **Step 5.1: Per-marker worker function**

  ```python
  def _metric_worker(args: tuple) -> dict:
      """Process a single marker: MSA features + freq + tree + pseudo-tree + saturation.
      Executed in a ProcessPoolExecutor worker."""
  ```

  - One worker call per marker
  - Reads MSA file (if exists) → MSA metrics + freq stats
  - Reads tree file (if exists) → tree metrics + DVMC + RF + saturation (if MSA also exists)
  - Reads pseudo-tree via FastTree (if `--pseudo-tree-metrics` and MSA exists)
  - Returns flat dict of all metrics + `loci` key
  - Exceptions caught, returned as error dict

- [ ] **Step 5.2: Main orchestrator**

  ```python
  def run_metrics(
      msa_dir: Path | None,
      tree_dir: Path | None,
      seq_type: str,
      threads: int,
      output_dir: Path,
      decimal_places: int = 6,
      skip_freq: bool = False,
      pseudo_tree: bool = False,
      skip_pairwise_identity: bool = False,
      outgroup_list: Path | None = None,
      ref_tree: Path | None = None,
      overwrite: bool = False,
      dry_run: bool = False,
      quiet: bool = False,
  ) -> dict:
      """Orchestrate the full metrics pipeline: compute → plot → correlate."""
  ```

  Flow:
  1. Validate inputs (at least one of msa_dir/tree_dir)
  2. Handle `--dry-run`: print plan and exit
  3. Handle `--overwrite`: clear output_dir
  4. Handle output_dir conflict (non-empty → exit 1)
  5. Pass 1: scan MSA headers for taxa pool
  6. Pair files by stem
  7. Check FastTree if `--pseudo-tree-metrics`
  8. **Pairwise identity warning:** if any marker >200 taxa and `--skip-pairwise-identity` is not set, print stderr warning once
  9. Load ref tree if `--ref-tree`
  10. Parallel computation via `ProcessPoolExecutor`
  11. Collect results → `metrics.csv` (rounded to `decimal_places`)
  11. Generate distribution plots (all numeric metrics)
  12. Generate `metrics.basic_statistics.csv`
  13. Generate correlation heatmap + matrix CSV
  14. Write `metrics.log` (resolved command, parameters, per-marker stderr, wall time, exit code)
  15. Write `result.json`
  16. Terminal output (unless `--quiet`): Rich progress bar during computation, Rich summary table, "Metrics saved to ..." / "Plots saved to ..." / "Correlation heatmap saved to ..." / "Results saved to ..."
  17. Return result dict for CLI layer

- [ ] **Step 5.3: CSV field ordering**

  ```python
  _METRICS_CSV_ORDER = [
      "loci", "DataType",
      # MSA: counts
      "num_taxa", "taxa_occupancy", "num_sites",
      "num_patterns", "proportion_patterns",
      "num_parsimony_sites", "proportion_parsimony",
      "num_singletons", "proportion_singletons",
      "num_sites/num_taxa", "num_patterns/num_taxa",
      "num_parsimony_sites/num_taxa", "num_singletons/num_taxa",
      # MSA: proportions
      "proportion_gaps", "proportion_invariant",
      # MSA: diversity
      "entropy", "bollback", "pattern_entropy",
      "rcfv", "nrcfv", "average_pairwise_identity", "GC_content",
      # Tree
      "average_BS", "sd_BS", "total_tree_length",
      "average_internal_branch_length", "sd_internal_branch_length",
      "average_terminal_branch_length", "sd_terminal_branch_length",
      "tree_diameter", "average_patristic_distance", "sd_patristic_distance",
      "evo_rate", "treeness", "dvmc", "saturation", "RF_distance",
  ]
  ```

  Frequency columns + pseudo-tree columns appended dynamically after these.
  `DataType` = `"AA"` or `"NT"` when `--msa-dir` is provided, `""` for tree-only mode.
  MSA-only metrics = `""` when no MSA; tree-only metrics = `""` when no tree.

- [ ] **Step 5.4: Orchestration tests**

  - Integration test: temp dir with 2 MSA files + 2 tree files → verify `metrics.csv` has correct shape
  - Integration test: only `--tree-dir` → verify only tree metrics present
  - Integration test: `--dry-run` → verify no files created, correct plan output
  - Integration test: empty output dir conflict → exit code 1

### Verification for Task 5

```bash
python -m pytest tests/pretree/test_metrics.py -k "test_run_metrics" -v
```

---

## Task 6: Distribution Plot Generation

**Files:** Modify `phyloai/pretree/metrics.py` (plot section)

- [ ] **Step 6.1: Single metric plot function**

  ```python
  def _plot_single_metric(
      data: np.ndarray,
      metric_name: str,
      output_path: Path,
      bins: int = 50,
      xmin: float | None = None,
      xmax: float | None = None,
      tukey_k: float | None = None,
      title: str | None = None,
      xlabel: str | None = None,
      ylabel: str = "Density",
      color: str = "#2E86AB",
      fig_width: float = 10.0,
      fig_height: float = 8.0,
      dpi: int = 150,
      font_size: int = 12,
  ) -> None:
      """Generate a density histogram PDF for a single metric."""
  ```

  - Use `matplotlib` (no seaborn dependency for simple histograms)
  - Clean data: remove NaN
  - Apply Tukey filter if `tukey_k` is set
  - Histogram with density normalization + KDE density curve overlay
  - Green-ish fill (`#D4EDDA` default), orange density line (`#FF8C00`) matching R script aesthetics
  - X-axis lims: auto or user-specified
  - Title, labels, font sizes configurable
  - Save as PDF via `fig.savefig(output_path, dpi=dpi, bbox_inches="tight")`
  - Close figure to avoid memory leaks

- [ ] **Step 6.2: Batch plot generation (all metrics)**

  ```python
  def _generate_all_plots(
      rows: list[dict],
      numeric_cols: list[str],
      plots_dir: Path,
      bins: int = 50,
  ) -> int:
      """Generate one PDF per numeric metric. Return count of generated plots."""
  ```

  - Iterate `numeric_cols` (exclude `loci`, `DataType`)
  - Extract values: `np.array([float(r[col]) for r in rows if r.get(col) not in (None, "")], dtype=float)`
  - Call `_plot_single_metric` for each
  - Output: `plots/<metric_name>.pdf`
  - Skip columns where all values are NaN

- [ ] **Step 6.3: Basic statistics CSV generation**

  ```python
  def _generate_basic_statistics(
      rows: list[dict],
      numeric_cols: list[str],
      output_path: Path,
  ) -> None:
  ```

  For each numeric column, compute: mean, median, min, max, q25, q75, std, n_ex_NA, n_total using `numpy`.
  Write as CSV with `metric` as the key column using `csv.DictWriter`.

- [ ] **Step 6.4: Plot tests**

  - Test `_plot_single_metric` with mock data → verify PDF file created
  - Test with `tukey_k=1.5` → verify outliers excluded from visualization
  - Test `_generate_basic_statistics` → verify correct values in CSV
  - Test with empty data → graceful handling

### Verification for Task 6

```bash
python -m pytest tests/pretree/test_metrics.py -k "test_plot" -v
```

---

## Task 7: Correlation Analysis

**Files:** Modify `phyloai/pretree/metrics.py` (correlation section)

- [ ] **Step 7.1: Correlation computation**

  ```python
  def _compute_correlation(
      rows: list[dict],
      columns: list[str],
      method: str = "spearman",
  ) -> tuple[np.ndarray, list[str]]:
      """Compute pairwise correlation matrix for selected columns.
      Returns (M×M ndarray, column_names)."""
  ```

  - Extract numeric values per column: `data = np.array([[float(r[col]) for r in rows] for col in columns], dtype=float).T`
  - Drop rows with any NaN (complete case analysis) via `np.isnan(data).any(axis=1)`
  - Compute correlation matrix using `scipy.stats.spearmanr` pairwise
  - Return M×M `np.ndarray` + ordered column names

- [ ] **Step 7.2: Ordered heatmap generation**

  ```python
  def _generate_correlation_heatmap(
      corr_matrix: np.ndarray,
      col_names: list[str],
      output_path: Path,
      triangle: str = "full",
      cluster_rectangles: int | None = None,
      cmap: str = "RdBu_r",
      annot: bool = True,
      fmt: str = ".2f",
      fig_width: float = 12.0,
      fig_height: float = 10.0,
      dpi: int = 150,
      font_size: int = 10,
      title: str = "Correlation Heatmap",
  ) -> None:
  ```

  - Compute Ward leaf order from distance matrix = `1 - |correlation|` for readability, but do not draw dendrogram axes.
  - Draw the heatmap/corrplot with plain Matplotlib axes so the output has no top/left dendrogram whitespace.
  - Use a single colorbar: left side for `--triangle upper`, right side for `full` and `lower`, with enough padding to avoid overlap with metric labels.
  - Apply `--triangle` mask: for "lower" or "upper", omit the opposite triangle cells and grid segments after ordering.
  - For `--triangle upper`, place x labels on top and y labels on the right; for `lower`, keep labels on left/bottom.
  - Expose `--label-angle` to control x-axis metric label rotation; default 45 degrees.
  - Hide rectangular spines in triangle modes and draw a stepped triangular border outside the diagonal cells.
  - Apply `--cluster-rectangles` only when `--triangle full`; ignore it for lower/upper triangle modes and warn when ignored.
  - Save as PDF

  ```python
  # Magnitude-based distance for clustering
  from scipy.cluster.hierarchy import linkage
  from scipy.spatial.distance import squareform

  dist = 1.0 - np.abs(corr_matrix)
  # Ensure diagonal is exactly 0
  np.fill_diagonal(dist, 0.0)
  condensed = squareform(dist, checks=False)
  linkage_matrix = linkage(condensed, method="ward")
  ```

  - Use Matplotlib axes directly for corrplot-style circle cells; no seaborn dependency.

- [ ] **Step 7.3: Correlation matrix CSV output**

  ```python
  def _write_correlation_csv(corr_matrix: np.ndarray, col_names: list[str], output_path: Path) -> None:
  ```

  Write the full correlation matrix to CSV using `csv.writer`, with variable names as the first row and first column.

- [ ] **Step 7.4: Correlation tests**

  - Test with known data → verify Spearman values match `scipy.stats.spearmanr`
  - Test `triangle="lower"` → verify upper triangle masked
  - Test magnitude-based clustering → verify strong positive and negative correlations are adjacent
  - **Edge case: constant column** (all values identical) → Spearman returns NaN for that pair → exclude from clustering, mark as NaN in matrix, print warning
  - **Edge case: all-NA column** → excluded before correlation computation, logged
  - **Edge case: only 1 valid column** → exit with error "Need at least 2 numeric columns with non-NA values"
  - **Edge case: only 2 variables** → generate matrix + heatmap without clustered reordering
  - **Edge case: all pairwise NaN** (no complete cases) → exit with error "No complete cases for correlation analysis"

### Verification for Task 7

```bash
python -m pytest tests/pretree/test_metrics.py -k "test_correlate" -v
```

---

## Task 8: CLI Registration

**Files:** Modify `phyloai/cli/commands/pretree.py`

- [ ] **Step 8.1: Update `_PretreeGroup.list_commands`**

  ```python
  return ["convert", "stats", "align", "trim", "metrics", "concat"]
  ```

- [ ] **Step 8.2: Register `metrics` command**

  ```python
  @pretree.command(
      "metrics",
      help="Compute MSA and tree metrics, generate distribution plots and correlation heatmap.",
  )
  @click.option("--msa-dir", type=click.Path(exists=True, file_okay=False, path_type=Path), help="MSA directory.")
  @click.option("--tree-dir", type=click.Path(exists=True, file_okay=False, path_type=Path), help="Gene tree directory.")
  @click.option("--seq-type", type=click.Choice(["AA", "NT", "auto"]), default="auto", show_default=True)
  @click.option("--outgroup-list", type=click.Path(exists=True, dir_okay=False, path_type=Path))
  @click.option("--ref-tree", type=click.Path(exists=True, dir_okay=False, path_type=Path))
  @click.option("--skip-freq-statistics", is_flag=True, default=False)
  @click.option("--pseudo-tree-metrics", is_flag=True, default=False)
  @click.option("--skip-pairwise-identity", is_flag=True, default=False, help="Skip average_pairwise_identity (slow for >200 taxa).")
  @click.option("--round", "decimal_places", type=click.IntRange(0, 12), default=6, show_default=True, help="Decimal places for metric values in CSV.")
  @click.option("--table-format", type=click.Choice(["csv", "tsv"]), default="csv", show_default=True, help="Table format for auxiliary tabular outputs.")
  @click.option("--output-dir", "-o", type=click.Path(file_okay=False, path_type=Path), default=Path("runs/pretree/metrics"), show_default=True)
  @click.option("--threads", "-t", type=int, default=4, show_default=True)
  @click.option("--dry-run", is_flag=True, default=False, help="Validate inputs and show plan without computing.")
  @click.option("--overwrite", is_flag=True, default=False)
  @click.option("--quiet", "-q", is_flag=True, default=False)
  def metrics_command(...):
  ```

  Validation:
  - At least one of `--msa-dir`, `--tree-dir` must be provided
  - `--pseudo-tree-metrics` requires `--msa-dir`
  - `--outgroup-list` requires `--tree-dir`
  - `--ref-tree` requires `--tree-dir`
  - `--threads` ≥ 1
  - **Output dir conflict:** applies the standard policy (non-empty → exit 1 unless `--overwrite`)

- [ ] **Step 8.3: Register `plot` and `correlate` as Click group subcommands**

  `metrics` is a Click group (`@click.group("metrics", invoke_without_command=True)`)
  registered on `pretree` via `pretree.add_command(metrics_group)`. The group's
  main callback (invoked when no subcommand is given) runs the full 3-step pipeline:
  compute → plot → correlation.

  Subcommands are decorated with `@metrics_group.command("plot")` and
  `@metrics_group.command("correlate")` — they are NOT flat `pretree` commands.

  **`metrics plot` options:**
  - `--csv` (required) — existing metrics.csv path
  - `--input-format csv|tsv|auto` (default auto) — input table format with content-based auto-detection
  - `--metric` (required) — column name to plot
  - `--bins` (default 50), `--xmin`, `--xmax`, `--tukey-k` (optional, saves CSV)
  - `--title`, `--xlabel`, `--ylabel`, `--color`, `--fig-width/height`, `--dpi`, `--font-size`
  - `--output-dir` — defaults to `<csv_parent>/plot_<metric>/`
  - When `--tukey-k` is set, filtered-out loci are saved as `<output_dir>/<metric>.tukey_filtered.csv` with columns `loci,value`.

  **`metrics correlate` options:**
  - `--csv` (required) — existing metrics.csv path
  - `--input-format csv|tsv|auto` (default auto) — input table format with content-based auto-detection
  - `--metrics` (optional) — comma-separated subset; `all` means every numeric column; omitted means automatic core-metric selection
  - `--include-freq` — include `freq*` columns in automatic selection
  - `--include-sd` — include `sd_*` columns in automatic selection
  - `--method` — `spearman` (default, rank-based) or `pearson` (z-score normalized)
  - `--triangle` — `full` (default), `lower`, `upper`
  - `--annot/--no-annot` — show correlation values in cells; default `--no-annot`
  - `--cluster-rectangles` — draw N rectangles on full matrices only; ignored with warning for lower/upper triangle modes
  - `--cmap`, `--fmt`, `--fig-width/height`, `--dpi`, `--font-size`, `--label-angle`, `--title`
  - `--output-dir` — defaults to `runs/pretree/metrics/correlate/`

  Circle cells use `matplotlib.patches.Circle` with radius proportional to √|correlation|,
  capped at 0.45 (cell boundary). Data is reordered by Ward leaf order
  (`scipy.cluster.hierarchy.leaves_list`) before drawing, but no dendrogram is rendered.
  Title is placed via the heatmap axis title to avoid excess whitespace above the matrix.

  Automatic column selection must exclude `loci`, `DataType`, all `freq*`, and all
  `sd_*` columns by default so the default PDF is readable. Explicit `--metrics`
  always wins, including frequency and standard-deviation columns when named.
  `--metrics all` must include all numeric columns.

- [ ] **Step 8.4: Import and wire-up**

  Import `metrics_group` into `phyloai/cli/commands/pretree.py` and register via
  `pretree.add_command(metrics_group)`.

### Verification for Task 8

```bash
python -m phyloai pretree metrics --help
python -m phyloai pretree metrics plot --help
python -m phyloai pretree metrics correlate --help
```

---

## Task 9: Complete Integration Test

**Files:** Modify `tests/pretree/test_metrics.py`

- [ ] **Step 9.1: End-to-end test with small real data**

  Create temp directory with:
  - 2 small AA MSA FASTA files (3–5 taxa, ~50 sites each)
  - 2 corresponding Newick tree files
  - Verify pipeline produces:
    - `metrics.csv` with correct row count and all expected columns
    - `plots/*.pdf` — one PDF per numeric metric
    - `metrics.basic_statistics.csv` — correct shape
    - `correlation_heatmap.pdf` — produced
    - `correlation_matrix.csv` — correct dimension
    - `result.json` — `status: success`, correct `key_results.n_markers`
    - `metrics.log` — exists

- [ ] **Step 9.2: Test `--dry-run`**

  Verify no output files created, terminal output shows planned actions.

- [ ] **Step 9.3: Test `--tree-dir` only mode**

  Verify only tree metrics in CSV, no MSA metrics.

- [ ] **Step 9.4: Test `--overwrite`**

  Create output dir with content → run with `--overwrite` → verify old content replaced.

- [ ] **Step 9.5: Test `--skip-freq-statistics`**

  Verify no `freq*` columns in CSV.

- [ ] **Step 9.6: Test taxon mismatch warnings**

  Create MSA with taxa {A,B,C}, tree with taxa {A,B,D} → verify warning in `result.json`.

### Verification for Task 9

```bash
python -m pytest tests/pretree/test_metrics.py -v
```

---

## Task 10: Command Documentation

**Files:** Create `docs/commands/pretree-metrics.md`, Modify `README.md`

Per total design Section 4.5, new or materially changed commands must update both their command document and the top-level README command index.

- [ ] **Step 10.1: Create `docs/commands/pretree-metrics.md`**

  Follow the required section structure from the total design:
  - **Purpose:** what `pretree metrics` does and does not do
  - **Usage:** minimal usage examples for `metrics`, `metrics plot`, `metrics correlate` + full parameter tables
  - **Inputs:** supported formats, directory scanning rules, tree-file suffix stripping logic
  - **Outputs:** `metrics.csv` schema (all column names with explanations), `plots/` directory, `correlation_heatmap.pdf`, `result.json` schema summary
  - **Examples:** AA-only, NT-with-trees, `--pseudo-tree-metrics`, re-plotting with `--tukey-k` (filtered loci in `<metric>.tukey_filtered.csv`), correlation with triangle display options
  - **Warnings and Errors:** taxon mismatches, unpaired files, FastTree not found, large-dataset pairwise identity warning
  - **Notes:** relationship to `pretree filter`, metric name glossary (what each metric means)

- [ ] **Step 10.2: Update README command index**

  Add `pretree metrics` entry to the top-level README command index table with brief description and link to `docs/commands/pretree-metrics.md`.

### Verification for Task 10

```bash
# Verify docs file exists and is valid markdown
head -30 docs/commands/pretree-metrics.md
```

---

## Implementation Order

```
Task 1  → Task 2  → Task 5 (orchestration)
Task 1  → Task 3  → Task 5
Task 1  → Task 4  → Task 5
                       ↓
                  Task 6 (plot)
                       ↓
                  Task 7 (correlate)
                       ↓
                  Task 8 (CLI)
                       ↓
                  Task 9 (integration tests)
                       ↓
                  Task 10 (docs)
```

Tasks 2, 3, 4 can be developed in parallel after Task 1 is complete.

---

## Commit policy

No per-step commits. A single commit is made at the end when all tests pass, on explicit user instruction.
