# pretree metrics Design Specification

**Date:** 2026-06-14  
**Status:** Approved for implementation  
**Reference:** `ref/scripts/extract_msa_tree_features.py`, `ref/scripts/MSA_Tree_metrics_exploration_20250822_data.table.R`, `docs/superpowers/specs/2026-06-07-phyloai-design.md`

---

## 1. Purpose

`phyloai pretree metrics` computes molecular marker attributes from MSAs and/or gene trees for statistical exploration and downstream filtering. It operates **after** `pretree trim` and **before** `pretree filter` in the standard pipeline.

The module has three layers:

1. **Compute** — parallel calculation of MSA + tree metrics → `metrics.csv`
2. **Plot** — distribution visualization per metric → one PDF per metric + basic statistics CSV
3. **Correlate** — Spearman correlation heatmap from `metrics.csv`

Steps 1–3 execute by default under `phyloai pretree metrics`. Steps 2 and 3 can also be re-run independently via `phyloai pretree metrics plot` and `phyloai pretree metrics correlate` against an existing `metrics.csv`.

What it does **not** do:
- UMAP clustering (moved to `pretree filter`)
- Marker-level filtering decisions (use `pretree filter`)
- Tree inference (use `tree genetree`)

---

## 2. Architecture

```
phyloai/pretree/metrics.py           # core library: compute, plot, correlate, helpers
phyloai/cli/commands/pretree.py      # CLI: metrics Click group, plot + correlate subcommands
```

**Reused components:**
- `core/runner.py` — external tool invocation (FastTree)
- `core/env.py` — tool detection (FastTree in `TOOL_REGISTRY`)
- `core/schema.py` — `MSACollection`, shared types
- `core/sequence_normalization.py` — `gap_chars`, `standard_chars`, `resolve_seq_type`
- `core/formats.py` — format detection

**External dependencies:**
- `BioPython` (AlignIO, Phylo, SeqIO) — already in project
- `numpy` — already in project
- `matplotlib` — for distribution plots and correlation heatmaps
- `scipy` — Spearman/Pearson helpers and Ward leaf ordering for correlation plots

**Not required:** `phykit` — all computations (RCFV, saturation, DVMC) are implemented directly in `pretree/metrics.py`.

---

## 3. Command Structure

```
phyloai pretree metrics
    --msa-dir <path>               # MSA directory
    [--tree-dir <path>]            # gene tree directory
    [--seq-type AA|NT|auto]        # molecule type, default auto
    [--outgroup-list <file>]       # outgroup taxa for DVMC pruning
    [--ref-tree <file>]            # reference species tree for nRF distance
    [--skip-freq-statistics]       # skip per-character frequency computation
    [--pseudo-tree-metrics]        # compute FastTree-derived pseudo-tree metrics
    [--skip-pairwise-identity]     # skip average_pairwise_identity to save time on large datasets
    [--round N]                    # decimal places in metrics.csv, default 6, range 0-12
    [--output-dir <dir>]           # output directory, default runs/pretree/metrics/
    [--threads N]                  # parallel workers, default 4
    [--dry-run]                    # validate inputs and show plan without computing
    [--overwrite]                  # overwrite existing output directory
    [--quiet]                      # suppress terminal output

At least one of --msa-dir or --tree-dir must be provided.
`--resume` is not supported for metrics — the module is a deterministic, read-only
analysis (same inputs → identical outputs), so re-running is simpler and produces
the same result as resuming from a checkpoint.

phyloai pretree metrics plot
    --csv <metrics.csv>            # metrics CSV file (required)
    --metric <col>                 # single metric column to plot (required)
    [--bins N]                     # histogram bins, default 50
    [--xmin FLOAT]                 # x-axis minimum
    [--xmax FLOAT]                 # x-axis maximum
    [--tukey-k FLOAT]              # Tukey's Fences k value (e.g. 1.5 = standard, 3.0 = conservative)
    [--title TEXT]                 # plot title (default: "Distribution of <metric>")
    [--xlabel TEXT]                # x-axis label (default: metric display name)
    [--ylabel TEXT]                # y-axis label, default "Density"
    [--color HEX]                  # bar fill colour, default "#2E86AB"
    [--fig-width FLOAT]            # figure width (inches), default 10.0
    [--fig-height FLOAT]           # figure height (inches), default 8.0
    [--dpi INT]                    # output resolution, default 150
    [--font-size INT]              # base font size, default 12
    [--output-dir <dir>]           # default <csv_parent>/plot_<metric>/
    [--overwrite] [--quiet]

When --tukey-k is set, outlier loci are saved to
`<output_dir>/<metric>.tukey_filtered.csv` with columns `loci,value`.

phyloai pretree metrics correlate
    --csv <metrics.csv>            # metrics CSV file (required)
    [--metrics <col1,col2,...>|all]# explicit subset; default auto-selects core numeric metrics
    [--include-freq]               # include freq* columns in automatic metric selection
    [--include-sd]                 # include sd_* columns in automatic metric selection
    [--method spearman|pearson]    # correlation method, default spearman
    [--triangle full|lower|upper]  # matrix display mode, default full
    [--annot / --no-annot]         # show numeric values in cells, default --no-annot
    [--cluster-rectangles N]       # number of cluster rectangles; full triangle only
    [--cmap NAME]                  # matplotlib colormap, default RdBu_r
    [--fmt FORMAT]                 # annotation format, default .2f
    [--fig-width FLOAT]            # figure width (inches), default 12.0
    [--fig-height FLOAT]           # figure height (inches), default 10.0
    [--dpi INT]                    # output resolution, default 150
    [--font-size INT]              # base font size, default 10
    [--label-angle FLOAT]          # x-axis label rotation angle, default 45
    [--title TEXT]                 # plot title (default: none)
    [--output-dir <dir>]           # default runs/pretree/metrics/correlate/
    [--overwrite] [--quiet]
```

Circle cells are drawn with `matplotlib.patches.Circle`, area proportional to
|correlation|, radius capped at 0.45 (cell boundary). Data is reordered by Ward
leaf order for readability, but no dendrogram is drawn. Clustering distance =
`1 - |correlation|` (magnitude-based, not sign-based). Pearson method applies
z-score normalization before computation.

---

## 4. Metrics Specification

### 4.1 MSA Metrics

Metrics are ordered logically: core counts first, then pairwise ratios (`num_XXX` / `proportion_XXX`), then single-entity proportions/statistics, then derived diversity metrics.

| #  | Metric                    | Formula / Source                                                                         | New?   |
|----|---------------------------|------------------------------------------------------------------------------------------|--------|
| 1  | `num_taxa`                | `len(msa)`                                                                                | —      |
| 2  | `taxa_occupancy`          | `num_taxa / total_taxa_pool` (union of all taxa across all input MSAs)                    | **NEW** |
| 3  | `num_sites`               | `msa.get_alignment_length()`                                                              | —      |
| 4  | `num_patterns`             | `len(set(normalized columns))` — non-standard chars → `-` before hashing (stats.py logic) | —      |
| 5  | `proportion_patterns`      | `num_patterns / num_sites`                                                                | —      |
| 6  | `num_parsimony_sites`      | sites where ≥2 chars each with count ≥2                                                   | —      |
| 7  | `proportion_parsimony`     | `num_parsimony_sites / num_sites`                                                         | —      |
| 8  | `num_singletons`           | variable but not parsimony-informative sites (stats.py logic)                             | **NEW** |
| 9  | `proportion_singletons`    | `num_singletons / num_sites`                                                              | **NEW** |
| 10 | `num_sites/num_taxa`      | —                                                                                        | —      |
| 11 | `num_patterns/num_taxa`   | —                                                                                        | —      |
| 12 | `num_parsimony_sites/num_taxa` | —                                                                                   | —      |
| 13 | `num_singletons/num_taxa` | —                                                                                        | **NEW** |
| 14 | `proportion_gaps`         | `total_gaps / (nsites × ntaxa)`                                                          | —      |
| 15 | `proportion_invariant`    | `invariant_sites / nsites`                                                                | —      |
| 16 | `entropy`                 | mean per-site Shannon entropy (valid chars only)                                          | —      |
| 17 | `bollback`                | Bollback multinomial statistic                                                            | —      |
| 18 | `pattern_entropy`         | pattern entropy                                                                           | —      |
| 19 | `rcfv`                    | RCFV (frequency-based, Fleming & Struck 2023 Eq. 1)                                       | —      |
| 20 | `nrcfv`                   | nRCFV = `RCFV / (p^(-0.5) · n^(0.01) · c) · 100` (Fleming & Struck 2023 Eq. 4)           | **NEW** |
| 21 | `average_pairwise_identity` | mean pairwise identity across all taxon pairs                                            | —      |
| 22 | `GC_content`              | `(G + C) / total_valid_chars` — **NT only**, `""` for AA                                  | **NEW** |

**Key adjustments from `extract_msa_tree_features.py`:**
- `num_patterns` (item 4): non-standard characters → gap `-` before hashing (stats.py logic). Since `pretree convert` normalises sequences, this primarily affects downstream consistency.
- `proportion_patterns`, `proportion_parsimony`, `proportion_singletons` (items 5, 7, 9): renamed from `num_XXX/num_sites` for clarity.
- `rcfv` (item 19): classic formula as in extract_msa_tree_features.py and Fleming & Struck 2023 Eq. 1. Has known biases vs sequence length and taxa count. **Gaps are excluded** from frequency calculations.
- `nrcfv` (item 20): the bias-corrected metric from Fleming & Struck 2023 Eq. 4, where *p* = num_sites, *n* = num_taxa, *c* = 4 (NT) or 20 (AA). This metric is comparable across datasets of different sizes.
- **"Valid characters" definition** applies to: `num_patterns` (non-standard→gap), `num_parsimony_sites`, `num_singletons`, `proportion_invariant`, `entropy`, `rcfv`, `nrcfv`, `average_pairwise_identity`. Gap `-` is never counted as a valid character state for the purpose of these metrics. `proportion_invariant` counts sites where ≤1 valid character state appears (a site of all gaps → 0 valid chars → invariant). Reference lines: `extract_msa_tree_features.py:251-267`, `stats.py:135-153`.

`taxa_occupancy` (item 2) requires two-pass processing: first pass gathers all taxon names across all MSAs to determine the total pool, second pass computes per-marker metrics. Pass 1 reads only FASTA headers.

### 4.2 Tree Metrics

| #  | Metric                            | Formula                                                              | New?   |
|----|-----------------------------------|----------------------------------------------------------------------|--------|
| 1  | `average_BS`                      | mean bootstrap support (internal nodes with `.confidence`)            | —      |
| 2  | `sd_BS`                           | stdev of bootstrap support                                            | —      |
| 3  | `total_tree_length`               | sum of all branch lengths                                             | —      |
| 4  | `average_internal_branch_length`  | mean of internal branch lengths                                       | —      |
| 5  | `sd_internal_branch_length`       | stdev of internal branch lengths                                      | —      |
| 6  | `average_terminal_branch_length`  | mean of terminal branch lengths                                       | —      |
| 7  | `sd_terminal_branch_length`       | stdev of terminal branch lengths                                      | —      |
| 8  | `tree_diameter`                   | max patristic distance between any two taxa                           | —      |
| 9  | `average_patristic_distance`      | mean patristic distance across all taxon pairs                        | —      |
| 10 | `sd_patristic_distance`           | stdev of patristic distance                                           | —      |
| 11 | `evo_rate`                        | `total_tree_length / n_taxa`                                          | —      |
| 12 | `treeness`                        | `sum(internal_brlen) / total_tree_length`                             | —      |
| 13 | `dvmc`                            | stddev of root-to-tip distances; outgroups pruned if `--outgroup-list` | —      |
| 14 | `saturation`                      | phykit saturation slope (1st output value) — requires MSA + tree     | **NEW** |
| 15 | `RF_distance`                     | normalised Robinson-Foulds distance — only when `--ref-tree` provided | —      |

**saturation** (item 14) requires both the MSA and the tree for each marker. Output value is `slope` (not `abs(1 - slope)`). Computation is O(n²) in taxa count; for large trees, this is the dominant cost. See Section 5.5 for implementation details.

**dvmc** (item 13) with `--outgroup-list`: outgroup taxa are pruned from the tree before DVMC computation. Without outgroups, DVMC is computed on the full tree. Both paths align with phykit's `dvmc.py` logic.

**RF_distance** (item 15) is only computed when `--ref-tree <file>` is provided. Trees are pruned to shared tips, rooted with the same tip, and normalized RF = `plain_rf / (2 × (tip_count - 3))`.

### 4.3 Frequency Statistics (optional)

When `--skip-freq-statistics` is **not** set, per-character frequencies are computed:

- **NT:** `freqA`, `freqC`, `freqG`, `freqT`
- **AA:** `freqA`, `freqR`, `freqN`, `freqD`, `freqC`, `freqQ`, `freqE`, `freqG`, `freqH`, `freqI`, `freqL`, `freqK`, `freqM`, `freqF`, `freqP`, `freqS`, `freqT`, `freqW`, `freqY`, `freqV`

Frequencies sum to 1.0 across valid characters (gaps excluded).

### 4.4 Pseudo-Tree Metrics (optional)

When `--pseudo-tree-metrics` is set:

| Metric                                 |
|----------------------------------------|
| `average_BS_FT`                        |
| `sd_BS_FT`                             |
| `total_tree_length_FT`                 |
| `average_internal_branch_length_FT`    |
| `sd_internal_branch_length_FT`         |
| `average_terminal_branch_length_FT`    |
| `sd_terminal_branch_length_FT`         |
| `tree_diameter_FT`                     |
| `average_patristic_distance_FT`        |
| `sd_patristic_distance_FT`             |

FastTree invocation: `FastTree [-nt -gtr | -lg] -noml -boot 500 < msa_file`. The `_FT` suffix distinguishes pseudo-tree metrics from actual gene-tree metrics.

**FastTree detection:** `phyloai doctor` must check for FastTree availability. Running `FastTree` without arguments prints a version line and help; the first line contains the version string. Added to `core/env.py` `TOOL_REGISTRY` under `fasttree`.

### 4.5 Data Type Column

The `metrics.csv` always includes a `DataType` column as the second field (after `loci`):
- `AA` or `NT` when `--msa-dir` is provided
- `""` (empty string) when only `--tree-dir` is provided (tree-only mode)

This enables downstream stratified analysis while guaranteeing column position stability regardless of input mode.

---

## 5. Step 1: Metric Computation

### 5.1 Workflow

```
Optional: --msa-dir          Optional: --tree-dir
        │                           │
        ▼                           ▼
   Pair by file prefix ────> Check taxon name consistency (warn on mismatch)
        │
        ▼
   Parallel computation (ProcessPoolExecutor, --threads workers)
        │
        ├── Per-marker: MSA features + freq stats + (pseudo-tree) + (tree) + (saturation)
        │
        ▼
   Aggregate → metrics.csv + result.json
```

**Input requirements:** At least one of `--msa-dir` or `--tree-dir` must be provided. If only `--tree-dir` is given, only tree-level metrics are computed (no MSA metrics, frequency stats, pseudo-tree metrics, saturation). If only `--msa-dir` is given, MSA metrics + optional freq stats + optional pseudo-tree metrics are computed.

### 5.2 File Pairing

This command follows the global file-matching policy in `docs/superpowers/specs/2026-06-07-phyloai-design.md`.

1. Scan `--msa-dir` for alignment inputs.
2. If `--tree-dir` provided, scan for tree inputs.
3. Derive MSA logical loci by taking the filename before the final `.`.
4. For each tree file, try one-suffix and two-suffix reduction candidates.
5. If exactly one candidate matches an available MSA locus, pair it.
6. If neither candidate matches, record the file under `result.json.data.unpaired` and continue.
7. If both candidates match different loci, raise an explicit ambiguity error and stop.
8. If `--msa-dir` + `--tree-dir` but zero pairs remain after matching, exit code 1.

Illustrative examples:

- `gene1.fa.treefile` -> candidates `gene1.fa`, `gene1`
- `gene2.tre` -> candidate `gene2`
- `gene3.FASTA.bestTree` -> candidates `gene3.FASTA`, `gene3`

### 5.3 MSA–Tree Taxon Consistency Check

For each paired marker, verify `set(msa_taxon_names) == set(tree_taxon_names)`:

- **Exact match** → no warning
- **MSA ⊂ tree** (tree has extra taxa) → `[WARN]` with extra taxa listed
- **Tree ⊂ MSA** (MSA has extra taxa) → `[WARN]` with extra taxa listed
- **Partial overlap** → `[WARN]` with symmetric difference listed

Warnings are collected in `result.json.data.taxon_mismatches`. Computation proceeds regardless — filtering decisions belong to the user.

### 5.4 Two-Pass Taxa Pool

`taxa_occupancy` requires the total taxa pool. Processing:

1. **Pass 1:** Quick scan of all MSAs to collect taxon names → union set `total_taxa_pool`
2. **Pass 2:** Parallel computation of all per-marker metrics, using `total_taxa_pool` for `taxa_occupancy`

Pass 1 reads only headers — cheap.

### 5.5 Saturation Computation

Saturation (item 14 in tree metrics) is implemented directly in `phyloai/pretree/metrics.py` — no dependency on phykit. The computation follows phykit's algorithm (`phykit/services/tree/saturation.py:43-70`):

1. Read the MSA and tree for the marker.
2. For all taxon pairs (i, j), compute:
   - Patristic distance from the tree (`tree.distance(taxon_i, taxon_j)`)
   - Uncorrected distance = `1 - (identity / alignment_length)` using ALL sites by default (`exclude_gaps=False`, matching phykit's default — not excluding gaps)
3. Fit a line through origin: `slope = Σ(x·y) / Σ(x²)` where `x = patristic, y = uncorrected`.
4. Output value: `slope` (not `abs(1 - slope)`).

When `exclude_gaps` is set (optional), positions where EITHER sequence has a gap character are excluded. Gap characters are detected via `core.sequence_normalization.gap_chars()` (includes `-`, `?`, `N` for NT; `-`, `?` for AA).

Only computed when both MSA and tree exist for a given marker. For large taxa counts (>~64), this is O(n²) pairwise and the dominant cost.

### 5.6 Parallel Execution

- `ProcessPoolExecutor` with `max_workers = --threads`
- Each worker processes one marker end-to-end (MSA read → metrics → tree read → metrics → saturation)
- Rich progress bar showing `[N/total] marker_name`
- Failed markers are recorded in `result.json.data.failed` with error message; successful markers proceed

### 5.7 Output Files

| File            | Content                                                                 |
|-----------------|-------------------------------------------------------------------------|
| `metrics.csv`   | One row per marker, columns = all computed metrics + `DataType`          |
| `metrics.log`   | Per-step log (see Section 5.8)                                          |
| `result.json`   | Structured result (see Section 8)                                       |

### 5.8 Terminal Display and Logging

Follows total design Section 9.6 requirements:

**Terminal output** (unless `--quiet`):
- Rich progress bar during parallel computation: `[N/total] marker_name`
- Rich summary table after completion (n_markers, n_metrics, n_taxon_mismatches)
- Explicit file paths: `Metrics saved to runs/pretree/metrics/metrics.csv`, `Plots saved to runs/pretree/metrics/plots/`, etc.
- Warnings printed for taxon mismatches and unpaired files

**Log file** (`metrics.log`):
```
# phyloai pretree metrics --msa-dir ... --threads 4
# Started: 2026-06-14T10:00:00
# Parameters: seq_type=AA, skip_freq=False, pseudo_tree=False
# Tool: FastTree 2.1.11 (if pseudo-tree enabled)
# --- per-marker stderr/warnings collected below ---
# Wall time: 142.3s
# Exit code: 0
```

**`--dry-run` output** (terminal only, no files):
```
[DRY RUN] phyloai pretree metrics
Would process 150 markers (MSA only)
  MSA dir: /path/to/trimmed
  Tree dir: (not provided)
  Output: runs/pretree/metrics/
  Metrics: 22 MSA + 20 frequency
  Plots: 42 PDFs → runs/pretree/metrics/plots/
  Correlation: Spearman heatmap → runs/pretree/metrics/correlate/correlation_heatmap.pdf
```

### 5.9 Numeric Precision

All metric values in `metrics.csv` are rounded to a configurable number of decimal places, controlled by `--round N` (default 6, range 0–12). This follows the project convention established in `stats.py` (`_round(value, 6)`) and `extract_msa_tree_features.py` (`format_float`). Six decimal places balance precision with readability; values like `proportion_gaps` are meaningful at 6 digits while very small numbers (e.g., `rcfv` ~0.001) remain precise.

### 5.10 Computational Bottlenecks and Optimizations

Two metrics dominate runtime due to O(n_taxa² × n_sites) pairwise operations:

| Metric | Complexity | Example cost (100 taxa × 1000 sites) |
|---|---|---|
| `average_pairwise_identity` | O(n² × L) | 4,950 pairs × 1000 = ~5M char comparisons |
| `saturation` | O(n² × L) | Same + `tree.distance()` per pair |

**Optimizations applied:**
1. **Numpy vectorization** — MSA sequences converted to `np.array([list(seq)])` of dtype `U1`, pairwise identity computed via broadcasted equality checks (`(seq_i == seq_j).sum()`). This is ~10–50× faster than Python string iteration.
2. **Gap mask precomputation** — For saturation with gap exclusion, gap positions are precomputed as boolean masks per sequence; pairwise `valid = ~(gap_i | gap_j)` is a single numpy operation.
3. **`tree.distance()` caching** — `Bio.Phylo` computes `tree.distance()` via `tree.depths()` which walks the tree each call. For repeated calls (all pairs), preload the distance matrix or use `Bio.Phylo`'s internal caching. If not available, compute all root-to-tip distances once and reconstruct pairwise distances via `dist(i,j) = rt_dist[i] + rt_dist[j] - 2×lca_dist`.

**`--skip-pairwise-identity`:** When set, `average_pairwise_identity` is not computed for any marker (value = `""` in CSV). For large datasets (>200 taxa per marker), this avoids ~50K+ pairwise comparisons per marker.

**Runtime warning:** Before computation begins, if any marker exceeds 200 taxa and `--skip-pairwise-identity` is NOT set, a warning is printed to stderr:

```
[WARN] 5 markers have >200 taxa (max: 350). average_pairwise_identity computes O(n²) pairs
       and may be slow. Use --skip-pairwise-identity to skip this metric.
```

This warning is printed once per run (not per marker), before the progress bar starts. Computation proceeds normally unless the user interrupts.

---

## 6. Step 2: Distribution Plots (`metrics plot`)

### 6.1 Default Behavior (in `phyloai pretree metrics`)

After computing `metrics.csv`, generate one PDF per numeric metric:

- **Filename:** `<metric_name>.pdf` under `runs/pretree/metrics/plots/`
- **Plot style:** histogram (density-normalised) + density curve overlay, matching the R script's style
- **Bins:** 50 (default)
- **Title:** `Distribution of <metric_name>`
- **Axes:** x = metric name, y = density

A `metrics.basic_statistics.csv` is also written to `runs/pretree/metrics/`.

### 6.2 Standalone Subcommand (`metrics plot`)

Used to re-generate a **single** metric's distribution plot with custom parameters, matching the interactive use case in the R script (pick a variable, tweak bins/x-range/Tukey filter, replot).

```
phyloai pretree metrics plot --csv metrics.csv --metric num_sites
    [--bins 60] [--xmin 0] [--xmax 50000]
    [--tukey-k 1.5]
    [--title "Custom Title"] [--xlabel "Sites"] [--ylabel "Density"]
    [--color "#2E86AB"] [--fig-width 10] [--fig-height 8] [--dpi 150]
    [--font-size 12]
    [--output-dir <dir>]
```

**Parameters:**

| Parameter     | Type  | Default                         | Description                                                     |
|---------------|-------|----------------------------------|-----------------------------------------------------------------|
| `--csv`       | Path  | (required)                       | Input metrics CSV                                               |
| `--metric`    | str   | (required)                       | A single metric column to plot                                  |
| `--bins`      | int   | 50                               | Histogram bins                                                  |
| `--xmin`      | float | auto                             | X-axis minimum                                                  |
| `--xmax`      | float | auto                             | X-axis maximum                                                  |
| `--tukey-k`   | float | None (disabled)                  | Tukey's Fences multiplier; when set, filter outliers on `--metric` |
| `--title`     | str   | `Distribution of <metric_name>` | Plot title                                                      |
| `--xlabel`    | str   | `--metric` value                 | X-axis label                                                    |
| `--ylabel`    | str   | `Density`                        | Y-axis label                                                    |
| `--color`     | str   | `#2E86AB`                        | Histogram fill color (hex)                                      |
| `--fig-width` | float | 10                               | Figure width in inches                                          |
| `--fig-height`| float | 8                                | Figure height in inches                                         |
| `--dpi`       | int   | 150                              | Output resolution (DPI)                                         |
| `--font-size` | int   | 12                               | Base font size for labels and title                             |
| `--output-dir`| Path  | `runs/pretree/metrics/`          | Output directory (one PDF placed here, not in `plots/`)          |
| `--overwrite` | flag  | False                            | Overwrite existing output                                       |
| `--quiet`     | flag  | False                            | Suppress terminal output                                        |

**Tukey filter:** When `--tukey-k` is set (default 1.5 when explicitly enabled), rows where the `--metric` value falls outside `[Q1 - k×IQR, Q3 + k×IQR]` are excluded before plotting. The default k=1.5 corresponds to the standard Tukey's Fences for "outer" outliers.

### 6.3 Basic Statistics CSV

Generated alongside plot PDFs as `metrics.basic_statistics.csv`, saved to `runs/pretree/metrics/` (in the main metrics command) or the specified `--output-dir` (in the standalone `plot` subcommand):

| Column      | Description                    |
|-------------|--------------------------------|
| `metric`    | Metric name                    |
| `mean`      | Mean value                     |
| `median`    | Median value                   |
| `min`       | Minimum value                  |
| `max`       | Maximum value                  |
| `q25`       | 25th percentile                |
| `q75`       | 75th percentile                |
| `std`       | Standard deviation             |
| `n_ex_NA`   | Count of non-NA observations   |
| `n_total`   | Total observations             |

Column naming follows the project convention of explicitness (e.g., `n_ex_NA` makes clear NA values are excluded). Detailed descriptions of each column and metric are provided in the command-specific README (`docs/commands/pretree-metrics.md`).

---

## 7. Step 3: Correlation Analysis (`metrics correlate`)

### 7.1 Default Behavior

In `phyloai pretree metrics`, after computing `metrics.csv`, generate:

- **Correlation heatmap PDF:** `correlation_heatmap.pdf` under `runs/pretree/metrics/correlate/`
- **Correlation matrix CSV:** `correlation_matrix.csv` under `runs/pretree/metrics/correlate/`

### 7.2 Standalone Subcommand

```
phyloai pretree metrics correlate --csv metrics.csv
    [--metrics num_sites,entropy,...]
    [--method spearman|pearson]
    [--triangle full|lower|upper]
    [--cluster-rectangles N]
    [--cmap RdBu_r] [--annot True|False] [--fmt ".2f"]
    [--fig-width 12] [--fig-height 10] [--dpi 150]
    [--font-size 10] [--label-angle 45] [--title "Correlation Heatmap"]
    [--output-dir <dir>] [--overwrite]
```

**Parameters:**

| Parameter              | Type   | Default                          | Description                                               |
|------------------------|--------|----------------------------------|-----------------------------------------------------------|
| `--csv`                | Path   | (required)                       | Input metrics CSV                                         |
| `--metrics`            | str    | core numeric cols                | Comma-separated list of columns to correlate; `all` means every numeric column |
| `--include-freq`       | flag   | false                            | Include `freq*` columns in automatic metric selection     |
| `--include-sd`         | flag   | false                            | Include `sd_*` columns in automatic metric selection      |
| `--method`             | choice | `spearman`                       | Correlation method: `spearman` or `pearson`               |
| `--triangle`           | choice | `full`                           | `full` / `lower` / `upper` triangle display               |
| `--cluster-rectangles` | int    | None (no rectangles)             | Number of cluster rectangles to draw on full heatmaps only; ignored for `lower`/`upper` |
| `--cmap`               | str    | `RdBu_r`                         | Matplotlib colormap for the heatmap                       |
| `--annot`              | bool   | `True`                           | Show correlation values on cells                          |
| `--fmt`                | str    | `".2f"`                          | Numeric format for annotations                            |
| `--fig-width`          | float  | 12                               | Figure width in inches                                    |
| `--fig-height`         | float  | 10                               | Figure height in inches                                   |
| `--dpi`                | int    | 150                              | Output resolution                                         |
| `--font-size`          | int    | 10                               | Base font size                                            |
| `--label-angle`        | float  | 45                               | X-axis metric label rotation angle                        |
| `--title`              | str    | `Correlation Heatmap`            | Plot title                                                |
| `--output-dir`         | Path   | `runs/pretree/metrics/correlate/`| Output directory                                          |
| `--overwrite`          | flag   | False                            | Overwrite existing output                                 |
| `--quiet`              | flag   | False                            | Suppress terminal output                                  |

### 7.3 Methodology

1. **Data preparation:** By default, select core numeric biological metrics and exclude identifier columns (`loci`, `DataType`), frequency columns (`freq*`), and standard-deviation columns (`sd_*`). Users can include frequency columns with `--include-freq`, include standard-deviation columns with `--include-sd`, provide an explicit comma-separated `--metrics` list, or use `--metrics all` for every numeric column. Drop rows with any NA in selected columns (complete case analysis).
2. **Correlation:** Spearman rank correlation (no normalization needed). Matrix of size M×M.
3. **Ordering:** Compute Ward leaf order from distance `1 - |correlation|` so that strong positive and strong negative correlations group together, but do not draw dendrograms.
4. **Visualisation:** Draw a compact Matplotlib heatmap/corrplot with one colorbar, diverging colormap (`RdBu_r` or similar), fixed range `[-1, 1]`, optional annotations, and triangle masks that remove both cells and grid lines from the hidden half. Lower triangles use left/bottom labels and a right-side colorbar; upper triangles use top/right labels and a left-side colorbar. Full matrices use a right-side colorbar. Triangle modes hide rectangular spines and draw a stepped triangular border that runs outside the diagonal cells.
5. **Cluster rectangles:** Optional, drawn via `matplotlib` patch overlays on full matrices only. They are ignored for `--triangle lower` and `--triangle upper`; the CLI prints a warning when the user requests both.

### 7.4 Key Differences from R Script

| Aspect                 | R Script                                                   | Python Implementation                              |
|------------------------|------------------------------------------------------------|----------------------------------------------------|
| Variable selection     | Manual via checkbox UI                                     | Auto core metrics, opt-in `freq*`/`sd_*`, or explicit `--metrics` |
| Clustering distance    | `1 - correlation` (sign-sensitive)                          | `1 - \|correlation\|` (magnitude-based)               |
| Rectangles             | Optional, pre-set count                                    | Optional via `--cluster-rectangles` in full mode only |
| Normalization          | Not needed (Spearman)                                       | Same                                               |
| Output format          | `corrplot` + CSV download                                   | PDF + CSV                                          |

---

## 8. Output Directory and result.json

### 8.1 Default Output Layout

```
runs/pretree/metrics/
├── metrics.csv                         # per-marker metrics table
├── metrics.log                         # log file
├── result.json                         # structured result
├── metrics.basic_statistics.csv        # per-metric summary statistics
├── plots/                              # distribution plots (all metrics)
│   ├── num_taxa.pdf
│   ├── num_sites.pdf
│   ├── ...                             # one PDF per metric
└── correlate/
    ├── correlation_heatmap.pdf          # Spearman correlation heatmap
    └── correlation_matrix.csv           # full correlation matrix
```

When `metrics plot` is used independently, the single PDF is written directly to `<csv_parent>/plot_<metric>/` by default, alongside its own `result.json`. When `metrics correlate` is used independently, its outputs go to `runs/pretree/metrics/correlate/` by default.

### 8.2 result.json Schema

```json
{
  "status": "success",
  "command": "phyloai pretree metrics --msa-dir ...",
  "wall_time": 142.3,
  "tool_versions": {},
  "params": {
    "msa_dir": "/path/to/msa",
    "tree_dir": null,
    "seq_type": "AA",
    "round": 6,
    "skip_freq_statistics": false,
    "pseudo_tree_metrics": false,
    "skip_pairwise_identity": false,
    "outgroup_list": null,
    "ref_tree": null,
    "threads": 4
  },
  "key_results": {
    "n_markers": 150,
    "n_msa_metrics": 22,
    "n_tree_metrics": 15,
    "n_markers_with_trees": 120,
    "n_saturation_computed": 115,
    "n_taxon_mismatches": 3
  },
  "error": null,
  "data": {
    "metrics_csv": "runs/pretree/metrics/metrics.csv",
    "taxon_mismatches": [
      {"marker": "gene_042", "msa_only": ["taxon_X"], "tree_only": ["taxon_Y"]}
    ],
    "unpaired": {
      "msa_only": [],
      "tree_only": ["gene_099"]
    },
    "failed": [],
    "skipped": []
  }
}
```

### 8.3 Separate subcommand result.jsons

`metrics plot` and `metrics correlate` each write their own `result.json`:

**plot/result.json:**
```json
{
  "status": "success",
  "params": {"csv": "...", "metric": "num_sites", "bins": 50, "tukey_k": null},
  "key_results": {"n_plots": 1},
  "data": {"plot_pdf": "runs/pretree/metrics/num_sites.pdf", "statistics_csv": "runs/pretree/metrics/metrics.basic_statistics.csv"}
}
```

**correlate/result.json:**
```json
{
  "status": "success",
  "params": {"csv": "...", "triangle": "full", "cluster_rectangles": null, "metrics": [...]},
  "key_results": {"n_variables": 21, "n_complete_cases": 145},
  "data": {"heatmap_pdf": "...", "matrix_csv": "..."}
}
```

---

## 9. CLI Registration

In `phyloai/cli/commands/pretree.py`:

- `_PretreeGroup.list_commands` updated to include `"metrics"` after `"trim"` and before `"concat"`.
- `metrics` registered as a Click group with nested `plot` and `correlate` subcommands.

---

## 10. FastTree in `phyloai doctor`

Add `fasttree` to `core/env.py` `TOOL_REGISTRY` dict:

```python
"fasttree": {
    "required": False,
    "version_flag": "",
    "install": "Download from http://www.microbesonline.org/fasttree/",
},
```

Detection: run `FastTree` without arguments. The first line of stdout contains the version (e.g., `FastTree Version 2.1.11`). `version_flag: ""` means the tool prints its version when invoked with no arguments — the existing detection code handles this.

---

## 11. Key Design Decisions

| Decision | Rationale |
|---|---|
| Both `rcfv` and `nrcfv` | `rcfv` for backward compatibility; `nrcfv` (Fleming & Struck 2023) corrects length/taxa bias |
| `num_patterns` uses stats.py logic | Consistent with `phyloai pretree stats` output; normalised sequences from `pretree convert` make this the correct approach |
| `--dry-run` supported, `--resume` not | `--dry-run` shows plan before saturating computation; no resume because metrics is a deterministic read-only analysis — re-running produces identical output, so checkpoint overhead is unwarranted |
| `taxa_occupancy` needs two-pass | Total taxa pool not known until all files are scanned; header-only first pass is trivial I/O |
| Separate `metrics plot` and `metrics correlate` | Enables iterative exploration without recomputing metrics; plot regenerates a single metric, correlate regenerates the heatmap |
| `metrics plot` outputs one metric only | Mirrors the R script's interactive use case: pick a variable, tweak bins/x-range/Tukey, replot |
| Standalone subcommands skip dir conflict | `metrics plot` and `metrics correlate` overwrite only their own outputs when requested; they are designed for iterative replotting against an existing metrics directory |
| `metrics.csv` as canonical intermediate | Decouples computation from visualisation; filter module reads the same CSV |
| No interactive dashboard | Too complex for CLI phase; static PDF output covers the essential use case |
| Clustering by `\|\r\|` not `r` | Groups strong positive and strong negative correlations together, matching user's intent |
| Core metrics by default | Keeps the default PDF readable; `freq*` and `sd_*` metrics remain available through `--include-freq`, `--include-sd`, explicit `--metrics`, or `--metrics all` |
| DVMC pruning before computation | extract_features.py approach; outgroups affect root-to-tip distance distribution |
| Saturation implemented directly | Avoids phykit dependency; phykit algorithm is straightforward to reimplement |
| Parallel per-marker | Each marker is independent; ProcessPoolExecutor avoids GIL contention |

---

## 12. Open Questions (resolved in discussion)

| Question | Resolution |
|---|---|
| RCFV vs phykit RCV? | Output both `rcfv` (classic) and `nrcfv` (bias-corrected, Fleming & Struck 2023); phykit RCV is a different metric entirely (count-based, not frequency-based) |
| `num_patterns` character handling? | stats.py approach (normalize → gap) |
| Interactive dashboard? | No, CLI subcommands instead |
| Correlation as separate module? | Yes, `metrics correlate` |
| MSA–tree taxon check? | Verify equality; warn on mismatch |
| DVMC formula consistency? | Identical without outgroup; prune first with outgroup |
