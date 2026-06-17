# phyloai pretree filter

## Purpose

`phyloai pretree filter` is the quality-control gate between tree inference and supermatrix concatenation. It provides four complementary workflows:

| Subcommand | Scope | Action |
|------------|-------|--------|
| `taper` | Site-level | Mask erroneous positions; retain the locus |
| `treeshrink` | Taxon-level | Prune outlier long-branch taxa from gene trees (and optionally MSAs) |
| `metrics` | Locus-level | Drop or retain whole loci by numeric/string rules on a metrics table |
| `cluster` | Population-level | Group loci by metric profiles (PCA/UMAP + hierarchical clustering); optionally drop outlier clusters |

**Position in the pipeline:**

```
phyloai pretree align    →  aligned MSAs
phyloai pretree trim     →  trimmed MSAs
phyloai pretree filter   →  quality-controlled MSAs & trees  ← YOU ARE HERE
phyloai pretree concat   →  supermatrix
```

Filter does **not** compute metrics (use `phyloai pretree metrics`), infer trees, or concatenate. It reads `pretree metrics` output and writes structured decisions and filtered files.

All subcommands write `result.json` and `filter.log` to their output directory. Terminal output uses Rich tables; suppress with `--quiet`.

---

## Architecture

```
phyloai/core/file_matching.py       ← shared logical locus-name parsing, MSA/tree scanning, pairing
phyloai/pretree/filter.py           ← core library: run_taper, run_treeshrink, run_metrics_filter, run_cluster_filter
phyloai/cli/commands/pretree.py     ← CLI: filter Click group with four ordered subcommands
```

**Shared infrastructure:**
- `core/checkpoint.py` + `pretree/checkpoint_helpers.py` — used by TAPER for resume (follows the same pattern as `pretree align` and `pretree trim`)
- `core/env.py` — tool detection for Julia and `run_treeshrink.py`
- `core/file_matching.py` — suffix-agnostic file scanning and locus-name matching (shared with `pretree metrics`)

**External tools:**
- TAPER: bundled `phyloai/bundled/TAPER-1.0.0/correction_multi.jl`, executed by Julia (must be on PATH or specified via `--julia-path`)
- TreeShrink: user-installed `run_treeshrink.py` (must be on PATH or specified via `--treeshrink-path`)

**Dependencies:** `scikit-learn` for PCA and clustering; `umap-learn` optional for UMAP reduction.

### Design decisions

Key architectural choices and their rationale:

| Decision | Why |
|----------|-----|
| Separate subcommands, not one combined command | Avoids many invalid option combinations; each subcommand has distinct required inputs, tool requirements, and output semantics |
| TAPER has no `--tree-dir` | Masked MSAs should be used for *new* tree inference; copying old trees would risk using stale topologies |
| TAPER runs per-locus in parallel | Matches existing `align`/`trim` patterns for process pool and checkpoint/resume |
| TreeShrink runs once across the dataset | TreeShrink's statistical model can use information from multiple trees jointly; per-locus invocation would break this |
| TreeShrink has no `--resume` or `--threads` in v1 | Keeps the first version simple and statistically correct; the external tool runs as one batch |
| `metrics` supports AND-only logic | Keeps parsing transparent and safe; OR logic can be achieved by running multiple filter passes |
| Cluster defaults to PCA | Stable, deterministic, dependency-light (`scikit-learn` only). UMAP is available for nonlinear exploration |
| Cluster auto-drop is opt-in | Cluster interpretation should remain user-guided; automatic removal is potentially risky |
| Cluster count uses lightweight validation voting (silhouette, Calinski-Harabasz, Davies-Bouldin) | Transparent, dependency-light alternative to R's NbClust; sufficient for exploratory outlier screening |
| No `pandas` dependency | Standard-library `csv` plus `numpy` |
| No `plotly` dependency | Existing `matplotlib` for all diagnostic plots |

---

## Common options

All four subcommands share these options:

| Option | Default | Purpose |
|--------|---------|---------|
| `--output-dir` / `-o` | `runs/pretree/filter/<subcommand>` | All outputs go here |
| `--table-format` | `csv` | Delimiter and suffix for auxiliary tables (`retained_loci`, `dropped_loci`, `filter_decisions`, etc.). Does not affect `result.json`. |
| `--overwrite` | off | Delete and recreate `--output-dir` if it exists |
| `--dry-run` | off | Validate inputs, show planned actions, write no files (not even `result.json`) |
| `--quiet` / `-q` | off | Suppress all terminal output except errors |

### File matching policy

All subcommands that accept `--msa-dir` or `--tree-dir` use suffix-agnostic logical locus-name matching rather than hard-coded extension whitelists:

| File | Logical locus name |
|------|--------------------|
| `gene1.fa` | `gene1` |
| `gene2.v1.ALI` | `gene2.v1` |
| `gene3.treefile` | `gene3` |
| `gene4.fa.treefile` | `gene4.fa`, then `gene4` |

For MSA directories, every regular non-empty file is scanned; file format is validated when parsed, not by extension. For tree directories, the same rule applies but with one- and two-suffix reduction candidates to handle `.treefile`, `.fa.treefile`, etc. Ambiguous matches (where a filename's candidates are already occupied) cause an error with details.

`phyloai pretree metrics` uses the same helpers, so `metrics` and `filter` behave identically for non-standard file naming.

---

## `filter taper` — TAPER Error-Site Masking

### When to use

After alignment and trimming, individual sites within otherwise-good loci may still contain sequencing errors or alignment artifacts. TAPER identifies stretches of amino acids (or nucleotides) that are unexpectedly divergent relative to the rest of the alignment and masks them to `X` (or `N`).

Use `taper` for **site-level** quality control. If you want to remove entire loci, use `filter metrics`. If you want to remove individual taxa, use `filter treeshrink`.

### CLI

```bash
phyloai pretree filter taper \
  --msa-dir <aa_or_nt_msa_dir> \
  [--nt-dir <codon_aligned_nt_msa_dir>] \
  [--seq-type AA|NT|auto] \
  [--cutoff 3] \
  [--taper-path <correction_multi.jl>] \
  [--julia-path <julia>] \
  [--tool-args "..."] \
  [--show-masked-sites] \
  [--output-dir runs/pretree/filter/taper] \
  [--threads 4] [--resume] [--dry-run] [--overwrite]
```

### Operating modes

| Mode | Input | Command | Output |
|------|-------|---------|--------|
| **AA-only** | AA MSA files in `--msa-dir` | `julia correction_multi.jl -c <cutoff> INPUT > OUTPUT` | Masked AA → `seqs/` |
| **NT-only** | NT MSA files in `--msa-dir`, `--seq-type NT` | `julia correction_multi.jl -m N -a N -c <cutoff> INPUT > OUTPUT` | Masked NT → `seqs/` |
| **AA+CDS** | AA MSA in `--msa-dir` + codon-aligned NT MSA in `--nt-dir` | Run on AA; project masks to CDS | Masked AA → `seqs/faa/`, projected CDS → `seqs/fna/` |

AA+CDS mode projection rules:
- Original `X` in the input → left unchanged (not counted as a TAPER mask)
- Standard AA → `X` by TAPER → corresponding codon replaced with `NNN`
- Gap (`-`) → `X` → warning; no CDS change (defensive check, not expected in normal output)

### Parameters

| Parameter | Default | Notes |
|-----------|---------|-------|
| `--msa-dir` | required | Directory of aligned MSA files (any suffix). All regular non-empty files scanned. |
| `--nt-dir` | — | AA+CDS mode. Codon-aligned NT MSA directory (one file per AA locus, NT length == 3 × AA length). Raw CDS not accepted. |
| `--seq-type` | `auto` | `AA`, `NT`, or `auto`. `auto` inspects the first file: presence of `EFILPQWYZ` → AA, otherwise NT. Override when auto-detection is unreliable. |
| `--cutoff` | 3 | TAPER `-c` error-correction cutoff (≥1). **Lower = more aggressive.** 1-10 typical; 1-2 for noisy data, 3 is TAPER's default, 5+ for conservative. |
| `--taper-path` | — | Explicit path to `correction_multi.jl`. Uses bundled copy by default. |
| `--julia-path` | — | Explicit Julia executable. Resolved via PATH; check with `phyloai doctor`. |
| `--tool-args` | — | Additional TAPER flags passed verbatim. Managed flags (`-m`, `-a`, `-c`, `-l`, input path, output redirection) cause an error. |
| `--threads` / `-t` | 4 | Worker processes (one locus per worker). Follows the same `ProcessPoolExecutor` + checkpoint pattern as `pretree align` and `pretree trim`. |
| `--show-masked-sites` | off | Include `masked_taxa_detail` column in `filter_decisions.csv` (`taxonA:3; taxonB:5`). Default off for compact output on large datasets. |
| `--resume` | off | Resume from `checkpoint.json`. Parameters must match; completed loci skipped only if their output files exist and pass validation. |
| `--overwrite` | off | Mutually exclusive with `--resume`. Delete and recreate `--output-dir`. |
| `--dry-run` | off | Show detected mode, paired loci count, output layout, tool command template, resume eligibility. No files written. |
| `--table-format` | `csv` | Format for `retained_loci`, `dropped_loci`, `filter_decisions`. |
| `--output-dir` / `-o` | `runs/pretree/filter/taper` | Output directory. |

### Outputs

```
runs/pretree/filter/taper/
├── seqs/                              (or seqs/faa/ + seqs/fna/ for AA+CDS)
├── retained_loci.csv|tsv
├── dropped_loci.csv|tsv               (locus, reason)
├── filter_decisions.csv|tsv           (locus, status, new_masked_sites, masked_taxa_count,
│                                       masked_taxa_detail when --show-masked-sites)
├── checkpoint.json                    (internal; for --resume)
├── filter.log
└── result.json
```

`filter_decisions.csv` columns:
- `locus`, `status` (success/failed)
- `new_masked_sites` — total AA sites newly masked to `X` (original input `X` not counted)
- `masked_taxa_count` — number of taxa with ≥1 masked site
- `masked_taxa_detail` — semicolon-separated `taxon:count` per locus (only with `--show-masked-sites`)

Terminal output: Filter Results table (input/retained/dropped/masked loci/masked taxa/masked sites) + Retained MSA Statistics table (MSA count, total/mean/min/max length, mean taxa). Julia version auto-detected and recorded in `result.json` and `filter.log`.

### Examples

```bash
# Default AA masking (most common starting point)
phyloai pretree filter taper --msa-dir ./trimmed

# Aggressive masking for noisy data
phyloai pretree filter taper --msa-dir ./trimmed --cutoff 1 --threads 8

# Conservative masking (only clearest errors)
phyloai pretree filter taper --msa-dir ./trimmed --cutoff 5

# NT-only mode
phyloai pretree filter taper --msa-dir ./trimmed_nt --seq-type NT

# AA+CDS: mask on AA, project to codon-aligned NT
phyloai pretree filter taper --msa-dir ./trimmed_aa --nt-dir ./trimmed_fna

# Resume after interruption
phyloai pretree filter taper --msa-dir ./trimmed --resume

# Inspect which taxa got masked
phyloai pretree filter taper --msa-dir ./trimmed --show-masked-sites
```

### Warnings and errors

| Condition | Behaviour |
|-----------|-----------|
| `--nt-dir` + `--seq-type NT` | Exit 1 — AA+CDS needs AA input |
| `--threads` < 1 | Exit 1 |
| `--resume` + `--overwrite` | Exit 1 — mutually exclusive |
| Julia not found | Exit 3 — install Julia or set `--julia-path` |
| Non-empty output directory without `--overwrite` or `--resume` | Exit 1 |
| No valid MSA files | Exit 1 |
| TAPER exits non-zero | Locus skipped; reason in `dropped_loci.csv` |
| TAPER output missing or fails FASTA validation | Locus skipped |
| All loci fail | Exit 2 |

---

## `filter treeshrink` — TreeShrink Taxon Pruning

### When to use

After gene tree inference, some taxa may appear on unusually long branches — typically artifacts of misassembly, paralogy, or contamination. TreeShrink uses a statistical test to detect outlier branch lengths across all trees simultaneously and prunes the offending taxa.

This is **taxon-level** filtering: taxa are removed from specific gene trees. The locus is retained. If you want to remove entire poorly-performing loci, use `filter metrics`.

**Design rationale:** TreeShrink is invoked once across the entire dataset (not per locus) because its statistical model can pool information from multiple trees. Per-locus invocation would change the statistical model. For the same reason, `--resume` and `--threads` are not provided — the external tool runs as one batch.

### CLI

```bash
phyloai pretree filter treeshrink \
  --tree-dir <gene_tree_dir> \
  [--msa-dir <msa_dir>] \
  [--threshold 0.05] \
  [--treeshrink-mode auto|per-gene|all-genes|per-species] \
  [--treeshrink-path <run_treeshrink.py>] \
  [--tool-args "..."] \
  [--output-dir runs/pretree/filter/treeshrink] \
  [--keep-work-dir] [--dry-run] [--overwrite]
```

### How it works

1. Scan `--tree-dir` for gene tree files; optionally scan `--msa-dir` for matching MSAs (paired by logical locus name)
2. Create a per-gene work layout in a temporary directory (`gene/input.tree`, optional `gene/input.fasta`)
3. Execute `run_treeshrink.py -i <work_dir> -t input.tree [-a input.fasta] -q <threshold> [-m <mode>]`
4. Collect `output.tree` (and `output.fasta` if MSA mode) for each gene
5. Compare input and output trees to track which taxa were removed per locus
6. Clean up work directory unless `--keep-work-dir` is set

### Parameters

| Parameter | Default | Notes |
|-----------|---------|-------|
| `--tree-dir` | required | Directory of gene tree files (any suffix). |
| `--msa-dir` | — | Optional MSA directory. When provided, MSAs are paired by locus name and shrunk alongside trees. |
| `--threshold` | 0.05 | TreeShrink `-q` false-positive threshold. **Smaller = more taxa removed.** 0.05 works for most datasets; 0.01 for aggressive, 0.1 for conservative. |
| `--treeshrink-mode` | `auto` | `auto` omits `-m` (TreeShrink default). `per-gene`: independent per gene. `all-genes` / `per-species`: pools cross-gene information. |
| `--treeshrink-path` | — | Explicit path to `run_treeshrink.py`. Resolved via PATH by default. |
| `--tool-args` | — | Additional TreeShrink flags. Managed flags (`-i`, `-t`, `-a`, `-q`, `-m`, `-o`, `-O`) cause an error. |
| `--keep-work-dir` | off | Retain per-gene work directory under `--output-dir/work/` for debugging failed loci. |
| `--output-dir` / `-o` | `runs/pretree/filter/treeshrink` | Output directory. |
| `--table-format` | `csv` | Format for auxiliary tables. |
| `--overwrite` | off | Delete and recreate output directory. |
| `--dry-run` | off | Print resolved command and locus count without executing. |
| `--quiet` / `-q` | off | Suppress terminal output. |

### Outputs

```
runs/pretree/filter/treeshrink/
├── trees/                              (shrunk gene trees)
├── seqs/                               (only when --msa-dir provided)
├── retained_loci.csv|tsv
├── modified_loci.csv|tsv               (loci where ≥1 taxon was pruned)
├── dropped_loci.csv|tsv                (loci with missing/invalid output)
├── removed_taxa.csv|tsv                (locus, taxon per row)
├── filter_decisions.csv|tsv            (locus, status, removed_count)
├── work/                               (only with --keep-work-dir)
├── filter.log
└── result.json
```

Terminal output: Filter Results table (input/retained/modified/dropped/taxa removed) + Retained MSA Statistics table (when `--msa-dir` provided) + contextual tip.

### Examples

```bash
# Basic taxon pruning
phyloai pretree filter treeshrink --tree-dir ./genetrees

# Trees + matching MSAs with conservative threshold
phyloai pretree filter treeshrink \
  --tree-dir ./genetrees \
  --msa-dir ./trimmed \
  --threshold 0.1

# Per-species mode (cross-gene information pooling)
phyloai pretree filter treeshrink \
  --tree-dir ./genetrees \
  --treeshrink-mode per-species

# Debug failed loci
phyloai pretree filter treeshrink --tree-dir ./genetrees --keep-work-dir
```

### Warnings and errors

| Condition | Behaviour |
|-----------|-----------|
| `run_treeshrink.py` not found | Exit 3 |
| No valid tree files | Exit 1 |
| Ambiguous locus matching (tree ↔ MSA) | Exit 1 with details |
| TreeShrink exits non-zero | All loci marked failed; Exit 2 if none succeed |
| Non-empty output directory without `--overwrite` | Exit 1 |

---

## `filter metrics` — Metric Rule Filtering

### When to use

After computing per-locus quality metrics with `phyloai pretree metrics`, use `filter metrics` to drop entire loci that fail explicit thresholds. This is **locus-level** filtering: the whole gene is kept or discarded.

`filter metrics` is deliberately separate from `pretree metrics` computation — you can explore multiple threshold combinations without re-computing metrics each time.

### CLI

```bash
phyloai pretree filter metrics \
  --table <metrics.csv|metrics.tsv> \
  --keep "dvmc>=0,dvmc<=0.3,average_BS>=0.8" \
  [--input-format auto|csv|tsv] \
  [--loci-column loci] \
  [--msa-dir <msa_dir>] [--tree-dir <tree_dir>] [--copy] \
  [--output-dir runs/pretree/filter/metrics] \
  [--table-format csv|tsv] [--dry-run] [--overwrite]
```

### Rule syntax

Comma-separated conditions, AND logic (all must pass):

```
column operator value
```

| Operator | Value type | Example |
|----------|-----------|---------|
| `>=`, `>`, `<=`, `<` | numeric only | `average_BS>=0.8` |
| `==`, `!=` | numeric or string | `DataType==AA` |

Using `>=`/`>`/`<=`/`<` on a string column exits with an error.

```bash
# Good: numeric thresholds
--keep "dvmc>=0,dvmc<=0.3,average_BS>=0.8"

# Good: mixed numeric + string
--keep "DataType==AA,num_sites>=300"

# Good: simple single condition
--keep "num_sites>=1000"

# Error: inequality on string column
--keep "DataType>=AA"
```

OR logic is not supported in v1. To achieve OR (e.g., "drop if dvmc > 0.3 OR average_BS < 0.5"), run the command twice with different `--output-dir` values.

### Choosing thresholds

Thresholds are dataset-dependent. As starting points for a typical phylogenomic dataset:

| Metric | Conservative | Moderate | Lenient |
|--------|-------------|----------|---------|
| `average_BS` | ≥ 0.8 | ≥ 0.5 | ≥ 0.3 |
| `dvmc` | ≤ 0.1 | ≤ 0.3 | ≤ 0.5 |
| `num_sites` | ≥ 200 | ≥ 100 | ≥ 50 |
| `gap_proportion` | ≤ 0.2 | ≤ 0.5 | ≤ 0.8 |

Start moderate, check `condition_failure_counts` in `result.json`, and relax or remove the most restrictive conditions.

### Parameters

| Parameter | Default | Notes |
|-----------|---------|-------|
| `--table` | required | Metrics CSV/TSV (typically from `phyloai pretree metrics`). Delimiter auto-detected. |
| `--keep` | required | Comma-separated AND conditions. |
| `--input-format` | `auto` | `csv`, `tsv`, or `auto`. Override when auto-detection fails. |
| `--loci-column` | `loci` | Column name for locus identifier. |
| `--msa-dir` | — | MSA directory for computing retained-MSA statistics (shown in terminal, stored in `result.json`). |
| `--tree-dir` | — | Tree directory for `--copy` mode. |
| `--copy` | off | Copy retained MSAs/trees to `--output-dir/seqs/` and `--output-dir/trees/`. Requires `--msa-dir` or `--tree-dir`. |
| `--output-dir` / `-o` | `runs/pretree/filter/metrics` | Output directory. |
| `--table-format` | `csv` | Format for auxiliary tables. |
| `--overwrite` | off | Delete and recreate output directory. |
| `--dry-run` | off | Parse rules, report pass/fail counts, write no files. |
| `--quiet` / `-q` | off | Suppress terminal output. |

### Outputs

```
runs/pretree/filter/metrics/
├── retained_loci.csv|tsv
├── dropped_loci.csv|tsv
├── filter_decisions.csv|tsv
├── seqs/                              (only with --copy --msa-dir)
├── trees/                             (only with --copy --tree-dir)
├── filter.log
└── result.json
```

`result.json.key_results.condition_failure_counts` maps each condition to how many loci it rejected — the key diagnostic for understanding which threshold is most restrictive.

### Examples

```bash
# Basic quality filter (moderate thresholds)
phyloai pretree filter metrics \
  --table ./metrics/metrics.csv \
  --keep "dvmc<=0.3,average_BS>=0.8"

# Mixed numeric + string with file copy
phyloai pretree filter metrics \
  --table ./metrics.csv \
  --keep "DataType==AA,num_sites>=300" \
  --copy --msa-dir ./trimmed

# Dry-run first to explore thresholds
phyloai pretree filter metrics \
  --table ./metrics.csv \
  --keep "num_sites>=500,average_BS>=0.7" \
  --dry-run

# Custom locus identifier column
phyloai pretree filter metrics \
  --table ./table.tsv \
  --keep "average_BS>=0.9" \
  --loci-column gene_id

# Compare two filtering strategies (different output dirs)
phyloai pretree filter metrics \
  --table ./metrics.csv \
  --keep "average_BS>=0.8" \
  -o ./runs/strategy_conservative

phyloai pretree filter metrics \
  --table ./metrics.csv \
  --keep "average_BS>=0.5" \
  -o ./runs/strategy_lenient
```

### Warnings and errors

| Condition | Behaviour |
|-----------|-----------|
| `--table` missing or empty | Exit 1 |
| Malformed `--keep` syntax | Exit 1 with parse detail |
| Unknown column in `--keep` | Exit 1 |
| `>=`/`>`/`<=`/`<` on string column | Exit 1 |
| `--copy` without `--msa-dir` or `--tree-dir` | Exit 1 |
| No loci pass all conditions | 0 retained (not an error) |
| Non-empty output directory without `--overwrite` | Exit 1 |

---

## `filter cluster` — Cluster-Based Exploration

### When to use

`filter cluster` is the **exploratory** arm of `filter`. It groups loci by their metric profiles using dimensionality reduction (PCA or UMAP) plus hierarchical clustering, then produces diagnostic visualizations. Use it to:

- Identify sub-populations of loci with systematically different properties
- Spot isolated outliers
- Understand whether single-metric filtering is missing correlated problems
- Decide on filtering strategy before applying `filter metrics`

By default, `cluster` writes diagnostics only — no loci are removed. Use `--drop-outlier-clusters auto` to optionally remove the worst-performing clusters.

### CLI

```bash
phyloai pretree filter cluster \
  --table <metrics.csv|metrics.tsv> \
  [--input-format auto|csv|tsv] \
  [--metrics all|col1,col2,...] \
  [--exclude-regex REGEX] [--exclude-regex REGEX] \
  [--reduction pca|umap] \
  [--n-clusters N] [--max-clusters N] \
  [--cluster-linkage ward|average|complete|single] \
  [--cluster-distance euclidean|cosine|manhattan] \
  [--drop-outlier-clusters none|auto] \
  [--outlier-metric average_BS] [--outlier-direction low|high] \
  [--max-drop-fraction 0.2] \
  [--plot-metrics-per-page auto|N] [--plot-label-angle 45] \
  [--umap-n-neighbors 15] [--umap-min-dist 0.001] \
  [--umap-replicates 1] [--umap-random-state 0] \
  [--msa-dir <msa_dir>] [--tree-dir <tree_dir>] [--copy] \
  [--output-dir runs/pretree/filter/cluster] \
  [--dry-run] [--overwrite]
```

### How it works

1. **Feature selection**: All numeric columns from the input table, excluding the locus ID column, `DataType`, constant columns, and columns matching `--exclude-regex`. The `features_used.csv` output is your audit trail.

2. **Scaling**: Features z-score standardized (mean 0, variance 1).

3. **Dimensionality reduction**: PCA (default, 3 components via `sklearn.decomposition.PCA`) or UMAP (2D embedding, optional `umap-learn`). PCA is preferred for determinism and speed; UMAP for detecting nonlinear manifold structure.

4. **Clustering**: Agglomerative hierarchical clustering (`sklearn.cluster.AgglomerativeClustering`). Default: Ward linkage + Euclidean distance (minimizes within-cluster variance).

5. **Cluster count selection** (when `--n-clusters` not set): Evaluate `k=2..max_clusters`. Each of three internal validation metrics votes for its best `k`:
   - Silhouette score (higher = better)
   - Calinski-Harabasz index (higher = better)
   - Davies-Bouldin index (lower = better)
   
   Ties broken by higher silhouette, then smaller `k`. Saved in `cluster_selection.csv`.

6. **Diagnostics**: 2D/3D scatter plots, per-cluster metric heatmap, per-cluster metric means table, per-metric boxplots by cluster.

7. **Outlier dropping** (only when `--drop-outlier-clusters auto`):
   - Rank clusters by mean `--outlier-metric`
   - `--outlier-direction low`: smaller values worse (e.g., `average_BS`)
   - `--outlier-direction high`: larger values worse (e.g., `dvmc`)
   - Drop worst clusters until cumulative fraction exceeds `--max-drop-fraction`

### Choosing the right method

| Aspect | PCA | UMAP |
|--------|-----|------|
| Deterministic | Yes | No (use `--umap-replicates` for stability) |
| Dependencies | scikit-learn (always installed) | `pip install umap-learn` |
| Best for | Linear structure, speed, reproducibility | Non-linear manifolds, complex relationships |
| Output | `PC1`, `PC2`, `PC3` | `UMAP1`, `UMAP2` |

Start with PCA (default). Try UMAP if PCA clusters look uninformative or if you suspect nonlinear metric relationships.

### Parameters

| Parameter | Default | Notes |
|-----------|---------|-------|
| `--table` | required | Metrics CSV/TSV. |
| `--input-format` | `auto` | `csv`, `tsv`, or `auto`. |
| `--metrics` | all | Comma-separated metric columns. Default uses all numeric except exclusions. |
| `--exclude-regex` | — | Repeatable. `--exclude-regex '^freq' --exclude-regex '^sd_'`. Combined regex also valid: `--exclude-regex '^(freq\|sd_)'`. |
| `--reduction` | `pca` | `pca` or `umap`. |
| `--n-clusters` | auto | Fixed cluster count. Auto-selected when omitted. |
| `--max-clusters` | auto | Upper bound for auto-selection. Default: `min(30, max(6, ceil(sqrt(n_loci)/3)))`. |
| `--cluster-linkage` | `ward` | `ward` (minimizes within-cluster variance), `average`, `complete`, or `single`. |
| `--cluster-distance` | `euclidean` | `euclidean`, `cosine`, or `manhattan`. Ward requires Euclidean. |
| `--drop-outlier-clusters` | `none` | `none`: diagnostics only. `auto`: remove worst clusters. |
| `--outlier-metric` | `average_BS` | Metric column for ranking clusters. |
| `--outlier-direction` | `low` | `low`: smaller values worse (BS, support). `high`: larger values worse (dvmc, gap). |
| `--max-drop-fraction` | 0.2 | Safety limit on fraction of loci dropped (0.0–1.0). |
| `--msa-dir` / `--tree-dir` | — | Input directories for `--copy` mode. |
| `--copy` | off | Copy retained MSAs/trees when outlier dropping is active. No-op if no loci dropped. |
| `--plot-metrics-per-page` | `auto` | Boxplots per PDF page. `auto`: ≤6 clusters → 12/page, ≤12 → 6, ≤20 → 4, >20 → 2. |
| `--plot-label-angle` | 45.0 | X-axis label rotation in diagnostic plots. |
| `--umap-n-neighbors` | 15 | UMAP local/global balance. |
| `--umap-min-dist` | 0.001 | UMAP point packing. |
| `--umap-replicates` | 1 | UMAP runs with different seeds; best selected by cluster-validation rank-sum scoring. |
| `--umap-random-state` | 0 | Base random seed. Each replicate uses `base + index`. |
| `--output-dir` / `-o` | `runs/pretree/filter/cluster` | Output directory. |
| `--table-format` | `csv` | Format for all output tables. |
| `--overwrite` | off | Delete and recreate output directory. |
| `--dry-run` | off | Show selected features, reduction, cluster count range, drop plan. No files. |
| `--quiet` / `-q` | off | Suppress terminal output. |

### Outputs

Core (always written):
```
runs/pretree/filter/cluster/
├── features_used.csv|tsv              (column, included, reason)
├── reduction.csv|tsv                  (coordinates per locus: PC1/PC2/PC3 or UMAP1/UMAP2)
├── cluster_selection.csv|tsv          (k, silhouette, calinski_harabasz, davies_bouldin)
├── clusters.csv|tsv                   (locus, cluster, distance to centroid)
├── cluster_summary.csv|tsv            (per-cluster size, metric means)
├── cluster_metric_means.csv|tsv       (per-cluster mean for every numeric metric)
├── cluster_metric_heatmap.pdf         (z-score heatmap: metrics × clusters)
├── cluster_2d.pdf                     (first 2 reduced dimensions, colored by cluster)
├── cluster_3d.pdf                     (3D scatter, PCA only)
├── cluster_metric_boxplots_*.pdf      (per-metric distributions by cluster)
├── cluster_loci/cluster_*.csv|tsv     (loci in each cluster)
├── filter.log
└── result.json
```

With `--drop-outlier-clusters auto`, additionally:
```
├── retained_loci.csv|tsv              (after dropping)
├── dropped_loci.csv|tsv
├── filter_decisions.csv|tsv
├── outlier_comparison.csv|tsv        (normal vs outlier: mean, median, sd, count per metric)
├── outlier_comparison_boxplots_*.pdf (per-metric distributions by outlier status)
├── seqs/                             (only with --copy --msa-dir)
└── trees/                            (only with --copy --tree-dir)
```

### Examples

```bash
# Exploratory: see how loci cluster by metric profiles
phyloai pretree filter cluster --table ./metrics/metrics.csv

# UMAP with fixed cluster count
phyloai pretree filter cluster \
  --table ./metrics.csv \
  --reduction umap --n-clusters 5

# Drop outlier clusters by average_BS, copy surviving files
phyloai pretree filter cluster \
  --table ./metrics.csv \
  --drop-outlier-clusters auto \
  --outlier-metric average_BS \
  --max-drop-fraction 0.15 \
  --copy --msa-dir ./trimmed

# Exclude frequency and standard deviation columns
phyloai pretree filter cluster \
  --table ./metrics.csv \
  --exclude-regex '^freq' --exclude-regex '^sd_'

# Specific metric subset
phyloai pretree filter cluster \
  --table ./metrics.csv \
  --metrics "average_BS,dvmc,gc_content,num_sites"
```

### Notes

- Without `--drop-outlier-clusters auto`, **no loci are removed**. All outputs are diagnostic.
- `features_used.csv` is your audit trail: shows every column, whether included, and why excluded.
- UMAP requires `pip install umap-learn`. Missing dependency exits with a clear install hint.
- `--resume` and `--threads` are not supported — clustering runs once in-memory.
- Cluster count selection uses lightweight validation voting, not R's `NbClust`. This is intentionally simpler and dependency-light, suitable for exploratory screening.

---

## Common workflows

### Full pipeline: taper → metrics → filter → concat

```bash
# 1. Mask error sites
phyloai pretree filter taper --msa-dir ./trimmed -o ./runs/taper

# 2. Compute metrics on masked MSAs and gene trees
phyloai pretree metrics \
  --msa-dir ./runs/taper/seqs \
  --tree-dir ./genetrees \
  -o ./runs/metrics

# 3. Filter by quality thresholds (copy surviving MSAs)
phyloai pretree filter metrics \
  --table ./runs/metrics/metrics.csv \
  --keep "average_BS>=0.7,dvmc<=0.3,num_sites>=100" \
  --copy --msa-dir ./runs/taper/seqs \
  -o ./runs/filtered

# 4. Concatenate into supermatrix
phyloai pretree concat \
  --msa-dir ./runs/filtered/seqs \
  -o ./runs/concat
```

### Add TreeShrink before metrics

```bash
# TreeShrink first to clean trees, then metrics on shrunk data
phyloai pretree filter treeshrink \
  --tree-dir ./genetrees \
  --msa-dir ./trimmed \
  --threshold 0.1 \
  -o ./runs/treeshrink

phyloai pretree metrics \
  --msa-dir ./runs/treeshrink/seqs \
  --tree-dir ./runs/treeshrink/trees \
  -o ./runs/metrics_shrunk
```

### Iterative threshold exploration

```bash
# Dry-run to explore without writing files
phyloai pretree filter metrics --table ./metrics.csv \
  --keep "average_BS>=0.8" --dry-run

phyloai pretree filter metrics --table ./metrics.csv \
  --keep "average_BS>=0.7" --dry-run

phyloai pretree filter metrics --table ./metrics.csv \
  --keep "average_BS>=0.5" --dry-run

# Write results for the chosen threshold
phyloai pretree filter metrics --table ./metrics.csv \
  --keep "average_BS>=0.7" \
  --copy --msa-dir ./taper/seqs \
  -o ./runs/filtered
```

### Cluster-guided filtering

```bash
# 1. Explore clustering
phyloai pretree filter cluster --table ./metrics.csv -o ./runs/cluster_explore

# 2. Inspect cluster_metric_means.csv, cluster_metric_heatmap.pdf

# 3. If a clear outlier cluster exists, drop it with copy
phyloai pretree filter cluster \
  --table ./metrics.csv \
  --drop-outlier-clusters auto \
  --outlier-metric average_BS \
  --max-drop-fraction 0.2 \
  --copy --msa-dir ./taper/seqs \
  -o ./runs/cluster_filtered
```

---

## result.json schema

All subcommands follow the same structured result format (matching main design §9.4):

```json
{
  "status": "success | error",
  "command": "phyloai pretree filter ...",
  "wall_time": 1.23,
  "tool_versions": {"julia": "julia version 1.12.6", ...},
  "params": {"msa_dir": "...", "cutoff": 3, ...},
  "key_results": {
    "n_input": 1066,
    "n_retained": 1064,
    "n_dropped": 2,
    "masked_loci": 42,
    "total_masked_taxa": 156,
    "total_masked_aa_sites": 1203
  },
  "error": null,
  "data": {
    "retained_msa_stats": {"n_msa": 1064, "total_length": 312345, ...},
    ...
  }
}
```

`filter.log` records the resolved command(s), tool versions, wall time, exit codes, and per-locus outcomes.
