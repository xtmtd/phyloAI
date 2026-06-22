# phyloai pretree filter

## Purpose

`phyloai pretree filter` is the quality-control gate between tree inference and supermatrix concatenation. It provides five complementary filtering workflows:

- `taper` — mask erroneous sites within individual loci (site-level)
- `treeshrink` — prune outlier long-branch taxa from gene trees (taxon-level)
- `metrics` — drop or retain whole loci by numeric/string rules on a metrics table (locus-level)
- `cluster` — group loci by metric profiles; optionally drop outlier clusters (population-level)

Filter does **not** compute marker metrics (use `phyloai pretree metrics`), infer gene trees, or concatenate retained MSAs (use `phyloai pretree concat`). It reads `pretree metrics` output tables where appropriate and writes structured decision files and optionally filtered MSA/tree directories.

The module keeps metric computation separate from filtering decisions so you can explore multiple threshold combinations without re-computing metrics.

Position in the pipeline:

```
phyloai pretree align    →  aligned MSAs
phyloai pretree trim     →  trimmed MSAs
phyloai pretree filter   →  quality-controlled MSAs & trees  ← YOU ARE HERE
phyloai pretree concat   →  supermatrix
```

### Recommended workflow

For typical phylogenomic MSAs, the five subcommands are designed to be applied in sequence:

1. **`taper`** — mask potential site-level errors in the trimmed MSAs. This produces cleaner alignments without discarding loci or taxa.

2. **Build gene trees** from the masked MSAs (using an external tree-inference tool). These trees reflect the corrected sequences.

3. **`treeshrink`** — feed the gene trees (and optionally the masked MSAs) into TreeShrink to identify and prune outlier long-branch taxa. The result is a set of shrunk trees and optionally shrunk MSAs with problematic taxa removed.

4. **(Optional) Re-infer gene trees** on the shrunk MSAs — TreeShrink ensures long-branch taxa are removed, but the tree topology may improve further with the pruned alignment.

5. **`metrics`**, **`symtest`** and/or **`cluster`** — only after site masking and taxon pruning should you evaluate per-locus quality. These subcommands are optional but complementary: `metrics` computes marker attributes for rule-based filtering, `symtest` tests phylogenetic symmetry assumptions (stationarity/homogeneity) via IQ-TREE3 and filters loci by p-value, and `cluster` groups loci by metric profiles for exploration and outlier removal. Use one or more of these on the cleaned dataset as needed.

The subcommands can also be used independently. For example, if you already have gene trees and only want to prune taxa, start from step 3. If you only need to apply metric thresholds, jump directly to `filter metrics`.

All subcommands write `result.json` to their output directory. Terminal output uses Rich tables; suppress with `--quiet`.

Any PhyloAI-authored FASTA-family outputs written by these subcommands wrap sequence lines at 60 characters.

### Shared options

| Option | Default | Purpose |
|--------|---------|---------|
| `--output-dir` / `-o` | `runs/pretree/filter/<subcommand>` | Output directory |
| `--table-format` | `csv` | Delimiter and suffix for auxiliary tables; does not affect `result.json` |
| `--overwrite` | off | Delete and recreate `--output-dir` if it exists |
| `--dry-run` | off | Validate inputs and show planned actions; no files written (no `result.json`) |
| `--quiet` / `-q` | off | Suppress terminal output except errors |

### File matching policy

When `--msa-dir` or `--tree-dir` are accepted, all subcommands use suffix-agnostic logical locus-name matching from `phyloai/core/file_matching.py`:

| File | Logical locus |
|------|---------------|
| `gene1.fa` | `gene1` |
| `gene2.v1.ALI` | `gene2.v1` |
| `gene3.treefile` | `gene3` |
| `gene4.fa.treefile` | `gene4.fa`, then `gene4` |

Every regular non-empty file is scanned; format is validated when parsed, not by extension. Ambiguous tree matches (where filename candidates are already occupied) cause an error. `phyloai pretree metrics` uses the same helpers, so `metrics` and `filter` behave identically for non-standard naming.

---

## `filter taper` — TAPER Error-Site Masking

### Purpose

Run TAPER to mask erroneous amino-acid or nucleotide sites within MSAs. TAPER identifies stretches of residues that are unexpectedly divergent relative to the rest of the alignment and replaces them with `X` (AA) or `N` (NT). Only newly introduced masks are counted; original ambiguity characters in the input are ignored.

This is site-level quality control: the locus is retained, but problematic positions are neutralized. TAPER does **not** remove loci (use `filter metrics`) or remove taxa (use `filter treeshrink`).

### Usage

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

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--msa-dir` | required | Directory of aligned MSA files. Every regular non-empty file is scanned. |
| `--nt-dir` | — | Codon-aligned NT MSA directory for AA+CDS mode (NT length == 3 × AA length). |
| `--seq-type` | `auto` | `AA`, `NT`, or `auto` (detects from first file: `EFILPQWYZ` → AA). |
| `--cutoff` | 3 | TAPER correction cutoff (≥1). Lower = more aggressive. 3 is TAPER's default; use 1-2 for noisy data, 5+ for conservative masking. |
| `--taper-path` | — | Explicit path to `correction_multi.jl`. Uses bundled copy by default. |
| `--julia-path` | — | Explicit Julia executable. Resolved via PATH; verify with `phyloai doctor`. |
| `--tool-args` | — | Additional TAPER flags passed verbatim. Must not include `-m`, `-a`, `-c`, `-l`, input path, or output redirection (these are managed; error if present). |
| `--threads` / `-t` | 4 | Worker processes, one locus per worker. Uses `ProcessPoolExecutor` with checkpoint/resume (same pattern as `pretree align` and `pretree trim`). |
| `--show-masked-sites` | off | Add `masked_taxa_detail` column to `filter_decisions.csv` (`taxonA:3; taxonB:5`). |
| `--table-format` | `csv` | Format for `retained_loci`, `dropped_loci`, `filter_decisions`. |
| `--resume` | off | Resume from `checkpoint.json`; parameters must match; completed loci skipped if output passes validation. |
| `--overwrite` | off | Mutually exclusive with `--resume`. |
| `--dry-run` | off | Show detected mode, paired loci, command template, output layout. No files. |
| `--output-dir` / `-o` | `runs/pretree/filter/taper` | Output directory. |
| `--quiet` / `-q` | off | Suppress terminal output except errors. |

### Inputs

Three operating modes:

| Mode | Input | Output |
|------|-------|--------|
| AA-only | AA MSA files in `--msa-dir` | Masked AA → `seqs/` |
| NT-only | NT MSA files + `--seq-type NT` | Masked NT → `seqs/` |
| AA+CDS | AA MSA + `--nt-dir` codon-aligned NT | Masked AA → `seqs/faa/`, projected CDS → `seqs/fna/` |

AA+CDS mode requires that: NT records form valid codon MSAs (equal length, divisible by 3); AA and NT taxa match exactly per locus; AA length == NT length / 3.

Projection rules for AA+CDS:
- Original `X` in input → unchanged (not counted as a TAPER mask)
- Standard AA → `X` by TAPER → corresponding codon replaced with `NNN`
- Gap `-` → `X` → warning; no CDS change (defensive check, not expected)

### Outputs

```
runs/pretree/filter/taper/
├── seqs/                              (or seqs/faa/ + seqs/fna/ for AA+CDS)
├── retained_loci.csv|tsv
├── dropped_loci.csv|tsv               (locus, reason)
├── filter_decisions.csv|tsv           (locus, status, new_masked_sites, masked_taxa_count,
│                                       masked_taxa_detail when --show-masked-sites)
├── checkpoint.json                    (internal; only with --resume)
├── logs/
│   ├── gene1.log
│   └── ...
├── result.json
```

Terminal output: two Rich tables — Filter Results (input/retained/dropped/masked loci/taxa/sites) and Retained MSA Statistics (MSA count, total/mean/min/max alignment length, mean taxa). Julia version auto-detected via `julia -v` and recorded in `result.json`.

### Examples

```bash
# Default AA masking
phyloai pretree filter taper --msa-dir ./trimmed

# Aggressive masking for noisy data
phyloai pretree filter taper --msa-dir ./trimmed --cutoff 1 --threads 8

# Conservative masking
phyloai pretree filter taper --msa-dir ./trimmed --cutoff 5

# NT-only mode
phyloai pretree filter taper --msa-dir ./trimmed_nt --seq-type NT

# AA+CDS: mask AA, project to codon-aligned NT
phyloai pretree filter taper --msa-dir ./trimmed_aa --nt-dir ./trimmed_fna

# Resume after interruption
phyloai pretree filter taper --msa-dir ./trimmed --resume

# Include per-taxon mask detail for inspection
phyloai pretree filter taper --msa-dir ./trimmed --show-masked-sites
```

### Warnings and Errors

| Condition | Behaviour |
|-----------|-----------|
| `--nt-dir` with `--seq-type NT` | Exit 1 |
| `--threads` < 1 | Exit 1 |
| `--resume` + `--overwrite` | Exit 1 |
| Julia not found | Exit 3 |
| Non-empty output directory without `--overwrite` or `--resume` | Exit 1 |
| No valid MSA files in `--msa-dir` | Exit 1 |
| TAPER exits non-zero for a locus | Locus skipped; reason in `dropped_loci.csv` |
| TAPER output missing or fails FASTA validation | Locus skipped |
| All loci fail | Exit 2 |

### Notes

TAPER is always the first filtering step because masking should happen before tree inference. After masking, compute per-locus metrics with `phyloai pretree metrics` on the masked MSAs (and optionally on re-inferred gene trees), then apply `filter metrics` or `filter cluster`.

`--resume` is supported because masking large numbers of loci is compute-intensive. The checkpoint follows the same pattern as `pretree align` and `pretree trim`.

---

## `filter treeshrink` — TreeShrink Taxon Pruning

### Purpose

Run TreeShrink to detect and remove outlier long-branch taxa from gene trees. TreeShrink uses a statistical test to identify taxa with unexpectedly long branches across multiple trees jointly. When `--msa-dir` is provided, matching MSAs are also shrunk to remove the same pruned taxa.

This is taxon-level filtering: taxa are removed from specific gene trees. The locus is retained. TreeShrink does **not** remove entire loci (use `filter metrics`).

### Usage

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

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--tree-dir` | required | Directory of gene tree files. Every regular non-empty file scanned. |
| `--msa-dir` | — | Optional MSA directory. Files paired by logical locus name; shrunk alongside trees. |
| `--threshold` | 0.05 | TreeShrink false-positive threshold (`-q`). Smaller = more taxa removed. 0.05 is TreeShrink's default. |
| `--treeshrink-mode` | `auto` | `auto` omits `-m` (TreeShrink default). `per-gene`: independent per gene. `all-genes` / `per-species`: pools cross-gene info. |
| `--treeshrink-path` | — | Explicit path to `run_treeshrink.py`. Resolved via PATH. |
| `--tool-args` | — | Additional TreeShrink flags. Must not include `-i`, `-t`, `-a`, `-q`, `-m`, `-o`, `-O`. |
| `--keep-work-dir` | off | Retain per-gene working directory under `output_dir/work/` for debugging. |
| `--output-dir` / `-o` | `runs/pretree/filter/treeshrink` | Output directory. |
| `--table-format` | `csv` | Format for auxiliary tables. |
| `--overwrite` | off | Delete and recreate output directory. |
| `--dry-run` | off | Print resolved command and locus count. |
| `--quiet` / `-q` | off | Suppress terminal output. |

### Inputs

`--tree-dir` is required. TreeShrink is invoked once across the entire dataset (not per locus) because its statistical model can pool information from multiple trees. PhyloAI creates a per-gene work layout in a temporary directory:

```
<work_dir>/input/
├── gene1/
│   ├── input.tree
│   └── input.fasta     (only when --msa-dir provided)
├── gene2/
│   ├── input.tree
│   └── input.fasta
```

### Outputs

```
runs/pretree/filter/treeshrink/
├── trees/                              (shrunk gene trees)
├── seqs/                               (only when --msa-dir provided)
├── retained_loci.csv|tsv
├── modified_loci.csv|tsv               (loci where ≥1 taxon pruned)
├── dropped_loci.csv|tsv                (loci with missing/invalid output)
├── removed_taxa.csv|tsv                (locus, taxon per row)
├── filter_decisions.csv|tsv            (locus, status, removed_count)
├── work/                               (only with --keep-work-dir)
├── logs/treeshrink.log                 (single shared tool stderr)
└── result.json
```

Decision categories: retained (incl. unmodified), modified (taxa removed), dropped (output missing).

Terminal output: Filter Results table (input/retained/modified/dropped/taxa removed) + Retained MSA Statistics table (when `--msa-dir` provided). A contextual tip reminds users that filtered alignments may be used to re-construct phylogenetic trees, which are possibly more accurate than those pruned by TreeShrink.

### Examples

```bash
# Basic taxon pruning
phyloai pretree filter treeshrink --tree-dir ./genetrees

# Trees + matching MSAs, conservative threshold
phyloai pretree filter treeshrink \
  --tree-dir ./genetrees --msa-dir ./trimmed --threshold 0.1

# Per-species mode (cross-gene pooling)
phyloai pretree filter treeshrink \
  --tree-dir ./genetrees --treeshrink-mode per-species

# Debug output
phyloai pretree filter treeshrink --tree-dir ./genetrees --keep-work-dir
```

### Warnings and Errors

| Condition | Behaviour |
|-----------|-----------|
| `run_treeshrink.py` not found | Exit 3 |
| No valid tree files in `--tree-dir` | Exit 1 |
| Ambiguous locus matching between trees and MSAs | Exit 1 with details |
| TreeShrink exits non-zero | All loci marked failed |
| All loci fail | Exit 2 |
| Non-empty output directory without `--overwrite` | Exit 1 |

### Notes

`--resume` and `--threads` are not supported: TreeShrink runs once across the entire dataset, and per-locus parallelization would change the statistical model. TreeShrink's `-q` threshold controls the false-positive rate; 0.05 means ~5% chance of incorrectly removing a taxon.

After TreeShrink, you may want to re-infer gene trees on the shrunk MSAs for more accurate topologies, then compute metrics and filter.

---

## `filter metrics` — Metric Rule Filtering

### Purpose

Filter whole loci by explicit numeric or string conditions on a metrics CSV/TSV table (typically the output of `phyloai pretree metrics`). All conditions in `--keep` are combined with AND logic: a locus must satisfy every condition to be retained.

This is locus-level filtering: the entire gene is kept or discarded. `filter metrics` does **not** compute metrics (use `pretree metrics`) or filter at the site/taxon level (use `filter taper` or `filter treeshrink`).

### Usage

```bash
phyloai pretree filter metrics \
  --table <metrics.csv|metrics.tsv> \
  --keep "col>=val,col<=val,..." \
  [--input-format auto|csv|tsv] \
  [--loci-column loci] \
  [--msa-dir <msa_dir>] [--tree-dir <tree_dir>] [--copy] \
  [--output-dir runs/pretree/filter/metrics] \
  [--table-format csv|tsv] [--dry-run] [--overwrite]
```

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--table` | required | Path to metrics CSV/TSV. Delimiter auto-detected unless `--input-format` specified. |
| `--keep` | required | Comma-separated AND conditions. Operators: `>=`, `>`, `<=`, `<`, `==`, `!=`. Only `==`/`!=` on string columns. |
| `--input-format` | `auto` | `csv`, `tsv`, or `auto`. Override when auto-detection is ambiguous. |
| `--loci-column` | `loci` | Column name holding the locus identifier. |
| `--msa-dir` | — | MSA directory for retained-MSA statistics (terminal + `result.json`). |
| `--tree-dir` | — | Tree directory for `--copy` mode. |
| `--copy` | off | Copy retained MSAs/trees to output dir. Requires `--msa-dir` or `--tree-dir`. |
| `--output-dir` / `-o` | `runs/pretree/filter/metrics` | Output directory. |
| `--table-format` | `csv` | Format for `retained_loci`, `dropped_loci`, `filter_decisions`. |
| `--overwrite` | off | Delete and recreate output directory. |
| `--dry-run` | off | Parse rules and report pass/fail counts without writing files. |
| `--quiet` / `-q` | off | Suppress terminal output. |

### Inputs

The `--table` file must be CSV or TSV with a header row. The locus identifier column (default `loci`) identifies each row. Delimiter auto-detection inspects the first 1024 bytes. Empty files cause an error.

Rule syntax:

```
column operator value
```

```bash
# Numeric thresholds (AND only)
--keep "dvmc>=0,dvmc<=0.3,average_BS>=0.8"

# Mixed numeric + string
--keep "DataType==AA,num_sites>=300"

# Simple single condition
--keep "num_sites>=1000"
```

Using `>=`/`>`/`<=`/`<` on a string column exits with an error. OR logic is not supported in v1.

### Outputs

```
runs/pretree/filter/metrics/
├── retained_loci.csv|tsv
├── dropped_loci.csv|tsv
├── filter_decisions.csv|tsv
├── seqs/                              (only with --copy --msa-dir)
├── trees/                             (only with --copy --tree-dir)
├── result.json
```

Terminal output: Filter Results table (total/retained/dropped) + Retained MSA Statistics table when `--msa-dir` is provided.

`result.json.key_results.condition_failure_counts` maps each condition to how many loci it rejected — use this to identify the most restrictive threshold.

### Examples

```bash
# Basic quality filter
phyloai pretree filter metrics \
  --table ./metrics/metrics.csv \
  --keep "dvmc<=0.3,average_BS>=0.8"

# Mixed numeric + string with file copy
phyloai pretree filter metrics \
  --table ./metrics.csv \
  --keep "DataType==AA,num_sites>=300" \
  --copy --msa-dir ./trimmed

# Dry-run to explore thresholds before writing
phyloai pretree filter metrics \
  --table ./metrics.csv \
  --keep "num_sites>=500,average_BS>=0.7" \
  --dry-run

# Custom locus column name
phyloai pretree filter metrics \
  --table ./table.tsv \
  --keep "average_BS>=0.9" \
  --loci-column gene_id

# Compare two strategies (different output dirs)
phyloai pretree filter metrics \
  --table ./metrics.csv --keep "average_BS>=0.8" \
  -o ./runs/strategy_conservative

phyloai pretree filter metrics \
  --table ./metrics.csv --keep "average_BS>=0.5" \
  -o ./runs/strategy_lenient
```

### Warnings and Errors

| Condition | Behaviour |
|-----------|-----------|
| `--table` does not exist | Exit 1 |
| `--table` is empty | Exit 1 |
| Malformed `--keep` syntax | Exit 1 with parse error detail |
| Unknown column referenced in `--keep` | Exit 1 |
| Numeric operator (`>=`, `>`, `<=`, `<`) on string column | Exit 1 |
| `--copy` without `--msa-dir` or `--tree-dir` | Exit 1 |
| No loci satisfy all conditions | Result reports 0 retained (not an error) |
| Non-empty output directory without `--overwrite` | Exit 1 |

### Notes

`filter metrics` is deliberately separate from `pretree metrics` computation so you can explore threshold combinations without re-computing metrics. Use `--dry-run` to iterate quickly, then apply the final thresholds with `--copy` to produce filtered files.

Without `--copy`, only decision tables are written — fast and disk-friendly for threshold exploration. The `condition_failure_counts` in `result.json` shows exactly which condition drops the most loci.

To implement OR logic, run the command twice with different `--keep` and `--output-dir` values.

---

## `filter symtest` — Symmetry Test Filtering

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
| `--iqtree-path` | — | Explicit path to iqtree binary. IQ-TREE3 >= 2.3.0 required. |
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

Temporary files are cleaned up after the run.

### Outputs

```
runs/pretree/filter/symtest/
├── seqs/                              (retained MSAs)
├── trees/                             (only when --tree-dir provided)
├── retained_loci.csv|tsv
├── dropped_loci.csv|tsv               (locus, reason)
├── filter_decisions.csv|tsv           (locus, status, p_value, symtest_type,
│                                       sym_pval, mar_pval, int_pval,
│                                       sym_sig, sym_non, mar_sig, mar_non,
│                                       int_sig, int_non)
├── logs/symtest.log                    (single shared tool stderr)
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

### Notes

Symmetry testing should be run after alignment and trimming but before supermatrix concatenation, since violations of stationarity or homogeneity can bias phylogenetic inference. The `--symtest-type` default (combined Sym test) is the most general and widely applicable.

References: Naser-Khdour et al. (2019) "Assessing the Goodness of Fit of Phylogenetic Models..." doi:10.1093/gbe/evz193.

---

## `filter cluster` — Cluster-Based Exploration

### Purpose

Group loci by their metric profiles using dimensionality reduction (PCA or UMAP) followed by hierarchical clustering. This is primarily an exploratory tool: by default it writes clusters, diagnostic plots, and per-cluster metric summaries without removing any loci.

Use `--drop-outlier-clusters auto` to optionally remove the worst-performing clusters. `filter cluster` does **not** apply rule-based filtering (use `filter metrics`) or mask/prune individual sites/taxa (use `filter taper` or `filter treeshrink`).

### Usage

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
  [--plot-metrics-cols N] [--plot-label-angle 45] \
  [--outlier-boxplot-cols N] \
  [--umap-n-neighbors 15] [--umap-min-dist 0.001] \
  [--umap-replicates 1] [--umap-random-state 42] \
  [--threads 1] \
  [--msa-dir <msa_dir>] [--tree-dir <tree_dir>] [--copy] \
  [--output-dir runs/pretree/filter/cluster] \
  [--dry-run] [--overwrite]
```

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--table` | required | Metrics CSV/TSV. |
| `--input-format` | `auto` | `csv`, `tsv`, or `auto`. |
| `--metrics` | all | Comma-separated metric columns. Default: all numeric except locus ID, `DataType`, constant columns, and `--exclude-regex` matches. |
| `--exclude-regex` | — | Repeatable. `--exclude-regex '^freq' --exclude-regex '^sd_'` or combined: `--exclude-regex '^(freq\|sd_)'`. |
| `--reduction` | `pca` | `pca` (3 components, deterministic, `scikit-learn`) or `umap` (3 components, stochastic, optional `umap-learn`). |
| `--n-clusters` | auto | Fixed cluster count. Auto-selected via multi-metric voting when omitted. |
| `--max-clusters` | auto | Upper bound for auto-selection. Default: `min(30, max(6, ceil(sqrt(n_loci)/3)))`. |
| `--cluster-linkage` | `ward` | `ward` (minimizes within-cluster variance), `average`, `complete`, `single`. |
| `--cluster-distance` | `euclidean` | `euclidean`, `cosine`, `manhattan`. Ward requires Euclidean. |
| `--drop-outlier-clusters` | `none` | `none`: diagnostics only. `auto`: remove worst clusters. |
| `--outlier-metric` | `average_BS` | Metric for ranking clusters when dropping. |
| `--outlier-direction` | `low` | `low`: smaller values worse. `high`: larger values worse. |
| `--max-drop-fraction` | 0.2 | Maximum fraction of loci removable (0.0–1.0). |
| `--plot-metrics-cols` | 2 | Boxplot columns per figure (cluster_metric_boxplots). Rows auto-calculated. |
| `--plot-label-angle` | 45.0 | X-axis label rotation in plots. |
| `--outlier-boxplot-cols` | 4 | Boxplot columns per figure (outlier_comparison_boxplots). Rows auto-calculated. |
| `--umap-n-neighbors` | 15 | UMAP local/global balance (PCA: ignored). |
| `--umap-min-dist` | 0.001 | UMAP point packing (PCA: ignored). |
| `--umap-replicates` | 1 | UMAP runs; best selected by cluster-validation rank-sum scoring (PCA: ignored). |
| `--umap-random-state` | 42 | Random seed for reproducible UMAP. Only applied when `--umap-replicates 1`; when >1, no seed is set so `--threads` can parallelize. |
| `--threads` | 1 | CPU threads for UMAP n_jobs. Only used when `--reduction umap` and `--umap-replicates > 1`. |
| `--msa-dir` / `--tree-dir` | — | Input directories for `--copy` mode. |
| `--output-dir` / `-o` | `runs/pretree/filter/cluster` | Output directory. |
| `--table-format` | `csv` | Format for all output tables. |
| `--overwrite` | off | Delete and recreate output directory. |
| `--dry-run` | off | Show selected features, reduction, cluster range, drop plan. |
| `--quiet` / `-q` | off | Suppress terminal output. |

### Inputs

Feature selection: all numeric columns from the input table, excluding locus ID, `DataType`, constant columns, and `--exclude-regex` matches. All features are z-score scaled before reduction.

PCA produces `PC1`/`PC2`/`PC3` (3 components via `sklearn.decomposition.PCA`). UMAP produces `UMAP1`/`UMAP2`/`UMAP3` (requires `pip install umap-learn`; missing dependency exits with install hint).

Cluster count selection (when `--n-clusters` not set): evaluate `k=2..max_clusters`. Three internal validation metrics vote — silhouette (higher), Calinski-Harabasz (higher), Davies-Bouldin (lower). Ties broken by higher silhouette, then smaller `k`. UMAP replicate selection uses rank-sum scoring across the same three metrics.

Outlier dropping (when `--drop-outlier-clusters auto`): rank clusters by mean `--outlier-metric`, drop worst until cumulative fraction exceeds `--max-drop-fraction`.

### Outputs

Core (always written):
```
runs/pretree/filter/cluster/
├── 01-input/
│   └── features_used.csv|tsv          (column, included, reason)
├── 02-reduction/
│   ├── reduction.csv|tsv              (PC1/PC2/PC3 or UMAP1/UMAP2/UMAP3 per locus)
│   ├── cluster_selection.csv|tsv      (k, silhouette, calinski_harabasz, davies_bouldin)
│   └── umap_replicates.csv|tsv        (only when UMAP replicates > 1)
├── 03-clustering/
│   ├── clusters.csv|tsv               (locus, cluster)
│   ├── cluster_summary.csv|tsv        (per-cluster size)
│   └── cluster_loci/cluster_*.csv|tsv  (loci in each cluster)
├── 04-diagnostics/
│   ├── cluster_metric_means.csv|tsv   (per-cluster mean for each numeric metric)
│   └── plots/
│       ├── cluster_2d.pdf             (first 2 reduced dimensions)
│       ├── cluster_3d.pdf             (3D scatter)
│       ├── cluster_metric_heatmap.pdf (z-score heatmap: metrics × clusters)
│       └── cluster_metric_boxplots.pdf (per-metric distributions by cluster)
├── result.json
```

With `--drop-outlier-clusters auto`, additionally:
```
├── 05-outlier-drop/
│   ├── retained_loci.csv|tsv
│   ├── dropped_loci.csv|tsv
│   ├── filter_decisions.csv|tsv
│   ├── outlier_comparison.csv|tsv     (normal vs outlier: mean, median, sd, count per metric)
│   ├── outlier_wilcoxon.csv|tsv       (Mann-Whitney U p-values per metric)
│   ├── plots/
│   │   └── outlier_comparison_boxplots.pdf  (all metrics, * p<0.05 ** p<0.01 *** p<0.001)
│   ├── seqs/                          (only with --copy --msa-dir)
│   └── trees/                         (only with --copy --tree-dir)
```

### Examples

```bash
# Exploratory: see how loci cluster by metric profiles
phyloai pretree filter cluster --table ./metrics/metrics.csv

# UMAP with fixed cluster count
phyloai pretree filter cluster \
  --table ./metrics.csv --reduction umap --n-clusters 5

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

### Warnings and Errors

| Condition | Behaviour |
|-----------|-----------|
| `--table` does not exist | Exit 1 |
| `--table` is empty or has no numeric columns | Exit 1 |
| `--reduction umap` and `umap-learn` not installed | Exit 1 with `pip install umap-learn` hint |
| `--cluster-linkage ward` + non-Euclidean `--cluster-distance` | Exit 1 |
| `--n-clusters` < 2 or > number of loci | Exit 1 |
| `--copy` without `--msa-dir` or `--tree-dir` | Exit 1 |
| `--drop-outlier-clusters auto` with `--copy` but no loci dropped | Copy is a no-op (warning) |
| Non-empty output directory without `--overwrite` | Exit 1 |

### Notes

Without `--drop-outlier-clusters auto`, the command is read-only — no loci are removed, only diagnostics are written. This is intentional: cluster interpretation should remain user-guided; automatic removal is potentially risky.

`features_used.csv` is the audit trail — it shows every column, whether included in the feature set, and why excluded.

`plots/cluster_metric_heatmap.pdf` is a z-score standardized heatmap where:
- **Rows** = clusters, **Columns** = metrics, **Color** = how far above (red) or below (blue) a cluster's mean is from the global mean across all clusters, measured in standard deviations.
- A deep-red cell means that cluster scores much higher than average on that metric; deep-blue means much lower.
- Each cell is annotated with the exact z-score value (e.g. `+1.53` or `−0.87`).
- Use this to quickly identify: which clusters score well/poorly on which metrics, whether outlier clusters have systematically low values on quality proxies like `average_BS`, and whether cluster separation is driven by a few dominant metrics.

PCA is the default reduction because it is deterministic, stable, and requires only `scikit-learn`. UMAP is available for exploring non-linear structure but adds the optional `umap-learn` dependency and stochasticity. With `--umap-replicates 1`, a fixed `--umap-random-state` seed ensures reproducibility. With `--umap-replicates > 1`, no seed is set and `--threads` can parallelize UMAP.

Terminal output shows step-by-step progress (feature selection → reduction → clustering → diagnostics → outlier detection) and describes each output file. After `--drop-outlier-clusters auto`, the screen also shows retained MSA statistics when `--msa-dir` is provided.

---

## result.json schema

All subcommands follow the same result format:

```json
{
  "status": "success | error",
  "command": "phyloai pretree filter ...",
  "wall_time": 1.23,
  "tool_versions": {},
  "params": {},
  "key_results": {},
  "error": null,
  "data": {}
}
```

Mode-specific `key_results` and `data`:

| Subcommand | key_results | data |
|------------|-------------|------|
| `taper` | n_input, n_retained, n_dropped, masked_loci, total_masked_taxa, total_masked_aa_sites | retained_msa_stats, dry_run_cmds |
| `treeshrink` | n_input, n_retained, n_modified, n_dropped, n_removed_taxa_total | retained_loci, modified_loci, dropped_loci, removed_taxa, retained_msa_stats |
| `metrics` | n_total, n_retained, n_dropped, condition_failure_counts | copied_msa, copied_tree, retained_msa_stats |
| `cluster` | n_loci, n_valid_loci, n_features, n_clusters, reduction, selected_umap_replicate, n_retained, n_dropped | features, cluster_sizes, drop_clusters, retained_loci, retained_msa_stats, plot_paths, umap_replicates |

Tool stderr for external-tool subcommands (`taper`, `treeshrink`, `symtest`) is written to per-locus `logs/<locus>.log` (taper) or a shared `logs/<tool>.log` file (treeshrink, symtest). Pure-Python subcommands (`metrics`, `cluster`) write no external log.
