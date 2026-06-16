# phyloai pretree filter Design Specification

**Date:** 2026-06-15  
**Status:** Draft for review  
**Reference:** `docs/superpowers/specs/2026-06-07-phyloai-design.md`, `docs/superpowers/specs/2026-06-14-pretree-metrics-design.md`, `docs/superpowers/specs/2026-06-12-pretree-trim-design.md`, `ref/scripts/loci_filtering_alignment-based.sh`, `ref/scripts/loci_filtering_tree-based.sh`, `ref/scripts/outlier_detection.R`, `ref/cluster.txt`, TAPER README, TreeShrink README

---

## 1. Purpose

`phyloai pretree filter` runs after `pretree metrics` and before `pretree concat`. It provides four filtering workflows:

1. **TAPER site masking**: mask erroneous sequence stretches within MSAs.
2. **TreeShrink taxon pruning**: remove outlier taxa from gene trees and optionally matching MSAs.
3. **Metric rule filtering**: remove whole loci by explicit conditions on `metrics.csv`-like tables.
4. **Cluster-based exploration/filtering**: group loci by metric profiles using PCA or UMAP plus hierarchical clustering, then optionally remove small worst outlier clusters.

The module keeps metric computation separate from filtering decisions. It reads `pretree metrics` output where appropriate, writes structured decisions, and produces filtered or copied MSA/tree directories only when the selected mode needs them.

What it does **not** do:
- compute marker metrics; use `pretree metrics`
- infer gene trees after masking; use later tree commands
- concatenate retained MSAs; use `pretree concat`
- implement inconsistent-gene likelihood workflows from the legacy scripts in the first version

---

## 2. Architecture

```
phyloai/core/file_matching.py       # shared logical locus-name parsing and MSA/tree pairing
phyloai/pretree/filter.py           # core library: taper, treeshrink, metrics rules, clustering, outputs
phyloai/cli/commands/pretree.py     # CLI: filter group and subcommands
docs/commands/pretree-filter.md     # user-facing command documentation
```

**Shared components:**
- `core/env.py`: tool detection for Julia, bundled `correction_multi.jl`, and `run_treeshrink.py`
- `core/checkpoint.py` and `pretree/checkpoint_helpers.py`: used by TAPER only
- `core/sequence_normalization.py`: sequence type detection and codon validation
- `core/sequence_output_validation.py`: FASTA/MSA validation
- `core/file_matching.py`: new shared helper used by both `pretree filter` and `pretree metrics`

**External tools:**
- TAPER: bundled `phyloai/bundled/TAPER-1.0.0/correction_multi.jl`, executed by Julia
- TreeShrink: user-installed `run_treeshrink.py`

**Python dependencies:**
- No new dependency for `taper`, `treeshrink`, or `metrics`
- Add `scikit-learn` for `cluster --reduction pca` and clustering validation
- Use optional `umap-learn` only when `cluster --reduction umap` is requested; missing dependency exits with an install hint
- Do not add `pandas`; use standard-library `csv` plus `numpy`
- Use existing `matplotlib` and `scipy` for plots and tests

---

## 3. Command Structure

`pretree filter` is a Click group with four subcommands:

```bash
phyloai pretree filter taper
phyloai pretree filter treeshrink
phyloai pretree filter metrics
phyloai pretree filter cluster
```

Subcommands have separate input validation because their required inputs, external tools, and output semantics differ.

### 3.1 Common Output Options

All subcommands support:

| Parameter | Type | Default | Notes |
|---|---|---|---|
| `--output-dir`, `-o` | path | `runs/pretree/filter` | Output directory |
| `--table-format` | `csv`\|`tsv` | `csv` | Format for auxiliary tables; file suffix follows the format |
| `--dry-run` | flag | off | Validate inputs and show planned work without writing output files |
| `--quiet`, `-q` | flag | off | Suppress terminal output except errors |
| `--overwrite` | flag | off | Replace a non-empty output directory |

All successful non-dry-run invocations write `result.json` and `filter.log`. `--dry-run` writes no files.

---

## 4. Shared File Matching

This command follows the global file-matching policy in `docs/superpowers/specs/2026-06-07-phyloai-design.md`. Flat MSA and tree directories are matched by logical locus name, using suffix-agnostic dot-boundary parsing rather than command-local suffix whitelists.

Illustrative examples under the global policy:

| File | Matching candidate(s) |
|---|---|
| `gene.fa` | `gene` |
| `gene.v1.ALI` | `gene.v1` |
| `gene.treefile` | `gene` |
| `gene.fa.treefile` | `gene.fa`, then `gene` |
| `gene.v1.fa.treefile` | `gene.v1.fa`, then `gene.v1` |

If one-dot and two-dot tree reductions would match different loci, the command must fail with an explicit ambiguity error instead of choosing automatically.

### 4.1 Shared Helper Adoption

Add `phyloai/core/file_matching.py` with helpers such as:
- `logical_msa_locus_name(path: Path) -> str`
- `logical_tree_locus_name(path: Path) -> str`
- `scan_msa_dir(path: Path) -> dict[str, Path]`
- `scan_tree_dir(path: Path) -> dict[str, Path]`
- `pair_loci(msa_dir: Path | None, tree_dir: Path | None) -> PairingResult`

`pretree metrics` should use the same helper so `metrics` and `filter` behave identically for files such as `gene.fa.treefile`, uppercase suffixes, and other non-standard naming patterns.

---

## 5. TAPER Site Masking

### 5.1 CLI

```bash
phyloai pretree filter taper \
  --msa-dir <aa_or_nt_msa_dir> \
  [--nt-dir <codon_aligned_nt_msa_dir>] \
  [--seq-type AA|NT|auto] \
  [--cutoff 3] \
  [--taper-path <correction_multi.jl>] \
  [--julia-path <julia>] \
  [--tool-args "..."] \
  [--threads 4] [--resume] [--dry-run] [--overwrite]
```

### 5.2 Inputs and Modes

`--msa-dir` is required. `--tree-dir` is not supported for TAPER because masking occurs before gene tree inference.

Operating modes:
- **AA-only**: AA MSA input, no `--nt-dir`; output masked AA MSAs to `seqs/`.
- **NT-only**: NT MSA input, no `--nt-dir`; output masked NT MSAs to `seqs/`.
- **AA+CDS**: AA MSA input plus `--nt-dir` containing codon-aligned NT MSAs; output masked AA to `seqs/faa/` and projected masked CDS to `seqs/fna/`.

`--nt-dir` must contain codon-aligned NT MSAs paired by logical locus name. Raw unaligned CDS is not accepted.

### 5.3 Commands

AA command:

```bash
julia correction_multi.jl -c <cutoff> INPUT > OUTPUT
```

NT command:

```bash
julia correction_multi.jl -m N -a N -c <cutoff> INPUT > OUTPUT
```

`--cutoff` maps to TAPER `-c`, default `3`. It must be `>=1`. Help text should say: lower values are more aggressive, values around `1-10` are typical, `10` is conservative, and larger values are allowed.

`--tool-args` can pass other TAPER strategy options but cannot include PhyloAI-managed flags: `-m`, `-a`, `-c`, `-l`, input path, or output redirection.

### 5.4 AA+CDS Projection

Run TAPER on the AA MSA. Only newly introduced AA `X` masks are projected back to CDS as `NNN`.

Validation before running:
- NT records form a valid codon MSA: equal lengths and length divisible by 3
- AA and NT taxa match exactly for each locus
- AA alignment length equals NT alignment length / 3

Projection rules:

| Original AA | TAPER AA | CDS action |
|---|---|---|
| `X` | `X` | no change; original ambiguity is not a TAPER mask |
| standard AA | `X` | replace corresponding codon with `NNN` |
| `-` | `X` | warning; no CDS change |
| `X` | non-`X` | warning or failure depending severity |

The last two cases are defensive checks; they are not expected in normal TAPER output.

### 5.5 Execution and Resume

TAPER runs one locus per worker, following existing `pretree align` and `pretree trim` patterns.

- `--threads` controls worker count.
- `--resume` uses per-locus checkpoint tasks keyed by logical locus name.
- A successful resumed task is skipped only if expected output files exist and pass FASTA/MSA validation.
- Per-locus failures are recorded and do not stop other loci.
- The command exits successfully if at least one locus succeeds; it exits with code 1 if all loci fail.

### 5.6 TAPER Outputs

Output directories:
- AA-only or NT-only: `seqs/`
- AA+CDS: `seqs/faa/`, `seqs/fna/`

Tables:
- `retained_loci.csv|tsv`: loci with successful masked output
- `dropped_loci.csv|tsv`: loci that failed or were skipped
- `filter_decisions.csv|tsv`: one row per locus with status, reason, masked site counts, and output paths

Terminal and `result.json` include retained MSA statistics for the generated masked MSAs.

---

## 6. TreeShrink Taxon Pruning

### 6.1 CLI

```bash
phyloai pretree filter treeshrink \
  --tree-dir <gene_tree_dir> \
  [--msa-dir <msa_dir>] \
  [--threshold 0.05] \
  [--treeshrink-mode auto|per-gene|all-genes|per-species] \
  [--treeshrink-path <run_treeshrink.py>] \
  [--tool-args "..."] \
  [--dry-run] [--overwrite] [--keep-work-dir]
```

### 6.2 Inputs and Execution Model

`--tree-dir` is required. `--msa-dir` is optional.

TreeShrink can use information from multiple gene trees jointly, so PhyloAI runs TreeShrink once across the dataset rather than per gene. The first version does not provide `--threads` or `--resume` for TreeShrink.

PhyloAI creates a per-gene working layout in a safe temporary directory by default. `--keep-work-dir` writes and retains the work directory under `output_dir/work/` for debugging.

Per-gene work layout:

```
work_input/
  gene1/
    input.tree
    input.fasta    # only when --msa-dir is provided
  gene2/
    input.tree
    input.fasta
```

Commands:

```bash
run_treeshrink.py -i <work_input> -t input.tree -q <threshold>
run_treeshrink.py -i <work_input> -t input.tree -a input.fasta -q <threshold>
```

`--threshold` maps to TreeShrink `-q`, default `0.05`, and PhyloAI accepts only one value. `--treeshrink-mode auto` omits `-m`; explicit values pass to `-m`.

`--tool-args` can pass TreeShrink strategy options but cannot include PhyloAI-managed flags: `-i`, `-t`, `-a`, `-q`, `-m`, `-o`, `-O`.

### 6.3 Outputs and Semantics

Output directories:
- `trees/`: shrunk gene trees
- `seqs/`: shrunk MSAs only when `--msa-dir` was provided

Decision categories:
- `retained_loci`: loci with valid shrunk outputs, including unmodified loci and loci with pruned taxa
- `modified_loci`: retained loci where TreeShrink removed at least one taxon
- `dropped_loci`: loci with missing or invalid outputs

Tables:
- `retained_loci.csv|tsv`
- `dropped_loci.csv|tsv`
- `modified_loci.csv|tsv`
- `removed_taxa.csv|tsv` with `locus,taxon`
- `filter_decisions.csv|tsv`

Terminal and `result.json` include retained MSA statistics only when `--msa-dir` was provided and shrunk MSAs were generated.

---

## 7. Metrics Rule Filtering

### 7.1 CLI

```bash
phyloai pretree filter metrics \
  --table <metrics.csv|metrics.tsv> \
  --keep "dvmc>=0,dvmc<=0.3,average_BS>=0.8" \
  [--input-format auto|csv|tsv] \
  [--loci-column loci] \
  [--msa-dir <msa_dir>] [--tree-dir <tree_dir>] [--copy] \
  [--dry-run] [--overwrite]
```

### 7.2 Rule Syntax

`--keep` is a comma-separated list of AND conditions. OR logic is not supported in the first version.

Supported operators:

```
>=, >, <=, <, ==, !=
```

Examples:

```bash
--keep "dvmc>=0,dvmc<=0.3,average_BS>=0.8"
--keep "DataType==AA,num_sites>=300,proportion_gaps<=0.5"
```

Rules:
- unknown columns are errors
- malformed conditions are errors
- missing values fail numeric comparisons
- `==` and `!=` can compare strings; numeric parsing is attempted first

### 7.3 Table Input and Output

`--table` accepts CSV or TSV. `--input-format auto` is the default:
- `.tsv` or `.tab` means TSV
- all other suffixes mean CSV

`--table-format csv|tsv` controls output table delimiters and suffixes.

### 7.4 Copy and Summaries

Metrics filtering always writes decision tables unless `--dry-run` is set. It copies files only when `--copy` is set.

Copy behavior:
- if `--copy` and `--msa-dir` are provided, retained MSAs are copied to `seqs/`
- if `--copy` and `--tree-dir` are provided, retained trees are copied to `trees/`
- if `--copy` is set without either input directory, exit with a clear error

When `--msa-dir` is provided, terminal output reports retained input MSA statistics even without `--copy`, because this is useful for threshold exploration.

Tables:
- `retained_loci.csv|tsv`
- `dropped_loci.csv|tsv`
- `filter_decisions.csv|tsv`

Terminal summary:
- total loci in table
- retained and dropped counts
- condition failure counts
- optional MSA/tree copy counts and missing matches
- retained MSA count, total retained alignment length, mean/min/max marker length, and mean taxa count when MSA data are available

---

## 8. Cluster-Based Exploration and Filtering

### 8.1 CLI

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
  [--outlier-metric average_BS] \
  [--outlier-direction low|high] \
  [--max-drop-fraction 0.2] \
  [--plot-metrics-per-page auto|N] \
  [--plot-label-angle FLOAT] \
  [--msa-dir <msa_dir>] [--tree-dir <tree_dir>] [--copy] \
  [--dry-run] [--overwrite]
```

### 8.2 Feature Selection

Default `--metrics all` uses all numeric columns except:
- locus identifier column
- `DataType`
- empty columns
- constant columns
- columns matching any repeatable `--exclude-regex`

`--exclude-regex` can be repeated:

```bash
--exclude-regex "^freq" --exclude-regex "^sd_"
```

A combined regex is also valid:

```bash
--exclude-regex "^(freq|sd_)"
```

Save `features_used.csv|tsv` with included columns and excluded columns plus exclusion reasons.

### 8.3 Reduction

Default reduction is PCA because it is stable, reproducible, and only requires `scikit-learn`. Documentation should recommend trying UMAP for nonlinear structure exploration.

All features are z-score scaled before reduction. Reduction produces three coordinates:
- PCA: `PC1`, `PC2`, `PC3`
- UMAP: `UMAP1`, `UMAP2`, `UMAP3`

PCA:
- `sklearn.decomposition.PCA(n_components=3)`

UMAP:
- optional dependency `umap-learn`
- missing dependency exits with `pip install umap-learn` guidance
- options:
  - `--umap-replicates 1`
  - `--umap-random-state 0`
  - `--umap-n-neighbors 15`
  - `--umap-min-dist 0.001`

For each UMAP replicate, use `random_state + replicate_index`.

### 8.4 Hierarchical Clustering

Default clustering:
- `--cluster-linkage ward`
- `--cluster-distance euclidean`

Ward linkage is recommended because it minimizes within-cluster variance and matches the legacy R workflow (`ward.D2` with Euclidean distance) closely enough for this purpose. Ward requires Euclidean distance; incompatible combinations exit with a clear error.

Supported linkage options:
- `ward`: recommended default; Euclidean only
- `complete`: conservative, compact clusters
- `average`: balanced alternative
- `single`: available for advanced use but prone to chaining

Implementation uses `sklearn.cluster.AgglomerativeClustering`.

### 8.5 Cluster Count Selection

If `--n-clusters` is provided, use it directly.

Otherwise evaluate `k=2..max_clusters`, where default `max_clusters` is:

```text
min(30, max(6, ceil(sqrt(n_loci) / 3)))
```

`--max-clusters` overrides this default.

Cluster count scoring uses three internal validation metrics:
- silhouette: higher is better
- Calinski-Harabasz: higher is better
- Davies-Bouldin: lower is better

Each metric votes for its best `k`. Ties are broken by higher silhouette, then smaller `k`. Save `cluster_selection.csv|tsv`.

This is intentionally not a full Python replacement for R `NbClust(index="all")`; it is a transparent, lightweight selector suitable for exploratory outlier cluster screening.

### 8.6 UMAP Replicate Selection

When `--reduction umap --umap-replicates N` with `N > 1`, each replicate runs reduction, clustering, and cluster-count scoring.

Select the best replicate by rank sum:

```text
rank_desc(silhouette) + rank_desc(Calinski-Harabasz) + rank_asc(Davies-Bouldin)
```

Lower rank sum is better. Ties are broken by higher silhouette, lower Davies-Bouldin, then smaller replicate index.

Save `umap_replicates.csv|tsv` with each replicate's selected `k`, validation metrics, and rank score.

### 8.7 Outlier Cluster Removal

Default `--drop-outlier-clusters none` means no loci are removed. The command only writes clusters, plots, and diagnostics.

When `--drop-outlier-clusters auto` is set:
- rank clusters by mean `--outlier-metric`
- `--outlier-direction low` treats low values as worse
- `--outlier-direction high` treats high values as worse
- add worst clusters to the dropped set until adding the next cluster would exceed `--max-drop-fraction`
- default `--outlier-metric average_BS`
- default `--outlier-direction low`
- default `--max-drop-fraction 0.2`

If `--copy` is set, retained MSAs and/or trees are copied from provided input directories only when outlier dropping is active. If no loci are dropped, copy is a no-op with a warning.

### 8.8 Cluster Outputs

Core tables:
- `features_used.csv|tsv`
- `reduction.csv|tsv`
- `cluster_selection.csv|tsv`
- `clusters.csv|tsv`
- `cluster_summary.csv|tsv`
- `cluster_loci/cluster_*.csv|tsv`

Core plots:
- `cluster_2d.pdf`: first two reduced dimensions, colored by cluster
- `cluster_3d.pdf`: static 3D reduced coordinates, colored by cluster

When outlier dropping is active, also write:
- `retained_loci.csv|tsv`
- `dropped_loci.csv|tsv`
- `filter_decisions.csv|tsv`
- optional copied `seqs/` and/or `trees/` when `--copy` is set

### 8.9 Cluster Metric Diagnostics

Cluster diagnostics are required because the scatter plots alone are insufficient for deciding whether clusters are biologically or analytically useful.

Always write:
- `cluster_metric_means.csv|tsv`: per-cluster `n_loci` and means for every numeric input metric
- `cluster_metric_heatmap.pdf`: standardized cluster mean heatmap across metrics
- `cluster_metric_boxplots_001.pdf`, etc.: per-metric distributions grouped by cluster

Boxplot pagination is adaptive:
- `--plot-metrics-per-page auto|N`, default `auto`
- `n_clusters <= 6`: up to 12 metrics per page
- `7 <= n_clusters <= 12`: up to 6 metrics per page
- `13 <= n_clusters <= 20`: up to 4 metrics per page
- `n_clusters > 20`: 1-2 metrics per page, depending on label readability
- `--plot-label-angle` defaults to `45`, but auto-increases to `60` or `90` when cluster labels are dense
- figure width increases moderately with cluster count but has an upper bound to avoid very large PDFs

When outlier dropping is active, also write:
- `outlier_comparison.csv|tsv`: normal vs outlier mean, median, standard deviation, and count for each metric
- `outlier_wilcoxon.csv|tsv`: Mann-Whitney U / Wilcoxon rank-sum p-values and direction for each metric, using `scipy.stats.mannwhitneyu`
- `outlier_comparison_boxplots_001.pdf`, etc.: per-metric distributions grouped by outlier status

The design intentionally avoids `plotly` in the first version to minimize dependencies.

---

## 9. Retained MSA Statistics

Any subcommand that produces, copies, or evaluates retained MSAs should report basic retained MSA statistics in terminal output and `result.json.data.retained_msa_stats`.

Metrics:
- retained MSA count
- retained concat total length, computed as the sum of retained alignment lengths
- mean marker length
- min/max marker length
- mean taxa count

Applicability:
- `taper`: always, because it generates masked MSAs
- `treeshrink`: when `--msa-dir` is provided and shrunk MSAs are generated
- `metrics`: when `--msa-dir` is provided, even without `--copy`, using retained input MSAs
- `cluster`: only when outlier dropping is active and retained MSAs can be identified; normally paired with `--copy` and `--msa-dir`

If some retained MSA files cannot be read, show a warning and record details in `result.json.data.warnings`.

---

## 10. Logging and Result Schema

All non-dry-run executions write:
- `filter.log`
- `result.json`

`filter.log` includes:
- resolved command(s)
- tool versions where applicable
- stderr
- wall time
- exit code
- stdout only when it is diagnostic text, not primary sequence output

`result.json` follows existing pipeline conventions:

```json
{
  "params": {},
  "data": {},
  "key_results": {},
  "tool_versions": {},
  "wall_time_seconds": 0.0
}
```

Mode-specific key results:
- `taper`: input loci, successful loci, failed loci, masked AA site count, projected CDS codon count, output paths
- `treeshrink`: input loci, retained loci, modified loci, dropped loci, removed taxon count, output paths
- `metrics`: retained/dropped counts, condition failure counts, copied MSA/tree counts, retained MSA statistics
- `cluster`: selected reduction, selected UMAP replicate if applicable, selected `k`, cluster sizes, outlier clusters if any, retained/dropped counts if dropping

---

## 11. Documentation

Implementation must add:
- `docs/commands/pretree-filter.md`
- README command index update

The command documentation should include:
- when to use each subcommand
- practical examples for AA TAPER, NT TAPER, AA+CDS TAPER, TreeShrink with and without MSA, metrics filtering, PCA clustering, UMAP clustering, and outlier cluster dropping
- table format behavior
- file matching rules and recommended suffixes
- warnings that `cluster` is exploratory and does not remove loci unless explicitly requested

The main design spec should be updated only where its old flat `pretree filter` examples conflict with the new subcommand interface.

---

## 12. Testing Strategy

Unit tests:
- logical locus matching, including `gene.fa.treefile`, `gene.fasta.bestTree`, and dotted locus names like `gene.v1.fa.treefile`
- table input delimiter auto-detection
- output table delimiter selection
- metrics condition parsing and evaluation
- retained MSA statistics
- TAPER AA+CDS projection, including original `X` preservation
- cluster feature selection and repeated `--exclude-regex`
- cluster-count selection on synthetic data
- UMAP replicate selection scoring, using mocked replicate metrics if needed

CLI tests:
- each subcommand validates required inputs
- `--dry-run` writes no files
- `metrics --copy` behavior with MSA and tree dirs
- `cluster` no-drop default writes clusters and diagnostics but no retained/dropped filtering tables
- `cluster --drop-outlier-clusters auto` writes retained/dropped decisions

External-tool tests:
- mock TAPER by a small executable/script that transforms FASTA in predictable ways
- mock TreeShrink by writing expected `output.tree`, `output.fasta`, and `output.txt` into the work layout
- do not require Julia, TAPER, or TreeShrink in CI

---

## 13. Design Decisions

| Decision | Rationale |
|---|---|
| Use filter subcommands | Avoids one command with many invalid option combinations |
| TAPER has no `--tree-dir` | Masked MSAs should be used before gene tree inference; copying old trees risks misuse |
| TAPER runs per locus | Matches existing `align/trim` parallel and resume patterns; simpler than chunked `-l` execution |
| TreeShrink runs once across dataset | TreeShrink can use multiple trees jointly; per-gene parallelization would change the statistical model |
| TreeShrink has no `--resume` or `--threads` in v1 | Keeps the first version simple and statistically correct |
| `metrics` and `cluster` use `--copy` | Users can explore thresholds/clusters without copying large directories |
| Rule filtering supports AND only | Keeps parsing transparent and safe for v1 |
| Cluster defaults to PCA | Stable, reproducible, fewer dependencies |
| UMAP is optional but recommended for exploration | Better nonlinear structure discovery, but adds dependency and stochasticity |
| Cluster count uses lightweight validation voting, not NbClust | Python has no lightweight direct NbClust equivalent; voting is transparent and dependency-light |
| Cluster auto-drop is opt-in | Cluster interpretation should remain user-guided; automatic removal is potentially risky |
| Diagnostic metric plots are required | Cluster scatter plots alone do not explain why clusters differ |
