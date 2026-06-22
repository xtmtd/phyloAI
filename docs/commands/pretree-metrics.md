# phyloai pretree metrics

## Purpose

`phyloai pretree metrics` computes molecular marker attributes from MSAs and/or gene trees for statistical exploration and downstream filtering. It runs after `pretree trim` and before `pretree filter` in the standard pipeline.

The command group has three entry points:

- `phyloai pretree metrics` computes metrics, writes `metrics.csv`, generates distribution plots, and creates a correlation heatmap.
- `phyloai pretree metrics plot` re-generates one metric distribution plot from an existing `metrics.csv`.
- `phyloai pretree metrics correlate` re-generates a correlation matrix and heatmap from an existing `metrics.csv`.

It does not run marker filtering, UMAP clustering, or tree inference.

## Usage

Compute metrics from MSA files:

```bash
phyloai pretree metrics --msa-dir ./runs/pretree/trim/seqs/faa
```

Compute MSA and gene-tree metrics together:

```bash
phyloai pretree metrics \
  --msa-dir ./runs/pretree/trim/seqs/faa \
  --tree-dir ./data/gene_trees \
  --output-dir ./runs/pretree/metrics
```

Preview work without writing files:

```bash
phyloai pretree metrics --msa-dir ./msa --dry-run
```

Re-plot one metric with custom styling and Tukey filtering:

```bash
phyloai pretree metrics plot \
  --csv ./runs/pretree/metrics/metrics.csv \
  --metric entropy \
  --tukey-k 1.5 \
  --color '#2E86AB' \
  --title 'Entropy distribution'
```

Re-generate a compact upper-triangle correlation heatmap:

```bash
phyloai pretree metrics correlate \
  --csv ./runs/pretree/metrics/metrics.csv \
  --triangle upper \
  --label-angle 60 \
  --output-dir ./runs/pretree/metrics/correlate
```

Include all numeric columns, including `freq*` and `sd_*`, and draw full-mode cluster rectangles:

```bash
phyloai pretree metrics correlate \
  --csv ./runs/pretree/metrics/metrics.csv \
  --metrics all \
  --triangle full \
  --cluster-rectangles 5
```

## Parameters

### `phyloai pretree metrics`

| Parameter | Default | Description |
|---|---|---|
| `--msa-dir` | none | Directory of aligned FASTA files (`.fa`, `.fasta`, `.fas`, `.fna`, `.faa`, `.aln`). At least one of `--msa-dir` or `--tree-dir` is required. |
| `--tree-dir` | none | Directory of Newick tree files (`.tre`, `.tree`, `.nwk`, `.newick`, `.treefile`, `.bestTree`, `.contree`). |
| `--seq-type` | `auto` | Molecule type: `AA`, `NT`, or per-marker auto-detection. |
| `--outgroup-list` | none | File with one outgroup taxon name per line for DVMC pruning; requires `--tree-dir`. |
| `--ref-tree` | none | Reference species tree for normalized RF distance; requires `--tree-dir`. |
| `--skip-freq-statistics` | off | Skip per-character frequency columns (`freqA`, `freqC`, etc.). |
| `--pseudo-tree-metrics` | off | Compute FastTree-derived pseudo-tree metrics with `_FT` suffix; requires `--msa-dir`. |
| `--fasttree-path` | `FastTree` | Explicit FastTree executable path. |
| `--skip-pairwise-identity` | off | Skip `average_pairwise_identity`; recommended for markers with many taxa. |
| `--round` | `6` | Decimal places for numeric CSV values; range 0-12. |
| `--table-format` | `csv` | Table format for auxiliary tabular outputs: `csv` or `tsv`. All auxiliary tables (`metrics`, `basic_statistics`, `correlation_matrix`) use the same format. |
| `--output-dir`, `-o` | `runs/pretree/metrics` | Output directory for the metrics table, plots, correlation outputs, and `result.json`. |
| `--threads`, `-t` | `4` | Worker process count; must be at least 1. |
| `--dry-run` | off | Validate inputs and show planned work without writing files. |
| `--overwrite` | off | Delete and recreate a non-empty output directory. |
| `--quiet`, `-q` | off | Suppress terminal progress and summary output. |

### `phyloai pretree metrics plot`

| Parameter | Default | Description |
|---|---|---|
| `--csv` | required | Existing `metrics.csv` (or `.tsv`). |
| `--input-format` | `auto` | Table format of the input file: `csv`, `tsv`, or `auto` (detects by content — tab/comma counts — falling back to file extension). |
| `--metric` | required | Exact metric column name to plot. |
| `--bins` | `50` | Histogram bin count; valid range is 1-500. |
| `--xmin` | auto | Force x-axis lower limit. |
| `--xmax` | auto | Force x-axis upper limit. |
| `--tukey-k` | disabled | Tukey's Fences multiplier. When set, filtered loci are written to `<output_dir>/<metric>.tukey_filtered.csv`. |
| `--title` | `Distribution of <metric>` | Plot title. |
| `--xlabel` | metric display name | X-axis label. |
| `--ylabel` | `Density` | Y-axis label. |
| `--color` | `#2E86AB` | Histogram bar fill color. |
| `--fig-width` | `10.0` | Figure width in inches. |
| `--fig-height` | `8.0` | Figure height in inches. |
| `--dpi` | `150` | Output resolution. |
| `--font-size` | `12` | Base font size. |
| `--output-dir`, `-o` | `<csv_parent>/plot_<metric>/` | Output directory for the PDF and `result.json`. |
| `--overwrite` | off | Replace an existing output PDF/directory. |
| `--quiet`, `-q` | off | Suppress terminal output. |

### `phyloai pretree metrics correlate`

| Parameter | Default | Description |
|---|---|---|
| `--csv` | required | Existing `metrics.csv` (or `.tsv`). |
| `--input-format` | `auto` | Table format of the input file: `csv`, `tsv`, or `auto` (detects by content — tab/comma counts — falling back to file extension). |
| `--metrics` | core numeric | Comma-separated metric columns. Use `all` for every numeric column. Omitted means automatic readable core-metric selection. |
| `--include-freq` | off | Include `freq*` columns in automatic metric selection. |
| `--include-sd` | off | Include `sd_*` columns in automatic metric selection. |
| `--method` | `spearman` | Correlation method: `spearman` or `pearson`. |
| `--triangle` | `full` | Matrix display: `full`, `lower`, or `upper`. Lower mode uses left/bottom labels; upper mode uses top/right labels. |
| `--annot` / `--no-annot` | `--no-annot` | Show numeric correlation values inside cells. |
| `--cluster-rectangles` | none | Draw N cluster rectangles on full matrices only. If used with `--triangle lower` or `upper`, PhyloAI warns and ignores it. |
| `--cmap` | `RdBu_r` | Matplotlib colormap. |
| `--fmt` | `.2f` | Numeric format for annotations. |
| `--fig-width` | `12.0` | Figure width in inches. |
| `--fig-height` | `10.0` | Figure height in inches. |
| `--dpi` | `150` | Output resolution. |
| `--font-size` | `10` | Base font size for metric labels. |
| `--label-angle` | `45.0` | Rotation angle for x-axis metric labels in degrees. Useful for dense upper/lower triangle plots. |
| `--title` | none | Optional plot title. |
| `--output-dir`, `-o` | `runs/pretree/metrics/correlate` | Directory for `correlation_heatmap.pdf`, `correlation_matrix.csv`, and `result.json`. |
| `--overwrite` | off | Replace existing correlation outputs. |
| `--quiet`, `-q` | off | Suppress terminal output and non-critical warnings. |

## Inputs

MSA and tree files are paired by logical locus name using the shared global matching policy.

- MSA logical locus: everything before the final `.` in the filename.
- Tree logical locus: try removing the final suffix segment, then the final two suffix segments.
- If exactly one tree candidate matches an MSA locus, that pair is used.
- If both tree candidates match different loci, PhyloAI exits with an ambiguity error instead of guessing.

Examples: `EOG090X002Z.fas` -> `EOG090X002Z`; `EOG090X002Z.fas.treefile` tries `EOG090X002Z.fas` and `EOG090X002Z`.

## Outputs

### Main `metrics` output directory

All auxiliary tables follow `--table-format` (default `csv`, producing `.csv` files). If `--table-format tsv` is given, `.tsv` files are written instead.

| File or directory | Description |
|---|---|
| `metrics.csv` (or `.tsv`) | One row per marker with identifiers, MSA metrics, tree metrics, optional frequency columns, and optional pseudo-tree metrics. |
| `plots/` | One density histogram PDF per numeric metric. |
| `metrics.basic_statistics.csv` (or `.tsv`) | Mean, median, min, max, q25, q75, standard deviation, non-NA count, and total count per metric. |
| `correlate/correlation_heatmap.pdf` | Default compact correlation heatmap from core numeric metrics. |
| `correlate/correlation_matrix.csv` (or `.tsv`) | Correlation matrix with metric names as row and column labels. |
| `result.json` | Structured status, parameters, key counts, warnings, and data paths. |

### Important `metrics.csv` columns

- `loci` and `DataType` identify the marker and molecule type.
- MSA metrics include `num_taxa`, `taxa_occupancy`, `num_sites`, `num_patterns`, `proportion_patterns`, `num_parsimony_sites`, `proportion_parsimony`, `num_singletons`, `proportion_singletons`, `proportion_gaps`, `proportion_invariant`, `entropy`, `bollback`, `pattern_entropy`, `rcfv`, `nrcfv`, `average_pairwise_identity`, and `GC_content`.
- Tree metrics include `average_BS`, `sd_BS`, `total_tree_length`, branch-length summaries, patristic-distance summaries, `evo_rate`, `treeness`, `dvmc`, `saturation`, and `RF_distance`.
- Frequency metrics use `freq*` column names and are omitted when `--skip-freq-statistics` is set.
- Pseudo-tree metrics use `_FT` suffixes and are only produced with `--pseudo-tree-metrics`.

## Correlation Notes

The default correlation plot intentionally excludes identifier columns, `freq*`, and `sd_*` columns to keep the PDF readable. Use `--include-freq`, `--include-sd`, explicit `--metrics`, or `--metrics all` when those columns should be visualized.

Variables are ordered by Ward clustering on magnitude-based distance `1 - |corr|`, but no dendrogram is drawn. Triangle modes draw only the visible half with a stepped triangular border outside the diagonal cells. `--triangle upper` places metric labels on top/right and moves the colorbar to the left; `--triangle lower` places labels on left/bottom and keeps the colorbar on the right. Full matrices also keep the colorbar on the right.

## Warnings and Errors

| Situation | Behavior |
|---|---|
| Missing both `--msa-dir` and `--tree-dir` | Exit with error. |
| `--pseudo-tree-metrics` without `--msa-dir` | Exit with error. |
| `--outgroup-list` or `--ref-tree` without `--tree-dir` | Exit with error. |
| Non-empty output directory without `--overwrite` | Main `metrics` exits with error. |
| MSA/tree taxa mismatch | Warning is recorded; computation continues. |
| Unpaired MSA/tree files | Warning is recorded for each unmatched stem. |
| Marker with >200 taxa and pairwise identity enabled | Warning suggests `--skip-pairwise-identity`. |
| FastTree unavailable | Pseudo-tree metrics are empty; core metrics continue. |
| `--cluster-rectangles` with `--triangle lower/upper` | Warning is printed unless `--quiet`; rectangles are ignored. |

## Notes

- `metrics.csv` is the canonical intermediate consumed by `pretree filter`.
- `rcfv` is the classic relative composition frequency variability; `nrcfv` is the bias-corrected metric from Fleming and Struck.
- `saturation` is the slope of patristic distance versus uncorrected sequence distance through the origin.
- `average_pairwise_identity` is O(n² x L); use `--skip-pairwise-identity` for many-taxon datasets.
