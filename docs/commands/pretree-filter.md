# `phyloai pretree filter`

Marker-level filtering with four subcommands for distinct filtering workflows.

## Subcommands

| Subcommand | Purpose |
|------------|---------|
| `taper` | Mask erroneous sequence stretches within MSAs using TAPER. |
| `treeshrink` | Remove outlier taxa from gene trees (and optionally MSAs) using TreeShrink. |
| `metrics` | Remove whole loci by explicit conditions on a `metrics.csv`-like table. |
| `cluster` | Group loci by metric profiles using PCA or UMAP + hierarchical clustering, with optional outlier cluster removal. |

## `filter taper` — TAPER Error-Site Masking

### Purpose

Run TAPER (bundled `correction_multi.jl`) to mask error sites in MSAs. Operates in three modes: AA-only, NT-only, AA+CDS (mask AA then project to codon-aligned NT).

### Usage

```bash
phyloai pretree filter taper --msa-dir <dir> [--nt-dir <dir>] [--seq-type AA|NT|auto]
```

### Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--msa-dir` | Path | required | Directory of MSA files. |
| `--nt-dir` | Path | — | Codon-aligned NT MSA directory for AA+CDS mode. |
| `--seq-type` | AA\|NT\|auto | auto | Molecule type. |
| `--cutoff` | int | 3 | TAPER `-c` cutoff; lower = more aggressive. |
| `--taper-path` | Path | — | Path to `correction_multi.jl`. |
| `--julia-path` | Path | — | Julia executable path. |
| `--tool-args` | str | — | Extra TAPER args (not -m,-a,-c,-l). |
| `-t, --threads` | int | 4 | Worker count. |
| `--resume` | flag | off | Resume from checkpoint. |

### Inputs

AA-only: AA MSA files in `--msa-dir`. Output masked AA to `seqs/`.
NT-only: NT MSA files + `--seq-type NT`. Output masked NT to `seqs/`.
AA+CDS: AA MSA files + `--nt-dir` with codon-aligned NT MSAs. Output masked AA to `seqs/faa/`, projected CDS to `seqs/fna/`.

### Outputs

```
runs/pretree/filter/taper/
  seqs/                         (or seqs/faa/ + seqs/fna/ for AA+CDS)
  retained_loci.csv|tsv
  dropped_loci.csv|tsv
  filter_decisions.csv|tsv
  filter.log
  result.json
  checkpoint.json               (internal; only with --resume)
```

### Examples

```bash
phyloai pretree filter taper --msa-dir ./trimmed
phyloai pretree filter taper --msa-dir ./trimmed --cutoff 5 --threads 8
phyloai pretree filter taper --msa-dir ./trimmed_aa --nt-dir ./trimmed_fna
```

---

## `filter treeshrink` — TreeShrink Taxon Pruning

### Purpose

Run TreeShrink to detect and remove outlier long-branch taxa from gene trees. Optionally prune matching MSAs.

### Usage

```bash
phyloai pretree filter treeshrink --tree-dir <dir> [--msa-dir <dir>]
```

### Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--tree-dir` | Path | required | Directory of gene tree files. |
| `--msa-dir` | Path | — | Optional MSA directory to also shrink. |
| `--threshold` | float | 0.05 | TreeShrink `-q` threshold. |
| `--treeshrink-mode` | Choice | auto | auto\|per-gene\|all-genes\|per-species. |
| `--keep-work-dir` | flag | off | Retain TreeShrink work directory for debugging. |

### Outputs

```
runs/pretree/filter/treeshrink/
  trees/
  seqs/                         (only when --msa-dir provided)
  retained_loci.csv|tsv
  modified_loci.csv|tsv
  dropped_loci.csv|tsv
  removed_taxa.csv|tsv
  filter_decisions.csv|tsv
  filter.log
  result.json
```

### Examples

```bash
phyloai pretree filter treeshrink --tree-dir ./genetrees
phyloai pretree filter treeshrink --tree-dir ./genetrees --msa-dir ./trimmed --threshold 0.1
```

---

## `filter metrics` — Metric Rule Filtering

### Purpose

Filter loci by explicit numeric or string conditions on a metrics CSV/TSV table. Conditions are AND-only.

### Usage

```bash
phyloai pretree filter metrics --table <file> --keep "col>=val,col<=val,..."
```

### Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--table` | Path | required | Metrics CSV/TSV. |
| `--keep` | str | required | Comma-separated conditions (AND-only). |
| `--loci-column` | str | loci | Column name for locus identifier. |
| `--copy` | flag | off | Copy retained MSAs/trees to output dir. |

### Rule Syntax

Supported operators: `>=`, `>`, `<=`, `<`, `==`, `!=`.

```bash
--keep "dvmc>=0,dvmc<=0.3,average_BS>=0.8"
--keep "DataType==AA,num_sites>=300"
```

### Outputs

```
runs/pretree/filter/metrics/
  retained_loci.csv|tsv
  dropped_loci.csv|tsv
  filter_decisions.csv|tsv
  seqs/                         (only with --copy --msa-dir)
  trees/                        (only with --copy --tree-dir)
  filter.log
  result.json
```

### Examples

```bash
phyloai pretree filter metrics --table ./metrics/metrics.csv --keep "dvmc<=0.3,average_BS>=0.8"
phyloai pretree filter metrics --table ./metrics.csv --keep "num_sites>=300" --copy --msa-dir ./trimmed
```

---

## `filter cluster` — Cluster-Based Exploration

### Purpose

Group loci by metric profiles using PCA (default) or UMAP dimensionality reduction + hierarchical clustering. Optionally remove outlier clusters.

### Usage

```bash
phyloai pretree filter cluster --table <file> [--reduction pca|umap]
```

### Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--table` | Path | required | Metrics CSV/TSV. |
| `--reduction` | pca\|umap | pca | Dimensionality reduction method. |
| `--n-clusters` | int | auto | Fixed cluster count. |
| `--cluster-linkage` | ward\|average\|complete\|single | ward | Linkage method. |
| `--drop-outlier-clusters` | none\|auto | none | Enable automatic outlier cluster removal. |
| `--outlier-metric` | str | average_BS | Metric for ranking clusters. |
| `--max-drop-fraction` | float | 0.2 | Max fraction of loci to drop. |
| `--copy` | flag | off | Copy retained MSAs/trees (only with auto drop). |

### Outputs

```
runs/pretree/filter/cluster/
  features_used.csv|tsv
  reduction.csv|tsv
  cluster_selection.csv|tsv
  clusters.csv|tsv
  cluster_summary.csv|tsv
  cluster_metric_means.csv|tsv
  cluster_metric_heatmap.pdf
  cluster_2d.pdf
  cluster_3d.pdf
  cluster_metric_boxplots_001.pdf ...
  cluster_loci/cluster_*.csv|tsv
  retained_loci.csv|tsv              (only with auto drop)
  dropped_loci.csv|tsv               (only with auto drop)
  filter_decisions.csv|tsv           (only with auto drop)
  outlier_comparison.csv|tsv         (only with auto drop)
  filter.log
  result.json
```

### Examples

```bash
phyloai pretree filter cluster --table ./metrics/metrics.csv
phyloai pretree filter cluster --table ./metrics.csv --reduction umap --n-clusters 5
phyloai pretree filter cluster --table ./metrics.csv --drop-outlier-clusters auto --copy --msa-dir ./trimmed
```

---

## Notes

- `cluster` is exploratory by default; it does not remove loci unless `--drop-outlier-clusters auto` is set.
- `taper` is the only filter subcommand supporting `--resume` (checkpoint-based recovery).
- File matching across `--msa-dir` and `--tree-dir` uses logical locus names, not filename suffix whitelists.
- All subcommands write `result.json` and `filter.log` to their output directory.
