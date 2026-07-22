# phyloai posttree signal

[English](posttree-signal.md) | [中文](posttree-signal.zh.md)


## Purpose

Performs phylogenetic signal distribution analysis through three independent
subcommands:

| Subcommand | Analysis | Core tool |
|------------|----------|-----------|
| `lnl` | Site-wise and gene-wise log-likelihood score distribution | IQ-TREE3 `-wslr` |
| `consistent` | Consistent gene identification via GLS + GQS (Shen et al. 2021) | IQ-TREE3 + wASTRAL |
| `fclm` | Four-cluster Likelihood Mapping | IQ-TREE3 `-lmap -lmclust` |

These analyses examine how phylogenetic signal is distributed across sites and
genes, identify outlier loci with disproportionate influence on topology, and
assess phylogenetic signal in contentious branches via likelihood mapping.

## Usage

```bash
# Site-wise lnL distribution (homogeneous model)
phyloai posttree signal lnl --matrix ./matrix.fa --candidate-trees trees --model-expr LG+F+R4

# Partitions model with gene-wise output
phyloai posttree signal lnl --matrix ./matrix.fa --candidate-trees trees --partitions partitions.txt

# Gene-wise lnL with locus ranges + outlier analysis
phyloai posttree signal lnl --matrix ./matrix.fa --candidate-trees trees --model-expr LG+F+R4 --locus-ranges partitions.txt --metrics metrics.csv

# Consistent gene identification (GLS + GQS)
phyloai posttree signal consistent --matrix ./matrix.fa --candidate-trees T1.tre,T2.tre --tree-dir ./gene_trees --model-expr LG+F+R4 --locus-ranges partitions.txt

# Four-cluster Likelihood Mapping
phyloai posttree signal fclm --matrix ./matrix.fa --taxset-csv taxsets.csv --model-expr LG+C60+F+R4

# Four-cluster Likelihood Mapping with partitions
phyloai posttree signal fclm --matrix ./matrix.fa --taxset-csv taxsets.csv --partitions matrix.best_model.nex
```

## signal lnl — Site-wise and Gene-wise lnL Distribution

### Purpose

Computes site-wise and gene-wise log-likelihood scores across multiple candidate
trees using IQ-TREE3's `-wslr`. Identifies genes with disproportionate
phylogenetic signal (ΔGLS) following Shen et al. (2017).

### Inputs

| Input | Description |
|-------|-------------|
| `--matrix` | Single supermatrix alignment (FASTA, PHYLIP, NEXUS). Required. Maps to IQ-TREE `-s`. |
| `--candidate-trees` | Tree-list file or comma-separated individual NEWICK files. Required. Maps to IQ-TREE `-z`. Same format as `posttree topology`. |
| `--model-expr` | Complete IQ-TREE `-m` expression (e.g. `LG+F+R4`, `C20+F+R4`). When combined with `--partitions`, the same model is applied independently to each partition. |
| `--partitions` | Partition file passed to IQ-TREE as `-p` or `-Q` (per `--partition-mode`). Also extracts locus boundaries for gene-wise calculation. Mutually exclusive with `--locus-ranges`. When combined with `--model-expr`, each partition independently estimates parameters using that model. |
| `--partition-mode` | `p` = `-p` (edge-linked proportional); `Q` = `-Q` (edge-unlinked). Default `p`. Only valid when `--partitions` is provided. |
| `--locus-ranges` | Partition file for locus boundary extraction only (not passed to IQ-TREE). Mutually exclusive with `--partitions`. |
| `--guide-tree` | Guide tree for PMSF-style models. Maps to IQ-TREE `-ft`. |
| `--metrics` | Metrics CSV from `phyloai pretree metrics` for outlier-vs-nonoutlier gene comparison. |
| `--threads` | IQ-TREE `-T` value (integer or `auto`, default `auto`). |
| `--tool-args` | Extra IQ-TREE flags. Blocked: `-s`, `-z`, `-wslr`, `--prefix`, `-p`, `-Q`. |
| `--prefix` | IQ-TREE output prefix (default: `lnl`). |
| `--iqtree-path` | Explicit path to iqtree3 executable. |
| `--resume` | Resume incomplete IQ-TREE run (native checkpoint). |
| `--output-dir` | Output directory (default: `runs/posttree/signal/lnl`). |
| `--overwrite` | Delete and recreate output directory before running. |
| `--dry-run` | Print the IQ-TREE command without executing. |
| `--quiet` | Suppress terminal output except errors. |

### Gene-wise ΔSLS/ΔGLS Formulas

- **2 candidate trees:** ΔSLS = lnL_T1 − lnL_T2 (signed); ΔGLS = same sum over gene
- **>2 trees:** ΔSLS/ΔGLS = mean of all pairwise |lnL_Ta − lnL_Tb|

`support_sig` column (|ΔGLS| ≥ 2) only appears in gene-wise tables with 2 trees.

### Outputs

```
runs/posttree/signal/lnl/
├── result.json
├── candidate.trees              # merged (only when >1 tree file provided)
├── site_lnl.csv                 # site-wise table, ΔSLS-sorted
├── site_support.pdf             # site support distribution bar chart
├── support_summary_sites.csv    # site counts supporting each topology
├── gene_lnl.csv                 # [if locus boundaries provided]
├── gene_support.pdf             # [if locus boundaries provided]
├── support_summary_genes.csv    # [if locus boundaries provided] gene counts per topology
├── outlier_genes.txt            # [if locus boundaries provided]
├── outlier_comparison.csv       # [if --metrics provided]
├── outlier_comparison.pdf       # [if --metrics provided]
└── iqtree/
    ├── <prefix>.sitelh          # IQ-TREE raw site log-likelihoods
    ├── <prefix>.iqtree          # IQ-TREE native report
    └── <prefix>.log             # IQ-TREE log
```

### Examples

```bash
# Site-wise only, no gene-wise breakdown
phyloai posttree signal lnl --matrix matrix.fa --candidate-trees trees --model-expr LG+F+R4

# With gene-wise output via locus ranges
phyloai posttree signal lnl --matrix matrix.fa --candidate-trees trees --model-expr LG+F+R4 --locus-ranges partitions.txt

# With outlier-vs-normal metrics comparison
phyloai posttree signal lnl --matrix matrix.fa --candidate-trees trees --model-expr LG+F+R4 --locus-ranges partitions.txt --metrics metrics.csv
```

---

## signal consistent — Consistent Gene Identification

### Purpose

Identifies genes where both likelihood-based (GLS) and quartet-based (GQS)
phylogenetic signal agree in supporting one of two candidate topologies.
Requires exactly 2 candidate trees. Uses IQ-TREE3 for GLS and wASTRAL for GQS.
Based on Shen et al. (2021).

### Inputs

| Input | Description |
|-------|-------------|
| `--matrix` | Single supermatrix alignment. Required. |
| `--candidate-trees` | Exactly 2 candidate trees (tree-list or comma-separated). Required. |
| `--tree-dir` | Directory of gene tree files for GQS calculation. Required. |
| `--model-expr` | IQ-TREE model expression. When combined with `--partitions`, the same model is applied independently to each partition. |
| `--partitions` | Partition file passed to IQ-TREE (`-p` or `-Q` per `--partition-mode`). Also extracts locus boundaries. |
| `--partition-mode` | `p` (edge-linked) or `Q` (edge-unlinked). Default `p`. Only with `--partitions`. |
| `--locus-ranges` | Partition file for locus boundary extraction only. Mutually exclusive with `--partitions`. |
| `--guide-tree` | Guide tree for PMSF models. |
| `--metrics` | Metrics CSV for consistent-vs-inconsistent gene comparison. |
| `--threads` | IQ-TREE `-T` (default `auto`). Also controls wASTRAL parallelism. |
| `--tool-args` | Extra IQ-TREE flags. Blocked: `-s`, `-z`, `-wslr`, `--prefix`, `-p`, `-Q`. |
| `--prefix` | IQ-TREE output prefix (default: `consistent`). |
| `--iqtree-path` | Explicit path to iqtree3 executable. |
| `--wastral-path` | Explicit path to wastral executable. |
| `--resume` | Resume incomplete IQ-TREE run (native checkpoint). |
| `--output-dir` | Output directory (default: `runs/posttree/signal/consistent`). |
| `--overwrite` | Delete and recreate output directory before running. |
| `--dry-run` | Print the IQ-TREE command without executing. |
| `--quiet` | Suppress terminal output except errors. |

Validation rules:
- Exactly 2 candidate trees required (hard error for 1 or >2).
- All loci from `--partitions`/`--locus-ranges` must have a matching gene
  tree file in `--tree-dir`. Extra gene tree files in `--tree-dir` are silently ignored.

Gene tree taxon handling:
- **Extra taxa** (taxon in gene tree but not in reference trees): hard error —
  indicates a mismatch between matrix and gene tree taxon sets.
- **Missing taxa** (taxon in reference trees but not in gene tree): reference
  trees are pruned to match the gene tree's taxon set before GQS computation.
- **Post-prune < 4 taxa**: if a pruned reference tree has fewer than 4 taxa,
  the locus is skipped for GQS (recorded in `gqs.csv` with `status: skipped`,
  `reason: pruned_tree_too_small`). Such loci are always considered inconsistent
  regardless of their GLS support.

### Outputs

```
runs/posttree/signal/consistent/
├── result.json
├── candidate.trees              # [if >1 tree file merged]
├── gls.csv                      # gene-wise lnL comparison
├── gqs.csv                      # gene-wise GQS comparison
├── consistent_genes.txt
├── inconsistent_genes.txt
├── gls_support.pdf
├── gqs_support.pdf
├── consistent_comparison.csv    # [if --metrics provided]
├── consistent_comparison.pdf    # [if --metrics provided]
└── iqtree/
    ├── <prefix>.sitelh
    ├── <prefix>.iqtree
    └── <prefix>.log
```

`gqs.csv` columns:

| Column | Description |
|--------|-------------|
| `locus` | Logical locus name matching partition file |
| `GQS_T1` | wASTRAL quartet score against Tree 1 (`null` if skipped) |
| `GQS_T2` | wASTRAL quartet score against Tree 2 (`null` if skipped) |
| `ΔGQS` | `GQS_T1 - GQS_T2` (`null` if skipped) |
| `support` | `"T1"`, `"T2"`, or `"ambiguous"` (tie within 1e-9 tolerance) |
| `status` | `"success"` or `"skipped"` |
| `reason` | `null` for success; `"pruned_tree_too_small"` when post-prune reference tree has fewer than 4 taxa |

A locus is **consistent** when `gls.support == gqs.support`, both are not
`"ambiguous"`, and `gqs.status` is not `"skipped"`. All other loci are
**inconsistent**.

### Example

```bash
phyloai posttree signal consistent \
  --matrix matrix.fa \
  --candidate-trees T1.tre,T2.tre \
  --tree-dir gene_trees/ \
  --model-expr LG+F+R4 \
  --locus-ranges partitions.txt
```

---

## signal fclm — Four-cluster Likelihood Mapping

### Purpose

Performs Four-cluster Likelihood Mapping (FcLM) to assess phylogenetic signal
supporting alternative hypotheses of relationship among four taxon clusters.
Uses IQ-TREE3's `-lmap` and `-lmclust` flags.

### Inputs

| Input | Description |
|-------|-------------|
| `--matrix` | Single supermatrix alignment. Required. |
| `--taxset-csv` | Two-column CSV (`taxon,taxset`) defining cluster membership. Minimum 4 taxsets. Required. |
| `--model-expr` | IQ-TREE model expression (e.g. `LG+C60+F+R4`). When combined with `--partitions`, the same model is applied independently to each partition. |
| `--partitions` | Partition file (e.g. `.best_model.nexus` from IQ-TREE). Passed to IQ-TREE as `-p` or `-Q` (per `--partition-mode`). When combined with `--model-expr`, each partition independently estimates parameters using that model. |
| `--partition-mode` | `p` = `-p` (edge-linked); `Q` = `-Q` (edge-unlinked). Default `p`. Only when `--partitions` is provided. |
| `--lmap` | Quartet count: `ALL` for all quartets, integer for fixed count, or omit for `50 × n_taxa`. Maps to IQ-TREE `-lmap`. |
| `--guide-tree` | Guide tree for PMSF models. |
| `--threads` | IQ-TREE `-T` (default `auto`). |
| `--tool-args` | Extra IQ-TREE flags. Blocked: `-s`, `-lmap`, `-lmclust`, `-n`, `-p`, `-Q`, `--prefix`. |
| `--prefix` | IQ-TREE output prefix (default: `fclm`). |
| `--iqtree-path` | Explicit path to iqtree3 executable. |
| `--resume` | Resume incomplete IQ-TREE run (native checkpoint). |
| `--output-dir` | Output directory (default: `runs/posttree/signal/fclm`). |
| `--overwrite` | Delete and recreate output directory before running. |
| `--dry-run` | Print the IQ-TREE command without executing. |
| `--quiet` | Suppress terminal output except errors. |

Validation rules:
- All taxa in CSV must match taxa in `--matrix` exactly.
- `taxset` assignments must be mutually exclusive (each taxon in exactly one taxset).
- Minimum 4 taxsets required (IQ-TREE FcLM requirement).

### Outputs

```
runs/posttree/signal/fclm/
├── result.json
├── cluster.nexus                # generated from --taxset-csv
└── iqtree/
    ├── <prefix>.lmap.eps        # IQ-TREE likelihood mapping figure
    ├── <prefix>.iqtree          # IQ-TREE native report (all lmap statistics)
    └── <prefix>.log
```

### Example

```bash
# Homogeneous model
phyloai posttree signal fclm --matrix matrix.fa --taxset-csv taxsets.csv --model-expr LG+C60+F+R4

# Partition model
phyloai posttree signal fclm --matrix matrix.fa --taxset-csv taxsets.csv --partitions matrix.best_model.nex
```

---

## Shared Notes

- All three subcommands are single-matrix only (no batch mode).
- `--model-expr` and `--partitions` can be combined: `--model-expr` specifies the
  substitution model and `--partitions` provides partition boundaries. Each
  partition independently estimates parameters using the same model formula.
  `--partition-mode` (default `p`) controls whether `--partitions` is passed to
  IQ-TREE as `-p` (edge-linked) or `-Q` (edge-unlinked).
- `--partitions` vs `--locus-ranges`: `--partitions` is passed to IQ-TREE AND
  used for locus boundaries; `--locus-ranges` is used only for boundaries.
- IQ-TREE output files (`.sitelh`, `.iqtree`, `.lmap.eps`, `.log`) are placed in
  an `iqtree/` subdirectory under the output directory. IQ-TREE stdout streams to
  the terminal during execution.
- Default `--output-dir` values: `runs/posttree/signal/lnl`,
  `runs/posttree/signal/consistent`, `runs/posttree/signal/fclm`.
- `--dry-run` prints the IQ-TREE command and validates inputs without running
  external tools. Still reads input files to validate partition ranges, tree
  matching, and taxset assignments; does NOT create output directories or write
  cluster.nexus.
- References: Shen et al. (2017) *Nature Ecology & Evolution*; Shen et al. (2021)
  *Systematic Biology*.

## Exit Codes

| Code | Meaning |
|------|---------|
| 0 | Success |
| 1 | User input error (missing files, invalid parameters, output conflict) |
| 2 | External tool execution failed |
| 3 | External tool executable not found |
