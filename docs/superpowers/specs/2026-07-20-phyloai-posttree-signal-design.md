# PhyloAI Posttree Signal Design Specification

**Date:** 2026-07-20
**Status:** Approved for implementation
**Parent spec:** `2026-06-07-phyloai-design.md`
**References:**
- Shen et al. 2017, Nature Ecology & Evolution (site/gene-wise lnL, ΔSLS/ΔGLS, outlier genes)
- Shen et al. 2021, Systematic Biology (consistent genes, GLS + GQS)
- IQ-TREE3 Command Reference: `--wslr`, `-lmap`, `-lmclust`

---

## 1. Purpose

`phyloai posttree signal` provides three phylogenetic signal distribution analyses, all implemented as independent subcommands under a shared group:

| Subcommand | Analysis | Core tool |
|------------|----------|-----------|
| `lnl` | Site-wise and gene-wise log-likelihood score distribution | IQ-TREE3 `-wslr` |
| `consistent` | Consistent gene identification via GLS + GQS | IQ-TREE3 `-wslr` + wASTRAL |
| `fclm` | Four-cluster Likelihood Mapping | IQ-TREE3 `-lmap -lmclust` |

All three share the same model parameter interface as `posttree topology` (`--model-expr`, `--partitions`, `--guide-tree`, `--tool-args`).

---

## 2. Design Principles

1. **Reuse topology infrastructure.** Model parameter handling, candidate-tree merging, IQ-TREE path resolution, and `--tool-args` validation are lifted directly from `posttree/topology.py`.
2. **Separate subcommands, not modes.** The three analyses have different required inputs, external tool dependencies, and output schemas. Merging them into one command would create parameter ghosts and unclear help text.
3. **`--partitions` vs `--locus-ranges` separation.** In `lnl` and `consistent`, `--partitions` is passed to IQ-TREE (`-p` or `-Q`) AND used for locus boundary extraction. `--locus-ranges` is a parallel input accepted in the same partition file format but used only for locus boundary extraction, not passed to IQ-TREE. The two flags are mutually exclusive.
4. **Gene-wise output is conditional.** Gene-wise tables, plots, and outlier analysis are only produced when locus boundaries are available via `--partitions` or `--locus-ranges`.
5. **ΔSLS/ΔGLS formulas are tree-count-aware.** For exactly 2 candidate trees: ΔSLS = lnL_T1 − lnL_T2 (signed, Shen 2017 eq. 1/2). For >2 trees: ΔSLS/ΔGLS = mean of all pairwise |lnL_Ta − lnL_Tb| (eq. 5/6). `support_sig` column (|ΔGLS| ≥ 2) only appears in gene-wise tables when exactly 2 trees are compared.
6. **`consistent` requires exactly 2 candidate trees.** GLS and GQS are both T1 vs T2 comparisons. >2 trees → hard error.

---

## 3. CLI Surface

```bash
# Site-wise and gene-wise lnL distribution
phyloai posttree signal lnl \
  --matrix matrix.aa.fa \
  --candidate-trees trees \
  [--model-expr LG+F+R4 | --partitions partitions.txt [--partition-mode p|Q]] \
  [--locus-ranges partitions.txt] \
  [--guide-tree guide.nwk] \
  [--tool-args "..."] \
  [--threads auto] \
  [--prefix lnl] \
  [--output-dir runs/posttree/signal/lnl] \
  [--metrics metrics.csv] \
  [--overwrite] [--dry-run] [--quiet] [--resume]

# Consistent gene identification
phyloai posttree signal consistent \
  --matrix matrix.aa.fa \
  --candidate-trees T1.tre,T2.tre \
  --tree-dir gene_trees/ \
  [--model-expr LG+F+R4 | --partitions partitions.txt [--partition-mode p|Q]] \
  [--locus-ranges partitions.txt] \
  [--guide-tree guide.nwk] \
  [--tool-args "..."] \
  [--threads 4] \
  [--prefix consistent] \
  [--output-dir runs/posttree/signal/consistent] \
  [--metrics metrics.csv] \
  [--overwrite] [--dry-run] [--quiet] [--resume]

# Four-cluster Likelihood Mapping
phyloai posttree signal fclm \
  --matrix matrix.aa.fa \
  --taxset-csv taxsets.csv \
  [--model-expr LG+C60+F+R4 | --partitions partitions.txt [--partition-mode p|Q]] \
  [--lmap ALL|<N>] \
  [--guide-tree guide.nwk] \
  [--tool-args "..."] \
  [--threads auto] \
  [--prefix fclm] \
  [--output-dir runs/posttree/signal/fclm] \
  [--overwrite] [--dry-run] [--quiet] [--resume]
```

### 3.1 Command Hierarchy

```
phyloai posttree
└── signal
    ├── lnl
    ├── consistent
    └── fclm
```

---

## 4. Parameters

### 4.1 Shared Parameters (all three subcommands)

| Flag | Type | Default | Description |
|------|------|---------|-------------|
| `--matrix` | Path | required | Single supermatrix alignment. Maps to IQ-TREE `-s`. |
| `--model-expr` | str | None | Complete IQ-TREE model expression (e.g. `LG+F+R4`, `C20+F+R4`). When combined with `--partitions`, the same model is applied independently to each partition. |
| `--partitions` | Path | None | Partition file passed to IQ-TREE as `-p` or `-Q` (per `--partition-mode`). In `lnl`/`consistent`, also extracts locus boundaries. When combined with `--model-expr`, each partition independently estimates parameters using that model. |
| `--partition-mode` | `p\|Q` | None | `p` = `-p` (edge-linked); `Q` = `-Q` (edge-unlinked). Defaults to `p` when `--partitions` is provided. Only valid with `--partitions`. |
| `--guide-tree` | Path | None | Guide tree for PMSF-style models. Maps to IQ-TREE `-ft`. |
| `--tool-args` | str | None | Extra IQ-TREE flags. Blocked: `-s`, `-z`, `-p`, `-Q`, `--prefix` (all three); `-wslr` (`lnl`/`consistent`); `-lmap`, `-lmclust`, `-n` (`fclm`). |
| `--threads` | str | `auto` | IQ-TREE `-T` value (integer or `auto`). For `consistent`, also controls wastral gene-level parallelism. |
| `--prefix` | str | subcmd name | IQ-TREE output prefix. Defaults: `lnl`, `consistent`, `fclm`. |
| `--resume` | flag | False | Resume incomplete IQ-TREE run from native checkpoint. |
| `--output-dir` | Path | see §6 | Output directory. |
| `--overwrite` | flag | False | Delete and recreate output directory. |
| `--dry-run` | flag | False | Print IQ-TREE command without executing. |
| `--quiet` | flag | False | Suppress terminal output except errors. |

**Note:** `--partitions` and `--partition-mode` are shared across all three subcommands. They map to IQ-TREE `-p`/`-Q` (model specification) in `lnl` and `consistent`, and additionally extract locus boundaries for gene-wise calculation in those two subcommands. In `fclm`, `--partitions` is used for model specification only (no locus boundary extraction).

### 4.2 `signal lnl` Additional Parameters

| Flag | Type | Default | Description |
|------|------|---------|-------------|
| `--candidate-trees` | str | required | Tree-list file or comma-separated individual NEWICK files. Maps to IQ-TREE `-z` after optional merge. Same format as `posttree topology`. |
| `--partitions` | Path | None | Partition file passed to IQ-TREE as `-p` or `-Q` (per `--partition-mode`). Also used to extract locus boundaries for gene-wise calculation. Mutually exclusive with `--locus-ranges`. When combined with `--model-expr`, the same model is applied to each partition. |
| `--partition-mode` | `p\|Q` | `p` | Controls how `--partitions` is passed to IQ-TREE: `p` = `-p` (edge-linked proportional, shared topology + rate multipliers per partition); `Q` = `-Q` (edge-unlinked, independent branch lengths per partition). Only valid when `--partitions` is provided. |
| `--locus-ranges` | Path | None | Partition file used only to define locus boundaries for gene-wise calculation. Not passed to IQ-TREE. Mutually exclusive with `--partitions`. |
| `--metrics` | Path | None | Metrics CSV from `phyloai pretree metrics`. When provided, generates (i) outlier vs non-outlier comparison, and (ii) pairwise comparisons of gene groups supporting different candidate trees. All loci in `gene_lnl.csv` must be present in this file. |

**Validation rules:**
- At least one model source required: `--model-expr`, `--partitions`, or `-m`/`-p` in `--tool-args`.
- `--partitions` and `--locus-ranges` mutually exclusive.

### 4.3 `signal consistent` Additional Parameters

| Flag | Type | Default | Description |
|------|------|---------|-------------|
| `--candidate-trees` | str | required | Exactly 2 candidate trees (tree-list file or comma-separated). >2 trees → hard error. |
| `--tree-dir` | Path | required | Directory of gene tree files for GQS calculation via wASTRAL. Logical locus name is resolved per the global shared file matching policy (§9.7 of parent design): suffix-agnostic, removing 1–2 dot segments as needed (e.g. `gene.fa.treefile` → `gene`, `gene.treefile` → `gene`). Ambiguous matches raise a hard error. |
| `--partitions` | Path | None | Partition file passed to IQ-TREE as `-p` or `-Q` (per `--partition-mode`). Also used to extract locus boundaries for GLS. Mutually exclusive with `--locus-ranges`. When combined with `--model-expr`, the same model is applied to each partition. |
| `--partition-mode` | `p\|Q` | `p` | Controls how `--partitions` is passed to IQ-TREE: `p` = `-p` (edge-linked proportional, shared topology + rate multipliers per partition); `Q` = `-Q` (edge-unlinked, independent branch lengths per partition). Only valid when `--partitions` is provided. |
| `--locus-ranges` | Path | None | Partition file used only to define locus boundaries for GLS calculation. Not passed to IQ-TREE. Mutually exclusive with `--partitions`. |
| `--metrics` | Path | None | Metrics CSV from `phyloai pretree metrics`. Generates consistent vs inconsistent gene comparison. |

**Validation rules:**
- Exactly 2 candidate trees required; error if 1 or >2.
- Must have at least one of `--partitions` or `--locus-ranges` (GLS is required for consistent gene identification).
- `--partitions` and `--locus-ranges` mutually exclusive.
- `--partition-mode` only valid when `--partitions` is provided.
- **Locus–gene tree matching (hard errors):** Every locus in `--partitions`/`--locus-ranges` must have a matching gene tree file in `--tree-dir` (per §9.7 global matching policy). Extra gene tree files in `--tree-dir` that have no matching partition locus are silently ignored.

### 4.4 `signal fclm` Additional Parameters

| Flag | Type | Default | Description |
|------|------|---------|-------------|
| `--taxset-csv` | Path | required | Two-column CSV (`taxon,taxset`) defining cluster membership. PhyloAI converts to NEXUS format. Maps to IQ-TREE `-lmclust`. |
| `--partitions` | Path | None | Partition file (e.g. `.best_model.nex` from IQ-TREE). Passed to IQ-TREE as `-p` or `-Q` (per `--partition-mode`) for the likelihood mapping model. When combined with `--model-expr`, the same model is applied to each partition. |
| `--partition-mode` | `p\|Q` | `p` | Controls how `--partitions` is passed to IQ-TREE: `p` = `-p`; `Q` = `-Q`. Only valid when `--partitions` is provided. |
| `--lmap` | str | None | Number of quartets for likelihood mapping. `ALL` = all quartets; integer = fixed count; unset = `50 * n_taxa`. Maps to IQ-TREE `-lmap`. |

**Validation rules for `--taxset-csv`:**
- All taxa in CSV must match taxa in `--matrix` exactly (no missing, no extra).
- `taxset` assignments must be mutually exclusive (each taxon in exactly one taxset).
- Minimum 4 taxsets required for FcLM (IQ-TREE requirement).

---

## 5. Computation

### 5.1 `signal lnl`

**IQ-TREE command:**
```
iqtree3 -s <matrix> -m <model>|-p|-Q <partitions> -z <candidate.trees> -wslr [-ft <guide>] --prefix <prefix> -T <threads> [tool_args]
```

**`.sitelh` parsing:**
- Line 1: `<n_trees> <n_sites>`
- Lines 2+: `TreeN  lnL_s1  lnL_s2 ... lnL_sN`

**`site_lnl.csv` construction:**
- Columns: `site, lnL_T1, lnL_T2, ..., ΔSLS, support`
- 2 trees: `ΔSLS = lnL_T1 − lnL_T2` (signed, eq. 1)
- >2 trees: `ΔSLS = mean of all pairwise |lnL_Ta − lnL_Tb|` (eq. 5)
- `support`: tree with highest lnL; ties → `ambiguous`
- Sorted by `ΔSLS` descending

**`gene_lnl.csv` construction (when locus boundaries available):**
- **Partition boundary parsing:** `_parse_partition_ranges` auto-detects RAxML
  (`LG, geneA = 1-235`) vs NEXUS (`charset geneA = 1-235;`) format from the
  first line. NEXUS `charpartition` (model assignment) lines are silently skipped;
  only `charset` lines are extracted. When both `--model-expr` and `--partitions`
  are given, `--model-expr` supplies the model formula while `--partitions`
  provides the boundary definitions — model info in a `.best_model.nex` file is
  ignored for boundary extraction.
- For each locus, sum site lnL values over `start:end` range per tree
- Columns: `locus, lnL_T1, lnL_T2, ..., ΔGLS, support [, support_sig]`
- 2 trees: `ΔGLS = lnL_T1 − lnL_T2` (signed, eq. 2); `support_sig` = `True` when `|ΔGLS| >= 2`
- >2 trees: `ΔGLS = mean of all pairwise |lnL_Ta − lnL_Tb|` (eq. 6); no `support_sig` column
- Sorted by `ΔGLS` descending

**Outlier gene identification (eq. 3/4):**
- Compute IQR of `|ΔGLS|` values across all loci
- Upper whisker = Q3 + 1.5 × IQR; lower whisker = Q1 − 1.5 × IQR
- Outliers: loci with `|ΔGLS|` > upper whisker or < lower whisker
- Write locus names to `outlier_genes.txt`

**`--metrics` comparison:**
- Validate all outlier loci present in metrics CSV
- Split loci into outlier vs non-outlier groups
- Per-metric: compute group means and Wilcoxon rank-sum p-value
- Write `outlier_comparison.csv` (columns: `metric, outlier_mean, outlier_n, nonoutlier_mean, nonoutlier_n, wilcoxon_p`)
- Write `outlier_comparison.pdf` (grid-layout boxplots with colored fill, alpha 0.6, significance annotations ***/**/* / p=0.xxx)

**`--metrics` support-group comparison (when >=2 trees, locus boundaries available):**
- Group loci by their `support` value (which tree they favor; ambiguous excluded)
- Write a single merged `support_comparison.csv` with per-group means, counts, and all pairwise Wilcoxon p-values
- Write a single merged `support_comparison.pdf` with side-by-side boxplots per metric, pairwise significance brackets between groups

**Support summary:**
- `support_summary_sites.csv`: tree-level site support counts (columns: `tree, n_sites`), including `ambiguous`
- `support_summary_genes.csv` [if boundaries]: tree-level gene support counts (columns: `tree, n_genes`), including `ambiguous`

### 5.2 `signal consistent`

**IQ-TREE command:** identical to `lnl` (same `-wslr` flag).

**GLS calculation:** identical to `lnl` gene-wise, always 2-tree mode (ΔGLS = lnL_T1 − lnL_T2).

**GQS calculation (parallelised over gene trees using `--threads`):**

For each gene tree file in `--tree-dir`:
1. Read gene tree taxa using `Bio.Phylo`
2. **Extra-taxa check:** if gene tree contains taxa not present in the reference trees → hard error (indicates matrix/gene tree mismatch)
3. Identify missing taxa: `ref_taxa − gene_tree_taxa`
4. Prune T1 and T2 reference trees using `Bio.Phylo.prune()` for each missing taxon → `T1_pruned`, `T2_pruned`
5. **Post-prune viability check:** if pruned reference tree has < 4 taxa → skip this locus, record `status: skipped` and `reason: "pruned_tree_too_small"` in `gqs.csv`; GLS support for this locus is also set to `ambiguous` for consistency determination
6. Write temporary `T1_pruned.nwk` and `T2_pruned.nwk`
7. Run `wastral -i <gene.tre> -C -c T1_pruned.nwk --mode 4`; extract `Score: XX` from log → `GQS_T1`
8. Repeat for T2 → `GQS_T2`
9. `ΔGQS = GQS_T1 − GQS_T2`; support comparison uses tolerance `|GQS_T1 − GQS_T2| < 1e-9` → `ambiguous`; otherwise T1 or T2

**Consistency determination:**
- A locus is **consistent** if `GLS.support == GQS.support` AND `support != "ambiguous"` AND `gqs.status != "skipped"`
- All other loci (including skipped GQS loci) are **inconsistent**

**`--metrics` comparison:** identical logic to `lnl`, but comparing consistent vs inconsistent gene groups. Output files: `consistent_comparison.csv`, `consistent_comparison.pdf`.

### 5.3 `signal fclm`

**CSV → NEXUS conversion:**
```nexus
#NEXUS
begin sets;
  taxset <taxset1> = <taxon1> <taxon2> ...;
  taxset <taxset2> = <taxon3> <taxon4> ...;
end;
```

**IQ-TREE command:**
```
iqtree3 -s <matrix> -m <model>|-p|-Q <partitions> -lmap <N|ALL> -lmclust cluster.nexus -n 0 -T <threads> --prefix <prefix> [-ft <guide>] [tool_args]
```

`-lmap` value resolution: `ALL` → `"ALL"`; user integer → pass through; unset → `50 * n_taxa` (count from matrix).

All likelihood mapping statistics are contained in the native `<prefix>.iqtree` report. No additional summary CSV is generated.

---

## 6. Output Structure

### 6.1 `signal lnl`

```
runs/posttree/signal/lnl/
├── result.json
├── candidate.trees              # merged (only when >1 tree file provided)
├── site_lnl.csv                 # site-wise table, ΔSLS-sorted
├── site_support.pdf             # site support distribution bar chart
├── support_summary_sites.csv    # tree-level site support counts
├── gene_lnl.csv                 # [if locus boundaries provided]
├── gene_support.pdf             # [if locus boundaries provided]
├── support_summary_genes.csv    # [if locus boundaries] tree-level gene support counts
├── outlier_genes.txt            # [if locus boundaries provided]
├── outlier_comparison.csv       # [if --metrics provided]
├── outlier_comparison.pdf       # [if --metrics provided]
├── support_comparison.csv       # [if --metrics + >=2 support groups]
├── support_comparison.pdf       # [if --metrics + >=2 support groups]
└── iqtree/
    ├── <prefix>.sitelh          # IQ-TREE raw site log-likelihoods
    ├── <prefix>.iqtree          # IQ-TREE native report
    └── <prefix>.log             # IQ-TREE log
```

### 6.2 `signal consistent`

```
runs/posttree/signal/consistent/
├── result.json
├── candidate.trees              # [if >1 tree file merged]
├── gls.csv                      # locus, lnL_T1, lnL_T2, ΔGLS, support, support_sig
├── gqs.csv                      # locus, GQS_T1, GQS_T2, ΔGQS, support, status, reason
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

### 6.3 `signal fclm`

```
runs/posttree/signal/fclm/
├── result.json
├── cluster.nexus                # generated from --taxset-csv
└── iqtree/
    ├── <prefix>.lmap.eps        # IQ-TREE likelihood mapping figure
    ├── <prefix>.iqtree          # IQ-TREE native report (contains all lmap statistics)
    └── <prefix>.log
```

Default `<prefix>` values: `lnl`, `consistent`, `fclm` respectively.

Default `--output-dir` values: `runs/posttree/signal/lnl`, `runs/posttree/signal/consistent`, `runs/posttree/signal/fclm`.

---

## 7. `result.json` Schema

All three subcommands use the **Single Pattern** (JSON standard §5.2).

### 7.1 `signal lnl`

```json
{
  "status": "success",
  "command": "phyloai posttree signal lnl --matrix /abs/matrix.aa.fa --candidate-trees trees --model-expr LG+F+R4 --locus-ranges /abs/partitions.txt --threads auto --output-dir /abs/runs/posttree/signal/lnl",
  "wall_time": 42.1,
  "tool_versions": {"iqtree3": "3.x.x"},
  "params": {
    "matrix": "/abs/matrix.aa.fa",
    "candidate_trees_raw": "trees",
    "model_expr": "LG+F+R4",
    "partitions": null,
    "locus_ranges": "/abs/partitions.txt",
    "guide_tree": null,
    "threads": "auto",
    "tool_args": null,
    "metrics": null,
    "output_dir": "/abs/runs/posttree/signal/lnl",
    "overwrite": false,
    "dry_run": false,
    "quiet": false
  },
  "key_results": {
    "n_trees": 2,
    "n_sites": 5604,
    "n_loci": 20,
    "n_outlier_genes": 3
  },
  "error": null,
  "data": {
    "cmd": ["iqtree3", "-s", "...", "-m", "LG+F+R4", "-z", "candidate.trees", "-wslr", "-T", "auto", "--prefix", "lnl"],
    "tool_stderr": "",
    "tool_log": "iqtree/lnl.log",
    "summary": {
      "n_trees": 2,
      "n_sites": 5604,
      "n_loci": 20,
      "n_outlier_genes": 3
    },
    "output_files": {
      "site_lnl": {"path": "/abs/.../site_lnl.csv", "description": "Site-wise lnL scores per tree, ΔSLS, support; sorted by ΔSLS descending"},
      "site_support_plot": {"path": "/abs/.../site_support.pdf", "description": "Site support distribution bar chart"},
      "support_summary_sites": {"path": "/abs/.../support_summary_sites.csv", "description": "Number of sites supporting each topology"},
      "gene_lnl": {"path": "/abs/.../gene_lnl.csv", "description": "Gene-wise lnL scores per tree, ΔGLS, support; sorted by ΔGLS descending"},
      "gene_support_plot": {"path": "/abs/.../gene_support.pdf", "description": "Gene support distribution bar chart"},
      "support_summary_genes": {"path": "/abs/.../support_summary_genes.csv", "description": "Number of genes supporting each topology"},
      "outlier_genes": {"path": "/abs/.../outlier_genes.txt", "description": "Loci with |ΔGLS| outside boxplot whiskers (Shen 2017 eq. 3/4)"},
      "outlier_comparison": {"path": "/abs/.../outlier_comparison.csv", "description": "Outlier vs non-outlier per-metric means and Wilcoxon p-values"},
      "outlier_comparison_plot": {"path": "/abs/.../outlier_comparison.pdf", "description": "Outlier vs non-outlier metric distribution boxplots"},
      "iqtree_report": {"path": "/abs/.../iqtree/lnl.iqtree", "description": "IQ-TREE native report"},
      "iqtree_sitelh": {"path": "/abs/.../iqtree/lnl.sitelh", "description": "IQ-TREE raw site log-likelihoods"},
      "iqtree_log": {"path": "/abs/.../iqtree/lnl.log", "description": "IQ-TREE console log"}
    }
  }
}
```

Conditional `output_files` keys (gene-wise: present only when locus boundaries provided; `outlier_comparison*`: present when `--metrics` is provided; `support_comparison`/`support_comparison_plot`: present when `--metrics` is provided and at least two non-ambiguous support groups exist) follow the additive pattern (§5.4): the CLI handler updates and rewrites `result.json` after generating optional files.

`n_loci` and `n_outlier_genes` in `key_results` and `summary` are absent when no locus boundaries are provided. `support_comparison_sig_metrics` is present only when at least one support-group pair has a metric with p < 0.05.

### 7.2 `signal consistent`

Same structure as `lnl` with the following differences:

- `params` adds: `tree_dir`, `partition_mode`
- `tool_versions` adds: `wastral`
- `key_results` replaces outlier fields with: `n_loci`, `n_consistent`, `n_inconsistent`, `n_gqs_skipped`
- `data.summary` adds: `wastral_n_gene_trees`, `wastral_threads_used`, `n_gqs_skipped`
- `data.output_files` replaces outlier keys with: `gls`, `gqs`, `consistent_genes`, `inconsistent_genes`, `gls_support_plot`, `gqs_support_plot`, `consistent_comparison`, `consistent_comparison_plot`

**`gqs.csv` column schema:** `locus, GQS_T1, GQS_T2, ΔGQS, support, status, reason`
- `status`: `"success"` or `"skipped"`
- `reason`: `null` for success; `"pruned_tree_too_small"` when post-prune reference tree has < 4 taxa
- `GQS_T1`, `GQS_T2`, `ΔGQS`: `null` for skipped rows
- `support`: `"ambiguous"` for skipped rows

### 7.3 `signal fclm`

- `params`: `matrix`, `taxset_csv`, `model_expr`, `lmap`, `guide_tree`, `threads`, `tool_args`, `output_dir`, `overwrite`, `dry_run`, `quiet`
- `tool_versions`: `iqtree3`
- `key_results`: `n_taxsets`, `n_quartets`
- `data.output_files`: `cluster_nexus`, `lmap_figure`, `iqtree_report`, `iqtree_log`

---

## 8. External Tool Dependencies

| Tool | Subcommands | Detection | Notes |
|------|-------------|-----------|-------|
| `iqtree3` | all three | Already registered in `doctor` | Reuse existing path resolution |
| `wastral` | `consistent` | Already registered in `doctor` (via `tree msc`/`tree cf`) | Reuse existing detection |
| `Bio.Phylo` | `consistent` | Built-in (biopython dependency) | Used for taxon pruning |

No new tools require `doctor` registration.

---

## 9. Associated Updates

### 9.1 Report Integration

Update `phyloai/report/collector.py`:
- Add `"signal"` to `_THIRD_LEVEL` with members `{"lnl", "consistent", "fclm"}`
- Add three entries to `STEP_ORDER`, replacing the single `"posttree.signal"` stub:
  - `"posttree.signal.lnl"`
  - `"posttree.signal.consistent"`
  - `"posttree.signal.fclm"`
- Add `generate_methods_posttree_signal_lnl`, `generate_methods_posttree_signal_consistent`, `generate_methods_posttree_signal_fclm` methods to `phyloai/report/templates.py`, following the pattern of `generate_methods_posttree_topology`.

### 9.2 MCP

Replace the existing `phyloai_posttree_signal` stub with three MCP tools:
- `phyloai_posttree_signal_lnl`
- `phyloai_posttree_signal_consistent`
- `phyloai_posttree_signal_fclm`

Parameter schema follows `phyloai_posttree_topology` conventions.

### 9.3 README

Replace the single `phyloai posttree signal` stub line with three example lines:

```bash
phyloai posttree signal lnl        --matrix ./matrix.fa --candidate-trees trees --model-expr LG+F+R4
phyloai posttree signal consistent --matrix ./matrix.fa --candidate-trees T1.tre,T2.tre --tree-dir ./gene_trees --model-expr LG+F+R4
phyloai posttree signal fclm       --matrix ./matrix.fa --taxset-csv taxsets.csv --model-expr LG+C60+F+R4
```

### 9.4 commands documentation

Add (English + Chinese):
- `docs/commands/posttree-signal.md`
- `docs/commands/posttree-signal.zh.md`

Structure follows existing `posttree-topology.md`: Purpose, Usage, Inputs, Outputs, Examples, Warnings/Errors, Notes. Covers all three subcommands in one document.

### 9.5 Skills

Update `skills/` PhyloAI workflow skill to describe the three `signal` subcommands: typical usage, input/output, and parameter mutual-exclusion rules (`--partitions` vs `--locus-ranges`).

---

## 10. Implementation Notes

- `signal lnl` and `signal consistent` share `.sitelh` parsing and gene-wise lnL summation logic. Extract to a shared helper (e.g., `_parse_sitelh`, `_sum_gene_lnl`) in `posttree/signal.py`.
- `signal consistent` GQS loop: use `concurrent.futures.ProcessPoolExecutor` with `max_workers=threads` (matching the pattern in `pretree/align.py`). Each worker handles one gene tree: prune → wastral → parse score. Skipped loci (post-prune < 4 taxa) are recorded in `gqs.csv` with `GQS_T1=None`, `GQS_T2=None`, `ΔGQS=None`, `support=ambiguous`, `status=skipped`.
- **Floating-point ambiguous threshold:** `support` in GQS is determined by `abs(GQS_T1 - GQS_T2) < 1e-9`; same tolerance applied to site-wise lnL ties in `support` column of `site_lnl.csv` and `gene_lnl.csv`. Document this constant and test it explicitly.
- Outlier/consistent comparison plots: reuse `pretree/filter/cluster.py` boxplot helper for style consistency.
- `--taxset-csv` validation (fclm): read matrix taxon list from the alignment header before running IQ-TREE, using the same format-detection helpers as `posttree/topology.py`.
