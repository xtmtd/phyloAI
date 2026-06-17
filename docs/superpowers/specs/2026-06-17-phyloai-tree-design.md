# PhyloAI Tree Module Design Specification

**Date:** 2026-06-17  
**Status:** Draft for user review  
**Parent spec:** `2026-06-07-phyloai-design.md`

---

## 1. Purpose

`phyloai tree` is the phylogenetic inference layer of PhyloAI. It consumes prepared alignments or gene trees from `pretree` and produces tree-level inference results and support summaries.

The tree layer has four submodules:

- `ml` for maximum-likelihood inference with IQ-TREE and FastTree
- `bi` for Bayesian inference with PhyloBayes
- `msc` for multispecies coalescent inference with wASTRAL
- `concordance` for concordance-factor calculations such as gCF and sCF

The module boundary is deliberate:

- ML and BI are supermatrix methods that consume a single concatenated alignment or matrix.
- MSC is a supertree-style method that consumes a set of gene trees.
- Concordance factors are node-support summaries, not downstream post-tree analytics, so they belong in `tree`.

---

## 2. CLI Surface

```bash
phyloai tree ml iqtree --matrix ./concat/matrix.fa --mode pmsf
phyloai tree ml fasttree --matrix ./concat/matrix.fa
phyloai tree bi phylobayes --matrix ./concat/matrix.phy --chains 3
phyloai tree msc wastral --gene-trees ./genetrees/
phyloai tree concordance --tree ./tree.nwk --gene-trees ./genetrees/
```

Each subcommand is a distinct inference or support workflow. The top-level `tree` group exists to keep the CLI organized; it does not imply a shared execution model beyond shared output conventions from the main design.

---

## 3. Submodule Responsibilities

### 3.1 `tree ml`

`tree ml` owns maximum-likelihood inference from a single matrix. It provides two backends:

- `iqtree` for partitioned, unpartitioned, mixture, and model-rich ML analyses
- `fasttree` for quick approximate tree inference and lightweight pseudo-tree use cases

The submodule is responsible for tool selection, matrix input validation, output verification, and result reporting. It is not responsible for building the matrix itself.

### 3.2 `tree bi`

`tree bi` owns Bayesian tree inference with PhyloBayes. It is responsible for chain setup, chain-level tracking, convergence-related outputs, and tool-native resume integration.

### 3.3 `tree msc`

`tree msc` owns multispecies coalescent inference with wASTRAL. It consumes gene trees and produces a species tree plus any local branch support or summary outputs associated with the run.

### 3.4 `tree concordance`

`tree concordance` computes gCF/sCF-style branch support from a tree and its supporting gene trees. It is a support-summary module, not a topology-search module.

---

## 4. Output Conventions

All tree subcommands follow the shared PhyloAI output conventions:

- write `result.json` in the command output directory
- write a command log alongside `result.json`
- use the standard output-directory conflict policy and `--resume` only where long-running behavior justifies it
- preserve any tool-native files needed for auditability, but always include a PhyloAI result wrapper

Tree commands do not define a separate FASTA-writing policy unless they themselves emit sequence files for intermediate or summarized outputs.

---

## 5. Documentation Requirements

Implementation must update or add:

- `docs/commands/tree-ml.md`
- `docs/commands/tree-bi.md`
- `docs/commands/tree-msc.md`
- `docs/commands/tree-concordance.md`
- `docs/commands` and main design references that still mention `tree genetree`, `tree iqtree`, `tree astral`, or `tree phylobayes`

---

## 6. Relationship to Posttree

`posttree` keeps topology tests, dating, signal analysis, systematic error diagnostics, and simulation.

`concordance` moves into `tree` because concordance factors are branch-support measures attached to tree inference, not post hoc analysis of finished phylogenetic results.
