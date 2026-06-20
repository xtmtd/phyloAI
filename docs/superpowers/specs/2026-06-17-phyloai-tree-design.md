# PhyloAI Tree Module Design Specification

**Date:** 2026-06-17  
**Last updated:** 2026-06-19 (added ml/iqtree spec; updated 9.9 tool-args semantics)  
**Status:** Approved  
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
# ML: batch gene trees
phyloai tree ml fasttree --msa-dir ./trimmed/seqs --seq-type AA --model lg \
    --mode normal --boot 1000 --cat 20 --gamma --threads 8 -o runs/tree/ml/fasttree
# ML: single supermatrix
phyloai tree ml fasttree --matrix ./concat/matrix.fa --seq-type NT --model gtr \
    --mode slow --boot 1000 -o runs/tree/ml/fasttree
# ML: IQ-TREE3
phyloai tree ml iqtree --matrix ./concat/matrix.fa --model LG --rate-heterogeneity +G4

# Bayesian
phyloai tree bi phylobayes --matrix ./concat/matrix.phy --chains 3

# MSC
phyloai tree msc wastral --gene-trees ./genetrees/

# Concordance factors
phyloai tree concordance --tree ./tree.nwk --gene-trees ./genetrees/
```

### CLI Hierarchy

```
phyloai tree (click.Group)
├── ml (click.Group)          # Maximum-likelihood tree inference
│   ├── fasttree              # FastTree backend
│   └── iqtree                # IQ-TREE3 backend
├── bi (click.Group)          # Bayesian inference
│   └── phylobayes            # PhyloBayes backend
├── msc (click.Group)         # Multispecies coalescent
│   └── wastral               # wASTRAL backend
└── concordance               # Concordance factors (gCF/sCF)
```

Each subcommand is a distinct inference or support workflow. The top-level `tree` group exists to keep the CLI organized; it does not imply a shared execution model beyond shared output conventions from the main design.

---

## 3. Submodule Responsibilities

### 3.1 `tree ml`

`tree ml` owns maximum-likelihood inference. It provides two backends:

- `fasttree` for quick approximate tree inference (batch gene trees or single supermatrix)
- `iqtree` for partitioned, unpartitioned, mixture, and model-rich ML analyses

**Two input modes:**
- `--msa-dir`: directory of MSA files → batch parallel gene tree inference via `ProcessPoolExecutor`
- `--matrix`: single concatenated matrix file → single-tree supermatrix inference
- `--msa-dir` and `--matrix` are mutually exclusive

**Shared parameters** (applicable to both `fasttree` and `iqtree`): `--msa-dir`, `--matrix`, `--seq-type` (AA|NT|auto), `--model` (domain varies by seq-type and backend), `--mode` (backend-specific: FastTree `normal|fastest|slow`, IQ-TREE `normal|fast`), `--boot` (int; optional, omit for no support), `--output-dir` / `-o` (default `runs/tree/ml/<backend>`), `--threads` (controls batch parallelism or tool threads per backend policy), `--overwrite`, `--resume`, `--dry-run`, `--quiet` / `-q`, `--tool-args`.

Detailed specification for FastTree: `docs/superpowers/specs/2026-06-18-phyloai-tree-ml-fasttree-design.md`.

Detailed specification for IQ-TREE: `docs/superpowers/specs/2026-06-19-phyloai-tree-ml-iqtree-design.md`.

### 3.2 `tree bi`

`tree bi` owns Bayesian tree inference with PhyloBayes. It is responsible for chain setup, chain-level tracking, convergence-related outputs, and tool-native resume integration.

### 3.3 `tree msc`

`tree msc` owns multispecies coalescent inference with wASTRAL. It consumes gene trees and produces a species tree plus any local branch support or summary outputs associated with the run.

### 3.4 `tree concordance`

`tree concordance` computes gCF/sCF-style branch support from a tree and its supporting gene trees. It is a support-summary module, not a topology-search module.

---

## 4. Input Format Policy

All tree subcommands that consume alignment or matrix files should support both **FASTA** (`.fa`, `.fas`, `.fasta`, `.faa`, `.fna`) and **phylip-relaxed** (`.phy`, `.phylip`) input formats wherever the underlying tool can consume them natively.

- **NEXUS** (`.nex`, `.nxs`, `.nexus`) is not required at the tree layer. Subcommands that need it must document the requirement individually and refer users to `phyloai pretree convert` for pre-conversion.
- If a tool cannot read a supported format natively, the subcommand spec must document the limitation and may either auto-convert via `core/formats.py` or reject the input with a clear message.
- Specific tools may impose additional constraints (e.g., IQ-TREE expects PHYLIP with a specific header form). Those constraints are documented in the tool-specific subcommand spec, not here.

---

## 5. Output Conventions

All tree subcommands follow the shared PhyloAI output conventions:

- write `result.json` in the command output directory
- write a command log alongside `result.json`
- use the standard output-directory conflict policy and `--resume` only where long-running behavior justifies it
- preserve any tool-native files needed for auditability, but always include a PhyloAI result wrapper

Tree commands do not define a separate FASTA-writing policy unless they themselves emit sequence files for intermediate or summarized outputs.

---

## 6. Documentation Requirements

Implementation must update or add:

- `docs/commands/tree-ml.md`
- `docs/commands/tree-bi.md`
- `docs/commands/tree-msc.md`
- `docs/commands/tree-concordance.md`
- `docs/commands` and main design references that still mention `tree genetree`, `tree iqtree`, `tree astral`, or `tree phylobayes`

---

## 7. Relationship to Posttree

`posttree` keeps topology tests, dating, signal analysis, systematic error diagnostics, and simulation.

`concordance` moves into `tree` because concordance factors are branch-support measures attached to tree inference, not post hoc analysis of finished phylogenetic results.
