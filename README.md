# PhyloAI

An AI-native modular phylogenomics analysis platform.

PhyloAI connects common pre-tree, tree, and post-tree tasks into a documented
command-line workflow, while keeping each step inspectable through structured
`result.json` outputs and optional HTML reports. It is designed to turn marker
sequence folders into traceable phylogenomic evidence: every major step records
what was run, which tools were used, which loci or taxa were retained, and what
diagnostics support the final tree.

With its MCP server and guided workflow Skill, PhyloAI can also be driven
conversationally: an AI assistant can check the environment, choose the right
command, review parameters, inspect run status, read `result.json`, diagnose
failures, and help explain the results. Finished analyses can be summarized into
readable HTML reports with embedded figures, tables, provenance, and draft
Methods text.

## Why PhyloAI

Modern phylogenomics is no longer only about inferring one best tree. Practical
analyses need data cleaning, alignment, trimming, marker diagnostics, gene-tree
and species-tree inference, topology tests, dating, concordance analysis, and
publication-ready reporting. Existing tools often excel at one stage but leave
users to stitch together parameters, logs, file formats, and diagnostics by
hand.

PhyloAI focuses on that integration layer and makes it stronger:

- **One framework for many study designs:** run supermatrix or supertree analyses from the same marker directory, then extend into topology tests, concordance factors, Bayesian inference, molecular dating, and reports without rebuilding the workflow by hand.
- **Transparent rather than black-box:** every command writes parameters, detected tool versions, logs, decision tables, output paths, and summaries to predictable output directories.
- **Quality-control centered:** marker statistics, symmetry tests, TAPER masking, TreeShrink pruning, cluster-based marker exploration, and correlation plots make data problems visible before final inference.
- **Best-practice backends without lock-in:** PhyloAI orchestrates established tools such as MAFFT, trimAl, BMGE, ClipKIT, TAPER, IQ-TREE3, FastTree, wASTRAL, PhyloBayes-MPI, TreeShrink, and MCMCtree while preserving their native behavior and citations.
- **Recoverable and auditable runs:** structured `result.json` files make interrupted analyses, failed loci, retained/dropped decisions, and tool-version differences easy to inspect or rerun.
- **AI-native by design:** the MCP server and guided workflow Skill let AI assistants inspect schemas, check run status, read results, explain analyses, and support conversational phylogenomics without hiding command-line details.
- **Readable reports, not just files:** `phyloai report` turns run directories into self-contained HTML reports with embedded plots, sortable tables, provenance, and draft Methods text for manuscript preparation.

The input boundary is deliberate: PhyloAI does not assemble reads, call targets,
infer ortholog groups, or extract BUSCO/UCE-style markers. Those upstream steps
should be completed before PhyloAI. Once marker files are available, PhyloAI
provides a lightweight, scriptable framework for alignment, filtering,
supermatrix or supertree inference, diagnostics, and reporting.

## Workflow Overview

PhyloAI commands follow the common shape of a phylogenomic study:

1. **Prepare marker sequences:** convert formats, inspect sequence statistics, align unaligned loci, trim MSAs, and remove problematic sites or taxa.
2. **Evaluate markers:** compute occupancy, entropy, pairwise identity, compositional bias, saturation, tree distance, and correlation summaries.
3. **Build matrices or gene trees:** concatenate retained loci for supermatrix analyses, or infer per-locus gene trees for coalescent workflows.
4. **Infer species relationships:** run ML supermatrix trees, Bayesian analyses, or wASTRAL species-tree inference.
5. **Diagnose conflict and robustness:** compute concordance factors, run topology tests, compare marker clusters, and check posterior/prior behavior in dating analyses.
6. **Report the run:** collect parameters, tool versions, tables, figures, and draft Methods text into reproducible JSON and HTML reports.

Each command writes its own output directory. The important files are usually
`result.json`, logs, decision tables, and final sequence/tree/report outputs.
This makes failed or partial runs easier to inspect, resume, or explain with an
AI assistant.

## Installation

```bash
git clone https://github.com/xtmtd/phyloAI.git
cd phyloAI
pip install -e .
```

PhyloAI does not bundle third-party phylogenetics executables. Install the external tools needed for your workflow, then verify them with:

```bash
phyloai doctor
```

See [docs/commands/installation.md](docs/commands/installation.md) for Python environment options, external tool groups, and operating-system notes.

## Quick Start

```bash
phyloai doctor
```

One-click phylogenomics pipeline from raw sequences to a species tree:

```bash
phyloai run --seq-dir ./markers
phyloai run --seq-dir ./markers --mode supertree --speed fast --threads 16
```

Show all available commands:

```bash
phyloai --help
```

## Shell Completion

PhyloAI can generate static shell completion scripts for Bash, Zsh, and Fish:

```bash
phyloai completion bash
phyloai completion zsh
phyloai completion fish
```

Generate the script once and configure your shell to load the saved file. See [docs/commands/completion.md](docs/commands/completion.md) for Bash, Zsh, and Fish setup examples.

## AI Integration

PhyloAI includes an MCP server and a guided-workflow Skill for conversational analysis. See [docs/commands/ai-integration.md](docs/commands/ai-integration.md) for setup and usage with OpenCode, Claude Code, or Codex.

## License

PhyloAI-authored code is free to use, copy, modify, and distribute for academic,
educational, and non-commercial research purposes. Commercial use, commercial
redistribution, sublicensing, sale, or integration into commercial products or
services requires prior written permission from the copyright holder. See
[LICENSE](LICENSE).

This repository also interoperates with third-party software, and each
third-party component keeps its own license. Tools listed in
[docs/commands/installation.md](docs/commands/installation.md) are external
dependencies and must be installed and used under their upstream licenses. This
section is a project-level license notice, not a replacement for those
third-party licenses.

## Author And Contact

Feng ZHANG  
Nanjing Agricultural University  
Email: <xtmtd.zf@gmail.com>

## Commands

`phyloai doctor` is the only command that supports `--output-format text|json`. Other commands write structured results to `result.json` inside their output directory and use Rich terminal output unless `--quiet` is set.

| Command | Purpose | Documentation |
|---------|---------|---------------|
| `phyloai doctor` | Inspect external tool availability. | [docs/commands/doctor.md](docs/commands/doctor.md) |
| Installation | Set up Python environments and external tools, then verify with `phyloai doctor`. | [docs/commands/installation.md](docs/commands/installation.md) |
| `phyloai completion` | Generate static Bash, Zsh, or Fish shell completion scripts. | [docs/commands/completion.md](docs/commands/completion.md) |
| `phyloai run`     | One-click phylogenomics pipeline from raw sequences to a species tree. | [docs/commands/run.md](docs/commands/run.md) |
| `phyloai pretree convert` | Normalize and convert sequence files before downstream analysis. | [docs/commands/pretree-convert.md](docs/commands/pretree-convert.md) |
| `phyloai pretree stats`   | Inspect one sequence/alignment file or summarize a directory of files. | [docs/commands/pretree-stats.md](docs/commands/pretree-stats.md)     |
| `phyloai pretree align`   | Align sequences with MAFFT or MAGUS. | [docs/commands/pretree-align.md](docs/commands/pretree-align.md)     |
| `phyloai pretree trim`    | Batch-trim aligned MSAs with PhyloAI using trimAl, BMGE, or ClipKIT backends. | [docs/commands/pretree-trim.md](docs/commands/pretree-trim.md)       |
| `phyloai pretree metrics` | Compute MSA/tree metrics, distribution plots, and compact correlation heatmaps for marker evaluation. | [docs/commands/pretree-metrics.md](docs/commands/pretree-metrics.md) |
| `phyloai pretree filter`  | Marker-level filtering: TAPER error-site masking, TreeShrink taxon pruning, metric-rule filtering, symmetry test filtering, cluster-based exploration. | [docs/commands/pretree-filter.md](docs/commands/pretree-filter.md) |
| `phyloai pretree concat`  | Concatenate multiple MSAs into a supermatrix with occupancy filtering, recoding, codon variants, and outgroup reordering. | [docs/commands/pretree-concat.md](docs/commands/pretree-concat.md) |
| `phyloai tree ml fasttree` | Infer ML gene trees or supermatrix trees using FastTree. | [docs/commands/tree-ml-fasttree.md](docs/commands/tree-ml-fasttree.md) |
| `phyloai tree ml iqtree`   | Infer ML trees with IQ-TREE3: homogeneous, heterogeneous, partitioned, and ModelFinder workflows. | [docs/commands/tree-ml-iqtree.md](docs/commands/tree-ml-iqtree.md) |
| `phyloai tree bi`    | Bayesian phylogenetic inference with PhyloBayes-MPI: multi-chain MCMC, real-time convergence monitoring, trace plots, and resume. | [docs/commands/tree-bi.md](docs/commands/tree-bi.md) |
| `phyloai tree msc`   | Multispecies coalescent species tree inference with wASTRAL. | [docs/commands/tree-msc.md](docs/commands/tree-msc.md) |
| `phyloai tree cf`    | Concordance factor computation: gCF, sCF, sCFl (IQ-TREE3) and qCF (wASTRAL). | [docs/commands/tree-cf.md](docs/commands/tree-cf.md) |
| `phyloai posttree topology` | Tree topology tests (AU / KH / SH / WKH / WSH / c-ELW) comparing candidate trees against a supermatrix. | [docs/commands/posttree-topology.md](docs/commands/posttree-topology.md) |
| `phyloai posttree dating`  | Bayesian molecular dating with MCMCtree: IQ-TREE Hessian computation + MCMC divergence time estimation with diagnostics. | [docs/commands/posttree-dating.md](docs/commands/posttree-dating.md) |
| `phyloai report`   | Generate a reproducible analysis report (JSON + self-contained HTML with embedded figures, sortable tables, and a draft Methods paragraph). Auto-generated methods text should be carefully verified before publication use. | [docs/commands/report.md](docs/commands/report.md) |
