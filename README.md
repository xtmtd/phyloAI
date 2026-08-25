# PhyloAI

[English](README.md) | [中文](README.zh.md)

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
3. **Build matrices or gene trees:** concatenate retained loci for supermatrix analyses, generate gene-jackknife pseudoreplicates when large matrices make downstream inference too costly, or infer per-locus gene trees for coalescent workflows.
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

Four-command simulation workflow: extract empirical per-locus parameters from
IQ-TREE reports, simulate a batch of alignments, then re-apply the original
gap mask and assess adequacy:

```bash
phyloai posttree simulate alisim params --iqtree-dir reports --tree-dir trees -o runs/params
phyloai posttree simulate alisim iqtree --model-params runs/params/params.tsv --strategy pdf --num-simulations 100 -o runs/sim
phyloai posttree simulate alisim transfergaps --original-msa markers/concat.aa.fa --simulated-dir runs/sim/MSAs -o runs/transfer
phyloai posttree simulate adequacy --original-msa markers/concat.aa.fa --simulated-dir runs/transfer -o runs/adequacy
```

Systematic-error branch-length screen across posterior/model trees (node map
identifies the same biological node across differing topologies):

```bash
# Across-taxon composition heterogeneity screen (alignment only)
phyloai posttree syserror taxcomp --matrix matrix.aa.fa --seq-type AA \
  -o runs/posttree/syserror/taxcomp

# Systematic-error branch-length screen across posterior/model trees
phyloai posttree syserror brlen --tree-dir posterior_trees --mode node-to-tip \
  --map nodes.map.txt --node1 Collembola -o runs/posttree/syserror/brlen

# Site-rate ranking/extraction sensitivity analysis
phyloai posttree syserror rate --iqtree-rate matrix.rate --matrix raw.fa \
  --subset slow --fraction 0.25,0.5,0.75 -o runs/posttree/syserror/rate

# Compositional-constraint diagnostic across two model analyses
phyloai posttree syserror cca --site-freq chain1.sitefreq \
  --site-lnl1 lnl_LG/site_lnl.csv --site-lnl2 lnl_C20/site_lnl.csv \
  --model1-name LG --model2-name C20 -o runs/posttree/syserror/cca
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
| `phyloai pretree concat`  | Concatenate multiple MSAs into a supermatrix with occupancy filtering, recoding, codon variants, outgroup reordering, and gene-jackknife pseudoreplicates for expensive downstream runs such as Bayesian or high-memory heterogeneous ML analyses. | [docs/commands/pretree-concat.md](docs/commands/pretree-concat.md) |
| `phyloai tree ml fasttree` | Infer ML gene trees or supermatrix trees using FastTree. | [docs/commands/tree-ml-fasttree.md](docs/commands/tree-ml-fasttree.md) |
| `phyloai tree ml iqtree`   | Infer ML trees with IQ-TREE3: homogeneous, heterogeneous, partitioned, ModelFinder, and custom exchangeability/site-profile workflows. | [docs/commands/tree-ml-iqtree.md](docs/commands/tree-ml-iqtree.md) |
| `phyloai tree bi pb`    | MCMC chain inference with PhyloBayes-MPI: multi-chain MCMC, real-time convergence monitoring, trace plots, and resume. | [docs/commands/tree-bi.md](docs/commands/tree-bi.md) |
| `phyloai tree bi bpcomp` | Topology convergence analysis with bpcomp using user-specified burn-in. | [docs/commands/tree-bi.md](docs/commands/tree-bi.md) |
| `phyloai tree bi tracecomp` | Parameter convergence analysis with tracecomp using user-specified burn-in. | [docs/commands/tree-bi.md](docs/commands/tree-bi.md) |
| `phyloai tree bi readpb` | Posterior analysis and predictive checks with readpb_mpi; `--mode ss,rr,r` also creates a posterior-parameterized PMSF simulation partition. | [docs/commands/tree-bi.md](docs/commands/tree-bi.md) |
| `phyloai tree msc`   | Multispecies coalescent species tree inference with wASTRAL. | [docs/commands/tree-msc.md](docs/commands/tree-msc.md) |
| `phyloai tree cf`    | Concordance factor computation: gCF, sCF, sCFl (IQ-TREE3) and qCF (wASTRAL). | [docs/commands/tree-cf.md](docs/commands/tree-cf.md) |
| `phyloai posttree topology` | Tree topology tests (AU / KH / SH / WKH / WSH / c-ELW) comparing candidate trees against a supermatrix. | [docs/commands/posttree-topology.md](docs/commands/posttree-topology.md) |
| `phyloai posttree dating`  | Bayesian molecular dating with MCMCtree: IQ-TREE Hessian computation + MCMC divergence time estimation with diagnostics. | [docs/commands/posttree-dating.md](docs/commands/posttree-dating.md) |
| `phyloai posttree signal lnl`        | Site-wise and gene-wise log-likelihood score distribution across candidate trees using IQ-TREE3 `-wslr`. Identifies outlier genes (Shen et al. 2017) and compares metrics across genes supporting different topologies when `--metrics` is provided. | [docs/commands/posttree-signal.md](docs/commands/posttree-signal.md) |
| `phyloai posttree signal consistent` | Consistent gene identification via GLS + GQS (Shen et al. 2021). Requires exactly 2 candidate trees and a gene tree directory. | [docs/commands/posttree-signal.md](docs/commands/posttree-signal.md) |
| `phyloai posttree signal fclm`       | Four-cluster Likelihood Mapping assessing phylogenetic signal supporting alternative hypotheses among taxon groups. | [docs/commands/posttree-signal.md](docs/commands/posttree-signal.md) |
| `phyloai posttree modelcompare iqtree` | Relative model comparison via IQ-TREE3 ModelFinder (BIC/AIC/AICc), including optional heterogeneous mixture model expansion via `-madd`. | [docs/commands/posttree-modelcompare.md](docs/commands/posttree-modelcompare.md) |
| `phyloai posttree modelcompare pb` | Relative model comparison via PhyloBayes LOO-CV / wAIC from `.sitelogl` site log-likelihood files (Lartillot 2023), pure Python. | [docs/commands/posttree-modelcompare.md](docs/commands/posttree-modelcompare.md) |
| `phyloai posttree simulate alisim` | IQ-TREE3 AliSim simulation preserving empirical dataset properties: `params` extracts per-locus parameters from IQ-TREE reports, `iqtree` simulates single or batch MSAs (complete/mixed/pdf strategies, resumable), `transfergaps` re-introduces the original gap mask onto one or many simulated MSAs. | [docs/commands/posttree-simulate-alisim.md](docs/commands/posttree-simulate-alisim.md) |
| `phyloai posttree simulate adequacy` | Compare observed PPA-DIV, PPA-CONV, PPA-VAR, and PPA-COMP statistics with simulated MSAs using a local pure-Python posterior predictive check. | [docs/commands/posttree-simulate-adequacy.md](docs/commands/posttree-simulate-adequacy.md) |
| `phyloai posttree syserror taxcomp` | Across-taxon composition screen using a Pearson common-composition X2 with nominal and Holm-adjusted per-taxon p-values plus observed PPA-COMP composition distances. Exploratory only; no taxon removal or recoding. | [docs/commands/posttree-syserror-taxcomp.md](docs/commands/posttree-syserror-taxcomp.md) |
| `phyloai posttree syserror brlen` | Extract branch-length statistics (total, terminal, internal, patristic, tip-to-tip, node-to-node, node-to-tip) from one tree, a directory, or Newick multi-tree files to diagnose rate heterogeneity across taxa; `label-nodes` generates labeled reference trees and map templates. | [docs/commands/posttree-syserror-brlen.md](docs/commands/posttree-syserror-brlen.md) |
| `phyloai posttree syserror rate` | Site-rate ranking/extraction sensitivity utility for IQ-TREE or PhyloBayes rates; optionally writes slow/fast alignment subsets. | [docs/commands/posttree-syserror-rate.md](docs/commands/posttree-syserror-rate.md) |
| `phyloai posttree syserror cca` | Composition-constraint diagnostic comparing site-wise topology preference across two model analyses. | [docs/commands/posttree-syserror-cca.md](docs/commands/posttree-syserror-cca.md) |
| `phyloai report`   | Generate a reproducible analysis report (JSON + self-contained HTML with embedded figures, sortable tables, and a draft Methods paragraph). Auto-generated methods text should be carefully verified before publication use. | [docs/commands/report.md](docs/commands/report.md) |
