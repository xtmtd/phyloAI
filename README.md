# PhyloAI

An AI-native modular phylogenomics analysis platform.

## Installation

```bash
pip install -e .
```

PhyloAI bundles TAPER 1.0.0 (`correction_multi.jl`) and BMGE 1.12 (`BMGE.jar`) inside the Python package, so these two tools do not need separate installation. Other external tools, such as IQ-TREE, MAFFT, PAML, PhyloBayes, ASTRAL/ASTER, TreeShrink, Java, and Julia, should still be installed by the user according to the operating system and local environment.

## Quick Start

```bash
phyloai doctor
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

Generate the script once from an environment where `phyloai` is installed, save it to a persistent file, and add a `source ...` line to your shell configuration so completion is available automatically in every new terminal.

Do not run `phyloai completion ...` dynamically from `.bashrc`, `.zshrc`, or other shell startup files. Generate the file once and source the saved script instead.

Example for Bash:

```bash
mkdir -p ~/.config/phyloai/completion
phyloai completion bash > ~/.config/phyloai/completion/phyloai.bash
```

Add this line to `~/.bashrc` so every new Bash shell loads the saved completion script automatically:

```bash
source ~/.config/phyloai/completion/phyloai.bash
```

If you only run `source ~/.config/phyloai/completion/phyloai.bash` manually in the current terminal, completion only works for that shell session.

Example for Zsh:

```bash
mkdir -p ~/.config/phyloai/completion
phyloai completion zsh > ~/.config/phyloai/completion/phyloai.zsh
```

Add this line to `~/.zshrc` so every new Zsh shell loads the saved completion script automatically:

```bash
source ~/.config/phyloai/completion/phyloai.zsh
```

If you only run `source ~/.config/phyloai/completion/phyloai.zsh` manually in the current terminal, completion only works for that shell session.

Example for Fish:

```bash
mkdir -p ~/.config/fish/completions
phyloai completion fish > ~/.config/fish/completions/phyloai.fish
```

Fish loads completions from that directory automatically in new shells, so you do not need to add an extra `source` line to `config.fish`.

## Commands

`phyloai doctor` is the only command that supports `--output-format text|json`. Other commands write structured results to `result.json` inside their output directory and use Rich terminal output unless `--quiet` is set.

| Command | Purpose | Documentation |
|---------|---------|---------------|
| `phyloai doctor` | Inspect external tool availability. | [docs/commands/doctor.md](docs/commands/doctor.md) |
| `phyloai pretree convert` | Normalize and convert sequence files before downstream analysis. | [docs/commands/pretree-convert.md](docs/commands/pretree-convert.md) |
| `phyloai pretree stats`   | Inspect one sequence/alignment file or summarize a directory of files. | [docs/commands/pretree-stats.md](docs/commands/pretree-stats.md)     |
| `phyloai pretree align`   | Align sequences with MAFFT or MAGUS. | [docs/commands/pretree-align.md](docs/commands/pretree-align.md)     |
| `phyloai pretree trim`    | Batch-trim aligned MSAs with PhyloAI using trimAl, BMGE, or ClipKIT backends. | [docs/commands/pretree-trim.md](docs/commands/pretree-trim.md)       |
| `phyloai pretree metrics` | Compute MSA/tree metrics, distribution plots, and compact correlation heatmaps for marker evaluation. | [docs/commands/pretree-metrics.md](docs/commands/pretree-metrics.md) |
| `phyloai pretree filter`  | Marker-level filtering: TAPER error-site masking, TreeShrink taxon pruning, metric-rule filtering, symmetry test filtering, cluster-based exploration. | [docs/commands/pretree-filter.md](docs/commands/pretree-filter.md) |
| `phyloai pretree concat`  | Concatenate multiple MSAs into a supermatrix with occupancy filtering, recoding, codon variants, and outgroup reordering. | [docs/commands/pretree-concat.md](docs/commands/pretree-concat.md) |
| `phyloai tree ml fasttree` | Infer ML gene trees or supermatrix trees using FastTree. | [docs/commands/tree-ml-fasttree.md](docs/commands/tree-ml-fasttree.md) |
| `phyloai tree ml iqtree`   | Infer ML trees with IQ-TREE3: homogeneous, heterogeneous, partitioned, and ModelFinder workflows. | [docs/commands/tree-ml-iqtree.md](docs/commands/tree-ml-iqtree.md) |
| `phyloai tree msc`   | Multispecies coalescent species tree inference with wASTRAL. | [docs/commands/tree-msc.md](docs/commands/tree-msc.md) |
| `phyloai tree cf`    | Concordance factor computation: gCF, sCF, sCFl (IQ-TREE3) and qCF (wASTRAL). | [docs/commands/tree-cf.md](docs/commands/tree-cf.md) |
