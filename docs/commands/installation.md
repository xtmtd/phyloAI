# PhyloAI Installation Guide

## Purpose

This guide explains how to install PhyloAI, make external tools visible to the active shell environment, and verify the setup with `phyloai doctor`.

PhyloAI does not install third-party phylogenetics tools automatically. Install those tools through your operating system, Conda/Mamba environment, cluster module system, or the upstream project instructions.

## Get The Source

```bash
git clone https://github.com/xtmtd/phyloAI.git
cd phyloAI
```

## Python Environment

Choose one environment style.

### uv

Recommended for local development and quick reproducibility.

```bash
uv venv
source .venv/bin/activate
uv pip install -e .
```

### Conda / Mamba

Recommended when Python packages and bioinformatics command-line tools need to live in the same environment.

```bash
mamba create -n phyloai python=3.11
mamba activate phyloai
pip install -e .
```

Use `conda` instead of `mamba` if that is what your system provides.

### venv

Suitable for a pure Python environment when external tools are installed elsewhere on `PATH`.

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -U pip
python -m pip install -e .
```

## Verify The Environment

```bash
phyloai doctor
phyloai doctor --output-format json
```

`phyloai doctor` checks the current shell environment. Activate the environment you intend to use before running it. It does not install tools or modify `PATH`.

## One-Click Pipeline Dependencies

| Command mode | Main external tools |
|--------------|---------------------|
| `phyloai run --mode supermatrix --speed normal` | MAFFT, trimAl, Julia for TAPER, IQ-TREE3 |
| `phyloai run --mode supermatrix --speed fast` | MAFFT, trimAl, FastTree |
| `phyloai run --mode supertree --speed normal` | MAFFT, trimAl, Julia for TAPER, IQ-TREE3 for gene trees, wASTRAL |
| `phyloai run --mode supertree --speed fast` | MAFFT, trimAl, FastTree, wASTRAL |

## External Tools

### Core Workflow

| Tool | Used by | Install entry point | Detection name | Verify |
|------|---------|---------------------|----------------|--------|
| IQ-TREE3 | `tree ml iqtree`, topology tests, dating Hessian, normal supermatrix `run` | https://github.com/iqtree/iqtree3/releases | `iqtree3` | `phyloai doctor` |
| MAFFT | `pretree align`, `phyloai run` | https://mafft.cbrc.jp/alignment/software/ | `mafft` | `phyloai doctor` |
| trimAl | `pretree trim`, backtranslation, `phyloai run` | https://github.com/inab/trimal | `trimal` | `phyloai doctor` |

Treat these as external tools. PhyloAI checks whether they are visible to the active shell, but does not redistribute them.

### Tree And Posttree Tools

| Tool | Used by | Install entry point | Detection name | Verify |
|------|---------|---------------------|----------------|--------|
| FastTree | `tree ml fasttree`, fast `phyloai run` | http://www.microbesonline.org/fasttree/ | `FastTree` | `phyloai doctor` |
| wASTRAL | `tree msc`, supertree `phyloai run` | https://github.com/chaoszhang/ASTER | `wastral` | `phyloai doctor` |
| MCMCtree / PAML | `posttree dating mcmc` | https://github.com/abacus-gene/paml/releases | `mcmctree` | `phyloai doctor` |

### Bayesian Inference

`phyloai tree bi` needs the PhyloBayes-MPI tool group.

| Tool | Purpose | Install entry point | Detection name | Verify |
|------|---------|---------------------|----------------|--------|
| pb_mpi | MCMC sampler | https://github.com/bayesiancook/pbmpi | `pb_mpi` | `phyloai doctor` |
| bpcomp | Topology convergence | https://github.com/bayesiancook/pbmpi | `bpcomp` | `phyloai doctor` |
| tracecomp | Parameter convergence | https://github.com/bayesiancook/pbmpi | `tracecomp` | `phyloai doctor` |
| mpirun | MPI launcher | https://www.open-mpi.org/ | `mpirun` | `phyloai doctor` |
| readpb_mpi | Optional chain reader | https://github.com/bayesiancook/pbmpi | `readpb_mpi` | `phyloai doctor` |

If the tools are installed outside `PATH`, use `phyloai tree bi --pb-path /path/to/pbmpi/bin`.

### Filtering And Trimming Extras

| Tool | Used by | Install entry point | Detection name | Verify |
|------|---------|---------------------|----------------|--------|
| TreeShrink | `pretree filter treeshrink` | https://github.com/uym2/TreeShrink | `run_treeshrink.py` | `phyloai doctor` |
| MAGUS | `pretree align --method magus` | https://github.com/vlasmirnov/MAGUS | `magus` | `phyloai doctor` |
| ClipKIT | `pretree trim --tool clipkit` | https://github.com/JLSteenwyk/ClipKIT | `clipkit` | `phyloai doctor` |
| BMGE | `pretree trim --tool bmge` | https://github.com/BMGE/BMGE or upstream BMGE distribution | `BMGE.jar` | `phyloai doctor` |
| TAPER | `pretree filter taper` | upstream TAPER distribution | `correction_multi.jl` | `phyloai doctor` |

### Runtime Dependencies

| Tool | Used by | Install entry point | Detection name | Verify |
|------|---------|---------------------|----------------|--------|
| Java | BMGE workflows | https://www.java.com/ | `java` | `phyloai doctor` |
| Julia | TAPER masking | https://julialang.org/downloads/ | `julia` | `phyloai doctor` |

If BMGE or TAPER are installed outside `PATH`, pass their file paths with `--bmge-path /path/to/BMGE.jar` or `--taper-path /path/to/correction_multi.jl`.

## Operating System Notes

- macOS: Homebrew or Conda/Mamba are usually the simplest way to make command-line tools visible on `PATH`.
- Linux: use Conda/Mamba, distribution packages, cluster modules, or upstream binaries depending on your environment.
- WSL: install tools inside the Linux distribution, not only on Windows, so `phyloai doctor` can see them from the WSL shell.
