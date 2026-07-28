# phyloai tree bi

[English](tree-bi.md) | [中文](tree-bi.zh.md)


Bayesian phylogenetic inference with [PhyloBayes-MPI](https://github.com/bayesiancook/pbmpi).

## Overview

`phyloai tree bi` is a Click Group with four subcommands:

| Subcommand | Purpose |
|---|---|
| `tree bi pb` | MCMC chain inference with `pb_mpi` |
| `tree bi bpcomp` | Topology convergence analysis with `bpcomp` |
| `tree bi tracecomp` | Parameter convergence analysis with `tracecomp` |
| `tree bi readpb` | Posterior analysis and predictive checks with `readpb_mpi` |

The default output root is `runs/tree/bi/`.

## Requirements

Each subcommand resolves only its required PhyloBayes-MPI tools from `PATH` (or from the directory passed with `--pb-path`):

| Tool | Required by | Purpose |
|------|-------------|---------|
| `pb_mpi` | `bi pb` | MCMC sampler |
| `bpcomp` | `bi pb`, `bi bpcomp` | Topology convergence |
| `tracecomp` | `bi pb`, `bi tracecomp` | Parameter convergence |
| `mpirun` | `bi pb`, `bi readpb` | Open MPI launcher |
| `readpb_mpi` | `bi readpb` | Posterior analysis |

Run `phyloai doctor` to confirm installation.

---

# phyloai tree bi pb

Run N independent MCMC chains in parallel, monitor convergence in real time, and produce a consensus tree.

## Usage

```bash
phyloai tree bi pb --matrix <alignment> [OPTIONS]
```

## Examples

```bash
# Default: 3 chains, CAT-GTR model, run forever
phyloai tree bi pb --matrix concat/matrix.phy

# Homogeneous LG+G4, stop after 10000 total MCMC cycles
phyloai tree bi pb --matrix concat/matrix.phy --model lg --mixture 1 --nsamples 10000

# Add two extra chains to an existing run
phyloai tree bi pb --matrix concat/matrix.phy --chain-names chain4,chain5 -o runs/tree/bi

# Resume all chains
phyloai tree bi pb -o runs/tree/bi --resume

# Resume selected chains
phyloai tree bi pb -o runs/tree/bi --resume chain1,chain3

# Resume and extend to a new target
phyloai tree bi pb -o runs/tree/bi --resume --nsamples 10000

# Custom PhyloBayes tool directory
phyloai tree bi pb --matrix concat/matrix.phy --pb-path /opt/pbmpi/bin

# Dry-run
phyloai tree bi pb --matrix concat/matrix.phy --dry-run
```

## Parameters

### Input / Output

| Flag | Type | Default | Description |
|------|------|---------|-------------|
| `--matrix / -m` | Path | required | Input alignment (PHYLIP or FASTA). Not required when `--resume` is used. |
| `--output-dir / -o` | Path | `runs/tree/bi` | Output directory. |
| `--overwrite` | flag | False | Delete and recreate output directory. Mutually exclusive with `--resume`. |

### Model

| Flag | Choice | Default | pb_mpi flag | Description |
|------|--------|---------|-------------|-------------|
| `--model` | `gtr`, `poisson`, `lg`, `wag`, `jtt`, `mtrev`, `mtzoa`, `mtart` | `gtr` | `-gtr`, `-poisson`, … | Rate matrix. |
| `--mixture` | `auto`, `1`, or integer N | `auto` | `-cat` / no mixture flag / `-ncat N` | `auto` = CAT Dirichlet process; `1` = homogeneous single matrix; integer N = fixed N-component mixture. |
| `--gamma-cats` | int ≥ 1 | 4 | `-dgam N` | Discrete Gamma rate categories. |
| `--start-tree` | Path | None | `-t <file>` | Starting tree. Mutually exclusive with `--fix-tree`. |
| `--fix-tree` | Path | None | `-T <file>` | Fixed topology. Mutually exclusive with `--start-tree`. |

### Chains & Parallelism

| Flag | Type | Default | Description |
|------|------|---------|-------------|
| `--chains` | int ≥ 1 | 3 | Number of independent chains. |
| `--chain-prefix` | str | `chain` | Prefix for auto-named chains. |
| `--chain-names` | str | None | Comma-separated names. Overrides `--chains`/`--chain-prefix`. |
| `--threads / -t` | int ≥ 2 | 4 | MPI processes per chain (`mpirun -np`). |

### Sampling

| Flag | Type | Default | Description |
|------|------|---------|-------------|
| `--sample-freq` | int ≥ 1 | 1 | Save one point every N cycles. |
| `--nsamples` | int | `-1` | Total MCMC cycles per chain before stopping. `-1` = run forever. With `--sample-freq N`, saved points = cycles / N. |

### Convergence Monitoring

| Flag | Type | Default | Description |
|------|------|---------|-------------|
| `--monitor-freq` | int ≥ 1 | 100 | Run bpcomp + tracecomp every N new samples. |
| `--burnin-frac` | float `[0.0, 1.0)` | 0.5 | Fraction discarded as burn-in during monitoring. |
| `--poll-interval` | int ≥ 1 | 60 | Seconds between .trace file polls. |

### Resume

| Flag | Type | Description |
|------|------|-------------|
| `--resume [CHAINS]` | optional str | Resume from `run_state.json`. Mutually exclusive with `--overwrite`. |

### Tool

| Flag | Type | Default | Description |
|------|------|---------|-------------|
| `--pb-path` | Path | None | Directory containing `pb_mpi`, `bpcomp`, `tracecomp`, and `mpirun`. |
| `--dry-run` | flag | False | Print commands without executing. |
| `-q, --quiet` | flag | False | Suppress terminal output. |

## How It Works

### Startup

1. Validate inputs (matrix exists, mutual exclusions, parameter ranges).
2. Detect tools (`pb_mpi`, `bpcomp`, `tracecomp`, `mpirun`). Missing → exit code 3.
3. Prepare output directory. Write `run_state.json`.
4. Auto-convert FASTA to PHYLIP if needed.
5. Launch all chain processes (`subprocess.Popen`, cwd `chains/`).
6. Enter monitoring loop.

### Chain Commands

Fresh run: `mpirun -np <threads> pb_mpi -d <matrix> [model_flags] -x <sample_freq> <nsamples> <chainname>`

Resume: `mpirun -np <threads> pb_mpi <chainname>`

### Monitoring Loop

- **Every `--poll-interval` sec:** poll `.trace` files, update progress display, check `--nsamples` target.
- **Every `--monitor-freq` new samples:** trigger convergence check (bpcomp + tracecomp).
- **Ctrl+C:** soft-stop all chains, final convergence check, write `result.json`.

## Convergence Thresholds

| Metric | Good | Acceptable | Not converged |
|--------|------|------------|---------------|
| `bpcomp maxdiff` | < 0.1 | < 0.3 | ≥ 0.3 |
| `tracecomp min effsize` | > 300 | > 50 | ≤ 50 |
| `tracecomp max rel_diff` | < 0.1 | < 0.3 | ≥ 0.3 |

## Resume Semantics

`run_state.json` is the source of truth. Adding chains with `--chain-names` validates model-parameter consistency. `--resume` reads stored state, optionally overrides `--nsamples`.

## Outputs

```
runs/tree/bi/
├── chains/
│   ├── chain1.{trace,treelist,chain,param,monitor,run,log}
│   └── ...
├── convergence/
│   ├── trace_plots.pdf
│   ├── bpcomp_all.{bpdiff,bplist,con.tre}
│   └── ...
├── run_state.json
└── result.json
```

## Exit Codes

| Code | Meaning |
|------|---------|
| 0 | Success (chains completed or soft-stopped) |
| 1 | Input validation error |
| 2 | pb_mpi chain non-zero exit |
| 3 | Required tool not found |

---

# phyloai tree bi bpcomp

Run `bpcomp` once with a user-specified integer burn-in for final topology convergence analysis.

## Usage

```bash
phyloai tree bi bpcomp --chain-dir <chains_dir> [OPTIONS]
```

## Examples

```bash
# Default burn-in 1000
phyloai tree bi bpcomp --chain-dir runs/tree/bi/chains --burnin 1000

# Sub-sampling every 10 trees
phyloai tree bi bpcomp --chain-dir runs/tree/bi/chains --burnin 5000 --sample-freq 10

# Specific chains only
phyloai tree bi bpcomp --chain-dir runs/tree/bi/chains --chain-names chain1,chain2 --burnin 5000
```

## Parameters

### Input / Output

| Flag | Type | Default | Description |
|------|------|---------|-------------|
| `--chain-dir` | Path | **required** | Directory with `.chain` files. |
| `--chain-names` | str | `all` | Comma-separated names. `all` = auto-discover from `--chain-dir`. |
| `--output-dir / -o` | Path | `runs/tree/bi/bpcomp` | Output directory. |
| `--overwrite` | flag | False | Delete and recreate output directory. |

### Analysis

| Flag | Type | Default | Description |
|------|------|---------|-------------|
| `--burnin` | int ≥ 0 | 0 | Saved samples to discard. |
| `--sample-freq` | int ≥ 1 | 1 | Sub-sampling frequency after burn-in. |
| `--until` | str | `all` | Stop at sample index. `all` = full chain. |
| `--cutoff` | float (0,1) | 0.5 | Majority-rule consensus cutoff. |

### Tool

| Flag | Type | Default | Description |
|------|------|---------|-------------|
| `--pb-path` | Path | None | Directory containing `bpcomp`. |
| `--dry-run` | flag | False | Print command without executing. |
| `-q, --quiet` | flag | False | Suppress terminal output. |

## bpcomp `-x` flag construction

| Parameters | Result |
|---|---|
| `--burnin 1000` | `-x 1000` |
| `--burnin 1000 --sample-freq 10` | `-x 1000 10` |
| `--burnin 1000 --until 5000` | `-x 1000 1 5000` |

## Outputs

```
runs/tree/bi/bpcomp/
├── bpcomp.bpdiff      # Summary statistics (maxdiff, meandiff)
├── bpcomp.bplist      # Bipartition list
├── bpcomp.con.tre     # Majority-rule consensus tree
└── result.json
```

## Exit Codes

| Code | Meaning |
|------|---------|
| 0 | Success |
| 1 | Input validation error |
| 2 | bpcomp non-zero exit |
| 3 | `bpcomp` not found |

---

# phyloai tree bi tracecomp

Run `tracecomp` once with a user-specified integer burn-in for final parameter convergence analysis.

## Usage

```bash
phyloai tree bi tracecomp --chain-dir <chains_dir> [OPTIONS]
```

## Examples

```bash
# Default
phyloai tree bi tracecomp --chain-dir runs/tree/bi/chains --burnin 5000
```

## Parameters

### Input / Output

| Flag | Type | Default | Description |
|------|------|---------|-------------|
| `--chain-dir` | Path | **required** | Directory with `.trace` files. |
| `--chain-names` | str | `all` | Comma-separated names. `all` = auto-discover. |
| `--output-dir / -o` | Path | `runs/tree/bi/tracecomp` | Output directory. |
| `--overwrite` | flag | False | Delete and recreate output directory. |

### Analysis

| Flag | Type | Default | Description |
|------|------|---------|-------------|
| `--burnin` | int ≥ 0 | 0 | Saved samples to discard. |

### Tool

| Flag | Type | Default | Description |
|------|------|---------|-------------|
| `--pb-path` | Path | None | Directory containing `tracecomp`. |
| `--dry-run` | flag | False | Print command without executing. |
| `-q, --quiet` | flag | False | Suppress terminal output. |

## Outputs

```
runs/tree/bi/tracecomp/
├── tracecomp.contdiff  # Raw tracecomp output (terminal display shows annotated version)
└── result.json
```

tracecomp output is captured and printed with per-line `[good]`/`[ok]`/`[no]` annotations.

## Exit Codes

| Code | Meaning |
|------|---------|
| 0 | Success |
| 1 | Input validation error |
| 2 | tracecomp non-zero exit |
| 3 | `tracecomp` not found |

---

# phyloai tree bi readpb

Run `readpb_mpi` for posterior analysis on a single chain. Supports multiple analysis modes with automatic format conversion.

## Usage

```bash
phyloai tree bi readpb --chain <chain_path> --mode <modes> [OPTIONS]
```

## Examples

```bash
# Posterior mean exchangeabilities + site frequencies
phyloai tree bi readpb --chain chains/chain1 --mode ss,rr --burnin 5000

# PMSF simulation partition for iqtree3 --alisim
phyloai tree bi readpb --chain chains/chain1 --mode ss,rr,r --burnin 5000

# Site rates only
phyloai tree bi readpb --chain chains/chain1 --mode r --burnin 1000

# All posterior predictive checks
phyloai tree bi readpb --chain chains/chain1 --mode allppred --burnin 2000

# Dry-run
phyloai tree bi readpb --chain chains/chain1 --mode ss,rr --dry-run
```

## Parameters

### Input / Output

| Flag | Type | Default | Description |
|------|------|---------|-------------|
| `--chain` | Path | **required** | Path to chain file without extension. |
| `--mode` | str | **required** | Comma-separated analysis modes. |
| `--output-dir / -o` | Path | `runs/tree/bi/readpb` | Output directory for readpb outputs and `result.json`. |
| `--overwrite` | flag | False | Delete and recreate `--output-dir`. |

### Analysis

| Flag | Type | Default | Description |
|------|------|---------|-------------|
| `--burnin` | int ≥ 0 | 0 | Saved samples to discard. |
| `--sample-freq` | int ≥ 1 | 1 | Sub-sampling frequency after burn-in. |
| `--until` | str | `all` | Stop at sample index. |
| `--threads / -t` | int ≥ 2 | 4 | MPI processes. |

### Tool

| Flag | Type | Default | Description |
|------|------|---------|-------------|
| `--pb-path` | Path | None | Directory containing `readpb_mpi` and `mpirun`. |
| `--dry-run` | flag | False | Print commands without executing. |
| `-q, --quiet` | flag | False | Suppress terminal output. |

## `--mode` Values

| Value | readpb_mpi flag | Output | Description |
|---|---|---|---|
| `rr` | `-rr` | `.meanrr` → `.exchangeabilities` | Posterior mean exchangeabilities (PAML format). |
| `ss` | `-ss` | `.siteprofiles` → `.sitefreq` | Site-specific frequencies (IQ-TREE format). |
| `r` | `-r` | `.meansiterates` | Posterior mean site rates. |
| `sitelogl` | `-sitelogl` | `.sitelogl`, `.cpo` | Site-wise marginal log-likelihoods and cross-validation. |
| `ppred` | `-ppred` | `.ppred` | MSA simulation from posterior predictive distribution. |
| `div` | `-div` | `.div` | Diversity test (PPA-DIV). |
| `sitecomp` | `-sitecomp` | `.sitecomp` | Compositional heterogeneity (PPA-VAR). |
| `siteconvprob` | `-siteconvprob` | `.siteconvprob` | Convergence probability (PPA-CONV). |
| `comp` | `-comp` | `.comp` | Compositional homogeneity test. |
| `allppred` | `-allppred` | `.ppred` | Combined posterior predictive checks. |

`allppred` is mutually exclusive with `div`, `sitecomp`, `siteconvprob`, `comp`.

## Post-Processing

### `rr` → exchangeabilities

`<chain>.meanrr` is converted to PAML lower-triangle format (`<chain>.exchangeabilities`) with IQ-TREE-compatible AA ordering and uniform prior frequencies.

### `ss` → sitefreq

`<chain>.siteprofiles` is converted to IQ-TREE `-fs` format (`<chain>.sitefreq`), reindexing from PhyloBayes AA order to IQ-TREE order, with a `1e-8` floor and re-normalization.

### `ss,rr,r` → PMSF simulation partition

The `r` output supplies headerless, zero-based `site rate` posterior means, the chain trace supplies the posterior mean alpha using the requested burn-in/subsampling window, and the chain log supplies the discrete Gamma category count. These are combined with each site's converted frequency profile and co-generated `<chain>.exchangeabilities` model, producing one-site `+Gk{alpha}` partitions in `partition.PMSF.nex` for `iqtree3 --alisim`. The example uses `-p` for edge-proportional partitions; use `-q` instead only for edge-equal branch lengths.

```bash
iqtree3 --alisim simulated.phy -t tree.nwk -p runs/tree/bi/readpb/partition.PMSF.nex
```

## Outputs

```
runs/tree/bi/readpb/
├── chain1.meanrr              # readpb_mpi -rr
├── chain1.exchangeabilities   # PhyloAI post-processing
├── chain1.siteprofiles        # readpb_mpi -ss
├── chain1.sitefreq            # PhyloAI post-processing
├── partition.PMSF.nex          # automatic ss,rr,r PMSF simulation partitions
├── chain1.meansiterates       # readpb_mpi -r
├── chain1.sitelogl            # readpb_mpi -sitelogl
├── chain1.cpo                 # readpb_mpi -sitelogl
├── chain1.ppred               # readpb_mpi -allppred
├── ppred/chain1_ppred*.ali    # readpb_mpi -ppred
└── result.json
```

## Exit Codes

| Code | Meaning |
|------|---------|
| 0 | Success |
| 1 | Input validation error |
| 2 | readpb_mpi non-zero exit |
| 3 | `readpb_mpi` or `mpirun` not found |
