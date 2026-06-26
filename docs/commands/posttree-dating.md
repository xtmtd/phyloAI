# phyloai posttree dating

## Purpose

Two-step Bayesian molecular dating with MCMCtree using approximate likelihood (`usedata=2`):

1. **`hessian`** — IQ-TREE3 `--dating mcmctree` to compute gradients and Hessian.
2. **`mcmc`** — MCMCtree with the IQ-TREE output, running independent posterior + prior chains in parallel with real-time progress and diagnostic plots.

Only approximate likelihood (`usedata=2`) is implemented. Exact likelihood (`usedata=1`, sequence data) and `usedata=3` (output `in.BV`) are not supported by PhyloAI.

## Quick Start

```bash
# Step 1: Compute Hessian (AA, auto model LG+F+G4)
phyloai posttree dating hessian --matrix concat.aa.fa --rooted-tree calib.tre

# Step 2: Run MCMCtree (2 independent runs, independent rates clock)
phyloai posttree dating mcmc --hessian-dir runs/posttree/dating/hessian

# Step 3: Review diagnostics in runs/posttree/dating/mcmc/diagnostics/
```

## Calibrated Tree Format

The `--rooted-tree` must be a NEWICK tree with MCMCtree-style calibration annotations and a mandatory root age constraint. Calibration units are 100 Mya (e.g. `'>3.1<3.8'` means 310–380 Mya).

```
# Single lower-bound constraint (310 Mya for the C-D clade)
(A,((B,(C,D)'>3.1<3.8'),(E,F)'>2.9<3.6'))'<4.2';

# Point calibration with soft bounds (uniform distribution)
(A,((B,(C,D)'>0.5<0.7'),(E,F)'>0.3<0.5'))'<1.0';

# Tip dating: node label marks extant taxon
(A_tip,((B,(C,D_tip)'>3.1<3.8'),E_tip)'<4.2')'<5.0';
```

Constraints use `>`, `<`, or both: `'>L<U'` defines lower `L` and upper `U` bounds. The root (outermost node) must have a constraint — PhyloAI rejects trees without one.

## Hessian Step

Computes the gradient and Hessian matrix IQ-TREE needs for `usedata=2` (approximate likelihood). Produces three files consumed by the `mcmc` step.

### Usage

```bash
phyloai posttree dating hessian --matrix concat.aa.fa --rooted-tree calib.tre [OPTIONS]
```

### Inputs

| Flag | Default | Description |
|------|---------|-------------|
| `--matrix PATH` | *(required)* | Alignment: FASTA (`.fa`/`.fas`/`.fasta`/`.faa`/`.fna`/`.aln`), PHYLIP, or NEXUS. `--seq-type auto` reads all three formats via PhyloAI's shared format detector. |
| `--rooted-tree PATH` | *(required)* | MCMCtree calibration tree (see above). Maps to IQ-TREE `-te`. |
| `--seq-type AA\|NT\|auto` | `auto` | Sequence type. `AA` → `LG+F+G4`; `NT` → `GTR+G4`. `auto` detects from alignment content. |
| `--model-expr STR` | — | IQ-TREE model string (e.g. `C10+F+G4`). Mutually exclusive with `--partitions`. |
| `--partitions PATH` | — | Partition file (RAxML-like, NEXUS `.best_model.nex`, or cluster file). `< 10` partitions run directly; `>= 10` are auto-merged with `--merge --rclusterf 10`. Maps to `-Q`. |
| `-o, --output-dir PATH` | `runs/posttree/dating/hessian` | Output directory. |
| `-t, --threads INT` | `4` | Threads for IQ-TREE (`-T`). |
| `--iqtree-path PATH` | — | Override auto-detected `iqtree3` binary. |
| `--tool-args STR` | — | Extra IQ-TREE arguments. Blocked: `-s`, `--dating`, `-te`, `--prefix`. |
| `--overwrite` | off | Delete and recreate output directory. |
| `--resume` | off | Resume interrupted IQ-TREE run (IQ-TREE native checkpoint). |
| `--dry-run` | off | Print IQ-TREE command without executing. |
| `-q, --quiet` | off | Suppress terminal output except errors. |

### Outputs

```
runs/posttree/dating/hessian/
├── result.json
├── iqtree.dummy.phy          → mcmc seqfile
├── iqtree.rooted.nwk         → mcmc treefile
└── iqtree.mcmctree.hessian   → renamed to in.BV in each run dir
```

The `iqtree.*` output files use a fixed prefix — PhyloAI blocks `--prefix` in `--tool-args` because the `mcmc` step depends on these exact filenames.

### Model Selection Logic

| Condition | Model |
|-----------|-------|
| `--model-expr` set | Uses that expression directly (`-m <expr>`) |
| `--partitions` set, AA | `-m MF -Q <file> --mset LG -mfreq F -mrate G` |
| `--partitions` set, NT | `-m MF -Q <file> --mset GTR -mrate G` |
| Unpartitioned, AA | `-m LG+F+G4` |
| Unpartitioned, NT | `-m GTR+G4` |

When `--partitions` is used with ModelFinder (`-m MF`), the model search space is constrained to the simplest model families appropriate for dating — no free-rate models (`+R`) or complex mixtures that would inflate computation time.

## MCMC Step

Runs MCMCtree Bayesian dating using the hessian output. Launches independent posterior runs in parallel, each paired with a prior-predictive run (`usedata=0`) using the same random seed.

### Usage

```bash
phyloai posttree dating mcmc --hessian-dir runs/posttree/dating/hessian [OPTIONS]
```

### Inputs

| Flag | Default | Description |
|------|---------|-------------|
| `--hessian-dir PATH` | *(required)* | Directory containing `iqtree.dummy.phy`, `iqtree.rooted.nwk`, `iqtree.mcmctree.hessian`. |
| `-o, --output-dir PATH` | `runs/posttree/dating/mcmc` | Output directory. |
| `--runs INT` | `2` | Number of independent posterior runs (each paired with a prior run). |
| `--clock 1\|2\|3` | `2` | Clock model: `1` = global, `2` = independent rates, `3` = correlated rates. Ignored when `--ctl` is provided. |
| `--burnin INT` | `100000` | MCMC burn-in iterations. Ignored when `--ctl` is provided. |
| `--sample-freq INT` | `10` | Sampling frequency (iterations per sample). Ignored when `--ctl` is provided. |
| `--nsamples INT` | `10000` | Samples kept post-burnin. Ignored when `--ctl` is provided. |
| `--ctl PATH` | — | Pre-configured `mcmctree.ctl`. Mutually exclusive with non-default `--clock`/`--burnin`/`--sample-freq`/`--nsamples`. |
| `--mcmctree-path PATH` | — | Override auto-detected `mcmctree` binary. |
| `--overwrite` | off | Delete and recreate output directory. |
| `--dry-run` | off | Print generated ctl and exit without running. |
| `-q, --quiet` | off | Suppress MCMC log tail and progress bar output. |

**Total iterations** = `--burnin` + (`--sample-freq` × `--nsamples`). Default: 200,000.

### Clock Model Guide

| Clock | Model | When to Use |
|-------|-------|-------------|
| `1` | Global | All lineages evolve at the same rate (rarely true). Fastest, simplest. |
| `2` | Independent (recommended) | Each branch has its own rate, drawn from a log-normal distribution. Good default for most datasets. |
| `3` | Correlated | Rates on neighbouring branches are correlated. More biologically realistic but computationally heavier. |

### Outputs

```
runs/posttree/dating/mcmc/
├── result.json
├── mcmctree.ctl                   # generated template (or copy of --ctl)
├── run1/
│   ├── mcmctree.ctl               # posterior ctl with injected seed
│   ├── mcmc.txt                   # MCMC parameter trace
│   ├── mcmctree.out               # divergence time summary
│   ├── mcmctree.log
│   ├── FigTree.tre                # annotated dated tree
│   ├── FigTree.node.tre           # node-labeled tree (parsed from mcmctree.out)
│   ├── SeedUsed                   # seed written by MCMCtree
│   ├── iqtree.dummy.phy -> <hessian-dir>/
│   ├── iqtree.rooted.nwk -> <hessian-dir>/
│   ├── in.BV -> <hessian-dir>/iqtree.mcmctree.hessian
│   └── prior/
│       ├── mcmctree.ctl           # usedata=0, same seed as posterior
│       ├── mcmc.txt
│       ├── mcmctree.out
│       └── FigTree.node.tre
├── run2/                          # identical structure
│   └── ...
└── diagnostics/
    ├── traces/                    # per-parameter MCMC trace PDFs
    ├── convergence/
    │   ├── posterior_times.csv    # all runs combined
    │   ├── prior_times.csv
    │   └── convergence_*_runX_vs_runY.pdf  # scatter plots + fit line
    ├── infinite_sites/            # mean age vs 95% CI width plots
    ├── posterior_vs_prior/        # posterior vs prior mean age per node
    └── spearman_correlations.csv
```

### Diagnostic Interpretations

| Diagnostic | What to Look For |
|------------|-----------------|
| **Traces** (`traces/`) | Parameter traces should be well-mixed ("hairy caterpillar"), not trending or stuck. Converge early, stay within a stable range. |
| **Convergence** (`convergence/`) | Mean divergence times from independent runs should align along y=x. Slope near 1.0, Spearman's ρ near 1.0, low RMSE. Systematic deviations suggest insufficient burnin or chain length. |
| **Infinite-sites** (`infinite_sites/`) | CI widths should not increase dramatically with node age. A strong positive slope means older nodes are poorly constrained — check fossil calibrations. |
| **Posterior vs Prior** | Posterior estimates that strongly deviate from the prior mean the data is informative. If posterior ≈ prior for most nodes, the data carries little temporal signal. |
| **Spearman correlations** (`spearman_correlations.csv`) | Reports ρ, p-value, linear fit (slope, intercept, RMSE) for every pairwise comparison. Use as a quick convergence check. |

## Examples

```bash
# ── Hessian ───────────────────────────────────────

# 1. Unpartitioned AA with default model (LG+F+G4)
phyloai posttree dating hessian --matrix concat.aa.fa --rooted-tree calib.tre

# 2. Custom mixture model
phyloai posttree dating hessian --matrix concat.aa.fa --rooted-tree calib.tre --model-expr C10+F+G4

# 3. Explicit NT sequence type
phyloai posttree dating hessian --matrix concat.nt.fa --rooted-tree calib.tre --seq-type NT

# 4. Partitioned analysis (< 10 partitions)
phyloai posttree dating hessian --matrix concat.aa.fa --rooted-tree calib.tre \
    --partitions partitions.nex -o runs/dating/hessian

# 5. Custom output directory
phyloai posttree dating hessian --matrix concat.aa.fa --rooted-tree calib.tre \
    -o runs/dating/hessian

# 6. Resume interrupted run (IQ-TREE native checkpoint)
phyloai posttree dating hessian --matrix concat.aa.fa --rooted-tree calib.tre --resume

# ── MCMC ──────────────────────────────────────────

# 7. Default: 2 runs, independent rates clock
phyloai posttree dating mcmc --hessian-dir runs/dating/hessian

# 8. Three runs, correlated clock, longer chain
phyloai posttree dating mcmc --hessian-dir runs/dating/hessian \
    --runs 3 --clock 3 --burnin 200000 --nsamples 20000

# 9. Dry-run: inspect generated ctl before running
phyloai posttree dating mcmc --hessian-dir runs/dating/hessian --dry-run

# 10. Use a pre-configured mcmctree.ctl for full control
phyloai posttree dating mcmc --hessian-dir runs/dating/hessian --ctl my_run.ctl

# 11. Single-run (no convergence diagnostics)
phyloai posttree dating mcmc --hessian-dir runs/dating/hessian --runs 1

# 12. Overwrite previous results
phyloai posttree dating mcmc --hessian-dir runs/dating/hessian --overwrite
```

## Exit Codes

| Code | Meaning |
|------|---------|
| `0` | Success |
| `1` | Input error (missing/empty files, conflicting flags, non-empty output dir) |
| `2` | Tool failure (IQ-TREE non-zero exit, MCMCtree posterior failure) |
| `3` | Environment error (`iqtree3` or `mcmctree` not found) |

## Warnings / Errors

| Condition | Behaviour |
|-----------|-----------|
| Empty or missing `hessian-dir` | Exit 1 with message listing which files are absent. |
| `--matrix` or `--rooted-tree` does not exist or is empty | Exit 1 before IQ-TREE starts. |
| `--rooted-tree` missing root age constraint | Exit 1; outermost node must have `'>L<U'` or `'<U'`. |
| `--model-expr` and `--partitions` together | Exit 1; mutually exclusive. |
| `--ctl` with non-default `--clock`/`--burnin`/`--sample-freq`/`--nsamples` | Exit 1; these flags only apply when generating ctl from scratch. |
| Non-empty output directory without `--overwrite` | Exit 1; use `--overwrite` to replace. |
| IQ-TREE returns non-zero | Status `error` in `result.json`; stderr captured in `data.warnings`. |
| IQ-TREE produced empty output files | Recorded as a warning; possible crash mid-write. |
| IQ-TREE report missing "Total CPU time used" | Warning; IQ-TREE may have been interrupted. |
| `mcmctree` binary not found | Exit 3 (env error). Run `phyloai doctor`. |
| Posterior exit code ≠ 0 | Top-level `status` becomes `"error"`. |
| Prior exit code ≠ 0 or missing output files | Recorded in `data.warnings`; posterior results and diagnostics still produced. |
| `mcmc.txt` or `mcmctree.out` empty | Recorded in `data.warnings`; diagnostics degraded gracefully. |
| `--runs=1` | Convergence plots skipped; other diagnostics still generated. |

## Notes

- **usedata=2 only.** `usedata=1` (exact likelihood via sequence data) and `usedata=3` (output `in.BV`) are not implemented in PhyloAI.
- **Seed injection.** PhyloAI generates a unique random seed per run (`random.randint(1, 2³¹-1)`), injected into the ctl before launch. Posterior and prior share the same seed within a run, enabling proper comparison.
- **ndata.** The number of data blocks is always counted directly from `iqtree.dummy.phy` — never from `result.json`. This ensures correctness when IQ-TREE merges partitions (`--merge --rclusterf 10`).
- **seqtype.** Read from `hessian/result.json` (`params.seq_type`), falling back to auto-detection from dummy.phy content if the result.json is missing or malformed.
- **Version detection.** The `mcmctree` binary is found via `shutil.which`; the version is extracted by running the binary with no arguments and matching `paml version (\d+(?:\.\d+)+)`. Falls back to parsing `mcmctree.log` after runs complete.
- **Prefix locking.** PhyloAI hardcodes `--prefix iqtree` for the hessian step. This ensures the mcmc step can find `iqtree.dummy.phy`, `iqtree.rooted.nwk`, and `iqtree.mcmctree.hessian` without requiring user to track naming. `--prefix` is blocked in `--tool-args`.
- **Custom ctl (`--ctl`).** PhyloAI always symlinks the standard hessian files (`iqtree.dummy.phy`, `iqtree.rooted.nwk`, `in.BV`) into each `runN/` directory. Custom relative `seqfile`/`treefile` paths in the ctl are resolved relative to the ctl file's directory and symlinked into each `runN/` as well.
- **Diagnostics.** Convergence plots require `--runs >= 2`. All scatter plots use a single dashed fit line with equation; node labels show `nXX` format. Posterior-vs-prior plots use equal-length axes for direct comparison.
- **`result.json`** follows the PhyloAI JSON Output Standard.
- **Prior-predictive run.** The prior `mcmctree.ctl` is derived from the posterior ctl by replacing `usedata = 2` with `usedata = 0` and preserving the same seed. All other parameters remain identical, ensuring the prior and posterior are directly comparable.
- **OMP_NUM_THREADS.** Always set to `1` for MCMCtree subprocesses. MCMC samplers do not benefit from multi-threading and can produce incorrect results with `OMP_NUM_THREADS > 1`.
