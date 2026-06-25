# PhyloAI posttree dating — Design Spec

**Date:** 2026-06-25
**Status:** Approved for implementation

---

## 1. Overview

`phyloai posttree dating` implements Bayesian molecular dating via MCMCtree
approximate likelihood calculation. It splits into two sub-subcommands that
mirror the natural two-phase workflow:

```
phyloai posttree
├── topology          (existing)
└── dating            (new subgroup)
    ├── hessian       IQ-TREE3: compute gradients/Hessian for MCMCtree
    └── mcmc          MCMCtree: Bayesian dating + full diagnostics
```

Splitting into two commands allows users to inspect and edit `mcmctree.ctl`
between phases, re-run MCMC with different parameters without re-running the
expensive IQ-TREE hessian step, and treat each step as an independently
checkpointable unit.

---

## 2. CLI Structure

### 2.1 `phyloai posttree dating hessian`

Runs IQ-TREE3 with `--dating mcmctree` to produce three files required by
MCMCtree:

| Output file                  | MCMCtree parameter   |
|------------------------------|----------------------|
| `iqtree.dummy.phy`           | `seqfile`            |
| `iqtree.rooted.nwk`          | `treefile`           |
| `iqtree.mcmctree.hessian`    | `usedata = 2` (in.BV) |

#### Parameters

| Flag | Type | Default | Description |
|------|------|---------|-------------|
| `--matrix` | Path | required | Supermatrix alignment (FASTA/PHYLIP/NEXUS) |
| `--rooted-tree` | Path | required | Rooted tree with fossil/tip calibrations in MCMCtree format |
| `--seq-type` | Choice[AA\|NT\|auto] | `auto` | Sequence type. AA→LG+F+G4, NT→HKY+G4. auto detects from alignment |
| `--model-expr` | str | None | Custom IQ-TREE model expression (e.g. `C10+F+G4`). Mutually exclusive with `--partitions` |
| `--partitions` | Path | None | Partition file (RAxML-like or NEXUS `.best_model.nex`). Mutually exclusive with `--model-expr` |
| `--prefix` | str | `iqtree` | IQ-TREE output prefix |
| `-o/--output-dir` | Path | `runs/posttree/dating/hessian` | Output directory |
| `-t/--threads` | int | 4 | IQ-TREE thread count |
| `--iqtree-path` | str | None | Explicit path to iqtree3 executable |
| `--tool-args` | str | None | Extra IQ-TREE arguments (override/extend managed flags) |
| `--overwrite` | flag | False | Delete and recreate output directory |
| `--resume` | flag | False | Resume interrupted IQ-TREE run (IQ-TREE native checkpoint) |
| `--dry-run` | flag | False | Print IQ-TREE command without executing |
| `-q/--quiet` | flag | False | Suppress terminal output except errors |

#### Model selection logic

Based on `--seq-type` (detected or explicit) and `--partitions`:

| Condition | IQ-TREE command |
|-----------|----------------|
| No partitions | `-s matrix -m <model-expr or LG+F+G4/HKY+G4> -te rooted.tre --dating mcmctree --prefix iqtree` |
| Partitions, AA, N < 10 | add `-m MF -Q partitions --mset LG -mfreq F -mrate G` |
| Partitions, NT, N < 10 | add `-m MF -Q partitions --mset HKY -mrate G` |
| Partitions, AA, N ≥ 10 | same as above + `--merge --rclusterf 10` |
| Partitions, NT, N ≥ 10 | same as above + `--merge --rclusterf 10` |

When `--partitions` is provided, `--model-expr` is ignored (they are mutually
exclusive). The partition count is read from the file before building the
command.

`--tool-args` tokens are appended last and may override any managed flag.
Blocked flags (rejected with error): `-s`, `--dating`.

#### `--rooted-tree` format

MCMCtree calibration newick format with age constraints on nodes (units: 100
Mya) and a mandatory root age constraint, e.g.:

```
(A,((B,(C,D)'>3.1<3.8'),(E,F)'>2.9<3.6'))'<4.2';
```

The root age constraint (`'<X'` or `'>X<Y'` on the outermost node) is
validated before running IQ-TREE. Missing root age → error with explanation.

---

### 2.2 `phyloai posttree dating mcmc`

Runs MCMCtree Bayesian dating using approximate likelihood (usedata=2).
Launches 4 processes in parallel: run1-posterior, run2-posterior,
run1-prior, run2-prior.

#### Parameters

| Flag | Type | Default | Description |
|------|------|---------|-------------|
| `--hessian-dir` | Path | required | Output directory from `hessian` step |
| `--clock` | Choice[1\|2\|3] | `2` | Clock model: 1=global, 2=independent rates, 3=correlated rates |
| `--burnin` | int | 100000 | MCMC burnin iterations |
| `--sample-freq` | int | 10 | Sampling frequency (record every N iterations) |
| `--nsamples` | int | 10000 | Number of samples to keep |
| `--runs` | int | 2 | Number of independent posterior runs (≥ 2 recommended) |
| `-o/--output-dir` | Path | `runs/posttree/dating/mcmc` | Output directory |
| `--mcmctree-path` | str | None | Explicit path to mcmctree executable |
| `--overwrite` | flag | False | Delete and recreate output directory |
| `--dry-run` | flag | False | Generate mcmctree.ctl and print commands without executing |
| `-q/--quiet` | flag | False | Suppress terminal output except errors |

Total MCMC iterations = `burnin + (sample_freq × nsamples)`.
Default: 100000 + (10 × 10000) = 200000 iterations, 10000 samples kept.

---

## 3. Output Directory Layout

### 3.1 `hessian` output

```
runs/posttree/dating/hessian/
├── iqtree.dummy.phy            # dummy alignment for mcmctree seqfile
├── iqtree.rooted.nwk           # rooted calibrated tree for mcmctree treefile
├── iqtree.mcmctree.hessian     # gradient/Hessian (→ rename to in.BV for mcmctree)
├── iqtree.iqtree               # IQ-TREE full report
├── iqtree.log                  # IQ-TREE log
├── iqtree.*                    # other IQ-TREE outputs (model, partition, etc.)
└── result.json
```

### 3.2 `mcmc` output

```
runs/posttree/dating/mcmc/
├── mcmctree.ctl                    # generated template; user may edit before run
├── run1/
│   ├── in.BV -> <hessian-dir>/iqtree.mcmctree.hessian
│   ├── iqtree.dummy.phy -> <hessian-dir>/iqtree.dummy.phy
│   ├── iqtree.rooted.nwk -> <hessian-dir>/iqtree.rooted.nwk
│   ├── mcmctree.ctl                # copy of top-level ctl (usedata=2, seed=-1)
│   ├── mcmctree.log                # mcmctree stdout
│   ├── mcmc.txt                    # MCMC parameter trace (progress source)
│   ├── mcmctree.out                # mcmctree summary
│   ├── SeedUsed                    # seed written by mcmctree at startup
│   ├── FigTree.tre                 # annotated dated tree
│   └── prior/
│       ├── in.BV -> <hessian-dir>/iqtree.mcmctree.hessian
│       ├── iqtree.dummy.phy -> <hessian-dir>/iqtree.dummy.phy
│       ├── iqtree.rooted.nwk -> <hessian-dir>/iqtree.rooted.nwk
│       ├── mcmctree.ctl            # usedata=0, seed fixed from ../SeedUsed
│       ├── mcmctree.log
│       ├── mcmc.txt
│       ├── mcmctree.out
│       └── FigTree.tre
├── run2/                           # identical structure to run1/
│   ├── ...
│   └── prior/
│       └── ...
└── diagnostics/
    ├── convergence/
    │   ├── posterior_times.txt     # node means + 95%CI, run1 & run2
    │   ├── prior_times.txt
    │   ├── convergence_posterior.pdf   # run1 vs run2 scatter + regression line
    │   └── convergence_prior.pdf
    ├── infinite_sites/
    │   ├── infinite_sites_run1_posterior.pdf   # mean age vs 95%CI width
    │   ├── infinite_sites_run2_posterior.pdf
    │   ├── infinite_sites_run1_prior.pdf
    │   └── infinite_sites_run2_prior.pdf
    ├── posterior_vs_prior/
    │   ├── posterior_vs_prior_run1.pdf   # posterior vs prior mean per node
    │   └── posterior_vs_prior_run2.pdf
    ├── traces/
    │   ├── mcmc_trace_run1_posterior.pdf
    │   ├── mcmc_trace_run2_posterior.pdf
    │   ├── mcmc_trace_run1_prior.pdf
    │   └── mcmc_trace_run2_prior.pdf
    └── spearman_correlations.csv       # columns: comparison, rho, pvalue
```

Symlinks in each `runN/` and `runN/prior/` point directly to `<hessian-dir>`.
`in.BV` is a symlink named `in.BV` pointing to `iqtree.mcmctree.hessian`.

---

## 4. `mcmctree.ctl` Template

Generated by Python directly (not copied from PAML examples). Template:

```
          seed = -1
       seqfile = iqtree.dummy.phy
      treefile = iqtree.rooted.nwk
       outfile = mcmctree.out

         ndata = {ndata}
       seqtype = {seqtype}   * 0: nucleotides; 1:codons; 2:AAs
       usedata = 2           * 0: no data; 1:seq like; 2:use in.BV; 3: out.BV
         clock = {clock}     * 1: global clock; 2: independent rates; 3: correlated rates
       RootAge =             * safe constraint on root age

     cleandata = 0

       BDparas = 1 1 0.1 M
   rgene_gamma = 2 20 1
  sigma2_gamma = 1 10 1

      finetune = 0: .1  .1  .1  .1 .1 .1

*** These parameters control the MCMC run
***  Note: Total number of MCMC iterations will be burnin + (sampfreq * nsample)

         print = 1
        burnin = {burnin}
      sampfreq = {sampfreq}
       nsample = {nsample}


*** The following parameters only needed to run MCMCtree with exact likelihood (usedata = 1)
*** no need to change anything for approximate likelihood (usedata = 2)

         model = 0
         alpha = 0.5
         ncatG = 4

   kappa_gamma = 6 2
   alpha_gamma = 1 1
```

- `ndata`: number of partitions (from `iqtree.dummy.phy` header, or 1 if unpartitioned)
- `seqtype`: 2 (AA) or 0 (NT) — inferred from `iqtree.dummy.phy` content.
  If sequences contain characters outside `ACGTN-` → seqtype=2 (AA),
  otherwise seqtype=0 (NT). `--seq-type` is not a parameter of `mcmc`.

- Prior run ctl: identical except `usedata = 0` and `seed = <value from SeedUsed>`

---

## 5. Run Flow

### 5.1 `hessian` flow

1. Validate `--matrix` exists and has supported extension
2. Validate `--rooted-tree` exists; check root age constraint present (regex
   on outermost node label)
3. Detect seq-type from alignment if `--seq-type auto`
4. Count partitions if `--partitions` provided; select IQ-TREE command variant
5. Build IQ-TREE command; apply `--tool-args` tokens last
6. Handle output dir lifecycle (overwrite/resume/error if non-empty)
7. Execute IQ-TREE; stream stdout to terminal (same pattern as `tree ml iqtree`)
8. Validate 3 output files exist: `iqtree.dummy.phy`, `iqtree.rooted.nwk`,
   `iqtree.mcmctree.hessian`
9. Write `result.json`

### 5.2 `mcmc` flow

```
Step 1:  Validate hessian-dir contains 3 required files
Step 2:  Infer seqtype from iqtree.dummy.phy
Step 3:  Count ndata from iqtree.dummy.phy header
Step 4:  Generate mcmctree.ctl at output-dir root
         If --dry-run: print ctl content + would-run commands, exit 0
Step 5:  Create run1/, run2/ and run1/prior/, run2/prior/ directories
         - Symlinks to hessian-dir files in each directory
         - Copy mcmctree.ctl to each runN/ (usedata=2, seed=-1)
         - Write placeholder runN/prior/mcmctree.ctl (will be completed
           once SeedUsed appears)
Step 6:  Launch run1-posterior and run2-posterior subprocesses
         stdout → runN/mcmctree.log
Step 7:  Background watcher threads poll for runN/SeedUsed
         - As soon as SeedUsed appears for runN:
           - Read seed value
           - Write runN/prior/mcmctree.ctl (usedata=0, seed=<value>)
           - Launch runN/prior subprocess
Step 8:  rich.Live progress bar monitors all 4 processes
         - One progress task per run
         - Progress = mcmc.txt line count (excluding header) / nsamples
         - Prior tasks show "waiting for seed..." until prior launched
         - Poll interval: 5 seconds
Step 9:  All 4 processes finish → stop Live display
Step 10: Generate diagnostics/ (see Section 6)
Step 11: Write result.json
```

**Ctrl+C**: hard-kill all subprocesses immediately via `proc.kill()`, re-raise
`KeyboardInterrupt`, exit non-zero. No `result.json` written.

**Progress bar design** (4 tasks in `rich.Progress`):

```
 run1-posterior  ████████████░░░░  6234/10000 samples  [02:14]
 run2-posterior  ████████████░░░░  6198/10000 samples  [02:14]
 run1-prior      waiting for seed...
 run2-prior      ████████████████  10000/10000 samples [01:58]
```

---

## 6. Diagnostics

All diagnostics generated after all 4 runs complete successfully.

### 6.1 Time tables

Parse `mcmctree.out` (lines matching `^t_`) for each run:

```
posterior_times.txt  columns: node, mean_run1, lower95_run1, upper95_run1,
                               ci_width_run1, mean_run2, lower95_run2,
                               upper95_run2, ci_width_run2
prior_times.txt      same structure
```

### 6.2 Convergence plots (`diagnostics/convergence/`)

- `convergence_posterior.pdf`: scatter of run1 vs run2 posterior node means +
  linear regression line. X=run1 mean, Y=run2 mean.
- `convergence_prior.pdf`: same for prior.
- If points fall tightly on the diagonal, runs have converged.
- Spearman ρ and p-value added to `spearman_correlations.csv`.

### 6.3 Infinite-sites plots (`diagnostics/infinite_sites/`)

- X = mean posterior (or prior) node age, Y = 95%CI width
- Points connected as line plot (ordered by X)
- Assesses data sufficiency: a straight line indicates approach to
  infinite-sites limit; scatter indicates limited molecular data
- One plot per run × posterior/prior = 4 plots total

### 6.4 Posterior vs prior plots (`diagnostics/posterior_vs_prior/`)

- X = mean posterior age, Y = mean prior age, per node
- Points connected as line plot (ordered by X)
- Diagnoses fossil calibration placement: large deviations suggest
  calibrations may be misplaced on the tree
- One plot per run = 2 plots total

### 6.5 Trace plots (`diagnostics/traces/`)

- Parse `mcmc.txt` columns (t_nX parameters = node times, mu = mean rate,
  sigma2 = rate variance, lnL = log-likelihood)
- Plot each parameter trace over samples as multi-panel figure
- One PDF per run × posterior/prior = 4 plots total

### 6.6 `spearman_correlations.csv`

Columns: `comparison`, `rho`, `pvalue`

Rows:
- `convergence_posterior` (run1 vs run2 posterior means)
- `convergence_prior`
- `infinite_sites_run1_posterior` (mean vs CI width)
- `infinite_sites_run2_posterior`
- `infinite_sites_run1_prior`
- `infinite_sites_run2_prior`
- `posterior_vs_prior_run1`
- `posterior_vs_prior_run2`

---

## 7. `doctor` Version Detection Fix

MCMCtree prints version on first line of stdout when called with no arguments:

```
MCMCTREE in paml version 4.10.10, 27 Jan 2026
```

**Fix in `env.py`:**

`TOOL_REGISTRY` entry for `mcmctree`:
```python
"mcmctree": {
    "required": False,
    "version_args": [],                              # run with no arguments
    "version_pattern": r"paml version (\d+(?:\.\d+)+)",  # custom regex
},
```

`_get_version()` gains an optional `version_pattern` override (falls back to
the existing `\d+(?:\.\d+)+` default if not set). Change is ~5 lines in
`env.py`.

---

## 8. Help Text

### `hessian` command docstring

```
Compute gradients and Hessian for MCMCtree approximate likelihood dating.

Runs IQ-TREE3 with --dating mcmctree to generate three files required by
MCMCtree:

  iqtree.dummy.phy        dummy alignment (seqfile in mcmctree.ctl)
  iqtree.rooted.nwk       rooted calibrated tree (treefile in mcmctree.ctl)
  iqtree.mcmctree.hessian gradient/Hessian matrix (rename to in.BV)

The rooted tree (--rooted-tree) must be in MCMCtree calibration format with
fossil/tip age constraints on nodes and a constrained root age, e.g.:

  (A,((B,(C,D)'>3.1<3.8'),(E,F)'>2.9<3.6'))'<4.2';

Calibration units are 100 Mya. The root age constraint is mandatory.

Model selection:
  --seq-type AA|NT|auto  detects sequence type from the alignment (default:
                         auto). Default models: LG+F+G4 (AA), HKY+G4 (NT).
  --model-expr           override with any IQ-TREE model string (e.g.
                         C10+F+G4). Mutually exclusive with --partitions.
  --partitions           partition file (RAxML-like or .best_model.nex from
                         phyloai tree ml iqtree). Recommended: <= 10
                         partitions for MCMCtree (too many partitions narrow
                         node age intervals). If > 10, PhyloAI automatically
                         merges them with --merge --rclusterf 10.

Examples:

  # Unpartitioned AA analysis (default model LG+F+G4)
  phyloai posttree dating hessian \
      --matrix concat.aa.fa --rooted-tree calib.tre

  # Custom mixture model
  phyloai posttree dating hessian \
      --matrix concat.aa.fa --rooted-tree calib.tre --model-expr C10+F+G4

  # Partitioned NT analysis (<= 10 partitions, fixed HKY+G4 per partition)
  phyloai posttree dating hessian \
      --matrix concat.nt.fa --rooted-tree calib.tre \
      --partitions loci.partitions

  # Resume interrupted IQ-TREE run
  phyloai posttree dating hessian \
      --matrix concat.aa.fa --rooted-tree calib.tre --resume
```

### `mcmc` command docstring

```
Run MCMCtree Bayesian molecular dating (approximate likelihood).

Reads the three IQ-TREE files from a completed `hessian` run and executes
MCMCtree to estimate divergence times under a Bayesian framework.

Two independent posterior runs (run1/, run2/) are launched in parallel,
each paired with a matching prior run (run1/prior/, run2/prior/) started
as soon as the posterior seed is available from SeedUsed. All four runs
use one CPU thread each (4 threads total).

A mcmctree.ctl control file is generated in the output directory before
any run starts. You may inspect and edit it freely — the file is copied
into each run directory at launch time.

MCMC settings:
  Total iterations = --burnin + (--sample-freq x --nsamples)
  Default: 100000 + (10 x 10000) = 200000 iterations, 10000 samples kept.
  Increase --nsamples (e.g. 20000) or --burnin for demanding datasets.

Clock models (--clock):
  1  Global clock (all lineages same rate)
  2  Independent rates (default; recommended for most datasets)
  3  Correlated rates (autocorrelated across branches)

Progress is tracked by polling mcmc.txt sample counts for all four runs.

Diagnostics generated after all runs complete:
  diagnostics/convergence/         run1 vs run2 scatter + regression line
  diagnostics/infinite_sites/      mean age vs 95%CI width (data sufficiency)
  diagnostics/posterior_vs_prior/  posterior vs prior age per node
  diagnostics/traces/              MCMC parameter trace plots
  diagnostics/spearman_correlations.csv

Examples:

  # Default 2-run analysis
  phyloai posttree dating mcmc \
      --hessian-dir runs/posttree/dating/hessian

  # Longer run with correlated clock
  phyloai posttree dating mcmc \
      --hessian-dir runs/posttree/dating/hessian \
      --clock 3 --burnin 200000 --nsamples 20000

  # Dry-run: inspect generated mcmctree.ctl without executing
  phyloai posttree dating mcmc \
      --hessian-dir runs/posttree/dating/hessian --dry-run
```

---

## 9. Implementation Notes

### File organisation

- CLI layer: `phyloai/cli/commands/posttree.py` — add `dating` subgroup +
  `hessian` and `mcmc` commands following existing `topology` pattern
- Library layer: `phyloai/posttree/dating_hessian.py` and
  `phyloai/posttree/dating_mcmc.py`
- Diagnostic helpers: `phyloai/posttree/dating_diagnostics.py`
- `env.py`: mcmctree version detection fix

### Dependencies

No new dependencies. Uses:
- `matplotlib` (already used elsewhere for plots)
- `scipy.stats.spearmanr` (scipy already a dependency)
- `rich` (already used)
- `subprocess`, `threading`, `symlink` from stdlib

### `result.json` schema

Follows standard PhyloAI schema. Key fields for `mcmc`:

```json
{
  "status": "success",
  "command": "phyloai posttree dating mcmc ...",
  "wall_time": 3421.5,
  "tool_versions": {"mcmctree": "4.10.10"},
  "params": { "clock": 2, "burnin": 100000, "nsamples": 10000, ... },
  "key_results": {
    "n_nodes": 12,
    "n_runs": 2,
    "convergence_rho_posterior": 0.998,
    "convergence_rho_prior": 0.997
  },
  "data": {
    "posterior_times": [...],
    "prior_times": [...],
    "diagnostics": { ... }
  }
}
```
