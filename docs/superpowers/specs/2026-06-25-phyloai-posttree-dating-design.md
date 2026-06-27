# PhyloAI posttree dating — Design Spec

**Date:** 2026-06-25
**Last updated:** 2026-06-26 (implemented; --prefix removed, ndata counting from dummy.phy, --seq-type auto supports FASTA/PHYLIP/NEXUS, help/examples cleaned up)
**Status:** Implemented

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
| `--seq-type` | Choice[AA\|NT\|auto] | `auto` | Sequence type. AA→LG+F+G4, NT→GTR+G4. auto detects from FASTA, PHYLIP, and NEXUS via shared format helpers. |
| `--model-expr` | str | None | Custom IQ-TREE model expression (e.g. `C10+F+G4`). Mutually exclusive with `--partitions` |
| `--partitions` | Path | None | Partition file (RAxML-like or NEXUS `.best_model.nex`, or clusters from `phyloai pretree filter cluster`). Mutually exclusive with `--model-expr` |
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
| No partitions | `-s matrix -m <model-expr or LG+F+G4/GTR+G4> -te rooted.tre --dating mcmctree --prefix iqtree` |
| Partitions, AA, N < 10 | add `-m MF -Q partitions --mset LG -mfreq F -mrate G` |
| Partitions, NT, N < 10 | add `-m MF -Q partitions --mset GTR -mrate G` |
| Partitions, AA, N ≥ 10 | same as above + `--merge --rclusterf 10` |
| Partitions, NT, N ≥ 10 | same as above + `--merge --rclusterf 10` |

When `--partitions` is provided, `--model-expr` is ignored (they are mutually
exclusive). The partition count is read from the file before building the
command.

`--tool-args` tokens are appended last and may override any managed flag.
**Blocked flags** (rejected with exit code 1 and a clear message):
`-s`, `--dating`, `-te`, `--prefix`. These four flags define the
PhyloAI → mcmctree contract: the hessian step hardcodes IQ-TREE prefix
`iqtree` and must emit
`iqtree.dummy.phy`, `iqtree.rooted.nwk`, and `iqtree.mcmctree.hessian`
under the output directory, and use the
calibrated tree the user supplied via `--rooted-tree` (controlled by
`-te`). Letting `--tool-args` silently override any of these would
break the mcmc step's ability to find the required files or use the
correct tree. All other IQ-TREE flags remain strategy parameters and
pass through `--tool-args` unchanged.

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
Launches `2 × --runs` processes: a posterior + prior pair for each run
(default 4 processes: run1-posterior, run2-posterior, run1-prior,
run2-prior).

#### Parameters

| Flag | Type | Default | Description |
|------|------|---------|-------------|
| `--hessian-dir` | Path | required | Output directory from `hessian` step |
| `--ctl` | Path | None | Use this mcmctree.ctl as-is instead of generating one. PhyloAI still copies it into each runN/ and runs MCMCtree. The provided ctl's `seqfile`/`treefile` are honored; symlinks to hessian files under runN/ are placed alongside so MCMCtree can find `in.BV` and `iqtree.*` regardless of `seqfile`/`treefile` content. Passing any of `--clock`/`--burnin`/`--sample-freq`/`--nsamples` with a **non-default value** together with `--ctl` exits code 1 (those flags only affect the generated template; when `--ctl` is provided no template is generated). Passing the default value (e.g. `--clock 2`) is allowed but silently ignored. |
| `--clock` | Choice[1\|2\|3] | `2` | Clock model: 1=global, 2=independent rates, 3=correlated rates. Ignored when `--ctl` is provided. |
| `--burnin` | int | 100000 | MCMC burnin iterations. Ignored when `--ctl` is provided. |
| `--sample-freq` | int | 10 | Sampling frequency (record every N iterations). Ignored when `--ctl` is provided. |
| `--nsamples` | int | 10000 | Number of samples to keep. Ignored when `--ctl` is provided. |
| `--runs` | int | 2 | Number of independent posterior runs. Each run is paired with a matching prior run. **I7**: `--runs=1` is allowed for exploratory runs that skip convergence diagnostics; `--runs≥3` is supported by computing pairwise convergence (run1 vs run2, run1 vs run3, run2 vs run3, …) with the per-pair Spearman ρ in `spearman_correlations.csv`. |
| `-o/--output-dir` | Path | `runs/posttree/dating/mcmc` | Output directory |
| `--mcmctree-path` | str | None | Explicit path to mcmctree executable |
| `--overwrite` | flag | False | Delete and recreate output directory |
| `--dry-run` | flag | False | Generate mcmctree.ctl and print commands without executing (no-op when `--ctl` is provided: just print which ctl would be used) |
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
├── mcmctree.ctl                    # generated template copy; edit via --dry-run + --ctl workflow
├── run1/
│   ├── in.BV -> <hessian-dir>/iqtree.mcmctree.hessian
│   ├── iqtree.dummy.phy -> <hessian-dir>/iqtree.dummy.phy
│   ├── iqtree.rooted.nwk -> <hessian-dir>/iqtree.rooted.nwk
│   ├── mcmctree.ctl                # copy of top-level ctl, seed injected with random int
│   ├── mcmctree.log                # mcmctree stdout
│   ├── mcmc.txt                    # MCMC parameter trace (progress source)
│   ├── mcmctree.out                # mcmctree summary
│   ├── FigTree.tre                 # annotated dated tree (mcmctree native)
│   ├── FigTree.node.tre            # same tree, internal labels as bare integers
│   └── prior/
│       ├── in.BV -> <hessian-dir>/iqtree.mcmctree.hessian
│       ├── iqtree.dummy.phy -> <hessian-dir>/iqtree.dummy.phy
│       ├── iqtree.rooted.nwk -> <hessian-dir>/iqtree.rooted.nwk
│       ├── mcmctree.ctl            # usedata=0, same random seed as posterior
│       ├── mcmctree.log
│       ├── mcmctree.log
│       ├── mcmc.txt
│       ├── mcmctree.out
│       ├── FigTree.tre
│       └── FigTree.node.tre
├── run2/                           # identical structure to run1/
│   ├── ...
│   └── prior/
│       └── ...
└── diagnostics/
    ├── convergence/
    │   ├── posterior_times.csv                    # run1 vs run2 (canonical pair)
    │   ├── prior_times.csv
    │   ├── convergence_posterior_run1_vs_run2.pdf # scatter + regression + y=x
    │   └── convergence_prior_run1_vs_run2.pdf
    │   #   (extra pairs generated as *_runA_vs_runB.pdf when --runs >= 3)
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
    └── spearman_correlations.csv       # columns: comparison, rho, pvalue, slope, intercept, rmse
```

Symlinks in each `runN/` and `runN/prior/` point directly to `<hessian-dir>`.
`in.BV` is a symlink named `in.BV` pointing to `iqtree.mcmctree.hessian`.

---

## 4. `mcmctree.ctl` Template

Generated by Python directly (not copied from PAML examples). The template
includes every parameter that the upstream `mcmctree.sh` reference script
exposes, with the inline `* ...` comments preserved verbatim. These comments
are harmless if MCMCtree ignores the parameter (e.g. `model`/`alpha`/`ncatG`
when `usedata = 2`) and they help users edit the ctl confidently.

```
          seed = -1
       seqfile = iqtree.dummy.phy
      treefile = iqtree.rooted.nwk
       outfile = mcmctree.out

         ndata = {ndata}
       seqtype = {seqtype}  * 0: nucleotides; 1:codons; 2:AAs
       usedata = 2    * 0: no data; 1:seq like; 2:use in.BV; 3: out.BV
         clock = {clock}    * 1: global clock; 2: independent rates; 3: correlated rates
       RootAge =   * safe constraint on root age, used if no fossil for root.

       BDparas = 1 1 0.1 M   * birth, death, sampling
   rgene_gamma = 2 20 1   * gamma prior for overall rates for genes
  sigma2_gamma = 1 10 1    * gamma prior for sigma^2     (for clock=2 or 3)

      finetune = 0: .1  .1  .1  .1 .1 .1  * auto (0 or 1) : times, musigma2, rates, mixing, paras, FossilErr

*** These parameters control the MCMC run
***  Note: Total number of MCMC iterations will be burnin + (sampfreq * nsample)

         print = 1
        burnin = {burnin}
      sampfreq = {sampfreq}
       nsample = {nsample}


*** The following parameters only needed to run MCMCtree with exact likelihood (usedata = 1)
*** no need to change anything for approximate likelihood (usedata = 2)

         model = 0    * 0:JC69, 1:K80, 2:F81, 3:F84, 4:HKY85
         alpha = 0.5    * alpha for gamma rates at sites
         ncatG = 4    * No. categories in discrete gamma

     cleandata = 0  * remove sites with ambiguity data (1:yes, 0:no)?

   kappa_gamma = 6 2      * gamma prior for kappa
   alpha_gamma = 1 1      * gamma prior for alpha

*** Note: Make your window wider (100 columns) before running the program.
```

- `ndata`: number of data blocks in `iqtree.dummy.phy`. **Always counted
  directly from the file** — this is the ground truth, especially when
  `--merge --rclusterf 10` reduces ≥ 10 original partitions to fewer
  megapartitions. `mcmc` counts blocks matching `^\s+\d+\s+\d+\s*$`,
  minimum 1. The `hessian` step also stores the original partition count
  in `result.json` (`params.n_partitions`; 1 for unpartitioned), but
  `mcmc` does **not** read `ndata` from there — it reads `seq_type` from
  `params.seq_type` but always counts `ndata` from the actual dummy.phy.
- `seqtype`: 2 (AA) or 0 (NT). Sourced from `hessian`'s `result.json`
  (`params.seq_type`); fallback scans `iqtree.dummy.phy` for non-ACGT
  characters. `--seq-type` is not a parameter of `mcmc`.
- `cleandata = 0`: default for IQ-TREE-emitted dummy.phy, where MCMCtree's
  internal likelihood is not used. Leave at 0 unless running with
  `usedata = 1`/`usedata = 3` (exact likelihood paths, not in this spec's
  scope).
- Prior run ctl: identical except `usedata = 0` and the random seed
  injected by PhyloAI (shared with the posterior run).

---

## 5. Run Flow

### 5.1 `hessian` flow

1. Validate `--matrix` exists and is a supported format (FASTA, PHYLIP, NEXUS).
2. Validate `--rooted-tree` exists; check root age constraint present (regex
   on outermost node label)
3. Detect seq-type from FASTA alignment if `--seq-type auto`
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
Step 2:  Read seq_type and ndata from <hessian-dir>/result.json
         (hessian records them under params.seq_type / params.n_partitions).
         Fallback (older hessian output without these fields, or
         out-of-PhyloAI hessian dir): infer seqtype from iqtree.dummy.phy
         content and count ndata from header blocks.
Step 3:  If --ctl provided:
            - Validate --ctl exists
            - Reject if --clock/--burnin/--sample-freq/--nsamples are
              passed with non-default values (those flags only affect the
              generated template; Click cannot distinguish an implicit
              default from an explicit one, so only non-default values
              trigger the conflict check).
            - Resolve ctl_text = contents of --ctl file
            - Skip Steps 4 (no generation)
         Else:
            - Generate mcmctree.ctl at output-dir root from template,
              substituting {seqtype}/{ndata}/{clock}/{burnin}/{sampfreq}/{nsample}
            - ctl_text = generated template
Step 4:  If --dry-run:
            - Print ctl_text + would-run commands, exit 0
Step 5:  Generate a random seed per run (positive 32-bit int).
         For each run directory:
           - Inject the seed into ctl_text via regex
           - Write posterior runN/mcmctree.ctl
           - Derive prior runN/prior/mcmctree.ctl (usedata=0, same seed)
         All symlinks to hessian-dir files are created in both runN/
         and runN/prior/.
Step 6:  Launch all posterior and prior subprocesses in parallel
         stdout → runN/mcmctree.log
Step 7:  rich.Live progress bar monitors every launched process
         (2 × --runs total: one posterior + one prior task per run)
         - Progress = mcmc.txt line count (excluding header) / nsamples
         - Poll interval: 5 seconds
Step 8:  All processes finish → stop Live display
Step 9:  Collect return codes; if any posterior non-zero → status=error,
         error message names the failed run/dir, data.warnings lists
         every failure. Validation: mcmc.txt and mcmctree.out non-empty
         for every run; failures recorded as warnings even when
         returncode=0 (silent truncation).
Step 10: Generate diagnostics/ (see Section 6)
Step 11: Write result.json
```

**Ctrl+C**: hard-kill all subprocesses immediately via `proc.kill()`, re-raise
`KeyboardInterrupt`, exit non-zero. No `result.json` written.

**Progress bar design** (4 tasks in `rich.Progress`):

```
 run1-posterior  ████████████░░░░  6234/10000 samples  [02:14]
 run2-posterior  ████████████░░░░  6198/10000 samples  [02:14]
 run1-prior      ████████████████  10000/10000 samples [01:58]
 run2-prior      █████████░░░░░░░░  4231/10000 samples  [01:21]
```

---

## 6. Diagnostics

Diagnostics are generated after all runs complete. Diagnostics use
whatever output is available: a run whose prior failed is still included
(posterior node times and traces); a run with zero nodes is highlighted
as a warning in `spearman_correlations.csv`. Convergence plots and
pairwise metrics are generated for every pair of runs that both produced
valid `mcmctree.out` files (`--runs=1` skips convergence entirely). The
diagnostics section degrades gracefully with partial output.

### 6.1 Time tables

Parse `mcmctree.out` (lines matching `^t_`) for each run:

```
posterior_times.csv  columns: node, mean_run1, lower95_run1, upper95_run1,
                               ci_width_run1, mean_run2, lower95_run2,
                               upper95_run2, ci_width_run2
prior_times.csv      same structure
```

### 6.2 Convergence plots (`diagnostics/convergence/`)

- `convergence_posterior_<run_a>_vs_<run_b>.pdf`: scatter of two
  posterior runs' node means + linear regression line + y=x reference
  line. X = run_a mean, Y = run_b mean.
- `convergence_prior_<run_a>_vs_<run_b>.pdf`: same for priors.
- For `--runs=2` (default) only one pair is generated:
  `run1_vs_run2`.
- For `--runs≥3` all `C(n_runs, 2)` pairs are generated.
- For `--runs=1` (I7) no convergence plots are written; the
  `convergence` entry in `data.diagnostics` is `{"status": "skipped",
  "reason": "n_runs=1"}`.
- D3: each plot annotates the linear-regression slope/intercept and
  RMSE in the legend alongside Spearman ρ. A point cloud lying on the
  y=x diagonal (slope ≈ 1, intercept ≈ 0, ρ ≈ 1, low RMSE) indicates
  convergence; systematic offset between runs shows up as a slope ≠ 1
  or non-zero intercept that pure ρ would hide.
- Per-pair Spearman ρ and p-value added to `spearman_correlations.csv`.

### 6.3 Infinite-sites plots (`diagnostics/infinite_sites/`)

- X = mean node age, Y = 95%CI width
- Points connected as line plot (ordered by X)
- One plot per available run × posterior/prior pair. For the default
  `--runs=2` this is 4 plots; `--runs=1` produces 2; `--runs=3` produces
  6.
- Interpretation differs by distribution type:
  - **Posterior** (`usedata=2`): a straight line (CI width proportional to
    age) indicates the infinite-sites limit is approached, i.e. additional
    molecular data would not substantially improve precision. Scatter or
    non-linear pattern suggests limited molecular data.
  - **Prior** (`usedata=0`): reflects fossil calibration information alone.
    A straight line indicates fossil constraints are internally consistent
    and informative. Scatter or non-linearity suggests fossil calibrations
    may be insufficient or conflicting.
- Comparing posterior vs prior infinite-sites plots for the same run
  reveals how much the molecular data updates the fossil-only prior.

### 6.4 Posterior vs prior plots (`diagnostics/posterior_vs_prior/`)

- X = mean posterior age, Y = mean prior age, per node
- Points connected as line plot (ordered by X)
- Diagnoses fossil calibration placement: large deviations suggest
  calibrations may be misplaced on the tree
- One plot per available run. Default `--runs=2` → 2 plots total.

### 6.5 Trace plots (`diagnostics/traces/`)

- Parse `mcmc.txt` columns (t_nX parameters = node times, mu = mean rate,
  sigma2 = rate variance, lnL = log-likelihood)
- Plot each parameter trace over samples as multi-panel figure
- One PDF per available run × posterior/prior. Default `--runs=2` → 4 plots.

### 6.6 `spearman_correlations.csv`

Columns: `comparison`, `rho`, `pvalue`, `slope`, `intercept`, `rmse`

Rows:
- `convergence_posterior_run1_vs_run2`, `convergence_prior_run1_vs_run2`
  (and every extra pair when `--runs≥3`)
- `infinite_sites_run1_posterior` (mean vs CI width)
- `infinite_sites_run2_posterior`
- `infinite_sites_run1_prior`
- `infinite_sites_run2_prior`
- `posterior_vs_prior_run1`
- `posterior_vs_prior_run2`

D3: `slope`/`intercept` come from the per-plot linear fit (only for
rows that have an X/Y pair; infinite-sites and posterior-vs-prior rows
fill `slope`/`intercept`/`rmse` with the same fit). For runs with
fewer than 3 valid points (D6) all five numeric columns are written as
empty strings and a warning is added to `data.warnings`.

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
  --seq-type AA|NT|auto  Sequence type (default: auto — reads FASTA,
                         PHYLIP, and NEXUS through shared format helpers).
  --model-expr           Override with any IQ-TREE model string (e.g.
                         C10+F+G4). Mutually exclusive with --partitions.
  --partitions           Partition file (RAxML-like, NEXUS .best_model.nex,
                         or cluster file). < 10 partitions run directly;
                         >= 10 are auto-merged.

Examples:

  # Unpartitioned AA analysis (default model LG+F+G4)

  phyloai posttree dating hessian --matrix concat.aa.fa --rooted-tree calib.tre

  # Custom mixture model

  phyloai posttree dating hessian --matrix concat.aa.fa --rooted-tree calib.tre --model-expr C10+F+G4

  # Partitioned NT analysis (< 10 partitions, fixed GTR+G4 per partition)

  phyloai posttree dating hessian --matrix concat.nt.fa --rooted-tree calib.tre --partitions loci.partitions

  # Resume interrupted IQ-TREE run

  phyloai posttree dating hessian --matrix concat.aa.fa --rooted-tree calib.tre --resume
```

### `mcmc` command docstring

```
Run MCMCtree Bayesian molecular dating (approximate likelihood only).

Reads the three IQ-TREE files from a completed `hessian` run and executes
MCMCtree to estimate divergence times. Uses usedata=2 (gradient/Hessian
from IQ-TREE in.BV). Does NOT implement usedata=1 or usedata=3.

Two independent posterior runs (run1/, run2/) launch in parallel, each
with a matching prior run using the same seed.

Review the generated ctl with --dry-run, then edit the copy and re-run
with --ctl edited.ctl to customize parameters beyond the built-in flags.

MCMC settings: Total iterations = --burnin + (--sample-freq x --nsamples).
Default: 200000 iterations, 10000 samples kept.

Clock models:
  1  Global clock
  2  Independent rates (default; recommended)
  3  Correlated rates

Diagnostics generated after all runs complete (under diagnostics/):
  convergence/       posterior_times.csv + run1-vs-run2 scatter plots
  infinite_sites/    mean age vs 95%CI width (data-sufficiency check)
  posterior_vs_prior/  posterior-vs-prior mean age per node
  traces/            MCMC parameter trace plots
  spearman_correlations.csv

Examples:

  # Default 2-run analysis

  phyloai posttree dating mcmc --hessian-dir runs/posttree/dating/hessian

  # Longer run with correlated clock

  phyloai posttree dating mcmc --hessian-dir runs/posttree/dating/hessian --clock 3 --burnin 200000 --nsamples 20000

  # Use a custom mcmctree.ctl

  phyloai posttree dating mcmc --hessian-dir runs/posttree/dating/hessian --ctl my_run.ctl

  # Dry-run: inspect generated ctl without executing

  phyloai posttree dating mcmc --hessian-dir runs/posttree/dating/hessian --dry-run
```

---

## 9. Implementation Notes

### Tool detection

`run_mcmc` finds the mcmctree binary via `shutil.which("mcmctree")` (or
uses `--mcmctree-path` if supplied), then runs a `subprocess.run` with no
arguments and extracts the version from stdout via the regex
`paml version (\d+(?:\.\d+)+)`.  No `ToolEnv` or `TOOL_REGISTRY` is
involved — the detection is self-contained in ~15 lines.

### Seed handling

`run_mcmc` generates a `random.randint(1, 2**31 - 1)` seed for each run,
injects it into the ctl text via regex, and launches posterior + prior
subprocesses immediately with the same seed.  No `SeedUsed` file is read
or waited on.

When `--ctl` is provided, the user's ctl is copied into each `runN/`
after the seed line is replaced. The prior ctl is derived from the
posterior ctl (`usedata = 0`, same seed) via `_derive_prior_ctl()`.

### MCMCtree output parsing

`parse_mcmctree_out()` extracts node age estimates from `mcmctree.out`.
The regex allows optional whitespace after the opening parenthesis
(`\(\s*`) to handle PAML's `( 1.234, 5.678)` formatting.

`extract_node_tree()` returns the first Newick string after the
`Species tree for FigTree.` marker. The regex `\([\s\S]+?\)[\s\d]*;`
handles node-labeled trees (`) 7 ;`) where a bare-integer label sits
between `)` and `;`.

### Progress monitoring

A `_SampleCounter` helper tracks `(inode, byte_offset)` per `mcmc.txt`
for incremental reading — each 5-second poll reads only new bytes.

Log files are tailed to the terminal via daemon threads (unless
`--quiet`), prefixed with the run key (`run1:posterior`, etc.).

### hessian → mcmc contract

`seq_type` ("AA"|"NT") is read from `<hessian-dir>/result.json`, falling
back to scanning `iqtree.dummy.phy` content when the field is absent.
`ndata` is **always counted from `iqtree.dummy.phy`** data blocks — never
from `result.json`. The `hessian` step stores the original partition count
in `params.n_partitions` (1 for unpartitioned) for reporting only.

### result.json schema (hessian step)

```json
{
  "status": "success",
  "command": "phyloai posttree dating hessian ...",
  "wall_time": 3421.5,
  "tool_versions": {"iqtree3": "3.1.3"},
  "params": {
    "matrix": "/abs/path/matrix.fa",
    "rooted_tree": "/abs/path/calib.tre",
    "model_expr": "LG+F+R4",
    "partitions": null,
    "seq_type": "AA",
    "n_partitions": 1,
    "output_dir": "runs/posttree/dating/hessian",
    "threads": 4,
    "overwrite": false,
    "resume": false,
    "dry_run": false,
    "quiet": false
  },
  "key_results": {
    "n_partitions": 1,
    "seq_type": "AA",
    "ndata_raw": 1
  },
  "error": null,
  "data": {
    "output_files": {
      "iqtree_dummy_phy": {
        "path": "/abs/path/iqtree.dummy.phy",
        "description": "Dummy PHYLIP alignment for MCMCTree approximate likelihood calculation"
      },
      "iqtree_rooted_nwk": {
        "path": "/abs/path/iqtree.rooted.nwk",
        "description": "Rooted, calibrated tree in Newick format with fossil constraints for MCMCTree"
      },
      "iqtree_hessian": {
        "path": "/abs/path/iqtree.mcmctree.hessian",
        "description": "Gradient and Hessian matrix for approximate likelihood dating (renamed to in.BV for mcmc step)"
      }
    },
    "cmd": ["iqtree3", "-s", "/abs/path/matrix.fa", "-m", "LG+F+R4", "-te", "/abs/path/calib.tre", "--dating", "mcmctree", "--prefix", "iqtree", "-T", "4"],
    "tool_stderr": ""
  }
}
```

### result.json schema (mcmc step)

```json
{
  "status": "success",
  "command": "phyloai posttree dating mcmc ...",
  "wall_time": 3421.5,
  "tool_versions": {"mcmctree": "4.10.10"},
  "params": {
    "ctl": null, "ctl_source": "generated",
    "clock": 2, "burnin": 100000, "sample_freq": 10,
    "nsamples": 10000, "n_runs": 2,
    "seqtype": "AA", "ndata": 2, "seqtype_ndata_source": "hessian-result.json"
  },
  "key_results": {
    "n_runs": 2,
    "n_posterior_failures": 0,
    "convergence_rho_posterior": 0.998
  },
  "data": {
    "output_files": {
      "trace_run1_posterior": {
        "path": "/abs/path/diagnostics/traces/mcmc_trace_run1_posterior.pdf",
        "description": "MCMC trace plot for posterior run 1 showing parameter sampling over iterations"
      },
      "trace_run2_posterior": {
        "path": "/abs/path/diagnostics/traces/mcmc_trace_run2_posterior.pdf",
        "description": "MCMC trace plot for posterior run 2"
      },
      "trace_run1_prior": {
        "path": "/abs/path/diagnostics/traces/mcmc_trace_run1_prior.pdf",
        "description": "MCMC trace plot for prior run 1"
      },
      "trace_run2_prior": {
        "path": "/abs/path/diagnostics/traces/mcmc_trace_run2_prior.pdf",
        "description": "MCMC trace plot for prior run 2"
      },
      "convergence_posterior": {
        "path": "/abs/path/diagnostics/convergence/convergence_posterior_run1_vs_run2.pdf",
        "description": "Posterior convergence diagnostic: scatter plot with regression and y=x line for run 1 vs run 2 node ages"
      },
      "convergence_prior": {
        "path": "/abs/path/diagnostics/convergence/convergence_prior_run1_vs_run2.pdf",
        "description": "Prior convergence diagnostic: scatter plot of prior node ages between run 1 and run 2"
      },
      "infinite_sites_run1_posterior": {
        "path": "/abs/path/diagnostics/infinite_sites/infinite_sites_run1_posterior.pdf",
        "description": "Infinite-sites diagnostic: mean age vs 95% credible interval width for posterior run 1"
      },
      "infinite_sites_run2_posterior": {
        "path": "/abs/path/diagnostics/infinite_sites/infinite_sites_run2_posterior.pdf",
        "description": "Infinite-sites diagnostic for posterior run 2"
      },
      "posterior_vs_prior_run1": {
        "path": "/abs/path/diagnostics/posterior_vs_prior/posterior_vs_prior_run1.pdf",
        "description": "Posterior vs prior mean node age comparison for run 1"
      },
      "posterior_vs_prior_run2": {
        "path": "/abs/path/diagnostics/posterior_vs_prior/posterior_vs_prior_run2.pdf",
        "description": "Posterior vs prior mean node age comparison for run 2"
      },
      "convergence_summary": {
        "path": "/abs/path/diagnostics/spearman_correlations.csv",
        "description": "Spearman rank correlation and convergence statistics for posterior run comparisons"
      },
      "node_ages": {
        "path": "/abs/path/diagnostics/convergence/posterior_times.csv",
        "description": "Posterior node age estimates with 95% HPD intervals per node"
      }
    },
    "diagnostics": {"spearman": [...], "warnings": [], "generated": [...], "skipped": [...]},
    "warnings": [],
    "return_codes": {"run1:posterior": 0, "run1:prior": 0, ...}
  }
}
```
