# phyloai posttree dating

## Purpose

Two-step Bayesian molecular dating with MCMCtree using approximate
likelihood (usedata=2):

1. **`hessian`** — IQ-TREE3 `--dating mcmctree` to compute gradients and Hessian.
2. **`mcmc`** — MCMCtree with the IQ-TREE output, running independent
   posterior + prior chains in parallel with real-time progress and
   diagnostic plots.

## Usage

```bash
# hessian — compute gradients and Hessian
phyloai posttree dating hessian --matrix concat.aa.fa --rooted-tree calib.tre
phyloai posttree dating hessian --matrix concat.aa.fa --rooted-tree calib.tre --model-expr C10+F+G4
phyloai posttree dating hessian --matrix concat.aa.fa --rooted-tree calib.tre --partitions partitions.txt

# mcmc — run MCMCtree
phyloai posttree dating mcmc --hessian-dir runs/dating/hessian
phyloai posttree dating mcmc --hessian-dir runs/dating/hessian --clock 3 --burnin 200000 --nsamples 20000
phyloai posttree dating mcmc --hessian-dir runs/dating/hessian --ctl my_run.ctl
phyloai posttree dating mcmc --hessian-dir runs/dating/hessian --dry-run
```

## Inputs

### hessian

| Input | Description |
|-------|-------------|
| `--matrix` | Alignment: FASTA (`.fa`/`.fas`/`.fasta`/`.faa`/`.fna`/`.aln`), PHYLIP, or NEXUS. `--seq-type auto` reads all three formats via PhyloAI's shared format detector. |
| `--rooted-tree` | MCMCtree calibration tree with fossil constraints on internal nodes and a mandatory root age constraint, e.g. `(A,((B,(C,D)'>3.1<3.8'),(E,F)'>2.9<3.6'))'<4.2';`. Units: 100 Mya. |
| `--model-expr` | IQ-TREE model string (e.g. `C10+F+G4`). Mutually exclusive with `--partitions`. |
| `--partitions` | Partition file (RAxML-like, NEXUS `.best_model.nex`, or cluster file). `< 10` partitions run directly; `>= 10` are auto-merged with `--merge --rclusterf 10`. |
| `--seq-type` | `AA`, `NT`, or `auto` (default). Default models: `LG+F+G4` (AA), `GTR+G4` (NT). |

### mcmc

| Input | Default | Description |
|-------|---------|-------------|
| `--hessian-dir` | *(required)* | Directory containing `iqtree.dummy.phy`, `iqtree.rooted.nwk`, `iqtree.mcmctree.hessian`. |
| `-o`, `--output` | `runs/posttree/dating/mcmc` | Output directory. |
| `--runs` | `2` | Number of independent posterior runs. |
| `--clock` | `2` | `1` = global, `2` = independent rates, `3` = correlated rates. |
| `--burnin` | `100000` | MCMC burn-in iterations. |
| `--sample-freq` | `10` | Sampling frequency (iterations per sample). |
| `--nsamples` | `10000` | Samples kept post-burnin. |
| `--ctl` | — | Pre-configured `mcmctree.ctl`. Mutually exclusive with non-default `--clock`/`--burnin`/`--sample-freq`/`--nsamples`. PhyloAI injects a random seed per run. |
| `--mcmctree-path` | — | Override auto-detected mcmctree binary. |
| `--dry-run` | — | Print generated ctl and exit. |
| `--overwrite` | — | Delete existing output directory before running. |
| `--quiet` | — | Suppress MCMC log tail output. |

Total iterations = `--burnin` + (`--sample-freq` × `--nsamples`).

## Outputs

### hessian

```
<output-dir>/
├── result.json
├── iqtree.dummy.phy
├── iqtree.rooted.nwk
└── iqtree.mcmctree.hessian
```

### mcmc

```
<output-dir>/
├── result.json
├── mcmctree.ctl                      # generated template (or copy of --ctl)
├── run1/
│   ├── mcmctree.ctl                  # posterior ctl, injected seed
│   ├── mcmc.txt                      # MCMC parameter trace
│   ├── mcmctree.out                  # divergence time summary
│   ├── mcmctree.log
│   ├── FigTree.tre                   # annotated dated tree
│   ├── FigTree.node.tre              # node-labeled tree
│   ├── iqtree.dummy.phy -> <hessian-dir>/
│   ├── iqtree.rooted.nwk -> <hessian-dir>/
│   ├── in.BV -> <hessian-dir>/iqtree.mcmctree.hessian
│   └── prior/
│       ├── mcmctree.ctl              # usedata=0, same seed
│       ├── mcmc.txt
│       └── mcmctree.out
├── run2/                             # identical structure
│   └── ...
└── diagnostics/
    ├── traces/                       # per-parameter PDFs
    ├── convergence/
    │   ├── posterior_times.csv       # all runs combined
    │   ├── prior_times.csv
    │   └── convergence_*_runX_vs_runY.pdf
    ├── infinite_sites/
    ├── posterior_vs_prior/
    └── spearman_correlations.csv
```

## Examples

```bash
# 1. Compute Hessian (unpartitioned AA)
phyloai posttree dating hessian --matrix concat.aa.fa --rooted-tree calib.tre -o runs/dating/hessian

# 2. Run MCMCtree (default 2 runs)
phyloai posttree dating mcmc --hessian-dir runs/dating/hessian -o runs/dating/mcmc

# 3. Three runs, correlated clock
phyloai posttree dating mcmc --hessian-dir runs/dating/hessian --runs 3 --clock 3

# 4. Dry-run to inspect the generated ctl
phyloai posttree dating mcmc --hessian-dir runs/dating/hessian --dry-run

# 5. Use a pre-configured mcmctree.ctl
phyloai posttree dating mcmc --hessian-dir runs/dating/hessian --ctl my_run.ctl

# 6. Partitioned Hessian
phyloai posttree dating hessian --matrix concat.aa.fa --rooted-tree calib.tre --partitions partitions.nex -o runs/dating/hessian

# 7. NT data with explicit seq-type
phyloai posttree dating hessian --matrix concat.nt.fa --rooted-tree calib.tre --seq-type NT -o runs/dating/hessian
```

## Warnings / Errors

| Condition | Behaviour |
|-----------|-----------|
| Empty or missing `hessian-dir` | Exit 1 with message listing which files are absent. |
| `--matrix` or `--rooted-tree` not found | Exit 1 before IQ-TREE starts. |
| IQ-TREE returns non-zero | Status `error` in `result.json`; stderr captured in `data.warnings`. |
| `--ctl` with non-default `--clock`/`--burnin`/`--sample-freq`/`--nsamples` | Exit 1; these flags only apply when generating the ctl from scratch. |
| Non-empty output directory without `--overwrite` | Exit 1; use `--overwrite` to replace. |
| mcmctree binary not found | Exit 3 (env error). |
| Posterior exit code ≠ 0 | Top-level `status` becomes `"error"`. |
| Prior exit code ≠ 0 or missing output files | Recorded in `data.warnings`; posterior results and diagnostics still produced. |
| `mcmc.txt` or `mcmctree.out` empty | Recorded in `data.warnings`; diagnostics degraded gracefully. |
| `--runs=1` | Convergence plots skipped; other diagnostics still generated. |

## Notes

- **Approximate likelihood only (usedata=2).** `usedata=1` (exact likelihood
  via sequence data) and `usedata=3` (output in.BV) are not implemented in
  PhyloAI.
- **Seed injection.** PhyloAI generates a unique random seed per run
  (`random.randint(1, 2³¹-1)`), injected into the ctl before launch.
  Posterior and prior share the same seed within a run.
- **ndata.** The number of data blocks is always counted directly from
  `iqtree.dummy.phy` — never from `result.json`. This ensures correctness
  when IQ-TREE merges partitions (`--merge --rclusterf 10`).
- **seqtype.** Read from `hessian/result.json` (`params.seq_type`),
  falling back to auto-detection from dummy.phy content.
- **Version detection.** The mcmctree binary is found via `shutil.which`;
  the version is extracted by running the binary with no arguments and
  matching `paml version (\d+(?:\.\d+)+)`.
- **Custom ctl (`--ctl`).** PhyloAI always symlinks the standard hessian
  files (`iqtree.dummy.phy`, `iqtree.rooted.nwk`, `in.BV`) into each `runN/`
  directory. Custom relative `seqfile`/`treefile` paths in the ctl are
  resolved relative to the ctl file's directory and symlinked into each
  `runN/` as well.
- **Diagnostics.** Convergence plots require `--runs >= 2`. All plots use a
  single dashed fit line with equation; node labels show `nXX` format.
  Posterior-vs-prior plots use equal-length axes.
- **`result.json`** follows the PhyloAI JSON Output Standard.
