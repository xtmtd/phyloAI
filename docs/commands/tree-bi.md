# phyloai tree bi

Bayesian phylogenetic inference with PhyloBayes-MPI.

## Usage

```bash
phyloai tree bi --matrix concat/matrix.phy
```

## Common Examples

```bash
# Default: 3 chains, CAT-GTR model, run forever
phyloai tree bi --matrix concat/matrix.phy

# Homogeneous LG+G4, auto-stop after 11000 cycles
phyloai tree bi --matrix concat/matrix.phy --model lg --mixture 1 --nsamples 11000

# Add extra chains to an existing run
phyloai tree bi --matrix concat/matrix.phy --chain-names chain4,chain5 -o runs/tree/bi

# Resume all chains
phyloai tree bi -o runs/tree/bi --resume

# Resume selected chains
phyloai tree bi -o runs/tree/bi --resume chain1,chain3

# Resume and extend to a higher nsamples
phyloai tree bi -o runs/tree/bi --resume --nsamples 10000

# Resume and run forever (override previous target)
phyloai tree bi -o runs/tree/bi --resume --nsamples -1
```

## Burn-in and Convergence

`--burnin-frac` (default 0.5) controls the fraction of saved samples discarded for **convergence monitoring only**. It is not passed to pb_mpi. After the run, inspect chains with bpcomp and tracecomp to determine a suitable final burn-in:

```bash
bpcomp -x 5000 chains/chain1 chains/chain2 chains/chain3
tracecomp -x 5000 chains/chain1.trace chains/chain2.trace chains/chain3.trace
```

## Sampling

`--nsamples` sets the total MCMC cycles (pb_mpi `-x <until>`). The number of saved points is `nsamples / sample-freq`. Use `-1` to run forever (stop with Ctrl+C).

## Monitoring

`--poll-interval` (default 60 s) controls how often trace files are read to update progress and trigger convergence checks. `--monitor-freq` controls how many new samples must accumulate before the next bpcomp+tracecomp check. Progress bars pre-read existing trace files at startup; resume runs show correct counts immediately, not 0.

Convergence statistics are displayed as a 6-column table:

```
  All chains
  bpcomp    maxdiff  0.081   meandiff  0.006   [good]
  tracecomp  min effsize  312   max rel_diff  0.094   [good]

  Pairwise
    pair              maxdiff  min effsize  max rel_diff  bpcomp  tracecomp
    chain1 x chain2   0.073       340           0.094     good       good
    chain1 x chain3   0.432        76           0.210       no         ok
```

Per-column status: `bpcomp` uses `maxdiff` alone; `tracecomp` uses the worse of `min effsize` and `max rel_diff`. Thresholds: `< 0.1 good, < 0.3 ok, else no`. The `tracecomp` column in the pairwise table uses both `min effsize` and `max rel_diff` thresholds combined.

Tiered notifications below the table:
- All pairs `good` → `"*** All convergence criteria met ..."`
- All pairs at least `ok` → `"Convergence acceptable ... Consider stopping when ready."`
- Some converged → `"Some chain pairs agree (N good, M ok, K not converged)."`

## Resume

When `--resume` is used, the previous `result.json` is automatically backed up with a timestamp (e.g. `result_20260624_134500.json`) before the new result is written. Multiple resume cycles each produce a distinct backup.

Chains that have already reached their `--nsamples` target are skipped. Use `--nsamples` with `--resume` to extend a run: the new value overrides the stored target in `run_state.json` and chains continue from their current state. For example, after a 5000-cycle run, `--resume --nsamples 10000` resumes all chains and stops when they reach 10000 cycles.

## Safe Stopping

Use Ctrl+C. PhyloAI writes `0` to each `chains/<chain>.run` file and waits for pb_mpi to finish its current cycle. Direct interruption of pb_mpi can leave incomplete samples.

## Outputs

- `chains/<chain>.trace`: MCMC trace.
- `chains/<chain>.treelist`: sampled trees.
- `chains/<chain>.chain`: saved chain state.
- `chains/<chain>.log`: merged stdout and stderr.
- `convergence/bpcomp_all.bpdiff`: bpcomp summary parsed by PhyloAI.
- `convergence/tracecomp_all.contdiff`: tracecomp summary.
- `convergence/trace_plots.pdf`: trace plots (requires matplotlib; silently skipped if unavailable).
- `run_state.json`: resume metadata.
- `result.json`: structured PhyloAI result.
