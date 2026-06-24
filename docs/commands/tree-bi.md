# phyloai tree bi

Bayesian phylogenetic inference with [PhyloBayes-MPI](https://github.com/bayesiancook/pbmpi) (`pb_mpi`).

## Purpose

`phyloai tree bi` runs N independent MCMC chains in parallel under a PhyloBayes substitution model, monitors chain convergence in real time using `bpcomp` and `tracecomp`, and produces a consensus tree once chains are stopped.

Unlike `tree ml` and `tree msc`, `tree bi` is a long-running interactive command. Chains may run for hours or days. The command stays alive throughout, showing a live progress display and periodic convergence statistics, until the user terminates the chains with Ctrl+C (soft-stop) or the chain target is reached.

`tree bi` has no subcommand layer — it is invoked as `phyloai tree bi [OPTIONS]`.

The default output directory is `runs/tree/bi/`.

## Requirements

PhyloBayes-MPI tools must be installed and discoverable on `PATH` (or via `--pb-path`):

| Tool | Purpose |
|------|---------|
| `pb_mpi` | MCMC sampler (always required) |
| `bpcomp` | Topology convergence (always required) |
| `tracecomp` | Parameter convergence (always required) |
| `mpirun` | Open MPI launcher (always required) |
| `readpb_mpi` | Reading chain files (optional) |

Run `phyloai doctor` to confirm installation.

## Usage

```bash
phyloai tree bi --matrix <alignment> [OPTIONS]
```

The input alignment may be PHYLIP or FASTA. pb_mpi requires PHYLIP; FASTA is auto-converted before chain launch.

## Examples

```bash
# Default: 3 chains, CAT-GTR model, run forever
phyloai tree bi --matrix concat/matrix.phy

# Homogeneous LG+G4, stop after 10000 saved points
phyloai tree bi --matrix concat/matrix.phy --model lg --mixture 1 --nsamples 10000

# WAG+C20 mixture model
phyloai tree bi --matrix concat/matrix.phy --model wag --mixture 20

# Add two extra chains to an existing run
phyloai tree bi --matrix concat/matrix.phy --chain-names chain4,chain5 -o runs/tree/bi

# Resume all chains from their previous state
phyloai tree bi -o runs/tree/bi --resume

# Resume only chain1 and chain3
phyloai tree bi -o runs/tree/bi --resume chain1,chain3

# Resume and extend to a new target
phyloai tree bi -o runs/tree/bi --resume --nsamples 10000

# Resume and run forever (override previous target)
phyloai tree bi -o runs/tree/bi --resume --nsamples -1

# Custom PhyloBayes tool directory
phyloai tree bi --matrix concat/matrix.phy --pb-path /opt/pbmpi/bin

# Print the commands that would be run, then exit
phyloai tree bi --matrix concat/matrix.phy --dry-run
```

## Parameters

### Input / Output

| Flag | Type | Default | Description |
|------|------|---------|-------------|
| `--matrix / -m` | Path | required | Input alignment (PHYLIP or FASTA). FASTA is converted to PHYLIP before chain launch. Not required when `--resume` is used. |
| `--output-dir / -o` | Path | `runs/tree/bi` | Output directory. Holds `chains/`, `convergence/`, `run_state.json`, and `result.json`. |
| `--overwrite` | flag | False | Delete and recreate the output directory before starting. Mutually exclusive with `--resume`. |

### Model

| Flag | Choice | Default | pb_mpi flag | Description |
|------|--------|---------|-------------|-------------|
| `--model` | `gtr`, `poisson`, `lg`, `wag`, `jtt`, `mtrev`, `mtzoa`, `mtart` | `gtr` | `-gtr`, `-poisson`, … | Rate matrix. |
| `--mixture` | str | `auto` | `-cat` / `-ncat N` | `auto` = CAT Dirichlet process; `1` = homogeneous (e.g. LG+G4); integer N > 1 = fixed N-component mixture. |
| `--gamma-cats` | int ≥ 1 | 4 | `-dgam N` | Discrete Gamma rate categories. |
| `--start-tree` | Path | None | `-t <file>` | Starting tree (Newick). Topology is free to change. Mutually exclusive with `--fix-tree`. |
| `--fix-tree` | Path | None | `-T <file>` | Fixed topology (Newick). Only branch lengths and other parameters are sampled. Must be bifurcating. Mutually exclusive with `--start-tree`. |

**Shorthand:**

| PhyloAI invocation | pb_mpi equivalent | IQ-TREE analogy |
|---|---|---|
| (defaults) | `-cat -gtr -dgam 4` | CAT-GTR |
| `--model lg --mixture 1` | `-lg -ncat 1 -dgam 4` | LG+G4 |
| `--model poisson` | `-cat -poisson -dgam 4` | CAT-Poisson |
| `--model wag --mixture 20` | `-wag -ncat 20 -dgam 4` | WAG+C20 |

### Chains & Parallelism

| Flag | Type | Default | Description |
|------|------|---------|-------------|
| `--chains` | int ≥ 1 | 3 | Number of independent chains; auto-named `<prefix>1`, `<prefix>2`, … |
| `--chain-prefix` | str | `chain` | Prefix for auto-named chains. |
| `--chain-names` | str | None | Comma-separated names (e.g. `chain4,chain5`). Overrides `--chains` and `--chain-prefix`. Use to add chains to an existing run. |
| `--threads / -t` | int ≥ 2 | 4 | MPI processes per chain (`mpirun -np`). Minimum 2 (1 master + N-1 slaves). |

The effective names list is: if `--chain-names` is given, use it; else generate `[prefix+str(i) for i in range(1, chains+1)]`.

### Sampling

| Flag | Type | Default | Description |
|------|------|---------|-------------|
| `--sample-freq` | int ≥ 1 | 1 | Save one MCMC point every N cycles (pb_mpi `-x <every>`). |
| `--nsamples` | int | `-1` | Stop after N MCMC cycles per chain (pb_mpi `-x <until>`). `-1` = run forever; stop with Ctrl+C. The number of saved points is `nsamples / sample-freq`. |

To stop a forever-running chain: use Ctrl+C, or write `echo 0 > chains/<chainname>.run`.

### Convergence Monitoring

| Flag | Type | Default | Description |
|------|------|---------|-------------|
| `--monitor-freq` | int ≥ 1 | 100 | Run `bpcomp` + `tracecomp` every N new samples (min chain length increases by N). |
| `--burnin-frac` | float `[0.0, 1.0)` | 0.5 | Fraction of saved samples discarded as burn-in **during convergence monitoring only**. Applied dynamically: `burnin = floor(min_chain_length × burnin_frac)`. Minimum burn-in 10 samples; checks are skipped with a warning if chains are too short. NOT passed to pb_mpi. |
| `--poll-interval` | int ≥ 1 | 60 | Seconds between `.trace` file polls for progress display and convergence triggers. |

### Resume

| Flag | Type | Description |
|------|------|-------------|
| `--resume [CHAINS]` | optional str | Resume from `run_state.json`. No value = resume all chains. Comma-separated names = resume only those chains (e.g. `--resume chain1,chain3`). Mutually exclusive with `--overwrite`. |

Resume uses the native pb_mpi mechanism: `mpirun -np <threads> pb_mpi <chainname>` (no `-d` or model flags; pb_mpi reads the existing `.chain` file). When `--resume` is used together with `--nsamples`, the new value overrides the stored target in `run_state.json` and chains continue from their current state. Chains already at the resolved target are skipped. Running chains receive a soft-stop when they reach the target.

Click implementation: `@click.option('--resume', default=None, is_flag=False, flag_value='__ALL__', help='...')` — absent = `None`, bare `--resume` = `'__ALL__'`, `--resume chain1,chain2` = `'chain1,chain2'`.

### Tool

| Flag | Type | Default | Description |
|------|------|---------|-------------|
| `--pb-path` | Path | None | Directory containing `pb_mpi`, `bpcomp`, `tracecomp`. Overrides PATH lookup. `readpb_mpi` is also detected here if present, but is not required. |
| `--dry-run` | flag | False | Print all commands without executing. |
| `-q, --quiet` | flag | False | Suppress terminal output (progress display and convergence statistics). |

## How It Works

### Startup

1. Validate inputs (matrix exists, `--start-tree`/`--fix-tree` mutual exclusion, `--mixture` value, `--threads ≥ 2`).
2. Detect tools (`pb_mpi`, `bpcomp`, `tracecomp`, `mpirun`) via `ToolEnv`. Missing required tools → exit code 3.
3. Prepare output directory: delete and recreate if `--overwrite`; reject non-empty directory without `--overwrite` unless `--resume` or `--chain-names` adds new chains.
4. Create `chains/` and `convergence/` subdirectories.
5. Write `run_state.json` (fresh run only; updated, not replaced, when new chains are added).
6. Build and launch all chain processes (`subprocess.Popen`), working directory `chains/`.
7. Enter the monitoring loop.

### Chain Commands

Fresh run:

```
mpirun -np <threads> pb_mpi -d <abs_matrix_path> [model_flags] -x <sample_freq> <nsamples> <chainname>
```

Resume:

```
mpirun -np <threads> pb_mpi <chainname>
```

No `-d` or model flags; pb_mpi reads all settings from the existing `.chain` file.

### Monitoring Loop

A single main thread polls `.trace` files and triggers convergence checks while chain subprocesses run:

- **Every `--poll-interval` seconds**: poll each chain's `.trace` file, update the progress display, and check whether `--nsamples` has been reached (writing `0` to `chains/<chain>.run` for chains that have).
- **Every `--monitor-freq` new samples**: trigger a `bpcomp` + `tracecomp` check (Section "Convergence Monitoring"). Convergence output is rendered via `live_display.stop() → print() → live_display.start()` so the progress bar is preserved.
- **On Ctrl+C**: write `0` to all `chains/<chainname>.run` files, wait for subprocesses to finish their current cycle, then run a final convergence check and write `result.json`.
- **On subprocess non-zero exit**: stop the remaining chains and write `result.json` with `status: "error"`.

Progress bar tasks are pre-initialised with each chain's existing trace sample count, so resume runs display correct values immediately (not 0).

## Convergence Monitoring

### Trigger and Burn-in

A check fires when `min(chain_lengths) - last_check_min ≥ monitor-freq`. The dynamic burn-in is `floor(min_chain_length × burnin_frac)`. If `burnin < 10`, the check is skipped and a warning `"Skipping convergence check: chains too short (burnin < 10)"` is emitted. On `--resume`, the new `--nsamples` (if provided) becomes the chain target; chains already at or beyond the target are skipped.

### bpcomp / tracecomp Invocations

For N chains, PhyloAI runs:

- **All chains**: `bpcomp -x <burnin> -o bpcomp_all ../chains/chain1 ../chains/chain2 …` (cwd `convergence/`)
- **All pairwise combinations**: `bpcomp -x <burnin> -o bpcomp_chain1_chain2 ../chains/chain1 ../chains/chain2`, etc.
- **tracecomp** with the same structure, taking `.trace` files and writing `.contdiff` outputs.

Both commands run with working directory `convergence/`, so chain files are referenced as `../chains/<chain>` and trace files as `../chains/<chain>.trace`.

### Thresholds

| Metric | Good | Acceptable | Not converged |
|--------|------|------------|---------------|
| `bpcomp maxdiff` | < 0.1 | < 0.3 | ≥ 0.3 |
| `tracecomp min effsize` | > 300 | > 50 | ≤ 50 |
| `tracecomp max rel_diff` | < 0.1 | < 0.3 | ≥ 0.3 |

Per-column status: `bpcomp` uses `maxdiff` alone; `tracecomp` uses the worse of `min effsize` and `max rel_diff`. Overall status = worst metric across all three. The command does **not** auto-stop chains.

### Terminal Display

```
  All chains
  bpcomp    maxdiff  0.081   meandiff  0.006   [good]
  tracecomp  min effsize  312   max rel_diff  0.094   [good]

  Pairwise
    pair              maxdiff  min effsize  max rel_diff  bpcomp  tracecomp
    chain1 x chain2   0.073       340           0.094     good       good
    chain1 x chain3   0.432        76           0.210       no         ok
    chain2 x chain3   0.065       355           0.072     good       good
```

Tiered notifications below the table:

- All pairs `good` → `"*** All convergence criteria met (all pairs good). You may stop chains with Ctrl+C when ready. ***"`
- All pairs at least `ok` → `"Convergence acceptable across all chain pairs (N good, M ok). Consider stopping when ready."`
- Some pairs converged → `"Some chain pairs agree (N good, M ok, K not converged)."`
- Nothing converged → no notification (the table alone is sufficient).

### Trace Plots

On each convergence check, PhyloAI regenerates `convergence/trace_plots.pdf` using matplotlib. One page per trace parameter column (all columns except `iter` and `time`); one line per chain. A vertical dashed line marks the current burn-in position. If matplotlib is not installed, the PDF generation is silently skipped and a one-time `"matplotlib not available; trace plots disabled."` is printed.

## Resume Semantics

`run_state.json` is the source of truth for resume. It is created on fresh launch and updated (not replaced) when new chains are added via `--chain-names`. It is updated again when `--resume` is used with a new `--nsamples`. The file never exists before a fresh launch.

Schema:

```json
{
  "chain_names": ["chain1", "chain2", "chain3"],
  "matrix": "/abs/path/matrix.phy",
  "model_flags": ["-cat", "-gtr", "-dgam", "4"],
  "sample_freq": 1,
  "nsamples": 10000,
  "threads": 4
}
```

**Adding chains (`--chain-names chain4,chain5` to an existing directory):** PhyloAI validates that the new invocation's model parameters (`model_flags`, `sample_freq`, `nsamples`, `threads`) match the stored values. If any differ, it exits with code 1 and the message `"Model parameters conflict with existing run_state.json. Use --resume to continue existing chains or choose a different --output-dir."`. It also validates that none of the new names already exist.

**Resuming:**

1. Read `run_state.json` to obtain stored `nsamples` and `chain_names`.
2. If the user provided `--nsamples` and it differs from the stored value, use the new value and update `run_state.json` (this enables extending a completed run).
3. Resolve which chains to resume: `'__ALL__'` → all names; comma-separated → only those names.
4. For each chain to resume: read its `.trace`. If `nsamples != -1` and `current_length ≥ nsamples`, skip that chain.
5. Launch remaining chains with the resume command. Monitor loop issues a soft-stop when each chain reaches the target.

The previous `result.json` is automatically backed up with a timestamp (e.g. `result_20260624_134500.json`) before the new result is written. Multiple resume cycles each produce a distinct backup.

## Safe Stopping

Use **Ctrl+C**. PhyloAI writes `0` to each `chains/<chain>.run` file and waits for pb_mpi to finish its current cycle. Direct interruption of pb_mpi (e.g. `kill -9`) can leave incomplete samples in the trace file. If chains are not yet at the target, `data.interrupted` is set to `true` in `result.json`.

## Outputs

```
runs/tree/bi/
├── chains/
│   ├── chain1.trace        # MCMC trace (TSV: iter time topo loglik ...)
│   ├── chain1.treelist     # Sampled trees (Newick, one per line)
│   ├── chain1.chain        # Full parameter state (binary; used by readpb_mpi)
│   ├── chain1.param        # Current parameter snapshot (text)
│   ├── chain1.monitor      # Mixing statistics (text)
│   ├── chain1.run          # Run flag: 1=running, 0=stop
│   ├── chain1.log          # merged stdout+stderr (PhyloAI written)
│   ├── chain2.{trace,...}
│   └── chain3.{trace,...}
│
├── convergence/
│   ├── trace_plots.pdf              # All trace parameters, one page per column
│   ├── bpcomp_all.bpdiff            # All-chain bpcomp summary (parsed)
│   ├── bpcomp_all.bplist            # All-chain bipartition list
│   ├── bpcomp_all.con.tre           # Consensus tree from all chains
│   ├── bpcomp_chain1_chain2.bpdiff  # Pairwise bpcomp outputs
│   ├── bpcomp_chain1_chain2.con.tre
│   ├── ...
│   ├── tracecomp_all.contdiff             # All-chain tracecomp output
│   ├── tracecomp_chain1_chain2.contdiff
│   └── ...
│
├── run_state.json           # Resume metadata
└── result.json              # Structured PhyloAI result
```

`result.json` follows the standard schema (status / command / wall_time / tool_versions / params / key_results / error / data). The BI command uses a **multi-chain extension** of the single pattern:

- `data.chain_cmds` is a dict keyed by chain name; each value is the actual argv list executed (fresh or resume format).
- `data.tool_stderr` is a dict keyed by chain name; each value is the merged stdout+stderr content from `chains/<chain>.log`.
- `data.tool_logs` references the per-chain log file paths.
- `data.interrupted` is `true` if the run ended via Ctrl+C before reaching `--nsamples`.
- `status` is `"success"` for both normal exit and Ctrl+C soft-stop; `"error"` only if a pb_mpi chain exits with a non-zero return code.
- `key_results.consensus_tree` is set to `convergence/bpcomp_all.con.tre` only when a `bpcomp` run actually produced that file (i.e. convergence check ran and succeeded). It is `null` when checks were skipped, failed, or chains were too short.
- `tool_versions.pb_mpi`, `bpcomp`, `tracecomp` are detected from filenames (e.g. `pb_mpiManual1.9.pdf`, `VERSION`) when present, otherwise `null`. `mpirun` is detected from `--version`.

## Exit Codes

| Code | Meaning |
|------|---------|
| 0 | Success (chains completed or soft-stopped) |
| 1 | Input validation error (missing matrix, conflicting chain names, invalid parameter combination, output-directory conflict) |
| 2 | Tool execution failed (a pb_mpi chain exited with a non-zero return code) |
| 3 | Required tool not found (pb_mpi, bpcomp, tracecomp, or mpirun) |

Validation, tool-detection, and input-conflict failures always write a `result.json` with `status: "error"` and the error message.

## Tips and Warnings

- **Let the chain warm up.** Convergence checks are only triggered when `burnin = floor(min_chain_length × burnin_frac) ≥ 10`. Short chains skip checks with a warning.
- **Adjust `--burnin-frac` to your model.** Larger values are safer for slow-mixing models but delay the first check.
- **Use `--threads ≥ 2`.** pb_mpi requires 1 master + N-1 slaves; the CLI rejects `--threads 1`.
- **Add chains without re-running the whole set.** `--chain-names chain4,chain5` (with the same output dir) appends new chains to the existing `run_state.json` and launches them as fresh chains. Model parameters must match the stored values.
- **Extend a completed run.** Use `--resume --nsamples <new_target>` to continue from where chains stopped. Chains already at the new target are skipped silently.
- **Visualise mixing.** `convergence/trace_plots.pdf` is regenerated on every check. Open it in any PDF viewer to inspect chain mixing.
- **Slow runs.** Larger `--poll-interval` reduces I/O on network filesystems; larger `--monitor-freq` reduces `bpcomp`/`tracecomp` overhead.
- **Set up a screen / tmux session.** Long-running chains should be launched from a session manager so Ctrl+C and resume work correctly.
