# PhyloAI `tree bi` Design

**Date:** 2026-06-23
**Status:** Approved for implementation
**Parent spec:** `2026-06-07-phyloai-design.md`, `2026-06-17-phyloai-tree-design.md`
**JSON standard:** `2026-06-21-phyloai-json-output-standard.md`

---

## 1. Overview

`phyloai tree bi` performs Bayesian phylogenetic inference using PhyloBayes-MPI (`pb_mpi`). It runs N independent MCMC chains in parallel, monitors convergence in real time via `bpcomp` and `tracecomp`, and generates a consensus tree when chains are stopped.

Unlike `ml` and `msc`, `bi` is a long-running interactive command — chains may run for hours or days. The command stays alive throughout, providing a live progress display and periodic convergence statistics. The user terminates chains via Ctrl+C (soft-stop) or by pre-specifying a sample target.

`bi` has no subcommand layer. The command is:

```
phyloai tree bi [OPTIONS]
```

Default output directory: `runs/tree/bi/`.

---

## 2. Command-Line Interface

### 2.1 Parameter Groups

Parameters are displayed in grouped sections via `rich_click` `OPTION_GROUPS`.

#### Input / Output

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--matrix / -m` | Path | required | Input alignment (PHYLIP or FASTA). pb_mpi natively supports PHYLIP; FASTA inputs are converted automatically if needed. |
| `--output-dir / -o` | Path | `runs/tree/bi` | Output directory. |
| `--overwrite` | flag | False | Delete and recreate existing output directory before starting. |

#### Model

| Parameter | Type | Default | pb_mpi flag | Description |
|-----------|------|---------|-------------|-------------|
| `--model` | choice | `gtr` | `-gtr`, `-poisson`, `-lg`, `-wag`, `-jtt`, `-mtrev`, `-mtzoa`, `-mtart` | Relative exchangeabilities (rate matrix). |
| `--mixture` | str | `auto` | `-cat` / `-ncat N` | Profile mixture model. `auto` → `-cat` (Dirichlet process, recommended); `1` → `-ncat 1` (homogeneous, e.g. LG+G4 when `--model lg`); integer N > 1 → `-ncat N` (fixed N-component mixture). |
| `--gamma-cats` | int | `4` | `-dgam N` | Categories for discrete Gamma rate variation across sites. |
| `--start-tree` | Path | None | `-t <file>` | Starting tree in Newick format. Affects initialization only; topology is free to change. Mutually exclusive with `--fix-tree`. |
| `--fix-tree` | Path | None | `-T <file>` | Fix topology throughout MCMC. Only branch lengths and other continuous parameters are sampled. Must be a bifurcating tree. Mutually exclusive with `--start-tree`. |

**Model shorthand examples:**

| phyloai invocation | pb_mpi equivalent | IQ-TREE analogy |
|---|---|---|
| (defaults) | `-cat -gtr -dgam 4` | CAT-GTR (no IQ-TREE equivalent) |
| `--model lg --mixture 1` | `-lg -ncat 1 -dgam 4` | `LG+G4` |
| `--model poisson` | `-cat -poisson -dgam 4` | CAT-Poisson |
| `--model wag --mixture 20` | `-wag -ncat 20 -dgam 4` | `WAG+C20` |

#### Chains & Parallelism

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--chains` | int | `3` | Number of independent chains to launch, auto-named `<prefix>1`, `<prefix>2`, ... |
| `--chain-prefix` | str | `chain` | Prefix for auto-named chains. |
| `--chain-names` | str | None | Explicit chain names, comma-separated (e.g. `chain4,chain5`). Overrides `--chains` and `--chain-prefix`. Use to add chains to an existing run directory. |
| `--threads` | int | `4` | MPI processes per chain (`mpirun -np`). Minimum 2 (1 master + N-1 slaves). |

The effective names list is resolved as: if `--chain-names` is given, use it; else generate `[prefix+str(i) for i in range(1, chains+1)]`.

#### Sampling

| Parameter | Type | Default | pb_mpi flag | Description |
|-----------|------|---------|-------------|-------------|
| `--sample-freq` | int | `1` | `-x <every>` | Save one MCMC point every N cycles. |
| `--nsamples` | int | `-1` | `-x <every> <until>` | Total MCMC cycles per chain after which pb_mpi stops (pb_mpi `-x <until>`). The number of saved points is `nsamples / sample_freq` (since samples are taken every `sample_freq` cycles). `-1` = run forever (passed to pb_mpi as `<until>`; pb_mpi interprets `-1` as indefinite). |

To stop a forever-running chain: use Ctrl+C (phyloai sends soft-stop), or directly write `echo 0 > chains/<chainname>.run`.

#### Convergence Monitoring

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--monitor-freq` | int | `100` | Run `bpcomp` + `tracecomp` every N new samples (measured as: the minimum chain length across all active chains has increased by N since last check). |
| `--burnin-frac` | float | `0.5` | Fraction of saved samples to discard as burn-in **during convergence monitoring only**. Applied dynamically: `burnin = floor(min_chain_length × burnin_frac)`. Minimum burn-in of 10 samples required; checks are skipped with a warning if chains are too short. This value is NOT passed to pb_mpi (no `-x <burnin>`). Use bpcomp and tracecomp after the run to choose a final burn-in for summarisation. |
| `--poll-interval` | int | `60` | Seconds between `.trace` file polls for progress display and nsamples/check triggers. |

#### Resume

| Parameter | Type | Description |
|-----------|------|-------------|
| `--resume [CHAINS]` | optional str | Resume existing chains from their current state. No value = resume all chains listed in `run_state.json`. Comma-separated names = resume only those chains (e.g. `--resume chain1,chain3`). Resume uses the native pb_mpi mechanism: `mpirun -np <threads> pb_mpi <chainname>` (no `-d` or model flags; pb_mpi reads from `.chain` file). The `--nsamples` target from `run_state.json` is used unless overridden by a user-provided `--nsamples` value (e.g. `--resume --nsamples 10000` to extend a 5000-cycle run). Chains already at the resolved target are skipped; running chains receive a soft-stop when they reach the target. Click implementation: `@click.option('--resume', default=None, is_flag=False, flag_value='__ALL__', help='...')` — absent = `None`, bare `--resume` = `'__ALL__'`, `--resume chain1,chain2` = `'chain1,chain2'`. |

#### Tool

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--pb-path` | Path | None | Directory containing PhyloBayes tools. Overrides PATH lookup for `pb_mpi`, `bpcomp`, and `tracecomp` (required by this command) and optionally `readpb_mpi` if present. Does not require `readpb_mpi` to exist. |
| `--dry-run` | flag | False | Print all commands without executing. |
| `-q, --quiet` | flag | False | Suppress terminal output (progress display and convergence statistics). |
| `--help` | — | — | Show help and exit. |

### 2.2 Help Text

```
Usage: phyloai tree bi [OPTIONS]

  Bayesian phylogenetic inference with PhyloBayes-MPI (pb_mpi).

  Runs N independent MCMC chains in parallel using mpirun + pb_mpi, monitors
  convergence in real time via bpcomp and tracecomp, and generates a consensus
  tree. Chains run until stopped: use Ctrl+C for a safe soft-stop, or set
  --nsamples to stop automatically after N MCMC cycles.

  Examples:
    # Standard: 3 chains, CAT-GTR, 4 MPI processes each, run forever
    phyloai tree bi --matrix concat/matrix.phy

    # Homogeneous model LG+G4, stop after 10000 samples
    phyloai tree bi --matrix concat/matrix.phy --model lg --mixture 1 \
        --nsamples 10000

    # Add two extra chains to an existing run directory
    phyloai tree bi --matrix concat/matrix.phy --chain-names chain4,chain5 \
        --output-dir runs/tree/bi

    # Resume all chains in an existing directory
    phyloai tree bi --output-dir runs/tree/bi --resume

    # Resume only chain1 and chain3
    phyloai tree bi --output-dir runs/tree/bi --resume chain1,chain3

    # Resume and extend to a new nsamples target
    phyloai tree bi --output-dir runs/tree/bi --resume --nsamples 10000

    # Resume and run forever (was previously set to 5000)
    phyloai tree bi --output-dir runs/tree/bi --resume --nsamples -1
```

### 2.3 `rich_click` Option Groups

```python
click.rich_click.OPTION_GROUPS["phyloai tree bi"] = [
    {"name": "Input / Output",
     "options": ["--matrix", "--output-dir", "--overwrite"]},
    {"name": "Model",
     "options": ["--model", "--mixture", "--gamma-cats", "--start-tree", "--fix-tree"]},
    {"name": "Chains & Parallelism",
     "options": ["--chains", "--chain-prefix", "--chain-names", "--threads"]},
    {"name": "Sampling",
     "options": ["--sample-freq", "--nsamples"]},
     {"name": "Convergence Monitoring",
      "options": ["--monitor-freq", "--burnin-frac", "--poll-interval"]},
    {"name": "Resume",
     "options": ["--resume"]},
    {"name": "Tool",
     "options": ["--pb-path", "--dry-run", "--quiet", "--help"]},
]
```

---

## 3. Execution Model

### 3.1 Startup Sequence

1. Validate inputs: matrix path exists, `--start-tree` / `--fix-tree` mutual exclusion, `--mixture` value is `"auto"` or a positive integer string, `--threads >= 2`.
2. Detect tools: `pb_mpi`, `bpcomp`, `tracecomp`, `mpirun` via `ToolEnv`. `require()` raises `FileNotFoundError` (exit code 3) if any is missing. `readpb_mpi` is registered but not required at this stage.
3. Prepare output directory: if `--overwrite`, delete and recreate. If directory exists and is non-empty without `--overwrite`, raise `ValueError` (exit code 1) — unless `--resume` is set or `--chain-names` adds new chains only.
4. Create subdirectories: `output_dir/chains/`, `output_dir/convergence/`.
5. Write `output_dir/run_state.json` (fresh run only; not overwritten on resume). This file persists parameters needed for resume (see Section 3.4).
6. Build and launch all chain processes in parallel (`subprocess.Popen`), working directory = `output_dir/chains/`.
7. Enter monitoring loop.

### 3.2 Chain Commands

**Fresh run:**
```
mpirun -np <threads> pb_mpi -d <abs_matrix_path> [model_flags] -x <sample_freq> <nsamples> <chainname>
```

Model flags are assembled from CLI parameters:
- `--model gtr` → `-gtr`; `--model poisson` → `-poisson`; etc.
- `--mixture auto` → `-cat`; `--mixture 1` → `-ncat 1`; `--mixture N` → `-ncat N`
- `--gamma-cats N` → `-dgam N`
- `--start-tree <file>` → `-t <abs_path>`
- `--fix-tree <file>` → `-T <abs_path>`
- `--nsamples` is always passed: `-x <freq> <nsamples>` (pb_mpi requires both `<every>` and `<until>`; `-1` means forever)

**Resume:**
```
mpirun -np <threads> pb_mpi <chainname>
```

No `-d` or model flags. pb_mpi reads all settings from the existing `.chain` file.

### 3.3 run_state.json

`output_dir/run_state.json` records the run configuration. It is created on fresh launch and **updated (not replaced) when new chains are added via `--chain-names`**. It is never modified during `--resume`.

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

**Adding chains (`--chain-names chain4,chain5` to an existing directory):**

1. Read existing `run_state.json`.
2. Validate that the current invocation's model parameters (`model_flags`, `sample_freq`, `nsamples`, `threads`) match the stored values. If any differ, raise `ValueError` (exit code 1): `"Model parameters conflict with existing run_state.json. Use --resume to continue existing chains or choose a different --output-dir."`.
3. Validate that none of the new chain names already exist in `chain_names` (would overwrite existing chains without `--overwrite`).
4. Append new names to `chain_names` and rewrite `run_state.json`.
5. Launch only the new chains (fresh run command).

**On `--resume`:**

1. Read `run_state.json` to obtain stored `nsamples`, `chain_names`, and other parameters. If the user-provided `--nsamples` value differs from the stored value, use the user-provided value and update `run_state.json` (this enables extending a completed run, e.g. `--resume --nsamples 10000` after a 5000-cycle run). If `--nsamples` is not explicitly provided (default `-1`), the stored value is used.
2. Resolve which chains to resume: `'__ALL__'` → all names in `run_state.json`; comma-separated string → only those names (validate each exists in `run_state.json`).
3. For each chain to resume: read its `.trace` file to get current length. If `nsamples != -1` and `current_length >= nsamples`, skip that chain and print `"chain1: already at target (6420 >= 10000 samples), skipping"`.
4. Launch remaining chains with resume command.
5. If `nsamples != -1`: monitor loop issues a soft-stop for each chain once its trace length reaches `nsamples` (checked on each 60-second poll).

### 3.4 Monitoring Loop

The main thread runs a monitoring loop while chain subprocesses are alive:

**Before entering loop:** Progress bar tasks are created with pre-read trace file sample counts (so resume runs show existing data immediately, not 0). `last_trace_read` is set to `time.monotonic()` so the first poll waits a full `--poll-interval` period.

**Every `--poll-interval` seconds (default 60):** Poll each chain's `.trace` file to update the progress display. Trace length = number of complete, newline-terminated non-header lines. Return 0 if file is absent, header-only, or contains only blank/partial lines. Never return negative values.

**Every `--monitor-freq` new samples** (measured as: `min_chain_length - last_check_length >= monitor_freq`): trigger a convergence check (Section 4). Convergence output is printed via `live_display.stop()` → `console.print()` → `live_display.start()` to ensure visibility during active Rich Live rendering.

**On Ctrl+C:** capture `SIGINT`, immediately call `live_display.stop()`, write `0` to all `chains/<chainname>.run` files, print soft-stop messages, wait for all subprocesses to exit (they finish their current cycle), then run a final convergence check and write `result.json`.

**On subprocess exit:** if any chain exits with non-zero return code, mark it as failed. If all chains exit normally (zero return code or soft-stopped), write `result.json` with `status: "success"`.

---

## 4. Convergence Monitoring

### 4.1 Trigger Condition

Check fires when `min(chain_lengths) - last_check_min >= monitor_freq`.

`burnin = floor(min_chain_length × burnin_frac)`. If `burnin < 10`, skip and warn: `"Skipping convergence check: chains too short (burnin < 10)"`.

### 4.2 bpcomp Invocations

All convergence commands run with working directory `output_dir/convergence/` (not `output_dir/`). Chain files are referenced as `../chains/chain1`, trace files as `../chains/chain1.trace`, and outputs written with base name (no path prefix).

For N chains, run:
- **All chains:** `bpcomp -x <burnin> -o bpcomp_all ../chains/chain1 ../chains/chain2 ...` (cwd = convergence/)
- **All pairwise combinations:** `bpcomp -x <burnin> -o bpcomp_chain1_chain2 ../chains/chain1 ../chains/chain2` (etc.)

Output files per invocation (per pb_mpi manual section 6.2): `<basename>.bpdiff` (summary), `<basename>.bplist` (bipartition list), `<basename>.con.tre` (consensus tree). PhyloAI's parser reads `.bpdiff` for the maxdiff/meandiff summary.

Parse from `.bpdiff`: `maxdiff` and `meandiff` lines.

### 4.3 tracecomp Invocations

- **All chains:** `tracecomp -x <burnin> ../chains/chain1.trace ../chains/chain2.trace ...` (cwd = convergence/; output written to stdout, captured to `convergence/tracecomp_all.contdiff` via subprocess stdout pipe)
- **All pairwise combinations:** same pattern, output captured to `convergence/tracecomp_chain1_chain2.contdiff` (etc.)

Note: `tracecomp` takes `.trace` filenames (with extension), `bpcomp` takes chain names (without extension). Both run from `output_dir/convergence/`.

Parse from `.contdiff`: extract `effsize` (minimum across all parameters) and `rel_diff` (maximum across all parameters).

### 4.4 Convergence Thresholds

| Metric | Good | Acceptable | Not converged |
|--------|------|-----------|---------------|
| `bpcomp maxdiff` | < 0.1 | < 0.3 | ≥ 0.3 |
| `tracecomp min effsize` | > 300 | > 50 | ≤ 50 |
| `tracecomp max rel_diff` | < 0.1 | < 0.3 | ≥ 0.3 |

Per-metric status (used in display):
- `_bpcomp_status(maxdiff)`: `good` if < 0.1, `ok` if < 0.3, else `no`
- `_tracecomp_status(min_effsize, max_rel_diff)`: uses the worse of the two tracecomp metrics; same thresholds

Overall status = worst individual metric across all three. The command notifies via tiered messages (Section 5.1).

The command does **not** auto-stop chains.

### 4.5 trace_plots.pdf Generation

On each convergence check, read all `.trace` files and regenerate `convergence/trace_plots.pdf` using `matplotlib`. **All parameter columns** in the trace file are plotted — one page per column, derived dynamically from the header row (skipping `iter` and `time`). Each page shows one line per chain (distinct colors), x-axis = iteration, y-axis = parameter value. A vertical dashed line marks the current burn-in position.

The parameter set is determined at first check from the header of any available `.trace` file; all chains are expected to share the same columns.

PDF generation runs after convergence statistics are printed, wrapped in try-except to prevent convergence-check crashes. Uses `matplotlib.use("Agg")` non-interactive backend. File is overwritten on each update. Users can open it in any PDF viewer; refreshing shows the latest state.

If `matplotlib` is not installed, skip silently and print once: `"matplotlib not available; trace plots disabled."`.

---

## 5. Terminal Display

### 5.1 Layout

Uses `rich.live.Live` for the upper progress section; convergence statistics are appended below as scrolling text.

**Progress section (origin-place refreshed every 60 s):**

```
Bayesian Inference | model: CAT-GTR | chains/ | 4 threads/chain

 chain1  [============================================================]  6420 samples  12.3 s/sample
 chain2  [============================================================]  6415 samples  12.1 s/sample
 chain3  [============================================================]  6418 samples  12.2 s/sample
```

- If `--nsamples N`: bar shows `current/N` fraction (deterministic progress).
- If `--nsamples -1`: bar shows an indeterminate rolling animation; only absolute sample count and speed are shown.
- Speed unit: `s/sample` when < 60 s/sample; `min/sample` otherwise.
- On resume: sample count starts from the existing `.trace` line count, not from zero. Progress bar tasks pre-read trace files at startup; display never briefly shows 0 before the first poll.
- Rich Live rendering uses `stop()` → `print()` → `start()` + 4 blank lines buffer for convergence output to prevent Live refresh from overwriting table rows; raw `print()` bypasses Rich's Console output processing.

**Convergence statistics (appended below on each check):**

```
--- Convergence Check @ 6400 samples (burnin 50% = 3200) ---

  All chains
  bpcomp    maxdiff  0.081   meandiff  0.006   [good]
  tracecomp  min effsize  312   max rel_diff  0.094   [good]

  Pairwise
    pair              maxdiff  min effsize  max rel_diff  bpcomp  tracecomp
    chain1 x chain2   0.073       340           0.094     good       good
    chain1 x chain3   0.089       298           0.085     good       good
    chain2 x chain3   0.065       355           0.072     good       good
-------------------------------------------------------------
  *** All convergence criteria met (all pairs good). You may stop chains with Ctrl+C when ready. ***
```

The Pairwise section uses a 6-column table: `pair`, `maxdiff`, `min effsize`, `max rel_diff`, `bpcomp` (single-metric status from `_bpcomp_status`), `tracecomp` (combined status from `_tracecomp_status`). The All chains bpcomp and tracecomp lines each show their own per-column status label (not the combined overall status).

Notification messages below the separator are tiered:
- All pairs good: `"*** All convergence criteria met ..."`
- All pairs at least ok: `"Convergence acceptable across all chain pairs ... Consider stopping when ready."`
- Some pairs converged: `"Some chain pairs agree (N good, M ok, K not converged)."`
- Nothing converged: no message (table alone is sufficient)

Status labels are ASCII: `good`, `ok`, `no`. The old label `[not converged]` is replaced by per-column `no` for each tool, making it clear which metric is failing.

### 5.2 Ctrl+C Soft-Stop

The `KeyboardInterrupt` handler immediately calls `live_display.stop()` to release the console, then prints each step of the shutdown sequence:

```
^C  Caught interrupt -- sending soft-stop to all chains...
    Wrote 0 -> chains/chain1.run
    Wrote 0 -> chains/chain2.run
    Wrote 0 -> chains/chain3.run
    Waiting for chains to finish current cycle...
    chain1 stopped at 6423 samples.
    chain2 stopped at 6420 samples.
    chain3 stopped at 6421 samples.
    Running final convergence check...
    Writing result.json  (status: success)
```

---

## 6. Output Structure

```
runs/tree/bi/
├── chains/
│   ├── chain1.trace        # MCMC trace (TSV: iter time topo loglik ...)
│   ├── chain1.treelist     # Sampled trees (Newick, one per line)
│   ├── chain1.chain        # Full parameter state (binary; used by readpb_mpi)
│   ├── chain1.param        # Current parameter snapshot (text)
│   ├── chain1.monitor      # Mixing statistics (text)
│   ├── chain1.run          # Run flag: 1=running, 0=stop
│   ├── chain2.{trace,...}
│   └── chain3.{trace,...}
│
├── convergence/
│   ├── trace_plots.pdf              # All trace parameters, one page per column (updated each check)
│   ├── bpcomp_all.bpdiff            # All-chain bpcomp summary
│   ├── bpcomp_all.bplist            # All-chain bipartition list
│   ├── bpcomp_all.con.tre           # Consensus tree from all chains
│   ├── bpcomp_chain1_chain2.bpdiff  # Pairwise bpcomp outputs
│   ├── bpcomp_chain1_chain2.con.tre
│   ├── bpcomp_chain1_chain3.bpdiff
│   ├── bpcomp_chain1_chain3.con.tre
│   ├── bpcomp_chain2_chain3.bpdiff
│   ├── bpcomp_chain2_chain3.con.tre
│   ├── tracecomp_all.contdiff             # All-chain tracecomp output
│   ├── tracecomp_chain1_chain2.contdiff   # Pairwise tracecomp outputs
│   ├── tracecomp_chain1_chain3.contdiff
│   └── tracecomp_chain2_chain3.contdiff
│
└── result.json
```

---

## 7. result.json

Follows the JSON output standard (`2026-06-21-phyloai-json-output-standard.md`). `bi` uses an extended variant of the single pattern, with per-chain dicts instead of single `cmd`/`tool_stderr` scalars.

```json
{
  "status": "success | error",
  "command": "phyloai tree bi --matrix concat/matrix.phy --model gtr --mixture auto --gamma-cats 4 --chains 3 --chain-prefix chain --threads 4 --sample-freq 1 --nsamples -1 --monitor-freq 100 --burnin-frac 0.5 --poll-interval 60 --output-dir runs/tree/bi",
  "wall_time": 3600.5,
  "tool_versions": {
    "pb_mpi": null,
    "bpcomp": null,
    "tracecomp": null,
    "mpirun": "4.1.2"
  },
  "params": {
    "matrix": "concat/matrix.phy",
    "output_dir": "runs/tree/bi",
    "overwrite": false,
    "model": "gtr",
    "mixture": "auto",
    "gamma_cats": 4,
    "start_tree": null,
    "fix_tree": null,
    "chains": 3,
    "chain_prefix": "chain",
    "chain_names": null,
    "threads": 4,
    "sample_freq": 1,
    "nsamples": -1,
    "resume": null,
    "monitor_freq": 100,
    "burnin_frac": 0.5,
    "poll_interval": 60,
    "pb_path": null,
    "dry_run": false,
    "quiet": false
  },
  "key_results": {
    "chain_names": ["chain1", "chain2", "chain3"],
    "chain_lengths": {"chain1": 6423, "chain2": 6420, "chain3": 6421},
    "final_convergence": {
      "all_chains": {
        "bpcomp_maxdiff": 0.081,
        "bpcomp_meandiff": 0.006,
        "tracecomp_min_effsize": 312,
        "tracecomp_max_reldiff": 0.094,
        "status": "good"
      },
      "pairwise": {
        "chain1_chain2": {
          "bpcomp_maxdiff": 0.073,
          "tracecomp_min_effsize": 340,
          "status": "good"
        },
        "chain1_chain3": {
          "bpcomp_maxdiff": 0.089,
          "tracecomp_min_effsize": 298,
          "status": "good"
        },
        "chain2_chain3": {
          "bpcomp_maxdiff": 0.065,
          "tracecomp_min_effsize": 355,
          "status": "good"
        }
      }
    },
    "consensus_tree": "convergence/bpcomp_all.con.tre"
  },
  "error": null,
  "data": {
    "chain_cmds": {
      "chain1": ["mpirun", "-np", "4", "pb_mpi", "-d", "/abs/path/matrix.phy",
                 "-cat", "-gtr", "-dgam", "4", "-x", "1", "chain1"],
      "chain2": ["mpirun", "-np", "4", "pb_mpi", "-d", "/abs/path/matrix.phy",
                 "-cat", "-gtr", "-dgam", "4", "-x", "1", "chain2"],
      "chain3": ["mpirun", "-np", "4", "pb_mpi", "-d", "/abs/path/matrix.phy",
                 "-cat", "-gtr", "-dgam", "4", "-x", "1", "chain3"]
    },
    "tool_stderr": {
      "chain1": "",
      "chain2": "",
      "chain3": ""
    },
    "interrupted": true,
    "output_files": {
      "trace_plots": {
        "path": "/abs/path/trace_plots.pdf",
        "description": "MCMC trace plots showing parameter sampling over iterations for all chains"
      }
    },
    "warnings": []
  }
}
```

**Notes:**
- `data.chain_cmds` records the actual argv lists executed (fresh run format or resume format).
- `data.tool_stderr` contains merged stdout+stderr per chain (per JSON standard section 5.3: field name is legacy but content is merged output). pb_mpi writes diagnostics to stdout; both streams are captured and written to `chains/<chainname>.log`. The `tool_stderr` field in result.json holds this merged content (or references the log file if content is large). The JSON standard `tool_log` field may also be used: `"tool_logs": {"chain1": "chains/chain1.log", ...}`.
- `data.interrupted`: `true` if stopped via Ctrl+C before reaching `--nsamples`; `false` if chains exited normally.
- `status: "success"` for both normal exit and Ctrl+C soft-stop. `status: "error"` only if pb_mpi exits with non-zero return code.
- `tool_versions.pb_mpi/bpcomp/tracecomp`: detected by searching for version-bearing files (e.g. `pb_mpiManual*.pdf`, `VERSION`, `CHANGELOG`) in the tool's parent directory using glob + regex. Falls back to `null` if not found.

### 7.1 JSON Standard Compliance Note

`bi` does not fit the standard single-mode pattern (one `cmd` + one `tool_stderr`) because it invokes multiple parallel processes. It uses a **multi-chain variant** of the single pattern: `data.chain_cmds` (dict) and `data.tool_stderr` (dict), keyed by chain name. This is explicitly documented here as the canonical extension for multi-process commands.

---

## 8. Tool Detection & doctor Integration

### 8.1 TOOL_REGISTRY Additions (`core/env.py`)

```python
"bpcomp": {
    "required": False,
    "version_flag": "",
    "install": "https://github.com/bayesiancook/pbmpi",
},
"tracecomp": {
    "required": False,
    "version_flag": "",
    "install": "https://github.com/bayesiancook/pbmpi",
},
"readpb_mpi": {
    "required": False,
    "version_flag": "",
    "install": "https://github.com/bayesiancook/pbmpi",
},
"mpirun": {
    "required": False,
    "version_flag": "--version",
    "install": "https://www.open-mpi.org  (or: brew install open-mpi / apt install openmpi-bin)",
},
```

`mpirun` supports `--version` and returns a parseable version string. The four PhyloBayes tools have no version flag; version detection falls back to file-based heuristic (see Section 8.2).

### 8.2 Version Detection for PhyloBayes Tools

When a pb tool is found, attempt to detect version by:

1. Identifying the tool's directory (resolved executable path's parent).
2. Globbing for `pb_mpi*Manual*.pdf`, `pb_mpi*README*`, `VERSION`, `CHANGELOG`, `*.pdf` in that directory and its parent.
3. Applying regex `[Vv]ersion\s*(\d+\.\d+[\.\d]*)` or `(\d+\.\d+)` against filenames.
4. Return first match as version string. If none found, return `null`.

### 8.3 `--pb-path` Propagation

```python
if pb_path:
    tool_paths = {
        "pb_mpi":    pb_path / "pb_mpi",
        "bpcomp":    pb_path / "bpcomp",
        "tracecomp": pb_path / "tracecomp",
    }
    # readpb_mpi registered only if present; not required for tree bi
    readpb = pb_path / "readpb_mpi"
    if readpb.exists():
        tool_paths["readpb_mpi"] = readpb
else:
    tool_paths = {}
env = ToolEnv(tool_paths=tool_paths)
```

### 8.4 doctor Display

`phyloai doctor` adds a **PhyloBayes MPI** section:

```
PhyloBayes MPI
  pb_mpi      found  /opt/pbmpi/bin/pb_mpi        version: 1.9
  bpcomp      found  /opt/pbmpi/bin/bpcomp         version: 1.9
  tracecomp   found  /opt/pbmpi/bin/tracecomp      version: 1.9
  readpb_mpi  found  /opt/pbmpi/bin/readpb_mpi     version: 1.9
  mpirun      found  /usr/local/bin/mpirun          version: 4.1.2
```

---

## 9. Code Structure

### 9.1 CLI Layer (`cli/commands/tree.py`)

- Register `bi` as a `@tree.command()` (direct command, not a group).
- Add `"bi"` to `_TreeGroup.list_commands()`.
- Register `OPTION_GROUPS` for `"phyloai tree bi"`.
- Import and call `run_bi()` from `phyloai/tree/bi.py`.
- Handle exit codes: 3 = tool not found, 2 = execution error, 1 = input error.

### 9.2 Library Layer (`phyloai/tree/bi.py`)

Follows the same structural pattern as `msc.py` and `cf.py`:

```
_resolve_chain_names()
_detect_tools()               # returns pb_mpi, bpcomp, tracecomp, mpirun paths
_detect_pb_version()          # file-based heuristic
_detect_mpirun_version()      # --version parsing
_detect_tool_versions()       # combines above
_build_chain_cmd()            # builds argv list, always passes -x <freq> <nsamples>
_build_resume_cmd()           # builds resume argv list
_build_model_flags()          # assembles -cat/-ncat/-dgam/tree flags
_count_trace_samples()        # safe .trace parser
_prepare_matrix()             # FASTA→PHYLIP conversion
_state_payload()              # run_state.json dict
_write_run_state()            # persist run_state
_read_run_state()             # read run_state
_update_run_state_for_new_chains()
_resolve_resume_names()
_parse_bpcomp_bpdiff()        # parse .bpdiff → maxdiff, meandiff
_parse_tracecomp_contdiff()   # parse .contdiff → min_effsize, max_rel_diff
_status_from_metrics()        # combined 3-metric status (good/ok/not converged)
_bpcomp_status()              # per-column: maxdiff only
_tracecomp_status()           # per-column: min_effsize + max_rel_diff
_run_convergence_check()      # bpcomp + tracecomp, cwd=convergence/
_format_convergence_text()    # 6-column ASCII table + per-column status; also writes convergence_render.txt diagnostic
_generate_trace_plots()       # matplotlib PDF
_soft_stop_chains()           # writes 0 to .run files
_assemble_result()            # builds result.json dict
_run_bi_processes()           # Rich Live (stop→print→start) + convergence checks + Ctrl+C handler
run_bi()                      # main entry point
```

`run_bi()` signature:
```python
def run_bi(
    matrix: Path,
    output_dir: Path,
    overwrite: bool,
    model: str,
    mixture: str,
    gamma_cats: int,
    start_tree: Path | None,
    fix_tree: Path | None,
    chains: int,
    chain_prefix: str,
    chain_names: str | None,
    threads: int,
    sample_freq: int,
    nsamples: int,
    resume: str | None,
    monitor_freq: int,
    burnin_frac: float,
    poll_interval: int,
    pb_path: Path | None,
    dry_run: bool,
    quiet: bool,
) -> dict:
```

The `params` dict is built from these exact parameter names at the top of `run_bi()` and reused in all return paths.

---

## 10. Exit Codes

| Code | Meaning |
|------|---------|
| 0 | Success (chains completed or soft-stopped) |
| 1 | Input validation error |
| 2 | Tool execution failed (pb_mpi non-zero exit) |
| 3 | Required tool not found (pb_mpi, bpcomp, tracecomp, or mpirun) |

---

## 11. Relationship to Other Specs

- **Parent design** (`2026-06-07-phyloai-design.md`): output dir convention `runs/tree/bi/`, execution mode single, exit codes.
- **Tree design** (`2026-06-17-phyloai-tree-design.md`): `bi` is defined as owning Bayesian inference; this doc supersedes the placeholder section there.
- **JSON standard** (`2026-06-21-phyloai-json-output-standard.md`): `bi` conforms with a multi-chain extension of the single pattern.
- **doctor design** (`2026-06-18-phyloai-doctor-design.md`): add PhyloBayes MPI tool group.
- **posttree**: `readpb_mpi` is registered in `TOOL_REGISTRY` here but not required; its use will be defined in the posttree spec.
- **doctor design** (`2026-06-18-phyloai-doctor-design.md`): the registry table and doctor output section must be updated to include `bpcomp`, `tracecomp`, `readpb_mpi`, and `mpirun` in a new "PhyloBayes MPI" group. The existing `pb_mpi` entry moves into that group. Doctor tests must be updated accordingly.

### Stale References to Update During Implementation

The following references in existing specs conflict with this design and must be updated:

| Spec | Stale reference | Correct value |
|------|----------------|---------------|
| `2026-06-07-phyloai-design.md` line ~241 | output path `runs/tree/bi/phylobayes/` | `runs/tree/bi/` |
| `2026-06-07-phyloai-design.md` CLI example | `phyloai tree bi phylobayes --matrix ...` | `phyloai tree bi --matrix ...` |
| `2026-06-17-phyloai-tree-design.md` lines ~41-43 | `bi` as a Click Group with `phylobayes` subcommand | `bi` as a direct `@tree.command()` |
| `2026-06-17-phyloai-tree-design.md` lines ~61-62 | CLI hierarchy showing `bi phylobayes` | `tree bi [OPTIONS]` |
