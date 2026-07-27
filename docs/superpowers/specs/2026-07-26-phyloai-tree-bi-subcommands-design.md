# PhyloAI `tree bi` Subcommands Design

**Date:** 2026-07-26
**Status:** Approved for implementation
**Parent spec:** `2026-06-23-phyloai-tree-bi-design.md`, `2026-06-07-phyloai-design.md`
**JSON standard:** `2026-06-21-phyloai-json-output-standard.md`

---

## 1. Overview

This spec extends `tree bi` from a single flat command into a Click Group with four subcommands:

| Subcommand | Purpose |
|---|---|
| `tree bi pb` | MCMC chain inference (exact rename of the existing `tree bi` command) |
| `tree bi bpcomp` | Final topology convergence analysis with user-specified burn-in |
| `tree bi tracecomp` | Final parameter convergence analysis with user-specified burn-in |
| `tree bi readpb` | Posterior analysis via `readpb_mpi` with output conversion |

**Breaking change:** `phyloai tree bi [OPTIONS]` becomes `phyloai tree bi pb [OPTIONS]`. All existing scripts and documentation must be updated accordingly.

The motivation for the split: `bi pb` runs `bpcomp`/`tracecomp` with a rolling 50% burn-in for live monitoring. `bi bpcomp` and `bi tracecomp` are for final analysis after chains are complete, where the user supplies a carefully chosen integer burn-in. `bi readpb` provides posterior summaries and predictive checks not available during live monitoring.

---

## 2. CLI Structure

### 2.1 Click Hierarchy

```
phyloai tree
└── bi                   ← _BiGroup (Click Group)
    ├── pb               ← direct @bi.command(); original run_bi() logic
    ├── bpcomp           ← direct @bi.command(); run_bi_bpcomp()
    ├── tracecomp        ← direct @bi.command(); run_bi_tracecomp()
    └── readpb           ← direct @bi.command(); run_bi_readpb()
```

### 2.2 Changes to Existing Code

**`phyloai/cli/commands/tree.py`:**

- `_TreeGroup.list_commands()` returns `["ml", "bi", "msc", "cf"]` — unchanged.
- Convert `@tree.command("bi", ...)` to `@tree.group("bi", cls=_BiGroup, ...)`.
- Add `class _BiGroup(click.Group)` with `list_commands()` returning `["pb", "bpcomp", "tracecomp", "readpb"]`.
- Register the original `bi_command` as `@bi.command("pb", ...)`.
- Rename internal reference from `run_bi()` to `run_bi_pb()`.
- Register `OPTION_GROUPS` for `"phyloai tree bi pb"` (was `"phyloai tree bi"`).
- Add three new subcommand registrations (Sections 3–5).

**`phyloai/tree/bi.py`:**

- Rename `run_bi()` → `run_bi_pb()` (signature unchanged).
- All shared helpers remain in `bi.py` — no extraction to a separate file.

### 2.3 New Library Files

```
phyloai/tree/
├── bi.py           # existing; run_bi() renamed to run_bi_pb(); shared helpers stay here
├── bi_bpcomp.py    # new; run_bi_bpcomp()
├── bi_tracecomp.py # new; run_bi_tracecomp()
└── bi_readpb.py    # new; run_bi_readpb()
```

`bi_bpcomp.py` and `bi_tracecomp.py` import shared helpers from `bi.py`:

```python
from phyloai.tree.bi import (
    _parse_bpcomp_bpdiff,
    _parse_tracecomp_contdiff,
    _bpcomp_status,
    _tracecomp_status,
)
```

Each new module resolves only its required executable(s) with `ToolEnv`; it does not reuse the live-monitoring `_detect_tools()` helper.

### 2.4 MCP Tool Name Change

| Old | New |
|---|---|
| `phyloai_tree_bi` | `phyloai_tree_bi_pb` |
| — | `phyloai_tree_bi_bpcomp` (new) |
| — | `phyloai_tree_bi_tracecomp` (new) |
| — | `phyloai_tree_bi_readpb` (new) |

---

## 3. `tree bi bpcomp`

### 3.1 Purpose

Run `bpcomp` once with a user-specified integer burn-in on a completed or running chains directory. Produces a final consensus tree and bipartition statistics.

### 3.2 Command

```
phyloai tree bi bpcomp [OPTIONS]
```

### 3.3 Parameters

#### Input / Output

| Parameter | Type | Default | Description |
|---|---|---|---|
| `--chain-dir` | Path | **required** | Directory containing chain files (`<chain>.chain`, `<chain>.treelist`, etc.). |
| `--chain-names` | str | `all` | Comma-separated chain names to include (e.g. `chain1,chain2`). `all` = all chains found in `--chain-dir` (sorted lexicographically by filename stem). |
| `--output-dir` | Path | `runs/tree/bi/bpcomp` | Output directory for result files and `result.json`. |
| `--overwrite` | flag | False | Delete and recreate the output directory. |

#### Analysis

| Parameter | Type | Default | pb flag | Description |
|---|---|---|---|---|
| `--burnin` | int ≥ 0 | `0` | `-x <burn-in>` | Number of saved samples to discard as burn-in. `0` = no burn-in. |
| `--sample-freq` | int ≥ 1 | `1` | `-x <burn-in> <every>` | Sub-sampling frequency: take one tree every N saved samples after burn-in. |
| `--until` | str | `all` | `-x <burn-in> <every> <until>` | Stop at this sample index. `all` = use the entire chain. Integer = stop at that saved sample index. |
| `--cutoff` | float (0,1) | `0.5` | `-c <cutoff>` | Majority-rule consensus cutoff: nodes with posterior probability below this are collapsed. |

**bpcomp `-x` flag construction:**
- `--burnin 1000 --sample-freq 1 --until all` → `-x 1000`
- `--burnin 1000 --sample-freq 10 --until all` → `-x 1000 10`
- `--burnin 1000 --sample-freq 10 --until 5000` → `-x 1000 10 5000`
- Only append `<every>` when `sample-freq != 1` **or** `until != all`.
- Only append `<until>` when `until != all`.

#### Tool

| Parameter | Type | Default | Description |
|---|---|---|---|
| `--pb-path` | Path | None | Directory containing PhyloBayes tools. Overrides PATH lookup. |
| `--dry-run` | flag | False | Print the bpcomp command without executing. |
| `--quiet` | flag | False | Suppress non-error terminal output. |

### 3.4 Chain Discovery

When `--chain-names all`:
1. List all files in `--chain-dir` matching `*.chain`.
2. Extract stems (filename without extension).
3. Sort lexicographically.
4. Use those stems as chain names.

If no `.chain` files are found, exit code 1 with error: `"No .chain files found in <chain-dir>"`.

### 3.5 bpcomp Invocation

Working directory: `--output-dir`. Chain files referenced as relative paths from output-dir to chain-dir (e.g. `../chains/chain1`).

A single bpcomp call using all resolved chain names:

```
bpcomp -x <burnin> [<every> [<until>]] [-c <cutoff>] -o bpcomp <rel_chain1> <rel_chain2> ...
```

Output basename is always `bpcomp`. Files produced by bpcomp:
- `bpcomp.bpdiff` — summary statistics (maxdiff, meandiff)
- `bpcomp.bplist` — bipartition list
- `bpcomp.con.tre` — majority-rule consensus tree

### 3.6 Output Transparency

`subprocess.run(..., capture_output=False)` — bpcomp stdout/stderr stream directly to the terminal unchanged. PhyloAI does not suppress or buffer tool output. This preserves the full bpcomp screen output the user would see when running bpcomp manually.

### 3.7 Terminal Output

When not `--quiet`:
- Print the bpcomp command before execution.
- bpcomp output streams to terminal directly (all original content visible).
- After bpcomp completes, parse `bpcomp.bpdiff` and print a PhyloAI summary line:
  ```
  PhyloAI: maxdiff 0.043  meandiff 0.003  [good]  → bpcomp/bpcomp.con.tre
  ```

### 3.8 result.json

```json
{
  "status": "success | error",
  "command": "phyloai tree bi bpcomp ...",
  "wall_time": 12.3,
  "tool_versions": { "bpcomp": "1.9", "mpirun": null },
  "params": {
    "chain_dir": "runs/tree/bi/chains",
    "chain_names": "all",
    "output_dir": "runs/tree/bi/bpcomp",
    "overwrite": false,
    "burnin": 1000,
    "sample_freq": 1,
    "until": "all",
    "cutoff": 0.5,
    "pb_path": null,
    "dry_run": false,
    "quiet": false
  },
  "key_results": {
    "chains_used": ["chain1", "chain2", "chain3"],
    "bpcomp_maxdiff": 0.043,
    "bpcomp_meandiff": 0.003,
    "bpcomp_status": "good",
    "consensus_tree": "runs/tree/bi/bpcomp/bpcomp.con.tre"
  },
  "error": null,
  "data": {
    "cmd": ["bpcomp", "-x", "1000", "-c", "0.5", "-o", "bpcomp", "../chains/chain1", "..."],
    "output_files": {
      "bpdiff": { "path": "runs/tree/bi/bpcomp/bpcomp.bpdiff", "description": "bpcomp summary" },
      "bplist": { "path": "runs/tree/bi/bpcomp/bpcomp.bplist", "description": "Bipartition list" },
      "consensus_tree": { "path": "runs/tree/bi/bpcomp/bpcomp.con.tre", "description": "Majority-rule consensus tree" }
    }
  }
}
```

### 3.9 Exit Codes

| Code | Meaning |
|---|---|
| 0 | Success |
| 1 | Input validation error (no chains found, invalid parameters) |
| 2 | bpcomp returned non-zero exit code |
| 3 | `bpcomp` executable not found |

---

## 4. `tree bi tracecomp`

### 4.1 Purpose

Run `tracecomp` once with a user-specified integer burn-in to assess parameter convergence across chains.

### 4.2 Command

```
phyloai tree bi tracecomp [OPTIONS]
```

### 4.3 Parameters

#### Input / Output

| Parameter | Type | Default | Description |
|---|---|---|---|
| `--chain-dir` | Path | **required** | Directory containing chain `.trace` files. |
| `--chain-names` | str | `all` | Comma-separated chain names. `all` = all chains with `.trace` files in `--chain-dir` (sorted lexicographically). |
| `--output-dir` | Path | `runs/tree/bi/tracecomp` | Output directory for `tracecomp.contdiff` and `result.json`. |
| `--overwrite` | flag | False | Delete and recreate the output directory. |

#### Analysis

| Parameter | Type | Default | tc flag | Description |
|---|---|---|---|---|
| `--burnin` | int ≥ 0 | `0` | `-x <burn-in>` | Number of saved samples to discard as burn-in. `0` = no burn-in. |

tracecomp does not support `<every>`, `<until>`, or `-c` — no such parameters are exposed.

#### Tool

| Parameter | Type | Default | Description |
|---|---|---|---|
| `--pb-path` | Path | None | Directory containing PhyloBayes tools. Overrides PATH lookup. |
| `--dry-run` | flag | False | Print the tracecomp command without executing. |
| `--quiet` | flag | False | Suppress non-error terminal output. |

### 4.4 Chain Discovery

Same logic as `bi bpcomp` but matches `*.trace` files instead of `*.chain` files. Error if none found.

### 4.5 tracecomp Invocation

Working directory: `--output-dir`. Trace files referenced as relative paths (e.g. `../chains/chain1.trace`).

```
tracecomp -x <burnin> <rel_chain1.trace> <rel_chain2.trace> ...
```

tracecomp writes its tabular output (the content that becomes `tracecomp.contdiff`) to stdout. PhyloAI captures stdout and:
1. Writes captured stdout to `output-dir/tracecomp.contdiff`.
2. Prints each data line to terminal (same content the user would see running tracecomp manually), with a PhyloAI-appended `[good]`/`[ok]`/`[no]` status column based on per-row `effsize` and `rel_diff` values using the convergence thresholds from `bi pb` (effsize > 300 and rel_diff < 0.1 → `good`; effsize > 50 and rel_diff < 0.3 → `ok`; else `no`).

Implementation: `subprocess.run(stdout=subprocess.PIPE, stderr=None, ...)` — stderr (if any) streams to terminal, stdout is captured. After capture, replay with per-line status annotation.

**Example annotated output:**
```
name                effsize    rel_diff    status
loglik              1529       0.0524634   [good]
length              948        0.0360727   [good]
alpha               1863       0.0694708   [good]
Nmode               2385       0.0170518   [good]
statent             1059       0.0706501   [good]
statalpha           1605       0.0512962   [good]
rrent               2041       0.0318876   [good]
rrmean              2544       0.018511    [good]
```

After printing the annotated table, PhyloAI prints a summary line:
```
PhyloAI: min effsize 948  max rel_diff 0.0706501  [good]
```

### 4.6 result.json

```json
{
  "status": "success | error",
  "command": "phyloai tree bi tracecomp ...",
  "wall_time": 3.1,
  "tool_versions": { "tracecomp": "1.9" },
  "params": {
    "chain_dir": "runs/tree/bi/chains",
    "chain_names": "all",
    "output_dir": "runs/tree/bi/tracecomp",
    "overwrite": false,
    "burnin": 1000,
    "pb_path": null,
    "dry_run": false,
    "quiet": false
  },
  "key_results": {
    "chains_used": ["chain1", "chain2", "chain3"],
    "tracecomp_min_effsize": 312,
    "tracecomp_max_reldiff": 0.094,
    "tracecomp_status": "good"
  },
  "error": null,
  "data": {
    "cmd": ["tracecomp", "-x", "1000", "../chains/chain1.trace", "..."],
    "output_files": {
      "contdiff": { "path": "runs/tree/bi/tracecomp/tracecomp.contdiff", "description": "tracecomp output" }
    }
  }
}
```

### 4.7 Exit Codes

| Code | Meaning |
|---|---|
| 0 | Success |
| 1 | Input validation error |
| 2 | tracecomp returned non-zero exit code |
| 3 | `tracecomp` executable not found |

---

## 5. `tree bi readpb`

### 5.1 Purpose

Run `readpb_mpi` for posterior analysis on a single chain. Supports multiple analysis modes. Automatically converts `rr` output to IQ-TREE exchangeabilities format and `ss` output to IQ-TREE site frequencies format.

### 5.2 Command

```
phyloai tree bi readpb --chain <chain_path> --mode <mode> [OPTIONS]
```

### 5.3 Parameters

#### Input / Output

| Parameter | Type | Default | Description |
|---|---|---|---|
| `--chain` | Path | **required** | Path to chain file without extension (e.g. `runs/tree/bi/chains/chain1`). readpb_mpi initially writes beside the chain; PhyloAI relocates the outputs to `--output-dir`. |
| `--output-dir` | Path | `runs/tree/bi/readpb` | Output directory for `result.json` and all readpb outputs. PhyloAI moves each mode's files here immediately after that mode completes. |
| `--overwrite` | flag | False | Delete and recreate `--output-dir`. |

#### Analysis

| Parameter | Type | Default | Description |
|---|---|---|---|
| `--mode` | str | **required** | Comma-separated list of analysis modes (see table below). Executed in order; each mode is a separate `readpb_mpi` call. |
| `--burnin` | int ≥ 0 | `0` | `-x <burn-in>`. Number of saved samples to discard. |
| `--sample-freq` | int ≥ 1 | `1` | `-x <burn-in> <every>`. Sub-sampling frequency after burn-in. |
| `--until` | str | `all` | `-x <burn-in> <every> <until>`. `all` = to end of chain. Integer = stop at that saved sample index. |

**`-x` flag construction follows the same rule as `bi bpcomp`.**

#### Parallelism

| Parameter | Type | Default | Description |
|---|---|---|---|
| `--threads` | int ≥ 2 | `4` | `mpirun -np <threads>`. MPI processes for readpb_mpi. |

#### Tool

| Parameter | Type | Default | Description |
|---|---|---|---|
| `--pb-path` | Path | None | Directory containing `readpb_mpi` and `mpirun`. Overrides PATH lookup. |
| `--dry-run` | flag | False | Print commands without executing. |
| `--quiet` | flag | False | Suppress non-error terminal output. |

### 5.4 `--mode` Values

| Value | readpb_mpi flag | Output file | Description |
|---|---|---|---|
| `rr` | `-rr` | `<chain>.meanrr` → auto-converted to `<chain>.exchangeabilities` | Posterior mean relative exchangeabilities. Only valid if the model has free exchangeabilities (e.g. GTR). Outputs PAML lower-triangle format for use with IQ-TREE `-m <file>`. |
| `ss` | `-ss` | `<chain>.siteprofiles` → auto-converted to `<chain>.sitefreq` | Posterior mean site-specific state frequencies. Only valid under infinite mixture models (CAT). Outputs IQ-TREE site frequency format for use with `-fs <file>`. |
| `r` | `-r` | `<chain>.meansiterates` | Posterior mean rates across sites. |
| `sitelogl` | `-sitelogl` | `<chain>.sitelogl`, optional `<chain>.cpo` | Site-specific marginal log-likelihoods, plus wAIC/LOO cross-validation quantities. |
| `ppred` | `-ppred` | `<chain>_ppred*.ali` | Simulate data replicates from the posterior predictive distribution. One replicate per saved sample. |
| `div` | `-div` | `<chain>.div.*` | Posterior predictive diversity test (PPA-DIV): mean diversity per site. |
| `sitecomp` | `-sitecomp` | `<chain>.sitecomp.*` | Posterior predictive test of compositional heterogeneity across sites (PPA-VAR). |
| `siteconvprob` | `-siteconvprob` | `<chain>.siteconvprob.*` | Posterior predictive convergence probability test (PPA-CONV): mean squared empirical frequency. |
| `comp` | `-comp` | `<chain>.comp.*` | Posterior predictive test of compositional homogeneity across taxa. |
| `allppred` | `-allppred` | `<chain>.ppred` | Combined posterior predictive checks. Mutually exclusive with `div`, `sitecomp`, `siteconvprob`, `comp`. |

### 5.5 Mutual Exclusion and Validation

- `allppred` combined with any of `div`, `sitecomp`, `siteconvprob`, `comp` → exit code 1: `"--mode allppred is mutually exclusive with div, sitecomp, siteconvprob, comp"`.
- Any unrecognised mode value → exit code 1.
- Duplicate mode values (e.g. `--mode rr,rr`) → exit code 1.

### 5.6 Execution

Modes are executed in the order specified by the user. Each mode is a separate `readpb_mpi` invocation:

```
mpirun -np <threads> readpb_mpi -x <burnin> [<every> [<until>]] <mode_flag(s)> <chain_stem>
```

Working directory: parent directory of `--chain` (required because readpb_mpi derives output names from the chain stem). Immediately after each successful mode, PhyloAI moves that mode's files to `--output-dir` before starting the next mode. This keeps the chain directory free of readpb analysis outputs.

Example for `--mode ss,rr --burnin 1000 --threads 8`:
1. `mpirun -np 8 readpb_mpi -x 1000 -ss chain1` (cwd = `chains/`)
2. `mpirun -np 8 readpb_mpi -x 1000 -rr chain1` (cwd = `chains/`)

Output transparency: use `stdout=subprocess.PIPE, stderr=None` so stderr (progress) streams to the terminal. If stdout is non-empty, write it to `<output-dir>/<mode>.stdout` and, unless `--quiet` is set, immediately replay it unchanged with `sys.stdout.write(proc.stdout)`. This preserves the screen output of a manual `readpb_mpi` invocation while retaining an audit copy per mode.

### 5.7 Post-Processing

#### `rr` → `exchangeabilities`

After `readpb_mpi -rr` completes, PhyloAI converts `<chain>.meanrr` to `<chain>.exchangeabilities` using numpy (no pandas).

**`.meanrr` file format:**
```
A C D E F G H I K L M N P Q R S T V W Y
                                          (blank line)
A   C   1.37435
A   D   0.486164
...
```

**Algorithm (pure numpy, no pandas):**
1. Read first non-empty line → `order_of_aa` (list of 20 AA symbols, space-separated).
2. Build index map: `aa_to_idx = {aa: i for i, aa in enumerate(order_of_aa)}`.
3. Initialise `exch = np.zeros((20, 20), dtype=np.float64)`.
4. For each remaining non-empty line: parse `source target value`; set `exch[aa_to_idx[source], aa_to_idx[target]] = exch[aa_to_idx[target], aa_to_idx[source]] = float(value)`.
5. Reorder to PAML order: `PAML_ORDER = ['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']`. Build `paml_idx = [aa_to_idx[aa] for aa in PAML_ORDER]`. Reindex: `paml_exch = exch[np.ix_(paml_idx, paml_idx)]`.
6. Write lower triangle for rows `i=0..19`, columns `0..i-1`, one row per line, values formatted as `%08.6f` space-separated, trailing space, then newline. Row `i=0` has no lower-triangle values and is therefore the required leading blank line. (The diagonal is not included; PAML format is the strict lower triangle only.)
7. Append a blank line, then `0.050000 ` × 20 + newline (uniform prior state frequencies placeholder).

Note: the reference script `convert-exchangeabilities.py` produces correct output but uses pandas with deprecated chained-assignment (hence the `ChainedAssignmentError` warnings). The numpy reimplementation avoids pandas entirely and produces identical output without warnings.

The raw `<chain>.meanrr` file is moved to `--output-dir` immediately after `rr` completes. Conversion then writes `<chain>.exchangeabilities` beside it in `--output-dir`.

#### `ss` → `sitefreq`

After `readpb_mpi -ss` completes, PhyloAI converts `<chain>.siteprofiles` to `<chain>.sitefreq` using numpy.

**`.siteprofiles` file format:**
```
<header line 1>
<header line 2>
<site_index> <freq_A> <freq_C> <freq_D> ... (PhyloBayes AA order)
...
```

PhyloBayes AA order: `A C D E F G H I K L M N P Q R S T V W Y`
IQ-TREE AA order: `A R N D C Q E G H I L K M F P S T W Y V`

**Algorithm:**
1. Skip first two lines.
2. For each data line: parse site index + 20 float frequencies.
3. Reindex from PhyloBayes order to IQ-TREE order using a precomputed index map.
4. Replace zeros/near-zeros with `1e-8` (floor), then re-normalize to sum to 1.
5. Write to `<chain>.sitefreq`: `<site_index> <20 space-separated floats %.8f>`.

The raw `<chain>.siteprofiles` file is moved to `--output-dir` immediately after `ss` completes. Conversion then writes `<chain>.sitefreq` beside it in `--output-dir`.

### 5.8 Output Structure

`readpb_mpi` initially writes beside the chain, but PhyloAI moves files immediately after each mode completes. All final analysis outputs therefore reside directly under `--output-dir`. The chain directory retains only the input chain and MCMC files.

```
runs/tree/bi/
├── chains/
│   ├── chain1.chain
│   └── chain1.trace
│
└── readpb/
    ├── chain1.meanrr
    ├── chain1.exchangeabilities
    ├── chain1.siteprofiles
    ├── chain1.sitefreq
    ├── chain1.meansiterates
    ├── chain1.sitelogl
    ├── chain1.sitecomp
    ├── chain1.ppred
    └── result.json
```

### 5.9 result.json

```json
{
  "status": "success | error",
  "command": "phyloai tree bi readpb --chain runs/tree/bi/chains/chain1 --mode ss,rr ...",
  "wall_time": 42.1,
  "tool_versions": { "readpb_mpi": "1.9", "mpirun": "4.1.2" },
  "params": {
    "chain": "runs/tree/bi/chains/chain1",
    "output_dir": "runs/tree/bi/readpb",
    "overwrite": false,
    "mode": "ss,rr",
    "burnin": 1000,
    "sample_freq": 1,
    "until": "all",
    "threads": 4,
    "pb_path": null,
    "dry_run": false,
    "quiet": false
  },
  "key_results": {
    "modes_run": ["ss", "rr"],
    "output_files": {
      "siteprofiles": "runs/tree/bi/chains/chain1.siteprofiles",
      "sitefreq": "runs/tree/bi/chains/chain1.sitefreq",
      "meanrr": "runs/tree/bi/chains/chain1.meanrr",
      "exchangeabilities": "runs/tree/bi/chains/chain1.exchangeabilities"
    }
  },
  "error": null,
  "data": {
    "cmds": {
      "ss": ["mpirun", "-np", "4", "readpb_mpi", "-x", "1000", "-ss", "chain1"],
      "rr": ["mpirun", "-np", "4", "readpb_mpi", "-x", "1000", "-rr", "chain1"]
    },
    "post_processing": {
      "ss": { "input": "chain1.siteprofiles", "output": "chain1.sitefreq", "status": "success" },
      "rr": { "input": "chain1.meanrr", "output": "chain1.exchangeabilities", "status": "success" }
    }
  }
}
```

### 5.10 Exit Codes

| Code | Meaning |
|---|---|
| 0 | Success (all modes completed) |
| 1 | Input validation error (bad mode, allppred conflict, chain not found) |
| 2 | readpb_mpi returned non-zero exit code |
| 3 | `readpb_mpi` or `mpirun` not found |

---

## 6. `tree bi pb` (Renamed from `tree bi`)

No functional changes. Identical parameters, behaviour, and output as the original `tree bi` command. The only change:

- CLI registration: `@bi.command("pb", ...)` instead of `@tree.command("bi", ...)`.
- Library: `run_bi_pb()` in `bi.py` (renamed from `run_bi()`).
- OPTION_GROUPS key: `"phyloai tree bi pb"` (was `"phyloai tree bi"`).
- MCP tool name: `phyloai_tree_bi_pb` (was `phyloai_tree_bi`).
- Help text: update example commands from `phyloai tree bi ...` to `phyloai tree bi pb ...`.

The help text for the `bi` group itself briefly describes all four subcommands:

```
Usage: phyloai tree bi COMMAND [ARGS]...

  Bayesian phylogenetic inference with PhyloBayes-MPI.

  Subcommands:
    pb         Run MCMC chains (pb_mpi).
    bpcomp     Topology convergence analysis (bpcomp).
    tracecomp  Parameter convergence analysis (tracecomp).
    readpb     Posterior analysis and predictive checks (readpb_mpi).
```

---

## 7. Shared Helper Reuse

All shared parsing and status helpers remain in `bi.py`. The new modules import only what they need:

| Function in `bi.py` | Used by |
|---|---|
| `_parse_bpcomp_bpdiff()` | `bi_bpcomp` |
| `_parse_tracecomp_contdiff()` | `bi_tracecomp` |
| `_bpcomp_status()` | `bi_bpcomp` |
| `_tracecomp_status()` | `bi_tracecomp` |

`_resolve_chain_names()` from `bi.py` handles explicit names; chain auto-discovery from directory (`.chain` / `.trace` glob) is new logic in `bi_bpcomp.py` / `bi_tracecomp.py` respectively.

Each new module resolves its own tools with `ToolEnv` and calls `env.require()` only for the executable(s) it needs. It must not call `_detect_tools()`, because that helper intentionally requires the complete live-monitoring tool set (`pb_mpi`, `bpcomp`, `tracecomp`, and `mpirun`). Each new module also builds its own `tool_versions` dictionary from its resolved executable paths.

---

## 8. Tool Detection Requirements

`bi bpcomp` requires: `bpcomp` (exit 3 if missing).
`bi tracecomp` requires: `tracecomp` (exit 3 if missing).
`bi readpb` requires: `readpb_mpi`, `mpirun` (exit 3 if either missing).

`--pb-path` propagation in all three subcommands follows the same pattern as `bi pb`:

```python
if pb_path:
    tool_paths = {tool: pb_path / tool for tool in required_tools}
    # readpb_mpi: also check pb_path / "readpb_mpi"
else:
    tool_paths = {}
env = ToolEnv(tool_paths=tool_paths)
```

---

## 9. stdout Passthrough Strategy

### bpcomp

bpcomp writes output to both files and to stdout/stderr. Use `capture_output=False` — all output streams directly to terminal unchanged. No post-processing of bpcomp's own output; PhyloAI only adds a summary line after completion (parsed from `bpcomp.bpdiff`).

### tracecomp

tracecomp writes its tabular result to stdout. Capture stdout, then replay with per-line status annotations:

```python
proc = subprocess.run(
    cmd,
    stdout=subprocess.PIPE,
    stderr=None,   # stderr streams to terminal directly
    cwd=working_dir,
    text=True,
)
# save
(output_dir / "tracecomp.contdiff").write_text(proc.stdout)
# annotate and print
_print_annotated_tracecomp(proc.stdout)
```

### readpb_mpi

readpb_mpi writes progress to stderr and results to stdout. Capture stdout per mode, save it to `<output-dir>/<mode>.stdout` if non-empty, then replay it unchanged with `sys.stdout.write(proc.stdout)` unless `--quiet` is set. Stderr streams to the terminal directly.

---

## 10. Output Directory Structure (Full Picture)

```
runs/tree/bi/
├── chains/                    ← written by bi pb
│   ├── chain1.{trace,treelist,chain,param,monitor,run,log}
│   ├── chain1.meanrr          ← written by bi readpb --mode rr
│   ├── chain1.exchangeabilities ← written by PhyloAI post-processing
│   ├── chain1.siteprofiles    ← written by bi readpb --mode ss
│   ├── chain1.sitefreq        ← written by PhyloAI post-processing
│   └── ...
│
├── convergence/               ← written by bi pb (live monitoring)
│   ├── trace_plots.pdf
│   ├── bpcomp_all.{bpdiff,bplist,con.tre}
│   └── ...
│
├── bpcomp/                    ← written by bi bpcomp
│   ├── bpcomp.{bpdiff,bplist,con.tre}
│   └── result.json
│
├── tracecomp/                 ← written by bi tracecomp
│   ├── tracecomp.contdiff
│   └── result.json
│
├── readpb/                    ← written by bi readpb (result.json only)
│   └── result.json
│
├── run_state.json             ← written by bi pb
└── result.json                ← written by bi pb
```

---

## 11. Documentation Updates Required

The following documents must be updated as part of this implementation:

| Document | Change |
|---|---|
| `docs/commands/tree-bi.md` | Replace with new structure: `bi pb` (original content), `bi bpcomp`, `bi tracecomp`, `bi readpb` sections; update all `phyloai tree bi` examples to `phyloai tree bi pb` |
| `docs/commands/tree-bi.zh.md` | Same as above in Chinese |
| `docs/superpowers/specs/2026-06-23-phyloai-tree-bi-design.md` | Add note at top referencing this spec as superseding the `bi` command structure |
| `docs/superpowers/specs/2026-06-07-phyloai-design.md` | Update CLI table: `tree bi` → `tree bi pb` |
| `docs/superpowers/specs/2026-06-17-phyloai-tree-design.md` | Update `bi` entry to reflect subcommand group |
| `docs/superpowers/plans/2026-07-26-phyloai-tree-bi-subcommands.md` | Add the implementation plan for this approved spec |
| `phyloai/mcp/schema_gen.py` | No manual tool registry change: verify its Click-tree walk automatically replaces `phyloai_tree_bi` with `phyloai_tree_bi_pb` and adds `phyloai_tree_bi_bpcomp`, `phyloai_tree_bi_tracecomp`, and `phyloai_tree_bi_readpb` |
| `skills/phyloai-workflow/SKILL.md` | Update `tree bi` workflow and required-tool guidance for the four subcommands |
| `skills/phyloai-workflow/references/parameter-annotations.md` | Replace `### tree bi` with sections for `tree bi pb`, `tree bi bpcomp`, `tree bi tracecomp`, and `tree bi readpb` |
| `skills/phyloai-workflow/references/workflow.md` | Update Bayesian workflow command from `tree bi` to `tree bi pb`; add final diagnostic and readpb workflow guidance |
| `skills/phyloai-workflow/references/error-catalog.md` | Split PhyloBayes missing-tool guidance by subcommand requirements |
| `pyproject.toml` | Bump package version from `0.3.0` to `0.4.0` |
| `phyloai/__init__.py` | Bump `__version__` from `0.3.0` to `0.4.0` |
| README (if present) | Update `tree bi` usage example |

---

## 12. `phyloai report` Template for `tree bi` Subcommands

Each subcommand must have a report template registered in the `phyloai report` system. Templates must be detailed, covering methods, parameters, and results.

### `bi bpcomp` Report Template

**Methods paragraph:**
> Topology convergence was assessed using `bpcomp` (PhyloBayes-MPI v{version}) applied to {n_chains} independent MCMC chains ({chain_names}). A burn-in of {burnin} saved samples was discarded; trees were sub-sampled every {sample_freq} points{until_clause}. The majority-rule consensus cutoff was set to {cutoff}. The maximum bipartition frequency discrepancy (maxdiff) between chains was {maxdiff:.4f} and the mean discrepancy (meandiff) was {meandiff:.6f}, indicating {status} convergence. The consensus tree was written to `{consensus_tree}`.

**Key results table:** chains used, burnin, maxdiff, meandiff, status, consensus tree path.

### `bi tracecomp` Report Template

**Methods paragraph:**
> Continuous parameter convergence was assessed using `tracecomp` (PhyloBayes-MPI v{version}) applied to {n_chains} chains with a burn-in of {burnin} saved samples. The minimum effective sample size across all parameters was {min_effsize:.0f} and the maximum relative difference was {max_reldiff:.4f}, indicating {status} mixing. Per-parameter diagnostics are summarised in the table below.

**Key results table:** parameter name, effsize, rel_diff, status — one row per parameter from `tracecomp.contdiff`.

### `bi readpb` Report Template

**Methods paragraph (one paragraph per mode run):**

- `rr`: Posterior mean relative exchangeabilities were estimated from {n_samples} post-burnin samples of chain `{chain}` using `readpb_mpi -rr`. The resulting exchangeability matrix was converted to PAML lower-triangle format (`{chain}.exchangeabilities`) for use with IQ-TREE.
- `ss`: Posterior mean site-specific amino acid frequencies were estimated from {n_samples} post-burnin samples using `readpb_mpi -ss`. Site frequency profiles were converted to IQ-TREE `-fs` format (`{chain}.sitefreq`).
- `r`: Posterior mean site rates were estimated using `readpb_mpi -r` and written to `{chain}.meansiterates`.
- `sitelogl`: Site-specific marginal log-likelihoods were computed using `readpb_mpi -sitelogl` and written to `{chain}.sitelogl`. These values can be used to compute wAIC and leave-one-out cross-validation scores.
- `div`/`sitecomp`/`siteconvprob`/`comp`/`allppred`: Posterior predictive checks were performed using `readpb_mpi -{mode}`. Observed test statistics were compared to the posterior predictive null distribution.

**Key results table:** modes run, output files, post-processing status (for `rr` and `ss`).

---

## 13. Relationship to Parent Specs

- **`2026-06-23-phyloai-tree-bi-design.md`**: This spec supersedes the command structure section (Section 1, "no subcommand layer"). All other sections of that spec remain valid for `tree bi pb`.
- **`2026-06-07-phyloai-design.md`**: `tree bi` now refers to a group, not a leaf command. Output directory convention `runs/tree/bi/` remains unchanged.
- **`2026-06-21-phyloai-json-output-standard.md`**: `bi bpcomp`, `bi tracecomp`, `bi readpb` all conform to the single-execution pattern (one `cmd` key, or `cmds` dict for multi-mode readpb).
- **`2026-06-18-phyloai-doctor-design.md`**: `readpb_mpi` is already registered in `TOOL_REGISTRY` as optional. No registry changes required; `bi readpb` uses `env.require("readpb_mpi")` which raises `FileNotFoundError` (exit 3) if absent.
