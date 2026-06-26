# PhyloAI `run` Command Design Specification

**Date:** 2026-06-26  
**Status:** Approved for implementation  
**Parent spec:** `2026-06-07-phyloai-design.md`  
**Related specs:**
- `2026-06-12-checkpoint-resume-design.md` — checkpoint/resume infrastructure
- `2026-06-11-pretree-align-design.md`
- `2026-06-12-pretree-trim-design.md`
- `2026-06-15-phyloai-pretree-filter-design.md`
- `2026-06-13-phyloai-pretree-concat-design.md`
- `2026-06-18-phyloai-tree-ml-fasttree-design.md`
- `2026-06-19-phyloai-tree-ml-iqtree-design.md`
- `2026-06-20-phyloai-tree-msc-design.md`

---

## 1. Purpose

`phyloai run` is the one-click pipeline entry point. It orchestrates the full phylogenomics workflow from raw sequence files to a final species tree, using sensible defaults for all intermediate steps. Users who want fine-grained control over individual steps should use the constituent subcommands directly.

`phyloai run` calls the Python library layer of each step directly (no subprocess invocation of the CLI), consistent with the architecture defined in the main design: Library → CLI → MCP.

---

## 2. CLI Interface

```bash
# Supermatrix mode (default): concatenated matrix → ML species tree
phyloai run --seq-dir ./markers --mode supermatrix [--speed normal|fast]
            [--output-dir ./runs/run] [--threads 8] [--resume] [--dry-run]
            [--overwrite] [--quiet]

# Supertree mode: per-gene trees → coalescent species tree
phyloai run --seq-dir ./markers --mode supertree [--speed normal|fast]
            [--output-dir ./runs/run] [--threads 8] [--resume] [--dry-run]
            [--overwrite] [--quiet]
```

### 2.1 Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--seq-dir` | Path | required | Input sequence directory. Any format accepted; sequences are always run through `pretree convert` first to normalise encoding and line wrapping. |
| `--mode` | `supermatrix\|supertree` | `supermatrix` | Pipeline type. `supermatrix`: concatenated matrix → IQ-TREE (or FastTree in fast mode). `supertree`: per-gene trees → wASTRAL. |
| `--speed` | `normal\|fast` | `normal` | Speed/accuracy trade-off. Controls tool selection and whether the TAPER filter step is included. See Section 3. |
| `--output-dir` / `-o` | Path | `./runs/run` | Root output directory for all pipeline steps. |
| `--threads` / `-t` | int | 4 | Thread count passed to all steps that support parallelism. |
| `--resume` | flag | False | Resume from `run_checkpoint.json`. Skips completed steps; passes `--resume` to any interrupted step's subcommand. |
| `--overwrite` | flag | False | Delete and recreate the output directory, then run from scratch. |
| `--dry-run` | flag | False | Print what would be executed at each step without running any tools. |
| `--quiet` / `-q` | flag | False | Suppress all terminal output except errors. |

`--resume` and `--overwrite` are mutually exclusive. Supplying both exits with code 1.

### 2.2 `--help` text

The `--help` output must include:

1. A one-paragraph description of what the command does.
2. A table or list of the two modes with their step sequences spelled out explicitly.
3. A `--speed` explanation showing which tools are used in each mode.
4. An example for each of the two pipeline modes.
5. A note pointing to constituent subcommands for fine-grained control.

Example help excerpt:

```
Usage: phyloai run [OPTIONS]

  One-click phylogenomics pipeline from raw sequences to a species tree.

  Runs the full pipeline automatically using sensible defaults.
  For fine-grained control over individual steps, use the constituent
  subcommands (phyloai pretree align, phyloai tree ml iqtree, etc.).

  Modes:
    supermatrix  convert → align → trim → [filter] → concat → iqtree
    supertree    convert → align → trim → [filter] → gene trees → wastral

  The [filter] step (TAPER error-site masking) is included in --speed normal
  and skipped in --speed fast.

  Speed modes:
    normal  MAFFT linsi, trimAl -automated1, TAPER filter, IQ-TREE3 / FastTree
    fast    MAFFT auto, trimAl -automated1, no filter, FastTree

Options:
  --seq-dir PATH          Input sequence directory (any format)  [required]
  --mode [supermatrix|supertree]
                          Pipeline mode  [default: supermatrix]
  --speed [normal|fast]   Speed/accuracy trade-off  [default: normal]
  -o, --output-dir PATH   Output directory  [default: runs/run]
  -t, --threads INTEGER   Thread count  [default: 4]
  --resume                Resume from checkpoint
  --overwrite             Overwrite existing output directory
  --dry-run               Show steps without running
  -q, --quiet             Suppress non-error output
  --help                  Show this message and exit.

Examples:
  phyloai run --seq-dir ./markers --mode supermatrix
  phyloai run --seq-dir ./markers --mode supertree --speed fast --threads 16
  phyloai run --seq-dir ./markers --mode supermatrix --resume
```

---

## 3. Pipeline Steps and Tool Mapping

### 3.1 Supermatrix mode

| Step | Subdirectory | `--speed normal` | `--speed fast` |
|------|--------------|-----------------|----------------|
| 1. Convert | `1-convert/` | `pretree convert` | same |
| 2. Align | `2-align/` | `pretree align --method linsi` | `pretree align --method auto` |
| 3. Trim | `3-trim/` | `pretree trim --tool trimal --trimal-mode automated1` | same |
| 4. Filter | `4-filter/` | `pretree filter taper` | **skipped** (directory not created) |
| 5. Concat | `5-concat/` | `pretree concat` | same |
| 6. Tree | `6-tree/` | `tree ml iqtree --matrix ... (unpartitioned, no partition model)` | `tree ml fasttree --matrix ...` |

### 3.2 Supertree mode

| Step | Subdirectory | `--speed normal` | `--speed fast` |
|------|--------------|-----------------|----------------|
| 1. Convert | `1-convert/` | `pretree convert` | same |
| 2. Align | `2-align/` | `pretree align --method linsi` | `pretree align --method auto` |
| 3. Trim | `3-trim/` | `pretree trim --tool trimal --trimal-mode automated1` | same |
| 4. Filter | `4-filter/` | `pretree filter taper` | **skipped** |
| 5. Gene Trees | `5-genetrees/` | `tree ml fasttree --msa-dir ... --mode normal` | `tree ml fasttree --msa-dir ... --mode fast` |
| 6. Species Tree | `6-tree/` | `tree msc --mode 1` | same |

### 3.3 Step count and `--speed fast`

In `--speed fast`, the filter step is omitted. The step counter displayed to the user reflects the actual steps to be run, not the maximum possible steps. A fast-mode run shows `[1/5]` through `[5/5]`, not `[1/6]`.

### 3.4 Default tool parameters

`phyloai run` uses each step's own default parameters unless noted above. It does not expose per-step tool parameters — users needing non-default tool parameters should run steps individually.

The `--threads` value from `phyloai run` is passed through to all steps that accept `--threads`.

### 3.5 IQ-TREE3 unpartitioned mode

In `--speed normal` supermatrix mode, IQ-TREE3 is called without a partition file and without a user-specified substitution model. IQ-TREE3 performs automatic model selection (ModelFinder) and builds the ML tree in a single unpartitioned analysis. This is intentional: `phyloai run` is a quick pipeline for a first-pass result. Partitioned or mixture-model analyses should be run via `phyloai tree ml iqtree` directly.

---

## 4. Progress Display

Each step prints a header line when it starts and a completion line when it finishes. Batch steps (align, trim, filter, gene trees) show a Rich progress bar. Single steps (convert, concat, tree) show a spinner.

```
[1/6] Converting sequences ...
      ✓ 100 sequences converted  [0:00:08]

[2/6] Aligning sequences (MAFFT linsi) ...
      ████████████░░░░░░  68/100 genes  ETA 3:12
      ✓ 100 alignments  [0:28:41]

[3/6] Trimming alignments (trimAl -automated1) ...
      ████████████████████  100/100 genes
      ✓ 100 trimmed alignments  [0:01:22]

[4/6] Filtering error sites (TAPER) ...
      ████████████████████  100/100 genes
      ✓ 98 genes passed  [0:04:15]

[5/6] Concatenating matrix ...
      ✓ matrix: 52 taxa × 43820 sites  [0:00:03]

[6/6] Inferring species tree (IQ-TREE3, unpartitioned) ...
      ✓ Tree complete  [1:12:04]

✓ Pipeline complete  [total: 2h 14m]
  Species tree:  runs/run/6-tree/iqtree.treefile
  Results:       runs/run/result.json
```

For fast mode without the filter step, the step counter shows `[1/5]` through `[5/5]`.

`--quiet` suppresses all of the above except the final summary line with file paths, and errors.

---

## 5. Output Directory Structure

```
runs/run/
├── run_checkpoint.json         # run-level checkpoint (step status only)
├── result.json                 # final run result (written on completion)
│
├── 1-convert/
│   ├── result.json
│   └── seqs/
│
├── 2-align/
│   ├── result.json
│   ├── checkpoint.json         # align subcommand's own per-gene checkpoint
│   ├── seqs/
│   └── logs/
│
├── 3-trim/
│   ├── result.json
│   ├── checkpoint.json
│   ├── seqs/
│   └── logs/
│
├── 4-filter/                   # only present in --speed normal
│   ├── result.json
│   ├── seqs/
│   └── logs/
│
│   [supermatrix only]
├── 5-concat/
│   ├── result.json
│   ├── matrix.fa
│   └── matrix.partitions
│
│   [supertree only]
├── 5-genetrees/
│   ├── result.json
│   ├── checkpoint.json
│   ├── trees/
│   └── logs/
│
└── 6-tree/
    ├── result.json
    └── <tool-native output files>
```

Numeric prefixes ensure natural sort order in file managers and `ls`. The prefix matches the step number displayed during progress output.

`4-filter/` is not created at all when `--speed fast` is used. The downstream step reads from `3-trim/seqs/` directly.

---

## 6. `--resume` Mechanism

### 6.1 `run_checkpoint.json` schema

`run_checkpoint.json` is written at `<output-dir>/run_checkpoint.json`. It tracks step-level status only. Fine-grained per-gene progress within a step is managed by each subcommand's own `checkpoint.json`.

```json
{
  "schema_version": 1,
  "step": "run",
  "command": "phyloai run --seq-dir ./markers --mode supermatrix --speed normal --threads 8",
  "status": "running",
  "params_hash": "sha256:...",
  "params": {
    "seq_dir": "./markers",
    "mode": "supermatrix",
    "speed": "normal",
    "threads": 8,
    "output_dir": "runs/run"
  },
  "started_at": "2026-06-26T10:00:00",
  "updated_at": "2026-06-26T12:30:00",
  "completed_at": null,
  "steps": [
    {"name": "convert",      "status": "success",     "output_dir": "runs/run/1-convert"},
    {"name": "align",        "status": "success",     "output_dir": "runs/run/2-align"},
    {"name": "trim",         "status": "interrupted", "output_dir": "runs/run/3-trim"},
    {"name": "filter_taper", "status": "pending",     "output_dir": "runs/run/4-filter"},
    {"name": "concat",       "status": "pending",     "output_dir": "runs/run/5-concat"},
    {"name": "tree",         "status": "pending",     "output_dir": "runs/run/6-tree"}
  ]
}
```

### 6.2 Step status values

| Status | Resume behaviour |
|--------|-----------------|
| `pending` | Run step normally |
| `running` | Treat as interrupted; pass `--resume` to subcommand if it supports it; otherwise rerun from scratch |
| `success` | Verify `<output_dir>/result.json` exists and has `status: "success"`. Skip if valid; rerun if invalid. |
| `failed` | Rerun step |
| `interrupted` | Pass `--resume` to subcommand if it supports it; otherwise rerun from scratch |
| `skipped` | Keep skipped (used for filter step in `--speed fast`) |

### 6.3 Resume flow

1. Load `run_checkpoint.json`. Exit 1 if missing, malformed, or `schema_version` is unsupported.
2. Validate `params_hash` against current invocation. Exit 1 if mismatch, with a diff of changed parameters.
3. Walk the step list in order. For each step:
   - `success` → verify output, skip if valid
   - `interrupted` or `running` → call subcommand with `--resume` if it supports it, otherwise `--overwrite` that step's directory
   - `pending` or `failed` → run step normally
   - `skipped` → skip
4. On completion, update `run_checkpoint.json` top-level `status` to `"success"` and write `result.json`.

### 6.4 Which steps support `--resume`

| Step | Subcommand resume support |
|------|--------------------------|
| convert | No (fast utility; reruns in seconds) |
| align | Yes (`checkpoint.json` per gene) |
| trim | Yes (`checkpoint.json` per gene) |
| filter taper | No (fast; rerun) |
| concat | No (single-shot; rerun) |
| gene trees (fasttree batch) | Yes (`checkpoint.json` per gene) |
| iqtree | Yes (IQ-TREE3 native checkpoint) |
| wastral | No (one-shot; rerun) |

Steps without resume support: if their status is `interrupted` or `running`, `phyloai run --resume` reruns that step from scratch using `--overwrite` on that step's subdirectory only.

### 6.5 Atomic writes and `run_checkpoint.json`

`run_checkpoint.json` is updated using `save_checkpoint_atomic` from `phyloai/core/checkpoint.py`, the same helper used by subcommands. The run-level checkpoint is written:
- once after each step completes (step count is small; no throttle needed at this level)
- on `KeyboardInterrupt` (status set to `interrupted`, `fsync=True`)
- on final completion (status set to `success`, `fsync=True`)

---

## 7. `result.json` Schema

`result.json` is written at `<output-dir>/result.json` only on successful pipeline completion.

```json
{
  "status": "success",
  "command": "phyloai run --seq-dir ./markers --mode supermatrix --speed normal --threads 8",
  "wall_time": 8040.5,
  "tool_versions": {
    "mafft": "7.520",
    "trimal": "1.4.1",
    "iqtree": "3.0.0"
  },
  "params": {
    "seq_dir": "./markers",
    "mode": "supermatrix",
    "speed": "normal",
    "threads": 8,
    "output_dir": "runs/run"
  },
  "key_results": {
    "n_input_genes": 200,
    "n_genes_after_filter": 187,
    "matrix_length": 45230,
    "matrix_taxa": 52,
    "final_tree": "runs/run/6-tree/iqtree.treefile"
  },
  "error": null,
  "data": {
    "mode": "supermatrix",
    "speed": "normal",
    "steps": [
      {"name": "convert",      "status": "success", "wall_time": 12.3,   "result_json": "runs/run/1-convert/result.json"},
      {"name": "align",        "status": "success", "wall_time": 3200.1, "result_json": "runs/run/2-align/result.json"},
      {"name": "trim",         "status": "success", "wall_time": 82.4,   "result_json": "runs/run/3-trim/result.json"},
      {"name": "filter_taper", "status": "success", "wall_time": 255.0,  "result_json": "runs/run/4-filter/result.json"},
      {"name": "concat",       "status": "success", "wall_time": 3.1,    "result_json": "runs/run/5-concat/result.json"},
      {"name": "tree",         "status": "success", "wall_time": 4487.6, "result_json": "runs/run/6-tree/result.json"}
    ]
  }
}
```

`key_results` fields:

| Field | Description |
|-------|-------------|
| `n_input_genes` | Number of valid input sequence files after convert |
| `n_genes_after_filter` | Number of genes after the filter step (equals `n_input_genes` in fast mode) |
| `matrix_length` | Concatenated matrix length in sites (supermatrix only; omitted for supertree) |
| `matrix_taxa` | Number of taxa in the final matrix or gene tree set |
| `final_tree` | Path to the final species tree file |

If the pipeline fails mid-run, `result.json` is written with `status: "error"` and `error` populated. Steps not yet started have `status: "pending"` in `data.steps`.

---

## 8. Error Handling

| Situation | Behaviour |
|-----------|-----------|
| `--seq-dir` does not exist | Exit 1 before starting |
| `--resume` with no `run_checkpoint.json` | Exit 1: tell user to use `--overwrite` or check path |
| `--resume` + `--overwrite` | Exit 1 |
| params_hash mismatch on resume | Exit 1: show which params changed |
| A step fails (subcommand exits non-zero) | Write `run_checkpoint.json` with step `status: "failed"`, write `result.json` with `status: "error"`, exit 2 |
| All genes filtered out before concat | Exit 1 with clear message from the filter step |
| External tool not found | Exit 3 (propagated from subcommand) |

When a step fails, the error message displayed to the user includes:
- which step failed
- the step's output directory for inspection
- the path to the step's `result.json` for details
- the suggestion to use `--resume` after fixing the issue

---

## 9. Internal Architecture

`phyloai run` is implemented in `phyloai/cli/commands/run.py` and calls library functions directly. It does not invoke subcommands via subprocess.

Each step is dispatched by calling the corresponding module's entry function, e.g.:

```python
from phyloai.pretree.convert import run_convert
from phyloai.pretree.align import run_align
from phyloai.pretree.trim import run_trim
from phyloai.pretree.filter import run_filter_taper
from phyloai.pretree.concat import run_concat
from phyloai.tree.ml import run_fasttree, run_iqtree
from phyloai.tree.msc import run_wastral
```

Each module function accepts a parameter dataclass or kwargs that mirror the CLI parameters, and returns a result dict (or raises on failure). The `run` orchestrator builds the parameter objects for each step based on the run-level `--mode`, `--speed`, `--threads`, and resolved output subdirectory paths.

The `command` field in each step's `result.json` is the equivalent CLI command string (constructed by the library function, not from an actual subprocess), consistent with the JSON Output Standard.

---

## 10. Main Design Updates Required

The following updates to `2026-06-07-phyloai-design.md` are needed after this spec is accepted:

1. **Section 4.1 command examples:** Update `--mode coalescent` to `--mode supertree` in the `phyloai run` examples.
2. **Section 4.2 pipeline table:** Replace the coalescent row with supertree; update both rows to match the step sequences in this spec (adding `convert` as first step, specifying `--speed` variants).
3. **Section 4.3 universal flags:** Add `--speed normal|fast` to the shared parameter registry table (applicable to `run` only).
4. **Section 9.2 shared parameter registry:** Add `--speed` row.

---

## 11. Documentation Requirements

Implementation must add or update:

- `docs/commands/run.md` — full command reference following the standard template (Purpose, Usage, Inputs, Outputs, Examples, Warnings/Errors, Notes)
- `README.md` — add `phyloai run` to the commands table

---

## 12. Testing Requirements

- `--help` output includes step sequences for both modes
- default `--output-dir` is `./runs/run`
- `--resume` + `--overwrite` exits 1
- `--resume` without `run_checkpoint.json` exits 1
- params_hash mismatch on resume exits 1
- dry-run prints step list and estimated tool invocations without executing
- fast mode skips filter step and shows `[1/5]` counter
- resume skips completed steps and continues from interrupted step
- failed step writes `result.json` with `status: "error"` and correct step name in error message
- on success, `result.json` `key_results.final_tree` points to an existing file
- `4-filter/` directory is not created in fast mode
- supertree mode creates `5-genetrees/` not `5-concat/`
