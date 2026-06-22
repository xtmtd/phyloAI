# PhyloAI Tree ML FastTree Design Specification

**Date:** 2026-06-18
**Status:** Approved
**Parent spec:** `2026-06-07-phyloai-design.md`, `2026-06-17-phyloai-tree-design.md`

---

## 1. Purpose

`phyloai tree ml fasttree` performs maximum-likelihood tree inference using FastTree. It supports two input modes:

- **Batch gene tree mode** (`--msa-dir`): parallel gene tree inference from a directory of MSA files using `ProcessPoolExecutor`
- **Single supermatrix mode** (`--matrix`): single-tree inference from one concatenated matrix file

---

## 2. CLI Surface

```bash
# Batch gene trees: parallel inference from MSA directory
phyloai tree ml fasttree --msa-dir ./trimmed/seqs \
    --seq-type AA --model lg --mode normal --boot 1000 \
    --cat 20 --gamma --threads 8 -o runs/tree/ml/fasttree

# Single supermatrix tree
phyloai tree ml fasttree --matrix ./concat/matrix.fa \
    --seq-type NT --model gtr --mode slow --boot 1000 \
    -o runs/tree/ml/fasttree

# Disable bootstrap
phyloai tree ml fasttree --msa-dir ./trimmed --boot 0

# Fast mode, no gamma, JTT model (AA default)
phyloai tree ml fasttree --msa-dir ./trimmed --mode fastest --no-gamma
```

### Command Hierarchy

```
phyloai tree (click.Group)
└── ml (click.Group)          # "Maximum-likelihood tree inference"
    ├── fasttree              # FastTree backend
    └── iqtree                # IQ-TREE3 backend (future)
```

---

## 3. Parameter Specification

### 3.1 Shared Parameters (`tree ml` level)

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--msa-dir` | Path | — | Directory of MSA files for batch gene tree inference. Mutually exclusive with `--matrix`. |
| `--matrix` | Path | — | Single concatenated matrix file for supermatrix tree inference. Mutually exclusive with `--msa-dir`. |
| `--seq-type` | `AA\|NT\|auto` | auto | Molecule type. |
| `--model` | domain varies | lg (AA) / gtr (NT) | Substitution model. AA: lg\|wag\|jtt. NT: gtr\|jc. |
| `--mode` | `normal\|fastest\|slow` | normal | Speed/accuracy trade-off. |
| `--boot` | int ≥ 0 | 1000 | Bootstrap replicates for node support. 0 disables support (`-nosupport`). |
| `--output-dir` / `-o` | Path | `runs/tree/ml/fasttree` | Output directory. |
| `--threads` / `-t` | int ≥ 1 | 4 | Parallel gene tree workers. Only used in `--msa-dir` mode. |
| `--overwrite` | flag | False | Delete and recreate output directory. |
| `--resume` | flag | False | Resume from checkpoint. `--msa-dir` mode only. |
| `--dry-run` | flag | False | Show commands without executing. |
| `--quiet` / `-q` | flag | False | Suppress terminal output except errors. |
| `--tool-args` | str | None | Extra flags passed verbatim to FastTree. Strategy-only; managed flags blocked (see §4.5). |

### 3.2 FastTree-Specific Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--cat` | int ≥ 1 | 20 | Number of rate categories (`-cat N`). |
| `--gamma` | bool | True | Enable gamma-distributed rate heterogeneity (`-gamma`). False omits the flag. |
| `--fasttree-path` | Path | None | Explicit FastTree executable path. None resolves via ToolEnv (see §12.4). |

---

## 4. FastTree Command Mapping

### 4.1 Model → Flags

| seq_type | model | FastTree flags |
|----------|-------|----------------|
| AA | jtt | *(none — FastTree default)* |
| AA | lg | `-lg` |
| AA | wag | `-wag` |
| NT | jc | `-nt` |
| NT | gtr | `-nt -gtr` |

### 4.2 Mode → Flags

| mode | FastTree flags |
|------|----------------|
| normal | *(none)* |
| fastest | `-fastest` |
| slow | `-slow` |

### 4.3 Other Parameters → Flags

| param | value | FastTree flag |
|-------|-------|---------------|
| `--gamma` | True (default) | `-gamma` |
| `--gamma` | False | *(omit)* |
| `--cat` | N (default 20) | `-cat N` |
| `--boot` | N > 0 (default 1000) | `-boot N` |
| `--boot` | 0 | `-nosupport` |

### 4.4 Full Command Examples

```bash
# AA, LG model, gamma, cat 20, 1000 bootstrap, normal mode
FastTree -lg -cat 20 -gamma -boot 1000 matrix.aa.fa

# NT, GTR model, gamma, cat 20, 1000 bootstrap, normal mode
FastTree -nt -gtr -cat 20 -gamma -boot 1000 matrix.nt.fa

# AA, JTT (default), no gamma, no bootstrap, fastest mode
FastTree -fastest -cat 20 -nosupport matrix.aa.fa
```

The tree output (stdout from FastTree) is captured and written to `<output_dir>/trees/<locus>.tre` or `<output_dir>/<matrix_stem>.tre`. FastTree's stderr output is saved as a per-task log file.

### 4.5 Managed Flags Blocklist for `--tool-args`

Per main design §9.9, PhyloAI manages input, output, work directory, data type, threads, logs, and codon/projection — all other parameters are **strategy parameters** free for `--tool-args` to control.

The following FastTree flags are managed by phyloai and BLOCKED in `--tool-args`:

| Managed flag | Reason |
|-------------|--------|
| `-nt` | Controlled by `--seq-type NT` |
| `-expert`, `-help` | Change interaction mode (always blocked) |
| `/` or `>` in any token | phyloai-managed I/O override |

All other FastTree flags (`-lg`, `-wag`, `-gtr`, `-jc`, `-cat`, `-gamma`, `-boot`, `-nosupport`, `-fastest`, `-slow`, `-noml`, `-spr`, `-mlacc`, `-slownni`, etc.) are **strategy parameters** and pass freely through `--tool-args`.

**Override semantics:** When `--tool-args` contains a flag for a category (model, mode, gamma, cat, boot), phyloai does NOT generate its own version of that flag — `--tool-args` fully takes over that category. For categories NOT mentioned in `--tool-args`, phyloai's defaults apply. This avoids duplicate/conflicting flags in the final FastTree command.

If `--tool-args` contains a blocked flag, exit code 1 with the blocked flag name.

---

## 5. Input Validation

### 5.1 Mutual Exclusivity

`--msa-dir` and `--matrix` are mutually exclusive. Providing both or neither exits with code 1.

### 5.2 `--msa-dir` Mode

- Scan directory for files matching FastTree-compatible extensions: `.fa, .fas, .fasta, .faa, .fna, .phy, .phylip`
- FastTree natively reads FASTA and phylip-relaxed formats (not NEXUS)
- NEXUS files (`.nex, .nxs, .nexus`): skip with warning; users must run `phyloai pretree convert` first
- Skip directories, empty files, unrecognized extensions; record as `skipped` in `data.skipped`
- If exactly 1 valid file: emit WARNING suggesting `--matrix` may be more appropriate, continue
- If zero valid files: exit code 1
- `--threads` controls `ProcessPoolExecutor(max_workers=threads)` for parallel gene tree computation

### 5.3 `--matrix` Mode

- Validate single file exists, is readable, and has a recognized extension (FASTA or PHYLIP)
- `--threads` has no effect
- No checkpoint/resume support

### 5.4 `--seq-type` Auto-Detection and Validation

When `--seq-type auto` (default):

- **`--msa-dir` mode**: scan ALL input files to detect sequence type per file using `detect_seq_type()` from `core/sequence_normalization.py`
- If all files share the same type, proceed with that type
- If mixed types are found: exit code 1 with counts and list of offending files (max 10 shown)
- **`--matrix` mode**: detect from the single file; failure exits code 1

When `--seq-type` is explicit (AA or NT):

- **`--msa-dir` mode**: validate all input files match the declared type
- On mismatch: exit code 1, listing offending files
- **`--matrix` mode**: validate the single file matches; mismatch exits code 1

### 5.5 `--model` Validation

| seq_type | valid values | default |
|----------|-------------|---------|
| AA | jtt, lg, wag | lg |
| NT | jc, gtr | gtr |
| auto | deferred until seq_type resolved | — |

Invalid model for resolved seq_type exits with code 1.

---

## 6. Output Directory Structure

### 6.1 `--msa-dir` (Batch Gene Tree Mode)

```
runs/tree/ml/fasttree/
├── result.json
├── checkpoint.json           # resume support
├── trees/
│   ├── gene001.tre
│   ├── gene002.tre
│   └── ...
└── logs/
    ├── gene001.log           # per-task FastTree stderr
    ├── gene002.log
    └── ...
```

### 6.2 `--matrix` (Single Supermatrix Mode)

```
runs/tree/ml/fasttree/matrix/
├── result.json               # data.tool_stderr inlined (single pattern)
└── <matrix_stem>.tre
```

---

## 7. result.json Schema

```json
{
  "status": "success",
  "command": "phyloai tree ml fasttree --msa-dir trimmed/seqs --seq-type AA --model lg ...",
  "wall_time": 23.5,
  "tool_versions": {"FastTree": "2.1.11"},
  "params": {
    "msa_dir": "/path/to/trimmed/seqs",
    "matrix": null,
    "seq_type": "AA",
    "model": "lg",
    "mode": "normal",
    "boot": 1000,
    "cat": 20,
    "gamma": true,
    "output_dir": "runs/tree/ml/fasttree",
    "threads": 4,
    "overwrite": false,
    "fasttree_path": null,
    "tool_args": null
  },
  "key_results": {
    "n_input": 1066,
    "n_trees": 1050,
    "n_failed": 15,
    "n_skipped": 1,
    "seq_type": "AA",
    "model": "lg",
    "mode": "normal",
    "boot": 1000
  },
  "error": null,
  "data": {
    "summary": {
      "n_input_files": 1066,
      "n_trees": 1050,
      "n_failed": 15,
      "n_skipped": 1,
      "mean_n_taxa": 8.5,
      "mean_wall_time": 0.3,
      "mode": "--msa-dir"
    },
    "cmd": ["FastTree", "-lg", "-cat", "20", "-gamma", "-boot", "1000", "matrix.aa.fa"],
    "tool_stderr": "# single mode: stderr inlined; null for batch",
    "files": [
      {
        "input": "runs/pretree/trim/seqs/gene001.fa",
        "output_tree": "runs/tree/ml/fasttree/trees/gene001.tre",
        "log_file": "runs/tree/ml/fasttree/logs/gene001.log",
        "cmd": ["FastTree", "-lg", "-cat", "20", "-gamma", "-boot", "1000", "gene001.fa"],
        "n_taxa": 8,
        "wall_time": 0.3,
        "warnings": []
      }
    ],
    "failed": [
      {
        "input": "runs/pretree/trim/seqs/bad_gene.fa",
        "reason": "FastTree returned exit code 1",
        "tool_stderr": "Error: ...",
        "wall_time": 0.1
      }
    ],
    "skipped": [
      {"path": "trimmed/seqs/empty.fa", "reason": "empty file"},
      {"path": "trimmed/seqs/data.nex", "reason": "NEXUS format not supported by FastTree; use pretree convert first"}
    ],
    "warnings": []
  }
}
```

---

## 8. Batch Partial-Failure Semantics

### 8.1 Exit Behavior

| Outcome | Exit Code |
|---------|-----------|
| All tasks successful | 0 |
| Partial success (≥1 tree produced, some failed/skipped) | 0 with warnings |
| All tasks failed (0 trees produced after attempted runs) | 2 |
| No valid input files at all | 1 |

### 8.2 Task Status Classification

| Status | Meaning | Criteria |
|--------|---------|----------|
| `success` | Tree produced successfully | FastTree exit code 0, output parses as Newick |
| `failed` | Tool execution error | FastTree exit code non-zero; recorded in `data.failed` |
| `skipped` | Input validation failure | Empty file, unrecognized format, wrong seq_type, NEXUS; recorded in `data.skipped` |

Failed FastTree tasks are stored with `status: "failed"` in checkpoint tasks (not `"skipped"`), so `--resume` retries them. Input-validation skips are always terminal and not retried.

---

## 9. Checkpoint and Resume

Applies to `--msa-dir` mode only (`--matrix` single-file mode does not support resume).

| Field | Value |
|-------|-------|
| `step` | `"tree.ml.fasttree"` |
| `schema_version` | 1 |
| `task_id` | input file stem (logical locus name) |
| `outputs` | `{"tree": "trees/<stem>.tre", "log": "logs/<stem>.log"}` |
| `status` values | `pending`, `running`, `success`, `failed` |

Resume behavior:
- Load `checkpoint.json` → validate params hash → `plan_resume` → re-run `pending|running|failed` tasks
- `success` tasks verified via: (1) output tree file exists and is non-empty, AND (2) file parses as valid Newick via `Bio.Phylo.read()`
- Newick validation failure: task re-marked as `failed` and re-run
- `--resume` and `--overwrite` are mutually exclusive
- Checkpoint flushed at most every 2 seconds during execution; final flush with `fsync=True`

---

## 10. Exit Codes

| Code | Meaning |
|------|---------|
| 0 | Success (partial-failure with warnings also exits 0 per §8.1) |
| 1 | User input error (mutual exclusivity, invalid model, no valid inputs, output dir exists, mixed seq_type, blocked tool-args) |
| 2 | All FastTree runs failed (0 trees produced) |
| 3 | FastTree not found (environment error) |

---

## 11. Warnings

| Condition | Behavior |
|-----------|----------|
| `--msa-dir` contains exactly 1 valid MSA file | WARN: suggest `--matrix` mode, continue |
| `--msa-dir` contains NEXUS files | WARN per file: "NEXUS not supported by FastTree; use pretree convert"; record in `data.skipped` |
| `--msa-dir` contains non-MSA files (unrecognized ext) | Skip silently, record in `data.skipped` |
| Mixed AA/NT types detected in `--msa-dir` | Exit code 1 (no WARN — hard error) |
| FastTree returns non-zero for a gene | Record in `data.failed`, continue batch |
| Partial failure (some genes failed, some succeeded) | WARN at end of run with counts |
| `--threads` passed with `--matrix` | WARN: `--threads` has no effect in single-file mode |

---

## 12. Logging

- `logs/<locus>.log`: per-task FastTree stderr output for batch (`--msa-dir`) mode. Tree output (FastTree stdout) is captured and written to `.tre` files, not duplicated in logs.
- **Single mode (`--matrix`):** FastTree stderr is inlined in `result.json` as `data.tool_stderr` (single pattern, JSON Output Standard Section 5.2). No external log file is written.
- No top-level `fasttree.log` — wall time, exit code, and summary counts are in `result.json`.

---

## 13. Implementation Notes

### 13.1 Files to Create

| File | Purpose |
|------|---------|
| `phyloai/cli/commands/tree.py` | CLI click.Group for `tree` → `ml` → `fasttree`/`iqtree` |
| `phyloai/tree/__init__.py` | Package init |
| `phyloai/tree/ml.py` | `run_fasttree()` library function |
| `phyloai/tree/checkpoint_helpers.py` | Checkpoint build/mark/plan_resume for tree step |
| `docs/commands/tree-ml.md` | User-facing command documentation |

### 13.2 Files to Modify

| File | Change |
|------|--------|
| `phyloai/cli/main.py` | Add `cli.add_command(tree)` |
| `docs/superpowers/specs/2026-06-17-phyloai-tree-design.md` | Update CLI examples, parameter names (`--msa` → `--matrix`), default output dir |

### 13.3 Key Patterns to Follow

- **CLI layer**: thin wrapper — validates params, delegates to library, writes result.json, renders Rich summary
- **Library layer**: `run_fasttree()` accepts all params + `progress_callback`, validates preconditions, scans inputs, runs with `ProcessPoolExecutor`, returns payload dict
- **Tool resolution**: see §13.4
- **Checkpoint**: reuse `core/checkpoint.py` dataclasses; build initial checkpoint, mark tasks, flush throttled
- **Progress bar**: Rich `Progress` with `transient=True`, total from input scan or checkpoint resume count
- **Batch parallelism**: `ProcessPoolExecutor(max_workers=threads)`, one worker per gene, worker runs FastTree as subprocess with 1 internal thread

### 13.4 FastTree Executable Resolution

Resolution order:

1. If `--fasttree-path` is provided: validate the path exists and is executable. If not: exit code 1.
2. If `--fasttree-path` is None: construct `ToolEnv(tool_paths={})` and call `require("FastTree")`. This resolves via: bundled (`phyloai/bundled/`) → PATH (`shutil.which("FastTree")` → `shutil.which("fasttree")`).
3. If ToolEnv also fails to resolve: exit code 3 (environment error).
4. Version detection: attempt `FastTree -version 2>&1` or `FastTree 2>&1` and extract version string via regex.

### 13.5 Tree Output Validation (Resume)

Resume mode validates completed tree outputs using:
1. File exists and is non-empty
2. `Bio.Phylo.read(path, "newick")` succeeds (raises no exception)
3. If parsing fails → task re-marked as `failed` and re-run on resume

---

## 14. Acceptance Criteria

Before merging, verify the following:

### 14.1 CLI Validation
- [ ] `--msa-dir` and `--matrix` together → exit 1
- [ ] Neither `--msa-dir` nor `--matrix` → exit 1
- [ ] Invalid `--model` for seq_type → exit 1
- [ ] Invalid `--seq-type` value → exit 1
- [ ] `--threads < 1` → exit 1
- [ ] `--overwrite` and `--resume` together → exit 1

### 14.2 Input Scanning
- [ ] `--msa-dir` with 0 valid files → exit 1
- [ ] `--msa-dir` with exactly 1 valid file → WARNING, continue
- [ ] `--msa-dir` with NEXUS files → WARNING per file, recorded in `data.skipped`
- [ ] `--msa-dir` with mixed AA/NT types in auto mode → exit 1 with file list
- [ ] `--matrix` with nonexistent file → exit 1

### 14.3 `--tool-args` Blocking
- [ ] `--tool-args "-nt"` → exit 1, blocked managed data-type flag `-nt`
- [ ] `--tool-args "-expert"` → exit 1
- [ ] `--tool-args "-lg"` → accepted (strategy parameter, overrides `--model` default)
- [ ] `--tool-args "-boot 500"` → accepted (strategy parameter, overrides `--boot` default)
- [ ] `--tool-args "-fastest"` → accepted (strategy parameter, overrides `--mode` default)
- [ ] Valid strategy arg (e.g., `--tool-args "-spr 4"`) → appended to command

### 14.4 Batch Execution
- [ ] All genes succeed → exit 0, all trees in `data.files`
- [ ] Some genes fail → exit 0 with warnings, `data.failed` populated, `data.files` for successes
- [ ] All genes fail → exit 2
- [ ] Per-gene stderr captured to `logs/<locus>.log`

### 14.5 Single Supermatrix
- [ ] `--matrix` mode produces single `.tre` file
- [ ] `--threads` with `--matrix` → WARNING, no effect

### 14.6 Dry Run
- [ ] `--dry-run` prints all FastTree commands without executing
- [ ] `--dry-run` creates no files

### 14.7 Resume
- [ ] Successfully completed tasks not re-run
- [ ] Failed/pending tasks re-run
- [ ] Truncated tree file (unparseable Newick) re-marked `failed` and re-run
- [ ] `result.json` matches final run state after resume

### 14.8 Output
- [ ] `result.json` written with correct schema
- [ ] `tool_versions` populated with key `FastTree`
- [ ] `key_results` includes `n_input`, `n_trees`, `n_failed`, `n_skipped`, `seq_type`, `model`, `mode`, `boot`

---

## 15. Relationship to IQ-TREE3 Backend

`phyloai tree ml iqtree` will share the same CLI group (`tree ml`) with shared parameters (`--msa-dir`, `--matrix`, `--seq-type`, `--model`, `--mode`, `--boot`, `--output-dir`, `--threads`, `-q`, `--overwrite`, `-h`, `--tool-args`). IQ-TREE3-specific parameters (e.g., `--partition`, `--merit`, `-bb`) will be documented in a separate spec.
