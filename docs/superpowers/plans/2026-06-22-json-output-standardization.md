# PhyloAI JSON/Log Output Standardization Plan

**Date:** 2026-06-22
**Parent spec:** `docs/superpowers/specs/2026-06-21-phyloai-json-output-standard.md`
**Design ref:** `docs/superpowers/specs/2026-06-07-phyloai-design.md` Section 9.4

## Goal

Align all 13 modules' `result.json` payload and tool stderr handling with the JSON Output Standard. Primarily payload/log standardization — changes are concentrated in the result assembly phase. Limited execution-path changes occur where stderr is not currently captured (e.g., `treeshrink`, `symtest` currently discard tool stderr and need subprocess stderr capture added). Core computation logic is untouched.

## Key Rules (from spec)

| Aspect | Batch (external tools) | Batch (pure Python) | Single |
|--------|------|-----|--------|
| `command` | Full CLI string | Full CLI string | Full CLI string |
| stderr location | `logs/<locus>.log` | N/A | Inlined `data.tool_stderr` |
| `files[].cmd` | ✓ (argv list) | Omit | N/A |
| `files[].log_file` | ✓ (path to log) | Omit | N/A |
| `data.cmd` | N/A | N/A | ✓ (argv list, `[]` if none) |
| `data.tool_stderr` | Omit | Omit | ✓ (raw stderr str) |
| Top-level `.log` | **DELETE** (PhyloAI-authored summary) | **DELETE** (PhyloAI-authored summary) | **DELETE** (PhyloAI-authored summary) |

## Changes Summary

| # | Module | File(s) changed | Changes |
|---|--------|-----------------|---------|
| 1 | `convert` | `pretree/convert.py`, `cli/commands/pretree.py` | Full `command` |
| 2 | `stats` | `cli/commands/pretree.py` | Full `command` |
| 3 | `metrics` | `pretree/metrics.py` | Full `command`, delete `metrics.log` |
| 4 | `concat` | `pretree/concat.py` | Full `command`, delete `concat.log` |
| 5 | `msc` | `tree/msc.py` | `data.tool_stderr` (success=`""`, error=real stderr), `data.tool_log`="wastral.log", delete `msc.log` |
| 6 | `cf` | `tree/cf.py` | `data.tool_stderr` (qCF success=`""`, else real stderr), `data.tool_log` (IQTREE: native `.log`; qCF: `wastral.log`) |
| 7 | `fasttree` | `tree/ml.py` | Batch: per-locus `logs/<locus>.log` + `files[].cmd`/`log_file`. Single: `data.cmd` + `data.tool_stderr`=`""` (log in `output/<stem>.log`). Delete `fasttree.log`. |
| 8 | `iqtree` | `tree/ml_iqtree.py` | Same as fasttree. Delete `iqtree.log`. |
| 9 | `align` | `pretree/align.py`, `cli/commands/pretree.py` | Full `command`, per-locus `logs/<locus>.log`, `files[].cmd`/`log_file`, delete `align.log` |
| 10 | `trim` | `pretree/trim.py`, `cli/commands/pretree.py` | Full `command`, `per_gene`→`files`, per-locus `logs/<locus>.log`, `files[].cmd`/`log_file`, delete `trim.log` |
| 11 | `filter taper` | `pretree/filter.py` | Full `command`, restructure to `files[]`, per-locus `logs/<locus>.log` |
| 12 | `filter treeshrink` | `pretree/filter.py` | Full `command`, single-invocation: `data.cmd` + `data.tool_stderr` + `data.results[]` |
| 13 | `filter symtest` | `pretree/filter.py` | Full `command`, single-invocation: `data.cmd` + `data.tool_stderr` + `data.results[]` |
| 14 | `filter metrics/cluster` | `pretree/filter.py` | Full `command`, delete `filter.log` |

---

## Task 1: `pretree convert` — Full `command` field

**Files:**
- Modify: `phyloai/pretree/convert.py` (payload assembly, lines 60–91)
- Modify: `phyloai/cli/commands/pretree.py` (convert_command, lines ~73–130)

**Current state:**
- `command` is partial (e.g., `"phyloai pretree convert"`)
- `_write_result_json` at CLI layer writes the payload

**Action:**
1. In `convert_input()`, build full `command` string from all resolved function parameters (input_path, output_dir, to_format, overwrite, seq_type, etc.)
2. Verify CLI layer `convert_command` passes enough context to reconstruct the string — if not, build it in CLI and include in payload

**No log changes** — convert already has no .log files (correct spec).

---

## Task 2: `pretree stats` — Full `command` field

**Files:**
- Modify: `phyloai/cli/commands/pretree.py` (`_run_stats_command`, lines ~275–381)

**Current state:**
- `command` is partial or missing
- All payload built inline in CLI function

**Action:**
1. Build full `command` string from all resolved Click parameters
2. Include in result.json payload before writing

**No log changes** — stats already has no .log files (correct spec).

---

## Task 3: `pretree metrics` — Full `command`, delete `metrics.log`

**Files:**
- Modify: `phyloai/pretree/metrics.py`

**Current state:**
- `command` is partial
- `_write_log()` writes `<output_dir>/metrics.log` at line 810–813
- `_write_result_json()` writes result.json at line 804–807
- Already uses batch `files[]` pattern; no `cmd` or `log_file` per spec (pure Python)

**Action:**
1. Build full `command` string from all resolved parameters (`--msa-dir`, `--output-dir`, `--seq-type`, etc.)
2. Delete the `_write_log()` function call (or the whole function) — NO replacement; pure Python batch has no log
3. Verify `files[]` entries do NOT include `cmd` or `log_file` (should already be correct)

---

## Task 4: `pretree concat` — Full `command`, delete `concat.log`

**Files:**
- Modify: `phyloai/pretree/concat.py` (payload assembly, lines ~751–820)

**Current state:**
- `command` = `"phyloai pretree concat --msa-dir <path>"` (very minimal)
- Writes `concat.log` with key=value summary (lines 804–819)
- Single pattern: `data` contains flat stats dict, no `data.cmd`, no `data.tool_stderr`

**Action:**
1. Build full `command` string from all resolved parameters (include `--output-dir`, `--prefix`, `--seq-type`, `--taxa-occupancy`, `--recoding`, `--outgroup`, `--to`, `--translate-codon`, `--exclude-codon3`, `--dry-run`, `--overwrite`)
2. Delete the `concat.log` writing block (lines 804–819)
3. Add `"cmd": []` to `data` (concat is pure Python, single pattern → empty argv list per spec)
4. Add `"tool_stderr": ""` to `data` (no external tool → empty string)

---

## Task 5: `tree msc` — Inline `data.tool_stderr`, delete log files

**Files:**
- Modify: `phyloai/tree/msc.py` (`_assemble_wastral_result`, lines ~273–389)
- Modify: `phyloai/cli/commands/tree.py` (msc_command, lines ~844–924)

**Current state:**
- `command` already near-complete (reconstructed from args) ✓
- `proc.stderr` written to `<output_dir>/wastral.log` (line 596)
- `msc.log` written for non-wastral errors (lines ~659–671)
- `data` has `cmd` (subprocess argv list) but NO `data.tool_stderr` inlined
- Single pattern — `data` contains `input`, `output_tree`, etc.

**Action:**
1. Add `"tool_stderr": ""` to `data` in `_assemble_wastral_result()` on success (diagnostic output is in `wastral.log`); keep real stderr on error paths.
2. Write `proc.stderr` to `<output_dir>/wastral.log` and reference via `"tool_log": "wastral.log"` in `data`. This follows the same pattern as IQ-TREE cf modes which reference their native `.log` files.
3. Remove the `msc.log` write (lines ~659–671) — any error info goes into `data.warnings`.
4. Verify `data.cmd` is present (already is) with full argv list.

---

## Task 6: `tree cf` — Inline `data.tool_stderr`, add `data.tool_log`, delete log files

**Files:**
- Modify: `phyloai/tree/cf.py` (`_assemble_cf_result`, lines ~516–627)

**Current state:**
- Single pattern for all modes (gCF, sCF, sCFl use IQ-TREE; qCF uses wASTRAL)
- `proc.stderr` saved to `wastral.log` for qCF (line ~971)
- `_write_cf_log()` writes `cf.log` with status summary (lines ~480–511)
- `data` has `cmd` (subprocess argv list) but NO `data.tool_stderr` inlined
- IQ-TREE writes native `.log` file (e.g., `gCF.log`, `sCF.log`, `sCFl.log`) as side effect

**Action:**
1. Add `"tool_stderr": proc.stderr` to `data` in `_assemble_cf_result()` for IQTREE modes; for qCF mode, `"tool_stderr": ""` (diagnostic output is in `wastral.log`).
2. For gCF/sCF/sCFl modes: add `"tool_log": "<prefix>.log"` to `data` (references IQ-TREE's own `.log` output — e.g., `"tool_log": "gCF.log"`). This file is saved by IQ-TREE, not by PhyloAI.
3. For qCF mode: add `"tool_log": "wastral.log"` (wASTRAL captured stderr saved as external log, same pattern as msc.py).
4. Write `proc.stderr` to `<output_dir>/wastral.log` for qCF mode only.
5. Delete `_write_cf_log()` function and its call — NO separate `cf.log`

---

## Task 7: `tree ml fasttree` — Per-locus logs, batch/single `data` restructure

**Files:**
- Modify: `phyloai/tree/ml.py` (payload assembly `_assemble_result`, lines ~646–783)

**Current state:**
- `command` already near-complete ✓
- **Batch mode (`--msa-dir`):**
  - `data.files[]` entries currently include `tool_stderr` inline in JSON — WRONG per spec
  - `data.files[]` entries include `cmd` but as a string, not list — needs to be list
  - One `fasttree.log` summary file written at lines 771–782
  - Per-locus stderr already saved to `logs/<locus>.log` at line 215 in `_run_fasttree()` — BUT `tool_stderr` still duplicated in the file result dict at line 222/235
  - Has `n_skipped`, `n_resume_skipped`, `n_failed` — fine
- **Single mode (`--matrix`):**
  - Result dict built inline (not via `_assemble_result()`) — need to check where
  - `data` should follow single pattern with `data.cmd`, `data.tool_stderr`

**Action (batch):**
1. In `_assemble_result()`: remove `tool_stderr` from each `files[]` entry — it's already in `logs/<locus>.log`
2. Add `"log_file": f"logs/{stem}.log"` to each `files[]` entry (path matches where `_run_fasttree` writes it at line 215)
3. Convert `cmd` from string to list (`shlex.split(cmd_str)`) in each `files[]` entry
4. Delete the `fasttree.log` summary write block (lines 771–782)
5. Verify `skipped` and `failed` sections are in `data`, not in `files[]`

**Action (single):**
1. For single-mode (`--matrix`): ensure `data` contains `"cmd": [...argv...]`, `"tool_stderr": ""` (per-locus log already written to `output/<stem>.log`), `"output": "<path>"`, `"warnings": [...]`
2. Per spec §5.2, `data.tool_stderr` MAY be `""` when diagnostic output is captured in an external log file (`output/<stem>.log`). No `tool_log` reference needed here — the log path is the output tree's sibling file.

---

## Task 8: `tree ml iqtree` — Same pattern as fasttree

**Files:**
- Modify: `phyloai/tree/ml_iqtree.py` (`_assemble_iqtree_result`, lines ~1505–1715)

**Current state:**
- `command` already near-complete ✓
- Same dual batch/single pattern as fasttree
- `iqtree.log` written at lines 1697–1715
- Batch: `data.files[]` entries may include `tool_stderr` inline — needs removal
- Batch: Per-locus stderr written to `logs/<locus>.log` at line 738
- Has `failed`, `skipped` sections ✓

**Action (batch):**
1. Remove `tool_stderr` from each `files[]` entry
2. Add `"log_file": f"logs/{stem}.log"` to each `files[]` entry
3. Convert `cmd` to list via `shlex.split()`
4. Delete the `iqtree.log` write block

**Action (single):**
1. Ensure single-mode `data` follows: `"cmd": [...]`, `"tool_stderr": ""` (per-locus log already written to `output/<stem>.log`), `"output": "<path>"`, `"warnings": [...]`
2. Same rationale as fasttree — diagnostic output in external log file, `tool_stderr` is `""`.

---

## Task 9: `pretree align` — Per-locus logs, `files[].cmd`/`log_file`, delete `align.log`

**Files:**
- Modify: `phyloai/pretree/align.py`
  - `reconstruct_align_result()` (lines ~236–310) — resume path payload
  - Fresh run payload assembly (lines ~926–987)
  - `_write_align_log()` (lines ~990–1005)
- Modify: `phyloai/cli/commands/pretree.py` (align_command, writes result.json, lines ~618–621)

**Current state:**
- `command` is partial (resume path uses `checkpoint.command`, fresh path hardcodes minimal flags)
- `data.files[]` structure exists but `tool_cmd`/`tool_stderr` are stripped out; NOT in result.json — GOOD
- Per-gene stderr written to monolithic `align.log` via `_write_align_log()` — WRONG per spec
- `files[]` entries have `input`, `output_aa`, `output_nt`, `n_taxa`, `alignment_length`, `wall_time`, `warnings` — needs `cmd` and `log_file` added
- `_align_one()` returns dict with `tool_cmd` (list), `tool_stderr` (str) — available for per-locus log writing

**Action:**
1. Build full `command` string: include `--output-dir`, `--backtrans`, `--overwrite`, `--tool-args`, `--mafft-path`, `--magus-path`, `--resume`, `--seq-type` in addition to `--seq-dir`, `--method`, `--threads` that are already there
2. Update `_write_align_log()` → `_write_per_locus_logs()`: iterate `file_results`, write each `tool_stderr` to `logs/<locus>.log`, create the `logs/` directory if needed
3. Add `"cmd": shlex.split(tool_cmd)` to each `files[]` entry (the `tool_cmd` list from `_align_one()` is already a string, split it)
4. Add `"log_file": f"logs/{locus}.log"` to each `files[]` entry
5. Update both `reconstruct_align_result()` and the fresh-run path consistently
6. Delete `_write_align_log()` (or repurpose it — rename + change behavior) — NO top-level `align.log`

---

## Task 10: `pretree trim` — Rename `per_gene`→`files`, per-locus logs, delete `trim.log`

**Files:**
- Modify: `phyloai/pretree/trim.py`
  - `_build_trim_payload()` (lines ~686–704)
  - `_write_trim_log()` (lines ~749–757)
- Modify: `phyloai/cli/commands/pretree.py` (trim_command, writes result.json, lines ~764–767)

**Current state:**
- `command` is partial (missing `--threads`, `--backtrans`, `--overwrite`, custom paths)
- Batch pattern uses `data.per_gene[]` — should be `data.files[]`
- Per-gene stderr written to monolithic `trim.log` via `_write_trim_log()`
- `per_gene[]` entries have `gene`, `length_before`, `length_after`, `columns_removed`, `outputs` — missing `cmd` and `log_file`
- `_make_success_result()` includes `tool_cmd` (string) and `tool_stderr` (str) — available for per-locus log writing

**Action:**
1. Build full `command` string: add `--threads`, `--backtrans`, `--overwrite`, `--trimAl-path`, `--bmge-path`, `--tool-args`, `--dry-run` to the existing `--msa-dir`, `--tool`, `--seq-type`
2. Rename `data.per_gene` → `data.files` throughout
3. Replace `_write_trim_log()` with per-locus log writing: iterate results, write each `tool_stderr` to `logs/<locus>.log`
4. Add `"cmd": shlex.split(tool_cmd)` to each `files[]` entry
5. Add `"log_file": f"logs/{locus}.log"` to each `files[]` entry
6. Rename `per_gene` fields to standard names: `gene` → `locus` or keep as-is if it's the locus stem name

---

## Task 11: `pretree filter taper` — Restructure to `files[]`, per-locus logs

**Files:**
- Modify: `phyloai/pretree/filter.py`
  - `run_taper()` payload assembly (lines ~440–473)
  - `_write_filter_log()` (lines ~88–98)
  - `_write_result_json()` (lines ~71–74)

**Current state:**
- `command` is partial (`--msa-dir`, `--seq-type`, `--cutoff` only)
- `data` structure is non-standard: `retained_loci[]`, `dropped_loci[]`, `file_results` (full per-locus dicts with inline `tool_stderr`)
- `filter.log` written with header-style format

**Action:**
1. Build full `command`: include all resolved parameters (listing is not exhaustive — every flag in the function's parameter list must appear; e.g., `--output-dir`, `--overwrite`, `--tool-args`)
2. Restructure `data`: 
   - `data.summary`: `n_input`, `n_retained`, `n_dropped`, `total_masked_aa_sites`, `total_masked_taxa`, `masked_loci`
   - `data.files[]`: one entry per locus with `input`, `locus`, `status` (`"retained"`/`"dropped"`), `log_file`, `cmd`, `warnings`, `n_taxa_before`, `n_taxa_after`, `length_before`, `length_after`, `masked_sites`, `wall_time`
   - `data.skipped[]`: entries for loci skipped before TAPER (invalid input, etc.)
3. Per-locus log: write each `tool_stderr` to `logs/<locus>.log` (create `logs/` dir)
4. Each `files[]` entry gets `"cmd": shlex.split(tool_cmd)` and `"log_file": f"logs/{locus}.log"`
5. Delete `_write_filter_log()` call — no top-level `filter.log`
6. Remove old `retained_loci`/`dropped_loci`/`file_results`/`retained_msa_stats` from `data`

---

## Task 12: `pretree filter treeshrink` — Single-invocation pattern with inlined stderr

**Files:**
- Modify: `phyloai/pretree/filter.py` — `run_treeshrink()` payload assembly

**Current state:**
- `command` is full ✓ (all resolved params via `_build_treeshrink_command`)
- Single-invocation pattern per spec §6.1: `data.cmd` (argv list), `data.tool_stderr` (merged stdout+stderr), `data.results[]` (per-locus decisions)
- No `files[]`, no shared log file — correct per spec

**Action:**
1. Already DONE. Treeshrink uses the Single pattern (§5.2 / §6.1): one external tool invocation, merged stderr inlined into `data.tool_stderr`, per-locus filtering decisions in `data.results[]`.
2. No per-locus `cmd` or `log_file` — tool is invoked once for the entire dataset.
3. `data.results[]` entries contain `locus`, `status`, `output_tree`, `removed_taxa`, `warnings`.

Note: Plan originally described a batch-with-shared-log pattern, but spec §6.1 correctly classifies treeshrink as Single-invocation. Code follows spec.

---

## Task 13: `pretree filter symtest` — Same single-invocation pattern as treeshrink

**Files:**
- Modify: `phyloai/pretree/filter.py` — `run_symtest()` payload assembly

**Current state:**
- Single IQ-TREE `--symtest-only` invocation, exactly one stderr output
- `command` is full ✓ (all resolved params via `_build_symtest_command`)
- Single-invocation pattern per spec §6.1: `data.cmd`, `data.tool_stderr` (merged), `data.results[]`

**Action:**
1. Already DONE. Symtest uses the Single pattern (§5.2 / §6.1): one IQ-TREE invocation, merged stdout+stderr inlined into `data.tool_stderr`, per-locus filtering decisions in `data.results[]`.
2. No per-locus `cmd` or `log_file` — tool is invoked once for the entire dataset.
3. `data.results[]` entries contain `locus`, `status`, `symtest_type`, `p_value`, `warnings`.

Note: Plan originally described a batch-with-shared-log pattern, but spec §6.1 correctly classifies symtest as Single-invocation. Code follows spec.

---

## Task 14: `pretree filter metrics` and `pretree filter cluster` — Full `command`, delete `filter.log`

**Files:**
- Modify: `phyloai/pretree/filter.py` — `run_metrics_filter()` and `run_cluster()` payload assembly

**Current state:**
- `command` is partial
- Pure Python — no external tool invocations
- `filter.log` written by shared `_write_filter_log()`

**Action:**
1. Build full `command` strings for each subcommand
2. Remove `_write_filter_log()` call (or conditionally skip it — but since ALL subcommands now drop `filter.log`, remove entirely)
3. Verify `files[]` entries omit `cmd` and `log_file` (pure Python per spec)

---

## Verification

After all modules are updated, run:

```bash
# Unit tests for each module
pytest tests/pretree/test_align.py -v
pytest tests/pretree/test_trim.py -v
pytest tests/pretree/test_metrics.py -v
pytest tests/pretree/test_filter.py -v
pytest tests/pretree/test_concat.py -v
pytest tests/pretree/test_convert.py -v
pytest tests/pretree/test_stats.py -v
pytest tests/tree/test_ml.py -v
pytest tests/tree/test_ml_iqtree.py -v
pytest tests/tree/test_msc.py -v
pytest tests/tree/test_cf.py -v

# Structural compliance (per JSON Output Standard Section 8)
# Every test that parses result.json should verify:
#   - command starts with "phyloai " and has >=3 tokens
#   - batch: files[].log_file starts with "logs/"
#   - single: data.cmd is list, data.tool_stderr is str
#   - single: if "tool_log" in data, assert isinstance(data["tool_log"], str)
#   - no tool_stderr contains wall_time, exit_code, or command

# No PhyloAI-authored summary .log files in output directory root
# (Tool-native .log files referenced via data.tool_log are preserved — e.g., gCF.log, wastral.log)
find runs/ -maxdepth 3 -name "align.log" -o -name "trim.log" -o -name "metrics.log" \
  -o -name "filter.log" -o -name "concat.log" -o -name "fasttree.log" \
  -o -name "iqtree.log" -o -name "msc.log" -o -name "cf.log"
```

## Execution Order

Recommended order (minimizes cross-module coupling risk):

1. **convert, stats** — simplest, no log changes
2. **metrics** — pure Python, no per-locus logs
3. **concat** — single mode, pure Python
4. **msc** — single mode, external tool, tool_stderr="" with external wastral.log
5. **cf** — single mode, external tool, + tool_log
6. **fasttree** — batch + single, per-locus logs
7. **iqtree** — batch + single, per-locus logs
8. **align** — batch, per-locus logs, complex
9. **trim** — batch, per-locus logs, rename per_gene
10. **filter taper** — batch, per-locus logs
11. **filter treeshrink/symtest** — single-invocation, inlined/merged tool_stderr *or* external log via tool_log
12. **filter metrics/cluster** — pure Python cleanup, delete filter.log

The filter module (tasks 11–14) should be done together since they share `_write_filter_log()` and `_write_result_json()`.
