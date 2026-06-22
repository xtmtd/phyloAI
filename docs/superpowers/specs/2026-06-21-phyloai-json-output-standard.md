# PhyloAI JSON Output Standard

**Date:** 2026-06-21  
**Status:** Approved for implementation  
**Parent spec:** `2026-06-07-phyloai-design.md` Section 9.4

---

## 1. Purpose

This document defines the binding `result.json` conventions for all PhyloAI non-`doctor` commands. It expands the schema skeleton in the parent design (Section 9.4) with exact field semantics, structural patterns, per-module requirements, and testing assertions.

---

## 2. Top-level Schema (recap)

From parent design Section 9.4:

```json
{
  "status": "success | error",
  "command": "phyloai pretree align --seq-dir ./raw --method linsi --threads 8",
  "wall_time": 142.3,
  "tool_versions": {"mafft": "7.520"},
  "params": {},
  "key_results": {},
  "error": null,
  "data": {}
}
```

All fields are required and genre-neutral except `data`, which is command-specific.

---

## 3. `command` Field

`command` MUST be the **full CLI invocation string** that can be re-executed verbatim to reproduce the identical run. It includes ALL parameters, including those with default values, to ensure reproducibility.

**Hard requirements:**
1. Starts with `"phyloai "` followed by the full subcommand path
2. Includes every parameter that affected the run (user-specified and resolved defaults)
3. Uses **exact CLI flag names** as defined in `click.option` decorators — no invented flag names
4. All non-flag arguments use their resolved (post-validation) values

Correct:  `"phyloai pretree align --seq-dir ./raw --method linsi --seq-type AA --threads 8"`
Incorrect: `"phyloai pretree align"` (missing arguments)
Incorrect: `"phyloai pretree concat --to-format fasta ..."` (flag is `--to`, not `--to-format`)

This field is the primary input for `report methods` paragraph generation.

---

## 4. `params` Field

`params` contains **every parameter** accepted by the command, including those with default values, in their resolved form. This ensures `result.json` is a complete record of the invocation.

**Hard requirements:**
1. **Completeness**: MUST include every Python parameter name from the `run_*` function signature. Omitted parameters indicate a bug.
2. **Naming**: Keys MUST match the Python parameter names in the `run_*` function signature (e.g., `seq_type`, not `seqType` or `sequence-type`).
3. **Resolved values**: Store the post-validation value (e.g., `"AA"` after auto-detection, not `"auto"` for `seq_type`). For commands that resolve parameters per-file (e.g., `pretree convert` where `seq_type` is auto-detected per input, or `pretree metrics` with per-locus metrics), the run-level `params` stores the user-supplied invocation value (`None` or `"auto"` is acceptable); per-file resolved values belong in `data.files[]`.
4. **Null handling**: Parameters not used in this invocation (e.g., mutually exclusive options) MUST be `null`, not absent.
5. **Single definition**: Each module MUST define its `params` dict in exactly one place (the main `run_*` function). CLI handler error paths and sub-functions MUST reuse this dict, not create parallel copies.

**Pattern for implementation:**
```python
# Each run_* function builds params from its own local variables — once.
params = {
    "msa_dir": str(msa_dir),
    "output_dir": str(output_dir),
    "seq_type": resolved_seq_type,
    "threads": threads,
    "tool_args": tool_args,
    "overwrite": overwrite,
    "dry_run": dry_run,
    "quiet": quiet,
}
# The same dict is used in all success/error/dry_run return paths.
# The _build_*_command function accepts the same parameters to guarantee flag-name consistency.
```

Parameter order is unspecified; consumers must key by name.

---

## 5. `data` Field — Standard Patterns

While `data` is command-specific, two structural patterns are used.

### 5.1 Batch Pattern

Used by: `align`, `trim`, `metrics`, `filter taper`, `filter metrics`, `filter cluster`, `fasttree` (per-gene), `iqtree` (per-gene).

```json
{
  "data": {
    "summary": { "n_input_files": 100, "n_success": 96, "n_skipped": 4 },
    "files": [
      {
        "input": "raw/gene1.fa",
        "output": "seqs/gene1.fa",
        "cmd": ["mafft", "--maxiterate", "1000", "gene1.fa"],
        "log_file": "logs/gene1.log",
        "status": "success",
        "wall_time": 1.2,
        "warnings": []
      }
    ],
    "failed": [],
    "skipped": []
  }
}
```

| Field | Required? | Description |
|-------|-----------|-------------|
| `files[].cmd` | Required when external tool invoked per task | Exact argv list. Pure-Python batch may omit. |
| `files[].log_file` | Required when `cmd` is present | Path relative to output dir, e.g. `logs/gene1.log` |
| `files[].wall_time` | **Required** when `cmd` is present | Per-task wall time in seconds. MUST be a measured positive float (`> 0`); `0.0` is only valid for dry-run tasks. |

- Do NOT inline tool stderr in `files[]` — stderr lives in `logs/<locus>.log`.
- Pure-Python batch commands (no external tool per task) may omit both `cmd` and `log_file`.

### 5.2 Single Pattern

Used by: `concat`, `msc`, `cf` (all sub-modes), `filter treeshrink`, `filter symtest`, `fasttree` (`--matrix`), `iqtree` (`--matrix`).

```json
{
  "data": {
    "cmd": ["wastral", "-i", "merged.trees", "-o", "wastral.tre"],
    "tool_stderr": "wASTRAL version 1.25.4 ...\nQuartet score: 0.95\n",
    "output": "wastral.tre",
    "warnings": []
  }
}
```

| Field | Required? | Description |
|-------|-----------|-------------|
| `data.cmd` | **MUST** | Exact argv list. If no external tool is invoked, MUST be `[]`. |
| `data.tool_stderr` | **MUST** | Raw tool stderr. If no external tool, MUST be `""`. May also be `""` when full diagnostic output is captured in an external log file referenced by `data.tool_log`. |
| `data.tool_log` | Optional | Path to tool-native report (e.g., IQ-TREE `.iqtree` or `.log`). |

- Single-mode PhyloAI-authored summary log files MUST NOT be written. Tool-native diagnostic files produced as a side effect of execution (e.g., IQ-TREE `.log`, wASTRAL captured output saved as `wastral.log`) are permitted and MUST be referenced via `data.tool_log` when present.
- `data.tool_stderr` stores raw stderr only. Summaries belong in `warnings` or `data.summary`.
- Top-level or sidecar log files MUST NOT duplicate `wall_time`, `exit_code`, `command`, or other fields already in `result.json`.

### 5.3 Tool Diagnostic Output Model (recap from parent Section 6.2)

| Mode | Diagnostic output location | Referenced in JSON? |
|------|---------------------------|---------------------|
| Batch | `<output-dir>/logs/<locus>.log` | `files[].log_file` |
| Single | Inlined in `data.tool_stderr`, or external log referenced by `data.tool_log` when `data.tool_stderr` is `""` | `data.tool_log` |
| Utility | No log written | N/A |

**Output stream capture**: Many tools write diagnostic/progress output to stdout rather than stderr. Per-locus log files and `data.tool_stderr` MUST capture **both stdout and stderr** from the subprocess, concatenated (stdout first, then stderr, separated by a newline when both are non-empty). The field name `tool_stderr` is retained for backward compatibility with the JSON schema; its content is the merged diagnostic output. |

---

## 6. Additional Rules

### 6.1 Filter Subcommands

Filter subcommands fall into three categories:

| Category | Subcommands | Pattern | `data` structure |
|----------|-------------|---------|-----------------|
| **True batch** | `taper` | Batch (§5.1) | `data.files[]` with per-locus `cmd`, `log_file`, `wall_time` |
| **Single-invocation** | `treeshrink`, `symtest` | Single (§5.2) | `data.cmd`, `data.tool_stderr` (inlined or `""` when external `tool_log` present), `data.results[]` for per-locus decisions |
| **Pure-Python** | `metrics`, `cluster` | Batch (§5.1) | `data.files[]` (lightweight, no `cmd`/`log_file`/`wall_time`) |

**Rationale for single-invocation pattern:** `treeshrink` and `symtest` each invoke their external tool exactly once for the entire dataset. The tool produces per-locus results, but there is only one tool stderr stream and one wall-time measurement. Using single pattern with `data.cmd` + `data.tool_stderr` + `data.results` avoids:
- Identical `wall_time` values across all per-locus entries (total / N approximation is misleading)
- Shared `log_file` references duplicated across every entry
- Missing `files[].cmd` (there is no per-locus command to record)

The `data.results[]` array holds per-locus filtering decisions (`locus`, `status`, and optional locus-specific metadata like `output_tree`, `removed_taxa`, `reason`), without per-task timing or command fields.

### 6.2 Pure-Python Batch

Commands without external tool invocations per task (e.g., `pretree metrics`, `pretree filter metrics`, `pretree filter cluster`) omit `files[].cmd` and `files[].log_file`. For `pretree metrics`, `data.files[]` includes lightweight per-locus entries (`input`, `status`) while metric values live in the `metrics.csv` sidecar file. `data.summary` is the canonical location for per-run aggregate stats (`n_markers`, `n_success`, `n_errors`, `warnings`).

### 6.3 `tool_stderr` Integrity

`data.tool_stderr` stores the merged diagnostic output (stdout + stderr) from external tools. Do not embed key=value summaries in this field. If a human-readable summary is warranted, place it in `warnings` or `data.summary`. Log files (`logs/<locus>.log`) follow the same merged-output convention.

---

## 7. Per-Module Field Requirements

| Module | `command` full? | `params` full? | `data.cmd` | `data.tool_stderr` | `data.tool_log` | `files[].cmd` | `files[].log_file` | `files[].wall_time` |
|--------|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| convert | ✓ | ✓ | — | — | — | — | — | — |
| stats | ✓ | ✓ | — | — | — | — | — | — |
| align | ✓ | ✓ | — | — | — | ✓ | ✓ | ✓ |
| trim | ✓ | ✓ | — | — | — | ✓ | ✓ | ✓ |
| metrics | ✓ | ✓ | — | — | — | — | — | — |
| filter taper | ✓ | ✓ | — | — | — | ✓ | ✓ | ✓ |
| filter treeshrink | ✓ | ✓ | ✓ | ✓ | — | — | — | — |
| filter symtest | ✓ | ✓ | ✓ | ✓ | — | — | — | — |
| filter metrics | ✓ | ✓ | — | — | — | — | — | — |
| filter cluster | ✓ | ✓ | — | — | — | — | — | — |
| concat | ✓ | ✓ | ✓ ([] ok) | ✓ | — | — | — | — |
| fasttree (batch) | ✓ | ✓ | — | — | — | ✓ | ✓ | ✓ |
| fasttree (single) | ✓ | ✓ | ✓ | ✓¹ | — | — | — | — |
| iqtree (batch) | ✓ | ✓ | — | — | — | ✓ | ✓ | ✓ |
| iqtree (single) | ✓ | ✓ | ✓ | ✓¹ | — | — | — | — |
| msc | ✓ | ✓ | ✓ | ✓² | ✓ | — | — | — |
| cf (gCF/sCF/sCFl) | ✓ | ✓ | ✓ | ✓ | ✓ (IQ-TREE native .log) | — | — | — |
| cf (qCF) | ✓ | ✓ | ✓ | ✓² | ✓ | — | — | — |

`—` = not required. `opt` = optional (pure-Python batch). `[] ok` = empty array valid when no external tool invoked.
`✓¹` = MAY be `""` when per-locus diagnostic output is captured in `output/<locus>.log` (single-mode FastTree/IQ-TREE).
`✓²` = MAY be `""` when full diagnostic output is in an external log file referenced by `data.tool_log` (wASTRAL-based commands).

**`params` full?** means the `params` dict MUST include every Python parameter name from the `run_*` function signature, with its resolved value. No parameter may be silently omitted. This is the single most common source of bugs in PhyloAI JSON output — manual dict construction without a shared source of truth inevitably leads to drift.

---

## 8. Testing Assertions

Every module's tests MUST verify structural compliance. Minimally:

```python
# command is a full CLI string with args
assert result["command"].startswith("phyloai ")
assert len(result["command"].split()) >= 3

# params completeness: NO param from run_* signature may be absent
# (module-specific: enumerate expected keys from function signature)
expected_params = {"msa_dir", "output_dir", "seq_type", "threads", "overwrite", "dry_run", "quiet", ...}
absent = expected_params - set(result["params"].keys())
assert not absent, f"params missing keys: {absent}"

# command reproducibility: flag names must match click.option definitions
# (re-executing command must not produce "no such option" error)
# Example: concat uses --to, never --to-format

# single-mode commands
assert isinstance(result["data"].get("cmd", []), list)
assert isinstance(result["data"].get("tool_stderr", ""), str)
if "tool_log" in result["data"]:
    assert isinstance(result["data"]["tool_log"], str)

# batch-mode commands with external tools
for f in result["data"].get("files", []):
    if "cmd" in f:
        assert isinstance(f["cmd"], list)
        assert f.get("wall_time", 0) > 0 or f.get("status") == "dry_run", \
            f"wall_time must be > 0 for non-dry-run tasks, got {f.get('wall_time')}"
    if "log_file" in f:
        assert f["log_file"].startswith("logs/")

# tool_stderr must not duplicate fields already in result.json
stderr = result["data"].get("tool_stderr", "")
assert "wall_time" not in stderr and "exit_code" not in stderr and "command" not in stderr
```

**Params regression prevention:** The `validate_result_json` helper in `tests/helpers.py` should be extended per-module to verify params completeness. Each test module defines its expected params key set from the `run_*` function signature and passes it to the validator.

---

## 9. Relationship to Other Specs

- **Parent design** (`2026-06-07-phyloai-design.md`): Defines the schema skeleton and output directory conventions. This doc interprets and expands those rules.
- **Per-subcommand specs**: Each defines its own `key_results` and command-specific `data` content, but must conform to the batch or single structural pattern defined here.
- **checkpoint spec** (`2026-06-12-checkpoint-resume-design.md`): Checkpoint format is independent; this doc does not modify it.
