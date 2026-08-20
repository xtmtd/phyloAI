# phyloai report

[English](report.md) | [中文](report.zh.md)


## Purpose

Generate a reproducible, auditable analysis report from a PhyloAI run directory. Produces two files:

- **`report.json`** — machine-readable record of every step (`params`, `key_results`, `tool_versions`, `methods_text`, output file paths). Serves as the entry point for AI/MCP diagnostics and reproducibility audits.
- **`report.html`** — human-readable report with embedded figures, sortable tables, inline CSV data, and a journal-ready Methods paragraph draft with one-click copy.

A single invocation covers the entire run directory; no sub-commands are needed.

## Usage

```bash
phyloai report --run-dir <run-dir> [OPTIONS]
```

## Inputs

`--run-dir` must contain one or more PhyloAI `result.json` files.

## Outputs

The command writes `report.json` and `report.html` under the report output directory.

## Quick Start

```bash
# Report on a single pipeline run (phyloai run output)
phyloai report --run-dir ./runs/run/faa

# Report on a module-level run (e.g. all pretree steps)
phyloai report --run-dir ./runs/pretree

# Overwrite an existing report
phyloai report --run-dir ./runs/run/faa --overwrite

# Custom output location
phyloai report --run-dir ./runs/pretree -o ./my-report
```

## Input Requirements

`--run-dir` must be a directory containing one or more `result.json` files produced by PhyloAI commands. The report auto-detects two structures:

| Structure | Detection | Typical use |
|-----------|-----------|-------------|
| **pipeline** | `run-dir/result.json` exists AND subdirectories also contain `result.json` | `phyloai run` output |
| **module** | No top-level `result.json`, but subdirectories contain `result.json` | `phyloai pretree`, `phyloai tree`, etc. |

Step discovery is purely filesystem-based — directory scanning with `report/`, `logs/`, and hidden directories excluded.

## Report Content

The HTML report has five panels:

### Panel A — Run Summary

Summary cards showing step counts (succeeded/failed), total wall time, and a progress bar for pipeline runs. Failed steps are listed by name.

### Panel B — Methods

A journal-ready Methods paragraph with one sentence group per analytical step. Each step has a clickable `[step_id]` badge linking to its Step Detail card. A "Copy to clipboard" button copies the plain text for manuscript use.

### Panel C — Steps Detail

One collapsible card per step showing:
- Status indicator and wall time
- `↑ Methods` back-link
- Scientific parameters table
- Key results table
- Embedded CSV tables (≤200 rows, ≤500 KB) with sortable columns
- Warnings and error messages
- Full CLI command

### Panel D — Figures

All PDF/PNG figures produced by the analysis commands, embedded natively. Vector PDF figures preserve quality. Numbered by analytical phase: `Fig-3.x` (pretree), `Fig-4.x` (tree), `Fig-5.x` (posttree).

### Panel E — Output Files Index

A sortable table listing every output file across all steps. Large tables (>20 rows) are collapsible. Each entry links to the actual file.

## Incomplete Runs

The report always succeeds, even when steps have failed:

- Failed steps are included with full error details and expanded by default
- Failed steps have empty `methods_text` and are excluded from the Methods paragraph
- Pipeline status is `"partial"` when some steps failed, `"failed"` when all failed

## Options

| Flag | Default | Description |
|------|---------|-------------|
| `--run-dir PATH` | *(required)* | Run directory to report on. |
| `-o, --output-dir PATH` | `<run-dir>/report` | Output directory for report files. |
| `--overwrite` | off | Overwrite existing report files. |
| `-q, --quiet` | off | Suppress terminal output except errors. |

## Output Structure

```
<run-dir>/report/
├── report.json    # Machine-readable source of truth
└── report.html    # Self-contained HTML report (no external dependencies)
```

## report.json Schema

```json
{
  "phyloai_version": "0.1.0",
  "generated_at": "2026-06-27T14:23:00Z",
  "run_dir": "/abs/path/runs/pretree",
  "run_mode": "module",
  "status": "complete",
  "pipeline_summary": {
    "status": "complete",
    "n_steps_total": 5,
    "n_steps_success": 5,
    "n_steps_failed": 0,
    "n_steps_skipped": 0,
    "total_wall_time": 142.3
  },
  "steps": [
    {
      "step_id": "pretree.align",
      "command": "phyloai pretree align --seq-dir ./raw --method linsi",
      "status": "success",
      "wall_time": 31.4,
      "tool_versions": {"mafft": "7.526"},
      "params": {"method": "linsi", "seq_type": "AA", "threads": 8},
      "key_results": {"n_aligned": 100, "n_skipped": 0},
      "methods_text": "Multiple sequence alignments were performed...",
      "output_files": {},
      "warnings": [],
      "error": null
    }
  ],
  "methods_paragraph": "Multiple sequence alignments were performed...",
  "methods_blocks": [
    {"step_id": "pretree.align", "text": "Multiple sequence alignments...", "step_index": 0}
  ],
  "figures_index": [
    {
      "figure_id": "Fig-3.1",
      "step_id": "pretree.metrics",
      "label": "correlation_heatmap",
      "caption": "Correlation heatmap",
      "path": "/abs/path/correlation_heatmap.pdf",
      "type": "pdf"
    }
  ],
  "tables_index": [
    {
      "table_id": "Table-3.1",
      "step_id": "pretree.metrics",
      "label": "metrics_table",
      "caption": "Phylogenetic informativeness metrics per locus",
      "path": "/abs/path/metrics.csv",
      "type": "csv"
    }
  ]
}
```

## Key Results Enrichment

The report automatically enriches `key_results` from data that some commands place outside the standard `key_results` field:

- `data.summary` scalars (int, float, str) are merged into `key_results`
- `data.*` top-level scalars are merged (e.g., stats single-file mode)
- Nested numeric dicts are flattened: `{length_before: {mean: 10}}` → `length_before_mean`
- Concat-specific metrics (`gap_ratio`, `pi_ratio`) are extracted from `data.variant_stats[0]`

This ensures templates always have complete data regardless of which module produced the `result.json`.

## Exit Codes

| Code | Meaning |
|------|---------|
| `0` | Report generated successfully (even if run had failed steps). |
| `1` | User input error — invalid `--run-dir`, no `result.json` found, or report files exist without `--overwrite`. |

## Warnings and Errors

| Condition | Behaviour |
|-----------|-----------|
| `--run-dir` not a valid run directory | Exit 1; "No result.json found." |
| Report files already exist | Exit 1; use `--overwrite`. |
| Step `result.json` contains `status: "error"` | Step included with error details; Methods paragraph excludes it. |
| Corrupt or unreadable `result.json` | Step recorded as error; report continues. |

## Examples

```bash
# Report on a completed pipeline run
phyloai report --run-dir ./runs/run/faa

# Report on pretree module (multi-step)
phyloai report --run-dir ./runs/pretree

# Report on a single command (e.g. 1-convert only)
phyloai report --run-dir ./runs/pretree/1-convert

# Regenerate report with updated code
phyloai report --run-dir ./runs/run/faa --overwrite

# Quiet mode (only errors to stderr)
phyloai report --run-dir ./runs/pretree -q

# Custom output path
phyloai report --run-dir ./runs/tree -o ./documents/methods
```

## Notes

- `report.html` is fully derived from `report.json` and can be re-rendered at any time without re-scanning the run directory.
- Methods text is deterministic (Python template functions, no LLM). All scientifically meaningful parameters are described; technical parameters (threads, paths, flags) are omitted.
- PDF figures are embedded via `<object>` tags preserving vector quality. No external dependencies (fonts, CDNs, JavaScript libraries) are required.
- The report is designed to be the primary entry point for AI/MCP diagnostics — `report.json` aggregates all step records, parameters, key results, and figure paths into a single queryable document.
- Small CSV/TSV tables (≤200 rows, ≤500 KB) are embedded inline in Step Detail cards as sortable HTML tables.
