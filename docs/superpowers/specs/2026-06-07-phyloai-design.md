# PhyloAI Design Specification

**Date:** 2026-06-07  
**Last updated:** 2026-06-26 (updated report module design; Section 3, 4.1, 6, 7 revised to match `2026-06-26-phyloai-report-design.md`)  
**Status:** Approved for implementation

---

## 1. Project Identity

**Name:** PhyloAI  
**Tagline:** *An AI-native modular phylogenomics analysis platform*  
**Scope:** MSA-first. Covers alignment through post-tree analysis. Does not handle upstream data acquisition (genome assembly, marker extraction, orthology inference).  
**Target users:** Intermediate-level users with limited theoretical background. AI interaction depth scales with analytical module complexity, not with user expertise level.

---

## 2. Architecture Overview

Three independent layers with unidirectional dependencies:

```
Python Library   ←  core algorithms, importable
      ↓
CLI              ←  built on library, single entry point
      ↓
MCP Server + Skill  ←  built on CLI, added after CLI stabilizes
```

**Development order:** Library + CLI first → MCP Server + Skill after CLI is stable. `phyloai run` (pipeline orchestration) and `report/` (cross-module aggregation) are separate phases: `run` depends only on the analysis modules it orchestrates (Phases 2–3), while `report` must wait for all analysis phases (2–4) to finalize their JSON output schemas.  
**MCP pre-requisite:** All non-`doctor` commands write a standard `result.json` file and use standard exit codes from day one, so MCP wrapping reads structured results from files rather than command-specific stdout formats. `phyloai doctor` is the only command with `--output-format text|json` because it is primarily a human-facing environment diagnostic. Doctor design details are in `docs/superpowers/specs/2026-06-18-phyloai-doctor-design.md`.

---

## 3. Module Structure

```
phyloai/
├── core/
│   ├── env.py          # tool detection, path resolution, bundled tool management
│   ├── runner.py       # unified external tool call interface, timeout, retry
│   ├── formats.py      # format detection & conversion (FASTA/Nexus/Phylip/Phylip-PAML)
│   ├── schema.py       # shared data structures (MSACollection, TreeSet, RunRecord)
│   └── logger.py       # per-task log writer under output-dir/logs/; single-mode stderr to result.json
│
├── pretree/
│   ├── convert.py      # format conversion between FASTA/Phylip/Nexus/Phylip-PAML
│   ├── stats.py        # sequence/alignment statistics: format, length distribution,
│   │                   # taxon count, gap ratio, min/max/mean/median length, total length
│   ├── align.py        # MAFFT (external), MAGUS (pip), batch AA + NT;
│   │                   # --backtrans: AA-alignment → NT codon alignment via trimAl
│   ├── trim.py         # trimAl (bundled), BMGE (bundled), ClipKIT (pip)
│   ├── metrics.py      # MSA + tree attribute extraction, metric distributions,
│   │                   # and correlation summaries for marker evaluation
│   ├── filter.py       # gene/marker-level filtering; reads metrics output as input
│   │                   #   Error-site: TAPER; MSA-based: PIS, length, GC, nRCFV,
│   │                   #   symtest, API, likelihood mapping; Tree-based: TreeShrink,
│   │                   #   ABS, treeness, DVMC, evo_rate, saturation, inconsistent genes;
│   │                   #   UMAP cluster/outlier filtering
│   └── concat.py       # matrix generation at multiple occupancy levels,
│                       # recoding (Dayhoff6 etc.), format conversion for downstream
│
├── tree/
│   ├── ml.py           # ML tree inference: IQ-TREE and FastTree
│   ├── bi.py           # Bayesian tree inference: PhyloBayes
│   ├── msc.py          # multispecies coalescent inference: wASTRAL
│   └── cf.py           # concordance factors: gCF, sCF, sCFl, qCF
│
├── posttree/
│   ├── topology.py     # AU / WKH / WSH tree topology tests
│   ├── dating.py       # MCMCTree (PAML): fossil calibration, convergence diagnostics
│   ├── signal.py       # phylogenetic signal distribution, Four-cluster Likelihood Mapping (FcLM)
│   ├── syserror.py     # systematic error diagnosis (iterative, atom operations)
│   └── simulate.py     # AliSim (MSA simulation), gene-jackknife resampling
│
├── report/
│   ├── collector.py    # directory scanning, result.json discovery, step ordering
│   ├── templates.py    # per-command methods text generation (Python functions)
│   ├── schema.py       # ReportRecord dataclass; report.json structure
│   ├── renderer.py     # report.json → report.html via Jinja2
│   └── html/
│       └── report.html.j2  # Jinja2 HTML template
│
└── cli/
    ├── main.py         # `phyloai` entry point
    ├── doctor.py       # environment detection command
    └── commands/
        ├── pretree.py
        ├── tree.py
        ├── posttree.py
        ├── report.py
        └── run.py      # one-click pipeline
```

---

## 4. CLI Design

### 4.1 Command Structure

```bash
# Environment check (separate spec: docs/superpowers/specs/2026-06-18-phyloai-doctor-design.md)
phyloai doctor

# Pre-tree
phyloai pretree convert  --input ./raw --output-dir ./runs/pretree/convert --to fasta
phyloai pretree stats    --seq-dir ./runs/pretree/convert
phyloai pretree align    --seq-dir ./raw --method magus [--tool-args "--maxsubsetsize 50"]
phyloai pretree trim     --msa-dir ./aligned --tool bmge [--bmge-matrix BLOSUM90] \
                         [--tool-args "-g 0.5"]
phyloai pretree metrics  --msa-dir ./trimmed [--tree-dir ./genetrees]
phyloai pretree filter taper     --msa-dir ./trimmed ...
phyloai pretree filter treeshrink --tree-dir ./genetrees ...
phyloai pretree filter metrics   --table ./metrics/metrics.csv --keep "..."
phyloai pretree filter symtest  --msa-dir ./trimmed [--tree-dir ./genetrees]
phyloai pretree filter cluster   --table ./metrics/metrics.csv ...
phyloai pretree concat   --msa-dir ./filtered --taxa-occupancy 0.75

# Tree
phyloai tree ml iqtree    --matrix ./concat/matrix.fa --model C20 --guide-tree guide.nwk
phyloai tree ml fasttree  --matrix ./concat/matrix.fa
phyloai tree bi pb --matrix ./concat/matrix.phy --chains 3 --threads 8
phyloai tree msc --tree-dir ./genetrees/
phyloai tree cf --cf gcf --ref-tree species.nwk --tree-dir ./genetrees/

# Post-tree
phyloai posttree topology    --matrix ./matrix.fa --candidate-trees candidate.trees
phyloai posttree dating hessian --matrix ./matrix.fa --rooted-tree calib.tre
phyloai posttree dating mcmc  --hessian-dir ./hessian
phyloai posttree signal      --matrix ./matrix.fa --hypotheses h1.nwk,h2.nwk
phyloai posttree simulate    --tree ./tree.nwk --replicates 100 --tool alisim
phyloai posttree syserror brlen  --tree ./tree.nwk
phyloai posttree syserror cca    --matrix ./matrix.fa --t1 lg.nwk --t2 pmsf.nwk
phyloai posttree syserror sites  --matrix ./matrix.fa --tree ./tree.nwk

# Report (see 2026-06-26-phyloai-report-design.md)
phyloai report --run-dir ./runs/run/faa      # single pipeline run
phyloai report --run-dir ./runs/pretree      # single module run

# One-click pipeline
phyloai run --seq-dir ./raw --output-dir ./runs/run --mode supermatrix
phyloai run --seq-dir ./raw --output-dir ./runs/run --mode supertree
```

### 4.2 One-click Pipeline (`phyloai run`)

> **Note:** Detailed specification superseded by `2026-06-26-phyloai-run-design.md`. This section provides a summary only.

Two modes, both start with: `convert → align → trim → [filter taper]`, then diverge.

`--speed normal` (default) uses MAFFT linsi + TAPER filter + IQ-TREE3 (unpartitioned) / FastTree (gene trees). `--speed fast` uses MAFFT auto, skips TAPER, and uses FastTree for all tree inference.

| Mode | Steps | Notes |
|------|-------|-------|
| `--mode supermatrix` | convert → align → trim → [filter taper] → concat → iqtree (unpartitioned) | `--speed fast`: FastTree instead of IQ-TREE3 |
| `--mode supertree` | convert → align → trim → [filter taper] → gene trees → wastral | `--speed fast`: FastTree fast mode for gene trees |

The filter step in `phyloai run` uses `phyloai pretree filter taper` (TAPER error-site masking only). Full marker filtering is an explicit manual step via the other `phyloai pretree filter` subcommands.

### 4.3 Universal CLI Flags

| Flag | Purpose |
|------|---------|
| `--output-format json\|text` | `doctor` only; human-readable or JSON diagnostic output |
| `--dry-run` | Show what would be executed without running |
| `--threads` / `-t` | Parallelism control |
| `--run-dir` | Override default run directory |
| `--quiet` / `-q` | Suppress all output except errors |
| `--overwrite` | Overwrite existing output directory |
| `--resume` | Resume long-running commands from `checkpoint.json` |

Commands that read alignment files support `--input-format` (auto-detect by default). Commands that invoke external tools support `--tool-args` for tool strategy parameters only. PhyloAI always manages command-level input paths, output directories, work directories, and structured result/log collection. Tool strategy parameters, including data type and thread flags, should be exposed as high-level PhyloAI options when useful but may be overridden by `--tool-args` unless a subcommand spec explicitly blocks them for safety.

**Format handling policy:** Each module uses the format required by its underlying tool. `pretree convert` and `core/formats.py` provide conversion as needed. There is no global FASTA-only mandate.

Full parameter naming rules and exit code definitions are in **Section 9**.

### 4.4 Shell Completion

Provides first-party shell completion (Bash/Zsh/Fish) via `phyloai completion <shell>`, wrapping Click's built-in support. Users generate a static script once from the phyloai environment and source it from shell config. See `docs/superpowers/plans/2026-06-10-cli-completion.md` for details.

### 4.5 Documentation Layout

Top-level `README.md` is a concise entry point. Detailed command documentation lives under `docs/commands/`, one file per command (e.g., `docs/commands/doctor.md`). Each command document covers: Purpose, Usage, Inputs, Outputs, Examples, Warnings/Errors, Notes. Detailed subcommand designs and implementation plans remain under `docs/superpowers/specs/` and `docs/superpowers/plans/`.

---

## 5. Dependency Management

External tool detection, bundling strategy, and the `phyloai doctor` command are specified in the separate doctor design doc: `docs/superpowers/specs/2026-06-18-phyloai-doctor-design.md`.

Summary:
- **Bundled:** TAPER 1.0.0, BMGE 1.12 (jar); IQ-TREE3, trimAl (planned)
- **pip-installable:** MAGUS, ClipKIT
- **User-installed:** PhyloBayes-MPI, wASTRAL, MAFFT, MCMCTree (PAML), TreeShrink
- **Runtime:** Java, Julia

---

## 6. Output Directory Convention

```
runs/
├── pretree/
│   ├── 1-convert/                 # utility: no log, no checkpoint
│   │   ├── result.json
│   │   └── seqs/
│   ├── 2-align/                   # batch: per-locus logs under logs/
│   │   ├── result.json
│   │   ├── checkpoint.json
│   │   ├── seqs/
│   │   └── logs/
│   │       └── <locus>.log
│   ├── 3-stats/                   # utility: no log, no checkpoint
│   │   └── result.json
│   ├── 4-trim/                    # batch: per-locus logs under logs/
│   │   ├── result.json
│   │   ├── checkpoint.json
│   │   ├── seqs/
│   │   └── logs/
│   │       └── <locus>.log
│   ├── 5-metrics/                 # batch: no per-locus logs (pure Python)
│   │   ├── result.json
│   │   ├── metrics.csv
│   │   └── logs/
│   │       └── <marker>.log
│   ├── 6-filter/                  # subcommands vary; batch write logs/
│   │   ├── taper/
│   │   │   ├── result.json
│   │   │   ├── seqs/
│   │   │   └── logs/
│   │   ├── treeshrink/
│   │   ├── metrics/
│   │   ├── symtest/
│   │   └── cluster/
│   └── concat/                    # single: stderr in result.json
│       ├── result.json
│       ├── matrix.fa
│       └── matrix.partitions
├── tree/
│   ├── ml/
│   │   ├── fasttree/              # batch: per-locus logs under gene_trees/logs/
│   │   │   ├── result.json
│   │   │   ├── gene_trees/
│   │   │   │   ├── trees/
│   │   │   │   └── logs/
│   │   │   └── ml/                # single (--matrix): stderr in result.json
│   │   └── iqtree/                # batch or single; logs/ for batch
│   ├── bi/                        # single: tool_stderr in result.json
│   ├── msc/                       # single: tool_stderr in result.json
│   └── cf/                        # single: tool_stderr in result.json
├── posttree/
│   ├── topology/
│   ├── dating/
│   ├── signal/
│   ├── syserror/
│   └── simulate/
└── <run-dir>/report/               # written alongside the run being reported
    ├── report.json                 # machine-readable source of truth; AI/MCP entry point
    └── report.html                 # human-readable; embedded PDF figures, sortable tables
```

### 6.1 `result.json` — the single source of truth

Every non-`doctor` command writes exactly one `result.json` at its output directory root. The `result.json` is the only file the report module and MCP server read. There is no separate top-level `<step>.log` file.

### 6.2 Tool stderr handling

How tool stderr is preserved depends on the command's execution mode:

**Batch mode** (align, trim, filter, fasttree, iqtree — commands that invoke external tools per task): Each batch task writes its tool stderr to `<output-dir>/logs/<locus>.log`. The `result.json` references these files via `data.files[].log_file` but does not inline stderr. This keeps `result.json` compact and lets users `grep` individual locus logs when debugging. Pure-Python batch commands (`metrics`, `filter metrics`, `filter cluster`) invoke no external tools and omit `cmd`/`log_file` per the JSON Output Standard.

**Single mode** (concat, msc, cf): The tool's full stderr is inlined in `result.json` as `data.tool_stderr`. There is no external log file. The single-invocation stderr volume is bounded and benefits from being directly queryable in JSON.

**Utility commands** (convert, stats): No log files are written. These commands either do not invoke external tools or are read-only.

### 6.3 Checkpoint files

Batch pipeline commands (`align`, `trim`, `fasttree`, `iqtree`) write `checkpoint.json` alongside `result.json` for resume support. Single-mode and utility commands do not use checkpoints.

---

## 7. Report Module

> **Full specification:** `docs/superpowers/specs/2026-06-26-phyloai-report-design.md`

`phyloai report --run-dir <path>` is a single command that produces two output files:

- **`report.json`** — machine-readable source of truth; structured entry point for AI/MCP diagnostics, reproducibility audits, and archival. Contains all step records, parameters, key results, methods text per step, and indices of all figures and tables.
- **`report.html`** — human-readable report; embeds PDF figures natively, renders sortable/collapsible tables, and includes a copyable Methods paragraph draft suitable for journal submission.

`report.html` is fully derived from `report.json` and can be re-rendered at any time without re-scanning the run directory.

**Directory auto-detection:** `report` identifies two run structures automatically:
- `pipeline` — `phyloai run` two-layer output (top-level `result.json` + per-step subdirectories)
- `module` — single-module output (one or more step subdirectories, with or without a top-level `result.json`)

**Methods paragraph:** Generated from per-command Python template functions in `templates.py`. Each template describes the tool used, key scientific parameters and their meaning, inputs, and quantitative outcomes. Templates are deterministic (no LLM involvement) and cover all scientific parameters; technical parameters (threads, paths) are omitted.

**Figures and tables:** Every command records all persistent output files it produces under `data.output_files` in its `result.json` (see JSON Output Standard Section 5.4). Each entry is an object with a required `"path"` and an optional `"description"` describing the file's content and analytical role. `phyloai report` reads `data.output_files` from each step to build `report.json:figures_index` (`.pdf`/`.png`) and `tables_index` (`.csv`/`.tsv`) — no hardcoded paths are needed. PDF figures are embedded directly in `report.html` via `<object>` tags, preserving vector quality.

**Incomplete runs:** Report generation always succeeds regardless of step failures. Failed steps are included in `report.json` with full error details; their `methods_text` is empty and they are excluded from the methods paragraph. `report.html` marks failed steps visually and expands their detail cards by default.

All steps including `pretree convert`, `pretree stats`, and `pretree metrics` contribute `methods_text` to the report.

---

## 8. AI Integration (Post-CLI Phase)

> **Full specification:** `docs/superpowers/specs/2026-06-27-phyloai-ai-integration-design.md`

Two components built on top of the stable CLI:

```
用户 ←→ Skill (对话/决策/解读) ←→ MCP Server (执行桥) ←→ phyloai CLI ←→ 文件系统
```

**MCP Server (Phase 7):** One tool per CLI subcommand. Schemas generated dynamically from the Click command tree at startup — zero manual sync. All commands fire-and-forget; `output_dir` is the persistent job handle across sessions. Transport: stdio.

**Skill `phyloai-workflow` (Phase 8):** Guided workflow with parameter confirmation cards, result interpretation, session recovery via `report.json`, demo mode, and error handling. Lives in `skills/phyloai-workflow/` inside this repo, version-coupled to CLI. Future: `phyloai-syserror` sub-workflow Skill and AI-assisted report review (separate specs).

---

## 9. CLI Conventions and Parameter Standards

This section defines binding conventions for all subcommand implementations. Per-subcommand designs must follow these rules. Deviations require explicit justification in the subcommand spec.

### 9.1 Naming Style

- All parameter names use **kebab-case**: `--msa-dir`, `--seq-type`, `--output-dir`
- Directories always use the `--xxx-dir` suffix; single files omit the suffix: `--tree`, `--matrix`, `--calibrations`
- Short aliases (`-t`, `-o`) are reserved only for the highest-frequency parameters listed in Section 9.2; all other parameters use long form only

### 9.2 Shared Parameter Registry

| Parameter | Short | Type | Default | Applicable commands |
|-----------|-------|------|---------|---------------------|
| `--msa-dir` | | Path | — | all commands reading an alignment directory |
| `--tree-dir` | | Path | — | commands requiring a gene tree directory |
| `--output-dir` | `-o` | Path | auto under `runs/` | all commands producing output files |
| `--threads` | `-t` | int | 4 | all commands invoking multi-threaded tools |
| `--seq-type` | | AA \| NT \| CODON \| auto | auto where safe | commands where molecule type matters |
| `--speed` | | `normal\|fast` | `normal` | `run` only |
| `--tool` | | str | tool-specific | commands offering multiple tool choices |
| `--input-format` | | str | auto-detect | commands reading alignment or tabular files |
| `--output-format` | | text \| json | text | `doctor` only |
| `--table-format` | | `csv\|tsv` | `csv` | commands writing auxiliary tabular outputs |
| `--tool-args` | | str | — | all commands invoking external tools; strategy parameters only |
| `--run-dir` | | Path | `./runs/` | all commands except `stats` and `convert` (root-level) |
| `--dry-run` | | flag | False | all commands |
| `--quiet` | `-q` | flag | False | all commands |
| `--overwrite` | | flag | False | all commands producing output directories |
| `--resume` | | flag | False | long-running pipeline commands only |

`--run-dir` is defined once on the `phyloai` root group and passed via Click context. `--input-format` keeps the same flag name across alignment-reading and table-reading commands, but the value domain is command-specific.

### 9.3 Exit Codes

| Code | Meaning |
|------|---------|
| `0` | Success |
| `1` | User input error (bad parameter, file not found, invalid value) |
| `2` | External tool execution failed (non-zero returncode from subprocess) |
| `3` | Environment error (required tool not installed or not detectable) |

### 9.4 JSON Output Schema

Every command writes a JSON result file to its output directory. All fields up to and including `error` are required and genre-neutral. `data` is command-specific.

```json
{
  "status": "success | error",
  "command": "phyloai pretree align --seq-dir ./raw --method linsi --threads 8",
  "wall_time": 142.3,
  "tool_versions": {"mafft": "7.520"},
  "params": { },
  "key_results": { },
  "error": null,
  "data": { }
}
```

**Field definitions:**

| Field | Required | Description |
|-------|----------|-------------|
| `status` | yes | `"success"` or `"error"` |
| `command` | yes | **Full CLI invocation string** with all resolved arguments (see 9.4.1) |
| `wall_time` | yes | Total wall-clock seconds (float) |
| `tool_versions` | yes | `{tool_name: version}` for every external tool invoked; `{}` if none |
| `params` | yes | All resolved input parameters (see 9.4.2) |
| `key_results` | yes | Quantitative outputs for report integration; `{}` for utility commands |
| `error` | yes | `null` on success; error message string on failure |
| `data` | yes | Command-specific detailed results (see `2026-06-21-phyloai-json-output-standard.md`) |

All commands follow the same output convention:
1. `result.json` always at output directory root
2. Data files alongside `result.json`, or in `seqs/` subdirectory for sequence outputs
3. Auxiliary files alongside `result.json`

Detailed field semantics, batch/single structural patterns, per-module requirements, and testing assertions are in the **JSON Output Standard**: `docs/superpowers/specs/2026-06-21-phyloai-json-output-standard.md`. Key rules in brief:

- `command`: full CLI invocation string with all resolved arguments.
- `params`: all parameters in resolved form. Order unspecified; key by name.
- `data`: follows either the **batch pattern** (`data.files[]` with per-task `cmd` and `log_file`) or the **single pattern** (`data.cmd`, `data.tool_stderr` inlined). Batch external-tool stderr is written to `logs/<locus>.log`, never inlined. Single-mode stderr is inlined; no external log file.
- Tool stderr stores raw text only; summaries belong in `warnings` or `data.summary`.

### 9.5 Output Directory Conflict and Resume Policy

Directory-producing commands use this policy:

- Default: if output directory exists and is non-empty, exit with code 1
- `--overwrite`: delete and recreate the output directory
- `--resume`: for long-running pipeline commands, load `checkpoint.json`, require exact parameter match, skip verified successful tasks
- `--overwrite` and `--resume` are mutually exclusive

Artifact-producing commands (for example `report` and metrics plotting commands) apply the same conflict check to each declared output artifact. Their `--overwrite` option replaces only those declared artifacts and does not delete unrelated files in the output directory.

Short utility commands (`pretree convert`, `pretree stats`) do not need resume support. A dry run never deletes an existing output directory, even when `--overwrite` is present.

Detailed checkpoint schema: `docs/superpowers/specs/2026-06-12-checkpoint-resume-design.md`.

### 9.6 Display and Logging

- Terminal output uses **Rich**: progress bars for batch operations, tables for summary results, colored status indicators
- Resume-aware progress bars display remaining work, not historical work
- Non-`doctor` commands always write structured output to `result.json`; no text/json stdout switch
- `--quiet` suppresses all terminal output except errors
- Every CLI command provides high-readability `--help` text
- When a command writes output files, terminal output states what was saved and where

**Tool stderr model** (see also Section 6.2):

- **Batch commands** (align, trim, filter, fasttree, iqtree per-gene) that invoke external tools: each task's tool stderr is written to `<output-dir>/logs/<locus>.log`. The `result.json` references these via `data.files[].log_file`. Per-task stderr is NOT duplicated in `result.json`. Pure-Python batch commands (`metrics`, `filter metrics`, `filter cluster`) invoke no external tools per task; they omit both `files[].cmd` and `files[].log_file`.
- **Single commands** (concat, msc, cf): tool stderr is inlined in `result.json` as `data.tool_stderr`. No external log file is written.
- **Utility commands** (convert, stats): no log files are written.
- On `--resume`: per-task log files are appended with a `=== RESUME ... ===` separator. On `--overwrite`: the `logs/` directory is deleted and recreated.

### 9.7 Shared File Matching Policy

Commands pairing flat MSA and tree directories match files by **logical locus name** (filename before the final `.`), not by fixed suffix whitelist. For tree inputs, suffix-agnostic matching attempts removing 1–2 dot segments. Ambiguity raises an explicit error. This policy does not depend on a hard-coded list of recognized suffixes.

### 9.8 Tabular Input and Output Policy

Commands reading CSV/TSV tables expose `--input-format csv|tsv|auto` (default `auto`). Commands writing auxiliary tabular outputs expose `--table-format csv|tsv` (default `csv`). `auto` is input-only. Auto-detection prefers delimiter inspection; if uncertain, exit with an error rather than guessing.

### 9.9 `--tool-args` Strategy-Only Semantics

PhyloAI uses a two-tier model for `--tool-args` interaction:

**Tier 1 — BLOCKED flags (hard-rejected):**

PhyloAI always controls command-level input file paths that define the PhyloAI operation and must not allow those paths to be redirected through `--tool-args`. Blocked flags are defined per subcommand and must be minimal. Common blocked items:
- The tool input flags corresponding to required PhyloAI input parameters (e.g., `-s` for IQ-TREE when `--matrix` is the PhyloAI input)
- Shell I/O redirect tokens (`>`, `<`, `|`, etc.)

If `--tool-args` contains any blocked flag, exit code 1 with the blocked flag name.

**Tier 2 — OVERRIDEABLE parameters (suppress-if-present):**

PhyloAI generates tool flags for its own managed parameters (model, bootstrap, partitions, etc.). If the user provides the same flag via `--tool-args`, PhyloAI suppresses its own version and lets `--tool-args` win. This allows users to bypass PhyloAI's parameter interface for edge cases while keeping the structured API for common use.

**Assembly order:**

1. PhyloAI assembles its base command (input file, output control, threads)
2. For each managed parameter, PhyloAI checks whether `--tool-args` already contains that flag; if so, skip PhyloAI's version; otherwise append it
3. `--tool-args` is tokenized with `shlex.split` and remaining tokens appended
4. No per-tool parser reimplementation; flag-name overlap check only

### 9.10 Generated Sequence and Alignment Validation

- Commands that generate sequence or alignment files in bulk (`pretree convert`, `pretree align`, `pretree trim`, and related future steps) must validate generated files before counting them as successful outputs.
- Validation is implemented through shared `core` helpers, not one-off command-local checks. FASTA validation must detect missing/empty files, unparsable FASTA, zero FASTA records, and empty sequences. MSA validation must additionally detect unequal sequence lengths.
- Per-file validation failures are recorded in `result.json` under `data.skipped` or per-file warnings, depending on whether the command can continue. If all generated files fail validation, the command exits with code 1.
- Input scanning follows the same principle: empty files, directories, unrecognized file types, and unparsable sequence files are skipped with explicit reasons; if no valid inputs remain, the command exits with code 1.

### 9.11 FASTA Line Wrapping Policy

All PhyloAI-authored FASTA-family outputs must wrap sequence lines at 60 characters. Applies to `pretree convert`, `pretree align`, `pretree trim`, `pretree concat`, and `pretree filter` (all subcommands) outputs. PHYLIP and NEXUS outputs keep their format-specific serialization rules. Externally-tool-generated output files (e.g., raw TAPER stdout) are reformatted to 60-char wrapping by PhyloAI before saving.

---

## 10. Platform Support

- **Primary:** Linux, macOS
- **Secondary:** WSL (Windows Subsystem for Linux)
- **Not supported:** Native Windows
- Tool-specific exceptions allowed when upstream binaries are not cross-platform (e.g., `pretree align --method magus` is Linux-only in Phase 2)

---

## 11. Development Phases

| Phase | Scope | Deliverable | Pre-requisites |
|-------|-------|-------------|----------------|
| 1 | `core/` infrastructure | env, runner, formats, logger | — |
| 2 | `pretree/` modules | stats, convert, align, trim, metrics, filter, concat | Phase 1 |
| 3 | `tree/` modules | ml, bi, msc, cf | Phase 1 |
| 4 | `posttree/` modules | topology, dating, signal, syserror, simulate | Phases 2–3 |
| 5 | `phyloai run` | one-click supermatrix and supertree pipelines | Phases 2–3 |
| 6 | `report/` module | collector, templates, schema, renderer; outputs report.json + report.html | Phases 2–4 |
| 7 | MCP Server | All CLI tools wrapped (fine-grained, one tool per subcommand); check_status / read_result / read_report / get_command_schema utilities; stub tools for future commands; stdio transport | Phases 1–6 |
| 8 | `phyloai-workflow` Skill | Full guided workflow; parameter cards with runtime schema; result interpretation; session recovery via report.json; demo mode; error handling (catalog + AI) | Phase 7 |
| 9 | `phyloai-syserror` Skill | Results-driven syserror orchestration sub-workflow (brlen → cca → sites) | syserror CLI (Phase 4) + Phase 7 |
| 10 | Report AI review | `polish_methods` MCP tool + Skill integration; scientific accuracy verification of methods text | Separate spec required |

**Spec granularity:** Phase 1 has one spec+plan. Phases 2–4 have one spec+plan per subcommand under `docs/superpowers/specs/` and `docs/superpowers/plans/`. Phases 5–8 have one spec+plan each. Every subcommand spec must be consistent with Section 9 conventions.

**Completion requirements for every subcommand:**
1. `result.json` output for all non-`doctor` commands
2. Output directory follows Section 6
3. Conflict policy follows Section 9.5
4. Tool stderr handled per Section 6.2 and 9.6 (batch: `logs/<locus>.log`; single: `data.tool_stderr` inlined; utility: no logs)
5. JSON schema follows Section 9.4
6. Exit codes follow Section 9.3

Testing for `result.json` structural compliance follows the assertions in `docs/superpowers/specs/2026-06-21-phyloai-json-output-standard.md` Section 8.

Modules within each phase can be developed in parallel. Phases are strictly sequential.

---

## 12. Key Design Decisions

| Decision | Rationale |
|----------|-----------|
| No GUI | CLI + AI interaction covers all use cases |
| stats/convert as pretree subcommands | Pre-analysis utilities, keeps top-level CLI clean |
| filter uses subcommands | TAPER, TreeShrink, metric-rule, and cluster filtering have different inputs and tool dependencies |
| metrics is prerequisite for filter | Decouples computation from filtering decisions |
| Filtering split by evidence source | Error-site, taxon-pruning, rule-based, and cluster-based each separate |
| --tool-args strategy-only model | Deterministic batch I/O while exposing tool-specific knobs |
| Format handling per-module | Different tools need different formats; per-module via core/formats.py |
| backtrans in align | Direct post-processing of alignment, uses trimAl already a dependency |
| syserror exposed as atomic ops only | Full diagnosis needs iterative human decisions; CLI atomics + Skill orchestration |
| genetree in tree/ not pretree/ | Gene trees are tree inference results, not preprocessing steps |
| JSON result.json for non-doctor commands | One stable machine-readable result path for MCP wrapping |
| JSON key_results in all pipeline modules | Enables report summary and methods text without post-hoc data extraction |
| report.json as session recovery entry point | Aggregates all step records; Skill calls read_report at session start to reconstruct run state without user explanation |
| report.json not auto-generated after each step | report generation scans all steps; too expensive for routine post-command use; per-step result.json used for immediate interpretation |
| report single command, not sub-commands | collector/templates/renderer are internal; user interface is one invocation |
| report.html embeds PDF figures natively | Preserves vector quality without format conversion; requires no extra dependencies |
| AI integration design decisions | See `docs/superpowers/specs/2026-06-27-phyloai-ai-integration-design.md` |
| Logical locus matching ignores suffix vocabularies | Real datasets use inconsistent suffixes |
| Table I/O uniform flags | Users don't relearn different table flags per command |
| Per-subcommand design + plan docs | Each subcommand complex enough to warrant own spec; main doc stays stable |
| Lightweight README + command docs | README maintainable; detailed behavior in `docs/commands/*.md` |
| Doctor in separate spec | Doctor is a standalone diagnostic concerned with tool detection, bundling, and registry; separate spec reduces main doc weight |
