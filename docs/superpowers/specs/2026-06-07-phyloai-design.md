# PhyloAI Design Specification

**Date:** 2026-06-07  
**Last updated:** 2026-06-19 (section 9.9 tool-args semantics clarified to two-tier model)  
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
│   └── logger.py       # per-step log file writer for command output dirs
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
│   └── concordance.py  # concordance factors: gCF, sCF, combined summaries
│
├── posttree/
│   ├── topology.py     # AU / WKH / WSH tests, Four-cluster Likelihood Mapping (FcLM)
│   ├── dating.py       # MCMCTree (PAML): fossil calibration, convergence diagnostics
│   ├── signal.py       # phylogenetic signal distribution
│   ├── syserror.py     # systematic error diagnosis (iterative, atom operations)
│   └── simulate.py     # AliSim (MSA simulation), gene-jackknife resampling
│
├── report/
│   ├── collector.py    # cross-module JSON aggregation
│   ├── methods.py      # Methods paragraph generation from run record
│   └── summary.py      # statistical summary, run_record.yaml output
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
phyloai tree bi phylobayes --matrix ./concat/matrix.phy --chains 3 --threads 8
phyloai tree msc wastral  --gene-trees ./genetrees/

# Post-tree
phyloai posttree topology    --matrix ./matrix.fa --hypotheses h1.nwk,h2.nwk,h3.nwk
phyloai posttree dating      --tree ./tree.nwk --matrix ./matrix.fa \
                             --calibrations calibrations.txt
phyloai posttree signal      --matrix ./matrix.fa --hypotheses h1.nwk,h2.nwk
phyloai posttree simulate    --tree ./tree.nwk --replicates 100 --tool alisim
phyloai posttree syserror brlen  --tree ./tree.nwk
phyloai posttree syserror cca    --matrix ./matrix.fa --t1 lg.nwk --t2 pmsf.nwk
phyloai posttree syserror sites  --matrix ./matrix.fa --tree ./tree.nwk

# Report
phyloai report generate --run-dir ./runs

# One-click pipeline
phyloai run --msa-dir ./raw --output ./runs --mode supermatrix
phyloai run --msa-dir ./raw --output ./runs --mode coalescent
```

### 4.2 One-click Pipeline (`phyloai run`)

Two modes, both include: `align → trim → filter taper → concat → [tree inference]`

| Mode | Steps | Notes |
|------|-------|-------|
| `--mode supermatrix` | align → trim → filter taper → concat → iqtree (unpartitioned) | Fast, single-step ML tree |
| `--mode coalescent` | align → trim → filter taper → concat → genetree → wastral | MSC-based species tree |

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

Commands that read alignment files support `--input-format` (auto-detect by default). Commands that invoke external tools support `--tool-args` for tool strategy parameters only. PhyloAI always manages input, output, work directory, data type, threads, logs, and codon/projection arguments.

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
│   ├── align/
│   │   ├── seqs/
│   │   ├── align.log
│   │   └── result.json
│   ├── trim/
│   │   ├── seqs/
│   │   ├── trim.log
│   │   └── result.json
│   ├── metrics/
│   │   ├── metrics.log
│   │   └── result.json
│   ├── filter/
│   │   ├── taper/
│   │   │   ├── seqs/
│   │   │   ├── filter.log
│   │   │   └── result.json
│   │   ├── treeshrink/
│   │   ├── metrics/
│   │   └── cluster/
│   └── concat/
│       ├── concat.log
│       └── result.json
├── tree/
│   ├── ml/
│   │   ├── iqtree/
│   │   └── fasttree/
│   ├── bi/phylobayes/
│   ├── msc/wastral/
│   └── concordance/
├── posttree/
│   ├── topology/
│   ├── dating/
│   ├── signal/
│   ├── syserror/
│   └── simulate/
└── report/
    ├── summary.json
    ├── methods.txt
    └── run_record.yaml
```

Log file content per step: tool version, full command, stderr, wall time, exit code, and stdout only when stdout is diagnostic text. Commands must not duplicate large primary data streams in logs. Pipelines writing multiple outputs under a subdirectory (e.g., `seqs/`) place the log in the output directory root alongside `result.json`.

---

## 7. Report Module

- Each module writes `result.json` with `params` and `key_results` on completion.
- `report collector` aggregates all JSON files from a run directory tree.
- `report methods` renders a Methods paragraph from the aggregated record.
- `report summary` outputs `run_record.yaml`.
- `report figures` renders tables and plots from `key_results`.

`pretree stats` and `pretree convert` are utility commands and do not contribute to `report`.

---

## 8. AI Integration (Post-CLI Phase)

**MCP Server:** Expose all CLI commands as MCP tools with JSON I/O, reading `result.json` for structured results.

**AI Coding Assistant Skill:** Workflow orchestration on top of MCP tools — parameter recommendation, result interpretation, Methods paragraph generation trigger, dynamic next-step suggestion.

**Systematic error diagnosis** (`posttree syserror`) exposes atomic operations individually via CLI. The Skill layer orchestrates them into a guided interactive sub-workflow. This cannot be encoded as a single CLI command.

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

Every command writes a JSON result file to its output directory:

```json
{
  "status": "success | error",
  "command": "phyloai pretree align ...",
  "wall_time": 142.3,
  "tool_versions": {"mafft": "7.520"},
  "params": { },
  "key_results": { },
  "error": null,
  "data": { }
}
```

- `params`: all resolved input parameters
- `key_results`: quantitative outputs for report integration; empty `{}` for utility commands
- `tool_versions`: version string for every external tool invoked
- `error`: null on success; error message string on failure
- `data`: command-specific detailed results

All commands follow the same output convention:
1. `result.json` always at output directory root
2. Data files alongside `result.json`, or in `seqs/` subdirectory for sequence outputs
3. Auxiliary files alongside `result.json`

### 9.5 Output Directory Conflict and Resume Policy

- Default: if output directory exists and is non-empty, exit with code 1
- `--overwrite`: delete and recreate the output directory
- `--resume`: for long-running pipeline commands, load `checkpoint.json`, require exact parameter match, skip verified successful tasks
- `--overwrite` and `--resume` are mutually exclusive
- Short utility commands (`pretree convert`, `pretree stats`) do not need resume support

Detailed checkpoint schema: `docs/superpowers/specs/2026-06-12-checkpoint-resume-design.md`.

### 9.6 Display and Logging

- Terminal output uses **Rich**: progress bars for batch operations, tables for summary results, colored status indicators
- Resume-aware progress bars display remaining work, not historical work
- Non-`doctor` commands always write structured output to `result.json`; no text/json stdout switch
- `--quiet` suppresses all terminal output except errors
- Pipeline commands write `<step>.log` to their output directory alongside `result.json`. Logs are appended on retry/resume. On `--overwrite`, the log is recreated
- Every CLI command provides high-readability `--help` text
- When a command writes output files, terminal output states what was saved and where

### 9.7 Shared File Matching Policy

Commands pairing flat MSA and tree directories match files by **logical locus name** (filename before the final `.`), not by fixed suffix whitelist. For tree inputs, suffix-agnostic matching attempts removing 1–2 dot segments. Ambiguity raises an explicit error. This policy does not depend on a hard-coded list of recognized suffixes.

### 9.8 Tabular Input and Output Policy

Commands reading CSV/TSV tables expose `--input-format csv|tsv|auto` (default `auto`). Commands writing auxiliary tabular outputs expose `--table-format csv|tsv` (default `csv`). `auto` is input-only. Auto-detection prefers delimiter inspection; if uncertain, exit with an error rather than guessing.

### 9.9 `--tool-args` Strategy-Only Semantics

PhyloAI uses a two-tier model for `--tool-args` interaction:

**Tier 1 — BLOCKED flags (hard-rejected):**

PhyloAI always controls the tool's input file path and must not allow I/O redirection through `--tool-args`. Blocked flags are defined per subcommand and must be minimal. Common blocked items:
- The tool's input-file flag (e.g., `-s` for IQ-TREE)
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

All PhyloAI-authored FASTA-family outputs must wrap sequence lines at 60 characters. Applies to `pretree convert`, `pretree align`, `pretree trim`, and `pretree concat` outputs. PHYLIP and NEXUS outputs keep their format-specific serialization rules.

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
| 3 | `tree/` modules | ml, bi, msc, concordance | Phase 1 |
| 4 | `posttree/` modules | topology, dating, signal, syserror, simulate | Phases 2–3 |
| 5 | `phyloai run` | one-click supermatrix and coalescent pipelines | Phases 2–3 |
| 6 | `report/` module | collector, methods, figures, run_record.yaml | Phases 2–4 |
| 7 | MCP Server | JSON tool interface | Phases 2–6 |
| 8 | AI Coding Assistant Skill | workflow orchestration, syserror sub-workflow | Phase 7 |

**Spec granularity:** Phase 1 has one spec+plan. Phases 2–4 have one spec+plan per subcommand under `docs/superpowers/specs/` and `docs/superpowers/plans/`. Phases 5–8 have one spec+plan each. Every subcommand spec must be consistent with Section 9 conventions.

**Completion requirements for every subcommand:**
1. `result.json` output for all non-`doctor` commands
2. Output directory follows Section 6
3. Conflict policy follows Section 9.5
4. Log file written as `<step>.log` for pipeline commands
5. JSON schema follows Section 9.4
6. Exit codes follow Section 9.3

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
| JSON key_results in all pipeline modules | Enables report figures/summary without post-hoc data extraction |
| Logical locus matching ignores suffix vocabularies | Real datasets use inconsistent suffixes |
| Table I/O uniform flags | Users don't relearn different table flags per command |
| Per-subcommand design + plan docs | Each subcommand complex enough to warrant own spec; main doc stays stable |
| Lightweight README + command docs | README maintainable; detailed behavior in `docs/commands/*.md` |
| Doctor in separate spec | Doctor is a standalone diagnostic concerned with tool detection, bundling, and registry; separate spec reduces main doc weight |
