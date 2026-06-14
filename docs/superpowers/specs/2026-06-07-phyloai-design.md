# PhyloAI Design Specification

**Date:** 2026-06-07  
**Last updated:** 2026-06-12 (checkpoint/resume policy added for long-running commands; pretree convert order clarified; README split into command-level docs; convert default output aligned with run layout; TAPER/BMGE bundled paths clarified)  
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
**MCP pre-requisite:** All non-`doctor` commands write a standard `result.json` file and use standard exit codes from day one, so MCP wrapping reads structured results from files rather than command-specific stdout formats. `phyloai doctor` is the only command with `--output-format text|json` because it is primarily a human-facing environment diagnostic.

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
│   ├── convert.py      # format conversion between FASTA/Phylip/Nexus/Phylip-PAML;
│   │                   # wraps core/formats.py, exposes as CLI subcommand
│   ├── stats.py        # sequence/alignment statistics: format, length distribution,
│   │                   # taxon count, gap ratio, min/max/mean/median length, total length
│   ├── align.py        # MAFFT (external), MAGUS (pip), batch AA + NT;
│   │                   # --backtrans: AA-alignment → NT codon alignment via trimAl
│   ├── trim.py         # trimAl (bundled), BMGE (bundled), ClipKIT (pip)
│   ├── metrics.py      # MSA + tree attribute extraction (16 MSA + 14 tree metrics),
│   │                   # layered: core metrics always computed; advanced (UMAP,
│   │                   # correlation matrix) via --advanced flag;
│   │                   # reference: github.com/xtmtd/MSA-and-tree-metrics-exploration
│   ├── filter.py       # gene/marker-level filtering; reads metrics output as input;
│   │                   #   MSA-based:  PIS, length, GC content, nRCFV,
│   │                   #               symtest, API, likelihood mapping
│   │                   #   Tree-based: TreeShrink, ABS, treeness, DVMC,
│   │                   #               evo_rate, saturation, inconsistent genes
│   │                   #   Error-site: TAPER (fast, non-destructive masking option)
│   └── concat.py       # matrix generation at multiple occupancy levels,
│                       # recoding (Dayhoff6 etc.), format conversion for downstream
│
├── tree/
│   ├── genetree.py     # batch gene tree inference (IQ-TREE, parallel)
│   ├── iqtree.py       # partitioned / unpartitioned / mixture (C20-C60, EX_EHO) /
│   │                   # GHOST / PMSF / recoding (Dayhoff6) / model selection
│   ├── astral.py       # astral-hybrid species-tree inference, local branch support
│   └── phylobayes.py   # CAT-GTR, convergence checking (bpcomp, tracecomp),
│                       # chain management (start, stop, resume)
│
├── posttree/
│   ├── concordance.py  # gCF, sCF, combined gCF+sCF
│   ├── topology.py     # AU / WKH / WSH tests, Four-cluster Likelihood Mapping (FcLM)
│   ├── dating.py       # MCMCTree (PAML): fossil calibration, convergence diagnostics,
│   │                   # infinite-sites plots, prior vs posterior comparison
│   ├── signal.py       # phylogenetic signal distribution:
│   │                   # site-wise likelihood support, gene-wise likelihood analysis,
│   │                   # FcLM-based signal, consistent gene detection
│   ├── syserror.py     # systematic error diagnosis (iterative, atom operations):
│   │                   #   - branch length heterogeneity (LBA detection, brlen stats)
│   │                   #   - rates across sites (fast/slow site analysis)
│   │                   #   - compositional bias (CCA, Keff analysis)
│   │                   #   - substitution pattern heterogeneity
│   │                   #   - model fit assessment (LOO-CV/PPC, internal use only)
│   │                   # note: full diagnosis workflow orchestrated in Skill layer
│   └── simulate.py     # AliSim (MSA simulation, empirical PDF parameters),
│                       # gene-jackknife resampling
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
# Environment check
phyloai doctor

# Pre-tree
phyloai pretree convert  --input ./raw --output-dir ./runs/pretree/convert --to fasta
phyloai pretree stats    --seq-dir ./runs/pretree/convert
phyloai pretree align    --seq-dir ./raw --method magus [--tool-args "--maxsubsetsize 50"]
phyloai pretree trim     --msa-dir ./aligned --tool bmge [--bmge-matrix BLOSUM90] \
                         [--tool-args "-g 0.5"]
phyloai pretree metrics  --msa-dir ./trimmed [--tree-dir ./genetrees] [--advanced]
phyloai pretree filter   --msa-dir ./trimmed --metrics-dir ./metrics \
                         [--tree-dir ./genetrees] \
                         --strategy outlier [--filter-by pis,abs,treeness] \
                         [--taper]
phyloai pretree concat   --msa-dir ./filtered --taxa-occupancy 0.75

# Tree
phyloai tree genetree    --msa-dir ./filtered --model LG --threads 4
phyloai tree iqtree      --matrix ./concat/matrix.fa --mode pmsf \
                         [--guide-tree ./tree.nwk]
phyloai tree astral      --gene-trees ./genetrees/ --mode wastral
phyloai tree phylobayes  --matrix ./concat/matrix.phy --chains 3 --threads 8

# Post-tree
phyloai posttree concordance --tree ./tree.nwk --gene-trees ./genetrees/
phyloai posttree topology    --matrix ./matrix.fa --hypotheses h1.nwk,h2.nwk,h3.nwk
phyloai posttree dating      --tree ./tree.nwk --matrix ./matrix.fa \
                             --calibrations calibrations.txt
phyloai posttree signal      --matrix ./matrix.fa --hypotheses h1.nwk,h2.nwk
phyloai posttree simulate    --tree ./tree.nwk --replicates 100 --tool alisim

# Systematic error diagnosis: individual atomic operations via CLI,
# full iterative workflow via Skill (post-CLI-stable)
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

Two modes, both include: `align → trim → filter (TAPER only) → concat → [tree inference]`

| Mode | Steps | Notes |
|------|-------|-------|
| `--mode supermatrix` | align → trim → filter (TAPER) → concat → iqtree (unpartitioned) | Fast, single-step ML tree |
| `--mode coalescent` | align → trim → filter (TAPER) → concat → genetree → astral-hybrid | MSC-based species tree |

The filter step in `phyloai run` uses TAPER only (fast error site masking). Full marker filtering — including MSA-based and tree-based criteria — is an explicit manual step via `phyloai pretree filter`.

### 4.3 Universal CLI Flags

All commands support the shared parameters defined in Section 9.2. Key universal flags:

| Flag | Purpose |
|------|---------|
| `--output-format json\|text` | `doctor` only; choose human-readable text or machine-readable JSON diagnostic output |
| `--dry-run` | Show what would be executed without running |
| `--config FILE` | Load parameters from YAML (HPC batch use); see Section 13 for template |
| `--threads` / `-t` | Parallelism control |
| `--run-dir` | Override default run directory |
| `--quiet` / `-q` | Suppress all output except errors |
| `--overwrite` | Overwrite existing output directory |
| `--resume` | Resume long-running commands from `checkpoint.json` |

Commands that read alignment files also support `--input-format` (auto-detect by default; accepted values match `AlignmentFormat` enum). Commands that invoke external tools also support `--tool-args` for tool strategy parameters only. PhyloAI always manages input, output, work directory, data type, threads, logs, and codon/projection arguments.

**Format handling policy:** Each module uses the format required by its underlying tool. `pretree convert` and `core/formats.py` provide conversion as needed. There is no global FASTA-only mandate; modules handle format requirements internally.

Full parameter naming rules and exit code definitions are in **Section 9**.

### 4.4 Shell Completion

PhyloAI should provide first-party shell completion for the full CLI surface without adding any new Python dependency beyond the existing `click>=8.1` requirement.

The supported scope is:

- command and subcommand completion for the entire `phyloai` CLI tree
- option-name completion for all commands
- fixed enumerated parameter values exposed through `click.Choice(...)`
- path completion delegated to Click and the user's shell integration

PhyloAI should not implement custom runtime-dependent completion values at this stage. In particular, completion should not inspect project configuration, plugin state, or filesystem-derived domain values beyond what Click and the shell already provide automatically.

The public interface should be a top-level command group:

```bash
phyloai completion bash
phyloai completion zsh
phyloai completion fish
```

Each subcommand should print the shell-specific completion script to standard output and exit successfully without modifying the user's shell configuration files. This command is intentionally a thin, stable wrapper around Click's built-in shell completion support. The wrapper exists to hide Click's internal environment-variable interface from end users and to give PhyloAI a stable, documented command surface.

The recommended installation model is static script generation, not shell-startup-time dynamic evaluation. Users should run the completion export command once from an environment where `phyloai` is installed, save the generated script to a persistent location, and source that static file from their shell configuration. This avoids startup failures when the user's default shell environment is not the Conda environment that contains `phyloai` and `click`.

The documentation must explicitly avoid recommending unconditional shell startup snippets that execute `phyloai` on every new shell session, such as dynamic `eval` patterns. Those patterns are fragile in Conda-based workflows because the default shell may not have access to the `phyloai` executable or its Python dependencies. Static script loading is the default documented path for Bash, Zsh, and Fish.

Example usage pattern:

```bash
# Generate once from the phyloai environment
phyloai completion zsh > ~/.config/phyloai/completion/phyloai.zsh

# Then source the static script from shell configuration
source ~/.config/phyloai/completion/phyloai.zsh
```

This feature is CLI-only. It does not require changes to the library API, MCP design, or run-record schema.

### 4.5 Documentation Layout

The top-level `README.md` should remain a concise project entry point, not a full command manual. It should contain only the project summary, installation, quick start, command index, shell completion summary, and links to command-specific documentation.

Detailed command documentation lives under `docs/commands/`, one file per user-facing command or subcommand. Examples:

- `docs/commands/doctor.md`
- `docs/commands/pretree-stats.md`
- `docs/commands/pretree-convert.md`

Each command document must include these sections:

- Purpose: what the command does and explicitly what it does not do
- Usage: minimal usage and full parameter table
- Inputs: supported formats, input modes, and directory scanning rules where applicable
- Outputs: terminal output, saved files, JSON schema summary, and default output paths
- Examples: common single-file and batch workflows where applicable
- Warnings and Errors: skipped inputs, overwrite behavior, format detection failures, and command-specific warnings
- Notes: relationship to adjacent commands, such as `pretree convert` before `pretree stats`

New or materially changed commands must update both their command document and the top-level README command index. Detailed subcommand designs and implementation plans remain separate under `docs/superpowers/specs/` and `docs/superpowers/plans/`.

---

## 5. Dependency Management

### 5.1 Strategy

- **Bundled in the Python package:** TAPER 1.0.0 (`phyloai/bundled/TAPER-1.0.0/correction_multi.jl`), BMGE 1.12 (`phyloai/bundled/BMGE-1.12/BMGE.jar`) with corresponding license/source files retained in the bundled directories
- **Bundled or auto-downloaded later:** trimAl — exact packaging strategy to be finalized before `pretree trim`
- **pip-installable (auto-installed):** MAGUS (`pip install magus-msa`), ClipKIT, PhyKIT
- **User-installed (detected by `doctor`):** IQ-TREE3, PhyloBayes-MPI, astral-hybrid, MAFFT, MCMCTree (PAML), TreeShrink

### 5.2 `phyloai doctor` Output Format

`phyloai doctor` defaults to `text` output. This is the **only** CLI command that supports `--output-format`; all other commands write structured results to `result.json`. Help text must make this explicit.

```
PhyloAI Environment Check
==========================
[OK]   iqtree3       3.0.1     /usr/local/bin/iqtree3
[OK]   mafft         7.520     /usr/local/bin/mafft
[OK]   magus         1.1.0     /usr/local/bin/magus
[OK]   trimal        1.4.1     bundled
[OK]   java          21.0.2    /usr/bin/java
[WARN] PhyloBayes              not found — CAT-GTR module unavailable
                                install: https://github.com/bayesiancook/pbmpi
[OK]   correction_multi.jl 1.0.0     bundled
[WARN] julia                   —         not found
                                install: https://julialang.org/downloads/
[WARN] MAGUS                   not found — large-dataset alignment
                                falls back to MAFFT
[MISS] MCMCTree                not found — dating module unavailable
                                install: https://github.com/abacus-gene/paml/releases
```

---

## 6. Output Directory Convention

```
runs/
├── pretree/
│   ├── align/
│   │   ├── seqs/
│   │   ├── align.log       # log co-located with outputs
│   │   └── result.json
│   ├── trim/
│   │   ├── seqs/
│   │   ├── trim.log
│   │   └── result.json
│   ├── metrics/
│   │   ├── metrics.log
│   │   └── result.json
│   ├── filter/
│   │   ├── seqs/
│   │   ├── filter.log
│   │   └── result.json
│   └── concat/
│       ├── concat.log
│       └── result.json
│   # note: convert defaults to runs/pretree/convert;
│   # stats writes reports to caller-specified files; neither writes a log
├── tree/
│   ├── genetree/
│   │   ├── genetree.log
│   │   └── result.json
│   ├── iqtree/
│   │   ├── iqtree.log
│   │   └── result.json
│   ├── astral/
│   │   ├── astral.log
│   │   └── result.json
│   └── phylobayes/
│       ├── phylobayes.log
│       └── result.json
├── posttree/
│   ├── concordance/
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

Log file content per step: tool version, full command, stderr, wall time, exit code, and stdout only when stdout is diagnostic text. Commands must not duplicate large primary data streams in logs when stdout is the primary output and that output is already saved to a dedicated file or subdirectory (e.g., aligned FASTA files under `seqs/`). Commands that produce multiple output files under a subdirectory (e.g., `seqs/`) place the log directly in the output directory root alongside `result.json`, without a separate `logs/` subdirectory.

---

## 7. Report Module

### 7.1 Design

- Each module writes a structured JSON result to its output directory on completion. The JSON contains two sections: `params` (all inputs and resolved parameters) and `key_results` (quantitative outputs and conclusions relevant for reporting).
- `report collector` aggregates all JSON files from a run directory.
- `report methods` renders a Methods paragraph from the aggregated record using a template engine.
- `report summary` outputs `run_record.yaml` — a complete, reproducible parameter snapshot plus key results from each step.
- `report figures` renders tables and plots from `key_results` data into `runs/report/figures/`. Only steps with meaningful visual output produce figures.

### 7.2 `key_results` Examples by Step

| Step | key_results content |
|------|---------------------|
| `pretree metrics` | MSA metric distributions, PIS vs length scatter data |
| `pretree filter` | retained/removed gene counts, removal reason breakdown |
| `pretree concat` | taxon × gene occupancy matrix, retained/dropped MSA counts, per-variant matrix stats |
| `pretree align` / `trim` | alignment length distribution before/after |
| `tree iqtree` / `genetree` | model selected, log-likelihood, tree file path |

`pretree stats` and `pretree convert` are utility commands and do not contribute to `report`.

### 7.3 Methods Paragraph Example

> Protein sequences were aligned using MAFFT 7.520 with the L-INS-i strategy. Alignments were trimmed using BMGE 1.12 with stringent parameters (−m BLOSUM90 −h 0.4). Error sites were detected and masked using TAPER. Trimmed alignments were concatenated into supermatrices at 75%, 90%, and 100% occupancy thresholds using PhyKIT. Phylogenetic trees were inferred using IQ-TREE3 3.0.1 under the PMSF model (LG+C60+F+R4), with node support estimated from 1,000 UFBoot2 replicates and 1,000 SH-aLRT replicates. Gene concordance factors (gCF) and site concordance factors (sCF) were calculated in IQ-TREE3 to quantify branch-level topological support.

---

## 8. AI Integration (Post-CLI Phase)

### 8.1 MCP Server

Expose all CLI commands as MCP tools with JSON I/O. Any MCP-compatible AI client (Claude Desktop, Cursor, OpenCode) can invoke phyloai operations directly.

### 8.2 AI Coding Assistant Skill

A generic skill (SKILL.md) compatible with any AI coding assistant that supports the superpowers plugin system (Claude Code, OpenCode, Cursor, Codex, etc.).

Provides workflow orchestration logic on top of MCP tools:
- parameter recommendation based on data characteristics
- result interpretation and biological explanation
- Methods paragraph generation trigger
- dynamic next-step suggestion

### 8.3 `posttree diagnose` as Skill Sub-workflow

Systematic error diagnosis is an iterative human-in-the-loop process:

```
observe brlen distribution
    → suspect LBA?
        → run CCA
            → confirm compositional bias?
                → run simulation to validate model effect
                    → decide: change model / remove taxa / recode
```

Each step's outcome determines the next step's strategy. This cannot be encoded as a single CLI command. Implementation plan:

- `syserror.py` exposes all atomic operations individually (each callable from CLI)
- Skill layer orchestrates them into a guided interactive sub-workflow
- CLI users can invoke any atomic operation manually

---

## 9. CLI Conventions and Parameter Standards

This section defines binding conventions for all subcommand implementations. Per-subcommand designs must follow these rules. Deviations require explicit justification in the subcommand spec.

### 9.1 Naming Style

- All parameter names use **kebab-case**: `--msa-dir`, `--seq-type`, `--output-dir`
- Directories always use the `--xxx-dir` suffix; single files omit the suffix: `--tree`, `--matrix`, `--calibrations`
- Short aliases (`-t`, `-o`) are reserved only for the highest-frequency parameters listed in Section 9.2; all other parameters use long form only to avoid cross-command conflicts

### 9.2 Shared Parameter Registry

Parameters that appear in multiple commands must use exactly these names and types:

| Parameter | Short | Type | Default | Applicable commands |
|-----------|-------|------|---------|---------------------|
| `--msa-dir` | | Path | — | all commands reading an alignment directory |
| `--tree-dir` | | Path | — | commands requiring a gene tree directory |
| `--output-dir` | `-o` | Path | auto under `runs/` | all commands producing output files |
| `--threads` | `-t` | int | 4 | all commands invoking multi-threaded tools |
| `--seq-type` | | AA \| NT \| auto | auto where safe; command-specific otherwise | commands where molecule type affects behavior |
| `--tool` | | str | tool-specific | commands offering multiple tool choices |
| `--input-format` | | str | auto-detect | all commands reading alignment files |
| `--output-format` | | text \| json | text | `doctor` only |
| `--tool-args` | | str | — | all commands invoking external tools; strategy parameters only |
| `--config` | | Path | — | all commands (global, inherited from CLI root) |
| `--run-dir` | | Path | `./runs/` | all commands except `stats` and `convert` (global, inherited from CLI root) |
| `--dry-run` | | flag | False | all commands |
| `--quiet` | `-q` | flag | False | all commands |
| `--overwrite` | | flag | False | all commands producing output directories |
| `--resume` | | flag | False | long-running pipeline commands only; resume from `checkpoint.json` |

Parameters marked "global, inherited from CLI root" are defined once on the `phyloai` root group and passed via Click context; subcommands do not re-declare them.

### 9.3 Exit Codes

All commands must exit with one of these codes:

| Code | Meaning |
|------|---------|
| `0` | Success |
| `1` | User input error (bad parameter, file not found, invalid value) |
| `2` | External tool execution failed (non-zero returncode from subprocess) |
| `3` | Environment error (required tool not installed or not detectable) |

### 9.4 JSON Output Schema

Every command writes a JSON result file to its output directory. The top-level structure is:

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

- `params`: all resolved input parameters (after config merge, after defaults applied)
- `key_results`: quantitative outputs and conclusions for report integration; empty `{}` for utility commands (stats, convert)
- `tool_versions`: version string for every external tool invoked
- `error`: null on success; error message string on failure
- `data`: command-specific detailed results (summary, files, skipped entries, warnings)

#### Unified Output Directory Structure

All commands follow the same output convention:

1. **JSON result** is always written to `result.json` inside the output directory
2. **Data files** (sequences, alignments, etc.) are written alongside `result.json`, or in a `seqs/` subdirectory for commands that produce sequence files
3. **Auxiliary files** (per-gene tables, etc.) are written alongside `result.json`

Example directory structure for `pretree convert`:
```
runs/pretree/convert/
├── seqs/                 # converted sequence files
│   ├── gene1.fa
│   ├── gene2.fa
│   └── ...
└── result.json           # JSON result
```

Example directory structure for `pretree stats`:
```
runs/pretree/stats/
├── result.json           # JSON result
└── per-gene.csv          # per-gene table (when --per-gene is used)
```

Example directory structure for pipeline commands (`tree iqtree`, etc.):
```
runs/tree/iqtree/
├── result.json           # JSON result
└── ...                   # tool output files
```

Example directory structure for `pretree align` (uses `seqs/` subdirectory, consistent with `pretree convert`):
```
runs/pretree/align/
├── seqs/                 # aligned sequence files (Mode 1/2: flat; Mode 3: faa/ and fna/)
│   ├── gene1.fa
│   └── ...
├── align.log             # per-gene tool log (co-located with outputs, not under runs/logs/)
└── result.json           # JSON result
```

This unified structure ensures that:
- Every output directory is self-contained with `result.json` for report aggregation
- `report collector` can aggregate results by scanning the run directory tree for `result.json` files
- Users can easily find all outputs from a command in one directory

### 9.5 Output Directory Conflict and Resume Policy

- Default behavior: if the output directory already exists and is non-empty, **exit with code 1** and a clear message
- `--overwrite`: delete and recreate the output directory before running
- `--resume`: for long-running pipeline commands, load `checkpoint.json` from the output directory, require exact parameter match, skip verified successful tasks, and rerun failed, interrupted, pending, or invalid-output tasks
- `--overwrite` and `--resume` are mutually exclusive
- Short utility commands such as `pretree convert` and `pretree stats` do not need resume support in the first checkpoint implementation

The detailed checkpoint schema and command adoption rules are defined in `docs/superpowers/specs/2026-06-12-checkpoint-resume-design.md`.

### 9.6 Display and Logging

- Terminal output always uses **Rich**: progress bars for batch operations, tables for summary results, colored status indicators. Non-`doctor` commands always write machine-readable structured output to `result.json`; they do not expose a text/json structured stdout switch.
- `--quiet` suppresses all terminal output except errors; useful for scripting and HPC batch jobs
- Every pipeline command (align, trim, metrics, filter, concat, genetree, iqtree, astral, phylobayes, and all posttree commands) writes a log file to its own output directory alongside `result.json`, named `<step>.log` (e.g., `runs/pretree/align/align.log`). The log contains: resolved command, tool versions, stderr, wall time, exit code, and stdout only when stdout is diagnostic text rather than the primary data output. Logs are never written to a separate `runs/<run>/logs/` directory. Utility commands (`stats`, `convert`) do not write run logs; `stats` is read-only, while `convert` writes normalized converted files to its `--output-dir` under the standard run layout by default.
- Log files are appended (not overwritten) on retry or resume, with a timestamp separator between runs. Resumed log entries should include resume context when useful. On `--overwrite` runs the log file is deleted and recreated together with the output directory.
- Every CLI command must provide **high-readability `--help` text**. Command help should explain what the command is for, when to use it, and what the major output means. Option help should explain practical intent, not just restate the flag name or type. Bare placeholders such as `TEXT`, missing descriptions, or one-line vague summaries are not acceptable for released commands.
- When a command writes one or more output files, terminal output must explicitly state what was saved and where, using concrete wording such as `Summary saved to <path>` or `Per-gene table saved to <path>`. Users should not need to infer output destinations from arguments alone.

### 9.7 Tabular Output Format

Commands that produce auxiliary tabular output (per-gene tables, per-taxon tables, metric tables) default to **CSV** format. TSV is available via command-specific table format flags such as `--per-gene-format tsv`, or by explicit output-file conventions when a command defines them. Structured command results always use `result.json`.

### 9.8 `--tool-args` Strategy-Only Semantics

1. PhyloAI builds managed input, output, work directory, data type, threads, logs, and codon/projection arguments first.
2. `--tool-args` is tokenized with standard shell splitting (respects quoted strings).
3. If `--tool-args` includes a PhyloAI-managed flag, the command exits with code 1 and names the blocked flag.
4. Non-managed strategy parameters are appended unchanged to the command and logged before execution.
5. Format of `--tool-args` must match the target tool's own CLI conventions; PhyloAI validates managed-flag conflicts but does not reimplement each tool's full parser.
6. Only `--tool-args` is supported for user-supplied external-tool strategy parameters.

### 9.9 Generated Sequence and Alignment Validation

- Commands that generate sequence or alignment files in bulk (`pretree convert`, `pretree align`, `pretree trim`, and related future steps) must validate generated files before counting them as successful outputs.
- Validation is implemented through shared `core` helpers, not one-off command-local checks. FASTA validation must detect missing/empty files, unparsable FASTA, zero FASTA records, and empty sequences. MSA validation must additionally detect unequal sequence lengths.
- Per-file validation failures are recorded in `result.json` under `data.skipped` or per-file warnings, depending on whether the command can continue. If all generated files fail validation, the command exits with code 1.
- Input scanning follows the same principle: empty files, directories, unrecognized file types, and unparsable sequence files are skipped with explicit reasons; if no valid inputs remain, the command exits with code 1.

---

## 10. Platform Support

- **Primary:** Linux, macOS
- **Secondary:** WSL (Windows Subsystem for Linux)
- **Not supported:** Native Windows
- Tool-specific exceptions are allowed when upstream binaries are not cross-platform. In Phase 2, `pretree align --method magus` is Linux-only because the pip-distributed MAGUS dependency bundle includes Linux binaries.

---

## 11. Development Phases

| Phase | Scope | Deliverable | Pre-requisites |
|-------|-------|-------------|----------------|
| 1 | `core/` infrastructure | env, runner, formats, logger | — |
| 2 | `pretree/` modules | stats, convert, align, trim, metrics, filter, concat | Phase 1 |
| 3 | `tree/` modules | genetree, iqtree, astral, phylobayes | Phase 1 |
| 4 | `posttree/` modules | concordance, topology, dating, signal, syserror, simulate | Phases 2–3 |
| 5 | `phyloai run` | one-click supermatrix and coalescent pipelines | Phases 2–3 (align, trim, filter/TAPER, concat, genetree, iqtree, astral) |
| 6 | `report/` module | collector, methods generation, figures rendering, run_record.yaml | Phases 2–4 |
| 7 | MCP Server | JSON tool interface | Phases 2–6 (full CLI stable) |
| 8 | AI Coding Assistant Skill | workflow orchestration, syserror guided sub-workflow | Phase 7 |

**Spec and plan granularity by phase:**

- **Phase 1** — one spec + plan covering the entire `core/` infrastructure (already complete).
- **Phases 2–4** — every subcommand (`pretree stats`, `pretree align`, `tree iqtree`, `posttree concordance`, etc.) has its own design spec and implementation plan under `docs/superpowers/specs/` and `docs/superpowers/plans/` respectively. Subcommands are designed and implemented one at a time in pipeline order within each phase.
- **Phase 5** — one spec + plan for `phyloai run` (pipeline orchestration only; does not include report).
- **Phase 6** — one spec + plan for the `report/` module (collector, methods, figures, run_record); depends on Phases 2–4 being stable so JSON output schemas are finalized.
- **Phase 7** — one spec + plan for the MCP Server.
- **Phase 8** — one spec + plan for the AI Coding Assistant Skill.

Each subcommand spec and plan must be consistent with the conventions in Section 9; any deviation requires explicit justification in the subcommand spec.

Every subcommand implementation must comply with the following binding requirements from the main design before it can be considered complete:

1. **Structured output** for all non-`doctor` commands is always saved as `result.json`; terminal Rich display is always on unless `--quiet` is set.
2. **Output directory structure** follows Section 6; commands with an `--output-dir` parameter default to the appropriate subdirectory under `runs/`.
3. **Conflict policy** follows Section 9.5: non-empty output directory → exit 1 by default; `--overwrite` replaces; supported long-running commands may use `--resume` with a valid `checkpoint.json`.
4. **Log file** written to the command output directory as `<step>.log` for all pipeline commands (Section 9.6).
5. **JSON output schema** follows Section 9.4; `result.json` written to the output directory on success.
6. **Exit codes** follow Section 9.3.

Modules within each phase can be developed in parallel. Phases are strictly sequential.

---

## 12. Key Design Decisions

| Decision | Rationale |
|----------|-----------|
| No GUI | CLI + AI interaction covers all use cases; GUI adds maintenance cost without proportional benefit in this era |
| stats/convert as pretree subcommands | Semantically these are pre-analysis utilities, consistent with concat being a data preparation step; keeps top-level CLI clean |
| filter.py unifies MSA + tree + TAPER filtering | Single entry point reduces cognitive load; TAPER as lightweight option within filter avoids an orphan command |
| metrics is prerequisite for filter | Decouples metric computation from filtering decisions; users can inspect metrics before committing to a filter strategy |
| metrics layered (core + --advanced) | Core metrics are fast and always needed; UMAP/correlation are expensive and optional |
| --tool-args strategy-only model | Keeps batch input/output and result schemas deterministic while still exposing tool-specific strategy knobs |
| Format handling per-module, not global | Different tools require different formats; per-module handling via core/formats.py is more correct than forcing a global FASTA mandate |
| backtrans in align, not a separate command | It is a direct post-processing of the alignment step and uses trimAl which is already a dependency |
| syserror.py exposed as atomic ops only | Full diagnosis requires iterative human decisions; CLI atomics + Skill orchestration is the correct separation |
| genetree.py in tree/ not pretree/ | Gene trees are tree inference results, not preprocessing steps |
| model_eval not exposed as standalone | Model evaluation logic is internal to syserror; exposing it separately creates redundant interface surface |
| JSON result files for non-`doctor` commands | Ensures MCP wrapping has one stable machine-readable result path without maintaining parallel text/json command outputs; `doctor` is the sole command with `--output-format` because it is a human-oriented diagnostic |
| JSON key_results in all pipeline modules | Enables report figures and summary without post-hoc data extraction; schema defined at design time |
| YAML for config, JSON for output | YAML supports comments and is human-writable; JSON is strict and machine-parseable |
| `--input-format` on alignment-reading commands | Real datasets often use inconsistent suffixes; explicit user intent must override guessing |
| Per-subcommand design + plan docs | Each pretree subcommand is complex enough to warrant its own spec; keeps main design doc stable while allowing detailed iteration |
| Lightweight top-level README + command docs | Keeps `README.md` maintainable as the CLI grows; detailed command behavior belongs in `docs/commands/*.md` with a consistent section structure |

---

## 13. YAML Config Files

All `--config FILE` parameters should map 1-to-1 to CLI flags. Unknown keys should be ignored with a warning, and CLI flags should override config values.

Example YAML templates should live under `examples/` in the repository rather than being embedded in the design document. Those example files should cover global options and representative command-level inputs, including `run_dir` and future `input_format` fields for alignment-reading commands.
