# PhyloAI Design Specification

**Date:** 2026-06-07  
**Last updated:** 2026-06-09 (pretree module design finalized)  
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

**Development order:** Library + CLI first → MCP Server + Skill after CLI is stable.  
**MCP pre-requisite:** All CLI commands support `--output-format json` and standard exit codes from day one, so MCP wrapping requires no interface redesign.

---

## 3. Module Structure

```
phyloai/
├── core/
│   ├── env.py          # tool detection, path resolution, bundled tool management
│   ├── runner.py       # unified external tool call interface, timeout, retry
│   ├── formats.py      # format detection & conversion (FASTA/Nexus/Phylip/Phylip-PAML)
│   ├── schema.py       # shared data structures (MSACollection, TreeSet, RunRecord)
│   └── logger.py       # per-step log files under runs/runNNN/logs/
│
├── pretree/
│   ├── stats.py        # sequence/alignment statistics: format, length distribution,
│   │                   # taxon count, gap ratio, min/max/mean/median length, total length
│   ├── convert.py      # format conversion between FASTA/Phylip/Nexus/Phylip-PAML;
│   │                   # wraps core/formats.py, exposes as CLI subcommand
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
phyloai pretree stats    --msa-dir ./aligned
phyloai pretree convert  --input ./file.phy --output ./file.fa --format fasta
phyloai pretree align    --msa-dir ./raw --method linsi [--nt-dir ./raw_nt] \
                         [--backtrans] [--extra-args "--op 1 --bl"]
phyloai pretree trim     --msa-dir ./aligned --tool bmge [--model BLOSUM90] \
                         [--extra-args "-h 0.4"]
phyloai pretree metrics  --msa-dir ./trimmed [--tree-dir ./genetrees] [--advanced]
phyloai pretree filter   --msa-dir ./trimmed --metrics-dir ./metrics \
                         [--tree-dir ./genetrees] \
                         --strategy outlier [--filter-by pis,abs,treeness] \
                         [--taper]
phyloai pretree concat   --msa-dir ./filtered --occupancy 75,80,90,100

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
phyloai report generate --run-dir ./runs/run001

# One-click pipeline
phyloai run --msa-dir ./raw --output ./runs/run001 --mode supermatrix
phyloai run --msa-dir ./raw --output ./runs/run001 --mode coalescent
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
| `--output-format json\|text` | Machine-readable output (MCP pre-requisite); default `text` |
| `--dry-run` | Show what would be executed without running |
| `--config FILE` | Load parameters from YAML (HPC batch use); see Section 13 for template |
| `--threads` / `-t` | Parallelism control |
| `--run-dir` | Override default run directory |
| `--quiet` / `-q` | Suppress all output except errors |
| `--overwrite` | Overwrite existing output directory |

Commands that read alignment files also support `--input-format` (auto-detect by default; accepted values match `AlignmentFormat` enum). Commands that invoke external tools also support `--extra-args` (see Section 9.7 for merge semantics).

**Format handling policy:** Each module uses the format required by its underlying tool. `pretree convert` and `core/formats.py` provide conversion as needed. There is no global FASTA-only mandate; modules handle format requirements internally.

Full parameter naming rules and exit code definitions are in **Section 9**.

---

## 5. Dependency Management

### 5.1 Strategy

- **Bundled (auto-downloaded on install):** trimAl, BMGE — permissive licenses, small binaries
- **pip-installable (auto-installed):** MAGUS (`pip install magus-msa`), ClipKIT, PhyKIT
- **User-installed (detected by `doctor`):** IQ-TREE3, PhyloBayes-MPI, astral-hybrid, MAFFT, MCMCTree (PAML), TreeShrink, TAPER (`correction_multi.jl`)

### 5.2 `phyloai doctor` Output Format

`phyloai doctor` should default to text output. Help text must make that default explicit.

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
[WARN] correction_multi.jl     —         not found
                                install: https://github.com/chaoszhang/TAPER
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
└── run001/
    ├── pretree/
    │   ├── align/
    │   ├── trim/
    │   ├── metrics/
    │   ├── filter/
    │   └── concat/
    │   # note: stats and convert write to caller-specified paths,
    │   # not to a fixed pretree subdirectory
    ├── tree/
    │   ├── genetree/
    │   ├── iqtree/
    │   ├── astral/
    │   └── phylobayes/
    ├── posttree/
    │   ├── concordance/
    │   ├── topology/
    │   ├── dating/
    │   ├── signal/
    │   ├── syserror/
    │   └── simulate/
    ├── report/
    │   ├── summary.json
    │   ├── methods.txt
    │   └── run_record.yaml
    └── logs/
        ├── align.log
        ├── trim.log
        ├── filter.log
        ├── concat.log
        ├── iqtree.log
        └── ...         # one log per step, appended on retry
```

Log file content per step: tool version, full command, stdout/stderr, wall time, exit code.

---

## 7. Report Module

### 7.1 Design

- Each module writes a structured JSON result to its output directory on completion. The JSON contains two sections: `params` (all inputs and resolved parameters) and `key_results` (quantitative outputs and conclusions relevant for reporting).
- `report collector` aggregates all JSON files from a run directory.
- `report methods` renders a Methods paragraph from the aggregated record using a template engine.
- `report summary` outputs `run_record.yaml` — a complete, reproducible parameter snapshot plus key results from each step.
- `report figures` renders tables and plots from `key_results` data into `runs/runNNN/report/figures/`. Only steps with meaningful visual output produce figures.

### 7.2 `key_results` Examples by Step

| Step | key_results content |
|------|---------------------|
| `pretree metrics` | MSA metric distributions, PIS vs length scatter data |
| `pretree filter` | retained/removed gene counts, removal reason breakdown |
| `pretree concat` | taxon × gene occupancy matrix, per-occupancy-level matrix stats |
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
| `--seq-type` | | AA \| NT | AA | commands where molecule type affects behavior |
| `--tool` | | str | tool-specific | commands offering multiple tool choices |
| `--input-format` | | str | auto-detect | all commands reading alignment files |
| `--output-format` | | text \| json | text | all commands (MCP pre-requisite) |
| `--extra-args` | | str | — | all commands invoking external tools |
| `--config` | | Path | — | all commands (global, inherited from CLI root) |
| `--run-dir` | | Path | `./runs/` | all commands (global, inherited from CLI root) |
| `--dry-run` | | flag | False | all commands |
| `--quiet` | `-q` | flag | False | all commands |
| `--overwrite` | | flag | False | all commands producing output directories |

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

When `--output-format json` is set, every command outputs a single JSON object to stdout with this top-level structure:

```json
{
  "status": "success | error",
  "command": "phyloai pretree align ...",
  "wall_time": 142.3,
  "tool_versions": {"mafft": "7.520"},
  "params": { },
  "key_results": { },
  "error": null
}
```

- `params`: all resolved input parameters (after config merge, after defaults applied)
- `key_results`: quantitative outputs and conclusions for report integration; empty `{}` for utility commands (stats, convert)
- `tool_versions`: version string for every external tool invoked
- `error`: null on success; error message string on failure

Each module also writes this same JSON to its output directory as `result.json` for `report collector` to aggregate.

### 9.5 Output Directory Conflict Policy

- Default behavior: if the output directory already exists and is non-empty, **exit with code 1** and a clear message
- `--overwrite`: delete and recreate the output directory before running
- `--resume` is not in scope for Phase 2; may be added for long-running commands (align, genetree) in a later phase

### 9.6 Display and Logging

- Terminal output uses **Rich**: progress bars for batch operations, tables for summary results, colored status indicators
- `--quiet` suppresses all terminal output except errors; useful for scripting and HPC batch jobs
- Every command writes a log file to `runs/runNNN/logs/<step>.log` containing: resolved command, tool versions, full stdout/stderr, wall time, exit code
- Log files are appended (not overwritten) on retry, with a timestamp separator between runs

### 9.7 Tabular Output Format

Commands that produce tabular output (per-gene tables, per-taxon tables, metric tables) default to **CSV** format. TSV is available via `--output-format tsv` or by using a `.tsv` file extension with `--output`. JSON output always uses `--output-format json` or a `.json` extension.

### 9.8 `--extra-args` Merge Semantics

1. PhyloAI builds its internal parameter set using tool-native argument format
2. `--extra-args` string is tokenized with standard shell splitting (respects quoted strings)
3. Any parameter in `--extra-args` that conflicts with an internal parameter **replaces** it (extra-args win)
4. The fully merged command is logged before execution
5. Format of `--extra-args` must match the target tool's own CLI conventions; PhyloAI does not validate tool-specific argument semantics

---

## 10. Platform Support

- **Primary:** Linux, macOS
- **Secondary:** WSL (Windows Subsystem for Linux)
- **Not supported:** Native Windows

---

## 11. Development Phases

| Phase | Scope | Deliverable |
|-------|-------|-------------|
| 1 | `core/` infrastructure | env, runner, formats, logger |
| 2 | `pretree/` modules | stats, convert, align, trim, metrics, filter, concat |
| 3 | `tree/` modules | genetree, iqtree, astral, phylobayes |
| 4 | `posttree/` modules | concordance, topology, dating, signal, syserror, simulate |
| 5 | `report/` + `phyloai run` | methods generation, figures rendering, one-click pipeline |
| 6 | MCP Server | JSON tool interface |
| 7 | AI Coding Assistant Skill | workflow orchestration, syserror guided sub-workflow |

Each `pretree` subcommand has its own design spec and implementation plan under `docs/superpowers/specs/` and `docs/superpowers/plans/` respectively. Subcommands within Phase 2 are designed and implemented one at a time in pipeline order.

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
| --extra-args with extra-wins merge | Avoids duplicate parameters in final command; behavior is deterministic and logged; users follow tool's own documentation |
| Format handling per-module, not global | Different tools require different formats; per-module handling via core/formats.py is more correct than forcing a global FASTA mandate |
| backtrans in align, not a separate command | It is a direct post-processing of the alignment step and uses trimAl which is already a dependency |
| syserror.py exposed as atomic ops only | Full diagnosis requires iterative human decisions; CLI atomics + Skill orchestration is the correct separation |
| genetree.py in tree/ not pretree/ | Gene trees are tree inference results, not preprocessing steps |
| model_eval not exposed as standalone | Model evaluation logic is internal to syserror; exposing it separately creates redundant interface surface |
| `--output-format json` from day one | Ensures MCP wrapping requires zero interface changes later |
| JSON key_results in all pipeline modules | Enables report figures and summary without post-hoc data extraction; schema defined at design time |
| YAML for config, JSON for output | YAML supports comments and is human-writable; JSON is strict and machine-parseable |
| `--input-format` on alignment-reading commands | Real datasets often use inconsistent suffixes; explicit user intent must override guessing |
| Per-subcommand design + plan docs | Each pretree subcommand is complex enough to warrant its own spec; keeps main design doc stable while allowing detailed iteration |

---

## 13. YAML Config Files

All `--config FILE` parameters should map 1-to-1 to CLI flags. Unknown keys should be ignored with a warning, and CLI flags should override config values.

Example YAML templates should live under `examples/` in the repository rather than being embedded in the design document. Those example files should cover global options and representative command-level inputs, including `output_format`, `run_dir`, and future `input_format` fields for alignment-reading commands.
