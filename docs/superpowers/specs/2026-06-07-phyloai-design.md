# PhyloAI Design Specification

**Date:** 2026-06-07  
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
│   ├── align.py        # MAFFT (external), MAGUS (pip), batch AA + NT
│   ├── trim.py         # trimAl (bundled), BMGE (bundled), ClipKIT (pip)
│   ├── filter.py       # unified MSA-based + gene-tree-based marker filtering
│   │                   #   MSA-based:  PIS, length, GC content, nRCFV,
│   │                   #               symtest, API, likelihood mapping, TAPER
│   │                   #   Tree-based: TreeShrink, ABS, treeness, DVMC,
│   │                   #               evo_rate, saturation, inconsistent genes
│   ├── metrics.py      # MSA + tree attribute extraction (16 MSA + 14 tree metrics),
│   │                   # correlation analysis, outlier detection (UMAP clustering)
│   │                   # reference: github.com/xtmtd/MSA-and-tree-metrics-exploration
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
phyloai pretree align    --msa-dir ./raw --method linsi [--nt-dir ./raw_nt] \
                         [--input-format fasta]
phyloai pretree trim     --msa-dir ./aligned --tool bmge [--model BLOSUM90] \
                         [--input-format phylip-relaxed]
phyloai pretree metrics  --msa-dir ./trimmed [--tree-dir ./genetrees]
phyloai pretree filter   --msa-dir ./trimmed [--tree-dir ./genetrees] \
                         --strategy outlier [--filter-by pis,abs,treeness]
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

Two modes, both include: `align → trim → TAPER → concat → [tree inference]`

| Mode | Steps | Notes |
|------|-------|-------|
| `--mode supermatrix` | align → trim → TAPER → concat → iqtree (unpartitioned) | Fast, single-step ML tree |
| `--mode coalescent` | align → trim → TAPER → concat → genetree → astral-hybrid | MSC-based species tree |

Filter step in `phyloai run` is TAPER only (error site detection). Full marker filtering is an explicit manual step via `phyloai pretree filter`.

### 4.3 Universal CLI Flags

All commands support:

| Flag | Purpose |
|------|---------|
| `--output-format json` | Machine-readable output (MCP pre-requisite) |
| `--dry-run` | Show what would be executed without running |
| `--config FILE` | Load parameters from YAML (HPC batch use); see Section 12 for template |
| `--threads N` | Parallelism control |
| `--run-dir DIR` | Override default run directory |

Commands that read alignment files should also support `--input-format`.

- Accepted values should match the core `AlignmentFormat` enum
- If `--input-format` is given, PhyloAI should trust it and skip guessing
- If it is omitted, format detection should use filename suffix first, then fall back to conservative content inspection
- Common suffixes to support include `.fa`, `.fasta`, `.fas`, `.faa`, `.fna`, `.phy`, `.phylip`, `.nex`, `.nxs`, and `.nexus`

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

- Each module writes a structured JSON result to its output directory on completion.
- `report collector` aggregates all JSON files from a run directory.
- `report methods` renders a Methods paragraph from the aggregated record using a template engine.
- `report summary` outputs `run_record.yaml` — a complete, reproducible parameter snapshot.

### 7.2 Methods Paragraph Example

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

## 9. Platform Support

- **Primary:** Linux, macOS
- **Secondary:** WSL (Windows Subsystem for Linux)
- **Not supported:** Native Windows

---

## 10. Development Phases

| Phase | Scope | Deliverable |
|-------|-------|-------------|
| 1 | `core/` infrastructure | env, runner, formats, logger |
| 2 | `pretree/` modules | align, trim, filter, metrics, concat |
| 3 | `tree/` modules | genetree, iqtree, astral, phylobayes |
| 4 | `posttree/` modules | concordance, topology, dating, signal, syserror, simulate |
| 5 | `report/` + `phyloai run` | methods generation, one-click pipeline |
| 6 | MCP Server | JSON tool interface |
| 7 | AI Coding Assistant Skill | workflow orchestration, syserror guided sub-workflow |

Modules within each phase can be developed in parallel. Phases are strictly sequential.

---

## 11. Key Design Decisions

| Decision | Rationale |
|----------|-----------|
| No GUI | CLI + AI interaction covers all use cases; GUI adds maintenance cost without proportional benefit in this era |
| filter.py unifies MSA + tree filtering | Gene trees are a prerequisite for tree-based filtering anyway; single entry point reduces user cognitive load |
| syserror.py exposed as atomic ops only | Full diagnosis requires iterative human decisions; CLI atomics + Skill orchestration is the correct separation |
| genetree.py in tree/ not pretree/ | Gene trees are tree inference results, not preprocessing steps |
| model_eval not exposed as standalone | Model evaluation logic is internal to syserror; exposing it separately creates redundant interface surface |
| TAPER only in `phyloai run` | TAPER is fast, non-destructive, and universally beneficial; full filtering requires data-specific decisions |
| `--output-format json` from day one | Ensures MCP wrapping requires zero interface changes later |
| YAML for config, JSON for output | YAML supports comments and is human-writable; JSON is strict and machine-parseable |
| `--input-format` on alignment-reading commands | Real datasets often use inconsistent suffixes; explicit user intent must override guessing |

---

## 12. YAML Config Files

All `--config FILE` parameters should map 1-to-1 to CLI flags. Unknown keys should be ignored with a warning, and CLI flags should override config values.

Example YAML templates should live under `examples/` in the repository rather than being embedded in the design document. Those example files should cover global options and representative command-level inputs, including `output_format`, `run_dir`, and future `input_format` fields for alignment-reading commands.
