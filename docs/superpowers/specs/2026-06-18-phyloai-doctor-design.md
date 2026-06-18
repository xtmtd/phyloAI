# PhyloAI Doctor Design Specification

**Date:** 2026-06-18  
**Status:** Draft  
**Parent spec:** `2026-06-07-phyloai-design.md`

---

## 1. Purpose

`phyloai doctor` checks whether external tools expected by PhyloAI can be found, where they are located, and whether a version string can be detected. Project-bundled tools are reported from the PhyloAI package path.

`phyloai doctor` is the **only** CLI command that supports `--output-format text|json`. It defaults to `text` because it is primarily a human-facing environment diagnostic. All other commands write structured results to `result.json`.

---

## 2. CLI

```bash
phyloai doctor [--output-format text|json]
```

### Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--output-format text\|json` | `text` | Human-readable Rich table or machine-readable JSON mapping tool name to structured status fields. |

---

## 3. Dependency Strategy

- **Bundled in the Python package:** TAPER 1.0.0 (`phyloai/bundled/TAPER-1.0.0/correction_multi.jl`), BMGE 1.12 (`phyloai/bundled/BMGE-1.12/BMGE.jar`) with corresponding license/source files retained in the bundled directories
- **Bundled or auto-downloaded later (planned):** IQ-TREE3, trimAl — exact packaging strategy to be finalized before `tree ml` and `pretree trim` respectively
- **pip-installable (auto-installed):** MAGUS (`pip install magus-msa`), ClipKIT
- **User-installed (detected by `doctor`):** PhyloBayes-MPI, wASTRAL, MAFFT, MCMCTree (PAML), TreeShrink
- **Runtime dependencies (detected by `doctor`):** Java (`java`), Julia

---

## 4. Tool Registry

The current registry includes required tools such as `iqtree3`, `mafft`, and `trimal`, plus optional tools such as `wastral`, `pb_mpi`, `mcmctree`, `run_treeshrink.py`, `magus`, `clipkit`, `java`, `julia`, and `FastTree`.

Bundled tools include `phyloai/bundled/TAPER-1.0.0/correction_multi.jl` and `phyloai/bundled/BMGE-1.12/BMGE.jar`. The bundled directories retain license and source files needed for distribution compliance. `bmge` is displayed as `BMGE.jar` in the table because PhyloAI calls the bundled jar file directly.

### Tool registry table (in `phyloai/core/env.py`)

| Tool | Required | Bundled | Install hint |
|------|----------|---------|-------------|
| `iqtree3` | Yes | Planned | github.com/iqtree/iqtree3/releases |
| `mafft` | Yes | No | mafft.cbrc.jp |
| `trimal` | Yes | Planned | github.com/inab/trimal/releases |
| `wastral` | No | No | github.com/chaoszhang/ASTER |
| `pb_mpi` | No | No | github.com/bayesiancook/pbmpi |
| `mcmctree` | No | No | github.com/abacus-gene/paml/releases |
| `correction_multi.jl` | No | Yes (TAPER-1.0.0) | bundled |
| `run_treeshrink.py` | No | No | github.com/uym2/TreeShrink |
| `magus` | No | No | pip install magus-msa |
| `clipkit` | No | No | pip install clipkit |
| `java` | No | No | java.com |
| `julia` | No | No | julialang.org |
| `bmge` | No | Yes (BMGE-1.12) | bundled |
| `FastTree` | No | No | microbesonline.org/fasttree |

---

## 5. Output Format

### Text output (default)

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
[MISS] MCMCTree                not found — dating module unavailable
                                install: https://github.com/abacus-gene/paml/releases
```

### JSON output

A mapping from tool name to structured status fields (`status`, `path`, `version`, `note`).

---

## 6. Cross-Cutting Conventions

Doctor inherits and implements the following conventions from the main design:

- **Universal flags:** `--output-format text|json` — `doctor` only command with this flag
- **Exit codes:** 0 on success, 3 if required tool is missing
- **Display:** Rich table for text output, compact JSON for machine consumption

---

## 7. Warnings and Errors

A missing optional tool is reported but does not mean the whole installation is unusable. Later workflow steps that require that tool may fail or become unavailable.

If a tool is missing in `doctor`, PhyloAI is seeing the same environment that later CLI commands will see. Activate the intended Conda or virtual environment before re-running the command.

`TAPER` (`correction_multi.jl`) and `BMGE` (`BMGE.jar`) are bundled with PhyloAI and do not need to be installed separately. If either is reported missing, it indicates a packaging or installation problem rather than a missing user-installed dependency.

---

## 8. Implementation

- **`phyloai/core/env.py`:** `TOOL_REGISTRY` dict, `ToolEnv` class for detection, path resolution, bundled tool management
- **`phyloai/cli/doctor.py`:** Click command, Rich table rendering, JSON serialization
- **`phyloai/cli/main.py`:** Registers `doctor` subcommand on the root `phyloai` group

Detailed implementation plan: `docs/superpowers/plans/2026-06-07-core-module.md` (Task 3: env.py, Task 7: CLI doctor)
