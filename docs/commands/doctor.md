# phyloai doctor

## Purpose

`phyloai doctor` checks whether external tools expected by the local PhyloAI environment can be found, where they are located, and whether a version string can be detected. Project-bundled tools are reported from the PhyloAI package path instead of requiring users to install them separately.

It does not install tools, modify the environment, or validate input datasets.

## Usage

```bash
phyloai doctor [--output-format text|json]
```

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--output-format text|json` | `text` | Choose human-readable terminal output or machine-readable JSON. `doctor` is the only command that defaults to `text`. |

## Inputs

`doctor` has no sequence or tree input files. It reads the current shell environment, including `PATH`, and checks whether known external tools are visible from that environment.

TAPER 1.0.0 (`correction_multi.jl`) and BMGE 1.12 (`BMGE.jar`) are bundled inside the PhyloAI package and do not need separate installation. Other external tools remain environment-dependent and should be installed by the user for their operating system, package manager, cluster, or Conda environment.

For practical setup guidance, see [installation.md](installation.md). It lists Python environment options, external tool groups, `phyloai run` dependency modes, and operating-system notes.

## Outputs

Text output is a Rich table showing each tool, status, detected version, path, and install note when relevant.

JSON output is a mapping from tool name to structured status fields such as `status`, `path`, `version`, and `note`.

## Examples

Human-readable environment check:

```bash
phyloai doctor
```

JSON for scripting or CI:

```bash
phyloai doctor --output-format json
```

## Warnings And Errors

A missing optional tool is reported but does not mean the whole installation is unusable. Later workflow steps that require that tool may fail or become unavailable.

If a tool is missing in `doctor`, PhyloAI is seeing the same environment that later CLI commands will see. Activate the intended Conda or virtual environment before re-running the command.

`TAPER` (`correction_multi.jl`) and `BMGE` (`BMGE.jar`) are bundled with PhyloAI and do not need to be installed separately. If either is reported missing, it indicates a packaging or installation problem rather than a missing user-installed dependency.

## Notes

The current registry includes required tools such as `iqtree3`, `mafft`, and `trimal`, plus optional tools such as `wastral`, `pb_mpi`, `mcmctree`, `run_treeshrink.py`, `magus`, `clipkit`, `java`, `julia`, and `FastTree`.

Bundled tools include `phyloai/bundled/TAPER-1.0.0/correction_multi.jl` and `phyloai/bundled/BMGE-1.12/BMGE.jar`. The bundled directories retain license and source files needed for distribution compliance. `bmge` is displayed as `BMGE.jar` in the table because PhyloAI calls the bundled jar file directly.
