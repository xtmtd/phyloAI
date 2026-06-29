# phyloai doctor

[English](doctor.md) | [中文](doctor.zh.md)


## Purpose

`phyloai doctor` checks whether external tools expected by the local PhyloAI environment can be found, where they are located, and whether a version string can be detected.

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

All third-party tools are environment-dependent and should be installed by the user for their operating system, package manager, cluster, or Conda environment. BMGE is detected as `BMGE.jar`; TAPER is detected as `correction_multi.jl`.

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

If `TAPER` (`correction_multi.jl`) or `BMGE` (`BMGE.jar`) is reported missing, install it externally or pass its explicit path to the command that needs it.

## Notes

The current registry includes required tools such as `iqtree3`, `mafft`, and `trimal`, plus optional tools such as `wastral`, `pb_mpi`, `mcmctree`, `run_treeshrink.py`, `magus`, `clipkit`, `BMGE.jar`, `correction_multi.jl`, `java`, `julia`, and `FastTree`.
