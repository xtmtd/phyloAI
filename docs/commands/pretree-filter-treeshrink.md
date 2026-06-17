# phyloai pretree filter treeshrink

## Purpose

`phyloai pretree filter treeshrink` detects and removes outlier long-branch taxa from gene trees using TreeShrink. When `--msa-dir` is provided, matching MSAs are also shrunk to remove the same pruned taxa.

TreeShrink is run once across the entire gene-tree dataset (not per gene) because it can use information from multiple trees jointly. PhyloAI creates a per-gene working layout (`input.tree`, optional `input.fasta`) in a temporary directory, invokes `run_treeshrink.py`, then collects the shrunk outputs.

## Usage

Minimal:
```bash
phyloai pretree filter treeshrink --tree-dir ./genetrees
```

With MSA pruning:
```bash
phyloai pretree filter treeshrink --tree-dir ./genetrees --msa-dir ./trimmed --threshold 0.1
```

## Parameters

| Parameter | Default | Notes |
|-----------|---------|-------|
| `--tree-dir` | required | Directory of input gene tree files (any suffix) |
| `--output-dir` / `-o` | `runs/pretree/filter/treeshrink` | Output directory |
| `--msa-dir` | — | Optional MSA directory for paired alignment pruning |
| `--threshold` | 0.05 | TreeShrink `-q` false-positive threshold; smaller = more taxa removed |
| `--treeshrink-mode` | `auto` | `auto` (omit `-m`), `per-gene`, `all-genes`, `per-species` |
| `--treeshrink-path` | — | Explicit path to `run_treeshrink.py`; resolved via PATH by default |
| `--tool-args` | — | Additional TreeShrink flags; PhyloAI manages `-i,-t,-a,-q,-m,-o,-O` |
| `--keep-work-dir` | off | Retain per-gene working directory under `--output-dir/work/` |
| `--table-format` | `csv` | `csv` or `tsv` for auxiliary tables |
| `--overwrite` | off | Delete and recreate non-empty output directory |
| `--dry-run` | off | Print command, create no files |
| `--quiet` / `-q` | off | Suppress terminal output except errors |

## Inputs

`--tree-dir` is scanned for gene tree files (any suffix, any format TreeShrink accepts). When `--msa-dir` is provided, MSAs are paired with trees by logical locus name.

## Outputs

```
runs/pretree/filter/treeshrink/
├── trees/
├── seqs/                         (only when --msa-dir provided)
├── retained_loci.csv|tsv
├── modified_loci.csv|tsv
├── dropped_loci.csv|tsv
├── removed_taxa.csv|tsv
├── filter_decisions.csv|tsv
├── work/                         (only with --keep-work-dir)
├── filter.log
└── result.json
```

Decision categories:
- `retained_loci`: loci with valid shrunk outputs (including unmodified)
- `modified_loci`: retained loci where TreeShrink removed at least one taxon
- `dropped_loci`: loci with missing or invalid outputs
- `removed_taxa`: per-taxon record of pruned taxa (`locus`, `taxon`)

Terminal summary includes input/retained/modified/dropped counts, taxa removed total, and retained MSA statistics when `--msa-dir` was provided. An indeterminate progress bar is shown during execution.

A tip is displayed reminding users that filtered alignments may be used to re-construct phylogenetic trees, which are possibly more accurate than those pruned by TreeShrink.

## Examples

```bash
# Basic taxon pruning
phyloai pretree filter treeshrink --tree-dir ./genetrees

# With MSA pruning and stricter threshold
phyloai pretree filter treeshrink --tree-dir ./genetrees --msa-dir ./trimmed --threshold 0.1

# Per-species mode with debugging work dir kept
phyloai pretree filter treeshrink --tree-dir ./genetrees --treeshrink-mode per-species --keep-work-dir
```

## Warnings and Errors

| Condition | Behaviour |
|-----------|-----------|
| `run_treeshrink.py` not found | Exit 3 |
| Non-empty output directory without `--overwrite` | Exit 1 |
| No valid tree input files | Exit 1 |
| Ambiguous locus matching | Exit 1 with details |
| TreeShrink exits non-zero | All loci marked failed |
| All loci failed | Exit 2 |

## Notes

- `--tool-args` passes strategy options. Do not include PhyloAI-managed flags: `-i`, `-t`, `-a`, `-q`, `-m`, `-o`, `-O`.
- TreeShrink mode `auto` omits `-m` entirely (TreeShrink default). `per-gene` runs independently per gene; `all-genes` and `per-species` use cross-gene information.
- `filter.log` records the resolved command, status, wall time, and per-locus outcomes.
