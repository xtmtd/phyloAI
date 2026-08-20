# phyloai posttree topology

[English](posttree-topology.md) | [中文](posttree-topology.zh.md)


## Purpose

Performs IQ-TREE tree topology tests (AU / KH / SH / WKH / WSH / c-ELW)
comparing a set of candidate trees against a supermatrix alignment. This
command tests whether alternative topologies are significantly worse than
the best-scoring candidate -- it does **not** infer new trees.

## Usage

```bash
# Homogeneous model
phyloai posttree topology --matrix matrix.fa --candidate-trees candidates.trees --model-expr LG+F+R4

# PMSF model with guide tree
phyloai posttree topology --matrix matrix.fa --candidate-trees candidates.trees --model-expr LG+C20+F+R4 --guide-tree guide.nwk

# Previously optimized partition model
phyloai posttree topology --matrix matrix.fa --candidate-trees candidates.trees --partitions matrix.best_model.nex

# Multiple individual tree files (comma-separated, merged by PhyloAI)
phyloai posttree topology --matrix matrix.fa --candidate-trees h1.nwk,h2.nwk,h3.nwk --model-expr LG+F+R4

# Custom exchangeabilities + site frequencies via --tool-args
phyloai posttree topology --matrix matrix.fa --candidate-trees trees --model-expr custom.exchangeabilities+R4 --tool-args "-fs custom.sitefreq" -t 30

# Heterogeneous model
phyloai posttree topology --matrix matrix.fa --candidate-trees trees --model-expr C20+F+R4
```

## Examples

```bash
phyloai posttree topology --matrix matrix.fa --candidate-trees candidates.trees --model-expr LG+F+R4
phyloai posttree topology --matrix matrix.fa --candidate-trees candidates.trees --partitions matrix.best_model.nex
```

## Inputs

| Input | Description |
|-------|-------------|
| `--matrix` | Single supermatrix alignment (FASTA, PHYLIP, NEXUS). Maps to IQ-TREE `-s`. |
| `--candidate-trees` | One tree-list file (one NEWICK tree per line) or multiple individual NEWICK files separated by commas (e.g. `h1.nwk,h2.nwk`). Multiple files are merged in order by PhyloAI into `candidate.trees`. Maps to IQ-TREE `-z`. |
| `--input-format` | PhyloAI-side matrix format hint (`auto|fasta|phylip-relaxed|nexus`, default `auto`). Not passed to IQ-TREE. |

## Model Source

Provide exactly one model source. PhyloAI does **not** re-run ModelFinder -- use
`phyloai tree ml iqtree` for model selection.

| Option | Description |
|--------|-------------|
| `--model-expr` | Complete IQ-TREE `-m` expression. Examples: `LG+F+R4`, `C20+F+R4`, `LG+C20+F+R4`, `custom.exchangeabilities+R4`. |
| `--partitions` | Previously optimized partition file (e.g., `.best_model.nex` from IQ-TREE). Maps to IQ-TREE `-p`. |

`--guide-tree` is used with PMSF models (e.g., `LG+C20+F+R4`). Maps to IQ-TREE `-ft`.

## Default Tests

PhyloAI generates the standard topology-test flags:

```
-n 0 -zb <replicates> -zw -au
```

| Test | Description |
|------|-------------|
| bp-RELL | Bootstrap proportion (RELL) |
| KH | Kishino-Hasegawa test |
| SH | Shimodaira-Hasegawa test |
| WKH | Weighted KH test |
| WSH | Weighted SH test |
| c-ELW | Expected likelihood weight |
| AU | Approximately unbiased test |

KH, SH, WKH, WSH, and AU are **p-values**. Trees with p < 0.05 are rejected
by that test. bp-RELL and c-ELW are **weights**, not p-values.
Recommended: AU, WSH, WKH.

## Advanced IQ-TREE Args

| Flag | Description |
|------|-------------|
| `--tool-args` | Additional IQ-TREE strategy parameters. **Blocked flags:** `-s` (matrix), `-z` (candidate trees). Shell I/O redirects (`<`, `>`, `|`) rejected. |
| `--iqtree-path` | Explicit path to `iqtree3` executable. |
| `--prefix` | IQ-TREE output prefix (default: matrix file stem). |

PhyloAI-built flags are suppressed when the same flag appears in `--tool-args`
(suppress-if-present). Overrideable flags: `-m`, `-p`, `-ft`, `-n`, `-zb`,
`-zw`, `-au`, `-T`, `--prefix`.

IQ-TREE stdout is streamed to the terminal for progress visibility. Results
are parsed from the `.iqtree` report and displayed as a formatted table.

## Outputs

IQ-TREE native files:
- `<prefix>.iqtree` -- full IQ-TREE report with topology test table
- `<prefix>.log` -- IQ-TREE log
- `<prefix>.treels.trees` -- IQ-TREE optimized candidate trees (suffix may vary)

PhyloAI files:
- `result.json` -- structured result with parsed test table
- `candidate.trees` -- merged tree file (only when multiple individual files were provided)

## Exit Codes

| Code | Meaning |
|------|---------|
| 0 | Success |
| 1 | User input error (missing files, invalid parameters, output conflict) |
| 2 | IQ-TREE execution failed |
| 3 | IQ-TREE executable not found |

## Warnings and Errors

- `--overwrite` and `--resume` are mutually exclusive.
- A non-empty output directory requires `--overwrite` or `--resume`.
- Exactly one of `--model-expr` and `--partitions` is required.

## Notes

- This command is single-matrix only (no batch mode).
- `--replicates` defaults to 10000. Very large values can make RELL resampling slow.
- IQ-TREE native resume (via `.ckp.gz`) is supported with `--resume`.
- All file paths (matrix, partitions, guide-tree) are resolved to absolute paths before IQ-TREE invocation.
