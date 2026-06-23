# phyloai posttree topology

## Purpose

Performs IQ-TREE tree topology tests (AU / KH / SH / WKH / WSH / c-ELW)
comparing a set of candidate trees against a supermatrix alignment. This
command tests whether alternative topologies are significantly worse than
the best-scoring candidate -- it does **not** infer new trees.

## Usage

```bash
# Homogeneous model
phyloai posttree topology \
  --matrix matrix.fa \
  --candidate-trees candidates.trees \
  --model-expr LG+F+R4

# PMSF model with guide tree
phyloai posttree topology \
  --matrix matrix.fa \
  --candidate-trees candidates.trees \
  --model-expr LG+C20+F+R4 \
  --guide-tree guide.nwk

# Previously optimized partition model
phyloai posttree topology \
  --matrix matrix.fa \
  --candidate-trees candidates.trees \
  --partitions matrix.best_model.nex

# Multiple individual tree files (merged by PhyloAI)
phyloai posttree topology \
  --matrix matrix.fa \
  --candidate-trees h1.nwk \
  --candidate-trees h2.nwk \
  --candidate-trees h3.nwk \
  --model-expr LG+F+R4

# Custom exchangeabilities + site frequencies via --tool-args
phyloai posttree topology \
  --matrix matrix.fa \
  --candidate-trees trees \
  --model-expr custom.exchangeabilities+R4 \
  --tool-args "-fs custom.sitefreq" -t 30

# Heterogeneous model
phyloai posttree topology \
  --matrix matrix.fa \
  --candidate-trees trees \
  --model-expr C20+F+R4
```

## Inputs

| Input | Description |
|-------|-------------|
| `--matrix` | Single supermatrix alignment (FASTA, PHYLIP, NEXUS, or CLUSTAL). Maps to IQ-TREE `-s`. |
| `--candidate-trees` | One tree-list file (one NEWICK tree per line) or multiple individual NEWICK files. Multiple files are merged in order by PhyloAI into `candidate.trees`. Maps to IQ-TREE `-z`. |
| `--input-format` | Optional PhyloAI-side matrix format hint (`auto` by default). Not passed to IQ-TREE. |

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
by that test. bp-RELL and c-ELW are **weights**, not p-values. The AU test is
generally considered the most reliable.

## Advanced IQ-TREE Args

| Flag | Description |
|------|-------------|
| `--tool-args` | Additional IQ-TREE strategy parameters. **Blocked flags:** `-s` (matrix), `-z` (candidate trees). Shell I/O redirects (`<`, `>`, `|`) rejected. |
| `--iqtree-path` | Explicit path to `iqtree3` executable. |
| `--prefix` | IQ-TREE output prefix (default: matrix file stem). |

PhyloAI-built flags are suppressed when the same flag appears in `--tool-args`
(suppress-if-present). Overrideable flags: `-m`, `-p`, `-ft`, `-n`, `-zb`,
`-zw`, `-au`, `-T`, `--prefix`.

## Input Format and Sequence Type

`--input-format` (`auto|fasta|phylip-relaxed|nexus|clustal`) only affects
PhyloAI's own matrix preflight validation; it is **not** passed to IQ-TREE.
Explicit IQ-TREE `--seqtype` belongs in `--tool-args` when needed
(e.g., `--tool-args "--seqtype AA"`).

## Outputs

IQ-TREE native files:
- `<prefix>.iqtree` -- full IQ-TREE report with topology test table
- `<prefix>.log` -- IQ-TREE log
- `<prefix>.treels.trees` -- IQ-TREE optimized candidate trees (suffix may vary)

PhyloAI files:
- `result.json` -- structured result with parsed test table
- `candidate.trees` -- merged tree file (only when multiple `--candidate-trees` were provided)

## Exit Codes

| Code | Meaning |
|------|---------|
| 0 | Success |
| 1 | User input error (missing files, invalid parameters, output conflict) |
| 2 | IQ-TREE execution failed |
| 3 | IQ-TREE executable not found |

## Notes

- This command is single-matrix only (no batch mode).
- `--replicates` defaults to 10000. Very large values can make RELL resampling slow.
- IQ-TREE native resume (via `.ckp.gz`) is supported with `--resume`.
