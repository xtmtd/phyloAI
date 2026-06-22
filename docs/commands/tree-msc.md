# phyloai tree msc

Multispecies coalescent species tree inference with [wASTRAL](https://github.com/chaoszhang/ASTER) (ASTER).

## Purpose

`phyloai tree msc` consumes gene trees and produces a species tree with local posterior probability branch support using wASTRAL. wASTRAL is a re-implementation of ASTRAL for species tree inference under the multispecies coalescent model.

wASTRAL computation is one-shot -- there is no `--resume` support.

## Usage

```bash
# Single gene tree file input
phyloai tree msc --tree gene_trees.trees -o runs/tree/msc

# Directory of gene tree files (auto-merged)
phyloai tree msc --tree-dir ./genetrees/

# Traditional unweighted Astral with exhaustive search
phyloai tree msc --tree-dir ./genetrees/ --mode 4 -R

# Bootstrap input support with custom range
phyloai tree msc --tree-dir ./genetrees/ \
    --mode 1 --boot 2 -R \
    --tree-boot-type bootstrap --tree-boot-min 10 --tree-boot-max 95 \
    -t 8 -o runs/tree/msc

# Override via --tool-args
phyloai tree msc --tree input.trees --tool-args "-r 32 -s 32"

# Root with outgroup species
phyloai tree msc --tree-dir ./genetrees/ --outgroup Oryza_sativa
```

## Inputs

| Option | Description |
|--------|-------------|
| `--tree` | Single gene tree file (newick, one tree per line). Mutually exclusive with `--tree-dir`. |
| `--tree-dir` | Directory of gene tree files. Scanned for `.nwk`, `.tre`, `.tree`, `.nw`, `.trees`, `.newick` extensions, merged into one input. Mutually exclusive with `--tree`. |

## Parameters

| Option | Default | Description |
|--------|---------|-------------|
| `--mode` | 1 | 1=hybrid, 2=branch support weighting, 3=branch length weighting, 4=traditional unweighted |
| `--boot` | 1 | wastral -u. 0=topology only, 1=local posterior probability, 2=quartet+local-PP, 3=2 + freqQuad.csv |
| `-R` / `--extra-rounds` | off | Enable exhaustive search (wastral -R). |
| `--tree-boot-type` | auto | Gene tree branch support type preset: `auto` (detect), `likelihood` (wastral -L/--lrt), `abayes` (wastral -B/--bayes), `bootstrap` (wastral -S/--bootstrap). Sets preset -x/-n values. |
| `--tree-boot-min` | -- | Min support value (wastral -n). Overrides preset default. |
| `--tree-boot-max` | -- | Max support value (wastral -x). Overrides preset default. |
| `--outgroup` | -- | Outgroup species name for rooting (wastral --root). |

## Outputs

```
runs/tree/msc/
├── result.json            # PhyloAI structured results (stderr inline in data.tool_stderr)
├── wastral.tre            # Species tree output (newick)
├── merged.trees           # Merged input (--tree-dir mode only)
└── freqQuad.csv           # Quartet frequency data (--boot 3 only)
```

## Exit Codes

| Code | Meaning |
|------|---------|
| 0 | Success |
| 1 | User input error |
| 2 | wastral execution failed |
| 3 | wastral not found |

## Notes

- wASTRAL must be installed and on PATH (or use `--wastral-path`).
- wASTRAL stderr is inlined in `result.json` as `data.tool_stderr`. No separate log file is written.
- `--boot 2` computes quartet support + local-PP and embeds values in the output tree; no separate data file is written. Use `--boot 3` for `freqQuad.csv`.
- `--tree-dir` mode merges all valid gene tree files into one input file saved as `merged.trees`.
- `--tool-args` passes extra flags verbatim to wastral. `-i` and `-o` are blocked. Strategy flags override phyloAI defaults.
- `--outgroup` specifies a single species name to root the tree (wastral `--root`).
