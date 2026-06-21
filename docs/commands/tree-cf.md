# phyloai tree cf

Concordance factor computation — quantify branch support using gene concordance (gCF), site concordance (sCF, sCFl), and quartet concordance (qCF).

## Purpose

Concordance factors measure the proportion of gene trees, sites, or quartets supporting each bipartition in a reference species tree. They provide complementary branch-support information beyond standard bootstrap.

## Usage

```bash
phyloai tree cf --cf MODE --ref-tree REF_TREE [INPUTS...] [OPTIONS]
```

## Modes

| `--cf`    | Index  | Tool    | Description | Origin |
|-----------|--------|---------|-------------|--------|
| `gcf`     | gCF    | IQ-TREE3 | Gene concordance factor | Minh et al. 2020, *Mol Biol Evol* 37(5):1530–1534 |
| `scf`     | sCF    | IQ-TREE3 | Site concordance factor (parsimony) | Minh et al. 2020, *Mol Biol Evol* 37(5):1530–1534 |
| `scfl`    | sCFl   | IQ-TREE3 | Site concordance factor (likelihood) | Mo et al. 2023, *Syst Biol* 72(3):559–574 |
| `gcf+scf` | gCF+sCF | IQ-TREE3 | Combined gCF + sCF in one run | Minh et al. 2020 |
| `qcf`     | qCF    | wASTRAL  | Quartet concordance factor | Mirarab et al. 2014, *Science* 346(6215):1250463 |

### CF Index Details

- **gCF** (gene concordance factor): Proportion of gene trees that contain a given bipartition. Assesses gene-tree heterogeneity.
- **sCF** (site concordance factor, parsimony): Proportion of alignment sites supporting a bipartition under the parsimony criterion. Computationally fast.
- **sCFl** (site concordance factor, likelihood): Proportion of alignment sites supporting a bipartition using maximum likelihood. More accurate than sCF but requires model selection.
- **gCF+sCF**: Computes both gCF and sCF in a single IQ-TREE3 invocation, saving compute time.
- **qCF** (quartet concordance factor): Proportion of quartets (from gene trees) supporting each bipartition. Computed by wASTRAL using its quartet scoring engine.

## Input Requirements by Mode

| Mode     | `--ref-tree` | `--tree`/`--tree-dir` | `--matrix` | `--model`/`--partitions` |
|----------|-------------|----------------------|-----------|-------------------------|
| `gcf`      | Required    | Required             | —         | —                       |
| `scf`      | Required    | —                    | Required  | —                       |
| `scfl`     | Required    | —                    | Required  | Optional (speedup)      |
| `gcf+scf`  | Required    | Required             | Required  | —                       |
| `qcf`      | Required    | Required             | —         | —                       |

## Options

| Option | Default | Description |
|--------|---------|-------------|
| `--cf MODE` | *required* | Concordance factor type (gcf/scf/scfl/gcf+scf/qcf) |
| `--ref-tree FILE` | *required* | Reference species tree (NEWICK) |
| `--tree FILE` | — | Single gene tree file (mutually exclusive with `--tree-dir`) |
| `--tree-dir DIR` | — | Directory of gene tree files |
| `--matrix FILE` | — | Multiple sequence alignment (required for scf/scfl/gcf+scf) |
| `--model TEXT` | — | Substitution model for scfl (e.g., `LG+F+R3`) |
| `--partitions FILE` | — | Partition file for scfl (e.g., `*.best_model.nex`) |
| `--scf-quartets N` | 100 | Number of quartets for sCF/sCFl (recommend >= 100) |
| `--prefix TEXT` | auto | Output prefix (default: gCF/sCF/sCFl/gCFsCF/qCF) |
| `-o, --output-dir DIR` | `runs/tree/cf` | Output directory |
| `-t, --threads N` | 4 | CPU threads |
| `--iqtree-path PATH` | auto | Explicit iqtree3 executable |
| `--wastral-path PATH` | auto | Explicit wastral executable |
| `--lpp` | off | Also append local posterior probability to qCF labels |
| `--overwrite` | off | Remove existing output directory |
| `--dry-run` | off | Show command without executing |
| `-q, --quiet` | off | Suppress non-error output |

## Examples

```bash
# gCF: gene trees + reference tree
phyloai tree cf --cf gcf --ref-tree species.nwk --tree-dir ./genetrees/

# gCF with single file
phyloai tree cf --cf gcf --ref-tree species.nwk --tree merged.trees

# sCF: alignment + reference tree (ideally gCF-annotated)
phyloai tree cf --cf scf --ref-tree gCF.cf.tree --matrix msa.fa

# sCFl (likelihood) with model for speedup
phyloai tree cf --cf scfl --ref-tree gCF.cf.tree --matrix msa.fa --model LG+F+R3

# sCFl with pre-computed partition model
phyloai tree cf --cf scfl --ref-tree gCF.cf.tree --matrix msa.fa \
    --partitions msa.best_model.nex

# Combined gCF + sCF
phyloai tree cf --cf gcf+scf --ref-tree species.nwk --tree-dir ./genetrees/ \
    --matrix msa.fa

# qCF via wASTRAL
phyloai tree cf --cf qcf --ref-tree species.nwk --tree merged.trees

# qCF with local posterior probability (pp1) appended
phyloai tree cf --cf qcf --ref-tree species.nwk --tree merged.trees --lpp

# Custom output prefix and threads
phyloai tree cf --cf gcf --ref-tree species.nwk --tree merged.trees \
    --prefix myCF -t 8
```

## Output Files

### IQ-TREE3 modes (gcf, scf, scfl, gcf+scf)

| File | Description |
|------|-------------|
| `<prefix>.cf.stat`  | Concordance factor statistics table |
| `<prefix>.cf.branch` | Tree with branch IDs |
| `<prefix>.cf.tree`  | Tree annotated with CF values |
| `<prefix>.cf.tree.nex` | NEXUS annotated tree for FigTree |
| `<prefix>.log`       | IQ-TREE3 log |
| `cf.log`             | PhyloAI module log (command, versions, timing) |
| `result.json`        | PhyloAI structured result |
| `merged.trees`       | Merged gene trees (if `--tree-dir` used) |

### qCF mode

| File | Description |
|------|-------------|
| `<prefix>.cf.tree` | Reference tree annotated with qCF (and optionally pp1) |
| `wastral.tre`      | Raw wASTRAL output (intermediate) |
| `wastral.log`      | wASTRAL log |
| `cf.log`           | PhyloAI module log (command, versions, timing) |
| `result.json`      | PhyloAI structured result |
| `merged.trees`     | Merged gene trees (if `--tree-dir` used) |

## qCF Output Format

qCF values are kept as raw decimals in [0,1] (not multiplied by 100). Trailing zeros
are stripped for readability (e.g., `1` not `1.0000`, `0.95` not `0.9500`). When appended
to existing support values, the format is:

- Without `--lpp`: `<support>/<q1>`  (e.g., `100/0.4221`)
- With `--lpp`: `<support>/<q1>/<pp1>`  (e.g., `100/0.4221/0.95`)

If no existing support exists, the qCF value becomes the sole label: `0.75`.

## Notes

- For best sCF/sCFl results, use a gCF-annotated tree as `--ref-tree` (e.g., run `--cf gcf` first).
- `--cf scfl` without `--model` or `--partitions` auto-computes the best-fit model — this is slow. Provide `--model` or `--partitions` for speedup.
- `--scf-quartets` should be >= 100 for reliable results. Higher values improve accuracy at the cost of runtime.
- gCF+sCF computation runs IQ-TREE3 once with both modes, saving significant time compared to two separate runs.
- qCF uses wASTRAL's calibrated quartet scoring (`-u 2 -C --mode 4`), which handles gene tree estimation error.

## References

- Minh BQ, Hahn MW, Lanfear R (2020) New methods to calculate concordance factors for phylogenomic datasets. *Molecular Biology and Evolution* **37**(5):1530–1534.
- Mo YK, Lanfear R, Hahn MW, Minh BQ (2023) Updated site concordance factors minimize effects of homoplasy and taxon sampling. *Systematic Biology* **72**(3):559–574.
- Mirarab S, Reaz R, Bayzid MS, Zimmermann T, Swenson MS, Warnow T (2014) ASTRAL: genome-scale coalescent-based species tree estimation. *Bioinformatics* **30**(17):i541–i548.
- Mirarab S, Reaz R, Bayzid MS, Zimmermann T, Swenson MS, Warnow T (2014) ASTRAL: genome-scale coalescent-based species tree estimation. *Science* **346**(6215):1250463.
- Zhang C, Rabiee M, Sayyari E, Mirarab S (2018) ASTRAL-III: polynomial time species tree reconstruction from partially resolved gene trees. *BMC Bioinformatics* **19**(Suppl 6):153.
