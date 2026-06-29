# phyloai tree ml iqtree

[English](tree-ml-iqtree.md) | [中文](tree-ml-iqtree.zh.md)


## Purpose

Infer maximum-likelihood phylogenetic trees using IQ-TREE3. Supports homogeneous models (with or without ModelFinder and partitions), heterogeneous AA models (mixture and PMSF), NT heterogeneous models (MIX+MF), and branch support computation (UFBoot, SH-aLRT).

## Usage

```bash
# Batch gene trees (homogeneous workflows only)
phyloai tree ml iqtree --msa-dir <dir> [OPTIONS]

# Single supermatrix (all workflows)
phyloai tree ml iqtree --matrix <file> [OPTIONS]
```

## Examples

```bash
# Batch: 20 gene trees, AA, LG model, 4 parallel jobs (default)
phyloai tree ml iqtree --msa-dir msas/ --seq-type AA

# Single matrix: fixed model, UFBoot + SH-aLRT
phyloai tree ml iqtree --matrix matrix.fa --model LG --boot 1000 --alrt 1000

# ModelFinder only (model selection, no tree)
phyloai tree ml iqtree --matrix matrix.fa --modelfinder MF --mset LG,WAG

# ModelFinder + tree + partitions + merge
phyloai tree ml iqtree --matrix matrix.fa --modelfinder MFP --partitions parts.nex

# Disable branch support
phyloai tree ml iqtree --matrix matrix.fa --model LG --boot 0

# AA mixture model (direct)
phyloai tree ml iqtree --matrix matrix.fa --model C20

# PMSF AA mixture
phyloai tree ml iqtree --matrix matrix.fa --model C20 --guide-tree guide.nwk

# NT heterogeneous
phyloai tree ml iqtree --matrix matrix.fa --seq-type NT --model MIX+MF
```

## Parameters

### Input (mutually exclusive)

| Flag | Type | Description |
|------|------|-------------|
| `--msa-dir` | Path | Directory of MSA files for batch gene trees |
| `--matrix` | Path | Single concatenated matrix for supermatrix inference |

Supported formats: `.fa`, `.fas`, `.fasta`, `.faa`, `.fna`, `.phy`, `.phylip`, `.nex`, `.nxs`, `.nexus`, `.aln`.

### Data Type

| Flag | Type | Default | Description |
|------|------|---------|-------------|
| `--seq-type` | `AA\|NT\|auto` | `auto` | Molecule type. `NT` → `--seqtype DNA` |

### Model (when `--modelfinder none`)

| Flag | Type | Default | Description |
|------|------|---------|-------------|
| `--model` | str | `LG` (AA) / `GTR` (NT) | Substitution model |
| `--state-freq` | `+F\|+FO\|+FQ\|+FU\|none` | `+F` | State frequency type |
| `--rate-heterogeneity` | `+I\|+G4\|+I+G4\|+R4\|+I+R4\|none` | `+R4` | Rate heterogeneity |

These combine to form `-m` (e.g., `LG+F+R4`). Ignored when `--modelfinder` is `MF` or `MFP`.

### ModelFinder

| Flag | Type | Default | Description |
|------|------|---------|-------------|
| `--modelfinder` | `MF\|MFP\|none` | `none` | `MF` = model-only; `MFP` = model + tree |
| `--mset` | str | `LG,WAG` / `GTR,HKY` | Model search space restriction |
| `--msub` | `nuclear\|mitochondrial\|chloroplast\|viral` | — | AA model source (AA only) |

### Partitions (`--matrix` only)

| Flag | Type | Default | Description |
|------|------|---------|-------------|
| `--partitions` | Path | — | Partition file (NEXUS or RAxML-style) |
| `--rclusterf` | int (1–100) | `10` | Merge percentage for MF/MFP |
| `--rcluster-max` | int | — | Max merge pairs (mutually exclusive with `--rclusterf`) |

### Heterogeneous Models (`--matrix` only)

| Flag | Type | Default | Description |
|------|------|---------|-------------|
| `--pmsf-base-model` | str | `LG` | Base AA model for PMSF (C10–C60 only) |
| `--guide-tree` | Path | — | Guide tree for PMSF (NEWICK) |
| `--qmax` | int | `10` | Rate categories for MIX+MF |

### Tree Search

| Flag | Type | Default | Description |
|------|------|---------|-------------|
| `--mode` | `normal\|fast` | `normal` | `fast` → `--fast` |

### Branch Support

| Flag | Type | Default | Description |
|------|------|---------|-------------|
| `--boot` | int (>=0) | `1000` | UFBoot replicates. `0` = skip |
| `--alrt` | int (>=0) | — | SH-aLRT replicates. `0` = parametric |
| `--bnni` | flag | `False` | NNI-optimize UFBoot trees |

### Output

| Flag | Type | Default | Description |
|------|------|---------|-------------|
| `--rate` | flag | `False` | Write site rates to `.rate` |
| `--wslr` | flag | `False` | Write site log-likelihoods to `.sitelh` |
| `--constraint` | Path | — | Constraint tree (NEWICK) |
| `--outgroup` | str | — | Outgroup taxa (comma-separated) |
| `--prefix` | str | — | Output prefix (`--matrix` only) |

### Execution

| Flag | Type | Default | Description |
|------|------|---------|-------------|
| `-o`, `--output-dir` | Path | `runs/tree/ml/iqtree` | Output directory |
| `--threads`, `-t` | int\|auto | `4` / `auto` | Batch: parallel jobs. Single: IQ-TREE threads |
| `--overwrite` | flag | `False` | Remove existing output dir first |
| `--resume` | flag | `False` | Resume from checkpoint |
| `--dry-run` | flag | `False` | Print commands without executing |
| `--keep-extra` | flag | `False` | Keep extra IQ-TREE files in `logs/` |
| `-q`, `--quiet` | flag | `False` | Suppress output except errors |
| `--iqtree-path` | Path | — | Custom path to `iqtree3` |
| `--tool-args` | str | — | Extra IQ-TREE flags |

## Output

### Single mode (`--matrix`)

```
runs/tree/ml/iqtree/
├── <prefix>.iqtree
├── <prefix>.treefile
├── <prefix>.log
├── ... (other IQ-TREE outputs)
├── result.json
```

### Batch mode (`--msa-dir`)

```
runs/tree/ml/iqtree/
├── trees/
│   ├── <gene1>.treefile
│   └── <gene2>.treefile
├── logs/
│   ├── <gene1>.iqtree
│   ├── <gene1>.log
│   └── ... (extra files only with --keep-extra)
├── checkpoint.json
├── result.json
```

## Exit Codes

| Code | Meaning |
|------|---------|
| `0` | Success |
| `1` | User input error |
| `2` | All IQ-TREE runs failed |
| `3` | `iqtree3` not found |

## Warnings & Errors

| Condition | Behavior |
|-----------|----------|
| `--modelfinder MF` with `--boot`/`--alrt`/`--bnni` | Warning: branch support flags ignored in model-only mode |
| `--bnni` without `--boot > 0` | Warning: `--bnni` has no effect |
| `--prefix` in `--msa-dir` batch mode | Warning: `--prefix` ignored; gene names used |
| `--rclusterf`/`--rcluster-max` without `--partitions` | Warning: flags have no effect |
| `--qmax` without `--model MIX+MF` | Warning: `--qmax` only takes effect with MIX+MF |
| Heterogeneous model with `--msa-dir` | Error: heterogeneous workflows require `--matrix` |
| `--partitions` with `--msa-dir` | Error: `--partitions` requires `--matrix` |
| `--overwrite` and `--resume` together | Error: mutually exclusive |
| Non-empty output dir without `--overwrite` | Error: directory already exists |

## Notes

- `--boot` defaults to `1000` (enabled), matching `phyloai tree ml fasttree`. Use `--boot 0` to skip branch support.
- `--threads` defaults to `4` parallel jobs in batch mode, `auto` in single mode (IQ-TREE determines optimal thread count).
- Single `--matrix` mode streams IQ-TREE stdout to the terminal for progress visibility.
- Batch `--msa-dir` mode keeps only `.iqtree` and `.log` files in `logs/` by default. Use `--keep-extra` to preserve all IQ-TREE output files.
- When `--modelfinder` is `MF` or `MFP`, `--model`, `--state-freq`, and `--rate-heterogeneity` are not passed to IQ-TREE and are recorded as `null` in `result.json`.
- IQ-TREE input and user-provided path arguments (`--partitions`, `--guide-tree`, `--constraint`) are resolved to absolute paths internally.
- `--resume` in `--matrix` mode re-runs the IQ-TREE command; IQ-TREE natively handles checkpoint/resume via its own mechanism (`--redo`). In `--msa-dir` batch mode, PhyloAI manages checkpoint state to skip completed gene trees.
