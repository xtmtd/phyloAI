# phyloai tree ml fasttree

[English](tree-ml-fasttree.md) | [中文](tree-ml-fasttree.zh.md)


## Purpose

Infer maximum-likelihood phylogenetic trees using FastTree.

## Usage

```bash
# Batch gene trees from MSA directory (parallel)
phyloai tree ml fasttree --msa-dir ./trimmed/seqs \
    --seq-type AA --model lg --mode normal --boot 1000 \
    --cat 20 --gamma --threads 8 -o runs/tree/ml/fasttree

# Single supermatrix tree
phyloai tree ml fasttree --matrix ./concat/matrix.fa \
    --seq-type NT --model gtr --mode slow --boot 1000 \
    -o runs/tree/ml/fasttree

# Disable bootstrap (no node support)
phyloai tree ml fasttree --msa-dir ./trimmed --boot 0

# Fast mode, JTT model (AA default)
phyloai tree ml fasttree --msa-dir ./trimmed --mode fastest --model jtt
```

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--msa-dir` | — | Directory of MSA files. Mutually exclusive with `--matrix`. |
| `--matrix` | — | Single concatenated matrix file. Mutually exclusive with `--msa-dir`. |
| `--seq-type` | auto | AA, NT, or auto (detect from input). |
| `--model` | lg (AA) / gtr (NT) | Substitution model. AA: jtt, lg, wag. NT: jc, gtr. |
| `--mode` | normal | Speed/accuracy: normal, fastest, slow. |
| `--boot` | 1000 | Bootstrap replicates. 0 = no support (-nosupport). |
| `--cat` | 20 | Number of rate categories. |
| `--gamma` | on | Gamma-distributed rate heterogeneity. Always on by default; use --tool-args to disable. |
| `--output-dir` / `-o` | runs/tree/ml/fasttree | Output directory. |
| `--threads` / `-t` | 4 | Parallel workers (--msa-dir only). |
| `--fasttree-path` | — | Explicit path to FastTree. |
| `--tool-args` | — | Extra strategy flags for FastTree. |
| `--overwrite` | — | Overwrite existing output dir. |
| `--resume` | — | Resume from checkpoint (--msa-dir only). |
| `--dry-run` | — | Print commands without executing. |
| `--quiet` / `-q` | — | Suppress output except errors. |

## Outputs

- `result.json`: structured results (trees, failed, skipped)
- `trees/`: Newick tree files (one per input, --msa-dir mode)
- `logs/`: per-gene FastTree logs (--msa-dir mode)
- `checkpoint.json`: resume state (--msa-dir mode)
- Single `.tre` file (--matrix mode)

## Supported Formats

FastTree natively reads FASTA (.fa, .fas, .fasta, .faa, .fna) and phylip-relaxed (.phy, .phylip).

NEXUS files (.nex, .nxs, .nexus) are not supported. Convert them first:
```bash
phyloai pretree convert --input data.nex --to fasta
```

## Warnings & Errors

| Condition | Behavior |
|-----------|----------|
| `--msa-dir` and `--matrix` both or neither provided | Error: exactly one is required |
| `--overwrite` and `--resume` together | Error: mutually exclusive |
| `--resume` in `--matrix` mode | Error: resume only in batch mode |
| `--threads` with `--matrix` mode | Warning: `--threads` has no effect in single mode |
| `--msa-dir` does not exist | Error: directory not found |
| No valid inputs in `--msa-dir` | Error: no valid input files |

## Notes

- `--threads` only controls parallel gene tree inference in `--msa-dir` mode. FastTree itself is single-threaded.
- `--resume` is only available in `--msa-dir` batch mode.
- Model default: LG for amino acid (AA), GTR for nucleotide (NT).
