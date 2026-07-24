# Demo Data

Demo paths are resolved from the installed package, not hardcoded.

## End-to-End Dataset

| Path | Content |
|---|---|
| `demo_data/end_to_end/raw/faa/` | 20 protein (AA) FASTA files, 6 taxa each |
| `demo_data/end_to_end/raw/fna/` | 20 nucleotide (NT) FASTA files, 6 taxa each, matched to `faa/` by gene name |

Full pipeline: raw → convert → align → trim → concat → ML tree → report. Runs in 2–5 minutes on a modern laptop.

## Per-Step Entry Points

Directories under `demo_data/per_step/` contain pre-computed outputs so you can jump into any step without running the full pipeline.

| Directory | Command | Notes |
|---|---|---|
| `aligned/` | `pretree align` result | Pre-aligned AA sequences in `seqs/faa/`. Use `--seq-dir` pointing here for `pretree stats` or `pretree trim`. |
| `trimmed/` | `pretree trim` result | Pre-trimmed AA sequences in `seqs/faa/`. Use for `pretree metrics`, `pretree concat`. |
| `concat/` | `pretree concat` result | Supermatrices in `faa/` (protein) and `fna/` (nucleotide). Use `--matrix concat/faa/matrix.fa` for tree inference, topology tests, dating hessian. |
| `gen_trees/` | `tree ml fasttree` result | Pre-inferred gene trees in `trees/`. Use `--tree-dir` for `pretree metrics --tree-dir`, `tree msc`, `tree cf`. |
| `topology_test/` | Candidate trees for `posttree topology` | Contains `candidate.trees` (one NEWICK tree per line). Pair with `--matrix concat/faa/matrix.fa` from the `concat/` entry point. |
| `dating/` | Inputs for `posttree dating` | Contains `input.tre` (rooted tree with MCMCtree calibrations, root age constraint in 100 Mya units) and `mcmctree.ctl` (sample control file). Pair with `--matrix concat/faa/matrix.fa` from `concat/` for `dating hessian`. Use `--ctl dating/mcmctree.ctl` for `dating mcmc`. |

Demo runs must write outputs to a user run directory, never into the demo data directory.
