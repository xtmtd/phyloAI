# phyloai posttree simulate alisim

[English](posttree-simulate-alisim.md) | [中文](posttree-simulate-alisim.zh.md)

Three-command workflow for IQ-TREE3 AliSim sequence simulation that
**preserves empirical properties of a real dataset** instead of using
arbitrary model parameters.

## Workflow

1. **params** — extract per-locus simulation parameters from existing IQ-TREE
   reports (model string components, frequencies, proportion of invariant
   sites, rate heterogeneity, tree) into a `params.tsv` table.
2. **iqtree** — simulate alignments with AliSim, either a single IQ-TREE call
   from explicit inputs, or a resumable batch that samples rows from the
   `params.tsv`.
3. **transfergaps** — re-introduce the original per-taxon gap mask onto a
   simulated alignment (AliSim produces gap-free MSAs).

```bash
# 1. Extract parameters from previous IQ-TREE runs
phyloai posttree simulate alisim params --iqtree-dir runs/tree/ml/iqtree --tree-dir runs/tree/ml/iqtree

# 2. Batch simulation sampling from the extracted table
phyloai posttree simulate alisim iqtree --model-params params.tsv --strategy complete --num-simulations 100

# 3. Transfer the original gap pattern onto one simulated alignment
phyloai posttree simulate alisim transfergaps --original-msa original.fa --simulated-msa sim001.fa
```

## `alisim params`

| Option | Description |
|--------|-------------|
| `--iqtree-dir` | Directory containing `.iqtree` report files (any nesting depth, globbed as `**/*.iqtree`). Required. |
| `--tree-dir` | Directory containing tree files. Matched to `.iqtree` files by logical locus name (suffix-agnostic). Required. |

Writes `params.tsv` (columns `id, seqtype, length, subs_model, subs_rate,
freq, prop_inv, rate_heterogeneity, rate_categories, rate_param, tree_path`)
plus `result.json` with `n_loci_parsed`, `n_loci_matched`, `n_loci_unmatched`,
and `seq_types`. Loci whose report has no matching tree are reported as
unmatched and skipped.

## `alisim iqtree`

Two mutually exclusive modes.

### Single mode (one IQ-TREE call)

Provide a reference tree, model string or partition file, sequence type, and
alignment length:

```bash
phyloai posttree simulate alisim iqtree --ref-tree ref.nwk --model LG+G4 --seq-type AA --length 2000
phyloai posttree simulate alisim iqtree --ref-tree ref.nwk --model-partitions matrix.best_model.nex --length 2000
```

| Option | Description |
|--------|-------------|
| `--ref-tree` | Reference tree (Newick). Maps to IQ-TREE `-t`. |
| `--model` | IQ-TREE model string (e.g. `GTR{XXX}+F{XXX}+G4{XXX}`). Maps to `-m`. Mutually exclusive with `--model-partitions`. |
| `--model-partitions` | NEXUS partition model file. Maps to `-p`. Mutually exclusive with `--model`. |
| `--seq-type` | `AA` or `DNA`. Maps to `--seqtype`. |
| `--length` | Alignment length. Maps to `--length`. |
| `--num-alignments` | Number of MSAs per IQ-TREE call (default 1). Single mode only. |
| `--msa-prefix` | Output MSA file prefix (default `sim`). |
| `--out-format` | Output MSA format: `fasta` (default) or `phy`. Maps to `--out-format`. |
| `--iqtree-threads` | Threads per IQ-TREE invocation (default 1). Maps to `-T`. |
| `--seed` | Random seed. Maps to `--seed`. |

Output: `MSAs/<prefix>.*`, `logs/<prefix>.iqtree`, `logs/<prefix>.log`,
`result.json`.

### Batch mode (resumable, one AliSim call per MSA)

Provide a `--model-params` table (from `alisim params`) plus a sampling
strategy and target count:

```bash
phyloai posttree simulate alisim iqtree --model-params params.tsv --strategy pdf --num-simulations 100 --seed 42
```

| Option | Description |
|--------|-------------|
| `--model-params` | TSV table from `alisim params`. Activates batch mode. |
| `--strategy` | `complete` (sample full rows uniformly), `mixed` (randomize model class + length), or `pdf` (histogram-density resampling). Required in batch mode. |
| `--num-simulations` | Total number of MSAs to simulate. Required in batch mode. |
| `--override` | Comma-separated `key=value` pairs fixed across all simulations, e.g. `length=500,prop_inv=0.1`. Valid keys: `length`, `prop_inv`. |
| `--noise-scale` | PDF resampling noise (0.0-1.0, default 1.0): 0 = bin centers, 1 = full within-bin jitter. Requires `--strategy pdf`. |
| `--pdf-params` | Comma-separated parameters sampled via density resampling (default `length,prop_inv,rate_param`). Valid: `length`, `prop_inv`, `rate_param`. Requires `--strategy pdf`. |
| `--seed` | Master seed; per-task seeds = master + task index. |
| `-t, --threads` | Parallel simulation tasks (default 4). |

Sampling notes:
- Complete rows are sampled as atomic units: the model core and its rate group
  stay together, and `+I` presence is decided before its value is resampled.
- Only `""` counts as absent for `prop_inv`; a non-empty value such as `"0"`
  is preserved and reconstructed as `+I{0}`.
- PDF density plots are written only for the selected PDF parameters
  (non-overridden); `length` uses Freedman-Diaconis bins, and `prop_inv` is
  clamped to `[0, 1)`.

Output: `MSAs/sim001.fa, ...`, `logs/<simulation_id>.log`, `params_sampled.tsv`
(every actual row used), `plots/*_density.pdf` (pdf strategy), `checkpoint.json`,
`result.json` with `source_loci`, `n_simulations_completed`,
`n_simulations_failed`.

Resume rules:
- `--resume` resumes a batch run from `checkpoint.json`; completed simulations
  are skipped and unfinished ones are retried. Batch mode only.
- `--overwrite` and `--resume` are mutually exclusive.
- A non-empty existing output directory requires `--overwrite` (or `--resume`).
- Resume requires the same parameters as the original run.

## `alisim transfergaps`

| Option | Description |
|--------|-------------|
| `--original-msa` | Single original (gapped) MSA file. Required. |
| `--simulated-msa` | Single simulated (gap-free) MSA file from `alisim iqtree`. Required. |
| `--seq-type` | `AA`, `NT`, or `auto` (default). Determines the valid character set. |
| `--exclude-ambiguity` | When set, only real gap characters (`-`, `.`) are transferred; ambiguity codes are left as simulated. Default masks every character outside the standard alphabet (incl. ambiguity codes). |

Validations: inputs must parse and be non-empty, no duplicate taxon IDs, the
taxon sets must match, and original/simulated alignments must have equal
length. Mask positions are replaced (never inserted); output order follows the
original MSA.

Output: `<original_stem>_transferred.<ext>` (60-column FASTA) and `result.json`
with `n_sequences`, `alignment_length`, `n_positions_masked`,
`mean_positions_masked_per_taxon`, `detected_seq_type`.

## Common Flags

All three commands accept `-o, --output-dir` (default
`runs/posttree/simulate/alisim/<sub>`), `--overwrite`, `--dry-run`, and
`-q, --quiet`. `dry-run` validates inputs and prints the plan without writing
any files.

## Exit Codes

| Code | Meaning |
|------|---------|
| 0 | Success |
| 1 | User input error (missing files, invalid parameters, output conflict) |
| 2 | IQ-TREE execution failed |
| 3 | IQ-TREE executable not found / environment failure |

## Notes

- `adequacy` and `phybase` subcommands are reserved for future work and
  currently return a not-implemented message.
- `--tool-args` on `alisim iqtree` passes extra IQ-TREE flags; managed flags
  (`--alisim`, `-t`, `-m`, `-p`, `-q`, `-Q`, `--seqtype`, `--length`,
  `--out-format`, `-af`, `--num-alignments`, `-T`, `--seed`, `--prefix`) are
  blocked.
- IQ-TREE3 (`iqtree3`) must be on `PATH` (or passed via `--iqtree-path`);
  `phyloai doctor` reports its detection status.
