# phyloai posttree simulate alisim

[English](posttree-simulate-alisim.md) | [中文](posttree-simulate-alisim.zh.md)

## Purpose

Three-command IQ-TREE3 AliSim workflow that **preserves the empirical
properties of a real dataset** instead of using arbitrary model parameters:

1. **params** — extract per-locus simulation parameters (model string
   components, frequencies, proportion of invariant sites, rate heterogeneity,
   tree) from existing IQ-TREE reports into a `params.tsv` table.
2. **iqtree** — simulate alignments with AliSim, either a single IQ-TREE call
   from explicit inputs, or a resumable batch that samples rows from
   `params.tsv`.
3. **transfergaps** — re-introduce the original per-taxon gap mask onto one
   simulated alignment (`--simulated-msa`) or a batch of simulated alignments
   (`--simulated-dir`). AliSim produces gap-free MSAs.

## Usage

```bash
# 1. Extract parameters from previous IQ-TREE runs
phyloai posttree simulate alisim params --iqtree-dir runs/tree/ml/iqtree --tree-dir runs/tree/ml/iqtree

# 2. Batch simulation sampling from the extracted table
phyloai posttree simulate alisim iqtree --model-params params.tsv --strategy complete --num-simulations 100

# 3. Transfer the original gap pattern onto one simulated alignment
phyloai posttree simulate alisim transfergaps --original-msa original.fa --simulated-msa sim001.fa

# 4. ...or onto every simulated alignment in a directory
phyloai posttree simulate alisim transfergaps --original-msa original.fa --simulated-dir MSAs/
```

## alisim params

Extracts one row per successfully parsed `.iqtree` report and pairs it with a
tree file by logical locus name.

### Parameters

| Option | Description |
|--------|-------------|
| `--iqtree-dir` | Directory containing `.iqtree` report files (any nesting depth, globbed as `**/*.iqtree`). Required. |
| `--tree-dir` | Directory containing tree files. Matched to `.iqtree` files by logical locus name (suffix-agnostic). Required. |

### Outputs

- `params.tsv` — columns `id, seqtype, length, subs_model, subs_rate, freq,
  prop_inv, rate_heterogeneity, rate_categories, rate_param, tree_path`.
  Multi-value columns (`subs_rate`, `freq`, `rate_param`) use IQ-TREE's comma
  delimiter, matching the `.iqtree` reports.
- `result.json` — `n_loci_parsed`, `n_loci_matched`, `n_loci_unmatched`,
  `seq_types`, and the unmatched-loci list.

Loci whose report has no matching tree are reported as unmatched and skipped.
Ambiguous tree matching (multiple tree files with the same logical locus name)
is a hard error.

## alisim iqtree

Two mutually exclusive modes.

### Single mode (one IQ-TREE call)

Provide a reference tree, model string or partition file, sequence type, and
alignment length:

```bash
phyloai posttree simulate alisim iqtree --ref-tree ref.nwk --model LG+G4 --seq-type AA --length 2000
phyloai posttree simulate alisim iqtree --ref-tree ref.nwk --model-partitions matrix.best_model.nex --seq-type AA
```

#### Parameters

| Option | Description |
|--------|-------------|
| `--ref-tree` | Reference tree (Newick). Maps to IQ-TREE `-t`. Required. |
| `--model` | IQ-TREE model string (e.g. `GTR{XXX}+F{XXX}+G4{XXX}`). Maps to `-m`. Mutually exclusive with `--model-partitions`. |
| `--model-partitions` | NEXUS partition model file. Maps to `-p`. Mutually exclusive with `--model`. When used, `--length` is inferred from the partition definitions and must be omitted. |
| `--seq-type` | `AA` or `DNA` (case-insensitive). Maps to IQ-TREE `--seqtype`. Required. |
| `--length` | Alignment length. Maps to `--length`. Required unless `--model-partitions` is used. |
| `--msa-prefix` | Output MSA file prefix (default `sim`). |
| `--num-alignments` | Number of MSAs per IQ-TREE call (default 1). Single mode only. |
| `--out-format` | Output MSA format: `fasta` (default) or `phy`. Maps to `--out-format`. |
| `--iqtree-threads` | Threads per IQ-TREE invocation (default 1). Maps to `-T`. |
| `--seed` | Random seed. Maps to `--seed`. |

#### Outputs

- `MSAs/<prefix>.*` — simulated MSAs (extension follows `--out-format`).
- `logs/<prefix>.log` — captured IQ-TREE console output.
- `result.json` — `n_msas_generated`, the executed IQ-TREE command, and output
  file paths.

AliSim produces no `.iqtree` report.

### Batch mode (resumable, one AliSim call per MSA)

Provide a `--model-params` table (from `alisim params`) plus a sampling
strategy and target count:

```bash
phyloai posttree simulate alisim iqtree --model-params params.tsv --strategy pdf --num-simulations 100 --seed 42
```

#### Parameters

| Option | Description |
|--------|-------------|
| `--model-params` | TSV table from `alisim params`. Activates batch mode. Required. |
| `--strategy` | `complete` (default), `mixed`, or `pdf`. See [Sampling strategies](#sampling-strategies). |
| `--num-simulations` | Total number of MSAs to simulate. Required. |
| `--override` | Comma-separated `key=value` pairs fixed across all simulations, e.g. `length=500,prop_inv=0.1`. Applies to all strategies. Valid keys: `length`, `prop_inv`. |
| `--noise-scale` | PDF resampling noise (0.0-1.0, default 1.0): 0 = bin centers, 1 = full within-bin jitter. Requires `--strategy pdf`. |
| `--pdf-params` | Comma-separated parameters resampled via density estimation (default `length,prop_inv,rate_param`). Valid: `length`, `prop_inv`, `rate_param`. Requires `--strategy pdf`. |
| `--seed` | Master seed; each simulation gets an independent random seed drawn from a master-seeded generator. |
| `-t, --threads` | Parallel simulation tasks (default 4). |
| `--out-format` | Output MSA format for every simulation: `fasta` (default) or `phy`. |
| `--iqtree-threads` | Threads per IQ-TREE invocation (default 1). Maps to `-T`. |

#### Sampling strategies

- **complete** — each simulated alignment replicates the full parameter set of
  a single source gene model (model core, rate heterogeneity, alignment
  length, invariant-site proportion, and tree all taken together from one row).
- **mixed** — the model core group, rate heterogeneity group, alignment length,
  invariant-site proportion, and reference tree are each sampled independently
  from the empirical gene-model distribution, preserving individual parameter
  distributions and their presence/absence ratios.
- **pdf** (probability density function) — built on mixed sampling; the
  parameters in `--pdf-params` are resampled from histogram-based estimates of
  the empirical probability density (Freedman-Diaconis binning) with noise
  scale `--noise-scale`, instead of drawn from the empirical column directly.

Sampling notes:
- Only `""` counts as absent for `prop_inv`; a non-empty value such as `"0"`
  is preserved and reconstructed as `+I{0}`.
- `rate_param` is density-resampled only when the sampled rate group is `G`
  (Gamma); FreeRate (`R`) groups are always sampled empirically.
- `--override` fixes parameters for all strategies and overrides both table
  values and density resampling.
- PDF density plots (empirical vs simulated Gaussian-KDE curves, `server.R`
  palette `#2E86AB`/`#A23B72`, curves only) are written only for `--strategy
  pdf`; `complete`/`mixed` runs produce no `plots/` directory.

#### Outputs

- `MSAs/sim001.<ext>, ...` — one MSA per simulation.
- `logs/<simulation_id>.log` — captured IQ-TREE console output per simulation.
- `params_sampled.tsv` — every actual row used. The `source_id` column is
  present for the `complete` strategy only; each row's `seed` is independent.
- `plots/*_density.pdf` — pdf strategy only.
- `checkpoint.json`, `result.json` — `result.json` carries `source_loci`,
  `n_simulations_completed`, `n_simulations_failed`; for `complete`/`mixed`,
  `noise_scale` and `pdf_params` are `null`.

Resume rules:
- `--resume` resumes a batch run from `checkpoint.json`; completed simulations
  are skipped and unfinished ones are retried. Batch mode only.
- `--overwrite` and `--resume` are mutually exclusive.
- A non-empty existing output directory requires `--overwrite` (or `--resume`).
- Resume requires the same parameters as the original run.

## alisim transfergaps

Exactly one of `--simulated-msa` (single mode) or `--simulated-dir` (batch
mode) is required; they are mutually exclusive.

### Parameters

| Option | Description |
|--------|-------------|
| `--original-msa` | Single original (gapped) MSA file. Required. |
| `--simulated-msa` | Single simulated (gap-free) MSA file from `alisim iqtree`. Mutually exclusive with `--simulated-dir`. |
| `--simulated-dir` | Directory of simulated (gap-free) MSA files from `alisim iqtree` (alignment extensions only). One transferred file is written per input as `<stem>.gaps.fa`. Mutually exclusive with `--simulated-msa`. |
| `--seq-type` | `AA`, `NT`, or `auto` (default, case-insensitive). Determines the valid character set. |
| `--exclude-ambiguity` | When set, only real gap characters (`-`, `.`) are transferred; ambiguity codes are left as simulated. Default masks every character outside the standard alphabet (incl. ambiguity codes). |

### Outputs

- Single mode: `<original_stem>.gaps.fa`.
- Batch mode: `<simulated_stem>.gaps.fa` per input.
- `result.json` — `n_sequences`, `alignment_length`, `n_positions_masked`,
  `mean_positions_masked_per_taxon`, `detected_seq_type`, `n_msas`, and the
  transferred-file list under `data.output_files`.

Output is always 60-column FASTA regardless of input format.

## Common Flags

All three commands accept `-o, --output-dir` (default
`runs/posttree/simulate/alisim/<sub>`), `--overwrite`, `--dry-run`, and
`-q, --quiet`. `dry-run` validates inputs and prints the plan without writing
any files.

## Warnings & Errors

| Condition | Behavior |
|-----------|----------|
| Output directory exists and is non-empty, no `--overwrite`/`--resume` | Hard error naming the directory |
| `params`: no matching tree for a report | Warning; locus listed under `data.unmatched` and skipped |
| `params`: multiple trees match one locus | Hard error listing the ambiguous files |
| `transfergaps`: taxon-set or length mismatch, duplicate taxa, unparsable/empty input | Hard error (batch mode names the failing file) |
| `transfergaps`: simulated length != original length | Hard error explaining that AliSim `--length` must equal the original column count |
| `iqtree`: blocked I/O flag in `--tool-args` | Hard error listing the flag |
| `iqtree`: IQ-TREE non-zero exit | Exit code 2 with stderr |

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
- `--tool-args` on `alisim iqtree` passes extra IQ-TREE flags; only
  PhyloAI-managed I/O flags (`--alisim`, `-t`, `--prefix`, `--out-format`,
  `-af`, in either `--flag` or `--flag=value` form) are blocked. Other flags
  (e.g. `--seqtype`, `--length`, `--num-alignments`, `-T`) override PhyloAI
  defaults, and PhyloAI suppresses its own copy of an overridden flag so the
  final IQ-TREE command contains each flag exactly once.
- AliSim does not produce a `.iqtree` report; `logs/` holds the captured
  IQ-TREE console output.
- IQ-TREE3 (`iqtree3`) must be on `PATH` (or passed via `--iqtree-path`);
  `phyloai doctor` reports its detection status.
