# phyloai posttree simulate adequacy

[English](posttree-simulate-adequacy.md) | [中文](posttree-simulate-adequacy.zh.md)

## Purpose

Compares PPA-DIV, PPA-CONV, PPA-VAR, and PPA-COMP statistics from an observed MSA with an empirical null distribution from simulated MSAs. The command is local-only and requires no external executable.

## Usage

```bash
phyloai posttree simulate adequacy \
  --original-msa matrix.fa \
  --simulated-dir runs/sim/MSAs \
  --threads 4 \
  --table-format csv \
  --output-dir runs/adequacy
```

## Parameters

| Option | Description |
|--------|-------------|
| `--original-msa` | Observed alignment. Required. FASTA, PHYLIP-relaxed, PHYLIP-PAML, and NEXUS are auto-detected; taxon IDs must be unique. |
| `--simulated-dir` | Directory of simulated alignments. Required. Each non-empty regular file is independently auto-detected and must have the same taxon set and alignment length as the observed MSA. |
| `--seq-type` | `AA`, `NT`, or `auto` (default; case-insensitive). `auto` resolves the type from the observed MSA before workers process simulations. |
| `--threads` | Parallel workers for simulated-MSA statistics (default `4`). |
| `--table-format` | Delimiter and suffix for all three output tables: `csv` (default) or `tsv`. |
| `-o, --output-dir` | Output directory (default `runs/posttree/simulate/adequacy`). |
| `--overwrite` | Delete and recreate a non-empty output directory. Mutually exclusive with `--resume`. |
| `--resume` | Resume from `checkpoint.json`. The observed MSA and required resume parameters must match the original run; changed simulated files are recomputed. |
| `--dry-run` | Validate the observed MSA, resolve sequence type, and count simulated files without writing output. |
| `-q, --quiet` | Suppress terminal output except errors. |

## Examples

```bash
# All-amino-acid data
phyloai posttree simulate adequacy --original-msa concat.aa.fa \
  --simulated-dir runs/sim/MSAs --seq-type AA --threads 4

# Resume a tab-delimited run
phyloai posttree simulate adequacy --original-msa concat.aa.fa \
  --simulated-dir runs/sim/MSAs --table-format tsv -o runs/adequacy --resume
```

## Outputs

```
runs/posttree/simulate/adequacy/
├── adequacy_summary.csv      # scalar PPA-DIV/CONV/VAR/COMP results
├── adequacy_taxon_comp.csv   # per-taxon PPA-COMP results
├── per_simulation_stats.csv  # raw scalar statistics per valid simulation
├── checkpoint.json           # resumable per-simulation statistics
└── result.json               # machine-readable result
```

With `--table-format tsv`, the three tables use `.tsv` suffixes and tab delimiters. `result.json` records their resolved paths under `data.output_files`.

`adequacy_summary` reports observed value, simulated mean and population SD, empirical 95% interval, z-score, posterior predictive p-value, and the number of valid simulations for `div`, `siteconvprob`, `sitecomp`, `comp_max`, and `comp_mean`. `adequacy_taxon_comp` provides the same comparison for each observed taxon.

## Warnings & Errors

| Condition | Behavior |
|-----------|----------|
| Fewer than 10 valid simulations | Hard error after processing; no adequacy summaries are produced. |
| Duplicate taxon ID, unequal sequence lengths, all-missing taxon, or no informative sites in observed MSA | Hard input error. |
| Simulated MSA has duplicate IDs, taxon/length mismatch, all-missing taxon, or cannot be parsed | File is skipped, recorded as failed, and included in warnings. |
| Output directory is non-empty without `--overwrite` or `--resume` | Preflight error; no output is changed. |
| `--resume` checkpoint is missing, incompatible, lacks the original-MSA fingerprint, or observed MSA changed | Preflight error; use `--overwrite` for a new run. |

## Exit Codes

| Code | Meaning |
|------|---------|
| `0` | Success. |
| `1` | User input, validation, output-conflict, or resume error. |

## Notes

- `div` pp is `P(sim <= obs)`; all other pp values are `P(sim > obs)`. A low pp (`< 0.05`) or `|z| > 2` indicates potential model inadequacy.
- Observed and simulated MSAs may use mixed supported formats.
- Apply `phyloai posttree simulate alisim transfergaps` before adequacy when the observed MSA has substantial missing data and simulations are gap-free.
- All statistics exclude non-standard characters: `ACDEFGHIKLMNPQRSTVWY` for AA and `ACGT` for NT. Zero-SD simulated distributions report `pp=null` in JSON and an empty pp table cell.
- This command compares one observed MSA with a replicate distribution; it is
  not an adequacy test run independently on each replicate. For strict
  `readpb --mode ppred` versus AliSim plug-in simulation choices, see the
  [systematic-error workflow reference](../../skills/phyloai-workflow/references/syserror-workflow.md).
