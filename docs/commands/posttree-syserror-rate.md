# phyloai posttree syserror rate

[English](posttree-syserror-rate.md) | [中文](posttree-syserror-rate.zh.md)

## Purpose

Rank empirical per-site rates from IQ-TREE3 or PhyloBayes and optionally
extract the slowest or fastest fractions of an existing alignment. This is a
site-rate ranking/extraction sensitivity utility: compare downstream analyses
across fractions; it does not prove a topology, correct systematic error, or
automatically choose a fraction.

This is a local pure-Python command. It invokes no external executable, so no
`phyloai doctor` check or resume support is needed.

## Usage

```bash
phyloai posttree syserror rate \
  (--iqtree-rate matrix.rate | --pb-rate chain.meansiterates) \
  [--matrix raw.fa] [--subset slow|fast] [--fraction 0.25,0.5,0.75] [options]
```

## Inputs

| Option | Required | Default | Description |
|---|---|---|---|
| `--iqtree-rate` | XOR | -- | IQ-TREE3 `--rate` table with `Site` and `Rate` columns. |
| `--pb-rate` | XOR | -- | PhyloBayes `readpb -r` headerless `<site> <rate>` table. |
| `--matrix` | no | -- | Original MSA; enables subset extraction. Supports FASTA, relaxed PHYLIP, PAML PHYLIP, and NEXUS. |
| `--subset` | with matrix | `slow` | Extract the `slow` or `fast` ranked sites. Invalid without `--matrix`. |
| `--fraction` | with matrix | -- | One or more comma-separated fractions in `(0, 1]`. Required with `--matrix`; invalid without it. |
| `-o`, `--output-dir` | no | `runs/posttree/syserror/rate` | Output directory. |
| `--overwrite` | no | false | Delete and recreate a non-empty output directory. |
| `--dry-run` | no | false | Validate inputs without writing files. |
| `-q`, `--quiet` | no | false | Suppress terminal output except errors. |

Exactly one rate source is required. A matrix must have `--fraction`; selection
defaults to `slow` only when a matrix is supplied.

## Rate Inputs And Indexing

IQ-TREE site identifiers must be strict consecutive 1-based indices. PhyloBayes
identifiers must be strict consecutive 0-based indices and are normalized by
adding one. Both sources therefore produce exactly `1..N` normalized sites.
Empty, malformed, duplicate, non-integer, non-finite, negative, or
non-consecutive rates are rejected.

`rates.csv` is deterministic and always uses 1-based normalized sites, sorted
slowest to fastest by `(rate, site)`:

```csv
site,rate
21,0.19145
```

## Outputs

Without `--matrix`:

```text
runs/posttree/syserror/rate/
├── rates.csv
└── result.json
```

With `--matrix` and `--subset slow --fraction 0.25,0.5`:

```text
runs/posttree/syserror/rate/
├── rates.csv
├── slow25/
│   ├── positions.txt
│   └── matrix.fa
├── slow50/
│   ├── positions.txt
│   └── matrix.fa
└── result.json
```

Each fraction retains `ceil(N * fraction)` sites. Ties at the boundary are not
expanded. `positions.txt` contains one 1-based original site position per line
in original alignment order. Generated subsets are FASTA, wrapped at 60
characters per sequence line. `result.json` records inputs, resolved settings,
summary statistics, subset counts, and every generated file.

## Examples

```bash
# Rank IQ-TREE rates only
phyloai posttree syserror rate --iqtree-rate matrix.rate

# Slow-site sensitivity analysis
phyloai posttree syserror rate --iqtree-rate matrix.rate --matrix raw.fa \
  --subset slow --fraction 0.25,0.5,0.75 -o runs/posttree/syserror/rate

# Fast-site extraction from PhyloBayes rates
phyloai posttree syserror rate --pb-rate chain.meansiterates --matrix raw.phy \
  --subset fast --fraction 0.1
```

## Warnings / Errors

- Exactly one of `--iqtree-rate` and `--pb-rate` is required.
- `--subset` and `--fraction` require `--matrix`; `--matrix` requires
  `--fraction`.
- Matrix length must equal the normalized rate count.
- Fractions must be valid, unique, and produce unique directory labels.
- A non-empty output directory requires `--overwrite`; no resume/checkpoint is
  available.

## Notes

- Use several defensible fractions, commonly `0.25,0.5,0.75`, as a sensitivity
  analysis rather than treating any one threshold as biologically privileged.
- Slow-site subsets reduce the contribution of rapidly evolving sites; fast-site
  subsets isolate that contribution. Neither interpretation establishes the
  true topology.
- Tree inference and comparison of extracted matrices are never automatic.
