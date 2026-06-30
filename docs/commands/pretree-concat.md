# phyloai pretree concat

[English](pretree-concat.md) | [中文](pretree-concat.zh.md)


## Purpose

`phyloai pretree concat` concatenates multiple MSA files into a supermatrix
for downstream phylogenetic inference.  Supports occupancy filtering, character
recoding, codon-specific variants (translation, 3rd-codon exclusion), outgroup
reordering, and multi-format output.

Inputs must be aligned MSA files (FASTA or Phylip).  Run `phyloai pretree
align` first, and use `phyloai pretree trim` beforehand if you need to remove
poorly aligned columns.

## Usage

Minimal:
```bash
phyloai pretree concat --msa-dir ./aligned
```

With recoding and CODON variants:
```bash
phyloai pretree concat \
  --msa-dir ./runs/pretree/trim/seqs/fna \
  --seq-type CODON \
  --recoding RY-nucleotide \
  --translate-codon \
  --exclude-codon3 \
  --outgroup Homo_sapiens \
  --taxa-occupancy 0.7 \
  --to fasta \
  --output-dir ./runs/pretree/concat
```

Gene-jackknife pseudoreplicates from an existing matrix:
```bash
phyloai pretree concat jackknife \
  --matrix runs/pretree/concat/matrix.fa \
  --partitions runs/pretree/concat/matrix.partitions \
  --replicates 100 \
  --target-length 50000 \
  --to fasta \
  --table-format csv \
  --seed 42 \
  -o runs/pretree/concat/jackknife
```

## Parameters

| Parameter | Default | Notes |
|-----------|---------|-------|
| `--msa-dir` | required | Directory of input MSA files |
| `--output-dir` / `-o` | `runs/pretree/concat` | Output directory |
| `--prefix` | `matrix` | Prefix for output filenames |
| `--seq-type` | `auto` | `AA`, `NT`, `CODON`, or `auto`; `CODON` never auto-detected |
| `--taxa-occupancy` | `0.5` | Min taxon ratio for MSA inclusion (0.0–1.0) |
| `--recoding` | — | `RY-nucleotide` (NT only; A/G→R, C/T/U→Y), `Dayhoff-6/9/12/15/18`, `SandR-6`, `KGB-6` (AA only) |
| `--outgroup` | — | Single taxon name to move to first position |
| `--to` | `fasta` | Output format: `fasta`, `phylip-relaxed`, `phylip-paml`, `nexus` |
| `--translate-codon` | off | Also produce CDS→AA translated matrix (CODON only) |
| `--exclude-codon3` | off | Also produce codon1+2 matrix (CODON only) |
| `--dry-run` | off | Validate inputs and report planned actions without writing, deleting, or replacing files |
| `--quiet` / `-q` | off | Suppress terminal output except errors |
| `--overwrite` | off | Delete and recreate non-empty output directory |

## Inputs

Scans `--msa-dir` one level deep for `.fa`, `.fas`, `.fasta`, `.faa`, `.fna`,
`.phy`.  Subdirectories, empty files, and unrecognized extensions are skipped.

Sequences are normalized per-gene via `core/sequence_normalization.py` before
concatenation.  Missing taxa in a given gene are filled with `?`.

## Outputs

```
runs/pretree/concat/
├── matrix.fa                   # (or .phy/.nex)
├── matrix.partitions           # RAxML-style partition file
├── matrix.recoded.fa           # if --recoding
├── matrix.recoded.partitions   # if --recoding
├── matrix.translated.fa        # if --translate-codon
├── matrix.translated.partitions # if --translate-codon
├── matrix.cds12.fa             # if --exclude-codon3
├── matrix.cds12.partitions     # if --exclude-codon3
├── dropped_alignments.csv      # if any MSA dropped
├── result.json
```

## Gene-Jackknife Pseudoreplicates

`phyloai pretree concat jackknife` creates pseudoreplicate matrices from an
existing concatenated matrix and partition file. It does not infer trees.

The main use case is making very large supermatrices cheaper to analyze by
sampling many smaller matrices with preserved gene boundaries. This is useful
not only for Bayesian workflows, but also for other expensive downstream runs
such as high-memory heterogeneous ML analyses.

Outputs are written as one directory per replicate:

```text
jackknife/
├── rep001/
│   ├── rep001.fa
│   └── rep001.partitions
├── rep002/
├── jackknife_summary.csv
└── result.json
```

The selected loci for each replicate are recorded in the corresponding
`repXXX.partitions` file and in `result.json`.

### Variants

| Variant | Condition | `seq_type` |
|---------|-----------|-----------|
| Original | Always | Resolved from input |
| Recoded | `--recoding` | `other` |
| Translated | `--translate-codon` | `AA` |
| Codon1+2 | `--exclude-codon3` | `NT` |

### Partition files

Each matrix is accompanied by a `.partitions` file in RAxML-style format,
describing gene boundaries in the supermatrix for partitioned analysis with
tools like IQ-TREE (`-p` option).

Each line in the file has the form:

```
TYPE, gene_name = start-end
```

**Prefix rule:**

| Matrix variant | `TYPE` |
|---|---|
| Original (NT / CODON) | `DNA` |
| Original (AA) | `LG` |
| Recoded (any) | `AUTO` |
| Translated (CODON→AA) | `LG` |
| Codon1+2 (CODON→NT) | `DNA` |

Gene names use the input file basename (without extension).  Positions are
1-indexed inclusive ranges.  For translated/cds12 variants, positions are
recomputed to match variant matrix lengths.  Not written under `--dry-run`.

Example:
```
DNA, COI = 1-654
DNA, 16S = 655-1203
```

### result.json

Generated and planned variant outputs are recorded as full paths in both
`key_results.variants_produced` and `data.variants[].path`.  All PhyloAI-authored
FASTA-family outputs from this command wrap sequence lines at 60 characters.

## Screen Display (Rich)

Three panels:

1. **Overview** — prefix, to_format, n_taxa, n_msa_* counts,
   taxon_occupancy_threshold, recoding, outgroup, variant files produced.

2. **Character Summary** — per-variant table: seq_type, total_length,
   gap_ratio, ambiguous_ratio, gap_ambiguous_ratio, standard_ratio.  Recoded
   variant (`seq_type = "other"`) shows meaningful gap_ratio and `—` for
   uncertain metrics.

3. **Site Patterns** — per-variant table: alignment_length,
   distinct_patterns (count + ratio), constant_sites, parsimony_informative,
   singleton_sites.  Ratios to 4 decimal places.  Distinct-pattern counting
   collapses all non-standard characters into the gap symbol, matching
   IQ-TREE's convention.

## Examples

```bash
# Basic NT concatenation
phyloai pretree concat --msa-dir ./aligned_nt --seq-type NT --to fasta

# CODON with all variants
phyloai pretree concat --msa-dir ./aligned_codon --seq-type CODON \
  --translate-codon --exclude-codon3

# Recoding + outgroup
phyloai pretree concat --msa-dir ./aligned_aa --seq-type AA \
  --recoding Dayhoff-6 --outgroup Sp_A

# Strict occupancy, dry-run first
phyloai pretree concat --msa-dir ./aligned --taxa-occupancy 1.0 --dry-run
```

## Warnings and Errors

| Condition | Behaviour |
|-----------|-----------|
| Missing `--msa-dir` or no MSA files found | Exit 1 |
| Non-empty output directory without `--overwrite` | Exit 1 |
| `--translate-codon` / `--exclude-codon3` on non-CODON seq_type | Exit 1 |
| AA-only recoding scheme on NT input (or vice versa) | Exit 1 |
| `--outgroup` taxon not found in matrix | Exit 1 |
| No MSAs pass `--taxa-occupancy` filter | Exit 1 |
| `--taxa-occupancy` outside 0.0–1.0 | Exit 1 |

## Notes

- Two-pass approach for memory efficiency: Pass 1 is a header-only scan
  (FASTA) to collect taxon sets and filter by occupancy; dropped files are
  never fully read.  Pass 2 is a streaming concatenation.
- Translation is per-gene (before concatenation) to preserve codon frame
  across gene boundaries.
- Statistics for the original variant are computed in-memory before format
  conversion to avoid name truncation issues with Phylip-PAML output.
- `result.json` includes `variant_stats` with per-variant character summary
  and site patterns.
- `.partitions` files are generated alongside each matrix for partitioned
  phylogenetic analysis (e.g., `iqtree -s matrix.fa -p matrix.partitions`).
- `--dry-run --overwrite` still leaves any existing output directory untouched;
  `--overwrite` only removes files during a real run.
