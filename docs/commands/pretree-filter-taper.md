# phyloai pretree filter taper

## Purpose

`phyloai pretree filter taper` masks erroneous amino-acid or nucleotide sites within multiple sequence alignments using the TAPER error-correction tool (bundled `correction_multi.jl`, executed by Julia).

It operates in three modes: AA-only, NT-only, and AA+CDS (mask AA then project masks to codon-aligned nucleotide). Only newly introduced `X` masks are counted; original ambiguities are not.

Parallel execution is supported via `--threads`, and `--resume` enables checkpoint-based recovery from interrupted runs.

## Usage

Minimal:
```bash
phyloai pretree filter taper --msa-dir ./trimmed
```

AA+CDS with custom cutoff:
```bash
phyloai pretree filter taper \
  --msa-dir ./trimmed_aa \
  --nt-dir ./trimmed_fna \
  --cutoff 5 \
  --threads 8
```

## Parameters

| Parameter | Default | Notes |
|-----------|---------|-------|
| `--msa-dir` | required | Directory of input MSA files (any suffix) |
| `--output-dir` / `-o` | `runs/pretree/filter/taper` | Output directory |
| `--seq-type` | `auto` | `AA`, `NT`, or `auto`; `auto` detects from first file |
| `--nt-dir` | — | AA+CDS mode: codon-aligned NT MSA directory |
| `--cutoff` | 3 | TAPER `-c` error-correction cutoff (1-10, lower = more aggressive) |
| `--taper-path` | — | Explicit path to `correction_multi.jl`; uses bundled copy by default |
| `--julia-path` | — | Explicit Julia executable path; resolved via PATH by default |
| `--tool-args` | — | Additional TAPER flags; PhyloAI manages `-m,-a,-c,-l`, input path, output redirection |
| `--threads` / `-t` | 4 | Parallel worker processes (one locus per worker) |
| `--table-format` | `csv` | `csv` or `tsv` for auxiliary tables |
| `--show-masked-sites` | off | Include per-taxon masked-site counts in `filter_decisions.csv` |
| `--resume` | off | Resume from `checkpoint.json`; parameters must match |
| `--overwrite` | off | Delete and recreate non-empty output directory |
| `--dry-run` | off | Print commands, create no files |
| `--quiet` / `-q` | off | Suppress terminal output except errors |

## Inputs

**AA-only:** AA MSA files in `--msa-dir`. Output masked AA to `seqs/`.

**NT-only:** NT MSA files + `--seq-type NT`. Output masked NT to `seqs/`.

**AA+CDS:** AA MSA files in `--msa-dir`, codon-aligned NT MSAs in `--nt-dir` (one per AA locus, length == 3 * AA length). Output masked AA to `seqs/faa/`, projected CDS to `seqs/fna/`.

## Outputs

```
runs/pretree/filter/taper/
├── seqs/                         (or seqs/faa/ + seqs/fna/ for AA+CDS)
├── retained_loci.csv|tsv
├── dropped_loci.csv|tsv
├── filter_decisions.csv|tsv
├── checkpoint.json               (internal; only with --resume)
├── filter.log
└── result.json
```

`filter_decisions.csv` columns: `locus`, `status`, `new_masked_sites`, `masked_taxa_count`. When `--show-masked-sites` is set, an additional `masked_taxa_detail` column lists `taxon:count` per locus.

Terminal summary includes input/retained/dropped counts, masked loci/taxa/sites, and retained MSA statistics (total length, mean/min/max length, mean taxa count).

## Examples

```bash
# Default AA masking
phyloai pretree filter taper --msa-dir ./trimmed

# Aggressive masking with 8 workers
phyloai pretree filter taper --msa-dir ./trimmed --cutoff 2 --threads 8

# NT-only masking
phyloai pretree filter taper --msa-dir ./trimmed_nt --seq-type NT

# AA+CDS with per-taxon detail
phyloai pretree filter taper --msa-dir ./trimmed_aa --nt-dir ./trimmed_fna --show-masked-sites
```

## Warnings and Errors

| Condition | Behaviour |
|-----------|-----------|
| `--nt-dir` with `--seq-type NT` | Exit 1 |
| `--threads` < 1 | Exit 1 |
| `--resume` + `--overwrite` | Exit 1 |
| Julia not found | Exit 3 |
| Non-empty output directory without `--overwrite` or `--resume` | Exit 1 |
| No valid input files | Exit 1 |
| TAPER exits non-zero | Locus skipped with stderr-derived reason |
| TAPER output file missing or validation fails | Locus skipped |
| All loci skipped | Exit 2 |

## Notes

- `--tool-args` is strategy-only. Do not pass `-m`, `-a`, `-c`, `-l`, input path, or output redirection. If any of these appear, the command exits with an error.
- Only newly introduced `X` masks are counted for the `new_masked_sites` and `masked_taxa_count` fields; original ambiguity characters in the input MSA are not counted.
- Julia version is detected via `julia -v` and recorded in `result.json` and `filter.log` under `tool_versions`.
- Progress bar: shows total locus count (or `N total` label under `--resume`) with per-locus advancement.
