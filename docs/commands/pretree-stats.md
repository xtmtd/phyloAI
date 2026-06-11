# phyloai pretree stats

## Purpose

`phyloai pretree stats` is a read-only QC and inspection command for sequence and alignment files. It summarizes format, sequence type, alignment status, taxon counts, length distributions, gap ratios, ambiguous-character ratios, and site-pattern statistics for aligned files.

It does not normalize, convert, align, trim, or modify input files. Use `phyloai pretree convert` first when raw inputs may contain mixed formats or characters that should be standardized before inspection.

## Usage

```bash
phyloai pretree stats [OPTIONS]
```

Exactly one of `--seq` or `--seq-dir` is required.

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--seq-dir DIRECTORY` | none | Directory mode. Scan supported sequence/alignment files in one folder and compute a dataset-level summary. |
| `--seq FILE` | none | Single-file mode. Inspect one sequence or alignment file in detail. |
| `--per-gene` | `False` | Directory mode only. Show per-gene results in terminal output when no `--output` is used, or save an adjacent per-gene table when `--output` is provided. |
| `--per-gene-format csv|tsv` | `csv` | Directory mode only. Format for the adjacent per-gene table written by `--per-gene --output`. |
| `--output PATH`, `-o PATH` | `runs/pretree-stats.json` | Save results to a file. Parent directories are created automatically. |
| `--output-format text|json` | `json` | Format for the file written with `--output`. This does not affect Rich terminal output. |
| `--input-format fasta|phylip-relaxed|nexus` | auto | Override automatic format detection. |
| `--seq-type AA|NT` | auto | Override automatic molecule-type detection. |
| `--threads INTEGER`, `-t INTEGER` | `4` | Directory mode only. Number of worker processes. Must be at least `1`. |
| `--quiet`, `-q` | `False` | Suppress terminal output except for errors. |

## Inputs

Single-file mode reads one file with `--seq`. Directory mode scans one directory with `--seq-dir` and processes supported sequence/alignment extensions.

Supported input formats are FASTA, Phylip-relaxed, and Nexus. Classic `phylip` is not exposed as a separate `stats` choice to avoid ambiguity with Phylip-relaxed.

Alignment status is detected per file. A file is aligned when it has more than one sequence and all sequences have identical length. Single-sequence files are treated as unaligned.

## Outputs

Terminal output uses Rich tables and panels unless `--quiet` is set. `--output-format` controls saved file content only.

Single-file mode shows overview, character summary, per-taxon statistics, and site-pattern statistics when the file is aligned.

Directory mode shows a summary table. With `--per-gene`, it also shows a per-gene table in the terminal when no output file is requested.

When `--output` is used, the command saves the structured result payload. With `--per-gene --output`, the per-gene table is saved next to the summary file as `<stem>.per-gene.csv` by default or `<stem>.per-gene.tsv` when `--per-gene-format tsv` is set.

## Examples

Inspect one aligned amino-acid file:

```bash
phyloai pretree stats --seq ref/phylogenomics_examples/test/EOG090X0971.faa
```

Inspect one unaligned nucleotide file:

```bash
phyloai pretree stats --seq ref/phylogenomics_examples/2-loci_filter/fna/EOG090X0971.fna
```

Summarize a directory:

```bash
phyloai pretree stats --seq-dir ref/phylogenomics_examples/3-align/faa
```

Save a directory summary plus per-gene CSV:

```bash
phyloai pretree stats \
  --seq-dir ref/phylogenomics_examples/2-loci_filter/fna \
  --per-gene \
  --output results/summary.json
```

Recommended order after raw input normalization:

```bash
phyloai pretree convert --input ./raw --output-dir ./runs/run001/pretree/convert --to fasta
phyloai pretree stats --seq-dir ./runs/run001/pretree/convert
```

## Warnings And Errors

`--seq` and `--seq-dir` are mutually exclusive. Supplying neither or both is an input error.

`--per-gene` is directory mode only. Using it with `--seq` is an input error.

`--threads` must be at least `1`.

If sequence type detection is ambiguous, the command defaults to `AA` and emits a warning.

If `*` appears in any sequence, the command emits a warning because it may indicate stop codons or upstream processing problems.

## Notes

Character classes are reported as standard, gap/missing, or ambiguous. Site-pattern statistics exclude gap/missing/ambiguous characters when determining parsimony-informative and singleton sites.

Per-gene CSV/TSV output removes columns that are completely empty for that export. Aligned-only datasets keep `alignment_length` and drop unaligned length columns; unaligned-only datasets keep `seq_length_*` columns and drop `alignment_length`.
