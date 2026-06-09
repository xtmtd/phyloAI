# PhyloAI

An AI-native modular phylogenomics analysis platform.

## Installation

```bash
pip install -e .
```

## Quick Start

```bash
phyloai doctor
```

The top-level CLI currently provides:

```bash
phyloai --help
```

- `phyloai doctor`: inspect external tool availability.
- `phyloai pretree stats`: inspect one sequence/alignment file or summarize a directory of files before downstream analysis.

## `phyloai doctor`

`phyloai doctor` checks whether external tools expected by the local project environment can be found, where they are located, and whether a version string can be detected.

### Usage

```bash
phyloai doctor [--output-format text|json]
```

### Parameters

- `--output-format text|json`: choose human-readable terminal output or machine-readable JSON. Default is `text`.

### Examples

Human-readable table:

```bash
phyloai doctor
```

JSON for scripting:

```bash
phyloai doctor --output-format json
```

### What it checks

The current tool registry includes these names:

- Required: `iqtree3`, `mafft`, `trimal`
- Optional: `astral-hybrid`, `pb_mpi`, `mcmctree`, `correction_multi.jl`, `run_treeshrink.py`, `magus`, `clipkit`, `phykit`, `java`, `julia`, `bmge`

`bmge` is displayed as `BMGE.jar` in the table because the local project may provide it as a bundled jar file.

### Notes

- `doctor` reports what the current shell environment can see. If a tool is missing in `doctor`, later workflow steps that depend on it will also fail.
- A missing optional tool is reported, but does not mean the whole installation is unusable.
- JSON output is useful for wrappers or CI checks that want structured status, path, version, and install notes.

## `phyloai pretree stats`

`phyloai pretree stats` is a read-only QC and inspection command for sequence and alignment files. It supports:

- Single-file mode with `--seq`
- Directory mode with `--seq-dir`
- Automatic aligned vs unaligned detection
- Automatic AA vs NT detection, with manual override when needed
- Text or JSON terminal output
- Optional saved output files for summaries and tables

Exactly one of `--seq` or `--seq-dir` is required.

### Usage

```bash
phyloai pretree stats [OPTIONS]
```

### Parameters

- `--seq-dir DIRECTORY`: directory mode. Scan all supported sequence/alignment files in one folder and compute a dataset-level summary.
- `--seq FILE`: single-file mode. Inspect one sequence or alignment file in detail.
- `--per-gene`: directory mode only. Show per-gene results in terminal output when no `--output` is used, or save an adjacent per-gene table when `--output` is provided.
- `--per-gene-format csv|tsv`: directory mode only. Format for the adjacent per-gene table written by `--per-gene --output`. Default is `csv`.
- `--output PATH`, `-o PATH`: save results to a file.
- `--output-format text|json`: terminal output format. Default is `text`.
- `--input-format fasta|phylip-relaxed|nexus`: override automatic format detection.
- `--seq-type AA|NT`: override automatic molecule type detection.
- `--threads INTEGER`, `-t INTEGER`: directory mode only. Number of worker processes. Default is `4`.
- `--quiet`, `-q`: suppress terminal output except for errors.

### Single-file examples

Aligned amino-acid file:

```bash
phyloai pretree stats --seq ref/phylogenomics_examples/test/EOG090X0971.faa
```

Unaligned nucleotide file:

```bash
phyloai pretree stats --seq ref/phylogenomics_examples/2-loci_filter/fna/EOG090X0971.fna
```

Force JSON on stdout:

```bash
phyloai pretree stats \
  --seq ref/phylogenomics_examples/test/EOG090X0971.faa \
  --output-format json
```

Save JSON even if the file suffix is `.txt`:

```bash
phyloai pretree stats \
  --seq ref/phylogenomics_examples/test/EOG090X0971.faa \
  --output out.txt \
  --output-format json
```

### Directory examples

Summarize an aligned directory and print only the summary table:

```bash
phyloai pretree stats --seq-dir ref/phylogenomics_examples/3-align/faa
```

Show terminal summary plus per-gene table:

```bash
phyloai pretree stats \
  --seq-dir ref/phylogenomics_examples/3-align/faa \
  --per-gene
```

Use more workers for a larger dataset:

```bash
phyloai pretree stats \
  --seq-dir ref/phylogenomics_examples/2-loci_filter/fna \
  --threads 8
```

Save summary plus adjacent per-gene CSV:

```bash
phyloai pretree stats \
  --seq-dir ref/phylogenomics_examples/2-loci_filter/fna \
  --per-gene \
  --output out.txt
```

Save summary plus adjacent per-gene TSV:

```bash
phyloai pretree stats \
  --seq-dir ref/phylogenomics_examples/2-loci_filter/fna \
  --per-gene \
  --per-gene-format tsv \
  --output out.txt
```

### Output behavior

Single-file mode:

- Text terminal output uses Rich panels and a per-taxon table.
- `--output file.csv` or `file.tsv` saves the per-taxon table.
- `--output file.txt` saves key-value summary content plus a `[per_taxon]` section.
- `--output-format json` prints JSON to stdout and also saves JSON to `--output`, even if the output suffix is not `.json`.

Directory mode:

- Text terminal output shows a Rich summary table.
- If `--per-gene` is set and `--output` is not set, a terminal per-gene table is also shown.
- In text mode, a transient progress bar is shown while files are processed unless `--quiet` is used.
- `--output file.csv`, `file.tsv`, or `file.txt` saves the dataset summary.
- If `--per-gene --output file` is used, the per-gene table is saved separately as:
  - `file.per-gene.csv` by default
  - `file.per-gene.tsv` when `--per-gene-format tsv` is set
- Per-gene CSV/TSV output dynamically removes columns that are completely empty for that specific export. This keeps aligned-only and unaligned-only tables compact.

### Per-gene columns

Common per-gene columns may include:

- `gene`
- `n_taxa`
- `n_taxa_ratio`
- `length_type`
- `alignment_length`
- `seq_length_min`
- `seq_length_max`
- `seq_length_mean`
- `seq_length_median`
- `seq_length_stdev`
- `gap_ratio`
- `ambiguous_ratio`
- `gap_ambiguous_ratio`
- `missing_taxa`
- `missing_taxa_ratio`

Not every export includes every column. For example:

- Fully aligned datasets keep `alignment_length` and drop the empty unaligned length columns.
- Fully unaligned datasets keep the `seq_length_*` columns and drop `alignment_length`.
- Mixed datasets may legitimately contain both groups.

### Notes and caveats

- `--per-gene` is directory mode only. Using it with `--seq` is an input error.
- `--threads` must be at least `1`.
- Sequence type detection defaults to `AA` if the command cannot distinguish AA vs NT from observed characters.
- `*` characters trigger a warning because they may indicate stop codons or upstream processing problems.
- Site-pattern statistics such as `distinct_patterns`, `constant_sites`, `parsimony_informative`, and `singleton_sites` are only available for aligned files.
- Single-sequence files are treated as unaligned.

### Help

For the latest CLI wording and option help:

```bash
phyloai pretree stats --help
```
