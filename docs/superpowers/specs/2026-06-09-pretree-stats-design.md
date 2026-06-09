# pretree stats — Design Specification

**Date:** 2026-06-09  
**Status:** Approved for implementation  
**Parent spec:** `2026-06-07-phyloai-design.md`

---

## 1. Purpose

`phyloai pretree stats` provides read-only batch or single-file statistics for sequence and alignment files. It does not modify any files. It serves two purposes:

1. **Quality inspection** before entering the analysis pipeline (unaligned or aligned)
2. **Diagnostic tool** at any pipeline stage to understand data characteristics

---

## 2. Two Operating Modes

### 2.1 Directory mode (`--seq-dir`)

Scans all sequence files in a directory. Files may be aligned or unaligned — auto-detected per file. Computes per-file statistics in parallel, then aggregates a directory-level summary.

```bash
phyloai pretree stats --seq-dir ./aligned [--per-gene] [--output report.csv] [--threads 4]
```

### 2.2 Single-file mode (`--seq`)

Detailed statistics for one file. Supports both aligned and unaligned input. Auto-detects alignment status.

```bash
phyloai pretree stats --seq ./EOG090X0971.faa [--output report.txt]
```

`--seq-dir` and `--seq` are mutually exclusive.

---

## 3. Input Format Detection

### 3.1 Alignment format

Uses `core/formats.py` `FormatConverter.detect()`. Supports FASTA, Phylip-relaxed, Nexus. `--input-format` overrides auto-detection.

### 3.2 Is-aligned detection

A file is considered **aligned** if all sequences have identical length after reading. Unequal lengths → unaligned. Single-sequence files are treated as unaligned.

### 3.3 Sequence type (AA vs NT) auto-detection

Inspect the set of standard characters (excluding gap/ambiguous) across all sequences:

- If any of `E F I L P Q W Y Z` are present → **AA**
- If all standard characters are subset of `A C G T U` → **NT**
- If ambiguous → emit `[WARN] Cannot determine seq_type, defaulting to AA` and proceed

`--seq-type AA|NT` overrides auto-detection. Affects which characters are classified as ambiguous.

---

## 4. Character Classification

Three classes, applied uniformly across both modes:

| Class       | AA characters                                     | NT characters                          | Treatment in site-pattern computation |
|-------------|---------------------------------------------------|----------------------------------------|----------------------------------------|
| standard    | `ACDEFGHIKLMNPQRSTVWY`                            | `ACGT`                                 | participate normally                   |
| gap/missing | `-` `?`                                           | `-` `?` `N`                            | excluded                               |
| ambiguous   | `BZJUO` + `*` + any other non-standard character  | `RYSWKMBDHVU` + `*` + any other        | excluded                               |

**Special warning:** if `*` (stop codon) is detected in any sequence, emit:
```
[WARN] Stop codon (*) found in <filename>. This may indicate upstream processing errors.
```

`gap_ratio` = (gap/missing chars) / total chars  
`ambiguous_ratio` = ambiguous chars / total chars  
`standard_ratio` = standard chars / total chars  
(three ratios sum to 1.0)

---

## 5. Statistics Computed

### 5.1 Directory mode — per-file metrics (basis for summary and per-gene table)

| Metric            | Aligned | Unaligned | Description                                               |
|-------------------|---------|-----------|-----------------------------------------------------------|
| `gene`            | yes     | yes       | filename without extension                                |
| `n_taxa`          | yes     | yes       | number of sequences                                       |
| `n_taxa_ratio`    | yes     | yes       | n_taxa / max_taxa_in_dataset                              |
| `is_aligned`      | yes     | yes       | boolean                                                   |
| `length`          | yes     | yes       | alignment length (aligned) or median seq length (unaligned) |
| `gap_ratio`       | yes     | yes       | mean gap_ratio across sequences in this file              |
| `ambiguous_ratio` | yes     | yes       | mean ambiguous_ratio across sequences in this file        |
| `missing_taxa`    | yes     | yes       | max_taxa_in_dataset − n_taxa                              |
| `missing_taxa_ratio` | yes  | yes       | missing_taxa / max_taxa_in_dataset                        |

### 5.2 Directory mode — summary (aggregated across all files)

| Metric                  | Description                                              |
|-------------------------|----------------------------------------------------------|
| `n_genes`               | total number of files processed                          |
| `format`                | detected format (or "mixed" if inconsistent)             |
| `is_aligned`            | True if all files aligned; False if all unaligned; "mixed" otherwise |
| `seq_type`              | AA / NT / mixed                                          |
| `total_taxa`            | number of unique taxon names across all files            |
| `taxa_per_gene`         | min / max / mean / median                                |
| `length`                | min / max / mean / median / total (alignment or seq length) |
| `gap_ratio`             | mean / median across all files                           |
| `ambiguous_ratio`       | mean / median across all files                           |
| `missing_taxa_ratio`    | mean / median across all files                           |
| `warnings`              | list of any stop-codon or unexpected character warnings  |

### 5.3 Single-file mode — metrics

**Common (aligned and unaligned):**

| Metric              | Description                                                     |
|---------------------|-----------------------------------------------------------------|
| `filename`          | input file path                                                 |
| `format`            | detected format                                                 |
| `seq_type`          | AA / NT                                                         |
| `is_aligned`        | boolean                                                         |
| `n_taxa`            | number of sequences                                             |
| `taxon_names`       | list of all taxon names                                         |
| `character_summary` | `{standard_ratio, gap_ratio, ambiguous_ratio}` across all chars |
| `per_taxon`         | per-sequence: `{name, raw_length, ungapped_length, gap_ratio, ambiguous_ratio}` |

**Aligned-only additions:**

| Metric                    | Description                                                     |
|---------------------------|-----------------------------------------------------------------|
| `alignment_length`        | number of columns                                               |
| `constant_sites`          | count + ratio                                                   |
| `variable_sites`          | count + ratio (all sites where ≥2 distinct standard chars)      |
| `parsimony_informative`   | count + ratio (≥2 standard chars each with count ≥2)           |
| `singleton_sites`         | count + ratio (variable but not parsimony-informative)          |
| `gap_only_sites`          | count + ratio (columns where all chars are gap/missing/ambiguous) |

Site pattern computation excludes gap/missing/ambiguous characters per site. Sites where fewer than 2 sequences have standard characters are excluded from variable/PIS/singleton counts and counted as `gap_only_sites`.

**Unaligned-only additions:**

| Metric          | Description                                        |
|-----------------|----------------------------------------------------|
| `seq_length`    | min / max / mean / median across all sequences     |
| `total_length`  | sum of all ungapped sequence lengths               |

---

## 6. Parameters

| Parameter         | Short | Type     | Default     | Notes                                                    |
|-------------------|-------|----------|-------------|----------------------------------------------------------|
| `--seq-dir`       |       | Path     | —           | directory mode; mutually exclusive with `--seq`          |
| `--seq`           |       | Path     | —           | single-file mode; mutually exclusive with `--seq-dir`    |
| `--per-gene`      |       | flag     | False       | directory mode only: include per-gene table in output    |
| `--output`        | `-o`  | Path     | —           | write results to file; format inferred from extension    |
| `--output-format` |       | text\|json | text      | terminal output format                                   |
| `--input-format`  |       | str      | auto-detect | override format detection                                |
| `--seq-type`      |       | AA\|NT   | auto-detect | override sequence type detection                         |
| `--threads`       | `-t`  | int      | 4           | directory mode only: files processed in parallel (ProcessPoolExecutor) |
| `--quiet`         | `-q`  | flag     | False       | suppress all output except errors                        |

No `--extra-args` (no external tool invoked).  
No `--overwrite` (read-only command, no output directory created).

---

## 7. Output Formats

### 7.1 Terminal output (default `--output-format text`)

**Directory mode:**
- Rich table: summary statistics
- Rich table: per-gene table (if `--per-gene`)
- `[WARN]` lines for any stop codons detected

**Single-file mode:**
- Rich panels: character summary, site patterns (aligned), per-taxon table

### 7.2 File output (`--output FILE`)

Extension determines format:

| Extension      | Content                                                          |
|----------------|------------------------------------------------------------------|
| `.csv`         | per-gene table (directory mode); per-taxon table (single-file)  |
| `.tsv`         | same as CSV but tab-separated                                    |
| `.json`        | full structured output (summary + per-gene or full single-file stats) |
| `.txt`         | plain text, same as terminal output without Rich formatting      |

Directory mode without `--per-gene`: `.csv`/`.tsv` output contains summary only (single-row or key-value format).

### 7.3 JSON output schema (`--output-format json` or `--output FILE.json`)

```json
{
  "status": "success",
  "command": "phyloai pretree stats --seq-dir ./aligned",
  "wall_time": 4.2,
  "tool_versions": {},
  "params": {
    "seq_dir": "./aligned",
    "threads": 4,
    "seq_type": "auto"
  },
  "key_results": {},
  "error": null,
  "data": {
    "summary": { },
    "per_gene": [ ]
  }
}
```

`key_results` is empty (`{}`): `stats` is a utility command and does not contribute to `report`.

---

## 8. Parallelism

Directory mode uses `concurrent.futures.ProcessPoolExecutor` with `--threads` workers. Each worker receives a file path and returns a stats dict. The main process collects results and computes summary aggregations.

`--threads 1` disables parallelism (single-process, useful for debugging).  
Default: 4. Recommended: match available CPU cores for large datasets (1000+ files).

Worker function must be a module-level function (not a lambda or closure) for pickling compatibility.

---

## 9. Integration with `core/`

- Format detection: `core/formats.py` `FormatConverter.detect()`
- No `Runner` invoked (pure Python, no external tools)
- No `RunRecord` entry (utility command)
- Warnings emitted via `core/logger.py` if `--run-dir` is set; otherwise printed to stderr

---

## 10. Test Data

| Path | Type | Notes |
|------|------|-------|
| `ref/phylogenomics_examples/2-loci_filter/faa/` | unaligned AA | 1066 files, 4–6 taxa |
| `ref/phylogenomics_examples/2-loci_filter/fna/` | unaligned NT | 1066 files, 4–6 taxa |
| `ref/phylogenomics_examples/3-align/faa/` | aligned AA | 1066 files, 4–6 taxa |
| `ref/phylogenomics_examples/3-align/fna/` | aligned NT | 1066 files, 4–6 taxa |
| `ref/phylogenomics_examples/test/EOG090X0971.faa` | aligned AA single file | 6 taxa, length 1042 |
| `ref/phylogenomics_examples/test/EOG090X0971.fna` | aligned NT single file | 6 taxa, length 3126 |

Expected values for `EOG090X0971.faa` (verified manually):
- alignment_length: 1042
- constant_sites: 242
- variable_sites: 560
- parsimony_informative: 87
- singleton_sites: 473
- gap_only_sites: 240 (sites with < 2 standard chars)
