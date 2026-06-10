# pretree stats — Design Specification

**Date:** 2026-06-09  
**Last updated:** 2026-06-10 (output-format default corrected to json; log policy clarified; --output renamed to --output-dir for consistency)  
**Status:** Approved for implementation — pending implementation update  
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

`--per-gene` controls whether per-gene results are displayed in terminal output when no output file is requested. When `--output FILE` is also provided, per-gene results are written to an adjacent table file named `<stem>.per-gene.<format>` instead of being mixed into the summary file. The adjacent table defaults to CSV; `--per-gene-format tsv` writes TSV.

### 2.2 Single-file mode (`--seq`)

Detailed statistics for one file. Supports both aligned and unaligned input. Auto-detects alignment status.

```bash
phyloai pretree stats --seq ./EOG090X0971.faa [--output report.txt]
```

`--seq-dir` and `--seq` are mutually exclusive.

Help text must make this mutual exclusivity explicit so users can see the mode split directly from `--help`.

---

## 3. Input Format Detection

### 3.1 Alignment format

Uses `core/formats.py` `FormatConverter.detect()`. Supports FASTA, Phylip-relaxed, Nexus. `--input-format` overrides auto-detection. For `pretree stats`, classic `phylip` is not exposed as a separate CLI choice to avoid ambiguity with `phylip-relaxed`.

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
`gap_ambiguous_ratio` = (`gap_ratio` + `ambiguous_ratio`)  
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
| `length_type`     | yes     | yes       | `alignment_length` for aligned files, `seq_length` for unaligned files |
| `alignment_length` | yes    | no        | alignment columns for aligned files                       |
| `seq_length_min`  | no      | yes       | minimum ungapped sequence length in the file              |
| `seq_length_max`  | no      | yes       | maximum ungapped sequence length in the file              |
| `seq_length_mean` | no      | yes       | mean ungapped sequence length in the file                 |
| `seq_length_median` | no    | yes       | median ungapped sequence length in the file               |
| `seq_length_stdev` | no     | yes       | sample standard deviation of ungapped sequence lengths    |
| `gap_ratio`       | yes     | yes       | mean gap_ratio across sequences in this file              |
| `ambiguous_ratio` | yes     | yes       | mean ambiguous_ratio across sequences in this file        |
| `gap_ambiguous_ratio` | yes | yes       | combined missing/problematic-character ratio (`gap_ratio + ambiguous_ratio`) |
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
| `gap_ambiguous_ratio`   | mean / median across all files                           |
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
| `character_summary` | `{standard_ratio, gap_ratio, ambiguous_ratio, gap_ambiguous_ratio}` across all chars |
| `per_taxon`         | per-sequence: `{name, raw_length, ungapped_length, gap_ratio, ambiguous_ratio, gap_ambiguous_ratio}` |

**Aligned-only additions:**

| Metric                    | Description                                                     |
|---------------------------|-----------------------------------------------------------------|
| `alignment_length`        | number of columns                                               |
| `distinct_patterns`       | count + ratio of unique full-column site patterns after normalizing `?` to `-` (IQ-TREE-compatible `distinct patterns`) |
| `constant_sites`          | count + ratio (`alignment_length - parsimony_informative - singleton_sites`; IQ-TREE-compatible summary) |
| `parsimony_informative`   | count + ratio (≥2 standard chars each with count ≥2)           |
| `singleton_sites`         | count + ratio (variable but not parsimony-informative)          |

Site pattern computation excludes gap/missing/ambiguous characters when determining whether a site is parsimony-informative or singleton. `distinct_patterns` follows the IQ-TREE summary convention and is computed from unique full-column patterns in the alignment after normalizing `?` and `-` to the same missing-state symbol. `constant_sites` follows the IQ-TREE summary convention and is reported as `alignment_length - parsimony_informative - singleton_sites`.

**Unaligned-only additions:**

| Metric          | Description                                        |
|-----------------|----------------------------------------------------|
| `seq_length`    | min / max / mean / median / stdev across all sequences |
| `total_length`  | sum of all ungapped sequence lengths               |

---

## 6. Parameters

| Parameter         | Short | Type     | Default     | Notes                                                    |
|-------------------|-------|----------|-------------|----------------------------------------------------------|
| `--seq-dir`       |       | Path     | —           | directory mode; mutually exclusive with `--seq`          |
| `--seq`           |       | Path     | —           | single-file mode; mutually exclusive with `--seq-dir`    |
| `--per-gene`      |       | flag     | False       | directory mode only: include per-gene results in terminal output when no `--output` is used; with `--output`, write an adjacent per-gene table |
| `--per-gene-format` |     | csv\|tsv | csv         | directory mode only: format for adjacent per-gene table written with `--per-gene --output` |
| `--output`        | `-o`  | Path     | —           | write results to file; extension controls text/csv/tsv/json unless `--output-format json` is used |
| `--output-format` |       | text\|json | json      | output format for stdout structured output and saved files; Rich terminal display is always on unless `--quiet`; `json` also saves JSON to `--output` regardless of file suffix |
| `--input-format`  |       | fasta\|phylip-relaxed\|nexus | auto-detect | override format detection; accepted values: `fasta`, `phylip-relaxed`, `nexus` |
| `--seq-type`      |       | AA\|NT   | auto-detect | override sequence type detection                         |
| `--threads`       | `-t`  | int      | 4           | directory mode only: files processed in parallel (ProcessPoolExecutor) |
| `--quiet`         | `-q`  | flag     | False       | suppress all output except errors                        |

No `--extra-args` (no external tool invoked).  
No `--overwrite` (read-only command, no output directory created).

---

## 7. Output Formats

### 7.1 Terminal output (always Rich, independent of `--output-format`)

`--output-format` controls the format of structured output written to stdout and saved files. It does **not** affect Rich terminal display. Rich tables, panels, and progress bars are always rendered to the terminal unless `--quiet` is set.

**Directory mode:**
- Rich progress bar while processing files, unless `--quiet` is set
- Rich table: summary statistics
- Rich table: per-gene table (if `--per-gene` and no `--output` path is given)
- `[WARN]` lines for any stop codons detected
- If `--output FILE` is used, terminal output must print the resolved output path after the summary so users can immediately locate the saved file
- If `--per-gene` and `--output FILE` are used together, the saved-file message must explicitly say that the per-gene table was written there

**Single-file mode:**
- Rich panels: character summary, site patterns (aligned), per-taxon table

### 7.2 File output (`--output FILE`)

Extension determines format:

| Extension      | Content                                                          |
|----------------|------------------------------------------------------------------|
| `.csv`         | summary table (directory mode); per-taxon table (single-file)  |
| `.tsv`         | same as CSV but tab-separated                                    |
| `.json`        | full structured output (summary + per-gene or full single-file stats) |
| `.txt`         | plain text; directory mode writes `[summary]`; single-file mode writes key-value summary plus `[per_taxon]` table |

Directory mode `.csv`/`.tsv` output contains summary only. If `--per-gene --output FILE` is used, the per-gene table is written next to `FILE` as `<stem>.per-gene.csv` by default or `<stem>.per-gene.tsv` when `--per-gene-format tsv` is set.

When `--output-format json` is used with `--output FILE`, the saved file is JSON even if `FILE` has a non-JSON suffix such as `.txt`.

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

In text terminal output mode, the CLI displays a Rich progress bar and advances it as each file result is collected. Progress is suppressed for `--quiet` and JSON terminal output.

`--threads 1` disables parallelism (single-process, useful for debugging).  
Default: 4. Recommended: match available CPU cores for large datasets (1000+ files).

Worker function must be a module-level function (not a lambda or closure) for pickling compatibility.

---

## 9. Integration with `core/`

- Format detection: `core/formats.py` `FormatConverter.detect()`
- No `Runner` invoked (pure Python, no external tools)
- No `RunRecord` entry (utility command)
- No log file is written. `stats` is a read-only utility that does not require or create a run directory. This is a documented exception to the general pipeline log policy (Section 9.6 of the main spec), justified by `stats` being a read-only utility with no output directory.

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
- distinct_patterns: 624
- constant_sites: 482
- parsimony_informative: 87
- singleton_sites: 473
