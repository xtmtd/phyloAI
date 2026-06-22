# `phyloai pretree concat` Design Specification

**Date:** 2026-06-13
**Status:** Draft

---

## 1. Purpose

Concatenate multiple MSA files into a supermatrix for downstream phylogenetic
inference. Provides occupancy filtering, recoding, codon-specific variants
(translation, 3rd-codon exclusion), outgroup reordering, and multi-format output.

## 2. CLI

```bash
phyloai pretree concat --msa-dir <dir> [OPTIONS]
```

### 2.1 Parameters

| Parameter | Type | Default | Description |
|---|---|---|---|
| `--msa-dir` | `PATH` (required) | — | Directory of input MSA files |
| `--output-dir` / `-o` | `PATH` | `runs/pretree/concat` | Output directory |
| `--prefix` | `str` | `matrix` | Prefix for output filenames |
| `--seq-type` | `AA|NT|CODON|auto` | `auto` | Sequence type |
| `--taxa-occupancy` | `float` (0.0–1.0) | `0.5` | Min taxon ratio for MSA inclusion |
| `--recoding` | `Choice[RY-nucleotide, Dayhoff-6/9/12/15/18, SandR-6, KGB-6]` | None | RY-nucleotide: NT only (A/G→R, C/T/U→Y). Others: AA only |
| `--outgroup` | `str` | None | Single taxon name to move to first position |
| `--to` | `fasta|phylip-relaxed|phylip-paml|nexus` | `fasta` | Output format |
| `--translate-codon` | flag | `False` | Also produce CDS→AA translated matrix (CODON only) |
| `--exclude-codon3` | flag | `False` | Also produce codon1+2 matrix (CODON only) |

Universal flags: `--quiet`, `--overwrite`, `--dry-run`. No `--config`, no
`--resume`, no parallelism.

### 2.2 CODON variants

When `--seq-type CODON`:

- **Original** codon concatenation is always produced.
- `--translate-codon` additionally translates CDS→AA (standard genetic code,
  Biopython) then concatenates to a separate matrix.
- `--exclude-codon3` additionally drops every 3rd position (indices 2,5,8,…)
  then concatenates the codon1+2 positions to a separate matrix.

Both variants share the same `--outgroup` and `--to` settings as the original.

### 2.3 Recoding

Applied to the **normalized and concatenated matrix** (not per-MSA).  Tables are
embedded in `concat.py` as `dict[str,str]` — no bundled data files.
Schemes and their letter codes follow PhyKIT
(<https://jlsteenwyk.com/PhyKIT/usage/index.html#cmd-alignment-recoding>).

Sequences are normalized (via `core/sequence_normalization.py`) *before*
recoding, so non-standard characters (`U/O/B/Z/J` etc.) are converted to
standard equivalents.  The recoding tables include `"X": "?"` (AA) and all
IUPAC ambiguous codes (NT) so that normalized input should not normally
trigger warnings.

Supported schemes:

| NT | AA |
|----|----|
| `RY-nucleotide` | `Dayhoff-6`, `Dayhoff-9`, `Dayhoff-12`, `Dayhoff-15`, `Dayhoff-18`, `SandR-6`, `KGB-6` |

Recoded matrices have `seq_type` set to `"other"` since the character alphabet
no longer matches AA or NT.  Statistics are limited to `gap_ratio` and site
patterns (distinct, constant, parsimony-informative, singleton).
`ambiguous_ratio` is always `0.0` and `standard_ratio` is `1 − gap_ratio`.

### 2.4 Occupancy filtering

`--taxa-occupancy THRESHOLD` filters **MSAs** (not taxa).  An MSA is retained
iff `len(msa_taxa) / len(total_taxa) >= THRESHOLD`.  Dropped MSAs are recorded
in `dropped_alignments.csv` inside the output directory.

### 2.5 Outgroup reordering

If `--outgroup TAXON` is given, that taxon's record is moved to position 0
(first line) in every output matrix.  Exit code 1 if the taxon is not found.

---

## 3. Data flow

```
[1]  Scan --msa-dir → list of MSA paths
[2]  Pass 1 (header-first): _read_msa_headers() on every file
     → taxon names per file → total_taxa (union) + per-file taxon counts
     (For FASTA inputs this is true header-only scanning; other formats may
     fall back to normal parsing when cheap header-only parsing is unavailable.
     Memory: O(T × N_filenames).)
[3]  Auto-detect seq_type (if --seq-type auto):
     _read_msa() on first 3 files → sample sequences → detect_seq_type()
[4]  Validate: CODON flags, recoding×seq_type compatibility
[5]  Occupancy filter: _filter_by_occupancy() → kept_paths + dropped list
[6]  Pass 2 (streaming concat, kept files only):
     for each kept file:
       _read_msa() → normalize_sequences() → stream-append to matrix[taxon]
     (Per-gene normalized seqs kept in memory only if CODON variants needed.)
[7]  Generate variants:
     variant_original  ← always (from step 6)
     variant_recoded   ← _apply_recoding() on normalized matrix (if --recoding)
     [CODON only, from cached normalized data]:
       variant_translated ← translate per-gene → concat → normalize(AA)
       variant_codon12    ← exclude_codon3 per-gene → concat → normalize(NT)
[8]  Apply --outgroup reordering to each variant
[9]  Write each variant matrix via core/formats.py + corresponding .partitions file
[10] Compute stats on every variant: gap_ratio, character summary, site patterns
[11] Rich display: Overview + Character Summary + Site Patterns
[12] Write result.json (tool_stderr inlined per JSON Output Standard Section 5.2)
```

---

## 4. Output files

| File | Condition | Content |
|---|---|---|
| `<prefix>.<ext>` | Always | Original concatenated matrix |
| `<prefix>.partitions` | Always | RAxML-style partition file for original matrix |
| `<prefix>.recoded.<ext>` | `--recoding` | Recoded matrix |
| `<prefix>.recoded.partitions` | `--recoding` | Partition file for recoded matrix |
| `<prefix>.translated.<ext>` | `--translate-codon` | CDS→AA matrix |
| `<prefix>.translated.partitions` | `--translate-codon` | Partition file for translated matrix |
| `<prefix>.cds12.<ext>` | `--exclude-codon3` | Codon1+2 matrix |
| `<prefix>.cds12.partitions` | `--exclude-codon3` | Partition file for cds12 matrix |
| `dropped_alignments.csv` | Any MSA dropped | CSV with columns: `filename,n_taxa,occupancy_ratio,total_taxa` |
| `result.json` | `not --dry-run` | Structured result |

Output directory structure:

```
runs/pretree/concat/
├── matrix.fa                   # (or .phy/.nex)
├── matrix.partitions           # partition file
├── matrix.recoded.fa           # if --recoding
├── matrix.recoded.partitions   # if --recoding
├── matrix.translated.fa        # if --translate-codon
├── matrix.translated.partitions # if --translate-codon
├── matrix.cds12.fa             # if --exclude-codon3
├── matrix.cds12.partitions     # if --exclude-codon3
├── dropped_alignments.csv      # if any dropped
└── result.json                # if not --dry-run
```

### 4.1 Partition files (RAxML-style)

Each matrix file is accompanied by a `.partitions` file describing gene
boundaries within the concatenated supermatrix. Format follows the RAxML
convention as used by IQ-TREE:

```
DNA, gene1 = 1-654
DNA, gene2 = 655-1203
DNA, gene3 = 1204-1980
```

**Prefix rule** (first column per line):

| Matrix variant | seq_type | Prefix |
|---|---|---|
| Original (NT/CODON) | nt / codon | `DNA` |
| Original (AA) | aa | `LG` |
| Recoded (any) | other | `AUTO` |
| Translated (CODON→AA) | aa | `LG` |
| Codon1+2 (CODON→NT) | nt | `DNA` |

**Partition naming:** Each partition is named after the gene file's basename
(without extension), e.g., `gene1.fa` → partition name `gene1`.

**File naming:** The partitions file shares the matrix file's stem. For
`<prefix>.<ext>`, the partitions file is `<prefix>.partitions`. For variant
matrices, the variant infix is preserved: `<prefix>.recoded.partitions`,
`<prefix>.translated.partitions`, `<prefix>.cds12.partitions`.

**Position tracking:** Gene start/end positions are tracked during Pass 2
(streaming concat) as 1-indexed inclusive ranges. For variant matrices
(recoded, translated, cds12), gene lengths may differ from the original, so
positions are recomputed from the per-gene variant data before re-concatenation.

**Edge cases:**
- Single gene kept → partitions file has one line (valid).
- Zero genes kept (all dropped by occupancy) → no matrix, no partitions.
- Partitions file only written if the corresponding matrix write succeeds.

---

## 5. Screen display (Rich)

Three panels displayed on the terminal:

1. **Overview** — prefix, to_format, n_taxa, n_msa_input,
   n_msa_used, n_msa_dropped, taxon_occupancy_threshold, recoding (if any),
   outgroup (if any), variant matrices produced.  (`seq_type` and
   `total_length` are per-variant and appear in Character Summary.)

2. **Character Summary** — per-variant table with columns for each variant
   matrix.  Rows: `seq_type`, `total_length`, `gap_ratio`, `ambiguous_ratio`,
   `gap_ambiguous_ratio`, `standard_ratio`.  Recoded variant (`seq_type =
   "other"`) shows `—` for meaningless metrics.

3. **Site Patterns** — per-variant table.  Rows: `alignment_length`,
   `distinct_patterns` (count + ratio), `constant_sites`,
   `parsimony_informative`, `singleton_sites`.  Ratios displayed to 4 decimal
   places.  Distinct-pattern counting collapses all non-standard characters
   (gaps, ambiguous codes) into the gap symbol, matching IQ-TREE's convention.

The following are **not shown on screen** (saved to `result.json` only when
files are written):
- Dropped MSAs list
- Per-taxon summary
- Per-gene occupancy

---

## 6. `result.json` schema

```json
{
  "status": "success",
  "command": "phyloai pretree concat --msa-dir ...",
  "wall_time": 12.3,
  "tool_versions": {},
  "params": {
    "msa_dir": "...",
    "output_dir": "...",
    "prefix": "matrix",
    "seq_type": "AA",
    "taxa_occupancy": 0.5,
    "recoding": "Dayhoff-6",
    "outgroup": "Sp_A",
    "to_format": "fasta",
    "translate_codon": false,
    "exclude_codon3": false
  },
  "key_results": {
    "n_taxa": 100,
    "n_msa_input": 500,
    "n_msa_used": 350,
    "n_msa_dropped": 150,
    "total_length": 250000,
    "variants_produced": ["runs/pretree/concat/matrix.fa", "runs/pretree/concat/matrix.recoded.fa"]
  },
  "error": null,
  "data": {
    "cmd": [],
    "tool_stderr": "",
    "character_summary": {
      "gap_ratio": 0.12,
      "ambiguous_ratio": 0.01,
      "gap_ambiguous_ratio": 0.13,
      "standard_ratio": 0.87
    },
    "site_patterns": {
      "alignment_length": 250000,
      "distinct_patterns": {"count": 120000, "ratio": 0.48},
      "constant_sites": {"count": 80000, "ratio": 0.32},
      "parsimony_informative": {"count": 30000, "ratio": 0.12},
      "singleton_sites": {"count": 20000, "ratio": 0.08}
    },
    "dropped_alignments": [
      {"filename": "gene123.fa", "n_taxa": 30, "occupancy_ratio": 0.30, "total_taxa": 100}
    ],
    "per_taxon": [
      {"name": "Sp_A", "n_msa_present": 350, "occupancy_ratio": 1.0,
       "raw_length": 250000, "ungapped_length": 220000, "gap_ratio": 0.12}
    ],
    "per_gene_occupancy": [
      {"gene": "gene001.fa", "n_present": 95, "n_missing": 5, "occupancy_ratio": 0.95}
    ],
    "variants": [
      {"variant": "original", "path": "runs/pretree/concat/matrix.fa", "seq_type": "AA", "length": 250000},
      {"variant": "recoded", "path": "runs/pretree/concat/matrix.recoded.fa", "seq_type": "other", "length": 250000}
    ],
    "variant_stats": [
      {
        "variant": "original",
        "seq_type": "AA",
        "total_length": 250000,
        "character_summary": { "gap_ratio": 0.12, ... },
        "site_patterns": { "alignment_length": 250000, "distinct_patterns": {...}, ... }
      },
      {
        "variant": "recoded",
        "seq_type": "other",
        "total_length": 250000,
        "character_summary": { "gap_ratio": 0.12, "ambiguous_ratio": 0.0, ... },
        "site_patterns": { ... }
      }
    ],
    "recoding_warnings": [],
    "normalization_replacements": {"u_to_t": 15, "aa_special_to_x": 3}
  }
}
```

`key_results.variants_produced` and each `data.variants[].path` entry are full
output paths constructed from `params.output_dir`, not bare filenames.  This
keeps downstream report collection and MCP wrappers independent of the caller's
current working directory.

---

## 7. Implementation notes

### 7.1 Module: `phyloai/pretree/concat.py`

New file. Main function `run_concat(...)` returns the result payload dict.
No PhyKIT dependency.  Internal helper functions:

| Function | Purpose |
|---|---|---|
| `_scan_msa_files(msa_dir)` | Glob for alignment files |
| `_read_msa_headers(path)` | Taxon-name scan; header-only for FASTA, parser fallback for other formats |
| `_read_msa(path)` | Parse one MSA → `dict[taxon_id, seq_str]` + length |
| `_concat_alignments(msa_paths, msa_data, total_taxa)` | Build supermatrix dict from pre-loaded data |
| `_filter_by_occupancy(msa_paths, total_taxa, threshold)` | Split kept/dropped |
| `_apply_recoding(matrix, scheme)` | Character-level recoding |
| `_translate_codon(seq)` | Codon→AA translation (standard code) |
| `_exclude_codon3(seq)` | Drop every 3rd position |
| `_reorder_outgroup(matrix, outgroup)` | Move taxon to index 0 |
| `_write_matrix(matrix, path, fmt, seq_type)` | Via `core/formats.py` |
| `_write_partitions(partitions_path, genes, prefix_type)` | Write RAxML-style `.partitions` file |
| `_compute_concat_stats(matrix, seq_type)` | Per-variant stats (gap_ratio, site patterns) |
| `_render_concat_panels(overview, variant_stats)` | Rich panels for terminal display |

### 7.2 Recoding tables

Reference: PhyKIT documentation on `alignment recoding`
<https://jlsteenwyk.com/PhyKIT/usage/index.html#cmd-alignment-recoding>.
This is the canonical source for the published recoding schemes; tables
embedded in code must match PhyKIT exactly so downstream users get
equivalent results.

Embedded as module-level `dict` in `phyloai/pretree/concat.py`.  Encoded
state labels (digits `0`–`9`, letters `A`–`H`) match PhyKIT exactly so
that recoded matrices are byte-for-byte equivalent to `phykit recode`.

```python
AA_RECODING_TABLES: dict[str, dict[str, str]] = {
    "Dayhoff-6":  {"A": "0", "G": "0", "P": "0", "S": "0", "T": "0",
                   "D": "1", "E": "1", "N": "1", "Q": "1",
                   "H": "2", "K": "2", "R": "2",
                   "I": "3", "L": "3", "M": "3", "V": "3",
                   "F": "4", "W": "4", "Y": "4",
                   "C": "5", "X": "?"},
    ...
}
NT_RECODING_TABLES: dict[str, dict[str, str]] = {
    "RY-nucleotide": {"A": "R", "G": "R", "C": "Y", "T": "Y", "U": "Y",
                      "N": "?", "X": "?",
                      "R": "R", "Y": "Y",
                      "S": "?", "W": "?", "K": "?", "M": "?",
                      "B": "?", "D": "?", "H": "?", "V": "?",
                      "-": "-", "?": "?", ".": "."},
}
```

Notes on embedding:

- Sequences are **normalized before recoding** (via
  `core/sequence_normalization.py`, mirroring `pretree convert` behavior).
  Non-standard characters (`U/O/B/Z/J`, invalid bases) are replaced with
  standard equivalents or `N`/`X`.  Recoding thus operates on clean input
  with a known character set.
- Each dict covers the **post-normalization character set** for its seq_type:
  - AA: 20 standard amino acids + `X`, plus gap chars (`-`, `?`, `.`, `*`).
  - NT: `A,C,G,T` + all IUPAC ambiguity codes (`R,Y,S,W,K,M,B,D,H,V,N`) +
    `X`, plus gap chars.
- Gap characters (`-`, `?`, `.`, `*`) are **preserved** (never recoded) via a
  dedicated gap check that runs before table lookup.
- Characters not in the table are passed through unchanged and a warning is
  recorded in `result.json.data.recoding_warnings`.  After normalization,
  this should only occur for truly unexpected edge cases.
- `--recoding` is **rejected at CLI validation time** (exit code 1) if the
  scheme is incompatible with the resolved `seq_type`:
  - AA-only schemes (`Dayhoff-*`, `SandR-6`, `KGB-6`) rejected for NT.
  - NT-only scheme (`RY-nucleotide`) rejected for AA.

### 7.3 Reuse existing modules

- `phyloai/core/formats.py` — `FormatConverter.write_alignment()` for
  writing matrices in all target formats.
- `phyloai/pretree/stats.py` — `file_stats_aligned()`,
  `compute_site_patterns()`, `render_single_file_panels()` for Rich display.
- `phyloai/core/sequence_normalization.py` — `detect_seq_type()`,
  `expand_dots_from_first_sequence()`, `normalize_sequences()`.
- `phyloai/core/schema.py` — `COMMON_ALIGNMENT_EXTENSIONS` for file scanning.
- `phyloai/cli/commands/pretree.py` — Add `concat_command()` following the
  existing pattern (`convert_command`, `trim_command`).

For character/site statistics, `CODON` matrices are treated as nucleotide
matrices (`NT`) because the underlying alphabet is still nucleotide-based.

### 7.4 CODON translation

Use Biopython's `Bio.Seq.Seq.translate()` with the standard genetic code
(table 1). Gaps (`-`) in codon sequences pass through as gaps in the
translated sequence (translated as `-`). For codon sequences whose length is
not a multiple of 3, the remainder is trimmed: only complete codons are
translated.

### 7.5 CODON3 exclusion

For each sequence of length N, keep positions `[0,1,3,4,6,7,…]` — remove
every 3rd position (index 3k+2 for k≥0).

### 7.6 Error handling

Exit codes follow Section 9.3:
- Code 1: missing/invalid input, empty MSA directory, outgroup not found,
  CODON variant requested for non-CODON seq_type, invalid recoding scheme,
  output directory conflict.
- Code 3: no external tool dependencies for this command.

### 7.7 CLI integration

Add `"concat"` to `_PretreeGroup.list_commands()` and register
`concat_command()` in `phyloai/cli/commands/pretree.py`. The CLI handler
validates parameters then delegates to `run_concat()`.

---

## 8. Acceptance criteria

- [ ] `--seq-type auto` resolves to `AA` or `NT` via `detect_seq_type()`; `CODON` is never auto-detected (must be explicit)
- [ ] `--seq-type auto` detection samples up to 3 files; resolved type recorded in `result.json` `params.seq_type`
- [ ] `--taxa-occupancy 0.5` drops MSAs below 50% total-taxon ratio
- [ ] `--taxa-occupancy 0.0` keeps all MSAs
- [ ] `--taxa-occupancy 1.0` keeps only MSAs containing 100% of taxa
- [ ] Dropped MSAs written to `dropped_alignments.csv`
- [ ] For FASTA inputs, Pass 1 uses header-only scanning; for other formats it may fall back to normal parsing; dropped files are skipped in Pass 2
- [ ] Streaming concat (Pass 2) normalizes and appends per-gene data without storing all raw sequences
- [ ] Original matrix written in FASTA format by default
- [ ] `--to phylip-relaxed` / `phylip-paml` / `nexus` produce correct output
- [ ] `--recoding Dayhoff-6` produces recoded matrix file
- [ ] `--recoding RY-nucleotide` with NT sequences works
- [ ] `--seq-type CODON --translate-codon` produces AA-translated matrix
- [ ] `--seq-type CODON --exclude-codon3` produces codon1+2 matrix
- [ ] `--outgroup Sp_X` moves that taxon to position 0 in all variant matrices
- [ ] Rich screen display shows Overview + Character Summary + Site Patterns
- [ ] `result.json` includes dropped MSAs, per-taxon, per-gene occupancy
- [ ] `data.tool_stderr` recorded in `result.json` (single mode, JSON Output Standard Section 5.2)
- [ ] `--dry-run` validates inputs and computes an in-memory summary without writing, deleting, or replacing any files or directories
- [ ] `--quiet` suppresses terminal output but still writes `result.json`
- [ ] Non-empty output directory without `--overwrite` → exit code 1
- [ ] `--overwrite` recreates the output directory
- [ ] Each matrix file has a corresponding `.partitions` file with RAxML-style gene boundaries
- [ ] NT/CODON matrices use `DNA` prefix; AA matrices use `LG`; recoded matrices use `AUTO`
- [ ] Partition names use gene file basenames (without extension)
- [ ] Partitions file not written under `--dry-run`
