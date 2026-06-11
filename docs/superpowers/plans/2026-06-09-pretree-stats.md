# pretree stats — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Implement `phyloai pretree stats` — a read-only batch and single-file sequence statistics command covering directory mode (`--seq-dir`) and single-file mode (`--seq`), with parallel computation and Rich terminal output.

**Design spec:** `docs/superpowers/specs/2026-06-09-pretree-stats-design.md`

**Tech stack:** Python ≥3.10, `biopython`, `rich`, `click`, `concurrent.futures.ProcessPoolExecutor`

**Test data:**
- `ref/phylogenomics_examples/2-loci_filter/faa/` — unaligned AA, 1066 files
- `ref/phylogenomics_examples/3-align/faa/` — aligned AA, 1066 files
- `ref/phylogenomics_examples/3-align/fna/` — aligned NT, 1066 files
- `ref/phylogenomics_examples/test/EOG090X0971.faa` — single aligned AA file
- `ref/phylogenomics_examples/test/EOG090X0971.fna` — single aligned NT file

---

## File Map

| File | Responsibility |
|------|----------------|
| `phyloai/pretree/__init__.py` | Package init |
| `phyloai/pretree/stats.py` | Core stats logic: character classification, site pattern computation, summary aggregation |
| `phyloai/cli/commands/pretree.py` | `pretree` Click group + `stats` subcommand |
| `tests/pretree/test_stats.py` | Unit and integration tests |

---

## Task 1: pretree package scaffold

**Files:**
- Create: `phyloai/pretree/__init__.py`
- Create: `phyloai/cli/commands/__init__.py`
- Create: `phyloai/cli/commands/pretree.py` (stub only — `pretree` group, no commands yet)
- Modify: `phyloai/cli/main.py` — register `pretree` group

**Steps:**
- [ ] Create `phyloai/pretree/__init__.py` (empty)
- [ ] Create `phyloai/cli/commands/__init__.py` (empty)
- [ ] Create `phyloai/cli/commands/pretree.py` with a bare `pretree` Click group:
  ```python
  import click
  @click.group()
  def pretree():
      """Pre-tree data preparation commands."""
      pass
  ```
- [ ] In `phyloai/cli/main.py`, import and register the `pretree` group:
  ```python
  from phyloai.cli.commands.pretree import pretree
  cli.add_command(pretree)
  ```
- [ ] Verify: `phyloai pretree --help` shows the group with no subcommands yet

---

## Task 2: Character classification and sequence type detection

**File:** `phyloai/pretree/stats.py`

**Steps:**
- [ ] Define character sets as module-level constants:
  ```python
  AA_STANDARD  = set("ACDEFGHIKLMNPQRSTVWY")
  AA_AMBIGUOUS = set("BZJUO*")  # plus any unrecognised char
  NT_STANDARD  = set("ACGT")
  NT_AMBIGUOUS = set("RYSWKMBDHVU*")  # plus any unrecognised char
  GAP_CHARS    = set("-?")
  NT_MISSING   = set("-?N")  # N is missing in NT only
  ```
- [ ] Implement `detect_seq_type(sequences: list[str]) -> str`:
  - Collect all characters from all sequences (upper-cased), exclude gap/ambiguous
  - If any of `EFILPQWYZ` present → return `"AA"`
  - If all standard chars subset of `ACGTU` → return `"NT"`
  - Else → emit warning, return `"AA"`
- [ ] Implement `classify_char(char: str, seq_type: str) -> str`:
  - Returns `"standard"`, `"gap"`, or `"ambiguous"`
  - Any char not in the standard or gap sets falls into `"ambiguous"`
- [ ] Implement `check_stop_codons(sequences: list[str], filename: str) -> list[str]`:
  - Returns list of warning strings if `*` found in any sequence
- [ ] Write unit tests in `tests/pretree/test_stats.py`:
  - [ ] `test_detect_seq_type_aa` — standard AA input
  - [ ] `test_detect_seq_type_nt` — standard NT input
  - [ ] `test_detect_seq_type_ambiguous_falls_back_to_aa`
  - [ ] `test_classify_char_standard`
  - [ ] `test_classify_char_gap`
  - [ ] `test_classify_char_ambiguous`
  - [ ] `test_classify_char_unexpected_is_ambiguous`
  - [ ] `test_stop_codon_warning`

---

## Task 3: Per-sequence and per-file statistics (unaligned)

**File:** `phyloai/pretree/stats.py`

**Steps:**
- [ ] Implement `per_taxon_stats(record, seq_type: str) -> dict`:
  - Inputs: a `SeqRecord`
  - Returns: `{name, raw_length, ungapped_length, gap_ratio, ambiguous_ratio, standard_ratio}`
  - `raw_length` = len(seq)
  - `ungapped_length` = count of non-gap characters
  - ratios computed over `raw_length`
- [ ] Implement `file_stats_unaligned(path: Path, seq_type: str | None, input_format) -> dict`:
  - Read sequences with `FormatConverter.read()`
  - Detect seq_type if not provided
  - Check stop codons
  - Compute `per_taxon_stats` for each sequence
  - Aggregate: `n_taxa`, `seq_length` (min/max/mean/median/stdev of ungapped lengths), `total_length`, `gap_ratio` (mean across taxa), `ambiguous_ratio` (mean across taxa), `gap_ambiguous_ratio` (mean gap + ambiguous), `standard_ratio` (mean), `warnings`
  - Return full dict including `is_aligned: False`
- [ ] Write tests:
  - [ ] `test_per_taxon_stats_basic`
  - [ ] `test_file_stats_unaligned` using `2-loci_filter/faa/EOG090X0007.faa`
  - [ ] Verify n_taxa=5, is_aligned=False

---

## Task 4: Per-file statistics (aligned) + site pattern computation

**File:** `phyloai/pretree/stats.py`

**Steps:**
- [ ] Implement `compute_site_patterns(sequences: list[str], seq_type: str) -> dict`:
  - Compute `distinct_patterns` as the number of unique full-column patterns across the alignment after normalizing `?` to `-` (IQ-TREE-compatible)
  - Iterate over each column index
  - For each column, collect characters classified as `standard` only
  - If fewer than 2 standard chars in column → skip variable/PIS/singleton classification for that site
  - Else classify column:
    - `parsimony_informative`: ≥2 distinct standard chars each appearing ≥2 times
    - `singleton`: variable but not parsimony-informative
    - `constant` = alignment_length − parsimony_informative − singleton_sites (IQ-TREE-compatible summary)
  - Return counts and ratios (denominator = alignment_length)
- [ ] Implement `file_stats_aligned(path: Path, seq_type: str | None, input_format) -> dict`:
  - Read alignment with `FormatConverter.read()`
  - Verify all sequences same length (is_aligned=True)
  - Detect seq_type if not provided
  - Check stop codons
  - Compute `per_taxon_stats` for each sequence
  - Call `compute_site_patterns`
  - Aggregate gap/ambiguous/gap+ambiguous/standard ratios (mean across taxa)
  - Return full dict including `alignment_length`, site pattern counts/ratios, `per_taxon`
- [ ] Write tests:
  - [ ] `test_compute_site_patterns_basic` — small hand-crafted alignment
  - [ ] `test_file_stats_aligned` using `test/EOG090X0971.faa`:
    - alignment_length = 1042
    - distinct_patterns = 624
    - constant_sites = 482
    - parsimony_informative = 87
    - singleton_sites = 473
  - [ ] `test_file_stats_aligned_nt` using `test/EOG090X0971.fna`

---

## Task 5: Unified single-file entry point

**File:** `phyloai/pretree/stats.py`

**Steps:**
- [ ] Implement `stats_single_file(path: Path, seq_type: str | None = None, input_format=None) -> dict`:
  - Detect format via `FormatConverter.detect()`
  - Read sequences to detect is_aligned
  - Dispatch to `file_stats_aligned` or `file_stats_unaligned`
  - Return unified dict with all fields (aligned-only fields absent when unaligned)
- [ ] Write tests:
  - [ ] `test_stats_single_file_aligned`
  - [ ] `test_stats_single_file_unaligned`
  - [ ] `test_stats_single_file_format_override`

---

## Task 6: Directory mode — per-file stats worker + parallel execution

**File:** `phyloai/pretree/stats.py`

**Steps:**
- [ ] Implement module-level `_worker(args: tuple) -> dict`:
  - Unpacks `(path, seq_type, input_format)` — must be module-level for `ProcessPoolExecutor` pickling
  - Calls `stats_single_file(path, seq_type, input_format)`
  - Catches exceptions: on error returns `{gene: path.stem, error: str(e)}`
- [ ] Implement `collect_seq_files(directory: Path) -> list[Path]`:
  - Glob for all files with extensions in `COMMON_ALIGNMENT_EXTENSIONS` from `core/schema.py`
  - Sort by filename for deterministic ordering
- [ ] Implement `stats_directory(directory: Path, seq_type: str | None, input_format, threads: int, progress_callback=None) -> tuple[list[dict], list[str]]`:
  - Returns `(per_file_results, warnings)`
  - `threads=1`: run serially (no executor) for debuggability
  - `threads>1`: use `ProcessPoolExecutor(max_workers=threads)` with `map(_worker, args_list)`
  - If `progress_callback` is provided, call it once for each completed file path
  - Collect errors from workers; include failed files in results with `error` field
  - Collect all warnings across files
- [ ] Write tests:
  - [ ] `test_collect_seq_files` — verify correct file discovery
  - [ ] `test_stats_directory_serial` using `3-align/faa/` (small subset, 5 files)
  - [ ] `test_stats_directory_parallel` same subset with threads=2
  - [ ] `test_stats_directory_error_handling` — one malformed file in otherwise valid directory

---

## Task 7: Directory mode — summary aggregation

**File:** `phyloai/pretree/stats.py`

**Steps:**
- [ ] Implement `aggregate_summary(per_file_results: list[dict]) -> dict`:
  - Skip results with `error` field for numeric aggregation; count them separately as `n_errors`
  - Compute:
    - `n_genes` (total files including errors)
    - `n_genes_ok` (successfully processed)
    - `format` — "mixed" if inconsistent across files
    - `is_aligned` — True/False/"mixed"
    - `seq_type` — "AA"/"NT"/"mixed"
    - `total_taxa` — count of unique taxon names across all files
    - `taxa_per_gene`: min/max/mean/median of `n_taxa`
    - `length`: min/max/mean/median/total of per-file length values
    - `gap_ratio`: mean/median across files
    - `ambiguous_ratio`: mean/median across files
    - `gap_ambiguous_ratio`: mean/median across files
    - `missing_taxa_ratio`: mean/median across files (requires knowing max_taxa)
    - `n_taxa_ratio`: derived from max_taxa detected in dataset
  - `max_taxa` = maximum `n_taxa` across all files (used for ratio computations)
- [ ] Write tests:
  - [ ] `test_aggregate_summary_aligned`
  - [ ] `test_aggregate_summary_mixed_alignment_status`
  - [ ] `test_aggregate_summary_with_errors`

---

## Task 8: Output rendering — Rich terminal display

**File:** `phyloai/pretree/stats.py` (rendering helpers) and `phyloai/cli/commands/pretree.py` (display calls)

**Steps:**
- [ ] Implement `render_summary_table(summary: dict) -> rich.table.Table`:
  - Two-column table: Metric / Value
  - Sections: Dataset Overview, Sequence Length, Gap & Ambiguous, Taxa Coverage
- [ ] Implement `render_per_gene_table(per_file: list[dict]) -> rich.table.Table`:
  - Columns: gene, n_taxa, n_taxa_ratio, length_type, alignment_length, seq_length_min, seq_length_max, seq_length_mean, seq_length_median, seq_length_stdev, gap_ratio, ambiguous_ratio, gap_ambiguous_ratio, missing_taxa, missing_taxa_ratio
  - Rows sorted by gene name
- [ ] Implement `render_single_file_panels(stats: dict) -> list[rich.panel.Panel]`:
  - Panel 1: Overview (filename, format, seq_type, is_aligned, n_taxa)
  - Panel 2: Character summary (standard/gap/ambiguous ratios)
  - Panel 3: Site patterns (aligned only) — include MSA length plus distinct-patterns/constant/PIS/singleton counts and ratios
  - Panel 4: Per-taxon table — name, raw_length, ungapped_length, gap_ratio, ambiguous_ratio
- [ ] Write visual smoke test (not automated):
  - [ ] Run `phyloai pretree stats --seq ref/phylogenomics_examples/test/EOG090X0971.faa` and verify output looks correct

---

## Task 9: Output serialisation — result.json plus per-gene table writing

**File:** `phyloai/pretree/stats.py`

**Steps:**
- [ ] Implement `write_output(data: dict, path: Path, mode: str, force_json: bool = False)` so the main command result is always written as structured JSON matching design spec Section 7.3
- [ ] Implement `per_gene_output_path(output_dir, output_format="csv")` to place `per-gene.csv` or `per-gene.tsv` inside `--output-dir`
- [ ] Implement `write_per_gene_output(data, path)` to write the optional per-gene CSV/TSV table
- [ ] Write tests:
  - [ ] `test_write_json_output`
  - [ ] `test_write_per_gene_csv_output`
  - [ ] `test_write_per_gene_tsv_output`

---

## Task 10: CLI command wiring

**File:** `phyloai/cli/commands/pretree.py`

**Steps:**
- [ ] Add `stats` subcommand to `pretree` group with all parameters from design spec Section 6:
  - `--seq-dir` (Path, mutually exclusive with `--seq`)
  - `--seq` (Path, mutually exclusive with `--seq-dir`)
  - `--per-gene` (flag)
  - `--per-gene-format` (choice: csv/tsv, default csv)
  - `--output-dir` / `-o` (Path)
  - `--input-format` (choice: `fasta`, `phylip-relaxed`, `nexus`; optional)
  - `--seq-type` (choice: AA/NT, optional)
  - `--threads` / `-t` (int, default 4)
  - `--quiet` / `-q` (flag)
- [ ] Implement mutual exclusivity check for `--seq-dir` / `--seq`: if both or neither provided, exit code 1 with clear message
- [ ] Implement command body:
  - Dispatch to `stats_single_file` or `stats_directory` + `aggregate_summary`
  - Render terminal output via Rich unless `--quiet` is set
  - Always write structured results to `result.json` inside `--output-dir`
  - When `--per-gene` is provided, write per-gene results under `--output-dir` with CSV default and TSV override via `--per-gene-format`
  - Print saved output paths in terminal output with content-aware wording (`Results saved to ...`, `Per-gene table saved to ...`, etc.)
  - Show a Rich progress bar during directory processing, suppressed under `--quiet`
  - Help text must explicitly state that exactly one of `--seq` or `--seq-dir` is required
  - Exit code 0 on success; 1 on input error; 2 on processing failure (all files errored)
- [ ] Write CLI integration tests:
  - [ ] `test_cli_stats_seq_single_file` — invoke via Click test runner
  - [ ] `test_cli_stats_seq_dir` — small directory subset
  - [ ] `test_cli_stats_writes_result_json`
  - [ ] `test_directory_per_gene_defaults_to_csv`
  - [ ] `test_directory_per_gene_format_can_write_tsv`
  - [ ] `test_unaligned_per_gene_output_includes_sequence_length_summary`
  - [ ] `test_single_file_result_json_includes_per_taxon_table`
  - [ ] `test_cli_stats_mutual_exclusivity_error`
  - [ ] `test_cli_stats_no_args_error`

---

## Task 11: End-to-end verification

**Steps:**
- [ ] Run full test suite: `pytest tests/pretree/ -v`
- [ ] Manual smoke tests:
  - [ ] `phyloai pretree stats --seq ref/phylogenomics_examples/test/EOG090X0971.faa`
  - [ ] `phyloai pretree stats --seq ref/phylogenomics_examples/test/EOG090X0971.fna`
  - [ ] `phyloai pretree stats --seq ref/phylogenomics_examples/2-loci_filter/faa/EOG090X0007.faa` (unaligned)
  - [ ] `phyloai pretree stats --seq-dir ref/phylogenomics_examples/3-align/faa/ --per-gene --output-dir /tmp/stats`
  - [ ] `phyloai pretree stats --seq-dir ref/phylogenomics_examples/2-loci_filter/faa/ --threads 8`
  - [ ] `phyloai pretree stats --seq-dir ref/phylogenomics_examples/2-loci_filter/fna --output-dir /tmp/stats-fna --threads 8 --per-gene` and verify `/tmp/stats-fna/per-gene.csv` contains length summary columns
  - [ ] `phyloai pretree stats --seq ref/phylogenomics_examples/2-loci_filter/fna/EOG090X0971.fna --output-dir /tmp/stats-single` and verify `per_taxon` is saved in `result.json`
- [ ] Verify site pattern values for `EOG090X0971.faa` match expected: distinct_patterns=624, constant=482, PIS=87, singleton=473
- [ ] Verify site pattern values for `raw.fa` match IQ-TREE distinct patterns after missing-state normalization: distinct_patterns=1053473
- [ ] Verify parallel and serial results are identical for same input directory
- [ ] Check all exit codes: success=0, missing file=1, bad format=1
