# pretree trim Design Specification

**Date:** 2026-06-12  
**Status:** Approved for implementation  
**Reference:** `ref/scripts/trimming_alignments.sh`, `docs/superpowers/specs/2026-06-07-phyloai-design.md`

---

## 1. Purpose

`phyloai pretree trim` trims multiple sequence alignments (MSAs) in batch, supporting three tools (trimAl, BMGE, ClipKIT), four sequence modes (AA-only, NT-only, CODON, AA+NT), and parallel execution. It is the step after `pretree align` and before `pretree metrics` in the standard phylogenomics pipeline.

What it does **not** do:
- Sequence alignment (use `pretree align`)
- Format conversion (use `pretree convert`)
- Marker filtering based on metrics (use `pretree filter`)

---

## 2. Architecture

```
phyloai/pretree/trim.py              # core library: tool builders, workers, orchestration
phyloai/cli/commands/pretree.py      # CLI: new trim subcommand registration
```

**Reused components:**
- `core/env.py` — `TOOL_REGISTRY` (trimal, bmge, clipkit already registered), `ToolEnv._detect_tool`, `_resolve_tool_paths` pattern
- `core/sequence_normalization.py` — `detect_seq_type`, extended with `_validate_codon_msa()`
- `core/sequence_output_validation.py` — `validate_fasta_output`
- `core/checkpoint.py` — checkpoint schema and helpers

`pretree/align.py` is referenced for implementation patterns (backtrans, stem pairing, CDS validation) but **not imported** by trim.py; trim.py implements its own versions to avoid cross-module coupling.

---

## 3. Operating Modes

Mode is determined automatically from `--seq-type` and whether `--nt-dir` is provided. No explicit mode flag is required.

### Mode 1 — AA-only
- `--seq-type AA` (or auto-detected), no `--nt-dir`
- Input: AA MSA directory (`--msa-dir`)
- Output: `seqs/` flat — trimmed AA files
- All three tools supported

### Mode 2 — NT-only
- `--seq-type NT` (or auto-detected), no `--nt-dir`
- Input: NT MSA directory (`--msa-dir`)
- Output: `seqs/` flat — trimmed NT files
- All three tools supported; BMGE uses `-t DNA`

### Mode 3 — CODON
- `--seq-type CODON` (must be explicit; cannot be auto-detected), no `--nt-dir`
- Input: codon-aligned NT MSA directory (`--msa-dir`)
- Output: `seqs/faa/` (Python-translated AA) + `seqs/fna/` (trimmed codon)
- All three tools supported; this is also the recommended path for BMGE AA+NT output
- trimAl: Python translate codon MSA → temp AA MSA FASTA; Python strip gaps from codon MSA → temp unaligned CDS FASTA; `trimal -in <temp_aa_msa> -out <faa/gene.fa> -<method> -backtrans <temp_unaligned_cds> -ignorestopcodon`; both temp files in `tempfile.TemporaryDirectory`
- BMGE: `java -jar BMGE.jar -i <codon_msa> -t CODON -m <matrix> -h <entropy> -of <fna/gene.fa>` → Python translate trimmed codon MSA → `faa/`
- ClipKIT: `clipkit <codon_msa> -o <fna/gene.fa> -m <method> --codon` → Python translate → `faa/`

### Mode 4 — AA+NT
- `--seq-type AA` (or auto-detected) + `--nt-dir <path>`
- Input: AA MSA directory (`--msa-dir`) + NT directory (`--nt-dir`)
- Output: `seqs/faa/` + `seqs/fna/`
- **`--nt-dir` content differs by tool** (documented in CLI help):
  - trimAl: unaligned CDS sequences (same as `pretree align --nt-dir`)
  - ClipKIT: codon-aligned NT MSA (must be pre-aligned)
  - BMGE: codon-aligned NT MSA — auto-downgrade to CODON mode (see below)
- trimAl: `trimal -in <aa_msa> -out <faa/gene.fa> -<method> -backtrans <nt_cds> -ignorestopcodon`; AA MSA separately trimmed for `faa/`
- ClipKIT: `clipkit <aa_msa> -o <faa/gene.fa> -m <method> -l` → parse `.log` → Python project kept columns onto codon-aligned NT MSA → `fna/`
- BMGE + Mode 4: **auto-downgrade with warning** — BMGE does not support AA+NT directly; when `--tool bmge --seq-type AA --nt-dir` is specified, a `[WARN]` is emitted and the command automatically uses `--nt-dir` files in CODON mode (`-t CODON`); `params.effective_seq_type` in `result.json` records `CODON`

#### Tool × Mode support matrix

| `--seq-type` | `--nt-dir` | trimAl | BMGE | ClipKIT | Output |
|---|---|---|---|---|---|
| `AA` (or auto) | no | ✓ | ✓ (`-t AA`) | ✓ | `seqs/` flat |
| `NT` (or auto) | no | ✓ | ✓ (`-t DNA`) | ✓ | `seqs/` flat |
| `CODON` | no | ✓ (Python+backtrans) | ✓ (`-t CODON`) | ✓ (`--codon`) | `seqs/faa/` + `seqs/fna/` |
| `AA` (or auto) | unaligned CDS | ✓ (`-backtrans`) | WARN→CODON | ✗ | `seqs/faa/` + `seqs/fna/` |
| `AA` (or auto) | codon-aligned NT MSA | ✗ | WARN→CODON | ✓ (log→Python) | `seqs/faa/` + `seqs/fna/` |

---

## 4. CLI Parameters

### Universal parameters (per main design Section 9.2)

| Parameter | Short | Type | Default | Notes |
|---|---|---|---|---|
| `--msa-dir` | | Path | required | Input MSA directory |
| `--output-dir` | `-o` | Path | `runs/pretree/trim` | Output directory |
| `--seq-type` | | `AA`\|`NT`\|`CODON`\|`auto` | `auto` | `auto` detects AA vs NT only; CODON must be explicit |
| `--nt-dir` | | Path | — | NT directory; content differs by tool — trimAl: unaligned CDS; ClipKIT: codon-aligned NT MSA; BMGE: codon-aligned NT MSA (auto-downgrade to CODON mode) |
| `--tool` | | `trimal`\|`bmge`\|`clipkit` | `trimal` | Trimming tool |
| `--threads` | `-t` | int | 4 | Parallel workers |
| `--extra-args` | | str | — | Extra arguments passed to tool (extra-wins merge semantics per main design Section 9.8) |
| `--dry-run` | | flag | False | Print commands without executing |
| `--quiet` | `-q` | flag | False | Suppress Rich output |
| `--overwrite` | | flag | False | Delete and recreate output directory |
| `--resume` | | flag | False | Resume from `checkpoint.json` |

### trimAl-specific parameters

| Parameter | Type | Default | Notes |
|---|---|---|---|
| `--trimal-method` | `automated1`\|`gappyout`\|`strict`\|`strictplus` | `automated1` | Automated trimming strategy; `automated1` is conservative (recommended for most datasets), `gappyout` is more aggressive |
| `--trimal-path` | Path | — | Override bundled/PATH trimal executable |

### BMGE-specific parameters

| Parameter | Type | Default | Notes |
|---|---|---|---|
| `--bmge-matrix` | str | `BLOSUM62` (AA/CODON); `DNAPAM100:2` (NT) | `-m` substitution matrix; AA/CODON: BLOSUM30–BLOSUM95 (higher = stricter); NT: DNAPAMx:y (lower first number = stricter); default is dynamic based on `--seq-type` |
| `--bmge-entropy` | float | `0.5` | `-h` entropy cutoff; lower = stricter; 0.2–0.4 is considered stringent |
| `--bmge-path` | Path | — | Override bundled BMGE.jar path |

### ClipKIT-specific parameters

| Parameter | Type | Default | Notes |
|---|---|---|---|
| `--clipkit-method` | str | `smart-gap` | ClipKIT trimming mode; common values: `smart-gap` (recommended, dynamic gap threshold), `kpi-smart-gap` (keeps parsimony-informative sites + smart-gap), `kpic-smart-gap` (keeps parsimony-informative and constant sites + smart-gap), `gappy`, `kpi-gappy`, `kpic-gappy`, `kpi`, `kpic`; full list and descriptions in `clipkit -h` |
| `--clipkit-path` | Path | — | Override PATH clipkit executable |

---

## 5. Core Logic

### File scanning

- Accepted extensions: `.fa`, `.fas`, `.fasta`, `.faa`, `.fna`
- Skips: subdirectories, empty files, unrecognized extensions
- Files sorted by name for reproducible ordering
- `--nt-dir` scanned with same rules; paired to `--msa-dir` files by stem (`gene1.faa` pairs with `gene1.fna`)
- No valid input files → exit 1

### Tool execution — per tool

**trimAl:**
- Mode 1/2: `trimal -in <msa> -out <out> -<method>`
- Mode 3 (CODON): Python translate codon MSA → temp AA FASTA (with gaps, normal MSA); Python strip gaps from codon MSA → temp unaligned CDS FASTA; then `trimal -in <temp_aa_msa> -out <faa/gene.fa> -<method> -backtrans <temp_unaligned_cds> -ignorestopcodon`; both temp files managed via `tempfile.TemporaryDirectory`
- Mode 4 (AA + `--nt-dir`): `trimal -in <aa_msa> -out <faa/gene.fa> -<method> -backtrans <nt_cds> -ignorestopcodon`; AA MSA also separately trimmed → `faa/`

**BMGE:**
- Mode 1 (AA): `java -jar BMGE.jar -i <msa> -t AA -m <matrix> -h <entropy> -of <out>`
- Mode 2 (NT): `java -jar BMGE.jar -i <msa> -t DNA -m <matrix> -h <entropy> -of <out>`
- Mode 3 (CODON): `java -jar BMGE.jar -i <codon_msa> -t CODON -m <matrix> -h <entropy> -of <fna/gene.fa>` → Python translate trimmed codon MSA → `faa/`
- Mode 4 auto-downgrade: emit `[WARN]`, switch to CODON logic using `--nt-dir` files

**ClipKIT:**
- Mode 1/2: `clipkit <msa> -o <out> -m <method>`
- Mode 3 (CODON): `clipkit <codon_msa> -o <fna/gene.fa> -m <method> --codon` → Python translate trimmed codon MSA → `faa/`
- Mode 4 (AA + `--nt-dir`): `clipkit <aa_msa> -o <faa/gene.fa> -m <method> -l` → parse `.log` → Python project kept columns onto codon-aligned NT MSA → `fna/`

### ClipKIT log parsing (Mode 4 only)

ClipKIT `-l` generates `<output>.log`, a space-separated file where column 1 is the **1-based** site index and column 2 is `trim` or `keep`. Empirically confirmed behavior:

- When run on an AA MSA (no `--codon`): each log row = one AA column
- Log column indices are 1-based; Python projection uses 0-based

```python
# Parse log: column 1 = 1-based site index, column 2 = trim/keep
log_rows = [line.strip().split() for line in open(log_path) if line.strip()]
kept_aa_cols = [int(row[0]) - 1 for row in log_rows if row[1] == "keep"]

# Project onto codon-aligned NT MSA: each AA col i -> NT cols i*3, i*3+1, i*3+2
nt_kept = []
for i in kept_aa_cols:
    nt_kept.extend([i*3, i*3+1, i*3+2])

trimmed_seq = "".join(str(record.seq)[j] for j in nt_kept)
```

This projection was empirically verified: output matches direct `clipkit --codon` trimming on the same codon MSA.

`.log` files are deleted after processing; they are not written to the output directory.

### CODON validation (`_validate_codon_msa` in `core/sequence_normalization.py`)

Applied to all inputs when `--seq-type CODON` or when BMGE/ClipKIT/trimAl process a codon MSA:

| Issue | Handling |
|---|---|
| Alignment length not divisible by 3 | Skip gene; reason: `codon_length_not_multiple_of_3` |
| Internal stop codons (TAA/TAG/TGA before last codon) | Warning + continue (lenient, consistent with trimAl `-ignorestopcodon`) |
| Terminal stop codon | Remove last codon before processing |
| Gap-only codon columns (`---`) | Treated as `-` in translation; no special handling needed |

This function reuses the stop codon detection logic from the existing `_validate_cds` in `pretree/align.py`, adapted for MSA input (gap-aware).

### Temporary file management

Workers that require intermediate files (trimAl CODON mode temp AA FASTA; ClipKIT `.log` files) use `tempfile.TemporaryDirectory` scoped to the worker function. All temporary files are cleaned up automatically when the worker exits, regardless of success or failure.

### `--seq-type auto` detection

When `--seq-type auto` (the default):
- Call `detect_seq_type` from `core/sequence_normalization.py` (samples up to 3 files, 10 sequences)
- Resolves to `AA` or `NT` only; `CODON` is never auto-detected
- If resolved type is unexpected, emit a warning and continue
- Detection result is recorded in `params.effective_seq_type` in `result.json`

### Alignment length statistics

`length_before` and `length_after` in `key_results` refer to **alignment column count** (including gap columns), not ungapped sequence length. This is the standard measure for trimming effectiveness.

---

## 6. Parallelism

Follows `pretree align` exactly:

- `ProcessPoolExecutor(max_workers=threads)`
- Worker function (`_trim_one_worker`) must be a **module-level function** (required for macOS spawn semantics with `ProcessPoolExecutor`)
- Each gene is one atomic task: AA trim + NT output (backtrans / projection / translation) are all handled within a single worker call
- Results collected via `as_completed`
- Checkpoint flushed at most every `CHECKPOINT_FLUSH_INTERVAL = 2.0` seconds (throttled, not per-gene)
- `KeyboardInterrupt`: mark checkpoint status as `interrupted` → force flush → re-raise

---

## 7. Parameter Validation

All checks run before any file processing. Failures exit immediately.

| Check | Condition | Behavior |
|---|---|---|
| `--seq-type CODON` + `--nt-dir` | CODON mode does not use `--nt-dir` | exit 1: error |
| BMGE + `--seq-type AA` + `--nt-dir` | BMGE cannot do AA+NT directly | `[WARN]` + auto-downgrade to CODON mode using `--nt-dir` |
| `--overwrite` + `--resume` | Mutually exclusive | exit 1: error |
| `--resume` without `checkpoint.json` | File not found | exit 1: error |
| `--resume` params mismatch | tool / seq-type / method changed vs checkpoint | exit 1: error |
| `--threads < 1` | Invalid value | exit 1: error |
| `--trimal-path` / `--bmge-path` / `--clipkit-path` exists | File does not exist | exit 1: error |
| Tool not available | Not in PATH, not bundled, no explicit path | exit 3: environment error |
| `--msa-dir` missing or empty | No valid input files after scanning | exit 1: error |
| NT pairing missing (Mode 4) | Gene in `--msa-dir` has no match in `--nt-dir` | Record in `data.skipped`; continue with remaining genes |

### Per-gene error handling

- Tool exits non-zero → gene skipped; stderr recorded in `result.json` and `trim.log`
- Output file validation fails (empty alignment, zero records, unequal lengths) → gene skipped
- Trimming removes all columns (output sequence length = 0) → gene skipped; reason: `all_columns_trimmed`
- All genes skipped → exit 2
- Partial success (some genes skipped) → exit 0; full detail in `result.json`

---

## 8. Checkpoint Design

Follows `docs/superpowers/specs/2026-06-12-checkpoint-resume-design.md` exactly.

```json
{
  "command": "pretree trim",
  "status": "running | completed | interrupted",
  "params": {
    "msa_dir": "./aligned",
    "nt_dir": null,
    "seq_type": "AA",
    "effective_seq_type": "AA",
    "tool": "trimal",
    "trimal_method": "automated1",
    "bmge_matrix": null,
    "bmge_entropy": null,
    "clipkit_method": null,
    "threads": 4,
    "extra_args": null
  },
  "tasks": {
    "gene1": {
      "status": "success",
      "outputs": ["seqs/faa/gene1.fa", "seqs/fna/gene1.fa"]
    },
    "gene2": {
      "status": "failed",
      "reason": "trimal exited with code 1: ..."
    },
    "gene3": {
      "status": "pending"
    }
  }
}
```

`--resume` skips tasks with `status: success`; reruns `failed`, `pending`, `interrupted`.

---

## 9. Output Directory Structure

**Mode 1/2 (AA-only or NT-only):**
```
runs/pretree/trim/
├── seqs/
│   ├── gene1.fa
│   ├── gene2.fa
│   └── ...
├── trim.log
├── checkpoint.json
└── result.json
```

**Mode 3/4 (CODON or AA+NT dual output):**
```
runs/pretree/trim/
├── seqs/
│   ├── faa/
│   │   ├── gene1.fa
│   │   └── ...
│   └── fna/
│       ├── gene1.fa
│       └── ...
├── trim.log
├── checkpoint.json
└── result.json
```

All output files use `.fa` suffix regardless of input suffix, consistent with `pretree align`.

---

## 10. Logging

`trim.log` content per gene entry:
- Fully merged command (including `--extra-args` tokens, after extra-wins merge)
- Tool version
- Wall time
- Exit code
- Full stderr
- Stdout only when diagnostic (trimAl/BMGE/ClipKIT primary output is always a file)

On `--resume`: append to existing log with timestamp separator:
```
=== RESUME 2026-06-12T14:32:01 ===
```

On `--overwrite`: log file deleted and recreated with the output directory.

---

## 11. result.json Schema

```json
{
  "status": "success | error",
  "command": "phyloai pretree trim --msa-dir ./aligned --tool trimal ...",
  "wall_time": 87.4,
  "tool_versions": {
    "trimal": "1.4.1"
  },
  "params": {
    "msa_dir": "./aligned",
    "nt_dir": null,
    "seq_type": "auto",
    "effective_seq_type": "AA",
    "tool": "trimal",
    "trimal_method": "automated1",
    "bmge_matrix": null,
    "bmge_entropy": null,
    "clipkit_method": null,
    "threads": 4,
    "extra_args": null,
    "output_dir": "runs/pretree/trim"
  },
  "key_results": {
    "total_genes": 350,
    "trimmed_genes": 342,
    "skipped_genes": 8,
    "skipped_reasons": {
      "all_columns_trimmed": 6,
      "validation_failed": 2
    },
    "length_before": {"mean": 412.3, "min": 89, "max": 1203},
    "length_after":  {"mean": 287.1, "min": 54, "max": 891},
    "columns_removed_pct": {"mean": 30.4, "min": 5.2, "max": 78.1}
  },
  "error": null,
  "data": {
    "mode": "AA-only",
    "skipped": [
      {"gene": "gene42", "reason": "all_columns_trimmed"},
      {"gene": "gene87", "reason": "nt_pairing_missing"}
    ],
    "warnings": [
      "BMGE does not support AA+NT mode directly; automatically switched to CODON mode using --nt-dir"
    ],
    "per_gene": [
      {
        "gene": "gene1",
        "length_before": 412,
        "length_after": 287,
        "columns_removed": 125,
        "outputs": ["seqs/faa/gene1.fa"]
      }
    ]
  }
}
```

`length_before` / `length_after` / `columns_removed` are **alignment column counts** (including gap columns).

---

## 12. doctor Integration

`phyloai doctor` already registers trimal, bmge (BMGE.jar), and clipkit in `TOOL_REGISTRY`. No changes required. Tool path resolution in `pretree trim` follows the same `ToolEnv` priority:
1. Explicit `--trimal-path` / `--bmge-path` / `--clipkit-path`
2. Bundled path (`phyloai/bundled/`)
3. `shutil.which()` on PATH

---

## 13. Key Design Decisions

| Decision | Rationale |
|---|---|
| Single tool per run (`--tool`) | Consistent with `pretree align --method`; users run multiple tools by running command multiple times with different `--output-dir` |
| No `--seq-type auto` for CODON | CODON sequences are NT-based; cannot be distinguished from NT by character inspection alone |
| BMGE AA+NT via CODON auto-downgrade | BMGE `-t CODON` is the correct mechanism; warning + auto-downgrade avoids user confusion without blocking the use case |
| ClipKIT NT via Python log parsing (not PhyKIT) | Avoids PhyKIT dependency; log column extraction is simple and reliable with BioPython |
| trimAl CODON via Python translate + `-backtrans` | Single trimal call handles both AA and NT output; simpler than column-index projection |
| `_validate_codon_msa` in `core/sequence_normalization.py` | Reuses stop codon detection from `_validate_cds`; shared location ensures consistency |
| Temporary files in `tempfile.TemporaryDirectory` | Auto-cleanup on worker exit; no residual files in output directory regardless of failure |
| `effective_seq_type` in `params` | Captures BMGE auto-downgrade transparently; ensures `result.json` accurately reflects what was executed |
| `length_before/after` = alignment column count | Standard measure for trimming; ungapped length varies per taxon and is not meaningful for MSA-level statistics |
| Checkpoint per gene | Atomic unit matches parallelism model; consistent with `pretree align` resume design |
