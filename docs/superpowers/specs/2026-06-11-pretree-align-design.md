# pretree align - Design Specification

**Date:** 2026-06-11
**Status:** Draft for user review
**Parent spec:** `2026-06-07-phyloai-design.md`

---

## 1. Purpose

`phyloai pretree align` aligns a directory of unaligned sequence files using MAFFT or MAGUS. It supports amino-acid sequences, nucleotide sequences, and a combined AA-alignment + backtranslation mode that produces both protein and codon-level alignments from a single run.

It is the first pipeline command in the `pretree` phase. Its outputs feed directly into `phyloai pretree trim`.

Recommended workflow:

```bash
# AA only
phyloai pretree align --seq-dir ./raw_aa --method linsi --seq-type AA \
  --output-dir ./runs/run001/pretree/align

# AA + codon backtranslation
phyloai pretree align --seq-dir ./raw_aa --method linsi --seq-type AA \
  --backtrans --nt-dir ./raw_nt \
  --output-dir ./runs/run001/pretree/align

# NT only
phyloai pretree align --seq-dir ./raw_nt --method linsi --seq-type NT \
  --output-dir ./runs/run001/pretree/align
```

---

## 2. Operating Model

### 2.1 Input scanning

`--seq-dir` is scanned one level deep. Files with extensions `.fa`, `.fas`, `.fasta`, `.faa`, `.fna` are treated as input sequences. Subdirectories, empty files, and unrecognized extensions are silently skipped and recorded in `data.skipped`.

`--input-format` overrides format detection for all input files (default: `auto`). Format detection follows `core/formats.py`.

### 2.2 Three operating modes

**Mode 1 — AA only** (`--seq-type AA`, no `--backtrans`):

```
seq-dir/gene*.fa  →  [MAFFT or MAGUS]  →  seqs/gene*.fa
```

**Mode 2 — NT only** (`--seq-type NT`, no `--backtrans`):

```
seq-dir/gene*.fa  →  [MAFFT or MAGUS]  →  seqs/gene*.fa
```

**Mode 3 — AA + backtrans** (`--seq-type AA`, `--backtrans`, `--nt-dir`):

```
seq-dir/gene*.fa  (AA)      →  [MAFFT or MAGUS]  →  seqs/faa/gene*.fa
nt-dir/gene*.fa   (NT CDS)  ─┘  trimAl -backtrans →  seqs/fna/gene*.fa
```

### 2.3 File pairing in backtrans mode

AA and CDS files are matched by stem (filename without extension). For a given AA file `gene1.fa`, PhyloAI searches `--nt-dir` for any file whose stem is `gene1` regardless of extension. If no matching CDS file is found, the gene is aligned (AA result written to `seqs/faa/`) but backtranslation is skipped with a warning.

### 2.4 CDS validation (backtrans mode only)

Before invoking trimAl, PhyloAI performs lightweight pre-checks on the CDS sequences. These checks produce per-gene warnings and skip backtranslation for affected genes; the AA alignment is always retained.

| Check | Executed by | Disposition |
|-------|-------------|-------------|
| Sequence length not a multiple of 3 | PhyloAI pre-check | Skip backtrans, record warning `CDS length not multiple of 3` |
| Taxon count mismatch between AA alignment and CDS file | PhyloAI pre-check | Skip backtrans, record warning `taxon count mismatch` |
| Trailing stop codon (TAA/TAG/TGA at sequence end) | trimAl (`-ignorestopcodon`) | Auto-handled by trimAl; no warning emitted unless trimAl reports one |
| Internal stop codon | trimAl detection | trimAl exits non-zero; PhyloAI captures stderr, records as warning, skips backtrans for that gene |

PhyloAI always passes `-ignorestopcodon` to trimAl for backtranslation. If trimAl exits with a non-zero code for any other reason, its stderr is captured verbatim as the warning message.

### 2.5 Parallelism

- `concurrent.futures.ProcessPoolExecutor(max_workers=threads)` runs alignment tasks in parallel
- Each task covers the full pipeline for one gene: alignment + optional backtranslation
- Each task uses a single thread for the underlying tool (`--thread 1` for MAFFT; `-np 1` for MAGUS)
- `--threads` controls the number of concurrent tasks (default: 4)
- MAGUS tasks each use an independent `tempfile.mkdtemp()` working directory (`-d`), cleaned up after task completion regardless of success or failure

### 2.6 Partial failure policy

- At least one gene aligned successfully → exit 0
- All genes failed or skipped → exit 1
- Each failure recorded in `data.skipped` with `path` and `reason`

---

## 3. Alignment Methods

`--method` selects the alignment strategy. Default is `linsi`.

| Value | Command | Notes |
|-------|---------|-------|
| `fftns1` | `mafft --retree 1 --thread 1` | Fastest, lowest accuracy |
| `fftns2` | `mafft --retree 2 --thread 1` | Fast, slightly lower accuracy |
| `auto` | `mafft --auto --thread 1` | MAFFT auto-selects strategy |
| `linsi` | `mafft --maxiterate 1000 --localpair --thread 1` | **Default.** High accuracy, recommended for most datasets |
| `einsi` | `mafft --maxiterate 1000 --genafpair --thread 1` | High accuracy, suitable for sequences with non-homologous regions |
| `ginsi` | `mafft --maxiterate 1000 --globalpair --thread 1` | High accuracy, global alignment |
| `magus` | `magus -i ... -o ... -d ... --datatype protein\|dna` | Highest accuracy, slowest; suited for large or difficult datasets |

For MAGUS, `--seq-type AA` maps to `--datatype protein`; `--seq-type NT` maps to `--datatype dna`.

`--extra-args` is only meaningful for `--method magus`. If passed with any MAFFT method, a warning is printed and the argument is ignored. The `--extra-args` string is tokenized with `shlex.split()` and appended to the MAGUS command. Conflicting parameters in `--extra-args` override PhyloAI's internal defaults, following the merge semantics in main design Section 9.8. The fully merged command is written to the log before execution.

---

## 4. Output Directory Structure

### 4.1 Mode 1 and Mode 2 (single sequence type)

```
runs/run001/pretree/align/
├── seqs/
│   ├── gene1.fa
│   ├── gene2.fa
│   └── ...
├── align.log
└── result.json
```

### 4.2 Mode 3 (AA + backtrans)

```
runs/run001/pretree/align/
├── seqs/
│   ├── faa/
│   │   ├── gene1.fa
│   │   └── ...
│   └── fna/
│       ├── gene1.fa
│       └── ...
├── align.log
└── result.json
```

All output files use `.fa` suffix regardless of input suffix or sequence type.

### 4.3 Output directory conflict policy

Follows main design Section 9.5: if `--output-dir` exists and is non-empty, exit 1 with a clear message. `--overwrite` deletes and recreates the output directory before running.

---

## 5. Logging

`align` is a pipeline command and writes a log file `align.log` to the output directory (alongside `result.json`). This deviates from the main design's `runs/runNNN/logs/` convention; the rationale is that the log is most useful when co-located with the outputs it describes.

Each gene's log entry appended to `align.log` contains:
- Tool name and version
- Full resolved command (after `--extra-args` merge)
- stdout and stderr
- Wall time
- Exit code

Log entries are separated by a timestamp header. On `--overwrite` runs the log file is also deleted and recreated.

---

## 6. CLI Parameters

| Parameter | Short | Type | Default | Notes |
|-----------|-------|------|---------|-------|
| `--seq-dir` | | Path | required | Input directory of unaligned sequences |
| `--method` | | `fftns1\|fftns2\|auto\|linsi\|einsi\|ginsi\|magus` | `linsi` | Alignment strategy |
| `--seq-type` | | `AA\|NT` | `AA` | Molecule type |
| `--backtrans` | | flag | `False` | Produce NT codon alignment via trimAl backtranslation |
| `--nt-dir` | | Path | — | Required when `--backtrans` is set; directory of unaligned CDS sequences |
| `--output-dir` | `-o` | Path | `runs/run001/pretree/align` | Output directory |
| `--threads` | `-t` | int | `4` | Number of concurrent alignment tasks |
| `--extra-args` | | str | — | Extra arguments passed to MAGUS only; ignored with warning for MAFFT methods |
| `--input-format` | | str | `auto` | Override format detection for all input files |
| `--overwrite` | | flag | `False` | Overwrite existing non-empty output directory |
| `--dry-run` | | flag | `False` | Print commands without executing; exit 0 |
| `--quiet` | `-q` | flag | `False` | Suppress terminal output except errors |

### 6.1 Parameter constraints (validated before execution)

- `--backtrans` without `--nt-dir` → exit 1
- `--nt-dir` without `--backtrans` → warning only (not an error; `--nt-dir` is ignored)
- `--seq-type NT` with `--backtrans` → exit 1 (backtrans requires AA alignment as input)
- `--threads` < 1 → exit 1

---

## 7. JSON Result Schema

```json
{
  "status": "success | error",
  "command": "phyloai pretree align --seq-dir ./raw_aa --method linsi ...",
  "wall_time": 142.3,
  "tool_versions": {
    "mafft": "7.526",
    "magus": "1.1.0",
    "trimal": "1.4.1"
  },
  "params": {
    "seq_dir": "./raw_aa",
    "method": "linsi",
    "seq_type": "AA",
    "backtrans": false,
    "nt_dir": null,
    "output_dir": "./runs/run001/pretree/align",
    "threads": 4,
    "extra_args": null,
    "input_format": "auto",
    "overwrite": false
  },
  "key_results": {
    "n_aligned": 96,
    "n_skipped": 4,
    "method": "linsi",
    "backtrans": false,
    "mean_alignment_length": 412.3,
    "mean_n_taxa": 48.2
  },
  "error": null,
  "data": {
    "summary": {
      "n_input_files": 100,
      "n_aligned": 96,
      "n_backtrans": 0,
      "n_skipped": 4
    },
    "files": [
      {
        "input": "raw_aa/gene1.fa",
        "output_aa": "seqs/gene1.fa",
        "output_nt": null,
        "n_taxa": 50,
        "alignment_length": 423,
        "wall_time": 1.2,
        "warnings": []
      }
    ],
    "skipped": [
      {
        "path": "raw_aa/gene5.fa",
        "reason": "mafft exited with code 1: <stderr excerpt>"
      }
    ]
  }
}
```

`key_results.mean_alignment_length` is in amino-acid positions for AA/backtrans mode and nucleotide positions for NT mode. `key_results.n_aligned` counts genes with at least one successful alignment output (AA counts in backtrans mode even if NT backtrans was skipped).

`tool_versions` includes only tools actually invoked. For MAFFT methods without backtrans, only `mafft` appears. For backtrans runs, `trimal` is added. For MAGUS, `magus` appears instead of `mafft`.

---

## 8. Library API

```python
# phyloai/pretree/align.py

def run_align(
    seq_dir: Path,
    output_dir: Path,
    method: str,
    seq_type: str,
    backtrans: bool = False,
    nt_dir: Path | None = None,
    threads: int = 4,
    extra_args: str | None = None,
    input_format: str = "auto",
    overwrite: bool = False,
    dry_run: bool = False,
    progress_callback: Callable[[Path], None] | None = None,
) -> dict[str, Any]: ...
```

The return value matches the JSON schema in Section 7. The CLI layer writes this dict to `result.json` and renders terminal output from it. The library function itself does not write `result.json`; that responsibility belongs to the CLI layer, consistent with `convert.py` and `stats.py`.

Internal helpers (not part of the public API):

- `_build_mafft_cmd(input_file, output_file, method)` → `list[str]`
- `_build_magus_cmd(input_file, output_file, work_dir, seq_type, extra_args)` → `list[str]`
- `_align_one(gene_path, output_dir, method, seq_type, extra_args, dry_run)` → per-gene result dict
- `_backtrans_one(aa_aln_path, nt_path, output_nt_path, dry_run)` → per-gene backtrans result dict
- `_validate_cds(nt_sequences)` → list of warning strings
- `_scan_input(seq_dir, input_format)` → list of valid input Paths + skipped list

---

## 9. Environment Requirements

Checked via `core/env.py` (`ToolEnv.require()`) before execution:

| Tool | When required |
|------|--------------|
| `mafft` | All MAFFT methods (`fftns1`, `fftns2`, `auto`, `linsi`, `einsi`, `ginsi`) |
| `magus` | `--method magus` |
| `trimal` | `--backtrans` |

Missing tool → exit 3 (environment error), consistent with main design Section 9.3.

---

## 10. Integration with Adjacent Commands

- **Upstream:** `phyloai pretree convert` or user-provided raw sequences. `pretree align` does not perform format conversion; inputs must be FASTA or another auto-detectable format.
- **Downstream:** `phyloai pretree trim` reads `seqs/` (or `seqs/faa/` and `seqs/fna/` in backtrans mode). The `--msa-dir` argument of `pretree trim` should point to the appropriate subdirectory.
- **Report:** `key_results` supplies alignment count, method, and mean length for the Methods paragraph (e.g., "X genes were aligned using MAFFT L-INS-i; mean alignment length Y aa").

---

## 11. `--dry-run` Behavior

With `--dry-run`:
- PhyloAI resolves all inputs and builds all commands
- Commands are printed in a Rich table (gene | tool | full command)
- No files are created; no tools are executed
- `result.json` is not written
- Exit 0

---

## 12. Terminal Output

Unless `--quiet` is set:
- Rich progress bar during batch alignment (one advance per completed gene)
- Summary table on completion: input files scanned, aligned, backtranslated, skipped, method, seq type, mean alignment length
- Per-gene warnings printed after the summary table
- Explicit save path messages: `Alignments saved to <output_dir>/seqs/` and `Result saved to <output_dir>/result.json`

---

## 13. Updates to Main Design

The main design (`2026-06-07-phyloai-design.md`) requires one correction:

**Section 9.4** example directory structure for `pretree align` currently shows aligned files directly in the top-level output directory:

```
runs/run001/pretree/align/
├── gene1.aln
├── gene2.aln
└── result.json
```

This should be updated to match the unified `seqs/` convention:

```
runs/run001/pretree/align/
├── seqs/
│   ├── gene1.fa
│   └── ...
├── align.log
└── result.json
```

The log file location also deviates from the main design's `runs/runNNN/logs/<step>.log` convention (Section 9.6). This spec places `align.log` inside the output directory. The main design's convention remains the default for other pipeline commands; this deviation is specific to `pretree align` and should be noted in the main design.

---

## 14. Documentation Requirements

Implementation must add or update:

- `docs/commands/pretree-align.md` with the standard sections: Purpose, Usage, Inputs, Outputs, Examples, Warnings and Errors, Notes
- `README.md` command index entry linking to `docs/commands/pretree-align.md`

---

## 15. Testing Requirements

Tests should cover:

- Directory input with mixed valid files, subdirectories, empty files, and unrecognized extensions
- All six MAFFT methods produce correct command strings
- MAGUS command construction with and without `--extra-args`
- `--extra-args` ignored with warning for MAFFT methods
- `--seq-type NT` maps to `--datatype dna` for MAGUS
- Backtrans mode: correct pairing, correct `seqs/faa/` and `seqs/fna/` output layout
- CDS pre-checks: length-not-multiple-of-3 and taxon-count-mismatch skip with warning
- trimAl non-zero exit captured as warning, gene skipped
- `--backtrans` without `--nt-dir` exits 1
- `--seq-type NT` with `--backtrans` exits 1
- Partial failure: some genes fail, exit 0 with skipped list
- All genes fail: exit 1
- `--dry-run` prints commands, creates no files
- Output directory conflict with and without `--overwrite`
- `result.json` schema shape and `key_results` field correctness
- MAGUS temp directory created per gene and cleaned up after task
