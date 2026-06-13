# phyloai pretree trim

## Purpose

`phyloai pretree trim` is PhyloAI's batch MSA-trimming command for already aligned FASTA files. It trims one gene alignment at a time with the selected backend tool (`trimAl`, `BMGE`, or `ClipKIT`) and supports AA-only, NT-only, CODON, and AA+NT dual-output workflows.

This command does not perform sequence alignment or format conversion. Inputs must already be aligned FASTA files. Run `phyloai pretree align` first, and use `phyloai pretree convert --to fasta` beforehand when the source files are not FASTA or when you want an extra normalization and validation pass on FASTA inputs.

## Usage

Minimal:
```bash
phyloai pretree trim --msa-dir ./aligned
```

AA + NT with trimAl manual thresholds:
```bash
phyloai pretree trim \
  --msa-dir ./runs/pretree/align/seqs/faa \
  --nt-dir ./runs/pretree/align/seqs/fna \
  --tool trimal \
  --tool-args "-gt 0.9 -cons 60" \
  --output-dir ./runs/pretree/trim
```

## Parameters

| Parameter | Default | Notes |
|-----------|---------|-------|
| `--msa-dir` | required | Directory of aligned input MSA files |
| `--output-dir` / `-o` | `runs/pretree/trim` | Output directory |
| `--tool` | `trimal` | `trimal`, `bmge`, or `clipkit` |
| `--seq-type` | `auto` | `AA`, `NT`, `CODON`, or `auto`; `auto` only resolves AA vs NT |
| `--nt-dir` | — | AA+NT mode only; trimAl accepts raw CDS or gapped codon-aligned NT because PhyloAI strips NT gaps before backtranslation; BMGE and ClipKIT expect codon-aligned NT MSAs |
| `--trimal-method` | `automated1` | trimAl automatic preset; ignored when `--tool-args` includes manual trimAl thresholds such as `-gt` or `-cons` |
| `--bmge-matrix` | dynamic | AA/CODON default `BLOSUM62`; NT default `DNAPAM100:2` |
| `--bmge-entropy` | `0.5` | Lower is more stringent |
| `--clipkit-method` | `smart-gap` | ClipKIT mode |
| `--trimal-path` | — | Explicit trimAl executable path |
| `--bmge-path` | — | Explicit BMGE.jar path |
| `--clipkit-path` | — | Explicit clipkit executable path |
| `--threads` / `-t` | 4 | Concurrent trimming tasks |
| `--tool-args` | — | Tool strategy args only; PhyloAI manages input/output/log/codon/thread flags |
| `--resume` | off | Resume from `checkpoint.json` |
| `--overwrite` | off | Delete and recreate non-empty output directory |
| `--dry-run` | off | Print commands, create no files |
| `--quiet` / `-q` | off | Suppress terminal output except errors |

## Inputs

Scans `--msa-dir` one level deep for `.fa`, `.fas`, `.fasta`, `.faa`, `.fna`. Subdirectories, empty files, and unrecognized extensions are skipped.

For AA+NT mode:

- `trimAl`: `--msa-dir` is aligned AA MSA; `--nt-dir` supplies matching NT files by stem. The NT input may be raw CDS or gapped codon-aligned NT. PhyloAI removes NT gaps before calling trimAl `-backtrans`.
- `BMGE`: `--msa-dir` is aligned AA MSA; `--nt-dir` is matching codon-aligned NT MSA. BMGE trims AA, then PhyloAI projects kept columns to NT.
- `ClipKIT`: `--msa-dir` is aligned AA MSA; `--nt-dir` is matching codon-aligned NT MSA. ClipKIT trims AA, then PhyloAI projects kept columns to NT.

## Outputs

**AA-only or NT-only:**
```
runs/pretree/trim/
├── seqs/
│   ├── gene1.fa
│   └── ...
├── trim.log
├── checkpoint.json
└── result.json
```

**CODON or AA+NT:**
```
runs/pretree/trim/
├── seqs/
│   ├── faa/
│   │   └── gene1.fa
│   └── fna/
│       └── gene1.fa
├── trim.log
├── checkpoint.json
└── result.json
```

`result.json` reports trimmed/skipped counts, skipped reasons, before/after alignment lengths, resolved parameters, and warnings.

## Examples

```bash
# Default trimAl automatic trimming
phyloai pretree trim --msa-dir ./aligned_aa --tool trimal

# trimAl manual thresholds via --tool-args
phyloai pretree trim --msa-dir ./aligned_aa --tool trimal \
  --tool-args "-gt 0.9 -cons 60"

# BMGE on codon-aligned NT input
phyloai pretree trim --msa-dir ./aligned_codon --seq-type CODON --tool bmge

# ClipKIT AA+NT projection mode
phyloai pretree trim --msa-dir ./aligned_aa --nt-dir ./aligned_nt \
  --tool clipkit --clipkit-method gappy
```

## Warnings and Errors

| Condition | Behaviour |
|-----------|-----------|
| `--seq-type CODON` with `--nt-dir` | Exit 1 |
| Missing external tool | Exit 3 |
| Non-empty output directory without `--overwrite` or `--resume` | Exit 1 |
| No valid input files | Exit 1 |
| Missing NT pair in AA+NT mode | Gene skipped |
| trimAl backtrans receives gapped NT input | PhyloAI strips gaps first and continues |
| Tool exits non-zero | Gene skipped with stderr-derived reason |
| All genes skipped | Exit 2 |

## Notes

- `--tool-args` is strategy-only. Do not pass tool-managed flags such as trimAl `-in/-out/-backtrans`, BMGE `-i/-of/-t`, or ClipKIT `-o/--codon`.
- For trimAl, if `--tool-args` includes manual thresholds like `-gt` or `-cons`, PhyloAI does not add the default automatic preset.
- `trim.log` records the resolved command, status, wall time, reason, and stderr for each gene.
