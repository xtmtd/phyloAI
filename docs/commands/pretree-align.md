# phyloai pretree align

## Purpose

Align a directory of unaligned sequence files using MAFFT or MAGUS. Produces one aligned FASTA file per input gene. Optionally produces codon-level NT alignments via trimAl backtranslation.

This command does not perform format conversion. Inputs must be FASTA. Run `phyloai pretree convert --to fasta` first for PHYLIP, Nexus, or other formats.

## Usage

Minimal:
```bash
phyloai pretree align --seq-dir ./raw_aa
```

Full:
```bash
phyloai pretree align \
  --seq-dir ./raw_aa \
  --method linsi \
  --seq-type AA \
  --output-dir ./runs/pretree/align \
  --threads 4
```

With backtranslation:
```bash
phyloai pretree align \
  --seq-dir ./raw_aa \
  --method linsi \
  --seq-type AA \
  --backtrans \
  --nt-dir ./raw_nt \
  --output-dir ./runs/pretree/align \
  --threads 4
```

## Parameters

| Parameter | Default | Notes |
|-----------|---------|-------|
| `--seq-dir` | required | Directory of unaligned sequence files |
| `--method` | `linsi` | fftns1, fftns2, auto, linsi, einsi, ginsi, magus |
| `--seq-type` | `auto` | AA, NT, or auto (auto-detects from first few genes) |
| `--backtrans` | off | Produce NT codon alignment; requires --nt-dir |
| `--nt-dir` | — | Unaligned CDS directory for backtrans |
| `--output-dir` / `-o` | `runs/pretree/align` | Output directory |
| `--threads` / `-t` | 4 | Concurrent alignment tasks (each uses 1 thread) |
| `--extra-args` | — | Extra args for MAGUS only; ignored for MAFFT methods |
| `--mafft-path` | — | Explicit MAFFT executable path for MAFFT methods |
| `--magus-path` | — | Explicit MAGUS executable path for `--method magus` |
| `--trimal-path` | — | Explicit trimAl executable path for `--backtrans` |
| `--resume` | off | Resume from `checkpoint.json`; requires exact same resolved parameters |
| `--overwrite` | off | Delete and recreate non-empty output directory |
| `--dry-run` | off | Print commands, create no files |
| `--quiet` / `-q` | off | Suppress terminal output except errors |

## Inputs

Scans `--seq-dir` one level deep for files with extensions: `.fa`, `.fas`, `.fasta`, `.faa`, `.fna`. Subdirectories, empty files, and unrecognized extensions are skipped.

## Outputs

**Mode AA or NT only:**
```
runs/pretree/align/
├── seqs/
│   ├── gene1.fa
│   └── ...
├── align.log
├── checkpoint.json
└── result.json
```

**Mode AA + backtrans:**
```
runs/pretree/align/
├── seqs/
│   ├── faa/
│   │   └── gene1.fa
│   └── fna/
│       └── gene1.fa
├── align.log
├── checkpoint.json
└── result.json
```

`result.json` contains `key_results` with `n_aligned`, `method`, `mean_alignment_length`, and `mean_n_taxa` for report integration.

`align.log` records commands, timing, exit status, and stderr. MAFFT alignment stdout is saved as FASTA files under `seqs/` and is not duplicated in the log.

## Examples

```bash
# Fast alignment for large dataset
phyloai pretree align --seq-dir ./raw_aa --method fftns2 --threads 8

# High-accuracy protein alignment + codon NT alignment
phyloai pretree align --seq-dir ./raw_aa --seq-type AA \
  --backtrans --nt-dir ./raw_nt --method linsi --threads 4

# NT direct alignment
phyloai pretree align --seq-dir ./raw_nt --seq-type NT --method linsi

# MAGUS with extra options
phyloai pretree align --seq-dir ./raw_aa --method magus \
  --extra-args "--maxsubsetsize 200" --threads 4

# Preview commands without running
phyloai pretree align --seq-dir ./raw_aa --method linsi --dry-run

# Resume an interrupted run
phyloai pretree align --seq-dir ./raw_aa --method linsi --seq-type AA \
  --output-dir ./runs/pretree/align --resume
```

## Resume behavior

`pretree align` supports `--resume` to recover from interruption, power loss, or external tool failure without redoing completed work.

- The output directory must already contain `checkpoint.json`.
- The current invocation's resolved parameters must match the checkpoint exactly. This includes analysis parameters and run-control parameters such as `--threads` and `--quiet`.
- Tasks with status `success` and still-valid output files are skipped.
- Tasks with status `failed`, `pending`, `running`, or `success` whose outputs are now missing or invalid are rerun.
- `--resume` and `--overwrite` are mutually exclusive.
- Resume appends to `align.log` and rewrites `result.json` on completion.

## Warnings and Errors

| Condition | Behaviour |
|-----------|-----------|
| `--backtrans` without `--nt-dir` | Exit 1 |
| `--seq-type NT` with `--backtrans` | Exit 1 |
| `--seq-type auto` detects NT with `--backtrans` | Exit 1 |
| `mafft` or `magus` not found | Exit 3 |
| `trimal` not found with `--backtrans` | Exit 3 |
| `--method magus` on non-Linux | Exit 1 (MAGUS bundled binaries are Linux-only) |
| Non-empty output directory | Exit 1 (use `--overwrite`) |
| `--resume` without checkpoint | Exit 1 |
| Resume parameter mismatch | Exit 1 |
| CDS length not multiple of 3 | Backtrans skipped for that gene, warning in result.json |
| Internal stop codon in CDS | Backtrans skipped for that gene, warning in result.json |
| trimAl exits non-zero | Backtrans skipped for that gene, stderr captured as warning |
| Generated MSA is empty, unparsable, or has unequal sequence lengths | Gene skipped with reason in result.json |
| All genes fail | Exit 1 |
| `--extra-args` used with MAFFT method | Warning printed, args ignored |

## Notes

- Downstream: pass `--msa-dir` of `phyloai pretree trim` to `seqs/` (Mode 1/2) or `seqs/faa/` (Mode 3 AA) as appropriate.
- `result.json` `key_results` feeds the Methods paragraph: "X genes were aligned using MAFFT L-INS-i; mean alignment length Y aa."
- Run `phyloai doctor` to verify MAFFT, MAGUS, and trimAl are detected.
