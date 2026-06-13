# phyloai pretree convert

## Purpose

`phyloai pretree convert` normalizes sequence characters and converts supported sequence/alignment formats for the PhyloAI workflow. It is intended for FASTA, Phylip-relaxed, Phylip-PAML, and Nexus files only; it is not a general-purpose format conversion tool.

## Usage

```bash
phyloai pretree convert --input ./raw --output-dir ./runs/pretree/convert --to fasta
```

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--input` | required | Input directory or single file |
| `--output-dir` | `runs/pretree/convert` | Directory for converted files and result.json |
| `--to` | `fasta` | Target format: `fasta`, `phylip-relaxed`, `phylip-paml`, `nexus` |
| `--input-format` | `auto` | Override format detection |
| `--seq-type` | `auto` | Override molecule type detection |
| `--aa-special` | `x` | Convert `B/Z/J/X/U/O` to `X`, or preserve with `keep` |
| `--threads` | `4` | Directory-mode worker count |
| `--quiet` | false | Suppress Rich terminal output except errors |
| `--overwrite` | false | Delete and recreate a non-empty output directory |

## Terminal Output

By default, the command displays a Rich progress bar and summary table in the terminal. The full JSON result is written to `result.json` inside `--output-dir`. Use `--quiet` to suppress terminal output.

## Inputs

`--input` may be a directory or a single file. Directory mode scans one level only and skips subdirectories, empty files, non-sequence files, and files that cannot be parsed.

## Outputs

Converted sequence files are written to `seqs/` subdirectory inside `--output-dir`. Target suffixes are `.fa`, `.phy`, `.paml.phy`, and `.nex`.

The JSON result is written to `result.json` inside `--output-dir`. The payload contains `summary`, `files`, `skipped`, and `warnings` under `data`. `key_results` is empty because `convert` is a utility command.

For `phylip-paml`, output records use the PAML sequential header form with `S`, write the normalized taxon name in a 30-character field, and place at least two spaces between the name field and the sequence.

Example output structure:
```
runs/pretree/convert/
├── seqs/
│   ├── gene1.fa
│   ├── gene2.fa
│   └── ...
└── result.json
```

## Examples

```bash
phyloai pretree convert --input ./raw
phyloai pretree stats --seq-dir ./runs/pretree/convert/seqs
phyloai pretree convert --input ./gene.phy --output-dir ./converted --to fasta --seq-type NT
phyloai pretree convert --input ./aligned --to phylip-paml --overwrite
```

## Warnings and Errors

If some files are invalid, they are skipped and listed in the output. If all inputs fail or are skipped, the command exits with code 1. If the output directory exists and is non-empty, use `--overwrite` to replace it.

## Notes

Use `pretree convert` before `pretree stats` when raw input files may contain mixed formats or non-standard characters.

Even if the source files are already FASTA, running `pretree convert --to fasta` is still recommended before downstream steps when the data source is uncertain. It re-checks parsing, normalizes sequence characters, and can catch malformed records, unexpected symbols, or other input problems early.
