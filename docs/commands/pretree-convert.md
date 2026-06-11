# phyloai pretree convert

## Purpose

`phyloai pretree convert` normalizes sequence characters and converts supported sequence/alignment formats for the PhyloAI workflow. It is intended for FASTA, Phylip-relaxed, Phylip-PAML, and Nexus files only; it is not a general-purpose format conversion tool.

## Usage

```bash
phyloai pretree convert --input ./raw --output-dir ./runs/run001/pretree/convert --to fasta
```

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--input` | required | Input directory or single file |
| `--output-dir`, `-o` | `runs/run001/pretree/convert` | Directory for converted files |
| `--to` | `fasta` | Target format: `fasta`, `phylip-relaxed`, `phylip-paml`, `nexus` |
| `--input-format` | `auto` | Override format detection |
| `--seq-type` | `auto` | Override molecule type detection |
| `--aa-special` | `x` | Convert `B/Z/J/X/U/O` to `X`, or preserve with `keep` |
| `--threads`, `-t` | `4` | Directory-mode worker count |
| `--output-format` | `json` | Structured output format |
| `--quiet`, `-q` | false | Suppress Rich terminal output except errors |
| `--overwrite` | false | Delete and recreate a non-empty output directory |

## Inputs

`--input` may be a directory or a single file. Directory mode scans one level only and skips subdirectories, empty files, non-sequence files, and files that cannot be parsed.

## Outputs

Converted files are written to `--output-dir`. Target suffixes are `.fa`, `.phy`, `.paml.phy`, and `.nex`.

The JSON payload contains `summary`, `files`, `skipped`, and `warnings` under `data`. `key_results` is empty because `convert` is a utility command.

## Examples

```bash
phyloai pretree convert --input ./raw
phyloai pretree stats --seq-dir ./runs/run001/pretree/convert
phyloai pretree convert --input ./gene.phy --output-dir ./converted --to fasta --seq-type NT
phyloai pretree convert --input ./aligned --to phylip-paml --overwrite
```

## Warnings and Errors

If some files are invalid, they are skipped and listed in the output. If all inputs fail or are skipped, the command exits with code 1. If the output directory exists and is non-empty, use `--overwrite` to replace it.

## Notes

Use `pretree convert` before `pretree stats` when raw input files may contain mixed formats or non-standard characters.
