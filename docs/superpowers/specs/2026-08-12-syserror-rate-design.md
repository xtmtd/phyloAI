# Systematic-Error Site-Rate Analysis Design

**Date:** 2026-08-12
**Status:** Implemented (2026-08-13)

## Purpose

`phyloai posttree syserror rate` is a pure-Python atomic diagnostic for site-rate heterogeneity. It normalizes and ranks site rates from IQ-TREE3 or PhyloBayes, then optionally retains the slowest or fastest site fractions from an input MSA. It supports sensitivity analyses such as retaining 25%, 50%, 75%, and 100% of slow-evolving sites without adding dependencies or invoking an external tool.

This command does not infer rates, choose a scientifically preferred fraction, infer trees, or interpret resulting topology changes. Those decisions remain iterative user/Skill workflow steps.

## Command Interface

```bash
phyloai posttree syserror rate \
  (--iqtree-rate matrix.rate | --pb-rate chain1.meansiterates) \
  [--matrix raw.fa] \
  [--subset slow|fast] \
  [--fraction 0.25,0.5,0.75] \
  [-o runs/posttree/syserror/rate] [--overwrite] [--dry-run] [-q]
```

| Parameter | Type | Default | Rules |
|---|---|---|---|
| `--iqtree-rate` | file | -- | Mutually exclusive and jointly required with `--pb-rate`; reads IQ-TREE3 `--rate` output. |
| `--pb-rate` | file | -- | Mutually exclusive and jointly required with `--iqtree-rate`; reads PhyloBayes `readpb -r` `.meansiterates` output. |
| `--matrix` | file | -- | Optional original MSA. Its presence enables site extraction. Existing supported MSA input formats are accepted. |
| `--subset` | `slow\|fast` | `slow` when extracting | Selection direction. Valid only with `--matrix`; omission resolves to `slow` after a matrix is supplied. |
| `--fraction` | comma-separated floats | -- | One or more fractions in `(0, 1]`; required with `--matrix`, forbidden without it. |
| `--output-dir`, `-o` | directory | `runs/posttree/syserror/rate` | Standard output conflict policy. |
| `--overwrite` | flag | false | Delete and recreate a non-empty output directory. |
| `--dry-run` | flag | false | Parse and validate all inputs, but write no files. |
| `--quiet`, `-q` | flag | false | Suppress terminal output except errors. |

The comma-list form makes single and multi-fraction runs use one interface: `--fraction 0.25` and `--fraction 0.25,0.5,0.75`. Internally the CLI represents an omitted `--subset` as unset, then resolves it to `slow` only for an extraction run; this preserves rejection of explicitly supplied extraction options without `--matrix`.

## Rate Inputs and Normalization

The command parses only site identifier and rate, retaining no input-specific metadata.

| Source | Input layout | Input index | Normalized index |
|---|---|---:|---:|
| IQ-TREE3 | Tabular output headed `Site` and `Rate`; comment lines begin with `#` | 1-based | unchanged |
| PhyloBayes | Whitespace-delimited, headerless `<site> <rate>` rows | 0-based | incremented by one |

Both parsers reject an empty input, malformed rows, duplicate site identifiers, non-integer site identifiers, non-finite rates, negative rates, or indices that are not consecutive from their required origin. Normalized sites must consequently be exactly `1..N`.

The canonical sequence is sorted by `(rate ascending, site ascending)`. The secondary key makes tie handling deterministic. `rates.csv` always has precisely these columns and row order:

```csv
site,rate
21,0.19145
...
```

## Optional MSA Extraction

When `--matrix` is supplied, its parsed alignment must contain records, have equal sequence lengths, and have a length exactly equal to the number of normalized rates.

For every requested fraction `f`, the selected site count is `ceil(N * f)`, where `N` is the alignment length. A non-empty alignment therefore always retains at least one site. Tied rates at the selection boundary are not expanded: each output has a predictable requested length.

- `slow` selects the first `ceil(N * f)` canonical sorted rows.
- `fast` selects the last `ceil(N * f)` canonical sorted rows.
- Selected sites are re-sorted numerically before serialization, preserving their original alignment order.

The command accepts existing PhyloAI-supported FASTA, relaxed PHYLIP, PAML PHYLIP, and NEXUS inputs through `core.formats.FormatConverter`. All generated submatrices are FASTA with lines wrapped at 60 characters.

## Output Layout

Without `--matrix`, the output directory contains only `rates.csv` and `result.json`.

With `--matrix`, the output directory additionally has one directory for every requested fraction:

```text
runs/posttree/syserror/rate/
├── rates.csv
├── slow25/
│   ├── positions.txt
│   └── matrix.fa
├── slow50/
│   ├── positions.txt
│   └── matrix.fa
└── result.json
```

Directory labels use `slow` or `fast` plus the fraction as a percentage with trailing zeroes removed (`0.25 -> slow25`, `0.125 -> slow12.5`, `1 -> slow100`). Fractions which would produce the same label after normalization are rejected. `positions.txt` has one 1-based original site position per line, in ascending alignment order.

## Result and Reporting

The command writes one standard `result.json` with no external tool versions. `params` stores both rate source paths, matrix/extraction parameters, and resolved output settings. `key_results` contains source type, total sites, minimum and maximum rate, and per-subset requested fraction, selected site count, actual fraction, and output directory.

`data.output_files` records `rates.csv`, every `positions.txt`, and every generated `matrix.fa`, each with a description. `data.warnings` is available for non-fatal notices. No checkpoint or external-tool log is created.

The report method template states the rate source (IQ-TREE empirical-Bayes estimate or PhyloBayes posterior mean), selection direction, requested fractions, and actual retained site counts. `rates.csv` participates in the report table index. Positions files and matrix FASTA outputs remain recorded as persistent files in the step output list but are not tables.

## Errors and Lifecycle

The default output-directory conflict policy applies. `--overwrite` and `--dry-run` follow existing utility-command conventions; no resume support is needed.

Validation rejects:

- zero or two rate input options;
- malformed, empty, duplicate, non-consecutive, non-finite, or negative rate input;
- `--subset` or `--fraction` without `--matrix`, and a matrix without `--fraction`;
- invalid, repeated, or label-colliding fractions;
- unparsable, empty, unaligned, or length-mismatched MSA input;
- a non-empty existing output directory without `--overwrite`.

Click accepts path strings without existence prevalidation so `run_rate()` can validate missing or non-file inputs and the CLI can write a standard error `result.json` when the output directory can safely be claimed. CLI validation errors use exit code 1. `--dry-run` writes nothing and displays the validated payload.

When validation fails, PhyloAI never deletes an existing output directory merely to write an error record. If the directory was already non-empty and `--overwrite` was not supplied, it is left entirely untouched. If `--overwrite` was supplied, the CLI may replace or create root `result.json` to record the validation error, but it preserves all other existing files because analysis execution and overwrite cleanup never began.

## Integration

- Add `phyloai/posttree/syserror_rate.py` and register `rate` under the existing `posttree syserror` Click group.
- MCP discovers the command dynamically from the Click tree; no hand-written MCP tool is added.
- Update the English and Chinese command documentation, README command index, `skills/phyloai-workflow/` guidance, the report template, and the high-level design's systematic-error examples and atomic-operation sequence.
- Include detailed Click help with source-specific examples, sorting-only mode, slow/fast single-fraction examples, and multi-fraction sensitivity analysis.

## Verification

Tests cover IQ-TREE and PhyloBayes parsing, normalized indices, deterministic tie sorting, sorting-only mode, slow/fast single and multi-fraction selection, ceiling rounding, ordered positions, FASTA serialization, supported input formats, missing-path and length/parameter validation, output conflicts, dry runs, standard `result.json`, CLI help, report methods text, report table indexing, and generated command/MCP discovery.

No generic site-filtering abstraction, dependency, rate inference, tree inference, or automated choice of biological thresholds is added. A future command may extract sites from other per-site scores only after a second concrete use case establishes a shared input contract.
