# phyloai pretree filter metrics

## Purpose

`phyloai pretree filter metrics` filters whole loci by explicit numeric or string conditions on a metrics CSV/TSV table (typically the output of `phyloai pretree metrics`).

All conditions in `--keep` are combined with AND logic — a locus must satisfy every condition to be retained. When `--msa-dir` is provided, retained-MSA statistics are computed and displayed. Use `--copy` to copy retained MSA and/or tree files into the output directory.

## Usage

Minimal:
```bash
phyloai pretree filter metrics --table ./metrics.csv --keep "dvmc<=0.3,average_BS>=0.8"
```

With file copy:
```bash
phyloai pretree filter metrics \
  --table ./metrics.csv \
  --keep "num_sites>=300,DataType==AA" \
  --copy --msa-dir ./trimmed
```

## Parameters

| Parameter | Default | Notes |
|-----------|---------|-------|
| `--table` | required | Path to metrics CSV/TSV |
| `--keep` | required | Comma-separated AND conditions (see Rule Syntax below) |
| `--output-dir` / `-o` | `runs/pretree/filter/metrics` | Output directory |
| `--input-format` | `auto` | `csv`, `tsv`, or `auto`; auto-detects delimiter from file content |
| `--loci-column` | `loci` | Column name holding the logical locus identifier |
| `--msa-dir` | — | Directory of MSA files for computing retained-MSA statistics |
| `--tree-dir` | — | Directory of tree files for copy mode |
| `--copy` | off | Copy retained MSAs/trees to output directory; requires `--msa-dir` or `--tree-dir` |
| `--table-format` | `csv` | `csv` or `tsv` for auxiliary tables |
| `--overwrite` | off | Delete and recreate non-empty output directory |
| `--dry-run` | off | Parse rules and report counts without writing files |
| `--quiet` / `-q` | off | Suppress terminal output except errors |

## Rule Syntax

Supported operators: `>=`, `>`, `<=`, `<`, `==`, `!=`.

Numeric comparisons (`>=`, `>`, `<=`, `<`) require numeric values. String comparisons (`==`, `!=`) work on any column. Using `>=`/`>`/`<=`/`<` with a string value exits with an error.

```bash
# Numeric conditions (AND only)
--keep "dvmc>=0,dvmc<=0.3,average_BS>=0.8"

# Mixed numeric + string
--keep "DataType==AA,num_sites>=300"

# Simple threshold
--keep "num_sites>=1000"
```

## Inputs

The `--table` file must be a CSV or TSV with a header row. The locus identifier column (default `loci`) identifies each row. Empty files cause an error. All other columns become available for filtering.

## Outputs

```
runs/pretree/filter/metrics/
├── retained_loci.csv|tsv
├── dropped_loci.csv|tsv
├── filter_decisions.csv|tsv
├── seqs/                         (only with --copy --msa-dir)
├── trees/                        (only with --copy --tree-dir)
├── filter.log
└── result.json
```

Terminal summary includes total/retained/dropped counts, condition failure counts, and retained MSA statistics when `--msa-dir` is provided (MSA count, total/mean/min/max alignment length, mean taxa count).

## Examples

```bash
# Basic numeric filter
phyloai pretree filter metrics --table ./metrics.csv --keep "dvmc<=0.3,average_BS>=0.8"

# String + numeric filter with copy
phyloai pretree filter metrics \
  --table ./metrics.csv \
  --keep "DataType==AA,num_sites>=300" \
  --copy --msa-dir ./trimmed

# Dry-run to explore thresholds
phyloai pretree filter metrics --table ./metrics.csv --keep "num_sites>=500" --dry-run

# Custom loci column
phyloai pretree filter metrics --table ./table.tsv --keep "average_BS>=0.9" --loci-column gene_id
```

## Warnings and Errors

| Condition | Behaviour |
|-----------|-----------|
| `--table` does not exist | Exit 1 |
| `--table` is empty | Exit 1 |
| `--keep` is malformed | Exit 1 with syntax detail |
| Unknown column referenced in `--keep` | Exit 1 |
| Non-equal operator on string column | Exit 1 |
| `--copy` without `--msa-dir` or `--tree-dir` | Exit 1 |
| No loci match conditions | Result reports 0 retained |
| Non-empty output directory without `--overwrite` | Exit 1 |

## Notes

- Only AND logic is supported. OR logic is not available in this version.
- Without `--copy` only decision tables are written, which is useful for threshold exploration without duplicating large files.
- Delimiter auto-detection inspects the first 1024 bytes of the input file; use `--input-format csv|tsv` to override when auto-detection is ambiguous.
- `result.json` includes `condition_failure_counts` mapping each condition to how many loci it rejected, useful for understanding which condition is most restrictive.
