# PhyloAI `pretree concat jackknife` Design

**Date:** 2026-06-30
**Status:** Draft for user review
**Parent spec:** `2026-06-07-phyloai-design.md`, `2026-06-13-phyloai-pretree-concat-design.md`
**JSON standard:** `2026-06-21-phyloai-json-output-standard.md`
**Reference:** `ref/scripts/jackknife_sampling.py`

---

## 1. Purpose

`phyloai pretree concat jackknife` creates gene-jackknife pseudoreplicate matrices from an existing concatenated matrix and its partition file. The main use case is reducing Bayesian inference cost by generating many smaller matrices, for example 100 pseudoreplicates of roughly 50,000 sites each, before running `phyloai tree bi` manually on selected replicates.

This command does not infer trees and does not read the original per-locus MSA directory. It operates on the output of `phyloai pretree concat`, so the same mechanism works for original, recoded, translated, and codon-position-excluded matrix variants.

---

## 2. CLI

```bash
phyloai pretree concat jackknife \
  --matrix runs/pretree/concat/matrix.fa \
  --partitions runs/pretree/concat/matrix.partitions \
  --replicates 100 \
  --target-length 50000 \
  --to fasta \
  --table-format csv \
  -o runs/pretree/concat/jackknife
```

### 2.1 Parameters

| Parameter | Type | Default | Description |
|---|---|---|---|
| `--matrix` | Path, required | - | Existing concatenated matrix. Supports formats readable by `core.formats.FormatConverter`. |
| `--partitions` | Path, required | - | RAxML-style partition file matching `--matrix`. |
| `--replicates` | int >= 1 | `100` | Number of pseudoreplicates to generate. |
| `--target-length` | int >= 1 | `50000` | Minimum total sampled site length per pseudoreplicate. Sampling stops once this length is reached or exceeded. |
| `--prefix` | str | `rep` | Output replicate prefix. Files are named `rep001`, `rep002`, etc. |
| `--to` | `fasta\|phylip-relaxed\|phylip-paml\|nexus` | `fasta` | Output matrix format for pseudoreplicates. |
| `--table-format` | `csv\|tsv` | `csv` | Table format for `jackknife_summary`. |
| `--seed` | int | `42` | Random seed for reproducible sampling. |
| `--output-dir` / `-o` | Path | `<matrix_parent>/jackknife` | Output directory. |
| `--overwrite` | flag | False | Delete and recreate non-empty output directory. |
| `--dry-run` | flag | False | Validate inputs and report planned outputs without writing files. |
| `--quiet` / `-q` | flag | False | Suppress terminal output except errors. |

No `--resume` is supported. The default fixed seed makes runs reproducible; users can pass a different `--seed` to generate a different pseudoreplicate set.

---

## 3. Sampling Semantics

Each pseudoreplicate samples loci from the partition file, not columns directly.

1. Parse all partition records as `(model, locus, start, end)`.
2. Read the concatenated matrix into memory.
3. Validate every partition range is within matrix bounds.
4. For each replicate, shuffle the partition list with the configured random generator.
5. Add loci without replacement until cumulative length is `>= --target-length`.
6. Slice the selected source ranges from every taxon sequence and concatenate them in sampled order.
7. Rewrite partition coordinates from 1 for the replicate matrix.

If one locus is longer than `--target-length`, it may form a one-locus replicate. If total available partition length is less than `--target-length`, the command exits with code 1 because the requested target cannot be reached.

Because sampling slices from an existing concatenated matrix, all taxa present in the input matrix are present in every pseudoreplicate. The command does not separately test taxon completeness per original locus.

---

## 4. Output Layout

```text
runs/pretree/concat/jackknife/
├── rep001/
│   ├── rep001.fa
│   └── rep001.partitions
├── rep002/
│   ├── rep002.fa
│   └── rep002.partitions
├── ...
├── jackknife_summary.csv
└── result.json
```

The matrix extension follows `--to`:

| `--to` | Extension |
|---|---|
| `fasta` | `.fa` |
| `phylip-relaxed` | `.phy` |
| `phylip-paml` | `.phy` |
| `nexus` | `.nex` |

No `loci.txt` is written. The selected locus names and rewritten coordinates are recoverable from each `repXXX.partitions` file.

FASTA outputs are written through `core.formats.FormatConverter.write_alignment()` and follow the project-wide 60-character line wrapping rule.

### 4.1 `jackknife_summary.csv` / `.tsv`

The summary table suffix follows `--table-format`: `jackknife_summary.csv` for `csv`, `jackknife_summary.tsv` for `tsv`.

Columns:

| Column | Description |
|---|---|
| `replicate` | Replicate name, e.g. `rep001`. |
| `matrix` | Replicate matrix path relative to `--output-dir`. |
| `partitions` | Replicate partition path relative to `--output-dir`. |
| `n_loci` | Number of sampled loci. |
| `total_length` | Total pseudoreplicate length after sampling. |
| `target_length` | Requested target length. |
| `seed` | Seed used for the run. |

---

## 5. `result.json`

The command follows the standard PhyloAI result schema. It is a pure-Python utility command, so `tool_versions` is `{}`, `data.cmd` is `[]`, and `data.tool_stderr` is `""`.

```json
{
  "status": "success",
  "command": "phyloai pretree concat jackknife --matrix ... --partitions ... --replicates 100 --target-length 50000 --prefix rep --to fasta --table-format csv --seed 42 --output-dir ...",
  "wall_time": 1.23,
  "tool_versions": {},
  "params": {
    "matrix": "runs/pretree/concat/matrix.fa",
    "partitions": "runs/pretree/concat/matrix.partitions",
    "replicates": 100,
    "target_length": 50000,
    "prefix": "rep",
    "to": "fasta",
    "table_format": "csv",
    "seed": 42,
    "output_dir": "runs/pretree/concat/jackknife",
    "overwrite": false,
    "dry_run": false,
    "quiet": false
  },
  "key_results": {
    "n_replicates": 100,
    "target_length": 50000,
    "min_length": 50012,
    "max_length": 53280,
    "mean_length": 51140.5,
    "min_loci": 18,
    "max_loci": 27
  },
  "error": null,
  "data": {
    "cmd": [],
    "tool_stderr": "",
    "output_files": {
      "summary": {
        "path": "/abs/path/jackknife_summary.csv",
        "description": "Summary table for generated gene-jackknife pseudoreplicates"
      },
      "rep001_matrix": {
        "path": "/abs/path/rep001/rep001.fa",
        "description": "Gene-jackknife pseudoreplicate matrix rep001"
      },
      "rep001_partitions": {
        "path": "/abs/path/rep001/rep001.partitions",
        "description": "Partition file for gene-jackknife pseudoreplicate rep001"
      }
    },
    "replicates": [
      {
        "name": "rep001",
        "matrix": "/abs/path/rep001/rep001.fa",
        "partitions": "/abs/path/rep001/rep001.partitions",
        "n_loci": 23,
        "total_length": 50721,
        "loci": ["geneA", "geneQ", "geneB"]
      }
    ],
    "warnings": []
  }
}
```

`data.replicates[].loci` is kept in `result.json` for machine-readable provenance even though no separate `loci.txt` file is written.

Every persistent matrix and partition file is also listed in `data.output_files` with labels `repXXX_matrix` and `repXXX_partitions`. This keeps report generation and MCP consumers aligned with the JSON Output Standard; `data.replicates[]` provides structured provenance but is not the discovery mechanism for output files.

---

## 6. CLI Integration

`phyloai pretree concat` becomes a Click group with `invoke_without_command=True` so the existing concat behavior remains unchanged while adding a nested command:

```text
phyloai pretree concat [OPTIONS]          # existing full concatenation
phyloai pretree concat jackknife [OPTIONS] # new pseudoreplicate generation
```

The current `concat --msa-dir` behavior remains unchanged. It can read PhyloAI-supported alignment formats via `FormatConverter`; users who need explicit normalization should run `phyloai pretree convert` first.

Implementation keeps the existing `run_concat()` entry point and adds `run_concat_jackknife()` in `phyloai/pretree/concat.py`.

The parent concat spec (`2026-06-13-phyloai-pretree-concat-design.md`) must be updated to record this backward-compatible CLI structure change.

---

## 7. Downstream Integration

### 7.1 Report

Report support must be updated so jackknife runs appear as first-class steps:

- `report.collector.parse_step_id()` maps `phyloai pretree concat jackknife ...` to `pretree.concat.jackknife`.
- `report.collector.STEP_ORDER` includes `pretree.concat.jackknife` immediately after `pretree.concat` and before `tree.*` steps.
- `report.templates` adds `generate_methods_pretree_concat_jackknife()` and registers it in `METHOD_GENERATORS`.
- Report file discovery relies on `data.output_files`, so every replicate matrix and partition file must be listed there.

### 7.2 MCP and Skill

The MCP schema generator is expected to discover the new nested Click command dynamically, but tests must confirm the generated tool name and schema. Expected tool name: `pretree_concat_jackknife`.

Command documentation and workflow guidance must be updated:

- Add jackknife usage to `docs/commands/pretree-concat.md` and Chinese docs if present.
- Update `phyloai-workflow` command guidance if it maintains command examples, parameter cards, or step descriptions outside dynamic MCP schemas.

---

## 8. Validation and Errors

Exit code 1 for:

- Missing or unreadable `--matrix` or `--partitions`.
- Empty partition file.
- Unparseable partition line.
- Partition range outside matrix length.
- Total partition length less than `--target-length`.
- Non-empty output directory without `--overwrite`.
- Invalid `--replicates`, `--target-length`, or `--prefix`.

The command writes an error `result.json` when the output directory can be created safely. If `--dry-run` is set, no files are written.

---

## 9. Acceptance Criteria

- [ ] `phyloai pretree concat` still works with the existing options and behavior.
- [ ] `phyloai pretree concat jackknife --help` lists all jackknife options.
- [ ] Given a small FASTA matrix and matching partitions, `--replicates 2 --target-length N` writes `rep001/rep001.fa`, `rep001/rep001.partitions`, `rep002/...`, `jackknife_summary.csv`, and `result.json`.
- [ ] Replicate matrices retain all input taxa in the same order.
- [ ] Replicate partition coordinates are rewritten from 1 and match the output matrix length.
- [ ] Sampling is reproducible with the default `--seed 42` and with user-provided seed values.
- [ ] `--to phylip-relaxed`, `--to phylip-paml`, and `--to nexus` write the expected output extension and readable matrix format.
- [ ] FASTA outputs are wrapped at 60 characters per line.
- [ ] `--table-format tsv` writes `jackknife_summary.tsv` and records `params.table_format = "tsv"`.
- [ ] `--dry-run` writes nothing and reports the planned replicate count and target length unless `--quiet` is set.
- [ ] Non-empty output directory without `--overwrite` exits code 1.
- [ ] Total partition length below target exits code 1.
- [ ] `result.json` passes JSON Output Standard Section 8 structural assertions, including full command, params completeness, `data.cmd`, `data.tool_stderr`, and `data.output_files` coverage.
- [ ] Report collector maps the command to `pretree.concat.jackknife`, report templates generate methods text, and report output indexes all replicate files from `data.output_files`.
- [ ] MCP schema generation exposes a `pretree_concat_jackknife` tool.

---

## 10. Deliberate Non-Goals

- No automatic `tree bi` execution.
- No separate `loci.txt` files.
- No bootstrap sampling with replacement.
- No per-locus taxon completeness checks from original MSA files.
- No changes to existing `concat --msa-dir` input format support.
