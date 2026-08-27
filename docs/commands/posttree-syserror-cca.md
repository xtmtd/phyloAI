# phyloai posttree syserror cca

[English](posttree-syserror-cca.md) | [中文](posttree-syserror-cca.zh.md)

## Purpose

Perform Compositional Constraint Analysis (CCA), a local
composition-constraint diagnostic that compares site-wise topology preference
across two model analyses. For each site, CCA calculates the effective number
of amino acids (Keff) from 20 site-specific amino-acid frequencies and the
Tree2-minus-Tree1 log-likelihood difference independently for each model.

CCA is diagnostic only. It does not generate site frequencies or likelihoods,
infer trees, filter sites, select a model, or establish the true topology.
Prepare the inputs separately, normally using IQ-TREE PMSF or
`phyloai tree bi readpb --mode ss` for site frequencies and two
`phyloai posttree signal lnl` runs for likelihood tables.

This is a pure-Python command with no external executable; `phyloai doctor`,
checkpointing, and resume support are not needed.

## Usage

```bash
phyloai posttree syserror cca \
  --site-freq chain1.sitefreq \
  --site-lnl1 lnl_LG/site_lnl.csv \
  --site-lnl2 lnl_C20/site_lnl.csv \
  [--model1-name model1] [--model2-name model2] \
  [--title TEXT] \
  [--xlabel "Effective number of amino acids"] \
  [--ylabel "Log-likelihood difference"] \
  [--fig-width 10] [--fig-height 6] [--dpi 300] [--font-size 16] \
  [-o runs/posttree/syserror/cca] [--overwrite] [--dry-run] [-q]
```

## Inputs

| Option | Required | Default | Description |
|---|---|---|---|
| `--site-freq` | yes | -- | IQ-TREE PMSF or PhyloAI-converted `readpb --mode ss` `.sitefreq` table. Raw PhyloBayes `.siteprofiles` is not accepted. |
| `--site-lnl1` | yes | -- | `site_lnl.csv` for model analysis 1. Must contain `site`, `lnL_Tree1`, and `lnL_Tree2`. |
| `--site-lnl2` | yes | -- | `site_lnl.csv` for model analysis 2. Must contain `site`, `lnL_Tree1`, and `lnL_Tree2`. |
| `--model1-name` | no | `model1` | Non-empty, distinct label for model analysis 1 in the CSV, legend, JSON, and report. |
| `--model2-name` | no | `model2` | Non-empty, distinct label for model analysis 2 in the CSV, legend, JSON, and report. |
| `--title` | no | empty | Optional PDF title. |
| `--xlabel` | no | `Effective number of amino acids` | X-axis label. |
| `--ylabel` | no | `Log-likelihood difference` | Y-axis label. |
| `--fig-width` | no | `10` | Positive PDF width in inches. |
| `--fig-height` | no | `6` | Positive PDF height in inches. |
| `--dpi` | no | `300` | Positive PDF rasterization DPI metadata. |
| `--font-size` | no | `16` | Positive figure-legend text size in points. |
| `-o`, `--output-dir` | no | `runs/posttree/syserror/cca` | Output directory. |
| `--overwrite` | no | false | Delete and recreate a non-empty output directory. |
| `--dry-run` | no | false | Validate and calculate the result payload without writing files. |
| `-q`, `--quiet` | no | false | Suppress terminal output except errors. |

Each `.sitefreq` data row contains a one-based site identifier followed by
exactly 20 finite, non-negative frequencies summing to one within `1e-6`.
Site identifiers must be unique and exactly consecutive `1..N`. Both LNL CSVs
require the literal headers `site`, `lnL_Tree1`, and `lnL_Tree2`; extra columns,
including `ΔSLS` and `support`, are ignored. Their rows may be unordered, but
the complete one-based consecutive site set in every input must match exactly.

## Calculation And CSV

For each site:

```text
Keff = 1 / sum(p_i^2), for i = 1..20
delta_lnl_tree2_tree1 = lnl_tree2 - lnl_tree1
```

CCA deliberately recalculates the second expression; it never reuses
`site_lnl.csv:ΔSLS`, whose signal-LNL convention is Tree1 minus Tree2.
Positive CCA values support Tree2; negative values support Tree1.

`cca.csv` has exactly these ASCII snake-case columns, sorted by site with
model 1 then model 2 at each site:

```csv
model,site,keff,lnl_tree1,lnl_tree2,delta_lnl_tree2_tree1
LG,1,11.974845235298696,-14.2296,-14.3580,-0.1284
C20,1,11.974845235298696,-13.8521,-13.9077,-0.0556
```

The CSV writes `lnl_tree1`, `lnl_tree2`, and `delta_lnl_tree2_tree1` to four decimal places; `keff` retains its full floating-point precision. The example below is the historical training `cca.txt` reference. Current bundled `site_lnl.csv` inputs instead yield site-1 deltas of `-0.0999` (LG) and `-0.0436` (C20).

These names map directly to the training `cca.txt` fields: `keff` → `Keff`,
`lnl_tree1` → `LnL_T1`, `lnl_tree2` → `LnL_T2`, and
`delta_lnl_tree2_tree1` → `deltaLnL_T2_T1`.

## Figure

CCA bins each model's rows by `floor(Keff)`, sums the CCA contrast per bin, and
writes a publication-oriented PDF. Bins are fixed at 1 through 20; missing
model/bin combinations are zero. Keff 20 belongs to bin 20.

The plot follows the training ggplot semantics: paired bars are centred at
`bin + 0.5` with width/dodge 1. The plotting range extends to 21 while ticks
remain 1–20 so the valid Keff=20 bin is fully visible. As in ggplot's default discrete fill scale,
model labels are sorted alphabetically before assigning `#F8766D` then
`#00BFC4`: for the training `C20`/`LG` labels, C20 is orange-red and LG is
blue-green. Vertical bin boundaries at 1–20 are solid grey 0.1-point lines;
a black zero baseline is drawn. The x scale has ticks 1–20, no expansion, and
no x major grid; the y major grid remains. The y
limits are `min(0, 1.1 * minimum_bin_sum)` and
`max(0, 1.1 * maximum_bin_sum)`, with `#ffdab9` / light-blue 50%-opacity
negative / positive backgrounds. The upper-right legend is anchored at
`(0.99, 0.9)`, has no title, a semi-transparent
white background, a black 0.5-point border, and the configured legend text size; remaining
text follows the approximately 11-point `theme_bw()` baseline.

## Outputs

```text
runs/posttree/syserror/cca/
├── cca.csv
├── cca.pdf
└── result.json
```

No PNG or aggregate duplicate table is created. `result.json` records all
resolved parameters, Keff and per-model delta summaries, the 20-bin plot
summaries, and absolute paths/descriptions for the CSV and PDF.

## Examples

```bash
# Compare LG and C20 likelihood analyses using shared site frequencies
phyloai posttree syserror cca \
  --site-freq chain1.sitefreq \
  --site-lnl1 lnl_LG/site_lnl.csv --site-lnl2 lnl_C20/site_lnl.csv \
  --model1-name LG --model2-name C20

# Validate inputs and calculation without creating an output directory
phyloai posttree syserror cca --site-freq chain1.sitefreq \
  --site-lnl1 lnl1/site_lnl.csv --site-lnl2 lnl2/site_lnl.csv --dry-run
```

## Warnings / Errors

- All input paths must be readable regular files.
- Site frequencies must have exactly 20 valid amino-acid frequencies per site;
  raw `.siteprofiles` and zero-based site IDs fail the format validation.
- Each LNL table must have the exact required headers; custom topology labels
  in header names are not accepted.
- Empty, malformed, duplicate, non-integer, non-finite, or non-consecutive
  site identifiers/values fail, as does any site-set mismatch across inputs.
- Model labels must be non-empty and distinct; figure dimensions, DPI, and
  legend font size must be positive.
- A non-empty output directory requires `--overwrite`. Dry-run writes no files;
  validation happens before overwrite deletion. If validation fails with
  `--overwrite`, existing files are retained but root `result.json` is replaced
  with the input-error record. No resume/checkpoint exists.

## Notes

- Keff measures the effective number of amino acids from inverse homozygosity,
  not entropy; CCA implements no entropy alternative.
- CCA can show whether topology preference changes with compositional
  constraint, but it cannot alone demonstrate that either model or topology is
  biologically correct.
- Each `site_lnl.csv` must come from `posttree signal lnl` on the same ordered
  Tree1/Tree2 candidate pair; CCA compares two model analyses, not one
  likelihood run per tree. For `.sitefreq` preparation, CAT/PMSF sensitivity
  choices, and the manual/future boundary for Keff filtering, see the
  [systematic-error workflow reference](../../skills/phyloai-workflow/references/syserror-workflow.md).
