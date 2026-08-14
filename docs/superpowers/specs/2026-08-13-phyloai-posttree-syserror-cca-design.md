# Compositional Constraint Analysis (CCA) Design

**Date:** 2026-08-13  
**Status:** Implemented (2026-08-13)  
**Parent spec:** `2026-06-07-phyloai-design.md` Sections 3, 4, 6, 7, 9, and 11

## Purpose

`phyloai posttree syserror cca` is a local, pure-Python atomic diagnostic for compositional heterogeneity across alignment sites. It combines site-specific stationary amino-acid frequencies with site-wise likelihoods for two fixed candidate topologies under two substitution-model analyses.

For each site and model, it calculates the effective number of amino acids:

```text
Keff = 1 / sum(p_i^2), for i = 1..20
```

and the topology-specific likelihood contrast:

```text
delta_lnl_tree2_tree1 = lnl_tree2 - lnl_tree1
```

It then sums that contrast within integer `floor(Keff)` bins and renders a composition-constraint plot. Positive values support Tree2 and negative values support Tree1. This follows the training workflow for comparing, for example, a homogeneous model with a more compositionally heterogeneous model.

The command does not generate site frequencies, calculate site likelihoods, infer trees, select a biologically correct topology/model, or filter sites. Users prepare inputs independently, normally with `phyloai tree bi readpb --mode ss` / IQ-TREE PMSF and two `phyloai posttree signal lnl` runs.

## Command Interface

```bash
phyloai posttree syserror cca \
  --site-freq chain1.sitefreq \
  --site-lnl1 lnl_LG/site_lnl.csv \
  --site-lnl2 lnl_CAT-PMSF/site_lnl.csv \
  [--model1-name model1] [--model2-name model2] \
  [--title TEXT] \
  [--xlabel "Effective number of amino acids"] \
  [--ylabel "Log-likelihood difference"] \
  [--fig-width 10] [--fig-height 6] [--dpi 300] [--font-size 16] \
  [-o runs/posttree/syserror/cca] [--overwrite] [--dry-run] [-q]
```

| Parameter | Type | Default | Rules |
|---|---|---|---|
| `--site-freq` | file | — | Required IQ-TREE PMSF or PhyloAI-converted `readpb --mode ss` `.sitefreq` file. |
| `--site-lnl1` | file | — | Required `site_lnl.csv` for the first model analysis. |
| `--site-lnl2` | file | — | Required `site_lnl.csv` for the second model analysis. |
| `--model1-name` | text | `model1` | Model label for `--site-lnl1`; used in `cca.csv`, figure legend, result JSON, and report. Must be non-empty and distinct from model 2. |
| `--model2-name` | text | `model2` | Model label for `--site-lnl2`; used in `cca.csv`, figure legend, result JSON, and report. Must be non-empty and distinct from model 1. |
| `--title` | text | empty | Optional figure title. |
| `--xlabel` | text | `Effective number of amino acids` | Figure x-axis label. |
| `--ylabel` | text | `Log-likelihood difference` | Figure y-axis label. |
| `--fig-width` | positive float | `10` | PDF figure width in inches. |
| `--fig-height` | positive float | `6` | PDF figure height in inches. |
| `--dpi` | positive integer | `300` | Figure rasterization DPI metadata. |
| `--font-size` | positive float | `16` | Figure-legend text size in points. |
| `--output-dir`, `-o` | directory | `runs/posttree/syserror/cca` | Standard output conflict policy. |
| `--overwrite` | flag | false | Delete and recreate a non-empty output directory. |
| `--dry-run` | flag | false | Parse, validate, calculate summaries, and return a payload without writing files. |
| `--quiet`, `-q` | flag | false | Suppress normal terminal output. |

This command invokes no external tool; no `doctor` check, checkpoint, or resume option is required.

## Inputs and Validation

### Site frequencies

`--site-freq` accepts **only** the common `.sitefreq` format emitted by IQ-TREE PMSF or by `phyloai tree bi readpb --mode ss`. It does not accept raw PhyloBayes `.siteprofiles` or any alternative indexing convention.

Each non-empty, non-comment row must contain a 1-based integer site identifier followed by exactly 20 finite, non-negative amino-acid frequencies. Site identifiers must be unique and exactly consecutive from `1` through `N`. Frequency rows must sum to 1 within a documented floating-point tolerance (`1e-6`).

### Site likelihood tables

Each `--site-lnl*` input should normally be the `site_lnl.csv` produced by `phyloai posttree signal lnl`. CCA requires columns named exactly:

```text
site,lnL_Tree1,lnL_Tree2
```

Other columns, including the standard `ΔSLS` and `support` columns, are permitted and ignored. Rows can occur in any order. Each table must have unique, consecutive 1-based sites `1..N`; its complete site set must match the `.sitefreq` site set and the other likelihood table exactly.

CCA never reuses `ΔSLS`: signal LNL defines that column as `lnL_Tree1 - lnL_Tree2`, whereas CCA deliberately uses the opposite sign, `lnL_Tree2 - lnL_Tree1`.

All likelihood values must be finite numbers. Empty files, missing columns, malformed CSV, duplicate/non-integer/non-consecutive site IDs, and any mismatch among the three inputs are hard errors.

## Calculation and Outputs

For every site `s`:

1. Compute `keff[s] = 1 / sum(p[s,i]^2)` from its 20 site frequencies.
2. Read `lnl_tree1[s]` and `lnl_tree2[s]` from each model's LNL table.
3. Compute `delta_lnl_tree2_tree1[s] = lnl_tree2[s] - lnl_tree1[s]` independently for each model.
4. Append one long-format row per model to `cca.csv`, in ascending site order and model-1-then-model-2 order.

`cca.csv` always has exactly these stable, ASCII column names:

```csv
model,site,keff,lnl_tree1,lnl_tree2,delta_lnl_tree2_tree1
LG,1,11.974845235298696,-14.2296,-14.3580,-0.1284
C20,1,11.974845235298696,-13.8521,-13.9077,-0.0556
```

The snake-case schema follows PhyloAI's tabular-output convention. `lnl_tree1`, `lnl_tree2`, and `delta_lnl_tree2_tree1` are serialized to four decimal places; `keff` retains full floating-point precision. Its direct correspondence to the historical training `cca.txt` header is `keff → Keff`, `lnl_tree1 → LnL_T1`, `lnl_tree2 → LnL_T2`, and `delta_lnl_tree2_tree1 → deltaLnL_T2_T1`; only column names differ, not their values or signs. The historical reference's site-1 deltas are `LG=-0.1284` and `C20=-0.0556`; current bundled `site_lnl.csv` inputs instead yield `LG=-0.0999` and `C20=-0.0436`.

`model` uses the corresponding `--model*-name` value, including the defaults `model1` and `model2` when no name is supplied.

The command additionally writes:

```text
runs/posttree/syserror/cca/
├── cca.csv
├── cca.pdf
└── result.json
```

No PNG or redundant aggregate CSV is written.

## Figure Semantics

`cca.pdf` is generated using Matplotlib without a new dependency. It reproduces the training R/ggplot intent:

- Bin each row by `floor(keff)` and sum `delta_lnl_tree2_tree1` by `(model, bin)`.
- Use fixed bins 1 through 20. Missing model/bin combinations are materialized as zero-height bars.
- A mathematical `Keff` of 20 remains in bin 20; normal numerical round-off at the upper boundary is handled safely before validation/binning.
- Plot models as side-by-side bars centred at `bin + 0.5`, matching the training `ggplot` expression. Use `width=1` and `position_dodge(width=1)` semantics. Match ggplot's default discrete fill mapping by sorting model labels alphabetically before assigning `#F8766D` then `#00BFC4`; thus training `C20` is orange-red and `LG` is blue-green.
- Draw vertical bin boundaries at integers 1 through 20, suppress major x grid lines, use a zero baseline, and set ticks to `1..20` without expansion. Extend the plotting range to 21 solely to keep the valid bin-20 bars fully visible.
- Derive plot limits exactly as the R command: `ymin = min(0, min_bin_sum * 1.1)` and `ymax = max(0, max_bin_sum * 1.1)`. Shade `[ymin, 0]` light orange (`#ffdab9`, alpha 0.5) and `[0, ymax]` light blue (alpha 0.5).
- Place a legend with no title at normalized plot coordinates `(0.99, 0.9)`, right/top justified, with white 50%-alpha fill, a black 0.5-point border, and `--font-size` legend text (default 16 pt). The remaining visual baseline follows `theme_bw()` defaults (approximately 11 pt), deliberately preserving the training R command's 11/16 mixed typography rather than treating `--font-size` as a global font setting.
- Use the configured title/labels/size/DPI settings and publication-oriented PDF layout.

The plot diagnoses whether topology preference varies systematically with compositional constraint. It is not by itself evidence that either model or topology is true.

## Result JSON and Reporting

`result.json` follows the pure-Python single-command pattern:

- `tool_versions` is `{}`.
- `params` includes every `run_cca()` parameter in resolved form, including all labels and figure settings.
- `key_results` contains total site count; Keff minimum, maximum, and mean; model labels; total delta log likelihood per model; and the 20-bin summaries used for the plot.
- `data.cmd` is `[]`, `data.tool_stderr` is `""`, and `data.warnings` is present.
- `data.output_files` contains `cca_table` for `cca.csv` and `cca_figure` for `cca.pdf`, with absolute paths and descriptions.

The report method template states that Keff was calculated as the inverse homozygosity of 20 site-specific amino-acid frequencies and that site-wise Tree2-minus-Tree1 likelihood differences were binned by `floor(Keff)` under the two named model analyses. It includes model labels and the number of sites but does not claim a preferred topology. The generic report collector will index `cca.csv` as a table and embed `cca.pdf` as a figure.

## Errors and Output Lifecycle

The default non-empty output-directory conflict policy applies. `--overwrite` and `--dry-run` follow the existing local-utility convention: dry run writes no files; overwrite applies only after validation succeeds. There is no resume support.

The CLI writes a standard error `result.json` for validation errors when the target directory can safely be claimed. With `--overwrite`, it replaces only root `result.json` after validation failure while preserving other existing files; without `--overwrite`, it never alters an existing non-empty directory merely to record an error. An existing output path that is a file remains untouched.

## Integration

- Update the parent design's Section 4.1 CCA example, replacing the obsolete matrix/tree interface (`--matrix --t1 --t2`) with the three prepared-input options (`--site-freq --site-lnl1 --site-lnl2`). Update the Section 11, Phase 9 systematic-error atomic sequence (`brlen → rate → cca → sites`) to describe CCA's prepared inputs and its position after site-frequency and site-likelihood production; update the Section 8 Skill overview if needed for consistency. This is a deliberate public-interface replacement, not an additional CCA mode.
- Add `phyloai/posttree/syserror_cca.py` and register `cca` below `phyloai posttree syserror` with detailed, grouped Click help. The help explicitly requires LNL columns named `site`, `lnL_Tree1`, and `lnL_Tree2`.
- Replace the current `posttree_syserror_cca` MCP stub by the dynamically generated Click-based execution tool; do not add a hand-written MCP schema.
- Update the English and Chinese CCA command documentation, README command/example index, high-level design examples and systematic-error atomic workflow, the PhyloAI workflow Skill, and its parameter annotations.
- Replace the existing placeholder CCA report template with an input-aware CCA methods template. `STEP_ORDER` already contains the correct command identifier.

## Verification

Tests will cover:

- Keff calculation and fixed historical-training anchors: the supplied CCA frequency fixture's site 1 has `keff = 11.974845235298696`; the historical `cca.txt` reference gives `LG=-0.1284` and `C20=-0.0556`, while current bundled likelihood fixtures give `LG=-0.0999` and `C20=-0.0436`; no entropy alternative;
- standard `.sitefreq` parsing, site/frequency validation, and upper-bound binning;
- LNL CSV required-column parsing while ignoring extra columns and arbitrary input row order;
- exact site-set agreement among all three inputs;
- independent re-calculation of the reverse-sign contrast;
- CSV naming, row count/order, requested/default model labels, bin completion, and aggregate values;
- PDF creation plus basic figure configuration;
- all plotting-parameter validation;
- output conflict, overwrite, dry-run, error-result, and standard `result.json` behavior;
- CLI help and dynamic MCP command discovery, including removal of the old stub;
- report methods text and figure/table indexing.

No generic per-site table abstraction, R runtime, automatic upstream IQ-TREE/PhyloBayes execution, entropy Keff calculation, topology/model verdict, or automatic site filtering is added.
