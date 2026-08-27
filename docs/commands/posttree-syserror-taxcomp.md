# phyloai posttree syserror taxcomp

[English](posttree-syserror-taxcomp.md) | [中文](posttree-syserror-taxcomp.zh.md)

## Purpose

Screen compositional heterogeneity **across taxa** in one nucleotide or
amino-acid alignment. `taxcomp` builds a taxon-by-state count table and
reports two complementary views from the same table:

1. A **Pearson common-composition chi-square test**: an overall homogeneity
   statistic (P4-style omnibus) plus each taxon's row contribution, with
   nominal per-taxon p-values (the IQ-TREE/P4-style per-taxon screening
   values) and per-taxon Holm-adjusted exploratory p-values.
2. The **observed PPA-COMP descriptive statistics** already used by
   `phyloai posttree simulate adequacy`: squared composition distance per
   taxon, `comp_max`, and `comp_mean`.

`taxcomp` is a pure-Python diagnostic. It uses no tree, no substitution
model, and no external executable, so `phyloai doctor`, checkpointing, and
resume support are not needed.

**Non-goals.** This command does not delete taxa, recode data, select a
model or topology, recommend a threshold, or run posterior-predictive
simulation. Its p-values are nominal exploratory results, not
phylogenetically calibrated model-adequacy tests. A low p-value or a large
composition distance is a prompt for inspection, not a verdict.

## Usage

```bash
phyloai posttree syserror taxcomp \
  --matrix matrix.aa.fa \
  [--seq-type AA|NT|auto] \
  [--table-format csv|tsv] \
  [-o runs/posttree/syserror/taxcomp] \
  [--overwrite] [--dry-run] [-q]
```

## Examples

Screen an amino-acid alignment with explicit sequence type:

```bash
phyloai posttree syserror taxcomp --matrix matrix.aa.fa --seq-type AA
```

Screen a nucleotide alignment to TSV summaries in a custom directory:

```bash
phyloai posttree syserror taxcomp --matrix matrix.nt.fa --seq-type NT --table-format tsv -o runs/posttree/syserror/taxcomp-nt
```

## Inputs

| Option | Required | Default | Description |
|---|---|---|---|
| `--matrix` | yes | -- | One aligned FASTA, PHYLIP, PHYLIP-PAML, or Nexus MSA. Format is auto-detected with the existing PhyloAI readers; Clustal is not supported. No tree or model input is needed. |
| `--seq-type` | no | `auto` | `AA`, `NT`, or automatic detection. AA counts only standard amino acids and NT only `ACGT`; every other character (gaps, unknowns, ambiguity codes, stops) is treated as missing and contributes no fractional count. |
| `--table-format` | no | `csv` | Delimiter and suffix for both summary tables (`csv` or `tsv`). |
| `-o`, `--output-dir` | no | `runs/posttree/syserror/taxcomp` | Output directory. Standard non-empty-directory conflict policy applies. |
| `--overwrite` | no | false | Delete and recreate a non-empty output directory only after validation succeeds. |
| `--dry-run` | no | false | Validate, parse, and calculate all summaries without writing any file. |
| `-q`, `--quiet` | no | false | Suppress terminal output except errors. |

The alignment must contain at least two uniquely named taxa and at least two
globally observed standard states; all sequences must be aligned to the same
length. Duplicate taxon identifiers, an empty or unreadable alignment, a
taxon with zero valid characters, and unequal sequence lengths are hard
errors.

## Outputs

```
<output_dir>/
├── overall_summary.csv|tsv
├── taxon_summary.csv|tsv
└── result.json
```

### overall_summary

One row with stable columns:

| Column | Meaning |
|---|---|
| `n_taxa` | Number of taxa |
| `n_states` | Effective Pearson `K`: number of globally observed standard states entering the expected counts and degrees of freedom |
| `x2` | Overall Pearson chi-square (sum of per-taxon contributions) |
| `df` | `(n_taxa - 1) * (n_states - 1)` |
| `p_nominal` | Nominal chi-square survival probability for the overall test |
| `sparse_count_check` | `triggered` or `not_triggered` |
| `expected_cells_total` | Number of expected cells (`n_taxa * n_states`) |
| `expected_cells_below_1` | Expected cells `< 1` |
| `expected_cells_below_5` | Expected cells `< 5` |
| `expected_cells_below_5_fraction` | Fraction of expected cells `< 5` |
| `comp_max` | Maximum per-taxon squared composition distance |
| `comp_mean` | Mean per-taxon squared composition distance |

### taxon_summary

One row per taxon, in input order:

| Column | Meaning |
|---|---|
| `taxon` | Taxon identifier, preserved exactly |
| `x2_contribution` | That taxon's row contribution to the overall X2 |
| `df` | `n_states - 1` (per-taxon screening degrees of freedom) |
| `p_nominal` | Nominal per-taxon chi-square p-value (IQ-TREE/P4-style screening value) |
| `p_holm` | Holm step-down adjusted nominal p-value, restored to input order |
| `squared_composition_distance` | Unitless squared Euclidean frequency discrepancy from the equal-taxon mean composition; this is the per-taxon observed PPA-COMP value (the `obs comp` reported by PhyloBayes `chain1.comp` and the `obs` column of `posttree simulate adequacy`'s `adequacy_taxon_comp.csv`) |

## Interpretation

- **Overall X2** pools evidence against one common taxon composition under
  the conventional chi-square screen. It is a PhyloAI/P4-style extension; it
  is not an IQ-TREE output.
- **Taxon X2 contribution** is each row's contribution to the overall X2. Its
  p-value is a screening value, **not** an independent one-taxon-versus-rest
  contingency test.
- **`p_holm`** is Holm's multiplicity adjustment applied to the nominal
  per-taxon p-values. When the underlying marginal p-values are valid, Holm
  controls family-wise error under arbitrary dependence, but here it inherits
  the nominal chi-square calibration limitations.
- **`sparse_count_check`** reports whether the conventional sparse-cell rule
  fired (any expected cell `< 1`, or more than 20% of expected cells `< 5`).
  `not_triggered` does **not** establish independence, validate the
  phylogenetic null model, or turn the screening p-value into a
  posterior-predictive test; `triggered` warns that the nominal p-values are
  especially unreliable.
- **Squared composition distance** is a unitless squared Euclidean
  frequency-space discrepancy from the equal-taxon mean. It is **not** an
  evolutionary distance and has **no universal cutoff**; `comp_max` and
  `comp_mean` describe the largest and average observed departures. Values
  across different alphabets, taxon sets, and missing-data patterns need not
  be directly comparable.

Phylogenetic dependence limits conventional chi-square p-value
interpretation even when the sparse-cell rule is not triggered. The report
template repeats these caveats; it never classifies a taxon as significant,
failed, outlier, biased, or removable.

## Sensitivity follow-up

A large distance or low p-value is a prompt for inspecting annotation,
coverage, contamination, lineage composition, and model adequacy. Recoding is
a separate, user-approved sensitivity analysis using the existing concat
workflow and must be reported alongside the original analysis:

- AA: `phyloai pretree concat --recoding Dayhoff-6`
- NT: `phyloai pretree concat --recoding RY-nucleotide`

Taxon removal is a manual data-curation decision requiring independent
evidence and, where applicable, comparison with model-calibrated
`phyloai posttree simulate adequacy` results. This command never produces a
removal list.

## Warnings / Errors

- A non-empty output directory is refused unless `--overwrite` is given.
- Validation failures before the directory is claimed write no files.
- `--overwrite` deletes and recreates the directory only after validation
  succeeds; with an error, it preserves existing files and may replace only
  the root `result.json`.
- `--dry-run` writes nothing and prints the computed summaries as JSON.
- `result.json` follows the standard PhyloAI schema with
  `error_category: null` on success and a standard category on error.

## Notes

- `taxcomp` is a screen only: it never deletes taxa, recodes data, selects a
  model or topology, or recommends a threshold.
- For theory, evidence limits, recoding/taxon-sampling sensitivity choices, and
  the boundary with nonstationary composition models, see the
  [systematic-error workflow reference](../../skills/phyloai-workflow/references/syserror-workflow.md).
- All p-values are nominal and exploratory; phylogenetic dependence limits
  conventional chi-square interpretation even when the sparse-cell rule is
  not triggered.
- `squared_composition_distance` is the per-taxon observed PPA-COMP value
  (`obs comp`); model-calibrated interpretation requires
  `posttree simulate adequacy`, which adds simulation-based z-scores and
  posterior-predictive pp values.
- Values across different alphabets, taxon sets, and missing-data patterns
  are not directly comparable.
