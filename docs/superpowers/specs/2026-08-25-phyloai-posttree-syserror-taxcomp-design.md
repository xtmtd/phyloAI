# Taxon Composition Heterogeneity (`taxcomp`) Design

**Date:** 2026-08-25  
**Status:** Draft - pending approval  
**Parent spec:** `2026-06-07-phyloai-design.md` Sections 3, 4, 6, 7, 9, and 11  
**JSON standard:** `2026-06-21-phyloai-json-output-standard.md`  
**Related design:** `2026-08-03-phyloai-posttree-simulate-adequacy-design.md`

## 1. Purpose

`phyloai posttree syserror taxcomp` is a local, pure-Python atomic diagnostic
for compositional heterogeneity across taxa in one nucleotide or amino-acid
alignment. It provides two complementary views from the same taxon-by-state
count table:

1. A Pearson chi-square screening test: a P4-compatible overall homogeneity
   test plus per-taxon row contributions and IQ-TREE/P4-style nominal
   per-taxon screening p-values. IQ-TREE itself reports the per-taxon test;
   the overall Pearson result is a PhyloAI/P4-style extension.
2. The observed PPA-COMP descriptive statistics already implemented by
   `phyloai posttree simulate adequacy`: squared composition distance per
   taxon, `comp_max`, and `comp_mean`.

The command does not use a tree or substitution model. Its chi-square
p-values are nominal exploratory results, not phylogenetically calibrated
model-adequacy tests. It does not delete taxa, recode data, select a model or
topology, recommend a threshold, or invoke posterior-predictive simulation.

## 2. Design Principles

1. **One count table.** Pearson and PPA-COMP statistics use the same parsed
   alignment, standard-state alphabet, and missing-character policy.
2. **No duplicated pretree metrics.** The command does not reproduce general
   completeness, gap, ambiguity, state-frequency, RCFV, or nRCFV tables.
3. **Separate statistic from calibration.** Chi-square values remain valid
   descriptive discrepancies when asymptotic p-value assumptions are weak.
   The output therefore reports the sparse-count diagnostic separately.
4. **No automatic biological decision.** Neither raw nor Holm-adjusted
   p-values identify contamination, justify taxon deletion, or establish that
   a topology is correct.
5. **Reuse existing formulas.** PPA-COMP values must call or share the same
   calculation used by `posttree simulate adequacy`; no second definition is
   introduced.
6. **Minimal output.** Two tables and standard `result.json`; no plot, copied
   alignment, frequency matrix, or redundant summary file.

## 3. Command Interface

```bash
phyloai posttree syserror taxcomp \
  --matrix matrix.aa.fa \
  [--seq-type AA|NT|auto] \
  [--table-format csv|tsv] \
  [-o runs/posttree/syserror/taxcomp] \
  [--overwrite] [--dry-run] [-q]
```

| Parameter | Type | Default | Rules |
|---|---|---|---|
| `--matrix` | file | - | Required aligned FASTA, PHYLIP, PHYLIP-PAML, or Nexus MSA. Format is detected using existing PhyloAI alignment readers. Clustal is not currently supported. |
| `--seq-type` | choice | `auto` | `AA`, `NT`, or automatic detection using existing sequence-normalization logic. The resolved type is recorded. |
| `--table-format` | choice | `csv` | Controls both tabular outputs and their suffixes. |
| `--output-dir`, `-o` | directory | `runs/posttree/syserror/taxcomp` | Standard output-conflict policy. |
| `--overwrite` | flag | false | Delete and recreate a non-empty output directory after input validation succeeds. |
| `--dry-run` | flag | false | Parse, validate, and calculate all summaries without writing files. |
| `--quiet`, `-q` | flag | false | Suppress normal terminal output. |

The command invokes no external program. It needs no `doctor`, thread,
checkpoint, or resume option.

## 4. Inputs and Validation

The alignment must contain at least two uniquely named taxa and at least two
globally observed standard states. All sequences must be aligned to the same
length. Duplicate taxon identifiers, an empty alignment, an invalid or
unreadable format, and a taxon with zero valid characters are hard errors.

The valid alphabets and ordering are:

```text
AA: ACDEFGHIKLMNPQRSTVWY
NT: ACGT
```

Input is converted to uppercase. Only standard states are counted. Gaps,
unknowns, stop characters, and ambiguity codes, including `-`, `.`, `?`, `X`,
`N`, nucleotide IUPAC ambiguity codes, `B`, `Z`, and `J`, are treated as
missing and contribute no fractional counts. This matches the existing
`simulate adequacy` PPA-COMP contract and keeps all cells integer-valued.

States absent from the entire alignment are removed before computing Pearson
expected counts and degrees of freedom. This avoids zero expected cells. The
effective Pearson state count `K` is therefore the number of globally observed
standard states, not necessarily 20 for AA or 4 for NT. PPA-COMP continues to
use the existing adequacy helper's resolved full alphabet; its globally zero
coordinates contribute zero, so retaining or omitting them is mathematically
equivalent for the distance.

## 5. Statistics

### 5.1 Shared count table

For taxon `i` and retained state `j`:

```text
O_ij = observed count
R_i  = sum_j O_ij                 taxon row total
C_j  = sum_i O_ij                 state column total
N    = sum_i sum_j O_ij           grand total
E_ij = R_i * C_j / N              expected count under common composition
```

This null model states that all taxa share one composition vector while
allowing taxa to have different valid-character totals.

### 5.2 Overall Pearson chi-square

```text
X2_overall = sum_i sum_j (O_ij - E_ij)^2 / E_ij
df_overall = (T - 1) * (K - 1)
p_nominal  = P(ChiSquare(df_overall) >= X2_overall)
```

`T` is the number of taxa and `K` the number of retained states. No Yates
correction is applied. The nominal survival probability is calculated with
the installed SciPy chi-square distribution. This overall result is not an
IQ-TREE output; it is the conventional P4-style omnibus homogeneity statistic
calculated from the same count table.

This is an exploratory homogeneity test. Homologous sequences are related by
a phylogeny, so the usual independent-sampling chi-square null distribution
does not fully describe these data. A small nominal p-value means the observed
taxon-by-state table is inconsistent with the conventional common-composition
screen under its asymptotic reference distribution; it is not a model-based
adequacy result.

### 5.3 Per-taxon chi-square contribution

For each taxon:

```text
X2_i = sum_j (O_ij - E_ij)^2 / E_ij
df_i = K - 1
p_i  = P(ChiSquare(df_i) >= X2_i)
```

The row contributions satisfy, within floating-point tolerance:

```text
sum_i X2_i = X2_overall
```

The per-taxon p-values are conventional IQ-TREE/P4-style screening values,
not independent one-taxon-versus-rest contingency tests. The command must not
describe them as such.

### 5.4 Holm adjustment

The family of `T` nominal per-taxon p-values is adjusted by Holm's step-down
method. Sort raw values `p_(1) <= ... <= p_(T)` and calculate:

```text
adjusted_(i) = min(1, max_{j<=i} ((T - j + 1) * p_(j)))
```

Adjusted values are returned to original taxon order as `p_holm`. Holm is
applied as a multiplicity adjustment to the nominal per-taxon p-values. When
the underlying marginal p-values are valid, Holm controls family-wise error
under arbitrary dependence. Here `p_holm` inherits the nominal chi-square
calibration limitations and remains exploratory; it does not correct
phylogenetic-dependence or sparse-count limitations. BH/FDR q-values are not
added.

### 5.5 Expected-count diagnostic

The expected-count diagnostic evaluates whether the ordinary asymptotic
chi-square approximation is additionally weakened by a sparse contingency
table. It is not a heterogeneity statistic and does not modify X2.

Across all `T * K` expected cells:

```text
sparse_count_check = triggered
    if any E_ij < 1
    or more than 20% of E_ij are < 5
otherwise sparse_count_check = not_triggered
```

The comparison is strict: exactly 1 is not below 1, exactly 5 is not below 5,
and exactly 20% is not more than 20%. The result also records the cell counts
and fraction so the status is auditable.

`not_triggered` means only that this conventional sparse-cell rule did not
fire. It does not establish independence, validate the phylogenetic null
distribution, or turn the screening p-value into a posterior-predictive test.
`triggered` leaves X2 and p-values in the output but adds a warning that the
nominal p-values are especially unreliable.

### 5.6 PPA-COMP descriptive statistics

For each taxon, convert its standard-state counts to a frequency vector:

```text
f_ij = O_ij / R_i
```

Following the existing `simulate adequacy` implementation, calculate the
unweighted mean of taxon frequency vectors:

```text
mean_f_j = (1 / T) * sum_i f_ij
```

and each taxon's squared composition distance:

```text
squared_composition_distance_i = sum_j (f_ij - mean_f_j)^2
```

Then:

```text
comp_max  = max_i squared_composition_distance_i
comp_mean = mean_i squared_composition_distance_i
```

The per-taxon distance measures the magnitude of that taxon's frequency-vector
departure from the equal-taxon mean composition. It is zero only when the
taxon's observed state frequencies equal the mean. It uses equal state weights
and, because it is frequency-based, does not grow mechanically with alignment
length. `comp_max` describes the largest departure; `comp_mean` describes the
average across-taxon dispersion.

These are descriptive observed-data statistics. Each value is a unitless
squared Euclidean discrepancy: no square root is taken, and it is neither an
evolutionary sequence distance nor a Pearson chi-square contribution. The
statistics have no universal cutoff or standalone p-value. They may be used to
rank taxa or compare otherwise comparable sensitivity analyses, but values
across different alphabets, taxon sets, and missing-data patterns need not be
directly comparable. Model adequacy requires comparing them with a simulated
reference distribution through `posttree simulate adequacy`.

The Pearson expectation uses the pooled, valid-character-weighted composition
`C_j/N`, whereas PPA-COMP uses the equal-taxon mean `mean_f_j`. This deliberate
difference reproduces the respective definitions. It can matter when taxa
have unequal valid-character totals and must be stated in help and reports.

## 6. Output Files

```text
<output_dir>/
├── overall_summary.csv|tsv
├── taxon_summary.csv|tsv
└── result.json
```

### 6.1 Overall summary

`overall_summary` contains exactly one row with stable columns. `n_states` is
the effective Pearson `K`: the number of globally observed states that enters
its expected counts and degrees of freedom. The resolved full AA/NT alphabet
is recorded separately in `data.character_policy`.

```text
n_taxa
n_states
x2
df
p_nominal
sparse_count_check
expected_cells_total
expected_cells_below_1
expected_cells_below_5
expected_cells_below_5_fraction
comp_max
comp_mean
```

### 6.2 Taxon summary

`taxon_summary` contains one row per taxon in input order:

```text
taxon
x2_contribution
df
p_nominal
p_holm
squared_composition_distance
```

The table deliberately omits valid-character count, gap/ambiguity fractions,
state frequencies, residuals, verdict labels, and deletion recommendations.
Those either duplicate pretree metrics, expand the interface without changing
the decision, or invite unsupported automated filtering.

Floating-point values are written with sufficient precision for numerical
round-tripping. Taxon identifiers are preserved exactly. CSV/TSV quoting uses
the standard Python `csv` module.

## 7. Result JSON and Terminal Output

`result.json` follows the pure-Python single-command pattern:

- `tool_versions` is `{}`.
- `params` contains every resolved `run_taxcomp()` parameter.
- `key_results` contains all fields from the one-row overall summary plus the
  resolved `seq_type`, the taxon with the largest X2 contribution, and the
  taxon with the largest squared composition distance. It contains no
  threshold-based taxon counts or pass/fail decisions.
- `error_category` is `null` for success and uses the standard category such
  as `input` or `output` for an error result.
- `data.cmd` is `[]`, `data.tool_stderr` is `""`, and `data.warnings` is always
  present.
- For a non-dry run, `data.output_files` contains `overall_summary` and
  `taxon_summary`, each with an absolute path and description. A dry run
  returns `data.output_files: {}` because no persistent files exist.
- `data.character_policy` records the resolved standard alphabet and that all
  non-standard characters were excluded as missing.

Normal terminal output prints the overall X2, df, nominal p-value, the phrase
`sparse-cell rule triggered` or `sparse-cell rule not triggered`, `comp_max`,
`comp_mean`, and output paths. It does not print a pass/fail verdict or a list
of taxa to remove. `--quiet` suppresses normal output. `--dry-run` returns all
computed summaries to the CLI for JSON printing and writes no files.

## 8. Interpretation and Reporting

Detailed CLI help and the report methods/interpretation text must distinguish:

- overall X2: pooled evidence against one common taxon composition under the
  conventional chi-square screen;
- taxon X2 contribution: each row's contribution to overall X2, with nominal
  and Holm-adjusted nominal exploratory p-values;
- expected-count diagnostic: whether a conventional sparse-cell rule was
  triggered, not a validation of assumptions, effect size, or biological
  result;
- squared composition distance: an unweighted frequency-space departure from
  the mean taxon composition;
- `comp_max` / `comp_mean`: maximum and average observed PPA-COMP discrepancy,
  without model-based calibration.

The report must state that phylogenetic dependence limits conventional
chi-square p-value interpretation even when the sparse-cell rule is not
triggered. It must not use `significant`, `failed`, `outlier`, `biased`, or
`remove` as an automatic taxon classification. A low p-value or large distance
is a prompt for inspecting annotation, coverage, contamination, lineage
composition, recoding sensitivity, and model adequacy. A recoding sensitivity
analysis uses the existing `phyloai pretree concat --recoding SCHEME` workflow
and must be reported alongside the original analysis. The scheme depends on
sequence type: `Dayhoff-6` is AA-only, whereas `RY-nucleotide` is NT-only.
Taxon removal remains a manual data curation decision requiring independent
evidence and, where applicable, comparison with model-calibrated `posttree
simulate adequacy` results.

### 8.1 Report methods acceptance baseline

The generated report methods text must include the input taxon count, alignment
length, and resolved sequence type. It must identify Pearson X2 as an observed
common-composition screening statistic, distinguish each taxon's X2 row
contribution from an independent test, name `p_nominal` and `p_holm` as
exploratory nominal values, and state that Holm addresses multiplicity only.
It must describe the sparse-cell rule as a check of the asymptotic reference
rather than a biological result, and state that phylogenetic dependence remains
unaddressed whether or not the rule is triggered. It must define squared
composition distance as a unitless squared Euclidean frequency discrepancy,
not an evolutionary distance or a p-value; it must describe `comp_max` and
`comp_mean` as the maximum and mean of those observed distances. It must state
that model-calibrated adequacy requires simulated comparison, and that the
command makes no automatic recoding, taxon-removal, model, or topology choice.

For example, given an AA alignment with 6 taxa and 5,604 sites:

> Taxon compositional heterogeneity was screened in a 6-taxon, 5,604-site AA
> alignment using a Pearson common-composition X2 statistic calculated from
> the taxon-by-amino-acid count table. Overall X2 and the contribution of each
> taxon row were reported with nominal chi-square reference p-values; the
> per-taxon nominal p-values were additionally adjusted by Holm's method for
> multiplicity, but both remain exploratory because homologous taxa are
> phylogenetically dependent. A conventional sparse-cell rule was evaluated
> from expected counts (any expected count below 1 or more than 20% below 5);
> its result concerns the asymptotic chi-square reference only and does not
> validate the phylogenetic null model. For each taxon, a unitless squared
> Euclidean composition discrepancy from the equal-taxon mean was calculated;
> its maximum (`comp_max`) and mean (`comp_mean`) summarize observed taxon
> heterogeneity. These observed distances have no universal cutoff and require
> simulated model comparison for adequacy assessment. No taxon was removed,
> data recoded, or topology/model selected automatically.

The rendered text may substitute the actual inputs and values, but it must
preserve the statistical distinctions and non-decision language above.

### 8.2 CLI help acceptance baseline

Detailed `phyloai posttree syserror taxcomp --help` text must cover these points
without requiring users to read the report:

| Flag or section | Required help content |
|---|---|
| Command summary | Screens composition heterogeneity across taxa from one aligned MSA; provides exploratory diagnostics and makes no taxon-removal, recoding, model, or topology decision. |
| `--matrix` | Requires one aligned FASTA, PHYLIP, PHYLIP-PAML, or Nexus MSA; no tree or model input; Clustal is unsupported. |
| `--seq-type` | `auto` resolves AA/NT from the alignment; AA counts only standard amino acids and NT only `ACGT`; gaps and ambiguity codes are excluded as missing. |
| `--table-format` | Chooses CSV or TSV for both summaries. |
| `--output-dir`, `--overwrite`, `--dry-run`, `--quiet` | State normal output lifecycle, destructive overwrite behavior, no-file dry run, and terminal suppression respectively. |
| Interpretation block | Explain that overall X2 is a common-composition screen; taxon p-values and Holm p-values are nominal/exploratory; sparse-cell status is not an assumption pass; squared composition distance, `comp_max`, and `comp_mean` describe observed frequency departure and have no universal cutoff; phylogenetic adequacy requires simulation. |

## 9. Errors and Output Lifecycle

The standard non-empty output-directory conflict policy applies. The command
performs preflight validation, including output-conflict checks, before it
claims the directory. `--overwrite` and `--dry-run` follow the existing
local-utility convention: dry run writes no files; a non-dry run creates or
recreates the output directory only after preflight validation succeeds; later
calculation errors follow the standard error-`result.json` behavior. There is
no resume behavior.

The CLI writes a standard error `result.json` for validation errors when the
target directory can safely be claimed. With `--overwrite`, a validation
failure may replace only root `result.json` while preserving existing other
files; without `--overwrite`, an existing non-empty directory is untouched.
An output path that is a file is never modified.

## 10. Integration

- Add `phyloai/posttree/syserror_taxcomp.py` and register `taxcomp` under
  `phyloai posttree syserror` with detailed Click help and examples.
- Extract only the taxon-composition portion of `simulate_adequacy._compute_statistics()`
  into a private shared helper returning taxon frequencies, per-taxon squared
  distances, `comp_max`, and `comp_mean`; both commands call it. Preserve
  `simulate adequacy`'s public results and do not introduce a public utility
  abstraction or duplicate the three-line distance formula.
- Expose the stable CLI dynamically as MCP tool
  `posttree_syserror_taxcomp`; do not add a hand-written MCP schema.
- Add `posttree.syserror.taxcomp` to report collection order and provide an
  input-aware methods/interpretation template.
- Update the parent design, English and Chinese README/command documentation,
  and `phyloai-workflow` Skill only after this spec and a separate
  implementation plan are approved.

## 11. Verification

Tests will cover:

- exact Pearson X2, overall df, nominal p-value, and row-contribution sum on a
  hand-calculated small alignment;
- per-taxon `df = K - 1` and exact Holm step-down adjustment, including tied
  p-values and restoration of input order;
- removal of globally absent states and the resulting effective `K`;
- expected-count boundaries at 1, 5, and exactly 20%, plus both sparse warning
  triggers;
- identity of `squared_composition_distance`, `comp_max`, and `comp_mean` with
  the existing `simulate adequacy` calculation;
- unequal taxon valid-character totals, verifying pooled Pearson expectations
  versus the equal-taxon PPA-COMP mean;
- exclusion of gaps and ambiguity codes without fractional counts;
- AA, NT, and auto detection; supported alignment formats; input taxon order;
- duplicate taxa, unequal lengths, fewer than two taxa/states, all-missing
  taxa, malformed alignments, and missing input errors;
- exact output columns, table delimiter/suffix, JSON fields, warnings, and
  output-file registration;
- output conflict, overwrite, dry-run, quiet, and error-result behavior;
- detailed CLI help, dynamic MCP discovery, report collection, and report text;
- before freezing the regression fixture, directly compare a known
  no-ambiguity alignment against IQ-TREE's per-sequence composition-test
  output to verify the expected-count convention; the fixture then covers
  per-taxon numerical compatibility, while the P4-style omnibus result is
  separately hand-checked;

No one-taxon-versus-rest test, BH/FDR correction, permutation/bootstrap,
posterior-predictive simulation, residual table, state-frequency table,
plotting, automatic taxon deletion, recoding, topology comparison, GFmix, or
NDCH2 integration is added.
