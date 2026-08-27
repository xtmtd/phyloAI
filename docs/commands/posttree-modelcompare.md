# phyloai posttree modelcompare

[English](posttree-modelcompare.md) | [中文](posttree-modelcompare.zh.md)


## Purpose

Performs relative substitution-model comparison and selection through two
independent subcommands:

| Subcommand | Analysis | Core tool |
|------------|----------|-----------|
| `iqtree` | ModelFinder BIC/AIC/AICc model comparison | IQ-TREE3 `-m MF` |
| `pb` | Leave-one-out cross-validation (LOO-CV) / wAIC model comparison | PhyloBayes `.sitelogl` (pure Python) |

`iqtree` compares substitution models on a single alignment using
ModelFinder (Kalyaanamoorthy et al. 2017) and reports BIC/AIC/AICc scores,
weights, and 95% confidence-set membership. `pb` compares fitted models by
their predictive accuracy using LOO-CV and the widely applicable information
criterion (wAIC) computed from site log-likelihood files (Lartillot 2023).

## Usage

```bash
# ModelFinder on a supermatrix with homogeneous models (BIC/AIC/AICc)
phyloai posttree modelcompare iqtree --matrix ./matrix.aa.fa --homogeneous-model LG,WAG --mrate E,G,R

# With heterogeneous mixture models (AA only)
phyloai posttree modelcompare iqtree --matrix ./matrix.aa.fa --homogeneous-model LG --heterogeneous-model C10,C20 --het-mrate G,R

# LOO-CV / wAIC comparison of two fitted models (each directory = one model)
phyloai posttree modelcompare pb --sitelogl-dir ./cat_sitelogl,./gtr_sitelogl

# Explicit chain groups (repeated --sitelogl; each occurrence = one model)
phyloai posttree modelcompare pb --sitelogl model1/c1.sitelogl,model1/c2.sitelogl --sitelogl model2/c1.sitelogl,model2/c2.sitelogl

# Single-model fit reporting (one directory)
phyloai posttree modelcompare pb --sitelogl-dir ./gtr_sitelogl --model-names GTR
```

## modelcompare iqtree — ModelFinder Model Comparison

### Purpose

Runs IQ-TREE3 `-m MF` with a homogeneous model search space (`-mset`) and
rate heterogeneity types (`-mrate`), plus optional heterogeneous mixture
models expanded via `-madd`. Parses the "List of models sorted by BIC scores:"
section of the `.iqtree` report into a comparison table.

### Inputs

| Input | Description |
|-------|-------------|
| `--matrix` | Single supermatrix alignment (FASTA, PHYLIP, NEXUS). Required. Maps to IQ-TREE `-s`. |
| `--homogeneous-model` | Comma-separated standard models for the homogeneous search space. Required. Maps to IQ-TREE `-mset`. |
| `--mrate` | Rate heterogeneity types for homogeneous models. Valid: any subset of `E`, `G`, `R` (comma-separated). Default `E,G`. Maps to IQ-TREE `-mrate`. |
| `--heterogeneous-model` | Comma-separated AA mixture models (`C10`–`C60`, `EX*`, `EHO`, `UL*`, `EX_EHO`, `LG4M`, `LG4X`) evaluated via `-madd`. AA only; rejected for NT data. |
| `--het-mrate` | Rate heterogeneity for heterogeneous expansion. Each token selects a variant family: `E` = base models (`C10, C10+F`), `G` = `+G4`, `R` = `+R4`. Valid: any subset of `E`, `G`, `R` (comma-separated). Default `E,G`. AA only. |
| `--seq-type` | `AA`, `NT`, or `auto` (default `auto`). When `auto`, the alignment is read to detect AA vs NT before model validation. |
| `--prefix` | IQ-TREE output prefix (default: `modelcompare`). Must be a single filename — path separators, `..`, and absolute paths are rejected so output stays inside the run directory. |
| `--threads` | IQ-TREE `-T` value (integer or `auto`, default `auto`). |
| `--iqtree-path` | Explicit path to iqtree3 executable. |
| `--tool-args` | Extra IQ-TREE flags. Blocked: `-s`, `--prefix`. When a flag PhyloAI also manages (`-m`, `-mset`, `-mrate`, `-madd`, `-cmin`, `-cmax`, `-T`) is given, the `--tool-args` value overrides the PhyloAI-generated one (no duplicates). |
| `--output-dir` | Output directory (default: `runs/posttree/modelcompare/iqtree`). |
| `--overwrite` | Delete and recreate output directory before running. |
| `--resume` | Resume incomplete IQ-TREE run (native checkpoint). |
| `--dry-run` | Print the IQ-TREE command without executing. |
| `--quiet` | Suppress terminal output except errors. |

### Heterogeneous model expansion

Each token in `--het-mrate` selects a variant family for each model M in
`--heterogeneous-model`, mirroring `--mrate` semantics:

- `E` → `M, M+F` (empirical state frequencies, no rate category)
- `G` → `M+G4, M+F+G4`
- `R` → `M+R4, M+F+R4`

Only the requested families are produced. For example,
`--heterogeneous-model C10 --het-mrate E,G` produces:

```
C10, C10+F, C10+G4, C10+F+G4
```

while `--het-mrate G` alone produces just `C10+G4, C10+F+G4` (no base `C10`).

The expanded list is joined comma-separated and passed to `-madd`.

### Outputs

```
runs/posttree/modelcompare/iqtree/
├── result.json
├── model_fit.csv                 # Rank,Model,LogL,AIC,w_AIC,In_AIC_95,AICc,w_AICc,In_AICc_95,BIC,w_BIC,In_BIC_95
└── iqtree/
    ├── modelcompare.iqtree       # IQ-TREE native report (BIC/AIC/AICc table)
    ├── modelcompare.log
    ├── modelcompare.model.gz
    └── modelcompare.treefile
```

`model_fit.csv` is sorted by BIC. `In_AIC_95` / `In_AICc_95` / `In_BIC_95`
indicate membership in each criterion's 95% confidence set (`+` in IQ-TREE
output).

### Example

```bash
phyloai posttree modelcompare iqtree --matrix concat.aa.fa --homogeneous-model LG,WAG --mrate E,G,R --heterogeneous-model C10,C20 --het-mrate G,R
```

---

## modelcompare pb — LOO-CV / wAIC Model Comparison

### Purpose

Compares fitted models by predictive accuracy using leave-one-out
cross-validation (LOO-CV) and the widely applicable information criterion
(wAIC, Watanabe 2009), following Lartillot (2023). Computed entirely in Python
from PhyloBayes `.sitelogl` site log-likelihood files (no external tool).

### Inputs

| Input | Description |
|-------|-------------|
| `--sitelogl-dir` | Comma-separated directories; each directory represents one model, globbing `*.sitelogl` (≥2 files per directory). Mutually exclusive with `--sitelogl`. Shell path completion supported. |
| `--sitelogl` | Repeatable option; each occurrence is comma-separated `.sitelogl` file paths for one model (≥2 files per model). Mutually exclusive with `--sitelogl-dir`. Shell path completion supported. |
| `--model-names` | Comma-separated model labels matching the number of model groups. If omitted, models are labeled `model_1`, `model_2`, etc. Each label must be a single path component (no `/`, `..`) since it names the output subdirectory under `sitelogl/`. |
| `--output-dir` | Output directory (default: `runs/posttree/modelcompare/pb`). |
| `--overwrite` | Delete and recreate output directory before running. |
| `--quiet` | Suppress terminal output except errors. |

Validation rules:
- Exactly one of `--sitelogl-dir` or `--sitelogl` must be provided.
- At least 1 model group is required; at least 2 `.sitelogl` files per model group.
- All `.sitelogl` files within a model group must have the same number of data rows.
- **Cross-model site validation:** when ≥2 model groups are provided, all groups
  must have the same site count and identical ordered `site` identifiers. Scores
  from different alignments or site orders are not comparable; mismatch → hard error.
- If `--model-names` is provided, the label count must match the number of model
  groups and labels must be unique.
- Duplicate basenames within a group are disambiguated with a numeric suffix
  (`chain1_1.sitelogl`, `chain1_2.sitelogl`).

### Computation

Per model (Lartillot 2023):
- **LOO-CV:** per-run mean of site `logcpo`, then cross-run mean; debiased as
  `score − 0.5 × mean(site-wise variance of logcpo across runs)`.
- **wAIC:** per-run `mean(logpostmeanl) − mean(var)`, then cross-run mean;
  debiased as `score + 0.5 × mean(site-wise variance of logpostmeanl across runs)`.
- **ESS quality:** `%(ess<10)` (fraction of sites with ESS < 10) and
  `f(ess<10)` (fraction of score from such sites). Quality is `good`
  (max < 0.1), `ok` (max < 0.3), or `no` (max ≥ 0.3).
- **Confidence interval:** Student's t critical value for `n_runs − 1` df
  (exact table values for df 1–30; linear interpolation toward 1.96 beyond).
- **Δ values (≥2 models):** each metric independently selects its best model
  (highest debiased score); Δ = model score − best score, so the best model has
  Δ = 0 and all others ≤ 0.

### Outputs

```
runs/posttree/modelcompare/pb/
├── result.json
├── model_fit.csv
└── sitelogl/
    ├── model_1/
    │   ├── chain1.sitelogl
    │   └── chain2.sitelogl
    └── model_2/
        ├── chain1.sitelogl
        └── chain2.sitelogl
```

`model_fit.csv` uses a per-metric format (Metric, Score, Bias, StDev,
CI95_min, CI95_max, ESS, Pct_ESS_lt10, Frac_ESS_lt10, Quality) for a single
model, or a per-model wide format including Delta_LOOCV and Delta_wAIC for
multiple models.

`result.json` `key_results` exposes both `best_loocv_quality` and
`best_waic_quality` (each: `good` / `ok` / `no`).

### Examples

```bash
# Two model directories, each with >= 2 chain files
phyloai posttree modelcompare pb --sitelogl-dir model1,model2 --model-names CAT,GTR

# Explicit chain groups
phyloai posttree modelcompare pb \
  --sitelogl model1/c1.sitelogl,model1/c2.sitelogl \
  --sitelogl model2/c1.sitelogl,model2/c2.sitelogl
```

---

## Shared Notes

- `--seq-type auto` reads the alignment to detect AA vs NT *before* model
  validation; `--heterogeneous-model` is rejected for NT data. An explicit
  `--seq-type AA|NT` is cross-checked against the detected type (mismatch →
  error).
- IQ-TREE output files for `iqtree` are placed in an `iqtree/` subdirectory;
  IQ-TREE stdout streams to the terminal during execution.
- `--dry-run` prints the IQ-TREE command and validates inputs without running
  external tools.
- Output directories must be empty (or use `--overwrite`); `pb` refuses a
  non-empty directory without `--overwrite`.
- Relative model fit identifies the better candidate among those supplied; it
  does not establish absolute adequacy or topology correctness. For model
  adequacy and optional posterior-predictive follow-up, see the
  [systematic-error workflow reference](../../skills/phyloai-workflow/references/syserror-workflow.md).
- References: Kalyaanamoorthy et al. (2017) *Nature Methods*; Lartillot (2023)
  *Systematic Biology* 72(3):616–638; Watanabe (2009) *JMLR*.

## Exit Codes

| Code | Meaning |
|------|---------|
| 0 | Success |
| 1 | User input error (missing files, invalid parameters, output conflict) |
| 2 | External tool execution failed |
| 3 | External tool executable not found |
