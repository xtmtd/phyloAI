# PhyloAI Posttree Model-Compare Design Specification

**Date:** 2026-08-02
**Status:** Draft — pending approval
**Parent spec:** `2026-06-07-phyloai-design.md`
**References:**
- Kalyaanamoorthy et al. 2017, Nature Methods (ModelFinder: BIC/AIC model selection)
- Lartillot N. 2023. Identifying the Best Approximating Model in Bayesian Phylogenetics: Bayes Factors, Cross-Validation or wAIC? Systematic Biology 72(3):616–638.
- IQ-TREE3 Command Reference: `-m MF`, `-mset`, `-mrate`, `-madd`, `-cmin`, `-cmax`
- PhyloBayes-MPI 1.9 Manual §4 (Model comparison and selection: LOO-CV / wAIC)

---

## 1. Purpose

`phyloai posttree modelcompare` estimates the relative fit of alternative substitution models via two independent subcommands:

| Subcommand | Analysis | Core tool |
|------------|----------|-----------|
| `iqtree` | ModelFinder BIC/AIC/AICc comparison of homogeneous + heterogeneous models | IQ-TREE3 `-m MF` |
| `pb` | Leave-one-out cross-validation (LOO-CV) and wAIC from PhyloBayes sitelogl | Pure Python (integrated) |

This module answers "which model fits relatively better?" by comparing multiple candidate models. It is strictly about relative model comparison, not model adequacy.

**Future module:** `posttree simulate` will be a group with at least two subcommands: `alisim` (AliSim sequence simulation and gene-jackknife) and `adequacy` (posterior predictive adequacy tests: PPA-DIV/VAR/CONV). The `adequacy` subcommand answers "does the model fit well enough?" — complementary to this module's "which model fits relatively better?"

---

## 2. Design Principles

1. **Reuse `tree ml iqtree` infrastructure.** IQ-TREE path resolution (`_resolve_iqtree_path`), version detection (`_detect_iqtree_version`), file extension validation (`IQTREE_COMPATIBLE_EXTENSIONS`) are shared from `phyloai.core.iqtree`.
2. **Separate subcommands.** `iqtree` requires external tool execution; `pb` is pure computation on existing `.sitelogl` files. Different input requirements, different outputs.
3. **Heterogeneous model expansion is PhyloAI's job.** IQ-TREE's `-madd` does not interact with `-mrate`/`-mfreq`/`-cmin`/`-cmax`. PhyloAI expands `--heterogeneous-model` × `--het-mrate` × `{∅, +F}` into the full `-madd` string, ensuring rate categories lock at 4.
4. **Real-time terminal output.** The IQ-TREE subprocess stdout is printed directly to terminal (`stdout=None`, inheriting the parent stdout) so users can observe ModelFinder progress (which can be slow). This is the same pattern used by `posttree topology` and `posttree signal lnl`.
5. **LOO-CV/wAIC math is integrated.** No external Python script or subprocess; the `readwaic()` logic is embedded directly in PhyloAI, eliminating the `basicstat.py`/`studentC.py` dependency.
6. **Quality indicators use the same good/ok/no convention** as `tree bi pb` convergence diagnostics.
7. **Multi-model comparison for `pb`.** The `pb` subcommand supports comparing multiple models in a single invocation by accepting multiple sitelogl groups (one group per model).
8. **Seq-type detection before model validation.** When `--seq-type auto`, the alignment file is read first to determine molecule type; only then are `--homogeneous-model` and `--heterogeneous-model` validated against the detected type. An explicit `--seq-type` is cross-checked against the detected type (mismatch → hard error).

---

## 3. CLI Surface

```bash
# IQ-TREE ModelFinder model comparison
phyloai posttree modelcompare iqtree \
  --matrix concat.aa.fa \
  --homogeneous-model LG,WAG \
  [--mrate E,G] \
  [--heterogeneous-model C10,C20] \
  [--het-mrate G,R] \
  [--seq-type AA|NT|auto] \
  [--prefix modelcompare] \
  [--output-dir runs/posttree/modelcompare/iqtree] \
  [--threads auto] \
  [--iqtree-path /path/to/iqtree3] \
  [--tool-args "..."] \
  [--overwrite] [--resume] [--dry-run] [--quiet]

# PhyloBayes LOO-CV / wAIC — multi-model comparison
# Mode 1: multiple directories (one per model)
phyloai posttree modelcompare pb \
  --sitelogl-dir dir_CAT,dir_GTR,dir_LG \
  [--model-names CAT,GTR,LG] \
  [--output-dir runs/posttree/modelcompare/pb] \
  [--overwrite] [--quiet]

# Mode 2: multiple --sitelogl flags (each = one model's chains)
phyloai posttree modelcompare pb \
  --sitelogl chain1_cat.sitelogl,chain2_cat.sitelogl \
  --sitelogl chain1_gtr.sitelogl,chain2_gtr.sitelogl \
  --sitelogl chain1_lg.sitelogl,chain2_lg.sitelogl \
  [--model-names CAT,GTR,LG] \
  [--output-dir runs/posttree/modelcompare/pb] \
  [--overwrite] [--quiet]
```

### 3.1 Command Hierarchy

```
phyloai posttree
└── modelcompare
    ├── iqtree
    └── pb
```

---

## 4. Parameters

### 4.1 `modelcompare iqtree` Parameters

| Flag | Type | Required | Default | Description |
|------|------|----------|---------|-------------|
| `--matrix` | Path | yes | — | Single alignment file (FASTA/PHYLIP/NEXUS). Maps to IQ-TREE `-s`. |
| `--homogeneous-model` | str | yes | — | Comma-separated standard models for homogeneous model search space. Maps to IQ-TREE `-mset`. |
| `--mrate` | str | no | `E,G` | Rate heterogeneity types for homogeneous models. Maps to IQ-TREE `-mrate`. Valid: any combination of `E`, `G`, `R` (comma-separated). |
| `--heterogeneous-model` | str | no | — | Comma-separated AA mixture models for `-madd` expansion. AA only; rejected when seq-type is NT. |
| `--het-mrate` | str | no | `E,G` | Rate heterogeneity for heterogeneous model expansion. Each token selects a variant family: `E` = base models (`C10,C10+F`), `G` = `+G4`, `R` = `+R4`. Valid: any subset of `E`, `G`, `R`. AA only. Mirrors `--mrate` semantics. |
| `--seq-type` | choice | no | `auto` | Sequence type: `AA`, `NT`, `auto`. When `auto`, alignment is read to detect molecule type before model validation. |
| `--prefix` | str | no | `modelcompare` | IQ-TREE output prefix. Must be a single filename (no `/`, `..`, or absolute paths) so output stays inside the run directory. Maps to `--prefix`. |
| `--output-dir` | Path | no | `runs/posttree/modelcompare/iqtree` | Output directory. |
| `--threads` | str | no | `auto` | IQ-TREE `-T` value (integer or `auto`). |
| `--iqtree-path` | str | no | — | Custom path to iqtree3 executable. |
| `--tool-args` | str | no | — | Extra IQ-TREE flags. Blocked: `-s`, `--prefix`. Overrides PhyloAI-managed flags (e.g. `-madd`, `-mrate`, `-T`) when present. |
| `--overwrite` | flag | no | False | Delete and recreate output directory. |
| `--resume` | flag | no | False | Resume incomplete IQ-TREE run (native checkpoint). |
| `--dry-run` | flag | no | False | Print IQ-TREE command without executing. |
| `--quiet` | flag | no | False | Suppress terminal output except errors. |

**Valid `--homogeneous-model` values:**

- AA standard models: `LG`, `Poisson`, `cpREV`, `mtREV`, `Dayhoff`, `mtMAM`, `JTT`, `WAG`, `mtART`, `mtZOA`, `VT`, `rtREV`, `DCMut`, `PMB`, `HIVb`, `HIVw`, `JTTDCMut`, `FLU`, `Blosum62`, `GTR20`, `mtMet`, `mtVer`, `mtInv`, `FLAVI`, `Q.LG`, `Q.pfam`, `Q.pfam_gb`, `Q.bird`, `Q.mammal`, `Q.insect`, `Q.plant`, `Q.yeast`
- NT standard models: `GTR`, `HKY`, `JC`, `F81`, `K2P`, `K3P`, `K81uf`, `TN`, `TNef`, `TIM`, `TIMef`, `TVM`, `TVMef`, `SYM`

**Valid `--heterogeneous-model` values (AA only):**

`C10`, `C20`, `C30`, `C40`, `C50`, `C60`, `EX2`, `EX3`, `EHO`, `UL2`, `UL3`, `EX_EHO`, `LG4M`, `LG4X`

**Validation rules:**
- `--matrix` must exist and have a supported extension.
- `--homogeneous-model` is required.
- `--heterogeneous-model` requires detected/specified seq-type to be AA; hard error if NT.
- `--homogeneous-model` values are validated against the AA or NT model set based on detected/specified seq-type.
- `--mrate` valid values: any non-empty subset of {E, G, R}.
- `--het-mrate` valid values: any subset of `E` (base models only), `G` (`+G4`), `R` (`+R4`), comma-separated.
- `--threads` must be a positive integer or `auto`.
- `--prefix` must be a single filename (no `/`, `..`, or absolute paths) — it is passed to IQ-TREE and used to build output paths, so unsafe prefixes are rejected to keep reports/logs/checkpoints inside the `iqtree/` subdirectory.
- `--overwrite` and `--resume` are mutually exclusive.
- `--tool-args` must not contain blocked flags (`-s`, `--prefix`). When `--tool-args` supplies a flag PhyloAI also manages (`-m`, `-mset`, `-mrate`, `-madd`, `-cmin`, `-cmax`, `-T`), the `--tool-args` value overrides the PhyloAI-generated one (no duplicate flags in the command).
- **Seq-type cross-validation:** When `--seq-type` is explicitly `AA` or `NT`, the alignment is still read to detect the actual molecule type; a mismatch → hard error (`--seq-type NT` with an AA matrix).

**Validation order:**
1. Verify `--matrix` exists and has valid extension
2. Read alignment file to detect molecule type (AA vs NT). For `--seq-type auto`, this is the resolved type; for explicit `--seq-type`, it is cross-checked against the detected type (mismatch → hard error)
3. Validate `--homogeneous-model` against the resolved type's model set
4. If `--heterogeneous-model` provided and resolved type is NT → hard error
5. Proceed with remaining validation

### 4.2 `modelcompare pb` Parameters

| Flag | Type | Required | Default | Description |
|------|------|----------|---------|-------------|
| `--sitelogl-dir` | Path | mutex | — | Comma-separated directories; each directory represents one model, globbing `*.sitelogl` (≥2 files per directory). Click `Path` type → shell path completion. |
| `--sitelogl` | Path (multiple) | mutex | — | Repeatable option; each occurrence is comma-separated `.sitelogl` file paths for one model (≥2 files per model). Use multiple `--sitelogl` flags for multiple models. Click `Path` type → shell path completion. |
| `--model-names` | str | no | — | Comma-separated model labels matching the number of model groups. If omitted, models are labeled `model_1`, `model_2`, etc. Each label must be a single path component (no `/`, `..`) — it names the output subdirectory under `sitelogl/`, so unsafe labels are rejected to keep copies inside the run directory. |
| `--output-dir` | Path | no | `runs/posttree/modelcompare/pb` | Output directory. |
| `--overwrite` | flag | no | False | Delete and recreate output directory. |
| `--quiet` | flag | no | False | Suppress terminal output except errors. |

**Validation rules:**
- Exactly one of `--sitelogl-dir` or `--sitelogl` must be provided.
- At least 1 model group is required; at least 2 `.sitelogl` files per model group.
- All `.sitelogl` files within a model group must have the same number of data rows (sites).
- **Within-group site validation:** All `.sitelogl` files within a model group must share identical ordered `site` identifiers (not just the same count). Sites are merged by row index; mismatched ordering would silently corrupt LOO-CV/wAIC → hard error.
- **Cross-model site validation:** When ≥2 model groups are provided, all groups must have the same site count and identical ordered `site` identifiers. Scores from different alignments or site orders are not comparable; mismatch → hard error.
- If `--model-names` provided, count must match number of model groups.
- Model labels must be unique and be safe path components (no `/`, `..`, or absolute paths) — they become subdirectory names under `sitelogl/`, so unsafe labels are rejected to keep copied files inside the run directory.
- **Duplicate basenames:** Repeated `.sitelogl` basenames within a model group are allowed. The copy step disambiguates them with a numeric suffix (`chain1_1.sitelogl`, `chain1_2.sitelogl`).
- **Output-directory conflict:** Without `--overwrite`, a non-empty existing output directory → hard error (mirrors `iqtree` and `signal` commands). `--overwrite` deletes and recreates the directory.

---

## 5. Computation

### 5.1 `modelcompare iqtree`

**IQ-TREE command assembly:**
```
iqtree3 -s <matrix> -m MF -mset <homogeneous_models> -mrate <mrate> -cmin 4 -cmax 4 [-madd <expanded_heterogeneous>] --prefix <prefix> -T <threads> [tool_args]
```

**Auto-set flags:**
- `-m MF` — ModelFinder only (no tree search beyond BIC evaluation)
- `-cmin 4 -cmax 4` — lock FreeRate categories to exactly 4 (so `+Rn` becomes `+R4`)

**Heterogeneous model expansion algorithm:**

Given `--heterogeneous-model C10,C20 --het-mrate G,R`:

For each model M in heterogeneous list:
For each rate token R in het-mrate:
- If R == "E": append `M`, `M+F` (empirical frequencies, no rate category; mirrors `--mrate E`)
- If R == "G": append `M+G4`, `M+F+G4`
- If R == "R": append `M+R4`, `M+F+R4`

Only the requested token families are produced (default `E,G` yields all four).

Result for C10 with G,R: `C10+G4,C10+F+G4,C10+R4,C10+F+R4`
Result for C10 with E: `C10,C10+F`
Result for C10 with E,G: `C10,C10+F,C10+G4,C10+F+G4`

All expanded models are joined comma-separated → passed to `-madd`.

**Result parsing:**

Parse `<prefix>.iqtree` for "List of models sorted by BIC scores:" section:
- Header line: `Model  LogL  AIC  w-AIC  AICc  w-AICc  BIC  w-BIC`
- Data lines until blank line
- Each criterion (AIC, AICc, BIC) independently marks its 95% confidence set members with `+` and excluded models with `-` before the weight value
- Best model by BIC: first row in BIC-sorted table (lowest BIC)

Output: `model_fit.csv` with columns:
```
Rank,Model,LogL,AIC,w_AIC,In_AIC_95,AICc,w_AICc,In_AICc_95,BIC,w_BIC,In_BIC_95
```

`In_AIC_95`, `In_AICc_95`, `In_BIC_95`: boolean columns indicating membership in each criterion's 95% confidence set (`+` in IQ-TREE output).

### 5.2 `modelcompare pb`

**Workflow:**
1. Validate and collect `.sitelogl` file groups (one group per model, ≥2 files per group)
2. Site validation: within each group, all files must share identical ordered `site` identifiers; if ≥2 model groups, all groups must have identical site counts and ordered `site` identifiers (mismatch → hard error with descriptive message)
3. Output-directory conflict check: non-empty existing output dir without `--overwrite` → hard error
4. Copy to `<output_dir>/sitelogl/model_1/`, `model_2/`, etc. If duplicate basenames exist within a group, append numeric suffix (`chain1_1.sitelogl`, `chain1_2.sitelogl`)
5. For each model group: parse files, compute LOO-CV and wAIC using integrated logic (§5.3)
6. If ≥2 models: compute Δ values relative to best model
7. Write `model_fit.csv` + print terminal table

### 5.3 LOO-CV/wAIC computation (per Lartillot 2023)

**LOO-CV (per model):**
- For each run r: `loocv_r = mean(logcpo across sites)`
- Cross-run mean: `loocv = mean(loocv_r for all runs)`
- Cross-run bias: `bias = 0.5 * mean(site-wise variance of logcpo across runs)`
- Debiased score: `loocv - bias`

**wAIC (Watanabe 2009, per model):**
- For each run r: `waic_r = mean(logpostmeanl) - mean(var)`
- Cross-run mean: `waic = mean(waic_r for all runs)`
- Cross-run bias: `bias = -0.5 * mean(site-wise variance of logpostmeanl across runs)`
- Debiased score: `waic - bias`

**ESS quality indicators (per model):**
- `%(ess<10)`: fraction of sites with ESS < 10
- `f(ess<10)`: fraction of total score contributed by sites with ESS < 10

**Quality classification (per model):**
- `good`: `max(%(ess<10), f(ess<10)) < 0.1`
- `ok`: `max(%(ess<10), f(ess<10)) < 0.3`
- `no`: `max(%(ess<10), f(ess<10)) >= 0.3`

**Confidence interval:** Uses Student's t distribution critical value for `n_runs - 1` degrees of freedom. For df 1–30, exact table values are used. For df > 30, linear interpolation toward z_0.975 = 1.96: `t(df) = t30 + (1.96 - t30) * (1 - 30/df)`, which monotonically decreases from t30=2.042 toward 1.96 as df → ∞.

**Δ calculation (≥2 models):**
- LOO-CV and wAIC select their best models independently: each is the model with the highest (least negative) debiased score.
- `Δ_LOOCV = loocv_i - loocv_best` and `Δ_wAIC = waic_i - waic_best`; each metric's best model is 0 and all others are ≤ 0.

### 5.4 Output CSV formats

**Single model `model_fit.csv`:**
```
Metric,Score,Bias,StDev,CI95_min,CI95_max,ESS,Pct_ESS_lt10,Frac_ESS_lt10,Quality
LOO-CV,-27.6689,0.0107,0.0066,-27.7526,-27.5851,550.42,0.022,0.055,good
wAIC,-27.6800,-0.0016,0.0015,-27.6990,-27.6609,589.87,0.000,0.000,good
```

**Multi-model `model_fit.csv`:**
```
Model,LOO-CV,LOO-CV_Bias,LOO-CV_StDev,LOO-CV_CI95min,LOO-CV_CI95max,LOO-CV_ESS,LOO-CV_Pct_ESS_lt10,LOO-CV_Frac_ESS_lt10,LOO-CV_Quality,wAIC,wAIC_Bias,wAIC_StDev,wAIC_CI95min,wAIC_CI95max,wAIC_ESS,wAIC_Pct_ESS_lt10,wAIC_Frac_ESS_lt10,wAIC_Quality,Delta_LOOCV,Delta_wAIC
model_1,-27.64,0.035,0.009,-27.76,-27.52,58.4,0.091,0.200,ok,-27.66,-0.014,0.009,-27.77,-27.55,61.0,0.041,0.086,good,0.000,0.000
model_2,-28.31,0.005,0.009,-28.42,-28.20,69.7,0.003,0.007,good,-28.30,-0.003,0.007,-28.39,-28.21,72.9,0.002,0.003,good,-0.670,-0.643
```

---

## 6. Output Structure

### 6.1 `modelcompare iqtree`

```
runs/posttree/modelcompare/iqtree/
├── result.json
├── model_fit.csv
└── iqtree/
    ├── modelcompare.iqtree
    ├── modelcompare.log
    ├── modelcompare.model.gz
    └── modelcompare.treefile
```

### 6.2 `modelcompare pb`

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

---

## 7. result.json Schema

Conforms to the global JSON output standard (`2026-06-21-phyloai-json-output-standard.md`).

### 7.1 `modelcompare iqtree`

```json
{
  "status": "success|error",
  "command": "phyloai posttree modelcompare iqtree --matrix ... --homogeneous-model ...",
  "wall_time": 42.5,
  "tool_versions": {"iqtree3": "3.1.2"},
  "params": {
    "matrix": "/abs/path/to/matrix.fa",
    "homogeneous_model": "LG,WAG",
    "mrate": "E,G",
    "heterogeneous_model": "C10,C20",
    "het_mrate": "E,G",
    "seq_type": "AA",
    "detected_seq_type": "AA",
    "prefix": "modelcompare",
    "output_dir": "/abs/path/to/output",
    "threads": "auto",
    "iqtree_path": null,
    "tool_args": null,
    "overwrite": false,
    "resume": false,
    "dry_run": false,
    "quiet": false
  },
  "key_results": {
    "best_model_bic": "LG+G4",
    "best_model_aic": "LG+F+G4",
    "best_model_aicc": "LG+F+G4",
    "n_models_tested": 12,
    "madd_expanded": "C10,C10+F,C10+G4,C10+F+G4,C10+R4,C10+F+R4,C20,C20+F,C20+G4,C20+F+G4,C20+R4,C20+F+R4"
  },
  "error": null,
  "error_category": null,
  "data": {
    "cmd": ["iqtree3", "-s", "...", "-m", "MF", "..."],
    "tool_stderr": "",
    "tool_log": "/abs/path/to/iqtree/modelcompare.log",
    "output_files": {
      "model_fit_csv": {"path": "...", "description": "Model comparison table sorted by BIC"},
      "iqtree_report": {"path": "...", "description": "IQ-TREE native report"},
      "iqtree_log": {"path": "...", "description": "IQ-TREE console log"}
    },
    "models": [
      {"rank": 1, "model": "LG+G4", "logl": -2680.063, "bic": 5414.722, "w_bic": 1.0, "in_bic_95": true, "in_aic_95": false, "in_aicc_95": false},
      {"rank": 2, "model": "LG+R4", "logl": -2677.000, "bic": 5435.894, "w_bic": 2.53e-05, "in_bic_95": false, "in_aic_95": false, "in_aicc_95": false}
    ]
  }
}
```

> **Params linkage:** `params.het_mrate` is `null` when `params.heterogeneous_model` is `null` (no heterogeneous search requested). They are always set together.


### 7.2 `modelcompare pb`

```json
{
  "status": "success|error",
  "command": "phyloai posttree modelcompare pb --sitelogl-dir dir1,dir2",
  "wall_time": 0.3,
  "tool_versions": {},
  "params": {
    "sitelogl_dir": "dir1,dir2",
    "sitelogl": null,
    "model_names": "model_1,model_2",
    "output_dir": "/abs/path/to/output",
    "overwrite": false,
    "quiet": false
  },
  "key_results": {
    "n_models": 2,
    "n_sites": 235,
    "best_model_loocv": "model_1",
    "best_model_waic": "model_1",
    "best_loocv_score": -27.64,
    "best_loocv_quality": "ok",
    "best_waic_score": -27.66,
    "best_waic_quality": "good"
  },
  "error": null,
  "error_category": null,
  "data": {
    "output_files": {
      "model_fit_csv": {"path": "...", "description": "LOO-CV and wAIC scores per model with Δ values"}
    },
    "models": [
      {
        "model": "model_1",
        "n_runs": 3,
        "loocv": {"score": -27.64, "bias": 0.035, "stdev": 0.009, "ci95_min": -27.76, "ci95_max": -27.52, "ess": 58.4, "pct_ess_lt10": 0.091, "frac_ess_lt10": 0.200, "quality": "ok"},
        "waic": {"score": -27.66, "bias": -0.014, "stdev": 0.009, "ci95_min": -27.77, "ci95_max": -27.55, "ess": 61.0, "pct_ess_lt10": 0.041, "frac_ess_lt10": 0.086, "quality": "good"},
        "delta_loocv": 0.0,
        "delta_waic": 0.0
      },
      {
        "model": "model_2",
        "n_runs": 2,
        "loocv": {"score": -28.31, "bias": 0.005, "stdev": 0.009, "ci95_min": -28.42, "ci95_max": -28.20, "ess": 69.7, "pct_ess_lt10": 0.003, "frac_ess_lt10": 0.007, "quality": "good"},
        "waic": {"score": -28.30, "bias": -0.003, "stdev": 0.007, "ci95_min": -28.39, "ci95_max": -28.21, "ess": 72.9, "pct_ess_lt10": 0.002, "frac_ess_lt10": 0.003, "quality": "good"},
        "delta_loocv": -0.670,
        "delta_waic": -0.643
      }
    ]
  }
}
```

---

## 8. Terminal Output

### 8.1 `modelcompare iqtree`

IQ-TREE stdout printed to terminal in real-time via `Popen(stdout=None)`. After completion:

```
Model comparison (sorted by BIC):

  Rank  Model        LogL        BIC       w-BIC  In_BIC_95
     1  LG+G4      -2680.063   5414.722   1.000   +
     2  LG+R4      -2677.000   5435.894   0.000   -
     3  LG+F+G4    -2644.688   5447.704   0.000   -
     ...

Best model (BIC): LG+G4
Models tested: 12
Result written to runs/posttree/modelcompare/iqtree/result.json
```

### 8.2 `modelcompare pb` (multi-model)

```
LOO-CV / wAIC model comparison:

  Model      LOO-CV     LOO-CV_Q  wAIC       wAIC_Q    Δ_LOOCV   Δ_wAIC
  model_1    -27.6400   ok        -27.6600   good       0.000     0.000
  model_2    -28.3100   good      -28.3000   good      -0.670    -0.643

Best model (LOO-CV): model_1
Sites: 235 | Models: 2
Result written to runs/posttree/modelcompare/pb/result.json
```

(Full ESS detail columns are in `model_fit.csv`; terminal table shows condensed summary.)

### 8.3 `modelcompare pb` (single model)

```
LOO-CV / wAIC model fit:

  Metric   Score      Bias       StDev    CI95min    CI95max    ESS     %(ess<10)  f(ess<10)  Quality
  LOO-CV   -27.6689   0.0107     0.0066   -27.7526   -27.5851   550.42  0.022      0.055      good
  wAIC     -27.6800   -0.0016    0.0015   -27.6990   -27.6609   589.87  0.000      0.000      good

Sites: 627 | Runs: 2
Result written to runs/posttree/modelcompare/pb/result.json
```

---

## 9. Report Template

### 9.1 `posttree.modelcompare.iqtree`

> Relative model fit was assessed using ModelFinder (Kalyaanamoorthy et al. 2017) as implemented in IQ-TREE3 v{version}. The homogeneous model search space included {homogeneous_models} with rate heterogeneity types {mrate}. {het_sentence}A total of {n_models_tested} model configurations were evaluated. The best-fitting model according to BIC was {best_model} (BIC = {bic}, w-BIC = {w_bic}).

Where `{het_sentence}` (only if `--heterogeneous-model` provided):
> Heterogeneous mixture models ({heterogeneous_models}) were additionally evaluated {with rate variants {het_mrate} }{yielding {n_madd} expanded model configurations passed via -madd}.

### 9.2 `posttree.modelcompare.pb`

**Multi-model:**
> Relative model fit was evaluated using leave-one-out cross-validation (LOO-CV) and the widely applicable information criterion (wAIC) following Lartillot (2023), computed from site log-likelihood files ({n_sites} sites). {n_models} candidate models were compared. The best-fitting model according to LOO-CV was {best_model_loocv} (LOO-CV = {score:.4f}, Δ = 0; quality: {quality}). The best-fitting model according to wAIC was {best_model_waic} (wAIC = {score:.4f}).

**Single model:**
> Model fit was evaluated using leave-one-out cross-validation (LOO-CV) and the widely applicable information criterion (wAIC) following Lartillot (2023), computed from {n_runs} independent MCMC chain site log-likelihood files ({n_sites} sites). The debiased LOO-CV score was {loocv_score:.4f} (quality: {loocv_quality}; ESS = {ess:.1f}; %(ESS<10) = {pct:.3f}; f(ESS<10) = {frac:.3f}). The debiased wAIC score was {waic_score:.4f} (quality: {waic_quality}).

---

## 10. MCP Integration

Both subcommands auto-generate MCP tools via the existing `schema_gen.py` Click-tree walker:
- `posttree_modelcompare_iqtree`
- `posttree_modelcompare_pb`

No manual stub required.

---

## 11. Dependencies

- **`modelcompare iqtree`:** IQ-TREE3 (external binary)
- **`modelcompare pb`:** None (pure Python — math is stdlib `math` module only)
- **No new pip dependencies**

---

## 12. Test Data

- IQ-TREE: `runs/modelCompare/EOG090X0A0V.fa` (6 taxa, 235 AA sites)
- PhyloBayes: `runs/modelCompare/LOOCV_wAIC/chain{1,2,3}.sitelogl` (235 sites, 3 chains)
