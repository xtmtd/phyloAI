# PhyloAI Posttree Simulate Adequacy Design Specification

**Date:** 2026-08-03
**Status:** Draft — pending approval
**Parent spec:** `2026-06-07-phyloai-design.md`
**JSON standard:** `2026-06-21-phyloai-json-output-standard.md`
**References:**
- Lartillot et al. (2007) PhyloBayes posterior predictive checks
- Blanquart & Lartillot (2006) PPA-COMP
- PhyloBayes-MPI 1.9 manual section 5–6 (`ref/pb_mpiManual1.9.md`)
- PhyloBayes-MPI source: `~/tools/pbmpi-1.9/sources/SequenceAlignment.h`
- Reference script: `ref/scripts/calculate_zscore.py`

---

## 1. Purpose

`phyloai posttree simulate adequacy` evaluates model adequacy by comparing four
summary statistics computed on an observed MSA against the null distribution
derived from a directory of simulated MSAs (e.g. PhyloBayes `readpb_mpi -ppred`
replicates, or AliSim replicates from `phyloai posttree simulate alisim iqtree`).

The four statistics faithfully replicate the PhyloBayes `SequenceAlignment.h`
formulas and sign conventions:

| Short name | PB flag | Full name | What it captures |
|---|---|---|---|
| `div` | `-div` | Mean diversity | Mean distinct states per site (PPA-DIV) |
| `siteconvprob` | `-siteconvprob` | Convergence probability | Mean Σf² across states/sites (PPA-CONV) |
| `sitecomp` | `-sitecomp` | Across-site comp. heterogeneity | Mean variance of site-specific freqs (PPA-VAR) |
| `comp` | `-comp` | Compositional homogeneity across taxa | Squared deviation taxon vs. global freq (PPA-COMP) |

`comp` yields two scalars (max and mean heterogeneity across taxa) plus a
per-taxon breakdown. No external tools are required; the module is pure Python
+ BioPython.

---

## 2. Design Principles

1. **No new dependencies.** BioPython (MSA reading) + stdlib only (`statistics`,
   `math`, `concurrent.futures`).
2. **Faithful formulas.** All four statistics derived from `SequenceAlignment.h`
   with matching sign conventions for z-score and pp.
3. **Reuse existing infrastructure.** Format detection via
   `phyloai.core.formats.detect_alignment_format`; checkpoint/resume via
   `phyloai.core.checkpoint`; progress via Rich; output via standard
   `result.json` schema.
4. **Per-simulation transparency.** `per_simulation_stats.csv` records every
   simulated file's raw statistic values.
5. **Empirical interval.** 95% interval = bounded empirical p2.5 / p97.5
   percentiles via `statistics.quantiles(data, n=40, method="inclusive")`.
   This interpolation never extrapolates beyond simulated values, including at
   the minimum of 10 replicates. No normality assumption.
6. **Checkpoint/resume.** Each simulated file's computed statistics are
   persisted in `checkpoint.json` so an interrupted run skips already-completed
   files. Useful for large MSAs (>1000 taxa or >50k sites) with many replicates.

---

## 3. CLI Surface

```bash
phyloai posttree simulate adequacy \
   --original-msa PATH        # observed MSA (required)
   --simulated-dir PATH       # directory of simulated MSAs (required)
   [--seq-type AA|NT|auto]    # default: auto
   [--threads INT]            # default: 4
   [--table-format csv|tsv]   # default: csv
  [--output-dir PATH]        # default: runs/posttree/simulate/adequacy
  [--overwrite]
  [--resume]
  [--dry-run]
  [--quiet]
```

### Parameter details

| Option | Required | Default | Description |
|---|---|---|---|
| `--original-msa` | Yes | — | Observed alignment. Formats: `fasta`, `phylip-relaxed`, `phylip-paml`, `nexus`. Auto-detected by content if extension is unknown. |
| `--simulated-dir` | Yes | — | Directory of simulated MSAs. All regular non-empty files are scanned; format auto-detected per file (supports `.ali` via content sniffing). |
| `--seq-type` | No | `auto` | Override sequence-type detection. `auto` infers from character content of `--original-msa`. |
| `--threads` | No | `4` | Parallel workers for per-simulation statistics computation. `--original-msa` is always processed first in the main thread. |
| `--table-format` | No | `csv` | Delimiter and suffix for the three output tables: `csv` or `tsv`. |
| `--output-dir` | No | `runs/posttree/simulate/adequacy` | Output directory. |
| `--overwrite` | No | `False` | Delete and recreate output directory. Mutually exclusive with `--resume`. |
| `--resume` | No | `False` | Resume from `checkpoint.json`. Skips simulated files already computed. |
| `--dry-run` | No | `False` | Validate inputs, detect seq-type, count simulated files; write nothing. |
| `--quiet` | No | `False` | Suppress Rich terminal output except errors. |

---

## 4. Statistics: Formulas and Sign Conventions

**Valid characters:** PhyloBayes encodes all gap and ambiguity characters as
`unknown = -1` (`BiologicalSequences.h:50`) and excludes them uniformly in all
statistic functions via `if (state != unknown)`. This implementation matches
that behaviour: for AA the valid set is `ACDEFGHIKLMNPQRSTVWY`; for NT the
valid set is `ACGT`. All other characters (gaps, ambiguity codes, `X`, `N`,
etc.) are treated as missing. Sites with zero valid observations are excluded
from all averages (`n_informative_sites`).

**Standard deviation convention:** All SD values use the **population SD**
(`sqrt(E[x²] − E[x]²)`), matching PhyloBayes `AllPostPred`
(`PhyloProcess.cpp:2225`: `varstatarray[k] /= samplesize; varstatarray[k] -= mean²`).
Minimum valid simulations: **≥ 10** (fewer raises a `ValueError` before aggregation).
When `sd_sim == 0` (all simulated values identical), `z_score` is set to `0.0`
and `pp` is set to `null` (JSON) / empty string (CSV) to signal an
uninformative null distribution; both `ci_lower` and `ci_upper` equal
`mean_sim`. Implementations must not write `NaN` to JSON (not valid JSON) nor
to CSV without quoting; use the sentinel values above.

**Input validation:** Before any statistic is computed:

1. **Unique taxon IDs in original MSA** — if any taxon name appears more than
   once, raise `ValueError` immediately (before scanning simulated files).
2. **Each simulated MSA** must satisfy — otherwise that file is skipped with a
   warning and counted in `n_failed`:
   a. No duplicate taxon IDs within the file.
   b. Taxon set identical to `--original-msa` (same names, no extras, no missing).
   c. Alignment length identical to `--original-msa`.
3. Taxon order is re-matched by name before `comp` per-taxon statistics, so
   file order does not need to match.

Simulated files failing these checks are skipped with a warning and counted in
`n_failed`; if all files fail, the run exits with an error.

**All-missing taxon:** If any taxon has zero valid characters across all sites
in the original MSA, the run exits with an error (`ValueError`). If a taxon has
zero valid characters in a simulated MSA, that file is skipped with a warning.

### 4.1 Diversity (`div`) — PPA-DIV

Mean number of distinct valid states per site:

```
div = (1 / S) × Σ_i  |distinct_states(site_i)|
```

where `S` = number of informative sites.

- **Z-score:** `z = (mean_sim − obs) / sd_sim`
  Positive → observed diversity lower than simulated (model over-predicts).
- **pp:** `P(sim ≤ obs)` = fraction of replicates where `sim_val ≤ obs`.
  Source: `PhyloProcess.cpp:2237` `pp = 1 - ppstatarray[0]` where
  `ppstatarray[0]` accumulates `sim > obs`.
  Low pp → observed diversity unexpectedly low (model over-predicts).

### 4.2 Convergence probability (`siteconvprob`) — PPA-CONV

Mean squared empirical frequency summed over all states, averaged over sites:

```
siteconvprob = (1 / S) × Σ_i Σ_k  f_ik²
```

where `f_ik` = empirical frequency of state `k` at site `i` (valid chars only).

- **Z-score:** `z = (obs − mean_sim) / sd_sim`
  Positive → data more convergent than model predicts.
- **pp:** `P(sim > obs)` = fraction of replicates where `sim_val > obs`.
  Source: `PhyloProcess.cpp:2245` `pp = ppstatarray[1]` where
  `ppstatarray[1]` accumulates `sim > obs`.
  Low pp → observed convergence unexpectedly high (anomalously convergent data).

### 4.3 Across-site compositional heterogeneity (`sitecomp`) — PPA-VAR

Mean variance of site-specific empirical frequencies across sites, averaged
over all states:

```
sitecomp = (1 / K) × Σ_k  Var_sites(f_k)
         = (1 / K) × Σ_k  [ E[f_k²] − E[f_k]² ]
```

where `K` = number of states (20 for AA, 4 for NT).

- **Z-score:** `z = (obs − mean_sim) / sd_sim`
  Positive → data more heterogeneous across sites than model predicts.
- **pp:** `P(sim > obs)` = fraction of replicates where `sim_val > obs`.
  Source: `PhyloProcess.cpp:2253` `pp = ppstatarray[2]`.
  Low pp → observed across-site heterogeneity unexpectedly high.

### 4.4 Compositional homogeneity across taxa (`comp`) — PPA-COMP

For each taxon `j`, compute the squared L2 distance between the taxon-specific
empirical frequency vector and the global (mean across taxa) frequency vector:

```
dist_j = Σ_k  (taxfreq_j[k] − globalfreq[k])²
```

Two scalars are reported:

| Sub-statistic | Formula | Interpretation |
|---|---|---|
| `comp_max` | `max_j(dist_j)` | Maximum heterogeneity: worst-offending taxon |
| `comp_mean` | `mean_j(dist_j)` | Mean squared heterogeneity across all taxa |

Both scalars:
- **Z-score:** `z = (obs − mean_sim) / sd_sim`
  Positive → data more compositionally heterogeneous than model predicts.
- **pp:** `P(sim > obs)` = fraction of replicates where `sim_val > obs`.
  Source: `PhyloProcess.cpp:2264,2272` `pp = ppstatarray[3/4]`.
  Low pp → observed compositional heterogeneity unexpectedly high.

Per-taxon statistics: for each taxon `j`, the distribution of `dist_j` over
simulated replicates yields `obs`, `mean_pred`, `sd_pred`, `ci_lower`,
`ci_upper`, `z_score`, and `pp`.
- **Per-taxon pp:** `P(sim_dist_j > obs_dist_j)` = fraction of replicates
  where taxon `j`'s `dist_j` exceeds the observed value.
  Source: `PhyloProcess.cpp:2283` `pp = pptaxstat[j]` accumulating `sim > obs`.

---

## 5. Output Files

```
<output_dir>/
├── adequacy_summary.csv         # one row per scalar statistic (5 scalars: div, siteconvprob, sitecomp, comp_max, comp_mean)
├── adequacy_taxon_comp.csv      # per-taxon comp statistics
├── per_simulation_stats.csv     # raw statistic values for all simulated files
├── checkpoint.json              # resume state
└── result.json                  # standard PhyloAI output
```

The three tabular outputs above use `.csv`; with `--table-format tsv` the same
three tables are tab-delimited with `.tsv` suffixes, and the matching paths are
declared in `result.json:data.output_files`.

### 5.1 `adequacy_summary.csv`

One row per scalar statistic. Columns:

| Column | Description |
|---|---|
| `statistic` | One of 5 values: `div`, `siteconvprob`, `sitecomp`, `comp_max`, `comp_mean` |
| `obs` | Value computed on `--original-msa` |
| `mean_sim` | Mean over simulated replicates |
| `sd_sim` | Standard deviation over simulated replicates |
| `ci_lower` | Empirical p2.5 percentile of simulated distribution |
| `ci_upper` | Empirical p97.5 percentile of simulated distribution |
| `z_score` | As defined in Section 4 (sign convention matches PhyloBayes) |
| `pp` | Posterior predictive p-value (fraction of replicates more extreme than obs) |
| `n_simulations` | Number of valid simulated replicates used |

Example (mirrors `chain1.ppred` from test data):

```
statistic,obs,mean_sim,sd_sim,ci_lower,ci_upper,z_score,pp,n_simulations
div,2.78723,2.80834,0.076548,...,...,0.2757,0.40,100
siteconvprob,0.554364,0.553529,0.015318,...,...,0.0545,0.49,100
sitecomp,0.024619,0.024464,0.000769,...,...,0.2024,0.40,100
comp_max,0.003965,0.003140,0.000834,...,...,0.9896,0.15,100
comp_mean,0.002321,0.001978,0.000388,...,...,0.8856,0.22,100
```

### 5.2 `adequacy_taxon_comp.csv`

One row per taxon. Columns:

| Column | Description |
|---|---|
| `taxon` | Taxon name from `--original-msa` |
| `obs` | `dist_j` computed on original MSA |
| `mean_pred` | Mean `dist_j` across simulated replicates |
| `sd_pred` | SD of `dist_j` across simulated replicates |
| `ci_lower` | Empirical p2.5 of simulated `dist_j` distribution |
| `ci_upper` | Empirical p97.5 of simulated `dist_j` distribution |
| `z_score` | `(obs − mean_pred) / sd_pred` |
| `pp` | Fraction simulated `dist_j > obs` (lower pp = more extreme) |

Example (mirrors taxon table in `chain1.ppred`):

```
taxon,obs,mean_pred,sd_pred,ci_lower,ci_upper,z_score,pp
Drosophila_melanogaster,0.002774,0.002677,...,...,...,0.1047,0.43
Maconellicoccus_hirsutus,0.003965,0.002401,...,...,...,1.9408,0.04
...
```

### 5.3 `per_simulation_stats.csv`

One row per simulated file. Columns: `file`, `div`, `siteconvprob`, `sitecomp`,
`comp_max`, `comp_mean`. Enables manual inspection of the null distribution.

### 5.4 `result.json`

Standard PhyloAI schema. Key fields:

```json
{
  "status": "success",
  "command": "phyloai posttree simulate adequacy ...",
  "wall_time": 12.4,
  "tool_versions": {},
  "error": null,
  "params": { ... },
  "key_results": {
    "n_simulations": 100,
    "seq_type": "AA",
    "n_taxa": 6,
    "n_sites": 235,
    "statistics": {
      "div":          { "obs": 2.787, "mean_sim": 2.808, "sd_sim": 0.077, "z_score": 0.276, "pp": 0.40 },
      "siteconvprob": { "obs": 0.554, "mean_sim": 0.554, "sd_sim": 0.015, "z_score": 0.054, "pp": 0.49 },
      "sitecomp":     { "obs": 0.025, "mean_sim": 0.024, "sd_sim": 0.001, "z_score": 0.202, "pp": 0.40 },
      "comp": {
        "max":  { "obs": 0.004, "mean_sim": 0.003, "sd_sim": 0.001, "z_score": 0.990, "pp": 0.15 },
        "mean": { "obs": 0.002, "mean_sim": 0.002, "sd_sim": 0.000, "z_score": 0.886, "pp": 0.22 }
      }
    }
  },
  "data": {
    "cmd": [],
    "tool_stderr": "",
    "output_files": {
      "adequacy_summary": {
        "path": "/abs/path/adequacy_summary.csv",
        "description": "Model adequacy summary: obs, mean_sim, sd_sim, CI, z-score, pp for 5 statistics"
      },
      "adequacy_taxon_comp": {
        "path": "/abs/path/adequacy_taxon_comp.csv",
        "description": "Per-taxon PPA-COMP statistics: obs, mean_pred, sd_pred, CI, z-score, pp"
      },
      "per_simulation_stats": {
        "path": "/abs/path/per_simulation_stats.csv",
        "description": "Raw statistic values for all simulated replicates (null distribution)"
      }
    }
  }
}
```

---

## 6. Sequence Type Detection

`--seq-type auto` inspects the first 200 valid characters of `--original-msa`.
Characters matching `EFILPQWYZ` unambiguously indicate AA. Otherwise NT is
assumed. This matches the existing convention in `phyloai.core.formats`.

`--seq-type` is accepted case-insensitively and normalized with `str.upper()`
at the entry point. Only the normalized `AA` or `NT` value is ever passed to a
statistic worker; the literal `auto` is resolved before dispatch and is never
used as a state set selector.

**Format detection:** The command intentionally does not expose a single
`--input-format`: the observed MSA and individual simulated MSAs may use
different supported formats. `FormatConverter` therefore detects each file
independently by content/extension. This is a scoped exception to the parent
design's shared `--input-format` convention; separate original/simulated format
flags would add interface complexity without improving mixed-format support.

---

## 7. Checkpoint / Resume

The checkpoint (`checkpoint.json`) follows `phyloai.core.checkpoint` schema.
One `CheckpointTask` is created per simulated file with:

- `task_id`: absolute path of the simulated file (stable, unique).
- `input`: `"<abs_path>|<file_size>|<mtime_ns>"` — path + size + mtime
  fingerprint. On resume, if the stored fingerprint does not match the current
  file, the task is treated as `pending` (force-recompute). This prevents
  silent reuse of stale statistics when a simulated file is replaced.
- `outputs`: dict storing **all values needed to reconstruct both output CSVs**
  without re-reading any simulated file. Keys:
  - `div`, `siteconvprob`, `sitecomp`, `comp_max`, `comp_mean` — five scalar
    stats as strings.
  - `taxon_dist_j` — JSON-encoded dict `{taxon_name: dist_j_float}` for all
    taxa, serialised as a single string value. Example:
    `'{"Taxon_A": "0.00277", "Taxon_B": "0.00174"}'`.

  All values are stored **as strings** (per
  `CheckpointTask.outputs: dict[str, str | None]` in `checkpoint.py:59`).
  The aggregation step must:
  - Cast scalar stats with `float(v)`.
  - Parse `taxon_dist_j` with `json.loads(v)` then cast each leaf to `float`.
  Implementations must not assume numeric types from the checkpoint.
- `status`: `pending | success | failed`.

**Resume fingerprint and membership handling:**

- The checkpoint also stores an `original_msa` fingerprint
  (`"<abs_path>|<size>|<mtime_ns>"`) captured at run start. On resume the
  observed MSA is fingerprinted again and must match, otherwise the resume is
  rejected with a clean input error: the persisted simulated distributions
  would no longer correspond to the observed data, so reusing them is unsafe.

  This fingerprint is stored as a dedicated `Checkpoint.original_msa_fingerprint`
  field — a metadata attribute **outside `checkpoint.params`**. It therefore
  never enters `validate_resume_params`' `canonical_params_hash`, which hashes
   only the `params` dict (checkpoint.py:154-175), so persisting it cannot
   invalidate an otherwise-valid resume. The field is optional and read with
   `data.get("original_msa_fingerprint")` (default `None`) for backward
   compatibility with existing checkpoints. A missing fingerprint is unsafe to
   resume and requires a fresh `--overwrite` run.
- On resume the simulated directory is **rescanned** and reconciled against the
  checkpoint:
  - Files present on disk but absent from the checkpoint become new `pending`
    tasks.
  - `success` tasks whose file still exists with a matching fingerprint are
    kept.
  - A checkpoint task whose file no longer exists on disk is dropped from the
    run with a warning and excluded from aggregation. A per-file `stat()` must
    not raise on a missing path; existence is checked before fingerprinting.

On `--resume`:
1. `validate_resume_params` checks that `original_msa`, `simulated_dir`, and
   requested `seq_type` match the checkpoint.
2. Require a non-null matching `original_msa` fingerprint; absence or mismatch
   → `PreflightError`.
3. Rescan `simulated_dir` and reconcile membership as above.
4. Tasks with `status == success` AND matching fingerprint are skipped.
5. All other tasks (added, stale, or invalidated) are recomputed as described.
6. Aggregation always re-runs from all current `success` task outputs.

**Pre-flight refusals write no `result.json` (scoped exception to the parent
JSON Output Standard §6.1).** The module must complete output lifecycle
preflight before parsing inputs or raising ordinary validation errors: reject
`--overwrite --resume`; reject a non-empty fresh output directory; and, for
`--resume`, load a present, schema-compatible checkpoint, validate resume
parameters, and validate a non-null matching observed-MSA fingerprint. Only
after those checks may the invocation claim the output directory and use the
standard error-result path. Any earlier failure, including a missing or
unparseable observed MSA when the directory has not been claimed, exits with
code 1 after printing the error to stderr and writes no `result.json`.

Preflight failures use a dedicated `PreflightError(ValueError)` subclass; the
CLI catches it before the general `ValueError` handler and calls `_fail`
directly. This prevents writing into an unrelated directory or a prior run
that failed resume validation. A malformed, unsupported-schema, or
missing-fingerprint checkpoint is also a `PreflightError` and requires
`--overwrite`. For the unit-test contract, refusals are identified by exception
type, never by string-matching error messages.

Checkpoint written atomically via `save_checkpoint_atomic` after each
completed file.

---

## 8. Parallelism

`--original-msa` statistics are computed in the main thread (single call).
Simulated files are dispatched to a `ProcessPoolExecutor(max_workers=threads)`.
Each worker receives the file path and seq-type; returns a dict of five
statistic values or raises an exception (logged as a warning; file counted as
failed). Progress displayed via Rich `Progress` with a single task bar showing
completed / total simulated files.

---

## 9. Dry-run

`--dry-run`:
1. Detect and validate `--original-msa` (format + seq-type).
2. Scan `--simulated-dir` and count valid files.
3. Print plan: seq-type, n_taxa, n_sites, n_simulated_files, output_dir.
4. Write nothing.

---

## 10. File Map

| Action | File | Responsibility |
|---|---|---|
| Create | `phyloai/posttree/simulate_adequacy.py` | Core library: statistic functions, batch computation, checkpoint, result writing |
| Modify | `phyloai/cli/commands/posttree.py` | Replace stub `simulate_adequacy_command` with full Click command |
| Modify | `phyloai/report/collector.py` | (a) Add `"posttree.simulate.adequacy"` to `_STEP_ORDER`; (b) add `"adequacy"` to `_THIRD_LEVEL["simulate"]` set; (c) no `_FOURTH_LEVEL` entry needed |
| Modify | `phyloai/report/templates.py` | Add `generate_methods_posttree_simulate_adequacy` and register in `METHODS_GENERATORS` |
| Create | `docs/commands/posttree-simulate-adequacy.md` | English user docs |
| Create | `docs/commands/posttree-simulate-adequacy.zh.md` | Chinese user docs |
| Modify | `docs/superpowers/specs/2026-06-07-phyloai-design.md` | Update adequacy stub line; update decision table |
| Modify | `README.md` | Add adequacy workflow example |
| Modify | `README.zh.md` | Add adequacy workflow example (Chinese) |
| Modify | `skills/phyloai-workflow/SKILL.md` | Add adequacy approval card, execution/recovery guidance, result interpretation |
| Create | `tests/posttree/test_simulate_adequacy.py` | Unit tests: statistic functions, z-score/pp, CI, checkpoint, CLI |
| Modify | `tests/cli/test_posttree_simulate_alisim.py` | CLI help, dry-run, and input-error coverage for adequacy |
| Modify | `tests/report/test_collector.py` | Adequacy report step parsing and ordering |
| Modify | `tests/report/test_templates.py` | Adequacy methods-template coverage |
| Modify | `tests/mcp/test_schema_gen.py` | Click-derived adequacy MCP schema coverage |

**Note on MCP:** MCP tools are auto-generated by `cli_tools.py` walking the
Click tree (`walk_click_tree`). No manual registration in `cli_tools.py` is
required; replacing the stub Click command in `posttree.py` is sufficient for
the MCP tool to appear automatically.

---

## 11. Report Template (`generate_methods_posttree_simulate_adequacy`)

The generated methods text should describe:

1. Software and approach (pure Python, four statistics, no external tool)
2. Input: observed MSA (n_taxa × n_sites, seq-type)
3. Simulated MSAs: count (the command does not record simulation provenance,
   so the source is intentionally not stated)
4. Statistics computed: div, siteconvprob, sitecomp, comp — one sentence each
5. Null distribution: empirical mean ± population SD and 95% CI (p2.5/p97.5)
6. Comparison: z-score and posterior predictive p-value; significance threshold
   (|z| > 2 or pp < 0.05 indicate potential model inadequacy)

Example output:

> Model adequacy was assessed by comparing four summary statistics computed on
> the observed alignment (6 taxa × 235 sites, AA) against the empirical null
> distribution derived from 100 simulated replicates. The statistics evaluated
> were: mean diversity per site (PPA-DIV; Lartillot et al. 2007), mean squared
> empirical frequency (PPA-CONV), mean variance of site-specific frequencies
> (PPA-VAR), and maximum/mean squared compositional deviation across taxa
> (PPA-COMP; Blanquart & Lartillot 2006). For each statistic, the null
> distribution mean ± SD and empirical 95% confidence interval (p2.5–p97.5)
> were computed from the simulated replicates; the observed value was compared
> against this distribution using a z-score and a posterior predictive p-value
> (pp). Statistics with |z| > 2 or pp < 0.05 indicate potential model
> inadequacy. All computations were performed using phyloai posttree simulate
> adequacy (pure Python implementation).

---

## 12. Test Cases

| Test | Input | Expected |
|---|---|---|
| `test_div_known` | Tiny 4-taxon, 3-site AA matrix with known diversity | `div` = expected value ± 1e-9 |
| `test_siteconvprob_known` | Same matrix | `siteconvprob` = expected value |
| `test_sitecomp_known` | Same matrix | `sitecomp` = expected value |
| `test_comp_known` | Same matrix | `comp_max`, `comp_mean`, per-taxon `dist_j` = expected values |
| `test_gap_exclusion` | Site with all gaps → excluded from `S` | `n_informative_sites` decremented |
| `test_zscore_sign` | obs < mean_sim for `div` | z-score > 0 |
| `test_pp_div` | All 10 sim values > obs | `pp = 0.0` (P(sim ≤ obs) = 0/10) |
| `test_ci_percentile` | 100 sim values uniform [0,1] | ci_lower ≈ 0.025, ci_upper ≈ 0.975 |
| `test_duplicate_taxon_original` | Original MSA has duplicate taxon name | `ValueError` before run |
| `test_duplicate_taxon_simulated` | One simulated file has duplicate taxon name | That file skipped with warning |
| `test_taxon_mismatch` | Simulated MSA has different taxon set | File skipped with warning |
| `test_length_mismatch` | Simulated MSA has different alignment length | File skipped with warning |
| `test_all_missing_taxon_original` | One taxon all-gap in original MSA | `ValueError` before run |
| `test_all_missing_taxon_simulated` | One taxon all-gap in one simulated file | That file skipped with warning |
| `test_resume_taxon_comp` | Resume after 5/10 files, 3 taxa | `adequacy_taxon_comp.csv` identical to full run |
| `test_sd_zero` | All simulated values identical | `z_score=0.0`, `pp=null` (JSON) / empty (CSV) |
| `test_pp_direction_div` | 8/10 sim values > obs | `pp = 0.2` (P(sim ≤ obs) = 2/10) |
| `test_pp_direction_siteconvprob` | 8/10 sim values > obs | `pp = 0.8` (P(sim > obs) = 8/10) |
| `test_checkpoint_resume` | Interrupt after 5/10 files | Resume skips first 5, total = 10 |
| `test_checkpoint_fingerprint` | Replace one simulated file after partial run | Replaced file recomputed on resume |
| `test_checkpoint_original_fingerprint` | Replace the observed MSA in place after a partial run | Resume rejected with `ValueError` |
| `test_checkpoint_missing_original_fingerprint` | Resume a legacy adequacy checkpoint without the field | `PreflightError`; user must use `--overwrite` |
| `test_resume_membership_changes` | Delete one simulated file and add another after a partial run | Deleted task dropped with warning; added file recomputed; aggregation excludes the deleted task |
| `test_dry_run` | Valid inputs | No files written; prints plan |
| `test_phylip_relaxed_input` | `.ali` file (phylip-relaxed content) | Auto-detected, read correctly |
| `test_result_json_schema` | Full run on test data | All required keys present |
| `test_taxon_comp_csv` | 6-taxon matrix | `adequacy_taxon_comp.csv` has 6 rows with all columns |
| `test_table_format_tsv` | Valid inputs with `--table-format tsv` | All output tables are tab-delimited with `.tsv` suffixes and listed in result.json |
| `test_zscore_fixture_regression` | `runs/zscore/matrix.XX` plus `runs/zscore/simulated/` | Five observed/mean/SD/z/pp values match `runs/zscore/readpb/chain1.ppred` within documented rounding tolerance |

---

## 13. Interpretation Guide (for report and skill)

| z-score | pp | Interpretation |
|---|---|---|
| `|z| < 2` and `pp > 0.05` | — | Model adequately reproduces this statistic |
| `|z| > 2` or `pp < 0.05` | — | Potential model inadequacy; investigate further |
| `pp = 0` for `div` | `z >> 0` | Severe failure; all simulated diversity values exceed observed (model greatly over-predicts diversity) |
| `pp = 0` for `siteconvprob`/`sitecomp`/`comp` | `z >> 0` | Severe failure; all simulated values below observed (data far more extreme than model predicts) |
| Negative z for `div` | — | Observed diversity higher than model predicts (model under-estimates heterogeneity) |
| High z for `comp` per-taxon | pp < 0.05 | That taxon has unusually divergent composition; consider compositional heterogeneity models |

**Note on sign conventions and pp direction** (source: `PhyloProcess.cpp`):

| Statistic | Z-score formula | pp = P(…) | Low pp means |
|---|---|---|---|
| `div` | `(mean_sim − obs) / sd` | `P(sim ≤ obs)` | obs unusually low (model over-predicts diversity) |
| `siteconvprob` | `(obs − mean_sim) / sd` | `P(sim > obs)` | obs unusually high (data anomalously convergent) |
| `sitecomp` | `(obs − mean_sim) / sd` | `P(sim > obs)` | obs unusually high (data more heterogeneous across sites) |
| `comp_max` / `comp_mean` | `(obs − mean_sim) / sd` | `P(sim > obs)` | obs unusually high (data more compositionally heterogeneous) |

---

## 14. Known Limitations

- Statistics are computed from empirical MSA frequencies only, not from model
  parameters. This is correct for comparing observed vs. simulated empirical
  patterns, but differs from model-parameter-based adequacy checks.
- No missing-data re-masking step is applied to simulated MSAs before
  computing statistics (unlike PhyloBayes internal `readpb_mpi`, which masks
  simulated data to match the original gap pattern). Users running AliSim
  simulations should consider applying `phyloai posttree simulate alisim
  transfergaps` to their simulated MSAs before running adequacy checks if
  missing data is substantial.
- CI computation uses bounded inclusive quantiles; with fewer than 40 simulated
  files, the p2.5 / p97.5 estimates may be imprecise.
