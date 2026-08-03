# PhyloAI Posttree Simulate AliSim Design Specification

**Date:** 2026-08-02
**Status:** Draft — pending approval
**Parent spec:** `2026-06-07-phyloai-design.md`
**JSON standard:** `2026-06-21-phyloai-json-output-standard.md`
**References:**
- IQ-TREE3 `--alisim` documentation (AliSim: Ly-Trong et al. 2023, Molecular Biology and Evolution)
- PhyloAI reference implementation: `ref/scripts/server.R` (Shiny simulation parameter sampling)
- PhyloAI reference implementation: `ref/scripts/transfer_gaps.py` (gap transfer)

---

## 1. Purpose

`phyloai posttree simulate alisim` provides alignment simulation using IQ-TREE3 `--alisim`, organized as three subcommands:

| Subcommand | Purpose | Core tool |
|------------|---------|-----------|
| `params` | Extract simulation parameters from batch `.iqtree` files into a TSV table | Pure Python (parsing) |
| `iqtree` | Simulate MSAs using IQ-TREE3 `--alisim` (single or batch mode) | IQ-TREE3 |
| `transfergaps` | Transfer gap patterns from original MSAs to simulated gap-free MSAs | Pure Python |

**Context:** This module is part of the broader `posttree simulate` group. Future sibling subcommand groups include `adequacy` (model adequacy checks on simulated MSAs) and `phybase` (R script generation for gene tree simulation). This spec covers only the `alisim` subcommand group.

---

## 2. Design Principles

1. **Reuse existing infrastructure.** IQ-TREE path resolution via `phyloai.core.iqtree`, batch parallelism via `ProcessPoolExecutor` + Rich progress bars, checkpoint/resume via `checkpoint.json`, and file matching via the shared logical locus name policy.
2. **Separate parsing from simulation.** `params` produces a self-contained TSV table that can be inspected, edited, or used by external scripts before feeding into `iqtree`. This decoupling enables manual parameter tuning.
3. **Three sampling strategies in one command.** `iqtree` exposes `--strategy complete|mixed|pdf` as a single flag. PDF-specific options (`--noise-scale`, `--pdf-params`) are only valid when `strategy=pdf`.
4. **Comma-delimited compound values.** Multi-value fields in `params.tsv` (freq, subs_rate, rate_param) use `,` as internal delimiter, matching IQ-TREE's native `+F{}`, `+G{}` and `+R{}` format in `.iqtree` reports and the `--alisim -m` string. No conversion needed when building `--alisim` commands. (Legacy `/`-delimited tables remain readable: the command builder normalizes `/` to `,`.)
5. **Single-file gap transfer.** `transfergaps` operates on a single pair of files (one original, one simulated), matching the reference script's design. Batch gap transfer across directories is left to user scripting or shell loops.
6. **Suffix-agnostic file matching.** Tree file matching in `params` follows the project-wide logical locus name policy (Section 9.7 of the parent spec): strip 1–2 dot segments, no hardcoded suffix whitelist.

---

## 3. CLI Surface

```bash
# Extract parameters from .iqtree files
phyloai posttree simulate alisim params \
  --iqtree-dir runs/tree/ml/iqtree/logs \
  --tree-dir runs/tree/ml/iqtree/trees \
  [--output-dir runs/posttree/simulate/alisim/params] \
  [--overwrite] [--dry-run] [--quiet]

# Single-parameter simulation
phyloai posttree simulate alisim iqtree \
  --ref-tree ref.treefile \
  --model "GTR{XXX}+F{XXX}+G4{XXX}" \
  --seq-type DNA \
  --length 500 \
  [--msa-prefix sim] \
  [--num-alignments 2] \
  [--out-format fasta] \
  [--iqtree-threads 1] \
  [--seed 42] \
  [--iqtree-path /path/to/iqtree3] \
  [--tool-args "..."] \
  [--output-dir runs/posttree/simulate/alisim/iqtree] \
  [--overwrite] [--dry-run] [--quiet]

# Single-parameter simulation with partition model
phyloai posttree simulate alisim iqtree \
  --ref-tree ref.treefile \
  --model-partitions partition.PMSF.nex \
  --seq-type AA \
  [--msa-prefix sim] \
  [--num-alignments 1] \
  [--out-format fasta] \
  [--iqtree-threads 1] \
  [--seed 42] \
  [--iqtree-path /path/to/iqtree3] \
  [--tool-args "..."] \
  [--output-dir runs/posttree/simulate/alisim/iqtree] \
  [--overwrite] [--dry-run] [--quiet]

# Batch simulation from params table
phyloai posttree simulate alisim iqtree \
  --model-params params.tsv \
  --strategy complete \
  --num-simulations 100 \
  [--override length=500] \
  [--msa-prefix sim] \
  [--out-format fasta] \
  [--iqtree-threads 1] \
  [--threads 4] \
  [--seed 42] \
  [--iqtree-path /path/to/iqtree3] \
  [--tool-args "..."] \
  [--output-dir runs/posttree/simulate/alisim/iqtree] \
  [--overwrite] [--resume] [--dry-run] [--quiet]

# PDF sampling with noise control
phyloai posttree simulate alisim iqtree \
  --model-params params.tsv \
  --strategy pdf \
  --num-simulations 200 \
  --noise-scale 1.0 \
  --pdf-params length,prop_inv,rate_param \
  [--override length=500,prop_inv=0.1] \
  [--threads 8] \
  [--output-dir runs/posttree/simulate/alisim/iqtree] \
  [--overwrite] [--resume] [--dry-run] [--quiet]

# Transfer gap pattern from original to simulated MSA (single-file mode)
phyloai posttree simulate alisim transfergaps \
  --original-msa original.fa \
  --simulated-msa sim001.fa \
  [--seq-type auto] \
  [--exclude-ambiguity] \
  [--output-dir runs/posttree/simulate/alisim/transfergaps] \
  [--overwrite] [--dry-run] [--quiet]

# Transfer gap patterns from original MSA to a batch of simulated MSAs (batch mode)
phyloai posttree simulate alisim transfergaps \
  --original-msa original.fa \
  --simulated-dir MSAs/ \
  [--seq-type auto] \
  [--exclude-ambiguity] \
  [--output-dir runs/posttree/simulate/alisim/transfergaps] \
  [--overwrite] [--dry-run] [--quiet]
```

### 3.1 Command Hierarchy

```
phyloai posttree
└── simulate
    ├── alisim
    │   ├── params
    │   ├── iqtree
    │   └── transfergaps
    ├── adequacy          (future)
    └── phybase           (future)
```

### 3.2 Mutual Exclusions

- `alisim iqtree`: `--model` and `--model-partitions` are mutually exclusive.
- `alisim iqtree`: `--model-params` is mutually exclusive with `--ref-tree`, `--model`, `--model-partitions`, `--seq-type`, `--length` (single-parameter mode vs batch mode).
- `alisim iqtree`: `--strategy`, `--num-simulations`, `--noise-scale`, `--pdf-params`, `--override` require `--model-params` (batch mode only).
- `alisim iqtree`: `--num-alignments` is single-parameter mode only.
- `alisim iqtree`: `--resume` is batch mode only.
- `alisim iqtree`: `--threads` controls parallel simulation tasks (batch mode); `--iqtree-threads` controls per-task IQ-TREE `-T`.
- `alisim iqtree`: `--noise-scale` and `--pdf-params` require `--strategy pdf`.
- `alisim iqtree`: `--overwrite` and `--resume` are mutually exclusive.
- `alisim iqtree`: `--override` parameters are fixed and not sampled; they take precedence over both table values and PDF estimation.
- `alisim transfergaps`: `--simulated-msa` (single mode) and `--simulated-dir` (batch mode) are mutually exclusive; exactly one is required.

---

## 4. Parameters

### 4.1 `alisim params` Parameters

| Flag | Type | Required | Default | Description |
|------|------|----------|---------|-------------|
| `--iqtree-dir` | Path | yes | — | Directory containing `.iqtree` report files (any nesting depth; globbed as `**/*.iqtree`). |
| `--tree-dir` | Path | yes | — | Directory containing tree files. Matched to `.iqtree` files by logical locus name (suffix-agnostic per Section 9.7: strip 1–2 dot segments, no suffix whitelist). |
| `--output-dir` | Path | no | `runs/posttree/simulate/alisim/params` | Output directory for `params.tsv` and `result.json`. |
| `--overwrite` | flag | no | False | Delete and recreate output directory. |
| `--dry-run` | flag | no | False | Validate inputs, report file counts, but write nothing. |
| `--quiet` | flag | no | False | Suppress terminal output except errors. |

**Validation rules:**
- `--iqtree-dir` must exist and contain at least one `.iqtree` file.
- `--tree-dir` must exist. Tree files matched by logical locus name following the project-wide suffix-agnostic policy (Section 9.7): strip 1–2 dot segments from filenames to derive logical locus names. Ambiguity (multiple tree files with the same logical locus name) → hard error with explicit message.
- Unmatched `.iqtree` files (no corresponding tree): logged as warning, excluded from `params.tsv`, listed in `result.json` under `data.unmatched`.
- Unmatched tree files: ignored silently (trees without models are not actionable).

### 4.2 `alisim iqtree` Parameters — Single-Parameter Mode

| Flag | Type | Required | Default | Description |
|------|------|----------|---------|-------------|
| `--ref-tree` | Path | yes* | — | Reference tree (Newick). Maps to IQ-TREE `-t`. |
| `--model` | str | yes* | — | IQ-TREE model string (e.g. `"GTR{XXX}+F{XXX}+G4{XXX}"`). Maps to IQ-TREE `-m`. Mutually exclusive with `--model-partitions`. |
| `--model-partitions` | Path | yes* | — | NEXUS partition model file (e.g. `partition.PMSF.nex`). Maps to IQ-TREE `-p` (edge-proportional partition model: shared topology, per-partition rate multipliers). Mutually exclusive with `--model`. When used, `--length` is not required (inferred from partition definitions). |
| `--seq-type` | choice | yes* | — | `AA` or `DNA`. Maps to IQ-TREE `--seqtype`. |
| `--length` | int | yes* | — | Alignment length. Maps to IQ-TREE `--length`. Not required when `--model-partitions` is used. |
| `--msa-prefix` | str | no | `sim` | Output MSA file prefix. Maps to IQ-TREE `--alisim`. |
| `--num-alignments` | int | no | `1` | Number of MSAs to simulate in one IQ-TREE call. Maps to IQ-TREE `--num-alignments`. Single-parameter mode only. |
| `--out-format` | choice | no | `fasta` | `fasta` or `phy`. Maps to IQ-TREE `--out-format`. |
| `--iqtree-threads` | int | no | `1` | Threads per IQ-TREE invocation. Maps to IQ-TREE `-T`. |
| `--seed` | int | no | random | Random seed. Maps to IQ-TREE `--seed`. Default: `random.randint(1, 2**31 - 1)`. |
| `--iqtree-path` | str | no | — | Custom IQ-TREE3 executable path. |
| `--tool-args` | str | no | — | Extra IQ-TREE flags. Blocked (PhyloAI-managed I/O only): `--alisim`, `-t`, `--prefix`, `--out-format`, `-af`. All other flags (e.g. `--seqtype`, `--length`, `--num-alignments`, `-T`, `--seed`, `-m`, `-p`) may be re-specified to override PhyloAI's defaults; PhyloAI then suppresses its own copy of the overridden flag so the final IQ-TREE command contains each flag exactly once with the tool-args value. |
| `--output-dir` | Path | no | `runs/posttree/simulate/alisim/iqtree` | Output directory. |
| `--overwrite` | flag | no | False | Delete and recreate output directory. |
| `--dry-run` | flag | no | False | Print IQ-TREE command without executing. |
| `--quiet` | flag | no | False | Suppress terminal output except errors. |

*Required in single-parameter mode (when `--model-params` is not used). `--model` or `--model-partitions` (exactly one), plus `--ref-tree` and `--seq-type`, are required. `--length` is required unless `--model-partitions` is used.

### 4.3 `alisim iqtree` Parameters — Batch Mode

| Flag | Type | Required | Default | Description |
|------|------|----------|---------|-------------|
| `--model-params` | Path | yes | — | TSV table from `alisim params` (or manually constructed). Activates batch mode. |
| `--strategy` | choice | no | `complete` | Sampling strategy: `complete`, `mixed`, `pdf`. |
| `--num-simulations` | int | yes | — | Total number of MSAs to simulate. Each simulation is one IQ-TREE `--alisim` call producing 1 MSA. |
| `--override` | str | no | — | Comma-separated `key=value` pairs to fix specific parameters across all simulations (e.g. `length=500,prop_inv=0.1`). Overridden parameters are never sampled, for **all** strategies (complete, mixed, pdf). Valid keys: `length`, `prop_inv`. Coupled group columns (`seqtype`, `subs_model`, `subs_rate`, `freq`, `rate_heterogeneity`, `rate_categories`, `rate_param`) cannot be individually overridden because partial override would break model consistency. To fix model parameters, pre-filter `params.tsv` manually. |
| `--noise-scale` | float | no | `1.0` | Histogram resampling noise amplitude: 0 = no noise (bin centers), 1 = full within-bin uniform jitter. Range [0.0, 1.0]. Only valid with `--strategy pdf`. Reported as `null` in result.json for non-pdf strategies. |
| `--pdf-params` | str | no | `length,prop_inv,rate_param` | Comma-separated parameter columns to sample via histogram-based density resampling. Only valid with `--strategy pdf`. Valid values: `length`, `prop_inv`, `rate_param`. `rate_param` is only PDF-sampled when `rate_heterogeneity` is `G` (Gamma); FreeRate (`R`) parameters are always empirically sampled due to their multi-valued coupled structure. Reported as `null` in result.json for non-pdf strategies. |
| `--msa-prefix` | str | no | `sim` | Output MSA file prefix. |
| `--out-format` | choice | no | `fasta` | `fasta` or `phy`. |
| `--iqtree-threads` | int | no | `1` | Threads per IQ-TREE invocation. |
| `--threads` | int | no | `4` | Parallel simulation tasks (ProcessPoolExecutor). |
| `--seed` | int | no | random | Master random seed. Each simulation gets an independent random seed drawn from a master-seeded generator (reproducible: same master seed yields the same per-simulation seed sequence). |
| `--iqtree-path` | str | no | — | Custom IQ-TREE3 executable path. |
| `--tool-args` | str | no | — | Extra IQ-TREE flags (same blocked set as single mode). |
| `--output-dir` | Path | no | `runs/posttree/simulate/alisim/iqtree` | Output directory. |
| `--overwrite` | flag | no | False | |
| `--resume` | flag | no | False | Resume from checkpoint. |
| `--dry-run` | flag | no | False | Show sampling plan and example commands without executing. |
| `--quiet` | flag | no | False | |

### 4.4 `alisim transfergaps` Parameters

| Flag | Type | Required | Default | Description |
|------|------|----------|---------|-------------|
| `--original-msa` | Path | yes | — | Single original (gapped) MSA file. |
| `--simulated-msa` | Path | yes* | — | Single simulated (gap-free) MSA file from `alisim iqtree`. Mutually exclusive with `--simulated-dir`. |
| `--simulated-dir` | Path | yes* | — | Directory of simulated (gap-free) MSA files from `alisim iqtree` (alignment extensions only, e.g. `.fa`, `.fasta`, `.phy`, `.nex`). One transferred file is written per input as `<stem>.gaps.fa` in the output directory. Mutually exclusive with `--simulated-msa`. |
| `--seq-type` | choice | no | `auto` | `AA`, `NT`, or `auto`. Determines the valid character set for identifying non-standard positions. AA valid: `ACDEFGHIKLMNPQRSTVWY`; NT valid: `ACGT`. |
| `--exclude-ambiguity` | flag | no | False | When set, only real gap characters (`-`, `.`) are transferred; ambiguity codes (X, N, B, Z, ?, etc.) are left as simulated characters. Default behavior (without this flag): all non-standard characters (gaps + ambiguity codes) are replaced with `-`. |
| `--output-dir` | Path | no | `runs/posttree/simulate/alisim/transfergaps` | Output directory for result.json and the transferred MSA file(s). |
| `--overwrite` | flag | no | False | Delete and recreate output directory. |
| `--dry-run` | flag | no | False | Validate inputs (file existence, taxon name matching, length compatibility) and report what would be done without writing any output files or creating the output directory. |
| `--quiet` | flag | no | False | Suppress terminal output except errors. |

*Exactly one of `--simulated-msa` or `--simulated-dir` is required.

---

## 5. Computation

### 5.1 `alisim params` — Parameter Extraction

**Parsing algorithm per `.iqtree` file:**

1. **Seq type:** Locate the line starting with `"Input data:"`. Extract `amino-acid` → `AA`, `nucleotide` → `DNA` via regex.

2. **AliSim command line:** Locate the line starting with `"To simulate an alignment"`. The *next line* contains the full `--alisim` command string. Extract:
   - `--length <N>` → `length` (integer)
   - `-m "<model_string>"` → full model string for parsing

3. **Model string parsing** (regex on the `-m` value):
   - **subs_model:** The model name before any `{` or `+` (e.g. `GTR`, `LG`, `WAG`).
   - **subs_rate:** For DNA models with rate matrix (e.g. `GTR{a,b,c,d,e}`), the content inside the first `{}` after the model name. Empty for AA models (which have no explicit rate matrix). `,`-delimited.
   - **freq:** Content of `+F{...}` if present. `,`-delimited. If `+F` is present but `{}` is empty, this indicates empirical frequencies were used but values must be extracted from pi() lines (see step 4).
   - **prop_inv:** Content of `+I{...}` if present. Single float value.
   - **rate_heterogeneity:** `G` if `+G<N>` present, `R` if `+R<N>` present. Empty if neither.
   - **rate_categories:** The integer `N` from `+G<N>` or `+R<N>`.
   - **rate_param:** Content of `+G<N>{...}` (single alpha value for Gamma) or `+R<N>{...}` (alternating proportion/rate pairs for FreeRate). `,`-delimited.

4. **AA frequency extraction (when +F in model but no explicit freq values in AliSim line):**
   Parse `pi(X) = value` lines from the `.iqtree` file. Collect all 20 amino acid frequencies in canonical ARNDCQEGHILKMFPSTWYV order. Join with `,` delimiter. This handles AA models like `LG+F+G4` where IQ-TREE's AliSim command line uses `+F{values}` but older versions or certain models may require extraction from the report body.

   **Note:** In practice, IQ-TREE3 always includes explicit freq values in the AliSim command line (e.g. `+F{0.0819,0.0436,...}`). The pi() fallback is a safety net for edge cases. The freq values from the AliSim command line take priority; pi() lines are only used when `+F` appears in the model without explicit `{values}`.

5. **Tree matching:** For each `.iqtree` file, derive the logical locus name (remove `.iqtree` suffix). Search `--tree-dir` using the project-wide suffix-agnostic policy: for each file in the tree directory, strip 1–2 dot segments to derive its logical locus name; match against the `.iqtree` locus name. If multiple tree files share the same logical locus name → hard error. No match → warning + skip this locus. Store the absolute path of the matched tree file in `tree_path`.

**Output:** `params.tsv` (TSV), one row per successfully parsed locus:

```
id	seqtype	length	subs_model	subs_rate	freq	prop_inv	rate_heterogeneity	rate_categories	rate_param	tree_path
EOG090X002Z	AA	2082	LG				G	4	0.609836	/abs/path/to/EOG090X002Z.treefile
EOG090X005G	AA	1879	LG		0.0819,0.0436,...	0.1403	G	4	1.1873	/abs/path/to/EOG090X005G.treefile
EOG090X002Z	DNA	6246	GTR	2.804,3.677,1.359,1.875,6.215	0.3114,0.1972,0.2256,0.2657	0.2141	G	4	1.39694	/abs/path/to/EOG090X002Z.treefile
EOG090X0064	AA	1640	LG					R	2	0.8637,0.351,0.1363,5.1118	/abs/path/to/EOG090X0064.treefile
```

Empty fields indicate the parameter is absent from the model (e.g. no `+I`, no `+F`, no subs_rate for AA models).

### 5.2 `alisim iqtree` — Single-Parameter Mode

**Note:** The `length` value extracted from `params.tsv` is the **alignment length** (total columns including gap positions). AliSim's `--length` parameter produces a simulated MSA of exactly this many columns (all ungapped). This ensures simulated MSAs have the same column count as the original MSAs, which is a prerequisite for `transfergaps` to work correctly.

**IQ-TREE command assembly:**

```
iqtree3 --alisim <msa_prefix> --seqtype <seqtype> -t <ref_tree> -m "<model>" --length <length> --out-format <out_format> --num-alignments <num_alignments> -T <iqtree_threads> --seed <seed> [tool_args]
```

Flags re-specified in `--tool-args` replace (rather than duplicate) PhyloAI's managed flag, so each flag appears once in the final command. The full re-executable command recorded in `result.json` uses shell quoting (`shlex.join`) so multi-token `--tool-args` values are unambiguous.

When `--model-partitions` is used instead of `--model`:
```
iqtree3 --alisim <msa_prefix> --seqtype <seqtype> -t <ref_tree> -p <model_partitions> --out-format <out_format> --num-alignments <num_alignments> -T <iqtree_threads> --seed <seed> [tool_args]
```

With `--model-partitions`, `--length` is omitted (partition file defines site ranges and per-partition models). `-p` (edge-proportional) shares the tree topology with per-partition rate multipliers, which is correct for PMSF-style per-site model definitions.

**Output file handling:**
- IQ-TREE `--alisim` produces files like `sim_1.fa`, `sim_2.fa`, ... (when `--num-alignments > 1`) or `sim.fa` (when `--num-alignments 1`).
- AliSim does **not** write a `{prefix}.iqtree` report or a `{prefix}.log` in the working directory; its console/screen log is written to `<ref_tree>.log` next to the tree. PhyloAI therefore captures IQ-TREE's stdout/stderr and writes it to `<output-dir>/logs/<msa_prefix>.log`, then removes the stray `<ref_tree>.log` (only if it did not pre-exist before the run).
- Move generated MSA files to `<output-dir>/MSAs/` (all files matching `<msa_prefix>*<ext>`, so a `--tool-args --num-alignments` override is honored).
- **Validation (per Section 9.10):** Each generated MSA file is validated before counting as successful: must be non-empty, parsable as FASTA/PHYLIP, contain ≥1 sequence, and all sequences must be equal length. Failed validation → recorded in result.json as failed.
- **FASTA line wrapping (per Section 9.11):** If `--out-format fasta`, generated FASTA files are reformatted to 60-character line wrapping before saving to `MSAs/`.

### 5.3 `alisim iqtree` — Batch Mode

Each of `--num-simulations` tasks:
1. Apply `--override` values (if any) to fix specified parameters. Overrides are applied after sampling for **all** strategies (complete, mixed, pdf).
2. Sample remaining parameters from `params.tsv` according to `--strategy`.
3. Build one IQ-TREE `--alisim` command producing 1 MSA.
4. Execute with an independent random seed per simulation (drawn from a master-seeded generator; reproducible given `--seed`).
5. Name output: `<msa_prefix><NNN>.<ext>` (zero-padded index, e.g. `sim001.fa`).
6. Validate output MSA (per Section 9.10): non-empty, parsable, sequences equal length. Failed → recorded as failed in checkpoint and result.json.

**Sampling strategies:**

#### Strategy 1: `complete` — Complete empirical match

All parameters come from the same row of `params.tsv`. A row is randomly selected (uniform, with replacement). The entire row's parameters are used as-is.

Rows with empty `tree_path` are excluded (no tree → cannot simulate).

#### Strategy 2: `mixed` — Mixed empirical sampling

Each parameter is independently sampled from its column in `params.tsv`, with these coupling constraints:

- **Model core group** (`seqtype` + `subs_model` + `subs_rate` + `freq`): sampled together from the same row. These four columns are structurally interdependent: the substitution model defines the rate matrix structure (GTR has 5 rates, AA models have none), and frequency dimensionality depends on sequence type (4 for DNA, 20 for AA). Mixing across sequence types would produce invalid models. Empty values within the group are valid states (e.g. `subs_rate` is empty for AA models; `freq` is empty for models without `+F`).
- **Rate heterogeneity group** (`rate_heterogeneity` + `rate_categories` + `rate_param`): sampled together from the same row. Gamma alpha vs FreeRate proportion/rate pairs are structurally different and cannot be mixed. All three columns may be empty (model without `+G`/`+R`) — this is a valid sampled state meaning no rate heterogeneity.
- **Independent parameters:** `length`, `prop_inv`, `tree_path` — each independently sampled from their respective columns.

When a coupled group is sampled, all columns within the group come from the same row. Empty values within a group are valid (they represent legitimate model configurations such as AA models without explicit rates, or models without `+F`/`+G`/`+R`). The only hard requirement is that `tree_path` must be non-empty for the simulation to execute.

**`prop_inv` two-step sampling:** `prop_inv` is sampled in two steps to preserve the empirical presence/absence ratio:
1. **Presence decision:** Sample whether `+I` is present by drawing from the empirical presence/absence ratio in the `prop_inv` column (i.e., fraction of rows with non-empty values).
2. **Value (if present):** If the presence decision is "yes", sample the actual `prop_inv` value from the non-empty values in the column (empirically in `mixed`, or via density resampling in `pdf` mode).
3. **If absent:** `prop_inv` is empty — the model does not include `+I`.

This parallels the rate heterogeneity group: the coupled group first determines whether `+G`/`+R`/none is present (and which type), then the `rate_param` value is sampled only from rows matching that decision.

**Tree selection in mixed mode:** `tree_path` is independently sampled from rows with non-empty `tree_path`.

**Cross-group consistency:** After sampling all groups independently, no additional cross-group validation is needed because `seqtype` is locked to the model core group, ensuring the final assembled model string is always self-consistent.

#### Strategy 3: `pdf` — Histogram-based density resampling

Built on top of `mixed` sampling. For parameters listed in `--pdf-params`, values are drawn from histogram-based density resampling (Freedman-Diaconis binning) instead of directly from the empirical column:

1. Compute histogram with Freedman-Diaconis bin edges on the non-empty values.
2. Sample a bin weighted by density (probability proportional to bin count).
3. Within the selected bin, draw a uniform value in `[bin_left, bin_right)` scaled by `--noise-scale`:
   - `noise-scale=0`: use bin midpoint (no jitter).
   - `noise-scale=1`: uniform jitter across full bin width (default).
   - Values between 0 and 1: proportional jitter around midpoint.
4. **Post-processing:** `length` is rounded to the nearest positive integer (min 1). `prop_inv` is clamped to [0.0, 1.0) — values resampled near boundaries are kept within valid range.

**Note:** This is histogram-based resampling with within-bin jitter, not kernel density estimation (KDE). The approach smooths the empirical distribution without requiring bandwidth selection.

**Eligibility per parameter:**
- `length`: always eligible.
- `prop_inv`: eligible for density resampling of the **value** (step 2 above). The presence/absence decision (step 1) is always empirically sampled regardless of `--pdf-params`. Only non-empty `prop_inv` values are used for the density estimation.
- `rate_param`: eligible only when the rate heterogeneity group sampled `rate_heterogeneity=G` (Gamma, single alpha value). When `rate_heterogeneity=R` (FreeRate), the multi-valued coupled structure makes independent density estimation inappropriate; the entire rate group is empirically sampled. When `rate_heterogeneity` is empty (no rate heterogeneity), `rate_param` is also empty — no density resampling occurs.

Parameters not in `--pdf-params` use `mixed` sampling logic.

**`--override` interaction:** Parameters listed in `--override` are never sampled (neither empirically nor via density resampling). They use the fixed value from `--override` regardless of `--strategy` or `--pdf-params`.

**Diagnostic plots (PDF mode only):**
For each PDF-sampled parameter, generate a density comparison plot saved to `<output-dir>/plots/<param>_density.pdf`:
- Empirical distribution: density curve in `#2E86AB` (matching `ref/scripts/server.R`).
- Simulated distribution: density curve in `#A23B72`.
- Both drawn as Gaussian-KDE density curves only (no histogram, no fill), matching the R script's `geom_density`.
- Plots are generated **only** for `--strategy pdf`; `complete` and `mixed` runs produce no `plots/` directory (the empirical-vs-sampled density comparison is only meaningful for density-resampled parameters).

### 5.4 Batch Output Naming and Tracking

Each batch simulation produces one MSA named `<msa_prefix><NNN>.<ext>` (e.g. `sim001.fa`). The `params_sampled.tsv` file records the full provenance for each simulation:

| Column | Description |
|--------|-------------|
| `simulation_id` | Output filename stem (e.g. `sim001`) |
| `source_id` | For `complete` strategy only: the `id` of the source row. Omitted entirely for `mixed`/`pdf` (rows are assembled from multiple sources, so the column would be meaningless). |
| `seqtype` | Actual seqtype used |
| `length` | Actual length used |
| `subs_model` | Actual model used |
| `subs_rate` | Actual rates used |
| `freq` | Actual frequencies used |
| `prop_inv` | Actual prop_inv used (empty if none) |
| `rate_heterogeneity` | Actual rate het type |
| `rate_categories` | Actual categories |
| `rate_param` | Actual rate parameter |
| `tree_path` | Actual tree used |
| `seed` | Per-task seed |

This enables full reproducibility and post-hoc analysis of which parameters were used for each simulation.

### 5.5 Model String Reconstruction

Given sampled (or overridden) parameters, reconstruct the IQ-TREE `-m` string:

```
<subs_model>{<subs_rate>}+F{<freq>}+I{<prop_inv>}+<rate_het><rate_cats>{<rate_param>}
```

Rules:
- `subs_rate`: include `{values}` only if non-empty (DNA models).
- `+F{values}`: include only if `freq` is non-empty.
- `+I{value}`: include only if `prop_inv` is non-empty.
- `+G<N>{value}` or `+R<N>{values}`: include only if `rate_heterogeneity` is non-empty.
- Values use `,` as delimiter within `{}` (IQ-TREE format). The table stores values already `,`-delimited; `/`-delimited legacy values are normalized for compatibility.

### 5.6 `alisim transfergaps` — Gap Pattern Transfer

`transfergaps` supports two mutually exclusive modes selected by the required input flag:

- **Single-file mode** (`--simulated-msa`): one original MSA + one simulated MSA → one output MSA, mirroring the reference script `ref/scripts/transfer_gaps.py`.
- **Batch mode** (`--simulated-dir`): one original MSA + a directory of simulated MSAs → one output MSA per simulated input. Each simulated file is validated and transferred independently; the first failing file raises a hard error naming that file.

In batch mode, alignment files are discovered by extension using `COMMON_ALIGNMENT_EXTENSIONS` from `phyloai.core.schema` (same set used by `pretree stats`), so non-alignment files are skipped. Output naming per simulated input is `<simulated_stem>.gaps.fa` (e.g. `sim001.fa` → `sim001.gaps.fa`, `sim002.phy` → `sim002.gaps.fa`).

**Input formats:** All four formats handled by the shared `FormatConverter` (the same reader used by `pretree convert`) are accepted — FASTA, PHYLIP, NEXUS, and PHYLIP-PAML — for both the original and each simulated MSA. Format is auto-detected per file (extension-first, content fallback).

**Output format:** Always FASTA with 60-character line wrapping, regardless of input format. Output filenames always end in `.gaps.fa`.

**Prerequisite:** Each simulated MSA must have the same alignment length as the original MSA. `transfergaps` validates this at runtime (hard error if mismatched). In a typical workflow, this is satisfied because `alisim params` extracts `--length` from the `.iqtree` file's AliSim command — that value equals the original alignment column count — and the subsequent `alisim iqtree` simulation produces MSAs of exactly that length.

**Gap transfer algorithm** (per simulated MSA):
1. Read original MSA and simulated MSA.
2. **Taxon matching by name:** Build a mapping from taxon names in the original to taxon names in the simulated MSA. Taxon sets must be identical (same names, order-independent). Taxon names from the original are used in the output; output order follows the original MSA order.
3. For each matched taxon pair, identify positions to mask in the original sequence:
   - **Default behavior (without `--exclude-ambiguity`):** All positions where the character is NOT in the valid alphabet (`ACDEFGHIKLMNPQRSTVWY` for AA, `ACGT` for NT) are masked. This includes gap characters (`-`, `.`) AND ambiguity codes (`X`, `N`, `B`, `Z`, `?`, etc.).
   - **With `--exclude-ambiguity`:** Only real gap characters (`-`, `.`) are masked. Ambiguity codes are left untouched (simulated character is preserved at those positions).
4. **Replace** (not insert) characters at masked positions in the simulated sequence with `-`.
5. Validate: simulated sequence length equals original sequence length (should always hold since AliSim `--length` uses the full alignment length).
6. Write the output MSA:
   - Single mode: `<output-dir>/<original_filename_stem>.gaps.fa`.
   - Batch mode: `<output-dir>/<simulated_filename_stem>.gaps.fa` for each simulated input.
   - Both always written as FASTA with 60-character line wrapping (per Section 9.11), regardless of input format.

**Error handling:**
- Taxon name mismatch (original and simulated have different taxon sets) → hard error (exit code 1) with list of mismatched names, naming the failing simulated file in batch mode.
- Length mismatch between original and simulated sequences → hard error with message explaining the expected length relationship.
- Empty or unparsable input files → hard error (exit code 1).
- `--simulated-dir` that does not exist, or contains no alignment files → hard error.
- Both `--simulated-msa` and `--simulated-dir` provided (or neither) → hard error requiring exactly one.

---

## 6. Output Structure

### 6.1 `alisim params`

```
runs/posttree/simulate/alisim/params/
├── result.json
└── params.tsv
```

### 6.2 `alisim iqtree` (single-parameter mode)

```
runs/posttree/simulate/alisim/iqtree/
├── result.json
├── MSAs/
│   ├── sim.fa              # (--num-alignments 1)
│   └── sim_1.fa ... sim_N.fa  # (--num-alignments > 1)
└── logs/
    └── sim.log             # captured IQ-TREE console output
```

Note: AliSim does not produce a `.iqtree` report; `logs/` contains only the captured console log. The output MSA count follows what IQ-TREE actually generated (a `--tool-args --num-alignments` override is honored).

### 6.3 `alisim iqtree` (batch mode)

```
runs/posttree/simulate/alisim/iqtree/
├── result.json
├── checkpoint.json
├── params_sampled.tsv          # actual parameters used per simulation (see §5.4)
├── MSAs/
│   ├── sim001.fa
│   ├── sim002.fa
│   └── ...
├── logs/
│   ├── sim001.iqtree
│   ├── sim001.log
│   └── ...
└── plots/                      # PDF mode only
    ├── length_density.pdf
    ├── prop_inv_density.pdf
    └── rate_param_density.pdf
```

`params_sampled.tsv`: Full provenance per simulation (see Section 5.4 for column definitions).

### 6.4 `alisim transfergaps`

Single mode (`--simulated-msa`):
```
runs/posttree/simulate/alisim/transfergaps/
├── result.json
├── <original_stem>.gaps.fa
```

Batch mode (`--simulated-dir`):
```
runs/posttree/simulate/alisim/transfergaps/
├── result.json
├── sim001.gaps.fa
├── sim002.gaps.fa
└── ...
```

Output MSA filenames: both modes use `<stem>.gaps.fa` — single mode derives from the original MSA filename stem, batch mode derives from each simulated MSA filename stem. Output is always FASTA regardless of input format.

---

## 7. result.json Schema

### 7.1 `alisim params`

```json
{
  "status": "success|error",
  "command": "phyloai posttree simulate alisim params --iqtree-dir ... --tree-dir ... -o ... [--overwrite] [--dry-run] [--quiet]",
  "wall_time": 1.2,
  "tool_versions": {},
  "params": {
    "iqtree_dir": "/abs/path/to/iqtree_dir",
    "tree_dir": "/abs/path/to/tree_dir",
    "output_dir": "/abs/path/to/output",
    "overwrite": false,
    "dry_run": false,
    "quiet": false
  },
  "key_results": {
    "n_loci_parsed": 780,
    "n_loci_matched": 780,
    "n_loci_unmatched": 20,
    "seq_types": {"AA": 780}
  },
  "error": null,
  "error_category": null,
  "data": {
    "output_files": {
      "params_tsv": {"path": "/abs/path/to/output/params.tsv", "description": "Simulation parameters extracted from IQ-TREE reports (TSV)"}
    },
    "unmatched": [
      {"id": "EOG090X0ZZZ", "reason": "no matching tree file found"}
    ]
  }
}
```

### 7.2 `alisim iqtree` (single-parameter mode)

```json
{
  "status": "success|error",
  "command": "phyloai posttree simulate alisim iqtree --ref-tree ... --model ... --seq-type DNA ...",
  "wall_time": 3.5,
  "tool_versions": {"iqtree3": "3.1.2"},
  "params": {
    "ref_tree": "/abs/path/to/ref.treefile",
    "model": "GTR{...}+F{...}+G4{...}",
    "model_partitions": null,
    "seq_type": "DNA",
    "length": 500,
    "msa_prefix": "sim",
    "num_alignments": 2,
    "out_format": "fasta",
    "iqtree_threads": 1,
    "seed": 42,
    "output_dir": "/abs/path/to/output",
    "iqtree_path": null,
    "tool_args": null,
    "overwrite": false,
    "dry_run": false,
    "quiet": false
  },
  "key_results": {
    "n_msas_generated": 2,
    "seq_type": "DNA",
    "length": 500,
    "model": "GTR{...}+F{...}+G4{...}"
  },
  "error": null,
  "error_category": null,
  "data": {
    "cmd": ["iqtree3", "--alisim", "sim", ...],
    "tool_stderr": "...",
    "output_files": {
      "msas": {"path": "/abs/path/to/output/MSAs/", "description": "Simulated MSA files"},
      "iqtree_log": {"path": "/abs/path/to/output/logs/sim.log", "description": "IQ-TREE console log"}
    }
  }
}
```

### 7.3 `alisim iqtree` (batch mode)

```json
{
  "status": "success|error",
  "command": "phyloai posttree simulate alisim iqtree --model-params ... --strategy pdf --num-simulations ... [--override ...] [--noise-scale ... --pdf-params ...] --msa-prefix ... --out-format ... --iqtree-threads ... -t ... --seed ... [--tool-args ...] -o ... [--overwrite] [--resume] [--dry-run] [--quiet]",
  "wall_time": 120.5,
  "tool_versions": {"iqtree3": "3.1.2"},
  "params": {
    "model_params": "/abs/path/to/params.tsv",
    "strategy": "pdf",
    "num_simulations": 200,
    "override": null,
    "noise_scale": 1.0,
    "pdf_params": "length,prop_inv,rate_param",
    "msa_prefix": "sim",
    "out_format": "fasta",
    "iqtree_threads": 1,
    "threads": 8,
    "seed": 42,
    "output_dir": "/abs/path/to/output",
    "iqtree_path": null,
    "tool_args": null,
    "overwrite": false,
    "resume": false,
    "dry_run": false,
    "quiet": false
  },
  "key_results": {
    "n_simulations_requested": 200,
    "n_simulations_completed": 200,
    "n_simulations_failed": 0,
    "strategy": "pdf",
    "source_loci": 780
  },
  "error": null,
  "error_category": null,
  "data": {
    "output_files": {
      "msas_dir": {"path": "/abs/path/to/output/MSAs/", "description": "Simulated MSA files (one per simulation)"},
      "params_sampled": {"path": "/abs/path/to/output/params_sampled.tsv", "description": "Actual parameters used per simulation (TSV)"},
      "length_density": {"path": "/abs/path/to/output/plots/length_density.pdf", "description": "Empirical vs simulated density: alignment length"},
      "prop_inv_density": {"path": "/abs/path/to/output/plots/prop_inv_density.pdf", "description": "Empirical vs simulated density: proportion of invariable sites"},
      "rate_param_density": {"path": "/abs/path/to/output/plots/rate_param_density.pdf", "description": "Empirical vs simulated density: Gamma alpha parameter"}
    },
    "files": [
      {"simulation_id": "sim001", "status": "success", "wall_time": 0.6, "cmd": ["iqtree3", "--alisim", "sim001", ...], "log_file": "logs/sim001.log", "output_file": "MSAs/sim001.fa"},
      {"simulation_id": "sim002", "status": "success", "wall_time": 0.5, "cmd": ["iqtree3", "--alisim", "sim002", ...], "log_file": "logs/sim002.log", "output_file": "MSAs/sim002.fa"}
    ]
  }
}
```

For `complete`/`mixed` strategies: `params.noise_scale` and `params.pdf_params` are `null`, `data.output_files` contains no density plots, and `--noise-scale`/`--pdf-params` are absent from `command`.

### 7.4 `alisim transfergaps`

Single mode (`--simulated-msa`):
```json
{
  "status": "success|error",
  "command": "phyloai posttree simulate alisim transfergaps --original-msa ... --simulated-msa ... --seq-type auto -o ... [--exclude-ambiguity] [--overwrite] [--dry-run] [--quiet]",
  "wall_time": 0.3,
  "tool_versions": {},
  "params": {
    "original_msa": "/abs/path/to/original.fa",
    "simulated_msa": "/abs/path/to/sim001.fa",
    "simulated_dir": null,
    "seq_type": "auto",
    "exclude_ambiguity": false,
    "output_dir": "/abs/path/to/output",
    "overwrite": false,
    "quiet": false
  },
  "key_results": {
    "n_msas": 1,
    "n_sequences": 6,
    "alignment_length": 2082,
    "mean_positions_masked_per_taxon": 142,
    "detected_seq_type": "AA"
  },
  "error": null,
  "error_category": null,
  "data": {
    "output_files": {
      "transferred_msa": {"path": "/abs/path/to/output/original.gaps.fa", "description": "Simulated MSA with gap pattern transferred from original"}
    }
  }
}
```

Batch mode (`--simulated-dir`): same shape, except `params.simulated_msa` is `null`, `params.simulated_dir` is the directory path, `key_results.n_msas` is the number of simulated MSAs processed, and `data.output_files.transferred_msas` maps each simulated filename stem to its transferred file:
```json
{
  "key_results": {
    "n_msas": 12,
    "n_sequences": 6,
    "alignment_length": 2082,
    "mean_positions_masked_per_taxon": 142,
    "detected_seq_type": "AA"
  },
  "data": {
    "output_files": {
      "transferred_msas": {
        "sim001": {"path": "/abs/path/to/output/sim001.gaps.fa", "description": "Simulated MSA with gap pattern transferred from original"},
        "sim002": {"path": "/abs/path/to/output/sim002.gaps.fa", "description": "Simulated MSA with gap pattern transferred from original"}
      }
    }
  }
}
```

---

## 8. Terminal Output

### 8.1 `alisim params`

```
Extracting simulation parameters...
  .iqtree files found: 800
  Tree files found:    800
  Matched pairs:       780
  Unmatched (no tree): 20

Sequence types: AA=780

Params written to runs/posttree/simulate/alisim/params/params.tsv
Result written to runs/posttree/simulate/alisim/params/result.json
```

### 8.2 `alisim iqtree` (single)

```
IQ-TREE AliSim command:
  iqtree3 --alisim sim --seqtype DNA -t ref.treefile -m "GTR{...}+F{...}+G4{...}" --length 500 --out-format fasta --num-alignments 2 -T 1 --seed 42

[IQ-TREE stdout]

Simulated 2 MSAs → runs/posttree/simulate/alisim/iqtree/MSAs/
Result written to runs/posttree/simulate/alisim/iqtree/result.json
```

### 8.3 `alisim iqtree` (batch)

```
Batch simulation: 200 tasks (strategy: pdf, noise_scale: 1.0)
Source: 780 loci in params.tsv

Simulating... ━━━━━━━━━━━━━━━━━━━━ 100% 200/200

Completed: 200 | Failed: 0 | Skipped: 0

Density plots → runs/posttree/simulate/alisim/iqtree/plots/
MSAs → runs/posttree/simulate/alisim/iqtree/MSAs/
Sampled params → runs/posttree/simulate/alisim/iqtree/params_sampled.tsv
Result written to runs/posttree/simulate/alisim/iqtree/result.json
```

### 8.4 `alisim transfergaps`

Single mode (`--simulated-msa`):
```
Gap transfer: original.fa → sim001.fa
  Sequences: 6
  Non-standard positions masked: 142 (mean per taxon)
  Output: runs/posttree/simulate/alisim/transfergaps/original.gaps.fa
Result written to runs/posttree/simulate/alisim/transfergaps/result.json
```

Batch mode (`--simulated-dir`):
```
Gap transfer: original.fa → sims
  MSAs processed: 12
  Sequences: 6
  Non-standard positions masked: 142 (mean per taxon)
  Outputs (12 files): sim001, sim002, ...
Result written to runs/posttree/simulate/alisim/transfergaps/result.json
```

---

## 9. Report Template

### 9.1 `posttree.simulate.alisim.params`

> Simulation parameters were extracted from {n_loci_parsed} IQ-TREE report files. {n_loci_matched} loci were successfully matched with corresponding tree files. {unmatched_sentence}Extracted parameters include substitution model, state frequencies, proportion of invariable sites, rate heterogeneity model, and alignment length for each locus.

Where `{unmatched_sentence}` (only if unmatched > 0):
> {n_loci_unmatched} loci could not be matched to tree files and were excluded.

### 9.2 `posttree.simulate.alisim.iqtree` (single)

> Sequence alignment was simulated using IQ-TREE3 v{version} AliSim (Ly-Trong et al. 2023). The simulation used a {seqtype} {model} model with {length} sites. {n_msas_generated} replicate alignment(s) were generated from the reference tree with random seed {seed}.{partition_sentence}

Where `{partition_sentence}` (only if `--model-partitions` used):
> A per-site partition model was used with edge-proportional branch lengths (-p).

### 9.3 `posttree.simulate.alisim.iqtree` (batch)

> Sequence alignments were simulated using IQ-TREE3 v{version} AliSim (Ly-Trong et al. 2023). Simulation parameters were drawn from an empirical distribution of {source_loci} gene models using {strategy_sentence}. {override_sentence}A total of {n_simulations_completed} alignments were generated ({n_simulations_failed} failed).{tree_sentence}

Where `{strategy_sentence}` depends on the strategy:
- `complete`:
  > the complete sampling strategy, in which each simulated alignment replicates the full parameter set of a single source gene model (substitution model, rate heterogeneity, alignment length, invariant-site proportion, and reference tree all taken together from one empirical row of the source table)
- `mixed`:
  > the mixed sampling strategy, in which the model core, rate heterogeneity group, alignment length, invariant-site proportion, and reference tree were each sampled independently from the empirical gene-model distribution, preserving the empirical distributions of individual parameters and their presence/absence ratios
- `pdf`:
  > the probability density function (PDF) sampling strategy, built on mixed sampling; for the parameters {pdf_params}, values were resampled from histogram-based estimates of the empirical probability density (Freedman-Diaconis binning) with noise scale {noise_scale}

Where `{override_sentence}` (only if `--override` used):
> The following parameters were fixed across all simulations: {override_params}.

Where `{tree_sentence}`:
> Reference trees were sampled {independently from the model parameters | together with all model parameters (complete strategy)} from the empirical gene tree set.

### 9.4 `posttree.simulate.alisim.transfergaps`

Single mode (`n_msas == 1`):
> Gap patterns from the original alignment were transferred to the simulated alignment to restore biologically realistic indel patterns. Positions containing non-standard characters (gaps and ambiguity codes) in the original sequences were replaced with gap characters (-) at the corresponding positions in the simulated sequences. The valid character set was {AA: ACDEFGHIKLMNPQRSTVWY | NT: ACGT}; all other characters were treated as gaps.

Batch mode (`n_msas > 1`):
> Gap patterns from the original alignment were transferred to {n_msas} simulated alignments to restore biologically realistic indel patterns. Positions containing non-standard characters (gaps and ambiguity codes) in the original sequences were replaced with gap characters (-) at the corresponding positions in the simulated sequences. The valid character set was {AA: ACDEFGHIKLMNPQRSTVWY | NT: ACGT}; all other characters were treated as gaps.

---

## 10. MCP Integration

Three MCP tools auto-generated via the Click-tree walker in `schema_gen.py`:
- `posttree_simulate_alisim_params`
- `posttree_simulate_alisim_iqtree`
- `posttree_simulate_alisim_transfergaps`

The existing `posttree_simulate` stub in `stubs.py` is removed once the `simulate` group is registered.

---

## 11. Click Hierarchy Changes

### 11.1 New Groups

```python
class _SimulateGroup(click.Group):
    """Simulation commands."""
    def list_commands(self, ctx):
        return ["alisim", "adequacy", "phybase"]

class _AlisimGroup(click.Group):
    """IQ-TREE AliSim simulation commands."""
    def list_commands(self, ctx):
        return ["params", "iqtree", "transfergaps"]
```

### 11.2 Registration Changes

- Add `"simulate"` to `_PosttreeGroup.list_commands()`.
- Wire `_SimulateGroup` under `posttree`.
- Wire `_AlisimGroup` under `simulate`.
- `adequacy` and `phybase` remain as stub subcommands (print "not yet implemented") until their specs are written.

### 11.3 New Library Files

```
phyloai/posttree/
├── simulate_alisim_params.py     # parameter extraction logic
├── simulate_alisim_iqtree.py     # simulation execution logic
└── simulate_alisim_transfergaps.py  # gap transfer logic
```

CLI commands added to `phyloai/cli/commands/posttree.py`.

### 11.4 MCP Stub Update

Remove `"posttree_simulate"` from `STUB_TOOL_NAMES` in `stubs.py`. The three new leaf commands auto-generate MCP tools.

---

## 12. Dependencies

- **`alisim params`:** None (pure Python parsing)
- **`alisim iqtree`:** IQ-TREE3 (external binary) — same dependency as `tree ml iqtree`
- **`alisim transfergaps`:** None (pure Python)
- **PDF density estimation:** `numpy` (already a dependency for histogram computation)
- **Density plots:** `matplotlib` (already a dependency for metrics plots)
- **No new pip dependencies**

---

## 13. Test Data

- AA `.iqtree` files: `runs/simulate/faa/logs/*.iqtree` (800 files)
- NT `.iqtree` files: `runs/simulate/fna/logs/*.iqtree` (800 files)
- FreeRate example: `runs/simulate/faa/logs/EOG090X0064.iqtree` (LG+R2)
- FreeRate+I example: `runs/simulate/fna/logs/EOG090X01AN.iqtree` (GTR+F+I+R2)
- PMSF partition file: `runs/tree/bi/readpb/partition.PMSF.nex`
- Gap transfer reference: `ref/scripts/transfer_gaps.py`
- Sampling strategy reference: `ref/scripts/server.R`
