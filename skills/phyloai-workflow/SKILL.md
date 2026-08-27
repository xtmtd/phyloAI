---
name: phyloai-workflow
description: >-
  Guide users through PhyloAI CLI analyses through the local MCP server.
  Use for PhyloAI analysis, command execution, run recovery, systematic-error
  and heterogeneity workflows, `doctor`, environment checks, missing external
  tools, installation, and external-tool failures.
---

# PhyloAI Workflow

## Core Rules

- Use `doctor` before commands that invoke external tools, on first run, or when the environment is unknown.
- Use this Skill for environment and installation requests too, including `doctor failed`, `missing pb_mpi`, `install iqtree`, `缺少 MAFFT`, `环境检查失败`, and similar external-tool setup questions.
- Read-only tools (`check_status`, `read_result`, `read_report`, `get_command_schema`) do not require `doctor` first.
- **Execution tools vs workflow:** The MCP execution tools (e.g. `phyloai_posttree_signal_lnl`, etc.) are thin wrappers that launch CLI subprocesses directly — they DO NOT enforce schema review, `doctor`, or user approval on their own. This Skill is the process layer: before calling an execution tool, always run `get_command_schema`, render a parameter card, run `doctor` if external tools are needed, and wait for explicit user approval. Do NOT call an execution tool just because its MCP schema is available.
- Before executing a CLI command, call `get_command_schema`, render a parameter card that lists **every** parameter from the schema, and wait for explicit user approval. Do not filter out parameters — annotations in `references/parameter-annotations.md` are decorations, not a display filter. Parameters without annotations must still be shown with their CLI `--help` text. For every parameter, show both the current value and the schema default (e.g. `--threads  4  (默认: 4)  ...`). If the schema marks a parameter as required, it MUST have an explicit value before approval — do not launch with an unset required parameter, including conditionally required ones like `--matrix` for `tree bi pb`.
- After launching a fire-and-forget command, call `check_status` to verify the job actually started before declaring success. Do not claim "已启动" based solely on the launch response — the subprocess may have exited immediately. Report the `check_status` result to the user; if the status is `error` or `unknown`, show the error details and suggest next steps.
- Treat `--overwrite` as destructive. When the target `--output-dir` already exists and the user has not explicitly requested overwrite, prefer suggesting a new `--output-dir` or `--resume` when available before offering `--overwrite`. If a parameter card sets `--overwrite true`, ask for separate explicit confirmation naming the affected `--output-dir`; general command approval is not enough.
- Never invent parameter names, aliases, defaults, or enum values. Unknown parameters block execution.
- After a command completes, summarize `key_results`, warnings, and next steps. Do not auto-run the next step.
- When a user asks about progress for a running job, call `check_status` and summarize using the checkpoint or result state. For `tree bi pb`, `convergence/convergence_render.txt` contains human-readable convergence diagnostics (pairwise chain status, after-burnin sample counts); this is the primary progress indicator once the first convergence check completes.

## Entry Modes

- New task: ask for input data path, run `doctor` if needed, then start pretree workflow.
- Resume task: call `read_report(run_dir)`; if missing, ask whether to run `report` or inspect a specific step with `read_result`.
- Recovery is verified, not assumed: when resuming or inspecting an existing `output_dir`, smoke-test the read path with `check_status` → `read_result` → (when requested) `read_report`, and report the observed state. A completed atomic run (`result.json` present, `status: success`) is NOT re-run merely because the Skill conversation was interrupted — show the existing result and only relaunch if the user explicitly asks for a fresh run or grants `--overwrite`.
- Single-step task: render the parameter card for the requested command and wait for confirmation.

## Language Policy

- Parameter cards use English parameter names with Chinese annotations and recommendations.
- Conversation and interpretation follow the user's language.
- CLI commands are shown in English exactly as executable commands.

## Workflow

- Pretree: `convert -> align -> trim -> metrics / filter -> concat` (supermatrix) or `... -> gene trees` (supertree). `stats` inspects results at any step.
- Tree: `tree ml iqtree` + `tree msc` as primary, `tree ml fasttree` for fast exploration, `tree bi pb` for Bayesian MCMC, `tree bi bpcomp`/`tree bi tracecomp` for final convergence diagnostics with user-chosen burn-in, `tree bi readpb` for posterior summaries. For custom CAT-PMSF-style ML, pass an AA exchangeability file with `tree ml iqtree --model`, a profile with `--site-freq-file`, and `--state-freq none`; raw `--tool-args -fs` overrides the structured profile. For PMSF simulation input, use `tree bi readpb --mode ss,rr,r`; it writes `partition.PMSF.nex` from posterior site rates, alpha, Gamma category count, and the co-generated `.exchangeabilities` model. Use `cf` on species trees.
- Posttree: `topology`, `dating hessian`, `dating mcmc`, `signal lnl`, `signal consistent`, `signal fclm`, `modelcompare iqtree`, `modelcompare pb`, `simulate alisim`, `simulate adequacy`.
- Systematic error / heterogeneity: when users ask about LBA, long branches, rate/composition heterogeneity, heterotachy/GHOST, CCA, model adequacy, or systematic-error correction, load `references/syserror-workflow.md` before selecting a path. Explain the relevant theory, offer only the applicable subset of the six optional analyses, and do not run a full diagnostic battery by default. The reference governs prerequisites, evidence limits, sensitivity options, and advanced simulations; existing schema review, doctor, explicit approval, overwrite, status, and recovery rules still apply one command at a time.
- Report: run `report` only when the user requests a report/methods draft or recovery needs `report.json`.

### posttree signal

Three subcommands for phylogenetic signal distribution analysis. All are
single-matrix (no batch mode). Model source: `--model-expr` and `--partitions`
can be combined — `--model-expr` specifies the model formula and `--partitions`
provides partition boundaries; each partition independently estimates parameters.

#### signal lnl

Purpose: Site-wise and gene-wise log-likelihood score distribution across
candidate trees using IQ-TREE3 `-wslr`. Identifies outlier genes with
disproportionate signal (Shen et al. 2017, *Nature Ecology & Evolution*).
When `--metrics` is provided with locus boundaries, also compares metrics
across gene groups supporting different candidate trees.

Required inputs: `--matrix`, `--candidate-trees`, plus at least one model source
(`--model-expr`, `--partitions`, or `-m`/`-p` in `--tool-args`).

Optional: `--locus-ranges` for gene-wise breakdown + outlier detection;
`--metrics` for outlier-vs-nonoutlier and tree-support-group pairwise comparisons.
`--guide-tree` for PMSF models.

Mutual exclusions: `--partitions` vs `--locus-ranges`.

Outputs: `site_lnl.csv`, `support_summary_sites.csv`, `gene_lnl.csv` (if boundaries), `support_summary_genes.csv` (if boundaries), `outlier_genes.txt`,
plots, `support_comparison.csv/pdf` (if `--metrics` + >=2 support groups), `result.json`. IQ-TREE files go in `iqtree/` subdirectory.

#### signal consistent

Purpose: Consistent gene identification where both GLS (likelihood-based) and
GQS (quartet-based) agree on supporting one of two candidate topologies.
Uses IQ-TREE3 for GLS + wASTRAL for GQS (Shen et al. 2021, *Systematic Biology*).

Required inputs: `--matrix`, `--candidate-trees` (exactly 2 trees), `--tree-dir`.
At least one of `--partitions` or `--locus-ranges` for GLS.

Optional: `--metrics` for consistent-vs-inconsistent comparison;
`--partition-mode` (`p` or `Q`) when `--partitions` is provided.

Mutual exclusions: `--partitions` vs
`--locus-ranges`. Exactly 2 candidate trees only. Partition loci must all
have matching gene tree files in `--tree-dir` (extra trees ignored).

Outputs: `gls.csv`, `gqs.csv`, `consistent_genes.txt`, `inconsistent_genes.txt`,
plots, `result.json`. IQ-TREE files go in `iqtree/` subdirectory.

#### signal fclm

Purpose: Four-cluster Likelihood Mapping (FcLM) to assess signal supporting
alternative hypotheses of relationship among four taxon clusters.
Uses IQ-TREE3 `-lmap -lmclust`.

Required inputs: `--matrix`, `--taxset-csv` (at least 4 taxsets). Model source:
`--model-expr` and/or `--partitions`. `--partition-mode` (`p`
or `Q`, default `p`) controls how `--partitions` is passed to IQ-TREE.

Optional: `--lmap` (quartet count; default `50 * n_taxa`); `--guide-tree` for
PMSF models.

Validation: all taxa in CSV must match matrix taxa; taxset assignments must be
mutually exclusive; minimum 4 taxsets.

Outputs: IQ-TREE native `.iqtree` report (contains all lmap statistics),
`<prefix>.lmap.eps` figure, `result.json`. IQ-TREE files in `iqtree/` subdirectory.

### posttree modelcompare

Two subcommands for relative substitution-model comparison and selection.

#### modelcompare iqtree

Purpose: ModelFinder (BIC/AIC/AICc) model comparison on a single alignment
using IQ-TREE3 `-m MF` (Kalyaanamoorthy et al. 2017). Reports each model's
BIC/AIC/AICc scores, weights, and 95% confidence-set membership.

Required inputs: `--matrix`, `--homogeneous-model` (comma-separated standard
models, maps to `-mset`).

Optional: `--mrate` (E/G/R, default `E,G`), `--heterogeneous-model` (AA mixture
models via `-madd`; AA only), `--het-mrate` (any subset of E/G/R, default `E,G`;
each token selects a variant family — E = base models `C10,C10+F`, G = `+G4`,
R = `+R4`, mirroring `--mrate`), `--seq-type`
(`AA`/`NT`/`auto`, default `auto` — alignment is read to detect AA vs NT before
model validation; an explicit value is cross-checked against the detected
type), `--prefix` (must be a single filename — no `/`, `..`, or absolute
paths), `--threads` (positive integer or `auto`), `--tool-args`
(blocked: `-s`, `--prefix`; overrides PhyloAI-managed flags like `-madd`,
`-mrate`, `-T` when present), `--resume`.

Mutual exclusions: `--overwrite` vs `--resume`. Heterogeneous models rejected
for NT data. Non-empty existing output dir requires `--overwrite`.

Outputs: `model_fit.csv` (sorted by BIC, with `In_*_95` confidence-set columns),
`result.json`. IQ-TREE files in `iqtree/` subdirectory.

#### modelcompare pb

Purpose: LOO-CV / wAIC model comparison from PhyloBayes `.sitelogl` site
log-likelihood files (Lartillot 2023). Pure Python — no external tool.

Required inputs: exactly one of `--sitelogl-dir` (comma-separated dirs, each
globbing `*.sitelogl`) or repeated `--sitelogl` (each occurrence = one model's
chain files). At least 2 `.sitelogl` files per model group. Both options use
the Click `Path` type, so shell path completion works.

Optional: `--model-names` (labels matching group count; unique; each must be a
safe path component — no `/`, `..` — as it names an output subdirectory).

Validation: within each group all files must share identical ordered `site`
identifiers; across groups all must share identical site count AND ordered
`site` identifiers (mismatch → hard error). Scores from different alignments
are not comparable. Non-empty existing output dir requires `--overwrite`.

Outputs: `model_fit.csv` (LOO-CV/wAIC score, bias, stdev, CI95, ESS quality
columns, Delta values for multi-model), copied `.sitelogl` files under
`sitelogl/model_N/`, `result.json`.

### posttree simulate alisim

Three-subcommand IQ-TREE3 AliSim workflow that preserves empirical dataset
properties instead of using arbitrary model parameters:

- `alisim params` — extract per-locus simulation parameters (model string
  components, frequencies, proportion of invariant sites, rate heterogeneity,
  tree) from existing `.iqtree` reports into `params.tsv`.
- `alisim iqtree` — simulate MSAs. Single mode takes `--ref-tree` + a model
  source (`--model` or `--model-partitions`). Batch mode takes `--model-params`
  + `--strategy complete|mixed|pdf` (default `complete`) +
  `--num-simulations`, samples rows as atomic units (model core + rate group
  stay together; `+I` presence decided before its value is resampled), is
  resumable via `--resume`, and writes every used row to `params_sampled.tsv`
  (`source_id` column only for the `complete` strategy; each simulation gets an
  independent random seed from the master-seeded generator). `--override` (keys
  `length`, `prop_inv`) fixes values across all simulations. PDF density plots
  (`server.R` palette) are generated only under `--strategy pdf`.
  `--tool-args` may override non-I/O IQ-TREE flags (e.g. `--length`,
  `--num-alignments`); only `--alisim`, `-t`, `--prefix`, `--out-format`,
  `-af` are blocked.
- `alisim transfergaps` — transfer the original per-taxon gap mask onto one
  simulated (gap-free) MSA (`--simulated-msa`) or a directory of simulated
  MSAs (`--simulated-dir`, exactly one required, mutually exclusive). Taxon
  sets and lengths must match. Batch mode writes `<stem>.gaps.fa` per input
  (output is always FASTA).
  `--exclude-ambiguity` transfers only `-`/`.` instead of all non-standard
  characters.

`simulate adequacy` compares observed PPA-DIV, PPA-CONV, PPA-VAR, and
PPA-COMP statistics with at least 10 simulated MSAs. It independently detects
the supported alignment format for every input file, writes CSV or TSV tables,
and can resume from its checkpoint. Recommend `transfergaps` first when
observed missing data are substantial. A low pp (`< 0.05`) or `|z| > 2`
indicates potential inadequacy; `div` pp is `P(sim <= obs)`, while the other
statistics use `P(sim > obs)`. `phybase` remains a future placeholder.

### posttree syserror taxcomp

Across-taxon composition heterogeneity screen, pure Python (no external tool,
so no `doctor` needed — but still review parameters and get explicit user
approval before running).

- **Input selection:** exactly one aligned MSA (`--matrix`); no tree or model
  input required. `--seq-type auto` resolves AA/NT from the alignment;
  standard states are counted and gaps/ambiguity codes are excluded.
- **Approval:** requires the full schema card and explicit execution approval;
  it is a read-only diagnostic that writes only two summaries plus
  `result.json`.
- **Results interpretation:** inspect `sparse_count_check` first — a
  `triggered` sparse-cell rule warns that the nominal chi-square p-values are
  especially unreliable. Treat the overall and per-taxon `p_nominal` and
  `p_holm` values as nominal/exploratory, never as a pass/fail verdict. Use
  the taxon `x2_contribution` and `squared_composition_distance` only to
  prioritize which taxa to inspect (annotation, coverage, contamination,
  lineage composition). Do not delete taxa from this screen alone; model-
  calibrated PPA-COMP requires `posttree simulate adequacy`. Optional recoding
  is a separate, user-approved sensitivity analysis with `Dayhoff-6` for AA or
  `RY-nucleotide` for NT.

### posttree syserror brlen

Branch-length extraction for rate-heterogeneity / LBA diagnosis, pure Python
(no external tool, so no `doctor` needed — but still review parameters and get
explicit user approval before running).

- Batch modes (`total`, `terminal`, `internal`, `patristic`, `all`) extract
  per-tree branch-length statistics. Endpoint modes (`tip-to-tip`,
  `node-to-node`, `node-to-tip`) compute distances between specified
  tips/nodes; endpoint modes are mutually exclusive with batch modes.
- **Node identification:** for batches whose topology or taxon set differs
  across trees (posterior/model tree batches), always use a `--map` file —
  rooted maps require an exact monophyletic clade, unrooted maps an exact split
  side. Use `label-nodes` only on an inspected, stable reference topology; `Nxx`
  labels are not portable across topologies.
- `patristic` is O(n²) per tree; review the `--max-rows` estimate (default
  5,000,000, `0` disables) for large posterior batches.
- Rootedness is structural (root children: 2 = rooted, 3+ = unrooted); state
  this caveat when interpreting `internal`/`label-nodes`/`node-to-tip` results.
- In batch endpoint mode, a tree whose tip/node/map endpoint cannot be resolved
  is skipped with a warning and counted in `n_trees_skipped`. `--dry-run` runs
  the same endpoint resolution, so it rejects single-tree unresolved endpoints
  before any real run. Batch processing shows a transient progress bar (visible
  on large posterior batches; `--threads` parallelizes `--tree-dir`).
- **Results interpretation:** distinguish terminal vs internal vs patristic
  measurements, and node-to-node vs node-to-tip distances. Compare these
  distributions across model runs (e.g. LG vs CAT-PMSF posterior trees). Do not
  declare model superiority from branch lengths alone; combine with `cca` and
  `taxcomp` diagnostics (for across-taxon composition) and `rate` (for
  site-rate heterogeneity).

### posttree syserror rate

Local pure-Python site-rate ranking/extraction sensitivity utility. No
`doctor` check is needed, but before execution still call
`get_command_schema`, render the full parameter card, and obtain explicit user
approval.

- **Input selection:** require exactly one rate source: IQ-TREE3 `--rate`
  output via `--iqtree-rate`, or PhyloBayes `readpb -r` `.meansiterates` via
  `--pb-rate`. IQ-TREE indices are 1-based; PhyloBayes indices are normalized
  from 0-based to 1-based. Both sources must be consecutive, producing exactly
  `1..N` sites.
- **Extraction:** `--matrix` is optional for ranking only. With a matrix,
  `--fraction` is required and `--subset slow|fast` selects the direction
  (`slow` is the extraction default). Recommend a sensitivity series such as
  `0.25,0.5,0.75`, not an asserted optimal threshold.
- **Interpretation:** slow subsets reduce rapidly evolving-site contribution;
  fast subsets isolate it. Neither establishes a topology or an automatic
  systematic-error correction. Downstream tree inference/comparison is always
  a separate user-approved action.

### posttree syserror cca

Local pure-Python compositional-constraint diagnostic. No `doctor` check is
needed, but before execution call `get_command_schema`, render every runtime
parameter in the parameter card, and obtain explicit user approval.

- **Prepared inputs:** require one IQ-TREE PMSF or `readpb --mode ss`
  `.sitefreq` input and two independent model-specific `site_lnl.csv` files.
  The frequency table is one-based `1..N` with 20 frequencies per site; each
  LNL table must use literal `site`, `lnL_Tree1`, and `lnL_Tree2` headers. The
  three site sets must exactly match. Raw `.siteprofiles` and custom LNL-header
  labels are not accepted.
- **Model comparison:** `--model1-name` and `--model2-name` label the two LNL
  analyses; they are non-empty and distinct. CCA calculates Keff as inverse
  homozygosity and independently calculates Tree2-minus-Tree1 likelihood
  differences, rather than using signal LNL's reverse-sign `ΔSLS` column.
- **Interpretation:** positive bars support Tree2 and negative bars support
  Tree1 within a model; the plot can diagnose composition-associated changes
  in preference but must not be used alone to declare a preferred topology or
  model. Tree/model inference and follow-up analyses remain separate,
  user-approved actions.

**Approval:** `alisim iqtree` executes IQ-TREE3, so run `doctor`, review the
full parameter card, and get explicit user approval before invoking it —
including `--dry-run` to show the sampling plan/command first. `alisim params`,
`alisim transfergaps`, and `simulate adequacy` are local-only (no external tool)
but still get the standard parameter review. `posttree syserror taxcomp`,
`posttree syserror brlen`, `label-nodes`, `posttree syserror rate`, and
`posttree syserror cca` are
likewise local-only: run `get_command_schema`, render a parameter card, and
get explicit user approval before running. `check_status`,
`read_result`, `read_report`, and `get_command_schema` remain directly
inspectable without approval.

## Demo Data

- Bundled demo data contains 20 genes and 6 species, with AA (`faa/`) and NT (`fna/`) raw inputs.
- Per-step demo directories include aligned, trimmed, concatenated, gene-tree, topology-test, and dating entry points when present.
- Demo runs must write new outputs to a user run directory, not back into `demo_data/`.

## Error Handling

- Exit 1/3: use `references/error-catalog.md`, then show a fix card.
- Exit 2: diagnose tool stderr. For batch commands, inspect only failed loci logs, capped at about 10 loci, and state when truncated.

## References

- `references/parameter-annotations.md`
- `references/error-catalog.md`
- `references/dialog-templates.md`
- `references/demo-data.md`
- `references/workflow.md`
- `references/syserror-workflow.md` for optional systematic-error and heterogeneity workflows
- `docs/commands/installation.md` for external-tool setup guidance
