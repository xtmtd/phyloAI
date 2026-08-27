# PhyloAI Systematic Error Workflow Design

**Date:** 2026-08-28  
**Last updated:** 2026-08-28  
**Status:** Draft — pending approval  
**Parent spec:** `2026-06-07-phyloai-design.md`  
**Related command specs:**

- `2026-08-11-phyloai-posttree-syserror-brlen-design.md`
- `2026-08-12-syserror-rate-design.md`
- `2026-08-13-phyloai-posttree-syserror-cca-design.md`
- `2026-08-25-phyloai-posttree-syserror-taxcomp-design.md`
- `2026-08-02-phyloai-posttree-modelcompare-design.md`
- `2026-08-02-phyloai-posttree-simulate-alisim-design.md`
- `2026-08-03-phyloai-posttree-simulate-adequacy-design.md`
- `2026-06-19-phyloai-tree-ml-iqtree-design.md`
- `2026-07-26-phyloai-tree-bi-subcommands-design.md`

---

## 1. Purpose

Define an AI-guided workflow for exploring, diagnosing, testing the sensitivity
of, and addressing common substitution-process systematic errors in
phylogenomics. Complete analysis normally requires several atomic PhyloAI
commands plus user decisions; no single `phyloai` command should claim to make
the complete diagnosis.

This design composes existing atomic commands. It does not add a one-click
`diagnose-all` command, a new top-level Skill, or new MCP execution tools.

All six analysis types are optional. The Skill first explains the relevant
theory, then helps the user choose analyses based on the dataset, focal clades,
competing topologies, and available compute. It must not run all six analyses by
default.

---

## 2. Scope and Evidence Policy

Most phylogenetic substitution models make simplifying assumptions about how
rates, equilibrium compositions, and substitution processes vary across sites,
taxa, lineages, and time. Violations can produce systematic error whose support
may increase rather than disappear as alignments become longer.

The workflow covers six traditional categories:

| # | Error type | Main dimension | Primary PhyloAI route |
|---|---|---|---|
| 1 | Heterogeneity of rates across taxa | Across branches/lineages | `posttree syserror brlen` |
| 2 | Heterogeneity of rates across sites | Across alignment columns | `posttree syserror rate` |
| 3 | Heterotachy | Site × lineage/taxon | `tree ml iqtree` with reviewed GHOST `--tool-args` |
| 4 | Compositions across taxa | Across taxa/lineages | `posttree syserror taxcomp` |
| 5 | Compositions across sites | Across alignment columns | `posttree syserror cca` |
| 6 | Substitution patterns across sites | General model fit/adequacy | `posttree modelcompare` and posterior prediction |

These names follow common systematic-error and heterogeneity literature. The
sixth workflow operationally tests broader or residual substitution-model
misspecification. It must not imply that CAT or C10–C60 assign a different
exchangeability matrix to every site: these models principally mix
site-frequency profiles while sharing exchangeabilities within an analysis.

### 2.1 Evidence vocabulary

The Skill must distinguish the following evidence levels:

1. **Screen:** identifies unusual taxa, sites, branch lengths, compositions, or
   model-fit patterns.
2. **Consistent with:** the observation matches a proposed systematic-error
   mechanism but does not establish causality.
3. **Sensitivity-supported:** a controlled model, taxon, site, recoding, or
   missing-data comparison changes the result in the predicted direction.
4. **Simulation-supported:** simulations under explicit assumptions reproduce
   the predicted behavior.

No single branch-length, topology, composition, likelihood, or simulation
result establishes a biologically correct topology. Avoid causal wording such
as `confirmed`, `proved`, `the true tree`, or `adequately handled` unless the
specific evidence justifies it.

### 2.2 Adjacent but separate domains

ILS, HGT, introgression, recombination, gene duplication/loss, natural
selection, saturation, missing data, alignment/orthology errors, and other
processing errors may generate similar conflicts. They are not silently folded
into this substitution-systematic-error workflow. Missing data is considered
here only as a controlled simulation factor, not as a seventh systematic-error
type.

---

## 3. Guided Workflows

Before proposing any command, the Skill gives the short theoretical explanation
in the corresponding section, asks what biological contrast or topology is
being tested, obtains the command schema, shows every parameter, and waits for
explicit approval under the existing `phyloai-workflow` execution rules.

### 3.1 Heterogeneity of rates across taxa

#### Theory

Branch length is the expected number of substitutions per site. A long branch
may reflect elapsed time, a high evolutionary rate, model-dependent branch
length estimation, or a combination of these factors. Unequal rates across
lineages can increase long-branch-attraction risk, but long branches alone do
not demonstrate that LBA caused a topology.

Terminal branches describe individual taxa. Internal branches describe clade
separation and short internodes. Map-defined node-to-tip or node-to-node
distances can represent biologically meaningful quantities such as tip-to-stem,
tip-to-clade-root, or stem-to-stem distances. Pairwise patristic distances are
not part of the standard rates-across-taxa workflow.

#### Screening and sensitivity flow

1. Define focal taxa/clades and the branch-length quantity before execution.
2. Use `phyloai posttree syserror brlen` to extract relevant terminal, internal,
   node-to-tip, or node-to-node lengths.
3. Compare the same quantity under homogeneous and heterogeneous models.
   Heterogeneous examples include CAT-LG and CAT-GTR in the Bayesian framework,
   and CAT-PMSF-style downstream ML analyses. CAT-PMSF is an example, not the
   only heterogeneous model.
4. Prefer distributions from posterior trees over a comparison of only one ML
   tree per model. Use `phyloai tree bi pb` with a fixed, explicitly chosen
   topology and model, then apply `brlen` to the resulting posterior trees.
   Remind the user to inspect convergence and normally discard an appropriate
   burn-in before preparing the posterior-tree input; optional thinning and its
   value remain user decisions. The Skill does not filter posterior trees or
   select burn-in/thinning values. A fixed-topology CAT-LG or CAT-GTR analysis
   is the Bayesian estimation stage used by CAT-PMSF-style workflows; the
   terminology must not imply that all heterogeneous analyses are CAT-PMSF.
5. Interpret differences as model-dependent branch-length behavior and LBA
   risk. They do not by themselves establish model superiority or topology
   correctness.

#### Resolution options

Present, but never auto-select:

- infer trees with an appropriate heterogeneous model, such as CAT-LG,
  CAT-GTR, or a CAT-PMSF-style approximation;
- perform a user-defined taxon-sampling sensitivity analysis with and without
  suspected long-branched taxa;
- combine model and taxon-sampling sensitivity analyses.

Taxon removal requires biological justification. A branch-length ratio can be
an exploratory threshold, not an automatic deletion rule.

#### Optional advanced simulation

This step is displayed only when the user asks for simulation or accepts it
after the Skill explains its computational cost.

**Primary posterior-predictive route:**

1. Run fixed-topology homogeneous and/or CAT-LG/CAT-GTR PhyloBayes analyses.
2. Use `phyloai tree bi readpb --mode ppred` to simulate posterior-predictive
   replicate MSAs from the selected chain. This is the preferred posterior
   predictive route because parameters are sampled through PhyloBayes rather
   than replaced by posterior means.
3. Run `phyloai posttree simulate adequacy` with the observed MSA and the full
   replicate set. It compares one observed MSA with a distribution across
   simulated MSAs; it is not a separate adequacy test for each replicate.
4. Re-estimate trees or fixed-topology branch lengths under the selected models
   and compare branch-length distributions with `syserror brlen`.

**Optional AliSim plug-in route:**

1. For a compatible CAT analysis, use
   `phyloai tree bi readpb --mode ss,rr,r` to generate posterior-mean site
   frequencies, exchangeabilities, rates, and `partition.PMSF.nex`.
2. Simulate gapless replicates with `phyloai posttree simulate alisim iqtree`.
   This is a posterior-mean plug-in parametric simulation, not a strict PPC.
3. Use `phyloai posttree simulate alisim transfergaps` to apply the observed
   per-taxon gap mask to a second replicate set.
4. Analyze both gapless and gap-transferred replicates. Their contrast can test
   whether the observed missing-data pattern materially changes the simulated
   inference.

The `ss,rr,r` AliSim route requires a chain/model for which all three outputs
and `partition.PMSF.nex` are valid. Do not request free exchangeabilities from
a fixed-LG chain or invent unavailable parameters.

---

### 3.2 Heterogeneity of rates across sites

#### Theory

Sites can evolve at different rates because of structural, functional, or
selective constraints. A single-rate model may underestimate repeated changes
at fast sites. Gamma (`+G`) and FreeRate (`+R`) models approximate this
variation, but rate estimates remain conditional on the fitted model and tree.
Slow sites may contain less information, while fast sites may contain more
homoplasy; neither class is automatically preferable.

#### Screening and sensitivity flow

1. Produce site rates with IQ-TREE `--rate` or PhyloBayes `readpb -r`.
2. Rank and extract subsets with an explicit rate source and matrix, for example:

   ```bash
   phyloai posttree syserror rate \
     --iqtree-rate matrix.rate \
     --matrix matrix.fa \
     --subset slow \
     --fraction 0.25,0.5,0.75
   ```

   Repeat with `--subset fast` when isolating fast-site behavior is relevant.
3. Rebuild trees from a user-chosen series of retained slow and/or fast subsets.
4. Separately compare inference with and without an across-site rate model, for
   example `LG+F` versus `LG+F+G4` or `LG+F+R4`.
5. Compare topology, support, focal branch lengths, and likelihood behavior.

A topology change after subsetting is consistent with rate-associated
sensitivity, but may also reflect reduced information, composition shifts, or
model-conditional rate estimates. A stable topology indicates no conspicuous
topological effect at the tested fractions; it does not exclude smaller effects
on support, branch lengths, or likelihood.

#### Resolution options

- use an appropriate `+G` or `+R` rate model;
- retain or exclude a user-justified slow or fast site class in a sensitivity
  analysis;
- use biologically defined partitions where appropriate.

The workflow must not assume that only fast sites may be removed. Thresholds
are sensitivity choices rather than universally optimal cutoffs.

---

### 3.3 Heterotachy

#### Theory

Heterotachy is rate variation across both sites and lineages/taxa: the same site
may evolve rapidly in one lineage and slowly in another. Standard `+G` and `+R`
models allow rates to differ among sites but generally keep each site's relative
rate constant across the tree. GHOST (General Heterogeneous evolution On a
Single Topology) provides an IQ-TREE route for modeling this site-by-lineage
rate variation.

#### Guided flow

1. Consider this workflow when the rates-across-taxa and rates-across-sites
   screens, prior knowledge, or the study design make heterotachy plausible.
2. Run GHOST through `phyloai tree ml iqtree --tool-args` using an exact,
   user-reviewed IQ-TREE model expression appropriate for the sequence type and
   installed IQ-TREE version.
3. The raw `-m` supplied through `--tool-args` overrides PhyloAI's managed model
   string. The parameter card must display the complete raw expression. The
   Skill must not guess GHOST syntax, component count, state-frequency option,
   or defaults.
4. Interpret the result as an explicit heterotachy-aware inference, not as a
   standalone proof that heterotachy caused an alternative topology.

No AU-test or `posttree modelcompare iqtree` step is prescribed for this
workflow. Current structured model validation and `modelcompare iqtree` do not
provide a dedicated GHOST comparison contract.

GHOST does not solve across-site profile heterogeneity or lineage-specific
composition heterogeneity. CAT-LG and CAT-GTR likewise do not model
heterotachy. Models intended to combine these dimensions are outside the
current design.

---

### 3.4 Compositions across taxa

#### Theory

Across-taxon composition heterogeneity occurs when lineages differ in observed
nucleotide or amino-acid composition. A stationary homogeneous model assumes a
shared equilibrium composition across the tree; violation can group unrelated
lineages with similar compositions. Observed differences may also reflect
contamination, annotation problems, missing data, or genuine lineage-specific
evolution.

#### Screening and sensitivity flow

1. Run `phyloai posttree syserror taxcomp --matrix <alignment>`.
2. Inspect `sparse_count_check` before interpreting nominal chi-square
   p-values.
3. Use overall and per-taxon chi-square contributions plus PPA-COMP descriptive
   distances to prioritize taxa for biological and data-quality inspection.
4. Treat nominal and Holm-adjusted p-values as exploratory screens, not as
   phylogenetically calibrated pass/fail decisions.
5. If justified, perform recoding or taxon-sampling sensitivity analyses and
   rebuild the tree.

#### Resolution options

- inspect annotation, contamination, coverage, and missing data for flagged
  taxa;
- perform a biologically justified taxon-sampling sensitivity analysis;
- recode with `phyloai pretree concat --recoding Dayhoff-6` for amino acids or
  `--recoding RY-nucleotide` for nucleotides, then rebuild the tree;
- consider nonstationary composition models outside the current PhyloAI
  implementation.

Profile-mixture models must not be presented as a direct solution for
composition heterogeneity across taxa: their site-frequency profiles primarily
address heterogeneity across sites and are shared across lineages. P4 and GFmix
may be mentioned as planned/external approaches for more explicit across-lineage
composition modeling, but this Skill and the current PhyloAI tools do not
implement or deeply guide them.

---

### 3.5 Compositions across sites

#### Theory

Different alignment columns can favor different amino acids or nucleotides due
to structural and functional constraints. A global stationary-frequency vector
can underestimate convergent substitutions at strongly constrained sites.
Profile-mixture models approximate this across-site variation. CCA asks whether
site-wise preference between two fixed topologies changes with effective amino
acid diversity (`Keff`) and with the fitted model.

CCA requires one `.sitefreq` table and two model-specific `site_lnl.csv` files.
Each likelihood table must contain the same two candidate trees in the same
order and therefore include `lnL_Tree1` and `lnL_Tree2` in one file.

#### Prepare `.sitefreq`

Use either source:

1. `phyloai tree bi readpb --mode ss` from a PhyloBayes CAT analysis; or
2. an IQ-TREE PMSF run that writes `.sitefreq`, including the classic PMSF
   preparation pattern:

   ```bash
   iqtree -s <alignment> -m LG+C20+F+G -ft <guide_tree> -n 0
   ```

A `.sitefreq` retained from an earlier CAT-PMSF or IQ-TREE PMSF analysis can be
reused when its alignment and site order exactly match the CCA likelihood
inputs.

#### Prepare two likelihood tables

Use the same ordered `T1,T2` candidate-tree set in two independent
`phyloai posttree signal lnl` runs:

```bash
phyloai posttree signal lnl \
  --matrix matrix.fa \
  --candidate-trees T1.tre,T2.tre \
  --model-expr LG+F+R4 \
  --output-dir runs/syserror/cca-input/lnl-LG

phyloai posttree signal lnl \
  --matrix matrix.fa \
  --candidate-trees T1.tre,T2.tre \
  --model-expr <reviewed-heterogeneous-model-expression> \
  --output-dir runs/syserror/cca-input/lnl-heterogeneous \
  <reviewed profile options when required>
```

For a custom CAT-PMSF-style analysis, the second run may require the same custom
exchangeability model and `.sitefreq` profile through reviewed supported
options, including `--tool-args "-fs <sitefreq>"` where appropriate. The Skill
must obtain the exact command schema and reuse the model/profile combination
that generated the corresponding tree; it must not invent a generic command.

#### Run and interpret CCA

```bash
phyloai posttree syserror cca \
  --site-freq analysis.sitefreq \
  --site-lnl1 runs/syserror/cca-input/lnl-LG/site_lnl.csv \
  --site-lnl2 runs/syserror/cca-input/lnl-heterogeneous/site_lnl.csv \
  --model1-name LG \
  --model2-name heterogeneous
```

- Low `Keff` denotes strongly composition-constrained sites.
- Positive `delta_lnl_tree2_tree1` supports Tree2 within that model; negative
  values support Tree1.
- A model-associated change concentrated in low-`Keff` bins is consistent with
  composition-constrained sites affecting topology preference.
- A flatter or reduced pattern under a heterogeneous model suggests reduced
  sensitivity; it does not prove that the model is adequate or that one
  topology is correct.

#### Resolution options

- infer with an across-site profile-mixture model such as CAT-LG, CAT-GTR,
  C10–C60/PMSF, or a CAT-PMSF-style analysis;
- perform a user-defined low-`Keff` site-removal sensitivity analysis;
- combine profile-mixture inference with a site-filtering sensitivity analysis.

Current PhyloAI has no atomic command that filters arbitrary sites by a `Keff`
cutoff. Low-`Keff` extraction must be described as manual preprocessing or
future command work, not as an already implemented atomic workflow. Cutoffs
such as 5, 6, or 7 are exploratory and dataset-dependent.

---

### 3.6 Substitution patterns across sites

#### Theory

Substitution-model misspecification is broader than rate or composition alone.
Partitioning and mixture models can describe different aspects of across-site
heterogeneity, while relative model fit asks which candidate model is better
and absolute adequacy asks whether a model can reproduce important properties
of the observed data. A model can be best among candidates and still be
inadequate.

This sixth workflow focuses on general or residual substitution-model
misspecification. It overlaps biologically with compositions across sites but
uses model comparison and posterior prediction rather than attributing the
signal to `Keff`.

#### Relative model fit

For IQ-TREE AA analyses, compare explicitly named homogeneous and heterogeneous
candidates, for example:

```bash
phyloai posttree modelcompare iqtree \
  --matrix concat.aa.fa \
  --homogeneous-model LG \
  --heterogeneous-model C10,C20 \
  --het-mrate G,R
```

Use `phyloai posttree modelcompare pb` for PhyloBayes LOO-CV/wAIC comparisons
when the required site-log-likelihood inputs are available.

Lower AIC/AICc/BIC or better LOO-CV/wAIC supports relative preference among the
specified candidates only. It does not establish absolute adequacy or prove
that a topology is correct.

#### Optional absolute adequacy / simulation

**Primary posterior-predictive route:**

1. Fit the selected PhyloBayes model and verify chain sampling/convergence as
   required by its own workflow.
2. Use `phyloai tree bi readpb --mode ppred` to create posterior-predictive MSA
   replicates.
3. Run `phyloai posttree simulate adequacy` once with the observed MSA and the
   replicate directory to evaluate PPA-DIV, PPA-CONV, PPA-VAR, and PPA-COMP.
4. Interpret posterior predictive p-values and z-scores as statistic-specific
   model checks, not a single global pass/fail score.

**Optional AliSim plug-in route:**

For a compatible CAT model, use `readpb --mode ss,rr,r` to generate
`partition.PMSF.nex`, simulate with AliSim, and optionally apply
`simulate alisim transfergaps`. This uses posterior-mean plug-in parameters and
must not be labeled strict PPC. Comparing gapless and gap-transferred replicate
sets can expose missing-data sensitivity.

#### Resolution options

- use a better-supported homogeneous, partition, or profile-mixture model;
- use CAT-LG/CAT-GTR in PhyloBayes or C10–C60/PMSF/CAT-PMSF-style inference in
  IQ-TREE where scientifically appropriate;
- if all tested models remain inadequate, report that limitation and combine
  model-fit evidence with the other systematic-error sensitivity analyses.

---

## 4. Skill Organization

Do not create a new top-level `phyloai-syserror` Skill or a separately
registered subskill. Extend the existing version-coupled workflow Skill:

```text
skills/phyloai-workflow/
├── SKILL.md
└── references/
    └── syserror-workflow.md
```

### 4.1 `SKILL.md` responsibilities

Keep the main Skill concise. Add:

- triggers for systematic-error, heterogeneity, LBA, composition, CCA,
  heterotachy/GHOST, model adequacy, and related analysis requests;
- the six-choice selection prompt and the rule that all analyses are optional;
- the requirement to explain theory before parameter review;
- a link instructing the agent to load `references/syserror-workflow.md`;
- existing schema review, `doctor`, explicit approval, overwrite confirmation,
  status, and result-reading rules;
- recovery instructions based on existing structured outputs.

### 4.2 Reference responsibilities

`references/syserror-workflow.md` contains:

- the scientific introductions and evidence vocabulary;
- input-preparation dependencies;
- per-type diagnosis, sensitivity, simulation, and resolution choices;
- interpretation cautions;
- optional advanced-compute branches.

This structure avoids turning `SKILL.md` into a long scientific manual while
keeping one discoverable workflow Skill.

### 4.3 Superseded master-design decision

Update the parent design so that its future standalone `phyloai-syserror` Skill
entry is superseded by this reference module inside `phyloai-workflow`.

---

## 5. Execution, Provenance, and Recovery

No new workflow database, manifest, or orchestration command is introduced.

1. The Skill proposes one command at a time and obtains approval for that
   command under the existing workflow contract.
2. For a multi-step analysis, propose deterministic output paths such as:

   ```text
   <workflow-root>/
   ├── rates-across-taxa/<model-or-variant>/
   ├── rates-across-sites/<fraction-or-model>/
   ├── heterotachy/<ghost-variant>/
   ├── compositions-across-taxa/<variant>/
   ├── compositions-across-sites/<model-or-input>/
   └── substitution-patterns/<model-or-simulation>/
   ```

   The exact paths remain user-reviewed parameters; no directory is created
   merely by loading the Skill.
3. Each atomic command's `result.json` remains the execution and provenance
   source of truth. Upstream files must be taken from completed result records
   or explicitly supplied by the user.
4. Recover by inspecting `check_status(output_dir)`, `read_result(output_dir)`,
   and, when the user requested aggregation, `read_report(run_dir)`.
5. `report.json` aggregates completed atomic analyses only. It is not real-time
   workflow state and does not replace status checks.

---

## 6. Documentation and Report Boundaries

### 6.1 Required documentation work

After this design is approved:

- update the parent design and its phase table;
- add the concise trigger/selection section to `SKILL.md`;
- add `skills/phyloai-workflow/references/syserror-workflow.md`;
- add a README overview linking to the detailed workflow;
- cross-link relevant command documentation (`brlen`, `rate`, `taxcomp`, `cca`,
  `signal lnl`, `modelcompare`, `simulate`, `tree bi readpb`, and IQ-TREE);
- audit CLI help examples and Skill parameter annotations for consistency with
  actual schemas, especially CCA preparation, GHOST raw `-m`, recoding names,
  and model-comparison parameters.

Help and reference text must be detailed enough to explain prerequisites,
output meaning, evidence limits, and the next user decision. No command
interface is changed merely to document this workflow.

### 6.2 Report boundary

No new report templates or report-module behavior are required. Interpretation
of the multi-step workflow is AI-driven from atomic `result.json` outputs and
the Skill reference.

Implementation acceptance still verifies that existing report templates and
output indices for completed atomic commands are not broken. The report must
not claim to persist or resume the live systematic-error decision process.

---

## 7. No New Commands or MCP Tools

The workflow uses existing commands and their generated MCP wrappers:

- `posttree syserror brlen`
- `posttree syserror rate`
- `posttree syserror taxcomp`
- `posttree syserror cca`
- `posttree signal lnl`
- `posttree modelcompare iqtree` / `pb`
- `posttree simulate alisim iqtree` / `transfergaps`
- `posttree simulate adequacy`
- `tree ml iqtree`
- `tree bi pb` / `readpb`
- `pretree concat`

GHOST remains a reviewed raw IQ-TREE model expression supplied through the
existing `tree ml iqtree --tool-args` escape hatch. It does not justify a new
MCP tool or an unreviewed default.

---

## 8. Implementation and Acceptance Checklist

- [ ] Approve this workflow design.
- [ ] Update `2026-06-07-phyloai-design.md` to supersede the standalone
      `phyloai-syserror` Skill plan.
- [ ] Add the concise systematic-error entry to `skills/phyloai-workflow/SKILL.md`.
- [ ] Add `skills/phyloai-workflow/references/syserror-workflow.md`.
- [ ] Update README and relevant command-document cross-links.
- [ ] Verify all documented CLI examples against current `--help` output.
- [ ] Verify `get_command_schema` exposes every referenced parameter.
- [ ] Verify GHOST uses an explicitly reviewed raw `-m` in `--tool-args` and no
      guessed model expression.
- [ ] Verify CCA preparation produces two `site_lnl.csv` files, each containing
      the same ordered Tree1/Tree2 pair.
- [ ] Verify strict posterior prediction is labeled `readpb --mode ppred`, while
      AliSim posterior-mean simulation is labeled plug-in parametric simulation.
- [ ] Verify gapless and gap-transferred AliSim paths remain distinguishable.
- [ ] Verify existing report collection/templates still index completed atomic
      outputs without treating the workflow as live report state.
- [ ] Do not commit without separate user approval.
