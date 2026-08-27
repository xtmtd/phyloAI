# Systematic Error and Heterogeneity Workflow Reference

Use this reference for systematic-error, heterogeneity, LBA, long-branch,
composition-bias, CCA, heterotachy/GHOST, model-adequacy, or related requests.
It guides selection and interpretation of existing atomic commands. It does not
create a one-click diagnostic battery, make an automatic correction, or replace
the standard parameter-card and approval procedure.

## Core rules

1. **All six paths are optional.** Ask which biological contrast, candidate
topologies, taxa, or data properties motivate the analysis. Do not run every
path by default.
2. **Explain the relevant theory first.** Then state what a positive/negative
result can and cannot show before proposing a command.
3. **Use evidence language precisely.** A screen identifies a pattern;
`consistent with` means it matches a mechanism; `sensitivity-supported` means
a user-approved comparison changes in the predicted direction;
`simulation-supported` means an explicit simulation reproduces the behavior.
Do not say a diagnostic proves LBA, establishes a true topology, or fully
validates a model.
4. **One command per approval.** Obtain its full runtime schema, render every
parameter and default, run `doctor` if an external executable is required, and
wait for explicit approval. Do not automatically launch a downstream tree,
simulation, filtering, recoding, or report command.
5. **Never automate irreversible curation.** Taxon/site removal and recoding
are separate user-approved sensitivity analyses. A p-value, branch length,
Keff bin, or rate cutoff is not an automatic deletion criterion.
6. **Preserve alternatives.** Use separate, deterministic, user-reviewed output
paths such as `<workflow-root>/rates-across-sites/slow-50` or
`<workflow-root>/composition-across-taxa/dayhoff6`. Recover completed atomic
steps from `check_status`, `read_result`, or user-requested `read_report`; do
not rerun a completed output merely because the conversation was interrupted.

## Choose an analysis

| User concern or observation | Relevant optional path |
|---|---|
| Long branches, suspected LBA, model-dependent branch lengths | [Rates across taxa](#1-rates-across-taxa) |
| Rapidly/conservatively evolving sites, +G/+R question, site-rate subsets | [Rates across sites](#2-rates-across-sites) |
| A site's rate may change among lineages | [Heterotachy](#3-heterotachy) |
| GC/AT or AA composition differs among taxa | [Compositions across taxa](#4-compositions-across-taxa) |
| Highly constrained site profiles, Keff, CCA, CAT/PMSF topology behavior | [Compositions across sites](#5-compositions-across-sites) |
| Candidate substitution models, mixture/partition choice, relative fit, PPC | [Substitution patterns across sites](#6-substitution-patterns-across-sites) |

ILS, HGT, introgression, recombination, natural selection, gene duplication/loss,
saturation, alignment/orthology errors, and processing errors can also explain
conflict.
They are adjacent but separate diagnostic domains. Do not silently attribute all
conflict to substitution systematic error. Missing data is considered here only
as an optional controlled simulation factor.

## 1. Rates across taxa

### Theory

Branch length is expected substitutions per site. A long branch can result from
more elapsed time, a higher substitution rate, or model-dependent estimation;
it does not alone prove long-branch attraction (LBA). Terminal branches describe
individual taxa; internal branches describe clade separation. Map-defined
node-to-tip and node-to-node distances can represent tip-to-stem,
tip-to-clade-root, or stem-to-stem quantities. Pairwise patristic distances are
not part of the standard rates-across-taxa workflow.

### Screen and compare

1. Ask the user to define focal taxa/clades and the biological branch-length
   quantity before choosing a mode.
2. Use `phyloai posttree syserror brlen` with relevant `terminal`, `internal`,
   `node-to-tip`, or `node-to-node` modes.
3. When comparing topologies/taxa sets, use a `--map` file. Labels created by
   `label-nodes` are only portable to an inspected stable reference topology.
4. Compare the same quantity under the user's chosen homogeneous and
   heterogeneous analyses. Examples include CAT-LG and CAT-GTR Bayesian
   analyses and CAT-PMSF-style ML analyses; CAT-PMSF is an example, not the
   only heterogeneous model.
5. Prefer a posterior-tree distribution to a comparison of one ML tree per
   model. For `tree bi pb` posterior trees, remind the user to inspect
   convergence and normally remove an appropriate burn-in before preparing the
   tree input. Thinning and its value remain user decisions. This Skill does
   not filter treelists or choose burn-in/thinning values.

Interpret model-associated branch-length shifts as model-dependent behavior or
LBA risk. They do not by themselves prove that a model is superior or that one
topology is correct.

### Resolution choices

Offer user-approved options only:

- infer with an appropriate heterogeneous model, for example CAT-LG, CAT-GTR,
  or a CAT-PMSF-style approximation;
- perform a taxon-sampling sensitivity analysis with and without biologically
  justified suspected long-branched taxa;
- combine model and taxon-sampling comparisons.

### Optional advanced posterior prediction

Display only after explaining its cost and the user asks to proceed.

**Primary posterior-predictive route:** run a fixed-topology PhyloBayes analysis
under a user-selected homogeneous and/or CAT-LG/CAT-GTR model; use
`phyloai tree bi readpb --mode ppred` to generate replicate MSAs; then use
`phyloai posttree simulate adequacy` once to compare the observed MSA against
the replicate distribution. Adequacy is a comparison of one observed MSA with
a set of replicates, not an adequacy test independently performed on every
replicate. Re-estimate user-selected trees or branch lengths and compare their
distributions.

**Optional AliSim plug-in route:** for a compatible CAT analysis, run
`phyloai tree bi readpb --mode ss,rr,r`, which writes posterior-mean simulation
inputs including `partition.PMSF.nex`; simulate gapless replicates with
`phyloai posttree simulate alisim iqtree`; then optionally use
`phyloai posttree simulate alisim transfergaps` to restore the observed gap
mask. This is posterior-mean plug-in parametric simulation, **not** a strict
posterior predictive check. Comparing gapless and gap-transferred replicates
can test missing-data sensitivity. Do not request free exchangeabilities from a
fixed-LG chain or invent unavailable posterior parameters.

## 2. Rates across sites

### Theory

Sites evolve at different rates because of structural, functional, and selective
constraints. A single-rate analysis can underestimate repeated changes at fast
sites. Gamma (`+G`) and FreeRate (`+R`) models approximate across-site rate
variation, but inferred rates remain conditional on the fitted model and tree.
Slow sites can have less information; fast sites can have more homoplasy;
neither is universally preferable.

### Screen and sensitivity analysis

1. Use rates produced by IQ-TREE `--rate` or PhyloBayes `readpb -r`.
2. Run `phyloai posttree syserror rate` with exactly one rate source. Include
   `--matrix`, `--subset`, and `--fraction` when extracting alignments:

   ```bash
   phyloai posttree syserror rate \
     --iqtree-rate matrix.rate \
     --matrix matrix.fa \
     --subset slow \
     --fraction 0.25,0.5,0.75
   ```

3. Use `--subset fast` when isolating fast-site behavior is relevant. The
   command retains the selected slow **or** fast fraction; it does not imply
   that only fast sites may be removed.
4. Rebuild user-selected subset matrices separately, then compare topology,
   support, focal branch lengths, and likelihood behavior.
5. Separately compare analysis with and without a rate model, for example
   `LG+F` against `LG+F+G4` or `LG+F+R4`.

A topology change is consistent with rate-associated sensitivity but may also
reflect information loss, composition shifts, or rate estimates conditional on
the model. Stable topology means no conspicuous topology change at the tested
fractions; it does not exclude subtler effects on support, branch lengths, or
likelihood.

### Resolution choices

- use an appropriate `+G` or `+R` model;
- retain or exclude a user-justified slow or fast class for sensitivity;
- use biologically defined partitions where appropriate.

Thresholds are exploratory sensitivity choices, not universal optima.

## 3. Heterotachy

### Theory

Heterotachy is rate variation across both sites and lineages/taxa: a site can be
fast in one lineage and slow in another. Standard `+G`/`+R` analyses normally
allow rates to differ among sites while retaining a site's relative rate across
the tree. IQ-TREE GHOST (General Heterogeneous evolution On a Single Topology)
provides a route for site-by-lineage rate variation.

### Guided flow

1. Use this path when both rates-across-taxa and rates-across-sites effects,
   prior biology, or study design make heterotachy plausible.
2. Run GHOST only through `phyloai tree ml iqtree --tool-args` with an exact,
   user-reviewed raw IQ-TREE `-m` expression appropriate to the sequence type
   and installed IQ-TREE version.
3. The raw `-m` overrides PhyloAI's managed model string. Show the full raw
   expression in the parameter card; never guess GHOST syntax, component count,
   state-frequency option, or defaults.
4. Interpret the tree as heterotachy-aware inference, not proof that
   heterotachy caused an alternative topology.

Do not prescribe AU testing or `posttree modelcompare iqtree` as a GHOST
workflow step. Current structured validation/model-comparison contracts do not
provide a dedicated GHOST comparison route. GHOST does not solve across-site
profile heterogeneity or lineage-specific composition heterogeneity. CAT-LG and
CAT-GTR also do not model heterotachy.

## 4. Compositions across taxa

### Theory

Across-taxon composition heterogeneity means taxa/lineages differ in observed
nucleotide or amino-acid composition. A stationary homogeneous model assumes a
shared equilibrium composition and may group unrelated lineages with similar
compositions. Observed differences can also arise from contamination,
annotation problems, coverage/missing data, or genuine lineage-specific
processes.

### Screen and sensitivity analysis

1. Run `phyloai posttree syserror taxcomp --matrix <alignment>`.
2. Inspect `sparse_count_check` before interpreting nominal chi-square
   p-values.
3. Use per-taxon chi-square contributions and squared composition distances to
   prioritize annotation, contamination, coverage, and lineage-composition
   inspection.
4. Treat nominal and Holm-adjusted p-values as exploratory screens, not
   phylogenetically calibrated pass/fail decisions.

### Resolution choices

- inspect data quality for taxa with unusual composition statistics;
- perform a biologically justified taxon-sampling sensitivity analysis;
- recode using `phyloai pretree concat --recoding Dayhoff-6` (AA) or
  `--recoding RY-nucleotide` (NT), then rebuild the tree;
- consider external/planned nonstationary models such as P4 or GFmix when
  explicit lineage-specific composition modeling is required.

Do not present profile-mixture models as a direct correction for
across-taxon composition heterogeneity. Their profiles primarily model
across-site variation and are shared across lineages. P4/GFmix are not
implemented or deeply guided by the current PhyloAI workflow.

## 5. Compositions across sites

### Theory

Alignment columns can prefer different amino acids/nucleotides because of
structural and functional constraints. A global stationary-frequency vector can
underestimate convergent substitutions at strongly constrained sites.
Profile-mixture models address this across-site profile variation. CCA tests
whether site-wise preference between two fixed topologies changes with the
effective number of amino acids (`Keff`) and model analysis.

CCA needs one `.sitefreq` table plus two model-specific `site_lnl.csv` files.
Every likelihood CSV must contain the same ordered Tree1/Tree2 pair and thus
contain `site`, `lnL_Tree1`, and `lnL_Tree2` columns.

### Prepare site frequencies

Use either compatible source:

1. `phyloai tree bi readpb --mode ss` from a PhyloBayes CAT analysis; or
2. a prior IQ-TREE PMSF analysis that generated a matching `.sitefreq`; the
   classic direct IQ-TREE preparation pattern is:

   ```bash
   iqtree -s <alignment> -m LG+C20+F+G -ft <guide_tree> -n 0
   ```

A prior CAT-PMSF or PMSF `.sitefreq` can be reused only when alignment and site
order exactly match the CCA likelihood inputs.

### Prepare likelihood tables and run CCA

Run `phyloai posttree signal lnl` twice using the same ordered candidate trees,
but distinct model analyses:

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
  --output-dir runs/syserror/cca-input/lnl-heterogeneous
```

For CAT-PMSF-style likelihood analysis, use the exact reviewed model/profile
combination required by the existing `signal lnl` schema; do not invent a
generic profile command. Then run:

```bash
phyloai posttree syserror cca \
  --site-freq analysis.sitefreq \
  --site-lnl1 runs/syserror/cca-input/lnl-LG/site_lnl.csv \
  --site-lnl2 runs/syserror/cca-input/lnl-heterogeneous/site_lnl.csv \
  --model1-name LG \
  --model2-name heterogeneous
```

### Interpretation and choices

Low `Keff` indicates strongly constrained sites. Positive Tree2-minus-Tree1
likelihood differences support Tree2 within that model; negative values support
Tree1. A model-associated pattern concentrated in low-Keff bins is consistent
with constrained-site sensitivity. A weaker/flatter pattern under a
heterogeneous model suggests reduced sensitivity, not proof of adequacy or a
correct topology.

Offer profile-mixture inference (CAT-LG, CAT-GTR, C10–C60/PMSF, or
CAT-PMSF-style), a user-defined low-Keff filtering sensitivity analysis, or
both. PhyloAI has no atomic arbitrary-Keff filtering command: position
extraction is manual preprocessing/future work, and cutoffs such as 5, 6, or 7
are exploratory.

## 6. Substitution patterns across sites

### Theory

Substitution-model misspecification extends beyond rate or composition alone.
Partitioning and mixture models capture different across-site behaviors.
Relative model fit asks which candidate is preferable; absolute adequacy asks
whether a model can reproduce selected observed properties. A model may be best
among the candidates and still be inadequate.

This path addresses broader/residual model misspecification. It overlaps with
compositions across sites but uses model comparison and posterior prediction,
not a Keff attribution. Do not describe CAT/C10–C60 as using a separate
exchangeability matrix at each site: they principally mix site-frequency
profiles with shared exchangeabilities within an analysis.

### Relative fit

For an AA IQ-TREE comparison, provide both model families explicitly:

```bash
phyloai posttree modelcompare iqtree \
  --matrix concat.aa.fa \
  --homogeneous-model LG \
  --heterogeneous-model C10,C20 \
  --het-mrate G,R
```

Use `phyloai posttree modelcompare pb` for PhyloBayes LOO-CV/wAIC when the
needed site-log-likelihood inputs are available. Better AIC/AICc/BIC or
LOO-CV/wAIC is a relative preference among nominated candidates, not proof of
absolute adequacy or topology correctness.

### Optional absolute adequacy / simulation

**Primary route:** fit the chosen PhyloBayes model, check the relevant chain
sampling/convergence diagnostics, run `phyloai tree bi readpb --mode ppred`,
and then run `phyloai posttree simulate adequacy` once using the observed MSA
and complete replicate set. Interpret PPA-DIV, PPA-CONV, PPA-VAR, and PPA-COMP
as separate statistic-specific checks, not one global pass/fail score.

**Optional AliSim route:** for a compatible CAT analysis, generate
`partition.PMSF.nex` with `readpb --mode ss,rr,r`, simulate with AliSim, and
optionally transfer gaps. Label this posterior-mean plug-in simulation, not
strict PPC. Gapless versus gap-transferred replicates can reveal missing-data
sensitivity.

### Resolution choices

- choose a better-supported homogeneous, partition, or profile-mixture model;
- use CAT-LG/CAT-GTR in PhyloBayes or C10–C60/PMSF/CAT-PMSF-style inference in
  IQ-TREE when scientifically appropriate;
- if all tested models remain inadequate, report that limitation and interpret
  it with other sensitivity evidence.

## Report and recovery boundary

No new report template or workflow-state database is created. Existing
`result.json` files remain the provenance source for atomic commands; a
user-requested `report.json` aggregates completed atomic records only. It does
not preserve live decision state or replace `check_status`. Do not claim that
this reference repairs lower-level command provenance, posterior-tree filtering,
or tool-specific model support.
