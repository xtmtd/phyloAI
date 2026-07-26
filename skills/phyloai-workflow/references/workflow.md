# Workflow Reference

## Phase 0: Doctor

Run before external-tool commands when environment is unknown. Summarize missing tools and installation actions.

## Phase 1: Pretree

- **Supermatrix**: `convert -> align -> trim -> metrics / filter -> concat`, then Phase 2 tree inference on the concatenated supermatrix.
- **Supertree**: `convert -> align -> trim -> metrics / filter -> gene trees`, then `tree msc` for species tree from gene trees.
- `stats` is a utility for inspecting sequences/alignments at any step; it is not part of the mandatory pipeline flow.

Check file counts, skipped records, gap ratios, taxon occupancy, and output paths after each step.

## Phase 2: Tree

- **iqtree + msc** as primary tree inference: `tree ml iqtree` for ML trees (gene trees or supermatrix), `tree msc` for coalescent species tree from gene trees.
- **fasttree** for quick exploration: `tree ml fasttree` for fast gene trees or preliminary supermatrix trees.
- **bi** as optional: `tree bi pb` for Bayesian MCMC chain inference. Use `tree bi bpcomp` and `tree bi tracecomp` for final convergence diagnostics with user-chosen burn-in; use `tree bi readpb` for posterior summaries and predictive checks.
- `cf` for concordance factors on any species tree.

## Phase 3: Posttree

Use `topology` for candidate tree tests. Use `dating hessian` then `dating mcmc` for MCMCtree dating.

## Phase 4: Report

Run `report` for `report.json`, `report.html`, methods text, and session recovery.
