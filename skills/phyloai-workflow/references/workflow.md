# Workflow Reference

## Phase 0: Doctor

Run before external-tool commands when environment is unknown. Summarize missing tools and installation actions.

## Phase 1: Pretree

`convert -> stats -> align -> trim -> metrics -> filter -> concat`

Check file counts, skipped records, gap ratios, taxon occupancy, and output paths.

## Phase 2: Tree

Use `tree ml iqtree` for final ML inference or `tree ml fasttree` for fast exploration. Use `msc` for species tree from gene trees and `cf` for concordance factors.

## Phase 3: Posttree

Use `topology` for candidate tree tests. Use `dating hessian` then `dating mcmc` for MCMCtree dating.

## Phase 4: Report

Run `report` for `report.json`, `report.html`, methods text, and session recovery.
