# phyloai posttree syserror brlen

[English](posttree-syserror-brlen.md) | [中文](posttree-syserror-brlen.zh.md)


## Purpose

Extract branch-length statistics from phylogenetic trees to diagnose
heterogeneity of rates across taxa (branch length heterogeneity) — a major
source of systematic error in phylogenomics. Branch lengths are compared
across trees inferred under different substitution models to assess
long-branch attraction (LBA) and model-dependent branch length estimation.
This command provides atomic branch-length extraction measurements; it does
not itself identify LBA causality.

Pure Python (Bio.Phylo) — no external executable is invoked, so no
`phyloai doctor` check is needed.

## Usage

```bash
phyloai posttree syserror brlen [options]
phyloai posttree syserror brlen label-nodes --tree <tree.nwk> [options]
```

## Inputs

| Option | Type | Required | Default | Description |
|--------|------|----------|---------|-------------|
| `--tree` | Path | mutex | — | Single tree file (Newick). Mutually exclusive with `--tree-dir`. |
| `--tree-dir` | Path | mutex | — | Directory of tree files. Mutually exclusive with `--tree`. |
| `--mode` | str | yes | — | Comma-separated mode list (see Modes). |
| `--map` | Path | no | — | Node-species map file for node identification. |
| `--node1` | str | no | — | First node name (endpoint modes). |
| `--node2` | str | no | — | Second node name (node-to-node mode). |
| `--tip1` | str | no | — | First tip taxon name. |
| `--tip2` | str | no | — | Second tip taxon name (tip-to-tip mode). |
| `-o`, `--output-dir` | Path | no | `runs/posttree/syserror/brlen` | Output directory. |
| `--table-format` | csv\|tsv | no | csv | Output table delimiter. |
| `-t`, `--threads` | int | no | 4 | Parallel workers for batch mode. |
| `--max-rows` | int | no | 5000000 | Safety limit for patristic output rows. 0 = unlimited. |
| `--overwrite` | flag | no | False | Delete and recreate output directory. |
| `--dry-run` | flag | no | False | Validate inputs (including endpoint resolution) without writing files. |
| `-q`, `--quiet` | flag | no | False | Suppress terminal output except errors. |

## Modes

### Batch modes (combinable via commas, no endpoint parameters required)

| Mode | Description | Output columns |
|------|-------------|----------------|
| `total` | Sum of all branch lengths per tree | tree_file, total_branch_length |
| `terminal` | Each terminal (tip) branch | tree_file, taxon, branch_length |
| `internal` | Each non-root internal branch | tree_file, representation, edge_taxa, branch_length |
| `patristic` | All pairwise tip-to-tip distances | tree_file, tip1, tip2, distance |
| `all` | = total + terminal + internal + patristic | All four CSVs |

`internal.csv` uses a stable schema for mixed rooted/unrooted batches:
`representation` is `rooted` or `unrooted`, and `edge_taxa` lists the
descendant taxa for a rooted edge or the **canonical split side** (the smaller
leaf set, ties broken by lexicographic order) for an unrooted edge. The root is
excluded because it has no incoming branch.

### Endpoint modes (one at a time, not combinable with batch modes)

| Mode | Required parameters | Output columns |
|------|--------------------|----------------|
| `tip-to-tip` | `--tip1`, `--tip2` | tree_file, tip1, tip2, distance |
| `node-to-node` | `--node1`, `--node2`, (`--map` or labeled tree) | tree_file, node1, node2, node1_type, node2_type, distance |
| `node-to-tip` | `--node1`, (`--map` or labeled tree; `--tip1` optional) | tree_file, node, node_type, tip, distance |

`node_type` is `internal` or `tip` (a single-taxon map entry resolves to that
tip; the column makes this explicit). Endpoint modes are mutually exclusive
with batch modes and with each other.

## Node Identification

Node-based modes (`node-to-node`, `node-to-tip`) identify internal nodes by:

1. A `--map` file, if given — it always takes priority over labels.
2. Otherwise, internal node labels (e.g. `N1` from `label-nodes`).

Map file format (colon separates node name from the comma-separated species
list; whitespace is stripped; empty lines and lines without `:` are skipped):

```
NodeName:sp1,sp2,sp3
Outgroup:spA,spB
```

`--map` uses the taxa present in each tree. A present subset can resolve when
it is an exact rooted clade or unrooted split; an empty overlap or incompatible
group emits a warning and skips that tree's endpoint calculation. `Nxx` labels
are suitable only for the reference topology that produced them.

- **Rooted trees** (root has exactly 2 children): a map group must equal the
  exact descendant leaf set of the MRCA (monophyletic clade).
- **Unrooted trees** (root has 3+ children): a map group must equal one side of
  an internal split. Both sides are tested; if only the complement (non-intuitive)
  side matches, a warning records this.

`node-to-tip` without `--tip1`:
- With `--map`: one row per taxon in (map group ∩ tree tips).
- Without `--map`, rooted labeled tree: one row per descendant of the labeled
  node.
- Without `--map`, unrooted tree: error — descendant inference is ambiguous.

## Rooted vs Unrooted Representation

Representation is detected structurally from the Newick root child count:
2 children = rooted, 3+ children = unrooted. This follows the convention of
IQ-TREE, RAxML, FastTree, wASTRAL, and gotree. It is a structural heuristic,
**not a proof of biological rooting** — the same unrooted tree can be written
with a bifurcating root, and a unary root or rooted multifurcation is treated
as unrooted. Users requiring rooted semantics (rooted `edge_taxa`, rooted
labeling, descendant-based `node-to-tip`) must supply a bifurcating rooted
Newick representation.

## Outputs

```
runs/posttree/syserror/brlen/
├── result.json
└── tables/
    ├── total.csv
    ├── terminal.csv
    ├── internal.csv
    ├── patristic.csv
    ├── tip_to_tip.csv
    ├── node_to_node.csv
    └── node_to_tip.csv
```

Only tables for requested modes are created; the file extension follows
`--table-format` (`.csv` or `.tsv`). All successful non-dry runs write exactly
one root `result.json` with status, resolved params, key results (tree counts,
modes, per-mode `n_values`/mean/population SD/min/max), warnings, and
`data.output_files`. `data.warnings` records skipped trees, the patristic row
estimate, and endpoint skip reasons. `tool_versions` is `{}` (pure Python).
No `--resume` and no checkpoints: this is a one-shot utility.

`terminal` reports **all** terminal branches per tree; filter its table for a
single taxon rather than requesting a single-tip terminal mode.

## Examples

```bash
# Batch: all branch lengths from a directory of posterior trees (multi-tree files supported)
phyloai posttree syserror brlen --tree-dir ./posterior_trees --mode all --max-rows 0

# Single: terminal and internal branch lengths
phyloai posttree syserror brlen --tree LG.tre --mode terminal,internal

# Single: tip-to-tip distance between two taxa
phyloai posttree syserror brlen --tree LG.tre --mode tip-to-tip --tip1 Neelus_murinus --tip2 Folsomia_candida

# Batch: node-to-node distance across model trees using map
phyloai posttree syserror brlen --tree-dir ./model_trees --mode node-to-node \
    --map nodes.map.txt --node1 Collembola --node2 Outgroup

# Batch: node-to-tip (map-defined taxa) across posterior trees
phyloai posttree syserror brlen --tree-dir ./posterior_trees --mode node-to-tip \
    --map nodes.map.txt --node1 Collembola

# Node-to-tip with specific tip on a rooted labeled tree (no map needed)
phyloai posttree syserror brlen --tree species.labeled.nwk --mode node-to-tip \
    --node1 N5 --tip1 Folsomia_candida

# Node-to-tip all descendants on a rooted labeled tree (no map needed)
phyloai posttree syserror brlen --tree species.labeled.nwk --mode node-to-tip --node1 N5

# Generate labeled tree and map template for a reference tree
phyloai posttree syserror brlen label-nodes --tree species.nwk
```

### label-nodes

Labels internal nodes `N1..Nxx` (preorder) of a single tree and writes:

| File | Description |
|------|-------------|
| `<stem>.labeled.nwk` | Newick with internal node labels |
| `<stem>.map.txt` | Node-species map template for the main command |
| `<stem>.labeled.pdf` | Tree visualization (matplotlib Agg) |
| `result.json` | Standard PhyloAI result |

Rooted trees label every internal node including the root; unrooted trees
exclude the artificial root. Default output directory:
`runs/posttree/syserror/brlen/label_nodes`.

Labels are unpadded (`N1`, `N2`, ... `Nxx`). Labeling replaces each node's
numeric support value (the `confidence` field is cleared), so `labeled.nwk`
contains clean `Nxx:length` tokens with no mixed label/support text, and the
PDF shows only labels. Branch lengths are written losslessly
(`format_branch_length="%r"`) — the parsed float values are preserved without
the default 5-decimal truncation.

## Warnings / Errors

- Exit code **0**: success. **1**: input validation error (no valid trees,
  invalid mode combination, missing required params, patristic row limit
  exceeded, unresolvable endpoint in single mode). No external tools are
  invoked, so exit codes 2/3 do not apply.
- Single `--tree` input: a parse failure is an exit-1 error. In `--tree-dir`
  mode, invalid files/trees are skipped with warnings; failure occurs only when
  no valid tree remains.
- Single endpoint mode: an unresolved tip/node/map endpoint is an exit-1 error.
  Batch endpoint mode: warn and skip only that tree, and count it in
  `n_trees_skipped`.
- Trees with fewer than two tips are skipped with a warning. A tree whose
  branch lengths are all missing is processed as 0.0 with a warning.
- Patristic output is O(n²) per tree. Before writing, the estimated row count
  is checked against `--max-rows` (default 5,000,000; 0 disables the limit),
  and a warning records the estimate.
- `Bio.Phylo` may parse arbitrary text as a degenerate one-tip tree; such input
  emits the fewer-than-two-tips warning rather than a parse-failure warning.

## Notes

- Branch lengths are substitutions per site; they do not independently
  distinguish elapsed time from substitution rate.
- Multi-tree files are identified as `filename:index` (zero-based); a one-tree
  file keeps the bare filename.
- Batch processing shows a transient Rich progress bar ("Processing trees")
  unless `--quiet`/`--dry-run`. Fast runs finish before it is perceptible;
  large posterior batches (especially patristic) show steady progress.
  `--threads` parallelizes `--tree-dir` (non-patristic) processing; a single
  `--tree` file executes serially.
- The `command` recorded in `result.json` is the reproducible invocation: it
  always includes the required inputs and `-o`, plus explicitly provided
  options and non-default values; default-valued flags are omitted. `params`
  always carries the full resolved values.
- This command extracts diagnostic measurements only; it does not declare which
  model is superior. Compare distributions across model runs.
- For systematic-error interpretation, model/taxon sensitivity choices, and
  optional posterior-predictive simulation, see the
  [systematic-error workflow reference](../../skills/phyloai-workflow/references/syserror-workflow.md).
  When a posterior tree distribution is used, inspect convergence and normally
  prepare a post-burn-in tree input before running `brlen`; this command does
  not select burn-in or thin/filter a treelist.
