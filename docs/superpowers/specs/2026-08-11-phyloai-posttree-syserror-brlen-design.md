# PhyloAI `posttree syserror brlen` Design Specification

**Date:** 2026-08-11
**Status:** Implemented (2026-08-11)
**Parent spec:** `2026-06-07-phyloai-design.md` (Section 4.1, Phase 4)

---

## 1. Purpose

Compute branch length statistics from phylogenetic trees to diagnose **heterogeneity of rates across taxa** (branch length heterogeneity) — a major source of systematic error in phylogenomics. Long-branch attraction (LBA) and model-dependent branch length estimation are assessed by comparing branch lengths across trees inferred under different substitution models.

This command provides atomic branch length extraction operations.

---

## 2. Command Structure

```
phyloai posttree syserror brlen [options]
phyloai posttree syserror brlen label-nodes --tree <tree.nwk> [options]
```

CLI hierarchy addition:

```
posttree
└── syserror (new click.Group)
    └── brlen (click.Group)
        ├── (default command — branch length calculation)
        └── label-nodes (sub-command)
```

`_PosttreeGroup.list_commands` adds `"syserror"`.

MCP tools auto-generated from Click tree:
- `posttree_syserror_brlen` — main branch length calculation
- `posttree_syserror_brlen_label_nodes` — node labeling helper

Upon implementation, remove `posttree_syserror_brlen` from `phyloai/mcp/tools/stubs.py`.

---

## 3. Main Command: Branch Length Calculation

### 3.1 Input Parameters

| Parameter | Type | Required | Default | Description |
|-----------|------|----------|---------|-------------|
| `--tree` | Path | mutex | — | Single tree file (Newick). Mutually exclusive with `--tree-dir`. |
| `--tree-dir` | Path | mutex | — | Directory of tree files. Mutually exclusive with `--tree`. |
| `--mode` | str | yes | — | Comma-separated mode list (see 3.2). |
| `--map` | Path | no | — | Node-species map file for node identification. |
| `--node1` | str | no | — | First node name (endpoint modes). |
| `--node2` | str | no | — | Second node name (node-to-node mode). |
| `--tip1` | str | no | — | First tip taxon name. |
| `--tip2` | str | no | — | Second tip taxon name (tip-to-tip mode). |
| `--output-dir` / `-o` | Path | no | `runs/posttree/syserror/brlen` | Output directory. |
| `--table-format` | csv\|tsv | no | csv | Output table delimiter. |
| `--threads` / `-t` | int | no | 4 | Parallel workers for batch mode. |
| `--max-rows` | int | no | 5000000 | Safety limit for patristic output rows. 0 = unlimited. |
| `--overwrite` | flag | no | False | Delete and recreate output directory. |
| `--dry-run` | flag | no | False | Validate inputs (including endpoint resolution) without writing files. |
| `--quiet` / `-q` | flag | no | False | Suppress terminal output except errors. |

### 3.2 Mode Classification

#### Batch Modes (combinable, no endpoint parameters required)

| Mode | Description | Output CSV columns |
|------|-------------|-------------------|
| `total` | Sum of all branch lengths per tree | tree_file, total_branch_length |
| `terminal` | Each terminal (tip) branch | tree_file, taxon, branch_length |
| `internal` | Each non-root internal branch | tree_file, representation, edge_taxa, branch_length |
| `patristic` | All pairwise tip-to-tip distances | tree_file, tip1, tip2, distance |
| `all` | = total + terminal + internal + patristic | All four CSVs |

**`internal` mode column behavior:**
- Every row uses `representation` (`rooted` or `unrooted`) and `edge_taxa`, so a mixed rooted/unrooted batch has one stable CSV schema.
- **Rooted tree**: `edge_taxa` lists all descendant taxa of the internal node (unambiguous direction).
- **Unrooted tree**: `edge_taxa` lists the **canonical side** of the bipartition defined by that internal edge. Canonical side = the smaller leaf set; ties broken by lexicographic order of the first taxon. This ensures stable identifiers across trees with differing root placements or Newick representations, enabling cross-tree matching of the same biological edge.
- The root is excluded because it has no incoming branch.

#### Endpoint Modes (one at a time, not combinable with batch modes)

| Mode | Required parameters | Output CSV columns |
|------|--------------------|--------------------|
| `tip-to-tip` | --tip1, --tip2 | tree_file, tip1, tip2, distance |
| `node-to-node` | --node1, --node2, (--map or labeled tree) | tree_file, node1, node2, node1_type, node2_type, distance |
| `node-to-tip` | --node1, [--tip1 optional], (--map or labeled tree; see 3.6) | tree_file, node, node_type, tip, distance |

`node_type` column: `internal` or `tip` (when a map entry defines only one taxon, the "node" is actually a tip; this column makes the distinction explicit to avoid misinterpretation).

#### Constraint Rules

- `all` is mutually exclusive with other batch modes.
- Endpoint modes cannot be combined with batch modes.
- Only one endpoint mode per invocation.
- Batch modes can be freely combined (e.g., `--mode terminal,internal,total`).

### 3.3 Rooted vs Unrooted Representation Detection

PhyloAI determines tree representation from the Newick structure:

- **Rooted representation**: `tree.root` has exactly 2 children (bifurcating root).
- **Unrooted representation**: `tree.root` has 3+ children (trifurcating/multifurcating root).

This follows the standard Newick convention used by IQ-TREE, RAxML, FastTree, wASTRAL, gotree, and other mainstream phylogenetic tools: unrooted binary trees are written with a trifurcating root; explicitly rooted trees (via outgroup rooting, midpoint rooting, or time-calibration) have a bifurcating root.

Detection is automatic; there is no `--rooted` flag.

**User responsibility caveat**: PhyloAI cannot determine biological rootedness from Newick topology alone — the same unrooted tree can technically be written with a bifurcating root. Users must ensure their input tree is rooted as intended (e.g., via outgroup rooting) before relying on rooted-representation behaviors. A unary root or rooted multifurcation is also classified as unrooted by this structural heuristic; users requiring rooted behavior must supply a bifurcating rooted Newick representation. Help text states these limits explicitly.

The rooted/unrooted representation status affects:
- `internal` mode `edge_taxa` interpretation (Section 3.2)
- label-nodes: whether the root node is labeled (Section 4)
- node-to-tip without --tip1: whether descendant inference is allowed without --map (Section 3.7)

### 3.4 Node Identification Logic

For modes requiring node identification (node-to-node, node-to-tip):

1. Check if tree internal nodes already have labels (e.g., N1, N2... from `label-nodes`).
2. If a label matches `--node1`/`--node2`, use the labeled node directly.
3. If no matching label exists, `--map` must be provided.
4. If both --map and labeled tree are present, **--map takes priority** (overrides tree labels).

### 3.5 Map File Format

```
NodeName:sp1,sp2,sp3
Outgroup:spA,spB
```

Parsing rules:
- Colon separates node name from species list.
- Commas separate species names.
- Leading/trailing whitespace around all tokens is stripped (tolerant parsing).
- Both `Node1:sp1,sp2` and `Node1: sp1, sp2, sp3` are valid.
- Empty lines and lines without `:` are skipped.
- Single-taxon entries are allowed (e.g., `OutgroupTip:spA`); the "node" resolves to that tip.

### 3.6 Monophyly Validation — Rooted Clade vs Unrooted Split Matching

Node identification semantics differ by tree rootedness:

#### Rooted trees (bifurcating root)

Standard clade-based matching. For each map-defined node with species set S:

1. Compute the intersection T = S ∩ {tree tips}.
2. If T is empty → warn and skip.
3. If |T| == 1 → the "node" is that single tip (node_type = `tip`).
4. If |T| >= 2 → find the MRCA of T. Check whether the set of leaves under MRCA equals T.
   - **If yes** → monophyletic. Use MRCA.
   - **If no** → not monophyletic. Warn and skip this node for this tree.

#### Unrooted trees (trifurcating root)

Map semantics operate on **unrooted bipartitions** (splits). A map-defined set S matches if it corresponds to either side of any internal edge in the tree.

Algorithm (edge-traversal split matching):

1. Compute the intersection T = S ∩ {tree tips}.
2. If T is empty → warn and skip.
3. If |T| == 1 → the "node" is that single tip (node_type = `tip`).
4. If |T| >= 2 → enumerate all internal edges of the tree. Each edge defines a bipartition: removing the edge splits the tree into two leaf sets (side_A, side_B).
   - If T == side_A or T == side_B → match found. The matched edge's internal endpoint (the node on the T side) is the target node.
   - If no edge matches → S is not a valid split on this tree. Warn and skip this node for this tree.

This approach is independent of root placement and correctly identifies splits regardless of Newick representation or ladderization. For node-to-tip distance computation, the matched edge endpoint provides the unambiguous starting node for distance calculation.

**Implementation note:** Edge enumeration is O(n) per tree (one pass over all clades). For each internal node, the leaf set "below" it vs. the complement forms one bipartition. This is equivalent to traversing all internal edges without explicitly removing edges.

The `data.warnings` records which side matched when the complement (non-intuitive) side was used.

### 3.7 node-to-tip Behavior Without --tip1

| --map | --tip1 | labeled tree | rooted tree | Behavior |
|-------|--------|--------------|-------------|----------|
| yes | yes | any | any | Compute node → tip1 distance |
| yes | no | any | any | Compute node → each tip in (map species set ∩ tree tips) |
| no | yes | yes | any | Compute node → tip1 distance (node by label) |
| no | no | yes | **yes** | Compute node → all descendants of labeled node |
| no | no | yes | **no** | **Error**: `--tip1 or --map required for unrooted trees (descendant ambiguity)` |
| no | any | no | any | **Error**: `--map or labeled tree required for node identification` |

When --map is provided and --tip1 is omitted, the output contains one row per tip in (map species set ∩ tree tips) — not all descendants of the MRCA (which may include extra taxa in subset-match scenarios).

### 3.8 Tree File Scanning and Validation

**Directory scanning** (`--tree-dir`):
- Scans all non-empty regular files in the directory (non-recursive).
- No extension filtering — posterior tree files (.treelist, .trees, etc.) have no standard suffix.
- Each file is attempted with `Bio.Phylo.parse(path, "newick")` to support **multi-tree files** (e.g., PhyloBayes posterior treelists).
- Files that fail parsing entirely are skipped with a warning.
- Sorting: alphabetical by filename for deterministic output order.

**Multi-tree files**:
- `Bio.Phylo.parse()` yields multiple trees from a single file.
- In CSV output, multi-tree files are identified as `filename:index` (0-based), e.g., `posterior.treelist:0`, `posterior.treelist:1`, etc.
- Single-tree files use just the filename (no `:0` suffix).

**Per-tree validation**:
1. Tree has >= 2 terminals (otherwise skip with warning).
2. If all branch lengths are None → warn (calculations yield 0.0).

**Failure behavior**:
- **Single mode** (`--tree`): parse failure → error exit (code 1).
- **Batch mode** (`--tree-dir`): skip invalid files/trees; if no valid trees remain → exit code 1.
- **Endpoint modes**: in single mode, a missing tip/node/map match or incompatible clade/split → error exit (code 1). In batch mode, warn and skip only the affected tree.
- **Skip counting**: `n_trees_skipped` counts skipped processing units. Each degenerate parsed tree counts as one; a file that fails parsing or that is non-empty but yields zero Newick trees counts as one; in batch endpoint mode, each tree whose tip/node/map endpoint cannot be resolved also counts as one. Each skipped unit has a warning.
- **Parser tolerance**: `Bio.Phylo` may parse arbitrary text as a degenerate one-tip tree. Such input emits the fewer-than-two-tips warning rather than a parse-failure warning; structurally invalid Newick emits the parse-failure warning.
- **`--dry-run`** runs the same endpoint resolution as a real run (so an unresolvable single-tree endpoint fails preflight with exit 1 and batch endpoint skips are reported), but writes no files.

### 3.9 Patristic Mode Safety Limit

Patristic distance is O(n²) per tree. For large datasets (e.g., 100 taxa × 1000 posterior trees = ~5M rows), unbounded output can be problematic.

- Before computing, estimate total rows: `sum(n_tips*(n_tips-1)/2 for each tree)`.
- If estimated rows > `--max-rows` (default 5,000,000) and `--max-rows != 0`: exit with error message showing the estimate and suggesting `--max-rows 0` to force.
- Help text explicitly warns about O(n²) scaling and the structural rootedness heuristic.
- `result.json` `data.warnings` records the row count estimate.
- Patristic rows are streamed directly to their CSV/TSV and summarized online; they are never accumulated as in-memory row dictionaries.

### 3.10 Batch Mode Parallelism

- `--threads` controls `concurrent.futures.ProcessPoolExecutor` worker count for non-patristic calculations.
- Each tree file (including a multi-tree file, to avoid excessive task granularity) is one parallel unit.
- A transient Rich progress bar ("Processing trees") is shown when not `--quiet`/`--dry-run`: it advances per file for batch non-patristic modes and per tree for patristic streaming. Fast runs finish before the bar is perceptible; large posterior batches show steady progress.
- Requested patristic output is written serially in deterministic tree order to preserve streaming memory bounds; in `--mode all`, only the patristic stage is serial.
- Results merged into unified CSVs with `tree_file` column identifying source.
- Default: 4 threads.

---

## 4. Sub-command: `label-nodes`

```bash
phyloai posttree syserror brlen label-nodes --tree tree.nwk [-o output_dir] [--overwrite] [--quiet]
```

### 4.1 Behavior

- **Single mode only** (one tree). No `--tree-dir` support.
- Assigns sequential labels N1, N2, ..., Nxx to internal nodes.
- Labels replace the node's numeric support value: the `confidence` attribute of
  each labeled node is cleared so `labeled.nwk` contains only `Nxx:length` and
  the PDF shows only labels (no mixed label/support tokens).
- Branch lengths are written losslessly (`format_branch_length="%r"`), so
  `labeled.nwk` preserves the parsed float values without the default 5-decimal
  truncation.
- **Rooted tree** (bifurcating root): all internal nodes including root are labeled.
- **Unrooted tree** (trifurcating root): the artificial root node created by Bio.Phylo is **excluded** from labeling (it has no biological meaning and its "descendants" depend on arbitrary Newick representation). Only true internal nodes below root are labeled.
- Numbering follows preorder traversal (excluding artificial root for unrooted trees).
- **Intended for a single reference tree** (e.g., species tree). Nxx labels are NOT stable across different topologies and should NOT be reused across trees. For batch analysis of trees with differing topologies, always use map files with split-based node definitions.

### 4.2 Parameters

| Parameter | Type | Required | Default | Description |
|-----------|------|----------|---------|-------------|
| `--tree` | Path | yes | — | Single Newick tree file. |
| `--output-dir` / `-o` | Path | no | `runs/posttree/syserror/brlen/label_nodes` | Output directory. |
| `--overwrite` | flag | no | False | Overwrite output. |
| `--quiet` / `-q` | flag | no | False | Suppress output. |

### 4.3 Outputs

All files placed directly under `<output-dir>`:

| File | Description |
|------|-------------|
| `<stem>.labeled.nwk` | Newick with internal nodes labeled N1..Nxx |
| `<stem>.map.txt` | Map file: each node and its descendant taxa |
| `<stem>.labeled.pdf` | Tree visualization via `Bio.Phylo.draw()` + matplotlib |
| `result.json` | Standard PhyloAI result (at `<output-dir>/result.json`) |

Map file format (generated without spaces for compactness):
```
N1:sp1,sp2,sp3,sp4,sp5
N2:sp1,sp2,sp3
N3:sp1,sp2
N4:sp4,sp5
```

For rooted trees, "descendant taxa" is unambiguous. For unrooted trees (artificial root excluded), each labeled node's descendants are defined by the parse-tree structure — the map file serves as a starting template that users should review and edit before use with the main command.

PDF generated via:
```python
import matplotlib.pyplot as plt
from Bio import Phylo

fig, ax = plt.subplots(figsize=(12, 8))
Phylo.draw(tree, axes=ax, do_show=False)
fig.savefig(output_path, format="pdf", bbox_inches="tight")
plt.close(fig)
```

---

## 5. Output Directory Structure

### Main command:
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

Only files for computed modes are present.

### label-nodes:
```
runs/posttree/syserror/brlen/label_nodes/
├── result.json
├── <stem>.labeled.nwk
├── <stem>.map.txt
└── <stem>.labeled.pdf
```

---

## 6. result.json Schema

### 6.1 Main command (success):

```json
{
  "status": "success",
  "command": "phyloai posttree syserror brlen --tree-dir ./trees --mode terminal,internal -o runs/posttree/syserror/brlen",
  "wall_time": 2.3,
  "tool_versions": {},
  "params": {
    "tree": null,
    "tree_dir": "./trees",
    "mode": "terminal,internal",
    "map": null,
    "node1": null,
    "node2": null,
    "tip1": null,
    "tip2": null,
    "table_format": "csv",
    "threads": 4,
    "max_rows": 5000000,
    "output_dir": "runs/posttree/syserror/brlen",
    "overwrite": false,
    "dry_run": false,
    "quiet": false
  },
  "key_results": {
    "n_trees": 50,
    "n_trees_skipped": 2,
    "modes": ["terminal", "internal"],
    "summary": {
      "terminal": {
        "n_values": 2400,
        "mean": 0.0523,
        "sd": 0.0312,
        "min": 0.0001,
        "max": 0.4512
      },
      "internal": {
        "n_values": 2350,
        "mean": 0.0089,
        "sd": 0.0045,
        "min": 0.0001,
        "max": 0.0321
      }
    }
  },
  "error": null,
  "data": {
    "summary": {
      "n_trees_processed": 50,
      "n_trees_skipped": 2,
      "n_multi_tree_files": 3
    },
    "warnings": [
      "tree file gene42.tre: failed to parse as Newick, skipped",
      "tree gene17.tre: node 'Collembola' is not monophyletic in this tree"
    ],
    "output_files": {
      "terminal_table": {
        "path": "/abs/path/runs/posttree/syserror/brlen/tables/terminal.csv",
        "description": "Terminal branch lengths per taxon per tree"
      },
      "internal_table": {
        "path": "/abs/path/runs/posttree/syserror/brlen/tables/internal.csv",
        "description": "Internal branch lengths per tree (edge_taxa with rooted/unrooted representation)"
      }
    }
  }
}
```

The `command` field is the reproducible invocation: it always records the
required inputs (`--tree`/`--tree-dir`, `--mode`) and `-o`, plus any explicitly
provided optional flags (`--map`, `--node1`...`--tip2`) and non-default values
(`--table-format`, `-t`, `--max-rows`) or set flags (`--overwrite`, `--dry-run`,
`-q`). Default-valued flags are omitted; `params` always carries the full
resolved values for verification.

### 6.2 label-nodes (success):

```json
{
  "status": "success",
  "command": "phyloai posttree syserror brlen label-nodes --tree species.nwk -o runs/posttree/syserror/brlen/label_nodes",
  "wall_time": 0.5,
  "tool_versions": {},
  "params": {
    "tree": "species.nwk",
    "output_dir": "runs/posttree/syserror/brlen/label_nodes",
    "overwrite": false,
    "quiet": false
  },
  "key_results": {
    "n_internal_nodes_labeled": 47,
    "n_terminals": 48,
    "rooted": true
  },
  "error": null,
  "data": {
    "output_files": {
      "labeled_tree": {
        "path": "/abs/path/runs/posttree/syserror/brlen/label_nodes/species.labeled.nwk",
        "description": "Newick tree with internal node labels N1..N47"
      },
      "map_file": {
        "path": "/abs/path/runs/posttree/syserror/brlen/label_nodes/species.map.txt",
        "description": "Node-species map template for brlen node modes"
      },
      "tree_figure": {
        "path": "/abs/path/runs/posttree/syserror/brlen/label_nodes/species.labeled.pdf",
        "description": "Tree visualization with labeled internal nodes"
      }
    }
  }
}
```

**JSON Output Standard compliance notes:**
- This command is a **pure-Python utility** (no external tool invocations). Per JSON Output Standard Section 6.2, it does NOT use `data.cmd`, `data.tool_stderr`, `data.files[]`, or `data.tool_log`.
- `data.output_files` uses the `{label: {path, description}}` dict format per Section 5.4.
- `data.summary` provides aggregate statistics per Section 6.2 (pure-Python batch pattern).

---

## 7. Report Template

Update `phyloai/report/templates.py`. The report template covers the **main command only**; `label-nodes` is an auxiliary helper that does not produce analytical conclusions and is excluded from methods text generation.

```python
_MODE_MEAN_PHRASES = {
    "total": "Mean total branch length per tree",
    "terminal": "Mean terminal branch length",
    "internal": "Mean internal branch length",
    "patristic": "Mean pairwise tip-to-tip distance",
    "tip-to-tip": "Mean tip-to-tip distance",
    "node-to-node": "Mean node-to-node distance",
    "node-to-tip": "Mean node-to-tip distance",
}


def generate_methods_posttree_syserror_brlen(
    params: dict[str, Any],
    key_results: dict[str, Any],
    tool_versions: dict[str, Any],
) -> str:
    n_trees = key_results.get("n_trees", 0)
    modes = key_results.get("modes", [])
    modes_str = ", ".join(modes)
    has_map = bool(params.get("map"))

    text = (
        f"Branch length heterogeneity was assessed by extracting "
        f"{modes_str} branch lengths from "
        f"{_describe_n(n_trees, 'phylogenetic tree', 'phylogenetic trees')}"
    )
    if has_map:
        text += " with internal nodes identified via a monophyletic group map file"
    text += " using PhyloAI (Bio.Phylo)."

    summary = key_results.get("summary", {})
    for mode in modes:
        phrase = _MODE_MEAN_PHRASES.get(mode)
        if phrase is None or mode not in summary:
            continue
        s = summary[mode]
        if s.get("mean") is None:
            continue
        text += (
            f" {phrase} was {_safe_fmt(s.get('mean'), '.4f')}"
            f" (SD = {_safe_fmt(s.get('sd'), '.4f')})."
        )

    n_skipped = key_results.get("n_trees_skipped", 0)
    if n_skipped > 0:
        text += (
            f" {_describe_n(n_skipped, 'tree', 'trees')} were skipped due to "
            f"parsing or monophyly validation failures."
        )

    return text
```

---

## 8. Dependencies

- **Bio.Phylo** (Biopython, existing) — tree parsing, distance computation, MRCA, traversal
- **matplotlib** (existing) — PDF tree visualization in label-nodes
- **concurrent.futures** (stdlib) — batch parallelism
- **itertools.combinations** (stdlib) — patristic pairwise enumeration
- **No new dependencies**

Custom implementation required:
- Rooted/unrooted detection (Section 3.3): root children count (2 = rooted, 3+ = unrooted)
- Rooted clade matching (Section 3.6): MRCA leaf set == map set intersection
- Unrooted split matching (Section 3.6): edge-traversal bipartition enumeration, match T against each edge's two sides
- Canonical split representation (Section 3.2): smaller leaf set of bipartition; ties by lexicographic order

---

## 9. Integration Points

| Component | Change |
|-----------|--------|
| `phyloai/posttree/syserror_brlen.py` | New file: core computation logic |
| `phyloai/cli/commands/posttree.py` | Add `syserror` group → `brlen` group + `label-nodes` sub-command |
| `phyloai/report/templates.py` | Update `generate_methods_posttree_syserror_brlen` (main cmd only) |
| `phyloai/mcp/tools/stubs.py` | Remove `posttree_syserror_brlen` from stubs |
| MCP schema | Auto-generates `posttree_syserror_brlen` + `posttree_syserror_brlen_label_nodes` from Click |
| `docs/commands/posttree-syserror-brlen.md` | New command documentation |
| `README.md` | Add syserror brlen to command list |
| `README.zh.md` | Add syserror brlen to command list (Chinese) |
| `docs/commands/ai-integration.md` | Replace the brlen stub entry with generated brlen and label-nodes MCP tools |
| `docs/commands/ai-integration.zh.md` | Replace the brlen stub entry with generated brlen and label-nodes MCP tools (Chinese) |
| Main design doc | No change needed (already lists `syserror brlen`) |

---

## 10. Exit Codes

| Code | Condition |
|------|-----------|
| 0 | Success |
| 1 | Input validation error (no valid trees, invalid mode combination, missing required params, patristic row limit exceeded) |

No external tools are invoked, so exit code 2 (tool failure) and 3 (environment error) do not apply.

---

## 11. Examples

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

---

## 12. Help Text

The `--help` output must be detailed and self-contained, covering:
- All modes with descriptions
- Mode classification (batch vs endpoint) and mutual exclusivity rules
- Node identification methods (map file vs labeled tree) with priority rules
- Rooted vs unrooted tree behavior differences
- Map file format with inline example
- Patristic O(n²) warning and --max-rows explanation
- Multi-tree file support (Bio.Phylo.parse)
- Input/output descriptions
- Workflow examples showing typical syserror diagnosis patterns

---

## 13. Terminology

- **node-to-tip**: always this form; never "tip-to-node" or "tip-to-stem" in code/output. Documentation may note that "tip-to-stem" in literature corresponds to `node-to-tip` mode with the stem node defined via map.
- **tip-to-tip** (endpoint mode): compute distance between a specified pair of tips. Distinct from **patristic** (batch mode) which computes ALL pairwise tip distances.
- **split matching**: the unrooted bipartition semantics used for map validation (Section 3.6).

---

## 14. Future Considerations

- **Comparison plots**: The `phyloai-syserror` Skill (Phase 9) will orchestrate branch length comparison across models (e.g., LG vs CAT-PMSF posterior trees) and generate boxplot visualizations. This is workflow logic, not part of the atomic `brlen` command.
- **Statistical tests**: Wilcoxon rank-sum tests comparing branch length distributions between groups/models will be a Skill-level operation reading the CSV outputs.
- **Root-to-tip mode**: For rooted trees, a dedicated `root-to-tip` batch mode computing distances from root to every terminal could be useful for clock-likeness assessment. Deferred until needed.
