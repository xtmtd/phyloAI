# PhyloAI Tree CF (Concordance Factors) Design Specification

**Date:** 2026-06-21
**Status:** Draft
**Parent spec:** `2026-06-07-phyloai-design.md`, `2026-06-17-phyloai-tree-design.md`
**Reference:** Minh et al. (2024) *Mol Biol Evol* 41(11):msae214 — gCF, sCF (parsimony), sCFl (likelihood), qCF naming conventions

---

## 1. Purpose

`phyloai tree cf` computes concordance factors — branch-support measures that quantify the proportion of gene trees or sites supporting each bipartition in a reference species tree. It implements four measures originating from three publications:

| Measure | Origin | Tool | Description |
|---------|--------|------|-------------|
| gCF | Minh et al. (2020) *Mol Biol Evol* 37(5):1530–1534 | IQ-TREE3 | Gene concordance factor: proportion of input gene trees concordant with each reference-tree bipartition |
| sCF (parsimony) | Minh et al. (2020) *Mol Biol Evol* 37(5):1530–1534 | IQ-TREE3 | Site concordance factor using parsimony: average of random quartets per branch |
| sCFl (likelihood) | Mo et al. (2023) *Syst Biol* 72(3):559–574 | IQ-TREE3 | Site concordance factor using maximum likelihood |
| qCF | Mirarab et al. (2014) *Science* 346(6215):1250463 | wASTRAL | Quartet concordance factor: proportion of quartets supporting each bipartition |

These are exposed through five invocation modes:

| Mode      | Measure       | Tool   |
|-----------|---------------|--------|
| `gcf`       | gCF           | IQ-TREE3 |
| `scf`       | sCF (parsimony) | IQ-TREE3 |
| `scfl`      | sCFl (likelihood) | IQ-TREE3 |
| `gcf+scf`   | gCF + sCF     | IQ-TREE3 |
| `qcf`       | qCF           | wASTRAL |

In summary, four measures of CFs are implemented: gCF (Minh et al. 2020), qCF (Mirarab et al. 2014), and the sCF calculated with parsimony (Minh et al. 2020) and likelihood (Mo et al. 2023).

`cf` is a direct `click.Command` (not a Group) — there are two backends (IQ-TREE3 + wASTRAL) but they serve disjoint CF types controlled by a single `--cf` mode selector. CF computation is one-shot — no `--resume` support.

---

## 2. CLI Surface

```bash
# gCF: gene trees only
phyloai tree cf --cf gcf --ref-tree species.nwk --tree-dir ./genetrees/

# gCF: single gene tree file
phyloai tree cf --cf gcf --ref-tree species.nwk --tree merged.trees

# sCF (parsimony): alignment + reference tree (ideally gCF-annotated)
phyloai tree cf --cf scf --ref-tree gCF.cf.tree --matrix msa.fa

# sCFl (likelihood): alignment + reference tree
phyloai tree cf --cf scfl --ref-tree gCF.cf.tree --matrix msa.fa
phyloai tree cf --cf scfl --ref-tree gCF.cf.tree --matrix msa.fa --model-expr LG+F+R4
phyloai tree cf --cf scfl --ref-tree gCF.cf.tree --matrix msa.fa --partitions msa.best_model.nex

# gCF + sCF combined: all inputs in one IQ-TREE3 call
phyloai tree cf --cf gcf+scf --ref-tree species.nwk --tree-dir ./genetrees/ --matrix msa.fa

# qCF: gene trees + reference tree via wASTRAL
phyloai tree cf --cf qcf --ref-tree species.nwk --tree merged.trees

# Custom prefix
phyloai tree cf --cf gcf --ref-tree species.nwk --tree merged.trees --prefix myCF
```

### Command Hierarchy

```
phyloai tree (click.Group)
├── ml (click.Group)
│   ├── fasttree
│   └── iqtree
├── bi                        # Bayesian inference (PhyloBayes-MPI, direct command)
├── msc
└── cf                 # Direct command (click.Command), mode selector: --cf
```

---

## 3. Parameter Specification

### 3.1 Shared Parameters (from main design §9.2)

| Parameter | Short | Type | Default | Description |
|-----------|-------|------|---------|-------------|
| `--output-dir` | `-o` | Path | `runs/tree/cf` | Output directory. |
| `--threads` | `-t` | int >= 1 | 4 | Maps to IQ-TREE3 `-T` or wASTRAL `-t`. |
| `--quiet` | `-q` | flag | False | Suppress terminal output except errors. |
| `--overwrite` | | flag | False | Delete and recreate output directory. |
| `--dry-run` | | flag | False | Show commands without executing. |
| `-h` / `--help` | — | — | — | Auto-handled by Click via root `CONTEXT_SETTINGS`. |

**No `--resume`** — CF computation is one-shot.

**No `--tool-args`** — the CF subcommand covers all known options; no extra tool flags are needed.

### 3.2 Executable Path Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--iqtree-path` | Path | None | Explicit IQ-TREE3 executable path. Resolved via ToolEnv if None. |
| `--wastral-path` | Path | None | Explicit wASTRAL executable path. Resolved via ToolEnv if None. |

### 3.3 CF-Specific Parameters

#### Mode Selector

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--cf` | `gcf`\|`scf`\|`scfl`\|`gcf+scf`\|`qcf` | (required) | Concordance factor type to compute. Determines which input parameters are required and which backend tool is invoked. |

#### Input

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--ref-tree` | Path | (required) | Reference species tree (NEWICK). For gcf/gcf+scf/qcf: topology constraint. For scf/scfl: fixed tree topology (`-te`). For best sCF/sCFl results, provide a gCF-annotated tree. |
| `--tree` | Path | None | Single gene tree file (NEWICK, one tree per line). Mutually exclusive with `--tree-dir`. Required for gcf, gcf+scf, qcf. |
| `--tree-dir` | Path | None | Directory of gene tree files. Merged into one input file (merged.trees). Mutually exclusive with `--tree`. Required for gcf, gcf+scf, qcf. |
| `--matrix` | Path | None | Multiple sequence alignment. Required for scf, scfl, gcf+scf. Maps to IQ-TREE3 `-s`. |
| `--partitions` | Path | None | Partition file for scfl model reuse (e.g., `*.best_model.nex` from IQ-TREE3). NEXUS or RAxML format. Optional, only valid with `--cf scfl`. Mutually exclusive with `--model-expr`. |

#### Options

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
 | `--scf-quartets` | int >= 1 | 100 | Number of quartets for sCF/sCFl. Maps to `--scf N` or `--scfl N`. Recommend >= 100. |
| `--lpp` | flag | False | Append local posterior probabilities (pp1) to qCF support labels. Only valid with `--cf qcf`. |
| `--model-expr` | str | None | Substitution model for scfl (e.g., `LG+F+R4`). Maps to IQ-TREE3 `-m`. Optional, only valid with `--cf scfl`. Mutually exclusive with `--partitions`. |
| `--prefix` | str | auto | Output file prefix. Maps to IQ-TREE3 `--prefix`. Default auto-derived from `--cf` value (see §3.4). |

### 3.4 Default Prefix per Mode

| `--cf` | Default `--prefix` |
|--------|-------------------|
| `gcf`    | `gCF`               |
| `scf`    | `sCF`               |
| `scfl`   | `sCFl`              |
| `gcf+scf` | `gCFsCF`           |
| `qcf`    | `qCF`               |

When user provides explicit `--prefix X`, that value is used for all modes.

---

## 4. CF Mode to Tool Command Mapping

### 4.1 gCF

```
iqtree3 -t <ref-tree> --gcf <gene-trees> --prefix <P> -T <N>
```

IQ-TREE3 generates: `<P>.cf.stat`, `<P>.cf.branch`, `<P>.cf.tree`, `<P>.cf.tree.nex`, `<P>.log`.

### 4.2 sCF (parsimony-based)

```
iqtree3 -s <matrix> -te <ref-tree> --scf <Q> --prefix <P> -T <N>
```

- `-te`: fixed tree topology (unlike gcf which uses `-t`)
- `--scf <Q>`: number of quartets (default 100)

IQ-TREE3 generates: `<P>.cf.stat`, `<P>.cf.branch`, `<P>.cf.tree`, `<P>.cf.tree.nex`, `<P>.log`.

### 4.3 sCFl (likelihood-based)

```
iqtree3 -s <matrix> -te <ref-tree> --scfl <Q> [ -m <model> | -p <partitions> ] --prefix <P> -T <N>
```

- `--scfl <Q>`: number of quartets (default 100)
- `-m <model>` or `-p <partitions>`: optional, for computational speedup using pre-computed model

When `--model-expr` and `--partitions` are both None, IQ-TREE3 computes the best-fit model internally (slow).

IQ-TREE3 generates: `<P>.cf.stat`, `<P>.cf.branch`, `<P>.cf.tree`, `<P>.cf.tree.nex`, `<P>.log` (and model files `.iqtree`, `.ckp.gz`, `.model.gz` when model is auto-computed).

### 4.4 gCF + sCF (combined)

```
iqtree3 -t <ref-tree> --gcf <gene-trees> -s <matrix> --scf <Q> --prefix <P> -T <N>
```

IQ-TREE3 handles the combined workflow internally: compute gCF first, then use the gCF-annotated tree as input for sCF. This is the only combination that can run in a single tool invocation.

IQ-TREE3 generates: `<P>.cf.stat` (containing both gCF and sCF columns), `<P>.cf.branch`, `<P>.cf.tree`, `<P>.cf.tree.nex`, `<P>.log`.

### 4.5 qCF (wASTRAL)

```
wastral -i <gene-trees> -o wastral.tre -u 2 -c <ref-tree> -C --mode 4 -t <N>
```

Fixed flags (not user-configurable):
- `-u 2`: quartet support (required for qCF extraction)
- `-C`: use constraint tree mode
- `--mode 4`: traditional unweighted ASTRAL (required for qCF)

wASTRAL generates raw `wastral.tre`. PhyloAI post-processes this to map qCF values onto the reference tree topology (§4.6).

Output files: `<P>.cf.tree` (mapped qCF tree), `wastral.tre` (raw wASTRAL output, intermediate).

### 4.6 qCF Value Mapping

After wASTRAL runs, phyloAI maps qCF (q1, and optionally pp1) values from `wastral.tre` onto the reference tree topology:

1. Parse both the reference tree (`--ref-tree`) and the wASTRAL output tree (`wastral.tre`) using Bio.Phylo (root the reference tree at its last internal node for consistency).
2. For each internal node in the reference tree, find the corresponding node in the wASTRAL tree by matching the **bipartition** (not rooted clade). To do this:
   - Compute the leaf set under the reference node: `L`
   - In the wASTRAL tree, find the node whose descendant leaf set equals **either** `L` **or** the complement of `L` (all taxa minus `L`). Use the canonical form `min(L, complement(L))` for comparison.
   - This handles unrooted tree equivalence: a bipartition may appear as either a clade or its complement depending on rooting.
3. Extract the `q1` value (and optionally `pp1`) from the matched wASTRAL node. wASTRAL stores q-values inside single-quoted internal node labels such as `'[...;q1=0.422083;pp1=0.95;...]'`. The parser checks `clade.name`, `clade.comment`, and/or the raw Newick label to reliably extract values. Fallback: if Bio.Phylo strips the label, parse the raw Newick string with regex `q1=([0-9.]+)` and `pp1=([0-9.]+)` on each internal node's label.
4. Keep q1 as a raw decimal in [0,1] (no multiplication/rounding to integer). Format: up to
   4 decimal places, trailing zeros stripped (e.g., `0.422083` → `0.4221`, `0.95` → `0.95`,
   `1.0` → `1`).
5. Append the qCF value(s) to the reference node's existing support value, separated by `/`:
   - Without `--lpp`: `<support>/<q1>` (e.g., `100/90` → `100/90/0.4221`)
   - With `--lpp`: `<support>/<q1>/<pp1>` (e.g., `100/90` → `100/90/0.4221/0.95`)
   - If the reference node has no existing support, the result is just the qCF value(s) (e.g., `0.4221` or `1`).
6. Write the result directly to `<P>.cf.tree` by injecting the annotations into the **raw ref-tree Newick string** (not via Bio.Phylo round-trip). This preserves the original branch-length precision and existing support-label formatting exactly — only the new qCF/pp1 values are appended.

This follows IQ-TREE3's convention of appending new support values after existing ones separated by `/`. No additional package dependencies are required beyond Bio.Phylo and the Python standard library.

---

## 5. Input Validation

### 5.1 Parameter Dependencies by Mode

| Check | gcf | scf | scfl | gcf+scf | qcf |
|-------|-----|-----|------|---------|-----|
| `--ref-tree` exists and is readable | ✓ | ✓ | ✓ | ✓ | ✓ |
| `--tree` xor `--tree-dir` required | ✓ | — | — | ✓ | ✓ |
| `--tree` and `--tree-dir` not both set | ✓ | — | — | ✓ | ✓ |
| `--tree`/`--tree-dir` must NOT be set | — | ✓ | ✓ | — | — |
| `--matrix` required | — | ✓ | ✓ | ✓ | — |
| `--matrix` must NOT be set | ✓ | — | — | — | ✓ |
| `--partitions` only valid for scfl | — | — | ✓ | — | — |
| `--model-expr` only valid for scfl | — | — | ✓ | — | — |
| `--model-expr` and `--partitions` not both set | — | — | ✓ | — | — |
| `--scf-quartets` only valid for scf/scfl/gcf+scf | — | ✓ | ✓ | ✓ | — |
| `--scf-quartets` must NOT be set | ✓ | — | — | — | ✓ |
| `--model` must NOT be set | ✓ | — | — | ✓ | ✓ |
| `--partitions` must NOT be set | ✓ | — | — | ✓ | ✓ |

### 5.2 `--tree` (Single File) Mode

- Validate file exists and is readable
- Content validation deferred to backend tool (IQ-TREE3 or wASTRAL)

### 5.3 `--tree-dir` (Directory) Mode

- Scan directory for gene tree files with extensions: `.nwk`, `.tre`, `.tree`, `.nw`, `.trees`, `.newick`
- Skip directories, empty files, files with unrecognized extensions
- Record skipped files in `data.skipped`
- If zero valid files found: exit code 1
- Read each valid file, concatenate into a single file (one newick tree per line)
- Save merged file to `<output_dir>/merged.trees`
- Pass merged file to the backend tool
- If exactly 1 valid file: emit WARNING suggesting `--tree` may be more appropriate, continue

This is the same `_merge_gene_trees()` pattern shared with `tree msc`.

### 5.4 `--scf-quartets` Validation

- Must be >= 1 (Click `IntRange(1, None)`)
- Below 100: WARNING recommending >= 100, continue

### 5.5 Tool Availability

- gcf/scf/scfl/gcf+scf: IQ-TREE3 must be available. If not found: exit code 3.
- qcf: wASTRAL must be available. If not found: exit code 3.

---

## 6. Output Directory Structure

### 6.1 `--cf gcf` mode

```
runs/tree/cf/
├── result.json
├── gCF.cf.stat
├── gCF.cf.branch
├── gCF.cf.tree
├── gCF.cf.tree.nex
├── gCF.log
└── merged.trees             # if --tree-dir used
```

### 6.2 `--cf scf` mode

```
runs/tree/cf/
├── result.json
├── sCF.cf.stat
├── sCF.cf.branch
├── sCF.cf.tree
├── sCF.cf.tree.nex
└── sCF.log
```

### 6.3 `--cf scfl` mode

```
runs/tree/cf/
├── result.json
├── sCFl.cf.stat
├── sCFl.cf.branch
├── sCFl.cf.tree
├── sCFl.cf.tree.nex
├── sCFl.log
├── sCFl.iqtree              # if model auto-computed
├── sCFl.ckp.gz              # if model auto-computed
└── sCFl.model.gz            # if model auto-computed
```

### 6.4 `--cf gcf+scf` mode

```
runs/tree/cf/
├── result.json
├── gCFsCF.cf.stat
├── gCFsCF.cf.branch
├── gCFsCF.cf.tree
├── gCFsCF.cf.tree.nex
├── gCFsCF.log
└── merged.trees             # if --tree-dir used
```

### 6.5 `--cf qcf` mode

```
runs/tree/cf/
├── result.json
├── qCF.cf.tree              # mapped qCF tree
├── wastral.tre              # raw wASTRAL output (intermediate)
└── merged.trees             # if --tree-dir used
```

### 6.6 Output Directory Conflict Policy

- Default: if output directory exists and is non-empty, exit with code 1
- `--overwrite`: delete and recreate the output directory before running
- No `--resume`: CF computation is one-shot

---

## 7. result.json Schema

### 7.1 gCF / sCF / sCFl Example

```json
{
  "status": "success",
  "command": "phyloai tree cf --cf gcf --ref-tree species.nwk --tree merged.trees -t 5",
  "wall_time": 12.34,
  "tool_versions": {"iqtree3": "2.3.6"},
  "params": {
    "cf": "gcf",
    "ref_tree": "/path/to/species.nwk",
    "tree": "/path/to/merged.trees",
    "tree_dir": null,
    "matrix": null,
    "partitions": null,
    "model": null,
    "scf_quartets": null,
    "prefix": "gCF",
    "output_dir": "runs/tree/cf",
    "threads": 5,
    "overwrite": false,
    "dry_run": false,
    "iqtree_path": null,
    "wastral_path": null
  },
  "key_results": {
    "cf_type": "gcf",
    "cf_stat": "runs/tree/cf/gCF.cf.stat",
    "cf_tree": "runs/tree/cf/gCF.cf.tree",
    "prefix": "gCF"
  },
  "error": null,
  "data": {
    "input_mode": "--tree",
    "input": {"path": "/path/to/merged.trees"},
    "cmd": ["iqtree3", "-t", "/path/to/species.nwk", "--gcf", "/path/to/merged.trees",
            "--prefix", "gCF", "-T", "5"],
    "tool_log": "gCF.log",
    "tool_stderr": "# IQ-TREE3 stderr (single pattern, JSON Output Standard Section 5.2)"
  }
}
```

### 7.2 gCF+sCF Combined Example

```json
{
  "key_results": {
    "cf_type": "gcf+scf",
    "cf_stat": "runs/tree/cf/gCFsCF.cf.stat",
    "cf_tree": "runs/tree/cf/gCFsCF.cf.tree",
    "prefix": "gCFsCF"
  }
}
```

### 7.3 qCF Example

```json
{
  "key_results": {
    "cf_type": "qcf",
    "cf_tree": "runs/tree/cf/qCF.cf.tree",
    "prefix": "qCF"
  },
  "data": {
    "input_mode": "--tree",
    "input": {"path": "/path/to/merged.trees"},
    "cmd": ["wastral", "-i", "/path/to/merged.trees", "-o", "wastral.tre",
            "-u", "2", "-c", "/path/to/species.nwk", "-C", "--mode", "4", "-t", "4"],
    "tool_stderr": "# wASTRAL stderr (single pattern, JSON Output Standard Section 5.2)"
  }
}
```

### 7.4 `--tree-dir` Mode Extras (any applicable mode)

```json
{
  "data": {
    "input_mode": "--tree-dir",
    "input": {
      "path": "runs/tree/cf/merged.trees",
      "n_trees": 800
    },
    "skipped": [],
    "warnings": []
  }
}
```

### 7.5 Key Results Fields

| Field | Description |
|-------|-------------|
| `cf_type` | `--cf` mode value |
| `cf_stat` | Path to `.cf.stat` file (IQ-TREE3 modes only) |
| `cf_tree` | Path to `.cf.tree` annotated tree |
| `prefix` | Resolved prefix value |

`data.tool_stderr` captures the tool's raw stderr (single pattern). For IQ-TREE3 modes, `data.tool_log` references the tool's native `.log` file (e.g., `gCF.log`), which is preserved in the output directory as a tool-native artifact and is NOT part of the PhyloAI log model.

---

## 8. Exit Codes

| Code | Meaning |
|------|---------|
| 0 | Success |
| 1 | User input error (mutual exclusivity, invalid param, output dir exists, no valid inputs) |
| 2 | External tool execution failed (IQ-TREE3 or wASTRAL non-zero exit) |
| 3 | Environment error (IQ-TREE3 or wASTRAL not found) |

---

## 9. Warnings

| Condition | Behavior |
|-----------|----------|
| `--scf-quartets` < 100 | WARN: recommend >= 100, continue |
| `--tree-dir` contains exactly 1 valid gene tree file | WARN: suggest `--tree` mode, continue |
| `--tree-dir` contains non-newick files | WARN per file: "unrecognized file extension", record in `data.skipped` |
| `--tree-dir` contains empty files | Skip silently, record in `data.skipped` |
| `--tree-dir` contains 0 valid files | Exit code 1 (not a warning) |
| `--cf scfl` with auto-model computation | INFO: model computation may be slow; suggest `--model-expr` or `--partitions` for speedup |

---

## 10. Logging

- **IQ-TREE3 modes:** the tool's native `.log` file (e.g., `gCF.log`) is preserved in the output directory. IQ-TREE3 writes its own log containing version info, command, and progress.
- **qCF mode:** wASTRAL stderr is inlined in `result.json` as `data.tool_stderr` (single pattern, JSON Output Standard Section 5.2).
- **All modes:** `result.json` is written alongside tool output files. No separate `cf.log` or `wastral.log` is written — summary info is in `result.json`.
- **Terminal output:** During execution, subprocess stdout/stderr is streamed to the terminal in real-time via `subprocess.Popen` + `select`-based line reading, so users can monitor long-running computations. Output is also captured for `data.tool_stderr`. When `--quiet` is set, only errors are shown.

### 10.1 Path Resolution

To ensure correct behavior when subprocesses run with `cwd=output_dir`, all input file paths are resolved to absolute before being passed to tool commands:
- `--ref-tree`, `--matrix`, `--partitions`: resolved via `Path.resolve()` after existence validation.
- `--tree`: resolved via `Path.resolve()` (validated to exist).
- `merged.trees` (from `--tree-dir`): resolved via `Path.resolve()` (live) or `os.path.abspath()` (dry-run, file may not exist yet).
- `--iqtree-path`, `--wastral-path`: resolved via `Path.resolve()` after validation, then returned as absolute `str`.


---

## 11. Implementation Notes

### 11.1 Files to Create

| File | Purpose |
|------|---------|
| `phyloai/tree/cf.py` | `run_cf()` library function and internal helpers |
| `docs/commands/tree-cf.md` | User-facing command documentation |

### 11.2 Files to Modify

| File | Change |
|------|--------|
| `phyloai/cli/commands/tree.py` | Add `"cf"` to `_TreeGroup.list_commands()`, register `@tree.command("cf", cls=_GroupedHelpCommand)` |
| `docs/superpowers/specs/2026-06-17-phyloai-tree-design.md` | Rename `concordance` → `cf` throughout, update CLI examples and hierarchy |
| `docs/superpowers/specs/2026-06-07-phyloai-design.md` | Rename `concordance.py` → `cf.py`, update output dir `concordance/` → `cf/` |

### 11.3 Test Files to Create

| File | Purpose |
|------|---------|
| `tests/tree/test_cf.py` | Library-level tests for `run_cf()`, command builders, qCF mapper |
| `tests/cli/test_tree_cf.py` | CLI integration tests (help output, validation, dry-run) |

### 11.4 Module Structure

```
phyloai/tree/cf.py
├── run_cf()                  # Main entry point
├── _dispatch_mode()          # Select backend and validate mode-specific params
├── _build_iqtree_cf_cmd()    # Build IQ-TREE3 command for any IQ-TREE3 mode
├── _run_iqtree_cf()          # Execute IQ-TREE3, collect output metadata
├── _build_wastral_qcf_cmd()  # Build wASTRAL command for qCF
├── _run_wastral_qcf()        # Execute wASTRAL, capture stderr to data.tool_stderr, map qCF to ref tree
├── _map_qcf_to_tree()        # Map qCF values from wastral.tre onto ref-tree
└── _merge_gene_trees()       # Merge --tree-dir into merged.trees (shared with msc.py)
```

### 11.5 Key Patterns to Follow

- **CLI layer**: thin wrapper — validates params, resolves tools, delegates to library, writes `result.json`, renders Rich summary
- **Library layer**: `run_cf()` accepts all params, validates preconditions, resolves inputs, builds command, runs subprocess, returns payload dict
- **No checkpoint**: CF is one-shot; no `checkpoint.json`, no `--resume`
- **No `--tool-args`**: the `cf` command covers all known use cases; no extra tool flags needed
- **qCF mapping**: `_map_qcf_to_tree()` uses Bio.Phylo to parse both trees, match nodes by canonical bipartition (leaf set vs complement, taking `min(L, complement(L))`), extract q1 (and optionally pp1 when `--lpp`) from `clade.name`/`clade.comment`/raw Newick label (regex `q1=([0-9.]+)` and `pp1=([0-9.]+)` fallback), keep as raw decimal [0,1] with up to 4 decimal places and trailing zeros stripped, then append to existing support value separated by `/`. No additional packages required.
- **Gene tree merging**: `_merge_gene_trees()` follows the same pattern as `phyloai/tree/msc.py`. Extract to a shared helper if both modules stabilize.

### 11.6 Tool Resolution

**IQ-TREE3** (for gcf/scf/scfl/gcf+scf):
1. If `--iqtree-path` provided: validate exists and executable
2. Otherwise: `ToolEnv.require("iqtree3")` — resolves via bundled path or `shutil.which`
3. If not found: exit code 3

Follows the same `_resolve_iqtree_path()` pattern from `phyloai/tree/ml_iqtree.py`.

**wASTRAL** (for qcf):
1. If `--wastral-path` provided: validate exists and executable
2. Otherwise: `ToolEnv.require("wastral")` — resolves via `shutil.which("wastral")` and `shutil.which("aster")`
3. If not found: exit code 3

Follows the same pattern from `phyloai/tree/msc.py`.

### 11.7 CLI Registration in `tree.py`

```python
# cf is a direct command, not a Group (mode selector: --cf)
@tree.command("cf", cls=_GroupedHelpCommand)
@click.option("--cf", type=click.Choice(["gcf", "scf", "scfl", "gcf+scf", "qcf"]),
              required=True, help="Concordance factor type to compute.")
@click.option("--ref-tree", type=click.Path(exists=True, path_type=Path),
              required=True, help="Reference species tree (NEWICK).")
# ... (all remaining options)
def cf_command(...):
    """Concordance factor computation (gCF, sCF, sCFl, qCF)."""
    ...
```

The `tree` group in `_TreeGroup.list_commands()` should include `"cf"` in its return list.

---

## 12. Acceptance Criteria

### 12.1 CLI Validation
- [ ] `--cf` not provided → help shown with error
- [ ] `--cf scf` without `--matrix` → exit 1
- [ ] `--cf scfl --model-expr X --partitions Y` → exit 1 (mutually exclusive)
- [ ] `--cf gcf --matrix msa.fa` → exit 1 (matrix not valid for gcf)
- [ ] `--cf qcf --matrix msa.fa` → exit 1 (matrix not valid for qcf)
- [ ] `--cf scf --tree merged.trees` → exit 1 (tree not valid for scf)
- [ ] `--cf gcf --model LG` → exit 1 (model not valid for gcf)
- [ ] `--cf qcf --partitions p.nex` → exit 1 (partitions not valid for qcf)
- [ ] `--cf gcf --scf-quartets 200` → exit 1 (scf-quartets not valid for gcf)
- [ ] `--tree` and `--tree-dir` together → exit 1

### 12.2 Command Building
- [ ] `--cf gcf` → correct `iqtree3 -t ... --gcf ... --prefix gCF -T 4`
- [ ] `--cf scf` → correct `iqtree3 -s ... -te ... --scf 100 --prefix sCF -T 4`
- [ ] `--cf scfl --model-expr LG+F+R4` → correct `iqtree3 -s ... -te ... --scfl 100 -m LG+F+R4 --prefix sCFl -T 4`
- [ ] `--cf scfl --partitions p.nex` → correct `iqtree3 -s ... -te ... --scfl 100 -p p.nex --prefix sCFl -T 4`
- [ ] `--cf gcf+scf` → correct combined command
- [ ] `--cf qcf` → correct `wastral -i ... -o wastral.tre -u 2 -c ... -C --mode 4 -t 4`

### 12.3 Default Prefix
- [ ] `--cf gcf` (no `--prefix`) → prefix = `gCF`
- [ ] `--cf scf` (no `--prefix`) → prefix = `sCF`
- [ ] `--cf scfl` (no `--prefix`) → prefix = `sCFl`
- [ ] `--cf gcf+scf` (no `--prefix`) → prefix = `gCFsCF`
- [ ] `--cf qcf` (no `--prefix`) → prefix = `qCF`
- [ ] `--cf gcf --prefix myCF` → prefix = `myCF`

### 12.4 Output
- [ ] `result.json` written with correct schema
- [ ] `.cf.stat`, `.cf.tree`, `.cf.branch`, `.cf.tree.nex` produced (IQ-TREE3 modes)
- [ ] `.log` preserved (IQ-TREE3 modes)
- [ ] `wastral.tre` + `<prefix>.cf.tree` produced (qcf mode); stderr inlined in `data.tool_stderr`
- [ ] `merged.trees` produced (`--tree-dir` mode)
- [ ] `tool_versions` populated correctly per mode
- [ ] qCF tree annotations follow `support/q1` convention

### 12.5 Exit Codes
- [ ] Successful run → exit 0
- [ ] Tool non-zero exit → exit 2
- [ ] Tool not found → exit 3
- [ ] Invalid input → exit 1

### 12.6 Warnings
- [ ] `--scf-quartets 50` → WARN "recommend >= 100", continue
- [ ] `--tree-dir` with 1 valid file → WARN "suggest --tree", continue
- [ ] `--tree-dir` with 0 valid files → exit 1

---

## 13. Relationship to Other Modules

- **`tree ml iqtree`**: Consumes MSA, produces gene trees or species tree. CF input data (species tree, gene trees, MSA) typically comes from `tree ml`.
- **`tree msc`**: Consumes gene trees, produces species tree. qCF is conceptually related to wASTRAL's quartet scoring (both use wASTRAL), but qCF is a branch-support measure while MSC is topology inference.
- **`pretree`**: Produces the aligned MSA consumed by CF (via `--matrix`).
- **`posttree`**: Will perform downstream topology tests, dating, and signal analysis on CF-annotated trees.

---

## 14. Design Updates to Parent Specs

### 14.1 `docs/superpowers/specs/2026-06-17-phyloai-tree-design.md`

1. Rename `concordance` → `cf` throughout
2. CLI example: `phyloai tree cf --cf gcf --ref-tree species.nwk --tree-dir ./genetrees/`
3. Hierarchy: `concordance` → `cf` with mode selector `--cf`
4. Module file: `concordance.py` → `cf.py`
5. Documentation file: `docs/commands/tree-concordance.md` → `docs/commands/tree-cf.md`

### 14.2 `docs/superpowers/specs/2026-06-07-phyloai-design.md`

1. Module listing §3: `concordance.py` → `cf.py`
2. Output directory §6: `runs/tree/concordance/` → `runs/tree/cf/`
3. Add CLI example: `phyloai tree cf --cf gcf --ref-tree ... --tree-dir .../`
