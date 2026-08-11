# PhyloAI Posttree Syserror Brlen Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Implement `phyloai posttree syserror brlen` and `label-nodes` to extract branch-length statistics from one tree, a tree directory, or Newick multi-tree files without new dependencies.

**Architecture:** Put all parsing, topology-aware node resolution, calculations, CSV writing, and `result.json` construction in one pure-Python library module. Keep the Click group as a thin wrapper. A labeled reference tree is convenient for one stable topology; map files are the required durable node-identification interface for batch trees whose topology or taxon set can differ.

**Tech Stack:** Python 3.10+, Biopython `Bio.Phylo`, matplotlib, Click, Rich, pytest, and stdlib (`csv`, `itertools`, `concurrent.futures`, `json`, `shlex`, `statistics`).

## Global Constraints

- Implement the approved design in `docs/superpowers/specs/2026-08-11-phyloai-posttree-syserror-brlen-design.md`; do not add features outside it.
- Add no pip dependency. Use existing Biopython and matplotlib only.
- Accept exactly one of `--tree` or `--tree-dir`; scan a directory non-recursively, alphabetically, with no suffix filter, and parse every file with `Bio.Phylo.parse(..., "newick")`.
- A single input parse failure is an exit-code-1 error. In directory mode, skip invalid files/trees with warnings and fail only when no valid tree remains.
- Treat `None` branch lengths as `0.0`. Validate at least two terminals. Support multiple Newick trees per file as `filename:index` (zero-based); a one-tree file remains `filename`.
- Add a per-tree warning when every branch length is `None`, and record the estimated total patristic row count in `data.warnings` whenever patristic output is requested.
- Keep `--mode terminal` as the sole terminal-branch operation. Do not add a redundant single-tip terminal mode; users can filter `terminal.csv`.
- Batch modes are `total`, `terminal`, `internal`, `patristic`, and `all`; endpoint modes are `tip-to-tip`, `node-to-node`, and `node-to-tip`. Enforce every mutual-exclusion and required-parameter rule from spec section 3.2 before output creation.
- Detect representation only by root child count: 2 is rooted, 3+ is unrooted. State in help/docs that this is structural representation, not proof of biological rooting.
- A unary root or rooted multifurcation is classified as unrooted by this heuristic; help/docs must state that users needing rooted semantics must provide a bifurcating rooted Newick representation.
- Resolve map entries against the taxa present in each tree: empty intersection or non-monophyly/non-split means warning plus skipped endpoint for that tree; a one-taxon match is a tip endpoint.
- For rooted maps, accept only an MRCA whose descendant leaf set equals the present mapped taxon set. For unrooted maps, accept only an internal-edge split whose one side equals the present mapped taxon set.
- Labeled internal nodes are accepted when no `--map` is given; `--map` overrides labels. `label-nodes` is single-tree only and labels all rooted internal nodes but excludes an unrooted artificial root.
- Write requested CSV/TSV files in `<output-dir>/tables/`; all successful non-dry runs write exactly one root `result.json`. No checkpoints or `--resume`: this is a pure-Python utility, not a resumable pipeline.
- `label-nodes` defaults to the separate `runs/posttree/syserror/brlen/label_nodes` directory and writes its three artifacts plus `result.json` directly there, never below the main command output root.
- Enforce `--max-rows` before writing any patristic table; default 5,000,000, and `0` disables the limit. Use `n * (n - 1) / 2` across all accepted trees.
- A non-empty output directory fails unless `--overwrite`; `--dry-run` validates (including endpoint resolution, so an unresolvable single-tree endpoint fails preflight) and returns a payload without creating files. `--overwrite` and `--resume` are not applicable because the command exposes no `--resume`.
- `result.json` must contain the shared top-level fields, `tool_versions: {}`, all resolved params, aggregate key results, warnings, and `data.output_files`; do not emit external-tool fields.
- In single endpoint mode, an unresolved tip/node/map endpoint is exit code 1. In batch endpoint mode, warn and skip only that processing unit. `n_trees_skipped` counts every degenerate parsed tree, every parse-failed or non-empty zero-tree file, and (in batch endpoint mode) every tree whose endpoint cannot be resolved as one skipped processing unit.
- Stream patristic rows in deterministic order and summarize them online; do not collect up to `--max-rows` dictionaries in worker results. Use process workers only for non-patristic modes.
- MCP leaves come only from the Click tree. Remove the legacy `posttree_syserror_brlen` stub; do not hand-register an MCP tool.
- Report methods text exists only for the main analysis command. `label-nodes` remains report-visible as a result but produces no methods prose.
- Do not commit any implementation, documentation, or plan change unless the user explicitly approves a commit.

---

## File Map

| Action | File | Responsibility |
|---|---|---|
| Create | `phyloai/posttree/syserror_brlen.py` | Newick reading, map/labeled-node resolution, all branch-length modes, writers, output lifecycle, label rendering, payloads. |
| Modify | `phyloai/cli/commands/posttree.py` | Register `syserror`, the invoke-without-command `brlen` group callback, and `label-nodes` Click wrapper. |
| Create | `tests/posttree/test_syserror_brlen.py` | Core, runner, map/split, multi-tree, table, result, and labeling regression tests. |
| Create | `tests/cli/test_posttree_syserror_brlen.py` | Help, option validation, dry run, command hierarchy, and CLI error-result tests. |
| Modify | `phyloai/mcp/tools/stubs.py` | Remove only `posttree_syserror_brlen` from stub names and descriptions. |
| Modify | `tests/mcp/test_stubs.py` | Assert brlen is no longer a stub. |
| Modify | `tests/mcp/test_schema_gen.py` | Assert generated main and `label_nodes` MCP schemas. |
| Modify | `phyloai/report/templates.py` | Replace the brlen placeholder with the approved deterministic method generator. |
| Modify | `phyloai/report/collector.py` | Parse `brlen label-nodes` as a distinct fourth-level auxiliary step. |
| Modify | `tests/report/test_collector.py` | Verify main and label-nodes step-ID parsing and ordering. |
| Modify | `tests/report/test_templates.py` | Assert brlen methods text uses tree count, modes, map, and statistics. |
| Create | `docs/commands/posttree-syserror-brlen.md` | Complete English command reference. |
| Create | `docs/commands/posttree-syserror-brlen.zh.md` | Chinese command reference with matching sections and behavior. |
| Modify | `docs/superpowers/specs/2026-08-11-phyloai-posttree-syserror-brlen-design.md` | Change the approved draft status to implemented only after all acceptance checks pass. |
| Modify | `README.md`, `README.zh.md`, `docs/commands/ai-integration.md`, `docs/commands/ai-integration.zh.md` | Add command-table entries, concise examples, and replace the obsolete brlen MCP stub listings. |
| Modify | `skills/phyloai-workflow/SKILL.md` | Add parameter-review, confirmation, and results-interpretation guidance. |
| Modify | `skills/phyloai-workflow/references/parameter-annotations.md` | Add Chinese annotations for every brlen and label-nodes parameter. |

## Task 1: Tree Input, Topology, And Primitive Calculations

**Files:**
- Create: `phyloai/posttree/syserror_brlen.py`
- Create: `tests/posttree/test_syserror_brlen.py`

**Interfaces:**
- Produces `TreeRecord(tree_id: str, tree: Tree)`, `_read_tree_file(path: Path) -> tuple[list[TreeRecord], list[str], int]`, `_is_rooted_representation(tree: Tree) -> bool`, `_branch_length(clade: Clade) -> float`, `_canonical_split(side: frozenset[str], all_tips: frozenset[str]) -> tuple[str, ...]`, and batch row helpers used by Tasks 2-3.

- [ ] **Step 1: Write failing primitive tests**

```python
from Bio import Phylo
from io import StringIO

from phyloai.posttree.syserror_brlen import (
    _branch_length, _canonical_split, _is_rooted_representation,
)


def _tree(text: str):
    return Phylo.read(StringIO(text), "newick")


def test_representation_and_none_length_rules() -> None:
    assert _is_rooted_representation(_tree("((A:1,B:2):3,C:4);"))
    assert not _is_rooted_representation(_tree("(A:1,B:2,C:3);"))
    assert _branch_length(_tree("(A,B);").find_any(name="A")) == 0.0


def test_canonical_split_uses_smaller_then_lexical_side() -> None:
    assert _canonical_split(frozenset({"A", "B"}), frozenset({"A", "B", "C", "D", "E"})) == ("A", "B")
    assert _canonical_split(frozenset({"C", "D"}), frozenset({"A", "B", "C", "D"})) == ("A", "B")
```

- [ ] **Step 2: Run the primitive tests and verify collection fails**

Run: `pytest tests/posttree/test_syserror_brlen.py -q`

Expected: FAIL because `phyloai.posttree.syserror_brlen` does not exist.

- [ ] **Step 3: Implement the minimal tree primitives and multi-tree reader**

```python
@dataclass
class TreeRecord:
    tree_id: str
    tree: Tree


def _is_rooted_representation(tree: Tree) -> bool:
    return len(tree.root.clades) == 2


def _branch_length(clade: Clade) -> float:
    return float(clade.branch_length or 0.0)


def _canonical_split(side: frozenset[str], all_tips: frozenset[str]) -> tuple[str, ...]:
    other = all_tips - side
    candidates = (tuple(sorted(side)), tuple(sorted(other)))
    return min(candidates, key=lambda values: (len(values), values))


def _read_tree_file(path: Path) -> tuple[list[TreeRecord], list[str], int]:
    trees = list(Phylo.parse(path, "newick"))
    records: list[TreeRecord] = []
    warnings: list[str] = []
    skipped = 0
    if not trees:
        return records, [f"tree file {path.name}: no Newick trees found, skipped"], 1
    for index, tree in enumerate(trees):
        if len(tree.get_terminals()) < 2:
            warnings.append(f"tree {path.name}:{index}: fewer than two tips, skipped")
            skipped += 1
            continue
        tree_id = path.name if len(trees) == 1 else f"{path.name}:{index}"
        records.append(TreeRecord(tree_id, tree))
    return records, warnings, skipped
```

Use `Bio.Phylo.parse`, preserve parse order within each file, and convert parser errors to a caller-visible `ValueError` that names the file. Add helpers that produce total `{"tree_file", "total_branch_length"}`, terminal `{"tree_file", "taxon", "branch_length"}`, internal `{"tree_file", "representation", "edge_taxa", "branch_length"}`, and patristic `{"tree_file", "tip1", "tip2", "distance"}` rows. `edge_taxa` is rooted descendants when `representation="rooted"`, and the canonical split side when `representation="unrooted"`. Internal rows must iterate only non-root internal clades: the root has no incoming edge and must never appear as a zero-length internal branch. Sort taxa in identifiers and use `itertools.combinations(sorted(terminals, key=name), 2)` for deterministic patristic output.

- [ ] **Step 4: Add and pass row-level tests**

```python
def test_all_branch_rows_are_deterministic() -> None:
    tree = _tree("((B:2,A:1):3,C:4);")
    assert _total_rows(TreeRecord("x", tree)) == [{"tree_file": "x", "total_branch_length": 10.0}]
    assert _terminal_rows(TreeRecord("x", tree)) == [
        {"tree_file": "x", "taxon": "A", "branch_length": 1.0},
        {"tree_file": "x", "taxon": "B", "branch_length": 2.0},
        {"tree_file": "x", "taxon": "C", "branch_length": 4.0},
    ]


def test_internal_rows_exclude_root_and_all_missing_lengths_warn() -> None:
    record = TreeRecord("x", _tree("((A,B),C);"))
    assert len(_internal_rows(record)) == 1
    assert _missing_length_warning(record) == "tree x: all branch lengths are missing; treating them as 0.0"


def test_internal_rows_use_one_schema_for_mixed_representations() -> None:
    rooted = _internal_rows(TreeRecord("rooted", _tree("((A:1,B:1):2,C:1);")))[0]
    unrooted = _internal_rows(TreeRecord("unrooted", _tree("((A:1,B:1):2,C:1,D:1);")))[0]
    assert list(rooted) == list(unrooted) == ["tree_file", "representation", "edge_taxa", "branch_length"]
    assert rooted["representation"] == "rooted"
    assert unrooted["representation"] == "unrooted"
```

Run: `pytest tests/posttree/test_syserror_brlen.py -q`

Expected: PASS.

- [ ] **Step 5: Review Task 1 diff; do not commit**

Run: `git diff --check && git diff -- phyloai/posttree/syserror_brlen.py tests/posttree/test_syserror_brlen.py`

Expected: no whitespace errors and only the core/test files changed.

## Task 2: Map And Labeled Node Resolution

**Files:**
- Modify: `phyloai/posttree/syserror_brlen.py`
- Modify: `tests/posttree/test_syserror_brlen.py`

**Interfaces:**
- Consumes `TreeRecord` and topology helpers from Task 1.
- Produces `_parse_map(path: Path) -> dict[str, frozenset[str]]`, `_resolve_endpoint(tree: Tree, name: str, node_map: dict[str, frozenset[str]] | None) -> ResolvedEndpoint`, `_endpoint_distance(tree: Tree, left: ResolvedEndpoint, right: ResolvedEndpoint) -> float`, and endpoint rows for Task 3.

- [ ] **Step 1: Write failing resolution tests**

```python
import pytest


def test_rooted_map_allows_missing_taxa_but_requires_present_clade() -> None:
    tree = _tree("((A:1,B:1):2,C:1);")
    endpoint = _resolve_endpoint(tree, "AB", {"AB": frozenset({"A", "B", "Missing"})})
    assert endpoint.node_type == "internal"
    assert endpoint.present_taxa == frozenset({"A", "B"})
    with pytest.raises(ValueError, match="not monophyletic"):
        _resolve_endpoint(tree, "AC", {"AC": frozenset({"A", "C"})})


def test_unrooted_map_matches_either_side_of_an_internal_edge() -> None:
    tree = _tree("((A:1,B:1):2,C:1,D:1);")
    endpoint = _resolve_endpoint(tree, "CD", {"CD": frozenset({"C", "D"})})
    assert endpoint.node_type == "internal"


def test_map_overrides_internal_node_label_and_single_taxon_is_tip(tmp_path) -> None:
    tree = _tree("((A:1,B:1)N1:2,C:1);")
    endpoint = _resolve_endpoint(tree, "N1", {"N1": frozenset({"A"})})
    assert endpoint.node_type == "tip"
```

- [ ] **Step 2: Run the tests and verify failure**

Run: `pytest tests/posttree/test_syserror_brlen.py -q`

Expected: FAIL because `_resolve_endpoint` and `ResolvedEndpoint` are undefined.

- [ ] **Step 3: Implement map parsing and exact rooted/unrooted semantics**

```python
@dataclass
class ResolvedEndpoint:
    clade: Clade
    node_type: str
    present_taxa: frozenset[str]


def _parse_map(path: Path) -> dict[str, frozenset[str]]:
    mapping: dict[str, frozenset[str]] = {}
    for raw_line in path.read_text().splitlines():
        if ":" not in raw_line:
            continue
        label, raw_taxa = raw_line.split(":", 1)
        taxa = frozenset(value.strip() for value in raw_taxa.split(",") if value.strip())
        if label.strip() and taxa:
            mapping[label.strip()] = taxa
    return mapping
```

When no map is supplied, match `name` first against internal `clade.name`, then against terminal `clade.name`; otherwise raise `ValueError("node '<name>' not found in labeled tree")`. When map is supplied, never fall back to labels. For rooted trees use `tree.common_ancestor(sorted(present_taxa))`, then compare its terminal-name set exactly. For unrooted trees, build `{child: parent for parent in tree.find_clades() for child in parent.clades}` once, then examine every non-root internal clade: its descendants and the complement define a split. A match must equal either side; use the clade for a descendant-side match and `parents[clade]` for a complement-side match. Store a warning that names the complement-side match. Do not reroot or mutate the input tree.

Use `tree.distance(left.clade, right.clade)` for node-to-node and `tree.distance(node.clade, tip.clade)` for node-to-tip. For map-driven `node-to-tip` without `--tip1`, emit one row for each `present_taxa` taxon, sorted; do not expand to MRCA descendants. Without a map, allow omitted `--tip1` only for a rooted labeled internal node and emit its sorted descendants; reject the same request for unrooted input.

- [ ] **Step 4: Add endpoint validation and row tests, then run them**

```python
def test_node_to_tip_without_map_is_rooted_labeled_descendants_only() -> None:
    tree = _tree("((A:1,B:2)N1:3,C:4);")
    rows = _node_to_tip_rows(TreeRecord("x", tree), "N1", None, None)
    assert [row["tip"] for row in rows] == ["A", "B"]


def test_unrooted_labeled_node_to_tip_requires_tip_or_map() -> None:
    with pytest.raises(ValueError, match="--tip1 or --map required"):
        _node_to_tip_rows(TreeRecord("x", _tree("(A:1,B:1,C:1)N1;")), "N1", None, None)
```

Run: `pytest tests/posttree/test_syserror_brlen.py -q`

Expected: PASS.

- [ ] **Step 5: Review Task 2 diff; do not commit**

Run: `git diff --check && git diff -- phyloai/posttree/syserror_brlen.py tests/posttree/test_syserror_brlen.py`

Expected: no whitespace errors.

## Task 3: Main Runner, Files, And Result Payload

**Files:**
- Modify: `phyloai/posttree/syserror_brlen.py`
- Modify: `tests/posttree/test_syserror_brlen.py`

**Interfaces:**
- Produces `run_brlen(tree: Path | None, tree_dir: Path | None, mode: str, map_file: Path | None = None, node1: str | None = None, node2: str | None = None, tip1: str | None = None, tip2: str | None = None, table_format: str = "csv", threads: int = 4, max_rows: int = 5_000_000, output_dir: Path = Path("runs/posttree/syserror/brlen"), overwrite: bool = False, dry_run: bool = False, quiet: bool = False) -> dict[str, Any]`.

- [ ] **Step 1: Write failing runner tests**

```python
import csv
import json
import pytest


def test_batch_all_writes_four_csv_tables_and_json(tmp_path) -> None:
    trees = tmp_path / "trees"
    trees.mkdir()
    (trees / "one.tre").write_text("((A:1,B:2):3,C:4);")
    result = run_brlen(tree=None, tree_dir=trees, mode="all", output_dir=tmp_path / "out", quiet=True)
    assert result["status"] == "success"
    assert set(result["key_results"]["modes"]) == {"total", "terminal", "internal", "patristic"}
    assert (tmp_path / "out" / "tables" / "total.csv").exists()
    assert json.loads((tmp_path / "out" / "result.json").read_text())["error"] is None


def test_patristic_limit_is_checked_before_writing(tmp_path) -> None:
    tree_file = tmp_path / "tree.nwk"
    tree_file.write_text("(A:1,B:1,C:1,D:1);")
    with pytest.raises(ValueError, match="--max-rows"):
        run_brlen(tree=tree_file, tree_dir=None, mode="patristic", max_rows=5, output_dir=tmp_path / "out", quiet=True)
    assert not (tmp_path / "out" / "tables").exists()
```

- [ ] **Step 2: Run the runner tests and verify failure**

Run: `pytest tests/posttree/test_syserror_brlen.py -q`

Expected: FAIL because `run_brlen` is undefined.

- [ ] **Step 3: Implement validation, scanning, writing, and aggregation**

Implement this control flow exactly:

```python
# 1. Validate one input source, table_format, threads >= 1, max_rows >= 0,
#    map readability, and mode grammar before claiming output_dir.
# 2. Expand mode="all" to ["total", "terminal", "internal", "patristic"].
# 3. Read --tree exactly once or scan sorted non-empty regular --tree-dir files.
#    A directory parse failure or non-empty zero-tree file increments skipped once
#    and becomes warnings.append(...); a single failure raises.
# 4. Estimate all requested patristic rows, append the exact estimate to warnings,
#    and reject over limit before mkdir/write.
# 5. If non-dry-run, enforce output lifecycle then mkdir tables/.
# 6. Process every record; endpoint resolution ValueError is a single-mode error
#    or a batch warning that skips only that record. Merge rows in input-tree order.
# 7. Stream patristic rows directly to its writer and maintain online count/sum/
#    sum-of-squares/min/max; write other requested table rows normally.
# 8. Summarize each numeric mode as n_values, mean, population SD, min, max;
#    build result.json and return it.
```

Use `ProcessPoolExecutor(max_workers=threads)` only for non-patristic directory calculations and only after parsing/row-limit validation. Each worker receives one file path and returns serializable non-patristic rows/warnings; preserve alphabetical file order when merging futures. A multi-tree file is one worker unit. For `patristic` or `all`, process files serially in that same deterministic order and write distances as generated, so no worker result contains millions of row dictionaries. Single mode executes directly. Processing (including endpoint resolution) runs in `--dry-run` too, so an unresolvable single-tree endpoint fails preflight and batch endpoint skips are reported; `dry_run` returns a success payload with validated `params`, resolved mode list, empty `output_files`, and no filesystem effects. A transient Rich progress bar ("Processing trees") is shown when not `--quiet`/`--dry-run`: it advances per file for batch non-patristic modes and per tree for patristic streaming.

Payload invariants:

```python
key_results = {
    "n_trees": n_processed,
    "n_trees_skipped": n_skipped,
    "modes": modes,
    "summary": per_mode_summary,
}

data = {
    "summary": {"n_trees_processed": n, "n_trees_skipped": skipped, "n_multi_tree_files": multi},
    "warnings": warnings,
    "output_files": {
        "terminal_table": {"path": str(path), "description": "Terminal branch lengths per taxon per tree"}
    },
}
```

For endpoint outputs, use `tables/tip_to_tip.<suffix>`, `node_to_node.<suffix>`, or `node_to_tip.<suffix>` and the columns in design section 3.2. Include endpoint node-type columns exactly where specified.

- [ ] **Step 4: Add failure-mode and format tests, then run them**

```python
def test_batch_counts_parse_failed_and_degenerate_processing_units(tmp_path) -> None:
    trees = tmp_path / "trees"
    trees.mkdir()
    (trees / "bad.txt").write_text("((A,B),C;")
    (trees / "degenerate.txt").write_text("not a tree")
    (trees / "posterior.trees").write_text("(A:1,B:1);\n(A:2,B:2);")
    result = run_brlen(tree=None, tree_dir=trees, mode="total", output_dir=tmp_path / "out", quiet=True)
    rows = list(csv.DictReader((tmp_path / "out" / "tables" / "total.csv").open()))
    assert [row["tree_file"] for row in rows] == ["posterior.trees:0", "posterior.trees:1"]
    assert result["key_results"]["n_trees_skipped"] == 2
    assert any("failed to parse" in warning for warning in result["data"]["warnings"])
    assert any("fewer than two tips" in warning for warning in result["data"]["warnings"])


def test_tsv_and_endpoint_missing_tip_warns_not_fails_batch(tmp_path) -> None:
    trees = tmp_path / "trees"
    trees.mkdir()
    (trees / "a.tre").write_text("(A:1,B:1);")
    result = run_brlen(tree=None, tree_dir=trees, mode="tip-to-tip", tip1="A", tip2="Missing", table_format="tsv", output_dir=tmp_path / "out", quiet=True)
    assert "Missing" in result["data"]["warnings"][0]
    assert (tmp_path / "out" / "tables" / "tip_to_tip.tsv").read_text().startswith("tree_file\ttip1\ttip2\tdistance")


def test_single_endpoint_missing_tip_is_an_error(tmp_path) -> None:
    tree_file = tmp_path / "one.tre"
    tree_file.write_text("(A:1,B:1);")
    with pytest.raises(ValueError, match="Missing"):
        run_brlen(tree=tree_file, tree_dir=None, mode="tip-to-tip", tip1="A", tip2="Missing", output_dir=tmp_path / "out", quiet=True)


def test_missing_lengths_and_patristic_estimate_are_warnings(tmp_path) -> None:
    tree_file = tmp_path / "tree.nwk"
    tree_file.write_text("(A,B,C);")
    result = run_brlen(tree=tree_file, tree_dir=None, mode="patristic", output_dir=tmp_path / "out", quiet=True)
    warnings = result["data"]["warnings"]
    assert any("all branch lengths are missing" in warning for warning in warnings)
    assert any("estimated patristic rows: 3" in warning for warning in warnings)
    assert result["key_results"]["n_trees"] == 1
    assert result["key_results"]["n_trees_skipped"] == 0
```

Run: `pytest tests/posttree/test_syserror_brlen.py -q`

Expected: PASS.

- [ ] **Step 5: Review Task 3 diff; do not commit**

Run: `git diff --check && git diff --stat`

Expected: no whitespace errors; only planned files changed.

## Task 4: `label-nodes` Helper

**Files:**
- Modify: `phyloai/posttree/syserror_brlen.py`
- Modify: `tests/posttree/test_syserror_brlen.py`

**Interfaces:**
- Produces `run_label_nodes(tree: Path, output_dir: Path = Path("runs/posttree/syserror/brlen/label_nodes"), overwrite: bool = False, quiet: bool = False) -> dict[str, Any]`.

- [ ] **Step 1: Write failing label-node tests**

```python
def test_label_nodes_writes_labeled_tree_map_pdf_and_result(tmp_path) -> None:
    tree_file = tmp_path / "species.nwk"
    tree_file.write_text("((A:1,B:1):2,C:1);")
    result = run_label_nodes(tree_file, output_dir=tmp_path / "out", quiet=True)
    labeled = tmp_path / "out" / "species.labeled.nwk"
    assert "N1" in labeled.read_text()
    assert (tmp_path / "out" / "species.map.txt").exists()
    assert (tmp_path / "out" / "species.labeled.pdf").exists()
    assert result["key_results"]["n_internal_nodes_labeled"] == 2


def test_unrooted_artificial_root_is_not_labeled(tmp_path) -> None:
    tree_file = tmp_path / "unrooted.nwk"
    tree_file.write_text("((A:1,B:1):1,C:1,D:1);")
    run_label_nodes(tree_file, output_dir=tmp_path / "out", quiet=True)
    map_rows = (tmp_path / "out" / "unrooted.map.txt").read_text().splitlines()
    assert len(map_rows) == 1


def test_label_output_isolated_from_main_result(tmp_path) -> None:
    tree_file = tmp_path / "species.nwk"
    tree_file.write_text("((A:1,B:1):2,C:1);")
    main_dir = tmp_path / "brlen"
    run_brlen(tree=tree_file, tree_dir=None, mode="total", output_dir=main_dir, quiet=True)
    run_label_nodes(tree_file, output_dir=main_dir / "label_nodes", quiet=True)
    assert (main_dir / "result.json").exists()
    assert (main_dir / "label_nodes" / "result.json").exists()
```

- [ ] **Step 2: Run label-node tests and verify failure**

Run: `pytest tests/posttree/test_syserror_brlen.py -q`

Expected: FAIL because `run_label_nodes` is undefined.

- [ ] **Step 3: Implement labeling and rendering**

```python
def _label_internal_nodes(tree: Tree) -> list[Clade]:
    rooted = _is_rooted_representation(tree)
    nodes = [node for node in tree.find_clades(order="preorder") if node.clades]
    if not rooted:
        nodes = [node for node in nodes if node is not tree.root]
    for index, node in enumerate(nodes, 1):
        node.name = f"N{index}"
        node.confidence = None
    return nodes
```

Reject non-Newick/unreadable input and a non-empty output directory without `--overwrite`. The default output directory is `runs/posttree/syserror/brlen/label_nodes`; output `<stem>.labeled.nwk`, `<stem>.map.txt`, `<stem>.labeled.pdf`, and `result.json` directly in that directory. Generate map rows in label/preorder order as `Nxx:taxon1,taxon2`; sort taxa. Labels are unpadded (`N1`, `N2`, ... `Nxx`) and clear each labeled node's support `confidence` so `labeled.nwk` never mixes labels with support values; branch lengths are written losslessly with `format_branch_length="%r"`. Use matplotlib's noninteractive `Agg` backend before importing pyplot, call `Phylo.draw(tree, axes=ax, do_show=False)`, write a PDF, and close the figure. The payload has `n_internal_nodes_labeled`, `n_terminals`, and `rooted`; `data.output_files` describes all three files; the recorded `command` always includes `-o <output-dir>`.

- [ ] **Step 4: Run labeling tests and inspect outputs**

Run: `pytest tests/posttree/test_syserror_brlen.py -q`

Expected: PASS.

- [ ] **Step 5: Review Task 4 diff; do not commit**

Run: `git diff --check && git diff -- phyloai/posttree/syserror_brlen.py tests/posttree/test_syserror_brlen.py`

Expected: no whitespace errors.

## Task 5: CLI, MCP, And Report Integration

**Files:**
- Modify: `phyloai/cli/commands/posttree.py`
- Create: `tests/cli/test_posttree_syserror_brlen.py`
- Modify: `phyloai/mcp/tools/stubs.py`
- Modify: `tests/mcp/test_stubs.py`
- Modify: `tests/mcp/test_schema_gen.py`
- Modify: `phyloai/report/templates.py`
- Modify: `phyloai/report/collector.py`
- Modify: `tests/report/test_collector.py`
- Modify: `tests/report/test_templates.py`

**Interfaces:**
- Consumes `run_brlen()` and `run_label_nodes()` from Tasks 3-4; existing `_fail()` and `_write_error_result_json()`.
- Produces Click/MCP leaves `posttree_syserror_brlen` and `posttree_syserror_brlen_label_nodes`, and full brlen report methods text.

- [ ] **Step 1: Write failing CLI/MCP/report tests**

```python
from click.testing import CliRunner


def test_brlen_help_lists_modes_map_and_patristic_guard() -> None:
    result = CliRunner().invoke(cli, ["posttree", "syserror", "brlen", "--help"])
    assert result.exit_code == 0
    for text in ("--tree-dir", "tip-to-tip", "node-to-tip", "--map", "--max-rows", "O(n²)"):
        assert text in result.output


def test_brlen_click_leaves_replace_stub() -> None:
    names = {item["tool_name"] for item in walk_click_tree(cli)}
    assert {"posttree_syserror_brlen", "posttree_syserror_brlen_label_nodes"} <= names
    assert "posttree_syserror_brlen" not in STUB_TOOL_NAMES


def test_brlen_methods_text_is_quantitative() -> None:
    text = generate_all_methods("posttree.syserror.brlen", {"map": "nodes.map"}, {"n_trees": 20, "modes": ["terminal"], "summary": {"terminal": {"mean": 0.1, "sd": 0.02}}}, {})
    assert "20" in text and "map" in text.lower() and "0.1000" in text


def test_label_nodes_has_distinct_report_step_without_methods_text() -> None:
    command = "phyloai posttree syserror brlen label-nodes --tree species.nwk"
    assert parse_step_id(command) == "posttree.syserror.brlen.label-nodes"
    assert "posttree.syserror.brlen.label-nodes" in STEP_ORDER
    assert generate_all_methods("posttree.syserror.brlen.label-nodes", {}, {}, {}) == ""
```

- [ ] **Step 2: Run integration tests and verify failure**

Run: `pytest tests/cli/test_posttree_syserror_brlen.py tests/mcp/test_stubs.py tests/mcp/test_schema_gen.py tests/report/test_collector.py tests/report/test_templates.py -q`

Expected: FAIL because `syserror` is absent and brlen remains an MCP stub.

- [ ] **Step 3: Add the Click hierarchy and thin wrappers**

```python
class _SyserrorGroup(click.Group):
    def list_commands(self, ctx: click.Context) -> list[str]:
        return ["brlen"]


@posttree.group("syserror", cls=_SyserrorGroup)
def syserror() -> None:
    """Atomic systematic-error diagnostic operations."""


@syserror.group("brlen", invoke_without_command=True)
@click.pass_context
def brlen(ctx: click.Context, ...all main options...) -> None:
    if ctx.invoked_subcommand is None:
        run_brlen(...)
```

Add every option and default in design section 3.1, including `click.Choice` for `table-format`, `click.IntRange(1)` for threads, and `click.IntRange(0)` for max rows. The brlen docstring must contain all mode rules, map example, rooted/unrooted caveat, multi-tree support, O(n²) guard, and the seven design examples. Register `@brlen.command("label-nodes")` with only `--tree`, `-o/--output-dir`, `--overwrite`, and `-q/--quiet`.

Update `_PosttreeGroup.list_commands()` to include `syserror`. Convert library `ValueError` into `_write_error_result_json(..., "input")` followed by `_fail(..., 1)` only after the output directory is valid to write; preserve the runner's preflight no-write behavior for a rejected preexisting directory. Construct a complete reproducible command string using `shlex.join`.

- [ ] **Step 4: Remove the obsolete stub and implement report step handling**

Remove only the brlen item and description from `STUB_TOOL_NAMES` / `_DESCRIPTIONS`; retain `cca` and `sites`. Replace `generate_methods_posttree_syserror_brlen()` with the design-section-7 generator, using `_describe_n` and `_safe_fmt` already in the module. It must read `key_results["n_trees"]`, describe modes, map use, and skipped trees, and emit one mean/SD sentence per analyzed mode in `key_results["summary"]` (total, terminal, internal, patristic, and endpoint modes) via the `_MODE_MEAN_PHRASES` map.

Add `"brlen": {"label-nodes"}` to `phyloai/report/collector.py`'s `_FOURTH_LEVEL`, add `"posttree.syserror.brlen.label-nodes"` directly after the main brlen step in `STEP_ORDER`, and register:

```python
def generate_methods_posttree_syserror_brlen_label_nodes(
    params: dict[str, Any], key_results: dict[str, Any], tool_versions: dict[str, Any],
) -> str:
    return ""
```

under the corresponding fourth-level key in `METHODS_GENERATORS`. This ensures the helper is collected and displayed but never contributes analytical methods prose.

- [ ] **Step 5: Run focused integration tests**

Run: `pytest tests/cli/test_posttree_syserror_brlen.py tests/mcp/test_stubs.py tests/mcp/test_schema_gen.py tests/report/test_collector.py tests/report/test_templates.py -q`

Expected: PASS.

- [ ] **Step 6: Review Task 5 diff; do not commit**

Run: `git diff --check && git diff --stat`

Expected: no whitespace errors and no hand-authored MCP tool schema.

## Task 6: Documentation, Workflow Guidance, And Final Verification

**Files:**
- Create: `docs/commands/posttree-syserror-brlen.md`
- Create: `docs/commands/posttree-syserror-brlen.zh.md`
- Modify: `docs/superpowers/specs/2026-08-11-phyloai-posttree-syserror-brlen-design.md`
- Modify: `README.md`
- Modify: `README.zh.md`
- Modify: `docs/commands/ai-integration.md`
- Modify: `docs/commands/ai-integration.zh.md`
- Modify: `skills/phyloai-workflow/SKILL.md`
- Modify: `skills/phyloai-workflow/references/parameter-annotations.md`

**Interfaces:**
- Consumes the final Click help and output schemas.
- Produces accurate bilingual user documentation and workflow guidance. The parent design already lists `syserror brlen`; do not modify it merely to restate this implemented command.

- [ ] **Step 1: Write the English command reference**

Create `docs/commands/posttree-syserror-brlen.md` with the established sections: Purpose, Usage, Inputs, Modes, Node Identification, Outputs, Examples, Warnings/Errors, and Notes. Include every main and label-nodes option with its default. Include both map parsing examples and this explicit behavior:

```markdown
`--map` uses the taxa present in each tree. A present subset can resolve when it
is an exact rooted clade or unrooted split; an empty overlap or incompatible
group emits a warning and skips that tree's endpoint calculation. `Nxx` labels
are suitable only for the reference topology that produced them.
```

List the seven tables and clarify that only requested modes create tables. State that `terminal` reports all terminal branches and users should filter its table for one taxon. Document `internal.csv` as `tree_file,representation,edge_taxa,branch_length`, and distinguish rooted descendant taxa from unrooted canonical split taxa. Explain the bifurcating-root requirement, single-mode endpoint exit-1 versus batch warning/skip behavior, and that output branch lengths are substitutions per site and do not independently distinguish elapsed time from rate.

- [ ] **Step 2: Write the Chinese reference and update README entries**

Create `posttree-syserror-brlen.zh.md` with the same section order, option coverage, warnings, and examples in Chinese. Add this example to both READMEs beside other post-tree commands:

```bash
phyloai posttree syserror brlen --tree-dir posterior_trees --mode node-to-tip \
  --map nodes.map.txt --node1 Collembola -o runs/posttree/syserror/brlen
```

Add a concise table entry linking each README to the new appropriate-language command document. In both AI integration documents, replace `posttree_syserror_brlen` in the stub list with `posttree_syserror_brlen` and `posttree_syserror_brlen_label_nodes` entries in the generated CLI-tools table; retain `cca` and `sites` as stubs. Do not claim this command identifies LBA causality; it extracts diagnostic branch-length measurements.

- [ ] **Step 3: Add workflow and parameter-review guidance**

In `skills/phyloai-workflow/SKILL.md`, require the standard parameter review plus explicit user approval before running `posttree syserror brlen`; no `doctor` is needed because it has no external executable. Guide users to choose `--map` for posterior/model tree batches and `label-nodes` only for an inspected stable reference topology. In results interpretation, distinguish terminal, internal, patristic, node-to-node, and node-to-tip measurements; advise comparing distributions across model runs without declaring model superiority from branch lengths alone.

In `parameter-annotations.md`, add Chinese annotations for `--tree`, `--tree-dir`, `--mode`, `--map`, `--node1`, `--node2`, `--tip1`, `--tip2`, `--output-dir`, `--table-format`, `--threads`, `--max-rows`, `--overwrite`, `--dry-run`, `--quiet`, and label-nodes options. State mutex input rules, mode dependencies, map/split behavior, row-limit semantics, and no-resume behavior.

- [ ] **Step 4: Update implementation status only after the preceding tests pass**

Change the brlen-specific design front matter from `Status: Draft — pending user review` to `Status: Implemented`, and append the implementation date. Do not modify the parent `2026-06-07` design: it already reserves the `syserror brlen` command and requires no semantic correction.

- [ ] **Step 5: Verify help, documentation references, and full regression suite**

Run: `python -m phyloai.cli.main posttree syserror brlen --help`

Expected: exit 0; detailed modes, map format, rooted/unrooted caveat, and O(n²) warning are visible.

Run: `python -m phyloai.cli.main posttree syserror brlen label-nodes --help`

Expected: exit 0; only label-node options are visible.

Run: `pytest tests/posttree/test_syserror_brlen.py tests/cli/test_posttree_syserror_brlen.py tests/mcp/test_stubs.py tests/mcp/test_schema_gen.py tests/report/test_collector.py tests/report/test_templates.py -q`

Expected: PASS.

Run: `pytest -q`

Expected: PASS.

- [ ] **Step 6: Final review; do not commit**

Run: `git diff --check && git status --short && git diff -- docs/commands README.md README.zh.md skills/phyloai-workflow`

Expected: no whitespace errors; only intended files are changed. Do not commit unless the user separately approves a commit.

## Final Acceptance Checklist

- [ ] One tree, a directory, and a multi-tree Newick file produce deterministic CSV/TSV output and spec-compliant `result.json`.
- [ ] All batch and endpoint modes validate their exclusive options and required parameters.
- [ ] Rooted maps validate exact present clades; unrooted maps validate exact splits; missing/incompatible endpoints warn and skip the affected batch tree.
- [ ] `label-nodes` creates `Nxx` tree, map, and PDF with the correct rooted/unrooted root behavior.
- [ ] The generated MCP main and helper leaves replace the brlen stub; no manual MCP schema is added.
- [ ] Reports produce quantitative brlen methods text; bilingual docs and workflow guidance match CLI help.
- [ ] Targeted tests and `pytest -q` pass; changes remain uncommitted until explicit user approval.
