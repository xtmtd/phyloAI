# phyloai tree cf Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Implement `phyloai tree cf` — a concordance factor computation command supporting gCF, sCF, sCFl, gCF+sCF (via IQ-TREE3), and qCF (via wASTRAL) with a `--cf` mode selector.

**Architecture:** Single-file library (`phyloai/tree/cf.py`) with `run_cf()` entry point, dispatched by `--cf` mode. No `--tool-args`, no `--resume`. Follows `tree msc` conventions for CLI registration, input validation, `result.json` schemas, and output directory management. Gene tree merging reuses the pattern from `phyloai/tree/msc.py`.

**Tech Stack:** Python 3.10+, Click, Bio.Phylo, Rich, subprocess (IQ-TREE3, wASTRAL)

**Spec:** `docs/superpowers/specs/2026-06-21-phyloai-tree-cf-design.md`

---

## File Structure

| File | Action | Responsibility |
|------|--------|----------------|
| `phyloai/tree/cf.py` | Create | Library: `run_cf()`, command builders, qCF mapper, gene tree merger |
| `tests/tree/test_cf.py` | Create | Library-level unit tests |
| `tests/cli/test_tree_cf.py` | Create | CLI integration tests |
| `phyloai/cli/commands/tree.py` | Modify | Register `cf` command, update `_TreeGroup.list_commands` |
| `docs/commands/tree-cf.md` | Create | User-facing command documentation |

---

### Task 1: Create library test file with qCF mapper tests (TDD)

**Files:**
- Create: `tests/tree/test_cf.py`

- [ ] **Step 1: Write test file skeleton**

```python
"""Tests for phyloai/tree/cf.py — concordance factor computation."""

from __future__ import annotations

from pathlib import Path

import pytest

# All functions tested will be imported from phyloai.tree.cf once written.
# For now, we write tests that will fail until implementation exists.
```

- [ ] **Step 2: Write test for _map_qcf_to_tree with existing support**

```python
def test_map_qcf_to_tree_appends_qcf_to_existing_support(tmp_path: Path) -> None:
    """q1=0.422083 -> 42, appended to '100/90' as '100/90/42'."""
    from phyloai.tree.cf import _map_qcf_to_tree
    from Bio import Phylo
    from io import StringIO

    ref_nwk = tmp_path / "ref.nwk"
    ref_nwk.write_text("((A:0.1,B:0.1)100/90:0.05,C:0.15);")

    wastral_nwk = tmp_path / "wastral.tre"
    wastral_nwk.write_text(
        "((A:0.1,B:0.1)'[q1=0.422083;q2=0.288958;q3=0.288958]':0.05,C:0.15);"
    )

    output_nwk = tmp_path / "qCF.cf.tree"
    _map_qcf_to_tree(ref_nwk, wastral_nwk, output_nwk)

    # Parse result and check the internal (A,B) node label exactly
    tree = Phylo.read(str(output_nwk), "newick")
    # Select the clade with exactly {A, B} as terminals (not root with {A,B,C})
    ab_node = None
    for c in tree.find_clades():
        if not c.is_terminal():
            leaves = frozenset(t.name for t in c.get_terminals())
            if leaves == frozenset({"A", "B"}):
                ab_node = c
                break
    assert ab_node is not None, "Could not find (A,B) clade"
    # Support is stored in clade.name by our mapper
    assert ab_node.name == "100/90/42"
```

- [ ] **Step 3: Write test for _map_qcf_to_tree with no existing support**

```python
def test_map_qcf_to_tree_no_existing_support(tmp_path: Path) -> None:
    """When ref tree has no support, qCF becomes the sole label."""
    from phyloai.tree.cf import _map_qcf_to_tree
    from Bio import Phylo
    from io import StringIO

    ref_nwk = tmp_path / "ref.nwk"
    ref_nwk.write_text("((A:0.1,B:0.1):0.05,C:0.15);")

    wastral_nwk = tmp_path / "wastral.tre"
    wastral_nwk.write_text(
        "((A:0.1,B:0.1)'[q1=0.750000;q2=0.125000;q3=0.125000]':0.05,C:0.15);"
    )

    output_nwk = tmp_path / "qCF.cf.tree"
    _map_qcf_to_tree(ref_nwk, wastral_nwk, output_nwk)

    tree = Phylo.read(str(output_nwk), "newick")
    ab_node = None
    for c in tree.find_clades():
        if not c.is_terminal():
            leaves = frozenset(t.name for t in c.get_terminals())
            if leaves == frozenset({"A", "B"}):
                ab_node = c
                break
    assert ab_node is not None
    assert ab_node.name == "75"
```

- [ ] **Step 4: Write test for complementary rooting — 4-taxon tree, qCF on opposite side**

```python
def test_map_qcf_to_tree_handles_complementary_rooting(tmp_path: Path) -> None:
    """4-taxon: ref annotates (A,B) side but wASTRAL labels (C,D) side.
    Canonical bipartition matching must handle the complement."""
    from phyloai.tree.cf import _map_qcf_to_tree
    from Bio import Phylo
    from io import StringIO

    # Reference: ((A,B),(C,D)) with support 100 on (A,B) clade
    ref_nwk = tmp_path / "ref.nwk"
    ref_nwk.write_text("((A:0.1,B:0.1)100:0.05,(C:0.1,D:0.1):0.05);")

    # wASTRAL: rooted differently so the internal node labels (C,D) not (A,B)
    # Same bipartition (A,B)|(C,D), just the complement side is labeled
    wastral_nwk = tmp_path / "wastral.tre"
    wastral_nwk.write_text(
        "((A:0.1,B:0.1):0.05,(C:0.1,D:0.1)'[q1=0.6;q2=0.2;q3=0.2]':0.05);"
    )

    output_nwk = tmp_path / "qCF.cf.tree"
    _map_qcf_to_tree(ref_nwk, wastral_nwk, output_nwk)

    # Both internal nodes exist; the (A,B) node should get '100/60'
    tree = Phylo.read(str(output_nwk), "newick")
    ab_clade = None
    for c in tree.find_clades():
        if not c.is_terminal():
            leaves = frozenset(t.name for t in c.get_terminals())
            if leaves == frozenset({"A", "B"}):
                ab_clade = c
                break
    assert ab_clade is not None
    assert ab_clade.name == "100/60"
```

- [ ] **Step 5: Write test for _map_qcf_to_tree with multiple internal nodes**

```python
def test_map_qcf_to_tree_multiple_internal_nodes(tmp_path: Path) -> None:
    """4-taxon tree with two internal nodes, both should get qCF annotations."""
    from phyloai.tree.cf import _map_qcf_to_tree

    ref_nwk = tmp_path / "ref.nwk"
    ref_nwk.write_text("(((A:0.1,B:0.1)100:0.05,C:0.1)80:0.02,D:0.15);")

    wastral_nwk = tmp_path / "wastral.tre"
    wastral_nwk.write_text(
        "(((A:0.1,B:0.1)'[q1=0.9;]':0.05,C:0.1)'[q1=0.7;]':0.02,D:0.15);"
    )

    output_nwk = tmp_path / "qCF.cf.tree"
    _map_qcf_to_tree(ref_nwk, wastral_nwk, output_nwk)

    result = output_nwk.read_text().strip()
    # Both nodes should have qCF values appended
    assert result.count("/") >= 3  # at least 100/90 for one node and 80/70 for the other
```

- [ ] **Step 6: Write test for _map_qcf_to_tree when wastral node has no q1 (graceful skip)**

```python
def test_map_qcf_to_tree_missing_q1_skips_node(tmp_path: Path) -> None:
    """Node without q1 annotation should be left unchanged."""
    from phyloai.tree.cf import _map_qcf_to_tree

    ref_nwk = tmp_path / "ref.nwk"
    ref_nwk.write_text("((A:0.1,B:0.1)100:0.05,C:0.15);")

    # No q1 annotation in the label
    wastral_nwk = tmp_path / "wastral.tre"
    wastral_nwk.write_text(
        "((A:0.1,B:0.1)'[pp1=0.9]':0.05,C:0.15);"
    )

    output_nwk = tmp_path / "qCF.cf.tree"
    _map_qcf_to_tree(ref_nwk, wastral_nwk, output_nwk)

    result = output_nwk.read_text().strip()
    # Node should still have 100 (no qCF appended since no q1)
    assert "100" in result
```

- [ ] **Step 7: Run tests to confirm they fail**

```bash
pytest tests/tree/test_cf.py -v
```

Expected: all tests FAIL with `ModuleNotFoundError: No module named 'phyloai.tree.cf'`.

---

### Task 2: Write command builder tests (TDD)

**Files:**
- Modify: `tests/tree/test_cf.py` (append)

- [ ] **Step 1: Write test for _build_iqtree_cf_cmd — gcf mode**

```python
def test_build_iqtree_cf_cmd_gcf(tmp_path: Path) -> None:
    from phyloai.tree.cf import _build_iqtree_cf_cmd

    ref_tree = tmp_path / "species.nwk"
    ref_tree.write_text("(A,B);")
    gene_trees = tmp_path / "merged.trees"
    gene_trees.write_text("(A,B);")

    cmd = _build_iqtree_cf_cmd(
        cf_mode="gcf",
        executable="iqtree3",
        ref_tree=ref_tree,
        gene_trees=gene_trees,
        matrix=None,
        scf_quartets=100,
        prefix="gCF",
        threads=4,
    )
    assert cmd == [
        "iqtree3", "-t", str(ref_tree), "--gcf", str(gene_trees),
        "--prefix", "gCF", "-T", "4",
    ]
```

- [ ] **Step 2: Write test for _build_iqtree_cf_cmd — scf mode**

```python
def test_build_iqtree_cf_cmd_scf(tmp_path: Path) -> None:
    from phyloai.tree.cf import _build_iqtree_cf_cmd

    ref_tree = tmp_path / "gCF.cf.tree"
    ref_tree.write_text("(A,B);")
    matrix = tmp_path / "msa.fa"
    matrix.write_text(">A\nACGT\n>B\nACGT\n")

    cmd = _build_iqtree_cf_cmd(
        cf_mode="scf",
        executable="iqtree3",
        ref_tree=ref_tree,
        gene_trees=None,
        matrix=matrix,
        scf_quartets=200,
        prefix="sCF",
        threads=8,
    )
    assert cmd == [
        "iqtree3", "-s", str(matrix), "-te", str(ref_tree),
        "--scf", "200", "--prefix", "sCF", "-T", "8",
    ]
```

- [ ] **Step 3: Write test for _build_iqtree_cf_cmd — scfl with model**

```python
def test_build_iqtree_cf_cmd_scfl_with_model(tmp_path: Path) -> None:
    from phyloai.tree.cf import _build_iqtree_cf_cmd

    ref_tree = tmp_path / "gCF.cf.tree"; ref_tree.write_text("(A,B);")
    matrix = tmp_path / "msa.fa"; matrix.write_text(">A\nACGT\n>B\nACGT\n")

    cmd = _build_iqtree_cf_cmd(
        cf_mode="scfl",
        executable="iqtree3",
        ref_tree=ref_tree,
        gene_trees=None,
        matrix=matrix,
        scf_quartets=100,
        prefix="sCFl",
        threads=4,
        model="LG+F+R3",
    )
    assert cmd == [
        "iqtree3", "-s", str(matrix), "-te", str(ref_tree),
        "--scfl", "100", "-m", "LG+F+R3", "--prefix", "sCFl", "-T", "4",
    ]
```

- [ ] **Step 4: Write test for _build_iqtree_cf_cmd — scfl with partitions**

```python
def test_build_iqtree_cf_cmd_scfl_with_partitions(tmp_path: Path) -> None:
    from phyloai.tree.cf import _build_iqtree_cf_cmd

    ref_tree = tmp_path / "gCF.cf.tree"; ref_tree.write_text("(A,B);")
    matrix = tmp_path / "msa.fa"; matrix.write_text(">A\nACGT\n>B\nACGT\n")
    partitions = tmp_path / "msa.best_model.nex"; partitions.write_text("#nexus")

    cmd = _build_iqtree_cf_cmd(
        cf_mode="scfl",
        executable="iqtree3",
        ref_tree=ref_tree,
        gene_trees=None,
        matrix=matrix,
        scf_quartets=100,
        prefix="sCFl",
        threads=4,
        partitions=str(partitions),
    )
    assert cmd == [
        "iqtree3", "-s", str(matrix), "-te", str(ref_tree),
        "--scfl", "100", "-p", str(partitions), "--prefix", "sCFl", "-T", "4",
    ]
```

- [ ] **Step 5: Write test for _build_iqtree_cf_cmd — scfl without model/partitions (auto-model)**

```python
def test_build_iqtree_cf_cmd_scfl_auto_model(tmp_path: Path) -> None:
    from phyloai.tree.cf import _build_iqtree_cf_cmd

    ref_tree = tmp_path / "gCF.cf.tree"; ref_tree.write_text("(A,B);")
    matrix = tmp_path / "msa.fa"; matrix.write_text(">A\nACGT\n>B\nACGT\n")

    cmd = _build_iqtree_cf_cmd(
        cf_mode="scfl",
        executable="iqtree3",
        ref_tree=ref_tree,
        gene_trees=None,
        matrix=matrix,
        scf_quartets=100,
        prefix="sCFl",
        threads=4,
    )
    # No -m or -p — IQ-TREE computes model internally
    assert cmd == [
        "iqtree3", "-s", str(matrix), "-te", str(ref_tree),
        "--scfl", "100", "--prefix", "sCFl", "-T", "4",
    ]
```

- [ ] **Step 6: Write test for _build_iqtree_cf_cmd — gcf+scf combined**

```python
def test_build_iqtree_cf_cmd_gcf_scf_combined(tmp_path: Path) -> None:
    from phyloai.tree.cf import _build_iqtree_cf_cmd

    ref_tree = tmp_path / "species.nwk"; ref_tree.write_text("(A,B);")
    gene_trees = tmp_path / "merged.trees"; gene_trees.write_text("(A,B);")
    matrix = tmp_path / "msa.fa"; matrix.write_text(">A\nACGT\n>B\nACGT\n")

    cmd = _build_iqtree_cf_cmd(
        cf_mode="gcf+scf",
        executable="iqtree3",
        ref_tree=ref_tree,
        gene_trees=gene_trees,
        matrix=matrix,
        scf_quartets=150,
        prefix="gCFsCF",
        threads=4,
    )
    assert cmd == [
        "iqtree3", "-t", str(ref_tree), "--gcf", str(gene_trees),
        "-s", str(matrix), "--scf", "150", "--prefix", "gCFsCF", "-T", "4",
    ]
```

- [ ] **Step 7: Write test for _build_wastral_qcf_cmd**

```python
def test_build_wastral_qcf_cmd(tmp_path: Path) -> None:
    from phyloai.tree.cf import _build_wastral_qcf_cmd

    ref_tree = tmp_path / "species.nwk"; ref_tree.write_text("(A,B);")
    gene_trees = tmp_path / "merged.trees"; gene_trees.write_text("(A,B);")
    output_dir = tmp_path / "out"

    cmd = _build_wastral_qcf_cmd(
        executable="wastral",
        gene_trees=gene_trees,
        ref_tree=ref_tree,
        output_dir=output_dir,
        threads=8,
    )
    assert cmd == [
        "wastral", "-i", str(gene_trees.resolve()),
        "-o", "wastral.tre",
        "-u", "2", "-c", str(ref_tree.resolve()), "-C", "--mode", "4",
        "-t", "8",
    ]
```

- [ ] **Step 8: Run tests to confirm failure**

```bash
pytest tests/tree/test_cf.py -v -k "build"
```

Expected: all FAIL (module not found).

---

### Task 3: Write validation tests for run_cf (TDD)

**Files:**
- Modify: `tests/tree/test_cf.py` (append)

- [ ] **Step 1: Write validation test — scf without matrix**

```python
def test_run_cf_scf_without_matrix_raises(tmp_path: Path) -> None:
    from phyloai.tree.cf import run_cf

    ref_tree = tmp_path / "ref.nwk"; ref_tree.write_text("(A,B);")

    with pytest.raises(ValueError, match="--matrix is required"):
        run_cf(
            cf_mode="scf", ref_tree=ref_tree,
            output_dir=tmp_path / "out",
            threads=4, dry_run=True,
        )
```

- [ ] **Step 2: Write validation test — gcf with matrix**

```python
def test_run_cf_gcf_with_matrix_raises(tmp_path: Path) -> None:
    from phyloai.tree.cf import run_cf

    ref_tree = tmp_path / "ref.nwk"; ref_tree.write_text("(A,B);")
    gene_trees = tmp_path / "trees"; gene_trees.write_text("(A,B);")
    matrix = tmp_path / "msa.fa"; matrix.write_text(">A\nA\n")

    with pytest.raises(ValueError, match="not valid"):
        run_cf(
            cf_mode="gcf", ref_tree=ref_tree,
            tree=gene_trees, matrix=matrix,
            output_dir=tmp_path / "out",
            threads=4, dry_run=True,
        )
```

- [ ] **Step 3: Write validation test — scfl with model AND partitions**

```python
def test_run_cf_scfl_model_and_partitions_mutually_exclusive(tmp_path: Path) -> None:
    from phyloai.tree.cf import run_cf

    ref_tree = tmp_path / "ref.nwk"; ref_tree.write_text("(A,B);")
    matrix = tmp_path / "msa.fa"; matrix.write_text(">A\nA\n")
    partitions = tmp_path / "p.nex"; partitions.write_text("#nexus")

    with pytest.raises(ValueError, match="mutually exclusive"):
        run_cf(
            cf_mode="scfl", ref_tree=ref_tree,
            matrix=matrix, model="LG", partitions=partitions,
            output_dir=tmp_path / "out",
            threads=4, dry_run=True,
        )
```

- [ ] **Step 4: Write validation test — qcf with matrix**

```python
def test_run_cf_qcf_with_matrix_raises(tmp_path: Path) -> None:
    from phyloai.tree.cf import run_cf

    ref_tree = tmp_path / "ref.nwk"; ref_tree.write_text("(A,B);")
    gene_trees = tmp_path / "trees"; gene_trees.write_text("(A,B);")
    matrix = tmp_path / "msa.fa"; matrix.write_text(">A\nA\n")

    with pytest.raises(ValueError, match="not valid"):
        run_cf(
            cf_mode="qcf", ref_tree=ref_tree,
            tree=gene_trees, matrix=matrix,
            output_dir=tmp_path / "out",
            threads=4, dry_run=True,
        )
```

- [ ] **Step 5: Write validation test — scf with tree**

```python
def test_run_cf_scf_with_tree_raises(tmp_path: Path) -> None:
    from phyloai.tree.cf import run_cf

    ref_tree = tmp_path / "ref.nwk"; ref_tree.write_text("(A,B);")
    gene_trees = tmp_path / "trees"; gene_trees.write_text("(A,B);")
    matrix = tmp_path / "msa.fa"; matrix.write_text(">A\nA\n")

    with pytest.raises(ValueError, match="not needed"):
        run_cf(
            cf_mode="scf", ref_tree=ref_tree,
            tree=gene_trees, matrix=matrix,
            output_dir=tmp_path / "out",
            threads=4, dry_run=True,
        )
```

- [ ] **Step 6: Write validation test — gcf with model**

```python
def test_run_cf_gcf_with_model_raises(tmp_path: Path) -> None:
    from phyloai.tree.cf import run_cf

    ref_tree = tmp_path / "ref.nwk"; ref_tree.write_text("(A,B);")
    gene_trees = tmp_path / "trees"; gene_trees.write_text("(A,B);")

    with pytest.raises(ValueError, match="not valid"):
        run_cf(
            cf_mode="gcf", ref_tree=ref_tree,
            tree=gene_trees, model="LG",
            output_dir=tmp_path / "out",
            threads=4, dry_run=True,
        )
```

- [ ] **Step 7: Write validation test — qcf with scf_quartets**

```python
def test_run_cf_qcf_with_scf_quartets_raises(tmp_path: Path) -> None:
    from phyloai.tree.cf import run_cf

    ref_tree = tmp_path / "ref.nwk"; ref_tree.write_text("(A,B);")
    gene_trees = tmp_path / "trees"; gene_trees.write_text("(A,B);")

    with pytest.raises(ValueError, match="not valid"):
        run_cf(
            cf_mode="qcf", ref_tree=ref_tree,
            tree=gene_trees, scf_quartets=200,
            output_dir=tmp_path / "out",
            threads=4, dry_run=True,
        )
```

- [ ] **Step 8: Write test — tree and tree-dir together**

```python
def test_run_cf_tree_and_tree_dir_mutually_exclusive(tmp_path: Path) -> None:
    from phyloai.tree.cf import run_cf

    ref_tree = tmp_path / "ref.nwk"; ref_tree.write_text("(A,B);")
    gene_trees = tmp_path / "trees"; gene_trees.write_text("(A,B);")
    tree_dir = tmp_path / "tdir"; tree_dir.mkdir()

    with pytest.raises(ValueError, match="mutually exclusive"):
        run_cf(
            cf_mode="gcf", ref_tree=ref_tree,
            tree=gene_trees, tree_dir=tree_dir,
            output_dir=tmp_path / "out",
            threads=4, dry_run=True,
        )
```

- [ ] **Step 9: Write test — gcf with neither tree nor tree-dir**

```python
def test_run_cf_gcf_without_gene_trees_raises(tmp_path: Path) -> None:
    from phyloai.tree.cf import run_cf

    ref_tree = tmp_path / "ref.nwk"; ref_tree.write_text("(A,B);")

    with pytest.raises(ValueError, match="--tree or --tree-dir"):
        run_cf(
            cf_mode="gcf", ref_tree=ref_tree,
            output_dir=tmp_path / "out",
            threads=4, dry_run=True,
        )
```

- [ ] **Step 10: Write test — scf_quartets warning when < 100 (only for scf mode)**

```python
def test_run_cf_scf_quartets_below_100_warns(tmp_path: Path) -> None:
    from phyloai.tree.cf import run_cf

    ref_tree = tmp_path / "ref.nwk"; ref_tree.write_text("(A,B);")
    matrix = tmp_path / "msa.fa"; matrix.write_text(">A\nACGT\n>B\nACGT\n")

    result = run_cf(
        cf_mode="scf", ref_tree=ref_tree,
        matrix=matrix, scf_quartets=50,
        output_dir=tmp_path / "out",
        threads=4, dry_run=True,
    )
    assert any(">= 100" in w for w in result.get("data", {}).get("warnings", []))
```

- [ ] **Step 11: Write test — default prefix per mode (each with valid inputs)**

```python
def test_run_cf_default_prefix_per_mode(tmp_path: Path) -> None:
    from phyloai.tree.cf import run_cf

    ref_tree = tmp_path / "ref.nwk"; ref_tree.write_text("(A,B);")
    gene_trees = tmp_path / "trees"; gene_trees.write_text("(A,B);\n")
    matrix = tmp_path / "msa.fa"; matrix.write_text(">A\nACGT\n>B\nACGT\n")

    # Modes needing gene trees
    for mode, expected in [("gcf", "gCF"), ("qcf", "qCF")]:
        result = run_cf(
            cf_mode=mode, ref_tree=ref_tree,
            tree=gene_trees,
            output_dir=tmp_path / "out",
            threads=4, dry_run=True,
        )
        assert result["params"]["prefix"] == expected

    # Modes needing matrix (no gene trees)
    for mode, expected in [("scf", "sCF"), ("scfl", "sCFl")]:
        result = run_cf(
            cf_mode=mode, ref_tree=ref_tree,
            matrix=matrix,
            output_dir=tmp_path / "out",
            threads=4, dry_run=True,
        )
        assert result["params"]["prefix"] == expected

    # gcf+scf needs both
    result = run_cf(
        cf_mode="gcf+scf", ref_tree=ref_tree,
        tree=gene_trees, matrix=matrix,
        output_dir=tmp_path / "out",
        threads=4, dry_run=True,
    )
    assert result["params"]["prefix"] == "gCFsCF"
```

- [ ] **Step 12: Run tests to confirm failure**

```bash
pytest tests/tree/test_cf.py -v -k "run_cf"
```

Expected: all FAIL (module not found).

---

### Task 4: Write dry-run flow tests (TDD)

**Files:**
- Modify: `tests/tree/test_cf.py` (append)

- [ ] **Step 1: Write test — dry_run gcf produces correct payload shape**

```python
def test_run_cf_dry_run_gcf_produces_payload(tmp_path: Path) -> None:
    from phyloai.tree.cf import run_cf

    ref_tree = tmp_path / "ref.nwk"; ref_tree.write_text("(A,B);")
    gene_trees = tmp_path / "trees"; gene_trees.write_text("(A,B);")

    result = run_cf(
        cf_mode="gcf", ref_tree=ref_tree,
        tree=gene_trees,
        output_dir=tmp_path / "out",
        threads=4, dry_run=True,
    )

    assert result["status"] == "success"
    assert result["params"]["cf"] == "gcf"
    assert result["key_results"]["cf_type"] == "gcf"
    assert result["key_results"]["prefix"] == "gCF"
    assert "iqtree3" in result["data"]["cmd"]
    assert "--gcf" in result["data"]["cmd"]
```

- [ ] **Step 2: Write test — dry_run qcf produces correct payload shape**

```python
def test_run_cf_dry_run_qcf_produces_payload(tmp_path: Path) -> None:
    from phyloai.tree.cf import run_cf

    ref_tree = tmp_path / "ref.nwk"; ref_tree.write_text("(A,B);")
    gene_trees = tmp_path / "trees"; gene_trees.write_text("(A,B);")

    result = run_cf(
        cf_mode="qcf", ref_tree=ref_tree,
        tree=gene_trees,
        output_dir=tmp_path / "out",
        threads=4, dry_run=True,
    )

    assert result["status"] == "success"
    assert result["params"]["cf"] == "qcf"
    assert result["key_results"]["cf_type"] == "qcf"
    assert result["key_results"]["prefix"] == "qCF"
    assert "wastral" in result["data"]["cmd"]
    assert "-u" in result["data"]["cmd"] and "2" in result["data"]["cmd"]
    assert "-C" in result["data"]["cmd"]
    assert "--mode" in result["data"]["cmd"] and "4" in result["data"]["cmd"]
```

- [ ] **Step 3: Write test — dry_run scfl with model produces correct command**

```python
def test_run_cf_dry_run_scfl_with_model(tmp_path: Path) -> None:
    from phyloai.tree.cf import run_cf

    ref_tree = tmp_path / "ref.nwk"; ref_tree.write_text("(A,B);")
    matrix = tmp_path / "msa.fa"; matrix.write_text(">A\nACGT\n>B\nACGT\n")

    result = run_cf(
        cf_mode="scfl", ref_tree=ref_tree,
        matrix=matrix, model="LG+F+R3",
        output_dir=tmp_path / "out",
        threads=4, dry_run=True,
    )

    cmd = result["data"]["cmd"]
    assert "--scfl" in cmd
    assert "-m" in cmd and "LG+F+R3" in cmd
    assert "-te" in cmd
```

- [ ] **Step 4: Write test — dry_run scfl with partitions produces correct command**

```python
def test_run_cf_dry_run_scfl_with_partitions(tmp_path: Path) -> None:
    from phyloai.tree.cf import run_cf

    ref_tree = tmp_path / "ref.nwk"; ref_tree.write_text("(A,B);")
    matrix = tmp_path / "msa.fa"; matrix.write_text(">A\nACGT\n>B\nACGT\n")
    partitions = tmp_path / "best_model.nex"; partitions.write_text("#nexus")

    result = run_cf(
        cf_mode="scfl", ref_tree=ref_tree,
        matrix=matrix, partitions=partitions,
        output_dir=tmp_path / "out",
        threads=4, dry_run=True,
    )

    cmd = result["data"]["cmd"]
    assert "-p" in cmd
    assert str(partitions) in cmd
```

- [ ] **Step 5: Write test — dry_run with explicit prefix overrides default**

```python
def test_run_cf_explicit_prefix_overrides_default(tmp_path: Path) -> None:
    from phyloai.tree.cf import run_cf

    ref_tree = tmp_path / "ref.nwk"; ref_tree.write_text("(A,B);")
    gene_trees = tmp_path / "trees"; gene_trees.write_text("(A,B);")

    result = run_cf(
        cf_mode="gcf", ref_tree=ref_tree,
        tree=gene_trees, prefix="myCF",
        output_dir=tmp_path / "out",
        threads=4, dry_run=True,
    )

    assert result["params"]["prefix"] == "myCF"
    assert result["key_results"]["prefix"] == "myCF"
```

- [ ] **Step 6: Write test — dry_run with --tree-dir**

```python
def test_run_cf_dry_run_tree_dir_merges_trees(tmp_path: Path) -> None:
    from phyloai.tree.cf import run_cf

    ref_tree = tmp_path / "ref.nwk"; ref_tree.write_text("(A,B);")
    td = tmp_path / "tdir"; td.mkdir()
    (td / "gene1.nwk").write_text("(A,B);\n")
    (td / "gene2.tre").write_text("(B,A);\n")

    result = run_cf(
        cf_mode="gcf", ref_tree=ref_tree,
        tree_dir=td,
        output_dir=tmp_path / "out",
        threads=4, dry_run=True,
    )

    assert result["data"]["input_mode"] == "--tree-dir"
    assert result["data"]["input"]["n_trees"] == 2
```

- [ ] **Step 7: Write test — dry_run tree_dir with 0 valid files raises**

```python
def test_run_cf_tree_dir_zero_valid_files_raises(tmp_path: Path) -> None:
    from phyloai.tree.cf import run_cf

    ref_tree = tmp_path / "ref.nwk"; ref_tree.write_text("(A,B);")
    td = tmp_path / "tdir"; td.mkdir()
    (td / "data.txt").write_text("not a tree")

    with pytest.raises(ValueError, match="No valid gene tree files"):
        run_cf(
            cf_mode="gcf", ref_tree=ref_tree,
            tree_dir=td,
            output_dir=tmp_path / "out",
            threads=4, dry_run=True,
        )
```

- [ ] **Step 8: Write test — dry_run tree_dir with 1 valid file warns**

```python
def test_run_cf_tree_dir_one_file_warns(tmp_path: Path) -> None:
    from phyloai.tree.cf import run_cf

    ref_tree = tmp_path / "ref.nwk"; ref_tree.write_text("(A,B);")
    td = tmp_path / "tdir"; td.mkdir()
    (td / "gene1.nwk").write_text("(A,B);\n")

    result = run_cf(
        cf_mode="gcf", ref_tree=ref_tree,
        tree_dir=td,
        output_dir=tmp_path / "out",
        threads=4, dry_run=True,
    )

    warnings = result.get("data", {}).get("warnings", [])
    assert any("--tree" in w for w in warnings)
```

- [ ] **Step 9: Verify tests fail**

```bash
pytest tests/tree/test_cf.py -v -k "dry_run"
```

Expected: all FAIL.

---

### Task 5: Create cf.py skeleton with constants and _merge_gene_trees

**Files:**
- Create: `phyloai/tree/cf.py`

- [ ] **Step 1: Write module skeleton**

```python
"""Concordance factor computation (gCF, sCF, sCFl, qCF) via IQ-TREE3 and wASTRAL."""

from __future__ import annotations

import json
import os
import re as _re
import shlex
import shutil
import subprocess
import time as _time
from pathlib import Path
from typing import Any

from phyloai.core.env import ToolEnv

# Gene tree file extensions (newick variants)
_CF_TREE_EXTENSIONS = frozenset({
    ".nwk", ".tre", ".tree", ".nw", ".trees", ".newick",
})

# Default prefix per CF mode
_DEFAULT_PREFIX: dict[str, str] = {
    "gcf": "gCF",
    "scf": "sCF",
    "scfl": "sCFl",
    "gcf+scf": "gCFsCF",
    "qcf": "qCF",
}

# CF modes that require gene trees
_CF_MODES_NEED_GENE_TREES = frozenset({"gcf", "gcf+scf", "qcf"})

# CF modes that require a matrix
_CF_MODES_NEED_MATRIX = frozenset({"scf", "scfl", "gcf+scf"})

# CF modes that use iqtree3
_CF_MODES_IQTREE = frozenset({"gcf", "scf", "scfl", "gcf+scf"})

# Parameters only valid for scfl mode
_SCFL_ONLY_PARAMS = frozenset({"model", "partitions"})

# Parameters that MUST NOT be set for non-IQ-TREE modes
# Each tuple: (param_name, error_fragment_for_match)
_NON_IQTREE_BLOCKED_PARAMS = frozenset({"model", "partitions"})
```

- [ ] **Step 2: Write _scan_input_cf**

```python
def _scan_input_cf(
    tree_dir: Path,
) -> tuple[list[Path], list[dict[str, str]]]:
    """Scan a directory for valid gene tree files.

    Returns:
        (valid_files, skipped_entries)
    """
    if not tree_dir.exists():
        return [], []

    found: list[Path] = []
    skipped: list[dict[str, str]] = []

    for entry in sorted(tree_dir.iterdir()):
        if entry.is_dir():
            skipped.append({"path": str(entry), "reason": "directory"})
            continue
        if not entry.is_file():
            skipped.append({"path": str(entry), "reason": "not a regular file"})
            continue
        if entry.stat().st_size == 0:
            skipped.append({"path": str(entry), "reason": "empty file"})
            continue

        ext = entry.suffix.lower()
        if ext in _CF_TREE_EXTENSIONS:
            found.append(entry)
        else:
            skipped.append({"path": str(entry), "reason": f"unrecognized extension: {ext}"})

    return found, skipped
```

- [ ] **Step 3: Write _merge_gene_trees**

```python
def _merge_gene_trees(
    tree_dir: Path,
    output_path: Path,
) -> tuple[int, list[dict[str, str]]]:
    """Scan tree_dir for newick files, merge into one file (one tree per line).

    Returns:
        (count_of_trees_merged, skipped_entries)
    """
    found, skipped = _scan_input_cf(tree_dir)

    count = 0
    with open(output_path, "w") as out:
        for f in found:
            content = f.read_text().strip()
            if not content:
                continue
            for line in content.splitlines():
                line = line.strip()
                if line:
                    out.write(line + "\n")
                    count += 1

    return count, skipped
```

- [ ] **Step 4: Run tests to see which pass**

```bash
pytest tests/tree/test_cf.py -v
```

Expected: `_map_qcf_to_tree` and `_build_*` tests still fail (not implemented). Validation tests may partially pass.

---

### Task 6: Implement _map_qcf_to_tree

**Files:**
- Modify: `phyloai/tree/cf.py` (append)

- [ ] **Step 1: Write _map_qcf_to_tree with clade.name storage and raw-Newick fallback**

```python
def _map_qcf_to_tree(
    ref_tree_path: Path,
    wastral_output_path: Path,
    output_path: Path,
) -> None:
    """Map qCF (q1) values from wastral tree onto reference tree topology.

    For each internal node in the reference tree, finds the corresponding
    bipartition in the wastral tree (comparing leaf set vs complement with
    canonical min-form), extracts q1, multiplies by 100, rounds to integer,
    and appends after existing support value separated by '/'.

    Support strings are stored in clade.name (not clade.confidence) to avoid
    Bio.Phylo serialization issues with non-numeric confidence values.
    """
    from Bio import Phylo
    from io import StringIO

    # Parse both trees via raw text first for label extraction fallback
    ref_raw = ref_tree_path.read_text().strip()
    wastral_raw = wastral_output_path.read_text().strip()

    ref_tree = Phylo.read(StringIO(ref_raw), "newick")
    wastral_tree = Phylo.read(StringIO(wastral_raw), "newick")

    # --- Raw Newick q1 extraction (fallback if Bio.Phylo strips labels) ---
    def _extract_raw_q1_map(raw_nwk: str, all_leaves: frozenset[str]) -> dict[frozenset[str], float]:
        """Parse raw Newick to extract q1 per bipartition by matching node order.

        Strategy: collect leaf names via a stack as we walk the string, then at each
        ')' look ahead for the label. Build a list of (leaf_set, label) pairs in
        post-order. Map these to Bio.Phylo nodes by matching the post-order position.
        """
        result: dict[frozenset[str], float] = {}
        q1_pattern = _re.compile(r"q1=([0-9.]+)")

        # First pass: collect all terminal labels
        leaf_names: list[str] = []
        i = 0
        while i < len(raw_nwk):
            c = raw_nwk[i]
            if c == '(':
                i += 1
                continue
            if c in (')', ',', ';'):
                i += 1
                continue
            if c == ':':
                # Skip branch length
                i += 1
                while i < len(raw_nwk) and raw_nwk[i] not in (',', ')', ';'):
                    i += 1
                continue
            # Leaf name — read until : , ) or ;
            start = i
            while i < len(raw_nwk) and raw_nwk[i] not in (':', ',', ')', ';'):
                i += 1
            name = raw_nwk[start:i].strip()
            if name:
                leaf_names.append(name)
            if i < len(raw_nwk) and raw_nwk[i] == ':':
                i += 1
                while i < len(raw_nwk) and raw_nwk[i] not in (',', ')', ';'):
                    i += 1
            continue

        # Second pass: post-order walk collecting internal node labels
        # Use a simpler strategy: split on ')', look at label before ':'
        # and map to wASTRAL leaf sets via Bio.Phylo tree topology
        # Since wastral tree and ref tree share topology, we skip raw leaf-set
        # collection and instead use the Bio.Phylo wastral tree to get leaf sets,
        # then try raw label extraction for each wastral clade by position

        # Build post-order list of wastral clades
        wastral_clades_postorder = [
            c for c in wastral_tree.find_clades(order="postorder")
            if not c.is_terminal()
        ]

        # Extract internal node labels from raw Newick in post-order
        # Each ')' marks an internal node; its label precedes ':'
        internal_labels: list[str] = []
        i = 0
        while i < len(raw_nwk):
            if raw_nwk[i] == ')':
                i += 1
                label_start = i
                while i < len(raw_nwk) and raw_nwk[i] not in (':', ',', ';', ')'):
                    i += 1
                label = raw_nwk[label_start:i].strip()
                internal_labels.append(label)
                if i < len(raw_nwk) and raw_nwk[i] == ':':
                    i += 1
                    while i < len(raw_nwk) and raw_nwk[i] not in (',', ';', ')'):
                        i += 1
                continue
            i += 1

        # Match post-order clades with extracted labels
        for j, clade in enumerate(wastral_clades_postorder):
            if j >= len(internal_labels):
                break
            label = internal_labels[j]
            m = q1_pattern.search(label)
            if m:
                ls = _leaf_set(clade)
                complement = all_leaves - ls
                canonical = ls if sorted(ls) < sorted(complement) else complement
                result[canonical] = float(m.group(1))

        return result

    # Collect all leaf names from Bio.Phylo tree
    all_leaves = frozenset(
        clade.name for clade in ref_tree.find_clades()
        if clade.is_terminal() and clade.name is not None
    )

    # Try raw Newick extraction first
    raw_q1_map = _extract_raw_q1_map(wastral_raw, all_leaves)

    def _leaf_set(clade):
        return frozenset(
            t.name for t in clade.find_clades(terminal=True)
            if t.name is not None
        )

    # Build canonical bipartition -> wastral clade map from Bio.Phylo
    wastral_bip_map: dict[frozenset[str], Any] = {}
    for clade in wastral_tree.find_clades():
        if clade.is_terminal():
            continue
        ls = _leaf_set(clade)
        complement = all_leaves - ls
        canonical = ls if sorted(ls) < sorted(complement) else complement
        wastral_bip_map[canonical] = clade

    # Extract q1 from Bio.Phylo clade (name / comment)
    def _extract_q1(clade) -> float | None:
        for attr in (clade.name, getattr(clade, "comment", None)):
            if attr and isinstance(attr, str) and "q1=" in attr:
                m = _re.search(r"q1=([0-9.]+)", attr)
                if m:
                    return float(m.group(1))
        return None

    # Walk reference tree and annotate
    for clade in ref_tree.find_clades():
        if clade.is_terminal():
            continue
        ls = _leaf_set(clade)
        complement = all_leaves - ls
        canonical = ls if sorted(ls) < sorted(complement) else complement

        # Try Bio.Phylo match first, then raw Newick fallback
        q1 = None
        matched = wastral_bip_map.get(canonical)
        if matched is not None:
            q1 = _extract_q1(matched)
        if q1 is None:
            q1 = raw_q1_map.get(canonical)

        if q1 is None:
            continue

        qcf_str = str(round(q1 * 100))

        # Existing support: clade.name (if string label) or clade.confidence (if numeric)
        existing = None
        if clade.name is not None:
            existing = clade.name
        elif clade.confidence is not None:
            # Render numeric support without trailing .0 for integer values
            conf = clade.confidence
            if isinstance(conf, float) and conf == int(conf):
                existing = str(int(conf))
            else:
                existing = str(conf)

        new_label = f"{existing}/{qcf_str}" if existing else qcf_str

        # Store as clade.name (string) to avoid Bio.Phylo confidence serialization issues
        clade.name = new_label
        clade.confidence = None

    Phylo.write(ref_tree, str(output_path), "newick")
```

- [ ] **Step 2: Run qCF mapper tests**

```bash
pytest tests/tree/test_cf.py -v -k "map_qcf"
```

Expect some to pass, some may need adjustment. If Bio.Phylo strips the `[q1=...]` label, adjust `_extract_q1` to also parse raw Newick string directly with regex on the wastral_raw text.

If Bio.Phylo fails to parse wastral's quoted labels, replace the wastral tree parsing with a regex-based approach that extracts q1 per node using the raw Newick string alongside Bio.Phylo for topology comparison only.

- [ ] **Step 3: Debug and fix until all 6 qCF mapper tests pass**

```bash
pytest tests/tree/test_cf.py -v -k "map_qcf"
```

Expected: all 6 PASS.

- [ ] **Step 4: Commit**

```bash
git add phyloai/tree/cf.py tests/tree/test_cf.py
git commit -m "feat(tree/cf): add _map_qcf_to_tree with bipartition matching"
```

---

### Task 7: Implement command builders

**Files:**
- Modify: `phyloai/tree/cf.py` (append)

- [ ] **Step 1: Write _build_iqtree_cf_cmd**

```python
def _build_iqtree_cf_cmd(
    *,
    cf_mode: str,
    executable: str,
    ref_tree: Path,
    gene_trees: Path | None,
    matrix: Path | None,
    scf_quartets: int,
    prefix: str,
    threads: int,
    model: str | None = None,
    partitions: str | None = None,
) -> list[str]:
    """Build the IQ-TREE3 command for a given CF mode."""
    cmd = [executable]

    if cf_mode in ("gcf", "gcf+scf"):
        cmd.extend(["-t", str(ref_tree)])
    else:
        # scf, scfl: use -te (fixed tree topology)
        cmd.extend(["-te", str(ref_tree)])

    if cf_mode in ("gcf", "gcf+scf"):
        assert gene_trees is not None
        cmd.extend(["--gcf", str(gene_trees)])

    if cf_mode in ("scf", "scfl", "gcf+scf"):
        assert matrix is not None
        cmd.extend(["-s", str(matrix)])

    if cf_mode in ("scf", "gcf+scf"):
        cmd.extend(["--scf", str(scf_quartets)])
    elif cf_mode == "scfl":
        cmd.extend(["--scfl", str(scf_quartets)])

    if cf_mode == "scfl":
        if partitions is not None:
            cmd.extend(["-p", partitions])
        elif model is not None:
            cmd.extend(["-m", model])

    cmd.extend(["--prefix", prefix])
    cmd.extend(["-T", str(threads)])

    return cmd
```

- [ ] **Step 2: Write _build_wastral_qcf_cmd**

```python
def _build_wastral_qcf_cmd(
    *,
    executable: str,
    gene_trees: Path,
    ref_tree: Path,
    output_dir: Path,   # cwd for the subprocess
    threads: int,
) -> list[str]:
    """Build the wASTRAL command for qCF computation.

    Uses absolute paths for all inputs; output uses 'wastral.tre' (relative
    to cwd) since the subprocess runs with cwd=output_dir.
    """
    return [
        executable,
        "-i", str(gene_trees.resolve()),
        "-o", "wastral.tre",
        "-u", "2",
        "-c", str(ref_tree.resolve()),
        "-C",
        "--mode", "4",
        "-t", str(threads),
    ]
```

- [ ] **Step 3: Run command builder tests**

```bash
pytest tests/tree/test_cf.py -v -k "build"
```

Expected: all 7 PASS.

- [ ] **Step 4: Commit**

```bash
git add phyloai/tree/cf.py
git commit -m "feat(tree/cf): add command builders for iqtree3 and wastral"
```

---

### Task 8: Implement tool resolution and version detection

**Files:**
- Modify: `phyloai/tree/cf.py` (append)

- [ ] **Step 1: Write _resolve_iqtree_path**

```python
def _resolve_iqtree_path(iqtree_path: str | None, dry_run: bool) -> str:
    """Resolve iqtree3 executable path."""
    if iqtree_path:
        p = Path(iqtree_path)
        if not p.exists():
            raise ValueError(f"--iqtree-path does not exist: {iqtree_path}")
        if not os.access(str(p), os.X_OK):
            raise ValueError(f"--iqtree-path is not executable: {iqtree_path}")
        return iqtree_path
    if dry_run:
        return "iqtree3"
    try:
        env = ToolEnv()
        return str(env.require("iqtree3"))
    except FileNotFoundError:
        raise FileNotFoundError(
            "iqtree3 not found. Install from https://github.com/iqtree/iqtree3/releases "
            "or use --iqtree-path."
        )
```

- [ ] **Step 2: Write _detect_iqtree_version**

```python
def _detect_iqtree_version(executable: str) -> dict[str, str]:
    """Detect iqtree3 version via --version."""
    try:
        proc = subprocess.run(
            [executable, "--version"],
            capture_output=True, text=True, timeout=10,
        )
        output = proc.stdout + proc.stderr
        m = _re.search(r"version\s+(\S+)", output)
        if m:
            return {"iqtree3": m.group(1)}
    except Exception:
        pass
    return {"iqtree3": "unknown"}
```

- [ ] **Step 3: Write _resolve_wastral_path**

```python
def _resolve_wastral_path(wastral_path: str | None, dry_run: bool) -> str:
    """Resolve wastral executable path."""
    if wastral_path:
        p = Path(wastral_path)
        if not p.exists():
            raise ValueError(f"--wastral-path does not exist: {wastral_path}")
        if not os.access(str(p), os.X_OK):
            raise ValueError(f"--wastral-path is not executable: {wastral_path}")
        return wastral_path
    if dry_run:
        return "wastral"
    for name in ("wastral", "aster"):
        env = ToolEnv()
        try:
            return str(env.require(name))
        except FileNotFoundError:
            continue
    raise FileNotFoundError(
        "wastral not found. Install from https://github.com/chaoszhang/ASTER "
        "or use --wastral-path."
    )
```

- [ ] **Step 4: Write _detect_wastral_version**

```python
def _detect_wastral_version(executable: str) -> dict[str, str]:
    """Detect wastral version via -v or -h fallback."""
    combined = ""
    for flag in ("-v", "-h"):
        try:
            proc = subprocess.run(
                [executable, flag],
                capture_output=True, text=True, timeout=10,
            )
            combined = proc.stdout + "\n" + proc.stderr
            if proc.returncode == 0 or (proc.stdout.strip() or proc.stderr.strip()):
                break
        except Exception:
            continue

    # Try to extract version from output
    m = _re.search(r"version\s+([\d.]+)", combined, _re.IGNORECASE)
    if m:
        return {"wastral": m.group(1)}
    m = _re.search(r"ASTER[-\s]+([\d.]+)", combined, _re.IGNORECASE)
    if m:
        return {"wastral": m.group(1)}
    return {"wastral": "unknown"}
```

- [ ] **Step 5: Commit**

```bash
git add phyloai/tree/cf.py
git commit -m "feat(tree/cf): add tool resolution and version detection"
```

---

### Task 9: Implement run_cf entry point

**Files:**
- Modify: `phyloai/tree/cf.py` (append)

- [ ] **Step 1: Write _assemble_cf_result**

```python
def _assemble_cf_result(
    *,
    run_start: float,
    cf_mode: str,
    ref_tree: Path,
    tree: Path | None,
    tree_dir: Path | None,
    matrix: Path | None,
    partitions: Path | None,
    model: str | None,
    scf_quartets: int,
    prefix: str,
    output_dir: Path,
    threads: int,
    iqtree_path: str | None,
    wastral_path: str | None,
    overwrite: bool,
    dry_run: bool,
    input_path: Path,
    n_input_trees: int,
    cmd: list[str],
    wall_time: float,
    skipped: list[dict[str, str]],
    warnings_list: list[str],
    is_error: bool,
    error_msg: str | None,
    versions: dict[str, str],
) -> dict[str, Any]:
    """Build the result.json payload."""
    if tree_dir is not None:
        input_mode = "--tree-dir"
    elif tree is not None:
        input_mode = "--tree"
    else:
        input_mode = "--matrix"

    # Reconstruct CLI invocation string
    cmd_parts = ["phyloai", "tree", "cf", "--cf", cf_mode]
    cmd_parts.extend(["--ref-tree", str(ref_tree)])
    if tree is not None:
        cmd_parts.extend(["--tree", str(tree)])
    elif tree_dir is not None:
        cmd_parts.extend(["--tree-dir", str(tree_dir)])
    if matrix is not None:
        cmd_parts.extend(["--matrix", str(matrix)])
    if partitions is not None:
        cmd_parts.extend(["--partitions", str(partitions)])
    if model is not None:
        cmd_parts.extend(["--model", model])
    if cf_mode not in ("gcf", "qcf"):
        cmd_parts.extend(["--scf-quartets", str(scf_quartets)])
    cmd_parts.extend(["--prefix", prefix])
    cmd_parts.extend(["-o", str(output_dir)])
    cmd_parts.extend(["-t", str(threads)])
    if overwrite:
        cmd_parts.append("--overwrite")
    cmd_str = " ".join(cmd_parts)

    input_data: dict[str, Any] = {"path": str(input_path)}
    if input_mode == "--tree-dir":
        input_data["n_trees"] = n_input_trees

    key_results: dict[str, Any] = {
        "cf_type": cf_mode,
        "prefix": prefix,
    }
    if cf_mode in _CF_MODES_IQTREE:
        key_results["cf_stat"] = str(output_dir / f"{prefix}.cf.stat")
        key_results["cf_tree"] = str(output_dir / f"{prefix}.cf.tree")
    else:
        # qcf
        key_results["cf_tree"] = str(output_dir / f"{prefix}.cf.tree")
        key_results["wastral_log"] = str(output_dir / "wastral.log")

    return {
        "status": "error" if is_error else "success",
        "command": cmd_str,
        "wall_time": _time.monotonic() - run_start,
        "tool_versions": versions,
        "params": {
            "cf": cf_mode,
            "ref_tree": str(ref_tree),
            "tree": str(tree) if tree else None,
            "tree_dir": str(tree_dir) if tree_dir else None,
            "matrix": str(matrix) if matrix else None,
            "partitions": str(partitions) if partitions else None,
            "model": model,
            "scf_quartets": scf_quartets if cf_mode not in ("gcf", "qcf") else None,
            "prefix": prefix,
            "output_dir": str(output_dir),
            "threads": threads,
            "overwrite": overwrite,
            "dry_run": dry_run,
            "iqtree_path": iqtree_path,
            "wastral_path": wastral_path,
        },
        "key_results": key_results,
        "error": error_msg,
        "data": {
            "input_mode": input_mode,
            "input": input_data,
            "cmd": cmd,
            "skipped": skipped,
            "warnings": warnings_list,
        },
    }
```

- [ ] **Step 2: Write run_cf**

```python
def run_cf(
    *,
    cf_mode: str,
    ref_tree: Path,
    tree: Path | None = None,
    tree_dir: Path | None = None,
    matrix: Path | None = None,
    partitions: Path | None = None,
    model: str | None = None,
    scf_quartets: int = 100,
    prefix: str | None = None,
    output_dir: Path = Path("runs/tree/cf"),
    threads: int = 4,
    iqtree_path: str | None = None,
    wastral_path: str | None = None,
    overwrite: bool = False,
    dry_run: bool = False,
    quiet: bool = False,
) -> dict[str, Any]:
    """Run concordance factor computation.

    Returns a result.json-compatible payload dict.
    Raises ValueError for invalid inputs.
    Raises FileNotFoundError for missing tools.
    """
    # --- Validate cf_mode ---
    valid_modes = frozenset({"gcf", "scf", "scfl", "gcf+scf", "qcf"})
    if cf_mode not in valid_modes:
        raise ValueError(
            f"Invalid --cf mode: {cf_mode}. Valid: {', '.join(sorted(valid_modes))}"
        )

    # --- Resolve prefix ---
    if prefix is None:
        prefix = _DEFAULT_PREFIX[cf_mode]

    # --- Validate ref_tree ---
    if not ref_tree.exists():
        raise ValueError(f"--ref-tree does not exist: {ref_tree}")

    # --- Validate gene trees (tree / tree-dir) ---
    needs_gene_trees = cf_mode in _CF_MODES_NEED_GENE_TREES
    if needs_gene_trees:
        if (tree is None and tree_dir is None) or (tree is not None and tree_dir is not None):
            raise ValueError(
                "--tree or --tree-dir must be provided (mutually exclusive) "
                f"for --cf {cf_mode}."
            )
        if tree is not None and not tree.exists():
            raise ValueError(f"--tree does not exist: {tree}")
        if tree_dir is not None and not tree_dir.exists():
            raise ValueError(f"--tree-dir does not exist: {tree_dir}")
    else:
        # scf, scfl: gene trees not needed
        if tree is not None or tree_dir is not None:
            raise ValueError(
                f"--tree/--tree-dir is not needed for --cf {cf_mode}."
            )

    # --- Validate matrix ---
    needs_matrix = cf_mode in _CF_MODES_NEED_MATRIX
    if needs_matrix:
        if matrix is None:
            raise ValueError(f"--matrix is required for --cf {cf_mode}.")
        if not matrix.exists():
            raise ValueError(f"--matrix does not exist: {matrix}")
    else:
        if matrix is not None:
            raise ValueError(
                f"--matrix is not valid for --cf {cf_mode}."
            )

    # --- Validate scfl-only params ---
    if cf_mode != "scfl":
        if model is not None:
            raise ValueError(f"--model is not valid for --cf {cf_mode}.")
        if partitions is not None:
            raise ValueError(f"--partitions is not valid for --cf {cf_mode}.")

    if cf_mode == "scfl":
        if model is not None and partitions is not None:
            raise ValueError(
                "--model and --partitions are mutually exclusive for --cf scfl."
            )

    # --- Validate scf_quartets ---
    if cf_mode in ("gcf", "qcf"):
        if scf_quartets != 100:  # user explicitly set it
            raise ValueError(
                f"--scf-quartets is not valid for --cf {cf_mode}."
            )
    else:
        if scf_quartets < 1:
            raise ValueError(f"--scf-quartets must be >= 1. Got: {scf_quartets}")

    if threads < 1:
        raise ValueError(f"--threads must be >= 1. Got: {threads}")

    # --- Output directory conflict ---
    if not dry_run:
        if overwrite and output_dir.exists():
            shutil.rmtree(output_dir)
        if not overwrite and output_dir.exists() and any(output_dir.iterdir()):
            raise ValueError(
                f"Output directory {output_dir} already exists and is non-empty. "
                "Use --overwrite to replace."
            )

    run_start = _time.monotonic()
    skipped: list[dict[str, str]] = []
    warnings_list: list[str] = []
    n_input_trees = 0

    # --- Input resolution ---
    if tree is not None:
        input_path = tree
    elif tree_dir is not None:
        valid_files, scanned_skipped = _scan_input_cf(tree_dir)
        skipped = scanned_skipped
        n_valid_files = len(valid_files)

        if n_valid_files == 0:
            raise ValueError(
                f"No valid gene tree files found in --tree-dir: {tree_dir}"
            )

        if n_valid_files == 1 and not quiet:
            warnings_list.append(
                "Exactly 1 valid gene tree file in --tree-dir. "
                "Consider using --tree mode directly."
            )

        unrecognized = [
            s for s in skipped
            if s.get("reason", "").startswith("unrecognized")
        ]
        if unrecognized and not quiet:
            warnings_list.append(
                f"--tree-dir contains {len(unrecognized)} non-newick file(s); "
                "skipped. See result.json data.skipped for details."
            )

        if dry_run:
            merged_path = output_dir / "merged.trees"
            n_input_trees = 0
            for f in valid_files:
                content = f.read_text().strip()
                if content:
                    n_input_trees += len([l for l in content.splitlines() if l.strip()])
            input_path = merged_path
        else:
            output_dir.mkdir(parents=True, exist_ok=True)
            merged_path = output_dir / "merged.trees"
            n_input_trees, _ = _merge_gene_trees(tree_dir, merged_path)
            input_path = merged_path
    else:
        # scf/scfl modes: no gene tree input
        input_path = matrix  # type: ignore[assignment]

    # --- Warnings ---
    if cf_mode not in ("gcf", "qcf") and scf_quartets < 100 and not quiet:
        warnings_list.append(
            f"--scf-quartets is {scf_quartets}; recommend >= 100 for reliable results."
        )

    if cf_mode == "scfl" and model is None and partitions is None and not quiet:
        warnings_list.append(
            "--cf scfl without --model or --partitions: IQ-TREE3 will "
            "auto-compute the best-fit model (slow). Consider providing "
            "--model or --partitions for speedup."
        )

    # --- Resolve executables ---
    if cf_mode in _CF_MODES_IQTREE:
        iqtree_exe = _resolve_iqtree_path(iqtree_path, dry_run)
    else:
        iqtree_exe = None
        wastral_exe = _resolve_wastral_path(wastral_path, dry_run)

    # --- Command building ---
    if cf_mode in _CF_MODES_IQTREE:
        cmd = _build_iqtree_cf_cmd(
            cf_mode=cf_mode,
            executable=iqtree_exe,  # type: ignore[arg-type]
            ref_tree=ref_tree,
            gene_trees=input_path if cf_mode in ("gcf", "gcf+scf") else None,
            matrix=matrix,
            scf_quartets=scf_quartets,
            prefix=prefix,
            threads=threads,
            model=model,
            partitions=str(partitions) if partitions else None,
        )
    else:
        cmd = _build_wastral_qcf_cmd(
            executable=wastral_exe,  # type: ignore[arg-type]
            gene_trees=input_path,
            ref_tree=ref_tree,
            output_dir=output_dir,
            threads=threads,
        )

    # --- Dry run: return payload without execution ---
    if dry_run:
        versions = {}
        if cf_mode in _CF_MODES_IQTREE:
            versions = {"iqtree3": "unknown"}
        else:
            versions = {"wastral": "unknown"}
        return _assemble_cf_result(
            run_start=run_start,
            cf_mode=cf_mode, ref_tree=ref_tree,
            tree=tree, tree_dir=tree_dir,
            matrix=matrix, partitions=partitions,
            model=model, scf_quartets=scf_quartets,
            prefix=prefix, output_dir=output_dir,
            threads=threads,
            iqtree_path=iqtree_path, wastral_path=wastral_path,
            overwrite=overwrite, dry_run=dry_run,
            input_path=input_path, n_input_trees=n_input_trees,
            cmd=cmd, wall_time=0.0,
            skipped=skipped, warnings_list=warnings_list,
            is_error=False, error_msg=None,
            versions=versions,
        )

    # --- Execution ---
    output_dir.mkdir(parents=True, exist_ok=True)

    versions: dict[str, str] = {}
    try:
        proc = subprocess.run(
            cmd,
            cwd=str(output_dir),
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
    except Exception as exc:
        if cf_mode in _CF_MODES_IQTREE:
            versions = {"iqtree3": "unknown"}
        else:
            versions = {"wastral": "unknown"}
        return _assemble_cf_result(
            run_start=run_start,
            cf_mode=cf_mode, ref_tree=ref_tree,
            tree=tree, tree_dir=tree_dir,
            matrix=matrix, partitions=partitions,
            model=model, scf_quartets=scf_quartets,
            prefix=prefix, output_dir=output_dir,
            threads=threads,
            iqtree_path=iqtree_path, wastral_path=wastral_path,
            overwrite=overwrite, dry_run=dry_run,
            input_path=input_path, n_input_trees=n_input_trees,
            cmd=cmd, wall_time=0.0,
            skipped=skipped, warnings_list=warnings_list,
            is_error=True, error_msg=str(exc),
            versions=versions,
        )

    if cf_mode in _CF_MODES_IQTREE:
        versions = _detect_iqtree_version(iqtree_exe)  # type: ignore[arg-type]
    else:
        versions = _detect_wastral_version(wastral_exe)  # type: ignore[arg-type]

    wall_time = _time.monotonic() - run_start

    # For qCF mode, always save wastral log before checking return code (diagnostics)
    if cf_mode == "qcf":
        (output_dir / "wastral.log").write_text(proc.stderr + "\n" + proc.stdout)

    if proc.returncode != 0:
        error_msg = (
            f"IQ-TREE3 exited with code {proc.returncode}"
            if cf_mode in _CF_MODES_IQTREE
            else f"wASTRAL exited with code {proc.returncode}"
        )
        error_msg += f": {proc.stderr[:500]}"
        return _assemble_cf_result(
            run_start=run_start,
            cf_mode=cf_mode, ref_tree=ref_tree,
            tree=tree, tree_dir=tree_dir,
            matrix=matrix, partitions=partitions,
            model=model, scf_quartets=scf_quartets,
            prefix=prefix, output_dir=output_dir,
            threads=threads,
            iqtree_path=iqtree_path, wastral_path=wastral_path,
            overwrite=overwrite, dry_run=dry_run,
            input_path=input_path, n_input_trees=n_input_trees,
            cmd=cmd, wall_time=wall_time,
            skipped=skipped, warnings_list=warnings_list,
            is_error=True, error_msg=error_msg,
            versions=versions,
        )

    # Post-process for qCF: map values from wastral.tre to ref tree
    if cf_mode == "qcf":
        wastral_tre_path = output_dir / "wastral.tre"

        if wastral_tre_path.exists():
            try:
                _map_qcf_to_tree(ref_tree, wastral_tre_path, output_dir / f"{prefix}.cf.tree")
            except Exception as exc:
                # qCF mapping failure is fatal — primary output not produced
                error_msg = f"qCF mapping failed: {exc}"
                return _assemble_cf_result(
                    run_start=run_start,
                    cf_mode=cf_mode, ref_tree=ref_tree,
                    tree=tree, tree_dir=tree_dir,
                    matrix=matrix, partitions=partitions,
                    model=model, scf_quartets=scf_quartets,
                    prefix=prefix, output_dir=output_dir,
                    threads=threads,
                    iqtree_path=iqtree_path, wastral_path=wastral_path,
                    overwrite=overwrite, dry_run=dry_run,
                    input_path=input_path, n_input_trees=n_input_trees,
                    cmd=cmd, wall_time=wall_time,
                    skipped=skipped, warnings_list=warnings_list,
                    is_error=True, error_msg=error_msg,
                    versions=versions,
                )
        else:
            error_msg = "wASTRAL completed but did not produce wastral.tre"
            return _assemble_cf_result(
                run_start=run_start,
                cf_mode=cf_mode, ref_tree=ref_tree,
                tree=tree, tree_dir=tree_dir,
                matrix=matrix, partitions=partitions,
                model=model, scf_quartets=scf_quartets,
                prefix=prefix, output_dir=output_dir,
                threads=threads,
                iqtree_path=iqtree_path, wastral_path=wastral_path,
                overwrite=overwrite, dry_run=dry_run,
                input_path=input_path, n_input_trees=n_input_trees,
                cmd=cmd, wall_time=wall_time,
                skipped=skipped, warnings_list=warnings_list,
                is_error=True, error_msg=error_msg,
                versions=versions,
            )

    # For IQ-TREE modes, the log is already written by iqtree3 as <prefix>.log

    return _assemble_cf_result(
        run_start=run_start,
        cf_mode=cf_mode, ref_tree=ref_tree,
        tree=tree, tree_dir=tree_dir,
        matrix=matrix, partitions=partitions,
        model=model, scf_quartets=scf_quartets,
        prefix=prefix, output_dir=output_dir,
        threads=threads,
        iqtree_path=iqtree_path, wastral_path=wastral_path,
        overwrite=overwrite, dry_run=dry_run,
        input_path=input_path, n_input_trees=n_input_trees,
        cmd=cmd, wall_time=wall_time,
        skipped=skipped, warnings_list=warnings_list,
        is_error=False, error_msg=None,
        versions=versions,
    )
```

- [ ] **Step 3: Run all library tests**

```bash
pytest tests/tree/test_cf.py -v
```

Expected: all tests PASS. Fix any failures before proceeding.

- [ ] **Step 4: Commit**

```bash
git add phyloai/tree/cf.py
git commit -m "feat(tree/cf): implement run_cf entry point with validation"
```

---

### Task 10: Register cf command in CLI

**Files:**
- Modify: `phyloai/cli/commands/tree.py`

- [ ] **Step 1: Add "cf" to _TreeGroup.list_commands**

```python
class _TreeGroup(click.Group):
    def list_commands(self, ctx: click.Context) -> list[str]:
        return ["ml", "msc", "cf"]
```

- [ ] **Step 2: Register the cf command (append to end of tree.py before any final lines)**

```python
@tree.command(
    "cf",
    cls=_GroupedHelpCommand,
    help=(
        "Concordance factor computation (gCF, sCF, sCFl, qCF).\n\n"
        "  --cf gcf     : gene concordance factor (IQ-TREE3)\n"
        "  --cf scf     : site concordance factor, parsimony-based (IQ-TREE3)\n"
        "  --cf scfl    : site concordance factor, likelihood-based (IQ-TREE3)\n"
        "  --cf gcf+scf : combined gCF + sCF in one IQ-TREE3 invocation\n"
        "  --cf qcf     : quartet concordance factor (wASTRAL)\n\n"
        "Input requirements by mode:\n"
        "  gcf, gcf+scf, qcf: --ref-tree + (--tree or --tree-dir)\n"
        "  scf, scfl        : --ref-tree + --matrix\n"
        "  gcf+scf           : --ref-tree + (--tree or --tree-dir) + --matrix\n"
        "  scfl              : optionally --model or --partitions for speedup\n\n"
        "CF computation is one-shot (no --resume)."
    ),
)
@click.option(
    "--cf",
    type=click.Choice(["gcf", "scf", "scfl", "gcf+scf", "qcf"]),
    required=True,
    help="Concordance factor type to compute.",
)
@click.option(
    "--ref-tree",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    required=True,
    help="Reference species tree (NEWICK).",
)
@click.option(
    "--tree",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    default=None,
    help="Single gene tree file (NEWICK, one tree per line). Mutually exclusive with --tree-dir.",
)
@click.option(
    "--tree-dir",
    type=click.Path(exists=True, file_okay=False, path_type=Path),
    default=None,
    help="Directory of gene tree files (merged into merged.trees). Mutually exclusive with --tree.",
)
@click.option(
    "--matrix",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    default=None,
    help="Multiple sequence alignment (required for scf/scfl/gcf+scf).",
)
@click.option(
    "--partitions",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    default=None,
    help="Partition file for scfl model reuse (e.g., *.best_model.nex from IQ-TREE3).",
)
@click.option(
    "--model",
    type=str,
    default=None,
    help="Substitution model for scfl speedup (e.g., LG+F+R3). Mutually exclusive with --partitions.",
)
@click.option(
    "--scf-quartets",
    type=click.IntRange(1, None),
    default=100,
    show_default=True,
    help="Number of quartets for sCF/sCFl (recommend >= 100).",
)
@click.option(
    "--prefix",
    type=str,
    default=None,
    help="Output file prefix (default: auto-derived from --cf, e.g., gCF, sCF, sCFl, gCFsCF, qCF).",
)
@click.option(
    "--output-dir", "-o",
    type=click.Path(file_okay=False, path_type=Path),
    default=Path("runs/tree/cf"),
    show_default=True,
    help="Output directory.",
)
@click.option(
    "--threads", "-t",
    type=int,
    default=4,
    show_default=True,
    help="Thread count (IQ-TREE3 -T or wASTRAL -t).",
)
@click.option(
    "--iqtree-path",
    type=Path,
    default=None,
    help="Explicit path to iqtree3 executable.",
)
@click.option(
    "--wastral-path",
    type=Path,
    default=None,
    help="Explicit path to wastral executable.",
)
@click.option("--overwrite", is_flag=True, default=False, help="Overwrite existing output directory.")
@click.option("--dry-run", is_flag=True, default=False, help="Show commands without executing.")
@click.option("--quiet", "-q", is_flag=True, default=False, help="Suppress terminal output except errors.")
def cf_command(
    cf: str,
    ref_tree: Path,
    tree: Path | None,
    tree_dir: Path | None,
    matrix: Path | None,
    partitions: Path | None,
    model: str | None,
    scf_quartets: int,
    prefix: str | None,
    output_dir: Path,
    threads: int,
    iqtree_path: Path | None,
    wastral_path: Path | None,
    overwrite: bool,
    dry_run: bool,
    quiet: bool,
) -> None:
    """Concordance factor computation (gCF, sCF, sCFl, qCF)."""
    from phyloai.tree.cf import run_cf

    # CLI-layer early validation for executable paths
    if iqtree_path is not None:
        if not iqtree_path.exists():
            _fail(f"--iqtree-path does not exist: {iqtree_path}", 1)
        if not os.access(str(iqtree_path), os.X_OK):
            _fail(f"--iqtree-path is not executable: {iqtree_path}", 1)
    if wastral_path is not None:
        if not wastral_path.exists():
            _fail(f"--wastral-path does not exist: {wastral_path}", 1)
        if not os.access(str(wastral_path), os.X_OK):
            _fail(f"--wastral-path is not executable: {wastral_path}", 1)

    error_msg: str | None = None

    try:
        payload = run_cf(
            cf_mode=cf,
            ref_tree=ref_tree,
            tree=tree,
            tree_dir=tree_dir,
            matrix=matrix,
            partitions=partitions,
            model=model,
            scf_quartets=scf_quartets,
            prefix=prefix,
            output_dir=output_dir,
            threads=threads,
            iqtree_path=str(iqtree_path) if iqtree_path else None,
            wastral_path=str(wastral_path) if wastral_path else None,
            overwrite=overwrite,
            dry_run=dry_run,
            quiet=quiet,
        )
    except (ValueError, FileNotFoundError) as exc:
        error_msg = str(exc)
    except SystemExit:
        raise
    except Exception as exc:
        error_msg = str(exc)

    if error_msg is not None:
        if "iqtree3 not found" in error_msg.lower() or "wastral not found" in error_msg.lower():
            exit_code = 3
        else:
            exit_code = 1
        _fail(error_msg, exit_code)

    if dry_run:
        if not quiet:
            cf_type = payload["key_results"]["cf_type"]
            click.echo(f"Dry run: --cf {cf_type} would be executed.")
            click.echo(" ".join(payload["data"]["cmd"]))
        return

    # Write result.json
    result_path = output_dir / "result.json"
    result_path.parent.mkdir(parents=True, exist_ok=True)
    with open(result_path, "w") as fh:
        json.dump(payload, fh, indent=2)

    if payload["status"] == "error":
        _fail(payload.get("error", "CF computation failed"), 2)

    if not quiet:
        prefix_val = payload["key_results"]["prefix"]
        if cf in ("gcf", "scf", "scfl", "gcf+scf"):
            click.echo(f"CF results saved to {output_dir}/{prefix_val}.cf.*")
        else:
            click.echo(f"qCF tree saved to {output_dir}/{prefix_val}.cf.tree")
```

- [ ] **Step 3: Verify help text renders**

```bash
python -m phyloai tree cf --help
```

Expected: help text shows all options, organized logically.

- [ ] **Step 4: Commit**

```bash
git add phyloai/cli/commands/tree.py
git commit -m "feat(cli): register tree cf command with --cf mode selector"
```

---

### Task 11: Write CLI integration tests

**Files:**
- Create: `tests/cli/test_tree_cf.py`

- [ ] **Step 1: Write CLI test file**

```python
"""CLI integration tests for phyloai tree cf."""

from __future__ import annotations

from pathlib import Path

from click.testing import CliRunner

from phyloai.cli.main import cli


def test_tree_cf_help_shows_all_flags() -> None:
    """--help lists all expected parameters."""
    result = CliRunner().invoke(cli, ["tree", "cf", "--help"])
    assert result.exit_code == 0
    for flag in [
        "--cf", "--ref-tree", "--tree", "--tree-dir",
        "--matrix", "--partitions", "--model",
        "--scf-quartets", "--prefix",
        "--output-dir", "--threads", "--iqtree-path",
        "--wastral-path", "--overwrite", "--dry-run", "--quiet",
    ]:
        assert flag in result.output, f"Missing flag: {flag}"


def test_tree_cf_gcf_dry_run(tmp_path: Path) -> None:
    """--cf gcf --dry-run with valid inputs produces exit 0."""
    ref_tree = tmp_path / "ref.nwk"
    ref_tree.write_text("(A,B);")
    gene_trees = tmp_path / "trees.nwk"
    gene_trees.write_text("(A,B);\n")
    out_dir = tmp_path / "out"

    result = CliRunner().invoke(cli, [
        "tree", "cf", "--cf", "gcf",
        "--ref-tree", str(ref_tree),
        "--tree", str(gene_trees),
        "-o", str(out_dir), "--dry-run",
    ])
    assert result.exit_code == 0
    assert "Dry run" in result.output


def test_tree_cf_scf_dry_run(tmp_path: Path) -> None:
    """--cf scf --dry-run with ref-tree + matrix produces exit 0."""
    ref_tree = tmp_path / "ref.nwk"
    ref_tree.write_text("(A,B);")
    matrix = tmp_path / "msa.fa"
    matrix.write_text(">A\nACGT\n>B\nACGT\n")
    out_dir = tmp_path / "out"

    result = CliRunner().invoke(cli, [
        "tree", "cf", "--cf", "scf",
        "--ref-tree", str(ref_tree),
        "--matrix", str(matrix),
        "-o", str(out_dir), "--dry-run",
    ])
    assert result.exit_code == 0


def test_tree_cf_qcf_dry_run(tmp_path: Path) -> None:
    """--cf qcf --dry-run with valid inputs produces exit 0."""
    ref_tree = tmp_path / "ref.nwk"
    ref_tree.write_text("(A,B);")
    gene_trees = tmp_path / "trees.nwk"
    gene_trees.write_text("(A,B);\n")
    out_dir = tmp_path / "out"

    result = CliRunner().invoke(cli, [
        "tree", "cf", "--cf", "qcf",
        "--ref-tree", str(ref_tree),
        "--tree", str(gene_trees),
        "-o", str(out_dir), "--dry-run",
    ])
    assert result.exit_code == 0


def test_tree_cf_missing_cf_flag_shows_error() -> None:
    """Missing --cf flag produces exit 2 (Click's missing-option code)."""
    result = CliRunner().invoke(cli, [
        "tree", "cf", "--ref-tree", "/fake/path",
    ])
    assert result.exit_code != 0


def test_tree_cf_scf_without_matrix_exits_1(tmp_path: Path) -> None:
    """--cf scf without --matrix exits 1."""
    ref_tree = tmp_path / "ref.nwk"
    ref_tree.write_text("(A,B);")
    out_dir = tmp_path / "out"

    result = CliRunner().invoke(cli, [
        "tree", "cf", "--cf", "scf",
        "--ref-tree", str(ref_tree),
        "-o", str(out_dir), "--dry-run",
    ])
    assert result.exit_code == 1


def test_tree_cf_gcf_with_matrix_exits_1(tmp_path: Path) -> None:
    """--cf gcf with --matrix exits 1."""
    ref_tree = tmp_path / "ref.nwk"
    ref_tree.write_text("(A,B);")
    gene_trees = tmp_path / "trees.nwk"
    gene_trees.write_text("(A,B);\n")
    matrix = tmp_path / "msa.fa"
    matrix.write_text(">A\nA\n")
    out_dir = tmp_path / "out"

    result = CliRunner().invoke(cli, [
        "tree", "cf", "--cf", "gcf",
        "--ref-tree", str(ref_tree),
        "--tree", str(gene_trees),
        "--matrix", str(matrix),
        "-o", str(out_dir), "--dry-run",
    ])
    assert result.exit_code == 1


def test_tree_cf_scfl_model_and_partitions_exits_1(tmp_path: Path) -> None:
    """--cf scfl with --model and --partitions exits 1."""
    ref_tree = tmp_path / "ref.nwk"
    ref_tree.write_text("(A,B);")
    matrix = tmp_path / "msa.fa"
    matrix.write_text(">A\nA\n")
    partitions = tmp_path / "p.nex"
    partitions.write_text("#nexus")
    out_dir = tmp_path / "out"

    result = CliRunner().invoke(cli, [
        "tree", "cf", "--cf", "scfl",
        "--ref-tree", str(ref_tree),
        "--matrix", str(matrix),
        "--model", "LG", "--partitions", str(partitions),
        "-o", str(out_dir), "--dry-run",
    ])
    assert result.exit_code == 1


def test_tree_cf_explicit_prefix(tmp_path: Path) -> None:
    """--prefix myCF overrides default."""
    ref_tree = tmp_path / "ref.nwk"
    ref_tree.write_text("(A,B);")
    gene_trees = tmp_path / "trees.nwk"
    gene_trees.write_text("(A,B);\n")
    out_dir = tmp_path / "out"

    result = CliRunner().invoke(cli, [
        "tree", "cf", "--cf", "gcf",
        "--ref-tree", str(ref_tree),
        "--tree", str(gene_trees),
        "--prefix", "myCF",
        "-o", str(out_dir), "--dry-run",
    ])
    assert result.exit_code == 0


def test_tree_cf_tree_and_tree_dir_mutually_exclusive(tmp_path: Path) -> None:
    """--tree and --tree-dir together exits 1."""
    ref_tree = tmp_path / "ref.nwk"
    ref_tree.write_text("(A,B);")
    gene_trees = tmp_path / "trees.nwk"
    gene_trees.write_text("(A,B);\n")
    td = tmp_path / "tdir"
    td.mkdir()
    out_dir = tmp_path / "out"

    result = CliRunner().invoke(cli, [
        "tree", "cf", "--cf", "gcf",
        "--ref-tree", str(ref_tree),
        "--tree", str(gene_trees),
        "--tree-dir", str(td),
        "-o", str(out_dir), "--dry-run",
    ])
    assert result.exit_code == 1
```

- [ ] **Step 2: Run CLI tests**

```bash
pytest tests/cli/test_tree_cf.py -v
```

Expected: all PASS.

- [ ] **Step 3: Commit**

```bash
git add tests/cli/test_tree_cf.py
git commit -m "test(cli): add integration tests for tree cf command"
```

---

### Task 12: Create user-facing documentation

**Files:**
- Create: `docs/commands/tree-cf.md`

- [ ] **Step 1: Write docs/commands/tree-cf.md**

```markdown
# phyloai tree cf

Concordance factor computation — quantify branch support using gene concordance (gCF), site concordance (sCF, sCFl), and quartet concordance (qCF).

## Purpose

Concordance factors measure the proportion of gene trees, sites, or quartets supporting each bipartition in a reference species tree. They provide complementary branch-support information beyond standard bootstrap.

## Usage

```bash
phyloai tree cf --cf MODE --ref-tree REF_TREE [INPUTS...] [OPTIONS]
```

## Modes

| `--cf`    | Index  | Tool    | Description |
|-----------|--------|---------|-------------|
| `gcf`       | gCF    | IQ-TREE3 | Gene concordance factor |
| `scf`       | sCF    | IQ-TREE3 | Site concordance factor (parsimony) |
| `scfl`      | sCFl   | IQ-TREE3 | Site concordance factor (likelihood) |
| `gcf+scf`   | gCF+sCF | IQ-TREE3 | Combined gCF + sCF in one run |
| `qcf`       | qCF    | wASTRAL  | Quartet concordance factor |

## Input Requirements by Mode

| Mode     | `--ref-tree` | `--tree`/`--tree-dir` | `--matrix` | `--model`/`--partitions` |
|----------|-------------|----------------------|-----------|-------------------------|
| `gcf`      | Required    | Required             | —         | —                       |
| `scf`      | Required    | —                    | Required  | —                       |
| `scfl`     | Required    | —                    | Required  | Optional (speedup)      |
| `gcf+scf`  | Required    | Required             | Required  | —                       |
| `qcf`      | Required    | Required             | —         | —                       |

## Examples

```bash
# gCF: gene trees + reference tree
phyloai tree cf --cf gcf --ref-tree species.nwk --tree-dir ./genetrees/

# gCF with single file
phyloai tree cf --cf gcf --ref-tree species.nwk --tree merged.trees

# sCF: alignment + reference tree (ideally gCF-annotated)
phyloai tree cf --cf scf --ref-tree gCF.cf.tree --matrix msa.fa

# sCFl (likelihood) with model for speedup
phyloai tree cf --cf scfl --ref-tree gCF.cf.tree --matrix msa.fa --model LG+F+R3

# sCFl with pre-computed partition model
phyloai tree cf --cf scfl --ref-tree gCF.cf.tree --matrix msa.fa \
    --partitions msa.best_model.nex

# Combined gCF + sCF
phyloai tree cf --cf gcf+scf --ref-tree species.nwk --tree-dir ./genetrees/ \
    --matrix msa.fa

# qCF via wASTRAL
phyloai tree cf --cf qcf --ref-tree species.nwk --tree merged.trees

# Custom output prefix and threads
phyloai tree cf --cf gcf --ref-tree species.nwk --tree merged.trees \
    --prefix myCF -t 8
```

## Output Files

### IQ-TREE3 modes (gcf, scf, scfl, gcf+scf)

| File | Description |
|------|-------------|
| `<prefix>.cf.stat`  | Concordance factor statistics table |
| `<prefix>.cf.branch` | Tree with branch IDs |
| `<prefix>.cf.tree`  | Tree annotated with CF values |
| `<prefix>.cf.tree.nex` | NEXUS annotated tree for FigTree |
| `<prefix>.log`       | IQ-TREE3 log |
| `result.json`        | PhyloAI structured result |
| `merged.trees`       | Merged gene trees (if `--tree-dir` used) |

### qCF mode

| File | Description |
|------|-------------|
| `<prefix>.cf.tree` | Reference tree with qCF values appended |
| `wastral.tre`      | Raw wASTRAL output (intermediate) |
| `wastral.log`      | wASTRAL log |
| `result.json`      | PhyloAI structured result |
| `merged.trees`     | Merged gene trees (if `--tree-dir` used) |

## Notes

- For best sCF/sCFl results, use a gCF-annotated tree as `--ref-tree` (e.g., run `--cf gcf` first).
- `--cf scfl` without `--model` or `--partitions` auto-computes the best-fit model — this is slow. Provide `--model` or `--partitions` for speedup.
- qCF values are multiplied by 100 and rounded to integers before appending to reference tree support values.
- `--scf-quartets` should be >= 100 for reliable results.
```

- [ ] **Step 2: Commit**

```bash
git add docs/commands/tree-cf.md
git commit -m "docs: add user documentation for tree cf command"
```

---

### Task 13: Final verification

- [ ] **Step 1: Run full test suite**

```bash
pytest tests/tree/test_cf.py tests/cli/test_tree_cf.py -v
```

Expected: all tests PASS.

- [ ] **Step 2: Verify help output**

```bash
python -m phyloai tree cf --help
python -m phyloai tree --help
```

Expected: both show `cf` in the command list.

- [ ] **Step 3: Run lint/typecheck**

```bash
ruff check phyloai/tree/cf.py phyloai/cli/commands/tree.py
```

- [ ] **Step 4: Commit any remaining lint/test fixes**

```bash
git add -A && git diff --cached --stat
```

Review staged changes, then:

```bash
git commit -m "chore: final lint and test fixes for tree cf"
```

---

## Self-Review

**Spec coverage:**
- §1 Purpose: covered by Task 9 (run_cf entry point, mode dispatch) ✓
- §2 CLI Surface: covered by Task 10 (CLI registration) ✓
- §3 Parameter Specification: covered by Task 10 (click options) + Task 9 (validation) ✓
- §3.4 Default Prefix: covered by Task 3 (test) + Task 5 (constants) + Task 9 (resolution) ✓
- §4.1-4.4 IQ-TREE commands: covered by Task 2 (tests) + Task 7 (_build_iqtree_cf_cmd) ✓
- §4.5 wASTRAL command: covered by Task 2 (test) + Task 7 (_build_wastral_qcf_cmd) ✓
- §4.6 qCF mapping: covered by Task 1 (tests) + Task 6 (_map_qcf_to_tree) ✓
- §5 Input Validation: covered by Task 3 (tests) + Task 9 (validation in run_cf) ✓
- §6 Output Dir Structure: handled by subprocess cwd + tool-native output ✓
- §7 result.json Schema: covered by Task 4 (test) + Task 9 (_assemble_cf_result) ✓
- §8 Exit Codes: covered by Task 10 (CLI error handling) ✓
- §9 Warnings: covered by Task 3/4 (tests) + Task 9 (warnings in run_cf) ✓
- §10 Logging: covered by Task 9 (wastral.log save, iqtree native log) ✓
- §12 Acceptance Criteria: all major criteria covered by test tasks ✓

**Placeholder scan:** No TBD, TODO, or incomplete sections. All code is shown inline.

**Type consistency:** `cf_mode: str` used consistently. `prefix: str | None` with default resolution. All function signatures match across tasks.

**Addendum — Implementation deviations from plan:**

After implementation, the following design refinements were applied based on code review and user feedback:

1. **qCF values kept as decimal [0,1]** (not integer ×100): The original plan specified `round(q1 * 100)` producing integer labels (e.g., `42`). Per user request, qCF values are now raw decimals formatted to 4 decimal places (e.g., `0.4221`). This preserves the original wASTRAL precision. Tests updated accordingly.

2. **`--lpp` flag for local posterior probabilities**: New CLI flag between `--scf-quartets` and `--prefix`. When set, pp1 is extracted from wASTRAL labels and appended after q1: `<support>/<q1>/<pp1>`. Without `--lpp`, only q1 is appended.

3. **`cf.log` (module log)**: Added project-level log file (matching `msc.log` pattern from `tree msc`). Written to output directory containing timestamp, command, tool versions, wall time, input tree count, and output file paths.

4. **Real-time subprocess streaming**: Changed from `subprocess.run` with `PIPE` to `subprocess.Popen` with `select`-based streaming. Tools' stdout/stderr is written to terminal in real-time while also captured for log files. When `--quiet` is set, output is suppressed.

5. **Input path resolution**: All input paths (`--ref-tree`, `--matrix`, `--tree`, `--partitions`, `merged.trees`) are resolved to absolute before being passed to subprocess commands. Executable paths (`--iqtree-path`, `--wastral-path`) also resolved via `Path.resolve()`. This prevents failures when subprocesses run with `cwd=output_dir`.

6. **Help text formatting**: Created `_CFCommand` class (subclass of `_GroupedHelpCommand`) that overrides `format_help_text` to preserve newline formatting, ensuring modes and input requirements render one-per-line.

7. **`--partitions` path resolved**: Added validation and `.resolve()` for `--partitions` Path, fixing the relative-path not-found error when IQ-TREE runs from cwd=output_dir.

8. **Raw Newick qCF mapping (no Bio.Phylo round-trip)**: `_map_qcf_to_tree` no longer uses `Phylo.write` to output the annotated tree. Instead, it walks the raw ref-tree Newick string character-by-character and injects qCF/pp1 annotations directly into the original labels. This preserves branch-length precision (e.g., 10 decimal places) and existing support-label formatting exactly — only the new `/q1[/pp1]` values are appended. The `_fmt_val` helper strips unnecessary trailing zeros (e.g., `1` not `1.0000`, `0.95` not `0.9500`).
