"""Tests for phyloai/tree/cf.py — concordance factor computation."""

from __future__ import annotations

from pathlib import Path

import pytest


# ── qCF mapper tests ────────────────────────────────────────────────


def test_map_qcf_to_tree_appends_qcf_to_existing_support(tmp_path: Path) -> None:
    """q1=0.422083 -> 42, appended to '100/90' as '100/90/42'."""
    from phyloai.tree.cf import _map_qcf_to_tree
    from Bio import Phylo

    ref_nwk = tmp_path / "ref.nwk"
    ref_nwk.write_text("((A:0.1,B:0.1)100/90:0.05,C:0.15);")

    wastral_nwk = tmp_path / "wastral.tre"
    wastral_nwk.write_text(
        "((A:0.1,B:0.1)'[q1=0.422083;q2=0.288958;q3=0.288958]':0.05,C:0.15);"
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
    assert ab_node is not None, "Could not find (A,B) clade"
    assert ab_node.name == "100/90/0.4221"


def test_map_qcf_to_tree_no_existing_support(tmp_path: Path) -> None:
    """When ref tree has no support, qCF becomes the sole label."""
    from phyloai.tree.cf import _map_qcf_to_tree
    from Bio import Phylo

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
    assert ab_node.confidence == 0.75


def test_map_qcf_to_tree_handles_complementary_rooting(tmp_path: Path) -> None:
    """4-taxon: ref annotates (A,B) side but wASTRAL labels (C,D) side.
    Canonical bipartition matching must handle the complement."""
    from phyloai.tree.cf import _map_qcf_to_tree
    from Bio import Phylo

    ref_nwk = tmp_path / "ref.nwk"
    ref_nwk.write_text("((A:0.1,B:0.1)100:0.05,(C:0.1,D:0.1):0.05);")

    wastral_nwk = tmp_path / "wastral.tre"
    wastral_nwk.write_text(
        "((A:0.1,B:0.1):0.05,(C:0.1,D:0.1)'[q1=0.6;q2=0.2;q3=0.2]':0.05);"
    )

    output_nwk = tmp_path / "qCF.cf.tree"
    _map_qcf_to_tree(ref_nwk, wastral_nwk, output_nwk)

    tree = Phylo.read(str(output_nwk), "newick")
    ab_clade = None
    for c in tree.find_clades():
        if not c.is_terminal():
            leaves = frozenset(t.name for t in c.get_terminals())
            if leaves == frozenset({"A", "B"}):
                ab_clade = c
                break
    assert ab_clade is not None
    assert ab_clade.name == "100/0.6"


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
    assert result.count("/") >= 2


def test_map_qcf_to_tree_missing_q1_skips_node(tmp_path: Path) -> None:
    """Node without q1 annotation should be left unchanged."""
    from phyloai.tree.cf import _map_qcf_to_tree

    ref_nwk = tmp_path / "ref.nwk"
    ref_nwk.write_text("((A:0.1,B:0.1)100:0.05,C:0.15);")

    wastral_nwk = tmp_path / "wastral.tre"
    wastral_nwk.write_text(
        "((A:0.1,B:0.1)'[pp1=0.9]':0.05,C:0.15);"
    )

    output_nwk = tmp_path / "qCF.cf.tree"
    _map_qcf_to_tree(ref_nwk, wastral_nwk, output_nwk)

    result = output_nwk.read_text().strip()
    assert "100" in result


# ── Command builder tests ────────────────────────────────────────────


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


def test_build_iqtree_cf_cmd_scfl_with_model(tmp_path: Path) -> None:
    from phyloai.tree.cf import _build_iqtree_cf_cmd

    ref_tree = tmp_path / "gCF.cf.tree"
    ref_tree.write_text("(A,B);")
    matrix = tmp_path / "msa.fa"
    matrix.write_text(">A\nACGT\n>B\nACGT\n")

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


def test_build_iqtree_cf_cmd_scfl_with_partitions(tmp_path: Path) -> None:
    from phyloai.tree.cf import _build_iqtree_cf_cmd

    ref_tree = tmp_path / "gCF.cf.tree"
    ref_tree.write_text("(A,B);")
    matrix = tmp_path / "msa.fa"
    matrix.write_text(">A\nACGT\n>B\nACGT\n")
    partitions = tmp_path / "msa.best_model.nex"
    partitions.write_text("#nexus")

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


def test_build_iqtree_cf_cmd_scfl_auto_model(tmp_path: Path) -> None:
    from phyloai.tree.cf import _build_iqtree_cf_cmd

    ref_tree = tmp_path / "gCF.cf.tree"
    ref_tree.write_text("(A,B);")
    matrix = tmp_path / "msa.fa"
    matrix.write_text(">A\nACGT\n>B\nACGT\n")

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
    assert cmd == [
        "iqtree3", "-s", str(matrix), "-te", str(ref_tree),
        "--scfl", "100", "--prefix", "sCFl", "-T", "4",
    ]


def test_build_iqtree_cf_cmd_gcf_scf_combined(tmp_path: Path) -> None:
    from phyloai.tree.cf import _build_iqtree_cf_cmd

    ref_tree = tmp_path / "species.nwk"
    ref_tree.write_text("(A,B);")
    gene_trees = tmp_path / "merged.trees"
    gene_trees.write_text("(A,B);")
    matrix = tmp_path / "msa.fa"
    matrix.write_text(">A\nACGT\n>B\nACGT\n")

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
        "iqtree3", "-s", str(matrix), "-t", str(ref_tree), "--gcf", str(gene_trees),
        "--scf", "150", "--prefix", "gCFsCF", "-T", "4",
    ]


def test_build_wastral_qcf_cmd(tmp_path: Path) -> None:
    from phyloai.tree.cf import _build_wastral_qcf_cmd

    ref_tree = tmp_path / "species.nwk"
    ref_tree.write_text("(A,B);")
    gene_trees = tmp_path / "merged.trees"
    gene_trees.write_text("(A,B);")
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


# ── Validation tests ──────────────────────────────────────────────────


def test_run_cf_scf_without_matrix_raises(tmp_path: Path) -> None:
    from phyloai.tree.cf import run_cf

    ref_tree = tmp_path / "ref.nwk"
    ref_tree.write_text("(A,B);")

    with pytest.raises(ValueError, match="--matrix is required"):
        run_cf(
            cf_mode="scf", ref_tree=ref_tree,
            output_dir=tmp_path / "out",
            threads=4, dry_run=True,
        )


def test_run_cf_gcf_with_matrix_raises(tmp_path: Path) -> None:
    from phyloai.tree.cf import run_cf

    ref_tree = tmp_path / "ref.nwk"
    ref_tree.write_text("(A,B);")
    gene_trees = tmp_path / "trees"
    gene_trees.write_text("(A,B);")
    matrix = tmp_path / "msa.fa"
    matrix.write_text(">A\nA\n")

    with pytest.raises(ValueError, match="not valid"):
        run_cf(
            cf_mode="gcf", ref_tree=ref_tree,
            tree=gene_trees, matrix=matrix,
            output_dir=tmp_path / "out",
            threads=4, dry_run=True,
        )


def test_run_cf_scfl_model_and_partitions_mutually_exclusive(tmp_path: Path) -> None:
    from phyloai.tree.cf import run_cf

    ref_tree = tmp_path / "ref.nwk"
    ref_tree.write_text("(A,B);")
    matrix = tmp_path / "msa.fa"
    matrix.write_text(">A\nA\n")
    partitions = tmp_path / "p.nex"
    partitions.write_text("#nexus")

    with pytest.raises(ValueError, match="mutually exclusive"):
        run_cf(
            cf_mode="scfl", ref_tree=ref_tree,
            matrix=matrix, model="LG", partitions=partitions,
            output_dir=tmp_path / "out",
            threads=4, dry_run=True,
        )


def test_run_cf_qcf_with_matrix_raises(tmp_path: Path) -> None:
    from phyloai.tree.cf import run_cf

    ref_tree = tmp_path / "ref.nwk"
    ref_tree.write_text("(A,B);")
    gene_trees = tmp_path / "trees"
    gene_trees.write_text("(A,B);")
    matrix = tmp_path / "msa.fa"
    matrix.write_text(">A\nA\n")

    with pytest.raises(ValueError, match="not valid"):
        run_cf(
            cf_mode="qcf", ref_tree=ref_tree,
            tree=gene_trees, matrix=matrix,
            output_dir=tmp_path / "out",
            threads=4, dry_run=True,
        )


def test_run_cf_scf_with_tree_raises(tmp_path: Path) -> None:
    from phyloai.tree.cf import run_cf

    ref_tree = tmp_path / "ref.nwk"
    ref_tree.write_text("(A,B);")
    gene_trees = tmp_path / "trees"
    gene_trees.write_text("(A,B);")
    matrix = tmp_path / "msa.fa"
    matrix.write_text(">A\nA\n")

    with pytest.raises(ValueError, match="not needed"):
        run_cf(
            cf_mode="scf", ref_tree=ref_tree,
            tree=gene_trees, matrix=matrix,
            output_dir=tmp_path / "out",
            threads=4, dry_run=True,
        )


def test_run_cf_gcf_with_model_raises(tmp_path: Path) -> None:
    from phyloai.tree.cf import run_cf

    ref_tree = tmp_path / "ref.nwk"
    ref_tree.write_text("(A,B);")
    gene_trees = tmp_path / "trees"
    gene_trees.write_text("(A,B);")

    with pytest.raises(ValueError, match="not valid"):
        run_cf(
            cf_mode="gcf", ref_tree=ref_tree,
            tree=gene_trees, model="LG",
            output_dir=tmp_path / "out",
            threads=4, dry_run=True,
        )


def test_run_cf_qcf_with_scf_quartets_raises(tmp_path: Path) -> None:
    from phyloai.tree.cf import run_cf

    ref_tree = tmp_path / "ref.nwk"
    ref_tree.write_text("(A,B);")
    gene_trees = tmp_path / "trees"
    gene_trees.write_text("(A,B);")

    with pytest.raises(ValueError, match="not valid"):
        run_cf(
            cf_mode="qcf", ref_tree=ref_tree,
            tree=gene_trees, scf_quartets=200,
            output_dir=tmp_path / "out",
            threads=4, dry_run=True,
        )


def test_run_cf_tree_and_tree_dir_mutually_exclusive(tmp_path: Path) -> None:
    from phyloai.tree.cf import run_cf

    ref_tree = tmp_path / "ref.nwk"
    ref_tree.write_text("(A,B);")
    gene_trees = tmp_path / "trees"
    gene_trees.write_text("(A,B);")
    tree_dir = tmp_path / "tdir"
    tree_dir.mkdir()

    with pytest.raises(ValueError, match="mutually exclusive"):
        run_cf(
            cf_mode="gcf", ref_tree=ref_tree,
            tree=gene_trees, tree_dir=tree_dir,
            output_dir=tmp_path / "out",
            threads=4, dry_run=True,
        )


def test_run_cf_gcf_without_gene_trees_raises(tmp_path: Path) -> None:
    from phyloai.tree.cf import run_cf

    ref_tree = tmp_path / "ref.nwk"
    ref_tree.write_text("(A,B);")

    with pytest.raises(ValueError, match="--tree or --tree-dir"):
        run_cf(
            cf_mode="gcf", ref_tree=ref_tree,
            output_dir=tmp_path / "out",
            threads=4, dry_run=True,
        )


def test_run_cf_scf_quartets_below_100_warns(tmp_path: Path) -> None:
    from phyloai.tree.cf import run_cf

    ref_tree = tmp_path / "ref.nwk"
    ref_tree.write_text("(A,B);")
    matrix = tmp_path / "msa.fa"
    matrix.write_text(">A\nACGT\n>B\nACGT\n")

    result = run_cf(
        cf_mode="scf", ref_tree=ref_tree,
        matrix=matrix, scf_quartets=50,
        output_dir=tmp_path / "out",
        threads=4, dry_run=True,
    )
    assert any(">= 100" in w for w in result.get("data", {}).get("warnings", []))


def test_run_cf_default_prefix_per_mode(tmp_path: Path) -> None:
    from phyloai.tree.cf import run_cf

    ref_tree = tmp_path / "ref.nwk"
    ref_tree.write_text("(A,B);")
    gene_trees = tmp_path / "trees"
    gene_trees.write_text("(A,B);\n")
    matrix = tmp_path / "msa.fa"
    matrix.write_text(">A\nACGT\n>B\nACGT\n")

    for mode, expected in [("gcf", "gCF"), ("qcf", "qCF")]:
        result = run_cf(
            cf_mode=mode, ref_tree=ref_tree,
            tree=gene_trees,
            output_dir=tmp_path / "out",
            threads=4, dry_run=True,
        )
        assert result["params"]["prefix"] == expected

    for mode, expected in [("scf", "sCF"), ("scfl", "sCFl")]:
        result = run_cf(
            cf_mode=mode, ref_tree=ref_tree,
            matrix=matrix,
            output_dir=tmp_path / "out",
            threads=4, dry_run=True,
        )
        assert result["params"]["prefix"] == expected

    result = run_cf(
        cf_mode="gcf+scf", ref_tree=ref_tree,
        tree=gene_trees, matrix=matrix,
        output_dir=tmp_path / "out",
        threads=4, dry_run=True,
    )
    assert result["params"]["prefix"] == "gCFsCF"


# ── Dry-run flow tests ────────────────────────────────────────────────


def test_run_cf_dry_run_gcf_produces_payload(tmp_path: Path) -> None:
    from phyloai.tree.cf import run_cf
    from tests.helpers import validate_params_completeness, validate_result_json

    ref_tree = tmp_path / "ref.nwk"
    ref_tree.write_text("(A,B);")
    gene_trees = tmp_path / "trees"
    gene_trees.write_text("(A,B);")

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

    validate_result_json(result)
    validate_params_completeness(result, {
        "cf", "ref_tree", "tree", "tree_dir", "matrix", "partitions",
        "model", "scf_quartets", "lpp", "prefix", "output_dir", "threads",
        "overwrite", "dry_run", "iqtree_path", "wastral_path", "quiet",
    })


def test_run_cf_dry_run_qcf_produces_payload(tmp_path: Path) -> None:
    from phyloai.tree.cf import run_cf

    ref_tree = tmp_path / "ref.nwk"
    ref_tree.write_text("(A,B);")
    gene_trees = tmp_path / "trees"
    gene_trees.write_text("(A,B);")

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


def test_run_cf_dry_run_scfl_with_model(tmp_path: Path) -> None:
    from phyloai.tree.cf import run_cf

    ref_tree = tmp_path / "ref.nwk"
    ref_tree.write_text("(A,B);")
    matrix = tmp_path / "msa.fa"
    matrix.write_text(">A\nACGT\n>B\nACGT\n")

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


def test_run_cf_dry_run_scfl_with_partitions(tmp_path: Path) -> None:
    from phyloai.tree.cf import run_cf

    ref_tree = tmp_path / "ref.nwk"
    ref_tree.write_text("(A,B);")
    matrix = tmp_path / "msa.fa"
    matrix.write_text(">A\nACGT\n>B\nACGT\n")
    partitions = tmp_path / "best_model.nex"
    partitions.write_text("#nexus")

    result = run_cf(
        cf_mode="scfl", ref_tree=ref_tree,
        matrix=matrix, partitions=partitions,
        output_dir=tmp_path / "out",
        threads=4, dry_run=True,
    )

    cmd = result["data"]["cmd"]
    assert "-p" in cmd
    assert str(partitions) in cmd


def test_run_cf_explicit_prefix_overrides_default(tmp_path: Path) -> None:
    from phyloai.tree.cf import run_cf

    ref_tree = tmp_path / "ref.nwk"
    ref_tree.write_text("(A,B);")
    gene_trees = tmp_path / "trees"
    gene_trees.write_text("(A,B);")

    result = run_cf(
        cf_mode="gcf", ref_tree=ref_tree,
        tree=gene_trees, prefix="myCF",
        output_dir=tmp_path / "out",
        threads=4, dry_run=True,
    )

    assert result["params"]["prefix"] == "myCF"
    assert result["key_results"]["prefix"] == "myCF"


def test_run_cf_dry_run_tree_dir_merges_trees(tmp_path: Path) -> None:
    from phyloai.tree.cf import run_cf

    ref_tree = tmp_path / "ref.nwk"
    ref_tree.write_text("(A,B);")
    td = tmp_path / "tdir"
    td.mkdir()
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


def test_run_cf_tree_dir_zero_valid_files_raises(tmp_path: Path) -> None:
    from phyloai.tree.cf import run_cf

    ref_tree = tmp_path / "ref.nwk"
    ref_tree.write_text("(A,B);")
    td = tmp_path / "tdir"
    td.mkdir()
    (td / "data.txt").write_text("not a tree")

    with pytest.raises(ValueError, match="No valid gene tree files"):
        run_cf(
            cf_mode="gcf", ref_tree=ref_tree,
            tree_dir=td,
            output_dir=tmp_path / "out",
            threads=4, dry_run=True,
        )


def test_run_cf_tree_dir_one_file_warns(tmp_path: Path) -> None:
    from phyloai.tree.cf import run_cf

    ref_tree = tmp_path / "ref.nwk"
    ref_tree.write_text("(A,B);")
    td = tmp_path / "tdir"
    td.mkdir()
    (td / "gene1.nwk").write_text("(A,B);\n")

    result = run_cf(
        cf_mode="gcf", ref_tree=ref_tree,
        tree_dir=td,
        output_dir=tmp_path / "out",
        threads=4, dry_run=True,
    )

    warnings = result.get("data", {}).get("warnings", [])
    assert any("--tree" in w for w in warnings)


def test_run_cf_lpp_only_valid_for_qcf(tmp_path: Path) -> None:
    from phyloai.tree.cf import run_cf

    ref_tree = tmp_path / "ref.nwk"
    ref_tree.write_text("(A,B);")
    gene_trees = tmp_path / "genes.nwk"
    gene_trees.write_text("(A,B);\n")

    # gcf should reject --lpp
    with pytest.raises(ValueError, match="--lpp is only valid"):
        run_cf(
            cf_mode="gcf", ref_tree=ref_tree,
            tree=gene_trees, lpp=True,
            output_dir=tmp_path / "out",
            threads=4, dry_run=True,
        )

    # scf should reject --lpp
    matrix = tmp_path / "msa.fa"
    matrix.write_text(">A\nACGT\n>B\nACGT\n")
    with pytest.raises(ValueError, match="--lpp is only valid"):
        run_cf(
            cf_mode="scf", ref_tree=ref_tree,
            matrix=matrix, lpp=True,
            output_dir=tmp_path / "out",
            threads=4, dry_run=True,
        )

    # qcf should accept --lpp
    result = run_cf(
        cf_mode="qcf", ref_tree=ref_tree,
        tree=gene_trees, lpp=True,
        output_dir=tmp_path / "out",
        threads=4, dry_run=True,
    )
    assert result["params"]["lpp"] is True
    assert "--lpp" in result["command"]


def test_map_qcf_raw_newick_preserves_branch_lengths(tmp_path: Path) -> None:
    """Verify _map_qcf_to_tree preserves original branch lengths exactly."""
    from phyloai.tree.cf import _map_qcf_to_tree

    ref_nwk = tmp_path / "ref.nwk"
    ref_nwk.write_text("((A:0.1234567890,B:0.0987654321)100:0.0012345678,C:0.5555555555);")

    wastral_nwk = tmp_path / "wastral.tre"
    wastral_nwk.write_text("((A:0.1,B:0.1)'[q1=0.7500]',C:0.1);")

    output_nwk = tmp_path / "qCF.cf.tree"
    _map_qcf_to_tree(ref_nwk, wastral_nwk, output_nwk)

    out = output_nwk.read_text()
    assert "0.1234567890" in out
    assert "0.0987654321" in out
    assert "0.0012345678" in out
    assert "0.5555555555" in out
    assert "100/0.75" in out


def test_map_qcf_raw_newick_with_lpp(tmp_path: Path) -> None:
    """Verify _map_qcf_to_tree with --lpp appends pp1 in raw output."""
    from phyloai.tree.cf import _map_qcf_to_tree

    ref_nwk = tmp_path / "ref.nwk"
    ref_nwk.write_text("((A:0.1,B:0.1)100:0.05,C:0.15);")

    wastral_nwk = tmp_path / "wastral.tre"
    wastral_nwk.write_text("((A:0.1,B:0.1)'[q1=0.5;pp1=0.9]':0.05,C:0.15);")

    output_nwk = tmp_path / "qCF.cf.tree"
    _map_qcf_to_tree(ref_nwk, wastral_nwk, output_nwk, lpp=True)

    out = output_nwk.read_text()
    assert "100/0.5/0.9" in out


def test_map_qcf_canonical_bipartition_matching(tmp_path: Path) -> None:
    """Canonical bipartition matching: both (A,B) and (C,D) sides get annotated
    regardless of rooting differences between ref and wastral trees."""
    from phyloai.tree.cf import _map_qcf_to_tree

    ref_nwk = tmp_path / "ref.nwk"
    ref_nwk.write_text("((A:0.1,B:0.1)100:0.05,(C:0.1,D:0.1):0.05);")

    # wASTRAL labels the (C,D) side with q1=0.6 (different rooting)
    wastral_nwk = tmp_path / "wastral.tre"
    wastral_nwk.write_text(
        "((A:0.1,B:0.1):0.05,(C:0.1,D:0.1)'[q1=0.6;q2=0.2;q3=0.2]':0.05);"
    )

    output_nwk = tmp_path / "qCF_out.tre"
    _map_qcf_to_tree(ref_nwk, wastral_nwk, output_nwk)

    out = output_nwk.read_text()
    assert "100/0.6" in out
    # Verify branch lengths untouched
    assert "0.05" in out
