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
        "--wastral-path", "--overwrite", "--dry-run", "--quiet", "--lpp",
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
    """Missing --cf flag produces non-zero exit."""
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


def test_tree_cf_lpp_non_qcf_writes_error_outputs(tmp_path: Path) -> None:
    """--lpp with --cf gcf exits 1 and writes result.json."""
    import json

    ref_tree = tmp_path / "ref.nwk"
    ref_tree.write_text("(A,B);")
    gene_trees = tmp_path / "trees.nwk"
    gene_trees.write_text("(A,B);\n")
    out_dir = tmp_path / "out"

    result = CliRunner().invoke(cli, [
        "tree", "cf", "--cf", "gcf",
        "--ref-tree", str(ref_tree),
        "--tree", str(gene_trees),
        "--lpp",
        "-o", str(out_dir),
    ])
    assert result.exit_code == 1

    result_json = out_dir / "result.json"
    assert result_json.exists()
    payload = json.loads(result_json.read_text())
    assert payload["status"] == "error"
    assert payload["error"]
    assert "wall_time" in payload
    assert "tool_versions" in payload
    assert "key_results" in payload
    assert "data" in payload
    assert payload["params"]["lpp"] is True


def test_tree_cf_validation_error_writes_error_outputs(tmp_path: Path) -> None:
    """Validation errors (e.g. scf with tree) exit !=0 and write result.json."""
    import json

    ref_tree = tmp_path / "ref.nwk"
    ref_tree.write_text("(A,B);")
    gene_trees = tmp_path / "trees.nwk"
    gene_trees.write_text("(A,B);\n")
    out_dir = tmp_path / "out"

    # scf mode does not accept --tree
    result = CliRunner().invoke(cli, [
        "tree", "cf", "--cf", "scf",
        "--ref-tree", str(ref_tree),
        "--tree", str(gene_trees),
        "-o", str(out_dir),
    ])
    assert result.exit_code == 1

    result_json = out_dir / "result.json"
    assert result_json.exists()
    payload = json.loads(result_json.read_text())
    assert payload["status"] == "error"
    assert "not needed" in payload["error"]
    assert "wall_time" in payload
    assert "key_results" in payload
    assert "data" in payload
