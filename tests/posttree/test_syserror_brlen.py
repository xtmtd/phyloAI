from __future__ import annotations

import csv
import json
from io import StringIO
from pathlib import Path

import pytest
from Bio import Phylo

from phyloai.posttree.syserror_brlen import (
    TreeRecord,
    _branch_length,
    _canonical_split,
    _internal_rows,
    _is_rooted_representation,
    _missing_length_warning,
    _node_to_tip_rows,
    _resolve_endpoint,
    _terminal_rows,
    _total_rows,
    run_brlen,
    run_label_nodes,
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


def test_unrooted_complement_side_match_is_warned(tmp_path) -> None:
    tree_file = tmp_path / "unrooted.nwk"
    tree_file.write_text("((A:1,B:1):2,C:1,D:1);")
    map_file = tmp_path / "nodes.map.txt"
    map_file.write_text("CD:C,D\n")
    result = run_brlen(
        tree=tree_file, tree_dir=None, mode="node-to-node",
        map_file=map_file, node1="CD", node2="CD",
        output_dir=tmp_path / "out", quiet=True,
    )
    assert any("matched complement side of a split" in warning for warning in result["data"]["warnings"])


def test_map_overrides_internal_node_label_and_single_taxon_is_tip(tmp_path) -> None:
    tree = _tree("((A:1,B:1)N01:2,C:1);")
    endpoint = _resolve_endpoint(tree, "N01", {"N01": frozenset({"A"})})
    assert endpoint.node_type == "tip"


def test_node_to_tip_without_map_is_rooted_labeled_descendants_only() -> None:
    tree = _tree("((A:1,B:2)N01:3,C:4);")
    rows = _node_to_tip_rows(TreeRecord("x", tree), "N01", None, None)
    assert [row["tip"] for row in rows] == ["A", "B"]


def test_unrooted_labeled_node_to_tip_requires_tip_or_map() -> None:
    with pytest.raises(ValueError, match="--tip1 or --map required"):
        _node_to_tip_rows(TreeRecord("x", _tree("(A:1,B:1,C:1)N01;")), "N01", None, None)


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
    result = run_brlen(
        tree=None, tree_dir=trees, mode="tip-to-tip", tip1="A", tip2="Missing",
        table_format="tsv", output_dir=tmp_path / "out", quiet=True,
    )
    assert "Missing" in result["data"]["warnings"][0]
    assert (tmp_path / "out" / "tables" / "tip_to_tip.tsv").read_text().startswith("tree_file\ttip1\ttip2\tdistance")


def test_single_endpoint_missing_tip_is_an_error(tmp_path) -> None:
    tree_file = tmp_path / "one.tre"
    tree_file.write_text("(A:1,B:1);")
    with pytest.raises(ValueError, match="Missing"):
        run_brlen(tree=tree_file, tree_dir=None, mode="tip-to-tip", tip1="A", tip2="Missing", output_dir=tmp_path / "out", quiet=True)


def test_batch_endpoint_failure_counts_into_skipped(tmp_path) -> None:
    trees = tmp_path / "trees"
    trees.mkdir()
    (trees / "a.tre").write_text("(A:1,B:1);")
    (trees / "bad.txt").write_text("((A,B),C;")
    result = run_brlen(
        tree=None, tree_dir=trees, mode="tip-to-tip", tip1="A", tip2="Missing",
        output_dir=tmp_path / "out", quiet=True,
    )
    assert any("Missing" in warning for warning in result["data"]["warnings"])
    assert result["key_results"]["n_trees_skipped"] == 2


def test_dry_run_validates_and_writes_nothing(tmp_path) -> None:
    tree_file = tmp_path / "tree.nwk"
    tree_file.write_text("(A:1,B:1,C:1);")
    result = run_brlen(
        tree=tree_file, tree_dir=None, mode="patristic",
        output_dir=tmp_path / "out", dry_run=True, quiet=True,
    )
    assert result["status"] == "success"
    assert "estimated patristic rows: 3" in result["data"]["warnings"]
    assert result["data"]["output_files"] == {}
    assert not (tmp_path / "out").exists()


def test_dry_run_single_endpoint_error_is_rejected(tmp_path) -> None:
    tree_file = tmp_path / "one.tre"
    tree_file.write_text("(A:1,B:1);")
    with pytest.raises(ValueError, match="Missing"):
        run_brlen(
            tree=tree_file, tree_dir=None, mode="tip-to-tip", tip1="A", tip2="Missing",
            output_dir=tmp_path / "out", dry_run=True, quiet=True,
        )


def test_dry_run_batch_endpoint_warns_and_counts_skipped(tmp_path) -> None:
    trees = tmp_path / "trees"
    trees.mkdir()
    (trees / "a.tre").write_text("(A:1,B:1);")
    result = run_brlen(
        tree=None, tree_dir=trees, mode="tip-to-tip", tip1="A", tip2="Missing",
        output_dir=tmp_path / "out", dry_run=True, quiet=True,
    )
    assert result["status"] == "success"
    assert any("Missing" in warning for warning in result["data"]["warnings"])
    assert result["key_results"]["n_trees_skipped"] == 1
    assert not (tmp_path / "out").exists()


def test_dry_run_still_validates_mode_grammar(tmp_path) -> None:
    tree_file = tmp_path / "tree.nwk"
    tree_file.write_text("(A:1,B:1);")
    with pytest.raises(ValueError, match="unknown mode"):
        run_brlen(tree=tree_file, tree_dir=None, mode="bogus", dry_run=True, output_dir=tmp_path / "out", quiet=True)


def test_missing_lengths_and_patristic_estimate_are_warnings(tmp_path) -> None:
    tree_file = tmp_path / "tree.nwk"
    tree_file.write_text("(A,B,C);")
    result = run_brlen(tree=tree_file, tree_dir=None, mode="patristic", output_dir=tmp_path / "out", quiet=True)
    warnings = result["data"]["warnings"]
    assert any("all branch lengths are missing" in warning for warning in warnings)
    assert any("estimated patristic rows: 3" in warning for warning in warnings)
    assert result["key_results"]["n_trees"] == 1
    assert result["key_results"]["n_trees_skipped"] == 0


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


def test_label_nodes_strips_support_values(tmp_path) -> None:
    tree_file = tmp_path / "support.nwk"
    tree_file.write_text("((A:1,B:1)100:2,C:1);")
    run_label_nodes(tree_file, output_dir=tmp_path / "out", quiet=True)
    labeled = (tmp_path / "out" / "support.labeled.nwk").read_text()
    assert "N1" in labeled and "N2" in labeled
    assert "N1100" not in labeled and "N2100" not in labeled


def test_label_nodes_command_always_records_output_dir(tmp_path, monkeypatch) -> None:
    tree_file = tmp_path / "species.nwk"
    tree_file.write_text("((A:1,B:1):2,C:1);")
    monkeypatch.chdir(tmp_path)
    result = run_label_nodes(tree_file, output_dir=Path("runs/posttree/syserror/brlen/label_nodes"), quiet=True)
    assert "-o" in result["command"]
    assert "runs/posttree/syserror/brlen/label_nodes" in result["command"]


def test_label_nodes_preserves_branch_length_precision(tmp_path) -> None:
    tree_file = tmp_path / "precise.nwk"
    tree_file.write_text("((A:1.23456789,B:0.000001):0.123456789012,C:1);")
    run_label_nodes(tree_file, output_dir=tmp_path / "out", quiet=True)
    labeled = (tmp_path / "out" / "precise.labeled.nwk").read_text()
    assert "1.23456789" in labeled
    assert "0.123456789012" in labeled


def test_command_omits_defaults_but_keeps_required(tmp_path) -> None:
    tree_file = tmp_path / "one.tre"
    tree_file.write_text("((A:1,B:2):3,C:4);")
    result = run_brlen(
        tree=tree_file, tree_dir=None, mode="total",
        output_dir=tmp_path / "out", overwrite=True, quiet=True,
    )
    command = result["command"]
    assert "--mode total" in command
    assert "-o" in command
    assert "--max-rows" not in command
    assert "--table-format" not in command
    assert "-t " not in command
    assert "--overwrite" in command


def test_command_records_non_default_flags(tmp_path) -> None:
    tree_file = tmp_path / "one.tre"
    tree_file.write_text("(A:1,B:1,C:1,D:1);")
    result = run_brlen(
        tree=tree_file, tree_dir=None, mode="patristic", max_rows=0, table_format="tsv",
        output_dir=tmp_path / "out", quiet=True,
    )
    command = result["command"]
    assert "--max-rows 0" in command
    assert "--table-format tsv" in command


def test_label_output_isolated_from_main_result(tmp_path) -> None:
    tree_file = tmp_path / "species.nwk"
    tree_file.write_text("((A:1,B:1):2,C:1);")
    main_dir = tmp_path / "brlen"
    run_brlen(tree=tree_file, tree_dir=None, mode="total", output_dir=main_dir, quiet=True)
    run_label_nodes(tree_file, output_dir=main_dir / "label_nodes", quiet=True)
    assert (main_dir / "result.json").exists()
    assert (main_dir / "label_nodes" / "result.json").exists()


def test_map_file_that_is_a_directory_is_rejected(tmp_path) -> None:
    tree_file = tmp_path / "one.tre"
    tree_file.write_text("((A:1,B:1):1,C:1);")
    with pytest.raises(ValueError, match="map file not readable"):
        run_brlen(
            tree=tree_file, tree_dir=None, mode="node-to-node",
            map_file=tmp_path, node1="N1", node2="N2",
            output_dir=tmp_path / "out", quiet=True,
        )


def test_map_file_with_non_utf8_bytes_is_an_input_error(tmp_path) -> None:
    tree_file = tmp_path / "one.tre"
    tree_file.write_text("((A:1,B:1):1,C:1);")
    map_file = tmp_path / "bad.map"
    map_file.write_bytes(b"N1:\xff\xfe\n")
    with pytest.raises(ValueError, match="cannot read"):
        run_brlen(
            tree=tree_file, tree_dir=None, mode="node-to-node",
            map_file=map_file, node1="N1", node2="N2",
            output_dir=tmp_path / "out", quiet=True,
        )


def test_batch_warning_identifies_the_offending_tree(tmp_path) -> None:
    trees = tmp_path / "trees"
    trees.mkdir()
    (trees / "ok.tre").write_text("((A:1,B:1):1,C:1);")
    (trees / "bad.tre").write_text("(A:1,B:1,C:1);")
    map_file = tmp_path / "nodes.map"
    map_file.write_text("N1:A,B\n")
    result = run_brlen(
        tree=None, tree_dir=trees, mode="node-to-node",
        map_file=map_file, node1="N1", node2="Missing",
        output_dir=tmp_path / "out", quiet=True,
    )
    warnings = result["data"]["warnings"]
    assert any("tree bad.tre:" in warning for warning in warnings)
    assert any("tree ok.tre:" in warning for warning in warnings)
    assert result["key_results"]["n_trees_skipped"] == 2
