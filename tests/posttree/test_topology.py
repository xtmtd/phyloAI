"""Tests for phyloai.posttree.topology."""
from __future__ import annotations

import shutil
from pathlib import Path

import pytest

from phyloai.posttree.topology import (
    _build_topology_cmd,
    _merge_candidate_trees,
    _parse_user_trees_table,
    _validate_inputs,
    run_topology,
)


# ------------------------------------------------------------------
# _validate_inputs
# ------------------------------------------------------------------

class TestValidateInputs:
    def test_matrix_does_not_exist(self, tmp_path: Path) -> None:
        errs = _validate_inputs(
            matrix=tmp_path / "nope.fa",
            candidate_trees=[tmp_path / "t.nwk"],
            replicates=10000, threads=4,
            overwrite=False, resume=False,
            model_expr="LG+F+R4", partitions=None,
            tool_args=None, guide_tree=None,
        )
        assert any("matrix" in e.lower() for e in errs)

    def test_candidate_tree_empty(self, tmp_path: Path) -> None:
        matrix = tmp_path / "m.fa"
        matrix.write_text(">a\nMKT\n")
        ct = tmp_path / "empty.nwk"
        ct.write_text("")
        errs = _validate_inputs(
            matrix=matrix, candidate_trees=[ct],
            replicates=10000, threads=4,
            overwrite=False, resume=False,
            model_expr="LG+F+R4", partitions=None,
            tool_args=None, guide_tree=None,
        )
        assert any("empty" in e for e in errs)

    def test_replicates_too_low(self, tmp_path: Path) -> None:
        matrix = tmp_path / "m.fa"
        matrix.write_text(">a\nMKT\n")
        ct = tmp_path / "t.nwk"
        ct.write_text("(a,b);\n")
        errs = _validate_inputs(
            matrix=matrix, candidate_trees=[ct],
            replicates=999, threads=4,
            overwrite=False, resume=False,
            model_expr="LG+F+R4", partitions=None,
            tool_args=None, guide_tree=None,
        )
        assert any("replicates" in e.lower() for e in errs)

    def test_overwrite_and_resume_mutually_exclusive(self, tmp_path: Path) -> None:
        matrix = tmp_path / "m.fa"
        matrix.write_text(">a\nMKT\n")
        ct = tmp_path / "t.nwk"
        ct.write_text("(a,b);\n")
        errs = _validate_inputs(
            matrix=matrix, candidate_trees=[ct],
            replicates=10000, threads=4,
            overwrite=True, resume=True,
            model_expr="LG+F+R4", partitions=None,
            tool_args=None, guide_tree=None,
        )
        assert any("mutually exclusive" in e.lower() for e in errs)

    def test_no_model_source(self, tmp_path: Path) -> None:
        matrix = tmp_path / "m.fa"
        matrix.write_text(">a\nMKT\n")
        ct = tmp_path / "t.nwk"
        ct.write_text("(a,b);\n")
        errs = _validate_inputs(
            matrix=matrix, candidate_trees=[ct],
            replicates=10000, threads=4,
            overwrite=False, resume=False,
            model_expr=None, partitions=None,
            tool_args=None, guide_tree=None,
        )
        assert any("model" in e.lower() for e in errs)

    def test_both_model_expr_and_partitions(self, tmp_path: Path) -> None:
        matrix = tmp_path / "m.fa"
        matrix.write_text(">a\nMKT\n")
        ct = tmp_path / "t.nwk"
        ct.write_text("(a,b);\n")
        part = tmp_path / "m.nex"
        part.write_text("#NEXUS\nbegin sets;\ncharset p1 = 1-3;\nend;\n")
        errs = _validate_inputs(
            matrix=matrix, candidate_trees=[ct],
            replicates=10000, threads=4,
            overwrite=False, resume=False,
            model_expr="LG", partitions=str(part),
            tool_args=None, guide_tree=None,
        )
        assert errs == []

    def test_all_valid(self, tmp_path: Path) -> None:
        matrix = tmp_path / "m.fa"
        matrix.write_text(">a\nMKT\n")
        ct = tmp_path / "t.nwk"
        ct.write_text("(a,b);\n")
        errs = _validate_inputs(
            matrix=matrix, candidate_trees=[ct],
            replicates=10000, threads=4,
            overwrite=False, resume=False,
            model_expr="LG+F+R4", partitions=None,
            tool_args=None, guide_tree=None,
        )
        assert errs == []


# ------------------------------------------------------------------
# _merge_candidate_trees
# ------------------------------------------------------------------

class TestMergeCandidateTrees:
    def test_merge_two_files(self, tmp_path: Path) -> None:
        (tmp_path / "h1.nwk").write_text("(A,B);\n")
        (tmp_path / "h2.nwk").write_text("(A,C);\n")
        merged = _merge_candidate_trees(
            [tmp_path / "h1.nwk", tmp_path / "h2.nwk"], tmp_path,
        )
        assert merged.name == "candidate.trees"
        content = merged.read_text()
        assert "(A,B);" in content
        assert "(A,C);" in content

    def test_empty_file_raises(self, tmp_path: Path) -> None:
        (tmp_path / "empty.nwk").write_text("")
        with pytest.raises(ValueError, match="Empty candidate tree file"):
            _merge_candidate_trees([tmp_path / "empty.nwk"], tmp_path)


# ------------------------------------------------------------------
# _build_topology_cmd
# ------------------------------------------------------------------

class TestBuildTopologyCmd:
    def test_basic_command(self, tmp_path: Path) -> None:
        cmd = _build_topology_cmd(
            executable="iqtree3",
            matrix=tmp_path / "matrix.fa",
            candidate_trees=tmp_path / "trees",
            prefix="matrix",
            model_expr="LG+F+R4",
            partitions=None,
            guide_tree=None,
            replicates=10000, threads=20,
            tool_args=None,
        )
        assert cmd[0] == "iqtree3"
        assert "-s" in cmd
        assert "-z" in cmd
        assert "-m" in cmd and "LG+F+R4" in cmd
        assert "-n" in cmd and "0" in cmd
        assert "-zb" in cmd and "10000" in cmd
        assert "-zw" in cmd
        assert "-au" in cmd
        assert "-T" in cmd and "20" in cmd

    def test_partitions_mode(self, tmp_path: Path) -> None:
        cmd = _build_topology_cmd(
            executable="iqtree3",
            matrix=tmp_path / "matrix.fa",
            candidate_trees=tmp_path / "trees",
            prefix="m",
            model_expr=None, partitions="m.best_model.nex",
            guide_tree=None,
            replicates=10000, threads=4,
            tool_args=None,
        )
        assert "-p" in cmd
        assert "-m" not in cmd

    def test_guide_tree(self, tmp_path: Path) -> None:
        cmd = _build_topology_cmd(
            executable="iqtree3",
            matrix=tmp_path / "matrix.fa",
            candidate_trees=tmp_path / "trees",
            prefix="m",
            model_expr="LG+C20+F+R4", partitions=None,
            guide_tree="guide.nwk",
            replicates=10000, threads=4,
            tool_args=None,
        )
        assert "-ft" in cmd

    def test_suppress_threads_via_tool_args(self, tmp_path: Path) -> None:
        cmd = _build_topology_cmd(
            executable="iqtree3",
            matrix=tmp_path / "matrix.fa",
            candidate_trees=tmp_path / "trees",
            prefix="m",
            model_expr="LG+F+R4", partitions=None, guide_tree=None,
            replicates=10000, threads=20,
            tool_args="-T 30",
        )
        t_indices = [i for i, t in enumerate(cmd) if t == "-T"]
        assert len(t_indices) == 1
        assert cmd[t_indices[0] + 1] == "30"

    def test_suppress_zb_via_tool_args(self, tmp_path: Path) -> None:
        cmd = _build_topology_cmd(
            executable="iqtree3",
            matrix=tmp_path / "matrix.fa",
            candidate_trees=tmp_path / "trees",
            prefix="m",
            model_expr="LG+F+R4", partitions=None, guide_tree=None,
            replicates=10000, threads=4,
            tool_args="-zb 5000",
        )
        zb_indices = [i for i, t in enumerate(cmd) if t == "-zb"]
        assert len(zb_indices) == 1
        assert cmd[zb_indices[0] + 1] == "5000"

    def test_blocked_s_flag_raises(self, tmp_path: Path) -> None:
        with pytest.raises(ValueError, match="-s"):
            _build_topology_cmd(
                executable="iqtree3",
                matrix=tmp_path / "matrix.fa",
                candidate_trees=tmp_path / "trees",
                prefix="m",
                model_expr="LG+F+R4", partitions=None, guide_tree=None,
                replicates=10000, threads=4,
                tool_args="-s other.fa",
            )

    def test_blocked_z_flag_raises(self, tmp_path: Path) -> None:
        with pytest.raises(ValueError, match="-z"):
            _build_topology_cmd(
                executable="iqtree3",
                matrix=tmp_path / "matrix.fa",
                candidate_trees=tmp_path / "trees",
                prefix="m",
                model_expr="LG+F+R4", partitions=None, guide_tree=None,
                replicates=10000, threads=4,
                tool_args="-z other.trees",
            )

    def test_replicates_value(self, tmp_path: Path) -> None:
        cmd = _build_topology_cmd(
            executable="iqtree3",
            matrix=tmp_path / "matrix.fa",
            candidate_trees=tmp_path / "trees",
            prefix="m",
            model_expr="LG+F+R4", partitions=None, guide_tree=None,
            replicates=2000, threads=4,
            tool_args=None,
        )
        zb_idx = cmd.index("-zb")
        assert cmd[zb_idx + 1] == "2000"


# ------------------------------------------------------------------
# _parse_user_trees_table
# ------------------------------------------------------------------

class TestParseUserTreesTable:
    def test_parse_standard_iqtree_output(self, tmp_path: Path) -> None:
        content = """\
USER TREES:

See http://www.iqtree.org/doc/Topology-Tests

Tree      logL    deltaL  bp-RELL    p-KH     p-SH    p-WKH    p-WSH    c-ELW       p-AU
------------------------------------------------------------------------------------------
 1  -21152.617    0.000  0.7110 +  0.7400 +  1.0000 +  0.7380 +  1.0000 +  0.6954 +  0.7939 +
 2  -21158.123    5.506  0.2299 +  0.2590 +  0.1260 +  0.2690 +  0.1330 +  0.2275 +  0.2336 +
 3  -21162.987   10.370  0.0392 +  0.0080 -  0.0070 -  0.0080 -  0.0060 -  0.0404 +  0.0140 -
 4  -21168.456   15.839  0.0199 -  0.0010 -  0.0000 -  0.0010 -  0.0000 -  0.0367 -  0.0030 -

------------------------------------------------------------------------------------------
"""
        p = tmp_path / "test.iqtree"
        p.write_text(content)
        tests, warnings = _parse_user_trees_table(p)
        assert len(tests) == 4
        t1 = tests[0]
        assert t1["tree_id"] == 1
        assert t1["log_likelihood"] == -21152.617
        assert t1["bp_rell"] == 0.7110
        assert t1["bp_rell_sign"] == "+"
        assert t1["p_au"] == 0.7939
        # raw_line preserved
        assert "raw_line" in t1
        assert "1  -21152.617" in t1["raw_line"]

    def test_missing_file(self, tmp_path: Path) -> None:
        tests, warnings = _parse_user_trees_table(tmp_path / "nope.iqtree")
        assert tests == []
        assert len(warnings) > 0

    def test_no_user_trees_section(self, tmp_path: Path) -> None:
        (tmp_path / "test.iqtree").write_text("OTHER SECTION\nNo trees\n")
        tests, warnings = _parse_user_trees_table(tmp_path / "test.iqtree")
        assert tests == []
        assert len(warnings) > 0

    def test_missing_column_stores_none(self, tmp_path: Path) -> None:
        # IQ-TREE may omit some columns; absent values must be None
        content = """\
USER TREES:
Tree      logL    deltaL
 1  -100.0  0.0
 2  -105.0  5.0
"""
        (tmp_path / "test.iqtree").write_text(content)
        tests, _ = _parse_user_trees_table(tmp_path / "test.iqtree")
        assert len(tests) == 2
        assert tests[0]["tree_id"] == 1
        assert tests[0]["p_au"] is None  # column absent entirely


# ------------------------------------------------------------------
# run_topology (dry-run + integration)
# ------------------------------------------------------------------

class TestRunTopology:
    def test_dry_run_single_tree_file(self, tmp_path: Path) -> None:
        from tests.helpers import validate_result_json

        (tmp_path / "matrix.fa").write_text(">a\nMKTLLL\n>b\nMKTLLL\n")
        (tmp_path / "trees").write_text("(a,b);\n")

        result = run_topology(
            matrix=tmp_path / "matrix.fa",
            candidate_trees=[tmp_path / "trees"],
            model_expr="LG+F+R4",
            output_dir=tmp_path / "out",
            dry_run=True,
        )
        validate_result_json(result)
        assert result["status"] == "success"
        assert result["params"]["candidate_trees_mode"] == "tree-list"
        assert isinstance(result["data"]["cmd"], list)

    def test_dry_run_multiple_tree_files(self, tmp_path: Path) -> None:
        (tmp_path / "matrix.fa").write_text(">a\nMKTLLL\n>b\nMKTLLL\n")
        (tmp_path / "h1.nwk").write_text("(a,b);\n")
        (tmp_path / "h2.nwk").write_text("(a,c);\n")

        result = run_topology(
            matrix=tmp_path / "matrix.fa",
            candidate_trees=[tmp_path / "h1.nwk", tmp_path / "h2.nwk"],
            model_expr="LG+F+R4",
            output_dir=tmp_path / "out",
            dry_run=True,
        )
        assert result["params"]["candidate_trees_mode"] == "individual-files"
        assert len(result["params"]["candidate_trees"]) == 2

    def test_validation_error_returns_error_payload(self, tmp_path: Path) -> None:
        from tests.helpers import validate_result_json

        result = run_topology(
            matrix=tmp_path / "nope.fa",
            candidate_trees=[],
            model_expr="LG+F+R4",
            dry_run=True,
        )
        assert result["status"] == "error"
        assert result["error"] is not None
        validate_result_json(result)
        assert result["command"].startswith("phyloai ")
        assert isinstance(result["params"], dict) and result["params"]
        assert isinstance(result["data"], dict)

    @pytest.mark.skipif(
        not shutil.which("iqtree3"),
        reason="iqtree3 not found in PATH",
    )
    def test_run_topology_real_iqtree(self, tmp_path: Path) -> None:
        from tests.helpers import validate_result_json

        matrix = tmp_path / "matrix.fa"
        matrix.write_text(
            ">t1\nMKTLLLTLWVV\n>t2\nMKTLLLTLWVI\n>t3\nMKTLLLSLWVI\n>t4\nMKTLLLTLWVA\n"
        )
        (tmp_path / "t1.nwk").write_text("(t1,t2,(t3,t4));\n")
        (tmp_path / "t2.nwk").write_text("(t1,t3,(t2,t4));\n")
        (tmp_path / "t3.nwk").write_text("(t1,t4,(t2,t3));\n")
        (tmp_path / "t4.nwk").write_text("(t2,t3,(t1,t4));\n")

        out = tmp_path / "out"
        result = run_topology(
            matrix=matrix,
            candidate_trees=[
                tmp_path / "t1.nwk", tmp_path / "t2.nwk",
                tmp_path / "t3.nwk", tmp_path / "t4.nwk",
            ],
            model_expr="LG",
            replicates=1000, output_dir=out, threads=1,
        )

        validate_result_json(result)
        assert result["status"] == "success"
        assert "iqtree3" in result["tool_versions"]
        assert len(result["data"]["tests"]) == 4
        assert result["data"]["merged_candidate_trees"] == "candidate.trees"
        assert (out / "candidate.trees").exists()
        assert (out / "matrix.iqtree").exists()
        assert (out / "matrix.log").exists()
