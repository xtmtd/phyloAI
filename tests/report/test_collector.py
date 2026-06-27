"""Tests for phyloai.report.collector."""
from __future__ import annotations

import json

from phyloai.report.collector import (
    STEP_ORDER,
    discover_steps,
    parse_step_id,
)


class TestParseStepId:
    def test_basic_two_part(self):
        assert parse_step_id("phyloai pretree align --seq-dir ./raw --method linsi") == "pretree.align"

    def test_three_part_filter(self):
        assert parse_step_id("phyloai pretree filter taper --msa-dir ./trimmed") == "pretree.filter.taper"

    def test_tree_ml_iqtree(self):
        assert parse_step_id("phyloai tree ml iqtree --matrix ./concat/matrix.fa --model C20") == "tree.ml.iqtree"

    def test_pipeline_run(self):
        assert parse_step_id("phyloai run --seq-dir ./raw --mode supermatrix") == "run"

    def test_dating_mcmc(self):
        assert parse_step_id("phyloai posttree dating mcmc --tree ./tree.nwk --matrix ./matrix.fa") == "posttree.dating.mcmc"

    def test_with_paths_and_flags(self):
        assert parse_step_id("phyloai pretree trim --msa-dir /abs/path/to/trimmed --tool bmge --bmge-matrix BLOSUM90 --tool-args '-g 0.5'") == "pretree.trim"

    def test_global_flag_before_subcommand(self):
        assert parse_step_id("phyloai --run-dir ./runs pretree align --method linsi") == "pretree.align"
        assert parse_step_id("phyloai --quiet tree ml iqtree --model LG") == "tree.ml.iqtree"
        assert parse_step_id("phyloai --overwrite --quiet pretree filter taper --cutoff 0.1") == "pretree.filter.taper"

    def test_empty_command(self):
        assert parse_step_id("") == "unknown"

    def test_not_phyloai(self):
        assert parse_step_id("some other command") == "unknown"


class TestStepOrder:
    def test_step_order_is_list(self):
        assert isinstance(STEP_ORDER, list)
        assert len(STEP_ORDER) > 10

    def test_common_steps_in_order(self):
        assert "pretree.convert" in STEP_ORDER
        assert "pretree.align" in STEP_ORDER
        assert "tree.ml.iqtree" in STEP_ORDER
        assert "posttree.topology" in STEP_ORDER
        assert STEP_ORDER.index("pretree.convert") < STEP_ORDER.index("pretree.align")
        assert STEP_ORDER.index("pretree.align") < STEP_ORDER.index("tree.ml.iqtree")


class TestDiscoverStepsModule:
    def test_single_step_module(self, tmp_path):
        sub = tmp_path / "2-align"
        sub.mkdir()
        (sub / "result.json").write_text(json.dumps({
            "status": "success",
            "command": "phyloai pretree align --seq-dir ./raw --method linsi --threads 8",
            "wall_time": 31.4,
            "tool_versions": {"mafft": "7.526"},
            "params": {"seq_dir": "./raw", "method": "linsi", "threads": 8},
            "key_results": {"n_aligned": 100},
            "error": None,
            "data": {"output_files": {}},
        }))

        result = discover_steps(sub)
        assert result["run_mode"] == "module"
        assert len(result["steps"]) == 1
        assert result["steps"][0]["step_id"] == "pretree.align"

    def test_multi_step_module(self, tmp_path):
        for name, step_cmd in [("2-align", "phyloai pretree align --some-flag"),
                                ("4-trim", "phyloai pretree trim --some-flag"),
                                ("5-metrics", "phyloai pretree metrics --some-flag")]:
            d = tmp_path / name
            d.mkdir()
            (d / "result.json").write_text(json.dumps({
                "status": "success",
                "command": step_cmd,
                "wall_time": 10.0,
                "tool_versions": {},
                "params": {},
                "key_results": {},
                "error": None,
                "data": {"output_files": {}},
            }))

        result = discover_steps(tmp_path)
        assert result["run_mode"] == "module"
        assert len(result["steps"]) == 3

    def test_no_result_json_error(self, tmp_path):
        try:
            discover_steps(tmp_path)
            assert False, "should have raised"
        except ValueError as e:
            assert "result.json" in str(e).lower()

    def test_excludes_report_and_logs(self, tmp_path):
        sub = tmp_path / "2-align"
        sub.mkdir()
        (sub / "result.json").write_text(json.dumps({
            "status": "success",
            "command": "phyloai pretree align --some-flag",
            "wall_time": 10.0,
            "tool_versions": {},
            "params": {},
            "key_results": {},
            "error": None,
            "data": {"output_files": {}},
        }))
        report_dir = tmp_path / "report"
        report_dir.mkdir()
        (report_dir / "result.json").write_text(json.dumps({"status": "success"}))

        result = discover_steps(tmp_path)
        assert len(result["steps"]) == 1
        assert result["steps"][0]["step_id"] == "pretree.align"


class TestDiscoverStepsPipeline:
    def test_pipeline_detection(self, tmp_path):
        (tmp_path / "result.json").write_text(json.dumps({
            "status": "success",
            "command": "phyloai run --seq-dir ./raw --mode supermatrix",
            "wall_time": 8040.5,
            "tool_versions": {},
            "params": {"mode": "supermatrix", "speed": "normal"},
            "key_results": {"n_input_genes": 200},
            "error": None,
            "data": {"mode": "supermatrix", "speed": "normal"},
        }))

        for name, step_cmd in [("1-convert", "phyloai pretree convert --input ./raw --to fasta"),
                                ("2-align", "phyloai pretree align --seq-dir ./raw --method linsi")]:
            d = tmp_path / name
            d.mkdir()
            (d / "result.json").write_text(json.dumps({
                "status": "success",
                "command": step_cmd,
                "wall_time": 1.0,
                "tool_versions": {},
                "params": {},
                "key_results": {},
                "error": None,
                "data": {"output_files": {}},
            }))

        result = discover_steps(tmp_path)
        assert result["run_mode"] == "pipeline"
        assert result["pipeline_summary"] is not None
        assert result["pipeline_summary"]["mode"] == "supermatrix"
        assert len(result["steps"]) == 2

    def test_pipeline_reads_metadata_only(self, tmp_path):
        """data.steps is NOT required — detection is purely filesystem-based."""
        (tmp_path / "result.json").write_text(json.dumps({
            "status": "success",
            "command": "phyloai run --mode supermatrix",
            "wall_time": 1.0,
            "tool_versions": {},
            "params": {},
            "key_results": {},
            "error": None,
            "data": {"mode": "supermatrix"},  # no data.steps
        }))
        sub = tmp_path / "1-convert"
        sub.mkdir()
        (sub / "result.json").write_text(json.dumps({
            "status": "success",
            "command": "phyloai pretree convert --to fasta",
            "wall_time": 1.0,
            "tool_versions": {},
            "params": {},
            "key_results": {},
            "error": None,
            "data": {"output_files": {}},
        }))

        result = discover_steps(tmp_path)
        assert result["run_mode"] == "pipeline"
        assert len(result["steps"]) == 1
