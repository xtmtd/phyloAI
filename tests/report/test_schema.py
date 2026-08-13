"""Tests for phyloai.report.schema."""
from __future__ import annotations

from pathlib import Path

from phyloai.report.schema import (
    ReportRecord,
    ReportStep,
    assemble_report,
    build_figures_index,
    build_tables_index,
)


class TestReportStep:
    def test_create_minimal(self):
        step = ReportStep(
            step_id="pretree.align",
            command="phyloai pretree align --seq-dir ./raw",
            status="success",
            wall_time=31.4,
            tool_versions={"mafft": "7.526"},
            params={"method": "linsi", "threads": 8},
            key_results={"n_aligned": 100},
            methods_text="Multiple sequence alignments were performed...",
        )
        assert step.step_id == "pretree.align"
        assert step.status == "success"
        assert step.error is None

    def test_failed_step(self):
        step = ReportStep(
            step_id="pretree.align",
            command="phyloai pretree align ...",
            status="error",
            wall_time=0.5,
            tool_versions={},
            params={},
            key_results={},
            methods_text="",
            error="MAFFT not found",
        )
        assert step.status == "error"
        assert step.error == "MAFFT not found"
        assert step.methods_text == ""

    def test_to_dict(self):
        step = ReportStep(
            step_id="pretree.align",
            command="phyloai pretree align --seq-dir ./raw --method linsi",
            status="success",
            wall_time=31.4,
            tool_versions={"mafft": "7.526"},
            params={"method": "linsi", "seq_dir": "./raw", "threads": 8},
            key_results={"n_aligned": 100, "n_skipped": 0},
            methods_text="Alignments were performed...",
            output_files={"aligned_sequences": {"path": "/tmp/out.fa", "description": "Aligned sequences"}},
            warnings=[],
        )
        d = step.to_dict()
        assert d["step_id"] == "pretree.align"
        assert d["status"] == "success"
        assert d["methods_text"] == "Alignments were performed..."
        assert "aligned_sequences" in d["output_files"]


class TestReportRecord:
    def test_create_minimal(self):
        record = ReportRecord(
            run_dir=Path("/tmp/runs/pretree"),
            run_mode="module",
            status="complete",
        )
        assert record.run_mode == "module"
        assert record.status == "complete"
        assert record.steps == []
        assert record.figures_index == []
        assert record.tables_index == []

    def test_to_dict(self):
        record = ReportRecord(
            run_dir=Path("/tmp/runs/pretree"),
            run_mode="module",
            status="complete",
            steps=[
                ReportStep(
                    step_id="pretree.align",
                    command="phyloai pretree align ...",
                    status="success",
                    wall_time=10.0,
                    tool_versions={"mafft": "7.0"},
                    params={},
                    key_results={"n_aligned": 5},
                    methods_text="Test methods.",
                )
            ],
            methods_paragraph="Test methods.",
            pipeline_summary=None,
        )
        d = record.to_dict()
        assert d["run_mode"] == "module"
        assert d["status"] == "complete"
        assert len(d["steps"]) == 1
        assert d["methods_paragraph"] == "Test methods."
        assert "phyloai_version" in d
        assert "generated_at" in d
        # methods_blocks should be generated from successful steps
        assert len(d["methods_blocks"]) == 1
        assert d["methods_blocks"][0] == {"step_id": "pretree.align", "text": "Test methods.", "step_index": 0}

    def test_to_dict_excludes_failed_from_blocks(self):
        record = ReportRecord(
            run_dir=Path("/tmp/runs"),
            run_mode="module",
            status="partial",
            steps=[
                ReportStep(step_id="a", command="...", status="success", wall_time=1.0,
                           tool_versions={}, params={}, key_results={}, methods_text="Good."),
                ReportStep(step_id="b", command="...", status="error", wall_time=0.1,
                           tool_versions={}, params={}, key_results={}, methods_text="Bad.",
                           error="failed"),
            ],
        )
        d = record.to_dict()
        assert len(d["methods_blocks"]) == 1
        assert d["methods_blocks"][0]["step_id"] == "a"


class TestBuildFiguresIndex:
    def test_extracts_pdf_and_png(self):
        steps = [
            {
                "step_id": "pretree.metrics",
                "output_files": {
                    "correlation_heatmap": {"path": "/tmp/corr.pdf", "description": "Correlation heatmap"},
                    "metrics_table": {"path": "/tmp/metrics.csv", "description": "Metrics table"},
                    "distribution_plot": {"path": "/tmp/dist.png", "description": "Distribution plot"},
                },
            },
        ]
        figures = build_figures_index(steps)
        assert len(figures) == 2
        paths = {f["path"] for f in figures}
        assert "/tmp/corr.pdf" in paths
        assert "/tmp/dist.png" in paths

    def test_skips_non_figure_types(self):
        steps = [
            {
                "step_id": "pretree.align",
                "output_files": {
                    "aligned": {"path": "/tmp/aligned.fa"},
                    "log": {"path": "/tmp/run.log"},
                },
            },
        ]
        figures = build_figures_index(steps)
        assert figures == []

    def test_figure_numbering_by_phase(self):
        steps = [
            {"step_id": "pretree.metrics", "output_files": {"a": {"path": "/tmp/a.pdf"}}},
            {"step_id": "tree.ml.iqtree", "output_files": {"b": {"path": "/tmp/b.pdf"}}},
            {"step_id": "posttree.topology", "output_files": {"c": {"path": "/tmp/c.pdf"}}},
        ]
        figures = build_figures_index(steps)
        assert figures[0]["figure_id"] == "Fig-3.1"
        assert figures[1]["figure_id"] == "Fig-4.1"
        assert figures[2]["figure_id"] == "Fig-5.1"

    def test_sequential_within_phase(self):
        steps = [
            {"step_id": "pretree.metrics", "output_files": {"a": {"path": "/tmp/a.pdf"}, "b": {"path": "/tmp/b.pdf"}}},
        ]
        figures = build_figures_index(steps)
        assert figures[0]["figure_id"] == "Fig-3.1"
        assert figures[1]["figure_id"] == "Fig-3.2"

    def test_missing_description_fallback(self):
        steps = [
            {"step_id": "pretree.metrics", "output_files": {"heatmap": {"path": "/tmp/hm.pdf"}}},
        ]
        figures = build_figures_index(steps)
        assert figures[0]["description"] == "heatmap"

    def test_reads_from_data_output_files(self):
        steps = [
            {
                "step_id": "pretree.metrics",
                "output_files": {},
                "data": {"output_files": {"fig": {"path": "/tmp/f.pdf", "description": "A figure"}}},
            },
        ]
        figures = build_figures_index(steps)
        assert len(figures) == 1
        assert figures[0]["caption"] == "A figure"


class TestBuildTablesIndex:
    def test_extracts_csv_and_tsv(self):
        steps = [
            {
                "step_id": "pretree.metrics",
                "output_files": {
                    "metrics_table": {"path": "/tmp/metrics.csv", "description": "Metrics"},
                    "results": {"path": "/tmp/results.tsv", "description": "Results"},
                    "plot": {"path": "/tmp/plot.pdf"},
                },
            },
        ]
        tables = build_tables_index(steps)
        assert len(tables) == 2
        paths = {t["path"] for t in tables}
        assert "/tmp/metrics.csv" in paths
        assert "/tmp/results.tsv" in paths

    def test_rate_output_indexes_only_rates_csv(self, tmp_path):
        discovered = {
            "run_mode": "module",
            "steps": [
                {
                    "step_id": "posttree.syserror.rate",
                    "command": "phyloai posttree syserror rate --iqtree-rate matrix.rate",
                    "status": "success",
                    "wall_time": 0.0,
                    "tool_versions": {},
                    "params": {},
                    "key_results": {},
                    "error": None,
                    "data": {"output_files": {
                        "rates": {"path": "/tmp/rates.csv"},
                        "slow25_positions": {"path": "/tmp/slow25/positions.txt"},
                        "slow25_matrix": {"path": "/tmp/slow25/matrix.fa"},
                    }},
                },
            ],
            "pipeline_summary": None,
        }

        assert [table["label"] for table in build_tables_index(discovered["steps"])] == ["rates"]
        report = assemble_report(discovered, tmp_path)
        assert set(report["steps"][0]["output_files"]) == {
            "rates", "slow25_positions", "slow25_matrix",
        }


class TestAssembleReport:
    def test_complete_run(self, tmp_path):
        discovered = {
            "run_mode": "module",
            "steps": [
                {
                    "step_id": "pretree.align",
                    "command": "phyloai pretree align --seq-dir ./raw --method linsi",
                    "status": "success",
                    "wall_time": 31.4,
                    "tool_versions": {"mafft": "7.526"},
                    "params": {"method": "linsi", "seq_dir": "./raw", "threads": 4},
                    "key_results": {"n_aligned": 100, "n_skipped": 0},
                    "error": None,
                    "data": {"output_files": {"aligned": {"path": str(tmp_path / "aligned.fa")}}},
                    "methods_text": "MSA was performed using MAFFT v7.526 with L-INS-i.",
                },
            ],
            "pipeline_summary": None,
        }

        report = assemble_report(discovered, tmp_path)
        assert report["status"] == "complete"
        assert report["run_mode"] == "module"
        assert len(report["steps"]) == 1
        assert report["steps"][0]["methods_text"] == "MSA was performed using MAFFT v7.526 with L-INS-i."
        assert "methods_paragraph" in report
        assert "figures_index" in report
        assert "tables_index" in report

    def test_partial_failure(self, tmp_path):
        discovered = {
            "run_mode": "module",
            "steps": [
                {
                    "step_id": "pretree.align",
                    "command": "phyloai pretree align ...",
                    "status": "success",
                    "wall_time": 10.0,
                    "tool_versions": {},
                    "params": {},
                    "key_results": {},
                    "error": None,
                    "data": {"output_files": {}},
                    "methods_text": "Methods for align.",
                },
                {
                    "step_id": "pretree.trim",
                    "command": "phyloai pretree trim ...",
                    "status": "error",
                    "wall_time": 0.1,
                    "tool_versions": {},
                    "params": {},
                    "key_results": {},
                    "error": "trimAl exited with code 1",
                    "data": {},
                },
            ],
            "pipeline_summary": None,
        }

        report = assemble_report(discovered, tmp_path)
        assert report["status"] == "partial"
        assert report["pipeline_summary"]["n_steps_success"] == 1
        assert report["pipeline_summary"]["n_steps_failed"] == 1

    def test_methods_paragraph_excludes_failed(self, tmp_path):
        discovered = {
            "run_mode": "module",
            "steps": [
                {
                    "step_id": "pretree.align",
                    "command": "phyloai pretree align ...",
                    "status": "success",
                    "wall_time": 10.0,
                    "tool_versions": {},
                    "params": {},
                    "key_results": {},
                    "error": None,
                    "data": {"output_files": {}},
                    "methods_text": "Align methods.",
                },
                {
                    "step_id": "pretree.trim",
                    "command": "phyloai pretree trim ...",
                    "status": "error",
                    "wall_time": 0.1,
                    "tool_versions": {},
                    "params": {},
                    "key_results": {},
                    "error": "trimAl failed",
                    "data": {},
                    "methods_text": "Trim methods.",
                },
            ],
            "pipeline_summary": None,
        }

        report = assemble_report(discovered, tmp_path)
        assert "Trim methods" not in report["methods_paragraph"]
        assert "Align methods" in report["methods_paragraph"]

    def test_same_step_different_params_kept(self, tmp_path):
        """Same step_id with different seq_type → both texts kept."""
        discovered = {
            "run_mode": "module",
            "steps": [
                {
                    "step_id": "pretree.convert",
                    "command": "phyloai pretree convert --to fasta --seq-type NT",
                    "status": "success", "wall_time": 1.0,
                    "tool_versions": {}, "params": {"to": "fasta", "seq_type": "NT"},
                    "key_results": {"n_converted": 50}, "error": None,
                    "data": {"output_files": {}},
                    "methods_text": "Converted 50 NT files to FASTA.",
                },
                {
                    "step_id": "pretree.convert",
                    "command": "phyloai pretree convert --to fasta --seq-type AA",
                    "status": "success", "wall_time": 1.0,
                    "tool_versions": {}, "params": {"to": "fasta", "seq_type": "AA"},
                    "key_results": {"n_converted": 30}, "error": None,
                    "data": {"output_files": {}},
                    "methods_text": "Converted 30 AA files to FASTA.",
                },
            ],
            "pipeline_summary": None,
        }
        report = assemble_report(discovered, tmp_path)
        assert "NT files" in report["methods_paragraph"]
        assert "AA files" in report["methods_paragraph"]
