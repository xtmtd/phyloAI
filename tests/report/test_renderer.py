"""Smoke tests for phyloai.report.renderer."""
from __future__ import annotations

import json

from phyloai.report.renderer import render_html


class TestRenderHtml:
    def test_produces_valid_html(self, tmp_path):
        report = {
            "phyloai_version": "0.1.0",
            "generated_at": "2026-06-27T00:00:00Z",
            "run_dir": "/tmp/runs/pretree",
            "run_mode": "module",
            "status": "complete",
            "pipeline_summary": {
                "n_steps_total": 1, "n_steps_success": 1,
                "n_steps_failed": 0, "n_steps_skipped": 0, "total_wall_time": 100.0,
            },
            "steps": [{
                "step_id": "pretree.align",
                "command": "phyloai pretree align --seq-dir ./raw --method linsi",
                "status": "success", "wall_time": 31.4,
                "tool_versions": {"mafft": "7.526"},
                "params": {"method": "linsi", "seq_dir": "./raw", "threads": 8, "seq_type": "AA", "backtrans": False},
                "key_results": {"n_aligned": 100, "n_skipped": 0, "mean_alignment_length": 500.0},
                "methods_text": "Multiple sequence alignments were performed using MAFFT v7.526...",
                "output_files": {}, "warnings": [], "error": None,
            }],
            "methods_paragraph": "Multiple sequence alignments were performed using MAFFT v7.526...",
            "figures_index": [], "tables_index": [],
        }
        output_dir = tmp_path / "report"
        output_dir.mkdir()
        (output_dir / "report.json").write_text(json.dumps(report))

        result_path = render_html(report, output_dir)
        assert result_path.exists()
        content = result_path.read_text()
        assert "<!DOCTYPE html>" in content
        assert "PhyloAI Report" in content
        assert "pretree.align" in content
        assert "MAFFT" in content

    def test_with_figures(self, tmp_path):
        pdf = tmp_path / "corr.pdf"
        pdf.write_text("fake pdf")
        report = {
            "phyloai_version": "0.1.0", "generated_at": "2026-06-27T00:00:00Z",
            "run_dir": str(tmp_path), "run_mode": "module", "status": "complete",
            "pipeline_summary": {"n_steps_total": 1, "n_steps_success": 1, "n_steps_failed": 0, "n_steps_skipped": 0, "total_wall_time": 10.0},
            "steps": [{"step_id": "pretree.metrics", "command": "...", "status": "success", "wall_time": 10.0,
                        "tool_versions": {}, "params": {}, "key_results": {}, "methods_text": "Metrics.",
                        "output_files": {"correlation_heatmap": {"path": str(pdf)}}, "warnings": [], "error": None}],
            "methods_paragraph": "Metrics.",
            "figures_index": [{"figure_id": "Fig-3.1", "step_id": "pretree.metrics", "caption": "Correlation heatmap",
                               "path": str(pdf), "type": "pdf"}],
            "tables_index": [],
        }
        output_dir = tmp_path / "report"
        output_dir.mkdir()
        (output_dir / "report.json").write_text(json.dumps(report))
        result_path = render_html(report, output_dir)
        content = result_path.read_text()
        assert "Fig-3.1" in content
        assert '<object' in content

    def test_with_failed_step(self, tmp_path):
        report = {
            "phyloai_version": "0.1.0", "generated_at": "2026-06-27T00:00:00Z",
            "run_dir": str(tmp_path), "run_mode": "module", "status": "partial",
            "pipeline_summary": {"n_steps_total": 2, "n_steps_success": 1, "n_steps_failed": 1, "n_steps_skipped": 0, "total_wall_time": 10.0},
            "steps": [
                {"step_id": "pretree.align", "command": "...", "status": "success", "wall_time": 3.0,
                 "tool_versions": {}, "params": {}, "key_results": {}, "methods_text": "Align.", "output_files": {}, "warnings": [], "error": None},
                {"step_id": "pretree.trim", "command": "...", "status": "error", "wall_time": 0.1,
                 "tool_versions": {}, "params": {}, "key_results": {}, "methods_text": "", "output_files": {}, "warnings": [], "error": "trimAl failed"},
            ],
            "methods_paragraph": "Align.", "figures_index": [], "tables_index": [],
        }
        output_dir = tmp_path / "report"
        output_dir.mkdir()
        (output_dir / "report.json").write_text(json.dumps(report))
        result_path = render_html(report, output_dir)
        content = result_path.read_text()
        assert 'class="failed"' in content
        assert "trimAl failed" in content

    def test_output_files_index_includes_non_figure_files(self, tmp_path):
        fasta_path = tmp_path / "matrix.fa"
        fasta_path.write_text(">seq\nATCG")
        nwk_path = tmp_path / "tree.nwk"
        nwk_path.write_text("(A,B);")
        report = {
            "phyloai_version": "0.1.0", "generated_at": "2026-06-27T00:00:00Z",
            "run_dir": str(tmp_path), "run_mode": "module", "status": "complete",
            "pipeline_summary": {"n_steps_total": 1, "n_steps_success": 1, "n_steps_failed": 0, "n_steps_skipped": 0, "total_wall_time": 10.0},
            "steps": [{"step_id": "pretree.concat", "command": "...", "status": "success", "wall_time": 1.0,
                        "tool_versions": {}, "params": {}, "key_results": {},
                        "methods_text": "Concat.", "warnings": [], "error": None,
                        "output_files": {
                            "matrix": {"path": str(fasta_path), "description": "Supermatrix"},
                            "tree": {"path": str(nwk_path), "description": "Species tree"},
                        }}],
            "methods_paragraph": "Concat.", "figures_index": [], "tables_index": [],
        }
        output_dir = tmp_path / "report"
        output_dir.mkdir()
        (output_dir / "report.json").write_text(json.dumps(report))
        result_path = render_html(report, output_dir)
        content = result_path.read_text()
        assert "matrix.fa" in content
        assert "tree.nwk" in content
        assert "Supermatrix" in content
        assert "Species tree" in content

    def test_tsv_embedded_with_tabs(self, tmp_path):
        tsv = tmp_path / "data.tsv"
        tsv.write_text("col1\tcol2\nval1\tval2\n")
        report = {
            "phyloai_version": "0.1.0", "generated_at": "2026-06-27T00:00:00Z",
            "run_dir": str(tmp_path), "run_mode": "module", "status": "complete",
            "pipeline_summary": {"n_steps_total": 1, "n_steps_success": 1, "n_steps_failed": 0, "n_steps_skipped": 0, "total_wall_time": 1.0},
            "steps": [{"step_id": "x", "command": "...", "status": "success", "wall_time": 1.0,
                        "tool_versions": {}, "params": {}, "key_results": {}, "methods_text": "",
                        "output_files": {"t": {"path": str(tsv)}}, "warnings": [], "error": None}],
            "methods_paragraph": "", "figures_index": [], "tables_index": [],
        }
        output_dir = tmp_path / "report"
        output_dir.mkdir()
        result_path = render_html(report, output_dir)
        content = result_path.read_text()
        # TSV should parse as 2 columns, not 1
        assert "col1" in content
        assert "col2" in content
