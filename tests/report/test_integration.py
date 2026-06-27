"""End-to-end integration tests for phyloai report."""
from __future__ import annotations

import json
from pathlib import Path

from click.testing import CliRunner

from phyloai.cli.commands.report import report


def _make_result_json(path: Path, step_command: str, extra: dict | None = None) -> None:
    data = {
        "status": "success",
        "command": step_command,
        "wall_time": 10.0,
        "tool_versions": {"mafft": "7.526"},
        "params": {"method": "linsi", "seq_type": "AA", "threads": 8, "seq_dir": "./raw"},
        "key_results": {"n_aligned": 100, "n_skipped": 0},
        "error": None,
        "data": {"output_files": {}},
    }
    if extra:
        data.update(extra)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data))


class TestEndToEnd:
    def test_module_run_single_step(self, tmp_path):
        run_dir = tmp_path / "runs" / "pretree"
        step_dir = run_dir / "2-align"
        _make_result_json(
            step_dir / "result.json",
            "phyloai pretree align --seq-dir ./raw --method linsi --threads 4",
        )

        runner = CliRunner()
        result = runner.invoke(report, ["--run-dir", str(run_dir)])
        assert result.exit_code == 0, f"CLI failed: {result.output}"

        report_dir = run_dir / "report"
        assert (report_dir / "report.json").exists()
        assert (report_dir / "report.html").exists()

        rj = json.loads((report_dir / "report.json").read_text())
        assert rj["status"] == "complete"
        assert rj["run_mode"] == "module"
        assert len(rj["steps"]) == 1
        assert rj["steps"][0]["step_id"] == "pretree.align"
        assert "methods_paragraph" in rj
        assert len(rj["methods_paragraph"]) > 0

    def test_module_run_multi_step(self, tmp_path):
        run_dir = tmp_path / "runs" / "pretree"
        _make_result_json(
            run_dir / "2-align" / "result.json",
            "phyloai pretree align --seq-dir ./raw --method linsi",
        )
        _make_result_json(
            run_dir / "4-trim" / "result.json",
            "phyloai pretree trim --msa-dir ./aligned --tool bmge",
            extra={"tool_versions": {"bmge": "1.12"}},
        )

        runner = CliRunner()
        result = runner.invoke(report, ["--run-dir", str(run_dir)])
        assert result.exit_code == 0

        rj = json.loads((run_dir / "report" / "report.json").read_text())
        assert len(rj["steps"]) == 2
        step_ids = {s["step_id"] for s in rj["steps"]}
        assert step_ids == {"pretree.align", "pretree.trim"}

    def test_failed_step_included(self, tmp_path):
        run_dir = tmp_path / "runs" / "pretree"
        align_dir = run_dir / "2-align"
        align_dir.mkdir(parents=True)
        (align_dir / "result.json").write_text(json.dumps({
            "status": "error",
            "command": "phyloai pretree align --seq-dir ./raw",
            "wall_time": 0.5,
            "tool_versions": {},
            "params": {},
            "key_results": {},
            "error": "MAFFT returned exit code 1",
            "data": {"output_files": {}},
        }))

        runner = CliRunner()
        result = runner.invoke(report, ["--run-dir", str(run_dir)])
        assert result.exit_code == 0

        rj = json.loads((run_dir / "report" / "report.json").read_text())
        assert rj["status"] == "failed"
        assert rj["steps"][0]["error"] == "MAFFT returned exit code 1"
        assert rj["steps"][0]["methods_text"] == ""
        assert rj["methods_paragraph"] == ""

    def test_overwrite_protection(self, tmp_path):
        run_dir = tmp_path / "runs" / "pretree"
        _make_result_json(
            run_dir / "2-align" / "result.json",
            "phyloai pretree align --seq-dir ./raw",
        )

        runner = CliRunner()
        r1 = runner.invoke(report, ["--run-dir", str(run_dir)])
        assert r1.exit_code == 0

        r2 = runner.invoke(report, ["--run-dir", str(run_dir)])
        assert r2.exit_code != 0
        assert "overwrite" in r2.output.lower()

        r3 = runner.invoke(report, ["--run-dir", str(run_dir), "--overwrite"])
        assert r3.exit_code == 0

    def test_custom_output_dir(self, tmp_path):
        run_dir = tmp_path / "runs" / "pretree"
        _make_result_json(
            run_dir / "2-align" / "result.json",
            "phyloai pretree align --seq-dir ./raw",
        )
        out_dir = tmp_path / "my-reports" / "run1"

        runner = CliRunner()
        result = runner.invoke(report, ["--run-dir", str(run_dir), "-o", str(out_dir)])
        assert result.exit_code == 0
        assert (out_dir / "report.json").exists()
        assert (out_dir / "report.html").exists()
