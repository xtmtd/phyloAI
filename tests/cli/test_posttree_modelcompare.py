"""CLI tests for phyloai posttree modelcompare."""
from __future__ import annotations

import json
from pathlib import Path

from click.testing import CliRunner

from phyloai.cli.main import cli


def _make_sitelogl_dir(tmp_path: Path) -> Path:
    model_dir = tmp_path / "model1"
    model_dir.mkdir()
    for i in (1, 2):
        p = model_dir / f"chain{i}.sitelogl"
        p.write_text(
            "site\tlogl\tvar\tlogcpo\tess\tlogpostmeanl\tess\n"
            "1\t-1.0\t0.1\t-1.0\t60.0\t-1.0\t65.0\n"
            "2\t-1.5\t0.2\t-1.5\t60.0\t-1.5\t65.0\n"
        )
    return model_dir


def test_pb_output_conflict_keeps_existing_result(tmp_path: Path) -> None:
    out = tmp_path / "out"
    out.mkdir()
    existing = out / "result.json"
    existing.write_text(json.dumps({"status": "success", "key_results": {"old": True}}))

    result = CliRunner().invoke(cli, [
        "posttree", "modelcompare", "pb",
        "--sitelogl-dir", str(_make_sitelogl_dir(tmp_path)),
        "--output-dir", str(out), "--quiet",
    ])
    assert result.exit_code != 0
    assert json.loads(existing.read_text())["key_results"]["old"] is True


def test_pb_overwrite_replaces_existing_result(tmp_path: Path) -> None:
    out = tmp_path / "out"
    out.mkdir()
    existing = out / "result.json"
    existing.write_text(json.dumps({"status": "success", "old": True}))

    result = CliRunner().invoke(cli, [
        "posttree", "modelcompare", "pb",
        "--sitelogl-dir", str(_make_sitelogl_dir(tmp_path)),
        "--output-dir", str(out), "--overwrite", "--quiet",
    ])
    assert result.exit_code == 0
    assert "old" not in json.loads(existing.read_text())
