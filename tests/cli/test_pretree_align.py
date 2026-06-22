from __future__ import annotations

import json
import shutil
from pathlib import Path

import pytest
from click.testing import CliRunner

from phyloai.cli.main import cli


def _make_seq_dir(tmp_path: Path) -> Path:
    d = tmp_path / "seqs"
    d.mkdir()
    (d / "gene1.fa").write_text(">sp1\nMKTLL\n>sp2\nMKTAA\n>sp3\nMKTVV\n")
    (d / "gene2.fa").write_text(">sp1\nGHTLL\n>sp2\nGHTAA\n>sp3\nGHTVV\n")
    return d


def test_align_help() -> None:
    runner = CliRunner()
    result = runner.invoke(cli, ["pretree", "align", "--help"])
    assert result.exit_code == 0
    assert "--method" in result.output
    assert "--seq-type" in result.output
    assert "--backtrans" in result.output
    assert "--nt-dir" in result.output
    assert "--threads" in result.output
    assert "--mafft-path" in result.output
    assert "--magus-path" in result.output
    assert "--trimal-path" in result.output
    assert "--input-format" not in result.output


def test_align_dry_run_exits_zero(tmp_path: Path) -> None:
    runner = CliRunner()
    seq_dir = _make_seq_dir(tmp_path)
    out_dir = tmp_path / "out"

    result = runner.invoke(cli, [
        "pretree", "align",
        "--seq-dir", str(seq_dir),
        "--output-dir", str(out_dir),
        "--method", "linsi",
        "--seq-type", "AA",
        "--dry-run",
    ])

    assert result.exit_code == 0
    assert not out_dir.exists()
    assert "mafft" in result.output
    assert "gene1.fa" in result.output


def test_align_backtrans_without_nt_dir_exits_1(tmp_path: Path) -> None:
    runner = CliRunner()
    seq_dir = _make_seq_dir(tmp_path)
    out_dir = tmp_path / "out"

    result = runner.invoke(cli, [
        "pretree", "align",
        "--seq-dir", str(seq_dir),
        "--output-dir", str(out_dir),
        "--method", "linsi",
        "--seq-type", "AA",
        "--backtrans",
    ])

    assert result.exit_code == 1


def test_align_nt_seq_type_with_backtrans_exits_1(tmp_path: Path) -> None:
    runner = CliRunner()
    seq_dir = _make_seq_dir(tmp_path)
    nt_dir = tmp_path / "nt"
    nt_dir.mkdir()
    out_dir = tmp_path / "out"

    result = runner.invoke(cli, [
        "pretree", "align",
        "--seq-dir", str(seq_dir),
        "--output-dir", str(out_dir),
        "--method", "linsi",
        "--seq-type", "NT",
        "--backtrans",
        "--nt-dir", str(nt_dir),
    ])

    assert result.exit_code == 1


def test_align_dry_run_ignores_output_dir_conflict(tmp_path: Path) -> None:
    runner = CliRunner()
    seq_dir = _make_seq_dir(tmp_path)
    out_dir = tmp_path / "out"
    out_dir.mkdir()
    (out_dir / "old.txt").write_text("old")

    result = runner.invoke(cli, [
        "pretree", "align",
        "--seq-dir", str(seq_dir),
        "--output-dir", str(out_dir),
        "--method", "linsi",
        "--dry-run",
    ])

    assert result.exit_code == 0


def test_align_writes_result_json(tmp_path: Path) -> None:
    if not shutil.which("mafft"):
        pytest.skip("mafft not found")
    runner = CliRunner()
    seq_dir = _make_seq_dir(tmp_path)
    out_dir = tmp_path / "out"

    result = runner.invoke(cli, [
        "pretree", "align",
        "--seq-dir", str(seq_dir),
        "--output-dir", str(out_dir),
        "--method", "linsi",
        "--seq-type", "AA",
        "--threads", "2",
    ])

    assert result.exit_code == 0
    result_json = out_dir / "result.json"
    assert result_json.exists()
    payload = json.loads(result_json.read_text())
    assert payload["status"] == "success"
    assert payload["key_results"]["n_aligned"] == 2
    assert (out_dir / "seqs" / "gene1.fa").exists()
    assert (out_dir / "seqs" / "gene2.fa").exists()


def test_align_quiet_suppresses_rich_output(tmp_path: Path) -> None:
    if not shutil.which("mafft"):
        pytest.skip("mafft not found")
    runner = CliRunner()
    seq_dir = _make_seq_dir(tmp_path)
    out_dir = tmp_path / "out"

    result = runner.invoke(cli, [
        "pretree", "align",
        "--seq-dir", str(seq_dir),
        "--output-dir", str(out_dir),
        "--method", "linsi",
        "--quiet",
    ])

    assert result.exit_code == 0
    assert "pretree align summary" not in result.output
