from __future__ import annotations

import json
from pathlib import Path

from click.testing import CliRunner

from phyloai.cli.main import cli


def test_align_cli_rejects_overwrite_with_resume(tmp_path: Path) -> None:
    runner = CliRunner()
    seq_dir = tmp_path / "seqs"
    seq_dir.mkdir()
    (seq_dir / "gene1.fa").write_text(">a\nMKT\n>b\nMKA\n")
    out_dir = tmp_path / "out"
    out_dir.mkdir()
    (out_dir / "old.txt").write_text("old")

    result = runner.invoke(
        cli,
        [
            "pretree", "align",
            "--seq-dir", str(seq_dir),
            "--method", "linsi",
            "--seq-type", "AA",
            "--output-dir", str(out_dir),
            "--resume",
            "--overwrite",
        ],
    )

    assert result.exit_code == 1
    assert "--overwrite and --resume" in result.output


def test_align_cli_resume_requires_checkpoint(tmp_path: Path) -> None:
    runner = CliRunner()
    seq_dir = tmp_path / "seqs"
    seq_dir.mkdir()
    (seq_dir / "gene1.fa").write_text(">a\nMKT\n>b\nMKA\n")
    out_dir = tmp_path / "out"
    out_dir.mkdir()

    result = runner.invoke(
        cli,
        [
            "pretree", "align",
            "--seq-dir", str(seq_dir),
            "--method", "linsi",
            "--seq-type", "AA",
            "--output-dir", str(out_dir),
            "--resume",
        ],
    )

    assert result.exit_code == 1
    assert "No checkpoint" in result.output


def test_align_cli_resume_dry_run_validates_params(tmp_path: Path) -> None:
    from phyloai.pretree.align import _resolved_align_params
    from phyloai.pretree.checkpoint_helpers import build_initial_checkpoint

    runner = CliRunner()
    seq_dir = tmp_path / "seqs"
    seq_dir.mkdir()
    (seq_dir / "gene1.fa").write_text(">a\nMKT\n>b\nMKA\n")
    out_dir = tmp_path / "out"
    out_dir.mkdir()
    aa_dir = out_dir / "seqs"
    aa_dir.mkdir()

    params = _resolved_align_params(
        seq_dir=seq_dir,
        output_dir=out_dir,
        method="linsi",
        resolved_seq_type="AA",
        backtrans=False,
        nt_dir=None,
        threads=4,
        tool_args=None,
        mafft_path="mafft",
        magus_path="magus",
        trimal_path="trimal",
        quiet=False,
    )
    checkpoint = build_initial_checkpoint(
        step="pretree.align",
        command="phyloai pretree align --seq-dir /data --output-dir /out --method linsi --seq-type AA --threads 4",
        params=params,
        inputs=[seq_dir / "gene1.fa"],
        output_for=lambda p: aa_dir / f"{p.stem}.fa",
        nt_output_for=lambda p: None,
    )
    (out_dir / "checkpoint.json").write_text(json.dumps(checkpoint.to_dict()))

    result = runner.invoke(
        cli,
        [
            "pretree", "align",
            "--seq-dir", str(seq_dir),
            "--method", "fftns1",
            "--seq-type", "AA",
            "--output-dir", str(out_dir),
            "--resume",
            "--dry-run",
        ],
    )

    assert result.exit_code == 1
    assert "Resume parameter mismatch" in result.output


def test_align_cli_resume_dry_run_reports_summary(tmp_path: Path) -> None:
    from phyloai.pretree.align import _resolved_align_params
    from phyloai.pretree.checkpoint_helpers import build_initial_checkpoint, mark_task

    runner = CliRunner()
    seq_dir = tmp_path / "seqs"
    seq_dir.mkdir()
    (seq_dir / "gene1.fa").write_text(">a\nMKT\n>b\nMKA\n")
    (seq_dir / "gene2.fa").write_text(">a\nGHT\n>b\nGHA\n")
    out_dir = tmp_path / "out"
    out_dir.mkdir()
    aa_dir = out_dir / "seqs"
    aa_dir.mkdir()
    (aa_dir / "gene1.fa").write_text(">a\nMKT\n>b\nMKA\n")

    params = _resolved_align_params(
        seq_dir=seq_dir,
        output_dir=out_dir,
        method="linsi",
        resolved_seq_type="AA",
        backtrans=False,
        nt_dir=None,
        threads=4,
        tool_args=None,
        mafft_path="mafft",
        magus_path="magus",
        trimal_path="trimal",
        quiet=False,
    )
    checkpoint = build_initial_checkpoint(
        step="pretree.align",
        command="phyloai pretree align --seq-dir /data --output-dir /out --method linsi --seq-type AA --threads 4",
        params=params,
        inputs=[seq_dir / "gene1.fa", seq_dir / "gene2.fa"],
        output_for=lambda p: aa_dir / f"{p.stem}.fa",
        nt_output_for=lambda p: None,
    )
    mark_task(checkpoint, "gene1", status="success", reason=None)
    mark_task(checkpoint, "gene2", status="failed", reason="previous run error")
    (out_dir / "checkpoint.json").write_text(json.dumps(checkpoint.to_dict()))

    result = runner.invoke(
        cli,
        [
            "pretree", "align",
            "--seq-dir", str(seq_dir),
            "--method", "linsi",
            "--seq-type", "AA",
            "--output-dir", str(out_dir),
            "--resume",
            "--dry-run",
        ],
    )

    assert result.exit_code == 0
    assert "Resume dry-run:" in result.output
    assert "skip 1 tasks" in result.output
    assert "rerun 1 tasks" in result.output
