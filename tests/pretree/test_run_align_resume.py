from __future__ import annotations

import json
import os
import shutil
import time
from pathlib import Path

import pytest


def test_run_align_resume_requires_checkpoint(tmp_path: Path) -> None:
    from phyloai.pretree.align import run_align

    seq_dir = tmp_path / "seqs"
    seq_dir.mkdir()
    (seq_dir / "gene1.fa").write_text(">a\nMKT\n>b\nMKA\n")
    out_dir = tmp_path / "out"
    out_dir.mkdir()

    with pytest.raises(ValueError, match="No checkpoint"):
        run_align(
            seq_dir=seq_dir,
            output_dir=out_dir,
            method="linsi",
            seq_type="AA",
            resume=True,
        )


def test_run_align_resume_rejects_overwrite(tmp_path: Path) -> None:
    from phyloai.pretree.align import run_align

    seq_dir = tmp_path / "seqs"
    seq_dir.mkdir()
    (seq_dir / "gene1.fa").write_text(">a\nMKT\n>b\nMKA\n")
    out_dir = tmp_path / "out"
    out_dir.mkdir()

    with pytest.raises(ValueError, match="--overwrite and --resume"):
        run_align(
            seq_dir=seq_dir,
            output_dir=out_dir,
            method="linsi",
            seq_type="AA",
            overwrite=True,
            resume=True,
        )


def test_run_align_resume_detects_param_mismatch(tmp_path: Path) -> None:
    from phyloai.pretree.align import _resolved_align_params, run_align
    from phyloai.pretree.checkpoint_helpers import build_initial_checkpoint

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
        method="fftns1",
        resolved_seq_type="AA",
        backtrans=False,
        nt_dir=None,
        threads=2,
        extra_args=None,
        mafft_executable="mafft",
        magus_executable="magus",
        trimal_executable="trimal",
        quiet=False,
    )
    checkpoint = build_initial_checkpoint(
        step="pretree.align",
        command="phyloai pretree align",
        params=params,
        inputs=[seq_dir / "gene1.fa"],
        output_for=lambda p: aa_dir / f"{p.stem}.fa",
        nt_output_for=lambda p: None,
    )
    (out_dir / "checkpoint.json").write_text(json.dumps(checkpoint.to_dict()))

    with pytest.raises(ValueError, match="Resume parameter mismatch"):
        run_align(
            seq_dir=seq_dir,
            output_dir=out_dir,
            method="linsi",
            seq_type="AA",
            threads=2,
            resume=True,
            dry_run=True,
        )


def test_run_align_resume_success_checkpoint_still_validates_params(tmp_path: Path) -> None:
    from phyloai.pretree.align import _resolved_align_params, run_align
    from phyloai.pretree.checkpoint_helpers import build_initial_checkpoint, mark_task

    seq_dir = tmp_path / "seqs"
    seq_dir.mkdir()
    (seq_dir / "gene1.fa").write_text(">a\nMKT\n>b\nMKA\n")
    out_dir = tmp_path / "out"
    out_dir.mkdir()
    aa_dir = out_dir / "seqs"
    aa_dir.mkdir()
    (aa_dir / "gene1.fa").write_text(">a\nMKT\n>b\nMKA\n")

    params = _resolved_align_params(
        seq_dir=seq_dir,
        output_dir=out_dir,
        method="fftns1",
        resolved_seq_type="AA",
        backtrans=False,
        nt_dir=None,
        threads=2,
        extra_args=None,
        mafft_executable="mafft",
        magus_executable="magus",
        trimal_executable="trimal",
        quiet=False,
    )
    checkpoint = build_initial_checkpoint(
        step="pretree.align",
        command="phyloai pretree align",
        params=params,
        inputs=[seq_dir / "gene1.fa"],
        output_for=lambda p: aa_dir / f"{p.stem}.fa",
        nt_output_for=lambda p: None,
    )
    mark_task(checkpoint, "gene1", status="success", reason=None)
    checkpoint.status = "success"
    checkpoint.completed_at = checkpoint.updated_at
    (out_dir / "checkpoint.json").write_text(json.dumps(checkpoint.to_dict()))

    with pytest.raises(ValueError, match="Resume parameter mismatch"):
        run_align(
            seq_dir=seq_dir,
            output_dir=out_dir,
            method="linsi",
            seq_type="AA",
            threads=2,
            resume=True,
            dry_run=True,
        )


def test_run_align_resume_skips_successful_tasks(tmp_path: Path) -> None:
    if not shutil.which("mafft"):
        pytest.skip("mafft not found")

    from phyloai.pretree.align import _align_one, _resolved_align_params, run_align
    from phyloai.pretree.checkpoint_helpers import build_initial_checkpoint, mark_task

    seq_dir = tmp_path / "seqs"
    seq_dir.mkdir()
    (seq_dir / "gene1.fa").write_text(">a\nMKTLLLTLVVVTIVC\n>b\nMKTLLLTLAAVTIVC\n>c\nMKTLLLTLVVVTIVC\n")
    (seq_dir / "gene2.fa").write_text(">a\nGHTLLLTLVVVTIVC\n>b\nGHTLLLTLAAVTIVC\n>c\nGHTLLLTLVVVTIVC\n")
    out_dir = tmp_path / "out"
    out_dir.mkdir()
    aa_dir = out_dir / "seqs"
    aa_dir.mkdir()

    from phyloai.pretree.align import _resolve_tool_paths

    mafft_exe, magus_exe, trimal_exe = _resolve_tool_paths(
        method="linsi",
        backtrans=False,
        mafft_path=None,
        magus_path=None,
        trimal_path=None,
        dry_run=False,
    )

    params = _resolved_align_params(
        seq_dir=seq_dir,
        output_dir=out_dir,
        method="linsi",
        resolved_seq_type="AA",
        backtrans=False,
        nt_dir=None,
        threads=2,
        extra_args=None,
        mafft_executable=mafft_exe,
        magus_executable=magus_exe,
        trimal_executable=trimal_exe,
        quiet=False,
    )
    checkpoint = build_initial_checkpoint(
        step="pretree.align",
        command="phyloai pretree align",
        params=params,
        inputs=[seq_dir / "gene1.fa", seq_dir / "gene2.fa"],
        output_for=lambda p: aa_dir / f"{p.stem}.fa",
        nt_output_for=lambda p: None,
    )

    first = _align_one(
        seq_dir / "gene1.fa",
        aa_dir,
        method="linsi",
        seq_type="AA",
        extra_args=None,
        dry_run=False,
    )
    assert first["status"] == "success"
    mark_task(checkpoint, "gene1", status="success", reason=None)
    mark_task(checkpoint, "gene2", status="failed", reason="previous run error")
    (out_dir / "checkpoint.json").write_text(json.dumps(checkpoint.to_dict()))

    payload = run_align(
        seq_dir=seq_dir,
        output_dir=out_dir,
        method="linsi",
        seq_type="AA",
        threads=2,
        resume=True,
    )

    assert payload["status"] == "success"
    assert payload["key_results"]["n_aligned"] == 2
    assert (aa_dir / "gene1.fa").exists()
    assert (aa_dir / "gene2.fa").exists()


def test_run_align_backtrans_resume_skips_completed_genes(tmp_path: Path) -> None:
    """Regression: with --backtrans, an interrupted run must not re-align genes
    whose AA+NT outputs already exist. Previously the checkpoint marked AA tasks
    success before backtrans ran, so resume saw missing NT files and reran all."""
    if not shutil.which("mafft") or not shutil.which("trimal"):
        pytest.skip("mafft and trimal required")

    from phyloai.core.checkpoint import load_checkpoint
    from phyloai.pretree.align import run_align

    faa = tmp_path / "faa"
    fna = tmp_path / "fna"
    faa.mkdir()
    fna.mkdir()

    genes = {
        "g1": ("MKTLLLT", "ATGAAAACTTTGCTTTTGACT"),
        "g2": ("GHTLLLT", "GGTCATACTTTGCTTTTGACT"),
        "g3": ("MKAVVVT", "ATGAAAGCTGTTGTTGTTACT"),
        "g4": ("WYTLLLT", "TGGTATACTTTGCTTTTGACT"),
    }
    for name, (aa, nt) in genes.items():
        (faa / f"{name}.fa").write_text(f">a\n{aa}\n>b\n{aa}\n>c\n{aa}\n")
        (fna / f"{name}.fa").write_text(f">a\n{nt}\n>b\n{nt}\n>c\n{nt}\n")

    out_dir = tmp_path / "out"

    # Full run: all genes get AA + NT outputs and a success checkpoint.
    run_align(
        seq_dir=faa, output_dir=out_dir, method="fftns1", seq_type="AA",
        backtrans=True, nt_dir=fna, threads=2,
    )

    ckpt = load_checkpoint(out_dir / "checkpoint.json")
    assert all(t.status == "success" for t in ckpt.tasks)
    for t in ckpt.tasks:
        assert t.outputs.get("nt") and Path(t.outputs["nt"]).exists()

    # Simulate an interrupted run: g3, g4 incomplete (status running, outputs gone).
    for t in ckpt.tasks:
        if t.task_id in {"g3", "g4"}:
            t.status = "running"
            Path(t.outputs["aa"]).unlink(missing_ok=True)
            Path(t.outputs["nt"]).unlink(missing_ok=True)
    ckpt.status = "interrupted"
    (out_dir / "checkpoint.json").write_text(json.dumps(ckpt.to_dict()))

    aa_dir = out_dir / "seqs" / "faa"
    mtimes_before = {g: os.path.getmtime(aa_dir / f"{g}.fa") for g in ("g1", "g2")}
    time.sleep(0.05)

    payload = run_align(
        seq_dir=faa, output_dir=out_dir, method="fftns1", seq_type="AA",
        backtrans=True, nt_dir=fna, threads=2, resume=True,
    )

    assert payload["status"] == "success"
    # Completed genes must NOT be re-aligned (mtime unchanged).
    for g in ("g1", "g2"):
        assert os.path.getmtime(aa_dir / f"{g}.fa") == mtimes_before[g], (
            f"{g} was re-aligned on resume but should have been skipped"
        )
    # All genes complete with NT outputs after resume.
    ckpt2 = load_checkpoint(out_dir / "checkpoint.json")
    assert all(t.status == "success" for t in ckpt2.tasks)
    for t in ckpt2.tasks:
        assert Path(t.outputs["nt"]).exists()
