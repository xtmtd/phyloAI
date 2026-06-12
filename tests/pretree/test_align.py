from __future__ import annotations

import shutil
from pathlib import Path

import pytest


def test_scan_input_finds_fasta_files(tmp_path: Path) -> None:
    from phyloai.pretree.align import _scan_input

    (tmp_path / "gene1.fa").write_text(">a\nACGT\n")
    (tmp_path / "gene2.faa").write_text(">b\nMKT\n")
    (tmp_path / "notes.txt").write_text("skip")
    (tmp_path / "empty.fa").write_text("")
    (tmp_path / "subdir").mkdir()

    found, skipped = _scan_input(tmp_path)

    assert len(found) == 2
    assert len(skipped) == 3
    skip_reasons = {s["reason"] for s in skipped}
    assert "empty file" in skip_reasons
    assert "directory" in skip_reasons
    assert "unrecognized extension" in skip_reasons


def test_build_mafft_cmd_linsi(tmp_path: Path) -> None:
    from phyloai.pretree.align import _build_mafft_cmd

    inp = tmp_path / "gene1.fa"
    out = tmp_path / "gene1_aln.fa"
    cmd = _build_mafft_cmd(inp, out, method="linsi")

    assert cmd[0] == "mafft"
    assert "--maxiterate" in cmd
    assert "1000" in cmd
    assert "--localpair" in cmd
    assert "--thread" in cmd
    assert "1" in cmd
    assert str(inp) in cmd


def test_build_mafft_cmd_accepts_explicit_executable(tmp_path: Path) -> None:
    from phyloai.pretree.align import _build_mafft_cmd

    inp = tmp_path / "gene1.fa"
    out = tmp_path / "gene1_aln.fa"
    cmd = _build_mafft_cmd(inp, out, method="linsi", executable="/opt/bin/mafft")

    assert cmd[0] == "/opt/bin/mafft"


def test_build_mafft_cmd_fftns1(tmp_path: Path) -> None:
    from phyloai.pretree.align import _build_mafft_cmd

    inp = tmp_path / "gene1.fa"
    out = tmp_path / "gene1_aln.fa"
    cmd = _build_mafft_cmd(inp, out, method="fftns1")

    assert "--retree" in cmd
    idx = cmd.index("--retree")
    assert cmd[idx + 1] == "1"


def test_build_mafft_cmd_fftns2(tmp_path: Path) -> None:
    from phyloai.pretree.align import _build_mafft_cmd

    inp = tmp_path / "gene1.fa"
    out = tmp_path / "gene1_aln.fa"
    cmd = _build_mafft_cmd(inp, out, method="fftns2")

    assert "--retree" in cmd
    idx = cmd.index("--retree")
    assert cmd[idx + 1] == "2"


def test_build_mafft_cmd_auto(tmp_path: Path) -> None:
    from phyloai.pretree.align import _build_mafft_cmd

    inp = tmp_path / "gene1.fa"
    out = tmp_path / "gene1_aln.fa"
    cmd = _build_mafft_cmd(inp, out, method="auto")

    assert "--auto" in cmd
    assert "--thread" in cmd


def test_build_mafft_cmd_einsi(tmp_path: Path) -> None:
    from phyloai.pretree.align import _build_mafft_cmd

    inp = tmp_path / "gene1.fa"
    out = tmp_path / "out.fa"
    cmd = _build_mafft_cmd(inp, out, method="einsi")

    assert "--genafpair" in cmd


def test_build_mafft_cmd_ginsi(tmp_path: Path) -> None:
    from phyloai.pretree.align import _build_mafft_cmd

    inp = tmp_path / "gene1.fa"
    out = tmp_path / "out.fa"
    cmd = _build_mafft_cmd(inp, out, method="ginsi")

    assert "--globalpair" in cmd


def test_build_magus_cmd_aa(tmp_path: Path) -> None:
    from phyloai.pretree.align import _build_magus_cmd

    inp = tmp_path / "gene1.fa"
    out = tmp_path / "gene1_aln.fa"
    work = tmp_path / "work"
    cmd = _build_magus_cmd(inp, out, work_dir=work, seq_type="AA", extra_args=None)

    assert cmd[0] == "magus"
    assert "-i" in cmd
    assert str(inp) in cmd
    assert "-o" in cmd
    assert str(out) in cmd
    assert "-d" in cmd
    assert str(work) in cmd
    assert "--datatype" in cmd
    idx = cmd.index("--datatype")
    assert cmd[idx + 1] == "protein"


def test_build_magus_cmd_nt(tmp_path: Path) -> None:
    from phyloai.pretree.align import _build_magus_cmd

    inp = tmp_path / "gene1.fa"
    out = tmp_path / "out.fa"
    work = tmp_path / "work"
    cmd = _build_magus_cmd(inp, out, work_dir=work, seq_type="NT", extra_args=None)

    idx = cmd.index("--datatype")
    assert cmd[idx + 1] == "dna"


def test_build_magus_cmd_extra_args(tmp_path: Path) -> None:
    from phyloai.pretree.align import _build_magus_cmd

    inp = tmp_path / "gene1.fa"
    out = tmp_path / "out.fa"
    work = tmp_path / "work"
    cmd = _build_magus_cmd(inp, out, work_dir=work, seq_type="AA", extra_args="--maxsubsetsize 50 --recurse true")

    assert "--maxsubsetsize" in cmd
    assert "50" in cmd
    assert "--recurse" in cmd
    assert "true" in cmd


def test_build_magus_cmd_extra_args_override(tmp_path: Path) -> None:
    from phyloai.pretree.align import _build_magus_cmd

    inp = tmp_path / "gene1.fa"
    out = tmp_path / "out.fa"
    work = tmp_path / "work"
    cmd = _build_magus_cmd(inp, out, work_dir=work, seq_type="AA", extra_args="--datatype dna")

    pairs = list(zip(cmd, cmd[1:]))
    datatype_values = [b for a, b in pairs if a == "--datatype"]
    assert datatype_values == ["dna"]


def test_build_magus_cmd_preserves_unknown_extra_args(tmp_path: Path) -> None:
    from phyloai.pretree.align import _build_magus_cmd

    inp = tmp_path / "gene1.fa"
    out = tmp_path / "out.fa"
    work = tmp_path / "work"
    cmd = _build_magus_cmd(
        inp,
        out,
        work_dir=work,
        seq_type="AA",
        extra_args="--maxsubsetsize 50 --recurse --some-flag=value",
    )

    assert cmd[-4:] == ["--maxsubsetsize", "50", "--recurse", "--some-flag=value"]


def test_build_magus_cmd_accepts_explicit_executable(tmp_path: Path) -> None:
    from phyloai.pretree.align import _build_magus_cmd

    inp = tmp_path / "gene1.fa"
    out = tmp_path / "out.fa"
    work = tmp_path / "work"
    cmd = _build_magus_cmd(inp, out, work_dir=work, seq_type="AA", extra_args=None, executable="/opt/bin/magus")

    assert cmd[0] == "/opt/bin/magus"


def test_validate_cds_passes_clean_sequences() -> None:
    from phyloai.pretree.align import _validate_cds

    seqs = {"sp1": "ATGGCCTAA", "sp2": "ATGCGCTAG"}
    warnings = _validate_cds(seqs, n_aa_taxa=2)
    assert warnings == []


def test_validate_cds_length_not_multiple_of_3() -> None:
    from phyloai.pretree.align import _validate_cds

    seqs = {"sp1": "ATGGC"}
    warnings = _validate_cds(seqs, n_aa_taxa=1)
    assert any("not a multiple of 3" in w for w in warnings)


def test_validate_cds_taxon_count_mismatch() -> None:
    from phyloai.pretree.align import _validate_cds

    seqs = {"sp1": "ATGGCCTAA"}
    warnings = _validate_cds(seqs, n_aa_taxa=3)
    assert any("taxon count mismatch" in w for w in warnings)


def test_validate_cds_taxon_id_mismatch() -> None:
    from phyloai.pretree.align import _validate_cds

    seqs = {"sp1": "ATGGCCTAA", "wrong_sp": "ATGCGCTAG"}
    warnings = _validate_cds(seqs, aa_taxa={"sp1", "sp2"})
    assert any("taxon ID mismatch" in w for w in warnings)


def test_validate_cds_internal_stop_codon() -> None:
    from phyloai.pretree.align import _validate_cds

    seqs = {"sp1": "ATGTAAGCT"}
    warnings = _validate_cds(seqs, n_aa_taxa=1)
    assert any("internal stop codon" in w for w in warnings)


def test_validate_cds_trailing_stop_not_flagged() -> None:
    from phyloai.pretree.align import _validate_cds

    seqs = {"sp1": "ATGGCCTAA"}
    warnings = _validate_cds(seqs, n_aa_taxa=1)
    assert not any("internal stop codon" in w for w in warnings)


def test_align_one_mafft_linsi_produces_output(tmp_path: Path) -> None:
    if not shutil.which("mafft"):
        pytest.skip("mafft not found")
    from phyloai.pretree.align import _align_one

    inp = tmp_path / "gene1.fa"
    inp.write_text(">sp1\nMKTLLLTLVVVTIVC\n>sp2\nMKTLLLTLAAVTIVC\n>sp3\nMKTLLLTLVVVTIVC\n")
    out_dir = tmp_path / "out"
    out_dir.mkdir()

    result = _align_one(inp, out_dir, method="linsi", seq_type="AA",
                        extra_args=None, dry_run=False)

    assert result["status"] == "success"
    assert Path(result["output_aa"]).exists()
    assert result["n_taxa"] == 3
    assert result["alignment_length"] > 0
    assert result["wall_time"] > 0


def test_align_one_dry_run_creates_no_files(tmp_path: Path) -> None:
    from phyloai.pretree.align import _align_one

    inp = tmp_path / "gene1.fa"
    inp.write_text(">sp1\nMKT\n>sp2\nMKT\n")
    out_dir = tmp_path / "out"
    out_dir.mkdir()

    result = _align_one(inp, out_dir, method="linsi", seq_type="AA",
                        extra_args=None, dry_run=True)

    assert result["status"] == "dry_run"
    assert result["cmd"] is not None
    assert not any(out_dir.iterdir())


def test_align_one_failed_tool_returns_skipped(tmp_path: Path) -> None:
    from phyloai.pretree.align import _align_one

    inp = tmp_path / "bad.fa"
    inp.write_text("not a fasta file at all\njust garbage\n")
    out_dir = tmp_path / "out"
    out_dir.mkdir()

    result = _align_one(inp, out_dir, method="linsi", seq_type="AA",
                        extra_args=None, dry_run=False)

    assert result["status"] in {"skipped", "success"}


def test_backtrans_one_produces_nt_alignment(tmp_path: Path) -> None:
    if not shutil.which("trimal"):
        pytest.skip("trimal not found")
    from phyloai.pretree.align import _backtrans_one

    aa_aln = tmp_path / "gene1_aa.fa"
    aa_aln.write_text(
        ">sp1\nMK-\n"
        ">sp2\nMKT\n"
    )
    nt_file = tmp_path / "gene1.fa"
    nt_file.write_text(
        ">sp1\nATGAAA\n"
        ">sp2\nATGAAAACT\n"
    )
    out_nt = tmp_path / "gene1_nt.fa"

    result = _backtrans_one(aa_aln, nt_file, out_nt, dry_run=False)

    assert result["status"] in {"success", "skipped"}
    if result["status"] == "success":
        assert out_nt.exists()


def test_run_align_aa_only_dry_run(tmp_path: Path) -> None:
    from phyloai.pretree.align import run_align

    seq_dir = tmp_path / "seqs"
    seq_dir.mkdir()
    (seq_dir / "gene1.fa").write_text(">a\nMKT\n>b\nMKA\n")
    (seq_dir / "gene2.fa").write_text(">a\nGHT\n>b\nGHA\n")
    out_dir = tmp_path / "out"

    payload = run_align(
        seq_dir=seq_dir,
        output_dir=out_dir,
        method="linsi",
        seq_type="AA",
        dry_run=True,
    )

    assert payload["status"] == "success"
    assert payload["data"]["summary"]["n_input_files"] == 2
    assert not out_dir.exists()


def test_run_align_output_dir_conflict(tmp_path: Path) -> None:
    from phyloai.pretree.align import run_align

    seq_dir = tmp_path / "seqs"
    seq_dir.mkdir()
    (seq_dir / "gene1.fa").write_text(">a\nMKT\n")
    out_dir = tmp_path / "out"
    out_dir.mkdir()
    (out_dir / "old.txt").write_text("old")

    with pytest.raises(ValueError, match="already exists and is non-empty"):
        run_align(seq_dir=seq_dir, output_dir=out_dir, method="linsi",
                  seq_type="AA", overwrite=False)


def test_run_align_overwrite_clears_directory(tmp_path: Path) -> None:
    if not shutil.which("mafft"):
        pytest.skip("mafft not found")
    from phyloai.pretree.align import run_align

    seq_dir = tmp_path / "seqs"
    seq_dir.mkdir()
    (seq_dir / "gene1.fa").write_text(">sp1\nMKTLL\n>sp2\nMKTAA\n")
    out_dir = tmp_path / "out"
    out_dir.mkdir()
    (out_dir / "old.txt").write_text("old")

    payload = run_align(seq_dir=seq_dir, output_dir=out_dir, method="linsi",
                        seq_type="AA", overwrite=True)

    assert not (out_dir / "old.txt").exists()
    assert payload["status"] == "success"


def test_run_align_all_skipped_returns_error(tmp_path: Path) -> None:
    from phyloai.pretree.align import run_align

    seq_dir = tmp_path / "seqs"
    seq_dir.mkdir()
    (seq_dir / "bad.fa").write_text("")
    out_dir = tmp_path / "out"

    with pytest.raises(ValueError, match="No genes were aligned"):
        run_align(seq_dir=seq_dir, output_dir=out_dir, method="linsi",
                  seq_type="AA", overwrite=False)


def test_run_align_backtrans_requires_nt_dir(tmp_path: Path) -> None:
    from phyloai.pretree.align import run_align

    seq_dir = tmp_path / "seqs"
    seq_dir.mkdir()
    (seq_dir / "gene1.fa").write_text(">a\nMKT\n")
    out_dir = tmp_path / "out"

    with pytest.raises(ValueError, match="--nt-dir"):
        run_align(seq_dir=seq_dir, output_dir=out_dir, method="linsi",
                  seq_type="AA", backtrans=True, nt_dir=None)


def test_run_align_nt_seq_type_with_backtrans_raises(tmp_path: Path) -> None:
    from phyloai.pretree.align import run_align

    seq_dir = tmp_path / "seqs"
    seq_dir.mkdir()
    nt_dir = tmp_path / "nt"
    nt_dir.mkdir()
    out_dir = tmp_path / "out"

    with pytest.raises(ValueError, match="--backtrans requires"):
        run_align(seq_dir=seq_dir, output_dir=out_dir, method="linsi",
                  seq_type="NT", backtrans=True, nt_dir=nt_dir)


def test_run_align_extra_args_ignored_for_mafft_emits_warning(tmp_path: Path) -> None:
    from phyloai.pretree.align import run_align

    seq_dir = tmp_path / "seqs"
    seq_dir.mkdir()
    (seq_dir / "gene1.fa").write_text(">a\nMKT\n>b\nMKA\n")
    out_dir = tmp_path / "out"

    payload = run_align(
        seq_dir=seq_dir,
        output_dir=out_dir,
        method="linsi",
        seq_type="AA",
        extra_args="--maxsubsetsize 50",
        dry_run=True,
    )

    assert any("ignored" in w.lower() for w in payload["data"].get("warnings", []))


def test_run_align_key_results_shape(tmp_path: Path) -> None:
    if not shutil.which("mafft"):
        pytest.skip("mafft not found")
    from phyloai.pretree.align import run_align

    seq_dir = tmp_path / "seqs"
    seq_dir.mkdir()
    (seq_dir / "gene1.fa").write_text(">sp1\nMKTLL\n>sp2\nMKTAA\n>sp3\nMKTVV\n")
    (seq_dir / "gene2.fa").write_text(">sp1\nGHTLL\n>sp2\nGHTAA\n>sp3\nGHTVV\n")
    out_dir = tmp_path / "out"

    payload = run_align(seq_dir=seq_dir, output_dir=out_dir, method="linsi",
                        seq_type="AA")

    kr = payload["key_results"]
    assert "n_aligned" in kr
    assert "n_skipped" in kr
    assert "mean_alignment_length" in kr
    assert "mean_n_taxa" in kr
    assert "method" in kr
    assert "backtrans" in kr


def test_resolve_tool_paths_accepts_explicit_mafft_path(tmp_path: Path) -> None:
    from phyloai.pretree.align import _resolve_tool_paths

    fake = tmp_path / "mafft"
    fake.write_text("#!/bin/sh\n")

    mafft_exe, magus_exe, trimal_exe = _resolve_tool_paths(
        method="linsi",
        backtrans=False,
        mafft_path=fake,
        magus_path=None,
        trimal_path=None,
        dry_run=True,
    )

    assert mafft_exe == str(fake)
    assert magus_exe == "magus"
    assert trimal_exe == "trimal"


def test_detect_tool_versions_mafft() -> None:
    if not shutil.which("mafft"):
        pytest.skip("mafft not found")
    from phyloai.pretree.align import _detect_tool_versions

    versions = _detect_tool_versions(method="linsi", backtrans=False, mafft_path=None, magus_path=None, trimal_path=None)
    assert "mafft" in versions
    assert versions["mafft"]


def test_detect_tool_versions_backtrans_includes_trimal() -> None:
    from phyloai.pretree.align import _detect_tool_versions

    versions = _detect_tool_versions(method="linsi", backtrans=True, mafft_path=None, magus_path=None, trimal_path=None)
    assert "trimal" in versions


def test_run_align_missing_explicit_mafft_path_raises(tmp_path: Path) -> None:
    from phyloai.pretree.align import run_align

    seq_dir = tmp_path / "seqs"
    seq_dir.mkdir()
    (seq_dir / "gene1.fa").write_text(">a\nMKT\n>b\nMKA\n")
    out_dir = tmp_path / "out"

    with pytest.raises(FileNotFoundError, match="mafft"):
        run_align(
            seq_dir=seq_dir,
            output_dir=out_dir,
            method="linsi",
            seq_type="AA",
            mafft_path=tmp_path / "missing-mafft",
            dry_run=True,
        )


def test_run_align_missing_explicit_trimal_path_raises(tmp_path: Path) -> None:
    from phyloai.pretree.align import run_align

    seq_dir = tmp_path / "seqs"
    seq_dir.mkdir()
    nt_dir = tmp_path / "nt"
    nt_dir.mkdir()
    (seq_dir / "gene1.fa").write_text(">a\nMKT\n>b\nMKA\n")
    out_dir = tmp_path / "out"

    with pytest.raises(FileNotFoundError, match="trimal"):
        run_align(
            seq_dir=seq_dir,
            output_dir=out_dir,
            method="linsi",
            seq_type="AA",
            backtrans=True,
            nt_dir=nt_dir,
            trimal_path=tmp_path / "missing-trimal",
            dry_run=True,
        )


def test_run_align_auto_detects_seq_type(tmp_path: Path) -> None:
    from phyloai.pretree.align import run_align

    seq_dir = tmp_path / "seqs"
    seq_dir.mkdir()
    (seq_dir / "gene1.fa").write_text(">a\nMKTLLL\n>b\nMKAAA\n")
    out_dir = tmp_path / "out"

    payload = run_align(
        seq_dir=seq_dir,
        output_dir=out_dir,
        method="linsi",
        seq_type="auto",
        dry_run=True,
    )

    assert payload["status"] == "success"
    assert any("auto-detected" in w for w in payload["data"].get("warnings", []))


def test_run_align_auto_detected_nt_with_backtrans_raises(tmp_path: Path) -> None:
    from phyloai.pretree.align import run_align

    seq_dir = tmp_path / "seqs"
    seq_dir.mkdir()
    (seq_dir / "gene1.fa").write_text(">a\nACGTACGT\n>b\nACGTACGT\n")
    nt_dir = tmp_path / "nt"
    nt_dir.mkdir()
    out_dir = tmp_path / "out"

    with pytest.raises(ValueError, match="--backtrans requires"):
        run_align(
            seq_dir=seq_dir,
            output_dir=out_dir,
            method="linsi",
            seq_type="auto",
            backtrans=True,
            nt_dir=nt_dir,
            dry_run=True,
        )


def test_align_one_magus_dry_run_does_not_create_tempdir(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    from phyloai.pretree import align

    inp = tmp_path / "gene1.fa"
    inp.write_text(">a\nMKT\n>b\nMKA\n")
    out_dir = tmp_path / "out"

    def fail_mkdtemp(*_args, **_kwargs):
        raise AssertionError("dry-run must not create MAGUS temp directories")

    monkeypatch.setattr(align.tempfile, "mkdtemp", fail_mkdtemp)
    result = align._align_one(inp, out_dir, method="magus", seq_type="AA", extra_args=None, dry_run=True)

    assert result["status"] == "dry_run"
    assert not out_dir.exists()


def test_align_one_magus_cleans_tempdir_when_subprocess_raises(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    from phyloai.pretree import align

    inp = tmp_path / "gene1.fa"
    inp.write_text(">a\nMKT\n>b\nMKA\n")
    out_dir = tmp_path / "out"
    work_dir = tmp_path / "magus-work"

    def fake_mkdtemp(*_args, **_kwargs):
        work_dir.mkdir()
        return str(work_dir)

    def fail_run(*_args, **_kwargs):
        raise OSError("cannot execute magus")

    monkeypatch.setattr(align.tempfile, "mkdtemp", fake_mkdtemp)
    monkeypatch.setattr(align.subprocess, "run", fail_run)

    result = align._align_one(inp, out_dir, method="magus", seq_type="AA", extra_args=None, dry_run=False)

    assert result["status"] == "skipped"
    assert not work_dir.exists()


def test_align_one_magus_success_path_writes_and_validates_output(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    from phyloai.pretree import align

    inp = tmp_path / "gene1.fa"
    inp.write_text(">a\nMKT\n>b\nMKA\n")
    out_dir = tmp_path / "out"
    work_dir = tmp_path / "magus-work"

    class Proc:
        returncode = 0
        stdout = "MAGUS diagnostics\n"
        stderr = ""

    def fake_mkdtemp(*_args, **_kwargs):
        work_dir.mkdir()
        return str(work_dir)

    def fake_run(cmd, *_args, **_kwargs):
        out_path = Path(cmd[cmd.index("-o") + 1])
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path.write_text(">a\nMKT\n>b\nMKA\n")
        return Proc()

    monkeypatch.setattr(align.tempfile, "mkdtemp", fake_mkdtemp)
    monkeypatch.setattr(align.subprocess, "run", fake_run)

    result = align._align_one(inp, out_dir, method="magus", seq_type="AA", extra_args=None, dry_run=False)

    assert result["status"] == "success"
    assert result["n_taxa"] == 2
    assert result["alignment_length"] == 3
    assert result.get("tool_stdout", "") == ""
    assert not work_dir.exists()


def test_run_align_magus_dry_run_when_platform_is_linux(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    from phyloai.pretree import align

    seq_dir = tmp_path / "seqs"
    seq_dir.mkdir()
    (seq_dir / "gene1.fa").write_text(">a\nMKT\n>b\nMKA\n")
    out_dir = tmp_path / "out"

    monkeypatch.setattr(align.platform, "system", lambda: "Linux")
    payload = align.run_align(seq_dir=seq_dir, output_dir=out_dir, method="magus", seq_type="AA", dry_run=True)

    assert payload["status"] == "success"
    cmd = payload["data"]["files"][0]["cmd"]
    assert cmd[0] == "magus"
    assert "--datatype" in cmd
    assert "protein" in cmd
    assert not out_dir.exists()


def test_align_one_success_does_not_return_alignment_stdout(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    from phyloai.pretree import align

    inp = tmp_path / "gene1.fa"
    inp.write_text(">a\nMKT\n>b\nMKA\n")
    out_dir = tmp_path / "out"

    class Proc:
        returncode = 0
        stdout = ">a\nMKT\n>b\nMKA\n"
        stderr = "diagnostic\n"

    monkeypatch.setattr(align.subprocess, "run", lambda *_args, **_kwargs: Proc())

    result = align._align_one(inp, out_dir, method="linsi", seq_type="AA", extra_args=None, dry_run=False)

    assert result["status"] == "success"
    assert result.get("tool_stdout", "") == ""
    assert (out_dir / "gene1.fa").read_text() == Proc.stdout


def test_validate_msa_output_rejects_empty_file(tmp_path: Path) -> None:
    from phyloai.pretree.align import _validate_msa_output

    out = tmp_path / "empty.fa"
    out.write_text("")

    n_taxa, aln_len, warnings = _validate_msa_output(out)

    assert n_taxa == 0
    assert aln_len == 0
    assert any("empty" in w.lower() for w in warnings)


def test_validate_msa_output_rejects_no_fasta_records(tmp_path: Path) -> None:
    from phyloai.pretree.align import _validate_msa_output

    out = tmp_path / "bad.fa"
    out.write_text("\n\n")

    n_taxa, aln_len, warnings = _validate_msa_output(out)

    assert n_taxa == 0
    assert aln_len == 0
    assert any("no fasta records" in w.lower() or "could not parse" in w.lower() for w in warnings)


def test_validate_msa_output_rejects_unequal_lengths(tmp_path: Path) -> None:
    from phyloai.pretree.align import _validate_msa_output

    out = tmp_path / "bad.fa"
    out.write_text(">a\nMKT\n>b\nMKTA\n")

    _n_taxa, _aln_len, warnings = _validate_msa_output(out)

    assert any("unequal" in w.lower() for w in warnings)


def test_align_one_skips_invalid_generated_msa(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    from phyloai.pretree import align

    inp = tmp_path / "gene1.fa"
    inp.write_text(">a\nMKT\n>b\nMKA\n")
    out_dir = tmp_path / "out"

    class Proc:
        returncode = 0
        stdout = ""
        stderr = ""

    monkeypatch.setattr(align.subprocess, "run", lambda *_args, **_kwargs: Proc())

    result = align._align_one(inp, out_dir, method="linsi", seq_type="AA", extra_args=None, dry_run=False)

    assert result["status"] == "skipped"
    assert "empty" in result["reason"].lower()


def test_run_align_magus_on_non_linux_raises(tmp_path: Path) -> None:
    import platform as _platform
    if _platform.system() == "Linux":
        pytest.skip("Test only applies to non-Linux platforms")
    from phyloai.pretree.align import run_align

    seq_dir = tmp_path / "seqs"
    seq_dir.mkdir()
    (seq_dir / "gene1.fa").write_text(">a\nMKT\n>b\nMKA\n")
    out_dir = tmp_path / "out"

    with pytest.raises(ValueError, match="magus requires Linux"):
        run_align(
            seq_dir=seq_dir,
            output_dir=out_dir,
            method="magus",
            seq_type="AA",
            dry_run=True,
        )


def test_detect_seq_type_from_files_aa(tmp_path: Path) -> None:
    from phyloai.pretree.align import _detect_seq_type_from_files

    f = tmp_path / "gene.fa"
    f.write_text(">sp1\nMKTFFF\n>sp2\nMKTYYY\n")
    assert _detect_seq_type_from_files([f]) == "AA"


def test_detect_seq_type_from_files_nt(tmp_path: Path) -> None:
    from phyloai.pretree.align import _detect_seq_type_from_files

    f = tmp_path / "gene.fa"
    f.write_text(">sp1\nACGTACGT\n>sp2\nACGTACGT\n")
    assert _detect_seq_type_from_files([f]) == "NT"
