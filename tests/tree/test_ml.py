from __future__ import annotations

from pathlib import Path


def test_scan_input_finds_fasta_phylip_files(tmp_path: Path) -> None:
    from phyloai.tree.ml import _scan_input

    (tmp_path / "gene1.fa").write_text(">a\nACGT\n")
    (tmp_path / "gene2.faa").write_text(">b\nMKT\n")
    (tmp_path / "gene3.phy").write_text("2 10\na  ACGT\nb  ACGT\n")
    (tmp_path / "gene4.phylip").write_text("2 10\na  ACGT\nb  ACGT\n")
    (tmp_path / "gene5.nex").write_text("#NEXUS\n")
    (tmp_path / "notes.txt").write_text("skip")
    (tmp_path / "empty.fa").write_text("")
    (tmp_path / "subdir").mkdir()

    found, skipped = _scan_input(tmp_path)

    assert len(found) == 4
    assert len(skipped) == 4
    skip_reasons = {s["reason"] for s in skipped}
    assert "NEXUS format not supported by FastTree; use pretree convert first" in skip_reasons
    assert "empty file" in skip_reasons
    assert "directory" in skip_reasons
    assert "unrecognized extension: .txt" in skip_reasons


import pytest


def test_build_fasttree_cmd_aa_lg_full(tmp_path: Path) -> None:
    from phyloai.tree.ml import _build_fasttree_cmd

    inp = tmp_path / "gene.fa"
    out = tmp_path / "gene.tre"
    cmd = _build_fasttree_cmd(inp, out, seq_type="AA", model="lg", mode="normal",
                              boot=1000, cat=20, gamma=True)

    assert cmd[0] == "FastTree"
    assert "-lg" in cmd
    assert "-gamma" in cmd
    assert "-cat" in cmd and "20" in cmd
    assert "-boot" in cmd and "1000" in cmd
    assert "-nosupport" not in cmd
    assert str(inp) == cmd[-1]


def test_build_fasttree_cmd_nt_gtr(tmp_path: Path) -> None:
    from phyloai.tree.ml import _build_fasttree_cmd

    inp = tmp_path / "gene.fa"
    out = tmp_path / "gene.tre"
    cmd = _build_fasttree_cmd(inp, out, seq_type="NT", model="gtr")

    assert "-nt" in cmd
    assert "-gtr" in cmd
    assert "-lg" not in cmd


def test_build_fasttree_cmd_aa_jtt_default_no_flags(tmp_path: Path) -> None:
    from phyloai.tree.ml import _build_fasttree_cmd

    inp = tmp_path / "gene.fa"
    out = tmp_path / "gene.tre"
    cmd = _build_fasttree_cmd(inp, out, seq_type="AA", model="jtt")

    assert "-lg" not in cmd
    assert "-wag" not in cmd


def test_build_fasttree_cmd_fastest_mode(tmp_path: Path) -> None:
    from phyloai.tree.ml import _build_fasttree_cmd

    inp = tmp_path / "gene.fa"
    out = tmp_path / "gene.tre"
    cmd = _build_fasttree_cmd(inp, out, mode="fastest")

    assert "-fastest" in cmd
    assert "-slow" not in cmd


def test_build_fasttree_cmd_boot_zero_gives_nosupport(tmp_path: Path) -> None:
    from phyloai.tree.ml import _build_fasttree_cmd

    inp = tmp_path / "gene.fa"
    out = tmp_path / "gene.tre"
    cmd = _build_fasttree_cmd(inp, out, boot=0)

    assert "-nosupport" in cmd
    assert "-boot" not in cmd


def test_build_fasttree_cmd_no_gamma(tmp_path: Path) -> None:
    from phyloai.tree.ml import _build_fasttree_cmd

    inp = tmp_path / "gene.fa"
    out = tmp_path / "gene.tre"
    cmd = _build_fasttree_cmd(inp, out, gamma=False)

    assert "-gamma" not in cmd


def test_build_fasttree_cmd_with_tool_args(tmp_path: Path) -> None:
    from phyloai.tree.ml import _build_fasttree_cmd

    inp = tmp_path / "gene.fa"
    out = tmp_path / "gene.tre"
    cmd = _build_fasttree_cmd(inp, out, boot=1000, tool_args="-spr 4 -mlacc 2")

    assert "-spr" in cmd
    assert "4" in cmd
    assert "-mlacc" in cmd
    assert "2" in cmd


def test_check_managed_flag_conflict_blocks_lg() -> None:
    from phyloai.tree.ml import _check_managed_flag_conflict

    with pytest.raises(ValueError, match="Blocked managed flag.*-lg"):
        _check_managed_flag_conflict("-lg")


def test_check_managed_flag_conflict_blocks_boot() -> None:
    from phyloai.tree.ml import _check_managed_flag_conflict

    with pytest.raises(ValueError, match="Blocked managed flag.*-boot"):
        _check_managed_flag_conflict("-boot 500")


def test_check_managed_flag_conflict_allows_strategy_args() -> None:
    from phyloai.tree.ml import _check_managed_flag_conflict

    _check_managed_flag_conflict("-spr 4 -mlacc 2 -slownni")


def test_build_fasttree_cmd_with_explicit_executable(tmp_path: Path) -> None:
    from phyloai.tree.ml import _build_fasttree_cmd

    inp = tmp_path / "gene.fa"
    out = tmp_path / "gene.tre"
    cmd = _build_fasttree_cmd(inp, out, executable="/opt/bin/FastTree")

    assert cmd[0] == "/opt/bin/FastTree"


def test_run_one_fasttree_success(tmp_path: Path) -> None:
    from phyloai.tree.ml import _run_one_fasttree

    inp = tmp_path / "gene.fa"
    inp.write_text(">a\nMKTLLL\n>b\nMKTLLL\n")
    log_dir = tmp_path / "logs"
    log_dir.mkdir()

    result = _run_one_fasttree(
        gene_path=inp, seq_type="AA", model="lg", mode="normal",
        boot=1000, cat=20, gamma=True, tool_args=None,
        log_dir=log_dir, fasttree_executable="FastTree",
    )

    assert result["status"] == "success"
    assert "output_tree" in result
    assert "log_file" in result


def test_run_one_fasttree_dry_run(tmp_path: Path) -> None:
    from phyloai.tree.ml import _run_one_fasttree

    inp = tmp_path / "gene.fa"
    inp.write_text(">a\nMKT\n")
    log_dir = tmp_path / "logs"
    log_dir.mkdir()

    result = _run_one_fasttree(
        gene_path=inp, seq_type="AA", model="lg", mode="normal",
        boot=1000, cat=20, gamma=True, tool_args=None,
        log_dir=log_dir, fasttree_executable="FastTree",
        dry_run=True,
    )

    assert result["status"] == "dry_run"
    assert "cmd" in result


def test_run_one_fasttree_missing_input(tmp_path: Path) -> None:
    from phyloai.tree.ml import _run_one_fasttree

    inp = tmp_path / "missing.fa"
    log_dir = tmp_path / "logs"
    log_dir.mkdir()

    result = _run_one_fasttree(
        gene_path=inp, seq_type="AA", model="lg", mode="normal",
        boot=1000, cat=20, gamma=True, tool_args=None,
        log_dir=log_dir, fasttree_executable="FastTree",
    )

    assert result["status"] == "failed"
    assert "reason" in result


def test_validate_seq_types_homogeneous_aa(tmp_path: Path) -> None:
    from phyloai.tree.ml import _validate_seq_types

    (tmp_path / "g1.fa").write_text(">a\nMKTLLL\n")
    (tmp_path / "g2.fa").write_text(">b\nMKTLLL\n")
    files = sorted(tmp_path.glob("*.fa"))

    resolved, offending = _validate_seq_types(files, declared_type=None)

    assert resolved == "AA"
    assert len(offending) == 0


def test_validate_seq_types_mixed_raises(tmp_path: Path) -> None:
    from phyloai.tree.ml import _validate_seq_types

    (tmp_path / "g1.fa").write_text(">a\nMKTLLL\n")
    (tmp_path / "g2.fa").write_text(">b\nACGTAC\n")
    files = sorted(tmp_path.glob("*.fa"))

    resolved, offending = _validate_seq_types(files, declared_type=None)

    assert resolved is None
    assert len(offending) >= 1


def test_validate_seq_types_explicit_mismatch(tmp_path: Path) -> None:
    from phyloai.tree.ml import _validate_seq_types

    (tmp_path / "g1.fa").write_text(">a\nMKTLLL\n")
    files = sorted(tmp_path.glob("*.fa"))

    resolved, offending = _validate_seq_types(files, declared_type="NT")

    assert resolved == "NT"
    assert len(offending) == 1


def test_validate_seq_types_no_files(tmp_path: Path) -> None:
    from phyloai.tree.ml import _validate_seq_types

    resolved, offending = _validate_seq_types([], declared_type=None)

    assert resolved == "AA"
    assert len(offending) == 0


def test_run_fasttree_batch_dry_run(tmp_path: Path) -> None:
    from phyloai.tree.ml import run_fasttree

    msa_dir = tmp_path / "msas"
    msa_dir.mkdir()
    (msa_dir / "g1.fa").write_text(">a\nMKTLLL\n>b\nMKTLLL\n")
    (msa_dir / "g2.fa").write_text(">c\nMKTLLL\n>d\nMKTLLL\n")

    out_dir = tmp_path / "out"

    payload = run_fasttree(
        msa_dir=msa_dir, output_dir=out_dir,
        seq_type="AA", model="lg", dry_run=True, quiet=True,
    )

    assert payload["status"] == "success"
    assert payload["data"]["summary"]["n_input_files"] >= 2
    assert "files" in payload["data"]
    for f in payload["data"]["files"]:
        assert "cmd" in f


def test_run_fasttree_single_dry_run(tmp_path: Path) -> None:
    from phyloai.tree.ml import run_fasttree

    mat = tmp_path / "matrix.fa"
    mat.write_text(">a\nMKTLLL\n>b\nMKTLLL\n")

    out_dir = tmp_path / "out"

    payload = run_fasttree(
        matrix=mat, output_dir=out_dir,
        seq_type="AA", model="lg", dry_run=True, quiet=True,
    )

    assert payload["data"]["summary"]["mode"] == "--matrix"
    assert len(payload["data"]["files"]) >= 1
    if payload["data"]["files"]:
        assert "cmd" in payload["data"]["files"][0]


def test_run_fasttree_batch_auto_detects_mixed_and_fails(tmp_path: Path) -> None:
    from phyloai.tree.ml import run_fasttree

    msa_dir = tmp_path / "mixed"
    msa_dir.mkdir()
    (msa_dir / "aa.fa").write_text(">a\nMKTLLL\n")
    (msa_dir / "nt.fa").write_text(">b\nACGTAC\n")

    out_dir = tmp_path / "out"

    with pytest.raises(ValueError, match="Mixed sequence types"):
        run_fasttree(msa_dir=msa_dir, output_dir=out_dir, seq_type="auto", quiet=True)


def test_run_fasttree_batch_explicit_mismatch_fails(tmp_path: Path) -> None:
    from phyloai.tree.ml import run_fasttree

    msa_dir = tmp_path / "dir"
    msa_dir.mkdir()
    (msa_dir / "aa.fa").write_text(">a\nMKTLLL\n")

    out_dir = tmp_path / "out"

    with pytest.raises(ValueError, match="Files with wrong --seq-type"):
        run_fasttree(msa_dir=msa_dir, output_dir=out_dir, seq_type="NT", quiet=True)


def test_run_fasttree_invalid_model_for_aa(tmp_path: Path) -> None:
    from phyloai.tree.ml import run_fasttree

    msa_dir = tmp_path / "dir"
    msa_dir.mkdir()
    (msa_dir / "g1.fa").write_text(">a\nMKTLLL\n")

    out_dir = tmp_path / "out"

    with pytest.raises(ValueError, match="Invalid model for AA.*gtr"):
        run_fasttree(msa_dir=msa_dir, output_dir=out_dir, seq_type="AA", model="gtr", quiet=True)


def test_run_fasttree_neither_input_raises() -> None:
    from phyloai.tree.ml import run_fasttree
    from pathlib import Path

    out_dir = Path("/tmp/out")
    with pytest.raises(ValueError, match="Either --msa-dir or --matrix"):
        run_fasttree(output_dir=out_dir, seq_type="AA", model="lg", quiet=True)
