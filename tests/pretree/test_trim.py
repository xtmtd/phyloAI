from __future__ import annotations

import shutil
from pathlib import Path

import pytest
from click.testing import CliRunner

from phyloai.cli.main import cli


def test_scan_input_finds_fasta_files(tmp_path: Path) -> None:
    from phyloai.pretree.trim import _scan_input

    (tmp_path / "gene1.fa").write_text(">a\nMKT\n")
    (tmp_path / "gene2.faa").write_text(">b\nMKT\n")
    (tmp_path / "notes.txt").write_text("skip")
    (tmp_path / "empty.fa").write_text("")
    (tmp_path / "subdir").mkdir()

    found, skipped = _scan_input(tmp_path)
    assert len(found) == 2
    reasons = {s["reason"] for s in skipped}
    assert "empty file" in reasons
    assert "directory" in reasons
    assert "unrecognized extension" in reasons


def test_build_trimal_cmd_aa_only(tmp_path: Path) -> None:
    from phyloai.pretree.trim import _build_trimal_cmd

    inp = tmp_path / "gene1.fa"
    out = tmp_path / "gene1_trim.fa"
    cmd = _build_trimal_cmd(inp, out, method="automated1", executable="trimal")

    assert cmd[0] == "trimal"
    assert "-in" in cmd
    assert str(inp) in cmd
    assert "-out" in cmd
    assert str(out) in cmd
    assert "-automated1" in cmd
    assert "-backtrans" not in cmd


def test_build_trimal_cmd_backtrans(tmp_path: Path) -> None:
    from phyloai.pretree.trim import _build_trimal_cmd

    inp = tmp_path / "gene1.fa"
    out = tmp_path / "gene1_trim.fa"
    nt = tmp_path / "gene1.fna"
    cmd = _build_trimal_cmd(inp, out, method="gappyout", executable="trimal", backtrans_path=nt)

    assert "-gappyout" in cmd
    assert "-backtrans" in cmd
    assert str(nt) in cmd
    assert "-ignorestopcodon" in cmd


def test_build_bmge_cmd_aa(tmp_path: Path) -> None:
    from phyloai.pretree.trim import _build_bmge_cmd

    inp = tmp_path / "gene1.fa"
    out = tmp_path / "gene1_trim.fa"
    cmd = _build_bmge_cmd(inp, out, seq_type="AA", matrix="BLOSUM62", entropy=0.5, java_executable="java", bmge_jar="/path/to/BMGE.jar")

    assert "java" in cmd
    assert "-jar" in cmd
    assert "/path/to/BMGE.jar" in cmd
    assert "-i" in cmd
    assert str(inp) in cmd
    assert "-t" in cmd
    assert "AA" in cmd
    assert "-m" in cmd
    assert "BLOSUM62" in cmd
    assert "-h" in cmd
    assert "0.5" in cmd
    assert "-of" in cmd
    assert str(out) in cmd


def test_build_bmge_cmd_codon(tmp_path: Path) -> None:
    from phyloai.pretree.trim import _build_bmge_cmd

    cmd = _build_bmge_cmd(
        tmp_path / "gene1.fna",
        tmp_path / "gene1_trim.fna",
        seq_type="CODON",
        matrix="BLOSUM90",
        entropy=0.4,
        java_executable="java",
        bmge_jar="/path/to/BMGE.jar",
    )

    assert "CODON" in cmd
    assert "BLOSUM90" in cmd


def test_build_bmge_cmd_nt(tmp_path: Path) -> None:
    from phyloai.pretree.trim import _build_bmge_cmd

    cmd = _build_bmge_cmd(
        tmp_path / "gene1.fna",
        tmp_path / "gene1_trim.fna",
        seq_type="NT",
        matrix="DNAPAM100:2",
        entropy=0.5,
        java_executable="java",
        bmge_jar="/path/to/BMGE.jar",
    )

    assert "DNA" in cmd
    assert "DNAPAM100:2" in cmd


def test_split_tool_args_rejects_managed_bmge_io_args() -> None:
    from phyloai.pretree.trim import _split_tool_args

    with pytest.raises(ValueError, match="PhyloAI-managed bmge"):
        _split_tool_args("-i custom.fa -m BLOSUM90", "bmge")

    assert _split_tool_args("-m BLOSUM90 -h 0.4", "bmge") == ["-m", "BLOSUM90", "-h", "0.4"]


def test_detect_trim_tool_versions_uses_registry_detection() -> None:
    from unittest.mock import patch

    from phyloai.core.env import ToolInfo, ToolStatus
    from phyloai.pretree.trim import _detect_trim_tool_versions

    def fake_detect(self, name, **_kwargs):
        versions = {"bmge": "1.12", "java": "17.0.15", "clipkit": "2.0.1"}
        return ToolInfo(name=name, status=ToolStatus.OK, path=Path(f"/fake/{name}"), version=versions[name])

    with patch("phyloai.core.env.ToolEnv._detect_tool", fake_detect), patch("shutil.which", return_value="/fake/java"):
        assert _detect_trim_tool_versions("bmge", None, None, None) == {"bmge": "1.12", "java": "17.0.15"}
        assert _detect_trim_tool_versions("clipkit", None, None, None) == {"clipkit": "2.0.1"}


def test_build_clipkit_cmd_aa_only(tmp_path: Path) -> None:
    from phyloai.pretree.trim import _build_clipkit_cmd

    inp = tmp_path / "gene1.fa"
    out = tmp_path / "gene1_trim.fa"
    cmd = _build_clipkit_cmd(inp, out, mode="smart-gap", codon=False, log_path=None, executable="clipkit")

    assert cmd[0] == "clipkit"
    assert str(inp) in cmd
    assert "-o" in cmd
    assert str(out) in cmd
    assert "-m" in cmd
    assert "smart-gap" in cmd
    assert "--codon" not in cmd
    assert "-l" not in cmd


def test_build_clipkit_cmd_codon_mode(tmp_path: Path) -> None:
    from phyloai.pretree.trim import _build_clipkit_cmd

    cmd = _build_clipkit_cmd(tmp_path / "gene1.fna", tmp_path / "gene1_trim.fna", mode="smart-gap", codon=True, log_path=None, executable="clipkit")
    assert "--codon" in cmd


def test_build_clipkit_cmd_with_log(tmp_path: Path) -> None:
    from phyloai.pretree.trim import _build_clipkit_cmd

    log = tmp_path / "gene1_trim.fa.log"
    cmd = _build_clipkit_cmd(tmp_path / "gene1.fa", tmp_path / "gene1_trim.fa", mode="smart-gap", codon=False, log_path=log, executable="clipkit")
    assert "-l" in cmd


def test_parse_clipkit_log(tmp_path: Path) -> None:
    from phyloai.pretree.trim import _parse_clipkit_log

    log_path = tmp_path / "test.log"
    log_path.write_text("1 keep constant 0.0\n2 trim other 0.9\n3 trim other 0.9\n4 keep constant 0.0\n5 keep constant 0.0\n")

    assert _parse_clipkit_log(log_path) == [0, 3, 4]


def test_project_columns_onto_nt_msa() -> None:
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord

    from phyloai.pretree.trim import _project_columns_onto_nt_msa

    codon_records = [
        SeqRecord(Seq("ATGGCTTCTACTAAA"), id="seq1", description=""),
        SeqRecord(Seq("ATGGCTTCTACT---"), id="seq2", description=""),
    ]
    result = _project_columns_onto_nt_msa(codon_records, [0, 3, 4])

    assert len(result) == 2
    assert str(result[0].seq) == "ATGACTAAA"
    assert str(result[1].seq) == "ATGACT---"


def test_translate_codon_msa() -> None:
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord

    from phyloai.pretree.trim import _translate_codon_msa

    codon_records = [
        SeqRecord(Seq("ATGGCTTCT"), id="seq1", description=""),
        SeqRecord(Seq("ATG------"), id="seq2", description=""),
    ]
    aa_records = _translate_codon_msa(codon_records)

    assert str(aa_records[0].seq) == "MAS"
    assert str(aa_records[1].seq) == "M--"
    assert aa_records[0].id == "seq1"


def test_translate_codon_msa_maps_invalid_or_ambiguous_codons_to_x() -> None:
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord

    from phyloai.pretree.trim import _translate_codon_msa

    aa_records = _translate_codon_msa([SeqRecord(Seq("ATGXXXNNN"), id="seq1", description="")])

    assert str(aa_records[0].seq) == "MXX"


def test_trim_one_trimal_dry_run(tmp_path: Path) -> None:
    from phyloai.pretree.trim import _trim_one_trimal

    msa = tmp_path / "gene1.faa"
    msa.write_text(">seq1\nMKTPQ\n>seq2\nMKTPQ\n")
    out_dir = tmp_path / "seqs"
    out_dir.mkdir()

    result = _trim_one_trimal(msa, out_dir, None, None, "automated1", "AA", None, True, "trimal")

    assert result["status"] == "dry_run"
    assert "cmd" in result


def test_trim_one_trimal_codon_dry_run(tmp_path: Path) -> None:
    from phyloai.pretree.trim import _trim_one_trimal

    msa = tmp_path / "gene1.fna"
    msa.write_text(">seq1\nATGGCTTCT\n>seq2\nATGGCT---\n")
    faa_dir = tmp_path / "seqs" / "faa"
    fna_dir = tmp_path / "seqs" / "fna"
    faa_dir.mkdir(parents=True)
    fna_dir.mkdir(parents=True)

    result = _trim_one_trimal(msa, faa_dir, fna_dir, None, "automated1", "CODON", None, True, "trimal")

    assert result["status"] == "dry_run"


def test_trim_one_bmge_dry_run(tmp_path: Path) -> None:
    from phyloai.pretree.trim import _trim_one_bmge

    msa = tmp_path / "gene1.faa"
    msa.write_text(">seq1\nMKTPQ\n>seq2\nMKTPQ\n")
    out_dir = tmp_path / "seqs"
    out_dir.mkdir()

    result = _trim_one_bmge(msa, out_dir, None, "AA", "BLOSUM62", 0.5, None, True, "java", "/fake/BMGE.jar")

    assert result["status"] == "dry_run"
    assert "cmd" in result


def test_trim_one_bmge_codon_dry_run(tmp_path: Path) -> None:
    from phyloai.pretree.trim import _trim_one_bmge

    msa = tmp_path / "gene1.fna"
    msa.write_text(">seq1\nATGGCTTCT\n>seq2\nATGGCT---\n")
    faa_dir = tmp_path / "seqs" / "faa"
    fna_dir = tmp_path / "seqs" / "fna"
    faa_dir.mkdir(parents=True)
    fna_dir.mkdir(parents=True)

    result = _trim_one_bmge(msa, faa_dir, fna_dir, "CODON", "BLOSUM62", 0.5, None, True, "java", "/fake/BMGE.jar")

    assert result["status"] == "dry_run"


def test_trim_one_clipkit_dry_run(tmp_path: Path) -> None:
    from phyloai.pretree.trim import _trim_one_clipkit

    msa = tmp_path / "gene1.faa"
    msa.write_text(">seq1\nMKTPQ\n>seq2\nMKTPQ\n")
    out_dir = tmp_path / "seqs"
    out_dir.mkdir()

    result = _trim_one_clipkit(msa, out_dir, None, None, "smart-gap", "AA", None, True, "clipkit")

    assert result["status"] == "dry_run"
    assert "cmd" in result


def test_trim_one_clipkit_method4_dry_run(tmp_path: Path) -> None:
    from phyloai.pretree.trim import _trim_one_clipkit

    msa = tmp_path / "gene1.faa"
    msa.write_text(">seq1\nMKTPQ\n>seq2\nMKTPQ\n")
    nt_msa = tmp_path / "gene1.fna"
    nt_msa.write_text(">seq1\nATGAAGACCCCTCAA\n>seq2\nATGAAGACCCCTCAA\n")
    faa_dir = tmp_path / "seqs" / "faa"
    fna_dir = tmp_path / "seqs" / "fna"
    faa_dir.mkdir(parents=True)
    fna_dir.mkdir(parents=True)

    result = _trim_one_clipkit(msa, faa_dir, fna_dir, nt_msa, "smart-gap", "AA", None, True, "clipkit")

    assert result["status"] == "dry_run"


def test_verify_trim_outputs_both_exist(tmp_path: Path) -> None:
    from phyloai.pretree.trim import verify_trim_outputs

    aa = tmp_path / "gene1.fa"
    nt = tmp_path / "gene1_nt.fa"
    aa.write_text(">seq1\nMKT\n>seq2\nMKT\n")
    nt.write_text(">seq1\nATGAAGACT\n>seq2\nATGAAGACT\n")

    assert verify_trim_outputs(aa, nt) is True


def test_verify_trim_outputs_missing_aa(tmp_path: Path) -> None:
    from phyloai.pretree.trim import verify_trim_outputs

    assert verify_trim_outputs(tmp_path / "gene1.fa", None) is False


def test_verify_trim_outputs_empty_aa(tmp_path: Path) -> None:
    from phyloai.pretree.trim import verify_trim_outputs

    aa = tmp_path / "gene1.fa"
    aa.write_text("")
    assert verify_trim_outputs(aa, None) is False


def test_detect_trim_seq_type_auto(tmp_path: Path) -> None:
    from phyloai.pretree.trim import _detect_seq_type_from_files

    f = tmp_path / "gene1.fa"
    f.write_text(">seq1\nMKTPQWER\n>seq2\nMKTPQWER\n")
    assert _detect_seq_type_from_files([f]) == "AA"


def test_run_trim_validation_codon_with_nt_dir_raises(tmp_path: Path) -> None:
    from phyloai.pretree.trim import run_trim

    msa_dir = tmp_path / "msa"
    msa_dir.mkdir()
    (msa_dir / "gene1.fa").write_text(">a\nATGGCT\n")

    with pytest.raises(ValueError, match="CODON mode does not use --nt-dir"):
        run_trim(msa_dir=msa_dir, output_dir=tmp_path / "out", tool="trimal", seq_type="CODON", nt_dir=tmp_path / "nt", threads=1)


def test_run_trim_validation_overwrite_resume_mutual_exclusive(tmp_path: Path) -> None:
    from phyloai.pretree.trim import run_trim

    msa_dir = tmp_path / "msa"
    msa_dir.mkdir()
    (msa_dir / "gene1.fa").write_text(">a\nMKT\n")

    with pytest.raises(ValueError, match="mutually exclusive"):
        run_trim(msa_dir=msa_dir, output_dir=tmp_path / "out", tool="trimal", seq_type="AA", overwrite=True, resume=True, threads=1)


def test_run_trim_dry_run_returns_payload(tmp_path: Path) -> None:
    from phyloai.pretree.trim import run_trim

    msa_dir = tmp_path / "msa"
    msa_dir.mkdir()
    (msa_dir / "gene1.fa").write_text(">a\nMKT\n>b\nMKT\n")

    payload = run_trim(msa_dir=msa_dir, output_dir=tmp_path / "out", tool="trimal", seq_type="AA", dry_run=True, threads=1)

    assert payload["status"] == "success"
    assert payload["data"]["summary"]["n_input_files"] == 1


def test_run_trim_resume_completed_checkpoint_returns_existing_success(tmp_path: Path) -> None:
    from phyloai.core.checkpoint import save_checkpoint_atomic
    from phyloai.pretree.checkpoint_helpers import build_initial_checkpoint, mark_task
    from phyloai.pretree.trim import run_trim

    msa_dir = tmp_path / "msa"
    out_dir = tmp_path / "out"
    msa_dir.mkdir()
    gene = msa_dir / "gene1.fa"
    gene.write_text(">a\nMKT\n>b\nMKT\n")
    aa_out = out_dir / "seqs" / "gene1.fa"
    aa_out.parent.mkdir(parents=True)
    aa_out.write_text(">a\nMKT\n>b\nMKT\n")
    params = {
        "msa_dir": str(msa_dir), "nt_dir": None, "seq_type": "AA",
        "tool": "trimal", "trimal_method": "automated1", "trimal_path": None,
        "bmge_matrix": None, "bmge_entropy": None, "bmge_path": None,
        "clipkit_method": None, "clipkit_path": None,
        "threads": 1,
        "tool_args": None, "output_dir": str(out_dir),
        "overwrite": False, "dry_run": False,
        "resume": True, "quiet": False,
    }
    checkpoint = build_initial_checkpoint(
        step="pretree.trim",
        command="phyloai pretree trim --msa-dir /data --output-dir /out --tool trimal --seq-type AA --threads 4",
        params=params,
        inputs=[gene],
        output_for=lambda _path: aa_out,
        nt_output_for=lambda _path: None,
    )
    mark_task(checkpoint, "gene1", status="success")
    checkpoint.status = "success"
    checkpoint.completed_at = checkpoint.touch()
    save_checkpoint_atomic(checkpoint, out_dir / "checkpoint.json")

    payload = run_trim(msa_dir=msa_dir, output_dir=out_dir, tool="trimal", seq_type="AA", resume=True, threads=1)

    assert payload["status"] == "success"
    assert payload["data"]["summary"]["n_trimmed"] == 1


def test_run_trim_resume_reconstructs_skipped_successes(tmp_path: Path) -> None:
    from phyloai.core.checkpoint import save_checkpoint_atomic
    from phyloai.pretree.checkpoint_helpers import build_initial_checkpoint, mark_task
    from phyloai.pretree.trim import run_trim

    msa_dir = tmp_path / "msa"
    out_dir = tmp_path / "out"
    msa_dir.mkdir()
    gene1 = msa_dir / "gene1.fa"
    gene2 = msa_dir / "gene2.fa"
    gene1.write_text(">a\nMKT\n>b\nMKT\n")
    gene2.write_text(">a\nMKT\n>b\nMKT\n")
    aa1 = out_dir / "seqs" / "gene1.fa"
    aa2 = out_dir / "seqs" / "gene2.fa"
    aa1.parent.mkdir(parents=True)
    aa1.write_text(">a\nMKT\n>b\nMKT\n")
    params = {
        "msa_dir": str(msa_dir), "nt_dir": None, "seq_type": "AA",
        "tool": "trimal", "trimal_method": "automated1", "trimal_path": None,
        "bmge_matrix": None, "bmge_entropy": None, "bmge_path": None,
        "clipkit_method": None, "clipkit_path": None,
        "threads": 1,
        "tool_args": None, "output_dir": str(out_dir),
        "overwrite": False, "dry_run": True,
        "resume": True, "quiet": False,
    }
    checkpoint = build_initial_checkpoint(
        step="pretree.trim",
        command="phyloai pretree trim --msa-dir /data --output-dir /out --tool trimal --seq-type AA --threads 4",
        params=params,
        inputs=[gene1, gene2],
        output_for=lambda path: aa1 if path.stem == "gene1" else aa2,
        nt_output_for=lambda _path: None,
    )
    mark_task(checkpoint, "gene1", status="success")
    mark_task(checkpoint, "gene2", status="failed", reason="previous failure")
    save_checkpoint_atomic(checkpoint, out_dir / "checkpoint.json")

    payload = run_trim(msa_dir=msa_dir, output_dir=out_dir, tool="trimal", seq_type="AA", resume=True, dry_run=True, threads=1)

    assert payload["data"]["summary"]["n_input_files"] == 2


def test_codon_workers_write_terminal_stop_stripped_temp_inputs(tmp_path: Path) -> None:
    import phyloai.pretree.trim as trim

    msa = tmp_path / "gene1.fna"
    msa.write_text(">seq1\nATGGCTTAA\n>seq2\nATGGCTTCT\n")
    faa_dir = tmp_path / "faa"
    fna_dir = tmp_path / "fna"
    faa_dir.mkdir()
    fna_dir.mkdir()
    captured: list[list[str]] = []
    captured_input_text: list[str] = []

    def fake_run(cmd: list[str]):
        captured.append(cmd)
        captured_input_text.append(Path(cmd[cmd.index("-i") + 1]).read_text())
        out = Path(cmd[cmd.index("-of") + 1] if "-of" in cmd else cmd[cmd.index("-o") + 1])
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(">seq1\nM\n>seq2\nM\n")
        class Proc:
            returncode = 0
            stdout = ""
            stderr = ""
        return Proc(), 0.01, None

    old = trim._run_cmd
    trim._run_cmd = fake_run
    try:
        result = trim._trim_one_bmge(msa, faa_dir, fna_dir, "CODON", "BLOSUM62", 0.5, None, False, "java", "BMGE.jar")
    finally:
        trim._run_cmd = old

    assert result["status"] == "success"
    assert captured
    assert "MA" in captured_input_text[0]
    assert "*" not in captured_input_text[0]


def test_run_trim_missing_nt_pair_marks_checkpoint_failed(tmp_path: Path) -> None:
    from phyloai.core.checkpoint import load_checkpoint
    from phyloai.pretree.trim import run_trim

    msa_dir = tmp_path / "msa"
    nt_dir = tmp_path / "nt"
    msa_dir.mkdir()
    nt_dir.mkdir()
    (msa_dir / "gene1.fa").write_text(">a\nMKT\n>b\nMKT\n")

    with pytest.raises(ValueError):
        run_trim(msa_dir=msa_dir, output_dir=tmp_path / "out", tool="clipkit", seq_type="AA", nt_dir=nt_dir, threads=1)

    checkpoint = load_checkpoint(tmp_path / "out" / "checkpoint.json")
    assert checkpoint.tasks[0].status == "failed"
    assert checkpoint.tasks[0].reason == "nt_pairing_missing"


def test_trim_cli_all_skipped_exits_two(tmp_path: Path) -> None:
    msa_dir = tmp_path / "msa"
    msa_dir.mkdir()
    (msa_dir / "gene1.fa").write_text(">a\nMKT\n")
    fake_trimal = tmp_path / "fake-trimal"
    fake_trimal.write_text("#!/bin/sh\nexit 1\n")
    fake_trimal.chmod(0o755)

    result = CliRunner().invoke(cli, [
        "pretree", "trim",
        "--msa-dir", str(msa_dir),
        "--output-dir", str(tmp_path / "out"),
        "--tool", "trimal",
        "--trimal-path", str(fake_trimal),
    ])

    assert result.exit_code == 2


def test_trim_cli_uses_manual_trimal_thresholds_from_tool_args(tmp_path: Path) -> None:
    msa_dir = tmp_path / "msa"
    msa_dir.mkdir()
    (msa_dir / "gene1.fa").write_text(">a\nMKT\n>b\nMKT\n")

    result = CliRunner().invoke(cli, [
        "pretree", "trim",
        "--msa-dir", str(msa_dir),
        "--output-dir", str(tmp_path / "out"),
        "--tool", "trimal",
        "--tool-args", "-gt 0.9 -cons 60",
        "--dry-run",
    ])

    assert result.exit_code == 0, result.output
    assert "-gt 0.9" in result.output
    assert "-cons 60" in result.output
    assert "-automated1" not in result.output


def test_trim_one_trimal_omits_auto_method_when_manual_thresholds_in_tool_args(tmp_path: Path) -> None:
    import phyloai.pretree.trim as trim

    msa = tmp_path / "gene1.faa"
    msa.write_text(">a\nMKT\n>b\nMKT\n")
    out_dir = tmp_path / "out"
    out_dir.mkdir()
    captured: list[list[str]] = []

    def fake_run(cmd: list[str]):
        captured.append(cmd)
        out = Path(cmd[cmd.index("-out") + 1])
        out.write_text(">a\nMKT\n>b\nMKT\n")
        class Proc:
            returncode = 0
            stdout = ""
            stderr = ""
        return Proc(), 0.01, None

    old = trim._run_cmd
    trim._run_cmd = fake_run
    try:
        result = trim._trim_one_trimal(
            msa, out_dir, None, None, "automated1", "AA", "-gt 0.9 -cons 60", False, "trimal"
        )
    finally:
        trim._run_cmd = old

    assert result["status"] == "success"
    assert "-automated1" not in captured[0]
    assert "-gt" in captured[0]
    assert "-cons" in captured[0]


def test_trim_one_trimal_backtrans_strips_gaps_from_nt_input(tmp_path: Path) -> None:
    import phyloai.pretree.trim as trim

    msa = tmp_path / "gene1.faa"
    msa.write_text(">a\nMKT\n>b\nMKT\n")
    nt_msa = tmp_path / "gene1.fna"
    nt_msa.write_text(">a\nATGAAA---ACT\n>b\nATGAAA---ACT\n")
    faa_dir = tmp_path / "faa"
    fna_dir = tmp_path / "fna"
    faa_dir.mkdir()
    fna_dir.mkdir()
    captured_backtrans_input: list[str] = []

    def fake_run(cmd: list[str]):
        out = Path(cmd[cmd.index("-out") + 1])
        out.parent.mkdir(parents=True, exist_ok=True)
        if "-backtrans" in cmd:
            backtrans_input = Path(cmd[cmd.index("-backtrans") + 1])
            captured_backtrans_input.append(backtrans_input.read_text())
            out.write_text(">a\nATGAAAACT\n>b\nATGAAAACT\n")
        else:
            out.write_text(">a\nMKT\n>b\nMKT\n")

        class Proc:
            returncode = 0
            stdout = ""
            stderr = ""

        return Proc(), 0.01, None

    old = trim._run_cmd
    trim._run_cmd = fake_run
    try:
        result = trim._trim_one_trimal(msa, faa_dir, fna_dir, nt_msa, "automated1", "AA", None, False, "trimal")
    finally:
        trim._run_cmd = old

    assert result["status"] == "success"
    assert captured_backtrans_input
    assert "---" not in captured_backtrans_input[0]
    assert "ATGAAAACT" in captured_backtrans_input[0]


def test_bmge_codon_output_is_normalized_to_three_times_aa(tmp_path: Path) -> None:
    import phyloai.pretree.trim as trim

    msa = tmp_path / "gene1.fna"
    msa.write_text(">seq1\nATGGCTTCT\n>seq2\nATGGCTTCT\n")
    faa_dir = tmp_path / "faa"
    fna_dir = tmp_path / "fna"
    faa_dir.mkdir()
    fna_dir.mkdir()

    def fake_run(cmd: list[str]):
        out = Path(cmd[cmd.index("-of") + 1])
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(">seq1\nMA\n>seq2\nMA\n")
        class Proc:
            returncode = 0
            stdout = ""
            stderr = ""
        return Proc(), 0.01, None

    old = trim._run_cmd
    trim._run_cmd = fake_run
    try:
        result = trim._trim_one_bmge(msa, faa_dir, fna_dir, "CODON", "BLOSUM62", 0.5, None, False, "java", "BMGE.jar")
    finally:
        trim._run_cmd = old

    assert result["status"] == "success"
    nt_len = len((fna_dir / "gene1.fa").read_text().splitlines()[1])
    aa_len = len((faa_dir / "gene1.fa").read_text().splitlines()[1])
    assert nt_len == aa_len * 3


def test_bmge_codon_trims_translated_aa_then_projects_nt(tmp_path: Path) -> None:
    import phyloai.pretree.trim as trim

    msa = tmp_path / "gene1.fna"
    msa.write_text(">seq1\nATGGCTTCT\n>seq2\nATGGCTTCT\n")
    faa_dir = tmp_path / "faa"
    fna_dir = tmp_path / "fna"
    faa_dir.mkdir()
    fna_dir.mkdir()
    captured: list[list[str]] = []

    def fake_run(cmd: list[str]):
        captured.append(cmd)
        inp = Path(cmd[cmd.index("-i") + 1])
        assert inp.read_text().splitlines()[1] == "MAS"
        out = Path(cmd[cmd.index("-of") + 1])
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(">seq1\nMS\n>seq2\nMS\n")
        class Proc:
            returncode = 0
            stdout = ""
            stderr = ""
        return Proc(), 0.01, None

    old = trim._run_cmd
    trim._run_cmd = fake_run
    try:
        result = trim._trim_one_bmge(msa, faa_dir, fna_dir, "CODON", "BLOSUM62", 0.5, None, False, "java", "BMGE.jar")
    finally:
        trim._run_cmd = old

    assert result["status"] == "success"
    assert "AA" in captured[0]
    assert "CODON" not in captured[0]
    assert (fna_dir / "gene1.fa").read_text().splitlines()[1] == "ATGTCT"


def test_bmge_success_return_without_output_is_skipped(tmp_path: Path) -> None:
    import phyloai.pretree.trim as trim

    msa = tmp_path / "gene1.faa"
    msa.write_text(">seq1\nMKT\n>seq2\nMKT\n")
    out_dir = tmp_path / "out"
    out_dir.mkdir()

    def fake_run(_cmd: list[str]):
        class Proc:
            returncode = 0
            stdout = ""
            stderr = ""
        return Proc(), 0.01, None

    old = trim._run_cmd
    trim._run_cmd = fake_run
    try:
        result = trim._trim_one_bmge(msa, out_dir, None, "AA", "BLOSUM62", 0.5, None, False, "java", "BMGE.jar")
    finally:
        trim._run_cmd = old

    assert result["status"] == "skipped"
    assert "did not create output" in result["reason"]


def test_trim_help_documents_tool_options() -> None:
    result = CliRunner().invoke(cli, ["pretree", "trim", "--help"])
    normalized = " ".join(result.output.split())

    assert result.exit_code == 0
    assert "AA+NT mode: --msa-dir is aligned AA MSAs" in normalized
    assert "BLOSUM30" in normalized
    assert "more stringent" in normalized
    assert "smart-gap|entropy|gappy|block-gappy" in normalized
    assert "-g 0.8" in normalized
    assert "--tool-args" in normalized


@pytest.mark.skipif(shutil.which("trimal") is None, reason="trimal not available")
def test_trim_trimal_aa_integration(tmp_path: Path) -> None:
    from Bio import SeqIO

    from phyloai.pretree.trim import run_trim

    msa_dir = tmp_path / "msa"
    msa_dir.mkdir()
    aa_seqs = {
        "seq1": "MASTKLIVDE", "seq2": "M----LIVDE", "seq3": "M----LIVDE",
        "seq4": "M----LIVDE", "seq5": "M----LIVDE", "seq6": "M----LIVDE",
        "seq7": "M----LIVDE", "seq8": "M----LIVDE", "seq9": "M----LIVDE",
        "seq10": "M----LIVDE",
    }
    (msa_dir / "gene1.faa").write_text("\n".join(f">{name}\n{seq}" for name, seq in aa_seqs.items()))

    payload = run_trim(msa_dir=msa_dir, output_dir=tmp_path / "out", tool="trimal", seq_type="AA", trimal_method="gappyout", threads=1)

    assert payload["status"] == "success"
    assert payload["data"]["summary"]["n_trimmed"] == 1
    out_file = tmp_path / "out" / "seqs" / "gene1.fa"
    records = list(SeqIO.parse(str(out_file), "fasta"))
    assert len(records) == 10
    assert len(str(records[0].seq)) < 10


@pytest.mark.skipif(shutil.which("clipkit") is None, reason="clipkit not available")
def test_trim_clipkit_aa_integration(tmp_path: Path) -> None:
    from phyloai.pretree.trim import run_trim

    msa_dir = tmp_path / "msa"
    msa_dir.mkdir()
    aa_seqs = {
        "seq1": "MASTKLIVDE", "seq2": "M----LIVDE", "seq3": "M----LIVDE",
        "seq4": "M----LIVDE", "seq5": "M----LIVDE", "seq6": "M----LIVDE",
        "seq7": "M----LIVDE", "seq8": "M----LIVDE", "seq9": "M----LIVDE",
        "seq10": "M----LIVDE",
    }
    (msa_dir / "gene1.faa").write_text("\n".join(f">{name}\n{seq}" for name, seq in aa_seqs.items()))

    payload = run_trim(msa_dir=msa_dir, output_dir=tmp_path / "out", tool="clipkit", seq_type="AA", clipkit_method="gappy", threads=1)

    assert payload["status"] == "success"
    assert payload["data"]["summary"]["n_trimmed"] == 1


@pytest.mark.skipif(shutil.which("clipkit") is None, reason="clipkit not available")
def test_trim_clipkit_method4_integration(tmp_path: Path) -> None:
    from Bio import SeqIO as _SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord

    from phyloai.pretree.trim import run_trim

    msa_dir = tmp_path / "msa"
    nt_dir = tmp_path / "nt"
    msa_dir.mkdir()
    nt_dir.mkdir()
    codon_map = {"M": "ATG", "A": "GCT", "S": "TCT", "T": "ACT", "K": "AAA", "L": "CTT", "I": "ATT", "V": "GTT", "D": "GAT", "E": "GAA", "-": "---"}
    aa_seqs = {
        "seq1": "MASTKLIVDE", "seq2": "M----LIVDE", "seq3": "M----LIVDE",
        "seq4": "M----LIVDE", "seq5": "M----LIVDE", "seq6": "M----LIVDE",
        "seq7": "M----LIVDE", "seq8": "M----LIVDE", "seq9": "M----LIVDE",
        "seq10": "M----LIVDE",
    }
    aa_recs = [SeqRecord(Seq(seq), id=name, description="") for name, seq in aa_seqs.items()]
    codon_recs = [SeqRecord(Seq("".join(codon_map[aa] for aa in seq)), id=name, description="") for name, seq in aa_seqs.items()]
    _SeqIO.write(aa_recs, str(msa_dir / "gene1.faa"), "fasta")
    _SeqIO.write(codon_recs, str(nt_dir / "gene1.fna"), "fasta")

    payload = run_trim(msa_dir=msa_dir, output_dir=tmp_path / "out", tool="clipkit", seq_type="AA", nt_dir=nt_dir, clipkit_method="gappy", threads=1)

    assert payload["status"] == "success"
    aa_out = list(_SeqIO.parse(str(tmp_path / "out" / "seqs" / "faa" / "gene1.fa"), "fasta"))
    nt_out = list(_SeqIO.parse(str(tmp_path / "out" / "seqs" / "fna" / "gene1.fa"), "fasta"))
    assert len(str(nt_out[0].seq)) == len(str(aa_out[0].seq)) * 3
