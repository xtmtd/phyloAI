"""Tests for dating_hessian pure helpers."""
from __future__ import annotations
from pathlib import Path
import pytest
from phyloai.posttree.dating_hessian import (
    detect_seqtype_from_alignment,
    count_partitions,
    validate_root_age,
    build_iqtree_dating_cmd,
    HESSIAN_OUTPUT_FILES,
)


def test_detect_aa(tmp_path):
    fa = tmp_path / "aa.fa"
    fa.write_text(">sp1\nMKTVFLGEI\n>sp2\nMLTVFLGEI\n")
    assert detect_seqtype_from_alignment(fa) == "AA"


def test_detect_nt(tmp_path):
    fa = tmp_path / "nt.fa"
    fa.write_text(">sp1\nACGTACGT\n>sp2\nACGTACGT\n")
    assert detect_seqtype_from_alignment(fa) == "NT"


def test_detect_auto_defaults_aa_for_mixed(tmp_path):
    fa = tmp_path / "m.fa"
    fa.write_text(">sp1\nACGTMKLWI\n")
    assert detect_seqtype_from_alignment(fa) == "AA"


def test_detect_auto_reads_phylip(tmp_path):
    phy = tmp_path / "aa.phy"
    phy.write_text("2 9\nsp1  MKTVFLGEI\nsp2  MLTVFLGEI\n")
    assert detect_seqtype_from_alignment(phy) == "AA"


def test_detect_auto_reads_nexus(tmp_path):
    nex = tmp_path / "nt.nex"
    nex.write_text(
        "#NEXUS\n"
        "begin data;\n"
        "dimensions ntax=2 nchar=8;\n"
        "format datatype=dna gap=- missing=?;\n"
        "matrix\n"
        "sp1 ACGTACGT\n"
        "sp2 ACGTACGT\n"
        ";\n"
        "end;\n"
    )
    assert detect_seqtype_from_alignment(nex) == "NT"


def test_count_raxml_partitions(tmp_path):
    pf = tmp_path / "parts.txt"
    pf.write_text("LG, p1 = 1-100\nLG, p2 = 101-200\nLG, p3 = 201-300\n")
    assert count_partitions(pf) == 3


def test_count_nexus_partitions(tmp_path):
    pf = tmp_path / "parts.nex"
    pf.write_text(
        "#NEXUS\nbegin sets;\n"
        "  charset p1 = 1-100;\n"
        "  charset p2 = 101-200;\n"
        "end;\n"
    )
    assert count_partitions(pf) == 2


def test_valid_root_age_upper_only():
    tree = "(A,(B,C))'<4.2';"
    assert validate_root_age(tree) is True


def test_valid_root_age_range():
    tree = "(A,(B,C))'>3.1<4.2';"
    assert validate_root_age(tree) is True


def test_missing_root_age():
    tree = "(A,(B,C));"
    assert validate_root_age(tree) is False


def test_unpartitioned_aa_default_model(tmp_path):
    matrix = tmp_path / "m.fa"
    matrix.touch()
    tree = tmp_path / "t.nwk"
    tree.touch()
    cmd = build_iqtree_dating_cmd(
        iqtree_path=Path("/usr/bin/iqtree3"),
        matrix=matrix,
        rooted_tree=tree,
        seq_type="AA",
        model_expr=None,
        partitions=None,
        n_partitions=0,
        threads=4,
        tool_args=None,
    )
    assert "-m" in cmd
    idx = cmd.index("-m")
    assert cmd[idx + 1] == "LG+F+G4"
    assert "--dating" in cmd
    assert "mcmctree" in cmd


def test_unpartitioned_nt_default_model(tmp_path):
    matrix = tmp_path / "m.fa"
    matrix.touch()
    tree = tmp_path / "t.nwk"
    tree.touch()
    cmd = build_iqtree_dating_cmd(
        iqtree_path=Path("/usr/bin/iqtree3"),
        matrix=matrix,
        rooted_tree=tree,
        seq_type="NT",
        model_expr=None,
        partitions=None,
        n_partitions=0,

        threads=4,
        tool_args=None,
    )
    idx = cmd.index("-m")
    assert cmd[idx + 1] == "GTR+G4"


def test_partitioned_aa_few(tmp_path):
    matrix = tmp_path / "m.fa"
    matrix.touch()
    tree = tmp_path / "t.nwk"
    tree.touch()
    parts = tmp_path / "parts.nex"
    parts.touch()
    cmd = build_iqtree_dating_cmd(
        iqtree_path=Path("/usr/bin/iqtree3"),
        matrix=matrix,
        rooted_tree=tree,
        seq_type="AA",
        model_expr=None,
        partitions=parts,
        n_partitions=5,

        threads=4,
        tool_args=None,
    )
    assert "-Q" in cmd
    assert "--merge" not in cmd
    assert "--mset" in cmd
    idx = cmd.index("--mset")
    assert cmd[idx + 1] == "LG"


def test_partitioned_aa_many_merges(tmp_path):
    matrix = tmp_path / "m.fa"
    matrix.touch()
    tree = tmp_path / "t.nwk"
    tree.touch()
    parts = tmp_path / "parts.nex"
    parts.touch()
    cmd = build_iqtree_dating_cmd(
        iqtree_path=Path("/usr/bin/iqtree3"),
        matrix=matrix,
        rooted_tree=tree,
        seq_type="AA",
        model_expr=None,
        partitions=parts,
        n_partitions=12,

        threads=4,
        tool_args=None,
    )
    assert "--merge" in cmd
    assert "--rclusterf" in cmd


def test_tool_args_appended_last(tmp_path):
    matrix = tmp_path / "m.fa"
    matrix.touch()
    tree = tmp_path / "t.nwk"
    tree.touch()
    cmd = build_iqtree_dating_cmd(
        iqtree_path=Path("/usr/bin/iqtree3"),
        matrix=matrix,
        rooted_tree=tree,
        seq_type="AA",
        model_expr=None,
        partitions=None,
        n_partitions=0,

        threads=4,
        tool_args="--redo",
    )
    assert cmd[-1] == "--redo"


def test_hessian_output_files_constant():
    assert "iqtree.dummy.phy" in HESSIAN_OUTPUT_FILES
    assert "iqtree.rooted.nwk" in HESSIAN_OUTPUT_FILES
    assert "iqtree.mcmctree.hessian" in HESSIAN_OUTPUT_FILES


def test_run_hessian_fails_when_output_files_missing(tmp_path):
    from unittest.mock import patch
    from phyloai.posttree.dating_hessian import run_hessian
    import subprocess

    matrix = tmp_path / "m.fa"
    matrix.write_text(">sp1\nMKTV\n>sp2\nMLTV\n")
    tree = tmp_path / "t.nwk"
    tree.write_text("(sp1,sp2)'<4.2';\n")
    output = tmp_path / "out"
    output.mkdir()

    def _fake_run(cmd, **kwargs):
        return subprocess.CompletedProcess(cmd, returncode=0, stdout="", stderr="")

    with patch("subprocess.run", side_effect=_fake_run), \
         patch("phyloai.posttree.dating_hessian._resolve_iqtree_path",
               return_value="/usr/bin/iqtree3"), \
         patch("phyloai.posttree.dating_hessian._detect_iqtree_version",
               return_value={"iqtree3": "2.0.0"}):
        payload = run_hessian(
            matrix=matrix, rooted_tree=tree,
            output_dir=output, dry_run=False, quiet=True,
        )
    assert payload["status"] == "error"


def test_run_hessian_warns_on_empty_hessian_file(tmp_path):
    from unittest.mock import patch
    from phyloai.posttree.dating_hessian import run_hessian, HESSIAN_OUTPUT_FILES
    import subprocess

    matrix = tmp_path / "m.fa"
    matrix.write_text(">sp1\nMKTV\n>sp2\nMLTV\n")
    tree = tmp_path / "t.nwk"
    tree.write_text("(sp1,sp2)'<4.2';\n")
    output = tmp_path / "out"
    output.mkdir()

    def _fake_run(cmd, **kwargs):
        for fn in HESSIAN_OUTPUT_FILES:
            p = output / fn
            if fn == "iqtree.mcmctree.hessian":
                p.write_text("")
            else:
                p.write_text("x")
        return subprocess.CompletedProcess(cmd, returncode=0, stdout="", stderr="")

    with patch("subprocess.run", side_effect=_fake_run), \
         patch("phyloai.posttree.dating_hessian._resolve_iqtree_path",
               return_value="/usr/bin/iqtree3"), \
         patch("phyloai.posttree.dating_hessian._detect_iqtree_version",
               return_value={"iqtree3": "2.0.0"}):
        payload = run_hessian(
            matrix=matrix, rooted_tree=tree,
            output_dir=output, dry_run=False, quiet=True,
        )
    assert payload["status"] == "error"
    assert any("empty" in str(w).lower() for w in payload.get("data", {}).get("warnings", []))
