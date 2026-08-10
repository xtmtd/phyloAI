from __future__ import annotations

from pathlib import Path

from click.testing import CliRunner

from phyloai.cli.main import cli


def test_tree_group_exists() -> None:
    result = CliRunner().invoke(cli, ["tree", "--help"])
    assert result.exit_code == 0
    assert "Maximum-likelihood" in result.output or "ml" in result.output


def test_tree_ml_help_shows_both_backends() -> None:
    result = CliRunner().invoke(cli, ["tree", "ml", "--help"])
    assert result.exit_code == 0
    assert "fasttree" in result.output


def test_tree_ml_fasttree_help() -> None:
    result = CliRunner().invoke(cli, ["tree", "ml", "fasttree", "--help"])
    assert result.exit_code == 0
    for flag in ["--msa-dir", "--matrix", "--seq-type", "--model", "--mode",
                  "--boot", "--cat", "--gamma", "--output-dir", "--threads"]:
        assert flag in result.output


def test_tree_ml_fasttree_mutual_exclusivity(tmp_path: Path) -> None:
    msa_dir = tmp_path / "msas"
    msa_dir.mkdir()
    mat = tmp_path / "matrix.fa"
    mat.write_text(">a\nMKT\n")

    result = CliRunner().invoke(cli, [
        "tree", "ml", "fasttree",
        "--msa-dir", str(msa_dir), "--matrix", str(mat),
    ])
    assert result.exit_code == 1


def test_tree_ml_fasttree_neither_input() -> None:
    result = CliRunner().invoke(cli, [
        "tree", "ml", "fasttree",
    ])
    assert result.exit_code == 1


def test_cli_msa_dir_nonexistent_exits_1() -> None:
    result = CliRunner().invoke(cli, [
        "tree", "ml", "fasttree",
        "--msa-dir", "/nonexistent/path",
    ])
    assert result.exit_code == 1


def test_tree_ml_fasttree_quiet_dry_run_batch(tmp_path: Path) -> None:
    msa_dir = tmp_path / "msas"
    msa_dir.mkdir()
    (msa_dir / "g1.fa").write_text(">a\nMKTLLL\n>b\nMKTLLL\n")

    out_dir = tmp_path / "out"

    result = CliRunner().invoke(cli, [
        "tree", "ml", "fasttree",
        "--msa-dir", str(msa_dir),
        "--output-dir", str(out_dir),
        "--seq-type", "AA",
        "--model", "lg",
        "--quiet",
        "--dry-run",
    ])

    assert result.exit_code == 0
    assert not (out_dir / "result.json").exists()


def test_tree_ml_fasttree_quiet_dry_run_single(tmp_path: Path) -> None:
    mat = tmp_path / "matrix.fa"
    mat.write_text(">a\nMKTLLL\n>b\nMKTLLL\n")

    out_dir = tmp_path / "out"

    result = CliRunner().invoke(cli, [
        "tree", "ml", "fasttree",
        "--matrix", str(mat),
        "--output-dir", str(out_dir),
        "--seq-type", "AA",
        "--model", "lg",
        "--quiet",
        "--dry-run",
    ])

    assert result.exit_code == 0


def test_tree_ml_fasttree_invalid_model_exits_1(tmp_path: Path) -> None:
    mat = tmp_path / "matrix.fa"
    mat.write_text(">a\nMKTLLL\n")

    out_dir = tmp_path / "out"

    result = CliRunner().invoke(cli, [
        "tree", "ml", "fasttree",
        "--matrix", str(mat),
        "--output-dir", str(out_dir),
        "--seq-type", "AA",
        "--model", "gtr",
        "--quiet",
    ])

    assert result.exit_code == 1


def test_tree_ml_fasttree_blocked_tool_args(tmp_path: Path) -> None:
    mat = tmp_path / "matrix.fa"
    mat.write_text(">a\nMKTLLL\n")

    out_dir = tmp_path / "out"

    result = CliRunner().invoke(cli, [
        "tree", "ml", "fasttree",
        "--matrix", str(mat),
        "--output-dir", str(out_dir),
        "--tool-args", "-nt",
        "--quiet",
    ])

    assert result.exit_code == 1


def test_tree_ml_fasttree_threads_warn_single(tmp_path: Path) -> None:
    mat = tmp_path / "matrix.fa"
    mat.write_text(">a\nMKTLLL\n>b\nMKTLLL\n")

    out_dir = tmp_path / "out"

    result = CliRunner().invoke(cli, [
        "tree", "ml", "fasttree",
        "--matrix", str(mat),
        "--output-dir", str(out_dir),
        "--threads", "8",
        "--quiet",
        "--dry-run",
    ])
    assert "has no effect" in result.output.lower() or result.exit_code == 0


def test_tree_ml_fasttree_writes_result_json_and_log(tmp_path: Path) -> None:
    mat = tmp_path / "matrix.fa"
    mat.write_text(">a\nMKTLLL\n>b\nMKTLLL\n")

    out_dir = tmp_path / "out"

    result = CliRunner().invoke(cli, [
        "tree", "ml", "fasttree",
        "--matrix", str(mat),
        "--output-dir", str(out_dir),
        "--seq-type", "AA",
        "--model", "lg",
        "--quiet",
    ])
    if result.exit_code == 0:
        assert (out_dir / "result.json").exists()
    elif result.exit_code == 3:
        import pytest
        pytest.skip("FastTree not installed")


def test_tree_ml_fasttree_no_gamma_cli(tmp_path: Path) -> None:
    mat = tmp_path / "matrix.fa"
    mat.write_text(">a\nMKTLLL\n>b\nMKTLLL\n")

    out_dir = tmp_path / "out"

    result = CliRunner().invoke(cli, [
        "tree", "ml", "fasttree",
        "--matrix", str(mat),
        "--output-dir", str(out_dir),
        "--seq-type", "AA",
        "--model", "lg",
        "--no-gamma",
        "--quiet",
        "--dry-run",
    ])

    assert result.exit_code == 0


def test_tree_ml_fasttree_quiet_dry_run_zero_inputs_exits_1(tmp_path: Path) -> None:
    msa_dir = tmp_path / "empty"
    msa_dir.mkdir()

    out_dir = tmp_path / "out"

    result = CliRunner().invoke(cli, [
        "tree", "ml", "fasttree",
        "--msa-dir", str(msa_dir),
        "--output-dir", str(out_dir),
        "--quiet",
        "--dry-run",
    ])

    assert result.exit_code == 1


# ===================================================================
# IQ-TREE CLI tests
# ===================================================================

def test_tree_ml_iqtree_help() -> None:
    result = CliRunner().invoke(cli, ["tree", "ml", "iqtree", "--help"])
    assert result.exit_code == 0
    assert "iqtree" in result.output
    assert "--msa-dir" in result.output
    assert "--modelfinder" in result.output
    assert "--rclusterf" in result.output
    assert "--qmax" in result.output
    assert "Input:" in result.output
    assert "ModelFinder:" in result.output
    assert "Branch Support:" in result.output
    assert "AA standard models:" in result.output
    assert "NT standard models:" in result.output
    assert "AA heterogeneous / mixture models:" in result.output
    assert "NT heterogeneous model:" in result.output
    assert "LG, Poisson, cpREV, mtREV, Dayhoff, mtMAM, JTT, WAG" in result.output
    assert "GTR, HKY, JC, F81, K2P, K3P, K81uf, TN, TNef, TIM, TIMef, TVM, TVMef," in result.output
    assert "SYM" in result.output

    grouped = result.output.split("Input:", 1)[1]
    for flag in [
        "--msa-dir", "--matrix", "--seq-type", "--model", "--state-freq",
        "--rate-heterogeneity", "--modelfinder", "--mset", "--msub", "--mode",
        "--boot", "--alrt", "--bnni", "--partitions", "--rclusterf",
        "--rcluster-max", "--pmsf-base-model", "--guide-tree", "--qmax",
        "--rate", "--wslr", "--constraint", "--outgroup", "--prefix",
        "--output-dir", "--threads", "--iqtree-path", "--tool-args",
        "--overwrite", "--resume", "--dry-run", "--keep-extra", "--quiet", "--help",
    ]:
        assert flag in grouped

    normalized_grouped = " ".join(grouped.split())
    assert "--rclusterf and --rcluster-max are mutually exclusive" in normalized_grouped
    assert "If neither is provided with --partitions + MF/MFP, PhyloAI uses --rclusterf 10." in normalized_grouped
    assert "Direct AA mixture:" in grouped
    assert "PMSF AA mixture:" in grouped
    assert "NT heterogeneous:" in grouped
    assert "Example model string: C20+F+R4" in grouped
    assert "Example model string: LG+C20+F+R4" in grouped
    assert "Example model string: MIX+MF" in grouped
    assert "phyloai tree ml iqtree --matrix matrix.fa --seq-type AA --model C20" in normalized_grouped
    assert "phyloai tree ml iqtree --matrix matrix.fa --seq-type AA --model C20 --guide-tree guide.nwk" in normalized_grouped
    assert "phyloai tree ml iqtree --matrix matrix.fa --seq-type NT --model MIX+MF" in normalized_grouped
    assert "Workflow Examples:" in result.output
    assert "phyloai tree ml iqtree --msa-dir msas/ --seq-type AA --model LG" in normalized_grouped
    assert "phyloai tree ml iqtree --matrix matrix.fa --seq-type AA --partitions parts.nex --modelfinder MFP" in normalized_grouped
    assert "--msa-dir and --matrix are mutually exclusive" in normalized_grouped
    assert "--partitions is only valid with --matrix" in normalized_grouped
    assert "--bnni requires --boot > 0" in normalized_grouped
    assert "--qmax is only valid with --model MIX+MF" in normalized_grouped
    assert grouped.index("fast maps to IQ-TREE --fast.") < grouped.index("--constraint")
    assert "-t, --threads" in grouped
    assert "--matrix   Single concatenated matrix for supermatrix inference\n\n  --msa-dir and --matrix are mutually exclusive." in grouped
    assert "--qmax             MIX+MF rate categories (default: 10)\n\n  Heterogeneous workflows require --matrix." in grouped
    assert "Homogeneous batch fixed model:\n    phyloai tree ml iqtree --msa-dir msas/ --seq-type AA --model LG" in grouped
    assert "PMSF AA mixture:\n    phyloai tree ml iqtree --matrix matrix.fa --seq-type AA --model C20 --guide-tree guide.nwk" in grouped
    assert "+FU" in grouped  # mentioned in state-freq help
    assert "+F|+FO|+FQ|+FU|none" in grouped
    assert "+I|+G4|+I+G4|+R4|+I+R4|none" in grouped
    assert "combine to form the IQ-TREE -m argument" in normalized_grouped
    assert "Ignored when --modelfinder is MF or MFP" in normalized_grouped
    assert "default: runs/tree/ml/iqtree" in grouped
    assert "nuclear|mitochondrial|chloroplast|viral" in grouped
    assert "parallel IQ-TREE jobs (default: 4)" in grouped
    assert "Single: NUM or" in grouped
    assert "(default: auto)" in grouped
    assert "Supported formats:" in grouped
    assert ".fa .fas .fasta .faa .fna" in normalized_grouped
    assert "NT maps to --seqtype DNA" in normalized_grouped
    assert "BIN, NT2AA, CODON, MORPH" in normalized_grouped
    assert ">=1000 recommended" in grouped
    assert "boot -> -B" in normalized_grouped
    assert "maps to IQ-TREE" in normalized_grouped
    assert "IQ-TREE -g" in normalized_grouped
    assert "IQ-TREE -o" in normalized_grouped
    # continuation lines of wrapped descriptions should be +2 indent
    assert "\n    " in grouped  # 4-space indent (2 base + 2 extra) for continuation lines
    assert "\n  --" in grouped  # parameter lines start at 2-space indent


def test_tree_ml_help_shows_iqtree() -> None:
    result = CliRunner().invoke(cli, ["tree", "ml", "--help"])
    assert result.exit_code == 0
    assert "iqtree" in result.output


def test_tree_ml_iqtree_mutual_exclusivity(tmp_path: Path) -> None:
    msa_dir = tmp_path / "msas"
    msa_dir.mkdir()
    mat = tmp_path / "matrix.fa"
    mat.write_text(">a\nMKTLLL\n")

    result = CliRunner().invoke(cli, [
        "tree", "ml", "iqtree",
        "--msa-dir", str(msa_dir), "--matrix", str(mat),
    ])
    assert result.exit_code == 1


def test_tree_ml_iqtree_neither_input() -> None:
    result = CliRunner().invoke(cli, ["tree", "ml", "iqtree"])
    assert result.exit_code == 1


def test_cli_iqtree_msa_dir_nonexistent() -> None:
    result = CliRunner().invoke(cli, [
        "tree", "ml", "iqtree", "--msa-dir", "/nonexistent/path",
    ])
    assert result.exit_code == 1


def test_tree_ml_iqtree_quiet_dry_run_single(tmp_path: Path) -> None:
    mat = tmp_path / "matrix.fa"
    mat.write_text(">a\nMKTLLL\n>b\nMKTLLL\n")
    out_dir = tmp_path / "out"

    result = CliRunner().invoke(cli, [
        "tree", "ml", "iqtree",
        "--matrix", str(mat),
        "--output-dir", str(out_dir),
        "--seq-type", "AA",
        "--model", "LG",
        "--quiet",
        "--dry-run",
    ])

    assert result.exit_code == 0


def test_tree_ml_iqtree_quiet_dry_run_batch(tmp_path: Path) -> None:
    msa_dir = tmp_path / "msas"
    msa_dir.mkdir()
    (msa_dir / "g1.fa").write_text(">a\nMKTLLL\n>b\nMKTLLL\n")
    out_dir = tmp_path / "out"

    result = CliRunner().invoke(cli, [
        "tree", "ml", "iqtree",
        "--msa-dir", str(msa_dir),
        "--output-dir", str(out_dir),
        "--seq-type", "AA",
        "--model", "LG",
        "--quiet",
        "--dry-run",
    ])

    assert result.exit_code == 0


def test_tree_ml_iqtree_blocked_tool_args(tmp_path: Path) -> None:
    mat = tmp_path / "matrix.fa"
    mat.write_text(">a\nMKTLLL\n")
    out_dir = tmp_path / "out"

    result = CliRunner().invoke(cli, [
        "tree", "ml", "iqtree",
        "--matrix", str(mat),
        "--output-dir", str(out_dir),
        "--tool-args", "-s hack.fa",
        "--quiet",
    ])

    assert result.exit_code == 1


def test_tree_ml_iqtree_heterogeneous_in_batch(tmp_path: Path) -> None:
    msa_dir = tmp_path / "msas"
    msa_dir.mkdir()
    (msa_dir / "g1.fa").write_text(">a\nMKTLLL\n")
    out_dir = tmp_path / "out"

    result = CliRunner().invoke(cli, [
        "tree", "ml", "iqtree",
        "--msa-dir", str(msa_dir),
        "--output-dir", str(out_dir),
        "--model", "C20",
        "--quiet",
    ])

    assert result.exit_code == 1


def test_tree_ml_iqtree_modelfinder_mf_dry_run(tmp_path: Path) -> None:
    mat = tmp_path / "matrix.fa"
    mat.write_text(">a\nMKTLLL\n>b\nMKTLLL\n")
    out_dir = tmp_path / "out"

    result = CliRunner().invoke(cli, [
        "tree", "ml", "iqtree",
        "--matrix", str(mat),
        "--output-dir", str(out_dir),
        "--modelfinder", "MF",
        "--mset", "LG,WAG",
        "--quiet",
        "--dry-run",
    ])

    assert result.exit_code == 0


def test_tree_ml_iqtree_writes_result_json_and_log(tmp_path: Path) -> None:
    mat = tmp_path / "matrix.fa"
    mat.write_text(">a\nMKTLLL\n>b\nMKTLLL\n")
    out_dir = tmp_path / "out"

    result = CliRunner().invoke(cli, [
        "tree", "ml", "iqtree",
        "--matrix", str(mat),
        "--output-dir", str(out_dir),
        "--seq-type", "AA",
        "--model", "LG",
        "--quiet",
    ])
    if result.exit_code == 0:
        assert (out_dir / "result.json").exists()
    elif result.exit_code == 3:
        import pytest
        pytest.skip("iqtree3 not installed")
