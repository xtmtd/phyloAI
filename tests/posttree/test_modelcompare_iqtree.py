"""Tests for phyloai.posttree.modelcompare_iqtree."""
from __future__ import annotations

from pathlib import Path

import pytest

from phyloai.posttree.modelcompare_iqtree import (
    _build_cmd,
    _expand_heterogeneous_models,
    _parse_modelfinder_results,
    _validate_inputs,
    run_modelcompare_iqtree,
)

FIXTURE_DIR = Path("runs/modelCompare")


class TestExpandHeterogeneousModels:
    def test_single_model_g(self) -> None:
        assert _expand_heterogeneous_models(["C10"], "G") == [
            "C10+G4", "C10+F+G4",
        ]

    def test_single_model_e(self) -> None:
        assert _expand_heterogeneous_models(["C10"], "E") == [
            "C10", "C10+F",
        ]

    def test_single_model_e_g(self) -> None:
        assert _expand_heterogeneous_models(["C10"], "E,G") == [
            "C10", "C10+F", "C10+G4", "C10+F+G4",
        ]

    def test_single_model_g_r(self) -> None:
        assert _expand_heterogeneous_models(["C10"], "G,R") == [
            "C10+G4", "C10+F+G4", "C10+R4", "C10+F+R4",
        ]

    def test_two_models_e_g_r(self) -> None:
        assert len(_expand_heterogeneous_models(["C10", "C20"], "E,G,R")) == 12

    def test_f_then_rate_order(self) -> None:
        models = _expand_heterogeneous_models(["C60"], "E,G,R")
        assert models[0] == "C60"
        assert models[1] == "C60+F"
        assert models[2] == "C60+G4"
        assert models[3] == "C60+F+G4"
        assert models[4] == "C60+R4"
        assert models[5] == "C60+F+R4"


class TestValidateInputs:
    def test_missing_matrix(self, tmp_path: Path) -> None:
        errors = _validate_inputs(
            matrix=tmp_path / "nope.fa",
            homogeneous_models=["LG"],
            mrate="E,G",
            heterogeneous_models=None,
            het_mrate="G",
            seq_type="AA",
            threads="auto",
            tool_args=None,
            overwrite=False,
            resume=False,
        )
        assert any("does not exist" in e for e in errors)

    def test_rejects_nt_with_heterogeneous(self, tmp_path: Path) -> None:
        matrix = tmp_path / "m.fa"
        matrix.write_text(">a\nACGT\n")
        errors = _validate_inputs(
            matrix=matrix,
            homogeneous_models=["GTR"],
            mrate="E,G",
            heterogeneous_models=["C10"],
            het_mrate="G",
            seq_type="NT",
            threads="auto",
            tool_args=None,
            overwrite=False,
            resume=False,
        )
        assert any("amino-acid" in e for e in errors)

    def test_rejects_invalid_mrate(self, tmp_path: Path) -> None:
        matrix = tmp_path / "m.fa"
        matrix.write_text(">a\nACGT\n")
        errors = _validate_inputs(
            matrix=matrix,
            homogeneous_models=["GTR"],
            mrate="E,X",
            heterogeneous_models=None,
            het_mrate="G",
            seq_type="NT",
            threads="auto",
            tool_args=None,
            overwrite=False,
            resume=False,
        )
        assert any("invalid token" in e for e in errors)

    def test_rejects_bad_homogeneous_model(self, tmp_path: Path) -> None:
        matrix = tmp_path / "m.fa"
        matrix.write_text(">a\nACGT\n")
        errors = _validate_inputs(
            matrix=matrix,
            homogeneous_models=["LG"],  # LG is AA-only
            mrate="E,G",
            heterogeneous_models=None,
            het_mrate="G",
            seq_type="NT",
            threads="auto",
            tool_args=None,
            overwrite=False,
            resume=False,
        )
        assert any("not valid for NT" in e for e in errors)

    def test_overwrite_resume_mutex(self, tmp_path: Path) -> None:
        matrix = tmp_path / "m.fa"
        matrix.write_text(">a\nACGT\n")
        errors = _validate_inputs(
            matrix=matrix,
            homogeneous_models=["GTR"],
            mrate="E,G",
            heterogeneous_models=None,
            het_mrate="G",
            seq_type="NT",
            threads="auto",
            tool_args=None,
            overwrite=True,
            resume=True,
        )
        assert any("mutually exclusive" in e for e in errors)

    def test_rejects_invalid_threads(self, tmp_path: Path) -> None:
        matrix = tmp_path / "m.fa"
        matrix.write_text(">a\nACGT\n")
        for bad in ("abc", "0", "-1"):
            errors = _validate_inputs(
                matrix=matrix,
                homogeneous_models=["GTR"],
                mrate="E,G",
                heterogeneous_models=None,
                het_mrate="G",
                seq_type="NT",
                threads=bad,
                tool_args=None,
                overwrite=False,
                resume=False,
            )
            assert any("--threads" in e for e in errors), bad

    def test_accepts_valid_threads(self, tmp_path: Path) -> None:
        matrix = tmp_path / "m.fa"
        matrix.write_text(">a\nACGT\n")
        for good in ("auto", "1", "8"):
            errors = _validate_inputs(
                matrix=matrix,
                homogeneous_models=["GTR"],
                mrate="E,G",
                heterogeneous_models=None,
                het_mrate="G",
                seq_type="NT",
                threads=good,
                tool_args=None,
                overwrite=False,
                resume=False,
            )
            assert not any("--threads" in e for e in errors), good

    def test_rejects_blocked_tool_args(self, tmp_path: Path) -> None:
        matrix = tmp_path / "m.fa"
        matrix.write_text(">a\nACGT\n")
        errors = _validate_inputs(
            matrix=matrix,
            homogeneous_models=["GTR"],
            mrate="E,G",
            heterogeneous_models=None,
            het_mrate="G",
            seq_type="NT",
            threads="auto",
            tool_args="-s other.fa",
            overwrite=False,
            resume=False,
        )
        assert any("blocked flag" in e for e in errors)

    def test_rejects_equals_form_blocked_tool_args(self, tmp_path: Path) -> None:
        matrix = tmp_path / "m.fa"
        matrix.write_text(">a\nACGT\n")
        for tool_args in ("--prefix=other", "-s=other.fa"):
            errors = _validate_inputs(
                matrix=matrix,
                homogeneous_models=["GTR"],
                mrate="E,G",
                heterogeneous_models=None,
                het_mrate="G",
                seq_type="NT",
                threads="auto",
                tool_args=tool_args,
                overwrite=False,
                resume=False,
            )
            assert any("blocked flag" in e for e in errors), tool_args

    def test_allows_nonblocked_tool_args(self, tmp_path: Path) -> None:
        matrix = tmp_path / "m.fa"
        matrix.write_text(">a\nACGT\n")
        errors = _validate_inputs(
            matrix=matrix,
            homogeneous_models=["GTR"],
            mrate="E,G",
            heterogeneous_models=None,
            het_mrate="G",
            seq_type="NT",
            threads="auto",
            tool_args="-cmin 2 -cmax 8 -pers 0.5",
            overwrite=False,
            resume=False,
        )
        assert not any("blocked flag" in e for e in errors)

    def test_e_het_mrate_allowed(self, tmp_path: Path) -> None:
        matrix = tmp_path / "m.fa"
        matrix.write_text(">a\nMKTLLL\n")
        errors = _validate_inputs(
            matrix=matrix,
            homogeneous_models=["LG"],
            mrate="E,G",
            heterogeneous_models=["C10"],
            het_mrate="E",
            seq_type="AA",
            threads="auto",
            tool_args=None,
            overwrite=False,
            resume=False,
        )
        assert not errors

    def test_invalid_het_mrate_token_rejected(self, tmp_path: Path) -> None:
        matrix = tmp_path / "m.fa"
        matrix.write_text(">a\nMKTLLL\n")
        errors = _validate_inputs(
            matrix=matrix,
            homogeneous_models=["LG"],
            mrate="E,G",
            heterogeneous_models=["C10"],
            het_mrate="X",
            seq_type="AA",
            threads="auto",
            tool_args=None,
            overwrite=False,
            resume=False,
        )
        assert any("--het-mrate" in e for e in errors)

    def test_empty_het_mrate_rejected_when_heterogeneous(self, tmp_path: Path) -> None:
        matrix = tmp_path / "m.fa"
        matrix.write_text(">a\nMKTLLL\n")
        errors = _validate_inputs(
            matrix=matrix,
            homogeneous_models=["LG"],
            mrate="E,G",
            heterogeneous_models=["C10"],
            het_mrate="",
            seq_type="AA",
            threads="auto",
            tool_args=None,
            overwrite=False,
            resume=False,
        )
        assert any("at least one of E, G, R" in e for e in errors)

    def test_empty_het_mrate_ignored_without_heterogeneous(self, tmp_path: Path) -> None:
        matrix = tmp_path / "m.fa"
        matrix.write_text(">a\nMKTLLL\n")
        errors = _validate_inputs(
            matrix=matrix,
            homogeneous_models=["LG"],
            mrate="E,G",
            heterogeneous_models=None,
            het_mrate="",
            seq_type="AA",
            threads="auto",
            tool_args=None,
            overwrite=False,
            resume=False,
        )
        assert not errors

    def test_unsafe_prefix_rejected(self, tmp_path: Path) -> None:
        matrix = tmp_path / "m.fa"
        matrix.write_text(">a\nMKTLLL\n")
        for bad in ("../out", "a/b", "..", ".", "/abs"):
            errors = _validate_inputs(
                matrix=matrix,
                homogeneous_models=["LG"],
                mrate="E,G",
                heterogeneous_models=None,
                het_mrate="E,G",
                seq_type="AA",
                threads="auto",
                prefix=bad,
                tool_args=None,
                overwrite=False,
                resume=False,
            )
            assert any("--prefix" in e for e in errors), bad

    def test_safe_prefix_allowed(self, tmp_path: Path) -> None:
        matrix = tmp_path / "m.fa"
        matrix.write_text(">a\nMKTLLL\n")
        errors = _validate_inputs(
            matrix=matrix,
            homogeneous_models=["LG"],
            mrate="E,G",
            heterogeneous_models=None,
            het_mrate="E,G",
            seq_type="AA",
            threads="auto",
            prefix="modelcompare",
            tool_args=None,
            overwrite=False,
            resume=False,
        )
        assert not errors


class TestBuildCmd:
    def test_base_command(self) -> None:
        cmd = _build_cmd(
            executable="iqtree3",
            matrix=Path("m.fa"),
            homogeneous_models=["LG", "WAG"],
            mrate="E,G",
            expanded_het=None,
            prefix="modelcompare",
            threads="auto",
            tool_args=None,
        )
        assert cmd == [
            "iqtree3", "-s", "m.fa", "-m", "MF", "-mset", "LG,WAG",
            "-mrate", "E,G", "-cmin", "4", "-cmax", "4",
            "--prefix", "modelcompare", "-T", "auto",
        ]

    def test_with_het_and_tool_args(self) -> None:
        cmd = _build_cmd(
            executable="iqtree3",
            matrix=Path("m.fa"),
            homogeneous_models=["LG"],
            mrate="E,G",
            expanded_het=["C10", "C10+F"],
            prefix="mc",
            threads="4",
            tool_args="-n 0",
        )
        assert "-madd" in cmd
        assert cmd[cmd.index("-madd") + 1] == "C10,C10+F"
        assert cmd[-2:] == ["-n", "0"]

    def test_tool_args_override_managed_flags(self) -> None:
        cmd = _build_cmd(
            executable="iqtree3",
            matrix=Path("m.fa"),
            homogeneous_models=["LG"],
            mrate="E,G",
            expanded_het=["C10", "C10+F"],
            prefix="mc",
            threads="4",
            tool_args="-madd C20 -T 8 -mrate G",
        )
        assert cmd.count("-madd") == 1
        assert cmd[cmd.index("-madd") + 1] == "C20"
        assert cmd.count("-T") == 1
        assert cmd[cmd.index("-T") + 1] == "8"
        assert cmd.count("-mrate") == 1
        assert cmd[cmd.index("-mrate") + 1] == "G"

    def test_tool_args_equals_form_override(self) -> None:
        cmd = _build_cmd(
            executable="iqtree3",
            matrix=Path("m.fa"),
            homogeneous_models=["LG"],
            mrate="E,G",
            expanded_het=["C10", "C10+F"],
            prefix="mc",
            threads="auto",
            tool_args="-T=8",
        )
        assert "-T=8" in cmd
        assert "-T" not in cmd or cmd.count("-T=8") == 1


class TestParseModelFinderResults:
    def test_fixture_parses_12_models(self) -> None:
        iqtree_file = FIXTURE_DIR / "EOG090X0A0V.fa.iqtree"
        if not iqtree_file.exists():
            pytest.skip("ModelCompare IQ-TREE fixture not present")
        models = _parse_modelfinder_results(iqtree_file)
        assert len(models) == 12
        best = models[0]
        assert best["model"] == "LG+G4"
        assert best["bic"] < models[1]["bic"]
        assert best["in_bic_95"] is True
        assert "w_bic" in best
        assert "in_aic_95" in best
        assert "in_aicc_95" in best

    def test_missing_section_raises(self, tmp_path: Path) -> None:
        bad = tmp_path / "no_models.iqtree"
        bad.write_text("IQ-TREE header\nNo modelfinder section here.\n")
        with pytest.raises(ValueError, match="List of models sorted by BIC"):
            _parse_modelfinder_results(bad)


class TestRunModelcompareIqtree:
    def test_dry_run(self, tmp_path: Path) -> None:
        matrix = FIXTURE_DIR / "EOG090X0A0V.fa"
        if not matrix.exists():
            pytest.skip("ModelCompare fixture not present")
        result = run_modelcompare_iqtree(
            matrix=matrix,
            homogeneous_model="LG",
            mrate="E,G,R",
            heterogeneous_model="C10",
            het_mrate="G",
            seq_type="auto",
            output_dir=tmp_path / "mc_out",
            dry_run=True,
        )
        assert result["status"] == "success"
        assert result["params"]["detected_seq_type"] == "AA"
        assert result["data"]["cmd"][0:3] == ["iqtree3", "-s", str(matrix.resolve())]
        assert "-m" in result["data"]["cmd"]
        assert result["command"].startswith("phyloai posttree modelcompare iqtree")

    def test_params_het_mrate_none_when_no_heterogeneous(self, tmp_path: Path) -> None:
        matrix = FIXTURE_DIR / "EOG090X0A0V.fa"
        if not matrix.exists():
            pytest.skip("ModelCompare fixture not present")
        result = run_modelcompare_iqtree(
            matrix=matrix,
            homogeneous_model="LG",
            mrate="E,G",
            seq_type="auto",
            output_dir=tmp_path / "mc_hom",
            dry_run=True,
        )
        assert result["status"] == "success"
        assert result["params"]["heterogeneous_model"] is None
        assert result["params"]["het_mrate"] is None

    def test_params_het_mrate_linked_when_heterogeneous(self, tmp_path: Path) -> None:
        matrix = FIXTURE_DIR / "EOG090X0A0V.fa"
        if not matrix.exists():
            pytest.skip("ModelCompare fixture not present")
        result = run_modelcompare_iqtree(
            matrix=matrix,
            homogeneous_model="LG",
            mrate="E,G",
            heterogeneous_model="C10",
            het_mrate="E",
            seq_type="auto",
            output_dir=tmp_path / "mc_het",
            dry_run=True,
        )
        assert result["status"] == "success"
        assert result["params"]["heterogeneous_model"] == "C10"
        assert result["params"]["het_mrate"] == "E"
        assert "-madd" in result["data"]["cmd"]
        assert result["data"]["cmd"][result["data"]["cmd"].index("-madd") + 1] == "C10,C10+F"

    def test_missing_matrix_error(self, tmp_path: Path) -> None:
        result = run_modelcompare_iqtree(
            matrix=tmp_path / "nope.fa",
            homogeneous_model="LG",
            mrate="E,G",
            output_dir=tmp_path / "mc_out",
        )
        assert result["status"] == "error"
        assert result["error_category"] == "input"

    def test_seq_type_conflict_with_matrix(self, tmp_path: Path) -> None:
        matrix = FIXTURE_DIR / "EOG090X0A0V.fa"  # AA
        if not matrix.exists():
            pytest.skip("ModelCompare fixture not present")
        result = run_modelcompare_iqtree(
            matrix=matrix,
            homogeneous_model="GTR",
            mrate="E,G",
            seq_type="NT",  # wrong: matrix is AA
            output_dir=tmp_path / "mc_nt",
            dry_run=True,
        )
        assert result["status"] == "error"
        assert result["error_category"] == "input"
        assert "conflicts" in result["error"]

    def test_seq_type_match_passes(self, tmp_path: Path) -> None:
        matrix = FIXTURE_DIR / "EOG090X0A0V.fa"  # AA
        if not matrix.exists():
            pytest.skip("ModelCompare fixture not present")
        result = run_modelcompare_iqtree(
            matrix=matrix,
            homogeneous_model="LG",
            mrate="E,G",
            seq_type="AA",
            output_dir=tmp_path / "mc_aa",
            dry_run=True,
        )
        assert result["status"] == "success"
        assert result["params"]["detected_seq_type"] == "AA"
