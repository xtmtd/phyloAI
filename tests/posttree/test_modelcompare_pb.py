"""Tests for phyloai.posttree.modelcompare_pb."""
from __future__ import annotations

import csv as csv_mod
from pathlib import Path

import pytest

from phyloai.posttree.modelcompare_pb import (
    _classify_quality,
    _compute_loocv_waic,
    _parse_sitelogl,
    run_modelcompare_pb,
)

LOOCV_DIR = Path("runs/modelCompare/LOOCV_wAIC")


class TestParseSitelogl:
    def test_fixture_235_rows(self) -> None:
        path = LOOCV_DIR / "chain1.sitelogl"
        if not path.exists():
            pytest.skip("LOOCV fixture not present")
        rows = _parse_sitelogl(path)
        assert len(rows) == 235
        assert len(rows[0]) == 7
        assert rows[0][0] == 1.0  # site id

    def test_missing_file_raises(self, tmp_path: Path) -> None:
        with pytest.raises(OSError):
            _parse_sitelogl(tmp_path / "nope.sitelogl")

    def test_malformed_row_raises(self, tmp_path: Path) -> None:
        path = tmp_path / "bad.sitelogl"
        path.write_text("site\tlogl\tvar\n1\t2\n")
        with pytest.raises(ValueError, match="expected 7 columns"):
            _parse_sitelogl(path)

    def test_empty_file_raises(self, tmp_path: Path) -> None:
        path = tmp_path / "empty.sitelogl"
        path.write_text("site\tlogl\tvar\n")
        with pytest.raises(ValueError, match="no data rows"):
            _parse_sitelogl(path)

    def test_duplicate_site_id_raises(self, tmp_path: Path) -> None:
        path = tmp_path / "dup.sitelogl"
        path.write_text(
            "site\tlogl\tvar\tlogcpo\tess\tlogpostmeanl\tess\n"
            "1\t-1.0\t0.1\t-1.0\t60.0\t-1.0\t65.0\n"
            "1\t-1.5\t0.2\t-1.5\t60.0\t-1.5\t65.0\n"
        )
        with pytest.raises(ValueError, match="Duplicate site identifier"):
            _parse_sitelogl(path)


class TestClassifyQuality:
    def test_good(self) -> None:
        assert _classify_quality(0.05, 0.05) == "good"

    def test_ok(self) -> None:
        assert _classify_quality(0.15, 0.05) == "ok"

    def test_no(self) -> None:
        assert _classify_quality(0.35, 0.05) == "no"
        assert _classify_quality(0.05, 0.5) == "no"


class TestComputeLoocvWaic:
    def test_three_chain_fixture(self) -> None:
        files = [LOOCV_DIR / f"chain{i}.sitelogl" for i in (1, 2, 3)]
        if not files[0].exists():
            pytest.skip("LOOCV fixture not present")
        runs = [[_parse_sitelogl(p) for p in files]]
        stats = _compute_loocv_waic(runs[0])
        # Reference readwaic.py: loo=-11.5443, waic=-11.5330
        assert stats["loocv"]["score"] == pytest.approx(-11.5443, abs=1e-3)
        assert stats["waic"]["score"] == pytest.approx(-11.5330, abs=1e-3)
        assert stats["loocv"]["ess"] == pytest.approx(62.86, abs=1e-2)
        assert stats["waic"]["ess"] == pytest.approx(69.43, abs=1e-2)
        assert stats["loocv"]["quality"] == "good"
        assert stats["waic"]["quality"] == "good"


def _synthetic_sitelogl(tmp_path: Path, site_ids: list[int], logcpo: float,
                        logpostmeanl: float, var: float) -> Path:
    name = "_".join(f"s{s}" for s in site_ids)
    path = tmp_path / f"{name}.sitelogl"
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = ["site\tlogl\tvar\tlogcpo\tess\tlogpostmeanl\tess\n"]
    for s in site_ids:
        lines.append(
            f"{s}\t{logpostmeanl + var:.6f}\t{var:.6f}\t{logcpo:.6f}\t60.0\t{logpostmeanl:.6f}\t65.0\n"
        )
    path.write_text("".join(lines))
    return path


class TestRunModelcomparePb:
    def test_single_model_from_dir(self, tmp_path: Path) -> None:
        src = LOOCV_DIR / "model1"
        if not src.exists():
            pytest.skip("LOOCV fixture not present")
        result = run_modelcompare_pb(
            sitelogl_dir=str(src),
            model_names="model_1",
            output_dir=tmp_path / "pb_out",
        )
        assert result["status"] == "success"
        kr = result["key_results"]
        assert kr["n_models"] == 1
        assert kr["n_sites"] == 235
        assert kr["n_runs"] == 2
        assert result["command"].startswith("phyloai posttree modelcompare pb")
        assert (tmp_path / "pb_out" / "result.json").exists()
        assert (tmp_path / "pb_out" / "model_fit.csv").exists()
        assert (tmp_path / "pb_out" / "sitelogl" / "model_1" / "chain1.sitelogl").exists()

    def test_two_models_delta(self, tmp_path: Path) -> None:
        src = LOOCV_DIR
        if not (src / "model1").exists() or not (src / "model2").exists():
            pytest.skip("LOOCV fixture not present")
        result = run_modelcompare_pb(
            sitelogl_dir=f"{src / 'model1'},{src / 'model2'}",
            model_names="model_1,model_2",
            output_dir=tmp_path / "pb_multi",
        )
        assert result["status"] == "success"
        kr = result["key_results"]
        assert kr["n_models"] == 2
        assert kr["n_sites"] == 235
        assert "best_loocv_quality" in kr
        assert "best_waic_quality" in kr
        assert kr["best_waic_quality"] in ("good", "ok", "no")
        models = {m["model"]: m for m in result["data"]["models"]}
        assert "model_1" in models and "model_2" in models
        # The best model has delta 0; the other is <= 0
        best_loocv = models[kr["best_model_loocv"]]
        assert best_loocv["delta_loocv"] == pytest.approx(0.0)
        other_loocv = models["model_2" if kr["best_model_loocv"] == "model_1" else "model_1"]
        assert other_loocv["delta_loocv"] <= 0.0
        assert best_loocv["delta_waic"] == pytest.approx(0.0)

    def test_synthetic_two_model_delta_sign(self, tmp_path: Path) -> None:
        # model_1 scores better on both metrics; model_2 must have delta < 0
        m1a = _synthetic_sitelogl(tmp_path / "m1", [1, 2], -10.0, -10.0, 0.5)
        m1b = _synthetic_sitelogl(tmp_path / "m1", [1, 2], -10.1, -10.1, 0.6)
        m2a = _synthetic_sitelogl(tmp_path / "m2", [1, 2], -15.0, -15.0, 0.5)
        m2b = _synthetic_sitelogl(tmp_path / "m2", [1, 2], -15.1, -15.1, 0.6)
        result = run_modelcompare_pb(
            sitelogl=[
                f"{m1a},{m1b}",
                f"{m2a},{m2b}",
            ],
            model_names="model_1,model_2",
            output_dir=tmp_path / "pb_syn",
        )
        assert result["status"] == "success"
        assert result["key_results"]["best_model_loocv"] == "model_1"
        assert result["key_results"]["best_model_waic"] == "model_1"
        assert result["key_results"]["best_waic_quality"] in ("good", "ok", "no")
        models = {m["model"]: m for m in result["data"]["models"]}
        assert models["model_1"]["delta_loocv"] == pytest.approx(0.0)
        assert models["model_2"]["delta_loocv"] < 0.0
        assert models["model_2"]["delta_waic"] < 0.0

    def test_within_group_site_mismatch_hard_error(self, tmp_path: Path) -> None:
        # Same model group, same site count, but different site identifiers/order
        m1a = _synthetic_sitelogl(tmp_path / "m1", [1, 2], -10.0, -10.0, 0.5)
        m1b = _synthetic_sitelogl(tmp_path / "m1", [2, 1], -10.1, -10.1, 0.6)
        result = run_modelcompare_pb(
            sitelogl=[f"{m1a},{m1b}"],
            model_names="model_1",
            output_dir=tmp_path / "pb_within",
        )
        assert result["status"] == "error"
        assert result["error_category"] == "input"
        assert "site identifiers" in result["error"]

    def test_cross_model_site_id_mismatch_hard_error(self, tmp_path: Path) -> None:
        # model_2 has a different site identifier set => not comparable
        m1a = _synthetic_sitelogl(tmp_path / "m1", [1, 2], -10.0, -10.0, 0.5)
        m1b = _synthetic_sitelogl(tmp_path / "m1", [1, 2], -10.1, -10.1, 0.6)
        m2a = _synthetic_sitelogl(tmp_path / "m2", [1, 3], -15.0, -15.0, 0.5)
        m2b = _synthetic_sitelogl(tmp_path / "m2", [1, 3], -15.1, -15.1, 0.6)
        result = run_modelcompare_pb(
            sitelogl=[f"{m1a},{m1b}", f"{m2a},{m2b}"],
            model_names="model_1,model_2",
            output_dir=tmp_path / "pb_mismatch",
        )
        assert result["status"] == "error"
        assert result["error_category"] == "input"
        assert "site identifiers" in result["error"]

    def test_cross_model_site_count_mismatch_hard_error(self, tmp_path: Path) -> None:
        # model_2 has 3 sites vs model_1's 2 => count mismatch
        m1a = _synthetic_sitelogl(tmp_path / "m1", [1, 2], -10.0, -10.0, 0.5)
        m1b = _synthetic_sitelogl(tmp_path / "m1", [1, 2], -10.1, -10.1, 0.6)
        m2a = _synthetic_sitelogl(tmp_path / "m2", [1, 2, 3], -15.0, -15.0, 0.5)
        m2b = _synthetic_sitelogl(tmp_path / "m2", [1, 2, 3], -15.1, -15.1, 0.6)
        result = run_modelcompare_pb(
            sitelogl=[f"{m1a},{m1b}", f"{m2a},{m2b}"],
            model_names="model_1,model_2",
            output_dir=tmp_path / "pb_count",
        )
        assert result["status"] == "error"
        assert result["error_category"] == "input"
        assert "site count" in result["error"]

    def test_nonempty_output_dir_rejected(self, tmp_path: Path) -> None:
        m1a = _synthetic_sitelogl(tmp_path / "m1", [1, 2], -10.0, -10.0, 0.5)
        m1b = _synthetic_sitelogl(tmp_path / "m1", [1, 2], -10.1, -10.1, 0.6)
        out = tmp_path / "pb_occupied"
        out.mkdir()
        (out / "stale.txt").write_text("old")
        result = run_modelcompare_pb(
            sitelogl=[f"{m1a},{m1b}"],
            model_names="model_1",
            output_dir=out,
        )
        assert result["status"] == "error"
        assert result["error_category"] == "input"
        assert "--overwrite" in result["error"]

    def test_overwrite_replaces_nonempty_output_dir(self, tmp_path: Path) -> None:
        m1a = _synthetic_sitelogl(tmp_path / "m1", [1, 2], -10.0, -10.0, 0.5)
        m1b = _synthetic_sitelogl(tmp_path / "m1", [1, 2], -10.1, -10.1, 0.6)
        out = tmp_path / "pb_overwrite"
        out.mkdir()
        (out / "stale.txt").write_text("old")
        result = run_modelcompare_pb(
            sitelogl=[f"{m1a},{m1b}"],
            model_names="model_1",
            output_dir=out,
            overwrite=True,
        )
        assert result["status"] == "success"
        assert not (out / "stale.txt").exists()
        assert (out / "result.json").exists()

    def test_both_input_modes_rejected(self, tmp_path: Path) -> None:
        result = run_modelcompare_pb(
            sitelogl_dir=str(LOOCV_DIR / "model1"),
            sitelogl=["a,b"],
            output_dir=tmp_path / "pb_both",
        )
        assert result["status"] == "error"
        assert result["error_category"] == "input"
        assert "Exactly one" in result["error"]

    def test_model_names_count_mismatch(self, tmp_path: Path) -> None:
        result = run_modelcompare_pb(
            sitelogl_dir=str(LOOCV_DIR / "model1"),
            model_names="a,b",
            output_dir=tmp_path / "pb_names",
        )
        assert result["status"] == "error"
        assert "label" in result["error"]

    def test_model_names_unsafe_rejected(self, tmp_path: Path) -> None:
        src = LOOCV_DIR / "model1"
        if not src.exists():
            pytest.skip("LOOCV fixture not present")
        for bad in ("../../target", "a/b", "..", ".", "/abs"):
            result = run_modelcompare_pb(
                sitelogl_dir=str(src),
                model_names=bad,
                output_dir=tmp_path / "pb_unsafe",
            )
            assert result["status"] == "error"
            assert "unsafe" in result["error"]
            assert not (tmp_path / "pb_unsafe").exists()

    def test_model_names_safe_copy_stays_inside(self, tmp_path: Path) -> None:
        src = LOOCV_DIR / "model1"
        if not src.exists():
            pytest.skip("LOOCV fixture not present")
        out = tmp_path / "pb_safe"
        result = run_modelcompare_pb(
            sitelogl_dir=str(src), model_names="CAT-GTR", output_dir=out,
        )
        assert result["status"] == "success"
        assert (out / "sitelogl" / "CAT-GTR").is_dir()

    def test_single_model_csv_columns(self, tmp_path: Path) -> None:
        src = LOOCV_DIR / "model1"
        if not src.exists():
            pytest.skip("LOOCV fixture not present")
        out = tmp_path / "pb_cols"
        run_modelcompare_pb(sitelogl_dir=str(src), output_dir=out)
        with open(out / "model_fit.csv") as handle:
            columns = csv_mod.DictReader(handle).fieldnames or []
        assert "Metric" in columns
        assert "Score" in columns
        assert "Quality" in columns
        rows = list(csv_mod.DictReader(open(out / "model_fit.csv")))
        assert {r["Metric"] for r in rows} == {"LOO-CV", "wAIC"}
