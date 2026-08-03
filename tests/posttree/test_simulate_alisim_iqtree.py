"""Tests for phyloai.posttree.simulate_alisim_iqtree."""
from __future__ import annotations

import random
from pathlib import Path

import pytest

from phyloai.posttree.simulate_alisim_iqtree import (
    build_model_string,
    run_alisim_iqtree,
    sample_batch_rows,
)


def _row(**overrides: str) -> dict[str, str]:
    base = {
        "id": "g1", "seqtype": "AA", "length": "2082", "subs_model": "LG",
        "subs_rate": "", "freq": "", "prop_inv": "", "rate_heterogeneity": "G",
        "rate_categories": "4", "rate_param": "0.6", "tree_path": "/t/g1.tre",
    }
    base.update(overrides)
    return base


def _rows() -> list[dict[str, str]]:
    return [
        _row(id="g1"),
        _row(id="g2", length="1000", prop_inv="0.1"),
        _row(id="g3", seqtype="DNA", subs_model="GTR",
             subs_rate="1/2/3/4/5", freq=".2/.3/.2/.3",
             rate_heterogeneity="G", rate_categories="4", rate_param="1.2",
             length="3000", tree_path="/t/g3.tre"),
    ]


def _fake_worker(args: tuple) -> dict:
    """Module-level stand-in for _run_simulation_worker (picklable).

    Writes a tiny FASTA so the batch completion path (move/wrap/validate)
    runs end to end without an IQ-TREE binary.
    """
    from pathlib import Path

    simulation_id = args[0]
    work_dir = Path(f"./_work_{simulation_id}")
    work_dir.mkdir(parents=True, exist_ok=True)
    (work_dir / f"{simulation_id}.fa").write_text(">taxa1\nMKT\n>taxa2\nMKA\n")
    return {
        "simulation_id": simulation_id,
        "status": "success",
        "wall_time": 0.01,
        "cmd": ["iqtree3", "--alisim", simulation_id],
        "log_file": f"logs/{simulation_id}.log",
        "log_text": "ok",
        "output_file": f"MSAs/{simulation_id}.fa",
        "output_path": work_dir / f"{simulation_id}.fa",
        "reason": None,
    }


class TestBuildModelString:
    def test_omits_absent_components(self) -> None:
        row = {"subs_model": "LG", "subs_rate": "", "freq": "", "prop_inv": "",
               "rate_heterogeneity": "R", "rate_categories": "2",
               "rate_param": "0.5/1.5"}
        assert build_model_string(row) == "LG+R2{0.5,1.5}"

    def test_full_dna_model(self) -> None:
        row = {"subs_model": "GTR", "subs_rate": "1/2/3/4/5",
               "freq": ".1/.2/.3/.4", "prop_inv": ".2",
               "rate_heterogeneity": "G", "rate_categories": "4",
               "rate_param": ".7"}
        assert build_model_string(row) == "GTR{1,2,3,4,5}+F{.1,.2,.3,.4}+I{.2}+G4{.7}"

    def test_retains_zero_invariable_proportion(self) -> None:
        row = {"subs_model": "LG", "subs_rate": "", "freq": "", "prop_inv": "0",
               "rate_heterogeneity": "", "rate_categories": "", "rate_param": ""}
        assert build_model_string(row) == "LG+I{0}"

    def test_bare_model(self) -> None:
        row = {"subs_model": "LG", "subs_rate": "", "freq": "", "prop_inv": "",
               "rate_heterogeneity": "", "rate_categories": "", "rate_param": ""}
        assert build_model_string(row) == "LG"


class TestSampleBatchRows:
    def test_complete_keeps_rows_intact(self) -> None:
        sampled = sample_batch_rows(
            _rows(), strategy="complete", n=10, rng=random.Random(1),
            pdf_params=(), noise_scale=1.0, overrides={},
        )
        assert len(sampled) == 10
        for row in sampled:
            assert row["subs_model"] in {"LG", "GTR"}
            assert row["tree_path"]

    def test_mixed_keeps_model_core_and_rate_group_together(self) -> None:
        sampled = sample_batch_rows(
            _rows(), strategy="mixed", n=20, rng=random.Random(3),
            pdf_params=(), noise_scale=1.0, overrides={},
        )
        for row in sampled:
            aa = row["seqtype"] == "AA" and row["subs_model"] == "LG" and row["subs_rate"] == ""
            dna = (row["seqtype"] == "DNA" and row["subs_model"] == "GTR"
                   and row["subs_rate"] == "1/2/3/4/5")
            assert aa or dna

    def test_mixed_prop_inv_presence_ratio_is_empirical(self) -> None:
        sampled = sample_batch_rows(
            _rows(), strategy="mixed", n=200, rng=random.Random(6),
            pdf_params=(), noise_scale=1.0, overrides={},
        )
        present = sum(1 for row in sampled if row["prop_inv"])
        assert present > 0
        assert present < 200

    def test_pdf_preserves_i_presence_before_resampling_value(self) -> None:
        sampled = sample_batch_rows(
            _rows(), strategy="pdf", n=100, rng=random.Random(4),
            pdf_params=("prop_inv",), noise_scale=0.0, overrides={},
        )
        assert any(row["prop_inv"] == "" for row in sampled)
        assert any(row["prop_inv"] for row in sampled)

    def test_pdf_resamples_rate_param_only_for_gamma(self) -> None:
        rows = [
            _row(id="a", rate_heterogeneity="G", rate_categories="4", rate_param="0.5"),
            _row(id="b", rate_heterogeneity="G", rate_categories="4", rate_param="1.5"),
            _row(id="c", rate_heterogeneity="R", rate_categories="2",
                 rate_param="0.5/1.5/0.5/1.5"),
            _row(id="d", rate_heterogeneity="", rate_categories="", rate_param=""),
        ]
        empirical_r = "0.5/1.5/0.5/1.5"
        sampled = sample_batch_rows(
            rows, strategy="pdf", n=100, rng=random.Random(5),
            pdf_params=("rate_param",), noise_scale=0.0, overrides={},
        )
        for row in sampled:
            if row["rate_heterogeneity"] == "R":
                assert row["rate_param"] == empirical_r
            elif row["rate_heterogeneity"] == "G":
                assert 0.5 <= float(row["rate_param"]) <= 1.5
            else:
                assert row["rate_param"] == ""

    def test_pdf_resamples_length(self) -> None:
        rows = [
            _row(id="a", length="100"),
            _row(id="b", length="200"),
            _row(id="c", length="300"),
        ]
        sampled = sample_batch_rows(
            rows, strategy="pdf", n=50, rng=random.Random(7),
            pdf_params=("length",), noise_scale=0.0, overrides={},
        )
        for row in sampled:
            assert int(row["length"]) >= 1

    def test_override_fixes_parameters(self) -> None:
        sampled = sample_batch_rows(
            _rows(), strategy="mixed", n=10, rng=random.Random(8),
            pdf_params=(), noise_scale=1.0, overrides={"length": "500", "prop_inv": "0.1"},
        )
        assert all(row["length"] == "500" for row in sampled)
        assert all(row["prop_inv"] == "0.1" for row in sampled)

    def test_sample_batch_rows_pdf_deterministic_seed(self) -> None:
        sampled = sample_batch_rows(
            _rows(), strategy="mixed", n=5, rng=random.Random(9),
            pdf_params=(), noise_scale=1.0, overrides={},
        )
        assert len(sampled) == 5


class TestRunAlisimIqtree:
    def test_batch_dry_run_records_deterministic_task_seeds(self, tmp_path: Path) -> None:
        table = tmp_path / "params.tsv"
        table.write_text(
            "id\tseqtype\tlength\tsubs_model\tsubs_rate\tfreq\tprop_inv\t"
            "rate_heterogeneity\trate_categories\trate_param\ttree_path\n"
            "g1\tAA\t2082\tLG\t\t\t\tG\t4\t0.6\t/t/g1.tre\n"
        )
        result = run_alisim_iqtree(
            model_params=table, strategy="complete", num_simulations=2,
            seed=7, output_dir=tmp_path / "out", dry_run=True,
        )
        assert [row["seed"] for row in result["data"]["sampled_rows"]] == [7, 8]

    def test_single_dry_run_uses_out_format_and_partition_p(
        self, tmp_path: Path,
    ) -> None:
        tree = tmp_path / "ref.tre"
        tree.write_text("(A,B);\n")
        parts = tmp_path / "p.nex"
        parts.write_text("charset c1 = 1-100;\n")
        result = run_alisim_iqtree(
            ref_tree=tree, model=None, model_partitions=parts, seq_type="AA",
            length=None, output_dir=tmp_path / "out", dry_run=True,
        )
        cmd = " ".join(result["data"]["cmd"])
        assert "-p" in cmd
        assert "--out-format" in cmd

    def test_single_dry_run_with_model_and_length(self, tmp_path: Path) -> None:
        tree = tmp_path / "ref.tre"
        tree.write_text("(A,B);\n")
        result = run_alisim_iqtree(
            ref_tree=tree, model="LG+G4{0.6}", seq_type="AA", length=100,
            output_dir=tmp_path / "out", dry_run=True,
        )
        cmd = result["data"]["cmd"]
        assert cmd[0] == "iqtree3"
        assert "--alisim" in cmd
        assert "-m" in cmd
        assert "--length" in cmd

    def test_batch_requires_model_params_options(self, tmp_path: Path) -> None:
        with pytest.raises(ValueError, match="strategy"):
            run_alisim_iqtree(
                model_params=tmp_path / "params.tsv", num_simulations=5,
                output_dir=tmp_path / "out", dry_run=True,
            )

    def test_noise_scale_requires_pdf(self, tmp_path: Path) -> None:
        with pytest.raises(ValueError, match="pdf"):
            run_alisim_iqtree(
                model_params=tmp_path / "params.tsv", strategy="complete",
                num_simulations=2, noise_scale=0.5,
                output_dir=tmp_path / "out", dry_run=True,
            )

    def test_invalid_override_key_rejected(self, tmp_path: Path) -> None:
        table = tmp_path / "params.tsv"
        table.write_text(
            "id\tseqtype\tlength\tsubs_model\tsubs_rate\tfreq\tprop_inv\t"
            "rate_heterogeneity\trate_categories\trate_param\ttree_path\n"
            "g1\tAA\t2082\tLG\t\t\t\tG\t4\t0.6\t/t/g1.tre\n"
        )
        with pytest.raises(ValueError, match="length|prop_inv"):
            run_alisim_iqtree(
                model_params=table, strategy="complete", num_simulations=2,
                override="seqtype=DNA", output_dir=tmp_path / "out", dry_run=True,
            )

    def test_single_and_batch_mutually_exclusive(self, tmp_path: Path) -> None:
        tree = tmp_path / "ref.tre"
        tree.write_text("(A,B);\n")
        table = tmp_path / "params.tsv"
        table.write_text(
            "id\tseqtype\tlength\tsubs_model\tsubs_rate\tfreq\tprop_inv\t"
            "rate_heterogeneity\trate_categories\trate_param\ttree_path\n"
            "g1\tAA\t2082\tLG\t\t\t\tG\t4\t0.6\t/t/g1.tre\n"
        )
        with pytest.raises(ValueError, match="model-params"):
            run_alisim_iqtree(
                ref_tree=tree, model="LG", seq_type="AA", length=100,
                model_params=table, strategy="complete", num_simulations=2,
                output_dir=tmp_path / "out", dry_run=True,
            )

    def test_blocked_tool_args_rejected(self, tmp_path: Path) -> None:
        tree = tmp_path / "ref.tre"
        tree.write_text("(A,B);\n")
        with pytest.raises(ValueError, match="Blocked"):
            run_alisim_iqtree(
                ref_tree=tree, model="LG", seq_type="AA", length=100,
                tool_args="-m LG+F", output_dir=tmp_path / "out", dry_run=True,
            )

    def test_overwrite_and_resume_mutually_exclusive(self, tmp_path: Path) -> None:
        with pytest.raises(ValueError, match="overwrite"):
            run_alisim_iqtree(
                model_params=tmp_path / "p.tsv", strategy="complete",
                num_simulations=2, overwrite=True, resume=True,
                output_dir=tmp_path / "out", dry_run=True,
            )

    def test_batch_dry_run_result_shape(self, tmp_path: Path) -> None:
        table = tmp_path / "params.tsv"
        table.write_text(
            "id\tseqtype\tlength\tsubs_model\tsubs_rate\tfreq\tprop_inv\t"
            "rate_heterogeneity\trate_categories\trate_param\ttree_path\n"
            "g1\tAA\t2082\tLG\t\t\t\tG\t4\t0.6\t/t/g1.tre\n"
        )
        result = run_alisim_iqtree(
            model_params=table, strategy="complete", num_simulations=2,
            seed=1, output_dir=tmp_path / "out", dry_run=True,
        )
        assert result["status"] == "success"
        assert result["key_results"]["n_simulations_requested"] == 2
        assert result["tool_versions"] == {"iqtree3": "dry-run"}
        assert result["data"]["sampled_rows"]
        assert result["data"]["files"] == []

    def test_batch_completes_with_worker_and_resumes(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch,
    ) -> None:
        from phyloai.posttree import simulate_alisim_iqtree as sim

        table = tmp_path / "params.tsv"
        table.write_text(
            "id\tseqtype\tlength\tsubs_model\tsubs_rate\tfreq\tprop_inv\t"
            "rate_heterogeneity\trate_categories\trate_param\ttree_path\n"
            "g1\tAA\t100\tLG\t\t\t\tG\t4\t0.6\t/t/g1.tre\n"
        )

        monkeypatch.setattr(sim, "_run_simulation_worker", _fake_worker)
        out = tmp_path / "out"
        first = run_alisim_iqtree(
            model_params=table, strategy="complete", num_simulations=2,
            seed=5, output_dir=out, quiet=True,
        )
        assert first["key_results"]["n_simulations_completed"] == 2
        assert first["key_results"]["n_simulations_failed"] == 0
        assert (out / "checkpoint.json").exists()
        assert (out / "params_sampled.tsv").exists()
        assert len(list((out / "MSAs").glob("*.fa"))) == 2
        assert (out / "result.json").exists()

        second = run_alisim_iqtree(
            model_params=table, strategy="complete", num_simulations=2,
            seed=5, output_dir=out, resume=True, quiet=True,
        )
        assert second["key_results"]["n_simulations_completed"] == 2
        assert second["key_results"]["n_simulations_failed"] == 0

    def test_single_mode_output_conflict(self, tmp_path: Path) -> None:
        tree = tmp_path / "ref.tre"
        tree.write_text("(A,B);\n")
        out = tmp_path / "out"
        out.mkdir()
        (out / "existing.txt").write_text("x")
        with pytest.raises(ValueError, match="non-empty"):
            run_alisim_iqtree(
                ref_tree=tree, model="LG", seq_type="AA", length=100,
                output_dir=out, dry_run=False, quiet=True,
            )
