from __future__ import annotations

import csv
import json
from pathlib import Path

import pytest

from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from phyloai.posttree.simulate_adequacy import (
    PreflightError,
    _compute_statistics,
    _compute_taxon_composition,
    _summarize_distribution,
    run_simulate_adequacy,
)


def _msa(*rows: tuple[str, str]) -> MultipleSeqAlignment:
    return MultipleSeqAlignment([
        SeqRecord(Seq(sequence), id=taxon, description="")
        for taxon, sequence in rows
    ])


def test_statistics_known_values_and_gap_exclusion() -> None:
    stats = _compute_statistics(
        _msa(("A", "AA-"), ("B", "AC-"), ("C", "CC-"), ("D", "CA-")),
        "AA",
    )

    assert stats["n_informative_sites"] == 2
    assert stats["div"] == 2.0
    assert stats["siteconvprob"] == 0.5
    assert stats["sitecomp"] == 0.0
    assert stats["comp_max"] == 0.5
    assert stats["comp_mean"] == 0.25


def test_taxon_composition_helper_matches_existing_statistics() -> None:
    alignment = _msa(
        ("A", "AA--"),
        ("B", "AC--"),
        ("C", "CC--"),
        ("D", "CA--"),
    )

    comp = _compute_taxon_composition(alignment, "AA")
    stats = _compute_statistics(alignment, "AA")

    assert comp["taxon_dist_j"] == pytest.approx(stats["taxon_dist_j"])
    assert comp["comp_max"] == pytest.approx(stats["comp_max"])
    assert comp["comp_mean"] == pytest.approx(stats["comp_mean"])
    assert comp["taxon_freqs"]["A"] == pytest.approx([1.0, 0.0] + [0.0] * 18)


def test_pp_directions_match_phylobayes() -> None:
    div = _summarize_distribution([2.0] * 8 + [1.0] * 2, 1.5, "div")
    conv = _summarize_distribution([2.0] * 8 + [1.0] * 2, 1.5, "high")

    assert div["pp"] == 0.2
    assert conv["pp"] == 0.8


def test_zero_sd_uses_json_safe_undefined_pp() -> None:
    summary = _summarize_distribution([1.0] * 10, 1.0, "high")

    assert summary["sd_sim"] == 0.0
    assert summary["z_score"] == 0.0
    assert summary["pp"] is None


def test_nt_ambiguity_and_all_missing_sites_are_ignored() -> None:
    stats = _compute_statistics(_msa(("A", "AN-"), ("B", "CG-")), "NT")

    assert stats["n_informative_sites"] == 2
    assert stats["div"] == 1.5


def test_fewer_than_ten_replicates_is_rejected() -> None:
    with pytest.raises(ValueError, match="at least 10"):
        _summarize_distribution([1.0] * 9, 1.0, "high")


def test_empirical_interval_never_extrapolates_beyond_simulations() -> None:
    summary = _summarize_distribution([float(value) for value in range(10)], 4.5, "high")

    assert 0.0 <= summary["ci_lower"] <= 9.0
    assert 0.0 <= summary["ci_upper"] <= 9.0


def _fasta(path: Path, rows: list[tuple[str, str]]) -> None:
    path.write_text("".join(f">{name}\n{sequence}\n" for name, sequence in rows))


def _write_ten_simulations(directory: Path) -> None:
    directory.mkdir()
    for index in range(10):
        _fasta(directory / f"sim{index}.fa", [("A", "AC"), ("B", "CA")])


def test_preflight_refusal_writes_no_result_json(tmp_path: Path) -> None:
    original = tmp_path / "original.fa"
    simulations = tmp_path / "simulations"
    _fasta(original, [("A", "AC"), ("B", "CA")])
    _write_ten_simulations(simulations)
    out = tmp_path / "out"
    out.mkdir()
    (out / "unrelated.txt").write_text("keep me")

    with pytest.raises(PreflightError):
        run_simulate_adequacy(original_msa=original, simulated_dir=simulations,
                              output_dir=out, quiet=True)

    assert not (out / "result.json").exists()
    assert (out / "unrelated.txt").exists()


def test_unreadable_original_does_not_claim_nonempty_output(tmp_path: Path) -> None:
    simulations = tmp_path / "simulations"
    _write_ten_simulations(simulations)
    out = tmp_path / "out"
    out.mkdir()
    (out / "unrelated.txt").write_text("keep me")

    with pytest.raises(PreflightError):
        run_simulate_adequacy(original_msa=tmp_path / "missing.fa",
                              simulated_dir=simulations, output_dir=out, quiet=True)

    assert not (out / "result.json").exists()
    assert (out / "unrelated.txt").exists()


def test_run_writes_all_tables_and_json_safe_result(tmp_path: Path) -> None:
    original = tmp_path / "original.fa"
    simulations = tmp_path / "simulations"
    _fasta(original, [("A", "AC"), ("B", "CA")])
    _write_ten_simulations(simulations)

    result = run_simulate_adequacy(original_msa=original, simulated_dir=simulations, output_dir=tmp_path / "out", quiet=True)

    assert result["status"] == "success"
    assert result["command"] == (
        f"phyloai posttree simulate adequacy --original-msa {original.resolve()} "
        f"--simulated-dir {simulations.resolve()} --seq-type auto --threads 4 "
        f"-o {(tmp_path / 'out').resolve()} --quiet"
    )
    assert isinstance(result["wall_time"], float)
    assert result["tool_versions"] == {}
    assert result["error"] is None
    assert result["key_results"]["seq_type"] == "NT"
    assert set(result["key_results"]["statistics"]) == {"div", "siteconvprob", "sitecomp", "comp"}
    assert set(result["key_results"]["statistics"]["comp"]) == {"max", "mean"}
    assert result["data"]["cmd"] == []
    assert result["data"]["tool_stderr"] == ""
    assert (tmp_path / "out" / "adequacy_summary.csv").exists()
    assert (tmp_path / "out" / "adequacy_taxon_comp.csv").exists()
    assert (tmp_path / "out" / "per_simulation_stats.csv").exists()
    assert json.loads((tmp_path / "out" / "result.json").read_text())["error"] is None


def test_duplicate_original_taxon_is_a_hard_error(tmp_path: Path) -> None:
    original = tmp_path / "original.fa"
    simulations = tmp_path / "simulations"
    _fasta(original, [("A", "AC"), ("A", "CA")])
    _write_ten_simulations(simulations)

    with pytest.raises(ValueError, match="duplicate taxon"):
        run_simulate_adequacy(original_msa=original, simulated_dir=simulations, output_dir=tmp_path / "out", quiet=True)


def test_invalid_original_does_not_modify_nonempty_output(tmp_path: Path) -> None:
    original = tmp_path / "original.fa"
    simulations = tmp_path / "simulations"
    output_dir = tmp_path / "out"
    _fasta(original, [("A", "AC"), ("A", "CA")])
    _write_ten_simulations(simulations)
    output_dir.mkdir()
    (output_dir / "result.json").write_text('{"keep": true}')

    with pytest.raises(PreflightError, match="already exists and is non-empty"):
        run_simulate_adequacy(
            original_msa=original, simulated_dir=simulations,
            output_dir=output_dir, quiet=True,
        )

    assert (output_dir / "result.json").read_text() == '{"keep": true}'


def test_invalid_original_is_rejected_before_simulations_run(tmp_path: Path) -> None:
    original = tmp_path / "original.fa"
    simulations = tmp_path / "simulations"
    _fasta(original, [("A", "--"), ("B", "--")])
    _write_ten_simulations(simulations)

    with pytest.raises(ValueError, match="no informative sites"):
        run_simulate_adequacy(
            original_msa=original, simulated_dir=simulations,
            output_dir=tmp_path / "out", quiet=True,
        )

    assert not (tmp_path / "out" / "checkpoint.json").exists()


def test_bad_simulation_is_skipped_and_valid_taxa_are_remapped(tmp_path: Path) -> None:
    original = tmp_path / "original.fa"
    simulations = tmp_path / "simulations"
    _fasta(original, [("A", "AC"), ("B", "CA")])
    _write_ten_simulations(simulations)
    _fasta(simulations / "duplicate.fa", [("A", "AC"), ("A", "CA")])
    _fasta(simulations / "reordered.fa", [("B", "CA"), ("A", "AC")])

    result = run_simulate_adequacy(original_msa=original, simulated_dir=simulations, output_dir=tmp_path / "out", quiet=True)

    assert result["key_results"]["n_failed"] == 1
    assert any("duplicate" in warning for warning in result["data"]["warnings"])
    with open(tmp_path / "out" / "adequacy_taxon_comp.csv", newline="") as handle:
        assert {row["taxon"] for row in csv.DictReader(handle)} == {"A", "B"}


def test_resume_uses_saved_taxon_values_without_reprocessing(tmp_path: Path) -> None:
    original = tmp_path / "original.fa"
    simulations = tmp_path / "simulations"
    _fasta(original, [("A", "AC"), ("B", "CA")])
    _write_ten_simulations(simulations)
    output_dir = tmp_path / "out"
    first = run_simulate_adequacy(original_msa=original, simulated_dir=simulations, output_dir=output_dir, quiet=True)
    checkpoint = json.loads((output_dir / "checkpoint.json").read_text())
    assert "taxon_dist_j" in checkpoint["tasks"][0]["outputs"]
    expected_taxon_csv = (output_dir / "adequacy_taxon_comp.csv").read_text()

    resumed = run_simulate_adequacy(original_msa=original, simulated_dir=simulations, output_dir=output_dir, resume=True, quiet=True)

    assert resumed["key_results"]["n_simulations"] == first["key_results"]["n_simulations"]
    assert (output_dir / "adequacy_taxon_comp.csv").read_text() == expected_taxon_csv


def test_nonempty_output_requires_overwrite_or_resume(tmp_path: Path) -> None:
    original = tmp_path / "original.fa"
    simulations = tmp_path / "simulations"
    output_dir = tmp_path / "out"
    _fasta(original, [("A", "AC"), ("B", "CA")])
    _write_ten_simulations(simulations)
    output_dir.mkdir()
    (output_dir / "existing.txt").write_text("keep")

    with pytest.raises(PreflightError, match="already exists and is non-empty"):
        run_simulate_adequacy(original_msa=original, simulated_dir=simulations, output_dir=output_dir, quiet=True)


def test_resume_requires_checkpoint(tmp_path: Path) -> None:
    original = tmp_path / "original.fa"
    simulations = tmp_path / "simulations"
    output_dir = tmp_path / "out"
    _fasta(original, [("A", "AC"), ("B", "CA")])
    _write_ten_simulations(simulations)
    output_dir.mkdir()

    with pytest.raises(PreflightError, match="No checkpoint found"):
        run_simulate_adequacy(original_msa=original, simulated_dir=simulations, output_dir=output_dir, resume=True, quiet=True)


def test_replaced_simulation_is_recomputed_on_resume(tmp_path: Path) -> None:
    original = tmp_path / "original.fa"
    simulations = tmp_path / "simulations"
    _fasta(original, [("A", "AC"), ("B", "CA")])
    _write_ten_simulations(simulations)
    output_dir = tmp_path / "out"
    run_simulate_adequacy(original_msa=original, simulated_dir=simulations, output_dir=output_dir, quiet=True)
    _fasta(simulations / "sim0.fa", [("A", "AA"), ("B", "CC")])

    run_simulate_adequacy(original_msa=original, simulated_dir=simulations, output_dir=output_dir, resume=True, quiet=True)

    checkpoint = json.loads((output_dir / "checkpoint.json").read_text())
    task = next(task for task in checkpoint["tasks"] if task["task_id"].endswith("sim0.fa"))
    assert task["input"].split("|")[1:] == [str((simulations / "sim0.fa").stat().st_size), str((simulations / "sim0.fa").stat().st_mtime_ns)]


def test_dry_run_writes_nothing(tmp_path: Path) -> None:
    original = tmp_path / "original.fa"
    simulations = tmp_path / "simulations"
    _fasta(original, [("A", "AC"), ("B", "CA")])
    _write_ten_simulations(simulations)
    output_dir = tmp_path / "out"

    result = run_simulate_adequacy(original_msa=original, simulated_dir=simulations, output_dir=output_dir, dry_run=True, quiet=True)

    assert result["status"] == "success"
    assert not output_dir.exists()


def test_table_format_tsv_uses_tsv_suffixes_and_delimiters(tmp_path: Path) -> None:
    original = tmp_path / "original.fa"
    simulations = tmp_path / "simulations"
    _fasta(original, [("A", "AC"), ("B", "CA")])
    _write_ten_simulations(simulations)

    result = run_simulate_adequacy(
        original_msa=original, simulated_dir=simulations,
        table_format="tsv", output_dir=tmp_path / "out", quiet=True,
    )

    summary = tmp_path / "out" / "adequacy_summary.tsv"
    assert "\t" in summary.read_text().splitlines()[0]
    assert result["data"]["output_files"]["adequacy_summary"]["path"] == str(summary)


def test_progress_callback_counts_only_pending_resume_tasks(tmp_path: Path) -> None:
    original = tmp_path / "original.fa"
    simulations = tmp_path / "simulations"
    _fasta(original, [("A", "AC"), ("B", "CA")])
    _write_ten_simulations(simulations)
    output_dir = tmp_path / "out"
    run_simulate_adequacy(original_msa=original, simulated_dir=simulations, output_dir=output_dir, quiet=True)
    _fasta(simulations / "sim0.fa", [("A", "AA"), ("B", "CC")])
    updates: list[tuple[int, int]] = []

    run_simulate_adequacy(
        original_msa=original, simulated_dir=simulations, output_dir=output_dir,
        resume=True, quiet=True, progress_callback=lambda done, total: updates.append((done, total)),
    )

    assert updates == [(1, 1)]


def test_changed_original_blocks_resume(tmp_path: Path) -> None:
    original = tmp_path / "original.fa"
    simulations = tmp_path / "simulations"
    _fasta(original, [("A", "AC"), ("B", "CA")])
    _write_ten_simulations(simulations)
    output_dir = tmp_path / "out"
    run_simulate_adequacy(original_msa=original, simulated_dir=simulations, output_dir=output_dir, quiet=True)
    _fasta(original, [("A", "AA"), ("B", "CC")])

    with pytest.raises(PreflightError, match="original.*changed"):
        run_simulate_adequacy(original_msa=original, simulated_dir=simulations, output_dir=output_dir, resume=True, quiet=True)


def test_legacy_checkpoint_without_original_fingerprint_is_rejected(tmp_path: Path) -> None:
    original = tmp_path / "original.fa"
    simulations = tmp_path / "simulations"
    _fasta(original, [("A", "AC"), ("B", "CA")])
    _write_ten_simulations(simulations)
    output_dir = tmp_path / "out"
    run_simulate_adequacy(original_msa=original, simulated_dir=simulations,
                          output_dir=output_dir, quiet=True)
    checkpoint_path = output_dir / "checkpoint.json"
    checkpoint = json.loads(checkpoint_path.read_text())
    checkpoint.pop("original_msa_fingerprint")
    checkpoint_path.write_text(json.dumps(checkpoint))

    with pytest.raises(PreflightError, match="fingerprint.*overwrite"):
        run_simulate_adequacy(original_msa=original, simulated_dir=simulations,
                              output_dir=output_dir, resume=True, quiet=True)


def test_resume_adds_new_and_drops_deleted_simulations(tmp_path: Path) -> None:
    original = tmp_path / "original.fa"
    simulations = tmp_path / "simulations"
    _fasta(original, [("A", "AC"), ("B", "CA")])
    _write_ten_simulations(simulations)
    output_dir = tmp_path / "out"
    run_simulate_adequacy(original_msa=original, simulated_dir=simulations, output_dir=output_dir, quiet=True)
    (simulations / "sim9.fa").unlink()
    _fasta(simulations / "sim10.fa", [("A", "AC"), ("B", "CA")])

    result = run_simulate_adequacy(original_msa=original, simulated_dir=simulations, output_dir=output_dir, resume=True, quiet=True)

    assert result["key_results"]["n_simulations"] == 10
    assert len(json.loads((output_dir / "checkpoint.json").read_text())["tasks"]) == 10


def test_zscore_fixture_regression(tmp_path: Path) -> None:
    repo_root = Path(__file__).resolve().parents[2]
    original = repo_root / "runs/zscore/matrix.XX"
    if not original.exists():
        original = repo_root / "runs/zscore/matrix.fa"
    simulated = repo_root / "runs/zscore/simulated/fasta"
    output_dir = tmp_path / "out"

    result = run_simulate_adequacy(original_msa=original, simulated_dir=simulated, output_dir=output_dir, quiet=True)

    stats = result["key_results"]["statistics"]
    expected = {
        "div": (2.8208936170212775, 0.0712155216268796, 0.4726437958910017, 0.32),
        "siteconvprob": (0.5508971158392437, 0.014775408637274967, 0.2346432805833419, 0.38),
        "sitecomp": (0.024310778713344406, 0.0007377783766110406, 0.41786635622153806, 0.31),
        "comp_max": (0.0030999442049923775, 0.0007365073338196676, 1.174913271897109, 0.15),
        "comp_mean": (0.0019202401184298704, 0.00031487384860307124, 1.2742692228464887, 0.1),
    }
    mapping = {
        "div": stats["div"],
        "siteconvprob": stats["siteconvprob"],
        "sitecomp": stats["sitecomp"],
        "comp_max": stats["comp"]["max"],
        "comp_mean": stats["comp"]["mean"],
    }
    for name, (mean_sim, sd_sim, z_score, pp) in expected.items():
        got = mapping[name]
        assert got["mean_sim"] == pytest.approx(mean_sim, abs=1e-12)
        assert got["sd_sim"] == pytest.approx(sd_sim, abs=1e-12)
        assert got["z_score"] == pytest.approx(z_score, abs=1e-12)
        assert got["pp"] == pytest.approx(pp, abs=1e-12)
