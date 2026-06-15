"""Integration and unit tests for phyloai.pretree.metrics."""

from __future__ import annotations

import csv
import itertools
import json
import tempfile
from pathlib import Path

import numpy as np
import pytest
from Bio import Phylo

from phyloai.pretree.metrics import (
    _check_taxon_consistency,
    _compute_correlation,
    _compute_frequencies,
    _compute_msa_metrics,
    _compute_saturation,
    _compute_tree_metrics,
    _generate_all_plots,
    _generate_basic_statistics,
    _generate_correlation_heatmap,
    _pair_files,
    _plot_single_metric,
    _select_correlation_columns,
    _scan_msa_headers,
    _strip_tree_suffixes,
    _write_correlation_csv,
    run_metrics,
)

# ---------------------------------------------------------------------------
# Helpers for test data
# ---------------------------------------------------------------------------

_AA_SEQS = [
    "A-RCD--E",
    "AR-CD-LE",
    "A-RCDT-E",
]

_NT_SEQS = [
    "A-CT--G",
    "A-CT-NG",
    "AGCT--G",
]

_NEWICK_TREE = "(A:0.1,B:0.2,(C:0.3,D:0.4):0.5);"


def _write_fasta(path: Path, sequences: list[str], seq_ids: list[str] | None = None) -> None:
    if seq_ids is None:
        seq_ids = [f"seq{i}" for i in range(len(sequences))]
    with open(path, "w") as fh:
        for sid, seq in zip(seq_ids, sequences):
            fh.write(f">{sid}\n{seq}\n")


def _write_newick(path: Path, tree_str: str) -> None:
    with open(path, "w") as fh:
        fh.write(tree_str + "\n")


# ---------------------------------------------------------------------------
# Task 1 tests
# ---------------------------------------------------------------------------


class TestStripTreeSuffixes:
    def test_iqtree_compound(self):
        assert _strip_tree_suffixes("EOG090X002Z.fas.treefile") == "EOG090X002Z"

    def test_contree(self):
        assert _strip_tree_suffixes("gene2.contree") == "gene2"

    def test_simple_tre(self):
        assert _strip_tree_suffixes("simple.tre") == "simple"

    def test_no_suffix_change(self):
        assert _strip_tree_suffixes("genome") == "genome"

    def test_treefile_suffix(self):
        assert _strip_tree_suffixes("gene1.fa.treefile") == "gene1"


class TestScanMsaHeaders:
    def test_basic(self, tmp_path):
        for i in range(3):
            _write_fasta(tmp_path / f"gene{i}.fa", _AA_SEQS, [f"taxon{j}" for j in range(3)])
        per_marker, total_pool = _scan_msa_headers(tmp_path)
        assert len(per_marker) == 3
        assert len(total_pool) == 3
        assert per_marker["gene0"] == {"taxon0", "taxon1", "taxon2"}


class TestPairFiles:
    def test_msa_tree_pairing(self, tmp_path):
        msa_dir = tmp_path / "msa"
        tree_dir = tmp_path / "trees"
        msa_dir.mkdir()
        tree_dir.mkdir()
        _write_fasta(msa_dir / "gene1.fa", _AA_SEQS, ["A", "B", "C"])
        _write_newick(tree_dir / "gene1.fas.treefile", "(A,B,C);")
        _write_fasta(msa_dir / "gene2.fa", _AA_SEQS, ["D", "E", "F"])
        _write_newick(tree_dir / "gene2.tre", "(D,E,F);")

        paired, warnings = _pair_files(msa_dir, tree_dir)
        assert len(paired) == 2
        assert paired["gene1"][0] is not None  # msa
        assert paired["gene1"][1] is not None  # tree
        assert paired["gene2"][0] is not None
        assert paired["gene2"][1] is not None

    def test_tree_only_warns(self, tmp_path):
        tree_dir = tmp_path / "trees"
        tree_dir.mkdir()
        _write_newick(tree_dir / "orphan.tre", "(A,B);")
        paired, warnings = _pair_files(None, tree_dir)
        assert "orphan" in paired
        assert len(warnings) > 0


class TestCheckTaxonConsistency:
    def test_match(self, tmp_path):
        msa_p = tmp_path / "test.fa"
        tree_p = tmp_path / "test.tre"
        _write_fasta(msa_p, ["ACG", "TGC", "CAT"], ["A", "B", "C"])
        _write_newick(tree_p, "(A,B,C);")
        assert _check_taxon_consistency(msa_p, tree_p) is None

    def test_mismatch(self, tmp_path):
        msa_p = tmp_path / "test.fa"
        tree_p = tmp_path / "test.tre"
        _write_fasta(msa_p, ["ACG", "TGC"], ["A", "B"])
        _write_newick(tree_p, "(A,C);")
        result = _check_taxon_consistency(msa_p, tree_p)
        assert result is not None
        assert "B" in result["msa_only"]
        assert "C" in result["tree_only"]


# ---------------------------------------------------------------------------
# Task 2 tests: MSA metrics
# ---------------------------------------------------------------------------


class TestComputeMsaMetrics:
    @staticmethod
    def _msa_data(sequences: list[str]) -> list[list[str]]:
        return [[c.upper() for c in s] for s in sequences]

    def test_basic_counts(self):
        msa = self._msa_data(_AA_SEQS)
        result = _compute_msa_metrics(msa, "AA", total_taxa_pool=3, skip_freq=True)
        assert result["num_taxa"] == 3
        assert result["num_sites"] == 8
        assert result["taxa_occupancy"] == 1.0

    def test_proportion_gaps(self):
        msa = self._msa_data(_AA_SEQS)
        result = _compute_msa_metrics(msa, "AA", total_taxa_pool=3, skip_freq=True)
        assert 0.0 < result["proportion_gaps"] < 1.0

    def test_num_patterns(self):
        msa = self._msa_data(_AA_SEQS)
        result = _compute_msa_metrics(msa, "AA", total_taxa_pool=3, skip_freq=True)
        assert result["num_patterns"] >= 1

    def test_num_singletons(self):
        msa = self._msa_data(_NT_SEQS)
        result = _compute_msa_metrics(msa, "NT", total_taxa_pool=3, skip_freq=True)
        assert result["num_singletons"] >= 0

    def test_rcfv(self):
        msa = self._msa_data(_AA_SEQS)
        result = _compute_msa_metrics(msa, "AA", total_taxa_pool=3, skip_freq=True)
        assert result["rcfv"] >= 0.0

    def test_nrcfv(self):
        msa = self._msa_data(_AA_SEQS)
        result = _compute_msa_metrics(msa, "AA", total_taxa_pool=3, skip_freq=True)
        assert result["nrcfv"] >= 0.0

    def test_gc_content_nt(self):
        msa = self._msa_data(_NT_SEQS)
        result = _compute_msa_metrics(msa, "NT", total_taxa_pool=3, skip_freq=True)
        assert isinstance(result["GC_content"], float)

    def test_gc_content_aa_empty(self):
        msa = self._msa_data(_AA_SEQS)
        result = _compute_msa_metrics(msa, "AA", total_taxa_pool=3, skip_freq=True)
        assert result["GC_content"] == ""

    def test_frequency_statistics(self):
        msa = self._msa_data(_AA_SEQS)
        result = _compute_msa_metrics(msa, "AA", total_taxa_pool=3, skip_freq=False)
        for label in "ARNDCQEGHILKMFPSTWYV":
            assert f"freq{label}" in result

    def test_empty_msa(self):
        result = _compute_msa_metrics([], "AA", total_taxa_pool=0, skip_freq=True)
        assert result["num_sites"] == 0

    def test_proportion_invariant(self):
        # All sequences identical
        msa = self._msa_data(["AAAA", "AAAA", "AAAA"])
        result = _compute_msa_metrics(msa, "NT", total_taxa_pool=3, skip_freq=True)
        assert result["proportion_invariant"] == 1.0


# ---------------------------------------------------------------------------
# Task 3 tests: Tree metrics
# ---------------------------------------------------------------------------


class TestComputeTreeMetrics:
    def test_basic_tree(self, tmp_path):
        tree_path = tmp_path / "test.tre"
        _write_newick(tree_path, "(A:0.1,B:0.2,(C:0.3,D:0.4):0.5);")
        result = _compute_tree_metrics(tree_path, None, None)
        assert result["total_tree_length"] > 0
        assert result["treeness"] > 0

    def test_bs_parsing(self, tmp_path):
        tree_path = tmp_path / "test.tre"
        _write_newick(tree_path, "(A:0.1,B:0.2)0.95:0.1;")
        result = _compute_tree_metrics(tree_path, None, None)
        assert 0 < result["average_BS"] <= 100

    def test_dvmc(self, tmp_path):
        tree_path = tmp_path / "test.tre"
        _write_newick(tree_path, "(A:0.1,B:0.2,(C:0.3,D:0.4):0.5);")
        result = _compute_tree_metrics(tree_path, None, None)
        assert result["dvmc"] >= 0.0

    def test_rf_distance(self, tmp_path):
        tree_path = tmp_path / "test.tre"
        ref_path = tmp_path / "ref.tre"
        _write_newick(tree_path, "(A:0.1,B:0.2,(C:0.3,D:0.4):0.5);")
        _write_newick(ref_path, "(A,B,(C,D));")
        result = _compute_tree_metrics(tree_path, None, ref_path)
        assert result["RF_distance"] >= 0.0

class TestComputeSaturation:
    def test_slope_perfect(self, tmp_path):
        msa_p = tmp_path / "test.fa"
        tree_p = tmp_path / "test.tre"
        _write_fasta(msa_p, ["AAA", "AAA", "AAA", "AAA"], ["A", "B", "C", "D"])
        _write_newick(tree_p, "(A:0.5,B:0.5,(C:0.3,D:0.3):0.2);")
        slope = _compute_saturation(msa_p, tree_p)
        assert 0.0 <= slope <= 1.0


# ---------------------------------------------------------------------------
# Task 5 tests: Orchestration
# ---------------------------------------------------------------------------


class TestRunMetrics:
    def test_full_pipeline(self, tmp_path):
        msa_dir = tmp_path / "msa"
        tree_dir = tmp_path / "trees"
        msa_dir.mkdir()
        tree_dir.mkdir()
        out_dir = tmp_path / "output"
        out_dir.mkdir()

        _write_fasta(msa_dir / "gene1.fa", _AA_SEQS, ["A", "B", "C"])
        _write_fasta(msa_dir / "gene2.fa", _AA_SEQS, ["A", "B", "C"])
        _write_newick(tree_dir / "gene1.tre", "(A,B,C);")
        _write_newick(tree_dir / "gene2.tre", "(A,B,C);")

        result = run_metrics(
            msa_dir=msa_dir,
            tree_dir=tree_dir,
            output_dir=out_dir,
            skip_freq=True,
            skip_pairwise_identity=True,
            quiet=True,
        )
        assert result["status"] == "success"
        assert result["key_results"]["n_markers"] == 2
        assert result["key_results"]["n_success"] == 2

    def test_tree_only(self, tmp_path):
        tree_dir = tmp_path / "trees"
        tree_dir.mkdir()
        out_dir = tmp_path / "output"

        _write_newick(tree_dir / "gene1.tre", "(A,B,C);")

        result = run_metrics(
            tree_dir=tree_dir,
            output_dir=out_dir,
            skip_freq=True,
            skip_pairwise_identity=True,
            quiet=True,
        )
        assert result["status"] == "success"

    def test_dry_run(self, tmp_path):
        msa_dir = tmp_path / "msa"
        msa_dir.mkdir()
        out_dir = tmp_path / "output"

        _write_fasta(msa_dir / "gene1.fa", _AA_SEQS, ["A", "B", "C"])

        result = run_metrics(
            msa_dir=msa_dir,
            output_dir=out_dir,
            skip_pairwise_identity=True,
            dry_run=True,
            quiet=True,
        )
        assert result["key_results"]["dry_run"] is True
        assert not out_dir.exists()

    def test_no_input_dirs(self, tmp_path):
        result = run_metrics(output_dir=tmp_path / "output", quiet=True)
        assert result["status"] == "error"

    def test_csv_columns(self, tmp_path):
        msa_dir = tmp_path / "msa"
        msa_dir.mkdir()
        out_dir = tmp_path / "output"

        _write_fasta(msa_dir / "gene1.fa", _AA_SEQS, ["A", "B", "C"])

        run_metrics(
            msa_dir=msa_dir,
            output_dir=out_dir,
            skip_freq=True,
            skip_pairwise_identity=True,
            quiet=True,
        )
        csv_path = out_dir / "metrics.csv"
        assert csv_path.exists()
        with open(csv_path) as f:
            reader = csv.DictReader(f)
            row = next(reader)
            assert "loci" in row
            assert "DataType" in row
            assert "num_taxa" in row
            assert "rcfv" in row

    def test_taxon_mismatch_warning(self, tmp_path):
        msa_dir = tmp_path / "msa"
        tree_dir = tmp_path / "trees"
        msa_dir.mkdir()
        tree_dir.mkdir()
        out_dir = tmp_path / "output"

        _write_fasta(msa_dir / "gene1.fa", _AA_SEQS, ["A", "B", "C"])
        _write_newick(tree_dir / "gene1.tre", "(A,B,D);")

        result = run_metrics(
            msa_dir=msa_dir,
            tree_dir=tree_dir,
            output_dir=out_dir,
            skip_freq=True,
            skip_pairwise_identity=True,
            quiet=True,
        )
        warnings = result["key_results"].get("warnings") or []
        assert any("taxa" in w.lower() for w in warnings)


# ---------------------------------------------------------------------------
# Task 6 tests: Plots
# ---------------------------------------------------------------------------


class TestPlotSingleMetric:
    def test_creates_pdf(self, tmp_path):
        data = np.random.normal(0, 1, 50)
        out = tmp_path / "test.pdf"
        _plot_single_metric(data, "test", out)
        assert out.exists()
        assert out.stat().st_size > 100

    def test_tukey_filter(self, tmp_path):
        data = np.array([1.0, 2.0, 3.0, 2.5, 2.0, 1.5, 100.0])
        out = tmp_path / "test_tukey.pdf"
        _plot_single_metric(data, "test", out, tukey_k=1.5)
        assert out.exists()

    def test_empty_data(self, tmp_path):
        out = tmp_path / "empty.pdf"
        _plot_single_metric(np.array([]), "empty", out)
        # Should not crash

    def test_uses_requested_bar_color(self, tmp_path, monkeypatch):
        captured = {}

        import matplotlib.axes
        original_hist = matplotlib.axes.Axes.hist

        def spy_hist(self, *args, **kwargs):
            captured["color"] = kwargs.get("color")
            return original_hist(self, *args, **kwargs)

        monkeypatch.setattr(matplotlib.axes.Axes, "hist", spy_hist)

        out = tmp_path / "colored.pdf"
        _plot_single_metric(np.array([1.0, 2.0, 3.0]), "test", out, color="#123456")

        assert captured["color"] == "#123456"


class TestGenerateAllPlots:
    def test_multiple_metrics(self, tmp_path):
        rows = [
            {"loci": "g1", "DataType": "AA", "a": "1.0", "b": "2.0"},
            {"loci": "g2", "DataType": "AA", "a": "3.0", "b": "4.0"},
        ]
        plots_dir = tmp_path / "plots"
        count = _generate_all_plots(rows, ["a", "b"], plots_dir)
        assert count == 2
        assert (plots_dir / "a.pdf").exists()
        assert (plots_dir / "b.pdf").exists()


class TestGenerateBasicStatistics:
    def test_stats_csv(self, tmp_path):
        rows = [
            {"a": "1.0", "b": "2.0"},
            {"a": "3.0", "b": "4.0"},
            {"a": "5.0", "b": "6.0"},
        ]
        out = tmp_path / "stats.csv"
        _generate_basic_statistics(rows, ["a", "b"], out)
        assert out.exists()
        with open(out) as f:
            reader = csv.DictReader(f)
            rows_read = list(reader)
            assert len(rows_read) == 2


# ---------------------------------------------------------------------------
# Task 7 tests: Correlation
# ---------------------------------------------------------------------------


class TestComputeCorrelation:
    def test_default_selection_excludes_frequency_and_sd_columns(self):
        rows = [
            {"loci": "g1", "DataType": "AA", "freqA": "0.1", "num_sites": "10", "sd_BS": "1", "entropy": "0.2"},
            {"loci": "g2", "DataType": "AA", "freqA": "0.2", "num_sites": "20", "sd_BS": "2", "entropy": "0.4"},
        ]
        columns = list(rows[0].keys())

        selected = _select_correlation_columns(rows, columns)

        assert selected == ["num_sites", "entropy"]

    def test_include_flags_add_frequency_and_sd_columns(self):
        rows = [
            {"loci": "g1", "DataType": "AA", "freqA": "0.1", "num_sites": "10", "sd_BS": "1", "entropy": "0.2"},
            {"loci": "g2", "DataType": "AA", "freqA": "0.2", "num_sites": "20", "sd_BS": "2", "entropy": "0.4"},
        ]
        columns = list(rows[0].keys())

        selected = _select_correlation_columns(rows, columns, include_freq=True, include_sd=True)

        assert selected == ["freqA", "num_sites", "sd_BS", "entropy"]

    def test_metrics_all_keeps_all_numeric_columns(self):
        rows = [
            {"loci": "g1", "DataType": "AA", "freqA": "0.1", "num_sites": "10", "sd_BS": "1"},
            {"loci": "g2", "DataType": "AA", "freqA": "0.2", "num_sites": "20", "sd_BS": "2"},
        ]
        columns = list(rows[0].keys())

        selected = _select_correlation_columns(rows, columns, requested="all")

        assert selected == ["freqA", "num_sites", "sd_BS"]

    def test_explicit_metrics_are_preserved(self):
        rows = [
            {"loci": "g1", "DataType": "AA", "freqA": "0.1", "num_sites": "10", "sd_BS": "1"},
            {"loci": "g2", "DataType": "AA", "freqA": "0.2", "num_sites": "20", "sd_BS": "2"},
        ]
        columns = list(rows[0].keys())

        selected = _select_correlation_columns(rows, columns, requested="freqA,sd_BS")

        assert selected == ["freqA", "sd_BS"]

    def test_spearman(self):
        rows = [
            {"a": "1", "b": "2", "c": "5"},
            {"a": "2", "b": "4", "c": "4"},
            {"a": "3", "b": "6", "c": "3"},
            {"a": "4", "b": "8", "c": "2"},
            {"a": "5", "b": "10", "c": "1"},
        ]
        corr, names = _compute_correlation(rows, ["a", "b", "c"], method="spearman")
        assert corr.shape == (3, 3)
        assert abs(corr[0, 1]) > 0.99  # a vs b near-perfect positive
        assert abs(corr[0, 2] + 1.0) < 0.01  # a vs c near-perfect negative

    def test_constant_column_excluded(self):
        rows = [
            {"a": "1", "b": "5", "c": "9"},
            {"a": "1", "b": "4", "c": "10"},
            {"a": "1", "b": "3", "c": "8"},
        ]
        corr, names = _compute_correlation(rows, ["a", "b", "c"], method="spearman")
        assert len(names) == 2  # 'b' and 'c' - 'a' is constant so excluded

    def test_all_na_excluded(self):
        rows = [
            {"a": "1", "b": "3"},
            {"a": "2", "b": "4"},
            {"a": "3", "b": "1"},
        ]
        corr, names = _compute_correlation(rows, ["a", "b", "x"], method="spearman")
        assert len(names) == 2  # 'a' and 'b' - 'x' all NA so excluded

    def test_insufficient_columns(self):
        rows = [{"a": "1"}, {"a": "2"}]
        corr, names = _compute_correlation(rows, ["a"], method="spearman")
        assert corr.size == 0


class TestGenerateCorrelationHeatmap:
    def test_creates_pdf(self, tmp_path):
        corr = np.array([[1.0, 0.8], [0.8, 1.0]])
        out = tmp_path / "heatmap.pdf"
        _generate_correlation_heatmap(corr, ["a", "b"], out)
        assert out.exists()

    def test_two_variables_no_dendrogram(self, tmp_path):
        corr = np.array([[1.0, 0.5], [0.5, 1.0]])
        out = tmp_path / "heatmap2.pdf"
        _generate_correlation_heatmap(corr, ["a", "b"], out, triangle="lower")
        assert out.exists()

    def test_circle_annotations_create_pdf(self, tmp_path):
        corr = np.array([[1.0, -0.7], [-0.7, 1.0]])
        out = tmp_path / "annotated_circle.pdf"
        _generate_correlation_heatmap(corr, ["a", "b"], out, annot=True)
        assert out.exists()
        assert out.stat().st_size > 100

    def test_circle_heatmap_uses_single_colorbar_without_clustermap_axes(self, tmp_path, monkeypatch):
        axes_counts = []

        import matplotlib.figure
        original_savefig = matplotlib.figure.Figure.savefig

        def spy_savefig(self, *args, **kwargs):
            axes_counts.append(len(self.axes))
            return original_savefig(self, *args, **kwargs)

        monkeypatch.setattr(matplotlib.figure.Figure, "savefig", spy_savefig)

        corr = np.eye(4)
        corr[corr == 0] = 0.4
        out = tmp_path / "single_colorbar.pdf"
        _generate_correlation_heatmap(corr, ["a", "b", "c", "d"], out)

        assert axes_counts[-1] == 2

    def test_triangle_circle_heatmap_skips_masked_grid_segments(self, tmp_path, monkeypatch):
        captured_segments = []

        import matplotlib.collections
        original_init = matplotlib.collections.LineCollection.__init__

        def spy_init(self, segments, *args, **kwargs):
            captured_segments.extend(segments)
            original_init(self, segments, *args, **kwargs)

        monkeypatch.setattr(matplotlib.collections.LineCollection, "__init__", spy_init)

        corr = np.eye(3)
        corr[0, 1] = corr[1, 0] = 0.5
        corr[0, 2] = corr[2, 0] = -0.2
        corr[1, 2] = corr[2, 1] = 0.8
        out = tmp_path / "lower_triangle.pdf"
        _generate_correlation_heatmap(corr, ["a", "b", "c"], out, triangle="lower")

        assert captured_segments
        assert all(not (x1 == 2 and x2 == 2 and y1 == 0 and y2 == 1)
                   for ((x1, y1), (x2, y2)) in captured_segments)

    def test_cluster_rectangles_are_full_mode_only(self, tmp_path, monkeypatch):
        rectangles = []

        import matplotlib.axes
        original_add_patch = matplotlib.axes.Axes.add_patch

        def spy_add_patch(self, patch):
            if patch.__class__.__name__ == "Rectangle":
                rectangles.append(patch)
            return original_add_patch(self, patch)

        monkeypatch.setattr(matplotlib.axes.Axes, "add_patch", spy_add_patch)

        corr = np.array([
            [1.0, 0.9, 0.1, 0.1],
            [0.9, 1.0, 0.1, 0.1],
            [0.1, 0.1, 1.0, 0.9],
            [0.1, 0.1, 0.9, 1.0],
        ])

        _generate_correlation_heatmap(corr, ["a", "b", "c", "d"], tmp_path / "lower_rect.pdf",
                                      triangle="lower", cluster_rectangles=2)
        lower_count = len(rectangles)

        _generate_correlation_heatmap(corr, ["a", "b", "c", "d"], tmp_path / "full_rect.pdf",
                                      triangle="full", cluster_rectangles=2)

        assert lower_count == 0
        assert len(rectangles) > lower_count

    def test_cluster_rectangles_warn_when_ignored_for_triangle_mode(self, tmp_path):
        warnings = []
        corr = np.eye(3)

        _generate_correlation_heatmap(
            corr,
            ["a", "b", "c"],
            tmp_path / "ignored_rectangles.pdf",
            triangle="upper",
            cluster_rectangles=2,
            warn=warnings.append,
        )

        assert any("cluster-rectangles" in message and "full" in message for message in warnings)

    def test_upper_triangle_moves_ticks_to_top_and_right(self, tmp_path, monkeypatch):
        positions = {}

        import matplotlib.figure
        original_savefig = matplotlib.figure.Figure.savefig

        def spy_savefig(self, *args, **kwargs):
            heatmap_ax = self.axes[0]
            positions["x_ticks"] = heatmap_ax.xaxis.get_ticks_position()
            positions["x_label"] = heatmap_ax.xaxis.get_label_position()
            positions["y_ticks"] = heatmap_ax.yaxis.get_ticks_position()
            positions["y_label"] = heatmap_ax.yaxis.get_label_position()
            positions["x_rotation"] = heatmap_ax.get_xticklabels()[0].get_rotation()
            positions["x_ha"] = heatmap_ax.get_xticklabels()[0].get_ha()
            positions["x_va"] = heatmap_ax.get_xticklabels()[0].get_va()
            return original_savefig(self, *args, **kwargs)

        monkeypatch.setattr(matplotlib.figure.Figure, "savefig", spy_savefig)

        corr = np.eye(3)
        out = tmp_path / "upper_ticks.pdf"
        _generate_correlation_heatmap(corr, ["a", "b", "c"], out, triangle="upper", label_angle=30)

        assert positions["x_ticks"] == "top"
        assert positions["x_label"] == "top"
        assert positions["y_ticks"] == "right"
        assert positions["y_label"] == "right"
        assert positions["x_rotation"] == 30
        assert positions["x_ha"] == "left"
        assert positions["x_va"] == "bottom"

    def test_triangle_heatmap_uses_triangle_border_not_rectangular_spines(self, tmp_path, monkeypatch):
        state = {}

        import matplotlib.figure
        original_savefig = matplotlib.figure.Figure.savefig

        def spy_savefig(self, *args, **kwargs):
            heatmap_ax = self.axes[0]
            state["spines_visible"] = [spine.get_visible() for spine in heatmap_ax.spines.values()]
            state["lines"] = [line.get_xydata().tolist() for line in heatmap_ax.lines]
            return original_savefig(self, *args, **kwargs)

        monkeypatch.setattr(matplotlib.figure.Figure, "savefig", spy_savefig)

        corr = np.eye(3)
        out = tmp_path / "triangle_border.pdf"
        _generate_correlation_heatmap(corr, ["a", "b", "c"], out, triangle="lower")

        assert state["spines_visible"] == [False, False, False, False]
        assert len(state["lines"]) >= 3
        assert all(line[0][0] == line[1][0] or line[0][1] == line[1][1] for line in state["lines"])

    def test_title_uses_axes_title_without_suptitle_gap(self, tmp_path, monkeypatch):
        titles = {}

        import matplotlib.figure
        original_suptitle = matplotlib.figure.Figure.suptitle

        def spy_suptitle(self, *args, **kwargs):
            titles["suptitle_called"] = True
            return original_suptitle(self, *args, **kwargs)

        def fake_savefig(self, *args, **kwargs):
            titles["axes_title"] = self.axes[0].get_title()

        monkeypatch.setattr(matplotlib.figure.Figure, "suptitle", spy_suptitle)
        monkeypatch.setattr(matplotlib.figure.Figure, "savefig", fake_savefig)

        corr = np.eye(3)
        _generate_correlation_heatmap(corr, ["a", "b", "c"], tmp_path / "title.pdf", title="My Title")

        assert titles.get("suptitle_called") is None
        assert titles["axes_title"] == "My Title"

    def test_label_angle_controls_bottom_axis_rotation(self, tmp_path, monkeypatch):
        rotations = []

        import matplotlib.figure
        original_savefig = matplotlib.figure.Figure.savefig

        def spy_savefig(self, *args, **kwargs):
            rotations.extend(label.get_rotation() for label in self.axes[0].get_xticklabels())
            return original_savefig(self, *args, **kwargs)

        monkeypatch.setattr(matplotlib.figure.Figure, "savefig", spy_savefig)

        corr = np.eye(3)
        _generate_correlation_heatmap(corr, ["a", "b", "c"], tmp_path / "angle.pdf", label_angle=60)

        assert rotations
        assert set(rotations) == {60.0}

    def test_upper_triangle_places_colorbar_on_left(self, tmp_path, monkeypatch):
        positions = []

        import matplotlib.figure
        original_savefig = matplotlib.figure.Figure.savefig

        def spy_savefig(self, *args, **kwargs):
            heatmap_ax, colorbar_ax = self.axes[0], self.axes[1]
            positions.append((heatmap_ax.get_position().x0, colorbar_ax.get_position().x0))
            return original_savefig(self, *args, **kwargs)

        monkeypatch.setattr(matplotlib.figure.Figure, "savefig", spy_savefig)

        corr = np.eye(3)
        _generate_correlation_heatmap(corr, ["a", "b", "c"], tmp_path / "upper_cbar.pdf", triangle="upper")

        heatmap_left, colorbar_left = positions[-1]
        assert colorbar_left < heatmap_left

    def test_lower_and_full_place_colorbar_on_right(self, tmp_path, monkeypatch):
        positions = []

        import matplotlib.figure
        original_savefig = matplotlib.figure.Figure.savefig

        def spy_savefig(self, *args, **kwargs):
            heatmap_ax, colorbar_ax = self.axes[0], self.axes[1]
            positions.append((heatmap_ax.get_position().x1, colorbar_ax.get_position().x0))
            return original_savefig(self, *args, **kwargs)

        monkeypatch.setattr(matplotlib.figure.Figure, "savefig", spy_savefig)

        corr = np.eye(3)
        _generate_correlation_heatmap(corr, ["a", "b", "c"], tmp_path / "lower_cbar.pdf", triangle="lower")
        _generate_correlation_heatmap(corr, ["a", "b", "c"], tmp_path / "full_cbar.pdf", triangle="full")

        assert all(colorbar_left > heatmap_right for heatmap_right, colorbar_left in positions)



class TestWriteCorrelationCsv:
    def test_csv_output(self, tmp_path):
        corr = np.array([[1.0, 0.5], [0.5, 1.0]])
        out = tmp_path / "corr.csv"
        _write_correlation_csv(corr, ["a", "b"], out)
        with open(out) as f:
            reader = csv.reader(f)
            header = next(reader)
            assert header == ["", "a", "b"]
            rows = list(reader)
            assert rows[0][0] == "a"


# ---------------------------------------------------------------------------
# Task 9: CLI integration tests
# ---------------------------------------------------------------------------


@pytest.fixture
def e2e_data(tmp_path):
    msa_dir = tmp_path / "msa"
    tree_dir = tmp_path / "trees"
    msa_dir.mkdir()
    tree_dir.mkdir()

    _write_fasta(msa_dir / "gene1.fa", _AA_SEQS, ["A", "B", "C"])
    _write_fasta(msa_dir / "gene2.fa", _AA_SEQS, ["A", "B", "C"])
    _write_newick(tree_dir / "gene1.tre", "(A,B,C);")
    _write_newick(tree_dir / "gene2.tre", "(A,B,C);")

    return tmp_path, msa_dir, tree_dir


class TestIntegrationCLI:
    def test_metrics_csv_output(self, e2e_data):
        tmp_path, msa_dir, tree_dir = e2e_data
        out = tmp_path / "output"
        result = run_metrics(
            msa_dir=msa_dir,
            tree_dir=tree_dir,
            output_dir=out,
            skip_freq=True,
            skip_pairwise_identity=True,
            quiet=True,
        )
        assert result["status"] == "success"
        assert (out / "metrics.csv").exists()
        assert (out / "result.json").exists()
        assert (out / "metrics.log").exists()

    def test_overwrite_flag(self, e2e_data):
        tmp_path, msa_dir, tree_dir = e2e_data
        out = tmp_path / "output"
        out.mkdir()
        (out / "old.txt").write_text("old")

        result = run_metrics(
            msa_dir=msa_dir,
            tree_dir=tree_dir,
            output_dir=out,
            overwrite=True,
            skip_freq=True,
            skip_pairwise_identity=True,
            quiet=True,
        )
        assert result["status"] == "success"

    def test_tree_only_mode(self, e2e_data):
        tmp_path, _, tree_dir = e2e_data
        out = tmp_path / "output_tree"
        result = run_metrics(
            tree_dir=tree_dir,
            output_dir=out,
            skip_freq=True,
            skip_pairwise_identity=True,
            quiet=True,
        )
        assert result["status"] == "success"
        with open(out / "metrics.csv") as f:
            reader = csv.DictReader(f)
            row = next(reader)
            assert row["DataType"] == ""

    def test_skip_freq_statistics(self, e2e_data):
        tmp_path, msa_dir, _ = e2e_data
        out = tmp_path / "output_nofreq"
        result = run_metrics(
            msa_dir=msa_dir,
            output_dir=out,
            skip_freq=True,
            skip_pairwise_identity=True,
            quiet=True,
        )
        with open(out / "metrics.csv") as f:
            reader = csv.DictReader(f)
            fieldnames = reader.fieldnames
            assert "freqA" not in fieldnames

    def test_with_freq_statistics(self, e2e_data):
        tmp_path, msa_dir, _ = e2e_data
        out = tmp_path / "output_freq"
        result = run_metrics(
            msa_dir=msa_dir,
            output_dir=out,
            skip_freq=False,
            skip_pairwise_identity=True,
            quiet=True,
        )
        with open(out / "metrics.csv") as f:
            reader = csv.DictReader(f)
            fieldnames = reader.fieldnames
            assert "freqA" in fieldnames

    def test_taxon_mismatch_in_pipeline(self, e2e_data):
        tmp_path, msa_dir, tree_dir = e2e_data
        _write_newick(tree_dir / "gene1.tre", "(A,B,D);")
        out = tmp_path / "output_mismatch"
        result = run_metrics(
            msa_dir=msa_dir,
            tree_dir=tree_dir,
            output_dir=out,
            skip_freq=True,
            skip_pairwise_identity=True,
            quiet=True,
        )
        assert result["status"] == "partial" or result["status"] == "success"

    def test_dry_run_no_files(self, e2e_data):
        tmp_path, msa_dir, _ = e2e_data
        out = tmp_path / "output_dry"
        result = run_metrics(
            msa_dir=msa_dir,
            output_dir=out,
            dry_run=True,
            skip_pairwise_identity=True,
            quiet=True,
        )
        assert not out.exists()
