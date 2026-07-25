"""Tests for phyloai.posttree.signal helpers."""
from __future__ import annotations

import csv as csv_mod
from pathlib import Path
from types import SimpleNamespace

import pytest

from phyloai.posttree.signal import (
    _compare_groups,
    _delta_score,
    _outlier_loci,
    _parse_partition_ranges,
    _parse_sitelh,
    _sum_gene_lnl,
    _support_label,
)


class TestParsePartitionRanges:
    def test_basic(self, tmp_path: Path) -> None:
        path = tmp_path / "p.txt"
        path.write_text("LG, geneA = 1-235\nLG, geneB = 236-461\n")

        assert _parse_partition_ranges(path) == [
            {"locus": "geneA", "start": 1, "end": 235},
            {"locus": "geneB", "start": 236, "end": 461},
        ]

    def test_empty_raises(self, tmp_path: Path) -> None:
        path = tmp_path / "empty.txt"
        path.write_text("")

        with pytest.raises(ValueError, match="empty"):
            _parse_partition_ranges(path)

    def test_bad_line_raises(self, tmp_path: Path) -> None:
        path = tmp_path / "bad.txt"
        path.write_text("not a partition line\n")

        with pytest.raises(ValueError, match="Unparseable"):
            _parse_partition_ranges(path)

    def test_nexus_charset(self, tmp_path: Path) -> None:
        path = tmp_path / "charset.nex"
        path.write_text("#NEXUS\nbegin sets;\n  charset geneA = 1-235;\n  charset geneB = 236-461;\nend;\n")

        assert _parse_partition_ranges(path) == [
            {"locus": "geneA", "start": 1, "end": 235},
            {"locus": "geneB", "start": 236, "end": 461},
        ]

    def test_nexus_uppercase(self, tmp_path: Path) -> None:
        path = tmp_path / "upper.nex"
        path.write_text("#NEXUS\nBEGIN SETS;\n  CHARSET gene = 1-10;\nEND;\n")

        assert _parse_partition_ranges(path) == [
            {"locus": "gene", "start": 1, "end": 10},
        ]

    def test_nexus_skips_charpartition(self, tmp_path: Path) -> None:
        path = tmp_path / "model.nex"
        path.write_text("#NEXUS\nbegin sets;\n  charset geneA = 1-100;\n  charpartition mymodels = LG+I+G: geneA;\nend;\n")

        assert _parse_partition_ranges(path) == [
            {"locus": "geneA", "start": 1, "end": 100},
        ]

    def test_nexus_no_charset_raises(self, tmp_path: Path) -> None:
        path = tmp_path / "nocharset.nex"
        path.write_text("#NEXUS\nbegin sets;\n  charpartition mymodels = LG+I+G: geneA;\nend;\n")

        with pytest.raises(ValueError, match="empty"):
            _parse_partition_ranges(path)


class TestParseSitelh:
    def test_two_trees(self, tmp_path: Path) -> None:
        path = tmp_path / "site.sitelh"
        path.write_text("2 3\nTree1 -1.0 -2.0 -3.0\nTree2 -1.5 -2.5 -3.5\n")

        labels, scores = _parse_sitelh(path)

        assert labels == ["Tree1", "Tree2"]
        assert scores[0] == pytest.approx([-1.0, -2.0, -3.0])
        assert scores[1] == pytest.approx([-1.5, -2.5, -3.5])

    def test_wrong_column_count_raises(self, tmp_path: Path) -> None:
        path = tmp_path / "bad.sitelh"
        path.write_text("2 3\nTree1 -1.0 -2.0\nTree2 -1.5 -2.5 -3.5\n")

        with pytest.raises(ValueError, match="columns"):
            _parse_sitelh(path)


class TestSumGeneLnl:
    def test_sums_inclusive_gene_range(self) -> None:
        assert _sum_gene_lnl([[-1.0, -2.0, -3.0], [-4.0, -5.0, -6.0]], 2, 3) == [-5.0, -11.0]


class TestDeltaScore:
    def test_two_trees_delta_is_signed(self) -> None:
        assert _delta_score([-5.0, -8.0]) == pytest.approx(3.0)

    def test_three_trees_delta_is_mean_pairwise_difference(self) -> None:
        assert _delta_score([-1.0, -2.0, -3.0]) == pytest.approx(4 / 3)


class TestSupportLabel:
    def test_support_label_uses_tree_number(self) -> None:
        assert _support_label([-5.0, -8.0], ["Tree1", "Tree2"]) == "Tree1"

    def test_support_label_marks_ties_ambiguous(self) -> None:
        assert _support_label([-5.0, -5.0 + 5e-10], ["Tree1", "Tree2"]) == "ambiguous"


class TestOutlierLoci:
    def test_no_outliers_all_same(self) -> None:
        assert not any(_outlier_loci([1.0] * 10))

    def test_extreme_value_is_outlier(self) -> None:
        flags = _outlier_loci([1.0] * 9 + [100.0])

        assert flags[-1] is True
        assert not any(flags[:-1])


class TestCompareGroups:
    def test_writes_comparison_files(self, tmp_path: Path) -> None:
        metrics = tmp_path / "metrics.csv"
        metrics.write_text("loci,entropy\na,1.0\nb,3.0\n")

        csv_path, pdf_path, sig_info = _compare_groups(["a"], ["b"], metrics, "supported", "other", tmp_path)

        assert csv_path == tmp_path / "supported_comparison.csv"
        assert pdf_path == tmp_path / "supported_comparison.pdf"
        assert csv_path.read_text().splitlines() == [
            "metric,supported_mean,supported_n,other_mean,other_n,wilcoxon_p",
            "entropy,1.0,1,3.0,1,1.0",
        ]
        assert pdf_path.is_file()

    def test_empty_group_a_writes_na_csv_and_explanatory_pdf(self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
        metrics = tmp_path / "metrics.csv"
        metrics.write_text("loci,entropy\na,1.0\nb,3.0\n")
        from matplotlib.axes import Axes

        text_calls: list[str] = []
        original_text = Axes.text

        def record_text(axis, x, y, text, *args, **kwargs):
            text_calls.append(text)
            return original_text(axis, x, y, text, *args, **kwargs)

        monkeypatch.setattr(Axes, "boxplot", lambda *args, **kwargs: pytest.fail("empty group must not create boxplots"))
        monkeypatch.setattr(Axes, "text", record_text)

        csv_path, pdf_path, sig_info = _compare_groups([], ["a", "b"], metrics, "supported", "other", tmp_path)

        assert csv_path.read_text().splitlines() == [
            "metric,supported_mean,supported_n,other_mean,other_n,wilcoxon_p",
            "entropy,NA,0,2.0,2,NA",
        ]
        assert pdf_path.is_file()
        assert pdf_path.stat().st_size > 0
        assert text_calls == ["No supported loci"]


class TestRunSignalLnlDryRun:
    def test_dry_run_returns_cmd(self, tmp_path: Path) -> None:
        from phyloai.posttree.signal import run_signal_lnl

        matrix = tmp_path / "m.fa"
        matrix.write_text(">A\nMKT\n>B\nMKA\n")
        t1 = tmp_path / "T1.nwk"
        t1.write_text("(A,B);\n")
        t2 = tmp_path / "T2.nwk"
        t2.write_text("(B,A);\n")

        result = run_signal_lnl(
            matrix=matrix,
            candidate_trees=[t1, t2],
            model_expr="LG+F+R4",
            partitions=None,
            locus_ranges=None,
            guide_tree=None,
            threads="auto",
            iqtree_path=None,
            tool_args=None,
            metrics=None,
            output_dir=tmp_path / "out",
            overwrite=False,
            dry_run=True,
            quiet=True,
        )

        assert result["status"] == "success"
        assert "iqtree3" in result["data"]["cmd"][0]
        assert "-wslr" in result["data"]["cmd"]
        assert result["params"]["model_expr"] == "LG+F+R4"
        assert result["params"]["partitions"] is None
        assert result["params"]["locus_ranges"] is None


class TestRunSignalLnlReviewFindings:
    @staticmethod
    def _run(
        tmp_path: Path,
        monkeypatch: pytest.MonkeyPatch,
        *,
        partitions: Path | None = None,
        locus_ranges: Path | None = None,
        metrics: Path | None = None,
        sitelh: str = "2 2\nT1 -1 -1\nT2 -1 -1\n",
    ) -> tuple[dict, Path]:
        from phyloai.posttree import signal

        matrix = tmp_path / "matrix.fa"
        trees = tmp_path / "trees.nwk"
        matrix.write_text(">A\nAA\n>B\nAA\n")
        trees.write_text("(A,B);\n")
        output_dir = tmp_path / "out"

        def succeed(*args, **kwargs):
            work_dir = Path(kwargs.get("cwd", output_dir))
            work_dir.mkdir(parents=True, exist_ok=True)
            (work_dir / "lnl.sitelh").write_text(sitelh)
            return SimpleNamespace(returncode=0, stdout="", stderr="")

        monkeypatch.setattr(signal, "_resolve_iqtree_path", lambda *args: "iqtree3")
        monkeypatch.setattr(signal, "_detect_iqtree_version", lambda *args: {"iqtree3": "test"})
        monkeypatch.setattr(signal.subprocess, "run", succeed)
        return signal.run_signal_lnl(
            matrix=matrix, candidate_trees=[trees], model_expr="LG", partitions=partitions,
            locus_ranges=locus_ranges, guide_tree=None, threads="1", iqtree_path=None,
            tool_args=None, metrics=metrics, output_dir=output_dir, overwrite=False,
            dry_run=False, quiet=True,
        ), output_dir

    def test_metrics_comparison_is_written_when_no_loci_are_outliers(self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
        ranges = tmp_path / "ranges.txt"
        metrics = tmp_path / "metrics.csv"
        ranges.write_text("LG, a = 1-1\nLG, b = 2-2\n")
        metrics.write_text("loci,entropy\na,1\nb,2\n")

        result, output_dir = self._run(tmp_path, monkeypatch, locus_ranges=ranges, metrics=metrics)

        assert result["status"] == "success"
        assert (output_dir / "outlier_comparison.csv").is_file()

    def test_support_group_pairwise_comparison_written(self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
        ranges = tmp_path / "ranges.txt"
        metrics = tmp_path / "metrics.csv"
        ranges.write_text("LG, a = 1-2\nLG, b = 3-4\n")
        metrics.write_text("loci,entropy\na,1\nb,2\n")
        sitelh = "2 4\nT1 -1.0 -1.1 -10.0 -10.1\nT2 -10.0 -10.1 -1.0 -1.1\n"

        result, output_dir = self._run(tmp_path, monkeypatch, locus_ranges=ranges, metrics=metrics, sitelh=sitelh)

        assert result["status"] == "success"
        assert (output_dir / "outlier_comparison.csv").is_file()
        assert (output_dir / "support_comparison.csv").is_file()
        assert (output_dir / "support_comparison.pdf").is_file()
        assert "support_comparison" in result["data"]["output_files"]
        assert "support_comparison_plot" in result["data"]["output_files"]

    def test_support_group_three_trees_merged_output_and_ambiguous_excluded(self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
        ranges = tmp_path / "ranges.txt"
        metrics = tmp_path / "metrics.csv"
        ranges.write_text("LG, a = 1-1\nLG, b = 2-2\nLG, c = 3-3\nLG, d = 4-4\n")
        metrics.write_text("loci,entropy\na,1\nb,2\nc,3\nd,0\n")
        sitelh = "3 4\nT1 -1.0 -10.0 -10.0 -1.0\nT2 -10.0 -1.0 -10.0 -1.0\nT3 -10.0 -10.0 -1.0 -1.0\n"

        result, output_dir = self._run(tmp_path, monkeypatch, locus_ranges=ranges, metrics=metrics, sitelh=sitelh)

        assert result["status"] == "success"
        assert (output_dir / "support_comparison.csv").is_file()
        assert (output_dir / "support_comparison.pdf").is_file()
        assert "support_comparison" in result["data"]["output_files"]
        assert "support_comparison_plot" in result["data"]["output_files"]
        with open(output_dir / "gene_lnl.csv") as fh:
            gene_rows = list(csv_mod.DictReader(fh))
        d_row = next(r for r in gene_rows if r["locus"] == "d")
        assert d_row["support"] == "ambiguous"
        with open(output_dir / "support_comparison.csv") as fh:
            comp = fh.read()
            assert "T1_mean" in comp
            assert "T2_mean" in comp
            assert "T3_mean" in comp
            assert "T1_vs_T2_wilcoxon_p" in comp
            assert "T1_vs_T3_wilcoxon_p" in comp
            assert "T2_vs_T3_wilcoxon_p" in comp

    def test_multi_group_preserves_insertion_order_not_lexicographic(self, tmp_path: Path) -> None:
        from phyloai.posttree.signal import _compare_multiple_groups

        metrics = tmp_path / "metrics.csv"
        metrics.write_text("loci,entropy,GC\na,1.0,0.4\nb,3.0,0.5\nc,2.0,0.3\n")

        groups = {"T2": ["b"], "T10": ["a"], "T3": ["c"]}
        csv_path, pdf_path, _ = _compare_multiple_groups(groups, metrics, tmp_path, prefix="order_test")

        with open(csv_path) as fh:
            header = fh.readline()
        cols = header.strip().split(",")
        assert cols[0] == "metric"
        assert cols[1] == "T2_mean"
        assert cols[2] == "T2_n"
        assert cols[3] == "T10_mean"
        assert cols[4] == "T10_n"
        assert cols[5] == "T3_mean"
        assert cols[6] == "T3_n"
        assert "T2_vs_T10_wilcoxon_p" in cols
        assert "T2_vs_T3_wilcoxon_p" in cols
        assert "T10_vs_T3_wilcoxon_p" in cols
        assert cols.index("T2_mean") < cols.index("T10_mean") < cols.index("T3_mean")

    def test_range_past_sitelh_site_count_returns_structured_error(self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
        ranges = tmp_path / "ranges.txt"
        ranges.write_text("LG, a = 1-3\n")

        result, output_dir = self._run(tmp_path, monkeypatch, locus_ranges=ranges)

        assert result["status"] == "error"
        assert "site count" in result["error"]
        assert (output_dir / "result.json").is_file()

    @pytest.mark.parametrize(
        ("sitelh", "ranges", "metrics"),
        [
            ("", None, None),
            ("2 2\nT1 -1 -1\nT2 -1 -1\n", "bad boundary\n", None),
            ("2 2\nT1 -1 -1\nT2 -1 -1\n", "LG, a = 1-1\nLG, b = 2-2\n", "loci,entropy\na,1\n"),
        ],
        ids=["missing-sitelh", "invalid-boundary", "invalid-metrics"],
    )
    def test_post_iqtree_output_failures_write_structured_error(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch, sitelh: str, ranges: str | None, metrics: str | None,
    ) -> None:
        range_path = tmp_path / "ranges.txt" if ranges else None
        metrics_path = tmp_path / "metrics.csv" if metrics else None
        if range_path:
            range_path.write_text(ranges)
        if metrics_path:
            metrics_path.write_text(metrics)

        result, output_dir = self._run(tmp_path, monkeypatch, locus_ranges=range_path, metrics=metrics_path, sitelh=sitelh)

        assert result["status"] == "error"
        assert result["error_category"] == "output"
        assert (output_dir / "result.json").is_file()

    @pytest.mark.parametrize("option", ["partitions", "locus_ranges", "metrics"])
    def test_optional_files_must_be_regular_readable_files(self, tmp_path: Path, option: str) -> None:
        from phyloai.posttree.signal import run_signal_lnl

        matrix = tmp_path / "matrix.fa"
        trees = tmp_path / "trees.nwk"
        invalid = tmp_path / option
        matrix.write_text(">A\nAA\n>B\nAA\n")
        trees.write_text("(A,B);\n")
        invalid.mkdir()
        kwargs = {"partitions": None, "locus_ranges": None, "metrics": None, option: invalid}

        result = run_signal_lnl(
            matrix=matrix, candidate_trees=[trees], model_expr=None, guide_tree=None, threads="1",
            iqtree_path=None, tool_args="-m LG", output_dir=tmp_path / "out", overwrite=False,
            dry_run=True, quiet=True, **kwargs,
        )

        assert result["status"] == "error"
        assert f"--{option.replace('_', '-')}" in result["error"]

    def test_tool_args_cannot_override_lnl_prefix(self, tmp_path: Path) -> None:
        from phyloai.posttree.signal import run_signal_lnl

        matrix = tmp_path / "matrix.fa"
        trees = tmp_path / "trees.nwk"
        matrix.write_text(">A\nAA\n>B\nAA\n")
        trees.write_text("(A,B);\n")

        result = run_signal_lnl(
            matrix=matrix, candidate_trees=[trees], model_expr="LG", partitions=None,
            locus_ranges=None, guide_tree=None, threads="1", iqtree_path=None,
            tool_args="--prefix custom", metrics=None, output_dir=tmp_path / "out",
            overwrite=False, dry_run=True, quiet=True,
        )

        assert result["status"] == "error"
        assert "Blocked flag in --tool-args: --prefix" in result["error"]

    def test_empty_candidate_trees_are_rejected_before_output_directory_handling(self, tmp_path: Path) -> None:
        from phyloai.posttree.signal import run_signal_lnl

        matrix = tmp_path / "matrix.fa"
        output_dir = tmp_path / "out"
        matrix.write_text(">A\nAA\n>B\nAA\n")
        output_dir.mkdir()
        (output_dir / "existing").write_text("keep")

        result = run_signal_lnl(
            matrix=matrix, candidate_trees=[], model_expr="LG", partitions=None,
            locus_ranges=None, guide_tree=None, threads="1", iqtree_path=None,
            tool_args=None, metrics=None, output_dir=output_dir, overwrite=False,
            dry_run=False, quiet=True,
        )

        assert result["status"] == "error"
        assert "--candidate-trees must not be empty" in result["error"]

    def test_guide_tree_must_be_a_readable_regular_file(self, tmp_path: Path) -> None:
        from phyloai.posttree.signal import run_signal_lnl

        matrix = tmp_path / "matrix.fa"
        trees = tmp_path / "trees.nwk"
        guide_tree = tmp_path / "guide"
        matrix.write_text(">A\nAA\n>B\nAA\n")
        trees.write_text("(A,B);\n")
        guide_tree.mkdir()

        result = run_signal_lnl(
            matrix=matrix, candidate_trees=[trees], model_expr="LG", partitions=None,
            locus_ranges=None, guide_tree=guide_tree, threads="1", iqtree_path=None,
            tool_args=None, metrics=None, output_dir=tmp_path / "out", overwrite=False,
            dry_run=True, quiet=True,
        )

        assert result["status"] == "error"
        assert "--guide-tree must be a readable regular file" in result["error"]


class TestRunSignalLnlWithFixture:
    """Uses runs/signal/ test data."""

    SIGNAL_DIR = Path("runs/signal")

    def test_site_lnl_csv_columns(self, tmp_path: Path) -> None:
        from phyloai.posttree.signal import run_signal_lnl

        if not self.SIGNAL_DIR.exists():
            pytest.skip("Signal test data not present")
        matrix = self.SIGNAL_DIR / "matrix.aa.fa"
        trees = self.SIGNAL_DIR / "trees"
        result = run_signal_lnl(
            matrix=matrix,
            candidate_trees=[trees],
            model_expr="LG+F+R4",
            partitions=None,
            locus_ranges=self.SIGNAL_DIR / "matrix.aa.partitions",
            guide_tree=None,
            threads="auto",
            iqtree_path=None,
            tool_args=None,
            metrics=None,
            output_dir=tmp_path / "lnl_out",
            overwrite=False,
            dry_run=False,
            quiet=True,
        )

        if result["status"] == "error" and result.get("error_category") == "env":
            pytest.skip("iqtree3 not available")
        assert result["status"] == "success"
        assert result["key_results"]["n_trees"] == 3
        assert result["key_results"]["n_sites"] == 5604
        assert result["key_results"]["n_loci"] == 20
        import csv as csv_mod

        with open(tmp_path / "lnl_out" / "site_lnl.csv") as handle:
            columns = csv_mod.DictReader(handle).fieldnames or []
        assert "ΔSLS" in columns
        assert "support" in columns
        assert "lnL_Tree1" in columns
