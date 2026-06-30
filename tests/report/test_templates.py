"""Tests for phyloai.report.templates."""
from __future__ import annotations

from phyloai.report.templates import (
    METHODS_GENERATORS,
    _safe_fmt,
    generate_all_methods,
    generate_methods_pretree_align,
)


class TestSafeFmt:
    def test_float(self):
        assert _safe_fmt(3.14159, ".2f") == "3.14"

    def test_int(self):
        assert _safe_fmt(1000, ",") == "1,000"

    def test_string_fallback(self):
        assert _safe_fmt("N/A", ".1f") == "?"
        assert _safe_fmt(None, ",") == "?"

    def test_percent(self):
        assert _safe_fmt(0.156, ".1%") == "15.6%"


class TestTemplatesRegistry:
    def test_has_all_implemented_steps(self):
        """Every step_id in STEP_ORDER MUST have a template generator."""
        from phyloai.report.collector import STEP_ORDER

        for step_id in STEP_ORDER:
            assert step_id in METHODS_GENERATORS, (
                f"Missing methods template for {step_id}. "
                f"Add generate_methods_{step_id.replace('.', '_')}() to templates.py."
            )
            gen = METHODS_GENERATORS[step_id]
            assert callable(gen), f"{step_id} generator is not callable"


class TestGenerateAllMethods:
    def test_dispatches(self):
        result = generate_all_methods(
            "pretree.convert",
            params={"to": "fasta", "input": "./raw"},
            key_results={"n_converted": 100},
            tool_versions={},
        )
        assert isinstance(result, str)
        assert len(result) > 0

    def test_unknown_step_returns_empty(self):
        result = generate_all_methods(
            "nonexistent.step",
            params={},
            key_results={},
            tool_versions={},
        )
        assert result == ""

    def test_failed_step_returns_empty(self):
        result = generate_all_methods(
            "pretree.align",
            params={},
            key_results={},
            tool_versions={},
            status="error",
        )
        assert result == ""


class TestGenerateMethodsPretreeAlign:
    def test_linsi(self):
        text = generate_methods_pretree_align(
            params={"method": "linsi", "seq_type": "AA", "backtrans": False},
            key_results={"n_aligned": 100, "n_skipped": 0, "mean_alignment_length": 500.0},
            tool_versions={"mafft": "7.526"},
        )
        assert "L-INS-i" in text
        assert "MAFFT" in text
        assert "7.526" in text

    def test_magus(self):
        text = generate_methods_pretree_align(
            params={"method": "magus", "seq_type": "AA", "backtrans": False},
            key_results={"n_aligned": 200, "n_skipped": 0, "mean_alignment_length": 300.0},
            tool_versions={"magus": "unknown version"},
        )
        assert "MAGUS" in text

    def test_backtrans(self):
        text = generate_methods_pretree_align(
            params={"method": "linsi", "seq_type": "AA", "backtrans": True},
            key_results={"n_aligned": 100, "n_skipped": 0, "mean_alignment_length": 500.0},
            tool_versions={"mafft": "7.526", "trimal": "1.4.1"},
        )
        assert "back-translation" in text.lower() or "codon-aware" in text.lower()

    def test_skipped_clause(self):
        text = generate_methods_pretree_align(
            params={"method": "auto", "seq_type": "NT", "backtrans": False},
            key_results={"n_aligned": 90, "n_skipped": 10, "mean_alignment_length": 200.0},
            tool_versions={"mafft": "7.526"},
        )
        assert "10" in text


class TestTaper:
    def test_basic(self):
        text = generate_all_methods(
            "pretree.filter.taper",
            params={"cutoff": 0.1},
            key_results={"n_input": 100, "n_retained": 95, "n_dropped": 5, "n_masked_sites": 1200},
            tool_versions={"taper": "1.0.0", "julia": "1.9"},
        )
        assert "TAPER" in text
        assert "1.0.0" in text
        assert "5" in text


class TestConcat:
    def test_basic(self):
        text = generate_all_methods(
            "pretree.concat",
            params={"taxa_occupancy": 0.75, "seq_type": "AA", "recoding": None},
            key_results={
                "n_msa_used": 187, "n_msa_dropped": 13, "n_taxa": 52,
                "total_length": 45000, "gap_ratio": 0.15, "pi_ratio": 0.42,
            },
            tool_versions={},
        )
        assert "75%" in text or "0.75" in text
        assert "187" in text
        assert "52" in text

    def test_recoding(self):
        text = generate_all_methods(
            "pretree.concat",
            params={"taxa_occupancy": 0.75, "seq_type": "AA", "recoding": "Dayhoff6"},
            key_results={
                "n_msa_used": 100, "n_msa_dropped": 0, "n_taxa": 50,
                "total_length": 30000, "gap_ratio": 0.1, "pi_ratio": 0.4,
            },
            tool_versions={},
        )
        assert "Dayhoff6" in text
        assert "recod" in text.lower()

    def test_jackknife(self):
        text = generate_all_methods(
            "pretree.concat.jackknife",
            params={"replicates": 100, "target_length": 50000, "seed": 42, "to": "fasta"},
            key_results={"n_replicates": 100, "min_length": 50012, "max_length": 53280, "mean_length": 51140.5},
            tool_versions={},
        )
        assert "gene-jackknife" in text
        assert "100 pseudoreplicates" in text
        assert "50,000" in text


class TestIqtree:
    def test_unpartitioned(self):
        text = generate_all_methods(
            "tree.ml.iqtree",
            params={"partitioned": False, "modelfinder": "MFP", "mset": "LG+C20+C60", "boot": 1000},
            key_results={"log_likelihood": -12345.67},
            tool_versions={"iqtree": "3.0.0"},
        )
        assert "IQ-TREE" in text
        assert "3.0.0" in text

    def test_partitioned(self):
        text = generate_all_methods(
            "tree.ml.iqtree",
            params={
                "partitioned": True, "merged_partitions": True, "rclusterf": 10,
                "modelfinder": "MFP", "mset": "LG+C20+C60", "boot": 1000,
            },
            key_results={"log_likelihood": -5000.0},
            tool_versions={"iqtree": "3.0.0"},
        )
        assert "partition" in text.lower()
        assert "rclusterf" in text.lower()
