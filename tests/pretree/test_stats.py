from __future__ import annotations

import json
from pathlib import Path

import pytest
from click.testing import CliRunner

from phyloai.cli.main import cli


TEST_DATA = Path("ref/phylogenomics_examples")


def test_pretree_group_is_registered() -> None:
    runner = CliRunner()

    result = runner.invoke(cli, ["pretree", "--help"])

    assert result.exit_code == 0
    assert "Pre-tree" in result.output


def test_detect_seq_type_aa() -> None:
    from phyloai.pretree.stats import detect_seq_type

    assert detect_seq_type(["MTEYKLVVVG", "ACDEFGHIKL"]) == "AA"


def test_detect_seq_type_nt() -> None:
    from phyloai.pretree.stats import detect_seq_type

    assert detect_seq_type(["ACGTACGT", "AUGCUUAA"]) == "NT"


def test_detect_seq_type_nt_with_iupac_ambiguity() -> None:
    from phyloai.pretree.stats import detect_seq_type

    assert detect_seq_type(["ACGT", "ACGR", "ACGN"]) == "NT"


def test_detect_seq_type_ambiguous_falls_back_to_aa() -> None:
    from phyloai.pretree.stats import detect_seq_type

    assert detect_seq_type(["ACGT", "ACGX"]) == "AA"


@pytest.mark.parametrize(
    ("char", "seq_type", "expected"),
    [
        ("A", "AA", "standard"),
        ("-", "AA", "gap"),
        ("B", "AA", "ambiguous"),
        ("X", "AA", "ambiguous"),
        ("N", "NT", "gap"),
    ],
)
def test_classify_char(char: str, seq_type: str, expected: str) -> None:
    from phyloai.pretree.stats import classify_char

    assert classify_char(char, seq_type) == expected


def test_normalize_pattern_char_treats_question_mark_as_gap() -> None:
    from phyloai.pretree.stats import normalize_pattern_char

    assert normalize_pattern_char("?") == "-"
    assert normalize_pattern_char("A") == "A"


def test_stop_codon_warning() -> None:
    from phyloai.pretree.stats import check_stop_codons

    warnings = check_stop_codons(["MTE*KL", "ACDEFG"], "example.faa")

    assert warnings == [
        "[WARN] Stop codon (*) found in example.faa. This may indicate upstream processing errors."
    ]


def test_compute_site_patterns_basic() -> None:
    from phyloai.pretree.stats import compute_site_patterns

    stats = compute_site_patterns(["AAAAA", "AAAAC", "AACCC", "AACCC"], "AA")

    assert stats["alignment_length"] == 5
    assert stats["constant_sites"] == {"count": 2, "ratio": 0.4}
    assert stats["distinct_patterns"] == {"count": 3, "ratio": 0.6}
    assert stats["parsimony_informative"] == {"count": 2, "ratio": 0.4}
    assert stats["singleton_sites"] == {"count": 1, "ratio": 0.2}


def test_stats_single_file_aligned() -> None:
    from phyloai.pretree.stats import stats_single_file

    stats = stats_single_file(TEST_DATA / "test" / "EOG090X0971.faa")

    assert stats["format"] == "fasta"
    assert stats["seq_type"] == "AA"
    assert stats["is_aligned"] is True
    assert stats["n_taxa"] == 6
    assert stats["alignment_length"] == 1042
    assert stats["constant_sites"]["count"] == 482
    assert stats["distinct_patterns"]["count"] == 624
    assert stats["parsimony_informative"]["count"] == 87
    assert stats["singleton_sites"]["count"] == 473
    assert stats["gap_ambiguous_ratio"] == 0.409949
    assert len(stats["per_taxon"]) == 6


def test_stats_single_file_unaligned() -> None:
    from phyloai.pretree.stats import stats_single_file

    stats = stats_single_file(TEST_DATA / "2-loci_filter" / "faa" / "EOG090X0007.faa")

    assert stats["format"] == "fasta"
    assert stats["seq_type"] == "AA"
    assert stats["is_aligned"] is False
    assert stats["n_taxa"] == 5
    assert stats["seq_length"]["median"] > 0
    assert stats["total_length"] > 0
    assert stats["gap_ambiguous_ratio"] == 0.0
    assert len(stats["per_taxon"]) == 5


def test_stats_single_file_raw_alignment_matches_iqtree_patterns() -> None:
    from phyloai.pretree.stats import stats_single_file

    stats = stats_single_file(TEST_DATA / "test" / "raw.fa")

    assert stats["is_aligned"] is True
    assert stats["distinct_patterns"]["count"] == 1053473


def test_stats_single_file_reports_seq_type_fallback_warning(tmp_path: Path) -> None:
    from phyloai.pretree.stats import stats_single_file

    path = tmp_path / "ambiguous.fa"
    path.write_text(">tax1\nACGTX\n>tax2\nACGTX\n")

    stats = stats_single_file(path)

    assert stats["seq_type"] == "AA"
    assert "Cannot determine seq_type" in stats["warnings"][0]


def test_collect_seq_files() -> None:
    from phyloai.pretree.stats import collect_seq_files

    files = collect_seq_files(TEST_DATA / "test")

    assert [path.name for path in files] == ["EOG090X0971.faa", "EOG090X0971.fna", "raw.fa"]


def test_stats_directory_serial() -> None:
    from phyloai.pretree.stats import aggregate_summary, stats_directory

    results, warnings = stats_directory(TEST_DATA / "test", seq_type=None, input_format=None, threads=1)
    summary = aggregate_summary(results)

    assert len(results) == 3
    assert warnings == []
    assert summary["n_genes"] == 3
    assert summary["n_genes_ok"] == 3
    assert summary["is_aligned"] is True
    assert summary["seq_type"] == "mixed"


def test_stats_directory_parallel_matches_serial() -> None:
    from phyloai.pretree.stats import stats_directory

    serial_results, serial_warnings = stats_directory(TEST_DATA / "test", None, None, threads=1)
    parallel_results, parallel_warnings = stats_directory(TEST_DATA / "test", None, None, threads=2)

    assert parallel_results == serial_results
    assert parallel_warnings == serial_warnings


def test_stats_directory_reports_progress_for_each_file() -> None:
    from phyloai.pretree.stats import collect_seq_files, stats_directory

    seen: list[str] = []

    results, warnings = stats_directory(
        TEST_DATA / "test",
        seq_type=None,
        input_format=None,
        threads=2,
        progress_callback=lambda path: seen.append(path.name),
    )

    assert len(results) == 3
    assert warnings == []
    assert sorted(seen) == sorted(path.name for path in collect_seq_files(TEST_DATA / "test"))


def test_cli_stats_json_output(tmp_path: Path) -> None:
    """Stats writes result.json to output-dir; terminal still shows Rich panels."""
    runner = CliRunner()
    output_dir = tmp_path / "stats_out"

    result = runner.invoke(
        cli,
        [
            "pretree",
            "stats",
            "--seq",
            str(TEST_DATA / "test" / "EOG090X0971.faa"),
            "--output-dir",
            str(output_dir),
        ],
    )
    assert result.exit_code == 0
    result_path = output_dir / "result.json"
    assert result_path.exists()
    saved = json.loads(result_path.read_text())
    assert saved["status"] == "success"
    assert saved["data"]["alignment_length"] == 1042
    assert "constant_sites" in saved["data"]


def test_cli_stats_mutual_exclusivity_error() -> None:
    runner = CliRunner()

    result = runner.invoke(
        cli,
        [
            "pretree",
            "stats",
            "--seq",
            str(TEST_DATA / "test" / "EOG090X0971.faa"),
            "--seq-dir",
            str(TEST_DATA / "test"),
        ],
    )

    assert result.exit_code == 1
    assert "mutually exclusive" in result.output.lower()


def test_cli_stats_no_args_error() -> None:
    runner = CliRunner()

    result = runner.invoke(cli, ["pretree", "stats"])

    assert result.exit_code == 1
    assert "one of --seq or --seq-dir must be provided" in result.output.lower()


def test_cli_stats_missing_file_exits_one() -> None:
    runner = CliRunner()

    result = runner.invoke(cli, ["pretree", "stats", "--seq", "missing.faa"])

    assert result.exit_code == 1


def test_cli_stats_all_failed_directory_exits_two(tmp_path: Path) -> None:
    runner = CliRunner()
    bad = tmp_path / "bad.fa"
    bad.write_text(">broken\n")

    result = runner.invoke(cli, ["pretree", "stats", "--seq-dir", str(tmp_path), "--output-dir", str(tmp_path / "out")])

    assert result.exit_code == 2


def test_directory_output_writes_summary_and_per_gene(tmp_path: Path) -> None:
    """Directory mode writes result.json and per-gene CSV in output-dir."""
    runner = CliRunner()
    output_dir = tmp_path / "stats_out"

    result = runner.invoke(
        cli,
        [
            "pretree",
            "stats",
            "--seq-dir",
            str(TEST_DATA / "test"),
            "--per-gene",
            "--output-dir",
            str(output_dir),
        ],
    )

    assert result.exit_code == 0
    result_path = output_dir / "result.json"
    per_gene_path = output_dir / "per-gene.csv"
    assert result_path.exists()
    assert per_gene_path.exists()
    saved = json.loads(result_path.read_text())
    assert saved["data"]["summary"]["n_genes"] > 0
    per_gene_header = per_gene_path.read_text().splitlines()[0]
    assert "gene" in per_gene_header
    assert "n_taxa" in per_gene_header
    assert "alignment_length" in per_gene_header


def test_directory_table_format_can_write_tsv(tmp_path: Path) -> None:
    runner = CliRunner()
    output_dir = tmp_path / "stats_out"

    result = runner.invoke(
        cli,
        [
            "pretree",
            "stats",
            "--seq-dir",
            str(TEST_DATA / "test"),
            "--per-gene",
            "--table-format",
            "tsv",
            "--output-dir",
            str(output_dir),
        ],
    )

    assert result.exit_code == 0
    per_gene_path = output_dir / "per-gene.tsv"
    assert per_gene_path.exists()
    assert "\t" in per_gene_path.read_text().splitlines()[0]


def test_cli_help_explains_options() -> None:
    runner = CliRunner()

    result = runner.invoke(cli, ["pretree", "stats", "--help"])

    assert result.exit_code == 0
    assert "Inspect one sequence file or summarize a directory of sequence files." in result.output
    assert "Exactly" in result.output
    assert "--seq or --seq-dir" in result.output
    assert "is required" in result.output
    assert "Override auto-detection when a file suffix is" in result.output
    assert "misleading" in result.output
    assert "[fasta|phylip-relaxed|nexus]" in result.output
    assert "Directory where result.json" in result.output
    assert "--table-format" in result.output
    assert "--per-gene-format" not in result.output


def test_cli_rejects_phylip_input_format_choice() -> None:
    runner = CliRunner()

    result = runner.invoke(
        cli,
        [
            "pretree",
            "stats",
            "--seq",
            str(TEST_DATA / "test" / "EOG090X0971.faa"),
            "--input-format",
            "phylip",
        ],
    )

    assert result.exit_code != 0
    assert "invalid value for '--input-format'" in result.output.lower()


def test_directory_per_gene_writes_to_output_dir(tmp_path: Path) -> None:
    """Per-gene data goes to per-gene.csv in output-dir."""
    runner = CliRunner()
    output_dir = tmp_path / "stats_out"

    result = runner.invoke(
        cli,
        [
            "pretree",
            "stats",
            "--seq-dir",
            str(TEST_DATA / "test"),
            "--per-gene",
            "--output-dir",
            str(output_dir),
        ],
    )

    assert result.exit_code == 0
    result_path = output_dir / "result.json"
    per_gene_path = output_dir / "per-gene.csv"
    assert result_path.exists()
    assert per_gene_path.exists()
    saved = json.loads(result_path.read_text())
    assert "per_gene" in saved["data"]
    per_gene_content = per_gene_path.read_text()
    per_gene_header = per_gene_content.splitlines()[0]
    assert "gene" in per_gene_header
    assert "n_taxa" in per_gene_header
    assert "EOG090X0971" in per_gene_content


def test_cli_stats_rejects_per_gene_in_single_file_mode() -> None:
    runner = CliRunner()

    result = runner.invoke(
        cli,
        [
            "pretree",
            "stats",
            "--seq",
            str(TEST_DATA / "test" / "EOG090X0971.faa"),
            "--per-gene",
        ],
    )

    assert result.exit_code == 1
    assert "--per-gene is directory mode only" in result.output


def test_single_file_output_includes_per_taxon_data(tmp_path: Path) -> None:
    """Default JSON output for single file includes per_taxon array."""
    runner = CliRunner()
    output_dir = tmp_path / "stats_out"

    result = runner.invoke(
        cli,
        [
            "pretree",
            "stats",
            "--seq",
            str(TEST_DATA / "2-loci_filter" / "fna" / "EOG090X0971.fna"),
            "--output-dir",
            str(output_dir),
        ],
    )

    assert result.exit_code == 0
    result_path = output_dir / "result.json"
    assert result_path.exists()
    saved = json.loads(result_path.read_text())
    assert "per_taxon" in saved["data"]
    assert len(saved["data"]["per_taxon"]) > 0
