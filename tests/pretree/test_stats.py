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


def test_cli_stats_json_output() -> None:
    """--output-format json with --output writes a JSON file; terminal still shows Rich panels."""
    runner = CliRunner()
    import tempfile, os
    with tempfile.NamedTemporaryFile(suffix=".json", delete=False) as f:
        out_path = Path(f.name)
    try:
        result = runner.invoke(
            cli,
            [
                "pretree",
                "stats",
                "--seq",
                str(TEST_DATA / "test" / "EOG090X0971.faa"),
                "--output",
                str(out_path),
                "--output-format",
                "json",
            ],
        )
        assert result.exit_code == 0
        saved = json.loads(out_path.read_text())
        assert saved["status"] == "success"
        assert saved["data"]["alignment_length"] == 1042
        assert "constant_sites" in saved["data"]
    finally:
        out_path.unlink(missing_ok=True)


def test_cli_stats_json_output_format_writes_json_file_ignoring_extension(tmp_path: Path) -> None:
    """--output-format json forces JSON in the output file regardless of file extension."""
    runner = CliRunner()
    output = tmp_path / "out.txt"

    result = runner.invoke(
        cli,
        [
            "pretree",
            "stats",
            "--seq",
            str(TEST_DATA / "test" / "EOG090X0971.faa"),
            "--output",
            str(output),
            "--output-format",
            "json",
        ],
    )

    assert result.exit_code == 0
    # file contains JSON despite .txt extension (--output-format overrides extension)
    saved_data = json.loads(output.read_text())
    assert saved_data["status"] == "success"
    assert saved_data["data"]["alignment_length"] == 1042


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


def test_cli_stats_output_any_extension_accepted(tmp_path: Path) -> None:
    """Output format is controlled by --output-format, not file extension; any path is valid."""
    runner = CliRunner()
    output = tmp_path / "stats.xyz"

    result = runner.invoke(
        cli,
        [
            "pretree",
            "stats",
            "--seq",
            str(TEST_DATA / "test" / "EOG090X0971.faa"),
            "--output",
            str(output),
        ],
    )

    assert result.exit_code == 0
    saved = json.loads(output.read_text())
    assert saved["status"] == "success"


def test_cli_stats_all_failed_directory_exits_two(tmp_path: Path) -> None:
    runner = CliRunner()
    bad = tmp_path / "bad.fa"
    bad.write_text(">broken\n")

    result = runner.invoke(cli, ["pretree", "stats", "--seq-dir", str(tmp_path)])

    assert result.exit_code == 2


def test_directory_output_writes_summary_and_adjacent_per_gene(tmp_path: Path) -> None:
    """Directory mode writes a JSON summary and an adjacent per-gene CSV."""
    runner = CliRunner()
    output = tmp_path / "stats.json"
    per_gene_output = tmp_path / "stats.per-gene.csv"

    result = runner.invoke(
        cli,
        [
            "pretree",
            "stats",
            "--seq-dir",
            str(TEST_DATA / "test"),
            "--per-gene",
            "--output",
            str(output),
        ],
    )

    assert result.exit_code == 0
    saved = json.loads(output.read_text())
    assert saved["data"]["summary"]["n_genes"] > 0
    per_gene_header = per_gene_output.read_text().splitlines()[0]
    assert per_gene_header == (
        "gene,n_taxa,n_taxa_ratio,length_type,alignment_length,"
        "gap_ratio,ambiguous_ratio,gap_ambiguous_ratio,missing_taxa,missing_taxa_ratio"
    )


def test_directory_per_gene_adjacent_csv_always_created(tmp_path: Path) -> None:
    """Per-gene adjacent CSV is created regardless of the main output format."""
    runner = CliRunner()
    output = tmp_path / "stats.json"
    per_gene_output = tmp_path / "stats.per-gene.csv"

    result = runner.invoke(
        cli,
        [
            "pretree",
            "stats",
            "--seq-dir",
            str(TEST_DATA / "test"),
            "--per-gene",
            "--output",
            str(output),
        ],
    )

    assert result.exit_code == 0
    assert per_gene_output.exists()
    assert f"Per-gene table saved to {per_gene_output}" in result.output


def test_directory_per_gene_format_can_write_tsv(tmp_path: Path) -> None:
    runner = CliRunner()
    output = tmp_path / "stats.txt"
    per_gene_output = tmp_path / "stats.per-gene.tsv"

    result = runner.invoke(
        cli,
        [
            "pretree",
            "stats",
            "--seq-dir",
            str(TEST_DATA / "test"),
            "--per-gene",
            "--per-gene-format",
            "tsv",
            "--output",
            str(output),
        ],
    )

    assert result.exit_code == 0
    assert per_gene_output.exists()
    assert "\t" in per_gene_output.read_text().splitlines()[0]


def test_unaligned_per_gene_output_includes_sequence_length_summary(tmp_path: Path) -> None:
    from phyloai.pretree.stats import write_per_gene_output

    per_gene_output = tmp_path / "stats.per-gene.csv"
    payload = {
        "data": {
            "per_gene": [
                {
                    "gene": "g1",
                    "n_taxa": 4,
                    "n_taxa_ratio": 1.0,
                    "length_type": "seq_length",
                    "alignment_length": "",
                    "seq_length_min": 100,
                    "seq_length_max": 120,
                    "seq_length_mean": 110.0,
                    "seq_length_median": 110.0,
                    "seq_length_stdev": 8.165,
                    "gap_ratio": 0.0,
                    "ambiguous_ratio": 0.0,
                    "gap_ambiguous_ratio": 0.0,
                    "missing_taxa": 0,
                    "missing_taxa_ratio": 0.0,
                }
            ]
        }
    }

    write_per_gene_output(payload, per_gene_output)

    header = per_gene_output.read_text().splitlines()[0].split(",")
    assert "alignment_length" not in header
    for column in [
        "length_type",
        "seq_length_min",
        "seq_length_max",
        "seq_length_mean",
        "seq_length_median",
        "seq_length_stdev",
    ]:
        assert column in header


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
    assert "Write full results to" in result.output


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


def test_directory_output_file_keeps_terminal_to_summary_only(tmp_path: Path) -> None:
    runner = CliRunner()
    output = tmp_path / "stats.txt"
    per_gene_output = tmp_path / "stats.per-gene.csv"

    result = runner.invoke(
        cli,
        [
            "pretree",
            "stats",
            "--seq-dir",
            str(TEST_DATA / "test"),
            "--per-gene",
            "--output",
            str(output),
            "--output-format",
            "text",
        ],
    )

    assert result.exit_code == 0
    assert "pretree stats summary" in result.output
    assert "Per-gene statistics" not in result.output
    assert str(output) in result.output
    assert str(per_gene_output) in result.output
    assert "Summary saved to" in result.output
    assert "Per-gene table saved to" in result.output


def test_directory_per_gene_writes_adjacent_table_not_inline(tmp_path: Path) -> None:
    """Per-gene data goes to the adjacent CSV, not embedded in the main output file."""
    runner = CliRunner()
    output = tmp_path / "stats.json"
    per_gene_output = tmp_path / "stats.per-gene.csv"

    result = runner.invoke(
        cli,
        [
            "pretree",
            "stats",
            "--seq-dir",
            str(TEST_DATA / "test"),
            "--per-gene",
            "--output",
            str(output),
        ],
    )

    assert result.exit_code == 0
    saved = json.loads(output.read_text())
    # per_gene key is present in payload (included because --per-gene was set)
    assert "per_gene" in saved["data"]
    per_gene_content = per_gene_output.read_text()
    assert per_gene_content.splitlines()[0] == (
        "gene,n_taxa,n_taxa_ratio,length_type,alignment_length,"
        "gap_ratio,ambiguous_ratio,gap_ambiguous_ratio,missing_taxa,missing_taxa_ratio"
    )
    assert "EOG090X0971" in per_gene_content


def test_directory_per_gene_csv_output_writes_adjacent_csv(tmp_path: Path) -> None:
    """With --per-gene, a .per-gene.csv file is written adjacent to the main output."""
    runner = CliRunner()
    output = tmp_path / "stats.json"
    per_gene_output = tmp_path / "stats.per-gene.csv"

    result = runner.invoke(
        cli,
        [
            "pretree",
            "stats",
            "--seq-dir",
            str(TEST_DATA / "test"),
            "--per-gene",
            "--output",
            str(output),
        ],
    )

    assert result.exit_code == 0
    assert output.exists()
    assert per_gene_output.exists()
    saved = json.loads(output.read_text())
    assert saved["data"]["summary"]["n_genes"] > 0
    assert per_gene_output.read_text().splitlines()[0] == (
        "gene,n_taxa,n_taxa_ratio,length_type,alignment_length,"
        "gap_ratio,ambiguous_ratio,gap_ambiguous_ratio,missing_taxa,missing_taxa_ratio"
    )


def test_directory_csv_output_reports_per_gene_destination(tmp_path: Path) -> None:
    runner = CliRunner()
    output = tmp_path / "stats.csv"
    per_gene_output = tmp_path / "stats.per-gene.csv"

    result = runner.invoke(
        cli,
        [
            "pretree",
            "stats",
            "--seq-dir",
            str(TEST_DATA / "test"),
            "--per-gene",
            "--output",
            str(output),
        ],
    )

    assert result.exit_code == 0
    assert f"Per-gene table saved to {per_gene_output}" in result.output


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


def test_single_file_text_output_mentions_msa_length() -> None:
    runner = CliRunner()

    result = runner.invoke(
        cli,
        [
            "pretree",
            "stats",
            "--seq",
            str(TEST_DATA / "test" / "EOG090X0971.faa"),
            "--output-format",
            "text",
        ],
    )

    assert result.exit_code == 0
    assert "MSA length" in result.output
    assert "1042" in result.output
    assert "distinct_patterns" in result.output


def test_single_file_text_output_shows_gap_ambiguous_per_taxon() -> None:
    runner = CliRunner()

    result = runner.invoke(
        cli,
        ["pretree", "stats", "--seq", str(TEST_DATA / "test" / "EOG090X0971.faa")],
    )

    assert result.exit_code == 0
    assert "gap_ambiguous_ratio" in result.output


def test_single_file_json_output_file_includes_per_taxon_data(tmp_path: Path) -> None:
    """Default JSON output for single file includes per_taxon array."""
    runner = CliRunner()
    output = tmp_path / "out.json"

    result = runner.invoke(
        cli,
        [
            "pretree",
            "stats",
            "--seq",
            str(TEST_DATA / "2-loci_filter" / "fna" / "EOG090X0971.fna"),
            "--output",
            str(output),
        ],
    )

    assert result.exit_code == 0
    saved = json.loads(output.read_text())
    assert "per_taxon" in saved["data"]
    assert len(saved["data"]["per_taxon"]) > 0
