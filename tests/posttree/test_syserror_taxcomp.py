from __future__ import annotations

from pathlib import Path

import pytest

from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from scipy.stats import chi2

from phyloai.core.formats import AlignmentFormat, FormatConverter
from phyloai.posttree.simulate_adequacy import _compute_taxon_composition
from phyloai.posttree.syserror_taxcomp import (
    compute_taxcomp_statistics,
    holm_adjust,
    run_taxcomp,
    sparse_count_check,
)


def _msa(*rows: tuple[str, str]) -> MultipleSeqAlignment:
    return MultipleSeqAlignment([
        SeqRecord(Seq(sequence), id=taxon, description="")
        for taxon, sequence in rows
    ])


def test_compute_taxcomp_statistics_known_nt_table() -> None:
    overall, rows = compute_taxcomp_statistics(
        _msa(("A", "AAAA"), ("B", "CCCC"), ("C", "AACC")),
        "NT",
    )

    assert overall["n_taxa"] == 3
    assert overall["n_states"] == 2
    assert overall["df"] == 2
    assert overall["x2"] == pytest.approx(8.0)
    assert sum(row["x2_contribution"] for row in rows) == pytest.approx(overall["x2"])
    assert [row["df"] for row in rows] == [1, 1, 1]
    assert [row["taxon"] for row in rows] == ["A", "B", "C"]


def test_p_nominal_matches_scipy_survival() -> None:
    overall, rows = compute_taxcomp_statistics(
        _msa(("A", "AAAA"), ("B", "CCCC"), ("C", "AACC")),
        "NT",
    )

    assert overall["p_nominal"] == pytest.approx(chi2.sf(8.0, 2))
    assert rows[0]["p_nominal"] == pytest.approx(chi2.sf(4.0, 1))


def test_holm_adjusts_in_original_order_and_is_monotone() -> None:
    adjusted = holm_adjust([0.03, 0.001, 0.02, 0.20])

    assert adjusted == pytest.approx([0.06, 0.004, 0.06, 0.20])


def test_holm_equal_p_values_get_equal_adjusted_values() -> None:
    adjusted = holm_adjust([0.01, 0.01, 0.50])

    assert adjusted == pytest.approx([0.03, 0.03, 0.50])


def test_sparse_count_check_has_strict_boundaries() -> None:
    assert sparse_count_check([1.0, 5.0, 5.0, 5.0, 5.0]) == {
        "sparse_count_check": "not_triggered",
        "expected_cells_total": 5,
        "expected_cells_below_1": 0,
        "expected_cells_below_5": 1,
        "expected_cells_below_5_fraction": 0.2,
    }


def test_one_expected_cell_below_one_triggers() -> None:
    result = sparse_count_check([0.999, 5.0, 5.0, 5.0, 5.0])

    assert result["sparse_count_check"] == "triggered"
    assert result["expected_cells_below_1"] == 1


def test_more_than_twenty_percent_below_five_triggers() -> None:
    cells = [4.0] * 6 + [5.0] * 19

    result = sparse_count_check(cells)

    assert result["sparse_count_check"] == "triggered"
    assert result["expected_cells_total"] == 25
    assert result["expected_cells_below_5"] == 6


def test_globally_absent_states_reduce_n_states() -> None:
    overall, _ = compute_taxcomp_statistics(
        _msa(("A", "AA"), ("B", "AC"), ("C", "CA")),
        "NT",
    )

    assert overall["n_states"] == 2
    assert overall["df"] == (3 - 1) * (2 - 1)


def test_gaps_and_ambiguity_codes_are_ignored() -> None:
    with_gaps = _msa(("A", "AA--N"), ("B", "CC--N"), ("C", "AC--N"))
    clean = _msa(("A", "AA"), ("B", "CC"), ("C", "AC"))

    overall_g, _ = compute_taxcomp_statistics(with_gaps, "NT")
    overall_c, _ = compute_taxcomp_statistics(clean, "NT")

    assert overall_g["x2"] == pytest.approx(overall_c["x2"])
    assert overall_g["n_states"] == overall_c["n_states"]


def test_ppa_comp_matches_adequacy_helper() -> None:
    alignment = _msa(("A", "AAAC"), ("B", "CCCA"), ("C", "AACC"))

    overall, rows = compute_taxcomp_statistics(alignment, "NT")
    comp = _compute_taxon_composition(alignment, "NT")

    assert overall["comp_max"] == pytest.approx(comp["comp_max"])
    assert overall["comp_mean"] == pytest.approx(comp["comp_mean"])
    by_taxon = {row["taxon"]: row["squared_composition_distance"] for row in rows}
    assert by_taxon == pytest.approx(comp["taxon_dist_j"])


def test_ppa_comp_uses_equal_taxon_mean_not_pooled() -> None:
    alignment = _msa(("A", "AAAAC"), ("B", "CC---"), ("C", "AC---"))

    _, rows = compute_taxcomp_statistics(alignment, "NT")

    freqs = {"A": [4 / 5, 1 / 5], "B": [0.0, 1.0], "C": [0.5, 0.5]}
    mean = [sum(freq[i] for freq in freqs.values()) / 3 for i in range(2)]
    dist_a = sum((freqs["A"][i] - mean[i]) ** 2 for i in range(2))
    by_taxon = {row["taxon"]: row["squared_composition_distance"] for row in rows}

    assert by_taxon["A"] == pytest.approx(dist_a)


def test_all_missing_taxon_raises() -> None:
    with pytest.raises(ValueError):
        compute_taxcomp_statistics(_msa(("A", "AAAA"), ("B", "----")), "NT")


def test_duplicate_taxon_raises() -> None:
    with pytest.raises(ValueError):
        compute_taxcomp_statistics(_msa(("A", "AAAA"), ("A", "CCCC")), "NT")


def test_unequal_lengths_raise() -> None:
    with pytest.raises(ValueError):
        compute_taxcomp_statistics(_msa(("A", "AAA"), ("B", "AAAA")), "NT")


def test_fewer_than_two_taxa_raises() -> None:
    with pytest.raises(ValueError):
        compute_taxcomp_statistics(_msa(("A", "AAAA")), "NT")


def test_fewer_than_two_global_states_raises() -> None:
    with pytest.raises(ValueError):
        compute_taxcomp_statistics(_msa(("A", "AAAA"), ("B", "AAAA")), "NT")


def test_auto_resolution_and_supported_formats_preserve_taxon_order(tmp_path: Path) -> None:
    """Design 11: auto AA/NT detection and FASTA/PHYLIP/PHYLIP-PAML/Nexus
    ingestion must resolve the type and preserve input taxon order."""
    import csv as _csv

    converter = FormatConverter()
    suffix = {
        AlignmentFormat.FASTA: ".fa",
        AlignmentFormat.PHYLIP: ".phy",
        AlignmentFormat.PHYLIP_PAML: ".paml.phy",
        AlignmentFormat.NEXUS: ".nex",
    }
    aa_aln = _msa(("TaxW", "WFWIYPWYVWG"), ("TaxF", "FFWWPFFYVWW"), ("TaxM", "MMPWMMFYVWM"))
    nt_aln = _msa(("N1", "ACGTACGTAC"), ("N2", "AACCGGTTAA"), ("N3", "GTACGTACGT"))

    for label, aln, expected_type, mol in (
        ("aa", aa_aln, "AA", "protein"),
        ("nt", nt_aln, "NT", "DNA"),
    ):
        for fmt in (AlignmentFormat.FASTA, AlignmentFormat.PHYLIP,
                    AlignmentFormat.PHYLIP_PAML, AlignmentFormat.NEXUS):
            path = tmp_path / f"{label}_{fmt.name}{suffix[fmt]}"
            converter.write_alignment(aln, path, target=fmt, molecule_type=mol)
            output = tmp_path / f"out_{label}_{fmt.name}"

            payload = run_taxcomp(path, seq_type="auto", output_dir=output, quiet=True)

            assert payload["params"]["detected_seq_type"] == expected_type
            assert payload["key_results"]["seq_type"] == expected_type
            rows = list(_csv.DictReader((output / "taxon_summary.csv").open()))
            assert [row["taxon"] for row in rows] == [rec.id for rec in aln]


def test_iqtree_per_taxon_composition_fixture_regression() -> None:
    """Freeze the IQ-TREE3 per-sequence composition chi2 test on a
    no-ambiguity AA alignment, verifying per-taxon df and nominal p-values
    match IQ-TREE's pooled expected-count convention."""
    alignment = _msa(
        ("TaxA", "ILVIIKDLMSRERTFAVVFMNALADQAIITLMIFIEVVTFEQEASHTHVTACYIVHAMSQREVSRQRASDLIDPNMRHCQVTEWCHWMRVCQSEEARKQKGIVIKWDNYTQDPVDVWKHRAVQFQTKEDPAQLWMEYCAQIADGKNGWVFFKNAAYFVAM"),
        ("TaxB", "ISGVRLSYSPTYHFYWHTWKMHKKRMVWPHSKIWRSAKIQFLKVEAGFQDYMIYGFCITYWHKGDLWWPVEDIEQHMLGWVIRAWAPKGTADIATQDAMIMMTVFCIDAIPDEGFLKSHTYRVPNIMYHLWTVIMQPTMPQYVVLTGYWGRILSIATRTD"),
        ("TaxC", "GSDDLPTFVLVRQLEPMLTKGHTYQWVLKQLLMAVLHTKLFRLHVVEINITYIGCMVNSACQFGGMNPNLECLIHHLDQMMLWMWQSTYGGVQVWDQPERTISNHFAQFYVNLSLFCRNAKLKDCQTVVRRVISGILIMAFTCIPHCDLAWTFWMHKLVT"),
        ("TaxD", "QYSSTHWWFSSHPRGIAGMMQVFGFKLHVFFSVAMPQHADHRQWKLVKVMTRKSAVDIDLMAVHLKGEFKEPEAQADSHWSAVVRWHIVQLVPIVTMNMICALNIFIEKHSLVVPQNWSMGEFKKMSDQVYWAAEHMNLRRFKMHICGSTRNTMSCHQTG"),
        ("TaxE", "FVGKSYHGRWRRKIFESRCFREPFEDGPGNPKYVYITEEHVWQGGGHHDKTKTKDSWEALKKVFTRKSRRSWRMDSFQQSLSHGVFLLPPYHTFWTTWVFLPTHHDSVHDCHEVAENFSDRIKWANDCQPKQDEDSEKTIQWYGKVIQAVTHLLEVFQAI"),
        ("TaxF", "HTPMDDGQWLMVELDSATWTREPRLIVKWESIVRTCKLDTGGREIPHFRHIVTFAHPQKIAVTCMEHLEGHFQDSGVKPGTRNRFPDHTMRQDNHRVILFSTHGIAHFVSFMAGGFRHLSTMDSYKHLNPVLPTDVMDRQKFGVHINVNYKKKEYVTMNT"),
    )
    # IQ-TREE3 3.1.2 per-sequence composition chi2 output (percentages, df=19)
    iqtree_pct = {"TaxA": 41.85, "TaxB": 32.37, "TaxC": 29.47,
                  "TaxD": 88.94, "TaxE": 52.83, "TaxF": 84.89}

    overall, rows = compute_taxcomp_statistics(alignment, "AA")

    assert overall["n_taxa"] == 6
    assert overall["n_states"] == 20
    for row in rows:
        assert row["df"] == 19
        assert row["p_nominal"] * 100 == pytest.approx(iqtree_pct[row["taxon"]], abs=0.01)


def test_run_taxcomp_writes_only_tables_and_result_json(tmp_path: Path) -> None:
    matrix = tmp_path / "matrix.fa"
    matrix.write_text(">A\nAAAA\n>B\nCCCC\n>C\nAACC\n")
    output = tmp_path / "out"

    payload = run_taxcomp(matrix, seq_type="NT", output_dir=output, quiet=True)

    assert {path.name for path in output.iterdir()} == {
        "overall_summary.csv", "taxon_summary.csv", "result.json",
    }
    assert payload["error_category"] is None
    assert set(payload["data"]["output_files"]) == {"overall_summary", "taxon_summary"}
    assert payload["key_results"]["seq_type"] == "NT"


def test_run_taxcomp_dry_run_writes_no_files(tmp_path: Path) -> None:
    matrix = tmp_path / "matrix.fa"
    matrix.write_text(">A\nAAAA\n>B\nCCCC\n>C\nAACC\n")

    payload = run_taxcomp(
        matrix, seq_type="NT", output_dir=tmp_path / "out",
        dry_run=True, quiet=True,
    )

    assert payload["data"]["output_files"] == {}
    assert not (tmp_path / "out").exists()


def test_run_taxcomp_tsv_uses_tsv_suffix_and_delimiter(tmp_path: Path) -> None:
    matrix = tmp_path / "matrix.fa"
    matrix.write_text(">A\nAAAA\n>B\nCCCC\n>C\nAACC\n")
    output = tmp_path / "out"

    payload = run_taxcomp(
        matrix, seq_type="NT", table_format="tsv", output_dir=output, quiet=True,
    )

    assert payload["data"]["output_files"]["overall_summary"]["path"].endswith(".tsv")
    header = (output / "overall_summary.tsv").read_text().splitlines()[0]
    assert "\t" in header
    assert "," not in header.split("\t")[0]


def test_run_taxcomp_nonempty_output_refused_without_overwrite(tmp_path: Path) -> None:
    matrix = tmp_path / "matrix.fa"
    matrix.write_text(">A\nAAAA\n>B\nCCCC\n>C\nAACC\n")
    output = tmp_path / "out"
    output.mkdir()
    (output / "keep.txt").write_text("existing")

    with pytest.raises(ValueError):
        run_taxcomp(matrix, seq_type="NT", output_dir=output, quiet=True)

    assert (output / "keep.txt").exists()


def test_run_taxcomp_overwrite_replaces_after_valid_preflight(tmp_path: Path) -> None:
    matrix = tmp_path / "matrix.fa"
    matrix.write_text(">A\nAAAA\n>B\nCCCC\n>C\nAACC\n")
    output = tmp_path / "out"
    output.mkdir()
    (output / "keep.txt").write_text("existing")

    run_taxcomp(matrix, seq_type="NT", output_dir=output, overwrite=True, quiet=True)

    assert not (output / "keep.txt").exists()
    assert (output / "result.json").exists()


def test_run_taxcomp_invalid_input_preserves_existing_output(tmp_path: Path) -> None:
    bad = tmp_path / "bad.fa"
    bad.write_text(">A\nAAAA\n>B\n----\n")
    output = tmp_path / "out"
    output.mkdir()
    (output / "keep.txt").write_text("existing")

    with pytest.raises(ValueError):
        run_taxcomp(bad, seq_type="NT", output_dir=output, quiet=True)

    assert (output / "keep.txt").exists()
    assert not (output / "result.json").exists()


def test_run_taxcomp_sparse_trigger_adds_warning(tmp_path: Path) -> None:
    matrix = tmp_path / "matrix.fa"
    matrix.write_text(">A\nA\n>B\nC\n>C\nA\n")
    output = tmp_path / "out"

    payload = run_taxcomp(matrix, seq_type="NT", output_dir=output, quiet=True)

    assert payload["key_results"]["sparse_count_check"] == "triggered"
    assert any("sparse" in warning for warning in payload["data"]["warnings"])
