import csv
import json
import re
import shlex
from pathlib import Path

import pytest
from Bio import AlignIO

from phyloai.posttree.syserror_rate import (
    RateRow,
    canonical_rates,
    parse_iqtree_rate,
    parse_fractions,
    parse_pb_rate,
    rate_source,
    read_matrix,
    run_rate,
    select_sites,
    subset_label,
    write_subset_fasta,
)


def test_iqtree_rate_keeps_one_based_sites_and_sorts_ties(tmp_path: Path) -> None:
    source = tmp_path / "matrix.rate"
    source.write_text("# comment\nSite\tRate\tCat\n2\t0.5\t1\n1\t0.5\t1\n3\t1.0\t2\n")
    assert canonical_rates(parse_iqtree_rate(source)) == [
        RateRow(site=1, rate=0.5),
        RateRow(site=2, rate=0.5),
        RateRow(site=3, rate=1.0),
    ]


def test_pb_rate_converts_zero_based_sites_to_one_based(tmp_path: Path) -> None:
    source = tmp_path / "chain.meansiterates"
    source.write_text("0 1.2\n1 0.2\n2 0.6\n")
    assert canonical_rates(parse_pb_rate(source)) == [
        RateRow(site=2, rate=0.2),
        RateRow(site=3, rate=0.6),
        RateRow(site=1, rate=1.2),
    ]


@pytest.mark.parametrize("iqtree,pb", [(None, None), (Path("a.rate"), Path("b.meansiterates"))])
def test_exactly_one_source_is_required(iqtree: Path | None, pb: Path | None) -> None:
    with pytest.raises(ValueError, match="exactly one"):
        rate_source(iqtree, pb)


@pytest.mark.parametrize(
    "text,match",
    [
        ("Site\tRate\n1\t-1\n", "negative"),
        ("Site\tRate\n1\tnan\n", "finite"),
        ("Site\tRate\n1\t0.1\n1\t0.2\n", "duplicate"),
        ("Site\tRate\n1\t0.1\n3\t0.2\n", "consecutive"),
    ],
)
def test_iqtree_rate_rejects_invalid_rows(tmp_path: Path, text: str, match: str) -> None:
    source = tmp_path / "bad.rate"
    source.write_text(text)
    with pytest.raises(ValueError, match=match):
        parse_iqtree_rate(source)


@pytest.mark.parametrize(
    "text,line_number",
    [
        ("Site\tRate\n1\t0.1\n1\t0.2\n", 3),
        ("Site\tRate\n1\t0.1\n3\t0.2\n", 3),
    ],
)
def test_iqtree_rate_duplicate_and_nonconsecutive_sites_name_path_and_line(
    tmp_path: Path, text: str, line_number: int
) -> None:
    source = tmp_path / "bad.rate"
    source.write_text(text)

    with pytest.raises(ValueError, match=rf"{re.escape(str(source))}.*line {line_number}"):
        parse_iqtree_rate(source)


def test_slow_and_fast_selection_round_up_and_restore_site_order() -> None:
    ranked = [RateRow(4, 0.1), RateRow(2, 0.2), RateRow(1, 0.3), RateRow(3, 0.4)]
    assert select_sites(ranked, "slow", 0.26) == [2, 4]
    assert select_sites(ranked, "fast", 0.5) == [1, 3]


def test_subset_fasta_is_aligned_and_wrapped(tmp_path: Path) -> None:
    matrix = tmp_path / "matrix.fa"
    matrix.write_text(">A\n" + "A" * 61 + "\n>B\n" + "C" * 61 + "\n")
    output = tmp_path / "subset.fa"
    write_subset_fasta(read_matrix(matrix), [1, 61], output)
    assert output.read_text() == ">A\nAA\n>B\nCC\n"
    assert AlignIO.read(output, "fasta").get_alignment_length() == 2


def test_subset_fasta_wraps_selected_sequences_at_60_columns(tmp_path: Path) -> None:
    matrix = tmp_path / "matrix.fa"
    matrix.write_text(">A\n" + "A" * 61 + "\n>B\n" + "C" * 61 + "\n")
    output = tmp_path / "subset.fa"

    write_subset_fasta(read_matrix(matrix), list(range(1, 62)), output)

    assert output.read_text() == ">A\n" + "A" * 60 + "\nA\n>B\n" + "C" * 60 + "\nC\n"


@pytest.mark.parametrize("value", [None, "", "0", "1.1", "0.25,0.25", "0.25,"])
def test_fraction_parser_rejects_invalid_lists(value: str | None) -> None:
    with pytest.raises(ValueError):
        parse_fractions(value)


def test_fraction_parser_and_labels_reject_label_collisions() -> None:
    assert parse_fractions("0.25,0.125,1") == [0.25, 0.125, 1.0]
    assert [subset_label("slow", value) for value in parse_fractions("0.25,0.125,1")] == [
        "slow25", "slow12.5", "slow100"
    ]
    with pytest.raises(ValueError):
        parse_fractions("0.1,0.10000000000001")


def test_read_matrix_accepts_relaxed_phylip(tmp_path: Path) -> None:
    matrix = tmp_path / "matrix.phy"
    matrix.write_text("2 3\nA  ACG\nB  TTT\n")
    assert read_matrix(matrix).get_alignment_length() == 3


@pytest.mark.parametrize(
    "name,text",
    [
        ("matrix.paml", "2 3\nA\nACG\nB\nTTT\n"),
        (
            "matrix.nex",
            "#NEXUS\nbegin data;\ndimensions ntax=2 nchar=3;\nformat datatype=dna gap=-;\nmatrix\nA ACG\nB TTT\n;\nend;\n",
        ),
    ],
)
def test_read_matrix_accepts_documented_paml_phylip_and_nexus(
    tmp_path: Path, name: str, text: str
) -> None:
    matrix = tmp_path / name
    matrix.write_text(text)

    assert read_matrix(matrix).get_alignment_length() == 3


@pytest.mark.parametrize(
    "text,match",
    [
        ("0 0.1 extra\n", "malformed"),
        ("0 0.1\n0 0.2\n", "duplicate"),
        ("0 0.1\n2 0.2\n", "consecutive"),
        ("0 nan\n", "finite"),
        ("0 inf\n", "finite"),
        ("0 -0.1\n", "negative"),
        ("0 -0.0\n", "negative"),
    ],
)
def test_pb_rate_rejects_invalid_rows(tmp_path: Path, text: str, match: str) -> None:
    source = tmp_path / "chain.meansiterates"
    source.write_text(text)

    with pytest.raises(ValueError, match=match):
        parse_pb_rate(source)


@pytest.mark.parametrize(
    "text",
    [">A\n\n", ">A\nAC\n>B\nA\n", ">A\nAC\n>A\nAC\n"],
)
def test_read_matrix_rejects_invalid_alignments(tmp_path: Path, text: str) -> None:
    matrix = tmp_path / "matrix.fa"
    matrix.write_text(text)
    with pytest.raises(ValueError):
        read_matrix(matrix)


def test_rate_only_writes_sorted_two_column_csv_and_result(tmp_path: Path) -> None:
    rate = tmp_path / "matrix.rate"
    rate.write_text("Site\tRate\n1\t0.4\n2\t0.1\n3\t0.2\n")

    result = run_rate(rate, None, output_dir=tmp_path / "out", quiet=True)

    assert list(csv.DictReader((tmp_path / "out" / "rates.csv").open())) == [
        {"site": "2", "rate": "0.1"},
        {"site": "3", "rate": "0.2"},
        {"site": "1", "rate": "0.4"},
    ]
    assert result["key_results"]["n_sites"] == 3
    assert json.loads((tmp_path / "out" / "result.json").read_text())["error"] is None


def test_multi_fraction_fast_outputs_positions_and_matrices(tmp_path: Path) -> None:
    rate = tmp_path / "rates"
    rate.write_text("0 0.1\n1 0.2\n2 0.3\n3 0.4\n")
    matrix = tmp_path / "matrix.fa"
    matrix.write_text(">A\nABCD\n>B\nWXYZ\n")

    result = run_rate(
        None, rate, matrix, subset="fast", fraction="0.25,0.5",
        output_dir=tmp_path / "out", quiet=True,
    )

    assert (tmp_path / "out" / "fast25" / "positions.txt").read_text() == "4\n"
    assert (tmp_path / "out" / "fast50" / "matrix.fa").read_text() == ">A\nCD\n>B\nYZ\n"
    assert [item["selected_sites"] for item in result["key_results"]["subsets"]] == [1, 2]


def test_matrix_length_mismatch_fails_before_output_creation(tmp_path: Path) -> None:
    rate = tmp_path / "rates"
    rate.write_text("0 0.1\n1 0.2\n")
    matrix = tmp_path / "matrix.fa"
    matrix.write_text(">A\nABC\n>B\nABC\n")

    with pytest.raises(ValueError, match="length"):
        run_rate(None, rate, matrix, fraction="0.5", output_dir=tmp_path / "out", quiet=True)
    assert not (tmp_path / "out").exists()


def test_dry_run_validates_extraction_but_writes_nothing(tmp_path: Path) -> None:
    rate = tmp_path / "rates"
    rate.write_text("0 0.1\n1 0.2\n")
    matrix = tmp_path / "matrix.fa"
    matrix.write_text(">A\nAC\n>B\nTG\n")

    result = run_rate(
        None, rate, matrix, fraction="0.5", output_dir=tmp_path / "out", dry_run=True, quiet=True,
    )

    assert result["key_results"]["subsets"] == [{
        "subset": "slow",
        "requested_fraction": 0.5,
        "selected_sites": 1,
        "actual_fraction": 0.5,
        "output_dir": str((tmp_path / "out" / "slow50").resolve()),
    }]
    assert not (tmp_path / "out").exists()


def test_success_result_command_is_a_complete_shell_quoted_call(tmp_path: Path) -> None:
    rate = tmp_path / "rate file"
    rate.write_text("Site\tRate\n1\t0.1\n")
    output = tmp_path / "output dir"

    result = run_rate(rate, None, output_dir=output, overwrite=True, quiet=True)

    assert result["command"] == shlex.join([
        "phyloai", "posttree", "syserror", "rate",
        "--iqtree-rate", str(rate.resolve()),
        "--output-dir", str(output.resolve()),
        "--overwrite", "--quiet",
    ])


def test_extraction_options_without_matrix_are_rejected(tmp_path: Path) -> None:
    rate = tmp_path / "rates"
    rate.write_text("0 0.1\n")

    with pytest.raises(ValueError, match="--matrix"):
        run_rate(None, rate, subset="fast", fraction="0.5", output_dir=tmp_path / "out", quiet=True)


def test_matrix_without_fraction_is_rejected(tmp_path: Path) -> None:
    rate = tmp_path / "rates"
    rate.write_text("0 0.1\n")
    matrix = tmp_path / "matrix.fa"
    matrix.write_text(">A\nA\n>B\nC\n")

    with pytest.raises(ValueError, match="--fraction"):
        run_rate(None, rate, matrix, output_dir=tmp_path / "out", quiet=True)


def test_nonempty_output_is_preserved_without_overwrite_and_replaced_with_it(tmp_path: Path) -> None:
    rate = tmp_path / "rates"
    rate.write_text("0 0.1\n")
    output = tmp_path / "out"
    output.mkdir()
    keep = output / "keep.txt"
    keep.write_text("existing data")

    with pytest.raises(ValueError, match="not empty"):
        run_rate(None, rate, output_dir=output, quiet=True)
    assert keep.read_text() == "existing data"

    run_rate(None, rate, output_dir=output, overwrite=True, quiet=True)
    assert not keep.exists()


def test_library_validation_never_deletes_existing_output_with_overwrite(tmp_path: Path) -> None:
    output = tmp_path / "out"
    output.mkdir()
    keep = output / "keep.txt"
    keep.write_text("existing data")

    with pytest.raises(ValueError, match="does not exist"):
        run_rate(tmp_path / "missing.rate", None, output_dir=output, overwrite=True, quiet=True)
    assert keep.read_text() == "existing data"


def test_output_dir_file_is_rejected_without_deleting_it(tmp_path: Path) -> None:
    rate = tmp_path / "rates"
    rate.write_text("0 0.1\n")
    output = tmp_path / "out"
    output.write_text("i am a file")

    with pytest.raises(ValueError, match="not a directory"):
        run_rate(None, rate, output_dir=output, quiet=True)
    assert output.read_text() == "i am a file"


def test_persisted_result_has_shared_fields_and_output_descriptions(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    rate = tmp_path / "rates"
    rate.write_text("0 0.1\n1 0.2\n")
    matrix = tmp_path / "matrix.fa"
    matrix.write_text(">A\nAC\n>B\nTG\n")

    output = tmp_path / "out"
    run_rate(None, rate, matrix, subset="slow", fraction="0.25", output_dir=output)

    result = json.loads((output / "result.json").read_text())
    assert {"status", "command", "wall_time", "error", "params", "key_results", "tool_versions", "data"} <= result.keys()
    assert result["params"]["quiet"] is False
    assert result["params"]["fraction"] == [0.25]
    files = result["data"]["output_files"]
    assert set(files) == {"rates", "slow25_positions", "slow25_matrix"}
    assert all(Path(item["path"]).is_absolute() and item["description"] for item in files.values())
    assert capsys.readouterr().out.splitlines() == [
        str(output / "rates.csv"),
        str(output / "slow25" / "matrix.fa"),
        str(output / "result.json"),
    ]
