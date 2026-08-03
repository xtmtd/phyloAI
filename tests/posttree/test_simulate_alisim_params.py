"""Tests for phyloai.posttree.simulate_alisim_params."""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from phyloai.posttree.simulate_alisim_params import (
    PARAM_COLUMNS,
    parse_iqtree_report,
    run_alisim_params,
)


def _write_report(path: Path, text: str) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)
    return path


def test_parse_gtr_f_i_g_and_tree_pair(tmp_path: Path) -> None:
    report = _write_report(
        tmp_path / "gene.iqtree",
        "Input data: DNA\nTo simulate an alignment of the same\n"
        'iqtree3 --alisim sim -m "GTR{1,2,3,4,5}+F{.1,.2,.3,.4}+I{.2}+G4{.7}" --length 100\n',
    )
    parsed = parse_iqtree_report(report)
    assert parsed["seqtype"] == "DNA"
    assert parsed["length"] == "100"
    assert parsed["subs_model"] == "GTR"
    assert parsed["subs_rate"] == "1/2/3/4/5"
    assert parsed["freq"] == ".1/.2/.3/.4"
    assert parsed["prop_inv"] == ".2"
    assert parsed["rate_heterogeneity"] == "G"
    assert parsed["rate_categories"] == "4"
    assert parsed["rate_param"] == ".7"


def test_parse_legal_empty_model_components(tmp_path: Path) -> None:
    report = _write_report(
        tmp_path / "aa.iqtree",
        "Input data: amino-acid\nTo simulate an alignment of the same\n"
        'iqtree3 --alisim sim -m "LG" --length 100\n',
    )
    parsed = parse_iqtree_report(report)
    assert parsed["seqtype"] == "AA"
    assert parsed["length"] == "100"
    assert parsed["subs_model"] == "LG"
    assert parsed["subs_rate"] == ""
    assert parsed["freq"] == ""
    assert parsed["prop_inv"] == ""
    assert parsed["rate_heterogeneity"] == ""
    assert parsed["rate_categories"] == ""
    assert parsed["rate_param"] == ""


def test_parse_freerate_pairs(tmp_path: Path) -> None:
    report = _write_report(
        tmp_path / "fr.iqtree",
        "Input data: amino-acid\nTo simulate an alignment of the same\n"
        'iqtree3 --alisim sim -m "LG+R2{0.863674,0.350977,0.136326,5.11179}" --length 1640\n',
    )
    parsed = parse_iqtree_report(report)
    assert parsed["seqtype"] == "AA"
    assert parsed["rate_heterogeneity"] == "R"
    assert parsed["rate_categories"] == "2"
    assert parsed["rate_param"] == "0.863674/0.350977/0.136326/5.11179"


def test_parse_aa_f_falls_back_to_pi_lines(tmp_path: Path) -> None:
    report = _write_report(
        tmp_path / "aaf.iqtree",
        "Input data: amino-acid\n"
        "pi(A) = 0.08\npi(R) = 0.06\npi(N) = 0.04\npi(D) = 0.05\npi(C) = 0.01\n"
        "pi(Q) = 0.03\npi(E) = 0.05\npi(G) = 0.07\npi(H) = 0.02\npi(I) = 0.06\n"
        "pi(L) = 0.09\npi(K) = 0.06\npi(M) = 0.02\npi(F) = 0.04\npi(P) = 0.04\n"
        "pi(S) = 0.07\npi(T) = 0.06\npi(W) = 0.01\npi(Y) = 0.03\npi(V) = 0.07\n"
        "To simulate an alignment of the same\n"
        'iqtree3 --alisim sim -m "LG+F" --length 100\n',
    )
    parsed = parse_iqtree_report(report)
    assert parsed["freq"].startswith("0.08/0.06/0.04/0.05")
    assert parsed["freq"].count("/") == 19


def test_parse_nt_no_f_leaves_freq_empty(tmp_path: Path) -> None:
    report = _write_report(
        tmp_path / "dna.iqtree",
        "Input data: nucleotide\nTo simulate an alignment of the same\n"
        'iqtree3 --alisim sim -m "GTR{1,2,3,4,5}+G4{0.5}" --length 300\n',
    )
    parsed = parse_iqtree_report(report)
    assert parsed["seqtype"] == "DNA"
    assert parsed["freq"] == ""
    assert parsed["subs_rate"] == "1/2/3/4/5"


def test_run_skips_unmatched_and_rejects_ambiguous_tree(tmp_path: Path) -> None:
    reports = tmp_path / "reports"
    trees = tmp_path / "trees"
    (reports / "nested").mkdir(parents=True)
    trees.mkdir()
    _write_report(
        reports / "nested" / "EOG090X002Z.iqtree",
        "Input data: amino-acid\nTo simulate an alignment of the same\n"
        'iqtree3 --alisim sim -m "LG+G4{0.6}" --length 100\n',
    )
    (trees / "EOG090X002Z.treefile").write_text("(A,B);\n")

    result = run_alisim_params(
        iqtree_dir=reports, tree_dir=trees, output_dir=tmp_path / "out",
    )
    assert result["status"] == "success"
    assert result["key_results"]["n_loci_parsed"] == 1
    assert result["key_results"]["n_loci_matched"] == 1
    assert result["key_results"]["n_loci_unmatched"] == 0
    assert result["data"]["unmatched"] == []
    assert (tmp_path / "out" / "params.tsv").exists()
    assert (tmp_path / "out" / "result.json").exists()


def test_run_reports_unmatched_locus(tmp_path: Path) -> None:
    reports = tmp_path / "reports"
    trees = tmp_path / "trees"
    reports.mkdir()
    trees.mkdir()
    _write_report(
        reports / "AAA.iqtree",
        "Input data: amino-acid\nTo simulate an alignment of the same\n"
        'iqtree3 --alisim sim -m "LG" --length 100\n',
    )
    result = run_alisim_params(
        iqtree_dir=reports, tree_dir=trees, output_dir=tmp_path / "out",
    )
    assert result["key_results"]["n_loci_unmatched"] == 1
    assert result["data"]["unmatched"][0]["id"] == "AAA"
    rows = (tmp_path / "out" / "params.tsv").read_text().splitlines()
    assert len(rows) == 1  # header only


def test_run_ambiguity_is_hard_error(tmp_path: Path) -> None:
    reports = tmp_path / "reports"
    trees = tmp_path / "trees"
    reports.mkdir()
    trees.mkdir()
    _write_report(
        reports / "AAA.iqtree",
        "Input data: amino-acid\nTo simulate an alignment of the same\n"
        'iqtree3 --alisim sim -m "LG" --length 100\n',
    )
    (trees / "AAA.treefile").write_text("(A,B);\n")
    (trees / "AAA.nwk").write_text("(A,B);\n")
    with pytest.raises(ValueError, match="[Aa]mbiguous"):
        run_alisim_params(
            iqtree_dir=reports, tree_dir=trees, output_dir=tmp_path / "out",
        )


def test_run_dry_run_writes_nothing(tmp_path: Path) -> None:
    reports = tmp_path / "reports"
    trees = tmp_path / "trees"
    reports.mkdir()
    trees.mkdir()
    _write_report(
        reports / "AAA.iqtree",
        "Input data: amino-acid\nTo simulate an alignment of the same\n"
        'iqtree3 --alisim sim -m "LG" --length 100\n',
    )
    (trees / "AAA.treefile").write_text("(A,B);\n")
    result = run_alisim_params(
        iqtree_dir=reports, tree_dir=trees,
        output_dir=tmp_path / "out", dry_run=True,
    )
    assert result["status"] == "success"
    assert not (tmp_path / "out").exists()


def test_run_writes_absolute_paths_and_tab_delimited(tmp_path: Path) -> None:
    reports = tmp_path / "reports"
    trees = tmp_path / "trees"
    reports.mkdir()
    trees.mkdir()
    _write_report(
        reports / "AAA.iqtree",
        "Input data: amino-acid\nTo simulate an alignment of the same\n"
        'iqtree3 --alisim sim -m "LG+G4{0.5}" --length 42\n',
    )
    (trees / "AAA.treefile").write_text("(A,B);\n")
    out = tmp_path / "out"
    run_alisim_params(iqtree_dir=reports, tree_dir=trees, output_dir=out)

    header = out / "params.tsv"
    first = header.read_text().splitlines()[0].split("\t")
    assert first == list(PARAM_COLUMNS)
    row = dict(zip(PARAM_COLUMNS, header.read_text().splitlines()[1].split("\t")))
    assert Path(row["tree_path"]).is_absolute()
    assert row["tree_path"].endswith("AAA.treefile")
    result = json.loads((out / "result.json").read_text())
    assert Path(result["data"]["output_files"]["params_tsv"]["path"]).is_absolute()
    assert result["params"]["iqtree_dir"] == str(reports.resolve())
