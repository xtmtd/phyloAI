import csv
import json
import shlex

from pathlib import Path

import pytest
from matplotlib import colors

from phyloai.posttree.syserror_cca import (
    CcaRow,
    SiteFrequency,
    SiteLikelihood,
    build_cca_rows,
    parse_site_freq,
    parse_site_lnl,
    run_cca,
    summarize_bins,
    plot_cca,
    validate_matching_sites,
)


def _frequency_row(site: int) -> str:
    return f"{site} " + " ".join(["0.05"] * 20) + "\n"


def test_parse_site_freq_accepts_one_based_20_state_rows(tmp_path: Path) -> None:
    source = tmp_path / "chain.sitefreq"
    source.write_text(_frequency_row(1))

    row = parse_site_freq(source)[0]

    assert row.site == 1
    assert row.frequencies == (0.05,) * 20


@pytest.mark.parametrize("row,match", [
    ("0 " + " ".join(["0.05"] * 20), "consecutive from 1"),
    ("1 " + " ".join(["0.05"] * 19), "exactly 20"),
    ("1 " + " ".join(["0.10"] * 20), "sum to 1"),
    ("1 " + " ".join(["-0.05"] + ["0.05"] * 19), "non-negative"),
])
def test_parse_site_freq_rejects_invalid_rows(tmp_path: Path, row: str, match: str) -> None:
    source = tmp_path / "bad.sitefreq"
    source.write_text(row + "\n")

    with pytest.raises(ValueError, match=match):
        parse_site_freq(source)


def test_site_freq_parser_sorts_input_rows_by_site(tmp_path: Path) -> None:
    source = tmp_path / "chain.sitefreq"
    source.write_text(_frequency_row(2) + _frequency_row(1))

    assert [row.site for row in parse_site_freq(source)] == [1, 2]


def test_lnl_parser_requires_named_columns_but_ignores_extra_columns(tmp_path: Path) -> None:
    source = tmp_path / "site_lnl.csv"
    source.write_text(
        "site,lnL_Tree1,lnL_Tree2,ΔSLS,support\n"
        "2,-2.0,-1.0,-1.0,Tree1\n"
        "1,-4.0,-3.5,-0.5,Tree1\n"
    )

    assert [(row.site, row.lnl_tree1, row.lnl_tree2) for row in parse_site_lnl(source)] == [
        (1, -4.0, -3.5),
        (2, -2.0, -1.0),
    ]


@pytest.mark.parametrize("text,match", [
    ("site,lnL_Tree1\n1,-1\n", "header requires"),
    ("site,lnL_Tree1,lnL_Tree2\n1,-1,nan\n", "finite"),
    ("site,lnL_Tree1,lnL_Tree2\n1,-1,-2\n1,-2,-3\n", "duplicate"),
    ("site,lnL_Tree1,lnL_Tree2\n1,-1,-2,unexpected\n", "too many fields"),
    ("site,lnL_Tree1,lnL_Tree2\n1,-1,-2\n3,-2,-3\n", "consecutive from 1"),
])
def test_lnl_parser_rejects_invalid_rows(tmp_path: Path, text: str, match: str) -> None:
    source = tmp_path / "bad.csv"
    source.write_text(text)

    with pytest.raises(ValueError, match=match):
        parse_site_lnl(source)


def test_matching_inputs_reject_different_site_sets(tmp_path: Path) -> None:
    freq = tmp_path / "x.sitefreq"
    freq.write_text(_frequency_row(1))
    first = tmp_path / "first.csv"
    first.write_text("site,lnL_Tree1,lnL_Tree2\n1,-1,-2\n")
    second = tmp_path / "second.csv"
    second.write_text("site,lnL_Tree1,lnL_Tree2\n1,-1,-2\n2,-3,-4\n")

    with pytest.raises(ValueError, match="site sets must match"):
        validate_matching_sites(parse_site_freq(freq), parse_site_lnl(first), parse_site_lnl(second))


def test_cca_rows_and_bins_use_approved_formula_and_sign() -> None:
    rows = build_cca_rows(
        [SiteFrequency(1, (0.2,) * 5 + (0.0,) * 15)],
        [SiteLikelihood(1, -14.2296, -14.3580)],
        [SiteLikelihood(1, -13.8521, -13.9077)],
        "LG", "C20",
    )
    assert rows[0].keff == pytest.approx(5.0)
    assert [row.delta_lnl_tree2_tree1 for row in rows] == pytest.approx([-0.1284, -0.0556])

    sums = summarize_bins([
        CcaRow("LG", 1, 1.9, -1, 1, 2),
        CcaRow("LG", 2, 20.0, -1, 2, 3),
        CcaRow("C20", 1, 1.9, -1, 0, 1),
    ], ("LG", "C20"))
    assert len(sums["LG"]) == len(sums["C20"]) == 20
    assert sums["LG"][0] == 2
    assert sums["LG"][19] == 3


def test_plot_retains_bin_twenty_bars_with_training_x_limits(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    captured = {}
    import matplotlib.pyplot as plt
    original = plt.savefig

    def capture(*args, **kwargs):
        captured["axis"] = plt.gca()
        return original(*args, **kwargs)

    monkeypatch.setattr(plt, "savefig", capture)
    plot_cca({"LG": [0.0] * 19 + [1.0], "C20": [0.0] * 19 + [2.0]}, ("LG", "C20"), "", "X", "Y", 10, 6, 300, 16, tmp_path / "cca.pdf")

    assert captured["axis"].get_xlim() == (1.0, 21.0)
    assert all(len(container) == 20 for container in captured["axis"].containers)
    bars = [container[-1] for container in captured["axis"].containers]
    assert [(bar.get_x(), bar.get_x() + bar.get_width()) for bar in bars] == pytest.approx([(20.0, 20.5), (20.5, 21.0)])
    boundaries = [line for line in captured["axis"].lines if line.get_xdata()[0] == line.get_xdata()[1]]
    assert all(line.get_color() == "grey" and line.get_linewidth() == 0.1 for line in boundaries)

def _write_inputs(tmp_path: Path) -> tuple[Path, Path, Path]:
    frequencies = " ".join(["0.05"] * 20)
    freq = tmp_path / "chain.sitefreq"
    freq.write_text(f"1 {frequencies}\n2 {frequencies}\n")
    first = tmp_path / "first.csv"
    first.write_text("site,lnL_Tree1,lnL_Tree2\n2,-2,-1\n1,-4,-3.5\n")
    second = tmp_path / "second.csv"
    second.write_text("site,lnL_Tree1,lnL_Tree2\n2,-3,-4\n1,-5,-4\n")
    return freq, first, second


def test_plot_writes_pdf_with_documented_bar_colors(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    captured = {}
    import matplotlib.pyplot as plt
    original = plt.savefig

    def capture(*args, **kwargs):
        captured["figure"] = plt.gcf()
        return original(*args, **kwargs)

    monkeypatch.setattr(plt, "savefig", capture)
    output = tmp_path / "cca.pdf"
    plot_cca({"LG": [0.0] * 20, "C20": [1.0] + [0.0] * 19}, ("LG", "C20"), "", "X", "Y", 10, 6, 300, 16, output)

    axis = captured["figure"].axes[0]
    bar_colors = {colors.to_hex(bar.get_facecolor()) for container in axis.containers for bar in container}
    assert output.read_bytes().startswith(b"%PDF")
    bar_colours = {
        container.get_label(): colors.to_hex(container[0].get_facecolor())
        for container in axis.containers
    }
    assert bar_colours == {"LG": "#00bfc4", "C20": "#f8766d"}
    assert axis.get_xlim() == (1.0, 21.0)
    assert axis.get_legend().get_title().get_text() == ""
    assert any(line.get_visible() for line in axis.get_ygridlines())
    assert axis.xaxis.label.get_size() == axis.yaxis.label.get_size() == axis.title.get_size() == 11


def test_bundled_cca_fixture_matches_current_input_anchors() -> None:
    root = Path("runs/cca")
    rows = build_cca_rows(
        parse_site_freq(root / "chain1.sitefreq"),
        parse_site_lnl(root / "lnl_LG/site_lnl.csv"),
        parse_site_lnl(root / "lnl_C20/site_lnl.csv"),
        "LG", "C20",
    )

    assert rows[0].keff == pytest.approx(11.974845235298696)
    assert rows[0].delta_lnl_tree2_tree1 == pytest.approx(-0.0999)
    assert rows[1].delta_lnl_tree2_tree1 == pytest.approx(-0.0436)


def test_build_cca_command_omits_default_model_labels(tmp_path: Path) -> None:
    from phyloai.posttree.syserror_cca import build_cca_command

    command = build_cca_command(
        tmp_path / "freq", tmp_path / "first", tmp_path / "second", "model1", "model2",
        "", "Effective number of amino acids", "Log-likelihood difference", 10, 6, 300, 16,
        tmp_path / "out", False, False, False,
    )

    assert "--model1-name" not in command
    assert "--model2-name" not in command


def test_run_cca_writes_csv_pdf_and_result(tmp_path: Path) -> None:
    freq, first, second = _write_inputs(tmp_path)
    output = tmp_path / "out"

    payload = run_cca(freq, first, second, "LG", "C20", output_dir=output, quiet=True)

    assert {path.name for path in output.iterdir()} == {"cca.csv", "cca.pdf", "result.json"}
    rows = list(csv.DictReader((output / "cca.csv").open()))
    assert list(rows[0]) == ["model", "site", "keff", "lnl_tree1", "lnl_tree2", "delta_lnl_tree2_tree1"]
    assert [row["model"] for row in rows] == ["LG", "C20", "LG", "C20"]
    assert rows[0]["lnl_tree1"] == "-4.0000"
    assert rows[0]["lnl_tree2"] == "-3.5000"
    assert rows[0]["delta_lnl_tree2_tree1"] == "0.5000"
    assert float(rows[0]["keff"]) == pytest.approx(20.0)
    assert payload["key_results"]["n_sites"] == 2
    persisted = json.loads((output / "result.json").read_text())
    assert set(persisted["data"]["output_files"]) == {"cca_table", "cca_figure"}


def test_run_cca_dry_run_writes_nothing(tmp_path: Path) -> None:
    freq, first, second = _write_inputs(tmp_path)
    output = tmp_path / "out"

    result = run_cca(freq, first, second, output_dir=output, dry_run=True, quiet=True)

    assert result["data"]["output_files"] == {}
    assert not output.exists()


def test_run_cca_preserves_existing_directory_on_invalid_input(tmp_path: Path) -> None:
    output = tmp_path / "out"
    output.mkdir()
    keep = output / "keep.txt"
    keep.write_text("keep")

    with pytest.raises(ValueError, match="does not exist"):
        run_cca(tmp_path / "missing", tmp_path / "one", tmp_path / "two", output_dir=output, overwrite=True)
    assert keep.read_text() == "keep"
