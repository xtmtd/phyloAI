from __future__ import annotations

from pathlib import Path

import pytest


def test_dayhoff6_recoding_maps_all_standard_amino_acids() -> None:
    from phyloai.pretree.concat import AA_RECODING_TABLES

    table = AA_RECODING_TABLES["Dayhoff-6"]
    assert table["A"] == "0"
    assert table["G"] == "0"
    assert table["P"] == "0"
    assert table["S"] == "0"
    assert table["T"] == "0"
    assert table["D"] == "1"
    assert table["E"] == "1"
    assert table["N"] == "1"
    assert table["Q"] == "1"
    assert table["H"] == "2"
    assert table["K"] == "2"
    assert table["R"] == "2"
    assert table["I"] == "3"
    assert table["L"] == "3"
    assert table["M"] == "3"
    assert table["V"] == "3"
    assert table["F"] == "4"
    assert table["W"] == "4"
    assert table["Y"] == "4"
    assert table["C"] == "5"


def test_ry_nucleotide_recoding_includes_ambiguous() -> None:
    from phyloai.pretree.concat import NT_RECODING_TABLES

    table = NT_RECODING_TABLES["RY-nucleotide"]
    assert table["A"] == "R"
    assert table["G"] == "R"
    assert table["C"] == "Y"
    assert table["T"] == "Y"
    assert table["U"] == "Y"
    assert table["N"] == "?"
    assert table["X"] == "?"
    assert table["-"] == "-"
    assert table["?"] == "?"
    assert table["."] == "."


def test_apply_recoding_dayhoff6() -> None:
    from phyloai.pretree.concat import _apply_recoding

    matrix = {"tax1": "ACDEFG", "tax2": "HIKLMN"}
    recoded, warnings = _apply_recoding(matrix, "Dayhoff-6")
    assert recoded["tax1"] == "051140"
    assert recoded["tax2"] == "232331"
    assert warnings == []


def test_apply_recoding_preserves_gaps_and_special_chars() -> None:
    from phyloai.pretree.concat import _apply_recoding

    matrix = {"tax1": "A-H.?*X"}
    recoded, warnings = _apply_recoding(matrix, "Dayhoff-6")
    assert recoded["tax1"] == "0-2.?*X"
    assert warnings == []


def test_apply_recoding_unknown_scheme_raises() -> None:
    from phyloai.pretree.concat import _apply_recoding

    with pytest.raises(ValueError, match="Unknown recoding scheme"):
        _apply_recoding({"tax1": "ACG"}, "FakeScheme")


def test_apply_recoding_ry_nucleotide() -> None:
    from phyloai.pretree.concat import _apply_recoding

    matrix = {"tax1": "ACGTN--?."}
    recoded, warnings = _apply_recoding(matrix, "RY-nucleotide")
    assert recoded["tax1"] == "RYRY?--?."
