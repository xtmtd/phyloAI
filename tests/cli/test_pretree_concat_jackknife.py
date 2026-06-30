from __future__ import annotations

import json
from pathlib import Path

from click.testing import CliRunner

from phyloai.cli.main import cli


def test_cli_pretree_concat_jackknife_help() -> None:
    result = CliRunner().invoke(cli, ["pretree", "concat", "jackknife", "--help"])
    assert result.exit_code == 0
    for flag in ["--matrix", "--partitions", "--replicates", "--target-length", "--table-format", "--seed"]:
        assert flag in result.output


def test_cli_pretree_concat_jackknife_writes_outputs(tmp_path: Path) -> None:
    matrix = tmp_path / "matrix.fa"
    matrix.write_text(">A\nAAAACCCC\n>B\nAAAACCCC\n")
    parts = tmp_path / "matrix.partitions"
    parts.write_text("LG, gene1 = 1-4\nLG, gene2 = 5-8\n")
    out = tmp_path / "jackknife"

    result = CliRunner().invoke(
        cli,
        [
            "pretree", "concat", "jackknife",
            "--matrix", str(matrix),
            "--partitions", str(parts),
            "--replicates", "1",
            "--target-length", "4",
            "--table-format", "tsv",
            "--output-dir", str(out),
            "--quiet",
        ],
    )

    assert result.exit_code == 0, result.output
    assert (out / "rep001" / "rep001.fa").exists()
    assert (out / "jackknife_summary.tsv").exists()
    payload = json.loads((out / "result.json").read_text())
    assert payload["params"]["table_format"] == "tsv"
