from __future__ import annotations

import json
from pathlib import Path

from click.testing import CliRunner

from phyloai.cli.main import cli


def test_pretree_convert_help_is_registered_before_stats() -> None:
    result = CliRunner().invoke(cli, ["pretree", "--help"])

    assert result.exit_code == 0
    assert result.output.index("convert") < result.output.index("stats")


def test_cli_pretree_convert_single_file_json(tmp_path: Path) -> None:
    src = tmp_path / "gene.fa"
    src.write_text(">tax one\nacgu?\n")
    out_dir = tmp_path / "converted"

    result = CliRunner().invoke(
        cli,
        [
            "pretree",
            "convert",
            "--input",
            str(src),
            "--output-dir",
            str(out_dir),
            "--to",
            "fasta",
            "--seq-type",
            "NT",
            "--quiet",
        ],
    )

    assert result.exit_code == 0
    payload = json.loads(result.output)
    assert payload["status"] == "success"
    assert payload["data"]["summary"]["n_converted"] == 1
    assert (out_dir / "gene.fa").exists()


def test_cli_pretree_convert_all_failed_exits_one(tmp_path: Path) -> None:
    src_dir = tmp_path / "raw"
    src_dir.mkdir()
    (src_dir / "empty.fa").write_text("")

    result = CliRunner().invoke(cli, ["pretree", "convert", "--input", str(src_dir), "--output-dir", str(tmp_path / "out")])

    assert result.exit_code == 1
    assert "All input entries failed" in result.output
