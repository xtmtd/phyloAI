from __future__ import annotations

import json
import shlex

import pytest
from click.testing import CliRunner

from phyloai.cli.main import cli
from phyloai.mcp.schema_gen import walk_click_tree


def test_rate_help_and_generated_mcp_leaf() -> None:
    result = CliRunner().invoke(cli, ["posttree", "syserror", "rate", "--help"])

    assert result.exit_code == 0
    for option in ("--iqtree-rate", "--pb-rate", "--matrix", "--subset", "--fraction", "--dry-run"):
        assert option in result.output
    assert "posttree_syserror_rate" in {item["tool_name"] for item in walk_click_tree(cli)}


def test_rate_cli_writes_error_result_but_dry_run_does_not(tmp_path) -> None:
    source = tmp_path / "rates"
    source.write_text("0 0.1\n")
    output = tmp_path / "out"

    failed = CliRunner().invoke(cli, [
        "posttree", "syserror", "rate", "--pb-rate", str(source), "--fraction", "0.5", "-o", str(output),
    ])
    assert failed.exit_code == 1
    assert (output / "result.json").exists()

    dry_output = tmp_path / "dry"
    dry = CliRunner().invoke(cli, [
        "posttree", "syserror", "rate", "--pb-rate", str(source), "--fraction", "0.5", "-o", str(dry_output), "--dry-run",
    ])
    assert dry.exit_code == 1
    assert not dry_output.exists()


def test_rate_validation_error_result_records_the_supplied_command(tmp_path) -> None:
    source = tmp_path / "rate file"
    source.write_text("0 0.1\n")
    output = tmp_path / "output dir"

    result = CliRunner().invoke(cli, [
        "posttree", "syserror", "rate", "--pb-rate", str(source),
        "--fraction", "0.5", "--output-dir", str(output), "--quiet",
    ])

    assert result.exit_code == 1
    payload = json.loads((output / "result.json").read_text())
    assert payload["command"] == shlex.join([
        "phyloai", "posttree", "syserror", "rate",
        "--pb-rate", str(source.resolve()), "--fraction", "0.5",
        "--output-dir", str(output.resolve()), "--quiet",
    ])


@pytest.mark.parametrize("subset", ["slow", "fast"])
def test_rate_validation_error_result_keeps_subset_without_matrix(tmp_path, subset) -> None:
    source = tmp_path / "rates"
    source.write_text("0 0.1\n")
    output = tmp_path / "out"

    result = CliRunner().invoke(cli, [
        "posttree", "syserror", "rate", "--pb-rate", str(source),
        "--subset", subset, "--output-dir", str(output),
    ])

    assert result.exit_code == 1
    payload = json.loads((output / "result.json").read_text())
    assert f"--subset {subset}" in payload["command"]


def test_rate_output_dir_file_fails_cleanly_without_result_json(tmp_path) -> None:
    source = tmp_path / "rates"
    source.write_text("0 0.1\n")
    output = tmp_path / "out"
    output.write_text("i am a file")

    result = CliRunner().invoke(cli, [
        "posttree", "syserror", "rate", "--pb-rate", str(source), "-o", str(output),
    ])

    assert result.exit_code == 1
    assert "not a directory" in result.output
    assert "Traceback" not in result.output
    assert output.read_text() == "i am a file"


def test_rate_invalid_paths_write_error_results(tmp_path) -> None:
    for option, value in (
        ("--pb-rate", tmp_path / "missing"),
        ("--pb-rate", tmp_path),
        ("--matrix", tmp_path / "missing"),
        ("--matrix", tmp_path),
    ):
        source = tmp_path / "rates"
        source.write_text("0 0.1\n")
        args = ["posttree", "syserror", "rate", "--pb-rate", str(source)]
        if option == "--pb-rate":
            args[args.index("--pb-rate") + 1] = str(value)
        else:
            args.extend([option, str(value), "--fraction", "0.5"])
        output = tmp_path / f"out-{option[2:]}-{value.name or 'dir'}"
        result = CliRunner().invoke(cli, [*args, "-o", str(output)])
        assert result.exit_code == 1
        assert (output / "result.json").exists()


def test_rate_overwrite_invalid_input_preserves_existing_files(tmp_path) -> None:
    output = tmp_path / "out"
    output.mkdir()
    (output / "keep.txt").write_text("existing data")

    result = CliRunner().invoke(cli, [
        "posttree", "syserror", "rate", "--pb-rate", str(tmp_path / "missing"),
        "--overwrite", "-o", str(output),
    ])

    assert result.exit_code == 1
    assert (output / "keep.txt").read_text() == "existing data"
    assert json.loads((output / "result.json").read_text())["status"] == "error"


def test_rate_cli_extracts_subset_and_registers_outputs(tmp_path) -> None:
    source = tmp_path / "rates"
    source.write_text("0 0.1\n1 0.2\n")
    matrix = tmp_path / "matrix.fa"
    matrix.write_text(">A\nAC\n>B\nGT\n")
    output = tmp_path / "out"

    result = CliRunner().invoke(cli, [
        "posttree", "syserror", "rate", "--pb-rate", str(source), "--matrix", str(matrix),
        "--fraction", "0.5", "-o", str(output),
    ])

    assert result.exit_code == 0
    assert (output / "slow50" / "matrix.fa").exists()
    payload = json.loads((output / "result.json").read_text())
    assert "rates" in payload["data"]["output_files"]
    assert "slow50_positions" in payload["data"]["output_files"]
