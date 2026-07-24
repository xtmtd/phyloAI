"""CLI tests for phyloai update."""

from __future__ import annotations

from unittest.mock import patch

from click.testing import CliRunner


def test_update_check_up_to_date_exits_0() -> None:
    from phyloai.cli.main import cli

    with patch("phyloai.core.update.check_update", return_value={
        "status": "up_to_date", "current": "0.3.0", "latest": "v0.3.0",
    }):
        result = CliRunner().invoke(cli, ["update", "--check"])
    assert result.exit_code == 0
    assert "up to date" in result.output


def test_update_check_available_exits_1() -> None:
    from phyloai.cli.main import cli

    with patch("phyloai.core.update.check_update", return_value={
        "status": "available", "current": "0.3.0", "latest": "v0.4.0",
    }):
        result = CliRunner().invoke(cli, ["update", "--check"])
    assert result.exit_code == 1
    assert "Update available" in result.output
    assert "0.3.0" in result.output
    assert "v0.4.0" in result.output


def test_update_check_error_exits_1() -> None:
    from phyloai.cli.main import cli

    with patch("phyloai.core.update.check_update", return_value={
        "status": "error", "message": "network error",
    }):
        result = CliRunner().invoke(cli, ["update", "--check"])
    assert result.exit_code == 1
    assert "Error" in result.output or "network error" in result.output


def test_update_yes_no_prompt() -> None:
    from phyloai.cli.main import cli

    with patch("phyloai.core.update.check_update", return_value={
        "status": "up_to_date", "current": "0.3.0", "latest": "v0.3.0",
    }):
        result = CliRunner().invoke(cli, ["update", "--yes"])
    assert result.exit_code == 0


def test_update_help_shows_options() -> None:
    from phyloai.cli.main import cli

    result = CliRunner().invoke(cli, ["update", "--help"])
    assert result.exit_code == 0
    assert "--check" in result.output
    assert "--yes" in result.output
