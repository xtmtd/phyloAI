from __future__ import annotations

from click.testing import CliRunner

from phyloai.cli.main import cli


def test_root_help_does_not_expose_config_option() -> None:
    result = CliRunner().invoke(cli, ["-h"])

    assert result.exit_code == 0, result.output
    assert "--config" not in result.output
    assert "YAML config" not in result.output
