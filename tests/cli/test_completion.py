from click.testing import CliRunner

from phyloai.cli.main import cli


def test_completion_group_is_registered() -> None:
    runner = CliRunner()

    result = runner.invoke(cli, ["completion", "--help"])

    assert result.exit_code == 0
    assert "Generate shell completion scripts" in result.output
    assert "bash" in result.output
    assert "zsh" in result.output
    assert "fish" in result.output


def test_completion_bash_outputs_script() -> None:
    runner = CliRunner()

    result = runner.invoke(cli, ["completion", "bash"])

    assert result.exit_code == 0
    assert "complete -F" in result.output
    assert "phyloai" in result.output


def test_completion_zsh_outputs_script() -> None:
    runner = CliRunner()

    result = runner.invoke(cli, ["completion", "zsh"])

    assert result.exit_code == 0
    assert "#compdef phyloai" in result.output
    assert "compdef" in result.output


def test_completion_fish_outputs_script() -> None:
    runner = CliRunner()

    result = runner.invoke(cli, ["completion", "fish"])

    assert result.exit_code == 0
    assert "complete --command phyloai" in result.output


def test_completion_help_explains_static_usage() -> None:
    runner = CliRunner()

    result = runner.invoke(cli, ["completion", "bash", "--help"])

    assert result.exit_code == 0
    assert "Print a Bash completion script for static sourcing." in result.output
