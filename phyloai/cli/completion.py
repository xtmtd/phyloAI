"""Shell completion helpers for the PhyloAI CLI."""

from __future__ import annotations

import click
from click.shell_completion import get_completion_class


def _completion_script(shell: str) -> str:
    from phyloai.cli.main import cli

    completion_class = get_completion_class(shell)
    if completion_class is None:
        raise click.ClickException(f"Unsupported shell '{shell}'.")

    return completion_class(
        cli=cli,
        ctx_args={},
        prog_name="phyloai",
        complete_var="_PHYLOAI_COMPLETE",
    ).source()


@click.group(help="Generate shell completion scripts for static installation.")
def completion() -> None:
    """Generate shell completion scripts for PhyloAI."""


@completion.command(help="Print a Bash completion script for static sourcing.")
def bash() -> None:
    click.echo(_completion_script("bash"), nl=False)


@completion.command(help="Print a Zsh completion script for static sourcing.")
def zsh() -> None:
    click.echo(_completion_script("zsh"), nl=False)


@completion.command(help="Print a Fish completion script for static sourcing.")
def fish() -> None:
    click.echo(_completion_script("fish"), nl=False)
