"""PhyloAI CLI entry point."""

from __future__ import annotations

import click

from phyloai.cli.commands.pretree import pretree
from phyloai.cli.commands.tree import tree
from phyloai.cli.completion import completion
from phyloai.cli.doctor import doctor

CONTEXT_SETTINGS = {"help_option_names": ["-h", "--help"]}


@click.group(context_settings=CONTEXT_SETTINGS)
@click.pass_context
def cli(ctx: click.Context) -> None:
    """PhyloAI — modular phylogenomics analysis platform.

    Run 'phyloai doctor' to check your environment before starting.
    """
    ctx.ensure_object(dict)


cli.add_command(completion)
cli.add_command(doctor)
cli.add_command(pretree)
cli.add_command(tree)


if __name__ == "__main__":
    cli()
