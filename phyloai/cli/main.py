"""PhyloAI CLI entry point."""

from __future__ import annotations

import click

from phyloai.cli.commands.posttree import posttree
from phyloai.cli.commands.pretree import pretree
from phyloai.cli.commands.report import report
from phyloai.cli.commands.run import run
from phyloai.cli.commands.tree import tree
from phyloai.cli.commands.mcp_server import mcp_server
from phyloai.cli.completion import completion
from phyloai.cli.doctor import doctor

CONTEXT_SETTINGS = {"help_option_names": ["-h", "--help"]}


class _RootGroup(click.Group):
    def list_commands(self, ctx: click.Context) -> list[str]:
        return ["completion", "doctor", "run", "pretree", "tree", "posttree", "report", "mcp-server"]


@click.group(context_settings=CONTEXT_SETTINGS, cls=_RootGroup)
@click.pass_context
def cli(ctx: click.Context) -> None:
    """PhyloAI — modular phylogenomics analysis platform.

    Run 'phyloai doctor' to check your environment before starting.
    """
    ctx.ensure_object(dict)


cli.add_command(completion)
cli.add_command(doctor)
cli.add_command(run)
cli.add_command(pretree)
cli.add_command(tree)
cli.add_command(posttree)
cli.add_command(report)
cli.add_command(mcp_server)


if __name__ == "__main__":
    cli()
