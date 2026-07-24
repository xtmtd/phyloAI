"""PhyloAI CLI entry point."""

from __future__ import annotations

import click

from phyloai import __version__
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
        return ["completion", "doctor", "update", "run", "pretree", "tree", "posttree", "report", "mcp-server"]


@click.group(context_settings=CONTEXT_SETTINGS, cls=_RootGroup)
@click.version_option(__version__, "-v", "--version")
@click.pass_context
def cli(ctx: click.Context) -> None:
    """PhyloAI — modular phylogenomics analysis platform.

    Run 'phyloai doctor' to check your environment before starting.
    """
    ctx.ensure_object(dict)


@cli.command("update")
@click.option("--check", is_flag=True, default=False, help="Only check for updates; do not install.")
@click.option("--yes", "-y", is_flag=True, default=False, help="Update without confirmation prompt.")
def update_command(check: bool, yes: bool) -> None:
    """Check for and install PhyloAI updates from GitHub releases."""
    from phyloai.core.update import check_update, run_update

    if check:
        result = check_update()
        if result["status"] == "error":
            click.echo(f"Error: {result['message']}", err=True)
            raise SystemExit(1)
        if result["status"] == "up_to_date":
            click.echo(f"PhyloAI is up to date (v{result['current']}).")
        else:
            click.echo(f"Update available: v{result['current']} -> v{result['latest']}")
            raise SystemExit(1)
    else:
        ret = run_update(confirm=yes)
        if ret != 0:
            raise SystemExit(ret)


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
