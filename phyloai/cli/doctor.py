"""phyloai doctor — environment detection command."""

from __future__ import annotations
import json

import click
from rich.console import Console
from rich.table import Table

from phyloai.core.env import ToolEnv, ToolStatus

console = Console()


@click.command("doctor")
@click.option(
    "--output-format",
    type=click.Choice(["text", "json"]),
    default="text",
    help="Output format.",
)
def doctor(output_format: str) -> None:
    """Check that required external tools are available."""
    env = ToolEnv()
    tools = env.check_all()

    if output_format == "json":
        out = {
            name: {
                "status": info.status.value,
                "path": str(info.path) if info.path else None,
                "version": info.version,
                "note": info.note,
            }
            for name, info in tools.items()
        }
        click.echo(json.dumps(out, indent=2))
        return

    table = Table(title="PhyloAI Environment Check", show_header=True)
    table.add_column("Status", width=8)
    table.add_column("Tool", width=14)
    table.add_column("Version", width=12)
    table.add_column("Path / Note")

    status_icon = {
        ToolStatus.OK:      "[green]OK[/green]",
        ToolStatus.WARN:    "[yellow]WARN[/yellow]",
        ToolStatus.MISSING: "[red]MISSING[/red]",
    }

    for name, info in tools.items():
        table.add_row(
            status_icon[info.status],
            name,
            info.version or "—",
            str(info.path) if info.path else info.note,
        )

    console.print(table)
