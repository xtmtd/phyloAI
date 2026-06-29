"""phyloai doctor — environment detection command."""

from __future__ import annotations
import json

import click
from rich.console import Console
from rich.table import Table

from phyloai.core.env import TOOL_GROUPS, TOOL_REGISTRY, ToolEnv, ToolStatus

console = Console()
DISPLAY_NAMES = {"bmge": "BMGE.jar"}

STATUS_ICON = {
    ToolStatus.OK:      "[green]OK[/green]",
    ToolStatus.WARN:    "[yellow]WARN[/yellow]",
    ToolStatus.MISSING: "[red]MISSING[/red]",
}


def _build_table(title: str, tools: dict) -> Table:
    table = Table(title=title, show_header=True)
    table.add_column("Status", width=8)
    table.add_column("Tool", width=14)
    table.add_column("Version", width=12)
    table.add_column("Path / Note")
    for name, info in tools.items():
        table.add_row(
            STATUS_ICON[info.status],
            DISPLAY_NAMES.get(name, name),
            info.version or "\u2014",
            str(info.path) if info.path else info.note,
        )
    return table


@click.command("doctor")
@click.option(
    "--output-format",
    type=click.Choice(["text", "json"]),
    default="text",
    help="Output format. Default: text.",
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
        if any(
            info.status == ToolStatus.MISSING and TOOL_REGISTRY.get(name, {}).get("required", False)
            for name, info in tools.items()
        ):
            ctx = click.get_current_context()
            ctx.exit(3)
        return

    # Text output: render grouped tables in TOOL_GROUPS order.
    # Tools not in any group appear last under "Other".
    covered: set[str] = set()
    for group_name, member_names in TOOL_GROUPS:
        group_tools = {name: tools[name] for name in member_names if name in tools}
        covered.update(group_tools.keys())
        if group_tools:
            console.print(_build_table(group_name, group_tools))

    other = {name: info for name, info in tools.items() if name not in covered}
    if other:
        console.print(_build_table("Other", other))

    if any(
        info.status == ToolStatus.MISSING and TOOL_REGISTRY.get(name, {}).get("required", False)
        for name, info in tools.items()
    ):
        ctx = click.get_current_context()
        ctx.exit(3)
