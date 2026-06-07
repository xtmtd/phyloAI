"""PhyloAI CLI entry point."""

from __future__ import annotations
from pathlib import Path
from typing import Optional

import click
import yaml

from phyloai.cli.doctor import doctor

CONTEXT_SETTINGS = {"help_option_names": ["-h", "--help"]}


@click.group(context_settings=CONTEXT_SETTINGS)
@click.option(
    "--config", "config_file",
    type=click.Path(exists=True, path_type=Path),
    default=None,
    help="Load parameters from a YAML config file. CLI flags override config values.",
)
@click.pass_context
def cli(ctx: click.Context, config_file: Optional[Path]) -> None:
    """PhyloAI — modular phylogenomics analysis platform.

    Run 'phyloai doctor' to check your environment before starting.
    """
    ctx.ensure_object(dict)
    if config_file:
        with open(config_file) as fh:
            ctx.obj["config"] = yaml.safe_load(fh)
    else:
        ctx.obj["config"] = {}


cli.add_command(doctor)
