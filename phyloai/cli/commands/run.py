"""One-click phylogenomics pipeline."""

from __future__ import annotations

from pathlib import Path

import click
from rich.console import Console

console = Console()


@click.command(
    "run",
    help=(
        "One-click phylogenomics pipeline from raw sequences to a species tree.\n\n"
        "Runs the full pipeline automatically using sensible defaults. "
        "For fine-grained control over individual steps, use the constituent "
        "subcommands (phyloai pretree align, phyloai tree ml iqtree, etc.).\n\n"
        "Modes:\n\n"
        "  supermatrix  convert -> align -> trim -> [filter] -> concat -> iqtree\n\n"
        "  supertree    convert -> align -> trim -> [filter] -> gene trees -> wastral\n\n"
        "The [filter] step (TAPER error-site masking) is included in --speed normal "
        "and skipped in --speed fast.\n\n"
        "Speed modes:\n\n"
        "  normal  MAFFT linsi, trimAl -automated1, TAPER filter, IQ-TREE3 / FastTree\n\n"
        "  fast    MAFFT auto, trimAl -automated1, no filter, FastTree\n\n"
        "Examples:\n\n"
        "  phyloai run --seq-dir ./markers --mode supermatrix\n\n"
        "  phyloai run --seq-dir ./markers --mode supertree --speed fast --threads 16\n\n"
        "  phyloai run --seq-dir ./markers --mode supermatrix --resume"
    ),
)
@click.option(
    "--seq-dir",
    required=True,
    type=click.Path(exists=True, file_okay=False, path_type=Path),
    help="Input sequence directory (any format; always converted first).",
)
@click.option(
    "--mode",
    type=click.Choice(["supermatrix", "supertree"]),
    default="supermatrix",
    show_default=True,
    help="Pipeline mode: supermatrix (concat -> iqtree) or supertree (gene trees -> wastral).",
)
@click.option(
    "--speed",
    type=click.Choice(["normal", "fast"]),
    default="normal",
    show_default=True,
    help=(
        "Speed/accuracy trade-off. normal: MAFFT linsi + TAPER + IQ-TREE3. "
        "fast: MAFFT auto, no TAPER, FastTree."
    ),
)
@click.option(
    "--output-dir",
    "-o",
    type=click.Path(path_type=Path),
    default=Path("runs/run"),
    show_default=True,
    help="Root output directory for all pipeline steps.",
)
@click.option(
    "--threads",
    "-t",
    type=int,
    default=4,
    show_default=True,
    help="Thread count passed to all steps.",
)
@click.option("--resume", is_flag=True, default=False, help="Resume from run_checkpoint.json.")
@click.option("--overwrite", is_flag=True, default=False, help="Delete and recreate output directory.")
@click.option("--dry-run", is_flag=True, default=False, help="Show steps without running.")
@click.option("--quiet", "-q", is_flag=True, default=False, help="Suppress non-error output.")
def run(
    seq_dir: Path,
    mode: str,
    speed: str,
    output_dir: Path,
    threads: int,
    resume: bool,
    overwrite: bool,
    dry_run: bool,
    quiet: bool,
) -> None:
    from phyloai.cli.commands._run_pipeline import execute_pipeline
    execute_pipeline(
        seq_dir=seq_dir,
        mode=mode,
        speed=speed,
        output_dir=output_dir,
        threads=threads,
        resume=resume,
        overwrite=overwrite,
        dry_run=dry_run,
        quiet=quiet,
    )
