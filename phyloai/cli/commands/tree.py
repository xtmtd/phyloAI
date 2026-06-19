"""Tree inference CLI commands."""

from __future__ import annotations

import json
import os
from pathlib import Path

import click
from rich.console import Console
from rich.progress import Progress

from phyloai.tree.ml import run_fasttree

console = Console()


def _fail(message: str, exit_code: int) -> None:
    click.echo(f"Error: {message}", err=True)
    raise click.exceptions.Exit(exit_code)


class _MLGroup(click.Group):
    def list_commands(self, ctx: click.Context) -> list[str]:
        return ["fasttree"]


class _TreeGroup(click.Group):
    def list_commands(self, ctx: click.Context) -> list[str]:
        return ["ml"]


@click.group(cls=_TreeGroup)
def tree() -> None:
    """Phylogenetic tree inference commands."""


@tree.group(cls=_MLGroup)
def ml() -> None:
    """Maximum-likelihood tree inference (FastTree / IQ-TREE3)."""


@ml.command(
    "fasttree",
    help=(
        "Infer ML trees using FastTree.\n\n"
        "  --msa-dir : batch gene trees from an MSA directory\n"
        "  --matrix  : single supermatrix tree from one file\n\n"
        "FastTree natively reads FASTA and phylip-relaxed formats.\n"
        "NEXUS files must be converted first via 'phyloai pretree convert'."
    ),
)
@click.option(
    "--msa-dir",
    type=click.Path(file_okay=False, path_type=Path),
    default=None,
    help="Directory of MSA files for batch gene tree inference.",
)
@click.option(
    "--matrix",
    type=click.Path(dir_okay=False, path_type=Path),
    default=None,
    help="Single concatenated matrix file for supermatrix tree inference.",
)
@click.option(
    "--seq-type",
    type=click.Choice(["AA", "NT", "auto"]),
    default="auto",
    show_default=True,
    help="Molecule type.",
)
@click.option(
    "--model",
    type=click.Choice(["jtt", "lg", "wag", "jc", "gtr"]),
    default=None,
    help="Substitution model. AA: jtt/lg/wag. NT: jc/gtr. AA default: lg. NT default: gtr.",
)
@click.option(
    "--mode",
    type=click.Choice(["normal", "fastest", "slow"]),
    default="normal",
    show_default=True,
    help="Speed/accuracy trade-off.",
)
@click.option(
    "--boot",
    type=click.IntRange(min=0),
    default=1000,
    show_default=True,
    help="Bootstrap replicates. 0 disables node support.",
)
@click.option(
    "--cat",
    type=click.IntRange(min=1),
    default=20,
    show_default=True,
    help="Number of rate categories for FastTree (-cat N).",
)
@click.option(
    "--gamma",
    is_flag=True,
    default=True,
    help="Enable gamma-distributed rate heterogeneity (default: on).",
)
@click.option(
    "--output-dir", "-o",
    type=click.Path(file_okay=False, path_type=Path),
    default=Path("runs/tree/ml/fasttree"),
    show_default=True,
    help="Output directory.",
)
@click.option(
    "--threads", "-t",
    type=int,
    default=4,
    show_default=True,
    help="Parallel gene tree workers. Only used in --msa-dir mode.",
)
@click.option(
    "--fasttree-path",
    type=Path,
    default=None,
    help="Explicit path to FastTree executable.",
)
@click.option(
    "--tool-args",
    type=str,
    default=None,
    help="Extra FastTree flags appended verbatim (e.g. '-spr 4 -mlacc 2 -slownni'). "
    "Managed I/O and data-type flags (-nt) are blocked. "
    "Strategy flags like -lg/-wag/-cat/-gamma/-boot may override PhyloAI defaults.",
)
@click.option("--overwrite", is_flag=True, default=False, help="Overwrite existing output directory.")
@click.option("--resume", is_flag=True, default=False, help="Resume from checkpoint (--msa-dir only).")
@click.option("--dry-run", is_flag=True, default=False, help="Show commands without executing.")
@click.option("--quiet", "-q", is_flag=True, default=False, help="Suppress terminal output except errors.")
def fasttree_command(
    msa_dir: Path | None,
    matrix: Path | None,
    seq_type: str,
    model: str | None,
    mode: str,
    boot: int,
    cat: int,
    gamma: bool,
    output_dir: Path,
    threads: int,
    fasttree_path: Path | None,
    tool_args: str | None,
    overwrite: bool,
    resume: bool,
    dry_run: bool,
    quiet: bool,
) -> None:
    batch_mode = msa_dir is not None
    single_mode = matrix is not None

    if batch_mode == single_mode:
        if not batch_mode and not single_mode:
            _fail("Either --msa-dir or --matrix is required.", 1)
        else:
            _fail("--msa-dir and --matrix are mutually exclusive.", 1)

    if threads < 1:
        _fail("--threads must be at least 1.", 1)

    if resume and overwrite:
        _fail("--overwrite and --resume are mutually exclusive.", 1)

    if resume and single_mode:
        _fail("--resume is only supported in --msa-dir mode.", 1)

    # Validate input paths exist
    if msa_dir is not None and not msa_dir.exists():
        _fail(f"--msa-dir does not exist: {msa_dir}", 1)
    if matrix is not None and not matrix.exists():
        _fail(f"--matrix does not exist: {matrix}", 1)

    # Validate fasttree-path
    if fasttree_path is not None:
        if not fasttree_path.exists():
            _fail(f"--fasttree-path does not exist: {fasttree_path}", 1)
        if not os.access(str(fasttree_path), os.X_OK):
            _fail(f"--fasttree-path is not executable: {fasttree_path}", 1)

    # Warn about --threads in single mode
    if single_mode and threads != 4:
        if not quiet:
            click.echo("Warning: --threads has no effect in single --matrix mode.", err=True)

    fasttree_path_str = str(fasttree_path) if fasttree_path else None

    def _invoke(progress_callback=None):
        return run_fasttree(
            msa_dir=msa_dir,
            matrix=matrix,
            output_dir=output_dir,
            seq_type=seq_type,
            model=model,
            mode=mode,
            boot=boot,
            cat=cat,
            gamma=gamma,
            threads=threads,
            fasttree_path=fasttree_path_str,
            tool_args=tool_args,
            overwrite=overwrite,
            resume=resume,
            dry_run=dry_run,
            quiet=quiet,
            progress_callback=progress_callback,
        )

    error_msg: str | None = None

    try:
        if not quiet and not dry_run and batch_mode:
            from phyloai.tree.ml import _scan_input

            found, _ = _scan_input(msa_dir)
            total = len(found)

            if total == 0:
                _fail("No valid input files found in --msa-dir.", 1)

            with Progress(console=console, transient=True) as progress:
                task = progress.add_task(
                    "[cyan]Inferring gene trees with FastTree", total=total
                )
                try:
                    payload = _invoke(
                        progress_callback=lambda _: progress.advance(task)
                    )
                    n_resume = payload["data"]["summary"].get("n_resume_skipped", 0)
                    if n_resume > 0:
                        progress.advance(task, advance=n_resume)
                except (ValueError, FileNotFoundError) as exc:
                    error_msg = str(exc)
        else:
            try:
                payload = _invoke()
            except (ValueError, FileNotFoundError) as exc:
                error_msg = str(exc)
    except SystemExit:
        raise
    except Exception as exc:
        error_msg = str(exc)

    if error_msg is not None:
        exit_code = 3 if "not found" in error_msg.lower() else 1
        _fail(error_msg, exit_code)

    if dry_run:
        if not quiet:
            click.echo(
                f"Dry run: {payload['data']['summary']['n_input_files']} input(s) would be processed."
            )
            for item in payload["data"].get("files", []):
                if "cmd" in item:
                    click.echo(" ".join(item["cmd"]))
        return

    result_path = output_dir / "result.json"
    result_path.parent.mkdir(parents=True, exist_ok=True)
    with open(result_path, "w") as fh:
        json.dump(payload, fh, indent=2)

    summary = payload["data"]["summary"]
    n_failed = summary.get("n_failed", 0)
    n_trees = summary.get("n_trees", 0)
    n_skipped = summary.get("n_skipped", 0)

    if not quiet:
        click.echo(
            f"Trees: {n_trees} | Failed: {n_failed} | Skipped: {n_skipped}"
        )
        if batch_mode:
            click.echo(f"Trees saved to {output_dir / 'trees'}", err=True)
            click.echo(f"Logs saved to {output_dir / 'logs'}", err=True)
        click.echo(f"Results saved to {result_path}", err=True)

        if n_failed > 0:
            click.echo(
                f"Warning: {n_failed} gene(s) failed. Check result.json data.failed for details.",
                err=True,
            )

    if n_trees == 0 and (n_failed > 0):
        _fail("All FastTree runs failed.", 2)
