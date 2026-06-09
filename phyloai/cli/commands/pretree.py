"""Pre-tree CLI commands."""

from __future__ import annotations

import json
from pathlib import Path

import click
from rich.console import Console
from rich.progress import Progress

from phyloai.pretree.stats import (
    aggregate_summary,
    collect_seq_files,
    per_gene_output_path,
    render_per_gene_table,
    render_single_file_panels,
    render_summary_table,
    stats_directory,
    stats_single_file,
    write_output,
    write_per_gene_output,
)

console = Console()


def _fail(message: str, exit_code: int) -> None:
    click.echo(f"Error: {message}", err=True)
    raise click.exceptions.Exit(exit_code)


@click.group()
def pretree() -> None:
    """Pre-tree data preparation commands."""


@pretree.command(
    "stats",
    help=(
        "Inspect one sequence file or summarize a directory of sequence files. "
        "Use this command for quick pre-analysis QC on aligned or unaligned data. "
        "Exactly one of --seq or --seq-dir is required. Directory-only options: "
        "--per-gene, --threads. Shared options: --output, --output-format, "
        "--input-format, --seq-type, --quiet."
    ),
)
@click.option(
    "--seq-dir",
    type=click.Path(file_okay=False, path_type=Path),
    default=None,
    help="Directory mode. Scan all supported sequence/alignment files in a folder and compute a dataset summary.",
)
@click.option(
    "--seq",
    type=click.Path(dir_okay=False, path_type=Path),
    default=None,
    help="Single-file mode. Inspect one sequence or alignment file in detail.",
)
@click.option(
    "--per-gene",
    is_flag=True,
    default=False,
    help="Directory mode only. Include per-gene results in terminal output, or write them to an adjacent .per-gene.csv/.per-gene.tsv table when --output is used.",
)
@click.option(
    "--per-gene-format",
    type=click.Choice(["csv", "tsv"]),
    default="csv",
    show_default=True,
    help="Directory mode only. Table format for the adjacent per-gene file written with --per-gene --output.",
)
@click.option(
    "--output",
    "output_path",
    type=click.Path(path_type=Path),
    default=None,
    help="Write full results to a file. In directory mode, --per-gene writes an adjacent table file instead of mixing it into this file.",
)
@click.option(
    "--output-format",
    type=click.Choice(["text", "json"]),
    default="text",
    show_default=True,
    help="Terminal output format. Use json for machine-readable stdout.",
)
@click.option(
    "--input-format",
    type=click.Choice(["fasta", "phylip-relaxed", "nexus"]),
    default=None,
    help="Override auto-detection when a file suffix is misleading, for example a FASTA file without a standard extension. Choices: [fasta|phylip-relaxed|nexus].",
)
@click.option(
    "--seq-type",
    type=click.Choice(["AA", "NT"]),
    default=None,
    help="Override automatic molecule-type detection and force amino-acid or nucleotide rules.",
)
@click.option(
    "--threads",
    "threads",
    type=int,
    default=4,
    show_default=True,
    help="Directory mode only. Number of worker processes used for per-file statistics.",
)
@click.option(
    "--quiet",
    is_flag=True,
    default=False,
    help="Suppress terminal output except for errors.",
)
def stats_command(
    seq_dir: Path | None,
    seq: Path | None,
    per_gene: bool,
    per_gene_format: str,
    output_path: Path | None,
    output_format: str,
    input_format: str | None,
    seq_type: str | None,
    threads: int,
    quiet: bool,
) -> None:
    """Compute statistics for sequence or alignment files."""
    if bool(seq_dir) == bool(seq):
        if seq_dir and seq:
            _fail("--seq-dir and --seq are mutually exclusive.", 1)
        _fail("One of --seq or --seq-dir must be provided.", 1)

    if threads < 1:
        _fail("--threads must be at least 1.", 1)

    if seq is not None and per_gene:
        _fail("--per-gene is directory mode only; use --seq-dir to request per-gene output.", 1)

    if seq is not None and not seq.exists():
        _fail(f"Input file '{seq}' does not exist.", 1)
    if seq_dir is not None and not seq_dir.exists():
        _fail(f"Input directory '{seq_dir}' does not exist.", 1)

    try:
        _run_stats_command(
            seq_dir,
            seq,
            per_gene,
            per_gene_format,
            output_path,
            output_format,
            input_format,
            seq_type,
            threads,
            quiet,
        )
    except ValueError as exc:
        _fail(str(exc), 1)


def _run_stats_command(
    seq_dir: Path | None,
    seq: Path | None,
    per_gene: bool,
    per_gene_format: str,
    output_path: Path | None,
    output_format: str,
    input_format: str | None,
    seq_type: str | None,
    threads: int,
    quiet: bool,
) -> None:
    """Run the stats command after CLI validation."""

    if seq is not None:
        stats = stats_single_file(seq, seq_type=seq_type, input_format=input_format)
        payload = {
            "status": "success",
            "command": "phyloai pretree stats",
            "wall_time": 0.0,
            "tool_versions": {},
            "params": {
                "seq": str(seq),
                "seq_type": seq_type or "auto",
                "input_format": input_format or "auto",
            },
            "key_results": {},
            "error": None,
            "data": stats,
        }
        if output_path is not None:
            write_output(payload, output_path, mode="single", force_json=output_format == "json")
        if quiet:
            return
        if output_format == "json":
            click.echo(json.dumps(payload, indent=2, sort_keys=True))
            if output_path is not None:
                click.echo(f"Results written to {output_path}", err=True)
            return
        for panel in render_single_file_panels(stats):
            console.print(panel)
        if output_path is not None:
            click.echo(f"Results written to {output_path}")
        return

    if not quiet and output_format == "text":
        files = collect_seq_files(seq_dir)
        with Progress(console=console, transient=True) as progress:
            task = progress.add_task("Processing sequence files", total=len(files))
            results, warnings = stats_directory(
                seq_dir,
                seq_type=seq_type,
                input_format=input_format,
                threads=threads,
                progress_callback=lambda _path: progress.advance(task),
            )
    else:
        results, warnings = stats_directory(seq_dir, seq_type=seq_type, input_format=input_format, threads=threads)
    summary = aggregate_summary(results)
    summary["warnings"] = warnings
    payload = {
        "status": "success" if summary["n_genes_ok"] > 0 else "error",
        "command": "phyloai pretree stats",
        "wall_time": 0.0,
        "tool_versions": {},
        "params": {
            "seq_dir": str(seq_dir),
            "seq_type": seq_type or "auto",
            "input_format": input_format or "auto",
            "threads": threads,
            "per_gene": per_gene,
        },
        "key_results": {},
        "error": None if summary["n_genes_ok"] > 0 else "All files failed during processing.",
        "data": {
            "summary": summary,
            "per_gene": results,
        },
    }
    if output_path is not None:
        write_output(payload, output_path, mode="directory", per_gene=per_gene, force_json=output_format == "json")
        per_gene_path = per_gene_output_path(output_path, per_gene_format) if per_gene else None
        if per_gene_path is not None:
            write_per_gene_output(payload, per_gene_path)
    else:
        per_gene_path = None
    if summary["n_genes_ok"] == 0:
        _fail("All files failed during processing.", 2)
    if quiet:
        return
    if output_format == "json":
        click.echo(json.dumps(payload, indent=2, sort_keys=True))
        if output_path is not None:
            if per_gene:
                click.echo(f"Summary saved to {output_path}", err=True)
                click.echo(f"Per-gene table saved to {per_gene_path}", err=True)
            else:
                click.echo(f"Summary saved to {output_path}", err=True)
        return
    console.print(render_summary_table(summary))
    if per_gene and output_path is None:
        console.print(render_per_gene_table(results))
    if output_path is not None:
        if per_gene:
            click.echo(f"Summary saved to {output_path}")
            click.echo(f"Per-gene table saved to {per_gene_path}")
        else:
            click.echo(f"Summary saved to {output_path}")
    for warning in warnings:
        click.echo(warning, err=True)
