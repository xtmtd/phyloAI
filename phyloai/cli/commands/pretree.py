"""Pre-tree CLI commands."""

from __future__ import annotations

import json
from pathlib import Path

import click
from rich.console import Console
from rich.progress import Progress

from phyloai.pretree.convert import convert_input, render_convert_summary_table
from phyloai.pretree.stats import (
    aggregate_summary,
    collect_seq_files,
    render_per_gene_table,
    render_single_file_panels,
    render_summary_table,
    stats_directory,
    stats_single_file,
)

console = Console()


def _fail(message: str, exit_code: int) -> None:
    click.echo(f"Error: {message}", err=True)
    raise click.exceptions.Exit(exit_code)


@click.group()
def pretree() -> None:
    """Pre-tree data preparation commands."""


@pretree.command(
    "convert",
    help=(
        "Normalize and convert one sequence file or a directory of sequence files. "
        "Directory input is the primary mode; --input may also be a single file. "
        "Invalid directory entries are skipped and summarized."
    ),
)
@click.option("--input", "input_path", type=click.Path(path_type=Path), required=True, help="Input directory or single sequence/alignment file.")
@click.option("--output-dir", "output_dir", type=click.Path(file_okay=False, path_type=Path), default=Path("runs/run001/pretree/convert"), show_default=True, help="Directory where converted files and result.json are written.")
@click.option("--to", "target_format", type=click.Choice(["fasta", "phylip-relaxed", "phylip-paml", "nexus"]), default="fasta", show_default=True, help="Target output format.")
@click.option("--input-format", type=click.Choice(["auto", "fasta", "phylip-relaxed", "phylip-paml", "nexus"]), default="auto", show_default=True, help="Override input format detection for all input files.")
@click.option("--seq-type", type=click.Choice(["AA", "NT", "auto"]), default="auto", show_default=True, help="Override sequence type detection.")
@click.option("--aa-special", type=click.Choice(["x", "keep"]), default="x", show_default=True, help="Convert B/Z/J/X/U/O to X, or preserve them with keep.")
@click.option("--threads", "threads", type=int, default=4, show_default=True, help="Directory mode worker count.")
@click.option("--quiet", "quiet", is_flag=True, default=False, help="Suppress Rich terminal output except errors.")
@click.option("--overwrite", "overwrite", is_flag=True, default=False, help="Delete and recreate a non-empty output directory before conversion.")
def convert_command(
    input_path: Path,
    output_dir: Path,
    target_format: str,
    input_format: str,
    seq_type: str,
    aa_special: str,
    threads: int,
    quiet: bool,
    overwrite: bool,
) -> None:
    if threads < 1:
        _fail("--threads must be at least 1.", 1)
    if not input_path.exists():
        _fail(f"Input path '{input_path}' does not exist.", 1)
    entries = [input_path] if input_path.is_file() else list(input_path.iterdir())
    conversion_error: str | None = None
    payload: dict | None = None
    if not quiet:
        with Progress(console=console, transient=True) as progress:
            task = progress.add_task("Converting sequence files", total=len(entries))
            try:
                payload = convert_input(
                    input_path,
                    output_dir,
                    target_format=target_format,
                    input_format=None if input_format == "auto" else input_format,
                    seq_type=None if seq_type == "auto" else seq_type,
                    aa_special=aa_special,
                    threads=threads,
                    overwrite=overwrite,
                    progress_callback=lambda _path: progress.advance(task),
                )
            except ValueError as exc:
                conversion_error = str(exc)
    else:
        try:
            payload = convert_input(
                input_path,
                output_dir,
                target_format=target_format,
                input_format=None if input_format == "auto" else input_format,
                seq_type=None if seq_type == "auto" else seq_type,
                aa_special=aa_special,
                threads=threads,
                overwrite=overwrite,
            )
        except ValueError as exc:
            conversion_error = str(exc)
    if conversion_error is not None:
        _fail(conversion_error, 1)
    if not quiet:
        console.print(render_convert_summary_table(payload["data"]["summary"]))
    result_path = output_dir / "result.json"
    result_path.parent.mkdir(parents=True, exist_ok=True)
    with open(result_path, "w") as fh:
        json.dump(payload, fh, indent=2)
    if not quiet:
        click.echo(f"Converted files saved to {output_dir / 'seqs'}", err=True)
        click.echo(f"Results saved to {result_path}", err=True)


@pretree.command(
    "stats",
    help=(
        "Inspect one sequence file or summarize a directory of sequence files. "
        "Use this command for quick pre-analysis QC on aligned or unaligned data. "
        "Exactly one of --seq or --seq-dir is required. Directory-only options: "
        "--per-gene, --threads. Shared options: --output-dir, --input-format, "
        "--seq-type, --quiet."
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
    help="Directory mode only. Include per-gene results in terminal output, or write them to per-gene.csv in the output directory.",
)
@click.option(
    "--per-gene-format",
    type=click.Choice(["csv", "tsv"]),
    default="csv",
    show_default=True,
    help="Directory mode only. Table format for the per-gene file written with --per-gene.",
)
@click.option(
    "--output-dir",
    "output_dir",
    type=click.Path(file_okay=False, path_type=Path),
    default=Path("runs/run001/pretree/stats"),
    show_default=True,
    help="Directory where result.json and per-gene files are written.",
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
@click.option(
    "--overwrite",
    "overwrite",
    is_flag=True,
    default=False,
    help="Delete and recreate a non-empty output directory before writing results.",
)
def stats_command(
    seq_dir: Path | None,
    seq: Path | None,
    per_gene: bool,
    per_gene_format: str,
    output_dir: Path,
    input_format: str | None,
    seq_type: str | None,
    threads: int,
    quiet: bool,
    overwrite: bool,
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

    if output_dir.exists() and any(output_dir.iterdir()):
        if not overwrite:
            _fail(f"Output directory '{output_dir}' already exists and is non-empty. Use --overwrite to replace it.", 1)
        import shutil
        shutil.rmtree(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    try:
        _run_stats_command(
            seq_dir,
            seq,
            per_gene,
            per_gene_format,
            output_dir,
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
    output_dir: Path,
    input_format: str | None,
    seq_type: str | None,
    threads: int,
    quiet: bool,
) -> None:
    """Run the stats command after CLI validation."""
    result_path = output_dir / "result.json"

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
        if not quiet:
            for panel in render_single_file_panels(stats):
                console.print(panel)
        with open(result_path, "w") as fh:
            json.dump(payload, fh, indent=2)
        click.echo(f"Results saved to {result_path}", err=True)
        return

    if not quiet:
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
    data: dict = {"summary": summary}
    if per_gene:
        data["per_gene"] = results
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
        "data": data,
    }
    if summary["n_genes_ok"] == 0:
        _fail("All files failed during processing.", 2)
    if not quiet:
        console.print(render_summary_table(summary))
        if per_gene:
            console.print(render_per_gene_table(results))
    with open(result_path, "w") as fh:
        json.dump(payload, fh, indent=2)
    click.echo(f"Summary saved to {result_path}", err=True)
    if per_gene:
        per_gene_path = output_dir / f"per-gene.{per_gene_format}"
        _write_per_gene_csv(results, per_gene_path, per_gene_format)
        click.echo(f"Per-gene table saved to {per_gene_path}", err=True)
    for warning in warnings:
        click.echo(warning, err=True)


def _write_per_gene_csv(results: list[dict], path: Path, fmt: str) -> None:
    """Write per-gene results to CSV or TSV."""
    import csv
    if not results:
        return
    columns = list(results[0].keys())
    delimiter = "\t" if fmt == "tsv" else ","
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=columns, delimiter=delimiter)
        writer.writeheader()
        writer.writerows(results)
