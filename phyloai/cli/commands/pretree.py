"""Pre-tree CLI commands."""

from __future__ import annotations

import json
import shlex
import time
from pathlib import Path

import click
from rich.console import Console
from rich.progress import Progress

from phyloai.pretree.align import render_align_summary_table, run_align
from phyloai.pretree.convert import convert_input, render_convert_summary_table
from phyloai.pretree.stats import (
    aggregate_summary,
    collect_seq_files,
    render_single_file_panels,
    render_summary_table,
    stats_directory,
    stats_single_file,
)
from phyloai.pretree.trim import render_trim_summary_table, run_trim, _scan_input as _trim_scan_input
from phyloai.pretree.concat import run_concat, _render_concat_panels, _build_concat_command
from phyloai.pretree.metrics import (
    _compute_correlation,
    _detect_input_delimiter,
    _generate_all_plots,
    _generate_basic_statistics,
    _generate_correlation_heatmap,
    _get_delimiter,
    _plot_single_metric,
    _select_correlation_columns,
    _table_suffix,
    _write_correlation_csv,
    run_metrics,
)

console = Console()


def _fail(message: str, exit_code: int) -> None:
    click.echo(f"Error: {message}", err=True)
    raise click.exceptions.Exit(exit_code)


class _PretreeGroup(click.Group):
    def list_commands(self, ctx: click.Context) -> list[str]:
        return ["convert", "stats", "align", "trim", "metrics", "filter", "concat"]


@click.group(cls=_PretreeGroup)
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
@click.option("--output-dir", "-o", "output_dir", type=click.Path(file_okay=False, path_type=Path), default=Path("runs/pretree/convert"), show_default=True, help="Directory where converted files and result.json are written.")
@click.option("--to", "target_format", type=click.Choice(["fasta", "phylip-relaxed", "phylip-paml", "nexus"]), default="fasta", show_default=True, help="Target output format.")
@click.option("--input-format", type=click.Choice(["auto", "fasta", "phylip-relaxed", "phylip-paml", "nexus"]), default="auto", show_default=True, help="Override input format detection for all input files.")
@click.option("--seq-type", type=click.Choice(["AA", "NT", "auto"]), default="auto", show_default=True, help="Override sequence type detection.")
@click.option("--aa-special", type=click.Choice(["x", "keep"]), default="x", show_default=True, help="Convert B/Z/J/X/U/O to X, or preserve them with keep.")
@click.option("--threads", "-t", "threads", type=int, default=4, show_default=True, help="Directory mode worker count.")
@click.option("--quiet", "-q", "quiet", is_flag=True, default=False, help="Suppress Rich terminal output except errors.")
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
                    quiet=quiet,
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
                quiet=quiet,
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
        "--per-gene, --table-format, --threads. Shared options: --output-dir, "
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
    "--unaligned",
    "unaligned",
    is_flag=True,
    default=False,
    help="Treat input as unaligned sequences (excludes alignment_length and site-pattern columns from per-gene CSV).",
)
@click.option(
    "--per-gene",
    is_flag=True,
    default=False,
    help="Directory mode only. Write per-gene results to a per-gene CSV/TSV file in the output directory.",
)
@click.option(
    "--table-format",
    type=click.Choice(["csv", "tsv"]),
    default="csv",
    show_default=True,
    help="Directory mode only. Table format for the per-gene file written with --per-gene.",
)
@click.option(
    "--output-dir",
    "-o",
    "output_dir",
    type=click.Path(file_okay=False, path_type=Path),
    default=Path("runs/pretree/stats"),
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
    "-t",
    "threads",
    type=int,
    default=4,
    show_default=True,
    help="Directory mode only. Number of worker processes used for per-file statistics.",
)
@click.option(
    "--quiet",
    "-q",
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
    table_format: str,
    output_dir: Path,
    input_format: str | None,
    seq_type: str | None,
    threads: int,
    quiet: bool,
    overwrite: bool,
    unaligned: bool,
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
            table_format,
            output_dir,
            input_format,
            seq_type,
            threads,
            quiet,
            is_aligned=not unaligned,
            unaligned=unaligned,
            overwrite=overwrite,
        )
    except ValueError as exc:
        _fail(str(exc), 1)


def _build_stats_command(
    seq_dir: Path | None,
    seq: Path | None,
    per_gene: bool,
    table_format: str,
    output_dir: Path,
    input_format: str | None,
    seq_type: str | None,
    threads: int,
    unaligned: bool,
    overwrite: bool,
) -> str:
    parts = ["phyloai", "pretree", "stats"]
    if seq is not None:
        parts.extend(["--seq", str(seq)])
    else:
        parts.extend(["--seq-dir", str(seq_dir)])
    if unaligned:
        parts.append("--unaligned")
    if per_gene:
        parts.append("--per-gene")
    parts.extend(["--table-format", table_format])
    parts.extend(["--output-dir", str(output_dir)])
    if input_format:
        parts.extend(["--input-format", input_format])
    if seq_type:
        parts.extend(["--seq-type", seq_type])
    parts.extend(["--threads", str(threads)])
    if overwrite:
        parts.append("--overwrite")
    return shlex.join(parts)


def _run_stats_command(
    seq_dir: Path | None,
    seq: Path | None,
    per_gene: bool,
    table_format: str,
    output_dir: Path,
    input_format: str | None,
    seq_type: str | None,
    threads: int,
    quiet: bool,
    is_aligned: bool = True,
    unaligned: bool = False,
    overwrite: bool = False,
) -> None:
    """Run the stats command after CLI validation."""
    start = time.monotonic()
    output_dir = output_dir.resolve()
    result_path = output_dir / "result.json"
    params: dict = {
        "seq_dir": str(seq_dir) if seq_dir is not None else None,
        "seq": str(seq) if seq is not None else None,
        "per_gene": per_gene,
        "table_format": table_format,
        "output_dir": str(output_dir),
        "input_format": input_format or "auto",
        "seq_type": seq_type or "auto",
        "threads": threads,
        "quiet": quiet,
        "is_aligned": is_aligned,
        "unaligned": unaligned,
        "overwrite": overwrite,
    }

    if seq is not None:
        cmd_str = _build_stats_command(seq_dir, seq, per_gene, table_format, output_dir, input_format, seq_type, threads, unaligned, overwrite)
        stats = stats_single_file(seq, seq_type=seq_type, input_format=input_format)
        elapsed = time.monotonic() - start
        payload = {
            "status": "success",
            "command": cmd_str,
            "wall_time": round(elapsed, 3),
            "tool_versions": {},
            "params": params,
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

    cmd_str = _build_stats_command(seq_dir, seq, per_gene, table_format, output_dir, input_format, seq_type, threads, unaligned, overwrite)
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

    # Warn when the declared alignment mode disagrees with per-file detection.
    if seq_dir is not None:
        ok_results = [r for r in results if "error" not in r]
        n_detected_aligned = sum(1 for r in ok_results if r.get("is_aligned"))
        n_detected_unaligned = len(ok_results) - n_detected_aligned
        if is_aligned and n_detected_unaligned > 0:
            warnings.append(
                f"--unaligned is NOT set, but {n_detected_unaligned} of "
                f"{len(ok_results)} files were detected as unaligned "
                "(unequal sequence lengths).  Use --unaligned to write "
                "unaligned-specific columns to the per-gene table."
            )
        elif not is_aligned and n_detected_aligned > 0:
            warnings.append(
                f"--unaligned is set, but {n_detected_aligned} of "
                f"{len(ok_results)} files were detected as aligned "
                "(equal sequence lengths).  Drop --unaligned to include "
                "alignment-specific columns in the per-gene table."
            )

    summary = aggregate_summary(results)
    summary["warnings"] = warnings
    data: dict = {"summary": summary}
    output_files: dict[str, dict[str, str]] = {}
    if per_gene:
        per_gene_path = output_dir / f"per-gene.{table_format}"
        _write_per_gene_csv(results, per_gene_path, table_format, is_aligned=is_aligned)
        data["per_gene"] = results
        output_files["per_gene_table"] = {"path": str(per_gene_path), "description": "Per-locus sequence statistics: length, taxon count, gap ratio for each input file"}
        if not quiet:
            click.echo(f"Per-gene table saved to {per_gene_path}", err=True)
    data["output_files"] = output_files
    payload = {
        "status": "success" if summary["n_genes_ok"] > 0 else "error",
        "command": cmd_str,
        "wall_time": round(time.monotonic() - start, 3),
        "tool_versions": {},
        "params": params,
        "key_results": {},
        "error": None if summary["n_genes_ok"] > 0 else "All files failed during processing.",
        "data": data,
    }
    if summary["n_genes_ok"] == 0:
        _fail("All files failed during processing.", 2)
    if not quiet:
        console.print(render_summary_table(summary))
    with open(result_path, "w") as fh:
        json.dump(payload, fh, indent=2)
    click.echo(f"Results saved to {result_path}", err=True)
    for warning in warnings:
        click.echo(warning, err=True)


def _write_per_gene_csv(
    results: list[dict],
    path: Path,
    fmt: str,
    is_aligned: bool = True,
) -> None:
    """Write per-gene results to CSV or TSV.

    Uses the curated ``PER_GENE_COLUMNS`` set (from ``pretree.stats``) so
    that only well-known, scalar metrics appear in the output.  When
    *is_aligned* is ``True`` (default) the output includes alignment-
    specific columns (``alignment_length``, site patterns) and omits
    per-sequence length statistics; when ``False`` the reverse applies.
    This avoids the ``ValueError`` that ``csv.DictWriter`` raises when
    mixed aligned / unaligned rows contain differing key sets.
    """
    import csv

    from phyloai.pretree.stats import per_gene_columns_for_rows

    if not results:
        return

    columns = per_gene_columns_for_rows(results, is_aligned=is_aligned)

    delimiter = "\t" if fmt == "tsv" else ","
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=columns, delimiter=delimiter)
        writer.writeheader()
        normalized = [
            {col: row.get(col, "") for col in columns} for row in results
        ]
        writer.writerows(normalized)


@pretree.command(
    "align",
    help=(
        "Align a directory of unaligned sequences using MAFFT or MAGUS.\n\n"
        "Strategies: fftns1/fftns2 are fastest but least accurate; "
        "linsi/einsi/ginsi offer high accuracy; "
        "magus is slowest but best for large or difficult datasets (Linux only).\n\n"
        "Use --backtrans with --nt-dir to also produce codon-level NT alignments "
        "from a protein alignment using trimAl backtranslation.\n\n"
        "--threads controls how many genes are aligned in parallel; each "
        "individual alignment uses a single thread."
    ),
)
@click.option("--seq-dir", type=click.Path(file_okay=False, path_type=Path),
              required=True, help="Input directory of unaligned sequence files.")
@click.option("--method",
              type=click.Choice(["fftns1", "fftns2", "auto", "linsi", "einsi", "ginsi", "magus"]),
              default="linsi", show_default=True,
              help=(
                  "Alignment strategy. "
                  "fftns1/fftns2: fast, lower accuracy. "
                  "auto: MAFFT chooses strategy automatically. "
                  "linsi/einsi/ginsi: high accuracy. "
                  "magus: highest accuracy, slowest, best for large datasets (Linux only)."
              ))
@click.option("--seq-type", type=click.Choice(["AA", "NT", "auto"]), default="auto",
              show_default=True, help="Molecule type of input sequences. Auto-detects from first gene if 'auto'.")
@click.option("--backtrans", is_flag=True, default=False,
              help="Produce codon NT alignments via trimAl -backtrans. Requires --nt-dir.")
@click.option("--nt-dir", type=click.Path(file_okay=False, path_type=Path), default=None,
              help="Directory of unaligned CDS sequences for --backtrans mode.")
@click.option("--output-dir", "-o", type=click.Path(file_okay=False, path_type=Path),
              default=Path("runs/pretree/align"), show_default=True,
              help="Output directory; contains seqs/, logs/<locus>.log, result.json.")
@click.option("--threads", "-t", type=int, default=4, show_default=True,
              help="Number of genes to align in parallel (each uses 1 thread).")
@click.option("--tool-args", type=str, default=None,
              help="MAGUS strategy arguments only. PhyloAI manages input, output, work dir, datatype, and threads.")
@click.option("--mafft-path", type=click.Path(dir_okay=False, path_type=Path), default=None,
              help="Explicit MAFFT executable path for MAFFT methods; PATH lookup is used when omitted.")
@click.option("--magus-path", type=click.Path(dir_okay=False, path_type=Path), default=None,
              help="Explicit MAGUS executable path for --method magus; PATH lookup is used when omitted.")
@click.option("--trimal-path", type=click.Path(dir_okay=False, path_type=Path), default=None,
              help="Explicit trimAl executable path for --backtrans; PATH lookup is used when omitted.")
@click.option("--resume", is_flag=True, default=False,
              help=(
                  "Resume from checkpoint.json in the output directory. "
                  "Requires the same parameters as the original run. "
                  "Mutually exclusive with --overwrite."
              ))
@click.option("--overwrite", is_flag=True, default=False,
              help="Delete and recreate a non-empty output directory before running.")
@click.option("--dry-run", is_flag=True, default=False,
              help="Print commands without executing; creates no files.")
@click.option("--quiet", "-q", is_flag=True, default=False,
              help="Suppress Rich terminal output except errors.")
def align_command(
    seq_dir: Path,
    method: str,
    seq_type: str,
    backtrans: bool,
    nt_dir: Path | None,
    output_dir: Path,
    threads: int,
    tool_args: str | None,
    mafft_path: Path | None,
    magus_path: Path | None,
    trimal_path: Path | None,
    resume: bool,
    overwrite: bool,
    dry_run: bool,
    quiet: bool,
) -> None:
    if threads < 1:
        _fail("--threads must be at least 1.", 1)
    if not seq_dir.exists():
        _fail(f"--seq-dir '{seq_dir}' does not exist.", 1)
    if nt_dir is not None and not nt_dir.exists():
        _fail(f"--nt-dir '{nt_dir}' does not exist.", 1)

    payload: dict | None = None
    error_msg: str | None = None

    def _invoke(progress_callback=None):
        return run_align(
            seq_dir=seq_dir,
            output_dir=output_dir,
            method=method,
            seq_type=seq_type,
            backtrans=backtrans,
            nt_dir=nt_dir,
            threads=threads,
            tool_args=tool_args,
            mafft_path=mafft_path,
            magus_path=magus_path,
            trimal_path=trimal_path,
            overwrite=overwrite,
            resume=resume,
            dry_run=dry_run,
            quiet=quiet,
            progress_callback=progress_callback,
        )

    if not quiet and not dry_run:
        if resume and overwrite:
            _fail("--overwrite and --resume are mutually exclusive.", 1)
        if resume:
            ckpt_path = output_dir / "checkpoint.json"
            if ckpt_path.exists():
                from phyloai.core.checkpoint import load_checkpoint
                from phyloai.pretree.checkpoint_helpers import plan_resume
                from phyloai.pretree.align import verify_align_outputs
                ckpt = load_checkpoint(ckpt_path)
                to_run, _ = plan_resume(ckpt, verify_align_outputs)
                total = len(to_run)
            else:
                total = 1
        else:
            from phyloai.pretree.align import _scan_input
            found, _ = _scan_input(seq_dir)
            total = len(found)
        if total == 0:
            click.echo("Run already complete — all tasks verified.", err=True)
            payload = _invoke()
        else:
            with Progress(console=console, transient=True) as progress:
                task = progress.add_task("Aligning sequences", total=total)
                try:
                    payload = _invoke(progress_callback=lambda _: progress.advance(task))
                except (ValueError, FileNotFoundError) as exc:
                    error_msg = str(exc)
    else:
        try:
            payload = _invoke()
        except (ValueError, FileNotFoundError) as exc:
            error_msg = str(exc)

    if error_msg is not None:
        exit_code = 3 if "not found" in error_msg.lower() else 1
        _fail(error_msg, exit_code)

    if dry_run:
        click.echo(f"Dry run: {payload['data']['summary']['n_input_files']} genes would be aligned.")
        if resume:
            from phyloai.core.checkpoint import (
                load_checkpoint,
                summarize_resume_tasks,
                validate_resume_params,
            )
            from phyloai.pretree.align import (
                _detect_seq_type_from_files,
                _resolved_align_params,
                _resolve_tool_paths,
                _scan_input,
                verify_align_outputs,
            )
            from phyloai.pretree.checkpoint_helpers import resume_verifier

            ckpt_path = output_dir / "checkpoint.json"
            found, _ = _scan_input(seq_dir)
            resolved_seq_type = seq_type
            if resolved_seq_type == "auto":
                resolved_seq_type = _detect_seq_type_from_files(found) if found else "AA"
            mafft_exe, magus_exe, trimal_exe = _resolve_tool_paths(
                method=method,
                backtrans=backtrans,
                mafft_path=mafft_path,
                magus_path=magus_path,
                trimal_path=trimal_path,
                dry_run=True,
            )
            resolved = _resolved_align_params(
                seq_dir=seq_dir,
                output_dir=output_dir,
                method=method,
                resolved_seq_type=resolved_seq_type,
                backtrans=backtrans,
                nt_dir=nt_dir,
                threads=threads,
                tool_args=tool_args,
                mafft_path=mafft_exe,
                magus_path=magus_exe,
                trimal_path=trimal_exe,
                quiet=quiet,
            )
            checkpoint = load_checkpoint(ckpt_path)
            validate_resume_params(checkpoint, resolved, step="pretree.align")
            summary = summarize_resume_tasks(checkpoint, resume_verifier(verify_align_outputs))
            click.echo(
                f"Resume dry-run: skip {summary['skip']} tasks, "
                f"rerun {summary['rerun']} tasks, "
                f"invalidate {summary['invalid']} recorded successes.",
                err=True,
            )
        for item in payload["data"].get("files", []):
            if item.get("cmd"):
                click.echo(" ".join(item["cmd"]))
        return

    result_path = output_dir / "result.json"
    result_path.parent.mkdir(parents=True, exist_ok=True)
    with open(result_path, "w") as fh:
        json.dump(payload, fh, indent=2)

    if not quiet:
        console.print(render_align_summary_table(payload["data"]["summary"]))
        click.echo(
            f"Alignments saved to {output_dir / 'seqs'}", err=True
        )
        click.echo(f"Results saved to {result_path}", err=True)
        for w in payload["data"].get("warnings", []):
            click.echo(f"Warning: {w}", err=True)


@pretree.command(
    "trim",
    help=(
        "PhyloAI batch-trims aligned FASTA MSAs using trimAl, BMGE, or ClipKIT.\n\n"
        "Modes are inferred from --seq-type and --nt-dir: AA-only, NT-only, "
        "CODON, or AA+NT dual output. CODON produces seqs/faa and seqs/fna. "
        "AA+NT mode: --msa-dir is aligned AA MSAs. trimAl accepts raw CDS or "
        "gapped codon-aligned NT in --nt-dir because PhyloAI strips NT gaps "
        "before backtranslation. ClipKIT and BMGE expect codon-aligned NT MSAs."
    ),
)
@click.option("--msa-dir", type=click.Path(file_okay=False, path_type=Path), required=True, help="Input directory of aligned MSA files. In AA+NT mode, this is the aligned AA MSA directory; --nt-dir supplies matching NT/CDS files.")
@click.option("--output-dir", "-o", type=click.Path(file_okay=False, path_type=Path), default=Path("runs/pretree/trim"), show_default=True, help="Output directory; contains seqs/, logs/, checkpoint.json, result.json.")
@click.option("--tool", type=click.Choice(["trimal", "bmge", "clipkit"]), default="trimal", show_default=True, help="Trimming tool to use.")
@click.option("--seq-type", type=click.Choice(["AA", "NT", "CODON", "auto"]), default="auto", show_default=True, help="Molecule type. 'auto' detects AA vs NT only; CODON must be explicit.")
@click.option("--nt-dir", type=click.Path(file_okay=False, path_type=Path), default=None, help="NT directory for AA+NT dual output. trimAl accepts raw CDS or gapped codon-aligned NT because PhyloAI strips NT gaps before backtranslation; ClipKIT/BMGE expect codon-aligned NT MSAs.")
@click.option("--trimal-method", type=click.Choice(["automated1", "gappyout", "strict", "strictplus"]), default="automated1", show_default=True, help="trimAl automated trimming strategy.")
@click.option("--bmge-matrix", type=str, default=None, help="BMGE substitution matrix (-m). Defaults: AA/CODON=BLOSUM62; NT=DNAPAM100:2. Common AA options include BLOSUM30, BLOSUM62, BLOSUM90; common NT options include DNAPAM100:2.")
@click.option("--bmge-entropy", type=float, default=0.5, show_default=True, help="BMGE entropy cutoff (-h). Lower values are more stringent.")
@click.option("--clipkit-method", type=click.Choice(["smart-gap", "entropy", "gappy", "block-gappy", "gappyout", "composition-bias", "heterotachy", "kpic", "kpic-smart-gap", "kpic-gappy", "kpi", "kpi-smart-gap", "kpi-gappy", "cst", "c3"]), default="smart-gap", show_default=True, help="ClipKIT trimming mode (-m).")
@click.option("--trimal-path", type=click.Path(dir_okay=False, path_type=Path), default=None, help="Explicit trimAl executable path.")
@click.option("--bmge-path", type=click.Path(dir_okay=False, path_type=Path), default=None, help="Explicit BMGE.jar path.")
@click.option("--clipkit-path", type=click.Path(dir_okay=False, path_type=Path), default=None, help="Explicit clipkit executable path.")
@click.option("--threads", "-t", type=int, default=4, show_default=True, help="Number of genes to trim in parallel.")
@click.option("--tool-args", type=str, default=None, help='Tool strategy arguments only. PhyloAI manages input/output/log/codon/threads, e.g. --tool-args "-g 0.8" for ClipKIT or --tool-args "-m BLOSUM90 -h 0.4" for BMGE.')
@click.option("--resume", is_flag=True, default=False, help="Resume from checkpoint.json in the output directory.")
@click.option("--overwrite", is_flag=True, default=False, help="Delete and recreate a non-empty output directory before running.")
@click.option("--dry-run", is_flag=True, default=False, help="Print commands without executing; creates no files.")
@click.option("--quiet", "-q", is_flag=True, default=False, help="Suppress Rich terminal output except errors.")
def trim_command(
    msa_dir: Path,
    output_dir: Path,
    tool: str,
    seq_type: str,
    nt_dir: Path | None,
    trimal_method: str,
    bmge_matrix: str | None,
    bmge_entropy: float,
    clipkit_method: str,
    trimal_path: Path | None,
    bmge_path: Path | None,
    clipkit_path: Path | None,
    threads: int,
    tool_args: str | None,
    resume: bool,
    overwrite: bool,
    dry_run: bool,
    quiet: bool,
) -> None:
    if threads < 1:
        _fail("--threads must be at least 1.", 1)
    if not msa_dir.exists():
        _fail(f"--msa-dir '{msa_dir}' does not exist.", 1)
    if nt_dir is not None and not nt_dir.exists():
        _fail(f"--nt-dir '{nt_dir}' does not exist.", 1)
    for flag, path in [("--trimal-path", trimal_path), ("--bmge-path", bmge_path), ("--clipkit-path", clipkit_path)]:
        if path is not None and not path.exists():
            _fail(f"{flag} '{path}' does not exist.", 1)

    payload: dict | None = None
    error_msg: str | None = None

    def _invoke(progress_callback=None):
        return run_trim(
            msa_dir=msa_dir,
            output_dir=output_dir,
            tool=tool,
            seq_type=seq_type,
            nt_dir=nt_dir,
            trimal_method=trimal_method,
            bmge_matrix=bmge_matrix,
            bmge_entropy=bmge_entropy,
            clipkit_method=clipkit_method,
            trimal_path=trimal_path,
            bmge_path=bmge_path,
            clipkit_path=clipkit_path,
            threads=threads,
            tool_args=tool_args,
            overwrite=overwrite,
            resume=resume,
            dry_run=dry_run,
            quiet=quiet,
            progress_callback=progress_callback,
        )

    if not quiet and not dry_run:
        if resume and overwrite:
            _fail("--overwrite and --resume are mutually exclusive.", 1)
        if resume:
            ckpt_path = output_dir / "checkpoint.json"
            if ckpt_path.exists():
                from phyloai.core.checkpoint import load_checkpoint
                from phyloai.pretree.checkpoint_helpers import plan_resume
                from phyloai.pretree.trim import verify_trim_outputs
                ckpt = load_checkpoint(ckpt_path)
                to_run, _ = plan_resume(ckpt, verify_trim_outputs)
                total = len(to_run)
            else:
                total = 1
        else:
            found, _ = _trim_scan_input(msa_dir)
            total = len(found)
        if total == 0:
            click.echo("Run already complete — all tasks verified.", err=True)
            payload = _invoke()
        else:
            with Progress(console=console, transient=True) as progress:
                task = progress.add_task("Trimming alignments", total=total)
                try:
                    payload = _invoke(progress_callback=lambda _path: progress.advance(task))
                except (ValueError, FileNotFoundError) as exc:
                    error_msg = str(exc)
    else:
        try:
            payload = _invoke()
        except (ValueError, FileNotFoundError) as exc:
            error_msg = str(exc)

    if error_msg is not None:
        if "No genes were trimmed" in error_msg:
            exit_code = 2
        else:
            exit_code = 3 if "not found" in error_msg.lower() else 1
        _fail(error_msg, exit_code)

    if dry_run:
        click.echo(f"Dry run: {payload['data']['summary']['n_input_files']} genes would be trimmed.")
        for cmd_str in payload["data"].get("dry_run_cmds", []):
            click.echo(cmd_str)
        return

    result_path = output_dir / "result.json"
    result_path.parent.mkdir(parents=True, exist_ok=True)
    with open(result_path, "w") as fh:
        json.dump(payload, fh, indent=2)

    if not quiet:
        console.print(render_trim_summary_table(payload["data"]["summary"]))
        click.echo(f"Trimmed alignments saved to {output_dir / 'seqs'}", err=True)
        click.echo(f"Results saved to {result_path}", err=True)
        for warning in payload["data"].get("warnings", []):
            click.echo(f"Warning: {warning}", err=True)


@pretree.command(
    "concat",
    help=(
        "Concatenate multiple MSA files into a supermatrix for phylogenetic inference. "
        "Supports occupancy filtering, recoding, codon variants, outgroup reordering, "
        "and multi-format output."
    ),
)
@click.option(
    "--msa-dir", type=click.Path(file_okay=False, path_type=Path),
    required=True, help="Directory of input MSA files.",
)
@click.option(
    "--output-dir", "-o", type=click.Path(file_okay=False, path_type=Path),
    default=Path("runs/pretree/concat"), show_default=True, help="Output directory.",
)
@click.option(
    "--prefix", type=str, default="matrix", show_default=True,
    help="Prefix for output filenames.",
)
@click.option(
    "--seq-type", type=click.Choice(["AA", "NT", "CODON", "auto"]),
    default="auto", show_default=True, help="Sequence type.",
)
@click.option(
    "--taxa-occupancy", type=float, default=0.5, show_default=True,
    help="Min taxon ratio for MSA inclusion (0.0-1.0).",
)
@click.option(
    "--recoding", type=click.Choice([
        "RY-nucleotide",
        "Dayhoff-6", "Dayhoff-9", "Dayhoff-12", "Dayhoff-15", "Dayhoff-18",
        "SandR-6", "KGB-6",
    ]), default=None,
    help=(
        "Character recoding scheme. "
        "RY-nucleotide: NT only (A/G->R, C/T/U->Y). "
        "Dayhoff-6/9/12/15/18, SandR-6, KGB-6: AA only."
    ),
)
@click.option(
    "--outgroup", type=str, default=None,
    help="Single taxon name to move to first position in each output matrix.",
)
@click.option(
    "--to", "to",
    type=click.Choice(["fasta", "phylip-relaxed", "phylip-paml", "nexus"]),
    default="fasta", show_default=True, help="Output format.",
)
@click.option(
    "--translate-codon", is_flag=True, default=False,
    help="Also produce CDS-->AA translated matrix (CODON only).",
)
@click.option(
    "--exclude-codon3", is_flag=True, default=False,
    help="Also produce codon1+2 matrix (CODON only).",
)
@click.option(
    "--dry-run", is_flag=True, default=False,
    help="Validate inputs and report planned actions without writing files.",
)
@click.option(
    "--quiet", "-q", is_flag=True, default=False,
    help="Suppress Rich terminal output.",
)
@click.option(
    "--overwrite", is_flag=True, default=False,
    help="Delete and recreate non-empty output directory.",
)
def concat_command(
    msa_dir: Path,
    output_dir: Path,
    prefix: str,
    seq_type: str,
    taxa_occupancy: float,
    recoding: str | None,
    outgroup: str | None,
    to: str,
    translate_codon: bool,
    exclude_codon3: bool,
    dry_run: bool,
    quiet: bool,
    overwrite: bool,
) -> None:
    if not msa_dir.exists():
        _fail(f"MSA directory '{msa_dir}' does not exist.", 1)
    if not (0.0 <= taxa_occupancy <= 1.0):
        _fail("--taxa-occupancy must be between 0.0 and 1.0.", 1)

    payload: dict | None = None
    error_msg: str | None = None
    try:
        payload = run_concat(
            msa_dir=msa_dir,
            output_dir=output_dir,
            prefix=prefix,
            seq_type=seq_type,
            taxa_occupancy=taxa_occupancy,
            recoding=recoding,
            outgroup=outgroup,
            to=to,
            translate_codon=translate_codon,
            exclude_codon3=exclude_codon3,
            dry_run=dry_run,
            overwrite=overwrite,
            quiet=quiet,
        )
    except ValueError as exc:
        error_msg = str(exc)
        if not dry_run and output_dir.exists():
            import json

            err_payload = {
                "status": "error",
                "command": _build_concat_command(msa_dir, output_dir, prefix, seq_type, taxa_occupancy, recoding, outgroup, to, translate_codon, exclude_codon3, dry_run, overwrite, quiet=quiet),
                "wall_time": 0.0,
                "tool_versions": {},
                # NOTE: params must be kept in sync with concat.py run_concat()
                "params": {
                    "msa_dir": str(msa_dir),
                    "output_dir": str(output_dir),
                    "prefix": prefix,
                    "seq_type": seq_type,
                    "taxa_occupancy": taxa_occupancy,
                    "recoding": recoding,
                    "outgroup": outgroup,
                    "to": to,
                    "translate_codon": translate_codon,
                    "exclude_codon3": exclude_codon3,
                    "dry_run": dry_run,
                    "overwrite": overwrite,
                    "quiet": quiet,
                },
                "key_results": {},
                "error": error_msg,
                "data": {"cmd": [], "tool_stderr": ""},
            }
            result_path = output_dir / "result.json"
            with open(result_path, "w") as fh:
                json.dump(err_payload, fh, indent=2)

    if error_msg is not None:
        _fail(error_msg, 1)

    if not quiet and payload is not None:
        overview = {
            "prefix": prefix,
            "to": payload["params"]["to"],
            "n_taxa": payload["key_results"]["n_taxa"],
            "n_msa_input": payload["key_results"]["n_msa_input"],
            "n_msa_used": payload["key_results"]["n_msa_used"],
            "n_msa_dropped": payload["key_results"]["n_msa_dropped"],
            "taxon_occupancy_threshold": payload["params"]["taxa_occupancy"],
            "recoding": payload["params"].get("recoding"),
            "outgroup": payload["params"].get("outgroup"),
            "variants_produced": payload["key_results"]["variants_produced"],
        }
        variant_stats = payload["data"]["variant_stats"]
        panels = _render_concat_panels(overview, variant_stats)
        for panel in panels:
            console.print(panel)

    if payload is not None and not dry_run and not quiet:
        click.echo(f"Results saved to {output_dir / 'result.json'}", err=True)
    elif dry_run and not quiet:
        click.echo("[dry-run] No files written.", err=True)


# ---------------------------------------------------------------------------
# metrics group
# ---------------------------------------------------------------------------


@click.group(
    "metrics",
    invoke_without_command=True,
    help=(
        "Compute MSA and tree metrics, distribution plots, and correlation "
        "analysis for molecular marker evaluation."
    ),
)
@click.option("--msa-dir", type=click.Path(exists=True, file_okay=False, path_type=Path), default=None,
              help="Directory of aligned FASTA files (.fa/.fasta/.fas/.fna/.faa/.aln).")
@click.option("--tree-dir", type=click.Path(exists=True, file_okay=False, path_type=Path), default=None,
              help="Directory of Newick tree files (.tre/.tree/.nwk/.newick/.treefile/.bestTree/.contree).")
@click.option("--seq-type", type=click.Choice(["AA", "NT", "auto"]), default="auto", show_default=True,
              help="Molecule type: AA, NT, or auto-detect per marker.")
@click.option("--outgroup-list", type=click.Path(exists=True, dir_okay=False, path_type=Path), default=None,
              help="File with one outgroup taxon name per line for DVMC pruning.")
@click.option("--ref-tree", type=click.Path(exists=True, dir_okay=False, path_type=Path), default=None,
              help="Reference species tree for normalized Robinson-Foulds distance.")
@click.option("--skip-freq-statistics", is_flag=True, default=False,
              help="Skip per-character frequency columns (freqA, freqC, ...).")
@click.option("--pseudo-tree-metrics", is_flag=True, default=False,
              help="Compute FastTree-derived pseudo-tree metrics (_FT suffix).")
@click.option("--fasttree-path", type=click.Path(dir_okay=False, path_type=Path), default=None,
              help="Explicit path to FastTree executable.")
@click.option("--skip-pairwise-identity", is_flag=True, default=False,
              help="Skip average_pairwise_identity (O(n^2 x L); recommended for >200 taxa).")
@click.option("--round", "decimal_places", type=click.IntRange(0, 12), default=6, show_default=True,
              help="Decimal places for numeric values in the metrics table (0=integer).")
@click.option("--table-format", type=click.Choice(["csv", "tsv"]), default="csv", show_default=True,
              help="Table format for auxiliary tabular outputs (metrics table, basic statistics, correlation matrix).")
@click.option("--output-dir", "-o", type=click.Path(file_okay=False, path_type=Path),
              default=Path("runs/pretree/metrics"), show_default=True,
              help="Output directory for the metrics table, plots/, correlation_heatmap.pdf, result.json.")
@click.option("--threads", "-t", type=int, default=4, show_default=True,
              help="Number of worker processes.")
@click.option("--dry-run", is_flag=True, default=False,
              help="Validate inputs and show plan without writing files.")
@click.option("--overwrite", is_flag=True, default=False,
              help="Delete and recreate a non-empty output directory.")
@click.option("--quiet", "-q", is_flag=True, default=False,
              help="Suppress terminal output.")
@click.pass_context
def metrics_group(
    ctx: click.Context,
    msa_dir: Path | None,
    tree_dir: Path | None,
    seq_type: str,
    outgroup_list: Path | None,
    ref_tree: Path | None,
    skip_freq_statistics: bool,
    pseudo_tree_metrics: bool,
    fasttree_path: Path | None,
    skip_pairwise_identity: bool,
    decimal_places: int,
    table_format: str,
    output_dir: Path,
    threads: int,
    dry_run: bool,
    overwrite: bool,
    quiet: bool,
) -> None:
    if ctx.invoked_subcommand is not None:
        return

    if not msa_dir and not tree_dir:
        _fail("At least one of --msa-dir or --tree-dir must be provided.", 1)
    if pseudo_tree_metrics and not msa_dir:
        _fail("--pseudo-tree-metrics requires --msa-dir.", 1)
    if outgroup_list and not tree_dir:
        _fail("--outgroup-list requires --tree-dir.", 1)
    if ref_tree and not tree_dir:
        _fail("--ref-tree requires --tree-dir.", 1)
    if threads < 1:
        _fail("--threads must be at least 1.", 1)

    if output_dir.exists() and any(output_dir.iterdir()) and not overwrite and not dry_run:
        _fail(
            f"Output directory '{output_dir}' already exists and is non-empty. "
            "Use --overwrite to replace it.",
            1,
        )

    progress = None
    if not quiet:
        progress = Progress(console=console, transient=True)
        progress.start()

    # Step 1: compute metrics
    payload = run_metrics(
        msa_dir=msa_dir, tree_dir=tree_dir, seq_type=seq_type,
        threads=threads, output_dir=output_dir, decimal_places=decimal_places,
        skip_freq_statistics=skip_freq_statistics, pseudo_tree_metrics=pseudo_tree_metrics,
        fasttree_path=str(fasttree_path) if fasttree_path else "FastTree",
        skip_pairwise_identity=skip_pairwise_identity,
        outgroup_list=outgroup_list, ref_tree=ref_tree,
        overwrite=overwrite, dry_run=dry_run, quiet=quiet,
        progress=progress, console=console,
        table_format=table_format,
    )

    if payload["status"] == "error":
        if progress is not None:
            progress.stop()
        _fail(payload.get("error", "Unknown error"), 1)

    if dry_run:
        if progress is not None:
            progress.stop()
        if not quiet:
            click.echo("[dry-run] No files written; no plots or correlation generated.", err=True)
        return

    # Step 2: generate distribution plots
    n_plots = 0
    try:
        import csv as _csv
        metrics_file = output_dir / f"metrics{_table_suffix(table_format)}"
        rows = []
        delimiter = _get_delimiter(table_format)
        with open(metrics_file, newline="") as fh:
            for row in _csv.DictReader(fh, delimiter=delimiter):
                rows.append(row)
        numeric_cols = [k for k in rows[0].keys() if k not in ("loci", "DataType")] if rows else []
        plots_dir = output_dir / "plots"
        n_plots = _generate_all_plots(rows, numeric_cols, plots_dir)
        _generate_basic_statistics(rows, numeric_cols, output_dir / f"metrics.basic_statistics{_table_suffix(table_format)}", table_format=table_format)
    except Exception as exc:
        if not quiet:
            click.echo(f"\n[WARN] Plot generation failed: {exc}", err=True)

    # Step 3: correlation heatmap
    try:
        corr_cols = _select_correlation_columns(rows, list(rows[0].keys()) if rows else [])
        corr_matrix, col_names = _compute_correlation(rows, corr_cols, method="spearman")
        if corr_matrix.size > 0:
            corr_dir = output_dir / "correlate"
            corr_dir.mkdir(parents=True, exist_ok=True)
            _generate_correlation_heatmap(
                corr_matrix, col_names, corr_dir / "correlation_heatmap.pdf",
                annot=False,
            )
            _write_correlation_csv(corr_matrix, col_names, corr_dir / f"correlation_matrix{_table_suffix(table_format)}", table_format=table_format)
            corr_matrix_path = corr_dir / f"correlation_matrix{_table_suffix(table_format)}"
            corr_heatmap_path = corr_dir / "correlation_heatmap.pdf"
            payload["data"]["output_files"]["correlation_matrix"] = {"path": str(corr_matrix_path), "description": "Pairwise Spearman correlation matrix for all phylogenetic metrics"}
            payload["data"]["output_files"]["correlation_heatmap"] = {"path": str(corr_heatmap_path), "description": "Spearman correlation heatmap of all computed phylogenetic informativeness metrics"}
    except Exception as exc:
        if not quiet:
            click.echo(f"\n[WARN] Correlation generation failed: {exc}", err=True)

    basic_stats_path = output_dir / f"metrics.basic_statistics{_table_suffix(table_format)}"
    payload["data"]["output_files"]["basic_statistics"] = {"path": str(basic_stats_path), "description": "Per-metric summary statistics: min, max, mean, median, standard deviation"}
    payload["data"]["output_files"]["plots_dir"] = {"path": str(plots_dir), "description": "Directory containing distribution plots for each computed metric"}
    payload["data"]["output_files"]["n_plots"] = n_plots

    # Rewrite result.json with updated output_files
    with open(output_dir / "result.json", "w") as fh:
        json.dump(payload, fh, indent=2)

    if progress:
        progress.stop()

    if not quiet:
        summary = payload["data"].get("summary", {})
        suffix = _table_suffix(table_format)
        click.echo(f"Metrics table  → {output_dir / f'metrics{suffix}'}", err=True)
        click.echo(f"Plots        → {output_dir / 'plots'} ({n_plots} PDFs)", err=True)
        click.echo(f"Basic stats  → {output_dir / f'metrics.basic_statistics{suffix}'}", err=True)
        click.echo(f"Correlation  → {output_dir / 'correlate' / 'correlation_heatmap.pdf'}", err=True)
        click.echo(f"Results      → {output_dir / 'result.json'}", err=True)
        click.echo(
            f"n_markers={summary.get('n_markers', '?')}, n_success={summary.get('n_success', '?')}, "
            f"n_errors={summary.get('n_errors', '?')}",
            err=True,
        )


@metrics_group.command(
    "plot",
    help=(
        "Re-generate a single metric's distribution density histogram from "
        "an existing metrics.csv.\n\n"
        "Reads the specified --metric column, filters non-numeric values, "
        "optionally applies Tukey's Fences outlier removal (saving filtered "
        "loci names), and writes a density-normalised histogram PDF with "
        "KDE overlay.\n\n"
        "Useful for iterative plot styling without re-running the full "
        "metrics computation."
    ),
)
@click.option("--csv", "csv_path", type=click.Path(exists=True, dir_okay=False, path_type=Path), required=True,
              help="Path to an existing metrics.csv (required).")
@click.option("--input-format", "input_format", type=click.Choice(["csv", "tsv", "auto"]), default="auto", show_default=True,
              help="Input format of the metrics table (csv, tsv, or auto-detect).")
@click.option("--metric", type=str, required=True,
              help="Exact column name to plot, e.g. entropy, rcfv, saturation (required).")
@click.option("--bins", type=int, default=50, show_default=True,
              help="Number of histogram bins (1-500).")
@click.option("--xmin", type=float, default=None, show_default=False,
              help="Force X-axis lower limit (auto-detected if omitted).")
@click.option("--xmax", type=float, default=None, show_default=False,
              help="Force X-axis upper limit (auto-detected if omitted).")
@click.option("--tukey-k", type=float, default=None, show_default=False,
              help="Tukey's Fences multiplier for outlier removal (e.g. 1.5 = standard, 3.0 = conservative). "
                   "Outlier loci names are saved to <output_dir>/<metric>.tukey_filtered.csv.")
@click.option("--title", type=str, default=None, show_default=False,
              help="Plot title (default: 'Distribution of <metric>').")
@click.option("--xlabel", type=str, default=None, show_default=False,
              help="X-axis label (default: metric display name).")
@click.option("--ylabel", type=str, default="Density", show_default=True,
              help="Y-axis label.")
@click.option("--color", type=str, default="#2E86AB", show_default=True,
              help="Bar fill colour (hex or named colour).")
@click.option("--fig-width", type=float, default=10.0, show_default=True,
              help="Figure width in inches.")
@click.option("--fig-height", type=float, default=8.0, show_default=True,
              help="Figure height in inches.")
@click.option("--dpi", type=int, default=150, show_default=True,
              help="Output resolution (72-600).")
@click.option("--font-size", type=int, default=12, show_default=True,
              help="Base font size in points.")
@click.option("--output-dir", "-o", type=click.Path(file_okay=False, path_type=Path),
              default=None, show_default=False,
              help="Directory for the PDF. Default: <csv_parent>/plot_<metric>/")
@click.option("--overwrite", is_flag=True, default=False)
@click.option("--quiet", "-q", is_flag=True, default=False)
def metrics_plot_command(
    csv_path: Path, input_format: str, metric: str, bins: int, xmin: float | None, xmax: float | None,
    tukey_k: float | None, title: str | None, xlabel: str | None, ylabel: str,
    color: str, fig_width: float, fig_height: float, dpi: int, font_size: int,
    output_dir: Path | None, overwrite: bool, quiet: bool,
) -> None:
    start = time.monotonic()
    import numpy as _np
    import csv as _csv_mod

    if output_dir is None:
        output_dir = csv_path.parent / f"plot_{metric}"
    output_dir = output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    if not overwrite and (output_dir / f"{metric}.pdf").exists():
        _fail(
            f"Output '{output_dir / f'{metric}.pdf'}' already exists. Use --overwrite to replace it.",
            1,
        )
    if overwrite and any(output_dir.iterdir()):
        import shutil
        shutil.rmtree(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

    values = []
    loci_labels = []
    delimiter = _detect_input_delimiter(csv_path, input_format)
    with open(csv_path, newline="") as fh:
        for row in _csv_mod.DictReader(fh, delimiter=delimiter):
            v = row.get(metric)
            lbl = row.get("loci", "")
            if v not in (None, "", "NA"):
                try:
                    values.append(float(v))
                    loci_labels.append(lbl)
                except (ValueError, TypeError):
                    continue
    data = _np.array(values, dtype=float)
    clean = data[_np.isfinite(data)]
    if len(clean) == 0:
        _fail(f"No valid numeric values found for metric '{metric}'.", 1)

    out_path = output_dir / f"{metric}.pdf"
    n_filtered, fltr_pairs = _plot_single_metric(
        clean, metric, out_path, bins=bins, xmin=xmin, xmax=xmax,
        tukey_k=tukey_k, raw_labels=loci_labels, title=title,
        xlabel=xlabel, ylabel=ylabel,
        color=color, fig_width=fig_width, fig_height=fig_height, dpi=dpi,
        font_size=font_size,
    )

    # Save filtered loci as CSV when tukey-k is used
    if tukey_k is not None and fltr_pairs:
        fltr_path = output_dir / f"{metric}.tukey_filtered.csv"
        with open(fltr_path, "w", newline="") as fh:
            writer = _csv_mod.writer(fh)
            writer.writerow(["loci", "value"])
            for loci_name, val in fltr_pairs:
                writer.writerow([loci_name, round(val, 6)])
        if not quiet:
            click.echo(f"Filtered loci list  → {fltr_path} ({n_filtered} filtered)", err=True)

    cmd_parts = ["phyloai", "pretree", "metrics", "plot", "--csv", str(csv_path)]
    if input_format != "auto":
        cmd_parts.extend(["--input-format", input_format])
    cmd_parts.extend(["--metric", metric, "--bins", str(bins)])
    if xmin is not None:
        cmd_parts.extend(["--xmin", str(xmin)])
    if xmax is not None:
        cmd_parts.extend(["--xmax", str(xmax)])
    if tukey_k is not None:
        cmd_parts.extend(["--tukey-k", str(tukey_k)])
    if title is not None:
        cmd_parts.extend(["--title", title])
    if xlabel is not None:
        cmd_parts.extend(["--xlabel", xlabel])
    cmd_parts.extend(["--ylabel", ylabel, "--color", color])
    cmd_parts.extend(["--fig-width", str(fig_width), "--fig-height", str(fig_height), "--dpi", str(dpi), "--font-size", str(font_size)])
    cmd_parts.extend(["--output-dir", str(output_dir)])
    if overwrite:
        cmd_parts.append("--overwrite")
    payload = {
        "status": "success",
        "command": " ".join(cmd_parts),
        "wall_time": round(time.monotonic() - start, 3), "tool_versions": {},
        "params": {"csv": str(csv_path), "input_format": input_format, "metric": metric, "bins": bins,
                   "tukey_k": tukey_k, "n_filtered": n_filtered,
                   "fig_width": fig_width, "fig_height": fig_height, "dpi": dpi, "font_size": font_size,
                   "xmin": xmin, "xmax": xmax, "title": title, "xlabel": xlabel, "ylabel": ylabel,
                   "color": color, "output_dir": str(output_dir), "overwrite": overwrite},
        "key_results": {"n_filtered": n_filtered}, "error": None,
        "data": {"cmd": [], "tool_stderr": "", "output_files": {"plot": {"path": str(out_path), "description": f"Histogram of {metric} values"}}},
    }
    with open(output_dir / "result.json", "w") as fh:
        json.dump(payload, fh, indent=2)
    if not quiet:
        click.echo(f"Plot saved to {out_path}", err=True)


@metrics_group.command(
    "correlate",
    help=(
        "Re-generate Spearman/Pearson correlation heatmap from an existing "
        "metrics.csv.\n\n"
        "By default core numeric columns are correlated; freq* and sd_* "
        "columns are omitted for readability. Use --include-freq, --include-sd, "
        "--metrics, or --metrics all to broaden the selection. "
        "Variables are ordered by Ward clustering on magnitude-based "
        "distance (1 - |correlation|), but no dendrogram is drawn. "
        "Writes correlation_heatmap.pdf, "
        "correlation_matrix.csv, and result.json.\n\n"
        "Cells use an R-style circle corrplot layout."
    ),
)
@click.option("--csv", "csv_path", type=click.Path(exists=True, dir_okay=False, path_type=Path), required=True,
              help="Path to an existing metrics.csv (required).")
@click.option("--input-format", "input_format", type=click.Choice(["csv", "tsv", "auto"]), default="auto", show_default=True,
              help="Input format of the metrics table (csv, tsv, or auto-detect).")
@click.option("--metrics", type=str, default=None,
              help="Comma-separated metric columns to correlate. Use 'all' for every numeric column; default uses core metrics only.")
@click.option("--include-freq", is_flag=True, default=False,
              help="Include freq* columns in automatic metric selection.")
@click.option("--include-sd", is_flag=True, default=False,
              help="Include sd_* columns in automatic metric selection.")
@click.option("--method", type=click.Choice(["spearman", "pearson"]), default="spearman", show_default=True,
              help="Spearman (rank-based, non-parametric) or Pearson (linear, z-score normalized).")
@click.option("--triangle", type=click.Choice(["full", "lower", "upper"]), default="full", show_default=True,
              help="Matrix display: full, lower with left/bottom labels, or upper with top/right labels.")
@click.option("--annot/--no-annot", default=False, show_default=True,
              help="Show numeric correlation values inside cells.")
@click.option("--cluster-rectangles", type=int, default=None,
              help="Draw N cluster rectangles on full matrices only; warns and ignores for --triangle lower/upper.")
@click.option("--cmap", type=str, default="RdBu_r", show_default=True,
              help="Matplotlib colormap (e.g. RdBu_r, coolwarm, viridis).")
@click.option("--fmt", type=str, default=".2f", show_default=True,
              help="Numeric format for cell annotations (only when --annot).")
@click.option("--fig-width", type=float, default=12.0, show_default=True,
              help="Figure width in inches.")
@click.option("--fig-height", type=float, default=10.0, show_default=True,
              help="Figure height in inches.")
@click.option("--dpi", type=int, default=150, show_default=True,
              help="Output resolution (72-600).")
@click.option("--font-size", type=int, default=10, show_default=True,
              help="Base font size for axis labels.")
@click.option("--label-angle", type=float, default=45.0, show_default=True,
              help="Rotation angle for x-axis metric labels in degrees.")
@click.option("--title", type=str, default=None,
              help="Plot title (default: no title).")
@click.option("--output-dir", "-o", type=click.Path(file_okay=False, path_type=Path),
              default=Path("runs/pretree/metrics/correlate"), show_default=True,
              help="Directory for heatmap PDF, correlation matrix CSV, and result.json.")
@click.option("--overwrite", is_flag=True, default=False)
@click.option("--quiet", "-q", is_flag=True, default=False)
def metrics_correlate_command(
    csv_path: Path, input_format: str, metrics: str | None, include_freq: bool, include_sd: bool, method: str, triangle: str,
    annot: bool, cluster_rectangles: int | None, cmap: str, fmt: str,
    fig_width: float, fig_height: float, dpi: int, font_size: int, label_angle: float,
    title: str | None, output_dir: Path, overwrite: bool, quiet: bool,
) -> None:
    start = time.monotonic()
    output_dir = output_dir.resolve()
    import csv as _csv_mod

    output_dir.mkdir(parents=True, exist_ok=True)

    heatmap_path = output_dir / "correlation_heatmap.pdf"
    if not overwrite and heatmap_path.exists():
        _fail(
            f"Output '{heatmap_path}' already exists. Use --overwrite to replace it.",
            1,
        )
    if overwrite and any(output_dir.iterdir()):
        import shutil
        shutil.rmtree(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

    rows = []
    delimiter = _detect_input_delimiter(csv_path, input_format)
    with open(csv_path, newline="") as fh:
        for row in _csv_mod.DictReader(fh, delimiter=delimiter):
            rows.append(row)
    columns = _select_correlation_columns(
        rows,
        list(rows[0].keys()) if rows else [],
        requested=metrics,
        include_freq=include_freq,
        include_sd=include_sd,
    )
    if not columns:
        _fail("No metric columns found for correlation.", 1)

    corr_matrix, col_names = _compute_correlation(rows, columns, method=method)
    if corr_matrix.size == 0:
        _fail("Not enough valid data for correlation analysis.", 1)

    _generate_correlation_heatmap(
        corr_matrix, col_names, heatmap_path,
        triangle=triangle, cluster_rectangles=cluster_rectangles,
        cmap=cmap, annot=annot, fmt=fmt,
        fig_width=fig_width, fig_height=fig_height, dpi=dpi,
        font_size=font_size, title=title,
        label_angle=label_angle,
        warn=(lambda message: click.echo(f"[WARN] {message}", err=True)) if not quiet else None,
    )
    _write_correlation_csv(corr_matrix, col_names, output_dir / "correlation_matrix.csv")

    cmd_parts = ["phyloai", "pretree", "metrics", "correlate", "--csv", str(csv_path)]
    if input_format != "auto":
        cmd_parts.extend(["--input-format", input_format])
    if metrics:
        cmd_parts.extend(["--metrics", metrics])
    if include_freq:
        cmd_parts.append("--include-freq")
    if include_sd:
        cmd_parts.append("--include-sd")
    cmd_parts.extend(["--method", method, "--triangle", triangle])
    if annot:
        cmd_parts.append("--annot")
    if cluster_rectangles is not None:
        cmd_parts.extend(["--cluster-rectangles", str(cluster_rectangles)])
    cmd_parts.extend(["--cmap", cmap, "--fmt", fmt])
    cmd_parts.extend(["--fig-width", str(fig_width), "--fig-height", str(fig_height), "--dpi", str(dpi), "--font-size", str(font_size), "--label-angle", str(label_angle)])
    if title:
        cmd_parts.extend(["--title", title])
    cmd_parts.extend(["--output-dir", str(output_dir)])
    if overwrite:
        cmd_parts.append("--overwrite")
    payload = {
        "status": "success",
        "command": " ".join(cmd_parts),
        "wall_time": round(time.monotonic() - start, 3), "tool_versions": {},
        "params": {"csv": str(csv_path), "input_format": input_format, "metrics": metrics,
                   "include_freq": include_freq, "include_sd": include_sd, "method": method,
                   "triangle": triangle,
                   "annot": annot, "cmap": cmap, "fmt": fmt,
                   "fig_width": fig_width, "fig_height": fig_height, "dpi": dpi, "font_size": font_size,
                   "label_angle": label_angle, "title": title, "output_dir": str(output_dir),
                   "overwrite": overwrite, "cluster_rectangles": cluster_rectangles},
        "key_results": {"n_variables": len(col_names) if col_names else 0}, "error": None,
        "data": {"cmd": [], "tool_stderr": "", "output_files": {
            "correlation_heatmap": {"path": str(heatmap_path), "description": "Spearman correlation heatmap of selected phylogenetic metrics"},
            "correlation_matrix": {"path": str(output_dir / "correlation_matrix.csv"), "description": "Pairwise Spearman correlation matrix in tabular format"},
        }},
    }
    with open(output_dir / "result.json", "w") as fh:
        json.dump(payload, fh, indent=2)
    if not quiet:
        click.echo(f"Heatmap saved to {heatmap_path}", err=True)
        click.echo(f"Correlation matrix saved to {output_dir / 'correlation_matrix.csv'}", err=True)


# ---------------------------------------------------------------------------
# filter group
# ---------------------------------------------------------------------------

from phyloai.pretree.filter import render_filter_summary_table, run_taper, run_treeshrink, run_metrics_filter, run_symtest, run_cluster_filter  # noqa: E402


class _FilterGroup(click.Group):
    def list_commands(self, ctx: click.Context) -> list[str]:
        return ["taper", "treeshrink", "metrics", "symtest", "cluster"]


@click.group(
    "filter",
    cls=_FilterGroup,
    help="TAPER site masking, TreeShrink taxa pruning, "
    "metric-rule loci filtering, symmetry test filtering, "
    "cluster-based exploration.",
)
def filter_group() -> None:
    pass


# ---- filter taper ----

_TAPER_HELP = (
    "Mask erroneous amino-acid or nucleotide sites within multiple sequence "
    "alignments using the TAPER error-correction tool (bundled "
    "correction_multi.jl, executed by Julia).\n"
    "\n"
    "Operating modes (one per paragraph):\n"
    "\n"
    "  AA-only\n"
    "    --msa-dir with AA alignments.\n"
    "    Output: masked AA to seqs/\n"
    "\n"
    "  NT-only\n"
    "    --msa-dir with NT alignments, --seq-type NT.\n"
    "    Output: masked NT to seqs/\n"
    "\n"
    "  AA+CDS\n"
    "    --msa-dir with AA alignments, --nt-dir with codon-aligned NT MSAs.\n"
    "    Output: masked AA to seqs/faa/, projected CDS to seqs/fna/\n"
    "\n"
    "TAPER is run per locus in parallel.  Only newly introduced 'X' masks "
    "(not original ambiguity) are counted.  --resume skips loci whose output "
    "files already exist and pass validation."
)

@filter_group.command("taper", help=_TAPER_HELP)
@click.option(
    "--msa-dir", type=click.Path(exists=True, file_okay=False, path_type=Path),
    required=True,
    help="Directory containing input MSA files (any suffix).  All regular non-empty "
    "files are scanned; format is validated when parsed.",
)
@click.option(
    "--nt-dir", type=click.Path(exists=True, file_okay=False, path_type=Path),
    default=None,
    help="AA+CDS mode only: directory of codon-aligned nucleotide MSAs "
    "(one per AA locus, length == 3 * AA length).  Requires amino-acid --msa-dir input.",
)
@click.option(
    "--seq-type", type=click.Choice(["AA", "NT", "auto"]),
    default="auto", show_default=True,
    help="Expected molecule type of the MSA files.  'auto' detects from the first "
    "file's sequence characters (presence of EFILPQWYZ -> AA, otherwise NT).  "
    "AA+CDS mode (--nt-dir) requires AA input; --seq-type NT is rejected with --nt-dir.",
)
@click.option(
    "--cutoff", type=click.IntRange(1), default=3, show_default=True,
    help="TAPER -c error-correction cutoff.  Lower values mask more aggressively "
    "(1 = most aggressive); values 1-10 are typical.  The TAPER default is 3.",
)
@click.option(
    "--taper-path", type=click.Path(exists=True, dir_okay=False, path_type=Path),
    default=None,
    help="Explicit path to correction_multi.jl.  When omitted, the bundled copy "
    "inside the PhyloAI package is used.",
)
@click.option(
    "--julia-path", type=click.Path(exists=True, dir_okay=False, path_type=Path),
    default=None,
    help="Explicit path to the Julia executable.  When omitted, PhyloAI resolves "
    "it via PATH (see 'phyloai doctor' for detection status).",
)
@click.option(
    "--tool-args", type=str, default=None,
    help="Additional TAPER command-line flags passed through verbatim.  "
    "PhyloAI-managed flags (-m, -a, -c, -l, input path, output redirection) "
    "must not appear here; the command exits with an error if they do.",
)
@click.option(
    "--threads", "-t", type=int, default=4, show_default=True,
    help="Number of parallel worker processes (one locus per worker).",
)
@click.option(
    "--output-dir", "-o", type=click.Path(file_okay=False, path_type=Path),
    default=Path("runs/pretree/filter/taper"), show_default=True,
    help="Directory for masked output files, decision tables, result.json, and logs.",
)
@click.option(
    "--table-format", type=click.Choice(["csv", "tsv"]),
    default="csv", show_default=True,
    help="Delimiter and file suffix for auxiliary tables (retained_loci, dropped_loci, "
    "filter_decisions).  Does not affect result.json.",
)
@click.option(
    "--show-masked-sites", is_flag=True, default=False,
    help="Include per-taxon masked-site counts in filter_decisions.csv.  "
    "Default off to keep the output table compact.",
)
@click.option(
    "--resume", is_flag=True, default=False,
    help="Resume a previous run from the checkpoint.json inside --output-dir.  "
    "Parameters must match exactly; successfully completed loci are skipped.",
)
@click.option(
    "--overwrite", is_flag=True, default=False,
    help="Delete and recreate --output-dir if it already exists.  "
    "Mutually exclusive with --resume.",
)
@click.option(
    "--dry-run", is_flag=True, default=False,
    help="Validate inputs and show the planned TAPER commands and output layout "
    "without executing anything or writing files.",
)
@click.option(
    "--quiet", "-q", is_flag=True, default=False,
    help="Suppress all terminal output except errors.",
)
def filter_taper_command(msa_dir, nt_dir, seq_type, cutoff, taper_path, julia_path, tool_args, threads, output_dir, table_format, show_masked_sites, resume, overwrite, dry_run, quiet):
    if threads < 1:
        _fail("--threads must be at least 1.", 1)
    if nt_dir is not None and seq_type == "NT":
        _fail("--nt-dir (AA+CDS mode) requires amino-acid input.  --seq-type NT is incompatible with --nt-dir.", 1)

    def _invoke(progress_callback=None):
        return run_taper(msa_dir=msa_dir, output_dir=output_dir, seq_type=seq_type, nt_dir=nt_dir, cutoff=cutoff, taper_path=taper_path, julia_path=julia_path, threads=threads, tool_args=tool_args, resume=resume, overwrite=overwrite, dry_run=dry_run, quiet=quiet, table_format=table_format, show_masked_sites=show_masked_sites, progress_callback=progress_callback)

    error_msg = None
    if not quiet and not dry_run:
        if resume and overwrite:
            _fail("--overwrite and --resume are mutually exclusive.", 1)
        from phyloai.core.file_matching import scan_msa_dir

        if resume:
            ckpt_path = output_dir / "checkpoint.json"
            if ckpt_path.exists():
                from phyloai.core.checkpoint import load_checkpoint as _load_ckpt
                from phyloai.pretree.checkpoint_helpers import plan_resume
                from phyloai.pretree.filter import _verify_taper_outputs
                ckpt = _load_ckpt(ckpt_path)
                to_run, _ = plan_resume(ckpt, _verify_taper_outputs)
                total = len(to_run)
                label = f"TAPER masking (resume, {len(ckpt.tasks)} total)"
            else:
                total = 1
                label = "TAPER masking"
        else:
            msa_map = scan_msa_dir(msa_dir)
            total = max(len(msa_map), 1)
            label = "TAPER masking"
        if total == 0:
            click.echo("Run already complete — all tasks verified.", err=True)
            payload = _invoke()
        else:
            with Progress(console=console, transient=True) as progress:
                task = progress.add_task(label, total=total)
                try:
                    payload = _invoke(progress_callback=lambda _: progress.advance(task))
                except (ValueError, FileNotFoundError) as exc:
                    error_msg = str(exc)
    else:
        try:
            payload = _invoke()
        except (ValueError, FileNotFoundError) as exc:
            error_msg = str(exc)

    if error_msg is not None:
        exit_code = 3 if "not found" in error_msg.lower() else 1
        _fail(error_msg, exit_code)
    if dry_run:
        click.echo(f"Dry run: {payload['key_results']['n_input']} loci would be processed.")
        for cmd in payload["data"]["dry_run_cmds"]:
            click.echo(cmd)
        return
    if not quiet:
        console.print(render_filter_summary_table({
            "Input": payload["key_results"]["n_input"],
            "Retained": payload["key_results"]["n_retained"],
            "Dropped": payload["key_results"]["n_dropped"],
            "Masked loci": payload["key_results"]["masked_loci"],
            "Masked taxa": payload["key_results"]["total_masked_taxa"],
            "Masked sites": payload["key_results"]["total_masked_aa_sites"],
        }))
        summary_data = payload["data"].get("summary", {})
        msa_stats = summary_data.get("retained_msa_stats", {})
        if msa_stats and msa_stats.get("n_msa", 0) > 0:
            console.print(render_filter_summary_table({
                "Retained MSAs": msa_stats["n_msa"],
                "Total length": msa_stats["total_length"],
                "Mean length": msa_stats["mean_length"],
                "Min length": msa_stats["min_length"],
                "Max length": msa_stats["max_length"],
                "Mean taxa": msa_stats["mean_taxa"],
            }))
        click.echo(f"Masked MSAs saved to {output_dir / 'seqs'}", err=True)
        click.echo(f"Results saved to {output_dir / 'result.json'}", err=True)
    if payload["status"] == "error":
        _fail(payload.get("error", "All loci failed."), 1)


# ---- filter treeshrink ----

_TREESHRINK_HELP = (
    "Detect and prune outlier long-branch taxa from gene trees using TreeShrink.\n\n"
    "TreeShrink is run once across the entire gene-tree dataset (not per gene) "
    "because it can use information from multiple trees jointly.  PhyloAI creates "
    "a per-gene working layout (input.tree, optional input.fasta) in a temporary "
    "directory, invokes run_treeshrink.py, then collects the shrunk outputs.\n\n"
    "When --msa-dir is provided, matching MSAs are also shrunk to remove the "
    "same pruned taxa, and retained-MSA statistics are included in the output."
)

@filter_group.command("treeshrink", help=_TREESHRINK_HELP)
@click.option(
    "--tree-dir", type=click.Path(exists=True, file_okay=False, path_type=Path),
    required=True,
    help="Directory of input gene tree files (any suffix).  All regular non-empty "
    "files are scanned; logical locus names are derived from filename stems.",
)
@click.option(
    "--msa-dir", type=click.Path(exists=True, file_okay=False, path_type=Path),
    default=None,
    help="Optional directory of MSA files paired with the gene trees by logical "
    "locus name.  When provided, TreeShrink is invoked with both tree and "
    "alignment input per locus, and shrunk MSAs are written to seqs/.",
)
@click.option(
    "--threshold", type=click.FloatRange(0.0), default=0.05, show_default=True,
    help="TreeShrink -q false-positive threshold.  Smaller values remove more taxa; "
    "0.05 is the TreeShrink default.  Must be >=0.",
)
@click.option(
    "--treeshrink-mode", type=click.Choice(["auto", "per-gene", "all-genes", "per-species"]),
    default="auto", show_default=True,
    help="TreeShrink -m operating mode.  'auto' omits -m (TreeShrink default).  "
    "'per-gene' runs independently per gene.  'all-genes' and 'per-species' "
    "use cross-gene information.",
)
@click.option(
    "--treeshrink-path", type=click.Path(exists=True, dir_okay=False, path_type=Path),
    default=None,
    help="Explicit path to run_treeshrink.py.  When omitted, PhyloAI resolves it "
    "via PATH (see 'phyloai doctor').",
)
@click.option(
    "--tool-args", type=str, default=None,
    help="Additional TreeShrink flags passed through verbatim.  "
    "PhyloAI-managed flags (-i, -t, -a, -q, -m, -o, -O) must not appear here.",
)
@click.option(
    "--keep-work-dir", is_flag=True, default=False,
    help="Retain the per-gene working directory under --output-dir/work/ for "
    "debugging.  By default the temporary work layout is deleted after the run.",
)
@click.option(
    "--output-dir", "-o", type=click.Path(file_okay=False, path_type=Path),
    default=Path("runs/pretree/filter/treeshrink"), show_default=True,
    help="Directory for shrunk trees, optional shrunk MSAs, decision tables, "
    "result.json, and logs.",
)
@click.option(
    "--table-format", type=click.Choice(["csv", "tsv"]),
    default="csv", show_default=True,
    help="Delimiter and file suffix for auxiliary tables (retained_loci, "
    "modified_loci, dropped_loci, removed_taxa, filter_decisions).  "
    "Does not affect result.json.",
)
@click.option(
    "--overwrite", is_flag=True, default=False,
    help="Delete and recreate --output-dir if it already exists.",
)
@click.option(
    "--dry-run", is_flag=True, default=False,
    help="Validate inputs and show the resolved TreeShrink command and number of "
    "loci without executing or writing files.",
)
@click.option(
    "--quiet", "-q", is_flag=True, default=False,
    help="Suppress all terminal output except errors.",
)
def filter_treeshrink_command(tree_dir, msa_dir, threshold, treeshrink_mode, treeshrink_path, tool_args, keep_work_dir, output_dir, table_format, overwrite, dry_run, quiet):
    if not quiet and not dry_run:
        with Progress(console=console, transient=True) as progress:
            progress.add_task("TreeShrink running...", total=None)
            try:
                payload = run_treeshrink(tree_dir=tree_dir, output_dir=output_dir, msa_dir=msa_dir, threshold=threshold, treeshrink_mode=treeshrink_mode, treeshrink_path=treeshrink_path, tool_args=tool_args, keep_work_dir=keep_work_dir, overwrite=overwrite, dry_run=dry_run, quiet=quiet, table_format=table_format)
            except (ValueError, FileNotFoundError) as exc:
                _fail(str(exc), 3 if "not found" in str(exc).lower() else 1)
    else:
        try:
            payload = run_treeshrink(tree_dir=tree_dir, output_dir=output_dir, msa_dir=msa_dir, threshold=threshold, treeshrink_mode=treeshrink_mode, treeshrink_path=treeshrink_path, tool_args=tool_args, keep_work_dir=keep_work_dir, overwrite=overwrite, dry_run=dry_run, quiet=quiet, table_format=table_format)
        except (ValueError, FileNotFoundError) as exc:
            _fail(str(exc), 3 if "not found" in str(exc).lower() else 1)
    if dry_run:
        click.echo(f"Dry run: would process {payload['key_results']['n_input']} loci.")
        click.echo(payload["data"]["dry_run_cmd"])
        return
    if not quiet:
        console.print(render_filter_summary_table({"Input": payload["key_results"]["n_input"], "Retained": payload["key_results"]["n_retained"], "Modified": payload["key_results"]["n_modified"], "Dropped": payload["key_results"]["n_dropped"], "Taxa removed": payload["key_results"]["n_removed_taxa_total"]}))
        summary_data = payload["data"].get("summary", {})
        msa_stats = summary_data.get("retained_msa_stats", {})
        if msa_stats and msa_stats.get("n_msa", 0) > 0:
            console.print(render_filter_summary_table({
                "Retained MSAs": msa_stats["n_msa"],
                "Total length": msa_stats["total_length"],
                "Mean length": msa_stats["mean_length"],
                "Min length": msa_stats["min_length"],
                "Max length": msa_stats["max_length"],
                "Mean taxa": msa_stats["mean_taxa"],
            }))
        console.print("[dim]Tip: filtered alignments may be used to re-construct "
                      "phylogenetic trees, which are possibly more accurate than "
                      "those pruned by TreeShrink.[/dim]")
        click.echo(f"Shrunk trees saved to {output_dir / 'trees'}", err=True)
        click.echo(f"Results saved to {output_dir / 'result.json'}", err=True)
    if payload["status"] == "error":
        _fail(payload.get("error", "All loci failed."), 1)


# ---- filter metrics ----

_METRICS_HELP = (
    "Filter whole loci by explicit numeric or string conditions on a metrics "
    "CSV/TSV table (typically the output of 'phyloai pretree metrics').\n\n"
    "All conditions in --keep are combined with AND logic (a locus must satisfy "
    "every condition to be retained).  OR logic is not supported in this version.\n\n"
    "Use --copy to copy retained MSA and/or tree files into the output directory.  "
    "Without --copy only decision tables are written, which is useful for "
    "threshold exploration without duplicating large files."
)

@filter_group.command("metrics", help=_METRICS_HELP)
@click.option(
    "--table", "table_path", type=click.Path(exists=True, dir_okay=False, path_type=Path),
    required=True,
    help="Path to a CSV or TSV metrics table.  Delimiter is auto-detected from "
    "file content unless --input-format csv|tsv is specified explicitly.",
)
@click.option(
    "--keep", type=str, required=True,
    help="Comma-separated list of AND conditions.  Supported operators: "
    ">=, >, <=, <, ==, !=.  Numeric comparisons require numeric values; "
    "only == and != are valid for string columns.  "
    "Example: 'dvmc>=0,dvmc<=0.3,average_BS>=0.8,DataType==AA'.",
)
@click.option(
    "--input-format", type=click.Choice(["csv", "tsv", "auto"]),
    default="auto", show_default=True,
    help="Delimiter format of the input --table.  'auto' inspects file content; "
    "use csv or tsv to override when auto-detection is ambiguous.",
)
@click.option(
    "--loci-column", type=str, default="loci", show_default=True,
    help="Name of the column that holds the logical locus identifier.",
)
@click.option(
    "--msa-dir", type=click.Path(exists=True, file_okay=False, path_type=Path),
    default=None,
    help="Directory of MSA files for computing retained-MSA statistics (terminal "
    "summary and result.json).  When combined with --copy, retained MSAs are "
    "copied to --output-dir/seqs/.",
)
@click.option(
    "--tree-dir", type=click.Path(exists=True, file_okay=False, path_type=Path),
    default=None,
    help="Directory of tree files.  When combined with --copy, retained trees are "
    "copied to --output-dir/trees/.",
)
@click.option(
    "--copy", is_flag=True, default=False,
    help="Copy retained MSA and/or tree files into the output directory.  "
    "Requires at least one of --msa-dir or --tree-dir.",
)
@click.option(
    "--output-dir", "-o", type=click.Path(file_okay=False, path_type=Path),
    default=Path("runs/pretree/filter/metrics"), show_default=True,
    help="Directory for decision tables, optional copied files, result.json, and logs.",
)
@click.option(
    "--table-format", type=click.Choice(["csv", "tsv"]),
    default="csv", show_default=True,
    help="Delimiter and file suffix for auxiliary tables (retained_loci, "
    "dropped_loci, filter_decisions).  Does not affect result.json.",
)
@click.option(
    "--overwrite", is_flag=True, default=False,
    help="Delete and recreate --output-dir if it already exists.",
)
@click.option(
    "--dry-run", is_flag=True, default=False,
    help="Parse --keep rules and show how many loci would be retained/dropped "
    "without writing any files.",
)
@click.option(
    "--quiet", "-q", is_flag=True, default=False,
    help="Suppress all terminal output except errors.",
)
def filter_metrics_command(table_path, keep, input_format, loci_column, msa_dir, tree_dir, copy, output_dir, table_format, overwrite, dry_run, quiet):
    try:
        payload = run_metrics_filter(table_path=table_path, output_dir=output_dir, keep=keep, input_format=input_format, loci_column=loci_column, msa_dir=msa_dir, tree_dir=tree_dir, copy=copy, overwrite=overwrite, dry_run=dry_run, quiet=quiet, table_format=table_format)
    except (ValueError, FileNotFoundError) as exc:
        _fail(str(exc), 1)
    if dry_run:
        click.echo(f"Dry run: {payload['key_results']['n_total']} loci -> {payload['key_results']['n_retained']} retained, {payload['key_results']['n_dropped']} dropped")
        return
    if not quiet:
        console.print(render_filter_summary_table({"Total": payload["key_results"]["n_total"], "Retained": payload["key_results"]["n_retained"], "Dropped": payload["key_results"]["n_dropped"]}))
        fm_summary = payload["data"].get("summary", {})
        msa_stats = fm_summary.get("retained_msa_stats", {})
        if msa_stats and msa_stats.get("n_msa", 0) > 0:
            console.print(render_filter_summary_table({
                "Retained MSAs": msa_stats["n_msa"],
                "Total length": msa_stats["total_length"],
                "Mean length": msa_stats["mean_length"],
                "Min length": msa_stats["min_length"],
                "Max length": msa_stats["max_length"],
                "Mean taxa": msa_stats["mean_taxa"],
            }))
        click.echo(f"Decision tables saved to {output_dir}", err=True)
        click.echo(f"Results saved to {output_dir / 'result.json'}", err=True)
    if payload["status"] == "error":
        _fail(payload.get("error", "Filtering failed."), 1)


# ---- filter symtest ----

_SYMTEST_HELP = (
    "Test phylogenetic symmetry assumptions per locus using IQ-TREE3's "
    "--symtest-only, then filter loci by p-value.\n\n"
    "The p-value column used depends on --symtest-type:\n"
    "  (default)  SymPval  combined stationarity + homogeneity\n"
    "  MAR        MarPval  marginal / stationarity\n"
    "  INT        IntPval  internal / homogeneity\n\n"
    "References: Naser-Khdour et al. (2019) doi:10.1093/gbe/evz193"
)


def _validate_symtest_pval(ctx, param, value):
    if value <= 0 or value > 1:
        raise click.BadParameter("must be > 0 and <= 1")
    return value


@filter_group.command("symtest", help=_SYMTEST_HELP)
@click.option(
    "--msa-dir", type=click.Path(exists=True, file_okay=False, path_type=Path),
    required=True,
    help="Directory containing per-locus MSA files (any suffix).",
)
@click.option(
    "--symtest-type", type=click.Choice(["MAR", "INT"]),
    default=None,
    help="Which symmetry test to use for filtering.  When omitted (default), "
    "the combined Sym test is used (SymPval column).  MAR uses marginal "
    "(stationarity) test.  INT uses internal (homogeneity) test.",
)
@click.option(
    "--symtest-pval", type=float, default=0.05, show_default=True,
    callback=_validate_symtest_pval,
    help="P-value threshold.  Loci with p >= threshold are retained.",
)
@click.option(
    "--symtest-keep-zero", is_flag=True, default=False,
    help="Pass --symtest-keep-zero to IQ-TREE (keep NAs in the tests).",
)
@click.option(
    "--iqtree-path", type=click.Path(exists=True, dir_okay=False, path_type=Path),
    default=None,
    help="Explicit path to iqtree binary.  When omitted, resolved via "
    "PATH ('phyloai doctor' for detection status).",
)
@click.option(
    "--threads", "-t", type=int, default=4, show_default=True,
    help="Number of threads for IQ-TREE (-T).",
)
@click.option(
    "--tree-dir", type=click.Path(exists=True, file_okay=False, path_type=Path),
    default=None,
    help="Optional directory of gene tree files.  Trees matching retained "
    "loci (by logical locus name) are copied to trees/.",
)
@click.option(
    "--output-dir", "-o", type=click.Path(file_okay=False, path_type=Path),
    default=Path("runs/pretree/filter/symtest"), show_default=True,
    help="Directory for retained MSAs, optional trees, decision tables, "
    "result.json, and logs.",
)
@click.option(
    "--table-format", type=click.Choice(["csv", "tsv"]),
    default="csv", show_default=True,
    help="Delimiter and file suffix for auxiliary tables.",
)
@click.option(
    "--overwrite", is_flag=True, default=False,
    help="Delete and recreate --output-dir if it already exists.",
)
@click.option(
    "--dry-run", is_flag=True, default=False,
    help="Validate inputs and show the planned IQ-TREE command without writing files.",
)
@click.option(
    "--quiet", "-q", is_flag=True, default=False,
    help="Suppress all terminal output except errors.",
)
def filter_symtest_command(msa_dir, symtest_type, symtest_pval, symtest_keep_zero,
                           iqtree_path, threads, tree_dir, output_dir, table_format,
                           overwrite, dry_run, quiet):
    if threads < 1:
        _fail("--threads must be at least 1.", 1)

    if not quiet and not dry_run:
        from phyloai.core.file_matching import scan_msa_dir
        msa_map = scan_msa_dir(msa_dir)
        with Progress(console=console, transient=True) as progress:
            progress.add_task("IQ-TREE symmetry test running...", total=None)
            try:
                payload = run_symtest(
                    msa_dir=msa_dir, output_dir=output_dir,
                    symtest_type=symtest_type, symtest_pval=symtest_pval,
                    symtest_keep_zero=symtest_keep_zero,
                    iqtree_path=iqtree_path, threads=threads,
                    tree_dir=tree_dir, msa_map=msa_map,
                    table_format=table_format,
                    dry_run=dry_run, overwrite=overwrite, quiet=quiet,
                )
            except (ValueError, FileNotFoundError, RuntimeError) as exc:
                msg = str(exc)
                exit_code = 3 if "not found" in msg.lower() else (
                    2 if "exited with code" in msg.lower() else 1)
                _fail(msg, exit_code)
    else:
        try:
            payload = run_symtest(
                msa_dir=msa_dir, output_dir=output_dir,
                symtest_type=symtest_type, symtest_pval=symtest_pval,
                symtest_keep_zero=symtest_keep_zero,
                iqtree_path=iqtree_path, threads=threads,
                tree_dir=tree_dir, table_format=table_format,
                dry_run=dry_run, overwrite=overwrite, quiet=quiet,
            )
        except (ValueError, FileNotFoundError, RuntimeError) as exc:
            msg = str(exc)
            exit_code = 3 if "not found" in msg.lower() else (
                2 if "exited with code" in msg.lower() else 1)
            _fail(msg, exit_code)

    if dry_run:
        click.echo(f"Dry run: {payload['key_results']['n_input']} loci would be processed.")
        click.echo(payload["data"]["dry_run_cmd"])
        return

    if not quiet:
        console.print(render_filter_summary_table({
            "Input": payload["key_results"]["n_input"],
            "Retained": payload["key_results"]["n_retained"],
            "Dropped": payload["key_results"]["n_dropped"],
            "P-value threshold": payload["key_results"]["p_value_threshold"],
            "Symtest type": payload["key_results"]["symtest_type"],
        }))
        summary_data = payload["data"].get("summary", {})
        msa_stats = summary_data.get("retained_msa_stats", {})
        if msa_stats and msa_stats.get("n_msa", 0) > 0:
            console.print(render_filter_summary_table({
                "Retained MSAs": msa_stats["n_msa"],
                "Total length": msa_stats["total_length"],
                "Mean length": msa_stats["mean_length"],
                "Min length": msa_stats["min_length"],
                "Max length": msa_stats["max_length"],
                "Mean taxa": msa_stats["mean_taxa"],
            }))
        if payload["key_results"].get("retained_trees_copied", 0) > 0:
            mt = summary_data.get("missed_tree_count", 0)
            console.print(render_filter_summary_table({
                "Trees copied": payload["key_results"]["retained_trees_copied"],
                "Trees missed": mt,
            }))
        click.echo(f"Retained MSAs saved to {output_dir / 'seqs'}", err=True)
        if payload["key_results"].get("retained_trees_copied", 0) > 0:
            click.echo(f"Retained trees saved to {output_dir / 'trees'}", err=True)
        click.echo(f"Results saved to {output_dir / 'result.json'}", err=True)


# ---- filter cluster ----

_CLUSTER_HELP = (
    "Group loci by their metric profiles using dimensionality reduction "
    "(PCA or UMAP) followed by hierarchical clustering.\n\n"
    "This is primarily an exploratory tool: by default it only writes clusters, "
    "diagnostic plots, and per-cluster metric summaries.  Use "
    "--drop-outlier-clusters auto to optionally remove the worst-performing "
    "clusters based on a specified metric (e.g. average_BS).\n\n"
    "Clustering is run once in-memory; --resume and --threads are not supported "
    "in this version.  PCA is the default reduction method because it is "
    "deterministic and dependency-light; UMAP is recommended for exploring "
    "non-linear structure but adds the optional umap-learn dependency."
)

@filter_group.command("cluster", help=_CLUSTER_HELP)
@click.option(
    "--table", "table_path", type=click.Path(exists=True, dir_okay=False, path_type=Path),
    required=True,
    help="Path to a CSV or TSV metrics table.  Delimiter auto-detected unless "
    "--input-format is specified.",
)
@click.option(
    "--input-format", type=click.Choice(["csv", "tsv", "auto"]),
    default="auto", show_default=True,
    help="Delimiter format of the input --table.",
)
@click.option(
    "--metrics", type=str, default="all", show_default=True,
    help="Comma-separated list of metric column names to use as features.  "
    "Default 'all' uses all numeric columns except the locus identifier, "
    "DataType column, constant columns, and columns matching --exclude-regex.",
)
@click.option(
    "--exclude-regex", type=str, multiple=True, default=None,
    help="Regular expression for excluding metric columns.  Repeatable: "
    "--exclude-regex '^freq' --exclude-regex '^sd_'.  Matched columns are "
    "dropped from the feature set (logged in features_used output table).",
)
@click.option(
    "--reduction", type=click.Choice(["pca", "umap"]),
    default="pca", show_default=True,
    help="Dimensionality reduction method.  'pca' uses sklearn PCA (n_components=3).  "
    "'umap' requires 'pip install umap-learn'; supports --umap-replicates for "
    "stability assessment.",
)
@click.option(
    "--n-clusters", type=int, default=None,
    help="Fixed number of clusters.  When omitted, the best k in [2..max_clusters] "
    "is selected by multi-metric voting (silhouette, Calinski-Harabasz, "
    "Davies-Bouldin).",
)
@click.option(
    "--max-clusters", type=int, default=None,
    help="Upper bound for automatic cluster-count search.  Defaults to "
    "min(30, max(6, ceil(sqrt(n_loci)/3))).  Ignored when --n-clusters is set.",
)
@click.option(
    "--cluster-linkage", type=click.Choice(["ward", "average", "complete", "single"]),
    default="ward", show_default=True,
    help="Agglomerative clustering linkage criterion.  'ward' minimizes "
    "within-cluster variance (requires Euclidean distance).  'average' and "
    "'complete' are conservative alternatives.  'single' is prone to chaining.",
)
@click.option(
    "--cluster-distance", type=click.Choice(["euclidean", "cosine", "manhattan"]),
    default="euclidean", show_default=True,
    help="Distance metric for clustering.  'ward' linkage requires 'euclidean'.",
)
@click.option(
    "--drop-outlier-clusters", type=click.Choice(["none", "auto"]),
    default="none", show_default=True,
    help="When 'auto', rank clusters by mean --outlier-metric and remove the "
    "worst-performing clusters up to --max-drop-fraction of total loci.  "
    "When 'none', no loci are removed; only diagnostics are written.",
)
@click.option(
    "--outlier-metric", type=str, default="average_BS", show_default=True,
    help="Metric column used to rank clusters for outlier removal.  "
    "Clusters with low (or high, per --outlier-direction) mean values are "
    "dropped first.  Ignored unless --drop-outlier-clusters auto.",
)
@click.option(
    "--outlier-direction", type=click.Choice(["low", "high"]),
    default="low", show_default=True,
    help="Interpretation of --outlier-metric for ranking.  'low' means smaller "
    "values are worse (e.g. average_BS).  'high' means larger values are worse.",
)
@click.option(
    "--max-drop-fraction", type=click.FloatRange(0.0, 1.0),
    default=0.2, show_default=True,
    help="Maximum fraction of total loci that can be removed by outlier-cluster "
    "dropping.  Range [0.0, 1.0].",
)
@click.option(
    "--plot-metrics-cols", type=int, default=2, show_default=True,
    help="Number of metric boxplot columns per figure (for cluster_metric_boxplots).  "
    "All metrics are placed in a single figure with this many columns; rows are auto-calculated.",
)
@click.option(
    "--plot-label-angle", type=float, default=45.0, show_default=True,
    help="Rotation angle in degrees for x-axis labels in diagnostic plots.",
)
@click.option(
    "--outlier-boxplot-cols", type=int, default=4, show_default=True,
    help="Number of boxplot columns per figure (for outlier_comparison_boxplots).  "
    "All metrics are placed in a single figure with this many columns; rows are auto-calculated.",
)
@click.option(
    "--umap-n-neighbors", type=int, default=15, show_default=True,
    help="UMAP n_neighbors parameter.  Controls the balance between local and "
    "global structure.  Ignored for PCA reduction.",
)
@click.option(
    "--umap-min-dist", type=float, default=0.001, show_default=True,
    help="UMAP min_dist parameter.  Controls how tightly points are packed.  "
    "Ignored for PCA reduction.",
)
@click.option(
    "--umap-replicates", type=int, default=1, show_default=True,
    help="Number of UMAP replicates with different random seeds.  The best "
    "replicate is selected by cluster-validation rank-sum scoring.  "
    "Ignored for PCA reduction.",
)
@click.option(
    "--umap-random-state", type=int, default=42, show_default=True,
    help="Random seed for reproducible UMAP.  Only applied when --umap-replicates 1; "
    "when replicates > 1, no seed is set so --threads can parallelize across CPUs.  "
    "Ignored for PCA reduction.",
)
@click.option(
    "--threads", type=int, default=1, show_default=True,
    help="Number of CPU threads for UMAP (n_jobs).  Only takes effect when "
    "--reduction umap and --umap-replicates > 1.  Ignored for PCA reduction.",
)
@click.option(
    "--msa-dir", type=click.Path(exists=True, file_okay=False, path_type=Path),
    default=None,
    help="Directory of MSA files.  Used with --copy to copy retained MSAs into "
    "--output-dir/seqs/ when outlier dropping is active.",
)
@click.option(
    "--tree-dir", type=click.Path(exists=True, file_okay=False, path_type=Path),
    default=None,
    help="Directory of tree files.  Used with --copy to copy retained trees into "
    "--output-dir/trees/ when outlier dropping is active.",
)
@click.option(
    "--copy", is_flag=True, default=False,
    help="Copy retained MSA and/or tree files into the output directory.  "
    "Only takes effect when --drop-outlier-clusters auto actually drops loci.",
)
@click.option(
    "--output-dir", "-o", type=click.Path(file_okay=False, path_type=Path),
    default=Path("runs/pretree/filter/cluster"), show_default=True,
    help="Directory for cluster assignments, diagnostic plots, decision tables, "
    "result.json, and logs.",
)
@click.option(
    "--table-format", type=click.Choice(["csv", "tsv"]),
    default="csv", show_default=True,
    help="Delimiter and file suffix for all auxiliary tables (features_used, "
    "reduction, clusters, cluster_summary, cluster_metric_means, and "
    "optionally retained/dropped/filter_decisions).  Does not affect result.json.",
)
@click.option(
    "--overwrite", is_flag=True, default=False,
    help="Delete and recreate --output-dir if it already exists.",
)
@click.option(
    "--dry-run", is_flag=True, default=False,
    help="Validate inputs, resolve the feature set, and show planned reduction, "
    "cluster-count range, and whether outlier dropping would be applied.",
)
@click.option(
    "--quiet", "-q", is_flag=True, default=False,
    help="Suppress all terminal output except errors.",
)
def filter_cluster_command(table_path, input_format, metrics, exclude_regex, reduction, n_clusters, max_clusters, cluster_linkage, cluster_distance, drop_outlier_clusters, outlier_metric, outlier_direction, max_drop_fraction, plot_metrics_cols, plot_label_angle, outlier_boxplot_cols, umap_n_neighbors, umap_min_dist, umap_replicates, umap_random_state, threads, msa_dir, tree_dir, copy, output_dir, table_format, overwrite, dry_run, quiet):
    try:
        payload = run_cluster_filter(table_path=table_path, output_dir=output_dir, input_format=input_format, metrics=metrics, exclude_regex=list(exclude_regex) if exclude_regex else None, reduction=reduction, n_clusters=n_clusters, max_clusters=max_clusters, cluster_linkage=cluster_linkage, cluster_distance=cluster_distance, drop_outlier_clusters=drop_outlier_clusters, outlier_metric=outlier_metric, outlier_direction=outlier_direction, max_drop_fraction=max_drop_fraction, plot_metrics_cols=plot_metrics_cols, plot_label_angle=plot_label_angle, outlier_boxplot_cols=outlier_boxplot_cols, umap_n_neighbors=umap_n_neighbors, umap_min_dist=umap_min_dist, umap_replicates=umap_replicates, umap_random_state=umap_random_state, threads=threads, msa_dir=msa_dir, tree_dir=tree_dir, copy=copy, overwrite=overwrite, dry_run=dry_run, quiet=quiet, table_format=table_format)
    except (ValueError, FileNotFoundError, ImportError) as exc:
        _fail(str(exc), 1)
    if dry_run:
        click.echo(f"Dry run: {payload['key_results']['n_loci']} loci, {payload['key_results']['n_features']} features")
        return
    if not quiet:
        console.print(render_filter_summary_table({
            "Loci": payload["key_results"]["n_loci"],
            "Valid loci": payload["key_results"]["n_valid_loci"],
            "Features": payload["key_results"]["n_features"],
            "Reduction": payload["key_results"]["reduction"],
            "Clusters": payload["key_results"]["n_clusters"],
            "Dropped": payload["key_results"]["n_dropped"],
        }))
        fc_summary = payload["data"].get("summary", {})
        if fc_summary.get("drop_clusters"):
            drop_list = fc_summary["drop_clusters"]
            console.print(f"[yellow]Dropped clusters: {drop_list} "
                          f"({payload['key_results']['n_dropped']} loci removed)[/yellow]")
        msa_stats = fc_summary.get("retained_msa_stats", {})
        if msa_stats and msa_stats.get("n_msa", 0) > 0:
            console.print(render_filter_summary_table({
                "Retained MSAs": msa_stats["n_msa"],
                "Total length": msa_stats["total_length"],
                "Mean length": msa_stats["mean_length"],
                "Min length": msa_stats["min_length"],
                "Max length": msa_stats["max_length"],
                "Mean taxa": msa_stats["mean_taxa"],
            }))
        console.print(f"Results saved to {output_dir}")


pretree.add_command(filter_group)

# Register the group on the pretree instance
pretree.add_command(metrics_group)
