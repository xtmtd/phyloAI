"""Pre-tree CLI commands."""

from __future__ import annotations

import json
from pathlib import Path

import click
from rich.console import Console
from rich.progress import Progress

from phyloai.pretree.align import render_align_summary_table, run_align
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
from phyloai.pretree.trim import render_trim_summary_table, run_trim, _scan_input as _trim_scan_input
from phyloai.pretree.concat import run_concat, _render_concat_panels

console = Console()


def _fail(message: str, exit_code: int) -> None:
    click.echo(f"Error: {message}", err=True)
    raise click.exceptions.Exit(exit_code)


class _PretreeGroup(click.Group):
    def list_commands(self, ctx: click.Context) -> list[str]:
        return ["convert", "stats", "align", "trim", "concat"]


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
@click.option("--output-dir", "output_dir", type=click.Path(file_okay=False, path_type=Path), default=Path("runs/pretree/convert"), show_default=True, help="Directory where converted files and result.json are written.")
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
              help="Output directory; contains seqs/, align.log, result.json.")
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
        from phyloai.pretree.align import _scan_input
        found, _ = _scan_input(seq_dir)
        with Progress(console=console, transient=True) as progress:
            task = progress.add_task("Aligning sequences", total=len(found))
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
                mafft_executable=mafft_exe,
                magus_executable=magus_exe,
                trimal_executable=trimal_exe,
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
@click.option("--output-dir", "-o", type=click.Path(file_okay=False, path_type=Path), default=Path("runs/pretree/trim"), show_default=True, help="Output directory; contains seqs/, trim.log, checkpoint.json, result.json.")
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
            clipkit_mode=clipkit_method,
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
        found, _ = _trim_scan_input(msa_dir)
        with Progress(console=console, transient=True) as progress:
            task = progress.add_task("Trimming alignments", total=len(found))
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
    "--recoding", type=str, default=None,
    help="Recoding scheme: RY-nucleotide, Dayhoff-6/9/12/15/18, SandR-6, KGB-6.",
)
@click.option(
    "--outgroup", type=str, default=None,
    help="Taxon name to move to first position.",
)
@click.option(
    "--to", "to_format",
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
    to_format: str,
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
            to_format=to_format,
            translate_codon=translate_codon,
            exclude_codon3=exclude_codon3,
            dry_run=dry_run,
            overwrite=overwrite,
        )
    except ValueError as exc:
        error_msg = str(exc)

    if error_msg is not None:
        _fail(error_msg, 1)

    if not quiet and payload is not None:
        display_stats = {
            "prefix": prefix,
            "seq_type": payload["params"]["seq_type"],
            "to_format": payload["params"]["to_format"],
            "n_taxa": payload["key_results"]["n_taxa"],
            "n_msa_input": payload["key_results"]["n_msa_input"],
            "n_msa_used": payload["key_results"]["n_msa_used"],
            "n_msa_dropped": payload["key_results"]["n_msa_dropped"],
            "total_length": payload["key_results"]["total_length"],
            "taxon_occupancy_threshold": payload["params"]["taxa_occupancy"],
            "recoding": payload["params"].get("recoding"),
            "outgroup": payload["params"].get("outgroup"),
            "variants_produced": payload["key_results"]["variants_produced"],
            "character_summary": payload["data"]["character_summary"],
            "site_patterns": payload["data"]["site_patterns"],
        }
        panels = _render_concat_panels(display_stats)
        for panel in panels:
            console.print(panel)

    if payload is not None and not dry_run:
        click.echo(f"Results saved to {output_dir / 'result.json'}", err=True)
    elif dry_run:
        click.echo("[dry-run] No files written.", err=True)
