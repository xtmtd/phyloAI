"""phyloai report — generate reproducible analysis reports."""

from __future__ import annotations

import json
from pathlib import Path

import click
from rich.console import Console

from phyloai.report.collector import discover_steps
from phyloai.report.renderer import render_html
from phyloai.report.schema import assemble_report
from phyloai.report.templates import generate_all_methods

console = Console()


def _fail(message: str, exit_code: int = 1) -> None:
    click.echo(f"Error: {message}", err=True)
    raise click.exceptions.Exit(exit_code)


@click.command()
@click.option(
    "--run-dir",
    required=True,
    type=click.Path(exists=True, file_okay=False, path_type=Path),
    help="Run directory to report on (pipeline or module output).",
)
@click.option(
    "-o", "--output-dir",
    type=click.Path(path_type=Path),
    help="Output directory for report files. Default: <run-dir>/report",
)
@click.option(
    "--overwrite",
    is_flag=True,
    help="Overwrite existing report files.",
)
@click.option(
    "-q", "--quiet",
    is_flag=True,
    help="Suppress terminal output except errors.",
)
def report(
    run_dir: Path,
    output_dir: Path | None,
    overwrite: bool,
    quiet: bool,
) -> None:
    """Generate a reproducible analysis report from a PhyloAI run directory.

    Produces report.json (machine-readable, AI/MCP diagnostic entry point)
    and report.html (human-readable, with embedded figures and methods draft).

    \b
    Examples:
      phyloai report --run-dir ./runs/run/faa
      phyloai report --run-dir ./runs/pretree -o ./my-report
    """
    run_dir = run_dir.resolve()

    if output_dir is None:
        output_dir = run_dir / "report"
    output_dir = output_dir.resolve()

    report_json_path = output_dir / "report.json"
    report_html_path = output_dir / "report.html"

    if not overwrite and (report_json_path.exists() or report_html_path.exists()):
        _fail(
            f"Report files already exist in {output_dir}. "
            f"Use --overwrite to replace them."
        )

    output_dir.mkdir(parents=True, exist_ok=True)

    if not quiet:
        console.print(f"[bold]Scanning[/bold] {run_dir}")
    try:
        discovered = discover_steps(run_dir)
    except ValueError as e:
        _fail(str(e))
    except Exception as e:
        _fail(f"Failed to scan run directory: {e}")

    if not quiet:
        console.print(
            f"  Run mode: [bold]{discovered['run_mode']}[/bold] "
            f"({len(discovered['steps'])} steps found)"
        )

    for raw_step in discovered["steps"]:
        step_id = raw_step["step_id"]
        status = raw_step.get("status", "error")

        # Enrich key_results by merging values that some commands put
        # outside key_results: data.summary (convert dir, stats dir) or
        # flat data.* (stats single-file).
        key_results = dict(raw_step.get("key_results", {}))
        raw_data = raw_step.get("data", {})
        data_summary = raw_data.get("summary", {})
        for k, v in data_summary.items():
            if isinstance(v, (int, float, str)) and k not in key_results:
                key_results[k] = v
            elif isinstance(v, dict) and all(isinstance(x, (int, float)) for x in v.values()):
                for sk, sv in v.items():
                    fk = f"{k}_{sk}"
                    if fk not in key_results:
                        key_results[fk] = sv
            elif isinstance(v, list) and k not in key_results:
                key_results[k] = v
        # Fallback: stats single-file puts scalars directly in data
        _STRUCTURAL_KEYS = {
            "output_files", "files", "per_gene", "cmd", "tool_stderr",
            "tool_log", "summary", "variant_stats", "dropped_alignments",
            "per_taxon", "per_gene_occupancy", "skipped", "warnings",
            "character_summary", "site_patterns", "recoding_warnings",
            "normalization_replacements",
        }
        for k, v in raw_data.items():
            if k not in _STRUCTURAL_KEYS and isinstance(v, (int, float, str, bool)) and k not in key_results:
                key_results[k] = v
        # Merge concat-specific metrics from variant_stats[0] (original variant)
        if "gap_ratio" not in key_results or "pi_ratio" not in key_results:
            variants = raw_data.get("variant_stats", [])
            if variants:
                orig = variants[0]
                cs = orig.get("character_summary", {})
                sp = orig.get("site_patterns", {})
                if "gap_ratio" not in key_results and "gap_ratio" in cs:
                    key_results["gap_ratio"] = cs["gap_ratio"]
                if "pi_ratio" not in key_results:
                    pi = sp.get("parsimony_informative")
                    if isinstance(pi, dict):
                        key_results["pi_ratio"] = pi.get("ratio", 0)
        # Flatten nested scalar dicts already in key_results (e.g. trim's
        # length_before: {mean, min, max}).
        for k in list(key_results.keys()):
            v = key_results[k]
            if isinstance(v, dict) and all(isinstance(x, (int, float)) for x in v.values()):
                for sk, sv in v.items():
                    fk = f"{k}_{sk}"
                    if fk not in key_results:
                        key_results[fk] = sv

        text = generate_all_methods(
            step_id,
            params=raw_step.get("params", {}),
            key_results=key_results,
            tool_versions=raw_step.get("tool_versions", {}),
            status=status,
        )
        raw_step["methods_text"] = text

    if not quiet:
        console.print("[bold]Assembling[/bold] report.json")

    report_dict = assemble_report(discovered, run_dir)

    with open(report_json_path, "w") as fh:
        json.dump(report_dict, fh, indent=2, ensure_ascii=False)

    if not quiet:
        n_ok = sum(1 for s in report_dict["steps"] if s["status"] == "success")
        n_fail = sum(1 for s in report_dict["steps"] if s["status"] == "error")
        status_color = "green" if n_fail == 0 else "yellow"
        console.print(
            f"  Status: [{status_color}]{report_dict['status']}[/{status_color}] "
            f"({n_ok} success, {n_fail} failed)"
        )

    if not quiet:
        console.print("[bold]Rendering[/bold] report.html")

    report_dict["run_dir"] = str(run_dir)
    html_path = render_html(report_dict, output_dir)
    if not quiet:
        console.print(f"  [green]report.html[/green] → {html_path}")

    if not quiet:
        console.print("\n[bold green]Report generated:[/bold green]")
        console.print(f"  {report_json_path}")
        console.print(f"  {report_html_path}")
