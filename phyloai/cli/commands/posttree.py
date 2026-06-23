"""Post-tree analysis CLI commands."""
from __future__ import annotations

import json
import shlex
import shutil
from pathlib import Path

import click
from rich.console import Console
from rich.table import Table

from phyloai.core.iqtree import IQTREE_COMPATIBLE_EXTENSIONS

console = Console()


class _PosttreeGroup(click.Group):
    def list_commands(self, ctx: click.Context) -> list[str]:
        return ["topology"]


@click.group(cls=_PosttreeGroup)
def posttree() -> None:
    """Post-tree analysis commands."""


def _fail(message: str, exit_code: int = 1) -> None:
    click.echo(f"Error: {message}", err=True)
    raise click.exceptions.Exit(exit_code)


def _write_error_result_json(
    output_dir: Path,
    command: str,
    error_msg: str,
    error_category: str = "input",
) -> None:
    """Write a spec-compliant error result.json before exiting on validation failure."""
    result_path = output_dir / "result.json"
    result_path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "status": "error",
        "command": command,
        "wall_time": 0.0,
        "tool_versions": {},
        "params": {},
        "key_results": {},
        "error": error_msg,
        "error_category": error_category,
        "data": {"cmd": [], "tool_stderr": "", "tests": [], "warnings": [error_msg]},
    }
    with open(result_path, "w") as fh:
        json.dump(payload, fh, indent=2)


def _display_tests_table(tests: list[dict], quiet: bool = False) -> None:
    """Render the USER TREES results as a Rich table."""
    if not tests or quiet:
        return

    # Collect all available columns across rows
    score_cols = ["bp_rell", "p_kh", "p_sh", "p_wkh", "p_wsh", "c_elw", "p_au"]
    available = []
    for col in score_cols:
        if any(t.get(col) is not None for t in tests):
            available.append(col)

    table = Table(title="Topology Test Results")
    table.add_column("Tree", justify="right", style="cyan")
    table.add_column("logL", justify="right")
    table.add_column("deltaL", justify="right")

    col_styles = {
        "bp_rell": ("bp-RELL", "magenta"),
        "p_kh": ("p-KH", "green"),
        "p_sh": ("p-SH", "green"),
        "p_wkh": ("p-WKH", "green"),
        "p_wsh": ("p-WSH", "green"),
        "c_elw": ("c-ELW", "yellow"),
        "p_au": ("p-AU", "red"),
    }
    for col in available:
        label, style = col_styles[col]
        table.add_column(label, justify="right", style=style)

    for t in tests:
        tid = t.get("tree_id")
        row = [
            str(tid) if tid is not None else "?",
            f"{t.get('log_likelihood', '-'):.3f}" if t.get("log_likelihood") is not None else "-",
            f"{t.get('delta_likelihood', '-'):.3f}" if t.get("delta_likelihood") is not None else "-",
        ]
        for col in available:
            val = t.get(col)
            sign = t.get(col + "_sign", "")
            if val is not None and sign:
                row.append(f"{val:.4f} {sign}")
            elif val is not None:
                row.append(f"{val:.4f}")
            else:
                row.append("-")
        table.add_row(*row)

    console.print()
    console.print(table)


@posttree.command("topology")
@click.option(
    "--matrix", type=click.Path(path_type=Path), default=None,
    help="Single supermatrix alignment (FASTA/PHYLIP/NEXUS/CLUSTAL).  Maps to IQ-TREE -s.",
)
@click.option(
    "--candidate-trees", "candidate_trees_raw",
    type=str,
    help=(
        "Candidate tree input. One tree-list file (one NEWICK tree per line), "
        "or multiple individual NEWICK tree files separated by commas "
        "(e.g. h1.nwk,h2.nwk,h3.nwk). Multiple files are merged in order by PhyloAI."
    ),
)
@click.option(
    "--input-format",
    type=click.Choice(["auto", "fasta", "phylip-relaxed", "nexus"]),
    default="auto", show_default=True,
    help="PhyloAI-side matrix format hint for preflight validation. Not passed to IQ-TREE.",
)
@click.option(
    "--model-expr", type=str, default=None,
    help="Complete IQ-TREE -m model expression (e.g. LG+F+R4, C20+F+R4).",
)
@click.option(
    "--partitions", type=click.Path(path_type=Path), default=None,
    help="Previously optimized partition model. Maps to IQ-TREE -p.",
)
@click.option(
    "--guide-tree", type=click.Path(path_type=Path), default=None,
    help="Guide tree for PMSF-style model expressions. Maps to IQ-TREE -ft.",
)
@click.option(
    "--replicates", type=int, default=10000,
    help="RELL replicates (min 1000, default 10000). Maps to IQ-TREE -zb.",
)
@click.option(
    "--prefix", type=str, default=None,
    help="IQ-TREE output prefix (default: matrix file stem).",
)
@click.option(
    "-o", "--output-dir", type=click.Path(path_type=Path),
    default=Path("runs/posttree/topology"),
    help="Output directory.",
)
@click.option(
    "-t", "--threads", type=int, default=4,
    help="Thread count. Maps to IQ-TREE -T unless overridden by --tool-args.",
)
@click.option(
    "--iqtree-path", type=str, default=None,
    help="Explicit path to iqtree3 executable.",
)
@click.option(
    "--tool-args", type=str, default=None,
    help="Additional IQ-TREE strategy parameters. Blocked flags: -s, -z.",
)
@click.option("--overwrite", is_flag=True, default=False,
              help="Delete and recreate output directory.")
@click.option("--resume", is_flag=True, default=False,
              help="Reuse existing output directory with IQ-TREE native resume.")
@click.option("--dry-run", is_flag=True, default=False,
              help="Print the IQ-TREE command without executing.")
@click.option("-q", "--quiet", is_flag=True, default=False,
              help="Suppress terminal output except errors.")
def topology_command(
    matrix: Path | None,
    candidate_trees_raw: str | None,
    input_format: str,
    model_expr: str | None,
    partitions: Path | None,
    guide_tree: Path | None,
    replicates: int,
    prefix: str | None,
    output_dir: Path,
    threads: int,
    iqtree_path: str | None,
    tool_args: str | None,
    overwrite: bool,
    resume: bool,
    dry_run: bool,
    quiet: bool,
) -> None:
    """Tree topology tests (AU / KH / SH / WKH / WSH / c-ELW).

    Tests whether alternative topologies are significantly worse than the
    best-scoring candidate. Does NOT infer new trees — use
    `phyloai tree ml iqtree` for ML tree inference and model selection.

    PhyloAI does NOT expose ModelFinder here because topology tests are run
    after model inference — you should already have a preferred model.

    Default tests: -n 0 -zb <replicates> -zw -au
    (bp-RELL, KH, SH, weighted KH, weighted SH, c-ELW, AU).

    KH / SH / WKH / WSH / AU are p-values (< 0.05 = rejected).
    bp-RELL and c-ELW are weights (not p-values). Recommended: AU, WSH, WKH.

    Examples:

      # Homogeneous model

      phyloai posttree topology --matrix raw.fa --candidate-trees trees --model-expr LG+F+R4

      # Heterogeneous model

      phyloai posttree topology --matrix raw.fa --candidate-trees trees --model-expr C20+F+R4

      # PMSF model with guide tree

      phyloai posttree topology --matrix raw.fa --candidate-trees trees --model-expr LG+C20+F+R4 --guide-tree guide.tree

      # Partition model

      phyloai posttree topology --matrix raw.fa --candidate-trees trees --partitions raw.best_model.nex

      # Multiple individual tree files (comma-separated, merged by PhyloAI)

      phyloai posttree topology --matrix raw.fa --candidate-trees h1.nwk,h2.nwk,h3.nwk --model-expr LG+F+R4

      # Custom exchangeabilities + site frequencies via --tool-args

      phyloai posttree topology --matrix raw.fa --candidate-trees trees --model-expr custom.exchangeabilities+R4 --tool-args "-fs custom.sitefreq"
    """
    from phyloai.posttree.topology import run_topology

    try:
        _run_topology_impl(
            matrix=matrix,
            candidate_trees_raw=candidate_trees_raw,
            input_format=input_format,
            model_expr=model_expr,
            partitions=partitions,
            guide_tree=guide_tree,
            replicates=replicates,
            prefix=prefix,
            output_dir=output_dir,
            threads=threads,
            iqtree_path=iqtree_path,
            tool_args=tool_args,
            overwrite=overwrite,
            resume=resume,
            dry_run=dry_run,
            quiet=quiet,
        )
    except click.exceptions.Exit as e:
        # Build a best-effort command string for the error result.json
        err_cmd_parts = ["phyloai", "posttree", "topology"]
        if matrix:
            err_cmd_parts.extend(["--matrix", str(matrix)])
        if candidate_trees_raw:
            err_cmd_parts.extend(["--candidate-trees", candidate_trees_raw])
        err_cmd_parts.extend(["--input-format", input_format])
        if model_expr:
            err_cmd_parts.extend(["--model-expr", model_expr])
        if partitions:
            err_cmd_parts.extend(["--partitions", str(partitions)])
        if guide_tree:
            err_cmd_parts.extend(["--guide-tree", str(guide_tree)])
        err_cmd_parts.extend(["--replicates", str(replicates)])
        err_cmd_parts.extend(["-o", str(output_dir)])
        err_cmd_parts.extend(["-t", str(threads)])
        if iqtree_path:
            err_cmd_parts.extend(["--iqtree-path", iqtree_path])
        if tool_args:
            err_cmd_parts.extend(["--tool-args", tool_args])
        if overwrite:
            err_cmd_parts.append("--overwrite")
        if resume:
            err_cmd_parts.append("--resume")
        if dry_run:
            err_cmd_parts.append("--dry-run")
        if quiet:
            err_cmd_parts.append("-q")
        err_cmd = shlex.join(err_cmd_parts)
        _write_error_result_json(output_dir.resolve(), err_cmd, str(e), "input")
        raise


def _run_topology_impl(
    *,
    matrix: Path | None,
    candidate_trees_raw: str | None,
    input_format: str,
    model_expr: str | None,
    partitions: Path | None,
    guide_tree: Path | None,
    replicates: int,
    prefix: str | None,
    output_dir: Path,
    threads: int,
    iqtree_path: str | None,
    tool_args: str | None,
    overwrite: bool,
    resume: bool,
    dry_run: bool,
    quiet: bool,
) -> None:
    from phyloai.posttree.topology import run_topology

    if matrix is None:
        _fail("--matrix is required", exit_code=1)
    if not candidate_trees_raw:
        _fail("At least one --candidate-trees is required", exit_code=1)

    matrix_path: Path = matrix

    # Matrix extension
    ext = matrix_path.suffix.lower()
    if ext not in IQTREE_COMPATIBLE_EXTENSIONS:
        _fail(
            f"Unsupported matrix extension: {ext}. "
            f"Supported: {', '.join(sorted(IQTREE_COMPATIBLE_EXTENSIONS))}",
            exit_code=1,
        )

    # Candidate trees: parse comma-separated individual tree files, or a single
    # tree-list file.  Strip whitespace from each segment.
    raw_parts = [p.strip() for p in candidate_trees_raw.split(",")]
    if not raw_parts or all(p == "" for p in raw_parts):
        _fail("--candidate-trees must contain at least one non-empty path", exit_code=1)
    candidate_trees_list = [Path(p) for p in raw_parts]
    for i, ct in enumerate(candidate_trees_list):
        if not ct.is_file():
            _fail(f"--candidate-trees #{i + 1} is not a regular file: {ct}", exit_code=1)
        if ct.stat().st_size == 0:
            _fail(f"--candidate-trees #{i + 1} is empty: {ct}", exit_code=1)

    # Model source
    has_explicit = model_expr is not None or partitions is not None
    has_tool_args_model = False
    if tool_args:
        tokens = shlex.split(tool_args)
        has_tool_args_model = "-m" in tokens or "-p" in tokens
    if not has_explicit and not has_tool_args_model:
        _fail(
            "Neither --model-expr, --partitions, nor -m/-p in --tool-args provided. "
            "Must specify one model source.",
            exit_code=1,
        )
    if model_expr and partitions:
        _fail("--model-expr and --partitions are mutually exclusive.", exit_code=1)

    # Cross-source model conflict
    if tool_args:
        tokens = shlex.split(tool_args)
        if model_expr and "-p" in tokens:
            _fail(
                "--model-expr is set but --tool-args contains -p. "
                "Remove --model-expr if you want -p from --tool-args to take effect.",
                exit_code=1,
            )
        if partitions and "-m" in tokens:
            _fail(
                "--partitions is set but --tool-args contains -m. "
                "Remove --partitions if you want -m from --tool-args to take effect.",
                exit_code=1,
            )

    # partitions / guide-tree existence and readability
    import os as _os
    if partitions:
        if not partitions.exists():
            _fail(f"--partitions does not exist: {partitions}", exit_code=1)
        if not partitions.is_file():
            _fail(f"--partitions is not a regular file: {partitions}", exit_code=1)
        if not _os.access(str(partitions), _os.R_OK):
            _fail(f"--partitions is not readable: {partitions}", exit_code=1)
    if guide_tree:
        if not guide_tree.exists():
            _fail(f"--guide-tree does not exist: {guide_tree}", exit_code=1)
        if not guide_tree.is_file():
            _fail(f"--guide-tree is not a regular file: {guide_tree}", exit_code=1)
        if not _os.access(str(guide_tree), _os.R_OK):
            _fail(f"--guide-tree is not readable: {guide_tree}", exit_code=1)

    # overwrite / resume
    if overwrite and resume:
        _fail("--overwrite and --resume are mutually exclusive.", exit_code=1)

    # Numeric bounds
    if replicates < 1000:
        _fail(f"--replicates must be >= 1000, got {replicates}", exit_code=1)
    if threads < 1:
        _fail(f"--threads must be >= 1, got {threads}", exit_code=1)

    # --tool-args blocked flags
    if tool_args:
        from phyloai.posttree.topology import _check_managed_flag_conflict
        try:
            _check_managed_flag_conflict(
                tool_args, blocked_flags=frozenset({"-s", "-z"}),
            )
        except ValueError as e:
            _fail(str(e), exit_code=1)

    # ---- Output directory lifecycle (CLI layer) ----

    if not dry_run:
        output_dir = output_dir.resolve()
        if not overwrite and not resume:
            if output_dir.exists() and any(output_dir.iterdir()):
                _fail(
                    f"Output directory exists and is not empty: {output_dir}\n"
                    "Use --overwrite to replace or --resume to reuse.",
                    exit_code=1,
                )
        if overwrite and output_dir.exists():
            shutil.rmtree(output_dir)

    # ---- Execute ----

    guide_tree_str = str(guide_tree) if guide_tree else None

    payload = run_topology(
        matrix=matrix_path,
        candidate_trees=candidate_trees_list,
        input_format=input_format,
        model_expr=model_expr,
        partitions=str(partitions) if partitions else None,
        guide_tree=guide_tree_str,
        replicates=replicates,
        prefix=prefix,
        output_dir=output_dir,
        threads=threads,
        iqtree_path=iqtree_path,
        tool_args=tool_args,
        overwrite=overwrite,
        resume=resume,
        dry_run=dry_run,
        quiet=quiet,
        stream_output=True,
    )

    # ---- Write / display ----

    if not dry_run:
        result_path = output_dir / "result.json"
        result_path.parent.mkdir(parents=True, exist_ok=True)
        with open(result_path, "w") as fh:
            json.dump(payload, fh, indent=2)
        if not quiet:
            click.echo(f"Result written to {result_path}")

    if dry_run:
        cmd_str = " ".join(payload["data"]["cmd"])
        click.echo(f"Would run: {cmd_str}")
    elif payload["status"] == "error":
        err_msg = payload.get("error") or "Unknown error"
        cat = payload.get("error_category")
        if cat == "input":
            code = 1
        elif cat == "env":
            code = 3
        else:
            code = 2
        _fail(err_msg, exit_code=code)

    if not quiet:
        kr = payload["key_results"]
        click.echo()
        click.echo(f"Status: {payload['status']}")
        click.echo(f"Wall time: {payload['wall_time']:.1f}s")
        if kr.get("n_candidate_trees"):
            click.echo(f"Candidate trees tested: {kr['n_candidate_trees']}")
        if kr.get("best_tree_id") is not None:
            click.echo(f"Best tree: #{kr['best_tree_id']}")
        if kr.get("n_rejected_au_0_05") is not None:
            click.echo(f"Rejected (AU < 0.05): {kr['n_rejected_au_0_05']}")
        _display_tests_table(payload["data"].get("tests", []))
