"""Post-tree analysis CLI commands."""
from __future__ import annotations

import json
import shlex
import shutil
from pathlib import Path

import click

from phyloai.core.iqtree import IQTREE_COMPATIBLE_EXTENSIONS


class _PosttreeGroup(click.Group):
    def list_commands(self, ctx: click.Context) -> list[str]:
        return ["topology"]


@click.group(cls=_PosttreeGroup)
def posttree() -> None:
    """Post-tree analysis commands."""


def _fail(message: str, exit_code: int = 1) -> None:
    click.echo(f"Error: {message}", err=True)
    raise click.exceptions.Exit(exit_code)


_TOPOLOGY_HELP = """
Tree topology tests (AU / KH / SH / WKH / WSH / c-ELW).

Compares a set of candidate trees against a supermatrix alignment using
IQ-TREE's built-in topology test framework.

\\b
PURPOSE
  This command tests whether alternative topologies are significantly worse
  than the best-scoring candidate. It does NOT infer new trees — use
  `phyloai tree ml iqtree` for ML tree inference and model selection.

\\b
INPUT
  --matrix             Single supermatrix alignment (FASTA/PHYLIP/NEXUS/CLUSTAL).
  --candidate-trees    Accepts either one tree-list file (one NEWICK tree per
                       line) or multiple individual NEWICK tree files passed
                       in order. Multiple files are merged by PhyloAI into
                       candidate.trees before invoking IQ-TREE.

\\b
MODEL SOURCE  (exactly one required)
  --model-expr          Complete IQ-TREE -m expression (e.g. LG+F+R4, C20+F+R4).
  --partitions PATH     Previously optimized partition model (e.g. .best_model.nex).
                        Maps to IQ-TREE -p.
  --guide-tree PATH     Guide tree for PMSF models. Maps to IQ-TREE -ft.

  PhyloAI does NOT expose ModelFinder here because topology tests are run
  after model inference — you should already have a preferred model.

\\b
DEFAULT TESTS
  PhyloAI generates standard topology-test flags:

      -n 0 -zb <replicates> -zw -au

  This produces: bp-RELL, KH, SH, weighted KH, weighted SH, c-ELW, and AU.
  Each test is individually suppressible via --tool-args.

\\b
ADVANCED IQ-TREE ARGS
  --tool-args TEXT    Additional IQ-TREE strategy parameters. Blocked flags:
                      -s, -z (managed by --matrix and --candidate-trees).
  --iqtree-path PATH  Explicit path to iqtree3 executable.

  PhyloAI-built flags are suppressed when the same flag appears in
  --tool-args (suppress-if-present).  Overrideable: -m, -p, -ft, -n, -zb,
  -zw, -au, -T, --prefix.

\\b
EXAMPLES

  # Homogeneous unpartitioned model
  phyloai posttree topology --matrix raw.fa --candidate-trees trees \\
      --model-expr LG+F+R4 --replicates 10000 -t 20

  # Heterogeneous model
  phyloai posttree topology --matrix raw.fa --candidate-trees trees \\
      --model-expr C20+F+R4 -t 20

  # PMSF model with guide tree
  phyloai posttree topology --matrix raw.fa --candidate-trees trees \\
      --model-expr LG+C20+F+R4 --guide-tree guide.tree -t 4

  # Previously optimized partition model
  phyloai posttree topology --matrix raw.fa --candidate-trees trees \\
      --partitions raw.best_model.nex -t 20

  # Multiple individual tree files (merged by PhyloAI)
  phyloai posttree topology --matrix raw.fa \\
      --candidate-trees h1.nwk --candidate-trees h2.nwk \\
      --candidate-trees h3.nwk --model-expr LG+F+R4 -t 20

  # Custom exchangeabilities + site frequencies via --tool-args
  phyloai posttree topology --matrix raw.fa --candidate-trees trees \\
      --model-expr custom.exchangeabilities+R4 \\
      --tool-args "-fs custom.sitefreq" -t 30

\\b
INPUT FORMAT AND SEQUENCE TYPE
  --input-format only affects PhyloAI's own matrix preflight validation;
  it is NOT passed to IQ-TREE.  IQ-TREE's --seqtype flag can be passed
  via --tool-args when needed (e.g. --tool-args "--seqtype AA").

\\b
INTERPRETATION
  KH / SH / WKH / WSH / AU are p-values.  Trees with p < 0.05 are rejected
  by that test.  bp-RELL and c-ELW are weights (not p-values).  The AU test
  is generally considered the most reliable.
"""


@posttree.command("topology", help=_TOPOLOGY_HELP)
@click.option(
    "--matrix", type=click.Path(path_type=Path), default=None,
    help="Single supermatrix alignment (FASTA/PHYLIP/NEXUS/CLUSTAL).  Maps to IQ-TREE -s.",
)
@click.option(
    "--candidate-trees", "candidate_trees_raw",
    multiple=True, type=click.Path(path_type=Path),
    help=(
        "Candidate tree input. Accepts either one tree-list file (one NEWICK tree "
        "per line) or multiple individual NEWICK tree files (merged in order by PhyloAI)."
    ),
)
@click.option(
    "--input-format",
    type=click.Choice(["auto", "fasta", "phylip-relaxed", "nexus", "clustal"]),
    default="auto",
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
    candidate_trees_raw: tuple[Path, ...],
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
    """Tree topology tests (AU / KH / SH / WKH / WSH / c-ELW)."""
    from phyloai.posttree.topology import run_topology

    # ---- Manual validation (exit code 1 for all user input errors) ----

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

    # Candidate trees (existence / non-empty / readability)
    candidate_trees_list = list(candidate_trees_raw)
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
        click.echo(f"Status: {payload['status']}")
        click.echo(f"Wall time: {payload['wall_time']:.1f}s")
        if kr.get("n_candidate_trees"):
            click.echo(f"Candidate trees tested: {kr['n_candidate_trees']}")
        if kr.get("best_tree_id") is not None:
            click.echo(f"Best tree: #{kr['best_tree_id']}")
        if kr.get("n_rejected_au_0_05") is not None:
            click.echo(f"Rejected (AU < 0.05): {kr['n_rejected_au_0_05']}")
