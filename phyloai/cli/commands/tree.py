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
        return ["fasttree", "iqtree"]


class _TreeGroup(click.Group):
    def list_commands(self, ctx: click.Context) -> list[str]:
        return ["ml"]


class _GroupedHelpCommand(click.Command):
    def format_epilog(self, ctx: click.Context, formatter: click.HelpFormatter) -> None:
        return None


class _IQTreeCommand(_GroupedHelpCommand):
    _HELP_GROUPS: list[tuple[str, list[str]]] = [
        (
            "Input",
            [
                "--msa-dir  Directory of MSA files for batch gene trees",
                "--matrix   Single concatenated matrix for supermatrix inference",
                "",
                "--msa-dir and --matrix are mutually exclusive.",
                "Exactly one of --msa-dir / --matrix is required.",
                "Supported formats: .fa .fas .fasta .faa .fna .phy .phylip .nex .nxs .nexus .aln",
            ],
        ),
        (
            "Data Type",
            [
                "--seq-type  AA | NT | auto (default: auto)",
                "",
                "auto validates input and lets IQ-TREE infer molecule type.",
                "NT maps to --seqtype DNA; other types (BIN, NT2AA, CODON, MORPH) via --tool-args.",
            ],
        ),
        (
            "Model",
            [
                "--model             Substitution model (omit for auto: AA->LG, NT->GTR)",
                "AA standard models:",
                "  LG, Poisson, cpREV, mtREV, Dayhoff, mtMAM, JTT, WAG, mtART, mtZOA, VT, rtREV, DCMut, PMB, HIVb, HIVw, JTTDCMut, FLU, Blosum62, GTR20, mtMet, mtVer, mtInv, FLAVI, Q.LG, Q.pfam, Q.pfam_gb, Q.bird, Q.mammal, Q.insect, Q.plant, Q.yeast",
                "NT standard models:",
                "  GTR, HKY, JC, F81, K2P, K3P, K81uf, TN, TNef, TIM, TIMef, TVM, TVMef, SYM",
                "AA heterogeneous / mixture models:",
                "  C10-C60, EX2, EX3, EHO, UL2, UL3, EX_EHO, LG4M, LG4X",
                "NT heterogeneous model:",
                "  MIX+MF",
                "--state-freq        +F|+FO|+FQ|+FU|none (default: +F)",
                "--rate-heterogeneity  +I|+G4|+I+G4|+R4|+I+R4|none (default: +R4)",
                "",
                "--model, --state-freq and --rate-heterogeneity combine to form the IQ-TREE -m argument string (e.g. LG+F+R4).",
                "+FU only valid for AA.  Ignored when --modelfinder is MF or MFP.",
                "--qmax is only valid with --model MIX+MF.",
            ],
        ),
        (
            "ModelFinder",
            [
                "--modelfinder  MF | MFP | none (default: none)",
                "  MF: model-only, no tree, no branch support",
                "  MFP: ModelFinder + tree search",
                "--mset  Restrict model search space (AA default: LG,WAG; NT default: GTR,HKY).",
                "  Use 'all' for unrestricted search.",
                "--msub  nuclear|mitochondrial|chloroplast|viral (default: nuclear; AA only).",
                "",
                "When --modelfinder is MF or MFP, --model/--state-freq/--rate-heterogeneity are ignored.",
                "Advanced freq/rate candidate control (--mfreq, --mrate) only via --tool-args.",
            ],
        ),
        (
            "Partitions",
            [
                "--partitions    Partition file (--matrix only)",
                "--rclusterf     Merge percentage for MF/MFP partition merging (default: 10)",
                "--rcluster-max  Max merge pairs; mutually exclusive with --rclusterf",
                "",
                "--partitions is only valid with --matrix.",
                "--rclusterf and --rcluster-max are mutually exclusive.",
                "If neither is provided with --partitions + MF/MFP, PhyloAI uses --rclusterf 10.",
            ],
        ),
        (
            "Heterogeneous",
            [
                "--pmsf-base-model  Base AA model for PMSF (default: LG)",
                "--guide-tree       Required to trigger PMSF with C10-C60",
                "--qmax             MIX+MF rate categories (default: 10)",
                "",
                "Heterogeneous workflows require --matrix.",
                "Direct AA mixture: choose a mixture model such as C20 or LG4X without --guide-tree. Example model string: C20+F+R4.",
                "PMSF AA mixture: choose C10-C60 and provide --guide-tree; PhyloAI defaults --pmsf-base-model to LG. Example model string: LG+C20+F+R4.",
                "NT heterogeneous: use --model MIX+MF; optional --mset restricts candidate models. Example model string: MIX+MF.",
            ],
        ),
        (
            "Tree Search",
            [
                "--mode  normal | fast",
                "",
                "fast maps to IQ-TREE --fast.",
            ],
        ),
        (
            "Branch Support",
            [
                "--boot  UFBoot replicates, >=1000 recommended (default: 1000). 0 = skip.",
                "--alrt  SH-aLRT replicates; 0 = parametric aLRT (optional)",
                "--bnni  NNI-optimize UFBoot trees; only when --boot is provided",
                "",
                "IQ-TREE mappings: boot -> -B, alrt -> --alrt, bnni -> --bnni.",
                "--bnni requires --boot > 0.  MF mode ignores boot/alrt/bnni.",
            ],
        ),
        (
            "Output",
            [
                "--rate        Write empirical Bayesian site rates to .rate file",
                "--wslr        Write site log-likelihoods per rate category to .sitelh file",
                "--constraint  Topological constraint tree (NEWICK). Maps to IQ-TREE -g.",
                "--outgroup    Comma-separated outgroup taxa. Maps to IQ-TREE -o.",
                "--prefix      Output prefix (--matrix only; ignored in batch).",
                "",
                "Single --matrix mode uses --prefix for all IQ-TREE native outputs.",
            ],
        ),
        (
            "Execution",
            [
                "-o, --output-dir  Output directory (default: runs/tree/ml/iqtree)",
                "--threads         Batch: parallel IQ-TREE jobs (default: 4); Single: NUM or auto (default: auto)",
                "--overwrite       Remove existing output directory before running",
                "--resume          Resume incomplete run",
                "--dry-run         Print commands without executing",
                "--keep-extra      Batch: keep extra IQ-TREE files (.ckp.gz, .bionj, etc.) in logs/",
                "-q, --quiet       Suppress all output except errors",
                "-h, --help        Show this message and exit",
                "--iqtree-path     Custom path to iqtree3 executable",
                "--tool-args       Extra IQ-TREE strategy flags (BLOCKED: -s, I/O redirects)",
                "",
                "Batch mode also blocks --prefix in --tool-args; other managed flags are overrideable.",
            ],
        ),
        (
            "Workflow Examples",
            [
                "Homogeneous batch fixed model:\n    phyloai tree ml iqtree --msa-dir msas/ --seq-type AA --model LG",
                "Homogeneous matrix fixed model:\n    phyloai tree ml iqtree --matrix matrix.fa --seq-type NT --model GTR",
                "ModelFinder only:\n    phyloai tree ml iqtree --matrix matrix.fa --seq-type AA --modelfinder MF --mset LG,WAG",
                "ModelFinder + tree:\n    phyloai tree ml iqtree --matrix matrix.fa --seq-type NT --modelfinder MFP --mset GTR,HKY",
                "Partitioned ModelFinder merge:\n    phyloai tree ml iqtree --matrix matrix.fa --seq-type AA --partitions parts.nex --modelfinder MFP",
                "Direct AA mixture:\n    phyloai tree ml iqtree --matrix matrix.fa --seq-type AA --model C20",
                "PMSF AA mixture:\n    phyloai tree ml iqtree --matrix matrix.fa --seq-type AA --model C20 --guide-tree guide.nwk",
                "NT heterogeneous:\n    phyloai tree ml iqtree --matrix matrix.fa --seq-type NT --model MIX+MF",
            ],
        ),
    ]

    def format_help(self, ctx: click.Context, formatter: click.HelpFormatter) -> None:
        import textwrap

        self.format_usage(ctx, formatter)
        self.format_help_text(ctx, formatter)
        formatter.write_paragraph()
        for section, rows in self._HELP_GROUPS:
            with formatter.section(section):
                formatter.write_paragraph()
                for row in rows:
                    if row == "":
                        formatter.write_paragraph()
                        continue
                    if section == "Workflow Examples":
                        formatter.write(f"  {row}\n")
                    else:
                        indent = formatter.current_indent
                        width = formatter.width - indent
                        if width < 40:
                            width = 70
                        wrapped = textwrap.fill(
                            row,
                            width=width,
                            initial_indent="",
                            subsequent_indent="  ",
                        )
                        indent_str = " " * indent
                        for line in wrapped.split("\n"):
                            formatter.buffer.append(indent_str + line + "\n")
        self.format_epilog(ctx, formatter)


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
        "  --msa-dir : batch gene trees from an MSA directory\n\n"
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
    "--gamma/--no-gamma",
    default=True,
    show_default=True,
    help="Enable gamma-distributed rate heterogeneity.",
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
        exit_code = 3 if "fasttree not found" in error_msg.lower() else 1
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


@ml.command(
    "iqtree",
    cls=_IQTreeCommand,
    help=(
        "Infer ML trees using IQ-TREE3.\n\n"
        "  --msa-dir : batch gene trees (homogeneous workflows only)\n\n"
        "  --matrix  : single supermatrix (all workflows: homogeneous, heterogeneous, partitioned)\n\n"
        "Reads FASTA, PHYLIP, NEXUS, CLUSTAL formats. "
        "Heterogeneous models (C10-C60, MIX+MF) require --matrix."
    ),
)
@click.option("--msa-dir", type=click.Path(file_okay=False, path_type=Path), default=None,
              help="Directory of MSA files for batch gene tree inference.")
@click.option("--matrix", type=click.Path(dir_okay=False, path_type=Path), default=None,
              help="Single concatenated matrix for supermatrix inference.")
@click.option("--seq-type", type=click.Choice(["AA", "NT", "auto"]), default="auto", show_default=True,
              help="Molecule type. NT maps to --seqtype DNA.")
@click.option("--model", type=str, default=None,
              help="Substitution model. Omit for auto-detect: AA->LG, NT->GTR. Heterogeneous: C10-C60, MIX+MF.")
@click.option("--state-freq", type=click.Choice(["+F", "+FO", "+FQ", "+FU", "none"]),
              default="+F", show_default=True, help="State frequency type.")
@click.option("--rate-heterogeneity", type=click.Choice(["+I", "+G4", "+I+G4", "+R4", "+I+R4", "none"]),
              default="+R4", show_default=True, help="Rate heterogeneity among sites.")
@click.option("--modelfinder", type=click.Choice(["MF", "MFP", "none"]), default="none", show_default=True,
              help="MF = ModelFinder only (no tree). MFP = ModelFinder + tree search.")
@click.option("--mset", type=str, default=None,
              help="Comma-separated model list for ModelFinder. AA default: LG,WAG. NT default: GTR,HKY.")
@click.option("--msub", type=click.Choice(["nuclear", "mitochondrial", "chloroplast", "viral"]),
              default=None, help="AA model source. AA only.")
@click.option("--mode", type=click.Choice(["normal", "fast"]), default="normal", show_default=True,
              help="Tree search mode. fast maps to --fast.")
@click.option("--boot", type=click.IntRange(min=0), default=1000, show_default=True,
              help="UFBoot replicates. 0 = skip branch support.")
@click.option("--alrt", type=click.IntRange(min=0), default=None,
              help="SH-aLRT replicates. 0 = parametric aLRT.")
@click.option("--bnni", is_flag=True, default=False,
              help="Optimize UFBoot trees by NNI.")
@click.option("--partitions", type=click.Path(dir_okay=False, path_type=Path), default=None,
              help="Partition file. --matrix only.")
@click.option("--rclusterf", type=click.IntRange(1, 100), default=None,
              help="Percent partition pairs for rclusterf merge. Partitions + MF/MFP default: 10.")
@click.option("--rcluster-max", type=int, default=None,
              help="Max partition pairs for rcluster merge. Mutually exclusive with --rclusterf.")
@click.option("--pmsf-base-model", type=str, default=None,
              help="Base AA model for PMSF (C10-C60 only). Default: LG. Requires --guide-tree.")
@click.option("--guide-tree", type=click.Path(dir_okay=False, path_type=Path), default=None,
              help="Guide tree for PMSF in NEWICK format. Only with --model C10-C60.")
@click.option("--qmax", type=int, default=None,
              help="Max rate categories for MIX+MF (default: 10).")
@click.option("--rate", is_flag=True, default=False,
              help="Write site rates to .rate file.")
@click.option("--wslr", is_flag=True, default=False,
              help="Write site log-likelihoods to .sitelh file.")
@click.option("--constraint", type=click.Path(dir_okay=False, path_type=Path), default=None,
              help="Topological constraint tree in NEWICK format.")
@click.option("--outgroup", type=str, default=None,
              help="Outgroup taxa, comma-separated.")
@click.option("--prefix", type=str, default=None,
              help="Prefix for output files. --matrix only; ignored in batch.")
@click.option("--output-dir", "-o", type=click.Path(file_okay=False, path_type=Path),
              default=Path("runs/tree/ml/iqtree"), show_default=True, help="Output directory.")
@click.option("--threads", "-t", type=str, default=None,
              help="Batch: parallel IQ-TREE jobs (default: 4). Single: NUM or 'auto' (default: auto).")
@click.option("--iqtree-path", type=Path, default=None,
              help="Explicit path to iqtree3 executable.")
@click.option("--tool-args", type=str, default=None,
              help="Extra IQ-TREE flags. I/O flags blocked; strategy flags override PhyloAI defaults.")
@click.option("--overwrite", is_flag=True, default=False, help="Overwrite existing output directory.")
@click.option("--resume", is_flag=True, default=False, help="Resume incomplete run.")
@click.option("--dry-run", is_flag=True, default=False, help="Show commands without executing.")
@click.option("--keep-extra", is_flag=True, default=False,
              help="Batch mode: keep extra IQ-TREE output files (.ckp.gz, .bionj, .mldist, etc.) in logs/.")
@click.option("--quiet", "-q", is_flag=True, default=False, help="Suppress terminal output except errors.")
def iqtree_command(
    msa_dir: Path | None,
    matrix: Path | None,
    seq_type: str,
    model: str,
    state_freq: str,
    rate_heterogeneity: str,
    modelfinder: str,
    mset: str | None,
    msub: str | None,
    mode: str,
    boot: int | None,
    alrt: int | None,
    bnni: bool,
    partitions: Path | None,
    rclusterf: int | None,
    rcluster_max: int | None,
    pmsf_base_model: str | None,
    guide_tree: Path | None,
    qmax: int | None,
    rate: bool,
    wslr: bool,
    constraint: Path | None,
    outgroup: str | None,
    prefix: str | None,
    output_dir: Path,
    threads: str | None,
    iqtree_path: Path | None,
    tool_args: str | None,
    overwrite: bool,
    resume: bool,
    dry_run: bool,
    keep_extra: bool,
    quiet: bool,
) -> None:
    from phyloai.tree.ml_iqtree import run_iqtree
    from phyloai.tree.ml_iqtree import _scan_input_iqtree as _scan_iqtree

    batch_mode = msa_dir is not None
    single_mode = matrix is not None

    if batch_mode == single_mode:
        if not batch_mode and not single_mode:
            _fail("Either --msa-dir or --matrix is required.", 1)
        else:
            _fail("--msa-dir and --matrix are mutually exclusive.", 1)

    if resume and overwrite:
        _fail("--overwrite and --resume are mutually exclusive.", 1)

    if msa_dir is not None and not msa_dir.exists():
        _fail(f"--msa-dir does not exist: {msa_dir}", 1)
    if matrix is not None and not matrix.exists():
        _fail(f"--matrix does not exist: {matrix}", 1)
    if partitions is not None and not partitions.exists():
        _fail(f"--partitions does not exist: {partitions}", 1)
    if guide_tree is not None and not guide_tree.exists():
        _fail(f"--guide-tree does not exist: {guide_tree}", 1)
    if constraint is not None and not constraint.exists():
        _fail(f"--constraint does not exist: {constraint}", 1)

    if iqtree_path is not None:
        if not iqtree_path.exists():
            _fail(f"--iqtree-path does not exist: {iqtree_path}", 1)
        if not os.access(str(iqtree_path), os.X_OK):
            _fail(f"--iqtree-path is not executable: {iqtree_path}", 1)

    if bnni and boot is None:
        if not quiet:
            click.echo("Warning: --bnni has no effect without --boot.", err=True)

    # Warn about MF mode branch support (CLI-level validation)
    if modelfinder == "MF" and (boot is not None or alrt is not None or bnni):
        if not quiet:
            click.echo(
                "Warning: --boot/--alrt/--bnni are ignored in MF (model-only) mode.", err=True
            )

    iqtree_path_str = str(iqtree_path) if iqtree_path else None
    partitions_str = str(partitions) if partitions else None
    guide_tree_str = str(guide_tree) if guide_tree else None
    constraint_str = str(constraint) if constraint else None

    # Resolve pmsf_base_model default: only when PMSF is triggered
    from phyloai.tree.ml_iqtree import AA_MIXTURE_MODELS
    pmsf_base_model_resolved = pmsf_base_model
    if pmsf_base_model_resolved is None and model in AA_MIXTURE_MODELS and guide_tree is not None:
        pmsf_base_model_resolved = "LG"

    def _invoke(progress_callback=None):
        return run_iqtree(
            msa_dir=msa_dir, matrix=matrix,
            output_dir=output_dir,
            seq_type=seq_type, model=model,
            state_freq=state_freq, rate_heterogeneity=rate_heterogeneity,
            modelfinder=modelfinder,
            mset=mset, msub=msub, mode=mode,
            boot=boot, alrt=alrt, bnni=bnni,
            partitions=partitions_str,
            rclusterf=rclusterf, rcluster_max=rcluster_max,
            pmsf_base_model=pmsf_base_model_resolved,
            guide_tree=guide_tree_str, qmax=qmax,
            rate=rate, wslr=wslr,
            constraint=constraint_str, outgroup=outgroup,
            prefix=prefix, threads=threads,
            iqtree_path=iqtree_path_str, tool_args=tool_args,
            overwrite=overwrite, resume=resume,
            dry_run=dry_run, keep_extra=keep_extra, quiet=quiet,
            progress_callback=progress_callback,
        )

    error_msg: str | None = None

    try:
        if not quiet and not dry_run and batch_mode:
            found, _ = _scan_iqtree(msa_dir)
            total = len(found)
            if total == 0:
                _fail("No valid input files found in --msa-dir.", 1)

            with Progress(console=console, transient=True) as progress:
                task = progress.add_task(
                    "[cyan]Inferring gene trees with IQ-TREE", total=total
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
        exit_code = 3 if "iqtree3 not found" in error_msg.lower() else 1
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
        click.echo(f"Trees: {n_trees} | Failed: {n_failed} | Skipped: {n_skipped}")
        if batch_mode:
            click.echo(f"Trees saved to {output_dir / 'trees'}", err=True)
            click.echo(f"Logs saved to {output_dir / 'logs'}", err=True)
        click.echo(f"Results saved to {result_path}", err=True)
        if n_failed > 0:
            click.echo(
                f"Warning: {n_failed} gene(s) failed. Check result.json data.failed for details.",
                err=True,
            )

    if n_trees == 0 and n_failed > 0:
        _fail("All IQ-TREE runs failed.", 2)
