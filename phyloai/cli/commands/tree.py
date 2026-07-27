"""Tree inference CLI commands."""

from __future__ import annotations

import json
import os
import shutil
from pathlib import Path

import click
from rich.console import Console
from rich.progress import Progress

from phyloai.tree.ml import run_fasttree

console = Console()


def _fail(message: str, exit_code: int) -> None:
    click.echo(f"Error: {message}", err=True)
    raise click.exceptions.Exit(exit_code)


def _write_error_result(kwargs: dict, message: str, exit_code: int) -> None:
    output_dir = kwargs.get("output_dir", Path("runs/tree/bi"))
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    result = {
        "status": "error",
        "command": " ".join(_build_command_tokens(kwargs)),
        "wall_time": 0.0,
        "tool_versions": {},
        "params": {
            k: str(v) if isinstance(v, Path) else v
            for k, v in kwargs.items()
        },
        "key_results": {},
        "error": message,
        "data": {},
    }
    result_path = output_dir / "result.json"
    result_path.write_text(json.dumps(result, indent=2))


def _build_command_tokens(kwargs: dict) -> list[str]:
    tokens = ["phyloai", "tree", "bi", "pb"]
    flag_map = {
        "matrix": "--matrix",
        "output_dir": "--output-dir",
        "model": "--model",
        "mixture": "--mixture",
        "gamma_cats": "--gamma-cats",
        "start_tree": "--start-tree",
        "fix_tree": "--fix-tree",
        "chains": "--chains",
        "chain_prefix": "--chain-prefix",
        "chain_names": "--chain-names",
        "threads": "--threads",
        "sample_freq": "--sample-freq",
        "nsamples": "--nsamples",
        "resume": "--resume",
        "monitor_freq": "--monitor-freq",
        "burnin_frac": "--burnin-frac",
        "poll_interval": "--poll-interval",
        "pb_path": "--pb-path",
    }
    for key, flag in flag_map.items():
        val = kwargs.get(key)
        if val is None:
            continue
        if isinstance(val, bool) and val:
            tokens.append(flag)
        elif isinstance(val, bool):
            continue
        else:
            tokens.append(flag)
            tokens.append(str(val))
    if kwargs.get("overwrite"):
        tokens.append("--overwrite")
    if kwargs.get("dry_run"):
        tokens.append("--dry-run")
    if kwargs.get("quiet"):
        tokens.append("--quiet")
    return tokens


_BPCOMP_FLAG_MAP: dict[str, str] = {
    "chain_dir": "--chain-dir", "chain_names": "--chain-names",
    "output_dir": "--output-dir", "burnin": "--burnin",
    "sample_freq": "--sample-freq", "until": "--until",
    "cutoff": "--cutoff", "pb_path": "--pb-path",
}
_TRACECOMP_FLAG_MAP: dict[str, str] = {
    "chain_dir": "--chain-dir", "chain_names": "--chain-names",
    "output_dir": "--output-dir", "burnin": "--burnin",
    "pb_path": "--pb-path",
}
_READPB_FLAG_MAP: dict[str, str] = {
    "chain": "--chain", "mode": "--mode",
    "output_dir": "--output-dir", "burnin": "--burnin",
    "sample_freq": "--sample-freq", "until": "--until",
    "threads": "--threads", "pb_path": "--pb-path",
}


def _build_bi_subcommand_tokens(subcommand: str, kwargs: dict, flag_map: dict[str, str]) -> list[str]:
    tokens = ["phyloai", "tree", "bi", subcommand]
    for key, flag in flag_map.items():
        val = kwargs.get(key)
        if val is None:
            continue
        if isinstance(val, bool) and val:
            tokens.append(flag)
        elif isinstance(val, bool):
            continue
        else:
            tokens.append(flag)
            tokens.append(str(val))
    if kwargs.get("overwrite"):
        tokens.append("--overwrite")
    if kwargs.get("dry_run"):
        tokens.append("--dry-run")
    if kwargs.get("quiet"):
        tokens.append("--quiet")
    return tokens


def _write_bi_error_result(kwargs: dict, message: str, exit_code: int, subcommand: str, flag_map: dict[str, str]) -> None:
    output_dir = Path(kwargs.get("output_dir", f"runs/tree/bi/{subcommand}"))
    output_dir.mkdir(parents=True, exist_ok=True)
    tokens = _build_bi_subcommand_tokens(subcommand, kwargs, flag_map)
    result = {
        "status": "error",
        "command": " ".join(tokens),
        "wall_time": 0.0,
        "tool_versions": {},
        "params": {k: (str(v) if isinstance(v, Path) else v) for k, v in kwargs.items()},
        "key_results": {},
        "error": message,
        "data": {"cmd": [], "tool_stderr": ""},
    }
    (output_dir / "result.json").write_text(json.dumps(result, indent=2))


class _MLGroup(click.Group):
    def list_commands(self, ctx: click.Context) -> list[str]:
        return ["fasttree", "iqtree"]


class _BiGroup(click.Group):
    def list_commands(self, ctx: click.Context) -> list[str]:
        return ["pb", "bpcomp", "tracecomp", "readpb"]

    def format_help_text(self, ctx: click.Context, formatter: click.HelpFormatter) -> None:
        text = self.help or ""
        if not text:
            return
        paragraphs = text.split("\n\n")
        for i, paragraph in enumerate(paragraphs):
            if i:
                formatter.write_paragraph()
            formatter.write(paragraph)


class _TreeGroup(click.Group):
    def list_commands(self, ctx: click.Context) -> list[str]:
        return ["ml", "bi", "msc", "cf"]


class _GroupedHelpCommand(click.Command):
    def format_epilog(self, ctx: click.Context, formatter: click.HelpFormatter) -> None:
        return None


class _CFCommand(_GroupedHelpCommand):
    """Preserve newline formatting in help text without Click rewrapping."""

    def format_help_text(self, ctx: click.Context, formatter: click.HelpFormatter) -> None:
        text = self.help or ""
        if not text:
            return
        paragraphs = text.split("\n\n")
        for i, paragraph in enumerate(paragraphs):
            if i:
                formatter.write_paragraph()
            formatter.write(paragraph)


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
            if batch_mode:
                click.echo(
                    f"Dry run: {payload['data']['summary']['n_input_files']} input(s) would be processed."
                )
                for item in payload["data"].get("files", []):
                    if "cmd" in item:
                        click.echo(" ".join(item["cmd"]))
            else:
                click.echo("Dry run: single matrix mode would be processed.")
                if payload["data"].get("cmd"):
                    click.echo(" ".join(payload["data"]["cmd"]))
        return

    result_path = output_dir / "result.json"
    result_path.parent.mkdir(parents=True, exist_ok=True)
    with open(result_path, "w") as fh:
        json.dump(payload, fh, indent=2)

    if batch_mode:
        summary = payload["data"]["summary"]
        n_failed = summary.get("n_failed", 0)
        n_trees = summary.get("n_trees", 0)
        n_skipped = summary.get("n_skipped", 0)
    else:
        n_trees = 1 if payload["data"].get("output") else 0
        n_failed = 0 if n_trees else 1
        n_skipped = 0

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
            if batch_mode:
                click.echo(
                    f"Dry run: {payload['data']['summary']['n_input_files']} input(s) would be processed."
                )
                for item in payload["data"].get("files", []):
                    if "cmd" in item:
                        click.echo(" ".join(item["cmd"]))
            else:
                click.echo("Dry run: single matrix mode would be processed.")
                if payload["data"].get("cmd"):
                    click.echo(" ".join(payload["data"]["cmd"]))
        return

    result_path = output_dir / "result.json"
    result_path.parent.mkdir(parents=True, exist_ok=True)
    with open(result_path, "w") as fh:
        json.dump(payload, fh, indent=2)

    if batch_mode:
        summary = payload["data"]["summary"]
        n_failed = summary.get("n_failed", 0)
        n_trees = summary.get("n_trees", 0)
        n_skipped = summary.get("n_skipped", 0)
    else:
        n_trees = 1 if payload["data"].get("output") else 0
        n_failed = 0 if n_trees else 1
        n_skipped = 0

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


@tree.command(
    "msc",
    cls=_GroupedHelpCommand,
    help=(
        "Multispecies coalescent species tree inference with wASTRAL.\n\n"
        "  --tree     : single gene tree file (newick, one tree per line)\n\n"
        "  --tree-dir : directory of gene tree files (merged into one input)\n\n"
        "--tree and --tree-dir are mutually exclusive. "
        "wASTRAL is one-shot computation (no --resume)."
    ),
)
@click.option(
    "--tree",
    type=click.Path(dir_okay=False, path_type=Path),
    default=None,
    help="Single gene tree file (newick, one tree per line).",
)
@click.option(
    "--tree-dir",
    type=click.Path(file_okay=False, path_type=Path),
    default=None,
    help="Directory of gene tree files for merging.",
)
@click.option(
    "--mode",
    type=click.IntRange(1, 4),
    default=1,
    show_default=True,
    help="wastral --mode. 1=hybrid, 2=branch support weighting, 3=branch length weighting, 4=traditional unweighted Astral.",
)
@click.option(
    "--boot",
    type=click.IntRange(0, 3),
    default=1,
    show_default=True,
    help="wastral -u/--support. 0=topology only, 1=local posterior probability, 2=quartet+local-PP, 3=2+freqQuad.csv.",
)
@click.option(
    "--extra-rounds", "-R",
    is_flag=True,
    default=False,
    help="Enable exhaustive search (wastral -R).",
)
@click.option(
    "--tree-boot-type",
    type=click.Choice(["auto", "likelihood", "abayes", "bootstrap"]),
    default="auto",
    show_default=True,
    help="Gene tree branch support type wastral preset: "
    "auto (detect from gene trees), "
    "likelihood (wastral -L/--lrt: alrt, -x 1 -n 0), "
    "abayes (wastral -B/--bayes, -x 1 -n 0.333), "
    "bootstrap (wastral -S/--bootstrap: default, -x 100 -n 0).",
)
@click.option(
    "--tree-boot-min",
    type=float,
    default=None,
    help="Minimum support threshold (wastral -n). Only with non-auto --tree-boot-type.",
)
@click.option(
    "--tree-boot-max",
    type=float,
    default=None,
    help="Maximum support value (wastral -x). Only with non-auto --tree-boot-type.",
)
@click.option(
    "--outgroup",
    type=str,
    default=None,
    help="Outgroup species for rooting (wastral --root).",
)
@click.option(
    "--output-dir", "-o",
    type=click.Path(file_okay=False, path_type=Path),
    default=Path("runs/tree/msc"),
    show_default=True,
    help="Output directory.",
)
@click.option(
    "--threads", "-t",
    type=int,
    default=4,
    show_default=True,
    help="Thread count for wastral -t.",
)
@click.option(
    "--wastral-path",
    type=Path,
    default=None,
    help="Explicit path to wastral executable.",
)
@click.option(
    "--tool-args",
    type=str,
    default=None,
    help="Extra wastral flags. -i/-o blocked; strategy flags override phyloAI defaults.",
)
@click.option("--overwrite", is_flag=True, default=False, help="Overwrite existing output directory.")
@click.option("--dry-run", is_flag=True, default=False, help="Show commands without executing.")
@click.option("--quiet", "-q", is_flag=True, default=False, help="Suppress terminal output except errors.")
def msc_command(
    tree: Path | None,
    tree_dir: Path | None,
    mode: int,
    boot: int,
    extra_rounds: bool,
    tree_boot_type: str,
    tree_boot_min: float | None,
    tree_boot_max: float | None,
    outgroup: str | None,
    output_dir: Path,
    threads: int,
    wastral_path: Path | None,
    tool_args: str | None,
    overwrite: bool,
    dry_run: bool,
    quiet: bool,
) -> None:
    from phyloai.tree.msc import run_wastral

    # Mutual exclusivity: CLI-layer early check for better error messages
    if (tree is None and tree_dir is None) or (tree is not None and tree_dir is not None):
        _fail("Either --tree or --tree-dir must be provided (mutually exclusive).", 1)

    if tree is not None and not tree.exists():
        _fail(f"--tree does not exist: {tree}", 1)
    if tree_dir is not None and not tree_dir.exists():
        _fail(f"--tree-dir does not exist: {tree_dir}", 1)

    if wastral_path is not None:
        if not wastral_path.exists():
            _fail(f"--wastral-path does not exist: {wastral_path}", 1)
        if not os.access(str(wastral_path), os.X_OK):
            _fail(f"--wastral-path is not executable: {wastral_path}", 1)

    error_msg: str | None = None

    try:
        payload = run_wastral(
            tree=tree,
            tree_dir=tree_dir,
            output_dir=output_dir,
            mode=mode,
            boot=boot,
            extra_rounds=extra_rounds,
            tree_boot_type=tree_boot_type,
            tree_boot_min=tree_boot_min,
            tree_boot_max=tree_boot_max,
            outgroup=outgroup,
            threads=threads,
            wastral_path=str(wastral_path) if wastral_path else None,
            tool_args=tool_args,
            overwrite=overwrite,
            dry_run=dry_run,
            quiet=quiet,
        )
    except (ValueError, FileNotFoundError) as exc:
        error_msg = str(exc)
    except SystemExit:
        raise
    except Exception as exc:
        error_msg = str(exc)

    if error_msg is not None:
        exit_code = 3 if "wastral not found" in error_msg.lower() else 1
        _fail(error_msg, exit_code)

    if dry_run:
        if not quiet:
            click.echo(
                f"Dry run: {payload['key_results']['n_input_trees']} gene tree(s) "
                f"would be processed."
            )
            click.echo(" ".join(payload["data"]["cmd"]))
        return

    # Write result.json (always, even on failure — design requires structured output)
    result_path = output_dir / "result.json"
    result_path.parent.mkdir(parents=True, exist_ok=True)
    with open(result_path, "w") as fh:
        json.dump(payload, fh, indent=2)

    # Check payload status for wastral execution failures (exit 2, after writing result.json)
    if payload["status"] == "error":
        _fail(payload.get("error", "wastral execution failed"), 2)

    if not quiet:
        click.echo(
            f"Species tree saved to {output_dir / 'wastral.tre'}"
        )
        click.echo(f"Results saved to {result_path}", err=True)


@tree.command(
    "cf",
    cls=_CFCommand,
    help=(
        "Compute concordance factors for a reference species tree.\n"
        "\n"
        "  Modes (--cf):\n"
        "    gcf      Gene concordance factor (IQ-TREE3)\n"
        "    scf      Site concordance factor, parsimony-based (IQ-TREE3)\n"
        "    scfl     Site concordance factor, likelihood-based (IQ-TREE3)\n"
        "    gcf+scf  Combined gCF + sCF in one IQ-TREE3 invocation\n"
        "    qcf      Quartet concordance factor (wASTRAL)\n"
        "\n"
        "  Input requirements by mode:\n"
        "    gcf / qcf       --ref-tree + (--tree or --tree-dir)\n"
        "    scf / scfl      --ref-tree + --matrix\n"
        "    gcf+scf         --ref-tree + (--tree or --tree-dir) + --matrix\n"
        "    scfl            optionally --model-expr or --partitions for speedup\n"
        "\n"
        "CF computation is one-shot (no --resume)."
    ),
)
@click.option(
    "--cf",
    type=click.Choice(["gcf", "scf", "scfl", "gcf+scf", "qcf"]),
    required=True,
    help="Concordance factor type to compute.",
)
@click.option(
    "--ref-tree",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    required=True,
    help="Reference species tree (NEWICK).",
)
@click.option(
    "--tree",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    default=None,
    help="Single gene tree file (NEWICK, one tree per line). Mutually exclusive with --tree-dir.",
)
@click.option(
    "--tree-dir",
    type=click.Path(exists=True, file_okay=False, path_type=Path),
    default=None,
    help="Directory of gene tree files (merged into merged.trees). Mutually exclusive with --tree.",
)
@click.option(
    "--matrix",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    default=None,
    help="Multiple sequence alignment (required for scf/scfl/gcf+scf).",
)
@click.option(
    "--partitions",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    default=None,
    help="Partition file for scfl model reuse (e.g., *.best_model.nex from IQ-TREE3).",
)
@click.option(
    "--model-expr",
    type=str,
    default=None,
    help="Substitution model for scfl speedup (e.g., LG+F+R4). Mutually exclusive with --partitions.",
)
@click.option(
    "--scf-quartets",
    type=click.IntRange(1, None),
    default=100,
    show_default=True,
    help="Number of quartets for sCF/sCFl (recommend >= 100).",
)
@click.option(
    "--lpp",
    is_flag=True,
    default=False,
    help="Append local posterior probabilities (pp1) to qCF support labels.",
)
@click.option(
    "--prefix",
    type=str,
    default=None,
    help="Output file prefix (default: auto-derived from --cf, e.g., gCF, sCF, sCFl, gCFsCF, qCF).",
)
@click.option(
    "--output-dir", "-o",
    type=click.Path(file_okay=False, path_type=Path),
    default=Path("runs/tree/cf"),
    show_default=True,
    help="Output directory.",
)
@click.option(
    "--threads", "-t",
    type=int,
    default=4,
    show_default=True,
    help="Thread count (IQ-TREE3 -T or wASTRAL -t).",
)
@click.option(
    "--iqtree-path",
    type=Path,
    default=None,
    help="Explicit path to iqtree3 executable.",
)
@click.option(
    "--wastral-path",
    type=Path,
    default=None,
    help="Explicit path to wastral executable.",
)
@click.option("--overwrite", is_flag=True, default=False, help="Overwrite existing output directory.")
@click.option("--dry-run", is_flag=True, default=False, help="Show commands without executing.")
@click.option("--quiet", "-q", is_flag=True, default=False, help="Suppress terminal output except errors.")
def cf_command(
    cf: str,
    ref_tree: Path,
    tree: Path | None,
    tree_dir: Path | None,
    matrix: Path | None,
    partitions: Path | None,
    model_expr: str | None,
    scf_quartets: int,
    prefix: str | None,
    output_dir: Path,
    threads: int,
    iqtree_path: Path | None,
    wastral_path: Path | None,
    overwrite: bool,
    dry_run: bool,
    quiet: bool,
    lpp: bool,
) -> None:
    """Concordance factor computation (gCF, sCF, sCFl, qCF)."""
    from phyloai.tree.cf import run_cf

    if iqtree_path is not None:
        if not iqtree_path.exists():
            _fail(f"--iqtree-path does not exist: {iqtree_path}", 1)
        if not os.access(str(iqtree_path), os.X_OK):
            _fail(f"--iqtree-path is not executable: {iqtree_path}", 1)
    if wastral_path is not None:
        if not wastral_path.exists():
            _fail(f"--wastral-path does not exist: {wastral_path}", 1)
        if not os.access(str(wastral_path), os.X_OK):
            _fail(f"--wastral-path is not executable: {wastral_path}", 1)

    error_msg: str | None = None

    try:
        payload = run_cf(
            cf_mode=cf,
            ref_tree=ref_tree,
            tree=tree,
            tree_dir=tree_dir,
            matrix=matrix,
            partitions=partitions,
            model_expr=model_expr,
            scf_quartets=scf_quartets,
            prefix=prefix,
            output_dir=output_dir,
            threads=threads,
            iqtree_path=str(iqtree_path) if iqtree_path else None,
            wastral_path=str(wastral_path) if wastral_path else None,
            overwrite=overwrite,
            dry_run=dry_run,
            quiet=quiet,
            lpp=lpp,
        )
    except (ValueError, FileNotFoundError) as exc:
        error_msg = str(exc)
    except SystemExit:
        raise
    except Exception as exc:
        error_msg = str(exc)

    if error_msg is not None:
        if "iqtree3 not found" in error_msg.lower() or "wastral not found" in error_msg.lower():
            exit_code = 3
        else:
            exit_code = 1
        try:
            import json as _json
            output_dir.mkdir(parents=True, exist_ok=True)

            # result.json
            _prefix = prefix or {"gcf": "gCF", "scf": "sCF", "scfl": "sCFl", "gcf+scf": "gCFsCF", "qcf": "qCF"}.get(cf, cf)
            _cmd_parts = ["phyloai", "tree", "cf", "--cf", cf]
            _cmd_parts.extend(["--ref-tree", str(ref_tree)])
            if tree is not None:
                _cmd_parts.extend(["--tree", str(tree)])
            elif tree_dir is not None:
                _cmd_parts.extend(["--tree-dir", str(tree_dir)])
            if matrix is not None:
                _cmd_parts.extend(["--matrix", str(matrix)])
            if partitions is not None:
                _cmd_parts.extend(["--partitions", str(partitions)])
            if model_expr is not None:
                _cmd_parts.extend(["--model-expr", model_expr])
            if cf not in ("gcf", "qcf"):
                _cmd_parts.extend(["--scf-quartets", str(scf_quartets)])
            if lpp:
                _cmd_parts.append("--lpp")
            _cmd_parts.extend(["--prefix", _prefix])
            _cmd_parts.extend(["-o", str(output_dir)])
            _cmd_parts.extend(["-t", str(threads)])
            if overwrite:
                _cmd_parts.append("--overwrite")
            if iqtree_path:
                _cmd_parts.extend(["--iqtree-path", str(iqtree_path)])
            if wastral_path:
                _cmd_parts.extend(["--wastral-path", str(wastral_path)])
            if dry_run:
                _cmd_parts.append("--dry-run")
            _cmd_str = " ".join(_cmd_parts)

            result = {
                "status": "error",
                "command": _cmd_str,
                "wall_time": 0.0,
                "tool_versions": {},
                "params": {
                    "cf": cf,
                    "ref_tree": str(ref_tree),
                    "tree": str(tree) if tree else None,
                    "tree_dir": str(tree_dir) if tree_dir else None,
                    "matrix": str(matrix) if matrix else None,
                    "partitions": str(partitions) if partitions else None,
                    "model_expr": model_expr,
                    "scf_quartets": scf_quartets if cf not in ("gcf", "qcf") else None,
                    "lpp": lpp,
                    "prefix": _prefix,
                    "output_dir": str(output_dir),
                    "threads": threads,
                    "overwrite": overwrite,
                    "dry_run": dry_run,
                    "iqtree_path": str(iqtree_path) if iqtree_path else None,
                    "wastral_path": str(wastral_path) if wastral_path else None,
                    "iqtree_exe": None,
                    "wastral_exe": None,
                },
                "key_results": {"cf_type": cf, "prefix": prefix or ""},
                "error": error_msg,
                "data": {
                    "input_mode": "--tree" if tree else "--tree-dir" if tree_dir else "--matrix",
                    "input": {},
                    "cmd": [],
                    "tool_stderr": "",
                    "skipped": [],
                    "warnings": [],
                },
            }
            (output_dir / "result.json").write_text(_json.dumps(result, indent=2))
        except Exception:
            pass
        _fail(error_msg, exit_code)

    if dry_run:
        if not quiet:
            cf_type = payload["key_results"]["cf_type"]
            click.echo(f"Dry run: --cf {cf_type} would be executed.")
            click.echo(" ".join(payload["data"]["cmd"]))
        return

    result_path = output_dir / "result.json"
    result_path.parent.mkdir(parents=True, exist_ok=True)
    with open(result_path, "w") as fh:
        json.dump(payload, fh, indent=2)

    if payload["status"] == "error":
        _fail(payload.get("error", "CF computation failed"), 2)

    if not quiet:
        prefix_val = payload["key_results"]["prefix"]
        if cf in ("gcf", "scf", "scfl", "gcf+scf"):
            click.echo(f"CF results saved to {output_dir}/{prefix_val}.cf.*")
        else:
            click.echo(f"qCF tree saved to {output_dir}/{prefix_val}.cf.tree")


@tree.group(
    "bi",
    cls=_BiGroup,
    help=(
        "Bayesian phylogenetic inference with PhyloBayes-MPI.\n\n"
        "Subcommands:\n"
        "  pb         Run MCMC chains (pb_mpi).\n"
        "  bpcomp     Topology convergence analysis (bpcomp).\n"
        "  tracecomp  Parameter convergence analysis (tracecomp).\n"
        "  readpb     Posterior analysis and predictive checks (readpb_mpi).\n"
    ),
)
def bi_group() -> None:
    """Bayesian phylogenetic inference with PhyloBayes-MPI."""


@bi_group.command(
    "pb",
    cls=_CFCommand,
    help=(
        "Run MCMC chains with PhyloBayes-MPI (pb_mpi).\n\n"
        "\n\n"
        "  Run N independent MCMC chains in parallel (mpirun + pb_mpi),\n"
        "  monitor convergence in real time (bpcomp + tracecomp), and\n"
        "  produce a consensus tree when chains are stopped.\n\n"
        "\n\n"
        "  Chains run until Ctrl+C (safe soft-stop) or --nsamples cycles\n"
        "  are reached. The command stays alive, showing a live progress\n"
        "  bar and periodic convergence statistics.\n\n"
        "\n\n"
        "  Examples:\n\n"
        "    Default: 3 chains, CAT-GTR, run forever\n"
        "      phyloai tree bi pb --matrix concat/matrix.phy\n\n"
        "    Homogeneous LG+G4 model, auto-stop after 11000 cycles\n"
        "      phyloai tree bi pb --matrix concat/matrix.phy --model lg --mixture 1 --nsamples 11000\n\n"
        "    Add extra chains to an existing run\n"
        "      phyloai tree bi pb --matrix concat/matrix.phy --chain-names chain4,chain5 -o runs/tree/bi\n\n"
        "    Resume all chains in a directory\n"
        "      phyloai tree bi pb -o runs/tree/bi --resume\n\n"
        "    Resume selected chains only\n"
        "      phyloai tree bi pb -o runs/tree/bi --resume chain1,chain3\n\n"
        "    Resume and extend to a higher nsamples target\n"
        "      phyloai tree bi pb -o runs/tree/bi --resume --nsamples 10000\n\n"
        "    Resume and run forever (override previous nsamples)\n"
        "      phyloai tree bi pb -o runs/tree/bi --resume --nsamples -1\n\n"
        "\n\n"
        "  After the run, determine burn-in and summarise results:\n"
        "    bpcomp -x 5000 chains/chain1 chains/chain2 chains/chain3\n"
        "    tracecomp -x 5000 chains/chain1.trace chains/chain2.trace chains/chain3.trace\n\n"
        "\n\n"
        "  Mutually exclusive:\n"
        "    --start-tree / --fix-tree\n"
        "    --overwrite / --resume\n"
        "    --chains/--chain-prefix vs --chain-names (names override count)\n\n"
    ),
)
@click.option(
    "--matrix",
    type=click.Path(dir_okay=False, path_type=Path),
    default=None,
    help="Input alignment (PHYLIP sequential or FASTA). FASTA is auto-converted to PHYLIP for pb_mpi. Required unless --resume is used.",
)
@click.option(
    "--output-dir", "-o",
    type=click.Path(file_okay=False, path_type=Path),
    default=Path("runs/tree/bi"),
    show_default=True,
    help="Output directory for chains/, convergence/, run_state.json, and result.json.",
)
@click.option(
    "--model",
    type=click.Choice(["gtr", "poisson", "lg", "wag", "jtt", "mtrev", "mtzoa", "mtart"]),
    default="gtr",
    show_default=True,
    help="Relative exchangeability matrix. gtr=general time reversible, lg/wag/jtt=empirical, poisson=equal rates, mtrev/mtzoa/mtart=mitochondrial.",
)
@click.option(
    "--mixture",
    type=str,
    default="auto",
    show_default=True,
    help="Profile mixture model. 'auto'=CAT Dirichlet process (recommended). '1'=homogeneous single matrix (e.g. LG+G4). Integer N=N-component fixed mixture (e.g. '20'=CAT20).",
)
@click.option(
    "--gamma-cats",
    type=click.IntRange(1, None),
    default=4,
    show_default=True,
    help="Number of discrete Gamma rate categories (pb_mpi -dgam). Typical: 4.",
)
@click.option(
    "--start-tree",
    type=click.Path(dir_okay=False, path_type=Path),
    default=None,
    help="Starting tree (Newick). Topology is free to change during MCMC. Mutually exclusive with --fix-tree.",
)
@click.option(
    "--fix-tree",
    type=click.Path(dir_okay=False, path_type=Path),
    default=None,
    help="Fix topology throughout MCMC (Newick). Only branch lengths and other parameters are sampled. Tree must be bifurcating. Mutually exclusive with --start-tree.",
)
@click.option(
    "--chains",
    type=click.IntRange(1, None),
    default=3,
    show_default=True,
    help="Number of independent MCMC chains. Auto-named as --chain-prefix + 1..N. Overridden by --chain-names.",
)
@click.option(
    "--chain-prefix",
    type=str,
    default="chain",
    show_default=True,
    help="Prefix for auto-generated chain names (chain1, chain2, ...). Ignored when --chain-names is given.",
)
@click.option(
    "--chain-names",
    type=str,
    default=None,
    help="Comma-separated explicit chain names (e.g. 'chain4,chain5'). Overrides --chains and --chain-prefix. Use to add new chains to an existing run.",
)
@click.option(
    "--threads", "-t",
    type=click.IntRange(2, None),
    default=4,
    show_default=True,
    help="MPI processes per chain (mpirun -np). Minimum 2: 1 master + N-1 slaves.",
)
@click.option(
    "--sample-freq",
    type=click.IntRange(1, None),
    default=1,
    show_default=True,
    help="Save one MCMC point every N cycles (pb_mpi -x <every>). Lower values give denser sampling but larger files.",
)
@click.option(
    "--nsamples",
    type=int,
    default=None,
    show_default=False,
    help="Total MCMC cycles per chain after which pb_mpi stops (pb_mpi -x <until>). -1 = run forever (Ctrl+C to stop). Default (not set) = -1 for fresh runs, or use the stored value on --resume. Note: with --sample-freq=N, the number of saved points is --nsamples / N. On --resume, a new --nsamples value overrides the stored target; use this to extend a completed run (e.g. --resume --nsamples 10000).",
)
@click.option(
    "--monitor-freq",
    type=click.IntRange(1, None),
    default=100,
    show_default=True,
    help="Run bpcomp + tracecomp convergence diagnostics every N new saved samples (minimum across all chains).",
)
@click.option(
    "--burnin-frac",
    type=click.FloatRange(0.0, None),
    default=0.5,
    show_default=True,
    help="Fraction of saved samples discarded for convergence monitoring only (0.0 <= x < 1.0). NOT passed to pb_mpi; use bpcomp/tracecomp after the run to choose a final burn-in.",
)
@click.option(
    "--poll-interval",
    type=click.IntRange(1, None),
    default=60,
    show_default=True,
    help="Seconds between .trace file reads for progress update and convergence triggers. Larger values reduce I/O on network filesystems.",
)
@click.option(
    "--pb-path",
    type=click.Path(file_okay=False, path_type=Path),
    default=None,
    help="Directory containing PhyloBayes-MPI tools (pb_mpi, bpcomp, tracecomp). Overrides PATH lookup.",
)
@click.option(
    "--resume",
    default=None,
    is_flag=False,
    flag_value="__ALL__",
    help="Resume chains from their .chain state. Bare --resume resumes all chains in run_state.json. --resume chain1,chain3 resumes only those. Mutually exclusive with --overwrite. Use --nsamples to extend to a new target (e.g. --resume --nsamples 10000 to continue from 5000).",
)
@click.option(
    "--overwrite",
    is_flag=True,
    default=False,
    help="Delete and recreate the output directory before starting. Mutually exclusive with --resume.",
)
@click.option(
    "--dry-run",
    is_flag=True,
    default=False,
    help="Print the mpirun + pb_mpi commands that would be run, then exit. No processes are started.",
)
@click.option(
    "--quiet", "-q",
    is_flag=True,
    default=False,
    help="Suppress all terminal output except errors.",
)
def bi_pb_command(**kwargs) -> None:
    from phyloai.tree.bi import run_bi_pb

    try:
        payload = run_bi_pb(**kwargs)
    except FileNotFoundError as exc:
        _write_error_result(kwargs, str(exc), 3)
        _fail(str(exc), 3)
    except ValueError as exc:
        _write_error_result(kwargs, str(exc), 1)
        _fail(str(exc), 1)
    if kwargs.get("dry_run"):
        if not kwargs.get("quiet"):
            for cmd in payload["data"]["chain_cmds"].values():
                click.echo(" ".join(cmd))
        return
    result_path = kwargs["output_dir"] / "result.json"
    result_path.parent.mkdir(parents=True, exist_ok=True)

    # Add convergence_render.txt to output_files if it exists
    conv_path = kwargs["output_dir"] / "convergence" / "convergence_render.txt"
    if conv_path.exists():
        payload.setdefault("data", {}).setdefault("output_files", {})["convergence_render"] = {
            "path": str(conv_path), "description": "Human-readable convergence diagnostic summary"
        }

    if result_path.exists():
        from datetime import datetime
        ts = datetime.now().strftime("%Y%m%d_%H%M%S")
        backup_path = result_path.with_name(f"result_{ts}.json")
        shutil.copy2(str(result_path), str(backup_path))
    with open(result_path, "w") as fh:
        json.dump(payload, fh, indent=2)
    if payload["status"] == "error":
        _fail(payload.get("error", "pb_mpi execution failed"), 2)
    if not kwargs.get("quiet"):
        click.echo(f"Results saved to {result_path}", err=True)


# ---------------------------------------------------------------------------
# tree bi bpcomp
# ---------------------------------------------------------------------------


@bi_group.command(
    "bpcomp",
    cls=_CFCommand,
    help=(
        "Run bpcomp for final topology convergence analysis.\n\n"
        "Runs bpcomp once with a user-specified integer burn-in on a "
        "completed or running chains directory. Produces a final "
        "consensus tree and bipartition statistics.\n\n"
        "Examples:\n\n"
        "  phyloai tree bi bpcomp --chain-dir runs/tree/bi/chains --burnin 5000\n\n"
        "  phyloai tree bi bpcomp --chain-dir runs/tree/bi/chains --chain-names chain1,chain3 --burnin 2000 --sample-freq 5 --until 3000 --cutoff 0.75\n\n"
    ),
)
@click.option(
    "--chain-dir",
    type=click.Path(file_okay=False, path_type=Path),
    required=True,
    help="Directory containing chain files (<chain>.chain, <chain>.treelist, etc.).",
)
@click.option(
    "--chain-names",
    type=str,
    default="all",
    show_default=True,
    help="Comma-separated chain names to include (e.g. 'chain1,chain2'). 'all' = all chains in --chain-dir.",
)
@click.option(
    "--output-dir", "-o",
    type=click.Path(file_okay=False, path_type=Path),
    default=Path("runs/tree/bi/bpcomp"),
    show_default=True,
    help="Output directory for result files and result.json.",
)
@click.option(
    "--burnin",
    type=click.IntRange(0),
    default=0,
    show_default=True,
    help="Number of saved samples to discard as burn-in. 0 = no burn-in.",
)
@click.option(
    "--sample-freq",
    type=click.IntRange(1),
    default=1,
    show_default=True,
    help="Sub-sampling frequency: take one tree every N saved samples after burn-in.",
)
@click.option(
    "--until",
    type=str,
    default="all",
    show_default=True,
    help="Stop at this sample index. 'all' = use the entire chain. Integer = stop at that sample index.",
)
@click.option(
    "--cutoff",
    type=click.FloatRange(0.0, 1.0),
    default=0.5,
    show_default=True,
    help="Majority-rule consensus cutoff: nodes below this are collapsed.",
)
@click.option(
    "--pb-path",
    type=click.Path(file_okay=False, path_type=Path),
    default=None,
    help="Directory containing PhyloBayes tools. Overrides PATH lookup.",
)
@click.option("--overwrite", is_flag=True, default=False, help="Delete and recreate the output directory.")
@click.option("--dry-run", is_flag=True, default=False, help="Print the bpcomp command without executing.")
@click.option("--quiet", "-q", is_flag=True, default=False, help="Suppress non-error terminal output.")
def bi_bpcomp_command(**kwargs) -> None:
    from phyloai.tree.bi_bpcomp import run_bi_bpcomp

    try:
        payload = run_bi_bpcomp(**kwargs)
    except FileNotFoundError as exc:
        _write_bi_error_result(kwargs, str(exc), 3, "bpcomp", _BPCOMP_FLAG_MAP)
        _fail(str(exc), 3)
    except ValueError as exc:
        if "already exists and is non-empty" not in str(exc):
            _write_bi_error_result(kwargs, str(exc), 1, "bpcomp", _BPCOMP_FLAG_MAP)
        _fail(str(exc), 1)
    if kwargs.get("dry_run"):
        if not kwargs.get("quiet"):
            click.echo(" ".join(payload["data"]["cmd"]))
        return
    result_path = kwargs["output_dir"] / "result.json"
    if result_path.exists():
        from datetime import datetime
        ts = datetime.now().strftime("%Y%m%d_%H%M%S")
        shutil.copy2(str(result_path), str(result_path.with_name(f"result_{ts}.json")))
    with open(result_path, "w") as fh:
        json.dump(payload, fh, indent=2)
    if payload["status"] == "error":
        _fail(payload.get("error", "bpcomp failed"), 2)
    if not kwargs.get("quiet"):
        click.echo(f"Results saved to {result_path}", err=True)


# ---------------------------------------------------------------------------
# tree bi tracecomp
# ---------------------------------------------------------------------------


@bi_group.command(
    "tracecomp",
    cls=_CFCommand,
    help=(
        "Run tracecomp for final parameter convergence analysis.\n\n"
        "Runs tracecomp once with a user-specified integer burn-in to "
        "assess parameter convergence across chains.\n\n"
        "Examples:\n\n"
        "  phyloai tree bi tracecomp --chain-dir runs/tree/bi/chains --burnin 5000\n\n"
        "  phyloai tree bi tracecomp --chain-dir runs/tree/bi/chains --chain-names chain1,chain2 --burnin 2000\n\n"
    ),
)
@click.option(
    "--chain-dir",
    type=click.Path(file_okay=False, path_type=Path),
    required=True,
    help="Directory containing chain .trace files.",
)
@click.option(
    "--chain-names",
    type=str,
    default="all",
    show_default=True,
    help="Comma-separated chain names (e.g. 'chain1,chain2'). 'all' = all chains with .trace files in --chain-dir.",
)
@click.option(
    "--output-dir", "-o",
    type=click.Path(file_okay=False, path_type=Path),
    default=Path("runs/tree/bi/tracecomp"),
    show_default=True,
    help="Output directory for tracecomp.contdiff and result.json.",
)
@click.option(
    "--burnin",
    type=click.IntRange(0),
    default=0,
    show_default=True,
    help="Number of saved samples to discard as burn-in. 0 = no burn-in.",
)
@click.option(
    "--pb-path",
    type=click.Path(file_okay=False, path_type=Path),
    default=None,
    help="Directory containing PhyloBayes tools. Overrides PATH lookup.",
)
@click.option("--overwrite", is_flag=True, default=False, help="Delete and recreate the output directory.")
@click.option("--dry-run", is_flag=True, default=False, help="Print the tracecomp command without executing.")
@click.option("--quiet", "-q", is_flag=True, default=False, help="Suppress non-error terminal output.")
def bi_tracecomp_command(**kwargs) -> None:
    from phyloai.tree.bi_tracecomp import run_bi_tracecomp

    try:
        payload = run_bi_tracecomp(**kwargs)
    except FileNotFoundError as exc:
        _write_bi_error_result(kwargs, str(exc), 3, "tracecomp", _TRACECOMP_FLAG_MAP)
        _fail(str(exc), 3)
    except ValueError as exc:
        if "already exists and is non-empty" not in str(exc):
            _write_bi_error_result(kwargs, str(exc), 1, "tracecomp", _TRACECOMP_FLAG_MAP)
        _fail(str(exc), 1)
    if kwargs.get("dry_run"):
        if not kwargs.get("quiet"):
            click.echo(" ".join(payload["data"]["cmd"]))
        return
    result_path = kwargs["output_dir"] / "result.json"
    if result_path.exists():
        from datetime import datetime
        ts = datetime.now().strftime("%Y%m%d_%H%M%S")
        shutil.copy2(str(result_path), str(result_path.with_name(f"result_{ts}.json")))
    with open(result_path, "w") as fh:
        json.dump(payload, fh, indent=2)
    if payload["status"] == "error":
        _fail(payload.get("error", "tracecomp failed"), 2)
    if not kwargs.get("quiet"):
        click.echo(f"Results saved to {result_path}", err=True)


# ---------------------------------------------------------------------------
# tree bi readpb
# ---------------------------------------------------------------------------


@bi_group.command(
    "readpb",
    cls=_CFCommand,
    help=(
        "Run readpb_mpi for posterior analysis on a single chain.\n\n"
        "Supports multiple analysis modes (-rr, -ss, -r, -sitelogl, etc.). "
        "Automatically converts rr output to IQ-TREE exchangeabilities format "
        "and ss output to IQ-TREE site frequencies format.\n\n"
        "Modes:\n\n"
        "  rr             Posterior mean relative exchangeabilities (-> exchangeabilities).\n"
        "  ss             Posterior mean site-specific state frequencies (-> sitefreq).\n"
        "  r              Posterior mean rates across sites.\n"
        "  sitelogl       Site-specific marginal log-likelihoods (wAIC / LOO).\n"
        "  ppred          Simulate data replicates from the posterior predictive distribution.\n"
        "  div            Posterior predictive diversity test (PPA-DIV).\n"
        "  sitecomp       Posterior predictive compositional heterogeneity test (PPA-VAR).\n"
        "  siteconvprob   Posterior predictive convergence probability test (PPA-CONV).\n"
        "  comp           Posterior predictive compositional homogeneity test.\n"
        "  allppred       All four predictive tests at once\n"
        "                 (mutually exclusive with div/sitecomp/siteconvprob/comp).\n\n"
        "  rr/ss convert to IQ-TREE-compatible exchangeabilities/sitefreq files.\n\n"
        "Examples:\n\n"
        "  phyloai tree bi readpb --chain runs/tree/bi/chains/chain1 --mode rr,ss --burnin 2000 --sample-freq 5\n\n"
        "  phyloai tree bi readpb --chain runs/tree/bi/chains/chain1 --mode allppred --burnin 5000\n\n"
    ),
)
@click.option(
    "--chain",
    type=click.Path(dir_okay=False, path_type=Path),
    required=True,
    help="Path to chain file without extension (e.g. runs/tree/bi/chains/chain1).",
)
@click.option(
    "--mode",
    type=str,
    required=True,
    help="Comma-separated list of analysis modes: rr, ss, r, sitelogl, ppred, div, sitecomp, siteconvprob, comp, allppred.",
)
@click.option(
    "--output-dir", "-o",
    type=click.Path(file_okay=False, path_type=Path),
    default=Path("runs/tree/bi/readpb"),
    show_default=True,
    help="Output directory for readpb outputs and result.json.",
)
@click.option(
    "--burnin",
    type=click.IntRange(0),
    default=0,
    show_default=True,
    help="Number of saved samples to discard as burn-in.",
)
@click.option(
    "--sample-freq",
    type=click.IntRange(1),
    default=1,
    show_default=True,
    help="Sub-sampling frequency after burn-in.",
)
@click.option(
    "--until",
    type=str,
    default="all",
    show_default=True,
    help="'all' = to end of chain. Integer = stop at that saved sample index.",
)
@click.option(
    "--threads", "-t",
    type=click.IntRange(2),
    default=4,
    show_default=True,
    help="MPI processes for readpb_mpi.",
)
@click.option(
    "--pb-path",
    type=click.Path(file_okay=False, path_type=Path),
    default=None,
    help="Directory containing PhyloBayes tools. Overrides PATH lookup.",
)
@click.option("--overwrite", is_flag=True, default=False, help="Delete and recreate --output-dir.")
@click.option("--dry-run", is_flag=True, default=False, help="Print commands without executing.")
@click.option("--quiet", "-q", is_flag=True, default=False, help="Suppress non-error terminal output.")
def bi_readpb_command(**kwargs) -> None:
    from phyloai.tree.bi_readpb import run_bi_readpb

    try:
        payload = run_bi_readpb(**kwargs)
    except FileNotFoundError as exc:
        _write_bi_error_result(kwargs, str(exc), 3, "readpb", _READPB_FLAG_MAP)
        _fail(str(exc), 3)
    except ValueError as exc:
        if "already exists and is non-empty" not in str(exc):
            _write_bi_error_result(kwargs, str(exc), 1, "readpb", _READPB_FLAG_MAP)
        _fail(str(exc), 1)
    if kwargs.get("dry_run"):
        if not kwargs.get("quiet"):
            for mode_name, cmd in payload["data"]["cmds"].items():
                click.echo(" ".join(cmd))
        return
    result_path = kwargs["output_dir"] / "result.json"
    if result_path.exists():
        from datetime import datetime
        ts = datetime.now().strftime("%Y%m%d_%H%M%S")
        shutil.copy2(str(result_path), str(result_path.with_name(f"result_{ts}.json")))
    with open(result_path, "w") as fh:
        json.dump(payload, fh, indent=2)
    if payload["status"] == "error":
        _fail(payload.get("error", "readpb failed"), 2)
    if not kwargs.get("quiet"):
        click.echo(f"Results saved to {result_path}", err=True)
