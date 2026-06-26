"""IQ-TREE3 hessian computation for MCMCtree approximate likelihood dating."""
from __future__ import annotations

import re
import shlex
import subprocess
import time
from pathlib import Path
from typing import Any

from phyloai.core.formats import FormatConverter
from phyloai.core.iqtree import _resolve_iqtree_path, _detect_iqtree_version
from phyloai.core.sequence_normalization import detect_seq_type

# ponytail: prefix hardcoded to "iqtree" — custom prefix breaks
# mcmc's ability to find hessian output files. Remove --prefix from
# CLI instead of making every consumer configurable.
HESSIAN_PREFIX = "iqtree"
HESSIAN_OUTPUT_FILES = (
    f"{HESSIAN_PREFIX}.dummy.phy",
    f"{HESSIAN_PREFIX}.rooted.nwk",
    f"{HESSIAN_PREFIX}.mcmctree.hessian",
)

def detect_seqtype_from_alignment(matrix: Path) -> str:
    converter = FormatConverter()
    alignment = converter.read(matrix)
    return detect_seq_type([str(record.seq) for record in alignment])


def count_partitions(partition_file: Path) -> int:
    text = partition_file.read_text(errors="ignore")
    if "charset" in text.lower():
        return len(re.findall(r"(?i)\bcharset\b", text))
    lines = [
        l.strip() for l in text.splitlines()
        if l.strip() and not l.strip().startswith("#")
    ]
    return len(lines)


def validate_root_age(newick: str) -> bool:
    return bool(re.search(r"\)['\"]?[^;,()]*[<>][^;,()]*['\"]?\s*;", newick))


def build_iqtree_dating_cmd(
    *,
    iqtree_path: Path,
    matrix: Path,
    rooted_tree: Path,
    seq_type: str,
    model_expr: str | None,
    partitions: Path | None,
    n_partitions: int,
    threads: int,
    tool_args: str | None,
) -> list[str]:
    cmd: list[str] = [str(iqtree_path), "-s", str(matrix), "-te", str(rooted_tree),
                      "--dating", "mcmctree", "--prefix", HESSIAN_PREFIX, "-T", str(threads)]

    if partitions is None:
        model = model_expr or ("LG+F+G4" if seq_type == "AA" else "GTR+G4")
        cmd += ["-m", model]
    else:
        cmd += ["-m", "MF", "-Q", str(partitions)]
        if seq_type == "AA":
            cmd += ["--mset", "LG", "-mfreq", "F", "-mrate", "G"]
        else:
            cmd += ["--mset", "GTR", "-mrate", "G"]
        if n_partitions >= 10:
            cmd += ["--merge", "--rclusterf", "10"]

    if tool_args:
        cmd += shlex.split(tool_args)

    return cmd



def _validate_hessian_inputs(
    *,
    matrix: Path,
    rooted_tree: Path,
    seq_type: str,
    model_expr: str | None,
    partitions: Path | None,
    threads: int,
    overwrite: bool,
    resume: bool,
    tool_args: str | None,
) -> list[str]:
    errors: list[str] = []
    import os as _os

    if not matrix.exists():
        errors.append(f"--matrix does not exist: {matrix}")
    elif not matrix.is_file():
        errors.append(f"--matrix is not a regular file: {matrix}")
    elif not _os.access(str(matrix), _os.R_OK):
        errors.append(f"--matrix is not readable: {matrix}")

    if not rooted_tree.exists():
        errors.append(f"--rooted-tree does not exist: {rooted_tree}")
    elif not rooted_tree.is_file():
        errors.append(f"--rooted-tree is not a regular file: {rooted_tree}")
    else:
        content = rooted_tree.read_text(errors="ignore")
        if not validate_root_age(content):
            errors.append(
                "--rooted-tree is missing a root age constraint. "
                "The outermost node must have a calibration label such as "
                "'<4.2' or '>3.1<4.2' (units: 100 Mya). "
                "Example: (A,(B,C))'<4.2';"
            )

    if model_expr and partitions:
        errors.append("--model-expr and --partitions are mutually exclusive.")

    if partitions:
        if not partitions.exists():
            errors.append(f"--partitions does not exist: {partitions}")
        elif not partitions.is_file():
            errors.append(f"--partitions is not a regular file: {partitions}")

    if threads < 1:
        errors.append(f"--threads must be >= 1, got {threads}")

    if overwrite and resume:
        errors.append("--overwrite and --resume are mutually exclusive.")

    if tool_args:
        blocked = {"-s", "--dating", "-te", "--prefix"}
        tokens = shlex.split(tool_args)
        for tok in tokens:
            if tok in blocked:
                errors.append(
                    f"--tool-args contains blocked flag '{tok}' "
                    f"(managed by PhyloAI). The hessian step must emit "
                    f"iqtree.dummy.phy, iqtree.rooted.nwk, and "
                    f"iqtree.mcmctree.hessian under the output "
                    f"directory using the supplied calibrated tree."
                )

    return errors


def run_hessian(
    *,
    matrix: Path,
    rooted_tree: Path,
    seq_type: str = "auto",
    model_expr: str | None = None,
    partitions: Path | None = None,
    output_dir: Path,
    threads: int = 4,
    iqtree_path: str | None = None,
    tool_args: str | None = None,
    overwrite: bool = False,
    resume: bool = False,
    dry_run: bool = False,
    quiet: bool = False,
    stream_output: bool = True,
) -> dict[str, Any]:
    t0 = time.time()

    errors = _validate_hessian_inputs(
        matrix=matrix, rooted_tree=rooted_tree,
        seq_type=seq_type,
        model_expr=model_expr, partitions=partitions,
        threads=threads, overwrite=overwrite, resume=resume,
        tool_args=tool_args,
    )
    if errors:
        return {
            "status": "error",
            "command": "",
            "wall_time": 0.0,
            "tool_versions": {},
            "params": {},
            "key_results": {},
            "error": errors[0],
            "error_category": "input",
            "data": {"cmd": [], "tool_stderr": "", "warnings": errors},
        }

    matrix = matrix.resolve()
    rooted_tree = rooted_tree.resolve()
    if partitions:
        partitions = partitions.resolve()

    if seq_type == "auto":
        try:
            seq_type = detect_seqtype_from_alignment(matrix)
        except Exception as e:
            return {
                "status": "error",
                "command": "",
                "wall_time": 0.0,
                "tool_versions": {},
                "params": {},
                "key_results": {},
                "error": f"Cannot read --matrix as FASTA/PHYLIP/NEXUS: {e}",
                "error_category": "input",
                "data": {"cmd": [], "tool_stderr": "", "warnings": [str(e)]},
            }

    n_partitions = 1  # unpartitioned = single ndata block
    if partitions:
        n_partitions = count_partitions(partitions)

    try:
        iqtree_exe = _resolve_iqtree_path(iqtree_path, dry_run)
    except (ValueError, FileNotFoundError) as e:
        return {
            "status": "error",
            "command": "",
            "wall_time": 0.0,
            "tool_versions": {},
            "params": {},
            "key_results": {},
            "error": str(e),
            "error_category": "env",
            "data": {"cmd": [], "tool_stderr": "", "warnings": [str(e)]},
        }
    tool_versions_iqtree = (
        _detect_iqtree_version(iqtree_exe) if not dry_run else {"iqtree3": "dry-run"}
    )

    cmd = build_iqtree_dating_cmd(
        iqtree_path=Path(iqtree_exe),
        matrix=matrix,
        rooted_tree=rooted_tree,
        seq_type=seq_type,
        model_expr=model_expr,
        partitions=partitions,
        n_partitions=n_partitions,
        threads=threads,
        tool_args=tool_args,
    )

    if dry_run:
        return {
            "status": "success",
            "command": " ".join(cmd),
            "wall_time": 0.0,
            "tool_versions": tool_versions_iqtree,
            "params": {"seq_type": seq_type, "n_partitions": n_partitions},
            "key_results": {},
            "error": None,
            "data": {"cmd": cmd, "tool_stderr": "", "warnings": []},
        }

    output_dir.mkdir(parents=True, exist_ok=True)

    proc = subprocess.run(
        cmd,
        cwd=output_dir,
        stdout=None if stream_output else subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )

    warnings: list[str] = []
    missing = [f for f in HESSIAN_OUTPUT_FILES if not (output_dir / f).exists()]
    if not missing:
        empty = [f for f in HESSIAN_OUTPUT_FILES
                 if (output_dir / f).stat().st_size == 0]
        if empty:
            missing = empty
            warnings.append(
                f"IQ-TREE produced empty output file(s): {empty}. "
                "Possible crash mid-write."
            )
    iqtree_report = output_dir / f"{HESSIAN_PREFIX}.iqtree"
    if proc.returncode == 0 and iqtree_report.exists():
        report_text = iqtree_report.read_text(errors="ignore")
        if "Total CPU time used" not in report_text:
            warnings.append(
                f"{iqtree_report.name} has no 'Total CPU time used' marker — "
                "IQ-TREE may have been interrupted before completion."
            )

    if missing or proc.returncode != 0:
        return {
            "status": "error",
            "command": " ".join(cmd),
            "wall_time": time.time() - t0,
            "tool_versions": tool_versions_iqtree,
            "params": {"seq_type": seq_type, "n_partitions": n_partitions},
            "key_results": {},
            "error": f"IQ-TREE failed (returncode={proc.returncode}). Missing: {missing}",
            "error_category": "tool",
            "data": {"cmd": cmd, "tool_stderr": getattr(proc, "stderr", ""), "warnings": warnings},
        }

    wall = time.time() - t0
    return {
        "status": "success",
        "command": " ".join(cmd),
        "wall_time": wall,
        "tool_versions": tool_versions_iqtree,
        "params": {
            "seq_type": seq_type,
            "n_partitions": n_partitions,
            "model_expr": model_expr,
            "partitions": str(partitions) if partitions else None,
            "threads": threads,
        },
        "key_results": {
            "hessian_file": str(output_dir / f"{HESSIAN_PREFIX}.mcmctree.hessian"),
        },
        "error": None,
        "data": {
            "cmd": cmd,
            "tool_stderr": getattr(proc, "stderr", ""),
            "warnings": warnings,
            "output_files": {f: str(output_dir / f) for f in HESSIAN_OUTPUT_FILES},
        },
    }
