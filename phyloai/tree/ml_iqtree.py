"""Maximum-likelihood tree inference with IQ-TREE3."""

from __future__ import annotations

import os
import re as _re
import shlex
import shutil
import subprocess
import tempfile
import time as _time
import warnings
from pathlib import Path
from typing import Any

from Bio import Phylo, SeqIO

from phyloai.core.env import ToolEnv
from phyloai.core.iqtree import (
    _resolve_iqtree_path,
    _detect_iqtree_version,
)
from phyloai.core.schema import COMMON_ALIGNMENT_EXTENSIONS
from phyloai.core.sequence_normalization import detect_seq_type

# ===================================================================
# Shared constants
# ===================================================================

IQTREE_COMPATIBLE_EXTENSIONS = frozenset({
    ".fa", ".fas", ".fasta", ".faa", ".fna",
    ".phy", ".phylip",
    ".nex", ".nxs", ".nexus",
    ".aln",
})

AA_STANDARD_MODELS = frozenset({
    "LG", "Poisson", "cpREV", "mtREV", "Dayhoff", "mtMAM", "JTT", "WAG",
    "mtART", "mtZOA", "VT", "rtREV", "DCMut", "PMB", "HIVb", "HIVw",
    "JTTDCMut", "FLU", "Blosum62", "GTR20", "mtMet", "mtVer", "mtInv",
    "FLAVI", "Q.LG", "Q.pfam", "Q.pfam_gb", "Q.bird", "Q.mammal",
    "Q.insect", "Q.plant", "Q.yeast",
})

NT_STANDARD_MODELS = frozenset({
    "GTR", "HKY", "JC", "F81", "K2P", "K3P", "K81uf", "TN", "TNef",
    "TIM", "TIMef", "TVM", "TVMef", "SYM",
})

AA_MIXTURE_MODELS = frozenset(
    {f"C{i}" for i in range(10, 61, 10)}
    | {"EX2", "EX3", "EHO", "UL2", "UL3", "EX_EHO", "LG4M", "LG4X"}
)

_IQTREE_MANAGED_LONG_FLAGS = frozenset({
    "--ufboot", "--alrt", "--bnni", "--fast", "--merge",
    "--rclusterf", "--rcluster-max", "--mset", "--msub",
    "--prefix", "--rate", "--qmax", "--seqtype",
    "--redo", "--redo-tree", "--undo",
})

_IQTREE_MANAGED_SHORT_FLAGS = frozenset({
    "-s", "-m", "-p", "-T", "-B", "-ft", "-g", "-o",
    "-q", "-Q", "-S", "-wslr",
})

_IQTREE_BLOCKED_FLAGS = frozenset({"-s"})
_IQTREE_BLOCKED_FLAGS_BATCH = frozenset({"-s", "--prefix"})
_IQTREE_BLOCKED_IO_CHARS = frozenset({"<", ">", "|"})

_IQTREE_FLAG_ALIASES: dict[str, frozenset[str]] = {
    "-B": frozenset({"-B", "--ufboot"}),
    "--ufboot": frozenset({"-B", "--ufboot"}),
}


# ===================================================================
# Input scanning
# ===================================================================

def _scan_input_iqtree(msa_dir: Path) -> tuple[list[Path], list[dict[str, str]]]:
    if not msa_dir.exists():
        return [], []

    found: list[Path] = []
    skipped: list[dict[str, str]] = []

    for entry in sorted(msa_dir.iterdir()):
        if entry.is_dir():
            skipped.append({"path": str(entry), "reason": "directory"})
            continue
        if not entry.is_file():
            skipped.append({"path": str(entry), "reason": "not a regular file"})
            continue
        if entry.stat().st_size == 0:
            skipped.append({"path": str(entry), "reason": "empty file"})
            continue

        ext = entry.suffix.lower()
        if ext in IQTREE_COMPATIBLE_EXTENSIONS:
            found.append(entry)
        elif ext in set(COMMON_ALIGNMENT_EXTENSIONS):
            skipped.append({"path": str(entry), "reason": f"unrecognized extension: {ext}"})
        else:
            skipped.append({"path": str(entry), "reason": f"unrecognized extension: {ext}"})

    return found, skipped


# ===================================================================
# Format-aware sequence type validation for IQ-TREE inputs
# ===================================================================

def _detect_file_format(ext: str) -> str:
    """Map file extension to Bio.SeqIO format name."""
    ext = ext.lower()
    if ext in {".nex", ".nxs", ".nexus"}:
        return "nexus"
    if ext in {".aln"}:
        return "clustal"
    if ext in {".phy", ".phylip"}:
        return "phylip-relaxed"
    return "fasta"


def _validate_seq_types_iqtree(
    files: list[Path],
    *,
    declared_type: str | None,
) -> tuple[str | None, list[dict[str, Any]]]:
    if not files:
        return (declared_type or "AA"), []

    all_types: dict[str, str] = {}
    offending: list[dict[str, Any]] = []

    for f in files:
        try:
            fmt = _detect_file_format(f.suffix)
            with open(str(f)) as fh:
                seqs = [str(r.seq) for r in SeqIO.parse(fh, fmt)]
            if not seqs:
                offending.append({"file": str(f), "reason": "no sequences found"})
                continue
            dt = detect_seq_type(seqs)
            all_types[str(f)] = dt
        except Exception:
            offending.append({"file": str(f), "reason": "failed to parse input file"})
            continue

    if declared_type:
        for f_str, dt in all_types.items():
            if dt != declared_type:
                offending.append({"file": f_str, "expected": declared_type, "detected": dt})
        return declared_type, offending

    type_counts: dict[str, int] = {}
    for dt in all_types.values():
        type_counts[dt] = type_counts.get(dt, 0) + 1

    if not type_counts:
        return None, offending

    if len(type_counts) == 1:
        resolved = next(iter(type_counts))
        return resolved, []

    majority = max(type_counts, key=type_counts.get)
    for f_str, dt in all_types.items():
        if dt != majority:
            offending.append({"file": f_str, "expected": majority, "detected": dt})

    return None, offending


# ===================================================================
# Two-tier --tool-args management
# ===================================================================

def _format_offender(o: dict[str, Any]) -> str:
    if "reason" in o:
        return f"{o['file']}: {o['reason']}"
    return f"{o['file']}: {o['detected']} (expected {o['expected']})"


def _check_managed_flag_conflict(tool_args: str, *, batch_mode: bool = False) -> None:
    """Reject BLOCKED flags and I/O redirects in --tool-args."""
    blocked = _IQTREE_BLOCKED_FLAGS_BATCH if batch_mode else _IQTREE_BLOCKED_FLAGS
    tokens = shlex.split(tool_args)
    for token in tokens:
        if token in blocked:
            raise ValueError(f"Blocked managed flag in --tool-args: {token}")
        if any(c in token for c in _IQTREE_BLOCKED_IO_CHARS):
            raise ValueError(f"Blocked I/O override in --tool-args: {token}")


def _is_flag_overridden(flag: str, tool_tokens: set[str]) -> bool:
    """Check whether a PhyloAI-generated flag is present in --tool-args tokens."""
    if flag in tool_tokens:
        return True
    aliases = _IQTREE_FLAG_ALIASES.get(flag)
    if aliases:
        for alias in aliases:
            if alias in tool_tokens:
                return True
    return False


def _get_tool_arg_value(flag: str, tool_args: str | None) -> str | None:
    """Return value immediately following a flag in --tool-args, if present."""
    if not tool_args:
        return None
    tokens = shlex.split(tool_args)
    try:
        idx = tokens.index(flag)
    except ValueError:
        return None
    if idx + 1 >= len(tokens):
        return None
    return tokens[idx + 1]


# ===================================================================
# Workflow classification
# ===================================================================

def _classify_workflow(
    *,
    modelfinder: str,
    model: str,
    seq_type: str,
    partitions: str | None = None,
    rclusterf: int | None = None,
    rcluster_max: int | None = None,
    guide_tree: str | None = None,
) -> str:
    if modelfinder in ("MF", "MFP"):
        if partitions is not None:
            return f"homogeneous-partition-{modelfinder}-merge"
        return f"homogeneous-no-partition-{modelfinder}"

    if model in AA_MIXTURE_MODELS and seq_type == "AA":
        if guide_tree:
            return "AA-heterogeneous-PMSF"
        return "AA-heterogeneous-direct"

    if model == "MIX+MF" and seq_type == "NT":
        return "NT-heterogeneous"

    if partitions is not None:
        return "homogeneous-partition-none"

    return "homogeneous-no-partition-none"


# ===================================================================
# Model helpers
# ===================================================================

def _build_model_string(
    *,
    model: str,
    state_freq: str,
    rate_heterogeneity: str,
    modelfinder: str,
    pmsf_base_model: str | None = None,
) -> str:
    if modelfinder in ("MF", "MFP"):
        return modelfinder

    base = model
    if pmsf_base_model and model in AA_MIXTURE_MODELS:
        base = f"{pmsf_base_model}+{model}"

    # MIX+MF returns bare model string (no freq/rate appended)
    if model == "MIX+MF":
        return base

    parts = [base]
    if state_freq != "none":
        parts.append(state_freq)
    if rate_heterogeneity != "none":
        parts.append(rate_heterogeneity)
    return "".join(parts)


def _resolve_custom_model_path(model: str) -> str | None:
    """Return an absolute custom-model path when *model* names a regular file."""
    candidate = Path(model).expanduser()
    if not candidate.exists():
        return None
    if not candidate.is_file():
        raise ValueError(f"Custom model path is not a regular file: {candidate}")
    return str(candidate.resolve())


def _resolve_site_freq_file(site_freq_file: str | Path | None) -> str | None:
    """Return an absolute per-site frequency-profile path."""
    if site_freq_file is None:
        return None
    candidate = Path(site_freq_file).expanduser()
    if not candidate.exists():
        raise ValueError(f"--site-freq-file does not exist: {candidate}")
    if not candidate.is_file():
        raise ValueError(f"--site-freq-file is not a regular file: {candidate}")
    return str(candidate.resolve())


def _validate_model(model: str, seq_type: str, modelfinder: str, *, custom_model: bool = False) -> None:
    if custom_model or modelfinder in ("MF", "MFP"):
        return
    if seq_type == "AA":
        valid = AA_STANDARD_MODELS | AA_MIXTURE_MODELS
    elif seq_type == "NT":
        valid = NT_STANDARD_MODELS | {"MIX+MF"}
    else:
        return
    if model not in valid:
        raise ValueError(f"Invalid model for {seq_type}: {model}. Valid: {sorted(valid)}")


def _validate_pmsf_base_model(base_model: str) -> None:
    if base_model not in AA_STANDARD_MODELS:
        raise ValueError(
            f"Invalid PMSF base model: {base_model}. "
            f"Must be a standard AA model: {sorted(AA_STANDARD_MODELS)}"
        )


def _is_heterogeneous_model(model: str, seq_type: str, modelfinder: str) -> bool:
    if modelfinder in ("MF", "MFP"):
        return False
    if seq_type == "AA" and model in AA_MIXTURE_MODELS:
        return True
    if seq_type == "NT" and model == "MIX+MF":
        return True
    return False


# ===================================================================
# Pre-flight validation
# ===================================================================

def _run_validations(
    *,
    batch_mode: bool,
    seq_type: str,
    modelfinder: str,
    model: str,
    partitions: str | None,
    guide_tree: str | None,
    boot: int | None = None,
    alrt: int | None = None,
    bnni: bool = False,
    prefix: str | None = None,
    pmsf_base_model: str | None = None,
    rclusterf: int | None = None,
    rcluster_max: int | None = None,
    qmax: int | None = None,
    state_freq: str = "+F",
    custom_model: bool = False,
    site_freq_file: str | None = None,
    tool_args: str | None = None,
    quiet: bool = False,
) -> None:
    has_raw_site_freq = _is_flag_overridden("-fs", set(shlex.split(tool_args))) if tool_args else False

    if custom_model:
        if batch_mode:
            raise ValueError("Custom model files are only supported with --matrix, not --msa-dir.")
        if seq_type != "AA":
            raise ValueError("Custom model files are only supported for AA data.")
        if modelfinder != "none":
            raise ValueError("Custom model files cannot be used with ModelFinder.")

    if site_freq_file or has_raw_site_freq:
        if state_freq != "none":
            raise ValueError("--site-freq-file and --tool-args -fs require --state-freq none.")
        if site_freq_file:
            if batch_mode:
                raise ValueError("--site-freq-file is only supported with --matrix, not --msa-dir.")
            if seq_type != "AA":
                raise ValueError("--site-freq-file is only supported for AA data.")
            if modelfinder != "none":
                raise ValueError("--site-freq-file cannot be used with ModelFinder.")
            if not custom_model:
                raise ValueError("--site-freq-file requires a custom model file supplied through --model.")

    # Heterogeneous workflows require --matrix
    if _is_heterogeneous_model(model, seq_type, modelfinder) and batch_mode:
        raise ValueError(
            f"Heterogeneous model '{model}' is only supported in --matrix mode, not --msa-dir."
        )

    # Partitions require --matrix
    if partitions and batch_mode:
        raise ValueError("--partitions is only valid with --matrix, not --msa-dir.")

    # rcluster mutual exclusive
    if rclusterf is not None and rcluster_max is not None:
        raise ValueError("--rclusterf and --rcluster-max are mutually exclusive.")

    # rcluster without partitions warn
    if (rclusterf is not None or rcluster_max is not None) and not partitions:
        if not quiet:
            warnings.warn(
                "--rclusterf/--rcluster-max have no effect without --partitions.", UserWarning
            )

    # MF mode branch support
    if modelfinder == "MF":
        if boot is not None or alrt is not None or bnni:
            if not quiet:
                warnings.warn(
                    "Branch support flags (--boot, --alrt, --bnni) are ignored in MF (model-only) mode.",
                    UserWarning,
                )

    # bnni without boot
    if bnni and (boot is None or boot == 0):
        if not quiet:
            warnings.warn("--bnni has no effect without --boot.", UserWarning)

    # Prefix in batch mode
    if prefix and batch_mode:
        if not quiet:
            warnings.warn("--prefix ignored in batch mode; gene names used as prefix.", UserWarning)

    # PMSF validations
    if pmsf_base_model:
        if model not in AA_MIXTURE_MODELS:
            raise ValueError(
                f"--pmsf-base-model is only valid with AA mixture models (C10-C60). Got: {model}"
            )
        if not guide_tree:
            raise ValueError(
                "PMSF mode requires --guide-tree with --model in C10-C60 range."
            )

    # Model validation
    _validate_model(model=model, seq_type=seq_type, modelfinder=modelfinder, custom_model=custom_model)

    # PMSF base model validation
    if pmsf_base_model and model in AA_MIXTURE_MODELS:
        _validate_pmsf_base_model(pmsf_base_model)

    # qmax only with MIX+MF
    if qmax is not None and model != "MIX+MF":
        warnings.warn("--qmax only takes effect with --model MIX+MF.", UserWarning)


# ===================================================================
# Threads parsing
# ===================================================================

def _parse_threads(value: str | None, batch_mode: bool) -> int | str:
    if value is None:
        return 4 if batch_mode else "auto"
    v = str(value).strip().lower()
    if v == "auto":
        if batch_mode:
            raise ValueError("--threads 'auto' is not valid in batch mode; use a numeric value.")
        return "auto"
    try:
        n = int(v)
    except ValueError:
        raise ValueError(f"--threads must be a positive integer or 'auto', got: {value!r}")
    if n < 1:
        raise ValueError(f"--threads must be >= 1, got: {n}")
    return n

# ===================================================================
# .iqtree report parser
# ===================================================================

_LOG_LIKE_RE = _re.compile(r"Log-likelihood(?: of the tree)?:\s+([-\d.eE+]+)", _re.IGNORECASE)
_MODEL_SELECTED_RE = _re.compile(r"Best-fit model(?: according to \w+)?:\s*(\S+)", _re.IGNORECASE)


def _parse_iqtree_report(iqtree_path: Path) -> dict[str, Any]:
    result: dict[str, Any] = {
        "log_likelihood": None,
        "model_selected": None,
    }
    if not iqtree_path.exists():
        return result
    try:
        text = iqtree_path.read_text()
    except Exception:
        return result

    m = _LOG_LIKE_RE.search(text)
    if m:
        try:
            result["log_likelihood"] = float(m.group(1))
        except ValueError:
            pass

    m = _MODEL_SELECTED_RE.search(text)
    if m:
        result["model_selected"] = m.group(1)

    return result


# ===================================================================
# Command builder
# ===================================================================

def _build_iqtree_cmd(
    input_path: Path,
    prefix: Path,
    *,
    model_string: str,
    seq_type: str,
    cmd_seq_type: str | None = None,
    boot: int | None,
    alrt: int | None,
    bnni: bool,
    mode: str,
    threads_arg: str,
    modelfinder: str = "none",
    executable: str = "iqtree3",
    mset: str | None = None,
    msub: str | None = None,
    partitions: str | None = None,
    rclusterf: int | None = None,
    rcluster_max: int | None = None,
    guide_tree: str | None = None,
    qmax: int | None = None,
    rate: bool = False,
    wslr: bool = False,
    constraint: str | None = None,
    outgroup: str | None = None,
    site_freq_file: str | None = None,
    tool_args: str | None = None,
    batch_mode: bool = False,
) -> list[str]:
    cmd = [executable]

    tool_tokens = set(shlex.split(tool_args)) if tool_args else set()

    if tool_args:
        _check_managed_flag_conflict(tool_args, batch_mode=batch_mode)

    # Backward-compatible: if cmd_seq_type not provided, use seq_type
    if cmd_seq_type is None:
        cmd_seq_type = seq_type

    # Always include: input, prefix (no --redo by default)
    cmd.extend(["-s", str(input_path)])
    # In single mode, --prefix is overridable; in batch mode always emit (blocked)
    if batch_mode or not _is_flag_overridden("--prefix", tool_tokens):
        cmd.extend(["--prefix", str(prefix)])

    # Model
    if not _is_flag_overridden("-m", tool_tokens):
        cmd.extend(["-m", model_string])

    # ModelFinder control: only add --mset/--msub when ModelFinder is active
    # or model_string is MIX+MF (NT heterogeneous with internal model selection)
    is_mf_mode = modelfinder in ("MF", "MFP")
    need_mset = is_mf_mode or model_string == "MIX+MF"
    if mset and not _is_flag_overridden("--mset", tool_tokens):
        if mset != "all" and need_mset:
            cmd.extend(["--mset", mset])
    if msub and not _is_flag_overridden("--msub", tool_tokens):
        if is_mf_mode:
            cmd.extend(["--msub", msub])

    # Partitions + merge (only for MF/MFP per design)
    if partitions and is_mf_mode and not _is_flag_overridden("-p", tool_tokens):
        cmd.extend(["-p", partitions])
    elif partitions and not is_mf_mode and not _is_flag_overridden("-p", tool_tokens):
        cmd.extend(["-p", partitions])
    if partitions and is_mf_mode and not _is_flag_overridden("--merge", tool_tokens):
        cmd.append("--merge")
    if rclusterf is not None and is_mf_mode and not _is_flag_overridden("--rclusterf", tool_tokens):
        cmd.extend(["--rclusterf", str(rclusterf)])
    if rcluster_max is not None and is_mf_mode and not _is_flag_overridden("--rcluster-max", tool_tokens):
        cmd.extend(["--rcluster-max", str(rcluster_max)])

    # PMSF guide tree
    if guide_tree and not _is_flag_overridden("-ft", tool_tokens):
        cmd.extend(["-ft", guide_tree])

    # MIX+MF qmax
    if qmax is not None and not _is_flag_overridden("-qmax", tool_tokens):
        cmd.extend(["-qmax", str(qmax)])

    # Tree search mode
    if mode == "fast" and not _is_flag_overridden("--fast", tool_tokens):
        cmd.append("--fast")

    # Seq type: only emit when user explicitly chose AA or NT
    if not _is_flag_overridden("--seqtype", tool_tokens):
        if cmd_seq_type == "NT":
            cmd.extend(["--seqtype", "DNA"])
        elif cmd_seq_type == "AA":
            cmd.extend(["--seqtype", "AA"])
        # "auto" → omit --seqtype entirely

    # Branch support: skip entirely in MF mode (model-only, no tree)
    # boot=0 means "no support" and is treated the same as None
    if modelfinder != "MF":
        if boot and boot > 0 and not _is_flag_overridden("-B", tool_tokens):
            cmd.extend(["-B", str(boot)])
        if alrt is not None and not _is_flag_overridden("--alrt", tool_tokens):
            cmd.extend(["--alrt", str(alrt)])
        if bnni and boot and boot > 0 and not _is_flag_overridden("--bnni", tool_tokens):
            cmd.append("--bnni")

    # Threads
    if not _is_flag_overridden("-T", tool_tokens):
        cmd.extend(shlex.split(threads_arg))

    # Output flags
    if site_freq_file and not _is_flag_overridden("-fs", tool_tokens):
        cmd.extend(["-fs", site_freq_file])
    if rate and not _is_flag_overridden("--rate", tool_tokens):
        cmd.append("--rate")
    if wslr and not _is_flag_overridden("-wslr", tool_tokens):
        cmd.append("-wslr")

    # Constraint / outgroup
    if constraint and not _is_flag_overridden("-g", tool_tokens):
        cmd.extend(["-g", constraint])
    if outgroup and not _is_flag_overridden("-o", tool_tokens):
        cmd.extend(["-o", outgroup])

    # Raw tool-args appended last
    if tool_args:
        cmd.extend(shlex.split(tool_args))

    return cmd


# ===================================================================
# Single MSA execution
# ===================================================================

def _run_one_iqtree(
    gene_path: Path,
    *,
    seq_type: str,
    cmd_seq_type: str | None = None,
    model_string: str,
    modelfinder: str,
    boot: int | None,
    alrt: int | None,
    bnni: bool,
    mode: str,
    threads_arg: str,
    log_dir: Path,
    output_dir: Path,
    executable: str = "iqtree3",
    mset: str | None = None,
    msub: str | None = None,
    partitions: str | None = None,
    rclusterf: int | None = None,
    rcluster_max: int | None = None,
    guide_tree: str | None = None,
    qmax: int | None = None,
    rate: bool = False,
    wslr: bool = False,
    constraint: str | None = None,
    outgroup: str | None = None,
    site_freq_file: str | None = None,
    prefix: str | None = None,
    tool_args: str | None = None,
    dry_run: bool = False,
    batch_mode: bool = False,
    keep_extra: bool = False,
    work_dir: Path | None = None,
) -> dict[str, Any]:
    result: dict[str, Any] = {
        "input": str(gene_path),
        "output_tree": None,
        "log_iqtree": None,
        "log_file": None,
        "n_taxa": 0,
        "log_likelihood": None,
        "model_selected": None,
    }

    if not gene_path.exists():
        return {
            **result,
            "status": "failed",
            "reason": f"input file not found: {gene_path}",
            "wall_time": 0,
            "warnings": [],
        }

    # Resolve user-provided paths to absolute so IQ-TREE can find them regardless of cwd
    gene_path = gene_path.resolve()
    if partitions:
        partitions = str(Path(partitions).resolve())
    if guide_tree:
        guide_tree = str(Path(guide_tree).resolve())
    if constraint:
        constraint = str(Path(constraint).resolve())

    stem = gene_path.stem
    is_mf_only = modelfinder == "MF"

    if cmd_seq_type is None:
        cmd_seq_type = seq_type

    try:
        n_taxa = sum(1 for _ in SeqIO.parse(str(gene_path), "fasta"))
        result["n_taxa"] = n_taxa
    except Exception:
        pass

    # --- Batch mode: temp work_dir with isolation ---
    own_work_dir = False
    if batch_mode:
        if work_dir is None:
            work_dir = Path(tempfile.mkdtemp(prefix=f"iqtree_{stem}_"))
            own_work_dir = True
        prefix_path = work_dir / stem
        effective_prefix_name = stem
    else:
        # Single mode: cwd is output_dir; use prefix name relative to it
        override_prefix = _get_tool_arg_value("--prefix", tool_args)
        effective_prefix_name = override_prefix if override_prefix else (prefix if prefix else stem)
        prefix_path = Path(effective_prefix_name)

    cmd = _build_iqtree_cmd(
        input_path=gene_path, prefix=prefix_path,
        model_string=model_string, seq_type=seq_type,
        cmd_seq_type=cmd_seq_type,
        boot=boot, alrt=alrt, bnni=bnni,
        mode=mode, threads_arg=threads_arg,
        modelfinder=modelfinder,
        executable=executable,
        mset=mset, msub=msub,
        partitions=partitions,
        rclusterf=rclusterf, rcluster_max=rcluster_max,
        guide_tree=guide_tree, qmax=qmax,
        rate=rate, wslr=wslr,
        constraint=constraint, outgroup=outgroup,
        site_freq_file=site_freq_file,
        tool_args=tool_args,
        batch_mode=batch_mode,
    )

    if batch_mode:
        out_tree = output_dir / f"{stem}.treefile"
    else:
        out_tree = output_dir / f"{effective_prefix_name}.treefile"
    out_iqtree = log_dir / f"{effective_prefix_name}.iqtree"
    out_log = log_dir / f"{effective_prefix_name}.log"

    result.update({
        "output_tree": str(out_tree) if not is_mf_only else None,
        "log_iqtree": str(out_iqtree),
        "log_file": str(out_log),
        "cmd": cmd,
    })

    if dry_run:
        return {**result, "status": "dry_run", "wall_time": 0, "warnings": []}

    warnings_list: list[str] = []
    start = _time.monotonic()

    try:
        log_dir.mkdir(parents=True, exist_ok=True)
        output_dir.mkdir(parents=True, exist_ok=True)
        if batch_mode and work_dir:
            work_dir.mkdir(parents=True, exist_ok=True)

        cwd = str(work_dir) if batch_mode else str(output_dir)

        if batch_mode:
            proc = subprocess.run(
                cmd, cwd=cwd,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
            )
            stderr_text = proc.stderr
            returncode = proc.returncode
        else:
            # Single mode: stream IQ-TREE stdout to terminal for progress visibility
            child = subprocess.Popen(
                cmd, cwd=cwd,
                stdout=None,  # inherit → terminal
                stderr=subprocess.PIPE,
                text=True,
            )
            _, stderr_text = child.communicate()
            returncode = child.returncode

        wall_time = _time.monotonic() - start

        if returncode != 0:
            return {
                **result,
                "status": "failed",
                "reason": f"iqtree3 exited with code {returncode}: {stderr_text[:200]}",
                "tool_stderr": stderr_text,
                "wall_time": wall_time,
                "warnings": warnings_list,
            }

        # Collect output files based on mode
        if batch_mode and work_dir:
            assert work_dir is not None
            prefix_base = work_dir / stem
            # Move .treefile
            src_tree = Path(f"{prefix_base}.treefile")
            if src_tree.exists() and not is_mf_only:
                output_dir.mkdir(parents=True, exist_ok=True)
                shutil.move(str(src_tree), str(out_tree))
            elif not is_mf_only:
                return {
                    **result,
                    "status": "failed",
                    "reason": "iqtree3 did not produce .treefile",
                    "wall_time": wall_time,
                    "warnings": warnings_list,
                }
            # Move .iqtree
            src_iqtree = Path(f"{prefix_base}.iqtree")
            if src_iqtree.exists():
                shutil.move(str(src_iqtree), str(out_iqtree))
            # Move .log
            src_log = Path(f"{prefix_base}.log")
            if src_log.exists():
                shutil.move(str(src_log), str(out_log))
            # Move extra generated files only if --keep-extra
            if keep_extra:
                for f in sorted(work_dir.iterdir()):
                    _name = f.name
                    if _name == f"{stem}.iqtree" or _name == f"{stem}.log":
                        continue  # already moved above
                    dest = log_dir / _name
                    if not dest.exists():
                        shutil.move(str(f), str(dest))
        else:
            # Single mode: files already in output_dir (IQ-TREE wrote there directly)
            # Just verify .treefile exists for non-MF modes
            if not is_mf_only and not out_tree.exists():
                return {
                    **result,
                    "status": "failed",
                    "reason": "iqtree3 did not produce .treefile",
                    "wall_time": wall_time,
                    "warnings": warnings_list,
                }

        # Parse .iqtree for metadata
        report = _parse_iqtree_report(out_iqtree)
        result["log_likelihood"] = report["log_likelihood"]
        result["model_selected"] = report["model_selected"]

        # MF mode: require non-empty .iqtree report
        if is_mf_only:
            if not out_iqtree.exists() or out_iqtree.stat().st_size == 0:
                return {
                    **result,
                    "status": "failed",
                    "reason": "MF mode: iqtree3 did not produce a valid .iqtree report",
                    "wall_time": wall_time,
                    "warnings": warnings_list,
                }

        # Some IQ-TREE models (for example GHOST) write one Newick tree per
        # heterotachy class. Require one or more parseable trees.
        if not is_mf_only and out_tree.exists():
            try:
                trees = Phylo.parse(str(out_tree), "newick")
                next(trees)
                for _ in trees:
                    pass
            except Exception as e:
                return {
                    **result,
                    "status": "failed",
                    "reason": f"iqtree3 produced unparseable Newick output: {e}",
                    "wall_time": wall_time,
                    "warnings": warnings_list,
                }

        return {
            **result,
            "status": "success",
            "wall_time": wall_time,
            "tool_stderr": stderr_text,
            "warnings": warnings_list,
        }

    except Exception as exc:
        return {
            **result,
            "status": "failed",
            "reason": str(exc),
            "wall_time": _time.monotonic() - start,
            "warnings": warnings_list,
        }
    finally:
        if own_work_dir and work_dir and work_dir.exists():
            shutil.rmtree(work_dir, ignore_errors=True)


# ===================================================================
# Main entry point
# ===================================================================

import json
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import Callable

from phyloai.core.checkpoint import load_checkpoint, save_checkpoint_atomic, validate_resume_params
from phyloai.tree.checkpoint_helpers import mark_task


CHECKPOINT_FLUSH_INTERVAL = 2.0


def _resolved_iqtree_params(
    *,
    msa_dir: Path | None,
    matrix: Path | None,
    output_dir: Path,
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
    partitions: str | None,
    rclusterf: int | None,
    rcluster_max: int | None,
    pmsf_base_model: str | None,
    guide_tree: str | None,
    site_freq_file: str | None,
    qmax: int | None,
    rate: bool,
    wslr: bool,
    constraint: str | None,
    outgroup: str | None,
    prefix: str | None,
    threads: str,
    overwrite: bool,
    iqtree_path: str | None,
    tool_args: str | None,
) -> dict[str, Any]:
    return {
        "msa_dir": str(msa_dir) if msa_dir else None,
        "matrix": str(matrix) if matrix else None,
        "seq_type": seq_type,
        "model": model,
        "state_freq": state_freq,
        "rate_heterogeneity": rate_heterogeneity,
        "modelfinder": modelfinder,
        "mset": mset,
        "msub": msub,
        "mode": mode,
        "boot": boot,
        "alrt": alrt,
        "bnni": bnni,
        "partitions": partitions,
        "rclusterf": rclusterf,
        "rcluster_max": rcluster_max,
        "pmsf_base_model": pmsf_base_model,
        "guide_tree": guide_tree,
        "site_freq_file": site_freq_file,
        "qmax": qmax,
        "rate": rate,
        "wslr": wslr,
        "constraint": constraint,
        "outgroup": outgroup,
        "prefix": prefix,
        "output_dir": str(output_dir),
        "threads": threads,
        "overwrite": overwrite,
        "iqtree_path": iqtree_path,
        "tool_args": tool_args,
    }


def run_iqtree(
    *,
    msa_dir: Path | None = None,
    matrix: Path | None = None,
    output_dir: Path,
    seq_type: str = "auto",
    model: str | None = None,
    state_freq: str = "+F",
    rate_heterogeneity: str = "+R4",
    modelfinder: str = "none",
    mset: str | None = None,
    msub: str | None = None,
    mode: str = "normal",
    boot: int | None = 1000,
    alrt: int | None = None,
    bnni: bool = False,
    partitions: str | None = None,
    rclusterf: int | None = None,
    rcluster_max: int | None = None,
    pmsf_base_model: str | None = None,
    guide_tree: str | None = None,
    site_freq_file: str | Path | None = None,
    qmax: int | None = None,
    rate: bool = False,
    wslr: bool = False,
    constraint: str | None = None,
    outgroup: str | None = None,
    prefix: str | None = None,
    threads: str | int | None = None,
    iqtree_path: str | None = None,
    tool_args: str | None = None,
    overwrite: bool = False,
    resume: bool = False,
    dry_run: bool = False,
    keep_extra: bool = False,
    quiet: bool = False,
    progress_callback: Callable[[Path], None] | None = None,
) -> dict[str, Any]:
    run_start = _time.monotonic()
    output_dir = output_dir.resolve()

    # --- Input mutual exclusivity ---
    if (msa_dir is None and matrix is None) or (msa_dir is not None and matrix is not None):
        raise ValueError("Either --msa-dir or --matrix must be provided (not both).")

    # --- overwrite/resume mutual exclusivity ---
    if overwrite and resume:
        raise ValueError("--overwrite and --resume are mutually exclusive.")

    # --- Resolve executable ---
    iqtree_exe = _resolve_iqtree_path(iqtree_path, dry_run)

    batch_mode = msa_dir is not None
    raw_site_freq_override = bool(tool_args and _is_flag_overridden("-fs", set(shlex.split(tool_args))))
    custom_model = False
    resolved_site_freq_file = _resolve_site_freq_file(site_freq_file)
    if model is not None:
        custom_model_path = _resolve_custom_model_path(model)
        if custom_model_path is not None:
            model = custom_model_path
            custom_model = True
    n_resume_skipped = 0
    resolved_seq_type = seq_type

    trees_dir = output_dir / "trees"
    logs_dir = output_dir / "logs"

    # --- Parse threads ---
    threads_str = str(threads) if threads is not None else None
    parsed_threads = _parse_threads(threads_str, batch_mode)
    threads_spec: str
    threads_int: int
    if isinstance(parsed_threads, int):
        threads_spec = str(parsed_threads)
        threads_int = parsed_threads
    else:
        threads_spec = "AUTO"
        threads_int = 1  # unused in single mode with AUTO

    # --- Single mode ---
    if not batch_mode:
        assert matrix is not None
        matrix_ext = matrix.suffix.lower()
        if matrix_ext not in IQTREE_COMPATIBLE_EXTENSIONS:
            raise ValueError(
                f"--matrix has unsupported extension: {matrix.suffix}. "
                f"Supported: {sorted(IQTREE_COMPATIBLE_EXTENSIONS)}"
            )
        fmt = _detect_file_format(matrix_ext)
        try:
            recs = list(SeqIO.parse(str(matrix), fmt))
        except Exception:
            recs = []
        if not recs:
            raise ValueError(f"Cannot parse --matrix: {matrix} (format: {fmt})")
        if seq_type == "auto":
            resolved_seq_type = detect_seq_type([str(r.seq) for r in recs]) if recs else "AA"
            cmd_seq_type = "auto"
        else:
            cmd_seq_type = seq_type
            if recs:
                sample = [str(r.seq) for r in recs[:10]]
                resolved_sample = detect_seq_type(sample)
                if resolved_sample != seq_type:
                    raise ValueError(
                        f"--seq-type {seq_type} but detected {resolved_sample} in {matrix}"
                    )
            resolved_seq_type = seq_type

        # Output directory handling
        if not dry_run:
            if overwrite and output_dir.exists():
                shutil.rmtree(output_dir)
            if not overwrite and not resume and output_dir.exists() and any(output_dir.iterdir()):
                raise ValueError(
                    f"Output directory {output_dir} already exists and is non-empty. "
                    "Use --overwrite to replace."
                )

        # Resolve model default
        if model is None:
            model = "GTR" if resolved_seq_type == "NT" else "LG"
        if mset is None:
            if modelfinder in ("MF", "MFP"):
                if resolved_seq_type == "NT":
                    mset = "GTR,HKY"
                elif resolved_seq_type == "AA":
                    mset = "all" if (msub and msub != "nuclear") else "LG,WAG"
            elif model == "MIX+MF":
                mset = "GTR,HKY"

        if pmsf_base_model is None and model in AA_MIXTURE_MODELS and guide_tree is not None:
            pmsf_base_model = "LG"

        # Resolve partition merge defaults
        resolved_rclusterf = rclusterf
        resolved_rcluster_max = rcluster_max
        if partitions is not None and modelfinder in ("MF", "MFP"):
            if resolved_rclusterf is None and resolved_rcluster_max is None:
                resolved_rclusterf = 10

        # Resolve MIX+MF qmax default
        resolved_qmax = qmax
        if model == "MIX+MF" and resolved_qmax is None:
            resolved_qmax = 10

        # MF mode: normalize branch support flags (model-only, no tree)
        resolved_boot = boot
        resolved_alrt = alrt
        resolved_bnni = bnni
        # Validations
        _run_validations(
            batch_mode=False, seq_type=resolved_seq_type,
            modelfinder=modelfinder, model=model,
            partitions=partitions, guide_tree=guide_tree,
            boot=resolved_boot, alrt=resolved_alrt, bnni=resolved_bnni,
            prefix=prefix, pmsf_base_model=pmsf_base_model,
            rclusterf=resolved_rclusterf, rcluster_max=resolved_rcluster_max,
            qmax=resolved_qmax,
            state_freq=state_freq, custom_model=custom_model,
            site_freq_file=resolved_site_freq_file, tool_args=tool_args,
            quiet=quiet,
        )

        model_string = _build_model_string(
            model=model, state_freq=state_freq, rate_heterogeneity=rate_heterogeneity,
            modelfinder=modelfinder,
            pmsf_base_model=pmsf_base_model if model in AA_MIXTURE_MODELS else None,
        )

        single_threads_arg = f"-T {threads_spec}"

        result = _run_one_iqtree(
            gene_path=matrix, seq_type=resolved_seq_type,
            cmd_seq_type=cmd_seq_type,
            model_string=model_string, modelfinder=modelfinder,
            boot=resolved_boot, alrt=resolved_alrt, bnni=resolved_bnni,
            mode=mode, threads_arg=single_threads_arg,
            log_dir=output_dir, output_dir=output_dir,
            executable=iqtree_exe,
            mset=mset, msub=msub,
            partitions=partitions,
            rclusterf=resolved_rclusterf, rcluster_max=resolved_rcluster_max,
            guide_tree=guide_tree, qmax=resolved_qmax,
            rate=rate, wslr=wslr,
            constraint=constraint, outgroup=outgroup,
            site_freq_file=resolved_site_freq_file,
            prefix=prefix,
            tool_args=tool_args,
            dry_run=dry_run, batch_mode=False, keep_extra=keep_extra,
        )

        # Single mode: separate result into success vs failed for correct assembly
        single_results: list[dict[str, Any]] = []
        single_failed: list[dict[str, Any]] = []
        if result["status"] in {"success", "dry_run"}:
            single_results.append(result)
        else:
            single_failed.append(result)

        return _assemble_iqtree_result(
            run_start=run_start, iqtree_exe=iqtree_exe,
            batch_mode=False, results=single_results,
            failed_results=single_failed,
            resolved_seq_type=resolved_seq_type, model=model,
            model_string=model_string,
            modelfinder=modelfinder, mode=mode,
            boot=boot, alrt=alrt, bnni=bnni,
            partitions=partitions,
            rclusterf=resolved_rclusterf, rcluster_max=resolved_rcluster_max,
            guide_tree=guide_tree,
            output_dir=output_dir,
            msa_dir=msa_dir, matrix=matrix,
            iqtree_path=iqtree_path, tool_args=tool_args,
            overwrite=overwrite, threads=str(threads_spec),
            skipped_input=[], dry_run=dry_run, resume=resume, keep_extra=keep_extra,
            state_freq=state_freq, rate_heterogeneity=rate_heterogeneity,
            mset=mset, msub=msub,
            pmsf_base_model=pmsf_base_model,
            site_freq_file=(None if raw_site_freq_override else resolved_site_freq_file),
            qmax=resolved_qmax, rate=rate, wslr=wslr,
            constraint=constraint, outgroup=outgroup,
            prefix=prefix,
            quiet=quiet,
        )

    # --- Batch mode ---
    assert msa_dir is not None
    found, skipped_input = _scan_input_iqtree(msa_dir)
    if not found:
        raise ValueError("No valid input files found in --msa-dir")

    stems = [p.stem for p in found]
    dupes = {s for s in stems if stems.count(s) > 1}
    if dupes:
        colliding = sorted(p.name for p in found if p.stem in dupes)
        raise ValueError(
            "Duplicate output stems detected in --msa-dir. Files with different extensions "
            "but the same base name collide:\n  " + "\n  ".join(colliding) +
            "\nRename files so each stem is unique within the directory."
        )

    declared = None if seq_type == "auto" else seq_type
    resolved_seq_type, offending = _validate_seq_types_iqtree(found, declared_type=declared)
    if resolved_seq_type is None:
        offending_strs = [_format_offender(o) for o in offending[:10]]
        raise ValueError("Mixed sequence types in --msa-dir:\n" + "\n".join(offending_strs))
    if offending:
        offending_strs = [_format_offender(o) for o in offending[:10]]
        raise ValueError(
            f"Files with wrong --seq-type ({declared}) in --msa-dir:\n" + "\n".join(offending_strs)
        )

    if model is None:
        model = "GTR" if resolved_seq_type == "NT" else "LG"
    if mset is None:
        if modelfinder in ("MF", "MFP"):
            if resolved_seq_type == "NT":
                mset = "GTR,HKY"
            elif resolved_seq_type == "AA":
                mset = "all" if (msub and msub != "nuclear") else "LG,WAG"
        elif model == "MIX+MF":
            mset = "GTR,HKY"

    if pmsf_base_model is None and model in AA_MIXTURE_MODELS and guide_tree is not None:
        pmsf_base_model = "LG"

    # Resolve partition merge defaults
    resolved_rclusterf = rclusterf
    resolved_rcluster_max = rcluster_max
    if partitions is not None and modelfinder in ("MF", "MFP"):
        if resolved_rclusterf is None and resolved_rcluster_max is None:
            resolved_rclusterf = 10

    # Resolve MIX+MF qmax default (batch: homogeneous only, but safe)
    resolved_qmax = qmax
    if model == "MIX+MF" and resolved_qmax is None:
        resolved_qmax = 10

    # MF mode: normalize branch support flags
    resolved_boot = boot
    resolved_alrt = alrt
    resolved_bnni = bnni
    if modelfinder == "MF":
        if boot is not None or alrt is not None or bnni:
            if not quiet:
                warnings.warn(
                    "Branch support flags ignored in MF (model-only) mode.", UserWarning
                )
        resolved_boot = None
        resolved_alrt = None
        resolved_bnni = False

    # cmd_seq_type for batch: always auto (IQ-TREE infers from input)
    cmd_seq_type = seq_type if seq_type != "auto" else "auto"

    _run_validations(
        batch_mode=True, seq_type=resolved_seq_type,
        modelfinder=modelfinder, model=model,
        partitions=partitions, guide_tree=guide_tree,
        boot=resolved_boot, alrt=resolved_alrt, bnni=resolved_bnni,
        prefix=prefix, pmsf_base_model=pmsf_base_model,
        rclusterf=resolved_rclusterf, rcluster_max=resolved_rcluster_max,
        qmax=resolved_qmax,
        state_freq=state_freq, custom_model=custom_model,
        site_freq_file=resolved_site_freq_file, tool_args=tool_args,
        quiet=quiet,
    )

    model_string = _build_model_string(
        model=model, state_freq=state_freq, rate_heterogeneity=rate_heterogeneity,
        modelfinder=modelfinder,
        pmsf_base_model=pmsf_base_model if model in AA_MIXTURE_MODELS else None,
    )

    is_mf_only = modelfinder == "MF"
    batch_threads_arg = "-T 1"

    # --- Output directory handling (batch) ---
    if not dry_run:
        if overwrite and output_dir.exists():
            shutil.rmtree(output_dir)
        if not overwrite and not resume and output_dir.exists() and any(output_dir.iterdir()):
            raise ValueError(
                f"Output directory {output_dir} already exists and is non-empty. "
                "Use --overwrite to replace."
            )
        trees_dir.mkdir(parents=True, exist_ok=True)
        logs_dir.mkdir(parents=True, exist_ok=True)

    # --- Checkpoint ---
    checkpoint: Any = None
    ckpt_path = output_dir / "checkpoint.json"

    if not dry_run:
        if resume:
            if not ckpt_path.exists():
                raise ValueError(f"--resume requires {ckpt_path}, not found")
            checkpoint = load_checkpoint(ckpt_path)
            resolved_params = _resolved_iqtree_params(
                msa_dir=msa_dir, matrix=matrix, output_dir=output_dir,
                seq_type=resolved_seq_type, model=model,
                state_freq=state_freq, rate_heterogeneity=rate_heterogeneity,
                modelfinder=modelfinder,
                mset=mset, msub=msub, mode=mode,
                boot=boot, alrt=alrt, bnni=bnni,
                partitions=partitions, rclusterf=rclusterf, rcluster_max=rcluster_max,
                pmsf_base_model=pmsf_base_model, guide_tree=guide_tree,
                site_freq_file=(None if raw_site_freq_override else resolved_site_freq_file), qmax=qmax,
                rate=rate, wslr=wslr, constraint=constraint, outgroup=outgroup,
                prefix=prefix, threads=str(threads_spec),
                overwrite=overwrite, iqtree_path=iqtree_path, tool_args=tool_args,
            )
            validate_resume_params(checkpoint, resolved_params, step="tree.ml.iqtree")
            if checkpoint.status == "success":
                return _reconstruct_result(output_dir, run_start)

            from phyloai.tree.checkpoint_helpers import plan_resume_iqtree
            to_run_ids, skipped_ids = plan_resume_iqtree(checkpoint, is_mf_only=is_mf_only)
            n_resume_skipped = len(skipped_ids)
            if not to_run_ids:
                checkpoint.status = "success"
                save_checkpoint_atomic(checkpoint, ckpt_path)
                return _reconstruct_result(output_dir, run_start)
            found = [Path(task.input) for task in checkpoint.tasks if task.task_id in to_run_ids]
        else:
            from phyloai.tree.checkpoint_helpers import build_initial_iqtree_checkpoint
            iqtree_cmd_parts = ["phyloai", "tree", "ml", "iqtree", "--msa-dir", str(msa_dir)]
            iqtree_cmd_parts.extend([
                "--seq-type", resolved_seq_type, "--model", model,
                "--state-freq", state_freq, "--rate-heterogeneity", rate_heterogeneity,
                "--mode", mode,
                "-o", str(output_dir),
                "-t", str(threads_spec),
            ])
            if modelfinder != "none":
                iqtree_cmd_parts.extend(["--modelfinder", modelfinder])
            if mset is not None:
                iqtree_cmd_parts.extend(["--mset", mset])
            if msub is not None:
                iqtree_cmd_parts.extend(["--msub", msub])
            if resolved_boot is not None:
                iqtree_cmd_parts.extend(["--boot", str(resolved_boot)])
            if resolved_alrt is not None:
                iqtree_cmd_parts.extend(["--alrt", str(resolved_alrt)])
            if resolved_bnni:
                iqtree_cmd_parts.append("--bnni")
            if partitions:
                iqtree_cmd_parts.extend(["--partitions", partitions])
            if resolved_rclusterf is not None:
                iqtree_cmd_parts.extend(["--rclusterf", str(resolved_rclusterf)])
            if resolved_rcluster_max is not None:
                iqtree_cmd_parts.extend(["--rcluster-max", str(resolved_rcluster_max)])
            if pmsf_base_model is not None:
                iqtree_cmd_parts.extend(["--pmsf-base-model", pmsf_base_model])
            if guide_tree is not None:
                iqtree_cmd_parts.extend(["--guide-tree", guide_tree])
            if resolved_site_freq_file is not None:
                iqtree_cmd_parts.extend(["--site-freq-file", resolved_site_freq_file])
            if resolved_qmax is not None:
                iqtree_cmd_parts.extend(["--qmax", str(resolved_qmax)])
            if rate:
                iqtree_cmd_parts.append("--rate")
            if wslr:
                iqtree_cmd_parts.append("--wslr")
            if constraint is not None:
                iqtree_cmd_parts.extend(["--constraint", constraint])
            if outgroup is not None:
                iqtree_cmd_parts.extend(["--outgroup", outgroup])
            if prefix is not None:
                iqtree_cmd_parts.extend(["--prefix", prefix])
            if iqtree_path:
                iqtree_cmd_parts.extend(["--iqtree-path", iqtree_path])
            if tool_args:
                if " " in tool_args:
                    iqtree_cmd_parts.append(f"--tool-args '{tool_args}'")
                else:
                    iqtree_cmd_parts.extend(["--tool-args", tool_args])
            if overwrite:
                iqtree_cmd_parts.append("--overwrite")
            checkpoint = build_initial_iqtree_checkpoint(
                step="tree.ml.iqtree",
                command=" ".join(iqtree_cmd_parts),
                params=_resolved_iqtree_params(
                    msa_dir=msa_dir, matrix=matrix, output_dir=output_dir,
                    seq_type=resolved_seq_type, model=model,
                    state_freq=state_freq, rate_heterogeneity=rate_heterogeneity,
                    modelfinder=modelfinder,
                    mset=mset, msub=msub, mode=mode,
                    boot=boot, alrt=alrt, bnni=bnni,
                    partitions=partitions, rclusterf=rclusterf, rcluster_max=rcluster_max,
                    pmsf_base_model=pmsf_base_model, guide_tree=guide_tree,
                    site_freq_file=(None if raw_site_freq_override else resolved_site_freq_file), qmax=qmax,
                    rate=rate, wslr=wslr, constraint=constraint, outgroup=outgroup,
                    prefix=prefix, threads=str(threads_spec),
                    overwrite=overwrite, iqtree_path=iqtree_path, tool_args=tool_args,
                ),
                inputs=found, trees_dir=trees_dir, logs_dir=logs_dir,
            )
            save_checkpoint_atomic(checkpoint, ckpt_path)

    # --- Build worker args (module-level function for ProcessPoolExecutor pickle-compat) ---
    worker_args: list[tuple] = [
        (
            p,
            resolved_seq_type,
            cmd_seq_type,
            model_string,
            modelfinder,
            resolved_boot,
            resolved_alrt,
            resolved_bnni,
            mode,
            batch_threads_arg,
            logs_dir,
            trees_dir,
            iqtree_exe,
            mset,
            msub,
            partitions,
            resolved_rclusterf,
            resolved_rcluster_max,
            guide_tree,
            resolved_qmax,
            rate,
            wslr,
            constraint,
            outgroup,
            tool_args,
            dry_run,
            keep_extra,
        )
        for p in found
    ]

    file_results: list[dict[str, Any]] = []
    failed_results: list[dict[str, Any]] = []

    _ckpt_write = checkpoint is not None and not dry_run
    _last_flush = _time.monotonic()

    def _maybe_flush(*, force: bool = False) -> None:
        nonlocal _last_flush
        if not _ckpt_write:
            return
        now = _time.monotonic()
        if force or (now - _last_flush) >= CHECKPOINT_FLUSH_INTERVAL:
            save_checkpoint_atomic(checkpoint, ckpt_path)
            _last_flush = now

    interrupted = False
    try:
        if dry_run:
            for args in worker_args:
                result = _run_one_iqtree(
                    gene_path=args[0], seq_type=args[1],
                    cmd_seq_type=args[2], model_string=args[3],
                    modelfinder=args[4], boot=args[5], alrt=args[6], bnni=args[7],
                    mode=args[8], threads_arg=args[9],
                    log_dir=args[10], output_dir=args[11],
                    executable=args[12], mset=args[13], msub=args[14],
                    partitions=args[15], rclusterf=args[16], rcluster_max=args[17],
                    guide_tree=args[18], qmax=args[19],
                    rate=args[20], wslr=args[21],
                    constraint=args[22], outgroup=args[23],
                    tool_args=args[24],
                    dry_run=args[25], batch_mode=True, keep_extra=args[26],
                )
                file_results.append(result)
                if progress_callback:
                    progress_callback(args[0])
        else:
            with ProcessPoolExecutor(max_workers=threads_int) as pool:
                futures = {
                    pool.submit(
                        _run_one_iqtree,
                        gene_path=args[0], seq_type=args[1],
                        cmd_seq_type=args[2], model_string=args[3],
                        modelfinder=args[4], boot=args[5], alrt=args[6], bnni=args[7],
                        mode=args[8], threads_arg=args[9],
                        log_dir=args[10], output_dir=args[11],
                        executable=args[12], mset=args[13], msub=args[14],
                        partitions=args[15], rclusterf=args[16], rcluster_max=args[17],
                        guide_tree=args[18], qmax=args[19],
                        rate=args[20], wslr=args[21],
                        constraint=args[22], outgroup=args[23],
                        tool_args=args[24],
                        dry_run=args[25], batch_mode=True, keep_extra=args[26],
                    ): args[0]
                    for args in worker_args
                }
                for future in as_completed(futures):
                    gene_path = futures[future]
                    result = future.result()
                    task_id = gene_path.stem

                    if result["status"] == "success":
                        file_results.append(result)
                        mark_task(checkpoint, task_id, status="success")
                    elif result["status"] == "failed":
                        failed_results.append(result)
                        mark_task(checkpoint, task_id, status="failed",
                                  reason=result.get("reason"))
                    else:
                        skipped_input.append({
                            "path": result.get("input", ""),
                            "reason": result.get("reason", "unknown"),
                        })

                    if progress_callback:
                        progress_callback(gene_path)
                    _maybe_flush()

    except KeyboardInterrupt:
        interrupted = True

    from datetime import datetime as _dt_cls, timezone as _tz
    if _ckpt_write:
        if interrupted:
            checkpoint.status = "interrupted"
        elif failed_results:
            checkpoint.status = "running"
        else:
            checkpoint.status = "success"
            checkpoint.completed_at = _dt_cls.now(_tz.utc).isoformat(timespec="seconds")
        save_checkpoint_atomic(checkpoint, ckpt_path, fsync=True)
    if interrupted:
        raise KeyboardInterrupt

    return _assemble_iqtree_result(
        run_start=run_start, iqtree_exe=iqtree_exe,
        batch_mode=True, results=file_results,
        failed_results=failed_results,
        resolved_seq_type=resolved_seq_type, model=model,
        model_string=model_string,
        modelfinder=modelfinder, mode=mode,
        boot=boot, alrt=alrt, bnni=bnni,
        partitions=partitions,
        rclusterf=resolved_rclusterf, rcluster_max=resolved_rcluster_max,
        guide_tree=guide_tree,
        output_dir=output_dir,
        msa_dir=msa_dir, matrix=matrix,
        iqtree_path=iqtree_path, tool_args=tool_args,
        overwrite=overwrite, threads=str(threads_spec),
        skipped_input=skipped_input,
        n_resume_skipped=n_resume_skipped,
        dry_run=dry_run, resume=resume, keep_extra=keep_extra,
        state_freq=state_freq, rate_heterogeneity=rate_heterogeneity,
        mset=mset, msub=msub,
        pmsf_base_model=pmsf_base_model,
        site_freq_file=(None if raw_site_freq_override else resolved_site_freq_file),
        qmax=resolved_qmax, rate=rate, wslr=wslr,
        constraint=constraint, outgroup=outgroup,
        prefix=prefix,
        quiet=quiet,
    )


def _assemble_iqtree_result(
    *,
    run_start: float,
    iqtree_exe: str,
    batch_mode: bool,
    results: list[dict[str, Any]],
    resolved_seq_type: str,
    model: str,
    model_string: str,
    modelfinder: str,
    mode: str,
    boot: int | None,
    alrt: int | None,
    bnni: bool,
    partitions: str | None,
    rclusterf: int | None,
    rcluster_max: int | None,
    guide_tree: str | None,
    output_dir: Path,
    msa_dir: Path | None,
    matrix: Path | None,
    iqtree_path: str | None,
    tool_args: str | None,
    overwrite: bool,
    threads: str,
    skipped_input: list[dict[str, str]],
    failed_results: list[dict[str, Any]] | None = None,
    n_resume_skipped: int = 0,
    dry_run: bool = False,
    resume: bool = False,
    keep_extra: bool = False,
    quiet: bool = False,
    state_freq: str = "+F",
    rate_heterogeneity: str = "+R4",
    mset: str | None = None,
    msub: str | None = None,
    pmsf_base_model: str | None = None,
    site_freq_file: str | None = None,
    qmax: int | None = None,
    rate: bool = False,
    wslr: bool = False,
    constraint: str | None = None,
    outgroup: str | None = None,
    prefix: str | None = None,
) -> dict[str, Any]:
    if failed_results is None:
        failed_results = []

    all_ok = [r for r in results if r["status"] in {"success", "dry_run"}]
    n_successful = len(all_ok) + n_resume_skipped
    n_failed = len(failed_results)
    n_skipped = len(skipped_input)

    # Status: error only if NO successful tasks AND something went wrong
    is_error = n_successful == 0 and (n_failed > 0 or n_skipped > 0)

    error_msg = "All IQ-TREE runs failed" if is_error else None

    mean_wall_time = 0.0
    if len(all_ok) > 0:
        total_wall = sum(r.get("wall_time", 0.0) for r in all_ok)
        mean_wall_time = total_wall / len(all_ok)

    try:
        versions = _detect_iqtree_version(iqtree_exe)
    except Exception:
        versions = {"iqtree3": "unknown"}

    cmd_parts = ["phyloai", "tree", "ml", "iqtree"]
    if batch_mode:
        cmd_parts.extend(["--msa-dir", str(msa_dir)])
    else:
        cmd_parts.extend(["--matrix", str(matrix)])
    cmd_parts.extend([
        "--seq-type", resolved_seq_type,
    ])
    if modelfinder == "none":
        cmd_parts.extend([
            "--model", model,
            "--state-freq", state_freq,
            "--rate-heterogeneity", rate_heterogeneity,
        ])
    cmd_parts.extend([
        "--mode", mode,
        "-o", str(output_dir),
    ])
    if threads != "4":
        cmd_parts.extend(["-t", threads])
    if modelfinder != "none":
        cmd_parts.extend(["--modelfinder", modelfinder])
    if mset is not None:
        cmd_parts.extend(["--mset", mset])
    if msub is not None:
        cmd_parts.extend(["--msub", msub])
    if boot is not None:
        cmd_parts.extend(["--boot", str(boot)])
    if alrt is not None:
        cmd_parts.extend(["--alrt", str(alrt)])
    if bnni:
        cmd_parts.append("--bnni")
    if partitions:
        cmd_parts.extend(["--partitions", partitions])
    if rclusterf is not None:
        cmd_parts.extend(["--rclusterf", str(rclusterf)])
    if rcluster_max is not None:
        cmd_parts.extend(["--rcluster-max", str(rcluster_max)])
    if pmsf_base_model is not None:
        cmd_parts.extend(["--pmsf-base-model", pmsf_base_model])
    if guide_tree is not None:
        cmd_parts.extend(["--guide-tree", guide_tree])
    if site_freq_file is not None:
        cmd_parts.extend(["--site-freq-file", site_freq_file])
    if qmax is not None:
        cmd_parts.extend(["--qmax", str(qmax)])
    if rate:
        cmd_parts.append("--rate")
    if wslr:
        cmd_parts.append("--wslr")
    if constraint is not None:
        cmd_parts.extend(["--constraint", constraint])
    if outgroup is not None:
        cmd_parts.extend(["--outgroup", outgroup])
    if prefix is not None:
        cmd_parts.extend(["--prefix", prefix])
    if iqtree_path:
        cmd_parts.extend(["--iqtree-path", iqtree_path])
    if tool_args:
        if " " in tool_args:
            cmd_parts.append(f"--tool-args '{tool_args}'")
        else:
            cmd_parts.extend(["--tool-args", tool_args])
    if overwrite:
        cmd_parts.append("--overwrite")
    if dry_run:
        cmd_parts.append("--dry-run")
    if resume:
        cmd_parts.append("--resume")
    if keep_extra:
        cmd_parts.append("--keep-extra")
    cmd_str = " ".join(cmd_parts)

    workflow = _classify_workflow(
        modelfinder=modelfinder, model=model, seq_type=resolved_seq_type,
        partitions=partitions, rclusterf=rclusterf, rcluster_max=rcluster_max,
        guide_tree=guide_tree,
    )

    single_log_likelihood = None
    if not batch_mode and len(all_ok) > 0:
        single_log_likelihood = all_ok[0].get("log_likelihood")

    model_selected = None
    if len(all_ok) > 0:
        model_selected = all_ok[0].get("model_selected")
    if modelfinder == "none":
        model_selected = model_string

    if batch_mode:
        ok_files = []
        for r in all_ok:
            entry = {k: v for k, v in r.items() if k != "tool_stderr"}
            for key in ("log_file", "log_iqtree"):
                log_path = Path(entry.get(key, ""))
                if key in entry and log_path.is_absolute():
                    try:
                        entry[key] = str(log_path.relative_to(output_dir))
                    except ValueError:
                        pass
            ok_files.append(entry)
        ok_failed = [
            {k: v for k, v in r.items() if k != "tool_stderr"}
            for r in failed_results
        ]
        for entry in ok_failed:
            for key in ("log_file", "log_iqtree"):
                log_path = Path(entry.get(key, ""))
                if key in entry and log_path.is_absolute():
                    try:
                        entry[key] = str(log_path.relative_to(output_dir))
                    except ValueError:
                        pass
        data_block: dict[str, Any] = {
            "summary": {
                "n_input_files": len(results) + n_failed + n_skipped + n_resume_skipped,
                "n_trees": n_successful if not (modelfinder == "MF") else 0,
                "n_failed": n_failed,
                "n_skipped": n_skipped,
                "n_resume_skipped": n_resume_skipped,
                "mean_wall_time": mean_wall_time,
                "mode": "--msa-dir" if batch_mode else "--matrix",
            },
            "files": ok_files,
            "failed": ok_failed,
            "skipped": skipped_input,
            "warnings": [],
        }
    else:
        first = results[0] if results else {}
        data_block = {
            "cmd": first.get("cmd", []),
            "tool_stderr": first.get("tool_stderr", ""),
            "output": first.get("output_tree", ""),
            "warnings": first.get("warnings", []),
        }

    payload: dict[str, Any] = {
        "status": "error" if is_error else "success",
        "command": cmd_str,
        "wall_time": _time.monotonic() - run_start,
        "tool_versions": versions,
        "params": {
            "msa_dir": str(msa_dir) if msa_dir else None,
            "matrix": str(matrix) if matrix else None,
            "seq_type": resolved_seq_type,
            "model": model if modelfinder == "none" else None,
            "state_freq": state_freq if modelfinder == "none" else None,
            "rate_heterogeneity": rate_heterogeneity if modelfinder == "none" else None,
            "modelfinder": modelfinder,
            "mset": mset,
            "msub": msub,
            "partitions": partitions,
            "rclusterf": rclusterf,
            "rcluster_max": rcluster_max,
            "mode": mode,
            "boot": boot if (boot and boot > 0) else None,
            "alrt": alrt,
            "bnni": bnni,
            "pmsf_base_model": pmsf_base_model,
            "guide_tree": guide_tree,
            "site_freq_file": site_freq_file,
            "qmax": qmax,
            "rate": rate,
            "wslr": wslr,
            "constraint": constraint,
            "outgroup": outgroup,
            "prefix": prefix,
            "output_dir": str(output_dir),
            "threads": threads,
            "overwrite": overwrite,
            "resume": resume,
            "dry_run": dry_run,
            "keep_extra": keep_extra,
            "quiet": quiet,
            "iqtree_path": iqtree_path,
            "tool_args": tool_args,
        },
        "key_results": {
            "n_input": len(results) + n_failed + n_skipped + n_resume_skipped,
            "n_trees": n_successful if not (modelfinder == "MF") else 0,
            "n_failed": n_failed,
            "n_skipped": n_skipped,
            "seq_type": resolved_seq_type,
            "workflow": workflow,
            "model_selected": model_selected,
            "log_likelihood": single_log_likelihood,
            "boot": boot if (boot and boot > 0) else None,
            "alrt": alrt,
            "partitioned": partitions is not None,
            "merged_partitions": (
                partitions is not None and modelfinder in ("MF", "MFP")
            ),
        },
        "error": error_msg,
        "data": data_block,
    }

    return payload


def _reconstruct_result(output_dir: Path, run_start: float) -> dict[str, Any]:
    result_path = output_dir / "result.json"
    if result_path.exists():
        return json.loads(result_path.read_text())
    ckpt_path = output_dir / "checkpoint.json"
    if ckpt_path.exists():
        ckpt = json.loads(ckpt_path.read_text())
        cmd = ckpt.get("command", "")
    else:
        cmd = ""
    if not cmd:
        return {
            "status": "error",
            "command": "phyloai tree ml iqtree",
            "wall_time": _time.monotonic() - run_start,
            "tool_versions": {},
            "params": {},
            "key_results": {},
            "error": "Cannot reconstruct result: result.json not found and checkpoint.json missing or has no command",
            "data": {"summary": {}, "files": [], "failed": [], "skipped": [], "warnings": []},
        }
    return {
        "status": "success",
        "command": cmd,
        "wall_time": _time.monotonic() - run_start,
        "tool_versions": {},
        "params": {},
        "key_results": {},
        "error": None,
        "data": {"summary": {}, "files": [], "failed": [], "skipped": [], "warnings": []},
    }
