"""Concordance factor computation (gCF, sCF, sCFl, qCF) via IQ-TREE3 and wASTRAL."""

from __future__ import annotations

import logging
import os
import re as _re
import shutil
import subprocess
import sys
import time as _time
from pathlib import Path
from typing import Any

from phyloai.core.env import ToolEnv
from phyloai.core.iqtree import (
    _resolve_iqtree_path,
    _detect_iqtree_version,
)

_logger = logging.getLogger(__name__)

# Default prefix per CF mode
_DEFAULT_PREFIX: dict[str, str] = {
    "gcf": "gCF",
    "scf": "sCF",
    "scfl": "sCFl",
    "gcf+scf": "gCFsCF",
    "qcf": "qCF",
}

# CF modes that require gene trees
_CF_MODES_NEED_GENE_TREES = frozenset({"gcf", "gcf+scf", "qcf"})

# CF modes that require a matrix
_CF_MODES_NEED_MATRIX = frozenset({"scf", "scfl", "gcf+scf"})

# CF modes that use iqtree3
_CF_MODES_IQTREE = frozenset({"gcf", "scf", "scfl", "gcf+scf"})


# ---- Input scanning ---------------------------------------------------


def _scan_input_cf(
    tree_dir: Path,
) -> tuple[list[Path], list[dict[str, str]]]:
    """Scan a directory for gene tree files.

    Suffix-agnostic (per Section 9.7 policy): every non-empty regular file
    is checked.  A file is accepted when the first line looks like a newick
    tree string (contains `(` and ends with `;`).

    Returns:
        (valid_files, skipped_entries)
    """
    if not tree_dir.exists():
        return [], []

    found: list[Path] = []
    skipped: list[dict[str, str]] = []

    for entry in sorted(tree_dir.iterdir()):
        if entry.is_dir():
            skipped.append({"path": str(entry), "reason": "directory"})
            continue
        if not entry.is_file():
            skipped.append({"path": str(entry), "reason": "not a regular file"})
            continue
        if entry.stat().st_size == 0:
            skipped.append({"path": str(entry), "reason": "empty file"})
            continue

        try:
            content = entry.read_text().strip()
            if not content:
                skipped.append({"path": str(entry), "reason": "empty file content"})
                continue
            if "(" in content and content.rstrip().endswith(";"):
                found.append(entry)
            else:
                skipped.append({"path": str(entry), "reason": "content does not look like newick"})
        except UnicodeDecodeError:
            skipped.append({"path": str(entry), "reason": "binary file"})

    return found, skipped


def _merge_gene_trees(
    tree_dir: Path,
    output_path: Path,
) -> tuple[int, list[dict[str, str]]]:
    """Scan tree_dir for newick files, merge into one file (one tree per line).

    Handles multi-line newick by splitting on tree terminators (``;``).

    Returns:
        (count_of_trees_merged, skipped_entries)
    """
    found, skipped = _scan_input_cf(tree_dir)

    count = 0
    with open(output_path, "w") as out:
        for f in found:
            content = f.read_text().strip()
            if not content:
                continue
            for tree in content.split(";"):
                tree = tree.strip()
                if not tree or not tree.strip():
                    continue
                out.write(tree.rstrip(";").strip() + ";\n")
                count += 1

    return count, skipped


# ---- qCF mapper -------------------------------------------------------


def _map_qcf_to_tree(
    ref_tree_path: Path,
    wastral_output_path: Path,
    output_path: Path,
    *,
    lpp: bool = False,
) -> None:
    """Map qCF values from wastral tree onto reference tree topology.

    Annotates the reference tree by injecting q1 (and optionally pp1) values
    into the raw Newick string at each internal node.  Branch lengths and
    existing support labels are preserved exactly — no Bio.Phylo round-trip
    distorts the original formatting.
    """
    from io import StringIO

    from Bio import Phylo

    ref_raw = ref_tree_path.read_text().strip()
    wastral_raw = wastral_output_path.read_text().strip()

    ref_tree = Phylo.read(StringIO(ref_raw), "newick")
    wastral_tree = Phylo.read(StringIO(wastral_raw), "newick")

    all_leaves = frozenset(
        clade.name for clade in ref_tree.find_clades()
        if clade.is_terminal() and clade.name is not None
    )

    def _fmt_val(v: float) -> str:
        """Format a CF value: 4 decimal places, strip trailing zeros."""
        s = f"{v:.4f}"
        s = s.rstrip("0").rstrip(".")
        return s if s else "0"

    def _leaf_set(clade):
        return frozenset(
            t.name for t in clade.find_clades(terminal=True)
            if t.name is not None
        )

    # ── Raw Newick fallback for wastral value extraction ──
    def _extract_raw_values(raw_nwk: str, target_clade) -> tuple[float | None, float | None]:
        """Extract q1/pp1 from raw wastral Newick for a given clade by post-order position."""
        q1_pattern = _re.compile(r"q1=([0-9.]+)")
        pp1_pattern = _re.compile(r"pp1=([0-9.]+)")

        wastral_clades_postorder = [
            c for c in wastral_tree.find_clades(order="postorder")
            if not c.is_terminal()
        ]

        # Find the target clade's position
        target_idx = None
        target_ls = _leaf_set(target_clade)
        for idx, c in enumerate(wastral_clades_postorder):
            if _leaf_set(c) == target_ls:
                target_idx = idx
                break
        if target_idx is None:
            return None, None

        # Extract internal node labels from raw Newick in post-order
        internal_labels: list[str] = []
        i = 0
        while i < len(raw_nwk):
            if raw_nwk[i] == ')':
                i += 1
                label_start = i
                while i < len(raw_nwk) and raw_nwk[i] not in (':', ',', ')'):
                    i += 1
                label = raw_nwk[label_start:i].strip()
                internal_labels.append(label)
                if i < len(raw_nwk) and raw_nwk[i] == ':':
                    i += 1
                    while i < len(raw_nwk) and raw_nwk[i] not in (',', ')'):
                        i += 1
                continue
            i += 1

        if target_idx >= len(internal_labels):
            return None, None

        label = internal_labels[target_idx]
        q1_val, pp1_val = None, None
        m = q1_pattern.search(label)
        if m:
            q1_val = float(m.group(1))
        if lpp:
            m = pp1_pattern.search(label)
            if m:
                pp1_val = float(m.group(1))
        return q1_val, pp1_val

    # ── Build canonical-bip → (q1_str, pp1_str) from wastral ──
    value_map: dict[frozenset[str], tuple[str, str]] = {}

    def _extract_clade_values(clade) -> tuple[float | None, float | None]:
        q1_val, pp1_val = None, None
        for attr in (clade.name, getattr(clade, "comment", None)):
            if attr and isinstance(attr, str):
                m = _re.search(r"q1=([0-9.]+)", attr)
                if m:
                    q1_val = float(m.group(1))
                if lpp:
                    m = _re.search(r"pp1=([0-9.]+)", attr)
                    if m:
                        pp1_val = float(m.group(1))
        return q1_val, pp1_val

    for clade in wastral_tree.find_clades():
        if clade.is_terminal():
            continue
        q1_val, pp1_val = _extract_clade_values(clade)
        if q1_val is None and (not lpp or pp1_val is None):
            # Fallback: raw Newick
            q1_val, pp1_val = _extract_raw_values(wastral_raw, clade)
        if q1_val is None:
            continue
        ls = _leaf_set(clade)
        complement = all_leaves - ls
        canonical = ls if sorted(ls) < sorted(complement) else complement
        if canonical not in value_map:
            qs = _fmt_val(q1_val)
            ps = _fmt_val(pp1_val) if lpp and pp1_val is not None else ""
            value_map[canonical] = (qs, ps)

    # ── Walk raw ref Newick, inject annotations at each internal node ──
    internal_nodes = [
        c for c in ref_tree.find_clades(order="postorder")
        if not c.is_terminal()
    ]

    buf: list[str] = []
    i = 0
    node_idx = 0
    while i < len(ref_raw):
        c = ref_raw[i]
        if c == ')':
            i += 1

            # Determine which clade this ')' corresponds to
            clade = internal_nodes[node_idx] if node_idx < len(internal_nodes) else None
            node_idx += 1

            # Extract the label between ')' and ':' (or ',' ';' ')')
            label_start = i
            while i < len(ref_raw) and ref_raw[i] not in (':', ',', ';', ')'):
                i += 1
            orig_label = ref_raw[label_start:i]
            next_char = ref_raw[i] if i < len(ref_raw) else ''

            # Compute annotation for this clade
            if clade is not None:
                ls = _leaf_set(clade)
                complement = all_leaves - ls
                canonical = ls if sorted(ls) < sorted(complement) else complement
                v = value_map.get(canonical)
                if v is not None:
                    q1_str, pp1_str = v
                    if lpp and pp1_str:
                        new_label = f"{orig_label}/{q1_str}/{pp1_str}" if orig_label.strip() else f"{q1_str}/{pp1_str}"
                    else:
                        new_label = f"{orig_label}/{q1_str}" if orig_label.strip() else q1_str
                    buf.append(')')
                    buf.append(new_label)
                else:
                    buf.append(')')
                    buf.append(orig_label)
            else:
                buf.append(')')
                buf.append(orig_label)

            # Append branch length and anything after the label
            if next_char == ':':
                i += 1  # skip ':'
                branch_start = i
                while i < len(ref_raw) and ref_raw[i] not in (',', ';', ')'):
                    i += 1
                buf.append(':')
                buf.append(ref_raw[branch_start:i])
        else:
            # Skip quoted strings (wASTRAL-style labels) to avoid misinterpreting ')' inside quotes
            if c in ("'", '"'):
                quote = c
                buf.append(c)
                i += 1
                while i < len(ref_raw) and ref_raw[i] != quote:
                    buf.append(ref_raw[i])
                    i += 1
                if i < len(ref_raw):
                    buf.append(ref_raw[i])
                    i += 1
            else:
                buf.append(c)
                i += 1

    output_path.write_text("".join(buf))


# ---- Command builders -------------------------------------------------


def _build_iqtree_cf_cmd(
    *,
    cf_mode: str,
    executable: str,
    ref_tree: Path,
    gene_trees: Path | None,
    matrix: Path | None,
    scf_quartets: int,
    prefix: str,
    threads: int,
    model_expr: str | None = None,
    partitions: str | None = None,
) -> list[str]:
    """Build the IQ-TREE3 command for a given CF mode."""
    cmd = [executable]

    if cf_mode in ("scf", "scfl", "gcf+scf"):
        assert matrix is not None
        cmd.extend(["-s", str(matrix)])

    if cf_mode in ("gcf", "gcf+scf"):
        cmd.extend(["-t", str(ref_tree)])
    else:
        cmd.extend(["-te", str(ref_tree)])

    if cf_mode in ("gcf", "gcf+scf"):
        assert gene_trees is not None
        cmd.extend(["--gcf", str(gene_trees)])

    if cf_mode in ("scf", "gcf+scf"):
        cmd.extend(["--scf", str(scf_quartets)])
    elif cf_mode == "scfl":
        cmd.extend(["--scfl", str(scf_quartets)])

    if cf_mode == "scfl":
        if partitions is not None:
            cmd.extend(["-p", partitions])
        elif model_expr is not None:
            cmd.extend(["-m", model_expr])

    cmd.extend(["--prefix", prefix])
    cmd.extend(["-T", str(threads)])

    return cmd


def _build_wastral_qcf_cmd(
    *,
    executable: str,
    gene_trees: Path,
    ref_tree: Path,
    output_dir: Path,
    threads: int,
) -> list[str]:
    """Build the wASTRAL command for qCF computation.

    Uses absolute paths for all inputs; output uses 'wastral.tre' (relative
    to cwd) since the subprocess runs with cwd=output_dir.
    """
    return [
        executable,
        "-i", str(gene_trees.resolve()),
        "-o", "wastral.tre",
        "-u", "2",
        "-c", str(ref_tree.resolve()),
        "-C",
        "--mode", "4",
        "-t", str(threads),
    ]


def _resolve_wastral_path(wastral_path: str | None, dry_run: bool) -> str:
    """Resolve wastral executable path (always returns absolute)."""
    if wastral_path:
        p = Path(wastral_path).resolve()
        if not p.exists():
            raise ValueError(f"--wastral-path does not exist: {wastral_path}")
        if not os.access(str(p), os.X_OK):
            raise ValueError(f"--wastral-path is not executable: {wastral_path}")
        return str(p)
    if dry_run:
        return "wastral"
    for name in ("wastral", "aster"):
        env = ToolEnv()
        try:
            return str(env.require(name))
        except FileNotFoundError:
            continue
    raise FileNotFoundError(
        "wastral not found. Install from https://github.com/chaoszhang/ASTER "
        "or use --wastral-path."
    )


def _detect_wastral_version(executable: str) -> dict[str, str]:
    """Detect wastral version via -v or -h fallback."""
    combined = ""
    for flag in ("-v", "-h"):
        try:
            proc = subprocess.run(
                [executable, flag],
                capture_output=True, text=True, timeout=10,
            )
            combined = proc.stdout + "\n" + proc.stderr
            if proc.returncode == 0 or (proc.stdout.strip() or proc.stderr.strip()):
                break
        except Exception:
            continue

    if not combined.strip():
        return {"wastral": "unknown"}

    m = _re.search(r"Version:\s*v?([\d.]+)", combined, _re.IGNORECASE)
    if m:
        return {"wastral": m.group(1)}
    m = _re.search(r"version\s+([\d.]+)", combined, _re.IGNORECASE)
    if m:
        return {"wastral": m.group(1)}
    m = _re.search(r"\bv?(\d+(?:\.\d+)+)\b", combined)
    if m:
        return {"wastral": m.group(1)}
    return {"wastral": "unknown"}


# ---- Logging helper --------------------------------------------------


# ---- Result assembly --------------------------------------------------


def _assemble_cf_result(
    *,
    run_start: float,
    cf_mode: str,
    ref_tree: Path,
    tree: Path | None,
    tree_dir: Path | None,
    matrix: Path | None,
    partitions: Path | None,
    model_expr: str | None,
    scf_quartets: int,
    prefix: str,
    output_dir: Path,
    threads: int,
    iqtree_path: str | None,
    wastral_path: str | None,
    overwrite: bool,
    dry_run: bool,
    quiet: bool = False,
    lpp: bool,
    input_path: Path,
    n_input_trees: int,
    cmd: list[str],
    wall_time: float,
    skipped: list[dict[str, str]],
    warnings_list: list[str],
    is_error: bool,
    error_msg: str | None,
    versions: dict[str, str],
    iqtree_exe: str | None = None,
    wastral_exe: str | None = None,
    tool_stderr: str = "",
) -> dict[str, Any]:
    """Build the result.json payload."""
    if tree_dir is not None:
        input_mode = "--tree-dir"
    elif tree is not None:
        input_mode = "--tree"
    else:
        input_mode = "--matrix"

    cmd_parts = ["phyloai", "tree", "cf", "--cf", cf_mode]
    cmd_parts.extend(["--ref-tree", str(ref_tree)])
    if tree is not None:
        cmd_parts.extend(["--tree", str(tree)])
    elif tree_dir is not None:
        cmd_parts.extend(["--tree-dir", str(tree_dir)])
    if matrix is not None:
        cmd_parts.extend(["--matrix", str(matrix)])
    if partitions is not None:
        cmd_parts.extend(["--partitions", str(partitions)])
    if model_expr is not None:
        cmd_parts.extend(["--model-expr", model_expr])
    if cf_mode not in ("gcf", "qcf"):
        cmd_parts.extend(["--scf-quartets", str(scf_quartets)])
    if lpp:
        cmd_parts.append("--lpp")
    cmd_parts.extend(["--prefix", prefix])
    cmd_parts.extend(["-o", str(output_dir)])
    cmd_parts.extend(["-t", str(threads)])
    if overwrite:
        cmd_parts.append("--overwrite")
    if iqtree_path:
        cmd_parts.extend(["--iqtree-path", iqtree_path])
    if wastral_path:
        cmd_parts.extend(["--wastral-path", wastral_path])
    if dry_run:
        cmd_parts.append("--dry-run")
    cmd_str = " ".join(cmd_parts)

    input_data: dict[str, Any] = {"path": str(input_path)}
    if input_mode == "--tree-dir":
        input_data["n_trees"] = n_input_trees

    key_results: dict[str, Any] = {
        "cf_type": cf_mode,
        "prefix": prefix,
    }
    if cf_mode in _CF_MODES_IQTREE:
        key_results["cf_stat"] = str(output_dir / f"{prefix}.cf.stat")
        key_results["cf_tree"] = str(output_dir / f"{prefix}.cf.tree")
    else:
        key_results["cf_tree"] = str(output_dir / f"{prefix}.cf.tree")

    return {
        "status": "error" if is_error else "success",
        "command": cmd_str,
        "wall_time": _time.monotonic() - run_start,
        "tool_versions": versions,
        "params": {
            "cf": cf_mode,
            "ref_tree": str(ref_tree),
            "tree": str(tree) if tree else None,
            "tree_dir": str(tree_dir) if tree_dir else None,
            "matrix": str(matrix) if matrix else None,
            "partitions": str(partitions) if partitions else None,
            "model_expr": model_expr,
            "scf_quartets": scf_quartets if cf_mode not in ("gcf", "qcf") else None,
            "lpp": lpp,
            "prefix": prefix,
            "output_dir": str(output_dir),
            "threads": threads,
            "overwrite": overwrite,
            "dry_run": dry_run,
            "iqtree_path": iqtree_path,
            "wastral_path": wastral_path,
            "quiet": quiet,
        },
        "key_results": key_results,
        "error": error_msg,
        "data": {
            "input_mode": input_mode,
            "input": input_data,
            "cmd": cmd,
            "tool_stderr": tool_stderr,
            **({"tool_log": f"{prefix}.log"} if cf_mode in _CF_MODES_IQTREE else {"tool_log": "wastral.log"} if cf_mode == "qcf" else {}),
            "skipped": skipped,
            "warnings": warnings_list,
        },
    }


# ---- Entry point ------------------------------------------------------


def run_cf(
    *,
    cf_mode: str,
    ref_tree: Path,
    tree: Path | None = None,
    tree_dir: Path | None = None,
    matrix: Path | None = None,
    partitions: Path | None = None,
    model_expr: str | None = None,
    scf_quartets: int = 100,
    prefix: str | None = None,
    output_dir: Path = Path("runs/tree/cf"),
    threads: int = 4,
    iqtree_path: str | None = None,
    wastral_path: str | None = None,
    overwrite: bool = False,
    dry_run: bool = False,
    quiet: bool = False,
    lpp: bool = False,
    **kwargs: Any,
) -> dict[str, Any]:
    """Run concordance factor computation.

    Returns a result.json-compatible payload dict.
    Raises ValueError for invalid inputs.
    Raises FileNotFoundError for missing tools.
    """
    if "model" in kwargs:
        import warnings
        warnings.warn(
            "run_cf(model=...) is deprecated, use model_expr= instead",
            DeprecationWarning, stacklevel=2,
        )
        if model_expr is None:
            model_expr = kwargs["model"]
        del kwargs["model"]
    if kwargs:
        raise TypeError(f"run_cf() got unexpected keyword arguments: {set(kwargs)}")
    # --- Validate cf_mode ---
    valid_modes = frozenset({"gcf", "scf", "scfl", "gcf+scf", "qcf"})
    if cf_mode not in valid_modes:
        raise ValueError(
            f"Invalid --cf mode: {cf_mode}. Valid: {', '.join(sorted(valid_modes))}"
        )

    if lpp and cf_mode != "qcf":
        raise ValueError("--lpp is only valid with --cf qcf.")

    # --- Resolve prefix ---
    if prefix is None:
        prefix = _DEFAULT_PREFIX[cf_mode]

    # --- Validate ref_tree ---
    if not ref_tree.exists():
        raise ValueError(f"--ref-tree does not exist: {ref_tree}")
    ref_tree = ref_tree.resolve()

    # --- Validate gene trees (tree / tree-dir) ---
    needs_gene_trees = cf_mode in _CF_MODES_NEED_GENE_TREES
    if needs_gene_trees:
        if (tree is None and tree_dir is None) or (tree is not None and tree_dir is not None):
            raise ValueError(
                "--tree or --tree-dir must be provided (mutually exclusive) "
                f"for --cf {cf_mode}."
            )
        if tree is not None and not tree.exists():
            raise ValueError(f"--tree does not exist: {tree}")
        if tree_dir is not None and not tree_dir.exists():
            raise ValueError(f"--tree-dir does not exist: {tree_dir}")
    else:
        if tree is not None or tree_dir is not None:
            raise ValueError(
                f"--tree/--tree-dir is not needed for --cf {cf_mode}."
            )

    # --- Validate matrix ---
    needs_matrix = cf_mode in _CF_MODES_NEED_MATRIX
    if needs_matrix:
        if matrix is None:
            raise ValueError(f"--matrix is required for --cf {cf_mode}.")
        if not matrix.exists():
            raise ValueError(f"--matrix does not exist: {matrix}")
        matrix = matrix.resolve()
    else:
        if matrix is not None:
            raise ValueError(
                f"--matrix is not valid for --cf {cf_mode}."
            )

    # --- Validate scfl-only params ---
    if cf_mode != "scfl":
        if model_expr is not None:
            raise ValueError(f"--model-expr is not valid for --cf {cf_mode}.")
        if partitions is not None:
            raise ValueError(f"--partitions is not valid for --cf {cf_mode}.")

    if cf_mode == "scfl":
        if model_expr is not None and partitions is not None:
            raise ValueError(
                "--model-expr and --partitions are mutually exclusive for --cf scfl."
            )
        if partitions is not None:
            if not partitions.exists():
                raise ValueError(f"--partitions does not exist: {partitions}")
            partitions = partitions.resolve()

    # --- Validate scf_quartets ---
    if cf_mode in ("gcf", "qcf"):
        if scf_quartets != 100:
            raise ValueError(
                f"--scf-quartets is not valid for --cf {cf_mode}."
            )
    else:
        if scf_quartets < 1:
            raise ValueError(f"--scf-quartets must be >= 1. Got: {scf_quartets}")

    if threads < 1:
        raise ValueError(f"--threads must be >= 1. Got: {threads}")

    # --- Output directory conflict ---
    if not dry_run:
        if overwrite and output_dir.exists():
            shutil.rmtree(output_dir)
        if not overwrite and output_dir.exists() and any(output_dir.iterdir()):
            raise ValueError(
                f"Output directory {output_dir} already exists and is non-empty. "
                "Use --overwrite to replace."
            )

    run_start = _time.monotonic()
    skipped: list[dict[str, str]] = []
    warnings_list: list[str] = []
    n_input_trees = 0

    # --- Input resolution ---
    if tree is not None:
        input_path = tree.resolve()
    elif tree_dir is not None:
        valid_files, scanned_skipped = _scan_input_cf(tree_dir)
        skipped = scanned_skipped
        n_valid_files = len(valid_files)

        if n_valid_files == 0:
            raise ValueError(
                f"No valid gene tree files found in --tree-dir: {tree_dir}"
            )

        if n_valid_files == 1:
            msg = (
                "Exactly 1 valid gene tree file in --tree-dir. "
                "Consider using --tree mode directly."
            )
            warnings_list.append(msg)
            if not quiet:
                _logger.warning(msg)

        unrecognized = [
            s for s in skipped
            if "newick" in s.get("reason", "") or s.get("reason", "").startswith("binary")
        ]
        if unrecognized:
            msg = (
                f"--tree-dir contains {len(unrecognized)} non-newick file(s); "
                "skipped. See result.json data.skipped for details."
            )
            warnings_list.append(msg)
            if not quiet:
                _logger.warning(msg)

        if dry_run:
            merged_path = output_dir / "merged.trees"
            n_input_trees = 0
            for f in valid_files:
                content = f.read_text().strip()
                if content:
                    n_input_trees += len([line for line in content.splitlines() if line.strip()])
            # merged.trees may not exist yet — use safe absolute path
            input_path = Path(os.path.abspath(str(merged_path)))
        else:
            output_dir.mkdir(parents=True, exist_ok=True)
            merged_path = output_dir / "merged.trees"
            n_input_trees, _ = _merge_gene_trees(tree_dir, merged_path)
            input_path = merged_path.resolve()
    else:
        input_path = matrix  # type: ignore[assignment]  (already resolved above)

    # --- Warnings ---
    if cf_mode not in ("gcf", "qcf") and scf_quartets < 100:
        msg = f"--scf-quartets is {scf_quartets}; recommend >= 100 for reliable results."
        warnings_list.append(msg)
        if not quiet:
            _logger.warning(msg)

    if cf_mode == "scfl" and model_expr is None and partitions is None:
        msg = (
            "--cf scfl without --model-expr or --partitions: IQ-TREE3 will "
            "auto-compute the best-fit model (slow). Consider providing "
            "--model-expr or --partitions for speedup."
        )
        warnings_list.append(msg)
        if not quiet:
            _logger.warning(msg)

    # --- Resolve executables ---
    if cf_mode in _CF_MODES_IQTREE:
        iqtree_exe = _resolve_iqtree_path(iqtree_path, dry_run)
        wastral_exe = None
    else:
        iqtree_exe = None
        wastral_exe = _resolve_wastral_path(wastral_path, dry_run)

    # --- Command building ---
    if cf_mode in _CF_MODES_IQTREE:
        cmd = _build_iqtree_cf_cmd(
            cf_mode=cf_mode,
            executable=iqtree_exe,  # type: ignore[arg-type]
            ref_tree=ref_tree,
            gene_trees=input_path if cf_mode in ("gcf", "gcf+scf") else None,
            matrix=matrix,
            scf_quartets=scf_quartets,
            prefix=prefix,
            threads=threads,
            model_expr=model_expr,
            partitions=str(partitions) if partitions else None,
        )
    else:
        cmd = _build_wastral_qcf_cmd(
            executable=wastral_exe,  # type: ignore[arg-type]
            gene_trees=input_path,
            ref_tree=ref_tree,
            output_dir=output_dir,
            threads=threads,
        )

    # --- Dry run: return payload without execution ---
    if dry_run:
        versions = {}
        if cf_mode in _CF_MODES_IQTREE:
            versions = {"iqtree3": "unknown"}
        else:
            versions = {"wastral": "unknown"}
        return _assemble_cf_result(
            run_start=run_start,
            cf_mode=cf_mode, ref_tree=ref_tree,
            tree=tree, tree_dir=tree_dir,
            matrix=matrix, partitions=partitions,
            model_expr=model_expr, scf_quartets=scf_quartets,
            prefix=prefix, output_dir=output_dir,
            threads=threads,
            iqtree_path=iqtree_path, wastral_path=wastral_path,
            overwrite=overwrite, dry_run=dry_run,
            quiet=quiet,
            lpp=lpp,
            input_path=input_path, n_input_trees=n_input_trees,
            cmd=cmd, wall_time=0.0,
            skipped=skipped, warnings_list=warnings_list,
            is_error=False, error_msg=None,
            versions=versions,
            tool_stderr="",
        )

    # --- Execution ---
    output_dir.mkdir(parents=True, exist_ok=True)

    versions: dict[str, str] = {}
    try:
        proc = subprocess.Popen(
            cmd,
            cwd=str(output_dir),
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
        # Stream output lines to terminal in real-time, also capture for diagnostics
        captured_stdout: list[str] = []
        captured_stderr: list[str] = []
        assert proc.stdout is not None
        assert proc.stderr is not None

        import select

        _stdout_fd = proc.stdout.fileno()
        _stderr_fd = proc.stderr.fileno()

        while True:
            rlist, _, _ = select.select([_stdout_fd, _stderr_fd], [], [], 0.1)
            done_reading = False
            for fd in rlist:
                if fd == _stdout_fd:
                    line = proc.stdout.readline()
                    if line:
                        captured_stdout.append(line)
                        if not quiet:
                            sys.stderr.write(line)
                    else:
                        done_reading = True
                elif fd == _stderr_fd:
                    line = proc.stderr.readline()
                    if line:
                        captured_stderr.append(line)
                        if not quiet:
                            sys.stderr.write(line)
                    else:
                        done_reading = True
            if done_reading and proc.poll() is not None:
                # Drain any remaining
                for remaining in proc.stdout.readlines():
                    captured_stdout.append(remaining)
                    if not quiet:
                        sys.stderr.write(remaining)
                for remaining in proc.stderr.readlines():
                    captured_stderr.append(remaining)
                    if not quiet:
                        sys.stderr.write(remaining)
                break

        proc.wait()
        proc.stdout = "".join(captured_stdout)  # type: ignore[assignment]
        proc.stderr = "".join(captured_stderr)  # type: ignore[assignment]
        if cf_mode == "qcf":
            (output_dir / "wastral.log").write_text(str(proc.stderr))
    except Exception as exc:
        if cf_mode in _CF_MODES_IQTREE:
            versions = {"iqtree3": "unknown"}
        else:
            versions = {"wastral": "unknown"}
        return _assemble_cf_result(
            run_start=run_start,
            cf_mode=cf_mode, ref_tree=ref_tree,
            tree=tree, tree_dir=tree_dir,
            matrix=matrix, partitions=partitions,
            model_expr=model_expr, scf_quartets=scf_quartets,
            prefix=prefix, output_dir=output_dir,
            threads=threads,
            iqtree_path=iqtree_path, wastral_path=wastral_path,
            overwrite=overwrite, dry_run=dry_run,
            quiet=quiet,
            lpp=lpp,
            input_path=input_path, n_input_trees=n_input_trees,
            cmd=cmd, wall_time=0.0,
            skipped=skipped, warnings_list=warnings_list,
            is_error=True, error_msg=str(exc),
            versions=versions,
            tool_stderr="",
        )

    if cf_mode in _CF_MODES_IQTREE:
        versions = _detect_iqtree_version(iqtree_exe)  # type: ignore[arg-type]
    else:
        versions = _detect_wastral_version(wastral_exe)  # type: ignore[arg-type]

    wall_time = _time.monotonic() - run_start

    if proc.returncode != 0:
        error_msg = (
            f"IQ-TREE3 exited with code {proc.returncode}"
            if cf_mode in _CF_MODES_IQTREE
            else f"wASTRAL exited with code {proc.returncode}"
        )
        error_msg += f": {str(proc.stderr)[:500]}"
        return _assemble_cf_result(
            run_start=run_start,
            cf_mode=cf_mode, ref_tree=ref_tree,
            tree=tree, tree_dir=tree_dir,
            matrix=matrix, partitions=partitions,
            model_expr=model_expr, scf_quartets=scf_quartets,
            prefix=prefix, output_dir=output_dir,
            threads=threads,
            iqtree_path=iqtree_path, wastral_path=wastral_path,
            overwrite=overwrite, dry_run=dry_run,
            quiet=quiet,
            lpp=lpp,
            input_path=input_path, n_input_trees=n_input_trees,
            cmd=cmd, wall_time=wall_time,
            skipped=skipped, warnings_list=warnings_list,
            is_error=True, error_msg=error_msg,
            versions=versions,
            tool_stderr=str(proc.stderr),
        )

    # Post-process for qCF: map values from wastral.tre to ref tree
    if cf_mode == "qcf":
        wastral_tre_path = output_dir / "wastral.tre"

        if wastral_tre_path.exists():
            try:
                _map_qcf_to_tree(ref_tree, wastral_tre_path, output_dir / f"{prefix}.cf.tree", lpp=lpp)
            except Exception as exc:
                error_msg = f"qCF mapping failed: {exc}"
                return _assemble_cf_result(
                    run_start=run_start,
                    cf_mode=cf_mode, ref_tree=ref_tree,
                    tree=tree, tree_dir=tree_dir,
                    matrix=matrix, partitions=partitions,
                    model_expr=model_expr, scf_quartets=scf_quartets,
                    prefix=prefix, output_dir=output_dir,
                    threads=threads,
                    iqtree_path=iqtree_path, wastral_path=wastral_path,
                    overwrite=overwrite, dry_run=dry_run,
                    quiet=quiet,
                    lpp=lpp,
                    input_path=input_path, n_input_trees=n_input_trees,
                    cmd=cmd, wall_time=wall_time,
                    skipped=skipped, warnings_list=warnings_list,
                    is_error=True, error_msg=error_msg,
                    versions=versions,
                    tool_stderr=str(proc.stderr),
                )
        else:
            error_msg = "wASTRAL completed but did not produce wastral.tre"
            return _assemble_cf_result(
                run_start=run_start,
                cf_mode=cf_mode, ref_tree=ref_tree,
                tree=tree, tree_dir=tree_dir,
                matrix=matrix, partitions=partitions,
                model_expr=model_expr, scf_quartets=scf_quartets,
                prefix=prefix, output_dir=output_dir,
                threads=threads,
                iqtree_path=iqtree_path, wastral_path=wastral_path,
                overwrite=overwrite, dry_run=dry_run,
                quiet=quiet,
                lpp=lpp,
                input_path=input_path, n_input_trees=n_input_trees,
                cmd=cmd, wall_time=wall_time,
                skipped=skipped, warnings_list=warnings_list,
                is_error=True, error_msg=error_msg,
                versions=versions,
                tool_stderr=str(proc.stderr),
            )

    payload = _assemble_cf_result(
        run_start=run_start,
        cf_mode=cf_mode, ref_tree=ref_tree,
        tree=tree, tree_dir=tree_dir,
        matrix=matrix, partitions=partitions,
        model_expr=model_expr, scf_quartets=scf_quartets,
        prefix=prefix, output_dir=output_dir,
        threads=threads,
        iqtree_path=iqtree_path, wastral_path=wastral_path,
        overwrite=overwrite, dry_run=dry_run,
        quiet=quiet,
        lpp=lpp,
        input_path=input_path, n_input_trees=n_input_trees,
        cmd=cmd, wall_time=wall_time,
        skipped=skipped, warnings_list=warnings_list,
        is_error=False, error_msg=None,
        versions=versions,
        iqtree_exe=iqtree_exe,
        wastral_exe=wastral_exe,
        tool_stderr=str(proc.stderr),
    )

    return payload
