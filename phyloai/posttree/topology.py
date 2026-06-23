"""Tree topology tests (AU / KH / SH / WKH / WSH / c-ELW) via IQ-TREE."""
from __future__ import annotations

import shlex
import subprocess
import time as _time
from pathlib import Path
from typing import Any

from phyloai.core.iqtree import (
    _detect_iqtree_version,
    _resolve_iqtree_path,
    IQTREE_COMPATIBLE_EXTENSIONS,
)


# ===================================================================
# Library-layer input validation
# ===================================================================

def _validate_inputs(
    *,
    matrix: Path,
    candidate_trees: list[Path],
    replicates: int,
    threads: int,
    overwrite: bool,
    resume: bool,
    model_expr: str | None,
    partitions: str | None,
    tool_args: str | None,
    guide_tree: str | None,
) -> list[str]:
    """Validate all spec-mandated inputs.  Returns a list of error messages.

    Covers: matrix existence, candidate-tree existence / non-empty / readable,
    numeric bounds, parameter mutual-exclusion, model-source completeness.
    Output-directory lifecycle (conflict / delete / resume-dir) is left to
    the CLI layer.
    """
    errors: list[str] = []

    import os as _os

    if not matrix.exists():
        errors.append(f"--matrix does not exist: {matrix}")
    elif not matrix.is_file():
        errors.append(f"--matrix is not a regular file: {matrix}")
    elif not _os.access(str(matrix), _os.R_OK):
        errors.append(f"--matrix is not readable: {matrix}")

    if not candidate_trees:
        errors.append("At least one --candidate-trees value is required")
    else:
        for i, ct in enumerate(candidate_trees):
            if not ct.exists():
                errors.append(f"--candidate-trees #{i + 1} does not exist: {ct}")
            elif not ct.is_file():
                errors.append(f"--candidate-trees #{i + 1} is not a regular file: {ct}")
            elif ct.stat().st_size == 0:
                errors.append(f"--candidate-trees #{i + 1} is empty: {ct}")
            elif not _os.access(str(ct), _os.R_OK):
                errors.append(f"--candidate-trees #{i + 1} is not readable: {ct}")

    if replicates < 1000:
        errors.append(f"--replicates must be >= 1000, got {replicates}")

    if threads < 1:
        errors.append(f"--threads must be >= 1, got {threads}")

    if overwrite and resume:
        errors.append("--overwrite and --resume are mutually exclusive")

    has_explicit_model = model_expr is not None or partitions is not None
    has_tool_args_model = False
    if tool_args:
        tokens = shlex.split(tool_args)
        has_tool_args_model = "-m" in tokens or "-p" in tokens
    if not has_explicit_model and not has_tool_args_model:
        errors.append(
            "Neither --model-expr, --partitions, nor -m/-p in --tool-args provided. "
            "Must specify one model source."
        )
    if model_expr and partitions:
        errors.append("--model-expr and --partitions are mutually exclusive")

    # Cross-source model conflict: if a high-level model-source flag is
    # given, --tool-args must not contain the OTHER model-source flag,
    # otherwise the effective command would have both -m and -p.
    if tool_args:
        tokens = shlex.split(tool_args)
        if model_expr and "-p" in tokens:
            errors.append(
                "--model-expr is set but --tool-args contains -p. "
                "Remove --model-expr if you want -p from --tool-args to take effect."
            )
        if partitions and "-m" in tokens:
            errors.append(
                "--partitions is set but --tool-args contains -m. "
                "Remove --partitions if you want -m from --tool-args to take effect."
            )

    if guide_tree:
        gt = Path(guide_tree)
        if not gt.exists():
            errors.append(f"--guide-tree does not exist: {guide_tree}")
        elif not gt.is_file():
            errors.append(f"--guide-tree is not a regular file: {guide_tree}")
        elif not _os.access(str(gt), _os.R_OK):
            errors.append(f"--guide-tree is not readable: {guide_tree}")

    return errors


# ===================================================================
# Candidate tree merging
# ===================================================================

def _merge_candidate_trees(
    tree_files: list[Path],
    output_dir: Path,
) -> Path:
    """Merge multiple NEWICK tree files into one candidate.trees file.

    Each input file should contain one NEWICK tree. Trees are concatenated
    in order, one per line.
    """
    merged_path = output_dir / "candidate.trees"
    with open(merged_path, "w") as out:
        for tf in tree_files:
            text = tf.read_text().strip()
            if not text:
                raise ValueError(f"Empty candidate tree file: {tf}")
            out.write(text)
            if not text.endswith("\n"):
                out.write("\n")
    return merged_path


# ===================================================================
# Command builder
# ===================================================================

_TOPOLOGY_BLOCKED_FLAGS = frozenset({"-s", "-z"})
_TOPOLOGY_BLOCKED_IO_CHARS = frozenset({"<", ">", "|"})


def _check_managed_flag_conflict(
    tool_args: str,
    *,
    blocked_flags: frozenset[str],
) -> None:
    """Reject BLOCKED flags and I/O redirects in --tool-args."""
    tokens = shlex.split(tool_args)
    for token in tokens:
        if token in blocked_flags:
            raise ValueError(f"Blocked managed flag in --tool-args: {token}")
        if any(c in token for c in _TOPOLOGY_BLOCKED_IO_CHARS):
            raise ValueError(f"Blocked I/O override in --tool-args: {token}")


def _is_flag_overridden(flag: str, tool_tokens: set[str]) -> bool:
    return flag in tool_tokens


def _build_topology_cmd(
    *,
    executable: str,
    matrix: Path,
    candidate_trees: Path,
    prefix: str,
    model_expr: str | None,
    partitions: str | None,
    guide_tree: str | None,
    replicates: int,
    threads: int,
    tool_args: str | None,
) -> list[str]:
    """Build the IQ-TREE topology test command line.

    Assembly order per design spec section 6.3:
    1. executable, -s <matrix>, -z <candidate-trees>
    2. --prefix (suppress-if-present)
    3. model source: -m or -p (suppress-if-present)
    4. -ft <guide-tree> (suppress-if-present)
    5. topology-test defaults: -n 0 -zb <replicates> -zw -au (each suppressible)
    6. -T <threads> (suppress-if-present)
    7. raw --tool-args appended last
    """
    cmd = [executable]
    tool_tokens = set(shlex.split(tool_args)) if tool_args else set()

    if tool_args:
        _check_managed_flag_conflict(
            tool_args, blocked_flags=_TOPOLOGY_BLOCKED_FLAGS,
        )

    # 1. Required inputs (never overridable)
    cmd.extend(["-s", str(matrix)])
    cmd.extend(["-z", str(candidate_trees)])

    # 2. Prefix
    if not _is_flag_overridden("--prefix", tool_tokens):
        cmd.extend(["--prefix", prefix])

    # 3. Model source
    if model_expr and not _is_flag_overridden("-m", tool_tokens):
        cmd.extend(["-m", model_expr])
    elif partitions and not _is_flag_overridden("-p", tool_tokens):
        cmd.extend(["-p", partitions])

    # 4. Guide tree
    if guide_tree and not _is_flag_overridden("-ft", tool_tokens):
        cmd.extend(["-ft", guide_tree])

    # 5. Topology test defaults
    if not _is_flag_overridden("-n", tool_tokens):
        cmd.extend(["-n", "0"])
    if not _is_flag_overridden("-zb", tool_tokens):
        cmd.extend(["-zb", str(replicates)])
    if not _is_flag_overridden("-zw", tool_tokens):
        cmd.append("-zw")
    if not _is_flag_overridden("-au", tool_tokens):
        cmd.append("-au")

    # 6. Threads
    if not _is_flag_overridden("-T", tool_tokens):
        cmd.extend(["-T", str(threads)])

    # 7. Raw tool-args
    if tool_args:
        cmd.extend(shlex.split(tool_args))

    return cmd


# ===================================================================
# USER TREES table parser (token-based, tolerant across IQ-TREE versions)
# ===================================================================

# Ordered list of columns that IQ-TREE typically emits in the USER TREES
# section.  Each entry is (canonical_name, header_token, takes_sign).
# The first three columns (Tree, logL, deltaL) are always scalars.
# Subsequent test columns carry an optional sign token (+/-) in the data row.
_USER_TREES_COLUMNS = [
    ("tree_id", "Tree", False),
    ("log_likelihood", "logL", False),
    ("delta_likelihood", "deltaL", False),
    ("bp_rell", "bp-RELL", True),
    ("p_kh", "p-KH", True),
    ("p_sh", "p-SH", True),
    ("p_wkh", "p-WKH", True),
    ("p_wsh", "p-WSH", True),
    ("c_elw", "c-ELW", True),
    ("p_au", "p-AU", True),
]


def _parse_user_trees_table(iqtree_path: Path) -> tuple[list[dict[str, Any]], list[str]]:
    """Parse the USER TREES section of an .iqtree report.

    Returns (tests, warnings).  Each test dict maps column names to parsed
    values; absent columns are stored as None.  Each test dict also contains
    a ``raw_line`` key with the original text of the data row.
    """
    if not iqtree_path.exists():
        return [], [f".iqtree file not found: {iqtree_path}"]

    text = iqtree_path.read_text()
    idx = text.find("USER TREES")
    if idx < 0:
        return [], [".iqtree file contains no USER TREES section"]

    lines = text[idx:].splitlines()

    # Find header row — the first line containing "Tree", "logL", "deltaL"
    header_idx = None
    for i, line in enumerate(lines):
        if "Tree" in line and "logL" in line and "deltaL" in line:
            header_idx = i
            break
    if header_idx is None:
        return [], ["Could not find USER TREES table header in .iqtree"]

    header_line = lines[header_idx]

    # Determine which columns are present by scanning header tokens in order.
    # Each column that appears in the header gets an index into the token stream
    # (scalar columns = 1 token; sign-bearing columns = 2 tokens).
    active_columns: list[tuple[str, bool, int]] = []  # (name, takes_sign, header_pos)
    for canonical, token, takes_sign in _USER_TREES_COLUMNS:
        pos = header_line.find(token)
        if pos >= 0:
            active_columns.append((canonical, takes_sign, pos))

    if not active_columns:
        return [], ["No known column tokens found in USER TREES header"]

    # Sort active columns by header position to match output order
    active_columns.sort(key=lambda x: x[2])

    # Parse data rows
    tests: list[dict[str, Any]] = []
    warnings: list[str] = []

    for line in lines[header_idx + 1:]:
        stripped = line.strip()
        if not stripped:
            break  # blank line ends the table
        if all(c in "-= " for c in stripped):
            continue  # separator line

        tokens = stripped.split()
        row: dict[str, Any] = {"raw_line": line}

        # Pre-fill with None
        for canonical, _token, _sign in _USER_TREES_COLUMNS:
            row[canonical] = None
            if canonical not in ("tree_id", "log_likelihood", "delta_likelihood"):
                row[canonical + "_sign"] = None

        # Walk active columns consuming tokens from the data line
        ti = 0  # token index
        for canonical, takes_sign, _hp in active_columns:
            if ti >= len(tokens):
                break

            if canonical == "tree_id":
                try:
                    row[canonical] = int(tokens[ti])
                except ValueError:
                    row[canonical] = None
                ti += 1
            elif canonical in ("log_likelihood", "delta_likelihood"):
                try:
                    row[canonical] = float(tokens[ti])
                except ValueError:
                    row[canonical] = None
                ti += 1
            elif takes_sign:
                try:
                    row[canonical] = float(tokens[ti])
                except ValueError:
                    row[canonical] = None
                ti += 1
                if ti < len(tokens) and tokens[ti] in ("+", "-"):
                    row[canonical + "_sign"] = tokens[ti]
                    ti += 1
            else:
                try:
                    row[canonical] = float(tokens[ti])
                except ValueError:
                    row[canonical] = None
                ti += 1

        tests.append(row)

    return tests, warnings


# ===================================================================
# Output-file discovery
# ===================================================================

def _discover_optimized_trees(output_dir: Path, prefix: str) -> str | None:
    """Find IQ-TREE-optimized candidate tree file by prefix and known suffix patterns."""
    for suffix in (".treels.trees", ".trees", ".treels"):
        path = output_dir / f"{prefix}{suffix}"
        if path.exists():
            return str(path.name)
    return None


# ===================================================================
# result.json assembly
# ===================================================================

def _assemble_topology_result(
    *,
    run_start: float,
    iqtree_exe: str,
    tool_versions: dict[str, str],
    params: dict[str, Any],
    matrix: Path,
    candidate_trees_raw: list[Path],
    candidate_trees_mode: str,
    candidate_trees_effective: Path,
    effective_cmd: list[str],
    tool_stderr: str,
    output_dir: Path,
    prefix: str,
    returncode: int,
    dry_run: bool,
    warnings: list[str],
    error_category: str | None = None,
) -> dict[str, Any]:
    """Assemble the result.json payload.

    error_category: None for success, 'input' for validation errors,
    'env' for missing executable, 'tool' for IQ-TREE execution failure.
    """
    wall_time = round(_time.time() - run_start, 2)

    # Reconstruct CLI command string from resolved params
    cmd_parts = ["phyloai", "posttree", "topology"]
    cmd_parts.extend(["--matrix", str(matrix)])
    for ct in candidate_trees_raw:
        cmd_parts.extend(["--candidate-trees", str(ct)])
    if params.get("input_format") and params["input_format"] != "auto":
        cmd_parts.extend(["--input-format", params["input_format"]])
    if params.get("model_expr"):
        cmd_parts.extend(["--model-expr", params["model_expr"]])
    if params.get("partitions"):
        cmd_parts.extend(["--partitions", params["partitions"]])
    if params.get("guide_tree"):
        cmd_parts.extend(["--guide-tree", params["guide_tree"]])
    if params.get("replicates") != 10000:
        cmd_parts.extend(["--replicates", str(params["replicates"])])
    if params.get("prefix") and params["prefix"] != matrix.stem:
        cmd_parts.extend(["--prefix", params["prefix"]])
    cmd_parts.extend(["-o", str(output_dir)])
    if params.get("threads") != 4:
        cmd_parts.extend(["-t", str(params["threads"])])
    if params.get("iqtree_path"):
        cmd_parts.extend(["--iqtree-path", params["iqtree_path"]])
    if params.get("tool_args"):
        cmd_parts.extend(["--tool-args", params["tool_args"]])
    if params.get("overwrite"):
        cmd_parts.append("--overwrite")
    if params.get("resume"):
        cmd_parts.append("--resume")
    if params.get("dry_run"):
        cmd_parts.append("--dry-run")
    if params.get("quiet"):
        cmd_parts.append("-q")
    command = " ".join(cmd_parts)

    iqtree_log = output_dir / f"{prefix}.log"
    iqtree_report = output_dir / f"{prefix}.iqtree"

    # Parse USER TREES table
    tests: list[dict[str, Any]] = []
    parse_warnings: list[str] = []
    if iqtree_report.exists() and not dry_run:
        tests, parse_warnings = _parse_user_trees_table(iqtree_report)
        if parse_warnings:
            warnings.extend(parse_warnings)
        if not tests and returncode == 0:
            warnings.append("USER TREES parsing produced no rows; IQ-TREE completed successfully")

    # Discover optimized trees
    optimized_trees = _discover_optimized_trees(output_dir, prefix)

    # Key results
    n_candidate = len(tests) if tests else None
    best_tree_id = None
    n_rejected_au = None
    if tests:
        candidates_with_delta = [
            t for t in tests if t.get("delta_likelihood") is not None
        ]
        if candidates_with_delta:
            best = min(candidates_with_delta, key=lambda t: t["delta_likelihood"])
            best_tree_id = best.get("tree_id")
        au_vals = [
            t.get("p_au") for t in tests
            if t.get("p_au") is not None
        ]
        if au_vals:
            n_rejected_au = sum(1 for v in au_vals if v < 0.05)

    # Model source tracking — reflect which source actually supplied the
    # -m/-p flag in the final IQ-TREE command, accounting for suppress-if-present.
    _model_expr = params.get("model_expr")
    _partitions = params.get("partitions")
    _tool_args = params.get("tool_args")
    tool_tokens_for_source = set(shlex.split(_tool_args)) if _tool_args else set()
    if _model_expr and _is_flag_overridden("-m", tool_tokens_for_source):
        model_source = "tool-args"
    elif _partitions and _is_flag_overridden("-p", tool_tokens_for_source):
        model_source = "tool-args"
    elif _model_expr:
        model_source = "model-expr"
    elif _partitions:
        model_source = "partitions"
    else:
        model_source = "tool-args"

    key_results: dict[str, Any] = {
        "n_candidate_trees": n_candidate,
        "best_tree_id": best_tree_id,
        "n_rejected_au_0_05": n_rejected_au,
        "replicates": params["replicates"],
        "model_source": model_source,
    }

    data: dict[str, Any] = {
        "cmd": effective_cmd,
        "tool_stderr": tool_stderr,
        "log_iqtree": str(iqtree_report.name) if iqtree_report.exists() else None,
        "tool_log": str(iqtree_log.name) if iqtree_log.exists() else None,
        "optimized_trees": optimized_trees,
        "merged_candidate_trees": (
            str(candidate_trees_effective.name)
            if candidate_trees_mode == "individual-files"
            else None
        ),
        "tests": tests,
        "warnings": warnings,
    }
    # Remove None values so validate_result_json does not fail on
    # isinstance checks when file paths are absent (e.g. dry-run).
    data = {k: v for k, v in data.items() if v is not None}

    return {
        "status": "success" if returncode == 0 else "error",
        "command": command,
        "wall_time": wall_time,
        "tool_versions": tool_versions,
        "params": params,
        "key_results": key_results,
        "error": None if returncode == 0 else f"IQ-TREE exited with code {returncode}",
        "error_category": error_category,
        "data": data,
    }


# ===================================================================
# Main entry point
# ===================================================================

def run_topology(
    *,
    matrix: Path,
    candidate_trees: list[Path],
    input_format: str = "auto",
    model_expr: str | None = None,
    partitions: str | None = None,
    guide_tree: str | None = None,
    replicates: int = 10000,
    prefix: str | None = None,
    output_dir: Path | None = None,
    threads: int = 4,
    iqtree_path: str | None = None,
    tool_args: str | None = None,
    overwrite: bool = False,
    resume: bool = False,
    dry_run: bool = False,
    quiet: bool = False,
) -> dict[str, Any]:
    """Run IQ-TREE tree topology tests.

    Returns a result.json payload dict.
    """
    run_start = _time.time()

    # --- Library-layer validation ---
    errors = _validate_inputs(
        matrix=matrix,
        candidate_trees=candidate_trees,
        replicates=replicates,
        threads=threads,
        overwrite=overwrite,
        resume=resume,
        model_expr=model_expr,
        partitions=partitions,
        tool_args=tool_args,
        guide_tree=guide_tree,
    )
    if errors:
        return {
            "status": "error",
            "command": "",
            "wall_time": 0.0,
            "tool_versions": {},
            "params": {},
            "key_results": {},
            "error": "; ".join(errors),
            "error_category": "input",
            "data": {"cmd": [], "tool_stderr": "", "tests": [], "warnings": errors},
        }

    # --- Resolve paths ---
    if output_dir is None:
        output_dir = Path("runs/posttree/topology")
    output_dir = output_dir.resolve()
    if not dry_run:
        output_dir.mkdir(parents=True, exist_ok=True)

    matrix = matrix.resolve()

    candidate_trees_raw = [ct.resolve() for ct in candidate_trees]
    if len(candidate_trees) == 1:
        candidate_trees_mode = "tree-list"
        candidate_trees_effective = candidate_trees_raw[0]
    else:
        candidate_trees_mode = "individual-files"
        if not dry_run:
            candidate_trees_effective = _merge_candidate_trees(
                candidate_trees_raw, output_dir,
            )
        else:
            candidate_trees_effective = output_dir / "candidate.trees"

    resolved_prefix = prefix if prefix else matrix.stem

    # --- Resolve IQ-TREE ---
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
            "data": {"cmd": [], "tool_stderr": "", "tests": [], "warnings": [str(e)]},
        }
    tool_versions = (
        _detect_iqtree_version(iqtree_exe)
        if not dry_run
        else {"iqtree3": "dry-run"}
    )

    # --- Params dict (complete, all parameters) ---
    params: dict[str, Any] = {
        "matrix": str(matrix),
        "candidate_trees": [str(ct) for ct in candidate_trees_raw],
        "candidate_trees_mode": candidate_trees_mode,
        "candidate_trees_effective": str(candidate_trees_effective),
        "input_format": input_format,
        "model_expr": model_expr,
        "partitions": partitions,
        "guide_tree": guide_tree,
        "replicates": replicates,
        "prefix": resolved_prefix,
        "output_dir": str(output_dir),
        "threads": threads,
        "iqtree_path": iqtree_path,
        "tool_args": tool_args,
        "overwrite": overwrite,
        "resume": resume,
        "dry_run": dry_run,
        "quiet": quiet,
    }

    # --- Build IQ-TREE command ---
    effective_cmd = _build_topology_cmd(
        executable=iqtree_exe,
        matrix=matrix,
        candidate_trees=candidate_trees_effective,
        prefix=resolved_prefix,
        model_expr=model_expr,
        partitions=partitions,
        guide_tree=guide_tree,
        replicates=replicates,
        threads=threads,
        tool_args=tool_args,
    )

    if dry_run:
        return _assemble_topology_result(
            run_start=run_start,
            iqtree_exe=iqtree_exe,
            tool_versions=tool_versions,
            params=params,
            matrix=matrix,
            candidate_trees_raw=candidate_trees_raw,
            candidate_trees_mode=candidate_trees_mode,
            candidate_trees_effective=candidate_trees_effective,
            effective_cmd=effective_cmd,
            tool_stderr="",
            output_dir=output_dir,
            prefix=resolved_prefix,
            returncode=0,
            dry_run=True,
            warnings=[],
        )

    # --- Execute IQ-TREE ---
    try:
        proc = subprocess.run(
            effective_cmd,
            capture_output=True,
            text=True,
            cwd=str(output_dir),
        )
    except FileNotFoundError:
        return _assemble_topology_result(
            run_start=run_start,
            iqtree_exe=iqtree_exe,
            tool_versions=tool_versions,
            params=params,
            matrix=matrix,
            candidate_trees_raw=candidate_trees_raw,
            candidate_trees_mode=candidate_trees_mode,
            candidate_trees_effective=candidate_trees_effective,
            effective_cmd=effective_cmd,
            tool_stderr=f"Executable not found: {iqtree_exe}",
            output_dir=output_dir,
            prefix=resolved_prefix,
            returncode=3,
            dry_run=False,
            warnings=[],
            error_category="env",
        )

    merged_output = ""
    if proc.stdout:
        merged_output += proc.stdout
    if proc.stderr:
        if merged_output and not merged_output.endswith("\n"):
            merged_output += "\n"
        merged_output += proc.stderr

    return _assemble_topology_result(
        run_start=run_start,
        iqtree_exe=iqtree_exe,
        tool_versions=tool_versions,
        params=params,
        matrix=matrix,
        candidate_trees_raw=candidate_trees_raw,
        candidate_trees_mode=candidate_trees_mode,
        candidate_trees_effective=candidate_trees_effective,
        effective_cmd=effective_cmd,
        tool_stderr=merged_output,
        output_dir=output_dir,
        prefix=resolved_prefix,
        returncode=proc.returncode,
        dry_run=False,
        warnings=[],
        error_category="tool" if proc.returncode != 0 else None,
    )
