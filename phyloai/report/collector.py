"""Step discovery, ordering, and run-mode detection for phyloai report."""

from __future__ import annotations

import json
import shlex
from pathlib import Path
from typing import Any

STEP_ORDER: list[str] = [
    "pretree.convert",
    "pretree.stats",
    "pretree.align",
    "pretree.trim",
    "pretree.metrics",
    "pretree.filter.taper",
    "pretree.filter.treeshrink",
    "pretree.filter.symtest",
    "pretree.filter.metrics",
    "pretree.filter.cluster",
    "pretree.concat",
    "pretree.concat.jackknife",
    "tree.ml.fasttree",
    "tree.ml.iqtree",
    "tree.msc",
    "tree.bi.pb",
    "tree.bi.bpcomp",
    "tree.bi.tracecomp",
    "tree.bi.readpb",
    "tree.cf",
    "posttree.topology",
    "posttree.dating.hessian",
    "posttree.dating.mcmc",
    "posttree.signal",
    "posttree.syserror.brlen",
    "posttree.syserror.cca",
    "posttree.syserror.sites",
    "posttree.simulate",
]

_EXCLUDE_DIRS = {"report", "logs"}


def parse_step_id(command: str) -> str:
    """Extract step_id from a full CLI command string.

    Drops all flag tokens (any token starting with '-'), then finds the
    first known root command (pretree/tree/posttree/run) in the remaining
    tokens and reads subsequent tokens as the subcommand path.  Flag VALUES
    (e.g. ``./runs`` from ``--run-dir ./runs``) are harmless noise that
    don't match any known root and are naturally skipped by the first-match
    lookup.

    Boolean flags before the root are handled correctly because only the
    flag token ``--quiet`` is dropped; the next token ``pretree`` is not
    consumed as a flag value.
    """
    _ROOT_SUBCOMMANDS: dict[str, set[str] | None] = {
        "pretree": {"convert", "stats", "align", "trim", "metrics", "filter", "concat"},
        "tree":    {"ml", "bi", "msc", "cf"},
        "posttree":{"topology", "dating", "signal", "syserror", "simulate"},
        "run":     None,
        "doctor":  None,
    }
    _THIRD_LEVEL: dict[str, set[str]] = {
        "filter": {"taper", "treeshrink", "symtest", "metrics", "cluster"},
        "concat": {"jackknife"},
        "ml":     {"fasttree", "iqtree"},
        "dating": {"hessian", "mcmc"},
        "syserror": {"brlen", "cca", "sites"},
        "bi": {"pb", "bpcomp", "tracecomp", "readpb"},
    }

    try:
        tokens = shlex.split(command)
    except ValueError:
        tokens = command.split()

    if not tokens or tokens[0] != "phyloai":
        return "unknown"

    parts: list[str] = [t for t in tokens[1:] if not t.startswith("-")]
    if not parts:
        return "unknown"

    root_index = -1
    root = ""
    for idx, p in enumerate(parts):
        if p in _ROOT_SUBCOMMANDS:
            root_index = idx
            root = p
            break

    if root_index < 0:
        return "unknown"

    subcmds = _ROOT_SUBCOMMANDS[root]
    if subcmds is None:
        return root

    after = parts[root_index + 1:]
    if not after:
        return root

    l2 = after[0]
    if l2 not in subcmds:
        return f"{root}.{l2}" if l2 else root
    if len(after) < 2:
        return f"{root}.{l2}"

    l3 = after[1]
    third = _THIRD_LEVEL.get(l2, set())
    if l3 in third:
        return f"{root}.{l2}.{l3}"
    return f"{root}.{l2}"


def _load_result_json(path: Path) -> dict[str, Any]:
    """Load and return parsed result.json, or error dict on failure."""
    try:
        with open(path) as fh:
            return json.load(fh)
    except (json.JSONDecodeError, OSError) as exc:
        return {"status": "error", "error": str(exc), "command": f"<unreadable: {path}>"}


def _find_result_json_files(base: Path, max_depth: int = 3) -> list[Path]:
    """Find result.json files up to max_depth levels below base.

    Skips directories named 'report' and 'logs' to avoid self-referential
    scanning and tool log directories.  Dot-prefixed dirs are also skipped.

    Depth=3 covers pipeline layouts with nested tool outputs
    (e.g. ``runs/run/faa/tree/ml/iqtree/result.json`` at depth 3).
    """
    found: list[Path] = []
    todo: list[tuple[Path, int]] = [(base, 0)]
    while todo:
        current, depth = todo.pop()
        if depth >= max_depth:
            continue
        try:
            entries = sorted(current.iterdir())
        except PermissionError:
            continue
        for entry in entries:
            if not entry.is_dir():
                continue
            if entry.name in _EXCLUDE_DIRS or entry.name.startswith("."):
                continue
            candidate = entry / "result.json"
            if candidate.exists():
                found.append(candidate)
            todo.append((entry, depth + 1))
    return found


def _detect_run_mode(run_dir: Path) -> tuple[str, list[Path], dict[str, Any] | None]:
    """Detect run_mode by filesystem structure only (Spec Section 5 priority rules).

    Priority:
    1. Top-level result.json + subdirs with result.json → pipeline
    2. Top-level result.json only → module (single command)
    3. No top-level, but subdirs → module (multi-step)
    4. Nothing found → error

    Pipeline step paths come purely from filesystem scan — NOT from
    top-level result.json:data.steps[], which is an implementation detail
    of phyloai run that may change.  The top-level result.json is only
    read for optional metadata (mode, speed) enrichment.
    """
    top_result = run_dir / "result.json"
    sub_results = _find_result_json_files(run_dir)

    has_top = top_result.exists()

    if has_top and sub_results:
        top_data = _load_result_json(top_result)
        pipeline_summary = {
            "mode": top_data.get("data", {}).get("mode", "unknown"),
            "speed": top_data.get("data", {}).get("speed", "normal"),
        }
        # Enrich with top-level key_results (input_genes, final_tree, etc.)
        top_kr = top_data.get("key_results", {})
        for k in ("n_input_genes", "n_genes_after_filter", "final_tree",
                  "matrix_length", "matrix_taxa", "input_genes",
                  "genes_after_filter"):
            if k in top_kr:
                pipeline_summary[k] = top_kr[k]
        return ("pipeline", sub_results, pipeline_summary)

    if has_top and not sub_results:
        return ("module", [top_result], None)

    if not has_top and sub_results:
        return ("module", sub_results, None)

    raise ValueError(
        f"No result.json found in {run_dir} or its subdirectories. "
        f"Ensure the directory contains PhyloAI run output."
    )


def _order_steps(steps: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Sort steps by STEP_ORDER, appending unknowns at end."""
    order_map = {sid: i for i, sid in enumerate(STEP_ORDER)}

    def sort_key(step: dict[str, Any]) -> tuple[int, str]:
        sid = step.get("step_id", "unknown")
        return (order_map.get(sid, 9999), sid)

    return sorted(steps, key=sort_key)


def discover_steps(run_dir: str | Path) -> dict[str, Any]:
    """Scan run_dir and return structured step data for report assembly.

    Returns:
        dict with keys: run_mode, steps, pipeline_summary
    """
    run_dir = Path(run_dir).resolve()

    run_mode, result_paths, pipeline_summary = _detect_run_mode(run_dir)

    steps: list[dict[str, Any]] = []
    for rp in result_paths:
        result = _load_result_json(rp)
        step_id = parse_step_id(result.get("command", ""))
        steps.append({
            "step_id": step_id,
            "path": str(rp.absolute()),
            "command": result.get("command", ""),
            "status": result.get("status", "error"),
            "wall_time": result.get("wall_time", 0.0),
            "tool_versions": result.get("tool_versions", {}),
            "params": result.get("params", {}),
            "key_results": result.get("key_results", {}),
            "error": result.get("error"),
            "data": result.get("data", {}),
        })

    steps = _order_steps(steps)

    return {
        "run_mode": run_mode,
        "steps": steps,
        "pipeline_summary": pipeline_summary,
    }
