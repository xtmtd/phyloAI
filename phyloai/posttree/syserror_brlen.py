"""Branch-length extraction for systematic-error diagnosis (pure Python).

Parses Newick trees with Bio.Phylo and computes branch-length statistics
(total, terminal, internal, patristic, and endpoint distances) with
topology-aware node resolution via node maps or labeled trees.

Tree representation is inferred structurally from the Newick root child
count: 2 children = rooted, 3+ = unrooted. This is a structural heuristic,
not a proof of biological rooting.
"""

from __future__ import annotations

import csv
import itertools
import json
import math
import shlex
import shutil
import statistics
import time
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

from Bio import Phylo
from Bio.Phylo.BaseTree import Clade, Tree

_BATCH_MODES = frozenset({"total", "terminal", "internal", "patristic"})
_ENDPOINT_MODES = frozenset({"tip-to-tip", "node-to-node", "node-to-tip"})
_ALL_MODES = _BATCH_MODES | _ENDPOINT_MODES | frozenset({"all"})

_TABLE_COLUMNS = {
    "total": ["tree_file", "total_branch_length"],
    "terminal": ["tree_file", "taxon", "branch_length"],
    "internal": ["tree_file", "representation", "edge_taxa", "branch_length"],
    "patristic": ["tree_file", "tip1", "tip2", "distance"],
    "tip-to-tip": ["tree_file", "tip1", "tip2", "distance"],
    "node-to-node": ["tree_file", "node1", "node2", "node1_type", "node2_type", "distance"],
    "node-to-tip": ["tree_file", "node", "node_type", "tip", "distance"],
}

_OUTPUT_FILE_DESCRIPTIONS = {
    "total": "Total branch length per tree",
    "terminal": "Terminal branch lengths per taxon per tree",
    "internal": "Internal branch lengths per tree (edge_taxa with rooted/unrooted representation)",
    "patristic": "Pairwise tip-to-tip distances per tree",
    "tip-to-tip": "Distance between the specified tip pair per tree",
    "node-to-node": "Distance between the specified node pair per tree",
    "node-to-tip": "Distance between a node and tips per tree",
}


@dataclass
class TreeRecord:
    tree_id: str
    tree: Tree


@dataclass
class ResolvedEndpoint:
    clade: Clade
    node_type: str
    present_taxa: frozenset[str]
    warning: str | None = None


def _is_rooted_representation(tree: Tree) -> bool:
    return len(tree.root.clades) == 2


def _branch_length(clade: Clade) -> float:
    return float(clade.branch_length or 0.0)


def _canonical_split(side: frozenset[str], all_tips: frozenset[str]) -> tuple[str, ...]:
    other = all_tips - side
    candidates = (tuple(sorted(side)), tuple(sorted(other)))
    return min(candidates, key=lambda values: (len(values), values))


def _read_tree_file(path: Path) -> tuple[list[TreeRecord], list[str], int]:
    try:
        trees = list(Phylo.parse(path, "newick"))
    except Exception as exc:  # NewickError and friends
        raise ValueError(f"tree file {path.name}: failed to parse as Newick, skipped") from exc
    records: list[TreeRecord] = []
    warnings: list[str] = []
    skipped = 0
    if not trees:
        return records, [f"tree file {path.name}: no Newick trees found, skipped"], 1
    for index, tree in enumerate(trees):
        if len(tree.get_terminals()) < 2:
            warnings.append(f"tree {path.name}:{index}: fewer than two tips, skipped")
            skipped += 1
            continue
        tree_id = path.name if len(trees) == 1 else f"{path.name}:{index}"
        records.append(TreeRecord(tree_id, tree))
    return records, warnings, skipped


def _read_tree_dir(path: Path) -> tuple[list[TreeRecord], list[str], int]:
    records: list[TreeRecord] = []
    warnings: list[str] = []
    skipped = 0
    try:
        entries = sorted(p for p in path.iterdir() if p.is_file())
    except OSError as exc:
        raise ValueError(f"cannot read tree directory {path}: {exc}") from exc
    for entry in entries:
        try:
            if entry.stat().st_size == 0:
                continue
        except OSError:
            continue
        try:
            file_records, file_warnings, file_skipped = _read_tree_file(entry)
        except ValueError as exc:
            warnings.append(str(exc))
            skipped += 1
            continue
        records.extend(file_records)
        warnings.extend(file_warnings)
        skipped += file_skipped
    return records, warnings, skipped


def _missing_length_warning(record: TreeRecord) -> str:
    return f"tree {record.tree_id}: all branch lengths are missing; treating them as 0.0"


def _terminal_names(clade: Clade) -> list[str]:
    return sorted(terminal.name for terminal in clade.get_terminals() if terminal.name)


def _total_rows(record: TreeRecord) -> list[dict[str, Any]]:
    total = sum(_branch_length(clade) for clade in record.tree.find_clades())
    return [{"tree_file": record.tree_id, "total_branch_length": total}]


def _terminal_rows(record: TreeRecord) -> list[dict[str, Any]]:
    rows = []
    for terminal in record.tree.get_terminals():
        if terminal.name is None:
            continue
        rows.append({
            "tree_file": record.tree_id,
            "taxon": terminal.name,
            "branch_length": _branch_length(terminal),
        })
    rows.sort(key=lambda row: (row["taxon"], row["branch_length"]))
    return rows


def _internal_rows(record: TreeRecord) -> list[dict[str, Any]]:
    rooted = _is_rooted_representation(record.tree)
    all_tips = frozenset(terminal.name for terminal in record.tree.get_terminals() if terminal.name)
    rows = []
    for clade in record.tree.get_nonterminals():
        if clade is record.tree.root:
            continue
        if rooted:
            edge_taxa = _terminal_names(clade)
        else:
            side = frozenset(terminal.name for terminal in clade.get_terminals() if terminal.name)
            edge_taxa = list(_canonical_split(side, all_tips))
        rows.append({
            "tree_file": record.tree_id,
            "representation": "rooted" if rooted else "unrooted",
            "edge_taxa": ",".join(edge_taxa),
            "branch_length": _branch_length(clade),
        })
    return rows


def _patristic_rows(record: TreeRecord) -> Iterable[dict[str, Any]]:
    terminals = [terminal.name for terminal in record.tree.get_terminals() if terminal.name]
    terminals.sort()
    for tip1, tip2 in itertools.combinations(terminals, 2):
        yield {
            "tree_file": record.tree_id,
            "tip1": tip1,
            "tip2": tip2,
            "distance": record.tree.distance(tip1, tip2),
        }


def _parse_map(path: Path) -> dict[str, frozenset[str]]:
    try:
        text = path.read_text()
    except (OSError, UnicodeDecodeError) as exc:
        raise ValueError(f"map file {path}: cannot read ({exc})") from exc
    mapping: dict[str, frozenset[str]] = {}
    for raw_line in text.splitlines():
        if ":" not in raw_line:
            continue
        label, raw_taxa = raw_line.split(":", 1)
        taxa = frozenset(value.strip() for value in raw_taxa.split(",") if value.strip())
        if label.strip() and taxa:
            mapping[label.strip()] = taxa
    return mapping


def _tip_clade(tree: Tree, name: str) -> Clade:
    for terminal in tree.get_terminals():
        if terminal.name == name:
            return terminal
    raise ValueError(f"tip '{name}' not found")


def _resolve_labeled_node(tree: Tree, name: str) -> ResolvedEndpoint:
    for clade in tree.get_nonterminals():
        if clade.name == name:
            return ResolvedEndpoint(clade, "internal", frozenset(_terminal_names(clade)))
    for terminal in tree.get_terminals():
        if terminal.name == name:
            return ResolvedEndpoint(terminal, "tip", frozenset({name}))
    raise ValueError(f"node '{name}' not found in labeled tree")


def _resolve_map_endpoint(
    tree: Tree, name: str, node_map: dict[str, frozenset[str]]
) -> ResolvedEndpoint:
    if name not in node_map:
        raise ValueError(f"map does not define node '{name}'")
    all_tips = frozenset(terminal.name for terminal in tree.get_terminals() if terminal.name)
    present = node_map[name] & all_tips
    if not present:
        raise ValueError(f"node '{name}': no taxa in map group present in this tree")
    if len(present) == 1:
        return ResolvedEndpoint(_tip_clade(tree, next(iter(present))), "tip", present)
    if _is_rooted_representation(tree):
        clade = tree.common_ancestor(sorted(present))
        if frozenset(_terminal_names(clade)) != present:
            raise ValueError(f"node '{name}' is not monophyletic in this tree")
        return ResolvedEndpoint(clade, "internal", present)
    parents = {child: parent for parent in tree.find_clades() for child in parent.clades}
    for clade in tree.get_nonterminals():
        if clade is tree.root:
            continue
        descendants = frozenset(terminal.name for terminal in clade.get_terminals() if terminal.name)
        if descendants == present:
            return ResolvedEndpoint(clade, "internal", present)
        if all_tips - descendants == present:
            parent = parents[clade]
            return ResolvedEndpoint(
                parent,
                "internal",
                present,
                warning=f"node '{name}': matched complement side of a split in this tree",
            )
    raise ValueError(f"node '{name}' is not a valid split in this tree")


def _resolve_endpoint(
    tree: Tree, name: str, node_map: dict[str, frozenset[str]] | None
) -> ResolvedEndpoint:
    if node_map is not None:
        return _resolve_map_endpoint(tree, name, node_map)
    return _resolve_labeled_node(tree, name)


def _tip_to_tip_rows(record: TreeRecord, tip1: str, tip2: str) -> list[dict[str, Any]]:
    tree = record.tree
    _tip_clade(tree, tip1)
    _tip_clade(tree, tip2)
    return [{
        "tree_file": record.tree_id,
        "tip1": tip1,
        "tip2": tip2,
        "distance": tree.distance(tip1, tip2),
    }]


def _node_to_node_rows(
    record: TreeRecord,
    node1: str,
    node2: str,
    node_map: dict[str, frozenset[str]] | None,
    warnings: list[str] | None = None,
) -> list[dict[str, Any]]:
    tree = record.tree
    left = _resolve_endpoint(tree, node1, node_map)
    right = _resolve_endpoint(tree, node2, node_map)
    if warnings is not None:
        for endpoint in (left, right):
            if endpoint.warning:
                warnings.append(f"tree {record.tree_id}: {endpoint.warning}")
    return [{
        "tree_file": record.tree_id,
        "node1": node1,
        "node2": node2,
        "node1_type": left.node_type,
        "node2_type": right.node_type,
        "distance": tree.distance(left.clade, right.clade),
    }]


def _node_to_tip_rows(
    record: TreeRecord,
    node1: str,
    tip1: str | None,
    node_map: dict[str, frozenset[str]] | None,
    warnings: list[str] | None = None,
) -> list[dict[str, Any]]:
    tree = record.tree
    endpoint = _resolve_endpoint(tree, node1, node_map)
    if warnings is not None and endpoint.warning:
        warnings.append(f"tree {record.tree_id}: {endpoint.warning}")
    if tip1 is not None:
        tip_clade = _tip_clade(tree, tip1)
        return [{
            "tree_file": record.tree_id,
            "node": node1,
            "node_type": endpoint.node_type,
            "tip": tip1,
            "distance": tree.distance(endpoint.clade, tip_clade),
        }]
    if node_map is not None:
        rows = []
        for tip in sorted(endpoint.present_taxa):
            tip_clade = _tip_clade(tree, tip)
            rows.append({
                "tree_file": record.tree_id,
                "node": node1,
                "node_type": endpoint.node_type,
                "tip": tip,
                "distance": tree.distance(endpoint.clade, tip_clade),
            })
        return rows
    if not _is_rooted_representation(tree):
        raise ValueError("--tip1 or --map required for unrooted trees (descendant ambiguity)")
    rows = []
    for tip in sorted(_terminal_names(endpoint.clade)):
        tip_clade = _tip_clade(tree, tip)
        rows.append({
            "tree_file": record.tree_id,
            "node": node1,
            "node_type": endpoint.node_type,
            "tip": tip,
            "distance": tree.distance(endpoint.clade, tip_clade),
        })
    return rows


def _parse_modes(mode: str) -> tuple[list[str], str | None]:
    tokens = [token.strip() for token in mode.split(",") if token.strip()]
    if not tokens:
        raise ValueError("--mode is required")
    unknown = [token for token in tokens if token not in _ALL_MODES]
    if unknown:
        raise ValueError(f"unknown mode(s): {', '.join(unknown)}")
    if "all" in tokens:
        if len(tokens) != 1:
            raise ValueError("--mode all is mutually exclusive with other modes")
        return ["total", "terminal", "internal", "patristic"], None
    endpoints = [token for token in tokens if token in _ENDPOINT_MODES]
    if endpoints:
        if len(tokens) != 1:
            raise ValueError("endpoint modes cannot be combined with batch modes or other modes")
        return tokens, endpoints[0]
    return tokens, None


def _validate_params(
    tree: Path | None,
    tree_dir: Path | None,
    mode: str,
    map_file: Path | None,
    table_format: str,
    threads: int,
    max_rows: int,
) -> None:
    if (tree is None) == (tree_dir is None):
        raise ValueError("exactly one of --tree or --tree-dir is required")
    if table_format not in ("csv", "tsv"):
        raise ValueError(f"invalid table format: {table_format}")
    if threads < 1:
        raise ValueError("--threads must be >= 1")
    if max_rows < 0:
        raise ValueError("--max-rows must be >= 0")
    if tree is not None and not tree.exists():
        raise ValueError(f"tree file not found: {tree}")
    if tree_dir is not None and not tree_dir.is_dir():
        raise ValueError(f"tree directory not found: {tree_dir}")
    if map_file is not None and not map_file.is_file():
        raise ValueError(f"map file not readable: {map_file}")
    _parse_modes(mode)


def _summarize(values: list[float]) -> dict[str, Any]:
    if not values:
        return {"n_values": 0, "mean": None, "sd": None, "min": None, "max": None}
    return {
        "n_values": len(values),
        "mean": statistics.fmean(values),
        "sd": statistics.pstdev(values),
        "min": min(values),
        "max": max(values),
    }


def _record_rows(
    record: TreeRecord,
    modes: list[str],
    endpoint: str | None,
    node_map: dict[str, frozenset[str]] | None,
    node1: str | None,
    node2: str | None,
    tip1: str | None,
    tip2: str | None,
    warnings: list[str] | None = None,
) -> dict[str, list[dict[str, Any]]]:
    rows: dict[str, list[dict[str, Any]]] = {}
    if "total" in modes:
        rows["total"] = _total_rows(record)
    if "terminal" in modes:
        rows["terminal"] = _terminal_rows(record)
    if "internal" in modes:
        rows["internal"] = _internal_rows(record)
    if endpoint == "tip-to-tip":
        rows["tip-to-tip"] = _tip_to_tip_rows(record, tip1, tip2)
    elif endpoint == "node-to-node":
        rows["node-to-node"] = _node_to_node_rows(record, node1, node2, node_map, warnings)
    elif endpoint == "node-to-tip":
        rows["node-to-tip"] = _node_to_tip_rows(record, node1, tip1, node_map, warnings)
    return rows


def _process_file(
    path_str: str,
    modes: list[str],
    endpoint: str | None,
    node_map: dict[str, frozenset[str]] | None,
    node1: str | None,
    node2: str | None,
    tip1: str | None,
    tip2: str | None,
) -> tuple[dict[str, list[dict[str, Any]]], list[str], int]:
    rows_by_mode: dict[str, list[dict[str, Any]]] = {}
    warnings: list[str] = []
    skipped = 0
    try:
        records, _, _ = _read_tree_file(Path(path_str))
    except ValueError:
        return rows_by_mode, warnings, skipped
    for record in records:
        try:
            record_rows = _record_rows(
                record, modes, endpoint, node_map, node1, node2, tip1, tip2, warnings
            )
        except ValueError as exc:
            warnings.append(f"tree {record.tree_id}: {exc}")
            skipped += 1
            continue
        for mode, rows in record_rows.items():
            rows_by_mode.setdefault(mode, []).extend(rows)
    return rows_by_mode, warnings, skipped


def _prepare_output_dir(output_dir: Path, overwrite: bool) -> None:
    output_dir = output_dir.resolve()
    if output_dir.exists() and any(output_dir.iterdir()):
        if not overwrite:
            raise ValueError(
                f"Output directory exists and is not empty: {output_dir}\n"
                "Use --overwrite to replace."
            )
        shutil.rmtree(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    return output_dir


def _build_command(
    tree: Path | None,
    tree_dir: Path | None,
    mode: str,
    map_file: Path | None,
    node1: str | None,
    node2: str | None,
    tip1: str | None,
    tip2: str | None,
    output_dir: Path,
    table_format: str,
    threads: int,
    max_rows: int,
    overwrite: bool,
    dry_run: bool,
    quiet: bool,
) -> str:
    parts = ["phyloai", "posttree", "syserror", "brlen"]
    if tree:
        parts += ["--tree", str(tree)]
    if tree_dir:
        parts += ["--tree-dir", str(tree_dir)]
    parts += ["--mode", mode]
    if map_file:
        parts += ["--map", str(map_file)]
    if node1:
        parts += ["--node1", node1]
    if node2:
        parts += ["--node2", node2]
    if tip1:
        parts += ["--tip1", tip1]
    if tip2:
        parts += ["--tip2", tip2]
    parts += ["-o", str(output_dir)]
    if table_format != "csv":
        parts += ["--table-format", table_format]
    if threads != 4:
        parts += ["-t", str(threads)]
    if max_rows != 5_000_000:
        parts += ["--max-rows", str(max_rows)]
    if overwrite:
        parts.append("--overwrite")
    if dry_run:
        parts.append("--dry-run")
    if quiet:
        parts.append("-q")
    return shlex.join(parts)


def _write_table(path: Path, rows: list[dict[str, Any]], columns: list[str], delimiter: str) -> None:
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=columns, delimiter=delimiter)
        writer.writeheader()
        writer.writerows(rows)


def run_brlen(
    tree: Path | None,
    tree_dir: Path | None,
    mode: str,
    map_file: Path | None = None,
    node1: str | None = None,
    node2: str | None = None,
    tip1: str | None = None,
    tip2: str | None = None,
    table_format: str = "csv",
    threads: int = 4,
    max_rows: int = 5_000_000,
    output_dir: Path = Path("runs/posttree/syserror/brlen"),
    overwrite: bool = False,
    dry_run: bool = False,
    quiet: bool = False,
) -> dict[str, Any]:
    """Run branch-length extraction and return a spec-compliant payload."""
    start = time.time()
    _validate_params(tree, tree_dir, mode, map_file, table_format, threads, max_rows)
    modes, endpoint = _parse_modes(mode)
    node_map = _parse_map(map_file) if map_file is not None else None

    if endpoint == "tip-to-tip":
        if tip1 is None or tip2 is None:
            raise ValueError("--mode tip-to-tip requires both --tip1 and --tip2")
    elif endpoint == "node-to-node":
        if node1 is None or node2 is None:
            raise ValueError("--mode node-to-node requires both --node1 and --node2")
    elif endpoint == "node-to-tip":
        if node1 is None:
            raise ValueError("--mode node-to-tip requires --node1")

    warnings: list[str] = []
    if tree is not None:
        records, read_warnings, skipped = _read_tree_file(tree)
        warnings.extend(read_warnings)
    else:
        records, read_warnings, skipped = _read_tree_dir(tree_dir)
        warnings.extend(read_warnings)
    if not records:
        raise ValueError("no valid trees found; nothing to compute")

    base_counts = Counter(record.tree_id.split(":")[0] for record in records)
    n_multi_tree_files = sum(1 for count in base_counts.values() if count > 1)

    if "patristic" in modes:
        estimate = sum(
            len(record.tree.get_terminals()) * (len(record.tree.get_terminals()) - 1) // 2
            for record in records
        )
        warnings.append(f"estimated patristic rows: {estimate}")
        if max_rows and estimate > max_rows:
            raise ValueError(
                f"patristic output would exceed --max-rows ({estimate} > {max_rows}); "
                "use --max-rows 0 to disable the limit"
            )

    for record in records:
        if all(clade.branch_length is None for clade in record.tree.find_clades()):
            warnings.append(_missing_length_warning(record))

    resolved_output_dir = output_dir.resolve()
    tables_dir = None
    if not dry_run:
        resolved_output_dir = _prepare_output_dir(resolved_output_dir, overwrite)
        tables_dir = resolved_output_dir / "tables"
        tables_dir.mkdir(parents=True, exist_ok=True)

    delimiter = "\t" if table_format == "tsv" else ","
    suffix = table_format

    mode_rows: dict[str, list[dict[str, Any]]] = {m: [] for m in modes if m != "patristic"}
    endpoint_skipped = 0
    files = sorted(p for p in tree_dir.iterdir() if p.is_file()) if tree is None else None
    non_patristic_modes = [m for m in modes if m != "patristic"]
    progress = None
    task = None
    if not quiet and not dry_run and (non_patristic_modes or "patristic" in modes):
        from rich.progress import Progress

        progress = Progress(transient=True)
        total_units = 0
        if non_patristic_modes:
            total_units += len(records) if tree is not None else len(files)
        if "patristic" in modes:
            total_units += len(records)
        task = progress.add_task("Processing trees", total=total_units)
        progress.start()
    if tree is not None:
        for record in records:
            try:
                record_rows = _record_rows(
                    record, modes, endpoint, node_map, node1, node2, tip1, tip2, warnings
                )
            except ValueError as exc:
                raise ValueError(str(exc)) from exc
            if not dry_run:
                for m, rows in record_rows.items():
                    mode_rows.setdefault(m, []).extend(rows)
            if progress is not None and non_patristic_modes:
                progress.advance(task)
    else:
        use_pool = not dry_run and non_patristic_modes and threads > 1
        if use_pool:
            with ProcessPoolExecutor(max_workers=threads) as pool:
                futures = {
                    path: pool.submit(
                        _process_file,
                        str(path),
                        modes,
                        endpoint,
                        node_map,
                        node1,
                        node2,
                        tip1,
                        tip2,
                    )
                    for path in files
                }
                for path in files:
                    file_rows, file_warnings, file_skipped = futures[path].result()
                    warnings.extend(file_warnings)
                    endpoint_skipped += file_skipped
                    for m, rows in file_rows.items():
                        mode_rows.setdefault(m, []).extend(rows)
                    if progress is not None:
                        progress.advance(task)
        else:
            for path in files:
                file_rows, file_warnings, file_skipped = _process_file(
                    str(path), modes, endpoint, node_map, node1, node2, tip1, tip2
                )
                warnings.extend(file_warnings)
                endpoint_skipped += file_skipped
                if not dry_run:
                    for m, rows in file_rows.items():
                        mode_rows.setdefault(m, []).extend(rows)
                if progress is not None and non_patristic_modes:
                    progress.advance(task)
    skipped += endpoint_skipped

    output_files: dict[str, dict[str, str]] = {}
    if not dry_run:
        for m in modes:
            if m == "patristic":
                table_path = tables_dir / f"patristic.{suffix}"
                patristic_count = 0
                patristic_sum = 0.0
                patristic_sumsq = 0.0
                patristic_min = None
                patristic_max = None
                with open(table_path, "w", newline="") as fh:
                    writer = csv.DictWriter(
                        fh, fieldnames=_TABLE_COLUMNS["patristic"], delimiter=delimiter
                    )
                    writer.writeheader()
                    for record in records:
                        for row in _patristic_rows(record):
                            writer.writerow(row)
                            value = row["distance"]
                            patristic_count += 1
                            patristic_sum += value
                            patristic_sumsq += value * value
                            if patristic_min is None or value < patristic_min:
                                patristic_min = value
                            if patristic_max is None or value > patristic_max:
                                patristic_max = value
                        if progress is not None:
                            progress.advance(task)
                if patristic_count:
                    mean = patristic_sum / patristic_count
                    sd = math.sqrt(max(patristic_sumsq / patristic_count - mean * mean, 0.0))
                else:
                    mean = None
                    sd = None
                mode_summary = {
                    "n_values": patristic_count,
                    "mean": mean,
                    "sd": sd,
                    "min": patristic_min,
                    "max": patristic_max,
                }
                output_files["patristic_table"] = {
                    "path": str(table_path),
                    "description": _OUTPUT_FILE_DESCRIPTIONS["patristic"],
                }
            else:
                table_path = tables_dir / f"{m.replace('-', '_')}.{suffix}"
                _write_table(table_path, mode_rows.get(m, []), _TABLE_COLUMNS[m], delimiter)
                output_files[f"{m.replace('-', '_')}_table"] = {
                    "path": str(table_path),
                    "description": _OUTPUT_FILE_DESCRIPTIONS[m],
                }

    if progress is not None:
        progress.stop()

    summary: dict[str, Any] = {}
    for m in modes:
        if m == "patristic":
            summary[m] = mode_summary if not dry_run else _summarize([])
        elif m == "total":
            summary[m] = _summarize([row["total_branch_length"] for row in mode_rows.get(m, [])])
        elif m in ("terminal", "internal"):
            summary[m] = _summarize([row["branch_length"] for row in mode_rows.get(m, [])])
        else:
            summary[m] = _summarize([row["distance"] for row in mode_rows.get(m, [])])

    key_results = {
        "n_trees": len(records),
        "n_trees_skipped": skipped,
        "modes": modes,
        "summary": summary,
    }
    data = {
        "summary": {
            "n_trees_processed": len(records),
            "n_trees_skipped": skipped,
            "n_multi_tree_files": n_multi_tree_files,
        },
        "warnings": warnings,
        "output_files": output_files,
    }
    payload = {
        "status": "success",
        "command": _build_command(
            tree, tree_dir, mode, map_file, node1, node2, tip1, tip2,
            output_dir, table_format, threads, max_rows, overwrite, dry_run, quiet,
        ),
        "wall_time": time.time() - start,
        "tool_versions": {},
        "params": {
            "tree": str(tree) if tree else None,
            "tree_dir": str(tree_dir) if tree_dir else None,
            "mode": mode,
            "map": str(map_file) if map_file else None,
            "node1": node1,
            "node2": node2,
            "tip1": tip1,
            "tip2": tip2,
            "table_format": table_format,
            "threads": threads,
            "max_rows": max_rows,
            "output_dir": str(output_dir),
            "overwrite": overwrite,
            "dry_run": dry_run,
            "quiet": quiet,
        },
        "key_results": key_results,
        "error": None,
        "data": data,
    }

    if not dry_run:
        result_path = resolved_output_dir / "result.json"
        with open(result_path, "w") as fh:
            json.dump(payload, fh, indent=2)

    if not quiet and not dry_run:
        click_echo = __import__("click").echo
        click_echo(f"Status:    {payload['status']}")
        click_echo(f"Trees:     {len(records)} processed, {skipped} skipped")
        click_echo(f"Result:    {resolved_output_dir / 'result.json'}")
    return payload


def _label_internal_nodes(tree: Tree) -> list[Clade]:
    rooted = _is_rooted_representation(tree)
    nodes = [node for node in tree.find_clades(order="preorder") if node.clades]
    if not rooted:
        nodes = [node for node in nodes if node is not tree.root]
    for index, node in enumerate(nodes, 1):
        node.name = f"N{index}"
        node.confidence = None
    return nodes


def _render_tree_pdf(tree: Tree, path: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(12, 8))
    Phylo.draw(tree, axes=ax, do_show=False)
    fig.savefig(path, format="pdf", bbox_inches="tight")
    plt.close(fig)


def run_label_nodes(
    tree: Path,
    output_dir: Path = Path("runs/posttree/syserror/brlen/label_nodes"),
    overwrite: bool = False,
    quiet: bool = False,
) -> dict[str, Any]:
    """Label internal nodes of a single tree and write tree/map/PDF artifacts."""
    start = time.time()
    if not tree.is_file():
        raise ValueError(f"tree file not found: {tree}")
    try:
        parsed = list(Phylo.parse(tree, "newick"))
    except Exception as exc:
        raise ValueError(f"tree file {tree.name}: failed to parse as Newick") from exc
    if len(parsed) != 1:
        raise ValueError("label-nodes accepts exactly one tree (single-mode only)")
    target = parsed[0]
    if len(target.get_terminals()) < 2:
        raise ValueError("tree must have at least two terminals")

    resolved_output_dir = output_dir.resolve()
    if resolved_output_dir.exists() and any(resolved_output_dir.iterdir()):
        if not overwrite:
            raise ValueError(
                f"Output directory exists and is not empty: {resolved_output_dir}\n"
                "Use --overwrite to replace."
            )
        shutil.rmtree(resolved_output_dir)
    resolved_output_dir.mkdir(parents=True, exist_ok=True)

    labeled_nodes = _label_internal_nodes(target)
    rooted = _is_rooted_representation(target)
    stem = tree.stem

    labeled_path = resolved_output_dir / f"{stem}.labeled.nwk"
    with open(labeled_path, "w") as fh:
        Phylo.write(target, fh, "newick", format_branch_length="%r")

    map_path = resolved_output_dir / f"{stem}.map.txt"
    map_lines = [f"{node.name}:{','.join(_terminal_names(node))}" for node in labeled_nodes]
    map_path.write_text("\n".join(map_lines) + "\n")

    pdf_path = resolved_output_dir / f"{stem}.labeled.pdf"
    _render_tree_pdf(target, pdf_path)

    payload = {
        "status": "success",
        "command": shlex.join(
            ["phyloai", "posttree", "syserror", "brlen", "label-nodes",
             "--tree", str(tree), "-o", str(output_dir)]
            + (["--overwrite"] if overwrite else [])
            + (["-q"] if quiet else [])
        ),
        "wall_time": time.time() - start,
        "tool_versions": {},
        "params": {
            "tree": str(tree),
            "output_dir": str(output_dir),
            "overwrite": overwrite,
            "quiet": quiet,
        },
        "key_results": {
            "n_internal_nodes_labeled": len(labeled_nodes),
            "n_terminals": len(target.get_terminals()),
            "rooted": rooted,
        },
        "error": None,
        "data": {
            "output_files": {
                "labeled_tree": {
                    "path": str(labeled_path),
                    "description": f"Newick tree with internal node labels N1..N{len(labeled_nodes)}",
                },
                "map_file": {
                    "path": str(map_path),
                    "description": "Node-species map template for brlen node modes",
                },
                "tree_figure": {
                    "path": str(pdf_path),
                    "description": "Tree visualization with labeled internal nodes",
                },
            }
        },
    }
    with open(resolved_output_dir / "result.json", "w") as fh:
        json.dump(payload, fh, indent=2)

    if not quiet:
        __import__("click").echo(
            f"Labeled {len(labeled_nodes)} internal nodes; wrote "
            f"{labeled_path}, {map_path}, {pdf_path}"
        )
    return payload
