"""Molecular marker MSA + tree metric computation for pretree workflows."""

from __future__ import annotations

import csv
import io
import json
import math
import shutil
import subprocess
import time
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from typing import Any

import numpy as np
from Bio import Phylo, SeqIO
from Bio.Phylo.BaseTree import Tree
from rich.progress import Progress

from phyloai.core.env import TOOL_REGISTRY, ToolEnv
from phyloai.core.sequence_normalization import (
    AA_STANDARD,
    NT_STANDARD,
    gap_chars,
    resolve_seq_type,
    standard_chars,
)

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

MSA_EXTENSIONS = {".fa", ".fasta", ".fas", ".fna", ".faa", ".aln"}
TREE_EXTENSIONS = {".tre", ".tree", ".nwk", ".newick", ".treefile", ".bestTree", ".contree"}
_TREE_SPECIFIC_SUFFIXES = (".treefile", ".contree", ".bestTree", ".iqtree", ".tree", ".tre", ".nwk", ".newick")
_GENERAL_SUFFIXES = (".fa", ".fasta", ".fas", ".fna", ".faa", ".aln")

_METRICS_CSV_ORDER = [
    "loci", "DataType",
    "num_taxa", "taxa_occupancy", "num_sites",
    "num_patterns", "proportion_patterns",
    "num_parsimony_sites", "proportion_parsimony",
    "num_singletons", "proportion_singletons",
    "num_sites/num_taxa", "num_patterns/num_taxa",
    "num_parsimony_sites/num_taxa", "num_singletons/num_taxa",
    "proportion_gaps", "proportion_invariant",
    "entropy", "bollback", "pattern_entropy",
    "rcfv", "nrcfv", "average_pairwise_identity", "GC_content",
    "average_BS", "sd_BS", "total_tree_length",
    "average_internal_branch_length", "sd_internal_branch_length",
    "average_terminal_branch_length", "sd_terminal_branch_length",
    "tree_diameter", "average_patristic_distance", "sd_patristic_distance",
    "evo_rate", "treeness", "dvmc", "saturation", "RF_distance",
]

_AA_FREQ_LABELS = list("ARNDCQEGHILKMFPSTWYV")
_NT_FREQ_LABELS = list("ACGT")

_PSEUDO_TREE_METRIC_NAMES = [
    "average_BS_FT", "sd_BS_FT", "total_tree_length_FT",
    "average_internal_branch_length_FT", "sd_internal_branch_length_FT",
    "average_terminal_branch_length_FT", "sd_terminal_branch_length_FT",
    "tree_diameter_FT", "average_patristic_distance_FT", "sd_patristic_distance_FT",
]


def _to_stdev(values: list[float]) -> float:
    """Sample standard deviation (ddof=1), via statistics.stdev."""
    import statistics
    return statistics.stdev(values) if len(values) > 1 else 0.0


def _to_mean(values: list[float]) -> float:
    return float(np.mean(values)) if values else 0.0


# ---------------------------------------------------------------------------
# Task 1: Helpers
# ---------------------------------------------------------------------------


def _read_msa(path: Path) -> list[list[str]]:
    records = list(SeqIO.parse(str(path), "fasta"))
    if not records:
        raise ValueError(f"No sequences found in '{path}'.")
    return [[c.upper() for c in str(r.seq)] for r in records]


def _read_msa_records(path: Path) -> list:
    records = list(SeqIO.parse(str(path), "fasta"))
    if not records:
        raise ValueError(f"No sequences found in '{path}'.")
    return records


def _scan_msa_headers(msa_dir: Path) -> tuple[dict[str, set[str]], set[str]]:
    per_marker: dict[str, set[str]] = {}
    total_pool: set[str] = set()
    for path in sorted(msa_dir.iterdir()):
        if path.suffix.lower() not in MSA_EXTENSIONS:
            continue
        try:
            records = _read_msa_records(path)
            taxa = {record.id for record in records}
            per_marker[path.stem] = taxa
            total_pool.update(taxa)
        except Exception:
            continue
    return per_marker, total_pool


def _strip_tree_suffixes(name: str) -> str:
    changed = True
    while changed:
        changed = False
        for suf in _TREE_SPECIFIC_SUFFIXES:
            if name.endswith(suf) and len(name) > len(suf):
                name = name[: -len(suf)]
                changed = True
                break
        else:
            for suf in _GENERAL_SUFFIXES:
                if name.endswith(suf) and len(name) > len(suf):
                    name = name[: -len(suf)]
                    changed = True
                    break
    return name


def _pair_files(
    msa_dir: Path | None,
    tree_dir: Path | None,
) -> tuple[dict[str, tuple[Path | None, Path | None]], list[str]]:
    warnings: list[str] = []
    msa_stems: dict[str, Path] = {}
    if msa_dir is not None:
        for path in sorted(msa_dir.iterdir()):
            if path.suffix.lower() in MSA_EXTENSIONS:
                msa_stems[path.stem] = path

    tree_map: dict[str, Path] = {}
    if tree_dir is not None:
        for path in sorted(tree_dir.iterdir()):
            if path.suffix.lower() in TREE_EXTENSIONS or path.name.endswith(_TREE_SPECIFIC_SUFFIXES):
                stem = _strip_tree_suffixes(path.name)
                tree_map[stem] = path

    all_stems = sorted(set(msa_stems.keys()) | set(tree_map.keys()))
    paired: dict[str, tuple[Path | None, Path | None]] = {}
    for stem in all_stems:
        paired[stem] = (msa_stems.get(stem), tree_map.get(stem))

    msa_only = sorted(set(msa_stems) - set(tree_map))
    tree_only = sorted(set(tree_map) - set(msa_stems))
    for stem in msa_only:
        warnings.append(f"[WARN] MSA file '{stem}' has no matching tree file.")
    for stem in tree_only:
        warnings.append(f"[WARN] Tree file '{stem}' has no matching MSA file.")

    return paired, warnings


def _check_taxon_consistency(msa_path: Path, tree_path: Path) -> dict | None:
    msa_taxa = {r.id for r in SeqIO.parse(str(msa_path), "fasta")}
    tree = Phylo.read(str(tree_path), "newick")
    tree_taxa = {term.name for term in tree.get_terminals()}
    msa_only = sorted(msa_taxa - tree_taxa)
    tree_only = sorted(tree_taxa - msa_taxa)
    if msa_only or tree_only:
        return {"msa_only": msa_only, "tree_only": tree_only}
    return None


def _resolve_seq_type(path: Path) -> str:
    records = list(SeqIO.parse(str(path), "fasta"))
    sequences = [str(r.seq).upper() for r in records]
    seq_type, _ = resolve_seq_type(sequences)
    return seq_type


# ---------------------------------------------------------------------------
# Task 2: MSA metrics — all formulas match extract_msa_tree_features.py
# ---------------------------------------------------------------------------


def _compute_msa_metrics(
    msa: list[list[str]],
    seq_type: str,
    total_taxa_pool: int,
    skip_freq: bool,
) -> dict[str, Any]:
    valid = standard_chars(seq_type)
    ntaxa = len(msa)
    nsites = len(msa[0]) if msa else 0
    result: dict[str, Any] = {
        "num_taxa": ntaxa,
        "taxa_occupancy": ntaxa / total_taxa_pool if total_taxa_pool > 0 else 0.0,
        "num_sites": nsites,
    }

    if nsites == 0:
        for key in _METRICS_CSV_ORDER:
            if key not in result:
                result[key] = 0.0 if "proportion" in key or "entropy" in key or "ratio" in key else ""
        return result

    std_codes = {ord(ch) for ch in valid}
    gap_code = ord("-")

    # --- Pass 1: site-level stats (patterns, parsimony, invariant, entropy using log2) ---
    pattern_set: set[bytes] = set()
    pattern_counter: Counter[bytes] = Counter()
    parsimony_info = 0
    singleton_sites = 0
    invariant_count = 0
    total_gaps = 0
    entropies: list[float] = []

    # Convert sequences to bytes once
    seq_bytes = ["".join(row).encode("ascii") for row in msa]

    for col in zip(*seq_bytes):
        norm = bytes(gap_code if ch not in std_codes else ch for ch in col)
        pattern_set.add(norm)
        pattern_counter[norm] += 1

        std_counts: dict[int, int] = {}
        for ch in col:
            if ch in std_codes:
                std_counts[ch] = std_counts.get(ch, 0) + 1
            else:
                total_gaps += 1

        n_std = sum(std_counts.values())
        unique_std = len(std_counts)

        # invariant: ≤1 unique valid char (matches extract_msa_tree_features.py)
        if unique_std <= 1:
            invariant_count += 1

        # parsimony-informative / singleton (matches stats.py / extract_msa_tree_features.py)
        if n_std >= 2 and unique_std > 1:
            repeated = sum(1 for v in std_counts.values() if v >= 2)
            if repeated >= 2:
                parsimony_info += 1
            else:
                singleton_sites += 1

        # Shannon entropy per site using log2 (matches extract_msa_tree_features.py)
        if n_std > 0:
            site_entropy = 0.0
            for cnt in std_counts.values():
                p = cnt / n_std
                if p > 0:
                    site_entropy -= p * math.log2(p)
            entropies.append(site_entropy)

    num_patterns = len(pattern_set)
    result["num_patterns"] = num_patterns
    result["proportion_patterns"] = num_patterns / nsites if nsites > 0 else 0.0
    result["num_parsimony_sites"] = parsimony_info
    result["proportion_parsimony"] = parsimony_info / nsites if nsites > 0 else 0.0
    result["num_singletons"] = singleton_sites
    result["proportion_singletons"] = singleton_sites / nsites if nsites > 0 else 0.0

    result["num_sites/num_taxa"] = nsites / ntaxa if ntaxa > 0 else 0.0
    result["num_patterns/num_taxa"] = num_patterns / ntaxa if ntaxa > 0 else 0.0
    result["num_parsimony_sites/num_taxa"] = parsimony_info / ntaxa if ntaxa > 0 else 0.0
    result["num_singletons/num_taxa"] = singleton_sites / ntaxa if ntaxa > 0 else 0.0

    total_cells = nsites * ntaxa
    result["proportion_gaps"] = total_gaps / total_cells if total_cells > 0 else 0.0
    result["proportion_invariant"] = invariant_count / nsites if nsites > 0 else 0.0

    # Mean per-site Shannon entropy (log2)
    result["entropy"] = _to_mean(entropies)

    # --- pattern_entropy / Bollback (matches extract_msa_tree_features.py) ---
    pe = 0.0
    for count in pattern_counter.values():
        if count > 0:
            pe += count * math.log(count)
    result["pattern_entropy"] = pe
    result["bollback"] = pe - nsites * math.log(nsites) if nsites > 0 else 0.0

    # --- RCFV (classic) — matches extract_msa_tree_features.py ---
    rcfv = 0.0
    for state in sorted(valid):
        if state == "-":
            continue
        state_freqs = []
        for seq in ["".join(row) for row in msa]:
            total_valid = sum(1 for c in seq if c in valid)
            if total_valid > 0:
                state_freqs.append(seq.count(state) / total_valid)
        if state_freqs and ntaxa > 0:
            mean_freq = sum(state_freqs) / ntaxa
            rcfv += sum(abs(f - mean_freq) for f in state_freqs) / ntaxa
    result["rcfv"] = rcfv

    # --- nRCFV (RCFV_Reader formula from https://github.com/JFFleming/RCFV_Reader) ---
    # DNA: nRCFV = RCFV * sqrt(L) / 400   (400 = 4 chars × 100)
    # Protein: nRCFV = RCFV * sqrt(L) / 2000  (2000 = 20 chars × 100)
    c_val = 4 if seq_type == "NT" else 20
    try:
        nrcfv = rcfv * math.sqrt(nsites) / (c_val * 100)
    except (ZeroDivisionError, OverflowError):
        nrcfv = 0.0
    result["nrcfv"] = nrcfv

    # --- average_pairwise_identity (matches extract_msa_tree_features.py) ---
    # All sites used: identity_count / alignment_length
    result["average_pairwise_identity"] = _compute_average_pairwise_identity_reference(
        ["".join(row) for row in msa]
    )

    # --- GC_content (NT only) ---
    if seq_type == "NT":
        gc = sum(s.count("G") + s.count("C") for s in ["".join(row) for row in msa])
        total_valid_nt = sum(sum(1 for c in s if c in valid) for s in ["".join(row) for row in msa])
        result["GC_content"] = gc / total_valid_nt if total_valid_nt > 0 else 0.0
    else:
        result["GC_content"] = ""

    # --- Frequency statistics (only relevant chars per seq_type) ---
    if not skip_freq:
        result.update(_compute_frequencies(["".join(row) for row in msa], seq_type))

    return result


def _compute_average_pairwise_identity_reference(sequences: list[str]) -> float:
    """Matches extract_msa_tree_features.py: identical positions / alignment_length."""
    ntaxa = len(sequences)
    if ntaxa < 2:
        return 0.0
    total_identity = 0.0
    pairs = 0
    for i in range(ntaxa):
        for j in range(i + 1, ntaxa):
            identical = sum(1 for a, b in zip(sequences[i], sequences[j]) if a == b)
            total_identity += identical / len(sequences[i])
            pairs += 1
    return total_identity / pairs if pairs > 0 else 0.0


def _compute_frequencies(sequences: list[str], seq_type: str) -> dict[str, float]:
    """Return only the frequency labels relevant for this seq_type."""
    valid = standard_chars(seq_type)
    labels = _NT_FREQ_LABELS if seq_type == "NT" else _AA_FREQ_LABELS
    combined = "".join(sequences)
    total_valid = sum(1 for c in combined if c in valid)
    freqs: dict[str, float] = {}
    for label in labels:
        cnt = combined.count(label)
        freqs[f"freq{label}"] = cnt / total_valid if total_valid > 0 else 0.0
    return freqs


# ---------------------------------------------------------------------------
# Task 3: Tree metrics — formulas match extract_msa_tree_features.py
# ---------------------------------------------------------------------------


def _compute_tree_metrics(
    tree_path: Path,
    outgroup_list: Path | None,
    ref_tree_path: Path | None,
) -> dict[str, Any]:
    tree = Phylo.read(str(tree_path), "newick")
    terminals = tree.get_terminals()
    n_taxa = len(terminals)
    internal = tree.get_nonterminals()

    # --- bootstrap support (population std = np.std, ddof=0) ---
    bs_values = []
    for clade in internal:
        if clade.confidence is not None and clade.confidence > 0:
            bs_values.append(float(clade.confidence))
    avg_bs = _to_mean(bs_values)
    sd_bs = _to_stdev(bs_values)

    # --- branch lengths ---
    all_brlen = [clade.branch_length for clade in tree.find_clades() if clade.branch_length is not None]
    terminal_brlen = [t.branch_length for t in terminals if t.branch_length is not None]
    internal_brlen = [i.branch_length for i in internal if i.branch_length is not None]
    total_len = sum(all_brlen)

    # --- patristic distances ---
    patristic = _patristic_distances(tree, terminals)

    result: dict[str, Any] = {
        "average_BS": avg_bs,
        "sd_BS": sd_bs,
        "total_tree_length": total_len,
        "average_internal_branch_length": _to_mean(internal_brlen),
        "sd_internal_branch_length": _to_stdev(internal_brlen),
        "average_terminal_branch_length": _to_mean(terminal_brlen),
        "sd_terminal_branch_length": _to_stdev(terminal_brlen),
        "tree_diameter": max(patristic) if patristic else 0.0,
        "average_patristic_distance": _to_mean(patristic),
        "sd_patristic_distance": _to_stdev(patristic),
        "evo_rate": total_len / n_taxa if n_taxa > 0 else 0.0,
        "treeness": sum(internal_brlen) / total_len if total_len > 0 else 0.0,
        "dvmc": _compute_dvmc(tree, outgroup_list),
        "saturation": "",
        "RF_distance": "",
    }

    if ref_tree_path is not None:
        ref_tree = Phylo.read(str(ref_tree_path), "newick")
        result["RF_distance"] = _compute_rf_distance(tree, ref_tree)

    return result


def _patristic_distances(tree: Tree, terminals: list) -> list[float]:
    dists = []
    n = len(terminals)
    for i in range(n):
        for j in range(i + 1, n):
            try:
                d = tree.distance(terminals[i], terminals[j])
                if d is not None and np.isfinite(d):
                    dists.append(float(d))
            except Exception:
                continue
    return dists


def _compute_dvmc(tree: Tree, outgroup_list: Path | None) -> float:
    """Matches extract_msa_tree_features.py and phykit: root-to-tip distance SD (population, ddof=0)."""
    import copy
    tree_copy = copy.deepcopy(tree)

    if outgroup_list is not None and outgroup_list.exists():
        outgroup_names = set(outgroup_list.read_text().strip().splitlines())
        terminal_names = {term.name for term in tree_copy.get_terminals()}
        outgroups_in_tree = outgroup_names & terminal_names
        for tip in list(tree_copy.get_terminals()):
            if tip.name in outgroups_in_tree:
                tree_copy.prune(tip)

    num_spp = tree_copy.count_terminals()
    if num_spp < 2:
        return 0.0

    try:
        sum_dist = 0.0
        sumi2N = 0.0
        for term in tree_copy.get_terminals():
            dist = tree_copy.distance(term)
            sum_dist += dist
            sumi2N += dist ** 2
        avg_dist = sum_dist / num_spp
        squared_diff_sum = sumi2N - num_spp * (avg_dist ** 2)
        return float(math.sqrt(squared_diff_sum / (num_spp - 1)))
    except Exception:
        return 0.0


def _compute_rf_distance(tree: Tree, ref_tree: Tree) -> float:
    tips1 = {t.name for t in tree.get_terminals()}
    tips2 = {t.name for t in ref_tree.get_terminals()}
    shared = tips1 & tips2
    if len(shared) < 4:
        return 0.0
    t1 = _root_and_prune(tree, shared)
    t2 = _root_and_prune(ref_tree, shared)
    clades1 = _get_clade_sets(t1)
    clades2 = _get_clade_sets(t2)
    plain_rf = len(clades1 - clades2) + len(clades2 - clades1)
    n_shared = len(shared)
    max_rf = 2 * (n_shared - 3)
    return plain_rf / max_rf if max_rf > 0 else 0.0


def _root_and_prune(tree: Tree, shared_tips: set[str]) -> Tree:
    import copy
    tree2 = copy.deepcopy(tree)
    for leaf in list(tree2.get_terminals()):
        if leaf.name not in shared_tips:
            tree2.prune(leaf)
    first_tip = sorted(shared_tips)[0]
    try:
        tree2.root_with_outgroup(first_tip)
    except Exception:
        pass
    return tree2


def _get_clade_sets(tree: Tree) -> set[frozenset[str]]:
    clades: set[frozenset[str]] = set()
    for clade in tree.find_clades():
        names = frozenset(t.name for t in clade.get_terminals())
        if 1 < len(names) < len(list(tree.get_terminals())):
            clades.add(names)
    return clades


# ---------------------------------------------------------------------------
# Saturation — matches phykit exactly with name-based lookup
# ---------------------------------------------------------------------------


def _compute_saturation(
    msa_path: Path,
    tree_path: Path,
    exclude_gaps: bool = False,
) -> float:
    """Slope of patristic vs uncorrected distance through origin (matches phykit)."""
    try:
        records = list(SeqIO.parse(str(msa_path), "fasta"))
        tree = Phylo.read(str(tree_path), "newick")
    except Exception:
        return 0.0

    tips = [term.name for term in tree.get_terminals()]
    n = len(tips)
    if n < 4:
        return 0.0

    # Build name-based lookup for MSA sequences (matches phykit order)
    name_to_seq: dict[str, str] = {}
    for rec in records:
        name_to_seq[rec.id] = str(rec.seq).upper()

    patristic = []
    uncorrected = []
    for i in range(n):
        for j in range(i + 1, n):
            try:
                pd = tree.distance(tips[i], tips[j])
                if pd is None or not np.isfinite(pd):
                    continue
            except Exception:
                continue
            patristic.append(float(pd))

            seq1 = name_to_seq.get(tips[i], "")
            seq2 = name_to_seq.get(tips[j], "")
            if not seq1 or not seq2:
                continue

            if exclude_gaps:
                # phykit: valid_positions = ~(gap_mask1 | gap_mask2)
                # i.e. positions where neither sequence has a gap character
                # Use full gap set like phykit (not just '-')
                try:
                    _st = _resolve_seq_type(msa_path)
                except Exception:
                    _st = "AA"
                gap_set = gap_chars(_st)
                total_len = 0
                matches = 0
                for a, b in zip(seq1, seq2):
                    if a not in gap_set and b not in gap_set:
                        total_len += 1
                        if a == b:
                            matches += 1
                ud = 1.0 - (matches / total_len) if total_len > 0 else 0.0
            else:
                matches = sum(1 for a, b in zip(seq1, seq2) if a == b)
                ud = 1.0 - matches / len(seq1) if len(seq1) > 0 else 0.0
            uncorrected.append(ud)

    if len(patristic) < 2 or len(uncorrected) < 2:
        return 0.0

    x = np.array(patristic, dtype=float)
    y = np.array(uncorrected, dtype=float)
    mask = np.isfinite(x) & np.isfinite(y)
    x, y = x[mask], y[mask]
    if len(x) < 2:
        return 0.0
    denom = float(np.dot(x, x))
    if denom == 0.0:
        return 0.0
    return float(np.dot(x, y) / denom)


# ---------------------------------------------------------------------------
# Task 4: Pseudo-tree (FastTree)
# ---------------------------------------------------------------------------


def _compute_pseudo_tree_metrics(
    msa_path: Path,
    seq_type: str,
    fasttree_path: str = "FastTree",
) -> dict[str, Any]:
    args = [fasttree_path]
    if seq_type == "AA":
        args.extend(["-lg", "-noml", "-boot", "500"])
    else:
        args.extend(["-nt", "-gtr", "-noml", "-boot", "500"])

    try:
        with open(msa_path) as fh:
            msa_content = fh.read()
        result = subprocess.run(
            args,
            input=msa_content.encode(),
            capture_output=True,
            timeout=120,
        )
        if result.returncode != 0 or not result.stdout.strip():
            return {k: "" for k in _PSEUDO_TREE_METRIC_NAMES}
        tree = Phylo.read(io.StringIO(result.stdout.decode()), "newick")
        ft_metrics = _extract_tree_features(tree)
        return {f"{k}_FT": v for k, v in ft_metrics.items()}
    except Exception:
        return {k: "" for k in _PSEUDO_TREE_METRIC_NAMES}


def _extract_tree_features(tree: Tree) -> dict[str, Any]:
    terminals = tree.get_terminals()
    internal = tree.get_nonterminals()

    bs_values = []
    for clade in internal:
        if clade.confidence is not None and clade.confidence > 0:
            bs_values.append(float(clade.confidence))

    all_brlen = [clade.branch_length for clade in tree.find_clades() if clade.branch_length is not None]
    terminal_brlen = [t.branch_length for t in terminals if t.branch_length is not None]
    internal_brlen = [i.branch_length for i in internal if i.branch_length is not None]
    total_len = sum(all_brlen)

    patristic = _patristic_distances(tree, terminals)

    return {
        "average_BS": _to_mean(bs_values),
        "sd_BS": _to_stdev(bs_values),
        "total_tree_length": total_len,
        "average_internal_branch_length": _to_mean(internal_brlen),
        "sd_internal_branch_length": _to_stdev(internal_brlen),
        "average_terminal_branch_length": _to_mean(terminal_brlen),
        "sd_terminal_branch_length": _to_stdev(terminal_brlen),
        "tree_diameter": max(patristic) if patristic else 0.0,
        "average_patristic_distance": _to_mean(patristic),
        "sd_patristic_distance": _to_stdev(patristic),
    }


# ---------------------------------------------------------------------------
# Task 5: Orchestration
# ---------------------------------------------------------------------------


def _metric_worker(args: tuple) -> dict:
    (
        stem, msa_path, tree_path, total_taxa_pool, seq_type, skip_freq,
        pseudo_tree, fasttree_bin, skip_pairwise_identity, outgroup_list,
        ref_tree_path, decimal_places,
    ) = args

    result: dict[str, Any] = {"loci": stem, "DataType": ""}

    try:
        if msa_path is not None:
            msa_data = _read_msa(msa_path)
            auto_type = _resolve_seq_type(msa_path)
            use_type = seq_type if seq_type != "auto" else auto_type
            result["DataType"] = use_type

            msa_metrics = _compute_msa_metrics(msa_data, use_type, total_taxa_pool, skip_freq)
            if skip_pairwise_identity:
                msa_metrics["average_pairwise_identity"] = ""
            result.update(msa_metrics)

            if pseudo_tree and fasttree_bin:
                ft = _compute_pseudo_tree_metrics(msa_path, use_type, fasttree_bin)
                result.update(ft)
            elif pseudo_tree:
                result.update({k: "" for k in _PSEUDO_TREE_METRIC_NAMES})

            if tree_path is not None:
                tree_metrics = _compute_tree_metrics(tree_path, outgroup_list, ref_tree_path)
                result.update(tree_metrics)
                sat = _compute_saturation(msa_path, tree_path)
                result["saturation"] = sat

        elif tree_path is not None:
            tree_metrics = _compute_tree_metrics(tree_path, outgroup_list, ref_tree_path)
            result.update(tree_metrics)
    except Exception as exc:
        result["_error"] = str(exc)

    for key, val in list(result.items()):
        if isinstance(val, float):
            result[key] = round(val, decimal_places)

    return result


def _build_csv_rows(
    results: list[dict],
    decimal_places: int,
    skip_freq: bool,
    pseudo_tree: bool,
) -> list[dict]:
    ordered = list(_METRICS_CSV_ORDER)

    if not skip_freq:
        # Determine which seq_types are actually present
        data_types_present = set()
        for row in results:
            dt = row.get("DataType", "")
            if dt in ("AA", "NT"):
                data_types_present.add(dt)

        if "AA" in data_types_present:
            for lbl in _AA_FREQ_LABELS:
                name = f"freq{lbl}"
                if name not in ordered:
                    ordered.append(name)
        if "NT" in data_types_present:
            for lbl in _NT_FREQ_LABELS:
                name = f"freq{lbl}"
                if name not in ordered:
                    ordered.append(name)

    if pseudo_tree:
        for name in _PSEUDO_TREE_METRIC_NAMES:
            if name not in ordered:
                ordered.append(name)

    seen_cols = set(ordered)
    for row in results:
        for key in sorted(row.keys()):
            if key not in seen_cols and not key.startswith("_"):
                ordered.append(key)
                seen_cols.add(key)

    rows = []
    for row in results:
        csv_row = {}
        for col in ordered:
            val = row.get(col)
            if val is None:
                csv_row[col] = ""
            elif isinstance(val, float):
                csv_row[col] = round(val, decimal_places)
            else:
                csv_row[col] = val
        rows.append(csv_row)
    return rows


def _drop_empty_columns(rows: list[dict]) -> list[dict]:
    """Remove columns where every row has an empty value ('' or None).
    Excludes 'loci' and 'DataType' from removal — these are identifier columns."""
    if not rows:
        return rows
    all_keys = list(rows[0].keys())
    to_drop = set()
    for key in all_keys:
        if key in ("loci", "DataType"):
            continue
        if all(r.get(key) in (None, "") for r in rows):
            to_drop.add(key)
    if not to_drop:
        return rows
    return [{k: v for k, v in r.items() if k not in to_drop} for r in rows]


def _write_metrics_csv(rows: list[dict], output_path: Path) -> None:
    if not rows:
        return
    rows = _drop_empty_columns(rows)
    fieldnames = list(rows[0].keys())
    with open(output_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def _write_result_json(payload: dict, output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    with open(output_dir / "result.json", "w") as fh:
        json.dump(payload, fh, indent=2)


def _write_log(log_data: dict, output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    with open(output_dir / "metrics.log", "w") as fh:
        json.dump(log_data, fh, indent=2)


def run_metrics(
    msa_dir: Path | None = None,
    tree_dir: Path | None = None,
    seq_type: str = "auto",
    threads: int = 4,
    output_dir: Path = Path("runs/pretree/metrics"),
    decimal_places: int = 6,
    skip_freq: bool = False,
    pseudo_tree: bool = False,
    fasttree_path: str = "FastTree",
    skip_pairwise_identity: bool = False,
    outgroup_list: Path | None = None,
    ref_tree: Path | None = None,
    overwrite: bool = False,
    dry_run: bool = False,
    quiet: bool = False,
    progress: Progress | None = None,
    console: Any = None,
) -> dict:
    t0 = time.monotonic()
    per_marker_stderr: list[str] = []

    if not msa_dir and not tree_dir:
        return {
            "status": "error",
            "command": "phyloai pretree metrics",
            "wall_time": 0.0,
            "tool_versions": {},
            "params": {},
            "key_results": {},
            "error": "At least one of --msa-dir or --tree-dir must be provided.",
            "data": {},
        }

    if output_dir.exists() and any(output_dir.iterdir()) and not overwrite:
        return {
            "status": "error",
            "command": "phyloai pretree metrics",
            "wall_time": time.monotonic() - t0,
            "tool_versions": {},
            "params": {},
            "key_results": {},
            "error": f"Output directory '{output_dir}' already exists and is non-empty. Use --overwrite to replace it.",
            "data": {},
        }

    paired, pair_warnings = _pair_files(msa_dir, tree_dir)
    per_marker_stderr.extend(pair_warnings)

    if not paired:
        return {
            "status": "error",
            "command": "phyloai pretree metrics",
            "wall_time": time.monotonic() - t0,
            "tool_versions": {},
            "params": {},
            "key_results": {},
            "error": "No paired MSA/tree files found.",
            "data": {},
        }

    total_taxa_pool = 0
    if msa_dir is not None:
        _, taxon_pool = _scan_msa_headers(msa_dir)
        total_taxa_pool = len(taxon_pool)

    if dry_run:
        stems = list(paired.keys())
        n_markers = len(stems)
        plan = {
            "n_markers": n_markers,
            "stems": stems[:10] if n_markers > 10 else stems,
            "pseudo_tree": pseudo_tree,
            "skip_freq": skip_freq,
            "skip_pairwise_identity": skip_pairwise_identity,
            "threads": threads,
            "output_dir": str(output_dir),
            "decimal_places": decimal_places,
        }
        return {
            "status": "success",
            "command": "phyloai pretree metrics",
            "wall_time": time.monotonic() - t0,
            "tool_versions": {},
            "params": {},
            "key_results": {"dry_run": True, **plan},
            "error": None,
            "data": {},
        }

    if overwrite and output_dir.exists():
        shutil.rmtree(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Pairwise identity warning
    if not skip_pairwise_identity and msa_dir is not None:
        large_markers = []
        for stem, (msa_p, _) in paired.items():
            if msa_p is not None:
                records = _read_msa_records(msa_p)
                if len(records) > 200:
                    large_markers.append(stem)
        if large_markers:
            msg = f"[WARN] {len(large_markers)} marker(s) have >200 taxa. Consider --skip-pairwise-identity."
            per_marker_stderr.append(msg)
            if console:
                console.print(f"[yellow]{msg}[/yellow]")

    # Taxon consistency checks
    for stem, (msa_p, tree_p) in paired.items():
        if msa_p is not None and tree_p is not None:
            mismatch = _check_taxon_consistency(msa_p, tree_p)
            if mismatch:
                msa_only = mismatch["msa_only"]
                tree_only = mismatch["tree_only"]
                if msa_only:
                    per_marker_stderr.append(f"[WARN] {stem}: taxa in MSA but not tree: {msa_only}")
                if tree_only:
                    per_marker_stderr.append(f"[WARN] {stem}: taxa in tree but not MSA: {tree_only}")

    worker_args = []
    for stem, (msa_p, tree_p) in paired.items():
        worker_args.append((
            stem, msa_p, tree_p, total_taxa_pool, seq_type, skip_freq,
            pseudo_tree, fasttree_path, skip_pairwise_identity, outgroup_list,
            ref_tree, decimal_places,
        ))

    results = []
    desc = "[cyan]Computing metrics..."
    if progress:
        task = progress.add_task(desc, total=len(worker_args))

    do_parallel = threads > 1 and len(worker_args) > 1
    if do_parallel:
        with ProcessPoolExecutor(max_workers=threads) as ex:
            futures = [ex.submit(_metric_worker, arg) for arg in worker_args]
            for future in futures:
                try:
                    results.append(future.result())
                except Exception as exc:
                    results.append({"_error": str(exc), "loci": "unknown"})
                if progress:
                    progress.advance(task)
    else:
        for arg in worker_args:
            results.append(_metric_worker(arg))
            if progress:
                progress.advance(task)

    if progress and not do_parallel:
        progress.update(task, completed=len(worker_args), visible=False)

    rows = _build_csv_rows(results, decimal_places, skip_freq, pseudo_tree)
    _write_metrics_csv(rows, output_dir / "metrics.csv")

    n_success = sum(1 for r in results if "_error" not in r)
    n_errors = sum(1 for r in results if "_error" in r)
    errors_list = [
        {"loci": r.get("loci", "?"), "error": r["_error"]}
        for r in results if "_error" in r
    ]

    wall_time = time.monotonic() - t0

    log_data = {
        "command": "phyloai pretree metrics",
        "params": {
            "msa_dir": str(msa_dir) if msa_dir else None,
            "tree_dir": str(tree_dir) if tree_dir else None,
            "seq_type": seq_type,
            "threads": threads,
            "decimal_places": decimal_places,
            "skip_freq_statistics": skip_freq,
            "pseudo_tree_metrics": pseudo_tree,
            "fasttree_path": fasttree_path if pseudo_tree else None,
            "skip_pairwise_identity": skip_pairwise_identity,
            "outgroup_list": str(outgroup_list) if outgroup_list else None,
            "ref_tree": str(ref_tree) if ref_tree else None,
            "overwrite": overwrite,
        },
        "wall_time": wall_time,
        "exit_code": 0,
        "per_marker_stderr": per_marker_stderr,
    }
    _write_log(log_data, output_dir)

    payload = {
        "status": "success" if n_errors == 0 else "partial",
        "command": "phyloai pretree metrics",
        "wall_time": wall_time,
        "tool_versions": {},
        "params": log_data["params"],
        "key_results": {
            "n_markers": len(paired),
            "n_success": n_success,
            "n_errors": n_errors,
            "errors": errors_list if errors_list else None,
            "warnings": per_marker_stderr if per_marker_stderr else None,
        },
        "error": None,
        "data": {},
    }
    _write_result_json(payload, output_dir)

    return payload


# ---------------------------------------------------------------------------
# Task 6: Distribution plots
# ---------------------------------------------------------------------------


def _plot_single_metric(
    data: np.ndarray,
    metric_name: str,
    output_path: Path,
    bins: int = 50,
    xmin: float | None = None,
    xmax: float | None = None,
    tukey_k: float | None = None,
    raw_labels: list[str] | None = None,
    title: str | None = None,
    xlabel: str | None = None,
    ylabel: str = "Density",
    color: str = "#2E86AB",
    fig_width: float = 10.0,
    fig_height: float = 8.0,
    dpi: int = 150,
    font_size: int = 12,
) -> tuple[int, list[tuple[str, float]]]:
    """Generate a density-distribution histogram PDF.

    Returns (n_filtered, [(label, value), ...]) for Tukey-filtered points.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from scipy.stats import gaussian_kde

    clean = data[np.isfinite(data)]
    if len(clean) == 0:
        return 0, []

    filtered_pairs: list[tuple[str, float]] = []

    if tukey_k is not None and raw_labels:
        q1 = np.percentile(clean, 25)
        q3 = np.percentile(clean, 75)
        iqr = q3 - q1
        lower = q1 - tukey_k * iqr
        upper = q3 + tukey_k * iqr
        for i, val in enumerate(clean):
            if val < lower or val > upper:
                lbl = raw_labels[i] if i < len(raw_labels) else ""
                filtered_pairs.append((lbl, float(val)))
        clean = clean[(clean >= lower) & (clean <= upper)]

    if len(clean) == 0:
        return len(filtered_pairs), filtered_pairs

    _title = title or f"Distribution of {_metric_display_name(metric_name)}"
    _xlabel = xlabel or _metric_display_name(metric_name)

    fig, ax = plt.subplots(figsize=(fig_width, fig_height))
    ax.hist(clean, bins=bins, density=True, alpha=0.6, color=color, edgecolor="white")
    try:
        kde = gaussian_kde(clean)
        x_range = np.linspace(clean.min(), clean.max(), 200)
        ax.plot(x_range, kde(x_range), color="#FF8C00", linewidth=2)
    except Exception:
        pass

    if xmin is not None:
        ax.set_xlim(left=xmin)
    if xmax is not None:
        ax.set_xlim(right=xmax)
    ax.set_title(_title, fontsize=font_size + 2)
    ax.set_xlabel(_xlabel, fontsize=font_size)
    ax.set_ylabel(ylabel, fontsize=font_size)
    ax.tick_params(labelsize=font_size - 2)

    fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return len(filtered_pairs), filtered_pairs


def _metric_display_name(name: str) -> str:
    """Convert snake_case metric name to human-readable display name."""
    display = {
        "num_taxa": "Number of Taxa",
        "taxa_occupancy": "Taxa Occupancy Ratio",
        "num_sites": "Number of Sites",
        "num_patterns": "Number of Patterns",
        "proportion_patterns": "Proportion of Patterns",
        "num_parsimony_sites": "Number of Parsimony-informative Sites",
        "proportion_parsimony": "Proportion of Parsimony-informative Sites",
        "num_singletons": "Number of Singleton Sites",
        "proportion_singletons": "Proportion of Singleton Sites",
        "num_sites/num_taxa": "Sites per Taxon",
        "num_patterns/num_taxa": "Patterns per Taxon",
        "num_parsimony_sites/num_taxa": "Parsimony Sites per Taxon",
        "num_singletons/num_taxa": "Singletons per Taxon",
        "proportion_gaps": "Gap Proportion",
        "proportion_invariant": "Invariant Site Proportion",
        "entropy": "Shannon Entropy",
        "bollback": "Bollback Multinomial Deviance",
        "pattern_entropy": "Pattern Entropy",
        "rcfv": "Relative Composition Frequency Variability (RCFV)",
        "nrcfv": "Normalized RCFV",
        "average_pairwise_identity": "Average Pairwise Identity",
        "GC_content": "GC Content",
        "average_BS": "Mean Bootstrap Support",
        "sd_BS": "SD Bootstrap Support",
        "total_tree_length": "Total Tree Length",
        "average_internal_branch_length": "Mean Internal Branch Length",
        "sd_internal_branch_length": "SD Internal Branch Length",
        "average_terminal_branch_length": "Mean Terminal Branch Length",
        "sd_terminal_branch_length": "SD Terminal Branch Length",
        "tree_diameter": "Tree Diameter",
        "average_patristic_distance": "Mean Patristic Distance",
        "sd_patristic_distance": "SD Patristic Distance",
        "evo_rate": "Evolutionary Rate",
        "treeness": "Treeness",
        "dvmc": "DVMC",
        "saturation": "Saturation Slope",
        "RF_distance": "Normalized RF Distance",
    }
    for freq_label in _AA_FREQ_LABELS + _NT_FREQ_LABELS:
        key = f"freq{freq_label}"
        if key not in display:
            display[key] = f"Frequency of {freq_label}"
    for pt_name in _PSEUDO_TREE_METRIC_NAMES:
        base = pt_name[:-3]  # strip _FT suffix
        display[pt_name] = display.get(base, base) + " (FastTree)"
    return display.get(name, name)



def _generate_all_plots(
    rows: list[dict],
    numeric_cols: list[str],
    plots_dir: Path,
    bins: int = 50,
) -> int:
    plots_dir.mkdir(parents=True, exist_ok=True)
    count = 0
    for col in numeric_cols:
        if col in ("loci", "DataType"):
            continue
        values = []
        labels = []
        for r in rows:
            v = r.get(col)
            lbl = r.get("loci", "")
            if v not in (None, "", "NA"):
                try:
                    values.append(float(v))
                    labels.append(lbl)
                except (ValueError, TypeError):
                    continue
        data = np.array(values, dtype=float)
        clean = data[np.isfinite(data)]
        if len(clean) == 0:
            continue
        safe_name = col.replace("/", "_")
        _plot_single_metric(clean, col, plots_dir / f"{safe_name}.pdf", bins=bins)
        count += 1
    return count


def _generate_basic_statistics(
    rows: list[dict],
    numeric_cols: list[str],
    output_path: Path,
) -> None:
    stats_rows = []
    for col in numeric_cols:
        values = []
        for r in rows:
            v = r.get(col)
            if v not in (None, "", "NA"):
                try:
                    values.append(float(v))
                except (ValueError, TypeError):
                    continue
        arr = np.array(values, dtype=float)
        clean = arr[np.isfinite(arr)]
        n_total = len(values)
        n_ex_na = len(clean)
        if n_ex_na == 0:
            stats_rows.append({
                "metric": col, "mean": "", "median": "", "min": "", "max": "",
                "q25": "", "q75": "", "std": "", "n_ex_NA": 0, "n_total": n_total,
            })
        else:
            stats_rows.append({
                "metric": col,
                "mean": float(np.mean(clean)),
                "median": float(np.median(clean)),
                "min": float(np.min(clean)),
                "max": float(np.max(clean)),
                "q25": float(np.percentile(clean, 25)),
                "q75": float(np.percentile(clean, 75)),
                "std": float(np.std(clean, ddof=1)) if len(clean) > 1 else 0.0,
                "n_ex_NA": n_ex_na,
                "n_total": n_total,
            })

    fieldnames = ["metric", "mean", "median", "min", "max", "q25", "q75", "std", "n_ex_NA", "n_total"]
    with open(output_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(stats_rows)


# ---------------------------------------------------------------------------
# Task 7: Correlation
# ---------------------------------------------------------------------------


def _compute_correlation(
    rows: list[dict],
    columns: list[str],
    method: str = "spearman",
) -> tuple[np.ndarray, list[str]]:
    from scipy.stats import pearsonr, spearmanr

    valid_cols = []
    data_cols = []
    for col in columns:
        values = []
        for r in rows:
            v = r.get(col)
            if v not in (None, "", "NA"):
                try:
                    values.append(float(v))
                except (ValueError, TypeError):
                    values.append(np.nan)
            else:
                values.append(np.nan)
        arr = np.array(values, dtype=float)
        if np.all(np.isnan(arr)):
            continue
        if np.nanstd(arr) == 0.0:
            continue
        valid_cols.append(col)
        data_cols.append(arr)

    if len(valid_cols) < 2:
        return np.array([]), []

    data = np.column_stack(data_cols)
    mask = ~np.any(np.isnan(data), axis=1)
    if not mask.any():
        return np.array([]), []
    data = data[mask, :]

    corr_func = spearmanr if method == "spearman" else pearsonr
    m = data.shape[1]
    corr_matrix = np.ones((m, m))

    if method == "pearson":
        mean = np.nanmean(data, axis=0)
        std = np.nanstd(data, axis=0, ddof=1)
        std[std == 0] = 1.0
        data = (data - mean) / std

    for i in range(m):
        for j in range(i + 1, m):
            try:
                r, _ = corr_func(data[:, i], data[:, j])
                corr_matrix[i, j] = corr_matrix[j, i] = float(r) if np.isfinite(r) else np.nan
            except Exception:
                corr_matrix[i, j] = corr_matrix[j, i] = np.nan

    return corr_matrix, valid_cols


def _select_correlation_columns(
    rows: list[dict],
    columns: list[str],
    requested: str | None = None,
    include_freq: bool = False,
    include_sd: bool = False,
) -> list[str]:
    """Select readable default columns for correlation plots."""
    if requested:
        if requested.strip().lower() == "all":
            return [col for col in columns if _column_has_numeric_values(rows, col)]
        return [col.strip() for col in requested.split(",") if col.strip()]

    selected: list[str] = []
    for col in columns:
        if col in ("loci", "DataType"):
            continue
        if col.startswith("freq") and not include_freq:
            continue
        if col.startswith("sd_") and not include_sd:
            continue
        if _column_has_numeric_values(rows, col):
            selected.append(col)
    return selected


def _column_has_numeric_values(rows: list[dict], col: str) -> bool:
    for row in rows:
        value = row.get(col)
        if value in (None, "", "NA"):
            continue
        try:
            numeric = float(value)
        except (TypeError, ValueError):
            continue
        if np.isfinite(numeric):
            return True
    return False


def _generate_correlation_heatmap(
    corr_matrix: np.ndarray,
    col_names: list[str],
    output_path: Path,
    triangle: str = "full",
    cluster_rectangles: int | None = None,
    cmap: str = "RdBu_r",
    annot: bool = False,
    fmt: str = ".2f",
    fig_width: float = 12.0,
    fig_height: float = 10.0,
    dpi: int = 150,
    font_size: int = 10,
    title: str | None = None,
    label_angle: float = 45.0,
    warn: Any | None = None,
) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    m = corr_matrix.shape[0]
    if m < 2:
        return

    linkage_matrix = None
    if m <= 2:
        order = np.arange(m)
    else:
        from scipy.cluster.hierarchy import linkage, leaves_list
        from scipy.spatial.distance import squareform

        clustering_matrix = np.nan_to_num(corr_matrix, nan=0.0)
        dist = 1.0 - np.abs(clustering_matrix)
        np.fill_diagonal(dist, 0.0)
        condensed = squareform(dist, checks=False)
        linkage_matrix = linkage(condensed, method="ward")
        order = leaves_list(linkage_matrix)

    # Reorder matrix and names by clustered leaf order without drawing dendrograms.
    reordered = corr_matrix[order, :][:, order]
    reordered_names = [col_names[i] for i in order]
    mask = None
    if triangle == "lower":
        mask = np.triu(np.ones_like(reordered, dtype=bool), k=1)
    elif triangle == "upper":
        mask = np.tril(np.ones_like(reordered, dtype=bool), k=-1)

    fig, ax = plt.subplots(figsize=(fig_width, fig_height))
    _draw_heatmap(ax, reordered, reordered_names, mask=mask, cmap=cmap,
                  annot=annot, fmt=fmt, font_size=font_size,
                  colorbar_side="left" if triangle == "upper" else "right")

    if title is not None:
        ax.set_title(title, fontsize=font_size + 4, pad=10)
    ax.set_xticklabels(ax.get_xticklabels(), fontsize=font_size, rotation=label_angle, ha="right")
    ax.set_yticklabels(ax.get_yticklabels(), fontsize=font_size, rotation=0)
    _style_triangle_axes(ax, triangle, reordered.shape[0], label_angle=label_angle)

    if triangle != "full" and cluster_rectangles is not None and cluster_rectangles > 0 and warn is not None:
        warn("--cluster-rectangles is ignored unless --triangle full is used.")

    if triangle == "full" and linkage_matrix is not None and cluster_rectangles is not None and cluster_rectangles > 0:
        from scipy.cluster.hierarchy import fcluster
        from matplotlib.patches import Rectangle
        clusters = fcluster(linkage_matrix, t=cluster_rectangles, criterion="maxclust")
        position_by_original = {original_idx: pos for pos, original_idx in enumerate(order)}
        unique_clusters = sorted(set(clusters))
        for cluster_id in unique_clusters:
            positions = sorted(position_by_original[i] for i, c in enumerate(clusters) if c == cluster_id)
            if len(positions) > 1:
                start = positions[0] - 0.5
                width = positions[-1] - positions[0] + 1.0
                rect = Rectangle((start, start), width, width, fill=False, edgecolor="black", linewidth=2)
                ax.add_patch(rect)

    fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def _draw_heatmap(ax, corr_matrix, col_names, *, mask=None, cmap="RdBu_r",
                   annot=False, fmt=".2f", font_size=10, colorbar_side="right"):
    """Draw a corrplot-style circle heatmap."""
    import matplotlib.pyplot as plt
    from matplotlib.collections import LineCollection

    m = corr_matrix.shape[0]

    # Circle area is proportional to |correlation|, capped at cell size.
    max_radius = 0.45
    cm = plt.get_cmap(cmap)
    norm = plt.Normalize(vmin=-1, vmax=1)
    for i in range(m):
        for j in range(m):
            if mask is not None and mask[i, j]:
                continue
            val = corr_matrix[i, j]
            if np.isnan(val):
                continue
            r = max_radius * np.sqrt(abs(val))
            color = cm(norm(val))
            circle = plt.Circle((j, i), r, facecolor=color, edgecolor="white",
                                linewidth=0.3, zorder=2)
            ax.add_patch(circle)
            if annot:
                ax.text(j, i, format(val, fmt), ha="center", va="center",
                        fontsize=max(font_size - 2, 6), color="black", zorder=3)
    ax.set_xlim(-0.5, m - 0.5)
    ax.set_ylim(m - 0.5, -0.5)
    ax.set_xticks(range(m))
    ax.set_yticks(range(m))
    ax.set_xticklabels(col_names, fontsize=font_size, rotation=45, ha="right")
    ax.set_yticklabels(col_names, fontsize=font_size, rotation=0)
    segments = []
    for i in range(m):
        for j in range(m):
            if mask is not None and mask[i, j]:
                continue
            left, right = j - 0.5, j + 0.5
            top, bottom = i - 0.5, i + 0.5
            segments.extend([
                [(left, top), (right, top)],
                [(right, top), (right, bottom)],
                [(right, bottom), (left, bottom)],
                [(left, bottom), (left, top)],
            ])
    ax.add_collection(LineCollection(segments, colors="#E5E5E5", linewidths=0.6, zorder=1))
    ax.tick_params(which="minor", bottom=False, left=False)
    sm = plt.cm.ScalarMappable(cmap=cm, norm=norm)
    sm.set_array([])
    location = "left" if colorbar_side == "left" else "right"
    ax.figure.colorbar(sm, ax=ax, fraction=0.035, pad=0.12,
                       location=location, label="Correlation")
    ax.set_aspect("equal")


def _style_triangle_axes(ax, triangle: str, matrix_size: int, label_angle: float = 45.0) -> None:
    """Align triangle labels and borders with the visible half of the matrix."""
    if triangle == "upper":
        ax.xaxis.tick_top()
        ax.xaxis.set_label_position("top")
        ax.yaxis.tick_right()
        ax.yaxis.set_label_position("right")
        ax.tick_params(axis="x", labeltop=True, labelbottom=False)
        ax.tick_params(axis="y", labelright=True, labelleft=False)
        for label in ax.get_xticklabels():
            label.set_rotation(label_angle)
            label.set_ha("left")
            label.set_va("bottom")
    elif triangle == "lower":
        ax.xaxis.tick_bottom()
        ax.xaxis.set_label_position("bottom")
        ax.yaxis.tick_left()
        ax.yaxis.set_label_position("left")
        ax.tick_params(axis="x", labelbottom=True, labeltop=False)
        ax.tick_params(axis="y", labelleft=True, labelright=False)
    else:
        return

    for spine in ax.spines.values():
        spine.set_visible(False)

    n = matrix_size
    if triangle == "lower":
        points = [(-0.5, -0.5), (-0.5, n - 0.5), (n - 0.5, n - 0.5)]
        for k in range(n - 1, -1, -1):
            points.append((k + 0.5, k - 0.5))
            points.append((k - 0.5, k - 0.5))
    else:
        points = [(-0.5, -0.5), (n - 0.5, -0.5), (n - 0.5, n - 0.5)]
        for k in range(n - 1, -1, -1):
            points.append((k - 0.5, k + 0.5))
            points.append((k - 0.5, k - 0.5))
    for (x1, y1), (x2, y2) in zip(points, points[1:]):
        ax.plot([x1, x2], [y1, y2], color="black", linewidth=0.8, clip_on=False, zorder=4)


def _write_correlation_csv(
    corr_matrix: np.ndarray,
    col_names: list[str],
    output_path: Path,
) -> None:
    if corr_matrix.size == 0:
        return
    with open(output_path, "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow([""] + col_names)
        for i, name in enumerate(col_names):
            writer.writerow([name] + [float(corr_matrix[i, j]) for j in range(corr_matrix.shape[1])])
