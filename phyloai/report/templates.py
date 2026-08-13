"""Per-command methods text generators for phyloai report.

Each step_id maps to a dedicated function that reads scientifically
meaningful parameters from params, key_results, and tool_versions,
producing 2-5 sentences of academic English suitable for journal Methods.
"""

from __future__ import annotations

import shlex
from pathlib import Path, PureWindowsPath
from typing import Any

_ALIGN_METHOD_MAP: dict[str, tuple[str, str]] = {
    "linsi": (
        "L-INS-i",
        "applies iterative local pairwise alignment refinement and is suited "
        "for sequences with conserved domains and insertions",
    ),
    "einsi": (
        "E-INS-i",
        "uses multiple local alignments and is suited for sequences with "
        "multiple conserved regions separated by unalignable regions",
    ),
    "ginsi": (
        "G-INS-i",
        "applies global pairwise alignment and is suited for sequences of "
        "similar length without large insertions",
    ),
    "fftns1": (
        "FFT-NS-1",
        "uses progressive alignment with single FFT iteration and is suited "
        "for large datasets where speed is prioritized",
    ),
    "fftns2": (
        "FFT-NS-2",
        "uses progressive alignment with two FFT iterations",
    ),
    "auto": (
        "auto-selected",
        "strategy selected automatically by MAFFT based on sequence length and count",
    ),
    "magus": (
        "MAGUS",
        "uses graph-based divide-and-conquer alignment and is suited for "
        "very large or highly divergent datasets",
    ),
}

_CLOCK_MAP: dict[str, str] = {
    "strict": "strict",
    "independent": "independent-rates",
    "correlated": "autocorrelated-rates",
}


def _describe_n(value: Any, singular: str, plural: str | None = None) -> str:
    try:
        n = int(value)
    except (TypeError, ValueError):
        return f"{value} {plural or singular + 's'}"
    if n == 1:
        return f"1 {singular}"
    return f"{n:,} {plural or singular + 's'}"


def _safe_fmt(value: Any, fmt_spec: str) -> str:
    """Format a value safely — returns '?' for non-numeric values."""
    if isinstance(value, (int, float)):
        return format(value, fmt_spec)
    return "?"


def _tool_arg_value(flag: str, tool_args: str | None) -> str | None:
    """Return the value immediately following a flag in --tool-args, if present.

    Uses shlex.split() to match ml_iqtree._get_tool_arg_value() so quoted values
    containing spaces (e.g. a custom matrix path like -m '/tmp/CUSTOM matrix')
    are parsed as one token instead of being truncated at the first space.
    """
    if not tool_args:
        return None
    try:
        tokens = shlex.split(tool_args)
    except ValueError:
        return None
    try:
        idx = tokens.index(flag)
    except ValueError:
        return None
    if idx + 1 >= len(tokens):
        return None
    return tokens[idx + 1]


# ---------------------------------------------------------------------------
# Per-step methods generators
# ---------------------------------------------------------------------------

def generate_methods_pretree_convert(
    params: dict[str, Any],
    key_results: dict[str, Any],
    tool_versions: dict[str, Any],
) -> str:
    to_fmt = params.get("to", "FASTA")
    n = key_results.get("n_converted", 0)
    n_skipped = key_results.get("n_skipped") or key_results.get("n_failed", 0)
    parts = [
        f"Raw sequence files were converted to {str(to_fmt).upper()} format "
        f"using phyloai pretree convert. "
        f"A total of {_describe_n(n, 'file')} were successfully converted",
    ]
    if n_skipped:
        parts.append(f" ({_describe_n(n_skipped, 'file')} skipped).")
    else:
        parts.append(".")
    return "".join(parts)


def generate_methods_pretree_stats(
    params: dict[str, Any],
    key_results: dict[str, Any],
    tool_versions: dict[str, Any],
) -> str:
    n = key_results.get("n_genes") or key_results.get("n_genes_ok") or key_results.get("n_taxa", 0)
    n_errors = key_results.get("n_errors", 0)
    format_name = key_results.get("format", params.get("input_format", "detected"))
    total_taxa = key_results.get("total_taxa") or key_results.get("n_taxa")
    seq_type = key_results.get("seq_type") or params.get("seq_type", "")
    unit = "sites" if str(seq_type).upper() == "AA" else "bp"

    parts = [
        f"Sequence statistics were computed using phyloai pretree stats "
        f"for {_describe_n(n, 'sequence file')} "
        f"in {str(format_name).upper() if isinstance(format_name, str) else str(format_name)} format.",
    ]
    if seq_type:
        parts.append(f" All sequences are {str(seq_type).upper()}.")
    if n_errors:
        parts.append(f" {_describe_n(n_errors, 'file')} encountered errors during processing.")
    if total_taxa:
        parts.append(f" A total of {_safe_fmt(total_taxa, ',')} unique taxa were identified across all files.")

    if key_results.get("length_mean"):
        parts.append(
            f" Mean alignment length was {_safe_fmt(key_results.get('length_mean'), '.1f')} {unit}; "
            f"taxa per locus ranged from {key_results.get('taxa_per_gene_min', '?')} "
            f"to {key_results.get('taxa_per_gene_max', '?')} "
        )
        if key_results.get("taxa_per_gene_mean"):
            parts[-1] += f"(mean {_safe_fmt(key_results.get('taxa_per_gene_mean'), '.1f')})."

    if key_results.get("gap_ratio_mean") is not None:
        parts[-1] += (
            f" Mean gap ratio was {_safe_fmt(key_results.get('gap_ratio_mean'), '.1%')}"
        )
        if key_results.get("gap_ratio_median") is not None:
            parts[-1] += f" (median {_safe_fmt(key_results.get('gap_ratio_median'), '.1%')})"
        parts[-1] += "."
    elif not parts[-1].endswith("."):
        parts[-1] += "."

    if key_results.get("ambiguous_ratio_mean") is not None:
        parts.append(
            f" Mean ambiguous character ratio was {_safe_fmt(key_results.get('ambiguous_ratio_mean'), '.1%')}"
        )
        if key_results.get("ambiguous_ratio_median") is not None:
            parts[-1] += f" (median {_safe_fmt(key_results.get('ambiguous_ratio_median'), '.1%')})"
        parts[-1] += "."

    return " ".join(parts)


def generate_methods_pretree_align(
    params: dict[str, Any],
    key_results: dict[str, Any],
    tool_versions: dict[str, Any],
) -> str:
    method = str(params.get("method", "auto")).lower()
    desc, rationale = _ALIGN_METHOD_MAP.get(
        method, (method.upper(), "performs multiple sequence alignment")
    )
    tool_name = "MAFFT" if method != "magus" else "MAGUS"
    version = tool_versions.get("mafft" if method != "magus" else "magus", "unknown version")

    n_aligned = key_results.get("n_aligned", 0)
    n_skipped = key_results.get("n_skipped", 0)
    seq_type = params.get("seq_type", "AA")
    mean_len = key_results.get("mean_alignment_length", 0)

    lines = [
        f"Multiple sequence alignments were performed using {tool_name} v{version} "
        f"with the {desc} algorithm, which {rationale}. "
    ]

    n_line = f"A total of {_describe_n(n_aligned, f'{seq_type} locus', f'{seq_type} loci')} were aligned"
    if n_skipped > 0:
        n_line += f" ({_describe_n(n_skipped, 'locus', 'loci')} skipped). "
    else:
        n_line += ". "
    lines.append(n_line)

    if mean_len:
        unit = "sites" if str(seq_type).upper() == "AA" else "bp"
        lines.append(
            f"Mean alignment length was {_safe_fmt(mean_len, '.1f')} {unit} "
            f"across a mean of {key_results.get('mean_n_taxa', '?')} taxa per locus."
        )

    if params.get("backtrans"):
        trimal_ver = tool_versions.get("trimal", "unknown version")
        lines.append(
            f" Codon-aware nucleotide alignments were produced via back-translation "
            f"using trimAl v{trimal_ver}, preserving reading frame in the nucleotide alignments."
        )

    return "".join(lines)


def generate_methods_pretree_trim(
    params: dict[str, Any],
    key_results: dict[str, Any],
    tool_versions: dict[str, Any],
) -> str:
    tool = params.get("tool", "trimal")
    tool_ver_key = {"trimal": "trimal", "bmge": "bmge", "clipkit": "clipkit"}.get(tool, tool)
    version = tool_versions.get(tool_ver_key, "unknown version")
    seq_type = str(key_results.get("seq_type") or params.get("seq_type", "")).upper()
    unit = "sites" if seq_type == "AA" else "bp"

    n_total = key_results.get("total_genes", 0)
    n_trimmed = key_results.get("trimmed_genes", 0)
    n_skipped = key_results.get("skipped_genes", 0)

    # Determine if backtrans was used (NT trimmed via AA results)
    backtrans = params.get("backtrans") or params.get("nt_dir") or key_results.get("mode") == "aa_then_nt"

    parts = [f"Alignments were trimmed using {str(tool).upper()} v{version}."]
    if backtrans:
        parts.append(
            f" Nucleotide alignments were trimmed based on the corresponding "
            f"amino-acid alignment trimming masks (back-translation mode), "
            f"ensuring codon reading frame is preserved."
        )
    else:
        parts.append(f" All alignments are {seq_type} sequences.")

    parts.append(f" Of {_describe_n(n_total, 'input alignment')}, "
                 f"{_describe_n(n_trimmed, 'was', 'were')} successfully trimmed")
    if n_skipped:
        parts.append(f" ({_describe_n(n_skipped, 'was', 'were')} skipped).")
    else:
        parts.append(".")

    # Length change summary
    before = key_results.get("length_before_mean")
    after = key_results.get("length_after_mean")
    pct = key_results.get("columns_removed_pct_mean")
    if before is not None and after is not None:
        parts.append(
            f" Mean alignment length decreased from {_safe_fmt(before, '.1f')} to {_safe_fmt(after, '.1f')} {unit}"
        )
        if pct is not None:
            parts[-1] += f" ({_safe_fmt(pct, '.1f')}% of columns removed on average)"
        parts[-1] += "."

    return " ".join(parts)


def generate_methods_pretree_metrics(
    params: dict[str, Any],
    key_results: dict[str, Any],
    tool_versions: dict[str, Any],
) -> str:
    n_markers = key_results.get("n_markers", key_results.get("n_files", 0))
    n_metrics = key_results.get("n_metrics", 0)
    tree_dir = params.get("tree_dir")
    seq_type = str(key_results.get("seq_type") or params.get("seq_type", "")).upper()
    text = (
        f"Phylogenetic informativeness metrics were computed using phyloai pretree metrics "
        f"for {_describe_n(n_markers, 'locus', 'loci')}"
    )
    if n_metrics:
        text += f" across {n_metrics} dimensions"
    text += (
        ". Evaluated metrics included alignment-based metrics including "
        "locus length, number of informative sites, gap percentage, "
        "and RCFV (relative composition frequency variability)"
    )
    if tree_dir is not None:
        text += (
            ", as well as gene tree-based metrics including treeness "
            "and average bootstrap support"
        )
    text += (
        ". Pairwise Spearman correlations were computed across all metrics "
        "and visualized as a heatmap for diagnostic evaluation of metric "
        "redundancy and complementarity."
    )
    return text


def generate_methods_pretree_filter_taper(
    params: dict[str, Any],
    key_results: dict[str, Any],
    tool_versions: dict[str, Any],
) -> str:
    taper_ver = tool_versions.get("taper") or tool_versions.get("correction_multi.jl", "unknown version")
    julia_ver = tool_versions.get("julia", "unknown version")
    cutoff = params.get("cutoff", 0.1)
    n_input = key_results.get("n_input", 0)
    n_retained = key_results.get("n_retained", 0)
    n_dropped = key_results.get("n_dropped", 0)
    n_masked = key_results.get("n_masked_sites", 0)

    text = (
        f"Aligned sequences were screened for compositional bias and systematic "
        f"sequencing errors using TAPER v{taper_ver} (correction_multi.jl, "
        f"executed via Julia v{julia_ver}). TAPER applies a moving-window approach "
        f"to identify and mask amino acid sites within individual sequences that "
        f"deviate from expected substitution patterns, without discarding entire loci. "
        f"The masking stringency cutoff was set to {cutoff} (`-c {cutoff}`). "
        f"Of {_describe_n(n_input, 'input locus', 'input loci')}, "
        f"{_describe_n(n_retained, 'was', 'were')} retained"
    )
    if n_dropped > 0:
        text += f" ({_describe_n(n_dropped, 'locus', 'loci')} dropped)"
    if n_masked > 0:
        text += f". A total of {_safe_fmt(n_masked, ',')} sites were masked across all retained loci."
    else:
        text += "."
    return text


def generate_methods_pretree_filter_treeshrink(
    params: dict[str, Any],
    key_results: dict[str, Any],
    tool_versions: dict[str, Any],
) -> str:
    ts_ver = tool_versions.get("run_treeshrink.py") or tool_versions.get("treeshrink", "unknown version")
    threshold = params.get("threshold", 0.05)
    n_input = key_results.get("n_input", 0)
    n_retained = key_results.get("n_retained", 0)
    n_modified = key_results.get("n_modified", 0)
    n_removed = key_results.get("n_removed_taxa_total", 0)
    text = (
        f"Gene trees were screened for outlier long branches using "
        f"TreeShrink v{ts_ver} (α = {threshold}), which removes taxa whose removal "
        f"disproportionately reduces tree diameter. "
        f"Of {_describe_n(n_input, 'input gene tree')}, "
        f"{_describe_n(n_retained, 'was', 'were')} retained unchanged"
    )
    if n_modified > 0:
        text += (
            f" and {_describe_n(n_modified, 'tree', 'trees')} were modified "
            f"by removing a total of {_safe_fmt(n_removed, ',')} outlier taxa."
        )
    else:
        text += "."
    return text


def generate_methods_pretree_filter_symtest(
    params: dict[str, Any],
    key_results: dict[str, Any],
    tool_versions: dict[str, Any],
) -> str:
    iqtree_ver = tool_versions.get("iqtree3") or tool_versions.get("iqtree", "unknown version")
    n_input = key_results.get("n_input", 0)
    n_retained = key_results.get("n_retained", 0)
    n_dropped = key_results.get("n_dropped", 0)
    pval = params.get("symtest_pval", key_results.get("p_value_threshold", 0.05))
    text = (
        f"Alignments were tested for substitutional symmetry using "
        f"IQ-TREE v{iqtree_ver} pairwise symmetry tests (p < {pval}). "
        f"Of {_describe_n(n_input, 'input locus', 'input loci')}, "
        f"{_describe_n(n_retained, 'was', 'were')} retained "
        f"({_describe_n(n_dropped, 'locus', 'loci')} failed the symmetry test)."
    )
    if key_results.get("retained_msa_stats_n_msa"):
        text += (
            f" Retained loci comprised "
            f"{_safe_fmt(key_results.get('retained_msa_stats_total_length'), ',')} "
            f"sites across {key_results.get('retained_msa_stats_n_msa')} loci "
            f"(mean {_safe_fmt(key_results.get('retained_msa_stats_mean_length'), '.1f')} sites)."
        )
    return text


def generate_methods_pretree_filter_metrics(
    params: dict[str, Any],
    key_results: dict[str, Any],
    tool_versions: dict[str, Any],
) -> str:
    n_input = key_results.get("n_total") or key_results.get("n_input", 0)
    n_retained = key_results.get("n_retained", 0)
    n_dropped = key_results.get("n_dropped", 0)
    keep_rule = params.get("keep", "default criteria")
    text = (
        f"Loci were filtered using phyloai pretree filter metrics "
        f"based on phylogenetic informativeness metrics "
        f"(keep rule: {keep_rule}). "
        f"Of {_describe_n(n_input, 'input locus', 'input loci')}, "
        f"{_describe_n(n_retained, 'was', 'were')} retained "
        f"({_describe_n(n_dropped, 'locus', 'loci')} excluded)."
    )
    # Add MSA stats if available
    if key_results.get("retained_msa_stats_n_msa"):
        text += (
            f" Retained loci comprised {_safe_fmt(key_results.get('retained_msa_stats_total_length'), ',')} "
            f"sites across {key_results.get('retained_msa_stats_n_msa')} loci "
            f"(mean {_safe_fmt(key_results.get('retained_msa_stats_mean_length'), '.1f')} sites, "
            f"taxa range {_safe_fmt(key_results.get('retained_msa_stats_mean_taxa'), '.1f')})."
        )
    return text


def generate_methods_pretree_filter_cluster(
    params: dict[str, Any],
    key_results: dict[str, Any],
    tool_versions: dict[str, Any],
) -> str:
    n_input = key_results.get("n_loci") or key_results.get("n_valid_loci") or key_results.get("n_input", 0)
    n_retained = key_results.get("n_retained", 0)
    n_dropped = key_results.get("n_dropped", 0)
    n_clusters = key_results.get("n_clusters", 0)
    n_features = key_results.get("n_features", 0)
    reduction = key_results.get("reduction", params.get("reduction", "umap"))
    drop_clusters = key_results.get("drop_clusters")
    linkage = params.get("cluster_linkage", "ward")
    distance = params.get("cluster_distance", "euclidean")
    outlier_metric = params.get("outlier_metric", "")
    max_drop = params.get("max_drop_fraction", 0.3)

    pkg = "umap-learn" if str(reduction).lower() == "umap" else "scikit-learn"
    text_parts = [
        f"Loci were projected into a low-dimensional space using {str(reduction).upper()} "
        f"(via {pkg}) based on {_describe_n(n_features, 'phylogenetic informativeness', 'phylogenetic informativeness features')}, "
        f"and partitioned into {n_clusters} groups using agglomerative hierarchical clustering "
        f"({linkage} linkage, {distance} distance) with the optimal number of clusters "
        f"determined by silhouette, Davies-Bouldin, and Calinski-Harabasz indices.",
    ]

    # Cluster sizes
    sizes = []
    for ci in range(n_clusters):
        sz = key_results.get(f"cluster_sizes_{ci}")
        if sz is not None:
            sizes.append(f"cluster {ci}: {int(sz)}")
    if sizes:
        text_parts.append(f" Cluster sizes were: {', '.join(sizes)}.")

    # Outlier detection
    if drop_clusters is not None and isinstance(drop_clusters, list) and drop_clusters:
        dropped_ids = [str(d) for d in drop_clusters]
        text_parts.append(
            f" Outlier clusters ({', '.join(dropped_ids)}) were identified "
            f"via Wilcoxon rank-sum tests on {outlier_metric} (max drop fraction: {_safe_fmt(max_drop, '.0%')}) "
            f"and removed, resulting in {_describe_n(n_dropped, 'locus', 'loci')} excluded."
        )
    elif drop_clusters is not None:
        text_parts.append(" No outlier clusters were detected.")

    text_parts.append(
        f" After filtering, {_describe_n(n_retained, 'locus', 'loci')} "
        f"were retained from {_describe_n(n_input, 'input locus', 'input loci')}."
    )

    if key_results.get("retained_msa_stats_n_msa"):
        text_parts.append(
            f" Retained loci comprised {_safe_fmt(key_results.get('retained_msa_stats_total_length'), ',')} "
            f"sites across {key_results.get('retained_msa_stats_n_msa')} loci "
            f"(mean {_safe_fmt(key_results.get('retained_msa_stats_mean_length'), '.1f')} sites, "
            f"mean {_safe_fmt(key_results.get('retained_msa_stats_mean_taxa'), '.1f')} taxa per locus)."
        )

    return " ".join(text_parts)


def generate_methods_pretree_concat(
    params: dict[str, Any],
    key_results: dict[str, Any],
    tool_versions: dict[str, Any],
) -> str:
    taxa_occ = params.get("taxa_occupancy", 0.75)
    seq_type = params.get("seq_type", "AA")
    n_used = key_results.get("n_msa_used", 0)
    n_dropped = key_results.get("n_msa_dropped", 0)
    n_taxa = key_results.get("n_taxa", 0)
    total_len = key_results.get("total_length", 0)
    gap_ratio = key_results.get("gap_ratio", 0)
    pi_ratio = key_results.get("pi_ratio", 0)

    text = (
        f"Trimmed alignments were concatenated into a supermatrix using "
        f"phyloai concat. Loci were included only if they met a minimum "
        f"taxon occupancy threshold of {_safe_fmt(taxa_occ, '.0%')} (`--taxa-occupancy {taxa_occ}`); "
    )
    if n_dropped:
        text += f"{_describe_n(n_dropped, 'locus', 'loci')} were excluded for failing this criterion. "
    text += (
        f"The final supermatrix comprised {_describe_n(n_used, 'locus', 'loci')} "
        f"across {_describe_n(n_taxa, 'taxon', 'taxa')} with a total alignment "
        f"length of {_safe_fmt(total_len, ',')} {seq_type} positions "
        f"(gap ratio: {_safe_fmt(gap_ratio, '.1%')}; parsimony-informative sites: {_safe_fmt(pi_ratio, '.1%')})."
    )

    recoding = params.get("recoding")
    if recoding:
        recoding_groups = {"Dayhoff6": 6, "SR4": 4}.get(recoding, 6)
        text += (
            f" To reduce the influence of substitution saturation and "
            f"compositional heterogeneity, sequences were recoded into "
            f"{recoding} categories (`--recoding {recoding}`), collapsing "
            f"the 20 standard amino acids into {recoding_groups} biochemically "
            f"similar groups; both the original and recoded matrices were retained."
        )

    return text


def generate_methods_pretree_concat_jackknife(
    params: dict[str, Any],
    key_results: dict[str, Any],
    tool_versions: dict[str, Any],
) -> str:
    reps = params.get("replicates", key_results.get("n_replicates", 0))
    target = params.get("target_length", key_results.get("target_length", 0))
    seed = params.get("seed", 42)
    fmt = params.get("to", "fasta")
    mean_len = key_results.get("mean_length")
    min_len = key_results.get("min_length")
    max_len = key_results.get("max_length")
    min_loci = key_results.get("min_loci")
    max_loci = key_results.get("max_loci")
    text = (
        f"A gene-jackknife set of pseudoreplicate matrices was generated from the concatenated "
        f"supermatrix using phyloai pretree concat jackknife. "
        f"A total of {_describe_n(reps, 'pseudoreplicate')} were sampled without replacement "
        f"from partition-defined loci until each reached at least {_safe_fmt(target, ',')} sites "
        f"(`--target-length {target}`), using random seed {seed}. "
        f"Pseudoreplicates were written in {fmt} format."
    )
    if mean_len is not None and min_len is not None and max_len is not None:
        text += (
            f" Final pseudoreplicate lengths ranged from {_safe_fmt(min_len, ',')} "
            f"to {_safe_fmt(max_len, ',')} sites "
            f"(mean {_safe_fmt(mean_len, ',.1f')})"
        )
        if min_loci is not None and max_loci is not None:
            text += (
                f", comprising {min_loci} to {max_loci} loci per replicate"
            )
        text += "."
    return text


def generate_methods_tree_ml_fasttree(
    params: dict[str, Any],
    key_results: dict[str, Any],
    tool_versions: dict[str, Any],
) -> str:
    version = tool_versions.get("FastTree") or tool_versions.get("fasttree", "unknown version")
    n_trees = key_results.get("n_trees", 1)
    model = key_results.get("model") or params.get("model", "LG")
    cat = params.get("cat")
    gamma = params.get("gamma", False)
    boot = params.get("boot")

    text = (
        f"Maximum likelihood phylogenetic trees were inferred using "
        f"FastTree v{version} under the {str(model).upper()} substitution model"
    )
    if cat:
        text += f" with {cat} discrete rate categories (`-cat {cat}`)"
    if gamma:
        text += " and gamma-distributed rates with branch-length rescaling (`-gamma`)"
    text += "."

    if boot:
        text += f" Node support was assessed using SH-like local support with {boot} pseudoreplicates."
    text += f" A total of {_describe_n(n_trees, 'gene tree')} were produced."
    return text


def generate_methods_tree_ml_iqtree(
    params: dict[str, Any],
    key_results: dict[str, Any],
    tool_versions: dict[str, Any],
) -> str:
    version = tool_versions.get("iqtree3") or tool_versions.get("iqtree", "unknown version")
    # partitioned: partitions file is present; merged: rclusterf is set.
    # Also accept legacy boolean keys for backward compatibility.
    partitioned = bool(params.get("partitions")) or params.get("partitioned", False)
    merged = params.get("rclusterf") is not None or params.get("merged_partitions", False)
    rclusterf = params.get("rclusterf", 10) if merged else None
    mfinder = params.get("modelfinder", "MFP")
    mset = params.get("mset")
    boot = params.get("boot")
    alrt = params.get("alrt")
    state_freq = params.get("state_freq", "")
    rate_het = params.get("rate_heterogeneity", "")
    log_lk = key_results.get("log_likelihood")

    text = (
        f"Maximum likelihood phylogenetic inference was performed using "
        f"IQ-TREE v{version}. "
    )

    if partitioned and merged:
        text += (
            f"A partitioned analysis was conducted with partition merging "
            f"enabled (`--merge`), using the rclusterf algorithm "
            f"(`--rclusterf {rclusterf}`) to identify the optimal merging "
            f"scheme by evaluating {rclusterf}% of candidate partition pairs. "
        )
    elif partitioned:
        text += "A partitioned analysis was conducted using the provided partition scheme. "

    # Model description: direct specification vs ModelFinder
    if mfinder and str(mfinder).lower() != "none":
        text += (
            f"Substitution models were selected using ModelFinder ({mfinder}) "
            f"from a candidate set comprising {mset} matrix models (`--mset {mset}`). "
        )
    else:
        # A raw `--tool-args -m` overrides the structured model (e.g. GHOST
        # `-m LG+H4`); describe the actually-executed model, not the structured
        # params spelling the user's run never used.
        raw_m = _tool_arg_value("-m", params.get("tool_args"))
        model_value = raw_m if raw_m is not None else str(params.get("model") or "LG")
        if Path(model_value).is_absolute() or PureWindowsPath(model_value).is_absolute():
            text += "A user-specified custom exchangeability matrix was used. "
            if raw_m is None and rate_het and rate_het != "none":
                text += f"Rate heterogeneity was modeled with {rate_het}. "
            if raw_m is None and params.get("site_freq_file"):
                text += "IQ-TREE site-specific state-frequency profiles (-fs) were used. "
        else:
            model_desc = model_value.upper()
            if raw_m is None:
                if state_freq and state_freq != "none":
                    model_desc += str(state_freq)
                if rate_het and rate_het != "none":
                    model_desc += str(rate_het)
            text += f"The {model_desc} substitution model was used. "

    if log_lk is not None:
        text += f"The final log-likelihood of the best tree was {_safe_fmt(log_lk, '.2f')}. "

    supports = []
    if boot:
        supports.append(f"{_safe_fmt(boot, ',')} ultrafast bootstrap replicates (`-B {boot}`)")
    if alrt:
        supports.append(f"{_safe_fmt(alrt, ',')} SH-aLRT replicates (`--alrt {alrt}`)")
    if supports:
        text += f"Branch support was assessed using {' and '.join(supports)}."
    else:
        text += "Branch support was assessed using ultrafast bootstrap."

    return text


def generate_methods_tree_msc(
    params: dict[str, Any],
    key_results: dict[str, Any],
    tool_versions: dict[str, Any],
) -> str:
    version = tool_versions.get("wastral", "unknown version")
    n_trees = key_results.get("n_input_trees") or key_results.get("n_gene_trees", 0)
    mode = key_results.get("mode", params.get("mode", ""))
    boot = key_results.get("boot", params.get("boot"))
    extra = key_results.get("extra_rounds", params.get("extra_rounds", False))
    boot_type = key_results.get("tree_boot_type", params.get("tree_boot_type", "auto"))

    _MODE_DESC = {1: "hybrid", 2: "branch support weighting", 3: "branch length weighting", 4: "traditional unweighted Astral"}
    _BOOT_DESC = {0: "topology only", 1: "local posterior probability", 2: "quartet + local posterior probability", 3: "quartet + local PP + freqQuad.csv"}
    mode_desc = _MODE_DESC.get(mode, str(mode))
    boot_desc = _BOOT_DESC.get(boot, str(boot))

    text = (
        f"Species tree inference was performed under the multispecies "
        f"coalescent model using wASTRAL v{version}. "
        f"A total of {_describe_n(n_trees, 'gene tree')} were used as input. "
        f"Tree search used {mode_desc} (`--mode {mode}`) "
        f"with {boot_desc} (`--boot {boot}`)"
    )
    if extra:
        text += " and exhaustive search enabled (`-R`)"
    text += "."

    if boot_type and boot_type != "auto":
        _BT_DESC = {"bootstrap": "bootstrap support", "likelihood": "likelihood-based", "abayes": "approximate Bayes"}
        btd = _BT_DESC.get(boot_type, boot_type)
        text += f" Input gene trees contained {btd} values."

    return text


def generate_methods_tree_bi_pb(
    params: dict[str, Any],
    key_results: dict[str, Any],
    tool_versions: dict[str, Any],
) -> str:
    version = tool_versions.get("pb_mpi") or tool_versions.get("phylobayes", "unknown version")
    chains = params.get("chains", 3)
    model = params.get("model", "CAT-GTR")
    mixture = params.get("mixture", "auto")
    gamma_cats = params.get("gamma_cats")
    sample_freq = params.get("sample_freq", 1)
    nsamples = params.get("nsamples", -1)
    threads = params.get("threads")

    text = (
        f"Bayesian phylogenetic inference was performed using "
        f"PhyloBayes-MPI v{version} (pb_mpi) under the "
    )
    # Build model name: CAT-GTR when mixture=auto, or explicit mixture name
    if mixture and str(mixture).lower() == "auto":
        text += f"CAT-{model.upper()} model"
    elif mixture and str(mixture).lower() != "none":
        text += f"{str(mixture).upper()}+{model.upper()} model"
    else:
        text += f"{model.upper()} model"
    if gamma_cats:
        text += f" with {gamma_cats} discrete gamma categories"
    if threads:
        text += f" using {threads} threads"
    text += ". "

    text += (
        f"{chains} independent MCMC chains were run, "
        f"sampling every {sample_freq} generation"
    )
    if nsamples and nsamples != -1:
        text += f" for {nsamples} samples"
    text += "."

    # Chain lengths
    chain_lengths = key_results.get("chain_lengths", {})
    if chain_lengths:
        text += f" Chains reached "
        parts = [f"{v} samples" for v in chain_lengths.values()]
        text += ", ".join(parts) + "."

    # Convergence
    conv = key_results.get("final_convergence", {})
    all_conv = conv.get("all_chains", {})
    pairwise = conv.get("pairwise", {})
    if all_conv:
        status = all_conv.get("status", "unknown")
        maxdiff = all_conv.get("bpcomp_maxdiff")
        effsize = all_conv.get("tracecomp_min_effsize")
        rel_diff = all_conv.get("tracecomp_max_reldiff")
        text += f" Final convergence status: {status}"
        if maxdiff is not None or effsize is not None:
            details = []
            if maxdiff is not None:
                details.append(f"bpcomp maxdiff {_safe_fmt(maxdiff, '.3f')}")
            if effsize is not None:
                details.append(f"min ESS {_safe_fmt(effsize, '.0f')}")
            if rel_diff is not None:
                details.append(f"max rel_diff {_safe_fmt(rel_diff, '.3f')}")
            if details:
                text += f" ({', '.join(details)})"
        text += "."
    if pairwise:
        pw_parts = []
        for pair_name, pair_data in sorted(pairwise.items()):
            s = pair_data.get("status", "?")
            pw_parts.append(f"{pair_name}: {s}")
        if pw_parts:
            text += f" Pairwise: {'; '.join(pw_parts)}."

    return text


def generate_methods_tree_bi_bpcomp(
    params: dict[str, Any],
    key_results: dict[str, Any],
    tool_versions: dict[str, Any],
) -> str:
    version = tool_versions.get("bpcomp", "unknown")
    chains = key_results.get("chains_used", [])
    burnin = params.get("burnin", 0)
    sample_freq = params.get("sample_freq", 1)
    cutoff = params.get("cutoff", 0.5)
    maxdiff = key_results.get("bpcomp_maxdiff")
    meandiff = key_results.get("bpcomp_meandiff")
    status = key_results.get("bpcomp_status", "unknown")
    consensus = key_results.get("consensus_tree", "")

    until = params.get("until", "all")
    until_clause = "" if until == "all" else f" up to sample {until}"

    text = (
        f"Topology convergence was assessed using bpcomp "
        f"(PhyloBayes-MPI v{version}) applied to {_describe_n(len(chains), 'independent MCMC chain')} "
        f"({', '.join(chains)}). "
        f"A burn-in of {_safe_fmt(burnin, ',')} saved samples was discarded; "
        f"trees were sub-sampled every {_safe_fmt(sample_freq, ',')} points{until_clause}. "
        f"The majority-rule consensus cutoff was set to {_safe_fmt(cutoff, '.2f')}. "
    )
    if maxdiff is not None and meandiff is not None:
        text += (
            f"The maximum bipartition frequency discrepancy (maxdiff) between chains "
            f"was {_safe_fmt(maxdiff, '.4f')} and the mean discrepancy (meandiff) was "
            f"{_safe_fmt(meandiff, '.6f')}, indicating {status} convergence. "
        )
    if consensus:
        text += f"The consensus tree was written to {consensus}."
    return text


def generate_methods_tree_bi_tracecomp(
    params: dict[str, Any],
    key_results: dict[str, Any],
    tool_versions: dict[str, Any],
) -> str:
    version = tool_versions.get("tracecomp", "unknown")
    chains = key_results.get("chains_used", [])
    burnin = params.get("burnin", 0)
    min_effsize = key_results.get("tracecomp_min_effsize")
    max_reldiff = key_results.get("tracecomp_max_reldiff")
    status = key_results.get("tracecomp_status", "unknown")

    text = (
        f"Continuous parameter convergence was assessed using tracecomp "
        f"(PhyloBayes-MPI v{version}) applied to {_describe_n(len(chains), 'chain')} "
        f"with a burn-in of {_safe_fmt(burnin, ',')} saved samples. "
    )
    if min_effsize is not None and max_reldiff is not None:
        text += (
            f"The minimum effective sample size across all parameters was "
            f"{_safe_fmt(min_effsize, '.0f')} and the maximum relative difference was "
            f"{_safe_fmt(max_reldiff, '.4f')}, indicating {status} mixing. "
        )
    text += "Per-parameter diagnostics are available in tracecomp.contdiff."
    return text


def generate_methods_tree_bi_readpb(
    params: dict[str, Any],
    key_results: dict[str, Any],
    tool_versions: dict[str, Any],
) -> str:
    version = tool_versions.get("readpb_mpi", "unknown")
    modes = key_results.get("modes_run", [])
    chain = params.get("chain", "")
    burnin = params.get("burnin", 0)
    output_files = key_results.get("output_files", {})
    pp = key_results.get("post_processing", {})

    parts: list[str] = [f"Posterior analysis was performed using readpb_mpi (PhyloBayes-MPI v{version}). "]

    mode_descriptions: dict[str, str] = {
        "rr": (
            f"Posterior mean relative exchangeabilities were estimated from "
            f"post-burnin samples of chain {chain} using readpb_mpi -rr. "
            f"The resulting exchangeability matrix was converted to PAML "
            f"lower-triangle format for use with IQ-TREE."
        ),
        "ss": (
            f"Posterior mean site-specific amino acid frequencies were estimated "
            f"from post-burnin samples using readpb_mpi -ss. "
            f"Site frequency profiles were converted to IQ-TREE -fs format."
        ),
        "r": f"Posterior mean site rates were estimated using readpb_mpi -r.",
        "sitelogl": (
            f"Site-specific marginal log-likelihoods were computed using "
            f"readpb_mpi -sitelogl. These values can be used to compute wAIC "
            f"and leave-one-out cross-validation scores."
        ),
        "ppred": f"Data replicates were simulated from the posterior predictive distribution using readpb_mpi -ppred.",
        "div": f"Posterior predictive diversity test (PPA-DIV) was performed using readpb_mpi -div.",
        "sitecomp": f"Posterior predictive test of compositional heterogeneity across sites (PPA-VAR) was performed using readpb_mpi -sitecomp.",
        "siteconvprob": f"Posterior predictive convergence probability test (PPA-CONV) was performed using readpb_mpi -siteconvprob.",
        "comp": f"Posterior predictive test of compositional homogeneity across taxa was performed using readpb_mpi -comp.",
        "allppred": f"All posterior predictive tests (PPA-DIV, PPA-VAR, PPA-CONV, taxon composition) were performed using readpb_mpi -allppred.",
    }
    for m in modes:
        desc = mode_descriptions.get(m, f"Analysis mode {m} was run using readpb_mpi -{m}.")
        parts.append(f" {desc}")

    partition = pp.get("pmsf_partition", {})
    if partition.get("status") == "success":
        parts.append(
            f" Posterior mean site rates, the selected posterior mean alpha, and "
            f"the discrete Gamma category count were combined with rr-derived "
            f"exchangeabilities and ss-derived site frequencies into "
            f"{partition.get('output', 'partition.PMSF.nex')} for iqtree3 --alisim simulation."
        )

    if pp:
        statuses = []
        for mode_name, info in pp.items():
            s = info.get("status", "error")
            statuses.append(f"{mode_name}: {s}")
        parts.append(f" Post-processing status: {'; '.join(statuses)}.")

    return "".join(parts)


def generate_methods_tree_cf(
    params: dict[str, Any],
    key_results: dict[str, Any],
    tool_versions: dict[str, Any],
) -> str:
    cf_type = str(params.get("cf", "")).lower()
    # Resolved type from key_results (may differ from params, e.g. gcf+scf)
    resolved_type = str(key_results.get("cf_type", cf_type)).lower()

    _CF_NAMES: dict[str, str] = {
        "gcf": "gene concordance factor (gCF)",
        "scf": "site concordance factor (sCF)",
        "scfl": "site concordance factor by likelihood (sCFl)",
        "qcf": "quartet concordance factor (qCF)",
        "gcf+scf": "gene and site concordance factors (gCF/sCF)",
    }
    cf_desc = _CF_NAMES.get(resolved_type, cf_type)

    # qCF uses wASTRAL, others use IQ-TREE
    if "qcf" in resolved_type:
        version = tool_versions.get("wastral", "unknown version")
        tool = "wASTRAL"
    else:
        version = tool_versions.get("iqtree3") or tool_versions.get("iqtree", "unknown version")
        tool = "IQ-TREE"

    return (
        f"{cf_desc} values were calculated using {tool} v{version} "
        f"to assess phylogenetic discordance across the dataset."
    )


def generate_methods_posttree_topology(
    params: dict[str, Any],
    key_results: dict[str, Any],
    tool_versions: dict[str, Any],
) -> str:
    version = tool_versions.get("iqtree3") or tool_versions.get("iqtree", "unknown version")
    n_trees = key_results.get("n_candidate_trees", 0)
    replicates = key_results.get("replicates", params.get("replicates"))
    model = params.get("model_expr") or params.get("model", "")
    best_id = key_results.get("best_tree_id")
    n_rejected = key_results.get("n_rejected_au_0_05", 0)

    text = "Topology hypothesis testing was performed using IQ-TREE"
    if version != "unknown version":
        text += f" v{version}"
    text += "."

    if model:
        text += f" The {model} substitution model was used"
        if replicates:
            text += f" with {_safe_fmt(replicates, ',')} RELL replicates"
        text += "."

    text += (
        f" Approximately unbiased (AU), weighted Kishino-Hasegawa (WKH), and "
        f"weighted Shimodaira-Hasegawa (WSH) tests were applied to "
        f"{_describe_n(n_trees, 'candidate topology', 'candidate topologies')}."
    )

    if best_id:
        text += f" Tree {best_id} was identified as the best topology"
        if n_rejected > 0:
            text += f", and {_describe_n(n_rejected, 'topology', 'topologies')} were significantly rejected (AU test, p < 0.05)"
        text += "."

    return text


def generate_methods_posttree_dating_hessian(
    params: dict[str, Any],
    key_results: dict[str, Any],
    tool_versions: dict[str, Any],
) -> str:
    version = tool_versions.get("iqtree3") or tool_versions.get("iqtree", "unknown version")
    model = params.get("model_expr")
    partitions = params.get("partitions")
    n_partitions = params.get("n_partitions")

    text = (
        f"The Hessian matrix and gradient vectors required for "
        f"approximate likelihood calculation in MCMCtree were computed "
        f"using IQ-TREE v{version} (`--dating mcmctree`)."
    )
    if partitions and n_partitions:
        text += f" A partitioned analysis was performed across {n_partitions} partitions"
        if model:
            text += f" using {model} substitution models"
        else:
            text += " with substitution models automatically selected (AA: LG+F+G4)"
        text += "; partitions were used directly without merging."
    elif model:
        text += f" The {model} substitution model was used."
    else:
        text += " The substitution model was automatically selected based on sequence type."

    return text


def generate_methods_posttree_dating_mcmc(
    params: dict[str, Any],
    key_results: dict[str, Any],
    tool_versions: dict[str, Any],
) -> str:
    paml_ver = tool_versions.get("mcmctree") or tool_versions.get("paml", "unknown version")
    clock = params.get("clock", 2)
    _CLOCK = {1: "strict", 2: "independent-rates", 3: "autocorrelated-rates"}
    clock_desc = _CLOCK.get(clock, str(clock))
    n_runs = params.get("n_runs", key_results.get("n_runs", 2))
    burnin = params.get("burnin", 0)
    nsamples = params.get("nsamples", 0)
    sample_freq = params.get("sample_freq", 1)
    # Guard against non-numeric params in arithmetic
    try:
        total_gen = int(burnin) + int(nsamples) * int(sample_freq)
    except (TypeError, ValueError):
        total_gen = 0

    text = (
        f"Divergence time estimation was performed using MCMCtree "
        f"(PAML v{paml_ver}) under a {clock_desc} molecular clock. "
        f"{n_runs} independent MCMC chains were run for "
        f"{_safe_fmt(total_gen, ',')} generations "
        f"(burn-in: {_safe_fmt(burnin, ',')}, sampling: {_safe_fmt(nsamples, ',')} samples every "
        f"{sample_freq} generations)."
    )

    rho = key_results.get("convergence_rho_posterior")
    n_fail = key_results.get("n_posterior_failures", 0)
    if rho is not None and n_runs > 1:
        text += (
            f" Convergence across independent runs was assessed using "
            f"Spearman rank correlation of posterior node age estimates "
            f"(ρ = {_safe_fmt(rho, '.3f')})."
        )
    if isinstance(n_fail, (int, float)) and n_fail > 0:
        text += f" {_describe_n(n_fail, 'MCMC run', 'MCMC runs')} failed to converge."

    text += (
        f" Posterior node age estimates and 95% highest posterior density "
        f"(HPD) intervals are provided in the output node age tables. "
        f"Diagnostic plots include: trace plots of posterior and prior "
        f"for visual inspection of chain mixing and stationarity; "
        f"convergence scatter plots comparing node ages between "
        f"independent runs; infinite-sites plots to assess "
        f"data sufficiency (posterior: evaluates whether the molecular "
        f"data contain sufficient information for divergence time estimation; "
        f"prior: evaluates whether fossil calibrations provide adequate "
        f"temporal constraints); and posterior-vs-prior comparisons to "
        f"evaluate the influence of the prior on divergence time estimates "
        f"and to detect potential fossil calibration errors."
    )

    return text


def generate_methods_posttree_signal_lnl(
    params: dict[str, Any],
    key_results: dict[str, Any],
    tool_versions: dict[str, Any],
) -> str:
    n_trees = key_results.get("n_trees", 2)
    n_sites = key_results.get("n_sites", "")
    model = params.get("model_expr") or "a partition model"
    iqtree_ver = tool_versions.get("iqtree3", "IQ-TREE3")
    if isinstance(iqtree_ver, str) and not iqtree_ver.lower().startswith("iqtree"):
        iqtree_ver = f"IQ-TREE3 v{iqtree_ver}"
    parts = [
        f"Site-wise log-likelihood scores were computed using {iqtree_ver} "
        f"({model}) for {_describe_n(n_trees, 'candidate topology', 'candidate topologies')} "
        f"across {_describe_n(n_sites, 'alignment site')}, "
        f"following Shen et al. (2017). "
        f"Per-site ΔSLS and per-gene ΔGLS values were derived from these scores by PhyloAI, "
        f"assigning each site (and gene, where boundaries were provided) to the best-supported "
        f"topology."
    ]
    site_counts = key_results.get("site_support_counts", {})
    if site_counts:
        tree_parts = []
        for label, count in sorted(site_counts.items(), key=lambda x: -x[1]):
            if label == "ambiguous":
                tree_parts.append(f"{count} ambiguous")
            else:
                tree_parts.append(f"{count} supported {label}")
        if tree_parts:
            parts.append(f" Among {_describe_n(n_sites, 'site')}, " + ", ".join(tree_parts) + ".")

    gene_counts = key_results.get("gene_support_counts", {})
    n_loci = key_results.get("n_loci")
    if gene_counts:
        gene_tree_parts = []
        for label, count in sorted(gene_counts.items(), key=lambda x: -x[1]):
            if label == "ambiguous":
                gene_tree_parts.append(f"{count} ambiguous")
            else:
                gene_tree_parts.append(f"{count} supported {label}")
        if gene_tree_parts:
            parts.append(
                f" At the locus level, across {_describe_n(n_loci, 'locus', 'loci')}, " + ", ".join(gene_tree_parts) + "."
            )
        n_outlier = key_results.get("n_outlier_genes", 0)
        n_sig_metrics = key_results.get("n_sig_metrics_outlier", 0)
        sig_names = key_results.get("sig_metric_names_outlier", [])
        if isinstance(n_outlier, (int, float)) and n_outlier > 0:
            parts.append(
                f" {_describe_n(n_outlier, 'gene')} with |ΔGLS| exceeding Tukey's 1.5×IQR criterion "
                f"were identified as outlier loci. "
                f"Mann–Whitney U tests comparing outlier vs non-outlier metric profiles "
                f"(outlier_comparison.csv) found {_describe_n(n_sig_metrics, 'metric')}"
            )
            if sig_names:
                parts[-1] += f" significantly different between groups ({', '.join(sig_names)})"
            parts[-1] += " (p < 0.05)."
        elif isinstance(n_outlier, (int, float)) and n_outlier == 0:
            parts.append(" No outlier genes were detected by the |ΔGLS| criterion.")
        support_sig_metrics = key_results.get("support_comparison_sig_metrics", {})
        if support_sig_metrics:
            comparisons = [
                f"{pair} ({', '.join(metrics)})"
                for pair, metrics in support_sig_metrics.items()
            ]
            parts.append(
                " Mann–Whitney U tests comparing metric profiles among genes supporting "
                f"different topologies (support_comparison.csv) found significant differences "
                f"for {', '.join(comparisons)} (p < 0.05)."
            )
        if n_trees == 2:
            n_sig = key_results.get("n_loci_support_sig")
            if isinstance(n_sig, (int, float)):
                parts.append(f" {_describe_n(n_sig, 'locus', 'loci')} had |ΔGLS| ≥ 2 (significant support).")
    elif n_trees == 2 and not gene_counts:
        parts.append(" No locus boundaries were provided; site-wise results only.")
    return " ".join(parts)


def generate_methods_posttree_signal_consistent(
    params: dict[str, Any],
    key_results: dict[str, Any],
    tool_versions: dict[str, Any],
) -> str:
    n_loci = key_results.get("n_loci", "")
    n_consistent = key_results.get("n_consistent", "")
    n_inconsistent = key_results.get("n_inconsistent", "")
    n_skipped = key_results.get("n_gqs_skipped", 0)
    model = params.get("model_expr") or "a partition model"
    iqtree_ver = tool_versions.get("iqtree3", "IQ-TREE3")
    if isinstance(iqtree_ver, str) and not iqtree_ver.lower().startswith("iqtree"):
        iqtree_ver = f"IQ-TREE3 v{iqtree_ver}"
    wastral_ver = tool_versions.get("wastral", "wASTRAL")
    if isinstance(wastral_ver, str) and not wastral_ver.lower().startswith("wastral"):
        wastral_ver = f"wASTRAL v{wastral_ver}"
    n_sig_metrics = key_results.get("n_sig_metrics_consistent", 0)
    sig_names = key_results.get("sig_metric_names_consistent", [])
    parts = [
        f"Consistent genes were identified across {_describe_n(n_loci, 'locus', 'loci')} "
        f"following Shen et al. (2021). Gene-wise log-likelihood scores (GLS) were "
        f"computed using {iqtree_ver} ({model}); gene-wise quartet scores (GQS) were "
        f"computed using {wastral_ver} under the same two candidate topologies. "
        f"{_describe_n(n_consistent, 'locus', 'loci')} were consistent (GLS and GQS agree on the supported topology)."
    ]
    if n_inconsistent:
        parts.append(f" {_describe_n(n_inconsistent, 'locus', 'loci')} were inconsistent (disagree or ambiguous).")
    if n_skipped:
        parts.append(f" {_describe_n(n_skipped, 'locus', 'loci')} were excluded from GQS due to < 4 taxa after pruning.")
    gls_counts = key_results.get("gls_support_counts", {})
    gqs_counts = key_results.get("gqs_support_counts", {})
    if gls_counts:
        tree_parts = []
        for label, count in sorted(gls_counts.items(), key=lambda x: -x[1]):
            if count > 0:
                tree_parts.append(f"{count} {label}")
        if tree_parts:
            parts.append(" By GLS: " + ", ".join(tree_parts) + ".")
    if gqs_counts:
        tree_parts = []
        for label, count in sorted(gqs_counts.items(), key=lambda x: -x[1]):
            if count > 0:
                tree_parts.append(f"{count} {label}")
        if tree_parts:
            parts.append(" By GQS: " + ", ".join(tree_parts) + ".")
    if n_sig_metrics:
        parts.append(
            f" Mann–Whitney U tests comparing consistent vs inconsistent metric profiles "
            f"(consistent_comparison.csv) found {_describe_n(n_sig_metrics, 'metric')}"
        )
        if sig_names:
            parts[-1] += f" significantly different between groups ({', '.join(sig_names)})"
        parts[-1] += " (p < 0.05)."
    return " ".join(parts)


def generate_methods_posttree_signal_fclm(
    params: dict[str, Any],
    key_results: dict[str, Any],
    tool_versions: dict[str, Any],
) -> str:
    n_taxsets = key_results.get("n_taxsets", 4)
    n_taxa = key_results.get("n_taxa", "?")
    n_quartets = key_results.get("n_quartets", "?")
    model = params.get("model_expr") or params.get("partitions") or "a substitution model"
    iqtree_ver = tool_versions.get("iqtree3", "IQ-TREE3")
    if isinstance(iqtree_ver, str) and not iqtree_ver.lower().startswith("iqtree"):
        iqtree_ver = f"IQ-TREE3 v{iqtree_ver}"
    lmap = params.get("lmap")
    if lmap:
        lmap_desc = f"{lmap} quartets"
    elif n_quartets and n_quartets != "?":
        lmap_desc = f"{n_quartets} quartets"
    else:
        lmap_desc = f"50 × {n_taxa} = {n_quartets} quartets"
    try:
        n_taxa_int = int(n_taxa)
    except (TypeError, ValueError):
        n_taxa_int = None
    taxa_span = f" (spanning {_describe_n(n_taxa_int, 'taxon', 'taxa')} total)" if n_taxa_int else ""
    return (
        f"Four-cluster Likelihood Mapping (FcLM; Strimmer & von Haeseler 1997) was performed "
        f"using {iqtree_ver} ({model}) across {_describe_n(n_taxsets, 'taxon cluster', 'taxon clusters')}"
        f"{taxa_span} with {lmap_desc} sampled. "
        f"Cluster definitions were provided via a taxset CSV file "
        f"and converted to NEXUS format for IQ-TREE3 -lmclust. "
        f"Results are reported in the native IQ-TREE .iqtree output, including the relative "
        f"frequencies of quartets supporting each of the three possible unrooted topologies "
        f"(fully resolved, partially resolved, and unresolved regions of the likelihood-mapping triangle). "
        f"The likelihood-mapping triangle is visualised in the .lmap.eps figure."
    )


def generate_methods_posttree_modelcompare_iqtree(
    params: dict[str, Any],
    key_results: dict[str, Any],
    tool_versions: dict[str, Any],
) -> str:
    version = tool_versions.get("iqtree3") or tool_versions.get("iqtree", "unknown version")
    hom_models = params.get("homogeneous_model", "")
    mrate = params.get("mrate", "")
    het_models = params.get("heterogeneous_model")
    het_mrate = params.get("het_mrate", "")
    n_models = key_results.get("n_models_tested", 0)
    best = key_results.get("best_model_bic", "")
    bic = key_results.get("best_model_bic_value")
    w_bic = key_results.get("best_model_bic_weight")

    text = (
        f"Relative model fit was assessed using ModelFinder "
        f"(Kalyaanamoorthy et al. 2017) as implemented in IQ-TREE3 v{version}. "
        f"The homogeneous model search space included {hom_models} "
        f"with rate heterogeneity types {mrate}. "
    )

    if het_models:
        n_madd = key_results.get("n_madd_expanded", 0)
        text += (
            f"Heterogeneous mixture models ({het_models}) were additionally evaluated "
            f"with rate-variation families {het_mrate}, yielding {_describe_n(n_madd, 'expanded model configuration', 'expanded model configurations')} "
            f"passed via -madd. "
        )

    text += f"A total of {_describe_n(n_models, 'model configuration', 'model configurations')} were evaluated. "
    if best:
        text += f"The best-fitting model according to BIC was {best}"
        if bic is not None:
            text += f" (BIC = {_safe_fmt(bic, '.3f')}"
            if w_bic is not None:
                text += f", w-BIC = {_safe_fmt(w_bic, '.4g')}"
            text += ")"
        text += "."

    return text


def generate_methods_posttree_modelcompare_pb(
    params: dict[str, Any],
    key_results: dict[str, Any],
    tool_versions: dict[str, Any],
) -> str:
    n_sites = key_results.get("n_sites", 0)
    n_models = key_results.get("n_models", 1)

    if n_models is not None and n_models > 1:
        best_loocv = key_results.get("best_model_loocv", "")
        best_waic = key_results.get("best_model_waic", "")
        loocv_score = key_results.get("best_loocv_score")
        loocv_quality = key_results.get("best_loocv_quality", "")
        waic_score = key_results.get("best_waic_score")
        waic_quality = key_results.get("best_waic_quality", "")
        text = (
            f"Relative model fit was evaluated using leave-one-out cross-validation (LOO-CV) "
            f"and the widely applicable information criterion (wAIC) following Lartillot (2023), "
            f"computed from site log-likelihood files ({n_sites} sites). "
            f"{_describe_n(n_models, 'candidate model', 'candidate models')} were compared. "
        )
        if best_loocv:
            text += f"The best-fitting model according to LOO-CV was {best_loocv}"
            if loocv_score is not None:
                text += f" (LOO-CV = {_safe_fmt(loocv_score, '.4f')}, Δ = 0"
                if loocv_quality:
                    text += f"; quality: {loocv_quality}"
                text += ")"
            text += ". "
        if best_waic:
            text += f"The best-fitting model according to wAIC was {best_waic}"
            if waic_score is not None:
                text += f" (wAIC = {_safe_fmt(waic_score, '.4f')}"
                if waic_quality:
                    text += f"; quality: {waic_quality}"
                text += ")"
            text += "."
    else:
        n_runs = key_results.get("n_runs", 0)
        loocv_score = key_results.get("best_loocv_score")
        loocv_quality = key_results.get("best_loocv_quality", "")
        loocv_ess = key_results.get("best_loocv_ess")
        loocv_pct = key_results.get("best_loocv_pct_ess_lt10")
        loocv_frac = key_results.get("best_loocv_frac_ess_lt10")
        waic_score = key_results.get("best_waic_score")
        waic_quality = key_results.get("best_waic_quality", "")
        text = (
            f"Model fit was evaluated using leave-one-out cross-validation (LOO-CV) "
            f"and the widely applicable information criterion (wAIC) following Lartillot (2023), "
            f"computed from {_describe_n(n_runs, 'independent MCMC chain', 'independent MCMC chains')} "
            f"site log-likelihood files ({n_sites} sites). "
        )
        if loocv_score is not None:
            text += f"The debiased LOO-CV score was {_safe_fmt(loocv_score, '.4f')}"
            if loocv_quality:
                text += f" (quality: {loocv_quality}"
                if loocv_ess is not None:
                    text += f"; ESS = {_safe_fmt(loocv_ess, '.1f')}"
                if loocv_pct is not None:
                    text += f"; %(ESS<10) = {_safe_fmt(loocv_pct, '.3f')}"
                if loocv_frac is not None:
                    text += f"; f(ESS<10) = {_safe_fmt(loocv_frac, '.3f')}"
                text += ")"
            text += ". "
        if waic_score is not None:
            text += f"The debiased wAIC score was {_safe_fmt(waic_score, '.4f')}"
            if waic_quality:
                text += f" (quality: {waic_quality})"
            text += "."

    return text


_MODE_MEAN_PHRASES = {
    "total": "Mean total branch length per tree",
    "terminal": "Mean terminal branch length",
    "internal": "Mean internal branch length",
    "patristic": "Mean pairwise tip-to-tip distance",
    "tip-to-tip": "Mean tip-to-tip distance",
    "node-to-node": "Mean node-to-node distance",
    "node-to-tip": "Mean node-to-tip distance",
}


def generate_methods_posttree_syserror_brlen(
    params: dict[str, Any],
    key_results: dict[str, Any],
    tool_versions: dict[str, Any],
) -> str:
    n_trees = key_results.get("n_trees", 0)
    modes = key_results.get("modes", [])
    modes_str = ", ".join(modes)
    has_map = bool(params.get("map"))

    text = (
        f"Branch length heterogeneity was assessed by extracting "
        f"{modes_str} branch lengths from "
        f"{_describe_n(n_trees, 'phylogenetic tree', 'phylogenetic trees')}"
    )
    if has_map:
        text += " with internal nodes identified via a monophyletic group map file"
    text += " using PhyloAI (Bio.Phylo)."

    summary = key_results.get("summary", {})
    for mode in modes:
        phrase = _MODE_MEAN_PHRASES.get(mode)
        if phrase is None or mode not in summary:
            continue
        s = summary[mode]
        if s.get("mean") is None:
            continue
        text += (
            f" {phrase} was {_safe_fmt(s.get('mean'), '.4f')}"
            f" (SD = {_safe_fmt(s.get('sd'), '.4f')})."
        )

    n_skipped = key_results.get("n_trees_skipped", 0)
    if n_skipped > 0:
        text += (
            f" {_describe_n(n_skipped, 'tree', 'trees')} were skipped due to "
            f"parsing or monophyly validation failures."
        )

    return text


def generate_methods_posttree_syserror_brlen_label_nodes(
    params: dict[str, Any],
    key_results: dict[str, Any],
    tool_versions: dict[str, Any],
) -> str:
    return ""


def generate_methods_posttree_syserror_rate(
    params: dict[str, Any],
    key_results: dict[str, Any],
    tool_versions: dict[str, Any],
) -> str:
    source = (
        "IQ-TREE empirical-Bayes site-rate estimates"
        if params.get("iqtree_rate") else "PhyloBayes posterior mean site rates"
    )
    text = (
        f"Site-rate heterogeneity was summarized from {source} across "
        f"{_describe_n(key_results.get('n_sites', 0), 'alignment site', 'alignment sites')} using PhyloAI."
    )
    subsets = key_results.get("subsets", [])
    if isinstance(subsets, list):
        details = []
        for item in subsets:
            if not isinstance(item, dict) or not {"subset", "requested_fraction", "selected_sites"} <= item.keys():
                continue
            details.append(
                f"{item['subset']} sites: {_safe_fmt(item['requested_fraction'], '.1%')} "
                f"({item['selected_sites']} sites)"
            )
        if details:
            text += f" Rate-ranked subsets retained {'; '.join(details)}."
    return text


def generate_methods_posttree_syserror_cca(
    params: dict[str, Any],
    key_results: dict[str, Any],
    tool_versions: dict[str, Any],
) -> str:
    return (
        "Cross-comparative analysis (CCA) of systematic error was performed "
        "by comparing tree topologies under different substitution models."
    )


def generate_methods_posttree_syserror_sites(
    params: dict[str, Any],
    key_results: dict[str, Any],
    tool_versions: dict[str, Any],
) -> str:
    return (
        "Site-wise systematic error was diagnosed by evaluating "
        "per-site phylogenetic signal contributions."
    )


def generate_methods_posttree_simulate_alisim_params(
    params: dict[str, Any],
    key_results: dict[str, Any],
    tool_versions: dict[str, Any],
) -> str:
    n_parsed = key_results.get("n_loci_parsed", 0)
    n_matched = key_results.get("n_loci_matched", 0)
    n_unmatched = key_results.get("n_loci_unmatched", 0)
    unmatched_sentence = ""
    if n_unmatched:
        unmatched_sentence = (
            f"{n_unmatched} loci could not be matched to tree files and were "
            "excluded. "
        )
    return (
        f"Simulation parameters were extracted from {n_parsed} IQ-TREE report "
        f"files. {n_matched} loci were successfully matched with corresponding "
        f"tree files. {unmatched_sentence}Extracted parameters include "
        "substitution model, state frequencies, proportion of invariable sites, "
        "rate heterogeneity model, and alignment length for each locus."
    )


def generate_methods_posttree_simulate_alisim_iqtree(
    params: dict[str, Any],
    key_results: dict[str, Any],
    tool_versions: dict[str, Any],
) -> str:
    version = (tool_versions or {}).get("iqtree3", "unknown version")
    if key_results.get("n_msas_generated") is not None:
        seqtype = params.get("seq_type", "?")
        model = params.get("model") or params.get("model_partitions") or "?"
        length = params.get("length", "?")
        n_msas = key_results["n_msas_generated"]
        seed = params.get("seed", "random")
        partition_sentence = ""
        if params.get("model_partitions"):
            partition_sentence = (
                " A per-site partition model was used with edge-proportional "
                "branch lengths (-p)."
            )
        return (
            f"Sequence alignment was simulated using IQ-TREE3 v{version} AliSim "
            f"(Ly-Trong et al. 2023). The simulation used a {seqtype} {model} "
            f"model with {length} sites. {n_msas} replicate alignment(s) were "
            f"generated from the reference tree with random seed {seed}."
            f"{partition_sentence}"
        )
    source_loci = key_results.get("source_loci", 0)
    strategy = params.get("strategy", "?")
    n_completed = key_results.get("n_simulations_completed", 0)
    n_failed = key_results.get("n_simulations_failed", 0)

    strategy_sentence = {
        "complete": (
            "the complete sampling strategy, in which each simulated alignment "
            "replicates the full parameter set of a single source gene model "
            "(substitution model, rate heterogeneity, alignment length, "
            "invariant-site proportion, and reference tree all taken together "
            "from one empirical row of the source table)"
        ),
        "mixed": (
            "the mixed sampling strategy, in which the model core, rate "
            "heterogeneity group, alignment length, invariant-site proportion, "
            "and reference tree were each sampled independently from the "
            "empirical gene-model distribution, preserving the empirical "
            "distributions of individual parameters and their "
            "presence/absence ratios"
        ),
        "pdf": (
            "the probability density function (PDF) sampling strategy, built "
            "on mixed sampling; for the parameters "
            f"{params.get('pdf_params', '')}, values were resampled from "
            "histogram-based estimates of the empirical probability density "
            f"(Freedman-Diaconis binning) with noise scale "
            f"{params.get('noise_scale', 1.0)}"
        ),
    }.get(strategy, f"the {strategy} sampling strategy")

    override_sentence = ""
    if params.get("override"):
        override_sentence = (
            f"The following parameters were fixed across all simulations: "
            f"{params['override']}. "
        )
    tree_sentence = (
        "Reference trees were sampled together with all model parameters "
        "(complete strategy)."
        if strategy == "complete"
        else "Reference trees were sampled independently from the model parameters."
    )
    return (
        f"Sequence alignments were simulated using IQ-TREE3 v{version} AliSim "
        f"(Ly-Trong et al. 2023). Simulation parameters were drawn from an "
        f"empirical distribution of {source_loci} gene models using "
        f"{strategy_sentence}. {override_sentence}A total of {n_completed} "
        f"alignments were generated ({n_failed} failed). {tree_sentence}"
    )


def generate_methods_posttree_simulate_alisim_transfergaps(
    params: dict[str, Any],
    key_results: dict[str, Any],
    tool_versions: dict[str, Any],
) -> str:
    detected = key_results.get("detected_seq_type", params.get("seq_type", "AA"))
    valid = "ACDEFGHIKLMNPQRSTVWY" if detected == "AA" else "ACGT"
    n_msas = key_results.get("n_msas", 1)
    target = "alignments" if n_msas > 1 else "alignment"
    count_sentence = (
        f"Gap patterns from the original alignment were transferred to {n_msas} "
        f"simulated {target} to restore biologically realistic indel patterns. "
    ) if n_msas > 1 else (
        "Gap patterns from the original alignment were transferred to the "
        "simulated alignment to restore biologically realistic indel patterns. "
    )
    exclude_ambiguity = bool(params.get("exclude_ambiguity", False))
    if exclude_ambiguity:
        masking_sentence = (
            "Positions marked as gaps (-) in the original sequences were replaced "
            "with gap characters (-) at the corresponding positions in the "
            "simulated sequences; ambiguity codes were left untouched."
        )
    else:
        masking_sentence = (
            "Positions containing non-standard characters (gaps and ambiguity "
            "codes) in the original sequences were replaced with gap characters "
            f"(-) at the corresponding positions in the simulated sequences. The "
            f"valid character set was {valid}; all other characters were treated "
            "as gaps."
        )
    return count_sentence + masking_sentence


def generate_methods_posttree_simulate_adequacy(
    params: dict[str, Any],
    key_results: dict[str, Any],
    tool_versions: dict[str, Any],
) -> str:
    return (
        "Model adequacy was assessed with phyloai posttree simulate adequacy by comparing four summary "
        f"statistics from the observed {key_results.get('n_taxa', '?')}-taxon "
        f"{key_results.get('seq_type', params.get('seq_type', '?'))} alignment "
        f"({key_results.get('n_sites', '?')} sites) against "
        f"{key_results.get('n_simulations', 0)} simulated replicates. Mean diversity "
        "per site (PPA-DIV), mean squared empirical state frequency (PPA-CONV), "
        "mean variance of site-specific frequencies (PPA-VAR), and maximum/mean "
        "squared compositional deviation across taxa (PPA-COMP) were calculated. "
        "For each statistic, the null distribution was summarized using its mean, "
        "population SD, and empirical 95% interval (p2.5-p97.5); observed values "
        "were assessed using z-scores and posterior predictive p-values. Values with "
        "|z| > 2 or pp < 0.05 were treated as potential model inadequacy."
    )


# ---------------------------------------------------------------------------
# Methods generator registry
# ---------------------------------------------------------------------------

METHODS_GENERATORS: dict[str, Any] = {
    "pretree.convert": generate_methods_pretree_convert,
    "pretree.stats": generate_methods_pretree_stats,
    "pretree.align": generate_methods_pretree_align,
    "pretree.trim": generate_methods_pretree_trim,
    "pretree.metrics": generate_methods_pretree_metrics,
    "pretree.filter.taper": generate_methods_pretree_filter_taper,
    "pretree.filter.treeshrink": generate_methods_pretree_filter_treeshrink,
    "pretree.filter.symtest": generate_methods_pretree_filter_symtest,
    "pretree.filter.metrics": generate_methods_pretree_filter_metrics,
    "pretree.filter.cluster": generate_methods_pretree_filter_cluster,
    "pretree.concat": generate_methods_pretree_concat,
    "pretree.concat.jackknife": generate_methods_pretree_concat_jackknife,
    "tree.ml.fasttree": generate_methods_tree_ml_fasttree,
    "tree.ml.iqtree": generate_methods_tree_ml_iqtree,
    "tree.msc": generate_methods_tree_msc,
    "tree.bi.pb": generate_methods_tree_bi_pb,
    "tree.bi.bpcomp": generate_methods_tree_bi_bpcomp,
    "tree.bi.tracecomp": generate_methods_tree_bi_tracecomp,
    "tree.bi.readpb": generate_methods_tree_bi_readpb,
    "tree.cf": generate_methods_tree_cf,
    "posttree.topology": generate_methods_posttree_topology,
    "posttree.dating.hessian": generate_methods_posttree_dating_hessian,
    "posttree.dating.mcmc": generate_methods_posttree_dating_mcmc,
    "posttree.signal.lnl": generate_methods_posttree_signal_lnl,
    "posttree.signal.consistent": generate_methods_posttree_signal_consistent,
    "posttree.signal.fclm": generate_methods_posttree_signal_fclm,
    "posttree.modelcompare.iqtree": generate_methods_posttree_modelcompare_iqtree,
    "posttree.modelcompare.pb": generate_methods_posttree_modelcompare_pb,
    "posttree.syserror.brlen": generate_methods_posttree_syserror_brlen,
    "posttree.syserror.brlen.label-nodes": generate_methods_posttree_syserror_brlen_label_nodes,
    "posttree.syserror.rate": generate_methods_posttree_syserror_rate,
    "posttree.syserror.cca": generate_methods_posttree_syserror_cca,
    "posttree.syserror.sites": generate_methods_posttree_syserror_sites,
    "posttree.simulate.alisim.params": generate_methods_posttree_simulate_alisim_params,
    "posttree.simulate.alisim.iqtree": generate_methods_posttree_simulate_alisim_iqtree,
    "posttree.simulate.alisim.transfergaps": generate_methods_posttree_simulate_alisim_transfergaps,
    "posttree.simulate.adequacy": generate_methods_posttree_simulate_adequacy,
}


def generate_all_methods(
    step_id: str,
    params: dict[str, Any],
    key_results: dict[str, Any],
    tool_versions: dict[str, Any],
    status: str = "success",
) -> str:
    """Dispatch to the appropriate methods generator for step_id.

    Returns empty string for unknown step_ids or failed steps.
    """
    if status != "success":
        return ""
    generator = METHODS_GENERATORS.get(step_id)
    if generator is None:
        return ""
    return generator(params, key_results, tool_versions)
