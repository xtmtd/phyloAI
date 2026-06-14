"""Concatenate multiple MSAs into a supermatrix for phylogenetic inference."""

from __future__ import annotations

import json
import shutil
import time
from pathlib import Path
from typing import Any

from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from rich.panel import Panel
from rich.table import Table

from phyloai.core.formats import AlignmentFormat, FormatConverter
from phyloai.core.schema import COMMON_ALIGNMENT_EXTENSIONS
from phyloai.core.sequence_normalization import detect_seq_type, normalize_sequences

AA_RECODING_TABLES: dict[str, dict[str, str]] = {
    "Dayhoff-6": {
        "A": "0", "G": "0", "P": "0", "S": "0", "T": "0",
        "D": "1", "E": "1", "N": "1", "Q": "1",
        "H": "2", "K": "2", "R": "2",
        "I": "3", "L": "3", "M": "3", "V": "3",
        "F": "4", "W": "4", "Y": "4",
        "C": "5",
        "X": "X",
    },
    "Dayhoff-9": {
        "D": "0", "E": "0", "H": "0", "N": "0", "Q": "0",
        "I": "1", "L": "1", "M": "1", "V": "1",
        "F": "2", "Y": "2",
        "A": "3", "S": "3", "T": "3",
        "K": "4", "R": "4",
        "G": "5",
        "P": "6",
        "C": "7",
        "W": "8",
        "X": "X",
    },
    "Dayhoff-12": {
        "D": "0", "E": "0", "Q": "0",
        "M": "1", "L": "1", "I": "1", "V": "1",
        "F": "2", "Y": "2",
        "K": "3", "H": "3", "R": "3",
        "G": "4",
        "A": "5",
        "P": "6",
        "S": "7",
        "T": "8",
        "N": "9",
        "W": "A",
        "C": "B",
        "X": "X",
    },
    "Dayhoff-15": {
        "D": "0", "E": "0", "Q": "0",
        "M": "1", "L": "1",
        "I": "2", "V": "2",
        "F": "3", "Y": "3",
        "G": "4",
        "A": "5",
        "P": "6",
        "S": "7",
        "T": "8",
        "N": "9",
        "K": "A",
        "H": "B",
        "R": "C",
        "W": "D",
        "C": "E",
        "X": "X",
    },
    "Dayhoff-18": {
        "F": "0", "Y": "0",
        "M": "1", "L": "1",
        "I": "2",
        "V": "3",
        "G": "4",
        "A": "5",
        "P": "6",
        "S": "7",
        "T": "8",
        "D": "9",
        "E": "A",
        "Q": "B",
        "N": "C",
        "K": "D",
        "H": "E",
        "R": "F",
        "W": "G",
        "C": "H",
        "X": "X",
    },
    "SandR-6": {
        "A": "0", "P": "0", "S": "0", "T": "0",
        "D": "1", "E": "1", "N": "1", "G": "1",
        "Q": "2", "K": "2", "R": "2",
        "M": "3", "I": "3", "V": "3", "L": "3",
        "W": "4", "C": "4",
        "F": "5", "Y": "5", "H": "5",
        "X": "X",
    },
    "KGB-6": {
        "A": "0", "G": "0", "P": "0", "S": "0",
        "D": "1", "E": "1", "N": "1", "Q": "1", "H": "1", "K": "1", "R": "1", "T": "1",
        "M": "2", "I": "2", "L": "2",
        "W": "3",
        "F": "4", "Y": "4",
        "C": "5", "V": "5",
        "X": "X",
    },
}

NT_RECODING_TABLES: dict[str, dict[str, str]] = {
    "RY-nucleotide": {
        "A": "R", "G": "R",
        "C": "Y", "T": "Y", "U": "Y",
        "N": "?",
        "X": "?",
        "R": "R", "Y": "Y",
        "S": "?", "W": "?", "K": "?", "M": "?",
        "B": "?", "D": "?", "H": "?", "V": "?",
        "-": "-",
        "?": "?",
        ".": ".",
    },
}

_GAP_CHARS = frozenset("-?.*")


def _apply_recoding(
    matrix: dict[str, str], scheme: str
) -> tuple[dict[str, str], list[str]]:
    if scheme in AA_RECODING_TABLES:
        table = AA_RECODING_TABLES[scheme]
    elif scheme in NT_RECODING_TABLES:
        table = NT_RECODING_TABLES[scheme]
    else:
        raise ValueError(f"Unknown recoding scheme: {scheme!r}")

    warnings_set: set[str] = set()
    result: dict[str, str] = {}
    for taxon, seq in matrix.items():
        chars: list[str] = []
        for ch in seq:
            if ch in _GAP_CHARS:
                chars.append(ch)
            elif ch in table:
                chars.append(table[ch])
            else:
                chars.append(ch)
                warnings_set.add(
                    f"Character '{ch}' not in recoding table '{scheme}', passed through unchanged"
                )
        result[taxon] = "".join(chars)
    return result, sorted(warnings_set)


def _translate_codon(seq: str) -> str:
    n_complete = (len(seq) // 3) * 3
    result: list[str] = []
    for i in range(0, n_complete, 3):
        codon = seq[i:i + 3]
        if "-" in codon:
            result.append("-")
        else:
            result.append(str(Seq(codon).translate()))
    return "".join(result)


def _exclude_codon3(seq: str) -> str:
    return "".join(ch for i, ch in enumerate(seq) if i % 3 != 2)


def _scan_msa_files(msa_dir: Path) -> list[Path]:
    found = []
    for ext in COMMON_ALIGNMENT_EXTENSIONS:
        found.extend(path for path in msa_dir.glob(f"*{ext}") if path.is_file())
    return sorted(set(found))


def _read_msa(path: Path) -> tuple[list[str], list[str], int]:
    converter = FormatConverter()
    fmt = converter.detect(path)
    alignment = converter.read(path, source_format=fmt)
    taxa = [record.id for record in alignment]
    seqs = [str(record.seq).upper() for record in alignment]
    length = alignment.get_alignment_length()
    return taxa, seqs, length


def _read_msa_headers(path: Path) -> list[str]:
    suffix = path.suffix.lower()
    if suffix in (".fa", ".faa", ".ffn", ".frn", ".fasta", ".fas"):
        ids = []
        with open(path) as fh:
            for line in fh:
                line = line.strip()
                if line.startswith(">"):
                    ids.append(line[1:].split(None, 1)[0])
        return ids
    else:
        converter = FormatConverter()
        fmt = converter.detect(path)
        alignment = converter.read(path, source_format=fmt)
        return [record.id for record in alignment]


def _filter_by_occupancy(
    msa_paths: list[Path],
    msa_taxa: dict[str, set[str]],
    total_taxa: set[str],
    threshold: float,
) -> tuple[list[Path], list[dict[str, Any]]]:
    n_total = len(total_taxa)
    kept: list[Path] = []
    dropped: list[dict[str, Any]] = []
    for path in msa_paths:
        taxa = msa_taxa[str(path)]
        n_taxa = len(taxa)
        ratio = n_taxa / n_total if n_total > 0 else 0.0
        if ratio >= threshold:
            kept.append(path)
        else:
            dropped.append({
                "filename": path.name,
                "n_taxa": n_taxa,
                "occupancy_ratio": round(ratio, 4),
                "total_taxa": n_total,
            })
    return kept, dropped


def _concat_alignments(
    msa_paths: list[Path],
    msa_data: dict[str, tuple[list[str], list[str], int]],
    total_taxa: set[str],
) -> tuple[dict[str, str], list[str]]:
    matrix: dict[str, list[str]] = {taxon: [] for taxon in total_taxa}
    for path in msa_paths:
        taxa, seqs, length = msa_data[str(path)]
        taxon_to_seq = dict(zip(taxa, seqs))
        for taxon in total_taxa:
            seq = taxon_to_seq.get(taxon, "?" * length)
            matrix[taxon].append(seq)
    concatenated = {taxon: "".join(parts) for taxon, parts in matrix.items()}
    taxon_order = sorted(total_taxa)
    return concatenated, taxon_order


def _reorder_outgroup(matrix: dict[str, str], outgroup: str | None) -> dict[str, str]:
    if outgroup is None:
        return matrix
    if outgroup not in matrix:
        raise ValueError(f"Outgroup taxon {outgroup!r} not found in matrix")
    reordered = {outgroup: matrix[outgroup]}
    for taxon, seq in matrix.items():
        if taxon != outgroup:
            reordered[taxon] = seq
    return reordered


def _write_matrix(
    matrix: dict[str, str],
    out_path: Path,
    target_format: str,
    seq_type: str,
) -> list[dict[str, str]]:
    fmt = AlignmentFormat(target_format)
    records = [SeqRecord(Seq(seq), id=taxon, description="") for taxon, seq in matrix.items()]
    alignment = MultipleSeqAlignment(records)
    molecule_type = "protein" if seq_type == "AA" else "DNA"
    converter = FormatConverter()
    return converter.write_alignment(alignment, out_path, fmt, molecule_type=molecule_type)


def _compute_concat_stats(matrix: dict[str, str], seq_type: str) -> dict[str, Any]:
    from phyloai.pretree.stats import compute_site_patterns, _summarize_per_taxon, per_taxon_stats

    records = [SeqRecord(Seq(seq), id=taxon) for taxon, seq in matrix.items()]
    sequences = [str(record.seq) for record in records]
    n_taxa = len(sequences)
    alignment_length = len(sequences[0]) if sequences else 0

    if seq_type == "other":
        gap_set = {"-", "?"}
        total_chars = sum(len(s) for s in sequences)
        total_gap = sum(sum(1 for ch in s if ch in gap_set) for s in sequences)
        gap_ratio = round(total_gap / total_chars, 4) if total_chars > 0 else 0.0

        per_taxon = []
        for record in records:
            seq = str(record.seq)
            raw = len(seq)
            gc = sum(1 for ch in seq if ch in gap_set)
            per_taxon.append({
                "name": record.id, "raw_length": raw, "ungapped_length": raw - gc,
                "gap_ratio": round(gc / raw, 4) if raw > 0 else 0.0,
                "ambiguous_ratio": 0.0,
                "standard_ratio": round(1 - gc / raw, 4) if raw > 0 else 1.0,
            })

        character_summary = {
            "gap_ratio": gap_ratio, "ambiguous_ratio": 0.0,
            "gap_ambiguous_ratio": gap_ratio,
            "standard_ratio": round(1 - gap_ratio, 4),
        }

        gap_codes = {ord("-"), ord("?")}
        distinct: set[bytes] = set()
        pi = 0
        ss = 0
        for col in zip(*(s.encode("ascii") for s in sequences)):
            distinct.add(bytes(ord("-") if ch in gap_codes else ch for ch in col))
            counts: dict[int, int] = {}
            for ch in col:
                if ch not in gap_codes:
                    counts[ch] = counts.get(ch, 0) + 1
            if sum(counts.values()) < 2:
                continue
            if len(counts) == 1:
                continue
            repeated = sum(1 for v in counts.values() if v >= 2)
            if repeated >= 2:
                pi += 1
            else:
                ss += 1
        var = pi + ss
        const = alignment_length - var

        site_patterns = {
            "alignment_length": alignment_length,
            "distinct_patterns": {"count": len(distinct), "ratio": round(len(distinct) / alignment_length, 4) if alignment_length > 0 else 0.0},
            "constant_sites": {"count": const, "ratio": round(const / alignment_length, 4) if alignment_length > 0 else 0.0},
            "parsimony_informative": {"count": pi, "ratio": round(pi / alignment_length, 4) if alignment_length > 0 else 0.0},
            "singleton_sites": {"count": ss, "ratio": round(ss / alignment_length, 4) if alignment_length > 0 else 0.0},
        }

        return {
            "n_taxa": n_taxa, "alignment_length": alignment_length,
            "n_msa_used": 0, "seq_type": seq_type,
            "character_summary": character_summary,
            "site_patterns": site_patterns, "per_taxon": per_taxon,
        }

    stats_seq_type = "NT" if seq_type == "CODON" else seq_type

    per_taxon = [per_taxon_stats(record, stats_seq_type) for record in records]
    summary = _summarize_per_taxon(per_taxon)
    site_patterns = compute_site_patterns(sequences, stats_seq_type)

    return {
        "n_taxa": n_taxa,
        "alignment_length": alignment_length,
        "n_msa_used": 0,
        "seq_type": seq_type,
        "character_summary": summary,
        "site_patterns": site_patterns,
        "per_taxon": per_taxon,
    }


def _render_concat_panels(overview: dict[str, Any], variant_stats: list[dict[str, Any]]) -> list[Panel]:
    overview_table = Table(show_header=False)
    overview_table.add_column("Metric")
    overview_table.add_column("Value")
    overview_table.add_row("prefix", str(overview.get("prefix", "")))
    overview_table.add_row("to_format", str(overview.get("to_format", "")))
    overview_table.add_row("n_taxa", str(overview.get("n_taxa", "")))
    overview_table.add_row("n_msa_input", str(overview.get("n_msa_input", "")))
    overview_table.add_row("n_msa_used", str(overview.get("n_msa_used", "")))
    overview_table.add_row("n_msa_dropped", str(overview.get("n_msa_dropped", "")))
    overview_table.add_row("taxon_occupancy_threshold", str(overview.get("taxon_occupancy_threshold", "")))
    if overview.get("recoding"):
        overview_table.add_row("recoding", str(overview["recoding"]))
    if overview.get("outgroup"):
        overview_table.add_row("outgroup", str(overview["outgroup"]))
    overview_table.add_row("variants_produced", str(overview.get("variants_produced", "")))

    variant_names = [v["variant"] for v in variant_stats]

    char_table = Table(title="Character Summary")
    char_table.add_column("")
    for name in variant_names:
        char_table.add_column(name)
    char_metrics = ["seq_type", "total_length", "gap_ratio", "ambiguous_ratio",
                    "gap_ambiguous_ratio", "standard_ratio"]
    for metric in char_metrics:
        row = [metric]
        for v in variant_stats:
            if metric in ("seq_type", "total_length"):
                row.append(str(v.get(metric, "")))
            else:
                cs = v.get("character_summary")
                if cs is None:
                    row.append("\u2014")
                else:
                    val = cs.get(metric, "")
                    row.append(f"{val:.4f}" if isinstance(val, float) else str(val))
        char_table.add_row(*row)

    site_table = Table(title="Site Patterns")
    site_table.add_column("")
    for name in variant_names:
        site_table.add_column(name)
    site_metrics = ["alignment_length", "distinct_patterns", "constant_sites",
                    "parsimony_informative", "singleton_sites"]
    for metric in site_metrics:
        row = [metric]
        for v in variant_stats:
            sp = v.get("site_patterns")
            if sp is None:
                row.append("\u2014")
            elif metric == "alignment_length":
                row.append(str(sp.get(metric, "")))
            else:
                entry = sp.get(metric, {})
                if isinstance(entry, dict):
                    row.append(f"{entry.get('count', 0)} ({entry.get('ratio', 0):.4f})")
                else:
                    row.append(str(entry))
        site_table.add_row(*row)

    return [
        Panel(overview_table, title="Overview"),
        Panel(char_table, title="Character Summary"),
        Panel(site_table, title="Site Patterns"),
    ]


def run_concat(
    msa_dir: Path,
    output_dir: Path,
    prefix: str = "matrix",
    seq_type: str = "auto",
    taxa_occupancy: float = 0.5,
    recoding: str | None = None,
    outgroup: str | None = None,
    to_format: str = "fasta",
    translate_codon: bool = False,
    exclude_codon3: bool = False,
    dry_run: bool = False,
    overwrite: bool = False,
) -> dict[str, Any]:
    start_time = time.time()

    if not msa_dir.exists():
        raise ValueError(f"MSA directory '{msa_dir}' does not exist")
    if not msa_dir.is_dir():
        raise ValueError(f"'{msa_dir}' is not a directory")

    if output_dir.exists() and any(output_dir.iterdir()):
        if not overwrite:
            raise ValueError(
                f"Output directory '{output_dir}' is non-empty. Use --overwrite to replace."
            )
        if not dry_run:
            shutil.rmtree(output_dir)
    if not dry_run:
        output_dir.mkdir(parents=True, exist_ok=True)

    msa_paths = _scan_msa_files(msa_dir)
    if not msa_paths:
        raise ValueError(f"No alignment files found in '{msa_dir}'")

    # --- Pass 1: Header-first scan ---------------------------------------------
    all_taxa: set[str] = set()
    msa_taxa_map: dict[str, set[str]] = {}
    for path in msa_paths:
        taxa = _read_msa_headers(path)
        all_taxa.update(taxa)
        msa_taxa_map[str(path)] = set(taxa)

    # --- Auto-detect seq_type (sample 3 files with full read) -------------------
    if seq_type == "auto":
        sample_seqs: list[str] = []
        for path in msa_paths[:3]:
            _, seqs, _ = _read_msa(path)
            sample_seqs.extend(seqs[:10])
        resolved_seq_type = detect_seq_type(sample_seqs)
    else:
        resolved_seq_type = seq_type

    # --- Validation -------------------------------------------------------------
    if resolved_seq_type != "CODON" and (translate_codon or exclude_codon3):
        raise ValueError(
            "--translate-codon and --exclude-codon3 require --seq-type CODON, "
            f"got {resolved_seq_type}"
        )
    if recoding:
        if recoding in AA_RECODING_TABLES and resolved_seq_type not in ("AA",):
            raise ValueError(
                f"Recoding scheme '{recoding}' requires AA seq_type, got {resolved_seq_type}"
            )
        if recoding in NT_RECODING_TABLES and resolved_seq_type not in ("NT", "CODON"):
            raise ValueError(
                f"Recoding scheme '{recoding}' requires NT or CODON seq_type, got {resolved_seq_type}"
            )

    # --- Occupancy filtering ----------------------------------------------------
    kept_paths, dropped = _filter_by_occupancy(msa_paths, msa_taxa_map, all_taxa, taxa_occupancy)
    if not kept_paths:
        raise ValueError("No MSAs passed occupancy filtering")

    # --- Pass 2: Streaming concat -----------------------------------------------
    def _accumulate_replacements(counts: dict[str, int]) -> None:
        for key, value in counts.items():
            all_normalization_replacements[key] = all_normalization_replacements.get(key, 0) + value

    all_normalization_replacements: dict[str, int] = {}
    norm_seq_type = "NT" if resolved_seq_type == "CODON" else resolved_seq_type
    needs_variant_data = resolved_seq_type == "CODON" and (translate_codon or exclude_codon3)
    msa_data: dict[str, tuple[list[str], list[str], int]] = {}

    matrix_parts: dict[str, list[str]] = {taxon: [] for taxon in all_taxa}
    for path in kept_paths:
        taxa, seqs, length = _read_msa(path)
        norm = normalize_sequences(seqs, norm_seq_type)
        _accumulate_replacements(norm.replacements)
        normalized_seqs = norm.sequences

        if needs_variant_data:
            msa_data[str(path)] = (taxa, normalized_seqs, length)

        taxon_to_seq = dict(zip(taxa, normalized_seqs))
        for taxon in all_taxa:
            seq = taxon_to_seq.get(taxon, "?" * length)
            matrix_parts[taxon].append(seq)

    matrix = {taxon: "".join(parts) for taxon, parts in matrix_parts.items()}

    # --- Variant generation -----------------------------------------------------
    variants: list[dict[str, Any]] = []
    ext_map = {"fasta": ".fa", "phylip-relaxed": ".phy", "phylip-paml": ".phy", "nexus": ".nex"}
    ext = ext_map.get(to_format, ".fa")
    recoding_warnings: list[str] = []

    matrix = _reorder_outgroup(matrix, outgroup)
    original_path = output_dir / f"{prefix}{ext}"
    if not dry_run:
        _write_matrix(matrix, original_path, to_format, resolved_seq_type)
    variants.append({
        "variant": "original", "path": str(original_path),
        "seq_type": resolved_seq_type,
        "length": len(list(matrix.values())[0]) if matrix else 0,
    })

    if recoding:
        recoded_seq_type = "other"
        recoded_matrix, rw = _apply_recoding(matrix, recoding)
        recoding_warnings = rw
        recoded_matrix = _reorder_outgroup(recoded_matrix, outgroup)
        recoded_path = output_dir / f"{prefix}.recoded{ext}"
        if not dry_run:
            _write_matrix(recoded_matrix, recoded_path, to_format, recoded_seq_type)
        variants.append({
            "variant": "recoded", "path": str(recoded_path),
            "seq_type": recoded_seq_type,
            "length": len(list(recoded_matrix.values())[0]) if recoded_matrix else 0,
        })

    if resolved_seq_type == "CODON" and translate_codon:
        translated_data: dict[str, tuple[list[str], list[str], int]] = {}
        for path in kept_paths:
            taxa, seqs, _ = msa_data[str(path)]
            translated_seqs = [_translate_codon(seq) for seq in seqs]
            translated_len = len(translated_seqs[0]) if translated_seqs else 0
            translated_data[str(path)] = (taxa, translated_seqs, translated_len)
        translated_matrix, _ = _concat_alignments(kept_paths, translated_data, all_taxa)
        translated_taxa = list(translated_matrix.keys())
        tnorm = normalize_sequences([translated_matrix[t] for t in translated_taxa], "AA")
        translated_matrix = dict(zip(translated_taxa, tnorm.sequences))
        _accumulate_replacements(tnorm.replacements)
        translated_matrix = _reorder_outgroup(translated_matrix, outgroup)
        translated_path = output_dir / f"{prefix}.translated{ext}"
        if not dry_run:
            _write_matrix(translated_matrix, translated_path, to_format, "AA")
        variants.append({
            "variant": "translated", "path": str(translated_path),
            "seq_type": "AA",
            "length": len(list(translated_matrix.values())[0]) if translated_matrix else 0,
        })

    if resolved_seq_type == "CODON" and exclude_codon3:
        cds12_data: dict[str, tuple[list[str], list[str], int]] = {}
        for path in kept_paths:
            taxa, seqs, _ = msa_data[str(path)]
            cds12_seqs = [_exclude_codon3(seq) for seq in seqs]
            cds12_len = len(cds12_seqs[0]) if cds12_seqs else 0
            cds12_data[str(path)] = (taxa, cds12_seqs, cds12_len)
        cds12_matrix, _ = _concat_alignments(kept_paths, cds12_data, all_taxa)
        cds12_taxa = list(cds12_matrix.keys())
        cnorm = normalize_sequences([cds12_matrix[t] for t in cds12_taxa], "NT")
        cds12_matrix = dict(zip(cds12_taxa, cnorm.sequences))
        _accumulate_replacements(cnorm.replacements)
        cds12_matrix = _reorder_outgroup(cds12_matrix, outgroup)
        cds12_path = output_dir / f"{prefix}.cds12{ext}"
        if not dry_run:
            _write_matrix(cds12_matrix, cds12_path, to_format, "NT")
        variants.append({
            "variant": "cds12", "path": str(cds12_path),
            "seq_type": "NT",
            "length": len(list(cds12_matrix.values())[0]) if cds12_matrix else 0,
        })

    # --- Compute per-variant stats ---
    variant_stats: list[dict[str, Any]] = []
    # always: original
    orig_stats = _compute_concat_stats(matrix, resolved_seq_type)
    variant_stats.append({
        "variant": "original",
        "seq_type": resolved_seq_type,
        "total_length": orig_stats["alignment_length"],
        "character_summary": orig_stats["character_summary"],
        "site_patterns": orig_stats["site_patterns"],
    })
    if recoding:
        rec_stats = _compute_concat_stats(recoded_matrix, recoded_seq_type)
        variant_stats.append({
            "variant": "recoded",
            "seq_type": recoded_seq_type,
            "total_length": rec_stats["alignment_length"],
            "character_summary": rec_stats["character_summary"],
            "site_patterns": rec_stats["site_patterns"],
        })
    if resolved_seq_type == "CODON" and translate_codon:
        tr_stats = _compute_concat_stats(translated_matrix, "AA")
        variant_stats.append({
            "variant": "translated",
            "seq_type": "AA",
            "total_length": tr_stats["alignment_length"],
            "character_summary": tr_stats["character_summary"],
            "site_patterns": tr_stats["site_patterns"],
        })
    if resolved_seq_type == "CODON" and exclude_codon3:
        c12_stats = _compute_concat_stats(cds12_matrix, "NT")
        variant_stats.append({
            "variant": "cds12",
            "seq_type": "NT",
            "total_length": c12_stats["alignment_length"],
            "character_summary": c12_stats["character_summary"],
            "site_patterns": c12_stats["site_patterns"],
        })

    if not dry_run and dropped:
        dropped_path = output_dir / "dropped_alignments.csv"
        with open(dropped_path, "w") as fh:
            fh.write("filename,n_taxa,occupancy_ratio,total_taxa\n")
            for entry in dropped:
                fh.write(f"{entry['filename']},{entry['n_taxa']},{entry['occupancy_ratio']},{entry['total_taxa']}\n")

    wall_time = time.time() - start_time
    payload = {
        "status": "success",
        "command": f"phyloai pretree concat --msa-dir {msa_dir}",
        "wall_time": round(wall_time, 3),
        "tool_versions": {},
        "params": {
            "msa_dir": str(msa_dir),
            "output_dir": str(output_dir),
            "prefix": prefix,
            "seq_type": resolved_seq_type,
            "taxa_occupancy": taxa_occupancy,
            "recoding": recoding,
            "outgroup": outgroup,
            "to_format": to_format,
            "translate_codon": translate_codon,
            "exclude_codon3": exclude_codon3,
            "dry_run": dry_run,
        },
        "key_results": {
            "n_taxa": len(all_taxa),
            "n_msa_input": len(msa_paths),
            "n_msa_used": len(kept_paths),
            "n_msa_dropped": len(dropped),
            "total_length": orig_stats["alignment_length"],
            "variants_produced": [v["path"] for v in variants],
        },
        "error": None,
        "data": {
            "character_summary": orig_stats["character_summary"],
            "site_patterns": orig_stats["site_patterns"],
            "variant_stats": variant_stats,
            "dropped_alignments": dropped,
            "per_taxon": orig_stats["per_taxon"],
            "per_gene_occupancy": [
                {
                    "gene": Path(path).name,
                    "n_present": len(msa_taxa_map[str(path)]),
                    "n_missing": len(all_taxa) - len(msa_taxa_map[str(path)]),
                    "occupancy_ratio": round(len(msa_taxa_map[str(path)]) / len(all_taxa), 4),
                }
                for path in kept_paths
            ],
            "variants": variants,
            "recoding_warnings": recoding_warnings,
            "normalization_replacements": all_normalization_replacements,
        },
    }

    if not dry_run:
        result_path = output_dir / "result.json"
        with open(result_path, "w") as fh:
            json.dump(payload, fh, indent=2)

        log_path = output_dir / "concat.log"
        log_lines = [
            "command=phyloai pretree concat",
            f"wall_time={round(wall_time, 3)}",
            "status=success",
            "exit_code=0",
            f"n_taxa={len(all_taxa)}",
            f"n_msa_input={len(msa_paths)}",
            f"n_msa_used={len(kept_paths)}",
            f"n_msa_dropped={len(dropped)}",
            f"total_length={orig_stats['alignment_length']}",
        ]
        for vs in variant_stats:
            log_lines.append(f"variant.{vs['variant']}.seq_type={vs['seq_type']}")
            log_lines.append(f"variant.{vs['variant']}.total_length={vs['total_length']}")
        log_path.write_text("\n".join(log_lines) + "\n")

    return payload
