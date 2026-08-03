"""Single-file gap-mask transfer from an original MSA onto a simulated MSA.

Transfers the gap (and optionally ambiguity) pattern of an original aligned
FASTA/PHYLIP/NEXUS onto a gap-free simulated MSA produced by
``phyloai posttree simulate alisim iqtree``, producing a single output MSA.
"""

from __future__ import annotations

import json
import shutil
import time as _time
from pathlib import Path
from typing import Any

from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from phyloai.core.formats import AlignmentFormat, detect_alignment_format
from phyloai.core.sequence_normalization import (
    AA_STANDARD,
    NT_STANDARD,
    detect_seq_type,
)

_BIOPYTHON_FORMAT = {
    AlignmentFormat.FASTA: "fasta",
    AlignmentFormat.PHYLIP: "phylip-relaxed",
    AlignmentFormat.NEXUS: "nexus",
    AlignmentFormat.PHYLIP_PAML: "phylip-relaxed",
}


def _read_msa(path: Path) -> MultipleSeqAlignment:
    """Parse one alignment file, raising ValueError on unparsable input."""
    fmt = detect_alignment_format(path)
    try:
        alignment = AlignIO.read(str(path), _BIOPYTHON_FORMAT[fmt])
    except Exception as exc:
        raise ValueError(f"unable to parse alignment file {path}: {exc}") from exc
    if len(alignment) == 0 or alignment.get_alignment_length() == 0:
        raise ValueError(f"alignment file {path} is empty")
    return alignment


def _validate_records(alignment: MultipleSeqAlignment, label: str) -> None:
    ids = [record.id for record in alignment]
    duplicates = sorted({name for name in ids if ids.count(name) > 1})
    if duplicates:
        raise ValueError(
            f"{label} alignment has duplicate taxon IDs: {', '.join(duplicates)}"
        )
    lengths = {len(record.seq) for record in alignment}
    if len(lengths) > 1:
        raise ValueError(
            f"{label} alignment sequences have unequal lengths: {sorted(lengths)}"
        )


def run_alisim_transfergaps(
    *,
    original_msa: Path,
    simulated_msa: Path,
    seq_type: str = "auto",
    exclude_ambiguity: bool = False,
    output_dir: Path = Path("runs/posttree/simulate/alisim/transfergaps"),
    overwrite: bool = False,
    dry_run: bool = False,
    quiet: bool = False,
) -> dict[str, Any]:
    """Replace selected simulated columns with the original per-taxon gap mask.

    Returns the standard result.json payload.  Raises ValueError on hard
    validation errors (missing/empty files, duplicate taxa, taxon-set or
    length mismatch).
    """
    run_start = _time.time()

    if seq_type not in {"auto", "AA", "NT"}:
        raise ValueError(f"--seq-type must be AA, NT, or auto, got {seq_type!r}")

    for label, path in (("--original-msa", original_msa),
                        ("--simulated-msa", simulated_msa)):
        if not path.exists():
            raise ValueError(f"{label} does not exist: {path}")
        if not path.is_file():
            raise ValueError(f"{label} is not a file: {path}")
        if path.stat().st_size == 0:
            raise ValueError(f"{label} is empty: {path}")

    original = _read_msa(original_msa)
    simulated = _read_msa(simulated_msa)
    _validate_records(original, "original")
    _validate_records(simulated, "simulated")

    original_ids = [record.id for record in original]
    simulated_ids = [record.id for record in simulated]
    missing = sorted(set(original_ids) - set(simulated_ids))
    extra = sorted(set(simulated_ids) - set(original_ids))
    if missing or extra:
        details = []
        if missing:
            details.append(f"missing from simulated: {', '.join(missing)}")
        if extra:
            details.append(f"only in simulated: {', '.join(extra)}")
        raise ValueError("taxon name mismatch between original and simulated "
                         f"MSAs ({'; '.join(details)})")

    original_length = original.get_alignment_length()
    simulated_length = simulated.get_alignment_length()
    if original_length != simulated_length:
        raise ValueError(
            f"length mismatch: original MSA has {original_length} columns but "
            f"simulated MSA has {simulated_length}. The simulated MSA must "
            "have the same alignment length as the original (AliSim --length "
            "should equal the original column count)."
        )

    if seq_type == "auto":
        resolved_seq_type = detect_seq_type([str(r.seq) for r in original])
    else:
        resolved_seq_type = seq_type
    standard = AA_STANDARD if resolved_seq_type == "AA" else NT_STANDARD

    simulated_by_id = {record.id: record for record in simulated}
    transferred: list[SeqRecord] = []
    n_masked_total = 0
    masked_per_taxon: list[int] = []
    for record in original:
        sim = simulated_by_id[record.id]
        original_chars = str(record.seq)
        sim_chars = list(str(sim.seq))
        n_masked = 0
        for index, char in enumerate(original_chars):
            if exclude_ambiguity:
                mask = char in ("-", ".")
            else:
                mask = char.upper() not in standard
            if mask:
                sim_chars[index] = "-"
                n_masked += 1
        masked_per_taxon.append(n_masked)
        n_masked_total += n_masked
        transferred.append(SeqRecord(Seq("".join(sim_chars)), id=record.id,
                                     name=record.id, description=""))

    output_dir = output_dir.resolve()
    output_name = f"{original_msa.stem}_transferred{original_msa.suffix or '.fa'}"
    output_path = output_dir / output_name

    if not dry_run:
        if output_dir.exists() and any(output_dir.iterdir()):
            if not overwrite:
                raise ValueError(
                    f"Output directory '{output_dir}' already exists and is "
                    "non-empty. Use --overwrite to replace it."
                )
            shutil.rmtree(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        _write_fasta(transferred, output_path)

    payload: dict[str, Any] = {
        "status": "success",
        "command": "phyloai posttree simulate alisim transfergaps "
                   f"--original-msa {original_msa} --simulated-msa {simulated_msa}",
        "wall_time": round(_time.time() - run_start, 3),
        "tool_versions": {},
        "params": {
            "original_msa": str(original_msa.resolve()),
            "simulated_msa": str(simulated_msa.resolve()),
            "seq_type": seq_type,
            "exclude_ambiguity": exclude_ambiguity,
            "output_dir": str(output_dir),
            "overwrite": overwrite,
            "dry_run": dry_run,
            "quiet": quiet,
        },
        "key_results": {
            "n_sequences": len(transferred),
            "alignment_length": original_length,
            "n_positions_masked": n_masked_total,
            "mean_positions_masked_per_taxon": (
                round(sum(masked_per_taxon) / len(masked_per_taxon))
                if masked_per_taxon else 0
            ),
            "detected_seq_type": resolved_seq_type,
        },
        "error": None,
        "error_category": None,
        "data": {
            "output_files": {
                "transferred_msa": {
                    "path": str(output_path),
                    "description": "Simulated MSA with gap pattern transferred "
                                   "from original",
                }
            }
        },
    }

    if not dry_run:
        with open(output_dir / "result.json", "w", encoding="utf-8") as fh:
            json.dump(payload, fh, indent=2)

    if not quiet:
        _print_summary(payload, dry_run)

    return payload


def _write_fasta(records: list[SeqRecord], path: Path) -> None:
    with open(path, "w", encoding="utf-8") as fh:
        for record in records:
            fh.write(f">{record.id}\n")
            seq = str(record.seq)
            for i in range(0, len(seq), 60):
                fh.write(seq[i:i + 60] + "\n")


def _print_summary(payload: dict[str, Any], dry_run: bool) -> None:
    click_echo = __import__("click").echo
    kr = payload["key_results"]
    data = payload["data"]["output_files"]["transferred_msa"]
    click_echo(
        "Gap transfer: "
        f"{Path(payload['params']['original_msa']).name} -> "
        f"{Path(payload['params']['simulated_msa']).name}"
    )
    click_echo(f"  Sequences: {kr['n_sequences']}")
    click_echo(
        f"  Non-standard positions masked: {kr['n_positions_masked']} "
        f"(mean {kr['mean_positions_masked_per_taxon']} per taxon)"
    )
    if not dry_run:
        click_echo(f"  Output: {data['path']}")
        click_echo(f"Result written to {Path(payload['params']['output_dir'])}/result.json")
