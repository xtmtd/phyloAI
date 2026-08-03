"""Gap-mask transfer from an original MSA onto one or more simulated MSAs.

Transfers the gap (and optionally ambiguity) pattern of an original aligned
FASTA/PHYLIP/NEXUS onto gap-free simulated MSAs produced by
``phyloai posttree simulate alisim iqtree``. Single mode takes one simulated
MSA file; batch mode takes a directory of simulated MSA files.
"""

from __future__ import annotations

import json
import shlex
import shutil
import time as _time
from pathlib import Path
from typing import Any

from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from phyloai.core.formats import FormatConverter
from phyloai.core.schema import COMMON_ALIGNMENT_EXTENSIONS
from phyloai.core.sequence_normalization import (
    AA_STANDARD,
    NT_STANDARD,
    detect_seq_type,
)

_FORMAT_CONVERTER = FormatConverter()


def _read_msa(path: Path) -> MultipleSeqAlignment:
    """Parse one alignment file (FASTA/PHYLIP/NEXUS/PHYLIP-PAML).

    Uses the shared FormatConverter so format detection and reading follow
    the same code path as `phyloai pretree convert`.
    """
    try:
        alignment = _FORMAT_CONVERTER.read(path)
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


def _list_alignment_files(directory: Path) -> list[Path]:
    """Sorted alignment files in a directory (non-recursive)."""
    if not directory.exists():
        raise ValueError(f"--simulated-dir does not exist: {directory}")
    if not directory.is_dir():
        raise ValueError(f"--simulated-dir is not a directory: {directory}")
    files = sorted(
        path for path in directory.iterdir()
        if path.is_file() and path.suffix.lower() in COMMON_ALIGNMENT_EXTENSIONS
    )
    if not files:
        raise ValueError(
            f"--simulated-dir contains no alignment files: {directory}"
        )
    return files


def _transfer_one(
    original: MultipleSeqAlignment,
    simulated: MultipleSeqAlignment,
    label: str,
    resolved_seq_type: str,
    exclude_ambiguity: bool,
) -> tuple[list[SeqRecord], int]:
    """Validate one simulated MSA against the original and build transferred records.

    Returns (records, total_masked_positions).  Raises ValueError on
    taxon-set or length mismatch.
    """
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
        raise ValueError(f"{label}: taxon name mismatch between original and "
                         f"simulated MSAs ({'; '.join(details)})")

    original_length = original.get_alignment_length()
    simulated_length = simulated.get_alignment_length()
    if original_length != simulated_length:
        raise ValueError(
            f"{label}: length mismatch: original MSA has {original_length} "
            f"columns but simulated MSA has {simulated_length}. The simulated "
            "MSA must have the same alignment length as the original (AliSim "
            "--length should equal the original column count)."
        )

    standard = AA_STANDARD if resolved_seq_type == "AA" else NT_STANDARD
    simulated_by_id = {record.id: record for record in simulated}
    transferred: list[SeqRecord] = []
    n_masked_total = 0
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
        n_masked_total += n_masked
        transferred.append(SeqRecord(Seq("".join(sim_chars)), id=record.id,
                                     name=record.id, description=""))
    return transferred, n_masked_total


def run_alisim_transfergaps(
    *,
    original_msa: Path,
    simulated_msa: Path | None = None,
    simulated_dir: Path | None = None,
    seq_type: str = "auto",
    exclude_ambiguity: bool = False,
    output_dir: Path = Path("runs/posttree/simulate/alisim/transfergaps"),
    overwrite: bool = False,
    dry_run: bool = False,
    quiet: bool = False,
) -> dict[str, Any]:
    """Replace selected simulated columns with the original per-taxon gap mask.

    Exactly one of ``simulated_msa`` (single mode) or ``simulated_dir``
    (batch mode) must be provided.  Returns the standard result.json payload.
    Raises ValueError on hard validation errors (missing/empty inputs,
    duplicate taxa, taxon-set or length mismatch).
    """
    run_start = _time.time()

    if seq_type not in {"auto", "AA", "NT"}:
        raise ValueError(f"--seq-type must be AA, NT, or auto, got {seq_type!r}")

    if (simulated_msa is None) == (simulated_dir is None):
        raise ValueError(
            "exactly one of --simulated-msa (single mode) or --simulated-dir "
            "(batch mode) must be provided"
        )

    if not original_msa.exists():
        raise ValueError(f"--original-msa does not exist: {original_msa}")
    if not original_msa.is_file():
        raise ValueError(f"--original-msa is not a file: {original_msa}")
    if original_msa.stat().st_size == 0:
        raise ValueError(f"--original-msa is empty: {original_msa}")

    if simulated_msa is not None:
        if not simulated_msa.exists():
            raise ValueError(f"--simulated-msa does not exist: {simulated_msa}")
        if not simulated_msa.is_file():
            raise ValueError(f"--simulated-msa is not a file: {simulated_msa}")
        if simulated_msa.stat().st_size == 0:
            raise ValueError(f"--simulated-msa is empty: {simulated_msa}")
        sim_tasks = [simulated_msa]
    else:
        sim_tasks = _list_alignment_files(simulated_dir)

    original = _read_msa(original_msa)
    _validate_records(original, "original")

    if seq_type == "auto":
        resolved_seq_type = detect_seq_type([str(r.seq) for r in original])
    else:
        resolved_seq_type = seq_type

    # Per simulated file: validate + transfer.  Fail fast on the first bad file.
    completed: list[tuple[Path, list[SeqRecord], int]] = []
    for sim_file in sim_tasks:
        sim = _read_msa(sim_file)
        _validate_records(sim, "simulated")
        transferred, n_masked = _transfer_one(
            original, sim, str(sim_file), resolved_seq_type, exclude_ambiguity,
        )
        completed.append((sim_file, transferred, n_masked))

    output_dir = output_dir.resolve()
    n_msas = len(completed)

    if simulated_dir is not None:
        output_files: dict[str, Any] = {
            "transferred_msas": {
                sim_file.stem: {
                    "path": str(output_dir / f"{sim_file.stem}.gaps.fa"),
                    "description": "Simulated MSA with gap pattern transferred "
                                   "from original",
                    "n_positions_masked": n_masked,
                }
                for sim_file, _transferred, n_masked in completed
            }
        }
        command_input: list[str] = [
            "--original-msa", str(original_msa),
            "--simulated-dir", str(simulated_dir),
        ]
        params_input = {
            "simulated_msa": None,
            "simulated_dir": str(simulated_dir.resolve()),
        }
    else:
        output_name = f"{original_msa.stem}.gaps.fa"
        output_path = output_dir / output_name
        output_files = {
            "transferred_msa": {
                "path": str(output_path),
                "description": "Simulated MSA with gap pattern transferred "
                               "from original",
            }
        }
        command_input = [
            "--original-msa", str(original_msa),
            "--simulated-msa", str(simulated_msa),
        ]
        params_input = {
            "simulated_msa": str(simulated_msa.resolve()),
            "simulated_dir": None,
        }

    command_parts = ["phyloai", "posttree", "simulate", "alisim", "transfergaps",
                     *command_input, "--seq-type", seq_type, "-o", str(output_dir)]
    if exclude_ambiguity:
        command_parts.append("--exclude-ambiguity")
    if overwrite:
        command_parts.append("--overwrite")
    if dry_run:
        command_parts.append("--dry-run")
    if quiet:
        command_parts.append("--quiet")
    full_command = shlex.join(command_parts)
    if not dry_run:
        if output_dir.exists() and any(output_dir.iterdir()):
            if not overwrite:
                raise ValueError(
                    f"Output directory '{output_dir}' already exists and is "
                    "non-empty. Use --overwrite to replace it."
                )
            shutil.rmtree(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        for sim_file, transferred, _n_masked in completed:
            if simulated_dir is not None:
                out_path = output_dir / f"{sim_file.stem}.gaps.fa"
            else:
                out_path = output_path
            _write_fasta(transferred, out_path)

    n_masked_total = sum(n_masked for _f, _t, n_masked in completed)
    n_taxa_total = sum(len(t) for _f, t, _n in completed)

    payload: dict[str, Any] = {
        "status": "success",
        "command": full_command,
        "wall_time": round(_time.time() - run_start, 3),
        "tool_versions": {},
        "params": {
            "original_msa": str(original_msa.resolve()),
            **params_input,
            "seq_type": seq_type,
            "exclude_ambiguity": exclude_ambiguity,
            "output_dir": str(output_dir),
            "overwrite": overwrite,
            "dry_run": dry_run,
            "quiet": quiet,
        },
        "key_results": {
            "n_msas": n_msas,
            "n_sequences": len(original),
            "alignment_length": original.get_alignment_length(),
            "n_positions_masked": n_masked_total,
            "mean_positions_masked_per_taxon": (
                round(n_masked_total / n_taxa_total) if n_taxa_total else 0
            ),
            "detected_seq_type": resolved_seq_type,
        },
        "error": None,
        "error_category": None,
        "data": {
            "output_files": output_files,
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
    data = payload["data"]["output_files"]
    sim_label = (Path(payload["params"]["simulated_msa"]).name
                 if payload["params"]["simulated_msa"]
                 else Path(payload["params"]["simulated_dir"]).name)
    click_echo(
        "Gap transfer: "
        f"{Path(payload['params']['original_msa']).name} -> {sim_label}"
    )
    click_echo(f"  MSAs processed: {kr['n_msas']}")
    click_echo(f"  Sequences: {kr['n_sequences']}")
    click_echo(
        f"  Non-standard positions masked: {kr['n_positions_masked']} "
        f"(mean {kr['mean_positions_masked_per_taxon']} per taxon)"
    )
    if not dry_run:
        if "transferred_msa" in data:
            click_echo(f"  Output: {data['transferred_msa']['path']}")
        else:
            outputs = data.get("transferred_msas", {})
            click_echo(f"  Outputs ({len(outputs)} files): {', '.join(outputs)}")
        click_echo(f"Result written to {Path(payload['params']['output_dir'])}/result.json")
