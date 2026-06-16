"""Shared helpers for matching MSA and tree files by logical locus names."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class PairingResult:
    paired: dict[str, tuple[Path | None, Path | None]]
    warnings: list[str]


def logical_msa_locus_name(path: Path) -> str:
    """Return the filename before its final suffix, or the full filename."""
    name = path.name
    stem = path.stem
    return stem if stem != name else name


def logical_tree_locus_candidates(path: Path) -> tuple[str | None, str | None]:
    """Return one-suffix and two-suffix logical reductions for a tree file."""
    name = path.name
    one_suffix = Path(name).stem
    if one_suffix == name:
        return name, None

    two_suffix = Path(one_suffix).stem
    if two_suffix == one_suffix:
        two_suffix = None

    return one_suffix, two_suffix


def pair_msa_and_tree_maps(msa_map: dict[str, Path], tree_paths: list[Path]) -> PairingResult:
    paired: dict[str, tuple[Path | None, Path | None]] = {
        locus: (msa_path, None) for locus, msa_path in msa_map.items()
    }
    warnings: list[str] = []

    for tree_path in sorted(tree_paths):
        candidates = [
            candidate for candidate in logical_tree_locus_candidates(tree_path)
            if candidate is not None and candidate in msa_map
        ]
        matched_loci = tuple(dict.fromkeys(candidates))

        if len(matched_loci) > 1:
            raise ValueError(
                f"ambiguous tree name {tree_path.name!r} matches loci "
                f"{', '.join(matched_loci)}"
            )

        if len(matched_loci) == 1:
            locus = matched_loci[0]
            if paired[locus][1] is not None:
                raise ValueError(f"duplicate tree file {tree_path.name!r} for locus {locus!r}")
            paired[locus] = (msa_map[locus], tree_path)
            continue

        locus = next(candidate for candidate in logical_tree_locus_candidates(tree_path) if candidate is not None)
        if locus in paired and paired[locus][1] is not None:
            raise ValueError(f"duplicate tree file {tree_path.name!r} for locus {locus!r}")
        paired[locus] = (None, tree_path)
        warnings.append(f"No matching MSA found for tree file {tree_path}")

    return PairingResult(paired=paired, warnings=warnings)
