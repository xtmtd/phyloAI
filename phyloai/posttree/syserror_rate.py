"""Parse and normalize site-rate output for systematic-error analysis."""

from __future__ import annotations

import csv
import shlex
import math
import shutil
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from Bio.Align import MultipleSeqAlignment

from phyloai.core.formats import FormatConverter
from phyloai.core.schema import write_result_json


@dataclass(frozen=True)
class RateRow:
    site: int
    rate: float


def canonical_rates(rows: list[RateRow]) -> list[RateRow]:
    return sorted(rows, key=lambda row: (row.rate, row.site))


def subset_label(subset: str, fraction: float) -> str:
    percentage = format(fraction * 100, ".12g")
    return f"{subset}{percentage.removesuffix('.0')}"


def parse_fractions(value: str | None) -> list[float]:
    if value is None or not value.strip():
        raise ValueError("at least one fraction is required")
    fractions: list[float] = []
    labels: set[str] = set()
    for token in value.split(","):
        if not token.strip():
            raise ValueError("fractions must not contain empty values")
        try:
            fraction = float(token)
        except ValueError as exc:
            raise ValueError(f"invalid fraction: {token}") from exc
        if not math.isfinite(fraction) or not 0 < fraction <= 1:
            raise ValueError("fractions must be within (0, 1]")
        label = subset_label("subset", fraction)
        if fraction in fractions or label in labels:
            raise ValueError("fractions must be unique")
        fractions.append(fraction)
        labels.add(label)
    return fractions


def select_sites(ranked: list[RateRow], subset: str, fraction: float) -> list[int]:
    if subset not in {"slow", "fast"}:
        raise ValueError("subset must be 'slow' or 'fast'")
    count = math.ceil(len(ranked) * fraction)
    selected = ranked[:count] if subset == "slow" else ranked[-count:]
    return sorted(row.site for row in selected)


def read_matrix(path: Path) -> MultipleSeqAlignment:
    try:
        alignment = FormatConverter().read(path)
    except Exception as exc:
        raise ValueError(f"failed to read alignment {path}: {exc}") from exc
    if not alignment:
        raise ValueError(f"alignment has no records: {path}")
    if any(not record.seq for record in alignment):
        raise ValueError(f"alignment has empty sequences: {path}")
    if len({record.id for record in alignment}) != len(alignment):
        raise ValueError(f"alignment has duplicate record IDs: {path}")
    lengths = {len(record.seq) for record in alignment}
    if len(lengths) != 1:
        raise ValueError(f"alignment sequences have unequal lengths: {path}")
    return alignment


def write_subset_fasta(alignment: MultipleSeqAlignment, sites: list[int], path: Path) -> None:
    with open(path, "w", encoding="utf-8") as handle:
        for record in alignment:
            sequence = "".join(str(record.seq)[site - 1] for site in sites)
            handle.write(f">{record.description}\n")
            for start in range(0, len(sequence), 60):
                handle.write(f"{sequence[start:start + 60]}\n")


def rate_source(iqtree_rate: Path | None, pb_rate: Path | None) -> tuple[str, Path]:
    if (iqtree_rate is None) == (pb_rate is None):
        raise ValueError("exactly one of --iqtree-rate or --pb-rate is required")
    return ("iqtree", iqtree_rate) if iqtree_rate is not None else ("pb", pb_rate)  # type: ignore[return-value]


def build_rate_command(
    iqtree_rate: Path | None,
    pb_rate: Path | None,
    matrix: Path | None,
    subset: str | None,
    fraction: str | None,
    output_dir: Path,
    overwrite: bool,
    dry_run: bool,
    quiet: bool,
) -> str:
    """Build the reproducible CLI command recorded in rate result payloads."""
    parts = ["phyloai", "posttree", "syserror", "rate"]
    if iqtree_rate:
        parts += ["--iqtree-rate", str(iqtree_rate.resolve())]
    if pb_rate:
        parts += ["--pb-rate", str(pb_rate.resolve())]
    if matrix:
        parts += ["--matrix", str(matrix.resolve())]
    if subset:
        parts += ["--subset", subset]
    if fraction is not None:
        parts += ["--fraction", fraction]
    parts += ["--output-dir", str(output_dir.resolve())]
    if overwrite:
        parts.append("--overwrite")
    if dry_run:
        parts.append("--dry-run")
    if quiet:
        parts.append("--quiet")
    return shlex.join(parts)


def _read_lines(path: Path, source: str) -> list[str]:
    if not path.exists():
        raise ValueError(f"{source} rate file does not exist: {path}")
    if not path.is_file():
        raise ValueError(f"{source} rate path is not a file: {path}")
    try:
        lines = path.read_text().splitlines()
    except OSError as exc:
        raise ValueError(f"could not read {source} rate file {path}: {exc}") from exc
    if not any(line.strip() for line in lines):
        raise ValueError(f"{source} rate file is empty: {path}")
    return lines


def _parse_row(path: Path, source: str, line_number: int, site_text: str, rate_text: str) -> RateRow:
    try:
        site = int(site_text)
    except ValueError as exc:
        raise ValueError(f"{source} rate file {path}, line {line_number}: site must be an integer") from exc
    try:
        rate = float(rate_text)
    except ValueError as exc:
        raise ValueError(f"{source} rate file {path}, line {line_number}: rate must be numeric") from exc
    if not math.isfinite(rate):
        raise ValueError(f"{source} rate file {path}, line {line_number}: rate must be finite")
    if math.copysign(1.0, rate) < 0:
        raise ValueError(f"{source} rate file {path}, line {line_number}: rate must not be negative")
    return RateRow(site, rate)


def _validate_rows(
    path: Path, source: str, rows: list[RateRow], line_numbers: list[int], start: int
) -> list[RateRow]:
    if not rows:
        raise ValueError(f"{source} rate file has no data rows: {path}")
    sites: set[int] = set()
    for row, line_number in zip(rows, line_numbers):
        if row.site in sites:
            raise ValueError(f"{source} rate file {path}, line {line_number}: duplicate site identifier")
        sites.add(row.site)
    expected = set(range(start, start + len(rows)))
    if sites != expected:
        line_number = next(
            line for row, line in zip(rows, line_numbers) if row.site not in expected
        )
        raise ValueError(
            f"{source} rate file {path}, line {line_number}: "
            f"site identifiers must be consecutive from {start}"
        )
    return rows


def parse_iqtree_rate(path: Path) -> list[RateRow]:
    lines = _read_lines(path, "IQ-TREE")
    header: list[str] | None = None
    rows: list[RateRow] = []
    line_numbers: list[int] = []
    for line_number, line in enumerate(lines, 1):
        if not line.strip() or line.startswith("#"):
            continue
        fields = line.split("\t")
        if header is None:
            header = fields
            if "Site" not in header or "Rate" not in header:
                raise ValueError(f"IQ-TREE rate file {path}, line {line_number}: header requires Site and Rate columns")
            continue
        site_index = header.index("Site")
        rate_index = header.index("Rate")
        if len(fields) <= max(site_index, rate_index):
            raise ValueError(f"IQ-TREE rate file {path}, line {line_number}: malformed row")
        rows.append(_parse_row(path, "IQ-TREE", line_number, fields[site_index], fields[rate_index]))
        line_numbers.append(line_number)
    if header is None:
        raise ValueError(f"IQ-TREE rate file {path}: missing Site and Rate header")
    return _validate_rows(path, "IQ-TREE", rows, line_numbers, 1)


def parse_pb_rate(path: Path) -> list[RateRow]:
    rows: list[RateRow] = []
    line_numbers: list[int] = []
    for line_number, line in enumerate(_read_lines(path, "PhyloBayes"), 1):
        if not line.strip():
            continue
        fields = line.split()
        if len(fields) != 2:
            raise ValueError(f"PhyloBayes rate file {path}, line {line_number}: malformed row")
        rows.append(_parse_row(path, "PhyloBayes", line_number, fields[0], fields[1]))
        line_numbers.append(line_number)
    return [
        RateRow(row.site + 1, row.rate)
        for row in _validate_rows(path, "PhyloBayes", rows, line_numbers, 0)
    ]


def run_rate(
    iqtree_rate: Path | None,
    pb_rate: Path | None,
    matrix: Path | None = None,
    subset: str | None = None,
    fraction: str | None = None,
    output_dir: Path = Path("runs/posttree/syserror/rate"),
    overwrite: bool = False,
    dry_run: bool = False,
    quiet: bool = False,
) -> dict[str, Any]:
    """Rank site rates and optionally write slow or fast alignment subsets."""
    start = time.monotonic()
    source_name, source_path = rate_source(iqtree_rate, pb_rate)
    rows = parse_iqtree_rate(source_path) if source_name == "iqtree" else parse_pb_rate(source_path)
    ranked = canonical_rates(rows)

    if matrix is None and (subset is not None or fraction is not None):
        raise ValueError("--subset and --fraction require --matrix")
    if matrix is not None and fraction is None:
        raise ValueError("--matrix requires --fraction")

    original_subset = subset
    fractions: list[float] = []
    alignment: MultipleSeqAlignment | None = None
    subset = subset or "slow"
    if matrix is not None:
        if subset not in {"slow", "fast"}:
            raise ValueError("subset must be 'slow' or 'fast'")
        alignment = read_matrix(matrix)
        if alignment.get_alignment_length() != len(ranked):
            raise ValueError("matrix length does not match the number of rates")
        fractions = parse_fractions(fraction)

    output_dir = Path(output_dir).resolve()
    params = {
        "iqtree_rate": str(iqtree_rate) if iqtree_rate else None,
        "pb_rate": str(pb_rate) if pb_rate else None,
        "matrix": str(matrix) if matrix else None,
        "subset": subset if matrix else None,
        "fraction": fractions if matrix else None,
        "output_dir": str(output_dir),
        "overwrite": overwrite,
        "dry_run": dry_run,
        "quiet": quiet,
    }
    key_results: dict[str, Any] = {
        "rate_source": source_name,
        "n_sites": len(ranked),
        "min_rate": ranked[0].rate,
        "max_rate": ranked[-1].rate,
        "subsets": [],
    }
    for value in fractions:
        sites = select_sites(ranked, subset, value)
        key_results["subsets"].append({
            "subset": subset, "requested_fraction": value, "selected_sites": len(sites),
            "actual_fraction": len(sites) / len(ranked),
            "output_dir": str(output_dir / subset_label(subset, value)),
        })
    payload: dict[str, Any] = {
        "status": "success",
        "command": build_rate_command(
            iqtree_rate, pb_rate, matrix, original_subset, fraction, output_dir,
            overwrite, dry_run, quiet,
        ),
        "wall_time": 0.0,
        "tool_versions": {},
        "params": params,
        "key_results": key_results,
        "error": None,
        "data": {"cmd": [], "tool_stderr": "", "warnings": [], "output_files": {}},
    }
    if dry_run:
        payload["wall_time"] = round(time.monotonic() - start, 3)
        return payload

    if output_dir.exists() and not output_dir.is_dir():
        raise ValueError(f"output path is not a directory: {output_dir}")
    if output_dir.is_dir() and any(output_dir.iterdir()):
        if not overwrite:
            raise ValueError(f"output directory is not empty: {output_dir}")
        shutil.rmtree(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    output_files: dict[str, dict[str, str]] = {}
    rates_path = output_dir / "rates.csv"
    with rates_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["site", "rate"])
        writer.writeheader()
        writer.writerows({"site": row.site, "rate": str(row.rate)} for row in ranked)
    output_files["rates"] = {
        "path": str(rates_path.resolve()), "description": "Site rates sorted from slowest to fastest",
    }

    for value in fractions:
        label = subset_label(subset, value)
        directory = output_dir / label
        directory.mkdir()
        sites = select_sites(ranked, subset, value)
        positions_path = directory / "positions.txt"
        positions_path.write_text("".join(f"{site}\n" for site in sites), encoding="utf-8")
        matrix_path = directory / "matrix.fa"
        write_subset_fasta(alignment, sites, matrix_path)  # type: ignore[arg-type]
        output_files[f"{label}_positions"] = {
            "path": str(positions_path.resolve()), "description": f"Selected {subset} site positions",
        }
        output_files[f"{label}_matrix"] = {
            "path": str(matrix_path.resolve()), "description": f"Alignment containing selected {subset} sites",
        }
    payload["data"]["output_files"] = output_files
    payload["wall_time"] = round(time.monotonic() - start, 3)
    write_result_json(payload, output_dir)
    if not quiet:
        for file in [rates_path, *(Path(output_files[f"{subset_label(subset, value)}_matrix"]["path"]) for value in fractions), output_dir / "result.json"]:
            print(file)
    return payload
