"""Compositional Constraint Analysis (CCA)."""

from __future__ import annotations

import csv
import math
import shlex
import shutil
import time
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt

from phyloai.core.schema import write_result_json


@dataclass(frozen=True)
class SiteFrequency:
    site: int
    frequencies: tuple[float, ...]


@dataclass(frozen=True)
class SiteLikelihood:
    site: int
    lnl_tree1: float
    lnl_tree2: float


@dataclass(frozen=True)
class CcaRow:
    model: str
    site: int
    keff: float
    lnl_tree1: float
    lnl_tree2: float
    delta_lnl_tree2_tree1: float


def _read_nonempty_file(path: Path, source: str) -> list[str]:
    if not path.exists():
        raise ValueError(f"{source} file does not exist: {path}")
    if not path.is_file():
        raise ValueError(f"{source} path is not a file: {path}")
    try:
        lines = path.read_text(encoding="utf-8").splitlines()
    except OSError as exc:
        raise ValueError(f"could not read {source} file {path}: {exc}") from exc
    if not any(line.strip() and not line.lstrip().startswith("#") for line in lines):
        raise ValueError(f"{source} file is empty: {path}")
    return lines


def _parse_site(value: str | None, path: Path, line_number: int) -> int:
    try:
        return int(value or "")
    except ValueError as exc:
        raise ValueError(f"file {path}, line {line_number}: site must be an integer") from exc


def _finite_float(value: str | None, path: Path, line_number: int) -> float:
    try:
        number = float(value or "")
    except ValueError as exc:
        raise ValueError(f"file {path}, line {line_number}: value must be numeric") from exc
    if not math.isfinite(number):
        raise ValueError(f"file {path}, line {line_number}: value must be finite")
    return number


def _validate_one_based_sites(rows: list[Any], path: Path, source: str) -> list[Any]:
    if not rows:
        raise ValueError(f"{source} file is empty: {path}")
    sites = [row.site for row in rows]
    if len(sites) != len(set(sites)):
        raise ValueError(f"{source} file {path}: duplicate site identifier")
    if set(sites) != set(range(1, len(rows) + 1)):
        raise ValueError(f"{source} file {path}: site identifiers must be consecutive from 1")
    return sorted(rows, key=lambda row: row.site)


def parse_site_freq(path: Path) -> list[SiteFrequency]:
    """Parse one-based, 20-state amino-acid site-frequency rows."""
    rows: list[SiteFrequency] = []
    for line_number, line in enumerate(_read_nonempty_file(path, "site-frequency"), 1):
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        fields = line.split()
        if len(fields) != 21:
            raise ValueError(f"site-frequency file {path}, line {line_number}: exactly 20 frequencies are required")
        try:
            frequencies = tuple(float(value) for value in fields[1:])
        except ValueError as exc:
            raise ValueError(f"site-frequency file {path}, line {line_number}: frequencies must be numeric") from exc
        if not all(math.isfinite(value) for value in frequencies):
            raise ValueError(f"site-frequency file {path}, line {line_number}: frequencies must be finite")
        if any(value < 0 for value in frequencies):
            raise ValueError(f"site-frequency file {path}, line {line_number}: frequencies must be non-negative")
        if not math.isclose(sum(frequencies), 1.0, abs_tol=1e-6):
            raise ValueError(f"site-frequency file {path}, line {line_number}: frequencies must sum to 1")
        rows.append(SiteFrequency(_parse_site(fields[0], path, line_number), frequencies))
    return _validate_one_based_sites(rows, path, "site-frequency")


def parse_site_lnl(path: Path) -> list[SiteLikelihood]:
    """Parse a signal-lnl CSV table, retaining the required three fields."""
    if not path.exists():
        raise ValueError(f"site-likelihood file does not exist: {path}")
    if not path.is_file():
        raise ValueError(f"site-likelihood path is not a file: {path}")
    try:
        with path.open(newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle)
            required = {"site", "lnL_Tree1", "lnL_Tree2"}
            if reader.fieldnames is None or not required <= set(reader.fieldnames):
                raise ValueError(f"site-likelihood file {path}: header requires site, lnL_Tree1, lnL_Tree2")
            rows: list[SiteLikelihood] = []
            for number, row in enumerate(reader, 2):
                if None in row:
                    raise ValueError(f"site-likelihood file {path}, line {number}: too many fields")
                rows.append(SiteLikelihood(
                    _parse_site(row.get("site"), path, number),
                    _finite_float(row.get("lnL_Tree1"), path, number),
                    _finite_float(row.get("lnL_Tree2"), path, number),
                ))
    except csv.Error as exc:
        raise ValueError(f"malformed site-likelihood CSV {path}: {exc}") from exc
    except OSError as exc:
        raise ValueError(f"could not read site-likelihood file {path}: {exc}") from exc
    return _validate_one_based_sites(rows, path, "site-likelihood")


def validate_matching_sites(site_freq: list[SiteFrequency], model1: list[SiteLikelihood], model2: list[SiteLikelihood]) -> None:
    expected = [row.site for row in site_freq]
    if [row.site for row in model1] != expected or [row.site for row in model2] != expected:
        raise ValueError("site sets must match across --site-freq, --site-lnl1, and --site-lnl2")


def _keff(frequencies: tuple[float, ...]) -> float:
    return 1.0 / sum(value * value for value in frequencies)


def _model_names(model1_name: str, model2_name: str) -> tuple[str, str]:
    names = (model1_name.strip(), model2_name.strip())
    if not all(names) or names[0] == names[1]:
        raise ValueError("model names must be non-empty and distinct")
    return names


def build_cca_rows(site_freq: list[SiteFrequency], model1: list[SiteLikelihood], model2: list[SiteLikelihood], model1_name: str, model2_name: str) -> list[CcaRow]:
    names = _model_names(model1_name, model2_name)
    validate_matching_sites(site_freq, model1, model2)
    rows: list[CcaRow] = []
    for freq, first, second in zip(site_freq, model1, model2):
        for name, likelihood in zip(names, (first, second)):
            rows.append(CcaRow(name, freq.site, _keff(freq.frequencies), likelihood.lnl_tree1, likelihood.lnl_tree2, likelihood.lnl_tree2 - likelihood.lnl_tree1))
    return rows


def _keff_bin(keff: float) -> int:
    bin_id = 20 if math.isclose(keff, 20.0, abs_tol=1e-9) else math.floor(keff)
    if not 1 <= bin_id <= 20:
        raise ValueError(f"Keff is outside the supported range [1, 20]: {keff}")
    return bin_id


def summarize_bins(rows: list[CcaRow], model_names: tuple[str, str]) -> dict[str, list[float]]:
    sums = {name: [0.0] * 20 for name in model_names}
    for row in rows:
        sums[row.model][_keff_bin(row.keff) - 1] += row.delta_lnl_tree2_tree1
    return sums


def plot_cca(bin_sums: dict[str, list[float]], model_names: tuple[str, str], title: str, xlabel: str, ylabel: str, fig_width: float, fig_height: float, dpi: int, font_size: float, path: Path) -> None:
    values = [value for name in model_names for value in bin_sums[name]]
    minimum, maximum = min(values), max(values)
    ymin, ymax = min(0.0, minimum * 1.1), max(0.0, maximum * 1.1)
    if ymin == ymax:
        ymin, ymax = -1.0, 1.0
    figure, axis = plt.subplots(figsize=(fig_width, fig_height), dpi=dpi)
    axis.axhspan(ymin, 0, color="#ffdab9", alpha=0.5, zorder=0)
    axis.axhspan(0, ymax, color="lightblue", alpha=0.5, zorder=0)
    centres = [index + 0.5 for index in range(1, 21)]
    colours = {name: color for name, color in zip(sorted(model_names), ("#F8766D", "#00BFC4"))}
    for offset, name in zip((-0.25, 0.25), model_names):
        axis.bar([centre + offset for centre in centres], bin_sums[name], width=0.5, label=name, color=colours[name], zorder=2)
    for boundary in range(1, 21):
        axis.axvline(boundary, color="grey", linewidth=0.1, zorder=1)
    axis.axhline(0, color="black", linewidth=0.8, zorder=3)
    axis.set(xlim=(1, 21), ylim=(ymin, ymax), xticks=list(range(1, 21)), title=title, xlabel=xlabel, ylabel=ylabel)
    for spine in axis.spines.values():
        spine.set_color("black")
        spine.set_linewidth(0.8)
    axis.tick_params(labelsize=11, colors="black")
    axis.xaxis.label.set_size(11)
    axis.yaxis.label.set_size(11)
    axis.title.set_size(11)
    axis.grid(True, axis="y", color="lightgrey", linewidth=0.8)
    axis.grid(False, axis="x")
    legend = axis.legend(title=None, loc="upper right", bbox_to_anchor=(0.99, 0.9), borderaxespad=0, framealpha=0.5, facecolor="white", edgecolor="black", fontsize=font_size)
    legend.get_frame().set_linewidth(0.5)
    figure.tight_layout()
    plt.savefig(path, format="pdf", dpi=dpi)
    plt.close(figure)


def build_cca_command(site_freq: Path, site_lnl1: Path, site_lnl2: Path, model1_name: str, model2_name: str, title: str, xlabel: str, ylabel: str, fig_width: float, fig_height: float, dpi: int, font_size: float, output_dir: Path, overwrite: bool, dry_run: bool, quiet: bool) -> str:
    parts = ["phyloai", "posttree", "syserror", "cca", "--site-freq", str(site_freq.resolve()), "--site-lnl1", str(site_lnl1.resolve()), "--site-lnl2", str(site_lnl2.resolve()), "--output-dir", str(output_dir.resolve())]
    if model1_name != "model1": parts += ["--model1-name", model1_name]
    if model2_name != "model2": parts += ["--model2-name", model2_name]
    if title: parts += ["--title", title]
    if xlabel != "Effective number of amino acids": parts += ["--xlabel", xlabel]
    if ylabel != "Log-likelihood difference": parts += ["--ylabel", ylabel]
    if fig_width != 10: parts += ["--fig-width", str(fig_width)]
    if fig_height != 6: parts += ["--fig-height", str(fig_height)]
    if dpi != 300: parts += ["--dpi", str(dpi)]
    if font_size != 16: parts += ["--font-size", str(font_size)]
    if overwrite: parts.append("--overwrite")
    if dry_run: parts.append("--dry-run")
    if quiet: parts.append("--quiet")
    return shlex.join(parts)


def run_cca(site_freq: Path, site_lnl1: Path, site_lnl2: Path, model1_name: str = "model1", model2_name: str = "model2", title: str = "", xlabel: str = "Effective number of amino acids", ylabel: str = "Log-likelihood difference", fig_width: float = 10, fig_height: float = 6, dpi: int = 300, font_size: float = 16, output_dir: Path = Path("runs/posttree/syserror/cca"), overwrite: bool = False, dry_run: bool = False, quiet: bool = False) -> dict[str, Any]:
    start = time.monotonic()
    if min(fig_width, fig_height, dpi, font_size) <= 0:
        raise ValueError("figure width, height, DPI, and font size must be positive")
    frequencies = parse_site_freq(site_freq)
    first, second = parse_site_lnl(site_lnl1), parse_site_lnl(site_lnl2)
    names = _model_names(model1_name, model2_name)
    rows = build_cca_rows(frequencies, first, second, *names)
    sums = summarize_bins(rows, names)
    output_dir = output_dir.resolve()
    params = {"site_freq": str(site_freq.resolve()), "site_lnl1": str(site_lnl1.resolve()), "site_lnl2": str(site_lnl2.resolve()), "model1_name": names[0], "model2_name": names[1], "title": title, "xlabel": xlabel, "ylabel": ylabel, "fig_width": fig_width, "fig_height": fig_height, "dpi": dpi, "font_size": font_size, "output_dir": str(output_dir), "overwrite": overwrite, "dry_run": dry_run, "quiet": quiet}
    totals = {name: sum(row.delta_lnl_tree2_tree1 for row in rows if row.model == name) for name in names}
    payload: dict[str, Any] = {"status": "success", "command": build_cca_command(site_freq, site_lnl1, site_lnl2, *names, title, xlabel, ylabel, fig_width, fig_height, dpi, font_size, output_dir, overwrite, dry_run, quiet), "wall_time": 0.0, "tool_versions": {}, "params": params, "key_results": {"n_sites": len(frequencies), "keff_min": min(row.keff for row in rows), "keff_max": max(row.keff for row in rows), "keff_mean": sum(row.keff for row in rows) / len(rows), "models": list(names), "total_delta_lnl_tree2_tree1": totals, "bin_summaries": {name: [{"bin": index + 1, "delta_lnl_tree2_tree1": value} for index, value in enumerate(sums[name])] for name in names}}, "error": None, "data": {"cmd": [], "tool_stderr": "", "warnings": [], "output_files": {}}}
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
    table, figure = output_dir / "cca.csv", output_dir / "cca.pdf"
    with table.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(CcaRow.__dataclass_fields__))
        writer.writeheader()
        writer.writerows({
            **asdict(row),
            "lnl_tree1": f"{row.lnl_tree1:.4f}",
            "lnl_tree2": f"{row.lnl_tree2:.4f}",
            "delta_lnl_tree2_tree1": f"{row.delta_lnl_tree2_tree1:.4f}",
        } for row in rows)
    plot_cca(sums, names, title, xlabel, ylabel, fig_width, fig_height, dpi, font_size, figure)
    payload["data"]["output_files"] = {"cca_table": {"path": str(table.resolve()), "description": "Per-site compositional constraint analysis values"}, "cca_figure": {"path": str(figure.resolve()), "description": "Likelihood difference by effective amino-acid count"}}
    payload["wall_time"] = round(time.monotonic() - start, 3)
    write_result_json(payload, output_dir)
    if not quiet:
        for file in (table, figure, output_dir / "result.json"):
            print(f"Output: {file}")
    return payload
