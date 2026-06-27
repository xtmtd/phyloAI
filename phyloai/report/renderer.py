"""HTML report renderer using Jinja2 templates."""

from __future__ import annotations

import csv
import os
from io import StringIO
from pathlib import Path
from typing import Any

from jinja2 import Environment, FileSystemLoader, select_autoescape

_TEMPLATE_DIR = Path(__file__).parent / "html"
_CSV_EMBED_MAX_ROWS = 200
_CSV_EMBED_MAX_BYTES = 500_000


def _relative_path(abs_path: str, base_dir: Path) -> str:
    try:
        return os.path.relpath(abs_path, str(base_dir))
    except ValueError:
        return abs_path


def _scientific_params(params: dict[str, Any]) -> dict[str, Any]:
    technical = {
        "threads", "output_dir", "run_dir", "overwrite", "resume",
        "dry_run", "quiet", "seq_dir", "msa_dir", "tree_dir",
        "mafft_path", "magus_path", "trimal_path", "iqtree_path",
        "fasttree_path", "wastral_path", "tool_args",
        "input", "output", "input_format", "table_format",
    }
    return {k: v for k, v in params.items() if k not in technical}


def render_html(report: dict[str, Any], output_dir: Path) -> Path:
    """Render report.json to report.html using Jinja2.

    All file paths in the output are pre-computed as relative to output_dir
    so the template never needs to call path logic.
    """
    env = Environment(
        loader=FileSystemLoader(str(_TEMPLATE_DIR)),
        autoescape=select_autoescape(["html"]),
    )
    template = env.get_template("report.html.j2")

    figures = []
    for fig in report.get("figures_index", []):
        fig_copy = dict(fig)
        fig_copy["rel_path"] = _relative_path(fig["path"], output_dir)
        figures.append(fig_copy)

    tables = []
    for tbl in report.get("tables_index", []):
        tbl_copy = dict(tbl)
        tbl_copy["rel_path"] = _relative_path(tbl["path"], output_dir)
        tables.append(tbl_copy)

    steps = []
    for step in report.get("steps", []):
        step_copy = dict(step)
        step_copy["scientific_params"] = _scientific_params(step.get("params", {}))
        step_copy["n_params"] = len(step_copy["scientific_params"])
        of = {}
        for label, file_obj in step.get("output_files", {}).items():
            if not isinstance(file_obj, dict) or "path" not in file_obj:
                continue
            of[label] = dict(file_obj)
            of[label]["rel_path"] = _relative_path(file_obj["path"], output_dir)
        step_copy["output_files_rel"] = of

        # Embed small CSV tables for inline viewing
        csv_tables: list[dict[str, Any]] = []
        for label, fo in of.items():
            ext = Path(fo.get("path", "")).suffix.lower()
            if ext not in (".csv", ".tsv"):
                continue
            try:
                fpath = Path(fo["path"])
                if not fpath.is_file():
                    continue
                if fpath.stat().st_size > _CSV_EMBED_MAX_BYTES:
                    continue
                content = fpath.read_text(encoding="utf-8", errors="replace")
                delimiter = "\t" if ext == ".tsv" else ","
                reader = csv.reader(StringIO(content), delimiter=delimiter)
                rows = list(reader)
                if len(rows) > _CSV_EMBED_MAX_ROWS + 1:  # +1 for header
                    continue
                csv_tables.append({
                    "label": label,
                    "description": fo.get("description", label),
                    "headers": rows[0] if rows else [],
                    "rows": rows[1:] if len(rows) > 1 else [],
                    "n_rows": len(rows) - 1 if len(rows) > 1 else 0,
                    "rel_path": fo["rel_path"],
                })
            except Exception:
                continue
        step_copy["csv_tables"] = csv_tables

        steps.append(step_copy)

    all_files: list[dict[str, Any]] = []
    for s in steps:
        for label, fo in s.get("output_files_rel", {}).items():
            all_files.append({
                "step_id": s["step_id"],
                "label": label,
                "description": fo.get("description", "\u2014"),
                "path": fo["path"],
                "rel_path": fo["rel_path"],
                "type": Path(fo["path"]).suffix.lstrip(".") if fo.get("path") else "?",
            })

    # Build step_index → methods block number mapping for back-links
    step_to_methods: dict[int, int] = {}
    for mi, mb in enumerate(report.get("methods_blocks", []), start=1):
        si = mb.get("step_index")
        if isinstance(si, int):
            step_to_methods[si] = mi
    for si, s in enumerate(steps):
        s["methods_link"] = step_to_methods.get(si)

    html = template.render(
        report=report,
        figures=figures,
        tables=tables,
        steps=steps,
        all_files=all_files,
    )

    out_path = output_dir / "report.html"
    out_path.write_text(html, encoding="utf-8")
    return out_path
