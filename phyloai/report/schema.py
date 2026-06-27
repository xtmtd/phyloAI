"""Report data structures and report.json assembly."""

from __future__ import annotations

from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from phyloai import __version__

_PHASE_PREFIX: dict[str, int] = {
    "pretree": 3,
    "tree": 4,
    "posttree": 5,
}

_FIGURE_EXTENSIONS = {".pdf", ".png"}
_TABLE_EXTENSIONS = {".csv", ".tsv"}


def _get_phase_prefix(step_id: str) -> int:
    top = step_id.split(".")[0]
    return _PHASE_PREFIX.get(top, 99)


@dataclass
class ReportStep:
    step_id: str
    command: str
    status: str
    wall_time: float
    tool_versions: dict[str, str]
    params: dict[str, Any]
    key_results: dict[str, Any]
    methods_text: str = ""
    output_files: dict[str, dict[str, str]] = field(default_factory=dict)
    warnings: list[str] = field(default_factory=list)
    error: str | None = None

    def to_dict(self) -> dict[str, Any]:
        return {
            "step_id": self.step_id,
            "command": self.command,
            "status": self.status,
            "wall_time": self.wall_time,
            "tool_versions": self.tool_versions,
            "params": self.params,
            "key_results": self.key_results,
            "methods_text": self.methods_text,
            "output_files": self.output_files,
            "warnings": self.warnings,
            "error": self.error,
        }


@dataclass
class ReportRecord:
    run_dir: Path
    run_mode: str
    status: str
    phyloai_version: str = __version__
    generated_at: str = field(default_factory=lambda: datetime.now(timezone.utc).isoformat())
    steps: list[ReportStep] = field(default_factory=list)
    methods_paragraph: str = ""
    pipeline_summary: dict[str, Any] | None = None
    figures_index: list[dict[str, Any]] = field(default_factory=list)
    tables_index: list[dict[str, Any]] = field(default_factory=list)

    def to_dict(self) -> dict[str, Any]:
        methods_blocks = []
        for i, s in enumerate(self.steps):
            if s.methods_text and s.status == "success":
                methods_blocks.append({"step_id": s.step_id, "text": s.methods_text, "step_index": i})
        return {
            "phyloai_version": self.phyloai_version,
            "generated_at": self.generated_at,
            "run_dir": str(self.run_dir.absolute()),
            "run_mode": self.run_mode,
            "status": self.status,
            "pipeline_summary": self.pipeline_summary or {},
            "steps": [s.to_dict() for s in self.steps],
            "methods_paragraph": self.methods_paragraph,
            "methods_blocks": methods_blocks,
            "figures_index": self.figures_index,
            "tables_index": self.tables_index,
        }


def build_figures_index(steps: list[dict[str, Any]]) -> list[dict[str, Any]]:
    figures: list[dict[str, Any]] = []
    counters: dict[int, int] = {}
    for step in steps:
        step_id = step.get("step_id", "unknown")
        phase = _get_phase_prefix(step_id)
        output_files = step.get("output_files") or step.get("data", {}).get("output_files", {})
        for label, file_obj in output_files.items():
            if not isinstance(file_obj, dict):
                continue
            path = file_obj.get("path", "")
            ext = Path(path).suffix.lower()
            if ext in _FIGURE_EXTENSIONS:
                counters[phase] = counters.get(phase, 0) + 1
                figures.append({
                    "figure_id": f"Fig-{phase}.{counters[phase]}",
                    "step_id": step_id,
                    "label": label,
                    "caption": file_obj.get("description", label),
                    "description": file_obj.get("description", label),
                    "path": path,
                    "type": ext.lstrip("."),
                })
    return figures


def build_tables_index(steps: list[dict[str, Any]]) -> list[dict[str, Any]]:
    tables: list[dict[str, Any]] = []
    counters: dict[int, int] = {}
    for step in steps:
        step_id = step.get("step_id", "unknown")
        phase = _get_phase_prefix(step_id)
        output_files = step.get("output_files") or step.get("data", {}).get("output_files", {})
        for label, file_obj in output_files.items():
            if not isinstance(file_obj, dict):
                continue
            path = file_obj.get("path", "")
            ext = Path(path).suffix.lower()
            if ext in _TABLE_EXTENSIONS:
                counters[phase] = counters.get(phase, 0) + 1
                tables.append({
                    "table_id": f"Table-{phase}.{counters[phase]}",
                    "step_id": step_id,
                    "label": label,
                    "caption": file_obj.get("description", label),
                    "description": file_obj.get("description", label),
                    "path": path,
                    "type": ext.lstrip("."),
                })
    return tables


def assemble_report(
    discovered: dict[str, Any],
    run_dir: Path,
) -> dict[str, Any]:
    steps: list[dict[str, Any]] = []
    n_success = 0
    n_failed = 0
    n_skipped = 0
    total_wall_time = 0.0

    for raw_step in discovered["steps"]:
        step_id = raw_step["step_id"]
        status = raw_step["status"]

        if status == "success":
            n_success += 1
        elif status == "error":
            n_failed += 1
        total_wall_time += raw_step.get("wall_time", 0.0)

        output_files = raw_step.get("data", {}).get("output_files", {})
        # Purge entries that are not {path, description} dicts (legacy
        # bare ints like n_plots) and known-redundant keys.
        _SKIP_OF_KEYS = {"n_plots"}
        output_files = {
            k: v for k, v in output_files.items()
            if isinstance(v, dict) and "path" in v and k not in _SKIP_OF_KEYS
        }

        # Enrich key_results by merging values that some commands put
        # outside key_results: data.summary (convert dir, stats dir) or
        # flat data.* (stats single-file).
        key_results = dict(raw_step.get("key_results", {}))
        raw_data = raw_step.get("data", {})
        data_summary = raw_data.get("summary", {})
        for k, v in data_summary.items():
            if isinstance(v, (int, float, str)) and k not in key_results:
                key_results[k] = v
            elif isinstance(v, dict) and all(isinstance(x, (int, float)) for x in v.values()):
                for sk, sv in v.items():
                    fk = f"{k}_{sk}"
                    if fk not in key_results:
                        key_results[fk] = sv
            elif isinstance(v, list) and k not in key_results:
                key_results[k] = v
        _STRUCTURAL_KEYS = {
            "output_files", "files", "per_gene", "cmd", "tool_stderr",
            "tool_log", "summary", "variant_stats", "dropped_alignments",
            "per_taxon", "per_gene_occupancy", "skipped", "warnings",
            "character_summary", "site_patterns", "recoding_warnings",
            "normalization_replacements",
        }
        for k, v in raw_data.items():
            if k not in _STRUCTURAL_KEYS and isinstance(v, (int, float, str, bool)) and k not in key_results:
                key_results[k] = v
        # Merge concat-specific metrics from variant_stats[0] (original)
        if "gap_ratio" not in key_results or "pi_ratio" not in key_results:
            variants = raw_data.get("variant_stats", [])
            if variants:
                orig = variants[0]
                cs = orig.get("character_summary", {})
                sp = orig.get("site_patterns", {})
                if "gap_ratio" not in key_results and "gap_ratio" in cs:
                    key_results["gap_ratio"] = cs["gap_ratio"]
                if "pi_ratio" not in key_results:
                    pi = sp.get("parsimony_informative")
                    if isinstance(pi, dict):
                        key_results["pi_ratio"] = pi.get("ratio", 0)
        # Flatten nested scalar dicts already in key_results (e.g. trim's
        # length_before: {mean, min, max}).
        for k in list(key_results.keys()):
            v = key_results[k]
            if isinstance(v, dict) and all(isinstance(x, (int, float)) for x in v.values()):
                for sk, sv in v.items():
                    fk = f"{k}_{sk}"
                    if fk not in key_results:
                        key_results[fk] = sv

        methods_text = raw_step.get("methods_text", "")
        if status != "success":
            methods_text = ""

        step_record = {
            "step_id": step_id,
            "command": raw_step.get("command", ""),
            "status": status,
            "wall_time": raw_step.get("wall_time", 0.0),
            "tool_versions": raw_step.get("tool_versions", {}),
            "params": raw_step.get("params", {}),
            "key_results": key_results,
            "methods_text": methods_text,
            "output_files": output_files,
            "warnings": raw_step.get("warnings") or raw_step.get("data", {}).get("warnings", []),
            "error": raw_step.get("error"),
        }
        steps.append(step_record)

    figures_index = build_figures_index(steps)
    tables_index = build_tables_index(steps)

    # Each step gets its own paragraph — even same step_id with different
    # data (e.g. fna vs faa) are separate analyses.
    methods_blocks: list[dict[str, Any]] = []
    for i, s in enumerate(steps):
        txt = s["methods_text"]
        if txt:
            methods_blocks.append({"step_id": s["step_id"], "text": txt, "step_index": i})
    methods_paragraph = " ".join(b["text"] for b in methods_blocks)

    if n_failed == 0 and n_success > 0:
        overall_status = "complete"
    elif n_success == 0 and n_failed > 0:
        overall_status = "failed"
    else:
        overall_status = "partial"

    pipeline_summary = {
        "status": overall_status,
        "n_steps_total": len(steps),
        "n_steps_success": n_success,
        "n_steps_failed": n_failed,
        "n_steps_skipped": n_skipped,
        "total_wall_time": total_wall_time,
    }
    if discovered.get("pipeline_summary"):
        pipeline_summary.update(discovered["pipeline_summary"])

    report = {
        "phyloai_version": __version__,
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "run_dir": str(run_dir.absolute()),
        "run_mode": discovered["run_mode"],
        "status": overall_status,
        "pipeline_summary": pipeline_summary,
        "steps": steps,
        "methods_paragraph": methods_paragraph,
        "methods_blocks": methods_blocks,
        "figures_index": figures_index,
        "tables_index": tables_index,
    }

    return report
