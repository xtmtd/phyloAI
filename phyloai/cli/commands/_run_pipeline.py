"""Pipeline orchestration for phyloai run."""

from __future__ import annotations

import json
import os
import shutil
import time as _time
from pathlib import Path
from typing import Any

import click
from rich.console import Console

from phyloai.core.checkpoint import (
    CHECKPOINT_SCHEMA_VERSION,
    canonical_params_hash,
)

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

_STEP_DEFINITIONS: dict[tuple[str, str], list[dict[str, Any]]] = {
    ("supermatrix", "normal"): [
        {"name": "convert",      "subdir": "1-convert"},
        {"name": "align",        "subdir": "2-align"},
        {"name": "trim",         "subdir": "3-trim"},
        {"name": "filter_taper", "subdir": "4-filter"},
        {"name": "concat",       "subdir": "5-concat"},
        {"name": "tree",         "subdir": "6-tree"},
    ],
    ("supermatrix", "fast"): [
        {"name": "convert",      "subdir": "1-convert"},
        {"name": "align",        "subdir": "2-align"},
        {"name": "trim",         "subdir": "3-trim"},
        {"name": "concat",       "subdir": "5-concat"},
        {"name": "tree",         "subdir": "6-tree"},
    ],
    ("supertree", "normal"): [
        {"name": "convert",      "subdir": "1-convert"},
        {"name": "align",        "subdir": "2-align"},
        {"name": "trim",         "subdir": "3-trim"},
        {"name": "filter_taper", "subdir": "4-filter"},
        {"name": "genetrees",    "subdir": "5-genetrees"},
        {"name": "tree",         "subdir": "6-tree"},
    ],
    ("supertree", "fast"): [
        {"name": "convert",      "subdir": "1-convert"},
        {"name": "align",        "subdir": "2-align"},
        {"name": "trim",         "subdir": "3-trim"},
        {"name": "genetrees",    "subdir": "5-genetrees"},
        {"name": "tree",         "subdir": "6-tree"},
    ],
}

_STEP_TOOL_LABELS: dict[tuple[str, str, str], str] = {
    ("supermatrix", "normal", "convert"):      "Converting sequences (pretree convert)",
    ("supermatrix", "normal", "align"):        "Aligning sequences (MAFFT linsi)",
    ("supermatrix", "normal", "trim"):         "Trimming alignments (trimAl -automated1)",
    ("supermatrix", "normal", "filter_taper"): "Filtering error sites (TAPER)",
    ("supermatrix", "normal", "concat"):       "Concatenating matrix (pretree concat)",
    ("supermatrix", "normal", "tree"):         "Inferring species tree (IQ-TREE3, unpartitioned)",
    ("supermatrix", "fast", "convert"):        "Converting sequences (pretree convert)",
    ("supermatrix", "fast", "align"):          "Aligning sequences (MAFFT auto)",
    ("supermatrix", "fast", "trim"):           "Trimming alignments (trimAl -automated1)",
    ("supermatrix", "fast", "concat"):         "Concatenating matrix (pretree concat)",
    ("supermatrix", "fast", "tree"):           "Inferring species tree (FastTree --matrix)",
    ("supertree", "normal", "convert"):        "Converting sequences (pretree convert)",
    ("supertree", "normal", "align"):          "Aligning sequences (MAFFT linsi)",
    ("supertree", "normal", "trim"):           "Trimming alignments (trimAl -automated1)",
    ("supertree", "normal", "filter_taper"):   "Filtering error sites (TAPER)",
    ("supertree", "normal", "genetrees"):      "Building gene trees (IQ-TREE3 --msa-dir)",
    ("supertree", "normal", "tree"):           "Inferring species tree (wASTRAL mode 1)",
    ("supertree", "fast", "convert"):          "Converting sequences (pretree convert)",
    ("supertree", "fast", "align"):            "Aligning sequences (MAFFT auto)",
    ("supertree", "fast", "trim"):             "Trimming alignments (trimAl -automated1)",
    ("supertree", "fast", "genetrees"):        "Building gene trees (FastTree --msa-dir, fast)",
    ("supertree", "fast", "tree"):             "Inferring species tree (wASTRAL mode 1)",
}


# ---------------------------------------------------------------------------
# Custom exception for step failures (exit code 2 per Section 9.3)
# ---------------------------------------------------------------------------

class _RunStepError(click.ClickException):
    """Step failure — exit code 2 (tool failure) per Section 9.3."""

    exit_code = 2


class _EnvError(click.ClickException):
    """Missing tool — exit code 3 (environment error) per Section 9.3."""

    exit_code = 3


# ---------------------------------------------------------------------------
# Checkpoint helpers (plain dicts, not Checkpoint dataclass)
# ---------------------------------------------------------------------------

def _build_run_params(
    seq_dir: Path,
    mode: str,
    speed: str,
    threads: int,
    output_dir: Path,
) -> dict[str, Any]:
    return {
        "seq_dir": str(seq_dir.resolve()),
        "mode": mode,
        "speed": speed,
        "threads": threads,
        "output_dir": str(output_dir.resolve()),
    }


def _build_run_checkpoint(
    command_str: str,
    params: dict[str, Any],
    *,
    mode: str,
    speed: str,
) -> dict[str, Any]:
    import datetime as _dt

    now = _dt.datetime.now(_dt.timezone.utc).isoformat(timespec="seconds")
    step_defs = _STEP_DEFINITIONS[(mode, speed)]
    steps = [
        {"name": defn["name"], "status": "pending", "output_dir": None}
        for defn in step_defs
    ]
    return {
        "schema_version": CHECKPOINT_SCHEMA_VERSION,
        "step": "run",
        "command": command_str,
        "status": "running",
        "params_hash": canonical_params_hash(params),
        "params": params,
        "started_at": now,
        "updated_at": now,
        "completed_at": None,
        "steps": steps,
    }


def _load_run_checkpoint(checkpoint_path: Path) -> dict[str, Any]:
    if not checkpoint_path.exists():
        raise click.ClickException(
            f"run_checkpoint.json not found at {checkpoint_path}. "
            "Use --overwrite to start a clean run."
        )
    try:
        data = json.loads(checkpoint_path.read_text(encoding="utf-8"))
    except Exception as exc:
        raise click.ClickException(
            f"Malformed run_checkpoint.json: {exc}"
        ) from exc
    if data.get("schema_version") != CHECKPOINT_SCHEMA_VERSION:
        raise click.ClickException(
            f"Unsupported run_checkpoint.json schema version: "
            f"{data.get('schema_version')}. "
            "Use --overwrite to start a clean run."
        )
    if data.get("step") != "run":
        raise click.ClickException(
            "run_checkpoint.json belongs to a different command. "
            "Use --overwrite to start a clean run."
        )
    return data


def _validate_run_resume(checkpoint: dict[str, Any], current_hash: str) -> None:
    if checkpoint["params_hash"] != current_hash:
        raise click.ClickException(
            "Parameter mismatch: current parameters differ from the "
            "original run. Use --overwrite to start a clean run."
        )


def _save_json_atomic(data: dict[str, Any], path: Path, *, fsync: bool = False) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    with open(tmp, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2, ensure_ascii=False)
        f.flush()
        if fsync:
            os.fsync(f.fileno())
    os.replace(str(tmp), str(path))


def _save_run_checkpoint(checkpoint: dict[str, Any], path: Path, *, fsync: bool = False) -> None:
    import datetime as _dt
    checkpoint["updated_at"] = _dt.datetime.now(_dt.timezone.utc).isoformat(timespec="seconds")
    _save_json_atomic(checkpoint, path, fsync=fsync)


def _write_step_result_json(result: dict[str, Any], step_out: Path) -> None:
    """Write result.json for steps whose library functions don't do it themselves.
    Skips if the file already exists (step was already completed from resume)."""
    dest = step_out / "result.json"
    if dest.exists():
        return
    dest.write_text(
        json.dumps(result, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )


# ---------------------------------------------------------------------------
# Step dispatch helper
# ---------------------------------------------------------------------------

def _dispatch_step(
    *,
    checkpoint: dict[str, Any],
    checkpoint_path: Path,
    step_index: int,
    step_name: str,
    step_label: str,
    step_number: int,
    total_steps: int,
    console: Console,
    quiet: bool,
    runner: Any,
) -> tuple[dict[str, Any], bool]:
    """Run one step, updating checkpoint. Returns (result, skipped)."""
    step = checkpoint["steps"][step_index]

    if step["status"] == "success":
        result_path = Path(step["output_dir"]) / "result.json"
        if result_path.exists():
            try:
                data = json.loads(result_path.read_text())
                if data.get("status") == "success":
                    if not quiet:
                        console.print(
                            f"  [{step_number}/{total_steps}] {step_label}  "
                            f"[dim](already done, skipping)[/dim]"
                        )
                    return data, True
            except Exception:
                pass

    if not quiet:
        console.print(f"\n[bold][{step_number}/{total_steps}][/bold] {step_label} ...")

    step["status"] = "running"
    _save_run_checkpoint(checkpoint, checkpoint_path)

    try:
        result = runner()
    except click.ClickException:
        step["status"] = "failed"
        _save_run_checkpoint(checkpoint, checkpoint_path, fsync=True)
        raise
    except FileNotFoundError as exc:
        step["status"] = "failed"
        _save_run_checkpoint(checkpoint, checkpoint_path, fsync=True)
        raise _EnvError(f"Required tool not found: {exc}") from exc
    except Exception as exc:
        step["status"] = "failed"
        _save_run_checkpoint(checkpoint, checkpoint_path, fsync=True)
        raise _RunStepError(f"Step '{step_name}' failed: {exc}") from exc

    if isinstance(result, dict) and result.get("status") == "error":
        step["status"] = "failed"
        _save_run_checkpoint(checkpoint, checkpoint_path, fsync=True)
        raise _RunStepError(f"Step '{step_name}' returned status=error: {result.get('error', 'unknown error')}")

    step["status"] = "success"
    _save_run_checkpoint(checkpoint, checkpoint_path)
    return result, False


# ---------------------------------------------------------------------------
# Pipeline executor
# ---------------------------------------------------------------------------

def execute_pipeline(
    *,
    seq_dir: Path,
    mode: str,
    speed: str,
    output_dir: Path,
    threads: int,
    resume: bool,
    overwrite: bool,
    dry_run: bool,
    quiet: bool,
) -> None:
    console = Console(quiet=quiet)

    if resume and overwrite:
        raise click.ClickException("--resume and --overwrite are mutually exclusive.")

    checkpoint_path = output_dir / "run_checkpoint.json"
    params = _build_run_params(seq_dir, mode, speed, threads, output_dir)
    params_hash = canonical_params_hash(params)
    parts = [
        "phyloai", "run",
        "--seq-dir", str(seq_dir),
        "--mode", mode,
        "--speed", speed,
        "--output-dir", str(output_dir),
        "--threads", str(threads),
    ]
    if resume:
        parts.append("--resume")
    if overwrite:
        parts.append("--overwrite")
    if dry_run:
        parts.append("--dry-run")
    if quiet:
        parts.append("--quiet")
    command_str = " ".join(parts)

    # --- Dry-run: print step list (before any directory checks) ---
    step_defs = _STEP_DEFINITIONS[(mode, speed)]
    total = len(step_defs)

    if dry_run:
        console.print(f"\n[bold]Dry run — {mode} / {speed}[/bold]  ({total} steps)\n")
        for i, defn in enumerate(step_defs, 1):
            label = _STEP_TOOL_LABELS.get((mode, speed, defn["name"]), defn["name"])
            console.print(f"  [{i}/{total}] {defn['name']:15s}  {label}")
        console.print()
        return

    if resume:
        checkpoint = _load_run_checkpoint(checkpoint_path)
        _validate_run_resume(checkpoint, params_hash)
    else:
        if output_dir.exists() and any(output_dir.iterdir()):
            if overwrite:
                shutil.rmtree(output_dir)
            else:
                raise click.ClickException(
                    f"Output directory '{output_dir}' is non-empty. "
                    "Use --overwrite to replace it or --resume to continue."
                )
        output_dir.mkdir(parents=True, exist_ok=True)
        checkpoint = _build_run_checkpoint(command_str, params, mode=mode, speed=speed)
        for i, defn in enumerate(step_defs):
            checkpoint["steps"][i]["output_dir"] = str(output_dir / defn["subdir"])
        _save_run_checkpoint(checkpoint, checkpoint_path)
    run_start = _time.monotonic()
    all_tool_versions: dict[str, str] = {}

    def _write_error_result(error_msg: str) -> None:
        wall_time = round(_time.monotonic() - run_start, 3)
        payload: dict[str, Any] = {
            "status": "error",
            "command": command_str,
            "wall_time": wall_time,
            "tool_versions": all_tool_versions,
            "params": params,
            "key_results": {},
            "error": error_msg,
            "data": {
                "mode": mode,
                "speed": speed,
                "steps": [
                    {
                        "name": s["name"],
                        "status": s["status"],
                        "output_dir": s["output_dir"],
                        "result_json": str(Path(s["output_dir"]) / "result.json") if s["output_dir"] else None,
                    }
                    for s in checkpoint["steps"]
                ],
            },
        }
        result_file = output_dir / "result.json"
        result_file.write_text(
            json.dumps(payload, indent=2, ensure_ascii=False),
            encoding="utf-8",
        )

    try:
        step_map = {s["name"]: (i, s) for i, s in enumerate(checkpoint["steps"])}

        from phyloai.pretree.convert import convert_input
        from phyloai.pretree.align import run_align
        from phyloai.pretree.trim import run_trim
        from phyloai.pretree.filter import run_taper
        from phyloai.pretree.concat import run_concat
        from phyloai.tree.ml import run_fasttree
        from phyloai.tree.ml_iqtree import run_iqtree
        from phyloai.tree.msc import run_wastral

        # ---- Step 1: Convert ----
        def _step_resume(step_name: str) -> tuple[bool, bool]:
            """Return (can_resume: bool, is_resuming: bool) for a step."""
            idx = step_map[step_name][0]
            st = checkpoint["steps"][idx]
            is_resuming = st["status"] in ("interrupted", "running")
            return (not is_resuming, is_resuming)  # overwrite, resume

        def _get_step_out(step_name: str) -> Path:
            idx = step_map[step_name][0]
            return Path(checkpoint["steps"][idx]["output_dir"])

        step_out = _get_step_out("convert")
        step_out.mkdir(parents=True, exist_ok=True)

        def _run_convert() -> dict[str, Any]:
            return convert_input(
                input_path=seq_dir,
                output_dir=step_out,
                target_format="fasta",
                threads=threads,
                overwrite=True,
                quiet=quiet,
            )

        convert_result, _ = _dispatch_step(
            checkpoint=checkpoint,
            checkpoint_path=checkpoint_path,
            step_index=step_map["convert"][0],
            step_name="convert",
            step_label=_STEP_TOOL_LABELS[(mode, speed, "convert")],
            step_number=1,
            total_steps=total,
            console=console,
            quiet=quiet,
            runner=_run_convert,
        )
        all_tool_versions.update(convert_result.get("tool_versions") or {})
        _write_step_result_json(convert_result, step_out)
        converted_seqs_dir = step_out / "seqs"

        # ---- Step 2: Align ----
        step_out = _get_step_out("align")
        step_out.mkdir(parents=True, exist_ok=True)
        align_method = "linsi" if speed == "normal" else "auto"
        align_overwrite, align_resume = _step_resume("align")

        def _run_align_step() -> dict[str, Any]:
            return run_align(
                seq_dir=converted_seqs_dir,
                output_dir=step_out,
                method=align_method,
                seq_type="auto",
                threads=threads,
                overwrite=align_overwrite,
                resume=align_resume,
                quiet=quiet,
            )

        align_result, align_skipped = _dispatch_step(
            checkpoint=checkpoint,
            checkpoint_path=checkpoint_path,
            step_index=step_map["align"][0],
            step_name="align",
            step_label=_STEP_TOOL_LABELS[(mode, speed, "align")],
            step_number=2,
            total_steps=total,
            console=console,
            quiet=quiet,
            runner=_run_align_step,
        )
        all_tool_versions.update(align_result.get("tool_versions") or {})
        _write_step_result_json(align_result, step_out)
        aligned_seqs_dir = step_out / "seqs"

        # ---- Step 3: Trim ----
        step_out = _get_step_out("trim")
        step_out.mkdir(parents=True, exist_ok=True)
        trim_overwrite, trim_resume = _step_resume("trim")

        def _run_trim_step() -> dict[str, Any]:
            return run_trim(
                msa_dir=aligned_seqs_dir,
                output_dir=step_out,
                tool="trimal",
                trimal_method="automated1",
                threads=threads,
                overwrite=trim_overwrite,
                resume=trim_resume,
                quiet=quiet,
            )

        trim_result, trim_skipped = _dispatch_step(
            checkpoint=checkpoint,
            checkpoint_path=checkpoint_path,
            step_index=step_map["trim"][0],
            step_name="trim",
            step_label=_STEP_TOOL_LABELS[(mode, speed, "trim")],
            step_number=3,
            total_steps=total,
            console=console,
            quiet=quiet,
            runner=_run_trim_step,
        )
        all_tool_versions.update(trim_result.get("tool_versions") or {})
        _write_step_result_json(trim_result, step_out)
        trimmed_seqs_dir = step_out / "seqs"

        # Input to concat or gene-trees: after filter (normal) or trim (fast)
        filtered_seqs_dir = trimmed_seqs_dir
        step_number = 4

        # ---- Step 4: Filter (taper) — normal mode only ----
        if "filter_taper" in step_map:
            step_out = _get_step_out("filter_taper")
            step_out.mkdir(parents=True, exist_ok=True)

            def _run_taper_step() -> dict[str, Any]:
                return run_taper(
                    msa_dir=trimmed_seqs_dir,
                    output_dir=step_out,
                    threads=threads,
                    overwrite=True,
                    quiet=quiet,
                )

            taper_result, _ = _dispatch_step(
                checkpoint=checkpoint,
                checkpoint_path=checkpoint_path,
                step_index=step_map["filter_taper"][0],
                step_name="filter_taper",
                step_label=_STEP_TOOL_LABELS[(mode, speed, "filter_taper")],
                step_number=step_number,
                total_steps=total,
                console=console,
                quiet=quiet,
                runner=_run_taper_step,
            )
            all_tool_versions.update(taper_result.get("tool_versions") or {})
            filtered_seqs_dir = step_out / "seqs"
            step_number += 1

        n_genes_after_filter: int = 0
        if filtered_seqs_dir.exists():
            n_genes_after_filter = len(list(filtered_seqs_dir.glob("*.fa")))

        final_tree_path: str = ""
        matrix_length: int | None = None
        matrix_taxa: int | None = None

        # ---- Branch: supermatrix vs supertree ----
        if mode == "supermatrix":
            # Step 5: Concat
            step_out = _get_step_out("concat")
            step_out.mkdir(parents=True, exist_ok=True)

            def _run_concat() -> dict[str, Any]:
                return run_concat(
                    msa_dir=filtered_seqs_dir,
                    output_dir=step_out,
                    overwrite=True,
                    quiet=quiet,
                )

            concat_result, _ = _dispatch_step(
                checkpoint=checkpoint,
                checkpoint_path=checkpoint_path,
                step_index=step_map["concat"][0],
                step_name="concat",
                step_label=_STEP_TOOL_LABELS[(mode, speed, "concat")],
                step_number=step_number,
                total_steps=total,
                console=console,
                quiet=quiet,
                runner=_run_concat,
            )
            all_tool_versions.update(concat_result.get("tool_versions") or {})
            concat_kr = concat_result.get("key_results", {})
            matrix_length = concat_kr.get("total_length")
            matrix_taxa = concat_kr.get("n_taxa")
            matrix_file = step_out / "matrix.fa"
            step_number += 1

            # Step 6: Tree (iqtree normal / fasttree fast)
            step_out = _get_step_out("tree")
            step_out.mkdir(parents=True, exist_ok=True)
            tree_overwrite, tree_resume = _step_resume("tree")

            if speed == "normal":
                def _run_tree() -> dict[str, Any]:
                    return run_iqtree(
                        matrix=matrix_file,
                        output_dir=step_out,
                        threads=threads,
                        overwrite=tree_overwrite,
                        resume=tree_resume,
                        modelfinder="MFP",
                        quiet=quiet,
                    )
            else:
                def _run_tree() -> dict[str, Any]:
                    return run_fasttree(
                        matrix=matrix_file,
                        output_dir=step_out,
                        mode="normal",
                        threads=threads,
                        overwrite=tree_overwrite,
                        quiet=quiet,
                    )

            tree_result, _ = _dispatch_step(
                checkpoint=checkpoint,
                checkpoint_path=checkpoint_path,
                step_index=step_map["tree"][0],
                step_name="tree",
                step_label=_STEP_TOOL_LABELS[(mode, speed, "tree")],
                step_number=step_number,
                total_steps=total,
                console=console,
                quiet=quiet,
                runner=_run_tree,
            )
            all_tool_versions.update(tree_result.get("tool_versions") or {})
            _write_step_result_json(tree_result, step_out)
            final_tree_path = str(tree_result.get("data", {}).get("output") or "")

        else:  # supertree
            # Step 5: Gene trees (iqtree normal / fasttree fast)
            step_out = _get_step_out("genetrees")
            step_out.mkdir(parents=True, exist_ok=True)
            gt_overwrite, gt_resume = _step_resume("genetrees")

            if speed == "normal":
                def _run_genetrees() -> dict[str, Any]:
                    return run_iqtree(
                        msa_dir=filtered_seqs_dir,
                        output_dir=step_out,
                        threads=threads,
                        overwrite=gt_overwrite,
                        resume=gt_resume,
                        quiet=quiet,
                    )
            else:
                def _run_genetrees() -> dict[str, Any]:
                    return run_fasttree(
                        msa_dir=filtered_seqs_dir,
                        output_dir=step_out,
                        mode="fast",
                        threads=threads,
                        overwrite=gt_overwrite,
                        resume=gt_resume,
                        quiet=quiet,
                    )

            gt_result, _ = _dispatch_step(
                checkpoint=checkpoint,
                checkpoint_path=checkpoint_path,
                step_index=step_map["genetrees"][0],
                step_name="genetrees",
                step_label=_STEP_TOOL_LABELS[(mode, speed, "genetrees")],
                step_number=step_number,
                total_steps=total,
                console=console,
                quiet=quiet,
                runner=_run_genetrees,
            )
            all_tool_versions.update(gt_result.get("tool_versions") or {})
            _write_step_result_json(gt_result, step_out)
            trees_dir = step_out / "trees"
            step_number += 1

            # Step 6: Species tree (wASTRAL)
            step_out = _get_step_out("tree")
            step_out.mkdir(parents=True, exist_ok=True)

            def _run_wastral() -> dict[str, Any]:
                return run_wastral(
                    tree_dir=trees_dir,
                    output_dir=step_out,
                    mode=1,
                    threads=threads,
                    overwrite=True,
                    quiet=quiet,
                )

            wastral_result, _ = _dispatch_step(
                checkpoint=checkpoint,
                checkpoint_path=checkpoint_path,
                step_index=step_map["tree"][0],
                step_name="tree",
                step_label=_STEP_TOOL_LABELS[(mode, speed, "tree")],
                step_number=step_number,
                total_steps=total,
                console=console,
                quiet=quiet,
                runner=_run_wastral,
            )
            all_tool_versions.update(wastral_result.get("tool_versions") or {})
            _write_step_result_json(wastral_result, step_out)
            final_tree_path = str(
                wastral_result.get("data", {}).get("output_tree") or ""
            )

        # --- Validate final tree ---
        if not final_tree_path or not Path(final_tree_path).exists():
            tree_step_idx = step_map["tree"][0]
            checkpoint["steps"][tree_step_idx]["status"] = "failed"
            _save_run_checkpoint(checkpoint, checkpoint_path)
            raise _RunStepError(
                f"Final tree file not found after tree step completed: "
                f"{final_tree_path or '<none>'}"
            )

        # --- Mark run complete ---
        wall_time = round(_time.monotonic() - run_start, 3)

        checkpoint["status"] = "success"
        import datetime as _dt

        checkpoint["completed_at"] = _dt.datetime.now(
            _dt.timezone.utc
        ).isoformat(timespec="seconds")
        _save_run_checkpoint(checkpoint, checkpoint_path, fsync=True)

        # --- Write result.json ---
        n_input: int = 0
        if (output_dir / "1-convert" / "seqs").exists():
            n_input = len(list((output_dir / "1-convert" / "seqs").glob("*.fa")))

        result_payload: dict[str, Any] = {
            "status": "success",
            "command": command_str,
            "wall_time": wall_time,
            "tool_versions": all_tool_versions,
            "params": params,
            "key_results": {
                "n_input_genes": n_input,
                "n_genes_after_filter": n_genes_after_filter,
                "final_tree": final_tree_path,
                **({"matrix_length": matrix_length, "matrix_taxa": matrix_taxa}
                   if mode == "supermatrix" and matrix_length is not None and matrix_taxa is not None
                   else {}),
            },
            "error": None,
            "data": {
                "mode": mode,
                "speed": speed,
                "steps": [
                    {
                        "name": s["name"],
                        "status": s["status"],
                        "output_dir": s["output_dir"],
                        "result_json": str(Path(s["output_dir"]) / "result.json")
                        if s["output_dir"]
                        else None,
                    }
                    for s in checkpoint["steps"]
                ],
            },
        }

        result_file = output_dir / "result.json"
        result_file.write_text(
            json.dumps(result_payload, indent=2, ensure_ascii=False),
            encoding="utf-8",
        )

        # Final summary always printed (even in --quiet mode)
        click.echo(
            f"\n✓ Pipeline complete  [total: {wall_time:.0f}s]"
        )
        if final_tree_path:
            click.echo(f"  Species tree:  {final_tree_path}")
        click.echo(f"  Results:       {result_file}")

    except click.ClickException as exc:
        _write_error_result(exc.format_message())
        raise
    except KeyboardInterrupt:
        checkpoint["status"] = "interrupted"
        for s in checkpoint["steps"]:
            if s["status"] == "running":
                s["status"] = "interrupted"
        _save_run_checkpoint(checkpoint, checkpoint_path, fsync=True)
        raise
    except FileNotFoundError as exc:
        _write_error_result(f"Required tool not found: {exc}")
        raise _EnvError(f"Required tool not found: {exc}") from exc
    except Exception as exc:
        _write_error_result(str(exc))
        raise _RunStepError(str(exc)) from exc
