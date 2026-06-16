# pretree filter Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Implement `phyloai pretree filter` as a Click group with four subcommands: `taper`, `treeshrink`, `metrics`, `cluster`.

**Architecture:** Core library `phyloai/pretree/filter.py` with four public entry points (`run_taper`, `run_treeshrink`, `run_metrics_filter`, `run_cluster_filter`). CLI layer in `phyloai/cli/commands/pretree.py` registers a `filter` Click group. TAPER follows the existing `trim.py` checkpoint/resume pattern with `ProcessPoolExecutor`. Shared infra: `core/file_matching.py`, `core/checkpoint.py`, `pretree/checkpoint_helpers.py`.

**Tech Stack:** Python stdlib, BioPython, numpy, scikit-learn, matplotlib, scipy. Optional `umap-learn`. External: Julia (bundled TAPER), `run_treeshrink.py` (user-installed).

**Output directories:**
```
runs/pretree/filter/taper/       # seqs/, filter.log, result.json
                                 # checkpoint.json (internal resume state; not user-facing report)
runs/pretree/filter/treeshrink/  # trees/, seqs/ (opt), filter.log, result.json
runs/pretree/filter/metrics/     # decision CSVs, filter.log, result.json
runs/pretree/filter/cluster/     # plots, CSVs, filter.log, result.json
```

Note: `checkpoint.json` is produced only by `filter taper` (which supports `--resume`). It is an internal recovery file following the existing `pretree align`/`trim` pattern, not a user-facing report artifact per design §3.1.

**Task ordering is important:** filter.py must exist before CLI subcommands are registered, because CLI decorators import from it.

---

## Phase 1: Foundation

### Task 1: Add `scan_msa_dir` and `scan_tree_dir` to `core/file_matching.py`

**Files:**
- Modify: `phyloai/core/file_matching.py`
- Test: `tests/core/test_file_matching.py`

**Policy note (per MAIN §9.7):** Per the global file-matching policy, filename suffixes are only dot-separated naming boundaries; hard-coded suffix whitelists must not gate file acceptance. `scan_msa_dir` and `scan_tree_dir` therefore scan ALL regular non-empty files and derive logical locus names from their dot-separated filename stems. Downstream format detection / parser validation (per §9.10) determines whether each entry is a valid alignment or tree.

- [ ] **Step 1: Write tests**

```python
# tests/core/test_file_matching.py — append after existing tests

from pathlib import Path


def test_scan_msa_dir_returns_locus_map(tmp_path):
    from phyloai.core.file_matching import scan_msa_dir
    (tmp_path / "gene1.fa").write_text(">a\nACGT\n")
    (tmp_path / "gene2.FASTA").write_text(">b\nACGT\n")
    (tmp_path / "gene3.tre").write_text("(a,b);")  # also scanned; validity decided downstream
    (tmp_path / "notes.txt").write_text("hello")    # also scanned; validity decided downstream
    (tmp_path / "subdir").mkdir()                   # directories skipped

    result = scan_msa_dir(tmp_path)
    # All regular files scanned; suffix is not the gate
    assert "gene1" in result
    assert "gene2" in result
    assert "gene3" in result
    assert "notes" in result
    assert len(result) == 4


def test_scan_tree_dir_returns_locus_map(tmp_path):
    from phyloai.core.file_matching import scan_tree_dir
    (tmp_path / "gene1.treefile").write_text("(a,b);")
    (tmp_path / "gene2.tre").write_text("(c,d);")
    (tmp_path / "gene3.fa").write_text(">a\nACGT\n")  # also scanned

    result = scan_tree_dir(tmp_path)
    assert "gene1" in result
    assert "gene2" in result
    assert "gene3" in result
    assert len(result) == 3


def test_scan_tree_dir_handles_ambiguous_one_two_suffix(tmp_path):
    from phyloai.core.file_matching import scan_tree_dir
    (tmp_path / "gene.v1.treefile").write_text("(a,b);")
    result = scan_tree_dir(tmp_path)
    assert len(result) == 1


def test_scan_msa_dir_empty_on_nonexistent():
    from phyloai.core.file_matching import scan_msa_dir
    assert scan_msa_dir(Path("/nonexistent")) == {}
```

- [ ] **Step 2: Run test — verify FAIL**

```bash
pytest tests/core/test_file_matching.py::test_scan_msa_dir_returns_locus_map -v
```

- [ ] **Step 3: Implement helpers (no extension whitelist)**

```python
# phyloai/core/file_matching.py — append after existing functions


def scan_msa_dir(path: Path) -> dict[str, Path]:
    """Scan a directory for potential MSA files, returning ``{logical_locus: path}``.

    Per the global file-matching policy (MAIN §9.7): filename suffixes are
    only dot-separated naming boundaries.  This helper does NOT use a
    hard-coded suffix whitelist — it scans every regular non-empty file.
    Downstream format detection and parser-level validation decide whether
    each file is a valid alignment.
    """
    if not path.exists() or not path.is_dir():
        return {}
    result: dict[str, Path] = {}
    for entry in sorted(path.iterdir(), key=lambda p: p.name):
        if not entry.is_file() or entry.stat().st_size == 0:
            continue
        locus = logical_msa_locus_name(entry)
        result[locus] = entry
    return result


def scan_tree_dir(path: Path) -> dict[str, Path]:
    """Scan a directory for potential tree files, returning ``{logical_locus: path}``.

    Same suffix-agnostic policy as :func:`scan_msa_dir`.  Logical locus
    names are derived from one- or two-suffix dot reductions of the
    filename.
    """
    if not path.exists() or not path.is_dir():
        return {}
    result: dict[str, Path] = {}
    for entry in sorted(path.iterdir(), key=lambda p: p.name):
        if not entry.is_file() or entry.stat().st_size == 0:
            continue
        locus, _ = logical_tree_locus_candidates(entry)
        if locus and locus not in result:
            result[locus] = entry
        elif locus and locus in result:
            _, candidate2 = logical_tree_locus_candidates(entry)
            if candidate2 and candidate2 not in result:
                result[candidate2] = entry
    return result
```

- [ ] **Step 4: Run tests — verify PASS**

```bash
pytest tests/core/test_file_matching.py -v
```

- [ ] **Step 5: Review checkpoint**

```bash
pytest tests/core/test_file_matching.py -v
```

---

### Task 2: Add `scikit-learn` dependency

**Files:**
- Modify: `pyproject.toml`

- [ ] **Step 1: Add to dependencies**

```toml
# pyproject.toml — in [project.dependencies], add:
"scikit-learn>=1.3.0",
```

- [ ] **Step 2: Install**

```bash
pip install -e ".[dev]"
```

- [ ] **Step 3: Review checkpoint**

```bash
python -c "import sklearn; print(sklearn.__version__)"
```

---

### Task 3: Create `phyloai/pretree/filter.py` skeleton

**Files:**
- Create: `phyloai/pretree/filter.py`

- [ ] **Step 1: Write skeleton with all shared helpers**

Full content of `phyloai/pretree/filter.py`:

```python
"""Marker-level filtering: TAPER, TreeShrink, metric rules, and clustering."""

from __future__ import annotations

import csv
import datetime
import json
import shlex
import shutil
import subprocess
import sys
import tempfile
import time
from collections.abc import Callable
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Any

from Bio import SeqIO
from rich.table import Table

from phyloai.core.checkpoint import (
    Checkpoint,
    load_checkpoint,
    save_checkpoint_atomic,
    validate_resume_params,
)
from phyloai.core.env import ToolEnv, ToolInfo
from phyloai.core.file_matching import (
    logical_msa_locus_name,
    pair_msa_and_tree_maps,
    scan_msa_dir,
    scan_tree_dir,
)
from phyloai.core.runner import Runner
from phyloai.core.sequence_normalization import (
    resolve_seq_type,
    validate_codon_msa,
)
from phyloai.core.sequence_output_validation import validate_fasta_output
from phyloai.pretree.checkpoint_helpers import (
    build_initial_checkpoint,
    mark_task,
    plan_resume,
)

_CHECKPOINT_FLUSH_INTERVAL = 2.0


# --- Shared output helpers ---

def _write_csv_table(rows: list[dict], path: Path, columns: list[str], delimiter: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=columns, delimiter=delimiter, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def _table_delimiter(table_format: str) -> str:
    return "\t" if table_format == "tsv" else ","


def _table_suffix(table_format: str) -> str:
    return ".tsv" if table_format == "tsv" else ".csv"


def _write_result_json(payload: dict, output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    with open(output_dir / "result.json", "w") as fh:
        json.dump(payload, fh, indent=2)


def _common_output_conflict(output_dir: Path, overwrite: bool) -> None:
    if output_dir.exists() and any(output_dir.iterdir()):
        if not overwrite:
            raise ValueError(
                f"Output directory '{output_dir}' already exists and is non-empty. "
                "Use --overwrite to replace it."
            )
        shutil.rmtree(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)


def _write_filter_log(output_dir: Path, command: str, wall_time: float,
                      tool_versions: dict, success: bool) -> None:
    log_path = output_dir / "filter.log"
    with open(log_path, "a") as fh:
        fh.write(f"# {command}\n")
        fh.write(f"# Started: {datetime.datetime.now().isoformat()}\n")
        fh.write(f"# Tool versions: {json.dumps(tool_versions)}\n")
        fh.write(f"# Wall time: {wall_time:.1f}s\n")
        fh.write(f"# Exit code: {'0' if success else '1'}\n")
        fh.write("---\n")


def render_filter_summary_table(summary: dict) -> Table:
    table = Table(title="Filter Results")
    table.add_column("Metric", style="cyan")
    table.add_column("Value", style="green")
    for key, value in summary.items():
        table.add_row(str(key), str(value))
    return table


def _compute_retained_msa_stats(msa_paths: list[Path]) -> dict:
    if not msa_paths:
        return {"n_msa": 0, "total_length": 0, "mean_length": 0,
                "min_length": 0, "max_length": 0, "mean_taxa": 0}
    lengths = []
    taxa_counts = []
    for p in msa_paths:
        try:
            records = list(SeqIO.parse(str(p), "fasta"))
            if records:
                lengths.append(len(records[0].seq))
                taxa_counts.append(len(records))
        except Exception:
            continue
    if not lengths:
        return {"n_msa": 0, "total_length": 0, "mean_length": 0,
                "min_length": 0, "max_length": 0, "mean_taxa": 0}
    return {
        "n_msa": len(lengths),
        "total_length": sum(lengths),
        "mean_length": round(sum(lengths) / len(lengths), 2),
        "min_length": min(lengths),
        "max_length": max(lengths),
        "mean_taxa": round(sum(taxa_counts) / len(taxa_counts), 2),
    }
```

- [ ] **Step 2: Verify it imports cleanly**

```bash
python -c "from phyloai.pretree.filter import _write_csv_table, render_filter_summary_table; print('OK')"
```

- [ ] **Step 3: Review checkpoint**

```bash
python -c "from phyloai.pretree.filter import _write_csv_table, render_filter_summary_table; print('OK')"
```

---

### Task 4: Register `filter` Click group in CLI (minimal, no subcommands yet)

**Files:**
- Modify: `phyloai/cli/commands/pretree.py`

- [ ] **Step 1: Update `_PretreeGroup.list_commands`**

```python
# phyloai/cli/commands/pretree.py:48 — change line
def list_commands(self, ctx: click.Context) -> list[str]:
    return ["convert", "stats", "align", "trim", "metrics", "filter", "concat"]
```

- [ ] **Step 2: Add filter Click group at end of file**

After the `metrics_group` registration, append:

```python
# ---------------------------------------------------------------------------
# filter group
# ---------------------------------------------------------------------------

# NOTE: Only import what exists now. run_* imports will be added in later tasks
# as each function is implemented. For now, just the shared summary renderer.
from phyloai.pretree.filter import render_filter_summary_table


@click.group("filter", help="Marker-level filtering: TAPER masking, TreeShrink pruning, metric-rule filtering, cluster-based exploration.")
def filter_group() -> None:
    pass


pretree.add_command(filter_group)
```

- [ ] **Step 3: Verify CLI help works**

```bash
python -m phyloai pretree filter --help
```

Expected: shows `Commands:` with no subcommands yet (or empty group help).

- [ ] **Step 4: Review checkpoint**

```bash
python -m phyloai pretree filter --help
```

---

## Phase 2: TAPER Subcommand

### Task 5: Implement TAPER command builder and per-locus worker

**Files:**
- Modify: `phyloai/pretree/filter.py`

- [ ] **Step 1: Append TAPER command builder and worker**

```python
# --- TAPER --- (append to phyloai/pretree/filter.py)

_TAPER_CUTOFF_DEFAULT = 3
_TAPER_MANAGED_FLAGS = {"-m", "-a", "-c", "-l"}
_TAPER_NT_CMD_EXTRA = ["-m", "N", "-a", "N"]


def _build_taper_cmd(
    input_file: Path, output_file: Path, seq_type: str, cutoff: int,
    julia_exe: str, taper_script: str, tool_args: str | None,
) -> list[str]:
    cmd = [julia_exe, taper_script, "-c", str(cutoff)]
    if seq_type == "NT":
        cmd.extend(_TAPER_NT_CMD_EXTRA)
    if tool_args:
        extra = shlex.split(tool_args)
        for flag in _TAPER_MANAGED_FLAGS:
            if flag in extra:
                raise ValueError(f"Flag {flag!r} is managed by PhyloAI; remove from --tool-args.")
        cmd.extend(extra)
    cmd.append(str(input_file))
    return cmd


def _run_taper_one(
    input_file: Path, output_file: Path, seq_type: str, cutoff: int,
    julia_exe: str, taper_script: str, tool_args: str | None, runner: Runner,
) -> dict:
    """Run TAPER on one MSA; count only newly-introduced X masks."""
    cmd = _build_taper_cmd(input_file, output_file, seq_type, cutoff, julia_exe, taper_script, tool_args)
    with open(output_file, "w") as fh:
        proc = subprocess.run(cmd, stdout=fh, stderr=subprocess.PIPE, text=True, timeout=86400)
    new_mask_count = 0
    if proc.returncode == 0 and output_file.exists() and seq_type == "AA":
        # Count only newly introduced X by comparing output vs input
        in_recs = {rec.id: str(rec.seq) for rec in SeqIO.parse(str(input_file), "fasta")}
        out_recs = {rec.id: str(rec.seq) for rec in SeqIO.parse(str(output_file), "fasta")}
        for taxon in in_recs:
            if taxon in out_recs:
                for i, (in_ch, out_ch) in enumerate(zip(in_recs[taxon], out_recs[taxon])):
                    if in_ch != "X" and out_ch == "X":
                        new_mask_count += 1
    return {
        "locus": input_file.stem,
        "status": "success" if proc.returncode == 0 else "failed",
        "returncode": proc.returncode,
        "cmd": " ".join(cmd),
        "stderr": proc.stderr[:500] if proc.stderr else "",
        "new_masked_sites": new_mask_count,
        "output": str(output_file),
    }


def _project_taper_masks_to_cds(
    aa_original_path: Path, aa_masked_path: Path,
    nt_input_path: Path, nt_output_path: Path,
) -> dict:
    """Project only TAPER-introduced AA 'X' masks to CDS as 'NNN'.

    Per design §5.4: original X → X = no CDS change;
    standard AA → X = replace codon with NNN;
    X → non-X = warning/failure.
    """
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord

    aa_orig = {rec.id: str(rec.seq) for rec in SeqIO.parse(str(aa_original_path), "fasta")}
    aa_masked = {rec.id: str(rec.seq) for rec in SeqIO.parse(str(aa_masked_path), "fasta")}
    nt_recs = {rec.id: str(rec.seq) for rec in SeqIO.parse(str(nt_input_path), "fasta")}

    if aa_orig.keys() != nt_recs.keys() or aa_masked.keys() != nt_recs.keys():
        raise ValueError(f"AA/NT taxa mismatch for {aa_masked_path.stem}")
    length = len(next(iter(aa_orig.values())))
    for nt_seq in nt_recs.values():
        if len(nt_seq) != length * 3:
            raise ValueError("NT alignment length != AA length * 3")

    projected = 0
    warnings_list: list[str] = []
    new_nt_records = []

    for taxon in aa_orig:
        orig_seq = aa_orig[taxon]
        masked_seq = aa_masked[taxon]
        nt_chars = list(nt_recs[taxon])
        for i, (orig_ch, mask_ch) in enumerate(zip(orig_seq, masked_seq)):
            codon_start = i * 3
            if orig_ch == "X" and mask_ch == "X":
                # Original ambiguity — no CDS change
                pass
            elif orig_ch in _STANDARD_AA and mask_ch == "X":
                # TAPER introduced a new mask — replace codon with NNN
                original_codon = "".join(nt_chars[codon_start:codon_start + 3])
                if original_codon not in ("---", "NNN"):
                    nt_chars[codon_start:codon_start + 3] = ["N", "N", "N"]
                    projected += 1
            elif orig_ch == "X" and mask_ch != "X":
                warnings_list.append(
                    f"Original X at pos {i} for {taxon} changed to {mask_ch!r}"
                )
            elif orig_ch == "-" and mask_ch == "X":
                warnings_list.append(
                    f"Gap at pos {i} for {taxon} became X — no CDS change"
                )
        new_nt_records.append(SeqRecord(Seq("".join(nt_chars)), id=taxon, description=""))
    nt_output_path.parent.mkdir(parents=True, exist_ok=True)
    SeqIO.write(new_nt_records, str(nt_output_path), "fasta")
    return {"projected_codons": projected, "warnings": warnings_list}


_STANDARD_AA = set("ARNDCQEGHILKMFPSTWYV")


def _verify_taper_outputs(aa_path: Path, nt_path: Path | None) -> bool:
    aa_ok = validate_fasta_output(aa_path, require_aligned=True).ok
    if nt_path is not None:
        return aa_ok and validate_fasta_output(nt_path, require_aligned=True).ok
    return aa_ok
```

- [ ] **Step 2: Verify import still works**

```bash
python -c "from phyloai.pretree.filter import _build_taper_cmd; print('OK')"
```

- [ ] **Step 3: Review checkpoint**

```bash
python -c "from phyloai.pretree.filter import _build_taper_cmd, _project_taper_masks_to_cds; print('OK')"
```

---

### Task 6: Implement `run_taper` main entry point with checkpoint/resume

**Files:**
- Modify: `phyloai/pretree/filter.py`

- [ ] **Step 1: Append `run_taper` function**

```python
# --- TAPER main entry point --- (append to phyloai/pretree/filter.py)


def run_taper(
    msa_dir: Path, output_dir: Path, *,
    seq_type: str = "auto", nt_dir: Path | None = None,
    cutoff: int = _TAPER_CUTOFF_DEFAULT,
    taper_path: Path | None = None, julia_path: Path | None = None,
    threads: int = 4, tool_args: str | None = None,
    resume: bool = False, overwrite: bool = False,
    dry_run: bool = False, quiet: bool = False,
    table_format: str = "csv",
    progress_callback: Callable[[Path], None] | None = None,
) -> dict[str, Any]:
    """Run TAPER error-site masking. Returns result.json-compatible dict."""
    start = time.monotonic()
    if overwrite and resume:
        raise ValueError("--overwrite and --resume are mutually exclusive.")

    env = ToolEnv()
    julia_exe = str(julia_path) if julia_path else str(env.require("julia"))
    taper_script = str(taper_path) if taper_path else str(env.require("correction_multi.jl"))

    delimiter = _table_delimiter(table_format)
    suffix = _table_suffix(table_format)

    msa_map = scan_msa_dir(msa_dir)
    if not msa_map:
        raise ValueError(f"No valid MSA files in {msa_dir}")
    if seq_type == "auto":
        first = list(msa_map.values())[0]
        sample = list(SeqIO.parse(str(first), "fasta"))
        seq_type = resolve_seq_type([str(r.seq) for r in sample])[0]

    nt_map: dict[str, Path] = {}
    if nt_dir is not None:
        nt_map = scan_msa_dir(nt_dir)
        if not nt_map:
            raise ValueError(f"No valid NT MSA files in {nt_dir}")
    is_aa_cds = nt_dir is not None

    params = {
        "msa_dir": str(msa_dir), "nt_dir": str(nt_dir) if nt_dir else None,
        "seq_type": seq_type, "cutoff": cutoff,
        "taper_path": taper_script, "julia_path": julia_exe,
        "threads": threads, "tool_args": tool_args, "table_format": table_format,
    }
    command = f"phyloai pretree filter taper --msa-dir {msa_dir} --seq-type {seq_type} --cutoff {cutoff}"

    ckpt_path = output_dir / "checkpoint.json"
    locus_list = sorted(msa_map.keys())
    input_files = [msa_map[l] for l in locus_list]

    def _output_for(inp: Path) -> Path:
        prefix = "seqs/faa" if is_aa_cds else "seqs"
        return output_dir / prefix / inp.name

    def _nt_output_for(inp: Path) -> Path | None:
        if not is_aa_cds:
            return None
        locus = logical_msa_locus_name(inp)
        if locus in nt_map:
            return output_dir / "seqs" / "fna" / nt_map[locus].name
        return None

    checkpoint: Checkpoint | None = None
    resume_success_results: list[dict] = []
    runner = Runner()

    if dry_run:
        cmds = [" ".join(_build_taper_cmd(inp, _output_for(inp), seq_type, cutoff, julia_exe, taper_script, tool_args)) for inp in input_files]
        return {
            "status": "success", "command": command, "wall_time": time.monotonic() - start,
            "tool_versions": {}, "params": params,
            "key_results": {"n_input": len(input_files)}, "error": None,
            "data": {"dry_run_cmds": cmds, "summary": {"n_input_files": len(input_files)}},
        }

    if resume:
        checkpoint = load_checkpoint(ckpt_path)
        validate_resume_params(checkpoint, params, step="pretree.filter.taper")
        to_run_ids, _ = plan_resume(checkpoint, _verify_taper_outputs)
        for task in checkpoint.tasks:
            if task.task_id in set(to_run_ids):
                continue
            if task.status == "success":
                resume_success_results.append({
                    "locus": task.task_id, "status": "success",
                    "output": task.outputs.get("aa", ""),
                    "nt_output": task.outputs.get("nt"),
                    "new_masked_sites": 0,
                })
        input_files = [Path(task.input) for task in checkpoint.tasks if task.task_id in set(to_run_ids)]
    else:
        _common_output_conflict(output_dir, overwrite)
        checkpoint = build_initial_checkpoint(
            step="pretree.filter.taper", command=command, params=params,
            inputs=input_files, output_for=_output_for, nt_output_for=_nt_output_for,
        )
        save_checkpoint_atomic(checkpoint, ckpt_path)

    to_run_ids = [logical_msa_locus_name(f) for f in input_files]
    if checkpoint and to_run_ids:
        for tid in to_run_ids:
            mark_task(checkpoint, tid, status="running")
        save_checkpoint_atomic(checkpoint, ckpt_path)

    file_results: list[dict] = list(resume_success_results)
    last_flush = time.monotonic()
    interrupted = False

    try:
        if to_run_ids:
            with ProcessPoolExecutor(max_workers=threads) as pool:
                futures = {}
                for inp in input_files:
                    fut = pool.submit(
                        _run_taper_one, input_file=inp, output_file=_output_for(inp),
                        seq_type=seq_type, cutoff=cutoff, julia_exe=julia_exe,
                        taper_script=taper_script, tool_args=tool_args, runner=runner,
                    )
                    futures[fut] = logical_msa_locus_name(inp)
                for fut in as_completed(futures):
                    task_id = futures[fut]
                    try:
                        result = fut.result()
                    except Exception as exc:
                        result = {"locus": task_id, "status": "failed", "reason": str(exc)[:200], "output": ""}
                    if is_aa_cds and result.get("status") == "success":
                        locus = result["locus"]
                        if locus in nt_map:
                            nt_out = _nt_output_for(msa_map[locus])
                            try:
                                proj = _project_taper_masks_to_cds(msa_map[locus], Path(result["output"]), nt_map[locus], nt_out)
                                result["nt_output"] = str(nt_out)
                                result["projected_codons"] = proj["projected_codons"]
                            except Exception as exc:
                                result["status"] = "failed"
                                result["nt_error"] = str(exc)[:200]
                    file_results.append(result)
                    if checkpoint:
                        mark_task(checkpoint, task_id, status=result.get("status", "failed"), reason=result.get("reason"))
                        now = time.monotonic()
                        if now - last_flush >= _CHECKPOINT_FLUSH_INTERVAL:
                            save_checkpoint_atomic(checkpoint, ckpt_path)
                            last_flush = now
                    if progress_callback:
                        progress_callback(Path(result.get("output", "")))
    except KeyboardInterrupt:
        interrupted = True
        if checkpoint:
            checkpoint.status = "interrupted"
            save_checkpoint_atomic(checkpoint, ckpt_path, fsync=True)
        raise

    if checkpoint:
        checkpoint.status = "success" if not interrupted else "interrupted"
        checkpoint.completed_at = None if interrupted else checkpoint.touch()
        save_checkpoint_atomic(checkpoint, ckpt_path, fsync=True)

    retained = [r for r in file_results if r.get("status") == "success"]
    dropped = [r for r in file_results if r.get("status") != "success"]

    if not dry_run:
        _write_csv_table([{"locus": r["locus"]} for r in retained], output_dir / f"retained_loci{suffix}", ["locus"], delimiter)
        _write_csv_table([{"locus": r["locus"], "reason": r.get("reason", "")} for r in dropped], output_dir / f"dropped_loci{suffix}", ["locus", "reason"], delimiter)
        decisions = [{"locus": r.get("locus", ""), "status": r.get("status", ""), "new_masked_sites": r.get("new_masked_sites", 0), "output": r.get("output", "")} for r in file_results]
        _write_csv_table(decisions, output_dir / f"filter_decisions{suffix}", ["locus", "status", "new_masked_sites", "output"], delimiter)

    wall_time = time.monotonic() - start
    payload = {
        "status": "success" if retained else "error",
        "command": command, "wall_time": round(wall_time, 2),
        "tool_versions": {"julia": "unknown", "correction_multi.jl": "1.0.0"},
        "params": params,
        "key_results": {
            "n_input": len(file_results), "n_retained": len(retained), "n_dropped": len(dropped),
            "total_masked_aa_sites": sum(r.get("new_masked_sites", 0) for r in file_results),
        },
        "error": None if retained else "All loci failed TAPER.",
        "data": {
            "retained_loci": [r["locus"] for r in retained],
            "dropped_loci": [r["locus"] for r in dropped],
            "file_results": file_results,
            "retained_msa_stats": _compute_retained_msa_stats(
                [Path(r["output"]) for r in retained if r.get("output")]),
        },
    }
    if not dry_run:
        _write_result_json(payload, output_dir)
        _write_filter_log(output_dir, command, wall_time, payload["tool_versions"], payload["status"] == "success")
    return payload
```

- [ ] **Step 2: Verify import cleanly**

```bash
python -c "from phyloai.pretree.filter import run_taper; print('OK')"
```

- [ ] **Step 3: Review checkpoint**

```bash
python -c "from phyloai.pretree.filter import run_taper; print('OK')"
```

---

### Task 7: Add TAPER CLI subcommand and unit tests

**Files:**
- Modify: `phyloai/cli/commands/pretree.py`
- Create: `tests/pretree/test_filter.py`

- [ ] **Step 1: Add `filter taper` CLI command**

Insert before `pretree.add_command(filter_group)` in `phyloai/cli/commands/pretree.py`:

```python
@filter_group.command(
    "taper",
    help=(
        "Mask error sites in MSAs using TAPER.\n\n"
        "Modes: AA-only, NT-only (--seq-type NT), AA+CDS (--nt-dir). "
        "--cutoff defaults to 3; lower values are more aggressive."
    ),
)
@click.option("--msa-dir", type=click.Path(exists=True, file_okay=False, path_type=Path), required=True, help="Directory of MSA files.")
@click.option("--nt-dir", type=click.Path(exists=True, file_okay=False, path_type=Path), default=None, help="Codon-aligned NT MSAs for AA+CDS mode.")
@click.option("--seq-type", type=click.Choice(["AA", "NT", "auto"]), default="auto", show_default=True)
@click.option("--cutoff", type=click.IntRange(1), default=3, show_default=True, help="TAPER -c cutoff; 1-10 typical, lower=more aggressive.")
@click.option("--taper-path", type=click.Path(exists=True, dir_okay=False, path_type=Path), default=None, help="Path to correction_multi.jl.")
@click.option("--julia-path", type=click.Path(exists=True, dir_okay=False, path_type=Path), default=None, help="Julia executable path.")
@click.option("--tool-args", type=str, default=None, help="Extra TAPER args (not -m,-a,-c,-l).")
@click.option("--threads", "-t", type=int, default=4, show_default=True)
@click.option("--output-dir", "-o", type=click.Path(file_okay=False, path_type=Path), default=Path("runs/pretree/filter/taper"), show_default=True)
@click.option("--table-format", type=click.Choice(["csv", "tsv"]), default="csv", show_default=True)
@click.option("--resume", is_flag=True, default=False)
@click.option("--overwrite", is_flag=True, default=False)
@click.option("--dry-run", is_flag=True, default=False)
@click.option("--quiet", "-q", is_flag=True, default=False)
def filter_taper_command(msa_dir, nt_dir, seq_type, cutoff, taper_path, julia_path, tool_args, threads, output_dir, table_format, resume, overwrite, dry_run, quiet):
    if threads < 1:
        _fail("--threads must be at least 1.", 1)
    if nt_dir is not None and seq_type == "NT":
        _fail("--nt-dir (AA+CDS mode) incompatible with --seq-type NT.", 1)
    try:
        payload = run_taper(msa_dir=msa_dir, output_dir=output_dir, seq_type=seq_type, nt_dir=nt_dir, cutoff=cutoff, taper_path=taper_path, julia_path=julia_path, threads=threads, tool_args=tool_args, resume=resume, overwrite=overwrite, dry_run=dry_run, quiet=quiet, table_format=table_format)
    except (ValueError, FileNotFoundError) as exc:
        _fail(str(exc), 3 if "not found" in str(exc).lower() else 1)
    if dry_run:
        click.echo(f"Dry run: {payload['key_results']['n_input']} loci would be processed.")
        for cmd in payload["data"]["dry_run_cmds"]:
            click.echo(cmd)
        return
    if not quiet:
        console.print(render_filter_summary_table({
            "Input": payload["key_results"]["n_input"],
            "Retained": payload["key_results"]["n_retained"],
            "Dropped": payload["key_results"]["n_dropped"],
            "Masked AA sites": payload["key_results"]["total_masked_aa_sites"],
        }))
        click.echo(f"Masked MSAs saved to {output_dir / 'seqs'}", err=True)
        click.echo(f"Results saved to {output_dir / 'result.json'}", err=True)
    if payload["status"] == "error":
        _fail(payload.get("error", "All loci failed."), 1)
```

- [ ] **Step 2: Write TAPER unit tests**

```python
# tests/pretree/test_filter.py

import pytest
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


class TestBuildTaperCmd:
    def test_aa_default(self):
        from phyloai.pretree.filter import _build_taper_cmd
        cmd = _build_taper_cmd(Path("a.fa"), Path("out.fa"), "AA", 3, "julia", "/t/c.jl", None)
        assert cmd[:3] == ["julia", "/t/c.jl", "-c"]
        assert "3" in cmd

    def test_nt_has_mode_flags(self):
        from phyloai.pretree.filter import _build_taper_cmd
        cmd = _build_taper_cmd(Path("a.fa"), Path("out.fa"), "NT", 3, "julia", "/t/c.jl", None)
        assert "-m" in cmd and "-a" in cmd

    def test_blocks_managed_flags(self):
        from phyloai.pretree.filter import _build_taper_cmd
        with pytest.raises(ValueError, match="managed"):
            _build_taper_cmd(Path("a.fa"), Path("out.fa"), "AA", 3, "julia", "/t/c.jl", "-c 5")


class TestTaperCDSProjection:
    def test_projection(self, tmp_path):
        from phyloai.pretree.filter import _project_taper_masks_to_cds
        aa_original = tmp_path / "aa_original.fa"
        aa_masked = tmp_path / "aa_masked.fa"
        SeqIO.write([SeqRecord(Seq("AA-"), id="t1", description=""), SeqRecord(Seq("AX-"), id="t2", description="")], str(aa_original), "fasta")
        SeqIO.write([SeqRecord(Seq("AX-"), id="t1", description=""), SeqRecord(Seq("AX-"), id="t2", description="")], str(aa_masked), "fasta")
        nt_path = tmp_path / "nt.fa"
        SeqIO.write([SeqRecord(Seq("GCANNN---"), id="t1", description=""), SeqRecord(Seq("GCANNN---"), id="t2", description="")], str(nt_path), "fasta")
        out = tmp_path / "out.fna"
        result = _project_taper_masks_to_cds(aa_original, aa_masked, nt_path, out)
        assert result["projected_codons"] == 1
        assert out.exists()


class TestRetainedMsaStats:
    def test_empty(self):
        from phyloai.pretree.filter import _compute_retained_msa_stats
        s = _compute_retained_msa_stats([])
        assert s["n_msa"] == 0

    def test_one_fasta(self, tmp_path):
        from phyloai.pretree.filter import _compute_retained_msa_stats
        p = tmp_path / "g.fa"
        SeqIO.write([SeqRecord(Seq("ACGT"), id="a", description=""), SeqRecord(Seq("ACGT"), id="b", description="")], str(p), "fasta")
        s = _compute_retained_msa_stats([p])
        assert s["n_msa"] == 1
        assert s["mean_taxa"] == 2.0
        assert s["total_length"] == 4
```

- [ ] **Step 3: Run tests**

```bash
pytest tests/pretree/test_filter.py -v
```

- [ ] **Step 4: Verify CLI**

```bash
python -m phyloai pretree filter taper --help
```

- [ ] **Step 5: Review checkpoint**

```bash
pytest tests/pretree/test_filter.py -v
python -m phyloai pretree filter taper --help
```

---

## Phase 3: TreeShrink Subcommand

### Task 8: Implement `run_treeshrink` and CLI

**Files:**
- Modify: `phyloai/pretree/filter.py`
- Modify: `phyloai/cli/commands/pretree.py`
- Modify: `tests/pretree/test_filter.py`

- [ ] **Step 1: Append `run_treeshrink` to filter.py**

```python
# --- TreeShrink --- (append to phyloai/pretree/filter.py)

_TREESHRINK_MANAGED_FLAGS = {"-i", "-t", "-a", "-q", "-m", "-o", "-O"}


def run_treeshrink(
    tree_dir: Path, output_dir: Path, *,
    msa_dir: Path | None = None, threshold: float = 0.05,
    treeshrink_mode: str = "auto", treeshrink_path: Path | None = None,
    tool_args: str | None = None, keep_work_dir: bool = False,
    overwrite: bool = False, dry_run: bool = False,
    quiet: bool = False, table_format: str = "csv",
) -> dict[str, Any]:
    """Run TreeShrink taxon pruning. Returns result.json-compatible dict."""
    start = time.monotonic()
    env = ToolEnv()
    treeshrink_exe = str(treeshrink_path) if treeshrink_path else str(env.require("run_treeshrink.py"))

    delimiter = _table_delimiter(table_format)
    suffix = _table_suffix(table_format)
    tree_map = scan_tree_dir(tree_dir)
    if not tree_map:
        raise ValueError(f"No valid tree files in {tree_dir}")
    msa_map: dict[str, Path] = scan_msa_dir(msa_dir) if msa_dir else {}
    pairing = pair_msa_and_tree_maps(msa_map, list(tree_map.values()))

    params = {"tree_dir": str(tree_dir), "msa_dir": str(msa_dir) if msa_dir else None,
              "threshold": threshold, "treeshrink_mode": treeshrink_mode, "table_format": table_format}
    command = f"phyloai pretree filter treeshrink --tree-dir {tree_dir} --threshold {threshold}"

    if dry_run:
        # Build cmd for display only; no work_dir needed
        work_dir_display = output_dir / "work" if keep_work_dir else Path("/tmp/treeshrink_tmp")
        cmd_display = [treeshrink_exe, "-i", str(work_dir_display / "input"), "-t", "input.tree", "-q", str(threshold)]
        if msa_dir:
            cmd_display.extend(["-a", "input.fasta"])
        if treeshrink_mode != "auto":
            cmd_display.extend(["-m", treeshrink_mode])
        return {"status": "success", "command": command, "wall_time": 0, "tool_versions": {},
                "params": params, "key_results": {"n_input": len(pairing.paired)}, "error": None,
                "data": {"dry_run_cmd": " ".join(cmd_display), "summary": {"n_input_files": len(pairing.paired)}}}

    # Conflict check BEFORE creating any work directories
    _common_output_conflict(output_dir, overwrite)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Now safe to create work layout
    work_dir = output_dir / "work" if keep_work_dir else Path(tempfile.mkdtemp(prefix="treeshrink_"))
    input_dir = work_dir / "input"

    for locus, (msa_path, tree_path) in pairing.paired.items():
        if tree_path is None:
            continue
        gene_dir = input_dir / locus
        gene_dir.mkdir(parents=True, exist_ok=True)
        shutil.copy2(tree_path, gene_dir / "input.tree")
        if msa_path is not None:
            shutil.copy2(msa_path, gene_dir / "input.fasta")

    cmd = [treeshrink_exe, "-i", str(input_dir), "-t", "input.tree", "-q", str(threshold)]
    if msa_dir:
        cmd.extend(["-a", "input.fasta"])
    if treeshrink_mode != "auto":
        cmd.extend(["-m", treeshrink_mode])
    if tool_args:
        extra = shlex.split(tool_args)
        for flag in _TREESHRINK_MANAGED_FLAGS:
            if flag in extra:
                raise ValueError(f"Flag {flag!r} is managed by PhyloAI; remove from --tool-args.")
        cmd.extend(extra)

    runner = Runner()
    result = runner.run(cmd, tool_name="run_treeshrink.py")

    trees_out = output_dir / "trees"
    seqs_out = output_dir / "seqs"
    trees_out.mkdir(parents=True, exist_ok=True)

    file_results, retained, dropped, modified_loci, removed_taxa = [], [], [], [], []
    for locus, (msa_path, tree_path) in pairing.paired.items():
        if tree_path is None:
            dropped.append({"locus": locus, "reason": "no tree input"})
            continue
        src_tree = input_dir / locus / "output.tree"
        if src_tree.exists():
            dst_tree = trees_out / f"{locus}.tre"
            shutil.copy2(src_tree, dst_tree)
            entry = {"locus": locus, "status": "success", "output_tree": str(dst_tree)}
            try:
                from Bio import Phylo
                in_tree = Phylo.read(str(tree_path), "newick")
                out_tree = Phylo.read(str(src_tree), "newick")
                in_taxa = {c.name for c in in_tree.get_terminals()}
                out_taxa = {c.name for c in out_tree.get_terminals()}
                removed = in_taxa - out_taxa
                if removed:
                    modified_loci.append({"locus": locus, "removed_count": len(removed)})
                    for t in sorted(removed):
                        removed_taxa.append({"locus": locus, "taxon": t})
            except Exception:
                pass
            if msa_path:
                src_fa = input_dir / locus / "output.fasta"
                if src_fa.exists():
                    seqs_out.mkdir(parents=True, exist_ok=True)
                    shutil.copy2(src_fa, seqs_out / f"{locus}.fa")
                    entry["output_msa"] = str(seqs_out / f"{locus}.fa")
            retained.append(entry)
            file_results.append(entry)
        else:
            dropped.append({"locus": locus, "reason": "output missing"})
            file_results.append({"locus": locus, "status": "failed", "reason": "output missing"})

    if not dry_run:
        _write_csv_table([{"locus": r["locus"]} for r in retained], output_dir / f"retained_loci{suffix}", ["locus"], delimiter)
        _write_csv_table(dropped, output_dir / f"dropped_loci{suffix}", ["locus", "reason"], delimiter)
        _write_csv_table(modified_loci, output_dir / f"modified_loci{suffix}", ["locus", "removed_count"], delimiter)
        _write_csv_table(removed_taxa, output_dir / f"removed_taxa{suffix}", ["locus", "taxon"], delimiter)
        decisions = [{"locus": r.get("locus", ""), "status": r.get("status", "failed"), "removed_count": sum(1 for t in removed_taxa if t["locus"] == r.get("locus", ""))} for r in file_results]
        _write_csv_table(decisions, output_dir / f"filter_decisions{suffix}", ["locus", "status", "removed_count"], delimiter)

    if not keep_work_dir:
        shutil.rmtree(work_dir, ignore_errors=True)

    msa_stats = _compute_retained_msa_stats(list(seqs_out.glob("*.fa"))) if msa_dir else {}
    wall_time = time.monotonic() - start
    payload = {"status": "success" if retained else "error", "command": command, "wall_time": round(wall_time, 2),
               "tool_versions": {"run_treeshrink.py": "unknown"}, "params": params,
               "key_results": {"n_input": len(pairing.paired), "n_retained": len(retained), "n_modified": len(modified_loci), "n_dropped": len(dropped), "n_removed_taxa_total": len(removed_taxa)},
               "error": None if retained else "All loci failed.", "data": {"retained_loci": [r["locus"] for r in retained], "modified_loci": modified_loci, "dropped_loci": dropped, "removed_taxa": removed_taxa, "retained_msa_stats": msa_stats}}
    _write_result_json(payload, output_dir)
    _write_filter_log(output_dir, command, wall_time, payload["tool_versions"], payload["status"] == "success")
    return payload
```

- [ ] **Step 2: Add `filter treeshrink` CLI command**

Insert before `pretree.add_command(filter_group)`:

```python
@filter_group.command("treeshrink", help="Prune outlier taxa from gene trees using TreeShrink.")
@click.option("--tree-dir", type=click.Path(exists=True, file_okay=False, path_type=Path), required=True)
@click.option("--msa-dir", type=click.Path(exists=True, file_okay=False, path_type=Path), default=None, help="Optional MSA directory to also shrink.")
@click.option("--threshold", type=click.FloatRange(0.0), default=0.05, show_default=True)
@click.option("--treeshrink-mode", type=click.Choice(["auto", "per-gene", "all-genes", "per-species"]), default="auto", show_default=True)
@click.option("--treeshrink-path", type=click.Path(exists=True, dir_okay=False, path_type=Path), default=None)
@click.option("--tool-args", type=str, default=None, help="Extra TreeShrink args (not -i,-t,-a,-q,-m,-o,-O).")
@click.option("--keep-work-dir", is_flag=True, default=False)
@click.option("--output-dir", "-o", type=click.Path(file_okay=False, path_type=Path), default=Path("runs/pretree/filter/treeshrink"), show_default=True)
@click.option("--table-format", type=click.Choice(["csv", "tsv"]), default="csv", show_default=True)
@click.option("--overwrite", is_flag=True, default=False)
@click.option("--dry-run", is_flag=True, default=False)
@click.option("--quiet", "-q", is_flag=True, default=False)
def filter_treeshrink_command(tree_dir, msa_dir, threshold, treeshrink_mode, treeshrink_path, tool_args, keep_work_dir, output_dir, table_format, overwrite, dry_run, quiet):
    try:
        payload = run_treeshrink(tree_dir=tree_dir, output_dir=output_dir, msa_dir=msa_dir, threshold=threshold, treeshrink_mode=treeshrink_mode, treeshrink_path=treeshrink_path, tool_args=tool_args, keep_work_dir=keep_work_dir, overwrite=overwrite, dry_run=dry_run, quiet=quiet, table_format=table_format)
    except (ValueError, FileNotFoundError) as exc:
        _fail(str(exc), 3 if "not found" in str(exc).lower() else 1)
    if dry_run:
        click.echo(f"Dry run: would process {payload['key_results']['n_input']} loci.")
        click.echo(payload["data"]["dry_run_cmd"])
        return
    if not quiet:
        console.print(render_filter_summary_table({"Input": payload["key_results"]["n_input"], "Retained": payload["key_results"]["n_retained"], "Modified": payload["key_results"]["n_modified"], "Dropped": payload["key_results"]["n_dropped"], "Taxa removed": payload["key_results"]["n_removed_taxa_total"]}))
        click.echo(f"Shrunk trees saved to {output_dir / 'trees'}", err=True)
        click.echo(f"Results saved to {output_dir / 'result.json'}", err=True)
    if payload["status"] == "error":
        _fail(payload.get("error", "All loci failed."), 1)
```

- [ ] **Step 3: Add TreeShrink unit test**

```python
# tests/pretree/test_filter.py — append

class TestTreeshrink:
    def test_empty_tree_dir_raises(self, tmp_path):
        from phyloai.pretree.filter import run_treeshrink
        (tmp_path / "trees").mkdir()
        with pytest.raises(ValueError, match="No valid tree files"):
            run_treeshrink(tree_dir=tmp_path / "trees", output_dir=tmp_path / "out", dry_run=True)
```

- [ ] **Step 4: Review checkpoint**

```bash
pytest tests/pretree/test_filter.py -v
python -m phyloai pretree filter treeshrink --help
```

---

## Phase 4: Metrics Rule Filtering

### Task 9: Implement condition parser and `run_metrics_filter` + CLI

**Files:**
- Modify: `phyloai/pretree/filter.py`
- Modify: `phyloai/cli/commands/pretree.py`
- Modify: `tests/pretree/test_filter.py`

- [ ] **Step 1: Append condition parser and `run_metrics_filter` to filter.py**

```python
# --- Metrics rule filtering --- (append to phyloai/pretree/filter.py)

import re as _re

_OP_PATTERN = _re.compile(r"^([^><=!]+?)(>=|<=|!=|==|>|<)(.+)$")
_NUMERIC_OPS = {">=": lambda a, b: a >= b, "<=": lambda a, b: a <= b, ">": lambda a, b: a > b, "<": lambda a, b: a < b, "==": lambda a, b: a == b, "!=": lambda a, b: a != b}


class FilterCondition:
    def __init__(self, col: str, op: str, value: float | str):
        self.col = col
        self.op = op
        self.value = value

    def __str__(self) -> str:
        value_repr = self.value if isinstance(self.value, float) else repr(self.value)
        return f"{self.col}{self.op}{value_repr}"

    def evaluate(self, row: dict) -> bool:
        raw = row.get(self.col, "")
        if raw in (None, "", "NA"):
            return False
        if isinstance(self.value, str):
            return _NUMERIC_OPS[self.op](str(raw).strip(), self.value)
        try:
            return _NUMERIC_OPS[self.op](float(raw), self.value)
        except (ValueError, TypeError):
            return False


def parse_keep_conditions(keep: str, known_columns: set[str]) -> list[FilterCondition]:
    conditions = []
    for part in keep.split(","):
        part = part.strip()
        if not part:
            continue
        m = _OP_PATTERN.match(part)
        if not m:
            raise ValueError(f"Malformed condition: {part!r}. Expected form: col>=val")
        col, op, val_str = m.group(1).strip(), m.group(2), m.group(3).strip()
        if col not in known_columns:
            raise ValueError(f"Unknown column {col!r} in --keep. Known: {sorted(known_columns)}")
        try:
            val: float | str = float(val_str)
        except ValueError:
            val = val_str.strip("\"'")
        conditions.append(FilterCondition(col, op, val))
    return conditions


def _apply_metric_filters(
    rows: list[dict], conditions: list[FilterCondition], loci_column: str = "loci"
) -> tuple[list[dict], list[dict], dict[str, int]]:
    retained, dropped = [], []
    failure_counts = {str(c): 0 for c in conditions}
    for row in rows:
        failed_conditions = [c for c in conditions if not c.evaluate(row)]
        if not failed_conditions:
            retained.append(row)
        else:
            failures = [str(c) for c in failed_conditions]
            for failure in failures:
                failure_counts[failure] += 1
            dropped.append({**row, "_filter_reason": "FAIL: " + ", ".join(failures)})
    return retained, dropped, failure_counts


def _detect_input_delimiter(path: Path, input_format: str) -> str:
    """Detect CSV/TSV delimiter; fail on ambiguity per MAIN §9.8."""
    if input_format == "csv":
        return ","
    if input_format == "tsv":
        return "\t"
    with open(path, newline="") as fh:
        sample = fh.read(4096)
    tabs = sample.count("\t")
    commas = sample.count(",")
    if tabs == 0 and commas == 0:
        raise ValueError(f"Cannot detect delimiter in {path}: no tabs or commas found.")
    if tabs > 0 and commas > 0:
        # Require at least 2:1 ratio for confident detection
        if max(tabs, commas) < 2 * min(tabs, commas):
            raise ValueError(
                f"Ambiguous delimiter in {path}: {tabs} tabs, {commas} commas. "
                "Use --input-format csv|tsv to specify explicitly."
            )
    return "\t" if tabs > commas else ","


def run_metrics_filter(
    table_path: Path, output_dir: Path, *, keep: str,
    input_format: str = "auto", loci_column: str = "loci",
    msa_dir: Path | None = None, tree_dir: Path | None = None,
    copy: bool = False, overwrite: bool = False,
    dry_run: bool = False, quiet: bool = False,
    table_format: str = "csv",
) -> dict[str, Any]:
    """Filter loci by metric conditions. Returns result.json-compatible dict."""
    start = time.monotonic()
    if copy and not msa_dir and not tree_dir:
        raise ValueError("--copy requires at least one of --msa-dir or --tree-dir.")
    delimiter_in = _detect_input_delimiter(table_path, input_format)
    delimiter_out = _table_delimiter(table_format)
    suffix = _table_suffix(table_format)
    rows = []
    with open(table_path, newline="") as fh:
        for row in csv.DictReader(fh, delimiter=delimiter_in):
            rows.append(row)
    if not rows:
        raise ValueError(f"No data rows in {table_path}")
    columns = list(rows[0].keys())
    conditions = parse_keep_conditions(keep, set(columns))
    retained, dropped, failure_counts = _apply_metric_filters(rows, conditions, loci_column)
    params = {"table": str(table_path), "keep": keep, "input_format": input_format, "loci_column": loci_column, "copy": copy, "table_format": table_format}
    command = f"phyloai pretree filter metrics --table {table_path} --keep {keep!r}"
    if dry_run:
        return {"status": "success", "command": command, "wall_time": 0, "tool_versions": {}, "params": params, "key_results": {"n_total": len(rows), "n_retained": len(retained), "n_dropped": len(dropped)}, "error": None, "data": {"condition_failure_counts": failure_counts}}
    _common_output_conflict(output_dir, overwrite)
    output_dir.mkdir(parents=True, exist_ok=True)
    _write_csv_table([{loci_column: r[loci_column]} for r in retained], output_dir / f"retained_loci{suffix}", [loci_column], delimiter_out)
    _write_csv_table([{loci_column: d[loci_column], "reason": d.get("_filter_reason", "")} for d in dropped], output_dir / f"dropped_loci{suffix}", [loci_column, "reason"], delimiter_out)
    decisions = [{loci_column: r[loci_column], "status": "retained", "reason": ""} for r in retained] + [{loci_column: d[loci_column], "status": "dropped", "reason": d.get("_filter_reason", "")} for d in dropped]
    _write_csv_table(decisions, output_dir / f"filter_decisions{suffix}", [loci_column, "status", "reason"], delimiter_out)
    copied_msa, copied_tree = 0, 0
    msa_map = scan_msa_dir(msa_dir) if msa_dir else {}
    tree_map = scan_tree_dir(tree_dir) if tree_dir else {}
    retained_set = {r[loci_column] for r in retained}

    if copy:
        if msa_map:
            (output_dir / "seqs").mkdir(parents=True, exist_ok=True)
            for locus in retained_set:
                if locus in msa_map:
                    shutil.copy2(msa_map[locus], output_dir / "seqs" / msa_map[locus].name)
                    copied_msa += 1
        if tree_map:
            (output_dir / "trees").mkdir(parents=True, exist_ok=True)
            for locus in retained_set:
                if locus in tree_map:
                    shutil.copy2(tree_map[locus], output_dir / "trees" / tree_map[locus].name)
                    copied_tree += 1

    msa_stats = _compute_retained_msa_stats(
        [msa_map[l] for l in retained_set if l in msa_map]
    ) if msa_map else {}
    wall_time = time.monotonic() - start
    payload = {"status": "success", "command": command, "wall_time": round(wall_time, 2), "tool_versions": {}, "params": params, "key_results": {"n_total": len(rows), "n_retained": len(retained), "n_dropped": len(dropped), "condition_failure_counts": failure_counts}, "error": None, "data": {"copied_msa": copied_msa, "copied_tree": copied_tree, "retained_msa_stats": msa_stats, "condition_failure_counts": failure_counts}}
    _write_result_json(payload, output_dir)
    _write_filter_log(output_dir, command, wall_time, {}, True)
    return payload
```

- [ ] **Step 2: Add `filter metrics` CLI subcommand**

Insert before `pretree.add_command(filter_group)`:

```python
@filter_group.command("metrics", help="Filter loci by metric conditions (AND-only rules).")
@click.option("--table", "table_path", type=click.Path(exists=True, dir_okay=False, path_type=Path), required=True, help="Metrics CSV/TSV table.")
@click.option("--keep", type=str, required=True, help="Comma-separated conditions, e.g. 'dvmc>=0,dvmc<=0.3'.")
@click.option("--input-format", type=click.Choice(["csv", "tsv", "auto"]), default="auto", show_default=True)
@click.option("--loci-column", type=str, default="loci", show_default=True)
@click.option("--msa-dir", type=click.Path(exists=True, file_okay=False, path_type=Path), default=None)
@click.option("--tree-dir", type=click.Path(exists=True, file_okay=False, path_type=Path), default=None)
@click.option("--copy", is_flag=True, default=False, help="Copy retained MSAs/trees to output dir.")
@click.option("--output-dir", "-o", type=click.Path(file_okay=False, path_type=Path), default=Path("runs/pretree/filter/metrics"), show_default=True)
@click.option("--table-format", type=click.Choice(["csv", "tsv"]), default="csv", show_default=True)
@click.option("--overwrite", is_flag=True, default=False)
@click.option("--dry-run", is_flag=True, default=False)
@click.option("--quiet", "-q", is_flag=True, default=False)
def filter_metrics_command(table_path, keep, input_format, loci_column, msa_dir, tree_dir, copy, output_dir, table_format, overwrite, dry_run, quiet):
    try:
        payload = run_metrics_filter(table_path=table_path, output_dir=output_dir, keep=keep, input_format=input_format, loci_column=loci_column, msa_dir=msa_dir, tree_dir=tree_dir, copy=copy, overwrite=overwrite, dry_run=dry_run, quiet=quiet, table_format=table_format)
    except (ValueError, FileNotFoundError) as exc:
        _fail(str(exc), 1)
    if dry_run:
        click.echo(f"Dry run: {payload['key_results']['n_total']} loci -> {payload['key_results']['n_retained']} retained, {payload['key_results']['n_dropped']} dropped")
        return
    if not quiet:
        console.print(render_filter_summary_table({"Total": payload["key_results"]["n_total"], "Retained": payload["key_results"]["n_retained"], "Dropped": payload["key_results"]["n_dropped"]}))
        click.echo(f"Decision tables saved to {output_dir}", err=True)
    if payload["status"] == "error":
        _fail(payload.get("error", "Filtering failed."), 1)
```

- [ ] **Step 3: Add metrics filter unit tests**

```python
# tests/pretree/test_filter.py — append

class TestFilterCondition:
    def test_numeric_gte(self):
        from phyloai.pretree.filter import FilterCondition
        c = FilterCondition("dvmc", ">=", 0.3)
        assert c.evaluate({"dvmc": "0.5"})
        assert not c.evaluate({"dvmc": "0.1"})
        assert not c.evaluate({"dvmc": ""})

    def test_string_eq(self):
        from phyloai.pretree.filter import FilterCondition
        c = FilterCondition("DataType", "==", "AA")
        assert c.evaluate({"DataType": "AA"})
        assert not c.evaluate({"DataType": "NT"})


class TestParseKeepConditions:
    def test_simple(self):
        from phyloai.pretree.filter import parse_keep_conditions
        conds = parse_keep_conditions("dvmc>=0,dvmc<=0.3,average_BS>=0.8", {"dvmc", "average_BS", "num_sites"})
        assert len(conds) == 3

    def test_unknown_column_raises(self):
        from phyloai.pretree.filter import parse_keep_conditions
        with pytest.raises(ValueError, match="Unknown column"):
            parse_keep_conditions("badcol>=0", {"dvmc"})

    def test_malformed_raises(self):
        from phyloai.pretree.filter import parse_keep_conditions
        with pytest.raises(ValueError, match="Malformed"):
            parse_keep_conditions("not_valid", {"dvmc"})
```

- [ ] **Step 4: Review checkpoint**

```bash
pytest tests/pretree/test_filter.py -v
python -m phyloai pretree filter metrics --help
```

---

## Phase 5: Cluster Subcommand

### Task 10: Implement cluster core (feature selection, PCA, clustering, k-selection)

**Files:**
- Modify: `phyloai/pretree/filter.py`

- [ ] **Step 1: Append cluster core functions**

```python
# --- Cluster-based filtering --- (append to phyloai/pretree/filter.py)

import numpy as np


def _select_features(rows: list[dict], columns: list[str], metrics: str | None, exclude_regex: list[str], loci_column: str) -> list[str]:
    """Select numeric feature columns from metrics table, excluding constant/empty columns."""
    exclude_patterns = [_re.compile(p) for p in (exclude_regex or [])]
    all_numeric = []
    for col in columns:
        if col == loci_column or col == "DataType":
            continue
        if any(pat.search(col) for pat in exclude_patterns):
            continue
        vals = set()
        for row in rows:
            v = row.get(col, "")
            if v not in (None, "", "NA"):
                try:
                    vals.add(float(v))
                except (ValueError, TypeError):
                    pass
        if len(vals) > 1:  # exclude constant columns
            all_numeric.append(col)
    if metrics and metrics != "all":
        requested = [m.strip() for m in metrics.split(",")]
        return [c for c in requested if c in all_numeric]
    return all_numeric


def _extract_feature_matrix(rows: list[dict], features: list[str], loci_column: str) -> tuple[np.ndarray, list[str], list[dict]]:
    data, labels, valid_rows = [], [], []
    for row in rows:
        vals = []
        all_valid = True
        for f in features:
            v = row.get(f, "")
            if v in (None, "", "NA"):
                all_valid = False
                break
            try:
                vals.append(float(v))
            except (ValueError, TypeError):
                all_valid = False
                break
        if all_valid:
            data.append(vals)
            labels.append(row.get(loci_column, ""))
            valid_rows.append(row)
    return np.array(data, dtype=float), labels, valid_rows


def _scale_features(matrix: np.ndarray) -> np.ndarray:
    from sklearn.preprocessing import StandardScaler
    return StandardScaler().fit_transform(matrix)


def _reduce_pca(scaled: np.ndarray) -> np.ndarray:
    from sklearn.decomposition import PCA
    return PCA(n_components=min(3, scaled.shape[1])).fit_transform(scaled)


def _hierarchical_clustering(
    reduced: np.ndarray, n_clusters: int | None, max_clusters: int | None,
    linkage: str, distance: str, n_loci: int,
) -> tuple[int, np.ndarray, list[dict]]:
    """Cluster and optionally auto-select k via multi-metric voting.

    Per design §8.5: silhouette (higher better), Calinski-Harabasz (higher better),
    Davies-Bouldin (lower better).  Each votes for best k.  Ties broken by
    higher silhouette, then smaller k.  Returns (selected_k, labels, selection_rows)
    where selection_rows is for ``cluster_selection.csv``.
    """
    from sklearn.cluster import AgglomerativeClustering
    from sklearn.metrics import silhouette_score, calinski_harabasz_score, davies_bouldin_score

    if max_clusters is None:
        max_clusters = min(30, max(6, int(np.ceil(np.sqrt(n_loci) / 3))))

    if n_clusters is not None:
        cl = AgglomerativeClustering(n_clusters=n_clusters, metric=distance, linkage=linkage)
        labels = cl.fit_predict(reduced)
        return n_clusters, labels, []

    k_min, k_max = 2, min(max_clusters, n_loci - 1)
    if k_max < k_min:
        return 1, np.zeros(n_loci, dtype=int), []

    results: list[dict] = []
    for k in range(k_min, k_max + 1):
        cl = AgglomerativeClustering(n_clusters=k, metric=distance, linkage=linkage)
        lbs = cl.fit_predict(reduced)
        sil = silhouette_score(reduced, lbs, metric=distance)
        ch = calinski_harabasz_score(reduced, lbs)
        db = davies_bouldin_score(reduced, lbs)
        results.append({"k": k, "silhouette": sil, "calinski_harabasz": ch, "davies_bouldin": db})

    # Vote: rank each k for each metric (higher=better for sil & ch; lower=better for db)
    ks = [r["k"] for r in results]
    sil_ordered = sorted(ks, key=lambda k: -next(r["silhouette"] for r in results if r["k"] == k))
    ch_ordered = sorted(ks, key=lambda k: -next(r["calinski_harabasz"] for r in results if r["k"] == k))
    db_ordered = sorted(ks, key=lambda k: next(r["davies_bouldin"] for r in results if r["k"] == k))

    rank_sums: dict[int, int] = {}
    for k in ks:
        rank_sums[k] = sil_ordered.index(k) + ch_ordered.index(k) + db_ordered.index(k)

    # Tie-break: lower rank-sum, then higher silhouette, then smaller k
    best_k = min(ks, key=lambda k: (rank_sums[k],
                                     -next(r["silhouette"] for r in results if r["k"] == k), k))

    cl = AgglomerativeClustering(n_clusters=best_k, metric=distance, linkage=linkage)
    labels = cl.fit_predict(reduced)
    return best_k, labels, results


def _select_best_umap_replicate(
    scaled: np.ndarray,
    n_replicates: int,
    base_random_state: int,
    n_neighbors: int,
    min_dist: float,
    n_clusters: int | None,
    max_clusters: int | None,
    linkage: str,
    distance: str,
) -> tuple[np.ndarray, int, int, list[dict], list[dict]]:
    from umap import UMAP

    replicate_rows: list[dict] = []
    best: tuple | None = None
    for replicate_index in range(n_replicates):
        random_state = base_random_state + replicate_index
        reducer = UMAP(n_components=3, n_neighbors=n_neighbors, min_dist=min_dist, random_state=random_state)
        reduced = reducer.fit_transform(scaled)
        selected_k, labels, selection_rows = _hierarchical_clustering(
            reduced, n_clusters, max_clusters, linkage, distance, len(scaled)
        )
        if selection_rows:
            selected_row = next(row for row in selection_rows if row["k"] == selected_k)
            silhouette = selected_row["silhouette"]
            ch = selected_row["calinski_harabasz"]
            db = selected_row["davies_bouldin"]
        else:
            silhouette = float("nan")
            ch = float("nan")
            db = float("nan")
        replicate_rows.append({
            "replicate": replicate_index,
            "random_state": random_state,
            "selected_k": selected_k,
            "silhouette": silhouette,
            "calinski_harabasz": ch,
            "davies_bouldin": db,
        })
        best = (reduced, selected_k, replicate_index, selection_rows) if best is None else best

    valid_metric_rows = [row for row in replicate_rows if row["silhouette"] == row["silhouette"]]
    if valid_metric_rows:
        sil_order = sorted(valid_metric_rows, key=lambda row: (-row["silhouette"], row["replicate"]))
        ch_order = sorted(valid_metric_rows, key=lambda row: (-row["calinski_harabasz"], row["replicate"]))
        db_order = sorted(valid_metric_rows, key=lambda row: (row["davies_bouldin"], row["replicate"]))
        for row in valid_metric_rows:
            row["rank_sum"] = sil_order.index(row) + ch_order.index(row) + db_order.index(row)
        best_row = min(valid_metric_rows, key=lambda row: (row["rank_sum"], -row["silhouette"], row["davies_bouldin"], row["replicate"]))
        best_index = best_row["replicate"]
    else:
        for row in replicate_rows:
            row["rank_sum"] = "NA"
        best_index = 0

    reducer = UMAP(n_components=3, n_neighbors=n_neighbors, min_dist=min_dist, random_state=base_random_state + best_index)
    best_reduced = reducer.fit_transform(scaled)
    best_k, best_labels, best_selection_rows = _hierarchical_clustering(
        best_reduced, n_clusters, max_clusters, linkage, distance, len(scaled)
    )
    return best_reduced, best_k, best_index, replicate_rows, best_selection_rows
```

- [ ] **Step 2: Review checkpoint**

```bash
python -c "from phyloai.pretree.filter import _select_features, _hierarchical_clustering; print('OK')"
```

---

### Task 11: Implement cluster diagnostics and outlier removal

**Files:**
- Modify: `phyloai/pretree/filter.py`

- [ ] **Step 1: Append diagnostic plots and outlier removal**

```python
# --- Cluster diagnostics --- (append after Task 10 code)


def _generate_cluster_plots(reduced: np.ndarray, labels: np.ndarray, output_dir: Path, n_clusters: int) -> list[str]:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D  # noqa: F401
    paths = []
    cmap = plt.get_cmap("tab20")
    # 2D scatter
    fig, ax = plt.subplots(figsize=(10, 8))
    for c in range(n_clusters):
        mask = labels == c
        ax.scatter(reduced[mask, 0], reduced[mask, 1], color=cmap(c % 20), label=f"Cluster {c}", alpha=0.7, s=30)
    ax.set_xlabel("Dim 1"); ax.set_ylabel("Dim 2"); ax.set_title(f"Cluster Scatter (2D) — k={n_clusters}")
    ax.legend(bbox_to_anchor=(1.05, 1), loc="upper left", fontsize=8)
    fig.tight_layout()
    pdf_2d = output_dir / "cluster_2d.pdf"
    fig.savefig(pdf_2d, dpi=150, bbox_inches="tight")
    plt.close(fig)
    paths.append(str(pdf_2d))
    # 3D scatter
    if reduced.shape[1] >= 3:
        fig = plt.figure(figsize=(10, 8))
        ax = fig.add_subplot(111, projection="3d")
        for c in range(n_clusters):
            mask = labels == c
            ax.scatter(reduced[mask, 0], reduced[mask, 1], reduced[mask, 2],
                       color=cmap(c % 20), label=f"Cluster {c}", alpha=0.7, s=20)
        ax.set_xlabel("Dim 1"); ax.set_ylabel("Dim 2"); ax.set_zlabel("Dim 3")
        ax.set_title(f"Cluster Scatter (3D) — k={n_clusters}")
        fig.tight_layout()
        pdf_3d = output_dir / "cluster_3d.pdf"
        fig.savefig(pdf_3d, dpi=150, bbox_inches="tight")
        plt.close(fig)
        paths.append(str(pdf_3d))
    return paths


def _generate_cluster_metric_means(
    valid_rows: list[dict], labels: np.ndarray, features: list[str],
    loci_column: str, valid_loci: list[str], output_dir: Path,
    table_format: str,
) -> str:
    """Write cluster_metric_means.<fmt> and cluster_metric_heatmap.pdf."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from sklearn.preprocessing import StandardScaler

    n_clusters = len(set(labels))
    cluster_means = np.zeros((n_clusters, len(features)))
    for c in range(n_clusters):
        indices = [i for i, lb in enumerate(labels) if lb == c]
        for j, f in enumerate(features):
            vals = [float(valid_rows[i].get(f, 0)) for i in indices if valid_rows[i].get(f, "") not in ("", "NA")]
            cluster_means[c, j] = np.mean(vals) if vals else 0.0

    # CSV
    means_rows = []
    for c in range(n_clusters):
        entry = {"cluster": c, "n_loci": int((labels == c).sum())}
        for j, f in enumerate(features):
            entry[f] = round(float(cluster_means[c, j]), 6)
        means_rows.append(entry)
    delimiter_out = _table_delimiter(table_format)
    suffix = _table_suffix(table_format)
    _write_csv_table(means_rows, output_dir / f"cluster_metric_means{suffix}",
                     ["cluster", "n_loci"] + features, delimiter_out)

    # Heatmap PDF
    scaler = StandardScaler()
    heat_data = scaler.fit_transform(cluster_means)
    fig, ax = plt.subplots(figsize=(max(8, len(features) * 0.5), max(3, n_clusters * 0.5)))
    im = ax.imshow(heat_data, aspect="auto", cmap="RdBu_r", interpolation="nearest")
    ax.set_xticks(range(len(features)))
    ax.set_xticklabels(features, rotation=45, ha="right", fontsize=8)
    ax.set_yticks(range(n_clusters))
    ax.set_yticklabels([f"Cluster {c}" for c in range(n_clusters)])
    ax.set_title("Standardized Cluster Metric Means")
    plt.colorbar(im, ax=ax)
    fig.tight_layout()
    pdf = output_dir / "cluster_metric_heatmap.pdf"
    fig.savefig(pdf, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return str(pdf)


def _auto_drop_outlier_clusters(labels: np.ndarray, rows: list[dict], n_loci: int, outlier_metric: str, outlier_direction: str, max_drop_fraction: float) -> tuple[set[int], list[dict]]:
    max_drop = int(np.floor(n_loci * max_drop_fraction))
    if max_drop == 0:
        return set(), []
    n_clusters = len(set(labels))
    cluster_means = {}
    cluster_sizes = {}
    for c in range(n_clusters):
        indices = [i for i, lb in enumerate(labels) if lb == c]
        cluster_sizes[c] = len(indices)
        vals = []
        for i in indices:
            v = rows[i].get(outlier_metric, "")
            try:
                vals.append(float(v))
            except (ValueError, TypeError):
                pass
        cluster_means[c] = np.mean(vals) if vals else 0.0
    reverse = outlier_direction == "high"
    sorted_clusters = sorted(cluster_means.items(), key=lambda x: x[1], reverse=reverse)
    drop = set()
    dropped_count = 0
    for c, _ in sorted_clusters:
        if dropped_count + cluster_sizes[c] <= max_drop:
            drop.add(c)
            dropped_count += cluster_sizes[c]
        else:
            break
    return drop, []


def _generate_cluster_metric_boxplots(
    valid_rows: list[dict], labels: np.ndarray, features: list[str],
    output_dir: Path, n_clusters: int, plot_metrics_per_page: str,
) -> list[str]:
    """Write cluster_metric_boxplots_*.pdf pages required by design §8.8."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    if plot_metrics_per_page != "auto":
        per_page = int(plot_metrics_per_page)
    else:
        per_page = 12 if n_clusters <= 6 else 6 if n_clusters <= 12 else 4 if n_clusters <= 20 else 2
    pdf_paths: list[str] = []
    for page_start in range(0, len(features), per_page):
        page_features = features[page_start:page_start + per_page]
        fig, axes = plt.subplots(len(page_features), 1, figsize=(10, max(3, 2.5 * len(page_features))))
        if len(page_features) == 1:
            axes = [axes]
        for ax, feature in zip(axes, page_features):
            grouped = []
            for c in range(n_clusters):
                values = []
                for i, lb in enumerate(labels):
                    if lb != c:
                        continue
                    raw = valid_rows[i].get(feature, "")
                    if raw in ("", "NA"):
                        continue
                    try:
                        values.append(float(raw))
                    except (TypeError, ValueError):
                        continue
                grouped.append(values)
            ax.boxplot(grouped, labels=[f"C{c}" for c in range(n_clusters)])
            ax.set_title(feature)
        fig.tight_layout()
        pdf_path = output_dir / f"cluster_metric_boxplots_{(page_start // per_page) + 1:03d}.pdf"
        fig.savefig(pdf_path, dpi=150, bbox_inches="tight")
        plt.close(fig)
        pdf_paths.append(str(pdf_path))
    return pdf_paths


def _write_outlier_diagnostics(
    valid_rows: list[dict], labels: np.ndarray, drop_clusters: set[int],
    features: list[str], output_dir: Path, table_format: str,
) -> list[str]:
    from scipy.stats import mannwhitneyu
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    delimiter_out = _table_delimiter(table_format)
    suffix = _table_suffix(table_format)
    outlier_flags = [label in drop_clusters for label in labels]

    comparison_rows = []
    wilcoxon_rows = []
    for feature in features:
        normal_vals = []
        outlier_vals = []
        for row, is_outlier in zip(valid_rows, outlier_flags):
            raw = row.get(feature, "")
            if raw in ("", "NA"):
                continue
            try:
                value = float(raw)
            except (TypeError, ValueError):
                continue
            if is_outlier:
                outlier_vals.append(value)
            else:
                normal_vals.append(value)
        comparison_rows.append({
            "metric": feature,
            "normal_mean": round(float(np.mean(normal_vals)), 6) if normal_vals else "NA",
            "normal_median": round(float(np.median(normal_vals)), 6) if normal_vals else "NA",
            "normal_std": round(float(np.std(normal_vals)), 6) if normal_vals else "NA",
            "normal_n": len(normal_vals),
            "outlier_mean": round(float(np.mean(outlier_vals)), 6) if outlier_vals else "NA",
            "outlier_median": round(float(np.median(outlier_vals)), 6) if outlier_vals else "NA",
            "outlier_std": round(float(np.std(outlier_vals)), 6) if outlier_vals else "NA",
            "outlier_n": len(outlier_vals),
        })
        if normal_vals and outlier_vals:
            stat = mannwhitneyu(normal_vals, outlier_vals, alternative="two-sided")
            direction = "outlier_higher" if np.mean(outlier_vals) > np.mean(normal_vals) else "outlier_lower"
            wilcoxon_rows.append({"metric": feature, "u_statistic": round(float(stat.statistic), 6), "p_value": round(float(stat.pvalue), 6), "direction": direction})
        else:
            wilcoxon_rows.append({"metric": feature, "u_statistic": "NA", "p_value": "NA", "direction": "insufficient_data"})

    _write_csv_table(comparison_rows, output_dir / f"outlier_comparison{suffix}", ["metric", "normal_mean", "normal_median", "normal_std", "normal_n", "outlier_mean", "outlier_median", "outlier_std", "outlier_n"], delimiter_out)
    _write_csv_table(wilcoxon_rows, output_dir / f"outlier_wilcoxon{suffix}", ["metric", "u_statistic", "p_value", "direction"], delimiter_out)

    pdf_paths: list[str] = []
    per_page = 12 if len(drop_clusters) <= 6 else 6
    for page_start in range(0, len(features), per_page):
        page_features = features[page_start:page_start + per_page]
        fig, axes = plt.subplots(len(page_features), 1, figsize=(10, max(3, 2.5 * len(page_features))))
        if len(page_features) == 1:
            axes = [axes]
        for ax, feature in zip(axes, page_features):
            normal_vals = []
            outlier_vals = []
            for row, is_outlier in zip(valid_rows, outlier_flags):
                raw = row.get(feature, "")
                if raw in ("", "NA"):
                    continue
                try:
                    value = float(raw)
                except (TypeError, ValueError):
                    continue
                (outlier_vals if is_outlier else normal_vals).append(value)
            ax.boxplot([normal_vals, outlier_vals], labels=["normal", "outlier"])
            ax.set_title(feature)
        fig.tight_layout()
        pdf_path = output_dir / f"outlier_comparison_boxplots_{(page_start // per_page) + 1:03d}.pdf"
        fig.savefig(pdf_path, dpi=150, bbox_inches="tight")
        plt.close(fig)
        pdf_paths.append(str(pdf_path))
    return pdf_paths
```

- [ ] **Step 2: Review checkpoint**

```bash
python -c "from phyloai.pretree.filter import _generate_cluster_plots, _generate_cluster_metric_means, _generate_cluster_metric_boxplots, _write_outlier_diagnostics, _select_best_umap_replicate; print('OK')"
```

---

### Task 12: Implement `run_cluster_filter` and CLI

**Files:**
- Modify: `phyloai/pretree/filter.py`
- Modify: `phyloai/cli/commands/pretree.py`
- Modify: `tests/pretree/test_filter.py`

- [ ] **Step 1: Append `run_cluster_filter`**

```python
# --- Cluster main entry point --- (append after Task 11 code)


def run_cluster_filter(
    table_path: Path, output_dir: Path, *,
    input_format: str = "auto", metrics: str | None = None,
    exclude_regex: list[str] | None = None, reduction: str = "pca",
    n_clusters: int | None = None, max_clusters: int | None = None,
    cluster_linkage: str = "ward", cluster_distance: str = "euclidean",
    drop_outlier_clusters: str = "none", outlier_metric: str = "average_BS",
    outlier_direction: str = "low", max_drop_fraction: float = 0.2,
    plot_metrics_per_page: str = "auto", plot_label_angle: float = 45.0,
    umap_n_neighbors: int = 15, umap_min_dist: float = 0.001,
    umap_replicates: int = 1, umap_random_state: int = 0,
    msa_dir: Path | None = None, tree_dir: Path | None = None,
    copy: bool = False, overwrite: bool = False,
    dry_run: bool = False, quiet: bool = False,
    table_format: str = "csv",
) -> dict[str, Any]:
    """Run cluster-based exploration and optional outlier filtering."""
    start = time.monotonic()
    if reduction == "umap":
        try:
            import umap  # noqa: F401
        except ImportError:
            raise ImportError("umap-learn required for --reduction umap. pip install umap-learn")
    if cluster_linkage == "ward" and cluster_distance != "euclidean":
        raise ValueError("Ward linkage requires Euclidean distance.")
    delimiter_in = _detect_input_delimiter(table_path, input_format)
    delimiter_out = _table_delimiter(table_format)
    suffix = _table_suffix(table_format)
    loci_column = "loci"
    rows = []
    with open(table_path, newline="") as fh:
        for row in csv.DictReader(fh, delimiter=delimiter_in):
            rows.append(row)
    if not rows:
        raise ValueError(f"No data rows in {table_path}")
    columns = list(rows[0].keys())
    features = _select_features(rows, columns, metrics, exclude_regex or [], loci_column)
    if len(features) < 2:
        raise ValueError(f"Need >=2 features; found {len(features)}.")
    params = {"table": str(table_path), "metrics": metrics, "reduction": reduction, "n_clusters": n_clusters, "max_clusters": max_clusters, "cluster_linkage": cluster_linkage, "cluster_distance": cluster_distance, "drop_outlier_clusters": drop_outlier_clusters, "table_format": table_format}
    command = f"phyloai pretree filter cluster --table {table_path} --reduction {reduction}"
    if dry_run:
        k_range = [n_clusters, n_clusters] if n_clusters is not None else [2, min(max_clusters or min(30, max(6, int(np.ceil(np.sqrt(len(rows)) / 3)))), max(2, len(rows) - 1))]
        return {"status": "success", "command": command, "wall_time": 0, "tool_versions": {}, "params": params, "key_results": {"n_loci": len(rows), "n_features": len(features)}, "error": None, "data": {"features": features, "reduction": reduction, "k_range": k_range, "drop_outlier_clusters": drop_outlier_clusters, "copy": copy}}
    _common_output_conflict(output_dir, overwrite)
    output_dir.mkdir(parents=True, exist_ok=True)
    matrix, valid_loci, valid_rows = _extract_feature_matrix(rows, features, loci_column)
    scaled = _scale_features(matrix)
    if reduction == "pca":
        reduced = _reduce_pca(scaled)
        selected_replicate = None
        umap_replicate_rows: list[dict] = []
        selected_k, labels, selection_rows = _hierarchical_clustering(reduced, n_clusters, max_clusters, cluster_linkage, cluster_distance, len(valid_loci))
    else:
        reduced, selected_k, selected_replicate, umap_replicate_rows, selection_rows = _select_best_umap_replicate(
            scaled,
            umap_replicates,
            umap_random_state,
            umap_n_neighbors,
            umap_min_dist,
            n_clusters,
            max_clusters,
            cluster_linkage,
            cluster_distance,
        )
        _, labels, _ = _hierarchical_clustering(reduced, selected_k, max_clusters, cluster_linkage, cluster_distance, len(valid_loci))
    # --- Write all diagnostic outputs (per design §8.8) ---
    coord_names = ["PC1", "PC2", "PC3"] if reduction == "pca" else ["UMAP1", "UMAP2", "UMAP3"]

    # features_used.csv
    _write_csv_table(
        [{"column": f, "included": True} for f in features],
        output_dir / f"features_used{suffix}", ["column", "included"], delimiter_out,
    )
    # reduction.csv (one row per locus with coordinates + cluster)
    red_rows = []
    for i, locus in enumerate(valid_loci):
        row = {loci_column: locus}
        for j, cname in enumerate(coord_names):
            if j < reduced.shape[1]:
                row[cname] = round(float(reduced[i, j]), 6)
        row["cluster"] = int(labels[i])
        red_rows.append(row)
    _write_csv_table(red_rows, output_dir / f"reduction{suffix}",
                     [loci_column] + coord_names[:reduced.shape[1]] + ["cluster"], delimiter_out)
    # cluster_selection.csv (validation metrics per k)
    if selection_rows:
        _write_csv_table(selection_rows, output_dir / f"cluster_selection{suffix}",
                         ["k", "silhouette", "calinski_harabasz", "davies_bouldin"], delimiter_out)
    if umap_replicate_rows:
        _write_csv_table(umap_replicate_rows, output_dir / f"umap_replicates{suffix}", ["replicate", "random_state", "selected_k", "silhouette", "calinski_harabasz", "davies_bouldin", "rank_sum"], delimiter_out)
    # clusters.csv (per-locus cluster assignments)
    cluster_assign_rows = [{loci_column: valid_loci[i], "cluster": int(labels[i])} for i in range(len(valid_loci))]
    _write_csv_table(cluster_assign_rows, output_dir / f"clusters{suffix}", [loci_column, "cluster"], delimiter_out)
    # cluster_summary.csv
    _write_csv_table([{"cluster": c, "n_loci": int((labels == c).sum())} for c in range(selected_k)],
                     output_dir / f"cluster_summary{suffix}", ["cluster", "n_loci"], delimiter_out)
    # cluster_loci/ per-cluster files
    cluster_loci_dir = output_dir / "cluster_loci"
    cluster_loci_dir.mkdir(parents=True, exist_ok=True)
    for c in range(selected_k):
        mask = labels == c
        loci_in = [valid_loci[i] for i in range(len(valid_loci)) if mask[i]]
        _write_csv_table([{loci_column: l} for l in loci_in], cluster_loci_dir / f"cluster_{c}{suffix}", [loci_column], delimiter_out)
    # cluster_metric_means.csv + heatmap
    means_path = _generate_cluster_metric_means(valid_rows, labels, features, loci_column, valid_loci, output_dir, table_format)
    boxplot_paths = _generate_cluster_metric_boxplots(valid_rows, labels, features, output_dir, selected_k, plot_metrics_per_page)
    # Diagnostic plots: 2D scatter + 3D scatter
    plot_paths = _generate_cluster_plots(reduced, labels, output_dir, selected_k)
    # --- Outlier removal (opt-in) ---
    drop_clusters: set[int] = set()
    if drop_outlier_clusters == "auto":
        drop_clusters, _ = _auto_drop_outlier_clusters(labels, valid_rows, len(valid_rows), outlier_metric, outlier_direction, max_drop_fraction)
        if drop_clusters:
            retained_set = [valid_loci[i] for i in range(len(valid_loci)) if labels[i] not in drop_clusters]
            dropped_set = [valid_loci[i] for i in range(len(valid_loci)) if labels[i] in drop_clusters]
            _write_csv_table([{loci_column: l} for l in retained_set], output_dir / f"retained_loci{suffix}", [loci_column], delimiter_out)
            _write_csv_table([{loci_column: l, "reason": f"outlier_cluster"} for l in dropped_set], output_dir / f"dropped_loci{suffix}", [loci_column, "reason"], delimiter_out)
            decisions = [{loci_column: valid_loci[i], "status": "dropped" if labels[i] in drop_clusters else "retained", "cluster": int(labels[i])} for i in range(len(valid_loci))]
            _write_csv_table(decisions, output_dir / f"filter_decisions{suffix}", [loci_column, "status", "cluster"], delimiter_out)
            outlier_plot_paths = _write_outlier_diagnostics(valid_rows, labels, drop_clusters, features, output_dir, table_format)
            if copy:
                retained_locus_names = set(retained_set)
                if msa_dir:
                    msa_map = scan_msa_dir(msa_dir)
                    (output_dir / "seqs").mkdir(parents=True, exist_ok=True)
                    for locus in retained_locus_names:
                        if locus in msa_map:
                            shutil.copy2(msa_map[locus], output_dir / "seqs" / msa_map[locus].name)
                if tree_dir:
                    tree_map = scan_tree_dir(tree_dir)
                    (output_dir / "trees").mkdir(parents=True, exist_ok=True)
                    for locus in retained_locus_names:
                        if locus in tree_map:
                            shutil.copy2(tree_map[locus], output_dir / "trees" / tree_map[locus].name)
        elif copy:
            # Per design §8.7: if no loci dropped, copy is a no-op with warning
            import sys
            print("[WARN] No outlier clusters dropped (all within max_drop_fraction). Copy skipped.",
                  file=sys.stderr)
    wall_time = time.monotonic() - start
    payload = {"status": "success", "command": command, "wall_time": round(wall_time, 2), "tool_versions": {}, "params": params, "key_results": {"n_loci": len(rows), "n_features": len(features), "n_clusters": selected_k, "reduction": reduction, "selected_umap_replicate": selected_replicate, "n_dropped": sum((labels == c).sum() for c in drop_clusters) if drop_clusters else 0}, "error": None, "data": {"features": features, "cluster_sizes": {c: int((labels == c).sum()) for c in range(selected_k)}, "drop_clusters": sorted(drop_clusters), "plot_paths": plot_paths + boxplot_paths + [means_path] + (outlier_plot_paths if drop_clusters else []), "umap_replicates": umap_replicate_rows}}
    _write_result_json(payload, output_dir)
    _write_filter_log(output_dir, command, wall_time, {}, True)
    return payload
```

- [ ] **Step 2: Add `filter cluster` CLI subcommand**

Insert before `pretree.add_command(filter_group)`:

```python
@filter_group.command("cluster", help="Group loci by metric profiles using PCA/UMAP + hierarchical clustering.")
@click.option("--table", "table_path", type=click.Path(exists=True, dir_okay=False, path_type=Path), required=True, help="Metrics CSV/TSV table.")
@click.option("--input-format", type=click.Choice(["csv", "tsv", "auto"]), default="auto", show_default=True)
@click.option("--metrics", type=str, default=None, help="Comma-separated columns; default=all numeric excl. freq/sd. 'all'=all.")
@click.option("--exclude-regex", type=str, multiple=True, default=None, help="Exclude columns matching regex (repeatable).")
@click.option("--reduction", type=click.Choice(["pca", "umap"]), default="pca", show_default=True)
@click.option("--n-clusters", type=int, default=None, help="Fixed cluster count; auto-selected if omitted.")
@click.option("--max-clusters", type=int, default=None)
@click.option("--cluster-linkage", type=click.Choice(["ward", "average", "complete", "single"]), default="ward", show_default=True)
@click.option("--cluster-distance", type=click.Choice(["euclidean", "cosine", "manhattan"]), default="euclidean", show_default=True)
@click.option("--drop-outlier-clusters", type=click.Choice(["none", "auto"]), default="none", show_default=True)
@click.option("--outlier-metric", type=str, default="average_BS", show_default=True)
@click.option("--outlier-direction", type=click.Choice(["low", "high"]), default="low", show_default=True)
@click.option("--max-drop-fraction", type=click.FloatRange(0.0, 1.0), default=0.2, show_default=True)
@click.option("--plot-metrics-per-page", type=str, default="auto", show_default=True)
@click.option("--plot-label-angle", type=float, default=45.0, show_default=True)
@click.option("--umap-n-neighbors", type=int, default=15)
@click.option("--umap-min-dist", type=float, default=0.001)
@click.option("--umap-replicates", type=int, default=1)
@click.option("--umap-random-state", type=int, default=0)
@click.option("--msa-dir", type=click.Path(exists=True, file_okay=False, path_type=Path), default=None)
@click.option("--tree-dir", type=click.Path(exists=True, file_okay=False, path_type=Path), default=None)
@click.option("--copy", is_flag=True, default=False)
@click.option("--output-dir", "-o", type=click.Path(file_okay=False, path_type=Path), default=Path("runs/pretree/filter/cluster"), show_default=True)
@click.option("--table-format", type=click.Choice(["csv", "tsv"]), default="csv", show_default=True)
@click.option("--overwrite", is_flag=True, default=False)
@click.option("--dry-run", is_flag=True, default=False)
@click.option("--quiet", "-q", is_flag=True, default=False)
def filter_cluster_command(table_path, input_format, metrics, exclude_regex, reduction, n_clusters, max_clusters, cluster_linkage, cluster_distance, drop_outlier_clusters, outlier_metric, outlier_direction, max_drop_fraction, plot_metrics_per_page, plot_label_angle, umap_n_neighbors, umap_min_dist, umap_replicates, umap_random_state, msa_dir, tree_dir, copy, output_dir, table_format, overwrite, dry_run, quiet):
    try:
        payload = run_cluster_filter(table_path=table_path, output_dir=output_dir, input_format=input_format, metrics=metrics, exclude_regex=list(exclude_regex) if exclude_regex else None, reduction=reduction, n_clusters=n_clusters, max_clusters=max_clusters, cluster_linkage=cluster_linkage, cluster_distance=cluster_distance, drop_outlier_clusters=drop_outlier_clusters, outlier_metric=outlier_metric, outlier_direction=outlier_direction, max_drop_fraction=max_drop_fraction, plot_metrics_per_page=plot_metrics_per_page, plot_label_angle=plot_label_angle, umap_n_neighbors=umap_n_neighbors, umap_min_dist=umap_min_dist, umap_replicates=umap_replicates, umap_random_state=umap_random_state, msa_dir=msa_dir, tree_dir=tree_dir, copy=copy, overwrite=overwrite, dry_run=dry_run, quiet=quiet, table_format=table_format)
    except (ValueError, FileNotFoundError, ImportError) as exc:
        _fail(str(exc), 1)
    if dry_run:
        click.echo(f"Dry run: {payload['key_results']['n_loci']} loci, {payload['key_results']['n_features']} features")
        return
    if not quiet:
        console.print(render_filter_summary_table({"Loci": payload["key_results"]["n_loci"], "Features": payload["key_results"]["n_features"], "Reduction": payload["key_results"]["reduction"], "Clusters": payload["key_results"]["n_clusters"], "Dropped": payload["key_results"]["n_dropped"]}))
        click.echo(f"Results saved to {output_dir}", err=True)
```

- [ ] **Step 3: Add cluster unit test**

```python
# tests/pretree/test_filter.py — append

class TestClusterFeatureSelection:
    def test_selects_numeric_excludes_loci(self):
        from phyloai.pretree.filter import _select_features
        rows = [{"loci": "g1", "DataType": "AA", "dvmc": "0.1", "bs": "0.9", "name": "abc"}]
        cols = list(rows[0].keys())
        feats = _select_features(rows, cols, None, [], "loci")
        assert "loci" not in feats
        assert "DataType" not in feats
        assert "name" not in feats  # non-numeric excluded
        assert "dvmc" in feats
        assert "bs" in feats
```

- [ ] **Step 4: Review checkpoint**

```bash
pytest tests/pretree/test_filter.py -v
python -m phyloai pretree filter cluster --help
```

---

## Phase 6: Documentation and Final Integration

### Task 13: Write user docs, update README, run final checks

**Files:**
- Create: `docs/commands/pretree-filter.md`
- Modify: `README.md`

- [ ] **Step 1: Create `docs/commands/pretree-filter.md`** following the structure from design spec §11.

- [ ] **Step 2: Update README command index** with four filter subcommand entries.

- [ ] **Step 3: Run full test suite**

```bash
pytest tests/pretree/test_filter.py tests/core/test_file_matching.py -v
```

- [ ] **Step 4: Run lint**

```bash
ruff check phyloai/pretree/filter.py phyloai/cli/commands/pretree.py tests/pretree/test_filter.py
ruff format --check phyloai/pretree/filter.py
```

- [ ] **Step 5: Verify all CLIs show help**

```bash
python -m phyloai pretree filter --help
python -m phyloai pretree filter taper --help
python -m phyloai pretree filter treeshrink --help
python -m phyloai pretree filter metrics --help
python -m phyloai pretree filter cluster --help
```

- [ ] **Step 6: Final review checkpoint**

```bash
pytest tests/pretree/test_filter.py tests/core/test_file_matching.py -v
ruff check phyloai/pretree/filter.py phyloai/cli/commands/pretree.py tests/pretree/test_filter.py
```

---

## File Summary

| File | Action | Tasks |
|------|--------|-------|
| `phyloai/core/file_matching.py` | Modify | Task 1 |
| `pyproject.toml` | Modify | Task 2 |
| `phyloai/pretree/filter.py` | Create | Tasks 3, 5, 6, 8, 9, 10, 11, 12 |
| `phyloai/cli/commands/pretree.py` | Modify | Tasks 4, 7, 8, 9, 12 |
| `tests/pretree/test_filter.py` | Create | Tasks 7, 8, 9, 12 |
| `tests/core/test_file_matching.py` | Modify | Task 1 |
| `docs/commands/pretree-filter.md` | Create | Task 13 |
| `README.md` | Modify | Task 13 |
