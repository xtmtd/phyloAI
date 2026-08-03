# Posttree Simulate Adequacy Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Implement `phyloai posttree simulate adequacy` to evaluate PPA-DIV, PPA-CONV, PPA-VAR, and PPA-COMP from an observed MSA and a directory of simulated MSAs.

**Architecture:** A single pure-Python module reads all supported MSA formats through the existing `FormatConverter`, computes PhyloBayes-compatible statistics, and persists each valid replicate's scalar and per-taxon values in the shared checkpoint schema. The Click command is a thin wrapper. The report, dynamically generated MCP tool, docs, and workflow skill consume the standard `result.json` output without bespoke integrations.

**Tech Stack:** Python 3.10+, BioPython, stdlib (`csv`, `json`, `math`, `statistics`, `concurrent.futures`), Click, Rich, pytest.

## Global Constraints

- Add no dependency; use `FormatConverter`, `detect_seq_type`, and `phyloai.core.checkpoint`.
- Support only FASTA, PHYLIP-relaxed, PHYLIP-PAML, and NEXUS alignment inputs; independently auto-detect the observed MSA and every simulated MSA through shared format logic.
- Treat only `ACDEFGHIKLMNPQRSTVWY` (AA) or `ACGT` (NT) as observed states; all other characters are missing.
- Match PhyloBayes `AllPostPred`: population SD, `div` pp = `P(sim <= obs)`, all other scalar and taxon pp = `P(sim > obs)`.
- Require at least 10 valid simulations; never serialize `NaN` (`null` in JSON and empty CSV cells for undefined pp).
- Reuse `CheckpointTask.outputs` string values: scalars as strings and `taxon_dist_j` as a JSON-encoded string; parse explicitly during aggregation.
- Use exact input validation and conflict policy from the approved design: unique taxa, equal lengths, name-based remapping, resume fingerprints, `--overwrite`/`--resume` mutual exclusion.
- Do not add a single `--input-format`: original and simulated MSAs may intentionally have mixed formats, so forcing one format would be incorrect. This is the approved scoped exception to the parent shared-flag convention.
- Resolve `--seq-type auto` from the original MSA before any simulated MSA is dispatched; use and report only the resolved `AA` or `NT` value thereafter.
- Expose `--table-format csv|tsv` (default `csv`) and use it consistently for all three output tables and their `result.json` file-object labels.
- For a non-dry run, reject a non-empty output directory unless `--overwrite` or `--resume` is set; `--overwrite` removes it, while `--resume` requires its `checkpoint.json`.
- Every successful payload includes a full resolved top-level `command`, `wall_time`, `tool_versions: {}`, `error: null`, and a nested `key_results.statistics` object with `comp.max` and `comp.mean`.
- Do not hand-register MCP tools. `walk_click_tree()` exposes completed Click leaf commands automatically.

---

## File Structure

| File | Responsibility |
|---|---|
| `phyloai/posttree/simulate_adequacy.py` | MSA validation, statistics, parallel replicate processing, checkpoint lifecycle, CSV/result construction, terminal summary. |
| `phyloai/cli/commands/posttree.py` | Replace the adequacy placeholder with Click options and error-result handling. |
| `tests/posttree/test_simulate_adequacy.py` | Core statistics, validation, aggregation, resume, JSON/CSV/TSV, fixture regression, and dry-run tests. |
| `tests/cli/test_posttree_simulate_alisim.py` | Replace the adequacy-stub assertion with command help, dry-run, and input-error assertions. |
| `phyloai/report/collector.py` | Recognize and order `posttree.simulate.adequacy`. |
| `phyloai/report/templates.py` | Produce deterministic adequacy methods text. |
| `tests/report/test_collector.py` | Test adequacy step parsing and ordering. |
| `tests/report/test_templates.py` | Test adequacy methods generation. |
| `tests/mcp/test_schema_gen.py` | Confirm the Click leaf is exposed as the correct generated MCP tool schema. |
| `docs/commands/posttree-simulate-adequacy.md` | English command documentation. |
| `docs/commands/posttree-simulate-adequacy.zh.md` | Chinese command documentation. |
| `README.md`, `README.zh.md` | Add the simulation-to-adequacy workflow and command-table entry. |
| `skills/phyloai-workflow/SKILL.md` | Describe local-only adequacy approval, outputs, recovery, and interpretation. |
| `skills/phyloai-workflow/references/parameter-annotations.md` | Add Chinese annotations for every adequacy parameter. |
| `docs/superpowers/specs/2026-06-07-phyloai-design.md` | Replace the adequacy future-stub references and update its module decision text. |

---

### Task 1: Build And Verify The Statistical Core

**Files:**
- Create: `phyloai/posttree/simulate_adequacy.py`
- Create: `tests/posttree/test_simulate_adequacy.py`

**Interfaces:**
- Consumes: `FormatConverter.read(path)`, `detect_seq_type(sequences)`, `MultipleSeqAlignment`, `CheckpointTask.outputs` string constraints.
- Produces: `_compute_statistics(alignment, seq_type) -> dict[str, Any]`, `_summarize_distribution(values, obs, direction) -> dict[str, float | int | None]`, `_validate_simulated_alignment(...) -> list[str]`, and `run_simulate_adequacy(...) -> dict[str, Any]` for Tasks 2-4.

- [ ] **Step 1: Write failing known-value and pp tests**

```python
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from phyloai.posttree.simulate_adequacy import (
    _compute_statistics,
    _summarize_distribution,
)


def _msa(*rows: tuple[str, str]) -> MultipleSeqAlignment:
    return MultipleSeqAlignment([
        SeqRecord(Seq(sequence), id=taxon, description="")
        for taxon, sequence in rows
    ])


def test_statistics_known_values_and_gap_exclusion() -> None:
    stats = _compute_statistics(
        _msa(("A", "AA-"), ("B", "AC-"), ("C", "CC-"), ("D", "CA-")),
        "AA",
    )

    assert stats["n_informative_sites"] == 2
    assert stats["div"] == 2.0
    assert stats["siteconvprob"] == 0.5
    assert stats["sitecomp"] == 0.0
    assert stats["comp_max"] == 0.5
    assert stats["comp_mean"] == 0.25


def test_pp_directions_match_phylobayes() -> None:
    div = _summarize_distribution([2.0] * 8 + [1.0] * 2, 1.5, "div")
    conv = _summarize_distribution([2.0] * 8 + [1.0] * 2, 1.5, "high")

    assert div["pp"] == 0.2
    assert conv["pp"] == 0.8


def test_zero_sd_uses_json_safe_undefined_pp() -> None:
    summary = _summarize_distribution([1.0] * 10, 1.0, "high")

    assert summary["sd_sim"] == 0.0
    assert summary["z_score"] == 0.0
    assert summary["pp"] is None
```

- [ ] **Step 2: Run the tests to verify failure**

Run: `pytest tests/posttree/test_simulate_adequacy.py -q`

Expected: FAIL during collection because `phyloai.posttree.simulate_adequacy` does not exist.

- [ ] **Step 3: Implement the minimal parser-independent statistic helpers**

```python
AA_STATES = "ACDEFGHIKLMNPQRSTVWY"
NT_STATES = "ACGT"
SCALAR_NAMES = ("div", "siteconvprob", "sitecomp", "comp_max", "comp_mean")


def _compute_statistics(alignment: MultipleSeqAlignment, seq_type: str) -> dict[str, Any]:
    states = AA_STATES if seq_type == "AA" else NT_STATES
    state_index = {state: index for index, state in enumerate(states)}
    n_states = len(states)
    names = [record.id for record in alignment]
    sequences = [str(record.seq).upper() for record in alignment]
    n_sites = alignment.get_alignment_length()

    site_freqs: list[list[float]] = []
    diversity_total = 0.0
    squared_freq_total = 0.0
    for site in range(n_sites):
        counts = [0] * n_states
        for sequence in sequences:
            index = state_index.get(sequence[site])
            if index is not None:
                counts[index] += 1
        observed = sum(counts)
        if not observed:
            continue
        freqs = [count / observed for count in counts]
        site_freqs.append(freqs)
        diversity_total += sum(count > 0 for count in counts)
        squared_freq_total += sum(freq * freq for freq in freqs)

    if not site_freqs:
        raise ValueError("alignment has no informative sites")
    n_informative = len(site_freqs)
    means = [sum(freq[k] for freq in site_freqs) / n_informative for k in range(n_states)]
    sitecomp = sum(
        sum(freq[k] * freq[k] for freq in site_freqs) / n_informative - means[k] * means[k]
        for k in range(n_states)
    ) / n_states

    taxon_freqs: dict[str, list[float]] = {}
    for name, sequence in zip(names, sequences):
        counts = [0] * n_states
        for state in sequence:
            index = state_index.get(state)
            if index is not None:
                counts[index] += 1
        total = sum(counts)
        if not total:
            raise ValueError(f"taxon {name!r} has no valid characters")
        taxon_freqs[name] = [count / total for count in counts]
    global_freq = [sum(freq[k] for freq in taxon_freqs.values()) / len(taxon_freqs) for k in range(n_states)]
    taxon_dist = {
        name: sum((freq[k] - global_freq[k]) ** 2 for k in range(n_states))
        for name, freq in taxon_freqs.items()
    }
    return {
        "div": diversity_total / n_informative,
        "siteconvprob": squared_freq_total / n_informative,
        "sitecomp": sitecomp,
        "comp_max": max(taxon_dist.values()),
        "comp_mean": sum(taxon_dist.values()) / len(taxon_dist),
        "taxon_dist_j": taxon_dist,
        "n_informative_sites": n_informative,
    }


def _summarize_distribution(values: list[float], obs: float, direction: str) -> dict[str, float | int | None]:
    if len(values) < 10:
        raise ValueError("at least 10 valid simulated MSAs are required")
    mean_sim = sum(values) / len(values)
    sd_sim = math.sqrt(sum(value * value for value in values) / len(values) - mean_sim * mean_sim)
    quantiles = statistics.quantiles(values, n=40)
    ci_lower, ci_upper = quantiles[0], quantiles[-1]
    if sd_sim == 0:
        return {"mean_sim": mean_sim, "sd_sim": 0.0, "ci_lower": mean_sim, "ci_upper": mean_sim, "z_score": 0.0, "pp": None}
    if direction == "div":
        z_score = (mean_sim - obs) / sd_sim
        pp = sum(value <= obs for value in values) / len(values)
    else:
        z_score = (obs - mean_sim) / sd_sim
        pp = sum(value > obs for value in values) / len(values)
    return {"mean_sim": mean_sim, "sd_sim": sd_sim, "ci_lower": ci_lower, "ci_upper": ci_upper, "z_score": z_score, "pp": pp}
```

- [ ] **Step 4: Add focused formula edge-case tests and make them pass**

```python
import pytest


def test_nt_ambiguity_and_all_missing_sites_are_ignored() -> None:
    stats = _compute_statistics(_msa(("A", "AN-"), ("B", "CG-")), "NT")

    assert stats["n_informative_sites"] == 2
    assert stats["div"] == 1.5


def test_fewer_than_ten_replicates_is_rejected() -> None:
    with pytest.raises(ValueError, match="at least 10"):
        _summarize_distribution([1.0] * 9, 1.0, "high")
```

Run: `pytest tests/posttree/test_simulate_adequacy.py -q`

Expected: PASS.

- [ ] **Step 5: Review Task 1 changes**

Do not commit without explicit user approval. Inspect the Task 1 diff and retain it for the next task.

### Task 2: Add Batch Validation, Checkpoint Resume, CSV, And Result Payload

**Files:**
- Modify: `phyloai/posttree/simulate_adequacy.py`
- Modify: `tests/posttree/test_simulate_adequacy.py`

**Interfaces:**
- Consumes: Task 1 `_compute_statistics()` and `_summarize_distribution()`; `FormatConverter`, `Checkpoint`, `CheckpointTask`, `save_checkpoint_atomic`, `load_checkpoint`, and `validate_resume_params`.
- Produces: `run_simulate_adequacy(original_msa: Path, simulated_dir: Path, seq_type: str = "auto", threads: int = 4, table_format: str = "csv", output_dir: Path = ..., overwrite: bool = False, resume: bool = False, dry_run: bool = False, quiet: bool = False, progress_callback: Callable[[int, int], None] | None = None) -> dict[str, Any]`.

- [ ] **Step 1: Write failing batch, validation, and resume tests**

```python
import csv
import json
from pathlib import Path

import pytest

from phyloai.posttree.simulate_adequacy import run_simulate_adequacy


def _fasta(path: Path, rows: list[tuple[str, str]]) -> None:
    path.write_text("".join(f">{name}\n{sequence}\n" for name, sequence in rows))


def _write_ten_simulations(directory: Path) -> None:
    directory.mkdir()
    for index in range(10):
        _fasta(directory / f"sim{index}.fa", [("A", "AC"), ("B", "CA")])


def test_run_writes_all_tables_and_json_safe_result(tmp_path: Path) -> None:
    original = tmp_path / "original.fa"
    simulations = tmp_path / "simulations"
    _fasta(original, [("A", "AC"), ("B", "CA")])
    _write_ten_simulations(simulations)

    result = run_simulate_adequacy(original_msa=original, simulated_dir=simulations, output_dir=tmp_path / "out", quiet=True)

    assert result["status"] == "success"
    assert result["command"] == (
        f"phyloai posttree simulate adequacy --original-msa {original.resolve()} "
        f"--simulated-dir {simulations.resolve()} --seq-type auto --threads 4 "
        f"-o {(tmp_path / 'out').resolve()} --quiet"
    )
    assert isinstance(result["wall_time"], float)
    assert result["tool_versions"] == {}
    assert result["error"] is None
    assert result["key_results"]["seq_type"] == "NT"
    assert set(result["key_results"]["statistics"]) == {"div", "siteconvprob", "sitecomp", "comp"}
    assert set(result["key_results"]["statistics"]["comp"]) == {"max", "mean"}
    assert result["data"]["cmd"] == []
    assert result["data"]["tool_stderr"] == ""
    assert (tmp_path / "out" / "adequacy_summary.csv").exists()
    assert (tmp_path / "out" / "adequacy_taxon_comp.csv").exists()
    assert (tmp_path / "out" / "per_simulation_stats.csv").exists()
    assert json.loads((tmp_path / "out" / "result.json").read_text())["error"] is None


def test_duplicate_original_taxon_is_a_hard_error(tmp_path: Path) -> None:
    original = tmp_path / "original.fa"
    simulations = tmp_path / "simulations"
    _fasta(original, [("A", "AC"), ("A", "CA")])
    _write_ten_simulations(simulations)

    with pytest.raises(ValueError, match="duplicate taxon"):
        run_simulate_adequacy(original_msa=original, simulated_dir=simulations, output_dir=tmp_path / "out", quiet=True)


def test_bad_simulation_is_skipped_and_valid_taxa_are_remapped(tmp_path: Path) -> None:
    original = tmp_path / "original.fa"
    simulations = tmp_path / "simulations"
    _fasta(original, [("A", "AC"), ("B", "CA")])
    _write_ten_simulations(simulations)
    _fasta(simulations / "duplicate.fa", [("A", "AC"), ("A", "CA")])
    _fasta(simulations / "reordered.fa", [("B", "CA"), ("A", "AC")])

    result = run_simulate_adequacy(original_msa=original, simulated_dir=simulations, output_dir=tmp_path / "out", quiet=True)

    assert result["key_results"]["n_failed"] == 1
    with open(tmp_path / "out" / "adequacy_taxon_comp.csv", newline="") as handle:
        assert {row["taxon"] for row in csv.DictReader(handle)} == {"A", "B"}


def test_resume_uses_saved_taxon_values_without_reprocessing(tmp_path: Path) -> None:
    original = tmp_path / "original.fa"
    simulations = tmp_path / "simulations"
    _fasta(original, [("A", "AC"), ("B", "CA")])
    _write_ten_simulations(simulations)
    output_dir = tmp_path / "out"
    first = run_simulate_adequacy(original_msa=original, simulated_dir=simulations, output_dir=output_dir, quiet=True)
    checkpoint = json.loads((output_dir / "checkpoint.json").read_text())
    assert "taxon_dist_j" in checkpoint["tasks"][0]["outputs"]
    expected_taxon_csv = (output_dir / "adequacy_taxon_comp.csv").read_text()

    resumed = run_simulate_adequacy(original_msa=original, simulated_dir=simulations, output_dir=output_dir, resume=True, quiet=True)

    assert resumed["key_results"]["n_simulations"] == first["key_results"]["n_simulations"]
    assert (output_dir / "adequacy_taxon_comp.csv").read_text() == expected_taxon_csv


def test_nonempty_output_requires_overwrite_or_resume(tmp_path: Path) -> None:
    original = tmp_path / "original.fa"
    simulations = tmp_path / "simulations"
    output_dir = tmp_path / "out"
    _fasta(original, [("A", "AC"), ("B", "CA")])
    _write_ten_simulations(simulations)
    output_dir.mkdir()
    (output_dir / "existing.txt").write_text("keep")

    with pytest.raises(ValueError, match="already exists and is non-empty"):
        run_simulate_adequacy(original_msa=original, simulated_dir=simulations, output_dir=output_dir, quiet=True)


def test_resume_requires_checkpoint(tmp_path: Path) -> None:
    original = tmp_path / "original.fa"
    simulations = tmp_path / "simulations"
    output_dir = tmp_path / "out"
    _fasta(original, [("A", "AC"), ("B", "CA")])
    _write_ten_simulations(simulations)
    output_dir.mkdir()

    with pytest.raises(ValueError, match="No checkpoint found"):
        run_simulate_adequacy(original_msa=original, simulated_dir=simulations, output_dir=output_dir, resume=True, quiet=True)
```

- [ ] **Step 2: Run the tests to verify failure**

Run: `pytest tests/posttree/test_simulate_adequacy.py -q`

Expected: FAIL because `run_simulate_adequacy` is not defined.

- [ ] **Step 3: Implement shared-format input validation and one-file processing**

```python
_FORMAT_CONVERTER = FormatConverter()


def _read_alignment(path: Path) -> MultipleSeqAlignment:
    try:
        alignment = _FORMAT_CONVERTER.read(path)
    except Exception as exc:
        raise ValueError(f"unable to parse alignment file {path}: {exc}") from exc
    if not len(alignment) or not alignment.get_alignment_length():
        raise ValueError(f"alignment file {path} is empty")
    return alignment


def _alignment_ids(alignment: MultipleSeqAlignment, label: str) -> list[str]:
    names = [record.id for record in alignment]
    duplicates = sorted({name for name in names if names.count(name) > 1})
    if duplicates:
        raise ValueError(f"{label} alignment has duplicate taxon IDs: {', '.join(duplicates)}")
    if len({len(record.seq) for record in alignment}) != 1:
        raise ValueError(f"{label} alignment sequences have unequal lengths")
    return names


def _process_simulation(path: Path, original_ids: list[str], original_length: int, seq_type: str) -> dict[str, Any]:
    alignment = _read_alignment(path)
    ids = _alignment_ids(alignment, "simulated")
    if set(ids) != set(original_ids):
        raise ValueError("taxon name mismatch between original and simulated MSAs")
    if alignment.get_alignment_length() != original_length:
        raise ValueError("length mismatch between original and simulated MSAs")
    records = {record.id: record for record in alignment}
    ordered = MultipleSeqAlignment([records[name] for name in original_ids])
    return _compute_statistics(ordered, seq_type)


def _resolved_seq_type(original: MultipleSeqAlignment, requested: str) -> str:
    if requested == "auto":
        return detect_seq_type([str(record.seq) for record in original])
    return requested
```

- [ ] **Step 4: Implement checkpoint persistence, parallel dispatch, aggregation, and writers**

```python
def _fingerprint(path: Path) -> str:
    stat = path.stat()
    return f"{path.resolve()}|{stat.st_size}|{stat.st_mtime_ns}"


def _task_outputs(stats: dict[str, Any]) -> dict[str, str]:
    return {
        **{name: repr(float(stats[name])) for name in SCALAR_NAMES},
        "taxon_dist_j": json.dumps({name: repr(float(value)) for name, value in stats["taxon_dist_j"].items()}, sort_keys=True),
    }


def _task_stats(task: CheckpointTask) -> dict[str, Any]:
    return {
        **{name: float(task.outputs[name]) for name in SCALAR_NAMES},
        "taxon_dist_j": {name: float(value) for name, value in json.loads(task.outputs["taxon_dist_j"]).items()},
    }


def _write_csv(path: Path, fieldnames: list[str], rows: list[dict[str, Any]]) -> None:
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="raise")
        writer.writeheader()
        writer.writerows(rows)
```

Implement `run_simulate_adequacy()` with this control flow:

```python
params = {
    "original_msa": str(original_msa.resolve()),
    "simulated_dir": str(simulated_dir.resolve()),
    "seq_type": seq_type,
}
# Validate `seq_type in {"AA", "NT", "auto"}`. Read and validate the
# original MSA first, then assign `resolved_seq_type = _resolved_seq_type(...)`.
# `params` retains requested `seq_type`; payload `key_results["seq_type"]` and
# all worker calls use `resolved_seq_type`, never the literal value "auto".
# Build a separate full payload params dict containing every resolved argument:
# input paths, requested and detected seq types, threads,
# table_format, output_dir, overwrite, resume, dry_run, and quiet. Construct
# `full_command` with shlex.join and every CLI option. Start timing before
# validation. The smaller `params` above is only the resume-compatibility dict,
# so changing threads/table format/output directory does not reject a valid
# resume.
# Reject overwrite and resume together. On a non-dry run, overwrite removes an
# existing output directory; resume loads its checkpoint. Catch FileNotFoundError
# from load_checkpoint and re-raise ValueError with the same message so the CLI
# emits a clean input error. Otherwise a non-empty output directory raises
# ValueError. Create the output directory only after this lifecycle decision.
# Fresh run creates one CheckpointTask per sorted regular non-empty input.
# Resume validates params, compares `_fingerprint(Path(task.task_id))` with
# `task.input`, and changes stale success tasks to pending before dispatch.
# Submit pending tasks to ProcessPoolExecutor(max_workers=threads).  Every
# completed future changes exactly one task to success (with `_task_outputs`)
# or failed (with `reason`), then calls save_checkpoint_atomic().  Aggregate
# every success task via `_task_stats`, require 10, and create scalar rows with
# `_summarize_distribution`: direction="div" for div and "high" for the other
# four statistics. Build per-taxon rows from each success task's `taxon_dist_j`
# values in original taxon order with the same helper using direction="high".
# For each original taxon, write exactly the columns `taxon`, `obs`,
# `mean_pred`, `sd_pred`, `ci_lower`, `ci_upper`, `z_score`, `pp`; use the
# original taxon order. Zero SD gives `ci_lower=ci_upper=mean_pred`,
# `z_score=0.0`, `pp=None`. Undefined pp remains None in key_results and is
# converted to "" only for either CSV. Write all three CSVs, with a delimiter
# and `.csv`/`.tsv` suffix from table_format. Build `key_results["statistics"]` as:
# {"div": summary, "siteconvprob": summary, "sitecomp": summary,
#  "comp": {"max": comp_max_summary, "mean": comp_mean_summary}}.
# Return/write a payload with top-level `command: full_command`, rounded
# `wall_time`, `tool_versions: {}`, and `error: None`, plus
# `data={"cmd": [], "tool_stderr": "", "warnings": warnings,
# "output_files": file_objects}`.
# When progress_callback is provided, call it once after each pending task
# completes with `(completed_remaining, total_remaining)`; the CLI owns one
# Rich Progress instance and passes a callback only when quiet is false. Its
# resume task total excludes already-valid success tasks.
```

- [ ] **Step 5: Add fingerprint, dry-run, and result-schema tests**

```python
def test_replaced_simulation_is_recomputed_on_resume(tmp_path: Path) -> None:
    original = tmp_path / "original.fa"
    simulations = tmp_path / "simulations"
    _fasta(original, [("A", "AC"), ("B", "CA")])
    _write_ten_simulations(simulations)
    output_dir = tmp_path / "out"
    run_simulate_adequacy(original_msa=original, simulated_dir=simulations, output_dir=output_dir, quiet=True)
    _fasta(simulations / "sim0.fa", [("A", "AA"), ("B", "CC")])

    run_simulate_adequacy(original_msa=original, simulated_dir=simulations, output_dir=output_dir, resume=True, quiet=True)

    checkpoint = json.loads((output_dir / "checkpoint.json").read_text())
    task = next(task for task in checkpoint["tasks"] if task["task_id"].endswith("sim0.fa"))
    assert task["input"].split("|")[1:] == [str((simulations / "sim0.fa").stat().st_size), str((simulations / "sim0.fa").stat().st_mtime_ns)]


def test_dry_run_writes_nothing(tmp_path: Path) -> None:
    original = tmp_path / "original.fa"
    simulations = tmp_path / "simulations"
    _fasta(original, [("A", "AC"), ("B", "CA")])
    _write_ten_simulations(simulations)
    output_dir = tmp_path / "out"

    result = run_simulate_adequacy(original_msa=original, simulated_dir=simulations, output_dir=output_dir, dry_run=True, quiet=True)

    assert result["status"] == "success"
    assert not output_dir.exists()


def test_table_format_tsv_uses_tsv_suffixes_and_delimiters(tmp_path: Path) -> None:
    original = tmp_path / "original.fa"
    simulations = tmp_path / "simulations"
    _fasta(original, [("A", "AC"), ("B", "CA")])
    _write_ten_simulations(simulations)

    result = run_simulate_adequacy(
        original_msa=original, simulated_dir=simulations,
        table_format="tsv", output_dir=tmp_path / "out", quiet=True,
    )

    summary = tmp_path / "out" / "adequacy_summary.tsv"
    assert "\t" in summary.read_text().splitlines()[0]
    assert result["data"]["output_files"]["adequacy_summary"]["path"] == str(summary)


def test_progress_callback_counts_only_pending_resume_tasks(tmp_path: Path) -> None:
    original = tmp_path / "original.fa"
    simulations = tmp_path / "simulations"
    _fasta(original, [("A", "AC"), ("B", "CA")])
    _write_ten_simulations(simulations)
    output_dir = tmp_path / "out"
    run_simulate_adequacy(original_msa=original, simulated_dir=simulations, output_dir=output_dir, quiet=True)
    _fasta(simulations / "sim0.fa", [("A", "AA"), ("B", "CC")])
    updates: list[tuple[int, int]] = []

    run_simulate_adequacy(
        original_msa=original, simulated_dir=simulations, output_dir=output_dir,
        resume=True, quiet=True, progress_callback=lambda done, total: updates.append((done, total)),
    )

    assert updates == [(1, 1)]
```

Add a fixture regression test using `runs/zscore/matrix.XX`,
`runs/zscore/simulated/`, and `runs/zscore/readpb/chain1.ppred`. Parse the PB
text file and compare `obs`, `mean_sim`, `sd_sim`, `z_score`, and `pp` for all
five scalar statistics after a full run in `tmp_path`, using `abs=1e-5`. Do not
hard-code rounded expected values.

Run: `pytest tests/posttree/test_simulate_adequacy.py -q`

Expected: PASS.

- [ ] **Step 6: Review Task 2 changes**

Do not commit without explicit user approval. Run the targeted test command and inspect the diff.

### Task 3: Expose The CLI And Integrate Reports

**Files:**
- Modify: `phyloai/cli/commands/posttree.py:1699-1720`
- Modify: `tests/cli/test_posttree_simulate_alisim.py:26-36`
- Modify: `phyloai/report/collector.py:10-45`
- Modify: `phyloai/report/collector.py:71-84`
- Modify: `phyloai/report/templates.py:1476-1515`
- Modify: `tests/report/test_collector.py`
- Modify: `tests/report/test_templates.py`

**Interfaces:**
- Consumes: `run_simulate_adequacy()` from Task 2, existing `_write_error_result_json()` and `_fail()` in the posttree CLI, report `parse_step_id()` and `generate_all_methods()` registries.
- Produces: the `posttree.simulate.adequacy` Click/MCP leaf and report methods text.

- [ ] **Step 1: Replace stub expectations with failing CLI/report tests**

```python
def test_adequacy_command_has_options_and_dry_run(tmp_path: Path) -> None:
    original = tmp_path / "original.fa"
    simulations = tmp_path / "simulations"
    original.write_text(">A\nAC\n>B\nCA\n")
    simulations.mkdir()
    (simulations / "sim1.fa").write_text(">A\nAC\n>B\nCA\n")

    help_result = CliRunner().invoke(cli, ["posttree", "simulate", "adequacy", "--help"])
    assert help_result.exit_code == 0
    assert "--original-msa" in help_result.output
    assert "--simulated-dir" in help_result.output

    result = CliRunner().invoke(cli, [
        "posttree", "simulate", "adequacy", "--original-msa", str(original),
        "--simulated-dir", str(simulations), "--dry-run", "--quiet",
    ])
    assert result.exit_code == 0


def test_adequacy_step_id_is_recognized() -> None:
    assert parse_step_id("phyloai posttree simulate adequacy --original-msa real.fa --simulated-dir sims") == "posttree.simulate.adequacy"
    assert "posttree.simulate.adequacy" in STEP_ORDER


def test_adequacy_methods_text() -> None:
    text = generate_all_methods(
        "posttree.simulate.adequacy",
        params={"seq_type": "AA"},
        key_results={"n_simulations": 100, "n_taxa": 6, "n_sites": 235},
        tool_versions={},
    )
    assert "100" in text
    assert "PPA-DIV" in text
    assert "PPA-COMP" in text
```

- [ ] **Step 2: Run targeted tests to verify failure**

Run: `pytest tests/cli/test_posttree_simulate_alisim.py tests/report/test_collector.py tests/report/test_templates.py -q`

Expected: FAIL because adequacy remains a no-options stub and report parsing returns `posttree.simulate`.

- [ ] **Step 3: Replace the Click stub with the minimal wrapper**

```python
@simulate.command("adequacy")
@click.option("--original-msa", type=click.Path(path_type=Path), required=True,
              help="Observed MSA (FASTA/PHYLIP-relaxed/PHYLIP-PAML/NEXUS).")
@click.option("--simulated-dir", type=click.Path(path_type=Path), required=True,
              help="Directory of simulated MSAs in supported alignment formats.")
@click.option("--seq-type", type=click.Choice(["AA", "NT", "auto"], case_sensitive=False), default="auto", show_default=True)
@click.option("--threads", type=click.IntRange(1), default=4, show_default=True)
@click.option("--table-format", type=click.Choice(["csv", "tsv"]), default="csv", show_default=True)
@click.option("-o", "--output-dir", type=click.Path(path_type=Path), default=Path("runs/posttree/simulate/adequacy"), show_default=True)
@click.option("--overwrite", is_flag=True, default=False)
@click.option("--resume", is_flag=True, default=False)
@click.option("--dry-run", is_flag=True, default=False)
@click.option("-q", "--quiet", is_flag=True, default=False)
def simulate_adequacy_command(original_msa: Path, simulated_dir: Path, seq_type: str, threads: int, table_format: str, output_dir: Path, overwrite: bool, resume: bool, dry_run: bool, quiet: bool) -> None:
    """Assess model adequacy from an observed MSA and simulated replicates."""
    from phyloai.posttree.simulate_adequacy import run_simulate_adequacy

    err_parts = [
        "phyloai", "posttree", "simulate", "adequacy",
        "--original-msa", str(original_msa.resolve()),
        "--simulated-dir", str(simulated_dir.resolve()), "--seq-type", seq_type,
        "--threads", str(threads), "--table-format", table_format,
        "-o", str(output_dir.resolve()),
    ]
    if overwrite:
        err_parts.append("--overwrite")
    if resume:
        err_parts.append("--resume")
    if dry_run:
        err_parts.append("--dry-run")
    if quiet:
        err_parts.append("--quiet")
    err_cmd = shlex.join(err_parts)
    try:
        run_simulate_adequacy(original_msa=original_msa, simulated_dir=simulated_dir, seq_type=seq_type, threads=threads, table_format=table_format, output_dir=output_dir, overwrite=overwrite, resume=resume, dry_run=dry_run, quiet=quiet)
    except ValueError as exc:
        _write_error_result_json(output_dir.resolve(), err_cmd, str(exc), "input")
        _fail(str(exc), exit_code=1)
```

- [ ] **Step 4: Register the report step and methods generator**

```python
# collector.py
STEP_ORDER.append("posttree.simulate.adequacy")
# In parse_step_id's _THIRD_LEVEL mapping:
"simulate": {"alisim", "adequacy"},

# templates.py
def generate_methods_posttree_simulate_adequacy(params: dict[str, Any], key_results: dict[str, Any], tool_versions: dict[str, Any]) -> str:
    source = Path(params["simulated_dir"]).name if params.get("simulated_dir") else "simulated MSA directory"
    return (
        f"Model adequacy was assessed in pure Python by comparing four summary statistics "
        f"from the observed {key_results.get('n_taxa', '?')}-taxon "
        f"{key_results.get('seq_type', params.get('seq_type', '?'))} alignment "
        f"({key_results.get('n_sites', '?')} sites) against {key_results.get('n_simulations', 0)} "
        f"simulated replicates from {source}. Mean diversity per site (PPA-DIV), mean squared "
        "empirical state frequency (PPA-CONV), mean variance of site-specific frequencies "
        "(PPA-VAR), and maximum/mean squared compositional deviation across taxa (PPA-COMP) "
        "were calculated. For each statistic, the null distribution was summarized using its mean, "
        "population SD, and empirical 95% interval (p2.5-p97.5); observed values were assessed "
        "using z-scores and posterior predictive p-values. Values with |z| > 2 or pp < 0.05 "
        "were treated as potential model inadequacy."
    )

METHODS_GENERATORS["posttree.simulate.adequacy"] = generate_methods_posttree_simulate_adequacy
```

- [ ] **Step 5: Verify Click-derived MCP exposure without manual registration**

```python
from phyloai.cli.main import cli
from phyloai.mcp.schema_gen import build_mcp_tool, walk_click_tree


def test_adequacy_mcp_tool_is_generated_from_click() -> None:
    tool = next(build_mcp_tool(item) for item in walk_click_tree(cli) if item["tool_name"] == "posttree_simulate_adequacy")

    assert tool["inputSchema"]["required"] == ["original_msa", "simulated_dir"]
    assert "threads" in tool["inputSchema"]["properties"]
```

Run: `pytest tests/posttree/test_simulate_adequacy.py tests/cli/test_posttree_simulate_alisim.py tests/report/test_collector.py tests/report/test_templates.py tests/mcp/test_schema_gen.py -q`

Expected: PASS.

- [ ] **Step 6: Review Task 3 changes**

Do not commit without explicit user approval. Run the targeted integration tests and inspect the diff.

### Task 4: Publish Documentation And Workflow Guidance

**Files:**
- Create: `docs/commands/posttree-simulate-adequacy.md`
- Create: `docs/commands/posttree-simulate-adequacy.zh.md`
- Modify: `README.md:91-99`
- Modify: `README.md:169-177`
- Modify: `README.zh.md` corresponding simulation workflow and command table
- Modify: `skills/phyloai-workflow/SKILL.md:32-35`
- Modify: `skills/phyloai-workflow/SKILL.md:154-190`
- Modify: `skills/phyloai-workflow/references/parameter-annotations.md`
- Modify: `docs/superpowers/specs/2026-06-07-phyloai-design.md:69-79`
- Modify: `docs/superpowers/specs/2026-06-07-phyloai-design.md:139-144`
- Modify: `docs/superpowers/specs/2026-06-07-phyloai-design.md:271-276`
- Modify: `docs/superpowers/specs/2026-06-07-phyloai-design.md:549`

**Interfaces:**
- Consumes: final CLI interface and output schema from Tasks 2-3.
- Produces: user-facing workflow instructions and accurate parent-design status.

- [ ] **Step 1: Write command documentation with the accepted workflow and limitations**

```markdown
# phyloai posttree simulate adequacy

[English](posttree-simulate-adequacy.md) | [中文](posttree-simulate-adequacy.zh.md)

## Purpose

Compares PPA-DIV, PPA-CONV, PPA-VAR, and PPA-COMP statistics from an observed MSA with an empirical null distribution from simulated MSAs. The command is local-only and requires no external executable.

## Usage

```bash
phyloai posttree simulate adequacy \
  --original-msa matrix.fa \
  --simulated-dir runs/sim/MSAs \
  --threads 4 \
  --table-format csv \
  --output-dir runs/adequacy
```

## Outputs

- `adequacy_summary.csv`: observed values, simulated mean/SD, empirical 95% interval, z-score, and pp for five scalar values.
- `adequacy_taxon_comp.csv`: per-taxon PPA-COMP null-distribution comparison.
- `per_simulation_stats.csv`: scalar values from every valid simulated MSA.
- `checkpoint.json`: resumable per-MSA statistics, including per-taxon values.
- `result.json`: machine-readable command result.

## Notes

- Simulations must have the same unique taxon set and alignment length as the observed MSA; invalid simulated files are skipped and reported.
- At least 10 valid simulations are required.
- The observed MSA and every simulated MSA are independently auto-detected. Mixed supported formats are allowed: FASTA, PHYLIP-relaxed, PHYLIP-PAML, and NEXUS.
- Use `--table-format tsv` when tab-delimited output is required; all three tables use the selected delimiter and suffix.
- Run `alisim transfergaps` before adequacy when the original MSA has substantial missing data, because AliSim output is gap-free.
```

- [ ] **Step 2: Update README, parent design, and Skill guidance**

```markdown
# README workflow addition after transfergaps
phyloai posttree simulate adequacy --original-msa markers/concat.aa.fa --simulated-dir runs/transfer -o runs/adequacy
```

```markdown
# Skill rule addition
- `simulate adequacy` is local-only: review all command parameters and get explicit approval, but `doctor` is not required. It independently auto-detects each observed/simulated MSA in the supported formats, writes CSV or TSV through `--table-format`, resumes from a checkpoint, and reports PPA-DIV/CONV/VAR/COMP results. Advise `alisim transfergaps` first when original missing data are substantial. Interpret low pp (<0.05) or |z| > 2 as potential inadequacy; `div` pp measures P(sim <= obs), while the other checks measure P(sim > obs).
```

Add parameter-annotation headings for `--original-msa`, `--simulated-dir`, `--seq-type`, `--threads`, `--output-dir`, `--overwrite`, `--resume`, `--dry-run`, and `--quiet`; each annotation must match Click help and state the format support, minimum 10 valid replicates, and resume constraint where applicable.

- [ ] **Step 3: Run documentation and full regression checks**

Run: `pytest tests/posttree/test_simulate_adequacy.py tests/cli/test_posttree_simulate_alisim.py tests/report/test_collector.py tests/report/test_templates.py tests/mcp -q`

Expected: PASS.

Run: `pytest -q`

Expected: PASS.

- [ ] **Step 4: Review Task 4 changes**

Do not commit without explicit user approval. Verify the documentation matches `phyloai posttree simulate adequacy --help` and retain the changes uncommitted.
