# PhyloAI Pretree Concat Jackknife Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add `phyloai pretree concat jackknife` to generate reproducible gene-jackknife pseudoreplicate matrices from an existing concatenated matrix and partition file.

**Architecture:** Keep existing `phyloai pretree concat` behavior intact by converting it to a backward-compatible Click group with `invoke_without_command=True`. Put jackknife parsing/sampling/writing in `phyloai/pretree/concat.py` beside existing concat helpers, then wire CLI, report, MCP tests, and command docs around that library entry point.

**Tech Stack:** Python, Click, Biopython via existing `FormatConverter`, pytest, PhyloAI report/MCP schema generation.

## Global Constraints

- Follow `docs/superpowers/specs/2026-06-30-phyloai-pretree-concat-jackknife-design.md`.
- Preserve existing `phyloai pretree concat [OPTIONS]` behavior.
- No automatic `tree bi` execution.
- No `loci.txt` outputs.
- Default `--replicates 100`, `--target-length 50000`, `--prefix rep`, `--to fasta`, `--table-format csv`, `--seed 42`.
- Every non-`doctor` command writes exactly one `result.json` at output directory root.
- `data.output_files` must list `summary`, every `repXXX_matrix`, and every `repXXX_partitions` with absolute paths.
- FASTA outputs must wrap sequence lines at 60 characters.
- No new dependencies.

---

## File Structure

- Modify `phyloai/pretree/concat.py`: add partition parsing, replicate sampling, matrix slicing, summary writing, command builder, and `run_concat_jackknife()`.
- Modify `phyloai/cli/commands/pretree.py`: convert `concat` command to an invoke-without-command group and add nested `jackknife` command.
- Modify `phyloai/report/collector.py`: add `pretree.concat.jackknife` step parsing and ordering.
- Modify `phyloai/report/templates.py`: add methods text generator for jackknife.
- Modify `docs/superpowers/specs/2026-06-13-phyloai-pretree-concat-design.md`: note backward-compatible nested `jackknife` command.
- Modify `docs/commands/pretree-concat.md` and `docs/commands/pretree-concat.zh.md`: document jackknife usage.
- Modify tests under `tests/pretree/`, `tests/cli/`, `tests/report/`, and `tests/mcp/`.

---

### Task 1: Partition Parsing and Sampling Helpers

**Files:**
- Modify: `phyloai/pretree/concat.py`
- Test: `tests/pretree/test_concat_jackknife.py`

**Interfaces:**
- Produces: `_parse_partitions(path: Path) -> list[dict[str, Any]]`
- Produces: `_sample_partition_replicate(partitions: list[dict[str, Any]], target_length: int, rng: random.Random) -> list[dict[str, Any]]`

- [ ] **Step 1: Write failing helper tests**

Add `tests/pretree/test_concat_jackknife.py`:

```python
from __future__ import annotations

import random
from pathlib import Path

import pytest


def test_parse_partitions_raxml_style(tmp_path: Path) -> None:
    from phyloai.pretree.concat import _parse_partitions

    path = tmp_path / "matrix.partitions"
    path.write_text("LG, geneA = 1-10\nLG, geneB = 11-25\n")

    parts = _parse_partitions(path)

    assert parts == [
        {"model": "LG", "locus": "geneA", "start": 1, "end": 10, "length": 10},
        {"model": "LG", "locus": "geneB", "start": 11, "end": 25, "length": 15},
    ]


def test_parse_partitions_rejects_bad_line(tmp_path: Path) -> None:
    from phyloai.pretree.concat import _parse_partitions

    path = tmp_path / "bad.partitions"
    path.write_text("not a partition\n")

    with pytest.raises(ValueError, match="Unparseable partition line 1"):
        _parse_partitions(path)


def test_sample_partition_replicate_without_replacement_reaches_target() -> None:
    from phyloai.pretree.concat import _sample_partition_replicate

    parts = [
        {"model": "LG", "locus": "a", "start": 1, "end": 10, "length": 10},
        {"model": "LG", "locus": "b", "start": 11, "end": 20, "length": 10},
        {"model": "LG", "locus": "c", "start": 21, "end": 30, "length": 10},
    ]

    sampled = _sample_partition_replicate(parts, 20, random.Random(42))

    assert sum(p["length"] for p in sampled) >= 20
    assert len({p["locus"] for p in sampled}) == len(sampled)
```

- [ ] **Step 2: Run tests to verify failure**

Run: `pytest tests/pretree/test_concat_jackknife.py -v`

Expected: FAIL with import errors for `_parse_partitions` and `_sample_partition_replicate`.

- [ ] **Step 3: Implement helper functions**

Add near existing partition helpers in `phyloai/pretree/concat.py`:

```python
import random
import re
```

Add functions:

```python
_PARTITION_RE = re.compile(r"^\s*([^,]+)\s*,\s*(.+?)\s*=\s*(\d+)\s*-\s*(\d+)\s*$")


def _parse_partitions(path: Path) -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []
    with open(path) as fh:
        for line_no, raw in enumerate(fh, start=1):
            line = raw.strip()
            if not line:
                continue
            match = _PARTITION_RE.match(line)
            if match is None:
                raise ValueError(f"Unparseable partition line {line_no}: {raw.rstrip()}")
            model, locus, start_s, end_s = match.groups()
            start = int(start_s)
            end = int(end_s)
            if start < 1 or end < start:
                raise ValueError(f"Invalid partition range on line {line_no}: {start}-{end}")
            records.append({
                "model": model.strip(),
                "locus": locus.strip(),
                "start": start,
                "end": end,
                "length": end - start + 1,
            })
    if not records:
        raise ValueError(f"Partition file is empty: {path}")
    return records


def _sample_partition_replicate(
    partitions: list[dict[str, Any]],
    target_length: int,
    rng: random.Random,
) -> list[dict[str, Any]]:
    shuffled = list(partitions)
    rng.shuffle(shuffled)
    selected: list[dict[str, Any]] = []
    total = 0
    for part in shuffled:
        selected.append(part)
        total += int(part["length"])
        if total >= target_length:
            break
    return selected
```

- [ ] **Step 4: Run helper tests**

Run: `pytest tests/pretree/test_concat_jackknife.py -v`

Expected: PASS.

---

### Task 2: Jackknife Library Function and Outputs

**Files:**
- Modify: `phyloai/pretree/concat.py`
- Test: `tests/pretree/test_concat_jackknife.py`

**Interfaces:**
- Consumes: `_parse_partitions()`, `_sample_partition_replicate()` from Task 1.
- Consumes existing `_read_msa(path: Path) -> tuple[list[str], list[str], int]`, returning taxa, sequence strings, and alignment length as implemented in `phyloai/pretree/concat.py`.
- Produces: `run_concat_jackknife(matrix: Path, partitions: Path, output_dir: Path | None = None, replicates: int = 100, target_length: int = 50000, prefix: str = "rep", to: str = "fasta", table_format: str = "csv", seed: int = 42, overwrite: bool = False, dry_run: bool = False, quiet: bool = False) -> dict[str, Any]`
- Produces: `_build_concat_jackknife_command(...) -> str`

- [ ] **Step 1: Add failing end-to-end library test**

Append to `tests/pretree/test_concat_jackknife.py`:

```python
import json


def test_run_concat_jackknife_writes_replicate_dirs_and_result_json(tmp_path: Path) -> None:
    from phyloai.pretree.concat import run_concat_jackknife

    matrix = tmp_path / "matrix.fa"
    matrix.write_text(">A\nAAAACCCCGGGGTTTT\n>B\nAAAACCCCGGGGTTTT\n")
    parts = tmp_path / "matrix.partitions"
    parts.write_text("LG, gene1 = 1-4\nLG, gene2 = 5-8\nLG, gene3 = 9-12\nLG, gene4 = 13-16\n")
    out = tmp_path / "jackknife"

    payload = run_concat_jackknife(
        matrix=matrix,
        partitions=parts,
        output_dir=out,
        replicates=2,
        target_length=8,
        prefix="rep",
        to="fasta",
        table_format="csv",
        seed=42,
        overwrite=False,
        dry_run=False,
        quiet=True,
    )

    assert payload["status"] == "success"
    assert (out / "rep001" / "rep001.fa").exists()
    assert (out / "rep001" / "rep001.partitions").exists()
    assert (out / "rep002" / "rep002.fa").exists()
    assert (out / "jackknife_summary.csv").exists()
    saved = json.loads((out / "result.json").read_text())
    assert saved["params"]["seed"] == 42
    assert "summary" in saved["data"]["output_files"]
    assert "rep001_matrix" in saved["data"]["output_files"]
    assert "rep001_partitions" in saved["data"]["output_files"]


def test_run_concat_jackknife_rewrites_partition_coordinates(tmp_path: Path) -> None:
    from phyloai.pretree.concat import run_concat_jackknife

    matrix = tmp_path / "matrix.fa"
    matrix.write_text(">A\nAAAACCCCGGGG\n>B\nAAAACCCCGGGG\n")
    parts = tmp_path / "matrix.partitions"
    parts.write_text("LG, gene1 = 1-4\nLG, gene2 = 5-8\nLG, gene3 = 9-12\n")

    run_concat_jackknife(
        matrix=matrix,
        partitions=parts,
        output_dir=tmp_path / "out",
        replicates=1,
        target_length=8,
        seed=42,
        quiet=True,
    )

    lines = (tmp_path / "out" / "rep001" / "rep001.partitions").read_text().splitlines()
    assert lines[0].endswith("= 1-4")
    assert lines[1].endswith("= 5-8")
```

- [ ] **Step 2: Run tests to verify failure**

Run: `pytest tests/pretree/test_concat_jackknife.py -v`

Expected: FAIL with import error for `run_concat_jackknife`.

- [ ] **Step 3: Implement matrix slicing helpers and command builder**

Add stdlib imports to `phyloai/pretree/concat.py`. `json`, `shlex`, `shutil`, and `time` already exist in this module; add the missing imports beside them:

```python
import csv
import random
import re
```

Then add these helpers:

```python
def _table_suffix(table_format: str) -> str:
    return ".tsv" if table_format == "tsv" else ".csv"


def _table_delimiter(table_format: str) -> str:
    return "\t" if table_format == "tsv" else ","


def _matrix_extension(target_format: str) -> str:
    return {"fasta": ".fa", "phylip-relaxed": ".phy", "phylip-paml": ".phy", "nexus": ".nex"}.get(target_format, ".fa")


def _build_concat_jackknife_command(
    matrix: Path,
    partitions: Path,
    output_dir: Path,
    replicates: int,
    target_length: int,
    prefix: str,
    to: str,
    table_format: str,
    seed: int,
    overwrite: bool,
    dry_run: bool,
    quiet: bool,
) -> str:
    parts = [
        "phyloai", "pretree", "concat", "jackknife",
        "--matrix", str(matrix),
        "--partitions", str(partitions),
        "--replicates", str(replicates),
        "--target-length", str(target_length),
        "--prefix", prefix,
        "--to", to,
        "--table-format", table_format,
        "--seed", str(seed),
        "--output-dir", str(output_dir),
    ]
    if overwrite:
        parts.append("--overwrite")
    if dry_run:
        parts.append("--dry-run")
    if quiet:
        parts.append("--quiet")
    return shlex.join(parts)


def _validate_partition_bounds(partitions: list[dict[str, Any]], matrix_length: int) -> None:
    for part in partitions:
        if int(part["end"]) > matrix_length:
            raise ValueError(
                f"Partition {part['locus']!r} range {part['start']}-{part['end']} "
                f"exceeds matrix length {matrix_length}"
            )


def _slice_matrix_by_partitions(
    source_matrix: dict[str, str],
    selected: list[dict[str, Any]],
) -> dict[str, str]:
    sliced: dict[str, str] = {}
    for taxon, seq in source_matrix.items():
        pieces = [seq[int(part["start"]) - 1:int(part["end"])] for part in selected]
        sliced[taxon] = "".join(pieces)
    return sliced


def _rewrite_selected_partitions(selected: list[dict[str, Any]]) -> list[tuple[str, int, int, str]]:
    rewritten: list[tuple[str, int, int, str]] = []
    pos = 1
    for part in selected:
        length = int(part["length"])
        start = pos
        end = pos + length - 1
        rewritten.append((str(part["model"]), start, end, str(part["locus"])))
        pos = end + 1
    return rewritten


def _write_jackknife_partitions(path: Path, rewritten: list[tuple[str, int, int, str]]) -> None:
    lines = [f"{model}, {locus} = {start}-{end}\n" for model, start, end, locus in rewritten]
    path.write_text("".join(lines))
```

- [ ] **Step 4: Implement `run_concat_jackknife()`**

Add to `phyloai/pretree/concat.py`:

```python
def run_concat_jackknife(
    matrix: Path,
    partitions: Path,
    output_dir: Path | None = None,
    replicates: int = 100,
    target_length: int = 50000,
    prefix: str = "rep",
    to: str = "fasta",
    table_format: str = "csv",
    seed: int = 42,
    overwrite: bool = False,
    dry_run: bool = False,
    quiet: bool = False,
) -> dict[str, Any]:
    start_time = time.time()
    output_dir = (output_dir or (matrix.parent / "jackknife")).resolve()
    params = {
        "matrix": str(matrix),
        "partitions": str(partitions),
        "replicates": replicates,
        "target_length": target_length,
        "prefix": prefix,
        "to": to,
        "table_format": table_format,
        "seed": seed,
        "output_dir": str(output_dir),
        "overwrite": overwrite,
        "dry_run": dry_run,
        "quiet": quiet,
    }
    command = _build_concat_jackknife_command(
        matrix, partitions, output_dir, replicates, target_length, prefix, to,
        table_format, seed, overwrite, dry_run, quiet,
    )

    def _error_payload(message: str) -> dict[str, Any]:
        return {
            "status": "error",
            "command": command,
            "wall_time": round(time.time() - start_time, 3),
            "tool_versions": {},
            "params": params,
            "key_results": {},
            "error": message,
            "data": {"cmd": [], "tool_stderr": "", "output_files": {}, "replicates": [], "warnings": []},
        }

    try:
        if replicates < 1:
            raise ValueError("--replicates must be at least 1")
        if target_length < 1:
            raise ValueError("--target-length must be at least 1")
        if not prefix:
            raise ValueError("--prefix must not be empty")
        if table_format not in {"csv", "tsv"}:
            raise ValueError("--table-format must be csv or tsv")
        if not matrix.exists():
            raise ValueError(f"--matrix does not exist: {matrix}")
        if not partitions.exists():
            raise ValueError(f"--partitions does not exist: {partitions}")
        if output_dir.exists() and any(output_dir.iterdir()) and not overwrite:
            raise ValueError(f"Output directory '{output_dir}' is non-empty. Use --overwrite to replace.")
        if output_dir.exists() and any(output_dir.iterdir()) and overwrite and not dry_run:
            shutil.rmtree(output_dir)

        partition_records = _parse_partitions(partitions)
        available_length = sum(int(part["length"]) for part in partition_records)
        if available_length < target_length:
            raise ValueError(
                f"Total partition length {available_length} is less than --target-length {target_length}"
            )

        taxa, seqs, matrix_length = _read_msa(matrix)
        source_matrix = dict(zip(taxa, seqs))
        _validate_partition_bounds(partition_records, matrix_length)

        if not dry_run:
            output_dir.mkdir(parents=True, exist_ok=True)

        rng = random.Random(seed)
        ext = _matrix_extension(to)
        output_files: dict[str, dict[str, str]] = {}
        replicate_rows: list[dict[str, Any]] = []
        lengths: list[int] = []
        n_loci_values: list[int] = []

        for idx in range(1, replicates + 1):
            name = f"{prefix}{idx:03d}"
            selected = _sample_partition_replicate(partition_records, target_length, rng)
            rep_matrix = _slice_matrix_by_partitions(source_matrix, selected)
            rewritten = _rewrite_selected_partitions(selected)
            total_length = sum(int(part["length"]) for part in selected)
            rep_dir = output_dir / name
            matrix_path = rep_dir / f"{name}{ext}"
            part_path = rep_dir / f"{name}.partitions"
            if not dry_run:
                rep_dir.mkdir(parents=True, exist_ok=True)
                _write_matrix(rep_matrix, matrix_path, to, "AA")
                _write_jackknife_partitions(part_path, rewritten)
            output_files[f"{name}_matrix"] = {
                "path": str(matrix_path.resolve()),
                "description": f"Gene-jackknife pseudoreplicate matrix {name}",
            }
            output_files[f"{name}_partitions"] = {
                "path": str(part_path.resolve()),
                "description": f"Partition file for gene-jackknife pseudoreplicate {name}",
            }
            replicate_rows.append({
                "name": name,
                "matrix": str(matrix_path.resolve()),
                "partitions": str(part_path.resolve()),
                "n_loci": len(selected),
                "total_length": total_length,
                "loci": [str(part["locus"]) for part in selected],
            })
            lengths.append(total_length)
            n_loci_values.append(len(selected))

        summary_path = output_dir / f"jackknife_summary{_table_suffix(table_format)}"
        output_files["summary"] = {
            "path": str(summary_path.resolve()),
            "description": "Summary table for generated gene-jackknife pseudoreplicates",
        }
        if not dry_run:
            with open(summary_path, "w", newline="") as fh:
                writer = csv.DictWriter(
                    fh,
                    fieldnames=["replicate", "matrix", "partitions", "n_loci", "total_length", "target_length", "seed"],
                    delimiter=_table_delimiter(table_format),
                )
                writer.writeheader()
                for row in replicate_rows:
                    writer.writerow({
                        "replicate": row["name"],
                        "matrix": str(Path(row["matrix"]).relative_to(output_dir)),
                        "partitions": str(Path(row["partitions"]).relative_to(output_dir)),
                        "n_loci": row["n_loci"],
                        "total_length": row["total_length"],
                        "target_length": target_length,
                        "seed": seed,
                    })

        payload = {
            "status": "success",
            "command": command,
            "wall_time": round(time.time() - start_time, 3),
            "tool_versions": {},
            "params": params,
            "key_results": {
                "n_replicates": replicates,
                "target_length": target_length,
                "min_length": min(lengths),
                "max_length": max(lengths),
                "mean_length": round(sum(lengths) / len(lengths), 3),
                "min_loci": min(n_loci_values),
                "max_loci": max(n_loci_values),
            },
            "error": None,
            "data": {
                "cmd": [],
                "tool_stderr": "",
                "output_files": output_files,
                "replicates": replicate_rows,
                "warnings": [],
            },
        }
        if not dry_run:
            with open(output_dir / "result.json", "w") as fh:
                json.dump(payload, fh, indent=2)
        return payload
    except ValueError as exc:
        payload = _error_payload(str(exc))
        if not dry_run:
            output_dir.mkdir(parents=True, exist_ok=True)
            with open(output_dir / "result.json", "w") as fh:
                json.dump(payload, fh, indent=2)
        raise
```

- [ ] **Step 5: Run library tests**

Run: `pytest tests/pretree/test_concat_jackknife.py -v`

Expected: PASS.

- [ ] **Step 6: Add JSON compliance test**

Append:

```python
def test_run_concat_jackknife_result_json_schema(tmp_path: Path) -> None:
    from phyloai.pretree.concat import run_concat_jackknife
    from tests.helpers import validate_params_completeness, validate_result_json

    matrix = tmp_path / "matrix.fa"
    matrix.write_text(">A\nAAAACCCC\n>B\nAAAACCCC\n")
    parts = tmp_path / "matrix.partitions"
    parts.write_text("LG, gene1 = 1-4\nLG, gene2 = 5-8\n")
    payload = run_concat_jackknife(matrix, parts, tmp_path / "out", replicates=1, target_length=4, quiet=True)

    validate_result_json(payload)
    validate_params_completeness(payload, {
        "matrix", "partitions", "output_dir", "replicates", "target_length",
        "prefix", "to", "table_format", "seed", "overwrite", "dry_run", "quiet",
    })
```

- [ ] **Step 7: Run JSON compliance test**

Run: `pytest tests/pretree/test_concat_jackknife.py::test_run_concat_jackknife_result_json_schema -v`

Expected: PASS.

---

### Task 3: CLI Backward-Compatible Group and Jackknife Command

**Files:**
- Modify: `phyloai/cli/commands/pretree.py`
- Test: `tests/cli/test_pretree_concat.py`
- Test: `tests/cli/test_pretree_concat_jackknife.py`

**Interfaces:**
- Consumes: `run_concat()` and `_render_concat_panels()` existing behavior.
- Consumes: `run_concat_jackknife()` from Task 2.
- Produces: `phyloai pretree concat [OPTIONS]` unchanged and `phyloai pretree concat jackknife [OPTIONS]`.

- [ ] **Step 1: Add failing CLI tests**

Append to `tests/cli/test_pretree_concat.py`:

```python
def test_cli_pretree_concat_still_works_after_group_conversion(tmp_path: Path) -> None:
    msa_dir = tmp_path / "msas"
    msa_dir.mkdir()
    (msa_dir / "gene1.fa").write_text(">A\nACGT\n>B\nACGT\n")

    result = CliRunner().invoke(
        cli,
        ["pretree", "concat", "--msa-dir", str(msa_dir), "--output-dir", str(tmp_path / "out"), "--seq-type", "NT", "--quiet"],
    )

    assert result.exit_code == 0, result.output
    assert (tmp_path / "out" / "result.json").exists()
```

Create `tests/cli/test_pretree_concat_jackknife.py`:

```python
from __future__ import annotations

import json
from pathlib import Path

from click.testing import CliRunner

from phyloai.cli.main import cli


def test_cli_pretree_concat_jackknife_help() -> None:
    result = CliRunner().invoke(cli, ["pretree", "concat", "jackknife", "--help"])
    assert result.exit_code == 0
    for flag in ["--matrix", "--partitions", "--replicates", "--target-length", "--table-format", "--seed"]:
        assert flag in result.output


def test_cli_pretree_concat_jackknife_writes_outputs(tmp_path: Path) -> None:
    matrix = tmp_path / "matrix.fa"
    matrix.write_text(">A\nAAAACCCC\n>B\nAAAACCCC\n")
    parts = tmp_path / "matrix.partitions"
    parts.write_text("LG, gene1 = 1-4\nLG, gene2 = 5-8\n")
    out = tmp_path / "jackknife"

    result = CliRunner().invoke(
        cli,
        [
            "pretree", "concat", "jackknife",
            "--matrix", str(matrix),
            "--partitions", str(parts),
            "--replicates", "1",
            "--target-length", "4",
            "--table-format", "tsv",
            "--output-dir", str(out),
            "--quiet",
        ],
    )

    assert result.exit_code == 0, result.output
    assert (out / "rep001" / "rep001.fa").exists()
    assert (out / "jackknife_summary.tsv").exists()
    payload = json.loads((out / "result.json").read_text())
    assert payload["params"]["table_format"] == "tsv"
```

- [ ] **Step 2: Run CLI tests to verify failure**

Run: `pytest tests/cli/test_pretree_concat.py tests/cli/test_pretree_concat_jackknife.py -v`

Expected: FAIL because `concat jackknife` is not registered.

- [ ] **Step 3: Convert concat command to group without changing callback behavior**

In `phyloai/cli/commands/pretree.py`, replace `@pretree.command("concat", ...)` with:

```python
@pretree.group(
    "concat",
    invoke_without_command=True,
    help=(
        "Concatenate multiple MSA files into a supermatrix for phylogenetic inference. "
        "Supports occupancy filtering, recoding, codon variants, outgroup reordering, "
        "multi-format output, and gene-jackknife pseudoreplicates."
    ),
)
```

Add `@click.pass_context` immediately before `def concat_command(...):`, add `ctx: click.Context` as the first parameter, and add this as the first body line:

```python
    if ctx.invoked_subcommand is not None:
        return
```

Keep every existing option and body branch unchanged after that guard.

- [ ] **Step 4: Add nested `jackknife` CLI command**

Add below `concat_command()` in `phyloai/cli/commands/pretree.py`:

```python
@concat_command.command(
    "jackknife",
    help="Generate gene-jackknife pseudoreplicate matrices from a concatenated matrix and partitions.",
)
@click.option("--matrix", type=click.Path(exists=True, dir_okay=False, path_type=Path), required=True, help="Existing concatenated matrix.")
@click.option("--partitions", type=click.Path(exists=True, dir_okay=False, path_type=Path), required=True, help="RAxML-style partition file matching --matrix.")
@click.option("--replicates", type=click.IntRange(1, None), default=100, show_default=True, help="Number of pseudoreplicates to generate.")
@click.option("--target-length", type=click.IntRange(1, None), default=50000, show_default=True, help="Minimum sampled site length per pseudoreplicate.")
@click.option("--prefix", type=str, default="rep", show_default=True, help="Replicate file and directory prefix.")
@click.option("--to", "to", type=click.Choice(["fasta", "phylip-relaxed", "phylip-paml", "nexus"]), default="fasta", show_default=True, help="Output matrix format.")
@click.option("--table-format", type=click.Choice(["csv", "tsv"]), default="csv", show_default=True, help="Table format for jackknife_summary.")
@click.option("--seed", type=int, default=42, show_default=True, help="Random seed for reproducible sampling.")
@click.option("--output-dir", "-o", type=click.Path(file_okay=False, path_type=Path), default=None, help="Output directory. Default: <matrix_parent>/jackknife.")
@click.option("--overwrite", is_flag=True, default=False, help="Delete and recreate non-empty output directory.")
@click.option("--dry-run", is_flag=True, default=False, help="Validate inputs and report planned outputs without writing files.")
@click.option("--quiet", "-q", is_flag=True, default=False, help="Suppress terminal output.")
def concat_jackknife_command(
    matrix: Path,
    partitions: Path,
    replicates: int,
    target_length: int,
    prefix: str,
    to: str,
    table_format: str,
    seed: int,
    output_dir: Path | None,
    overwrite: bool,
    dry_run: bool,
    quiet: bool,
) -> None:
    from phyloai.pretree.concat import run_concat_jackknife

    try:
        payload = run_concat_jackknife(
            matrix=matrix,
            partitions=partitions,
            output_dir=output_dir,
            replicates=replicates,
            target_length=target_length,
            prefix=prefix,
            to=to,
            table_format=table_format,
            seed=seed,
            overwrite=overwrite,
            dry_run=dry_run,
            quiet=quiet,
        )
    except ValueError as exc:
        _fail(str(exc), 1)

    if dry_run:
        if not quiet:
            click.echo(f"[dry-run] Would generate {replicates} pseudoreplicates with target length {target_length}.", err=True)
        return
    if not quiet:
        click.echo(f"Pseudoreplicates saved to {payload['params']['output_dir']}", err=True)
        click.echo(f"Results saved to {Path(payload['params']['output_dir']) / 'result.json'}", err=True)
```

- [ ] **Step 5: Run CLI tests**

Run: `pytest tests/cli/test_pretree_concat.py tests/cli/test_pretree_concat_jackknife.py -v`

Expected: PASS.

---

### Task 4: Report Integration

**Files:**
- Modify: `phyloai/report/collector.py`
- Modify: `phyloai/report/templates.py`
- Test: `tests/report/test_collector.py`
- Test: `tests/report/test_templates.py`

**Interfaces:**
- Produces step id: `pretree.concat.jackknife`
- Produces methods generator: `generate_methods_pretree_concat_jackknife(params, key_results, tool_versions) -> str`

- [ ] **Step 1: Add failing report tests**

Append to `tests/report/test_collector.py`:

```python
def test_parse_step_id_pretree_concat_jackknife():
    assert parse_step_id("phyloai pretree concat jackknife --matrix matrix.fa --partitions matrix.partitions") == "pretree.concat.jackknife"


def test_step_order_places_concat_jackknife_after_concat():
    assert STEP_ORDER.index("pretree.concat") < STEP_ORDER.index("pretree.concat.jackknife")
    assert STEP_ORDER.index("pretree.concat.jackknife") < STEP_ORDER.index("tree.ml.fasttree")
```

Append to `tests/report/test_templates.py`:

```python
def test_methods_generator_pretree_concat_jackknife():
    from phyloai.report.templates import generate_all_methods

    steps = [{
        "step_id": "pretree.concat.jackknife",
        "status": "success",
        "params": {"replicates": 100, "target_length": 50000, "seed": 42, "to": "fasta"},
        "key_results": {"n_replicates": 100, "min_length": 50012, "max_length": 53280, "mean_length": 51140.5},
        "tool_versions": {},
    }]

    out = generate_all_methods(steps)
    assert "gene-jackknife" in out[0]["methods_text"]
    assert "100 pseudoreplicates" in out[0]["methods_text"]
    assert "50,000" in out[0]["methods_text"]
```

- [ ] **Step 2: Run report tests to verify failure**

Run: `pytest tests/report/test_collector.py::test_parse_step_id_pretree_concat_jackknife tests/report/test_collector.py::test_step_order_places_concat_jackknife_after_concat tests/report/test_templates.py::test_methods_generator_pretree_concat_jackknife -v`

Expected: FAIL because step id/order/template are absent.

- [ ] **Step 3: Update report collector**

In `phyloai/report/collector.py`, insert in `STEP_ORDER` immediately after `"pretree.concat"`:

```python
    "pretree.concat.jackknife",
```

Update `_THIRD_LEVEL`:

```python
        "concat": {"jackknife"},
```

- [ ] **Step 4: Add methods generator**

In `phyloai/report/templates.py`, add after `generate_methods_pretree_concat()`:

```python
def generate_methods_pretree_concat_jackknife(
    params: dict[str, Any],
    key_results: dict[str, Any],
    tool_versions: dict[str, Any],
) -> str:
    reps = params.get("replicates", key_results.get("n_replicates", 0))
    target = params.get("target_length", key_results.get("target_length", 0))
    seed = params.get("seed", 42)
    fmt = params.get("to", "fasta")
    mean_len = key_results.get("mean_length")
    min_len = key_results.get("min_length")
    max_len = key_results.get("max_length")
    text = (
        f"Gene-jackknife pseudoreplicate matrices were generated from the concatenated "
        f"supermatrix using phyloai pretree concat jackknife. "
        f"A total of {_describe_n(reps, 'pseudoreplicate')} were sampled without replacement "
        f"from partition-defined loci until each reached at least {_safe_fmt(target, ',')} sites "
        f"(`--target-length {target}`), using random seed {seed}. "
        f"Pseudoreplicates were written in {fmt} format."
    )
    if mean_len is not None and min_len is not None and max_len is not None:
        text += (
            f" Final pseudoreplicate lengths ranged from {_safe_fmt(min_len, ',')} "
            f"to {_safe_fmt(max_len, ',')} sites "
            f"(mean {_safe_fmt(mean_len, ',.1f')})."
        )
    return text
```

Register in `METHOD_GENERATORS`:

```python
    "pretree.concat.jackknife": generate_methods_pretree_concat_jackknife,
```

- [ ] **Step 5: Run report tests**

Run: `pytest tests/report/test_collector.py tests/report/test_templates.py -v`

Expected: PASS.

---

### Task 5: MCP Schema Test

**Files:**
- Modify: `tests/mcp/test_schema_gen.py`

**Interfaces:**
- Consumes Click tree from Task 3.
- Verifies tool name `pretree_concat_jackknife` and schema defaults.

- [ ] **Step 1: Add MCP failing test**

Append to `tests/mcp/test_schema_gen.py`:

```python
def test_walk_click_tree_finds_pretree_concat_jackknife() -> None:
    descriptor = next(d for d in walk_click_tree(cli) if d["tool_name"] == "pretree_concat_jackknife")
    tool_def = build_mcp_tool(descriptor)

    assert tool_def["name"] == "pretree_concat_jackknife"
    props = tool_def["inputSchema"]["properties"]
    assert props["replicates"]["default"] == 100
    assert props["target_length"]["default"] == 50000
    assert props["seed"]["default"] == 42
    assert props["table_format"]["enum"] == ["csv", "tsv"]
```

- [ ] **Step 2: Run MCP test**

Run: `pytest tests/mcp/test_schema_gen.py::test_walk_click_tree_finds_pretree_concat_jackknife -v`

Expected: PASS after Task 3.

---

### Task 6: Command Docs and Parent Spec Sync

**Files:**
- Modify: `docs/commands/pretree-concat.md`
- Modify: `docs/commands/pretree-concat.zh.md`
- Modify: `docs/superpowers/specs/2026-06-13-phyloai-pretree-concat-design.md`

**Interfaces:**
- Consumes implemented CLI behavior from Task 3.
- Produces user-facing docs for jackknife.

- [ ] **Step 1: Update English command docs**

In `docs/commands/pretree-concat.md`, add a section after the main concat usage:

```markdown
## Gene-Jackknife Pseudoreplicates

`phyloai pretree concat jackknife` creates pseudoreplicate matrices from an existing concatenated matrix and partition file. It does not infer trees.

```bash
phyloai pretree concat jackknife \
  --matrix runs/pretree/concat/matrix.fa \
  --partitions runs/pretree/concat/matrix.partitions \
  --replicates 100 \
  --target-length 50000 \
  --to fasta \
  --table-format csv \
  --seed 42 \
  -o runs/pretree/concat/jackknife
```

Outputs are written as one directory per replicate:

```text
jackknife/
├── rep001/
│   ├── rep001.fa
│   └── rep001.partitions
├── rep002/
├── jackknife_summary.csv
└── result.json
```

The selected loci for each replicate are recorded in the corresponding `repXXX.partitions` file and in `result.json`.
```

- [ ] **Step 2: Update Chinese command docs**

In `docs/commands/pretree-concat.zh.md`, add equivalent content:

```markdown
## Gene-jackknife 伪重复矩阵

`phyloai pretree concat jackknife` 从已有的 concatenated matrix 和 partition 文件生成伪重复矩阵，不会自动推断树。

```bash
phyloai pretree concat jackknife \
  --matrix runs/pretree/concat/matrix.fa \
  --partitions runs/pretree/concat/matrix.partitions \
  --replicates 100 \
  --target-length 50000 \
  --to fasta \
  --table-format csv \
  --seed 42 \
  -o runs/pretree/concat/jackknife
```

每个 replicate 单独一个目录：

```text
jackknife/
├── rep001/
│   ├── rep001.fa
│   └── rep001.partitions
├── rep002/
├── jackknife_summary.csv
└── result.json
```

每个 replicate 抽到的 locus 可从对应的 `repXXX.partitions` 和 `result.json` 读取。
```

- [ ] **Step 3: Update parent concat spec**

In `docs/superpowers/specs/2026-06-13-phyloai-pretree-concat-design.md`, add under CLI section:

```markdown
`phyloai pretree concat` is implemented as a backward-compatible Click group: invoking `phyloai pretree concat [OPTIONS]` runs the original full concatenation, while `phyloai pretree concat jackknife [OPTIONS]` generates pseudoreplicates from an existing matrix and partition file. The jackknife design is specified in `2026-06-30-phyloai-pretree-concat-jackknife-design.md`.
```

- [ ] **Step 4: Run docs smoke checks**

Run: `pytest tests/report/test_templates.py::TestReportTemplates::test_all_step_ids_have_generators -v`

Expected: PASS.

---

### Task 7: Full Verification

**Files:**
- No code changes unless verification reveals failures.

**Interfaces:**
- Consumes all prior tasks.
- Produces verified feature branch state.

- [ ] **Step 1: Run focused test suite**

Run:

```bash
pytest \
  tests/pretree/test_concat_jackknife.py \
  tests/cli/test_pretree_concat.py \
  tests/cli/test_pretree_concat_jackknife.py \
  tests/report/test_collector.py \
  tests/report/test_templates.py \
  tests/mcp/test_schema_gen.py \
  -v
```

Expected: PASS.

- [ ] **Step 2: Run full existing concat tests**

Run: `pytest tests/pretree/test_concat.py tests/cli/test_pretree_concat.py -v`

Expected: PASS, confirming the group conversion did not break existing concat behavior.

- [ ] **Step 3: Inspect result output manually with CLI smoke test**

Run:

```bash
tmpdir=$(mktemp -d) && \
mkdir -p "$tmpdir/in" && \
printf '>A\nAAAACCCCGGGGTTTT\n>B\nAAAACCCCGGGGTTTT\n' > "$tmpdir/matrix.fa" && \
printf 'LG, gene1 = 1-4\nLG, gene2 = 5-8\nLG, gene3 = 9-12\nLG, gene4 = 13-16\n' > "$tmpdir/matrix.partitions" && \
phyloai pretree concat jackknife --matrix "$tmpdir/matrix.fa" --partitions "$tmpdir/matrix.partitions" --replicates 2 --target-length 8 --output-dir "$tmpdir/out" --quiet && \
test -f "$tmpdir/out/rep001/rep001.fa" && \
test -f "$tmpdir/out/rep001/rep001.partitions" && \
test -f "$tmpdir/out/jackknife_summary.csv" && \
test -f "$tmpdir/out/result.json"
```

Expected: command exits 0.

- [ ] **Step 4: Run broader regression suite if time allows**

Run: `pytest tests/cli tests/pretree tests/report tests/mcp -v`

Expected: PASS.

---

## Self-Review Notes

- Spec coverage: tasks cover core sampling, output layout, `--table-format`, seed 42, output files, FASTA writing through existing writer, report, MCP, docs, and parent spec sync.
- Deliberate simplification: jackknife writes matrices using `_write_matrix(..., seq_type="AA")`; this is safe for FASTA and standard alignment serialization, but if downstream requires molecule metadata in Nexus for NT/recoded matrices, add `--seq-type` later. Keep it out now because the spec does not require molecule-type inference for jackknife.
- No new dependencies are introduced.
