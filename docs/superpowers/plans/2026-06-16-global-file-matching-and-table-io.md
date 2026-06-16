# Global File Matching and Table I/O Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Implement the approved global logical-locus matching policy and unified CSV/TSV table interface across the current pretree codebase.

**Architecture:** Introduce a shared helper in `phyloai/core/file_matching.py` so filename parsing and MSA/tree pairing logic are defined once and reused by `pretree metrics` and future commands. Then migrate `pretree metrics` to that helper, replace the old `--per-gene-format` interface in `pretree stats` with the global `--table-format` interface, and update tests plus stale implementation-plan docs so design, code, and tests all agree.

**Tech Stack:** Python 3.10+, Click 8+, pathlib, csv, BioPython, pytest

---

## File Map

| Action | File | Responsibility |
|---|---|---|
| Create | `phyloai/core/file_matching.py` | Shared logical-locus parsing, pairing, ambiguity detection, and scan helpers |
| Modify | `phyloai/pretree/metrics.py` | Remove local suffix-based pairing and adopt the shared helper |
| Modify | `phyloai/pretree/stats.py` | Keep per-gene table writing aligned with the global `--table-format` interface |
| Modify | `phyloai/cli/commands/pretree.py` | Rename `pretree stats` flag from `--per-gene-format` to `--table-format` |
| Modify | `tests/pretree/test_metrics.py` | Replace old suffix-whitelist tests with logical-locus and ambiguity tests |
| Modify | `tests/pretree/test_stats.py` | Replace `--per-gene-format` CLI coverage with `--table-format` coverage |
| Modify | `docs/superpowers/plans/2026-06-14-pretree-metrics.md` | Remove obsolete iterative suffix-stripping guidance |
| Modify | `docs/superpowers/plans/2026-06-09-pretree-stats.md` | Keep the implementation plan text consistent with the final interface name |

---

### Task 1: Add shared logical-locus matching helpers

**Files:**
- Create: `phyloai/core/file_matching.py`
- Modify: `tests/pretree/test_metrics.py`

- [ ] **Step 1: Write failing tests for logical locus parsing and ambiguity handling**

```python
from pathlib import Path

import pytest

from phyloai.core.file_matching import (
    logical_msa_locus_name,
    logical_tree_locus_candidates,
    pair_msa_and_tree_maps,
)


def test_logical_msa_locus_name_uses_text_before_last_dot() -> None:
    assert logical_msa_locus_name(Path("gene.fa")) == "gene"
    assert logical_msa_locus_name(Path("gene.v1.ALI")) == "gene.v1"


def test_logical_tree_locus_candidates_try_one_or_two_suffixes() -> None:
    assert logical_tree_locus_candidates(Path("gene.fa.treefile")) == ("gene.fa", "gene")
    assert logical_tree_locus_candidates(Path("gene.tre")) == ("gene", None)


def test_pair_msa_and_tree_maps_raises_on_ambiguous_tree_name(tmp_path: Path) -> None:
    msa_map = {
        "gene": tmp_path / "gene.fa",
        "gene.fa": tmp_path / "gene.fa.fasta",
    }
    tree_map = {tmp_path / "gene.fa.treefile": None}

    with pytest.raises(ValueError, match="ambiguous"):
        pair_msa_and_tree_maps(msa_map, [tmp_path / "gene.fa.treefile"])
```

- [ ] **Step 2: Run the tests and verify they fail first**

Run: `pytest tests/pretree/test_metrics.py -k "logical_msa_locus_name or logical_tree_locus_candidates or ambiguous_tree_name" -v`
Expected: FAIL with import errors because `phyloai.core.file_matching` and its functions do not exist yet.

- [ ] **Step 3: Implement the shared helper module**

```python
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class PairingResult:
    paired: dict[str, tuple[Path | None, Path | None]]
    warnings: list[str]


def logical_msa_locus_name(path: Path) -> str:
    parts = path.name.rsplit(".", 1)
    return parts[0] if len(parts) == 2 else path.name


def logical_tree_locus_candidates(path: Path) -> tuple[str | None, str | None]:
    parts = path.name.split(".")
    if len(parts) == 1:
        return (path.name, None)
    one_suffix = ".".join(parts[:-1])
    two_suffix = ".".join(parts[:-2]) if len(parts) >= 3 else None
    return (one_suffix or None, two_suffix or None)


def pair_msa_and_tree_maps(
    msa_map: dict[str, Path],
    tree_paths: list[Path],
) -> PairingResult:
    paired = {locus: (path, None) for locus, path in msa_map.items()}
    warnings: list[str] = []

    for tree_path in sorted(tree_paths):
        one_suffix, two_suffix = logical_tree_locus_candidates(tree_path)
        one_hit = one_suffix in msa_map if one_suffix is not None else False
        two_hit = two_suffix in msa_map if two_suffix is not None else False

        if one_hit and two_hit and one_suffix != two_suffix:
            raise ValueError(
                f"Tree file '{tree_path.name}' is ambiguous: matches both '{one_suffix}' and '{two_suffix}'."
            )
        if one_hit:
            msa_path = msa_map[one_suffix]
            paired[one_suffix] = (msa_path, tree_path)
            continue
        if two_hit:
            msa_path = msa_map[two_suffix]
            paired[two_suffix] = (msa_path, tree_path)
            continue

        tree_locus = one_suffix or tree_path.name
        paired[tree_locus] = (None, tree_path)
        warnings.append(f"[WARN] Tree file '{tree_path.name}' has no matching MSA file.")

    return PairingResult(paired=paired, warnings=warnings)
```

- [ ] **Step 4: Run the helper tests and verify they pass**

Run: `pytest tests/pretree/test_metrics.py -k "logical_msa_locus_name or logical_tree_locus_candidates or ambiguous_tree_name" -v`
Expected: PASS

- [ ] **Step 5: Commit the helper module**

```bash
git add phyloai/core/file_matching.py tests/pretree/test_metrics.py
git commit -m "feat(core): add logical locus matching helpers"
```

---

### Task 2: Migrate `pretree metrics` to the shared matching policy

**Files:**
- Modify: `phyloai/pretree/metrics.py`
- Modify: `tests/pretree/test_metrics.py`

- [ ] **Step 1: Replace old suffix-stripping tests with shared-policy tests**

```python
def test_pair_files_matches_uppercase_and_nonstandard_msa_suffixes(tmp_path: Path) -> None:
    msa_dir = tmp_path / "msa"
    tree_dir = tmp_path / "trees"
    msa_dir.mkdir()
    tree_dir.mkdir()

    _write_fasta(msa_dir / "gene.v1.ALI", _AA_SEQS, ["A", "B", "C"])
    _write_newick(tree_dir / "gene.v1.ALI.treefile", "(A,B,C);")

    paired, warnings = _pair_files(msa_dir, tree_dir)

    assert warnings == []
    assert paired["gene.v1.ALI"][0] is not None
    assert paired["gene.v1.ALI"][1] is not None


def test_pair_files_raises_on_ambiguous_tree_candidate(tmp_path: Path) -> None:
    msa_dir = tmp_path / "msa"
    tree_dir = tmp_path / "trees"
    msa_dir.mkdir()
    tree_dir.mkdir()

    _write_fasta(msa_dir / "gene.fa", _AA_SEQS, ["A", "B", "C"])
    _write_fasta(msa_dir / "gene.fa.fasta", _AA_SEQS, ["A", "B", "C"])
    _write_newick(tree_dir / "gene.fa.treefile", "(A,B,C);")

    with pytest.raises(ValueError, match="ambiguous"):
        _pair_files(msa_dir, tree_dir)
```

- [ ] **Step 2: Run the pairing tests and verify they fail against the old implementation**

Run: `pytest tests/pretree/test_metrics.py -k "uppercase_and_nonstandard_msa_suffixes or ambiguous_tree_candidate" -v`
Expected: FAIL because `_pair_files()` still depends on hard-coded suffix logic and `_strip_tree_suffixes()`.

- [ ] **Step 3: Refactor `metrics.py` to use the shared helper**

```python
from phyloai.core.file_matching import (
    PairingResult,
    logical_msa_locus_name,
    pair_msa_and_tree_maps,
)


def _scan_msa_headers(msa_dir: Path) -> tuple[dict[str, set[str]], set[str]]:
    per_marker: dict[str, set[str]] = {}
    total_pool: set[str] = set()
    for path in sorted(msa_dir.iterdir()):
        if not path.is_file():
            continue
        try:
            records = _read_msa_records(path)
        except Exception:
            continue
        locus = logical_msa_locus_name(path)
        taxa = {record.id for record in records}
        per_marker[locus] = taxa
        total_pool.update(taxa)
    return per_marker, total_pool


def _pair_files(
    msa_dir: Path | None,
    tree_dir: Path | None,
) -> tuple[dict[str, tuple[Path | None, Path | None]], list[str]]:
    msa_map: dict[str, Path] = {}
    if msa_dir is not None:
        for path in sorted(msa_dir.iterdir()):
            if not path.is_file():
                continue
            msa_map[logical_msa_locus_name(path)] = path

    tree_paths = [] if tree_dir is None else [path for path in sorted(tree_dir.iterdir()) if path.is_file()]
    result = pair_msa_and_tree_maps(msa_map, tree_paths)

    for locus in sorted(set(msa_map) - {key for key, value in result.paired.items() if value[1] is not None}):
        result.warnings.append(f"[WARN] MSA file '{locus}' has no matching tree file.")

    return result.paired, result.warnings
```

- [ ] **Step 4: Delete the obsolete suffix-stripping helper and update imports**

```python
# Remove these constants and helper from metrics.py:
TREE_EXTENSIONS = ...
_TREE_SPECIFIC_SUFFIXES = ...
_GENERAL_SUFFIXES = ...

def _strip_tree_suffixes(name: str) -> str:
    ...
```

- [ ] **Step 5: Run the full metrics pairing test slice**

Run: `pytest tests/pretree/test_metrics.py -k "pair_files or scan_msa_headers or taxon_consistency" -v`
Expected: PASS

- [ ] **Step 6: Commit the metrics migration**

```bash
git add phyloai/pretree/metrics.py tests/pretree/test_metrics.py
git commit -m "refactor(metrics): use shared logical locus pairing"
```

---

### Task 3: Rename `pretree stats` table flag to `--table-format`

**Files:**
- Modify: `phyloai/cli/commands/pretree.py`
- Modify: `phyloai/pretree/stats.py`
- Modify: `tests/pretree/test_stats.py`

- [ ] **Step 1: Update CLI tests to use the new global flag name**

```python
def test_directory_table_format_can_write_tsv(tmp_path: Path) -> None:
    runner = CliRunner()
    output_dir = tmp_path / "stats_out"

    result = runner.invoke(
        cli,
        [
            "pretree",
            "stats",
            "--seq-dir",
            str(TEST_DATA / "test"),
            "--per-gene",
            "--table-format",
            "tsv",
            "--output-dir",
            str(output_dir),
        ],
    )

    assert result.exit_code == 0
    assert (output_dir / "per-gene.tsv").exists()
```

- [ ] **Step 2: Run the stats CLI tests and verify they fail first**

Run: `pytest tests/pretree/test_stats.py -k "table_format_can_write_tsv or directory_output_writes_summary_and_per_gene or directory_per_gene_writes_to_output_dir" -v`
Expected: FAIL because the CLI still exposes `--per-gene-format`.

- [ ] **Step 3: Rename the Click option and parameter in `pretree.py`**

```python
@click.option(
    "--table-format",
    type=click.Choice(["csv", "tsv"]),
    default="csv",
    show_default=True,
    help="Directory mode only. Table format for the per-gene file written with --per-gene.",
)
def stats_command(
    seq_dir: Path | None,
    seq: Path | None,
    per_gene: bool,
    table_format: str,
    output_dir: Path,
    input_format: str | None,
    seq_type: str | None,
    threads: int,
    quiet: bool,
    overwrite: bool,
) -> None:
    ...
    _run_stats_command(
        seq_dir,
        seq,
        per_gene,
        table_format,
        output_dir,
        input_format,
        seq_type,
        threads,
        quiet,
    )
```

- [ ] **Step 4: Keep `stats.py` signatures and path generation in sync with the renamed flag**

```python
def per_gene_output_path(summary_path: Path, output_format: str = "csv") -> Path:
    table_suffix = f".{output_format}"
    return summary_path.with_name(f"{summary_path.stem}.per-gene{table_suffix}")


def _run_stats_command(
    seq_dir: Path | None,
    seq: Path | None,
    per_gene: bool,
    table_format: str,
    output_dir: Path,
    input_format: str | None,
    seq_type: str | None,
    threads: int,
    quiet: bool,
) -> None:
    ...
    if per_gene:
        per_gene_path = output_dir / f"per-gene.{table_format}"
        write_per_gene_output(payload, per_gene_path)
```

- [ ] **Step 5: Run the targeted stats tests and then the CLI help test**

Run: `pytest tests/pretree/test_stats.py -k "table_format_can_write_tsv or directory_output_writes_summary_and_per_gene or directory_per_gene_writes_to_output_dir or cli_help_explains_options" -v`
Expected: PASS

- [ ] **Step 6: Commit the interface rename**

```bash
git add phyloai/cli/commands/pretree.py phyloai/pretree/stats.py tests/pretree/test_stats.py
git commit -m "refactor(stats): rename per-gene format flag"
```

---

### Task 4: Update stale implementation-plan docs to match the approved design

**Files:**
- Modify: `docs/superpowers/plans/2026-06-14-pretree-metrics.md`
- Modify: `docs/superpowers/plans/2026-06-09-pretree-stats.md`

- [ ] **Step 1: Update the old metrics plan so it no longer teaches suffix-whitelist stripping**

```md
- [ ] **Step 1.2: File pairing logic**

  ```python
  def _pair_files(
      msa_dir: Path | None,
      tree_dir: Path | None,
  ) -> tuple[dict[str, tuple[Path | None, Path | None]], list[str]]:
      """Match MSA and tree files by logical locus name using the shared helper."""
  ```

  - Derive MSA logical loci from the filename before the final `.`
  - For each tree file, try one-suffix and two-suffix reductions
  - If both tree candidates match different loci, raise an ambiguity error
  - Record unmatched files as warnings
```

- [ ] **Step 2: Update the old stats plan so it uses `--table-format` consistently**

```md
- `--table-format` (choice: csv/tsv, default csv)
- When `--per-gene` is provided, write per-gene results under `--output-dir` with CSV default and TSV override via `--table-format`
```

- [ ] **Step 3: Verify no stale plan text remains**

Run: `rg "per-gene-format|strip known suffixes|iterative suffix" docs/superpowers/plans`
Expected: no matches

- [ ] **Step 4: Commit the plan cleanup**

```bash
git add docs/superpowers/plans/2026-06-14-pretree-metrics.md docs/superpowers/plans/2026-06-09-pretree-stats.md
git commit -m "docs(plans): align pairing and table interfaces"
```

---

### Task 5: Final verification across code and docs

**Files:**
- Modify: none expected unless failures surface
- Test: `tests/pretree/test_metrics.py`, `tests/pretree/test_stats.py`

- [ ] **Step 1: Run the focused pretree test suites**

Run: `pytest tests/pretree/test_metrics.py tests/pretree/test_stats.py -v`
Expected: PASS

- [ ] **Step 2: Run a repo-wide grep for the retired CLI flag**

Run: `rg "per-gene-format"`
Expected: no matches

- [ ] **Step 3: Run a repo-wide grep for the removed suffix-stripping helper**

Run: `rg "_strip_tree_suffixes|TREE_EXTENSIONS|_TREE_SPECIFIC_SUFFIXES|_GENERAL_SUFFIXES" phyloai tests`
Expected: no matches for the removed pairing-only helpers; keep only any still-legitimate non-pairing constants if they remain for unrelated logic.

- [ ] **Step 4: Review `git diff` before handoff**

Run: `git diff -- phyloai/core/file_matching.py phyloai/pretree/metrics.py phyloai/pretree/stats.py phyloai/cli/commands/pretree.py tests/pretree/test_metrics.py tests/pretree/test_stats.py docs/superpowers/plans/2026-06-14-pretree-metrics.md docs/superpowers/plans/2026-06-09-pretree-stats.md`
Expected: only the planned logical-locus and table-format changes appear.

- [ ] **Step 5: Commit the final verification pass if any last fixes were needed**

```bash
git add phyloai/core/file_matching.py phyloai/pretree/metrics.py phyloai/pretree/stats.py phyloai/cli/commands/pretree.py tests/pretree/test_metrics.py tests/pretree/test_stats.py docs/superpowers/plans/2026-06-14-pretree-metrics.md docs/superpowers/plans/2026-06-09-pretree-stats.md
git commit -m "test: verify logical locus and table format migration"
```

---

## Self-Review

- Spec coverage: this plan covers the new global file matching policy, ambiguity-on-conflict behavior, and the unified `--input-format csv|tsv|auto` / `--table-format csv|tsv` interface where it currently impacts implemented pretree commands.
- Placeholder scan: no `TODO`, `TBD`, or “implement later” placeholders remain.
- Type consistency: the plan uses `logical_msa_locus_name`, `logical_tree_locus_candidates`, `pair_msa_and_tree_maps`, and `table_format` consistently across later tasks.
