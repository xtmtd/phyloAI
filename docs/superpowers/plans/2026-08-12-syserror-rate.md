# Systematic-Error Site-Rate Analysis Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add `phyloai posttree syserror rate` to normalize IQ-TREE3/PhyloBayes site rates and optionally create slow- or fast-site MSA subsets.

**Architecture:** Keep parsing, normalization, rate ranking, optional alignment slicing, output lifecycle, and result payload creation in `phyloai/posttree/syserror_rate.py`. Keep Click as a thin validation/error-result wrapper. Reuse the existing `FormatConverter` to read supported MSA formats and Biopython only to emit standard FASTA; MCP discovers the Click command automatically.

**Tech Stack:** Python 3.10+, stdlib (`csv`, `json`, `math`, `shlex`, `shutil`, `time`), existing Biopython, Click, Rich, pytest.

## Global Constraints

- Implement only `docs/superpowers/specs/2026-08-12-syserror-rate-design.md`; add no generic site-filtering framework or dependency.
- Exactly one of `--iqtree-rate` and `--pb-rate` is required. IQ-TREE sites are 1-based; PhyloBayes sites are 0-based and must normalize to 1-based.
- Canonical ranking is `(rate ascending, site ascending)` and `rates.csv` has exactly `site,rate` columns.
- `--matrix` enables extraction and requires `--fraction`; omitted `--subset` resolves to `slow` only for extraction; without `--matrix`, explicitly supplied extraction options are invalid.
- `--fraction` is a comma-separated unique list in `(0, 1]`; retain `ceil(n_sites * fraction)` sites without expanding tied boundary rates.
- Reorder selected columns by original 1-based site coordinate before writing `positions.txt` and `matrix.fa`.
- Read all existing supported alignment formats through `FormatConverter`; write 60-column wrapped FASTA only.
- Default output is `runs/posttree/syserror/rate`; no external tools, no `doctor`, no checkpoints, and no resume support.
- Click path options must not use `exists=True`; `run_rate()` validates `is_file()` so non-dry CLI validation failures can write an error `result.json`; dry-runs write nothing.
- Validation failures never delete an existing output directory. With `--overwrite`, the CLI may replace/create only root `result.json` while preserving every other pre-existing file.
- Only `rates.csv` enters the report table index. `positions.txt` and `matrix.fa` remain persistent `data.output_files` entries only.
- Do not commit any change unless the user explicitly approves a commit.

---

## File Map

| Action | File | Responsibility |
|---|---|---|
| Create | `phyloai/posttree/syserror_rate.py` | Rate parsers, validation, ranking, MSA extraction, output files, and result payload. |
| Modify | `phyloai/cli/commands/posttree.py` | Register and wrap the `rate` Click command under `syserror`. |
| Create | `tests/posttree/test_syserror_rate.py` | Unit and runner tests for both inputs, slicing, errors, dry-runs, and payloads. |
| Create | `tests/cli/test_posttree_syserror_rate.py` | Click help, end-to-end behavior, error-result, report-ID, and generated MCP discovery tests. |
| Modify | `phyloai/report/collector.py` | Recognize `posttree.syserror.rate`. |
| Modify | `phyloai/report/templates.py` | Generate deterministic rate-selection methods prose. |
| Modify | `tests/report/test_collector.py` | Assert rate command step-ID parsing. |
| Modify | `tests/report/test_templates.py` | Assert rate methods text. |
| Create | `docs/commands/posttree-syserror-rate.md` | Detailed English command reference. |
| Create | `docs/commands/posttree-syserror-rate.zh.md` | Matching Chinese command reference. |
| Modify | `README.md`, `README.zh.md` | Add rate command example and reference-table row. |
| Modify | `docs/superpowers/specs/2026-06-07-phyloai-design.md` | Add `syserror rate` CLI example and atom sequence. |
| Modify | `skills/phyloai-workflow/SKILL.md` | Add local-only parameter review and scientific interpretation guidance. |
| Modify | `skills/phyloai-workflow/references/parameter-annotations.md` | Add Chinese annotations for all rate-command parameters. |

### Task 1: Rate Parsing And Canonical Ranking

**Files:**
- Create: `phyloai/posttree/syserror_rate.py`
- Create: `tests/posttree/test_syserror_rate.py`

**Interfaces:**
- Produces `RateRow(site: int, rate: float)`, `parse_iqtree_rate(path: Path) -> list[RateRow]`, `parse_pb_rate(path: Path) -> list[RateRow]`, `canonical_rates(rows: list[RateRow]) -> list[RateRow]`, and `rate_source(iqtree_rate: Path | None, pb_rate: Path | None) -> tuple[str, Path]`.

- [ ] **Step 1: Write the failing parser tests**

```python
from pathlib import Path

import pytest

from phyloai.posttree.syserror_rate import (
    RateRow, canonical_rates, parse_iqtree_rate, parse_pb_rate, rate_source,
)


def test_iqtree_rate_keeps_one_based_sites_and_sorts_ties(tmp_path: Path) -> None:
    source = tmp_path / "matrix.rate"
    source.write_text("# comment\nSite\tRate\tCat\n2\t0.5\t1\n1\t0.5\t1\n3\t1.0\t2\n")
    assert canonical_rates(parse_iqtree_rate(source)) == [
        RateRow(site=1, rate=0.5), RateRow(site=2, rate=0.5), RateRow(site=3, rate=1.0),
    ]


def test_pb_rate_converts_zero_based_sites_to_one_based(tmp_path: Path) -> None:
    source = tmp_path / "chain.meansiterates"
    source.write_text("0 1.2\n1 0.2\n2 0.6\n")
    assert canonical_rates(parse_pb_rate(source)) == [
        RateRow(site=2, rate=0.2), RateRow(site=3, rate=0.6), RateRow(site=1, rate=1.2),
    ]


@pytest.mark.parametrize("iqtree,pb", [(None, None), (Path("a.rate"), Path("b.meansiterates"))])
def test_exactly_one_source_is_required(iqtree: Path | None, pb: Path | None) -> None:
    with pytest.raises(ValueError, match="exactly one"):
        rate_source(iqtree, pb)
```

- [ ] **Step 2: Run the tests and verify they fail**

Run: `pytest tests/posttree/test_syserror_rate.py -q`

Expected: FAIL because `phyloai.posttree.syserror_rate` does not exist.

- [ ] **Step 3: Implement the minimal strict parsers**

```python
@dataclass(frozen=True)
class RateRow:
    site: int
    rate: float


def canonical_rates(rows: list[RateRow]) -> list[RateRow]:
    return sorted(rows, key=lambda row: (row.rate, row.site))


def rate_source(iqtree_rate: Path | None, pb_rate: Path | None) -> tuple[str, Path]:
    if (iqtree_rate is None) == (pb_rate is None):
        raise ValueError("exactly one of --iqtree-rate or --pb-rate is required")
    return ("iqtree", iqtree_rate) if iqtree_rate is not None else ("pb", pb_rate)
```

Implement a shared path/line validator that rejects absent paths, non-files, unreadable/empty files, malformed rows, non-integer sites, duplicate sites, non-finite rates (`math.isfinite`), negative rates, and any source index sequence other than `1..N` for IQ-TREE or `0..N-1` for PhyloBayes. For IQ-TREE, ignore `#` lines, require a tabular header containing exact `Site` and `Rate` columns, and parse fields by header index. For PhyloBayes, require exactly two whitespace-separated fields per non-empty row, then add one to every parsed site. Return source-specific errors that name the path and offending line.

- [ ] **Step 4: Add malformed-input regression tests and pass them**

```python
@pytest.mark.parametrize("text,match", [
    ("Site\tRate\n1\t-1\n", "negative"),
    ("Site\tRate\n1\tnan\n", "finite"),
    ("Site\tRate\n1\t0.1\n1\t0.2\n", "duplicate"),
    ("Site\tRate\n1\t0.1\n3\t0.2\n", "consecutive"),
])
def test_iqtree_rate_rejects_invalid_rows(tmp_path: Path, text: str, match: str) -> None:
    source = tmp_path / "bad.rate"
    source.write_text(text)
    with pytest.raises(ValueError, match=match):
        parse_iqtree_rate(source)
```

Run: `pytest tests/posttree/test_syserror_rate.py -q`

Expected: PASS.

- [ ] **Step 5: Review Task 1; do not commit**

Run: `git diff --check && git diff -- phyloai/posttree/syserror_rate.py tests/posttree/test_syserror_rate.py`

Expected: no whitespace errors.

### Task 2: Fraction Selection And Alignment Serialization

**Files:**
- Modify: `phyloai/posttree/syserror_rate.py`
- Modify: `tests/posttree/test_syserror_rate.py`

**Interfaces:**
- Consumes `RateRow` and `canonical_rates()` from Task 1.
- Produces `parse_fractions(value: str | None) -> list[float]`, `select_sites(ranked: list[RateRow], subset: str, fraction: float) -> list[int]`, `read_matrix(path: Path) -> MultipleSeqAlignment`, and `write_subset_fasta(alignment: MultipleSeqAlignment, sites: list[int], path: Path) -> None`.

- [ ] **Step 1: Write failing selection and FASTA-output tests**

```python
from Bio import AlignIO


def test_slow_and_fast_selection_round_up_and_restore_site_order() -> None:
    ranked = [RateRow(4, 0.1), RateRow(2, 0.2), RateRow(1, 0.3), RateRow(3, 0.4)]
    assert select_sites(ranked, "slow", 0.26) == [2, 4]
    assert select_sites(ranked, "fast", 0.5) == [1, 3]


def test_subset_fasta_is_aligned_and_wrapped(tmp_path: Path) -> None:
    matrix = tmp_path / "matrix.fa"
    matrix.write_text(">A\n" + "A" * 61 + "\n>B\n" + "C" * 61 + "\n")
    output = tmp_path / "subset.fa"
    write_subset_fasta(read_matrix(matrix), [1, 61], output)
    assert output.read_text() == ">A\nAA\n>B\nCC\n"
    assert AlignIO.read(output, "fasta").get_alignment_length() == 2
```

- [ ] **Step 2: Run the new tests and verify failure**

Run: `pytest tests/posttree/test_syserror_rate.py -q`

Expected: FAIL because the selection and MSA helper interfaces are undefined.

- [ ] **Step 3: Implement selection and MSA helpers**

```python
def select_sites(ranked: list[RateRow], subset: str, fraction: float) -> list[int]:
    count = math.ceil(len(ranked) * fraction)
    selected = ranked[:count] if subset == "slow" else ranked[-count:]
    return sorted(row.site for row in selected)


def write_subset_fasta(alignment: MultipleSeqAlignment, sites: list[int], path: Path) -> None:
    with open(path, "w", encoding="utf-8") as handle:
        for record in alignment:
            sequence = "".join(str(record.seq)[site - 1] for site in sites)
            handle.write(f">{record.description}\n")
            for start in range(0, len(sequence), 60):
                handle.write(f"{sequence[start:start + 60]}\n")
```

`parse_fractions` must reject absent/empty values, empty tokens, values outside `(0, 1]`, repeats by numeric value, and percentage-label collisions after formatting. Use a helper `subset_label(subset: str, fraction: float) -> str` based on `format(fraction * 100, ".12g")`, stripping a trailing `.0`, so `0.25`, `0.125`, and `1` produce `slow25`, `slow12.5`, and `slow100`.

`read_matrix` must use the module-level existing `FormatConverter`, reject parsing errors, zero records, empty sequences, duplicate record IDs, and unequal sequence lengths. `run_rate` will separately compare alignment length to rate count. Preserve IDs/descriptions and use only existing Biopython types; do not add a writer abstraction.

- [ ] **Step 4: Add validation and multi-format tests, then pass them**

```python
@pytest.mark.parametrize("value", ["", "0", "1.1", "0.25,0.25", "0.25,"])
def test_fraction_parser_rejects_invalid_lists(value: str) -> None:
    with pytest.raises(ValueError):
        parse_fractions(value)


def test_read_matrix_accepts_relaxed_phylip(tmp_path: Path) -> None:
    matrix = tmp_path / "matrix.phy"
    matrix.write_text("2 3\nA  ACG\nB  TTT\n")
    assert read_matrix(matrix).get_alignment_length() == 3
```

Run: `pytest tests/posttree/test_syserror_rate.py -q`

Expected: PASS.

- [ ] **Step 5: Review Task 2; do not commit**

Run: `git diff --check && git diff -- phyloai/posttree/syserror_rate.py tests/posttree/test_syserror_rate.py`

Expected: no whitespace errors.

### Task 3: Runner, Outputs, And Result JSON

**Files:**
- Modify: `phyloai/posttree/syserror_rate.py`
- Modify: `tests/posttree/test_syserror_rate.py`

**Interfaces:**
- Produces `run_rate(iqtree_rate: Path | None, pb_rate: Path | None, matrix: Path | None = None, subset: str | None = None, fraction: str | None = None, output_dir: Path = Path("runs/posttree/syserror/rate"), overwrite: bool = False, dry_run: bool = False, quiet: bool = False) -> dict[str, Any]`.

- [ ] **Step 1: Write failing runner tests**

```python
import csv
import json


def test_rate_only_writes_sorted_two_column_csv_and_result(tmp_path: Path) -> None:
    rate = tmp_path / "matrix.rate"
    rate.write_text("Site\tRate\n1\t0.4\n2\t0.1\n3\t0.2\n")
    result = run_rate(rate, None, output_dir=tmp_path / "out", quiet=True)
    assert list(csv.DictReader((tmp_path / "out" / "rates.csv").open())) == [
        {"site": "2", "rate": "0.1"}, {"site": "3", "rate": "0.2"}, {"site": "1", "rate": "0.4"},
    ]
    assert result["key_results"]["n_sites"] == 3
    assert json.loads((tmp_path / "out" / "result.json").read_text())["error"] is None


def test_multi_fraction_fast_outputs_positions_and_matrices(tmp_path: Path) -> None:
    rate = tmp_path / "rates"
    rate.write_text("0 0.1\n1 0.2\n2 0.3\n3 0.4\n")
    matrix = tmp_path / "matrix.fa"
    matrix.write_text(">A\nABCD\n>B\nWXYZ\n")
    result = run_rate(None, rate, matrix, subset="fast", fraction="0.25,0.5", output_dir=tmp_path / "out", quiet=True)
    assert (tmp_path / "out" / "fast25" / "positions.txt").read_text() == "4\n"
    assert (tmp_path / "out" / "fast50" / "matrix.fa").read_text() == ">A\nCD\n>B\nYZ\n"
    assert [item["selected_sites"] for item in result["key_results"]["subsets"]] == [1, 2]
```

- [ ] **Step 2: Run the runner tests and verify failure**

Run: `pytest tests/posttree/test_syserror_rate.py -q`

Expected: FAIL because `run_rate` is undefined.

- [ ] **Step 3: Implement output lifecycle and payload**

Implement this order:

```python
# 1. Validate source XOR and existence; parse/normalize/rank rates.
# 2. Validate extraction-option relationships before output-dir creation.
# 3. If matrix is given, parse/validate it, compare exact length, parse fractions.
# 4. Build command and full params; dry-run returns payload now without mkdir.
# 5. Refuse non-empty output unless overwrite; overwrite removes it; mkdir output.
# 6. Write rates.csv; then one <subset><percent>/positions.txt and matrix.fa per fraction.
# 7. Write result.json and return the same payload.
```

Use CSV `fieldnames=["site", "rate"]`; serialize rates with `str(row.rate)`. `positions.txt` ends with a newline. Add `data.output_files` entries named `rates`, `slow25_positions`, and `slow25_matrix`, each containing absolute/fully resolved path strings and descriptions. Only `rates.csv` has a report-table-compatible extension; positions and FASTA files must remain output-file records without changing global report extension rules. Populate:

```python
key_results = {
    "rate_source": source_name,
    "n_sites": len(ranked),
    "min_rate": ranked[0].rate,
    "max_rate": ranked[-1].rate,
    "subsets": [
        {"subset": subset, "requested_fraction": fraction, "selected_sites": count,
         "actual_fraction": count / len(ranked), "output_dir": str(directory)}
    ],
}
```

Use `tool_versions: {}` and `data={"cmd": [], "tool_stderr": "", "warnings": [], "output_files": ...}`. Terminal output names `rates.csv`, each selected matrix, and `result.json` unless quiet.

- [ ] **Step 4: Add lifecycle, error, and dry-run tests**

```python
def test_matrix_length_mismatch_fails_before_output_creation(tmp_path: Path) -> None:
    rate = tmp_path / "rates"
    rate.write_text("0 0.1\n1 0.2\n")
    matrix = tmp_path / "matrix.fa"
    matrix.write_text(">A\nABC\n>B\nABC\n")
    with pytest.raises(ValueError, match="length"):
        run_rate(None, rate, matrix, fraction="0.5", output_dir=tmp_path / "out", quiet=True)
    assert not (tmp_path / "out").exists()


def test_dry_run_validates_extraction_but_writes_nothing(tmp_path: Path) -> None:
    rate = tmp_path / "rates"
    rate.write_text("0 0.1\n1 0.2\n")
    matrix = tmp_path / "matrix.fa"
    matrix.write_text(">A\nAC\n>B\nTG\n")
    result = run_rate(None, rate, matrix, fraction="0.5", output_dir=tmp_path / "out", dry_run=True, quiet=True)
    assert result["data"]["output_files"] == {}
    assert not (tmp_path / "out").exists()


def test_extraction_options_without_matrix_are_rejected(tmp_path: Path) -> None:
    rate = tmp_path / "rates"
    rate.write_text("0 0.1\n")
    with pytest.raises(ValueError, match="--matrix"):
        run_rate(None, rate, subset="fast", fraction="0.5", output_dir=tmp_path / "out", quiet=True)
```

Also assert: a matrix without fraction fails, non-empty output is preserved without overwrite, overwrite succeeds, and `result.json` contains all shared top-level fields and output descriptions.

Add this overwrite-validation regression:

```python
def test_library_validation_never_deletes_existing_output_with_overwrite(tmp_path: Path) -> None:
    output = tmp_path / "out"
    output.mkdir()
    keep = output / "keep.txt"
    keep.write_text("existing data")
    with pytest.raises(ValueError, match="does not exist"):
        run_rate(tmp_path / "missing.rate", None, output_dir=output, overwrite=True, quiet=True)
    assert keep.read_text() == "existing data"
```

- [ ] **Step 5: Run focused tests and review; do not commit**

Run: `pytest tests/posttree/test_syserror_rate.py -q && git diff --check`

Expected: PASS with no whitespace errors.

### Task 4: CLI, MCP Discovery, And Reporting

**Files:**
- Modify: `phyloai/cli/commands/posttree.py`
- Create: `tests/cli/test_posttree_syserror_rate.py`
- Modify: `phyloai/report/collector.py`
- Modify: `phyloai/report/templates.py`
- Modify: `tests/report/test_collector.py`
- Modify: `tests/report/test_templates.py`

**Interfaces:**
- Consumes `run_rate()` from Task 3 and existing `_output_dir_writable()`, `_write_error_result_json()`, `parse_step_id()`, `STEP_ORDER`, and `generate_all_methods()`.
- Produces Click/MCP leaf `posttree_syserror_rate`, report step ID `posttree.syserror.rate`, and deterministic rate methods prose.

- [ ] **Step 1: Write failing integration tests**

```python
from click.testing import CliRunner

from phyloai.cli.main import cli
from phyloai.mcp.schema_gen import walk_click_tree
from phyloai.report.collector import parse_step_id
from phyloai.report.templates import generate_all_methods


def test_rate_help_and_generated_mcp_leaf() -> None:
    result = CliRunner().invoke(cli, ["posttree", "syserror", "rate", "--help"])
    assert result.exit_code == 0
    for option in ("--iqtree-rate", "--pb-rate", "--matrix", "--subset", "--fraction", "--dry-run"):
        assert option in result.output
    assert "posttree_syserror_rate" in {item["tool_name"] for item in walk_click_tree(cli)}


def test_rate_step_id_and_methods_text() -> None:
    assert parse_step_id("phyloai posttree syserror rate --pb-rate x") == "posttree.syserror.rate"
    text = generate_all_methods(
        "posttree.syserror.rate", {"pb_rate": "chain.meansiterates", "subset": "slow", "fraction": "0.25"},
        {"n_sites": 100, "subsets": [{"requested_fraction": 0.25, "selected_sites": 25}]}, {},
    )
    assert "PhyloBayes" in text and "slow" in text and "25" in text
```

- [ ] **Step 2: Run tests and verify failure**

Run: `pytest tests/cli/test_posttree_syserror_rate.py tests/report/test_collector.py tests/report/test_templates.py -q`

Expected: FAIL because `rate` is not registered and no report template exists.

- [ ] **Step 3: Add the Click command and error wrapper**

Register `@syserror.command("rate")` with `click.Path(path_type=Path)` for both rate inputs and matrix, deliberately without `exists=True`, `file_okay=False`, or `dir_okay=False`, `click.Choice(["slow", "fast"])` for `--subset` with Click default `None`, and the approved defaults. `run_rate()` alone verifies that supplied paths are readable regular files, allowing every missing or directory-input error to pass through the standard error-result lifecycle. In help, state that its extraction default is `slow`. `run_rate()` resolves `None` to `slow` only after confirming a matrix is present, so `--subset slow` without `--matrix` remains an input error. The help text must explain source XOR, source index normalization, slow/fast choice, comma-list fractions, `ceil` selection, output layout, no external tools, and these examples:

```bash
phyloai posttree syserror rate --iqtree-rate matrix.rate
phyloai posttree syserror rate --iqtree-rate matrix.rate --matrix raw.fa --fraction 0.25,0.5,0.75
phyloai posttree syserror rate --pb-rate chain1.meansiterates --matrix raw.phy --subset fast --fraction 0.1
```

Call `run_rate()` directly. Before it runs, record `pre_existing = output_dir.resolve().exists() and any(output_dir.resolve().iterdir())`. On `ValueError`, if not dry-run and either the directory was not pre-existing or `--overwrite` is true, call `_write_error_result_json()` then `_fail(..., exit_code=1)`; otherwise call `_fail(..., exit_code=1)` without touching the existing directory. This matches brlen's error-result lifecycle while ensuring validation failure never triggers overwrite cleanup: when `--overwrite` is true, writing root `result.json` preserves all other existing files. Print a JSON payload on successful dry-run, matching `brlen` behavior. Update `_SyserrorGroup.list_commands()` and its help text to include `rate`.

- [ ] **Step 4: Add report discovery and prose**

Add `"rate"` to `collector.py`'s `syserror` third-level set and place `"posttree.syserror.rate"` after brlen in `STEP_ORDER`. Add and register:

```python
def generate_methods_posttree_syserror_rate(
    params: dict[str, Any], key_results: dict[str, Any], tool_versions: dict[str, Any],
) -> str:
    source = "IQ-TREE empirical-Bayes site-rate estimates" if params.get("iqtree_rate") else "PhyloBayes posterior mean site rates"
    text = f"Site-rate heterogeneity was summarized from {source} across {_describe_n(key_results.get('n_sites', 0), 'alignment site', 'alignment sites')} using PhyloAI."
    subsets = key_results.get("subsets", [])
    if subsets:
        details = "; ".join(
            f"{item['subset']} sites: {_safe_fmt(item['requested_fraction'], '.1%')} ({item['selected_sites']} sites)"
            for item in subsets
        )
        text += f" Rate-ranked subsets retained {details}."
    return text
```

Keep `tool_versions` unused. The implementation must tolerate absent/malformed optional subset entries rather than causing report generation to fail.

- [ ] **Step 5: Add CLI lifecycle tests and pass the focused suite**

```python
def test_rate_cli_writes_error_result_but_dry_run_does_not(tmp_path) -> None:
    source = tmp_path / "rates"
    source.write_text("0 0.1\n")
    output = tmp_path / "out"
    failed = CliRunner().invoke(cli, ["posttree", "syserror", "rate", "--pb-rate", str(source), "--fraction", "0.5", "-o", str(output)])
    assert failed.exit_code == 1
    assert (output / "result.json").exists()
    dry_output = tmp_path / "dry"
    dry = CliRunner().invoke(cli, ["posttree", "syserror", "rate", "--pb-rate", str(source), "--fraction", "0.5", "-o", str(dry_output), "--dry-run"])
    assert dry.exit_code == 1
    assert not dry_output.exists()
```

Add a missing-input regression using `--pb-rate <missing path>` and a missing-matrix regression using `--matrix <missing path>`: both non-dry runs must exit 1 and write a root error `result.json`. Add a directory-input regression using an existing directory as `--pb-rate`, then as `--matrix`: both non-dry runs must exit 1 and write a root error `result.json`, proving Click did not preempt `run_rate()` validation. Add a pre-existing-output regression with `--overwrite`, invalid input, and `keep.txt`: it must preserve `keep.txt` while creating/replacing only `result.json`. Also add a successful CLI extraction assertion for `slow50/matrix.fa` and assert `rates.csv`/`positions.txt` are registered in result JSON.

Run: `pytest tests/cli/test_posttree_syserror_rate.py tests/report/test_collector.py tests/report/test_templates.py -q`

Expected: PASS.

- [ ] **Step 6: Review Task 4; do not commit**

Run: `git diff --check && git diff --stat`

Expected: no whitespace errors and no manually authored MCP schema/tool.

### Task 5: Documentation, Workflow Guidance, And Final Verification

**Files:**
- Create: `docs/commands/posttree-syserror-rate.md`
- Create: `docs/commands/posttree-syserror-rate.zh.md`
- Modify: `README.md`
- Modify: `README.zh.md`
- Modify: `docs/superpowers/specs/2026-06-07-phyloai-design.md`
- Modify: `skills/phyloai-workflow/SKILL.md`
- Modify: `skills/phyloai-workflow/references/parameter-annotations.md`
- Modify: `tests/report/test_schema.py`

**Interfaces:**
- Consumes final Click help and `result.json` from Tasks 3-4.
- Produces matching bilingual user documentation and safe workflow guidance.

- [ ] **Step 1: Write the English command reference**

Create `docs/commands/posttree-syserror-rate.md` using the established sections: Purpose, Usage, Inputs, Rate Inputs and Indexing, Outputs, Examples, Warnings / Errors, and Notes. Cover all CLI parameters and defaults, exact output layout, 1-based normalized `rates.csv`, strict consecutive-index validation, no tie expansion, ceiling rounding, original-order `positions.txt`, supported input MSA formats, 60-character FASTA wrapping, and no `doctor`/resume requirement. Explicitly state that rate filtering is a sensitivity analysis and does not prove a topology or automatically choose a fraction.

- [ ] **Step 2: Write matching Chinese documentation and README entries**

Create the Chinese document with the identical behavioral coverage. Add this concise example to the post-tree examples in both READMEs:

```bash
phyloai posttree syserror rate --iqtree-rate matrix.rate --matrix raw.fa \
  --subset slow --fraction 0.25,0.5,0.75 -o runs/posttree/syserror/rate
```

Add command-table entries linking to the language-matched reference documents. The description must call this a site-rate ranking/extraction sensitivity utility, not an automatic systematic-error correction.

- [ ] **Step 3: Update project and workflow guidance**

In the main design, add `phyloai posttree syserror rate --iqtree-rate ./matrix.rate --matrix ./matrix.fa --fraction 0.25,0.5,0.75` beside other syserror examples and change the high-level atomic diagnostic sequence from `brlen → cca → sites` to `brlen → rate → cca → sites`.

In `skills/phyloai-workflow/SKILL.md`, add a `posttree syserror rate` section: it is local-only, so no doctor check is needed, but schema review, parameter card, and explicit approval remain mandatory. Explain input selection, index normalization, slow versus fast scientific interpretation, recommended sensitivity fractions, and that downstream tree inference is never automatic.

Add Chinese annotations for `--iqtree-rate`, `--pb-rate`, `--matrix`, `--subset`, `--fraction`, `--output-dir`, `--overwrite`, `--dry-run`, and `--quiet`; state the source XOR and matrix/fraction dependency.

Add a report-index regression in `tests/report/test_schema.py` using a `posttree.syserror.rate` step with `data.output_files` containing `rates.csv`, `slow25/positions.txt`, and `slow25/matrix.fa`. Assert that `build_tables_index()` includes only the `rates` entry, while the assembled report step retains all three `output_files`. Do not widen `_TABLE_EXTENSIONS` or alter renderer behavior for `.txt` files.

- [ ] **Step 4: Mark the feature spec implemented after tests pass**

Change `docs/superpowers/specs/2026-08-12-syserror-rate-design.md` from `Status: Approved for specification review` to `Status: Implemented`, adding the implementation date only after every command, report, documentation, and test check in this plan passes.

- [ ] **Step 5: Run final verification; do not commit**

Run: `python -m phyloai.cli.main posttree syserror rate --help`

Expected: exit 0 with both rate sources, extraction dependencies, output explanation, and all three examples.

Run: `pytest tests/posttree/test_syserror_rate.py tests/cli/test_posttree_syserror_rate.py tests/report/test_collector.py tests/report/test_schema.py tests/report/test_templates.py -q`

Expected: PASS.

Run: `pytest -q`

Expected: PASS.

Run: `git diff --check && git status --short`

Expected: no whitespace errors; only planned files are changed, plus any pre-existing user changes left untouched.

## Final Acceptance Checklist

- [ ] IQ-TREE and PhyloBayes inputs normalize to exactly `1..N` and produce deterministic slow-to-fast `rates.csv`.
- [ ] Matrix-free runs write only `rates.csv` and standard `result.json`.
- [ ] Matrix runs produce correctly named one or more slow/fast directories, 1-based original-order positions, and 60-column FASTA subsets.
- [ ] Invalid source/index/rate/fraction/matrix/lifecycle inputs fail before unsafe output writes; CLI error JSON and dry-run behavior follow the project convention.
- [ ] Generated Click/MCP discovery exposes `posttree_syserror_rate`; no manual MCP implementation or dependency exists.
- [ ] Report methods text, bilingual documentation, README, design, and workflow guidance match final behavior.
- [ ] Targeted tests and `pytest -q` pass; no commit is made without explicit user approval.
