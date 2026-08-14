# CCA Systematic-Error Diagnostic Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add `phyloai posttree syserror cca`, a local CCA diagnostic that combines one `.sitefreq` input and two site-likelihood tables into `cca.csv`, a ggplot-faithful `cca.pdf`, and standard `result.json`.

**Architecture:** Keep all parsing, site-set validation, Keff/delta calculations, fixed-bin aggregation, Matplotlib drawing, output lifecycle, and payload construction in one focused pure-Python module, `phyloai/posttree/syserror_cca.py`. Keep Click limited to declared options, command reconstruction, and the existing error-result lifecycle. The existing dynamic Click-tree MCP generator exposes the command after removing its stub; report integration remains a small template plus the existing generic output-file indexer.

**Tech Stack:** Python 3.10+, stdlib (`csv`, `math`, `shlex`, `shutil`, `time`), installed Matplotlib, Click, pytest, existing PhyloAI result-schema/report/MCP helpers.

**Spec:** `docs/superpowers/specs/2026-08-13-phyloai-posttree-syserror-cca-design.md` and `docs/superpowers/specs/2026-06-07-phyloai-design.md`

## Global Constraints

- Implement only the approved CCA design; do not add an R runtime, entropy Keff, automatic upstream IQ-TREE/PhyloBayes execution, automatic site filtering, or a topology/model verdict.
- `--site-freq` accepts only one-based, consecutive `.sitefreq` rows with exactly 20 finite non-negative frequencies that sum to one within `1e-6`; reject raw `.siteprofiles` by format rather than filename.
- Both likelihood CSVs require exactly named columns `site`, `lnL_Tree1`, and `lnL_Tree2`; extra columns, including `ΔSLS`, are ignored. All three inputs must use the same complete, one-based consecutive `1..N` site set.
- Compute `keff = 1 / sum(p_i ** 2)` and compute CCA's independent, reverse-sign contrast as `lnl_tree2 - lnl_tree1`; never reuse `ΔSLS`.
- `cca.csv` has exactly `model,site,keff,lnl_tree1,lnl_tree2,delta_lnl_tree2_tree1`, ordered by site then model 1 then model 2. Default labels remain literal `model1` / `model2`; do not infer names from input directories.
- Aggregate `floor(keff)` into fixed bins `1..20`, filling every model/bin absence with zero; a valid mathematical 20 must remain bin 20.
- Reproduce the documented R visual semantics: bars centred at `bin + 0.5`, width/dodge 1, colors `#F8766D` / `#00BFC4`, x breaks/limits 1–20 with no expansion, `ymin=min(0, min_sum*1.1)`, `ymax=max(0,max_sum*1.1)`, orange/blue half-plane shading, and the specified legend placement/style.
- Default output is `runs/posttree/syserror/cca` and contains only `cca.csv`, `cca.pdf`, and `result.json`; no PNG, checkpoint, resume, doctor check, or new dependency.
- `result.json` uses `{}` tool versions, resolved parameters, `data.cmd=[]`, `data.tool_stderr=""`, warnings, absolute output-file paths/descriptions, and command reconstruction with resolved paths.
- Click paths must not use `exists=True`; the library validates paths so non-dry CLI failures can write an error `result.json`. Dry runs validate/calculates but write nothing.
- Validation failures must precede overwrite deletion. A pre-existing non-empty directory stays untouched without overwrite; with `--overwrite`, an invalid CLI input may replace only root `result.json`, preserving other existing files.
- Update English and Chinese documentation, both READMEs, parent design Sections 4.1/8 (if needed)/11 Phase 9, MCP documentation in both languages, workflow Skill and parameter annotations, and report methods text.
- Do not commit unless the user explicitly approves a commit.

---

## File Map

| Action | File | Responsibility |
|---|---|---|
| Create | `phyloai/posttree/syserror_cca.py` | Strict parsers, CCA calculation/binning, PDF drawing, output lifecycle, and result payload. |
| Modify | `phyloai/cli/commands/posttree.py` | Register grouped `cca` Click command and standard error-result wrapper. |
| Create | `tests/posttree/test_syserror_cca.py` | Unit/runner tests for parsing, calculations, bins, PDF, lifecycle, and JSON. |
| Create | `tests/cli/test_posttree_syserror_cca.py` | Help, CLI lifecycle, dynamic MCP discovery, and stub-removal tests. |
| Modify | `phyloai/mcp/tools/stubs.py` | Remove only `posttree_syserror_cca` from unavailable stubs. |
| Modify | `phyloai/report/templates.py` | Replace placeholder CCA prose with input-aware deterministic methods text. |
| Modify | `tests/report/test_templates.py` | Verify CCA methods prose and malformed-data resilience. |
| Create | `docs/commands/posttree-syserror-cca.md` | English CCA reference. |
| Create | `docs/commands/posttree-syserror-cca.zh.md` | Matching Chinese CCA reference. |
| Modify | `README.md`, `README.zh.md` | Add CCA example and command-table row. |
| Modify | `docs/superpowers/specs/2026-06-07-phyloai-design.md` | Replace obsolete CCA interface and clarify Phase 9 input order. |
| Modify | `docs/commands/ai-integration.md`, `docs/commands/ai-integration.zh.md` | Move generated CCA MCP tool from the stub list to the fire-and-forget tool table. |
| Modify | `skills/phyloai-workflow/SKILL.md` | Add CCA parameter-review, input-preparation, and cautious interpretation guidance. |
| Modify | `skills/phyloai-workflow/references/parameter-annotations.md` | Add Chinese annotations for every CCA parameter. |
| Modify | `docs/superpowers/specs/2026-08-13-phyloai-posttree-syserror-cca-design.md` | Mark the approved CCA design implemented only after all verification passes. |

### Task 1: Strict CCA Input Parsers And Site Validation

**Files:**
- Create: `phyloai/posttree/syserror_cca.py`
- Create: `tests/posttree/test_syserror_cca.py`

**Interfaces:**
- Produces `SiteFrequency(site: int, frequencies: tuple[float, ...])`, `SiteLikelihood(site: int, lnl_tree1: float, lnl_tree2: float)`, strict site-frequency/LNL parsing functions, and a cross-input site-set validator.

- [ ] **Step 1: Write failing site-frequency parser tests**

```python
from pathlib import Path

import pytest

from phyloai.posttree.syserror_cca import parse_site_freq


def test_parse_site_freq_accepts_one_based_20_state_rows(tmp_path: Path) -> None:
    source = tmp_path / "chain.sitefreq"
    source.write_text("1 " + " ".join(["0.05"] * 20) + "\n")

    row = parse_site_freq(source)[0]

    assert row.site == 1
    assert row.frequencies == (0.05,) * 20


@pytest.mark.parametrize("row,match", [
    ("0 " + " ".join(["0.05"] * 20), "consecutive from 1"),
    ("1 " + " ".join(["0.05"] * 19), "exactly 20"),
    ("1 " + " ".join(["0.10"] * 20), "sum to 1"),
    ("1 " + " ".join(["-0.05"] + ["0.05"] * 19), "non-negative"),
])
def test_parse_site_freq_rejects_invalid_rows(tmp_path: Path, row: str, match: str) -> None:
    source = tmp_path / "bad.sitefreq"
    source.write_text(row + "\n")
    with pytest.raises(ValueError, match=match):
        parse_site_freq(source)
```

- [ ] **Step 2: Run parser tests and verify failure**

Run: `pytest tests/posttree/test_syserror_cca.py -q`

Expected: FAIL because `phyloai.posttree.syserror_cca` does not exist.

- [ ] **Step 3: Implement strict `.sitefreq` parsing**

```python
@dataclass(frozen=True)
class SiteFrequency:
    site: int
    frequencies: tuple[float, ...]


def parse_site_freq(path: Path) -> list[SiteFrequency]:
    rows = []
    for line_number, line in enumerate(_read_nonempty_file(path, "site-frequency"), 1):
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        fields = line.split()
        if len(fields) != 21:
            raise ValueError(f"site-frequency file {path}, line {line_number}: exactly 20 frequencies are required")
        rows.append(SiteFrequency(_parse_site(fields[0], path, line_number), _parse_frequencies(fields[1:], path, line_number)))
    return _validate_one_based_sites(rows, path, "site-frequency")
```

Use one private helper for regular readable-file checks and another for duplicate/consecutive one-based identifiers. `_validate_one_based_sites()` must return `sorted(rows, key=lambda row: row.site)`, so both parsers return ascending site order regardless of their input order. Include path and source-line number in malformed-row errors. Do not inspect extensions: the row contract, not `.sitefreq` suffix, rejects raw `.siteprofiles`.

- [ ] **Step 4: Write failing LNL and cross-input tests**

```python
from phyloai.posttree.syserror_cca import parse_site_lnl, validate_matching_sites


def test_lnl_parser_requires_named_columns_but_ignores_extra_columns(tmp_path: Path) -> None:
    source = tmp_path / "site_lnl.csv"
    source.write_text(
        "site,lnL_Tree1,lnL_Tree2,ΔSLS,support\n"
        "2,-2.0,-1.0,-1.0,Tree1\n"
        "1,-4.0,-3.5,-0.5,Tree1\n"
    )

    assert parse_site_lnl(source) == [
        SiteLikelihood(site=1, lnl_tree1=-4.0, lnl_tree2=-3.5),
        SiteLikelihood(site=2, lnl_tree1=-2.0, lnl_tree2=-1.0),
    ]


def test_matching_inputs_reject_different_site_sets(tmp_path: Path) -> None:
    freq = tmp_path / "x.sitefreq"
    freq.write_text("1 " + " ".join(["0.05"] * 20) + "\n")
    lnl = tmp_path / "x.csv"
    lnl.write_text("site,lnL_Tree1,lnL_Tree2\n2,-1,-2\n")

    with pytest.raises(ValueError, match="site sets must match"):
        validate_matching_sites(parse_site_freq(freq), parse_site_lnl(lnl), parse_site_lnl(lnl))
```

- [ ] **Step 5: Implement LNL parsing and matching validation**

```python
@dataclass(frozen=True)
class SiteLikelihood:
    site: int
    lnl_tree1: float
    lnl_tree2: float


def parse_site_lnl(path: Path) -> list[SiteLikelihood]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        required = {"site", "lnL_Tree1", "lnL_Tree2"}
        if reader.fieldnames is None or not required <= set(reader.fieldnames):
            raise ValueError(f"site-likelihood file {path}: header requires site, lnL_Tree1, lnL_Tree2")
        rows = [SiteLikelihood(_parse_site(row["site"], path, line), _finite_float(row["lnL_Tree1"], path, line), _finite_float(row["lnL_Tree2"], path, line)) for line, row in enumerate(reader, 2)]
    return _validate_one_based_sites(rows, path, "site-likelihood")
```

Use `csv.DictReader`, ignore every column except the required names, and validate CSV structural errors such as a missing header. The common `_validate_one_based_sites()` return value sorts LNL rows by site before `build_cca_rows()` zips the three sources. `validate_matching_sites()` compares the full ordered normalized IDs and raises a single clear `site sets must match across --site-freq, --site-lnl1, and --site-lnl2` error. Add tests for missing required header, duplicate/non-consecutive LNL IDs, `nan` likelihoods, malformed CSV rows, and input LNL rows in non-site order.

- [ ] **Step 6: Run focused tests and review; do not commit**

Run: `pytest tests/posttree/test_syserror_cca.py -q && git diff --check`

Expected: PASS with no whitespace errors.

### Task 2: CCA Rows, Fixed Bins, And R-Faithful Plot Helper

**Files:**
- Modify: `phyloai/posttree/syserror_cca.py`
- Modify: `tests/posttree/test_syserror_cca.py`

**Interfaces:**
- Consumes `SiteFrequency` and `SiteLikelihood` from Task 1.
- Produces `CcaRow(model: str, site: int, keff: float, lnl_tree1: float, lnl_tree2: float, delta_lnl_tree2_tree1: float)`, `build_cca_rows(site_freq: list[SiteFrequency], model1: list[SiteLikelihood], model2: list[SiteLikelihood], model1_name: str, model2_name: str) -> list[CcaRow]`, `summarize_bins(rows: list[CcaRow], model_names: tuple[str, str]) -> dict[str, list[float]]`, and `plot_cca(bin_sums: dict[str, list[float]], model_names: tuple[str, str], title: str, xlabel: str, ylabel: str, fig_width: float, fig_height: float, dpi: int, font_size: float, path: Path) -> None`.

- [ ] **Step 1: Write failing calculation and bin-completion tests**

```python
from phyloai.posttree.syserror_cca import (
    CcaRow, SiteFrequency, SiteLikelihood, build_cca_rows, summarize_bins,
)


def test_cca_rows_use_inverse_homozygosity_and_tree2_minus_tree1() -> None:
    frequencies = (0.2,) * 5 + (0.0,) * 15
    rows = build_cca_rows(
        [SiteFrequency(1, frequencies)],
        [SiteLikelihood(1, -14.2296, -14.3580)],
        [SiteLikelihood(1, -13.8521, -13.9077)],
        "LG", "C20",
    )

    assert [row.model for row in rows] == ["LG", "C20"]
    assert rows[0].keff == pytest.approx(5.0)
    assert rows[0].delta_lnl_tree2_tree1 == pytest.approx(-0.1284)
    assert rows[1].delta_lnl_tree2_tree1 == pytest.approx(-0.0556)


def test_bin_summary_fills_one_through_twenty_and_keff_20_stays_in_20() -> None:
    rows = [
        CcaRow("LG", 1, 1.9, -1, 1, 2),
        CcaRow("LG", 2, 20.0, -1, 2, 3),
        CcaRow("C20", 1, 1.9, -1, 0, 1),
    ]
    sums = summarize_bins(rows, ("LG", "C20"))

    assert len(sums["LG"]) == len(sums["C20"]) == 20
    assert sums["LG"][0] == 2
    assert sums["LG"][19] == 3
    assert sums["C20"][19] == 0
```

- [ ] **Step 2: Run calculation tests and verify failure**

Run: `pytest tests/posttree/test_syserror_cca.py -q`

Expected: FAIL because the calculation/bin interfaces are undefined.

- [ ] **Step 3: Implement the minimum calculation interfaces**

```python
def _keff(frequencies: tuple[float, ...]) -> float:
    return 1.0 / sum(value * value for value in frequencies)


def build_cca_rows(site_freq, model1, model2, model1_name, model2_name) -> list[CcaRow]:
    rows = []
    for freq, first, second in zip(site_freq, model1, model2):
        keff = _keff(freq.frequencies)
        for name, likelihood in ((model1_name, first), (model2_name, second)):
            rows.append(CcaRow(name, freq.site, keff, likelihood.lnl_tree1, likelihood.lnl_tree2, likelihood.lnl_tree2 - likelihood.lnl_tree1))
    return rows


def _keff_bin(keff: float) -> int:
    bin_id = 20 if math.isclose(keff, 20.0, abs_tol=1e-9) else math.floor(keff)
    if not 1 <= bin_id <= 20:
        raise ValueError(f"Keff is outside the supported range [1, 20]: {keff}")
    return bin_id
```

`_keff_bin()` must handle an upper-bound value within `1e-9` of 20 as bin 20, then explicitly reject a bin outside `1..20` with `ValueError`; this is defensive even though validated frequency vectors mathematically constrain Keff to that interval. `summarize_bins()` returns exactly 20 numeric values per supplied model in bins 1–20. Validate both names are non-empty after stripping and are distinct before rows are built.

- [ ] **Step 4: Write failing PDF semantics tests**

```python
from matplotlib import colors

from phyloai.posttree.syserror_cca import plot_cca


def test_plot_writes_pdf_with_documented_colors_and_legend(tmp_path: Path, monkeypatch) -> None:
    captured = {}
    original_savefig = __import__("matplotlib.pyplot", fromlist=["savefig"]).savefig

    def capture(path, *args, **kwargs):
        captured["figure"] = __import__("matplotlib.pyplot", fromlist=["gcf"]).gcf()
        return original_savefig(path, *args, **kwargs)

    monkeypatch.setattr("matplotlib.pyplot.savefig", capture)
    output = tmp_path / "cca.pdf"
    plot_cca({"LG": [0.0] * 20, "C20": [1.0] + [0.0] * 19}, ("LG", "C20"), "", "X", "Y", 10, 6, 300, 16, output)

    axis = captured["figure"].axes[0]
    assert output.read_bytes().startswith(b"%PDF")
    bar_colors = {colors.to_hex(bar.get_facecolor()) for container in axis.containers for bar in container}
    assert {"#f8766d", "#00bfc4"} <= bar_colors
    assert axis.get_legend().get_title().get_text() == ""
```

- [ ] **Step 5: Implement the Matplotlib plot helper**

Use only installed Matplotlib. Establish an approximately `theme_bw()`-like baseline with white axes, black spines, default-sized axis text, no major x grid, a horizontal zero line, vertical boundaries at integers 1–20, and `ax.margins(x=0)`. Compute the exact documented y limits, including a safe non-zero fallback when all bin sums are zero so Matplotlib never has identical limits.

For two models, emulate `geom_bar(width=1, position_dodge(width=1))`: category centres are `bin + 0.5`; each bar is width `0.5`, offset by `-0.25`/`+0.25`, so their combined group width is one. Use colors in model argument order: `#F8766D`, then `#00BFC4`. Draw the orange `#ffdab9` alpha-0.5 region from `ymin` to zero and a light-blue alpha-0.5 region from zero to `ymax`, behind bars.

Configure `ax.legend(title=None, loc="upper right", bbox_to_anchor=(0.99, 0.9), borderaxespad=0, framealpha=0.5, facecolor="white", edgecolor="black", fontsize=font_size)` and set its frame linewidth to `0.5`. Respect title/xlabel/ylabel/figure size/DPI, call `tight_layout()`, write only the requested PDF via `matplotlib.pyplot.savefig()`, then close the figure.

- [ ] **Step 6: Add numerical and visual regression coverage, then pass it**

Add the fixed fixture anchor using frequencies that yield `11.974845235298696`, and assert both reference deltas: `LG=-0.1284`, `C20=-0.0556`. Add plot assertions for 20-bin x ticks, x limits, `#00BFC4` second bar color, legend font size, figure size, and non-empty single-sign/all-zero cases. Run:

`pytest tests/posttree/test_syserror_cca.py -q`

Expected: PASS.

- [ ] **Step 7: Review Task 2; do not commit**

Run: `git diff --check && git diff -- phyloai/posttree/syserror_cca.py tests/posttree/test_syserror_cca.py`

Expected: no whitespace errors.

### Task 3: Runner, Stable Files, And Result JSON

**Files:**
- Modify: `phyloai/posttree/syserror_cca.py`
- Modify: `tests/posttree/test_syserror_cca.py`

**Interfaces:**
- Produces `build_cca_command()` for the fully resolved invocation and `run_cca()` with the named parameters/defaults in the approved command interface.

- [ ] **Step 1: Write failing runner/output tests**

```python
import csv
import json

from phyloai.posttree.syserror_cca import run_cca


def test_runner_writes_only_cca_csv_pdf_and_result(tmp_path: Path) -> None:
    freq, lnl1, lnl2 = write_cca_inputs(tmp_path)  # local test helper writes two valid sites
    output = tmp_path / "out"

    payload = run_cca(freq, lnl1, lnl2, "LG", "C20", output_dir=output, quiet=True)

    assert {path.name for path in output.iterdir()} == {"cca.csv", "cca.pdf", "result.json"}
    rows = list(csv.DictReader((output / "cca.csv").open()))
    assert list(rows[0]) == ["model", "site", "keff", "lnl_tree1", "lnl_tree2", "delta_lnl_tree2_tree1"]
    assert [row["model"] for row in rows] == ["LG", "C20", "LG", "C20"]
    assert payload["key_results"]["n_sites"] == 2
    persisted = json.loads((output / "result.json").read_text())
    assert set(persisted["data"]["output_files"]) == {"cca_table", "cca_figure"}


def test_dry_run_validates_but_writes_nothing(tmp_path: Path) -> None:
    freq, lnl1, lnl2 = write_cca_inputs(tmp_path)
    output = tmp_path / "out"

    result = run_cca(freq, lnl1, lnl2, "LG", "C20", output_dir=output, dry_run=True, quiet=True)

    assert result["data"]["output_files"] == {}
    assert not output.exists()
```

- [ ] **Step 2: Run runner tests and verify failure**

Run: `pytest tests/posttree/test_syserror_cca.py -q`

Expected: FAIL because `run_cca` is undefined.

- [ ] **Step 3: Implement `run_cca()` in validated-write order**

```python
def run_cca(
    site_freq: Path, site_lnl1: Path, site_lnl2: Path,
    model1_name: str = "model1", model2_name: str = "model2",
    title: str = "", xlabel: str = "Effective number of amino acids",
    ylabel: str = "Log-likelihood difference", fig_width: float = 10,
    fig_height: float = 6, dpi: int = 300, font_size: float = 16,
    output_dir: Path = Path("runs/posttree/syserror/cca"), overwrite: bool = False,
    dry_run: bool = False, quiet: bool = False,
) -> dict[str, Any]:
    site_frequencies = parse_site_freq(site_freq)
    first, second = parse_site_lnl(site_lnl1), parse_site_lnl(site_lnl2)
    validate_matching_sites(site_frequencies, first, second)
    rows = build_cca_rows(site_frequencies, first, second, model1_name, model2_name)
    # Construct payload after validating names/plot values; write only after dry-run return.
```

Validation/write sequence:

1. Parse all three inputs, validate same sites, validate names and all positive plot dimensions/DPI/font size, calculate rows/bins/key results.
2. Resolve output path; construct fully reproducible command and complete `params`.
3. Return the success payload immediately for dry run without creating a directory or listing future outputs.
4. Reject a file output path; reject a non-empty directory unless overwrite; remove a non-empty directory only after all prior validation succeeds; create it.
5. Write `cca.csv` with `csv.DictWriter`, exact field order, ascending-site/model1-then-model2 ordering, four fixed decimal places for `lnl_tree1`, `lnl_tree2`, and `delta_lnl_tree2_tree1`, and full floating-point precision for `keff`.
6. Write `cca.pdf`; set output records `cca_table` and `cca_figure` to resolved paths and scientific descriptions.
7. Set rounded elapsed time, call `write_result_json(payload, output_dir)`, print `cca.csv`, `cca.pdf`, then `result.json` unless quiet, and return payload.

Set `key_results` to `n_sites`, `keff_min`, `keff_max`, `keff_mean`, `models` (the two names), `total_delta_lnl_tree2_tree1` keyed by model, and `bin_summaries` keyed by model with 20 `{bin, delta_lnl_tree2_tree1}` records. Set `tool_versions={}` and an empty local command/stderr/warnings structure plus output records.

- [ ] **Step 4: Add lifecycle, validation, and JSON tests**

Add tests that reject blank/equal model labels and zero/negative figure dimensions, reject an output path that is a file without touching it, preserve a non-empty output directory without overwrite, replace it after successful overwrite, and preserve all non-`result.json` files when invalid input occurs with overwrite while replacing root `result.json` with the error record. Assert complete shell-quoted command reconstruction for paths containing spaces and nondefault names/plot settings. Assert result JSON has all standard top-level fields, resolved parameters, absolute output paths, 20 bin summary records for each model, `tool_versions == {}`, and empty local-only command/stderr fields.

Run: `pytest tests/posttree/test_syserror_cca.py -q`

Expected: PASS.

- [ ] **Step 5: Review Task 3; do not commit**

Run: `git diff --check && git diff --stat`

Expected: no whitespace errors.

### Task 4: Click Command, Error Results, Dynamic MCP, And Report Methods

**Files:**
- Modify: `phyloai/cli/commands/posttree.py`
- Create: `tests/cli/test_posttree_syserror_cca.py`
- Modify: `phyloai/mcp/tools/stubs.py`
- Modify: `phyloai/report/templates.py`
- Modify: `tests/report/test_templates.py`

**Interfaces:**
- Consumes `run_cca()`/`build_cca_command()` from Task 3 plus `_write_error_result_json()` and `walk_click_tree()`.
- Produces Click/MCP leaf `posttree_syserror_cca` and the updated `generate_methods_posttree_syserror_cca()`.

- [ ] **Step 1: Write failing command/discovery tests**

```python
import json

from click.testing import CliRunner

from phyloai.cli.main import cli
from phyloai.mcp.schema_gen import walk_click_tree
from phyloai.mcp.tools.stubs import STUB_TOOL_NAMES


def test_cca_help_and_generated_mcp_leaf() -> None:
    result = CliRunner().invoke(cli, ["posttree", "syserror", "cca", "--help"])

    assert result.exit_code == 0
    for option in ("--site-freq", "--site-lnl1", "--site-lnl2", "--model1-name", "--fig-width", "--dry-run"):
        assert option in result.output
    assert "site, lnL_Tree1, and lnL_Tree2" in result.output
    assert "posttree_syserror_cca" in {item["tool_name"] for item in walk_click_tree(cli)}
    assert "posttree_syserror_cca" not in STUB_TOOL_NAMES


def test_cca_cli_writes_error_result_but_dry_run_does_not(tmp_path) -> None:
    output = tmp_path / "out"
    failed = CliRunner().invoke(cli, [
        "posttree", "syserror", "cca", "--site-freq", str(tmp_path / "missing"),
        "--site-lnl1", str(tmp_path / "one.csv"), "--site-lnl2", str(tmp_path / "two.csv"), "-o", str(output),
    ])
    assert failed.exit_code == 1
    assert (output / "result.json").exists()

    dry = CliRunner().invoke(cli, [
        "posttree", "syserror", "cca", "--site-freq", str(tmp_path / "missing"),
        "--site-lnl1", str(tmp_path / "one.csv"), "--site-lnl2", str(tmp_path / "two.csv"),
        "-o", str(tmp_path / "dry"), "--dry-run",
    ])
    assert dry.exit_code == 1
    assert not (tmp_path / "dry").exists()
```

- [ ] **Step 2: Run integration tests and verify failure**

Run: `pytest tests/cli/test_posttree_syserror_cca.py tests/report/test_templates.py -q`

Expected: FAIL because `cca` is not registered and the MCP name remains a stub.

- [ ] **Step 3: Register a thin grouped Click command**

Add a `_GroupedCcaCommand` subclass whose `option_sections` split `site_freq`/`site_lnl1`/`site_lnl2` into **Required Inputs**, model labels into **Model Labels**, figure controls into **Figure Options**, and lifecycle flags into **Common Options**. Register `@syserror.command("cca", cls=_GroupedCcaCommand)` with required `click.Path(path_type=Path)` inputs and no `exists=True`; use `click.IntRange(1)` for DPI and `click.FloatRange(min=0, min_open=True)` for figure/font values. Include all approved defaults and alias/help text.

Update `_SyserrorGroup.list_commands()` and its subcommand listing to include `cca`. Its detailed command docstring/help must state accepted `.sitefreq` provenance/format; required exact LNL headers; the deliberate Tree2-minus-Tree1 sign; two model-specific LNL tables; fixed outputs; no external tools; and these examples:

```bash
phyloai posttree syserror cca \
  --site-freq chain1.sitefreq \
  --site-lnl1 lnl_LG/site_lnl.csv --site-lnl2 lnl_C20/site_lnl.csv \
  --model1-name LG --model2-name C20

phyloai posttree syserror cca --site-freq chain1.sitefreq \
  --site-lnl1 lnl1/site_lnl.csv --site-lnl2 lnl2/site_lnl.csv --dry-run
```

Build the error command using `build_cca_command()` and call `run_cca()`. Before invoking it, reject a file output path cleanly and record whether a directory was non-empty. On `ValueError`, write standard input-error JSON only when not dry-run and the directory was not already non-empty without overwrite; retain existing files when an overwrite run fails validation. Print JSON only on successful dry run.

- [ ] **Step 4: Remove stub and replace the report placeholder**

Remove `posttree_syserror_cca` and its description from `phyloai/mcp/tools/stubs.py`, leaving `posttree_syserror_sites` unchanged. Do not write MCP schema/handler code: dynamic discovery must be the only implementation.

Replace the placeholder CCA method generator with defensive prose such as:

```python
def generate_methods_posttree_syserror_cca(params, key_results, tool_versions) -> str:
    models = key_results.get("models", [params.get("model1_name", "model1"), params.get("model2_name", "model2")])
    n_sites = key_results.get("n_sites", 0)
    return (
        "Compositional constraint analysis (CCA) was performed using PhyloAI across "
        f"{_describe_n(n_sites, 'alignment site', 'alignment sites')}. "
        "The effective number of amino acids (Keff) was calculated as the inverse "
        "homozygosity of the 20 site-specific amino-acid frequencies, and site-wise "
        "Tree2-minus-Tree1 log-likelihood differences were summed within floor(Keff) bins "
        f"for the {models[0]} and {models[1]} model analyses."
    )
```

Normalize malformed/missing `models` safely to the two parameter defaults rather than allowing report generation to fail. Do not claim that a model or topology is preferred.

- [ ] **Step 5: Add lifecycle/report tests and pass focused suites**

Add a successful CLI invocation using tiny generated inputs; assert `cca.csv`, `cca.pdf`, result file registration, and dynamic discovery. Add directory-input cases for each required input, model-label/figure Click validation errors, complete error-command reconstruction, and invalid overwrite preservation (`keep.txt` survives while root error JSON is written). Add report template tests asserting `Keff`, `floor(Keff)`, both names, site count, and `Tree2-minus-Tree1` appear; assert malformed model data does not raise.

Run: `pytest tests/cli/test_posttree_syserror_cca.py tests/report/test_templates.py -q`

Expected: PASS.

- [ ] **Step 6: Review Task 4; do not commit**

Run: `git diff --check && git diff --stat`

Expected: no whitespace errors, and no authored MCP tool schema.

### Task 5: User Documentation, Parent Design, MCP Guide, And Workflow Skill

**Files:**
- Create: `docs/commands/posttree-syserror-cca.md`
- Create: `docs/commands/posttree-syserror-cca.zh.md`
- Modify: `README.md`
- Modify: `README.zh.md`
- Modify: `docs/superpowers/specs/2026-06-07-phyloai-design.md`
- Modify: `docs/commands/ai-integration.md`
- Modify: `docs/commands/ai-integration.zh.md`
- Modify: `skills/phyloai-workflow/SKILL.md`
- Modify: `skills/phyloai-workflow/references/parameter-annotations.md`

**Interfaces:**
- Consumes final Click help, standard result payload, and report methods behavior from Tasks 3–4.
- Produces consistent bilingual public documentation and safe CCA workflow guidance.

- [ ] **Step 1: Create the English CCA reference**

Create `docs/commands/posttree-syserror-cca.md`, linking to its Chinese counterpart, with the established sections **Purpose**, **Usage**, **Inputs**, **Calculation And CSV**, **Figure**, **Outputs**, **Examples**, **Warnings / Errors**, and **Notes**. Specify every option/default; `.sitefreq` only; 1-based exact site matching; 20 frequencies/sum tolerance; required literal LNL headers and ignored `ΔSLS`; Keff formula; reverse sign; exact CSV field schema and training-header mapping; bins 1–20; all fixed plotting semantics/colors; exactly three outputs; dry-run/overwrite lifecycle; and no doctor/resume. State clearly that the figure is diagnostic and cannot establish the true topology/model.

- [ ] **Step 2: Create matching Chinese reference and update READMEs**

Create the Chinese page with identical behavioral coverage. Add a concise CCA example after rate in both READMEs:

```bash
# Compositional-constraint diagnostic across two model analyses
phyloai posttree syserror cca --site-freq chain1.sitefreq \
  --site-lnl1 lnl_LG/site_lnl.csv --site-lnl2 lnl_C20/site_lnl.csv \
  --model1-name LG --model2-name C20 -o runs/posttree/syserror/cca
```

Add language-matched command-table entries linking to the CCA pages. The description must call CCA a composition-constraint/systematic-error diagnostic, not a correction or model-selection result.

- [ ] **Step 3: Correct parent-design and MCP user documentation**

In `docs/superpowers/specs/2026-06-07-phyloai-design.md` Section 4.1, replace only the obsolete CCA line using `--matrix --t1 --t2` with the approved prepared-input invocation. In Section 11 Phase 9, retain the atom sequence `brlen → rate → cca → sites` but clarify CCA consumes a prepared `.sitefreq` plus two `site_lnl.csv` inputs after site-frequency/site-likelihood generation. Update Section 8 only if needed to keep its future-Skill wording consistent; do not claim the future `phyloai-syserror` Skill has been implemented.

Move `posttree_syserror_cca` from the stub sentence to the fire-and-forget tool table in both AI-integration pages, and change the remaining stub list to only unavailable tools. Preserve the documentation’s explicit dynamic-schema/no-manual-sync architecture statement.

- [ ] **Step 4: Add workflow and parameter-card guidance**

In `skills/phyloai-workflow/SKILL.md`, insert a `### posttree syserror cca` section next to rate: it is local-only and needs no doctor, but still requires `get_command_schema`, a complete parameter card, and explicit approval before execution. Explain the three inputs, exact LNL headers, prepared-source provenance, two-model labels, calculation/sign, and that interpretation must not declare a preferred topology/model. Update the global local-only approval paragraph to include CCA.

Add Chinese annotations for `--site-freq`, `--site-lnl1`, `--site-lnl2`, `--model1-name`, `--model2-name`, `--title`, `--xlabel`, `--ylabel`, `--fig-width`, `--fig-height`, `--dpi`, `--font-size`, `--output-dir`, `--overwrite`, `--dry-run`, and `--quiet`. State the required columns, sign convention, defaults, only-PDF output, and destructive nature of overwrite. Do not omit visual controls: the Skill core requires a card for every runtime-schema parameter.

- [ ] **Step 5: Inspect documentation references and review; do not commit**

Run: `rg -n 'syserror cca|posttree_syserror_cca|--matrix ./matrix.fa --t1' README.md README.zh.md docs/commands skills/phyloai-workflow docs/superpowers/specs/2026-06-07-phyloai-design.md`

Expected: CCA appears in the README, CCA docs, MCP table, Skill, annotations, and corrected parent design; no obsolete CCA interface and no CCA stub listing remain.

Run: `git diff --check`

Expected: no whitespace errors.

### Task 6: Integrated Verification And Specification Status

**Files:**
- Modify: `docs/superpowers/specs/2026-08-13-phyloai-posttree-syserror-cca-design.md`

**Interfaces:**
- Consumes all completed behavior and documentation from Tasks 1–5.
- Produces an accurately marked implemented specification only after verification evidence exists.

- [ ] **Step 1: Verify real CLI help and generated MCP tool definitions**

Run:

```bash
python -m phyloai.cli.main posttree syserror cca --help
python - <<'PY'
from phyloai.mcp.tools.cli_tools import get_tool_definitions
from phyloai.mcp.tools.stubs import STUB_TOOL_NAMES
assert "posttree_syserror_cca" in get_tool_definitions()
assert "posttree_syserror_cca" not in STUB_TOOL_NAMES
PY
```

Expected: help exits 0 and documents all three required inputs; MCP exposes the generated CCA tool and no longer exposes a CCA stub.

- [ ] **Step 2: Run targeted verification**

Run:

```bash
pytest tests/posttree/test_syserror_cca.py \
  tests/cli/test_posttree_syserror_cca.py \
  tests/report/test_templates.py -q
```

Expected: PASS.

- [ ] **Step 3: Run the complete test suite**

Run: `pytest -q`

Expected: PASS. If it fails, use `superpowers:systematic-debugging` before changing code; do not weaken unrelated assertions.

- [ ] **Step 4: Mark the CCA design implemented only after successful checks**

Change the CCA design document status from `Proposed — awaiting user review` to `Implemented` and include the actual implementation date. Do not modify this status if any required test remains failing.

- [ ] **Step 5: Final diff and scope review; do not commit**

Run:

```bash
git diff --check
git status --short
git diff --stat
```

Expected: no whitespace errors; changed files match this plan; no accidental dependency/lockfile, R code, PNG output, generic abstraction, or unapproved commit.

## Final Acceptance Checklist

- [ ] `.sitefreq` and both LNL tables are strictly validated as matching one-based consecutive site sets; `ΔSLS` is ignored.
- [ ] Keff matches `11.974845235298696`; historical training `cca.txt` and current bundled likelihood fixtures are separately labeled rather than treated as one numerical fixture.
- [ ] `cca.csv` uses its exact snake-case schema, stable order, supplied/default labels, and no aggregate duplicate table.
- [ ] PDF binning/fill/color/layout/y-padding/legend semantics match the documented R plotting contract; only a PDF is written.
- [ ] `result.json`, dry-run, output conflicts, overwrite validation safety, error result behavior, and local-only fields comply with project standards.
- [ ] Click help is detailed and grouped; dynamic MCP exposes CCA; its stub is removed; report methods and generic figure/table indexing work.
- [ ] Parent design, bilingual command/MCP documentation, bilingual READMEs, and workflow Skill/annotations all match the final interface.
- [ ] Targeted tests and `pytest -q` pass. No commit occurs without explicit user approval.
