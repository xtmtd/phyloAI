# Taxon Composition Heterogeneity (`taxcomp`) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add `phyloai posttree syserror taxcomp`, a pure-Python diagnostic for across-taxon compositional heterogeneity, and remove the obsolete unimplemented `syserror sites` public placeholder.

**Architecture:** A small private helper extracted from `simulate_adequacy.py` computes the already-defined PPA-COMP taxon-frequency distances. A new `syserror_taxcomp.py` builds the Pearson taxon-by-state contingency table, overall and row statistics, Holm-adjusted nominal p-values, sparse-cell status, CSV/TSV summaries, and standard `result.json`. The Click command exposes this stable interface; MCP discovers it dynamically. Report collection, methods generation, user documentation, and the workflow Skill describe it as exploratory and remove the unimplemented `sites` interface.

**Tech Stack:** Python 3, BioPython `MultipleSeqAlignment`, SciPy `scipy.stats.chi2`, standard library (`csv`, `math`, `shlex`, `shutil`, `time`), Click, pytest.

**Spec:** `docs/superpowers/specs/2026-08-25-phyloai-posttree-syserror-taxcomp-design.md`

## Global Constraints

- Implement only the approved `taxcomp` scope and the removal of the never-implemented `syserror sites` placeholder.
- Use no new dependency; SciPy and BioPython are already installed project dependencies.
- `taxcomp` is a pure-Python atomic diagnostic. It invokes no external tool, needs no `doctor`, threads, checkpoint, or resume option.
- Accept only FASTA, PHYLIP, PHYLIP-PAML, and Nexus through existing `FormatConverter`; do not add Clustal support.
- Resolve `--seq-type auto|AA|NT`; count only `ACDEFGHIKLMNPQRSTVWY` for AA or `ACGT` for NT. Treat every other character as missing and do not assign fractional ambiguity counts.
- Overall Pearson p-values, per-taxon p-values, and `p_holm` are nominal exploratory values. Never classify a taxon as significant, failed, outlier, biased, or removable.
- `sparse_count_check` has only `triggered` and `not_triggered`; it is not an assumption pass/fail result.
- Preserve the existing public behavior of `posttree simulate adequacy`.
- Write `result.json` for persistent successful output and eligible claimed-directory errors. Successful payloads set `error_category: null`; dry runs have `data.output_files: {}`.
- Do not commit without separate user approval. Use focused tests after every task and the complete affected suite before requesting approval.

---

## File Map

| File | Change | Responsibility |
|---|---|---|
| `phyloai/posttree/simulate_adequacy.py` | Modify | Extract the private, reusable PPA-COMP taxon-composition helper without changing adequacy outputs. |
| `phyloai/posttree/syserror_taxcomp.py` | Create | Alignment validation, Pearson and Holm computations, output writing, command construction, and `run_taxcomp()`. |
| `phyloai/cli/commands/posttree.py` | Modify | Add the (ungrouped) Click `taxcomp` leaf command and correct CLI error-result behavior. |
| `tests/posttree/test_simulate_adequacy.py` | Modify | Lock the extracted helper to current PPA-COMP outputs. |
| `tests/posttree/test_syserror_taxcomp.py` | Create | Unit and output-lifecycle tests for the new module. |
| `tests/cli/test_posttree_syserror_taxcomp.py` | Create | Click help, dry run, result JSON, and dynamic MCP discovery tests. |
| `phyloai/report/collector.py` | Modify | Register `posttree.syserror.taxcomp` and remove `posttree.syserror.sites` from report parsing/order. |
| `phyloai/report/templates.py` | Modify | Add the constrained taxcomp methods generator; remove the sites placeholder generator. |
| `tests/report/test_collector.py` | Modify | Assert taxcomp parsing/order and sites removal. |
| `tests/report/test_templates.py` | Modify | Assert taxcomp methods wording and registry coverage without sites. |
| `phyloai/mcp/tools/stubs.py` | Modify | Remove `posttree_syserror_sites`; leave an empty, typed stub collection if no stubs remain. |
| `tests/mcp/test_schema_gen.py` | Modify | Assert generated `posttree_syserror_taxcomp` schema. |
| `tests/mcp/test_stubs.py` | Modify | Allow an empty stub registry and assert the removed sites tool is no longer advertised or handled. |
| `tests/mcp/test_cli_tools.py` | Modify | Remove the assumption that at least one MCP stub exists. |
| `README.md`, `README.zh.md` | Modify | Add taxcomp to the systematic-error command index/example in English and Chinese. |
| `docs/commands/posttree-syserror-taxcomp.md` | Create | English user command documentation. |
| `docs/commands/posttree-syserror-taxcomp.zh.md` | Create | Chinese user command documentation. |
| `docs/commands/ai-integration.md`, `docs/commands/ai-integration.zh.md` | Modify | Remove the obsolete sites MCP stub from public tool documentation. |
| `skills/phyloai-workflow/SKILL.md` | Modify | Add taxcomp workflow and interpretation guidance; replace remaining sites references. |
| `skills/phyloai-workflow/references/parameter-annotations.md` | Modify | Add Chinese annotations for the stable taxcomp parameters. |
| `docs/superpowers/specs/2026-06-07-phyloai-design.md` | Already modified | Do not edit during implementation unless a necessary contract discrepancy is found and separately approved. |
| `docs/superpowers/specs/2026-08-25-phyloai-posttree-syserror-taxcomp-design.md` | Already modified | Treat as the binding contract; do not alter silently during implementation. |

### Task 1: Extract the private PPA-COMP helper

**Files:**
- Modify: `phyloai/posttree/simulate_adequacy.py:508-566`
- Modify: `tests/posttree/test_simulate_adequacy.py`

**Interfaces:**
- Consumes: `MultipleSeqAlignment`, resolved `seq_type: str`, `AA_STATES`, and `NT_STATES`.
- Produces: `_compute_taxon_composition(alignment: MultipleSeqAlignment, seq_type: str) -> dict[str, Any]` with exactly `taxon_freqs`, `taxon_dist_j`, `comp_max`, and `comp_mean`.
- Preserves: `_compute_statistics()` keeps returning its current keys and values, including `taxon_dist_j`, `comp_max`, and `comp_mean`.

- [ ] **Step 1: Write the failing helper identity test**

Add this test to `tests/posttree/test_simulate_adequacy.py`:

```python
def test_taxon_composition_helper_matches_existing_statistics() -> None:
    alignment = _msa(
        ("A", "AA--"),
        ("B", "AC--"),
        ("C", "CC--"),
        ("D", "CA--"),
    )

    comp = _compute_taxon_composition(alignment, "AA")
    stats = _compute_statistics(alignment, "AA")

    assert comp["taxon_dist_j"] == pytest.approx(stats["taxon_dist_j"])
    assert comp["comp_max"] == pytest.approx(stats["comp_max"])
    assert comp["comp_mean"] == pytest.approx(stats["comp_mean"])
    assert comp["taxon_freqs"]["A"] == pytest.approx([1.0, 0.0] + [0.0] * 18)
```

Import `_compute_taxon_composition` in the test module.

- [ ] **Step 2: Run the focused test and verify it fails**

Run:

```bash
pytest tests/posttree/test_simulate_adequacy.py::test_taxon_composition_helper_matches_existing_statistics -v
```

Expected: FAIL during test collection because `_compute_taxon_composition` does not exist.

- [ ] **Step 3: Extract the minimal private helper**

In `phyloai/posttree/simulate_adequacy.py`, add this helper immediately before `_compute_statistics`:

```python
def _compute_taxon_composition(
    alignment: MultipleSeqAlignment,
    seq_type: str,
) -> dict[str, Any]:
    states = AA_STATES if seq_type == "AA" else NT_STATES
    state_index = {state: index for index, state in enumerate(states)}
    taxon_freqs: dict[str, list[float]] = {}

    for record in alignment:
        counts = [0] * len(states)
        for state in str(record.seq).upper():
            index = state_index.get(state)
            if index is not None:
                counts[index] += 1
        total = sum(counts)
        if not total:
            raise ValueError(f"taxon {record.id!r} has no valid characters")
        taxon_freqs[record.id] = [count / total for count in counts]

    global_freq = [
        sum(freq[index] for freq in taxon_freqs.values()) / len(taxon_freqs)
        for index in range(len(states))
    ]
    taxon_dist = {
        name: sum((freq[index] - global_freq[index]) ** 2 for index in range(len(states)))
        for name, freq in taxon_freqs.items()
    }
    return {
        "taxon_freqs": taxon_freqs,
        "taxon_dist_j": taxon_dist,
        "comp_max": max(taxon_dist.values()),
        "comp_mean": sum(taxon_dist.values()) / len(taxon_dist),
    }
```

Replace only the taxon-frequency/global-frequency/distance block at the end of `_compute_statistics()` with:

```python
    composition = _compute_taxon_composition(alignment, seq_type)
    return {
        "div": diversity_total / n_informative,
        "siteconvprob": squared_freq_total / n_informative,
        "sitecomp": sitecomp,
        "n_informative_sites": n_informative,
        **composition,
    }
```

Do not change AA/NT character sets, site-statistic calculations, or external function signatures.

- [ ] **Step 4: Run focused adequacy tests**

Run:

```bash
pytest tests/posttree/test_simulate_adequacy.py -v
```

Expected: PASS. Existing numerical PPA-COMP tests and the new identity test preserve the old outputs.

- [ ] **Step 5: Inspect the focused diff**

Run:

```bash
git diff --check
git diff -- phyloai/posttree/simulate_adequacy.py tests/posttree/test_simulate_adequacy.py
```

Expected: only the private helper extraction and its regression test; no public adequacy schema change.

### Task 2: Implement and test pure taxcomp statistics

**Files:**
- Create: `phyloai/posttree/syserror_taxcomp.py`
- Create: `tests/posttree/test_syserror_taxcomp.py`

**Interfaces:**
- Consumes: `Path` matrix input, `seq_type: str`, `table_format: str`, standard `FormatConverter`, `_compute_taxon_composition`, and `scipy.stats.chi2`.
- Produces:

```python
def compute_taxcomp_statistics(
    alignment: MultipleSeqAlignment,
    seq_type: str,
) -> tuple[dict[str, Any], list[dict[str, Any]]]: ...

def holm_adjust(p_values: list[float]) -> list[float]: ...

def sparse_count_check(expected_cells: list[float]) -> dict[str, int | float | str]: ...

def build_taxcomp_command(
    matrix: Path, seq_type: str, table_format: str, output_dir: Path,
    overwrite: bool, dry_run: bool, quiet: bool,
) -> str: ...

def run_taxcomp(
    matrix: Path, seq_type: str = "auto", table_format: str = "csv",
    output_dir: Path = Path("runs/posttree/syserror/taxcomp"),
    overwrite: bool = False, dry_run: bool = False, quiet: bool = False,
) -> dict[str, Any]: ...
```

- [ ] **Step 1: Write failing statistical unit tests**

Create `tests/posttree/test_syserror_taxcomp.py` with small direct tests. Use this 3-taxon NT alignment as the primary hand-check fixture:

```python
def _msa(*rows: tuple[str, str]) -> MultipleSeqAlignment:
    return MultipleSeqAlignment([
        SeqRecord(Seq(sequence), id=taxon, description="")
        for taxon, sequence in rows
    ])


def test_compute_taxcomp_statistics_known_nt_table() -> None:
    overall, rows = compute_taxcomp_statistics(
        _msa(("A", "AAAA"), ("B", "CCCC"), ("C", "AACC")),
        "NT",
    )

    assert overall["n_taxa"] == 3
    assert overall["n_states"] == 2
    assert overall["df"] == 2
    assert overall["x2"] == pytest.approx(8.0)
    assert sum(row["x2_contribution"] for row in rows) == pytest.approx(overall["x2"])
    assert [row["df"] for row in rows] == [1, 1, 1]
    assert [row["taxon"] for row in rows] == ["A", "B", "C"]


def test_holm_adjusts_in_original_order_and_is_monotone() -> None:
    adjusted = holm_adjust([0.03, 0.001, 0.02, 0.20])

    assert adjusted == pytest.approx([0.06, 0.004, 0.06, 0.20])


def test_sparse_count_check_has_strict_boundaries() -> None:
    assert sparse_count_check([1.0, 5.0, 5.0, 5.0, 5.0]) == {
        "sparse_count_check": "not_triggered",
        "expected_cells_total": 5,
        "expected_cells_below_1": 0,
        "expected_cells_below_5": 1,
        "expected_cells_below_5_fraction": 0.2,
    }
```

Add separate tests for: one expected cell `0.999` triggers; six of twenty-five expected cells below 5 triggers; globally absent states reduce `n_states`; gaps/ambiguities are ignored; all-missing taxon, duplicate taxon, fewer than two taxa, and fewer than two global states raise `ValueError`.

- [ ] **Step 2: Run the unit tests and verify they fail**

Run:

```bash
pytest tests/posttree/test_syserror_taxcomp.py -v
```

Expected: FAIL during collection because `phyloai.posttree.syserror_taxcomp` does not exist.

- [ ] **Step 3: Implement parsing and deterministic statistics**

Create `phyloai/posttree/syserror_taxcomp.py`. Reuse `FormatConverter` and the adequacy module’s `_resolved_seq_type` logic only if doing so does not create a circular import; otherwise place a small local `_resolve_seq_type()` using `phyloai.core.sequence_normalization.detect_seq_type` and the same AA/NT validation.

Implement these rules exactly:

```python
states = AA_STATES if resolved_seq_type == "AA" else NT_STATES
state_index = {state: index for index, state in enumerate(states)}
counts = [[0] * len(states) for _ in alignment]
```

- Uppercase each sequence and increment a cell only for a standard state.
- Require unique identifiers, aligned sequence lengths, at least two taxa, and a positive standard-character row total for every taxon.
- Remove globally zero state columns **only** for Pearson expected counts and df.
- Build `E_ij = row_total[i] * column_total[j] / grand_total`.
- Set `x2_contribution` to each row’s sum of `(observed - expected) ** 2 / expected`; set overall `x2` to the sum of row contributions.
- Calculate `p_nominal` with `chi2.sf(x2, df)` for the overall result and `chi2.sf(row_x2, n_states - 1)` for each row. Return Python finite floats.
- Implement Holm by sorting `(p_value, original_index)`. For the `i`th smallest p-value using one-based `i`, multiply by `T - i + 1`; apply a cumulative maximum, clip at `1.0`, then restore original order. Equal p-values must receive equal adjusted values.
- Implement the strict sparse-cell rule from the spec and return all five required fields.
- Call `_compute_taxon_composition()` using the full resolved alphabet. Add each `squared_composition_distance` to the taxon row and `comp_max`/`comp_mean` to the overall dict.

Use these exact stable output columns:

```python
OVERALL_FIELDS = [
    "n_taxa", "n_states", "x2", "df", "p_nominal",
    "sparse_count_check", "expected_cells_total", "expected_cells_below_1",
    "expected_cells_below_5", "expected_cells_below_5_fraction",
    "comp_max", "comp_mean",
]
TAXON_FIELDS = [
    "taxon", "x2_contribution", "df", "p_nominal", "p_holm",
    "squared_composition_distance",
]
```

- [ ] **Step 4: Run the unit tests and fix only the implementation defects they expose**

Run:

```bash
pytest tests/posttree/test_syserror_taxcomp.py -v
```

Expected: PASS for deterministic statistics and validation tests.

- [ ] **Step 5: Add output lifecycle tests before implementing `run_taxcomp()`**

Add failing tests that call `run_taxcomp()` with a tiny FASTA fixture and assert:

```python
def test_run_taxcomp_writes_only_tables_and_result_json(tmp_path: Path) -> None:
    matrix = tmp_path / "matrix.fa"
    matrix.write_text(">A\nAAAA\n>B\nCCCC\n>C\nAACC\n")
    output = tmp_path / "out"

    payload = run_taxcomp(matrix, seq_type="NT", output_dir=output, quiet=True)

    assert {path.name for path in output.iterdir()} == {
        "overall_summary.csv", "taxon_summary.csv", "result.json",
    }
    assert payload["error_category"] is None
    assert set(payload["data"]["output_files"]) == {"overall_summary", "taxon_summary"}
    assert payload["key_results"]["seq_type"] == "NT"


def test_run_taxcomp_dry_run_writes_no_files(tmp_path: Path) -> None:
    matrix = tmp_path / "matrix.fa"
    matrix.write_text(">A\nAAAA\n>B\nCCCC\n>C\nAACC\n")

    payload = run_taxcomp(matrix, seq_type="NT", output_dir=tmp_path / "out", dry_run=True, quiet=True)

    assert payload["data"]["output_files"] == {}
    assert not (tmp_path / "out").exists()
```

Add tests for TSV suffix/delimiter, non-empty output refusal without overwrite, overwrite replacement after valid preflight, and invalid input preserving existing output content.

- [ ] **Step 6: Implement result payload and output lifecycle**

Implement `build_taxcomp_command()` with all resolved values and exact CLI flag names. `run_taxcomp()` must:

1. Validate `table_format` and `seq_type` before claiming the output directory.
2. Validate `matrix.exists()` and `matrix.is_file()` before claiming the output directory.
3. Parse/validate the alignment and compute summaries before deleting or creating the output directory.
4. On `dry_run`, return a success payload with all computed `key_results`, `error: None`, `error_category: None`, `data.cmd: []`, `data.tool_stderr: ""`, `data.warnings: []`, `data.output_files: {}`, and no filesystem writes.
5. For a persistent run, refuse a non-empty directory unless `overwrite`; with `overwrite`, remove it only after successful preflight and statistics calculation.
6. Write `overall_summary.<csv|tsv>`, `taxon_summary.<csv|tsv>`, then `result.json` through `write_result_json()`.
7. Include full `params`, `key_results` fields, `data.character_policy`, and absolute output paths exactly as specified.
8. Add a warning only when `sparse_count_check == "triggered"`; the warning must say the nominal chi-square p-values are especially unreliable due to sparse expected counts and must not call any taxon failed or significant.

Use `csv.DictWriter` with `extrasaction="raise"`, standard newline handling, and `delimiter = "\t" if table_format == "tsv" else ","`.

- [ ] **Step 7: Run the new module test suite**

Run:

```bash
pytest tests/posttree/test_syserror_taxcomp.py tests/posttree/test_simulate_adequacy.py -v
```

Expected: PASS. The new taxcomp output behavior and previous adequacy outputs both remain correct.

### Task 3: Add Click command, dynamic MCP coverage, and error-result behavior

**Files:**
- Modify: `phyloai/cli/commands/posttree.py`
- Create: `tests/cli/test_posttree_syserror_taxcomp.py`
- Modify: `tests/mcp/test_schema_gen.py`

**Interfaces:**
- Consumes: `run_taxcomp()` and `build_taxcomp_command()` from Task 2.
- Produces: `phyloai posttree syserror taxcomp` and dynamically discovered MCP tool `posttree_syserror_taxcomp`.

- [ ] **Step 1: Write failing CLI and MCP tests**

Create `tests/cli/test_posttree_syserror_taxcomp.py` with tests for:

```python
def test_taxcomp_help_and_generated_mcp_leaf() -> None:
    result = CliRunner().invoke(cli, ["posttree", "syserror", "taxcomp", "--help"])

    assert result.exit_code == 0
    for option in ("--matrix", "--seq-type", "--table-format", "--overwrite", "--dry-run", "--quiet"):
        assert option in result.output
    assert "exploratory" in result.output.lower()
    assert "sparse-cell" in result.output.lower()
    assert "posttree_syserror_taxcomp" in {item["tool_name"] for item in walk_click_tree(cli)}
```

Add a CLI success test that verifies the two tables and `result.json`; a dry-run test that verifies no output directory; and an invalid-matrix test that verifies error `result.json` contains `status: "error"` and `error_category: "input"` when the directory can be claimed.

Add this test to `tests/mcp/test_schema_gen.py`:

```python
def test_taxcomp_mcp_tool_is_generated_from_click() -> None:
    descriptor = next(
        item for item in walk_click_tree(cli)
        if item["tool_name"] == "posttree_syserror_taxcomp"
    )
    tool = build_mcp_tool(descriptor)

    assert tool["inputSchema"]["required"] == ["matrix"]
    props = tool["inputSchema"]["properties"]
    assert props["seq_type"]["enum"] == ["AA", "NT", "auto"]
    assert props["table_format"]["enum"] == ["csv", "tsv"]
```

- [ ] **Step 2: Run the focused CLI/MCP tests and verify they fail**

Run:

```bash
pytest tests/cli/test_posttree_syserror_taxcomp.py tests/mcp/test_schema_gen.py -v
```

Expected: taxcomp tests FAIL because the Click leaf does not exist; existing MCP schema tests remain green.

- [ ] **Step 3: Register the (ungrouped) Click command**

In `phyloai/cli/commands/posttree.py`, register `taxcomp` as a **plain**
command (no `_GroupedCommand` grouping). The option set is small and simple,
so the `Required Inputs` / `Analysis Options` / `Common Options` sectioned-help
model was deliberately dropped on user approval; register it as
`@syserror.command("taxcomp")` (no `cls=`), which renders options flat.
Register it so it appears before `cca`. Use these decorators:

```python
@click.option("--matrix", type=click.Path(path_type=Path), required=True,
              help="Aligned FASTA, PHYLIP, PHYLIP-PAML, or Nexus MSA.")
@click.option("--seq-type", type=_CaseInsensitiveChoice(["AA", "NT", "auto"]),
              default="auto", show_default=True,
              help="Sequence type. Standard AA/NT states are counted; gaps and ambiguity codes are excluded.")
@click.option("--table-format", type=click.Choice(["csv", "tsv"]),
              default="csv", show_default=True,
              help="Delimiter and suffix for overall_summary and taxon_summary.")
@click.option("-o", "--output-dir", type=click.Path(path_type=Path),
              default=Path("runs/posttree/syserror/taxcomp"), show_default=True,
              help="Output directory containing two summaries and result.json.")
@click.option("--overwrite", is_flag=True, default=False,
              help="Delete and recreate a non-empty output directory after validation succeeds.")
@click.option("--dry-run", is_flag=True, default=False,
              help="Validate and calculate summaries without writing files.")
@click.option("-q", "--quiet", is_flag=True, default=False,
              help="Suppress terminal output except errors.")
```

The command docstring must include the spec’s help acceptance baseline: exploratory only; no taxon removal/recode/model/topology decision; p-values are nominal; the sparse-cell status is not an assumption pass; PPA-COMP distances have no universal cutoff and simulation is required for model adequacy.

Build the reproducible command with `build_taxcomp_command()`, call `run_taxcomp()`, print dry-run JSON, and render concise normal output from returned `key_results`. On `ValueError`, use `_write_error_result_json()` with `error_category="input"` only when the target can safely be claimed, mirroring CCA’s existing protected-output behavior.

- [ ] **Step 4: Verify shared CLI error payload category behavior**

Do not change `_write_error_result_json()`: it already serializes the provided top-level `error_category`. Add an assertion in the taxcomp invalid-input CLI test that its eligible error result contains:

```python
{"status": "error", "error_category": "input"}
```

This protects the established shared behavior while avoiding an unrelated no-op code change.

- [ ] **Step 5: Run CLI/MCP tests**

Run:

```bash
pytest tests/cli/test_posttree_syserror_taxcomp.py tests/cli/test_posttree_syserror_cca.py tests/mcp/test_schema_gen.py -v
```

Expected: PASS. CCA retains behavior, taxcomp is discoverable dynamically, and no manual MCP tool is added.

### Task 4: Integrate reports and remove `sites` placeholder code

**Files:**
- Modify: `phyloai/report/collector.py`
- Modify: `phyloai/report/templates.py`
- Modify: `phyloai/mcp/tools/stubs.py`
- Modify: `tests/report/test_collector.py`
- Modify: `tests/report/test_templates.py`
- Modify: `tests/mcp/test_stubs.py`
- Modify: `tests/mcp/test_cli_tools.py`

**Interfaces:**
- Consumes: result commands `phyloai posttree syserror taxcomp ...`; `key_results` fields from Task 2.
- Produces: report step ID `posttree.syserror.taxcomp`; no remaining public `posttree.syserror.sites` report/MCP placeholder.

- [ ] **Step 1: Write failing report and stub-removal tests**

Add tests that assert:

```python
def test_taxcomp_step_id_and_order() -> None:
    assert parse_step_id("phyloai posttree syserror taxcomp --matrix matrix.fa") == "posttree.syserror.taxcomp"
    assert "posttree.syserror.taxcomp" in STEP_ORDER
    assert "posttree.syserror.sites" not in STEP_ORDER
    assert STEP_ORDER.index("posttree.syserror.taxcomp") < STEP_ORDER.index("posttree.syserror.brlen")
```

Add a report-template test that calls `generate_all_methods("posttree.syserror.taxcomp", ...)` with a complete minimal payload and asserts the generated text contains all of:

```python
("Pearson", "Holm", "nominal", "sparse-cell", "phylogenetically", "comp_max", "comp_mean", "No taxon")
```

Replace the existing non-empty-stub assumption with tests that allow the
registry to be empty and assert:

```python
assert STUB_TOOL_NAMES == frozenset()
assert STUB_TOOLS == []
assert "posttree_syserror_sites" not in STUB_TOOL_NAMES
assert handle_stub("posttree_syserror_sites") is None
assert handle_stub("pretree_align") is None
```

In `tests/mcp/test_cli_tools.py`, replace any `next(iter(STUB_TOOL_NAMES))`
launch assertion with a check that dynamic CLI handlers still build when the
stub collection is empty and that `posttree_syserror_sites` is absent from the
handler mapping.

- [ ] **Step 2: Run the focused tests and verify they fail**

Run:

```bash
pytest tests/report/test_collector.py tests/report/test_templates.py tests/mcp/test_stubs.py tests/mcp/test_cli_tools.py -v
```

Expected: new taxcomp/removed-sites assertions FAIL.

- [ ] **Step 3: Update report step discovery**

In `phyloai/report/collector.py`:

- Insert `"posttree.syserror.taxcomp"` before `"posttree.syserror.brlen"` in `STEP_ORDER`.
- Remove `"posttree.syserror.sites"` from `STEP_ORDER`.
- Change `_THIRD_LEVEL["syserror"]` to exactly `{"taxcomp", "brlen", "rate", "cca"}`.

Do not change the generic parsing algorithm.

- [ ] **Step 4: Add the constrained taxcomp methods generator and remove sites generator**

In `phyloai/report/templates.py`, delete `generate_methods_posttree_syserror_sites()` and its registry entry. Add and register `generate_methods_posttree_syserror_taxcomp()`.

The function must interpolate `n_taxa`, `n_sites`, and `seq_type`, but retain this substance when values are absent:

```python
return (
    f"Taxon compositional heterogeneity was screened in a {n_taxa}-taxon, "
    f"{n_sites}-site {seq_type} alignment using a Pearson common-composition X2 "
    "statistic calculated from the taxon-by-state count table. Overall X2 and "
    "each taxon row contribution were reported with nominal chi-square reference "
    "p-values; per-taxon values were additionally adjusted by Holm's method for "
    "multiplicity, but remain exploratory because homologous taxa are "
    "phylogenetically dependent. A conventional sparse-cell rule was evaluated "
    "from expected counts; this concerns the asymptotic chi-square reference only "
    "and does not validate the phylogenetic null model. Per-taxon unitless squared "
    "Euclidean composition discrepancies from the equal-taxon mean, together with "
    "their maximum (comp_max) and mean (comp_mean), summarized observed taxon "
    "heterogeneity. These distances have no universal cutoff and require simulated "
    "comparison for model adequacy assessment. No taxon was removed, data recoded, "
    "or topology/model selected automatically."
)
```

Do not add threshold counts or labels such as significant, failed, outlier, biased, or remove.

- [ ] **Step 5: Remove the MCP sites stub**

In `phyloai/mcp/tools/stubs.py`, replace the obsolete populated structures with:

```python
STUB_TOOL_NAMES: frozenset[str] = frozenset()
_DESCRIPTIONS: dict[str, str] = {}
STUB_TOOLS: list[dict] = []
```

Keep `handle_stub()` unchanged except that it naturally returns `None` for every name. This is the smallest deletion compatible with `cli_tools.py` importing `STUB_TOOLS` and `handle_stub`.

- [ ] **Step 6: Run report and MCP stub tests**

Run:

```bash
pytest tests/report/test_collector.py tests/report/test_templates.py tests/mcp/test_stubs.py tests/mcp/test_cli_tools.py tests/mcp/test_schema_gen.py -v
```

Expected: PASS. Report generation knows taxcomp, the old sites step is absent, and the only taxcomp MCP tool comes from Click discovery.

### Task 5: Add end-user documentation and workflow guidance

**Files:**
- Create: `docs/commands/posttree-syserror-taxcomp.md`
- Create: `docs/commands/posttree-syserror-taxcomp.zh.md`
- Modify: `README.md`
- Modify: `README.zh.md`
- Modify: `docs/commands/ai-integration.md`
- Modify: `docs/commands/ai-integration.zh.md`
- Modify: `skills/phyloai-workflow/SKILL.md`
- Modify: `skills/phyloai-workflow/references/parameter-annotations.md`

**Interfaces:**
- Consumes: the stable CLI surface and result schema from Tasks 2-3.
- Produces: English/Chinese user documentation and user-safe workflow guidance.

- [ ] **Step 1: Write documentation assertions before documenting**

Add a lightweight documentation test in `tests/cli/test_posttree_syserror_taxcomp.py` that reads both command docs after they are created and asserts each contains the command path, `p_holm`, `sparse_count_check`, `comp_max`, `comp_mean`, and the phrase or Chinese equivalent that no automatic taxon removal occurs.

- [ ] **Step 2: Run the documentation assertion and verify it fails**

Run:

```bash
pytest tests/cli/test_posttree_syserror_taxcomp.py::test_taxcomp_command_docs_cover_interpretation_boundaries -v
```

Expected: FAIL because both command documents do not yet exist.

- [ ] **Step 3: Write both command documents from the stable contract**

Create the English and Chinese command documents with matching sections:

1. Purpose and explicit non-goals.
2. Usage, including one AA and one NT example.
3. Input formats and valid-character/missing-character policy.
4. Exact output filenames and columns.
5. Interpretation of overall X2, row contribution, `p_nominal`, `p_holm`, `sparse_count_check`, squared composition distance, `comp_max`, and `comp_mean`.
6. Clear caveats: no p-value is phylogenetically calibrated; `not_triggered` is not an assumption pass; distances are not evolutionary distances; no universal threshold.
7. Sensitivity follow-up: AA can use `phyloai pretree concat --recoding Dayhoff-6`; NT can use `phyloai pretree concat --recoding RY-nucleotide`; both must be run and reported beside original data. Taxon deletion requires independent curation evidence and is not produced by this command.
8. Error/output-directory behavior and dry-run behavior.

Use these executable examples:

```bash
phyloai posttree syserror taxcomp --matrix matrix.aa.fa --seq-type AA
phyloai posttree syserror taxcomp --matrix matrix.nt.fa --seq-type NT --table-format tsv -o runs/posttree/syserror/taxcomp-nt
```

Do not mention Clustal input, `dayhoff6`, BH/FDR, a fixed distance threshold, automatic recoding, automatic taxon removal, or `syserror sites`.

- [ ] **Step 4: Update indexes and MCP documentation**

In both README files, add `taxcomp` in the systematic-error command list with a short description: across-taxon composition screening using Pearson and PPA-COMP observed summaries. Do not describe it as a model-adequacy test.

In both AI integration command documents:

- Add `posttree_syserror_taxcomp` to the generated execution-tools table after the taxcomp CLI exists.
- Remove `posttree_syserror_sites` from the stub-tool sentence.
- Do not claim manual MCP registration; explain that the CLI leaf is dynamically discovered where that is already the local documentation pattern.

- [ ] **Step 5: Update the PhyloAI workflow Skill and annotations**

In `skills/phyloai-workflow/SKILL.md`:

- Add a `### posttree syserror taxcomp` section before `brlen`.
- State: pure Python/no doctor; one aligned MSA; requires full schema card and explicit execution approval; no tree/model required.
- Give a results interpretation sequence: inspect `sparse_count_check` first; treat p-values and `p_holm` as nominal exploratory values; use taxon row contribution and squared distance to prioritize inspection; use `simulate adequacy` for model-calibrated PPA-COMP; optional recoding is a separate, user-approved sensitivity analysis with `Dayhoff-6` for AA or `RY-nucleotide` for NT.
- Replace remaining `sites` references in the brlen/CCA context with `taxcomp`, `rate`, `cca`, or `simulate adequacy` as scientifically appropriate.

In `skills/phyloai-workflow/references/parameter-annotations.md`, add Chinese entries for `taxcomp` parameters, including why `--seq-type`, `--table-format`, `--output-dir`, and destructive `--overwrite` matter. Do not invent an annotation for an option not in the Click schema.

- [ ] **Step 6: Run documentation and targeted regression tests**

Run:

```bash
pytest tests/cli/test_posttree_syserror_taxcomp.py tests/report/test_collector.py tests/report/test_templates.py tests/mcp/test_stubs.py tests/mcp/test_cli_tools.py -v
rg -n "dayhoff6|posttree_syserror_sites|posttree\.syserror\.sites|syserror sites" README.md README.zh.md docs/commands skills/phyloai-workflow phyloai tests
rg -n "combine with .*sites diagnostics" skills/phyloai-workflow/SKILL.md && exit 1 || true
```

Expected: all pytest tests PASS. The first `rg` command returns no matches for removed `sites` names and wrong `dayhoff6`; the second returns no match, specifically guarding the previously bare `sites diagnostics` wording. Unrelated ordinary uses of the English word `sites` are not a failure.

### Task 6: Full verification and review handoff

**Files:**
- Review only: all files modified by Tasks 1-5.

**Interfaces:**
- Consumes: completed code, tests, CLI docs, and Skill updates.
- Produces: evidence for user approval to commit; no commit is made in this task.

- [ ] **Step 1: Run static and focused test suites**

Run:

```bash
python -m compileall -q phyloai
pytest tests/posttree/test_simulate_adequacy.py tests/posttree/test_syserror_taxcomp.py tests/cli/test_posttree_syserror_taxcomp.py tests/cli/test_posttree_syserror_cca.py tests/report/test_collector.py tests/report/test_templates.py tests/mcp/test_schema_gen.py tests/mcp/test_stubs.py tests/mcp/test_cli_tools.py -v
```

Expected: all commands succeed and every listed test passes.

- [ ] **Step 2: Run the command-level smoke checks**

Create a temporary three-taxon NT FASTA and run:

```bash
phyloai posttree syserror taxcomp --matrix /tmp/taxcomp.nt.fa --seq-type NT -o /tmp/taxcomp-out --overwrite
phyloai posttree syserror taxcomp --matrix /tmp/taxcomp.nt.fa --seq-type NT -o /tmp/taxcomp-dry --dry-run
python - <<'PY'
import json
from pathlib import Path
payload = json.loads(Path('/tmp/taxcomp-out/result.json').read_text())
assert payload['status'] == 'success'
assert payload['error_category'] is None
assert set(payload['data']['output_files']) == {'overall_summary', 'taxon_summary'}
assert payload['key_results']['seq_type'] == 'NT'
assert not Path('/tmp/taxcomp-dry').exists()
PY
```

Expected: persistent run produces exactly two summary tables plus `result.json`; dry run produces no directory; returned values are finite JSON-compatible values.

- [ ] **Step 3: Inspect interface discoverability and removed placeholder**

Run:

```bash
phyloai posttree syserror taxcomp --help
python - <<'PY'
from phyloai.cli.main import cli
from phyloai.mcp.schema_gen import walk_click_tree
from phyloai.mcp.tools.stubs import STUB_TOOL_NAMES
names = {item['tool_name'] for item in walk_click_tree(cli)}
assert 'posttree_syserror_taxcomp' in names
assert 'posttree_syserror_sites' not in names
assert 'posttree_syserror_sites' not in STUB_TOOL_NAMES
PY
```

Expected: help contains all six options and the interpretation boundary; taxcomp is dynamically discoverable and sites is neither a CLI-derived tool nor a stub.

- [ ] **Step 4: Review diff and report exact evidence**

Run:

```bash
git diff --check
git status --short
git diff --stat
```

Expected: no whitespace errors, only files named in this plan changed, and no unapproved unrelated change is reverted.

- [ ] **Step 5: Request code review and commit approval**

Summarize the changed files, command smoke results, and test results. State explicitly that no commit has been made. Ask the user for separate approval before creating one Conventional Commit such as:

```text
feat: add taxon composition diagnostics
```

## Plan Self-Review

### Spec coverage

- Command surface, supported formats, no external tool, and output lifecycle: Tasks 2-3.
- Shared PPA-COMP identity, full alphabet vs Pearson effective `K`, missing-character policy, X2, row contributions, nominal p-values, Holm, sparse-cell boundaries, and output schemas: Tasks 1-2.
- Detailed help and report scientific-language baseline: Tasks 3-5.
- Dynamic MCP integration: Task 3 and Task 6.
- Report parsing/template integration: Task 4.
- Removal of unimplemented `sites` from MCP, report, docs, Skill, and tests: Tasks 4-5.
- English/Chinese docs and correct AA/NT recoding guidance: Task 5.
- Regression, smoke, discoverability, whitespace, and diff evidence: Task 6.

### Placeholder scan

The plan contains no implementation placeholders. Every code task names concrete files, symbols, schema fields, tests, and commands. The only future-facing instruction is the explicit user approval gate before a commit.

### Type consistency

- Task 1 defines `_compute_taxon_composition()` used by Task 2.
- Task 2 defines `compute_taxcomp_statistics()`, `holm_adjust()`, `sparse_count_check()`, `build_taxcomp_command()`, and `run_taxcomp()` used by Tasks 3-6.
- Task 3 defines the stable Click leaf and dynamic MCP name used by Tasks 4-6.
- Task 4 defines the report step ID and generator used by Task 5 documentation and Task 6 verification.

## Execution Handoff

Plan complete and saved to `docs/superpowers/plans/2026-08-25-phyloai-posttree-syserror-taxcomp.md`.

Two execution options:

1. **Subagent-Driven (recommended)** - Dispatch a fresh subagent per task, review between tasks, and preserve the approval gate before commit.
2. **Inline Execution** - Execute tasks in this session using `executing-plans`, with checkpoints after each task.

Choose one approach after reviewing the plan.
