# Posttree AliSim Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Implement `phyloai posttree simulate alisim params`, `iqtree`, and `transfergaps` for empirical AliSim parameter extraction, reproducible single/batch MSA simulation, and single-file or directory-batch gap-mask transfer.

**Architecture:** Keep the three responsibilities in three focused `phyloai/posttree/` modules. The params and transfergaps commands are pure Python; iqtree delegates each simulation to IQ-TREE3, reuses shared checkpoint, executable-resolution, file-matching, validation, and result-json infrastructure, and records every sampled parameter for reproducibility. The Click tree, generated MCP schemas, report templates, documentation, and workflow skill remain thin adapters over these library entry points.

**Tech Stack:** Python 3.11+, Click, Rich, Biopython, NumPy, Matplotlib, IQ-TREE3, pytest.

## Global Constraints

- Follow `docs/superpowers/specs/2026-08-02-phyloai-posttree-simulate-alisim-design.md` exactly; this plan implements only the AliSim group, not `adequacy` or `phybase`.
- All non-`doctor` commands write one root `result.json` with full resolved parameters, full re-executable command, absolute `data.output_files` paths, and standard exit codes.
- Use kebab-case CLI options; use `--seq-type` in PhyloAI and translate it to IQ-TREE `--seqtype`.
- Reuse `phyloai.core.file_matching`, `phyloai.core.checkpoint`, `phyloai.core.runner`, IQ-TREE resolution/version helpers, and shared sequence-output validation rather than recreating them.
- Batch simulation uses `ProcessPoolExecutor`, Rich remaining-work progress, atomic checkpoint writes, per-simulation seeds drawn from a master-seeded generator (see Addendum), and per-task merged stdout/stderr logs.
- IQ-TREE commands use `--out-format`, use `-p` for partition files, and reject managed I/O/model/output flags in `--tool-args`.
- All generated FASTA outputs wrap at 60 characters; generated MSAs are validated before success is recorded.
- No new pip dependencies. Primary platforms are Linux/macOS, secondary WSL, native Windows unsupported.
- Do not implement directory-wide gap transfer, model adequacy, gene-tree simulation, or a custom density/KDE framework.

---

## File Map

| Action | File | Responsibility |
|---|---|---|
| Create | `phyloai/posttree/simulate_alisim_params.py` | Parse IQ-TREE reports, pair trees, and write `params.tsv`/result data. |
| Create | `phyloai/posttree/simulate_alisim_iqtree.py` | Validate models/tables, sample rows, construct AliSim commands, execute/resume, validate outputs, and plot PDF-mode diagnostics. |
| Create | `phyloai/posttree/simulate_alisim_transfergaps.py` | Name-match taxa and replace simulated sites using an original alignment mask. |
| Modify | `phyloai/cli/commands/posttree.py` | Register simulate/alisim Click groups and all command options. |
| Modify | `phyloai/mcp/tools/stubs.py` | Remove the top-level `posttree_simulate` stub; retain future leaf stubs through the Click group. |
| Modify | `phyloai/report/templates.py` | Add params, single/batch AliSim, and transfergaps methods generators. |
| Modify | `phyloai/report/collector.py` | Add the three step IDs and nested command depth to report ordering. |
| Modify | `README.md`, `README.zh.md` | Add compact AliSim workflow examples. |
| Create | `docs/commands/posttree-simulate-alisim.md`, `docs/commands/posttree-simulate-alisim.zh.md` | Document inputs, outputs, examples, warnings, and recovery. |
| Modify | `skills/phyloai-workflow/SKILL.md` | Add approval-card inputs, execution/recovery guidance, and result interpretation. |
| Create | `tests/posttree/test_simulate_alisim_params.py` | Parameter parser and tree pairing tests. |
| Create | `tests/posttree/test_simulate_alisim_iqtree.py` | Model reconstruction, sampling, dry-run, checkpoint, and result structure tests. |
| Create | `tests/posttree/test_simulate_alisim_transfergaps.py` | Name-based gap-mask transfer tests. |
| Create | `tests/cli/test_posttree_simulate_alisim.py` | Click hierarchy, option validation, and wrapper tests. |

## Task 1: Parameter Extraction Library

**Files:**
- Create: `phyloai/posttree/simulate_alisim_params.py`
- Create: `tests/posttree/test_simulate_alisim_params.py`

**Interfaces:**
- Consumes: `Path iqtree_dir`, `Path tree_dir`, `Path output_dir`, `overwrite`, `dry_run`, `quiet`.
- Produces: `parse_iqtree_report(path: Path) -> dict[str, str]`, `run_alisim_params(...) -> dict[str, Any]`, and a UTF-8 tab-delimited `params.tsv` with the 11 spec columns.

^- [x] **Step 1: Write failing parser and matching tests**

```python
def test_parse_gtr_f_i_g_and_tree_pair(tmp_path: Path) -> None:
    report = tmp_path / "gene.iqtree"
    report.write_text("Input data: DNA\nTo simulate an alignment of the same\n"
                      'iqtree3 --alisim sim -m "GTR{1,2,3,4,5}+F{.1,.2,.3,.4}+I{.2}+G4{.7}" --length 100\n')
    parsed = parse_iqtree_report(report)
    assert parsed["seqtype"] == "DNA"
    assert parsed["subs_rate"] == "1/2/3/4/5"
    assert parsed["rate_param"] == ".7"

def test_run_skips_unmatched_and_rejects_ambiguous_tree(tmp_path: Path) -> None:
    result = run_alisim_params(iqtree_dir=tmp_path / "reports", tree_dir=tmp_path / "trees", output_dir=tmp_path / "out")
    assert result["data"]["unmatched"]

def test_parse_legal_empty_model_components(tmp_path: Path) -> None:
    report = tmp_path / "aa.iqtree"
    report.write_text("Input data: amino-acid\nTo simulate an alignment of the same\n"
                      'iqtree3 --alisim sim -m "LG" --length 100\n')
    parsed = parse_iqtree_report(report)
    assert parsed["subs_rate"] == parsed["freq"] == ""
    assert parsed["rate_heterogeneity"] == parsed["rate_categories"] == parsed["rate_param"] == ""
```

^- [x] **Step 2: Run tests to verify they fail**

Run: `pytest tests/posttree/test_simulate_alisim_params.py -v`

Expected: FAIL because the module and entry point do not exist.

^- [x] **Step 3: Implement the minimal parser and result writer**

```python
PARAM_COLUMNS = ("id", "seqtype", "length", "subs_model", "subs_rate", "freq", "prop_inv", "rate_heterogeneity", "rate_categories", "rate_param", "tree_path")

def parse_iqtree_report(path: Path) -> dict[str, str]:
    """Return AliSim-compatible values extracted from one IQ-TREE report."""
    # Locate Input data and the line after the simulation heading; normalize comma lists to '/'.

def run_alisim_params(*, iqtree_dir: Path, tree_dir: Path, output_dir: Path, overwrite: bool = False, dry_run: bool = False, quiet: bool = False) -> dict[str, Any]:
    """Write matched empirical parameter rows and the standard result payload."""
```

Implement all model forms in the spec: AA/DNA, optional `+F` with ordered `pi()` fallback, optional `+I`, Gamma, FreeRate, and no rate heterogeneity. Preserve legal empty fields exactly: AA `subs_rate`, absent `+F`, and absent `+G/+R` are empty strings, not parse errors. Scan reports recursively; use suffix-agnostic tree matching with an ambiguity error. Use `csv.DictWriter(..., delimiter="\t")`; do not write files on dry-run or create `checkpoint.json`.

^- [x] **Step 4: Run focused tests**

Run: `pytest tests/posttree/test_simulate_alisim_params.py -v`

Expected: PASS, including AA pi fallback, FreeRate, unmatched warning, ambiguity, and absolute `tree_path`/output paths.

^- [x] **Step 5: Commit the extraction unit**

```bash
git add phyloai/posttree/simulate_alisim_params.py tests/posttree/test_simulate_alisim_params.py
git commit -m "feat: extract AliSim parameters from IQ-TREE reports"
```

## Task 2: AliSim Sampling, Commands, and Execution Library

**Files:**
- Create: `phyloai/posttree/simulate_alisim_iqtree.py`
- Create: `tests/posttree/test_simulate_alisim_iqtree.py`

**Interfaces:**
- Consumes: either single-mode tree/model/length inputs or a `params.tsv`; shared `Checkpoint`; IQ-TREE executable.
- Produces: `build_model_string(row: Mapping[str, str]) -> str`, `sample_batch_rows(...) -> list[dict[str, str]]`, `run_alisim_iqtree(...) -> dict[str, Any]`, `MSAs/`, `logs/`, and batch `params_sampled.tsv`/`checkpoint.json`.

^- [x] **Step 1: Write failing model/sampling tests**

```python
def test_build_model_string_omits_absent_components() -> None:
    assert build_model_string({"subs_model": "LG", "subs_rate": "", "freq": "", "prop_inv": "", "rate_heterogeneity": "R", "rate_categories": "2", "rate_param": "0.5/1.5"}) == "LG+R2{0.5,1.5}"

def test_build_model_string_retains_zero_invariable_proportion() -> None:
    row = {"subs_model": "LG", "subs_rate": "", "freq": "", "prop_inv": "0", "rate_heterogeneity": "", "rate_categories": "", "rate_param": ""}
    assert build_model_string(row) == "LG+I{0}"

def test_mixed_sampling_keeps_model_core_and_rate_group_together() -> None:
    sampled = sample_batch_rows(rows, strategy="mixed", n=20, rng=random.Random(3), pdf_params=(), noise_scale=1.0, overrides={})
    assert all(row["seqtype"] == "AA" and row["subs_model"] == "LG" or row["seqtype"] == "DNA" and row["subs_model"] == "GTR" for row in sampled)

def test_pdf_preserves_i_presence_before_resampling_value() -> None:
    sampled = sample_batch_rows(rows, strategy="pdf", n=100, rng=random.Random(4), pdf_params=("prop_inv",), noise_scale=0.0, overrides={})
    assert any(row["prop_inv"] == "" for row in sampled)
    assert any(row["prop_inv"] for row in sampled)

def test_pdf_resamples_rate_param_only_for_gamma() -> None:
    sampled = sample_batch_rows(rows_with_g_r_and_none, strategy="pdf", n=100, rng=random.Random(5), pdf_params=("rate_param",), noise_scale=0.0, overrides={})
    assert all(row["rate_heterogeneity"] != "R" or row["rate_param"] in empirical_r_rate_params for row in sampled)
    assert all(row["rate_heterogeneity"] or row["rate_param"] == "" for row in sampled)
```

^- [x] **Step 2: Run tests to verify they fail**

Run: `pytest tests/posttree/test_simulate_alisim_iqtree.py -v`

Expected: FAIL because the module and helpers do not exist.

^- [x] **Step 3: Implement table validation, sampling, and command construction**

```python
def build_model_string(row: Mapping[str, str]) -> str:
    model = row["subs_model"]
    if row["subs_rate"]:
        model += "{" + row["subs_rate"].replace("/", ",") + "}"
    if row["freq"]:
        model += "+F{" + row["freq"].replace("/", ",") + "}"
    if row["prop_inv"]:
        model += "+I{" + row["prop_inv"] + "}"
    if row["rate_heterogeneity"]:
        model += f'+{row["rate_heterogeneity"]}{row["rate_categories"]}' + "{" + row["rate_param"].replace("/", ",") + "}"
    return model
```

Validate all required TSV columns and types. Sample complete rows as units; in mixed/pdf keep model core and rate group atomic, select `+I` presence before sampling/resampling its value, and only PDF-resample Gamma alpha after Gamma is selected. Treat only `""` as absent: every non-empty `prop_inv`, including `"0"` or `"0.0"`, must reconstruct as `+I{value}`. Accept only `length` and `prop_inv` overrides. Use Freedman-Diaconis bins with the specified jitter semantics, safe fallback for constant/singleton inputs, positive rounded lengths, and `prop_inv` clamped to `[0, 1)`.

^- [x] **Step 4: Write failing execution/resume tests**

```python
def test_single_dry_run_uses_out_format_and_partition_p(tmp_path: Path) -> None:
    result = run_alisim_iqtree(ref_tree=tree, model=None, model_partitions=parts, seq_type="AA", length=None, output_dir=tmp_path / "out", dry_run=True)
    assert "-p" in result["data"]["cmd"]
    assert "--out-format" in result["data"]["cmd"]

def test_batch_dry_run_records_deterministic_task_seeds(tmp_path: Path) -> None:
    result = run_alisim_iqtree(model_params=table, strategy="complete", num_simulations=2, seed=7, output_dir=tmp_path / "out", dry_run=True)
    assert [row["seed"] for row in result["data"]["sampled_rows"]] == [7, 8]
```

^- [x] **Step 5: Implement execution, output handling, and resume**

```python
def run_alisim_iqtree(*, ref_tree: Path | None = None, model: str | None = None, model_partitions: Path | None = None, seq_type: str | None = None, length: int | None = None, model_params: Path | None = None, strategy: str | None = None, num_simulations: int | None = None, override: str | None = None, noise_scale: float = 1.0, pdf_params: str = "length,prop_inv,rate_param", msa_prefix: str = "sim", out_format: str = "fasta", num_alignments: int = 1, iqtree_threads: int = 1, threads: int = 4, seed: int | None = None, iqtree_path: str | None = None, tool_args: str | None = None, output_dir: Path = Path("runs/posttree/simulate/alisim/iqtree"), overwrite: bool = False, resume: bool = False, dry_run: bool = False, quiet: bool = False) -> dict[str, Any]:
    """Run one AliSim invocation or a resumable batch of one-alignment invocations."""
```

Use the shared resolver/version checker, `Runner`, checkpoint helpers, process pool, and FASTA validator; add a PHYLIP validator through Biopython rather than treating PHYLIP as FASTA. Assign `task_index` from zero so master seed `7` produces seeds `7, 8, ...`. Write merged stdout/stderr to `logs/<simulation_id>.log`, validate before checkpoint success, and write every actual row to `params_sampled.tsv`. Generate density PDFs only for selected non-overridden PDF parameters. In single mode keep merged diagnostics in `data.tool_stderr`; in batch mode use `data.files[]` with `cmd`, `log_file`, positive `wall_time`, and `output_file`.

^- [x] **Step 6: Run focused tests**

Run: `pytest tests/posttree/test_simulate_alisim_iqtree.py -v`

Expected: PASS for model reconstruction, validation, all three strategies, invalid overrides, output conflict/resume mismatch, blocked tool arguments, dry-run argv, result shape, and PDF plot selection.

^- [x] **Step 7: Commit the simulation unit**

```bash
git add phyloai/posttree/simulate_alisim_iqtree.py tests/posttree/test_simulate_alisim_iqtree.py
git commit -m "feat: add AliSim single and batch simulation"
```

## Task 3: Single-File Gap-Mask Transfer Library

**Files:**
- Create: `phyloai/posttree/simulate_alisim_transfergaps.py`
- Create: `tests/posttree/test_simulate_alisim_transfergaps.py`

**Interfaces:**
- Consumes: an original aligned sequence file plus (exactly one of) a single simulated MSA (`--simulated-msa`) or a directory of simulated MSAs (`--simulated-dir`), and output-directory controls.
- Produces: `run_alisim_transfergaps(...) -> dict[str, Any]`; single mode writes `<original_stem>.gaps.fa`, batch mode writes `<simulated_stem>.gaps.fa` per input (output is always FASTA regardless of input format).

^- [x] **Step 1: Write failing taxon-name and ambiguity tests**

```python
def test_transfer_matches_taxa_by_name_and_preserves_original_order(tmp_path: Path) -> None:
    result = run_alisim_transfergaps(original_msa=original, simulated_msa=reordered, output_dir=tmp_path / "out")
    records = list(SeqIO.parse(result["data"]["output_files"]["transferred_msa"]["path"], "fasta"))
    assert [record.id for record in records] == ["A", "B"]
    assert str(records[0].seq) == "A-C-"

def test_exclude_ambiguity_masks_only_dash_and_dot(tmp_path: Path) -> None:
    result = run_alisim_transfergaps(original_msa=original_with_x, simulated_msa=simulated, exclude_ambiguity=True, output_dir=tmp_path / "out")
    assert result["key_results"]["n_positions_masked"] == 1

def test_dry_run_validates_but_writes_nothing(tmp_path: Path) -> None:
    output_dir = tmp_path / "out"
    result = run_alisim_transfergaps(original_msa=original, simulated_msa=simulated, output_dir=output_dir, dry_run=True)
    assert result["status"] == "success"
    assert not output_dir.exists()
```

^- [x] **Step 2: Run tests to verify they fail**

Run: `pytest tests/posttree/test_simulate_alisim_transfergaps.py -v`

Expected: FAIL because the module and entry point do not exist.

^- [x] **Step 3: Implement minimal validated transfer**

```python
def run_alisim_transfergaps(*, original_msa: Path, simulated_msa: Path | None = None, simulated_dir: Path | None = None, seq_type: str = "auto", exclude_ambiguity: bool = False, output_dir: Path = Path("runs/posttree/simulate/alisim/transfergaps"), overwrite: bool = False, dry_run: bool = False, quiet: bool = False) -> dict[str, Any]:
    """Replace selected simulated columns with the original per-taxon mask."""
```

Reject unparsable/empty inputs, duplicate taxon IDs, taxon-set mismatch, and unequal original/simulated lengths. Require exactly one of `simulated_msa` / `simulated_dir`. Inputs may be FASTA/PHYLIP/NEXUS/PHYLIP-PAML (read via the shared `FormatConverter`). Batch mode discovers alignment files by extension, transfers each independently, and names outputs `<simulated_stem>.gaps.fa` (always FASTA). Detect molecule type when requested. Replace, never insert, mask positions; default masks non-standard symbols, `--exclude-ambiguity` masks only `-` and `.`. Write 60-column FASTA and a standard result payload (with `n_msas`); dry-run performs validation but writes nothing.

^- [x] **Step 4: Run focused tests**

Run: `pytest tests/posttree/test_simulate_alisim_transfergaps.py -v`

Expected: PASS for taxon reordering, set mismatch, length mismatch, both mask policies, 60-column output, dry-run, and result JSON.

^- [x] **Step 5: Commit the transfer unit**

```bash
git add phyloai/posttree/simulate_alisim_transfergaps.py tests/posttree/test_simulate_alisim_transfergaps.py
git commit -m "feat: add AliSim gap mask transfer"
```

## Task 4: CLI, MCP, Report, and Tests

**Files:**
- Modify: `phyloai/cli/commands/posttree.py`
- Modify: `phyloai/mcp/tools/stubs.py`
- Modify: `phyloai/report/templates.py`
- Modify: `phyloai/report/collector.py`
- Create: `tests/cli/test_posttree_simulate_alisim.py`
- Modify: `tests/report/test_templates.py`

**Interfaces:**
- Consumes: the three `run_alisim_*` functions from Tasks 1-3.
- Produces: the `posttree simulate alisim` Click hierarchy, three generated MCP leaf tools, report step IDs, and deterministic methods text.

^- [x] **Step 1: Write failing CLI and template tests**

```python
def test_simulate_alisim_commands_are_registered() -> None:
    result = CliRunner().invoke(cli, ["posttree", "simulate", "alisim", "--help"])
    assert result.exit_code == 0
    assert {"params", "iqtree", "transfergaps"} <= set(result.output.split())

def test_alisim_batch_methods_mentions_strategy_and_failures() -> None:
    text = generate_all_methods("posttree.simulate.alisim.iqtree", {"model_params": "x.tsv", "strategy": "pdf", "pdf_params": "length", "noise_scale": 1.0}, {"source_loci": 4, "n_simulations_completed": 3, "n_simulations_failed": 1}, {"iqtree3": "3.1"})
    assert "histogram" in text.lower()
    assert "1 failed" in text
```

^- [x] **Step 2: Run tests to verify they fail**

Run: `pytest tests/cli/test_posttree_simulate_alisim.py tests/report/test_templates.py -v`

Expected: FAIL because the command tree and methods generators are not registered.

^- [x] **Step 3: Register command tree and result consumers**

```python
class _SimulateGroup(click.Group):
    def list_commands(self, ctx: click.Context) -> list[str]:
        return ["alisim", "adequacy", "phybase"]

class _AlisimGroup(click.Group):
    def list_commands(self, ctx: click.Context) -> list[str]:
        return ["params", "iqtree", "transfergaps"]
```

Add wrappers with every option and validation from the spec. Add future `adequacy`/`phybase` Click stubs returning not-implemented. The existing Click-tree schema walker will auto-generate their MCP leaf schemas and `launch_cli` will run their harmless not-implemented CLI response; do not add manual entries to `STUB_TOOL_NAMES`. Remove only the obsolete top-level `posttree_simulate` stub. Register report templates under `posttree.simulate.alisim.params`, `.iqtree`, and `.transfergaps`; add their collector ordering/depth. Preserve all existing posttree command names.

^- [x] **Step 4: Run integration tests**

Run: `pytest tests/cli/test_posttree_simulate_alisim.py tests/report/test_templates.py tests/mcp/test_stubs.py tests/mcp/test_cli_tools.py -v`

Expected: PASS; generated MCP schemas include the three implemented leaves and retain not-implemented future leaves.

^- [x] **Step 5: Commit integration**

```bash
git add phyloai/cli/commands/posttree.py phyloai/mcp/tools/stubs.py phyloai/report/templates.py phyloai/report/collector.py tests/cli/test_posttree_simulate_alisim.py tests/report/test_templates.py
git commit -m "feat: expose AliSim simulation commands"
```

## Task 5: User Documentation and Workflow Guidance

**Files:**
- Modify: `README.md`
- Modify: `README.zh.md`
- Create: `docs/commands/posttree-simulate-alisim.md`
- Create: `docs/commands/posttree-simulate-alisim.zh.md`
- Modify: `skills/phyloai-workflow/SKILL.md`

**Interfaces:**
- Consumes: final CLI option names and `result.json` fields from Tasks 1-4.
- Produces: concise user-facing workflow documentation and approval/recovery guidance.

^- [x] **Step 1: Write documentation assertions as literal command checks**

```text
The documentation must contain:
phyloai posttree simulate alisim params --iqtree-dir ... --tree-dir ...
phyloai posttree simulate alisim iqtree --model-params params.tsv --strategy complete --num-simulations 100
phyloai posttree simulate alisim transfergaps --original-msa original.fa --simulated-msa sim001.fa
```

^- [x] **Step 2: Add the minimal documentation**

Keep each README addition to a 3-5 line workflow example and link to the command document for detail, preserving the project's lightweight README policy. Document the three-command workflow, mutually exclusive modes, `-p` partition behavior, all sampling strategies, PDF limitations, `--override`, 1:1 gap transfer, output layouts, resume rules, and key report fields. In the workflow skill, require an explicit parameter review and user approval before invoking IQ-TREE execution; allow read-only result/report inspection directly.

^- [x] **Step 3: Verify documented commands and prose**

Run: `rg -n "posttree simulate alisim|--model-params|--exclude-ambiguity" README.md README.zh.md docs/commands/posttree-simulate-alisim.md docs/commands/posttree-simulate-alisim.zh.md skills/phyloai-workflow/SKILL.md`

Expected: Each command document and both README files contain the implemented command path; the skill contains the approval requirement.

^- [x] **Step 4: Commit documentation**

```bash
git add README.md README.zh.md docs/commands/posttree-simulate-alisim.md docs/commands/posttree-simulate-alisim.zh.md skills/phyloai-workflow/SKILL.md
git commit -m "docs: add AliSim simulation workflow"
```

## Task 6: Full Verification and Release Documentation

**Files:**
- Modify: `docs/superpowers/specs/2026-06-07-phyloai-design.md`
- Modify: `docs/superpowers/specs/2026-08-02-phyloai-posttree-simulate-alisim-design.md` only if implementation exposes a verified specification contradiction.

**Interfaces:**
- Consumes: completed modules and test suites.
- Produces: verified current-state documentation and clean test evidence.

^- [x] **Step 1: Run the full relevant suite**

Run: `pytest tests/posttree/test_simulate_alisim_params.py tests/posttree/test_simulate_alisim_iqtree.py tests/posttree/test_simulate_alisim_transfergaps.py tests/cli/test_posttree_simulate_alisim.py tests/report/test_templates.py tests/mcp/test_stubs.py tests/mcp/test_cli_tools.py -v`

Expected: PASS.

^- [x] **Step 2: Run CLI and MCP smoke checks**

Run: `phyloai posttree simulate alisim --help`

Expected: exit 0 and lists `params`, `iqtree`, and `transfergaps`.

Run: `phyloai posttree simulate adequacy`

Expected: explicit not-implemented response without a traceback.

^- [x] **Step 3: Confirm parent design remains precise**

Verify that the parent design describes only three implemented AliSim library modules and identifies `adequacy`/`phybase` as future stubs; do not claim they are implemented.

^- [x] **Step 4: Commit final design alignment if needed**

```bash
git add docs/superpowers/specs/2026-06-07-phyloai-design.md docs/superpowers/specs/2026-08-02-phyloai-posttree-simulate-alisim-design.md
git commit -m "docs: align AliSim design with implementation"
```

---

## Addendum (2026-08-03) — post-implementation corrections

Changes applied after initial implementation, per user review. All deviate from the original plan/spec text above and supersede it.

1. **`params.tsv` delimiter changed from `/` to `,`.** Multi-value fields (`subs_rate`, `freq`, `rate_param`) now store IQ-TREE's native comma delimiter, matching the `.iqtree` reports and the AliSim `-m` string. `build_model_string` still normalizes legacy `/` to `,` for backwards compatibility.
2. **result.json `command` is now complete** for all three commands (`params`, `iqtree`, `transfergaps`): it includes `-o`, `--overwrite`, `--strategy`, `--num-simulations`, `--threads`, `--msa-prefix`, `--out-format`, `--iqtree-threads`, `--seed`, `--tool-args`, `--resume`, `--dry-run`, `--quiet`, etc. as applicable.
3. **`--tool-args` blocked set narrowed to I/O-only** (`--alisim`, `-t`, `--prefix`, `--out-format`, `-af`). Non-I/O flags (e.g. `--seqtype`, `--length`, `--num-alignments`, `-T`, `-m`) are appended after PhyloAI's managed flags and may override them. Single-mode MSA collection now globs `<msa_prefix>*<ext>` so a `--num-alignments` override is honored.
4. **Single-mode logs.** AliSim does not write `{prefix}.iqtree`/`{prefix}.log` in the work dir (its console log goes to `<ref_tree>.log`). Single mode now captures Runner stdout/stderr into `logs/<msa_prefix>.log` (matching batch mode) and removes the stray `<ref_tree>.log` when it did not pre-exist. The `iqtree_report` output-file entry was dropped.
5. **`--override` fixed for `complete` strategy** (previously only applied in `mixed`/`pdf`).
6. **Per-simulation seeds are now independent random** values drawn from a master-seeded generator instead of `master_seed + task_index`. `--seed` is the master seed; the seed column in `params_sampled.tsv` no longer increments sequentially.
7. **`source_id` column in `params_sampled.tsv` is complete-strategy only**; omitted for `mixed`/`pdf`.
8. **Density plots are PDF-mode only** (no `plots/` for `complete`/`mixed`), restyled as Gaussian-KDE density curves with `#2E86AB` (empirical) / `#A23B72` (simulated), matching `ref/scripts/server.R`.
9. **`--strategy` defaults to `complete`** in the CLI; `--seq-type` and transfergaps `--seq-type` help now display `[AA|DNA]` / `[AA|NT|auto]` (case-insensitive input still accepted).
10. **`noise_scale`/`pdf_params` are `null` in result.json** for non-pdf strategies.
11. **Resume param hashing excludes `_command`** (`core/checkpoint.py`) so the re-executable command string (which now includes `--resume`/`--overwrite`) no longer causes false resume mismatches.

### Addendum 2 (2026-08-03, second review round)

12. **`--seq-type` help display.** The case-insensitive choice metavar is now produced by overriding `Choice.get_metavar` (Click renders `get_metavar` in help, not the `metavar` property). `--seq-type` displays `[AA|DNA]` / `[AA|NT|auto]` (uppercase, consistent with the parent design and other modules) while still accepting lowercase input.
13. **No duplicate managed flags in the IQ-TREE command.** `_build_alisim_cmd` now parses `--tool-args` for re-specified managed flags (`--seqtype`, `-m`, `-p`, `--length`, `--num-alignments`, `-T`, `--seed`, supporting `=` form) and suppresses its own copy, so each flag appears exactly once with the tool-args value. Previously PhyloAI emitted its flag and tool-args appended a duplicate (IQ-TREE happened to use the last).
14. **`result.json` command uses `shlex.join`** for `params`, `iqtree`, and `transfergaps`, so multi-token `--tool-args` values (and paths containing spaces) are shell-quoted and unambiguous.
15. **PDF density plots draw curves only** (no histogram, no fill): Gaussian-KDE lines in `#2E86AB` (empirical) and `#A23B72` (simulated), matching `ref/scripts/server.R` `geom_density`.

### Addendum 3 (2026-08-03, review findings)

16. **MCP `output_dir` default (cli_tools.py).** The hardcoded default-output-dir allowlist (four Bayesian tools only) was replaced with `_default_output_dir()`, which reads the click `output_dir` default from the tool's own schema. All AliSim tools now accept calls that omit `output_dir`, matching their advertised CLI defaults.
17. **Blocked `--tool-args` flags could be bypassed with `=` syntax** (`--out-format=phy` passed the exact-match check). `_check_managed_flag_conflict` now normalizes tokens with `token.split("=")[0]`, closing the redirect/duplicate-flag bypass.
18. **Doc example fix:** single-mode `--model-partitions` examples in `docs/commands/posttree-simulate-alisim.{md,zh.md}` dropped the spurious `--length` (partition mode infers it).
19. **Parent design alignment:** `2026-06-07-phyloai-design.md` gap-transfer row now states single-file and directory-batch modes (implementation supports both via `--simulated-msa`/`--simulated-dir`).
20. **Plan smoke checks** use `phyloai` instead of the nonexistent `python -m phyloai` entry point.
21. **Report methods text (templates.py + design §9.3):** the batch AliSim methods now explain each strategy — `pdf` is introduced as "probability density function (PDF)" with a description of histogram-based density resampling; `complete` and `mixed` get per-strategy sentences. The previously generated `report.html` was regenerated.
22. **Command docs restructured to the parent-design convention.** `docs/commands/posttree-simulate-alisim.{md,zh.md}` previously used ad-hoc headings (`## alisim params` etc.). They now follow the shared structure required by `2026-06-07-phyloai-design.md` §9 (Purpose / Usage / Inputs / Outputs / Examples / Warnings & Errors / Exit Codes / Notes), with per-subcommand Parameters + Outputs sections, a dedicated Warnings & Errors table, and a Sampling-strategies subsection.
