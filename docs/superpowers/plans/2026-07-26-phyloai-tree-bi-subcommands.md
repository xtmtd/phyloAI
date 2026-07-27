# phyloai tree bi Subcommands Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

> **Commits are user-driven only.** No task in this plan runs `git add` or `git commit`. Each Task ends with a "Report done — no commit" step. The user reviews the diff and stages/commits when they ask.

**Goal:** Convert `phyloai tree bi` from a single flat command into a Click Group with four subcommands — `tree bi pb` (renamed original), `tree bi bpcomp`, `tree bi tracecomp`, and `tree bi readpb` — per the approved design spec `docs/superpowers/specs/2026-07-26-phyloai-tree-bi-subcommands-design.md`.

**Architecture:** Convert `@tree.command("bi", ...)` to `@tree.group("bi", cls=_BiGroup, ...)`, rename `run_bi()` to `run_bi_pb()` in `bi.py`, add three new library modules each containing a single `run_bi_*()` function that returns a payload dict (NOT writes result.json — the CLI layer is the sole writer), and register three new Click subcommands. MCP tool names auto-update via the existing Click-tree walk. All shared parsing helpers remain in `bi.py`. The CLI for the renamed original command (`bi pb`) is a pure rename — no behavior changes.

**Critical conventions applied throughout this plan:**
- Library `run_bi_*()` functions return payload dicts. CLI command functions are the sole writer of `result.json`, matching the existing `bi pb` pattern. No `write_result_json()` calls in library modules.
- `data.tool_stderr` is a **MUST** field per JSON standard §5.3. Per the approved design spec, tool output is streamed directly to terminal (not captured+replayed). Set `data.tool_stderr: ""` for bpcomp (all output streamed to terminal + written to tool output files). For tracecomp: set to captured stdout (stderr → terminal per design). For readpb_mpi: set to concatenated stdout across all modes (stderr → terminal per design).
- `os.path.relpath(chain_dir / name, output_dir)` for all bpcomp/tracecomp chain path references — never hardcoded `../chains/`.
- Dry-run must never call `env.require()` — use placeholder executable names like the existing `bi pb` does.
- Explicit `--chain-names` must validate each name has a corresponding file in `--chain-dir`.
- Error payloads MUST include a full, reproducible `command` field (not just `phyloai tree bi <subcommand>`). Build it from a subcommand-specific flag map, like the existing `_build_command_tokens()` does for `bi pb`.

**Tech Stack:** Python 3.12+, Click, subprocess, os, pathlib, json, numpy, click.rich_click, pytest, Click CliRunner.

**Spec References:**
- Design spec: `docs/superpowers/specs/2026-07-26-phyloai-tree-bi-subcommands-design.md`
- Parent spec: `docs/superpowers/specs/2026-06-23-phyloai-tree-bi-design.md`
- JSON standard: `docs/superpowers/specs/2026-06-21-phyloai-json-output-standard.md`

---

## File Map

| File | Action | Responsibility |
|------|--------|----------------|
| `phyloai/cli/commands/tree.py` | Modify | Convert `tree bi` command → `tree bi` group with 4 subcommands. Per-subcommand error handlers (NOT reusing the old `_build_command_tokens` / `_write_error_result`). |
| `phyloai/tree/bi.py` | Modify | Rename `run_bi()` → `run_bi_pb()`. No other changes. |
| `phyloai/tree/bi_bpcomp.py` | Create | `run_bi_bpcomp()`: chain discovery + existence validation, bpcomp invocation, bpdiff parsing, return payload dict (no file I/O beyond what bpcomp writes). |
| `phyloai/tree/bi_tracecomp.py` | Create | `run_bi_tracecomp()`: chain discovery + existence validation, tracecomp invocation, stdout annotation, return payload dict. |
| `phyloai/tree/bi_readpb.py` | Create | `run_bi_readpb()`: mode validation, readpb_mpi invocation per mode, rr→exchangeabilities post-processing, ss→sitefreq post-processing, return payload dict. |
| `tests/tree/test_bi_bpcomp.py` | Create | Unit + fake-tool tests for bpcomp. |
| `tests/tree/test_bi_tracecomp.py` | Create | Unit + fake-tool tests for tracecomp. |
| `tests/tree/test_bi_readpb.py` | Create | Unit tests for readpb helpers and post-processing. |
| `tests/tree/test_bi.py` | Modify | Update `run_bi` → `run_bi_pb` references. |
| `phyloai/report/collector.py` | Modify | Update `parse_step_id()` and `STEP_ORDER` for `tree.bi.pb`/`tree.bi.bpcomp`/`tree.bi.tracecomp`/`tree.bi.readpb`. |
| `phyloai/report/templates.py` | Modify | Add report templates for the four `tree.bi.*` step_ids; update dispatch map. |
| `tests/report/test_collector.py` | Modify | Add parse_step_id tests for new subcommands. |
| `docs/commands/tree-bi.md` | Modify | Replace with new four-subcommand structure. |
| `docs/commands/tree-bi.zh.md` | Modify | Same as above in Chinese. |
| Various spec docs | Modify | Add superseding notes, update CLI tables (see design spec §11). |
| Various skill files | Modify | Update workflow guidance, parameter annotations, error catalog (see design spec §11). |
| `pyproject.toml` | Modify | Bump version `0.3.0` → `0.4.0`. |
| `phyloai/__init__.py` | Modify | Bump `__version__` `0.3.0` → `0.4.0`. |

---

## Task 1: Refactor `tree bi` into `_BiGroup` and register `tree bi pb`

**Files:**
- Modify: `phyloai/cli/commands/tree.py`
- Modify: `phyloai/tree/bi.py`
- Modify: `tests/tree/test_bi.py`

**Goal:** Convert the existing `@tree.command("bi", ...)` into a Click group with the original command registered as `tree bi pb`. Zero behavior changes.

- [ ] **Step 1: Rename `run_bi` → `run_bi_pb` in `bi.py`**

Rename the `run_bi` function to `run_bi_pb`. Signature unchanged. No other changes.

- [ ] **Step 2: Add `_BiGroup` class**

In `phyloai/cli/commands/tree.py`, add after `_MLGroup`:

```python
class _BiGroup(click.Group):
    def list_commands(self, ctx: click.Context) -> list[str]:
        return ["pb", "bpcomp", "tracecomp", "readpb"]
```
- [ ] **Step 3: Convert `@tree.command("bi", ...)` to `@tree.group("bi", cls=_BiGroup, ...)`**

Change the decorator. Group help text describes all four subcommands (see design spec §6). Remove `cls=_GroupedHelpCommand` from the group decorator.

- [ ] **Step 4: Re-register the original command as `@bi.command("pb", cls=_GroupedHelpCommand, ...)`**

Move all Click decorators from the old `bi_command` to register as `bi pb`. Update help text examples from `phyloai tree bi` to `phyloai tree bi pb`. Change import to `from phyloai.tree.bi import run_bi_pb`.

- [ ] **Step 5: Update `OPTION_GROUPS` key**

Change from `"phyloai tree bi"` to `"phyloai tree bi pb"`. The option list is unchanged.

- [ ] **Step 6: Update `_build_command_tokens()`**

Change `tokens = ["phyloai", "tree", "bi"]` → `tokens = ["phyloai", "tree", "bi", "pb"]`.

- [ ] **Step 7: Update tests**

In `tests/tree/test_bi.py`, update all `run_bi` references to `run_bi_pb`.

- [ ] **Step 8: Run tests**

```bash
pytest tests/tree/test_bi.py -q
pytest tests/cli/ -q -k "bi"
```

---

## Task 2: Create `phyloai/tree/bi_bpcomp.py` — `run_bi_bpcomp()`

**Files:**
- Create: `phyloai/tree/bi_bpcomp.py`
- Create: `tests/tree/test_bi_bpcomp.py`

**Goal:** Implement `run_bi_bpcomp()` per design spec §3. Reuses `_parse_bpcomp_bpdiff()` and `_bpcomp_status()` from `bi.py`. Does NOT reuse `_detect_tools()`. **Returns a payload dict — does NOT write result.json.** The CLI handler writes result.json.

**Design decisions carried into implementation:**
- Dry-run: use placeholder string `"bpcomp"` for the executable; never call `env.require()`.
- Relative chain paths: `os.path.relpath(chain_dir / name, output_dir)`.
- Explicit `--chain-names`: validate every name has a corresponding `*.chain` file in `chain_dir`.
- `data.tool_stderr`: set to `""` since bpcomp stdout/stderr are streamed directly to the terminal (not captured).

- [ ] **Step 1: Create library module**

Implement these internal helpers and the main function:

**`_discover_chain_names(chain_dir: Path) -> list[str]`**: glob `*.chain`, extract stems, sort, raise `FileNotFoundError("No .chain files found in <chain_dir>")` if none.

**`_validate_chain_names(chain_dir: Path, names: list[str]) -> None`**: for each name, assert `<chain_dir>/<name>.chain` exists. Raise `FileNotFoundError` listing missing names if any are absent. Called when `chain_names != "all"`.

**`_build_bpcomp_x_flag(burnin: int, sample_freq: int, until: str) -> list[str]`**: per design spec §3.3 bpcomp `-x` flag rules:
- `--burnin 1000 --sample-freq 1 --until all` → `["-x", "1000"]`
- `--burnin 1000 --sample-freq 10 --until all` → `["-x", "1000", "10"]`
- `--burnin 1000 --sample-freq 1 --until 5000` → `["-x", "1000", "1", "5000"]`
- Only append `<every>` when `sample_freq != 1` OR `until != "all"`.
- Only append `<until>` when `until != "all"`.

**`run_bi_bpcomp(...)`** signature (full parameter list per design spec §3.3):

```python
def run_bi_bpcomp(
    chain_dir: Path = Path("runs/tree/bi/chains"),
    chain_names: str = "all",
    output_dir: Path = Path("runs/tree/bi/bpcomp"),
    overwrite: bool = False,
    burnin: int = 0,
    sample_freq: int = 1,
    until: str = "all",
    cutoff: float = 0.5,
    pb_path: Path | None = None,
    dry_run: bool = False,
    quiet: bool = False,
) -> dict[str, Any]:
```

Key behaviors:
1. Validation: `burnin >= 0`, `sample_freq >= 1`, `0 < cutoff < 1`. Raise `ValueError` on violation.
2. Chain discovery: if `chain_names == "all"`, call `_discover_chain_names()`. Else split + strip, call `_validate_chain_names()`.
3. Output dir: if `overwrite`, `shutil.rmtree` then `mkdir`. Else `mkdir(parents=True, exist_ok=True)`.
4. Tool resolution: if `dry_run`, use `"bpcomp"` as placeholder. Else: `env = ToolEnv(tool_paths={"bpcomp": pb_path / "bpcomp"})` if `pb_path` else `ToolEnv()`, call `env.require("bpcomp")`.
5. Command: `bpcomp -x <burnin> [<every> [<until>]] -c <cutoff> -o bpcomp <rel_paths...>` where each `<rel_path>` = `os.path.relpath(chain_dir / name, output_dir)`.
6. Dry-run: return payload with `key_results` containing placeholder `None` values and `data.cmd`.
7. Execute: `subprocess.run(cmd, cwd=output_dir)` (capture_output=False — terminal passthrough per design spec §3.6). On non-zero exit, return payload with `"status": "error"`.
8. Parse: read `output_dir/bpcomp.bpdiff` with `_parse_bpcomp_bpdiff()`.
9. Summary print (unless quiet): `PhyloAI: maxdiff <v>  meandiff <v>  [<status>]  -> bpcomp/bpcomp.con.tre`.
10. Return payload per design spec §3.8. Set `data.tool_stderr: ""` (all bpcomp output streamed to terminal per design §3.6; results are in the `.bpdiff`/`.bplist`/`.con.tre` output files).

- [ ] **Step 2: Create tests**

- `test_build_x_flag_default/sample_freq/until/full`: verify all four `-x` flag variants.
- `test_discover_chain_names`: success and empty-dir error.
- `test_validate_chain_names`: missing name → `FileNotFoundError`.
- `test_run_bi_bpcomp_dry_run`: verify result shape, chains_used populated, no tools needed.
- `test_run_bi_bpcomp_validation_errors`: negative burnin, zero sample_freq, invalid cutoff → `ValueError`.
- `test_run_bi_bpcomp_with_fake_tool`: create fake bpcomp script, chain .chain files, non-default `chain_dir`/`output_dir`, verify paths use `os.path.relpath`, verify parsed maxdiff/meandiff/status, verify `data.tool_stderr == ""`.

- [ ] **Step 3: Run tests**

```bash
pytest tests/tree/test_bi_bpcomp.py -q
```

---

## Task 3: Create `phyloai/tree/bi_tracecomp.py` — `run_bi_tracecomp()`

**Files:**
- Create: `phyloai/tree/bi_tracecomp.py`
- Create: `tests/tree/test_bi_tracecomp.py`

**Goal:** Implement `run_bi_tracecomp()` per design spec §4. Reuses `_tracecomp_status()` from `bi.py`. Does NOT reuse `_detect_tools()`. **Returns a payload dict — does NOT write result.json.**

- [ ] **Step 1: Create library module**

Implement:

**`_discover_trace_names(chain_dir: Path) -> list[str]`**: glob `*.trace`, extract stems, sort, raise `FileNotFoundError("No .trace files found in <chain_dir>")` if none.

**`_validate_trace_names(chain_dir: Path, names: list[str]) -> None`**: for each name, assert `<chain_dir>/<name>.trace` exists. Raise `FileNotFoundError` listing missing names.

**`_annotate_tracecomp_output(stdout: str) -> tuple[str, float | None, float | None]`**: parse each data line (skip header), compute per-row `_tracecomp_status(effsize, rel_diff)`, append `\t[good]`/`\t[ok]`/`\t[no]`. Return `(annotated_text, min_effsize, max_rel_diff)`.

**`run_bi_tracecomp(...)`** signature:

```python
def run_bi_tracecomp(
    chain_dir: Path = Path("runs/tree/bi/chains"),
    chain_names: str = "all",
    output_dir: Path = Path("runs/tree/bi/tracecomp"),
    overwrite: bool = False,
    burnin: int = 0,
    pb_path: Path | None = None,
    dry_run: bool = False,
    quiet: bool = False,
) -> dict[str, Any]:
```

Key behaviors:
1. Validation: `burnin >= 0`.
2. Chain discovery + validation: same pattern as bpcomp with `.trace` suffix.
3. Output dir: same pattern.
4. Tool resolution: dry-run → placeholder `"tracecomp"`. Else `env.require("tracecomp")`.
5. Command: `tracecomp -x <burnin> <rel_paths...>` where each `<rel_path>` = `os.path.relpath(chain_dir / f"{name}.trace", output_dir)`.
6. Dry-run: return payload with placeholder values.
7. Execute: `subprocess.run(cmd, cwd=output_dir, stdout=PIPE, stderr=None, text=True)`. Stderr → terminal directly per design spec §4.5.
8. Save raw stdout to `output_dir/tracecomp.contdiff`. Set `data.tool_stderr` to captured stdout.
9. Annotate stdout + print annotated table. Print PhyloAI summary line: `min effsize <v>  max rel_diff <v>  [<status>]`.
10. Return payload per design spec §4.6.

- [ ] **Step 2: Create tests**

- `test_discover_trace_names`: success and empty error.
- `test_validate_trace_names`: missing name → `FileNotFoundError`.
- `test_annotate_tracecomp_output`: input with good/ok/no values, verify annotations + min/max.
- `test_run_bi_tracecomp_dry_run`: verify result shape, no tools needed.
- `test_run_bi_tracecomp_validation_errors`: negative burnin.
- `test_run_bi_tracecomp_with_fake_tool`: fake tracecomp that prints to stdout, non-default dirs, verify relpath usage, verify annotation, verify `data.tool_stderr` contains captured stdout.

- [ ] **Step 3: Run tests**

```bash
pytest tests/tree/test_bi_tracecomp.py -q
```

---

## Task 4: Create `phyloai/tree/bi_readpb.py` — `run_bi_readpb()`

**Files:**
- Create: `phyloai/tree/bi_readpb.py`
- Create: `tests/tree/test_bi_readpb.py`

**Goal:** Implement `run_bi_readpb()` per design spec §5. **Returns a payload dict — does NOT write result.json.**

- [ ] **Step 1: Create library module — mode validation and execution**

Function signature:

```python
def run_bi_readpb(
    chain: Path,              # required, path to chain file without extension
    mode: str,                # required, comma-separated modes
    output_dir: Path = Path("runs/tree/bi/readpb"),
    overwrite: bool = False,
    burnin: int = 0,
    sample_freq: int = 1,
    until: str = "all",
    threads: int = 4,
    pb_path: Path | None = None,
    dry_run: bool = False,
    quiet: bool = False,
) -> dict[str, Any]:
```

Key behaviors:
1. Validate `--chain` exists as `<chain>.chain` file (e.g., `Path(str(chain) + ".chain").exists()`).
2. Mode validation per design spec §5.5: unrecognized mode → `ValueError`; duplicate modes → `ValueError`; `allppred` + any of `div`/`sitecomp`/`siteconvprob`/`comp` → `ValueError`.
3. Tool resolution: dry-run → placeholders. Else `env.require("readpb_mpi")` + `env.require("mpirun")`.
4. Per-mode execution: each mode = separate `readpb_mpi` invocation. Working directory = `chain.parent`. Mode flag mapping per design spec §5.4 table. Command: `mpirun -np <threads> readpb_mpi -x <burnin> [<every> [<until>]] <mode_flag> <chain_stem>`. Immediately after each successful invocation, detect that mode's generated files and move them to `output_dir` before starting the next mode. Move `ppred` mode files to `output_dir/ppred/`; move `allppred`'s `<chain>.ppred` and every other mode's files directly to `output_dir`.
5. `-x` flag: reuse `_build_bpcomp_x_flag()` logic from `bi_bpcomp.py` or make a shared helper.
6. Output capture per design spec §5.6: `stdout=PIPE, stderr=None` (stderr/progress → terminal directly). Save stdout per mode to `<output-dir>/<mode>.stdout`. Replay stdout to terminal unchanged unless `--quiet`. After all modes complete, build a single `data.tool_stderr` string by concatenating all modes' stdout, separated by `\n--- <mode> ---\n` headers. The `params` dict records the full mode string; `data.tool_stderr` is a single flat string (not per-mode keys).
7. Post-processing triggers: after moving the raw `rr` output to `output_dir`, call `_convert_meanrr_to_exchangeabilities()` there. After moving the raw `ss` output to `output_dir`, call `_convert_siteprofiles_to_sitefreq()` there. Both derived files therefore appear before the next mode begins.
8. Return payload per design spec §5.9.

- [ ] **Step 2: Implement `_convert_meanrr_to_exchangeabilities(meanrr_path: Path) -> Path`**

Algorithm (from design spec §5.7), using numpy, no pandas:
1. Read first non-empty line → `order_of_aa` (list of 20 AA symbols, space-separated).
2. Build `aa_to_idx = {aa: i for i, aa in enumerate(order_of_aa)}`.
3. Initialize `exch = np.zeros((20, 20), dtype=np.float64)`.
4. Parse each remaining non-empty line: `source target value` → set symmetric entries.
5. Reindex to PAML order: `PAML_ORDER = ['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']`. `paml_idx = [aa_to_idx[aa] for aa in PAML_ORDER]`. `paml_exch = exch[np.ix_(paml_idx, paml_idx)]`.
6. Write lower triangle for rows `i=0..19`, columns `0..i-1`: one row per line, values `%08.6f` space-separated, trailing space, then newline. **Row `i=0` has no lower-triangle values → the required leading blank line.**
7. Append blank line, then `"0.050000 "` × 20 + newline (uniform prior state frequencies placeholder).
8. Output file: `<chain>.exchangeabilities` (same directory as chain). Return the Path.

- [ ] **Step 3: Implement `_convert_siteprofiles_to_sitefreq(siteprofiles_path: Path) -> Path`**

Algorithm (from design spec §5.7):
1. Skip first 2 header lines.
2. PhyloBayes AA order: `A C D E F G H I K L M N P Q R S T V W Y`.
3. IQ-TREE AA order: `A R N D C Q E G H I L K M F P S T W Y V`.
4. Precompute index map: for each position in IQ-TREE order, find index in PhyloBayes order.
5. For each data line: parse site index + 20 float frequencies. Reindex to IQ-TREE order. Floor zeros/near-zeros to `1e-8`. Re-normalize to sum to 1.
6. Write `<chain>.sitefreq`: each line = `<site_index> <20 %.8f floats>`. Return Path.

- [ ] **Step 4: Create tests**

- `test_validate_modes_unrecognized/duplicates/allppred_conflict`: → `ValueError`.
- `test_convert_meanrr_to_exchangeabilities`: create sample `.meanrr`, verify output shape (21 rows, row 0 blank), PAML order correct, values symmetric, trailing spaces present.
- `test_convert_siteprofiles_to_sitefreq`: create sample `.siteprofiles`, verify reindex order, `1e-8` floor, row sum ≈ 1.0, format `%.8f`.
- `test_run_bi_readpb_dry_run`: verify result shape, modes_run, `data.cmds` dict per mode, no tools needed.
- `test_run_bi_readpb_ss_rr_roundtrip`: fake `readpb_mpi` + `mpirun` → end-to-end with actual `.meanrr` + `.siteprofiles` output files → verify raw and post-processing output files exist under `output_dir`, are absent from the chain directory, and are valid.
- `test_run_bi_readpb_moves_each_mode_before_next_starts`: fake `readpb_mpi` writes rr output, then its ss invocation asserts the rr raw and converted files already exist under `output_dir`; verify `ppred` files are in `output_dir/ppred/` and all other mode outputs are directly in `output_dir`.

- [ ] **Step 5: Run tests**

```bash
pytest tests/tree/test_bi_readpb.py -q
```

---

## Task 5: Register CLI commands for `tree bi bpcomp`, `tree bi tracecomp`, `tree bi readpb`

**Files:**
- Modify: `phyloai/cli/commands/tree.py`

**Goal:** Add three new Click command registrations. Each follows the `bi pb` pattern: validate kwargs → call library function → write result.json (with backup) → handle errors. **Do NOT reuse the old `_build_command_tokens()` or `_write_error_result()` — those are specific to `bi pb` and emit wrong command paths and have incomplete flag maps for the new subcommands.**

- [ ] **Step 1: Add per-subcommand error helper**

Add a helper that builds a full CLI command string from kwargs. Each subcommand has its own flag map:

```python
# Flag maps for bpcomp / tracecomp / readpb
_BPCOMP_FLAG_MAP = {
    "chain_dir": "--chain-dir", "chain_names": "--chain-names",
    "output_dir": "--output-dir", "burnin": "--burnin",
    "sample_freq": "--sample-freq", "until": "--until",
    "cutoff": "--cutoff", "pb_path": "--pb-path",
}
_TRACECOMP_FLAG_MAP = {
    "chain_dir": "--chain-dir", "chain_names": "--chain-names",
    "output_dir": "--output-dir", "burnin": "--burnin",
    "pb_path": "--pb-path",
}
_READPB_FLAG_MAP = {
    "chain": "--chain", "mode": "--mode",
    "output_dir": "--output-dir", "burnin": "--burnin",
    "sample_freq": "--sample-freq", "until": "--until",
    "threads": "--threads", "pb_path": "--pb-path",
}

def _build_bi_subcommand_tokens(subcommand: str, kwargs: dict, flag_map: dict[str, str]) -> list[str]:
    """Build full CLI token list for a tree bi subcommand error payload."""
    tokens = ["phyloai", "tree", "bi", subcommand]
    for key, flag in flag_map.items():
        val = kwargs.get(key)
        if val is None:
            continue
        if isinstance(val, bool) and val:
            tokens.append(flag)
        elif isinstance(val, bool):
            continue
        else:
            tokens.append(flag)
            tokens.append(str(val))
    if kwargs.get("overwrite"):
        tokens.append("--overwrite")
    if kwargs.get("dry_run"):
        tokens.append("--dry-run")
    if kwargs.get("quiet"):
        tokens.append("--quiet")
    return tokens

def _write_bi_error_result(kwargs: dict, message: str, exit_code: int, subcommand: str, flag_map: dict[str, str]) -> None:
    """Write error result.json for tree bi subcommands."""
    output_dir = Path(kwargs.get("output_dir", f"runs/tree/bi/{subcommand}"))
    output_dir.mkdir(parents=True, exist_ok=True)
    tokens = _build_bi_subcommand_tokens(subcommand, kwargs, flag_map)
    result = {
        "status": "error",
        "command": " ".join(tokens),
        "wall_time": 0.0,
        "tool_versions": {},
        "params": {k: (str(v) if isinstance(v, Path) else v) for k, v in kwargs.items()},
        "key_results": {},
        "error": message,
        "data": {"cmd": [], "tool_stderr": ""},
    }
    (output_dir / "result.json").write_text(json.dumps(result, indent=2))
```

The existing `_build_command_tokens()` and `_write_error_result()` remain unchanged for `bi pb`. Each subcommand's CLI command function calls the new helper with the appropriate flag map.

- [ ] **Step 2: Register `@bi.command("bpcomp", cls=_GroupedHelpCommand, ...)`**

Parameters per design spec §3.3. Help text references `phyloai tree bi bpcomp`. Command body:

```python
def bi_bpcomp_command(**kwargs):
    from phyloai.tree.bi_bpcomp import run_bi_bpcomp
    try:
        payload = run_bi_bpcomp(**kwargs)
    except FileNotFoundError as exc:
        _write_bi_error_result(kwargs, str(exc), 3, "bpcomp", _BPCOMP_FLAG_MAP)
        _fail(str(exc), 3)
    except ValueError as exc:
        _write_bi_error_result(kwargs, str(exc), 1, "bpcomp", _BPCOMP_FLAG_MAP)
        _fail(str(exc), 1)
    if kwargs.get("dry_run"):
        if not kwargs.get("quiet"):
            click.echo(" ".join(payload["data"]["cmd"]))
        return
    result_path = kwargs["output_dir"] / "result.json"
    # backup logic (same as bi pb)
    if result_path.exists():
        from datetime import datetime
        ts = datetime.now().strftime("%Y%m%d_%H%M%S")
        shutil.copy2(str(result_path), str(result_path.with_name(f"result_{ts}.json")))
    with open(result_path, "w") as fh:
        json.dump(payload, fh, indent=2)
    if payload["status"] == "error":
        _fail(payload.get("error", "bpcomp failed"), 2)
    if not kwargs.get("quiet"):
        click.echo(f"Results saved to {result_path}", err=True)
```

- [ ] **Step 3: Register `@bi.command("tracecomp", cls=_GroupedHelpCommand, ...)`**

Parameters per design spec §4.3. Same command body pattern as bpcomp, calling `run_bi_tracecomp()` and using `_TRACECOMP_FLAG_MAP` in error handler calls.

- [ ] **Step 4: Register `@bi.command("readpb", cls=_GroupedHelpCommand, ...)`**

Parameters per design spec §5.3. Same pattern, calling `run_bi_readpb()` and using `_READPB_FLAG_MAP` in error handler calls.

- [ ] **Step 5: Verify CLI structure**

```bash
phyloai tree bi --help
phyloai tree bi pb --help
phyloai tree bi bpcomp --help
phyloai tree bi tracecomp --help
phyloai tree bi readpb --help
```

- [ ] **Step 6: Run CLI tests**

Add CLI tests for help output, invalid parameters, and dry-run output for each new subcommand. Then:

```bash
pytest tests/cli/test_tree_bi.py -q
```

---

## Task 6: Update report collector and report templates

**Files:**
- Modify: `phyloai/report/collector.py`
- Modify: `phyloai/report/templates.py`
- Modify: `tests/report/test_collector.py`

**Goal:** Update `parse_step_id()` to correctly parse `tree bi pb`, `tree bi bpcomp`, `tree bi tracecomp`, `tree bi readpb` from command strings (e.g., `"phyloai tree bi bpcomp ..."`). Update `STEP_ORDER`. Add report template functions for each step_id. Add tests for new step_id parsing.

**IMPORTANT:** The plan in the original design spec §12 incorrectly references template keys like `bi_bpcomp`. The actual step_id format uses dots: `tree.bi.pb`, `tree.bi.bpcomp`, `tree.bi.tracecomp`, `tree.bi.readpb`. All registration must use these dotted names.

- [ ] **Step 1: Update `parse_step_id()`**

In `_THIRD_LEVEL`, add entry: `"bi": {"pb", "bpcomp", "tracecomp", "readpb"}`. This allows `phyloai tree bi bpcomp --chain-dir ...` to parse as `tree.bi.bpcomp`.

- [ ] **Step 2: Update `STEP_ORDER`**

Replace the single `"tree.bi"` entry with:
```python
"tree.bi.pb",
"tree.bi.bpcomp",
"tree.bi.tracecomp",
"tree.bi.readpb",
```

- [ ] **Step 3: Add report template functions**

In `phyloai/report/templates.py`, add four template functions (one per subcommand) per design spec §12. Register them in the methods dispatch map using the dotted step_id keys. Templates:
- `generate_methods_tree_bi_pb()` — wraps existing `generate_methods_tree_bi()`, or rename the existing one
- `generate_methods_tree_bi_bpcomp()` — per spec §12 bpcomp paragraph + table
- `generate_methods_tree_bi_tracecomp()` — per spec §12 tracecomp paragraph + table
- `generate_methods_tree_bi_readpb()` — per spec §12 readpb paragraphs + table

- [ ] **Step 4: Add tests for new parse_step_id outputs**

```python
def test_parse_step_id_tree_bi_bpcomp():
    assert parse_step_id("phyloai tree bi bpcomp --chain-dir chains --burnin 1000") == "tree.bi.bpcomp"

def test_parse_step_id_tree_bi_tracecomp():
    assert parse_step_id("phyloai tree bi tracecomp --chain-dir chains --burnin 5000") == "tree.bi.tracecomp"

def test_parse_step_id_tree_bi_readpb():
    assert parse_step_id("phyloai tree bi readpb --chain chain1 --mode ss,rr") == "tree.bi.readpb"

def test_parse_step_id_tree_bi_pb():
    assert parse_step_id("phyloai tree bi pb --matrix m.phy --chains 2") == "tree.bi.pb"
```

- [ ] **Step 5: Run tests**

```bash
pytest tests/report/test_collector.py -q
```

---

## Task 7: Verify MCP schema auto-generation

**Files:**
- Inspect: `phyloai/mcp/schema_gen.py` (no changes expected).

**Goal:** Confirm `walk_click_tree()` discovers `phyloai_tree_bi_pb`, `phyloai_tree_bi_bpcomp`, `phyloai_tree_bi_tracecomp`, `phyloai_tree_bi_readpb`. The old `phyloai_tree_bi` should not appear.

- [ ] **Step 1: Add MCP schema test**

```python
def test_tree_bi_subcommands_in_mcp_tools():
    from phyloai.mcp.tools.cli_tools import get_tool_definitions
    tools = get_tool_definitions()
    for name in ["phyloai_tree_bi_pb", "phyloai_tree_bi_bpcomp",
                 "phyloai_tree_bi_tracecomp", "phyloai_tree_bi_readpb"]:
        assert name in tools, f"{name} missing from MCP tools"
    assert "phyloai_tree_bi" not in tools, "old phyloai_tree_bi must not exist"
```

- [ ] **Step 2: Run test**

```bash
pytest tests/mcp/ -q -k "tree_bi_subcommands"
```

---

## Task 8: Documentation and version bump

**Files:**
- Various docs and skill files per design spec §11.

- [ ] **Step 1: Bump version to `0.4.0`**

`pyproject.toml`: `version = "0.4.0"`.
`phyloai/__init__.py`: `__version__ = "0.4.0"`.

- [ ] **Step 2: Update command docs**

- `docs/commands/tree-bi.md`: Replace with four subcommand sections. `tree bi pb` = original content. Add `tree bi bpcomp`/`tree bi tracecomp`/`tree bi readpb` with parameter tables, examples, output file descriptions.
- `docs/commands/tree-bi.zh.md`: Same in Chinese.

- [ ] **Step 3: Update design specs**

Per design spec §11 — add superseding notes, update CLI tables.

- [ ] **Step 4: Update skill files**

Per design spec §11 — update workflow guidance, parameter annotations, error catalog.

---

## Task 9: Final integration verification

- [ ] **Step 1: Run full test suite**

```bash
pytest tests/tree/ -q
pytest tests/cli/ -q
pytest tests/report/ -q
pytest tests/mcp/ -q
```

- [ ] **Step 2: Verify CLI help for all subcommands**

```bash
phyloai tree bi --help
phyloai tree bi pb --help
phyloai tree bi bpcomp --help
phyloai tree bi tracecomp --help
phyloai tree bi readpb --help
```

- [ ] **Step 3: Verify MCP tool listing**

```bash
python -c "from phyloai.mcp.tools.cli_tools import get_tool_definitions; tools = get_tool_definitions(); print([k for k in tools if 'tree_bi' in k])"
```

Expected output: `['phyloai_tree_bi_pb', 'phyloai_tree_bi_bpcomp', 'phyloai_tree_bi_tracecomp', 'phyloai_tree_bi_readpb']`.

---

## Implementation Order

```
Task 1 (refactor bi → group + pb)  ── required foundation
    │
    ├── Task 2 (bi_bpcomp.py)   ──┐
    ├── Task 3 (bi_tracecomp.py)──┤  can run in parallel
    └── Task 4 (bi_readpb.py)   ──┘
              │
    Task 5 (CLI registration)     ── depends on Tasks 2-4
    Task 6 (report collector)     ── depends on Task 1
    Task 7 (MCP verify)           ── depends on Task 5
              │
    Task 8 (docs + version bump)  ── independent, run anytime
    Task 9 (final verification)   ── runs last
```
