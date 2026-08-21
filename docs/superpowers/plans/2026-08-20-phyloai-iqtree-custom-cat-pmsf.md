# Custom IQ-TREE CAT-PMSF Inputs Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add structured custom IQ-TREE exchangeability and per-site frequency-profile inputs to `phyloai tree ml iqtree`, enabling CAT-PMSF-style ML analyses without copying input files.

**Architecture:** Treat an existing regular file supplied through `--model` as an AA custom exchangeability matrix, resolve it to an absolute path, and let IQ-TREE validate its content. Add `--site-freq-file` as the structured `-fs` counterpart, while preserving existing `--tool-args` override semantics: raw `-fs` suppresses the generated flag. Persist the resolved structured inputs in `result.json` and checkpoint data so reporting and resume remain reproducible.

**Tech Stack:** Python 3.10+, Click, Biopython, IQ-TREE3, pytest.

**Spec:** `docs/superpowers/specs/2026-06-19-phyloai-tree-ml-iqtree-design.md`

## Global Constraints

- Release version: `0.5.0` in both `pyproject.toml` and `phyloai/__init__.py`.
- Do not add dependencies, custom model parsers, new top-level CLI commands, or a separate MCP wrapper.
- Do not copy, move, mutate, or content-validate the custom exchangeability matrix or site-frequency profile; resolve structured input paths to absolute paths and let IQ-TREE validate file content.
- Structured custom-model and `--site-freq-file` workflows are AA, `--matrix`, and `--modelfinder none` only.
- `--site-freq-file` requires `--state-freq none`; raw `--tool-args "-fs PATH"` also requires `--state-freq none` and overrides structured `--site-freq-file`.
- Preserve `--pmsf-base-model` for built-in C10–C60 PMSF only; do not expand or reinterpret it.
- Preserve existing `--tool-args` strategy override semantics. `-s` and shell redirects remain blocked; `-m` and `-fs` are overrideable.
- Keep all code, comments, user-visible help, documentation, tests, and report prose in English, except the established Chinese documentation and Skill annotation files.

---

## File Structure

| File | Responsibility |
|---|---|
| `phyloai/tree/ml_iqtree.py` | Resolve and validate custom files, construct `-m`/`-fs`, preserve inputs in checkpoints and `result.json`. |
| `phyloai/cli/commands/tree.py` | Expose `--site-freq-file`, place it in the Heterogeneous help group, and pass it to the core runner. |
| `tests/tree/test_ml_iqtree.py` | Unit tests for model detection, validation, IQ-TREE command construction, dry-run payload, and result metadata. |
| `tests/cli/test_tree.py` | CLI help and invalid-combination coverage. |
| `tests/mcp/test_schema_gen.py` | Assert Click-derived MCP schema exposes `site_freq_file` as a path. |
| `tests/report/test_templates.py` | Verify custom-model methods prose. |
| `phyloai/report/templates.py` | Describe custom exchangeabilities and structured site-specific profiles accurately. |
| `docs/commands/tree-ml-iqtree.md` | English command reference and CAT-PMSF example. |
| `docs/commands/tree-ml-iqtree.zh.md` | Chinese command reference and CAT-PMSF example. |
| `README.md`, `README.zh.md` | Concise capability-table wording. |
| `skills/phyloai-workflow/SKILL.md` | Mention the supported structured custom CAT-PMSF ML route. |
| `skills/phyloai-workflow/references/parameter-annotations.md` | Chinese parameter guidance for custom `--model`, `--site-freq-file`, and raw `-fs` override. |
| `pyproject.toml`, `phyloai/__init__.py` | Release version `0.5.0`. |

## Task 1: Add core custom-model and profile semantics

**Files:**
- Modify: `phyloai/tree/ml_iqtree.py`
- Test: `tests/tree/test_ml_iqtree.py`

**Interfaces:**
- Consumes: `run_iqtree(..., model: str | None, state_freq: str, rate_heterogeneity: str, modelfinder: str, tool_args: str | None, ...)`.
- Produces: `run_iqtree(..., site_freq_file: str | Path | None = None, ...)` and `_build_iqtree_cmd(..., site_freq_file: str | None = None, ...)`.
- Produces: a resolved custom model string such as `/abs/chain1.exchangeabilities+R4`; generated `-fs /abs/chain1.sitefreq` unless raw `tool_args` contains `-fs`.

- [ ] **Step 1: Write failing core tests for custom models and `-fs` command construction**

Add to `tests/tree/test_ml_iqtree.py` beside the existing `_build_model_string`, `_run_validations`, and `_build_iqtree_cmd` tests:

```python
def test_build_model_string_custom_exchangeabilities(tmp_path: Path) -> None:
    from phyloai.tree.ml_iqtree import _build_model_string

    model = tmp_path / "chain1.exchangeabilities"
    model.write_text("0.5\n")

    assert _build_model_string(
        model=str(model.resolve()), state_freq="none",
        rate_heterogeneity="+R4", modelfinder="none",
    ) == f"{model.resolve()}+R4"


def test_build_iqtree_cmd_adds_site_freq_file(tmp_path: Path) -> None:
    from phyloai.tree.ml_iqtree import _build_iqtree_cmd

    profile = tmp_path / "chain1.sitefreq"
    cmd = _build_iqtree_cmd(
        input_path=tmp_path / "matrix.fa", prefix=tmp_path / "matrix",
        model_string="/models/chain1.exchangeabilities+R4", seq_type="AA",
        modelfinder="none", boot=0, alrt=None, bnni=False,
        mode="normal", threads_arg="-T 1", site_freq_file=str(profile),
    )

    assert cmd.count("-fs") == 1
    assert cmd[cmd.index("-fs") + 1] == str(profile)


def test_build_iqtree_cmd_tool_args_fs_overrides_structured_profile(tmp_path: Path) -> None:
    from phyloai.tree.ml_iqtree import _build_iqtree_cmd

    cmd = _build_iqtree_cmd(
        input_path=tmp_path / "matrix.fa", prefix=tmp_path / "matrix",
        model_string="/models/chain1.exchangeabilities+R4", seq_type="AA",
        modelfinder="none", boot=0, alrt=None, bnni=False,
        mode="normal", threads_arg="-T 1", site_freq_file="/managed.sitefreq",
        tool_args="-fs /override.sitefreq",
    )

    assert cmd.count("-fs") == 1
    assert cmd[cmd.index("-fs") + 1] == "/override.sitefreq"
```

Add validation tests covering each rejected structured case:

```python
@pytest.mark.parametrize(
    ("batch_mode", "seq_type", "modelfinder", "state_freq", "match"),
    [
        (True, "AA", "none", "none", "--matrix"),
        (False, "NT", "none", "none", "AA"),
        (False, "AA", "MF", "none", "ModelFinder"),
        (False, "AA", "none", "none", "custom model"),
        (False, "AA", "none", "+F", "--state-freq none"),
    ],
)
def test_validate_site_freq_file_rejects_invalid_context(...):
    ...
```

Also test that `--tool-args "-fs /override.sitefreq"` with `state_freq="+F"` raises an error mentioning `--state-freq none`, and that a custom model file is rejected in batch mode, for NT data, and with ModelFinder.

- [ ] **Step 2: Run the focused tests to verify they fail**

Run:

```bash
pytest tests/tree/test_ml_iqtree.py -k 'custom_exchangeabilities or site_freq_file or tool_args_fs' -v
```

Expected: FAIL because the new keyword argument and validation behavior do not exist.

- [ ] **Step 3: Implement minimal custom-file detection and validation**

In `phyloai/tree/ml_iqtree.py`:

1. Add a small helper that recognizes an existing custom model path, rejects an existing non-file path, and returns `str(path.resolve())` for a regular file. Standard model names continue through `_validate_model()` unchanged.
2. Resolve `site_freq_file` with the same regular-file requirement. Do not inspect file contents.
3. Thread `custom_model: bool`, `site_freq_file`, and `tool_args` into `_run_validations()`.
4. Enforce the Global Constraints exactly. For a raw `-fs`, detect the exact token with `shlex.split(tool_args)`; do not substring-match arbitrary values.
5. Keep `_build_model_string()` simple: custom paths follow the existing `base + optional state-frequency + optional rate` logic, and structured profile validation guarantees no state-frequency suffix is appended.

The helper shape should remain small and local:

```python
def _resolve_custom_model_path(model: str) -> str | None:
    candidate = Path(model).expanduser()
    if not candidate.exists():
        return None
    if not candidate.is_file():
        raise ValueError(f"Custom model path is not a regular file: {candidate}")
    return str(candidate.resolve())
```

Use the returned value as `model` for validation, model-string construction, reproducibility metadata, and command assembly. Retain the standard-model validator for non-file values.

- [ ] **Step 4: Implement managed `-fs` emission with raw override**

Add `site_freq_file` to `_build_iqtree_cmd()`. Use the existing `_is_flag_overridden()` pattern:

```python
if site_freq_file and not _is_flag_overridden("-fs", tool_tokens):
    cmd.extend(["-fs", site_freq_file])
```

Do not add `-fs` to `_IQTREE_BLOCKED_FLAGS`. It must be overrideable and raw tool args remain appended last. Ensure every caller of `_build_iqtree_cmd()` passes the new optional argument.

- [ ] **Step 5: Run focused tests to verify they pass**

Run:

```bash
pytest tests/tree/test_ml_iqtree.py -k 'custom_exchangeabilities or site_freq_file or tool_args_fs' -v
```

Expected: PASS.

- [ ] **Step 6: Run the entire IQ-TREE core test module**

Run:

```bash
pytest tests/tree/test_ml_iqtree.py -v
```

Expected: PASS with all existing IQ-TREE tests retained.

- [ ] **Step 7: Commit the core behavior**

```bash
git add phyloai/tree/ml_iqtree.py tests/tree/test_ml_iqtree.py
git commit -m "feat(iqtree): support custom exchangeability profiles"
```

## Task 2: Wire the CLI, reproducibility metadata, checkpoint data, and generated command

**Files:**
- Modify: `phyloai/cli/commands/tree.py`
- Modify: `phyloai/tree/ml_iqtree.py`
- Test: `tests/cli/test_tree.py`
- Test: `tests/tree/test_ml_iqtree.py`

**Interfaces:**
- Consumes: core `run_iqtree(..., site_freq_file=...)` from Task 1.
- Produces: Click option `--site-freq-file` as `Path | None`; `result.json.params.site_freq_file: str | None`; checkpoint/reproducible command includes the resolved structured profile.

- [ ] **Step 1: Write failing CLI and payload tests**

Add a CLI help assertion that `--site-freq-file` appears after `--guide-tree` in the `Heterogeneous:` group. Add a dry-run test:

```python
def test_tree_ml_iqtree_custom_profile_dry_run(tmp_path: Path) -> None:
    matrix = tmp_path / "matrix.fa"
    matrix.write_text(">a\nMKTLLL\n>b\nMKTLLL\n")
    model = tmp_path / "chain1.exchangeabilities"
    model.write_text("0.5\n")
    profile = tmp_path / "chain1.sitefreq"
    profile.write_text("1 " + " ".join(["0.05"] * 20) + "\n")

    result = CliRunner().invoke(cli, [
        "tree", "ml", "iqtree", "--matrix", str(matrix),
        "--seq-type", "AA", "--model", str(model),
        "--site-freq-file", str(profile), "--state-freq", "none",
        "--rate-heterogeneity", "+R4", "--boot", "0", "--threads", "1",
        "--dry-run", "--quiet",
    ])

    assert result.exit_code == 0
```

Add direct-API dry-run assertions that `payload["data"]["cmd"]` contains absolute `-m`, absolute `-fs`, and `payload["params"]["site_freq_file"] == str(profile.resolve())`. Add a second test where raw `tool_args="-fs /tmp/override.sitefreq"` leaves `params.site_freq_file` as `None` and emits only the raw profile.

- [ ] **Step 2: Run the new tests to verify they fail**

Run:

```bash
pytest tests/cli/test_tree.py -k 'site_freq or custom_profile' -v
pytest tests/tree/test_ml_iqtree.py -k 'site_freq or custom_model' -v
```

Expected: FAIL because Click has no option and result/checkpoint serialization lacks the field.

- [ ] **Step 3: Add Click option and Heterogeneous help text**

In `phyloai/cli/commands/tree.py`:

1. Add `@click.option("--site-freq-file", type=click.Path(dir_okay=False, path_type=Path), default=None, ...)` adjacent to `--guide-tree`.
2. Update the `--model` Click help to state that, besides supported built-in model names, it may be an existing IQ-TREE custom model file path.
3. Add `site_freq_file: Path | None` to `iqtree_command()` and pass its string value to `run_iqtree()`.
4. In `_IQTreeCommand._HELP_GROUPS`, put this row between `--guide-tree` and `--qmax`:

```text
--site-freq-file   Per-site AA frequency profile for a custom --model; maps to IQ-TREE -fs
```

4. Add one concise CAT-PMSF-style example to the Heterogeneous explanation and Workflow Examples. State that it needs `--state-freq none` and that `--tool-args "-fs PATH"` overrides the structured option.

Do not duplicate core validation in the CLI beyond the existing basic path-existence checks; direct library callers must receive the same protections from `run_iqtree()`.

- [ ] **Step 4: Persist effective values consistently**

In `phyloai/tree/ml_iqtree.py`, add `site_freq_file` to all effective parameter flows:

- `run_iqtree()` signature and both single/batch paths;
- `_resolved_iqtree_params()` so batch checkpoint resume detects a changed structured profile;
- `_assemble_iqtree_result()` command reconstruction and `params` object;
- initial batch checkpoint command and params;
- single-mode `result.json` payload.

Use the resolved absolute structured path. When raw `--tool-args -fs` overrides the structured value, pass `None` as the effective `site_freq_file` so `result.json` does not falsely describe the raw path as structured. Leave raw `tool_args` untouched.

- [ ] **Step 5: Run focused tests to verify they pass**

Run:

```bash
pytest tests/cli/test_tree.py -k 'site_freq or custom_profile or iqtree_help' -v
pytest tests/tree/test_ml_iqtree.py -k 'site_freq or custom_model or resume' -v
```

Expected: PASS.

- [ ] **Step 6: Commit CLI and metadata wiring**

```bash
git add phyloai/cli/commands/tree.py phyloai/tree/ml_iqtree.py \
  tests/cli/test_tree.py tests/tree/test_ml_iqtree.py
git commit -m "feat(iqtree): expose site frequency profile option"
```

## Task 3: Make report prose and MCP schema truthful

**Files:**
- Modify: `phyloai/report/templates.py`
- Test: `tests/report/test_templates.py`
- Test: `tests/mcp/test_schema_gen.py`

**Interfaces:**
- Consumes: `result.json.params.model`, `state_freq`, `rate_heterogeneity`, and `site_freq_file` from Task 2.
- Produces: methods text that distinguishes built-in models from custom exchangeabilities; automatically generated MCP `tree_ml_iqtree` schema with `site_freq_file` path property.

- [ ] **Step 1: Write failing methods and schema tests**

Add to `TestIqtree` in `tests/report/test_templates.py`:

```python
def test_custom_exchangeabilities_and_profile(self):
    text = generate_all_methods(
        "tree.ml.iqtree",
        params={
            "modelfinder": "none",
            "model": "/abs/chain1.exchangeabilities",
            "state_freq": "none",
            "rate_heterogeneity": "+R4",
            "site_freq_file": "/abs/chain1.sitefreq",
            "boot": None,
        },
        key_results={}, tool_versions={"iqtree3": "3.1.2"},
    )

    assert "custom exchangeability matrix" in text
    assert "site-specific state-frequency profiles" in text
    assert "+R4" in text
```

Add to `tests/mcp/test_schema_gen.py`:

```python
def test_tree_ml_iqtree_schema_exposes_site_freq_file() -> None:
    descriptor = next(d for d in walk_click_tree(cli) if d["tool_name"] == "tree_ml_iqtree")
    props = build_mcp_tool(descriptor)["inputSchema"]["properties"]

    assert props["site_freq_file"]["type"] == "string"
    assert props["site_freq_file"]["format"] == "path"
```

- [ ] **Step 2: Run the tests to verify they fail**

Run:

```bash
pytest tests/report/test_templates.py -k custom_exchangeabilities -v
pytest tests/mcp/test_schema_gen.py -k site_freq_file -v
```

Expected: report test FAILS because prose reports an uppercased path as a normal model; schema test FAILS until Task 2's Click option is present.

- [ ] **Step 3: Implement the narrow report branch**

In `generate_methods_tree_ml_iqtree()`:

1. Treat an absolute `params["model"]` path as the structured custom-model signal only when `modelfinder` is `none`. Do not call `Path.is_file()` here: report generation must remain accurate even if the original input file has since been moved or deleted.
2. For that branch, say a user-specified custom exchangeability matrix was used, append the non-`none` rate heterogeneity, and do not uppercase or fabricate a `+F` suffix.
3. If `site_freq_file` is non-empty, add that IQ-TREE site-specific state-frequency profiles (`-fs`) were used.
4. Keep current ModelFinder and built-in model wording unchanged.

The report must not label the run “standard PMSF”; use only neutral custom-input wording.

- [ ] **Step 4: Run the focused tests to verify they pass**

Run:

```bash
pytest tests/report/test_templates.py -v
pytest tests/mcp/test_schema_gen.py -v
```

Expected: PASS.

- [ ] **Step 5: Commit report/schema coverage**

```bash
git add phyloai/report/templates.py tests/report/test_templates.py tests/mcp/test_schema_gen.py
git commit -m "fix(report): describe custom IQ-TREE model inputs"
```

## Task 4: Update user documentation, workflow guidance, and release version

**Files:**
- Modify: `docs/commands/tree-ml-iqtree.md`
- Modify: `docs/commands/tree-ml-iqtree.zh.md`
- Modify: `README.md`
- Modify: `README.zh.md`
- Modify: `skills/phyloai-workflow/SKILL.md`
- Modify: `skills/phyloai-workflow/references/parameter-annotations.md`
- Modify: `pyproject.toml`
- Modify: `phyloai/__init__.py`
- Test: `tests/cli/test_main.py`

**Interfaces:**
- Consumes: final CLI contract from Tasks 1–2.
- Produces: accurate bilingual command docs, parameter-card context, and package/CLI version `0.5.0`.

- [ ] **Step 1: Write the version test**

Add a focused assertion to `tests/cli/test_main.py` using the existing CLI test conventions:

```python
def test_version_is_0_5_0() -> None:
    result = CliRunner().invoke(cli, ["--version"])

    assert result.exit_code == 0
    assert "0.5.0" in result.output
```

- [ ] **Step 2: Run the version test to verify it fails**

Run:

```bash
pytest tests/cli/test_main.py -k version -v
```

Expected: FAIL because the package reports `0.4.1`.

- [ ] **Step 3: Update English and Chinese command documentation**

In both `docs/commands/tree-ml-iqtree.md` and `.zh.md`:

1. Add a CAT-PMSF-style example using exactly the three `runs/test/` inputs, `--state-freq none`, `--rate-heterogeneity +R4`, `--boot 0`, and `--threads 1`.
2. Add `--site-freq-file` to the Heterogeneous Models table and clarify that it maps to IQ-TREE `-fs` for a user-supplied AA custom model.
3. Clarify that an existing regular file can be supplied to `--model`; PhyloAI resolves it to an absolute path and does not copy it.
4. Document all restrictions and raw override semantics: structured input is single-matrix AA/non-ModelFinder; `--state-freq none` is required; raw `--tool-args "-fs /absolute/profile"` overrides structured `--site-freq-file`, remains raw, and should use an absolute path.
5. Do not describe this generic facility as standard PMSF.

- [ ] **Step 4: Update README and workflow Skill text minimally**

- In both READMEs, expand only the `tree ml iqtree` table description to mention custom exchangeability matrices and site-specific profiles.
- In `skills/phyloai-workflow/SKILL.md`, add one short note in the Tree workflow about the structured CAT-PMSF-style ML route, separate from `tree bi readpb` simulation partitions.
- In `skills/phyloai-workflow/references/parameter-annotations.md`, update `--model`; add `--site-freq-file`; correct `--state-freq` descriptions to reflect the actual choices; explain the `--tool-args -fs` override and required `none` state frequency.

- [ ] **Step 5: Bump the release version**

Set both declarations exactly:

```toml
# pyproject.toml
version = "0.5.0"
```

```python
# phyloai/__init__.py
__version__ = "0.5.0"
```

- [ ] **Step 6: Run documentation and version checks**

Run:

```bash
pytest tests/cli/test_main.py -k version -v
python - <<'PY'
from phyloai import __version__
assert __version__ == "0.5.0"
print(__version__)
PY
rg -n -- '--site-freq-file|custom exchangeability|自定义交换率|0\.5\.0' \
  docs/commands/tree-ml-iqtree.md docs/commands/tree-ml-iqtree.zh.md \
  README.md README.zh.md skills/phyloai-workflow pyproject.toml phyloai/__init__.py
```

Expected: version test PASS, Python prints `0.5.0`, and every documentation surface has the intended terms.

- [ ] **Step 7: Commit docs and version**

```bash
git add docs/commands/tree-ml-iqtree.md docs/commands/tree-ml-iqtree.zh.md \
  README.md README.zh.md skills/phyloai-workflow/SKILL.md \
  skills/phyloai-workflow/references/parameter-annotations.md \
  pyproject.toml phyloai/__init__.py tests/cli/test_main.py
git commit -m "docs: document custom IQ-TREE profile inputs"
```

## Task 5: Run the complete regression suite and real CAT-PMSF acceptance analysis

**Files:**
- Verify: `runs/test/matrix.fa`
- Verify: `runs/test/chain1.exchangeabilities`
- Verify: `runs/test/chain1.sitefreq`
- Verify: repository test suite

**Interfaces:**
- Consumes: all completed implementation, local `iqtree3`, and user-provided ignored fixtures under `runs/test/`.
- Produces: an IQ-TREE run at `runs/test/cat-pmsf/` whose log proves both absolute custom files were read.

- [ ] **Step 1: Run static and full automated checks**

Run:

```bash
git diff --check
pytest -q
phyloai --version
phyloai tree ml iqtree --help | rg -- '--site-freq-file|Heterogeneous|custom'
```

Expected: no whitespace errors; pytest passes; CLI prints `0.5.0`; help contains the new Heterogeneous option and wording.

- [ ] **Step 2: Execute the agreed real acceptance analysis**

First remove only the test output directory so the command does not need `--overwrite`:

```bash
rm -rf runs/test/cat-pmsf
phyloai tree ml iqtree \
  --matrix runs/test/matrix.fa \
  --model runs/test/chain1.exchangeabilities \
  --site-freq-file runs/test/chain1.sitefreq \
  --state-freq none \
  --rate-heterogeneity +R4 \
  --boot 0 \
  --threads 1 \
  --output-dir runs/test/cat-pmsf
```

Expected: exit code 0 and `runs/test/cat-pmsf/result.json` exists. This command reads user-provided inputs but does not change them.

- [ ] **Step 3: Verify command provenance and IQ-TREE reads**

Run:

```bash
python - <<'PY'
import json
from pathlib import Path

result = json.loads(Path("runs/test/cat-pmsf/result.json").read_text())
params = result["params"]
cmd = result["data"]["cmd"]
assert params["model"] == str(Path("runs/test/chain1.exchangeabilities").resolve())
assert params["site_freq_file"] == str(Path("runs/test/chain1.sitefreq").resolve())
assert "-fs" in cmd
assert cmd[cmd.index("-fs") + 1] == params["site_freq_file"]
assert "+F" not in cmd[cmd.index("-m") + 1]
print("structured CAT-PMSF provenance verified")
PY
rg -n 'Reading site-specific state frequency file|Reading model parameters from file' \
  runs/test/cat-pmsf/*.log
```

Expected: provenance assertion prints success; IQ-TREE log shows both reads and their absolute paths.

- [ ] **Step 4: Test raw `-fs` override in a dry run**

Run:

```bash
phyloai tree ml iqtree \
  --matrix runs/test/matrix.fa \
  --model runs/test/chain1.exchangeabilities \
  --site-freq-file runs/test/chain1.sitefreq \
  --state-freq none \
  --rate-heterogeneity +R4 \
  --boot 0 --threads 1 --dry-run \
  --tool-args "-fs $(pwd)/runs/test/chain1.sitefreq" \
  --output-dir runs/test/cat-pmsf-override
```

Expected: printed IQ-TREE command has exactly one `-fs`, using the raw override path. Do not execute a second full inference.

- [ ] **Step 5: Review changed files and commit the design update if still uncommitted**

Run:

```bash
git status --short
git diff --check
git diff --stat
```

Confirm ignored `runs/test/cat-pmsf/` output is not staged. If the already-updated design file is still uncommitted, commit only that file:

```bash
git add docs/superpowers/specs/2026-06-19-phyloai-tree-ml-iqtree-design.md
git commit -m "docs: update IQ-TREE custom model design"
```

If it was committed earlier, make no empty commit.
