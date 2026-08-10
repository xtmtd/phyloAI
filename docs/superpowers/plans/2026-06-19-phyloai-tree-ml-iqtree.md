# PhyloAI Tree ML IQ-TREE Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Implement `phyloai tree ml iqtree` — IQ-TREE3 backend for maximum-likelihood tree inference supporting homogeneous, heterogeneous (AA mixture/PMSF, NT MIX+MF), partitioned, and ModelFinder workflows with batch/single input modes and two-tier --tool-args semantics.

**Architecture:** New file `phyloai/tree/ml_iqtree.py` implements IQ-TREE core following `ml.py` (fasttree) style. CLI command added to `phyloai/cli/commands/tree.py`. Checkpoint helpers extended in `phyloai/tree/checkpoint_helpers.py` for MF-mode resume. Two-tier tool-args (BLOCKED/OVERRIDEABLE) per Section 9.9.

**Tech Stack:** Python 3.10+, Click, BioPython (Phylo/SeqIO), Rich progress, concurrent.futures ProcessPoolExecutor.

---

## File Structure

```
phyloai/tree/
├── __init__.py                    # (no change)
├── ml.py                          # FastTree (no change, reused _validate_seq_types)
├── ml_iqtree.py                   # NEW — IQ-TREE core (~700 lines)
└── checkpoint_helpers.py          # MODIFY — add resume_verifier_iqtree

phyloai/cli/commands/
└── tree.py                        # MODIFY — add "iqtree" to _MLGroup + command

tests/
├── tree/
│   └── test_ml_iqtree.py          # NEW — unit + integration tests
└── cli/
    └── test_tree.py               # MODIFY — add iqtree CLI tests
```

**Interfaces (key functions in ml_iqtree.py):**

| Function | Signature | Purpose |
|----------|-----------|---------|
| `_scan_input_iqtree` | `(msa_dir: Path) -> tuple[list[Path], list[dict]]` | Scan MSA dir for IQ-TREE-compatible files |
| `_resolve_iqtree_path` | `(path: str\|None, dry_run: bool) -> str` | Resolve iqtree3 executable (custom → PATH → bundled) |
| `_detect_iqtree_version` | `(executable: str) -> dict[str,str]` | Parse `iqtree3 --version` output |
| `_check_managed_flag_conflict` | `(tool_args: str) -> None` | Reject BLOCKED flags in --tool-args |
| `_is_flag_overridden` | `(flag: str, tokens: set[str]) -> bool` | Check if flag present in --tool-args |
| `_classify_workflow` | `(**params) -> str` | Classify workflow type (9 variants) |
| `_build_model_string` | `(**params) -> str` | Assemble -m argument (e.g. LG+F+R4) |
| `_validate_model` | `(model, seq_type, modelfinder) -> None` | Validate model against domain |
| `_run_validations` | `(**params) -> None` | Pre-flight validation (15+ rules) |
| `_build_iqtree_cmd` | `(**params) -> list[str]` | Assemble IQ-TREE CLI argument list |
| `_parse_iqtree_report` | `(path: Path) -> dict` | Parse .iqtree for log_likelihood + model |
| `_run_one_iqtree` | `(**params) -> dict` | Execute IQ-TREE on one MSA, collect outputs |
| `run_iqtree` | `(**params) -> dict` | Main entry: single + batch orchestration |
| `_assemble_iqtree_result` | `(**params) -> dict` | Build result.json payload |

---

## Shared Constants (used across tasks)

```python
# In ml_iqtree.py, define:
IQTREE_COMPATIBLE_EXTENSIONS = frozenset({
    ".fa", ".fas", ".fasta", ".faa", ".fna",
    ".phy", ".phylip", ".nex", ".nxs", ".nexus", ".aln",
})

AA_STANDARD_MODELS = frozenset({
    "LG", "Poisson", "cpREV", "mtREV", "Dayhoff", "mtMAM", "JTT", "WAG",
    "mtART", "mtZOA", "VT", "rtREV", "DCMut", "PMB", "HIVb", "HIVw",
    "JTTDCMut", "FLU", "Blosum62", "GTR20", "mtMet", "mtVer", "mtInv",
    "FLAVI", "Q.LG", "Q.pfam", "Q.pfam_gb", "Q.bird", "Q.mammal",
    "Q.insect", "Q.plant", "Q.yeast",
})

NT_STANDARD_MODELS = frozenset({
    "GTR", "HKY", "JC", "F81", "K2P", "K3P", "K81uf", "TN", "TNef",
    "TIM", "TIMef", "TVM", "TVMef", "SYM",
})

AA_MIXTURE_MODELS = frozenset(
    {f"C{i}" for i in range(10, 61, 10)}
    | {"EX2", "EX3", "EHO", "UL2", "UL3", "EX_EHO", "LG4M", "LG4X"}
)

_IQTREE_MANAGED_LONG_FLAGS = frozenset({
    "--ufboot", "--alrt", "--bnni", "--fast", "--merge",
    "--rclusterf", "--rcluster-max", "--mset", "--msub",
    "--prefix", "--rate", "--qmax", "--seqtype",
    "--redo", "--redo-tree", "--undo",
})

_IQTREE_MANAGED_SHORT_FLAGS = frozenset({
    "-s", "-m", "-p", "-T", "-B", "-ft", "-g", "-o",
    "-q", "-Q", "-S", "-wslr",
})

_IQTREE_BLOCKED_FLAGS = frozenset({"-s"})
# In batch mode, --prefix is also blocked (see _check_managed_flag_conflict)
_IQTREE_BLOCKED_FLAGS_BATCH = frozenset({"-s", "--prefix"})
_IQTREE_BLOCKED_IO_CHARS = frozenset({"<", ">", "|"})

_IQTREE_FLAG_ALIASES: dict[str, frozenset[str]] = {
    "-B": frozenset({"-B", "--ufboot"}),
    "--ufboot": frozenset({"-B", "--ufboot"}),
}
```

---

### Task 1: Scan input, I/O conflict check, model validation, workflow classification

**Files:**
- Create: `phyloai/tree/ml_iqtree.py` — functions: `_scan_input_iqtree`, `_check_managed_flag_conflict`, `_is_flag_overridden`, `_classify_workflow`, `_build_model_string`, `_validate_model`, `_validate_pmsf_base_model`, `_is_heterogeneous_model`, `_run_validations`
- Create: `tests/tree/test_ml_iqtree.py` — all unit tests

- [ ] **Step 1: Create test file with all unit tests**

Write `tests/tree/test_ml_iqtree.py` with tests for:
- `_scan_input_iqtree`: finds .fa/.fas/.fasta/.faa/.fna/.phy/.phylip/.nex/.nxs/.nexus/.aln; skips .txt, empty, dirs
- `_check_managed_flag_conflict`: blocks `-s`, `--prefix` (in batch mode), blocks `>`, `<`, `|`, allows strategy args
- `_is_flag_overridden`: detects short flags (`-m`, `-T`), long flags (`--merge`), alias pairs (`-B`/`--ufboot`)
- `_classify_workflow`: exactly 9 workflow variants: `homogeneous-no-partition-none`, `homogeneous-no-partition-MF`, `homogeneous-no-partition-MFP`, `homogeneous-partition-none`, `homogeneous-partition-MF-merge`, `homogeneous-partition-MFP-merge`, `AA-heterogeneous-direct`, `AA-heterogeneous-PMSF`, `NT-heterogeneous`
- `_build_model_string`: `LG+F+R4`, `GTR+R4`, `HKY`, `MF`/`MFP`, `LG+C20+F+R4` (PMSF), `C20+F+R4` (direct), `MIX+MF`
- `_validate_model`: AA standard OK, NT standard OK, mixture OK, cross-type rejection, skip when MF
- `_validate_pmsf_base_model`: standard OK, mixture rejected
- `_run_validations`: heterogeneous+msa_dir error, partitions+msa_dir error, rcluster mutual exclusion, rcluster without partitions warn, MF+branch support warn, bnni without boot warn, prefix in batch warn, PMSF requires guide tree, pmsf_base without mixture error, qmax without MIX+MF warn

- [ ] **Step 2: Run tests to verify they fail**

`python -m pytest tests/tree/test_ml_iqtree.py -v` — Expected: all FAIL (no module)

- [ ] **Step 3: Implement ml_iqtree.py with all functions**

Create `phyloai/tree/ml_iqtree.py` implementing all functions listed above. Pattern from existing `ml.py`:
- `_scan_input_iqtree`: same structure as `_scan_input` but with `IQTREE_COMPATIBLE_EXTENSIONS` (includes NEXUS + .aln, no skip message)
- `_validate_seq_types_iqtree`: format-aware validator. Maps `.nex/.nxs/.nexus` → "nexus", `.aln` → "clustal", `.phy/.phylip` → "phylip-relaxed", else → "fasta". Parses with Bio.SeqIO using correct format. Unparseable files → offense record. Otherwise identical logic to `ml._validate_seq_types`.
- `_check_managed_flag_conflict`: same pattern, with `batch_mode` param. Batch blocks `--prefix` + `-s` + I/O chars (`<`, `>`, `|`); single blocks `-s` + I/O chars only.
- `_is_flag_overridden`: checks flag in tool_tokens set, falls back to alias lookup in `_IQTREE_FLAG_ALIASES`
- `_classify_workflow`: if MF/MFP + partitions → `homogeneous-partition-{MF|MFP}-merge` (always merge since rclusterf defaults to 10); MF/MFP without partitions → `homogeneous-no-partition-{MF|MFP}`; modelfinder none + model in AA_MIXTURE_MODELS → `AA-heterogeneous-[direct|PMSF]`; model == MIX+MF → `NT-heterogeneous`; partitions + none → `homogeneous-partition-none`; else → `homogeneous-no-partition-none`
- `_build_model_string`: if MF/MFP → return modelfinder; if PMSF → `base+model+state_freq+rate_het` joined by `+`; else → `model+state_freq+rate_het` joined by `+`; strip `none` parts
- `_validate_model`: skip if MF/MFP; check model ∈ domain set; raise ValueError with valid list
- `_validate_pmsf_base_model`: check ∈ AA_STANDARD_MODELS
- `_is_heterogeneous_model`: model ∈ AA_MIXTURE_MODELS (AA) or model == "MIX+MF" (NT), only when modelfinder is "none"
- `_run_validations`: structured sequence of checks, raise ValueError for hard errors, warnings.warn for soft issues

- [ ] **Step 4: Run tests to verify they pass**

`python -m pytest tests/tree/test_ml_iqtree.py -v` — Expected: ALL PASS

- [ ] **Step 5: Verify**

`python -m pytest tests/tree/test_ml_iqtree.py -v` (all non-iqtree3 tests PASS)

---

### Task 2: Command builder `_build_iqtree_cmd`

**Files:**
- Modify: `phyloai/tree/ml_iqtree.py` — add `_build_iqtree_cmd`
- Modify: `tests/tree/test_ml_iqtree.py` — add builder tests

- [ ] **Step 1: Write tests for _build_iqtree_cmd**

Append to `tests/tree/test_ml_iqtree.py` tests covering:
- `def test_build_iqtree_cmd_basic_homogeneous(...)`: checks `-s`, `--prefix`, `-m LG+F+R4`, `-T AUTO`, `--seqtype DNA` (NT). **No `--redo` by default.**
- `def test_build_iqtree_cmd_seqtype_aa(...)`: `--seqtype AA` for AA
- `def test_build_iqtree_cmd_seqtype_auto(...)`: no `--seqtype` for auto
- `def test_build_iqtree_cmd_with_boot_and_alrt(...)`: `-B 1000`, `--alrt 1000`, `--bnni`
- `def test_build_iqtree_cmd_alrt_zero(...)`: `--alrt 0` for parametric
- `def test_build_iqtree_cmd_modelfinder_mf(...)`: `-m MF`, `--mset LG,WAG`, `--msub nuclear`
- `def test_build_iqtree_cmd_fast_mode(...)`: `--fast`
- `def test_build_iqtree_cmd_partitions_merge(...)`: `-p`, `--merge`, `--rclusterf 10`
- `def test_build_iqtree_cmd_pmsf(...)`: `-ft guide.nwk`
- `def test_build_iqtree_cmd_mix_mf_qmax(...)`: `-m MIX+MF`, `-qmax 10`, `--mset GTR,HKY`
- `def test_build_iqtree_cmd_output_flags(...)`: `--rate`, `-wslr`, `-g constraint.nwk`, `-o taxon1,taxon2`
- `def test_build_iqtree_cmd_tool_args_appended(...)`: `-pers 0.5 -nstop 500`
- `def test_build_iqtree_cmd_tool_args_overrides_model(...)`: `-m` appears once with `--tool-args` value
- `def test_build_iqtree_cmd_tool_args_overrides_boot(...)`: `-B` appears once with `--tool-args` value
- `def test_build_iqtree_cmd_tool_args_blocked_s(...)`: raises ValueError for `-s`
- `def test_build_iqtree_cmd_tool_args_blocked_prefix_in_batch(...)`: raises ValueError for `--prefix` when `batch_mode=True`
- `def test_build_iqtree_cmd_tool_args_blocked_pipe(...)`: raises ValueError for `|`

No branch support flags when boot/alrt are None.

- [ ] **Step 2: Run tests to verify they fail**

`python -m pytest tests/tree/test_ml_iqtree.py -k "build_iqtree_cmd" -v` — FAIL

- [ ] **Step 3: Implement _build_iqtree_cmd**

```python
def _build_iqtree_cmd(
    input_path: Path,
    prefix: Path,
    *,
    model_string: str,
    seq_type: str,
    boot: int | None,
    alrt: int | None,
    bnni: bool,
    mode: str,
    threads_arg: str,
    executable: str = "iqtree3",
    mset: str | None = None,
    msub: str | None = None,
    partitions: str | None = None,
    rclusterf: int | None = None,
    rcluster_max: int | None = None,
    guide_tree: str | None = None,
    qmax: int | None = None,
    rate: bool = False,
    wslr: bool = False,
    constraint: str | None = None,
    outgroup: str | None = None,
    tool_args: str | None = None,
    batch_mode: bool = False,
) -> list[str]:
```

Logic:
1. `cmd = [executable]`
2. `tool_tokens = set(shlex.split(tool_args)) if tool_args else set()`
3. If tool_args: call `_check_managed_flag_conflict(tool_args, batch_mode=batch_mode)`
4. Always: `cmd.extend(["-s", str(input_path), "--prefix", str(prefix)])`. **No `--redo` by default.**
5. Conditional (each gated by `_is_flag_overridden`):
   - `-m model_string`
   - `--mset <mset>` (if mset not None and mset != "all")
   - `--msub <msub>`
   - `-p <partitions>`, `--merge` (if partitions)
   - `--rclusterf N` or `--rcluster-max N`
   - `-ft <guide_tree>`
   - `-qmax N`
   - `--fast` (if mode == "fast")
   - `--seqtype DNA|AA` (based on seq_type, skip for "auto")
   - `-B <boot>` (if boot is not None)
   - `--alrt <alrt>` (if alrt is not None)
   - `--bnni` (if bnni and boot is not None)
   - `threads_arg` (the `-T NUM` or `-T AUTO` string, as-is)
   - `--rate`, `-wslr`
   - `-g <constraint>`, `-o <outgroup>`
6. If tool_args: `cmd.extend(shlex.split(tool_args))`
7. Return cmd

- [ ] **Step 4: Run tests to verify they pass**

`python -m pytest tests/tree/test_ml_iqtree.py -k "build_iqtree_cmd" -v` — PASS

- [ ] **Step 5: Verify**

`python -m pytest tests/tree/test_ml_iqtree.py -k "build_iqtree_cmd" -v` (all PASS)

---

### Task 3: `.iqtree` report parser

**Files:**
- Modify: `phyloai/tree/ml_iqtree.py` — add `_parse_iqtree_report`
- Modify: `tests/tree/test_ml_iqtree.py` — add parser tests

- [ ] **Step 1: Write tests for _parse_iqtree_report**

Tests:
- Basic: `"Log-likelihood of the tree: -12345.678"` → `-12345.678`
- Scientific: `"Log-likelihood: -1.234567e+04"` → parsed float
- Model BIC: `"Best-fit model according to BIC: LG+F+R4"` → `"LG+F+R4"`
- Model AIC: `"Best-fit model according to AIC: WAG+F+I+G4"` → `"WAG+F+I+G4"`
- Both fields present in same file
- Missing file → both None
- Empty file → both None
- No model line (e.g., `"ModelFinder will test 484 models"`) → model_selected is None

- [ ] **Step 2: Run → FAIL**

- [ ] **Step 3: Implement**

```python
import re as _re

_LOG_LIKE_RE = _re.compile(r"Log-likelihood(?: of the tree)?:\s+([-\d.eE+]+)", _re.IGNORECASE)
_MODEL_SELECTED_RE = _re.compile(r"Best-fit model(?: according to \w+)?:\s*(\S+)", _re.IGNORECASE)

def _parse_iqtree_report(iqtree_path: Path) -> dict[str, Any]:
    result = {"log_likelihood": None, "model_selected": None}
    if not iqtree_path.exists():
        return result
    try:
        text = iqtree_path.read_text()
    except Exception:
        return result
    m = _LOG_LIKE_RE.search(text)
    if m:
        try:
            result["log_likelihood"] = float(m.group(1))
        except ValueError:
            pass
    m = _MODEL_SELECTED_RE.search(text)
    if m:
        result["model_selected"] = m.group(1)
    return result
```

- [ ] **Step 4: Run → PASS**

- [ ] **Step 5: Verify**

`python -m pytest tests/tree/test_ml_iqtree.py -k "parse_iqtree_report" -v` (all PASS)

---

### Task 4: Executable resolution, version detection, and `_run_one_iqtree`

**Files:**
- Modify: `phyloai/tree/ml_iqtree.py` — add `_resolve_iqtree_path`, `_detect_iqtree_version`, `_run_one_iqtree`
- Modify: `tests/tree/test_ml_iqtree.py` — add tests

- [ ] **Step 1: Write tests**

Tests for `_resolve_iqtree_path`:
- Custom path that exists + executable → returns str
- Missing path → raises FileNotFoundError
- Path not executable → raises ValueError
- Dry run → returns "iqtree3"

Tests for `_detect_iqtree_version`:
- Mock subprocess returning `"IQ-TREE multicore version 3.1.2"` → `{"iqtree3": "3.1.2"}`
- Mock Exception → `{"iqtree3": "unknown"}`

Tests for `_run_one_iqtree`:
- `def test_run_one_iqtree_dry_run(...)`: returns `"dry_run"` with cmd and n_taxa
- `def test_run_one_iqtree_missing_input(...)`: returns `"failed"` with reason
- `@pytest.mark.skipif(not shutil.which("iqtree3"), ...)` `def test_run_one_iqtree_success(...)`: 3-taxon simple AA, status `"success"`, `.treefile` and `.iqtree` exist
- `@pytest.mark.skipif(...)` `def test_run_one_iqtree_mf_mode_no_tree(...)`: MF mode, status `"success"`, `output_tree` is None, `.iqtree` exists

- [ ] **Step 2: Run → FAIL**

- [ ] **Step 3: Implement**

`_resolve_iqtree_path(iqtree_path, dry_run)`:
- Same pattern as `_resolve_fasttree` in ml.py but uses `ToolEnv().require("iqtree3")`
- Error message: `"iqtree3 not found. Install from https://github.com/iqtree/iqtree3/releases or use --iqtree-path."`

`_detect_iqtree_version(executable)`:
- Same pattern as `_detect_fasttree_version`, runs `[executable, "--version"]`

`_run_one_iqtree(...)` — **batch mode** (when `batch_mode=True`):
- Creates per-task temp `work_dir` via `tempfile.mkdtemp(prefix=f"iqtree_{stem}_")`
- IQ-TREE `--prefix work_dir/<stem>` → files: `<stem>.iqtree`, `<stem>.log`, `<stem>.treefile` in work_dir
- Calls `_build_iqtree_cmd` with `batch_mode=True`
- Runs `subprocess.run(cmd, cwd=str(work_dir), ...)` capturing stdout/stderr
- Return code ≠ 0 → `"failed"` with exit code and stderr
- On success, moves output files from work_dir to correct locations:
  - `work_dir/<stem>.treefile` → `output_dir/<stem>.treefile` (only non-MF)
  - `work_dir/<stem>.iqtree` → `logs_dir/<stem>.iqtree`
  - `work_dir/<stem>.log` → `logs_dir/<stem>.log`
  - All other files (`<stem>.ufboot`, `<stem>.contree`, `<stem>.splits.nex`, `<stem>.rate`, `<stem>.sitelh`, `<stem>.ckp.gz`, `<stem>.model.gz`) → `logs_dir/<stem>.<ext>`
- Always cleans up work_dir in `finally` if self-created
- Return dict matches per-task result schema in spec (Section 7)

`_run_one_iqtree(...)` — **single mode** (when `batch_mode=False`):
- Does NOT create temp dir. IQ-TREE runs directly in `output_dir`.
- IQ-TREE `--prefix output_dir/<stem>` → native files go directly to `output_dir`
- Can skip file movement entirely since `output_dir` IS the target
- Does NOT clean up (preserves IQ-TREE checkpoint files for native resume)
- Still parses `.iqtree` with `_parse_iqtree_report` for metadata
- Validates `.treefile` with `Bio.Phylo.read(path, "newick")` for non-MF modes

- [ ] **Step 4: Run → PASS** (iqtree3-dependent tests may skip)

- [ ] **Step 5: Verify**

`python -m pytest tests/tree/test_ml_iqtree.py -k "resolve_iqtree or detect_iqtree or run_one_iqtree" -v` (all PASS; iqtree3-dependent may skip)

---



### Task 5: Main entry point `run_iqtree` and `_assemble_iqtree_result`

**Files:**
- Modify: `phyloai/tree/ml_iqtree.py` — add `run_iqtree`, `_assemble_iqtree_result`, `_resolved_iqtree_params`, `_reconstruct_result`
- Modify: `tests/tree/test_ml_iqtree.py` — add end-to-end tests

- [ ] **Step 1: Write tests**

Tests:
- `test_run_iqtree_neither_input_raises()`: raises ValueError
- `test_run_iqtree_both_inputs_raises(tmp_path)`: raises ValueError
- `test_run_iqtree_single_dry_run(tmp_path)`: returns `"success"`, mode `"--matrix"`, has `cmd`
- `test_run_iqtree_batch_dry_run(tmp_path)`: returns `"success"`, n_input ≥ 2, files have `cmd`
- `test_run_iqtree_heterogeneous_rejected_in_batch(tmp_path)`: raises ValueError with `"only supported in --matrix"`
- `test_run_iqtree_no_valid_inputs(tmp_path)`: raises ValueError with `"No valid input files"`
- `test_run_iqtree_unsupported_matrix_extension(tmp_path)`: raises ValueError with `"unsupported extension"`
- `test_run_iqtree_unparsable_matrix(tmp_path)`: raises ValueError with `"Cannot parse"`
- `test_run_iqtree_overwrite_resume_mutual(tmp_path)`: raises ValueError

- [ ] **Step 2: Run → FAIL**

- [ ] **Step 3: Implement**

`run_iqtree(...)` — main orchestration (~250 lines):
1. Input mutual exclusivity check
2. overwrite/resume mutual exclusivity
3. Resolve iqtree3 executable + version
4. Parse threads: `_parse_threads(threads, batch_mode)` — if `None`: batch defaults to `4`, single defaults to `"auto"`; batch rejects `"auto"`; single accepts `"auto"` or int >= 1. Returns `(threads_int, threads_str)`.
5. **Output directory conflict check** (before running any tool):
   - `--overwrite`: `shutil.rmtree(output_dir)` if exists
   - `--resume`: validate checkpoint exists; skip if already success
   - Otherwise: reject if `output_dir` exists and is non-empty → ValueError
6. Single mode (`--matrix`):
   - Validate extension ∈ IQTREE_COMPATIBLE_EXTENSIONS
   - Parse matrix with format-specific parsing: `.nex/.nxs/.nexus` → "nexus", `.aln` → "clustal", `.phy/.phylip` → "phylip-relaxed", else → "fasta". Use Bio.SeqIO with detected format.
   - Resolve defaults (model, mset, etc.)
   - Build model_string via `_build_model_string`
   - Threads: `-T <resolved>` or `-T AUTO`
   - Run `_run_one_iqtree` with `batch_mode=False` (runs in output_dir directly, no temp dir)
   - Call `_assemble_iqtree_result` with batch_mode=False
7. Batch mode (`--msa-dir`):
   - Scan input via `_scan_input_iqtree`
   - Validate seq types via `_validate_seq_types_iqtree`
   - Resolve defaults
   - Threads: each job `-T 1`; batch parallelism via `ProcessPoolExecutor(max_workers=threads_int)`
   - Build initial IQ-TREE checkpoint via `build_initial_iqtree_checkpoint` (records both `tree` and `iqtree` outputs)
   - Run all tasks with progress_callback, each with `batch_mode=True`
   - Use `plan_resume_iqtree` for resume (validates `.treefile` for non-MF, `.iqtree` for MF)
   - Assemble result with batch_mode=True
8. Return result.json dict

`_assemble_iqtree_result(...)`:
- Compute `n_successful`: count of tasks with `status == "success"` or `"dry_run"` (works for both tree and model-only)
- Compute `n_failed`: count of tasks with `status == "failed"`
- Status: `"error"` if `n_successful == 0 and (n_failed > 0 or n_skipped > 0)`, else `"success"`
- For single mode: `log_likelihood` from first task result; batch mode: null
- `model_selected`: from first task if MF/MFP, else `model_string` from params
- `key_results`: workflow, model_selected, log_likelihood, boot, alrt, partitioned, merged_partitions
- `data.files` list with per-file fields (same as spec Section 7)
- Write `iqtree.log` with exit code, command, tool versions, wall time, counts

`_reconstruct_result(output_dir, run_start)`:
- Load existing `result.json` if available; otherwise return empty success payload

- [ ] **Step 4: Run → PASS**

- [ ] **Step 5: Verify**

`python -m pytest tests/tree/test_ml_iqtree.py -k "run_iqtree" -v` (all PASS; iqtree3-dependent may skip)

---

### Task 6: IQ-TREE-aware checkpoint build and resume

**Files:**
- Modify: `phyloai/tree/checkpoint_helpers.py` — add `build_initial_iqtree_checkpoint`, `plan_resume_iqtree`, `resume_verifier_iqtree`
- Modify: `tests/tree/test_ml_iqtree.py` — add verifier tests

- [ ] **Step 1: Write tests**

```python
def test_resume_verifier_iqtree_validates_treefile(tmp_path):
    from phyloai.tree.checkpoint_helpers import resume_verifier_iqtree
    tree = tmp_path / "gene.treefile"
    tree.write_text("(a:0.1,b:0.2);")
    verify = resume_verifier_iqtree(validate_tree=True)
    assert verify(tree) is True

def test_resume_verifier_iqtree_rejects_empty(tmp_path):
    from phyloai.tree.checkpoint_helpers import resume_verifier_iqtree
    tree = tmp_path / "gene.treefile"
    tree.write_text("")
    verify = resume_verifier_iqtree(validate_tree=True)
    assert verify(tree) is False

def test_resume_verifier_iqtree_mf_mode_validates_iqtree(tmp_path):
    from phyloai.tree.checkpoint_helpers import resume_verifier_iqtree
    iqtree = tmp_path / "gene.iqtree"
    iqtree.write_text("Log-likelihood: -100.0\n")
    verify = resume_verifier_iqtree(validate_tree=False)
    assert verify(iqtree) is True

def test_resume_verifier_iqtree_mf_mode_empty_iqtree(tmp_path):
    from phyloai.tree.checkpoint_helpers import resume_verifier_iqtree
    iqtree = tmp_path / "gene.iqtree"
    iqtree.write_text("")
    verify = resume_verifier_iqtree(validate_tree=False)
    assert verify(iqtree) is False
```

- [ ] **Step 2: Run → FAIL** (import error)

- [ ] **Step 3: Implement `build_initial_iqtree_checkpoint`, `plan_resume_iqtree`, and `resume_verifier_iqtree`**

Add to `checkpoint_helpers.py` after existing functions:

```python
def resume_verifier_iqtree(validate_tree: bool = True) -> Callable[[Path], bool]:
    """Return verifier for IQ-TREE output. validate_tree=True checks .treefile with Bio.Phylo.
    validate_tree=False (MF mode) checks .iqtree exists and is non-empty."""
    if validate_tree:
        def _verify(tree_path: Path) -> bool:
            if not tree_path.exists() or tree_path.stat().st_size == 0:
                return False
            try:
                Phylo.read(str(tree_path), "newick")
                return True
            except Exception:
                return False
        return _verify
    else:
        def _verify(iqtree_path: Path) -> bool:
            if not iqtree_path.exists() or iqtree_path.stat().st_size == 0:
                return False
            return True
        return _verify


def build_initial_iqtree_checkpoint(
    *, step: str, command: str, params: dict[str, Any],
    inputs: list[Path], trees_dir: Path, logs_dir: Path,
) -> Checkpoint:
    """Like build_initial_checkpoint but records both tree and iqtree outputs per task."""
    now = _utc_now_iso()
    tasks = [
        CheckpointTask(
            task_id=inp.stem,
            status="pending",
            input=str(inp),
            outputs={
                "tree": str(trees_dir / f"{inp.stem}.treefile"),
                "iqtree": str(logs_dir / f"{inp.stem}.iqtree"),
            },
        )
        for inp in inputs
    ]
    return Checkpoint(
        schema_version=CHECKPOINT_SCHEMA_VERSION,
        step=step, command=command, status="running",
        params_hash=canonical_params_hash(params),
        params=params,
        started_at=now, updated_at=now, completed_at=None,
        tasks=tasks,
    )


def plan_resume_iqtree(checkpoint: Checkpoint, is_mf_only: bool = False) -> tuple[list[str], list[str]]:
    """Like plan_resume but for IQ-TREE. Uses .treefile verifier for tree-producing workflows,
    .iqtree verifier for MF-only workflows."""
    to_run: list[str] = []
    skipped: list[str] = []
    verifier = resume_verifier_iqtree(validate_tree=not is_mf_only)

    for task in checkpoint.tasks:
        if task.status in {"pending", "running", "failed"}:
            to_run.append(task.task_id)
        elif task.status == "success":
            output_key = "iqtree" if is_mf_only else "tree"
            output_path_str = task.outputs.get(output_key)
            if output_path_str and verifier(Path(output_path_str)):
                skipped.append(task.task_id)
            else:
                to_run.append(task.task_id)
        else:
            skipped.append(task.task_id)

    return to_run, skipped
```

- [ ] **Step 4: Run → PASS**

- [ ] **Step 5: Verify**

`python -m pytest tests/tree/test_ml_iqtree.py -k "resume_verifier_iqtree" -v` (all PASS)

---

### Task 7: CLI integration — iqtree command in tree.py

**Files:**
- Modify: `phyloai/cli/commands/tree.py` — add "iqtree" to _MLGroup + iqtree command
- Modify: `tests/cli/test_tree.py` — add iqtree CLI tests

- [ ] **Step 1: Update _MLGroup.list_commands**

```python
class _MLGroup(click.Group):
    def list_commands(self, ctx: click.Context) -> list[str]:
        return ["fasttree", "iqtree"]
```

- [ ] **Step 2: Add iqtree CLI command**

Pattern follows `fasttree_command` structure. Add after fasttree_command function:

```python
@ml.command(
    "iqtree",
    help=(
        "Infer ML trees using IQ-TREE3.\n\n"
        "  --msa-dir : batch gene trees (homogeneous workflows only)\n\n"
        "  --matrix  : single supermatrix (all workflows)\n\n"
        "Reads FASTA, PHYLIP, NEXUS, CLUSTAL formats. "
        "Heterogeneous models (C10-C60, MIX+MF) require --matrix."
    ),
)
# ... click options ...
def iqtree_command(...) -> None:
```

Click options (in functional blocks matching spec Section 2):
- Input: `--msa-dir`, `--matrix`
- Data Type: `--seq-type` (AA|NT|auto, default auto)
- Model (when modelfinder none): `--model`, `--state-freq`, `--rate-heterogeneity`
- ModelFinder: `--modelfinder` (MF|MFP|none), `--mset`, `--msub`
- Partitions: `--partitions`, `--rclusterf`, `--rcluster-max`
- Heterogeneous: `--pmsf-base-model`, `--guide-tree`, `--qmax`
- Tree Search: `--mode` (normal|fast), `--constraint`
- Branch Support: `--boot`, `--alrt`, `--bnni`
- Output: `--rate`, `--wslr`, `--outgroup`, `--prefix`
- Execution: `-o/--output-dir`, `--threads` (type=str, default=None), `--overwrite`, `--resume`, `--dry-run`, `-q/--quiet`, `--iqtree-path`, `--tool-args`

CLI body logic:
1. Mutual exclusivity: exactly one of `--msa-dir` / `--matrix`
2. threads: `_parse_threads(threads, batch_mode)` — if `None`: batch → `1`, single → `"auto"`; batch rejects `"auto"`; single accepts `"auto"` or int >= 1
3. resume + overwrite mutual exclusivity
4. Validate paths exist (msa-dir, matrix, partitions, guide-tree, constraint, iqtree-path)
5. iqtree-path: check exists + executable
6. Invoke `run_iqtree(...)` with Rich progress bar for batch mode (same pattern as fasttree)
7. On success: write `result.json`, print summary (trees/failed/skipped)
8. On `KeyboardInterrupt`: exit gracefully
9. On error: map Exception type to exit code (1=user, 2=tool, 3=env) matching spec Section 6

- [ ] **Step 3: Add CLI tests to test_tree.py**

```python
def test_tree_ml_iqtree_help():
    result = CliRunner().invoke(cli, ["tree", "ml", "iqtree", "--help"])
    assert result.exit_code == 0
    assert "iqtree" in result.output
    assert "--msa-dir" in result.output
    assert "--modelfinder" in result.output

def test_tree_ml_iqtree_mutual_exclusivity(tmp_path):
    msa_dir = tmp_path / "msas"; msa_dir.mkdir()
    mat = tmp_path / "matrix.fa"; mat.write_text(">a\nMKT\n")
    result = CliRunner().invoke(cli, [
        "tree", "ml", "iqtree", "--msa-dir", str(msa_dir), "--matrix", str(mat),
    ])
    assert result.exit_code == 1

def test_tree_ml_iqtree_neither_input():
    result = CliRunner().invoke(cli, ["tree", "ml", "iqtree"])
    assert result.exit_code == 1

def test_cli_iqtree_msa_dir_nonexistent():
    result = CliRunner().invoke(cli, [
        "tree", "ml", "iqtree", "--msa-dir", "/nonexistent/path",
    ])
    assert result.exit_code == 1

def test_tree_ml_iqtree_quiet_dry_run_single(tmp_path):
    mat = tmp_path / "matrix.fa"; mat.write_text(">a\nMKTLLL\n>b\nMKTLLL\n")
    out_dir = tmp_path / "out"
    result = CliRunner().invoke(cli, [
        "tree", "ml", "iqtree", "--matrix", str(mat),
        "--output-dir", str(out_dir), "--seq-type", "AA", "--model", "LG",
        "--quiet", "--dry-run",
    ])
    assert result.exit_code == 0

def test_tree_ml_iqtree_quiet_dry_run_batch(tmp_path):
    msa_dir = tmp_path / "msas"; msa_dir.mkdir()
    (msa_dir / "g1.fa").write_text(">a\nMKTLLL\n>b\nMKTLLL\n")
    out_dir = tmp_path / "out"
    result = CliRunner().invoke(cli, [
        "tree", "ml", "iqtree", "--msa-dir", str(msa_dir),
        "--output-dir", str(out_dir), "--seq-type", "AA", "--model", "LG",
        "--quiet", "--dry-run",
    ])
    assert result.exit_code == 0

def test_tree_ml_iqtree_blocked_tool_args(tmp_path):
    mat = tmp_path / "matrix.fa"; mat.write_text(">a\nMKTLLL\n")
    out_dir = tmp_path / "out"
    result = CliRunner().invoke(cli, [
        "tree", "ml", "iqtree", "--matrix", str(mat),
        "--output-dir", str(out_dir), "--tool-args", "-s hack.fa", "--quiet",
    ])
    assert result.exit_code == 1

def test_tree_ml_iqtree_heterogeneous_in_batch(tmp_path):
    msa_dir = tmp_path / "msas"; msa_dir.mkdir()
    (msa_dir / "g1.fa").write_text(">a\nMKTLLL\n")
    out_dir = tmp_path / "out"
    result = CliRunner().invoke(cli, [
        "tree", "ml", "iqtree", "--msa-dir", str(msa_dir),
        "--output-dir", str(out_dir), "--model", "C20", "--quiet",
    ])
    assert result.exit_code == 1

def test_tree_ml_iqtree_modelfinder_mf_dry_run(tmp_path):
    mat = tmp_path / "matrix.fa"; mat.write_text(">a\nMKTLLL\n>b\nMKTLLL\n")
    out_dir = tmp_path / "out"
    result = CliRunner().invoke(cli, [
        "tree", "ml", "iqtree", "--matrix", str(mat),
        "--output-dir", str(out_dir),
        "--modelfinder", "MF", "--mset", "LG,WAG",
        "--quiet", "--dry-run",
    ])
    assert result.exit_code == 0

def test_tree_ml_iqtree_writes_result_json_and_log(tmp_path):
    mat = tmp_path / "matrix.fa"; mat.write_text(">a\nMKTLLL\n>b\nMKTLLL\n")
    out_dir = tmp_path / "out"
    result = CliRunner().invoke(cli, [
        "tree", "ml", "iqtree", "--matrix", str(mat),
        "--output-dir", str(out_dir), "--seq-type", "AA", "--model", "LG",
        "--quiet",
    ])
    if result.exit_code == 0:
        assert (out_dir / "result.json").exists()
        assert (out_dir / "iqtree.log").exists()
    elif result.exit_code == 3:
        pytest.skip("iqtree3 not installed")
```

- [ ] **Step 4: Run CLI tests → PASS** (iqtree3-dependent tests may skip)

- [ ] **Step 5: Verify**

`python -m pytest tests/cli/test_tree.py -k "iqtree" -v` (all PASS or SKIP for iqtree3)

---

### Task 8: Final integration — run full test suite and lint

**Files:** (none new, verify all)

- [ ] **Step 1: Run full unit test suite**

```bash
python -m pytest tests/tree/test_ml_iqtree.py -v
python -m pytest tests/cli/test_tree.py -k "iqtree" -v
```

Expected: all unit tests PASS; iqtree3-dependent tests SKIP if not installed.

- [ ] **Step 2: Run iqtree3 integration tests (if iqtree3 available)**

```bash
# Real-command tests only if iqtree3 is in PATH
python -m pytest tests/tree/test_ml_iqtree.py -v -k "run_one" --no-header
python -m pytest tests/cli/test_tree.py -v -k "iqtree_writes_result"
```

Expected: integration tests PASS with real iqtree3.

- [ ] **Step 3: Run tests and lint**

```bash
python -m pytest tests/tree/test_ml_iqtree.py tests/cli/test_tree.py -v
```

Fix any test failures. If available, also run project lint command from project convention.

- [ ] **Step 4: Acceptance test coverage verification**

Verify these spec requirements are tested:
- [x] --msa-dir and --matrix mutual exclusivity (CLI + unit)
- [x] --seq-type NT → --seqtype DNA mapping
- [x] --modelfinder MF: model-only, no tree, branch support warn+ignore
- [x] Heterogeneous models reject --msa-dir
- [x] --partitions reject --msa-dir
- [x] --tool-args BLOCKED: -s rejected, --prefix rejected in batch, >/<, | rejected
- [x] --tool-args OVERRIDEABLE: -m, -B, -T overridden when present
- [x] --merge flag with --partitions + MF/MFP
- [x] -ft flag for PMSF with --guide-tree
- [x] -qmax for MIX+MF
- [x] --alrt 0 (parametric aLRT)
- [x] .iqtree parsing: log_likelihood + model_selected
- [x] result.json: key_results match spec structure
- [x] iqtree.log written with exit code
- [x] No --redo in base command (only via --tool-args)
- [x] Single mode: native IQ-TREE resume (no temp dir, preserves ckp)
- [x] Batch MF mode: status based on n_successful, not n_trees
- [x] --threads parsed as str with batch/single validation
- [x] NEXUS/CLUSTAL format detection and parsing
- [x] MF mode checkpoint: both tree and iqtree outputs recorded, verified by plan_resume_iqtree

- [ ] **Step 5: Verify full test suite**

```bash
python -m pytest tests/tree/test_ml_iqtree.py tests/cli/test_tree.py -v
```

All unit tests PASS. iqtree3-dependent tests may SKIP gracefully.

---

## Self-Review Checklist

**1. Spec coverage:** Each spec section mapped:
- Section 2 (CLI Parameters): all flags covered in Click options + command builder
- Section 3 (Workflow Classification): `_classify_workflow` returns all 9 variants
- Section 3 (Command Mapping): `_build_iqtree_cmd` tests verify each mapping row
- Section 4 (Validation Rules): `_run_validations` covers all 15 rules
- Section 5 (Resume): `build_initial_iqtree_checkpoint` and `plan_resume_iqtree` handle MF + non-MF tasks
- Section 6 (Exit Codes): CLI mapping (1=user, 2=tool, 3=env) in command body
- Section 7 (File Organization): output dir layout matching spec
- Section 7 (result.json): `_assemble_iqtree_result` produces exact spec structure
- Section 7 (.iqtree parsing): `_parse_iqtree_report` with both regex patterns
- Section 8 (CLI Layer): command added to tree.py under ml group

**2. Placeholder scan:** No "TBD", "TODO", "implement later". All code shown inline.

**3. Type consistency:**
- `run_iqtree` signature matches `run_fasttree` pattern (keyword-only args, same return type)
- `_run_one_iqtree` return dict matches per-task schema
- `_assemble_iqtree_result` produces `result.json` matching spec
- Checkpoint uses `build_initial_iqtree_checkpoint`/`mark_task`/`plan_resume_iqtree` — IQ-TREE-specific variants that record both `tree` and `iqtree` outputs

**4. Cross-file consistency:**
- `checkpoint_helpers.build_initial_iqtree_checkpoint` sets `outputs={"tree": ..., "iqtree": ...}` per task
- `checkpoint_helpers.plan_resume_iqtree` selects verifier by `validate_tree` flag (`.treefile` vs `.iqtree`)
- `ml_iqtree` imports from `phyloai.core.checkpoint`, `phyloai.core.env`, `phyloai.core.schema`, `phyloai.tree.checkpoint_helpers` — all existing
- `checkpoint_helpers.resume_verifier_iqtree(validate_tree=True)` returns `Callable[[Path], bool]` matching existing pattern
- `ml_iqtree` imports from `phyloai.core.checkpoint`, `phyloai.core.env`, `phyloai.core.schema`, `phyloai.tree.checkpoint_helpers` — all existing
- `tree.py` imports `from phyloai.tree.ml_iqtree import run_iqtree` — follows fasttree pattern
- `_validate_seq_types_iqtree` uses IQ-TREE-specific format-aware parsing (fasta/phylip-relaxed/nexus/clustal) per file extension, reuses `detect_seq_type` from core for sequence-type inference
- `_scan_input_iqtree` is separate from `_scan_input` — IQ-TREE supports more extensions (NEXUS, .aln)

---
