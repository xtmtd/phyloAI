# Posttree Model-Compare Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Implement `phyloai posttree modelcompare iqtree` and `phyloai posttree modelcompare pb` subcommands for relative model comparison and selection.

**Architecture:** Two subcommands backed by separate library modules (`posttree/modelcompare_iqtree.py` and `posttree/modelcompare_pb.py`); CLI wiring in `cli/commands/posttree.py`. `modelcompare_iqtree.py` reuses shared IQ-TREE infrastructure from `phyloai.core.iqtree`. `modelcompare_pb.py` is pure Python with no external dependencies.

**Tech Stack:** Python 3.11+, click, math (stdlib), IQ-TREE3 (external).

## Global Constraints

- All result.json must conform to JSON output standard (`2026-06-21-phyloai-json-output-standard.md`)
- `params` dict must include EVERY parameter from `run_*` function signature; null for unused mutually-exclusive params
- `command` field must be full re-executable CLI string with all resolved defaults
- Reuse `phyloai.core.iqtree._resolve_iqtree_path`, `_detect_iqtree_version`, `IQTREE_COMPATIBLE_EXTENSIONS`
- Reuse `phyloai.core.formats` for seq-type auto-detection
- No new pip dependencies
- Tests live in `tests/posttree/test_modelcompare_iqtree.py`, `tests/posttree/test_modelcompare_pb.py`
- Spec: `docs/superpowers/specs/2026-08-02-phyloai-posttree-modelcompare-design.md`

---

## File Map

| Action | File | Responsibility |
|--------|------|----------------|
| Create | `phyloai/posttree/modelcompare_iqtree.py` | Library: validation, cmd build, het expansion, .iqtree parsing, run_modelcompare_iqtree() |
| Create | `phyloai/posttree/modelcompare_pb.py` | Library: sitelogl parsing, LOO-CV/wAIC math, quality classification, run_modelcompare_pb() |
| Modify | `phyloai/cli/commands/posttree.py` | Add `modelcompare` Click group + `iqtree`/`pb` subcommand wrappers; update _PosttreeGroup.list_commands |
| Modify | `phyloai/report/templates.py` | Add `generate_methods_posttree_modelcompare_iqtree` + `_pb`; register in METHODS_GENERATORS |
| Modify | `phyloai/report/collector.py` | Add `posttree.modelcompare.iqtree` and `posttree.modelcompare.pb` to STEP_ORDER and _THIRD_LEVEL |
| Modify | `docs/superpowers/specs/2026-06-07-phyloai-design.md` | Add modelcompare to posttree command tree / module structure |
| Create | `tests/posttree/test_modelcompare_iqtree.py` | Unit tests |
| Create | `tests/posttree/test_modelcompare_pb.py` | Unit tests |
| Modify | `README.md` | Add modelcompare example commands |
| Modify | `README.zh.md` | Add modelcompare example commands (Chinese) |
| Create | `docs/commands/posttree-modelcompare.md` | English command documentation |
| Create | `docs/commands/posttree-modelcompare.zh.md` | Chinese command documentation |
| Modify | `skills/phyloai-workflow/SKILL.md` | Add modelcompare usage, inputs/outputs, parameter rules |

---

## Task 1: Library `modelcompare_iqtree.py`

**File:** Create `phyloai/posttree/modelcompare_iqtree.py`

- [ ] **Step 1.1:** Module skeleton — imports, constants

```python
"""Model comparison via IQ-TREE ModelFinder (BIC/AIC/AICc)."""
from __future__ import annotations

import json
import re
import shlex
import shutil
import subprocess
import time as _time
from pathlib import Path
from typing import Any

from phyloai.core.iqtree import (
    _detect_iqtree_version,
    _resolve_iqtree_path,
    IQTREE_COMPATIBLE_EXTENSIONS,
)

_BLOCKED_FLAGS = frozenset({"-s", "-m", "-mset", "-mrate", "-madd", "-cmin", "-cmax", "--prefix"})

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

AA_HETEROGENEOUS_MODELS = frozenset({
    "C10", "C20", "C30", "C40", "C50", "C60",
    "EX2", "EX3", "EHO", "UL2", "UL3", "EX_EHO", "LG4M", "LG4X",
})

VALID_MRATE_TOKENS = frozenset({"E", "G", "R"})
VALID_HET_MRATE_TOKENS = frozenset({"E", "G", "R"})
```

- [ ] **Step 1.2:** `_expand_heterogeneous_models(models: list[str], het_mrate: str) -> list[str]`

Expansion logic: each token in het-mrate selects a variant family, mirroring `--mrate`: `E` → M, M+F; `G` → M+G4, M+F+G4; `R` → M+R4, M+F+R4. Only requested families are produced; default `E,G` yields all four. Return flat list.

- [ ] **Step 1.3:** `_detect_seq_type(matrix: Path) -> str`

Use `phyloai.core.formats` to read alignment and detect AA vs NT. Returns "AA" or "NT".

- [ ] **Step 1.4:** `_validate_inputs(...) -> list[str]`

Validation order:
1. Matrix exists + valid extension
2. Detect seq-type (if auto); if explicit `--seq-type`, cross-check against the detected type (mismatch → hard error)
3. Validate homogeneous models against detected type
4. Validate heterogeneous models (AA only)
5. Validate mrate / het-mrate tokens
6. Validate `--threads` is a positive integer or `auto`
7. Validate `--prefix` is a single filename (no `/`, `..`, or absolute paths)
8. Check blocked flags in tool-args (blocked: `-s`, `--prefix`)
9. Check overwrite/resume mutual exclusion

- [ ] **Step 1.5:** `_build_cmd(...) -> list[str]`

Assemble: executable, -s, -m MF, -mset, -mrate, -cmin 4, -cmax 4, [-madd], --prefix, -T, [tool_args]. If `--tool-args` supplies a flag PhyloAI also manages (`-m`, `-mset`, `-mrate`, `-madd`, `-cmin`, `-cmax`, `-T`), skip the PhyloAI-generated one so the flag appears once (tool-args overrides).

- [ ] **Step 1.6:** `_parse_modelfinder_results(iqtree_file: Path) -> list[dict]`

Parse "List of models sorted by BIC scores:" section. For each model row, extract: model, logl, aic, w_aic, in_aic_95, aicc, w_aicc, in_aicc_95, bic, w_bic, in_bic_95. Detect +/- before weight values.

- [ ] **Step 1.7:** `run_modelcompare_iqtree(...) -> dict[str, Any]`

Main entry: validate → build cmd → execute (stdout=None inheriting parent stdout, stderr=PIPE) → parse results → write model_fit.csv → write result.json → return payload. Enforce output-directory conflict: non-empty existing output dir without `--overwrite` → hard error.

---

## Task 2: Library `modelcompare_pb.py`

**File:** Create `phyloai/posttree/modelcompare_pb.py`

- [ ] **Step 2.1:** Module skeleton — imports, Student's t critical values

```python
"""Model comparison via PhyloBayes LOO-CV / wAIC."""
from __future__ import annotations

import json
import math
import shutil
import time as _time
from pathlib import Path
from typing import Any

# Student's t 95% critical values for df 1..30 (from standard tables)
# For df > 30: linear interpolation toward 1.96 (z_0.975).
# Formula: t30 + (1.96 - t30) * (1 - 30/df) — monotonically decreases from t30=2.042 to 1.96.
_STUDENT_T95 = {
    1: 12.706, 2: 4.303, 3: 3.182, 4: 2.776, 5: 2.571,
    6: 2.447, 7: 2.365, 8: 2.306, 9: 2.262, 10: 2.228,
    11: 2.201, 12: 2.179, 13: 2.160, 14: 2.145, 15: 2.131,
    16: 2.120, 17: 2.110, 18: 2.101, 19: 2.093, 20: 2.086,
    21: 2.080, 22: 2.074, 23: 2.069, 24: 2.064, 25: 2.060,
    26: 2.056, 27: 2.052, 28: 2.048, 29: 2.045, 30: 2.042,
}
```

- [ ] **Step 2.2:** `_mean(x)`, `_unbiased_var(x)`, `_student_t95(df)` — basic statistics helpers

- [ ] **Step 2.3:** `_parse_sitelogl(path: Path) -> list[list[float]]`

Parse one .sitelogl file: skip header, return list of rows (each row = [site, logl, var, logcpo, ess, logpostmeanl, ess]). Reject duplicate site identifiers within a single file.

- [ ] **Step 2.4:** `_compute_loocv_waic(runs: list[list[list[float]]]) -> dict`

Integrated readwaic() logic. Returns dict with loocv/waic scores, bias, stdev, CI95, ESS stats.

- [ ] **Step 2.5:** `_classify_quality(pct_ess_lt10: float, frac_ess_lt10: float) -> str`

`max(pct, frac) < 0.1` → "good"; `< 0.3` → "ok"; else → "no"

- [ ] **Step 2.6:** `run_modelcompare_pb(...) -> dict[str, Any]`

Main entry:
1. Validate inputs (mutual exclusion, ≥2 files per group, matching site counts within group)
2. **Site validation:** within each group, all files must share identical ordered `site` identifiers; if ≥2 model groups, all groups must have identical site counts and ordered `site` identifiers; mismatch → hard error
3. **Output-directory conflict:** non-empty existing output dir without `--overwrite` → hard error
4. Resolve model groups and labels (auto model_1/model_2/... for both input modes). With `--model-names`: count must match groups, labels must be unique and be safe path components (no `/`, `..`, or absolute paths) since they name output subdirectories
5. Copy sitelogl files to output/sitelogl/model_N/ (handle duplicate basenames by appending numeric suffix)
6. Compute LOO-CV/wAIC per model
7. If ≥2 models: independently select the best LOO-CV and wAIC models, then compute each metric's Δ relative to its own best
8. Write model_fit.csv (single-model or multi-model format; multi-model includes Pct_ESS_lt10 and Frac_ESS_lt10 columns for both metrics)
9. Write result.json
10. Return payload

---

## Task 3: CLI wiring in `posttree.py`

**File:** Modify `phyloai/cli/commands/posttree.py`

- [ ] **Step 3.1:** Add `_ModelcompareGroup(click.Group)` with `list_commands` returning `["iqtree", "pb"]`

- [ ] **Step 3.2:** Update `_PosttreeGroup.list_commands` to return `["topology", "dating", "signal", "modelcompare"]`

- [ ] **Step 3.3:** Add `@posttree.group("modelcompare", cls=_ModelcompareGroup)` with docstring

- [ ] **Step 3.4:** Add `modelcompare_iqtree_command` with Click options:
  - `--matrix` (required, type=Path)
  - `--homogeneous-model` (required, type=str)
  - `--mrate` (default="E,G")
  - `--heterogeneous-model` (optional)
  - `--het-mrate` (default="E,G")
  - `--seq-type` (choice: AA/NT/auto, default=auto)
  - `--prefix` (default="modelcompare"; must be a single filename)
  - `-o/--output-dir` (default=runs/posttree/modelcompare/iqtree)
  - `--threads` (default="auto")
  - `--iqtree-path`
  - `--tool-args`
  - `--overwrite`, `--resume`, `--dry-run`, `--quiet`
  
  After execution: print model comparison table to terminal.

- [ ] **Step 3.5:** Add `modelcompare_pb_command` with Click options:
  - `--sitelogl-dir` (optional, type=click.Path — comma-separated dirs; enables shell path completion)
  - `--sitelogl` (optional, type=click.Path, multiple=True — each occurrence = one model's chains; enables shell path completion)
  - `--model-names` (optional, type=str — comma-separated labels; each must be a safe path component)
  - `-o/--output-dir` (default=runs/posttree/modelcompare/pb)
  - `--overwrite`, `--quiet`
  
  After execution: print LOO-CV/wAIC table to terminal.

---

## Task 4: Report templates + collector

**File:** Modify `phyloai/report/templates.py`

- [ ] **Step 4.1:** Add `generate_methods_posttree_modelcompare_iqtree(params, key_results, tool_versions)` per spec §9.1
- [ ] **Step 4.2:** Add `generate_methods_posttree_modelcompare_pb(params, key_results, tool_versions)` per spec §9.2
- [ ] **Step 4.3:** Register both in `METHODS_GENERATORS` dict as `"posttree.modelcompare.iqtree"` and `"posttree.modelcompare.pb"`

**File:** Modify `phyloai/report/collector.py`

- [ ] **Step 4.4:** Add `"posttree.modelcompare.iqtree"` and `"posttree.modelcompare.pb"` to STEP_ORDER and _THIRD_LEVEL

---

## Task 5: Tests

- [ ] **Step 5.1:** Create `tests/posttree/test_modelcompare_iqtree.py`
  - Test `_expand_heterogeneous_models("C10", "G")` → `["C10", "C10+F", "C10+G4", "C10+F+G4"]`
  - Test `_expand_heterogeneous_models("C10", "G,R")` → 6 items
  - Test `_expand_heterogeneous_models("C10,C20", "G,R")` → 12 items
  - Test `_validate_inputs` rejects NT + heterogeneous
  - Test `_validate_inputs` rejects invalid mrate tokens
  - Test `_build_cmd` produces correct command list
  - Test `_parse_modelfinder_results` against `runs/modelCompare/EOG090X0A0V.fa.iqtree`
  - Test `run_modelcompare_iqtree` dry-run mode

- [ ] **Step 5.2:** Create `tests/posttree/test_modelcompare_pb.py`
  - Test `_parse_sitelogl` against `runs/modelCompare/LOOCV_wAIC/chain1.sitelogl` (235 rows)
  - Test `_compute_loocv_waic` with 3-chain test data (approximate expected values)
  - Test `_classify_quality(0.05, 0.05)` → "good"
  - Test `_classify_quality(0.15, 0.05)` → "ok"
  - Test `_classify_quality(0.35, 0.05)` → "no"
  - Test `run_modelcompare_pb` end-to-end with `runs/modelCompare/LOOCV_wAIC/` (single model)
  - Test multi-model Δ calculation with synthetic 2-model fixture (small, e.g. 10 sites, 2 chains per model; second model scores worse; verify Δ < 0 for worse model, Δ = 0 for best)
   - Test cross-model site count or ordered site-ID mismatch → hard error
  - Test that LOO-CV best ≠ wAIC best scenario produces correct per-metric Δ values

---

## Task 6: Documentation + external files

- [ ] **Step 6.1:** Modify `docs/superpowers/specs/2026-06-07-phyloai-design.md` — add `modelcompare` to posttree command tree and module structure
- [ ] **Step 6.2:** Update `README.md` — add modelcompare example commands in posttree section
- [ ] **Step 6.3:** Update `README.zh.md` — add modelcompare example commands (Chinese)
- [ ] **Step 6.4:** Create `docs/commands/posttree-modelcompare.md` — English command documentation
- [ ] **Step 6.5:** Create `docs/commands/posttree-modelcompare.zh.md` — Chinese command documentation
- [ ] **Step 6.6:** Modify `skills/phyloai-workflow/SKILL.md` — add modelcompare subcommands usage

---

## Task 7: Cleanup

- [x] **Step 7.1:** ~~Delete old draft files~~ (already deleted during design revision)

---

## Verification

After all tasks:
1. `python -m pytest tests/posttree/test_modelcompare_iqtree.py -v`
2. `python -m pytest tests/posttree/test_modelcompare_pb.py -v`
3. `phyloai posttree modelcompare iqtree --dry-run --matrix runs/modelCompare/EOG090X0A0V.fa --homogeneous-model LG --mrate E,G,R --heterogeneous-model C10`
4. `phyloai posttree modelcompare pb --sitelogl-dir runs/modelCompare/LOOCV_wAIC`
5. `phyloai posttree modelcompare --help` shows iqtree + pb subcommands
6. `phyloai posttree --help` shows modelcompare in subcommand list
7. MCP tool auto-generation: `posttree_modelcompare_iqtree` and `posttree_modelcompare_pb` appear
