# Tool Args and BMGE Projection Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace unsafe `--tool-args` merge semantics with strategy-only `--tool-args`, and make BMGE CODON/AA+NT trim via AA trimming plus NT projection.

**Architecture:** PhyloAI owns batch input/output, work directories, data type, threads, logs, and codon projection. Tool-specific strategy flags are accepted through `--tool-args` only after rejecting managed flags. BMGE paired AA/NT outputs run BMGE in AA mode and infer kept columns for NT projection.

**Tech Stack:** Python, Click, BioPython, pytest, trimAl, BMGE, ClipKIT, MAGUS.

---

### Task 1: Documentation Revisions

**Files:**
- Modify: `docs/superpowers/specs/2026-06-07-phyloai-design.md`
- Modify: `docs/superpowers/specs/2026-06-11-pretree-align-design.md`
- Modify: `docs/superpowers/specs/2026-06-12-pretree-trim-design.md`

- [x] **Step 1:** Replace global `--tool-args` merge semantics with strategy-only `--tool-args` semantics.
- [x] **Step 2:** Document `--tool-args` as a deprecated alias.
- [x] **Step 3:** Document BMGE CODON/AA+NT as AA trim plus NT projection, not direct BMGE `-t CODON` output.

### Task 2: Align Tool Args

**Files:**
- Modify: `phyloai/pretree/align.py`
- Modify: `phyloai/cli/commands/pretree.py`
- Test: `tests/pretree/test_align.py`

- [x] **Step 1:** Add tests that MAGUS `--tool-args` appends allowed strategy flags.
- [x] **Step 2:** Add tests that MAGUS managed flags (`-i`, `-o`, `-d`, `-np`, `--datatype`) are rejected.
- [x] **Step 3:** Implement `tool_args` in `run_align` and CLI, preserving `tool_args` only as a deprecated alias.

### Task 3: Trim Tool Args and ClipKIT Modes

**Files:**
- Modify: `phyloai/pretree/trim.py`
- Modify: `phyloai/cli/commands/pretree.py`
- Test: `tests/pretree/test_trim.py`

- [x] **Step 1:** Add tests that managed tool flags are rejected from `--tool-args`.
- [x] **Step 2:** Replace append/merge `tool_args` handling with strategy-only `tool_args` handling.
- [x] **Step 3:** Update `--clipkit-method` to the 15 modes exposed by `clipkit -h`.

### Task 4: BMGE AA Projection

**Files:**
- Modify: `phyloai/pretree/trim.py`
- Test: `tests/pretree/test_trim.py`

- [x] **Step 1:** Add a failing test proving BMGE CODON translates to AA, runs BMGE `-t AA`, and projects kept columns to NT.
- [x] **Step 2:** Implement kept-column inference from BMGE AA output and projection onto codon MSA.
- [x] **Step 3:** Verify on real `EOG090X0A0V` that BMGE CODON AA length matches BMGE AA length.

### Task 5: Verification

**Files:**
- Test: full repository tests

- [x] **Step 1:** Run targeted trim/align tests.
- [x] **Step 2:** Run representative CLI commands for help, BMGE CODON, and tool-args validation.
- [x] **Step 3:** Run `python -m pytest tests/ -v --tb=short`.
