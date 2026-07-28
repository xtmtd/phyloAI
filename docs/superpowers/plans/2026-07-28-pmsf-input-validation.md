# PMSF Input Validation Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Reject invalid PMSF inputs with actionable errors and document the actual AliSim and readpb formats.

**Architecture:** Keep validation in the existing `bi_readpb.py` conversion and partition-writing functions. Convert file-read failures to `ValueError` at the PMSF boundary so they become normal command errors; preserve `-p` as the documented edge-proportional AliSim option.

**Tech Stack:** Python, pytest, Click, IQ-TREE/AliSim documentation.

## Global Constraints

- Do not add dependencies or abstractions.
- Keep `.meansiterates` as headerless, zero-based `site rate` input.
- Use `-p` for edge-proportional AliSim partition simulation; mention `-q` only as the edge-equal alternative.

---

### Task 1: Add PMSF Validation Tests

**Files:**
- Modify: `tests/tree/test_bi_readpb.py`

**Interfaces:**
- Consumes: `_convert_siteprofiles_to_sitefreq(path)` and `_write_pmsf_partition(...)`.
- Produces: regression coverage for malformed rows, missing chain metadata, and invalid IDs/rates.

- [ ] **Step 1: Write failing tests**

```python
with pytest.raises(ValueError, match="Invalid siteprofiles row 3"):
    _convert_siteprofiles_to_sitefreq(siteprofiles)

with pytest.raises(ValueError, match="Invalid sitefreq site ID"):
    _write_pmsf_partition(...)

with pytest.raises(ValueError, match="Invalid meansiterates site index"):
    _write_pmsf_partition(...)

with pytest.raises(ValueError, match="Invalid meansiterates rate"):
    _write_pmsf_partition(...)
```

- [ ] **Step 2: Run the focused tests**

Run: `pytest tests/tree/test_bi_readpb.py -q`
Expected: failures for inputs currently accepted or misclassified.

### Task 2: Validate Inputs at Their Boundaries

**Files:**
- Modify: `phyloai/tree/bi_readpb.py:102-243`
- Test: `tests/tree/test_bi_readpb.py`

**Interfaces:**
- Consumes: readpb `.siteprofiles`, `.meansiterates`, `.trace`, and `.log` files.
- Produces: `ValueError` with an input-specific message before writing `partition.PMSF.nex`.

- [ ] **Step 1: Make the minimal implementation**

```python
if len(parts) != 21:
    raise ValueError(f"Invalid siteprofiles row {line_number}: expected site ID plus 20 frequencies")

if site <= 0:
    raise ValueError(f"Invalid sitefreq site ID on row {line_number}: {site}")
if site < 0:
    raise ValueError(f"Invalid meansiterates site index on row {line_number}: {site}")
if not isfinite(numeric_rate) or numeric_rate <= 0:
    raise ValueError(f"Invalid meansiterates rate on row {line_number}")
```

Catch `OSError` from trace/log reads inside `_write_partition_if_ready` and return a `pmsf_partition` validation error, rather than allowing the CLI to classify it as a missing executable.

- [ ] **Step 2: Run the focused tests**

Run: `pytest tests/tree/test_bi_readpb.py -q`
Expected: PASS.

### Task 3: Complete Error Metadata and Documentation

**Files:**
- Modify: `phyloai/tree/bi_readpb.py:549-562`
- Modify: `docs/commands/tree-bi.md`
- Modify: `docs/commands/tree-bi.zh.md`
- Modify: `docs/superpowers/plans/2026-07-26-phyloai-tree-bi-subcommands.md`

**Interfaces:**
- Consumes: PMSF error return payload and command documentation.
- Produces: reproducible error metadata and correct IQ-TREE/readpb instructions.

- [ ] **Step 1: Add omitted parameters to the PMSF error payload**

```python
"overwrite": overwrite,
"burnin": burnin,
```

- [ ] **Step 2: Correct documentation**

Document `-p partition.PMSF.nex` as AliSim edge-proportional partitions, and state that `-q` is the edge-equal alternative. Replace any planned `Site Rate C_Rate` header with headerless, zero-based `site rate` rows.

- [ ] **Step 3: Run regression tests**

Run: `pytest tests/tree/test_bi_readpb.py tests/cli/test_tree.py -q`
Expected: PASS.
