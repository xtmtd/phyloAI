# pretree align Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build `phyloai pretree align` as a FASTA-only batch MSA command supporting MAFFT (6 strategies) and MAGUS, with parallel execution, optional backtranslation to codon alignments, and full JSON result output.

**Architecture:** Add `phyloai/pretree/align.py` as the library layer (FASTA input scanning, command building, CDS validation, parallel dispatch, result aggregation). Wire a thin CLI layer in `phyloai/cli/commands/pretree.py`. Logging goes to `align.log` in the output directory, following the shared pipeline convention. MAFFT, MAGUS, and trimAl resolve from the environment by default, with optional per-tool path overrides (`--mafft-path`, `--magus-path`, `--trimal-path`).

**Tech Stack:** Python 3.10+, Click, Rich, Biopython SeqIO, concurrent.futures.ProcessPoolExecutor, shlex, tempfile, pytest, Click CliRunner.

---

## File Structure

- Create: `phyloai/pretree/align.py` — library: input scanning, command building, CDS validation, parallel alignment, backtranslation, result aggregation, Rich table renderer
- Modify: `phyloai/cli/commands/pretree.py` — register `align` subcommand
- Create: `tests/pretree/test_align.py` — module-level unit and integration tests
- Create: `tests/cli/test_pretree_align.py` — CLI-level tests via Click CliRunner
- Create: `docs/commands/pretree-align.md` — command documentation

## Required Design Decisions Applied In This Plan

- `pretree align` is FASTA-only. Do not expose `--input-format`; users must run `phyloai pretree convert --to fasta` first for PHYLIP/Nexus/other formats.
- MAFFT, MAGUS, and trimAl are user-provided external tools. Resolve each tool from PATH by default, or from `--mafft-path`, `--magus-path`, or `--trimal-path` when provided.
- Log files are written to the command output directory as `align.log`; do not write to a shared `runs/<run>/logs/` folder. Do not duplicate MAFFT alignment stdout in the log because it is primary output already saved under `seqs/`.
- `--seq-type` accepts `AA`, `NT`, or `auto`; default is `auto`. Auto-detection samples the first few genes and resolves to one molecule type before command construction. `--backtrans` still requires the resolved type to be `AA`.
- `--method magus` is Linux-only in Phase 2 because the pip-distributed MAGUS bundle includes Linux binaries. Non-Linux platforms should fail early with a user-facing error.
- Generated MSA files are validated through shared `core` sequence-output validation helpers before a gene is counted as aligned. Empty output, unparsable FASTA, zero FASTA records, empty sequences, or unequal sequence lengths are skipped with a recorded reason.
- MAGUS `--tool-args` are tokenized with `shlex.split()`. Known internal MAGUS options (`-i`, `-o`, `-d`, `-np`, `--datatype`) are replaced when supplied by the user; all other extra args are appended unchanged.
- Parallel workers must be top-level functions or use top-level callables only. Do not submit nested/local functions to `ProcessPoolExecutor`, because that fails under macOS spawn semantics.
- Backtranslation validation checks both CDS length and AA/CDS taxon identity. A matching taxon count alone is not sufficient.

---

## Task 1: Input Scanning and Command Building

**Files:**
- Create: `phyloai/pretree/align.py`
- Create: `tests/pretree/test_align.py`

- [ ] **Step 1: Write failing tests for input scanning**

Create `tests/pretree/test_align.py`:

```python
from __future__ import annotations

from pathlib import Path

import pytest


def test_scan_input_finds_fasta_files(tmp_path: Path) -> None:
    from phyloai.pretree.align import _scan_input

    (tmp_path / "gene1.fa").write_text(">a\nACGT\n")
    (tmp_path / "gene2.faa").write_text(">b\nMKT\n")
    (tmp_path / "notes.txt").write_text("skip")
    (tmp_path / "empty.fa").write_text("")
    (tmp_path / "subdir").mkdir()

    found, skipped = _scan_input(tmp_path)

    assert len(found) == 2
    assert len(skipped) == 3
    skip_reasons = {s["reason"] for s in skipped}
    assert "empty file" in skip_reasons
    assert "directory" in skip_reasons
    assert "unrecognized extension" in skip_reasons


def test_build_mafft_cmd_linsi(tmp_path: Path) -> None:
    from phyloai.pretree.align import _build_mafft_cmd

    inp = tmp_path / "gene1.fa"
    out = tmp_path / "gene1_aln.fa"
    cmd = _build_mafft_cmd(inp, out, method="linsi")

    assert cmd[0] == "mafft"
    assert "--maxiterate" in cmd
    assert "1000" in cmd
    assert "--localpair" in cmd
    assert "--thread" in cmd
    assert "1" in cmd
    assert str(inp) in cmd


def test_build_mafft_cmd_accepts_explicit_executable(tmp_path: Path) -> None:
    from phyloai.pretree.align import _build_mafft_cmd

    inp = tmp_path / "gene1.fa"
    out = tmp_path / "gene1_aln.fa"
    cmd = _build_mafft_cmd(inp, out, method="linsi", executable="/opt/bin/mafft")

    assert cmd[0] == "/opt/bin/mafft"


def test_build_mafft_cmd_fftns1(tmp_path: Path) -> None:
    from phyloai.pretree.align import _build_mafft_cmd

    inp = tmp_path / "gene1.fa"
    out = tmp_path / "gene1_aln.fa"
    cmd = _build_mafft_cmd(inp, out, method="fftns1")

    assert "--retree" in cmd
    idx = cmd.index("--retree")
    assert cmd[idx + 1] == "1"


def test_build_mafft_cmd_fftns2(tmp_path: Path) -> None:
    from phyloai.pretree.align import _build_mafft_cmd

    inp = tmp_path / "gene1.fa"
    out = tmp_path / "gene1_aln.fa"
    cmd = _build_mafft_cmd(inp, out, method="fftns2")

    assert "--retree" in cmd
    idx = cmd.index("--retree")
    assert cmd[idx + 1] == "2"


def test_build_mafft_cmd_auto(tmp_path: Path) -> None:
    from phyloai.pretree.align import _build_mafft_cmd

    inp = tmp_path / "gene1.fa"
    out = tmp_path / "gene1_aln.fa"
    cmd = _build_mafft_cmd(inp, out, method="auto")

    assert "--auto" in cmd
    assert "--thread" in cmd


def test_build_mafft_cmd_einsi(tmp_path: Path) -> None:
    from phyloai.pretree.align import _build_mafft_cmd

    inp = tmp_path / "gene1.fa"
    out = tmp_path / "out.fa"
    cmd = _build_mafft_cmd(inp, out, method="einsi")

    assert "--genafpair" in cmd


def test_build_mafft_cmd_ginsi(tmp_path: Path) -> None:
    from phyloai.pretree.align import _build_mafft_cmd

    inp = tmp_path / "gene1.fa"
    out = tmp_path / "out.fa"
    cmd = _build_mafft_cmd(inp, out, method="ginsi")

    assert "--globalpair" in cmd


def test_build_magus_cmd_aa(tmp_path: Path) -> None:
    from phyloai.pretree.align import _build_magus_cmd

    inp = tmp_path / "gene1.fa"
    out = tmp_path / "gene1_aln.fa"
    work = tmp_path / "work"
    cmd = _build_magus_cmd(inp, out, work_dir=work, seq_type="AA", tool_args=None)

    assert cmd[0] == "magus"
    assert "-i" in cmd
    assert str(inp) in cmd
    assert "-o" in cmd
    assert str(out) in cmd
    assert "-d" in cmd
    assert str(work) in cmd
    assert "--datatype" in cmd
    idx = cmd.index("--datatype")
    assert cmd[idx + 1] == "protein"


def test_build_magus_cmd_nt(tmp_path: Path) -> None:
    from phyloai.pretree.align import _build_magus_cmd

    inp = tmp_path / "gene1.fa"
    out = tmp_path / "out.fa"
    work = tmp_path / "work"
    cmd = _build_magus_cmd(inp, out, work_dir=work, seq_type="NT", tool_args=None)

    idx = cmd.index("--datatype")
    assert cmd[idx + 1] == "dna"


def test_build_magus_cmd_tool_args(tmp_path: Path) -> None:
    from phyloai.pretree.align import _build_magus_cmd

    inp = tmp_path / "gene1.fa"
    out = tmp_path / "out.fa"
    work = tmp_path / "work"
    cmd = _build_magus_cmd(inp, out, work_dir=work, seq_type="AA", tool_args="--maxsubsetsize 50 --recurse true")

    assert "--maxsubsetsize" in cmd
    assert "50" in cmd
    assert "--recurse" in cmd
    assert "true" in cmd


def test_build_magus_cmd_tool_args_override(tmp_path: Path) -> None:
    from phyloai.pretree.align import _build_magus_cmd

    inp = tmp_path / "gene1.fa"
    out = tmp_path / "out.fa"
    work = tmp_path / "work"
    # --datatype is a known internal option, so tool_args should replace it.
    cmd = _build_magus_cmd(inp, out, work_dir=work, seq_type="AA", tool_args="--datatype dna")

    pairs = list(zip(cmd, cmd[1:]))
    datatype_values = [b for a, b in pairs if a == "--datatype"]
    assert datatype_values == ["dna"]


def test_build_magus_cmd_preserves_unknown_tool_args(tmp_path: Path) -> None:
    from phyloai.pretree.align import _build_magus_cmd

    inp = tmp_path / "gene1.fa"
    out = tmp_path / "out.fa"
    work = tmp_path / "work"
    cmd = _build_magus_cmd(
        inp,
        out,
        work_dir=work,
        seq_type="AA",
        tool_args="--maxsubsetsize 50 --recurse --some-flag=value",
    )

    assert cmd[-4:] == ["--maxsubsetsize", "50", "--recurse", "--some-flag=value"]


def test_build_magus_cmd_accepts_explicit_executable(tmp_path: Path) -> None:
    from phyloai.pretree.align import _build_magus_cmd

    inp = tmp_path / "gene1.fa"
    out = tmp_path / "out.fa"
    work = tmp_path / "work"
    cmd = _build_magus_cmd(inp, out, work_dir=work, seq_type="AA", tool_args=None, executable="/opt/bin/magus")

    assert cmd[0] == "/opt/bin/magus"
```

- [ ] **Step 2: Run tests to verify they fail**

```bash
pytest tests/pretree/test_align.py -v 2>&1 | head -20
```

Expected: FAIL with `ModuleNotFoundError: No module named 'phyloai.pretree.align'`

- [ ] **Step 3: Implement `_scan_input`, `_build_mafft_cmd`, `_build_magus_cmd`**

Create `phyloai/pretree/align.py`:

```python
"""Batch sequence alignment using MAFFT or MAGUS."""

from __future__ import annotations

import shlex
from pathlib import Path
from typing import Any

from phyloai.core.schema import COMMON_ALIGNMENT_EXTENSIONS


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

MAFFT_METHODS = {"fftns1", "fftns2", "auto", "linsi", "einsi", "ginsi"}
MAGUS_METHODS = {"magus"}
ALL_METHODS = MAFFT_METHODS | MAGUS_METHODS

INPUT_EXTENSIONS = {ext for ext in COMMON_ALIGNMENT_EXTENSIONS if ext in
                    {".fa", ".fas", ".fasta", ".faa", ".fna"}}


# ---------------------------------------------------------------------------
# Input scanning
# ---------------------------------------------------------------------------

def _scan_input(
    seq_dir: Path,
) -> tuple[list[Path], list[dict[str, str]]]:
    """Scan seq_dir one level deep; return (found, skipped) lists."""
    found: list[Path] = []
    skipped: list[dict[str, str]] = []

    for entry in sorted(seq_dir.iterdir(), key=lambda p: p.name):
        if entry.is_dir():
            skipped.append({"path": str(entry), "reason": "directory"})
            continue
        if not entry.is_file():
            skipped.append({"path": str(entry), "reason": "not a file"})
            continue
        if entry.stat().st_size == 0:
            skipped.append({"path": str(entry), "reason": "empty file"})
            continue
        if entry.suffix.lower() not in INPUT_EXTENSIONS:
            skipped.append({"path": str(entry), "reason": "unrecognized extension"})
            continue
        found.append(entry)

    return found, skipped


# ---------------------------------------------------------------------------
# Command builders
# ---------------------------------------------------------------------------

def _build_mafft_cmd(
    input_file: Path,
    output_file: Path,
    method: str,
    executable: str = "mafft",
) -> list[str]:
    """Build a mafft command list. Output is redirected by the caller."""
    base = [executable]

    if method == "fftns1":
        base += ["--retree", "1", "--thread", "1"]
    elif method == "fftns2":
        base += ["--retree", "2", "--thread", "1"]
    elif method == "auto":
        base += ["--auto", "--thread", "1"]
    elif method == "linsi":
        base += ["--maxiterate", "1000", "--localpair", "--thread", "1"]
    elif method == "einsi":
        base += ["--maxiterate", "1000", "--genafpair", "--thread", "1"]
    elif method == "ginsi":
        base += ["--maxiterate", "1000", "--globalpair", "--thread", "1"]
    else:
        raise ValueError(f"Unknown MAFFT method: {method!r}")

    base.append(str(input_file))
    return base


def _build_magus_cmd(
    input_file: Path,
    output_file: Path,
    work_dir: Path,
    seq_type: str,
    tool_args: str | None,
    executable: str = "magus",
) -> list[str]:
    """Build a magus command list with controlled tool_args overrides."""
    datatype = "protein" if seq_type == "AA" else "dna"

    internal: dict[str, str] = {
        "-i": str(input_file),
        "-o": str(output_file),
        "-d": str(work_dir),
        "--datatype": datatype,
        "-np": "1",
    }
    known_internal = set(internal)

    # Replace only known internal options. Preserve all other MAGUS args unchanged.
    extra: list[str] = shlex.split(tool_args) if tool_args else []
    passthrough: list[str] = []
    i = 0
    while i < len(extra):
        token = extra[i]
        if token in known_internal:
            if i + 1 >= len(extra):
                raise ValueError(f"MAGUS option {token!r} requires a value in --tool-args")
            internal[token] = extra[i + 1]
            i += 2
        else:
            passthrough.append(token)
            i += 1

    cmd = [executable]
    for key, val in internal.items():
        cmd += [key, val]
    return cmd + passthrough
```

- [ ] **Step 4: Run tests to verify they pass**

```bash
pytest tests/pretree/test_align.py::test_scan_input_finds_fasta_files \
       tests/pretree/test_align.py::test_build_mafft_cmd_linsi \
       tests/pretree/test_align.py::test_build_mafft_cmd_accepts_explicit_executable \
       tests/pretree/test_align.py::test_build_mafft_cmd_fftns1 \
       tests/pretree/test_align.py::test_build_mafft_cmd_fftns2 \
       tests/pretree/test_align.py::test_build_mafft_cmd_auto \
       tests/pretree/test_align.py::test_build_mafft_cmd_einsi \
       tests/pretree/test_align.py::test_build_mafft_cmd_ginsi \
       tests/pretree/test_align.py::test_build_magus_cmd_aa \
       tests/pretree/test_align.py::test_build_magus_cmd_nt \
       tests/pretree/test_align.py::test_build_magus_cmd_tool_args \
       tests/pretree/test_align.py::test_build_magus_cmd_tool_args_override \
       tests/pretree/test_align.py::test_build_magus_cmd_preserves_unknown_tool_args \
       tests/pretree/test_align.py::test_build_magus_cmd_accepts_explicit_executable \
       -v
```

Expected: all PASS

- [ ] **Step 5: Commit**

```bash
git add phyloai/pretree/align.py tests/pretree/test_align.py
git commit -m "feat(align): add input scanning and command builders"
```

---

## Task 2: CDS Validation for Backtranslation

**Files:**
- Modify: `phyloai/pretree/align.py`
- Modify: `tests/pretree/test_align.py`

- [ ] **Step 1: Write failing tests for CDS validation**

Append to `tests/pretree/test_align.py`:

```python
def test_validate_cds_passes_clean_sequences() -> None:
    from phyloai.pretree.align import _validate_cds

    seqs = {"sp1": "ATGGCCTAA", "sp2": "ATGCGCTAG"}  # trailing stops only
    warnings = _validate_cds(seqs, n_aa_taxa=2)
    assert warnings == []


def test_validate_cds_length_not_multiple_of_3() -> None:
    from phyloai.pretree.align import _validate_cds

    seqs = {"sp1": "ATGGC"}  # 5 nt — not multiple of 3
    warnings = _validate_cds(seqs, n_aa_taxa=1)
    assert any("not a multiple of 3" in w for w in warnings)


def test_validate_cds_taxon_count_mismatch() -> None:
    from phyloai.pretree.align import _validate_cds

    seqs = {"sp1": "ATGGCCTAA"}  # 1 taxon
    warnings = _validate_cds(seqs, n_aa_taxa=3)  # AA alignment has 3
    assert any("taxon count mismatch" in w for w in warnings)


def test_validate_cds_taxon_id_mismatch() -> None:
    from phyloai.pretree.align import _validate_cds

    seqs = {"sp1": "ATGGCCTAA", "wrong_sp": "ATGCGCTAG"}
    warnings = _validate_cds(seqs, aa_taxa={"sp1", "sp2"})
    assert any("taxon ID mismatch" in w for w in warnings)


def test_validate_cds_internal_stop_codon() -> None:
    from phyloai.pretree.align import _validate_cds

    # TAA at positions 3-5 (internal), not at end
    seqs = {"sp1": "ATGTAAGCT"}
    warnings = _validate_cds(seqs, n_aa_taxa=1)
    assert any("internal stop codon" in w for w in warnings)


def test_validate_cds_trailing_stop_not_flagged() -> None:
    from phyloai.pretree.align import _validate_cds

    seqs = {"sp1": "ATGGCCTAA"}  # TAA is last codon
    warnings = _validate_cds(seqs, n_aa_taxa=1)
    assert not any("internal stop codon" in w for w in warnings)
```

- [ ] **Step 2: Run tests to verify they fail**

```bash
pytest tests/pretree/test_align.py::test_validate_cds_passes_clean_sequences \
       tests/pretree/test_align.py::test_validate_cds_length_not_multiple_of_3 \
       tests/pretree/test_align.py::test_validate_cds_taxon_count_mismatch \
       tests/pretree/test_align.py::test_validate_cds_internal_stop_codon \
       tests/pretree/test_align.py::test_validate_cds_trailing_stop_not_flagged \
       -v
```

Expected: FAIL with `ImportError: cannot import name '_validate_cds'`

- [ ] **Step 3: Implement `_validate_cds`**

Append to `phyloai/pretree/align.py`:

```python
# ---------------------------------------------------------------------------
# CDS validation
# ---------------------------------------------------------------------------

STOP_CODONS = {"TAA", "TAG", "TGA"}


def _validate_cds(
    sequences: dict[str, str],
    n_aa_taxa: int | None = None,
    aa_taxa: set[str] | None = None,
) -> list[str]:
    """
    Lightweight pre-checks on CDS sequences before backtranslation.

    Returns a list of warning strings. Empty list means all checks pass.
    Trailing stop codons are allowed (trimAl handles them with -ignorestopcodon).
    Internal stop codons are flagged.
    """
    warnings: list[str] = []

    if n_aa_taxa is None and aa_taxa is not None:
        n_aa_taxa = len(aa_taxa)

    if n_aa_taxa is not None and len(sequences) != n_aa_taxa:
        warnings.append(
            f"taxon count mismatch: CDS file has {len(sequences)} sequences "
            f"but AA alignment has {n_aa_taxa}"
        )

    if aa_taxa is not None:
        cds_taxa = set(sequences)
        if cds_taxa != aa_taxa:
            missing = sorted(aa_taxa - cds_taxa)
            extra = sorted(cds_taxa - aa_taxa)
            warnings.append(
                "taxon ID mismatch: "
                f"missing in CDS={missing}; extra in CDS={extra}"
            )

    for name, seq in sequences.items():
        seq_upper = seq.upper().replace("-", "")

        # Length check
        if len(seq_upper) % 3 != 0:
            warnings.append(
                f"{name}: CDS length {len(seq_upper)} is not a multiple of 3"
            )
            continue  # skip codon scan if length is wrong

        # Internal stop codon check (exclude last codon)
        codons = [seq_upper[i:i+3] for i in range(0, len(seq_upper) - 3, 3)]
        for pos, codon in enumerate(codons):
            if codon in STOP_CODONS:
                warnings.append(
                    f"{name}: internal stop codon '{codon}' at codon position {pos + 1}"
                )
                break  # one warning per sequence is enough

    return warnings
```

- [ ] **Step 4: Run tests to verify they pass**

```bash
pytest tests/pretree/test_align.py::test_validate_cds_passes_clean_sequences \
       tests/pretree/test_align.py::test_validate_cds_length_not_multiple_of_3 \
       tests/pretree/test_align.py::test_validate_cds_taxon_count_mismatch \
       tests/pretree/test_align.py::test_validate_cds_taxon_id_mismatch \
       tests/pretree/test_align.py::test_validate_cds_internal_stop_codon \
       tests/pretree/test_align.py::test_validate_cds_trailing_stop_not_flagged \
       -v
```

Expected: all PASS

- [ ] **Step 5: Commit**

```bash
git add phyloai/pretree/align.py tests/pretree/test_align.py
git commit -m "feat(align): add CDS pre-validation for backtranslation"
```

---

## Task 3: Single-Gene Alignment and Backtranslation

**Files:**
- Modify: `phyloai/pretree/align.py`
- Modify: `tests/pretree/test_align.py`

This task requires `mafft` and `trimal` available on PATH. Tests use `pytest.importorskip` pattern with `shutil.which` guards.

- [ ] **Step 1: Write failing tests for `_align_one` and `_backtrans_one`**

Append to `tests/pretree/test_align.py`:

```python
import shutil


def test_align_one_mafft_linsi_produces_output(tmp_path: Path) -> None:
    if not shutil.which("mafft"):
        pytest.skip("mafft not found")
    from phyloai.pretree.align import _align_one

    inp = tmp_path / "gene1.fa"
    inp.write_text(">sp1\nMKTLLLTLVVVTIVC\n>sp2\nMKTLLLTLAAVTIVC\n>sp3\nMKTLLLTLVVVTIVC\n")
    out_dir = tmp_path / "out"
    out_dir.mkdir()

    result = _align_one(inp, out_dir, method="linsi", seq_type="AA",
                        tool_args=None, dry_run=False)

    assert result["status"] == "success"
    assert Path(result["output_aa"]).exists()
    assert result["n_taxa"] == 3
    assert result["alignment_length"] > 0
    assert result["wall_time"] > 0


def test_align_one_dry_run_creates_no_files(tmp_path: Path) -> None:
    from phyloai.pretree.align import _align_one

    inp = tmp_path / "gene1.fa"
    inp.write_text(">sp1\nMKT\n>sp2\nMKT\n")
    out_dir = tmp_path / "out"
    out_dir.mkdir()

    result = _align_one(inp, out_dir, method="linsi", seq_type="AA",
                        tool_args=None, dry_run=True)

    assert result["status"] == "dry_run"
    assert result["cmd"] is not None
    assert not any(out_dir.iterdir())


def test_align_one_failed_tool_returns_skipped(tmp_path: Path) -> None:
    from phyloai.pretree.align import _align_one

    # Empty file will cause mafft to fail
    inp = tmp_path / "bad.fa"
    inp.write_text("")
    out_dir = tmp_path / "out"
    out_dir.mkdir()

    result = _align_one(inp, out_dir, method="linsi", seq_type="AA",
                        tool_args=None, dry_run=False)

    assert result["status"] == "skipped"
    assert "reason" in result


def test_backtrans_one_produces_nt_alignment(tmp_path: Path) -> None:
    if not shutil.which("trimal"):
        pytest.skip("trimal not found")
    from phyloai.pretree.align import _backtrans_one

    aa_aln = tmp_path / "gene1_aa.fa"
    aa_aln.write_text(
        ">sp1\nMK-\n"
        ">sp2\nMKT\n"
    )
    nt_file = tmp_path / "gene1.fa"
    nt_file.write_text(
        ">sp1\nATGAAA\n"
        ">sp2\nATGAAAACT\n"
    )
    out_nt = tmp_path / "gene1_nt.fa"

    result = _backtrans_one(aa_aln, nt_file, out_nt, dry_run=False)

    assert result["status"] in {"success", "skipped"}  # trimal may reject mismatched lengths
    if result["status"] == "success":
        assert out_nt.exists()
```

- [ ] **Step 2: Run tests to verify they fail**

```bash
pytest tests/pretree/test_align.py::test_align_one_mafft_linsi_produces_output \
       tests/pretree/test_align.py::test_align_one_dry_run_creates_no_files \
       tests/pretree/test_align.py::test_align_one_failed_tool_returns_skipped \
       tests/pretree/test_align.py::test_backtrans_one_produces_nt_alignment \
       -v
```

Expected: FAIL with `ImportError: cannot import name '_align_one'`

- [ ] **Step 3: Implement `_align_one` and `_backtrans_one`**

Append to `phyloai/pretree/align.py`:

```python
import subprocess
import tempfile
import time
from collections.abc import Callable

from Bio import SeqIO


# ---------------------------------------------------------------------------
# Single-gene alignment
# ---------------------------------------------------------------------------

def _align_one(
    gene_path: Path,
    output_dir: Path,
    method: str,
    seq_type: str,
    tool_args: str | None,
    dry_run: bool,
    mafft_executable: str = "mafft",
    magus_executable: str = "magus",
) -> dict[str, Any]:
    """Align one gene file. Returns a result dict with status/output_aa/etc."""
    out_aa = output_dir / f"{gene_path.stem}.fa"

    if method in MAFFT_METHODS:
        cmd = _build_mafft_cmd(gene_path, out_aa, method, executable=mafft_executable)
    else:
        work_dir = Path(tempfile.mkdtemp(prefix="phyloai_magus_"))
        cmd = _build_magus_cmd(gene_path, out_aa, work_dir, seq_type, tool_args, executable=magus_executable)

    if dry_run:
        return {"status": "dry_run", "input": str(gene_path), "cmd": cmd}

    start = time.monotonic()
    try:
        if method in MAFFT_METHODS:
            # mafft writes to stdout
            proc = subprocess.run(
                cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
            )
            wall_time = time.monotonic() - start
            if proc.returncode != 0:
                return {
                    "status": "skipped",
                    "input": str(gene_path),
                    "reason": f"mafft exited with code {proc.returncode}: {proc.stderr[:200]}",
                    "tool_stdout": proc.stdout,
                    "tool_stderr": proc.stderr,
                    "wall_time": wall_time,
                }
            out_aa.parent.mkdir(parents=True, exist_ok=True)
            out_aa.write_text(proc.stdout)
            stdout, stderr = proc.stdout, proc.stderr
        else:
            # magus writes directly to -o
            out_aa.parent.mkdir(parents=True, exist_ok=True)
            proc = subprocess.run(
                cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
            )
            wall_time = time.monotonic() - start
            _cleanup_magus_workdir(cmd)
            if proc.returncode != 0:
                return {
                    "status": "skipped",
                    "input": str(gene_path),
                    "reason": f"magus exited with code {proc.returncode}: {proc.stderr[:200]}",
                    "tool_stdout": proc.stdout,
                    "tool_stderr": proc.stderr,
                    "wall_time": wall_time,
                }
            stdout, stderr = proc.stdout, proc.stderr

    except Exception as exc:
        return {
            "status": "skipped",
            "input": str(gene_path),
            "reason": str(exc),
            "tool_stdout": "",
            "tool_stderr": "",
            "wall_time": time.monotonic() - start,
        }

    # Count taxa and alignment length from output
    try:
        records = list(SeqIO.parse(str(out_aa), "fasta"))
        n_taxa = len(records)
        alignment_length = len(records[0].seq) if records else 0
    except Exception:
        n_taxa = 0
        alignment_length = 0

    return {
        "status": "success",
        "input": str(gene_path),
        "output_aa": str(out_aa),
        "output_nt": None,
        "n_taxa": n_taxa,
        "alignment_length": alignment_length,
        "wall_time": wall_time,
        "tool_cmd": " ".join(cmd),
        "tool_stdout": stdout,
        "tool_stderr": stderr,
        "warnings": [],
    }


def _align_one_worker(args: tuple[Path, Path, str, str, str | None, bool, str, str]) -> dict[str, Any]:
    """Top-level worker for ProcessPoolExecutor; required for macOS spawn."""
    gene_path, output_dir, method, seq_type, tool_args, dry_run, mafft_exe, magus_exe = args
    return _align_one(
        gene_path,
        output_dir,
        method=method,
        seq_type=seq_type,
        tool_args=tool_args,
        dry_run=dry_run,
        mafft_executable=mafft_exe,
        magus_executable=magus_exe,
    )


def _cleanup_magus_workdir(cmd: list[str]) -> None:
    """Remove MAGUS working directory from a cmd list."""
    import shutil as _shutil
    try:
        idx = cmd.index("-d")
        work = Path(cmd[idx + 1])
        if work.exists():
            _shutil.rmtree(work, ignore_errors=True)
    except (ValueError, IndexError):
        pass


# ---------------------------------------------------------------------------
# Backtranslation
# ---------------------------------------------------------------------------

def _backtrans_one(
    aa_aln_path: Path,
    nt_path: Path,
    output_nt_path: Path,
    dry_run: bool,
    executable: str = "trimal",
) -> dict[str, Any]:
    """
    Run trimAl -backtrans to produce a codon alignment.
    Always passes -ignorestopcodon.
    """
    cmd = [
        executable,
        "-in", str(aa_aln_path),
        "-backtrans", str(nt_path),
        "-ignorestopcodon",
        "-out", str(output_nt_path),
        "-fasta",
    ]

    if dry_run:
        return {"status": "dry_run", "cmd": cmd}

    output_nt_path.parent.mkdir(parents=True, exist_ok=True)
    start = time.monotonic()
    try:
        proc = subprocess.run(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
        wall_time = time.monotonic() - start
    except Exception as exc:
        return {
            "status": "skipped",
            "reason": f"trimal error: {exc}",
            "tool_stdout": "",
            "tool_stderr": "",
            "wall_time": time.monotonic() - start,
        }

    if proc.returncode != 0:
        return {
            "status": "skipped",
            "reason": f"trimal -backtrans exited with code {proc.returncode}: {proc.stderr[:300]}",
            "tool_stdout": proc.stdout,
            "tool_stderr": proc.stderr,
            "wall_time": wall_time,
        }

    return {
        "status": "success",
        "output_nt": str(output_nt_path),
        "tool_cmd": " ".join(cmd),
        "tool_stdout": proc.stdout,
        "tool_stderr": proc.stderr,
        "wall_time": wall_time,
    }
```

- [ ] **Step 4: Run tests to verify they pass**

```bash
pytest tests/pretree/test_align.py::test_align_one_dry_run_creates_no_files \
       tests/pretree/test_align.py::test_align_one_failed_tool_returns_skipped \
       -v
# Run tool-dependent tests only if mafft/trimal available:
pytest tests/pretree/test_align.py::test_align_one_mafft_linsi_produces_output \
       tests/pretree/test_align.py::test_backtrans_one_produces_nt_alignment \
       -v
```

Expected: all available tests PASS (tool-dependent tests skip if tool absent)

- [ ] **Step 5: Commit**

```bash
git add phyloai/pretree/align.py tests/pretree/test_align.py
git commit -m "feat(align): implement _align_one and _backtrans_one"
```

---

## Task 4: `run_align` — Parallel Orchestration and Result Aggregation

**Files:**
- Modify: `phyloai/pretree/align.py`
- Modify: `tests/pretree/test_align.py`

- [ ] **Step 1: Write failing tests for `run_align`**

Append to `tests/pretree/test_align.py`:

```python
def test_run_align_aa_only_dry_run(tmp_path: Path) -> None:
    from phyloai.pretree.align import run_align

    seq_dir = tmp_path / "seqs"
    seq_dir.mkdir()
    (seq_dir / "gene1.fa").write_text(">a\nMKT\n>b\nMKA\n")
    (seq_dir / "gene2.fa").write_text(">a\nGHT\n>b\nGHA\n")
    out_dir = tmp_path / "out"

    payload = run_align(
        seq_dir=seq_dir,
        output_dir=out_dir,
        method="linsi",
        seq_type="AA",
        dry_run=True,
    )

    assert payload["status"] == "success"
    assert payload["data"]["summary"]["n_input_files"] == 2
    assert not out_dir.exists()  # dry_run: no files created


def test_run_align_output_dir_conflict(tmp_path: Path) -> None:
    from phyloai.pretree.align import run_align

    seq_dir = tmp_path / "seqs"
    seq_dir.mkdir()
    (seq_dir / "gene1.fa").write_text(">a\nMKT\n")
    out_dir = tmp_path / "out"
    out_dir.mkdir()
    (out_dir / "old.txt").write_text("old")

    import pytest as _pytest
    with _pytest.raises(ValueError, match="already exists and is non-empty"):
        run_align(seq_dir=seq_dir, output_dir=out_dir, method="linsi",
                  seq_type="AA", overwrite=False)


def test_run_align_overwrite_clears_directory(tmp_path: Path) -> None:
    if not shutil.which("mafft"):
        pytest.skip("mafft not found")
    from phyloai.pretree.align import run_align

    seq_dir = tmp_path / "seqs"
    seq_dir.mkdir()
    (seq_dir / "gene1.fa").write_text(">sp1\nMKTLL\n>sp2\nMKTAA\n")
    out_dir = tmp_path / "out"
    out_dir.mkdir()
    (out_dir / "old.txt").write_text("old")

    payload = run_align(seq_dir=seq_dir, output_dir=out_dir, method="linsi",
                        seq_type="AA", overwrite=True)

    assert not (out_dir / "old.txt").exists()
    assert payload["status"] == "success"


def test_run_align_all_skipped_returns_error(tmp_path: Path) -> None:
    from phyloai.pretree.align import run_align

    seq_dir = tmp_path / "seqs"
    seq_dir.mkdir()
    # Empty file will be skipped during scanning
    (seq_dir / "bad.fa").write_text("")
    out_dir = tmp_path / "out"

    import pytest as _pytest
    with _pytest.raises(ValueError, match="No genes were aligned"):
        run_align(seq_dir=seq_dir, output_dir=out_dir, method="linsi",
                  seq_type="AA", overwrite=False)


def test_run_align_backtrans_requires_nt_dir(tmp_path: Path) -> None:
    from phyloai.pretree.align import run_align

    seq_dir = tmp_path / "seqs"
    seq_dir.mkdir()
    (seq_dir / "gene1.fa").write_text(">a\nMKT\n")
    out_dir = tmp_path / "out"

    import pytest as _pytest
    with _pytest.raises(ValueError, match="--nt-dir"):
        run_align(seq_dir=seq_dir, output_dir=out_dir, method="linsi",
                  seq_type="AA", backtrans=True, nt_dir=None)


def test_run_align_nt_seq_type_with_backtrans_raises(tmp_path: Path) -> None:
    from phyloai.pretree.align import run_align

    seq_dir = tmp_path / "seqs"
    seq_dir.mkdir()
    nt_dir = tmp_path / "nt"
    nt_dir.mkdir()
    out_dir = tmp_path / "out"

    import pytest as _pytest
    with _pytest.raises(ValueError, match="--backtrans requires"):
        run_align(seq_dir=seq_dir, output_dir=out_dir, method="linsi",
                  seq_type="NT", backtrans=True, nt_dir=nt_dir)


def test_run_align_tool_args_ignored_for_mafft_emits_warning(tmp_path: Path) -> None:
    from phyloai.pretree.align import run_align

    seq_dir = tmp_path / "seqs"
    seq_dir.mkdir()
    (seq_dir / "gene1.fa").write_text(">a\nMKT\n>b\nMKA\n")
    out_dir = tmp_path / "out"

    payload = run_align(
        seq_dir=seq_dir,
        output_dir=out_dir,
        method="linsi",
        seq_type="AA",
        tool_args="--maxsubsetsize 50",
        dry_run=True,
    )

    assert any("ignored" in w.lower() for w in payload["data"].get("warnings", []))


def test_run_align_key_results_shape(tmp_path: Path) -> None:
    if not shutil.which("mafft"):
        pytest.skip("mafft not found")
    from phyloai.pretree.align import run_align

    seq_dir = tmp_path / "seqs"
    seq_dir.mkdir()
    (seq_dir / "gene1.fa").write_text(">sp1\nMKTLL\n>sp2\nMKTAA\n>sp3\nMKTVV\n")
    (seq_dir / "gene2.fa").write_text(">sp1\nGHTLL\n>sp2\nGHTAA\n>sp3\nGHTVV\n")
    out_dir = tmp_path / "out"

    payload = run_align(seq_dir=seq_dir, output_dir=out_dir, method="linsi",
                        seq_type="AA")

    kr = payload["key_results"]
    assert "n_aligned" in kr
    assert "n_skipped" in kr
    assert "mean_alignment_length" in kr
    assert "mean_n_taxa" in kr
    assert "method" in kr
    assert "backtrans" in kr
```

- [ ] **Step 2: Run tests to verify they fail**

```bash
pytest tests/pretree/test_align.py::test_run_align_aa_only_dry_run \
       tests/pretree/test_align.py::test_run_align_output_dir_conflict \
       tests/pretree/test_align.py::test_run_align_all_skipped_returns_error \
       tests/pretree/test_align.py::test_run_align_backtrans_requires_nt_dir \
       tests/pretree/test_align.py::test_run_align_nt_seq_type_with_backtrans_raises \
       tests/pretree/test_align.py::test_run_align_tool_args_ignored_for_mafft_emits_warning \
       -v
```

Expected: FAIL with `ImportError: cannot import name 'run_align'`

- [ ] **Step 3: Implement `run_align` and `render_align_summary_table`**

Append to `phyloai/pretree/align.py` (add missing imports at top of file: `import shutil`, `import json`, `import datetime`, `from concurrent.futures import ProcessPoolExecutor, as_completed`):

```python
# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def run_align(
    seq_dir: Path,
    output_dir: Path,
    method: str,
    seq_type: str,
    backtrans: bool = False,
    nt_dir: Path | None = None,
    threads: int = 4,
    tool_args: str | None = None,
    mafft_path: Path | None = None,
    magus_path: Path | None = None,
    trimal_path: Path | None = None,
    overwrite: bool = False,
    dry_run: bool = False,
    progress_callback: Callable[[Path], None] | None = None,
) -> dict[str, Any]:
    """
    Batch-align all sequences in seq_dir. Returns JSON-serialisable payload.
    Raises ValueError on parameter errors or total failure.
    """
    import shutil as _shutil
    from concurrent.futures import ProcessPoolExecutor, as_completed

    run_start = time.monotonic()

    # --- Parameter validation ---
    if backtrans and nt_dir is None:
        raise ValueError("--nt-dir is required when --backtrans is set.")
    if backtrans and seq_type == "NT":
        raise ValueError("--backtrans requires --seq-type AA (backtrans produces NT from AA alignment).")

    global_warnings: list[str] = []
    if tool_args and method in MAFFT_METHODS:
        global_warnings.append(
            f"--tool-args is ignored for MAFFT method '{method}'; "
            "it is only used with --method magus."
        )
        tool_args = None

    mafft_exe, magus_exe, trimal_exe = _resolve_tool_paths(
        method=method,
        backtrans=backtrans,
        mafft_path=mafft_path,
        magus_path=magus_path,
        trimal_path=trimal_path,
        dry_run=dry_run,
    )

    # --- Output directory setup ---
    if not dry_run:
        if output_dir.exists() and any(output_dir.iterdir()):
            if not overwrite:
                raise ValueError(
                    f"Output directory '{output_dir}' already exists and is non-empty. "
                    "Use --overwrite to replace it."
                )
            _shutil.rmtree(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

    # --- Input scanning ---
    found, scan_skipped = _scan_input(seq_dir)
    skipped: list[dict[str, str]] = list(scan_skipped)

    if not found and not dry_run:
        raise ValueError("No genes were aligned: no valid input files found.")

    # --- Determine output subdirectory for aligned files ---
    if backtrans:
        aa_out_dir = output_dir / "seqs" / "faa"
        nt_out_dir = output_dir / "seqs" / "fna"
    else:
        aa_out_dir = output_dir / "seqs"
        nt_out_dir = None

    if not dry_run:
        aa_out_dir.mkdir(parents=True, exist_ok=True)
        if nt_out_dir:
            nt_out_dir.mkdir(parents=True, exist_ok=True)

    # --- Parallel alignment ---
    file_results: list[dict[str, Any]] = []

    worker_args = [
        (g, aa_out_dir, method, seq_type, tool_args, dry_run, mafft_exe, magus_exe)
        for g in found
    ]
    all_tool_results: list[dict[str, Any]] = []
    with ProcessPoolExecutor(max_workers=threads) as pool:
        futures = {pool.submit(_align_one_worker, arg): arg[0] for arg in worker_args}
        for future in as_completed(futures):
            gene_path = futures[future]
            try:
                result = future.result()
            except Exception as exc:
                result = {
                    "status": "skipped",
                    "input": str(gene_path),
                    "reason": str(exc),
                }
            if result["status"] == "skipped":
                skipped.append({"path": result["input"], "reason": result.get("reason", "unknown")})
            else:
                file_results.append(result)
            all_tool_results.append(result)
            if progress_callback:
                progress_callback(gene_path)

    # --- Backtranslation pass ---
    n_backtrans = 0
    if backtrans and nt_dir and not dry_run:
        for res in file_results:
            gene_stem = Path(res["input"]).stem
            # Find matching CDS file by stem
            nt_candidates = [
                p for p in nt_dir.iterdir()
                if p.is_file() and p.stem == gene_stem
            ]
            if not nt_candidates:
                res["warnings"].append(
                    f"no matching CDS file found in --nt-dir for gene '{gene_stem}'"
                )
                continue

            nt_path = nt_candidates[0]
            aa_aln_path = Path(res["output_aa"])

            # Pre-check CDS sequences
            try:
                aa_records = list(SeqIO.parse(str(aa_aln_path), "fasta"))
                aa_taxa = {r.id for r in aa_records}
                nt_records = list(SeqIO.parse(str(nt_path), "fasta"))
                nt_seqs = {r.id: str(r.seq) for r in nt_records}
            except Exception as exc:
                res["warnings"].append(f"could not read CDS file: {exc}")
                continue

            cds_warnings = _validate_cds(nt_seqs, n_aa_taxa=res["n_taxa"], aa_taxa=aa_taxa)
            if cds_warnings:
                res["warnings"].extend(cds_warnings)
                res["warnings"].append("backtrans skipped due to CDS validation errors")
                continue

            out_nt = nt_out_dir / f"{gene_stem}.fa"
            bt_result = _backtrans_one(aa_aln_path, nt_path, out_nt, dry_run=False, executable=trimal_exe)
            all_tool_results.append({**bt_result, "input": str(nt_path)})
            if bt_result["status"] == "success":
                res["output_nt"] = bt_result["output_nt"]
                n_backtrans += 1
            else:
                res["warnings"].append(bt_result["reason"])

    # --- Write log ---
    if not dry_run and all_tool_results:
        _write_align_log(output_dir, all_tool_results)

    # --- Fail if nothing aligned ---
    if not dry_run and not file_results:
        raise ValueError("No genes were aligned: all input files failed or were skipped.")

    # --- Aggregate results ---
    n_aligned = len(file_results)
    aligned_lengths = [r["alignment_length"] for r in file_results if r.get("alignment_length")]
    aligned_taxa = [r["n_taxa"] for r in file_results if r.get("n_taxa")]
    mean_len = round(sum(aligned_lengths) / len(aligned_lengths), 1) if aligned_lengths else 0.0
    mean_taxa = round(sum(aligned_taxa) / len(aligned_taxa), 1) if aligned_taxa else 0.0

    all_warnings = list(global_warnings)
    for r in file_results:
        all_warnings.extend(r.get("warnings", []))

    payload: dict[str, Any] = {
        "status": "success",
        "command": (
            f"phyloai pretree align --seq-dir {seq_dir} --method {method} "
            f"--seq-type {seq_type} --threads {threads}"
        ),
        "wall_time": time.monotonic() - run_start,
        "tool_versions": {},
        "params": {
            "seq_dir": str(seq_dir),
            "method": method,
            "seq_type": seq_type,
            "backtrans": backtrans,
            "nt_dir": str(nt_dir) if nt_dir else None,
            "output_dir": str(output_dir),
            "threads": threads,
            "tool_args": tool_args,
            "mafft_path": str(mafft_path) if mafft_path else None,
            "magus_path": str(magus_path) if magus_path else None,
            "trimal_path": str(trimal_path) if trimal_path else None,
            "overwrite": overwrite,
        },
        "key_results": {
            "n_aligned": n_aligned,
            "n_skipped": len(skipped),
            "method": method,
            "backtrans": backtrans,
            "mean_alignment_length": mean_len,
            "mean_n_taxa": mean_taxa,
        },
        "error": None,
        "data": {
            "summary": {
                "n_input_files": len(found) + len(scan_skipped),
                "n_aligned": n_aligned,
                "n_backtrans": n_backtrans,
                "n_skipped": len(skipped),
            },
            "files": [
                {
                    "input": r["input"],
                    "output_aa": r.get("output_aa"),
                    "output_nt": r.get("output_nt"),
                    "n_taxa": r.get("n_taxa", 0),
                    "alignment_length": r.get("alignment_length", 0),
                    "wall_time": r.get("wall_time", 0.0),
                    "warnings": r.get("warnings", []),
                }
                for r in file_results
            ],
            "skipped": skipped,
            "warnings": all_warnings,
        },
    }
    return payload


# ---------------------------------------------------------------------------
# Log writer
# ---------------------------------------------------------------------------

def _write_align_log(output_dir: Path, file_results: list[dict[str, Any]]) -> None:
    """Append per-gene tool log entries to align.log in output_dir."""
    import datetime as _dt
    log_path = output_dir / "align.log"
    with open(log_path, "a") as fh:
        for res in file_results:
            ts = _dt.datetime.now().isoformat(timespec="seconds")
            fh.write(
                f"{'='*60}\n"
                f"timestamp:  {ts}\n"
                f"input:      {res.get('input')}\n"
                f"cmd:        {res.get('tool_cmd', '')}\n"
                f"returncode: {'0' if res['status'] == 'success' else 'non-zero'}\n"
                f"wall_time:  {res.get('wall_time', 0.0):.2f}s\n"
                f"--- stdout ---\n{res.get('tool_stdout', '')}\n"
                f"--- stderr ---\n{res.get('tool_stderr', '')}\n"
            )


# ---------------------------------------------------------------------------
# Rich table renderer
# ---------------------------------------------------------------------------

def render_align_summary_table(summary: dict[str, Any]) -> "Table":
    from rich.table import Table as _Table
    table = _Table(title="pretree align summary")
    table.add_column("Metric")
    table.add_column("Value")
    for key in ["n_input_files", "n_aligned", "n_backtrans", "n_skipped"]:
        table.add_row(key, str(summary.get(key, "")))
    return table
```

- [ ] **Step 4: Add missing imports at the top of `align.py`**

The top of `phyloai/pretree/align.py` should now read:

```python
"""Batch sequence alignment using MAFFT or MAGUS."""

from __future__ import annotations

import datetime
import shlex
import shutil
import subprocess
import tempfile
import time
from collections.abc import Callable
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Any

from Bio import SeqIO

from phyloai.core.schema import COMMON_ALIGNMENT_EXTENSIONS
```

- [ ] **Step 5: Run tests to verify they pass**

```bash
pytest tests/pretree/test_align.py::test_run_align_aa_only_dry_run \
       tests/pretree/test_align.py::test_run_align_output_dir_conflict \
       tests/pretree/test_align.py::test_run_align_all_skipped_returns_error \
       tests/pretree/test_align.py::test_run_align_backtrans_requires_nt_dir \
       tests/pretree/test_align.py::test_run_align_nt_seq_type_with_backtrans_raises \
       tests/pretree/test_align.py::test_run_align_tool_args_ignored_for_mafft_emits_warning \
       -v
```

Expected: all PASS

- [ ] **Step 6: Commit**

```bash
git add phyloai/pretree/align.py tests/pretree/test_align.py
git commit -m "feat(align): implement run_align parallel orchestration and result aggregation"
```

---

## Task 5: Tool Path Resolution and Version Detection

**Files:**
- Modify: `phyloai/pretree/align.py`
- Modify: `tests/pretree/test_align.py`

- [ ] **Step 1: Write failing test**

Append to `tests/pretree/test_align.py`:

```python
def test_resolve_tool_paths_accepts_explicit_mafft_path(tmp_path: Path) -> None:
    from phyloai.pretree.align import _resolve_tool_paths

    fake = tmp_path / "mafft"
    fake.write_text("#!/bin/sh\n")

    mafft_exe, magus_exe, trimal_exe = _resolve_tool_paths(
        method="linsi",
        backtrans=False,
        mafft_path=fake,
        magus_path=None,
        trimal_path=None,
        dry_run=True,
    )

    assert mafft_exe == str(fake)
    assert magus_exe == "magus"
    assert trimal_exe == "trimal"


def test_detect_tool_versions_mafft() -> None:
    if not shutil.which("mafft"):
        pytest.skip("mafft not found")
    from phyloai.pretree.align import _detect_tool_versions

    versions = _detect_tool_versions(method="linsi", backtrans=False, mafft_path=None, magus_path=None, trimal_path=None)
    assert "mafft" in versions
    assert versions["mafft"]  # non-empty string


def test_detect_tool_versions_backtrans_includes_trimal() -> None:
    from phyloai.pretree.align import _detect_tool_versions

    versions = _detect_tool_versions(method="linsi", backtrans=True, mafft_path=None, magus_path=None, trimal_path=None)
    assert "trimal" in versions
```

- [ ] **Step 2: Run tests to verify they fail**

```bash
pytest tests/pretree/test_align.py::test_resolve_tool_paths_accepts_explicit_mafft_path \
       tests/pretree/test_align.py::test_detect_tool_versions_mafft \
       tests/pretree/test_align.py::test_detect_tool_versions_backtrans_includes_trimal \
       -v
```

Expected: FAIL with `ImportError: cannot import name '_detect_tool_versions'`

- [ ] **Step 3: Implement `_resolve_tool_paths` and `_detect_tool_versions`**

Append to `phyloai/pretree/align.py`:

```python
def _resolve_tool_paths(
    method: str,
    backtrans: bool,
    mafft_path: Path | None,
    magus_path: Path | None,
    trimal_path: Path | None,
    dry_run: bool,
) -> tuple[str, str, str]:
    """Resolve executable paths while keeping dry-run usable on systems without tools."""
    from phyloai.core.env import ToolEnv

    env = ToolEnv()
    if method == "magus":
        magus_exe = str(_validate_executable_path(magus_path, "magus")) if magus_path else ("magus" if dry_run else str(env.require("magus")))
        mafft_exe = str(mafft_path) if mafft_path else "mafft"
    else:
        mafft_exe = str(_validate_executable_path(mafft_path, "mafft")) if mafft_path else ("mafft" if dry_run else str(env.require("mafft")))
        magus_exe = str(magus_path) if magus_path else "magus"

    if backtrans:
        trimal_exe = str(_validate_executable_path(trimal_path, "trimal")) if trimal_path else ("trimal" if dry_run else str(env.require("trimal")))
    else:
        trimal_exe = str(trimal_path) if trimal_path else "trimal"

    return mafft_exe, magus_exe, trimal_exe


def _validate_executable_path(path: Path, tool_name: str) -> Path:
    if not path.exists() or path.is_dir():
        raise FileNotFoundError(f"Required tool '{tool_name}' not found at explicit path: {path}")
    return path


def _detect_tool_versions(
    method: str,
    backtrans: bool,
    mafft_path: Path | None,
    magus_path: Path | None,
    trimal_path: Path | None,
) -> dict[str, str]:
    """Return version strings for tools that will be invoked."""
    from phyloai.core.env import TOOL_REGISTRY, ToolEnv, ToolStatus

    tool_paths = {}
    if mafft_path:
        tool_paths["mafft"] = mafft_path
    if magus_path:
        tool_paths["magus"] = magus_path
    if trimal_path:
        tool_paths["trimal"] = trimal_path

    env = ToolEnv(tool_paths=tool_paths)
    versions: dict[str, str] = {}
    names = ["magus" if method == "magus" else "mafft"]
    if backtrans:
        names.append("trimal")

    for name in names:
        meta = TOOL_REGISTRY[name]
        # align treats MAFFT, MAGUS, and trimAl as user-provided tools only:
        # explicit path override first, otherwise PATH lookup. Do not resolve bundled tools here.
        info = env._detect_tool(
            name,
            version_flag=meta.get("version_flag", ""),
            version_args=meta.get("version_args"),
            bundled=False,
            path_aliases=meta.get("path_aliases"),
        )
        if info.status == ToolStatus.OK and info.version:
            versions[name] = info.version
    return versions
```

Also update `run_align` to populate `tool_versions` in the payload. Find the `"tool_versions": {}` line and replace with:

```python
"tool_versions": _detect_tool_versions(
    method=method,
    backtrans=backtrans,
    mafft_path=mafft_path,
    magus_path=magus_path,
    trimal_path=trimal_path,
),
```

- [ ] **Step 4: Run tests to verify they pass**

```bash
pytest tests/pretree/test_align.py::test_resolve_tool_paths_accepts_explicit_mafft_path \
       tests/pretree/test_align.py::test_detect_tool_versions_mafft \
       tests/pretree/test_align.py::test_detect_tool_versions_backtrans_includes_trimal \
       -v
```

Expected: PASS (or SKIP if tools absent)

- [ ] **Step 5: Commit**

```bash
git add phyloai/pretree/align.py tests/pretree/test_align.py
git commit -m "feat(align): add tool version detection to result payload"
```

---

## Task 6: Explicit Tool Path Validation

**Files:**
- Modify: `phyloai/pretree/align.py`
- Modify: `tests/pretree/test_align.py`

- [ ] **Step 1: Write failing tests**

Append to `tests/pretree/test_align.py`:

```python
def test_run_align_missing_explicit_mafft_path_raises(tmp_path: Path) -> None:
    from phyloai.pretree.align import run_align

    seq_dir = tmp_path / "seqs"
    seq_dir.mkdir()
    (seq_dir / "gene1.fa").write_text(">a\nMKT\n>b\nMKA\n")
    out_dir = tmp_path / "out"

    import pytest as _pytest
    with _pytest.raises(FileNotFoundError, match="mafft"):
        run_align(
            seq_dir=seq_dir,
            output_dir=out_dir,
            method="linsi",
            seq_type="AA",
            mafft_path=tmp_path / "missing-mafft",
            dry_run=True,
        )


def test_run_align_missing_explicit_trimal_path_raises(tmp_path: Path) -> None:
    from phyloai.pretree.align import run_align

    seq_dir = tmp_path / "seqs"
    seq_dir.mkdir()
    nt_dir = tmp_path / "nt"
    nt_dir.mkdir()
    (seq_dir / "gene1.fa").write_text(">a\nMKT\n>b\nMKA\n")
    out_dir = tmp_path / "out"

    import pytest as _pytest
    with _pytest.raises(FileNotFoundError, match="trimal"):
        run_align(
            seq_dir=seq_dir,
            output_dir=out_dir,
            method="linsi",
            seq_type="AA",
            backtrans=True,
            nt_dir=nt_dir,
            trimal_path=tmp_path / "missing-trimal",
            dry_run=True,
        )
```

- [ ] **Step 2: Run test to verify it fails**

```bash
pytest tests/pretree/test_align.py::test_run_align_missing_explicit_mafft_path_raises \
       tests/pretree/test_align.py::test_run_align_missing_explicit_trimal_path_raises \
       -v
```

Expected: FAIL

- [ ] **Step 3: Confirm path validation is wired through `run_align`**

`run_align` should already call `_resolve_tool_paths(...)` from Task 5. Confirm that call happens before any output directory is created and passes all three explicit path values:

```python
    mafft_exe, magus_exe, trimal_exe = _resolve_tool_paths(
        method=method,
        backtrans=backtrans,
        mafft_path=mafft_path,
        magus_path=magus_path,
        trimal_path=trimal_path,
        dry_run=dry_run,
    )
```

Do not add a separate `shutil.which()`-based `_check_environment()` helper; use one path-resolution path through `ToolEnv` so explicit paths and PATH lookup behave consistently.

- [ ] **Step 4: Run test to verify it passes**

```bash
pytest tests/pretree/test_align.py::test_run_align_missing_explicit_mafft_path_raises \
       tests/pretree/test_align.py::test_run_align_missing_explicit_trimal_path_raises \
       -v
```

Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add phyloai/pretree/align.py tests/pretree/test_align.py
git commit -m "feat(align): validate explicit tool paths"
```

---

## Task 7: CLI Subcommand

**Files:**
- Modify: `phyloai/cli/commands/pretree.py`
- Create: `tests/cli/test_pretree_align.py`

- [ ] **Step 1: Write failing CLI tests**

Create `tests/cli/test_pretree_align.py`:

```python
from __future__ import annotations

import json
import shutil
from pathlib import Path

import pytest
from click.testing import CliRunner

from phyloai.cli.main import cli


def _make_seq_dir(tmp_path: Path) -> Path:
    d = tmp_path / "seqs"
    d.mkdir()
    (d / "gene1.fa").write_text(">sp1\nMKTLL\n>sp2\nMKTAA\n>sp3\nMKTVV\n")
    (d / "gene2.fa").write_text(">sp1\nGHTLL\n>sp2\nGHTAA\n>sp3\nGHTVV\n")
    return d


def test_align_help() -> None:
    runner = CliRunner()
    result = runner.invoke(cli, ["pretree", "align", "--help"])
    assert result.exit_code == 0
    assert "--method" in result.output
    assert "--seq-type" in result.output
    assert "--backtrans" in result.output
    assert "--nt-dir" in result.output
    assert "--threads" in result.output
    assert "--mafft-path" in result.output
    assert "--magus-path" in result.output
    assert "--trimal-path" in result.output
    assert "--input-format" not in result.output


def test_align_dry_run_exits_zero(tmp_path: Path) -> None:
    runner = CliRunner()
    seq_dir = _make_seq_dir(tmp_path)
    out_dir = tmp_path / "out"

    result = runner.invoke(cli, [
        "pretree", "align",
        "--seq-dir", str(seq_dir),
        "--output-dir", str(out_dir),
        "--method", "linsi",
        "--seq-type", "AA",
        "--dry-run",
    ])

    assert result.exit_code == 0
    assert not out_dir.exists()


def test_align_backtrans_without_nt_dir_exits_1(tmp_path: Path) -> None:
    runner = CliRunner()
    seq_dir = _make_seq_dir(tmp_path)
    out_dir = tmp_path / "out"

    result = runner.invoke(cli, [
        "pretree", "align",
        "--seq-dir", str(seq_dir),
        "--output-dir", str(out_dir),
        "--method", "linsi",
        "--seq-type", "AA",
        "--backtrans",
    ])

    assert result.exit_code == 1


def test_align_nt_seq_type_with_backtrans_exits_1(tmp_path: Path) -> None:
    runner = CliRunner()
    seq_dir = _make_seq_dir(tmp_path)
    nt_dir = tmp_path / "nt"
    nt_dir.mkdir()
    out_dir = tmp_path / "out"

    result = runner.invoke(cli, [
        "pretree", "align",
        "--seq-dir", str(seq_dir),
        "--output-dir", str(out_dir),
        "--method", "linsi",
        "--seq-type", "NT",
        "--backtrans",
        "--nt-dir", str(nt_dir),
    ])

    assert result.exit_code == 1


def test_align_dry_run_ignores_output_dir_conflict(tmp_path: Path) -> None:
    runner = CliRunner()
    seq_dir = _make_seq_dir(tmp_path)
    out_dir = tmp_path / "out"
    out_dir.mkdir()
    (out_dir / "old.txt").write_text("old")

    result = runner.invoke(cli, [
        "pretree", "align",
        "--seq-dir", str(seq_dir),
        "--output-dir", str(out_dir),
        "--method", "linsi",
        "--dry-run",
    ])

    assert result.exit_code == 0


def test_align_writes_result_json(tmp_path: Path) -> None:
    if not shutil.which("mafft"):
        pytest.skip("mafft not found")
    runner = CliRunner()
    seq_dir = _make_seq_dir(tmp_path)
    out_dir = tmp_path / "out"

    result = runner.invoke(cli, [
        "pretree", "align",
        "--seq-dir", str(seq_dir),
        "--output-dir", str(out_dir),
        "--method", "linsi",
        "--seq-type", "AA",
        "--threads", "2",
    ])

    assert result.exit_code == 0
    result_json = out_dir / "result.json"
    assert result_json.exists()
    payload = json.loads(result_json.read_text())
    assert payload["status"] == "success"
    assert payload["key_results"]["n_aligned"] == 2
    assert (out_dir / "seqs" / "gene1.fa").exists()
    assert (out_dir / "seqs" / "gene2.fa").exists()
    assert (out_dir / "align.log").exists()


def test_align_quiet_suppresses_rich_output(tmp_path: Path) -> None:
    if not shutil.which("mafft"):
        pytest.skip("mafft not found")
    runner = CliRunner()
    seq_dir = _make_seq_dir(tmp_path)
    out_dir = tmp_path / "out"

    result = runner.invoke(cli, [
        "pretree", "align",
        "--seq-dir", str(seq_dir),
        "--output-dir", str(out_dir),
        "--method", "linsi",
        "--quiet",
    ])

    assert result.exit_code == 0
    # No Rich tables in output
    assert "pretree align summary" not in result.output
```

- [ ] **Step 2: Run tests to verify they fail**

```bash
pytest tests/cli/test_pretree_align.py -v 2>&1 | head -20
```

Expected: FAIL — `align` subcommand not yet registered

- [ ] **Step 3: Register `align` subcommand in `pretree.py`**

Add import at the top of `phyloai/cli/commands/pretree.py`:

```python
from phyloai.pretree.align import render_align_summary_table, run_align
```

Add the subcommand after the `convert` command block:

```python
@pretree.command(
    "align",
    help=(
        "Align a directory of unaligned sequences using MAFFT or MAGUS.\n\n"
        "Strategies: fftns1/fftns2 are fastest but least accurate; "
        "linsi/einsi/ginsi offer high accuracy; "
        "magus is slowest but best for large or difficult datasets (Linux only).\n\n"
        "Use --backtrans with --nt-dir to also produce codon-level NT alignments "
        "from a protein alignment using trimAl backtranslation.\n\n"
        "--threads controls how many genes are aligned in parallel; each "
        "individual alignment uses a single thread."
    ),
)
@click.option("--seq-dir", type=click.Path(file_okay=False, path_type=Path),
              required=True, help="Input directory of unaligned sequence files.")
@click.option("--method",
              type=click.Choice(["fftns1", "fftns2", "auto", "linsi", "einsi", "ginsi", "magus"]),
              default="linsi", show_default=True,
              help=(
                  "Alignment strategy. "
                  "fftns1/fftns2: fast, lower accuracy. "
                  "auto: MAFFT chooses strategy automatically. "
                  "linsi/einsi/ginsi: high accuracy. "
                  "magus: highest accuracy, slowest, best for large datasets (Linux only)."
               ))
@click.option("--seq-type", type=click.Choice(["AA", "NT", "auto"]), default="auto",
              show_default=True, help="Molecule type of input sequences. Auto-detects from first gene if 'auto'.")
@click.option("--backtrans", is_flag=True, default=False,
              help="Produce codon NT alignments via trimAl -backtrans. Requires --nt-dir.")
@click.option("--nt-dir", type=click.Path(file_okay=False, path_type=Path), default=None,
              help="Directory of unaligned CDS sequences for --backtrans mode.")
@click.option("--output-dir", "-o", type=click.Path(file_okay=False, path_type=Path),
              default=Path("runs/pretree/align"), show_default=True,
              help="Output directory; contains seqs/, align.log, result.json.")
@click.option("--threads", "-t", type=int, default=4, show_default=True,
              help="Number of genes to align in parallel (each uses 1 thread).")
@click.option("--tool-args", type=str, default=None,
              help="Extra arguments passed to magus only; ignored with a warning for MAFFT methods.")
@click.option("--mafft-path", type=click.Path(dir_okay=False, path_type=Path), default=None,
              help="Explicit MAFFT executable path for MAFFT methods; PATH lookup is used when omitted.")
@click.option("--magus-path", type=click.Path(dir_okay=False, path_type=Path), default=None,
              help="Explicit MAGUS executable path for --method magus; PATH lookup is used when omitted.")
@click.option("--trimal-path", type=click.Path(dir_okay=False, path_type=Path), default=None,
              help="Explicit trimAl executable path for --backtrans; PATH lookup is used when omitted.")
@click.option("--overwrite", is_flag=True, default=False,
              help="Delete and recreate a non-empty output directory before running.")
@click.option("--dry-run", is_flag=True, default=False,
              help="Print commands without executing; creates no files.")
@click.option("--quiet", "-q", is_flag=True, default=False,
              help="Suppress Rich terminal output except errors.")
def align_command(
    seq_dir: Path,
    method: str,
    seq_type: str,
    backtrans: bool,
    nt_dir: Path | None,
    output_dir: Path,
    threads: int,
    tool_args: str | None,
    mafft_path: Path | None,
    magus_path: Path | None,
    trimal_path: Path | None,
    overwrite: bool,
    dry_run: bool,
    quiet: bool,
) -> None:
    if threads < 1:
        _fail("--threads must be at least 1.", 1)
    if not seq_dir.exists():
        _fail(f"--seq-dir '{seq_dir}' does not exist.", 1)
    if nt_dir is not None and not nt_dir.exists():
        _fail(f"--nt-dir '{nt_dir}' does not exist.", 1)

    payload: dict | None = None
    error_msg: str | None = None

    def _invoke(progress_callback=None):
        return run_align(
            seq_dir=seq_dir,
            output_dir=output_dir,
            method=method,
            seq_type=seq_type,
            backtrans=backtrans,
            nt_dir=nt_dir,
            threads=threads,
            tool_args=tool_args,
            mafft_path=mafft_path,
            magus_path=magus_path,
            trimal_path=trimal_path,
            overwrite=overwrite,
            dry_run=dry_run,
            progress_callback=progress_callback,
        )

    if not quiet and not dry_run:
        from rich.progress import Progress
        with Progress(console=console, transient=True) as progress:
            # Scan to get total count for progress bar
            from phyloai.pretree.align import _scan_input
            found, _ = _scan_input(seq_dir)
            task = progress.add_task("Aligning sequences", total=len(found))
            try:
                payload = _invoke(progress_callback=lambda _: progress.advance(task))
            except (ValueError, FileNotFoundError) as exc:
                error_msg = str(exc)
    else:
        try:
            payload = _invoke()
        except (ValueError, FileNotFoundError) as exc:
            error_msg = str(exc)

    if error_msg is not None:
        exit_code = 3 if "not found" in error_msg.lower() else 1
        _fail(error_msg, exit_code)

    if dry_run:
        click.echo(f"Dry run: {payload['data']['summary']['n_input_files']} genes would be aligned.")
        return

    result_path = output_dir / "result.json"
    result_path.parent.mkdir(parents=True, exist_ok=True)
    import json as _json
    with open(result_path, "w") as fh:
        _json.dump(payload, fh, indent=2)

    if not quiet:
        console.print(render_align_summary_table(payload["data"]["summary"]))
        click.echo(
            f"Alignments saved to {output_dir / 'seqs'}", err=True
        )
        click.echo(f"Results saved to {result_path}", err=True)
        for w in payload["data"].get("warnings", []):
            click.echo(f"Warning: {w}", err=True)
```

- [ ] **Step 4: Run CLI tests**

```bash
pytest tests/cli/test_pretree_align.py -v
```

Expected: all tests PASS (tool-dependent tests skip if mafft absent)

- [ ] **Step 5: Commit**

```bash
git add phyloai/cli/commands/pretree.py tests/cli/test_pretree_align.py
git commit -m "feat(align): register pretree align CLI subcommand"
```

---

## Task 8: Full Test Suite and Command Documentation

**Files:**
- Create: `docs/commands/pretree-align.md`
- Modify: `README.md`

- [ ] **Step 1: Run full test suite**

```bash
pytest tests/ -v --tb=short 2>&1 | tail -30
```

Expected: all existing tests still PASS; new align tests PASS or SKIP

- [ ] **Step 2: Create `docs/commands/pretree-align.md`**

Create the file with this content:

```markdown
# phyloai pretree align

## Purpose

Align a directory of unaligned sequence files using MAFFT or MAGUS. Produces one aligned FASTA file per input gene. Optionally produces codon-level NT alignments via trimAl backtranslation.

This command does not perform format conversion. Inputs must be FASTA. Run `phyloai pretree convert --to fasta` first for PHYLIP, Nexus, or other formats.

## Usage

Minimal:
```bash
phyloai pretree align --seq-dir ./raw_aa
```

Full:
```bash
phyloai pretree align \
  --seq-dir ./raw_aa \
  --method linsi \
  --seq-type AA \
  --output-dir ./runs/pretree/align \
  --threads 4
```

With backtranslation:
```bash
phyloai pretree align \
  --seq-dir ./raw_aa \
  --method linsi \
  --seq-type AA \
  --backtrans \
  --nt-dir ./raw_nt \
  --output-dir ./runs/pretree/align \
  --threads 4
```

## Parameters

| Parameter | Default | Notes |
|-----------|---------|-------|
| `--seq-dir` | required | Directory of unaligned sequence files |
| `--method` | `linsi` | fftns1, fftns2, auto, linsi, einsi, ginsi, magus |
| `--seq-type` | `auto` | AA, NT, or auto (auto-detects from first few genes) |
| `--backtrans` | off | Produce NT codon alignment; requires --nt-dir |
| `--nt-dir` | — | Unaligned CDS directory for backtrans |
| `--output-dir` / `-o` | `runs/pretree/align` | Output directory |
| `--threads` / `-t` | 4 | Concurrent alignment tasks (each uses 1 thread) |
| `--tool-args` | — | Extra args for MAGUS only; ignored for MAFFT methods |
| `--mafft-path` | — | Explicit MAFFT executable path for MAFFT methods |
| `--magus-path` | — | Explicit MAGUS executable path for `--method magus` |
| `--trimal-path` | — | Explicit trimAl executable path for `--backtrans` |
| `--overwrite` | off | Delete and recreate non-empty output directory |
| `--dry-run` | off | Print commands, create no files |
| `--quiet` / `-q` | off | Suppress terminal output except errors |

## Inputs

Scans `--seq-dir` one level deep for files with extensions: `.fa`, `.fas`, `.fasta`, `.faa`, `.fna`. Subdirectories, empty files, and unrecognized extensions are skipped.

## Outputs

**Mode AA or NT only:**
```
runs/pretree/align/
├── seqs/
│   ├── gene1.fa
│   └── ...
├── align.log
└── result.json
```

**Mode AA + backtrans:**
```
runs/pretree/align/
├── seqs/
│   ├── faa/
│   │   └── gene1.fa
│   └── fna/
│       └── gene1.fa
├── align.log
└── result.json
```

`result.json` contains `key_results` with `n_aligned`, `method`, `mean_alignment_length`, and `mean_n_taxa` for report integration.

## Examples

```bash
# Fast alignment for large dataset
phyloai pretree align --seq-dir ./raw_aa --method fftns2 --threads 8

# High-accuracy protein alignment + codon NT alignment
phyloai pretree align --seq-dir ./raw_aa --seq-type AA \
  --backtrans --nt-dir ./raw_nt --method linsi --threads 4

# NT direct alignment
phyloai pretree align --seq-dir ./raw_nt --seq-type NT --method linsi

# MAGUS with extra options
phyloai pretree align --seq-dir ./raw_aa --method magus \
  --tool-args "--maxsubsetsize 200" --threads 4

# Preview commands without running
phyloai pretree align --seq-dir ./raw_aa --method linsi --dry-run
```

## Warnings and Errors

| Condition | Behaviour |
|-----------|-----------|
| `--backtrans` without `--nt-dir` | Exit 1 |
| `--seq-type NT` with `--backtrans` | Exit 1 |
| `--seq-type auto` detects NT with `--backtrans` | Exit 1 |
| `mafft` or `magus` not found | Exit 3 |
| `trimal` not found with `--backtrans` | Exit 3 |
| `--method magus` on non-Linux | Exit 1 (MAGUS bundled binaries are Linux-only) |
| Non-empty output directory | Exit 1 (use `--overwrite`) |
| CDS length not multiple of 3 | Backtrans skipped for that gene, warning in result.json |
| Internal stop codon in CDS | Backtrans skipped for that gene, warning in result.json |
| trimAl exits non-zero | Backtrans skipped for that gene, stderr captured as warning |
| Generated MSA is empty, unparsable, or has unequal sequence lengths | Gene skipped with reason in result.json |
| All genes fail | Exit 1 |
| `--tool-args` used with MAFFT method | Warning printed, args ignored |

## Notes

- Downstream: pass `--msa-dir` of `phyloai pretree trim` to `seqs/` (Mode 1/2) or `seqs/faa/` (Mode 3 AA) as appropriate.
- `result.json` `key_results` feeds the Methods paragraph: "X genes were aligned using MAFFT L-INS-i; mean alignment length Y aa."
- Run `phyloai doctor` to verify MAFFT, MAGUS, and trimAl are detected.
```

- [ ] **Step 3: Update README command index**

Open `README.md` and add `pretree align` to the command index table under the pretree section:

```markdown
| `phyloai pretree align` | Align sequences with MAFFT or MAGUS | [docs](docs/commands/pretree-align.md) |
```

- [ ] **Step 4: Run full test suite one final time**

```bash
pytest tests/ -v --tb=short 2>&1 | tail -40
```

Expected: all tests PASS or SKIP (no FAIL)

- [ ] **Step 5: Commit**

```bash
git add docs/commands/pretree-align.md README.md
git commit -m "docs(align): add pretree-align command documentation and README entry"
```

- [ ] **Step 6: Push to GitHub**

```bash
git push
```
