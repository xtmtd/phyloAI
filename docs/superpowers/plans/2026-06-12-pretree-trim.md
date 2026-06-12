# pretree trim Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Implement `phyloai pretree trim` — batch MSA trimming with trimAl, BMGE, and ClipKIT supporting AA/NT/CODON/AA+NT modes with parallel execution and checkpoint/resume.

**Architecture:** Core library in `phyloai/pretree/trim.py` (tool builders, workers, orchestration); CLI registration appended to `phyloai/cli/commands/pretree.py`; CODON validation helper added to `phyloai/core/sequence_normalization.py`. Follows `pretree align` patterns throughout.

**Tech Stack:** Python 3.10+, BioPython, Click 8+, Rich, ProcessPoolExecutor, trimAl (bundled), BMGE (bundled jar), ClipKIT (pip)

---

## File Map

| Action | File | Responsibility |
|--------|------|----------------|
| Create | `phyloai/pretree/trim.py` | All trim logic: tool command builders, per-gene workers, `run_trim` orchestrator |
| Modify | `phyloai/core/sequence_normalization.py` | Add `validate_codon_msa()` helper |
| Modify | `phyloai/cli/commands/pretree.py` | Register `trim` subcommand; add to `_PretreeGroup` order |
| Create | `tests/pretree/test_trim.py` | Unit + integration tests for trim module |
| Create | `tests/core/test_sequence_normalization_codon.py` | Unit tests for `validate_codon_msa` |

---

## Task 1: `validate_codon_msa` in `core/sequence_normalization.py`

**Files:**
- Modify: `phyloai/core/sequence_normalization.py`
- Create: `tests/core/test_sequence_normalization_codon.py`

CODON MSA validation: checks alignment length divisible by 3, detects internal stop codons (warn + continue), strips terminal stop codons.

- [ ] **Step 1.1: Write failing tests**

Create `tests/core/test_sequence_normalization_codon.py`:

```python
from __future__ import annotations
import pytest
from phyloai.core.sequence_normalization import validate_codon_msa


def test_valid_codon_msa_no_issues():
    seqs = {
        "seq1": "ATGGCTTCT",  # M-A-S, 9 nt, 3 codons
        "seq2": "ATGGCTTCT",
    }
    result = validate_codon_msa(seqs)
    assert result.skip is False
    assert result.warnings == []
    # sequences unchanged
    assert result.sequences["seq1"] == "ATGGCTTCT"


def test_alignment_length_not_divisible_by_3_triggers_skip():
    # MSA alignment length (all seqs same length) not divisible by 3
    seqs = {
        "seq1": "ATGGCT",     # 6 nt alignment columns — ok
        "seq2": "ATGGCT",
    }
    result = validate_codon_msa(seqs)
    assert result.skip is False  # 6 % 3 == 0

    seqs_bad = {
        "seq1": "ATGGCTA",    # 7 nt alignment columns — not divisible by 3
        "seq2": "ATGGCTA",
    }
    result_bad = validate_codon_msa(seqs_bad)
    assert result_bad.skip is True
    assert any("multiple of 3" in w for w in result_bad.warnings)


def test_internal_stop_codon_warns_but_continues():
    # TAA is stop; placed before last codon → internal
    seqs = {
        "seq1": "ATGTAAGCT",  # ATG-TAA-GCT: TAA at position 1 (before last)
        "seq2": "ATGGCTTCT",
    }
    result = validate_codon_msa(seqs)
    assert result.skip is False
    assert any("internal stop" in w for w in result.warnings)


def test_terminal_stop_codon_is_removed():
    seqs = {
        "seq1": "ATGGCTTAA",  # ATG-GCT-TAA: TAA is terminal → stripped
        "seq2": "ATGGCTTAA",
    }
    result = validate_codon_msa(seqs)
    assert result.skip is False
    assert result.sequences["seq1"] == "ATGGCT"
    assert result.sequences["seq2"] == "ATGGCT"


def test_terminal_stop_with_trailing_gaps_preserved():
    # MSA has trailing gaps after stop codon; strip stop but preserve gap structure
    # seq1: ATG GCT --- TAA → after stripping terminal stop: ATG GCT ---
    # seq2: ATG GCT TCT ---  (no stop)
    seqs = {
        "seq1": "ATGGCT---TAA",  # 12 cols; ungapped = ATGGCTTAA; terminal TAA → strip last codon from ungapped, map back
        "seq2": "ATGGCTTCT---",
    }
    result = validate_codon_msa(seqs)
    assert result.skip is False
    # seq1 ungapped after strip: ATGGCT (6 nt); alignment cols = 12 - 3 = 9
    assert len(result.sequences["seq1"].replace("-", "")) == 6
    # seq2 unchanged
    assert result.sequences["seq2"] == "ATGGCTTCT---"


def test_gap_columns_in_msa_handled():
    # Gaps count toward alignment length; strip gaps for per-seq divisibility check
    seqs = {
        "seq1": "ATG---GCT",  # 9 nt alignment; ungapped = ATGGCT = 6 = divisible by 3
        "seq2": "ATGGCTTCT",
    }
    result = validate_codon_msa(seqs)
    assert result.skip is False


def test_all_gapped_sequence_does_not_raise():
    seqs = {
        "seq1": "---------",
        "seq2": "ATGGCTTCT",
    }
    result = validate_codon_msa(seqs)
    assert result.skip is False  # ungapped length 0 is divisible by 3


def test_empty_seqs_dict_returns_skip():
    result = validate_codon_msa({})
    assert result.skip is True
```

- [ ] **Step 1.2: Run tests to confirm they fail**

```bash
cd /Users/zf/data/coding/phyloAI
python -m pytest tests/core/test_sequence_normalization_codon.py -v 2>&1 | head -30
```

Expected: ImportError or AttributeError — `validate_codon_msa` does not exist yet.

- [ ] **Step 1.3: Implement `validate_codon_msa`**

Append to `phyloai/core/sequence_normalization.py`:

```python
from dataclasses import dataclass, field as _field


STOP_CODONS = frozenset({"TAA", "TAG", "TGA"})


@dataclass
class CodonMsaValidation:
    skip: bool                          # True = hard error, caller must skip gene
    sequences: dict[str, str]           # possibly modified (terminal stop stripped)
    warnings: list[str] = _field(default_factory=list)


def validate_codon_msa(sequences: dict[str, str]) -> CodonMsaValidation:
    """Validate a codon-aligned MSA (gap-aware).

    Rules:
    - Empty input → skip (hard error).
    - Alignment length (number of columns, shared by all seqs) not divisible by 3 → skip.
    - Terminal stop codon (in ungapped seq) → remove last codon from each seq that has one,
      preserving MSA column structure by removing the 3 rightmost non-gap characters.
    - Internal stop codon → warn and continue (lenient, mirrors trimAl -ignorestopcodon).
    - Gap-only sequences (ungapped length 0) → no codon-level validation needed.
    """
    if not sequences:
        return CodonMsaValidation(skip=True, sequences={}, warnings=["empty sequence dict"])

    # Check alignment-level column count (all seqs must be same length in a valid MSA)
    lengths = {len(seq) for seq in sequences.values()}
    if len(lengths) > 1:
        return CodonMsaValidation(
            skip=True, sequences=sequences,
            warnings=[f"MSA sequences have unequal lengths: {sorted(lengths)}"]
        )
    aln_len = next(iter(lengths))
    if aln_len % 3 != 0:
        return CodonMsaValidation(
            skip=True, sequences=sequences,
            warnings=[f"alignment length {aln_len} is not a multiple of 3 (codon_length_not_multiple_of_3)"]
        )

    warnings: list[str] = []
    result_seqs: dict[str, str] = {}

    for name, seq in sequences.items():
        ungapped = seq.replace("-", "").upper()

        if len(ungapped) == 0:
            result_seqs[name] = seq
            continue

        codons = [ungapped[i:i+3] for i in range(0, len(ungapped), 3)]

        # Check and strip terminal stop codon
        if codons and codons[-1] in STOP_CODONS:
            # Remove the last 3 non-gap characters from the raw (gapped) sequence
            # to preserve alignment column structure for other sequences in the MSA.
            chars = list(seq)
            removed = 0
            for i in range(len(chars) - 1, -1, -1):
                if chars[i] != "-":
                    chars[i] = "\x00"   # mark for deletion
                    removed += 1
                    if removed == 3:
                        break
            seq = "".join(c for c in chars if c != "\x00")
            ungapped = seq.replace("-", "").upper()
            codons = [ungapped[i:i+3] for i in range(0, len(ungapped), 3)] if ungapped else []

        # Check internal stops (all but the last codon)
        for pos, codon in enumerate(codons[:-1]):
            if codon in STOP_CODONS:
                warnings.append(
                    f"{name}: internal stop codon '{codon}' at codon position {pos + 1} "
                    "(proceeding with lenient handling)"
                )
                break

        result_seqs[name] = seq

    return CodonMsaValidation(skip=False, sequences=result_seqs, warnings=warnings)
```

- [ ] **Step 1.4: Run tests to confirm they pass**

```bash
python -m pytest tests/core/test_sequence_normalization_codon.py -v
```

Expected: all 7 tests PASS.

---

## Task 2: `trim.py` — constants, file scanning, command builders

**Files:**
- Create: `phyloai/pretree/trim.py`
- Create: `tests/pretree/test_trim.py`

- [ ] **Step 2.1: Write failing tests for constants and file scanning**

Create `tests/pretree/test_trim.py`:

```python
from __future__ import annotations
from pathlib import Path
import pytest


def test_scan_input_finds_fasta_files(tmp_path: Path) -> None:
    from phyloai.pretree.trim import _scan_input

    (tmp_path / "gene1.fa").write_text(">a\nMKT\n")
    (tmp_path / "gene2.faa").write_text(">b\nMKT\n")
    (tmp_path / "notes.txt").write_text("skip")
    (tmp_path / "empty.fa").write_text("")
    (tmp_path / "subdir").mkdir()

    found, skipped = _scan_input(tmp_path)
    assert len(found) == 2
    reasons = {s["reason"] for s in skipped}
    assert "empty file" in reasons
    assert "directory" in reasons
    assert "unrecognized extension" in reasons


def test_build_trimal_cmd_aa_only(tmp_path: Path) -> None:
    from phyloai.pretree.trim import _build_trimal_cmd

    inp = tmp_path / "gene1.fa"
    out = tmp_path / "gene1_trim.fa"
    cmd = _build_trimal_cmd(inp, out, method="automated1", executable="trimal")

    assert cmd[0] == "trimal"
    assert "-in" in cmd
    assert str(inp) in cmd
    assert "-out" in cmd
    assert str(out) in cmd
    assert "-automated1" in cmd
    assert "-backtrans" not in cmd


def test_build_trimal_cmd_backtrans(tmp_path: Path) -> None:
    from phyloai.pretree.trim import _build_trimal_cmd

    inp = tmp_path / "gene1.fa"
    out = tmp_path / "gene1_trim.fa"
    nt = tmp_path / "gene1.fna"
    cmd = _build_trimal_cmd(inp, out, method="gappyout", executable="trimal", backtrans_path=nt)

    assert "-gappyout" in cmd
    assert "-backtrans" in cmd
    assert str(nt) in cmd
    assert "-ignorestopcodon" in cmd


def test_build_bmge_cmd_aa(tmp_path: Path) -> None:
    from phyloai.pretree.trim import _build_bmge_cmd

    inp = tmp_path / "gene1.fa"
    out = tmp_path / "gene1_trim.fa"
    cmd = _build_bmge_cmd(
        inp, out,
        seq_type="AA",
        matrix="BLOSUM62",
        entropy=0.5,
        java_executable="java",
        bmge_jar="/path/to/BMGE.jar",
    )

    assert "java" in cmd
    assert "-jar" in cmd
    assert "/path/to/BMGE.jar" in cmd
    assert "-i" in cmd
    assert str(inp) in cmd
    assert "-t" in cmd
    assert "AA" in cmd
    assert "-m" in cmd
    assert "BLOSUM62" in cmd
    assert "-h" in cmd
    assert "0.5" in cmd
    assert "-of" in cmd
    assert str(out) in cmd


def test_build_bmge_cmd_codon(tmp_path: Path) -> None:
    from phyloai.pretree.trim import _build_bmge_cmd

    inp = tmp_path / "gene1.fna"
    out = tmp_path / "gene1_trim.fna"
    cmd = _build_bmge_cmd(
        inp, out,
        seq_type="CODON",
        matrix="BLOSUM90",
        entropy=0.4,
        java_executable="java",
        bmge_jar="/path/to/BMGE.jar",
    )

    assert "CODON" in cmd
    assert "BLOSUM90" in cmd


def test_build_bmge_cmd_nt(tmp_path: Path) -> None:
    from phyloai.pretree.trim import _build_bmge_cmd

    inp = tmp_path / "gene1.fna"
    out = tmp_path / "gene1_trim.fna"
    cmd = _build_bmge_cmd(
        inp, out,
        seq_type="NT",
        matrix="DNAPAM100:2",
        entropy=0.5,
        java_executable="java",
        bmge_jar="/path/to/BMGE.jar",
    )

    assert "DNA" in cmd  # -t DNA for NT
    assert "DNAPAM100:2" in cmd


def test_build_clipkit_cmd_aa_only(tmp_path: Path) -> None:
    from phyloai.pretree.trim import _build_clipkit_cmd

    inp = tmp_path / "gene1.fa"
    out = tmp_path / "gene1_trim.fa"
    cmd = _build_clipkit_cmd(
        inp, out,
        mode="smart-gap",
        codon=False,
        log_path=None,
        executable="clipkit",
    )

    assert cmd[0] == "clipkit"
    assert str(inp) in cmd
    assert "-o" in cmd
    assert str(out) in cmd
    assert "-m" in cmd
    assert "smart-gap" in cmd
    assert "--codon" not in cmd
    assert "-l" not in cmd


def test_build_clipkit_cmd_codon_mode(tmp_path: Path) -> None:
    from phyloai.pretree.trim import _build_clipkit_cmd

    inp = tmp_path / "gene1.fna"
    out = tmp_path / "gene1_trim.fna"
    cmd = _build_clipkit_cmd(
        inp, out,
        mode="smart-gap",
        codon=True,
        log_path=None,
        executable="clipkit",
    )

    assert "--codon" in cmd


def test_build_clipkit_cmd_with_log(tmp_path: Path) -> None:
    from phyloai.pretree.trim import _build_clipkit_cmd

    inp = tmp_path / "gene1.fa"
    out = tmp_path / "gene1_trim.fa"
    log = tmp_path / "gene1_trim.fa.log"
    cmd = _build_clipkit_cmd(
        inp, out,
        mode="smart-gap",
        codon=False,
        log_path=log,
        executable="clipkit",
    )

    assert "-l" in cmd


def test_parse_clipkit_log(tmp_path: Path) -> None:
    from phyloai.pretree.trim import _parse_clipkit_log

    log_content = (
        "1 keep constant 0.0\n"
        "2 trim other 0.9\n"
        "3 trim other 0.9\n"
        "4 keep constant 0.0\n"
        "5 keep constant 0.0\n"
    )
    log_path = tmp_path / "test.log"
    log_path.write_text(log_content)

    kept = _parse_clipkit_log(log_path)
    # 0-based indices of kept columns: positions 0, 3, 4
    assert kept == [0, 3, 4]


def test_project_columns_onto_nt_msa() -> None:
    from phyloai.pretree.trim import _project_columns_onto_nt_msa
    from Bio.SeqRecord import SeqRecord
    from Bio.Seq import Seq

    # AA kept cols [0, 3, 4] → NT cols [0,1,2, 9,10,11, 12,13,14]
    codon_records = [
        SeqRecord(Seq("ATGGCTTCTACTAAA"), id="seq1", description=""),  # 5 codons = 15 nt
        SeqRecord(Seq("ATGGCTTCTACT---"), id="seq2", description=""),
    ]
    kept_aa_cols = [0, 3, 4]
    result = _project_columns_onto_nt_msa(codon_records, kept_aa_cols)

    assert len(result) == 2
    # seq1: codons 0,3,4 → ATG + ACT + AAA = ATGACTAAA
    assert str(result[0].seq) == "ATGACTAAA"
    # seq2: codons 0,3,4 → ATG + ACT + --- = ATGACT---
    assert str(result[1].seq) == "ATGACT---"


def test_translate_codon_msa() -> None:
    from phyloai.pretree.trim import _translate_codon_msa
    from Bio.SeqRecord import SeqRecord
    from Bio.Seq import Seq

    codon_records = [
        SeqRecord(Seq("ATGGCTTCT"), id="seq1", description=""),  # M-A-S
        SeqRecord(Seq("ATG------"), id="seq2", description=""),  # M---
    ]
    aa_records = _translate_codon_msa(codon_records)

    assert str(aa_records[0].seq) == "MAS"
    assert str(aa_records[1].seq) == "M--"
    assert aa_records[0].id == "seq1"
```

- [ ] **Step 2.2: Run tests to confirm they fail**

```bash
python -m pytest tests/pretree/test_trim.py -v 2>&1 | head -20
```

Expected: ImportError — `trim` module does not exist.

- [ ] **Step 2.3: Create `phyloai/pretree/trim.py` with constants, scanning, builders**

```python
"""Batch MSA trimming using trimAl, BMGE, or ClipKIT."""

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
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from phyloai.core.checkpoint import Checkpoint
from phyloai.core.schema import COMMON_ALIGNMENT_EXTENSIONS
from phyloai.core.sequence_output_validation import validate_fasta_output
from phyloai.core.sequence_normalization import detect_seq_type, validate_codon_msa


INPUT_EXTENSIONS = {
    ext for ext in COMMON_ALIGNMENT_EXTENSIONS
    if ext in {".fa", ".fas", ".fasta", ".faa", ".fna"}
}

TRIMAL_METHODS = {"automated1", "gappyout", "strict", "strictplus"}
BMGE_SEQ_TYPE_MAP = {"AA": "AA", "NT": "DNA", "CODON": "CODON"}
CHECKPOINT_FLUSH_INTERVAL = 2.0


def _scan_input(seq_dir: Path) -> tuple[list[Path], list[dict[str, str]]]:
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


def _build_trimal_cmd(
    input_file: Path,
    output_file: Path,
    method: str,
    executable: str = "trimal",
    backtrans_path: Path | None = None,
) -> list[str]:
    cmd = [
        executable,
        "-in", str(input_file),
        "-out", str(output_file),
        f"-{method}",
    ]
    if backtrans_path is not None:
        cmd += ["-backtrans", str(backtrans_path), "-ignorestopcodon"]
    return cmd


def _build_bmge_cmd(
    input_file: Path,
    output_file: Path,
    seq_type: str,
    matrix: str,
    entropy: float,
    java_executable: str,
    bmge_jar: str,
) -> list[str]:
    bmge_type = BMGE_SEQ_TYPE_MAP.get(seq_type, "AA")
    return [
        java_executable, "-jar", bmge_jar,
        "-i", str(input_file),
        "-t", bmge_type,
        "-m", matrix,
        "-h", str(entropy),
        "-of", str(output_file),
    ]


def _build_clipkit_cmd(
    input_file: Path,
    output_file: Path,
    mode: str,
    codon: bool,
    log_path: Path | None,
    executable: str = "clipkit",
) -> list[str]:
    cmd = [executable, str(input_file), "-o", str(output_file), "-m", mode]
    if codon:
        cmd.append("--codon")
    if log_path is not None:
        cmd.append("-l")
    return cmd


def _parse_clipkit_log(log_path: Path) -> list[int]:
    """Return 0-based indices of columns marked 'keep' in a ClipKIT log file."""
    kept: list[int] = []
    with open(log_path) as fh:
        for line in fh:
            parts = line.strip().split()
            if len(parts) >= 2 and parts[1] == "keep":
                kept.append(int(parts[0]) - 1)  # 1-based → 0-based
    return kept


def _project_columns_onto_nt_msa(
    codon_records: list[SeqRecord],
    kept_aa_cols: list[int],
) -> list[SeqRecord]:
    """Project kept AA column indices onto a codon-aligned NT MSA.

    Each AA column i maps to NT columns i*3, i*3+1, i*3+2.
    """
    nt_kept: list[int] = []
    for i in kept_aa_cols:
        nt_kept.extend([i * 3, i * 3 + 1, i * 3 + 2])

    result: list[SeqRecord] = []
    for rec in codon_records:
        seq = str(rec.seq)
        new_seq = "".join(seq[j] for j in nt_kept if j < len(seq))
        result.append(SeqRecord(Seq(new_seq), id=rec.id, description=rec.description))
    return result


def _translate_codon_msa(codon_records: list[SeqRecord]) -> list[SeqRecord]:
    """Translate a codon-aligned NT MSA to an AA MSA.

    Gap triplets (---) are translated to '-'. Partial codons at end are dropped.
    """
    result: list[SeqRecord] = []
    for rec in codon_records:
        seq = str(rec.seq)
        aa_chars: list[str] = []
        for i in range(0, len(seq) - 2, 3):
            codon = seq[i:i+3]
            if set(codon) == {"-"}:
                aa_chars.append("-")
            elif "-" in codon:
                aa_chars.append("-")
            else:
                try:
                    aa = str(Seq(codon).translate())
                    aa_chars.append(aa if aa != "*" else "-")
                except Exception:
                    aa_chars.append("X")
        result.append(SeqRecord(Seq("".join(aa_chars)), id=rec.id, description=rec.description))
    return result
```

- [ ] **Step 2.4: Run tests to confirm they pass**

```bash
python -m pytest tests/pretree/test_trim.py -v -k "not run_trim"
```

Expected: all builder/scanner/utility tests PASS.

---

## Task 3: trimAl per-gene worker

**Files:**
- Modify: `phyloai/pretree/trim.py`
- Modify: `tests/pretree/test_trim.py`

- [ ] **Step 3.1: Write failing tests for trimAl worker**

Append to `tests/pretree/test_trim.py`:

```python
def test_trim_one_trimal_dry_run(tmp_path: Path) -> None:
    from phyloai.pretree.trim import _trim_one_trimal

    msa = tmp_path / "gene1.faa"
    msa.write_text(">seq1\nMKTPQ\n>seq2\nMKTPQ\n")
    out_dir = tmp_path / "seqs"
    out_dir.mkdir()

    result = _trim_one_trimal(
        msa_path=msa,
        aa_out_dir=out_dir,
        nt_out_dir=None,
        nt_path=None,
        method="automated1",
        seq_type="AA",
        extra_args=None,
        dry_run=True,
        executable="trimal",
    )

    assert result["status"] == "dry_run"
    assert "cmd" in result


def test_trim_one_trimal_codon_dry_run(tmp_path: Path) -> None:
    from phyloai.pretree.trim import _trim_one_trimal

    msa = tmp_path / "gene1.fna"
    msa.write_text(">seq1\nATGGCTTCT\n>seq2\nATGGCT---\n")
    faa_dir = tmp_path / "seqs" / "faa"
    fna_dir = tmp_path / "seqs" / "fna"
    faa_dir.mkdir(parents=True)
    fna_dir.mkdir(parents=True)

    result = _trim_one_trimal(
        msa_path=msa,
        aa_out_dir=faa_dir,
        nt_out_dir=fna_dir,
        nt_path=None,
        method="automated1",
        seq_type="CODON",
        extra_args=None,
        dry_run=True,
        executable="trimal",
    )

    assert result["status"] == "dry_run"
```

- [ ] **Step 3.2: Run tests to confirm they fail**

```bash
python -m pytest tests/pretree/test_trim.py -v -k "trimal" 2>&1 | head -20
```

Expected: ImportError on `_trim_one_trimal`.

- [ ] **Step 3.3: Implement `_trim_one_trimal`**

Append to `phyloai/pretree/trim.py`:

```python
def _trim_one_trimal(
    msa_path: Path,
    aa_out_dir: Path,
    nt_out_dir: Path | None,
    nt_path: Path | None,
    method: str,
    seq_type: str,
    extra_args: str | None,
    dry_run: bool,
    executable: str = "trimal",
) -> dict[str, Any]:
    """Trim one gene MSA with trimAl.

    Modes:
    - AA/NT-only: trimal -in <msa> -out <aa_out> -<method>
    - CODON: Python translate → temp AA MSA + strip gaps → temp unaligned CDS;
             trimal -in <temp_aa_msa> -out <faa/gene.fa> -<method> -backtrans <temp_cds>
    - AA+NT (Mode 4): trimal -in <msa> -out <faa/gene.fa> -<method> -backtrans <nt_path>
    """
    gene_stem = msa_path.stem
    aa_out = aa_out_dir / f"{gene_stem}.fa"
    nt_out = (nt_out_dir / f"{gene_stem}.fa") if nt_out_dir else None
    length_before = _read_msa_col_count(msa_path)

    # --- CODON mode: build temp AA MSA + temp unaligned CDS ---
    if seq_type == "CODON":
        try:
            codon_records = list(SeqIO.parse(str(msa_path), "fasta"))
        except Exception as exc:
            return {"status": "skipped", "input": str(msa_path),
                    "reason": f"could not parse codon MSA: {exc}"}

        seqs_dict = {r.id: str(r.seq) for r in codon_records}
        validation = validate_codon_msa(seqs_dict)
        if validation.skip:
            return {"status": "skipped", "input": str(msa_path),
                    "reason": "; ".join(validation.warnings)}

        # Rebuild records with (possibly terminal-stop-stripped) sequences
        codon_records = [
            SeqRecord(Seq(validation.sequences[r.id]), id=r.id, description=r.description)
            for r in codon_records
        ]
        aa_records = _translate_codon_msa(codon_records)

        if dry_run:
            return {
                "status": "dry_run",
                "input": str(msa_path),
                "cmd": f"[CODON] trimal -in <temp_aa_msa> -out {aa_out} -{method} "
                       f"-backtrans <temp_unaligned_cds> -ignorestopcodon",
                "codon_warnings": validation.warnings,
            }

        with tempfile.TemporaryDirectory(prefix="phyloai_trim_") as tmpdir:
            tmp_aa_msa = Path(tmpdir) / f"{gene_stem}_aa_msa.fa"
            tmp_cds = Path(tmpdir) / f"{gene_stem}_cds_unaligned.fa"

            # Write temp AA MSA
            SeqIO.write(aa_records, str(tmp_aa_msa), "fasta")

            # Write temp unaligned CDS (strip gaps from codon records)
            unaligned_cds = [
                SeqRecord(Seq(str(r.seq).replace("-", "")), id=r.id, description=r.description)
                for r in codon_records
            ]
            SeqIO.write(unaligned_cds, str(tmp_cds), "fasta")

            # Build command
            aa_out.parent.mkdir(parents=True, exist_ok=True)
            if nt_out:
                nt_out.parent.mkdir(parents=True, exist_ok=True)

            fna_target = nt_out if nt_out else aa_out.parent / f"{gene_stem}_nt.fa"
            cmd = _build_trimal_cmd(tmp_aa_msa, fna_target, method=method,
                                    executable=executable, backtrans_path=tmp_cds)
            _apply_extra_args(cmd, extra_args)

            return _run_trimal_cmd(cmd, msa_path=msa_path, aa_out=aa_out,
                                   nt_out=nt_out, codon_warnings=validation.warnings,
                                   fna_target=fna_target, mode="CODON",
                                   length_before=length_before)

    # --- Mode 4: AA + NT backtrans ---
    if nt_path is not None and seq_type == "AA":
        aa_out.parent.mkdir(parents=True, exist_ok=True)
        if nt_out:
            nt_out.parent.mkdir(parents=True, exist_ok=True)

        # Trim AA separately first
        cmd_aa = _build_trimal_cmd(msa_path, aa_out, method=method, executable=executable)
        _apply_extra_args(cmd_aa, extra_args)

        if dry_run:
            cmd_nt = _build_trimal_cmd(msa_path, nt_out or Path("/dev/null"), method=method,
                                       executable=executable, backtrans_path=nt_path)
            return {"status": "dry_run", "input": str(msa_path),
                    "cmd": " ".join(cmd_aa) + " && " + " ".join(cmd_nt)}

        # Run AA trim
        start = time.monotonic()
        proc_aa = subprocess.run(cmd_aa, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        if proc_aa.returncode != 0:
            return {"status": "skipped", "input": str(msa_path),
                    "reason": f"trimal (AA) exited {proc_aa.returncode}: {proc_aa.stderr[:300]}",
                    "tool_stderr": proc_aa.stderr, "wall_time": time.monotonic() - start}

        # Run NT backtrans
        cmd_nt = _build_trimal_cmd(msa_path, nt_out or aa_out.parent / f"{gene_stem}_nt.fa",
                                   method=method, executable=executable, backtrans_path=nt_path)
        _apply_extra_args(cmd_nt, extra_args)
        proc_nt = subprocess.run(cmd_nt, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        wall_time = time.monotonic() - start

        if proc_nt.returncode != 0:
            return {"status": "skipped", "input": str(msa_path),
                    "reason": f"trimal (NT backtrans) exited {proc_nt.returncode}: {proc_nt.stderr[:300]}",
                    "tool_stderr": proc_nt.stderr, "wall_time": wall_time}

        return _make_success_result(msa_path, aa_out, nt_out,
                                    cmd=" ".join(cmd_nt), wall_time=wall_time,
                                    tool_stderr=proc_nt.stderr, warnings=[],
                                    length_before=length_before)

    # --- Mode 1/2: AA-only or NT-only ---
    aa_out.parent.mkdir(parents=True, exist_ok=True)
    cmd = _build_trimal_cmd(msa_path, aa_out, method=method, executable=executable)
    _apply_extra_args(cmd, extra_args)

    if dry_run:
        return {"status": "dry_run", "input": str(msa_path), "cmd": cmd}

    start = time.monotonic()
    proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    wall_time = time.monotonic() - start
    if proc.returncode != 0:
        return {"status": "skipped", "input": str(msa_path),
                "reason": f"trimal exited {proc.returncode}: {proc.stderr[:300]}",
                "tool_stderr": proc.stderr, "wall_time": wall_time}

    return _make_success_result(msa_path, aa_out, None,
                                cmd=" ".join(cmd), wall_time=wall_time,
                                tool_stderr=proc.stderr, warnings=[],
                                length_before=length_before)


def _apply_extra_args(cmd: list[str], extra_args: str | None) -> None:
    """Tokenize extra_args and append to cmd in-place.

    Extra-args are appended after the internally-constructed command.  When a
    tool processes duplicate flags, later values typically win (trimAl, BMGE,
    ClipKIT all behave this way), which achieves the spec's "extra-wins merge"
    semantics without fragile token-level surgery.

    Scope: flag replacement via removal is NOT attempted because reliably
    distinguishing bool flags from valued flags across three different tools
    requires tool-specific knowledge.  Simple append is correct and safe.
    """
    if not extra_args:
        return
    cmd.extend(shlex.split(extra_args))


def _run_trimal_cmd(
    cmd: list[str],
    *,
    msa_path: Path,
    aa_out: Path,
    nt_out: Path | None,
    codon_warnings: list[str],
    fna_target: Path,
    mode: str,
    length_before: int = 0,
) -> dict[str, Any]:
    start = time.monotonic()
    try:
        proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    except Exception as exc:
        return {"status": "skipped", "input": str(msa_path),
                "reason": str(exc), "wall_time": time.monotonic() - start}
    wall_time = time.monotonic() - start
    if proc.returncode != 0:
        return {"status": "skipped", "input": str(msa_path),
                "reason": f"trimal exited {proc.returncode}: {proc.stderr[:300]}",
                "tool_stderr": proc.stderr, "wall_time": wall_time}

    # For CODON mode, fna_target is the trimmed NT output; aa_out must be derived
    # from it by Python translation. Copy or re-parse to produce AA.
    if mode == "CODON" and fna_target.exists():
        try:
            trimmed_codon = list(SeqIO.parse(str(fna_target), "fasta"))
            aa_records = _translate_codon_msa(trimmed_codon)
            aa_out.parent.mkdir(parents=True, exist_ok=True)
            SeqIO.write(aa_records, str(aa_out), "fasta")
        except Exception as exc:
            return {"status": "skipped", "input": str(msa_path),
                    "reason": f"codon→AA translation failed: {exc}"}

    return _make_success_result(
        msa_path, aa_out, nt_out,
        cmd=" ".join(cmd), wall_time=wall_time,
        tool_stderr=proc.stderr, warnings=list(codon_warnings),
        length_before=length_before,
    )


def _read_msa_col_count(path: Path) -> int:
    """Return alignment column count (length of first sequence) from a FASTA MSA.

    Returns 0 if the file cannot be read or is empty.
    """
    try:
        for rec in SeqIO.parse(str(path), "fasta"):
            return len(rec.seq)
    except Exception:
        pass
    return 0


def _make_success_result(
    msa_path: Path,
    aa_out: Path,
    nt_out: Path | None,
    *,
    cmd: str,
    wall_time: float,
    tool_stderr: str,
    warnings: list[str],
    length_before: int = 0,
) -> dict[str, Any]:
    length_after = _read_msa_col_count(aa_out) if aa_out.exists() else 0
    return {
        "status": "success",
        "input": str(msa_path),
        "output_aa": str(aa_out),
        "output_nt": str(nt_out) if nt_out else None,
        "tool_cmd": cmd,
        "tool_stderr": tool_stderr,
        "wall_time": wall_time,
        "warnings": warnings,
        "length_before": length_before,
        "length_after": length_after,
    }
```

- [ ] **Step 3.4: Run tests to confirm they pass**

```bash
python -m pytest tests/pretree/test_trim.py -v -k "trimal"
```

Expected: all trimAl dry-run tests PASS.

---

## Task 4: BMGE and ClipKIT per-gene workers

**Files:**
- Modify: `phyloai/pretree/trim.py`
- Modify: `tests/pretree/test_trim.py`

- [ ] **Step 4.1: Write failing tests for BMGE and ClipKIT workers**

Append to `tests/pretree/test_trim.py`:

```python
def test_trim_one_bmge_dry_run(tmp_path: Path) -> None:
    from phyloai.pretree.trim import _trim_one_bmge

    msa = tmp_path / "gene1.faa"
    msa.write_text(">seq1\nMKTPQ\n>seq2\nMKTPQ\n")
    out_dir = tmp_path / "seqs"
    out_dir.mkdir()

    result = _trim_one_bmge(
        msa_path=msa,
        aa_out_dir=out_dir,
        nt_out_dir=None,
        seq_type="AA",
        matrix="BLOSUM62",
        entropy=0.5,
        extra_args=None,
        dry_run=True,
        java_executable="java",
        bmge_jar="/fake/BMGE.jar",
    )

    assert result["status"] == "dry_run"
    assert "cmd" in result


def test_trim_one_bmge_codon_dry_run(tmp_path: Path) -> None:
    from phyloai.pretree.trim import _trim_one_bmge

    msa = tmp_path / "gene1.fna"
    msa.write_text(">seq1\nATGGCTTCT\n>seq2\nATGGCT---\n")
    faa_dir = tmp_path / "seqs" / "faa"
    fna_dir = tmp_path / "seqs" / "fna"
    faa_dir.mkdir(parents=True)
    fna_dir.mkdir(parents=True)

    result = _trim_one_bmge(
        msa_path=msa,
        aa_out_dir=faa_dir,
        nt_out_dir=fna_dir,
        seq_type="CODON",
        matrix="BLOSUM62",
        entropy=0.5,
        extra_args=None,
        dry_run=True,
        java_executable="java",
        bmge_jar="/fake/BMGE.jar",
    )

    assert result["status"] == "dry_run"


def test_trim_one_clipkit_dry_run(tmp_path: Path) -> None:
    from phyloai.pretree.trim import _trim_one_clipkit

    msa = tmp_path / "gene1.faa"
    msa.write_text(">seq1\nMKTPQ\n>seq2\nMKTPQ\n")
    out_dir = tmp_path / "seqs"
    out_dir.mkdir()

    result = _trim_one_clipkit(
        msa_path=msa,
        aa_out_dir=out_dir,
        nt_out_dir=None,
        nt_path=None,
        mode="smart-gap",
        seq_type="AA",
        extra_args=None,
        dry_run=True,
        executable="clipkit",
    )

    assert result["status"] == "dry_run"
    assert "cmd" in result


def test_trim_one_clipkit_mode4_dry_run(tmp_path: Path) -> None:
    from phyloai.pretree.trim import _trim_one_clipkit

    msa = tmp_path / "gene1.faa"
    msa.write_text(">seq1\nMKTPQ\n>seq2\nMKTPQ\n")
    nt_msa = tmp_path / "gene1.fna"
    nt_msa.write_text(">seq1\nATGAAGACCCCTCAA\n>seq2\nATGAAGACCCCTCAA\n")
    faa_dir = tmp_path / "seqs" / "faa"
    fna_dir = tmp_path / "seqs" / "fna"
    faa_dir.mkdir(parents=True)
    fna_dir.mkdir(parents=True)

    result = _trim_one_clipkit(
        msa_path=msa,
        aa_out_dir=faa_dir,
        nt_out_dir=fna_dir,
        nt_path=nt_msa,
        mode="smart-gap",
        seq_type="AA",
        extra_args=None,
        dry_run=True,
        executable="clipkit",
    )

    assert result["status"] == "dry_run"
```

- [ ] **Step 4.2: Run tests to confirm they fail**

```bash
python -m pytest tests/pretree/test_trim.py -v -k "bmge or clipkit" 2>&1 | head -20
```

Expected: ImportError on `_trim_one_bmge` and `_trim_one_clipkit`.

- [ ] **Step 4.3: Implement `_trim_one_bmge` and `_trim_one_clipkit`**

Append to `phyloai/pretree/trim.py`:

```python
def _trim_one_bmge(
    msa_path: Path,
    aa_out_dir: Path,
    nt_out_dir: Path | None,
    seq_type: str,
    matrix: str,
    entropy: float,
    extra_args: str | None,
    dry_run: bool,
    java_executable: str = "java",
    bmge_jar: str = "BMGE.jar",
) -> dict[str, Any]:
    """Trim one MSA with BMGE.

    Modes:
    - AA/NT: direct BMGE run → one output file in aa_out_dir
    - CODON: BMGE -t CODON → fna_out; Python translate → faa_out
    - Mode 4 auto-downgrade (seq_type==AA, nt_out_dir set): treated as CODON
      (caller ensures msa_path is the codon-aligned NT MSA from --nt-dir)
    """
    gene_stem = msa_path.stem
    is_dual = nt_out_dir is not None
    aa_out = aa_out_dir / f"{gene_stem}.fa"
    nt_out = (nt_out_dir / f"{gene_stem}.fa") if nt_out_dir else None
    length_before = _read_msa_col_count(msa_path)

    # Determine effective seq_type for BMGE: CODON mode is used for dual output
    effective_seq_type = "CODON" if is_dual else seq_type
    primary_out = nt_out if is_dual else aa_out

    cmd = _build_bmge_cmd(
        msa_path, primary_out,
        seq_type=effective_seq_type,
        matrix=matrix,
        entropy=entropy,
        java_executable=java_executable,
        bmge_jar=bmge_jar,
    )
    _apply_extra_args(cmd, extra_args)

    if dry_run:
        return {"status": "dry_run", "input": str(msa_path), "cmd": cmd}

    primary_out.parent.mkdir(parents=True, exist_ok=True)
    if aa_out_dir:
        aa_out.parent.mkdir(parents=True, exist_ok=True)

    # Validate CODON input
    if effective_seq_type == "CODON":
        try:
            codon_records = list(SeqIO.parse(str(msa_path), "fasta"))
        except Exception as exc:
            return {"status": "skipped", "input": str(msa_path),
                    "reason": f"could not parse codon MSA: {exc}"}
        seqs_dict = {r.id: str(r.seq) for r in codon_records}
        validation = validate_codon_msa(seqs_dict)
        if validation.skip:
            return {"status": "skipped", "input": str(msa_path),
                    "reason": "; ".join(validation.warnings)}
        codon_warnings = validation.warnings
    else:
        codon_warnings = []

    start = time.monotonic()
    try:
        proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    except Exception as exc:
        return {"status": "skipped", "input": str(msa_path),
                "reason": str(exc), "wall_time": time.monotonic() - start}
    wall_time = time.monotonic() - start

    if proc.returncode != 0:
        return {"status": "skipped", "input": str(msa_path),
                "reason": f"BMGE exited {proc.returncode}: {proc.stderr[:300]}",
                "tool_stderr": proc.stderr, "wall_time": wall_time}

    # For CODON mode: translate trimmed NT → AA
    if is_dual and nt_out and nt_out.exists():
        try:
            trimmed_codon = list(SeqIO.parse(str(nt_out), "fasta"))
            aa_records = _translate_codon_msa(trimmed_codon)
            SeqIO.write(aa_records, str(aa_out), "fasta")
        except Exception as exc:
            return {"status": "skipped", "input": str(msa_path),
                    "reason": f"codon→AA translation failed: {exc}"}

    return _make_success_result(
        msa_path, aa_out, nt_out,
        cmd=" ".join(cmd), wall_time=wall_time,
        tool_stderr=proc.stderr, warnings=list(codon_warnings),
        length_before=length_before,
    )


def _trim_one_clipkit(
    msa_path: Path,
    aa_out_dir: Path,
    nt_out_dir: Path | None,
    nt_path: Path | None,
    mode: str,
    seq_type: str,
    extra_args: str | None,
    dry_run: bool,
    executable: str = "clipkit",
) -> dict[str, Any]:
    """Trim one MSA with ClipKIT.

    Modes:
    - AA/NT-only: clipkit without -l, no NT output
    - CODON: clipkit --codon → fna_out; Python translate → faa_out
    - Mode 4 (AA + nt_path): clipkit -l on AA → parse log → project onto codon NT MSA
    """
    gene_stem = msa_path.stem
    aa_out = aa_out_dir / f"{gene_stem}.fa"
    nt_out = (nt_out_dir / f"{gene_stem}.fa") if nt_out_dir else None
    length_before = _read_msa_col_count(msa_path)
    is_mode4 = nt_path is not None and seq_type == "AA"
    is_codon = seq_type == "CODON"

    if is_mode4:
        # Mode 4: trim AA with -l, parse log, project onto codon NT MSA
        if dry_run:
            cmd = _build_clipkit_cmd(msa_path, aa_out, mode=mode, codon=False,
                                     log_path=Path("gene.log"), executable=executable)
            return {"status": "dry_run", "input": str(msa_path), "cmd": cmd}

        aa_out.parent.mkdir(parents=True, exist_ok=True)
        if nt_out:
            nt_out.parent.mkdir(parents=True, exist_ok=True)

        with tempfile.TemporaryDirectory(prefix="phyloai_trim_") as tmpdir:
            tmp_aa_out = Path(tmpdir) / f"{gene_stem}.fa"
            log_path = Path(str(tmp_aa_out) + ".log")
            cmd = _build_clipkit_cmd(msa_path, tmp_aa_out, mode=mode, codon=False,
                                     log_path=log_path, executable=executable)
            _apply_extra_args(cmd, extra_args)

            start = time.monotonic()
            proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            wall_time = time.monotonic() - start
            if proc.returncode != 0:
                return {"status": "skipped", "input": str(msa_path),
                        "reason": f"clipkit exited {proc.returncode}: {proc.stderr[:300]}",
                        "tool_stderr": proc.stderr, "wall_time": wall_time}

            # Copy trimmed AA to final location
            shutil.copy2(tmp_aa_out, aa_out)

            # Parse log and project onto codon NT MSA
            kept_cols = _parse_clipkit_log(log_path)
            try:
                codon_records = list(SeqIO.parse(str(nt_path), "fasta"))
                projected = _project_columns_onto_nt_msa(codon_records, kept_cols)
                if nt_out:
                    SeqIO.write(projected, str(nt_out), "fasta")
            except Exception as exc:
                return {"status": "skipped", "input": str(msa_path),
                        "reason": f"NT projection failed: {exc}"}

        return _make_success_result(
            msa_path, aa_out, nt_out,
            cmd=" ".join(cmd), wall_time=wall_time,
            tool_stderr=proc.stderr, warnings=[],
            length_before=length_before,
        )

    if is_codon:
        # CODON mode: clipkit --codon → fna; Python translate → faa
        nt_primary = nt_out if nt_out else aa_out_dir / f"{gene_stem}.fa"
        cmd = _build_clipkit_cmd(msa_path, nt_primary, mode=mode, codon=True,
                                 log_path=None, executable=executable)
        _apply_extra_args(cmd, extra_args)

        if dry_run:
            return {"status": "dry_run", "input": str(msa_path), "cmd": cmd}

        nt_primary.parent.mkdir(parents=True, exist_ok=True)
        aa_out.parent.mkdir(parents=True, exist_ok=True)

        # Validate codon MSA
        try:
            codon_records = list(SeqIO.parse(str(msa_path), "fasta"))
        except Exception as exc:
            return {"status": "skipped", "input": str(msa_path),
                    "reason": f"could not parse codon MSA: {exc}"}
        seqs_dict = {r.id: str(r.seq) for r in codon_records}
        validation = validate_codon_msa(seqs_dict)
        if validation.skip:
            return {"status": "skipped", "input": str(msa_path),
                    "reason": "; ".join(validation.warnings)}

        start = time.monotonic()
        proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        wall_time = time.monotonic() - start
        if proc.returncode != 0:
            return {"status": "skipped", "input": str(msa_path),
                    "reason": f"clipkit exited {proc.returncode}: {proc.stderr[:300]}",
                    "tool_stderr": proc.stderr, "wall_time": wall_time}

        # Translate trimmed codon → AA
        if nt_primary.exists():
            try:
                trimmed_codon = list(SeqIO.parse(str(nt_primary), "fasta"))
                aa_records = _translate_codon_msa(trimmed_codon)
                SeqIO.write(aa_records, str(aa_out), "fasta")
            except Exception as exc:
                return {"status": "skipped", "input": str(msa_path),
                        "reason": f"codon→AA translation failed: {exc}"}

        return _make_success_result(
            msa_path, aa_out, nt_out,
            cmd=" ".join(cmd), wall_time=wall_time,
            tool_stderr=proc.stderr, warnings=list(validation.warnings),
            length_before=length_before,
        )

    # Mode 1/2: AA-only or NT-only — no -l needed
    aa_out.parent.mkdir(parents=True, exist_ok=True)
    cmd = _build_clipkit_cmd(msa_path, aa_out, mode=mode, codon=False,
                             log_path=None, executable=executable)
    _apply_extra_args(cmd, extra_args)

    if dry_run:
        return {"status": "dry_run", "input": str(msa_path), "cmd": cmd}

    start = time.monotonic()
    proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    wall_time = time.monotonic() - start
    if proc.returncode != 0:
        return {"status": "skipped", "input": str(msa_path),
                "reason": f"clipkit exited {proc.returncode}: {proc.stderr[:300]}",
                "tool_stderr": proc.stderr, "wall_time": wall_time}

    return _make_success_result(
        msa_path, aa_out, None,
        cmd=" ".join(cmd), wall_time=wall_time,
        tool_stderr=proc.stderr, warnings=[],
        length_before=length_before,
    )
```

- [ ] **Step 4.4: Run tests to confirm they pass**

```bash
python -m pytest tests/pretree/test_trim.py -v -k "bmge or clipkit"
```

Expected: all BMGE and ClipKIT dry-run tests PASS.

---

## Task 5: Output validation, worker dispatch, `run_trim` orchestrator

**Files:**
- Modify: `phyloai/pretree/trim.py`
- Modify: `tests/pretree/test_trim.py`

- [ ] **Step 5.1: Write failing tests for validation and orchestration helpers**

Append to `tests/pretree/test_trim.py`:

```python
def test_verify_trim_outputs_both_exist(tmp_path: Path) -> None:
    from phyloai.pretree.trim import verify_trim_outputs

    aa = tmp_path / "gene1.fa"
    nt = tmp_path / "gene1_nt.fa"
    aa.write_text(">seq1\nMKT\n>seq2\nMKT\n")
    nt.write_text(">seq1\nATGAAGACT\n>seq2\nATGAAGACT\n")

    assert verify_trim_outputs(aa, nt) is True


def test_verify_trim_outputs_missing_aa(tmp_path: Path) -> None:
    from phyloai.pretree.trim import verify_trim_outputs

    aa = tmp_path / "gene1.fa"
    assert verify_trim_outputs(aa, None) is False


def test_verify_trim_outputs_empty_aa(tmp_path: Path) -> None:
    from phyloai.pretree.trim import verify_trim_outputs

    aa = tmp_path / "gene1.fa"
    aa.write_text("")
    assert verify_trim_outputs(aa, None) is False


def test_detect_trim_seq_type_auto(tmp_path: Path) -> None:
    from phyloai.pretree.trim import _detect_seq_type_from_files

    f = tmp_path / "gene1.fa"
    f.write_text(">seq1\nMKTPQWER\n>seq2\nMKTPQWER\n")
    result = _detect_seq_type_from_files([f])
    assert result == "AA"


def test_run_trim_validation_codon_with_nt_dir_raises(tmp_path: Path) -> None:
    from phyloai.pretree.trim import run_trim

    msa_dir = tmp_path / "msa"
    msa_dir.mkdir()
    (msa_dir / "gene1.fa").write_text(">a\nATGGCT\n")

    with pytest.raises(ValueError, match="CODON mode does not use --nt-dir"):
        run_trim(
            msa_dir=msa_dir,
            output_dir=tmp_path / "out",
            tool="trimal",
            seq_type="CODON",
            nt_dir=tmp_path / "nt",
            threads=1,
        )


def test_run_trim_validation_overwrite_resume_mutual_exclusive(tmp_path: Path) -> None:
    from phyloai.pretree.trim import run_trim

    msa_dir = tmp_path / "msa"
    msa_dir.mkdir()
    (msa_dir / "gene1.fa").write_text(">a\nMKT\n")

    with pytest.raises(ValueError, match="mutually exclusive"):
        run_trim(
            msa_dir=msa_dir,
            output_dir=tmp_path / "out",
            tool="trimal",
            seq_type="AA",
            overwrite=True,
            resume=True,
            threads=1,
        )


def test_run_trim_dry_run_returns_payload(tmp_path: Path) -> None:
    from phyloai.pretree.trim import run_trim

    msa_dir = tmp_path / "msa"
    msa_dir.mkdir()
    (msa_dir / "gene1.fa").write_text(">a\nMKT\n>b\nMKT\n")

    payload = run_trim(
        msa_dir=msa_dir,
        output_dir=tmp_path / "out",
        tool="trimal",
        seq_type="AA",
        dry_run=True,
        threads=1,
    )

    assert payload["status"] == "success"
    assert payload["data"]["summary"]["n_input_files"] == 1
```

- [ ] **Step 5.2: Run tests to confirm they fail**

```bash
python -m pytest tests/pretree/test_trim.py -v -k "verify or detect or run_trim" 2>&1 | head -20
```

Expected: ImportError on `verify_trim_outputs`, `_detect_seq_type_from_files`, `run_trim`.

- [ ] **Step 5.3: Implement `verify_trim_outputs`, `_detect_seq_type_from_files`, module-level worker, and `run_trim`**

Append to `phyloai/pretree/trim.py`:

```python
def verify_trim_outputs(aa_path: Path, nt_path: Path | None) -> bool:
    """Check that trim output files exist and are valid FASTA."""
    if not aa_path.exists() or aa_path.stat().st_size == 0:
        return False
    aa_result = validate_fasta_output(aa_path, require_aligned=True)
    if not aa_result.ok:
        return False
    if nt_path is None:
        return True
    if not nt_path.exists() or nt_path.stat().st_size == 0:
        return False
    nt_result = validate_fasta_output(nt_path, require_aligned=True)
    return nt_result.ok


def _detect_seq_type_from_files(files: list[Path], max_files: int = 3) -> str:
    sequences: list[str] = []
    for f in files[:max_files]:
        try:
            for rec in SeqIO.parse(str(f), "fasta"):
                sequences.append(str(rec.seq))
                if len(sequences) >= 10:
                    break
        except Exception:
            continue
        if len(sequences) >= 10:
            break
    return detect_seq_type(sequences) if sequences else "AA"


def _trim_one_worker(
    args: tuple[
        Path, Path, Path | None, Path | None, str, str, str,
        str | None, str | None, float, bool, str, str, str,
    ]
) -> dict[str, Any]:
    """Module-level worker function required for ProcessPoolExecutor on macOS (spawn).

    args = (msa_path, aa_out_dir, nt_out_dir, nt_path,
            tool, seq_type, trimal_method,
            clipkit_mode, bmge_matrix, bmge_entropy, extra_args,
            dry_run, trimal_exe, java_exe, bmge_jar, clipkit_exe)
    """
    (msa_path, aa_out_dir, nt_out_dir, nt_path,
     tool, seq_type, trimal_method,
     clipkit_mode, bmge_matrix, bmge_entropy, extra_args,
     dry_run, trimal_exe, java_exe, bmge_jar, clipkit_exe) = args

    if tool == "trimal":
        return _trim_one_trimal(
            msa_path=msa_path, aa_out_dir=aa_out_dir, nt_out_dir=nt_out_dir,
            nt_path=nt_path, method=trimal_method, seq_type=seq_type,
            extra_args=extra_args, dry_run=dry_run, executable=trimal_exe,
        )
    elif tool == "bmge":
        return _trim_one_bmge(
            msa_path=msa_path, aa_out_dir=aa_out_dir, nt_out_dir=nt_out_dir,
            seq_type=seq_type, matrix=bmge_matrix, entropy=bmge_entropy,
            extra_args=extra_args, dry_run=dry_run,
            java_executable=java_exe, bmge_jar=bmge_jar,
        )
    else:  # clipkit
        return _trim_one_clipkit(
            msa_path=msa_path, aa_out_dir=aa_out_dir, nt_out_dir=nt_out_dir,
            nt_path=nt_path, mode=clipkit_mode, seq_type=seq_type,
            extra_args=extra_args, dry_run=dry_run, executable=clipkit_exe,
        )


def _detect_trim_tool_versions(
    tool: str,
    trimal_path: Path | None,
    bmge_path: Path | None,
    clipkit_path: Path | None,
) -> dict[str, str]:
    """Return version strings for the active tool (and java for BMGE)."""
    from phyloai.core.env import TOOL_REGISTRY, ToolEnv, ToolStatus

    tool_paths: dict[str, Path] = {}
    if trimal_path:
        tool_paths["trimal"] = trimal_path
    if bmge_path:
        tool_paths["bmge"] = bmge_path
    if clipkit_path:
        tool_paths["clipkit"] = clipkit_path

    env = ToolEnv(tool_paths=tool_paths)
    names = [tool]
    if tool == "bmge":
        names.append("java")

    versions: dict[str, str] = {}
    for name in names:
        meta = TOOL_REGISTRY.get(name, {})
        info = env._detect_tool(
            name,
            version_flag=meta.get("version_flag", ""),
            version_args=meta.get("version_args"),
            bundled=meta.get("bundled", False),
            bundled_dir=meta.get("bundled_dir"),
            bundled_executable=meta.get("bundled_executable"),
            path_aliases=meta.get("path_aliases"),
        )
        if info.status == ToolStatus.OK and info.version:
            versions[name] = info.version
    return versions


def _resolve_trim_tool_paths(
    tool: str,
    trimal_path: Path | None,
    bmge_path: Path | None,
    clipkit_path: Path | None,
    dry_run: bool,
) -> tuple[str, str, str, str]:
    """Return (trimal_exe, java_exe, bmge_jar, clipkit_exe)."""
    from phyloai.core.env import ToolEnv

    env = ToolEnv(tool_paths={
        **({"trimal": trimal_path} if trimal_path else {}),
        **({"bmge": bmge_path} if bmge_path else {}),
        **({"clipkit": clipkit_path} if clipkit_path else {}),
    })

    if tool == "trimal":
        exe = str(trimal_path) if trimal_path else ("trimal" if dry_run else str(env.require("trimal")))
        return exe, "java", "BMGE.jar", "clipkit"

    if tool == "bmge":
        bmge_jar = str(bmge_path) if bmge_path else ("BMGE.jar" if dry_run else str(env.require("bmge")))
        # Also require java
        java_path = shutil.which("java")
        if not java_path and not dry_run:
            raise FileNotFoundError("java not found on PATH; required for BMGE")
        return "trimal", java_path or "java", bmge_jar, "clipkit"

    # clipkit
    exe = str(clipkit_path) if clipkit_path else ("clipkit" if dry_run else str(env.require("clipkit")))
    return "trimal", "java", "BMGE.jar", exe


def run_trim(
    msa_dir: Path,
    output_dir: Path,
    tool: str = "trimal",
    seq_type: str = "auto",
    nt_dir: Path | None = None,
    trimal_method: str = "automated1",
    bmge_matrix: str | None = None,
    bmge_entropy: float = 0.5,
    clipkit_mode: str = "smart-gap",
    trimal_path: Path | None = None,
    bmge_path: Path | None = None,
    clipkit_path: Path | None = None,
    threads: int = 4,
    extra_args: str | None = None,
    overwrite: bool = False,
    resume: bool = False,
    dry_run: bool = False,
    quiet: bool = False,
    progress_callback: Callable[[Path], None] | None = None,
) -> dict[str, Any]:
    from phyloai.core.checkpoint import (
        load_checkpoint,
        save_checkpoint_atomic,
        validate_resume_params,
    )
    from phyloai.pretree.checkpoint_helpers import (
        build_initial_checkpoint,
        mark_task,
        plan_resume,
    )

    run_start = time.monotonic()
    global_warnings: list[str] = []

    # --- Param validation ---
    if resume and overwrite:
        raise ValueError("--overwrite and --resume are mutually exclusive.")
    if seq_type == "CODON" and nt_dir is not None:
        raise ValueError(
            "CODON mode does not use --nt-dir. "
            "Place codon-aligned MSA in --msa-dir and omit --nt-dir."
        )
    if threads < 1:
        raise ValueError("--threads must be at least 1.")

    # BMGE + Mode 4 auto-downgrade
    bmge_mode4_downgrade = False
    if tool == "bmge" and seq_type in ("AA", "auto") and nt_dir is not None:
        bmge_mode4_downgrade = True
        global_warnings.append(
            "BMGE does not support AA+NT mode directly; "
            "automatically using --nt-dir files in CODON mode (-t CODON)."
        )

    # --- Scan inputs ---
    found, scan_skipped = _scan_input(msa_dir)
    if not found and not dry_run:
        raise ValueError("No valid input MSA files found in --msa-dir.")

    # --- Resolve seq_type ---
    if seq_type == "auto":
        resolved_seq_type = _detect_seq_type_from_files(found) if found else "AA"
        global_warnings.append(f"seq_type auto-detected as '{resolved_seq_type}'.")
    else:
        resolved_seq_type = seq_type

    # --- BMGE matrix default ---
    if bmge_matrix is None:
        bmge_matrix = "BLOSUM62" if resolved_seq_type in ("AA", "CODON") else "DNAPAM100:2"

    # --- Output directories ---
    is_dual = (nt_dir is not None) or resolved_seq_type == "CODON"
    if is_dual:
        aa_out_dir = output_dir / "seqs" / "faa"
        nt_out_dir = output_dir / "seqs" / "fna"
    else:
        aa_out_dir = output_dir / "seqs"
        nt_out_dir = None

    # --- Resolve tool paths ---
    trimal_exe, java_exe, bmge_jar, clipkit_exe = _resolve_trim_tool_paths(
        tool=tool,
        trimal_path=trimal_path,
        bmge_path=bmge_path,
        clipkit_path=clipkit_path,
        dry_run=dry_run,
    )

    # --- Resolved params for checkpoint ---
    resolved_params: dict[str, Any] = {
        "msa_dir": str(msa_dir),
        "nt_dir": str(nt_dir) if nt_dir else None,
        "seq_type": seq_type,
        "effective_seq_type": resolved_seq_type,
        "tool": tool,
        "trimal_method": trimal_method,
        "bmge_matrix": bmge_matrix,
        "bmge_entropy": bmge_entropy,
        "clipkit_mode": clipkit_mode,
        "threads": int(threads),
        "extra_args": extra_args,
        "output_dir": str(output_dir),
        "bmge_mode4_downgrade": bmge_mode4_downgrade,
    }

    # --- Output dir + checkpoint setup ---
    ckpt_path = output_dir / "checkpoint.json"
    checkpoint = None
    to_run_ids: list[str] | None = None

    if resume:
        try:
            checkpoint = load_checkpoint(ckpt_path)
        except FileNotFoundError as exc:
            raise ValueError(str(exc)) from exc
        validate_resume_params(checkpoint, resolved_params, step="pretree.trim")
        if checkpoint.status == "success":
            return _reconstruct_trim_result(
                checkpoint=checkpoint,
                params=resolved_params,
                global_warnings=global_warnings,
                scan_skipped=list(scan_skipped),
                wall_time=0.0,
            )
        to_run_ids, _ = plan_resume(checkpoint, verify_trim_outputs)
        if not to_run_ids:
            checkpoint.status = "success"
            checkpoint.completed_at = checkpoint.touch()
            save_checkpoint_atomic(checkpoint, ckpt_path)
            return _reconstruct_trim_result(
                checkpoint=checkpoint,
                params=resolved_params,
                global_warnings=global_warnings,
                scan_skipped=list(scan_skipped),
                wall_time=0.0,
            )
        found = [Path(t.input) for t in checkpoint.tasks if t.task_id in set(to_run_ids)]
    else:
        if not dry_run:
            if output_dir.exists() and any(output_dir.iterdir()):
                if not overwrite:
                    raise ValueError(
                        f"Output directory '{output_dir}' already exists and is non-empty. "
                        "Use --overwrite to replace it or --resume to continue."
                    )
                shutil.rmtree(output_dir)
            output_dir.mkdir(parents=True, exist_ok=True)

    if not dry_run:
        aa_out_dir.mkdir(parents=True, exist_ok=True)
        if nt_out_dir is not None:
            nt_out_dir.mkdir(parents=True, exist_ok=True)

    if not resume and not dry_run:
        checkpoint = build_initial_checkpoint(
            step="pretree.trim",
            command=(
                f"phyloai pretree trim --msa-dir {msa_dir} --tool {tool} "
                f"--seq-type {seq_type} --threads {threads}"
            ),
            params=resolved_params,
            inputs=found,
            output_for=lambda p: aa_out_dir / f"{p.stem}.fa",
            nt_output_for=(lambda p: None) if nt_out_dir is None
                           else (lambda p: nt_out_dir / f"{p.stem}.fa"),
        )
        save_checkpoint_atomic(checkpoint, ckpt_path)
        to_run_ids = [p.stem for p in found]

    _ckpt_write = checkpoint is not None and to_run_ids and not dry_run
    _to_run_set = set(to_run_ids) if to_run_ids else set()

    if _ckpt_write:
        for task_id in to_run_ids:
            mark_task(checkpoint, task_id, status="running", reason=None)
        save_checkpoint_atomic(checkpoint, ckpt_path)

    # --- Build worker args ---
    def _nt_path_for(msa_path: Path) -> Path | None:
        if nt_dir is None:
            return None
        candidates = [p for p in nt_dir.iterdir() if p.is_file() and p.stem == msa_path.stem]
        return candidates[0] if candidates else None

    # For BMGE mode4 downgrade: msa_path is the codon MSA from nt_dir, not msa_dir
    def _effective_msa(msa_path: Path) -> Path:
        if bmge_mode4_downgrade and nt_dir is not None:
            candidates = [p for p in nt_dir.iterdir() if p.is_file() and p.stem == msa_path.stem]
            return candidates[0] if candidates else msa_path
        return msa_path

    effective_seq_type_for_worker = "CODON" if bmge_mode4_downgrade else resolved_seq_type

    worker_args = [
        (
            _effective_msa(g), aa_out_dir, nt_out_dir,
            _nt_path_for(g) if not bmge_mode4_downgrade else None,
            tool, effective_seq_type_for_worker, trimal_method,
            clipkit_mode, bmge_matrix, bmge_entropy, extra_args,
            dry_run, trimal_exe, java_exe, bmge_jar, clipkit_exe,
        )
        for g in found
    ]

    file_results: list[dict[str, Any]] = []
    all_tool_results: list[dict[str, Any]] = []
    dry_run_cmds: list[str] = []
    skipped: list[dict[str, str]] = list(scan_skipped)
    _last_flush = time.monotonic()

    def _maybe_flush(*, force: bool = False) -> None:
        nonlocal _last_flush
        if not _ckpt_write:
            return
        now = time.monotonic()
        if force or (now - _last_flush) >= CHECKPOINT_FLUSH_INTERVAL:
            save_checkpoint_atomic(checkpoint, ckpt_path)
            _last_flush = now

    interrupted = False
    try:
        with ProcessPoolExecutor(max_workers=threads) as pool:
            futures = {pool.submit(_trim_one_worker, arg): arg[0] for arg in worker_args}
            for future in as_completed(futures):
                gene_path = futures[future]
                try:
                    result = future.result()
                except Exception as exc:
                    result = {"status": "skipped", "input": str(gene_path), "reason": str(exc)}

                task_id = Path(result.get("input", str(gene_path))).stem

                if result["status"] == "skipped":
                    skipped.append({"path": result["input"], "reason": result.get("reason", "unknown")})
                    all_tool_results.append(result)
                    if _ckpt_write and task_id in _to_run_set:
                        mark_task(checkpoint, task_id, status="failed", reason=result.get("reason"))
                    if progress_callback:
                        progress_callback(gene_path)
                    _maybe_flush()
                    continue

                if result["status"] == "dry_run":
                    all_tool_results.append(result)
                    cmd_val = result.get("cmd", "")
                    if isinstance(cmd_val, list):
                        dry_run_cmds.append(" ".join(cmd_val))
                    elif cmd_val:
                        dry_run_cmds.append(cmd_val)
                    if progress_callback:
                        progress_callback(gene_path)
                    continue

                # Validate outputs
                aa_out = Path(result["output_aa"])
                nt_out = Path(result["output_nt"]) if result.get("output_nt") else None
                if not verify_trim_outputs(aa_out, nt_out):
                    reason = "output validation failed (empty or unequal lengths)"
                    skipped.append({"path": result["input"], "reason": reason})
                    if _ckpt_write and task_id in _to_run_set:
                        mark_task(checkpoint, task_id, status="failed", reason=reason)
                    if progress_callback:
                        progress_callback(gene_path)
                    _maybe_flush()
                    continue

                file_results.append(result)
                all_tool_results.append(result)
                if _ckpt_write and task_id in _to_run_set:
                    mark_task(checkpoint, task_id, status="success", reason=None)
                if progress_callback:
                    progress_callback(gene_path)
                _maybe_flush()

    except KeyboardInterrupt:
        interrupted = True

    if _ckpt_write:
        if interrupted:
            checkpoint.status = "interrupted"
        save_checkpoint_atomic(checkpoint, ckpt_path, fsync=True)

    if interrupted:
        raise KeyboardInterrupt

    if not dry_run and all_tool_results:
        _write_trim_log(output_dir, all_tool_results)

    if not dry_run and not file_results:
        raise ValueError("No genes were trimmed: all input files failed or were skipped.")

    if checkpoint is not None and not dry_run:
        checkpoint.status = "success"
        checkpoint.completed_at = checkpoint.touch()
        save_checkpoint_atomic(checkpoint, ckpt_path, fsync=True)
        tool_versions = _detect_trim_tool_versions(
            tool=tool,
            trimal_path=trimal_path,
            bmge_path=bmge_path,
            clipkit_path=clipkit_path,
        )
        return _reconstruct_trim_result(
            checkpoint=checkpoint,
            params=resolved_params,
            global_warnings=global_warnings,
            scan_skipped=skipped,
            wall_time=time.monotonic() - run_start,
            tool_versions=tool_versions,
        )

    tool_versions = _detect_trim_tool_versions(
        tool=tool,
        trimal_path=trimal_path,
        bmge_path=bmge_path,
        clipkit_path=clipkit_path,
    ) if not dry_run else {}

    return _build_trim_payload(
        file_results=file_results,
        skipped=skipped,
        params=resolved_params,
        global_warnings=global_warnings,
        wall_time=time.monotonic() - run_start,
        tool_versions=tool_versions,
        dry_run_cmds=dry_run_cmds,
    )


def _build_trim_payload(
    *,
    file_results: list[dict[str, Any]],
    skipped: list[dict[str, str]],
    params: dict[str, Any],
    global_warnings: list[str],
    wall_time: float,
    tool_versions: dict[str, str],
    dry_run_cmds: list[str] | None = None,
) -> dict[str, Any]:
    lengths_before = [r["length_before"] for r in file_results if r.get("length_before")]
    lengths_after = [r["length_after"] for r in file_results if r.get("length_after")]

    def _stats(values: list[int]) -> dict[str, Any]:
        if not values:
            return {"mean": 0.0, "min": 0, "max": 0}
        return {"mean": round(sum(values) / len(values), 1), "min": min(values), "max": max(values)}

    cols_removed_pct: list[float] = []
    for r in file_results:
        lb = r.get("length_before", 0)
        la = r.get("length_after", 0)
        if lb and lb > 0:
            cols_removed_pct.append(round((lb - la) / lb * 100, 1))

    is_dry = bool(dry_run_cmds is not None and not file_results)
    return {
        "status": "success" if (file_results or is_dry) else "error",
        "command": (
            f"phyloai pretree trim --msa-dir {params['msa_dir']} "
            f"--tool {params['tool']} --seq-type {params['seq_type']}"
        ),
        "wall_time": wall_time,
        "tool_versions": tool_versions,
        "params": params,
        "key_results": {
            "total_genes": len(file_results) + len(skipped),
            "trimmed_genes": len(file_results),
            "skipped_genes": len(skipped),
            "skipped_reasons": _count_reasons(skipped),
            "length_before": _stats(lengths_before),
            "length_after": _stats(lengths_after),
            "columns_removed_pct": _stats([int(p) for p in cols_removed_pct]),
        },
        "error": None if (file_results or is_dry) else "No genes were trimmed.",
        "data": {
            "mode": _determine_mode(params),
            "skipped": [{"gene": Path(s["path"]).stem, "reason": s["reason"]} for s in skipped],
            "warnings": global_warnings,
            "dry_run_cmds": dry_run_cmds or [],
            "per_gene": [
                {
                    "gene": Path(r["input"]).stem,
                    "length_before": r.get("length_before", 0),
                    "length_after": r.get("length_after", 0),
                    "columns_removed": r.get("length_before", 0) - r.get("length_after", 0),
                    "outputs": [o for o in [r.get("output_aa"), r.get("output_nt")] if o],
                }
                for r in file_results
            ],
            "summary": {
                "n_input_files": len(file_results) + len(skipped),
                "n_trimmed": len(file_results),
                "n_skipped": len(skipped),
            },
        },
    }


def _reconstruct_trim_result(
    *,
    checkpoint: Checkpoint,
    params: dict[str, Any],
    global_warnings: list[str],
    scan_skipped: list[dict[str, str]],
    wall_time: float,
    tool_versions: dict[str, str],
) -> dict[str, Any]:
    file_results: list[dict[str, Any]] = []
    failed: list[dict[str, str]] = []

    for task in checkpoint.tasks:
        if task.status != "success":
            failed.append({"path": task.input, "reason": task.reason or task.status})
            continue
        aa_path = Path(task.outputs.get("aa", "")) if task.outputs.get("aa") else None
        if aa_path is None or not aa_path.exists():
            failed.append({"path": task.input, "reason": "missing AA output"})
            continue
        length_after = _read_msa_col_count(aa_path)
        file_results.append({
            "input": task.input,
            "output_aa": str(aa_path),
            "output_nt": task.outputs.get("nt"),
            "length_before": 0,  # not recoverable from checkpoint without re-reading original input
            "length_after": length_after,
        })

    all_skipped = list(scan_skipped) + failed
    return _build_trim_payload(
        file_results=file_results,
        skipped=all_skipped,
        params=params,
        global_warnings=global_warnings,
        wall_time=wall_time,
        tool_versions=tool_versions,
    )


def _count_reasons(skipped: list[dict[str, str]]) -> dict[str, int]:
    counts: dict[str, int] = {}
    for s in skipped:
        reason = s.get("reason", "unknown")
        counts[reason] = counts.get(reason, 0) + 1
    return counts


def _determine_mode(params: dict[str, Any]) -> str:
    seq_type = params.get("effective_seq_type", "AA")
    nt_dir = params.get("nt_dir")
    if seq_type == "CODON":
        return "CODON"
    if nt_dir:
        return "AA+NT"
    if seq_type == "NT":
        return "NT-only"
    return "AA-only"


def _write_trim_log(output_dir: Path, tool_results: list[dict[str, Any]]) -> None:
    log_path = output_dir / "trim.log"
    with open(log_path, "a") as fh:
        for res in tool_results:
            ts = datetime.datetime.now().isoformat(timespec="seconds")
            stderr = res.get("tool_stderr", "")
            fh.write(
                f"{'='*60}\n"
                f"timestamp:  {ts}\n"
                f"input:      {res.get('input')}\n"
                f"cmd:        {res.get('tool_cmd', res.get('cmd', ''))}\n"
                f"status:     {res.get('status')}\n"
                f"wall_time:  {res.get('wall_time', 0.0):.2f}s\n"
            )
            if stderr:
                fh.write(f"--- stderr ---\n{stderr}\n")


def render_trim_summary_table(summary: dict[str, Any]) -> "Table":
    from rich.table import Table as _Table
    table = _Table(title="pretree trim summary")
    table.add_column("Metric")
    table.add_column("Value")
    for key in ["n_input_files", "n_trimmed", "n_skipped"]:
        table.add_row(key, str(summary.get(key, "")))
    return table
```

- [ ] **Step 5.4: Run all trim tests**

```bash
python -m pytest tests/pretree/test_trim.py tests/core/test_sequence_normalization_codon.py -v
```

Expected: all tests PASS.

---

## Task 6: CLI registration

**Files:**
- Modify: `phyloai/cli/commands/pretree.py`

- [ ] **Step 6.1: Add trim import and subcommand**

In `phyloai/cli/commands/pretree.py`:

**1. Extend existing imports block** (the file already imports `json`, `Path`, `click`, `Console`, `Progress`; add the trim import alongside the existing pretree imports at the top):

```python
# Add this line after the existing pretree imports (align, convert, stats):
from phyloai.pretree.trim import run_trim, render_trim_summary_table, _scan_input as _trim_scan_input
```

**2. Update `_PretreeGroup.list_commands`** (line ~34):

```python
# Change:
return ["convert", "stats", "align"]
# To:
return ["convert", "stats", "align", "trim"]
```

No other imports are needed; `json`, `Path`, `click`, `Progress`, `console`, and `_fail` are already defined in the file and shared across all subcommands.

- [ ] **Step 6.2: Register the `trim` subcommand**

Append to `phyloai/cli/commands/pretree.py` (after the `align_command` function):

```python
@pretree.command(
    "trim",
    help=(
        "Trim multiple sequence alignments in batch using trimAl, BMGE, or ClipKIT.\n\n"
        "Supports four sequence modes:\n"
        "  AA-only:  --seq-type AA (or auto-detected)\n"
        "  NT-only:  --seq-type NT (or auto-detected)\n"
        "  CODON:    --seq-type CODON (codon-aligned NT MSA in --msa-dir; produces faa/ + fna/)\n"
        "  AA+NT:    --seq-type AA --nt-dir <path> (produces faa/ + fna/)\n\n"
        "For --nt-dir content by tool:\n"
        "  trimAl:  unaligned CDS sequences (same as pretree align --nt-dir)\n"
        "  ClipKIT: codon-aligned NT MSA\n"
        "  BMGE:    codon-aligned NT MSA (auto-switches to CODON mode; --seq-type AA is accepted)\n\n"
        "BMGE AA+NT: use --seq-type CODON with codon MSA in --msa-dir for explicit CODON trimming.\n\n"
        "--threads controls how many genes are trimmed in parallel."
    ),
)
@click.option("--msa-dir", type=click.Path(file_okay=False, path_type=Path),
              required=True,
              help="Input directory of aligned MSA files.")
@click.option("--output-dir", "-o", type=click.Path(file_okay=False, path_type=Path),
              default=Path("runs/pretree/trim"), show_default=True,
              help="Output directory; contains seqs/, trim.log, result.json.")
@click.option("--tool", type=click.Choice(["trimal", "bmge", "clipkit"]),
              default="trimal", show_default=True,
              help="Trimming tool to use.")
@click.option("--seq-type", type=click.Choice(["AA", "NT", "CODON", "auto"]),
              default="auto", show_default=True,
              help=(
                  "Molecule type. 'auto' detects AA vs NT from input; "
                  "CODON must be explicit (cannot be auto-detected)."
              ))
@click.option("--nt-dir", type=click.Path(file_okay=False, path_type=Path), default=None,
              help=(
                  "NT directory for AA+NT dual output. Content differs by tool: "
                  "trimAl expects unaligned CDS; ClipKIT/BMGE expect codon-aligned NT MSA."
              ))
@click.option("--trimal-method",
              type=click.Choice(["automated1", "gappyout", "strict", "strictplus"]),
              default="automated1", show_default=True,
              help=(
                  "trimAl automated trimming strategy. "
                  "automated1: conservative (recommended for most datasets). "
                  "gappyout: aggressive gap-based. "
                  "strict/strictplus: strictest heuristics."
              ))
@click.option("--bmge-matrix", type=str, default=None,
              help=(
                  "BMGE substitution matrix (-m). "
                  "AA/CODON: BLOSUM30 to BLOSUM95 (higher = stricter); default BLOSUM62. "
                  "NT: DNAPAMx:y (lower first number = stricter); default DNAPAM100:2. "
                  "Default is chosen automatically based on --seq-type."
              ))
@click.option("--bmge-entropy", type=float, default=0.5, show_default=True,
              help="BMGE entropy cutoff (-h). Lower = stricter; 0.2-0.4 is stringent.")
@click.option("--clipkit-method", type=str, default="smart-gap", show_default=True,
              help=(
                  "ClipKIT trimming mode (-m). Common values: "
                  "smart-gap (recommended, dynamic gap threshold), "
                  "kpi-smart-gap (keeps parsimony-informative sites + smart-gap), "
                  "kpic-smart-gap (keeps PI + constant sites + smart-gap), "
                  "gappy, kpi-gappy, kpic-gappy, kpi, kpic. "
                  "Full list: clipkit -h."
              ))
@click.option("--trimal-path", type=click.Path(dir_okay=False, path_type=Path), default=None,
              help="Explicit trimAl executable path; bundled version used when omitted.")
@click.option("--bmge-path", type=click.Path(dir_okay=False, path_type=Path), default=None,
              help="Explicit BMGE.jar path; bundled version used when omitted.")
@click.option("--clipkit-path", type=click.Path(dir_okay=False, path_type=Path), default=None,
              help="Explicit clipkit executable path; PATH lookup used when omitted.")
@click.option("--threads", "-t", type=int, default=4, show_default=True,
              help="Number of genes to trim in parallel.")
@click.option("--extra-args", type=str, default=None,
              help=(
                  "Extra arguments passed directly to the trimming tool. "
                  "Conflicts with internal arguments are resolved by extra-args winning."
              ))
@click.option("--resume", is_flag=True, default=False,
              help=(
                  "Resume from checkpoint.json in the output directory. "
                  "Requires the same parameters as the original run. "
                  "Mutually exclusive with --overwrite."
              ))
@click.option("--overwrite", is_flag=True, default=False,
              help="Delete and recreate a non-empty output directory before running.")
@click.option("--dry-run", is_flag=True, default=False,
              help="Print commands without executing; creates no files.")
@click.option("--quiet", "-q", is_flag=True, default=False,
              help="Suppress Rich terminal output except errors.")
def trim_command(
    msa_dir: Path,
    output_dir: Path,
    tool: str,
    seq_type: str,
    nt_dir: Path | None,
    trimal_method: str,
    bmge_matrix: str | None,
    bmge_entropy: float,
    clipkit_method: str,
    trimal_path: Path | None,
    bmge_path: Path | None,
    clipkit_path: Path | None,
    threads: int,
    extra_args: str | None,
    resume: bool,
    overwrite: bool,
    dry_run: bool,
    quiet: bool,
) -> None:
    if threads < 1:
        _fail("--threads must be at least 1.", 1)
    if not msa_dir.exists():
        _fail(f"--msa-dir '{msa_dir}' does not exist.", 1)
    if nt_dir is not None and not nt_dir.exists():
        _fail(f"--nt-dir '{nt_dir}' does not exist.", 1)
    if trimal_path is not None and not trimal_path.exists():
        _fail(f"--trimal-path '{trimal_path}' does not exist.", 1)
    if bmge_path is not None and not bmge_path.exists():
        _fail(f"--bmge-path '{bmge_path}' does not exist.", 1)
    if clipkit_path is not None and not clipkit_path.exists():
        _fail(f"--clipkit-path '{clipkit_path}' does not exist.", 1)

    payload: dict | None = None
    error_msg: str | None = None

    def _invoke(progress_callback=None):
        return run_trim(
            msa_dir=msa_dir,
            output_dir=output_dir,
            tool=tool,
            seq_type=seq_type,
            nt_dir=nt_dir,
            trimal_method=trimal_method,
            bmge_matrix=bmge_matrix,
            bmge_entropy=bmge_entropy,
            clipkit_mode=clipkit_method,   # CLI param is --clipkit-method; library param is clipkit_mode
            trimal_path=trimal_path,
            bmge_path=bmge_path,
            clipkit_path=clipkit_path,
            threads=threads,
            extra_args=extra_args,
            overwrite=overwrite,
            resume=resume,
            dry_run=dry_run,
            quiet=quiet,
            progress_callback=progress_callback,
        )

    if not quiet and not dry_run:
        found, _ = _trim_scan_input(msa_dir)
        with Progress(console=console, transient=True) as progress:
            task = progress.add_task("Trimming alignments", total=len(found))
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
        click.echo(f"Dry run: {payload['data']['summary']['n_input_files']} genes would be trimmed.")
        for cmd_str in payload["data"].get("dry_run_cmds", []):
            click.echo(cmd_str)
        return

    result_path = output_dir / "result.json"
    result_path.parent.mkdir(parents=True, exist_ok=True)
    with open(result_path, "w") as fh:
        json.dump(payload, fh, indent=2)

    if not quiet:
        console.print(render_trim_summary_table(payload["data"]["summary"]))
        seqs_path = output_dir / "seqs"
        click.echo(f"Trimmed alignments saved to {seqs_path}", err=True)
        click.echo(f"Results saved to {result_path}", err=True)
        for w in payload["data"].get("warnings", []):
            click.echo(f"Warning: {w}", err=True)
```

- [ ] **Step 6.3: Verify CLI is registered correctly**

```bash
python -m phyloai pretree --help
python -m phyloai pretree trim --help
```

Expected: `trim` appears in the pretree subcommand list; `trim --help` shows all options.

- [ ] **Step 6.4: Run existing tests to check no regressions**

```bash
python -m pytest tests/ -v --tb=short 2>&1 | tail -20
```

Expected: all previously passing tests still PASS.

---

## Task 7: Integration smoke tests (tools available)

**Files:**
- Modify: `tests/pretree/test_trim.py`

These tests are skipped when the tool is not available, so they run locally where tools are installed but do not block CI.

- [ ] **Step 7.1: Write integration tests**

Append to `tests/pretree/test_trim.py`:

```python
import shutil
import pytest


@pytest.mark.skipif(shutil.which("trimal") is None, reason="trimal not available")
def test_trim_trimal_aa_integration(tmp_path: Path) -> None:
    from phyloai.pretree.trim import run_trim

    msa_dir = tmp_path / "msa"
    msa_dir.mkdir()
    # 10 sequences, 10 AA positions, cols 2-5 are 90% gaps → will be trimmed by gappy
    aa_seqs = {
        "seq1": "MASTKLIVDE", "seq2": "M----LIVDE", "seq3": "M----LIVDE",
        "seq4": "M----LIVDE", "seq5": "M----LIVDE", "seq6": "M----LIVDE",
        "seq7": "M----LIVDE", "seq8": "M----LIVDE", "seq9": "M----LIVDE",
        "seq10": "M----LIVDE",
    }
    fasta = "\n".join(f">{n}\n{s}" for n, s in aa_seqs.items())
    (msa_dir / "gene1.faa").write_text(fasta)

    payload = run_trim(
        msa_dir=msa_dir,
        output_dir=tmp_path / "out",
        tool="trimal",
        seq_type="AA",
        trimal_method="gappyout",
        threads=1,
    )

    assert payload["status"] == "success"
    assert payload["data"]["summary"]["n_trimmed"] == 1
    out_file = tmp_path / "out" / "seqs" / "gene1.fa"
    assert out_file.exists()
    # Trimmed file should be shorter than 10 AA columns
    from Bio import SeqIO
    records = list(SeqIO.parse(str(out_file), "fasta"))
    assert len(records) == 10
    assert len(str(records[0].seq)) < 10


@pytest.mark.skipif(shutil.which("clipkit") is None, reason="clipkit not available")
def test_trim_clipkit_aa_integration(tmp_path: Path) -> None:
    from phyloai.pretree.trim import run_trim

    msa_dir = tmp_path / "msa"
    msa_dir.mkdir()
    aa_seqs = {
        "seq1": "MASTKLIVDE", "seq2": "M----LIVDE", "seq3": "M----LIVDE",
        "seq4": "M----LIVDE", "seq5": "M----LIVDE", "seq6": "M----LIVDE",
        "seq7": "M----LIVDE", "seq8": "M----LIVDE", "seq9": "M----LIVDE",
        "seq10": "M----LIVDE",
    }
    fasta = "\n".join(f">{n}\n{s}" for n, s in aa_seqs.items())
    (msa_dir / "gene1.faa").write_text(fasta)

    payload = run_trim(
        msa_dir=msa_dir,
        output_dir=tmp_path / "out",
        tool="clipkit",
        seq_type="AA",
        clipkit_mode="gappy",
        threads=1,
    )

    assert payload["status"] == "success"
    assert payload["data"]["summary"]["n_trimmed"] == 1


@pytest.mark.skipif(shutil.which("clipkit") is None, reason="clipkit not available")
def test_trim_clipkit_mode4_integration(tmp_path: Path) -> None:
    """ClipKIT Mode 4: trim AA with -l, project kept columns onto codon NT MSA."""
    from phyloai.pretree.trim import run_trim
    from Bio.SeqRecord import SeqRecord
    from Bio.Seq import Seq
    from Bio import SeqIO as _SeqIO

    msa_dir = tmp_path / "msa"
    nt_dir = tmp_path / "nt"
    msa_dir.mkdir()
    nt_dir.mkdir()

    codon_map = {
        "M": "ATG", "A": "GCT", "S": "TCT", "T": "ACT",
        "K": "AAA", "L": "CTT", "I": "ATT", "V": "GTT",
        "D": "GAT", "E": "GAA", "-": "---",
    }
    aa_seqs = {
        "seq1": "MASTKLIVDE", "seq2": "M----LIVDE", "seq3": "M----LIVDE",
        "seq4": "M----LIVDE", "seq5": "M----LIVDE", "seq6": "M----LIVDE",
        "seq7": "M----LIVDE", "seq8": "M----LIVDE", "seq9": "M----LIVDE",
        "seq10": "M----LIVDE",
    }
    aa_recs = [SeqRecord(Seq(s), id=n, description="") for n, s in aa_seqs.items()]
    codon_recs = [
        SeqRecord(Seq("".join(codon_map[aa] for aa in s)), id=n, description="")
        for n, s in aa_seqs.items()
    ]
    _SeqIO.write(aa_recs, str(msa_dir / "gene1.faa"), "fasta")
    _SeqIO.write(codon_recs, str(nt_dir / "gene1.fna"), "fasta")

    payload = run_trim(
        msa_dir=msa_dir,
        output_dir=tmp_path / "out",
        tool="clipkit",
        seq_type="AA",
        nt_dir=nt_dir,
        clipkit_mode="gappy",
        threads=1,
    )

    assert payload["status"] == "success"
    faa = tmp_path / "out" / "seqs" / "faa" / "gene1.fa"
    fna = tmp_path / "out" / "seqs" / "fna" / "gene1.fa"
    assert faa.exists()
    assert fna.exists()
    aa_out = list(_SeqIO.parse(str(faa), "fasta"))
    nt_out = list(_SeqIO.parse(str(fna), "fasta"))
    # NT should be 3x the length of AA
    assert len(str(nt_out[0].seq)) == len(str(aa_out[0].seq)) * 3
```

- [ ] **Step 7.2: Run integration tests**

```bash
python -m pytest tests/pretree/test_trim.py -v -k "integration"
```

Expected: tests run (or are skipped if tools not available); no errors.

- [ ] **Step 7.3: Run full test suite**

```bash
python -m pytest tests/ -v --tb=short 2>&1 | tail -30
```

Expected: all tests PASS.

- [ ] **Step 7.4: Commit when instructed**

When the user confirms all tasks are complete and tests pass, commit everything together:

```bash
git add phyloai/core/sequence_normalization.py \
        phyloai/pretree/trim.py \
        phyloai/cli/commands/pretree.py \
        tests/core/test_sequence_normalization_codon.py \
        tests/pretree/test_trim.py
git commit -m "feat: implement pretree trim command (trimAl/BMGE/ClipKIT, AA/NT/CODON/AA+NT modes)"
```

---

## Self-Review Checklist

**Spec coverage check:**
- [x] Section 1 (4 modes) → Tasks 3, 4, 5
- [x] Section 2 (all CLI params) → Task 6
- [x] Section 3 (file scanning, tool logic, log parsing, CODON validation) → Tasks 1, 2, 3, 4
- [x] Section 4 (param validation, per-gene errors, checkpoint) → Task 5
- [x] Section 5 (output dirs, log format, result.json) → Task 5
- [x] Section 6 (doctor integration) → No changes needed; TOOL_REGISTRY already has trimal/bmge/clipkit
- [x] Section 13 (key design decisions) → reflected in implementation choices

**Placeholder scan:** No TBDs, no "similar to" references. All code blocks are complete.

**Type consistency:** `_trim_one_worker` args tuple matches unpacking. `_make_success_result` used consistently in all three workers with `length_before` populated. `run_trim` passes `tool_versions` and `dry_run_cmds` through to `_build_trim_payload`. `_reconstruct_trim_result` passes `tool_versions` parameter.

**Commit policy:** No per-step commits. A single commit is made at the end of Task 7 when all tests pass, on explicit user instruction.
