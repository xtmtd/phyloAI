# pretree convert Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build `phyloai pretree convert` as a directory-first sequence normalization and format-conversion command with single-file compatibility.

**Architecture:** Add a focused core normalization module for shared character and taxon-name rules, extend format handling for `phylip-paml`, then build `phyloai.pretree.convert` as the batch orchestration layer. Wire the CLI as a thin Click wrapper that delegates all conversion work to the pretree module and renders Rich/JSON output.

**Tech Stack:** Python 3.10+, Click, Rich, Biopython `SeqIO`/`AlignIO`, pytest, Click `CliRunner`.

---

## File Structure

- Create `phyloai/core/sequence_normalization.py`: shared sequence-type detection, character normalization, dot expansion, replacement counters, and taxon-name normalization.
- Modify `phyloai/core/formats.py`: add public `phylip-paml` enum semantics, suffix mapping, format normalization helpers, and a custom PAML writer if Biopython cannot guarantee two-space separation and 30-character behavior.
- Create `phyloai/pretree/convert.py`: conversion orchestration for one file or one directory, skipped-entry reporting, output naming, summary aggregation, and Rich table rendering.
- Modify `phyloai/pretree/stats.py`: import shared character constants/classification from `core.sequence_normalization` instead of keeping local divergent rules.
- Modify `phyloai/cli/commands/pretree.py`: register `convert` before `stats`, validate CLI parameters, call `phyloai.pretree.convert`.
- Create `tests/core/test_sequence_normalization.py`: unit tests for NT/AA normalization, dot expansion, and taxon-name cleanup.
- Modify `tests/core/test_formats.py`: cover `phylip-paml` detection/writing semantics.
- Create `tests/pretree/test_convert.py`: module-level conversion tests for file, directory, skipped entries, default run-layout output directory, and summary payload.
- Modify `tests/pretree/test_stats.py`: update character classification expectations that now come from shared core rules.
- Create `tests/cli/test_pretree_convert.py`: CLI tests for help, success, JSON output, overwrite behavior, and all-failed error behavior.
- Create `docs/commands/pretree-convert.md`: command documentation required by the main design.
- Modify `README.md`: keep a lightweight command index and link to command docs.

---

### Task 1: Shared Sequence Normalization Core

**Files:**
- Create: `phyloai/core/sequence_normalization.py`
- Test: `tests/core/test_sequence_normalization.py`

- [ ] **Step 1: Write failing tests for NT normalization**

Add `tests/core/test_sequence_normalization.py`:

```python
from __future__ import annotations

from phyloai.core.sequence_normalization import normalize_sequences


def test_normalize_nt_preserves_iupac_and_converts_u_question_and_invalid() -> None:
    result = normalize_sequences(["acguryswkmbdhvn?.!"], seq_type="NT")

    assert result.sequences == ["ACGTRYSWKMBDHVNNNN"]
    assert result.seq_type == "NT"
    assert result.replacements["u_to_t"] == 1
    assert result.replacements["question_to_missing"] == 1
    assert result.replacements["invalid_to_missing"] == 2
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/core/test_sequence_normalization.py::test_normalize_nt_preserves_iupac_and_converts_u_question_and_invalid -v`

Expected: FAIL with `ModuleNotFoundError: No module named 'phyloai.core.sequence_normalization'`.

- [ ] **Step 3: Implement NT normalization dataclass and function**

Create `phyloai/core/sequence_normalization.py` with:

```python
"""Shared sequence character normalization for PhyloAI commands."""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass, field


AA_STANDARD = set("ARNDCQEGHILKMFPSTWYV")
AA_SPECIAL = set("BZJXUO")
NT_STANDARD = set("ACGT")
NT_AMBIGUOUS = set("RYSWKMBDHVN")
NT_AUTO_CHARS = set("ACGTURYSWKMBDHVN")
GAP_CHARS = set("-?")
NT_MISSING = set("-?N")


@dataclass
class NormalizationResult:
    sequences: list[str]
    seq_type: str
    replacements: dict[str, int] = field(default_factory=dict)
    warnings: list[str] = field(default_factory=list)


def normalize_sequences(
    sequences: list[str],
    seq_type: str,
    aa_special: str = "x",
) -> NormalizationResult:
    counters: Counter[str] = Counter()
    normalized = [_normalize_one(sequence, seq_type, aa_special, counters) for sequence in sequences]
    return NormalizationResult(sequences=normalized, seq_type=seq_type, replacements=dict(counters), warnings=[])


def _normalize_one(sequence: str, seq_type: str, aa_special: str, counters: Counter[str]) -> str:
    output: list[str] = []
    for raw_char in sequence:
        if raw_char.isspace():
            counters["whitespace_removed"] += 1
            continue
        char = raw_char.upper()
        if seq_type == "NT":
            output.append(_normalize_nt_char(char, counters))
        else:
            output.append(_normalize_aa_char(char, aa_special, counters))
    return "".join(output)


def _normalize_nt_char(char: str, counters: Counter[str]) -> str:
    if char in NT_STANDARD or char in NT_AMBIGUOUS or char == "-":
        return char
    if char == "U":
        counters["u_to_t"] += 1
        return "T"
    if char == "?":
        counters["question_to_missing"] += 1
        return "N"
    counters["invalid_to_missing"] += 1
    return "N"


def _normalize_aa_char(char: str, aa_special: str, counters: Counter[str]) -> str:
    if char in AA_STANDARD or char == "-":
        return char
    if char in AA_SPECIAL:
        if aa_special == "keep":
            return char
        counters["aa_special_to_x"] += 1
        return "X"
    if char == "?":
        counters["question_to_missing"] += 1
        return "X"
    if char == "*":
        counters["stop_to_x"] += 1
        return "X"
    counters["invalid_to_missing"] += 1
    return "X"
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest tests/core/test_sequence_normalization.py::test_normalize_nt_preserves_iupac_and_converts_u_question_and_invalid -v`

Expected: PASS.

- [ ] **Step 5: Write failing tests for AA normalization and `--aa-special keep`**

Append to `tests/core/test_sequence_normalization.py`:

```python
def test_normalize_aa_converts_special_question_stop_and_invalid_to_x_by_default() -> None:
    result = normalize_sequences(["arndbzjxuop?*1"], seq_type="AA")

    assert result.sequences == ["ARNDXXXXXXPXXX"]
    assert result.replacements["aa_special_to_x"] == 6
    assert result.replacements["question_to_missing"] == 1
    assert result.replacements["stop_to_x"] == 1
    assert result.replacements["invalid_to_missing"] == 1


def test_normalize_aa_special_keep_preserves_bzjxuO() -> None:
    result = normalize_sequences(["BZJXUO"], seq_type="AA", aa_special="keep")

    assert result.sequences == ["BZJXUO"]
    assert result.replacements == {}
```

- [ ] **Step 6: Run AA tests to verify they pass after current implementation**

Run: `pytest tests/core/test_sequence_normalization.py -v`

Expected: PASS.

- [ ] **Step 7: Write failing tests for PAML dot expansion**

Append:

```python
from phyloai.core.sequence_normalization import expand_dots_from_first_sequence


def test_expand_dots_from_first_sequence_when_lengths_match() -> None:
    expanded, counts = expand_dots_from_first_sequence(["ACGT", "A..T", "...."], missing_char="N")

    assert expanded == ["ACGT", "ACGT", "ACGT"]
    assert counts == {"dot_expanded": 6}


def test_expand_dots_uses_missing_char_when_lengths_do_not_match() -> None:
    expanded, counts = expand_dots_from_first_sequence(["ACGT", "A.."], missing_char="N")

    assert expanded == ["ACGT", "ANN"]
    assert counts == {"dot_to_missing": 2}
```

- [ ] **Step 8: Implement dot expansion helper**

Add to `phyloai/core/sequence_normalization.py`:

```python
def expand_dots_from_first_sequence(sequences: list[str], missing_char: str) -> tuple[list[str], dict[str, int]]:
    counters: Counter[str] = Counter()
    if not sequences:
        return [], {}
    reference = sequences[0]
    same_length = all(len(sequence) == len(reference) for sequence in sequences)
    expanded: list[str] = []
    for sequence in sequences:
        chars: list[str] = []
        for index, char in enumerate(sequence):
            if char != ".":
                chars.append(char)
                continue
            if same_length and index < len(reference) and reference[index] != ".":
                chars.append(reference[index])
                counters["dot_expanded"] += 1
            else:
                chars.append(missing_char)
                counters["dot_to_missing"] += 1
        expanded.append("".join(chars))
    return expanded, dict(counters)
```

- [ ] **Step 9: Run normalization tests**

Run: `pytest tests/core/test_sequence_normalization.py -v`

Expected: PASS.

- [ ] **Step 10: Commit shared normalization core**

Run:

```bash
git add phyloai/core/sequence_normalization.py tests/core/test_sequence_normalization.py
git commit -m "feat: add sequence normalization core"
```

Expected: commit succeeds.

---

### Task 2: Shared Character Classification for stats

**Files:**
- Modify: `phyloai/core/sequence_normalization.py`
- Modify: `phyloai/pretree/stats.py`
- Modify: `tests/pretree/test_stats.py`

- [ ] **Step 1: Add shared detection/classification tests**

Append to `tests/core/test_sequence_normalization.py`:

```python
from phyloai.core.sequence_normalization import classify_char, detect_seq_type, normalize_pattern_char


def test_shared_detect_seq_type_nt_with_iupac() -> None:
    assert detect_seq_type(["ACGT", "AUGCRY", "NNNN"]) == "NT"


def test_shared_detect_seq_type_defaults_to_aa_when_x_seen() -> None:
    assert detect_seq_type(["ACGTX"]) == "AA"


def test_shared_classify_char_matches_stats_terms() -> None:
    assert classify_char("A", "AA") == "standard"
    assert classify_char("B", "AA") == "ambiguous"
    assert classify_char("-", "AA") == "gap"
    assert classify_char("N", "NT") == "gap"
    assert classify_char("R", "NT") == "ambiguous"


def test_shared_normalize_pattern_char_treats_question_mark_as_gap() -> None:
    assert normalize_pattern_char("?") == "-"
    assert normalize_pattern_char("a") == "A"
```

- [ ] **Step 2: Implement shared detection/classification helpers**

Add to `phyloai/core/sequence_normalization.py`:

```python
def detect_seq_type(sequences: list[str]) -> str:
    seq_type, _warnings = resolve_seq_type(sequences)
    return seq_type


def resolve_seq_type(sequences: list[str]) -> tuple[str, list[str]]:
    observed: set[str] = set()
    for sequence in sequences:
        observed.update(sequence.upper())
    observed.difference_update(GAP_CHARS)
    observed.discard(".")
    if observed & set("EFILPQWYZ"):
        return "AA", []
    if observed and observed <= NT_AUTO_CHARS:
        return "NT", []
    return "AA", ["[WARN] Cannot determine seq_type, defaulting to AA"]


def normalize_pattern_char(char: str) -> str:
    upper = char.upper()
    return "-" if upper == "?" else upper


def standard_chars(seq_type: str) -> set[str]:
    return NT_STANDARD if seq_type == "NT" else AA_STANDARD


def gap_chars(seq_type: str) -> set[str]:
    return NT_MISSING if seq_type == "NT" else GAP_CHARS


def classify_char(char: str, seq_type: str) -> str:
    upper = char.upper()
    if upper in standard_chars(seq_type):
        return "standard"
    if upper in gap_chars(seq_type):
        return "gap"
    return "ambiguous"
```

- [ ] **Step 3: Refactor stats imports without changing public behavior**

Modify `phyloai/pretree/stats.py`:

```python
from phyloai.core.sequence_normalization import (
    AA_STANDARD,
    GAP_CHARS,
    NT_AUTO_CHARS,
    NT_MISSING,
    NT_STANDARD,
    classify_char,
    detect_seq_type,
    gap_chars,
    normalize_pattern_char,
    resolve_seq_type,
    standard_chars,
)
```

Remove local definitions for `AA_STANDARD`, `AA_AMBIGUOUS`, `NT_STANDARD`, `NT_AMBIGUOUS`, `GAP_CHARS`, `NT_MISSING`, `NT_AUTO_CHARS`, `normalize_pattern_char`, `_standard_chars`, `_gap_chars`, `detect_seq_type`, `_resolve_seq_type`, and `classify_char`.

Replace remaining calls:

```python
gap_count = sum(sequence.count(char) for char in gap_chars(seq_type))
standard_count = sum(sequence.count(char) for char in standard_chars(seq_type))
detected_seq_type, seq_type_warnings = (seq_type, []) if seq_type else resolve_seq_type(sequences)
standard_codes = {ord(char) for char in standard_chars(seq_type)}
```

- [ ] **Step 4: Run stats and normalization tests**

Run: `pytest tests/core/test_sequence_normalization.py tests/pretree/test_stats.py -v`

Expected: PASS.

- [ ] **Step 5: Commit shared stats classification**

Run:

```bash
git add phyloai/core/sequence_normalization.py phyloai/pretree/stats.py tests/core/test_sequence_normalization.py tests/pretree/test_stats.py
git commit -m "refactor: share sequence character rules"
```

Expected: commit succeeds.

---

### Task 3: Format Semantics and PAML Writer

**Files:**
- Modify: `phyloai/core/formats.py`
- Modify: `tests/core/test_formats.py`

- [ ] **Step 1: Add failing tests for `phylip-paml` enum, compound suffix detection, and writer semantics**

Append to `tests/core/test_formats.py`:

```python
def test_phylip_paml_format_value_is_public_name() -> None:
    assert AlignmentFormat.PHYLIP.value == "phylip-relaxed"
    assert AlignmentFormat.PHYLIP_PAML.value == "phylip-paml"


def test_detect_phylip_paml_compound_suffix(tmp_path):
    paml = tmp_path / "gene.paml.phy"
    paml.write_text("2 4\ntaxon1  ACGT\ntaxon2  ACGA\n")

    assert FormatConverter().detect(paml) == AlignmentFormat.PHYLIP_PAML


def test_write_phylip_paml_uses_two_spaces_and_truncates_names(tmp_path):
    from Bio.Align import MultipleSeqAlignment
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from phyloai.core.formats import write_phylip_paml

    alignment = MultipleSeqAlignment([
        SeqRecord(Seq("ACGT"), id="Taxon name with spaces and very long suffix"),
        SeqRecord(Seq("ACGA"), id="Taxon:bad#chars"),
    ])
    out = tmp_path / "out.paml.phy"

    name_changes = write_phylip_paml(alignment, out)

    content = out.read_text().splitlines()
    assert content[0] == "2 4"
    assert "  " in content[1]
    assert len(content[1].split("  ", 1)[0]) <= 30
    assert ":" not in content[2].split("  ", 1)[0]
    assert "#" not in content[2].split("  ", 1)[0]
    assert len(name_changes) == 2
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/core/test_formats.py::test_phylip_paml_format_value_is_public_name tests/core/test_formats.py::test_detect_phylip_paml_compound_suffix tests/core/test_formats.py::test_write_phylip_paml_uses_two_spaces_and_truncates_names -v`

Expected: FAIL because `PHYLIP_PAML.value` is still `phylip`, compound `.paml.phy` suffix detection is missing, and `write_phylip_paml` is missing.

- [ ] **Step 3: Update enum and add custom writer**

Modify `phyloai/core/formats.py`:

```python
from Bio.SeqRecord import SeqRecord


class AlignmentFormat(str, Enum):
    FASTA = "fasta"
    PHYLIP = "phylip-relaxed"
    PHYLIP_PAML = "phylip-paml"
    NEXUS = "nexus"


_BIOPYTHON_FORMATS = {
    AlignmentFormat.FASTA: "fasta",
    AlignmentFormat.PHYLIP: "phylip-relaxed",
    AlignmentFormat.NEXUS: "nexus",
}


_COMPOUND_EXT_MAP: dict[str, AlignmentFormat] = {
    ".paml.phy": AlignmentFormat.PHYLIP_PAML,
}


def detect_alignment_format(path: Path, declared_format: Optional[AlignmentFormat] = None) -> AlignmentFormat:
    if declared_format is not None:
        return declared_format
    name = path.name.lower()
    for suffix, fmt in _COMPOUND_EXT_MAP.items():
        if name.endswith(suffix):
            return fmt
    suffix = path.suffix.lower()
    if suffix in _EXT_MAP:
        return _EXT_MAP[suffix]
    if path.exists():
        first = path.read_text(errors="ignore")[:200]
        if first.strip().startswith(">"):
            return AlignmentFormat.FASTA
        if first.strip().upper().startswith("#NEXUS"):
            return AlignmentFormat.NEXUS
        lines = [line for line in first.splitlines() if line.strip()]
        if lines and all(part.isdigit() for part in lines[0].strip().split()[:2]):
            return AlignmentFormat.PHYLIP
    raise ValueError(
        f"Cannot detect alignment format for '{path}'. "
        f"Supported extensions: {list(_EXT_MAP.keys()) + list(_COMPOUND_EXT_MAP.keys())}"
    )


def write_phylip_paml(alignment: MultipleSeqAlignment, dst: Path) -> list[dict[str, str]]:
    dst.parent.mkdir(parents=True, exist_ok=True)
    name_changes = _normalize_paml_names([record.id for record in alignment])
    length = alignment.get_alignment_length()
    with open(dst, "w") as fh:
        fh.write(f"{len(alignment)} {length}\n")
        for record, change in zip(alignment, name_changes):
            fh.write(f"{change['output']}  {str(record.seq)}\n")
    return [change for change in name_changes if change["input"] != change["output"]]


def _normalize_paml_names(names: list[str]) -> list[dict[str, str]]:
    used: set[str] = set()
    changes: list[dict[str, str]] = []
    for raw_name in names:
        safe = _safe_paml_name(raw_name)
        base = safe[:30] or "taxon"
        candidate = base
        suffix_number = 2
        while candidate in used:
            suffix = f"_{suffix_number}"
            candidate = f"{base[:30 - len(suffix)]}{suffix}"
            suffix_number += 1
        used.add(candidate)
        changes.append({"input": raw_name, "output": candidate})
    return changes


def _safe_paml_name(name: str) -> str:
    bad = {'"', ',', ':', '#', '(', ')', '$', '='}
    safe = "".join("_" if char.isspace() or char in bad else char for char in name.strip())
    while "__" in safe:
        safe = safe.replace("__", "_")
    return safe.strip("_")
```

Update `FormatConverter.detect()` to delegate to `detect_alignment_format(...)` so compound suffixes such as `.paml.phy` are checked before the normal `Path.suffix` fallback.

Update `FormatConverter.read()` and `FormatConverter.convert()` so public enum values are never passed directly to Biopython. `PHYLIP_PAML.value` is `phylip-paml` for the PhyloAI public API, but Biopython does not know that format name.

```python
class FormatConverter:
    def detect(
        self,
        path: Path,
        declared_format: Optional[AlignmentFormat] = None,
    ) -> AlignmentFormat:
        return detect_alignment_format(path, declared_format=declared_format)

    def convert(
        self,
        src: Path,
        dst: Path,
        target: AlignmentFormat,
        source_format: Optional[AlignmentFormat] = None,
        molecule_type: str = "protein",
    ) -> Path:
        src_fmt = source_format or self.detect(src)
        alignment = self.read(src, source_format=src_fmt)
        self.write_alignment(alignment, dst, target=target, molecule_type=molecule_type)
        return dst

    def read(
        self,
        path: Path,
        source_format: Optional[AlignmentFormat] = None,
    ) -> MultipleSeqAlignment:
        fmt = source_format or self.detect(path)
        if fmt == AlignmentFormat.PHYLIP_PAML:
            fmt = AlignmentFormat.PHYLIP
        return AlignIO.read(str(path), _BIOPYTHON_FORMATS[fmt])
```

`FormatConverter.read(...)` returns a `MultipleSeqAlignment`; this object is iterable over `SeqRecord` instances, so later conversion code may safely use `list(converter.read(...))` when it needs records.

- [ ] **Step 4: Run format tests**

Run: `pytest tests/core/test_formats.py -v`

Expected: PASS.

- [ ] **Step 5: Commit format semantics**

Run:

```bash
git add phyloai/core/formats.py tests/core/test_formats.py
git commit -m "feat: define phylip paml writer"
```

Expected: commit succeeds.

---

### Task 4: Pretree Convert Module

**Files:**
- Create: `phyloai/pretree/convert.py`
- Test: `tests/pretree/test_convert.py`

- [ ] **Step 1: Write failing single-file conversion test**

Create `tests/pretree/test_convert.py`:

```python
from __future__ import annotations

from pathlib import Path


def test_convert_single_file_defaults_to_fasta_and_normalizes_nt(tmp_path: Path) -> None:
    from phyloai.pretree.convert import convert_input

    src = tmp_path / "gene.fna"
    src.write_text(">tax one\nacgu?ry!\n")
    out_dir = tmp_path / "runs" / "pretree" / "convert"

    payload = convert_input(src, out_dir, target_format="fasta", seq_type="NT", threads=1, overwrite=False)

    out = out_dir / "gene.fa"
    assert out.exists()
    assert out.read_text() == ">tax_one\nACGTNRYN\n"
    assert payload["data"]["summary"]["n_converted"] == 1
    assert payload["data"]["summary"]["total_replacements"] == 3
    assert payload["data"]["files"][0]["output"] == str(out)
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/pretree/test_convert.py::test_convert_single_file_defaults_to_fasta_and_normalizes_nt -v`

Expected: FAIL with `ModuleNotFoundError: No module named 'phyloai.pretree.convert'`.

- [ ] **Step 3: Implement minimal `convert_input`, dot expansion, and FASTA writing**

Create `phyloai/pretree/convert.py`:

```python
"""Format conversion and sequence normalization for pretree workflows."""

from __future__ import annotations

import shutil
from pathlib import Path
from typing import Any

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from rich.table import Table

from phyloai.core.formats import AlignmentFormat, FormatConverter
from phyloai.core.sequence_normalization import detect_seq_type, expand_dots_from_first_sequence, normalize_sequences


TARGET_SUFFIX = {
    "fasta": ".fa",
    "phylip-relaxed": ".phy",
    "phylip-paml": ".paml.phy",
    "nexus": ".nex",
}


def convert_input(
    input_path: Path,
    output_dir: Path,
    target_format: str = "fasta",
    input_format: str | None = None,
    seq_type: str | None = None,
    aa_special: str = "x",
    threads: int = 4,
    overwrite: bool = False,
) -> dict[str, Any]:
    if output_dir.exists() and any(output_dir.iterdir()):
        if not overwrite:
            raise ValueError(f"Output directory '{output_dir}' already exists and is non-empty. Use --overwrite to replace it.")
        shutil.rmtree(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    entries = [input_path] if input_path.is_file() else sorted(input_path.iterdir(), key=lambda path: path.name)
    files: list[dict[str, Any]] = []
    skipped: list[dict[str, str]] = []
    for entry in entries:
        result = _convert_one(entry, output_dir, target_format, input_format, seq_type, aa_special)
        if "skipped" in result:
            skipped.append(result["skipped"])
        else:
            files.append(result)
    if not files:
        raise ValueError("All input entries failed or were skipped during conversion.")
    total_replacements = sum(sum(item["replacements"].values()) for item in files)
    payload = {
        "status": "success",
        "command": "phyloai pretree convert",
        "wall_time": 0.0,
        "tool_versions": {},
        "params": {
            "input": str(input_path),
            "output_dir": str(output_dir),
            "to": target_format,
            "input_format": input_format or "auto",
            "seq_type": seq_type or "auto",
            "aa_special": aa_special,
            "threads": threads,
        },
        "key_results": {},
        "error": None,
        "data": {
            "summary": {
                "n_input_entries": len(entries),
                "n_converted": len(files),
                "n_skipped": len(skipped),
                "target_format": target_format,
                "seq_type_summary": _mixed({item["seq_type"] for item in files}),
                "total_replacements": total_replacements,
                "total_taxon_name_changes": sum(item["taxon_name_changes"] for item in files),
            },
            "files": files,
            "skipped": skipped,
            "warnings": [warning for item in files for warning in item["warnings"]],
        },
    }
    return payload


def _convert_one(entry: Path, output_dir: Path, target_format: str, input_format: str | None, seq_type: str | None, aa_special: str) -> dict[str, Any]:
    if entry.is_dir():
        return {"skipped": {"path": str(entry), "reason": "directory"}}
    if not entry.is_file():
        return {"skipped": {"path": str(entry), "reason": "not a file"}}
    if entry.stat().st_size == 0:
        return {"skipped": {"path": str(entry), "reason": "empty file"}}
    converter = FormatConverter()
    try:
        fmt = converter.detect(entry, declared_format=_alignment_format(input_format))
        records = list(SeqIO.parse(str(entry), "fasta")) if fmt == AlignmentFormat.FASTA else list(converter.read(entry, source_format=fmt))
        if not records:
            return {"skipped": {"path": str(entry), "reason": "no sequences found"}}
        raw_sequences = [str(record.seq) for record in records]
        detected_seq_type = seq_type or detect_seq_type(raw_sequences)
        missing_char = "N" if detected_seq_type == "NT" else "X"
        expanded_sequences, dot_counts = expand_dots_from_first_sequence(raw_sequences, missing_char=missing_char)
        normalized = normalize_sequences(expanded_sequences, detected_seq_type, aa_special=aa_special)
        replacements = {**normalized.replacements}
        for key, value in dot_counts.items():
            replacements[key] = replacements.get(key, 0) + value
        normalized_records = [SeqRecord(Seq(sequence), id=_safe_record_id(record.id), description="") for record, sequence in zip(records, normalized.sequences)]
        out = output_dir / f"{entry.stem}{TARGET_SUFFIX[target_format]}"
        _write_records(normalized_records, out, target_format, detected_seq_type)
        return {
            "input": str(entry),
            "output": str(out),
            "input_format": fmt.value,
            "target_format": target_format,
            "seq_type": detected_seq_type,
            "replacements": replacements,
            "taxon_name_changes": sum(1 for original, record in zip(records, normalized_records) if original.id != record.id),
            "warnings": normalized.warnings,
        }
    except Exception as exc:
        return {"skipped": {"path": str(entry), "reason": str(exc)}}


def _write_records(records: list[SeqRecord], out: Path, target_format: str, seq_type: str) -> None:
    out.parent.mkdir(parents=True, exist_ok=True)
    if target_format == "fasta":
        SeqIO.write(records, str(out), "fasta")
        return
    from Bio.Align import MultipleSeqAlignment
    alignment = MultipleSeqAlignment(records)
    converter = FormatConverter()
    if target_format == "nexus":
        molecule_type = "DNA" if seq_type == "NT" else "protein"
        for record in alignment:
            record.annotations["molecule_type"] = molecule_type
    converter.write_alignment(alignment, out, target=_alignment_format(target_format), molecule_type="DNA" if seq_type == "NT" else "protein")


def _alignment_format(value: str | None) -> AlignmentFormat | None:
    if value in {None, "auto"}:
        return None
    for fmt in AlignmentFormat:
        if value == fmt.value:
            return fmt
    raise ValueError(f"Unsupported alignment format '{value}'.")


def _safe_record_id(name: str) -> str:
    return "_".join(name.strip().split()) or "taxon"


def _mixed(values: set[str]) -> str | None:
    if not values:
        return None
    if len(values) == 1:
        return next(iter(values))
    return "mixed"
```

- [ ] **Step 4: Add `FormatConverter.write_alignment` helper if missing**

Modify `phyloai/core/formats.py`:

```python
    def write_alignment(
        self,
        alignment: MultipleSeqAlignment,
        dst: Path,
        target: AlignmentFormat,
        molecule_type: str = "protein",
    ) -> list[dict[str, str]]:
        dst.parent.mkdir(parents=True, exist_ok=True)
        if target == AlignmentFormat.PHYLIP_PAML:
            return write_phylip_paml(alignment, dst)
        if target == AlignmentFormat.NEXUS:
            for record in alignment:
                if "molecule_type" not in record.annotations:
                    record.annotations["molecule_type"] = molecule_type
        with open(dst, "w") as fh:
            AlignIO.write(alignment, fh, _BIOPYTHON_FORMATS[target])
        return []
```

Update `convert()` to call `self.write_alignment(...)` after reading.

- [ ] **Step 5: Run single-file convert test**

Run: `pytest tests/pretree/test_convert.py::test_convert_single_file_defaults_to_fasta_and_normalizes_nt -v`

Expected: PASS.

- [ ] **Step 6: Add directory skipped-entry, dot expansion, and overwrite tests**

Append to `tests/pretree/test_convert.py`:

```python
def test_convert_directory_skips_non_sequence_empty_and_subdirectory(tmp_path: Path) -> None:
    from phyloai.pretree.convert import convert_input

    src_dir = tmp_path / "raw"
    src_dir.mkdir()
    (src_dir / "ok.fa").write_text(">a\nACGT\n")
    (src_dir / "notes.txt").write_text("not sequence")
    (src_dir / "empty.fa").write_text("")
    (src_dir / "nested").mkdir()
    out_dir = tmp_path / "out"

    payload = convert_input(src_dir, out_dir, target_format="fasta", seq_type="NT", threads=1, overwrite=False)

    assert payload["data"]["summary"]["n_converted"] == 1
    assert payload["data"]["summary"]["n_skipped"] == 3
    reasons = {item["reason"] for item in payload["data"]["skipped"]}
    assert "empty file" in reasons
    assert "directory" in reasons


def test_convert_expands_paml_dots_before_normalization(tmp_path: Path) -> None:
    from phyloai.pretree.convert import convert_input

    src = tmp_path / "dots.fa"
    src.write_text(">ref\nACGT\n>second\nA..T\n")
    out_dir = tmp_path / "out"

    payload = convert_input(src, out_dir, target_format="fasta", seq_type="NT", threads=1, overwrite=False)

    assert (out_dir / "dots.fa").read_text() == ">ref\nACGT\n>second\nACGT\n"
    assert payload["data"]["files"][0]["replacements"]["dot_expanded"] == 2


def test_convert_output_dir_conflict_requires_overwrite(tmp_path: Path) -> None:
    from phyloai.pretree.convert import convert_input

    src = tmp_path / "gene.fa"
    src.write_text(">a\nACGT\n")
    out_dir = tmp_path / "out"
    out_dir.mkdir()
    (out_dir / "old.fa").write_text("old")

    try:
        convert_input(src, out_dir, target_format="fasta", seq_type="NT", threads=1, overwrite=False)
    except ValueError as exc:
        assert "already exists and is non-empty" in str(exc)
    else:
        raise AssertionError("Expected output directory conflict")

    payload = convert_input(src, out_dir, target_format="fasta", seq_type="NT", threads=1, overwrite=True)
    assert payload["data"]["summary"]["n_converted"] == 1
    assert not (out_dir / "old.fa").exists()
```

- [ ] **Step 7: Run pretree convert module tests**

Run: `pytest tests/pretree/test_convert.py -v`

Expected: PASS.

- [ ] **Step 8: Commit pretree convert module**

Run:

```bash
git add phyloai/pretree/convert.py phyloai/core/formats.py tests/pretree/test_convert.py
git commit -m "feat: add pretree convert module"
```

Expected: commit succeeds.

---

### Task 5: CLI Command Wiring

**Files:**
- Modify: `phyloai/cli/commands/pretree.py`
- Create: `tests/cli/test_pretree_convert.py`

- [ ] **Step 1: Write failing CLI tests**

Create `tests/cli/test_pretree_convert.py`:

```python
from __future__ import annotations

import json
from pathlib import Path

from click.testing import CliRunner

from phyloai.cli.main import cli


def test_pretree_convert_help_is_registered_before_stats() -> None:
    result = CliRunner().invoke(cli, ["pretree", "--help"])

    assert result.exit_code == 0
    assert result.output.index("convert") < result.output.index("stats")


def test_cli_pretree_convert_single_file_json(tmp_path: Path) -> None:
    src = tmp_path / "gene.fa"
    src.write_text(">tax one\nacgu?\n")
    out_dir = tmp_path / "converted"

    result = CliRunner().invoke(
        cli,
        [
            "pretree",
            "convert",
            "--input",
            str(src),
            "--output-dir",
            str(out_dir),
            "--to",
            "fasta",
            "--seq-type",
            "NT",
            "--quiet",
        ],
    )

    assert result.exit_code == 0
    payload = json.loads(result.output)
    assert payload["status"] == "success"
    assert payload["data"]["summary"]["n_converted"] == 1
    assert (out_dir / "gene.fa").exists()


def test_cli_pretree_convert_all_failed_exits_one(tmp_path: Path) -> None:
    src_dir = tmp_path / "raw"
    src_dir.mkdir()
    (src_dir / "empty.fa").write_text("")

    result = CliRunner().invoke(cli, ["pretree", "convert", "--input", str(src_dir), "--output-dir", str(tmp_path / "out")])

    assert result.exit_code == 1
    assert "All input entries failed" in result.output
```

- [ ] **Step 2: Run CLI tests to verify they fail**

Run: `pytest tests/cli/test_pretree_convert.py -v`

Expected: FAIL because `convert` is not registered.

- [ ] **Step 3: Add Click command before `stats`**

Modify `phyloai/cli/commands/pretree.py` by importing convert helpers:

```python
from phyloai.pretree.convert import convert_input, render_convert_summary_table
```

Add the command above `@pretree.command("stats", ...)`:

```python
@pretree.command(
    "convert",
    help=(
        "Normalize and convert one sequence file or a directory of sequence files. "
        "Directory input is the primary mode; --input may also be a single file. "
        "Invalid directory entries are skipped and summarized."
    ),
)
@click.option("--input", "input_path", type=click.Path(path_type=Path), required=True, help="Input directory or single sequence/alignment file.")
@click.option("--output-dir", "output_dir", type=click.Path(file_okay=False, path_type=Path), default=Path("runs/pretree/convert"), show_default=True, help="Directory where converted files are written.")
@click.option("--to", "target_format", type=click.Choice(["fasta", "phylip-relaxed", "phylip-paml", "nexus"]), default="fasta", show_default=True, help="Target output format.")
@click.option("--input-format", type=click.Choice(["auto", "fasta", "phylip-relaxed", "phylip-paml", "nexus"]), default="auto", show_default=True, help="Override input format detection for all input files.")
@click.option("--seq-type", type=click.Choice(["AA", "NT", "auto"]), default="auto", show_default=True, help="Override sequence type detection.")
@click.option("--aa-special", type=click.Choice(["x", "keep"]), default="x", show_default=True, help="Convert B/Z/J/X/U/O to X, or preserve them with keep.")
@click.option("--threads", "threads", type=int, default=4, show_default=True, help="Directory mode worker count.")
@click.option("--quiet", "quiet", is_flag=True, default=False, help="Suppress Rich terminal output except errors.")
@click.option("--overwrite", "overwrite", is_flag=True, default=False, help="Delete and recreate a non-empty output directory before conversion.")
def convert_command(
    input_path: Path,
    output_dir: Path,
    target_format: str,
    input_format: str,
    seq_type: str,
    aa_special: str,
    threads: int,
    quiet: bool,
    overwrite: bool,
) -> None:
    if threads < 1:
        _fail("--threads must be at least 1.", 1)
    if not input_path.exists():
        _fail(f"Input path '{input_path}' does not exist.", 1)
    try:
        payload = convert_input(
            input_path,
            output_dir,
            target_format=target_format,
            input_format=None if input_format == "auto" else input_format,
            seq_type=None if seq_type == "auto" else seq_type,
            aa_special=aa_special,
            threads=threads,
            overwrite=overwrite,
        )
    except ValueError as exc:
        _fail(str(exc), 1)
    if not quiet:
        console.print(render_convert_summary_table(payload["data"]["summary"]))
        click.echo(f"Converted files saved to {output_dir / 'seqs'}", err=True)
        click.echo(f"Results saved to {output_dir / 'result.json'}", err=True)
```

- [ ] **Step 4: Add Rich summary renderer**

Add to `phyloai/pretree/convert.py`:

```python
def render_convert_summary_table(summary: dict[str, Any]) -> Table:
    table = Table(title="pretree convert summary")
    table.add_column("Metric")
    table.add_column("Value")
    for key in [
        "n_input_entries",
        "n_converted",
        "n_skipped",
        "target_format",
        "seq_type_summary",
        "total_replacements",
        "total_taxon_name_changes",
    ]:
        table.add_row(key, str(summary[key]))
    return table
```

- [ ] **Step 5: Run CLI tests**

Run: `pytest tests/cli/test_pretree_convert.py -v`

Expected: PASS.

- [ ] **Step 6: Run full relevant test suite**

Run: `pytest tests/core tests/pretree tests/cli -v`

Expected: PASS.

- [ ] **Step 7: Commit CLI wiring**

Run:

```bash
git add phyloai/cli/commands/pretree.py phyloai/pretree/convert.py tests/cli/test_pretree_convert.py
git commit -m "feat: add pretree convert cli"
```

Expected: commit succeeds.

---

### Task 6: Documentation

**Files:**
- Create: `docs/commands/pretree-convert.md`
- Modify: `README.md`

- [ ] **Step 1: Create command documentation**

Create `docs/commands/pretree-convert.md`:

```markdown
# phyloai pretree convert

## Purpose

`phyloai pretree convert` normalizes sequence characters and converts supported sequence/alignment formats for the PhyloAI workflow. It is intended for FASTA, Phylip-relaxed, Phylip-PAML, and Nexus files only; it is not a general-purpose format conversion tool.

## Usage

```bash
phyloai pretree convert --input ./raw --output-dir ./runs/pretree/convert --to fasta
```

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--input` | required | Input directory or single file |
| `--output-dir`, `-o` | `runs/pretree/convert` | Directory for converted files |
| `--to` | `fasta` | Target format: `fasta`, `phylip-relaxed`, `phylip-paml`, `nexus` |
| `--input-format` | `auto` | Override format detection |
| `--seq-type` | `auto` | Override molecule type detection |
| `--aa-special` | `x` | Convert `B/Z/J/X/U/O` to `X`, or preserve with `keep` |
| `--threads`, `-t` | `4` | Directory-mode worker count |
| `--quiet`, `-q` | false | Suppress Rich terminal output except errors |
| `--overwrite` | false | Delete and recreate a non-empty output directory |

## Inputs

`--input` may be a directory or a single file. Directory mode scans one level only and skips subdirectories, empty files, non-sequence files, and files that cannot be parsed.

## Outputs

Converted files are written to `--output-dir`. Target suffixes are `.fa`, `.phy`, `.paml.phy`, and `.nex`.

The JSON payload contains `summary`, `files`, `skipped`, and `warnings` under `data`. `key_results` is empty because `convert` is a utility command.

## Examples

```bash
phyloai pretree convert --input ./raw
phyloai pretree stats --seq-dir ./runs/pretree/convert/seqs
phyloai pretree convert --input ./gene.phy --output-dir ./converted --to fasta --seq-type NT
phyloai pretree convert --input ./aligned --to phylip-paml --overwrite
```

## Warnings and Errors

If some files are invalid, they are skipped and listed in the output. If all inputs fail or are skipped, the command exits with code 1. If the output directory exists and is non-empty, use `--overwrite` to replace it.

## Notes

Use `pretree convert` before `pretree stats` when raw input files may contain mixed formats or non-standard characters.
```

- [ ] **Step 2: Update README command index**

Modify `README.md` to keep only the short command list and links:

```markdown
## Commands

| Command | Purpose | Documentation |
|---------|---------|---------------|
| `phyloai doctor` | Inspect external tool availability | `docs/commands/doctor.md` |
| `phyloai pretree convert` | Normalize and convert sequence files | `docs/commands/pretree-convert.md` |
| `phyloai pretree stats` | Inspect one sequence file or summarize a directory | `docs/commands/pretree-stats.md` |
```

Keep the existing installation, quick start, and shell completion sections. If the README still contains the full `doctor` or `pretree stats` manuals, move that content into matching command docs in a separate documentation cleanup task rather than expanding this implementation task.

- [ ] **Step 3: Verify docs references**

Run: `pytest tests/cli/test_pretree_convert.py -v`

Expected: PASS. Documentation changes do not require additional automated tests in this repository.

- [ ] **Step 4: Commit documentation**

Run:

```bash
git add README.md docs/commands/pretree-convert.md
git commit -m "docs: document pretree convert"
```

Expected: commit succeeds.

---

### Task 7: Final Verification

**Files:**
- Verify: all changed files

- [ ] **Step 1: Run full test suite**

Run: `pytest -v`

Expected: all tests pass.

- [ ] **Step 2: Run CLI smoke test on fixture data**

Run:

```bash
phyloai pretree convert --input ref/phylogenomics_examples/test --output-dir /tmp/phyloai-convert-smoke --to fasta --overwrite --quiet
```

Expected: exit code 0 and `/tmp/phyloai-convert-smoke/result.json` contains `"status": "success"` and `"n_converted": 3`.

- [ ] **Step 3: Verify stats can read converted output**

Run:

```bash
phyloai pretree stats --seq-dir /tmp/phyloai-convert-smoke/seqs --quiet --output-dir /tmp/phyloai-convert-stats
```

Expected: exit code 0 and `/tmp/phyloai-convert-stats/result.json` contains JSON with `"status": "success"`.

- [ ] **Step 4: Inspect git status**

Run: `git status --short`

Expected: no uncommitted changes after all task commits.

---

## Self-Review

- Spec coverage: Tasks cover directory-first conversion, single-file compatibility, skipped invalid entries, default `runs/pretree/convert`, target formats, Phylip-PAML semantics, NT/AA normalization, shared stats rules, CLI, JSON output, docs, and verification.
- Placeholder scan: No placeholder implementation steps remain; each code task includes concrete tests, code, and commands.
- Type consistency: Public names are `normalize_sequences`, `expand_dots_from_first_sequence`, `convert_input`, `render_convert_summary_table`, `AlignmentFormat.PHYLIP_PAML`, `--output-dir`, `--to`, and `--aa-special` throughout the plan.
