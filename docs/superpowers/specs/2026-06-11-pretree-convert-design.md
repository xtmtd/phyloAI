# pretree convert - Design Specification

**Date:** 2026-06-11  
**Status:** Draft for user review  
**Parent spec:** `2026-06-07-phyloai-design.md`

---

## 1. Purpose

`phyloai pretree convert` is the workflow-specific normalization and format-conversion command for pre-tree sequence data. It prepares user-provided sequence or alignment files for the rest of the PhyloAI workflow by converting supported formats and normalizing sequence characters into a predictable, downstream-friendly representation.

It is not a general-purpose bioinformatics format converter. It supports only formats used by the PhyloAI workflow: FASTA, Phylip-relaxed, Phylip-PAML, and Nexus.

The recommended workflow is:

```bash
phyloai pretree convert --input ./raw --output-dir ./runs/pretree/convert --to fasta
phyloai pretree stats --seq-dir ./runs/pretree/convert
```

`convert` modifies output files only. It never edits input files in place.

---

## 2. Operating Model

### 2.1 Directory-first behavior

The primary input mode is directory-to-directory conversion:

```bash
phyloai pretree convert --input ./raw --output-dir ./runs/pretree/convert --to fasta
```

Directory mode scans one directory level only. It does not recurse into subdirectories. Each recognized sequence/alignment file is converted independently and written to the output directory with the same stem and a target-format suffix.

### 2.2 Single-file compatibility

`--input` may also point to a single file:

```bash
phyloai pretree convert --input ./gene001.phy --output-dir ./runs/pretree/convert --to fasta
```

Single-file mode uses the same format detection, normalization, and reporting rules as directory mode. The output file is written inside `--output-dir` using the input stem and target-format suffix.

### 2.3 Invalid input handling

Directory mode skips invalid entries and reports them in the summary. Skipped entries include subdirectories, empty files, non-sequence files, files with unsupported or undetectable formats, and files that fail during parsing or writing.

If at least one file is converted successfully, the command exits with code 0. If all candidate files fail or are skipped, the command exits with code 1 and reports the skip reasons.

---

## 3. Supported Formats

### 3.1 Format names

The public format names are:

| Format | Meaning | Typical suffix |
|--------|---------|----------------|
| `fasta` | FASTA sequence/alignment file | `.fa` |
| `phylip-relaxed` | PHYLIP-style alignment without the classic 10-character taxon-name limit | `.phy` |
| `phylip-paml` | PAML-compatible PHYLIP-like sequential alignment | `.paml.phy` |
| `nexus` | Nexus alignment file | `.nex` |

`phylip-relaxed` is the default Phylip interpretation in PhyloAI. The project should avoid using bare `phylip` in user-facing documentation except when referring to external software terminology.

### 3.1.1 Phylip output format

All phylip output (`phylip-relaxed` and `phylip-paml`) uses **sequential** format (one sequence per taxon, all on one line). This is the most common format expected by downstream phylogenetic software.

### 3.2 Phylip-PAML semantics

`phylip-paml` is not classic strict PHYLIP. PAML uses a PHYLIP-like native sequence format with these relevant rules:

- The first line contains the number of species and sequence length, optionally followed by PAML option characters.
- PAML defaults to sequential format when no `I` or `S` option is provided.
- PAML's species-name length is controlled by `LSPNAME` in program source; the documented default is 30 characters.
- Two consecutive spaces mark the end of the species name, so `pretree convert` must write at least two spaces between the name and the sequence.
- Species names should not contain two consecutive spaces. Spaces inside names should be normalized to `_` for safety.
- The special symbols `"`, `,`, `:`, `#`, `(`, `)`, `$`, and `=` should not appear in species names because PAML uses them for special purposes.

When writing `phylip-paml`, PhyloAI uses sequential output with one sequence record per taxon. The header should include the PAML `S` option, and each record should use the 30-character PAML name field followed by at least two spaces before the sequence.

---

## 4. Character Normalization

`pretree convert` always normalizes sequence characters by default. This is part of the command's purpose, not an optional post-processing step.

The normalization principle is: preserve broadly meaningful biological ambiguity where downstream tools commonly support it, but remove or standardize characters that commonly break phylogenetic software.

### 4.1 Shared classification with stats

`pretree stats` and `pretree convert` must use the same core character classification rules. `stats` reports what it sees; `convert` applies the corresponding normalization rules and reports what it changed.

The implementation should place shared character rules in `core`, outside CLI code, so both commands call the same definitions.

### 4.2 Nucleotide normalization

For `--seq-type NT` or auto-detected nucleotide data:

| Input character | Output behavior |
|-----------------|-----------------|
| `A C G T` | Keep, uppercase |
| `U` | Convert to `T` |
| `R Y S W K M B D H V N` | Keep, uppercase |
| `-` | Keep as alignment gap |
| `?` | Convert to `N` |
| `.` | Expand to the first sequence's same-position character when alignment context allows; otherwise convert to `N` |
| whitespace inside sequence lines | Remove |
| other letters or symbols | Convert to `N` and count as invalid replacements |

IUPAC nucleotide ambiguity codes are preserved rather than collapsed to `N`. This avoids losing information that may be useful for later codon-aware translation or filtering.

### 4.3 Amino-acid normalization

For `--seq-type AA` or auto-detected amino-acid data:

| Input character | Output behavior |
|-----------------|-----------------|
| `A R N D C Q E G H I L K M F P S T W Y V` | Keep, uppercase |
| `B Z J X U O` | Convert to `X` by default |
| `B Z J X U O` with `--aa-special keep` | Keep, uppercase |
| `-` | Keep as alignment gap |
| `?` | Convert to `X` |
| `*` | Convert to `X` and count separately as stop/termination replacements |
| `.` | Expand to the first sequence's same-position character when alignment context allows; otherwise convert to `X` |
| whitespace inside sequence lines | Remove |
| other letters or symbols | Convert to `X` and count as invalid replacements |

Defaulting `B/Z/J/X/U/O` to `X` favors compatibility with common protein phylogenetic software and standard 20-state amino-acid models. Advanced users who know their downstream tools can preserve those symbols with `--aa-special keep`.

### 4.4 Alignment-aware dot expansion

PAML treats `.` as the same character as the first sequence at that site. `pretree convert` should honor this behavior when the input can be read as an aligned matrix and the first sequence has a character at the same position.

If dot expansion is not possible because sequences are unaligned, lengths differ, or no same-position reference exists, dots are converted to `N` for NT or `X` for AA and counted in the report.

---

## 5. Taxon Name Normalization

Taxon names should be preserved when safe, but normalized when a target format has stricter constraints.

General name normalization:

- Strip leading and trailing whitespace.
- Replace internal whitespace runs with `_`.
- Preserve case unless a target writer requires otherwise.
- Ensure output names remain unique after normalization; if collisions occur, append deterministic suffixes such as `_2`, `_3` and report them.

Additional `phylip-paml` rules:

- Remove or replace PAML-problematic symbols: `"`, `,`, `:`, `#`, `(`, `)`, `$`, `=`.
- Truncate names longer than 30 characters by default.
- Preserve uniqueness after truncation by deterministic suffixing within the 30-character limit.
- Report every name replacement or truncation.
- Write the normalized PAML name in a 30-character field, then at least two spaces before the sequence.

---

## 6. CLI Parameters

| Parameter | Short | Type | Default | Notes |
|-----------|-------|------|---------|-------|
| `--input` | | Path | required | Input directory or single file |
| `--output-dir` | `-o` | Path | `runs/pretree/convert` | Output directory for converted files; future run allocation may resolve `runNNN` automatically |
| `--to` | | `fasta\|phylip-relaxed\|phylip-paml\|nexus` | `fasta` | Target output format |
| `--input-format` | | `auto\|fasta\|phylip-relaxed\|phylip-paml\|nexus` | `auto` | Override detection for all input files |
| `--seq-type` | | `AA\|NT\|auto` | `auto` | Override molecule-type detection |
| `--aa-special` | | `x\|keep` | `x` | Convert or preserve `B/Z/J/X/U/O` in AA data |
| `--threads` | `-t` | int | `4` | Directory mode parallelism |
| `--quiet` | `-q` | flag | `False` | Suppress terminal display except errors |
| `--overwrite` | | flag | `False` | Allow writing into an existing non-empty output directory |

`--output-dir` follows the global directory parameter convention and defaults under the standard run layout. The command does not write run logs and does not contribute to `report`.

---

## 7. Output Behavior

### 7.1 Output file naming

Output file stems match input file stems after path normalization. Target suffixes are:

| Target format | Suffix |
|---------------|--------|
| `fasta` | `.fa` |
| `phylip-relaxed` | `.phy` |
| `phylip-paml` | `.paml.phy` |
| `nexus` | `.nex` |

Converted sequence files are written to `seqs/` subdirectory inside `--output-dir`. The command always writes one output file per converted input file.

All PhyloAI-authored FASTA-family outputs from `convert` wrap sequence lines at 60 characters.

Example output structure:
```
runs/pretree/convert/
├── seqs/
│   ├── gene1.fa
│   ├── gene2.fa
│   └── ...
└── result.json
```

### 7.2 Output directory conflict policy

If `--output-dir` exists and is non-empty, the command exits with code 1 by default. With `--overwrite`, the output directory is deleted and recreated before conversion, matching the main design's output conflict policy.

### 7.3 Terminal output

Unless `--quiet` is set, directory mode displays a Rich progress indicator and a summary table containing:

- input entries scanned
- files converted
- entries skipped
- target format
- sequence-type summary
- total character replacements
- total taxon-name changes

Warnings are summarized after the table. Large per-file details belong in JSON output, not in terminal spam.

### 7.4 JSON result file

The full JSON result is written to `result.json` inside `--output-dir`. The terminal never prints raw JSON; it shows Rich tables only. This matches the behavior of all PhyloAI commands.

### 7.5 JSON result schema

`convert` is a utility command, so `key_results` is `{}`. Detailed conversion data appears under `data`:

```json
{
  "status": "success",
  "command": "phyloai pretree convert --input ./raw --output-dir ./runs/pretree/convert --to fasta",
  "wall_time": 0.0,
  "tool_versions": {},
  "params": {
    "input": "./raw",
    "output_dir": "./runs/pretree/convert",
    "to": "fasta",
    "input_format": "auto",
    "seq_type": "auto",
    "aa_special": "x",
    "threads": 4
  },
  "key_results": {},
  "error": null,
  "data": {
    "summary": {
      "n_input_entries": 100,
      "n_converted": 96,
      "n_skipped": 4,
      "target_format": "fasta",
      "seq_type_summary": "AA",
      "total_replacements": 17,
      "total_taxon_name_changes": 2
    },
    "files": [],
    "skipped": [],
    "warnings": []
  }
}
```

Each `files` entry records input path, output path, detected input format, target format, detected sequence type, replacement counts by category, taxon-name changes, and warnings for that file.

Each `skipped` entry records path and reason.

---

## 8. Integration with Core

`core/formats.py` should remain responsible for format naming, detection, reading, and writing. Its current enum semantics need correction in comments and docs:

- `AlignmentFormat.PHYLIP` maps to `phylip-relaxed`.
- `AlignmentFormat.PHYLIP_PAML` maps to `phylip-paml` conceptually, even if a custom writer uses Biopython's `phylip` dialect internally.
- `PHYLIP_PAML` must not be documented as classic strict PHYLIP.

Shared sequence normalization and character classification should live in a core-level helper used by both `pretree stats` and `pretree convert`. The helper should expose enough structured counts for reports without making the CLI parse warning text.

---

## 9. Documentation Requirements

Implementation must add or update:

- `docs/commands/pretree-convert.md` with the command-specific documentation sections required by the main design.
- `README.md` command index entry linking to `docs/commands/pretree-convert.md`.
- Any existing `pretree stats` documentation references that should now show the recommended `convert -> stats` order.

The top-level README should not contain the full `pretree convert` manual.

---

## 10. Testing Requirements

Tests should cover:

- directory input with mixed valid files, non-sequence files, empty files, and subdirectories
- single-file input writing to the default `runs/pretree/convert` output directory
- FASTA, Phylip-relaxed, Phylip-PAML, and Nexus conversion paths where library support permits
- NT normalization including `U -> T`, IUPAC ambiguity preservation, `? -> N`, and invalid symbols to `N`
- AA normalization including default `B/Z/J/X/U/O -> X`, `--aa-special keep`, `* -> X`, and invalid symbols to `X`
- PAML dot expansion against the first sequence when alignment context allows
- `phylip-paml` taxon-name cleanup, 30-character truncation, uniqueness preservation, and two-space separation
- output directory conflict behavior with and without `--overwrite`
- JSON schema shape and skipped-file reporting

---

## 11. Open Implementation Notes

Biopython may not provide an exact PAML writer that satisfies all name-separation and 30-character behavior requirements. If so, implementation should use a small custom writer for `phylip-paml` rather than weakening the format semantics.

The implementation plan should verify how existing Biopython readers handle PAML-style dot notation. If reader behavior loses dot context before normalization, the plan should read raw sequence text for that case or explicitly document the limitation before implementation.
