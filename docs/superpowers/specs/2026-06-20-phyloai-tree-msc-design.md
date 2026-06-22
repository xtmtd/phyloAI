# PhyloAI Tree MSC (wASTRAL) Design Specification

**Date:** 2026-06-20
**Status:** Approved
**Parent spec:** `2026-06-07-phyloai-design.md`, `2026-06-17-phyloai-tree-design.md`

---

## 1. Purpose

`phyloai tree msc` performs multispecies coalescent species tree inference using wASTRAL (ASTER). It consumes gene trees and produces a species tree with local posterior probability branch support.

`msc` is a direct `click.Command` (not a Group) because there is only one MSC backend planned (wASTRAL/ASTER). No sub-backend selection is needed.

wASTRAL computation is one-shot — no `--resume` support. No checkpoint is produced.

---

## 2. CLI Surface

```bash
# Single gene tree file input
phyloai tree msc --tree gene_trees.trees -o runs/tree/msc

# Directory of gene tree files (auto-merge)
phyloai tree msc --tree-dir ./genetrees/

# Traditional unweighted Astral with exhaustive search
phyloai tree msc --tree-dir ./genetrees/ --mode 4 -R

# Full example: bootstrap input support, exhaustive search, all-three support output
phyloai tree msc --tree-dir ./genetrees/ \
    --mode 1 --boot 2 -R \
    --tree-boot-type bootstrap --tree-boot-min 10 --tree-boot-max 95 \
    -t 8 -o runs/tree/msc

# Override via --tool-args: custom search rounds
phyloai tree msc --tree input.trees --tool-args "-r 32 -s 32"

# Bayesian input support with custom range
phyloai tree msc --tree input.trees --tree-boot-type abayes \
    --tree-boot-min 0.5 --tree-boot-max 1.0
```

### Command Hierarchy

```
phyloai tree (click.Group)
├── ml (click.Group)
│   ├── fasttree
│   └── iqtree
├── bi (click.Group)
│   └── phylobayes
├── msc                # Direct command (click.Command), single backend: wASTRAL
└── concordance
```

---

## 3. Parameter Specification

### 3.1 Shared Parameters (from main design §9.2)

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--output-dir` / `-o` | Path | `runs/tree/msc` | Output directory. |
| `--threads` / `-t` | int ≥ 1 | 4 | Maps to wastral `-t` (thread count). |
| `--quiet` / `-q` | flag | False | Suppress terminal output except errors. |
| `--overwrite` | flag | False | Delete and recreate output directory. |
| `--tool-args` | str | None | Extra flags passed verbatim to wastral. Strategy-only; managed flags blocked (see §5). |
| `-h` / `--help` | — | — | Auto-handled by Click via root `CONTEXT_SETTINGS`. |

**No `--resume`** — wASTRAL computation is one-shot. Not applicable.

### 3.2 MSC-Specific Parameters

#### Input

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--tree` | Path | — | Single gene tree file (newick, one tree per line). Maps to wastral `-i`. Mutually exclusive with `--tree-dir`. |
| `--tree-dir` | Path | — | Directory of gene tree files. Merged into one input file. Mutually exclusive with `--tree`. |

`--tree` and `--tree-dir` are mutually exclusive. Providing both or neither exits with code 1.

#### Mode

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--mode` | int: 1/2/3/4 | 1 | wastral `--mode`. 1=hybrid (default), 2=branch support weighting, 3=branch length weighting, 4=traditional unweighted Astral. |

#### Branch Support (Output)

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--boot` | int: 0/1/2/3 | 1 | wastral `-u`. 0=species tree topology only (no support), 1=local posterior probability (default), 2=quartet support + local-PP for all three alternative topologies per branch, 3=same as 2 plus `freqQuad.csv` output. |

When `--boot` is overridden via `--tool-args` (contains `-u` or `--support`), phyloAI suppresses its own `-u` flag.

#### Extra Rounds

| Parameter | Short | Type | Default | Description |
|-----------|-------|------|---------|-------------|
| `--extra-rounds` | `-R` | flag | False | wastral `-R`: enables exhaustive search (`-r 16 -s 16`). |

When `-R` is set, phyloAI appends `-R` to the wastral command. If `--tool-args` already contains `-R`, `-r`, or `-s`, phyloAI suppresses its own `-R`.

#### Input Gene Tree Branch Support

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--tree-boot-type` | `auto` / `likelihood` / `abayes` / `bootstrap` | `auto` | Input gene tree branch support type. `auto`: wastral auto-detects (no flag passed). |
| `--tree-boot-min` | float | varies by type | wastral `-n` (minimum support threshold). Only valid when `--tree-boot-type` is non-auto. |
| `--tree-boot-max` | float | varies by type | wastral `-x` (maximum support value). Only valid when `--tree-boot-type` is non-auto. |

`--tree-boot-type` presets and their wastral flag mappings:

| Type | wastral flag(s) | default `--tree-boot-max` (`-x`) | default `--tree-boot-min` (`-n`) | `-d` (hardcoded) |
|------|----------------|--------------------------------|--------------------------------|---------------|
| `likelihood` | `--lrt` | 1 | 0 | 0 |
| `abayes` | `--bayes` | 1 | 0.333 | 0.333 |
| `bootstrap` | `--bootstrap` | 100 | 0 | 0 |
| `auto` | *(none)* | N/A | N/A | N/A |

`-d` is not a user-facing phyloAI parameter. Its default value per type is hardcoded and passed to wastral. Users can override `-d` via `--tool-args`.

When `--tree-boot-type` is `auto` (default): no `--lrt`/`--bayes`/`--bootstrap` flag is passed; wastral auto-detects. `--tree-boot-min` and `--tree-boot-max` are rejected with exit code 1 if provided.

When `--tree-boot-type` is explicitly set: phyloAI passes `--lrt`/`--bayes`/`--bootstrap` along with `-x VALUE -n VALUE -d VALUE`. If `--tool-args` already contains any of these flags, phyloAI suppresses its own corresponding output.

`--tree-boot-min` and `--tree-boot-max` are optional. When not provided, the type preset defaults are used (table above). When provided, they override the preset defaults.

#### Executable Path

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--wastral-path` | Path | None | Explicit wASTRAL executable path. None resolves via ToolEnv (see §11.4). |

#### Outgroup Rooting

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--outgroup` | str | None | Outgroup species name for rooting the output tree. Maps to wastral `--root`. |

---

## 4. wASTRAL Command Mapping

### 4.1 Mode → wastral flags

| `--mode` | wastral flag |
|----------|-------------|
| 1 (default) | `--mode 1` |
| 2 | `--mode 2` |
| 3 | `--mode 3` |
| 4 | `--mode 4` |

### 4.2 Boot → wastral flags

| `--boot` | wastral flag |
|----------|-------------|
| 0 | `-u 0` |
| 1 (default) | `-u 1` |
| 2 | `-u 2` |
| 3 | `-u 3` |

### 4.3 Extra Rounds → wastral flags

| `-R` / `--extra-rounds` | wastral flag |
|--------------------------|-------------|
| True | `-R` |
| False (default) | *(none)* |

### 4.4 Tree Boot Type → wastral flags

| `--tree-boot-type` | wastral flags (with defaults) |
|-------------------|-----------------------------|
| `auto` (default) | *(none)* |
| `likelihood` | `--lrt -x 1 -n 0 -d 0` |
| `abayes` | `--bayes -x 1 -n 0.333 -d 0.333` |
| `bootstrap` | `--bootstrap -x 100 -n 0 -d 0` |

When `--tree-boot-min` or `--tree-boot-max` is explicitly set, those values replace the defaults above.

### 4.5 Full Command Examples

```bash
# Default: hybrid mode, local-PP support, auto-detect input support
wastral -i input.trees -o output.tre -u 1 -t 4

# Unweighted Astral, exhaustive search, all-three resolutions
wastral -i input.trees -o output.tre --mode 4 -u 2 -t 8 -R

# Bootstrap input gene tree support, custom range
wastral -i input.trees -o output.tre -u 1 -t 4 --bootstrap -x 95 -n 10 -d 0

# Likelihood input support with custom range
wastral -i input.trees -o output.tre -u 1 -t 4 --lrt -x 1 -n 0.2 -d 0

# All-three resolutions + freqQuad.csv
wastral -i input.trees -o output.tre -u 3 -t 4
```

### 4.6 Command Assembly Order

```
wastral -i <input> -o <output_dir>/wastral.tre [--mode N] [-u N] [-t N] [--lrt/--bayes/--bootstrap -x MAX -n MIN -d D] [-R] [--tool-args tokens...]
```

1. wastral executable path
2. `-i <input>` (phyloAI-managed)
3. `-o <output_dir>/wastral.tre` (phyloAI-managed, blocked in tool-args)
4. Wastral managed parameters, each checked for `--tool-args` override:
   - `--mode N`
   - `-u N`
   - `-t N`
   - `--lrt`/`--bayes`/`--bootstrap` + `-x`/`-n`/`-d` (when applicable)
   - `-R`
5. `--tool-args` tokens appended last (via `shlex.split`)

---

## 5. `--tool-args` Two-Tier Model

### 5.1 Tier 1 — BLOCKED (Hard-Rejected)

PhyloAI always manages the tool's input and output file paths. The following are blocked in `--tool-args`:

| Flag | Reason |
|------|--------|
| `-i` | Input file, phyloAI-managed |
| `-o` | Output file, phyloAI-managed |
| Shell redirects (`>`, `<`, `|`) | phyloAI-managed I/O |

If `--tool-args` contains any blocked flag, exit code 1 with the blocked flag name.

### 5.2 Tier 2 — OVERRIDEABLE (Suppress-if-Present)

When `--tool-args` contains a flag that overlaps with a phyloAI-managed parameter, phyloAI suppresses its own version. Overlap detection is flag-name-only (no value parsing).

| phyloAI parameter | wastral flag(s) | Override behavior |
|-------------------|-----------------|-------------------|
| `--mode` | `--mode` | Suppress phyloAI `--mode` if present in tool-args |
| `--boot` | `-u`, `--support` | Suppress phyloAI `-u` if present in tool-args |
| `-R` / `--extra-rounds` | `-R`, `-r`, `-s` | Suppress phyloAI `-R` if any of `-R`/`-r`/`-s` present in tool-args |
| `--tree-boot-type` | `--lrt`, `--bayes`, `--bootstrap`, `-x`, `-n` | Suppress phyloAI boot-type flags if any overlap present in tool-args |
| `--threads` | `-t` | Suppress phyloAI `-t` if present in tool-args |

---

## 6. Input Validation

### 6.1 Mutual Exclusivity

`--tree` and `--tree-dir` are mutually exclusive. Providing both or neither exits with code 1.

### 6.2 `--tree` (Single File) Mode

- Validate file exists and is readable
- No format validation beyond existence (wastral will validate content)
- Pass directly to `wastral -i <file>`

### 6.3 `--tree-dir` (Directory) Mode

- Scan directory for gene tree files with extensions: `.nwk`, `.tre`, `.tree`, `.nw`, `.trees`, `.newick`
- Skip directories, empty files, files with unrecognized extensions
- Record unrecognized/empty files in `data.skipped`
- If zero valid files found: exit code 1
- Read each valid file, concatenate into a single file (one newick tree per line)
- Save merged file to `<output_dir>/merged.trees`
- Pass merged file to `wastral -i <output_dir>/merged.trees`
- If exactly 1 valid file: emit WARNING suggesting `--tree` may be more appropriate, continue

### 6.4 `--tree-boot-type` + min/max Validation

- `--tree-boot-min` and `--tree-boot-max` only valid when `--tree-boot-type` is non-auto
- If set with `--tree-boot-type auto`: exit code 1
- `--tree-boot-min` must be < `--tree-boot-max`: exit code 1 if violated

### 6.5 `--mode` Validation

- Must be one of 1, 2, 3, 4
- Click `IntRange(1, 4)` handles this at CLI layer

### 6.6 `--boot` Validation

- Must be one of 0, 1, 2, 3
- Click `IntRange(0, 3)` handles this at CLI layer

---

## 7. Output Directory Structure

```
runs/tree/msc/
├── result.json               # data.tool_stderr inlined (single pattern)
├── wastral.tre              # species tree output (wastral -o)
└── merged.trees             # merged gene tree input (--tree-dir mode only)
```

### 7.1 `--tree` mode

```
runs/tree/msc/
├── result.json               # data.tool_stderr inlined (single pattern)
└── wastral.tre
```

### 7.2 `--tree-dir` mode

```
runs/tree/msc/
├── result.json               # data.tool_stderr inlined (single pattern)
├── wastral.tre
└── merged.trees             # merged input for auditability
```

### 7.3 Output Directory Conflict Policy

- Default: if output directory exists and is non-empty, exit with code 1
- `--overwrite`: delete and recreate the output directory before running
- No `--resume`: wastral is one-shot computation

---

## 8. result.json Schema

```json
{
  "status": "success",
  "command": "phyloai tree msc --tree-dir ./genetrees/ --mode 1 --boot 2 -R ...",
  "wall_time": 45.2,
  "tool_versions": {"wastral": "1.25.4.8"},
  "params": {
    "tree": null,
    "tree_dir": "/path/to/genetrees",
    "mode": 1,
    "boot": 2,
    "extra_rounds": true,
    "tree_boot_type": "bootstrap",
    "tree_boot_min": 10,
    "tree_boot_max": 95,
    "outgroup": null,
    "output_dir": "runs/tree/msc",
    "threads": 8,
    "overwrite": false,
    "tool_args": null,
    "wastral_path": null
  },
  "key_results": {
    "mode": 1,
    "boot": 2,
    "extra_rounds": true,
    "tree_boot_type": "bootstrap",
    "outgroup": null,
    "n_input_trees": 1066,
    "input_mode": "--tree-dir"
  },
  "error": null,
  "data": {
    "input_mode": "--tree-dir",
    "input": {
      "path": "runs/tree/msc/merged.trees",
      "n_trees": 1066
    },
    "output_tree": "runs/tree/msc/wastral.tre",
    "tool_stderr": "# wastral diagnostic output (single pattern, JSON Output Standard Section 5.2)",
    "cmd": [
      "wastral", "-i", "runs/tree/msc/merged.trees",
      "-o", "runs/tree/msc/wastral.tre",
      "--mode", "1", "-u", "2", "-t", "8",
      "--bootstrap", "-x", "95", "-n", "10", "-d", "0",
      "-R"
    ],
    "skipped": [
      {"path": "genetrees/empty.tre", "reason": "empty file"},
      {"path": "genetrees/data.txt", "reason": "unrecognized file extension"}
    ],
    "warnings": [
      "genetrees directory contains non-newick files; skipped 2 file(s)"
    ]
  }
}
```

For `--tree` (single file) mode:

```json
{
  "data": {
    "input_mode": "--tree",
    "input": {
      "path": "/path/to/gene_trees.trees"
    },
    "output_tree": "runs/tree/msc/wastral.tre",
    "cmd": ["wastral", "-i", "/path/to/gene_trees.trees", "-o", "runs/tree/msc/wastral.tre", "-u", "1", "-t", "4"]
  }
}
```

### 8.1 Key Results Fields

| Field | Source | Description |
|-------|--------|-------------|
| `mode` | parameter | wastral `--mode` value |
| `boot` | parameter | wastral `-u` value |
| `extra_rounds` | parameter | `-R` flag value |
| `tree_boot_type` | parameter | input gene tree boot type |
| `n_input_trees` | computed | number of gene trees passed to wastral (`--tree-dir` mode) |
| `input_mode` | computed | `"--tree"` or `"--tree-dir"` |

---

## 9. Exit Codes

| Code | Meaning |
|------|---------|
| 0 | Success (species tree produced) |
| 1 | User input error (mutual exclusivity, invalid param, no valid inputs, output dir exists, blocked tool-args) |
| 2 | External tool execution failed (wastral non-zero exit) |
| 3 | Environment error (wastral not found or not executable) |

---

## 10. Warnings

| Condition | Behavior |
|-----------|----------|
| `--tree-dir` contains exactly 1 valid gene tree file | WARN: suggest `--tree` mode, continue |
| `--tree-dir` contains non-newick files | WARN per file: "unrecognized file extension", record in `data.skipped` |
| `--tree-dir` contains empty files | Skip silently, record in `data.skipped` |
| `--tree-dir` contains 0 valid files | Exit code 1 (not a warning) |

---

## 11. Logging

`msc` is a single-mode command. wASTRAL stderr is inlined in `result.json` as `data.tool_stderr` (single pattern, JSON Output Standard Section 5.2). No external log file is written.

wASTRAL stdout (the species tree in newick format) is captured to `wastral.tre` via `-o`.

---

## 12. Implementation Notes

### 12.1 Files to Create

| File | Purpose |
|------|---------|
| `phyloai/tree/msc.py` | `run_wastral()` library function |
| `docs/commands/tree-msc.md` | User-facing command documentation |

### 12.2 Files to Modify

| File | Change |
|------|--------|
| `phyloai/cli/commands/tree.py` | Add `msc` command (click.Command, not Group) to `tree` group |
| `phyloai/core/env.py` | Add `path_aliases: ["aster"]` to the `wastral` entry in TOOL_REGISTRY |
| `docs/superpowers/specs/2026-06-17-phyloai-tree-design.md` | Update msc section: direct command (not Group), `--tree`/`--tree-dir` input, output dir, CLI examples |
| `docs/superpowers/specs/2026-06-07-phyloai-design.md` | Update msc CLI examples and output dir examples (remove `wastral` subcommand and subdirectory) |

### 12.3 Key Patterns to Follow

- **CLI layer**: thin wrapper — validates params, resolves tool, delegates to library, writes result.json, renders Rich summary
- **Library layer**: `run_wastral()` accepts all params, validates preconditions, resolves input, builds command, runs subprocess, returns payload dict
- **No checkpoint**: wastral is one-shot; no `checkpoint.json`, no `--resume`
- **No batch/parallelism**: wastral handles its own internal multithreading via `-t`; no `ProcessPoolExecutor` needed at the phyloAI layer
- **Command override detection**: flag-name overlap check, same pattern as FastTree/IQ-TREE (see `_is_flag_overridden()` in `ml_iqtree.py`)
- **Log saving pattern**: wastral stderr captured via subprocess and inlined in `data.tool_stderr` (single pattern, JSON Output Standard Section 5.2)

### 12.4 wASTRAL Executable Resolution

Resolution order (same pattern as FastTree/IQ-TREE):

1. If `--wastral-path` is provided: validate the path exists and is executable. If not: exit code 1.
2. If `--wastral-path` is None: construct `ToolEnv(tool_paths={})` and call `require("wastral")`. This resolves via `shutil.which("wastral")` (and `shutil.which("aster")` if `path_aliases: ["aster"]` is added to the TOOL_REGISTRY entry in `core/env.py`).
3. If ToolEnv also fails to resolve: exit code 3 (environment error).
4. Version detection: attempt `wastral -v` and extract version string from output. Fallback: check `wastral -h` output for version information. The TOOL_REGISTRY entry in `core/env.py` already defines `version_args: [["-v"], ["-h"]]`.

### 12.5 Gene Tree Merging (`--tree-dir` mode)

```python
def _merge_gene_trees(tree_dir: Path, output_path: Path) -> int:
    """
    Scan tree_dir for newick files, merge into one file (one tree per line).
    Returns count of valid trees merged.
    """
    extensions = {".nwk", ".tre", ".tree", ".nw", ".trees", ".newick"}
    count = 0
    skipped = []
    with open(output_path, "w") as out:
        for f in sorted(tree_dir.iterdir()):
            if f.is_dir():
                continue
            if f.suffix.lower() not in extensions:
                skipped.append({"path": str(f), "reason": "unrecognized file extension"})
                continue
            content = f.read_text().strip()
            if not content:
                skipped.append({"path": str(f), "reason": "empty file"})
                continue
            # Write one tree per line (some files may have multiple trees)
            for line in content.splitlines():
                line = line.strip()
                if line:
                    out.write(line + "\n")
                    count += 1
    return count, skipped
```

### 12.6 Flag Override Detection for wastral

```python
_WASTRAL_MANAGED_FLAGS = {"-i", "-o"}  # Tier 1: BLOCKED

# Tier 2: OVERRIDEABLE flag groups
_WASTRAL_OVERRIDE_MAP = {
    "mode": {"--mode"},
    "boot": {"-u"},
    "extra_rounds": {"-R", "-r", "-s"},
    "tree_boot": {"--lrt", "--bayes", "--bootstrap", "-x", "-n"},
    "threads": {"-t"},
}
```

Managed flag checking follows the same pattern as IQ-TREE: tokenize `--tool-args` with `shlex.split`, then check token membership against blocked and override sets.

### 12.7 CLI Registration in `tree.py`

```python
# msc is a direct command, not a Group (single backend: wastral)
@tree.command("msc", cls=_GroupedHelpCommand)
@click.option("--tree", ...)
@click.option("--tree-dir", ...)
# ... (all options)
def msc_command(...):
    """Multispecies coalescent species tree inference with wASTRAL."""
    ...
```

The `tree` group in `_TreeGroup.list_commands()` should include `"msc"` in its return list.

---

## 13. Acceptance Criteria

Before merging, verify the following:

### 13.1 CLI Validation
- [ ] `--tree` and `--tree-dir` together → exit 1
- [ ] Neither `--tree` nor `--tree-dir` → exit 1
- [ ] Invalid `--mode` value → exit 1
- [ ] Invalid `--boot` value → exit 1
- [ ] `--tree-boot-min`/`--tree-boot-max` with `--tree-boot-type auto` → exit 1
- [ ] `--tree-boot-min` >= `--tree-boot-max` → exit 1

### 13.2 Input Scanning
- [ ] `--tree-dir` with 0 valid files → exit 1
- [ ] `--tree-dir` with exactly 1 valid file → WARNING, continue
- [ ] `--tree-dir` with non-newick files → WARNING per file, recorded in `data.skipped`
- [ ] `--tree` with nonexistent file → exit 1

### 13.3 `--tool-args` Blocking
- [ ] `--tool-args "-i other.trees"` → exit 1, blocked managed flag `-i`
- [ ] `--tool-args "-o output.tre"` → exit 1, blocked managed flag `-o`
- [ ] `--tool-args "-u 2"` → accepted (strategy parameter, overrides `--boot`)
- [ ] `--tool-args "-R"` → accepted (strategy parameter, overrides `--extra-rounds`)
- [ ] `--tool-args "--bootstrap -x 95 -n 10"` → accepted (strategy parameter, overrides `--tree-boot-type`)
- [ ] Valid strategy arg (e.g., `--tool-args "--root OUTGROUP"`) → appended to command

### 13.4 Tree Boot Type Presets
- [ ] `--tree-boot-type auto` (default) → no wastral boot flag
- [ ] `--tree-boot-type likelihood` → `--lrt -x 1 -n 0 -d 0`
- [ ] `--tree-boot-type abayes` → `--bayes -x 1 -n 0.333 -d 0.333`
- [ ] `--tree-boot-type bootstrap` → `--bootstrap -x 100 -n 0 -d 0`
- [ ] `--tree-boot-type bootstrap --tree-boot-min 10 --tree-boot-max 95` → `--bootstrap -x 95 -n 10 -d 0`

### 13.5 Mode Mapping
- [ ] `--mode 1` (default) → `--mode 1`
- [ ] `--mode 4` → `--mode 4` (unweighted Astral)

### 13.6 Boot Mapping
- [ ] `--boot 0` → `-u 0`
- [ ] `--boot 1` (default) → `-u 1`
- [ ] `--boot 2` → `-u 2`
- [ ] `--boot 3` → `-u 3`

### 13.7 Extra Rounds
- [ ] `--extra-rounds` / `-R` → wastral `-R`

### 13.8 Output
- [ ] `result.json` written with correct schema
- [ ] `wastral.tre` produced (valid newick species tree)
- [ ] `data.tool_stderr` populated with wastral stderr (single pattern)
- [ ] `merged.trees` produced (`--tree-dir` mode only)
- [ ] `tool_versions` populated with key `wastral`

### 13.9 Exit Codes
- [ ] Successful run → exit 0
- [ ] wastral non-zero exit → exit 2
- [ ] wastral not found → exit 3
- [ ] Invalid input → exit 1

---

## 14. Relationship to Other Modules

- **`tree ml`**: Supermatrix ML inference (FastTree, IQ-TREE). Consumes MSA, produces gene trees or species tree. Output gene trees from `tree ml fasttree --msa-dir` can be input to `tree msc --tree-dir`.
- **`tree bi`**: Bayesian inference (PhyloBayes). Alternative to ML and MSC.
- **`tree concordance`**: Computes gCF/sCF from a species tree and its supporting gene trees. Produced after `tree msc` if concordance analysis is desired.
- **`pretree`**: Produces the gene trees consumed by `tree msc`. The `phyloai run --mode coalescent` pipeline wires `pretree` → `tree ml fasttree --msa-dir` → `tree msc --tree-dir`.

---

## 15. Design Updates to Parent Specs

### 15.1 `docs/superpowers/specs/2026-06-17-phyloai-tree-design.md`

1. `msc` is a direct `click.Command` (not a `click.Group` with `wastral` subcommand)
2. CLI example: `phyloai tree msc --tree-dir ./genetrees/` (not `phyloai tree msc wastral --gene-trees ./genetrees/`)
3. Input parameters: `--tree` / `--tree-dir` (not `--gene-trees`)
4. Default output directory: `runs/tree/msc` (not `runs/tree/msc/wastral`)
5. No `--resume` support for msc

### 15.2 `docs/superpowers/specs/2026-06-07-phyloai-design.md`

1. CLI example at L120: change `phyloai tree msc wastral --gene-trees ./genetrees/` to `phyloai tree msc --tree-dir ./genetrees/`
2. Output directory tree at L223: change `msc/wastral/` to `msc/`
