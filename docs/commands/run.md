# phyloai run

## Purpose

One-click phylogenomics pipeline from raw sequence files to a species tree. Orchestrates all preprocessing and inference steps using sensible defaults. For fine-grained control over any individual step, use the constituent subcommands directly.

## Quick Start

```bash
# Typical usage: AA markers in a directory → species tree (supermatrix, normal speed)
phyloai run --seq-dir ./markers

# Fast exploratory run
phyloai run --seq-dir ./markers --speed fast

# Supertree (gene trees → species tree via wASTRAL)
phyloai run --seq-dir ./markers --mode supertree

# Resume an interrupted run
phyloai run --seq-dir ./markers --output-dir ./runs/run --resume

# Preview steps without executing
phyloai run --seq-dir ./markers --mode supertree --dry-run
```

## Input Requirements

`--seq-dir` must be a directory containing sequence files. Any format is accepted — all files are converted to normalized FASTA as the first step.

Supported input formats: `.fa`, `.fas`, `.fasta`, `.faa`, `.fna`, `.phy`, `.phylip`, `.nex`, `.nxs`, `.nexus`, `.aln`.

Sequence type (AA or NT) is auto-detected. For mixed directories, pre-convert with `phyloai pretree convert` directly.

## Pipeline Modes

### `--mode supermatrix` (default)

Aligns each gene, trims, optionally filters, concatenates into a single matrix, and infers a species tree from that matrix. Best for datasets where genes share a common history (no gene-tree/species-tree discordance expected).

```
convert → align → trim → [filter] → concat → species tree
```

### `--mode supertree`

Aligns each gene, trims, optionally filters, infers individual gene trees, and uses wASTRAL to build the species tree from the gene trees. Best when gene-tree heterogeneity is expected (ILS, HGT, duplication).

```
convert → align → trim → [filter] → gene trees → species tree
```

### Mode Comparison

| | Supermatrix | Supertree |
|---|---|---|
| Species tree method | IQ-TREE3 on concatenated matrix | wASTRAL from gene trees |
| Handles discordance | No | Yes |
| Computational cost | Lower (one tree) | Higher (N gene trees + coalescent) |
| Use when | Gene trees agree; genes share history | ILS, HGT, or gene duplication expected |
| Normal speed tree tool | IQ-TREE3 (unpartitioned) | IQ-TREE3 per gene + wASTRAL |
| Fast speed tree tool | FastTree on matrix | FastTree per gene + wASTRAL |

## Speed Modes

| Step | `--speed normal` | `--speed fast` |
|------|-----------------|---------------|
| Align | MAFFT `linsi` (highest accuracy) | MAFFT `auto` (heuristically chosen) |
| Trim | trimAl `-automated1` | trimAl `-automated1` |
| Filter | TAPER error-site masking | **skipped** |
| Gene trees | IQ-TREE3 with ModelFinder | FastTree |
| Species tree (supermatrix) | IQ-TREE3 (ModelFinder, unpartitioned) | FastTree |
| Species tree (supertree) | wASTRAL mode 1 | wASTRAL mode 1 |

`--speed fast` skips TAPER filtering entirely and trades IQ-TREE3 for FastTree. Use for quick exploratory analysis; use `--speed normal` for publication-quality results.

## Pipeline Steps in Detail

### Step 1: Convert (`1-convert/`)

Converts all input sequence files to normalized FASTA. Handles format detection, blank line removal, multi-line sequence joining, and illegal character stripping. Output files are written to `1-convert/seqs/*.fa`.

Uses `phyloai pretree convert`.

### Step 2: Align (`2-align/`)

Runs multiple sequence alignment on each FASTA file.

- `normal`: MAFFT with `--linsi` (iterative refinement, highest accuracy)
- `fast`: MAFFT with `--auto` (heuristic method selection)

Output: `2-align/seqs/*.fa`.

Uses `phyloai pretree align`.

### Step 3: Trim (`3-trim/`)

Trims poorly aligned columns using trimAl `-automated1` (heuristic selection of gap/similarity thresholds). Reduces noise from alignment uncertainty.

Output: `3-trim/seqs/*.fa`.

Uses `phyloai pretree trim`.

### Step 4: Filter (`4-filter/`, normal speed only)

Masks error-prone sites using TAPER (site-wise error probability estimation). Sites with high error probability are removed, improving downstream phylogenetic signal.

In `--speed fast`, this step is skipped and `4-filter/` is not created.

Output: `4-filter/seqs/*.fa`.

Uses `phyloai pretree filter`.

### Step 5: Concatenation or Gene Trees

**Supermatrix (`5-concat/`):** Concatenates all filtered/trimmed alignments into one supermatrix (`matrix.fa`). Outputs a partition file indicating gene boundaries.

Uses `phyloai pretree concat`.

**Supertree (`5-genetrees/`):** Infers a maximum-likelihood gene tree for each filtered/trimmed alignment.

- `normal`: IQ-TREE3 with ModelFinder per gene
- `fast`: FastTree per gene

Output: `5-genetrees/trees/*.treefile`.

Uses `phyloai tree ml iqtree --msa-dir` or `phyloai tree ml fasttree --msa-dir`.

### Step 6: Species Tree (`6-tree/`)

**Supermatrix:** Infers a species tree from the concatenated matrix.

- `normal`: IQ-TREE3 with automatic ModelFinder (unpartitioned; partitioned analyses require running `phyloai tree ml iqtree` directly)
- `fast`: FastTree

**Supertree:** Infers a species tree from gene trees using wASTRAL (mode 1, unrooted).

Uses `phyloai tree ml iqtree`, `phyloai tree ml fasttree`, or `phyloai tree msc`.

## Options

| Flag | Default | Description |
|------|---------|-------------|
| `--seq-dir PATH` | *(required)* | Input sequence directory. Any format; always converted first. |
| `--mode supermatrix\|supertree` | `supermatrix` | Pipeline mode. |
| `--speed normal\|fast` | `normal` | Speed/accuracy trade-off. `fast` skips TAPER and uses FastTree. |
| `-o, --output-dir PATH` | `runs/run` | Root output directory for all pipeline steps. |
| `-t, --threads INT` | `4` | Thread count passed to all steps. |
| `--resume` | off | Resume from `run_checkpoint.json`. |
| `--overwrite` | off | Delete and recreate output directory. Mutually exclusive with `--resume`. |
| `--dry-run` | off | Print step list without executing. |
| `-q, --quiet` | off | Suppress non-error output. |

## Output Structure

```
runs/run/
├── run_checkpoint.json        # pipeline-level checkpoint (resume support)
├── result.json                # overall pipeline result
├── 1-convert/
│   ├── result.json
│   └── seqs/                  # normalized FASTA files
├── 2-align/
│   ├── result.json
│   └── seqs/                  # aligned FASTA files
├── 3-trim/
│   ├── result.json
│   └── seqs/                  # trimmed alignments
├── 4-filter/                  # --speed normal only
│   ├── result.json
│   └── seqs/                  # TAPER-filtered alignments
├── 5-concat/                  # supermatrix mode
│   ├── result.json
│   └── matrix.fa              # concatenated supermatrix
├── 5-genetrees/               # supertree mode
│   ├── result.json
│   └── trees/                 # gene tree files
├── 6-tree/
│   └── result.json            # species tree result
```

Each step subdirectory contains its own `result.json` with detailed results (`tool_versions`, `key_results`, `data`).

## Checkpoint and Resume

`--resume` loads `run_checkpoint.json` from the output directory. The checkpoint records each step's status (`pending`, `running`, `success`, `failed`), output directory, and the parameter hash.

### Resume Behaviour

1. Validates that the checkpoint's `params_hash` matches current parameters. If parameters changed, exit 1 (use `--overwrite` for a clean run).
2. Skips steps marked `success` whose `result.json` contains `"status": "success"`.
3. Re-runs steps marked `running` or `interrupted` (restarts from scratch, not mid-step resume).
4. Steps that fail are marked `failed` and the pipeline halts.

Note: most individual steps have their own checkpoint/resume mechanisms. The pipeline-level resume is coarser — it re-runs interrupted steps entirely rather than resuming within a tool's native checkpoint.

### Checkpoint JSON Structure

```json
{
  "schema_version": 1,
  "step": "run",
  "command": "phyloai run --seq-dir ./markers --mode supermatrix",
  "status": "success",
  "params_hash": "abc123...",
  "params": { "seq_dir": "/abs/path/to/markers", "mode": "supermatrix", ... },
  "started_at": "2025-01-01T00:00:00Z",
  "updated_at": "2025-01-01T01:00:00Z",
  "completed_at": "2025-01-01T01:00:00Z",
  "steps": [
    { "name": "convert",      "status": "success", "output_dir": "/abs/path/1-convert" },
    { "name": "align",        "status": "success", "output_dir": "/abs/path/2-align" },
    ...
  ]
}
```

## Result JSON

The top-level `result.json` provides an overview of the entire pipeline:

```json
{
  "status": "success",
  "command": "phyloai run --seq-dir ./markers ...",
  "wall_time": 1234.5,
  "tool_versions": { "mafft": "7.520", "iqtree3": "2.3.6", ... },
  "params": { ... },
  "key_results": {
    "n_input_genes": 50,
    "n_genes_after_filter": 47,
    "final_tree": "/abs/path/6-tree/matrix.fa.treefile",
    "matrix_length": 35000,
    "matrix_taxa": 100
  },
  "data": {
    "mode": "supermatrix",
    "speed": "normal",
    "steps": [
      { "name": "convert", "status": "success", "output_dir": "...", "result_json": "..." },
      ...
    ]
  }
}
```

The `key_results` fields vary by mode — `matrix_length` and `matrix_taxa` appear only in supermatrix mode.

## Required Tools

| Tool | Used By | Check With |
|------|---------|------------|
| MAFFT | Step 2 (align) | `mafft --version` |
| trimAl | Step 3 (trim) | `trimal --version` |
| TAPER | Step 4 (filter, normal speed only) | `taper --version` or `taper -h` |
| IQ-TREE3 | Steps 5/6 (gene trees, species tree, normal speed) | `iqtree3 --version` or `iqtree2 --version` |
| FastTree | Steps 5/6 (gene trees, species tree, fast speed) | `fasttree -h` |
| wASTRAL | Step 6 (species tree, supertree mode) | `wastral -h` |

Run `phyloai doctor` to check your environment for all required tools.

## Exit Codes

| Code | Meaning |
|------|---------|
| `0` | Success — pipeline completed, final tree written to `6-tree/`. |
| `1` | Input error — missing `--seq-dir`, parameter mismatch on `--resume`, `--resume` + `--overwrite` together, non-empty output directory without `--overwrite`. |
| `2` | Step failure — a pipeline step (external tool error or data error) failed. Check the failing step's `result.json` and logs. |
| `3` | Environment error — a required tool is not installed. Run `phyloai doctor`. |

## Warnings / Errors

| Condition | Behaviour |
|-----------|-----------|
| `--seq-dir` does not exist or is empty | Exit 1 before any step starts. |
| `--resume` without existing checkpoint | Exit 1; use `--overwrite` for a clean run. |
| `--resume` + `--overwrite` together | Exit 1; mutually exclusive. |
| Parameter mismatch on `--resume` | Exit 1; checkpoint hash does not match current params. |
| Non-empty output dir without `--overwrite` | Exit 1; use `--overwrite` or `--resume`. |
| Step tool returns non-zero | Exit 2; step marked `failed` in checkpoint. |
| Required tool not found (`FileNotFoundError`) | Exit 3; e.g. `mafft`, `iqtree3`, `trimal`. |
| Final tree file not found after tree step | Exit 2; tree step claimed success but produced no output. |
| Intermediate result.json missing or not `"success"` on `--resume` | Step re-run (checkpoint state ignored). |

## Examples

```bash
# Default: supermatrix, normal speed
phyloai run --seq-dir ./markers

# Full command with explicit options
phyloai run --seq-dir ./markers --mode supermatrix --speed normal --threads 8

# Supertree with fast speed and 16 threads
phyloai run --seq-dir ./markers --mode supertree --speed fast --threads 16

# Custom output directory
phyloai run --seq-dir ./markers -o ./runs/my_analysis

# Resume a previously interrupted run
phyloai run --seq-dir ./markers --output-dir ./runs/run --resume

# Preview steps without running
phyloai run --seq-dir ./markers --mode supertree --dry-run

# Overwrite an existing run
phyloai run --seq-dir ./markers --overwrite
```

## Notes

- `phyloai run` uses each step's default parameters. For non-default settings (e.g. partitioned IQ-TREE, custom TAPER cutoff, specific MAFFT method), run steps individually via the constituent subcommands.
- In supermatrix normal mode, IQ-TREE3 runs automatic ModelFinder without a partition file. This is a first-pass unpartitioned result. Partitioned analyses require `phyloai tree ml iqtree` directly.
- In supertree normal mode, gene tree inference uses `phyloai tree ml iqtree --msa-dir` (batch IQ-TREE3 with ModelFinder per gene).
- `--threads` is passed to all steps. Each step may interpret it differently: MAFFT uses it as thread count, IQ-TREE batch mode uses it as parallel job count, trimAl uses it for parallelization.
- TAPER (`4-filter/`) runs only in `--speed normal`. In `--speed fast`, the directory is skipped and the post-trim alignments feed directly into concatenation or gene tree inference.
- Tool versions from every step are aggregated into the top-level `result.json` under `tool_versions`.
- The pipeline halts at the first failing step. Checkpoint state is saved before each step starts and after it completes, enabling `--resume` from the failed step.
- All convert/align/trim/filter steps use `overwrite=True` internally (forcing output regeneration within a pipeline run) since the pipeline manages directory state via its checkpoint.
- In `--speed fast`, both gene tree inference (supertree) and species tree inference (supermatrix) use FastTree, which is much faster but produces only approximate trees. FastTree uses a JTT+CAT model for AA and GTR+CAT for NT.
