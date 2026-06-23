# PhyloAI Posttree Topology Design Specification

**Date:** 2026-06-22
**Status:** Draft
**Parent spec:** `2026-06-07-phyloai-design.md`
**Reference:** IQ-TREE Advanced Tutorial, Tree topology tests

---

## 1. Purpose

`phyloai posttree topology` performs IQ-TREE tree topology tests for a single supermatrix and a user-provided set of candidate trees. It is a post-tree analysis command: by this stage, users usually have already inferred trees and selected or estimated a preferred model in `phyloai tree ml iqtree` or an equivalent workflow. Therefore, this command does not repeat the full `tree ml iqtree` model-selection interface.

The command focuses on one task:

```bash
iqtree3 -s <matrix> -z <candidate-trees> -n 0 -zb <replicates> -zw -au ...
```

It computes the standard IQ-TREE topology-test table for candidate topologies: log likelihoods, delta log likelihoods, bp-RELL, KH, SH, weighted KH, weighted SH, c-ELW, and AU where available in the `.iqtree` report.

Four-cluster Likelihood Mapping (FcLM) is not part of `posttree topology`. FcLM is a phylogenetic signal diagnostic and belongs in the future `posttree signal` module.

---

## 2. Design Principles

1. **Single supermatrix only.** No `--msa-dir` batch mode. Topology tests compare candidate trees against the same alignment and model.
2. **Fixed topology-test behavior by default.** PhyloAI generates `-n 0 -zb <replicates> -zw -au` so users get the full standard test set without choosing individual tests.
3. **Model expression, not model search.** Users provide a complete IQ-TREE model expression (`--model-expr`) or a previously optimized partition model (`--partitions`).
4. **Flexible advanced IQ-TREE usage.** `--tool-args` may override high-level defaults except matrix and candidate-tree inputs. This follows the project-wide strategy-only semantics: block only unsafe I/O overrides, not general tool strategy parameters.
5. **Reuse IQ-TREE infrastructure.** Reuse or extract helpers from `phyloai.tree.ml_iqtree` for IQ-TREE path resolution, version detection, input format checks, command assembly conventions, and result parsing where appropriate.

---

## 3. CLI Surface

```bash
phyloai posttree topology \
  --matrix raw.fa \
  --candidate-trees trees \
  --model-expr LG+F+R4 \
  --replicates 10000 \
  -t 20
```

### 3.1 Command Hierarchy

```
phyloai posttree
└── topology
```

`posttree` is a new Click group registered at the root CLI level. `topology` is a direct command, not a nested command group.

The parent design originally sketched this command with `--hypotheses`. This spec standardizes the name to `--candidate-trees` because IQ-TREE uses a tree-list file via `-z` and the command compares candidate topologies rather than naming statistical hypotheses directly.

---

## 4. Parameters

### 4.1 Required Input

| Flag | Type | Default | Description |
|------|------|---------|-------------|
| `--matrix` | Path | required | Single supermatrix alignment. Maps to IQ-TREE `-s`. |
| `--candidate-trees` | str | required | Candidate tree input. Accepts either one tree-list file path (one NEWICK tree per line), or multiple individual NEWICK tree file paths separated by commas (e.g. `h1.nwk,h2.nwk,h3.nwk`). Maps to IQ-TREE `-z` after optional merge. |
| `--input-format` | `auto\|fasta\|phylip-relaxed\|nexus` | `auto` | Optional PhyloAI-side matrix format hint for preflight validation. Not passed to IQ-TREE. |

Supported matrix extensions follow IQ-TREE-compatible formats already used by `tree ml iqtree`: `.fa`, `.fas`, `.fasta`, `.faa`, `.fna`, `.phy`, `.phylip`, `.nex`, `.nxs`, `.nexus`, `.aln`.

Candidate tree input has two forms:

1. A single `--candidate-trees candidate.trees` argument (no commas): treat it as the final IQ-TREE tree-list file and pass it directly to `-z`.
2. A comma-separated list `--candidate-trees h1.nwk,h2.nwk,...`: treat each segment as an individual NEWICK tree file, read them in order, write a merged `candidate.trees` under `--output-dir`, and pass that merged file to `-z`.

The merged `candidate.trees` file preserves user order exactly. Empty files, directories, and unreadable files are user input errors. The merged file contains one NEWICK tree per line.

`--seq-type` is intentionally not a high-level parameter in this command. The complete model expression already encodes the intended model family (for example, `LG` for amino acid and `GTR` for nucleotide), and IQ-TREE validates model/data compatibility. Advanced users who need an explicit IQ-TREE sequence type can pass it through `--tool-args`, for example `--tool-args "--seqtype AA"`.

### 4.2 Model Input

Provide one model source through either `--model-expr`, `--partitions`, or an explicit IQ-TREE model/partition flag in `--tool-args` (`-m` or `-p`). For ordinary use, prefer the high-level `--model-expr` or `--partitions` flags because they make `result.json` easier to interpret.

| Flag | Type | Default | Description |
|------|------|---------|-------------|
| `--model-expr` | str | None | Complete IQ-TREE `-m` model expression, e.g. `LG+F+R4`, `C20+F+R4`, `LG+C20+F+R4`, `custom.exchangeabilities+R4`. Mutually exclusive with high-level `--partitions`. |
| `--partitions` | Path | None | Previously optimized partition model, e.g. `*.best_model.nex`. Maps to IQ-TREE `-p`. Mutually exclusive with high-level `--model-expr`. |
| `--guide-tree` | Path | None | Guide tree for PMSF-style model expressions. Maps to IQ-TREE `-ft`. |

Rationale: topology tests are usually run after tree inference/model selection. Requiring users to restate a complete model expression is clearer than duplicating `tree ml iqtree` options such as `--model`, `--state-freq`, `--rate-heterogeneity`, `--modelfinder`, `--mset`, `--msub`, and heterogeneous workflow classification.

`--guide-tree` is retained because PMSF workflows commonly require `-ft` together with a complete model expression such as `LG+C20+F+R4`. Do not add `--fixed-tree`; the command does not expose IQ-TREE `-te` as a high-level parameter. Advanced users may pass IQ-TREE-specific strategy flags through `--tool-args`.

### 4.3 Topology Test Control

| Flag | Type | Default | Description |
|------|------|---------|-------------|
| `--replicates` | int >= 1000 | `10000` | RELL replicates. Maps to IQ-TREE `-zb`. IQ-TREE recommends at least 1000; PhyloAI recommends 10000. |

PhyloAI generates these defaults:

```bash
-n 0 -zb <replicates> -zw -au
```

This means: evaluate candidate trees without ML tree search (`-n 0`), run RELL resampling (`-zb`), include weighted KH/SH tests (`-zw`), and include AU test (`-au`). There is no `--tests`, `--weighted`, or `--au` option in the high-level interface.

### 4.4 Output and Execution

| Flag | Type | Default | Description |
|------|------|---------|-------------|
| `--prefix` | str | matrix stem | IQ-TREE output prefix. Maps to `--prefix`. |
| `-o`, `--output-dir` | Path | `runs/posttree/topology` | Output directory. |
| `-t`, `--threads` | int | `4` | Thread count. Maps to IQ-TREE `-T` unless overridden by `--tool-args`. |
| `--iqtree-path` | Path | None | Explicit path to `iqtree3`. |
| `--tool-args` | str | None | Additional IQ-TREE strategy parameters. See Section 6. |
| `--overwrite` | flag | False | Delete and recreate output directory. |
| `--resume` | flag | False | Allow IQ-TREE to continue/reuse its native output state in the existing output directory. No PhyloAI `checkpoint.json` is written. |
| `--dry-run` | flag | False | Print the IQ-TREE command and write no result files. |
| `-q`, `--quiet` | flag | False | Suppress terminal output except errors. |

---

## 5. Core Examples

### 5.1 Homogeneous Unpartitioned Model

```bash
phyloai posttree topology \
  --matrix raw.fa \
  --candidate-trees trees \
  --model-expr LG+F+R4 \
  --replicates 10000 \
  -t 20
```

Maps to:

```bash
iqtree3 -s raw.fa -m LG+F+R4 -z trees -n 0 -zb 10000 -zw -au -T 20
```

### 5.2 Heterogeneous Unpartitioned Model

```bash
phyloai posttree topology \
  --matrix raw.fa \
  --candidate-trees trees \
  --model-expr C20+F+R4 \
  -t 20
```

Maps to:

```bash
iqtree3 -s raw.fa -m C20+F+R4 -z trees -n 0 -zb 10000 -zw -au -T 20
```

### 5.3 PMSF Model

```bash
phyloai posttree topology \
  --matrix raw.fa \
  --candidate-trees trees \
  --model-expr LG+C20+F+R4 \
  --guide-tree guide.tree \
  -t 4
```

Maps to:

```bash
iqtree3 -s raw.fa -m LG+C20+F+R4 -ft guide.tree -z trees -n 0 -zb 10000 -zw -au -T 4
```

### 5.4 Previously Optimized Partition Model

```bash
phyloai posttree topology \
  --matrix raw.fa \
  --candidate-trees trees \
  --partitions raw.best_model.nex \
  -t 20
```

Maps to:

```bash
iqtree3 -s raw.fa -p raw.best_model.nex -z trees -n 0 -zb 10000 -zw -au -T 20
```

### 5.5 Multiple Individual Candidate Tree Files

```bash
phyloai posttree topology \
  --matrix raw.fa \
  --candidate-trees h1.nwk,h2.nwk,h3.nwk \
  --model-expr LG+F+R4 \
  -t 20
```

PhyloAI splits on commas, reads `h1.nwk`, `h2.nwk`, and `h3.nwk` in that order, merges them into `candidate.trees` under the output directory, then invokes:

```bash
iqtree3 -s raw.fa -m LG+F+R4 -z runs/posttree/topology/candidate.trees -n 0 -zb 10000 -zw -au -T 20
```

### 5.6 Custom Exchangeability Matrix and Site Frequencies

```bash
phyloai posttree topology \
  --matrix raw.fa \
  --candidate-trees trees \
  --model-expr custom.exchangeabilities+R4 \
  --tool-args "-fs custom.sitefreq" \
  -t 30
```

Maps to:

```bash
iqtree3 -s raw.fa -m custom.exchangeabilities+R4 -z trees -n 0 -zb 10000 -zw -au -T 30 -fs custom.sitefreq
```

---

## 6. `--tool-args` Semantics

The parent design uses a two-tier model for tool arguments. For `posttree topology`, the blocked set is intentionally minimal.

### 6.1 Blocked Flags

`--tool-args` must reject:

| Blocked token | Reason |
|---------------|--------|
| `-s` | Matrix input is managed by `--matrix`. |
| `-z` | Candidate trees input is managed by `--candidate-trees`. |
| `<`, `>`, `|` | Shell redirection/pipes are unsafe and break result capture. |

Only these flags are hard-blocked in the first implementation. `--prefix`, `-T`, `-m`, `-p`, `-ft`, `-n`, `-zb`, `-zw`, and `-au` are not blocked.

### 6.2 Suppress-if-Present Behavior

PhyloAI generates high-level defaults, but if the same IQ-TREE flag appears in `--tool-args`, PhyloAI suppresses its own version and lets `--tool-args` win.

Overrideable generated flags:

| PhyloAI option/default | IQ-TREE flag |
|------------------------|--------------|
| `--model-expr` | `-m` |
| `--partitions` | `-p` |
| `--guide-tree` | `-ft` |
| `--replicates` | `-zb` |
| fixed no-search default | `-n` |
| fixed weighted-test default | `-zw` |
| fixed AU-test default | `-au` |
| `--threads` | `-T` |
| `--prefix` | `--prefix` |

For ordinary users, the default command remains the complete recommended topology-test invocation. Advanced users can intentionally change IQ-TREE behavior with `--tool-args`, and `result.json` records the exact final command.

### 6.3 Assembly Order

1. Add executable, `-s <matrix>`, and `-z <candidate-trees>`.
2. Add `--prefix` unless `--tool-args` already contains `--prefix`.
3. Add model source: `-m <model-expr>` or `-p <partitions>` unless overridden.
4. Add `-ft <guide-tree>` unless overridden.
5. Add topology-test defaults `-n 0 -zb <replicates> -zw -au`, each suppressible by matching `--tool-args` tokens.
6. Add `-T <threads>` unless overridden.
7. Append tokenized `--tool-args` last.

---

## 7. Validation Rules

1. `--matrix` is required and must exist.
2. At least one `--candidate-trees` value is required.
3. A `--candidate-trees` value with no commas is treated as a tree-list file and must exist, be a regular non-empty file, and be readable.
4. A `--candidate-trees` value containing commas is treated as ordered individual tree files; each segment must exist, be a regular non-empty file, and be readable. PhyloAI writes `candidate.trees` under `--output-dir` before invoking IQ-TREE.
5. `--input-format` must be one of `auto`, `fasta`, `phylip-relaxed`, or `nexus`.
6. Exactly one of `--model-expr` and `--partitions` is required unless `--tool-args` contains `-m` or `-p`. If `--tool-args` provides the model source, both high-level model-source flags may be omitted.
7. `--model-expr` and `--partitions` are mutually exclusive unless one is intentionally overridden through `--tool-args`. First implementation should reject both high-level flags together.
8. `--guide-tree`, if provided, must exist.
9. `--partitions`, if provided, must exist.
10. `--replicates` must be >= 1000. No hard maximum is imposed; help text should warn that very large values can make RELL resampling slow.
11. `--threads` must be >= 1.
12. `--overwrite` and `--resume` are mutually exclusive.
13. Existing non-empty output directory without `--overwrite` or `--resume` exits with code 1.
14. `--tool-args` rejects only blocked flags and unsafe shell tokens listed in Section 6.1.

---

## 8. Output Directory

Default output directory:

```text
runs/posttree/topology/
├── result.json
├── <prefix>.iqtree
├── <prefix>.log
├── <prefix>.treels.trees       # IQ-TREE optimized candidate trees, exact suffix may vary by IQ-TREE version
├── <prefix>.ckp.gz             # if IQ-TREE writes checkpoint data
└── candidate.trees             # only when multiple --candidate-trees files were merged
```

All IQ-TREE-native output files are preserved in the output directory. `result.json` is the authoritative structured output for PhyloAI consumers. IQ-TREE diagnostic output follows the single-command JSON pattern: merged stdout/stderr is stored in `data.tool_stderr`, and IQ-TREE-native `.log` / `.iqtree` files are referenced via `data.tool_log` and `data.log_iqtree`.

The exact optimized candidate-tree suffix can differ between IQ-TREE versions (`*.treels.trees`, `*.trees`, or a prefix-derived equivalent). The implementation should discover candidate-tree output files by prefix and known suffix patterns, not by assuming only one exact filename.

---

## 9. `result.json` Schema

Example:

```json
{
  "status": "success",
  "command": "phyloai posttree topology --matrix raw.fa --candidate-trees trees --model-expr LG+F+R4 --replicates 10000 -t 20 -o runs/posttree/topology",
  "wall_time": 42.5,
  "tool_versions": {"iqtree3": "3.1.3"},
  "params": {
    "matrix": "/abs/raw.fa",
    "candidate_trees": ["/abs/trees"],
    "candidate_trees_mode": "tree-list",
    "candidate_trees_effective": "/abs/trees",
    "input_format": "auto",
    "model_expr": "LG+F+R4",
    "partitions": null,
    "guide_tree": null,
    "replicates": 10000,
    "prefix": "raw",
    "output_dir": "runs/posttree/topology",
    "threads": 20,
    "iqtree_path": null,
    "tool_args": null,
    "overwrite": false,
    "resume": false,
    "dry_run": false,
    "quiet": false
  },
  "key_results": {
    "n_candidate_trees": 4,
    "best_tree_id": 1,
    "n_rejected_au_0_05": 1,
    "replicates": 10000,
    "model_source": "model-expr"
  },
  "error": null,
  "data": {
    "cmd": ["iqtree3", "-s", "/abs/raw.fa", "-z", "/abs/trees", "--prefix", "raw", "-m", "LG+F+R4", "-n", "0", "-zb", "10000", "-zw", "-au", "-T", "20"],
    "tool_stderr": "# merged stdout/stderr from IQ-TREE",
    "log_iqtree": "raw.iqtree",
    "tool_log": "raw.log",
    "optimized_trees": "raw.treels.trees",
    "merged_candidate_trees": null,
    "tests": [
      {
        "tree_id": 1,
        "log_likelihood": -21152.617,
        "delta_likelihood": 0.0,
        "bp_rell": 0.711,
        "bp_rell_sign": "+",
        "p_kh": 0.74,
        "p_kh_sign": "+",
        "p_sh": 1.0,
        "p_sh_sign": "+",
        "c_elw": 0.6954,
        "c_elw_sign": "+",
        "p_au": 0.7939,
        "p_au_sign": "+"
      }
    ],
    "warnings": []
  }
}
```

`data.tests[]` is parsed from the `USER TREES` section of the `.iqtree` report. Parsing must be best-effort across IQ-TREE versions. If a column is absent, store `null`. If parsing fails entirely but IQ-TREE exits successfully and the `.iqtree` file exists, overall status remains `success` with a warning and an empty `tests` array.

`data.log_iqtree` points to the IQ-TREE `.iqtree` report file and follows the naming style used by `tree ml iqtree` per-file results. `data.tool_log` points to IQ-TREE's native `.log` file. `data.optimized_trees` points to the IQ-TREE-written candidate trees with optimized branch lengths when detected. `data.merged_candidate_trees` is `candidate.trees` only when multiple candidate tree files were merged by PhyloAI; otherwise it is `null`.

### 9.1 Key Results

| Field | Description |
|-------|-------------|
| `n_candidate_trees` | Number of parsed candidate rows, or `null` if parsing failed. |
| `best_tree_id` | Tree with minimum `delta_likelihood`, or `null`. |
| `n_rejected_au_0_05` | Count of trees with `p_au < 0.05` or `p_au_sign == "-"`, or `null` if AU absent. |
| `replicates` | Resolved RELL replicate count. |
| `model_source` | `model-expr`, `partitions`, or `tool-args` depending on where the final `-m` or `-p` source came from. |

---

## 10. Resume Behavior

`posttree topology` is one IQ-TREE invocation. It does not write a PhyloAI `checkpoint.json`; resume relies on IQ-TREE's native output/checkpoint behavior in the output directory, matching the `tree ml iqtree --matrix` approach.

Rules:

- `--resume` is allowed with an existing non-empty output directory.
- `--overwrite` and `--resume` are mutually exclusive.
- Without `--resume` or `--overwrite`, an existing non-empty output directory exits with code 1.
- On `--resume`, rerun the same IQ-TREE command in the same output directory and let IQ-TREE decide whether to reuse or continue native files.
- After IQ-TREE returns, rebuild `result.json` from current output files.
- No separate `checkpoint.json`, `topology_checkpoint.json`, or PhyloAI parameter-hash file is written in the first implementation.

---

## 11. Exit Codes

| Code | Meaning |
|------|---------|
| `0` | Success |
| `1` | User input error: missing files, invalid parameter combination, or output conflict |
| `2` | IQ-TREE execution failed |
| `3` | IQ-TREE executable not found |

---

## 12. CLI Help Requirements

The help text must not be terse. It should use grouped sections similar to `tree ml iqtree`:

1. **Purpose**: explain that this command tests existing candidate topologies, not tree inference.
2. **Input**: `--matrix`, `--candidate-trees` (comma-separated individual tree files, or a single tree-list file), and the two candidate-tree input forms.
3. **Model Source**: explain `--model-expr` vs `--partitions`, and why PhyloAI does not expose ModelFinder here.
4. **Default Tests**: explicitly state `-n 0 -zb <replicates> -zw -au`, list resulting tests, and explain `--replicates`.
5. **Advanced IQ-TREE Args**: document minimal blocked flags and override behavior.
6. **Examples**: include all six examples from Section 5.
7. **Input Format and Sequence Type**: explain that `--input-format` only affects PhyloAI preflight validation, while explicit IQ-TREE `--seqtype` belongs in `--tool-args` when needed.
8. **Interpretation**: explain that KH/SH/AU are p-values and trees are rejected at p < 0.05; bp-RELL and c-ELW are weights, not p-values.

---

## 13. Implementation Notes

### 13.1 Files to Create

| File | Purpose |
|------|---------|
| `phyloai/posttree/__init__.py` | Posttree package marker. |
| `phyloai/posttree/topology.py` | Library implementation: command builder, runner, parser, result assembly. |
| `phyloai/cli/commands/posttree.py` | Click group and topology command. |
| `docs/commands/posttree-topology.md` | User-facing command documentation. |
| `tests/posttree/test_topology.py` | Library-level tests. |
| `tests/cli/test_posttree_topology.py` | CLI validation/help/dry-run tests. |

### 13.2 Files to Modify

| File | Change |
|------|--------|
| `phyloai/cli/main.py` | Register `posttree` group. |
| `docs/superpowers/specs/2026-06-07-phyloai-design.md` | Already updated with `--candidate-trees`, FcLM under signal, and relaxed global `--tool-args` language. |

### 13.3 Reuse and Extraction

Prefer extracting shared IQ-TREE helpers instead of copying code from `phyloai.tree.ml_iqtree`:

| Existing helper | Use |
|-----------------|-----|
| `_resolve_iqtree_path()` | Resolve custom/PATH/bundled IQ-TREE executable. |
| `_detect_iqtree_version()` | Populate `tool_versions`. |
| `IQTREE_COMPATIBLE_EXTENSIONS` | Validate `--matrix` extension. |
| `_detect_file_format()` | Parse alignment for minimal validation if needed. |
| `_is_flag_overridden()` pattern | Implement suppress-if-present for generated IQ-TREE flags. |

If extraction is small, create `phyloai/tree/iqtree_common.py` or `phyloai/core/iqtree.py`. Keep `posttree/topology.py` focused on topology-test-specific behavior.

### 13.4 Parser Strategy

Parse `.iqtree` `USER TREES` table using a tolerant table parser:

1. Locate line matching `USER TREES`.
2. Locate header row containing `Tree`, `logL`, and `deltaL`.
3. Split rows by whitespace.
4. Recognize value/sign pairs such as `0.7110 +`.
5. Map known columns: `bp-RELL`, `p-KH`, `p-SH`, `c-ELW`, `p-AU`; include weighted columns if IQ-TREE labels them separately.
6. Stop at blank line or explanatory footer.

The parser must preserve raw row text in each parsed item as `raw_line` if practical. This helps users debug IQ-TREE version differences.

---

## 14. Acceptance Criteria

### 14.1 CLI Validation

- [ ] `--matrix` missing exits 1.
- [ ] `--candidate-trees` missing exits 1.
- [ ] one `--candidate-trees candidate.trees` (no comma) passes that file directly to IQ-TREE `-z`.
- [ ] `--candidate-trees h1.nwk,h2.nwk` (comma-separated) writes merged `candidate.trees` in that order and passes it to `-z`.
- [ ] `--input-format` accepts only `auto`, `fasta`, `phylip-relaxed`, `nexus`, or `clustal`.
- [ ] Neither `--model-expr`, `--partitions`, nor model source in `--tool-args` exits 1.
- [ ] Both high-level `--model-expr` and `--partitions` exits 1.
- [ ] `--replicates 999` exits 1.
- [ ] `--overwrite` with `--resume` exits 1.
- [ ] `--tool-args "-s other.fa"` exits 1.
- [ ] `--tool-args "-z other.trees"` exits 1.
- [ ] `--tool-args "--prefix custom -T 30 -fs custom.sitefreq"` is accepted.

### 14.2 Command Building

- [ ] `--model-expr LG+F+R4` builds `-m LG+F+R4`.
- [ ] `--partitions raw.best_model.nex` builds `-p raw.best_model.nex`.
- [ ] `--guide-tree guide.tree` builds `-ft guide.tree`.
- [ ] Default topology-test flags include `-n 0 -zb 10000 -zw -au`.
- [ ] `--replicates 2000` builds `-zb 2000`.
- [ ] `--tool-args "-T 30"` suppresses generated `-T`.
- [ ] `--tool-args "-zb 5000"` suppresses generated `-zb`.

### 14.3 Output and JSON

- [ ] Successful non-dry-run writes `result.json`.
- [ ] `params` includes every `run_topology()` argument.
- [ ] `data.cmd` is exact final argv.
- [ ] `data.tool_stderr` captures merged stdout/stderr.
- [ ] `.iqtree` and `.log` paths are referenced when present.
- [ ] Parsed `USER TREES` rows populate `data.tests[]`.
- [ ] Parsing failure leaves status `success` if IQ-TREE succeeded, with warning.

### 14.4 Resume

- [ ] `--resume` is accepted with an existing non-empty output directory.
- [ ] `--resume` does not require or write PhyloAI `checkpoint.json`.
- [ ] `--resume` reruns the same IQ-TREE command in the existing output directory and rebuilds `result.json` from current outputs.

---

## 15. Relationship to Other Modules

- **`tree ml iqtree`**: infers ML trees and selects/uses models. `posttree topology` consumes the alignment, candidate trees, and model decisions from that stage.
- **`tree cf`**: computes support/concordance on a reference tree. `posttree topology` tests whether alternative topologies are significantly worse.
- **`posttree signal`**: future location for FcLM and other phylogenetic signal diagnostics.
- **`report`**: will consume `result.json` fields such as rejected topology count and AU p-values.
