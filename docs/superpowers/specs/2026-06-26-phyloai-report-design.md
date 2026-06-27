# PhyloAI Report Module Design

**Date:** 2026-06-26  
**Last updated:** 2026-06-27 (implementation iteration: filter/tree/posttree template enrichment, module path fixes, key_results merge extensions, HTML polish, CSV embedding)  
**Status:** Implemented (phyloai 0.1.0)  
**Depends on:** All analysis phases (pretree, tree, posttree) finalized

---

## 1. Overview and Purpose

`phyloai report` generates a reproducible, auditable record of a PhyloAI analysis run. It serves three audiences:

1. **Authors** — produces a Methods paragraph draft suitable for journal submission, describing tools, versions, key parameters, and outcomes for each analytical step
2. **Users (self-audit)** — produces a structured record of every step, parameter, and result for reproducibility verification and archival
3. **AI/MCP diagnostics** — produces a machine-readable `report.json` that serves as a structured entry point for AI-assisted troubleshooting, quality assessment, and guided re-analysis

The command is a single invocation; all internal logic (collection, template rendering, HTML generation) is transparent to the user.

---

## 2. User Interface

Single command, two usage patterns (v1):

```bash
# Single pipeline run (phyloai run — two-layer structure)
phyloai report --run-dir ./runs/run/faa

# Single module run (pretree/tree/posttree — one or more steps)
phyloai report --run-dir ./runs/pretree
```



Output is always written to `<run-dir>/report/`:

```
<run-dir>/report/
├── report.json    # Machine-readable source of truth; AI/MCP diagnostic entry point
└── report.html    # Human-readable report; embeds PDF figures, sortable/collapsible tables
```

`report.json` is the source of truth. `report.html` is fully derived from `report.json` and can be re-rendered at any time without re-scanning the run directory.

---

## 3. Module Structure

```
phyloai/
├── report/
│   ├── __init__.py
│   ├── collector.py       # Directory scanning, result.json discovery, step ordering
│   ├── templates.py       # Per-command methods text generation (Python functions)
│   ├── schema.py          # ReportRecord dataclass and related types
│   ├── renderer.py        # report.json → report.html via Jinja2
│   └── html/
│       └── report.html.j2 # Jinja2 HTML template
cli/commands/report.py     # CLI entry point
```

---

## 4. Data Flow

```
--run-dir
    │
    ▼
collector.py
  Detect run_mode (pipeline / module)
  Discover all result.json files at correct depth
  Parse and order steps by STEP_ORDER
    │
    ▼
templates.py
  For each step: generate methods_text from params + key_results + tool_versions
  Concatenate ordered methods_text → methods_paragraph
    │
    ▼
schema.py
  Assemble ReportRecord
  Write report.json
    │
    ▼
renderer.py
  Read report.json
  Render via Jinja2 → report.html
  Write report.html
```

---

## 5. Directory Detection and run_mode

`collector.py` determines `run_mode` from the structure of `--run-dir` using the following **ordered** checks (earlier checks take priority):

| Priority | Condition | run_mode |
|----------|-----------|----------|
| 1 | `run-dir/result.json` exists **AND** at least one subdir also contains `result.json` | `pipeline` (`phyloai run` two-layer) |
| 2 | `run-dir/result.json` exists, no subdirs with `result.json` | `module` (single command) |
| 3 | No `run-dir/result.json`, but subdirs contain `result.json` | `module` (multi-step module, e.g. pretree) |
| 4 | No `result.json` found at any expected depth | Error: not a valid run directory |

**Rationale for priority 1 before priority 2:** Both `pipeline` and single-command `module` runs have a top-level `result.json`. The distinguishing feature of `pipeline` is the simultaneous presence of per-step subdirectory `result.json` files. Checking for subdirectory files first prevents a pipeline run directory from being misclassified as a single-command run.

**Pipeline step detection is purely filesystem-based.** The top-level `result.json:data.steps[]` is read only for optional metadata (mode, speed) enrichment; step `result.json` paths come from directory scanning, not from `data.steps[]`. This avoids coupling report discovery to the internal `phyloai run` data format.

Scan depth (implemented as a BFS walk): default `max_depth=2`, excluding `report/`, `logs/`, and dot-prefixed directories.

**`step_id` parsing** from the `command` field uses a known-root lookup table. Flag tokens (starting with `-`) are dropped; the first token matching a known root (`pretree`, `tree`, `posttree`, `run`, `doctor`) determines the root command, and subsequent known subcommand tokens build the full `step_id`. Flag values (e.g. `./runs` from `--run-dir ./runs`) are harmless noise that don't match any root. Boolean flags before the root (e.g. `--quiet pretree align`) are handled correctly because only the flag token is dropped, not the subsequent root token.

---

## 6. Step Ordering

Steps are sorted by `STEP_ORDER` to ensure `methods_paragraph` reads in logical analytical sequence. Steps not in the list (future commands) are appended at the end.

```python
STEP_ORDER = [
    "pretree.convert",
    "pretree.stats",
    "pretree.align",
    "pretree.trim",
    "pretree.metrics",
    "pretree.filter.taper",
    "pretree.filter.treeshrink",
    "pretree.filter.symtest",
    "pretree.filter.metrics",
    "pretree.filter.cluster",
    "pretree.concat",
    "tree.ml.fasttree",
    "tree.ml.iqtree",
    "tree.msc",
    "tree.bi",
    "tree.cf",
    "posttree.topology",
    "posttree.dating.hessian",
    "posttree.dating.mcmc",
    # Phase 4 commands not yet implemented — listed here so they are ordered
    # correctly when implemented rather than appended as unknowns
    "posttree.signal",
    "posttree.syserror.brlen",
    "posttree.syserror.cca",
    "posttree.syserror.sites",
    "posttree.simulate",
]
```

`step_id` is parsed from the `command` field of `result.json` by taking the first 2–3 positional tokens after `phyloai`, stripping flags and path arguments:

```
"phyloai pretree align ..."          → "pretree.align"
"phyloai pretree filter taper ..."   → "pretree.filter.taper"
"phyloai tree ml iqtree ..."         → "tree.ml.iqtree"
"phyloai posttree dating mcmc ..."   → "posttree.dating.mcmc"
"phyloai run ..."                    → "run"
```

---

## 7. Incomplete Run Handling

`report` operates in permissive mode by default:

- Steps with `status: "error"` are included in `report.json` with full error details
- Failed steps have `methods_text: ""` and are excluded from `methods_paragraph`
- `report.html` marks failed steps with a visible `[FAILED]` indicator and expands their cards by default
- `pipeline_summary.status` is set to `"partial"` if any step failed, `"complete"` if all succeeded
- Report generation never aborts due to step failures

---

## 8. report.json Schema

`report.json` is the authoritative record of the run. All fields are described below.

```json
{
  "phyloai_version": "0.3.1",
  "generated_at": "2026-06-26T14:23:00Z",
  "run_dir": "/abs/path/to/runs/run/faa",
  "run_mode": "pipeline",
  "status": "complete",

  "pipeline_summary": {
    "n_steps_total": 12,
    "n_steps_success": 11,
    "n_steps_failed": 1,
    "n_steps_skipped": 0,
    "total_wall_time": 3842.5,
    "input_genes": 312,
    "genes_after_filter": 187,
    "final_tree": "runs/run/faa/tree/ml/iqtree/ml/best.treefile"
  },

  "steps": [
    {
      "step_id": "pretree.align",
      "command": "phyloai pretree align --seq-dir ... --method linsi ...",
      "status": "success",
      "wall_time": 31.4,
      "tool_versions": {"mafft": "7.526", "trimal": "1.5.rev1"},
      "params": {"seq_dir": "runs/pretree/1-convert/faa/seqs", "method": "linsi", "seq_type": "AA", "backtrans": true, "nt_dir": "runs/pretree/1-convert/fna/seqs", "threads": 8, "tool_args": null, "mafft_path": null, "magus_path": null, "trimal_path": null, "resume": false, "overwrite": true, "dry_run": false, "quiet": false},
      "key_results": {"n_aligned": 1066, "n_skipped": 0, "mean_alignment_length": 591.5},
      "methods_text": "Multiple sequence alignments were performed using MAFFT v7.526...",
      "output_files": {},
      "warnings": []
    },
    {
      "step_id": "pretree.metrics",
      "command": "phyloai pretree metrics ...",
      "status": "success",
      "wall_time": 18.2,
      "tool_versions": {},
      "params": {"msa_dir": "runs/pretree/4-trim/seqs/faa", "tree_dir": null, "seq_type": "AA", "threads": 4, "output_dir": "runs/pretree/5-metrics/faa", "overwrite": false, "dry_run": false, "quiet": false},
      "key_results": {"errors": 0},
      "methods_text": "Phylogenetic informativeness metrics were computed...",
      "output_files": {
        "metrics_table": {
          "path": "/abs/path/runs/pretree/5-metrics/faa/metrics.csv",
          "description": "Phylogenetic informativeness metrics per locus"
        },
        "correlation_heatmap": {
          "path": "/abs/path/runs/pretree/5-metrics/faa/correlation_heatmap.pdf",
          "description": "Spearman correlation heatmap of all metrics"
        }
      },
      "warnings": []
    }
  ],

  "methods_paragraph": "Raw sequence files were converted... [full concatenated paragraph]",
  "methods_blocks": [
    {"step_id": "pretree.convert", "text": "Raw sequence files were converted...", "step_index": 0},
    {"step_id": "pretree.stats", "text": "Sequence statistics were computed...", "step_index": 1}
  ],

  "figures_index": [
    {
      "figure_id": "Fig-3.1",
      "step_id": "pretree.metrics",
      "label": "correlation_heatmap",
      "caption": "Correlation heatmap of phylogenetic informativeness metrics",
      "description": "Spearman correlation heatmap of all metrics",
      "path": "/abs/path/runs/pretree/5-metrics/faa/correlation_heatmap.pdf",
      "type": "pdf"
    },
    {
      "figure_id": "Fig-3.2",
      "step_id": "pretree.filter.cluster",
      "label": "umap_scatter",
      "caption": "UMAP clustering of loci by phylogenetic informativeness",
      "description": "UMAP projection of loci colored by cluster assignment",
      "path": "/abs/path/runs/pretree/6-filter/cluster/umap_scatter.pdf",
      "type": "pdf"
    }
  ],

  "tables_index": [
    {
      "table_id": "Table-3.1",
      "step_id": "pretree.metrics",
      "label": "metrics_table",
      "caption": "Phylogenetic informativeness metrics per locus",
      "description": "Phylogenetic informativeness metrics per locus",
      "path": "/abs/path/runs/pretree/5-metrics/faa/metrics.csv",
      "type": "csv"
    },
    {
      "table_id": "Table-3.2",
      "step_id": "pretree.filter.taper",
      "label": "filter_decisions",
      "caption": "Per-locus TAPER masking decisions",
      "description": "Per-locus TAPER masking decisions",
      "path": "/abs/path/runs/pretree/6-filter/taper/filter_decisions.csv",
      "type": "csv"
    }
  ]
}
```

**Field notes:**

- `run_mode`: `"pipeline"` | `"module"`
- `status`: `"complete"` | `"partial"` | `"failed"`
- `steps[].params`: the **complete** `params` dict copied verbatim from `result.json` (all parameters including `threads`, paths, flags); this preserves full reproducibility. Methods templates read from this complete dict but only describe scientifically meaningful parameters in the generated text — technical parameters are ignored inside the template function.
- `steps[].methods_text`: empty string `""` for failed steps; excluded from `methods_paragraph`
- `steps[].output_files`: copied from `result.json:data.output_files`. Non-dict entries (legacy bare ints, redundant keys like `n_plots`) are purged at assembly time. Each retained value is a `{path, description}` object.
- `steps[].key_results`: may include values merged from `data.summary` (e.g. convert module) or flat `data.*` scalars (e.g. stats single-file). Nested numeric dicts (e.g. `length_before: {mean, min, max}`) are flattened into `length_before_mean`, `length_before_min`, etc. This enrichment runs in both the CLI handler (before template generation) and `assemble_report` (for report.json fidelity).
- `methods_paragraph`: plain text concatenation of all successful step methods texts; intended for copy-to-clipboard.
- `methods_blocks`: annotated per-step methods list. Each entry has `step_id`, `text`, and `step_index` (index into `steps` for anchor linking). HTML Panel B renders each block as a separate paragraph with a clickable `[step_id]` tag linking to the Step Detail card. No text deduplication — each result.json is an independent analysis.
- `figure_id` and `table_id`: sequential numbering by section group (see Section 11), e.g. `Fig-3.1`, `Table-3.1`
- `figures_index`: built by filtering `output_files` for `.pdf`/`.png` extensions.  
- `tables_index`: built by filtering `output_files` for `.csv`/`.tsv` extensions. Small CSV files (≤200 rows, ≤500KB) are additionally embedded inline in the HTML Step Detail cards as sortable tables.
- `label`: the snake_case key from `data.output_files` in the source `result.json`

---

## 9. Methods Text Templates

Each `step_id` maps to a dedicated Python function in `templates.py`. Templates read `params`, `key_results`, and `tool_versions` from the step's `result.json` and produce 2–5 sentences of academic English.

**Design principles:**

- Each template function receives the **complete** `params` dict (as stored in `steps[].params`) but only reads the scientifically meaningful keys; technical parameters (`threads`, paths, `--quiet`, `--overwrite`, `--resume`, `--dry-run`) are ignored inside the function
- The scientific parameters and key results each template reads are defined in that command's own spec; templates follow the principle, not a central registry
- All scientifically meaningful parameters are described, whether or not they differ from defaults
- Parameter descriptions include the parameter's scientific meaning, not just its value
- Conditional branches handle parameter combinations (e.g. `backtrans=True`, `partitioned=True`, `modelfinder=MFP`)
- All placeholder values have fallbacks: `tool_versions.get("mafft", "unknown version")`
- New commands require adding one function to `templates.py`; no other files change

### 9.1 Template Examples

`pretree.align`:
> Multiple sequence alignments were performed using MAFFT v{mafft} with the {method_description} algorithm, which {method_rationale}. A total of {n_aligned} {seq_type} loci were aligned{skipped_clause}. Mean alignment length was {mean_alignment_length:.1f} bp across a mean of {mean_n_taxa:.1f} taxa per locus.{backtrans_clause}

Where `method_description` and `method_rationale` are mapped per strategy:
- `linsi` → "L-INS-i" / "applies iterative local pairwise alignment refinement and is suited for sequences with conserved domains and insertions"
- `einsi` → "E-INS-i" / "uses multiple local alignments and is suited for sequences with multiple conserved regions separated by unalignable regions"
- `ginsi` → "G-INS-i" / "applies global pairwise alignment and is suited for sequences of similar length without large insertions"
- `fftns1` → "FFT-NS-1" / "uses progressive alignment with single FFT iteration and is suited for large datasets where speed is prioritized"
- `fftns2` → "FFT-NS-2" / "uses progressive alignment with two FFT iterations"
- `auto` → "auto-selected" / "strategy selected automatically by MAFFT based on sequence length and count"
- `magus` → "MAGUS" / "uses graph-based divide-and-conquer alignment and is suited for very large or highly divergent datasets"

And `backtrans_clause`:
- `True` → " Codon-aware nucleotide alignments were produced via back-translation using trimAl v{trimal}, preserving reading frame in the nucleotide alignments."
- `False` → ""

`pretree.filter.taper`:
> Aligned sequences were screened for compositional bias and systematic sequencing errors using TAPER v{taper} (correction_multi.jl, executed via Julia v{julia}). TAPER applies a moving-window approach to identify and mask amino acid sites within individual sequences that deviate from expected substitution patterns, without discarding entire loci. The masking stringency cutoff was set to {cutoff} (`-c {cutoff}`){cutoff_note}. Of {n_input} input loci, {n_retained} were retained{dropped_clause}{masked_clause}.

`pretree.concat`:
> Trimmed alignments were concatenated into a supermatrix using phyloai concat. Loci were included only if they met a minimum taxon occupancy threshold of {taxa_occupancy:.0%} (`--taxa-occupancy {taxa_occupancy}`), requiring sequence data for at least {min_taxa} of {total_taxa} taxa; {n_msa_dropped} loci were excluded for failing this criterion. The final supermatrix comprised {n_msa_used} loci across {n_taxa} taxa with a total alignment length of {total_length:,} {seq_type} positions (gap ratio: {gap_ratio:.1%}; parsimony-informative sites: {pi_ratio:.1%}).{recoding_clause}{outgroup_clause}

Where `recoding_clause`:
- not None → " To reduce the influence of substitution saturation and compositional heterogeneity, sequences were recoded into {recoding} categories (`--recoding {recoding}`), collapsing the 20 standard amino acids into {n_groups} biochemically similar groups; both the original and recoded matrices were retained."
- None → ""

`tree.ml.iqtree` (supermatrix, partitioned, MFP):
> Maximum likelihood phylogenetic inference was performed using IQ-TREE v{iqtree}. {partition_clause} Substitution models were selected {modelfinder_clause} from a candidate set comprising {mset} matrix models (`--mset {mset}`). {model_result_clause} Branch support was assessed using {boot:,} ultrafast bootstrap replicates (`-B {boot}`){alrt_clause}. The final log-likelihood of the best tree was {log_likelihood:.2f}.

Where `partition_clause`:
- `partitioned=True, merged_partitions=True` → "A partitioned analysis was conducted with partition merging enabled (`--merge`), using the rclusterf algorithm (`--rclusterf {rclusterf}`) to identify the optimal merging scheme by evaluating {rclusterf}% of candidate partition pairs."
- `partitioned=True, merged_partitions=False` → "A partitioned analysis was conducted using the provided partition scheme."
- `partitioned=False` → ""

`pretree.metrics` (pure-Python, no external tool, representative of non-tool steps):
> Phylogenetic informativeness metrics were computed for {n_markers} loci across {n_metrics} dimensions. Evaluated metrics included locus length, number of informative sites, gap percentage, GC content, and RCFV (relative composition frequency variability){tree_metrics_clause}. Pairwise Spearman correlations were computed across all metrics and visualized as a heatmap for diagnostic evaluation of metric redundancy and complementarity.

Where `tree_metrics_clause`:
- `tree_dir is not None` → " , as well as tree-based metrics including treeness and Robinson-Foulds distance between gene trees and the reference species tree"
- `tree_dir is None` → ""

`posttree.dating.mcmc` (multi-run, convergence, representative of complex multi-output commands):
> Divergence time estimation was performed using MCMCTree (PAML v{paml}) under a {model_descr} substitution model with a {clock_descr} molecular clock. {n_runs} independent MCMC chains were run for {burnin:,} burn-in generations followed by {sample:,} sampling generations, with parameters sampled every {sample_freq} generations. Convergence was assessed using {diag_descr}: effective sample sizes (ESS) for all parameters exceeded 200, and the potential scale reduction factor (PSRF) approached 1.0. Posterior node age estimates and {n_cred_intervals}% highest posterior density (HPD) intervals are summarised in the node ages table; MCMC trace plots for posterior and prior are provided in the supplementary figures.

Where `model_descr` and `clock_descr` are mapped from params:
- model: `"JC69"`, `"HKY85"`, etc. → verbatim model name
- clock: `"strict"` → "strict" / `"independent"` → "independent-rates" / `"correlated"` → "autocorrelated-rates"
- diag_descr: combined from `params.diag` or key_results; e.g. "Gelman-Rubin diagnostic" or "trace inspection and ESS calculation"

---

## 10. HTML Report Structure

The HTML report is a self-contained single file (no external CDN dependencies) renderable offline in any modern browser. PDF figures are embedded natively; only minimal JavaScript is used (table sorting, clipboard copy).

### Page layout

The HTML report has five named **panels**. These panel names are distinct from the **figure/table numbering groups** defined in Section 11 (which use `Fig-3.x`, `Fig-4.x` etc. based on analytical phase, not HTML panel position).

```
┌──────────────────────────────────────────────┐
│ Header: PhyloAI Report                        │
│ run_dir · generated_at · phyloai_version      │
├──────────────────────────────────────────────┤
│ Panel A. Run Summary                          │
├──────────────────────────────────────────────┤
│ Panel B. Methods                              │
├──────────────────────────────────────────────┤
│ Panel C. Steps Detail                         │
├──────────────────────────────────────────────┤
│ Panel D. Figures                              │
├──────────────────────────────────────────────┤
│ Panel E. Output Files Index                   │
└──────────────────────────────────────────────┘
```

### Panel A — Run Summary

Summary cards:

```
Steps: 11/12 succeeded · 1 failed     Total wall time: 1h 03m 42s
Input genes: 1,066 → After filters: 1,039
Final tree: runs/tree/ml/iqtree/ml/matrix.aa.treefile  [link]
```

Failed steps listed by name in red below the cards. For `run_mode: pipeline`, a linear progress indicator shows which steps succeeded/failed.

### Panel B — Methods

Uses **`methods_blocks`** from `report.json`. Each block is rendered as its own paragraph with a clickable `[step_id]` badge linking to the corresponding Step Detail card via anchor (`#step-{index}`). A "Copy to clipboard" button copies the plain text `methods_paragraph` (without badges). Intended for direct use in manuscript Methods sections.

### Panel C — Steps Detail

One collapsible card per step, ordered by `STEP_ORDER`:

```html
<details id="step-{index}">  <!-- success: collapsed by default -->
<details id="step-{index}" open>  <!-- failed: expanded by default -->
  <summary>
    [✓|✗] {step_id} · {primary_tool} v{version} · {wall_time}s
  </summary>

  <!-- Link back to the corresponding Methods paragraph (Panel B) -->

  <!-- Parameters table (scientific params only; >10 rows wrapped) -->

  <!-- Key Results table -->

  <!-- Embedded CSV tables (≤200 rows, ≤500KB) -->

  <!-- Warnings (if any), styled amber -->

  <!-- Full CLI command, monospace, copyable -->
</details>
```

Steps do **not** repeat the methods text from Panel B; instead they provide a `↑ Methods` back-link. This keeps Panel B (the manuscript-readable Methods draft) and Panel C (the per-step audit trail) cleanly separated.

### Panel D — Figures

Each figure from `figures_index` is rendered as:

```html
<figure id="Fig-3.1">
  <object data="../../pretree/6-metrics/faa/correlation.pdf"
          type="application/pdf"
          width="100%" height="600px">
  </object>
  <figcaption>
    <strong>Figure 3.1</strong> Correlation heatmap of phylogenetic
    informativeness metrics across 1,066 loci. (pretree.metrics)<br>
    <a href="../../pretree/6-metrics/faa/correlation.pdf">
      runs/pretree/6-metrics/faa/correlation.pdf
    </a>
  </figcaption>
</figure>
```

Figure numbering uses the analytical-phase groups from Section 11 (e.g. all pretree figures are `Fig-3.x`, all tree figures `Fig-4.x`). These phase numbers are independent of the HTML panel labels (A–E).

### Panel E — Output Files Index

Sortable HTML table of all output files across all steps. Rows beyond 20 are collapsed under a `<details>` wrapper.

| # | Step | Label | Description | Path | Type |
|---|------|-------|-------------|------|------|
| 1 | pretree.metrics | metrics_table | Phylogenetic informativeness metrics per locus | runs/pretree/5-metrics/faa/metrics.csv | CSV |
| 2 | pretree.metrics | correlation_heatmap | Spearman correlation heatmap of all metrics | runs/pretree/5-metrics/faa/correlation_heatmap.pdf | PDF |
| 3 | pretree.concat | matrix_faa | Concatenated supermatrix of 52 taxa across 1,039 loci | runs/pretree/7-concat/matrix.fa | FASTA |
| 4 | tree.ml.iqtree | ml_tree | Maximum likelihood species tree in Newick format | runs/tree/ml/iqtree/ml/matrix.aa.treefile | treefile |

The `Description` column is populated from the `description` field in each `output_files` entry. Omitted descriptions render as `—`.

**Table caption:**

> **Table E.1** Output files produced across all steps.  
> Source: aggregated from `result.json:data.output_files` across all steps.

Table headers are clickable for client-side sorting (ascending/descending).

Tables within Step Detail cards that exceed 20 rows are similarly wrapped:

```html
<details>
  <summary>Table 3.2 — Show all 1,066 gene alignment records</summary>
  <table>...</table>
  <caption>
    <strong>Table 3.2</strong> Per-locus alignment results from pretree.align.<br>
    Source: <a href="../../pretree/2-align/result.json">
      runs/pretree/2-align/result.json
    </a> → data.files
  </caption>
</details>
```

---

## 11. Figure and Table Numbering Convention

Figure and table IDs use **analytical phase numbers**, which are independent of the HTML panel labels (A–E). Phase numbers reflect the analytical stage of the run, not the visual layout of the report.

| Phase | Analytical group | Figure/Table prefix |
|-------|-----------------|---------------------|
| 3 | pretree.* | Fig-3.x / Table-3.x |
| 4 | tree.* | Fig-4.x / Table-4.x |
| 5 | posttree.* | Fig-5.x / Table-5.x |

Within each phase, numbering is sequential in `STEP_ORDER` order. Phase numbers 1 and 2 are intentionally unused (reserved for run summary and methods, which produce no figures or tables).

---

## 12. Output Files by Module

Every command records all persistent output files it produces under `data.output_files` in its `result.json` (JSON Output Standard Section 5.4). Each entry has a required `"path"` and an optional `"description"` (see Section 8 field notes). `collector.py` reads `data.output_files` from each step to populate `figures_index` (filtering for `.pdf`/`.png` entries) and `tables_index` (filtering for `.csv`/`.tsv` entries) in `report.json`, preserving `description` verbatim. Files of other types (FASTA, Newick, TXT, etc.) appear in the Output Files Index panel but are not embedded as figures. No paths are hardcoded in the report module.

Known output files per module (for reference; authoritative source is always `result.json:data.output_files`):

| Module | `data.output_files` labels | Type | Description |
|--------|---------------------------|------|-------------|
| `pretree.stats` | `per_gene_table` | CSV | Per-locus length, taxon count, gap ratio statistics |
| `pretree.metrics` | `metrics_table` | CSV | All phylogenetic informativeness metrics per locus |
| `pretree.metrics` | `correlation_heatmap` | PDF | Spearman correlation heatmap of all metrics (produced by `metrics analyze`) |
| `pretree.metrics` | `correlation_matrix` | CSV | Pairwise Spearman correlation matrix in tabular format |
| `pretree.metrics` | `basic_statistics` | CSV | Per-metric summary statistics |
| `pretree.metrics` | `plots_dir` | directory | Distribution plots for each computed metric |
| `pretree.filter.taper` | `retained_loci`, `dropped_loci`, `filter_decisions` | CSV | Per-locus TAPER masking decisions |
| `pretree.filter.treeshrink` | `retained_loci`, `dropped_loci`, `modified_loci`, `removed_taxa`, `filter_decisions` | CSV | Per-locus TreeShrink pruning decisions |
| `pretree.filter.symtest` | `retained_loci`, `dropped_loci`, `filter_decisions` | CSV | Per-locus symmetry test results |
| `pretree.filter.metrics` | `retained_loci`, `dropped_loci`, `filter_decisions` | CSV | Per-locus metric-rule filtering decisions |
| `pretree.filter.cluster` | `features_used`, `reduction`, `clusters`, `cluster_summary`, `cluster_metric_means`, `cluster_metric_heatmap`, `cluster_selection`*, `umap_replicates`*, `cluster_scatter_*`, `cluster_metric_boxplots_*`, `cluster_*` | CSV/PDF | UMAP cluster assignments, diagnostics, and plots |
| `pretree.filter.cluster` | `outlier_retained_loci`, `outlier_dropped_loci`, `outlier_filter_decisions`, `outlier_comparison`, `outlier_wilcoxon`, `outlier_boxplots_*` | CSV/PDF | Outlier cluster diagnostic files (auto-drop mode only) |
| `pretree.concat` | `matrix_original`, `partitions_original`; `matrix_recoded`*, `partitions_recoded`*; `matrix_translated`*, `partitions_translated`*; `matrix_cds12`*, `partitions_cds12`* | FASTA/TXT | Supermatrix and partition files per variant; * = conditional on params |
| `posttree.topology` | `topology_test_results`, `iqtree_report`, `iqtree_log`, `optimized_trees` | CSV/TXT | AU/WKH/WSH test results and IQ-TREE output files |
| `posttree.dating.hessian` | `iqtree.dummy.phy`, `iqtree.rooted.nwk`, `iqtree.mcmctree.hessian` | PHY/NWK/TXT | Hessian computation outputs for MCMCTree |
| `posttree.dating.mcmc` | `trace_*`, `convergence_*`, `infinite_sites_*`, `posterior_vs_prior_*`, `spearman_correlations` | PDF/CSV | MCMC diagnostic plots and statistics |
| `tree.bi` | `trace_plots`, `bpcomp_*` (e.g. `bpcomp_all_bpdiff`, `bpcomp_all_con_tre`), `tracecomp_*` (e.g. `tracecomp_all_contdiff`) | PDF/TXT | PhyloBayes MCMC trace plots and convergence diagnostics |

---

## 13. Adding New Commands

When a new phyloai command is added:

1. Add its `step_id` to `STEP_ORDER` in `collector.py` at the correct position
2. Add a `generate_methods_<step_id>(params, key_results, tool_versions) -> str` function to `templates.py`
3. Ensure the command records all persistent user-facing output files under `data.output_files` as objects with `"path"` (required) and `"description"` (recommended) in its `result.json`; `collector.py` auto-populates `figures_index` and `tables_index` from these entries by file extension — no changes to the report module needed
4. No changes to `schema.py`, `renderer.py`, or `report.html.j2` are needed for standard cases

---

## 14. CLI Specification

```
Usage: phyloai report [OPTIONS]

  Generate a reproducible analysis report from a PhyloAI run directory.

  Produces report.json (machine-readable, AI/MCP diagnostic entry point)
  and report.html (human-readable, with embedded figures and methods draft).

Options:
  --run-dir PATH   Run directory to report on. Auto-detects pipeline run
                   (phyloai run output) or module run (pretree/tree/posttree).
                   [required]
  -o, --output-dir PATH
                   Output directory for report files.
                   [default: <run-dir>/report]
  --overwrite      Overwrite existing report files.
  -q, --quiet      Suppress terminal output except errors.
  --help           Show this message and exit.
```

---

## 15. Out of Scope

The following are explicitly not part of this design:

- **Tree visualization** — phylogenetic trees are listed as file links only; rendering is delegated to downstream tools (FigTree, iTOL, etc.)
- **Cross-run comparison** — comparing two runs with different parameters is not supported in this version
- **LLM-assisted methods polishing** — methods text is deterministic template output; no LLM involvement
- **Figure generation** — report does not generate new figures; only embeds figures produced by analysis commands
- **Sub-commands** — `report` is a single command, not a group with sub-commands (`report collect`, `report methods`, etc.)
