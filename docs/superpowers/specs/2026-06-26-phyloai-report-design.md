# PhyloAI Report Module Design

**Date:** 2026-06-26  
**Status:** Approved  
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

Single command, three usage patterns:

```bash
# Single pipeline run (phyloai run — two-layer structure)
phyloai report --run-dir ./runs/run/faa

# Single module run (pretree/tree/posttree — one-layer structure)
phyloai report --run-dir ./runs/pretree

# Top-level directory — auto-discovers all runs
phyloai report --run-dir ./runs
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
  Detect run_mode (pipeline / module / multi)
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

`collector.py` determines `run_mode` from the structure of `--run-dir`:

| Condition | run_mode |
|-----------|----------|
| `run-dir/result.json` exists directly | `module` (single command) |
| `run-dir/` has subdirs with `result.json`, and `run-dir/result.json` also exists | `pipeline` (`phyloai run` two-layer) |
| `run-dir/` has subdirs with `result.json`, no top-level `result.json` | `module` (multi-step module, e.g. pretree) |
| `run-dir/` has multiple subdirs each matching pipeline or module pattern | `multi` |

Scan depth:

| run_mode | Depth | Typical path |
|----------|-------|--------------|
| `module` | 1–2 levels | `runs/pretree/2-align/result.json` |
| `pipeline` | 2–3 levels | `runs/run/faa/pretree/2-align/result.json` |
| `multi` | recurse, each child treated independently | — |

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
      "params": {"method": "linsi", "backtrans": true, "seq_type": "AA"},
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
      "params": {"msa_dir": "...", "tree_dir": null},
      "key_results": {"errors": 0},
      "methods_text": "Phylogenetic informativeness metrics were computed...",
      "output_files": {
        "metrics_table": "/abs/path/runs/pretree/5-metrics/faa/metrics.csv",
        "correlation_heatmap": "/abs/path/runs/pretree/5-metrics/faa/correlation_heatmap.pdf"
      },
      "warnings": []
    }
  ],

  "methods_paragraph": "Raw sequence files were converted... [full concatenated paragraph]",

  "figures_index": [
    {
      "figure_id": "Fig-3.1",
      "step_id": "pretree.metrics",
      "label": "correlation_heatmap",
      "caption": "Correlation heatmap of phylogenetic informativeness metrics",
      "path": "/abs/path/runs/pretree/5-metrics/faa/correlation_heatmap.pdf",
      "type": "pdf"
    },
    {
      "figure_id": "Fig-3.2",
      "step_id": "pretree.filter.cluster",
      "label": "umap_scatter",
      "caption": "UMAP clustering of loci by phylogenetic informativeness",
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
      "path": "/abs/path/runs/pretree/5-metrics/faa/metrics.csv",
      "type": "csv"
    },
    {
      "table_id": "Table-3.2",
      "step_id": "pretree.filter.taper",
      "label": "filter_decisions",
      "caption": "Per-locus TAPER masking decisions",
      "path": "/abs/path/runs/pretree/6-filter/taper/filter_decisions.csv",
      "type": "csv"
    }
  ]
}
```

**Field notes:**

- `run_mode`: `"pipeline"` | `"module"` | `"multi"`
- `status`: `"complete"` | `"partial"` | `"failed"`
- `steps[].params`: only scientifically meaningful parameters; technical params (`threads`, executable paths) are excluded
- `steps[].methods_text`: empty string `""` for failed steps; excluded from `methods_paragraph`
- `steps[].output_files`: copied directly from `result.json:data.output_files` for that step (see JSON Output Standard Section 5.4); `{}` when the step produces no tables or figures
- `figures_index`: global index of all PDF/PNG figures across all steps, built by filtering `output_files` entries whose paths end in `.pdf` or `.png`; enables AI diagnostics and HTML renderer to locate all figures without traversing individual step records
- `tables_index`: global index of all CSV/TSV tables across all steps, built by filtering `output_files` entries whose paths end in `.csv` or `.tsv`
- `figure_id` and `table_id`: sequential numbering by section group (see Section 11), e.g. `Fig-3.1`, `Table-3.1`
- `label`: the snake_case key from `data.output_files` in the source `result.json`

---

## 9. Methods Text Templates

Each `step_id` maps to a dedicated Python function in `templates.py`. Templates read `params`, `key_results`, and `tool_versions` from the step's `result.json` and produce 2–5 sentences of academic English.

**Design principles:**

- All scientifically meaningful parameters are described, whether or not they differ from defaults
- Parameter descriptions include the parameter's scientific meaning, not just its value
- Conditional branches handle parameter combinations (e.g. `backtrans=True`, `partitioned=True`, `modelfinder=MFP`)
- All placeholder values have fallbacks: `tool_versions.get("mafft", "unknown version")`
- Technical parameters (threads, paths, `--quiet`, `--overwrite`) are never included
- New commands require adding one function to `templates.py`; no other files change

**Template examples (reference quality):**

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

---

## 10. HTML Report Structure

The HTML report is a self-contained single file (no external CDN dependencies) renderable offline in any modern browser. PDF figures are embedded natively; only minimal JavaScript is used (table sorting, clipboard copy).

### Page layout

```
┌──────────────────────────────────────────────┐
│ Header: PhyloAI Report                        │
│ run_dir · generated_at · phyloai_version      │
├──────────────────────────────────────────────┤
│ Section 1. Pipeline Summary                   │
├──────────────────────────────────────────────┤
│ Section 2. Methods                            │
├──────────────────────────────────────────────┤
│ Section 3. Steps Detail                       │
├──────────────────────────────────────────────┤
│ Section 4. Figures                            │
├──────────────────────────────────────────────┤
│ Section 5. Output Files Index                 │
└──────────────────────────────────────────────┘
```

### Section 1 — Pipeline Summary

Summary cards:

```
Steps: 11/12 succeeded · 1 failed     Total wall time: 1h 03m 42s
Input genes: 1,066 → After filters: 1,039
Final tree: runs/tree/ml/iqtree/ml/matrix.aa.treefile  [link]
```

Failed steps listed by name in red below the cards. For `run_mode: pipeline`, a linear progress indicator shows which steps succeeded/failed.

### Section 2 — Methods

Full `methods_paragraph` in a styled block with a "Copy to clipboard" button. Each sentence is linked (via anchor) to its corresponding Step Detail card. Intended for direct use in manuscript Methods sections.

### Section 3 — Steps Detail

One collapsible card per step, ordered by `STEP_ORDER`:

```html
<details>  <!-- success: collapsed by default -->
<details open>  <!-- failed: expanded by default -->
  <summary>
    [✓|✗] {step_id} · {primary_tool} v{version} · {wall_time}s
  </summary>

  <!-- methods_text for this step -->

  <!-- Parameters table (scientific params only) -->
  <!-- If > 10 rows, wrapped in inner <details> -->

  <!-- Key Results table -->

  <!-- Warnings (if any), styled amber -->

  <!-- Full CLI command, monospace, copyable -->
</details>
```

### Section 4 — Figures

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

Figure numbering is sequential within section groups derived from `STEP_ORDER` position (e.g. all pretree figures are `Fig-3.x`, all tree figures `Fig-4.x`).

### Section 5 — Output Files Index

Sortable HTML table of all output files across all steps. Rows beyond 20 are collapsed under a `<details>` wrapper.

| # | Step | Label | Path | Type |
|---|------|-------|------|------|
| 1 | pretree.align | Aligned AA MSAs | runs/pretree/2-align/seqs/faa/ | directory |
| 2 | pretree.concat | Supermatrix | runs/pretree/7-concat/matrix.fa | fasta |
| 3 | tree.ml.iqtree | ML tree | runs/tree/ml/iqtree/ml/matrix.aa.treefile | treefile |

**Table caption:**

> **Table 5.1** Output files produced across all steps.  
> Source: aggregated from `result.json:data.output` and `key_results` across all steps.

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

Section numbers follow `STEP_ORDER` groupings:

| Section | Steps | Figure/Table prefix |
|---------|-------|---------------------|
| 1 | pipeline_summary | — |
| 2 | methods | — |
| 3 | pretree.* | Fig-3.x / Table-3.x |
| 4 | tree.* | Fig-4.x / Table-4.x |
| 5 | posttree.* | Fig-5.x / Table-5.x |
| 6 | output files index | Table-6.x |

Within each section, numbering is sequential in `STEP_ORDER` order.

---

## 12. Output Files by Module

Every command records all CSV/TSV tables and PDF/PNG figures it produces under `data.output_files` in its `result.json` (JSON Output Standard Section 5.4). `collector.py` reads `data.output_files` from each step to populate `figures_index` (`.pdf`/`.png` entries) and `tables_index` (`.csv`/`.tsv` entries) in `report.json`. No paths are hardcoded in the report module.

Known output files per module (for reference; authoritative source is always `result.json:data.output_files`):

| Module | `data.output_files` labels | Type | Description |
|--------|---------------------------|------|-------------|
| `pretree.stats` | `per_gene_stats` | CSV | Per-locus length, taxon count, gap ratio statistics |
| `pretree.metrics` | `metrics_table` | CSV | All phylogenetic informativeness metrics per locus |
| `pretree.metrics` | `correlation_heatmap` | PDF | Spearman correlation heatmap of all metrics |
| `pretree.filter.taper` | `retained_loci`, `dropped_loci`, `filter_decisions` | CSV | Per-locus TAPER masking decisions |
| `pretree.filter.treeshrink` | `retained_loci`, `dropped_loci`, `modified_loci`, `removed_taxa`, `filter_decisions` | CSV | Per-locus TreeShrink pruning decisions |
| `pretree.filter.symtest` | `retained_loci`, `dropped_loci`, `filter_decisions` | CSV | Per-locus symmetry test results |
| `pretree.filter.metrics` | `retained_loci`, `dropped_loci`, `filter_decisions` | CSV | Per-locus metric-rule filtering decisions |
| `pretree.filter.cluster` | `cluster_assignments`, `reduction_summary`, `metric_means` | CSV | UMAP cluster assignments and statistics |
| `pretree.filter.cluster` | `umap_scatter`, `boxplots`, `outlier_diagnostics` | PDF | UMAP and cluster diagnostic plots |
| `posttree.topology` | `topology_test_results` | CSV | AU/WKH/WSH test p-values per candidate tree |
| `posttree.dating.mcmc` | `trace_run*_posterior`, `trace_run*_prior` | PDF | MCMCtree MCMC trace plots per run |
| `posttree.dating.mcmc` | `convergence_summary`, `node_ages` | CSV | Convergence statistics and posterior node age estimates |
| `tree.bi` | `trace_plots`, `convergence_render` | PDF/TXT | MrBayes MCMC trace and convergence diagnostics |

---

## 13. Adding New Commands

When a new phyloai command is added:

1. Add its `step_id` to `STEP_ORDER` in `collector.py` at the correct position
2. Add a `generate_methods_<step_id>(params, key_results, tool_versions) -> str` function to `templates.py`
3. Ensure the command records all CSV/TSV tables and PDF/PNG figures under `data.output_files` in its `result.json`; `collector.py` auto-populates `figures_index` and `tables_index` from these entries by file extension — no changes to the report module needed
4. No changes to `schema.py`, `renderer.py`, or `report.html.j2` are needed for standard cases

---

## 14. CLI Specification

```
Usage: phyloai report [OPTIONS]

  Generate a reproducible analysis report from a PhyloAI run directory.

  Produces report.json (machine-readable, AI/MCP diagnostic entry point)
  and report.html (human-readable, with embedded figures and methods draft).

Options:
  --run-dir PATH   Run directory to report on. Auto-detects single run,
                   pipeline run, or multi-run directory structure.
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
