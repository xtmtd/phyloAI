# PhyloAI AI Integration Design Specification

**Date:** 2026-06-27
**Status:** Approved for implementation
**Prerequisite:** All non-`doctor` CLI commands stable and writing `result.json` (Phases 1–6 complete)

---

## 1. Overview

The AI integration layer sits above the CLI and consists of two components with strict separation of concerns:

```
用户 ←→ Skill (对话/决策/解读) ←→ MCP Server (执行桥) ←→ phyloai CLI ←→ 文件系统
```

| Layer          | Responsibility                                                                 | Does NOT do                          |
| -------------- | ------------------------------------------------------------------------------ | ------------------------------------ |
| **Skill**      | Conversation flow, parameter confirmation, result interpretation, next-step suggestion, syserror orchestration | Direct tool invocation, path management |
| **MCP Server** | Execute CLI commands, read result.json, return structured results              | Decision-making, result interpretation |
| **CLI**        | Computation, writing result.json                                               | Conversation                         |

Development order: MCP Server (Phase 7) → Skill (Phase 8), after all CLI modules are stable.

---

## 2. MCP Server Design

### 2.1 Tool Granularity

One MCP tool per CLI subcommand. Example mapping:

| MCP Tool                    | CLI Command                         |
| --------------------------- | ----------------------------------- |
| `pretree_convert`           | `phyloai pretree convert`           |
| `pretree_stats`             | `phyloai pretree stats`             |
| `pretree_align`             | `phyloai pretree align`             |
| `pretree_trim`              | `phyloai pretree trim`              |
| `pretree_metrics`           | `phyloai pretree metrics`           |
| `pretree_filter_taper`      | `phyloai pretree filter taper`      |
| `pretree_filter_treeshrink` | `phyloai pretree filter treeshrink` |
| `pretree_filter_metrics`    | `phyloai pretree filter metrics`    |
| `pretree_filter_symtest`    | `phyloai pretree filter symtest`    |
| `pretree_filter_cluster`    | `phyloai pretree filter cluster`    |
| `pretree_concat`            | `phyloai pretree concat`            |
| `tree_ml_iqtree`            | `phyloai tree ml iqtree`            |
| `tree_ml_fasttree`          | `phyloai tree ml fasttree`          |
| `tree_bi_pb`               | `phyloai tree bi pb`                |
| `tree_bi_bpcomp`            | `phyloai tree bi bpcomp`            |
| `tree_bi_tracecomp`         | `phyloai tree bi tracecomp`         |
| `tree_bi_readpb`            | `phyloai tree bi readpb`            |
| `tree_msc`                  | `phyloai tree msc`                  |
| `tree_cf`                   | `phyloai tree cf`                   |
| `posttree_topology`           | `phyloai posttree topology`         |
| `posttree_dating_hessian`    | `phyloai posttree dating hessian`   |
| `posttree_dating_mcmc`       | `phyloai posttree dating mcmc`      |
| `posttree_signal`           | `phyloai posttree signal`           *(stub until CLI implemented)* |
| `posttree_simulate`         | `phyloai posttree simulate`         *(stub until CLI implemented)* |
| `posttree_syserror_brlen`   | `phyloai posttree syserror brlen`   *(stub until CLI implemented)* |
| `posttree_syserror_cca`     | `phyloai posttree syserror cca`     *(stub until CLI implemented)* |
| `posttree_syserror_sites`   | `phyloai posttree syserror sites`   *(stub until CLI implemented)* |
| `doctor`                    | `phyloai doctor --output-format json` |
| `report`                    | `phyloai report`                    |
| `run`                       | `phyloai run`                       |
| `check_status`              | *(reads checkpoint.json / result.json — no CLI call)* |
| `read_result`               | *(reads result.json at given path — no CLI call)* |
| `read_report`               | *(reads report.json at given path — no CLI call)* |
| `get_command_schema`        | *(introspects Click command tree — no CLI call)* |

**Stub vs real:** Stub tools exist only for commands listed in the public CLI but not yet implemented in the installed package. Once a command is implemented, the dynamic schema generation exposes the real tool automatically — no MCP-side change required. Stub tools return `{"status": "not_implemented", "message": "This command is not yet available in the installed version."}` immediately without invoking the CLI.

### 2.2 Schema Source

MCP tool parameter schemas are generated dynamically at server startup by directly importing the phyloai Click application and introspecting the command tree. No hand-written schemas, no separate export scripts, no static JSON files.

```python
# Conceptual implementation
import click
from phyloai.cli.main import cli

def build_mcp_tool(command: click.Command) -> dict:
    return {
        "name": ...,
        "description": command.help,
        "inputSchema": {
            "type": "object",
            "properties": {
                param.name: click_param_to_json_schema(param)
                for param in command.params
            },
            "required": [p.name for p in command.params if p.required]
        }
    }
```

`click_param_to_json_schema` must extract: JSON type, `required`, `default`, `choices` (for Click `Choice` params), `help` text, and path/file hints (derivable from `click.Path` / `click.File` param types). This is the complete set needed to render accurate parameter cards in the Skill.

CLI parameter changes are automatically reflected in MCP tool schemas without any manual update. This is the only source of truth for parameter structure.

### 2.3 Execution Model

All phylogenomics commands are treated as long-running. The MCP server uses a **fire-and-forget + output_dir as job handle** model for the majority of tools:

**Fire-and-forget tools** (all CLI commands that produce output):

1. **Pre-launch**: MCP resolves and validates `--output-dir` before starting the subprocess. The resolved absolute path is returned to the caller immediately as the job handle, so the user can track the job even in a new session. If the path is invalid or conflicts with an existing directory (and `--overwrite` is not set), fail before launching.
2. **Launch**: MCP starts the CLI command as a detached background subprocess. On launch, writes a minimal `job.json` to `output_dir` containing `{"pid": ..., "started_at": ..., "command": ...}`. If the subprocess fails to start (e.g. path error, permission denied), captures the launch stderr and returns it immediately as an error.
3. **Track**: `check_status(output_dir)` reads `job.json`, `checkpoint.json` (if present), and `result.json` to determine state. Returns one of: `not_started | running | success | error | unknown`.
4. **Complete**: When `result.json` appears with `status: "success"` or `status: "error"`, the job is done.

```json
// check_status response (in progress)
{
  "status": "running",
  "output_dir": "/abs/path/runs/tree/ml/iqtree",
  "checkpoint": {"completed": 42, "total": 200, "last_updated": "2026-06-27T14:23:00"}
}

// check_status response (complete)
{
  "status": "success",
  "output_dir": "/abs/path/runs/tree/ml/iqtree",
  "result": { ... }
}

// check_status response (process crashed, no result.json)
{
  "status": "unknown",
  "output_dir": "/abs/path/runs/tree/ml/iqtree",
  "message": "Process exited but result.json not found. Check logs/ for tool stderr."
}
```

The `output_dir` (always absolute) is the persistent job handle across sessions. MCP always passes `--output-dir` explicitly to the CLI; it never relies on CLI defaults.

**Synchronous exceptions:** The following tools run synchronously — they return results directly instead of using fire-and-forget:
- `doctor` — runs `phyloai doctor --output-format json` via synchronous subprocess, returns JSON output immediately (fast, seconds).
- `check_status`, `read_result`, `read_report`, `get_command_schema` — read-only file operations, no subprocess involved.
- `report` — runs `phyloai report --run-dir <dir>` synchronously. It is pure Python (no external tools) and completes quickly; the caller typically wants the result immediately. Output goes to `<run_dir>/report/report.json`.

### 2.4 Deployment

**Transport:** stdio (local process). No HTTP server, no port management, no authentication.

**Installation:** `pip install phyloai` includes the MCP server. Users add one entry to their AI client config:

```json
{
  "phyloai": {
    "command": "phyloai",
    "args": ["mcp-server"]
  }
}
```

HTTP transport is deferred. If HPC remote access becomes a requirement, it can be added as an optional transport without changing core server logic.

---

## 3. Skill Design

### 3.1 Scope and Structure

One primary Skill (`phyloai-workflow`) covers the full pipeline. A separate sub-Skill (`phyloai-syserror`) handles systematic error diagnosis orchestration and is deferred to a future release alongside the syserror CLI commands.

**`phyloai-workflow` covers:**
- Session entry and scenario detection
- `phyloai doctor` gate
- All pretree / tree / posttree (topology, dating) guided workflows
- Parameter confirmation cards
- Result interpretation and next-step suggestions
- Session recovery via `report.json`
- Demo mode with bundled dataset

**`phyloai-syserror` (future) covers:**
- Results-driven orchestration: `brlen` screening → conditional `cca` / `sites`
- Iterative human-in-the-loop diagnosis
- Final synthesis across atomic operations

### 3.2 Distribution

The Skill lives inside the phyloai repository at `skills/phyloai-workflow/`. It is installed automatically with `pip install phyloai`. Skill version is coupled to CLI version — parameter semantic annotations and error catalog must be updated in the same PR as any CLI parameter change.

```
skills/
└── phyloai-workflow/
    ├── SKILL.md                        # main skill document
    └── references/
        ├── parameter-annotations.md   # per-command Chinese annotations + recommended values
        ├── error-catalog.md           # exit 1/3 known error patterns + fixes
        ├── dialog-templates.md        # pre-run card and post-run card templates
        ├── demo-data.md               # demo dataset paths and per-step entry points
        └── workflow.md                # per-phase execution reference
```

### 3.3 Session Entry and Scenario Detection

`phyloai doctor` is required before any command that invokes external tools (align, trim, tree inference, etc.), on first run in a session, or when the environment is unknown. Read-only operations (`read_result`, `read_report`, `check_status`) do not require doctor first.

Skill detects scenario from the user's first message without asking a menu of questions:

| User says                                      | Detected scenario  | First action after doctor          |
| ---------------------------------------------- | ------------------ | ---------------------------------- |
| "I have protein sequences, want phylogenomics" | New task           | Ask for seq-dir, begin pretree     |
| "My iqtree finished, what's next"              | Resume task        | Ask for run_dir, call `read_report` |
| "Run pretree stats on this directory"          | Single-step mode   | Show parameter card, confirm, run  |
| Ambiguous                                      | Ask one question   | "Are you starting a new analysis or continuing an existing one?" |

### 3.4 Session Recovery

When resuming an existing run, the Skill attempts recovery in order:

1. **`read_report(run_dir)`** — if `report.json` exists, reconstruct full context (completed steps with key_results, running steps, failed steps, logical next step). Present status summary; user does not need to explain their state.
2. **Fallback — no `report.json`**: ask user "report.json not found. Run `phyloai report` to generate it, or tell me which step to check?" If the user points to a specific step directory, call `read_result(step_dir)` directly.
3. **Running job check**: for any step believed to be in progress, call `check_status(output_dir)` to get current state before suggesting next actions.

### 3.5 Parameter Confirmation Cards

Before executing any command, the Skill presents a parameter card and waits for explicit user confirmation or edits.

**Card structure:**

```
命令: phyloai pretree align
目的: 对原始序列进行多序列比对

参数:
  --seq-dir        ./raw              输入序列目录（FASTA格式）
  --method         linsi              比对算法。linsi: 最精确，适合<500序列；
                                      auto: 自动选择，适合大数据集 [推荐: linsi]
  --threads        8                  并行线程数
  --tool-args      (未设置)           传递给MAFFT的额外参数（可选）

Schema source: runtime CLI (phyloai mcp-server)

确认执行？还是需要调整参数？
```

**Rules:**
- Parameter structure (names, types, required, defaults) loaded from `get_command_schema` at card render time — never from memory.
- Scientific annotations (Chinese descriptions, recommended values, common pitfalls) come from `references/parameter-annotations.md`.
- If a parameter exists in the schema but has no annotation entry, show the CLI `--help` text only, without Chinese annotation.
- Never invent parameter names, aliases, or enum values.
- If the user requests a parameter not in the schema, explicitly mark it as unsupported and suggest the nearest valid parameter.

### 3.6 Language Policy

- **Parameter cards**: fixed format — English parameter name + Chinese annotation + recommended value.
- **Conversation and result interpretation**: follow the user's language. No forced language switching.
- **CLI commands**: always shown in English original form (copy-pasteable).

### 3.7 Result Interpretation

After each command completes, the Skill:

1. Reads `result.json` via `read_result(output_dir)` (not `report.json` — fresher, more precise)
2. Summarizes key results in plain language (drawing on `key_results` fields)
3. Flags any warnings or unexpected values
4. Proposes logical next steps — does not auto-run anything
5. Waits for explicit user selection before proceeding

Example post-align interpretation:
> "比对完成。共处理 200 个基因，成功 198 个，跳过 2 个（序列数不足）。平均比对长度 842 bp，间隙率中位数 12.3%。建议下一步运行 `pretree trim` 去除低质量位点，或先用 `pretree stats` 检查比对质量分布。"

### 3.8 Error Handling

**Exit code 1 / 3 (structural errors — catalog-driven):**

`references/error-catalog.md` covers known patterns:
- Missing required parameter → show which parameter and valid format
- Invalid enum value → show allowed values
- Input directory not found → show path resolution suggestion
- Required tool not installed (exit 3) → guide to `phyloai doctor`, show installation instructions

When a catalog match is found, show a fix card with exact invalid field and corrective action.

**Exit code 2 (external tool failure — AI-driven):**

Pass `tool_stderr` content to AI for free-form diagnosis. Tool error messages are themselves diagnostic — fixed catalog patterns cannot cover the combinatorial failure space of IQ-TREE / PhyloBayes / MAFFT / trimAl. For batch commands, do not dump all locus logs: pass only the stderr of failed loci (identifiable from `result.json data.files[].status`). Cap total stderr passed to AI at a reasonable size (e.g. first 10 failed loci); note truncation if applicable.

### 3.9 Report Integration

`phyloai report` serves two distinct roles in the Skill:

| Role                    | When triggered                                         | Tool called        |
| ----------------------- | ------------------------------------------------------ | ------------------ |
| Session recovery        | Start of new conversation with existing run_dir        | `read_report`      |
| Archive / Methods draft | User explicitly requests report or asks for Methods paragraph | `report` then `read_report` |

The deterministic `methods_text` generated by `templates.py` is a draft. Future capability (deferred): Skill calls `polish_methods(report_json)` to perform AI-assisted review — not merely language polish, but scientific accuracy verification: checking that parameter choices are appropriate for the data (e.g., flagging C20 model on a 500-site matrix), and that analytical logic is coherent. This requires a separate spec when implemented.

### 3.10 Demo Mode

The phyloai package includes a bundled demo dataset:
- **End-to-end dataset**: 20 genes × 6 species, with raw AA (`faa/`) and NT (`fna/`) sequence directories — supports fast supermatrix / supertree walkthroughs.
- **Per-step intermediate data**: pre-aligned, pre-trimmed, concatenated matrices, and gene trees — allows jumping directly to any step without running prior steps.
- **Post-tree demo inputs**: topology-test files (matrix, reference/candidate trees as needed) and dating files (matrix, rooted calibrated tree, and MCMCtree-ready inputs) — allows demonstrating `posttree topology`, `posttree dating hessian`, and `posttree dating mcmc` without asking users to construct specialized files first.

Demo mode entry points:
- After `doctor` passes: "如果需要，我可以用内置示例数据演示某个步骤，然后再用您的数据继续。"
- When blocked by data format / path issues: offer demo as unblocking option
- When user explicitly asks for a walkthrough or teaching run

Demo data paths are resolved dynamically from the installed package location, not hardcoded. After any demo run, the Skill explicitly switches back to user data paths.

Demo runs never write output into the demo data directory. All output goes to a user-specified or auto-generated `runs/` directory.

---

## 4. Maintenance Policy

### 4.1 CLI → MCP (zero maintenance)
MCP tool schemas are generated from the Click command tree at runtime. No action required when CLI parameters change.

### 4.2 CLI → Skill (manual, required)
`references/parameter-annotations.md` contains per-parameter Chinese annotations and recommended values. This file must be updated in the same PR as any CLI parameter change that:
- Renames a parameter
- Removes a parameter
- Adds a new parameter
- Changes a parameter's default value or enum set

This is a required item in the subcommand implementation checklist (see main design spec Section 11).

### 4.3 Version coupling
Skill version is coupled to phyloai package version. Breaking CLI changes require a coordinated Skill update before release.

---

## 5. Future Work (Not in Scope)

| Item                        | Description                                                                          | Trigger condition                        |
| --------------------------- | ------------------------------------------------------------------------------------ | ---------------------------------------- |
| `phyloai-syserror` Skill    | Results-driven syserror orchestration (brlen → cca → sites)                         | After syserror CLI commands implemented  |
| `polish_methods` MCP tool   | AI-assisted report review: scientific accuracy verification + methods text improvement | Separate spec required                   |
| `posttree signal` / `simulate` | Signal distribution and simulation workflows                                     | After CLI commands implemented           |
| HTTP MCP transport          | Remote access for HPC workflows                                                      | When HPC use cases emerge                |

---

## 6. Development Phases (update to main spec)

| Phase | Scope              | Deliverable                                              | Pre-requisites |
| ----- | ------------------ | -------------------------------------------------------- | -------------- |
| 7     | MCP Server         | All CLI tools wrapped; check_status / read_result / read_report / get_command_schema utilities; stub tools for future commands | Phases 1–6     |
| 8     | `phyloai-workflow` Skill | Full guided workflow; parameter cards; result interpretation; session recovery; demo mode; error handling | Phase 7        |
| 9     | `phyloai-syserror` Skill | syserror sub-workflow orchestration                  | syserror CLI + Phase 7 |
| 10    | Report AI review   | `polish_methods` tool + Skill integration                | Separate spec  |
