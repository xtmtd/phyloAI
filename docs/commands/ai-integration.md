# PhyloAI AI Integration

[English](ai-integration.md) | [中文](ai-integration.zh.md)


PhyloAI provides an AI interaction layer through **MCP Server** + **Skill**, enabling conversational phylogenomics analysis without memorizing CLI parameters.

```
用户  ←→  AI 客户端 (Skill 驱动对话)  ←→  MCP Server (执行桥)  ←→  phyloai CLI  ←→  文件系统
```

| Layer | Role |
|-------|------|
| **MCP Server** | Exposes analysis, report, doctor, and run CLI commands as callable tools with auto-generated schemas. Handles job launch, status tracking, and result reading. |
| **Skill** | Guides the conversation: parameter confirmation cards, result interpretation, next-step suggestions, session recovery, and error diagnosis. |

## Quick Start

### 1. Install phyloai

```bash
pip install -e .
```

Verify the MCP server starts correctly:

```bash
phyloai mcp-server --help
```

### 2. Configure your AI client

Add the MCP server to your client config, then make the Skill discoverable.

Do not enable only the MCP server and skip the Skill. MCP exposes execution tools, but `phyloai-workflow` is the guided layer that runs `doctor`, renders parameter cards, waits for approval, interprets results, and diagnoses missing tools.

`AGENTS.md` is an AI-agent convention file at the project root. AI clients (OpenCode, Claude Code, Codex) that open this repo as their workspace automatically load it and follow its rules — no explicit user action is needed. It does NOT activate when the working directory is outside this repo.

For Skill-first behavior outside the repo:
- **OpenCode:** copy the rules from `AGENTS.md` (or the whole file) into `~/.config/opencode/AGENTS.md` for global application across all sessions.
- **Claude Code:** copy into `~/.claude/CLAUDE.md`.
- Or install `skills/phyloai-workflow` in your AI client's skills directory (see the sections below).

---

## OpenCode (Primary)

### MCP Server

Add to `~/.config/opencode/opencode.jsonc`:

```jsonc
"mcp": {
  // ... existing servers ...
  "phyloai": {
    "type": "local",
    "command": ["phyloai", "mcp-server"],
    "enabled": true
  }
}
```

### Skill

Copy (or symlink) the phyloai-workflow skill into OpenCode's skills directory:

```bash
cp -r skills/phyloai-workflow ~/.config/opencode/skills/phyloai-workflow
```

Or symlink (keeps in sync with repo updates):

```bash
ln -s "$(pwd)/skills/phyloai-workflow" ~/.config/opencode/skills/phyloai-workflow
```

If you always work inside the phyloai repo, you can skip copying — OpenCode discovers the local `skills/` directory automatically when the workspace is the phyloai project root.

### Start a conversation

Restart OpenCode, then try any of these:

| You say | What happens |
|---------|--------------|
| "I have protein sequences, run phylogenomics" | Skill detects new task → runs `doctor` → asks for `seq-dir` → guides through pretree workflow |
| "My iqtree just finished, what's next?" | Skill calls `read_report` or `read_result` → summarizes completed steps → suggests next logical command |
| "Run pretree stats on ./raw" | Skill calls `get_command_schema("pretree_stats")` → renders parameter card with Chinese annotations → waits for your confirmation |
| "帮我做系统发育分析" | Same as above, conversation follows your language |
| "Show me a demo first" | Skill uses bundled demo data (20 genes × 6 species) → runs steps into a user/auto-generated run directory |

Key behaviors:
- **Every command requires explicit approval** — the Skill shows a parameter card and waits for your OK before executing.
- **Long-running jobs are fire-and-forget** — you get back an `output_dir` handle; call `check_status` anytime to see progress.
- **Session recovery** — start a new conversation, say "resume my analysis in runs/run/", the Skill reads `report.json` and reconstructs full context.

---

## Claude Code

### MCP Server

Add to `~/.claude/claude_desktop_config.json`:

```json
{
  "mcpServers": {
    "phyloai": {
      "command": "phyloai",
      "args": ["mcp-server"]
    }
  }
}
```

Or via `.mcp.json` in the project root for project-scoped servers.

### Skill

Copy to Claude Code's skills directory:

```bash
cp -r skills/phyloai-workflow ~/.claude/skills/phyloai-workflow
```

### Start a conversation

Same interaction patterns as OpenCode — the Skill defines the workflow regardless of which client you use.

---

## Codex (OpenAI Codex CLI)

### MCP Server

Add to `~/.codex/config.toml`:

```toml
[mcp_servers.phyloai]
command = "phyloai"
args = ["mcp-server"]
```

### Skill

Copy to Codex skills directory:

```bash
cp -r skills/phyloai-workflow ~/.codex/skills/phyloai-workflow
```

---

## Built-in Tools

The MCP server registers one tool per CLI subcommand plus four utility tools. All tool schemas are generated dynamically from the Click command tree — no manual sync needed.

### CLI tools

| Tool | CLI equivalent | Execution |
|------|---------------|
| `doctor` | `phyloai doctor --output-format json` | synchronous |
| `pretree_convert` | `phyloai pretree convert` | fire-and-forget |
| `pretree_stats` | `phyloai pretree stats` | fire-and-forget |
| `pretree_align` | `phyloai pretree align` | fire-and-forget |
| `pretree_trim` | `phyloai pretree trim` | fire-and-forget |
| `pretree_metrics` | `phyloai pretree metrics` | fire-and-forget |
| `pretree_filter_taper` | `phyloai pretree filter taper` | fire-and-forget |
| `pretree_filter_treeshrink` | `phyloai pretree filter treeshrink` | fire-and-forget |
| `pretree_filter_metrics` | `phyloai pretree filter metrics` | fire-and-forget |
| `pretree_filter_symtest` | `phyloai pretree filter symtest` | fire-and-forget |
| `pretree_filter_cluster` | `phyloai pretree filter cluster` | fire-and-forget |
| `pretree_concat` | `phyloai pretree concat` | fire-and-forget |
| `tree_ml_fasttree` | `phyloai tree ml fasttree` | fire-and-forget |
| `tree_ml_iqtree` | `phyloai tree ml iqtree` | fire-and-forget |
| `tree_bi_pb` | `phyloai tree bi pb` | fire-and-forget |
| `tree_bi_bpcomp` | `phyloai tree bi bpcomp` | fire-and-forget |
| `tree_bi_tracecomp` | `phyloai tree bi tracecomp` | fire-and-forget |
| `tree_bi_readpb` | `phyloai tree bi readpb` | fire-and-forget |
| `tree_msc` | `phyloai tree msc` | fire-and-forget |
| `tree_cf` | `phyloai tree cf` | fire-and-forget |
| `posttree_topology` | `phyloai posttree topology` | fire-and-forget |
| `posttree_dating_hessian` | `phyloai posttree dating hessian` | fire-and-forget |
| `posttree_dating_mcmc` | `phyloai posttree dating mcmc` | fire-and-forget |
| `posttree_syserror_brlen` | `phyloai posttree syserror brlen` | fire-and-forget |
| `posttree_syserror_brlen_label_nodes` | `phyloai posttree syserror brlen label-nodes` | fire-and-forget |
| `posttree_syserror_cca` | `phyloai posttree syserror cca` | fire-and-forget |
| `report` | `phyloai report` | synchronous |
| `run` | `phyloai run` | fire-and-forget |

Stub tools (return "not yet available"): `posttree_simulate`, `posttree_syserror_sites`.

### Utility tools (read-only, synchronous)

| Tool | Description |
|------|-------------|
| `check_status` | Check job status by `output_dir`. Returns `not_started` / `running` / `success` / `error` / `unknown`. |
| `read_result` | Read `result.json` from a step output directory. |
| `read_report` | Read `report.json` from a run directory (looks under `<run_dir>/report/`). |
| `get_command_schema` | Get runtime parameter schema (names, types, defaults, choices, help) for any tool. |

## Skill Reference Files

The Skill at `skills/phyloai-workflow/` includes:

| File | Content |
|------|---------|
| `SKILL.md` | Core rules, entry modes, workflow phases, error handling strategy |
| `references/parameter-annotations.md` | Chinese scientific annotations for commonly-used parameters; runtime schema from `get_command_schema` is the authoritative source; operational and visual parameters omitted |
| `references/error-catalog.md` | Exit 1/3 known error patterns and fix cards |
| `references/dialog-templates.md` | Parameter cards, result cards, recovery cards, demo prompts |
| `references/demo-data.md` | Bundled demo dataset paths and per-step entry points |
| `references/workflow.md` | Per-phase execution reference with input/output/check/next-step guides |

## Architecture Notes

- **Transport**: stdio (local process). No HTTP server, no ports, no authentication.
- **Schema source**: All MCP tool parameter schemas are generated at server startup by introspecting the Click command tree. CLI parameter changes are automatically reflected — zero maintenance.
- **Job handle**: All fire-and-forget commands return `output_dir` (absolute path) as the persistent job handle. Track across sessions with `check_status`.
- **Result.json**: Every CLI command writes structured JSON results. The Skill reads these for interpretation; the MCP server reads them for tool responses.
- **Skill versioning**: The Skill lives in-repo and is version-coupled to the CLI. Parameter annotation updates and CLI changes ship together.
