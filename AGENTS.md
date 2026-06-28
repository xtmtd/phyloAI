# PhyloAI Agent Instructions

When a user asks about PhyloAI analysis, command execution, run recovery, `doctor`, missing external tools, installation, environment checks, or external-tool failures, invoke the `phyloai-workflow` Skill before using PhyloAI MCP tools.

Read-only MCP tools may be used directly for inspection:

- `check_status`
- `read_result`
- `read_report`
- `get_command_schema`

Execution MCP tools must go through `phyloai-workflow` parameter review and explicit user approval. Do not guess defaults and launch an execution tool only because the MCP schema is available.

Scope: this file applies to AI sessions that read the repository root. AI clients that open this repo as their workspace automatically load `AGENTS.md` — no explicit user action is needed. It does NOT activate when the working directory is outside this repo.

For users who install PhyloAI elsewhere and want the same Skill-first behavior globally:
- OpenCode: copy the rules into `~/.config/opencode/AGENTS.md` (global rules, applied to all sessions).
- Claude Code: copy into `~/.claude/CLAUDE.md`.
- Or install `skills/phyloai-workflow` in the AI client's skills directory (see `docs/commands/ai-integration.md`).
