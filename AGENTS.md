# PhyloAI Agent Instructions

When a user asks about PhyloAI analysis, command execution, run recovery, `doctor`, missing external tools, installation, environment checks, or external-tool failures, invoke the `phyloai-workflow` Skill before using PhyloAI MCP tools.

Read-only MCP tools may be used directly for inspection:

- `check_status`
- `read_result`
- `read_report`
- `get_command_schema`

Execution MCP tools must go through `phyloai-workflow` parameter review and explicit user approval. Do not guess defaults and launch an execution tool only because the MCP schema is available.


