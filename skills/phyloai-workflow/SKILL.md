# PhyloAI Workflow

Guide users through PhyloAI CLI analyses through the local MCP server.

## Core Rules

- Use `doctor` before commands that invoke external tools, on first run, or when the environment is unknown.
- Read-only tools (`check_status`, `read_result`, `read_report`, `get_command_schema`) do not require `doctor` first.
- Before executing a CLI command, call `get_command_schema`, render a parameter card, and wait for explicit user approval.
- Never invent parameter names, aliases, defaults, or enum values. Unknown parameters block execution.
- After a command completes, summarize `key_results`, warnings, and next steps. Do not auto-run the next step.

## Entry Modes

- New task: ask for input data path, run `doctor` if needed, then start pretree workflow.
- Resume task: call `read_report(run_dir)`; if missing, ask whether to run `report` or inspect a specific step with `read_result`.
- Single-step task: render the parameter card for the requested command and wait for confirmation.

## Language Policy

- Parameter cards use English parameter names with Chinese annotations and recommendations.
- Conversation and interpretation follow the user's language.
- CLI commands are shown in English exactly as executable commands.

## Workflow

- Pretree: `convert -> stats -> align -> trim -> metrics -> filter -> concat`.
- Tree: `tree ml iqtree` or `tree ml fasttree`, then optional `msc` and `cf`.
- Posttree: `topology`, `dating hessian`, `dating mcmc`.
- Report: run `report` only when the user requests a report/methods draft or recovery needs `report.json`.

## Error Handling

- Exit 1/3: use `references/error-catalog.md`, then show a fix card.
- Exit 2: diagnose tool stderr. For batch commands, inspect only failed loci logs, capped at about 10 loci, and state when truncated.

## References

- `references/parameter-annotations.md`
- `references/error-catalog.md`
- `references/dialog-templates.md`
- `references/demo-data.md`
- `references/workflow.md`
