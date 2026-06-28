# PhyloAI Workflow

Guide users through PhyloAI CLI analyses through the local MCP server.

## Core Rules

- Use `doctor` before commands that invoke external tools, on first run, or when the environment is unknown.
- Use this Skill for environment and installation requests too, including `doctor failed`, `missing pb_mpi`, `install iqtree`, `缺少 MAFFT`, `环境检查失败`, and similar external-tool setup questions.
- Read-only tools (`check_status`, `read_result`, `read_report`, `get_command_schema`) do not require `doctor` first.
- Before executing a CLI command, call `get_command_schema`, render a parameter card, and wait for explicit user approval.
- Treat `--overwrite` as destructive. When the target `--output-dir` already exists and the user has not explicitly requested overwrite, prefer suggesting a new `--output-dir` or `--resume` when available before offering `--overwrite`. If a parameter card sets `--overwrite true`, ask for separate explicit confirmation naming the affected `--output-dir`; general command approval is not enough.
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

- Pretree: `convert -> align -> trim -> metrics / filter -> concat` (supermatrix) or `... -> gene trees` (supertree). `stats` inspects results at any step.
- Tree: `tree ml iqtree` + `tree msc` as primary, `tree ml fasttree` for fast exploration, `tree bi` optional, `cf` on species trees.
- Posttree: `topology`, `dating hessian`, `dating mcmc`.
- Report: run `report` only when the user requests a report/methods draft or recovery needs `report.json`.

## Demo Data

- Bundled demo data contains 20 genes and 6 species, with AA (`faa/`) and NT (`fna/`) raw inputs.
- Per-step demo directories include aligned, trimmed, concatenated, gene-tree, topology-test, and dating entry points when present.
- Demo runs must write new outputs to a user run directory, not back into `phyloai/demo_data/`.

## Error Handling

- Exit 1/3: use `references/error-catalog.md`, then show a fix card.
- Exit 2: diagnose tool stderr. For batch commands, inspect only failed loci logs, capped at about 10 loci, and state when truncated.

## References

- `references/parameter-annotations.md`
- `references/error-catalog.md`
- `references/dialog-templates.md`
- `references/demo-data.md`
- `references/workflow.md`
- `docs/commands/installation.md` for external-tool setup guidance
