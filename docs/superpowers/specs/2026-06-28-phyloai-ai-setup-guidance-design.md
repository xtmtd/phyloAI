# PhyloAI AI Setup Guidance Design

**Date:** 2026-06-28
**Status:** Draft for review

## Goal

Improve PhyloAI's AI-assisted user experience in two places:

1. When an AI client has the PhyloAI MCP server enabled, PhyloAI analysis should still be guided by `phyloai-workflow` instead of directly launching execution tools with guessed defaults.
2. When `doctor` or a workflow step reports missing external tools, the Skill should give concise install guidance and point users to a maintained installation command document.

## Scope

In scope:

- Add project-level AI instructions in `AGENTS.md`.
- Strengthen `phyloai-workflow` entry rules and missing-tool guidance.
- Require explicit destructive-action confirmation whenever a command would run with `--overwrite`.
- Add `docs/commands/installation.md` as the detailed installation reference.
- Add `docs/commands/completion.md` as the shell completion reference.
- Keep README installation content short and link to the detailed reference.
- Keep README shell completion content short and link to the detailed reference.
- Update `doctor` and AI integration docs to point users to the installation reference and the Skill-first MCP usage pattern.

Out of scope:

- Automatic tool installation.
- Changes to MCP execution behavior.
- Full, tool-by-tool compile tutorials that duplicate upstream documentation.

## 1. Skill-First AI Usage

Add `AGENTS.md` at the project root with a short, explicit rule set for AI agents working inside this repository.

Required behavior:

- If the user mentions PhyloAI analysis, running commands, resuming runs, `doctor`, missing tools, external-tool failures, or environment checks, invoke `phyloai-workflow` before using PhyloAI MCP tools.
- Direct use of read-only MCP tools is allowed for inspection: `check_status`, `read_result`, `read_report`, and `get_command_schema`.
- Execution MCP tools must go through the Skill's parameter card and explicit user approval.
- Do not guess defaults and run a workflow command just because the MCP schema is available.

This applies to AI sessions that read the repository root. External users who install PhyloAI with `pip` will not automatically load this project-level file; for those users, Skill-first behavior depends on the client exposing `phyloai-workflow`, as documented in `docs/commands/ai-integration.md`. For global application, OpenCode users can copy the rules into `~/.config/opencode/AGENTS.md` and Claude Code users into `~/.claude/CLAUDE.md`.

## 2. Installation Documentation

Add `docs/commands/installation.md` as the main installation reference.

README keeps only the short path:

- Python installation command.
- Note that TAPER and BMGE are bundled.
- `phyloai doctor` as the verification command.
- Link to `docs/commands/installation.md` for external tools and OS-specific notes.

The installation document should use a practical, maintainable structure:

### Get The Source

Show `git clone` from the PhyloAI repo, with the repository URL.

### Python Environments

Cover three options:

- `uv`: recommended for local development and quick reproducibility.
- `conda` / `mamba`: recommended when mixing Python with bioinformatics tools.
- `venv`: suitable for pure Python environments.

### Verification

Show:

```bash
phyloai doctor
phyloai doctor --output-format json
```

Explain that `doctor` checks the current shell environment and does not install anything.

### External Tool Groups

Group tools by workflow role:

- Core workflow: `iqtree3`, `mafft`, `trimal`.
- Tree and posttree: `FastTree`, `wastral`, `mcmctree`.
- Bayesian inference: `pb_mpi`, `bpcomp`, `tracecomp`, `mpirun`, optional `readpb_mpi`.
- Filtering and trimming extras: `run_treeshrink.py`, `magus`, `clipkit`.
- Runtime dependencies: `java`, `julia`.
- Bundled tools: TAPER and BMGE, marked as no separate install needed. Also note that IQ-TREE3 and trimAl are planned for bundling or auto-download later, but currently should still be treated as external unless `doctor` reports them as bundled.

Each tool entry should use the same compact fields:

- Purpose.
- PhyloAI commands that need it.
- Recommended installation entry point, preferably upstream official docs or project page.
- macOS / Linux / WSL notes when useful.
- Detection name used by PhyloAI.
- Verification command, usually `phyloai doctor`.

The document should prefer upstream links over copied build instructions to avoid stale installation advice.

For `phyloai run`, include a compact dependency map:

- `--speed normal --mode supermatrix`: MAFFT, trimAl, Julia for TAPER, IQ-TREE3.
- `--speed fast --mode supermatrix`: MAFFT, trimAl, FastTree.
- `--speed normal --mode supertree`: MAFFT, trimAl, Julia for TAPER, IQ-TREE3 for gene trees, wASTRAL.
- `--speed fast --mode supertree`: MAFFT, trimAl, FastTree, wASTRAL.

## 3. Shell Completion Documentation

Add `docs/commands/completion.md` and move the current README shell completion details there.

The completion document should cover:

- `phyloai completion bash`
- `phyloai completion zsh`
- `phyloai completion fish`
- Generate once, save to a persistent file, and source/load the saved script.
- Do not run `phyloai completion ...` dynamically from shell startup files.

README should keep only a short shell completion note and link to `docs/commands/completion.md`.

## 4. Doctor And AI Integration Docs

Update `docs/commands/doctor.md` only lightly:

- Keep `doctor` focused on checking tool availability.
- Add a link to `docs/commands/installation.md` for missing tools.

Update `docs/commands/ai-integration.md` with a short warning:

- MCP exposes execution tools, but the Skill is the guided workflow layer.
- Users should install or expose both MCP and `phyloai-workflow`.
- Inside this repository, AI agents should follow `AGENTS.md` and use the Skill first.

## 5. Skill Missing-Tool Guidance

Strengthen `skills/phyloai-workflow/SKILL.md` and `references/error-catalog.md`.

The Skill must treat `--overwrite` as destructive. Any parameter card with `--overwrite true` must ask for a separate confirmation that names the output directory that may be deleted or replaced. General command approval is not enough for `--overwrite`.

Required behavior:

- If the user did not explicitly ask to overwrite, prefer a new `--output-dir` or `--resume` when available.
- If the user or workflow requests `--overwrite`, show a warning before execution.
- The warning must include the affected `--output-dir` and ask the user to confirm overwrite explicitly.
- Do not execute until the user confirms the destructive overwrite action.

Skill trigger coverage should include:

- `doctor failed`
- `missing pb_mpi`
- `install iqtree`
- `缺少 MAFFT`
- `环境检查失败`
- Any user request about installing or fixing PhyloAI external tools

For exit code 3 or `doctor` missing-tool output, the Skill should show a short fix card rather than a long tutorial. This must be implemented in `references/error-catalog.md`, not only described in this design.

The card should distinguish required and optional tools:

- Missing required tools block the requested command and may make `doctor` exit 3.
- Missing optional tools make only the dependent module unavailable; `doctor` can still exit 0 with warnings.

Example card:

```text
缺少工具: pb_mpi
影响命令: phyloai tree bi
为什么需要: PhyloBayes-MPI MCMC sampler
下一步:
  1. 查看 docs/commands/installation.md#phylobayes-mpi
  2. 安装后运行 phyloai doctor
  3. 如果已安装但未检测到，检查 PATH 或使用 --pb-path
```

For grouped tools, the card should name the group. For example, `phyloai tree bi` needs `pb_mpi`, `bpcomp`, `tracecomp`, and `mpirun`; `readpb_mpi` is optional.

The Skill should not attempt to install tools automatically. It should diagnose, link to the installation reference, and provide the shortest verification command.

## Success Criteria

- In-repo AI sessions have a clear rule requiring `phyloai-workflow` before execution MCP tools.
- The docs state that `AGENTS.md` only affects repository-scoped AI sessions, not arbitrary external client installs.
- README remains concise.
- Shell completion details live in `docs/commands/completion.md` instead of README.
- Users can find a single detailed installation reference from README, doctor docs, AI integration docs, and Skill fix cards.
- The installation reference includes `phyloai run --mode/--speed` dependency notes.
- Missing-tool fix cards distinguish required tool failures from optional tool warnings.
- Any command with `--overwrite true` requires separate explicit confirmation naming the affected output directory.
- Missing `pb_mpi` and similar failures result in actionable guidance instead of only "tool not found".
- No MCP execution-layer changes are required.
