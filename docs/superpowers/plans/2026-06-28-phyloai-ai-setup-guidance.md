# PhyloAI AI Setup Guidance Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make PhyloAI AI sessions Skill-first and move installation/completion details out of README into focused command docs.

**Architecture:** This is documentation and Skill guidance only. `AGENTS.md` controls repository-scoped AI behavior, `docs/commands/installation.md` and `docs/commands/completion.md` hold detailed user docs, and existing README/doctor/AI docs link to those references. No MCP execution code changes are needed.

**Tech Stack:** Markdown documentation, in-repo OpenCode/Codex-style agent instructions, existing `phyloai-workflow` Skill files.

## Global Constraints

- Do not change MCP server execution behavior.
- Do not implement automatic external tool installation.
- Keep README concise; move detailed installation and shell completion guidance into `docs/commands/`.
- `AGENTS.md` only affects AI sessions that read this repository root. Outside the repo, AI clients do not load it automatically. For global Skill-first behavior, OpenCode users can copy the rules into `~/.config/opencode/AGENTS.md` and Claude Code users into `~/.claude/CLAUDE.md`; all users can install `skills/phyloai-workflow` in their AI client.
- Installation docs must prefer official upstream install links over copied long build tutorials.
- Distinguish required tool failures from optional tool warnings in missing-tool guidance.
- Any command with `--overwrite true` requires separate explicit confirmation naming the affected output directory.
- Do not commit changes unless the user explicitly asks for a commit.

---

## File Structure

- Create `AGENTS.md`: repository-scoped AI rules requiring `phyloai-workflow` before execution MCP tools.
- Create `docs/commands/installation.md`: practical installation reference with Python environments, `phyloai doctor`, external tools, and `phyloai run` dependency map.
- Create `docs/commands/completion.md`: detailed Bash/Zsh/Fish shell completion instructions moved from README.
- Modify `README.md`: shorten Installation and Shell Completion sections, link to new command docs, add new command-doc rows.
- Modify `docs/commands/doctor.md`: link missing-tool users to `installation.md` without expanding doctor into an install guide.
- Modify `docs/commands/ai-integration.md`: clarify MCP + Skill pairing and `AGENTS.md` scope.
- Modify `skills/phyloai-workflow/SKILL.md`: strengthen entry/trigger rules for environment and install requests.
- Modify `skills/phyloai-workflow/references/error-catalog.md`: replace one-line exit 3 fix with required/optional missing-tool fix-card guidance.

---

### Task 1: Add Repository Skill-First AI Rules

**Files:**
- Create: `AGENTS.md`

**Interfaces:**
- Consumes: Existing MCP tool names from `docs/commands/ai-integration.md`.
- Produces: Repo-scoped instruction file read by AI agents working in this checkout.

- [ ] **Step 1: Create `AGENTS.md`**

Use this exact content:

```markdown
# PhyloAI Agent Instructions

When a user asks about PhyloAI analysis, command execution, run recovery, `doctor`, missing external tools, installation, environment checks, or external-tool failures, invoke the `phyloai-workflow` Skill before using PhyloAI MCP tools.

Read-only MCP tools may be used directly for inspection:

- `check_status`
- `read_result`
- `read_report`
- `get_command_schema`

Execution MCP tools must go through `phyloai-workflow` parameter review and explicit user approval. Do not guess defaults and launch an execution tool only because the MCP schema is available.

When showing a parameter card, list **every** parameter from `get_command_schema`. Do not filter out parameters. Annotations are decorations, not a display filter. Parameters without annotations must still be shown with their CLI `--help` text.

Scope: this file applies to AI sessions that read the repository root. AI clients that open this repo as their workspace automatically load `AGENTS.md` — no explicit user action is needed. It does NOT activate when the working directory is outside this repo.

For users who install PhyloAI elsewhere and want the same Skill-first behavior globally:
- OpenCode: copy the rules into `~/.config/opencode/AGENTS.md` (global rules, applied to all sessions).
- Claude Code: copy into `~/.claude/CLAUDE.md`.
- Or install `skills/phyloai-workflow` in the AI client's skills directory (see `docs/commands/ai-integration.md`).
```

- [ ] **Step 2: Verify content exists**

Run: `rg -n "every.*parameter|output\.log|natural language|phyloai-workflow|Execution MCP tools|Scope:|opencode/AGENTS.md|CLAUDE.md" AGENTS.md`

Expected: eight matching lines, including the parameter completeness and progress inquiry rules, scope note, and global path options.

---

### Task 2: Add Practical Installation Reference

**Files:**
- Create: `docs/commands/installation.md`

**Interfaces:**
- Consumes: Tool list from `docs/commands/doctor.md` and `docs/superpowers/specs/2026-06-18-phyloai-doctor-design.md`.
- Produces: Single installation target linked by README, doctor docs, AI docs, and Skill fix cards.

- [ ] **Step 1: Create `docs/commands/installation.md`**

Use this structure and keep entries concise:

````markdown
# PhyloAI Installation Guide

## Purpose

This guide explains how to install PhyloAI, make external tools visible to the active shell environment, and verify the setup with `phyloai doctor`.

PhyloAI does not install most third-party phylogenetics tools automatically. Install those tools through your operating system, Conda/Mamba environment, cluster module system, or the upstream project instructions.

## Get The Source

```bash
git clone https://github.com/xtmtd/phyloAI.git
cd phyloAI
```

## Python Environment

Choose one environment style.

### uv

Recommended for local development and quick reproducibility.

```bash
uv venv
source .venv/bin/activate
uv pip install -e .
```

### Conda / Mamba

Recommended when Python packages and bioinformatics command-line tools need to live in the same environment.

```bash
mamba create -n phyloai python=3.11
mamba activate phyloai
pip install -e .
```

Use `conda` instead of `mamba` if that is what your system provides.

### venv

Suitable for a pure Python environment when external tools are installed elsewhere on `PATH`.

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -U pip
python -m pip install -e .
```

## Verify The Environment

```bash
phyloai doctor
phyloai doctor --output-format json
```

`phyloai doctor` checks the current shell environment. Activate the environment you intend to use before running it. It does not install tools or modify `PATH`.

## One-Click Pipeline Dependencies

| Command mode | Main external tools |
|--------------|---------------------|
| `phyloai run --mode supermatrix --speed normal` | MAFFT, trimAl, Julia for TAPER, IQ-TREE3 |
| `phyloai run --mode supermatrix --speed fast` | MAFFT, trimAl, FastTree |
| `phyloai run --mode supertree --speed normal` | MAFFT, trimAl, Julia for TAPER, IQ-TREE3 for gene trees, wASTRAL |
| `phyloai run --mode supertree --speed fast` | MAFFT, trimAl, FastTree, wASTRAL |

## External Tools

### Core Workflow

| Tool | Used by | Install entry point | Detection name | Verify |
|------|---------|---------------------|----------------|--------|
| IQ-TREE3 | `tree ml iqtree`, topology tests, dating Hessian, normal supermatrix `run` | https://github.com/iqtree/iqtree3/releases | `iqtree3` | `phyloai doctor` |
| MAFFT | `pretree align`, `phyloai run` | https://mafft.cbrc.jp/alignment/software/ | `mafft` | `phyloai doctor` |
| trimAl | `pretree trim`, backtranslation, `phyloai run` | https://github.com/inab/trimal | `trimal` | `phyloai doctor` |

IQ-TREE3 and trimAl are planned for future bundling or auto-download in PhyloAI. Until `phyloai doctor` reports them as bundled, treat them as external tools.

### Tree And Posttree Tools

| Tool | Used by | Install entry point | Detection name | Verify |
|------|---------|---------------------|----------------|--------|
| FastTree | `tree ml fasttree`, fast `phyloai run` | http://www.microbesonline.org/fasttree/ | `FastTree` | `phyloai doctor` |
| wASTRAL | `tree msc`, supertree `phyloai run` | https://github.com/chaoszhang/ASTER | `wastral` | `phyloai doctor` |
| MCMCtree / PAML | `posttree dating mcmc` | https://github.com/abacus-gene/paml/releases | `mcmctree` | `phyloai doctor` |

### Bayesian Inference

`phyloai tree bi` needs the PhyloBayes-MPI tool group.

| Tool | Purpose | Install entry point | Detection name | Verify |
|------|---------|---------------------|----------------|--------|
| pb_mpi | MCMC sampler | https://github.com/bayesiancook/pbmpi | `pb_mpi` | `phyloai doctor` |
| bpcomp | Topology convergence | https://github.com/bayesiancook/pbmpi | `bpcomp` | `phyloai doctor` |
| tracecomp | Parameter convergence | https://github.com/bayesiancook/pbmpi | `tracecomp` | `phyloai doctor` |
| mpirun | MPI launcher | https://www.open-mpi.org/ | `mpirun` | `phyloai doctor` |
| readpb_mpi | Optional chain reader | https://github.com/bayesiancook/pbmpi | `readpb_mpi` | `phyloai doctor` |

If the tools are installed outside `PATH`, use `phyloai tree bi --pb-path /path/to/pbmpi/bin`.

### Filtering And Trimming Extras

| Tool | Used by | Install entry point | Detection name | Verify |
|------|---------|---------------------|----------------|--------|
| TreeShrink | `pretree filter treeshrink` | https://github.com/uym2/TreeShrink | `run_treeshrink.py` | `phyloai doctor` |
| MAGUS | `pretree align --method magus` | https://github.com/vlasmirnov/MAGUS | `magus` | `phyloai doctor` |
| ClipKIT | `pretree trim --tool clipkit` | https://github.com/JLSteenwyk/ClipKIT | `clipkit` | `phyloai doctor` |

### Runtime Dependencies

| Tool | Used by | Install entry point | Detection name | Verify |
|------|---------|---------------------|----------------|--------|
| Java | BMGE workflows | https://www.java.com/ | `java` | `phyloai doctor` |
| Julia | TAPER masking | https://julialang.org/downloads/ | `julia` | `phyloai doctor` |

### Bundled Tools

TAPER 1.0.0 (`correction_multi.jl`) and BMGE 1.12 (`BMGE.jar`) are bundled inside the PhyloAI package. They do not need separate installation. If either is missing from `phyloai doctor`, treat it as a PhyloAI packaging or installation problem rather than a missing user-installed dependency.

## Operating System Notes

- macOS: Homebrew or Conda/Mamba are usually the simplest way to make command-line tools visible on `PATH`.
- Linux: use Conda/Mamba, distribution packages, cluster modules, or upstream binaries depending on your environment.
- WSL: install tools inside the Linux distribution, not only on Windows, so `phyloai doctor` can see them from the WSL shell.
````

- [ ] **Step 2: Verify required sections**

Run: `rg -n "uv|One-Click Pipeline Dependencies|PhyloBayes-MPI|Bundled Tools|IQ-TREE3 and trimAl" docs/commands/installation.md`

Expected: all five patterns match.

---

### Task 3: Move Shell Completion Details Into Command Doc

**Files:**
- Create: `docs/commands/completion.md`

**Interfaces:**
- Consumes: Current README Shell Completion section.
- Produces: Detailed shell completion doc for README to link.

- [ ] **Step 1: Create `docs/commands/completion.md`**

Use this exact content:

````markdown
# phyloai completion

## Purpose

`phyloai completion <shell>` generates static shell completion scripts for Bash, Zsh, and Fish.

Generate the script once from an environment where `phyloai` is installed, save it to a persistent file, and configure your shell to load that saved script.

Do not run `phyloai completion ...` dynamically from `.bashrc`, `.zshrc`, `config.fish`, or other shell startup files.

## Usage

```bash
phyloai completion bash
phyloai completion zsh
phyloai completion fish
```

## Bash

```bash
mkdir -p ~/.config/phyloai/completion
phyloai completion bash > ~/.config/phyloai/completion/phyloai.bash
```

Add this line to `~/.bashrc`:

```bash
source ~/.config/phyloai/completion/phyloai.bash
```

If you only run the `source` command manually in the current terminal, completion only works for that shell session.

## Zsh

```bash
mkdir -p ~/.config/phyloai/completion
phyloai completion zsh > ~/.config/phyloai/completion/phyloai.zsh
```

Add this line to `~/.zshrc`:

```bash
source ~/.config/phyloai/completion/phyloai.zsh
```

If you only run the `source` command manually in the current terminal, completion only works for that shell session.

## Fish

```bash
mkdir -p ~/.config/fish/completions
phyloai completion fish > ~/.config/fish/completions/phyloai.fish
```

Fish loads completion files from `~/.config/fish/completions/` automatically in new shells, so no extra `source` line is required.
````

- [ ] **Step 2: Verify completion commands**

Run: `rg -n "completion bash|completion zsh|completion fish|Do not run" docs/commands/completion.md`

Expected: all four patterns match.

---

### Task 4: Link README And Existing Command Docs

**Files:**
- Modify: `README.md`
- Modify: `docs/commands/doctor.md`
- Modify: `docs/commands/ai-integration.md`

**Interfaces:**
- Consumes: `docs/commands/installation.md`, `docs/commands/completion.md`, `AGENTS.md`.
- Produces: Public entry points that route users to the new references.

- [ ] **Step 1: Shorten README Installation section**

Replace the current `## Installation` section with:

````markdown
## Installation

```bash
git clone https://github.com/xtmtd/phyloAI.git
cd phyloAI
pip install -e .
```

PhyloAI bundles TAPER 1.0.0 (`correction_multi.jl`) and BMGE 1.12 (`BMGE.jar`). Other external tools should be installed for your operating system and workflow, then verified with:

```bash
phyloai doctor
```

See [docs/commands/installation.md](docs/commands/installation.md) for Python environment options, external tool groups, and operating-system notes.
````

- [ ] **Step 2: Shorten README Shell Completion section**

Replace the current detailed `## Shell Completion` section with:

````markdown
## Shell Completion

PhyloAI can generate static completion scripts for Bash, Zsh, and Fish:

```bash
phyloai completion bash
phyloai completion zsh
phyloai completion fish
```

Generate the script once and configure your shell to load the saved file. See [docs/commands/completion.md](docs/commands/completion.md) for Bash, Zsh, and Fish setup examples.
````

- [ ] **Step 3: Add command table rows to README**

Insert these rows near the top of the README command table, after `phyloai doctor`:

```markdown
| Installation | Set up Python environments and external tools, then verify with `phyloai doctor`. | [docs/commands/installation.md](docs/commands/installation.md) |
| `phyloai completion` | Generate static Bash, Zsh, or Fish shell completion scripts. | [docs/commands/completion.md](docs/commands/completion.md) |
```

- [ ] **Step 4: Add doctor doc link**

In `docs/commands/doctor.md`, add this paragraph after the paragraph ending with "should be installed by the user for their operating system, package manager, cluster, or Conda environment.":

```markdown
For practical setup guidance, see [installation.md](installation.md). It lists Python environment options, external tool groups, `phyloai run` dependency modes, and operating-system notes.
```

- [ ] **Step 5: Add AI integration warning**

In `docs/commands/ai-integration.md`, add this paragraph after line 30 (`Add the MCP server to your client config, then make the Skill discoverable.`):

```markdown
Do not enable only the MCP server and skip the Skill. MCP exposes execution tools, but `phyloai-workflow` is the guided layer that runs `doctor`, renders parameter cards, waits for approval, interprets results, and diagnoses missing tools.

`AGENTS.md` is an AI-agent convention file at the project root. AI clients (OpenCode, Claude Code, Codex) that open this repo as their workspace automatically load it and follow its rules — no explicit user action is needed. It does NOT activate when the working directory is outside this repo.

For Skill-first behavior outside the repo:
- **OpenCode:** copy the rules from `AGENTS.md` (or the whole file) into `~/.config/opencode/AGENTS.md` for global application across all sessions.
- **Claude Code:** copy into `~/.claude/CLAUDE.md`.
- Or install `skills/phyloai-workflow` in your AI client's skills directory (see the sections below).
```

- [ ] **Step 6: Verify links and README size reduction**

Run: `rg -n "installation.md|completion.md|AGENTS.md|Do not enable only the MCP server" README.md docs/commands/doctor.md docs/commands/ai-integration.md`

Expected: matches in all three files.

---

### Task 5: Strengthen Skill Entry, Overwrite Safety, And Missing-Tool Guidance

**Files:**
- Modify: `skills/phyloai-workflow/SKILL.md`
- Modify: `skills/phyloai-workflow/references/error-catalog.md`

**Interfaces:**
- Consumes: `docs/commands/installation.md` anchors and current Skill error-handling rules.
- Produces: Skill instructions that catch environment/install requests, require destructive overwrite confirmation, render complete parameter cards (all parameters from schema), and render concise missing-tool fix cards.

- [ ] **Step 1: Update Skill core rules**

In `skills/phyloai-workflow/SKILL.md`, add this bullet after the existing `doctor` core rule:

```markdown
- Use this Skill for environment and installation requests too, including `doctor failed`, `missing pb_mpi`, `install iqtree`, `缺少 MAFFT`, `环境检查失败`, and similar external-tool setup questions.
```

- [ ] **Step 2: Add parameter completeness and overwrite safety core rules**

In `skills/phyloai-workflow/SKILL.md`, add these bullets after the parameter-card approval core rule:

```markdown
- Before executing a CLI command, call `get_command_schema`, render a parameter card that lists **every** parameter from the schema, and wait for explicit user approval. Do not filter out parameters — annotations in `references/parameter-annotations.md` are decorations, not a display filter. Parameters without annotations must still be shown with their CLI `--help` text.
- Treat `--overwrite` as destructive. When the target `--output-dir` already exists and the user has not explicitly requested overwrite, prefer suggesting a new `--output-dir` or `--resume` when available before offering `--overwrite`. If a parameter card sets `--overwrite true`, ask for separate explicit confirmation naming the affected `--output-dir`; general command approval is not enough.
```

- [ ] **Step 2b: Add progress inquiry core rule**

In `skills/phyloai-workflow/SKILL.md`, add this bullet after the result-summary core rule:

```markdown
- When a user asks about progress for an MCP-launched running job, read the tail of `<output_dir>/output.log` alongside `check_status`. Summarize the progress in natural language. Do not tell the user to run `tail -f` unless they explicitly ask to see raw output. All MCP-launched PhyloAI commands write runtime output to `output.log`.
```

- [ ] **Step 3: Update Skill references list**

In `skills/phyloai-workflow/SKILL.md`, add this reference to the `## References` list:

```markdown
- `docs/commands/installation.md` for external-tool setup guidance
```

- [ ] **Step 3b: Update dialog templates**

Replace the `## Parameter Card` section in `skills/phyloai-workflow/references/dialog-templates.md` with the multi-parameter template. Add a new `## Progress Inquiry` section.

```text
## Parameter Card

命令: phyloai <command>
目的: <one-line purpose>

参数:
  --paramA      <value>    <中文说明> [推荐: <value>]
  --paramB      <value>    <中文说明 or --help text>
  ...

Schema source: runtime CLI via get_command_schema

确认执行？还是需要调整参数？

Rules:
- List every parameter from get_command_schema. Do not omit any.
- Parameters with annotations get Chinese descriptions.
- Parameters without annotations use CLI --help text verbatim.
```

```text
## Progress Inquiry

When user asks about a running job:
1. Call check_status for state.
2. Read tail of <output_dir>/output.log (~30 lines).
3. Summarize in natural language.
4. Do not tell user to run tail -f unless they explicitly ask.
```

- [ ] **Step 4: Expand overwrite conflict catalog guidance**

Replace the `## Exit 1: Output Directory Exists` section in `skills/phyloai-workflow/references/error-catalog.md` with:

````markdown
## Exit 1: Output Directory Exists
Pattern: `Output directory exists and is not empty`
Fix: Prefer a new `--output-dir` or `--resume` when available. Use `--overwrite` only after separate explicit confirmation.

Overwrite confirmation template:

```text
`--overwrite` 会删除或替换已有输出目录:
<output-dir>

请明确确认是否覆盖这个目录。未确认前不要执行。
```
````

- [ ] **Step 5: Replace exit 3 catalog guidance**

Replace the two exit 3 sections in `skills/phyloai-workflow/references/error-catalog.md` with:

````markdown
## Exit 3: Required Tool Missing
Pattern: `not installed`, `not detectable`, `not found`, `Missing required tool`
Fix: Show a concise missing-tool fix card. Include the missing tool name, whether it is required or optional for the requested command, the affected PhyloAI command, why the tool is needed, and a link to `docs/commands/installation.md`.

Template:

```text
缺少工具: <tool>
状态: required | optional
影响命令: <phyloai command>
为什么需要: <short purpose>
下一步:
  1. 查看 docs/commands/installation.md#<tool-or-group-anchor>
  2. 安装后运行 phyloai doctor
  3. 如果已安装但未检测到，检查 PATH 或使用该命令的显式工具路径参数
```

Notes:
- Missing required tools block the requested command and may make `phyloai doctor` exit 3.
- Missing optional tools only make the dependent module unavailable; `phyloai doctor` can still exit 0 with warnings.
- For `phyloai tree bi`, treat `pb_mpi`, `bpcomp`, `tracecomp`, and `mpirun` as the required PhyloBayes-MPI group; `readpb_mpi` is optional.

## Exit 3: Runtime Missing
Pattern: `Java`, `Julia`, `MPI`
Fix: Show the same missing-tool fix card. Link Java and Julia to the Runtime Dependencies section of `docs/commands/installation.md`; link MPI-related failures to the Bayesian Inference section.
````

- [ ] **Step 6: Verify Skill trigger, overwrite warning, parameter completeness, and progress rules**

Run: `rg -n "every.*parameter|output\.log|tail.*output\.log|natural language|missing pb_mpi|installation.md|状态: required|readpb_mpi|overwrite.*destructive|覆盖这个目录" skills/phyloai-workflow/SKILL.md skills/phyloai-workflow/references/error-catalog.md skills/phyloai-workflow/references/dialog-templates.md`

Expected: all ten patterns match across the three files.

---

### Task 6: Final Documentation Verification

**Files:**
- Verify only; modify files only if checks reveal broken links or missing required text.

**Interfaces:**
- Consumes: All files changed in Tasks 1-5.
- Produces: Evidence that docs are internally linked and the design requirements are covered.

- [ ] **Step 1: Check all new docs are reachable from README or command docs**

Run: `rg -n "docs/commands/installation.md|docs/commands/completion.md|installation.md\)|completion.md\)" README.md docs/commands/*.md skills/phyloai-workflow/**/*.md`

Expected: links from README, doctor docs, AI integration docs or Skill files.

- [ ] **Step 2: Check AGENTS scope is explicit**

Run: `rg -n "Scope:|opencode/AGENTS.md|CLAUDE.md|phyloai-workflow" AGENTS.md docs/commands/ai-integration.md`

Expected: `AGENTS.md` contains the scope limitation with global path options; AI integration mentions global OpenCode/Claude Code AGENTS.md paths.

- [ ] **Step 3: Check no stale README-only shell completion block remains**

Run: `rg -n "mkdir -p ~/.config/phyloai/completion|~/.bashrc|~/.zshrc|~/.config/fish/completions" README.md docs/commands/completion.md`

Expected: matches only in `docs/commands/completion.md`, not README.

- [ ] **Step 4: Check missing-tool guidance distinguishes required and optional**

Run: `rg -n "required|optional|doctor.*exit 3|warnings" skills/phyloai-workflow/references/error-catalog.md docs/commands/installation.md docs/commands/doctor.md`

Expected: error catalog distinguishes required vs optional; doctor docs still explain optional missing tools are warnings.

- [ ] **Step 5: Check overwrite confirmation guidance exists**

Run: `rg -n "overwrite|覆盖|separate explicit confirmation|output-dir" skills/phyloai-workflow/SKILL.md skills/phyloai-workflow/references/error-catalog.md`

Expected: Skill core rules require separate overwrite confirmation; error catalog includes the Chinese overwrite confirmation template.

- [ ] **Step 6: Review diff manually**

Run: `git diff -- README.md AGENTS.md docs/commands/installation.md docs/commands/completion.md docs/commands/doctor.md docs/commands/ai-integration.md skills/phyloai-workflow/SKILL.md skills/phyloai-workflow/references/error-catalog.md`

Expected: documentation-only changes; no MCP server or CLI execution code changed.

---

## Self-Review Notes

- Spec coverage: Task 1 covers Skill-first `AGENTS.md` with scope limitation and global `~/.config/opencode/AGENTS.md` / `~/.claude/CLAUDE.md` options; Task 2 covers installation with `git clone`, bundled/planned notes, uv, and `run` dependencies; Task 3 covers completion docs; Task 4 covers README/doctor/AI docs with AGENTS.md mechanism explanation; Task 5 covers Skill triggers, parameter completeness (all params from schema), overwrite safety, progress inquiry via `output.log`, dialog templates, and fix cards; Task 6 verifies links, scope with global paths, parameter completeness rules, `output.log` progress rules, and destructive overwrite confirmation.
- Placeholder scan: no task uses open-ended placeholders such as "fill in later" or unspecified tests.
- Type consistency: no code interfaces are introduced; all referenced paths and MCP tool names match existing docs.
