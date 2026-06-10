# CLI Shell Completion Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a first-party `phyloai completion` command that emits static Bash, Zsh, and Fish completion scripts for the full CLI without adding any new Python dependencies.

**Architecture:** Keep the implementation as a thin wrapper around Click's built-in shell completion classes. Add a dedicated `completion` command group under the top-level CLI, centralize script generation in one helper module, cover it with CLI regression tests, and document the static-script workflow in `README.md`.

**Tech Stack:** Python 3.10+, Click 8.1 shell completion APIs, pytest, `click.testing.CliRunner`

---

## File Map

| File | Responsibility |
|------|---------------|
| `phyloai/cli/main.py` | Register the new top-level `completion` command group alongside existing CLI groups |
| `phyloai/cli/completion.py` | Implement shell script generation and expose `bash`, `zsh`, and `fish` subcommands |
| `tests/cli/test_completion.py` | Cover help text, supported shells, script generation, and error behavior |
| `README.md` | Document one-time static completion script generation and shell setup |

---

### Task 1: Add failing CLI completion tests

**Files:**
- Create: `tests/cli/test_completion.py`
- Reference: `tests/cli/test_doctor.py`

- [ ] **Step 1: Write the failing test file**

```python
from click.testing import CliRunner

from phyloai.cli.main import cli


def test_completion_group_is_registered() -> None:
    runner = CliRunner()

    result = runner.invoke(cli, ["completion", "--help"])

    assert result.exit_code == 0
    assert "Generate shell completion scripts" in result.output
    assert "bash" in result.output
    assert "zsh" in result.output
    assert "fish" in result.output


def test_completion_bash_outputs_script() -> None:
    runner = CliRunner()

    result = runner.invoke(cli, ["completion", "bash"])

    assert result.exit_code == 0
    assert "complete -F" in result.output
    assert "phyloai" in result.output


def test_completion_zsh_outputs_script() -> None:
    runner = CliRunner()

    result = runner.invoke(cli, ["completion", "zsh"])

    assert result.exit_code == 0
    assert "#compdef phyloai" in result.output
    assert "_phyloai_completion" in result.output


def test_completion_fish_outputs_script() -> None:
    runner = CliRunner()

    result = runner.invoke(cli, ["completion", "fish"])

    assert result.exit_code == 0
    assert "complete --command phyloai" in result.output


def test_completion_help_explains_static_usage() -> None:
    runner = CliRunner()

    result = runner.invoke(cli, ["completion", "bash", "--help"])

    assert result.exit_code == 0
    assert "Print a Bash completion script" in result.output
    assert "static" in result.output.lower()
```

- [ ] **Step 2: Run the new tests and verify they fail**

Run: `pytest tests/cli/test_completion.py -v`

Expected: FAIL with a Click error indicating `No such command 'completion'`.

- [ ] **Step 3: Commit the failing tests**

```bash
git add tests/cli/test_completion.py
git commit -m "test(cli): add shell completion command coverage"
```

---

### Task 2: Implement the completion command

**Files:**
- Create: `phyloai/cli/completion.py`
- Modify: `phyloai/cli/main.py`
- Test: `tests/cli/test_completion.py`

- [ ] **Step 1: Create the completion module with a focused generator helper**

```python
"""Shell completion helpers for the PhyloAI CLI."""

from __future__ import annotations

import click
from click.shell_completion import get_completion_class


def _completion_script(shell: str) -> str:
    from phyloai.cli.main import cli

    complete_var = "_PHYLOAI_COMPLETE"
    completion_class = get_completion_class(shell)
    if completion_class is None:
        raise click.ClickException(f"Unsupported shell '{shell}'.")
    return completion_class(
        cli=cli,
        ctx_args={},
        prog_name="phyloai",
        complete_var=complete_var,
    ).source()


@click.group(help="Generate shell completion scripts for static installation.")
def completion() -> None:
    """Generate shell completion scripts for PhyloAI."""


@completion.command(help="Print a Bash completion script for static sourcing.")
def bash() -> None:
    click.echo(_completion_script("bash"), nl=False)


@completion.command(help="Print a Zsh completion script for static sourcing.")
def zsh() -> None:
    click.echo(_completion_script("zsh"), nl=False)


@completion.command(help="Print a Fish completion script for static sourcing.")
def fish() -> None:
    click.echo(_completion_script("fish"), nl=False)
```

- [ ] **Step 2: Register the new top-level command in `phyloai/cli/main.py`**

```python
from phyloai.cli.commands.pretree import pretree
from phyloai.cli.completion import completion
from phyloai.cli.doctor import doctor


cli.add_command(completion)
cli.add_command(doctor)
cli.add_command(pretree)
```

- [ ] **Step 3: Run the new completion tests and verify they pass**

Run: `pytest tests/cli/test_completion.py -v`

Expected: PASS for all completion command tests.

- [ ] **Step 4: Run the broader CLI test suite**

Run: `pytest tests/cli/ -q`

Expected: PASS with existing `doctor` tests still green.

- [ ] **Step 5: Commit the implementation**

```bash
git add phyloai/cli/main.py phyloai/cli/completion.py tests/cli/test_completion.py
git commit -m "feat(cli): add shell completion command"
```

---

### Task 3: Document the static completion workflow

**Files:**
- Modify: `README.md`
- Reference: `phyloai/cli/completion.py`
- Reference: `phyloai/cli/main.py`

- [ ] **Step 1: Add a README section for shell completion**

Insert a new section after the top-level CLI overview with content equivalent to:

```md
## Shell Completion

PhyloAI can generate static shell completion scripts for Bash, Zsh, and Fish:

```bash
phyloai completion bash
phyloai completion zsh
phyloai completion fish
```

The recommended setup is to generate the script once from an environment where `phyloai` is installed, save it to a persistent file, and source that file from your shell configuration.

Example for Zsh:

```bash
mkdir -p ~/.config/phyloai/completion
phyloai completion zsh > ~/.config/phyloai/completion/phyloai.zsh
```

Then add this to `~/.zshrc`:

```bash
source ~/.config/phyloai/completion/phyloai.zsh
```

Example for Bash:

```bash
mkdir -p ~/.config/phyloai/completion
phyloai completion bash > ~/.config/phyloai/completion/phyloai.bash
```

Then add this to `~/.bashrc`:

```bash
source ~/.config/phyloai/completion/phyloai.bash
```

Example for Fish:

```bash
mkdir -p ~/.config/fish/completions
phyloai completion fish > ~/.config/fish/completions/phyloai.fish
```

Avoid configuring shell startup to run `phyloai completion ...` dynamically on every new shell session. Static scripts are more robust for Conda-based workflows where the default shell environment may not contain `phyloai`.
```

- [ ] **Step 2: Sanity-check README command examples against the implemented CLI**

Run: `python -m phyloai.cli.main completion --help`

Expected: help output lists `bash`, `zsh`, and `fish` subcommands.

- [ ] **Step 3: Run focused CLI tests again after the docs update**

Run: `pytest tests/cli/test_completion.py tests/cli/test_doctor.py -q`

Expected: PASS.

- [ ] **Step 4: Commit the docs update**

```bash
git add README.md
git commit -m "docs: add CLI shell completion setup"
```

---

### Task 4: Final verification and integration

**Files:**
- Verify: `phyloai/cli/main.py`
- Verify: `phyloai/cli/completion.py`
- Verify: `tests/cli/test_completion.py`
- Verify: `README.md`

- [ ] **Step 1: Run end-to-end command verification for all supported shells**

Run: `python -m phyloai.cli.main completion bash >/tmp/phyloai.bash && python -m phyloai.cli.main completion zsh >/tmp/phyloai.zsh && python -m phyloai.cli.main completion fish >/tmp/phyloai.fish`

Expected: command exits 0 and writes three non-empty files.

- [ ] **Step 2: Run the full relevant test suite**

Run: `pytest tests/cli/ tests/pretree/test_stats.py -q`

Expected: PASS.

- [ ] **Step 3: Review working tree before final commit or push**

Run: `git status --short`

Expected: only the intended CLI completion files are modified or staged.

- [ ] **Step 4: Create the final integration commit if prior commits were squashed locally**

```bash
git add phyloai/cli/main.py phyloai/cli/completion.py tests/cli/test_completion.py README.md
git commit -m "feat(cli): add static shell completion support"
```

Use this step only if the work was not already committed incrementally in Tasks 1-3.

---

## Self-Review

- **Spec coverage:** The plan covers the approved command surface (`phyloai completion bash|zsh|fish`), zero new dependencies, static-script workflow, and README guidance that avoids dynamic shell-startup evaluation.
- **Placeholder scan:** No `TODO`/`TBD` placeholders remain. Each task names exact files, concrete test commands, and expected outcomes.
- **Type consistency:** The plan consistently uses a top-level `completion` Click group, a helper `_completion_script(shell: str) -> str`, and `CliRunner`-based regression tests.
