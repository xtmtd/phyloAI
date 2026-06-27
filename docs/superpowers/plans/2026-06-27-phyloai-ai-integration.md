# PhyloAI AI Integration Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use subagent-driven-development (recommended) or executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Implement the MCP Server (Phase 7) and Skill (Phase 8) for phyloAI AI integration.

**Architecture:** The MCP server dynamically generates tool schemas from the phyloai Click command tree at startup (zero hand-written schemas). All tools are fire-and-forget with `output_dir` as the persistent job handle. The Skill is a markdown document living in `skills/phyloai-workflow/` with parameter annotations, error catalog, dialog templates, and workflow reference.

**Tech Stack:** Python 3.10+, Click (existing CLI), `mcp` package for MCP protocol, stdio transport. No new language runtimes.

## Global Constraints

- All non-`doctor` commands write `result.json` — already enforced
- MCP server uses stdio transport only; HTTP deferred
- Skill lives inside phyloai repo at `skills/phyloai-workflow/`, version-coupled to CLI
- One MCP tool per CLI subcommand (fine-grained); stub tools for unimplemented commands
- All tools fire-and-forget; `output_dir` (absolute path) is the job handle
- MCP server entry point: `phyloai mcp-server`
- Relevant specs: `2026-06-07-phyloai-design.md` (main), `2026-06-27-phyloai-ai-integration-design.md` (AI)

---

## File Structure

```
phyloai/
├── mcp/                              # NEW: MCP server package
│   ├── __init__.py
│   ├── server.py                     # MCP stdio server, tool registration loop
│   ├── schema_gen.py                 # Click command tree → MCP JSON schemas
│   ├── tools/
│   │   ├── __init__.py
│   │   ├── cli_tools.py              # Fire-and-forget CLI wrappers
│   │   ├── utils.py                  # check_status, read_result, read_report, get_command_schema
│   │   └── stubs.py                  # Manual stub definitions for unimplemented commands
│   └── job.py                        # job.json write/read; launch helpers
├── cli/
│   ├── main.py                       # MODIFY: add `phyloai mcp-server` subcommand
│   └── commands/
│       └── mcp_server.py             # NEW: Click command that starts MCP server

skills/
└── phyloai-workflow/                 # NEW: Skill directory
    ├── SKILL.md                      # Main skill document
    └── references/
        ├── parameter-annotations.md  # Per-command Chinese annotations + recommended values
        ├── error-catalog.md          # Exit 1/3 error patterns + fixes
        ├── dialog-templates.md       # Pre-run/post-run card templates
        ├── demo-data.md              # Demo dataset paths and per-step entry points
        └── workflow.md               # Per-phase execution reference

phyloai/
└── demo_data/                        # NEW: bundled demo dataset
    ├── __init__.py                   # resolve_demo_path() helper
    ├── end_to_end/
    │   └── raw/                      # ~10 genes, ~20 taxa raw sequences
    │       ├── gene_001.faa
    │       ├── gene_002.faa
    │       └── ...
    └── per_step/
        ├── aligned/                  # pre-aligned sequences
        ├── trimmed/                  # pre-trimmed sequences
        └── gen_trees/               # pre-inferred gene trees

tests/
└── mcp/                              # NEW: MCP server tests
    ├── __init__.py
    ├── conftest.py                   # Fixtures: temp cli, sample data
    ├── test_schema_gen.py            # Schema generation tests
    ├── test_cli_tools.py             # CLI tool launch + job.json tests
    ├── test_utils.py                 # check_status, read_result, read_report tests
    ├── test_stubs.py                 # Stub tool tests
    └── test_job.py                   # job.json lifecycle tests
```

---

### Task 1: Infrastructure — Add dependency and package scaffolding

**Files:**
- Modify: `pyproject.toml`
- Create: `phyloai/mcp/__init__.py`
- Create: `phyloai/cli/commands/mcp_server.py`
- Modify: `phyloai/cli/main.py`

**Interfaces:**
- Consumes: (none — first task)
- Produces: `phyloai mcp-server` Click command (inert stub, starts server in later tasks)

- [ ] **Step 1: Add `mcp` dependency to pyproject.toml**

```toml
# pyproject.toml — add to dependencies:
"mcp>=1.0",
```

- [ ] **Step 2: Run pip install to verify dependency resolves**

Run: `pip install -e "."`
Expected: Success, `mcp` appears in `pip list | grep mcp`.

- [ ] **Step 3: Create `phyloai/mcp/__init__.py`**

```python
"""PhyloAI MCP Server — exposes CLI commands as MCP tools via stdio transport."""
```

- [ ] **Step 4: Create `phyloai/cli/commands/mcp_server.py` with inert Click command**

```python
"""phyloai mcp-server — start the MCP protocol server."""
from __future__ import annotations

import click


@click.command(name="mcp-server", hidden=False)
def mcp_server() -> None:
    """Start the PhyloAI MCP server (stdio transport)."""
    click.echo("MCP server not yet implemented.", err=True)
    raise click.exceptions.Exit(1)
```

- [ ] **Step 5: Register `mcp_server` command in `phyloai/cli/main.py`**

```python
# Add import at top, after existing imports:
from phyloai.cli.commands.mcp_server import mcp_server

# Add after cli.add_command(report):
cli.add_command(mcp_server)
```

- [ ] **Step 6: Verify the command registers correctly**

Run: `python -m phyloai.cli.main mcp-server`
Expected: `MCP server not yet implemented.` and exit code 1.

Run: `python -m phyloai.cli.main --help`
Expected: `mcp-server` appears in command list after `report`.

- [ ] **Step 7: Commit**

```bash
git add pyproject.toml phyloai/mcp/__init__.py phyloai/cli/commands/mcp_server.py phyloai/cli/main.py
git commit -m "feat(mcp): add phyloai mcp-server command stub and mcp dependency"
```

---

### Task 2: Schema Generation — Walk Click tree, produce MCP JSON schemas

**Files:**
- Create: `phyloai/mcp/schema_gen.py`
- Create: `tests/mcp/__init__.py`
- Create: `tests/mcp/conftest.py`
- Create: `tests/mcp/test_schema_gen.py`

**Interfaces:**
- Consumes: `phyloai.cli.main.cli` (Click root group)
- Produces:
  - `walk_click_tree(root: click.Group) -> list[dict]` — returns list of command descriptors, each with `{"tool_name": str, "command_path": list[str], "click_command": click.Command, "help": str}`
  - `click_param_to_json_schema(param: click.Parameter) -> dict` — maps one Click param to JSON Schema property
  - `build_mcp_tool(descriptor: dict) -> dict` — full MCP tool definition: `{name, description, inputSchema}`

- [ ] **Step 1: Write failing test for `click_param_to_json_schema` — required string option**

```python
# tests/mcp/test_schema_gen.py
import click
from phyloai.mcp.schema_gen import click_param_to_json_schema


def test_string_option_default():
    opt = click.Option(["--name"], type=str, default="default_val", help="A name")
    result = click_param_to_json_schema(opt)
    assert result["type"] == "string"
    assert result["default"] == "default_val"
    assert result["description"] == "A name"
    assert "enum" not in result
```

- [ ] **Step 2: Run test to verify it fails**

Run: `python -m pytest tests/mcp/test_schema_gen.py::test_string_option_default -v`
Expected: FAIL with ImportError (module doesn't exist).

- [ ] **Step 3: Write failing test for `click_param_to_json_schema` — Choice param**

```python
def test_choice_option():
    opt = click.Option(["--method"], type=click.Choice(["linsi", "auto"]), default="linsi")
    result = click_param_to_json_schema(opt)
    assert result["type"] == "string"
    assert result["default"] == "linsi"
    assert result["enum"] == ["linsi", "auto"]
```

- [ ] **Step 4: Write failing test for `click_param_to_json_schema` — Path param**

```python
import click
from pathlib import Path

def test_path_param():
    opt = click.Option(["--matrix"], type=click.Path(path_type=Path), required=True,
                       help="Supermatrix alignment")
    result = click_param_to_json_schema(opt)
    assert result["type"] == "string"
    assert result["description"] == "Supermatrix alignment"
    # click_path hints when type is click.Path
    assert result.get("format") == "path"
```

- [ ] **Step 5: Write failing test for `click_param_to_json_schema` — int with range**

```python
def test_int_option():
    opt = click.Option(["--threads"], type=int, default=4)
    result = click_param_to_json_schema(opt)
    assert result["type"] == "integer"
    assert result["default"] == 4
```

- [ ] **Step 6: Write failing test for `click_param_to_json_schema` — flag**

```python
def test_flag_option():
    opt = click.Option(["--overwrite"], is_flag=True, default=False)
    result = click_param_to_json_schema(opt)
    assert result["type"] == "boolean"
    assert result["default"] is False
```

- [ ] **Step 7: Write failing test for `walk_click_tree` — finds known commands**

```python
from phyloai.cli.main import cli
from phyloai.mcp.schema_gen import walk_click_tree


def test_walk_click_tree_finds_known_commands():
    descriptors = walk_click_tree(cli)
    tool_names = {d["tool_name"] for d in descriptors}
    # Leaf commands that must exist
    assert "pretree_convert" in tool_names
    assert "pretree_align" in tool_names
    assert "tree_ml_iqtree" in tool_names
    assert "doctor" in tool_names
    assert "report" in tool_names
    assert "posttree_dating_hessian" in tool_names
    # Groups (non-leaf) should not appear
    assert "pretree" not in tool_names
    assert "tree" not in tool_names
    assert "posttree_dating" not in tool_names


def test_walk_click_tree_each_descriptor_has_required_keys():
    descriptors = walk_click_tree(cli)
    assert len(descriptors) > 0
    for d in descriptors:
        assert "tool_name" in d
        assert "command_path" in d
        assert "click_command" in d
        assert "help" in d
        assert isinstance(d["command_path"], list)
```

- [ ] **Step 8: Write failing test for `build_mcp_tool` — produces valid MCP tool definition**

```python
from phyloai.mcp.schema_gen import build_mcp_tool, walk_click_tree
from phyloai.cli.main import cli


def test_build_mcp_tool_doctor():
    descriptors = walk_click_tree(cli)
    doctor_descriptor = next(d for d in descriptors if d["tool_name"] == "doctor")
    tool_def = build_mcp_tool(doctor_descriptor)
    assert tool_def["name"] == "doctor"
    assert "description" in tool_def
    assert "inputSchema" in tool_def
    schema = tool_def["inputSchema"]
    assert schema["type"] == "object"
    assert "properties" in schema
    # doctor has --output-format
    assert "output_format" in schema["properties"]
```

- [ ] **Step 9: Run all tests to verify they fail**

Run: `python -m pytest tests/mcp/test_schema_gen.py -v`
Expected: All 7 tests FAIL with ImportError.

- [ ] **Step 10: Implement `phyloai/mcp/schema_gen.py`**

```python
"""Generate MCP JSON schemas from the phyloai Click command tree."""
from __future__ import annotations

from pathlib import Path
from typing import Any

import click


def walk_click_tree(root: click.Group) -> list[dict[str, Any]]:
    """Walk the Click command tree and return a flat list of leaf command descriptors.

    Tool names are built by joining group names with underscores, skipping the
    root 'cli' group.  Groups with invoke_without_command=True are treated as
    leaf commands if they have a callback.
    """
    descriptors: list[dict[str, Any]] = []

    def _walk(group: click.Group, path_parts: list[str]) -> None:
        for cmd_name in group.list_commands(None):
            cmd = group.get_command(None, cmd_name)
            if cmd is None:
                continue
            new_parts = path_parts + [cmd_name]
            if isinstance(cmd, click.Group):
                # Sub-group: recurse, but also treat as leaf if it has a
                # callback and is invoke_without_command
                if cmd.invoke_without_command and cmd.callback is not None:
                    descriptors.append({
                        "tool_name": "_".join(new_parts),
                        "command_path": new_parts,
                        "click_command": cmd,
                        "help": cmd.help or "",
                    })
                _walk(cmd, new_parts)
            else:
                # Leaf command
                descriptors.append({
                    "tool_name": "_".join(new_parts),
                    "command_path": new_parts,
                    "click_command": cmd,
                    "help": cmd.help or "",
                })

    _walk(root, [])
    return descriptors


def click_param_to_json_schema(param: click.Parameter) -> dict[str, Any]:
    """Convert a Click parameter to a JSON Schema property descriptor."""
    schema: dict[str, Any] = {
        "description": param.help or "",
    }

    # Determine JSON type
    if param.is_flag:
        schema["type"] = "boolean"
    elif isinstance(param.type, click.IntRange):
        schema["type"] = "integer"
    elif isinstance(param.type, click.types.IntParamType):
        schema["type"] = "integer"
    elif isinstance(param.type, click.types.FloatParamType):
        schema["type"] = "number"
    elif isinstance(param.type, click.Choice):
        schema["type"] = "string"
        schema["enum"] = list(param.type.choices)
    elif isinstance(param.type, (click.Path, click.File)):
        schema["type"] = "string"
        schema["format"] = "path"
    else:
        schema["type"] = "string"

    # Default value
    if param.default is not None and param.default is not ...:
        # Convert Path defaults to string
        default_val = param.default
        if isinstance(default_val, Path):
            default_val = str(default_val)
        schema["default"] = default_val

    return schema


def build_mcp_tool(descriptor: dict[str, Any]) -> dict[str, Any]:
    """Build a full MCP tool definition from a command descriptor."""
    cmd: click.Command = descriptor["click_command"]
    properties: dict[str, dict[str, Any]] = {}
    required: list[str] = []

    for param in cmd.params:
        if param.hidden:
            continue
        prop = click_param_to_json_schema(param)
        properties[param.name] = prop
        if param.required:
            required.append(param.name)

    return {
        "name": descriptor["tool_name"],
        "description": descriptor["help"] or cmd.help or "",
        "inputSchema": {
            "type": "object",
            "properties": properties,
            "required": required,
        },
    }
```

- [ ] **Step 11: Run tests to verify they pass**

Run: `python -m pytest tests/mcp/test_schema_gen.py -v`
Expected: All 7 tests PASS.

- [ ] **Step 12: Commit**

```bash
git add phyloai/mcp/schema_gen.py tests/mcp/
git commit -m "feat(mcp): add Click-to-MCP schema generation"
```

---

### Task 3: Job Management — job.json write/read and launch helpers

**Files:**
- Create: `phyloai/mcp/job.py`
- Create: `tests/mcp/test_job.py`

**Interfaces:**
- Consumes: None (pure file I/O)
- Produces:
  - `write_job_json(output_dir: Path, pid: int, command: str, *, early_exit_stderr: str = "") -> dict` — writes job.json, returns payload dict
  - `read_job_json(output_dir: Path) -> dict | None` — reads job.json, returns None if missing
  - `build_cli_argv(descriptor: dict, params: dict[str, Any]) -> list[str]` — builds `["phyloai", "pretree", "align", "--seq-dir", ...]` list from command_path + resolved params. Requires descriptor to include `click_command` with `params` to resolve flag names.
  - `launch_cli(descriptor: dict, params: dict[str, Any], output_dir: Path) -> tuple[Path, int]` — validates path, starts subprocess.Popen, writes job.json, returns (output_dir, pid). Raises ValueError on invalid path or subprocess start failure.

- [ ] **Step 1: Write failing test for `write_job_json` and `read_job_json`**

```python
# tests/mcp/test_job.py
import json
import tempfile
from pathlib import Path
from phyloai.mcp.job import write_job_json, read_job_json


def test_write_and_read_job_json():
    with tempfile.TemporaryDirectory() as tmp:
        output_dir = Path(tmp)
        payload = write_job_json(output_dir, pid=12345, command="phyloai pretree align --seq-dir ./raw")
        assert (output_dir / "job.json").exists()
        assert payload["pid"] == 12345
        assert payload["command"] == "phyloai pretree align --seq-dir ./raw"
        assert "started_at" in payload

        loaded = read_job_json(output_dir)
        assert loaded is not None
        assert loaded["pid"] == 12345
        assert loaded["command"] == "phyloai pretree align --seq-dir ./raw"


def test_read_job_json_missing():
    with tempfile.TemporaryDirectory() as tmp:
        loaded = read_job_json(Path(tmp))
        assert loaded is None
```

- [ ] **Step 2: Write failing test for `build_cli_argv`**

```python
from phyloai.mcp.job import build_cli_argv


def test_build_cli_argv():
    import click
    cmd = click.Command(
        name="align",
        params=[
            click.Option(["--seq-dir"], type=str, required=True),
            click.Option(["--method"], type=click.Choice(["linsi", "auto"]), default="linsi"),
            click.Option(["-t", "--threads"], type=int, default=4),
            click.Option(["--tool-args"], type=str, default=None),
        ],
    )
    descriptor = {"command_path": ["pretree", "align"], "click_command": cmd}
    params = {"seq_dir": "./raw", "method": "linsi", "threads": 8, "tool_args": None}
    result = build_cli_argv(descriptor, params)
    assert result[0] == "phyloai"
    assert "--seq-dir" in result
    assert "./raw" in result
    assert "--method" in result
    assert "linsi" in result
    assert "--threads" in result
    assert "8" in result
    # tool_args=None should be omitted
    assert "--tool-args" not in result
    assert None not in result


def test_build_cli_argv_omits_none_and_false_flags():
    import click
    cmd = click.Command(
        name="report",
        params=[
            click.Option(["--run-dir"], type=str, required=True),
            click.Option(["-q", "--quiet"], is_flag=True, default=False),
            click.Option(["--overwrite"], is_flag=True, default=False),
        ],
    )
    descriptor = {"command_path": ["report"], "click_command": cmd}
    params = {"run_dir": "./runs/run", "quiet": False, "overwrite": True}
    result = build_cli_argv(descriptor, params)
    assert "--quiet" not in result  # False flag omitted
    assert "-q" not in result
    assert "--overwrite" in result   # True flag included
    assert "-q" not in result
```

- [ ] **Step 3: Run tests to verify they fail**

Run: `python -m pytest tests/mcp/test_job.py -v`
Expected: All 4 tests FAIL with ImportError.

- [ ] **Step 4: Implement `phyloai/mcp/job.py`**

```python
"""Job lifecycle: job.json, CLI argv building, subprocess launching."""
from __future__ import annotations

import json
import os
import subprocess
import sys
import time
from pathlib import Path
from typing import Any


def write_job_json(output_dir: Path, pid: int, command: str,
                   *, early_exit_stderr: str = "") -> dict[str, Any]:
    """Write job.json to output_dir and return the payload."""
    output_dir.mkdir(parents=True, exist_ok=True)
    payload = {
        "pid": pid,
        "started_at": time.strftime("%Y-%m-%dT%H:%M:%S"),
        "command": command,
    }
    if early_exit_stderr:
        payload["early_exit_stderr"] = early_exit_stderr
    with open(output_dir / "job.json", "w") as fh:
        json.dump(payload, fh, indent=2)
    return payload


def read_job_json(output_dir: Path) -> dict[str, Any] | None:
    """Read job.json from output_dir, return None if missing or unreadable."""
    job_path = output_dir / "job.json"
    if not job_path.exists():
        return None
    try:
        with open(job_path) as fh:
            return json.load(fh)
    except (json.JSONDecodeError, OSError):
        return None


def build_cli_argv(descriptor: dict[str, Any], params: dict[str, Any]) -> list[str]:
    """Build a CLI argv list from a command descriptor and resolved parameters.

    Converts param dict keys (underscore) to Click flag names (kebab-case).
    Boolean flags: include only when True.  None values: omitted.
    descriptor MUST include 'click_command' with Click params.
    """
    argv = ["phyloai"] + descriptor["command_path"]

    cmd = descriptor.get("click_command")
    if cmd is None:
        # Fallback: append all params as --key value (best-effort, no type info)
        for key, value in params.items():
            if value is None or value is False:
                continue
            flag = f"--{key.replace('_', '-')}"
            if value is True:
                argv.append(flag)
            else:
                argv.extend([flag, str(value)])
        return argv

    for param in cmd.params:
        if param.hidden:
            continue
        key = param.name
        value = params.get(key, param.default if param.default is not ... else None)

        if value is None and not param.is_flag:
            continue
        if param.is_flag:
            if value:
                flag = max(param.opts, key=len) if param.opts else f"--{key}"
                argv.append(flag)
            continue

        flag = max(param.opts, key=len) if param.opts else f"--{key.replace('_', '-')}"
        argv.append(flag)
        argv.append(str(value))

    return argv


def launch_cli(
    descriptor: dict[str, Any],
    params: dict[str, Any],
    output_dir: Path,
    *,
    env: dict[str, str] | None = None,
) -> tuple[Path, int]:
    """Build CLI argv, validate output_dir, launch subprocess, write job.json.

    Returns (output_dir, pid).  Raises ValueError on invalid output_dir or
    subprocess start failure.
    """
    output_dir = output_dir.resolve()

    if not output_dir.parent.exists():
        raise ValueError(f"Parent directory does not exist: {output_dir.parent}")

    argv = build_cli_argv(descriptor, params)
    command_str = " ".join(argv)

    proc_env = dict(os.environ)
    if env:
        proc_env.update(env)

    try:
        proc = subprocess.Popen(
            argv,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.PIPE,
            stdin=subprocess.DEVNULL,
            env=proc_env,
            start_new_session=True,
        )
    except FileNotFoundError as e:
        raise ValueError(f"Failed to start phyloai CLI: {e}")
    except OSError as e:
        raise ValueError(f"Failed to start subprocess: {e}")

    # Wait briefly to detect immediate launch failures (e.g. import errors)
    try:
        stdout_data, stderr_data = proc.communicate(timeout=2)
        raise ValueError(
            f"Process exited immediately with code {proc.returncode}. "
            f"stderr: {stderr_data.decode('utf-8', errors='replace')[:500]}"
        )
    except subprocess.TimeoutExpired:
        pass  # Process is still running — good

    write_job_json(output_dir, pid=proc.pid, command=command_str)
    return output_dir, proc.pid
```

- [ ] **Step 5: Run tests to verify they pass**

Run: `python -m pytest tests/mcp/test_job.py -v`
Expected: All 4 tests PASS.

- [ ] **Step 6: Commit**

```bash
git add phyloai/mcp/job.py tests/mcp/test_job.py
git commit -m "feat(mcp): add job.json lifecycle and CLI argv builder"
```

---

### Task 4: Utility Tools — check_status, read_result, read_report, get_command_schema

**Files:**
- Create: `phyloai/mcp/tools/__init__.py`
- Create: `phyloai/mcp/tools/utils.py`
- Create: `tests/mcp/test_utils.py`

**Interfaces:**
- Consumes: `phyloai.mcp.job` (read_job_json), `phyloai.mcp.schema_gen` (walk_click_tree, build_mcp_tool)
- Produces:
  - `check_status(output_dir: str) -> dict` — returns `{status: not_started|running|success|error|unknown, output_dir, ...}`
  - `read_result(output_dir: str) -> dict` — reads result.json, returns parsed dict or error
  - `read_report(run_dir: str) -> dict` — reads report.json, returns parsed dict or error
  - `get_command_schema(command_name: str) -> dict` — looks up a command by tool_name, returns its MCP tool definition

- [ ] **Step 1: Write failing test for `check_status` — all five states**

```python
# tests/mcp/test_utils.py
import json
import time
import tempfile
from pathlib import Path
from phyloai.mcp.tools.utils import check_status


def test_check_status_not_started():
    with tempfile.TemporaryDirectory() as tmp:
        output_dir = Path(tmp)
        result = check_status(str(output_dir))
        assert result["status"] == "not_started"
        assert result["output_dir"] == str(output_dir.resolve())


def test_check_status_running():
    with tempfile.TemporaryDirectory() as tmp:
        output_dir = Path(tmp)
        # Write job.json (no result.json, no checkpoint.json)
        job = {"pid": 12345, "started_at": time.strftime("%Y-%m-%dT%H:%M:%S"),
               "command": "phyloai pretree align"}
        with open(output_dir / "job.json", "w") as fh:
            json.dump(job, fh)
        result = check_status(str(output_dir))
        assert result["status"] == "running"


def test_check_status_running_with_checkpoint():
    with tempfile.TemporaryDirectory() as tmp:
        output_dir = Path(tmp)
        job = {"pid": 12345, "started_at": time.strftime("%Y-%m-%dT%H:%M:%S"),
               "command": "phyloai pretree align"}
        checkpoint = {"completed": 42, "total": 200, "last_updated": "2026-06-27T14:23:00"}
        with open(output_dir / "job.json", "w") as fh:
            json.dump(job, fh)
        with open(output_dir / "checkpoint.json", "w") as fh:
            json.dump(checkpoint, fh)
        result = check_status(str(output_dir))
        assert result["status"] == "running"
        assert result["checkpoint"] == checkpoint


def test_check_status_success():
    with tempfile.TemporaryDirectory() as tmp:
        output_dir = Path(tmp)
        job = {"pid": 12345, "started_at": time.strftime("%Y-%m-%dT%H:%M:%S"),
               "command": "phyloai pretree align"}
        with open(output_dir / "job.json", "w") as fh:
            json.dump(job, fh)
        with open(output_dir / "result.json", "w") as fh:
            json.dump({"status": "success", "key_results": {"n_success": 5}}, fh)
        result = check_status(str(output_dir))
        assert result["status"] == "success"
        assert result["result"]["key_results"]["n_success"] == 5


def test_check_status_error():
    with tempfile.TemporaryDirectory() as tmp:
        output_dir = Path(tmp)
        with open(output_dir / "result.json", "w") as fh:
            json.dump({"status": "error", "error": "Bad input"}, fh)
        result = check_status(str(output_dir))
        assert result["status"] == "error"
        assert result["result"]["error"] == "Bad input"


def test_check_status_unknown():
    with tempfile.TemporaryDirectory() as tmp:
        output_dir = Path(tmp)
        # job.json exists but no result.json, and no process
        # simulate this by writing a fake pid job
        job = {"pid": 99999, "started_at": "2020-01-01T00:00:00",
               "command": "phyloai pretree align"}
        with open(output_dir / "job.json", "w") as fh:
            json.dump(job, fh)
        result = check_status(str(output_dir))
        # pid 99999 unlikely exists; should be unknown
        assert result["status"] == "unknown"
```

- [ ] **Step 2: Write failing test for `read_result`**

```python
from phyloai.mcp.tools.utils import read_result


def test_read_result_success():
    with tempfile.TemporaryDirectory() as tmp:
        output_dir = Path(tmp)
        payload = {"status": "success", "command": "phyloai pretree stats",
                   "wall_time": 1.5}
        with open(output_dir / "result.json", "w") as fh:
            json.dump(payload, fh)
        result = read_result(str(output_dir))
        assert result["status"] == "success"
        assert result["wall_time"] == 1.5


def test_read_result_missing():
    with tempfile.TemporaryDirectory() as tmp:
        result = read_result(str(Path(tmp)))
        assert result["status"] == "error"
        assert "not found" in result["message"].lower()
```

- [ ] **Step 3: Write failing test for `read_report`**

```python
from phyloai.mcp.tools.utils import read_report


def test_read_report_success():
    with tempfile.TemporaryDirectory() as tmp:
        run_dir = Path(tmp)
        report = {"report_version": "1.0", "steps": []}
        with open(run_dir / "report" / "report.json", "w") as fh:
            json.dump(report, fh)
        result = read_report(str(run_dir))
        assert result["report_version"] == "1.0"


def test_read_report_missing():
    with tempfile.TemporaryDirectory() as tmp:
        result = read_report(str(Path(tmp)))
        assert result["status"] == "error"
        assert "not found" in result["message"].lower()
```

- [ ] **Step 4: Write failing test for `get_command_schema`**

```python
from phyloai.mcp.tools.utils import get_command_schema


def test_get_command_schema_doctor():
    schema = get_command_schema("doctor")
    assert schema["name"] == "doctor"
    assert "inputSchema" in schema
    assert "output_format" in schema["inputSchema"]["properties"]


def test_get_command_schema_unknown():
    schema = get_command_schema("nonexistent_command")
    assert "error" in schema
```

- [ ] **Step 5: Run tests to verify they fail**

Run: `python -m pytest tests/mcp/test_utils.py -v`
Expected: All tests FAIL with ImportError.

- [ ] **Step 6: Implement `phyloai/mcp/tools/utils.py`**

```python
"""MCP utility tools: check_status, read_result, read_report, get_command_schema."""
from __future__ import annotations

import json
import os
from pathlib import Path

from phyloai.mcp.job import read_job_json
from phyloai.mcp.schema_gen import build_mcp_tool, walk_click_tree


# Lazily populated on first call to get_command_schema
_schema_cache: dict[str, dict] | None = None


def _load_schema_cache() -> dict[str, dict]:
    global _schema_cache
    if _schema_cache is not None:
        return _schema_cache
    from phyloai.cli.main import cli
    descriptors = walk_click_tree(cli)
    _schema_cache = {d["tool_name"]: build_mcp_tool(d) for d in descriptors}
    return _schema_cache


def check_status(output_dir: str) -> dict:
    """Return the current state of a job by inspecting its output directory."""
    od = Path(output_dir).resolve()
    result_path = od / "result.json"
    checkpoint_path = od / "checkpoint.json"
    job = read_job_json(od)

    # Check for completed/error result first
    if result_path.exists():
        try:
            with open(result_path) as fh:
                result_data = json.load(fh)
            status = result_data.get("status", "unknown")
            return {
                "status": status,
                "output_dir": str(od),
                "result": result_data,
            }
        except (json.JSONDecodeError, OSError):
            return {"status": "unknown", "output_dir": str(od),
                    "message": "result.json exists but is unreadable"}

    # Read checkpoint.json for progress info
    checkpoint = None
    if checkpoint_path.exists():
        try:
            with open(checkpoint_path) as fh:
                checkpoint = json.load(fh)
        except (json.JSONDecodeError, OSError):
            pass

    # No result.json yet
    if job is None:
        return {"status": "not_started", "output_dir": str(od)}

    # job.json present, no result.json — check if process still alive
    pid = job.get("pid")
    if pid is not None:
        try:
            os.kill(int(pid), 0)
            response = {"status": "running", "output_dir": str(od)}
            if checkpoint:
                response["checkpoint"] = checkpoint
            return response
        except (OSError, ValueError):
            pass

    response = {"status": "unknown", "output_dir": str(od),
                "message": "Process exited but result.json not found. Check logs/ for tool stderr."}
    if checkpoint:
        response["checkpoint"] = checkpoint
    return response


def read_result(output_dir: str) -> dict:
    """Read result.json from a step output directory."""
    od = Path(output_dir).resolve()
    result_path = od / "result.json"
    if not result_path.exists():
        return {"status": "error", "message": f"result.json not found at {od}"}
    try:
        with open(result_path) as fh:
            return json.load(fh)
    except (json.JSONDecodeError, OSError) as e:
        return {"status": "error", "message": f"Failed to read result.json: {e}"}


def read_report(run_dir: str) -> dict:
    """Read report.json from a run directory (looks under <run_dir>/report/)."""
    rd = Path(run_dir).resolve()
    report_path = rd / "report" / "report.json"
    if not report_path.exists():
        return {"status": "error",
                "message": f"report.json not found at {report_path}. Run 'phyloai report --run-dir {rd}' to generate it."}
    try:
        with open(report_path) as fh:
            return json.load(fh)
    except (json.JSONDecodeError, OSError) as e:
        return {"status": "error", "message": f"Failed to read report.json: {e}"}


def get_command_schema(command_name: str) -> dict:
    """Return the MCP tool schema for a given command (by tool_name).

    Returns an error dict if the command is unknown.  Uses a lazy cache so
    the Click tree is walked only once.
    """
    cache = _load_schema_cache()
    if command_name in cache:
        return cache[command_name]
    return {"error": f"Unknown command: {command_name}",
            "available": sorted(cache.keys())}
```

- [ ] **Step 7: Run tests to verify they pass**

Run: `python -m pytest tests/mcp/test_utils.py -v`
Expected: All tests PASS.

- [ ] **Step 8: Commit**

```bash
git add phyloai/mcp/tools/__init__.py phyloai/mcp/tools/utils.py tests/mcp/test_utils.py
git commit -m "feat(mcp): add utility tools — check_status, read_result, read_report, get_command_schema"
```

---

### Task 5: Stub Tools — Manual definitions for unimplemented commands

**Files:**
- Create: `phyloai/mcp/tools/stubs.py`
- Create: `tests/mcp/test_stubs.py`

**Interfaces:**
- Consumes: None
- Produces:
  - `STUB_TOOL_NAMES: frozenset[str]` — set of tool names that are stubs
  - `STUB_TOOLS: list[dict]` — list of MCP tool definitions with stub implementations
  - `handle_stub(tool_name: str) -> dict` — returns the stub response if tool_name is a stub, else None

- [ ] **Step 1: Write failing test for stub tools**

```python
# tests/mcp/test_stubs.py
from phyloai.mcp.tools.stubs import STUB_TOOL_NAMES, STUB_TOOLS, handle_stub
from phyloai.mcp.schema_gen import walk_click_tree
from phyloai.cli.main import cli


def test_stub_tool_names():
    expected = {"posttree_signal", "posttree_simulate",
                "posttree_syserror_brlen", "posttree_syserror_cca",
                "posttree_syserror_sites"}
    assert STUB_TOOL_NAMES == expected


def test_stub_tools_have_valid_definitions():
    for tool_def in STUB_TOOLS:
        assert "name" in tool_def
        assert tool_def["name"] in STUB_TOOL_NAMES
        assert "description" in tool_def
        assert "inputSchema" in tool_def


def test_handle_stub_returns_response():
    result = handle_stub("posttree_syserror_brlen")
    assert result is not None
    assert result["status"] == "not_implemented"


def test_handle_stub_non_stub_returns_none():
    result = handle_stub("pretree_align")
    assert result is None


def test_stubs_do_not_overlap_with_dynamic_tools():
    """Stub names must not appear in the Click command tree — if they do,
    they are already real commands and should NOT be stubbed."""
    descriptors = walk_click_tree(cli)
    dynamic_names = {d["tool_name"] for d in descriptors}
    overlap = STUB_TOOL_NAMES & dynamic_names
    assert overlap == set(), (
        f"These are listed as stubs but actually exist in the CLI: {overlap}. "
        "Remove them from STUB_TOOL_NAMES."
    )
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `python -m pytest tests/mcp/test_stubs.py -v`
Expected: All 4 tests FAIL with ImportError.

- [ ] **Step 3: Implement `phyloai/mcp/tools/stubs.py`**

```python
"""Manual stub definitions for CLI commands not yet implemented.

Once a command is implemented in the CLI, the dynamic schema generation
(phyloai.mcp.schema_gen) picks it up automatically.  Stubs here serve as
placeholders that return a polite "not implemented" response, so the AI
Skill knows these commands exist in the design but are not yet available.
"""
from __future__ import annotations

STUB_TOOL_NAMES: frozenset[str] = frozenset({
    "posttree_signal",
    "posttree_simulate",
    "posttree_syserror_brlen",
    "posttree_syserror_cca",
    "posttree_syserror_sites",
})

_STUB_DESCRIPTIONS = {
    "posttree_signal": "Phylogenetic signal distribution analysis (future).",
    "posttree_simulate": "AliSim MSA simulation / gene-jackknife resampling (future).",
    "posttree_syserror_brlen": "Systematic error diagnosis: branch-length variation (future).",
    "posttree_syserror_cca": "Systematic error diagnosis: cross-composition analysis (future).",
    "posttree_syserror_sites": "Systematic error diagnosis: site-wise analysis (future).",
}

STUB_TOOLS: list[dict] = [
    {
        "name": name,
        "description": desc,
        "inputSchema": {"type": "object", "properties": {}, "required": []},
    }
    for name, desc in _STUB_DESCRIPTIONS.items()
]


def handle_stub(tool_name: str) -> dict | None:
    """Return stub response if tool_name is a stub, else None."""
    if tool_name in STUB_TOOL_NAMES:
        return {
            "status": "not_implemented",
            "message": f"The '{tool_name}' command is not yet available in the installed version.",
        }
    return None
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `python -m pytest tests/mcp/test_stubs.py -v`
Expected: All 4 tests PASS.

- [ ] **Step 5: Commit**

```bash
git add phyloai/mcp/tools/stubs.py tests/mcp/test_stubs.py
git commit -m "feat(mcp): add stub tools for unimplemented posttree commands"
```

---

### Task 6: CLI Tool Wrappers — Fire-and-forget execution per command

**Files:**
- Create: `phyloai/mcp/tools/cli_tools.py`
- Create: `tests/mcp/test_cli_tools.py`

**Interfaces:**
- Consumes: `phyloai.mcp.job` (launch_cli), `phyloai.mcp.tools.stubs` (handle_stub), `phyloai.mcp.schema_gen` (walk_click_tree, build_mcp_tool)
- Produces:
  - `register_all_tools(mcp_server) -> None` — registers all CLI + stub tools on the MCP server instance
  - Each tool handler: async function that calls launch_cli and returns `{output_dir, pid, message}`

- [ ] **Step 1: Write failing test for `register_all_tools` exhaustiveness**

```python
# tests/mcp/test_cli_tools.py
import asyncio

# We test tool registration logic without a real MCP server by inspecting
# the tool definitions themselves.
from phyloai.mcp.schema_gen import walk_click_tree
from phyloai.mcp.tools.stubs import STUB_TOOL_NAMES
from phyloai.cli.main import cli


def test_all_dynamic_commands_unless_stub():
    """Every leaf in the Click tree must become either a dynamic tool or a stub."""
    descriptors = walk_click_tree(cli)
    dynamic_names = {d["tool_name"] for d in descriptors}

    # Commands that ARE stubs should NOT appear in dynamic tree
    overlap = dynamic_names & STUB_TOOL_NAMES
    assert overlap == set(), f"Stub commands found in dynamic tree: {overlap}"
```

- [ ] **Step 2: Write failing test for launch_cli integration**

```python
import tempfile
from pathlib import Path
from phyloai.mcp.job import launch_cli, read_job_json


def test_launch_cli_writes_valid_job_json():
    import click
    import tempfile
    from pathlib import Path
    from phyloai.mcp.job import launch_cli, read_job_json

    cmd = click.Command(
        name="stats",
        params=[
            click.Option(["--seq-dir"], type=str, required=True),
            click.Option(["-o", "--output-dir"], type=click.Path(path_type=Path),
                         default=Path("runs/pretree/stats")),
        ],
    )
    descriptor = {"command_path": ["pretree", "stats"], "click_command": cmd}
    with tempfile.TemporaryDirectory() as tmp:
        output_dir = Path(tmp) / "stats"
        result_dir, pid = launch_cli(descriptor, {"seq_dir": "./data"}, output_dir)
        job = read_job_json(result_dir)
        assert job is not None
        assert job["pid"] == pid
        assert "stats" in job["command"]
        assert result_dir == output_dir.resolve()
```

- [ ] **Step 3: Run tests to verify they fail**

Run: `python -m pytest tests/mcp/test_cli_tools.py -v`
Expected: Both tests FAIL.

- [ ] **Step 4: Implement `phyloai/mcp/tools/cli_tools.py`**

```python
"""MCP CLI tool wrappers — fire-and-forget subprocess execution."""
from __future__ import annotations

import json
import os
from pathlib import Path
from typing import Any

from phyloai.mcp.job import launch_cli
from phyloai.mcp.schema_gen import build_mcp_tool, walk_click_tree
from phyloai.mcp.tools.stubs import STUB_TOOLS, handle_stub


def _resolve_output_dir(
    descriptor: dict[str, Any],
    kwargs: dict[str, Any],
) -> Path:
    """Resolve the output directory for fire-and-forget commands.

    doctor and report are handled synchronously, never call this.
    """
    tool_name = descriptor["tool_name"]
    output_dir_raw = kwargs.get("output_dir") or kwargs.get("output-dir")
    if output_dir_raw is None:
        raise ValueError(f"output_dir is required for {tool_name}")
    return Path(output_dir_raw).resolve()


def _make_launch_handler(
    descriptor: dict[str, Any],
) -> Any:
    """Create an async handler function for a single CLI tool."""
    tool_name = descriptor["tool_name"]

    async def handler(**kwargs: Any) -> str:
        # Check stub first
        stub_response = handle_stub(tool_name)
        if stub_response is not None:
            return json.dumps(stub_response)

        # doctor: synchronous subprocess.run, return raw JSON
        if tool_name == "doctor":
            result = subprocess.run(
                ["phyloai", "doctor", "--output-format", "json"],
                capture_output=True, text=True, timeout=30,
            )
            if result.returncode != 0:
                return json.dumps({
                    "status": "error",
                    "message": f"doctor failed (code {result.returncode}): {result.stderr[:500]}",
                    "stdout": result.stdout,
                })
            return result.stdout

        # report: synchronous subprocess.run (pure Python, fast, caller wants result immediately)
        if tool_name == "report":
            run_dir = kwargs.get("run_dir") or kwargs.get("run-dir")
            if run_dir is None:
                return json.dumps({"status": "error", "message": "run_dir is required"})
            result = subprocess.run(
                ["phyloai", "report", "--run-dir", str(run_dir)],
                capture_output=True, text=True, timeout=300,
            )
            report_path = Path(run_dir) / "report" / "report.json"
            if report_path.exists():
                try:
                    with open(report_path) as fh:
                        return json.dumps(json.load(fh))
                except (json.JSONDecodeError, OSError):
                    pass
            if result.returncode != 0:
                return json.dumps({
                    "status": "error",
                    "message": f"report failed (code {result.returncode}): {result.stderr[:500]}",
                })
            return json.dumps({
                "status": "error",
                "message": "report completed but report.json not found",
            })

        # All other commands: fire-and-forget
        try:
            output_dir = _resolve_output_dir(descriptor, kwargs)
        except ValueError as e:
            return json.dumps({"status": "error", "message": str(e)})

        try:
            result_dir, pid = launch_cli(descriptor, kwargs, output_dir)
            return json.dumps({
                "status": "launched",
                "output_dir": str(result_dir),
                "pid": pid,
                "message": (
                    f"Command launched with PID {pid}. "
                    f"Track progress with check_status('{result_dir}')."
                ),
            })
        except ValueError as e:
            return json.dumps({
                "status": "error",
                "message": str(e),
            })

    return handler


def register_all_tools(mcp_server: Any) -> None:
    """Register all CLI + stub tools on the MCP server instance."""
    from phyloai.cli.main import cli

    # Walk Click tree to get all implemented tools
    descriptors = walk_click_tree(cli)
    dynamic_names = {d["tool_name"] for d in descriptors}

    # Register dynamic tools from Click tree
    for desc in descriptors:
        tool_def = build_mcp_tool(desc)
        handler = _make_launch_handler(desc)
        mcp_server.tool(
            name=tool_def["name"],
            description=tool_def["description"],
        )(handler)

    # Register stub tools only for names NOT in the dynamic Click tree
    for stub_def in STUB_TOOLS:
        if stub_def["name"] in dynamic_names:
            continue  # Already registered as a real tool above
        async def _stub_handler(name=stub_def["name"], **kwargs: Any) -> str:
            result = handle_stub(name)
            return json.dumps(result or {})
        mcp_server.tool(
            name=stub_def["name"],
            description=stub_def["description"],
        )(_stub_handler)
```

- [ ] **Step 5: Run tests to verify they pass**

Run: `python -m pytest tests/mcp/test_cli_tools.py -v`
Expected: All tests PASS.

- [ ] **Step 6: Commit**

```bash
git add phyloai/mcp/tools/cli_tools.py tests/mcp/test_cli_tools.py
git commit -m "feat(mcp): add fire-and-forget CLI tool wrappers"
```

---

### Task 7: MCP Server — stdio server entry point and integration

**Files:**
- Create: `phyloai/mcp/server.py`
- Modify: `phyloai/cli/commands/mcp_server.py` (replace inert stub)

**Interfaces:**
- Consumes: `phyloai.mcp.tools.cli_tools` (register_all_tools), `phyloai.mcp.tools.utils` (check_status, read_result, read_report, get_command_schema)
- Produces: `main()` async function — creates MCP Server, registers all tools, runs stdio loop

- [ ] **Step 1: Implement `phyloai/mcp/server.py`**

```python
"""PhyloAI MCP Server — stdio transport, tool registration."""
from __future__ import annotations

import json
import sys

from mcp.server import Server
from mcp.server.stdio import stdio_server

from phyloai.mcp.tools.cli_tools import register_all_tools
from phyloai.mcp.tools.utils import check_status, get_command_schema, read_report, read_result


def create_server() -> Server:
    """Create and configure the PhyloAI MCP server with all tools registered."""
    server = Server("phyloai")

    # Register utility tools manually (they don't map to CLI commands)
    @server.tool(
        name="check_status",
        description="Check the progress of a running phyloai job by inspecting its output directory. Returns state (not_started|running|success|error|unknown) and partial results.",
    )
    async def _check_status(output_dir: str) -> str:
        return json.dumps(check_status(output_dir))

    @server.tool(
        name="read_result",
        description="Read result.json from a step output directory.",
    )
    async def _read_result(output_dir: str) -> str:
        return json.dumps(read_result(output_dir))

    @server.tool(
        name="read_report",
        description="Read report.json from a run directory. The report aggregates all step records into one document.",
    )
    async def _read_report(run_dir: str) -> str:
        return json.dumps(read_report(run_dir))

    @server.tool(
        name="get_command_schema",
        description="Get the MCP tool schema (parameter names, types, defaults, choices, help) for a phyloai command by its tool name.",
    )
    async def _get_command_schema(command_name: str) -> str:
        return json.dumps(get_command_schema(command_name))

    # Register dynamic CLI tools + stubs
    register_all_tools(server)

    return server


async def main() -> None:
    """Entry point: create server, connect stdio transport, run."""
    server = create_server()
    async with stdio_server() as (reader, writer):
        await server.run(reader, writer, server.create_initialization_options())


def entry() -> None:
    """Synchronous wrapper for the async main."""
    import asyncio
    try:
        asyncio.run(main())
    except KeyboardInterrupt:
        pass
    except Exception as e:
        print(f"MCP server error: {e}", file=sys.stderr)
        sys.exit(1)
```

- [ ] **Step 2: Replace inert stub in `phyloai/cli/commands/mcp_server.py`**

```python
"""phyloai mcp-server — start the MCP protocol server."""
from __future__ import annotations

import click


@click.command(name="mcp-server", hidden=False)
def mcp_server() -> None:
    """Start the PhyloAI MCP server (stdio transport).

    Configure your AI client to run:

        phyloai mcp-server

    The server registers tools for all CLI commands and utility operations.
    """
    from phyloai.mcp.server import entry
    entry()
```

- [ ] **Step 3: Write integration test — server object can be created**

```python
# tests/mcp/conftest.py — add to existing file or create
import pytest
from phyloai.mcp.server import create_server


@pytest.fixture
def mcp_server():
    return create_server()
```

```python
# tests/mcp/test_cli_tools.py — add test
def test_server_creates_without_error(mcp_server):
    """Server object creation should not raise exceptions."""
    assert mcp_server is not None
    assert hasattr(mcp_server, "tool")
```

- [ ] **Step 4: Run all MCP tests**

Run: `python -m pytest tests/mcp/ -v`
Expected: All tests PASS.

- [ ] **Step 5: Verify `phyloai mcp-server --help`**

Run: `python -m phyloai.cli.main mcp-server --help`
Expected: Help text displayed without errors.

- [ ] **Step 6: Commit**

```bash
git add phyloai/mcp/server.py phyloai/cli/commands/mcp_server.py tests/mcp/conftest.py tests/mcp/test_cli_tools.py
git commit -m "feat(mcp): implement sysentry point with utility and CLI tool registration"
```

---

### Task 8: Skill — Create phyloai-workflow SKILL.md and reference files

**Files:**
- Create: `skills/phyloai-workflow/SKILL.md`
- Create: `skills/phyloai-workflow/references/parameter-annotations.md`
- Create: `skills/phyloai-workflow/references/error-catalog.md`
- Create: `skills/phyloai-workflow/references/dialog-templates.md`
- Create: `skills/phyloai-workflow/references/demo-data.md`
- Create: `skills/phyloai-workflow/references/workflow.md`

**Note:** These are markdown documentation files, not code. They do not need tests, but must be complete and aligned with the design spec.

- [ ] **Step 1: Create `skills/phyloai-workflow/SKILL.md`**

Write the main Skill document with these sections:
- **Core Rules:** phyloai doctor gate (required before tool-invoking commands, not required for read-only ops like check_status/read_result/read_report). Every command requires explicit approval — show parameter card, wait for confirmation. After each result, summarize outputs first, wait for user before next step. Never auto-run. Dedicated run root under runs/.
- **Parameter Source-of-Truth Protocol:** Load schema from MCP `get_command_schema` at card render time (never from memory). If MCP is unreachable, state the fallback explicitly. Never invent parameter names/aliases/enum values.
- **Parameter Validation Gate:** Before execution, validate approved params against schema. Unknown parameter → block. Missing required → block. Enum out of range → block. On block, show fix card with invalid field and valid alternatives.
- **Data Source Decision:** User data always takes priority. Demo dataset available for teaching/troubleshooting. After demo, explicitly switch back to user paths.
- **Entry Modes:** New task, Resume task, Single-step mode. Auto-detect from user's first message.
- **Session Recovery:** `read_report(run_dir)` → if missing, ask user to run `phyloai report` or specify step dir → `read_result(step_dir)` as fallback.
- **Language Policy:** Parameter cards - English param name + Chinese annotation + recommended value. Conversation and result interpretation - follow user language.
- **Phases:** Phase 0: doctor. Phase 1: pretree (convert → stats → align → trim → metrics → filter → concat). Phase 2: tree (ml/bi/msc/cf). Phase 3: posttree (topology/dating). Phase 4: report.
- **Error Handling:** Exit 1/3 → catalog lookup → fix card. Exit 2 → pass tool_stderr to AI; for batch, only failed loci logs, capped at 10 loci.
- **References:** List all reference files.

- [ ] **Step 2: Create `skills/phyloai-workflow/references/parameter-annotations.md`**

For each CLI command, list all parameters with:
- English parameter name
- Chinese description of what it does and scientific meaning
- Recommended value (where applicable)
- Common pitfalls

Start with the most commonly-used commands: doctor, pretree align, pretree trim, pretree concat, pretree filter, tree ml iqtree, tree ml fasttree, tree msc, tree cf, posttree topology, posttree dating hessian, posttree dating mcmc, report, run.

Each entry format:
```
### pretree align --method
对比算法。
- linsi: 最精确，使用一致性迭代比对，适合<500个序列 [推荐: 小型数据集]
- auto: 自动选择策略，适合大型数据集 [推荐: 大型数据集]
- ginsi: 全局比对（序列长度差异大时避免使用）
```

- [ ] **Step 3: Create `skills/phyloai-workflow/references/error-catalog.md`**

Cover known error patterns for exit 1 (user input) and exit 3 (environment):

```
### Exit 1: Missing required parameter
Pattern: "Error: --matrix is required"
Fix: Provide the parameter. Show the parameter card for the command.

### Exit 1: Output directory exists
Pattern: "Output directory exists and is not empty"
Fix: Add --overwrite to replace, or use a different --output-dir.

### Exit 3: Required tool not installed
Pattern: "Error: iqtree3 not found"
Fix: Run `phyloai doctor` to check environment. Install the missing tool or specify its path via --xxx-path.
...
```

Minimum 10 catalog entries covering: missing required param, output dir conflict, invalid enum, file not found, empty file, unsupported format, tool not installed, insufficient threads, --overwrite/--resume conflict, --tool-args blocked flag.

- [ ] **Step 4: Create `skills/phyloai-workflow/references/dialog-templates.md`**

Templates for:
1. **Pre-run parameter card:** Command name, purpose (1 line), parameter table (English name + Chinese description + value), schema source line, confirmation question.
2. **Post-run result card:** Status, wall time, key results summary, warnings, next step suggestions.
3. **Doctor result card:** Tool table with status indicators.
4. **Session recovery card:** Steps completed / running / failed table, next step suggestion.
5. **Demo mode prompt:** One-line offer.

- [ ] **Step 5: Create `skills/phyloai-workflow/references/demo-data.md`**

Document the demo dataset structure and per-step entry points:
```
## End-to-end dataset
- Path: phyloai/demo_data/end_to_end/raw/
- Content: 10 genes, 20 taxa, small protein sequences
- Pipeline: raw → align → trim → concat → iqtree → report (runs in ~2-5 min)

## Per-step entry points
- Pre-aligned: phyloai/demo_data/per_step/aligned/
- Pre-trimmed: phyloai/demo_data/per_step/trimmed/
- Gene trees: phyloai/demo_data/per_step/gen_trees/
- Concatenated: phyloai/demo_data/per_step/concat/matrix.fa

Usage: copy demo data to a user run directory before starting.
```

- [ ] **Step 6: Create `skills/phyloai-workflow/references/workflow.md`**

Per-phase execution reference:
- Phase 0: `phyloai doctor` — check environment, present tool status table
- Phase 1 (pretree): `convert → stats → align → trim → metrics → filter taper → concat`
- Phase 2 (tree): `ml iqtree` or `ml fasttree`, then optional `msc`, `cf`
- Phase 3 (posttree): `topology` (if multiple trees), `dating hessian → dating mcmc`
- Phase 4 (report): `phyloai report --run-dir <run_dir>`

Each phase lists: input requirements, what to check after completion, common pitfalls, next steps.

- [ ] **Step 7: Commit**

```bash
git add skills/
git commit -m "docs(skill): add phyloai-workflow Skill and reference files"
```

---

### Task 9: Demo Dataset — Create and bundle small phylogenomics dataset

**Files:**
- Create: `phyloai/demo_data/__init__.py`
- Create: `phyloai/demo_data/end_to_end/raw/gene_001.faa` through `gene_010.faa`
- Create: `phyloai/demo_data/per_step/aligned/.gitkeep`
- Create: `phyloai/demo_data/per_step/trimmed/.gitkeep`
- Create: `phyloai/demo_data/per_step/gen_trees/.gitkeep`

- [ ] **Step 1: Create `phyloai/demo_data/__init__.py`**

```python
"""Demo dataset path resolution."""
from __future__ import annotations

from pathlib import Path


def resolve_demo_path(*parts: str) -> Path:
    """Return the absolute path to a file or directory within the bundled demo dataset."""
    base = Path(__file__).resolve().parent
    return base.joinpath(*parts)


def resolve_raw_dir() -> Path:
    return resolve_demo_path("end_to_end", "raw")


def resolve_per_step_dir(step: str) -> Path:
    return resolve_demo_path("per_step", step)
```

- [ ] **Step 2: Generate 10 minimal gene sequences (AA FASTA)**

Each file contains 15-20 representative taxa with 50-200 amino acids. Use varied sequence lengths (50-200 aa) to reflect real datasets. Ensure all FASTA files are valid (60-character line wrapping).  Create small between-gene variation so that downstream steps produce non-trivial results.

```python
# Example: gene_001.faa
>T1
MALWMRLLPL-LALLALWGPDPAAAFVNQHL-CGSHLVEALYLVCGERGFFYTPKT
>T2
MALWIRLLPL-LALLALWGPDPAAAFVNQHL-CGSHLVEALYLVCGERGFFYTPKT
>T3
MALWMRLLPL-VLLALWGPDPAAAFVNQHL-CGSHLVEALYLVCGERGFFYTPKT
...
```

Create 10 gene files (gene_001.faa through gene_010.faa) with 15 taxa each.

- [ ] **Step 3: Generate per-step intermediate data**

Run each step on the demo dataset to produce actual intermediate output:

```bash
# Convert
python -m phyloai.cli.main pretree convert --input phyloai/demo_data/end_to_end/raw -o phyloai/demo_data/per_step/converted --overwrite

# Align (use mafft auto for speed)
python -m phyloai.cli.main pretree align --seq-dir phyloai/demo_data/per_step/converted/seqs --method auto -o phyloai/demo_data/per_step/aligned --overwrite

# Trim
python -m phyloai.cli.main pretree trim --msa-dir phyloai/demo_data/per_step/aligned/seqs --tool clipkit -o phyloai/demo_data/per_step/trimmed --overwrite

# Concat
python -m phyloai.cli.main pretree concat --msa-dir phyloai/demo_data/per_step/trimmed/seqs -o phyloai/demo_data/per_step/concat --overwrite
```

Move only the output sequence files (seqs/) into the per_step directories for commit.  Keep result.json in each directory for completeness.  Remove large generated files (e.g. IQ-TREE checkpoints) — only keep aligned/trimmed sequences and concatenated matrix.

- [ ] **Step 4: Verify demo data structure**

Run: `python -c "from phyloai.demo_data import resolve_demo_path, resolve_raw_dir; print(resolve_raw_dir()); print(list(f.suffix for f in resolve_raw_dir().glob('*.faa')))"`

Expected: Lists the demo_data directory and shows 10 `.faa` files.

- [ ] **Step 5: Verify FASTA files are valid by running pretree stats on them**

Run: `python -m phyloai.cli.main pretree stats --seq-dir phyloai/demo_data/end_to_end/raw -o runs/demo/stats`
Expected: Success, `result.json` written with valid stats.

- [ ] **Step 6: Commit**

```bash
git add phyloai/demo_data/
git commit -m "feat: add bundled demo dataset with 10 genes, 20 taxa"
```

---

### Task 10: Final Verification — End-to-end pipeline test with demo dataset

**Files:**
- None (verification only)

- [ ] **Step 1: Run doctor on demo environment**

Run: `python -m phyloai.cli.main doctor`
Expected: All tools report OK or explain missing tools.

- [ ] **Step 2: Verify existing per-step demo data**

Run: `python -m phyloai.cli.main pretree stats --seq-dir phyloai/demo_data/per_step/aligned/seqs -o runs/demo/verify-stats --overwrite`
Expected: Success, valid stats for aligned demo sequences.

- [ ] **Step 3: Run full MCP test suite**

Run: `python -m pytest tests/mcp/ -v`
Expected: All MCP tests PASS.

- [ ] **Step 4: Commit if any changes**

```bash
git status
# If needed: git add ... && git commit -m "chore: verification fixes"
```

---

### Self-Review Checklist

1. **Spec coverage:** Verify each requirement maps to a task:
   - [x] MCP tool granularity (one per subcommand) → Task 2 (schema_gen), Task 6 (cli_tools)
   - [x] Schema from Click import → Task 2
   - [x] Fire-and-forget + output_dir as job handle → Task 3, Task 4 (check_status), Task 6
   - [x] job.json with pid, started_at, command → Task 3
   - [x] check_status with not_started|running|success|error|unknown → Task 4
   - [x] Pre-launch output_dir validation → Task 3 (launch_cli), Task 4
   - [x] MCP always passes --output-dir explicitly → Task 3 (build_cli_argv)
   - [x] Utility tools (check_status, read_result, read_report, get_command_schema) → Task 4
   - [x] get_command_schema exposes type, required, default, choices, help, path hints → Task 2
   - [x] Stub tools for unimplemented commands → Task 5
   - [x] Stubs: "not yet available in the installed version" message → Task 5
   - [x] stdio transport → Task 7
   - [x] Skill lives in skills/phyloai-workflow/ → Task 8
   - [x] Parameter cards: English name + Chinese annotation + recommended value → Task 8 (parameter-annotations.md)
   - [x] Session recovery with fallback (read_report → report → read_result) → Task 8 (SKILL.md)
   - [x] Doctor required before tool-invoking commands, not read-only ops → Task 8 (SKILL.md)
   - [x] Exit 1/3 catalog-driven, exit 2 AI free-form → Task 8 (error-catalog.md, SKILL.md)
   - [x] Batch failures: only failed loci logs, capped → Task 8 (SKILL.md)
   - [x] Demo dataset bundled in package → Task 9
   - [x] Skill version coupled to CLI, same repo → Task 8 (in-repo), Task 1 (same package)

2. **Placeholder scan:** No TBD, TODO, "implement later", "fill in details". All tasks have concrete code or content.

3. **Type consistency:** Tool names follow pattern `{group}_{subcommand}` (e.g., `pretree_align`, `tree_ml_iqtree`, `posttree_dating_hessian`). `launch_cli` returns `tuple[Path, int]`. `check_status` returns `dict` with `status` key. All interfaces defined.

4. **Missing spec requirements:**
   - `polish_methods` tool: Deferred to future, not in current scope. Listed in "Future Work" section of relevant spec.
   - HTTP transport: Deferred, not in scope.
   - `phyloai-syserror` Skill: Deferred to Phase 9, not in current scope.
   - Report AI review: Deferred to Phase 10, not in current scope.

All current-scope requirements are covered.
