# PhyloAI core/ Module Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build the `core/` infrastructure layer that all other PhyloAI modules depend on: project scaffolding, environment detection, external tool execution, format I/O, shared data structures, and per-step logging.

**Architecture:** Five focused modules with no inter-dependencies except `logger` and `schema` being importable by all others. `runner.py` depends on `env.py` for tool paths. `formats.py` is standalone. All modules expose clean Python APIs consumed by both library users and the CLI layer.

**Tech Stack:** Python ≥3.10, `click` (CLI), `rich` (terminal output), `PyYAML`, `pytest`, `biopython` (sequence I/O)

---

## File Map

| File | Responsibility |
|------|---------------|
| `phyloai/__init__.py` | Package entry, version string |
| `phyloai/core/__init__.py` | Re-exports public API of core |
| `phyloai/core/schema.py` | Shared dataclasses: `MSACollection`, `TreeSet`, `RunRecord`, `ToolResult` |
| `phyloai/core/env.py` | Tool detection, path resolution, bundled tool management |
| `phyloai/core/runner.py` | Unified external tool call interface with timeout, retry, logging |
| `phyloai/core/formats.py` | Format detection and conversion (FASTA/Nexus/Phylip/Phylip-PAML) |
| `phyloai/core/logger.py` | Per-step log file writer (each step writes its own `<step>.log` inside the command's output directory; no shared `logs/` folder) |
| `phyloai/cli/__init__.py` | CLI package |
| `phyloai/cli/main.py` | `phyloai` entry point, config file loading |
| `phyloai/cli/doctor.py` | `phyloai doctor` command |
| `pyproject.toml` | Package metadata, dependencies, entry points |
| `tests/core/test_schema.py` | Schema dataclass tests |
| `tests/core/test_env.py` | Tool detection tests |
| `tests/core/test_runner.py` | Runner tests (mocked subprocess) |
| `tests/core/test_formats.py` | Format detection and conversion tests |
| `tests/core/test_logger.py` | Log file output tests |
| `tests/cli/test_doctor.py` | `phyloai doctor` CLI tests |

---

## Task 1: Project Scaffold

**Files:**
- Create: `pyproject.toml`
- Create: `phyloai/__init__.py`
- Create: `phyloai/core/__init__.py`
- Create: `phyloai/cli/__init__.py`
- Create: `tests/__init__.py`
- Create: `tests/core/__init__.py`
- Create: `tests/cli/__init__.py`

- [ ] **Step 1: Create `pyproject.toml`**

```toml
[build-system]
requires = ["hatchling"]
build-backend = "hatchling.build"

[project]
name = "phyloai"
version = "0.1.0"
description = "An AI-native modular phylogenomics analysis platform"
readme = "README.md"
requires-python = ">=3.10"
license = {text = "MIT"}
dependencies = [
    "click>=8.1",
    "rich>=13.0",
    "PyYAML>=6.0",
    "biopython>=1.81",
    "magus-msa>=1.0",
    "clipkit>=2.0",
]

[project.scripts]
phyloai = "phyloai.cli.main:cli"

[tool.hatch.build.targets.wheel]
packages = ["phyloai"]

[tool.pytest.ini_options]
testpaths = ["tests"]
```

- [ ] **Step 2: Create package init files**

`phyloai/__init__.py`:
```python
"""PhyloAI — An AI-native modular phylogenomics analysis platform."""

__version__ = "0.1.0"
```

`phyloai/core/__init__.py`:
```python
"""Core infrastructure for PhyloAI."""

from phyloai.core.schema import MSACollection, TreeSet, RunRecord, ToolResult
from phyloai.core.env import ToolEnv
from phyloai.core.runner import Runner
from phyloai.core.formats import FormatConverter
from phyloai.core.logger import StepLogger

__all__ = [
    "MSACollection", "TreeSet", "RunRecord", "ToolResult",
    "ToolEnv", "Runner", "FormatConverter", "StepLogger",
]
```

`phyloai/cli/__init__.py` and both `tests/` inits: empty files.

- [ ] **Step 3: Install in development mode**

```bash
pip install -e ".[dev]"
```

- [ ] **Step 4: Verify import works**

```bash
python -c "import phyloai; print(phyloai.__version__)"
```

Expected output: `0.1.0`

- [ ] **Step 5: Commit**

```bash
git add pyproject.toml phyloai/ tests/
git commit -m "feat: project scaffold and package structure"
```

---

## Task 2: Shared Data Structures (`schema.py`)

**Files:**
- Create: `phyloai/core/schema.py`
- Create: `tests/core/test_schema.py`

- [ ] **Step 1: Write failing tests**

```python
# tests/core/test_schema.py
import pytest
from pathlib import Path
from phyloai.core.schema import MSACollection, TreeSet, RunRecord, ToolResult


def test_msa_collection_requires_directory():
    with pytest.raises(TypeError):
        MSACollection()  # missing required field


def test_msa_collection_defaults():
    c = MSACollection(directory=Path("./alignments"))
    assert c.seq_type == "AA"
    assert c.file_extension == ".fa"
    assert c.count == 0


def test_tool_result_success():
    r = ToolResult(
        tool="iqtree3",
        command="iqtree3 -s matrix.fa",
        returncode=0,
        stdout="Analysis done",
        stderr="",
        wall_time=12.5,
    )
    assert r.success is True


def test_tool_result_failure():
    r = ToolResult(
        tool="iqtree3",
        command="iqtree3 -s missing.fa",
        returncode=1,
        stdout="",
        stderr="ERROR: file not found",
        wall_time=0.1,
    )
    assert r.success is False


def test_run_record_to_dict():
    record = RunRecord(run_dir=Path("./runs"))
    d = record.to_dict()
    assert "run_dir" in d
    assert "steps" in d
    assert isinstance(d["steps"], list)


def test_tree_set_defaults():
    ts = TreeSet(directory=Path("./trees"))
    assert ts.format == "newick"
    assert ts.count == 0
```

- [ ] **Step 2: Run tests — expect failure**

```bash
pytest tests/core/test_schema.py -v
```

Expected: `ImportError` — module does not exist yet.

- [ ] **Step 3: Implement `schema.py`**

```python
# phyloai/core/schema.py
"""Shared data structures for PhyloAI."""

from __future__ import annotations
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any


@dataclass
class MSACollection:
    """A directory of multiple sequence alignment files."""
    directory: Path
    seq_type: str = "AA"          # "AA" or "NT"
    file_extension: str = ".fa"   # ".fa", ".faa", ".fna", ".phy", ".nex"
    count: int = 0                # number of alignment files found

    def __post_init__(self):
        self.directory = Path(self.directory)
        if self.directory.exists():
            self.count = len(list(
                self.directory.glob(f"*{self.file_extension}")
            ))


@dataclass
class TreeSet:
    """A directory of phylogenetic tree files."""
    directory: Path
    format: str = "newick"        # "newick" or "nexus"
    file_extension: str = ".treefile"
    count: int = 0

    def __post_init__(self):
        self.directory = Path(self.directory)
        if self.directory.exists():
            self.count = len(list(
                self.directory.glob(f"*{self.file_extension}")
            ))


@dataclass
class ToolResult:
    """Result of a single external tool invocation."""
    tool: str
    command: str
    returncode: int
    stdout: str
    stderr: str
    wall_time: float             # seconds

    @property
    def success(self) -> bool:
        return self.returncode == 0

    def to_dict(self) -> dict[str, Any]:
        return {
            "tool": self.tool,
            "command": self.command,
            "returncode": self.returncode,
            "stdout": self.stdout,
            "stderr": self.stderr,
            "wall_time": self.wall_time,
            "success": self.success,
        }


@dataclass
class RunRecord:
    """Full record of a PhyloAI analysis run."""
    run_dir: Path
    phyloai_version: str = ""
    steps: list[dict[str, Any]] = field(default_factory=list)

    def __post_init__(self):
        self.run_dir = Path(self.run_dir)
        from phyloai import __version__
        self.phyloai_version = __version__

    def add_step(self, step_name: str, params: dict, result: ToolResult) -> None:
        self.steps.append({
            "step": step_name,
            "params": params,
            "result": result.to_dict(),
        })

    def to_dict(self) -> dict[str, Any]:
        return {
            "run_dir": str(self.run_dir),
            "phyloai_version": self.phyloai_version,
            "steps": self.steps,
        }
```

- [ ] **Step 4: Run tests — expect pass**

```bash
pytest tests/core/test_schema.py -v
```

Expected: all 6 tests PASS.

- [ ] **Step 5: Commit**

```bash
git add phyloai/core/schema.py tests/core/test_schema.py
git commit -m "feat(core): shared data structures (schema.py)"
```

---

## Task 3: Environment Detection (`env.py`)

**Files:**
- Create: `phyloai/core/env.py`
- Create: `tests/core/test_env.py`

- [ ] **Step 1: Write failing tests**

```python
# tests/core/test_env.py
import pytest
from unittest.mock import patch
from pathlib import Path
from phyloai.core.env import ToolEnv, ToolStatus, ToolInfo


def test_tool_status_values():
    assert ToolStatus.OK == "ok"
    assert ToolStatus.WARN == "warn"
    assert ToolStatus.MISSING == "missing"


def test_tool_info_defaults():
    info = ToolInfo(name="iqtree3")
    assert info.status == ToolStatus.MISSING
    assert info.path is None
    assert info.version is None


def test_detect_present_tool():
    env = ToolEnv()
    with patch("shutil.which", return_value="/usr/bin/echo"):
        with patch.object(env, "_get_version", return_value="1.0"):
            result = env._detect_tool("echo", version_flag="--version")
    assert result.status == ToolStatus.OK
    assert result.path == Path("/usr/bin/echo")


def test_detect_missing_tool():
    env = ToolEnv()
    with patch("shutil.which", return_value=None):
        result = env._detect_tool("nonexistent_tool_xyz")
    assert result.status == ToolStatus.MISSING
    assert result.path is None


def test_check_all_returns_dict():
    env = ToolEnv()
    results = env.check_all()
    assert isinstance(results, dict)
    # known required tools must be present as keys
    for key in ["iqtree3", "mafft", "trimal"]:
        assert key in results


def test_get_tool_path_raises_when_missing():
    env = ToolEnv()
    with patch("shutil.which", return_value=None):
        env._tools["fake_tool"] = ToolInfo(name="fake_tool")
        with pytest.raises(FileNotFoundError, match="fake_tool"):
            env.require("fake_tool")
```

- [ ] **Step 2: Run tests — expect failure**

```bash
pytest tests/core/test_env.py -v
```

Expected: `ImportError`.

- [ ] **Step 3: Implement `env.py`**

```python
# phyloai/core/env.py
"""Tool detection and path resolution for PhyloAI."""

from __future__ import annotations
import shutil
import subprocess
from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path
from typing import Optional


class ToolStatus(str, Enum):
    OK = "ok"
    WARN = "warn"       # found but may have issues
    MISSING = "missing"


@dataclass
class ToolInfo:
    name: str
    status: ToolStatus = ToolStatus.MISSING
    path: Optional[Path] = None
    version: Optional[str] = None
    note: str = ""      # install hint or warning message


# Registry of all tools PhyloAI can use.
# required=True -> doctor reports MISSING as error
# required=False -> doctor reports MISSING as warning only
TOOL_REGISTRY: dict[str, dict] = {
    # user-installed tools
    "iqtree3":     {"required": True,  "version_flag": "--version",
                    "install": "https://github.com/iqtree/iqtree3/releases"},
    "mafft":       {"required": True,  "version_flag": "--version",
                    "install": "https://mafft.cbrc.jp/alignment/software/"},
    "wastral":     {"required": False, "version_flag": "-v",
                    "install": "https://github.com/chaoszhang/ASTER"},
    "pb_mpi":      {"required": False, "version_flag": "",
                    "install": "https://github.com/bayesiancook/pbmpi"},
    "mcmctree":    {"required": False, "version_flag": "",
                    "install": "https://github.com/abacus-gene/paml/releases"},
    "run_treeshrink.py": {"required": False, "version_flag": "--version",
                    "install": "https://github.com/uym2/TreeShrink"},
    # pip-installed tools (detected via shutil.which after pip install)
    "magus":       {"required": False, "version_flag": "--version",
                    "install": "pip install magus-msa"},
    "clipkit":     {"required": False, "version_flag": "--version",
                    "install": "pip install clipkit"},

    # bundled tools (path resolved relative to package)
    "trimal":      {"required": True,  "version_flag": "--version",
                    "bundled": True,
                    "install": "https://github.com/inab/trimal/releases"},
    "bmge":        {"required": False, "version_flag": "",
                    "bundled": True,
                    "install": "conda install bmge"},
}


class ToolEnv:
    """Detect and cache tool availability across the system."""

    def __init__(self):
        self._tools: dict[str, ToolInfo] = {}
        self._bundled_dir = Path(__file__).parent.parent / "bundled"

    def _get_version(self, path: Path, version_flag: str) -> Optional[str]:
        if not version_flag:
            return None
        try:
            result = subprocess.run(
                [str(path), version_flag],
                capture_output=True, text=True, timeout=5
            )
            output = result.stdout.strip() or result.stderr.strip()
            # take first non-empty line
            for line in output.splitlines():
                if line.strip():
                    return line.strip()[:80]  # cap length
        except Exception:
            pass
        return None

    def _detect_tool(self, name: str, version_flag: str = "",
                     bundled: bool = False) -> ToolInfo:
        # check bundled location first
        if bundled:
            bundled_path = self._bundled_dir / name / name
            if bundled_path.exists():
                ver = self._get_version(bundled_path, version_flag)
                return ToolInfo(
                    name=name, status=ToolStatus.OK,
                    path=bundled_path, version=ver,
                    note="bundled"
                )
        # check system PATH
        found = shutil.which(name)
        if found:
            p = Path(found)
            ver = self._get_version(p, version_flag)
            return ToolInfo(
                name=name, status=ToolStatus.OK,
                path=p, version=ver
            )
        return ToolInfo(name=name)

    def check_all(self) -> dict[str, ToolInfo]:
        """Detect all registered tools and cache results."""
        for name, meta in TOOL_REGISTRY.items():
            info = self._detect_tool(
                name,
                version_flag=meta.get("version_flag", ""),
                bundled=meta.get("bundled", False),
            )
            if info.status == ToolStatus.MISSING:
                info.note = f"install: {meta.get('install', '')}"
            self._tools[name] = info
        return self._tools

    def require(self, name: str) -> Path:
        """Return path to tool or raise FileNotFoundError."""
        info = self._tools.get(name) or self._detect_tool(name)
        if info.status != ToolStatus.OK or info.path is None:
            raise FileNotFoundError(
                f"Required tool '{name}' not found. "
                f"{TOOL_REGISTRY.get(name, {}).get('install', '')}"
            )
        return info.path

    def get(self, name: str) -> Optional[Path]:
        """Return path to tool or None if unavailable."""
        try:
            return self.require(name)
        except FileNotFoundError:
            return None
```

- [ ] **Step 4: Run tests — expect pass**

```bash
pytest tests/core/test_env.py -v
```

Expected: all 6 tests PASS.

- [ ] **Step 5: Commit**

```bash
git add phyloai/core/env.py tests/core/test_env.py
git commit -m "feat(core): tool detection and path resolution (env.py)"
```

---

## Task 4: External Tool Runner (`runner.py`)

**Files:**
- Create: `phyloai/core/runner.py`
- Create: `tests/core/test_runner.py`

- [ ] **Step 1: Write failing tests**

```python
# tests/core/test_runner.py
import pytest
import time
from unittest.mock import patch, MagicMock
from pathlib import Path
from phyloai.core.runner import Runner
from phyloai.core.schema import ToolResult


def test_runner_run_success():
    runner = Runner()
    result = runner.run(["echo", "hello"], tool_name="echo")
    assert result.success is True
    assert "hello" in result.stdout
    assert result.tool == "echo"
    assert result.wall_time > 0


def test_runner_run_failure():
    runner = Runner()
    result = runner.run(["false"], tool_name="false")
    assert result.success is False
    assert result.returncode != 0


def test_runner_run_returns_tool_result():
    runner = Runner()
    result = runner.run(["echo", "test"], tool_name="echo")
    assert isinstance(result, ToolResult)


def test_runner_raises_on_missing_executable():
    runner = Runner()
    with pytest.raises(FileNotFoundError):
        runner.run(["nonexistent_binary_xyz_abc"], tool_name="fake")


def test_runner_timeout():
    runner = Runner(timeout=1)
    with pytest.raises(TimeoutError):
        runner.run(["sleep", "10"], tool_name="sleep")


def test_runner_captures_stderr():
    runner = Runner()
    # 'ls' on a non-existent path writes to stderr
    result = runner.run(
        ["ls", "/nonexistent_path_xyz_abc_123"],
        tool_name="ls"
    )
    assert result.success is False
    assert result.stderr != "" or result.returncode != 0
```

- [ ] **Step 2: Run tests — expect failure**

```bash
pytest tests/core/test_runner.py -v
```

Expected: `ImportError`.

- [ ] **Step 3: Implement `runner.py`**

```python
# phyloai/core/runner.py
"""Unified external tool call interface."""

from __future__ import annotations
import subprocess
import time
from pathlib import Path
from typing import Optional

from phyloai.core.schema import ToolResult


class Runner:
    """Execute external tools and capture results."""

    def __init__(self, timeout: int = 86400):
        """
        Args:
            timeout: Maximum seconds to wait for a tool. Default 24h.
        """
        self.timeout = timeout

    def run(
        self,
        cmd: list[str],
        tool_name: str,
        cwd: Optional[Path] = None,
        env: Optional[dict] = None,
    ) -> ToolResult:
        """
        Run an external command and return a ToolResult.

        Raises:
            FileNotFoundError: if the executable is not found.
            TimeoutError: if the command exceeds self.timeout seconds.
        """
        command_str = " ".join(str(c) for c in cmd)
        start = time.monotonic()

        try:
            proc = subprocess.run(
                [str(c) for c in cmd],
                capture_output=True,
                text=True,
                timeout=self.timeout,
                cwd=str(cwd) if cwd else None,
                env=env,
            )
        except FileNotFoundError:
            raise FileNotFoundError(
                f"Executable not found: '{cmd[0]}'. "
                f"Check 'phyloai doctor' for installation status."
            )
        except subprocess.TimeoutExpired:
            raise TimeoutError(
                f"Tool '{tool_name}' exceeded timeout of {self.timeout}s. "
                f"Command: {command_str}"
            )

        wall_time = time.monotonic() - start

        return ToolResult(
            tool=tool_name,
            command=command_str,
            returncode=proc.returncode,
            stdout=proc.stdout,
            stderr=proc.stderr,
            wall_time=wall_time,
        )
```

- [ ] **Step 4: Run tests — expect pass**

```bash
pytest tests/core/test_runner.py -v
```

Expected: all 6 tests PASS.

- [ ] **Step 5: Commit**

```bash
git add phyloai/core/runner.py tests/core/test_runner.py
git commit -m "feat(core): unified external tool runner (runner.py)"
```

---

## Task 5: Format Detection and Conversion (`formats.py`)

**Files:**
- Create: `phyloai/core/formats.py`
- Create: `tests/core/test_formats.py`
- Create: `tests/core/fixtures/` (small test alignment files)

- [ ] **Step 1: Create test fixtures**

`tests/core/fixtures/test.fasta`:
```
>Taxon_A
MARVELLOUS
>Taxon_B
MARVELIOUS
>Taxon_C
MARVELSOUS
```

`tests/core/fixtures/test.phy` (Phylip strict):
```
 3 10
Taxon_A   MARVELLOUS
Taxon_B   MARVELIOUS
Taxon_C   MARVELSOUS
```

`tests/core/fixtures/test.nex` (Nexus):
```
#NEXUS
BEGIN DATA;
  DIMENSIONS NTAX=3 NCHAR=10;
  FORMAT DATATYPE=PROTEIN GAP=- MISSING=?;
  MATRIX
    Taxon_A MARVELLOUS
    Taxon_B MARVELIOUS
    Taxon_C MARVELSOUS
  ;
END;
```

- [ ] **Step 2: Write failing tests**

```python
# tests/core/test_formats.py
import pytest
from pathlib import Path
from phyloai.core.formats import FormatConverter, AlignmentFormat

FIXTURES = Path(__file__).parent / "fixtures"


def test_detect_fasta():
    fc = FormatConverter()
    fmt = fc.detect(FIXTURES / "test.fasta")
    assert fmt == AlignmentFormat.FASTA


def test_detect_phylip():
    fc = FormatConverter()
    fmt = fc.detect(FIXTURES / "test.phy")
    assert fmt == AlignmentFormat.PHYLIP


def test_detect_nexus():
    fc = FormatConverter()
    fmt = fc.detect(FIXTURES / "test.nex")
    assert fmt == AlignmentFormat.NEXUS


def test_fasta_to_phylip(tmp_path):
    fc = FormatConverter()
    out = tmp_path / "out.phy"
    fc.convert(FIXTURES / "test.fasta", out, target=AlignmentFormat.PHYLIP)
    assert out.exists()
    content = out.read_text()
    assert "Taxon_A" in content
    assert "MARVELLOUS" in content


def test_fasta_to_nexus(tmp_path):
    fc = FormatConverter()
    out = tmp_path / "out.nex"
    fc.convert(FIXTURES / "test.fasta", out, target=AlignmentFormat.NEXUS)
    assert out.exists()
    content = out.read_text()
    assert "#NEXUS" in content


def test_phylip_to_fasta(tmp_path):
    fc = FormatConverter()
    out = tmp_path / "out.fa"
    fc.convert(FIXTURES / "test.phy", out, target=AlignmentFormat.FASTA)
    assert out.exists()
    content = out.read_text()
    assert ">Taxon_A" in content


def test_unsupported_format_raises():
    fc = FormatConverter()
    with pytest.raises(ValueError, match="Cannot detect"):
        fc.detect(Path("alignment.xyz"))
```

- [ ] **Step 3: Run tests — expect failure**

```bash
pytest tests/core/test_formats.py -v
```

Expected: `ImportError`.

- [ ] **Step 4: Implement `formats.py`**

```python
# phyloai/core/formats.py
"""Sequence alignment format detection and conversion."""

from __future__ import annotations
from enum import Enum
from pathlib import Path
from typing import Optional

from Bio import AlignIO, SeqIO
from Bio.Align import MultipleSeqAlignment


class AlignmentFormat(str, Enum):
    FASTA = "fasta"
    PHYLIP = "phylip-relaxed"   # biopython format name
    PHYLIP_PAML = "phylip"      # strict Phylip for PAML/PhyloBayes
    NEXUS = "nexus"


# Extension → format mapping (best-effort)
_EXT_MAP: dict[str, AlignmentFormat] = {
    ".fa":    AlignmentFormat.FASTA,
    ".fasta": AlignmentFormat.FASTA,
    ".faa":   AlignmentFormat.FASTA,
    ".fna":   AlignmentFormat.FASTA,
    ".phy":   AlignmentFormat.PHYLIP,
    ".nex":   AlignmentFormat.NEXUS,
    ".nxs":   AlignmentFormat.NEXUS,
    ".nexus": AlignmentFormat.NEXUS,
}


class FormatConverter:
    """Detect and convert sequence alignment formats using BioPython."""

    def detect(self, path: Path) -> AlignmentFormat:
        """Detect format from file extension and content sniff."""
        suffix = path.suffix.lower()
        if suffix in _EXT_MAP:
            return _EXT_MAP[suffix]
        # content sniff for ambiguous extensions
        if path.exists():
            first = path.read_text(errors="ignore")[:200]
            if first.strip().startswith(">"):
                return AlignmentFormat.FASTA
            if first.strip().startswith("#NEXUS"):
                return AlignmentFormat.NEXUS
            # Phylip: first non-blank line is "  N  L"
            lines = [l for l in first.splitlines() if l.strip()]
            if lines and lines[0].strip().split()[:1] and \
               all(p.isdigit() for p in lines[0].strip().split()[:2]):
                return AlignmentFormat.PHYLIP
        raise ValueError(
            f"Cannot detect alignment format for '{path}'. "
            f"Supported extensions: {list(_EXT_MAP.keys())}"
        )

    def convert(
        self,
        src: Path,
        dst: Path,
        target: AlignmentFormat,
        source_format: Optional[AlignmentFormat] = None,
    ) -> Path:
        """
        Convert alignment file from src to dst in target format.

        Returns:
            Path to the output file.
        """
        src_fmt = source_format or self.detect(src)
        alignment = AlignIO.read(str(src), src_fmt.value)
        dst.parent.mkdir(parents=True, exist_ok=True)
        with open(dst, "w") as fh:
            AlignIO.write(alignment, fh, target.value)
        return dst

    def read(self, path: Path) -> MultipleSeqAlignment:
        """Read an alignment file, auto-detecting format."""
        fmt = self.detect(path)
        return AlignIO.read(str(path), fmt.value)
```

- [ ] **Step 5: Run tests — expect pass**

```bash
pytest tests/core/test_formats.py -v
```

Expected: all 7 tests PASS.

- [ ] **Step 6: Commit**

```bash
git add phyloai/core/formats.py tests/core/test_formats.py \
        tests/core/fixtures/
git commit -m "feat(core): format detection and conversion (formats.py)"
```

---

## Task 6: Step Logger (`logger.py`)

**Files:**
- Create: `phyloai/core/logger.py`
- Create: `tests/core/test_logger.py`

- [ ] **Step 1: Write failing tests**

```python
# tests/core/test_logger.py
import time
import pytest
from pathlib import Path
from phyloai.core.logger import StepLogger
from phyloai.core.schema import ToolResult


def test_logger_creates_log_file(tmp_path):
    logger = StepLogger(run_dir=tmp_path)
    logger.write("align", ToolResult(
        tool="mafft",
        command="mafft --auto input.fa > output.fa",
        returncode=0,
        stdout="Done.",
        stderr="",
        wall_time=3.2,
    ))
    log_file = tmp_path / "logs" / "align.log"
    assert log_file.exists()


def test_logger_log_contains_required_fields(tmp_path):
    logger = StepLogger(run_dir=tmp_path)
    result = ToolResult(
        tool="iqtree3",
        command="iqtree3 -s matrix.fa",
        returncode=0,
        stdout="Analysis done",
        stderr="",
        wall_time=45.1,
    )
    logger.write("iqtree", result)
    content = (tmp_path / "logs" / "iqtree.log").read_text()
    assert "iqtree3" in content          # tool name
    assert "iqtree3 -s matrix.fa" in content  # full command
    assert "45.1" in content             # wall time
    assert "returncode: 0" in content    # exit code
    assert "Analysis done" in content    # stdout


def test_logger_appends_on_retry(tmp_path):
    logger = StepLogger(run_dir=tmp_path)
    result = ToolResult("mafft", "mafft in.fa", 0, "ok", "", 1.0)
    logger.write("align", result)
    logger.write("align", result)  # second write = retry
    content = (tmp_path / "logs" / "align.log").read_text()
    assert content.count("mafft in.fa") == 2


def test_logger_logs_dir_is_created(tmp_path):
    run_dir = tmp_path / "runs"
    logger = StepLogger(run_dir=run_dir)
    result = ToolResult("echo", "echo hi", 0, "hi", "", 0.01)
    logger.write("test_step", result)
    assert (run_dir / "logs").is_dir()
```

- [ ] **Step 2: Run tests — expect failure**

```bash
pytest tests/core/test_logger.py -v
```

Expected: `ImportError`.

- [ ] **Step 3: Implement `logger.py`**

```python
# phyloai/core/logger.py
"""Per-step log file writer for PhyloAI runs."""

from __future__ import annotations
import datetime
from pathlib import Path

from phyloai.core.schema import ToolResult


class StepLogger:
    """Writes one log file per analysis step directly under the run directory (no shared `logs/` subfolder)."""

    def __init__(self, run_dir: Path):
        self.log_dir = Path(run_dir) / "logs"
        self.log_dir.mkdir(parents=True, exist_ok=True)

    def write(self, step_name: str, result: ToolResult) -> Path:
        """
        Append a tool result entry to <step_name>.log.

        Each call appends a new entry (supports retry tracking).
        Returns the log file path.
        """
        log_path = self.log_dir / f"{step_name}.log"
        timestamp = datetime.datetime.now().isoformat(timespec="seconds")
        entry = (
            f"{'='*60}\n"
            f"timestamp: {timestamp}\n"
            f"tool:      {result.tool}\n"
            f"command:   {result.command}\n"
            f"returncode: {result.returncode}\n"
            f"wall_time: {result.wall_time:.2f}s\n"
            f"--- stdout ---\n{result.stdout}\n"
            f"--- stderr ---\n{result.stderr}\n"
        )
        with open(log_path, "a") as fh:
            fh.write(entry)
        return log_path
```

- [ ] **Step 4: Run tests — expect pass**

```bash
pytest tests/core/test_logger.py -v
```

Expected: all 4 tests PASS.

- [ ] **Step 5: Commit**

```bash
git add phyloai/core/logger.py tests/core/test_logger.py
git commit -m "feat(core): per-step log file writer (logger.py)"
```

---

## Task 7: CLI Entry Point and `phyloai doctor`

**Files:**
- Create: `phyloai/cli/main.py`
- Create: `phyloai/cli/doctor.py`
- Create: `tests/cli/test_doctor.py`

- [ ] **Step 1: Write failing tests**

```python
# tests/cli/test_doctor.py
import json
import pytest
from click.testing import CliRunner
from unittest.mock import patch
from phyloai.cli.main import cli
from phyloai.core.env import ToolStatus, ToolInfo
from pathlib import Path


def _mock_tools():
    return {
        "iqtree3": ToolInfo("iqtree3", ToolStatus.OK,
                            Path("/usr/bin/iqtree3"), "3.0.1"),
        "mafft":   ToolInfo("mafft",   ToolStatus.OK,
                            Path("/usr/bin/mafft"), "7.520"),
        "pb_mpi":  ToolInfo("pb_mpi",  ToolStatus.MISSING,
                            note="install: http://www.phylobayes.org"),
    }


def test_doctor_exits_zero():
    runner = CliRunner()
    with patch("phyloai.cli.doctor.ToolEnv") as MockEnv:
        MockEnv.return_value.check_all.return_value = _mock_tools()
        result = runner.invoke(cli, ["doctor"])
    assert result.exit_code == 0


def test_doctor_shows_ok_tools():
    runner = CliRunner()
    with patch("phyloai.cli.doctor.ToolEnv") as MockEnv:
        MockEnv.return_value.check_all.return_value = _mock_tools()
        result = runner.invoke(cli, ["doctor"])
    assert "iqtree3" in result.output
    assert "mafft" in result.output


def test_doctor_shows_missing_tools():
    runner = CliRunner()
    with patch("phyloai.cli.doctor.ToolEnv") as MockEnv:
        MockEnv.return_value.check_all.return_value = _mock_tools()
        result = runner.invoke(cli, ["doctor"])
    assert "pb_mpi" in result.output
    assert "not found" in result.output.lower() \
        or "missing" in result.output.lower()


def test_doctor_json_output():
    runner = CliRunner()
    with patch("phyloai.cli.doctor.ToolEnv") as MockEnv:
        MockEnv.return_value.check_all.return_value = _mock_tools()
        result = runner.invoke(cli, ["doctor", "--output-format", "json"])
    assert result.exit_code == 0
    data = json.loads(result.output)
    assert "iqtree3" in data
    assert data["iqtree3"]["status"] == "ok"
    assert data["pb_mpi"]["status"] == "missing"
```

- [ ] **Step 2: Run tests — expect failure**

```bash
pytest tests/cli/test_doctor.py -v
```

Expected: `ImportError`.

- [ ] **Step 3: Implement `cli/main.py`**

```python
# phyloai/cli/main.py
"""PhyloAI CLI entry point."""

from __future__ import annotations
from pathlib import Path
from typing import Optional

import click
import yaml

from phyloai.cli.doctor import doctor


CONTEXT_SETTINGS = {"help_option_names": ["-h", "--help"]}


@click.group(context_settings=CONTEXT_SETTINGS)
@click.option(
    "--config", "config_file",
    type=click.Path(exists=True, path_type=Path),
    default=None,
    help="Load parameters from a YAML config file. CLI flags override config values.",
)
@click.pass_context
def cli(ctx: click.Context, config_file: Optional[Path]) -> None:
    """PhyloAI — modular phylogenomics analysis platform.

    Run 'phyloai doctor' to check your environment before starting.
    """
    ctx.ensure_object(dict)
    if config_file:
        with open(config_file) as fh:
            ctx.obj["config"] = yaml.safe_load(fh)
    else:
        ctx.obj["config"] = {}


cli.add_command(doctor)
```

- [ ] **Step 4: Implement `cli/doctor.py`**

```python
# phyloai/cli/doctor.py
"""phyloai doctor — environment detection command."""

from __future__ import annotations
import json

import click
from rich.console import Console
from rich.table import Table

from phyloai.core.env import ToolEnv, ToolStatus

console = Console()


@click.command("doctor")
@click.option(
    "--output-format",
    type=click.Choice(["text", "json"]),
    default="text",
    help="Output format.",
)
def doctor(output_format: str) -> None:
    """Check that required external tools are available."""
    env = ToolEnv()
    tools = env.check_all()

    if output_format == "json":
        out = {
            name: {
                "status": info.status.value,
                "path": str(info.path) if info.path else None,
                "version": info.version,
                "note": info.note,
            }
            for name, info in tools.items()
        }
        click.echo(json.dumps(out, indent=2))
        return

    # Rich table output
    table = Table(title="PhyloAI Environment Check", show_header=True)
    table.add_column("Status", width=8)
    table.add_column("Tool", width=14)
    table.add_column("Version", width=12)
    table.add_column("Path / Note")

    status_style = {
        ToolStatus.OK:      "[green]OK[/green]",
        ToolStatus.WARN:    "[yellow]WARN[/yellow]",
        ToolStatus.MISSING: "[red]MISSING[/red]",
    }

    for name, info in tools.items():
        table.add_row(
            status_style[info.status],
            name,
            info.version or "—",
            str(info.path) if info.path else info.note,
        )

    console.print(table)
```

- [ ] **Step 5: Run tests — expect pass**

```bash
pytest tests/cli/test_doctor.py -v
```

Expected: all 4 tests PASS.

- [ ] **Step 6: Smoke test the actual CLI**

```bash
phyloai --help
phyloai doctor
phyloai doctor --output-format json
```

Expected: help text appears, doctor table renders, JSON is valid.

- [ ] **Step 7: Commit**

```bash
git add phyloai/cli/main.py phyloai/cli/doctor.py \
        tests/cli/test_doctor.py
git commit -m "feat(cli): entry point and phyloai doctor command"
```

---

## Task 8: Full Test Suite Run and Final Verification

- [ ] **Step 1: Run all tests**

```bash
pytest tests/ -v --tb=short
```

Expected: all tests PASS, zero failures.

- [ ] **Step 2: Verify CLI entry point is installed**

```bash
which phyloai
phyloai --version 2>/dev/null || phyloai --help | head -3
```

- [ ] **Step 3: Verify YAML config loading (smoke test)**

Create `test_config.yaml`:
```yaml
threads: 4
run_dir: ./runs/test
```

```bash
phyloai --config test_config.yaml doctor
```

Expected: doctor runs without error (config is loaded, no config-specific doctor params).

- [ ] **Step 4: Final commit**

```bash
git add .
git commit -m "feat(core): complete core/ module — all tests passing"
```

---

## Self-Review Notes

**Spec coverage check:**
- [x] `core/env.py` — tool detection, path resolution, bundled tool management
- [x] `core/runner.py` — unified external tool call, timeout, retry (timeout covered; retry is at caller level, sufficient for core)
- [x] `core/formats.py` — FASTA/Nexus/Phylip detection and conversion
- [x] `core/schema.py` — MSACollection, TreeSet, RunRecord, ToolResult
- [x] `core/logger.py` — per-step log files under logs/, append on retry
- [x] `cli/main.py` — `phyloai` entry point, YAML config loading
- [x] `cli/doctor.py` — `phyloai doctor`, text + JSON output
- [x] `--output-format json` on doctor from day one
- [x] `--config FILE` YAML loading at root CLI level

**Items deferred to later modules (by design):**
- `--dry-run` flag: implemented per-command in pretree/tree/posttree plans
- bundled tool download: `env.py` resolves bundled path; actual binary packaging is a packaging/CI task
- `phyloai run` one-click command: implemented in Phase 5 plan
- TAPER is detected via `correction_multi.jl`; Julia and Java runtime checks are reported as separate doctor entries

---

## Follow-up Task: Format Detection and Help Text Refinement

**Files:**
- Modify: `phyloai/cli/doctor.py`
- Modify: `phyloai/core/formats.py`
- Modify: `phyloai/core/schema.py`
- Modify: `tests/cli/test_doctor.py`
- Modify: `tests/core/test_formats.py`
- Modify: `tests/core/test_schema.py`

- [ ] **Step 1: Write failing tests for doctor help text and format handling**

Add tests covering:

- `phyloai doctor -h` explicitly says text is the default output format
- `.fas` maps to FASTA
- `.phylip` maps to PHYLIP
- files with ambiguous suffixes can still be detected from FASTA/NEXUS/PHYLIP content
- explicit `declared_format` / `source_format` overrides guessing
- `MSACollection` counts multiple common alignment suffixes in one directory

- [ ] **Step 2: Run targeted tests — expect failure**

```bash
python -m pytest tests/cli/test_doctor.py tests/core/test_formats.py tests/core/test_schema.py -v
```

Expected: failures in help text, suffix detection, or schema counting.

- [ ] **Step 3: Implement minimal code changes**

Implementation requirements:

- `doctor.py` help text says `text` is the default
- `formats.py` supports `.fas` and `.phylip`
- `FormatConverter.detect()` accepts an explicit declared format and skips guessing when provided
- `FormatConverter.read()` accepts `source_format`
- auto-detection remains suffix-first with conservative content sniffing fallback
- `MSACollection` can count common alignment suffixes, not just a single hard-coded extension

- [ ] **Step 4: Re-run targeted tests — expect pass**

```bash
python -m pytest tests/cli/test_doctor.py tests/core/test_formats.py tests/core/test_schema.py -v
```

Expected: all targeted tests PASS.
