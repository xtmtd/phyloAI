"""Tool detection and path resolution for PhyloAI."""

from __future__ import annotations
import re
import shutil
import subprocess
from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from typing import Mapping, Optional


class ToolStatus(str, Enum):
    OK = "ok"
    WARN = "warn"
    MISSING = "missing"


@dataclass
class ToolInfo:
    name: str
    status: ToolStatus = ToolStatus.MISSING
    path: Optional[Path] = None
    version: Optional[str] = None
    note: str = ""


TOOL_REGISTRY: dict[str, dict] = {
    "iqtree3":    {"required": True,  "version_flag": "--version",
                   "install": "https://github.com/iqtree/iqtree3/releases"},
    "mafft":      {"required": True,  "version_flag": "--version",
                   "install": "https://mafft.cbrc.jp/alignment/software/"},
    "astral-hybrid": {"required": False, "version_args": [["-v"], ["-h"]],
                   "install": "https://github.com/chaoszhang/ASTER"},
    "pb_mpi":     {"required": False, "version_flag": "",
                   "install": "https://github.com/bayesiancook/pbmpi"},
    "mcmctree":   {"required": False, "version_flag": "",
                   "install": "https://github.com/abacus-gene/paml/releases"},
    "correction_multi.jl": {"required": False, "version_flag": "",
                   "install": "https://github.com/chaoszhang/TAPER"},
    "run_treeshrink.py": {"required": False, "version_flag": "--version",
                   "install": "https://github.com/uym2/TreeShrink"},
    "magus":      {"required": False, "version_flag": "--version",
                   "install": "pip install magus-msa"},
    "clipkit":    {"required": False, "version_flag": "--version",
                   "install": "pip install clipkit"},
    "phykit":     {"required": False, "version_args": [["-h"], ["version"]],
                   "install": "pip install phykit"},
    "java":       {"required": False, "version_args": [["-version"]],
                   "install": "https://www.java.com/download/"},
    "julia":      {"required": False, "version_args": [["--version"]],
                   "install": "https://julialang.org/downloads/"},
    "trimal":     {"required": True,  "version_flag": "--version",
                   "bundled": True,
                   "install": "https://github.com/inab/trimal/releases"},
    "bmge":       {"required": False, "version_flag": "",
                   "bundled": True,
                   "install": "conda install bmge"},
}


class ToolEnv:
    def __init__(self, tool_paths: Optional[Mapping[str, Path | str]] = None):
        self._tools: dict[str, ToolInfo] = {}
        self._bundled_dir = Path(__file__).parent.parent / "bundled"
        self._tool_paths = {
            name: Path(path) for name, path in (tool_paths or {}).items()
        }

    def _normalize_version(self, output: str) -> Optional[str]:
        match = re.search(r"\b(?:v)?(\d+(?:\.\d+)+(?:[-+._a-zA-Z0-9]*)?)\b", output)
        if not match:
            return None
        return match.group(1)

    def _get_version(self, path: Path, version_args: str | list[list[str]]) -> Optional[str]:
        if not version_args:
            return None
        candidates = [[version_args]] if isinstance(version_args, str) else version_args
        try:
            for args in candidates:
                result = subprocess.run(
                    [str(path), *args],
                    capture_output=True, text=True, timeout=5
                )
                output = result.stdout.strip() or result.stderr.strip()
                for line in output.splitlines():
                    if line.strip():
                        version = self._normalize_version(line.strip()[:200])
                        if version:
                            return version
        except Exception:
            pass
        return None

    def _detect_tool(self, name: str, version_flag: str = "",
                     version_args: Optional[list[list[str]]] = None,
                     bundled: bool = False) -> ToolInfo:
        version_probe = version_args if version_args is not None else version_flag
        override_path = self._tool_paths.get(name)
        if override_path is not None:
            if override_path.exists():
                ver = self._get_version(override_path, version_probe)
                return ToolInfo(name=name, status=ToolStatus.OK,
                                path=override_path, version=ver)
            return ToolInfo(name=name)
        if bundled:
            bundled_path = self._bundled_dir / name / name
            if bundled_path.exists():
                ver = self._get_version(bundled_path, version_probe)
                return ToolInfo(name=name, status=ToolStatus.OK,
                                path=bundled_path, version=ver, note="bundled")
        found = shutil.which(name)
        if found:
            p = Path(found)
            ver = self._get_version(p, version_probe)
            return ToolInfo(name=name, status=ToolStatus.OK, path=p, version=ver)
        return ToolInfo(name=name)

    def check_all(self) -> dict[str, ToolInfo]:
        for name, meta in TOOL_REGISTRY.items():
            info = self._detect_tool(
                name,
                version_flag=meta.get("version_flag", ""),
                version_args=meta.get("version_args"),
                bundled=meta.get("bundled", False),
            )
            if info.status == ToolStatus.MISSING:
                info.note = f"install: {meta.get('install', '')}"
            self._tools[name] = info
        return self._tools

    def require(self, name: str) -> Path:
        info = self._tools.get(name) or self._detect_tool(name)
        if info.status != ToolStatus.OK or info.path is None:
            raise FileNotFoundError(
                f"Required tool '{name}' not found. "
                f"{TOOL_REGISTRY.get(name, {}).get('install', '')}"
            )
        return info.path

    def get(self, name: str) -> Optional[Path]:
        try:
            return self.require(name)
        except FileNotFoundError:
            return None
