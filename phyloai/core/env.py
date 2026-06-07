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
    "iqtree2":    {"required": True,  "version_flag": "--version",
                   "install": "https://github.com/iqtree/iqtree2/releases"},
    "mafft":      {"required": True,  "version_flag": "--version",
                   "install": "https://mafft.cbrc.jp/alignment/software/"},
    "astral":     {"required": False, "version_flag": "-v",
                   "install": "https://github.com/chaoszhang/ASTER"},
    "pb_mpi":     {"required": False, "version_flag": "",
                   "install": "http://www.phylobayes.org"},
    "mcmctree":   {"required": False, "version_flag": "",
                   "install": "http://abacus.gene.ucl.ac.uk/software/paml.html"},
    "simphy":     {"required": False, "version_flag": "-v",
                   "install": "https://github.com/adamallo/SimPhy"},
    "treeshrink": {"required": False, "version_flag": "--version",
                   "install": "pip install treeshrink"},
    "magus":      {"required": False, "version_flag": "--version",
                   "install": "pip install magus"},
    "clipkit":    {"required": False, "version_flag": "--version",
                   "install": "pip install clipkit"},
    "phykit":     {"required": False, "version_flag": "version",
                   "install": "pip install phykit"},
    "trimal":     {"required": True,  "version_flag": "--version",
                   "bundled": True,
                   "install": "auto-bundled with phyloai"},
    "bmge":       {"required": False, "version_flag": "",
                   "bundled": True,
                   "install": "auto-bundled with phyloai"},
}


class ToolEnv:
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
            for line in output.splitlines():
                if line.strip():
                    return line.strip()[:80]
        except Exception:
            pass
        return None

    def _detect_tool(self, name: str, version_flag: str = "",
                     bundled: bool = False) -> ToolInfo:
        if bundled:
            bundled_path = self._bundled_dir / name / name
            if bundled_path.exists():
                ver = self._get_version(bundled_path, version_flag)
                return ToolInfo(name=name, status=ToolStatus.OK,
                                path=bundled_path, version=ver, note="bundled")
        found = shutil.which(name)
        if found:
            p = Path(found)
            ver = self._get_version(p, version_flag)
            return ToolInfo(name=name, status=ToolStatus.OK, path=p, version=ver)
        return ToolInfo(name=name)

    def check_all(self) -> dict[str, ToolInfo]:
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
