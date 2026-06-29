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
    "wastral":    {"required": False, "version_args": [["-v"], ["-h"]],
                    "path_aliases": ["aster"],
                    "install": "https://github.com/chaoszhang/ASTER"},
    "pb_mpi":     {"required": False, "version_flag": "",
                   "install": "https://github.com/bayesiancook/pbmpi"},
    "bpcomp":     {"required": False, "version_flag": "",
                   "install": "https://github.com/bayesiancook/pbmpi"},
    "tracecomp":  {"required": False, "version_flag": "",
                   "install": "https://github.com/bayesiancook/pbmpi"},
    "readpb_mpi": {"required": False, "version_flag": "",
                   "install": "https://github.com/bayesiancook/pbmpi"},
    "mpirun":     {"required": False, "version_flag": "--version",
                   "install": "https://www.open-mpi.org  (or: brew install open-mpi / apt install openmpi-bin)"},
    "mcmctree":   {"required": False, "version_args": [],
                    "version_pattern": r"paml version (\d+(?:\.\d+)+)",
                    "install": "https://github.com/abacus-gene/paml/releases"},
    "correction_multi.jl": {"required": False, "version_args": [["-h"]],
                   "version_pattern": r"Version\s+(\d+(?:\.\d+)+)",
                   "install": "Install TAPER and make correction_multi.jl visible on PATH, or pass --taper-path"},
    "run_treeshrink.py": {"required": False, "version_flag": "--version",
                   "install": "https://github.com/uym2/TreeShrink"},
    "magus":      {"required": False, "version_flag": "--version",
                   "install": "pip install magus-msa"},
    "clipkit":    {"required": False, "version_flag": "--version",
                   "install": "pip install clipkit"},

    "java":       {"required": False, "version_args": [["-version"]],
                   "install": "https://www.java.com/download/"},
    "julia":      {"required": False, "version_args": [["--version"]],
                   "install": "https://julialang.org/downloads/"},
    "trimal":     {"required": True,  "version_flag": "--version",
                    "install": "https://github.com/inab/trimal/releases"},
    "bmge":       {"required": False, "version_args": [["-?"]],
                     "path_aliases": ["BMGE.jar"],
                     "install": "Install BMGE and make BMGE.jar visible on PATH, or pass --bmge-path"},
    "FastTree":   {"required": False, "version_args": [[]],
                     "path_aliases": ["fasttree"],
                     "install": "Download from http://www.microbesonline.org/fasttree/"},
}


TOOL_GROUPS: list[tuple[str, list[str]]] = [
    ("Core", ["mafft", "trimal", "iqtree3"]),
    ("Pre-tree", ["magus", "clipkit", "bmge", "java", "julia", "correction_multi.jl", "run_treeshrink.py"]),
    ("Tree", ["FastTree", "wastral"]),
    ("PhyloBayes MPI", ["pb_mpi", "bpcomp", "tracecomp", "readpb_mpi", "mpirun"]),
    ("Post-tree", ["mcmctree"]),
]


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

    def _version_from_tool_dir(self, path: Path) -> Optional[str]:
        resolved = path.resolve() if path.exists() else path
        tool_dir = resolved.parent
        search_dirs = [tool_dir, tool_dir.parent]
        for sd in search_dirs:
            version = self._normalize_version(sd.name)
            if version:
                return version
        patterns = ["*Manual*.pdf", "*README*", "VERSION", "CHANGELOG", "*.pdf"]
        for sd in list(dict.fromkeys(search_dirs)):
            if not sd.is_dir():
                continue
            for pat in patterns:
                for fpath in sorted(sd.glob(pat), key=lambda p: str(p).lower()):
                    if not fpath.is_file():
                        continue
                    for re_pat in [r"[Vv]ersion\s*(\d+\.\d+(?:\.\d+)?)", r"(\d+\.\d+)"]:
                        m = re.search(re_pat, fpath.name)
                        if m:
                            return m.group(1)
        return None

    def _get_version(self, path: Path, version_args: str | list[list[str]],
                      version_pattern: Optional[str] = None) -> Optional[str]:
        if version_args is None:
            return None
        if isinstance(version_args, str):
            if version_args == "":
                return self._version_from_tool_dir(path)
            candidates = [[version_args]]
        else:
            candidates = version_args if version_args else [[]]
        try:
            for args in candidates:
                command = [str(path), *args]
                if path.suffix.lower() == ".jar":
                    command = ["java", "-jar", str(path), *args]
                elif path.suffix.lower() == ".jl":
                    command = [str(self._tool_paths.get("julia", "julia")), str(path), *args]
                result = subprocess.run(
                    command,
                    capture_output=True, text=True, timeout=5
                )
                output = "\n".join(part.strip() for part in [result.stdout, result.stderr] if part and part.strip())
                if version_pattern:
                    m = re.search(version_pattern, output)
                    if m:
                        return m.group(1)
                for line in output.splitlines():
                    if line.strip():
                        version = self._normalize_version(line.strip()[:200])
                        if version:
                            return version
        except Exception:
            pass
        return self._version_from_tool_dir(path)

    def _detect_tool(self, name: str, version_flag: str = "",
                      version_args: Optional[list[list[str]]] = None,
                      version_pattern: Optional[str] = None,
                      bundled: bool = False,
                      bundled_dir: Optional[str] = None,
                      bundled_executable: Optional[str] = None,
                      path_aliases: Optional[list[str]] = None) -> ToolInfo:
        version_probe = version_args if version_args is not None else version_flag
        override_path = self._tool_paths.get(name)
        if override_path is not None:
            if override_path.exists():
                ver = self._get_version(override_path, version_probe, version_pattern)
                return ToolInfo(name=name, status=ToolStatus.OK,
                                path=override_path, version=ver)
            return ToolInfo(name=name)
        if bundled:
            bundled_name = bundled_executable or name
            bundled_path = self._bundled_dir / (bundled_dir or name) / bundled_name
            if bundled_path.exists():
                ver = self._get_version(bundled_path, version_probe, version_pattern) or self._version_from_bundled_dir(bundled_dir)
                return ToolInfo(name=name, status=ToolStatus.OK,
                                path=bundled_path, version=ver, note="bundled")
        candidates = [name, *(path_aliases or [])]
        for candidate in candidates:
            found = shutil.which(candidate)
            if found:
                p = Path(found)
                ver = self._get_version(p, version_probe, version_pattern)
                return ToolInfo(name=name, status=ToolStatus.OK, path=p, version=ver)
        return ToolInfo(name=name)

    def _version_from_bundled_dir(self, bundled_dir: Optional[str]) -> Optional[str]:
        if not bundled_dir:
            return None
        return self._normalize_version(bundled_dir)

    def check_all(self) -> dict[str, ToolInfo]:
        for name, meta in TOOL_REGISTRY.items():
            info = self._detect_tool(
                name,
                version_flag=meta.get("version_flag", ""),
                version_args=meta.get("version_args"),
                version_pattern=meta.get("version_pattern"),
                bundled=meta.get("bundled", False),
                bundled_dir=meta.get("bundled_dir"),
                bundled_executable=meta.get("bundled_executable"),
                path_aliases=meta.get("path_aliases"),
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
