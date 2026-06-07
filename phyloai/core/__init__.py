"""Core infrastructure for PhyloAI."""

from phyloai.core.schema import MSACollection, TreeSet, RunRecord, ToolResult

try:
    from phyloai.core.env import ToolEnv
except ImportError:
    ToolEnv = None  # type: ignore[assignment,misc]

try:
    from phyloai.core.runner import Runner
except ImportError:
    Runner = None  # type: ignore[assignment,misc]

try:
    from phyloai.core.formats import FormatConverter
except ImportError:
    FormatConverter = None  # type: ignore[assignment,misc]

try:
    from phyloai.core.logger import StepLogger
except ImportError:
    StepLogger = None  # type: ignore[assignment,misc]

__all__ = [
    "MSACollection", "TreeSet", "RunRecord", "ToolResult",
    "ToolEnv", "Runner", "FormatConverter", "StepLogger",
]
