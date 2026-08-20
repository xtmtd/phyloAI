#!/usr/bin/env python3
"""Fail when existing CLI command docs omit required layout sections."""

from __future__ import annotations

import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
REQUIRED = ("purpose", "usage", "inputs", "outputs", "examples", "warnings", "notes")
REQUIRED_ZH = ("目的", "用法", "输入", "输出", "示例", "警告", "备注")
ALIASES = {
    "warnings": ("warnings", "warning"),
    "备注": ("备注", "注意"),
}
# These are setup/integration guides rather than one-document-per-CLI-command references.
EXCLUDED = {"ai-integration", "completion", "installation"}


def headings(path: Path) -> list[str]:
    result = []
    in_fence = False
    for line in path.read_text().splitlines():
        if line.strip().startswith("```"):
            in_fence = not in_fence
            continue
        if not in_fence:
            match = re.match(r"^#{1,6}\s+(.+)$", line)
            if match:
                result.append(match.group(1).strip().lower())
    return result


def missing_sections(path: Path) -> list[str]:
    found = headings(path)
    required = REQUIRED_ZH if path.name.endswith(".zh.md") else REQUIRED
    return [
        section for section in required
        if not any(alias in heading for alias in ALIASES.get(section, (section,)) for heading in found)
    ]


def test_command_docs_follow_layout() -> None:
    failures = []
    for path in sorted((ROOT / "docs/commands").glob("*.md")):
        stem = path.name.removesuffix(".zh.md").removesuffix(".md")
        if stem in EXCLUDED:
            continue
        missing = missing_sections(path)
        if missing:
            failures.append(f"{path}: missing {', '.join(missing)}")
    assert not failures, "\n".join(failures)


if __name__ == "__main__":
    test_command_docs_follow_layout()
