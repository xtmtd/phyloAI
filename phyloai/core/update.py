"""Self-update helpers via GitHub releases."""

from __future__ import annotations

import json
import re
import subprocess
import sys
from dataclasses import dataclass
from urllib.request import Request, urlopen

from phyloai import __version__


_GITHUB_API_TAGS = "https://api.github.com/repos/xtmtd/phyloai/tags"
_GITHUB_INSTALL = "git+https://github.com/xtmtd/phyloai.git"


@dataclass
class _SemVer:
    major: int
    minor: int
    patch: int
    prerelease: tuple[int | str, ...] = ()
    build: str = ""

    @staticmethod
    def parse(version: str) -> _SemVer:
        v = version.strip().lstrip("v")
        build = ""
        if "+" in v:
            v, build = v.split("+", 1)
        pre: tuple[int | str, ...] = ()
        if "-" in v:
            v, pre_str = v.split("-", 1)
            parts: list[int | str] = []
            for p in pre_str.split("."):
                try:
                    parts.append(int(p))
                except ValueError:
                    parts.append(p)
            pre = tuple(parts)
        parts = v.split(".")
        major = int(parts[0]) if len(parts) > 0 else 0
        minor = int(parts[1]) if len(parts) > 1 else 0
        patch = int(parts[2]) if len(parts) > 2 else 0
        return _SemVer(major=major, minor=minor, patch=patch, prerelease=pre, build=build)

    def __lt__(self, other: _SemVer) -> bool:
        if (self.major, self.minor, self.patch) != (other.major, other.minor, other.patch):
            return (self.major, self.minor, self.patch) < (other.major, other.minor, other.patch)
        if not self.prerelease and not other.prerelease:
            return False
        if not self.prerelease:
            return False
        if not other.prerelease:
            return True
        return _compare_prerelease(self.prerelease, other.prerelease)


def _compare_prerelease(
    a: tuple[int | str, ...], b: tuple[int | str, ...]
) -> bool:
    min_len = min(len(a), len(b))
    for i in range(min_len):
        if isinstance(a[i], int) and isinstance(b[i], int):
            if a[i] != b[i]:
                return a[i] < b[i]
        elif isinstance(a[i], int):
            return True
        elif isinstance(b[i], int):
            return False
        else:
            if a[i] != b[i]:
                return a[i] < b[i]
    return len(a) < len(b)


def _fetch_all_tags() -> list[str]:
    results: list[str] = []
    url = _GITHUB_API_TAGS + "?per_page=100"
    while url:
        req = Request(
            url,
            headers={
                "Accept": "application/vnd.github+json",
                "X-GitHub-Api-Version": "2022-11-28",
            },
        )
        try:
            with urlopen(req, timeout=10) as resp:
                data = json.loads(resp.read().decode())
        except Exception:
            return results

        if isinstance(data, list):
            for tag in data:
                name = tag.get("name", "")
                if name:
                    results.append(name)

        url = ""
        link_header = resp.headers.get("Link", "")
        for part in link_header.split(","):
            if 'rel="next"' in part:
                url = part.split(";")[0].strip(" <>")
                break

    return results


def _fetch_latest_version() -> str | None:
    tags = _fetch_all_tags()
    if not tags:
        return None

    best: _SemVer | None = None
    best_name: str | None = None
    for name in tags:
        try:
            sv = _SemVer.parse(name)
        except Exception:
            continue
        if best is None or best < sv:
            best = sv
            best_name = name

    return best_name


def check_update() -> dict[str, str | bool]:
    latest = _fetch_latest_version()
    if latest is None:
        return {"status": "error", "message": "Failed to fetch latest version from GitHub."}

    try:
        current = _SemVer.parse(__version__)
        remote = _SemVer.parse(latest)
    except Exception:
        return {"status": "error", "message": f"Could not parse versions (current={__version__}, remote={latest})."}

    if current < remote:
        return {
            "status": "available",
            "current": __version__,
            "latest": latest,
        }
    return {
        "status": "up_to_date",
        "current": __version__,
        "latest": latest,
    }


def run_update(confirm: bool = False) -> int:
    result = check_update()
    if result["status"] == "error":
        print(f"Error: {result['message']}", file=sys.stderr)
        return 1
    if result["status"] == "up_to_date":
        print(f"PhyloAI is up to date (v{result['current']}).")
        return 0

    print(f"Current version: v{result['current']}")
    print(f"Latest version:  {result['latest']}")

    if not confirm:
        try:
            answer = input("Proceed with update? [y/N] ").strip().lower()
        except (EOFError, KeyboardInterrupt):
            print("\nUpdate cancelled.", file=sys.stderr)
            return 1
        if answer not in ("y", "yes"):
            print("Update cancelled.")
            return 0

    print("Updating phyloai via pip ...")
    ret = subprocess.run(
        [sys.executable, "-m", "pip", "install", "--upgrade", f"{_GITHUB_INSTALL}@{result['latest']}"],
        check=False,
    )
    return ret.returncode
