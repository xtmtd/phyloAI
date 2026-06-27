"""Bundled demo dataset path helpers."""

from __future__ import annotations

from pathlib import Path


def resolve_demo_path(*parts: str) -> Path:
    """Return an absolute path inside the bundled demo dataset."""
    return Path(__file__).resolve().parent.joinpath(*parts)


def resolve_raw_dir() -> Path:
    """Return the end-to-end demo raw sequence directory."""
    return resolve_demo_path("end_to_end", "raw")


def resolve_per_step_dir(step: str) -> Path:
    """Return a per-step demo directory."""
    return resolve_demo_path("per_step", step)
