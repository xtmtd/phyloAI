"""Multispecies coalescent inference with wASTRAL (ASTER)."""

from __future__ import annotations

import json
import os
import re as _re
import shlex
import shutil
import subprocess
import time as _time
from pathlib import Path
from typing import Any

from phyloai.core.env import ToolEnv


# Gene tree file extensions wastral can read (newick variants)
_WASTRAL_TREE_EXTENSIONS = frozenset({
    ".nwk", ".tre", ".tree", ".nw", ".trees", ".newick", ".treefile",
})


# ---- Two-tier --tool-args management ----------------------------------

# Tier 1: BLOCKED -- phyloAI always manages these; reject hard
_WASTRAL_BLOCKED_FLAGS = frozenset({"-i", "-o"})
_WASTRAL_BLOCKED_IO_CHARS = frozenset({"<", ">", "|"})

# Tier 2: OVERRIDEABLE -- suppress phyloAI's own flag if present in --tool-args
# Maps phyloAI parameter group -> set of wastral flag strings
_WASTRAL_OVERRIDE_MAP: dict[str, frozenset[str]] = {
    "mode": frozenset({"--mode"}),
    "boot": frozenset({"-u", "--support"}),
    "extra_rounds": frozenset({"-R", "-r", "-s"}),
    "tree_boot": frozenset({"--lrt", "--bayes", "--bootstrap", "-x", "-n"}),
    "threads": frozenset({"-t"}),
}


def _check_managed_flag_conflict(tool_args: str) -> None:
    """Reject BLOCKED flags and I/O redirects in --tool-args."""
    tokens = shlex.split(tool_args)
    for token in tokens:
        if token in _WASTRAL_BLOCKED_FLAGS:
            raise ValueError(f"Blocked managed flag in --tool-args: {token}")
        if any(c in token for c in _WASTRAL_BLOCKED_IO_CHARS):
            raise ValueError(f"Blocked I/O override in --tool-args: {token}")


def _is_flag_overridden(group: str, tool_tokens: set[str]) -> bool:
    """Check whether any flag in an override group appears in --tool-args."""
    flags = _WASTRAL_OVERRIDE_MAP.get(group, frozenset())
    return bool(flags & tool_tokens)


# ---- Input scanning ---------------------------------------------------

def _scan_input_wastral(
    tree_dir: Path,
) -> tuple[list[Path], list[dict[str, str]]]:
    """Scan a directory for valid gene tree files.

    Returns:
        (valid_files, skipped_entries)
    """
    if not tree_dir.exists():
        return [], []

    found: list[Path] = []
    skipped: list[dict[str, str]] = []

    for entry in sorted(tree_dir.iterdir()):
        if entry.is_dir():
            skipped.append({"path": str(entry), "reason": "directory"})
            continue
        if not entry.is_file():
            skipped.append({"path": str(entry), "reason": "not a regular file"})
            continue
        if entry.stat().st_size == 0:
            skipped.append({"path": str(entry), "reason": "empty file"})
            continue

        ext = entry.suffix.lower()
        if ext in _WASTRAL_TREE_EXTENSIONS:
            found.append(entry)
        else:
            skipped.append({"path": str(entry), "reason": f"unrecognized extension: {ext}"})

    return found, skipped


def _merge_gene_trees(
    tree_dir: Path,
    output_path: Path,
) -> tuple[int, list[dict[str, str]]]:
    """Scan tree_dir for newick files, merge into one file (one tree per line).

    Returns:
        (count_of_trees_merged, skipped_entries)
    """
    found, skipped = _scan_input_wastral(tree_dir)

    count = 0
    with open(output_path, "w") as out:
        for f in found:
            content = f.read_text().strip()
            if not content:
                continue
            for line in content.splitlines():
                line = line.strip()
                if line:
                    out.write(line + "\n")
                    count += 1

    return count, skipped


# ---- Executable resolution and version detection ----------------------

def _resolve_wastral_path(wastral_path: str | None, dry_run: bool) -> str:
    """Resolve wastral executable path.

    Priority: 1) explicit --wastral-path  2) ToolEnv.require("wastral")
    """
    if wastral_path:
        p = Path(wastral_path)
        if not p.exists():
            raise ValueError(f"--wastral-path does not exist: {wastral_path}")
        if not os.access(str(p), os.X_OK):
            raise ValueError(f"--wastral-path is not executable: {wastral_path}")
        return wastral_path
    if dry_run:
        return "wastral"
    for name in ("wastral", "aster"):
        env = ToolEnv()
        try:
            return str(env.require(name))
        except FileNotFoundError:
            continue
    raise FileNotFoundError(
        "wastral not found. Install from https://github.com/chaoszhang/ASTER "
        "or use --wastral-path."
    )


def _detect_wastral_version(executable: str) -> dict[str, str]:
    """Detect wastral version via -v or -h fallback."""
    combined = ""
    for flag in ("-v", "-h"):
        try:
            proc = subprocess.run(
                [executable, flag],
                capture_output=True, text=True, timeout=10,
            )
            combined = proc.stdout + "\n" + proc.stderr
            if proc.returncode == 0 or (proc.stdout.strip() or proc.stderr.strip()):
                break
        except Exception:
            continue

    if not combined.strip():
        return {"wastral": "unknown"}

    m = _re.search(r"Version:\s*v?([\d.]+)", combined, _re.IGNORECASE)
    if m:
        return {"wastral": m.group(1)}
    m = _re.search(r"version\s+([\d.]+)", combined, _re.IGNORECASE)
    if m:
        return {"wastral": m.group(1)}
    m = _re.search(r"\bv?(\d+(?:\.\d+)+)\b", combined)
    if m:
        return {"wastral": m.group(1)}
    return {"wastral": "unknown"}


# ---- Command builder --------------------------------------------------

def _build_wastral_cmd(
    input_path: Path,
    output_path: Path,
    *,
    mode: int,
    boot: int,
    extra_rounds: bool,
    tree_boot_type: str,
    tree_boot_min: float | None,
    tree_boot_max: float | None,
    threads: int,
    outgroup: str | None = None,
    executable: str = "wastral",
    tool_args: str | None = None,
) -> list[str]:
    """Build the full wastral command line."""
    cmd = [executable]

    # Managed input/output (always emitted, blocked in --tool-args)
    cmd.extend(["-i", str(input_path)])
    cmd.extend(["-o", str(output_path)])

    tool_tokens = set(shlex.split(tool_args)) if tool_args else set()

    if tool_args:
        _check_managed_flag_conflict(tool_args)

    # --mode
    if not _is_flag_overridden("mode", tool_tokens):
        cmd.extend(["--mode", str(mode)])

    # -u (branch support)
    if not _is_flag_overridden("boot", tool_tokens):
        cmd.extend(["-u", str(boot)])

    # -t (threads)
    if not _is_flag_overridden("threads", tool_tokens):
        cmd.extend(["-t", str(threads)])

    # Tree boot type (input gene tree branch support)
    if tree_boot_type != "auto" and not _is_flag_overridden("tree_boot", tool_tokens):
        boot_flag_map = {
            "likelihood": "--lrt",
            "abayes": "--bayes",
            "bootstrap": "--bootstrap",
        }
        boot_defaults = {
            "likelihood":  ("-d", "0",   "-x", "1.0", "-n", "0.0"),
            "abayes":      ("-d", "0.333", "-x", "1.0", "-n", "0.333"),
            "bootstrap":   ("-d", "0",   "-x", "100", "-n", "0"),
        }
        flag = boot_flag_map[tree_boot_type]
        cmd.append(flag)

        # -d, -x, -n: preset defaults, overridable by explicit min/max
        defaults = boot_defaults[tree_boot_type]
        max_val = str(tree_boot_max) if tree_boot_max is not None else defaults[3]
        min_val = str(tree_boot_min) if tree_boot_min is not None else defaults[5]
        cmd.extend([
            defaults[0], defaults[1],  # -d VALUE (hardcoded per type)
            "-x", max_val,
            "-n", min_val,
        ])

    # -R (extra rounds / exhaustive search)
    if extra_rounds and not _is_flag_overridden("extra_rounds", tool_tokens):
        cmd.append("-R")

    # --root (outgroup)
    if outgroup is not None:
        cmd.extend(["--root", outgroup])

    # --tool-args appended last
    if tool_args:
        cmd.extend(shlex.split(tool_args))

    return cmd


# ---- Result assembly --------------------------------------------------

def _has_boot_3(boot: int, tool_args: str | None) -> bool:
    """Check whether wASTRAL is configured to produce freqQuad.csv."""
    if boot == 3:
        return True
    if tool_args:
        tokens = shlex.split(tool_args)
        for i, t in enumerate(tokens):
            if t in ("-u", "--support"):
                if i + 1 < len(tokens) and tokens[i + 1] == "3":
                    return True
    return False


def _assemble_wastral_result(
    *,
    run_start: float,
    wastral_exe: str,
    tree: Path | None,
    tree_dir: Path | None,
    output_dir: Path,
    mode: int,
    boot: int,
    extra_rounds: bool,
    tool_stderr: str = "",
    tree_boot_type: str,
    resolved_boot_type: str | None = None,
    tree_boot_min: float | None,
    tree_boot_max: float | None,
    outgroup: str | None,
    threads: int,
    wastral_path: str | None,
    tool_args: str | None,
    overwrite: bool,
    dry_run: bool,
    quiet: bool = False,
    n_input_trees: int,
    input_path: Path,
    cmd: list[str],
    wall_time: float,
    skipped: list[dict[str, str]],
    warnings_list: list[str],
    is_error: bool,
    error_msg: str | None,
) -> dict[str, Any]:
    """Build the result.json payload."""
    input_mode = "--tree-dir" if tree_dir is not None else "--tree"

    if dry_run:
        versions = {"wastral": "unknown"}
    else:
        try:
            versions = _detect_wastral_version(wastral_exe)
        except Exception:
            versions = {"wastral": "unknown"}

    # Reconstruct CLI invocation string
    cmd_parts = ["phyloai", "tree", "msc"]
    if tree is not None:
        cmd_parts.extend(["--tree", str(tree)])
    else:
        cmd_parts.extend(["--tree-dir", str(tree_dir)])
    cmd_parts.extend([
        "--mode", str(mode),
        "--boot", str(boot),
        "-t", str(threads),
        "-o", str(output_dir),
    ])
    if extra_rounds:
        cmd_parts.append("-R")
    if tree_boot_type != "auto":
        cmd_parts.extend(["--tree-boot-type", tree_boot_type])
        if tree_boot_min is not None:
            cmd_parts.extend(["--tree-boot-min", str(tree_boot_min)])
        if tree_boot_max is not None:
            cmd_parts.extend(["--tree-boot-max", str(tree_boot_max)])
    if wastral_path:
        cmd_parts.extend(["--wastral-path", str(wastral_path)])
    if tool_args:
        if " " in tool_args:
            cmd_parts.append(f"--tool-args '{tool_args}'")
        else:
            cmd_parts.extend(["--tool-args", tool_args])
    if overwrite:
        cmd_parts.append("--overwrite")
    if outgroup:
        cmd_parts.extend(["--outgroup", outgroup])
    if dry_run:
        cmd_parts.append("--dry-run")
    cmd_str = " ".join(cmd_parts)

    input_data: dict[str, Any] = {
        "path": str(input_path),
    }
    if input_mode == "--tree-dir":
        input_data["n_trees"] = n_input_trees

    return {
        "status": "error" if is_error else "success",
        "command": cmd_str,
        "wall_time": _time.monotonic() - run_start,
        "tool_versions": versions,
        "params": {
            "tree": str(tree) if tree else None,
            "tree_dir": str(tree_dir) if tree_dir else None,
            "mode": mode,
            "boot": boot,
            "extra_rounds": extra_rounds,
            "tree_boot_type": tree_boot_type,
            "tree_boot_min": tree_boot_min,
            "tree_boot_max": tree_boot_max,
            "outgroup": outgroup,
            "output_dir": str(output_dir),
            "threads": threads,
            "overwrite": overwrite,
            "tool_args": tool_args,
            "wastral_path": wastral_path,
            "dry_run": dry_run,
            "quiet": quiet,
        },
        "key_results": {
            "mode": mode,
            "boot": boot,
            "extra_rounds": extra_rounds,
            "tree_boot_type": resolved_boot_type if resolved_boot_type else tree_boot_type,
            "outgroup": outgroup,
            "n_input_trees": n_input_trees,
            "input_mode": input_mode,
        },
        "error": error_msg,
        "data": {
            "input_mode": input_mode,
            "input": input_data,
            "output_tree": str(output_dir / "wastral.tre"),
            "cmd": cmd,
            "tool_stderr": tool_stderr,
            "tool_log": "wastral.log",
            "skipped": skipped,
            "warnings": warnings_list,
            **({"freq_quad_csv": str(output_dir / "freqQuad.csv")} if _has_boot_3(boot, tool_args) else {}),
        },
    }


# ---- Main entry point -------------------------------------------------

def run_wastral(
    *,
    tree: Path | None = None,
    tree_dir: Path | None = None,
    output_dir: Path = Path("runs/tree/msc"),
    mode: int = 1,
    boot: int = 1,
    extra_rounds: bool = False,
    tree_boot_type: str = "auto",
    tree_boot_min: float | None = None,
    tree_boot_max: float | None = None,
    outgroup: str | None = None,
    threads: int = 4,
    wastral_path: str | None = None,
    tool_args: str | None = None,
    overwrite: bool = False,
    dry_run: bool = False,
    quiet: bool = False,
) -> dict[str, Any]:
    """Run wASTRAL multispecies coalescent species tree inference.

    Returns a result.json-compatible payload dict.
    """
    # --- Library-level parameter validation ---
    if mode not in (1, 2, 3, 4):
        raise ValueError(f"--mode must be 1, 2, 3, or 4. Got: {mode}")
    if boot not in (0, 1, 2, 3):
        raise ValueError(f"--boot must be 0, 1, 2, or 3. Got: {boot}")
    if threads < 1:
        raise ValueError(f"--threads must be >= 1. Got: {threads}")
    if outgroup is not None and outgroup.strip() == "":
        raise ValueError("--outgroup must not be empty")

    # --- Input mutual exclusivity ---
    if (tree is None and tree_dir is None) or (tree is not None and tree_dir is not None):
        raise ValueError(
            "Either --tree or --tree-dir must be provided (mutually exclusive)."
        )

    # --- Input path validation ---
    if tree is not None and not tree.exists():
        raise ValueError(f"--tree does not exist: {tree}")
    if tree_dir is not None and not tree_dir.exists():
        raise ValueError(f"--tree-dir does not exist: {tree_dir}")

    # --- Tree boot type + min/max validation ---
    if tree_boot_type != "auto" and tree_boot_type not in {"likelihood", "abayes", "bootstrap"}:
        raise ValueError(
            f"Invalid --tree-boot-type: {tree_boot_type}. "
            f"Valid: auto, likelihood, abayes, bootstrap"
        )

    if tree_boot_type == "auto":
        if tree_boot_min is not None or tree_boot_max is not None:
            raise ValueError(
                "--tree-boot-min and --tree-boot-max are only valid when "
                "--tree-boot-type is not 'auto'."
            )

    if tree_boot_min is not None and tree_boot_max is not None and tree_boot_min >= tree_boot_max:
        raise ValueError(
            f"--tree-boot-min ({tree_boot_min}) must be < --tree-boot-max ({tree_boot_max})."
        )

    # --- Output directory conflict ---
    if not dry_run:
        if overwrite and output_dir.exists():
            shutil.rmtree(output_dir)
        if not overwrite and output_dir.exists() and any(output_dir.iterdir()):
            raise ValueError(
                f"Output directory {output_dir} already exists and is non-empty. "
                "Use --overwrite to replace."
            )

    run_start = _time.monotonic()

    # --- Input resolution (before tool resolution to ensure input errors come first) ---
    skipped: list[dict[str, str]] = []
    warnings_list: list[str] = []
    n_input_trees = 0

    if tree is not None:
        # Single file mode
        input_path = tree
        n_input_trees = 0  # unknown, wastral counts internally
    else:
        assert tree_dir is not None
        # Scan for valid files (always; no file writes)
        valid_files, scanned_skipped = _scan_input_wastral(tree_dir)
        skipped = scanned_skipped
        n_valid_files = len(valid_files)

        if n_valid_files == 0:
            raise ValueError(
                f"No valid gene tree files found in --tree-dir: {tree_dir}"
            )

        if n_valid_files == 1 and not quiet:
            warnings_list.append(
                "Exactly 1 valid gene tree file in --tree-dir. "
                "Consider using --tree mode directly."
            )

        # Warn about non-newick files
        unrecognized = [
            s for s in skipped
            if s.get("reason", "").startswith("unrecognized")
        ]
        if unrecognized and not quiet:
            warnings_list.append(
                f"--tree-dir contains {len(unrecognized)} non-newick file(s); "
                "skipped. See result.json data.skipped for details."
            )

        if dry_run:
            # Use a not-yet-existing path for command display only
            merged_path = output_dir / "merged.trees"
            # Count tree lines from valid files (match _merge_gene_trees logic)
            n_input_trees = 0
            for f in valid_files:
                content = f.read_text().strip()
                if content:
                    n_input_trees += len([l for l in content.splitlines() if l.strip()])
            input_path = merged_path
        else:
            # Write merged.trees only during real execution
            output_dir.mkdir(parents=True, exist_ok=True)
            merged_path = output_dir / "merged.trees"
            n_input_trees, _ = _merge_gene_trees(tree_dir, merged_path)
            input_path = merged_path

    # --- Resolve executable (after input validation, skip for dry_run) ---
    wastral_exe = _resolve_wastral_path(wastral_path, dry_run)

    # --- Command building ---
    # Resolve paths to absolute so wastral sees correct paths regardless of cwd
    output_tree = (output_dir / "wastral.tre").resolve()
    cmd = _build_wastral_cmd(
        input_path=input_path.resolve(),
        output_path=output_tree,
        mode=mode,
        boot=boot,
        extra_rounds=extra_rounds,
        tree_boot_type=tree_boot_type,
        tree_boot_min=tree_boot_min,
        tree_boot_max=tree_boot_max,
        threads=threads,
        outgroup=outgroup,
        executable=wastral_exe,
        tool_args=tool_args,
    )

    if dry_run:
        return _assemble_wastral_result(
            run_start=run_start,
            wastral_exe=wastral_exe,
            tree=tree, tree_dir=tree_dir,
            output_dir=output_dir,
            mode=mode, boot=boot, extra_rounds=extra_rounds,
            tree_boot_type=tree_boot_type, tree_boot_min=tree_boot_min,
            tree_boot_max=tree_boot_max,
            outgroup=outgroup,
            threads=threads, wastral_path=wastral_path, tool_args=tool_args,
            overwrite=overwrite,
            dry_run=dry_run,
            quiet=quiet,
            n_input_trees=n_input_trees,
            input_path=input_path, cmd=cmd, wall_time=0.0,
            skipped=skipped, warnings_list=warnings_list,
            is_error=False, error_msg=None,
            tool_stderr="",
        )

    # --- Execution (cwd = output_dir so freqQuad.csv lands in output dir) ---
    output_dir.mkdir(parents=True, exist_ok=True)

    try:
        proc = subprocess.run(
            cmd,
            cwd=str(output_dir),
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
        (output_dir / "wastral.log").write_text(proc.stderr or "")
    except Exception as exc:
        return _assemble_wastral_result(
            run_start=run_start,
            wastral_exe=wastral_exe,
            tree=tree, tree_dir=tree_dir,
            output_dir=output_dir,
            mode=mode, boot=boot, extra_rounds=extra_rounds,
            tree_boot_type=tree_boot_type, tree_boot_min=tree_boot_min,
            tree_boot_max=tree_boot_max,
            outgroup=outgroup,
            threads=threads, wastral_path=wastral_path, tool_args=tool_args,
            overwrite=overwrite,
            dry_run=dry_run,
            quiet=quiet,
            n_input_trees=n_input_trees,
            input_path=input_path, cmd=cmd, wall_time=0.0,
            skipped=skipped, warnings_list=warnings_list,
            is_error=True, error_msg=str(exc),
            tool_stderr="",
        )

    wall_time = _time.monotonic() - run_start

    # Resolve tree_boot_type from wastral stderr when set to "auto".
    # Keep the original user value in params; only key_results gets the
    # resolved value.
    resolved_boot_type = tree_boot_type
    if tree_boot_type == "auto" and proc.stderr:
        import re
        stderr_lower = proc.stderr.lower()
        if "bootstrap-like" in stderr_lower:
            resolved_boot_type = "bootstrap"
        elif "posterior probability" in stderr_lower or "abayes" in stderr_lower:
            resolved_boot_type = "abayes"
        elif "likelihood" in stderr_lower and "support" in stderr_lower:
            resolved_boot_type = "likelihood"

    if proc.returncode != 0:
        return _assemble_wastral_result(
            run_start=run_start,
            wastral_exe=wastral_exe,
            tree=tree, tree_dir=tree_dir,
            output_dir=output_dir,
            mode=mode, boot=boot, extra_rounds=extra_rounds,
            tree_boot_type=tree_boot_type, tree_boot_min=tree_boot_min,
            tree_boot_max=tree_boot_max,
            outgroup=outgroup,
            threads=threads, wastral_path=wastral_path, tool_args=tool_args,
            overwrite=overwrite,
            dry_run=dry_run,
            quiet=quiet,
            n_input_trees=n_input_trees,
            input_path=input_path, cmd=cmd, wall_time=wall_time,
            skipped=skipped, warnings_list=warnings_list,
            is_error=True,
            error_msg=f"wastral exited with code {proc.returncode}: {proc.stderr[:200]}",
            tool_stderr=proc.stderr,
        )

    # Verify output tree was produced
    if not output_tree.exists() or output_tree.stat().st_size == 0:
        return _assemble_wastral_result(
            run_start=run_start,
            wastral_exe=wastral_exe,
            tree=tree, tree_dir=tree_dir,
            output_dir=output_dir,
            mode=mode, boot=boot, extra_rounds=extra_rounds,
            tree_boot_type=tree_boot_type, tree_boot_min=tree_boot_min,
            tree_boot_max=tree_boot_max,
            outgroup=outgroup,
            threads=threads, wastral_path=wastral_path, tool_args=tool_args,
            overwrite=overwrite,
            dry_run=dry_run,
            quiet=quiet,
            n_input_trees=n_input_trees,
            input_path=input_path, cmd=cmd, wall_time=wall_time,
            skipped=skipped, warnings_list=warnings_list,
            is_error=True,
            error_msg="wastral did not produce output tree",
            tool_stderr=proc.stderr,
        )

    payload = _assemble_wastral_result(
        run_start=run_start,
        wastral_exe=wastral_exe,
        tree=tree, tree_dir=tree_dir,
        output_dir=output_dir,
        mode=mode, boot=boot, extra_rounds=extra_rounds,
        tree_boot_type=tree_boot_type, tree_boot_min=tree_boot_min,
        tree_boot_max=tree_boot_max,
        resolved_boot_type=resolved_boot_type,
        outgroup=outgroup,
        threads=threads, wastral_path=wastral_path, tool_args=tool_args,
        overwrite=overwrite,
        dry_run=dry_run,
        quiet=quiet,
        n_input_trees=n_input_trees,
        input_path=input_path, cmd=cmd, wall_time=wall_time,
            skipped=skipped, warnings_list=warnings_list,
            is_error=False, error_msg=None,
            tool_stderr="",
        )

    return payload
