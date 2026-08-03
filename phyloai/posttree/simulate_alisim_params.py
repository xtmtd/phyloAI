"""Extract AliSim simulation parameters from IQ-TREE report files.

`phyloai posttree simulate alisim params` parses batch ``.iqtree`` reports,
pairs each report with a tree file by logical locus name, and writes a
self-contained ``params.tsv`` table that can be fed back to
``alisim iqtree --model-params`` (or edited / consumed by external scripts).
"""
from __future__ import annotations

import csv
import json
import re
import shlex
import shutil
import time as _time
from pathlib import Path
from typing import Any

PARAM_COLUMNS = (
    "id", "seqtype", "length", "subs_model", "subs_rate", "freq",
    "prop_inv", "rate_heterogeneity", "rate_categories", "rate_param",
    "tree_path",
)

_AA_PI_ORDER = "ARNDCQEGHILKMFPSTWYV"
_IQTREE_SUFFIX = ".iqtree"

_INPUT_DATA_RE = re.compile(r"Input data:")
_ALISIM_HEADER_RE = re.compile(r"To simulate an alignment")
_PI_RE = re.compile(r"^\s*pi\(\s*([A-Za-z])\s*\)\s*=\s*([0-9.]+)", re.MULTILINE)

# Model string components: +F{...}, +I{...}, +G<n>{...} / +R<n>{...}.
_F_RE = re.compile(r"\+F(?:{([^}]*)})?")
_I_RE = re.compile(r"\+I\{([^}]*)\}")
_RATE_RE = re.compile(r"\+([GR])(\d+)\{([^}]*)\}")


def _resolve_seqtype(input_line: str) -> str:
    text = input_line.lower()
    if "amino-acid" in text or "protein" in text:
        return "AA"
    if "nucleotide" in text or "dna" in text:
        return "DNA"
    raise ValueError(f"cannot determine sequence type from report line: {input_line.strip()}")


def _extract_pi_frequencies(report_text: str) -> str:
    by_aa: dict[str, str] = {}
    for aa, value in _PI_RE.findall(report_text):
        by_aa.setdefault(aa.upper(), value)
    values = [by_aa[aa] for aa in _AA_PI_ORDER if aa in by_aa]
    if not values:
        raise ValueError(
            "+F was specified without explicit values and no pi(X) lines "
            "were found in the report"
        )
    return "/".join(values)


def _parse_model_string(model_string: str, report_text: str) -> dict[str, str]:
    """Parse one IQ-TREE ``-m`` string into its AliSim table components."""
    subs_model = re.split(r"[\{+]", model_string, maxsplit=1)[0]

    subs_rate = ""
    freq = ""
    prop_inv = ""
    rate_heterogeneity = ""
    rate_categories = ""
    rate_param = ""

    body = model_string[len(subs_model):]
    if body.startswith("{"):
        end = body.find("}")
        if end != -1:
            subs_rate = body[1:end]
            body = body[end + 1:]
        else:
            raise ValueError(f"unbalanced substitution-rate braces in model: {model_string!r}")

    f_match = _F_RE.search(body)
    if f_match:
        inner = f_match.group(1)
        if inner:
            freq = inner
        else:
            freq = _extract_pi_frequencies(report_text)
        body = body[:f_match.start()] + body[f_match.end():]
    i_match = _I_RE.search(body)
    if i_match:
        prop_inv = i_match.group(1)
        body = body[:i_match.start()] + body[i_match.end():]

    rate_match = _RATE_RE.search(body)
    if rate_match:
        rate_heterogeneity = rate_match.group(1)
        rate_categories = rate_match.group(2)
        rate_param = rate_match.group(3)
        body = body[:rate_match.start()] + body[rate_match.end():]

    residue = body.replace("{", "").replace("}", "").strip()
    if residue and residue != "+":
        raise ValueError(f"unrecognized model components: {residue!r} in {model_string!r}")

    def _slash(value: str) -> str:
        return value.replace(",", "/")

    return {
        "subs_model": subs_model,
        "subs_rate": _slash(subs_rate),
        "freq": _slash(freq),
        "prop_inv": prop_inv,
        "rate_heterogeneity": rate_heterogeneity,
        "rate_categories": rate_categories,
        "rate_param": _slash(rate_param),
    }


def parse_iqtree_report(path: Path) -> dict[str, str]:
    """Return AliSim-compatible values extracted from one IQ-TREE report.

    Raises ValueError when the report cannot be parsed.
    """
    text = path.read_text(encoding="utf-8", errors="replace")

    input_line = next(
        (line for line in text.splitlines() if _INPUT_DATA_RE.search(line)), None
    )
    if input_line is None:
        raise ValueError(f"no 'Input data:' line in {path.name}")
    seqtype = _resolve_seqtype(input_line)

    lines = text.splitlines()
    alisim_index = next(
        (idx for idx, line in enumerate(lines) if _ALISIM_HEADER_RE.search(line)), None
    )
    if alisim_index is None:
        raise ValueError(f"no AliSim command in {path.name}")
    command_line = ""
    for line in lines[alisim_index + 1:]:
        if line.strip():
            command_line = line
            break

    tokens = shlex.split(command_line)
    if not tokens or "--alisim" not in tokens:
        raise ValueError(f"malformed AliSim command line in {path.name}: {command_line!r}")

    model_string: str | None = None
    length: str | None = None
    for idx, token in enumerate(tokens):
        if token == "-m" and idx + 1 < len(tokens):
            model_string = tokens[idx + 1]
        elif token == "--length" and idx + 1 < len(tokens):
            length = tokens[idx + 1]
    if model_string is None:
        raise ValueError(f"no -m model string in AliSim command of {path.name}")
    if length is None:
        raise ValueError(f"no --length in AliSim command of {path.name}")

    return {
        "seqtype": seqtype,
        "length": length,
        **_parse_model_string(model_string, text),
    }


def _tree_locus_candidates(filename: str) -> list[str]:
    candidates: list[str] = []
    one_suffix = Path(filename).stem
    if one_suffix != filename:
        candidates.append(one_suffix)
    two_suffix = Path(one_suffix).stem
    if two_suffix != one_suffix:
        candidates.append(two_suffix)
    if not candidates:
        candidates.append(filename)
    return candidates


def _scan_trees_by_locus(tree_dir: Path) -> dict[str, list[Path]]:
    """Map logical locus names to every candidate tree file under tree_dir."""
    mapping: dict[str, list[Path]] = {}
    if not tree_dir.exists() or not tree_dir.is_dir():
        return mapping
    for entry in sorted(tree_dir.rglob("*")):
        if not entry.is_file() or entry.stat().st_size == 0:
            continue
        seen_in_this_entry: set[str] = set()
        for candidate in _tree_locus_candidates(entry.name):
            if candidate in seen_in_this_entry:
                continue
            seen_in_this_entry.add(candidate)
            mapping.setdefault(candidate, []).append(entry)
    return mapping


def run_alisim_params(
    *,
    iqtree_dir: Path,
    tree_dir: Path,
    output_dir: Path,
    overwrite: bool = False,
    dry_run: bool = False,
    quiet: bool = False,
) -> dict[str, Any]:
    """Match empirical parameters from ``.iqtree`` reports to tree files.

    Returns the standard result.json payload.  Raises ValueError on hard
    validation errors (missing directories, tree-name ambiguity).
    """
    run_start = _time.time()

    errors: list[str] = []
    if not iqtree_dir.exists():
        errors.append(f"--iqtree-dir does not exist: {iqtree_dir}")
    elif not iqtree_dir.is_dir():
        errors.append(f"--iqtree-dir is not a directory: {iqtree_dir}")
    if not tree_dir.exists():
        errors.append(f"--tree-dir does not exist: {tree_dir}")
    elif not tree_dir.is_dir():
        errors.append(f"--tree-dir is not a directory: {tree_dir}")
    if errors:
        raise ValueError("; ".join(errors))

    report_paths = sorted(
        path for path in iqtree_dir.rglob(f"*{_IQTREE_SUFFIX}")
        if path.is_file() and path.stat().st_size > 0
    )
    if not report_paths:
        raise ValueError(
            f"no {_IQTREE_SUFFIX} files found under {iqtree_dir}"
        )

    trees_by_locus = _scan_trees_by_locus(tree_dir)

    rows: list[dict[str, str]] = []
    unmatched: list[dict[str, str]] = []
    seq_type_counts: dict[str, int] = {}

    for report_path in report_paths:
        locus = report_path.name[: -len(_IQTREE_SUFFIX)]
        try:
            parsed = parse_iqtree_report(report_path)
        except ValueError as exc:
            unmatched.append({"id": locus, "reason": str(exc)})
            continue

        tree_candidates = trees_by_locus.get(locus, [])
        if len(tree_candidates) > 1:
            names = ", ".join(sorted(p.name for p in tree_candidates))
            raise ValueError(
                f"ambiguous tree matching for locus {locus!r}: multiple tree files "
                f"share this logical locus name ({names}). "
                "Rename files so each locus has exactly one tree file."
            )
        if not tree_candidates:
            unmatched.append({"id": locus, "reason": "no matching tree file found"})
            continue

        seq_type_counts[parsed["seqtype"]] = seq_type_counts.get(parsed["seqtype"], 0) + 1
        rows.append({
            "id": locus,
            **parsed,
            "tree_path": str(tree_candidates[0].resolve()),
        })

    output_dir = output_dir.resolve()
    if not dry_run:
        if output_dir.exists() and any(output_dir.iterdir()):
            if not overwrite:
                raise ValueError(
                    f"Output directory '{output_dir}' already exists and is non-empty. "
                    "Use --overwrite to replace it."
                )
            shutil.rmtree(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        params_tsv = output_dir / "params.tsv"
        with open(params_tsv, "w", encoding="utf-8", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=PARAM_COLUMNS, delimiter="\t")
            writer.writeheader()
            for row in rows:
                writer.writerow(row)

    n_parsed = len(report_paths)
    n_matched = len(rows)
    n_unmatched = len(unmatched)

    payload: dict[str, Any] = {
        "status": "success",
        "command": "phyloai posttree simulate alisim params "
                   f"--iqtree-dir {iqtree_dir} --tree-dir {tree_dir}",
        "wall_time": round(_time.time() - run_start, 3),
        "tool_versions": {},
        "params": {
            "iqtree_dir": str(iqtree_dir.resolve()),
            "tree_dir": str(tree_dir.resolve()),
            "output_dir": str(output_dir),
            "overwrite": overwrite,
            "dry_run": dry_run,
            "quiet": quiet,
        },
        "key_results": {
            "n_loci_parsed": n_parsed,
            "n_loci_matched": n_matched,
            "n_loci_unmatched": n_unmatched,
            "seq_types": dict(sorted(seq_type_counts.items())),
        },
        "error": None,
        "error_category": None,
        "data": {
            "output_files": {
                "params_tsv": {
                    "path": str(output_dir / "params.tsv"),
                    "description": "Simulation parameters extracted from IQ-TREE reports (TSV)",
                }
            },
            "unmatched": unmatched,
        },
    }

    if not dry_run:
        with open(output_dir / "result.json", "w", encoding="utf-8") as fh:
            json.dump(payload, fh, indent=2)

    if not quiet:
        _print_summary(payload, output_dir, dry_run)

    return payload


def _print_summary(payload: dict[str, Any], output_dir: Path, dry_run: bool) -> None:
    kr = payload["key_results"]
    click_echo = __import__("click").echo
    click_echo("Extracting simulation parameters...")
    click_echo(f"  .iqtree files found: {kr['n_loci_parsed']}")
    click_echo(f"  Matched pairs:       {kr['n_loci_matched']}")
    click_echo(f"  Unmatched (no tree): {kr['n_loci_unmatched']}")
    click_echo("")
    click_echo(
        "  Sequence types: " + " ".join(
            f"{key}={value}" for key, value in kr["seq_types"].items()
        ) or "  Sequence types: (none)"
    )
    click_echo("")
    if dry_run:
        click_echo("Dry run: no files written.")
        return
    click_echo(f"Params written to {output_dir / 'params.tsv'}")
    click_echo(f"Result written to {output_dir / 'result.json'}")
