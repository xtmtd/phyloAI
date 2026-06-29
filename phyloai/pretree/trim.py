"""Batch MSA trimming using trimAl, BMGE, or ClipKIT."""

from __future__ import annotations

import shlex
import shutil
import subprocess
import tempfile
import time
from collections.abc import Callable
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Any

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from rich.table import Table

from phyloai.core.checkpoint import Checkpoint
from phyloai.core.schema import COMMON_ALIGNMENT_EXTENSIONS
from phyloai.core.sequence_normalization import detect_seq_type, validate_codon_msa
from phyloai.core.sequence_output_validation import validate_fasta_output


INPUT_EXTENSIONS = {ext for ext in COMMON_ALIGNMENT_EXTENSIONS if ext in {".fa", ".fas", ".fasta", ".faa", ".fna"}}
BMGE_SEQ_TYPE_MAP = {"AA": "AA", "NT": "DNA", "CODON": "CODON"}
CHECKPOINT_FLUSH_INTERVAL = 2.0
TRIMAL_AUTOMATIC_METHODS = {"automated1", "gappyout", "strict", "strictplus"}
TRIMAL_MANUAL_FLAGS = {"-gt", "-st", "-ct", "-cons", "-gapthreshold", "-simthreshold", "-conthreshold"}


def _scan_input(seq_dir: Path) -> tuple[list[Path], list[dict[str, str]]]:
    found: list[Path] = []
    skipped: list[dict[str, str]] = []
    for entry in sorted(seq_dir.iterdir(), key=lambda path: path.name):
        if entry.is_dir():
            skipped.append({"path": str(entry), "reason": "directory"})
        elif not entry.is_file():
            skipped.append({"path": str(entry), "reason": "not a file"})
        elif entry.stat().st_size == 0:
            skipped.append({"path": str(entry), "reason": "empty file"})
        elif entry.suffix.lower() not in INPUT_EXTENSIONS:
            skipped.append({"path": str(entry), "reason": "unrecognized extension"})
        else:
            found.append(entry)
    return found, skipped


def _build_trimal_cmd(input_file: Path, output_file: Path, method: str | None, executable: str = "trimal", backtrans_path: Path | None = None) -> list[str]:
    cmd = [executable, "-in", str(input_file), "-out", str(output_file)]
    if method is not None:
        cmd.append(f"-{method}")
    if backtrans_path is not None:
        cmd.extend(["-backtrans", str(backtrans_path), "-ignorestopcodon"])
    return cmd


def _build_bmge_cmd(input_file: Path, output_file: Path, seq_type: str, matrix: str, entropy: float, java_executable: str, bmge_jar: str) -> list[str]:
    return [
        java_executable, "-jar", bmge_jar,
        "-i", str(input_file),
        "-t", BMGE_SEQ_TYPE_MAP.get(seq_type, "AA"),
        "-m", matrix,
        "-h", str(entropy),
        "-of", str(output_file),
    ]


def _build_clipkit_cmd(input_file: Path, output_file: Path, mode: str, codon: bool, log_path: Path | None, executable: str = "clipkit") -> list[str]:
    cmd = [executable, str(input_file), "-o", str(output_file), "-m", mode]
    if codon:
        cmd.append("--codon")
    if log_path is not None:
        cmd.append("-l")
    return cmd


MANAGED_TOOL_ARGS = {
    "trimal": {"-in", "-out", "-backtrans", "-ignorestopcodon"},
    "bmge": {"-i", "-of", "-t"},
    "clipkit": {"-o", "--output", "-l", "--log", "-co", "--codon", "-t", "--threads"},
}


def _split_tool_args(tool_args: str | None, tool: str) -> list[str]:
    if not tool_args:
        return []
    args = shlex.split(tool_args)
    blocked = MANAGED_TOOL_ARGS.get(tool, set()).intersection(args)
    if blocked:
        blocked_list = ", ".join(sorted(blocked))
        raise ValueError(f"--tool-args cannot include PhyloAI-managed {tool} argument(s): {blocked_list}")
    return args


def _parse_clipkit_log(log_path: Path) -> list[int]:
    kept: list[int] = []
    with open(log_path) as fh:
        for line in fh:
            parts = line.strip().split()
            if len(parts) >= 2 and parts[1] == "keep":
                kept.append(int(parts[0]) - 1)
    return kept


def _project_columns_onto_nt_msa(codon_records: list[SeqRecord], kept_aa_cols: list[int]) -> list[SeqRecord]:
    nt_cols = [col for aa_col in kept_aa_cols for col in (aa_col * 3, aa_col * 3 + 1, aa_col * 3 + 2)]
    projected: list[SeqRecord] = []
    for rec in codon_records:
        seq = str(rec.seq)
        projected.append(SeqRecord(Seq("".join(seq[col] for col in nt_cols if col < len(seq))), id=rec.id, description=rec.description))
    return projected


def _translate_codon_msa(codon_records: list[SeqRecord]) -> list[SeqRecord]:
    translated: list[SeqRecord] = []
    for rec in codon_records:
        seq = str(rec.seq)
        aa_chars: list[str] = []
        for idx in range(0, len(seq) - 2, 3):
            codon = seq[idx:idx + 3]
            if "-" in codon:
                aa_chars.append("-")
            else:
                try:
                    aa = str(Seq(codon).translate())
                except Exception:
                    aa = "X"
                aa_chars.append("-" if aa == "*" else aa)
        translated.append(SeqRecord(Seq("".join(aa_chars)), id=rec.id, description=rec.description))
    return translated


def _read_msa_col_count(path: Path) -> int:
    try:
        for rec in SeqIO.parse(str(path), "fasta"):
            return len(rec.seq)
    except Exception:
        return 0
    return 0


def _make_success_result(msa_path: Path, aa_out: Path, nt_out: Path | None, *, cmd: str, wall_time: float, tool_stderr: str, tool_stdout: str = "", warnings: list[str], length_before: int) -> dict[str, Any]:
    combined = _merge_tool_output(tool_stdout, tool_stderr)
    return {
        "status": "success",
        "input": str(msa_path),
        "output_aa": str(aa_out),
        "output_nt": str(nt_out) if nt_out else None,
        "tool_cmd": cmd,
        "tool_stderr": combined,
        "tool_stdout": tool_stdout,
        "wall_time": wall_time,
        "warnings": warnings,
        "length_before": length_before,
        "length_after": _read_msa_col_count(aa_out) if aa_out.exists() else 0,
    }


def _validate_codon_records(msa_path: Path) -> tuple[list[SeqRecord] | None, list[str], str | None]:
    try:
        records = list(SeqIO.parse(str(msa_path), "fasta"))
    except Exception as exc:
        return None, [], f"could not parse codon MSA: {exc}"
    validation = validate_codon_msa({rec.id: str(rec.seq) for rec in records})
    if validation.skip:
        return None, validation.warnings, "; ".join(validation.warnings)
    return [SeqRecord(Seq(validation.sequences[rec.id]), id=rec.id, description=rec.description) for rec in records], validation.warnings, None


def _write_normalized_codon_msa(records: list[SeqRecord], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    SeqIO.write(records, str(path), "fasta")


def _write_gapless_fasta(input_path: Path, output_path: Path) -> None:
    records = []
    for rec in SeqIO.parse(str(input_path), "fasta"):
        records.append(SeqRecord(Seq(str(rec.seq).replace("-", "")), id=rec.id, description=rec.description))
    output_path.parent.mkdir(parents=True, exist_ok=True)
    SeqIO.write(records, str(output_path), "fasta")


def _normalize_bmge_codon_output(path: Path) -> tuple[list[SeqRecord] | None, list[str], str | None]:
    try:
        records = list(SeqIO.parse(str(path), "fasta"))
    except Exception as exc:
        return None, [], f"could not parse BMGE codon output: {exc}"
    lengths = {len(rec.seq) for rec in records}
    if len(lengths) != 1:
        return None, [], f"BMGE codon output sequences have unequal lengths: {sorted(lengths)}"
    warnings: list[str] = []
    length = next(iter(lengths), 0)
    remainder = length % 3
    if remainder:
        warnings.append(f"BMGE codon output length {length} is not a multiple of 3; removed {remainder} trailing column(s).")
        records = [SeqRecord(Seq(str(rec.seq)[:-remainder]), id=rec.id, description=rec.description) for rec in records]
    validation = validate_codon_msa({rec.id: str(rec.seq) for rec in records})
    if validation.skip:
        return None, validation.warnings + warnings, "; ".join(validation.warnings)
    return [SeqRecord(Seq(validation.sequences[rec.id]), id=rec.id, description=rec.description) for rec in records], warnings + validation.warnings, None


def _run_cmd(cmd: list[str]) -> tuple[subprocess.CompletedProcess[str] | None, float, str | None]:
    start = time.monotonic()
    try:
        proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    except Exception as exc:
        return None, time.monotonic() - start, str(exc)
    return proc, time.monotonic() - start, None


def _merge_tool_output(stdout: str, stderr: str) -> str:
    out = stdout.strip()
    err = stderr.strip()
    if out and err:
        return f"{out}\n{err}"
    return out or err


def _tool_failure_reason(tool_name: str, returncode: int, stderr: str) -> str:
    text = stderr.strip()
    if len(text) <= 600:
        detail = text
    else:
        tail = text.splitlines()[-1] if text.splitlines() else text[-300:]
        detail = f"{text[:300]} ... {tail}"
    return f"{tool_name} exited {returncode}: {detail}"


def _resolve_trimal_method(method: str | None, tool_args_tokens: list[str]) -> str | None:
    if method is None:
        return None
    if method in TRIMAL_AUTOMATIC_METHODS and TRIMAL_MANUAL_FLAGS.intersection(tool_args_tokens):
        return None
    return method


def _trim_one_trimal(msa_path: Path, aa_out_dir: Path, nt_out_dir: Path | None, nt_path: Path | None, method: str | None, seq_type: str, tool_args: str | None, dry_run: bool, executable: str = "trimal") -> dict[str, Any]:
    gene = msa_path.stem
    aa_out = aa_out_dir / f"{gene}.fa"
    nt_out = nt_out_dir / f"{gene}.fa" if nt_out_dir else None
    length_before = _read_msa_col_count(msa_path)
    tool_args_tokens = shlex.split(tool_args) if tool_args else []
    method = _resolve_trimal_method(method, tool_args_tokens)

    if seq_type == "CODON":
        records, warnings, reason = _validate_codon_records(msa_path)
        if reason:
            return {"status": "skipped", "input": str(msa_path), "reason": reason}
        if dry_run:
            method_part = f"-{method}" if method else ""
            return {"status": "dry_run", "input": str(msa_path), "cmd": f"[CODON] trimal -in <temp_aa_msa> -out {nt_out or aa_out} {method_part} -backtrans <temp_cds> -ignorestopcodon", "codon_warnings": warnings}
        with tempfile.TemporaryDirectory(prefix="phyloai_trim_") as tmpdir:
            tmp_aa = Path(tmpdir) / f"{gene}_aa.fa"
            tmp_cds = Path(tmpdir) / f"{gene}_cds.fa"
            SeqIO.write(_translate_codon_msa(records or []), str(tmp_aa), "fasta")
            SeqIO.write([SeqRecord(Seq(str(rec.seq).replace("-", "")), id=rec.id, description=rec.description) for rec in records or []], str(tmp_cds), "fasta")
            target_nt = nt_out or aa_out
            target_nt.parent.mkdir(parents=True, exist_ok=True)
            cmd = _build_trimal_cmd(tmp_aa, target_nt, method, executable, tmp_cds)
            cmd.extend(_split_tool_args(tool_args, "trimal"))
            proc, wall, err = _run_cmd(cmd)
            if err or proc is None or proc.returncode != 0:
                return {"status": "skipped", "input": str(msa_path), "reason": err or _tool_failure_reason("trimal", proc.returncode, proc.stderr), "tool_stderr": "" if proc is None else proc.stderr, "wall_time": wall, "tool_cmd": " ".join(cmd)}
            aa_out.parent.mkdir(parents=True, exist_ok=True)
            SeqIO.write(_translate_codon_msa(list(SeqIO.parse(str(target_nt), "fasta"))), str(aa_out), "fasta")
            return _make_success_result(msa_path, aa_out, nt_out, cmd=" ".join(cmd), wall_time=wall, tool_stdout=proc.stdout, tool_stderr=proc.stderr, warnings=warnings, length_before=length_before)

    if nt_path is not None and seq_type == "AA":
        aa_out.parent.mkdir(parents=True, exist_ok=True)
        cmd_aa = _build_trimal_cmd(msa_path, aa_out, method, executable)
        cmd_aa.extend(_split_tool_args(tool_args, "trimal"))
        if dry_run:
            cmd_nt = _build_trimal_cmd(msa_path, nt_out or aa_out_dir / f"{gene}_nt.fa", method, executable, nt_path)
            cmd_nt.extend(_split_tool_args(tool_args, "trimal"))
            return {"status": "dry_run", "input": str(msa_path), "cmd": " ".join(cmd_aa) + " && " + " ".join(cmd_nt)}
        with tempfile.TemporaryDirectory(prefix="phyloai_trim_") as tmpdir:
            gapless_nt = Path(tmpdir) / f"{gene}_cds.fa"
            _write_gapless_fasta(nt_path, gapless_nt)
            cmd_nt = _build_trimal_cmd(msa_path, nt_out or aa_out_dir / f"{gene}_nt.fa", method, executable, gapless_nt)
            cmd_nt.extend(_split_tool_args(tool_args, "trimal"))
            proc_aa, wall_aa, err_aa = _run_cmd(cmd_aa)
            if err_aa or proc_aa is None or proc_aa.returncode != 0:
                return {"status": "skipped", "input": str(msa_path), "reason": err_aa or _tool_failure_reason("trimal (AA)", proc_aa.returncode, proc_aa.stderr), "tool_stderr": "" if proc_aa is None else proc_aa.stderr, "wall_time": wall_aa, "tool_cmd": " ".join(cmd_aa)}
            if nt_out:
                nt_out.parent.mkdir(parents=True, exist_ok=True)
            proc_nt, wall_nt, err_nt = _run_cmd(cmd_nt)
            wall = wall_aa + wall_nt
            if err_nt or proc_nt is None or proc_nt.returncode != 0:
                return {"status": "skipped", "input": str(msa_path), "reason": err_nt or _tool_failure_reason("trimal (NT backtrans)", proc_nt.returncode, proc_nt.stderr), "tool_stderr": "" if proc_nt is None else proc_nt.stderr, "wall_time": wall, "tool_cmd": " ".join(cmd_nt)}
            return _make_success_result(msa_path, aa_out, nt_out, cmd=" ".join(cmd_nt), wall_time=wall, tool_stdout=proc_nt.stdout, tool_stderr=proc_nt.stderr, warnings=[], length_before=length_before)

    cmd = _build_trimal_cmd(msa_path, aa_out, method, executable)
    cmd.extend(_split_tool_args(tool_args, "trimal"))
    if dry_run:
        return {"status": "dry_run", "input": str(msa_path), "cmd": cmd}
    aa_out.parent.mkdir(parents=True, exist_ok=True)
    proc, wall, err = _run_cmd(cmd)
    if err or proc is None or proc.returncode != 0:
        return {"status": "skipped", "input": str(msa_path), "reason": err or _tool_failure_reason("trimal", proc.returncode, proc.stderr), "tool_stderr": "" if proc is None else proc.stderr, "wall_time": wall, "tool_cmd": " ".join(cmd)}
    return _make_success_result(msa_path, aa_out, None, cmd=" ".join(cmd), wall_time=wall, tool_stdout=proc.stdout, tool_stderr=proc.stderr, warnings=[], length_before=length_before)


def _infer_kept_columns(original_records: list[SeqRecord], trimmed_records: list[SeqRecord]) -> list[int]:
    original_by_id = {rec.id: str(rec.seq) for rec in original_records}
    trimmed_by_id = {rec.id: str(rec.seq) for rec in trimmed_records}
    if set(original_by_id) != set(trimmed_by_id):
        raise ValueError("trimmed output taxon IDs do not match input taxon IDs")
    if not trimmed_records:
        return []
    ids = [rec.id for rec in original_records]
    original_len = len(next(iter(original_by_id.values())))
    trimmed_len = len(next(iter(trimmed_by_id.values())))
    original_cols = [tuple(original_by_id[taxon][idx] for taxon in ids) for idx in range(original_len)]
    trimmed_cols = [tuple(trimmed_by_id[taxon][idx] for taxon in ids) for idx in range(trimmed_len)]
    kept: list[int] = []
    start = 0
    for col in trimmed_cols:
        try:
            idx = original_cols.index(col, start)
        except ValueError as exc:
            raise ValueError("could not map trimmed BMGE AA output back to original AA columns") from exc
        if idx == -1:
            raise ValueError("could not map trimmed BMGE AA output back to original AA columns")
        kept.append(idx)
        start = idx + 1
    return kept


def _trim_one_bmge(msa_path: Path, aa_out_dir: Path, nt_out_dir: Path | None, seq_type: str, matrix: str, entropy: float, tool_args: str | None, dry_run: bool, java_executable: str = "java", bmge_jar: str = "BMGE.jar") -> dict[str, Any]:
    gene = msa_path.stem
    aa_out = aa_out_dir / f"{gene}.fa"
    nt_out = nt_out_dir / f"{gene}.fa" if nt_out_dir else None
    is_dual = nt_out is not None
    effective = "AA" if is_dual else seq_type
    primary = aa_out
    length_before = _read_msa_col_count(msa_path)
    if dry_run:
        cmd = _build_bmge_cmd(msa_path, primary, effective, matrix, entropy, java_executable, bmge_jar)
        cmd.extend(_split_tool_args(tool_args, "bmge"))
        return {"status": "dry_run", "input": str(msa_path), "cmd": cmd}
    warnings: list[str] = []
    input_for_tool = msa_path
    codon_records: list[SeqRecord] | None = None
    if is_dual:
        records, warnings, reason = _validate_codon_records(msa_path)
        if reason:
            return {"status": "skipped", "input": str(msa_path), "reason": reason}
        codon_records = records or []
    primary.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="phyloai_trim_") as tmpdir:
        aa_records_for_projection: list[SeqRecord] | None = None
        if is_dual and codon_records is not None:
            aa_records_for_projection = _translate_codon_msa(codon_records)
            input_for_tool = Path(tmpdir) / f"{gene}_aa.fa"
            SeqIO.write(aa_records_for_projection, str(input_for_tool), "fasta")
        cmd = _build_bmge_cmd(input_for_tool, primary, effective, matrix, entropy, java_executable, bmge_jar)
        cmd.extend(_split_tool_args(tool_args, "bmge"))
        proc, wall, err = _run_cmd(cmd)
        if err or proc is None or proc.returncode != 0:
            return {"status": "skipped", "input": str(msa_path), "reason": err or _tool_failure_reason("BMGE", proc.returncode, proc.stderr), "tool_stderr": "" if proc is None else proc.stderr, "wall_time": wall, "tool_cmd": " ".join(cmd)}
        if not primary.exists():
            return {"status": "skipped", "input": str(msa_path), "reason": f"BMGE exited 0 but did not create output file: {primary}", "tool_stderr": proc.stderr, "wall_time": wall, "tool_cmd": " ".join(cmd)}
        if is_dual and nt_out and codon_records is not None and aa_records_for_projection is not None:
            try:
                kept_cols = _infer_kept_columns(aa_records_for_projection, list(SeqIO.parse(str(aa_out), "fasta")))
            except ValueError as exc:
                return {"status": "skipped", "input": str(msa_path), "reason": str(exc), "tool_stderr": proc.stderr, "wall_time": wall, "tool_cmd": " ".join(cmd)}
            nt_out.parent.mkdir(parents=True, exist_ok=True)
            SeqIO.write(_project_columns_onto_nt_msa(codon_records, kept_cols), str(nt_out), "fasta")
        return _make_success_result(msa_path, aa_out, nt_out, cmd=" ".join(cmd), wall_time=wall, tool_stdout=proc.stdout, tool_stderr=proc.stderr, warnings=warnings, length_before=length_before)


def _trim_one_clipkit(msa_path: Path, aa_out_dir: Path, nt_out_dir: Path | None, nt_path: Path | None, mode: str, seq_type: str, tool_args: str | None, dry_run: bool, executable: str = "clipkit") -> dict[str, Any]:
    gene = msa_path.stem
    aa_out = aa_out_dir / f"{gene}.fa"
    nt_out = nt_out_dir / f"{gene}.fa" if nt_out_dir else None
    length_before = _read_msa_col_count(msa_path)
    if nt_path is not None and seq_type == "AA":
        cmd = _build_clipkit_cmd(msa_path, aa_out, mode, False, Path(str(aa_out) + ".log"), executable)
        cmd.extend(_split_tool_args(tool_args, "clipkit"))
        if dry_run:
            return {"status": "dry_run", "input": str(msa_path), "cmd": cmd}
        aa_out.parent.mkdir(parents=True, exist_ok=True)
        if nt_out:
            nt_out.parent.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(prefix="phyloai_trim_") as tmpdir:
            tmp_aa = Path(tmpdir) / f"{gene}.fa"
            log = Path(str(tmp_aa) + ".log")
            cmd = _build_clipkit_cmd(msa_path, tmp_aa, mode, False, log, executable)
            cmd.extend(_split_tool_args(tool_args, "clipkit"))
            proc, wall, err = _run_cmd(cmd)
            if err or proc is None or proc.returncode != 0:
                return {"status": "skipped", "input": str(msa_path), "reason": err or _tool_failure_reason("clipkit", proc.returncode, proc.stderr), "tool_stderr": "" if proc is None else proc.stderr, "wall_time": wall, "tool_cmd": " ".join(cmd)}
            shutil.copy2(tmp_aa, aa_out)
            SeqIO.write(_project_columns_onto_nt_msa(list(SeqIO.parse(str(nt_path), "fasta")), _parse_clipkit_log(log)), str(nt_out), "fasta")
            return _make_success_result(msa_path, aa_out, nt_out, cmd=" ".join(cmd), wall_time=wall, tool_stdout=proc.stdout, tool_stderr=proc.stderr, warnings=[], length_before=length_before)
    if seq_type == "CODON":
        nt_primary = nt_out or aa_out
        if dry_run:
            cmd = _build_clipkit_cmd(msa_path, nt_primary, mode, True, None, executable)
            cmd.extend(_split_tool_args(tool_args, "clipkit"))
            return {"status": "dry_run", "input": str(msa_path), "cmd": cmd}
        records, warnings, reason = _validate_codon_records(msa_path)
        if reason:
            return {"status": "skipped", "input": str(msa_path), "reason": reason}
        nt_primary.parent.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(prefix="phyloai_trim_") as tmpdir:
            input_for_tool = Path(tmpdir) / f"{gene}_codon.fa"
            _write_normalized_codon_msa(records or [], input_for_tool)
            cmd = _build_clipkit_cmd(input_for_tool, nt_primary, mode, True, None, executable)
            cmd.extend(_split_tool_args(tool_args, "clipkit"))
            proc, wall, err = _run_cmd(cmd)
            if err or proc is None or proc.returncode != 0:
                return {"status": "skipped", "input": str(msa_path), "reason": err or _tool_failure_reason("clipkit", proc.returncode, proc.stderr), "tool_stderr": "" if proc is None else proc.stderr, "wall_time": wall, "tool_cmd": " ".join(cmd)}
            aa_out.parent.mkdir(parents=True, exist_ok=True)
            SeqIO.write(_translate_codon_msa(list(SeqIO.parse(str(nt_primary), "fasta"))), str(aa_out), "fasta")
            return _make_success_result(msa_path, aa_out, nt_out, cmd=" ".join(cmd), wall_time=wall, tool_stdout=proc.stdout, tool_stderr=proc.stderr, warnings=warnings, length_before=length_before)
    cmd = _build_clipkit_cmd(msa_path, aa_out, mode, False, None, executable)
    cmd.extend(_split_tool_args(tool_args, "clipkit"))
    if dry_run:
        return {"status": "dry_run", "input": str(msa_path), "cmd": cmd}
    aa_out.parent.mkdir(parents=True, exist_ok=True)
    proc, wall, err = _run_cmd(cmd)
    if err or proc is None or proc.returncode != 0:
        return {"status": "skipped", "input": str(msa_path), "reason": err or _tool_failure_reason("clipkit", proc.returncode, proc.stderr), "tool_stderr": "" if proc is None else proc.stderr, "wall_time": wall, "tool_cmd": " ".join(cmd)}
    return _make_success_result(msa_path, aa_out, None, cmd=" ".join(cmd), wall_time=wall, tool_stdout=proc.stdout, tool_stderr=proc.stderr, warnings=[], length_before=length_before)


def verify_trim_outputs(aa_path: Path, nt_path: Path | None) -> bool:
    aa_result = validate_fasta_output(aa_path, require_aligned=True)
    if not aa_result.ok or aa_result.length == 0:
        return False
    if nt_path is None:
        return True
    nt_result = validate_fasta_output(nt_path, require_aligned=True)
    return nt_result.ok and nt_result.length > 0


def _detect_seq_type_from_files(files: list[Path], max_files: int = 3) -> str:
    sequences: list[str] = []
    for file in files[:max_files]:
        try:
            for rec in SeqIO.parse(str(file), "fasta"):
                sequences.append(str(rec.seq))
                if len(sequences) >= 10:
                    break
        except Exception:
            continue
        if len(sequences) >= 10:
            break
    return detect_seq_type(sequences) if sequences else "AA"


def _trim_one_worker(args: tuple[Path, Path, Path | None, Path | None, str, str, str, str, str, float, str | None, bool, str, str, str, str]) -> dict[str, Any]:
    msa_path, aa_out_dir, nt_out_dir, nt_path, tool, seq_type, trimal_method, clipkit_method, bmge_matrix, bmge_entropy, tool_args, dry_run, trimal_exe, java_exe, bmge_jar, clipkit_exe = args
    if tool == "trimal":
        return _trim_one_trimal(msa_path, aa_out_dir, nt_out_dir, nt_path, trimal_method, seq_type, tool_args, dry_run, trimal_exe)
    if tool == "bmge":
        return _trim_one_bmge(msa_path, aa_out_dir, nt_out_dir, seq_type, bmge_matrix, bmge_entropy, tool_args, dry_run, java_exe, bmge_jar)
    return _trim_one_clipkit(msa_path, aa_out_dir, nt_out_dir, nt_path, clipkit_method, seq_type, tool_args, dry_run, clipkit_exe)


def _validate_executable_path(path: Path, tool_name: str) -> Path:
    if not path.exists() or path.is_dir():
        raise FileNotFoundError(f"Required tool '{tool_name}' not found at explicit path: {path}")
    return path


def _resolve_trim_tool_paths(tool: str, trimal_path: Path | None, bmge_path: Path | None, clipkit_path: Path | None, dry_run: bool) -> tuple[str, str, str, str]:
    from phyloai.core.env import TOOL_REGISTRY, ToolEnv, ToolStatus

    env = ToolEnv()
    def require_registered(name: str) -> str:
        meta = TOOL_REGISTRY[name]
        info = env._detect_tool(
            name,
            version_flag=meta.get("version_flag", ""),
            version_args=meta.get("version_args"),
            bundled=meta.get("bundled", False),
            bundled_dir=meta.get("bundled_dir"),
            bundled_executable=meta.get("bundled_executable"),
            path_aliases=meta.get("path_aliases"),
        )
        if info.status != ToolStatus.OK or info.path is None:
            raise FileNotFoundError(f"Required tool '{name}' not found. {meta.get('install', '')}")
        return str(info.path)

    if tool == "trimal":
        trimal = str(_validate_executable_path(trimal_path, "trimal")) if trimal_path else ("trimal" if dry_run else require_registered("trimal"))
        return trimal, "java", "BMGE.jar", "clipkit"
    if tool == "bmge":
        bmge = str(_validate_executable_path(bmge_path, "bmge")) if bmge_path else ("BMGE.jar" if dry_run else require_registered("bmge"))
        java = shutil.which("java") or "java"
        if not dry_run and shutil.which("java") is None:
            raise FileNotFoundError("Required tool 'java' not found on PATH")
        return "trimal", java, bmge, "clipkit"
    clipkit = str(_validate_executable_path(clipkit_path, "clipkit")) if clipkit_path else ("clipkit" if dry_run else require_registered("clipkit"))
    return "trimal", "java", "BMGE.jar", clipkit


def _detect_trim_tool_versions(tool: str, trimal_path: Path | None, bmge_path: Path | None, clipkit_path: Path | None) -> dict[str, str]:
    from phyloai.core.env import TOOL_REGISTRY, ToolEnv, ToolStatus

    overrides = {name: path for name, path in {"trimal": trimal_path, "bmge": bmge_path, "clipkit": clipkit_path}.items() if path is not None}
    env = ToolEnv(tool_paths=overrides)

    def detect(name: str) -> str:
        meta = TOOL_REGISTRY[name]
        info = env._detect_tool(
            name,
            version_flag=meta.get("version_flag", ""),
            version_args=meta.get("version_args"),
            bundled=meta.get("bundled", False),
            bundled_dir=meta.get("bundled_dir"),
            bundled_executable=meta.get("bundled_executable"),
            path_aliases=meta.get("path_aliases"),
        )
        return info.version or "unknown" if info.status == ToolStatus.OK else "unknown"

    if tool == "bmge":
        return {"bmge": detect("bmge"), "java": detect("java")}
    return {tool: detect(tool)}


def run_trim(msa_dir: Path, output_dir: Path, tool: str = "trimal", seq_type: str = "auto", nt_dir: Path | None = None, trimal_method: str = "automated1", bmge_matrix: str | None = None, bmge_entropy: float = 0.5, clipkit_method: str = "smart-gap", trimal_path: Path | None = None, bmge_path: Path | None = None, clipkit_path: Path | None = None, threads: int = 4, tool_args: str | None = None, overwrite: bool = False, resume: bool = False, dry_run: bool = False, quiet: bool = False, progress_callback: Callable[[Path], None] | None = None) -> dict[str, Any]:
    from phyloai.core.checkpoint import load_checkpoint, save_checkpoint_atomic, validate_resume_params
    from phyloai.pretree.checkpoint_helpers import build_initial_checkpoint, mark_task, plan_resume

    start = time.monotonic()
    if overwrite and resume:
        raise ValueError("--overwrite and --resume are mutually exclusive.")
    if seq_type == "CODON" and nt_dir is not None:
        raise ValueError("CODON mode does not use --nt-dir. Place codon-aligned MSA in --msa-dir and omit --nt-dir.")
    if threads < 1:
        raise ValueError("--threads must be at least 1.")
    if tool not in {"trimal", "bmge", "clipkit"}:
        raise ValueError(f"Unknown trim tool: {tool}")
    warnings: list[str] = []
    tool_args_tokens = shlex.split(tool_args) if tool_args else []
    effective_trimal_method: str | None = trimal_method
    if tool == "trimal" and trimal_method in TRIMAL_AUTOMATIC_METHODS and TRIMAL_MANUAL_FLAGS.intersection(tool_args_tokens):
        warnings.append(
            f"--tool-args contains manual trimAl thresholds; ignoring --trimal-method {trimal_method}."
        )
        effective_trimal_method = None
    _split_tool_args(tool_args, tool)
    found, scan_skipped = _scan_input(msa_dir)
    if not found and not dry_run:
        raise ValueError("No valid input MSA files found in --msa-dir.")
    resolved_seq_type = _detect_seq_type_from_files(found) if seq_type == "auto" and found else ("AA" if seq_type == "auto" else seq_type)
    if seq_type == "auto":
        warnings.append(f"seq_type auto-detected as '{resolved_seq_type}'.")
    bmge_mode4 = tool == "bmge" and resolved_seq_type == "AA" and nt_dir is not None
    if bmge_mode4:
        warnings.append("BMGE AA+NT mode trims AA MSAs with BMGE -t AA, then projects kept AA columns onto matching codon-aligned NT MSAs.")
    if bmge_matrix is None:
        bmge_matrix = "BLOSUM62" if resolved_seq_type in {"AA", "CODON"} else "DNAPAM100:2"
    is_dual = nt_dir is not None or resolved_seq_type == "CODON"
    aa_out_dir = output_dir / "seqs" / "faa" if is_dual else output_dir / "seqs"
    nt_out_dir = output_dir / "seqs" / "fna" if is_dual else None
    logs_dir = output_dir / "logs"
    trimal_exe, java_exe, bmge_jar, clipkit_exe = _resolve_trim_tool_paths(tool, trimal_path, bmge_path, clipkit_path, dry_run)
    params: dict[str, Any] = {
        "msa_dir": str(msa_dir), "nt_dir": str(nt_dir) if nt_dir else None, "seq_type": resolved_seq_type,
        "tool": tool, "threads": threads, "tool_args": tool_args, "output_dir": str(output_dir),
        "overwrite": overwrite, "dry_run": dry_run, "resume": resume, "quiet": quiet,
        "trimal_method": effective_trimal_method if tool == "trimal" else None,
        "trimal_path": str(trimal_path) if tool == "trimal" and trimal_path else None,
        "bmge_matrix": bmge_matrix if tool == "bmge" else None,
        "bmge_entropy": bmge_entropy if tool == "bmge" else None,
        "bmge_path": str(bmge_path) if tool == "bmge" and bmge_path else None,
        "clipkit_method": clipkit_method if tool == "clipkit" else None,
        "clipkit_path": str(clipkit_path) if tool == "clipkit" and clipkit_path else None,
    }

    checkpoint = None
    ckpt_path = output_dir / "checkpoint.json"
    to_run_ids: list[str] | None = None
    resume_success_results: list[dict[str, Any]] = []
    if resume:
        checkpoint = load_checkpoint(ckpt_path)
        validate_resume_params(checkpoint, params, step="pretree.trim")
        to_run_ids, _ = plan_resume(checkpoint, verify_trim_outputs)
        resume_success_results = _results_from_checkpoint_successes(checkpoint, exclude=set(to_run_ids))
        if not to_run_ids:
            return _build_trim_payload(
                file_results=resume_success_results,
                skipped=list(scan_skipped),
                params=params,
                global_warnings=warnings,
                wall_time=0.0,
                tool_versions=_detect_trim_tool_versions(tool, trimal_path, bmge_path, clipkit_path),
            )
        found = [Path(task.input) for task in checkpoint.tasks if task.task_id in set(to_run_ids)]
    elif not dry_run:
        if output_dir.exists() and any(output_dir.iterdir()):
            if not overwrite:
                raise ValueError(f"Output directory '{output_dir}' already exists and is non-empty. Use --overwrite to replace it or --resume to continue.")
            shutil.rmtree(output_dir)
        aa_out_dir.mkdir(parents=True, exist_ok=True)
        if nt_out_dir:
            nt_out_dir.mkdir(parents=True, exist_ok=True)
        checkpoint = build_initial_checkpoint(
            step="pretree.trim",
            command=_build_trim_command(params),
            params=params,
            inputs=found,
            output_for=lambda path: aa_out_dir / f"{path.stem}.fa",
            nt_output_for=(lambda path: None) if nt_out_dir is None else (lambda path: nt_out_dir / f"{path.stem}.fa"),
        )
        save_checkpoint_atomic(checkpoint, ckpt_path)
        to_run_ids = [path.stem for path in found]

    def nt_match(msa_path: Path) -> Path | None:
        if nt_dir is None:
            return None
        return next((path for path in nt_dir.iterdir() if path.is_file() and path.stem == msa_path.stem), None)

    worker_args = []
    skipped: list[dict[str, str]] = list(scan_skipped)
    for msa_path in found:
        effective_msa = nt_match(msa_path) if bmge_mode4 else msa_path
        nt_path_for_gene = None if bmge_mode4 else nt_match(msa_path)
        if nt_dir is not None and (effective_msa is None or (not bmge_mode4 and nt_path_for_gene is None)):
            skipped.append({"path": str(msa_path), "reason": "nt_pairing_missing"})
            if checkpoint is not None and not dry_run:
                mark_task(checkpoint, msa_path.stem, status="failed", reason="nt_pairing_missing")
            continue
        worker_args.append((effective_msa or msa_path, aa_out_dir, nt_out_dir, nt_path_for_gene, tool, "CODON" if bmge_mode4 else resolved_seq_type, effective_trimal_method, clipkit_method, bmge_matrix, bmge_entropy, tool_args, dry_run, trimal_exe, java_exe, bmge_jar, clipkit_exe))

    file_results: list[dict[str, Any]] = []
    all_results: list[dict[str, Any]] = []
    dry_cmds: list[str] = []
    last_flush = time.monotonic()
    ckpt_write = checkpoint is not None and to_run_ids and not dry_run
    if not dry_run:
        logs_dir.mkdir(parents=True, exist_ok=True)
    to_run_set = set(to_run_ids or [])
    if ckpt_write:
        runnable_ids = {Path(args[0]).stem for args in worker_args}
        for task_id in runnable_ids:
            mark_task(checkpoint, task_id, status="running", reason=None)
        save_checkpoint_atomic(checkpoint, ckpt_path)

    def maybe_flush(force: bool = False) -> None:
        nonlocal last_flush
        if ckpt_write and (force or time.monotonic() - last_flush >= CHECKPOINT_FLUSH_INTERVAL):
            save_checkpoint_atomic(checkpoint, ckpt_path)
            last_flush = time.monotonic()

    interrupted = False
    try:
        with ProcessPoolExecutor(max_workers=threads) as pool:
            futures = {pool.submit(_trim_one_worker, args): args[0] for args in worker_args}
            for future in as_completed(futures):
                gene_path = futures[future]
                try:
                    result = future.result()
                except Exception as exc:
                    result = {"status": "skipped", "input": str(gene_path), "reason": str(exc)}
                all_results.append(result)
                task_id = Path(result.get("input", str(gene_path))).stem
                if result["status"] == "dry_run":
                    cmd = result.get("cmd", "")
                    dry_cmds.append(" ".join(cmd) if isinstance(cmd, list) else str(cmd))
                elif result["status"] == "skipped":
                    skipped.append({"path": result.get("input", str(gene_path)), "reason": result.get("reason", "unknown")})
                    if ckpt_write and task_id in to_run_set:
                        mark_task(checkpoint, task_id, status="failed", reason=result.get("reason"))
                else:
                    aa_out = Path(result["output_aa"])
                    nt_out = Path(result["output_nt"]) if result.get("output_nt") else None
                    if verify_trim_outputs(aa_out, nt_out):
                        file_results.append(result)
                        if not dry_run:
                            (logs_dir / f"{gene_path.stem}.log").write_text(result.get("tool_stderr", ""))
                        if ckpt_write and task_id in to_run_set:
                            mark_task(checkpoint, task_id, status="success", reason=None)
                    else:
                        skipped.append({"path": result["input"], "reason": "output validation failed (empty or unequal lengths)"})
                        if ckpt_write and task_id in to_run_set:
                            mark_task(checkpoint, task_id, status="failed", reason="output validation failed")
                if progress_callback:
                    progress_callback(gene_path)
                maybe_flush()
    except KeyboardInterrupt:
        interrupted = True

    if not dry_run and not file_results and not resume_success_results:
        raise ValueError("No genes were trimmed: all input files failed or were skipped.")

    if ckpt_write:
        if file_results:
            checkpoint.status = "interrupted" if interrupted else "success"
        else:
            checkpoint.status = "error"
        checkpoint.completed_at = None if interrupted else checkpoint.touch()
        save_checkpoint_atomic(checkpoint, ckpt_path, fsync=True)

    if interrupted:
        raise KeyboardInterrupt
    file_results = resume_success_results + file_results
    if resume_success_results and not dry_run:
        log_dir = output_dir / "logs"
        log_dir.mkdir(parents=True, exist_ok=True)
        for res in resume_success_results:
            if not res.get("tool_stderr"):
                locus = Path(res["input"]).stem
                log_path = log_dir / f"{locus}.log"
                if not log_path.exists():
                    log_path.write_text("# resumed from checkpoint — original stderr unavailable\n")
    return _build_trim_payload(file_results=file_results, skipped=skipped, params=params, global_warnings=warnings, wall_time=time.monotonic() - start, tool_versions={} if dry_run else _detect_trim_tool_versions(tool, trimal_path, bmge_path, clipkit_path), dry_run_cmds=dry_cmds if dry_run else None)


def _build_trim_command(params: dict[str, Any]) -> str:
    parts = ["phyloai", "pretree", "trim"]
    parts.extend(["--msa-dir", str(params["msa_dir"])])
    parts.extend(["--output-dir", str(params["output_dir"])])
    parts.extend(["--tool", str(params["tool"])])
    parts.extend(["--seq-type", str(params["seq_type"])])
    parts.extend(["--threads", str(params["threads"])])
    if params.get("nt_dir"):
        parts.extend(["--nt-dir", str(params["nt_dir"])])
    tool = params.get("tool", "")
    if tool == "trimal":
        if params.get("trimal_method"):
            parts.extend(["--trimal-method", str(params["trimal_method"])])
        if params.get("trimal_path"):
            parts.extend(["--trimal-path", str(params["trimal_path"])])
    elif tool == "bmge":
        parts.extend(["--bmge-entropy", str(params["bmge_entropy"])])
        if params.get("bmge_matrix"):
            parts.extend(["--bmge-matrix", str(params["bmge_matrix"])])
        if params.get("bmge_path"):
            parts.extend(["--bmge-path", str(params["bmge_path"])])
    elif tool == "clipkit":
        if params.get("clipkit_method"):
            parts.extend(["--clipkit-method", str(params["clipkit_method"])])
        if params.get("clipkit_path"):
            parts.extend(["--clipkit-path", str(params["clipkit_path"])])
    if params.get("tool_args"):
        parts.extend(["--tool-args", params["tool_args"]])
    if params.get("overwrite"):
        parts.append("--overwrite")
    if params.get("resume"):
        parts.append("--resume")
    if params.get("dry_run"):
        parts.append("--dry-run")
    if params.get("quiet"):
        parts.append("--quiet")
    return " ".join(parts)


def _build_trim_payload(*, file_results: list[dict[str, Any]], skipped: list[dict[str, str]], params: dict[str, Any], global_warnings: list[str], wall_time: float, tool_versions: dict[str, str], dry_run_cmds: list[str] | None = None) -> dict[str, Any]:
    before = [r.get("length_before", 0) for r in file_results if r.get("length_before", 0)]
    after = [r.get("length_after", 0) for r in file_results if r.get("length_after", 0)]

    def stats(values: list[int | float]) -> dict[str, Any]:
        return {"mean": round(sum(values) / len(values), 1), "min": min(values), "max": max(values)} if values else {"mean": 0.0, "min": 0, "max": 0}

    removed = [round((r["length_before"] - r["length_after"]) / r["length_before"] * 100, 1) for r in file_results if r.get("length_before")]
    is_dry = dry_run_cmds is not None
    files = []
    for r in file_results:
        locus = Path(r["input"]).stem
        tool_cmd = r.get("tool_cmd", "")
        entry = {
            "gene": locus,
            "length_before": r.get("length_before", 0),
            "length_after": r.get("length_after", 0),
            "columns_removed": r.get("length_before", 0) - r.get("length_after", 0),
            "outputs": [out for out in [r.get("output_aa"), r.get("output_nt")] if out],
            "log_file": f"logs/{locus}.log",
        }
        if tool_cmd:
            entry["cmd"] = shlex.split(tool_cmd)
            entry["wall_time"] = r.get("wall_time", 0.0)
        files.append(entry)
    return {
        "status": "success" if file_results or is_dry else "error",
        "command": _build_trim_command(params),
        "wall_time": wall_time,
        "tool_versions": tool_versions,
        "params": params,
        "key_results": {"total_genes": len(file_results) + len(skipped), "trimmed_genes": len(file_results), "skipped_genes": len(skipped), "skipped_reasons": _count_reasons(skipped), "length_before": stats(before), "length_after": stats(after), "columns_removed_pct": stats(removed)},
        "error": None if file_results or is_dry else "No genes were trimmed.",
        "data": {"mode": _determine_mode(params), "summary": {"n_input_files": len(file_results) + len(skipped) if not is_dry else len(file_results) + len(skipped) + len(dry_run_cmds), "n_trimmed": len(file_results), "n_skipped": len(skipped)}, "skipped": [{"gene": Path(s["path"]).stem, "reason": s["reason"]} for s in skipped], "warnings": global_warnings, "dry_run_cmds": dry_run_cmds or [], "files": files},
    }


def _results_from_checkpoint_successes(checkpoint: Checkpoint, exclude: set[str]) -> list[dict[str, Any]]:
    results: list[dict[str, Any]] = []
    for task in checkpoint.tasks:
        if task.task_id in exclude or task.status != "success":
            continue
        aa_path = Path(task.outputs.get("aa") or "")
        nt_path = Path(task.outputs.get("nt")) if task.outputs.get("nt") else None
        if not aa_path or not verify_trim_outputs(aa_path, nt_path):
            continue
        results.append(
            {
                "status": "success",
                "input": task.input,
                "output_aa": str(aa_path),
                "output_nt": str(nt_path) if nt_path else None,
                "tool_cmd": "",
                "tool_stderr": "",
                "length_before": _read_msa_col_count(Path(task.input)),
                "length_after": _read_msa_col_count(aa_path),
                "wall_time": 0.0,
                "warnings": [],
            }
        )
    return results


def _count_reasons(skipped: list[dict[str, str]]) -> dict[str, int]:
    counts: dict[str, int] = {}
    for item in skipped:
        reason = item.get("reason", "unknown")
        counts[reason] = counts.get(reason, 0) + 1
    return counts


def _determine_mode(params: dict[str, Any]) -> str:
    if params.get("nt_dir") and params.get("bmge_matrix") is not None:
        return "CODON"
    if params.get("nt_dir"):
        return "AA+NT"
    if params.get("seq_type") == "NT":
        return "NT-only"
    return "AA-only"


def render_trim_summary_table(summary: dict[str, Any]) -> Table:
    table = Table(title="pretree trim summary")
    table.add_column("Metric")
    table.add_column("Value")
    for key in ["n_input_files", "n_trimmed", "n_skipped"]:
        table.add_row(key, str(summary.get(key, "")))
    return table
