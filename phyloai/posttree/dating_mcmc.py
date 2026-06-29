"""MCMCtree Bayesian dating run management."""
from __future__ import annotations

import os as _os
import random
import re
import shutil
import subprocess
import threading
import time
from pathlib import Path
from typing import Any

from rich.console import Console
from rich.progress import (
    BarColumn,
    Progress,
    TaskProgressColumn,
    TextColumn,
    TimeElapsedColumn,
)
from rich.live import Live

from phyloai.posttree.dating_hessian import HESSIAN_OUTPUT_FILES, HESSIAN_PREFIX

console = Console()

SEQTYPE_CODE = {"AA": 2, "NT": 0}


def detect_seqtype_from_phylip(phylip_path: Path) -> str:
    text = phylip_path.read_text(errors="ignore").upper()
    lines = [l.strip() for l in text.splitlines() if l.strip()]
    aa_only = set("ARNDCQEGHILKMFPSTWYV")
    nt_chars = set("ACGTURYSWKMBDHVN-?")
    for line in lines[1:]:
        parts = line.split(None, 1)
        seq = parts[1] if len(parts) > 1 else parts[0]
        seq = seq.replace(" ", "").replace("-", "").replace("?", "")
        extra = set(seq) - nt_chars
        if extra & aa_only:
            return "AA"
    return "NT"


def count_ndata_from_phylip(phylip_path: Path) -> int:
    text = phylip_path.read_text(errors="ignore")
    blocks = re.findall(r"^\s+\d+\s+\d+\s*$", text, re.MULTILINE)
    return max(1, len(blocks))


def validate_hessian_dir(hessian_dir: Path) -> list[str]:
    errors: list[str] = []
    for fname in HESSIAN_OUTPUT_FILES:
        p = hessian_dir / fname
        if not p.exists():
            errors.append(f"Missing required file in --hessian-dir: {fname}")
        elif p.stat().st_size == 0:
            errors.append(f"Empty required file in --hessian-dir: {fname}")
    return errors


def generate_mcmctree_ctl(
    *,
    seqtype_code: int,
    ndata: int,
    clock: int,
    burnin: int,
    sampfreq: int,
    nsample: int,
    usedata: int,
    seed: int,
) -> str:
    return f"""\
          seed = {seed}
       seqfile = {HESSIAN_PREFIX}.dummy.phy
      treefile = {HESSIAN_PREFIX}.rooted.nwk
       outfile = mcmctree.out

         ndata = {ndata}
       seqtype = {seqtype_code}  * 0: nucleotides; 1:codons; 2:AAs
       usedata = {usedata}    * 0: no data; 1:seq like; 2:use in.BV; 3: out.BV
         clock = {clock}    * 1: global clock; 2: independent rates; 3: correlated rates
       RootAge =   * safe constraint on root age, used if no fossil for root.

       BDparas = 1 1 0.1 M   * birth, death, sampling
   rgene_gamma = 2 20 1   * gamma prior for overall rates for genes
  sigma2_gamma = 1 10 1    * gamma prior for sigma^2     (for clock=2 or 3)

      finetune = 0: .1  .1  .1  .1 .1 .1  * auto (0 or 1) : times, musigma2, rates, mixing, paras, FossilErr

*** These parameters control the MCMC run
***  Note: Total number of MCMC iterations will be burnin + (sampfreq * nsample)

         print = 1
        burnin = {burnin}
      sampfreq = {sampfreq}
       nsample = {nsample}


*** The following parameters only needed to run MCMCtree with exact likelihood (usedata = 1)
*** no need to change anything for approximate likelihood (usedata = 2)

         model = 0    * 0:JC69, 1:K80, 2:F81, 3:F84, 4:HKY85
         alpha = 0.5    * alpha for gamma rates at sites
         ncatG = 4    * No. categories in discrete gamma

     cleandata = 0  * remove sites with ambiguity data (1:yes, 0:no)?

   kappa_gamma = 6 2      * gamma prior for kappa
   alpha_gamma = 1 1      * gamma prior for alpha

*** Note: Make your window wider (100 columns) before running the program.
"""


def count_mcmc_samples(mcmc_txt: Path) -> int:
    if not mcmc_txt.exists():
        return 0
    try:
        lines = mcmc_txt.read_bytes().splitlines()
        data_lines = [l for l in lines[1:] if l.strip()]
        return len(data_lines)
    except Exception:
        return 0


def _resolve_seqtype_and_ndata(
    hessian_dir: Path,
) -> tuple[str, int, str]:
    """Read seq_type from hessian's result.json (preferred) or detect from
    dummy.phy.  ndata is always counted from the actual dummy.phy data
    blocks — this is the ground truth, especially when IQ-TREE merged
    partitions with --merge --rclusterf.
    """
    dummy_phy = hessian_dir / f"{HESSIAN_PREFIX}.dummy.phy"
    if not dummy_phy.exists():
        raise FileNotFoundError(
            f"{dummy_phy} not found in --hessian-dir. "
            "Re-run `phyloai posttree dating hessian` to regenerate the hessian directory."
        )

    ndata = count_ndata_from_phylip(dummy_phy)

    import json
    result_json = hessian_dir / "result.json"
    if result_json.exists():
        try:
            data = json.loads(result_json.read_text(errors="ignore"))
            params = data.get("params") or {}
            seq_type = params.get("seq_type")
            if seq_type in ("AA", "NT"):
                return seq_type, ndata, "hessian-result.json"
        except Exception:
            pass

    seq_type = detect_seqtype_from_phylip(dummy_phy)
    return seq_type, ndata, "dummy.phy-fallback"


def _read_seed_used(run_dir: Path) -> int | None:
    seed_file = run_dir / "SeedUsed"
    if not seed_file.exists():
        return None
    try:
        return int(seed_file.read_text().strip())
    except Exception:
        return None


def _parse_seed_from_ctl(ctl_text: str) -> int | None:
    """Extract seed value from a mcmctree ctl file, or None if not found."""
    for line in ctl_text.splitlines():
        stripped = line.strip().split("*")[0].strip()
        m = re.match(r"^\s*seed\s*=\s*(-?\d+)\s*$", stripped, re.IGNORECASE)
        if m:
            return int(m.group(1))
    return None


def _detect_mcmctree_version(mcmctree_exe: Path) -> str:
    output = ""
    try:
        result = subprocess.run(
            [str(mcmctree_exe)],
            capture_output=True, text=True, timeout=5,
        )
        output = (result.stdout or "") + (result.stderr or "")
    except subprocess.TimeoutExpired as e:
        stdout = e.stdout or e.output or ""
        stderr = e.stderr or ""
        if isinstance(stdout, bytes):
            stdout = stdout.decode(errors="replace")
        if isinstance(stderr, bytes):
            stderr = stderr.decode(errors="replace")
        output = stdout + stderr
    except Exception:
        return "unknown"
    m = re.search(r"paml version (\d+(?:\.\d+)+)", output)
    return m.group(1) if m else "unknown"


def _detect_mcmctree_version_from_logs(run_dirs: list[Path]) -> str:
    for run_dir in run_dirs:
        for log_path in (run_dir / "mcmctree.log", run_dir / "prior" / "mcmctree.log"):
            if not log_path.exists():
                continue
            m = re.search(
                r"paml version (\d+(?:\.\d+)+)",
                log_path.read_text(errors="ignore"),
            )
            if m:
                return m.group(1)
    return "unknown"


class _SampleCounter:
    def __init__(self) -> None:
        self._inode: int | None = None
        self._offset: int = 0
        self._count: int = 0
        self._header_done: bool = False

    def count(self, path: Path) -> int:
        if not path.exists():
            return self._count
        try:
            st = path.stat()
        except OSError:
            return self._count
        if self._inode is not None and st.st_ino != self._inode:
            self._inode = st.st_ino
            self._offset = 0
            self._count = 0
            self._header_done = False
        with open(path, "rb") as fh:
            if not self._header_done:
                fh.seek(0)
                fh.readline()
                self._offset = fh.tell()
                self._header_done = True
            fh.seek(self._offset)
            new_bytes = fh.read()
            if not new_bytes:
                return self._count
            self._offset = fh.tell()
            for line in new_bytes.splitlines():
                stripped = line.strip()
                if not stripped:
                    continue
                first = stripped.split()[0]
                try:
                    float(first)
                except ValueError:
                    continue
                self._count += 1
        self._inode = st.st_ino
        return self._count


_log_tail_threads: list[Any] = []


def _start_log_tail(log_path: Path, *, prefix: str) -> None:
    import threading

    def _tail() -> None:
        pos = 0
        try:
            while True:
                if not log_path.exists():
                    time.sleep(1)
                    continue
                with open(log_path, "rb") as fh:
                    fh.seek(0, 2)
                    end = fh.tell()
                if end > pos:
                    with open(log_path, "rb") as fh:
                        fh.seek(pos)
                        chunk = fh.read(end - pos)
                    pos = end
                    try:
                        text = chunk.decode("utf-8", errors="replace")
                    except Exception:
                        text = ""
                    for line in text.splitlines():
                        if line.strip():
                            console.print(f"[dim]{prefix}[/dim] {line}")
                time.sleep(2)
        except Exception:
            pass

    t = threading.Thread(target=_tail, daemon=True, name=f"tail:{prefix}")
    t.start()
    _log_tail_threads.append(t)


def _derive_prior_ctl(posterior_ctl: str, *, seed: int) -> str:
    text = posterior_ctl
    usedata_repl = "      usedata = 0    * 0: no data; 1:seq like; 2:use in.BV; 3: out.BV"
    seed_repl = f"          seed = {seed}"

    if re.search(r"^\s*usedata\s*=", text, re.MULTILINE):
        text = re.sub(r"^\s*usedata\s*=.*$", usedata_repl, text, count=1, flags=re.MULTILINE)
    else:
        text = text.rstrip() + "\n\n" + usedata_repl + "\n"

    if re.search(r"^\s*seed\s*=", text, re.MULTILINE):
        text = re.sub(r"^\s*seed\s*=.*$", seed_repl, text, count=1, flags=re.MULTILINE)
    else:
        text = text.rstrip() + "\n" + seed_repl + "\n"

    return text


_RE_SEQFILE_CTL = re.compile(r"^\s*seqfile\s*=\s*(\S+)", re.IGNORECASE | re.MULTILINE)
_RE_TREEFILE_CTL = re.compile(r"^\s*treefile\s*=\s*(\S+)", re.IGNORECASE | re.MULTILINE)


def _stage_ctl_files(run_dir: Path, ctl_text: str, ctl_dir: Path) -> None:
    std_names = set(HESSIAN_OUTPUT_FILES)
    for pattern in (_RE_SEQFILE_CTL, _RE_TREEFILE_CTL):
        m = pattern.search(ctl_text)
        if not m:
            continue
        value = m.group(1).strip()
        if value in std_names:
            continue
        src = Path(value)
        if not src.is_absolute():
            src = (ctl_dir / value).resolve()
        if not src.exists():
            continue
        link = run_dir / value
        if link.exists() or link.is_symlink():
            link.unlink()
        _os.symlink(str(src), str(link))


def _setup_run_dir(
    run_dir: Path,
    hessian_dir: Path,
    ctl_text: str,
    ctl_dir: Path | None = None,
) -> None:
    run_dir.mkdir(parents=True, exist_ok=True)
    hessian_abs = hessian_dir.resolve()
    for fname in HESSIAN_OUTPUT_FILES:
        link = run_dir / fname
        target = hessian_abs / fname
        if link.exists() or link.is_symlink():
            link.unlink()
        _os.symlink(str(target), str(link))
    inbv = run_dir / "in.BV"
    if inbv.exists() or inbv.is_symlink():
        inbv.unlink()
    _os.symlink(str(hessian_abs / f"{HESSIAN_PREFIX}.mcmctree.hessian"), str(inbv))
    if ctl_dir is not None:
        _stage_ctl_files(run_dir, ctl_text, ctl_dir)
    (run_dir / "mcmctree.ctl").write_text(ctl_text)


def run_mcmc(
    *,
    hessian_dir: Path,
    ctl: Path | None = None,
    clock: int = 2,
    burnin: int = 100000,
    sample_freq: int = 10,
    nsamples: int = 10000,
    n_runs: int = 2,
    output_dir: Path,
    mcmctree_path: str | None = None,
    overwrite: bool = False,
    dry_run: bool = False,
    quiet: bool = False,
) -> dict[str, Any]:
    t0 = time.time()

    errors = validate_hessian_dir(hessian_dir)
    if errors:
        return _error_result(errors[0], "input")

    if n_runs < 1:
        return _error_result(f"--runs must be >= 1, got {n_runs}", "input")

    if ctl is not None:
        if not ctl.exists():
            return _error_result(f"--ctl does not exist: {ctl}", "input")
        if not ctl.is_file():
            return _error_result(f"--ctl is not a regular file: {ctl}", "input")
        if clock != 2 or burnin != 100000 or sample_freq != 10 or nsamples != 10000:
            return _error_result(
                "--ctl is mutually exclusive with --clock, --burnin, --sample-freq, and --nsamples. "
                "Remove those flags or drop --ctl to use generated settings.",
                "input",
            )

    if mcmctree_path:
        mcmctree_exe = Path(mcmctree_path)
    else:
        found = shutil.which("mcmctree")
        if found is None:
            return _error_result("mcmctree not found. Install PAML.", "env")
        mcmctree_exe = Path(found)

    mcmctree_version = _detect_mcmctree_version(mcmctree_exe)

    seqtype_str, ndata, src = _resolve_seqtype_and_ndata(hessian_dir)

    seqtype_code = SEQTYPE_CODE[seqtype_str]

    if ctl is not None:
        ctl_text = ctl.read_text(errors="ignore")
        ctl_source = "user-supplied"
    else:
        ctl_text = generate_mcmctree_ctl(
            seqtype_code=seqtype_code,
            ndata=ndata,
            clock=clock,
            burnin=burnin,
            sampfreq=sample_freq,
            nsample=nsamples,
            usedata=2,
            seed=-1,
        )
        ctl_source = "generated"

    if dry_run:
        return {
            "status": "success",
            "command": f"phyloai posttree dating mcmc --hessian-dir {hessian_dir}",
            "wall_time": 0.0,
            "tool_versions": {"mcmctree": mcmctree_version},
            "params": {
                "ctl": str(ctl) if ctl else None,
                "clock": clock, "burnin": burnin,
                "sample_freq": sample_freq, "nsamples": nsamples,
                "n_runs": n_runs, "seqtype": seqtype_str, "ndata": ndata,
            },
            "key_results": {},
            "error": None,
            "data": {"ctl": ctl_text, "ctl_source": ctl_source, "warnings": []},
        }

    if overwrite and output_dir.exists():
        shutil.rmtree(output_dir)
    elif not overwrite and not dry_run:
        if output_dir.exists() and any(output_dir.iterdir()):
            return _error_result(
                f"Output directory exists and is not empty: {output_dir}\n"
                "Use --overwrite to replace.",
                "input",
            )
    output_dir.mkdir(parents=True, exist_ok=True)

    output_dir = output_dir.resolve()
    hessian_dir = hessian_dir.resolve()

    (output_dir / "mcmctree.ctl").write_text(ctl_text)

    run_dirs = [output_dir / f"run{i+1}" for i in range(n_runs)]
    run_seeds: dict[str, int] = {}
    for run_dir in run_dirs:
        run_seeds[run_dir.name] = random.randint(1, 2**31 - 1)

    ctl_dir = ctl.parent.resolve() if ctl is not None else None

    _RE_SEED = re.compile(r"^\s*seed\s*=\s*\S+", re.IGNORECASE | re.MULTILINE)
    for run_dir in run_dirs:
        seed = run_seeds[run_dir.name]
        run_ctl = _RE_SEED.sub(f"seed = {seed}", ctl_text) if _RE_SEED.search(ctl_text) else ctl_text + f"\nseed = {seed}\n"
        _setup_run_dir(run_dir, hessian_dir, run_ctl, ctl_dir=ctl_dir)

        prior_dir = run_dir / "prior"
        prior_dir.mkdir(parents=True, exist_ok=True)
        hessian_abs = hessian_dir.resolve()
        for fname in HESSIAN_OUTPUT_FILES:
            link = prior_dir / fname
            target = hessian_abs / fname
            if link.exists() or link.is_symlink():
                link.unlink()
            _os.symlink(str(target), str(link))
        inbv = prior_dir / "in.BV"
        if inbv.exists() or inbv.is_symlink():
            inbv.unlink()
        _os.symlink(str(hessian_abs / f"{HESSIAN_PREFIX}.mcmctree.hessian"), str(inbv))

        prior_ctl = _derive_prior_ctl(run_ctl, seed=seed)
        if ctl_dir is not None:
            _stage_ctl_files(prior_dir, prior_ctl, ctl_dir)
        (prior_dir / "mcmctree.ctl").write_text(prior_ctl)

    mcmctree_env = {**_os.environ, "OMP_NUM_THREADS": "1"}
    proc_log_handles: dict[str, Any] = {}

    def _launch_run(run_dir: Path, *, phase: str) -> subprocess.Popen:
        log_path = run_dir / "mcmctree.log"
        log_fh = open(log_path, "w")
        run_key = f"{run_dir.name}:{phase}"
        proc_log_handles[run_key] = log_fh
        proc = subprocess.Popen(
            [str(mcmctree_exe), "mcmctree.ctl"],
            cwd=run_dir,
            stdout=log_fh,
            stderr=subprocess.STDOUT,
            env=mcmctree_env,
        )
        if not quiet:
            _start_log_tail(log_path, prefix=run_key)
        return proc

    procs: dict[str, subprocess.Popen] = {}
    prior_procs: dict[str, subprocess.Popen] = {}
    try:
        for run_dir in run_dirs:
            procs[run_dir.name] = _launch_run(run_dir, phase="posterior")
            prior_procs[run_dir.name] = _launch_run(run_dir / "prior", phase="prior")

        if quiet:
            while True:
                post_done = all(p.poll() is not None for p in procs.values())
                prior_done = all(p.poll() is not None for p in prior_procs.values())
                if post_done and prior_done:
                    break
                time.sleep(5)
        else:
            progress = Progress(
                TextColumn(" {task.description:<30}"),
                BarColumn(),
                TaskProgressColumn(),
                TextColumn(" {task.completed}/{task.total} samples"),
                TimeElapsedColumn(),
                console=console,
            )
            task_ids: dict[str, Any] = {}
            for run_dir in run_dirs:
                tid = progress.add_task(f"{run_dir.name}-posterior", total=nsamples)
                task_ids[f"{run_dir.name}-posterior"] = tid
                tid = progress.add_task(f"{run_dir.name}-prior", total=nsamples)
                task_ids[f"{run_dir.name}-prior"] = tid

            sample_counters: dict[Path, "_SampleCounter"] = {
                run_dir / "mcmc.txt": _SampleCounter() for run_dir in run_dirs
            }
            sample_counters.update({
                run_dir / "prior" / "mcmc.txt": _SampleCounter() for run_dir in run_dirs
            })

            with Live(progress, console=console, refresh_per_second=0.5):
                while True:
                    for run_dir in run_dirs:
                        n = sample_counters[run_dir / "mcmc.txt"].count(run_dir / "mcmc.txt")
                        progress.update(task_ids[f"{run_dir.name}-posterior"], completed=n)
                        prior_n = sample_counters[run_dir / "prior" / "mcmc.txt"].count(run_dir / "prior" / "mcmc.txt")
                        progress.update(task_ids[f"{run_dir.name}-prior"], completed=prior_n)

                    post_done = all(p.poll() is not None for p in procs.values())
                    prior_done = all(p.poll() is not None for p in prior_procs.values())
                    if post_done and prior_done:
                        break

                    time.sleep(5)

    except KeyboardInterrupt:
        for p in procs.values():
            try:
                p.kill()
            except Exception:
                pass
        for p in prior_procs.values():
            try:
                p.kill()
            except Exception:
                pass
        raise
    finally:
        for fh in proc_log_handles.values():
            try:
                fh.close()
            except Exception:
                pass

    run_failures: list[str] = []
    for run_dir in run_dirs:
        run_key = f"{run_dir.name}:posterior"
        rc = procs[run_dir.name].returncode
        if rc != 0:
            run_failures.append(
                f"{run_dir.name}-posterior: mcmctree exited with code {rc}"
            )
        for required in ("mcmctree.out", "mcmc.txt"):
            p = run_dir / required
            if not p.exists():
                run_failures.append(f"{run_key}: missing {required}")
            elif p.stat().st_size == 0:
                run_failures.append(f"{run_key}: empty {required}")
    for run_dir in run_dirs:
        run_name = run_dir.name
        prior_dir = run_dir / "prior"
        rc = prior_procs[run_name].returncode
        if rc != 0:
            run_failures.append(
                f"{run_name}-prior: mcmctree exited with code {rc}"
            )
        for required in ("mcmctree.out", "mcmc.txt"):
            p = prior_dir / required
            if not p.exists():
                run_failures.append(f"{run_name}-prior: missing {required}")
            elif p.stat().st_size == 0:
                run_failures.append(f"{run_name}-prior: empty {required}")
    all_warnings = run_failures

    from phyloai.posttree.dating_diagnostics import generate_all_diagnostics
    diag_dir = output_dir / "diagnostics"
    diag_summary = generate_all_diagnostics(
        run_dirs=run_dirs,
        diag_dir=diag_dir,
        n_runs=n_runs,
    )

    if mcmctree_version == "unknown":
        mcmctree_version = _detect_mcmctree_version_from_logs(run_dirs)

    if n_runs < 2:
        diag_summary["convergence"] = {"status": "skipped", "reason": "n_runs=1"}

    wall = time.time() - t0
    posterior_failed = any(
        procs[rd.name].returncode != 0 for rd in run_dirs
    ) or any(
        f"-posterior:" in w for w in run_failures
    )
    return {
        "status": "error" if posterior_failed else "success",
        "command": f"phyloai posttree dating mcmc --hessian-dir {hessian_dir}",
        "wall_time": wall,
        "tool_versions": {"mcmctree": mcmctree_version},
        "params": {
            "ctl": str(ctl) if ctl else None,
            "clock": clock, "burnin": burnin,
            "sample_freq": sample_freq, "nsamples": nsamples,
            "n_runs": n_runs, "seqtype": seqtype_str, "ndata": ndata,
            "seqtype_ndata_source": src,
            "ctl_source": ctl_source,
        },
        "key_results": {
            "n_runs": n_runs,
            "n_posterior_failures": sum(
                1 for w in run_failures if "-posterior:" in w
            ),
            "convergence_rho_posterior": next(
                (r["rho"] for r in diag_summary.get("spearman", [])
                 if r["comparison"].startswith("convergence_posterior_")), None
            ) if n_runs >= 2 else None,
        },
        "error": run_failures[0] if posterior_failed else None,
        "error_category": "tool" if posterior_failed else None,
        "data": {
            "diagnostics": diag_summary,
            "output_files": diag_summary.pop("output_files", {}),
            "warnings": all_warnings,
            "return_codes": {
                **{f"{rd.name}:posterior": procs[rd.name].returncode
                   for rd in run_dirs},
                **{f"{rd.name}:prior": prior_procs[rd.name].returncode
                   for rd in run_dirs},
            },
        },
    }


def _error_result(msg: str, category: str) -> dict[str, Any]:
    return {
        "status": "error",
        "command": "",
        "wall_time": 0.0,
        "tool_versions": {},
        "params": {},
        "key_results": {},
        "error": msg,
        "error_category": category,
        "data": {"warnings": [msg]},
    }
