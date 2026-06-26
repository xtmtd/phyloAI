"""Tests for dating_mcmc pure helpers."""
from __future__ import annotations
from pathlib import Path
import pytest
from phyloai.posttree.dating_mcmc import (
    detect_seqtype_from_phylip,
    count_ndata_from_phylip,
    generate_mcmctree_ctl,
    validate_hessian_dir,
    count_mcmc_samples,
    _resolve_seqtype_and_ndata,
    _derive_prior_ctl,
    _detect_mcmctree_version,
)


DUMMY_PHY_AA = " 3 50\nsp1  MKTV\nsp2  MLTV\nsp3  MSTV\n"
DUMMY_PHY_NT = " 3 50\nsp1  ACGT\nsp2  ACGT\nsp3  ACGT\n"
DUMMY_PHY_2PART = " 3 50\nsp1  MKTV\nsp2  MLTV\nsp3  MSTV\n\n 3 50\nsp1  MKTV\nsp2  MLTV\nsp3  MSTV\n"


def test_detect_seqtype_aa(tmp_path):
    p = tmp_path / "dummy.phy"
    p.write_text(DUMMY_PHY_AA)
    assert detect_seqtype_from_phylip(p) == "AA"


def test_detect_seqtype_nt(tmp_path):
    p = tmp_path / "dummy.phy"
    p.write_text(DUMMY_PHY_NT)
    assert detect_seqtype_from_phylip(p) == "NT"


def test_count_ndata_single(tmp_path):
    p = tmp_path / "dummy.phy"
    p.write_text(DUMMY_PHY_AA)
    assert count_ndata_from_phylip(p) == 1


def test_count_ndata_two_partitions(tmp_path):
    p = tmp_path / "dummy.phy"
    p.write_text(DUMMY_PHY_2PART)
    assert count_ndata_from_phylip(p) == 2


def test_generate_mcmctree_ctl_posterior():
    ctl = generate_mcmctree_ctl(
        seqtype_code=2,
        ndata=1,
        clock=2,
        burnin=100000,
        sampfreq=10,
        nsample=10000,
        usedata=2,
        seed=-1,
    )
    assert "seqtype = 2" in ctl
    assert "usedata = 2" in ctl
    assert "seed = -1" in ctl
    assert "burnin = 100000" in ctl
    assert "ndata = 1" in ctl


def test_generate_mcmctree_ctl_prior():
    ctl = generate_mcmctree_ctl(
        seqtype_code=2,
        ndata=1,
        clock=2,
        burnin=100000,
        sampfreq=10,
        nsample=10000,
        usedata=0,
        seed=42,
    )
    assert "usedata = 0" in ctl
    assert "seed = 42" in ctl


def test_validate_hessian_dir_ok(tmp_path):
    for f in ("iqtree.dummy.phy", "iqtree.rooted.nwk", "iqtree.mcmctree.hessian"):
        (tmp_path / f).write_text("x")
    errs = validate_hessian_dir(tmp_path)
    assert errs == []


def test_validate_hessian_dir_missing(tmp_path):
    (tmp_path / "iqtree.dummy.phy").write_text("x")
    errs = validate_hessian_dir(tmp_path)
    assert len(errs) > 0


def test_count_mcmc_samples_empty(tmp_path):
    p = tmp_path / "mcmc.txt"
    assert count_mcmc_samples(p) == 0


def test_count_mcmc_samples_with_header(tmp_path):
    p = tmp_path / "mcmc.txt"
    p.write_text("Gen\tt_n7\tmu\n1\t0.4\t0.01\n2\t0.41\t0.011\n")
    assert count_mcmc_samples(p) == 2


def test_resolve_seqtype_and_ndata_prefers_hessian_result_json(tmp_path):
    """seq_type from result.json; ndata always from dummy.phy (ground truth)."""
    import json
    (tmp_path / "result.json").write_text(json.dumps({
        "params": {"seq_type": "AA", "n_partitions": 3},
        "key_results": {},
    }))
    # DUMMY_PHY_NT has 1 data block — ndata=1, regardless of result.json
    (tmp_path / "iqtree.dummy.phy").write_text(DUMMY_PHY_NT)
    seq, n, src = _resolve_seqtype_and_ndata(tmp_path)
    assert seq == "AA"       # from result.json (preferred)
    assert n == 1            # from dummy.phy (ground truth, 1 block)
    assert src == "hessian-result.json"


def test_resolve_seqtype_and_ndata_fallback_to_dummy_phy(tmp_path):
    (tmp_path / "iqtree.dummy.phy").write_text(DUMMY_PHY_AA)
    seq, n, src = _resolve_seqtype_and_ndata(tmp_path)
    assert seq == "AA"
    assert n == 1
    assert src == "dummy.phy-fallback"


def test_resolve_seqtype_and_ndata_missing_both_raises(tmp_path):
    with pytest.raises(FileNotFoundError):
        _resolve_seqtype_and_ndata(tmp_path)


def test_resolve_seqtype_and_ndata_partial_result_json_falls_back(tmp_path):
    import json
    (tmp_path / "result.json").write_text(json.dumps({"params": {}}))
    (tmp_path / "iqtree.dummy.phy").write_text(DUMMY_PHY_NT)
    seq, n, src = _resolve_seqtype_and_ndata(tmp_path)
    assert seq == "NT"
    assert n == 1
    assert src == "dummy.phy-fallback"


def test_derive_prior_ctl_substitutes_usedata_and_seed():
    posterior_ctl = (
        "      seed = -1\n"
        "      usedata = 2    * 0: no data; 1:seq like; 2:use in.BV; 3: out.BV\n"
        "      model = 0    * 0:JC69, 1:K80, 2:F81, 3:F84, 4:HKY85\n"
    )
    prior = _derive_prior_ctl(posterior_ctl, seed=12345)
    assert "usedata = 0" in prior
    assert "usedata = 2" not in prior
    assert "seed = 12345" in prior
    assert "seed = -1" not in prior
    assert "model = 0" in prior
    assert "* 0:JC69, 1:K80" in prior


def test_derive_prior_ctl_appends_missing_lines():
    ctl_text = "      ndata = 1\n      clock = 2\n"
    prior = _derive_prior_ctl(ctl_text, seed=42)
    assert "usedata = 0" in prior
    assert "seed = 42" in prior
    assert "ndata = 1" in prior
    assert "clock = 2" in prior


def test_detect_mcmctree_version_reads_timeout_output(tmp_path):
    import subprocess
    from unittest.mock import patch

    err = subprocess.TimeoutExpired(
        cmd=[str(tmp_path / "mcmctree")],
        timeout=5,
        output=b"MCMCTREE in paml version 4.10.10, 27 Jan 2026\n",
    )
    with patch("subprocess.run", side_effect=err):
        assert _detect_mcmctree_version(tmp_path / "mcmctree") == "4.10.10"


def test_detect_mcmctree_version_ignores_nonzero_returncode(tmp_path):
    import subprocess
    from unittest.mock import patch

    result = subprocess.CompletedProcess(
        [str(tmp_path / "mcmctree")],
        returncode=255,
        stdout=(
            "MCMCTREE in paml version 4.10.10, 27 Jan 2026\n\n"
            "error when opening file mcmctree.ctl\n"
        ),
        stderr="",
    )
    with patch("subprocess.run", return_value=result):
        assert _detect_mcmctree_version(tmp_path / "mcmctree") == "4.10.10"


def test_run_mcmc_with_fake_mcmctree_exits_successfully(tmp_path):
    import json, stat as _stat

    hessian = tmp_path / "hessian"
    hessian.mkdir()
    for fname in ("iqtree.dummy.phy", "iqtree.rooted.nwk", "iqtree.mcmctree.hessian"):
        (hessian / fname).write_text("x")
    (hessian / "result.json").write_text(json.dumps({
        "params": {"seq_type": "AA", "n_partitions": 1},
        "key_results": {},
    }))

    fake_mcmctree = tmp_path / "fake_mcmctree"
    fake_mcmctree.write_text(
        r"""#!/usr/bin/env python3
import shutil, sys, time, os
ctl_file = sys.argv[1]
run_dir = os.path.dirname(os.path.abspath(ctl_file))

with open(ctl_file) as f:
    ctl_text = f.read()
import re
match = re.search(r'nsample\s*=\s*(\d+)', ctl_text)
nsample = int(match.group(1)) if match else 2

with open(os.path.join(run_dir, "mcmc.txt"), "w") as f:
    f.write("Gen\tt_n7\tmu\tsigma2\tlnL\n")
    for i in range(1, nsample + 1):
        f.write(f"{i}\t0.42\t0.01\t0.005\t-10.0\n")

with open(os.path.join(run_dir, "mcmctree.out"), "w") as f:
    f.write("Posterior means and 95% Equal-tail CIs\n")
    f.write("t_n7  0.4213 (0.3521, 0.4891)\n")
    f.write("\n")
    f.write("Species tree for FigTree.  Branch lengths = posterior mean times; 95% CIs = labels\n")
    f.write("(sp1,sp2) 7 ;\n")
    f.write("(sp1,sp2) 7 ;\n")
    f.write("(sp1,sp2);\n")
"""
    )
    fake_mcmctree.chmod(fake_mcmctree.stat().st_mode | _stat.S_IXUSR | _stat.S_IXGRP | _stat.S_IXOTH)

    from phyloai.posttree.dating_mcmc import run_mcmc
    output = tmp_path / "out"
    payload = run_mcmc(
        hessian_dir=hessian,
        mcmctree_path=str(fake_mcmctree),
        nsamples=3,
        burnin=0,
        sample_freq=1,
        n_runs=2,
        output_dir=output,
        quiet=True,
    )
    assert payload["status"] == "success", f"Expected success, got error: {payload.get('error')}"
    assert payload["key_results"]["n_runs"] == 2
    assert payload["key_results"]["n_posterior_failures"] == 0
    assert (output / "diagnostics" / "convergence" / "convergence_posterior_run1_vs_run2.pdf").exists()
    assert (output / "diagnostics" / "spearman_correlations.csv").exists()


def test_run_mcmc_backfills_version_from_log_when_probe_unknown(tmp_path, monkeypatch):
    import json, stat as _stat

    hessian = tmp_path / "hessian"
    hessian.mkdir()
    for fname in ("iqtree.dummy.phy", "iqtree.rooted.nwk", "iqtree.mcmctree.hessian"):
        (hessian / fname).write_text("x")
    (hessian / "result.json").write_text(json.dumps({"params": {"seq_type": "AA"}}))

    fake_mcmctree = tmp_path / "fake_mcmctree"
    fake_mcmctree.write_text(
        r'''#!/usr/bin/env python3
import os, re, sys
print("MCMCTREE in paml version 4.10.10, 27 Jan 2026")
ctl_file = sys.argv[1]
run_dir = os.getcwd()
with open(ctl_file) as f:
    ctl_text = f.read()
match = re.search(r'nsample\s*=\s*(\d+)', ctl_text)
nsample = int(match.group(1)) if match else 2
with open(os.path.join(run_dir, "mcmc.txt"), "w") as f:
    f.write("Gen\tt_n7\tt_n8\tt_n9\tmu\n")
    for i in range(1, nsample + 1):
        f.write(f"{i}\t0.42\t0.52\t0.62\t0.01\n")
with open(os.path.join(run_dir, "mcmctree.out"), "w") as f:
    f.write("Posterior means and 95% Equal-tail CIs\n")
    f.write("t_n7  0.4213 (0.3521, 0.4891)\n")
    f.write("t_n8  0.5213 (0.4521, 0.5891)\n")
    f.write("t_n9  0.6213 (0.5521, 0.6891)\n")
    f.write("Species tree for FigTree.\n(sp1,sp2) 7 ;\n")
'''
    )
    fake_mcmctree.chmod(fake_mcmctree.stat().st_mode | _stat.S_IXUSR | _stat.S_IXGRP | _stat.S_IXOTH)

    from phyloai.posttree import dating_mcmc
    monkeypatch.setattr(dating_mcmc, "_detect_mcmctree_version", lambda _: "unknown")

    payload = dating_mcmc.run_mcmc(
        hessian_dir=hessian,
        mcmctree_path=str(fake_mcmctree),
        nsamples=3,
        burnin=0,
        sample_freq=1,
        n_runs=1,
        output_dir=tmp_path / "out",
        quiet=True,
    )
    assert payload["tool_versions"]["mcmctree"] == "4.10.10"


def test_run_mcmc_fake_mcmctree_posterior_failure_is_error(tmp_path):
    import json, stat as _stat

    hessian = tmp_path / "hessian"
    hessian.mkdir()
    for fname in ("iqtree.dummy.phy", "iqtree.rooted.nwk", "iqtree.mcmctree.hessian"):
        (hessian / fname).write_text("x")
    (hessian / "result.json").write_text(json.dumps({
        "params": {"seq_type": "AA", "n_partitions": 1},
    }))

    fake_mcmctree = tmp_path / "fake_fail"
    fake_mcmctree.write_text(
        "#!/usr/bin/env python3\nimport sys\nsys.exit(1)\n"
    )
    fake_mcmctree.chmod(fake_mcmctree.stat().st_mode | _stat.S_IXUSR | _stat.S_IXGRP | _stat.S_IXOTH)

    from phyloai.posttree.dating_mcmc import run_mcmc
    output = tmp_path / "out"
    payload = run_mcmc(
        hessian_dir=hessian,
        mcmctree_path=str(fake_mcmctree),
        nsamples=3, burnin=0, sample_freq=1, n_runs=2,
        output_dir=output, quiet=True,
    )
    assert payload["status"] == "error"


def test_sample_counter_handles_growing_file(tmp_path):
    from phyloai.posttree.dating_mcmc import _SampleCounter
    f = tmp_path / "mcmc.txt"
    f.write_text("Gen\tt_n7\tmu\n")
    c = _SampleCounter()
    assert c.count(f) == 0
    with open(f, "a") as fh:
        fh.write("1\t0.4\t0.01\n2\t0.41\t0.011\n")
    assert c.count(f) == 2
    with open(f, "a") as fh:
        fh.write("3\t0.42\t0.012\n")
    assert c.count(f) == 3


def test_sample_counter_skips_partial_trailing_line(tmp_path):
    from phyloai.posttree.dating_mcmc import _SampleCounter
    f = tmp_path / "mcmc.txt"
    f.write_text("Gen\tt_n7\tmu\n1\t0.4\t0.01\nincomplete")
    c = _SampleCounter()
    assert c.count(f) == 1


def test_sample_counter_recovers_from_file_replacement(tmp_path):
    from phyloai.posttree.dating_mcmc import _SampleCounter
    f = tmp_path / "mcmc.txt"
    f.write_text("Gen\tt_n7\tmu\n1\t0.4\t0.01\n")
    c = _SampleCounter()
    assert c.count(f) == 1
    f.unlink()
    f.write_text("Gen\tt_n7\tmu\n")
    assert c.count(f) == 0
    with open(f, "a") as fh:
        fh.write("10\t0.5\t0.02\n")
    assert c.count(f) == 1
