from __future__ import annotations

from pathlib import Path

import pytest

from phyloai.tree.bi import (
    RESUME_ALL,
    _build_chain_cmd,
    _build_model_flags,
    _build_resume_cmd,
    _count_trace_samples,
    _parse_bpcomp_bpdiff,
    _parse_tracecomp_contdiff,
    _read_run_state,
    _resolve_chain_names,
    _resolve_resume_names,
    _state_payload,
    _status_from_metrics,
    _update_run_state_for_new_chains,
    _write_run_state,
    run_bi_pb,
)


# ---------------------------------------------------------------------------
# Task 2: pure helpers
# ---------------------------------------------------------------------------


def test_resolve_chain_names_auto():
    assert _resolve_chain_names(3, "chain", None) == ["chain1", "chain2", "chain3"]


def test_resolve_chain_names_default_chains_count():
    assert _resolve_chain_names(1, "chain", None) == ["chain1"]


def test_resolve_chain_names_explicit():
    assert _resolve_chain_names(3, "chain", "a,b") == ["a", "b"]


def test_resolve_chain_names_rejects_duplicates():
    with pytest.raises(ValueError, match="duplicate"):
        _resolve_chain_names(3, "chain", "a,a")


def test_resolve_chain_names_rejects_empty():
    with pytest.raises(ValueError, match="non-empty"):
        _resolve_chain_names(3, "chain", ",")


def test_resolve_chain_names_rejects_chains_zero():
    with pytest.raises(ValueError, match="at least 1"):
        _resolve_chain_names(0, "chain", None)


def test_build_model_flags_default():
    assert _build_model_flags("gtr", "auto", 4, None, None) == ["-cat", "-gtr", "-dgam", "4"]


def test_build_model_flags_homogeneous_lg():
    assert _build_model_flags("lg", "1", 4, None, None) == ["-ncat", "1", "-lg", "-dgam", "4"]


def test_build_model_flags_fixed_mixture():
    assert _build_model_flags("wag", "20", 4, None, None) == ["-ncat", "20", "-wag", "-dgam", "4"]


def test_build_model_flags_rejects_invalid_model():
    with pytest.raises(ValueError, match="Invalid --model"):
        _build_model_flags("blosum62", "auto", 4, None, None)


def test_build_model_flags_rejects_gamma_zero():
    with pytest.raises(ValueError, match="gamma-cats"):
        _build_model_flags("gtr", "auto", 0, None, None)


def test_build_model_flags_rejects_both_trees(tmp_path: Path):
    tree = tmp_path / "t.nwk"
    tree.write_text("(a,b);\n")
    with pytest.raises(ValueError, match="mutually exclusive"):
        _build_model_flags("gtr", "auto", 4, tree, tree)


def test_build_model_flags_with_start_tree(tmp_path: Path):
    tree = tmp_path / "t.nwk"
    tree.write_text("(a,b);\n")
    flags = _build_model_flags("gtr", "auto", 4, tree, None)
    assert "-t" in flags
    assert "-T" not in flags


def test_build_model_flags_with_fix_tree(tmp_path: Path):
    tree = tmp_path / "t.nwk"
    tree.write_text("(a,b);\n")
    flags = _build_model_flags("gtr", "auto", 4, None, tree)
    assert "-T" in flags
    assert "-t" not in flags


def test_build_model_flags_rejects_invalid_mixture_str():
    with pytest.raises(ValueError, match="must be 'auto' or a positive integer"):
        _build_model_flags("gtr", "hello", 4, None, None)


def test_build_model_flags_rejects_mixture_zero():
    with pytest.raises(ValueError, match="at least 1"):
        _build_model_flags("gtr", "0", 4, None, None)


def test_count_trace_samples_empty(tmp_path: Path):
    trace = tmp_path / "chain1.trace"
    trace.write_text("")
    assert _count_trace_samples(trace) == 0


def test_count_trace_samples_header_only(tmp_path: Path):
    trace = tmp_path / "chain1.trace"
    trace.write_text("iter time loglik\n")
    assert _count_trace_samples(trace) == 0


def test_count_trace_samples_with_data(tmp_path: Path):
    trace = tmp_path / "chain1.trace"
    trace.write_text("iter\ttime\tloglik\n1\t0\t-10\n2\t1\t-9\n")
    assert _count_trace_samples(trace) == 2


def test_count_trace_samples_ignores_partial_line(tmp_path: Path):
    trace = tmp_path / "chain1.trace"
    trace.write_bytes(b"iter\ttime\tloglik\n1\t0\t-10\n2\t1\t-9")
    assert _count_trace_samples(trace) == 1


def test_count_trace_samples_nonexistent(tmp_path: Path):
    assert _count_trace_samples(tmp_path / "nope.trace") == 0


# ---------------------------------------------------------------------------
# Task 3: command builders
# ---------------------------------------------------------------------------


def test_build_chain_cmd_forever(tmp_path: Path):
    matrix = tmp_path / "m.phy"
    matrix.write_text("2 3\na AAA\nb AAA\n")
    cmd = _build_chain_cmd("mpirun", "pb_mpi", 4, matrix, ["-cat", "-gtr", "-dgam", "4"], 1, -1, "chain1")
    assert cmd[-3:] == ["1", "-1", "chain1"]


def test_build_chain_cmd_with_target(tmp_path: Path):
    matrix = tmp_path / "m.phy"
    matrix.write_text("2 3\na AAA\nb AAA\n")
    cmd = _build_chain_cmd("mpirun", "pb_mpi", 4, matrix, ["-ncat", "1", "-lg", "-dgam", "4"], 1, 10000, "chain1")
    assert cmd[-3:] == ["1", "10000", "chain1"]


def test_build_chain_cmd_rejects_threads_one():
    with pytest.raises(ValueError, match="threads"):
        _build_chain_cmd("mpirun", "pb_mpi", 1, Path("m.phy"), [], 1, -1, "c1")


def test_build_chain_cmd_rejects_sample_freq_zero():
    with pytest.raises(ValueError, match="sample-freq"):
        _build_chain_cmd("mpirun", "pb_mpi", 2, Path("m.phy"), [], 0, -1, "c1")


def test_build_chain_cmd_rejects_nsamples_zero():
    with pytest.raises(ValueError, match="nsamples"):
        _build_chain_cmd("mpirun", "pb_mpi", 2, Path("m.phy"), [], 1, 0, "c1")


def test_build_resume_cmd():
    assert _build_resume_cmd("mpirun", "pb_mpi", 4, "chain1") == ["mpirun", "-np", "4", "pb_mpi", "chain1"]


def test_build_resume_cmd_rejects_threads_one():
    with pytest.raises(ValueError, match="threads"):
        _build_resume_cmd("mpirun", "pb_mpi", 1, "c1")


# ---------------------------------------------------------------------------
# Task 4: run state
# ---------------------------------------------------------------------------


def test_run_state_roundtrip(tmp_path: Path):
    matrix = tmp_path / "m.phy"
    matrix.write_text("2 3\na AAA\nb AAA\n")
    payload = _state_payload(["chain1"], matrix, ["-cat", "-gtr", "-dgam", "4"], 1, -1, 4)
    _write_run_state(tmp_path, payload)
    assert _read_run_state(tmp_path)["chain_names"] == ["chain1"]


def test_read_run_state_missing(tmp_path: Path):
    with pytest.raises(ValueError, match="Missing run_state.json"):
        _read_run_state(tmp_path)


def test_update_run_state_adds_new_chain(tmp_path: Path):
    matrix = tmp_path / "m.phy"
    matrix.write_text("2 3\na AAA\nb AAA\n")
    payload = _state_payload(["chain1"], matrix, ["-cat", "-gtr", "-dgam", "4"], 1, -1, 4)
    _write_run_state(tmp_path, payload)
    updated = _update_run_state_for_new_chains(tmp_path, ["chain2"], payload)
    assert updated["chain_names"] == ["chain1", "chain2"]


def test_update_run_state_rejects_duplicate_chain(tmp_path: Path):
    matrix = tmp_path / "m.phy"
    matrix.write_text("2 3\na AAA\nb AAA\n")
    payload = _state_payload(["chain1"], matrix, ["-cat", "-gtr", "-dgam", "4"], 1, -1, 4)
    _write_run_state(tmp_path, payload)
    with pytest.raises(ValueError, match="already exist"):
        _update_run_state_for_new_chains(tmp_path, ["chain1"], payload)


def test_update_run_state_rejects_mismatched_model(tmp_path: Path):
    matrix = tmp_path / "m.phy"
    matrix.write_text("2 3\na AAA\nb AAA\n")
    payload_1 = _state_payload(["chain1"], matrix, ["-cat", "-gtr", "-dgam", "4"], 1, -1, 4)
    payload_2 = _state_payload(["chain2"], matrix, ["-ncat", "1", "-lg", "-dgam", "4"], 1, -1, 4)
    _write_run_state(tmp_path, payload_1)
    with pytest.raises(ValueError, match="conflict"):
        _update_run_state_for_new_chains(tmp_path, ["chain2"], payload_2)


def test_resolve_resume_all():
    state = {"chain_names": ["chain1", "chain2"]}
    assert _resolve_resume_names(RESUME_ALL, state) == ["chain1", "chain2"]


def test_resolve_resume_subset():
    state = {"chain_names": ["chain1", "chain2", "chain3"]}
    assert _resolve_resume_names("chain1,chain3", state) == ["chain1", "chain3"]


def test_resolve_resume_missing():
    state = {"chain_names": ["chain1"]}
    with pytest.raises(ValueError, match="not found"):
        _resolve_resume_names("chain99", state)


def test_resolve_resume_none():
    state = {"chain_names": ["chain1"]}
    assert _resolve_resume_names(None, state) == []


# ---------------------------------------------------------------------------
# Task 5: convergence parsers
# ---------------------------------------------------------------------------


def test_parse_bpcomp_bpdiff(tmp_path: Path):
    path = tmp_path / "bpcomp_all.bpdiff"
    path.write_text("maxdiff 0.081\nmeandiff 0.006\n")
    assert _parse_bpcomp_bpdiff(path) == {"maxdiff": 0.081, "meandiff": 0.006}


def test_parse_bpcomp_bpdiff_missing(tmp_path: Path):
    assert _parse_bpcomp_bpdiff(tmp_path / "nope.bpdiff") == {"maxdiff": None, "meandiff": None}


def test_parse_tracecomp_contdiff(tmp_path: Path):
    path = tmp_path / "tracecomp_all.contdiff"
    path.write_text("name effsize rel_diff\nloglik 312 0.094\nlength 400 0.050\n")
    result = _parse_tracecomp_contdiff(path)
    assert result["min_effsize"] == 312.0
    assert result["max_rel_diff"] == 0.094


def test_parse_tracecomp_contdiff_missing(tmp_path: Path):
    result = _parse_tracecomp_contdiff(tmp_path / "nope.contdiff")
    assert result == {"min_effsize": None, "max_rel_diff": None}


def test_status_from_metrics_good():
    assert _status_from_metrics(0.081, 312, 0.094) == "good"


def test_status_from_metrics_ok():
    assert _status_from_metrics(0.25, 200, 0.25) == "ok"


def test_status_from_metrics_not_converged():
    assert _status_from_metrics(0.5, 40, 0.5) == "not converged"


def test_status_from_metrics_none():
    assert _status_from_metrics(None, 312, 0.094) == "not converged"


# ---------------------------------------------------------------------------
# Task 7: run_bi_pb dry-run
# ---------------------------------------------------------------------------


def test_run_bi_dry_run_result_json_shape(tmp_path: Path):
    from tests.helpers import validate_params_completeness, validate_result_json

    matrix = tmp_path / "m.phy"
    matrix.write_text("2 3\na AAA\nb AAA\n")
    payload = run_bi_pb(matrix=matrix, output_dir=tmp_path / "out", dry_run=True)
    assert payload["status"] == "success"
    assert payload["data"]["chain_cmds"]["chain1"][-1] == "chain1"
    assert isinstance(payload["data"]["tool_stderr"], dict)
    validate_result_json(payload)
    validate_params_completeness(payload, {
        "matrix", "output_dir", "overwrite", "model", "mixture", "gamma_cats",
        "start_tree", "fix_tree", "chains", "chain_prefix", "chain_names",
        "threads", "sample_freq", "nsamples", "resume", "monitor_freq",
        "burnin_frac", "poll_interval", "pb_path", "dry_run", "quiet",
    })


def test_run_bi_requires_matrix_without_resume():
    with pytest.raises(ValueError, match="--matrix is required"):
        run_bi_pb(matrix=None)


def test_run_bi_overwrite_and_resume_mutually_exclusive():
    with pytest.raises(ValueError, match="mutually exclusive"):
        run_bi_pb(matrix=Path("m.phy"), resume=RESUME_ALL, overwrite=True, dry_run=True)


def test_run_bi_rejects_missing_matrix(tmp_path: Path):
    with pytest.raises(ValueError, match="does not exist"):
        run_bi_pb(matrix=tmp_path / "nope.phy")


def test_run_bi_dry_run_fasta_produces_chain_cmds(tmp_path: Path):
    matrix = tmp_path / "m.fa"
    matrix.write_text(">a\nAAA\n>b\nAAA\n")
    payload = run_bi_pb(matrix=matrix, output_dir=tmp_path / "out", dry_run=True)
    assert payload["status"] == "success"
    assert "chain1" in payload["data"]["chain_cmds"]


# ---------------------------------------------------------------------------
# Resume --nsamples override tests
# ---------------------------------------------------------------------------


def test_resume_bare_preserves_stored_nsamples(tmp_path: Path):
    matrix = tmp_path / "m.phy"
    matrix.write_text("2 3\na AAA\nb AAA\n")
    out = tmp_path / "out"
    (out / "chains").mkdir(parents=True)
    (out / "convergence").mkdir()
    (out / "chains" / "chain1.trace").write_text("iter time loglik\n1 0 -10\n")
    (out / "run_state.json").write_text(
        '{"chain_names": ["chain1"], "matrix": "' + str(matrix) + '", '
        '"model_flags": ["-cat", "-gtr", "-dgam", "4"], '
        '"sample_freq": 1, "nsamples": 5000, "threads": 4, '
        '"model": "gtr", "mixture": "auto", "gamma_cats": 4, '
        '"start_tree": null, "fix_tree": null}'
    )
    payload = run_bi_pb(matrix=None, output_dir=out, resume=RESUME_ALL, dry_run=True)
    assert payload["params"]["nsamples"] == 5000


def test_resume_with_nsamples_override_fake_tools(tmp_path: Path):
    tool_dir = tmp_path / "tools"
    tool_dir.mkdir()
    mpirun = tool_dir / "mpirun"
    mpirun.write_text("#!/bin/sh\nshift 2\nexec \"$@\"\n")
    pb = tool_dir / "pb_mpi"
    pb.write_text("#!/bin/sh\nfor last do :; done\nname=\"$last\"\nprintf 'iter time loglik\\n1 0 -10\\n2 1 -9\\n' > ${name}.trace\nprintf '(a,b);\\n' > ${name}.treelist\nprintf 'state\\n' > ${name}.chain\nprintf 'stdout from pb\\n'\n")
    bp = tool_dir / "bpcomp"
    bp.write_text("#!/bin/sh\nwhile [ $# -gt 0 ]; do if [ \"$1\" = \"-o\" ]; then shift; out=$1; fi; shift; done\nprintf 'maxdiff 0.081\\nmeandiff 0.006\\n' > ${out}.bpdiff\nprintf '(a,b);\\n' > ${out}.con.tre\nprintf 'split\\n' > ${out}.bplist\n")
    tr = tool_dir / "tracecomp"
    tr.write_text("#!/bin/sh\nprintf 'name effsize rel_diff\\nloglik 312 0.094\\n'\n")
    for path in [mpirun, pb, bp, tr]:
        path.chmod(0o755)
    matrix = tmp_path / "m.phy"
    matrix.write_text("2 3\na AAA\nb AAA\n")
    out = tmp_path / "out"
    (out / "chains").mkdir(parents=True)
    (out / "convergence").mkdir()
    (out / "chains" / "chain1.trace").write_text("iter time loglik\n1 0 -10\n")
    (out / "run_state.json").write_text(
        '{"chain_names": ["chain1"], "matrix": "' + str(matrix) + '", '
        '"model_flags": ["-cat", "-gtr", "-dgam", "4"], '
        '"sample_freq": 1, "nsamples": 5000, "threads": 4, '
        '"model": "gtr", "mixture": "auto", "gamma_cats": 4, '
        '"start_tree": null, "fix_tree": null}'
    )
    payload = run_bi_pb(matrix=None, output_dir=out, resume=RESUME_ALL, nsamples=10000, pb_path=tool_dir, quiet=True)
    assert payload["params"]["nsamples"] == 10000
    import json
    reloaded = json.loads((out / "run_state.json").read_text())
    assert reloaded["nsamples"] == 10000


def test_resume_with_nsamples_minus_one(tmp_path: Path):
    matrix = tmp_path / "m.phy"
    matrix.write_text("2 3\na AAA\nb AAA\n")
    out = tmp_path / "out"
    (out / "chains").mkdir(parents=True)
    (out / "convergence").mkdir()
    (out / "chains" / "chain1.trace").write_text("iter time loglik\n1 0 -10\n")
    (out / "run_state.json").write_text(
        '{"chain_names": ["chain1"], "matrix": "' + str(matrix) + '", '
        '"model_flags": ["-cat", "-gtr", "-dgam", "4"], '
        '"sample_freq": 1, "nsamples": 5000, "threads": 4, '
        '"model": "gtr", "mixture": "auto", "gamma_cats": 4, '
        '"start_tree": null, "fix_tree": null}'
    )
    payload = run_bi_pb(matrix=None, output_dir=out, resume=RESUME_ALL, nsamples=-1, dry_run=True)
    assert payload["params"]["nsamples"] == -1


def test_resume_all_chains_already_at_target_no_crash(tmp_path: Path):
    matrix = tmp_path / "m.phy"
    matrix.write_text("2 3\na AAA\nb AAA\n")
    out = tmp_path / "out"
    (out / "chains").mkdir(parents=True)
    (out / "convergence").mkdir()
    (out / "chains" / "chain1.trace").write_text("iter time loglik\n1 0 -10\n2 1 -9\n")
    (out / "run_state.json").write_text(
        '{"chain_names": ["chain1"], "matrix": "' + str(matrix) + '", '
        '"model_flags": ["-cat", "-gtr", "-dgam", "4"], '
        '"sample_freq": 1, "nsamples": 2, "threads": 4, '
        '"model": "gtr", "mixture": "auto", "gamma_cats": 4, '
        '"start_tree": null, "fix_tree": null}'
    )
    payload = run_bi_pb(matrix=None, output_dir=out, resume=RESUME_ALL, quiet=True)
    assert payload["status"] == "success"
    assert "already at target" in payload["data"]["warnings"][0]


# ---------------------------------------------------------------------------
# Task 3-4: integration test with fake tools
# ---------------------------------------------------------------------------


def test_run_bi_fake_tools_executes_chains(tmp_path: Path):
    tool_dir = tmp_path / "tools"
    tool_dir.mkdir()
    mpirun = tool_dir / "mpirun"
    mpirun.write_text("#!/bin/sh\nshift 2\nexec \"$@\"\n")
    pb = tool_dir / "pb_mpi"
    pb.write_text("#!/bin/sh\nfor last do :; done\nname=\"$last\"\nprintf 'iter time loglik\\n1 0 -10\\n2 1 -9\\n' > ${name}.trace\nprintf '(a,b);\\n' > ${name}.treelist\nprintf 'state\\n' > ${name}.chain\nprintf 'stdout from pb\\n'\n")
    bp = tool_dir / "bpcomp"
    bp.write_text("#!/bin/sh\nwhile [ $# -gt 0 ]; do if [ \"$1\" = \"-o\" ]; then shift; out=$1; fi; shift; done\nprintf 'maxdiff 0.081\\nmeandiff 0.006\\n' > ${out}.bpdiff\nprintf '(a,b);\\n' > ${out}.con.tre\nprintf 'split\\n' > ${out}.bplist\n")
    tr = tool_dir / "tracecomp"
    tr.write_text("#!/bin/sh\nprintf 'name effsize rel_diff\\nloglik 312 0.094\\n'\n")
    for path in [mpirun, pb, bp, tr]:
        path.chmod(0o755)
    matrix = tmp_path / "m.phy"
    matrix.write_text("2 3\na AAA\nb AAA\n")
    payload = run_bi_pb(matrix=matrix, output_dir=tmp_path / "out", chains=2, nsamples=2, monitor_freq=1, burnin_frac=0.5, pb_path=tool_dir, quiet=True)
    assert payload["status"] == "success"
    assert payload["key_results"]["chain_lengths"] == {"chain1": 2, "chain2": 2}


# ---------------------------------------------------------------------------
# Task 6: convergence runner with fake tools
# ---------------------------------------------------------------------------


def test_run_convergence_check_with_fake_tools(tmp_path: Path):
    from phyloai.tree.bi import _run_convergence_check

    chains = tmp_path / "chains"
    chains.mkdir()
    for name in ["chain1", "chain2"]:
        (chains / f"{name}.trace").write_text("iter time loglik\n1 0 -10\n2 1 -9\n")
    bpcomp = tmp_path / "bpcomp"
    bpcomp.write_text("#!/bin/sh\nwhile [ $# -gt 0 ]; do if [ \"$1\" = \"-o\" ]; then shift; out=$1; fi; shift; done\nprintf 'maxdiff 0.081\\nmeandiff 0.006\\n' > ${out}.bpdiff\nprintf '(a,b);\\n' > ${out}.con.tre\nprintf 'split\\n' > ${out}.bplist\n")
    tracecomp = tmp_path / "tracecomp"
    tracecomp.write_text("#!/bin/sh\nprintf 'name effsize rel_diff\\nloglik 312 0.094\\n'\n")
    bpcomp.chmod(0o755)
    tracecomp.chmod(0o755)
    result = _run_convergence_check(tmp_path, ["chain1", "chain2"], {"bpcomp": str(bpcomp), "tracecomp": str(tracecomp)}, 1)
    assert result["all_chains"]["status"] == "good"
    assert (tmp_path / "convergence" / "bpcomp_all.bpdiff").exists()
    assert result["pairwise"]["chain1_chain2"]["status"] == "good"


def test_run_convergence_check_warns_on_failed_tool(tmp_path: Path):
    from phyloai.tree.bi import _run_convergence_check

    chains = tmp_path / "chains"
    chains.mkdir()
    for name in ["chain1"]:
        (chains / f"{name}.trace").write_text("iter time loglik\n1 0 -10\n")
    bpcomp = tmp_path / "bpcomp"
    bpcomp.write_text("#!/bin/sh\nexit 1\n")
    tracecomp = tmp_path / "tracecomp"
    tracecomp.write_text("#!/bin/sh\nprintf 'name effsize rel_diff\\nloglik 312 0.094\\n'\n")
    bpcomp.chmod(0o755)
    tracecomp.chmod(0o755)
    result = _run_convergence_check(tmp_path, ["chain1"], {"bpcomp": str(bpcomp), "tracecomp": str(tracecomp)}, 1)
    assert len(result["warnings"]) >= 1
    assert "exited with code" in result["warnings"][0]
