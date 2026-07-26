from __future__ import annotations

from pathlib import Path

import pytest

from phyloai.tree.bi_tracecomp import (
    _annotate_tracecomp_output,
    _discover_trace_names,
    _validate_trace_names,
    run_bi_tracecomp,
)


def test_discover_trace_names(tmp_path: Path):
    for name in ["chain1", "chain2"]:
        (tmp_path / f"{name}.trace").write_text("")
    (tmp_path / "chain1.chain").write_text("")
    result = _discover_trace_names(tmp_path)
    assert result == ["chain1", "chain2"]


def test_discover_trace_names_empty(tmp_path: Path):
    with pytest.raises(FileNotFoundError, match="No .trace files"):
        _discover_trace_names(tmp_path)


def test_validate_trace_names_happy(tmp_path: Path):
    (tmp_path / "chain1.trace").write_text("")
    _validate_trace_names(tmp_path, ["chain1"])


def test_validate_trace_names_missing(tmp_path: Path):
    with pytest.raises(FileNotFoundError, match="not found"):
        _validate_trace_names(tmp_path, ["chain1", "chain99"])


def test_annotate_tracecomp_output():
    stdout = "name effsize rel_diff\nloglik 1529 0.0524634\nlength 948 0.0360727\nalpha 1863 0.0694708\nnmode 2385 0.0170518\nbad_one 40 0.5\nok_one 200 0.25\n"
    annotated, min_effsize, max_rel_diff = _annotate_tracecomp_output(stdout)
    assert min_effsize == 40.0
    assert max_rel_diff == 0.5
    assert "[good]" in annotated
    assert "[no]" in annotated
    assert "[ok]" in annotated
    lines = annotated.splitlines()
    assert lines[0].strip().endswith("status")


def test_annotate_tracecomp_output_empty():
    annotated, min_eff, max_rel = _annotate_tracecomp_output("")
    assert annotated == ""
    assert min_eff is None
    assert max_rel is None


def test_run_bi_tracecomp_dry_run(tmp_path: Path):
    chain_dir = tmp_path / "chains"
    chain_dir.mkdir()
    (chain_dir / "chain1.trace").write_text("")
    (chain_dir / "chain2.trace").write_text("")
    out = tmp_path / "tracecomp"

    payload = run_bi_tracecomp(chain_dir=chain_dir, output_dir=out, burnin=1000, dry_run=True)
    assert payload["status"] == "success"
    assert payload["key_results"]["chains_used"] == ["chain1", "chain2"]
    assert payload["key_results"]["tracecomp_min_effsize"] is None
    assert isinstance(payload["data"]["cmd"], list)


def test_run_bi_tracecomp_validation_negative_burnin():
    with pytest.raises(ValueError, match="burnin"):
        run_bi_tracecomp(burnin=-1)


def test_run_bi_tracecomp_with_fake_tool(tmp_path: Path):
    chain_dir = tmp_path / "chains"
    chain_dir.mkdir()
    for name in ["chain1", "chain2"]:
        (chain_dir / f"{name}.trace").write_text("iter time loglik\n1 0 -10\n")

    tool_dir = tmp_path / "tools"
    tool_dir.mkdir()
    tracecomp = tool_dir / "tracecomp"
    tracecomp.write_text("""#!/bin/sh
printf 'name effsize rel_diff\\nloglik 1529 0.0524634\\nlength 948 0.0360727\\n'
""")
    tracecomp.chmod(0o755)

    out = tmp_path / "tracecomp"
    payload = run_bi_tracecomp(
        chain_dir=chain_dir, output_dir=out, burnin=500,
        pb_path=tool_dir, quiet=True,
    )
    assert payload["status"] == "success"
    assert payload["key_results"]["tracecomp_min_effsize"] == 948.0
    assert payload["key_results"]["tracecomp_max_reldiff"] == 0.0524634
    assert payload["key_results"]["tracecomp_status"] == "good"
    assert "loglik" in payload["data"]["tool_stderr"]
    assert (out / "tracecomp.contdiff").exists()
