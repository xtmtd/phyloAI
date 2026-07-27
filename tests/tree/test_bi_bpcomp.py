from __future__ import annotations

from pathlib import Path

import pytest

from phyloai.tree.bi_bpcomp import (
    _discover_chain_names,
    _validate_chain_names,
    run_bi_bpcomp,
)
from phyloai.tree.bi import _build_x_flag


def test_build_x_flag_default():
    assert _build_x_flag(1000, 1, "all") == ["-x", "1000"]


def test_build_x_flag_with_sample_freq():
    assert _build_x_flag(1000, 10, "all") == ["-x", "1000", "10"]


def test_build_x_flag_with_until():
    assert _build_x_flag(1000, 1, "5000") == ["-x", "1000", "1", "5000"]


def test_build_x_flag_full():
    assert _build_x_flag(1000, 10, "5000") == ["-x", "1000", "10", "5000"]


def test_discover_chain_names(tmp_path: Path):
    for name in ["chain1", "chain2", "chain3"]:
        (tmp_path / f"{name}.chain").write_text("")
    (tmp_path / "chain1.trace").write_text("")
    result = _discover_chain_names(tmp_path)
    assert result == ["chain1", "chain2", "chain3"]


def test_discover_chain_names_empty(tmp_path: Path):
    with pytest.raises(FileNotFoundError, match="No .chain files"):
        _discover_chain_names(tmp_path)


def test_validate_chain_names_happy(tmp_path: Path):
    (tmp_path / "chain1.chain").write_text("")
    _validate_chain_names(tmp_path, ["chain1"])


def test_validate_chain_names_missing(tmp_path: Path):
    with pytest.raises(FileNotFoundError, match="not found"):
        _validate_chain_names(tmp_path, ["chain1", "chain99"])


def test_run_bi_bpcomp_dry_run(tmp_path: Path):
    chain_dir = tmp_path / "chains"
    chain_dir.mkdir()
    (chain_dir / "chain1.chain").write_text("")
    (chain_dir / "chain2.chain").write_text("")
    out = tmp_path / "bpcomp"

    payload = run_bi_bpcomp(chain_dir=chain_dir, output_dir=out, burnin=1000, dry_run=True)
    assert payload["status"] == "success"
    assert payload["key_results"]["chains_used"] == ["chain1", "chain2"]
    assert payload["key_results"]["bpcomp_maxdiff"] is None
    assert isinstance(payload["data"]["cmd"], list)


def test_run_bi_bpcomp_validation_negative_burnin():
    with pytest.raises(ValueError, match="burnin"):
        run_bi_bpcomp(chain_dir=Path("."), burnin=-1)


def test_run_bi_bpcomp_validation_zero_sample_freq():
    with pytest.raises(ValueError, match="sample-freq"):
        run_bi_bpcomp(chain_dir=Path("."), sample_freq=0)


def test_run_bi_bpcomp_validation_invalid_cutoff():
    with pytest.raises(ValueError, match="cutoff"):
        run_bi_bpcomp(chain_dir=Path("."), cutoff=1.0)


def test_run_bi_bpcomp_validation_invalid_until():
    with pytest.raises(ValueError, match="until"):
        run_bi_bpcomp(chain_dir=Path("."), until="abc")


def test_run_bi_bpcomp_validation_until_zero():
    with pytest.raises(ValueError, match="until"):
        run_bi_bpcomp(chain_dir=Path("."), until="0")


def test_run_bi_bpcomp_validation_until_int_ok(tmp_path: Path):
    chain_dir = tmp_path / "chains"
    chain_dir.mkdir()
    (chain_dir / "chain1.chain").write_text("")
    result = run_bi_bpcomp(chain_dir=chain_dir, output_dir=tmp_path / "bpcomp", until="5000", dry_run=True)
    assert result["status"] == "success"


def test_run_bi_bpcomp_validation_empty_chain_names(tmp_path: Path):
    chain_dir = tmp_path / "chains"
    chain_dir.mkdir()
    with pytest.raises(ValueError, match="at least one"):
        run_bi_bpcomp(chain_dir=chain_dir, chain_names=",,")


def test_run_bi_bpcomp_with_fake_tool(tmp_path: Path):
    chain_dir = tmp_path / "chains"
    chain_dir.mkdir()
    for name in ["chain1", "chain2"]:
        (chain_dir / f"{name}.chain").write_text("fake state")
        (chain_dir / f"{name}.treelist").write_text("(a,b);\n")

    tool_dir = tmp_path / "tools"
    tool_dir.mkdir()
    bpcomp = tool_dir / "bpcomp"
    bpcomp.write_text("""#!/bin/sh
while [ $# -gt 0 ]; do
    if [ "$1" = "-o" ]; then shift; out="$1"; fi
    shift
done
printf 'maxdiff 0.043\\nmeandiff 0.003\\n' > "${out}.bpdiff"
printf '(a,b,c);\\n' > "${out}.bplist"
printf '(a,(b,c));\\n' > "${out}.con.tre"
""")
    bpcomp.chmod(0o755)

    out = tmp_path / "bpcomp"
    payload = run_bi_bpcomp(
        chain_dir=chain_dir, output_dir=out, burnin=1000,
        pb_path=tool_dir, quiet=True,
    )
    assert payload["status"] == "success"
    assert payload["key_results"]["bpcomp_maxdiff"] == 0.043
    assert payload["key_results"]["bpcomp_meandiff"] == 0.003
    assert payload["key_results"]["bpcomp_status"] == "good"
    assert payload["data"]["tool_stderr"] == ""
    assert "bpcomp_all" in str(payload["key_results"]["consensus_tree"]) or "bpcomp.con" in str(payload["key_results"]["consensus_tree"])
