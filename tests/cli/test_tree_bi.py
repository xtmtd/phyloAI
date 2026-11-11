from __future__ import annotations

import json
import signal
from pathlib import Path

from click.testing import CliRunner

from phyloai.cli.main import cli


def test_tree_bi_pb_help_shows_flags():
    result = CliRunner().invoke(cli, ["tree", "bi", "pb", "--help"])
    assert result.exit_code == 0
    for flag in ["--matrix", "--model", "--mixture", "--gamma-cats", "--chains", "--threads", "--resume", "--pb-path", "--poll-interval"]:
        assert flag in result.output


def test_tree_group_shows_bi():
    result = CliRunner().invoke(cli, ["tree", "--help"])
    assert result.exit_code == 0
    assert "bi" in result.output


def test_tree_bi_group_help():
    result = CliRunner().invoke(cli, ["tree", "bi", "--help"])
    assert result.exit_code == 0
    for sub in ["pb", "bpcomp", "tracecomp", "readpb"]:
        assert sub in result.output


def test_tree_bi_bpcomp_help():
    result = CliRunner().invoke(cli, ["tree", "bi", "bpcomp", "--help"])
    assert result.exit_code == 0
    for flag in ["--chain-dir", "--burnin", "--cutoff"]:
        assert flag in result.output


def test_tree_bi_tracecomp_help():
    result = CliRunner().invoke(cli, ["tree", "bi", "tracecomp", "--help"])
    assert result.exit_code == 0
    for flag in ["--chain-dir", "--burnin"]:
        assert flag in result.output


def test_tree_bi_readpb_help():
    result = CliRunner().invoke(cli, ["tree", "bi", "readpb", "--help"])
    assert result.exit_code == 0
    for flag in ["--chain", "--mode", "--threads"]:
        assert flag in result.output


def test_tree_bi_pb_dry_run(tmp_path: Path):
    matrix = tmp_path / "m.phy"
    matrix.write_text("2 3\na AAA\nb AAA\n")
    result = CliRunner().invoke(cli, ["tree", "bi", "pb", "--matrix", str(matrix), "--output-dir", str(tmp_path / "out"), "--dry-run"])
    assert result.exit_code == 0
    assert "pb_mpi" in result.output


def test_tree_bi_pb_resume_bare_parses(tmp_path: Path):
    out = tmp_path / "out"
    out.mkdir()
    (out / "run_state.json").write_text('{"chain_names": ["chain1"], "matrix": "m.phy", "model_flags": ["-cat", "-gtr", "-dgam", "4"], "sample_freq": 1, "nsamples": -1, "threads": 4}')
    result = CliRunner().invoke(cli, ["tree", "bi", "pb", "--output-dir", str(out), "--resume", "--dry-run", "--quiet"])
    assert result.exit_code == 0


def test_tree_bi_pb_dry_run_fails_missing_matrix():
    result = CliRunner().invoke(cli, ["tree", "bi", "pb", "--matrix", "/nonexistent.phy", "--dry-run"])
    assert result.exit_code == 1
    assert "does not exist" in result.output


def test_tree_bi_pb_dry_run_phylogeny_with_fasta(tmp_path: Path):
    matrix = tmp_path / "test.fa"
    matrix.write_text(">a\nAAA\n>b\nAAA\n")
    result = CliRunner().invoke(cli, ["tree", "bi", "pb", "--matrix", str(matrix), "--output-dir", str(tmp_path / "out"), "--dry-run"])
    assert result.exit_code == 0
    assert "pb_mpi" in result.output


def test_tree_bi_bpcomp_dry_run(tmp_path: Path):
    chain_dir = tmp_path / "chains"
    chain_dir.mkdir()
    (chain_dir / "chain1.chain").write_text("")
    (chain_dir / "chain2.chain").write_text("")
    result = CliRunner().invoke(cli, ["tree", "bi", "bpcomp", "--chain-dir", str(chain_dir), "--output-dir", str(tmp_path / "bpcomp"), "--burnin", "1000", "--dry-run"])
    assert result.exit_code == 0


def test_tree_bi_tracecomp_dry_run(tmp_path: Path):
    chain_dir = tmp_path / "chains"
    chain_dir.mkdir()
    (chain_dir / "chain1.trace").write_text("")
    (chain_dir / "chain2.trace").write_text("")
    result = CliRunner().invoke(cli, ["tree", "bi", "tracecomp", "--chain-dir", str(chain_dir), "--output-dir", str(tmp_path / "tracecomp"), "--dry-run"])
    assert result.exit_code == 0


def test_tree_bi_readpb_dry_run(tmp_path: Path):
    chain = tmp_path / "chain1"
    (tmp_path / "chain1.chain").write_text("fake state")
    result = CliRunner().invoke(cli, ["tree", "bi", "readpb", "--chain", str(chain), "--mode", "rr", "--output-dir", str(tmp_path / "readpb"), "--dry-run"])
    assert result.exit_code == 0


def test_previous_result_is_not_backed_up_before_the_run(tmp_path: Path):
    """Backup belongs to final publication, not run startup."""
    out = tmp_path / "out"
    out.mkdir()
    (out / "result.json").write_text(json.dumps({"command": "phyloai previous run", "sentinel": True}))
    matrix = tmp_path / "m.phy"
    matrix.write_text("2 3\na AAA\nb AAA\n")

    result = CliRunner().invoke(
        cli,
        ["tree", "bi", "pb", "--matrix", str(matrix), "--output-dir", str(out), "--burnin-frac", "1.5"],
    )

    assert result.exit_code == 1
    assert not list(out.glob("result_*.json"))
    payload = json.loads((out / "result.json").read_text())
    assert payload["status"] == "error"


def _fake_payload(*, status: str = "success", shutdown_state: str = "completed", force_signal=None) -> dict:
    data = {
        "chain_cmds": {},
        "tool_stderr": {},
        "tool_logs": {},
        "output_files": {},
        "interrupted": shutdown_state != "completed",
        "warnings": [],
        "shutdown_state": shutdown_state,
    }
    if force_signal is not None:
        data["stop_signal"] = force_signal
        data["force_signal"] = force_signal
    return {
        "status": status,
        "command": "phyloai tree bi pb --matrix m.phy --output-dir out",
        "wall_time": 0.0,
        "tool_versions": {},
        "params": {},
        "key_results": {},
        "error": "forced" if status == "error" else None,
        "data": data,
    }


def _invoke_pb(tmp_path: Path):
    matrix = tmp_path / "m.phy"
    matrix.write_text("2 3\na AAA\nb AAA\n")
    return CliRunner().invoke(
        cli,
        ["tree", "bi", "pb", "--matrix", str(matrix), "--output-dir", str(tmp_path / "out")],
    )


def test_cli_backs_up_old_result_only_when_publishing(tmp_path: Path, monkeypatch):
    import phyloai.tree.bi as bi_module

    out = tmp_path / "out"
    out.mkdir()
    (out / "result.json").write_text('{"old": true}')
    monkeypatch.setattr(bi_module, "run_bi_pb", lambda **kwargs: _fake_payload())

    result = _invoke_pb(tmp_path)

    assert result.exit_code == 0
    assert list(out.glob("result_*.json"))
    assert json.loads((out / "result.json").read_text())["status"] == "success"


def test_cli_publish_failure_returns_4(tmp_path: Path, monkeypatch):
    import phyloai.cli.commands.tree as tree_module
    import phyloai.tree.bi as bi_module

    monkeypatch.setattr(bi_module, "run_bi_pb", lambda **kwargs: _fake_payload())

    def fail_publish(*_args, **_kwargs):
        raise OSError("disk full")

    monkeypatch.setattr(tree_module, "write_result_json", fail_publish)
    result = _invoke_pb(tmp_path)

    assert result.exit_code == 4
    assert "disk full" in result.output


def test_cli_publish_type_error_returns_4(tmp_path: Path, monkeypatch):
    import phyloai.cli.commands.tree as tree_module
    import phyloai.tree.bi as bi_module

    monkeypatch.setattr(bi_module, "run_bi_pb", lambda **kwargs: _fake_payload())

    def fail_publish(*_args, **_kwargs):
        raise TypeError("Object is not JSON serializable")

    monkeypatch.setattr(tree_module, "write_result_json", fail_publish)
    result = _invoke_pb(tmp_path)

    assert result.exit_code == 4
    assert "not JSON serializable" in result.output


def test_cli_publish_circular_value_error_returns_4_and_preserves_old_result(tmp_path: Path, monkeypatch):
    import phyloai.cli.commands.tree as tree_module
    import phyloai.tree.bi as bi_module

    out = tmp_path / "out"
    out.mkdir()
    (out / "result.json").write_text("old-result")
    monkeypatch.setattr(bi_module, "run_bi_pb", lambda **kwargs: _fake_payload())

    def fail_publish(*_args, **_kwargs):
        raise ValueError("Circular reference detected")

    monkeypatch.setattr(tree_module, "write_result_json", fail_publish)
    result = _invoke_pb(tmp_path)

    assert result.exit_code == 4
    assert "Circular reference" in result.output
    assert (out / "result.json").read_text() == "old-result"


def test_tool_failure_with_publish_failure_preserves_exit_2(tmp_path: Path, monkeypatch):
    import phyloai.cli.commands.tree as tree_module
    import phyloai.tree.bi as bi_module

    payload = _fake_payload(status="error")
    payload["error"] = "pb_mpi failed"
    payload["data"].pop("shutdown_state", None)
    monkeypatch.setattr(bi_module, "run_bi_pb", lambda **kwargs: payload)

    def fail_publish(*_args, **_kwargs):
        raise OSError("disk full")

    monkeypatch.setattr(tree_module, "write_result_json", fail_publish)
    result = _invoke_pb(tmp_path)

    assert result.exit_code == 2
    assert "disk full" in result.output


def test_missing_tool_error_write_failure_preserves_exit_3(tmp_path: Path, monkeypatch):
    import phyloai.cli.commands.tree as tree_module
    import phyloai.tree.bi as bi_module

    def missing(**_kwargs):
        raise FileNotFoundError("pb_mpi not found")

    monkeypatch.setattr(bi_module, "run_bi_pb", missing)

    def fail_write(*_args, **_kwargs):
        raise OSError("cannot write error result")

    monkeypatch.setattr(tree_module, "write_result_json", fail_write)
    result = _invoke_pb(tmp_path)

    assert result.exit_code == 3
    assert "pb_mpi not found" in result.output
    assert "cannot write error result" in result.output


def test_cli_backup_failure_returns_4_and_preserves_old_result(tmp_path: Path, monkeypatch):
    import phyloai.cli.commands.tree as tree_module
    import phyloai.tree.bi as bi_module

    out = tmp_path / "out"
    out.mkdir()
    (out / "result.json").write_text("old-result")
    monkeypatch.setattr(bi_module, "run_bi_pb", lambda **kwargs: _fake_payload())

    def fail_backup(*_args, **_kwargs):
        raise OSError("backup failed")

    monkeypatch.setattr(tree_module.shutil, "copy2", fail_backup)
    result = _invoke_pb(tmp_path)

    assert result.exit_code == 4
    assert "backup failed" in result.output
    assert (out / "result.json").read_text() == "old-result"


def test_cli_publishes_final_result_with_fsync(tmp_path: Path, monkeypatch):
    import phyloai.cli.commands.tree as tree_module
    import phyloai.tree.bi as bi_module

    monkeypatch.setattr(bi_module, "run_bi_pb", lambda **kwargs: _fake_payload())
    real_write = tree_module.write_result_json
    calls: list[bool] = []

    def capture(payload, output_dir, *, fsync=False):
        calls.append(fsync)
        return real_write(payload, output_dir, fsync=fsync)

    monkeypatch.setattr(tree_module, "write_result_json", capture)
    result = _invoke_pb(tmp_path)

    assert result.exit_code == 0
    assert calls == [True]


def test_cli_sidecar_delete_failure_warns_without_rewriting_result(tmp_path: Path, monkeypatch):
    import phyloai.tree.bi as bi_module

    monkeypatch.setattr(bi_module, "run_bi_pb", lambda **kwargs: _fake_payload())

    def fail_clear(_output_dir):
        raise OSError("cannot unlink sidecar")

    monkeypatch.setattr(bi_module, "_clear_stop_state", fail_clear)
    result = _invoke_pb(tmp_path)

    assert result.exit_code == 0
    assert "could not remove stop_state.json" in result.output
    assert json.loads((tmp_path / "out" / "result.json").read_text())["status"] == "success"


def test_cli_forced_publish_failure_preserves_143(tmp_path: Path, monkeypatch):
    import phyloai.cli.commands.tree as tree_module
    import phyloai.tree.bi as bi_module

    monkeypatch.setattr(
        bi_module,
        "run_bi_pb",
        lambda **kwargs: _fake_payload(
            status="error", shutdown_state="forced", force_signal=signal.SIGTERM
        ),
    )

    def fail_publish(*_args, **_kwargs):
        raise OSError("disk full")

    monkeypatch.setattr(tree_module, "write_result_json", fail_publish)
    result = _invoke_pb(tmp_path)

    assert result.exit_code == 143
    assert "disk full" in result.output


def test_cli_late_second_signal_is_finalized_before_publish(tmp_path: Path, monkeypatch):
    import phyloai.tree.bi as bi_module

    def fake_run(**kwargs):
        stop = kwargs["stop"]
        stop.note_signal(signal.SIGTERM)
        stop.mark_handled(1)
        stop.note_signal(signal.SIGINT)
        return _fake_payload(shutdown_state="stopped")

    monkeypatch.setattr(bi_module, "run_bi_pb", fake_run)
    result = _invoke_pb(tmp_path)

    assert result.exit_code == 130
    payload = json.loads((tmp_path / "out" / "result.json").read_text())
    assert payload["status"] == "error"
    assert payload["data"]["shutdown_state"] == "forced"
    assert payload["data"]["force_signal"] == signal.SIGINT
