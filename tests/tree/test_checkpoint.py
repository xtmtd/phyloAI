from __future__ import annotations

from pathlib import Path


def test_build_initial_checkpoint_tree(tmp_path: Path) -> None:
    from phyloai.tree.checkpoint_helpers import build_initial_checkpoint

    inputs = [Path("/data/gene1.fa"), Path("/data/gene2.fa")]
    trees_dir = Path("/out/trees")
    logs_dir = Path("/out/logs")

    ck = build_initial_checkpoint(
        step="tree.ml.fasttree",
        command="phyloai tree ml fasttree --msa-dir /data",
        params={"seq_type": "AA", "model": "lg"},
        inputs=inputs,
        trees_dir=trees_dir,
        logs_dir=logs_dir,
    )

    assert ck.schema_version == 1
    assert ck.step == "tree.ml.fasttree"
    assert ck.status == "running"
    assert len(ck.tasks) == 2
    assert ck.tasks[0].task_id == "gene1"
    assert ck.tasks[1].task_id == "gene2"
    assert ck.tasks[0].status == "pending"
    assert ck.tasks[0].outputs["tree"] == str(trees_dir / "gene1.tre")
    assert ck.tasks[0].outputs["log"] == str(logs_dir / "gene1.log")


def test_mark_task_updates_checkpoint(tmp_path: Path) -> None:
    from phyloai.tree.checkpoint_helpers import build_initial_checkpoint, mark_task

    inputs = [Path("/data/gene1.fa")]
    ck = build_initial_checkpoint(
        step="tree.ml.fasttree",
        command="cmd",
        params={},
        inputs=inputs,
        trees_dir=Path("/out/trees"),
        logs_dir=Path("/out/logs"),
    )

    mark_task(ck, "gene1", status="success")
    assert ck.tasks[0].status == "success"
    assert ck.tasks[0].attempts == 1

    mark_task(ck, "gene1", status="failed", reason="FastTree error")
    assert ck.tasks[0].status == "failed"
    assert ck.tasks[0].reason == "FastTree error"


def test_resume_verifier_valid_newick(tmp_path: Path) -> None:
    from phyloai.tree.checkpoint_helpers import resume_verifier

    tree_path = tmp_path / "gene1.tre"
    tree_path.write_text("(a:0.1,b:0.2);\n")

    verify = resume_verifier()
    assert verify(tree_path) is True


def test_resume_verifier_invalid_newick(tmp_path: Path) -> None:
    from phyloai.tree.checkpoint_helpers import resume_verifier

    tree_path = tmp_path / "gene1.tre"
    tree_path.write_text("(a:1,b:2;)\n")

    verify = resume_verifier()
    assert verify(tree_path) is False


def test_resume_verifier_empty_file(tmp_path: Path) -> None:
    from phyloai.tree.checkpoint_helpers import resume_verifier

    tree_path = tmp_path / "gene1.tre"
    tree_path.write_text("")

    verify = resume_verifier()
    assert verify(tree_path) is False


def test_resume_verifier_nonexistent_file(tmp_path: Path) -> None:
    from phyloai.tree.checkpoint_helpers import resume_verifier

    tree_path = tmp_path / "missing.tre"

    verify = resume_verifier()
    assert verify(tree_path) is False


def test_plan_resume_splits_tasks(tmp_path: Path) -> None:
    from phyloai.tree.checkpoint_helpers import build_initial_checkpoint, mark_task, plan_resume

    inputs = [Path("/data/g1.fa"), Path("/data/g2.fa"), Path("/data/g3.fa")]
    ck = build_initial_checkpoint(
        step="tree.ml.fasttree",
        command="cmd",
        params={},
        inputs=inputs,
        trees_dir=tmp_path / "trees",
        logs_dir=tmp_path / "logs",
    )

    (tmp_path / "trees").mkdir(parents=True)
    (tmp_path / "trees" / "g1.tre").write_text("(a:0.1,b:0.2);\n")

    mark_task(ck, "g1", status="success")
    mark_task(ck, "g2", status="failed", reason="error")
    mark_task(ck, "g3", status="pending")

    to_run, skipped = plan_resume(ck)

    assert "g2" in to_run
    assert "g3" in to_run
    assert "g1" in skipped


def test_plan_resume_all_succeeded_skips_all(tmp_path: Path) -> None:
    from phyloai.tree.checkpoint_helpers import build_initial_checkpoint, mark_task, plan_resume

    inputs = [Path("/data/g1.fa"), Path("/data/g2.fa")]
    ck = build_initial_checkpoint(
        step="tree.ml.fasttree",
        command="cmd",
        params={},
        inputs=inputs,
        trees_dir=tmp_path / "trees",
        logs_dir=tmp_path / "logs",
    )

    (tmp_path / "trees").mkdir(parents=True)
    (tmp_path / "trees" / "g1.tre").write_text("(a:0.1,b:0.2);\n")
    (tmp_path / "trees" / "g2.tre").write_text("(c:0.3,d:0.4);\n")

    mark_task(ck, "g1", status="success")
    mark_task(ck, "g2", status="success")

    to_run, skipped = plan_resume(ck)

    assert to_run == []
    assert "g1" in skipped
    assert "g2" in skipped
