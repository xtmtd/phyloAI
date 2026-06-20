from __future__ import annotations

from pathlib import Path

from click.testing import CliRunner

from phyloai.cli.main import cli


def test_tree_msc_help_shows_all_flags() -> None:
    result = CliRunner().invoke(cli, ["tree", "msc", "--help"])
    assert result.exit_code == 0
    for flag in [
        "--tree", "--tree-dir", "--mode", "--boot",
        "--extra-rounds", "--tree-boot-type", "--tree-boot-min",
        "--tree-boot-max", "--output-dir", "--threads",
        "--wastral-path", "--tool-args", "--overwrite",
        "--dry-run", "--quiet",
    ]:
        assert flag in result.output


def test_tree_msc_mutual_exclusivity(tmp_path: Path) -> None:
    tree_file = tmp_path / "g.trees"
    tree_file.write_text("((a,b),c);\n")
    tree_dir = tmp_path / "genetrees"
    tree_dir.mkdir()

    result = CliRunner().invoke(cli, [
        "tree", "msc",
        "--tree", str(tree_file),
        "--tree-dir", str(tree_dir),
    ])
    assert result.exit_code == 1


def test_tree_msc_neither_input() -> None:
    result = CliRunner().invoke(cli, ["tree", "msc"])
    assert result.exit_code == 1


def test_tree_msc_tree_nonexistent() -> None:
    result = CliRunner().invoke(cli, [
        "tree", "msc", "--tree", "/nonexistent/file.trees",
    ])
    assert result.exit_code == 1


def test_tree_msc_tree_dir_nonexistent() -> None:
    result = CliRunner().invoke(cli, [
        "tree", "msc", "--tree-dir", "/nonexistent/dir",
    ])
    assert result.exit_code == 1


def test_tree_msc_dry_run_tree_single(tmp_path: Path) -> None:
    tree_file = tmp_path / "g.trees"
    tree_file.write_text("((a,b),c);\n")
    out_dir = tmp_path / "out"

    result = CliRunner().invoke(cli, [
        "tree", "msc",
        "--tree", str(tree_file),
        "--output-dir", str(out_dir),
        "--dry-run",
    ])

    assert result.exit_code == 0
    assert "Dry run" in result.output


def test_tree_msc_dry_run_tree_dir(tmp_path: Path) -> None:
    tree_dir = tmp_path / "genetrees"
    tree_dir.mkdir()
    (tree_dir / "a.nwk").write_text("((a,b),c);\n")
    (tree_dir / "b.tre").write_text("((x,y),z);\n")
    out_dir = tmp_path / "out"

    result = CliRunner().invoke(cli, [
        "tree", "msc",
        "--tree-dir", str(tree_dir),
        "--output-dir", str(out_dir),
        "--dry-run",
    ])

    assert result.exit_code == 0
    assert "Dry run" in result.output
    assert not (out_dir / "merged.trees").exists()  # dry-run must not write files


def test_tree_msc_invalid_mode(tmp_path: Path) -> None:
    tree_file = tmp_path / "g.trees"
    tree_file.write_text("((a,b),c);\n")

    result = CliRunner().invoke(cli, [
        "tree", "msc",
        "--tree", str(tree_file),
        "--mode", "5",
    ])
    assert result.exit_code != 0


def test_tree_msc_invalid_boot(tmp_path: Path) -> None:
    tree_file = tmp_path / "g.trees"
    tree_file.write_text("((a,b),c);\n")

    result = CliRunner().invoke(cli, [
        "tree", "msc",
        "--tree", str(tree_file),
        "--boot", "4",
    ])
    assert result.exit_code != 0


def test_tree_msc_tree_boot_min_max_with_auto(tmp_path: Path) -> None:
    tree_file = tmp_path / "g.trees"
    tree_file.write_text("((a,b),c);\n")
    out_dir = tmp_path / "out"

    result = CliRunner().invoke(cli, [
        "tree", "msc",
        "--tree", str(tree_file),
        "--output-dir", str(out_dir),
        "--tree-boot-type", "auto",
        "--tree-boot-min", "10",
    ])
    assert result.exit_code == 1


def test_tree_msc_tool_args_blocked_i(tmp_path: Path) -> None:
    tree_file = tmp_path / "g.trees"
    tree_file.write_text("((a,b),c);\n")
    out_dir = tmp_path / "out"

    result = CliRunner().invoke(cli, [
        "tree", "msc",
        "--tree", str(tree_file),
        "--output-dir", str(out_dir),
        "--tool-args", "-i other.trees",
        "--dry-run",
    ])
    assert result.exit_code == 1


def test_tree_msc_tool_args_override_u(tmp_path: Path) -> None:
    tree_file = tmp_path / "g.trees"
    tree_file.write_text("((a,b),c);\n")
    out_dir = tmp_path / "out"

    result = CliRunner().invoke(cli, [
        "tree", "msc",
        "--tree", str(tree_file),
        "--output-dir", str(out_dir),
        "--tool-args", "-u 3",
        "--dry-run",
        "--quiet",
    ])
    assert result.exit_code == 0


def test_tree_msc_outgroup_in_cmd(tmp_path: Path) -> None:
    tree_file = tmp_path / "g.trees"
    tree_file.write_text("((a,b),c);\n")
    out_dir = tmp_path / "out"

    result = CliRunner().invoke(cli, [
        "tree", "msc",
        "--tree", str(tree_file),
        "--output-dir", str(out_dir),
        "--outgroup", "Gallus_gallus",
        "--dry-run",
    ])

    assert result.exit_code == 0
    assert "--root" in result.output
    assert "Gallus_gallus" in result.output
