from __future__ import annotations

from pathlib import Path

import pytest


_WASTRAL_EXTENSIONS = frozenset({
    ".nwk", ".tre", ".tree", ".nw", ".trees", ".newick",
})


# ===================================================================
# _scan_input_wastral
# ===================================================================

def test_scan_input_wastral_finds_all_supported(tmp_path: Path) -> None:
    from phyloai.tree.msc import _scan_input_wastral

    (tmp_path / "gene1.nwk").write_text("((a,b),c);\n")
    (tmp_path / "gene2.tre").write_text("((x,y),z);\n")
    (tmp_path / "gene3.tree").write_text("((1,2),3);\n")
    (tmp_path / "gene4.nw").write_text("(A,B);\n")
    (tmp_path / "gene5.trees").write_text("((a,b),c);\n")
    (tmp_path / "gene6.newick").write_text("(X,Y);\n")
    (tmp_path / "notes.txt").write_text("skip")
    (tmp_path / "empty.nwk").write_text("")
    (tmp_path / "subdir").mkdir()

    found, skipped = _scan_input_wastral(tmp_path)

    assert len(found) == 6
    assert len(skipped) == 3
    skip_reasons = {s["reason"] for s in skipped}
    assert "empty file" in skip_reasons
    assert "directory" in skip_reasons
    assert "unrecognized extension: .txt" in skip_reasons


def test_scan_input_wastral_nonexistent_dir() -> None:
    from phyloai.tree.msc import _scan_input_wastral

    found, skipped = _scan_input_wastral(Path("/nonexistent/dir"))
    assert found == []
    assert skipped == []


# ===================================================================
# _merge_gene_trees
# ===================================================================

def test_merge_gene_trees_concatenates(tmp_path: Path) -> None:
    from phyloai.tree.msc import _merge_gene_trees

    (tmp_path / "a.nwk").write_text("((a,b),c);\n")
    (tmp_path / "b.tre").write_text("((x,y),z);\n")

    merged_path = tmp_path / "merged.trees"
    count, skipped = _merge_gene_trees(tmp_path, merged_path)

    assert count == 2
    assert merged_path.exists()
    lines = merged_path.read_text().strip().split("\n")
    assert len(lines) == 2
    assert lines[0] == "((a,b),c);"
    assert lines[1] == "((x,y),z);"


def test_merge_gene_trees_multi_line_file(tmp_path: Path) -> None:
    from phyloai.tree.msc import _merge_gene_trees

    (tmp_path / "multi.trees").write_text("((a,b),c);\n((x,y),z);\n(A,B);\n")

    merged_path = tmp_path / "merged.trees"
    count, _ = _merge_gene_trees(tmp_path, merged_path)

    assert count == 3


def test_merge_gene_trees_skips_non_newick(tmp_path: Path) -> None:
    from phyloai.tree.msc import _merge_gene_trees

    (tmp_path / "a.nwk").write_text("((a,b),c);\n")
    (tmp_path / "notes.txt").write_text("not a tree")
    (tmp_path / "empty.tre").write_text("")

    merged_path = tmp_path / "merged.trees"
    count, skipped = _merge_gene_trees(tmp_path, merged_path)

    assert count == 1
    assert len(skipped) == 2
    skip_reasons = {s["reason"] for s in skipped}
    assert "unrecognized extension: .txt" in skip_reasons
    assert "empty file" in skip_reasons


def test_merge_gene_trees_no_valid_files(tmp_path: Path) -> None:
    from phyloai.tree.msc import _merge_gene_trees

    (tmp_path / "notes.txt").write_text("not a tree")
    (tmp_path / "empty.tre").write_text("")

    merged_path = tmp_path / "merged.trees"
    count, skipped = _merge_gene_trees(tmp_path, merged_path)

    assert count == 0
    assert len(skipped) == 2


# ===================================================================
# _build_wastral_cmd
# ===================================================================

def test_build_wastral_cmd_defaults(tmp_path: Path) -> None:
    from phyloai.tree.msc import _build_wastral_cmd

    inp = tmp_path / "input.trees"
    out = tmp_path / "output.tre"
    cmd = _build_wastral_cmd(inp, out, mode=1, boot=1, extra_rounds=False,
                              tree_boot_type="auto", tree_boot_min=None,
                              tree_boot_max=None, threads=4)

    assert cmd[0] == "wastral"
    assert "-i" in cmd
    assert str(inp) in cmd
    assert "-o" in cmd
    assert str(out) in cmd
    assert "--mode" in cmd
    assert "1" in cmd
    assert "-u" in cmd
    assert "-t" in cmd
    assert "4" in cmd
    assert "-R" not in cmd
    assert "--lrt" not in cmd
    assert "--bayes" not in cmd
    assert "--bootstrap" not in cmd


def test_build_wastral_cmd_mode_4_exhaustive(tmp_path: Path) -> None:
    from phyloai.tree.msc import _build_wastral_cmd

    inp = tmp_path / "input.trees"
    out = tmp_path / "output.tre"
    cmd = _build_wastral_cmd(inp, out, mode=4, boot=2, extra_rounds=True,
                              tree_boot_type="auto", tree_boot_min=None,
                              tree_boot_max=None, threads=8)

    assert "--mode" in cmd and "4" in cmd
    assert "-u" in cmd and "2" in cmd
    assert "-R" in cmd
    assert "-t" in cmd and "8" in cmd


def test_build_wastral_cmd_bootstrap_type(tmp_path: Path) -> None:
    from phyloai.tree.msc import _build_wastral_cmd

    inp = tmp_path / "input.trees"
    out = tmp_path / "output.tre"
    cmd = _build_wastral_cmd(inp, out, mode=1, boot=1, extra_rounds=False,
                              tree_boot_type="bootstrap", tree_boot_min=10,
                              tree_boot_max=95, threads=4)

    assert "--bootstrap" in cmd
    assert "-x" in cmd
    assert "95" in cmd
    assert "-n" in cmd
    assert "10" in cmd
    assert "-d" in cmd
    assert "0" in cmd  # hardcoded -d for bootstrap


def test_build_wastral_cmd_likelihood_type(tmp_path: Path) -> None:
    from phyloai.tree.msc import _build_wastral_cmd

    inp = tmp_path / "input.trees"
    out = tmp_path / "output.tre"
    cmd = _build_wastral_cmd(inp, out, mode=1, boot=1, extra_rounds=False,
                              tree_boot_type="likelihood", tree_boot_min=None,
                              tree_boot_max=None, threads=4)

    assert "--lrt" in cmd
    assert "-x" in cmd
    assert "1.0" in cmd
    assert "-n" in cmd
    assert "0.0" in cmd
    assert "-d" in cmd
    assert "0" in cmd


def test_build_wastral_cmd_abayes_type(tmp_path: Path) -> None:
    from phyloai.tree.msc import _build_wastral_cmd

    inp = tmp_path / "input.trees"
    out = tmp_path / "output.tre"
    cmd = _build_wastral_cmd(inp, out, mode=1, boot=1, extra_rounds=False,
                              tree_boot_type="abayes", tree_boot_min=None,
                              tree_boot_max=None, threads=4)

    assert "--bayes" in cmd
    assert "-x" in cmd
    assert "1.0" in cmd
    assert "-n" in cmd
    assert "0.333" in cmd
    assert "-d" in cmd
    assert "0.333" in cmd


# ===================================================================
# --tool-args management
# ===================================================================

def test_build_wastral_cmd_tool_args_blocks_minus_i(tmp_path: Path) -> None:
    from phyloai.tree.msc import _build_wastral_cmd

    inp = tmp_path / "input.trees"
    out = tmp_path / "output.tre"
    with pytest.raises(ValueError, match="Blocked managed flag"):
        _build_wastral_cmd(inp, out, mode=1, boot=1, extra_rounds=False,
                            tree_boot_type="auto", tree_boot_min=None,
                            tree_boot_max=None, threads=4,
                            tool_args="-i other.trees")


def test_build_wastral_cmd_tool_args_blocks_minus_o(tmp_path: Path) -> None:
    from phyloai.tree.msc import _build_wastral_cmd

    inp = tmp_path / "input.trees"
    out = tmp_path / "output.tre"
    with pytest.raises(ValueError, match="Blocked managed flag"):
        _build_wastral_cmd(inp, out, mode=1, boot=1, extra_rounds=False,
                            tree_boot_type="auto", tree_boot_min=None,
                            tree_boot_max=None, threads=4,
                            tool_args="-o output.tre")


def test_build_wastral_cmd_tool_args_overrides_u(tmp_path: Path) -> None:
    from phyloai.tree.msc import _build_wastral_cmd

    inp = tmp_path / "input.trees"
    out = tmp_path / "output.tre"
    cmd = _build_wastral_cmd(inp, out, mode=1, boot=1, extra_rounds=False,
                              tree_boot_type="auto", tree_boot_min=None,
                              tree_boot_max=None, threads=4,
                              tool_args="-u 3")

    # phyloAI should suppress its own -u, and -u 3 from tool-args should appear
    assert cmd.count("-u") >= 1
    assert "3" in cmd


def test_build_wastral_cmd_tool_args_overrides_support_flag(tmp_path: Path) -> None:
    from phyloai.tree.msc import _build_wastral_cmd

    inp = tmp_path / "input.trees"
    out = tmp_path / "output.tre"
    cmd = _build_wastral_cmd(inp, out, mode=1, boot=1, extra_rounds=False,
                              tree_boot_type="auto", tree_boot_min=None,
                              tree_boot_max=None, threads=4,
                              tool_args="--support 3")

    # When --support is in tool-args, phyloAI must suppress its own -u
    # count check: -u should NOT appear (suppressed), --support 3 from tool-args
    assert "-u" not in cmd
    assert "--support" in cmd
    assert "3" in cmd


def test_build_wastral_cmd_tool_args_overrides_R(tmp_path: Path) -> None:
    from phyloai.tree.msc import _build_wastral_cmd

    inp = tmp_path / "input.trees"
    out = tmp_path / "output.tre"
    cmd = _build_wastral_cmd(inp, out, mode=1, boot=1, extra_rounds=True,
                              tree_boot_type="auto", tree_boot_min=None,
                              tree_boot_max=None, threads=4,
                              tool_args="-R")

    # -R should appear exactly once (from tool-args, phyloAI suppressed)
    assert cmd.count("-R") == 1


def test_build_wastral_cmd_with_outgroup(tmp_path: Path) -> None:
    from phyloai.tree.msc import _build_wastral_cmd

    inp = tmp_path / "input.trees"
    out = tmp_path / "output.tre"
    cmd = _build_wastral_cmd(inp, out, mode=1, boot=1, extra_rounds=False,
                              tree_boot_type="auto", tree_boot_min=None,
                              tree_boot_max=None, threads=4,
                              outgroup="Arabidopsis")

    assert "--root" in cmd
    assert "Arabidopsis" in cmd


def test_build_wastral_cmd_no_outgroup(tmp_path: Path) -> None:
    from phyloai.tree.msc import _build_wastral_cmd

    inp = tmp_path / "input.trees"
    out = tmp_path / "output.tre"
    cmd = _build_wastral_cmd(inp, out, mode=1, boot=1, extra_rounds=False,
                              tree_boot_type="auto", tree_boot_min=None,
                              tree_boot_max=None, threads=4)

    assert "--root" not in cmd


# ===================================================================
# run_wastral validation
# ===================================================================

def test_run_wastral_mutual_exclusivity_both_none() -> None:
    from phyloai.tree.msc import run_wastral

    with pytest.raises(ValueError, match="Either --tree or --tree-dir"):
        run_wastral(tree=None, tree_dir=None)


def test_run_wastral_mutual_exclusivity_both(tmp_path: Path) -> None:
    from phyloai.tree.msc import run_wastral

    tree_file = tmp_path / "g.trees"
    tree_file.write_text("((a,b),c);\n")
    tree_dir = tmp_path / "genetrees"
    tree_dir.mkdir()
    (tree_dir / "a.nwk").write_text("((a,b),c);\n")

    with pytest.raises(ValueError, match="mutually exclusive"):
        run_wastral(tree=tree_file, tree_dir=tree_dir)


def test_run_wastral_tree_file_not_found() -> None:
    from phyloai.tree.msc import run_wastral

    with pytest.raises(ValueError, match="does not exist"):
        run_wastral(tree=Path("/nonexistent/file.trees"))


def test_run_wastral_tree_dir_no_valid_files(tmp_path: Path) -> None:
    from phyloai.tree.msc import run_wastral

    tree_dir = tmp_path / "empty"
    tree_dir.mkdir()
    out_dir = tmp_path / "out"

    with pytest.raises(ValueError, match="No valid gene tree files"):
        run_wastral(tree_dir=tree_dir, output_dir=out_dir)


def test_run_wastral_tree_boot_min_ge_max(tmp_path: Path) -> None:
    from phyloai.tree.msc import run_wastral

    tree_file = tmp_path / "g.trees"
    tree_file.write_text("((a,b),c);\n")

    with pytest.raises(ValueError, match="tree-boot-min.*<.*tree-boot-max"):
        run_wastral(
            tree=tree_file,
            tree_boot_type="bootstrap",
            tree_boot_min=95,
            tree_boot_max=10,
        )


def test_run_wastral_tree_boot_min_max_with_auto(tmp_path: Path) -> None:
    from phyloai.tree.msc import run_wastral

    tree_file = tmp_path / "g.trees"
    tree_file.write_text("((a,b),c);\n")

    with pytest.raises(ValueError, match="tree-boot-min.*only valid"):
        run_wastral(
            tree=tree_file,
            tree_boot_type="auto",
            tree_boot_min=10,
        )


def test_run_wastral_output_dir_exists(tmp_path: Path) -> None:
    from phyloai.tree.msc import run_wastral

    tree_file = tmp_path / "g.trees"
    tree_file.write_text("((a,b),c);\n")
    out_dir = tmp_path / "out"
    out_dir.mkdir()
    (out_dir / "existing.txt").write_text("old")

    with pytest.raises(ValueError, match="already exists"):
        run_wastral(tree=tree_file, output_dir=out_dir)


def test_run_wastral_invalid_mode(tmp_path: Path) -> None:
    from phyloai.tree.msc import run_wastral

    tree_file = tmp_path / "g.trees"
    tree_file.write_text("((a,b),c);\n")

    with pytest.raises(ValueError, match="mode must be"):
        run_wastral(tree=tree_file, mode=5)


def test_run_wastral_invalid_boot(tmp_path: Path) -> None:
    from phyloai.tree.msc import run_wastral

    tree_file = tmp_path / "g.trees"
    tree_file.write_text("((a,b),c);\n")

    with pytest.raises(ValueError, match="boot must be"):
        run_wastral(tree=tree_file, boot=4)


def test_run_wastral_invalid_threads(tmp_path: Path) -> None:
    from phyloai.tree.msc import run_wastral

    tree_file = tmp_path / "g.trees"
    tree_file.write_text("((a,b),c);\n")

    with pytest.raises(ValueError, match="threads must be"):
        run_wastral(tree=tree_file, threads=0)


# ===================================================================
# run_wastral dry-run
# ===================================================================

def test_run_wastral_dry_run_produces_payload(tmp_path: Path) -> None:
    from phyloai.tree.msc import run_wastral

    tree_file = tmp_path / "g.trees"
    tree_file.write_text("((a,b),c);\n")
    out_dir = tmp_path / "out"

    result = run_wastral(tree=tree_file, output_dir=out_dir, dry_run=True)

    assert result["status"] == "success"
    assert "cmd" in result["data"]
    assert "phyloai tree msc" in result["command"]
    assert not (out_dir / "wastral.tre").exists()


def test_run_wastral_dry_run_tree_dir(tmp_path: Path) -> None:
    from phyloai.tree.msc import run_wastral

    tree_dir = tmp_path / "genetrees"
    tree_dir.mkdir()
    (tree_dir / "a.nwk").write_text("((a,b),c);\n")
    (tree_dir / "b.tre").write_text("((x,y),z);\n")
    out_dir = tmp_path / "out"

    result = run_wastral(tree_dir=tree_dir, output_dir=out_dir, dry_run=True)

    assert result["status"] == "success"
    assert result["key_results"]["input_mode"] == "--tree-dir"
    assert not (out_dir / "wastral.tre").exists()
    assert not (out_dir / "merged.trees").exists()  # dry-run must not write files


def test_run_wastral_dry_run_tree_dir_multi_line(tmp_path: Path) -> None:
    from phyloai.tree.msc import run_wastral

    tree_dir = tmp_path / "genetrees"
    tree_dir.mkdir()
    # One file with 3 trees (multi-line), one file with 1 tree
    (tree_dir / "multi.trees").write_text("((a,b),c);\n((d,e),f);\n((g,h),i);\n")
    (tree_dir / "single.nwk").write_text("((x,y),z);\n")
    out_dir = tmp_path / "out"

    result = run_wastral(tree_dir=tree_dir, output_dir=out_dir, dry_run=True)

    assert result["status"] == "success"
    assert result["key_results"]["n_input_trees"] == 4  # 3+1 tree lines


def test_run_wastral_dry_run_outgroup(tmp_path: Path) -> None:
    from phyloai.tree.msc import run_wastral

    tree_file = tmp_path / "g.trees"
    tree_file.write_text("((a,b),c);\n")
    out_dir = tmp_path / "out"

    result = run_wastral(tree=tree_file, output_dir=out_dir, outgroup="Oryza",
                          dry_run=True)

    assert result["status"] == "success"
    assert result["params"]["outgroup"] == "Oryza"
    assert result["key_results"]["outgroup"] == "Oryza"
    assert "--root" in " ".join(result["data"]["cmd"])
    assert "Oryza" in result["data"]["cmd"]


def test_run_wastral_outgroup_empty_raises() -> None:
    from phyloai.tree.msc import run_wastral

    with pytest.raises(ValueError, match="outgroup must not be empty"):
        run_wastral(tree=Path("/nonexistent/file.trees"), outgroup="")


def test_run_wastral_dry_run_boot_3_includes_freq_quad(tmp_path: Path) -> None:
    from phyloai.tree.msc import run_wastral

    tree_file = tmp_path / "g.trees"
    tree_file.write_text("((a,b),c);\n")
    out_dir = tmp_path / "out"

    result = run_wastral(tree=tree_file, output_dir=out_dir, boot=3, dry_run=True)

    assert result["status"] == "success"
    assert "freq_quad_csv" in result["data"]
    assert "freqQuad.csv" in result["data"]["freq_quad_csv"]


def test_run_wastral_dry_run_boot_1_no_freq_quad(tmp_path: Path) -> None:
    from phyloai.tree.msc import run_wastral

    tree_file = tmp_path / "g.trees"
    tree_file.write_text("((a,b),c);\n")
    out_dir = tmp_path / "out"

    result = run_wastral(tree=tree_file, output_dir=out_dir, boot=1, dry_run=True)

    assert result["status"] == "success"
    assert "freq_quad_csv" not in result["data"]


def test_run_wastral_dry_run_tool_args_u3_includes_freq_quad(tmp_path: Path) -> None:
    from phyloai.tree.msc import run_wastral

    tree_file = tmp_path / "g.trees"
    tree_file.write_text("((a,b),c);\n")
    out_dir = tmp_path / "out"

    result = run_wastral(tree=tree_file, output_dir=out_dir, boot=1,
                          tool_args="-u 3", dry_run=True)

    assert result["status"] == "success"
    assert "freq_quad_csv" in result["data"]
    assert "freqQuad.csv" in result["data"]["freq_quad_csv"]
