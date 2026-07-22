"""Tests for phyloai.posttree.signal consistent subcommand."""
from __future__ import annotations

from pathlib import Path

import pytest


class TestRunSignalConsistentValidation:
    def test_more_than_two_trees_errors(self, tmp_path: Path) -> None:
        from phyloai.posttree.signal import run_signal_consistent

        matrix = tmp_path / "m.fa"
        matrix.write_text(">A\nMKT\n>B\nMKA\n")
        trees = [tmp_path / f"T{i}.nwk" for i in range(3)]
        for t in trees:
            t.write_text("(A,B);\n")
        tree_dir = tmp_path / "gtrees"
        tree_dir.mkdir()
        locus_ranges = tmp_path / "locus_ranges.txt"
        locus_ranges.write_text("DNA, locus1 = 1-6\n")
        result = run_signal_consistent(
            matrix=matrix,
            candidate_trees=trees,
            tree_dir=tree_dir,
            model_expr="LG+F+R4",
            partitions=None,
            partition_mode=None,
            locus_ranges=locus_ranges,
            guide_tree=None,
            threads="auto",
            iqtree_path=None,
            wastral_path=None,
            tool_args=None,
            metrics=None,
            output_dir=tmp_path / "out",
            overwrite=False,
            dry_run=True,
            quiet=True,
        )
        assert result["status"] == "error"
        assert "2" in result["error"]

    def test_no_locus_boundaries_errors(self, tmp_path: Path) -> None:
        from phyloai.posttree.signal import run_signal_consistent

        matrix = tmp_path / "m.fa"
        matrix.write_text(">A\nMKT\n>B\nMKA\n")
        trees = [tmp_path / f"T{i}.nwk" for i in range(2)]
        for t in trees:
            t.write_text("(A,B);\n")
        tree_dir = tmp_path / "gtrees"
        tree_dir.mkdir()
        result = run_signal_consistent(
            matrix=matrix,
            candidate_trees=trees,
            tree_dir=tree_dir,
            model_expr="LG+F+R4",
            partitions=None,
            partition_mode=None,
            locus_ranges=None,
            guide_tree=None,
            threads="auto",
            iqtree_path=None,
            wastral_path=None,
            tool_args=None,
            metrics=None,
            output_dir=tmp_path / "out",
            overwrite=False,
            dry_run=True,
            quiet=True,
        )
        assert result["status"] == "error"
        assert "locus" in result["error"].lower() or "partition" in result["error"].lower()

    def test_threads_zero_errors(self, tmp_path: Path) -> None:
        from phyloai.posttree.signal import run_signal_consistent

        matrix = tmp_path / "m.fa"
        matrix.write_text(">A\nMKT\n>B\nMKA\n")
        trees = [tmp_path / f"T{i}.nwk" for i in range(2)]
        for t in trees:
            t.write_text("(A,B);\n")
        tree_dir = tmp_path / "gtrees"
        tree_dir.mkdir()
        locus_ranges = tmp_path / "locus_ranges.txt"
        locus_ranges.write_text("DNA, locus1 = 1-6\n")
        result = run_signal_consistent(
            matrix=matrix,
            candidate_trees=trees,
            tree_dir=tree_dir,
            model_expr="LG+F+R4",
            partitions=None,
            partition_mode=None,
            locus_ranges=locus_ranges,
            guide_tree=None,
            threads="0",
            iqtree_path=None,
            wastral_path=None,
            tool_args=None,
            metrics=None,
            output_dir=tmp_path / "out",
            overwrite=False,
            dry_run=True,
            quiet=True,
        )
        assert result["status"] == "error"
        assert "threads" in result["error"].lower()

    def test_threads_negative_errors(self, tmp_path: Path) -> None:
        from phyloai.posttree.signal import run_signal_consistent

        matrix = tmp_path / "m.fa"
        matrix.write_text(">A\nMKT\n>B\nMKA\n")
        trees = [tmp_path / f"T{i}.nwk" for i in range(2)]
        for t in trees:
            t.write_text("(A,B);\n")
        tree_dir = tmp_path / "gtrees"
        tree_dir.mkdir()
        locus_ranges = tmp_path / "locus_ranges.txt"
        locus_ranges.write_text("DNA, locus1 = 1-6\n")
        result = run_signal_consistent(
            matrix=matrix,
            candidate_trees=trees,
            tree_dir=tree_dir,
            model_expr="LG+F+R4",
            partitions=None,
            partition_mode=None,
            locus_ranges=locus_ranges,
            guide_tree=None,
            threads="-1",
            iqtree_path=None,
            wastral_path=None,
            tool_args=None,
            metrics=None,
            output_dir=tmp_path / "out",
            overwrite=False,
            dry_run=True,
            quiet=True,
        )
        assert result["status"] == "error"
        assert "threads" in result["error"].lower()

    def test_threads_non_numeric_errors(self, tmp_path: Path) -> None:
        from phyloai.posttree.signal import run_signal_consistent

        matrix = tmp_path / "m.fa"
        matrix.write_text(">A\nMKT\n>B\nMKA\n")
        trees = [tmp_path / f"T{i}.nwk" for i in range(2)]
        for t in trees:
            t.write_text("(A,B);\n")
        tree_dir = tmp_path / "gtrees"
        tree_dir.mkdir()
        locus_ranges = tmp_path / "locus_ranges.txt"
        locus_ranges.write_text("DNA, locus1 = 1-6\n")
        result = run_signal_consistent(
            matrix=matrix,
            candidate_trees=trees,
            tree_dir=tree_dir,
            model_expr="LG+F+R4",
            partitions=None,
            partition_mode=None,
            locus_ranges=locus_ranges,
            guide_tree=None,
            threads="abc",
            iqtree_path=None,
            wastral_path=None,
            tool_args=None,
            metrics=None,
            output_dir=tmp_path / "out",
            overwrite=False,
            dry_run=True,
            quiet=True,
        )
        assert result["status"] == "error"
        assert "threads" in result["error"].lower()

    def test_extra_trees_in_tree_dir_are_ignored(self, tmp_path: Path) -> None:
        from phyloai.posttree.signal import run_signal_consistent

        matrix = tmp_path / "m.fa"
        matrix.write_text(">A\nMKT\n>B\nMKA\n")
        trees = [tmp_path / f"T{i}.nwk" for i in range(2)]
        for t in trees:
            t.write_text("(A,B);\n")
        tree_dir = tmp_path / "gtrees"
        tree_dir.mkdir()
        (tree_dir / "gene1.nwk").write_text("(A,B);\n")
        (tree_dir / "extra_gene.nwk").write_text("(A,B);\n")
        locus_ranges = tmp_path / "locus_ranges.txt"
        locus_ranges.write_text("DNA, gene1 = 1-6\n")

        result = run_signal_consistent(
            matrix=matrix, candidate_trees=trees, tree_dir=tree_dir,
            model_expr="LG+F+R4", partitions=None, partition_mode=None,
            locus_ranges=locus_ranges, guide_tree=None, threads="auto",
            iqtree_path=None, wastral_path=None, tool_args=None, metrics=None,
            output_dir=tmp_path / "out", overwrite=False, dry_run=True, quiet=True,
        )
        assert result["status"] == "success"

    def test_multi_suffix_gene_tree_map_expansion(self, tmp_path: Path) -> None:
        from phyloai.core.file_matching import logical_tree_locus_candidates, scan_tree_dir

        tree_dir = tmp_path / "gtrees"
        tree_dir.mkdir()
        (tree_dir / "gene1.fa.treefile").write_text("(A,B);\n")
        (tree_dir / "gene2.raxml.bestTree").write_text("(A,B);\n")

        gene_tree_map = scan_tree_dir(tree_dir)
        # scan_tree_dir indexes one-suffix reduction, not two-suffix
        assert "gene1.fa" in gene_tree_map
        assert "gene1" not in gene_tree_map
        assert "gene2.raxml" in gene_tree_map
        assert "gene2" not in gene_tree_map

        # Apply the expansion (same code as run_signal_consistent)
        for tree_path in list(gene_tree_map.values()):
            _, cand2 = logical_tree_locus_candidates(tree_path)
            if cand2 and cand2 not in gene_tree_map:
                gene_tree_map[cand2] = tree_path

        # Now partition loci "gene1", "gene2" resolve via GQS gene_tree_map lookup
        assert "gene1" in gene_tree_map
        assert gene_tree_map["gene1"] == gene_tree_map["gene1.fa"]
        assert "gene2" in gene_tree_map
        assert gene_tree_map["gene2"] == gene_tree_map["gene2.raxml"]


class TestRunSignalConsistentFixture:
    SIGNAL_DIR = Path("runs/signal")

    def test_consistent_with_fixture(self, tmp_path: Path) -> None:
        from phyloai.posttree.signal import run_signal_consistent, _parse_partition_ranges

        if not self.SIGNAL_DIR.exists():
            pytest.skip("Signal test data not present")
        matrix = self.SIGNAL_DIR / "matrix.aa.fa"
        t1 = self.SIGNAL_DIR / "T1.tre"
        t2 = self.SIGNAL_DIR / "T2.tre"
        full_tree_dir = self.SIGNAL_DIR / "gene_trees1066"
        partitions = self.SIGNAL_DIR / "matrix.aa.partitions"
        if not matrix.exists() or not full_tree_dir.exists():
            pytest.skip("Fixture files missing")

        # Partition file has 20 loci; full tree dir has 1066. Only copy the
        # 20 matching gene trees so locus<->tree matching passes.
        partition_recs = _parse_partition_ranges(partitions)
        partition_loci = {r["locus"] for r in partition_recs}
        test_tree_dir = tmp_path / "gene_trees"
        test_tree_dir.mkdir()
        import shutil as _shutil
        for entry in sorted(full_tree_dir.iterdir()):
            if not entry.is_file():
                continue
            from phyloai.core.file_matching import logical_tree_locus_candidates
            for cand in logical_tree_locus_candidates(entry):
                if cand is not None and cand in partition_loci:
                    _shutil.copy2(entry, test_tree_dir / entry.name)
                    break

        assert len(list(test_tree_dir.iterdir())) == 20, "Subset gene tree dir missing entries"

        result = run_signal_consistent(
            matrix=matrix,
            candidate_trees=[t1, t2],
            tree_dir=test_tree_dir,
            model_expr="LG+F+R4",
            partitions=None,
            partition_mode=None,
            locus_ranges=partitions,
            guide_tree=None,
            threads="auto",
            iqtree_path=None,
            wastral_path=None,
            tool_args=None,
            metrics=None,
            output_dir=tmp_path / "consistent_out",
            overwrite=False,
            dry_run=False,
            quiet=True,
        )
        if result["status"] == "error" and result.get("error_category") == "env":
            pytest.skip("iqtree3 or wastral not available")
        assert result["status"] == "success"
        assert result["key_results"]["n_loci"] == 20
        assert (tmp_path / "consistent_out" / "consistent_genes.txt").exists()
        assert (tmp_path / "consistent_out" / "inconsistent_genes.txt").exists()
