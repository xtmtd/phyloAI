"""Tests for phyloai.posttree.simulate_alisim_transfergaps."""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from phyloai.posttree.simulate_alisim_transfergaps import run_alisim_transfergaps


def _fasta(path: Path, records: list[tuple[str, str]]) -> None:
    path.write_text("".join(f">{name}\n{seq}\n" for name, seq in records))


@pytest.fixture
def original(tmp_path: Path) -> Path:
    path = tmp_path / "original.fa"
    _fasta(path, [("A", "AC-GT-"), ("B", "ACG-TA")])
    return path


@pytest.fixture
def simulated(tmp_path: Path) -> Path:
    path = tmp_path / "sim001.fa"
    _fasta(path, [("A", "ACGTGT"), ("B", "ACGTAC")])
    return path


class TestRunAlisimTransfergaps:
    def test_transfer_matches_taxa_by_name_and_preserves_original_order(
        self, tmp_path: Path, original: Path, simulated: Path,
    ) -> None:
        reordered = tmp_path / "reordered.fa"
        _fasta(reordered, [("B", "ACGTAC"), ("A", "ACGTGT")])
        result = run_alisim_transfergaps(
            original_msa=original, simulated_msa=reordered,
            output_dir=tmp_path / "out",
        )
        output = Path(result["data"]["output_files"]["transferred_msa"]["path"])
        from Bio import SeqIO

        records = list(SeqIO.parse(str(output), "fasta"))
        assert [record.id for record in records] == ["A", "B"]
        assert str(records[0].seq) == "AC-TG-"
        assert str(records[1].seq) == "ACG-AC"
        assert result["status"] == "success"

    def test_exclude_ambiguity_masks_only_dash_and_dot(
        self, tmp_path: Path,
    ) -> None:
        original = tmp_path / "orig_with_x.fa"
        _fasta(original, [("A", "AC-G-X")])
        simulated = tmp_path / "sim.fa"
        _fasta(simulated, [("A", "ACGTGT")])
        result = run_alisim_transfergaps(
            original_msa=original, simulated_msa=simulated,
            exclude_ambiguity=True, output_dir=tmp_path / "out",
        )
        output = Path(result["data"]["output_files"]["transferred_msa"]["path"])
        from Bio import SeqIO

        record = next(SeqIO.parse(str(output), "fasta"))
        assert str(record.seq) == "AC-T-T"  # X (ambiguity) kept, '-' transferred
        assert result["key_results"]["n_positions_masked"] == 2

    def test_default_masks_ambiguity_codes(
        self, tmp_path: Path,
    ) -> None:
        original = tmp_path / "orig_with_x.fa"
        _fasta(original, [("A", "AC-G-X")])
        simulated = tmp_path / "sim.fa"
        _fasta(simulated, [("A", "ACGTGT")])
        result = run_alisim_transfergaps(
            original_msa=original, simulated_msa=simulated,
            output_dir=tmp_path / "out",
        )
        output = Path(result["data"]["output_files"]["transferred_msa"]["path"])
        from Bio import SeqIO

        record = next(SeqIO.parse(str(output), "fasta"))
        assert str(record.seq) == "AC-T--"  # both '-' and 'X' masked
        assert result["key_results"]["n_positions_masked"] == 3

    def test_dry_run_validates_but_writes_nothing(
        self, tmp_path: Path, original: Path, simulated: Path,
    ) -> None:
        output_dir = tmp_path / "out"
        result = run_alisim_transfergaps(
            original_msa=original, simulated_msa=simulated,
            output_dir=output_dir, dry_run=True,
        )
        assert result["status"] == "success"
        assert result["key_results"]["n_sequences"] == 2
        assert not output_dir.exists()

    def test_auto_detects_aa_when_unambiguous(
        self, tmp_path: Path, simulated: Path,
    ) -> None:
        original = tmp_path / "aa.fa"
        _fasta(original, [("A", "MKWY")])
        sim = tmp_path / "sim.fa"
        _fasta(sim, [("A", "MKAP")])
        result = run_alisim_transfergaps(
            original_msa=original, simulated_msa=sim,
            output_dir=tmp_path / "out", dry_run=True,
        )
        assert result["key_results"]["detected_seq_type"] == "AA"

    def test_taxon_set_mismatch_rejected(
        self, tmp_path: Path, original: Path,
    ) -> None:
        sim = tmp_path / "sim.fa"
        _fasta(sim, [("A", "ACGTGT"), ("C", "ACGTAC")])
        with pytest.raises(ValueError, match="taxon name mismatch"):
            run_alisim_transfergaps(
                original_msa=original, simulated_msa=sim,
                output_dir=tmp_path / "out",
            )

    def test_length_mismatch_rejected(
        self, tmp_path: Path, original: Path,
    ) -> None:
        sim = tmp_path / "sim.fa"
        _fasta(sim, [("A", "ACGT"), ("B", "ACGT")])
        with pytest.raises(ValueError, match="length mismatch"):
            run_alisim_transfergaps(
                original_msa=original, simulated_msa=sim,
                output_dir=tmp_path / "out",
            )

    def test_duplicate_taxon_rejected(
        self, tmp_path: Path, simulated: Path,
    ) -> None:
        orig = tmp_path / "orig.fa"
        _fasta(orig, [("A", "AC-GT-"), ("A", "ACG-TA")])
        with pytest.raises(ValueError, match="duplicate taxon"):
            run_alisim_transfergaps(
                original_msa=orig, simulated_msa=simulated,
                output_dir=tmp_path / "out",
            )

    def test_missing_file_rejected(self, tmp_path: Path) -> None:
        with pytest.raises(ValueError, match="does not exist"):
            run_alisim_transfergaps(
                original_msa=tmp_path / "nope.fa",
                simulated_msa=tmp_path / "nope2.fa",
                output_dir=tmp_path / "out",
            )

    def test_empty_file_rejected(self, tmp_path: Path, simulated: Path) -> None:
        orig = tmp_path / "orig.fa"
        orig.write_text("")
        with pytest.raises(ValueError, match="empty"):
            run_alisim_transfergaps(
                original_msa=orig, simulated_msa=simulated,
                output_dir=tmp_path / "out",
            )

    def test_60_column_output(
        self, tmp_path: Path, simulated: Path,
    ) -> None:
        long_seq = "M" * 130
        original = tmp_path / "long.fa"
        _fasta(original, [("A", long_seq)])
        sim = tmp_path / "sim.fa"
        _fasta(sim, [("A", "A" * 130)])
        result = run_alisim_transfergaps(
            original_msa=original, simulated_msa=sim,
            output_dir=tmp_path / "out",
        )
        output = Path(result["data"]["output_files"]["transferred_msa"]["path"])
        lines = output.read_text().strip().splitlines()
        assert len(lines) == 1 + 3  # header + 130/60 wrapped
        assert max(len(line) for line in lines[1:]) == 60

    def test_result_json_written(
        self, tmp_path: Path, original: Path, simulated: Path,
    ) -> None:
        output_dir = tmp_path / "out"
        run_alisim_transfergaps(
            original_msa=original, simulated_msa=simulated,
            output_dir=output_dir,
        )
        payload = json.loads((output_dir / "result.json").read_text())
        assert payload["status"] == "success"
        assert payload["key_results"]["alignment_length"] == 6
        assert payload["data"]["output_files"]["transferred_msa"]["path"].endswith(
            "original_transferred.fa"
        )
        assert payload["tool_versions"] == {}

    def test_output_conflict_rejected(
        self, tmp_path: Path, original: Path, simulated: Path,
    ) -> None:
        output_dir = tmp_path / "out"
        output_dir.mkdir()
        (output_dir / "old.txt").write_text("x")
        with pytest.raises(ValueError, match="non-empty"):
            run_alisim_transfergaps(
                original_msa=original, simulated_msa=simulated,
                output_dir=output_dir,
            )

    def test_seq_type_forced_nt(
        self, tmp_path: Path,
    ) -> None:
        original = tmp_path / "nt.fa"
        _fasta(original, [("A", "AC-N-")])
        sim = tmp_path / "sim.fa"
        _fasta(sim, [("A", "ACGTA")])
        result = run_alisim_transfergaps(
            original_msa=original, simulated_msa=sim,
            seq_type="NT", output_dir=tmp_path / "out",
        )
        output = Path(result["data"]["output_files"]["transferred_msa"]["path"])
        from Bio import SeqIO

        record = next(SeqIO.parse(str(output), "fasta"))
        assert str(record.seq) == "AC---"  # '-', 'N' and '-' masked (NT valid = ACGT)
        assert result["key_results"]["detected_seq_type"] == "NT"

    def test_invalid_seq_type_rejected(
        self, tmp_path: Path, original: Path, simulated: Path,
    ) -> None:
        with pytest.raises(ValueError, match="--seq-type"):
            run_alisim_transfergaps(
                original_msa=original, simulated_msa=simulated,
                seq_type="RNA", output_dir=tmp_path / "out",
            )
