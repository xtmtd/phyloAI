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
            "original.gaps.fa"
        )
        assert payload["tool_versions"] == {}
        assert "--original-msa" in payload["command"]
        assert "--simulated-msa" in payload["command"]
        assert "--seq-type auto" in payload["command"]
        assert "-o" in payload["command"]

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


class TestBatchSimulatedDir:
    def test_dir_mode_writes_one_gaps_file_per_simulated(
        self, tmp_path: Path, original: Path,
    ) -> None:
        sim_dir = tmp_path / "sims"
        sim_dir.mkdir()
        _fasta(sim_dir / "sim001.fa", [("A", "ACGTGT"), ("B", "ACGTAC")])
        _fasta(sim_dir / "sim002.fa", [("A", "TTTTTT"), ("B", "CCCCCC")])
        result = run_alisim_transfergaps(
            original_msa=original, simulated_dir=sim_dir,
            output_dir=tmp_path / "out",
        )
        out_files = result["data"]["output_files"]["transferred_msas"]
        assert set(out_files) == {"sim001", "sim002"}
        assert Path(out_files["sim001"]["path"]).name == "sim001.gaps.fa"
        assert Path(out_files["sim002"]["path"]).name == "sim002.gaps.fa"
        from Bio import SeqIO

        rec = next(SeqIO.parse(str(out_files["sim001"]["path"]), "fasta"))
        assert str(rec.seq) == "AC-TG-"
        assert result["key_results"]["n_msas"] == 2
        assert result["key_results"]["n_positions_masked"] == 6
        assert result["key_results"]["detected_seq_type"] == "NT"

    def test_dir_mode_only_considers_alignment_extensions(
        self, tmp_path: Path, original: Path,
    ) -> None:
        sim_dir = tmp_path / "sims"
        sim_dir.mkdir()
        _fasta(sim_dir / "sim001.fa", [("A", "ACGTGT"), ("B", "ACGTAC")])
        (sim_dir / "notes.txt").write_text("not an alignment")
        result = run_alisim_transfergaps(
            original_msa=original, simulated_dir=sim_dir,
            output_dir=tmp_path / "out",
        )
        out_files = result["data"]["output_files"]["transferred_msas"]
        assert set(out_files) == {"sim001"}

    def test_dir_mode_mismatch_names_failing_file(
        self, tmp_path: Path, original: Path,
    ) -> None:
        sim_dir = tmp_path / "sims"
        sim_dir.mkdir()
        _fasta(sim_dir / "sim001.fa", [("A", "ACGTGT"), ("B", "ACGTAC")])
        _fasta(sim_dir / "sim002.fa", [("A", "TTTTTT"), ("C", "CCCCCC")])
        with pytest.raises(ValueError, match=r"sim002\.fa: taxon name mismatch"):
            run_alisim_transfergaps(
                original_msa=original, simulated_dir=sim_dir,
                output_dir=tmp_path / "out",
            )

    def test_dir_mode_dry_run_writes_nothing(
        self, tmp_path: Path, original: Path,
    ) -> None:
        sim_dir = tmp_path / "sims"
        sim_dir.mkdir()
        _fasta(sim_dir / "sim001.fa", [("A", "ACGTGT"), ("B", "ACGTAC")])
        output_dir = tmp_path / "out"
        result = run_alisim_transfergaps(
            original_msa=original, simulated_dir=sim_dir,
            output_dir=output_dir, dry_run=True,
        )
        assert result["status"] == "success"
        assert result["key_results"]["n_msas"] == 1
        assert not output_dir.exists()

    def test_missing_dir_rejected(self, tmp_path: Path, original: Path) -> None:
        with pytest.raises(ValueError, match="--simulated-dir does not exist"):
            run_alisim_transfergaps(
                original_msa=original, simulated_dir=tmp_path / "nope",
                output_dir=tmp_path / "out",
            )

    def test_empty_dir_rejected(self, tmp_path: Path, original: Path) -> None:
        sim_dir = tmp_path / "sims"
        sim_dir.mkdir()
        with pytest.raises(ValueError, match="no alignment files"):
            run_alisim_transfergaps(
                original_msa=original, simulated_dir=sim_dir,
                output_dir=tmp_path / "out",
            )

    def test_neither_input_rejected(
        self, tmp_path: Path, original: Path,
    ) -> None:
        with pytest.raises(ValueError, match="exactly one of"):
            run_alisim_transfergaps(
                original_msa=original, output_dir=tmp_path / "out",
            )

    def test_both_inputs_rejected(
        self, tmp_path: Path, original: Path, simulated: Path,
    ) -> None:
        sim_dir = tmp_path / "sims"
        sim_dir.mkdir()
        with pytest.raises(ValueError, match="exactly one of"):
            run_alisim_transfergaps(
                original_msa=original, simulated_msa=simulated,
                simulated_dir=sim_dir, output_dir=tmp_path / "out",
            )

    def test_dir_mode_result_json_written(
        self, tmp_path: Path, original: Path,
    ) -> None:
        sim_dir = tmp_path / "sims"
        sim_dir.mkdir()
        _fasta(sim_dir / "sim001.fa", [("A", "ACGTGT"), ("B", "ACGTAC")])
        output_dir = tmp_path / "out"
        run_alisim_transfergaps(
            original_msa=original, simulated_dir=sim_dir,
            output_dir=output_dir,
        )
        payload = json.loads((output_dir / "result.json").read_text())
        assert payload["status"] == "success"
        assert payload["params"]["simulated_dir"].endswith("sims")
        assert payload["params"]["simulated_msa"] is None
        assert payload["key_results"]["n_msas"] == 1


class TestInputFormats:
    """Non-FASTA inputs are read via the shared FormatConverter; output is
    always FASTA with a .gaps.fa name."""

    def test_single_mode_accepts_phylip_and_writes_gaps_fa(
        self, tmp_path: Path,
    ) -> None:
        original_phy = tmp_path / "original.phy"
        original_phy.write_text("2 6\nA      AC-GT-\nB      ACG-TA\n")
        sim_phy = tmp_path / "sim.phy"
        sim_phy.write_text("2 6\nA      ACGTGT\nB      ACGTAC\n")
        result = run_alisim_transfergaps(
            original_msa=original_phy, simulated_msa=sim_phy,
            output_dir=tmp_path / "out",
        )
        out = Path(result["data"]["output_files"]["transferred_msa"]["path"])
        assert out.name == "original.gaps.fa"
        assert out.read_text().startswith(">A\n")
        assert ">B" in out.read_text()

    def test_single_mode_accepts_nexus(self, tmp_path: Path) -> None:
        original_nex = tmp_path / "original.nex"
        original_nex.write_text(
            "#NEXUS\nBEGIN DATA;\nDIMENSIONS NTAX=2 NCHAR=6;\n"
            "FORMAT DATATYPE=PROTEIN MISSING=? GAP=-;\nMATRIX\n"
            "A AC-GT-\nB ACG-TA\n;\nEND;\n"
        )
        sim_nex = tmp_path / "sim.nex"
        sim_nex.write_text(
            "#NEXUS\nBEGIN DATA;\nDIMENSIONS NTAX=2 NCHAR=6;\n"
            "FORMAT DATATYPE=PROTEIN MISSING=? GAP=-;\nMATRIX\n"
            "A ACGTGT\nB ACGTAC\n;\nEND;\n"
        )
        result = run_alisim_transfergaps(
            original_msa=original_nex, simulated_msa=sim_nex,
            output_dir=tmp_path / "out",
        )
        out = Path(result["data"]["output_files"]["transferred_msa"]["path"])
        assert out.name == "original.gaps.fa"

    def test_single_mode_accepts_phylip_paml(self, tmp_path: Path) -> None:
        original_paml = tmp_path / "original.paml.phy"
        original_paml.write_text("2 6 S\nA                             AC-GT-\nB                             ACG-TA\n")
        sim_paml = tmp_path / "sim.paml.phy"
        sim_paml.write_text("2 6 S\nA                             ACGTGT\nB                             ACGTAC\n")
        result = run_alisim_transfergaps(
            original_msa=original_paml, simulated_msa=sim_paml,
            output_dir=tmp_path / "out",
        )
        out = Path(result["data"]["output_files"]["transferred_msa"]["path"])
        assert out.name.endswith(".gaps.fa")
        assert ">A" in out.read_text()

    def test_batch_mixed_formats_all_write_gaps_fa(
        self, tmp_path: Path, original: Path,
    ) -> None:
        sim_dir = tmp_path / "sims"
        sim_dir.mkdir()
        _fasta(sim_dir / "sim001.fa", [("A", "ACGTGT"), ("B", "ACGTAC")])
        (sim_dir / "sim002.phy").write_text(
            "2 6\nA      ACGTGT\nB      ACGTAC\n"
        )
        (sim_dir / "sim003.nex").write_text(
            "#NEXUS\nBEGIN DATA;\nDIMENSIONS NTAX=2 NCHAR=6;\n"
            "FORMAT DATATYPE=PROTEIN MISSING=? GAP=-;\nMATRIX\n"
            "A ACGTGT\nB ACGTAC\n;\nEND;\n"
        )
        result = run_alisim_transfergaps(
            original_msa=original, simulated_dir=sim_dir,
            output_dir=tmp_path / "out",
        )
        out_files = result["data"]["output_files"]["transferred_msas"]
        assert set(out_files) == {"sim001", "sim002", "sim003"}
        names = {Path(info["path"]).name for info in out_files.values()}
        assert names == {"sim001.gaps.fa", "sim002.gaps.fa", "sim003.gaps.fa"}
        assert result["key_results"]["n_msas"] == 3
        for info in out_files.values():
            assert Path(info["path"]).read_text().startswith(">A\n")

    def test_unparsable_file_rejected(self, tmp_path: Path, original: Path) -> None:
        sim_dir = tmp_path / "sims"
        sim_dir.mkdir()
        _fasta(sim_dir / "sim001.fa", [("A", "ACGTGT"), ("B", "ACGTAC")])
        (sim_dir / "sim002.fa").write_text("this is not an alignment\n")
        with pytest.raises(ValueError, match="unable to parse"):
            run_alisim_transfergaps(
                original_msa=original, simulated_dir=sim_dir,
                output_dir=tmp_path / "out",
            )
