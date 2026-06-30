from __future__ import annotations

import json
import random
from pathlib import Path

import pytest


def test_parse_partitions_raxml_style(tmp_path: Path) -> None:
    from phyloai.pretree.concat import _parse_partitions

    path = tmp_path / "matrix.partitions"
    path.write_text("LG, geneA = 1-10\nLG, geneB = 11-25\n")

    parts = _parse_partitions(path)

    assert parts == [
        {"model": "LG", "locus": "geneA", "start": 1, "end": 10, "length": 10},
        {"model": "LG", "locus": "geneB", "start": 11, "end": 25, "length": 15},
    ]


def test_parse_partitions_rejects_bad_line(tmp_path: Path) -> None:
    from phyloai.pretree.concat import _parse_partitions

    path = tmp_path / "bad.partitions"
    path.write_text("not a partition\n")

    with pytest.raises(ValueError, match="Unparseable partition line 1"):
        _parse_partitions(path)


def test_sample_partition_replicate_without_replacement_reaches_target() -> None:
    from phyloai.pretree.concat import _sample_partition_replicate

    parts = [
        {"model": "LG", "locus": "a", "start": 1, "end": 10, "length": 10},
        {"model": "LG", "locus": "b", "start": 11, "end": 20, "length": 10},
        {"model": "LG", "locus": "c", "start": 21, "end": 30, "length": 10},
    ]

    sampled = _sample_partition_replicate(parts, 20, random.Random(42))

    assert sum(p["length"] for p in sampled) >= 20
    assert len({p["locus"] for p in sampled}) == len(sampled)


def test_run_concat_jackknife_writes_replicate_dirs_and_result_json(tmp_path: Path) -> None:
    from phyloai.pretree.concat import run_concat_jackknife

    matrix = tmp_path / "matrix.fa"
    matrix.write_text(">A\nAAAACCCCGGGGTTTT\n>B\nAAAACCCCGGGGTTTT\n")
    parts = tmp_path / "matrix.partitions"
    parts.write_text("LG, gene1 = 1-4\nLG, gene2 = 5-8\nLG, gene3 = 9-12\nLG, gene4 = 13-16\n")
    out = tmp_path / "jackknife"

    payload = run_concat_jackknife(
        matrix=matrix,
        partitions=parts,
        output_dir=out,
        replicates=2,
        target_length=8,
        prefix="rep",
        to="fasta",
        table_format="csv",
        seed=42,
        overwrite=False,
        dry_run=False,
        quiet=True,
    )

    assert payload["status"] == "success"
    assert (out / "rep001" / "rep001.fa").exists()
    assert (out / "rep001" / "rep001.partitions").exists()
    assert (out / "rep002" / "rep002.fa").exists()
    assert (out / "jackknife_summary.csv").exists()
    saved = json.loads((out / "result.json").read_text())
    assert saved["params"]["seed"] == 42
    assert "summary" in saved["data"]["output_files"]
    assert "rep001_matrix" in saved["data"]["output_files"]
    assert "rep001_partitions" in saved["data"]["output_files"]


def test_run_concat_jackknife_rewrites_partition_coordinates(tmp_path: Path) -> None:
    from phyloai.pretree.concat import run_concat_jackknife

    matrix = tmp_path / "matrix.fa"
    matrix.write_text(">A\nAAAACCCCGGGG\n>B\nAAAACCCCGGGG\n")
    parts = tmp_path / "matrix.partitions"
    parts.write_text("LG, gene1 = 1-4\nLG, gene2 = 5-8\nLG, gene3 = 9-12\n")

    run_concat_jackknife(
        matrix=matrix,
        partitions=parts,
        output_dir=tmp_path / "out",
        replicates=1,
        target_length=8,
        seed=42,
        quiet=True,
    )

    lines = (tmp_path / "out" / "rep001" / "rep001.partitions").read_text().splitlines()
    assert lines[0].endswith("= 1-4")
    assert lines[1].endswith("= 5-8")


def test_run_concat_jackknife_result_json_schema(tmp_path: Path) -> None:
    from phyloai.pretree.concat import run_concat_jackknife
    from tests.helpers import validate_params_completeness, validate_result_json

    matrix = tmp_path / "matrix.fa"
    matrix.write_text(">A\nAAAACCCC\n>B\nAAAACCCC\n")
    parts = tmp_path / "matrix.partitions"
    parts.write_text("LG, gene1 = 1-4\nLG, gene2 = 5-8\n")
    payload = run_concat_jackknife(matrix, parts, tmp_path / "out", replicates=1, target_length=4, quiet=True)

    validate_result_json(payload)
    validate_params_completeness(payload, {
        "matrix", "partitions", "output_dir", "replicates", "target_length",
        "prefix", "to", "table_format", "seed", "overwrite", "dry_run", "quiet",
    })
