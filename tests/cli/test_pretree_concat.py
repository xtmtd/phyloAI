from __future__ import annotations

import json
from pathlib import Path

from click.testing import CliRunner

from phyloai.cli.main import cli


def test_cli_pretree_concat_basic(tmp_path: Path) -> None:
    msa_dir = tmp_path / "msas"
    msa_dir.mkdir()
    (msa_dir / "gene1.fa").write_text(">A\nACGT\n>B\nACGT\n>C\nACGT\n")
    (msa_dir / "gene2.fa").write_text(">A\nGGCC\n>B\nGGCC\n>C\nGGCC\n")

    output_dir = tmp_path / "out"
    result = CliRunner().invoke(
        cli,
        [
            "pretree", "concat",
            "--msa-dir", str(msa_dir),
            "--output-dir", str(output_dir),
            "--seq-type", "NT",
            "--to", "fasta",
        ],
    )

    assert result.exit_code == 0, result.output
    result_path = output_dir / "result.json"
    assert result_path.exists()
    payload = json.loads(result_path.read_text())
    assert payload["status"] == "success"
    assert payload["key_results"]["n_taxa"] == 3
