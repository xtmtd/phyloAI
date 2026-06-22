"""Tests for symtest filtering in filter.py."""

import io

import pytest

from phyloai.pretree.filter import _parse_symtest_csv, _build_symtest_supermatrix


# --- _parse_symtest_csv ---

def test_parse_symtest_csv_all_columns():
    csv_content = (
        "Name,SymSig,SymNon,SymPval,MarSig,MarNon,MarPval,IntSig,IntNon,IntPval\n"
        "gene1,44,92,0.475,50,86,0.722,4,132,0.239\n"
        "gene2,43,93,0.142,49,87,0.205,5,131,0.170\n"
        "gene3,53,83,0.005,58,78,0.002,6,130,0.343\n"
    )
    fp = io.StringIO(csv_content)
    results = _parse_symtest_csv(fp)
    assert len(results) == 3
    assert results[0]["Name"] == "gene1"
    assert results[0]["SymPval"] == 0.475
    assert results[2]["SymPval"] == 0.005
    assert results[1]["MarPval"] == 0.205
    assert results[1]["IntPval"] == 0.170


def test_parse_symtest_csv_skips_comment_lines():
    csv_content = (
        "# Matched-pair tests of symmetry\n"
        "# comment\n"
        "Name,SymSig,SymNon,SymPval,MarSig,MarNon,MarPval,IntSig,IntNon,IntPval\n"
        "gene1,44,92,0.475,50,86,0.722,4,132,0.239\n"
    )
    fp = io.StringIO(csv_content)
    results = _parse_symtest_csv(fp)
    assert len(results) == 1
    assert results[0]["Name"] == "gene1"


def test_parse_symtest_csv_empty():
    csv_content = (
        "# comment only\n"
        "Name,SymSig,SymNon,SymPval,MarSig,MarNon,MarPval,IntSig,IntNon,IntPval\n"
    )
    fp = io.StringIO(csv_content)
    results = _parse_symtest_csv(fp)
    assert results == []


def test_parse_symtest_csv_missing_header_raises():
    csv_content = "bad,header,here\nval1,val2,val3\n"
    fp = io.StringIO(csv_content)
    with pytest.raises(ValueError, match="missing expected columns"):
        _parse_symtest_csv(fp)


def test_parse_symtest_csv_non_numeric_pval():
    csv_content = (
        "Name,SymSig,SymNon,SymPval,MarSig,MarNon,MarPval,IntSig,IntNon,IntPval\n"
        "gene1,44,92,NA,50,86,NA,4,132,NA\n"
    )
    fp = io.StringIO(csv_content)
    results = _parse_symtest_csv(fp)
    assert results[0]["SymPval"] is None


# --- _filter_by_symtest_pval ---

def test_filter_symtest_retain_above_threshold():
    from phyloai.pretree.filter import _filter_by_symtest_pval
    results = [
        {"Name": "gene1", "SymPval": 0.475, "SymSig": 44, "SymNon": 92,
         "MarSig": 50, "MarNon": 86, "MarPval": 0.722,
         "IntSig": 4, "IntNon": 132, "IntPval": 0.239},
        {"Name": "gene2", "SymPval": 0.005, "SymSig": 53, "SymNon": 83,
         "MarSig": 58, "MarNon": 78, "MarPval": 0.002,
         "IntSig": 6, "IntNon": 130, "IntPval": 0.343},
    ]
    retained, dropped, decisions = _filter_by_symtest_pval(results, "Sym", 0.05)
    assert len(retained) == 1
    assert retained[0]["locus"] == "gene1"
    assert len(dropped) == 1
    assert dropped[0]["locus"] == "gene2"
    assert decisions[0]["status"] == "retained"
    assert decisions[1]["status"] == "dropped"


def test_filter_symtest_uses_mar_column():
    from phyloai.pretree.filter import _filter_by_symtest_pval
    results = [
        {"Name": "gene1", "SymPval": 0.001, "SymSig": 44, "SymNon": 92,
         "MarSig": 50, "MarNon": 86, "MarPval": 0.722,
         "IntSig": 4, "IntNon": 132, "IntPval": 0.239},
    ]
    retained, dropped, decisions = _filter_by_symtest_pval(results, "MAR", 0.05)
    assert len(retained) == 1  # MarPval=0.722 >= 0.05
    assert retained[0]["p_value"] == 0.722


def test_filter_symtest_none_pval_dropped():
    from phyloai.pretree.filter import _filter_by_symtest_pval
    results = [
        {"Name": "gene1", "SymPval": None, "SymSig": 44, "SymNon": 92,
         "MarSig": 50, "MarNon": 86, "MarPval": None,
         "IntSig": 4, "IntNon": 132, "IntPval": None},
    ]
    retained, dropped, decisions = _filter_by_symtest_pval(results, "Sym", 0.05)
    assert len(retained) == 0
    assert dropped[0]["reason"] == "p_value is null"


def test_filter_symtest_uses_int_column():
    from phyloai.pretree.filter import _filter_by_symtest_pval
    results = [
        {"Name": "gene1", "SymPval": 0.001, "SymSig": 44, "SymNon": 92,
         "MarSig": 50, "MarNon": 86, "MarPval": 0.002,
         "IntSig": 4, "IntNon": 132, "IntPval": 0.343},
    ]
    retained, dropped, decisions = _filter_by_symtest_pval(results, "INT", 0.05)
    assert len(retained) == 1  # IntPval=0.343 >= 0.05
    assert retained[0]["p_value"] == 0.343


def test_filter_symtest_boundary_at_threshold():
    from phyloai.pretree.filter import _filter_by_symtest_pval
    results = [
        {"Name": "gene1", "SymPval": 0.05, "SymSig": 44, "SymNon": 92,
         "MarSig": 50, "MarNon": 86, "MarPval": 0.05,
         "IntSig": 4, "IntNon": 132, "IntPval": 0.05},
    ]
    retained, dropped, decisions = _filter_by_symtest_pval(results, "Sym", 0.05)
    assert len(retained) == 1  # p == threshold is retained (p >= 0.05)
    assert retained[0]["p_value"] == 0.05


def test_filter_symtest_decision_has_all_pval_columns():
    from phyloai.pretree.filter import _filter_by_symtest_pval
    results = [
        {"Name": "gene1", "SymPval": 0.475, "SymSig": 44, "SymNon": 92,
         "MarSig": 50, "MarNon": 86, "MarPval": 0.722,
         "IntSig": 4, "IntNon": 132, "IntPval": 0.239},
    ]
    _, _, decisions = _filter_by_symtest_pval(results, "Sym", 0.05)
    d = decisions[0]
    assert d["sym_pval"] == 0.475
    assert d["mar_pval"] == 0.722
    assert d["int_pval"] == 0.239
    assert d["p_value"] == 0.475  # selected column


# --- _build_symtest_supermatrix ---

def test_build_symtest_supermatrix(tmp_path):
    msa1 = tmp_path / "gene1.fa"
    msa1.write_text(">taxa1\nACGT\n>taxa2\nACGT\n")
    msa2 = tmp_path / "gene2.fa"
    msa2.write_text(">taxa1\nTTTT\n>taxa2\nGGGG\n>taxa3\nCCCC\n")
    msa_map = {"gene1": msa1, "gene2": msa2}

    matrix_str, genes, prefix_type = _build_symtest_supermatrix(msa_map)

    assert prefix_type in ("DNA", "LG")
    assert len(genes) == 2
    assert genes[0][0] == "gene1"
    assert genes[0][1] == 1
    assert genes[0][2] == 4
    assert genes[1][0] == "gene2"
    assert genes[1][1] == 5
    assert genes[1][2] == 8
    assert ">taxa1" in matrix_str
    assert ">taxa2" in matrix_str
    assert ">taxa3" in matrix_str


def test_build_symtest_supermatrix_empty_dir():
    with pytest.raises(ValueError, match="No valid MSA files"):
        _build_symtest_supermatrix({})


# --- run_symtest integration with mock IQ-TREE ---

def test_run_symtest_with_mock_iqtree(tmp_path, monkeypatch):
    """Full integration: build supermatrix, invoke mock iqtree, parse output, filter."""
    from phyloai.pretree.filter import run_symtest

    msa_dir = tmp_path / "msa"
    msa_dir.mkdir()
    (msa_dir / "gene1.fa").write_text(">t1\nACGT\n>t2\nACGT\n")
    (msa_dir / "gene2.fa").write_text(">t1\nTTTT\n>t2\nGGGG\n")
    (msa_dir / "gene3.fa").write_text(">t1\nCCCC\n>t2\nAAAA\n")

    mock_iqtree = tmp_path / "mock_iqtree"
    mock_script = (
        '#!/usr/bin/env bash\n'
        '# Find the partitions file argument (-p <path>)\n'
        'partfile=""\n'
        'while [ $# -gt 0 ]; do\n'
        '  if [ "$1" = "-p" ]; then shift; partfile="$1"; fi\n'
        '  shift\n'
        'done\n'
        'symcsv="${partfile}.symtest.csv"\n'
        "cat > \"$symcsv\" <<'CSVEOF'\n"
        '# comment\n'
        'Name,SymSig,SymNon,SymPval,MarSig,MarNon,MarPval,IntSig,IntNon,IntPval\n'
        'gene1,44,92,0.475,50,86,0.722,4,132,0.239\n'
        'gene2,43,93,0.004,49,87,0.003,5,131,0.170\n'
        'gene3,53,83,0.620,58,78,0.550,6,130,0.343\n'
        'CSVEOF\n'
        'exit 0\n'
    )
    mock_iqtree.write_text(mock_script)
    mock_iqtree.chmod(0o755)
    monkeypatch.setattr("shutil.which", lambda x, path=None: str(mock_iqtree))

    output_dir = tmp_path / "out"
    payload = run_symtest(
        msa_dir=msa_dir, output_dir=output_dir,
        symtest_type=None, symtest_pval=0.05,
        iqtree_path=mock_iqtree, threads=1,
        table_format="csv",
    )

    assert payload["status"] == "success"
    assert payload["key_results"]["n_input"] == 3
    assert payload["key_results"]["n_retained"] == 2
    assert payload["key_results"]["n_dropped"] == 1

    assert (output_dir / "result.json").exists()
    assert isinstance(payload["data"]["cmd"], list)
    assert isinstance(payload["data"]["tool_stderr"], str)
    assert len(payload["data"]["results"]) == 3
    assert (output_dir / "retained_loci.csv").exists()
    assert (output_dir / "dropped_loci.csv").exists()
    assert (output_dir / "filter_decisions.csv").exists()

    seqs = output_dir / "seqs"
    assert (seqs / "gene1.fa").exists()
    assert (seqs / "gene3.fa").exists()
    assert not (seqs / "gene2.fa").exists()

    import csv
    with open(output_dir / "filter_decisions.csv") as fh:
        rows = list(csv.DictReader(fh))
    assert len(rows) == 3
    statuses = {r["locus"]: r["status"] for r in rows}
    assert statuses["gene1"] == "retained"
    assert statuses["gene2"] == "dropped"
    assert statuses["gene3"] == "retained"


def test_run_symtest_with_tree_dir(tmp_path, monkeypatch):
    """Tree copy: tree dir provided, only retained-locus trees copied."""
    from phyloai.pretree.filter import run_symtest

    msa_dir = tmp_path / "msa"
    msa_dir.mkdir()
    (msa_dir / "gene1.fa").write_text(">t1\nACGT\n>t2\nACGT\n")
    (msa_dir / "gene2.fa").write_text(">t1\nTTTT\n>t2\nGGGG\n")

    tree_dir = tmp_path / "trees"
    tree_dir.mkdir()
    (tree_dir / "gene1.tre").write_text("(t1,t2);")
    (tree_dir / "gene2.tre").write_text("(t1,t2);")

    mock_iqtree = tmp_path / "mock_iqtree"
    mock_script = (
        '#!/usr/bin/env bash\n'
        'partfile=""\n'
        'while [ $# -gt 0 ]; do\n'
        '  if [ "$1" = "-p" ]; then shift; partfile="$1"; fi\n'
        '  shift\n'
        'done\n'
        'symcsv="${partfile}.symtest.csv"\n'
        "cat > \"$symcsv\" <<'CSVEOF'\n"
        'Name,SymSig,SymNon,SymPval,MarSig,MarNon,MarPval,IntSig,IntNon,IntPval\n'
        'gene1,44,92,0.475,50,86,0.722,4,132,0.239\n'
        'gene2,43,93,0.004,49,87,0.003,5,131,0.170\n'
        'CSVEOF\n'
        'exit 0\n'
    )
    mock_iqtree.write_text(mock_script)
    mock_iqtree.chmod(0o755)
    monkeypatch.setattr("shutil.which", lambda x, path=None: str(mock_iqtree))

    output_dir = tmp_path / "out"
    payload = run_symtest(
        msa_dir=msa_dir, output_dir=output_dir,
        symtest_type=None, symtest_pval=0.05,
        iqtree_path=mock_iqtree, threads=1,
        tree_dir=tree_dir, table_format="csv",
    )

    assert payload["status"] == "success"
    assert payload["key_results"]["retained_trees_copied"] == 1
    trees = output_dir / "trees"
    assert (trees / "gene1.tre").exists()
    assert not (trees / "gene2.tre").exists()


def test_run_symtest_iqtree_nonzero_exit(tmp_path, monkeypatch):
    """IQ-TREE non-zero exit code raises RuntimeError."""
    from phyloai.pretree.filter import run_symtest

    msa_dir = tmp_path / "msa"
    msa_dir.mkdir()
    (msa_dir / "gene1.fa").write_text(">t1\nACGT\n>t2\nACGT\n")

    mock_iqtree = tmp_path / "mock_iqtree"
    mock_iqtree.write_text("#!/usr/bin/env bash\necho 'SIMULATED ERROR' >&2\nexit 1\n")
    mock_iqtree.chmod(0o755)
    monkeypatch.setattr("shutil.which", lambda x, path=None: str(mock_iqtree))

    output_dir = tmp_path / "out"
    with pytest.raises(RuntimeError, match="exited with code 1"):
        run_symtest(
            msa_dir=msa_dir, output_dir=output_dir,
            iqtree_path=mock_iqtree, threads=1, table_format="csv",
        )
