"""Tests for dating_diagnostics helpers."""
from __future__ import annotations
from pathlib import Path
from phyloai.posttree.dating_diagnostics import (
    parse_mcmctree_out,
    extract_node_tree,
    build_time_table,
)


SAMPLE_OUT = """\
Posterior means and 95% Equal-tail CIs

t_n7          0.4213 (0.3521, 0.4891)
t_n8          0.3102 (0.2641, 0.3589)
t_n9          0.5521 (0.4822, 0.6198)

Species tree for FigTree.  Branch lengths = posterior mean times; 95% CIs = labels
((sp1:0.32,sp2:0.31) 7 :0.12,sp3:0.44) 8 ;

((sp1,sp2) 7 ,sp3) 8 ;

(sp1,sp2,sp3);
"""


def test_parse_mcmctree_out_extracts_times():
    rows = parse_mcmctree_out(SAMPLE_OUT)
    assert len(rows) == 3
    assert rows[0]["node"] == "t_n7"
    assert abs(rows[0]["mean"] - 0.4213) < 1e-4
    assert abs(rows[0]["lower"] - 0.3521) < 1e-4
    assert abs(rows[0]["upper"] - 0.4891) < 1e-4
    assert abs(rows[0]["ci_width"] - (0.4891 - 0.3521)) < 1e-4


def test_extract_node_tree_returns_first_tree_with_bare_integers():
    tree = extract_node_tree(SAMPLE_OUT)
    assert tree is not None
    assert " 7 " in tree or ")7" in tree.replace(" ", "")


def test_build_time_table_two_runs():
    run1 = [
        {"node": "t_n7", "mean": 0.42, "lower": 0.35, "upper": 0.49, "ci_width": 0.14},
        {"node": "t_n8", "mean": 0.31, "lower": 0.26, "upper": 0.36, "ci_width": 0.10},
    ]
    run2 = [
        {"node": "t_n7", "mean": 0.43, "lower": 0.36, "upper": 0.50, "ci_width": 0.14},
        {"node": "t_n8", "mean": 0.30, "lower": 0.25, "upper": 0.35, "ci_width": 0.10},
    ]
    table = build_time_table(run1, run2)
    assert len(table) == 2
    assert "mean_run1" in table[0]
    assert "mean_run2" in table[0]
    assert "ci_width_run1" in table[0]
