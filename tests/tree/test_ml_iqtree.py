from __future__ import annotations

import shutil
from pathlib import Path

import pytest


# --- Constants from spec ---

IQTREE_EXTENSIONS = frozenset({
    ".fa", ".fas", ".fasta", ".faa", ".fna",
    ".phy", ".phylip", ".nex", ".nxs", ".nexus", ".aln",
})

AA_STANDARD_MODELS = frozenset({
    "LG", "Poisson", "cpREV", "mtREV", "Dayhoff", "mtMAM", "JTT", "WAG",
    "mtART", "mtZOA", "VT", "rtREV", "DCMut", "PMB", "HIVb", "HIVw",
    "JTTDCMut", "FLU", "Blosum62", "GTR20", "mtMet", "mtVer", "mtInv",
    "FLAVI", "Q.LG", "Q.pfam", "Q.pfam_gb", "Q.bird", "Q.mammal",
    "Q.insect", "Q.plant", "Q.yeast",
})

NT_STANDARD_MODELS = frozenset({
    "GTR", "HKY", "JC", "F81", "K2P", "K3P", "K81uf", "TN", "TNef",
    "TIM", "TIMef", "TVM", "TVMef", "SYM",
})

AA_MIXTURE_MODELS = frozenset(
    {f"C{i}" for i in range(10, 61, 10)}
    | {"EX2", "EX3", "EHO", "UL2", "UL3", "EX_EHO", "LG4M", "LG4X"}
)


# ===================================================================
# _scan_input_iqtree
# ===================================================================

def test_iqtree_extensions_coverage() -> None:
    assert ".aln" in IQTREE_EXTENSIONS
    assert ".nex" in IQTREE_EXTENSIONS
    assert ".nxs" in IQTREE_EXTENSIONS
    assert ".nexus" in IQTREE_EXTENSIONS
    assert ".fa" in IQTREE_EXTENSIONS
    assert ".phy" in IQTREE_EXTENSIONS


def test_scan_input_iqtree_finds_all_supported(tmp_path: Path) -> None:
    from phyloai.tree.ml_iqtree import _scan_input_iqtree

    (tmp_path / "gene1.fa").write_text(">a\nACGT\n")
    (tmp_path / "gene2.faa").write_text(">b\nMKT\n")
    (tmp_path / "gene3.fna").write_text(">c\nACGT\n")
    (tmp_path / "gene4.phy").write_text("2 10\na  ACGT\nb  ACGT\n")
    (tmp_path / "gene5.phylip").write_text("2 10\na  ACGT\nb  ACGT\n")
    (tmp_path / "gene6.nex").write_text(
        "#NEXUS\nbegin data; dimensions ntax=2 nchar=5; matrix\na ACGT-\nb ACGT-\n;\nend;"
    )
    (tmp_path / "gene7.nxs").write_text(
        "#NEXUS\nbegin data; dimensions ntax=2 nchar=5; matrix\na ACGT-\nb ACGT-\n;\nend;"
    )
    (tmp_path / "gene8.nexus").write_text(
        "#NEXUS\nbegin data; dimensions ntax=2 nchar=5; matrix\na ACGT-\nb ACGT-\n;\nend;"
    )
    (tmp_path / "gene9.aln").write_text("CLUSTAL\n\na ACGT\nb ACGT\n")
    (tmp_path / "notes.txt").write_text("skip")
    (tmp_path / "empty.fa").write_text("")
    (tmp_path / "subdir").mkdir()

    found, skipped = _scan_input_iqtree(tmp_path)

    assert len(found) == 9
    assert len(skipped) == 3
    skip_reasons = {s["reason"] for s in skipped}
    assert "empty file" in skip_reasons
    assert "directory" in skip_reasons
    assert "unrecognized extension: .txt" in skip_reasons


# ===================================================================
# _check_managed_flag_conflict + _is_flag_overridden
# ===================================================================

def test_check_managed_flag_conflict_blocks_s_flag() -> None:
    from phyloai.tree.ml_iqtree import _check_managed_flag_conflict

    with pytest.raises(ValueError, match="Blocked managed flag"):
        _check_managed_flag_conflict("-s some.fa")


def test_check_managed_flag_conflict_blocks_io_redirect_gt() -> None:
    from phyloai.tree.ml_iqtree import _check_managed_flag_conflict

    with pytest.raises(ValueError, match="Blocked I/O override"):
        _check_managed_flag_conflict("> out.txt")


def test_check_managed_flag_conflict_blocks_io_redirect_lt() -> None:
    from phyloai.tree.ml_iqtree import _check_managed_flag_conflict

    with pytest.raises(ValueError, match="Blocked I/O override"):
        _check_managed_flag_conflict("< input.txt")


def test_check_managed_flag_conflict_blocks_pipe() -> None:
    from phyloai.tree.ml_iqtree import _check_managed_flag_conflict

    with pytest.raises(ValueError, match="Blocked I/O override"):
        _check_managed_flag_conflict("| tee log.txt")


def test_check_managed_flag_conflict_blocks_prefix_in_batch() -> None:
    from phyloai.tree.ml_iqtree import _check_managed_flag_conflict

    with pytest.raises(ValueError, match="Blocked managed flag.*--prefix"):
        _check_managed_flag_conflict("--prefix shared", batch_mode=True)


def test_check_managed_flag_conflict_allows_prefix_in_single() -> None:
    from phyloai.tree.ml_iqtree import _check_managed_flag_conflict

    _check_managed_flag_conflict("--prefix myname", batch_mode=False)


def test_check_managed_flag_conflict_allows_strategy_args() -> None:
    from phyloai.tree.ml_iqtree import _check_managed_flag_conflict

    _check_managed_flag_conflict("-pers 0.5 -nstop 500 -nm 2000")


def test_is_flag_overridden_detects_short_flag() -> None:
    from phyloai.tree.ml_iqtree import _is_flag_overridden

    tokens = {"-m", "MFP", "-T", "4"}
    assert _is_flag_overridden("-m", tokens) is True
    assert _is_flag_overridden("-T", tokens) is True
    assert _is_flag_overridden("-B", tokens) is False


def test_is_flag_overridden_detects_long_flag() -> None:
    from phyloai.tree.ml_iqtree import _is_flag_overridden

    tokens = {"--merge", "--rclusterf", "10"}
    assert _is_flag_overridden("--merge", tokens) is True
    assert _is_flag_overridden("--fast", tokens) is False


def test_is_flag_overridden_detects_ufboot_aliases() -> None:
    from phyloai.tree.ml_iqtree import _is_flag_overridden

    tokens = {"-B", "1000"}
    assert _is_flag_overridden("-B", tokens) is True
    assert _is_flag_overridden("--ufboot", tokens) is True

    tokens2 = {"--ufboot", "500"}
    assert _is_flag_overridden("-B", tokens2) is True
    assert _is_flag_overridden("--ufboot", tokens2) is True


# ===================================================================
# _classify_workflow
# ===================================================================

def test_classify_workflow_homogeneous_no_partition_none() -> None:
    from phyloai.tree.ml_iqtree import _classify_workflow

    result = _classify_workflow(
        modelfinder="none", model="LG", seq_type="AA",
        partitions=None,
    )
    assert result == "homogeneous-no-partition-none"


def test_classify_workflow_homogeneous_no_partition_mf() -> None:
    from phyloai.tree.ml_iqtree import _classify_workflow

    result = _classify_workflow(
        modelfinder="MF", model="LG", seq_type="AA",
        partitions=None,
    )
    assert result == "homogeneous-no-partition-MF"


def test_classify_workflow_homogeneous_no_partition_mfp() -> None:
    from phyloai.tree.ml_iqtree import _classify_workflow

    result = _classify_workflow(
        modelfinder="MFP", model="LG", seq_type="AA",
        partitions=None,
    )
    assert result == "homogeneous-no-partition-MFP"


def test_classify_workflow_homogeneous_partition_none() -> None:
    from phyloai.tree.ml_iqtree import _classify_workflow

    result = _classify_workflow(
        modelfinder="none", model="LG", seq_type="AA",
        partitions="/some/p.txt",
    )
    assert result == "homogeneous-partition-none"


def test_classify_workflow_homogeneous_partition_mf_merge() -> None:
    from phyloai.tree.ml_iqtree import _classify_workflow

    result = _classify_workflow(
        modelfinder="MF", model="LG", seq_type="AA",
        partitions="/some/p.txt", rclusterf=10,
    )
    assert result == "homogeneous-partition-MF-merge"


def test_classify_workflow_homogeneous_partition_mfp_merge() -> None:
    from phyloai.tree.ml_iqtree import _classify_workflow

    result = _classify_workflow(
        modelfinder="MFP", model="LG", seq_type="AA",
        partitions="/some/p.txt", rcluster_max=50,
    )
    assert result == "homogeneous-partition-MFP-merge"


def test_classify_workflow_aa_heterogeneous_direct() -> None:
    from phyloai.tree.ml_iqtree import _classify_workflow

    result = _classify_workflow(
        modelfinder="none", model="C20", seq_type="AA",
        partitions=None,
    )
    assert result == "AA-heterogeneous-direct"


def test_classify_workflow_aa_heterogeneous_pmsf() -> None:
    from phyloai.tree.ml_iqtree import _classify_workflow

    result = _classify_workflow(
        modelfinder="none", model="C20", seq_type="AA",
        partitions=None, guide_tree="/some/guide.nwk",
    )
    assert result == "AA-heterogeneous-PMSF"


def test_classify_workflow_nt_heterogeneous() -> None:
    from phyloai.tree.ml_iqtree import _classify_workflow

    result = _classify_workflow(
        modelfinder="none", model="MIX+MF", seq_type="NT",
        partitions=None,
    )
    assert result == "NT-heterogeneous"


# ===================================================================
# _build_model_string
# ===================================================================

def test_build_model_string_homogeneous() -> None:
    from phyloai.tree.ml_iqtree import _build_model_string

    result = _build_model_string(
        model="LG", state_freq="+F", rate_heterogeneity="+R4",
        modelfinder="none",
    )
    assert result == "LG+F+R4"


def test_build_model_string_no_freq() -> None:
    from phyloai.tree.ml_iqtree import _build_model_string

    result = _build_model_string(
        model="GTR", state_freq="none", rate_heterogeneity="+R4",
        modelfinder="none",
    )
    assert result == "GTR+R4"


def test_build_model_string_no_rate() -> None:
    from phyloai.tree.ml_iqtree import _build_model_string

    result = _build_model_string(
        model="GTR", state_freq="+F", rate_heterogeneity="none",
        modelfinder="none",
    )
    assert result == "GTR+F"


def test_build_model_string_no_freq_no_rate() -> None:
    from phyloai.tree.ml_iqtree import _build_model_string

    result = _build_model_string(
        model="HKY", state_freq="none", rate_heterogeneity="none",
        modelfinder="none",
    )
    assert result == "HKY"


def test_build_model_string_mf() -> None:
    from phyloai.tree.ml_iqtree import _build_model_string

    result = _build_model_string(
        model="LG", state_freq="+F", rate_heterogeneity="+R4",
        modelfinder="MF",
    )
    assert result == "MF"


def test_build_model_string_mfp() -> None:
    from phyloai.tree.ml_iqtree import _build_model_string

    result = _build_model_string(
        model="LG", state_freq="+F", rate_heterogeneity="+R4",
        modelfinder="MFP",
    )
    assert result == "MFP"


def test_build_model_string_pmsf() -> None:
    from phyloai.tree.ml_iqtree import _build_model_string

    result = _build_model_string(
        model="C20", state_freq="+F", rate_heterogeneity="+R4",
        modelfinder="none", pmsf_base_model="LG",
    )
    assert result == "LG+C20+F+R4"


def test_build_model_string_direct_mixture() -> None:
    from phyloai.tree.ml_iqtree import _build_model_string

    result = _build_model_string(
        model="C20", state_freq="+F", rate_heterogeneity="+R4",
        modelfinder="none",
    )
    assert result == "C20+F+R4"


def test_build_model_string_nt_heterogeneous() -> None:
    from phyloai.tree.ml_iqtree import _build_model_string

    result = _build_model_string(
        model="MIX+MF", state_freq="+F", rate_heterogeneity="+R4",
        modelfinder="none",
    )
    assert result == "MIX+MF"


def test_build_model_string_custom_exchangeabilities(tmp_path: Path) -> None:
    from phyloai.tree.ml_iqtree import _build_model_string

    model = tmp_path / "chain1.exchangeabilities"
    model.write_text("0.5\n")

    assert _build_model_string(
        model=str(model.resolve()), state_freq="none",
        rate_heterogeneity="+R4", modelfinder="none",
    ) == f"{model.resolve()}+R4"


# ===================================================================
# _validate_model + _validate_pmsf_base_model
# ===================================================================

def test_validate_model_aa_standard() -> None:
    from phyloai.tree.ml_iqtree import _validate_model
    _validate_model(model="LG", seq_type="AA", modelfinder="none")


def test_validate_model_nt_standard() -> None:
    from phyloai.tree.ml_iqtree import _validate_model
    _validate_model(model="GTR", seq_type="NT", modelfinder="none")


def test_validate_model_rejects_nt_model_for_aa() -> None:
    from phyloai.tree.ml_iqtree import _validate_model
    with pytest.raises(ValueError, match="Invalid model"):
        _validate_model(model="GTR", seq_type="AA", modelfinder="none")


def test_validate_model_rejects_aa_model_for_nt() -> None:
    from phyloai.tree.ml_iqtree import _validate_model
    with pytest.raises(ValueError, match="Invalid model"):
        _validate_model(model="LG", seq_type="NT", modelfinder="none")


def test_validate_model_allows_aa_mixture() -> None:
    from phyloai.tree.ml_iqtree import _validate_model
    _validate_model(model="C20", seq_type="AA", modelfinder="none")


def test_validate_model_allows_nt_heterogeneous() -> None:
    from phyloai.tree.ml_iqtree import _validate_model
    _validate_model(model="MIX+MF", seq_type="NT", modelfinder="none")


def test_validate_model_skips_when_mf_active() -> None:
    from phyloai.tree.ml_iqtree import _validate_model
    _validate_model(model="INVALID", seq_type="AA", modelfinder="MFP")


def test_validate_pmsf_base_model_allows_standard() -> None:
    from phyloai.tree.ml_iqtree import _validate_pmsf_base_model
    _validate_pmsf_base_model("LG")


def test_validate_pmsf_base_model_rejects_mixture() -> None:
    from phyloai.tree.ml_iqtree import _validate_pmsf_base_model
    with pytest.raises(ValueError, match="Invalid PMSF base model"):
        _validate_pmsf_base_model("C20")


# ===================================================================
# _run_validations
# ===================================================================

def test_validate_heterogeneous_requires_matrix() -> None:
    from phyloai.tree.ml_iqtree import _run_validations
    with pytest.raises(ValueError, match="only supported in --matrix"):
        _run_validations(
            batch_mode=True, seq_type="AA", modelfinder="none",
            model="C20", partitions=None, guide_tree=None,
        )


def test_validate_nt_heterogeneous_requires_matrix() -> None:
    from phyloai.tree.ml_iqtree import _run_validations
    with pytest.raises(ValueError, match="only supported in --matrix"):
        _run_validations(
            batch_mode=True, seq_type="NT", modelfinder="none",
            model="MIX+MF", partitions=None, guide_tree=None,
        )


def test_validate_partitions_requires_matrix() -> None:
    from phyloai.tree.ml_iqtree import _run_validations
    with pytest.raises(ValueError, match="only valid with --matrix"):
        _run_validations(
            batch_mode=True, seq_type="AA", modelfinder="MFP",
            model="LG", partitions="/some/p.txt", guide_tree=None,
        )


def test_validate_boot_bnni_without_boot_warns() -> None:
    from phyloai.tree.ml_iqtree import _run_validations
    with pytest.warns(UserWarning, match="--bnni has no effect"):
        _run_validations(
            batch_mode=False, seq_type="AA", modelfinder="none",
            model="LG", partitions=None, guide_tree=None,
            bnni=True, boot=None,
        )


def test_validate_mf_ignores_branch_support_warns() -> None:
    from phyloai.tree.ml_iqtree import _run_validations
    with pytest.warns(UserWarning, match="Branch support flags"):
        _run_validations(
            batch_mode=False, seq_type="AA", modelfinder="MF",
            model="LG", partitions=None, guide_tree=None,
boot=1000, alrt=None, bnni=False        )


def test_validate_prefix_in_batch_warns() -> None:
    from phyloai.tree.ml_iqtree import _run_validations
    with pytest.warns(UserWarning, match="--prefix ignored in batch"):
        _run_validations(
            batch_mode=True, seq_type="AA", modelfinder="none",
            model="LG", partitions=None, guide_tree=None,
            prefix="myprefix",
        )


def test_validate_pmsf_requires_guide_tree() -> None:
    from phyloai.tree.ml_iqtree import _run_validations
    with pytest.raises(ValueError, match="PMSF mode requires --guide-tree"):
        _run_validations(
            batch_mode=False, seq_type="AA", modelfinder="none",
            model="C20", partitions=None, guide_tree=None,
            pmsf_base_model="LG",
        )


def test_validate_pmsf_base_model_requires_mixture_model() -> None:
    from phyloai.tree.ml_iqtree import _run_validations
    with pytest.raises(ValueError, match="only valid with AA mixture models"):
        _run_validations(
            batch_mode=False, seq_type="AA", modelfinder="none",
            model="LG", partitions=None, guide_tree=None,
            pmsf_base_model="LG",
        )


@pytest.mark.parametrize(
    ("batch_mode", "seq_type", "modelfinder", "state_freq", "custom_model", "match"),
    [
        (True, "AA", "none", "none", True, "--matrix"),
        (False, "NT", "none", "none", True, "AA"),
        (False, "AA", "MF", "none", True, "ModelFinder"),
        (False, "AA", "none", "+F", True, "--state-freq none"),
        (False, "AA", "none", "none", False, "custom model"),
    ],
)
def test_validate_site_freq_file_rejects_invalid_context(
    batch_mode: bool,
    seq_type: str,
    modelfinder: str,
    state_freq: str,
    custom_model: bool,
    match: str,
) -> None:
    from phyloai.tree.ml_iqtree import _run_validations

    with pytest.raises(ValueError, match=match):
        _run_validations(
            batch_mode=batch_mode, seq_type=seq_type, modelfinder=modelfinder,
            model="/tmp/chain1.exchangeabilities" if custom_model else "LG",
            partitions=None, guide_tree=None, state_freq=state_freq,
            custom_model=custom_model, site_freq_file="/tmp/chain1.sitefreq",
        )


def test_validate_tool_args_fs_requires_no_state_frequency() -> None:
    from phyloai.tree.ml_iqtree import _run_validations

    with pytest.raises(ValueError, match="--state-freq none"):
        _run_validations(
            batch_mode=False, seq_type="AA", modelfinder="none",
            model="/tmp/chain1.exchangeabilities", partitions=None, guide_tree=None,
            state_freq="+F", custom_model=True, tool_args="-fs /tmp/override.sitefreq",
        )


def test_validate_rcluster_without_partitions_warns() -> None:
    from phyloai.tree.ml_iqtree import _run_validations
    with pytest.warns(UserWarning, match="no effect without --partitions"):
        _run_validations(
            batch_mode=False, seq_type="AA", modelfinder="none",
            model="LG", partitions=None, guide_tree=None,
            rclusterf=10,
        )


def test_validate_rcluster_mutually_exclusive() -> None:
    from phyloai.tree.ml_iqtree import _run_validations
    with pytest.raises(ValueError, match="mutually exclusive"):
        _run_validations(
            batch_mode=False, seq_type="AA", modelfinder="MFP",
            model="LG", partitions="/some/p.txt", guide_tree=None,
            rclusterf=10, rcluster_max=50,
        )


def test_validate_qmax_without_mix_mf_warns() -> None:
    from phyloai.tree.ml_iqtree import _run_validations
    with pytest.warns(UserWarning, match="--qmax only takes effect"):
        _run_validations(
            batch_mode=False, seq_type="NT", modelfinder="none",
            model="GTR", partitions=None, guide_tree=None,
            qmax=10,
        )


# ===================================================================
# _validate_seq_types_iqtree
# ===================================================================

def test_validate_seq_types_iqtree_fasta(tmp_path: Path) -> None:
    from phyloai.tree.ml_iqtree import _validate_seq_types_iqtree

    (tmp_path / "g1.fa").write_text(">a\nMKTLLL\n")
    (tmp_path / "g2.fa").write_text(">b\nMKTLLL\n")
    files = sorted(tmp_path.glob("*.fa"))

    resolved, offending = _validate_seq_types_iqtree(files, declared_type=None)
    assert resolved == "AA"
    assert len(offending) == 0


def test_validate_seq_types_iqtree_mixed_raises(tmp_path: Path) -> None:
    from phyloai.tree.ml_iqtree import _validate_seq_types_iqtree

    (tmp_path / "g1.fa").write_text(">a\nMKTLLL\n")
    (tmp_path / "g2.fa").write_text(">b\nACGTAC\n")
    files = sorted(tmp_path.glob("*.fa"))

    resolved, offending = _validate_seq_types_iqtree(files, declared_type=None)
    assert resolved is None
    assert len(offending) >= 1


def test_validate_seq_types_iqtree_explicit_mismatch(tmp_path: Path) -> None:
    from phyloai.tree.ml_iqtree import _validate_seq_types_iqtree

    (tmp_path / "g1.fa").write_text(">a\nMKTLLL\n")
    files = sorted(tmp_path.glob("*.fa"))

    resolved, offending = _validate_seq_types_iqtree(files, declared_type="NT")
    assert resolved == "NT"
    assert len(offending) == 1


def test_validate_seq_types_iqtree_no_files(tmp_path: Path) -> None:
    from phyloai.tree.ml_iqtree import _validate_seq_types_iqtree

    resolved, offending = _validate_seq_types_iqtree([], declared_type=None)
    assert resolved == "AA"
    assert len(offending) == 0


def test_validate_seq_types_iqtree_nexus_format(tmp_path: Path) -> None:
    from phyloai.tree.ml_iqtree import _validate_seq_types_iqtree

    (tmp_path / "g1.nex").write_text(
        "#NEXUS\nbegin data;\ndimensions ntax=2 nchar=10;\nformat datatype=protein;\nmatrix\na  MKTLLLMKTK\nb  MKTLLLMKTK\n;\nend;"
    )
    files = sorted(tmp_path.glob("*.nex"))

    resolved, offending = _validate_seq_types_iqtree(files, declared_type=None)
    assert resolved == "AA"
    assert len(offending) == 0


def test_validate_seq_types_iqtree_unparsable_skipped(tmp_path: Path) -> None:
    from phyloai.tree.ml_iqtree import _validate_seq_types_iqtree

    (tmp_path / "bad.nex").write_text("not a valid nexus file\n")
    files = sorted(tmp_path.glob("*.nex"))

    resolved, offending = _validate_seq_types_iqtree(files, declared_type=None)
    # unparsable nexus returns 0 records -> "no sequences found" or parse failure
    assert len(offending) == 1
    assert "parse" in offending[0]["reason"].lower() or "no sequences" in offending[0]["reason"].lower()


# ===================================================================
# _parse_threads
# ===================================================================

def test_parse_threads_none_batch() -> None:
    from phyloai.tree.ml_iqtree import _parse_threads
    val = _parse_threads(None, batch_mode=True)
    assert val == 4


def test_parse_threads_none_single() -> None:
    from phyloai.tree.ml_iqtree import _parse_threads
    val = _parse_threads(None, batch_mode=False)
    assert val == "auto"


def test_parse_threads_auto_single() -> None:
    from phyloai.tree.ml_iqtree import _parse_threads
    val = _parse_threads("auto", batch_mode=False)
    assert val == "auto"


def test_parse_threads_auto_batch_rejected() -> None:
    from phyloai.tree.ml_iqtree import _parse_threads
    with pytest.raises(ValueError, match="--threads"):
        _parse_threads("auto", batch_mode=True)


def test_parse_threads_int_single() -> None:
    from phyloai.tree.ml_iqtree import _parse_threads
    val = _parse_threads("4", batch_mode=False)
    assert val == 4


def test_parse_threads_int_batch() -> None:
    from phyloai.tree.ml_iqtree import _parse_threads
    val = _parse_threads("8", batch_mode=True)
    assert val == 8


def test_parse_threads_zero_rejected() -> None:
    from phyloai.tree.ml_iqtree import _parse_threads
    with pytest.raises(ValueError, match="--threads"):
        _parse_threads("0", batch_mode=False)


def test_parse_threads_negative_rejected() -> None:
    from phyloai.tree.ml_iqtree import _parse_threads
    with pytest.raises(ValueError, match="--threads"):
        _parse_threads("-1", batch_mode=False)


# ===================================================================
# _build_iqtree_cmd
# ===================================================================

def test_build_iqtree_cmd_basic_homogeneous_aa(tmp_path: Path) -> None:
    from phyloai.tree.ml_iqtree import _build_iqtree_cmd

    inp = tmp_path / "gene.fa"
    pref = tmp_path / "gene"

    cmd = _build_iqtree_cmd(
        input_path=inp, prefix=pref,
        model_string="LG+F+R4", seq_type="AA",
        modelfinder="none", boot=None, alrt=None, bnni=False, mode="normal", threads_arg="-T AUTO",
    )

    assert cmd[0] == "iqtree3"
    flat = " ".join(cmd)
    assert "-s" in cmd
    assert "--prefix" in cmd
    assert "-m LG+F+R4" in flat
    assert "--seqtype AA" in flat
    assert "-T AUTO" in flat
    assert "--redo" not in cmd  # no --redo by default


def test_build_iqtree_cmd_basic_homogeneous_nt(tmp_path: Path) -> None:
    from phyloai.tree.ml_iqtree import _build_iqtree_cmd

    inp = tmp_path / "gene.fa"
    pref = tmp_path / "gene"

    cmd = _build_iqtree_cmd(
        input_path=inp, prefix=pref,
        model_string="GTR+F+R4", seq_type="NT",
        modelfinder="none", boot=None, alrt=None, bnni=False, mode="normal", threads_arg="-T AUTO",
    )

    flat = " ".join(cmd)
    assert "--seqtype DNA" in flat


def test_build_iqtree_cmd_seqtype_auto(tmp_path: Path) -> None:
    from phyloai.tree.ml_iqtree import _build_iqtree_cmd

    inp = tmp_path / "gene.fa"
    pref = tmp_path / "gene"

    cmd = _build_iqtree_cmd(
        input_path=inp, prefix=pref,
        model_string="LG+F+R4", seq_type="auto",
        boot=None, alrt=None, bnni=False,
        modelfinder="none", mode="normal", threads_arg="-T AUTO",
    )

    flat = " ".join(cmd)
    assert "--seqtype" not in cmd


def test_build_iqtree_cmd_with_boot_and_alrt(tmp_path: Path) -> None:
    from phyloai.tree.ml_iqtree import _build_iqtree_cmd

    inp = tmp_path / "gene.fa"
    pref = tmp_path / "gene"

    cmd = _build_iqtree_cmd(
        input_path=inp, prefix=pref,
        model_string="GTR+F+R4", seq_type="NT",
        boot=1000, alrt=1000, bnni=True,
        modelfinder="none", mode="normal", threads_arg="-T AUTO",
    )

    flat = " ".join(cmd)
    assert "-B 1000" in flat
    assert "--alrt 1000" in flat
    assert "--bnni" in cmd


def test_build_iqtree_cmd_alrt_zero_parametric(tmp_path: Path) -> None:
    from phyloai.tree.ml_iqtree import _build_iqtree_cmd

    inp = tmp_path / "gene.fa"
    pref = tmp_path / "gene"

    cmd = _build_iqtree_cmd(
        input_path=inp, prefix=pref,
        model_string="LG+F+R4", seq_type="AA",
        boot=None, alrt=0, bnni=False,
        modelfinder="none", mode="normal", threads_arg="-T AUTO",
    )

    flat = " ".join(cmd)
    assert "--alrt 0" in flat
    assert "-B" not in cmd


def test_build_iqtree_cmd_modelfinder_mf(tmp_path: Path) -> None:
    from phyloai.tree.ml_iqtree import _build_iqtree_cmd

    inp = tmp_path / "gene.fa"
    pref = tmp_path / "gene"

    cmd = _build_iqtree_cmd(
        input_path=inp, prefix=pref,
        model_string="MF", seq_type="AA",
        modelfinder="MF", boot=None, alrt=None, bnni=False, mode="normal", threads_arg="-T AUTO",
        mset="LG,WAG", msub="nuclear",
    )

    flat = " ".join(cmd)
    assert "-m MF" in flat
    assert "--mset LG,WAG" in flat
    assert "--msub nuclear" in flat


def test_build_iqtree_cmd_mset_all_omitted(tmp_path: Path) -> None:
    from phyloai.tree.ml_iqtree import _build_iqtree_cmd

    inp = tmp_path / "gene.fa"
    pref = tmp_path / "gene"

    cmd = _build_iqtree_cmd(
        input_path=inp, prefix=pref,
        model_string="MF", seq_type="AA",
        modelfinder="none", boot=None, alrt=None, bnni=False, mode="normal", threads_arg="-T AUTO",
        mset="all",
    )

    assert "--mset" not in cmd


def test_build_iqtree_cmd_fast_mode(tmp_path: Path) -> None:
    from phyloai.tree.ml_iqtree import _build_iqtree_cmd

    inp = tmp_path / "gene.fa"
    pref = tmp_path / "gene"

    cmd = _build_iqtree_cmd(
        input_path=inp, prefix=pref,
        model_string="LG+F+R4", seq_type="AA",
        boot=None, alrt=None, bnni=False,
        modelfinder="none", mode="fast", threads_arg="-T AUTO",
    )
    assert "--fast" in cmd


def test_build_iqtree_cmd_partitions_and_merge(tmp_path: Path) -> None:
    from phyloai.tree.ml_iqtree import _build_iqtree_cmd

    inp = tmp_path / "gene.fa"
    pref = tmp_path / "gene"

    cmd = _build_iqtree_cmd(
        input_path=inp, prefix=pref,
        model_string="MFP", seq_type="AA",
        boot=None, alrt=None, bnni=False,
        modelfinder="MFP", mode="normal", threads_arg="-T AUTO",
        partitions="/some/p.txt",
        rclusterf=10,
    )

    flat = " ".join(cmd)
    assert "-p /some/p.txt" in flat
    assert "--merge" in cmd
    assert "--rclusterf 10" in flat


def test_build_iqtree_cmd_pmsf_guide_tree(tmp_path: Path) -> None:
    from phyloai.tree.ml_iqtree import _build_iqtree_cmd

    inp = tmp_path / "gene.fa"
    pref = tmp_path / "gene"

    cmd = _build_iqtree_cmd(
        input_path=inp, prefix=pref,
        model_string="LG+C20+F+R4", seq_type="AA",
        modelfinder="none", boot=None, alrt=None, bnni=False, mode="normal", threads_arg="-T AUTO",
        guide_tree="/some/guide.nwk",
    )

    flat = " ".join(cmd)
    assert "-ft /some/guide.nwk" in flat


def test_build_iqtree_cmd_mix_mf_with_qmax(tmp_path: Path) -> None:
    from phyloai.tree.ml_iqtree import _build_iqtree_cmd

    inp = tmp_path / "gene.fa"
    pref = tmp_path / "gene"

    cmd = _build_iqtree_cmd(
        input_path=inp, prefix=pref,
        model_string="MIX+MF", seq_type="NT",
        modelfinder="none", boot=None, alrt=None, bnni=False, mode="normal", threads_arg="-T AUTO",
        mset="GTR,HKY", qmax=10,
    )

    flat = " ".join(cmd)
    assert "-m MIX+MF" in flat
    assert "-qmax 10" in flat
    assert "--mset GTR,HKY" in flat


def test_build_iqtree_cmd_output_flags(tmp_path: Path) -> None:
    from phyloai.tree.ml_iqtree import _build_iqtree_cmd

    inp = tmp_path / "gene.fa"
    pref = tmp_path / "gene"

    cmd = _build_iqtree_cmd(
        input_path=inp, prefix=pref,
        model_string="LG+F+R4", seq_type="AA",
        modelfinder="none", boot=None, alrt=None, bnni=False, mode="normal", threads_arg="-T AUTO",
        rate=True, wslr=True,
        constraint="/some/constraint.nwk",
        outgroup="taxon1,taxon2",
    )

    flat = " ".join(cmd)
    assert "--rate" in cmd
    assert "-wslr" in cmd
    assert "-g /some/constraint.nwk" in flat
    assert "-o taxon1,taxon2" in flat


def test_build_iqtree_cmd_tool_args_appended(tmp_path: Path) -> None:
    from phyloai.tree.ml_iqtree import _build_iqtree_cmd

    inp = tmp_path / "gene.fa"
    pref = tmp_path / "gene"

    cmd = _build_iqtree_cmd(
        input_path=inp, prefix=pref,
        model_string="LG+F+R4", seq_type="AA",
        modelfinder="none", boot=None, alrt=None, bnni=False, mode="normal", threads_arg="-T AUTO",
        tool_args="-pers 0.5 -nstop 500",
    )

    assert "-pers" in cmd
    assert "0.5" in cmd
    assert "-nstop" in cmd
    assert "500" in cmd


def test_build_iqtree_cmd_tool_args_overrides_model(tmp_path: Path) -> None:
    from phyloai.tree.ml_iqtree import _build_iqtree_cmd

    inp = tmp_path / "gene.fa"
    pref = tmp_path / "gene"

    cmd = _build_iqtree_cmd(
        input_path=inp, prefix=pref,
        model_string="LG+F+R4", seq_type="AA",
        modelfinder="none", boot=None, alrt=None, bnni=False, mode="normal", threads_arg="-T AUTO",
        tool_args="-m WAG+F+R4",
    )

    assert cmd.count("-m") == 1
    flat = " ".join(cmd)
    assert "WAG+F+R4" in flat


def test_build_iqtree_cmd_tool_args_overrides_boot(tmp_path: Path) -> None:
    from phyloai.tree.ml_iqtree import _build_iqtree_cmd

    inp = tmp_path / "gene.fa"
    pref = tmp_path / "gene"

    cmd = _build_iqtree_cmd(
        input_path=inp, prefix=pref,
        model_string="LG+F+R4", seq_type="AA",
        modelfinder="none", boot=None, alrt=None, bnni=False, mode="normal", threads_arg="-T AUTO",
        tool_args="-B 500",
    )

    assert cmd.count("-B") == 1
    flat = " ".join(cmd)
    assert "500" in flat
    assert "1000" not in flat


def test_build_iqtree_cmd_tool_args_blocked_s(tmp_path: Path) -> None:
    from phyloai.tree.ml_iqtree import _build_iqtree_cmd

    inp = tmp_path / "gene.fa"
    pref = tmp_path / "gene"

    with pytest.raises(ValueError, match="Blocked managed flag"):
        _build_iqtree_cmd(
            input_path=inp, prefix=pref,
            model_string="LG+F+R4", seq_type="AA",
            boot=None, alrt=None, bnni=False,
            mode="normal", threads_arg="-T AUTO",
            tool_args="-s hack.fa",
        )


def test_build_iqtree_cmd_adds_site_freq_file(tmp_path: Path) -> None:
    from phyloai.tree.ml_iqtree import _build_iqtree_cmd

    profile = tmp_path / "chain1.sitefreq"
    cmd = _build_iqtree_cmd(
        input_path=tmp_path / "matrix.fa", prefix=tmp_path / "matrix",
        model_string="/models/chain1.exchangeabilities+R4", seq_type="AA",
        modelfinder="none", boot=0, alrt=None, bnni=False,
        mode="normal", threads_arg="-T 1", site_freq_file=str(profile),
    )

    assert cmd.count("-fs") == 1
    assert cmd[cmd.index("-fs") + 1] == str(profile)


def test_build_iqtree_cmd_tool_args_fs_overrides_structured_profile(tmp_path: Path) -> None:
    from phyloai.tree.ml_iqtree import _build_iqtree_cmd

    cmd = _build_iqtree_cmd(
        input_path=tmp_path / "matrix.fa", prefix=tmp_path / "matrix",
        model_string="/models/chain1.exchangeabilities+R4", seq_type="AA",
        modelfinder="none", boot=0, alrt=None, bnni=False,
        mode="normal", threads_arg="-T 1", site_freq_file="/managed.sitefreq",
        tool_args="-fs /override.sitefreq",
    )

    assert cmd.count("-fs") == 1
    assert cmd[cmd.index("-fs") + 1] == "/override.sitefreq"


def test_build_iqtree_cmd_tool_args_blocked_pipe(tmp_path: Path) -> None:
    from phyloai.tree.ml_iqtree import _build_iqtree_cmd

    inp = tmp_path / "gene.fa"
    pref = tmp_path / "gene"

    with pytest.raises(ValueError, match="Blocked I/O override"):
        _build_iqtree_cmd(
            input_path=inp, prefix=pref,
            model_string="LG+F+R4", seq_type="AA",
            boot=None, alrt=None, bnni=False,
            mode="normal", threads_arg="-T AUTO",
            tool_args="| tee log.txt",
        )


# ===================================================================
# _resolve_iqtree_path + _detect_iqtree_version
# ===================================================================

def test_resolve_iqtree_path_custom(tmp_path: Path) -> None:
    from phyloai.tree.ml_iqtree import _resolve_iqtree_path

    fake = tmp_path / "iqtree3"
    fake.write_text("#!/bin/sh\necho IQ-TREE 3.1.2\n")
    fake.chmod(0o755)

    result = _resolve_iqtree_path(str(fake), dry_run=False)
    assert result == str(fake)


def test_resolve_iqtree_path_missing_raises() -> None:
    from phyloai.tree.ml_iqtree import _resolve_iqtree_path

    with pytest.raises(ValueError, match="does not exist"):
        _resolve_iqtree_path("/nonexistent/iqtree3", dry_run=False)


def test_resolve_iqtree_path_not_executable_raises(tmp_path: Path) -> None:
    from phyloai.tree.ml_iqtree import _resolve_iqtree_path

    fake = tmp_path / "iqtree3"
    fake.write_text("not executable")

    with pytest.raises(ValueError, match="not executable"):
        _resolve_iqtree_path(str(fake), dry_run=False)


def test_resolve_iqtree_path_dry_run() -> None:
    from phyloai.tree.ml_iqtree import _resolve_iqtree_path

    result = _resolve_iqtree_path(None, dry_run=True)
    assert result == "iqtree3"


def test_detect_iqtree_version_parses_output() -> None:
    from unittest.mock import patch
    from phyloai.tree.ml_iqtree import _detect_iqtree_version

    with patch("subprocess.run") as mock_run:
        mock_run.return_value.stdout = "IQ-TREE multicore version 3.1.2 for Linux 64-bit built Mar 15 2025"
        mock_run.return_value.stderr = ""
        result = _detect_iqtree_version("/usr/bin/iqtree3")
        assert result["iqtree3"] == "3.1.2"


def test_detect_iqtree_version_unknown_on_failure() -> None:
    from unittest.mock import patch
    from phyloai.tree.ml_iqtree import _detect_iqtree_version

    with patch("subprocess.run", side_effect=Exception("fail")):
        result = _detect_iqtree_version("/usr/bin/iqtree3")
        assert result["iqtree3"] == "unknown"


# ===================================================================
# _run_one_iqtree
# ===================================================================

def test_run_one_iqtree_dry_run(tmp_path: Path) -> None:
    from phyloai.tree.ml_iqtree import _run_one_iqtree

    inp = tmp_path / "gene.fa"
    inp.write_text(">a\nMKTLLL\n>b\nMKTLLL\n")
    log_dir = tmp_path / "logs"
    log_dir.mkdir()

    result = _run_one_iqtree(
        gene_path=inp, seq_type="AA",
        model_string="LG+F+R4", modelfinder="none",
        boot=None, alrt=None, bnni=False,
        mode="normal", threads_arg="-T AUTO",
        log_dir=log_dir, output_dir=tmp_path,
        dry_run=True,
    )

    assert result["status"] == "dry_run"
    assert "cmd" in result
    assert result["n_taxa"] == 2


def test_run_one_iqtree_single_prefix_tracks_all_outputs(tmp_path: Path) -> None:
    from phyloai.tree.ml_iqtree import _run_one_iqtree

    inp = tmp_path / "gene.fa"
    inp.write_text(">a\nMKTLLL\n>b\nMKTLLL\n")

    result = _run_one_iqtree(
        gene_path=inp, seq_type="AA",
        model_string="LG+F+R4", modelfinder="none",
        boot=None, alrt=None, bnni=False,
        mode="normal", threads_arg="-T AUTO",
        log_dir=tmp_path, output_dir=tmp_path,
        prefix="custom", dry_run=True,
    )

    assert result["output_tree"] == str(tmp_path / "custom.treefile")
    assert result["log_iqtree"] == str(tmp_path / "custom.iqtree")
    assert result["log_file"] == str(tmp_path / "custom.log")


def test_run_one_iqtree_single_tool_args_prefix_tracks_override(tmp_path: Path) -> None:
    from phyloai.tree.ml_iqtree import _run_one_iqtree

    inp = tmp_path / "gene.fa"
    inp.write_text(">a\nMKTLLL\n>b\nMKTLLL\n")

    result = _run_one_iqtree(
        gene_path=inp, seq_type="AA",
        model_string="LG+F+R4", modelfinder="none",
        boot=None, alrt=None, bnni=False,
        mode="normal", threads_arg="-T AUTO",
        log_dir=tmp_path, output_dir=tmp_path,
        tool_args="--prefix override", dry_run=True,
    )

    flat = " ".join(result["cmd"])
    assert "--prefix override" in flat
    assert f"--prefix {tmp_path / 'gene'}" not in flat
    assert result["output_tree"] == str(tmp_path / "override.treefile")
    assert result["log_iqtree"] == str(tmp_path / "override.iqtree")
    assert result["log_file"] == str(tmp_path / "override.log")


def test_run_one_iqtree_missing_input(tmp_path: Path) -> None:
    from phyloai.tree.ml_iqtree import _run_one_iqtree

    inp = tmp_path / "missing.fa"
    log_dir = tmp_path / "logs"
    log_dir.mkdir()

    result = _run_one_iqtree(
        gene_path=inp, seq_type="AA",
        model_string="LG+F+R4", modelfinder="none",
        boot=None, alrt=None, bnni=False,
        mode="normal", threads_arg="-T AUTO",
        log_dir=log_dir, output_dir=tmp_path,
    )

    assert result["status"] == "failed"
    assert "reason" in result


def test_run_one_iqtree_dry_run_batch(tmp_path: Path) -> None:
    from phyloai.tree.ml_iqtree import _run_one_iqtree

    inp = tmp_path / "gene.fa"
    inp.write_text(">a\nMKTLLL\n>b\nMKTLLL\n")
    log_dir = tmp_path / "logs"
    log_dir.mkdir()
    trees_dir = tmp_path / "trees"
    trees_dir.mkdir()

    result = _run_one_iqtree(
        gene_path=inp, seq_type="AA",
        model_string="LG+F+R4", modelfinder="none",
        boot=None, alrt=None, bnni=False,
        mode="normal", threads_arg="-T 1",
        log_dir=log_dir, output_dir=trees_dir,
        dry_run=True, batch_mode=True,
    )

    assert result["status"] == "dry_run"
    assert result["output_tree"] == str(trees_dir / "gene.treefile")
    assert result["log_iqtree"] == str(log_dir / "gene.iqtree")


def test_run_one_iqtree_dry_run_mf_mode(tmp_path: Path) -> None:
    from phyloai.tree.ml_iqtree import _run_one_iqtree

    inp = tmp_path / "gene.fa"
    inp.write_text(">a\nMKTLLL\n>b\nMKTLLL\n")
    log_dir = tmp_path / "logs"
    log_dir.mkdir()

    result = _run_one_iqtree(
        gene_path=inp, seq_type="AA",
        model_string="MF", modelfinder="MF",
        boot=None, alrt=None, bnni=False,
        mode="normal", threads_arg="-T AUTO",
        log_dir=log_dir, output_dir=tmp_path,
        dry_run=True,
    )

    assert result["status"] == "dry_run"
    assert result["output_tree"] is None  # MF mode: no tree


@pytest.mark.skipif(not shutil.which("iqtree3"), reason="iqtree3 not found in PATH")
def test_run_one_iqtree_success_single(tmp_path: Path) -> None:
    from phyloai.tree.ml_iqtree import _run_one_iqtree

    inp = tmp_path / "gene.fa"
    inp.write_text(">a\nMKTLLLMKT\n>b\nMKTLLLMKT\n>c\nMKTLLLMKT\n")

    out_dir = tmp_path / "out"
    out_dir.mkdir()

    result = _run_one_iqtree(
        gene_path=inp, seq_type="AA",
        model_string="LG+F+R4", modelfinder="none",
        boot=None, alrt=None, bnni=False,
        mode="normal", threads_arg="-T 2",
        log_dir=out_dir, output_dir=out_dir,
        executable="iqtree3", batch_mode=False,
    )

    assert result["status"] == "success"
    assert result["output_tree"] is not None
    assert Path(result["output_tree"]).exists()
    assert Path(result["log_iqtree"]).exists()


def test_run_one_iqtree_single_accepts_multiple_trees(tmp_path: Path) -> None:
    from unittest.mock import patch
    from phyloai.tree.ml_iqtree import _run_one_iqtree

    inp = tmp_path / "matrix.fa"
    inp.write_text(">a\nMKTLLL\n>b\nMKTLLL\n")
    out_dir = tmp_path / "out"
    out_dir.mkdir()

    class _FakePopen:
        returncode = 0

        def __init__(self, *args, **kwargs):
            (out_dir / "matrix.iqtree").write_text("Log-likelihood: -100.0\n")
            (out_dir / "matrix.treefile").write_text("(a:0.1,b:0.1);\n(a:0.2,b:0.2);\n")
            (out_dir / "matrix.log").write_text("done\n")

        def communicate(self):
            return "", ""

    with patch("subprocess.Popen", _FakePopen):
        result = _run_one_iqtree(
            gene_path=inp, seq_type="AA",
            model_string="LG+H4", modelfinder="none",
            boot=None, alrt=None, bnni=False,
            mode="normal", threads_arg="-T 1",
            log_dir=out_dir, output_dir=out_dir,
        )

    assert result["status"] == "success"


@pytest.mark.skipif(not shutil.which("iqtree3"), reason="iqtree3 not found in PATH")
def test_run_one_iqtree_success_batch(tmp_path: Path) -> None:
    from phyloai.tree.ml_iqtree import _run_one_iqtree

    inp = tmp_path / "gene.fa"
    inp.write_text(">a\nMKTLLLMKT\n>b\nMKTLLLMKT\n>c\nMKTLLLMKT\n")

    trees_dir = tmp_path / "trees"
    trees_dir.mkdir()
    log_dir = tmp_path / "logs"
    log_dir.mkdir()

    result = _run_one_iqtree(
        gene_path=inp, seq_type="AA",
        model_string="LG+F+R4", modelfinder="none",
        boot=None, alrt=None, bnni=False,
        mode="normal", threads_arg="-T 1",
        log_dir=log_dir, output_dir=trees_dir,
        executable="iqtree3", batch_mode=True,
    )

    assert result["status"] == "success"
    assert Path(result["output_tree"]).exists()
    assert Path(result["log_iqtree"]).exists()
    # batch mode: tree in trees/, iqtree in logs/
    assert "trees" in str(result["output_tree"])


@pytest.mark.skipif(not shutil.which("iqtree3"), reason="iqtree3 not found in PATH")
def test_run_one_iqtree_mf_mode_no_tree(tmp_path: Path) -> None:
    from phyloai.tree.ml_iqtree import _run_one_iqtree

    inp = tmp_path / "gene.fa"
    inp.write_text(">a\nMKTLLLMKT\n>b\nMKTLLLMKT\n>c\nMKTLLLMKT\n")

    out_dir = tmp_path / "out"
    out_dir.mkdir()

    result = _run_one_iqtree(
        gene_path=inp, seq_type="AA",
        model_string="MF", modelfinder="MF",
        boot=None, alrt=None, bnni=False,
        mode="normal", threads_arg="-T 2",
        log_dir=out_dir, output_dir=out_dir,
        executable="iqtree3", batch_mode=False,
    )

    assert result["status"] == "success"
    assert result["output_tree"] is None
    assert Path(result["log_iqtree"]).exists()


# ===================================================================
# run_iqtree — main entry point
# ===================================================================

def test_run_iqtree_neither_input_raises(tmp_path: Path) -> None:
    from phyloai.tree.ml_iqtree import run_iqtree

    with pytest.raises(ValueError, match="Either --msa-dir or --matrix"):
        run_iqtree(output_dir=tmp_path / "out", quiet=True)


def test_run_iqtree_both_inputs_raises(tmp_path: Path) -> None:
    from phyloai.tree.ml_iqtree import run_iqtree

    msa_dir = tmp_path / "msas"
    msa_dir.mkdir()
    mat = tmp_path / "matrix.fa"
    mat.write_text(">a\nMKT\n")

    with pytest.raises(ValueError, match="Either --msa-dir or --matrix"):
        run_iqtree(msa_dir=msa_dir, matrix=mat, output_dir=tmp_path / "out", quiet=True)


def test_run_iqtree_single_dry_run(tmp_path: Path) -> None:
    from phyloai.tree.ml_iqtree import run_iqtree

    mat = tmp_path / "matrix.fa"
    mat.write_text(">a\nMKTLLL\n>b\nMKTLLL\n")
    out_dir = tmp_path / "out"

    payload = run_iqtree(
        matrix=mat, output_dir=out_dir,
        seq_type="AA", model="LG", dry_run=True, quiet=True,
    )

    assert payload["status"] == "success"
    assert isinstance(payload["data"]["cmd"], list)
    assert "tool_stderr" in payload["data"]
    assert "output" in payload["data"]


def test_run_iqtree_batch_dry_run(tmp_path: Path) -> None:
    from phyloai.tree.ml_iqtree import run_iqtree

    msa_dir = tmp_path / "msas"
    msa_dir.mkdir()
    (msa_dir / "g1.fa").write_text(">a\nMKTLLL\n>b\nMKTLLL\n")
    (msa_dir / "g2.fa").write_text(">c\nMKTLLL\n>d\nMKTLLL\n")

    out_dir = tmp_path / "out"

    payload = run_iqtree(
        msa_dir=msa_dir, output_dir=out_dir,
        seq_type="AA", model="LG", dry_run=True, quiet=True,
    )

    assert payload["status"] == "success"
    assert payload["data"]["summary"]["n_input_files"] >= 2
    for f in payload["data"]["files"]:
        assert "cmd" in f


def test_run_iqtree_heterogeneous_rejected_in_batch(tmp_path: Path) -> None:
    from phyloai.tree.ml_iqtree import run_iqtree

    msa_dir = tmp_path / "msas"
    msa_dir.mkdir()
    (msa_dir / "g1.fa").write_text(">a\nMKTLLL\n")

    with pytest.raises(ValueError, match="only supported in --matrix"):
        run_iqtree(msa_dir=msa_dir, output_dir=tmp_path / "out",
                   seq_type="AA", model="C20", quiet=True)


def test_run_iqtree_no_valid_inputs(tmp_path: Path) -> None:
    from phyloai.tree.ml_iqtree import run_iqtree

    msa_dir = tmp_path / "empty"
    msa_dir.mkdir()

    with pytest.raises(ValueError, match="No valid input files"):
        run_iqtree(msa_dir=msa_dir, output_dir=tmp_path / "out", quiet=True)


def test_run_iqtree_unsupported_matrix_extension(tmp_path: Path) -> None:
    from phyloai.tree.ml_iqtree import run_iqtree

    bad = tmp_path / "matrix.xyz"
    bad.write_text("anything\n")

    with pytest.raises(ValueError, match="unsupported extension"):
        run_iqtree(matrix=bad, output_dir=tmp_path / "out", quiet=True)


def test_run_iqtree_unparsable_matrix(tmp_path: Path) -> None:
    from phyloai.tree.ml_iqtree import run_iqtree

    bad = tmp_path / "matrix.fa"
    bad.write_text("not a valid fasta file\n")

    with pytest.raises(ValueError, match="Cannot parse"):
        run_iqtree(matrix=bad, output_dir=tmp_path / "out", quiet=True)


def test_run_iqtree_overwrite_resume_mutual(tmp_path: Path) -> None:
    from phyloai.tree.ml_iqtree import run_iqtree

    mat = tmp_path / "matrix.fa"
    mat.write_text(">a\nMKT\n")

    with pytest.raises(ValueError, match="mutually exclusive"):
        run_iqtree(matrix=mat, output_dir=tmp_path / "out",
                   overwrite=True, resume=True, quiet=True)


def test_run_iqtree_nonempty_output_dir(tmp_path: Path) -> None:
    from phyloai.tree.ml_iqtree import run_iqtree

    mat = tmp_path / "matrix.fa"
    mat.write_text(">a\nMKT\n>b\nMKT\n")

    out_dir = tmp_path / "out"
    out_dir.mkdir()
    (out_dir / "junk.txt").write_text("old data")

    with pytest.raises(ValueError, match="already exists"):
        run_iqtree(matrix=mat, output_dir=out_dir, quiet=True)


def test_run_iqtree_overwrite_removes_dir(tmp_path: Path) -> None:
    from phyloai.tree.ml_iqtree import run_iqtree

    mat = tmp_path / "matrix.fa"
    mat.write_text(">a\nMKTLLL\n>b\nMKTLLL\n")

    out_dir = tmp_path / "out"
    out_dir.mkdir()
    (out_dir / "junk.txt").write_text("old data")

    payload = run_iqtree(
        matrix=mat, output_dir=out_dir,
        overwrite=True, dry_run=True, quiet=True,
    )
    assert payload["status"] == "success"


def test_run_iqtree_auto_seq_type_detection(tmp_path: Path) -> None:
    from phyloai.tree.ml_iqtree import run_iqtree

    mat = tmp_path / "matrix.fa"
    mat.write_text(">a\nMKTLLL\n>b\nMKTLLL\n")

    out_dir = tmp_path / "out"

    payload = run_iqtree(
        matrix=mat, output_dir=out_dir,
        seq_type="auto", dry_run=True, quiet=True,
    )
    assert payload["status"] == "success"


def test_run_iqtree_pmsf_default_base_model_direct_api(tmp_path: Path) -> None:
    from phyloai.tree.ml_iqtree import run_iqtree

    mat = tmp_path / "matrix.fa"
    mat.write_text(">a\nMKTLLL\n>b\nMKTLLL\n")
    guide = tmp_path / "guide.nwk"
    guide.write_text("(a,b);\n")

    payload = run_iqtree(
        matrix=mat,
        output_dir=tmp_path / "out",
        seq_type="AA",
        model="C20",
        guide_tree=str(guide),
        dry_run=True,
        quiet=True,
    )

    cmd = " ".join(payload["data"]["cmd"])
    assert "-m LG+C20+F+R4" in cmd
    assert payload["params"]["pmsf_base_model"] == "LG"


def test_run_iqtree_nt_seq_type_mapping(tmp_path: Path) -> None:
    from phyloai.tree.ml_iqtree import run_iqtree

    mat = tmp_path / "matrix.fa"
    mat.write_text(">a\nACGTACGT\n>b\nACGTACGT\n")

    out_dir = tmp_path / "out"

    payload = run_iqtree(
        matrix=mat, output_dir=out_dir,
        seq_type="NT", model="GTR", dry_run=True, quiet=True,
    )
    assert payload["status"] == "success"
    assert payload["key_results"]["seq_type"] == "NT"


def test_run_iqtree_seq_type_mismatch(tmp_path: Path) -> None:
    from phyloai.tree.ml_iqtree import run_iqtree

    mat = tmp_path / "matrix.fa"
    mat.write_text(">a\nMKTLLL\n>b\nMKTLLL\n")

    out_dir = tmp_path / "out"

    with pytest.raises(ValueError, match="but detected"):
        run_iqtree(matrix=mat, output_dir=out_dir,
                   seq_type="NT", quiet=True)


# ===================================================================
# Checkpoint resume verifier (IQ-TREE)
# ===================================================================

def test_resume_verifier_iqtree_validates_treefile(tmp_path: Path) -> None:
    from phyloai.tree.checkpoint_helpers import resume_verifier_iqtree

    tree = tmp_path / "gene.treefile"
    tree.write_text("(a:0.1,b:0.2);")

    verify = resume_verifier_iqtree(validate_tree=True)
    assert verify(tree) is True


def test_resume_verifier_iqtree_accepts_multiple_trees(tmp_path: Path) -> None:
    from phyloai.tree.checkpoint_helpers import resume_verifier_iqtree

    tree = tmp_path / "gene.treefile"
    tree.write_text("(a:0.1,b:0.2);\n(a:0.3,b:0.4);\n")

    assert resume_verifier_iqtree()(tree) is True


def test_resume_verifier_iqtree_rejects_empty(tmp_path: Path) -> None:
    from phyloai.tree.checkpoint_helpers import resume_verifier_iqtree

    tree = tmp_path / "gene.treefile"
    tree.write_text("")

    verify = resume_verifier_iqtree(validate_tree=True)
    assert verify(tree) is False


def test_resume_verifier_iqtree_rejects_unparsable(tmp_path: Path) -> None:
    from phyloai.tree.checkpoint_helpers import resume_verifier_iqtree

    tree = tmp_path / "gene.treefile"
    tree.write_text("(a,b")

    verify = resume_verifier_iqtree(validate_tree=True)
    assert verify(tree) is False


def test_resume_verifier_iqtree_missing_file(tmp_path: Path) -> None:
    from phyloai.tree.checkpoint_helpers import resume_verifier_iqtree

    verify = resume_verifier_iqtree()
    assert verify(tmp_path / "missing.treefile") is False


def test_resume_verifier_iqtree_mf_mode_validates_iqtree(tmp_path: Path) -> None:
    from phyloai.tree.checkpoint_helpers import resume_verifier_iqtree

    iqtree = tmp_path / "gene.iqtree"
    iqtree.write_text("Log-likelihood: -100.0\n")

    verify = resume_verifier_iqtree(validate_tree=False)
    assert verify(iqtree) is True


def test_resume_verifier_iqtree_mf_mode_empty_iqtree(tmp_path: Path) -> None:
    from phyloai.tree.checkpoint_helpers import resume_verifier_iqtree

    iqtree = tmp_path / "gene.iqtree"
    iqtree.write_text("")

    verify = resume_verifier_iqtree(validate_tree=False)
    assert verify(iqtree) is False


# ===================================================================
# #15: --boot defaults (1000 by default, 0 suppresses)
# ===================================================================


def test_build_cmd_boot_default_1000(tmp_path: Path) -> None:
    from phyloai.tree.ml_iqtree import _build_iqtree_cmd

    inp = tmp_path / "gene.fa"
    pref = tmp_path / "gene"
    cmd = _build_iqtree_cmd(
        input_path=inp, prefix=pref,
        model_string="LG+F+R4", seq_type="AA",
        modelfinder="none", boot=1000, alrt=None, bnni=False,
        mode="normal", threads_arg="-T AUTO",
    )
    flat = " ".join(cmd)
    assert "-B 1000" in flat


def test_build_cmd_boot_zero_suppresses(tmp_path: Path) -> None:
    from phyloai.tree.ml_iqtree import _build_iqtree_cmd

    inp = tmp_path / "gene.fa"
    pref = tmp_path / "gene"
    cmd = _build_iqtree_cmd(
        input_path=inp, prefix=pref,
        model_string="LG+F+R4", seq_type="AA",
        modelfinder="none", boot=0, alrt=None, bnni=False,
        mode="normal", threads_arg="-T AUTO",
    )
    assert "-B" not in cmd


def test_build_cmd_boot_none_omits(tmp_path: Path) -> None:
    from phyloai.tree.ml_iqtree import _build_iqtree_cmd

    inp = tmp_path / "gene.fa"
    pref = tmp_path / "gene"
    cmd = _build_iqtree_cmd(
        input_path=inp, prefix=pref,
        model_string="LG+F+R4", seq_type="AA",
        modelfinder="none", boot=None, alrt=None, bnni=False,
        mode="normal", threads_arg="-T AUTO",
    )
    assert "-B" not in cmd


# ===================================================================
# #14: --rclusterf 10 default when --partitions + MF/MFP
# ===================================================================


def test_build_cmd_partitions_mfp_rclusterf_default(tmp_path: Path) -> None:
    from phyloai.tree.ml_iqtree import _build_iqtree_cmd

    inp = tmp_path / "gene.fa"
    pref = tmp_path / "gene"
    cmd = _build_iqtree_cmd(
        input_path=inp, prefix=pref,
        model_string="MFP", seq_type="AA",
        boot=None, alrt=None, bnni=False,
        modelfinder="MFP", mode="normal", threads_arg="-T AUTO",
        partitions="/p.nex", rclusterf=10,
    )
    flat = " ".join(cmd)
    assert "--merge" in cmd
    assert "--rclusterf 10" in flat


def test_build_cmd_partitions_mf_rclusterf_default(tmp_path: Path) -> None:
    from phyloai.tree.ml_iqtree import _build_iqtree_cmd

    inp = tmp_path / "gene.fa"
    pref = tmp_path / "gene"
    cmd = _build_iqtree_cmd(
        input_path=inp, prefix=pref,
        model_string="MF", seq_type="AA",
        boot=None, alrt=None, bnni=False,
        modelfinder="MF", mode="normal", threads_arg="-T AUTO",
        partitions="/p.nex", rclusterf=10,
    )
    flat = " ".join(cmd)
    assert "--rclusterf 10" in flat


# ===================================================================
# #16: --keep-extra behavior
# ===================================================================


def test_run_one_iqtree_keep_extra(tmp_path: Path) -> None:
    from unittest.mock import patch
    from phyloai.tree.ml_iqtree import _run_one_iqtree

    inp = tmp_path / "gene.fa"
    inp.write_text(">a\nMKTLLL\n>b\nMKTLLL\n")
    out_dir = tmp_path / "out"
    out_dir.mkdir()
    log_dir = tmp_path / "logs"
    log_dir.mkdir()
    work_dir = tmp_path / "work"
    work_dir.mkdir()

    def _fake_run(cmd, **kw):
        iqtree = work_dir / "gene.iqtree"
        iqtree.write_text("Log-likelihood: -100.0\n")
        (work_dir / "gene.treefile").write_text("(a:0.1,b:0.1);\n")
        (work_dir / "gene.log").write_text("done\n")
        (work_dir / "gene.ufboot").write_text("boot data\n")
        return type("R", (), {"returncode": 0, "stdout": "", "stderr": ""})()

    with patch("subprocess.run", _fake_run):
        result = _run_one_iqtree(
            gene_path=inp, seq_type="AA",
            model_string="LG+F+R4", modelfinder="none",
            boot=1000, alrt=None, bnni=False,
            mode="normal", threads_arg="-T 1",
            log_dir=log_dir, output_dir=out_dir,
            keep_extra=True, batch_mode=True, work_dir=work_dir,
        )
    assert result["status"] == "success"
    assert (log_dir / "gene.ufboot").exists()


def test_run_one_iqtree_accepts_multiple_trees(tmp_path: Path) -> None:
    from unittest.mock import patch
    from phyloai.tree.ml_iqtree import _run_one_iqtree

    inp = tmp_path / "gene.fa"
    inp.write_text(">a\nMKTLLL\n>b\nMKTLLL\n")
    out_dir = tmp_path / "out"
    log_dir = tmp_path / "logs"
    work_dir = tmp_path / "work"
    out_dir.mkdir()
    log_dir.mkdir()
    work_dir.mkdir()

    def _fake_run(cmd, **kw):
        (work_dir / "gene.iqtree").write_text("Log-likelihood: -100.0\n")
        (work_dir / "gene.treefile").write_text("(a:0.1,b:0.1);\n(a:0.2,b:0.2);\n")
        (work_dir / "gene.log").write_text("done\n")
        return type("R", (), {"returncode": 0, "stdout": "", "stderr": ""})()

    with patch("subprocess.run", _fake_run):
        result = _run_one_iqtree(
            gene_path=inp, seq_type="AA",
            model_string="LG+H8", modelfinder="none",
            boot=None, alrt=None, bnni=False,
            mode="normal", threads_arg="-T 1",
            log_dir=log_dir, output_dir=out_dir,
            batch_mode=True, work_dir=work_dir,
        )

    assert result["status"] == "success"


def test_run_one_iqtree_no_keep_extra(tmp_path: Path) -> None:
    from unittest.mock import patch
    from phyloai.tree.ml_iqtree import _run_one_iqtree

    inp = tmp_path / "gene.fa"
    inp.write_text(">a\nMKTLLL\n>b\nMKTLLL\n")
    out_dir = tmp_path / "out"
    log_dir = tmp_path / "logs"
    work_dir = tmp_path / "work"
    out_dir.mkdir()
    log_dir.mkdir()
    work_dir.mkdir()

    def _fake_run(cmd, **kw):
        (work_dir / "gene.iqtree").write_text("Log-likelihood: -100.0\n")
        (work_dir / "gene.treefile").write_text("(a:0.1,b:0.1);\n")
        (work_dir / "gene.log").write_text("done\n")
        (work_dir / "gene.ufboot").write_text("boot data\n")
        return type("R", (), {"returncode": 0, "stdout": "", "stderr": ""})()

    with patch("subprocess.run", _fake_run):
        result = _run_one_iqtree(
            gene_path=inp, seq_type="AA",
            model_string="LG+F+R4", modelfinder="none",
            boot=1000, alrt=None, bnni=False,
            mode="normal", threads_arg="-T 1",
            log_dir=log_dir, output_dir=out_dir,
            keep_extra=False, batch_mode=True, work_dir=work_dir,
        )

    assert result["status"] == "success"
    assert not (log_dir / "gene.ufboot").exists()


# ===================================================================
# #7: Duplicate stem detection in batch mode
# ===================================================================


def test_run_iqtree_duplicate_stems_raises(tmp_path: Path) -> None:
    from phyloai.tree.ml_iqtree import run_iqtree

    msa_dir = tmp_path / "msas"
    msa_dir.mkdir()
    (msa_dir / "gene.fa").write_text(">a\nMKT\n>b\nMKT\n")
    (msa_dir / "gene.phy").write_text("2 3\na MKT\nb MKT\n")

    with pytest.raises(ValueError, match="Duplicate output stems"):
        run_iqtree(msa_dir=msa_dir, output_dir=tmp_path / "out", quiet=True)


def test_run_iqtree_unique_stems_passes(tmp_path: Path) -> None:
    from phyloai.tree.ml_iqtree import run_iqtree

    msa_dir = tmp_path / "msas"
    msa_dir.mkdir()
    (msa_dir / "gene1.fa").write_text(">a\nMKT\n>b\nMKT\n")
    (msa_dir / "gene2.fa").write_text(">a\nMKT\n>b\nMKT\n")

    result = run_iqtree(
        msa_dir=msa_dir, output_dir=tmp_path / "out",
        dry_run=True, quiet=True,
    )
    assert result["key_results"]["n_input"] == 2


# ===================================================================
# #9: _format_offender handles both reason and detected/expected
# ===================================================================


def test_format_offender_with_reason() -> None:
    from phyloai.tree.ml_iqtree import _format_offender

    o = {"file": "bad.fa", "reason": "failed to parse input file"}
    assert "bad.fa" in _format_offender(o)
    assert "failed to parse" in _format_offender(o)


def test_format_offender_with_detected_expected() -> None:
    from phyloai.tree.ml_iqtree import _format_offender

    o = {"file": "gene.fa", "detected": "NT", "expected": "AA"}
    s = _format_offender(o)
    assert "gene.fa" in s
    assert "NT" in s
    assert "AA" in s


# ===================================================================
# #13: Single-mode resume via IQ-TREE native checkpoint
# ===================================================================


def test_run_iqtree_single_resume_no_checkpoint(tmp_path: Path) -> None:
    from phyloai.tree.ml_iqtree import run_iqtree

    mat = tmp_path / "matrix.fa"
    mat.write_text(">a\nMKTLLL\n>b\nMKTLLL\n")
    out = tmp_path / "out"
    out.mkdir()

    result = run_iqtree(
        matrix=mat, output_dir=out,
        resume=True, dry_run=True, quiet=True,
    )
    assert result["status"] == "success"


def test_run_iqtree_single_resume_skips_dir_check(tmp_path: Path) -> None:
    from phyloai.tree.ml_iqtree import run_iqtree

    mat = tmp_path / "matrix.fa"
    mat.write_text(">a\nMKTLLL\n>b\nMKTLLL\n")
    out = tmp_path / "out"
    out.mkdir()
    (out / "existing.txt").write_text("old")

    result = run_iqtree(
        matrix=mat, output_dir=out,
        resume=True, dry_run=True, quiet=True,
    )
    assert result["status"] == "success"


# ===================================================================
# #12: -qmax override detection
# ===================================================================


def test_build_cmd_qmax_overridable_via_single_dash(tmp_path: Path) -> None:
    from phyloai.tree.ml_iqtree import _build_iqtree_cmd

    inp = tmp_path / "gene.fa"
    pref = tmp_path / "gene"
    cmd = _build_iqtree_cmd(
        input_path=inp, prefix=pref,
        model_string="MIX+MF", seq_type="NT",
        boot=None, alrt=None, bnni=False,
        modelfinder="none", mode="normal", threads_arg="-T AUTO",
        qmax=10,
        tool_args="-qmax 20",
    )
    flat = " ".join(cmd)
    assert "-qmax 20" in flat
    assert "-qmax 10" not in flat


# ===================================================================
# #11: --quiet suppresses warnings
# ===================================================================


def test_run_validations_quiet_suppresses_warnings() -> None:
    import warnings
    from phyloai.tree.ml_iqtree import _run_validations

    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        _run_validations(
            batch_mode=False, seq_type="AA",
            modelfinder="MF", model="LG",
            partitions=None, guide_tree=None,
            boot=1000, alrt=None, bnni=False,
            rclusterf=5, rcluster_max=None,
            quiet=True,
        )
        assert len(w) == 0


def test_run_validations_non_quiet_emits_warnings() -> None:
    import warnings
    from phyloai.tree.ml_iqtree import _run_validations

    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        _run_validations(
            batch_mode=False, seq_type="AA",
            modelfinder="MF", model="LG",
            partitions=None, guide_tree=None,
            boot=1000, alrt=None, bnni=False,
            rclusterf=5, rcluster_max=None,
            quiet=False,
        )
        assert len(w) >= 1
