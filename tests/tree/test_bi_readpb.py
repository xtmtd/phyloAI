from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from phyloai.tree.bi_readpb import (
    _convert_meanrr_to_exchangeabilities,
    _convert_siteprofiles_to_sitefreq,
    _write_pmsf_partition,
    _validate_modes,
    run_bi_readpb,
    PAML_ORDER,
)


def test_write_pmsf_partition_uses_posterior_rates_alpha_and_gamma_categories(tmp_path: Path):
    sitefreq = tmp_path / "chain1.sitefreq"
    frequencies = " ".join(["0.05000000"] * 20)
    sitefreq.write_text(f"2 {frequencies}\n1 {frequencies}\n")
    exchangeabilities = tmp_path / "chain1.exchangeabilities"
    exchangeabilities.write_text("custom model\n")
    meansiterates = tmp_path / "chain1.meansiterates"
    meansiterates.write_text("1\t0.82152\n0\t1.42275\n")
    trace = tmp_path / "chain1.trace"
    trace.write_text("iter\talpha\n0\t1.0\n1\t2.0\n")
    log = tmp_path / "chain1.log"
    log.write_text("discrete gamma distribution of rates across sites (4 categories)\n")

    result, mean_alpha, categories = _write_pmsf_partition(
        sitefreq, exchangeabilities, meansiterates, trace, log, 0, 1, "all",
    )

    text = result.read_text()
    assert result.name == "partition.PMSF.nex"
    assert "charset site_1 = 1;" in text
    assert text.index("charset site_1") < text.index("charset site_2")
    assert "chain1.exchangeabilities+F{" in text
    assert "+G4{1.50000000}:site_1{1.42275}" in text
    assert "+G4{1.50000000}:site_2{0.82152};" in text
    assert mean_alpha == 1.5
    assert categories == 4


@pytest.mark.parametrize(
    ("meansiterates_rows", "error"),
    [
        (["0\t1.0", "0\t1.0"], "Duplicate meansiterates"),
        (["0\tbad", "1\t1.0"], "Invalid meansiterates rate"),
        (["-1\t1.0", "1\t1.0"], "Invalid meansiterates site index"),
        (["0\t0", "1\t1.0"], "Invalid meansiterates rate"),
        (["0\t1.0"], "site IDs do not match"),
        (["0\t1.0", "1\t1.0", "2\t1.0"], "site IDs do not match"),
    ],
)
def test_write_pmsf_partition_rejects_invalid_meansiterates_records(
    tmp_path: Path, meansiterates_rows: list[str], error: str,
):
    sitefreq = tmp_path / "chain1.sitefreq"
    frequencies = " ".join(["0.05"] * 20)
    sitefreq.write_text(f"1 {frequencies}\n2 {frequencies}\n")
    exchangeabilities = tmp_path / "chain1.exchangeabilities"
    exchangeabilities.write_text("model\n")
    meansiterates = tmp_path / "chain1.meansiterates"
    meansiterates.write_text("\n".join(meansiterates_rows) + "\n")
    trace = tmp_path / "chain1.trace"
    trace.write_text("iter\talpha\n0\t1.0\n")
    log = tmp_path / "chain1.log"
    log.write_text("discrete gamma distribution of rates across sites (4 categories)\n")

    with pytest.raises(ValueError, match=error):
        _write_pmsf_partition(sitefreq, exchangeabilities, meansiterates, trace, log, 0, 1, "all")
    assert not (tmp_path / "partition.PMSF.nex").exists()


def test_write_pmsf_partition_uses_requested_trace_sampling_window(tmp_path: Path):
    sitefreq = tmp_path / "chain1.sitefreq"
    frequencies = " ".join(["0.05"] * 20)
    sitefreq.write_text(f"1 {frequencies}\n")
    exchangeabilities = tmp_path / "chain1.exchangeabilities"
    exchangeabilities.write_text("model\n")
    meansiterates = tmp_path / "chain1.meansiterates"
    meansiterates.write_text("0\t1.0\n")
    trace = tmp_path / "chain1.trace"
    trace.write_text("iter\talpha\n0\t1.0\n1\t2.0\n2\t4.0\n3\t8.0\n")
    log = tmp_path / "chain1.log"
    log.write_text("discrete gamma distribution of rates across sites (4 categories)\n")

    result, mean_alpha, _ = _write_pmsf_partition(
        sitefreq, exchangeabilities, meansiterates, trace, log, 1, 2, "3",
    )

    assert mean_alpha == 5.0
    assert "+G4{5.00000000}" in result.read_text()


def test_write_pmsf_partition_rejects_empty_sitefreq(tmp_path: Path):
    sitefreq = tmp_path / "chain1.sitefreq"
    sitefreq.write_text("")
    exchangeabilities = tmp_path / "chain1.exchangeabilities"
    exchangeabilities.write_text("model\n")
    meansiterates = tmp_path / "chain1.meansiterates"
    meansiterates.write_text("0\t1.0\n")
    trace = tmp_path / "chain1.trace"
    trace.write_text("iter\talpha\n0\t1.0\n")
    log = tmp_path / "chain1.log"
    log.write_text("discrete gamma distribution of rates across sites (4 categories)\n")

    with pytest.raises(ValueError, match="no sites"):
        _write_pmsf_partition(sitefreq, exchangeabilities, meansiterates, trace, log, 0, 1, "all")


def test_write_pmsf_partition_rejects_nonpositive_sitefreq_id(tmp_path: Path):
    sitefreq = tmp_path / "chain1.sitefreq"
    sitefreq.write_text("0 " + " ".join(["0.05"] * 20) + "\n")
    exchangeabilities = tmp_path / "chain1.exchangeabilities"
    exchangeabilities.write_text("model\n")
    meansiterates = tmp_path / "chain1.meansiterates"
    meansiterates.write_text("0\t1.0\n")
    trace = tmp_path / "chain1.trace"
    trace.write_text("iter\talpha\n0\t1.0\n")
    log = tmp_path / "chain1.log"
    log.write_text("discrete gamma distribution of rates across sites (4 categories)\n")

    with pytest.raises(ValueError, match="Invalid sitefreq site ID"):
        _write_pmsf_partition(sitefreq, exchangeabilities, meansiterates, trace, log, 0, 1, "all")


def test_write_pmsf_partition_reports_missing_metadata_as_validation_error(tmp_path: Path):
    sitefreq = tmp_path / "chain1.sitefreq"
    sitefreq.write_text("1 " + " ".join(["0.05"] * 20) + "\n")
    exchangeabilities = tmp_path / "chain1.exchangeabilities"
    exchangeabilities.write_text("model\n")
    meansiterates = tmp_path / "chain1.meansiterates"
    meansiterates.write_text("0\t1.0\n")

    with pytest.raises(ValueError, match="Could not read PMSF metadata"):
        _write_pmsf_partition(
            sitefreq, exchangeabilities, meansiterates,
            tmp_path / "chain1.trace", tmp_path / "chain1.log", 0, 1, "all",
        )


def test_validate_modes_unrecognized():
    with pytest.raises(ValueError, match="Unrecognised"):
        _validate_modes(["bad_mode"])


def test_validate_modes_duplicates():
    with pytest.raises(ValueError, match="Duplicate"):
        _validate_modes(["rr", "rr"])


def test_validate_modes_allppred_conflict():
    with pytest.raises(ValueError, match="mutually exclusive"):
        _validate_modes(["allppred", "div"])
    with pytest.raises(ValueError, match="mutually exclusive"):
        _validate_modes(["allppred", "sitecomp", "rr"])


def test_validate_modes_happy():
    _validate_modes(["rr", "ss", "r"])
    _validate_modes(["allppred"])


def test_convert_meanrr_to_exchangeabilities(tmp_path: Path):
    header = "A C D E F G H I K L M N P Q R S T V W Y"
    lines = [header]
    # Add a few symmetric exchangeability values
    for (src, tgt, val) in [("A", "R", "1.5"), ("A", "N", "0.8"), ("C", "D", "2.1")]:
        lines.append(f"{src} {tgt} {val}")
    for (tgt, src, val) in [("R", "A", "1.5"), ("N", "A", "0.8"), ("D", "C", "2.1")]:
        lines.append(f"{src} {tgt} {val}")

    meanrr = tmp_path / "chain1.meanrr"
    meanrr.write_text("\n".join(lines))

    result = _convert_meanrr_to_exchangeabilities(meanrr)
    assert result.name == "chain1.exchangeabilities"
    assert result.exists()

    content = result.read_text()
    content_lines = content.splitlines()

    # Row 0 is a blank line (no lower-triangle values for A)
    assert content_lines[0].strip() == ""
    # Last non-empty line should be the uniform prior
    assert "0.050000" in content_lines[-1]
    # Check PAML order: A,R,N,D,C,Q,E,G,H,I,L,K,M,F,P,S,T,W,Y,V
    idx_r = PAML_ORDER.index("R")
    row_r = content_lines[idx_r].strip()
    values_r = [float(x) for x in row_r.split()]
    assert len(values_r) == 1
    assert abs(values_r[0] - 1.5) < 0.001


def test_convert_siteprofiles_to_sitefreq(tmp_path: Path):
    header1 = "# PhyloBayes siteprofiles"
    header2 = "# Site-specific frequency profiles"
    pb_order = "A C D E F G H I K L M N P Q R S T V W Y"
    # Create one site with non-uniform frequencies
    freqs = [0.05] * 20
    freqs[0] = 0.5  # A is high
    freq_str = " ".join(f"{f:.6f}" for f in freqs)
    data_line = f"1 {freq_str}"

    siteprof = tmp_path / "chain1.siteprofiles"
    siteprof.write_text("\n".join([header1, header2, data_line]))

    result = _convert_siteprofiles_to_sitefreq(siteprof)
    assert result.name == "chain1.sitefreq"
    assert result.exists()

    content = result.read_text().strip()
    parts = content.split()
    assert parts[0] == "1"
    reindexed = [float(x) for x in parts[1:]]
    assert len(reindexed) == 20
    assert abs(sum(reindexed) - 1.0) < 0.001
    # Values should be >= 1e-8
    for v in reindexed:
        assert v >= 1e-8


@pytest.mark.parametrize(
    "data_line",
    [
        "1 " + " ".join(["0.05"] * 19),
        "1 " + " ".join(["0.05"] * 20) + " extra",
        "1 bad " + " ".join(["0.05"] * 19),
    ],
)
def test_convert_siteprofiles_rejects_malformed_data_rows(tmp_path: Path, data_line: str):
    siteprof = tmp_path / "chain1.siteprofiles"
    siteprof.write_text("# header\n# header\n" + data_line + "\n")

    with pytest.raises(ValueError, match="Invalid siteprofiles row 3"):
        _convert_siteprofiles_to_sitefreq(siteprof)


def test_run_bi_readpb_dry_run(tmp_path: Path):
    chain = tmp_path / "chain1"
    (tmp_path / "chain1.chain").write_text("fake state")
    out = tmp_path / "readpb"

    payload = run_bi_readpb(
        chain=chain, mode="ss,rr", output_dir=out,
        burnin=1000, dry_run=True,
    )
    assert payload["status"] == "success"
    assert payload["key_results"]["modes_run"] == ["ss", "rr"]
    assert "ss" in payload["data"]["cmds"]
    assert "rr" in payload["data"]["cmds"]
    assert "-x" in payload["data"]["cmds"]["ss"]


def test_run_bi_readpb_missing_chain(tmp_path: Path):
    chain = tmp_path / "chain1"
    with pytest.raises(ValueError, match="not found"):
        run_bi_readpb(chain=chain, mode="rr")


def test_run_bi_readpb_ss_rr_roundtrip(tmp_path: Path):
    chain_dir = tmp_path / "chains"
    chain_dir.mkdir()
    chain = chain_dir / "chain1"
    (chain_dir / "chain1.chain").write_text("fake state")
    (chain_dir / "chain1.trace").write_text("iter\talpha\n0\t1.0\n")
    (chain_dir / "chain1.treelist").write_text("(a,b);\n")

    tool_dir = tmp_path / "tools"
    tool_dir.mkdir()
    mpirun = tool_dir / "mpirun"
    mpirun.write_text("#!/bin/sh\nshift 2\nexec \"$@\"\n")
    readpb = tool_dir / "readpb_mpi"
    readpb.write_text("""#!/bin/sh
mode=""
chain_name=""
while [ $# -gt 0 ]; do
    case "$1" in
        -ss) mode="ss";;
        -rr) mode="rr";;
        -r) mode="r";;
        -*) ;;
        *) chain_name="$1";;
    esac
    shift
done
if [ "$mode" = "ss" ]; then
    printf '# PhyloBayes siteprofiles\\n# Site profiles\\n' > "${chain_name}.siteprofiles"
    AA="A C D E F G H I K L M N P Q R S T V W Y"
    freqs=$(printf '0.050000 %.0s' {1..20})
    printf '1 %s\\n' "$freqs" >> "${chain_name}.siteprofiles"
elif [ "$mode" = "rr" ]; then
    printf 'A C D E F G H I K L M N P Q R S T V W Y\\nA R 1.5\\nA N 0.8\\nC D 2.1\\n' > "${chain_name}.meanrr"
fi
""")
    mpirun.chmod(0o755)
    readpb.chmod(0o755)

    out = tmp_path / "readpb"
    payload = run_bi_readpb(
        chain=chain, mode="ss,rr", output_dir=out,
        burnin=1000, threads=4,
        pb_path=tool_dir, quiet=True,
    )
    assert payload["status"] == "success"
    assert payload["key_results"]["modes_run"] == ["ss", "rr"]
    assert "chain1_siteprofiles" in payload["key_results"]["output_files"]
    assert "chain1_meanrr" in payload["key_results"]["output_files"]
    assert "sitefreq" in payload["key_results"]["output_files"]
    assert "exchangeabilities" in payload["key_results"]["output_files"]
    assert (out / "chain1.siteprofiles").exists()
    assert (out / "chain1.sitefreq").exists()
    assert (out / "chain1.meanrr").exists()
    assert (out / "chain1.exchangeabilities").exists()
    assert not (out / "partition.PMSF.nex").exists()
    # Check post-processing status
    pp = payload["data"]["post_processing"]
    assert pp["ss"]["status"] == "success"
    assert pp["rr"]["status"] == "success"


def test_run_bi_readpb_moves_each_mode_before_next_starts(tmp_path: Path):
    chain_dir = tmp_path / "chains"
    chain_dir.mkdir()
    chain = chain_dir / "chain1"
    (chain_dir / "chain1.chain").write_text("fake state")
    out = tmp_path / "readpb"

    tool_dir = tmp_path / "tools"
    tool_dir.mkdir()
    mpirun = tool_dir / "mpirun"
    mpirun.write_text("#!/bin/sh\nshift 2\nexec \"$@\"\n")
    readpb = tool_dir / "readpb_mpi"
    readpb.write_text(f"""#!/bin/sh
set -e
mode=""
chain_name=""
while [ $# -gt 0 ]; do
    case "$1" in
        -rr) mode="rr";;
        -ss) mode="ss";;
        -ppred) mode="ppred";;
        -*) ;;
        *) chain_name="$1";;
    esac
    shift
done
if [ "$mode" = "rr" ]; then
    printf 'A C D E F G H I K L M N P Q R S T V W Y\\n' > "${{chain_name}}.meanrr"
elif [ "$mode" = "ss" ]; then
    test -s "{out}/chain1.meanrr"
    test -s "{out}/chain1.exchangeabilities"
    printf '# header\\n# header\\n1 ' > "${{chain_name}}.siteprofiles"
    printf '0.050000 %.0s' {{1..20}} >> "${{chain_name}}.siteprofiles"
    printf '\\n' >> "${{chain_name}}.siteprofiles"
elif [ "$mode" = "ppred" ]; then
    test -s "{out}/chain1.siteprofiles"
    test -s "{out}/chain1.sitefreq"
    printf 'replicate\\n' > "${{chain_name}}.ppred.1.ali"
fi
""")
    mpirun.chmod(0o755)
    readpb.chmod(0o755)

    payload = run_bi_readpb(
        chain=chain, mode="rr,ss,ppred", output_dir=out,
        pb_path=tool_dir, quiet=True,
    )

    assert payload["status"] == "success"
    assert (out / "chain1.meanrr").exists()
    assert (out / "chain1.exchangeabilities").exists()
    assert (out / "chain1.siteprofiles").exists()
    assert (out / "chain1.sitefreq").exists()
    assert (out / "ppred" / "chain1.ppred.1.ali").exists()
    assert not (chain_dir / "chain1.meanrr").exists()
    assert not (chain_dir / "chain1.siteprofiles").exists()
    assert not (chain_dir / "chain1.ppred.1.ali").exists()


def test_run_bi_readpb_writes_pmsf_partition_before_later_modes(tmp_path: Path):
    chain_dir = tmp_path / "chains"
    chain_dir.mkdir()
    chain = chain_dir / "chain1"
    (chain_dir / "chain1.chain").write_text("fake state")
    (chain_dir / "chain1.trace").write_text("iter\talpha\n0\t1.0\n")
    (chain_dir / "chain1.log").write_text("discrete gamma distribution of rates across sites (4 categories)\n")
    out = tmp_path / "readpb"

    tool_dir = tmp_path / "tools"
    tool_dir.mkdir()
    mpirun = tool_dir / "mpirun"
    mpirun.write_text("#!/bin/sh\nshift 2\nexec \"$@\"\n")
    readpb = tool_dir / "readpb_mpi"
    readpb.write_text(f"""#!/bin/sh
set -e
mode=""
chain_name=""
while [ $# -gt 0 ]; do
    case "$1" in
        -ss) mode="ss";;
        -rr) mode="rr";;
        -r) mode="r";;
        -*) ;;
        *) chain_name="$1";;
    esac
    shift
done
if [ "$mode" = "ss" ]; then
    printf '# header\\n# header\\n1 ' > "${{chain_name}}.siteprofiles"
    printf '0.050000 %.0s' {{1..20}} >> "${{chain_name}}.siteprofiles"
    printf '\\n' >> "${{chain_name}}.siteprofiles"
elif [ "$mode" = "rr" ]; then
    printf 'A C D E F G H I K L M N P Q R S T V W Y\\n' > "${{chain_name}}.meanrr"
elif [ "$mode" = "r" ]; then
    printf '0\\t1.0\\n' > "${{chain_name}}.meansiterates"
fi
""")
    mpirun.chmod(0o755)
    readpb.chmod(0o755)

    payload = run_bi_readpb(
        chain=chain, mode="ss,rr,r", output_dir=out,
        pb_path=tool_dir, quiet=True,
    )

    assert payload["status"] == "success"
    assert (out / "partition.PMSF.nex").exists()


def test_run_bi_readpb_does_not_write_pmsf_partition_when_r_fails(tmp_path: Path):
    chain_dir = tmp_path / "chains"
    chain_dir.mkdir()
    chain = chain_dir / "chain1"
    (chain_dir / "chain1.chain").write_text("fake state")
    (chain_dir / "chain1.trace").write_text("iter\talpha\n0\t1.0\n")
    (chain_dir / "chain1.log").write_text("discrete gamma distribution of rates across sites (4 categories)\n")
    out = tmp_path / "readpb"

    tool_dir = tmp_path / "tools"
    tool_dir.mkdir()
    mpirun = tool_dir / "mpirun"
    mpirun.write_text("#!/bin/sh\nshift 2\nexec \"$@\"\n")
    readpb = tool_dir / "readpb_mpi"
    readpb.write_text("""#!/bin/sh
mode=""
chain_name=""
while [ $# -gt 0 ]; do
    case "$1" in
        -ss) mode="ss";;
        -rr) mode="rr";;
        -r) mode="r";;
        -*) ;;
        *) chain_name="$1";;
    esac
    shift
done
if [ "$mode" = "ss" ]; then
    printf '# header\n# header\n1 ' > "${chain_name}.siteprofiles"
    printf '0.050000 %.0s' {1..20} >> "${chain_name}.siteprofiles"
    printf '\n' >> "${chain_name}.siteprofiles"
elif [ "$mode" = "rr" ]; then
    printf 'A C D E F G H I K L M N P Q R S T V W Y\n' > "${chain_name}.meanrr"
elif [ "$mode" = "r" ]; then
    exit 7
fi
""")
    mpirun.chmod(0o755)
    readpb.chmod(0o755)

    payload = run_bi_readpb(
        chain=chain, mode="ss,rr,r", output_dir=out,
        pb_path=tool_dir, quiet=True,
    )

    assert payload["status"] == "error"
    assert not (out / "partition.PMSF.nex").exists()


def test_run_bi_readpb_reports_new_rate_output(tmp_path: Path):
    chain_dir = tmp_path / "chains"
    chain_dir.mkdir()
    chain = chain_dir / "chain1"
    (chain_dir / "chain1.chain").write_text("fake state")

    tool_dir = tmp_path / "tools"
    tool_dir.mkdir()
    mpirun = tool_dir / "mpirun"
    mpirun.write_text("#!/bin/sh\nshift 2\nexec \"$@\"\n")
    readpb = tool_dir / "readpb_mpi"
    readpb.write_text("#!/bin/sh\nprintf 'rate\\n' > chain1.meansiterates\n")
    mpirun.chmod(0o755)
    readpb.chmod(0o755)

    payload = run_bi_readpb(
        chain=chain, mode="r", output_dir=tmp_path / "readpb",
        pb_path=tool_dir, quiet=True,
    )

    assert payload["status"] == "success"
    assert "chain1_meansiterates" in payload["data"]["output_files"]


def test_run_bi_readpb_moves_sitelogl_auxiliary_cpo(tmp_path: Path):
    chain_dir = tmp_path / "chains"
    chain_dir.mkdir()
    chain = chain_dir / "chain1"
    (chain_dir / "chain1.chain").write_text("fake state")

    tool_dir = tmp_path / "tools"
    tool_dir.mkdir()
    mpirun = tool_dir / "mpirun"
    mpirun.write_text("#!/bin/sh\nshift 2\nexec \"$@\"\n")
    readpb = tool_dir / "readpb_mpi"
    readpb.write_text(
        "#!/bin/sh\nprintf 'logl\\n' > chain1.sitelogl\nprintf 'cpo\\n' > chain1.cpo\n"
    )
    mpirun.chmod(0o755)
    readpb.chmod(0o755)

    out = tmp_path / "readpb"
    payload = run_bi_readpb(
        chain=chain, mode="sitelogl", output_dir=out,
        pb_path=tool_dir, quiet=True,
    )

    assert payload["status"] == "success"
    assert (out / "chain1.sitelogl").exists()
    assert (out / "chain1.cpo").exists()
    assert not (chain_dir / "chain1.cpo").exists()


def test_run_bi_readpb_moves_allppred_outputs(tmp_path: Path):
    chain_dir = tmp_path / "chains"
    chain_dir.mkdir()
    chain = chain_dir / "chain1"
    (chain_dir / "chain1.chain").write_text("fake state")

    tool_dir = tmp_path / "tools"
    tool_dir.mkdir()
    mpirun = tool_dir / "mpirun"
    mpirun.write_text("#!/bin/sh\nshift 2\nexec \"$@\"\n")
    readpb = tool_dir / "readpb_mpi"
    readpb.write_text(
        "#!/bin/sh\nprintf 'ppred\\n' > chain1.ppred\n"
    )
    mpirun.chmod(0o755)
    readpb.chmod(0o755)

    out = tmp_path / "readpb"
    payload = run_bi_readpb(
        chain=chain, mode="allppred", output_dir=out,
        pb_path=tool_dir, quiet=True,
    )

    assert payload["status"] == "success"
    assert (out / "chain1.ppred").exists()
    assert payload["data"]["tool_stderr"] == ""
    assert not (chain_dir / "chain1.ppred").exists()


def test_run_bi_readpb_rejects_stale_output_directory(tmp_path: Path):
    chain_dir = tmp_path / "chains"
    chain_dir.mkdir()
    chain = chain_dir / "chain1"
    (chain_dir / "chain1.chain").write_text("fake state")
    (chain_dir / "chain1.div").mkdir()

    tool_dir = tmp_path / "tools"
    tool_dir.mkdir()
    mpirun = tool_dir / "mpirun"
    mpirun.write_text("#!/bin/sh\nshift 2\nexec \"$@\"\n")
    readpb = tool_dir / "readpb_mpi"
    readpb.write_text("#!/bin/sh\nexit 0\n")
    mpirun.chmod(0o755)
    readpb.chmod(0o755)

    payload = run_bi_readpb(
        chain=chain, mode="div", output_dir=tmp_path / "readpb",
        pb_path=tool_dir, quiet=True,
    )

    assert payload["status"] == "error"
    assert "div" in payload["error"]


def test_run_bi_readpb_reports_stat_failures_as_warnings(tmp_path: Path, monkeypatch):
    chain_dir = tmp_path / "chains"
    chain_dir.mkdir()
    chain = chain_dir / "chain1"
    (chain_dir / "chain1.chain").write_text("fake state")

    tool_dir = tmp_path / "tools"
    tool_dir.mkdir()
    mpirun = tool_dir / "mpirun"
    mpirun.write_text("#!/bin/sh\nshift 2\nexec \"$@\"\n")
    readpb = tool_dir / "readpb_mpi"
    readpb.write_text("#!/bin/sh\nprintf 'result\\n' > chain1.div\n")
    mpirun.chmod(0o755)
    readpb.chmod(0o755)

    original_stat = Path.stat

    def failing_stat(path: Path, *args, **kwargs):
        if path.name == "chain1.div":
            raise OSError("simulated concurrent removal")
        return original_stat(path, *args, **kwargs)

    monkeypatch.setattr(Path, "stat", failing_stat)
    payload = run_bi_readpb(
        chain=chain, mode="div", output_dir=tmp_path / "readpb",
        pb_path=tool_dir, quiet=True,
    )

    assert payload["status"] == "error"
    assert "div" in payload["error"]
