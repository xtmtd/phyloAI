import pytest
from pathlib import Path
from phyloai.core.formats import FormatConverter, AlignmentFormat

FIXTURES = Path(__file__).parent / "fixtures"


def test_detect_fasta():
    fc = FormatConverter()
    fmt = fc.detect(FIXTURES / "test.fasta")
    assert fmt == AlignmentFormat.FASTA


def test_detect_phylip():
    fc = FormatConverter()
    fmt = fc.detect(FIXTURES / "test.phy")
    assert fmt == AlignmentFormat.PHYLIP


def test_detect_nexus():
    fc = FormatConverter()
    fmt = fc.detect(FIXTURES / "test.nex")
    assert fmt == AlignmentFormat.NEXUS


def test_fasta_to_phylip(tmp_path):
    fc = FormatConverter()
    out = tmp_path / "out.phy"
    fc.convert(FIXTURES / "test.fasta", out, target=AlignmentFormat.PHYLIP)
    assert out.exists()
    content = out.read_text()
    assert "Taxon_A" in content
    assert "MARVELLOUS" in content


def test_fasta_to_nexus(tmp_path):
    fc = FormatConverter()
    out = tmp_path / "out.nex"
    fc.convert(FIXTURES / "test.fasta", out, target=AlignmentFormat.NEXUS)
    assert out.exists()
    content = out.read_text()
    assert "#NEXUS" in content or "NEXUS" in content


def test_phylip_to_fasta(tmp_path):
    fc = FormatConverter()
    out = tmp_path / "out.fa"
    fc.convert(FIXTURES / "test.phy", out, target=AlignmentFormat.FASTA)
    assert out.exists()
    content = out.read_text()
    assert ">Taxon_A" in content


def test_unsupported_format_raises():
    fc = FormatConverter()
    with pytest.raises(ValueError, match="Cannot detect"):
        fc.detect(Path("alignment.xyz"))
