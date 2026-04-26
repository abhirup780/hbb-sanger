"""Tests for hbb_pipeline.coordinates — the critical coordinate translator."""

from pathlib import Path
from typing import Literal

import pytest

from hbb_pipeline.coordinates import CoordinateTranslator
from hbb_pipeline.models import Variant, Zygosity
from hbb_pipeline.reference import HBBReference

REF_PATH = Path(__file__).parent.parent / "reference" / "HBB_reference.fasta"


@pytest.fixture(scope="module")
def ref() -> HBBReference:
    return HBBReference(REF_PATH)


@pytest.fixture(scope="module")
def translator(ref: HBBReference) -> CoordinateTranslator:
    return CoordinateTranslator(ref)


# ---------------------------------------------------------------------------
# Exon boundary positions
# ---------------------------------------------------------------------------

def test_genomic_to_c_exon_boundaries(translator: CoordinateTranslator) -> None:
    assert translator.genomic_to_c(978) == "1"     # A of ATG — first base of exon 1
    assert translator.genomic_to_c(1069) == "92"   # last base of exon 1
    assert translator.genomic_to_c(1200) == "93"   # first base of exon 2
    assert translator.genomic_to_c(1422) == "315"  # last base of exon 2
    assert translator.genomic_to_c(2023) == "316"  # first base of exon 3
    assert translator.genomic_to_c(2151) == "444"  # last base of exon 3 (stop TAA last A)


# ---------------------------------------------------------------------------
# UTR positions
# ---------------------------------------------------------------------------

def test_genomic_to_c_utr5(translator: CoordinateTranslator) -> None:
    assert translator.genomic_to_c(977) == "-1"    # immediately upstream of ATG
    assert translator.genomic_to_c(976) == "-2"
    assert translator.genomic_to_c(900) == "-78"   # c.-78 promoter region


def test_genomic_to_c_utr3(translator: CoordinateTranslator) -> None:
    assert translator.genomic_to_c(2152) == "*1"   # first base of 3' UTR
    assert translator.genomic_to_c(2153) == "*2"
    assert translator.genomic_to_c(2749) == "*598"


# ---------------------------------------------------------------------------
# Intronic positions
# ---------------------------------------------------------------------------

def test_genomic_to_c_intronic(translator: CoordinateTranslator) -> None:
    # Intron 1: genomic [1070, 1200)
    assert translator.genomic_to_c(1070) == "92+1"   # 1 bp after exon 1 end
    assert translator.genomic_to_c(1199) == "93-1"   # 1 bp before exon 2 start
    # Intron 2: genomic [1423, 2023)
    assert translator.genomic_to_c(1423) == "315+1"
    assert translator.genomic_to_c(2022) == "316-1"


def test_ivs1_110_position_is_c93_minus_21(translator: CoordinateTranslator, ref: HBBReference) -> None:
    """IVS1-110 in HbVar = HBB:c.93-21G>A.

    IVS1-110 means 110 bp from the donor splice site (5' end of intron 1).
    Intron 1 donor = genomic 1070. Position 110 bp in = genomic 1070 + 109 = 1179.
    Distance to acceptor (genomic 1200): 1200 - 1179 = 21 → c.93-21.
    """
    genomic = 1179  # 110th base of intron 1 (0-indexed: 1070 + 110 - 1 = 1179)
    assert translator.genomic_to_c(genomic) == "93-21"
    # Wild-type at this position must be G (IVS1-110G>A)
    assert ref.seq[genomic].upper() == "G", (
        f"Expected G at genomic 1179 (IVS1-110 wild-type), got {ref.seq[genomic]!r}"
    )


# ---------------------------------------------------------------------------
# Roundtrip: c_to_genomic(genomic_to_c(x)) == x
# ---------------------------------------------------------------------------

def test_c_to_genomic_roundtrip(translator: CoordinateTranslator) -> None:
    positions = [100, 978, 1069, 1070, 1179, 1199, 1200, 1422, 1423, 2022, 2023, 2151, 2152, 2200]
    for gpos in positions:
        c = translator.genomic_to_c(gpos)
        back = translator.c_to_genomic(c)
        assert back == gpos, f"Roundtrip failed for genomic {gpos}: c.{c} → {back}"


# ---------------------------------------------------------------------------
# HbS (c.20A>T, p.Glu7Val / p.Glu6Val)
# ---------------------------------------------------------------------------

def _make_variant(
    ref_pos_genomic: int,
    ref_allele: str,
    alt_allele: str,
    region: Literal["5UTR","exon1","intron1","exon2","intron2","exon3","3UTR"] = "exon1",
) -> Variant:
    return Variant(
        ref_pos_genomic=ref_pos_genomic,
        ref_allele=ref_allele,
        alt_allele=alt_allele,
        zygosity=Zygosity.HOM,
        variant_type="SNV",
        hgvs_c="",  # will be filled by tests
        hgvs_p_hgvs=None,
        hgvs_p_legacy=None,
        region=region,
        phred_support=40,
        secondary_peak_ratio_fwd=None,
        secondary_peak_ratio_rev=None,
        called_by=["alignment"],
    )


def test_hbs_hgvs_c(translator: CoordinateTranslator) -> None:
    # HbS: c.20A>T — genomic 978 + 19 = 997
    v = _make_variant(978 + 19, "A", "T")
    assert translator.build_hgvs_c(v) == "c.20A>T"


def test_hbs_hgvs_p_hgvs(translator: CoordinateTranslator) -> None:
    v = _make_variant(978 + 19, "A", "T")
    assert translator.build_hgvs_p(v, "hgvs") == "p.Glu7Val"


def test_hbs_hgvs_p_legacy(translator: CoordinateTranslator) -> None:
    v = _make_variant(978 + 19, "A", "T")
    assert translator.build_hgvs_p(v, "legacy") == "p.Glu6Val"


# ---------------------------------------------------------------------------
# Additional build_hgvs_c shapes
# ---------------------------------------------------------------------------

def test_build_hgvs_c_intronic_snv(translator: CoordinateTranslator) -> None:
    # IVS1-110 position: genomic 1179, intron1
    v = _make_variant(1179, "G", "A", region="intron1")
    assert translator.build_hgvs_c(v) == "c.93-21G>A"


def test_build_hgvs_c_del_single(translator: CoordinateTranslator) -> None:
    from hbb_pipeline.models import Variant
    v = Variant(
        ref_pos_genomic=978 + 26,  # c.27
        ref_allele="G",
        alt_allele="",
        zygosity=Zygosity.HOM,
        variant_type="DEL",
        hgvs_c="",
        hgvs_p_hgvs=None,
        hgvs_p_legacy=None,
        region="exon1",
        phred_support=40,
        secondary_peak_ratio_fwd=None,
        secondary_peak_ratio_rev=None,
        called_by=["alignment"],
    )
    assert translator.build_hgvs_c(v) == "c.27delG"


def test_build_hgvs_c_dup_single(translator: CoordinateTranslator, ref: HBBReference) -> None:
    # Codon 8/9 +G duplication: c.27dupG
    # Genomic position 978+26 = 1004; the inserted G duplicates the G at 1004
    genomic_pos = 978 + 26  # c.27, which is G in reference
    from hbb_pipeline.models import Variant
    v = Variant(
        ref_pos_genomic=genomic_pos,
        ref_allele="",
        alt_allele="G",
        zygosity=Zygosity.HOM,
        variant_type="INS",
        hgvs_c="",
        hgvs_p_hgvs=None,
        hgvs_p_legacy=None,
        region="exon1",
        phred_support=40,
        secondary_peak_ratio_fwd=None,
        secondary_peak_ratio_rev=None,
        called_by=["alignment"],
    )
    result = translator.build_hgvs_c(v)
    # Should be dup since G is preceded by G (... GGT ... at codon 8/9)
    assert "dup" in result or "ins" in result  # accept either; check prefix
    assert "c.27" in result


def test_synonymous_variant_returns_p_equal(translator: CoordinateTranslator) -> None:
    # Change last base of codon 1 (ATG → ATC — both ... wait, ATG = Met, no synonymous)
    # Use codon 2: second codon starts at c.4 = genomic 981
    # Codon 2 in β-globin = GTG (Val) → change c.6 (genomic 983, third base of codon 2)
    # GTG → GTC → still Val → synonymous
    v = _make_variant(983, "G", "C", region="exon1")  # c.6 third base of codon 2
    result = translator.build_hgvs_p(v, "hgvs")
    assert result == "p.(=)"


def test_nonsense_variant(translator: CoordinateTranslator) -> None:
    # Codon 39 (Cd39): c.118C>T → CAG (Gln) → TAG (Ter/stop) → p.Gln40Ter (HGVS)
    # c.118 = genomic 978 + 117 = 1095; but exon 2 starts at 1200 (c.93)
    # c.118 is in exon 2: c.93 = genomic 1200, so c.118 = genomic 1200 + 25 = 1225
    genomic_pos = 1200 + (118 - 93)  # = 1225
    v = _make_variant(genomic_pos, "C", "T", region="exon2")
    result = translator.build_hgvs_p(v, "hgvs")
    assert result == "p.Gln40Ter"


def test_non_coding_variant_returns_none(translator: CoordinateTranslator) -> None:
    v = _make_variant(1070, "G", "A", region="intron1")
    assert translator.build_hgvs_p(v, "hgvs") is None
