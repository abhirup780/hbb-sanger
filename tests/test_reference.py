"""Tests for hbb_pipeline.reference — safety gate for the HBB reference FASTA."""

from pathlib import Path

import pytest

from hbb_pipeline.reference import HBBReference, ReferenceValidationError

REF_PATH = Path(__file__).parent.parent / "reference" / "HBB_reference.fasta"


@pytest.fixture(scope="module")
def ref() -> HBBReference:
    return HBBReference(REF_PATH)


def test_reference_loads_with_expected_structure(ref: HBBReference) -> None:
    assert ref.length == 2750

    # Segment lengths
    assert ref.utr5 == (0, 978)
    assert ref.exons[0] == (978, 1070)
    assert ref.introns[0] == (1070, 1200)
    assert ref.exons[1] == (1200, 1423)
    assert ref.introns[1] == (1423, 2023)
    assert ref.exons[2] == (2023, 2152)
    assert ref.utr3 == (2152, 2750)

    # Lengths of each segment
    assert ref.exons[0][1] - ref.exons[0][0] == 92
    assert ref.introns[0][1] - ref.introns[0][0] == 130
    assert ref.exons[1][1] - ref.exons[1][0] == 223
    assert ref.introns[1][1] - ref.introns[1][0] == 600
    assert ref.exons[2][1] - ref.exons[2][0] == 129
    assert ref.utr3[1] - ref.utr3[0] == 598

    # CDS length
    assert len(ref.cds) == 444

    # ATG start
    assert ref.upper_seq[978:981] == "ATG"

    # TAA stop (last 3 bases of exon 3)
    assert ref.upper_seq[2149:2152] == "TAA"


def test_cds_translates_to_canonical_beta_globin(ref: HBBReference) -> None:
    from Bio.Seq import Seq

    protein = str(Seq(ref.cds).translate())
    assert protein.endswith("*"), "CDS should end with stop codon"
    translated = protein.rstrip("*")
    assert len(translated) == 147
    assert translated.startswith("MVH"), f"β-globin starts MVH, got {translated[:3]}"
    assert translated.endswith("KYH"), f"β-globin ends KYH, got {translated[-3:]}"


def test_region_of_boundaries(ref: HBBReference) -> None:
    assert ref.region_of(0) == "5UTR"
    assert ref.region_of(977) == "5UTR"
    assert ref.region_of(978) == "exon1"
    assert ref.region_of(1069) == "exon1"
    assert ref.region_of(1070) == "intron1"
    assert ref.region_of(1199) == "intron1"
    assert ref.region_of(1200) == "exon2"
    assert ref.region_of(1422) == "exon2"
    assert ref.region_of(1423) == "intron2"
    assert ref.region_of(2022) == "intron2"
    assert ref.region_of(2023) == "exon3"
    assert ref.region_of(2151) == "exon3"
    assert ref.region_of(2152) == "3UTR"
    assert ref.region_of(2749) == "3UTR"


def test_reference_rejects_wrong_length(tmp_path: Path) -> None:
    bad_fasta = tmp_path / "bad.fasta"
    bad_fasta.write_text(">bad\nACGT\n")
    with pytest.raises(ReferenceValidationError, match="length"):
        HBBReference(bad_fasta)
