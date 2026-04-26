"""Tests for hbb_pipeline.variants — variant calling and annotation."""

from pathlib import Path
from typing import Literal

import pytest

from hbb_pipeline.coordinates import CoordinateTranslator
from hbb_pipeline.models import AlignedRead, TraceData, Variant, Zygosity
from hbb_pipeline.reference import HBBReference
from hbb_pipeline.variants import annotate_known, call_variants_from_alignment, merge_variant_calls

REF_PATH = Path(__file__).parent.parent / "reference" / "HBB_reference.fasta"


@pytest.fixture(scope="module")
def ref() -> HBBReference:
    return HBBReference(REF_PATH)


@pytest.fixture(scope="module")
def translator(ref: HBBReference) -> CoordinateTranslator:
    return CoordinateTranslator(ref)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_dummy_trace(n: int = 200) -> TraceData:
    return TraceData(
        sequence="A" * n,
        phred_scores=[30] * n,
        channels={"A": [100] * n, "C": [0] * n, "G": [0] * n, "T": [0] * n},
        peak_positions=list(range(n)),
        sample_name="dummy",
    )


def _make_aligned(seq: str, ref_seq: str, ref_start: int, quals: list[int] | None = None) -> AlignedRead:
    n = len(seq)
    if quals is None:
        quals = [30] * n
    return AlignedRead(
        read_name="test",
        reference_start=ref_start,
        reference_end=ref_start + n,
        aligned_seq=seq,
        aligned_ref=ref_seq,
        per_base_quality=quals,
        trace_index_map=list(range(n)),
    )


def _make_variant_obj(
    ref_pos: int,
    ref_allele: str,
    alt_allele: str,
    hgvs_c: str,
    region: str = "exon1",
    zygosity: Zygosity = Zygosity.HOM,
    called_by: list = None,
) -> Variant:
    return Variant(
        ref_pos_genomic=ref_pos,
        ref_allele=ref_allele,
        alt_allele=alt_allele,
        zygosity=zygosity,
        variant_type="SNV",
        hgvs_c=hgvs_c,
        hgvs_p_hgvs=None,
        hgvs_p_legacy=None,
        region=region,  # type: ignore[arg-type]
        phred_support=30,
        secondary_peak_ratio_fwd=None,
        secondary_peak_ratio_rev=None,
        called_by=called_by or ["alignment"],
    )


# ---------------------------------------------------------------------------
# HbS roundtrip test (the REQUIRED test from the spec)
# ---------------------------------------------------------------------------

def test_hbs_roundtrip_on_synthetic_consensus(ref: HBBReference, translator: CoordinateTranslator) -> None:
    """Introduce A→T at position 978+19 = 997 in the consensus and verify HbS is called."""
    # Build a synthetic consensus = reference with one SNV at codon 6 (c.20)
    consensus = list(ref.upper_seq)
    quals = [30] * ref.length

    hbs_pos = 978 + 19  # c.20 = genomic 997
    assert consensus[hbs_pos].upper() == "A", f"Expected A at {hbs_pos}, got {consensus[hbs_pos]}"
    consensus[hbs_pos] = "T"

    # Build synthetic aligned reads covering the HbS region
    window = 50
    win_start = max(0, hbs_pos - window)
    win_end = min(ref.length, hbs_pos + window)
    seg = "".join(consensus[win_start:win_end])
    ref_seg = ref.upper_seq[win_start:win_end]

    fwd_aligned = _make_aligned(seg, ref_seg, win_start)
    rev_aligned = _make_aligned(seg, ref_seg, win_start)
    fwd_trace = _make_dummy_trace(win_end - win_start)
    rev_trace = _make_dummy_trace(win_end - win_start)

    variants = call_variants_from_alignment(
        "".join(consensus),
        quals,
        ref,
        translator,
        fwd_trace,
        rev_trace,
        fwd_aligned,
        rev_aligned,
    )

    # Should call exactly the HbS SNV
    hbs_calls = [v for v in variants if v.ref_pos_genomic == hbs_pos]
    assert len(hbs_calls) == 1, f"Expected 1 HbS call, got {len(hbs_calls)}: {variants}"
    v = hbs_calls[0]
    assert v.hgvs_c == "c.20A>T", f"Expected c.20A>T, got {v.hgvs_c}"
    assert v.known_variant_name is not None, "HbS should be annotated as known"
    assert "HbS" in v.known_variant_name or "Sickle" in v.known_variant_name


# ---------------------------------------------------------------------------
# merge_variant_calls
# ---------------------------------------------------------------------------

def test_merge_both_agree_uses_tracy_zygosity() -> None:
    aln = [_make_variant_obj(997, "A", "T", "c.20A>T", zygosity=Zygosity.UNKNOWN, called_by=["alignment"])]
    trc = [_make_variant_obj(997, "A", "T", "c.20A>T", zygosity=Zygosity.HET, called_by=["tracy"])]
    merged = merge_variant_calls(aln, trc)
    assert len(merged) == 1
    assert merged[0].zygosity == Zygosity.HET
    assert set(merged[0].called_by) == {"alignment", "tracy"}
    assert merged[0].requires_manual_review is False


def test_merge_zygosity_discordance_flags_review() -> None:
    aln = [_make_variant_obj(997, "A", "T", "c.20A>T", zygosity=Zygosity.HOM, called_by=["alignment"])]
    trc = [_make_variant_obj(997, "A", "T", "c.20A>T", zygosity=Zygosity.HET, called_by=["tracy"])]
    merged = merge_variant_calls(aln, trc)
    assert merged[0].requires_manual_review is True
    assert "discordance" in (merged[0].review_reason or "")


def test_merge_alignment_only_indel_flags_review_when_tracy_ran() -> None:
    # Tracy is an indel detector; an alignment-only DEL when tracy ran should be flagged.
    aln = [_make_variant_obj(997, "AT", "A", "c.20delT")]
    aln[0] = aln[0].model_copy(update={"variant_type": "DEL"})
    merged = merge_variant_calls(aln, [], tracy_was_run=True)
    assert merged[0].requires_manual_review is True
    assert "alignment only" in (merged[0].review_reason or "")


def test_merge_alignment_only_snv_passes_through() -> None:
    # SNVs are never called by tracy, so alignment-only SNVs should not be flagged.
    aln = [_make_variant_obj(997, "A", "T", "c.20A>T")]
    merged = merge_variant_calls(aln, [], tracy_was_run=True)
    assert merged[0].requires_manual_review is False


def test_merge_tracy_only_flags_review() -> None:
    trc = [_make_variant_obj(997, "A", "T", "c.20A>T", called_by=["tracy"])]
    merged = merge_variant_calls([], trc)
    assert merged[0].requires_manual_review is True
    assert "tracy only" in (merged[0].review_reason or "")


# ---------------------------------------------------------------------------
# annotate_known
# ---------------------------------------------------------------------------

def test_annotate_known_hbs() -> None:
    v = _make_variant_obj(997, "A", "T", "c.20A>T")
    annotate_known([v])
    assert v.known_variant_name is not None
    assert "HbS" in v.known_variant_name or "Sickle" in v.known_variant_name


def test_annotate_known_unknown_variant_leaves_none() -> None:
    v = _make_variant_obj(978, "A", "C", "c.1A>C")
    annotate_known([v])
    assert v.known_variant_name is None
