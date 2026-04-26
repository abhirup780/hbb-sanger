"""Tests for hbb_pipeline.alignment — reverse-complement and consensus building."""

from pathlib import Path

import pytest

from hbb_pipeline.alignment import build_consensus, reverse_complement_trace
from hbb_pipeline.models import AlignedRead, TraceData


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_trace(seq: str, channels: dict[str, list[int]] | None = None) -> TraceData:
    n = len(seq)
    if channels is None:
        channels = {b: [i * 10 for i in range(n)] for b in "ACGT"}
    return TraceData(
        sequence=seq,
        phred_scores=[30] * n,
        channels=channels,
        peak_positions=list(range(n)),
        sample_name="test",
    )


# ---------------------------------------------------------------------------
# Reverse-complement tests
# ---------------------------------------------------------------------------

def test_rc_sequence_complement_and_reverse() -> None:
    trace = _make_trace("GAATTC")  # EcoRI palindrome
    rc = reverse_complement_trace(trace)
    assert rc.sequence == "GAATTC"  # palindrome → same


def test_rc_non_palindrome() -> None:
    trace = _make_trace("ATCG")
    rc = reverse_complement_trace(trace)
    assert rc.sequence == "CGAT"


def test_rc_double_application_is_identity() -> None:
    seq = "ATCGATCGATCG"
    trace = _make_trace(seq)
    rc_rc = reverse_complement_trace(reverse_complement_trace(trace))
    assert rc_rc.sequence == seq
    assert rc_rc.phred_scores == trace.phred_scores
    assert rc_rc.peak_positions == trace.peak_positions


def test_rc_phred_reversed() -> None:
    trace = _make_trace("ATCG")
    trace2 = TraceData(
        sequence="ATCG",
        phred_scores=[10, 20, 30, 40],
        channels={"A": [0]*4, "C": [0]*4, "G": [0]*4, "T": [0]*4},
        peak_positions=[0, 1, 2, 3],
        sample_name="q",
    )
    rc = reverse_complement_trace(trace2)
    assert rc.phred_scores == [40, 30, 20, 10]


def test_rc_channel_swap() -> None:
    """After RC: channels are complement-swapped AND array-reversed.

    peak_positions are remapped as new_p = n_scans-1 - old_p (so they stay
    ascending for plotting).  For a detect_secondary_peaks lookup to return the
    correct intensity after this remapping, the channel arrays must also be
    reversed: channel[new_p] == reversed(complement_channel)[new_p]
                               == complement_channel[n_scans-1 - new_p]
                               == complement_channel[old_p]  ✓

    In this symmetric test n_scans==n==6 and peaks==[0..5], so the remapped
    peaks happen to equal [0..5] again (each new_p = 5-old_p, applied to the
    reversed list, cancels out).  The channel reversal is still required and
    visible: rc.channels["A"] = reversed(orig.T).
    """
    n = 6
    a_channel = [i * 10 for i in range(n)]  # [0,10,20,30,40,50]
    t_channel = [i * 5 for i in range(n)]   # [0, 5,10,15,20,25]
    channels = {
        "A": a_channel,
        "C": [0] * n,
        "G": [0] * n,
        "T": t_channel,
    }
    trace = TraceData(
        sequence="ATCGAT",
        phred_scores=[30] * n,
        channels=channels,
        peak_positions=list(range(n)),
        sample_name="chan",
    )
    rc = reverse_complement_trace(trace)

    # Channels are complement-swapped AND reversed so scan lookups stay correct
    assert rc.channels["A"] == list(reversed(t_channel)), "rc.A should equal reversed(orig.T)"
    assert rc.channels["T"] == list(reversed(a_channel)), "rc.T should equal reversed(orig.A)"
    # peak_positions remapped: new_p = n_scans-1 - old_p applied to reversed list.
    # For the symmetric [0..5] / n_scans=6 case this yields [0..5] again.
    n_scans = n
    expected_peaks = [n_scans - 1 - p for p in reversed(range(n))]
    assert rc.peak_positions == expected_peaks, f"peaks should be remapped, got {rc.peak_positions}"


def test_rc_sets_is_reverse_complemented() -> None:
    trace = _make_trace("ATCG")
    assert trace.is_reverse_complemented is False
    rc = reverse_complement_trace(trace)
    assert rc.is_reverse_complemented is True


# ---------------------------------------------------------------------------
# build_consensus tests
# ---------------------------------------------------------------------------

def _make_aligned(
    seq: str,
    ref_seq: str,
    ref_start: int,
    quals: list[int] | None = None,
) -> AlignedRead:
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


def test_consensus_fwd_only_region() -> None:
    """Positions covered only by fwd use fwd base."""
    from hbb_pipeline.reference import HBBReference
    ref = HBBReference(Path(__file__).parent.parent / "reference" / "HBB_reference.fasta")

    # Fwd covers [978, 988), rev covers [985, 995)
    fwd = _make_aligned("ATGGTGCATC", ref.upper_seq[978:988], 978)
    rev = _make_aligned("ATCTGACTCC", ref.upper_seq[985:995], 985)

    consensus, quals = build_consensus(fwd, rev, ref)

    # Position 978 is only in fwd
    assert consensus[978] == ref.upper_seq[978], f"Expected {ref.upper_seq[978]} at 978"
    assert quals[978] == 30


def test_consensus_agreement_boosts_quality() -> None:
    from hbb_pipeline.reference import HBBReference
    ref = HBBReference(Path(__file__).parent.parent / "reference" / "HBB_reference.fasta")

    seq = ref.upper_seq[978:998]
    fwd = _make_aligned(seq, seq, 978, quals=[30] * 20)
    rev = _make_aligned(seq, seq, 978, quals=[30] * 20)

    _, quals = build_consensus(fwd, rev, ref)
    # Agreement adds +3 up to Q60
    assert quals[978] == min(60, 30 + 3)


def test_consensus_disagreement_high_quality_diff() -> None:
    """When q_diff >= 10, use the higher-quality base."""
    from hbb_pipeline.reference import HBBReference
    ref = HBBReference(Path(__file__).parent.parent / "reference" / "HBB_reference.fasta")

    seq = ref.upper_seq[978:979]  # one base
    fwd = _make_aligned("A", seq, 978, quals=[35])
    rev = _make_aligned("T", seq, 978, quals=[20])  # diff = 15 >= 10

    consensus, quals = build_consensus(fwd, rev, ref)
    assert consensus[978] == "A"  # higher quality wins
    assert quals[978] == 15


def test_consensus_disagreement_low_quality_diff_emits_N() -> None:
    """When |q_fwd - q_rev| < 10, emit N with quality 0."""
    from hbb_pipeline.reference import HBBReference
    ref = HBBReference(Path(__file__).parent.parent / "reference" / "HBB_reference.fasta")

    seq = ref.upper_seq[978:979]
    fwd = _make_aligned("A", seq, 978, quals=[25])
    rev = _make_aligned("T", seq, 978, quals=[20])  # diff = 5 < 10

    consensus, quals = build_consensus(fwd, rev, ref)
    assert consensus[978] == "N"
    assert quals[978] == 0
