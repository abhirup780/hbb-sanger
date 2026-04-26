"""Tests for hbb_pipeline.heterozygosity."""

import pytest

from hbb_pipeline.heterozygosity import classify_zygosity, detect_secondary_peaks
from hbb_pipeline.models import TraceData, Zygosity


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_trace_with_channels(channels: dict[str, list[int]], peak_positions: list[int]) -> TraceData:
    n = len(peak_positions)
    return TraceData(
        sequence="A" * n,
        phred_scores=[30] * n,
        channels=channels,
        peak_positions=peak_positions,
        sample_name="test",
    )


# ---------------------------------------------------------------------------
# detect_secondary_peaks
# ---------------------------------------------------------------------------

def test_het_detection_clean_secondary_peak() -> None:
    """A=1000, T=500 → T is secondary at ratio 0.5."""
    scan_len = 10
    channels = {
        "A": [0] * 3 + [1000] + [0] * (scan_len - 4),
        "C": [0] * scan_len,
        "G": [0] * scan_len,
        "T": [0] * 3 + [500] + [0] * (scan_len - 4),
    }
    trace = _make_trace_with_channels(channels, peak_positions=[3])
    peaks = detect_secondary_peaks(trace, trace_index=0)

    assert "A" in peaks
    assert peaks["A"] == pytest.approx(1.0)
    assert "T" in peaks
    assert peaks["T"] == pytest.approx(0.5)
    assert "C" not in peaks
    assert "G" not in peaks


def test_het_not_called_below_threshold() -> None:
    """A=1000, T=200 (0.2 ratio, below 0.25 threshold) → T not included."""
    scan_len = 10
    channels = {
        "A": [0] * 3 + [1000] + [0] * (scan_len - 4),
        "C": [0] * scan_len,
        "G": [0] * scan_len,
        "T": [0] * 3 + [200] + [0] * (scan_len - 4),
    }
    trace = _make_trace_with_channels(channels, peak_positions=[3])
    peaks = detect_secondary_peaks(trace, trace_index=0)

    assert "T" not in peaks
    assert "A" in peaks


def test_detect_secondary_out_of_bounds_returns_empty() -> None:
    channels = {"A": [0] * 5, "C": [0] * 5, "G": [0] * 5, "T": [0] * 5}
    trace = _make_trace_with_channels(channels, peak_positions=[2])
    assert detect_secondary_peaks(trace, trace_index=-1) == {}
    assert detect_secondary_peaks(trace, trace_index=99) == {}


def test_detect_secondary_zero_intensity_returns_empty() -> None:
    channels = {"A": [0] * 5, "C": [0] * 5, "G": [0] * 5, "T": [0] * 5}
    trace = _make_trace_with_channels(channels, peak_positions=[2])
    assert detect_secondary_peaks(trace, trace_index=0) == {}


def test_custom_threshold() -> None:
    """With threshold=0.4, a 0.5-ratio peak is still detected; 0.3-ratio is not."""
    scan_len = 10
    channels = {
        "A": [0] * 2 + [1000] + [0] * (scan_len - 3),
        "C": [0] * 2 + [500] + [0] * (scan_len - 3),   # 0.5 ratio
        "G": [0] * 2 + [300] + [0] * (scan_len - 3),   # 0.3 ratio
        "T": [0] * scan_len,
    }
    trace = _make_trace_with_channels(channels, peak_positions=[2])
    peaks = detect_secondary_peaks(trace, trace_index=0, ratio_threshold=0.4)

    assert "C" in peaks      # 0.5 ≥ 0.4
    assert "G" not in peaks  # 0.3 < 0.4


# ---------------------------------------------------------------------------
# classify_zygosity
# ---------------------------------------------------------------------------

def test_classify_het_concordant_both_strands() -> None:
    fwd = {"A": 1.0, "T": 0.42}
    rev = {"A": 1.0, "T": 0.38}
    zyg, fr, rr = classify_zygosity(fwd, rev, "T")
    assert zyg == Zygosity.HET
    assert fr == pytest.approx(0.42)
    assert rr == pytest.approx(0.38)


def test_classify_hom_no_secondary_peaks() -> None:
    # Alt is the only base on both strands
    fwd = {"T": 1.0}
    rev = {"T": 1.0}
    zyg, fr, rr = classify_zygosity(fwd, rev, "T")
    assert zyg == Zygosity.HOM


def test_classify_unknown_missing_one_strand() -> None:
    fwd = {"A": 1.0, "T": 0.4}
    rev: dict = {}
    zyg, _, _ = classify_zygosity(fwd, rev, "T")
    assert zyg == Zygosity.UNKNOWN


def test_classify_unknown_both_missing() -> None:
    zyg, _, _ = classify_zygosity({}, {}, "T")
    assert zyg == Zygosity.UNKNOWN


def test_classify_unknown_discordant() -> None:
    # fwd shows T as secondary, rev does not
    fwd = {"A": 1.0, "T": 0.35}
    rev = {"A": 1.0}
    zyg, _, _ = classify_zygosity(fwd, rev, "T")
    assert zyg == Zygosity.UNKNOWN
