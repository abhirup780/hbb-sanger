"""Tests for hbb_pipeline.parsing — ABI parsing and Mott trimming."""

from pathlib import Path

import pytest

from hbb_pipeline.models import TraceData
from hbb_pipeline.parsing import InvalidTraceFileError, parse_abi, trim_by_quality

FIXTURES = Path(__file__).parent / "fixtures"


# ---------------------------------------------------------------------------
# Helpers to build synthetic TraceData
# ---------------------------------------------------------------------------

def _make_trace(quals: list[int], seq: str | None = None) -> TraceData:
    """Build a minimal TraceData for trimming tests."""
    n = len(quals)
    if seq is None:
        seq = "A" * n
    dummy_channel = list(range(n * 10))  # arbitrary intensities
    return TraceData(
        sequence=seq,
        phred_scores=quals,
        channels={"A": dummy_channel, "C": dummy_channel, "G": dummy_channel, "T": dummy_channel},
        peak_positions=list(range(n)),
        sample_name="synthetic",
    )


# ---------------------------------------------------------------------------
# Mott trimming tests
# ---------------------------------------------------------------------------

def test_mott_trim_removes_low_quality_ends() -> None:
    """High-quality core flanked by low-quality ends should be trimmed."""
    # Q5 = error prob ~0.316 → score ≈ 0.01 - 0.316 = negative (bad)
    # Q30 = error prob 0.001 → score ≈ 0.01 - 0.001 = +0.009 (good)
    low = [5] * 10
    high = [30] * 50
    quals = low + high + low
    trace = _make_trace(quals)

    trimmed, (s, e) = trim_by_quality(trace, threshold=0.01)

    assert s >= 5, f"Expected trim_start >= 5, got {s}"
    assert e <= len(quals) - 5, f"Expected trim_end <= {len(quals)-5}, got {e}"
    assert all(q >= 20 for q in trimmed.phred_scores), "Trimmed region should be high quality"


def test_mott_trim_preserves_full_channels() -> None:
    """Channels array is NOT sliced — the full trace is needed by the chromatogram."""
    quals = [5] * 5 + [30] * 20 + [5] * 5
    trace = _make_trace(quals)
    trimmed, _ = trim_by_quality(trace)

    # Full channel length preserved
    assert len(trimmed.channels["A"]) == len(trace.channels["A"])


def test_mott_trim_all_high_quality_keeps_all() -> None:
    """If all bases are high quality, trimming should keep almost everything."""
    quals = [35] * 100
    trace = _make_trace(quals)
    trimmed, (s, e) = trim_by_quality(trace)

    assert s == 0
    assert e == 100


def test_mott_trim_all_low_quality_falls_back() -> None:
    """All low-quality bases should not crash — fallback returns full span."""
    quals = [3] * 30
    trace = _make_trace(quals)
    trimmed, (s, e) = trim_by_quality(trace)
    # Should not raise; trim span should be within [0, 30]
    assert 0 <= s <= e <= 30


def test_mott_trim_returns_correct_indices() -> None:
    """trim_start and trim_end must index the original (untrimmed) sequence."""
    quals = [5] * 8 + [35] * 20 + [5] * 8
    trace = _make_trace(quals)
    trimmed, (s, e) = trim_by_quality(trace)

    # Trimmed sequence must equal original sliced by (s, e)
    assert trimmed.sequence == trace.sequence[s:e]
    assert trimmed.phred_scores == trace.phred_scores[s:e]


def test_mott_trim_threshold_sensitivity() -> None:
    """A stricter threshold (lower value) should trim more aggressively."""
    quals = [20] * 5 + [30] * 40 + [20] * 5
    trace = _make_trace(quals)

    _, (s_loose, e_loose) = trim_by_quality(trace, threshold=0.01)   # Q20 is fine
    _, (s_strict, e_strict) = trim_by_quality(trace, threshold=0.001)  # Q20 is bad

    span_loose = e_loose - s_loose
    span_strict = e_strict - s_strict
    assert span_strict <= span_loose, (
        f"Stricter threshold should trim more: loose={span_loose}, strict={span_strict}"
    )


# ---------------------------------------------------------------------------
# parse_abi: integration test (skips if no fixture file present)
# ---------------------------------------------------------------------------

@pytest.mark.skipif(
    not any(FIXTURES.glob("*.ab1")) and not any(FIXTURES.glob("*.abi")),
    reason="No .ab1 fixture files in tests/fixtures/ — place real trace files there to enable",
)
def test_parse_abi_smoke() -> None:
    """Smoke test: parse the first .ab1 fixture and check field types."""
    ab1_files = list(FIXTURES.glob("*.ab1")) + list(FIXTURES.glob("*.abi"))
    trace = parse_abi(ab1_files[0])

    assert isinstance(trace.sequence, str)
    assert len(trace.sequence) > 0
    assert len(trace.phred_scores) == len(trace.sequence)
    assert len(trace.peak_positions) == len(trace.sequence)
    assert set(trace.channels.keys()) == {"A", "C", "G", "T"}
    assert all(isinstance(v, list) for v in trace.channels.values())


def test_parse_abi_missing_file_raises() -> None:
    with pytest.raises(InvalidTraceFileError, match="not found"):
        parse_abi(Path("nonexistent_file.ab1"))
