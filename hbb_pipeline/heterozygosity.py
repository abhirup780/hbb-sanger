"""Secondary peak detection and zygosity classification.

Two-level API:
  1. get_two_peaks() — pre-alignment, per-position secondary peak
     detection using a rule-based slope filter.  Used by
     apply_iupac_symbols() to IUPAC-code the trace before alignment.
  2. detect_secondary_peaks() / classify_zygosity() — post-call, used to
     assign HOM/HET/UNKNOWN to variants already called from the consensus.
"""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING

from hbb_pipeline.models import Zygosity

if TYPE_CHECKING:
    from hbb_pipeline.models import TraceData

logger = logging.getLogger(__name__)

_HET_THRESHOLD = 0.25  # post-call zygosity threshold — do not raise
_PRE_ALIGN_CUTOFF = 0.30  # pre-alignment secondary peak threshold
_SMALL_INC = 10   # smallIncrement for slope check
_BIG_INC = 20     # bigIncrement for slope check


# ---------------------------------------------------------------------------
# Pre-alignment secondary peak detection
# ---------------------------------------------------------------------------

def _check_rising_right(
    channel: list[int],
    scan_pos: int,
    right_end: int,
) -> bool:
    """True if channel rises consistently from scan_pos to right_end.

    A True return means the secondary peak is on a monotonic rising slope,
    not a genuine local peak.
    """
    i = scan_pos + 1
    while i <= right_end:
        if i >= len(channel):
            break
        increment = channel[i] - channel[i - 1]
        if i == scan_pos + 2 or i == right_end - 1:
            if increment < _SMALL_INC:
                break
        elif scan_pos + 2 < i < right_end - 1:
            if increment < _BIG_INC:
                break
        i += 1
    return i > right_end


def _check_rising_left(
    channel: list[int],
    scan_pos: int,
    left_end: int,
) -> bool:
    """True if channel rises consistently going left from scan_pos to left_end.

    A True return means the secondary peak is on the falling edge of a
    peak to the left — not a genuine local peak at this position.
    """
    i = scan_pos - 1
    while i >= left_end:
        if i < 0:
            break
        increment = channel[i] - channel[i + 1]   # positive when value rises leftward
        if i == left_end + 1 or i == scan_pos - 2:
            if increment < _SMALL_INC:
                break
        elif left_end + 1 < i < scan_pos - 2:
            if increment < _BIG_INC:
                break
        i -= 1
    return i < left_end


def get_two_peaks(
    trace: "TraceData",
    trace_index: int,
    cutoff: float = _PRE_ALIGN_CUTOFF,
) -> tuple[str, str | None]:
    """Pre-alignment secondary peak detection with rule-based slope filtering.

    Reads raw channel intensities at trace.peak_positions[trace_index],
    applies the ratio threshold, then filters out false secondary peaks
    that lie on monotonic slopes (not genuine local maxima).

    Returns:
        (primary_base, secondary_base)
        secondary_base is None when no valid secondary peak is found.
        primary_base is always the channel-intensity winner (may differ
        from the basecaller's call at noisy positions).

    Falls back to (original_base, None) when channel data is missing.
    """
    bases = ("A", "C", "G", "T")
    fallback_base = trace.sequence[trace_index] if trace_index < len(trace.sequence) else "N"

    if trace_index < 0 or trace_index >= len(trace.peak_positions):
        return fallback_base, None

    scan_pos = trace.peak_positions[trace_index]

    # Read all four channel heights at the basecall scan position
    heights: dict[str, int] = {}
    for b in bases:
        ch = trace.channels.get(b, [])
        heights[b] = ch[scan_pos] if scan_pos < len(ch) else 0

    primary_h = max(heights.values())
    if primary_h == 0:
        return fallback_base, None

    sorted_by_h = sorted(bases, key=lambda b: heights[b], reverse=True)
    primary = sorted_by_h[0]
    secondary = sorted_by_h[1]
    secondary_h = heights[secondary]

    # Ratio threshold (default 0.30)
    if secondary_h / primary_h < cutoff:
        return primary, None

    # --- Rule-based slope filter: direction + continuation check ---
    sec_ch = trace.channels.get(secondary, [])
    trace_len = max((len(trace.channels.get(b, [])) for b in bases), default=0)

    if trace_len == 0 or scan_pos >= trace_len:
        return primary, secondary   # can't validate — accept

    # Direction: uses immediate neighbours (position ± 1 scan points)
    if scan_pos == 0:
        direction = 1   # edge → treat as rising right
    elif scan_pos == trace_len - 1:
        direction = -1  # edge → treat as falling from left
    else:
        left_val = sec_ch[scan_pos - 1] if scan_pos - 1 < len(sec_ch) else 0
        right_val = sec_ch[scan_pos + 1] if scan_pos + 1 < len(sec_ch) else 0
        center_val = secondary_h

        if left_val < center_val < right_val:
            direction = 1   # still rising rightward
        elif left_val > center_val > right_val:
            direction = -1  # falling (peak was to the left)
        elif left_val > center_val and center_val < right_val:
            direction = 2   # valley between two peaks
        else:
            direction = 0   # peak shape (L < C > R) — genuine local max

    if direction == 0:
        return primary, secondary   # unambiguous peak, no filtering needed

    peak_positions = trace.peak_positions
    n_seq = len(trace.sequence)

    right_peak_found = False
    left_peak_found = False

    if direction in (1, 2) and scan_pos + 7 < trace_len:
        right_end = (
            trace_len - 1
            if trace_index >= n_seq - 1
            else peak_positions[trace_index + 1]
            if trace_index + 1 < len(peak_positions)
            else trace_len - 1
        )
        right_peak_found = _check_rising_right(sec_ch, scan_pos, right_end)

    if direction in (-1, 2) and scan_pos - 7 >= 0:
        left_end = (
            0
            if trace_index <= 0
            else peak_positions[trace_index - 1]
            if trace_index - 1 >= 0
            else 0
        )
        left_peak_found = _check_rising_left(sec_ch, scan_pos, left_end)

    # Apply filter: monotonic slope → not a genuine peak
    if direction == 1 and right_peak_found:
        return primary, None
    if direction == -1 and left_peak_found:
        return primary, None
    if direction == 2 and right_peak_found and left_peak_found:
        return primary, None

    return primary, secondary


def detect_secondary_peaks(
    trace: "TraceData",
    trace_index: int,
    ratio_threshold: float = _HET_THRESHOLD,
) -> dict[str, float]:
    """Detect secondary peaks at a given base position in the trace.

    Reads all four channel intensities at the scan position corresponding to
    trace_index (i.e., trace.peak_positions[trace_index]).  The primary channel
    (highest intensity) is normalised to 1.0; any other channel ≥ ratio_threshold
    is included in the returned dict.

    Args:
        trace: Parsed TraceData.
        trace_index: Index into trace.sequence (and trace.peak_positions).
        ratio_threshold: Minimum intensity ratio relative to primary to count as
                         a secondary peak.  Default 0.25 — do not raise this.

    Returns:
        Dict mapping base letter to intensity ratio, e.g. {"A": 1.0, "T": 0.42}.
        Returns {} if trace_index is out of bounds or primary intensity is 0.

    Example:
        peaks = detect_secondary_peaks(trace, 45, ratio_threshold=0.25)
        if "T" in peaks:
            print(f"Secondary T at ratio {peaks['T']:.2f}")
    """
    if trace_index < 0 or trace_index >= len(trace.peak_positions):
        return {}

    scan_pos = trace.peak_positions[trace_index]

    intensities: dict[str, int] = {}
    for base in "ACGT":
        channel = trace.channels.get(base, [])
        if scan_pos < len(channel):
            intensities[base] = channel[scan_pos]
        else:
            intensities[base] = 0

    primary_intensity = max(intensities.values())
    if primary_intensity == 0:
        return {}

    result: dict[str, float] = {}
    for base, intensity in intensities.items():
        ratio = intensity / primary_intensity
        if ratio >= ratio_threshold:
            result[base] = round(ratio, 4)

    return result


def classify_zygosity(
    fwd_peaks: dict[str, float],
    rev_peaks: dict[str, float],
    variant_alt: str,
) -> tuple[Zygosity, float | None, float | None]:
    """Classify a variant as HOM, HET, or UNKNOWN from secondary peak evidence.

    Logic:
    - HET: both fwd and rev show the alt base as a secondary peak AND the
           reference base is also present (either strand) — strand concordance
           is the strongest evidence.
    - HOM: neither trace shows a secondary peak; the alt base IS the primary
           in both (i.e., only the alt base appears in both peak dicts at ratio 1.0).
    - UNKNOWN: discordant between strands, only one strand covers, or ambiguous.

    Args:
        fwd_peaks: Output of detect_secondary_peaks for the forward trace.
        rev_peaks: Output of detect_secondary_peaks for the reverse trace.
        variant_alt: The alt allele base (single character, e.g. "T").

    Returns:
        (zygosity, fwd_alt_ratio, rev_alt_ratio)

    Example:
        zyg, fr, rr = classify_zygosity({"A": 1.0, "T": 0.42}, {"A": 1.0, "T": 0.38}, "T")
        # → (Zygosity.HET, 0.42, 0.38)
    """
    alt = variant_alt.upper()
    fwd_ratio = fwd_peaks.get(alt)
    rev_ratio = rev_peaks.get(alt)

    # One or both strands missing
    if not fwd_peaks and not rev_peaks:
        return Zygosity.UNKNOWN, None, None

    if not fwd_peaks or not rev_peaks:
        # Only one strand covers this position.
        # alt is the primary call → HOM (reasonably reliable from one strand).
        # alt is a secondary peak → UNKNOWN: HET requires dual-strand concordance.
        single = fwd_peaks or rev_peaks
        primary = _primary_base(single)
        if primary == alt:
            return Zygosity.HOM, fwd_ratio, rev_ratio
        return Zygosity.UNKNOWN, fwd_ratio, rev_ratio

    fwd_alt_present = alt in fwd_peaks
    rev_alt_present = alt in rev_peaks

    # Determine primary base on each strand
    fwd_primary = _primary_base(fwd_peaks)
    rev_primary = _primary_base(rev_peaks)

    # HOM: alt is the dominant call on both strands with no other prominent peaks
    if (
        fwd_primary == alt
        and rev_primary == alt
        and len(fwd_peaks) == 1
        and len(rev_peaks) == 1
    ):
        return Zygosity.HOM, fwd_ratio, rev_ratio

    # HET: alt is a secondary peak on both strands
    if fwd_alt_present and rev_alt_present:
        # Also require that ref base (non-alt) appears somewhere
        fwd_has_ref = any(b != alt for b in fwd_peaks)
        rev_has_ref = any(b != alt for b in rev_peaks)
        if fwd_has_ref or rev_has_ref:
            return Zygosity.HET, fwd_ratio, rev_ratio

    # HOM with secondary noise: alt is primary on both, some noise channel present
    if fwd_primary == alt and rev_primary == alt:
        return Zygosity.HOM, fwd_ratio, rev_ratio

    # Discordant or ambiguous
    logger.debug(
        "Ambiguous zygosity for alt=%s: fwd_peaks=%s, rev_peaks=%s",
        alt,
        fwd_peaks,
        rev_peaks,
    )
    return Zygosity.UNKNOWN, fwd_ratio, rev_ratio


def _primary_base(peaks: dict[str, float]) -> str:
    """Return the base with the highest ratio (primary call)."""
    return max(peaks, key=lambda b: peaks[b])
