"""ABI trace file parsing and quality trimming.

Parses .ab1 / .abi files via Biopython SeqIO and extracts:
- sequence (PBAS1)
- Phred quality scores (PCON1)
- four-channel intensity arrays (DATA9-12, ordered via FWO_)
- peak positions (PLOC1)

Mott's modified trimming algorithm is implemented for read-end quality trimming.
"""

from __future__ import annotations

import logging
from pathlib import Path

from Bio import SeqIO

from hbb_pipeline.models import TraceData

logger = logging.getLogger(__name__)


class InvalidTraceFileError(Exception):
    """Raised when an ABI file cannot be parsed or is structurally invalid."""


# ABI internal channel order tag — value like b"GATC" means DATA9=G, DATA10=A, ...
_FWO_TAG = "FWO_"
_DATA_TAGS = ["DATA9", "DATA10", "DATA11", "DATA12"]
_FALLBACK_ORDER = "GATC"  # Sanger convention if FWO_ absent


def parse_abi(path: Path) -> TraceData:
    """Parse an ABI trace file into a TraceData model.

    Reads channel order from the FWO_ tag so no channel is hardcoded.

    Args:
        path: Path to the .ab1 or .abi trace file.

    Returns:
        TraceData with sequence, phred_scores, channels, peak_positions,
        and sample_name populated.

    Raises:
        InvalidTraceFileError: if the file cannot be read or required tags are missing.

    Example:
        trace = parse_abi(Path("sample_fwd.ab1"))
        print(trace.sequence[:10])
    """
    if not path.exists():
        raise InvalidTraceFileError(f"Trace file not found: {path}")

    try:
        record = SeqIO.read(str(path), "abi")
    except Exception as exc:
        raise InvalidTraceFileError(f"Failed to parse ABI file {path}: {exc}") from exc

    annotations = record.annotations.get("abif_raw", {})

    # --- Sequence ---
    sequence = str(record.seq)

    # --- Phred scores ---
    raw_qual = annotations.get("PCON1")
    if raw_qual is None:
        raise InvalidTraceFileError(f"PCON1 (quality scores) missing in {path}")
    if isinstance(raw_qual, (bytes, bytearray)):
        phred_scores = list(raw_qual)
    else:
        phred_scores = list(raw_qual)

    if len(phred_scores) != len(sequence):
        logger.warning(
            "%s: PCON1 length %d != sequence length %d; truncating to shorter",
            path.name,
            len(phred_scores),
            len(sequence),
        )
        min_len = min(len(phred_scores), len(sequence))
        phred_scores = phred_scores[:min_len]
        sequence = sequence[:min_len]

    # --- Channel order from FWO_ ---
    fwo = annotations.get(_FWO_TAG)
    if fwo is None:
        logger.warning("%s: FWO_ tag absent; assuming %s channel order", path.name, _FALLBACK_ORDER)
        channel_order = _FALLBACK_ORDER
    else:
        channel_order = fwo.decode("ascii") if isinstance(fwo, (bytes, bytearray)) else str(fwo)

    # Map base letter → raw intensity array
    channels: dict[str, list[int]] = {}
    for base, tag in zip(channel_order, _DATA_TAGS):
        raw = annotations.get(tag)
        if raw is None:
            raise InvalidTraceFileError(f"{tag} (channel {base}) missing in {path}")
        channels[base] = list(raw)

    # Ensure all four bases are present
    for base in "ACGT":
        if base not in channels:
            raise InvalidTraceFileError(
                f"Channel for base '{base}' not found in {path}. "
                f"FWO_ gave order: {channel_order!r}"
            )

    # --- Peak positions ---
    ploc = annotations.get("PLOC1")
    if ploc is None:
        raise InvalidTraceFileError(f"PLOC1 (peak positions) missing in {path}")
    peak_positions = list(ploc)

    if len(peak_positions) != len(sequence):
        logger.warning(
            "%s: PLOC1 length %d != sequence length %d; truncating",
            path.name,
            len(peak_positions),
            len(sequence),
        )
        min_len = min(len(peak_positions), len(sequence))
        peak_positions = peak_positions[:min_len]
        sequence = sequence[:min_len]
        phred_scores = phred_scores[:min_len]

    sample_name = record.name or path.stem

    logger.info(
        "Parsed %s: %d bases, mean Q=%.1f",
        sample_name,
        len(sequence),
        sum(phred_scores) / max(len(phred_scores), 1),
    )

    return TraceData(
        sequence=sequence,
        phred_scores=phred_scores,
        channels=channels,
        peak_positions=peak_positions,
        sample_name=sample_name,
    )


def trim_by_quality(
    trace: TraceData,
    window: int = 20,
    threshold: float = 0.01,
) -> tuple[TraceData, tuple[int, int]]:
    """Trim trace ends using Mott's modified quality trimming algorithm.

    The algorithm finds the highest-scoring contiguous subsequence where
    each base contributes (threshold - error_probability). Low-quality ends
    contribute negative scores and are trimmed off.

    The *full* trace (channels, peak_positions) is kept in the returned
    TraceData because the chromatogram viewer needs it; only sequence,
    phred_scores, and peak_positions are sliced to the trimmed window.

    Args:
        trace: Input TraceData (full, untrimmed).
        window: Minimum window size for initial smoothing (not used in the
                pure Mott implementation but kept for API compatibility).
        threshold: Error probability threshold (default 0.01 ≈ Q20).

    Returns:
        (trimmed_trace, (trim_start, trim_end)) where trim_start/end are
        0-based indices into the original trace sequence.

    Example:
        trimmed, (s, e) = trim_by_quality(trace)
        print(f"Kept bases {s}..{e} ({e-s} bp)")
    """
    quals = trace.phred_scores
    n = len(quals)

    if n == 0:
        return trace, (0, 0)

    # Mott's algorithm: score_i = threshold - 10^(-Q/10)
    # High-quality bases (small error prob) give positive scores; low-quality give negative.
    scores = [threshold - 10 ** (-q / 10.0) for q in quals]

    # Kadane-style maximum subarray to find trim window
    max_sum = 0.0
    current_sum = 0.0
    best_start = 0
    best_end = 0
    temp_start = 0

    for i, s in enumerate(scores):
        current_sum += s
        if current_sum < 0:
            current_sum = 0.0
            temp_start = i + 1
        if current_sum > max_sum:
            max_sum = current_sum
            best_start = temp_start
            best_end = i + 1  # half-open

    if best_end == 0:
        # No base above threshold — return an empty (zero-length) trim window
        logger.warning("%s: Mott trimming removed all bases", trace.sample_name)
        best_start = 0
        best_end = 0

    logger.debug(
        "%s: Mott trim [%d, %d) of %d total bases",
        trace.sample_name,
        best_start,
        best_end,
        n,
    )

    # If trim window is empty (all bases below threshold), keep at least 1 base
    # so downstream callers receive a valid (though tiny) TraceData.
    if best_end <= best_start:
        best_start = 0
        best_end = min(1, n)

    trimmed = TraceData(
        sequence=trace.sequence[best_start:best_end],
        phred_scores=trace.phred_scores[best_start:best_end],
        channels=trace.channels,          # full channels preserved for chromatogram
        peak_positions=trace.peak_positions[best_start:best_end],
        sample_name=trace.sample_name,
        is_reverse_complemented=trace.is_reverse_complemented,
    )

    return trimmed, (best_start, best_end)
