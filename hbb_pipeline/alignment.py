"""Reverse-complement, pairwise alignment, and consensus building.

Alignment scoring is tuned for Sanger error profiles:
    match=2, mismatch=-3, open_gap=-15, extend_gap=-2

Gap open penalty is deliberately high (-15) to strongly prefer a mismatch
over introducing a gap.  This improves gap placement in homopolymer
regions (e.g. Cd41/42 -CTTT).
"""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING

from Bio import Align
from Bio.Seq import Seq

from hbb_pipeline.models import AlignedRead, TraceData

if TYPE_CHECKING:
    from hbb_pipeline.reference import HBBReference

logger = logging.getLogger(__name__)

# IUPAC two-base ambiguity code lookup (frozenset of two bases → IUPAC letter)
_IUPAC_FROM_PAIR: dict[frozenset, str] = {
    frozenset("AG"): "R",
    frozenset("CT"): "Y",
    frozenset("CG"): "S",
    frozenset("AT"): "W",
    frozenset("GT"): "K",
    frozenset("AC"): "M",
}


def apply_iupac_symbols(trace: TraceData, cutoff: float = 0.30) -> TraceData:
    """Pre-alignment IUPAC ambiguity coding via secondary peak detection.

    For every base position in the trace, reads raw channel intensities and
    applies a rule-based slope filter to identify genuine secondary peaks.
    Where a valid secondary peak is found, the basecaller's single-base call is
    replaced with the appropriate IUPAC ambiguity code.  This overrides the
    basecaller, catching heterozygous positions the basecaller silently called
    as the primary allele.

    Should be called on the raw trace immediately before trim_by_quality(),
    so that het positions (boosted to Q20) survive Mott trimming.

    Args:
        trace:  Raw TraceData (pre-trim).
        cutoff: Secondary peak ratio threshold (default 0.30).

    Returns:
        New TraceData with IUPAC-coded sequence.  All other fields unchanged.
    """
    from hbb_pipeline.heterozygosity import get_two_peaks

    new_seq: list[str] = []
    new_phred: list[int] = list(trace.phred_scores)
    for i in range(len(trace.sequence)):
        primary, secondary = get_two_peaks(trace, i, cutoff)
        if secondary is not None:
            pair = frozenset({primary.upper(), secondary.upper()})
            iupac = _IUPAC_FROM_PAIR.get(pair)
            if iupac:
                new_seq.append(iupac)
                # Basecaller often assigns low phred at het positions (mixed signal).
                # Override with Q20 so downstream quality filters don't suppress
                # a valid IUPAC call that passed the slope filter.
                if i < len(new_phred):
                    new_phred[i] = max(new_phred[i], 20)
                continue
        new_seq.append(primary)

    return TraceData(
        sequence="".join(new_seq),
        phred_scores=new_phred,
        channels=trace.channels,
        peak_positions=trace.peak_positions,
        sample_name=trace.sample_name,
        is_reverse_complemented=trace.is_reverse_complemented,
    )


class AlignmentFailedError(Exception):
    """Raised when alignment produces no usable result."""


# Complement map for reverse-complement — covers plain bases AND all IUPAC
# ambiguity codes so that IUPAC-coded bases in trace sequences (e.g. 'R', 'Y')
# are properly complemented during reverse_complement_trace.
# Without this, 'Y' on the reverse strand stays 'Y' instead of becoming 'R',
# causing build_consensus to see R vs Y as a disagreement → N at Q0.
_COMPLEMENT = str.maketrans(
    "ACGTRYSWKMBVDHNacgtryswkmbvdhn",
    "TGCAYRSWMKVBHDNtgcayrswmkvbhdn",
)


_IUPAC_HET_CODES = frozenset("RYSWKM")  # unambiguous het codes (excludes N/B/D/H/V)


def _is_iupac_het(base: str) -> bool:
    return base.upper() in _IUPAC_HET_CODES


def reverse_complement_trace(trace: TraceData) -> TraceData:
    """Return a new TraceData that is the reverse complement of the input.

    Transformations applied:
    - sequence: reverse + complement
    - phred_scores: reversed
    - peak_positions: reversed
    - channels: A↔T, C↔G, each array reversed
    - is_reverse_complemented: set to True

    Args:
        trace: Input TraceData (forward orientation).

    Returns:
        New TraceData in reverse-complement orientation.

    Example:
        rc = reverse_complement_trace(trace)
        # rc.sequence == str(Seq(trace.sequence).reverse_complement())
    """
    rc_seq = trace.sequence.translate(_COMPLEMENT)[::-1]
    rc_quals = trace.phred_scores[::-1]

    # Mirror scan coordinates so peak_positions remain ascending.
    # plot_chromatogram expects ascending peaks; reversing without remapping
    # produces descending peaks and an empty channel slice (x_min > x_max).
    # New scan index = (n_scans - 1 - original_scan), which maps the last
    # physical scan to position 0 and the first to position n_scans-1.
    orig = trace.channels
    n_scans = max((len(ch) for ch in orig.values() if ch), default=0)
    rc_peaks = [n_scans - 1 - p for p in reversed(trace.peak_positions)]

    # Reverse AND complement the channel arrays so that
    # rc_channels[base][rc_scan] == orig[complement[base]][n_scans-1-rc_scan].
    rc_channels: dict[str, list[int]] = {
        "A": list(reversed(orig.get("T", []))),
        "T": list(reversed(orig.get("A", []))),
        "C": list(reversed(orig.get("G", []))),
        "G": list(reversed(orig.get("C", []))),
    }

    return TraceData(
        sequence=rc_seq,
        phred_scores=rc_quals,
        channels=rc_channels,
        peak_positions=rc_peaks,
        sample_name=trace.sample_name,
        is_reverse_complemented=True,
    )


def align_to_reference(
    trace: TraceData,
    ref: "HBBReference",
) -> "AlignedRead":
    """Local-align a (possibly trimmed) trace to the HBB reference.

    Uses Bio.Align.PairwiseAligner in local mode with Sanger-tuned scoring.
    If the top two alignments are within 5% score, a warning is logged.

    Args:
        trace: Trimmed TraceData (forward or already RC'd for reverse reads).
        ref: Loaded HBBReference.

    Returns:
        AlignedRead with genomic coordinates, gap strings, per-base quality,
        and trace_index_map (alignment column → index in trace.sequence).

    Raises:
        AlignmentFailedError: if no alignment is produced.

    Example:
        aligned = align_to_reference(trimmed_fwd, ref)
        print(aligned.reference_start, aligned.reference_end)
    """
    aligner = Align.PairwiseAligner()
    aligner.mode = "local"
    aligner.match_score = 2
    aligner.mismatch_score = -3
    aligner.open_gap_score = -15
    aligner.extend_gap_score = -2

    query = trace.sequence.upper()
    target = ref.upper_seq

    alignments = list(aligner.align(target, query))
    if not alignments:
        raise AlignmentFailedError(
            f"No alignment produced for {trace.sample_name} against HBB reference"
        )

    best = alignments[0]

    # Warn if top-2 scores are ambiguous
    if len(alignments) > 1:
        s0 = best.score
        s1 = alignments[1].score
        if s1 >= 0.95 * s0:
            logger.warning(
                "%s: top two alignments within 5%% (%.1f vs %.1f) — result may be ambiguous",
                trace.sample_name,
                s0,
                s1,
            )

    # Extract aligned strings and coordinate mapping
    aligned_ref, aligned_query = _extract_aligned_strings(best)

    # Genomic bounds from alignment on target (reference)
    ref_coords = best.aligned[0]   # list of (start, end) pairs on target
    query_coords = best.aligned[1] # list of (start, end) pairs on query

    if len(ref_coords) == 0:
        raise AlignmentFailedError(f"Empty alignment for {trace.sample_name}")

    ref_start = int(ref_coords[0][0])
    ref_end = int(ref_coords[-1][1])

    # Build trace_index_map: for each column in aligned_query, the index into
    # the original trace.sequence (-1 for gap columns)
    trace_index_map = _build_trace_index_map(aligned_query, query_coords)

    # Per-base quality aligned to the query (gaps get quality 0)
    per_base_quality = _build_per_base_quality(aligned_query, trace.phred_scores, trace_index_map)

    logger.info(
        "%s: aligned [%d, %d), score=%.1f, %d alignment columns",
        trace.sample_name,
        ref_start,
        ref_end,
        best.score,
        len(aligned_ref),
    )

    return AlignedRead(
        read_name=trace.sample_name,
        reference_start=ref_start,
        reference_end=ref_end,
        aligned_seq=aligned_query,
        aligned_ref=aligned_ref,
        per_base_quality=per_base_quality,
        trace_index_map=trace_index_map,
    )


def build_consensus(
    fwd: "AlignedRead | None",
    rev: "AlignedRead | None",
    ref: "HBBReference",
) -> tuple[str, list[int]]:
    """Merge forward and/or reverse alignments into a consensus sequence.

    Accepts either or both strands; pass None for a missing strand.

    Iterates over every genomic position in the reference.  At each:
    1. Only fwd covers → use fwd base + qual
    2. Only rev covers → use rev base + qual
    3. Both agree → base + boosted quality (conservative cap Q60)
    4. Both disagree:
       - |q_fwd - q_rev| ≥ 10: use higher-Q base, qual = (higher - lower)
       - Otherwise: emit 'N', qual = 0 (flag caller must check)

    Returns:
        (consensus_str, per_genomic_qual) indexed from 0 to ref.length - 1.
        Positions not covered by either read get '?' and quality 0.

    Raises:
        ValueError: if both fwd and rev are None.

    Example:
        consensus, quals = build_consensus(fwd_aligned, rev_aligned, ref)
        consensus, quals = build_consensus(fwd_aligned, None, ref)  # single-strand
    """
    if fwd is None and rev is None:
        raise ValueError("build_consensus requires at least one aligned read")

    # Pre-compute per-genomic-position base+qual for each read
    fwd_bases, fwd_quals = _unpack_to_genomic(fwd) if fwd is not None else ({}, {})
    rev_bases, rev_quals = _unpack_to_genomic(rev) if rev is not None else ({}, {})

    consensus_list: list[str] = []
    qual_list: list[int] = []

    for gpos in range(ref.length):
        f_base = fwd_bases.get(gpos)
        r_base = rev_bases.get(gpos)
        fq = fwd_quals.get(gpos, 0)
        rq = rev_quals.get(gpos, 0)

        if f_base is None and r_base is None:
            consensus_list.append("?")
            qual_list.append(0)
        elif f_base is not None and r_base is None:
            consensus_list.append(f_base)
            qual_list.append(fq)
        elif f_base is None and r_base is not None:
            consensus_list.append(r_base)
            qual_list.append(rq)
        else:
            # Both cover this position
            if f_base == "-" and r_base == "-":
                # Both strands agree on a deletion — quality 30 passes Q20 caller filter
                consensus_list.append("-")
                qual_list.append(30)
            elif f_base.upper() == r_base.upper():
                # Agreement: conservative quality boost
                combined = min(60, max(fq, rq) + 3)
                consensus_list.append(f_base.upper())
                qual_list.append(combined)
            elif _is_iupac_het(f_base) and r_base.upper() == "N":
                # One strand has a specific IUPAC het code; other called N (ambiguous).
                # Trust the specific IUPAC call — N is strictly less informative than R/Y/etc.
                consensus_list.append(f_base.upper())
                qual_list.append(max(5, fq))
            elif _is_iupac_het(r_base) and f_base.upper() == "N":
                consensus_list.append(r_base.upper())
                qual_list.append(max(5, rq))
            else:
                # Disagreement (includes one-strand deletion vs base on other strand)
                if abs(fq - rq) >= 10:
                    if fq >= rq:
                        consensus_list.append(f_base.upper() if f_base != "-" else r_base.upper())
                    else:
                        consensus_list.append(r_base.upper() if r_base != "-" else f_base.upper())
                    qual_list.append(abs(fq - rq))
                else:
                    consensus_list.append("N")
                    qual_list.append(0)

    return "".join(consensus_list), qual_list


# ---------------------------------------------------------------------------
# Private helpers
# ---------------------------------------------------------------------------

def _extract_aligned_strings(alignment: Align.Alignment) -> tuple[str, str]:
    """Return (aligned_target, aligned_query) strings with '-' for gaps."""
    fmt = format(alignment)
    lines = fmt.strip().split("\n")
    # Biopython alignment format: target line, match line, query line
    # May have multiple blocks; reconstruct full strings
    target_parts: list[str] = []
    query_parts: list[str] = []

    i = 0
    while i < len(lines):
        line = lines[i]
        # Lines starting with "target" or digit (coordinate lines) delimit blocks
        if line.startswith("target"):
            parts = line.split()
            # parts[2] is the sequence string for target in this block
            if len(parts) >= 3:
                target_parts.append(parts[2])
        elif line.startswith("query"):
            parts = line.split()
            if len(parts) >= 3:
                query_parts.append(parts[2])
        i += 1

    # Fallback: use Biopython's substituted string method if parsing failed
    if not target_parts or not query_parts:
        return _aligned_strings_fallback(alignment)

    return "".join(target_parts), "".join(query_parts)


def _aligned_strings_fallback(alignment: Align.Alignment) -> tuple[str, str]:
    """Build aligned strings from alignment.aligned coordinate arrays."""
    target_seq = alignment.target
    query_seq = alignment.query
    ref_blocks = alignment.aligned[0]
    qry_blocks = alignment.aligned[1]

    t_out: list[str] = []
    q_out: list[str] = []

    prev_t_end = ref_blocks[0][0]
    prev_q_end = qry_blocks[0][0]

    for (ts, te), (qs, qe) in zip(ref_blocks, qry_blocks):
        # Gap in target between blocks
        t_gap = ts - prev_t_end
        q_gap = qs - prev_q_end
        if t_gap > 0 and q_gap == 0:
            t_out.append(str(target_seq[prev_t_end:ts]))
            q_out.append("-" * t_gap)
        elif q_gap > 0 and t_gap == 0:
            t_out.append("-" * q_gap)
            q_out.append(str(query_seq[prev_q_end:qs]))

        t_out.append(str(target_seq[ts:te]))
        q_out.append(str(query_seq[qs:qe]))
        prev_t_end = te
        prev_q_end = qe

    return "".join(t_out), "".join(q_out)


def _build_trace_index_map(aligned_query: str, query_coords: list) -> list[int]:
    """Map each column in aligned_query to an index in the original trace (-1 = gap)."""
    # query_coords: list of (start, end) pairs showing which query indices are aligned
    # We reconstruct the sequential query index at each column
    col_to_qidx: list[int] = []
    # Build a flat list of query indices in column order
    qidx_sequence: list[int] = []
    for qs, qe in query_coords:
        qidx_sequence.extend(range(qs, qe))

    qi = 0
    for col_char in aligned_query:
        if col_char == "-":
            col_to_qidx.append(-1)
        else:
            if qi < len(qidx_sequence):
                col_to_qidx.append(qidx_sequence[qi])
                qi += 1
            else:
                col_to_qidx.append(-1)

    return col_to_qidx


def _build_per_base_quality(
    aligned_query: str,
    phred_scores: list[int],
    trace_index_map: list[int],
) -> list[int]:
    """Build per-column quality array (0 for gap columns)."""
    quals: list[int] = []
    for col_char, tidx in zip(aligned_query, trace_index_map):
        if col_char == "-" or tidx < 0:
            quals.append(0)
        elif tidx < len(phred_scores):
            quals.append(phred_scores[tidx])
        else:
            quals.append(0)
    return quals


def _unpack_to_genomic(
    aligned: AlignedRead,
) -> tuple[dict[int, str], dict[int, int]]:
    """Unpack an AlignedRead into {genomic_pos: base} and {genomic_pos: qual} dicts.

    Gap columns ("-") in aligned_seq are recorded with quality 0 so
    build_consensus can distinguish "deletion present here" from "no coverage".
    Insertion columns (aligned_ref == "-") are skipped — they have no reference
    position, and advancing gpos for them would shift every downstream position.
    """
    bases: dict[int, str] = {}
    quals: dict[int, int] = {}

    gpos = aligned.reference_start
    for ref_char, col_char, qual in zip(
        aligned.aligned_ref, aligned.aligned_seq, aligned.per_base_quality
    ):
        if ref_char == "-":
            # Insertion in the query — no reference coordinate; skip.
            continue
        if col_char == "-":
            bases[gpos] = "-"
            quals[gpos] = 0
        else:
            bases[gpos] = col_char
            quals[gpos] = qual
        gpos += 1

    return bases, quals
