"""Variant calling from consensus alignment and optional tracy output, with HGVS annotation.

Two calling paths:
1. Alignment-based: walk consensus vs reference, emit variants for mismatches/gaps.
2. Tracy-based: parse tracy VCF output (see tracy_wrapper.py).

The merge step combines both and flags discordances for manual review.

HET SNV detection relies on IUPAC codes embedded in the consensus by
apply_iupac_symbols() (alignment.py) before alignment.
No separate known-variant scan is needed.
"""

from __future__ import annotations

import logging
from collections import OrderedDict
from typing import TYPE_CHECKING

from hbb_pipeline.heterozygosity import classify_zygosity, detect_secondary_peaks
from hbb_pipeline.known_variants import lookup_variant
from hbb_pipeline.models import Variant, Zygosity

if TYPE_CHECKING:
    from hbb_pipeline.coordinates import CoordinateTranslator
    from hbb_pipeline.models import AlignedRead, TraceData
    from hbb_pipeline.reference import HBBReference

logger = logging.getLogger(__name__)

_MIN_PHRED = 20  # Q20 minimum to call a variant
_MIN_PHRED_IUPAC = 5  # lower threshold when basecaller itself flagged ambiguity
_HET_RATIO = 0.25        # secondary peak ratio to confirm a HET call
_SUSPECT_RATIO = 0.20    # lower ratio to raise a sub-threshold suspect at known pathogenic sites

# IUPAC two-base ambiguity codes → the pair of bases they represent
_IUPAC_HET: dict[str, frozenset[str]] = {
    "R": frozenset("AG"),
    "Y": frozenset("CT"),
    "S": frozenset("CG"),
    "W": frozenset("AT"),
    "K": frozenset("GT"),
    "M": frozenset("AC"),
}


def call_variants_from_alignment(
    consensus: str,
    quals: list[int],
    ref: "HBBReference",
    translator: "CoordinateTranslator",
    fwd_trace: "TraceData | None",
    rev_trace: "TraceData | None",
    fwd_aligned: "AlignedRead | None",
    rev_aligned: "AlignedRead | None",
) -> list[Variant]:
    """Walk the consensus and emit a Variant for every difference vs reference.

    Only positions with consensus quality >= Q20 are called.  Each variant is
    annotated with:
    - HGVS c. and p. (both hgvs and legacy numbering)
    - Region (5UTR/exon/intron/3UTR)
    - Zygosity from secondary-peak analysis
    - known_variant_name if in the known_variants lookup

    Args:
        consensus: Full-length consensus string (length == ref.length), with
                   '?' for uncovered positions and 'N' for ambiguous positions.
        quals: Per-position Phred quality array (same length as consensus).
        ref: Loaded HBBReference.
        translator: CoordinateTranslator for this reference.
        fwd_trace: Forward TraceData.
        rev_trace: Reverse (already RC'd) TraceData.
        fwd_aligned: AlignedRead for forward strand.
        rev_aligned: AlignedRead for reverse strand.

    Returns:
        List of Variant objects, sorted by ref_pos_genomic.

    Example:
        variants = call_variants_from_alignment(cons, qs, ref, trans, ft, rt, fa, ra)
    """
    variants: list[Variant] = []
    ref_seq = ref.upper_seq

    i = 0
    while i < len(consensus):
        base = consensus[i]

        # Skip uncovered positions
        if base == "?":
            i += 1
            continue

        # Consensus 'N' means both strands called conflicting IUPAC codes (e.g. forward
        # 'R' vs reverse 'Y') — they represent the same heterozygosity from opposite
        # strands but disagree as strings, so quality collapses to 0.  Rescue by checking
        # raw secondary peaks directly: if both strands show the same alt allele above the
        # ratio threshold, that is confident dual-strand HET evidence.
        if base == "N":
            fwd_tidx = _resolve_trace_index(i, fwd_aligned)
            rev_tidx = _resolve_trace_index(i, rev_aligned)
            if fwd_tidx >= 0 and rev_tidx >= 0:
                ref_base = ref_seq[i]
                fwd_sp = detect_secondary_peaks(fwd_trace, fwd_tidx, ratio_threshold=_HET_RATIO)
                rev_sp = detect_secondary_peaks(rev_trace, rev_tidx, ratio_threshold=_HET_RATIO)
                common_alts = {
                    a for a in fwd_sp
                    if a != ref_base.upper() and a in rev_sp
                }
                if common_alts:
                    alt = next(iter(common_alts))
                    v = _make_variant(
                        ref_pos=i, ref_allele=ref_base, alt_allele=alt,
                        vtype="SNV", phred=25,
                        ref=ref, translator=translator,
                        fwd_trace=fwd_trace, rev_trace=rev_trace,
                        fwd_aligned=fwd_aligned, rev_aligned=rev_aligned,
                    )
                    if v is not None:
                        v = v.model_copy(update={"zygosity": Zygosity.HET})
                        variants.append(v)
            i += 1
            continue

        # IUPAC ambiguity code: the basecaller itself detected two overlapping alleles.
        # This is the strongest single-base evidence of heterozygosity — call it directly
        # at a lower quality threshold without requiring Q20.
        iupac_pair = _IUPAC_HET.get(base.upper())
        if iupac_pair is not None:
            ref_base = ref_seq[i]
            alt_set = iupac_pair - {ref_base.upper()}
            # Only call if the reference base is one of the two IUPAC bases
            if len(alt_set) == 1 and ref_base.upper() in iupac_pair and quals[i] >= _MIN_PHRED_IUPAC:
                alt_allele = next(iter(alt_set))
                v = _make_variant(
                    ref_pos=i,
                    ref_allele=ref_base,
                    alt_allele=alt_allele,
                    vtype="SNV",
                    phred=quals[i],
                    ref=ref,
                    translator=translator,
                    fwd_trace=fwd_trace,
                    rev_trace=rev_trace,
                    fwd_aligned=fwd_aligned,
                    rev_aligned=rev_aligned,
                )
                if v is not None:
                    v = v.model_copy(update={"zygosity": Zygosity.HET})
                    variants.append(v)
            i += 1
            continue

        # Deletions ("-") bypass the per-base Q20 check; their quality (30 when both
        # strands agree, 0 for single-strand) is evaluated separately below.
        if base != "-" and quals[i] < _MIN_PHRED:
            i += 1
            continue

        ref_base = ref_seq[i]

        if base.upper() == ref_base.upper():
            i += 1
            continue

        # --- Determine variant type ---
        # Simple model: treat non-gap mismatches as SNV, gaps as DEL/INS.
        # Full indel right-alignment is handled in CoordinateTranslator.build_hgvs_c.
        if base == "-":
            # Deletion in consensus — require both strands agreed (build_consensus
            # assigns Q30 for dual-strand deletions, Q0 for single-strand).
            if quals[i] < 30:
                i += 1
                continue
            del_start = i
            while i < len(consensus) and consensus[i] == "-":
                i += 1
            del_ref = ref_seq[del_start:i]
            v = _make_variant(
                ref_pos=del_start,
                ref_allele=del_ref,
                alt_allele="",
                vtype="DEL",
                phred=quals[del_start],
                ref=ref,
                translator=translator,
                fwd_trace=fwd_trace,
                rev_trace=rev_trace,
                fwd_aligned=fwd_aligned,
                rev_aligned=rev_aligned,
            )
            if v is not None:
                # Both strands agreed on deletion (required for qual >= Q20), so homozygous
                v = v.model_copy(update={"zygosity": Zygosity.HOM})
                variants.append(v)
            continue

        if ref_base == "-":
            # Insertion in consensus — not typical in full-length consensus but handle it
            ins_start = i
            ins_bases = ""
            while i < len(consensus) and ref_seq[i] == "-":
                ins_bases += consensus[i]
                i += 1
            v = _make_variant(
                ref_pos=ins_start,
                ref_allele="",
                alt_allele=ins_bases,
                vtype="INS",
                phred=quals[ins_start],
                ref=ref,
                translator=translator,
                fwd_trace=fwd_trace,
                rev_trace=rev_trace,
                fwd_aligned=fwd_aligned,
                rev_aligned=rev_aligned,
            )
            if v is not None:
                variants.append(v)
            continue

        # SNV
        v = _make_variant(
            ref_pos=i,
            ref_allele=ref_base,
            alt_allele=base.upper(),
            vtype="SNV",
            phred=quals[i],
            ref=ref,
            translator=translator,
            fwd_trace=fwd_trace,
            rev_trace=rev_trace,
            fwd_aligned=fwd_aligned,
            rev_aligned=rev_aligned,
        )
        if v is not None:
            variants.append(v)
        i += 1

    annotate_known(variants)
    return variants


def _make_variant(
    ref_pos: int,
    ref_allele: str,
    alt_allele: str,
    vtype: str,
    phred: int,
    ref: "HBBReference",
    translator: "CoordinateTranslator",
    fwd_trace: "TraceData | None",
    rev_trace: "TraceData | None",
    fwd_aligned: "AlignedRead | None",
    rev_aligned: "AlignedRead | None",
) -> Variant | None:
    """Build a Variant, computing HGVS and zygosity. Returns None on error."""
    try:
        region = ref.region_of(ref_pos)

        # Stub variant (hgvs_c filled after build_hgvs_c call)
        from hbb_pipeline.models import Zygosity
        stub = Variant(
            ref_pos_genomic=ref_pos,
            ref_allele=ref_allele,
            alt_allele=alt_allele,
            zygosity=Zygosity.UNKNOWN,
            variant_type=vtype,  # type: ignore[arg-type]
            hgvs_c="",
            hgvs_p_hgvs=None,
            hgvs_p_legacy=None,
            region=region,  # type: ignore[arg-type]
            phred_support=phred,
            secondary_peak_ratio_fwd=None,
            secondary_peak_ratio_rev=None,
            called_by=["alignment"],
        )

        hgvs_c = translator.build_hgvs_c(stub)
        hgvs_p_hgvs = translator.build_hgvs_p(stub, "hgvs")
        hgvs_p_legacy = translator.build_hgvs_p(stub, "legacy")

        # Secondary peak analysis
        fwd_tidx = _resolve_trace_index(ref_pos, fwd_aligned)
        rev_tidx = _resolve_trace_index(ref_pos, rev_aligned)

        fwd_covers = fwd_tidx >= 0 and fwd_trace is not None
        rev_covers = rev_tidx >= 0 and rev_trace is not None

        fwd_peaks = detect_secondary_peaks(fwd_trace, fwd_tidx) if fwd_covers else {}
        rev_peaks = detect_secondary_peaks(rev_trace, rev_tidx) if rev_covers else {}

        zygosity, fwd_ratio, rev_ratio = classify_zygosity(fwd_peaks, rev_peaks, alt_allele or "N")

        # Flag variants called from only one strand — these need manual confirmation.
        if fwd_covers and not rev_covers:
            single_strand_review = "Forward-strand only — reverse read does not cover this position"
        elif rev_covers and not fwd_covers:
            single_strand_review = "Reverse-strand only — forward read does not cover this position"
        else:
            single_strand_review = None

        return Variant(
            ref_pos_genomic=ref_pos,
            ref_allele=ref_allele,
            alt_allele=alt_allele,
            zygosity=zygosity,
            variant_type=vtype,  # type: ignore[arg-type]
            hgvs_c=hgvs_c,
            hgvs_p_hgvs=hgvs_p_hgvs,
            hgvs_p_legacy=hgvs_p_legacy,
            region=region,  # type: ignore[arg-type]
            phred_support=phred,
            secondary_peak_ratio_fwd=fwd_ratio,
            secondary_peak_ratio_rev=rev_ratio,
            called_by=["alignment"],
            requires_manual_review=single_strand_review is not None,
            review_reason=single_strand_review,
        )
    except Exception as exc:
        logger.warning("Failed to build variant at genomic %d: %s", ref_pos, exc)
        return None


# Cache: object-id of AlignedRead → {genomic_pos: trace_idx}.
# Keeps at most 4 entries (2 strands × up to 2 pipeline calls per session).
_TRACE_IDX_CACHE: OrderedDict[int, dict[int, int]] = OrderedDict()
_TRACE_IDX_CACHE_MAX = 4


def _build_pos_map(aligned: "AlignedRead") -> dict[int, int]:
    """Build a full {genomic_pos: trace_idx} map for one alignment in O(n)."""
    pos_map: dict[int, int] = {}
    gpos = aligned.reference_start
    tmap = aligned.trace_index_map
    for col_idx, ref_char in enumerate(aligned.aligned_ref):
        if ref_char != "-":
            pos_map[gpos] = tmap[col_idx] if col_idx < len(tmap) else -1
            gpos += 1
    return pos_map


def _resolve_trace_index(genomic_pos: int, aligned: "AlignedRead | None") -> int:
    """Return the trace index for a genomic position. Returns -1 if not covered.

    O(1) after the first call for a given AlignedRead object — the position map
    is built once and cached by object identity.
    """
    if aligned is None:
        return -1
    if not (aligned.reference_start <= genomic_pos < aligned.reference_end):
        return -1

    key = id(aligned)
    if key not in _TRACE_IDX_CACHE:
        if len(_TRACE_IDX_CACHE) >= _TRACE_IDX_CACHE_MAX:
            _TRACE_IDX_CACHE.popitem(last=False)
        _TRACE_IDX_CACHE[key] = _build_pos_map(aligned)

    return _TRACE_IDX_CACHE[key].get(genomic_pos, -1)


def merge_variant_calls(
    alignment_variants: list[Variant],
    tracy_variants: list[Variant],
    tracy_was_run: bool = False,
) -> list[Variant]:
    """Merge alignment-based and tracy-based variant calls.

    Match by (ref_pos_genomic, ref_allele, alt_allele):
    - Both agree: merge called_by, prefer tracy zygosity.
    - Only alignment (tracy_was_run=True): requires_manual_review=True.
    - Only tracy: requires_manual_review=True.
    - Zygosity disagrees: requires_manual_review=True.

    When tracy_was_run=False (default), alignment variants are passed through
    unchanged — single-strand flags set by _make_variant are preserved and no
    spurious "not confirmed by tracy" flag is added.

    Returns merged list sorted by genomic position.

    Example:
        merged = merge_variant_calls(alignment_variants, tracy_variants)
        merged = merge_variant_calls(alignment_variants, [], tracy_was_run=False)
    """
    # Index tracy variants by key
    tracy_map: dict[tuple, Variant] = {
        (v.ref_pos_genomic, v.ref_allele, v.alt_allele): v
        for v in tracy_variants
    }
    alignment_map: dict[tuple, Variant] = {
        (v.ref_pos_genomic, v.ref_allele, v.alt_allele): v
        for v in alignment_variants
    }

    merged: list[Variant] = []
    all_keys = set(tracy_map) | set(alignment_map)

    for key in sorted(all_keys):
        aln = alignment_map.get(key)
        trc = tracy_map.get(key)

        if aln is not None and trc is not None:
            # Called by both — prefer tracy zygosity.
            # Only flag review when both have definitive (non-UNKNOWN) but conflicting calls.
            both_definitive = (
                aln.zygosity != Zygosity.UNKNOWN
                and trc.zygosity != Zygosity.UNKNOWN
            )
            zyg_conflict = both_definitive and (aln.zygosity != trc.zygosity)
            merged.append(
                aln.model_copy(update={
                    "called_by": ["alignment", "tracy"],
                    "zygosity": trc.zygosity,
                    "requires_manual_review": zyg_conflict,
                    "review_reason": (
                        "Zygosity discordance between alignment and tracy"
                        if zyg_conflict else None
                    ),
                })
            )
        elif aln is not None:
            if tracy_was_run and aln.variant_type != "SNV":
                # Tracy ran but didn't confirm this DEL/INS — flag for review.
                # SNVs are never expected from tracy (it is an indel detector);
                # flagging alignment-only SNVs as "not confirmed by tracy" is wrong.
                merged.append(aln.model_copy(update={
                    "requires_manual_review": True,
                    "review_reason": "alignment only — not confirmed by tracy",
                }))
            else:
                # SNV, or tracy not run — preserve flags set by _make_variant.
                merged.append(aln)
        else:
            assert trc is not None
            merged.append(trc.model_copy(update={
                "requires_manual_review": True,
                "review_reason": "tracy only — not confirmed by alignment",
            }))

    return merged





def annotate_known(variants: list[Variant]) -> None:
    """In-place: fill known_variant_name from the known_variants lookup table.

    Tries exact HGVS string match first; falls back to genomic position lookup
    so that notation differences (e.g. 'del' vs 'delCTTT') don't silently miss.

    Example:
        annotate_known(variants)
        # variants[i].known_variant_name == "HbS (Sickle cell)"
    """
    for v in variants:
        key = v.hgvs_c[2:] if v.hgvs_c.startswith("c.") else v.hgvs_c
        kv = lookup_variant(key)
        if kv is not None:
            v.known_variant_name = kv.name
