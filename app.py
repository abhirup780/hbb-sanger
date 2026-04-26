"""Streamlit UI for the HBB Sanger Sequencing Analysis Pipeline."""

from __future__ import annotations

import hashlib
import logging
import sys
import tempfile
from pathlib import Path

import streamlit as st

sys.path.insert(0, str(Path(__file__).parent))

from plots import plot_chromatogram, plot_coverage_map
from hbb_pipeline.alignment import (
    align_to_reference,
    apply_iupac_symbols,
    build_consensus,
    reverse_complement_trace,
)
from hbb_pipeline.coordinates import CoordinateTranslator
from hbb_pipeline.models import TraceData
from hbb_pipeline.parsing import parse_abi, trim_by_quality
from hbb_pipeline.reference import HBBReference
from hbb_pipeline.reporting import generate_report, render_markdown_report
from hbb_pipeline.tracy_wrapper import is_tracy_available, run_tracy_on_pair
from hbb_pipeline.variants import (
    call_variants_from_alignment,
    merge_variant_calls,
)

logger = logging.getLogger(__name__)
_DEFAULT_REF = Path(__file__).parent / "reference" / "HBB_reference.fasta"

# Bump this whenever pipeline code changes that affect cached results.
_PIPELINE_VERSION = "26"


# ---------------------------------------------------------------------------
# Pipeline helpers
# ---------------------------------------------------------------------------

def _slice_trace(trace: TraceData, s: int, e: int) -> TraceData:
    return TraceData(
        sequence=trace.sequence[s:e],
        phred_scores=trace.phred_scores[s:e],
        channels=trace.channels,
        peak_positions=trace.peak_positions[s:e],
        sample_name=trace.sample_name,
        is_reverse_complemented=trace.is_reverse_complemented,
    )


def _cds_coverage(consensus: str, ref: HBBReference) -> str:
    covered = sum(
        1 for gs, ge in ref.exons
        for pos in range(gs, ge)
        if pos < len(consensus) and consensus[pos] not in ("?", "N")
    )
    total = sum(ge - gs for gs, ge in ref.exons)
    return f"{100 * covered / total:.1f}%" if total else "N/A"


def _resolve_idx(genomic_pos: int, aligned) -> int:
    from hbb_pipeline.variants import _resolve_trace_index
    return _resolve_trace_index(genomic_pos, aligned)  # -1 = not covered


def _process_trace(
    raw_bytes: bytes,
    name: str,
    tmpdir: str,
    ref: HBBReference,
    is_reverse: bool,
) -> tuple[TraceData | None, TraceData | None, object | None, int, int, str | None]:
    """Parse, trim, (optionally RC), and align one trace.

    Returns (raw, trimmed, aligned, trim_s, trim_e, error_msg).
    On failure returns (None, None, None, 0, 0, error_message).
    """
    try:
        path = Path(tmpdir) / name
        path.write_bytes(raw_bytes)
        raw = parse_abi(path)

        # IUPAC coding BEFORE trimming: basecallers assign Q1 at heterozygous
        # positions (mixed signal). The phred boost in apply_iupac_symbols
        # raises those to Q20 so Mott's trimmer doesn't cut on het positions.
        iupac_raw = apply_iupac_symbols(raw)
        _, (s, e) = trim_by_quality(iupac_raw)
        trimmed = _slice_trace(iupac_raw, s, e)
        if is_reverse:
            trimmed = reverse_complement_trace(trimmed)

        aligned = align_to_reference(trimmed, ref)
        return raw, trimmed, aligned, s, e, None
    except Exception as exc:
        logger.warning("Trace %s failed: %s", name, exc)
        return None, None, None, 0, 0, f"{type(exc).__name__}: {exc}"


def _run_variant_pipeline(
    ref: HBBReference,
    translator: CoordinateTranslator,
    fwd_trimmed: TraceData | None,
    rev_rc: TraceData | None,
    fwd_aligned,
    rev_aligned,
    tracy_variants: list | None = None,
    tracy_was_run: bool = False,
) -> tuple[list, str]:
    """Run consensus → variant calling → merge.

    HET SNV detection relies on IUPAC codes already embedded in the traces
    by apply_iupac_symbols() (called in _process_trace before alignment).

    Returns (merged_variants, consensus_str).
    """
    consensus, cons_quals = build_consensus(fwd_aligned, rev_aligned, ref)

    alignment_variants = call_variants_from_alignment(
        consensus, cons_quals, ref, translator,
        fwd_trimmed, rev_rc, fwd_aligned, rev_aligned,
    )
    merged = merge_variant_calls(
        alignment_variants, tracy_variants or [], tracy_was_run=tracy_was_run,
    )
    return merged, consensus


# ---------------------------------------------------------------------------
# Cached pipeline entry points
# ---------------------------------------------------------------------------

@st.cache_data(show_spinner="Running analysis…")
def run_pipeline(
    fwd_bytes: bytes,
    rev_bytes: bytes,
    fwd_name: str,
    rev_name: str,
    ref_path_str: str,
    min_phred: int,
    het_ratio: float,
) -> dict:
    """Paired-trace pipeline with graceful single-strand degradation.

    If one trace fails to parse or align, the analysis continues on the
    remaining strand and a warning is surfaced in the result.
    """
    import traceback
    try:
        ref = HBBReference(Path(ref_path_str))
        translator = CoordinateTranslator(ref)
        warnings: list[str] = []

        with tempfile.TemporaryDirectory() as tmpdir:
            fwd_raw, fwd_trimmed, fwd_aligned, fwd_s, fwd_e, fwd_err = \
                _process_trace(fwd_bytes, fwd_name, tmpdir, ref, is_reverse=False)
            rev_raw, rev_rc, rev_aligned, rev_s, rev_e, rev_err = \
                _process_trace(rev_bytes, rev_name, tmpdir, ref, is_reverse=True)

            if fwd_err:
                warnings.append(f"Forward trace failed ({fwd_err}) — analysis uses reverse strand only")
            if rev_err:
                warnings.append(f"Reverse trace failed ({rev_err}) — analysis uses forward strand only")

            if fwd_aligned is None and rev_aligned is None:
                return {"error": "Both traces failed to process.\n" + "\n".join(warnings)}

            # --- Tracy paired decomposition (DEL/DUP) ---
            tracy_variants: list = []
            tracy_was_run = False
            if is_tracy_available():
                fwd_abi = Path(tmpdir) / fwd_name if not fwd_err else None
                rev_abi = Path(tmpdir) / rev_name if not rev_err else None
                try:
                    tracy_variants = run_tracy_on_pair(
                        fwd_abi, rev_abi,
                        Path(ref_path_str),
                        Path(tmpdir) / "tracy",
                        translator, ref,
                    )
                    tracy_was_run = True
                    logger.info("Tracy called %d variant(s)", len(tracy_variants))
                except Exception as exc:
                    logger.error("Tracy failed: %s", exc)
                    warnings.append(f"Tracy decompose failed ({exc}) — falling back to alignment-only calling")
            else:
                warnings.append(
                    "ℹ️ Heterozygous DEL/DUP detection is not available in this online version. "
                    "SNV calling and alignment-based indel detection are fully active. "
                    "For complete HET indel analysis, run the pipeline locally with tracy installed."
                )

            merged, consensus = _run_variant_pipeline(
                ref, translator,
                fwd_trimmed, rev_rc,
                fwd_aligned, rev_aligned,
                tracy_variants=tracy_variants,
                tracy_was_run=tracy_was_run,
            )

            qc: dict = {}
            if fwd_raw is not None:
                qc["mean_phred_fwd"] = round(
                    sum(fwd_raw.phred_scores) / max(len(fwd_raw.phred_scores), 1), 1)
                qc["usable_length_fwd"] = fwd_e - fwd_s
            else:
                qc["mean_phred_fwd"] = "N/A (trace failed)"
                qc["usable_length_fwd"] = "N/A"
            if rev_raw is not None:
                qc["mean_phred_rev"] = round(
                    sum(rev_raw.phred_scores) / max(len(rev_raw.phred_scores), 1), 1)
                qc["usable_length_rev"] = rev_e - rev_s
            else:
                qc["mean_phred_rev"] = "N/A (trace failed)"
                qc["usable_length_rev"] = "N/A"
            qc["cds_coverage_pct"] = _cds_coverage(consensus, ref)
            qc["tracy_available"] = str(is_tracy_available())
            if tracy_was_run:
                qc["tracy_variants_called"] = len(tracy_variants)
            if warnings:
                qc["analysis_warnings"] = "; ".join(warnings)

            sample_id = fwd_name.replace(".ab1", "").replace(".abi", "")
            report = generate_report(merged, qc, sample_id, consensus,
                                     tracy_used=tracy_was_run)
            md = render_markdown_report(report)

        return {
            "report": report,
            "markdown": md,
            "fwd_trace": fwd_trimmed,
            "rev_rc": rev_rc,
            "fwd_aligned": fwd_aligned,
            "rev_aligned": rev_aligned,
            "fwd_phred_raw": fwd_raw.phred_scores if fwd_raw else None,
            "rev_phred_raw": rev_raw.phred_scores if rev_raw else None,
            "warnings": warnings,
            "error": None,
        }
    except Exception as exc:
        return {"error": f"{type(exc).__name__}: {exc}\n\n{traceback.format_exc()}"}


@st.cache_data(show_spinner="Running analysis…")
def run_pipeline_single(
    trace_bytes: bytes,
    trace_name: str,
    is_reverse: bool,
    ref_path_str: str,
    min_phred: int,
    het_ratio: float,
) -> dict:
    """Single-trace pipeline — one forward or one reverse file."""
    import traceback
    try:
        ref = HBBReference(Path(ref_path_str))
        translator = CoordinateTranslator(ref)

        with tempfile.TemporaryDirectory() as tmpdir:
            raw, trace_proc, aligned, s, e, err = \
                _process_trace(trace_bytes, trace_name, tmpdir, ref, is_reverse=is_reverse)

            if err:
                return {"error": f"Trace failed to process: {err}"}

            if is_reverse:
                fwd_trimmed, rev_rc = None, trace_proc
                fwd_aligned, rev_aligned = None, aligned
                fwd_abi, rev_abi = None, Path(tmpdir) / trace_name
            else:
                fwd_trimmed, rev_rc = trace_proc, None
                fwd_aligned, rev_aligned = aligned, None
                fwd_abi, rev_abi = Path(tmpdir) / trace_name, None

            # --- Tracy single-strand decomposition (DEL/DUP) ---
            tracy_variants: list = []
            tracy_was_run = False
            single_warnings: list[str] = []
            if is_tracy_available():
                try:
                    tracy_variants = run_tracy_on_pair(
                        fwd_abi, rev_abi,
                        Path(ref_path_str),
                        Path(tmpdir) / "tracy",
                        translator, ref,
                    )
                    tracy_was_run = True
                    logger.info("Tracy (single-strand) called %d variant(s)", len(tracy_variants))
                except Exception as exc:
                    logger.error("Tracy failed: %s", exc)
                    single_warnings.append(f"Tracy decompose failed ({exc})")

            merged, consensus = _run_variant_pipeline(
                ref, translator,
                fwd_trimmed, rev_rc,
                fwd_aligned, rev_aligned,
                tracy_variants=tracy_variants,
                tracy_was_run=tracy_was_run,
            )

            strand_label = "reverse" if is_reverse else "forward"
            qc: dict = {
                f"mean_phred_{strand_label}": round(
                    sum(raw.phred_scores) / max(len(raw.phred_scores), 1), 1),
                f"usable_length_{strand_label}": e - s,
                "cds_coverage_pct": _cds_coverage(consensus, ref),
                "analysis_mode": f"Single strand ({strand_label})",
                "tracy_available": str(is_tracy_available()),
            }
            if tracy_was_run:
                qc["tracy_variants_called"] = len(tracy_variants)
            if single_warnings:
                qc["analysis_warnings"] = "; ".join(single_warnings)

            sample_id = trace_name.replace(".ab1", "").replace(".abi", "")
            report = generate_report(merged, qc, sample_id, consensus,
                                     tracy_used=tracy_was_run)
            md = render_markdown_report(report)

        return {
            "report": report,
            "markdown": md,
            "fwd_trace": fwd_trimmed,
            "rev_rc": rev_rc,
            "fwd_aligned": fwd_aligned,
            "rev_aligned": rev_aligned,
            "fwd_phred_raw": raw.phred_scores if not is_reverse else None,
            "rev_phred_raw": raw.phred_scores if is_reverse else None,
            "warnings": [
                f"Single-strand analysis ({strand_label} only). All variant calls require confirmation.",
                *single_warnings,
            ],
            "error": None,
        }
    except Exception as exc:
        return {"error": f"{type(exc).__name__}: {exc}\n\n{traceback.format_exc()}"}


# ---------------------------------------------------------------------------
# Main UI
# ---------------------------------------------------------------------------

def main() -> None:
    st.set_page_config(page_title="HBB Sanger", page_icon="🧬", layout="wide")
    st.title("🧬 HBB Sanger Analysis")
    st.caption("β-thalassemia & Sickle Cell")

    with st.sidebar:
        st.header("Trace Files")
        mode = st.radio("Mode", ["Paired (Fwd + Rev)", "Single Strand"],
                        horizontal=True, label_visibility="collapsed")

        if mode == "Paired (Fwd + Rev)":
            fwd_file = st.file_uploader("Forward (.ab1)", type=["ab1", "abi"])
            rev_file = st.file_uploader("Reverse (.ab1)", type=["ab1", "abi"])
            single_file = None
            is_reverse = False
            can_run = fwd_file is not None and rev_file is not None
        else:
            strand_dir = st.radio("Strand direction", ["Forward", "Reverse"],
                                  horizontal=True)
            is_reverse = strand_dir == "Reverse"
            single_file = st.file_uploader(
                f"{'Reverse' if is_reverse else 'Forward'} trace (.ab1)",
                type=["ab1", "abi"],
            )
            fwd_file = rev_file = None
            can_run = single_file is not None

        st.header("Reference")
        use_default = st.checkbox("Use bundled HBB reference", value=True)
        if use_default:
            ref_path_str = str(_DEFAULT_REF)
            st.caption(_DEFAULT_REF.name)
        else:
            ref_up = st.file_uploader("Custom FASTA", type=["fasta", "fa"])
            ref_path_str = None
            if ref_up:
                tmp = tempfile.NamedTemporaryFile(delete=False, suffix=".fasta")
                tmp.write(ref_up.read()); tmp.flush()
                ref_path_str = tmp.name

        st.header("QC Settings")
        min_phred = st.slider("Trim quality (Phred)", 10, 30, 20)
        het_ratio = st.slider("HET peak ratio", 0.15, 0.50, 0.25, step=0.01)

        st.divider()
        run_btn = st.button("▶ Run Analysis", type="primary",
                            use_container_width=True, disabled=not can_run)

    if not can_run:
        if mode == "Paired (Fwd + Rev)":
            st.info("Upload Forward and Reverse trace files, then click Run Analysis.")
        else:
            st.info(f"Upload a {'Reverse' if is_reverse else 'Forward'} trace file, then click Run Analysis.")
        return

    # Invalidate cache when files or pipeline version change
    if mode == "Paired (Fwd + Rev)":
        files_key = hashlib.md5(fwd_file.getvalue() + rev_file.getvalue()).hexdigest()
    else:
        files_key = hashlib.md5(single_file.getvalue() + bytes([int(is_reverse)])).hexdigest()

    if st.session_state.get("_files_key") != files_key:
        st.session_state["_files_key"] = files_key
        st.session_state.pop("result", None)
    if st.session_state.get("_pipeline_v") != _PIPELINE_VERSION:
        st.session_state["_pipeline_v"] = _PIPELINE_VERSION
        st.session_state.pop("result", None)
        run_pipeline.clear()
        run_pipeline_single.clear()

    if run_btn:
        if not ref_path_str:
            st.error("Select or upload a reference FASTA.")
        else:
            if mode == "Paired (Fwd + Rev)":
                result = run_pipeline(
                    fwd_file.getvalue(), rev_file.getvalue(),
                    fwd_file.name, rev_file.name,
                    ref_path_str, min_phred, het_ratio,
                )
            else:
                result = run_pipeline_single(
                    single_file.getvalue(), single_file.name,
                    is_reverse, ref_path_str, min_phred, het_ratio,
                )
            st.session_state["result"] = result

    if "result" not in st.session_state:
        return

    result = st.session_state["result"]

    if result.get("error"):
        st.error("Pipeline failed:")
        st.code(result["error"])
        return

    report    = result["report"]
    fwd_trace = result["fwd_trace"]   # None if fwd failed / not uploaded
    rev_rc    = result["rev_rc"]      # None if rev failed / not uploaded
    warnings  = result.get("warnings", [])

    # Surface pipeline warnings prominently
    for w in warnings:
        st.warning(w)

    tab_sum, tab_qc, tab_chrom, tab_trace, tab_report = st.tabs(
        ["Summary", "QC", "Variant Reviewer", "Trace Viewer", "Report"]
    )

    # ── Summary ──────────────────────────────────────────────────────────────
    with tab_sum:
        all_sorted = sorted(report.variants, key=lambda x: x.ref_pos_genomic)
        display_variants = [
            v for v in all_sorted
            if v.known_variant_name or v.region not in ("5UTR", "3UTR")
        ]
        utr_noise_count = len(all_sorted) - len(display_variants)
        n = len(display_variants)

        if n == 0:
            if utr_noise_count:
                st.success(
                    f"✓ No pathogenic variants detected. "
                    f"({utr_noise_count} background UTR haplotype variant(s) excluded from display)"
                )
            else:
                st.success("✓ No variants detected — wild type at all covered positions.")
        else:
            st.subheader(f"{n} variant{'s' if n > 1 else ''} detected")
            for v in display_variants:
                hgvs = f"HBB:{v.hgvs_c}"
                if v.hgvs_p_hgvs:
                    hgvs += f" {v.hgvs_p_hgvs}"
                zyg   = v.zygosity.value.capitalize()
                label = f"**{hgvs}** ({zyg})"
                if v.known_variant_name:
                    label += f"  —  {v.known_variant_name}"

                detail = f"Q{v.phred_support} · {v.region}"

                # Distinguish single-strand warnings from other review flags
                reason = v.review_reason or ""
                is_fwd_only = "Forward-strand only" in reason
                is_rev_only = "Reverse-strand only" in reason
                is_single   = is_fwd_only or is_rev_only

                if is_single:
                    badge = "⚠ Fwd only" if is_fwd_only else "⚠ Rev only"
                    detail += f" · {badge}"
                elif v.requires_manual_review:
                    detail += f" · ⚠ {reason or 'review required'}"

                if v.requires_manual_review:
                    st.warning(f"{label}  \n{detail}")
                elif v.known_variant_name:
                    st.error(f"🔴 {label}  \n{detail}")
                else:
                    st.info(f"{label}  \n{detail}")

        with st.expander("Full table"):
            import pandas as pd
            rows = [{
                "HGVS c.":   v.hgvs_c,
                "p. (HGVS)": v.hgvs_p_hgvs or "—",
                "Zygosity":  v.zygosity.value,
                "Region":    v.region,
                "Q":         v.phred_support,
                "Type":      v.variant_type,
                "Known":     v.known_variant_name or "—",
                "Called by": ", ".join(v.called_by),
                "Review":    (
                    "⚠ Fwd only" if "Forward-strand only" in (v.review_reason or "")
                    else "⚠ Rev only" if "Reverse-strand only" in (v.review_reason or "")
                    else "⚠" if v.requires_manual_review
                    else ""
                ),
            } for v in all_sorted]
            if utr_noise_count:
                st.caption(
                    f"{utr_noise_count} UTR haplotype variant(s) shown here but excluded from "
                    "the summary above — no clinical significance (reference haplotype mismatch)."
                )
            st.dataframe(pd.DataFrame(rows), width="stretch")


    # ── Trace Viewer ──────────────────────────────────────────────────────────
    with tab_trace:
        available = []
        if fwd_trace is not None:
            available.append("Forward")
        if rev_rc is not None:
            available.append("Reverse (RC)")

        if len(available) == 0:
            st.info("No trace data available.")
        else:
            if len(available) > 1:
                strand = st.radio("Strand", available, horizontal=True,
                                  label_visibility="collapsed")
            else:
                strand = available[0]
                st.caption(f"Showing {strand} strand only")

            trace_to_show = fwd_trace if strand == "Forward" else rev_rc

            import plotly.graph_objects as go
            _TV_COLORS = {"A": "#2ca02c", "C": "#1f77b4", "G": "#999999", "T": "#d62728"}
            channels = trace_to_show.channels
            n_scans  = max((len(ch) for ch in channels.values() if ch), default=1)
            stride   = max(1, n_scans // 3000)

            _all_sampled: list[int] = []
            for _ch in channels.values():
                if _ch:
                    _all_sampled.extend(_ch[i] for i in range(0, len(_ch), stride) if i < len(_ch))
            if _all_sampled:
                _all_sampled.sort()
                _y_cap = _all_sampled[min(int(len(_all_sampled) * 0.99), len(_all_sampled) - 1)]
            else:
                _y_cap = 1000
            _y_range = [0, max(_y_cap * 1.08, 100)]

            fig_tv = go.Figure()
            for base, color in _TV_COLORS.items():
                ch = channels.get(base, [])
                if not ch:
                    continue
                xs = list(range(0, n_scans, stride))
                ys = [ch[i] for i in xs if i < len(ch)]
                xs = xs[:len(ys)]
                fig_tv.add_trace(go.Scattergl(
                    x=xs, y=ys, mode="lines", name=base,
                    line=dict(color=color, width=1),
                ))

            peaks = trace_to_show.peak_positions
            seq   = trace_to_show.sequence
            label_step = max(1, len(seq) // 50)
            for i in range(0, len(seq), label_step):
                if i >= len(peaks):
                    break
                b = seq[i]
                fig_tv.add_annotation(
                    x=peaks[i], y=-0.06, xref="x", yref="paper",
                    text=f"<b>{b}</b>", showarrow=False,
                    font=dict(color=_TV_COLORS.get(b.upper(), "#888888"), size=9),
                )

            _AXIS = dict(color="#aaaaaa", gridcolor="rgba(120,120,120,0.2)",
                         zerolinecolor="rgba(120,120,120,0.3)")
            fig_tv.update_layout(
                height=500,
                margin=dict(l=10, r=10, t=36, b=50),
                xaxis=dict(
                    title="Scan",
                    rangeslider=dict(visible=True, thickness=0.08,
                                     bgcolor="rgba(30,30,30,0.95)",
                                     bordercolor="rgba(150,150,150,0.4)",
                                     borderwidth=1),
                    **_AXIS,
                ),
                yaxis=dict(title="Intensity", range=_y_range, fixedrange=False, **_AXIS),
                font=dict(color="#aaaaaa"),
                legend=dict(orientation="h", x=1, xanchor="right", y=1, yanchor="bottom",
                            bgcolor="rgba(0,0,0,0)", font=dict(color="#aaaaaa")),
                plot_bgcolor="rgba(0,0,0,0)",
                paper_bgcolor="rgba(0,0,0,0)",
                title=dict(text=trace_to_show.sample_name, x=0,
                           font=dict(size=12, color="#aaaaaa")),
            )
            st.plotly_chart(fig_tv, width="stretch")

    # ── Variant Reviewer ──────────────────────────────────────────────────────
    with tab_chrom:
        if not report.variants:
            st.info("No variants to inspect.")
        else:
            labels  = [f"{v.hgvs_c} ({v.region})" for v in report.variants]
            genomic = {f"{v.hgvs_c} ({v.region})": v.ref_pos_genomic for v in report.variants}

            sel  = st.selectbox("Variant to inspect", labels)
            gpos = genomic[sel]

            c1, c2 = st.columns(2)
            with c1:
                st.caption("Forward")
                if fwd_trace is None:
                    st.info("Forward trace not available.")
                else:
                    tidx = _resolve_idx(gpos, result["fwd_aligned"])
                    if tidx < 0:
                        st.info("Position not covered by forward read.")
                    else:
                        fig = plot_chromatogram(
                            fwd_trace,
                            max(0, tidx - 20),
                            min(len(fwd_trace.sequence), tidx + 21),
                            highlight_pos=tidx,
                        )
                        st.plotly_chart(fig, key="chrom_fwd", width="stretch")
            with c2:
                st.caption("Reverse (RC)")
                if rev_rc is None:
                    st.info("Reverse trace not available.")
                else:
                    tidx_r = _resolve_idx(gpos, result["rev_aligned"])
                    if tidx_r < 0:
                        st.info("Position not covered by reverse read.")
                    else:
                        fig = plot_chromatogram(
                            rev_rc,
                            max(0, tidx_r - 20),
                            min(len(rev_rc.sequence), tidx_r + 21),
                            highlight_pos=tidx_r,
                        )
                        st.plotly_chart(fig, key="chrom_rev", width="stretch")

    # ── QC ────────────────────────────────────────────────────────────────────
    with tab_qc:
        qc = report.qc_metrics

        fwd_al = result["fwd_aligned"]
        rev_al = result["rev_aligned"]
        fig_cov = plot_coverage_map(
            fwd_al.reference_start if fwd_al else None,
            fwd_al.reference_end   if fwd_al else None,
            rev_al.reference_start if rev_al else None,
            rev_al.reference_end   if rev_al else None,
        )
        st.plotly_chart(fig_cov, width="stretch")

        c1, c2, c3, c4, c5 = st.columns(5)
        c1.metric("Mean Q (Fwd)", qc.get("mean_phred_fwd", "—"))
        c2.metric("Mean Q (Rev)", qc.get("mean_phred_rev", "—"))
        c3.metric("CDS Coverage", qc.get("cds_coverage_pct", "—"))
        c4.metric("Usable (Fwd)", f"{qc.get('usable_length_fwd', '—')} bp")
        c5.metric("Usable (Rev)", f"{qc.get('usable_length_rev', '—')} bp")

        import plotly.graph_objects as go
        fig_q = go.Figure()
        fwd_phred = result.get("fwd_phred_raw")
        rev_phred = result.get("rev_phred_raw")
        if fwd_phred:
            fig_q.add_trace(go.Histogram(x=fwd_phred, name="Forward",
                                         opacity=0.7, nbinsx=40,
                                         marker_color="#1f77b4"))
        if rev_phred:
            fig_q.add_trace(go.Histogram(x=rev_phred, name="Reverse",
                                         opacity=0.7, nbinsx=40,
                                         marker_color="#d62728"))
        if fwd_phred or rev_phred:
            fig_q.update_layout(
                barmode="overlay", height=280, title="Phred distribution",
                xaxis_title="Phred", yaxis_title="Count",
                plot_bgcolor="rgba(0,0,0,0)", paper_bgcolor="rgba(0,0,0,0)",
                font=dict(color="#aaaaaa"),
                xaxis=dict(color="#aaaaaa", gridcolor="rgba(120,120,120,0.2)"),
                yaxis=dict(color="#aaaaaa", gridcolor="rgba(120,120,120,0.2)"),
            )
            st.plotly_chart(fig_q, width="stretch")

    # ── Report ────────────────────────────────────────────────────────────────
    with tab_report:
        st.markdown(result["markdown"])
        st.download_button(
            "⬇ Download report (.md)",
            data=result["markdown"],
            file_name=f"{report.sample_id}_HBB_report.md",
            mime="text/markdown",
        )


if __name__ == "__main__":
    main()
