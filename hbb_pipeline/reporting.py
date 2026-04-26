"""Clinical report generation and Markdown rendering.

All reports include a mandatory disclaimer.  Never auto-diagnose.
"""

from __future__ import annotations

import logging
from datetime import datetime
from typing import TYPE_CHECKING

from hbb_pipeline.models import ClinicalReport, Zygosity

if TYPE_CHECKING:
    from hbb_pipeline.models import Variant

logger = logging.getLogger(__name__)

_DISCLAIMER = (
    "For research use only. Not validated for clinical diagnostic use. "
    "All variants, particularly those marked for manual review, must be confirmed "
    "by a certified clinical laboratory before any clinical decision is made."
)


def generate_report(
    variants: list["Variant"],
    qc_metrics: dict,
    sample_id: str,
    consensus_seq: str,
    reference_version: str = "HBB_reference_2750bp",
    tracy_used: bool = False,
) -> ClinicalReport:
    """Assemble a ClinicalReport from pipeline outputs."""
    return ClinicalReport(
        sample_id=sample_id,
        analysis_timestamp=datetime.now(),
        variants=variants,
        consensus_seq=consensus_seq,
        qc_metrics=qc_metrics,
        reference_version=reference_version,
        tracy_used=tracy_used,
    )


def render_markdown_report(report: ClinicalReport) -> str:
    """Render a ClinicalReport to a Markdown string.

    Sections:
    1. Sample Info & Timestamp
    2. QC Summary
    3. Variants Table
    4. Clinically Significant Findings
    5. Variants Requiring Manual Review
    6. Methods
    7. Disclaimer

    Args:
        report: ClinicalReport to render.

    Returns:
        Full Markdown string.

    Example:
        md = render_markdown_report(report)
        Path("report.md").write_text(md)
    """
    lines: list[str] = []

    # --- Header ---
    lines.append("# HBB Sanger Sequencing Analysis Report\n")

    # --- 1. Sample Info ---
    lines.append("## 1. Sample Information\n")
    lines.append(f"| Field | Value |")
    lines.append(f"|---|---|")
    lines.append(f"| Sample ID | {report.sample_id} |")
    lines.append(f"| Analysis timestamp | {report.analysis_timestamp.strftime('%Y-%m-%d %H:%M:%S')} |")
    lines.append(f"| Reference | {report.reference_version} |")
    lines.append("")

    # --- 2. QC Summary ---
    lines.append("## 2. QC Summary\n")
    qc = report.qc_metrics
    lines.append("| Metric | Value |")
    lines.append("|---|---|")
    for key, val in qc.items():
        lines.append(f"| {key.replace('_', ' ').title()} | {val} |")
    lines.append("")

    # --- 3. Variants Table ---
    lines.append("## 3. All Variants\n")
    if not report.variants:
        lines.append("*No variants detected above Q20.*\n")
    else:
        lines.append(
            "| Genomic pos | HGVS c. | HGVS p. (HGVS) | HGVS p. (legacy) | "
            "Region | Zygosity | Called by | Fwd ratio | Rev ratio | Known variant |"
        )
        lines.append("|---|---|---|---|---|---|---|---|---|---|")
        for v in sorted(report.variants, key=lambda x: x.ref_pos_genomic):
            p_hgvs = v.hgvs_p_hgvs or "—"
            p_leg = v.hgvs_p_legacy or "—"
            fwd_r = f"{v.secondary_peak_ratio_fwd:.2f}" if v.secondary_peak_ratio_fwd is not None else "—"
            rev_r = f"{v.secondary_peak_ratio_rev:.2f}" if v.secondary_peak_ratio_rev is not None else "—"
            known = v.known_variant_name or "—"
            called = ", ".join(v.called_by)
            flag = " ⚠" if v.requires_manual_review else ""
            lines.append(
                f"| {v.ref_pos_genomic} | {v.hgvs_c} | {p_hgvs} | {p_leg} | "
                f"{v.region} | {v.zygosity.value} | {called} | "
                f"{fwd_r} | {rev_r} | {known}{flag} |"
            )
    lines.append("")

    # --- 4. Clinically Significant Findings ---
    lines.append("## 4. Clinically Significant Findings\n")
    known_vars = [v for v in report.variants if v.known_variant_name]
    if not known_vars:
        lines.append("*No known pathogenic variants detected.*\n")
    else:
        lines.append("> **The following known pathogenic variants were identified:**\n")
        for v in known_vars:
            zyg_str = v.zygosity.value.capitalize()
            lines.append(
                f"- **{v.hgvs_c}** ({v.known_variant_name}) — "
                f"{zyg_str} | HGVS p.: {v.hgvs_p_hgvs or '—'} "
                f"(legacy: {v.hgvs_p_legacy or '—'})"
            )
        lines.append("")
        lines.append(
            "> ⚠ These findings require confirmation by a certified clinical laboratory."
        )
    lines.append("")

    # --- 5. Variants Requiring Manual Review ---
    lines.append("## 5. Variants Requiring Manual Review\n")
    flagged = [v for v in report.variants if v.requires_manual_review]
    if not flagged:
        lines.append("*No variants require manual review.*\n")
    else:
        lines.append("| HGVS c. | Zygosity | Reason |")
        lines.append("|---|---|---|")
        for v in flagged:
            lines.append(f"| {v.hgvs_c} | {v.zygosity.value} | {v.review_reason or '—'} |")
    lines.append("")

    # --- 6. Methods ---
    lines.append("## 6. Methods\n")
    lines.append(
        "Sanger trace files (.ab1) were parsed using Biopython. "
        "Reads were quality-trimmed using Mott's modified algorithm (threshold Q20). "
        "Trimmed reads were aligned to the case-annotated HBB reference "
        f"({report.reference_version}, 2750 bp) using BioPython PairwiseAligner "
        "(local mode; match=2, mismatch=−3, gap open=−15, gap extend=−2). "
        "Forward and reverse alignments were merged into a consensus; positions of "
        "disagreement with |ΔQ| < 10 were called 'N' and flagged for review. "
        "Variants were annotated with HGVS c. and p. nomenclature using an internal "
        "coordinate translator against the reference coordinate map. "
        "Zygosity was determined from secondary peak ratios (threshold 0.25). "
        + (
            "Heterozygous indel detection (DEL/DUP) was performed using tracy decompose "
            "(TIDE-based decomposition, bioconda/tracy); forward and reverse traces were "
            "decomposed independently and results merged by strand concordance. "
            if report.tracy_used else
            "Tracy (TIDE-based DEL/DUP decomposition) was not available at runtime; "
            "indel detection relies on alignment-based calling only. "
        )
    )
    lines.append("")

    # --- 7. Disclaimer ---
    lines.append("## 7. Disclaimer\n")
    lines.append(f"*{_DISCLAIMER}*\n")

    return "\n".join(lines)
