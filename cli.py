"""CLI entry point for the HBB Sanger analysis pipeline.

Usage:
    hbb-analyze run fwd.ab1 rev.ab1 --reference reference/HBB_reference.fasta --out report.md
    hbb-analyze run fwd.ab1 rev.ab1 --json out.json
    hbb-analyze validate-reference reference/HBB_reference.fasta
"""

from __future__ import annotations

import json
import logging
import sys
import tempfile
from pathlib import Path

import typer

# Ensure hbb_pipeline importable when run from hbb_sanger/
sys.path.insert(0, str(Path(__file__).parent))

app = typer.Typer(
    name="hbb-analyze",
    help="HBB Sanger sequencing analysis — beta-thalassemia & sickle cell variant detection.",
    add_completion=False,
)

logging.basicConfig(
    level=logging.INFO,
    format="%(levelname)s %(name)s: %(message)s",
)
logger = logging.getLogger(__name__)

_DEFAULT_REF = Path(__file__).parent / "reference" / "HBB_reference.fasta"


@app.command()
def run(
    fwd: Path = typer.Argument(..., help="Forward .ab1 trace file"),
    rev: Path = typer.Argument(..., help="Reverse .ab1 trace file"),
    reference: Path = typer.Option(_DEFAULT_REF, "--reference", "-r", help="Reference FASTA"),
    out: Path | None = typer.Option(None, "--out", "-o", help="Output Markdown report path"),
    json_out: Path | None = typer.Option(None, "--json", help="Output JSON report path"),
    min_phred: int = typer.Option(20, "--min-phred", help="Minimum Phred for quality trimming"),
) -> None:
    """Run the full HBB analysis pipeline on a forward + reverse trace pair."""
    from hbb_pipeline.alignment import (
        align_to_reference,
        apply_iupac_symbols,
        build_consensus,
        reverse_complement_trace,
    )
    from hbb_pipeline.coordinates import CoordinateTranslator
    from hbb_pipeline.parsing import parse_abi, trim_by_quality
    from hbb_pipeline.reference import HBBReference, ReferenceValidationError
    from hbb_pipeline.reporting import generate_report, render_markdown_report
    from hbb_pipeline.tracy_wrapper import is_tracy_available, run_tracy_on_pair, tracy_version
    from hbb_pipeline.variants import call_variants_from_alignment, merge_variant_calls

    # Validate reference
    typer.echo(f"Loading reference: {reference}")
    try:
        ref = HBBReference(reference)
    except ReferenceValidationError as exc:
        typer.echo(f"ERROR: Reference validation failed: {exc}", err=True)
        raise typer.Exit(1)
    typer.echo(f"  ✓ Reference OK — {ref.length} bp, CDS {len(ref.cds)} bp")

    translator = CoordinateTranslator(ref)

    # Parse traces
    typer.echo(f"Parsing: {fwd.name} (fwd), {rev.name} (rev)")
    from hbb_pipeline.parsing import InvalidTraceFileError

    try:
        fwd_raw = parse_abi(fwd)
        rev_raw = parse_abi(rev)
    except InvalidTraceFileError as exc:
        typer.echo(f"ERROR: {exc}", err=True)
        raise typer.Exit(1)

    # IUPAC coding BEFORE trimming: basecallers assign Q1 at heterozygous
    # positions (mixed signal). The phred boost raises those to Q20 so Mott's
    # trimmer doesn't cut on het positions and discard genuine variant evidence.
    fwd_iupac_raw = apply_iupac_symbols(fwd_raw)
    rev_iupac_raw = apply_iupac_symbols(rev_raw)

    threshold = 10 ** (-min_phred / 10)
    fwd_iupac, (fs, fe) = trim_by_quality(fwd_iupac_raw, threshold=threshold)
    rev_iupac_trimmed, (rs, re) = trim_by_quality(rev_iupac_raw, threshold=threshold)
    typer.echo(
        f"  Fwd trimmed: {fe - fs} bp  |  Rev trimmed: {re - rs} bp"
    )

    # Alignment
    typer.echo("Aligning to reference…")
    rev_iupac = reverse_complement_trace(rev_iupac_trimmed)
    try:
        from hbb_pipeline.alignment import AlignmentFailedError
        fwd_aligned = align_to_reference(fwd_iupac, ref)
        rev_aligned = align_to_reference(rev_iupac, ref)
    except Exception as exc:
        typer.echo(f"ERROR: Alignment failed: {exc}", err=True)
        raise typer.Exit(1)

    consensus, cons_quals = build_consensus(fwd_aligned, rev_aligned, ref)

    # Variant calling
    typer.echo("Calling variants…")
    alignment_variants = call_variants_from_alignment(
        consensus, cons_quals, ref, translator,
        fwd_iupac, rev_iupac, fwd_aligned, rev_aligned,
    )

    tracy_used = False
    tracy_variants: list = []
    if is_tracy_available():
        typer.echo(f"  Running tracy ({tracy_version()}) on fwd + rev…")
        with tempfile.TemporaryDirectory() as tmpdir:
            tracy_variants = run_tracy_on_pair(
                fwd, rev,
                reference,
                Path(tmpdir) / "tracy",
                translator, ref,
            )
            tracy_used = True
        typer.echo(f"  Tracy called {len(tracy_variants)} variant(s)")
    else:
        typer.echo("  tracy not available — using pure-Python HET detection")

    merged = merge_variant_calls(alignment_variants, tracy_variants)
    typer.echo(f"  {len(merged)} variant(s) called")

    # Build QC
    qc = {
        "mean_phred_fwd": round(sum(fwd_raw.phred_scores) / max(len(fwd_raw.phred_scores), 1), 1),
        "mean_phred_rev": round(sum(rev_raw.phred_scores) / max(len(rev_raw.phred_scores), 1), 1),
        "usable_length_fwd": fe - fs,
        "usable_length_rev": re - rs,
        "tracy_available": is_tracy_available(),
        "tracy_version": tracy_version() or "N/A",
    }

    sample_id = fwd.stem
    report = generate_report(merged, qc, sample_id, consensus,
                             tracy_used=tracy_used)

    # Output
    md = render_markdown_report(report)

    if out:
        out.write_text(md, encoding="utf-8")
        typer.echo(f"  Markdown report: {out}")
    else:
        typer.echo("\n" + md)

    if json_out:
        json_out.write_text(
            report.model_dump_json(indent=2),
            encoding="utf-8",
        )
        typer.echo(f"  JSON report: {json_out}")

    # Print summary of known variants
    known = [v for v in merged if v.known_variant_name]
    if known:
        typer.echo("\n=== Clinically Significant Findings ===")
        for v in known:
            typer.echo(
                f"  {v.hgvs_c}  {v.known_variant_name}  [{v.zygosity.value}]  "
                f"p.: {v.hgvs_p_hgvs or '—'} / {v.hgvs_p_legacy or '—'}"
            )

    flagged = [v for v in merged if v.requires_manual_review]
    if flagged:
        typer.echo(f"\n⚠  {len(flagged)} variant(s) require manual review.")


@app.command(name="validate-reference")
def validate_reference(
    reference: Path = typer.Argument(..., help="Path to the case-annotated reference FASTA"),
) -> None:
    """Validate the structure of a case-annotated HBB reference FASTA."""
    from hbb_pipeline.reference import HBBReference, ReferenceValidationError

    typer.echo(f"Validating: {reference}")
    try:
        ref = HBBReference(reference)
    except ReferenceValidationError as exc:
        typer.echo(f"FAIL: {exc}", err=True)
        raise typer.Exit(1)

    typer.echo(f"✓ PASS: {ref.length} bp")
    typer.echo(f"  5' flank:  {ref.utr5[1] - ref.utr5[0]} bp")
    typer.echo(f"  Exon 1:    {ref.exons[0][1] - ref.exons[0][0]} bp  (starts ATG: {ref.upper_seq[978:981]})")
    typer.echo(f"  Intron 1:  {ref.introns[0][1] - ref.introns[0][0]} bp")
    typer.echo(f"  Exon 2:    {ref.exons[1][1] - ref.exons[1][0]} bp")
    typer.echo(f"  Intron 2:  {ref.introns[1][1] - ref.introns[1][0]} bp")
    typer.echo(f"  Exon 3:    {ref.exons[2][1] - ref.exons[2][0]} bp  (ends TAA: {ref.upper_seq[2149:2152]})")
    typer.echo(f"  3' flank:  {ref.utr3[1] - ref.utr3[0]} bp")
    typer.echo(f"  CDS:       {len(ref.cds)} bp → 147 aa β-globin")


if __name__ == "__main__":
    app()
