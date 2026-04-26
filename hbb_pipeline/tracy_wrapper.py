"""Wrapper for the tracy external binary (conda install -c bioconda tracy).

Tracy implements TIDE-like decomposition for Sanger traces, making it the
gold-standard for heterozygous SNV, DEL, and DUP detection.

If tracy is not available at runtime all functions return gracefully and the
pipeline falls back to alignment-only calling.  A warning is added to QC.

Install:  conda install -c bioconda tracy
On Windows (WSL):  wsl -d Ubuntu sudo install -m755 /tmp/tracy /usr/local/bin/tracy
  (binary from https://github.com/gear-genomics/tracy/releases)
"""

from __future__ import annotations

import logging
import shutil
import subprocess
import sys
from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from hbb_pipeline.coordinates import CoordinateTranslator
    from hbb_pipeline.models import Variant
    from hbb_pipeline.reference import HBBReference

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Windows / WSL path helpers
# ---------------------------------------------------------------------------

_IS_WINDOWS = sys.platform == "win32"


def _wsl_path(p: Path) -> str:
    """Convert a Windows Path to a WSL-accessible path string.

    e.g. C:\\Users\\foo\\bar.ab1  →  /mnt/c/Users/foo/bar.ab1
    Falls back to str(p) on non-Windows.
    """
    if not _IS_WINDOWS:
        return str(p)
    try:
        result = subprocess.run(
            ["wsl", "-d", "Ubuntu", "wslpath", str(p).replace("\\", "/")],
            capture_output=True, text=True, timeout=10,
        )
        if result.returncode == 0:
            return result.stdout.strip()
    except Exception:
        pass
    # Fallback: manual conversion C:\foo\bar → /mnt/c/foo/bar
    s = str(p).replace("\\", "/")
    if len(s) >= 2 and s[1] == ":":
        s = "/mnt/" + s[0].lower() + s[2:]
    return s


def _tracy_cmd() -> list[str]:
    """Return the base tracy command appropriate for this platform."""
    if not _IS_WINDOWS:
        return ["tracy"]
    # On Windows: native tracy.bat wrapper calls wsl internally, but paths
    # passed via that wrapper are raw Windows paths that the Linux binary
    # won't parse.  Call WSL directly so we can convert paths ourselves.
    return ["wsl", "-d", "Ubuntu", "tracy"]


# ---------------------------------------------------------------------------
# Availability helpers
# ---------------------------------------------------------------------------

def is_tracy_available() -> bool:
    """Return True if tracy is callable (native or via WSL on Windows)."""
    if shutil.which("tracy") is not None:
        return True
    if _IS_WINDOWS:
        # Check WSL Ubuntu has tracy
        try:
            r = subprocess.run(
                ["wsl", "-d", "Ubuntu", "which", "tracy"],
                capture_output=True, text=True, timeout=10,
            )
            return r.returncode == 0
        except Exception:
            return False
    return False


def tracy_version() -> str | None:
    """Return tracy version string, or None if not available."""
    if not is_tracy_available():
        return None
    try:
        result = subprocess.run(
            [*_tracy_cmd(), "--version"],
            capture_output=True, text=True, timeout=10,
        )
        return (result.stdout or result.stderr).strip().splitlines()[0]
    except Exception:
        return None


# ---------------------------------------------------------------------------
# Core subprocess runner
# ---------------------------------------------------------------------------

def _wsl_stage(src: Path, wsl_dir: str, dest_name: str) -> str:
    """Copy a Windows file into a WSL temp dir, return the WSL path.

    Tracy v0.8.9 rejects paths that contain spaces (boost::program_options
    invalid_option_value).  We stage all inputs into a WSL-native /tmp dir
    that is guaranteed to be space-free before calling the binary.
    """
    src_wsl = _wsl_path(src)
    dest = f"{wsl_dir}/{dest_name}"
    subprocess.run(
        ["wsl", "-d", "Ubuntu", "cp", src_wsl, dest],
        capture_output=True, timeout=30,
    )
    return dest


def _run_decompose(
    abi_path: Path,
    reference_fasta: Path,
    out_dir: Path,
) -> Path | None:
    """Run `tracy decompose` on one ABI file. Returns VCF path or None."""
    out_dir.mkdir(parents=True, exist_ok=True)
    out_prefix = out_dir / abi_path.stem

    if _IS_WINDOWS:
        # Tracy v0.8.9 rejects paths with spaces (boost program_options bug).
        # Stage reference and ABI into a WSL-native /tmp dir, run tracy there,
        # then cat the VCF back out to the Windows output directory.
        stage = subprocess.run(
            ["wsl", "-d", "Ubuntu", "mktemp", "-d"],
            capture_output=True, text=True, timeout=10,
        ).stdout.strip()
        if not stage:
            logger.error("Failed to create WSL staging dir")
            return None

        try:
            ref_wsl  = _wsl_stage(reference_fasta, stage, "ref.fasta")
            abi_wsl  = _wsl_stage(abi_path, stage, abi_path.name)
            out_wsl  = f"{stage}/out"

            cmd = [
                "wsl", "-d", "Ubuntu", "tracy", "decompose",
                "-v",          # --callVariants: required to produce VCF output
                "-t", "3",     # quality-based trimming (removes noisy trace ends)
                "-r", ref_wsl,
                "-o", out_wsl,
                abi_wsl,
            ]
            logger.info("Running (WSL staged): %s", " ".join(cmd))

            try:
                result = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
            except subprocess.TimeoutExpired:
                logger.error("tracy timed out for %s", abi_path.name)
                return None
            except OSError as exc:
                logger.error("Failed to launch tracy: %s", exc)
                return None

            if result.returncode != 0:
                logger.error(
                    "tracy decompose failed for %s (exit %d):\n%s",
                    abi_path.name, result.returncode, result.stderr,
                )
                return None

            # Tracy v0.8.9 outputs .bcf (binary); convert to VCF text via bcftools.
            # Also check plain .vcf in case a future version changes this.
            vcf_text: str | None = None
            wsl_out_vcf = out_wsl + ".vcf"

            # Try BCF → VCF conversion first
            bcf_path = out_wsl + ".bcf"
            bcftools = subprocess.run(
                ["wsl", "-d", "Ubuntu", "bcftools", "view", bcf_path],
                capture_output=True, text=True, timeout=30,
            )
            if bcftools.returncode == 0 and bcftools.stdout.strip():
                vcf_text = bcftools.stdout
                logger.info("tracy BCF converted to VCF via bcftools")
            else:
                # Fallback: try plain .vcf
                for suffix in (".vcf", "_decompose.vcf"):
                    cat = subprocess.run(
                        ["wsl", "-d", "Ubuntu", "cat", out_wsl + suffix],
                        capture_output=True, text=True, timeout=10,
                    )
                    if cat.returncode == 0 and cat.stdout.strip():
                        vcf_text = cat.stdout
                        break

            if not vcf_text:
                logger.error("tracy ran but produced no VCF/BCF at %s/out*", stage)
                return None

            win_vcf = Path(str(out_prefix) + ".vcf")
            win_vcf.write_text(vcf_text, encoding="utf-8")
            logger.info("tracy VCF (staged): %s", win_vcf)
            return win_vcf
        finally:
            # Clean up WSL staging dir
            subprocess.run(
                ["wsl", "-d", "Ubuntu", "rm", "-rf", stage],
                capture_output=True, timeout=10,
            )

    # --- Non-Windows path: direct call ---
    cmd = [
        "tracy", "decompose",
        "-v",          # --callVariants: required to produce VCF output
        "-t", "3",     # quality-based trimming (removes noisy trace ends)
        "-r", str(reference_fasta),
        "-o", str(out_prefix),
        str(abi_path),
    ]
    logger.info("Running: %s", " ".join(cmd))

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
    except subprocess.TimeoutExpired:
        logger.error("tracy timed out for %s", abi_path.name)
        return None
    except OSError as exc:
        logger.error("Failed to launch tracy: %s", exc)
        return None

    if result.returncode != 0:
        logger.error(
            "tracy decompose failed for %s (exit %d):\n%s",
            abi_path.name, result.returncode, result.stderr,
        )
        return None

    # Tracy v0.8.9+ outputs .bcf; convert via bcftools if available
    bcf_path = Path(str(out_prefix) + ".bcf")
    if bcf_path.exists() and shutil.which("bcftools"):
        vcf_path = Path(str(out_prefix) + ".vcf")
        conv = subprocess.run(
            ["bcftools", "view", str(bcf_path), "-o", str(vcf_path), "-O", "v"],
            capture_output=True, text=True, timeout=30,
        )
        if conv.returncode == 0 and vcf_path.exists():
            logger.info("tracy BCF → VCF: %s", vcf_path)
            return vcf_path

    for suffix in (".vcf", "_decompose.vcf"):
        vcf_path = Path(str(out_prefix) + suffix)
        if vcf_path.exists():
            logger.info("tracy VCF: %s", vcf_path)
            return vcf_path

    logger.error("tracy ran but produced no VCF at %s*.vcf", out_prefix)
    return None


# ---------------------------------------------------------------------------
# VCF parser (no pysam — plain text)
# ---------------------------------------------------------------------------

def _parse_vcf(
    vcf_path: Path,
    translator: "CoordinateTranslator",
    ref: "HBBReference",
) -> list["Variant"]:
    """Parse a tracy VCF file without pysam.

    Handles SNV, DEL, INS/DUP using standard VCF anchor-base notation.
    Skips ref/ref calls (GT 0/0) and FILTER != PASS lines.
    Annotates known HBB pathogenic variants.
    """
    from hbb_pipeline.known_variants import lookup_variant
    from hbb_pipeline.models import Variant, Zygosity

    try:
        text = vcf_path.read_text(encoding="utf-8", errors="replace")
    except OSError as exc:
        logger.error("Cannot read tracy VCF %s: %s", vcf_path, exc)
        return []

    variants: list[Variant] = []

    for line in text.splitlines():
        if line.startswith("#"):
            continue
        parts = line.split("\t")
        if len(parts) < 8:
            continue

        _chrom, pos_str, _id, vcf_ref, vcf_alt_field, qual_str, filter_, _info = parts[:8]
        fmt_str  = parts[8] if len(parts) > 8 else ""
        samp_str = parts[9] if len(parts) > 9 else ""

        # Skip filtered records
        if filter_ not in ("PASS", ".", ""):
            continue

        # Skip monomorphic / missing alt
        if vcf_alt_field in (".", ""):
            continue

        vcf_alt = vcf_alt_field.split(",")[0]  # take first alt only

        # Skip if ref equals alt (shouldn't happen)
        if vcf_ref.upper() == vcf_alt.upper():
            continue

        # --- VCF POS → 0-based genomic position of the anchor base ---
        try:
            pos_0 = int(pos_str) - 1  # VCF is 1-based
        except ValueError:
            continue

        if pos_0 < 0 or pos_0 >= ref.length:
            continue

        # --- Zygosity from GT field ---
        zygosity = Zygosity.UNKNOWN
        if fmt_str and samp_str:
            fmt_keys = fmt_str.split(":")
            fmt_vals = samp_str.split(":")
            fmt = dict(zip(fmt_keys, fmt_vals))
            gt_raw = fmt.get("GT", "")
            gt_alleles = [a for a in gt_raw.replace("|", "/").split("/") if a not in (".", "")]
            if len(gt_alleles) == 2:
                if all(a == "0" for a in gt_alleles):
                    continue  # ref/ref — skip entirely
                if gt_alleles[0] == gt_alleles[1]:
                    zygosity = Zygosity.HOM
                else:
                    zygosity = Zygosity.HET

        # --- QUAL → phred ---
        try:
            phred = max(0, int(float(qual_str))) if qual_str not in (".", "") else 0
        except ValueError:
            phred = 0

        if phred < 20:
            logger.debug("Tracy variant at %s skipped: phred %d < 20", pos_str, phred)
            continue

        # --- Determine variant type and strip VCF anchor base(s) ---
        # Standard VCF: anchor base(s) shared at the left are included in REF and ALT.
        # DEL example:  REF=ACTTT  ALT=A      → deleted sequence = CTTT
        # INS example:  REF=G      ALT=GG     → inserted sequence = G
        # SNV example:  REF=A      ALT=T      → simple substitution

        if len(vcf_ref) == 1 and len(vcf_alt) == 1:
            vtype = "SNV"
            actual_ref = vcf_ref
            actual_alt = vcf_alt
            genomic_pos = pos_0

        elif len(vcf_ref) > len(vcf_alt):
            # Deletion — strip shared prefix (anchor)
            vtype = "DEL"
            anchor = 0
            for i in range(min(len(vcf_ref), len(vcf_alt))):
                if vcf_ref[i].upper() == vcf_alt[i].upper():
                    anchor += 1
                else:
                    break
            actual_ref = vcf_ref[anchor:]
            actual_alt = ""
            genomic_pos = pos_0 + anchor

        else:
            # Insertion / duplication — strip shared prefix
            vtype = "INS"
            anchor = 0
            for i in range(min(len(vcf_ref), len(vcf_alt))):
                if vcf_ref[i].upper() == vcf_alt[i].upper():
                    anchor += 1
                else:
                    break
            actual_ref = ""
            actual_alt = vcf_alt[anchor:]
            genomic_pos = pos_0 + anchor

        if genomic_pos < 0 or genomic_pos >= ref.length:
            continue

        # --- Build HGVS via translator ---
        try:
            region = ref.region_of(genomic_pos)
            stub = Variant(
                ref_pos_genomic=genomic_pos,
                ref_allele=actual_ref,
                alt_allele=actual_alt,
                zygosity=zygosity,
                variant_type=vtype,  # type: ignore[arg-type]
                hgvs_c="",
                hgvs_p_hgvs=None,
                hgvs_p_legacy=None,
                region=region,  # type: ignore[arg-type]
                phred_support=phred,
                secondary_peak_ratio_fwd=None,
                secondary_peak_ratio_rev=None,
                called_by=["tracy"],
            )
            hgvs_c      = translator.build_hgvs_c(stub)
            hgvs_p_hgvs = translator.build_hgvs_p(stub, "hgvs")
            hgvs_p_legacy = translator.build_hgvs_p(stub, "legacy")
        except Exception as exc:
            logger.warning("HGVS build failed for tracy variant at %d: %s", genomic_pos, exc)
            continue

        # --- Annotate known pathogenic variants ---
        key = hgvs_c[2:] if hgvs_c.startswith("c.") else hgvs_c
        kv = lookup_variant(key)

        variants.append(Variant(
            ref_pos_genomic=genomic_pos,
            ref_allele=actual_ref,
            alt_allele=actual_alt,
            zygosity=zygosity,
            variant_type=vtype,  # type: ignore[arg-type]
            hgvs_c=hgvs_c,
            hgvs_p_hgvs=hgvs_p_hgvs,
            hgvs_p_legacy=hgvs_p_legacy,
            region=region,  # type: ignore[arg-type]
            phred_support=phred,
            secondary_peak_ratio_fwd=None,
            secondary_peak_ratio_rev=None,
            called_by=["tracy"],
            known_variant_name=kv.name if kv else None,
        ))

    logger.info("Parsed %d variants from %s", len(variants), vcf_path.name)
    return variants


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _indel_key(v: "Variant") -> tuple:
    """Merge key for concordance matching between fwd and rev tracy calls.

    SNVs: exact (pos, ref, alt).
    DEL/INS: (pos, type, ref_len, alt_len) — omits the allele string so that
    the same net indel matches even if tracy chose a slightly different sequence
    representation between strands (same position, same net change = same event).
    Anchor-position drift of >0 bases will still produce two separate entries
    (the safe outcome: both flagged for review rather than wrongly merged).
    """
    if v.variant_type == "SNV":
        return (v.ref_pos_genomic, v.ref_allele, v.alt_allele)
    return (v.ref_pos_genomic, v.variant_type, len(v.ref_allele), len(v.alt_allele))


# ---------------------------------------------------------------------------
# Paired-strand analysis — main public entry point
# ---------------------------------------------------------------------------

def run_tracy_on_pair(
    fwd_path: Path | None,
    rev_path: Path | None,
    reference_fasta: Path,
    out_dir: Path,
    translator: "CoordinateTranslator",
    ref: "HBBReference",
) -> list["Variant"]:
    """Run tracy decompose on forward and/or reverse traces, then merge.

    Both strands are decomposed independently.  Results are merged by
    (genomic_pos, ref_allele, alt_allele):

    - Called by BOTH strands → high-confidence, no review flag.
    - Called by one strand only → kept, flagged for manual confirmation.

    This is especially important for DEL/DUP: a deletion at the end of the
    forward read may only be seen by the reverse (or vice versa).

    Returns:
        Merged list of Variant objects sorted by genomic position.
        Empty list if tracy is unavailable or both runs fail.
    """
    if not is_tracy_available():
        logger.warning(
            "tracy not found on PATH — DEL/DUP detection unavailable. "
            "Install via: conda install -c bioconda tracy"
        )
        return []

    from hbb_pipeline.models import Zygosity

    fwd_variants: list[Variant] = []
    rev_variants: list[Variant] = []

    if fwd_path is not None and fwd_path.exists():
        vcf = _run_decompose(fwd_path, reference_fasta, out_dir / "fwd")
        if vcf:
            fwd_variants = _parse_vcf(vcf, translator, ref)
            logger.info("Tracy fwd: %d variant(s)", len(fwd_variants))

    if rev_path is not None and rev_path.exists():
        vcf = _run_decompose(rev_path, reference_fasta, out_dir / "rev")
        if vcf:
            rev_variants = _parse_vcf(vcf, translator, ref)
            logger.info("Tracy rev: %d variant(s)", len(rev_variants))

    if not fwd_variants and not rev_variants:
        return []

    fwd_map = {_indel_key(v): v for v in fwd_variants}
    rev_map = {_indel_key(v): v for v in rev_variants}
    all_keys = set(fwd_map) | set(rev_map)

    merged: list[Variant] = []
    for key in sorted(all_keys):
        fv = fwd_map.get(key)
        rv = rev_map.get(key)

        if fv is not None and rv is not None:
            # Both strands agree — use higher phred, reconcile zygosity
            zyg = (
                fv.zygosity if fv.zygosity == rv.zygosity
                else Zygosity.HET if Zygosity.HET in (fv.zygosity, rv.zygosity)
                else Zygosity.UNKNOWN
            )
            merged.append(fv.model_copy(update={
                "zygosity": zyg,
                "phred_support": max(fv.phred_support, rv.phred_support),
                "requires_manual_review": False,
                "review_reason": None,
            }))
            logger.info(
                "Tracy concordant: %s [%s] phred=%d",
                fv.hgvs_c, zyg.value, max(fv.phred_support, rv.phred_support),
            )
        elif fv is not None:
            # Single-strand only: only keep DEL/INS.
            # SNVs from one strand are noise-prone; the alignment pipeline already
            # handles SNVs reliably.  Known pathogenic variants are always kept.
            if fv.variant_type == "SNV" and not fv.known_variant_name:
                logger.warning(
                    "Tracy fwd-only SNV dropped (single-strand, not in known variants): %s",
                    fv.hgvs_c,
                )
                continue
            merged.append(fv.model_copy(update={
                "requires_manual_review": True,
                "review_reason": "Tracy: forward strand only — reverse trace did not call this variant",
            }))
            logger.info("Tracy fwd-only: %s", fv.hgvs_c)
        else:
            assert rv is not None
            # Same filter: single-strand SNVs dropped unless known pathogenic
            if rv.variant_type == "SNV" and not rv.known_variant_name:
                logger.warning(
                    "Tracy rev-only SNV dropped (single-strand, not in known variants): %s",
                    rv.hgvs_c,
                )
                continue
            merged.append(rv.model_copy(update={
                "requires_manual_review": True,
                "review_reason": "Tracy: reverse strand only — forward trace did not call this variant",
            }))
            logger.info("Tracy rev-only: %s", rv.hgvs_c)

    return merged
