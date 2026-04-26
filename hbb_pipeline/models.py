"""Pydantic v2 data models for the HBB Sanger sequencing pipeline."""

from __future__ import annotations

from datetime import datetime
from enum import Enum
from typing import Literal

from pydantic import BaseModel


class TraceData(BaseModel):
    """Parsed ABI trace data for one sequencing read."""

    sequence: str
    phred_scores: list[int]
    channels: dict[str, list[int]]  # keys: "A", "C", "G", "T"
    peak_positions: list[int]  # one entry per base in `sequence`
    sample_name: str
    is_reverse_complemented: bool = False


class AlignedRead(BaseModel):
    """Result of aligning a TraceData to the HBB reference."""

    read_name: str
    reference_start: int  # 0-based genomic
    reference_end: int  # 0-based half-open
    aligned_seq: str  # with gaps ("-")
    aligned_ref: str  # with gaps ("-")
    per_base_quality: list[int]  # same length as ungapped aligned_seq
    trace_index_map: list[int]  # alignment col → index in original trace


class Zygosity(str, Enum):
    HOM = "homozygous"
    HET = "heterozygous"
    UNKNOWN = "unknown"


class Variant(BaseModel):
    """A single variant call with HGVS annotation and zygosity."""

    ref_pos_genomic: int  # 0-based in reference FASTA
    ref_allele: str
    alt_allele: str
    zygosity: Zygosity
    variant_type: Literal["SNV", "INS", "DEL"]
    hgvs_c: str  # e.g. "c.20A>T"
    hgvs_p_hgvs: str | None  # HGVS protein (Met=1): "p.Glu7Val"
    hgvs_p_legacy: str | None  # Legacy mature-protein (Met removed): "p.Glu6Val"
    region: Literal["5UTR", "exon1", "intron1", "exon2", "intron2", "exon3", "3UTR"]
    phred_support: int
    secondary_peak_ratio_fwd: float | None
    secondary_peak_ratio_rev: float | None
    called_by: list[Literal["alignment", "tracy"]]
    known_variant_name: str | None = None  # e.g. "HbS", "IVS1-110", "Cd39"
    requires_manual_review: bool = False
    review_reason: str | None = None


class ClinicalReport(BaseModel):
    """Top-level output of the pipeline for one sample."""

    sample_id: str
    analysis_timestamp: datetime
    variants: list[Variant]
    consensus_seq: str
    qc_metrics: dict
    reference_version: str
    tracy_used: bool
