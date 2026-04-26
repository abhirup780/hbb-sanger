"""Loads and validates the case-annotated HBB reference FASTA.

Convention: UPPERCASE = CDS (exons), lowercase = non-coding (flanks/introns).
Raises ReferenceValidationError on any structural mismatch — this is a safety gate.
"""

from __future__ import annotations

import logging
from pathlib import Path

from Bio import SeqIO

logger = logging.getLogger(__name__)

# Expected structure (0-based half-open intervals)
_EXPECTED = {
    "total_length": 2750,
    "utr5": (0, 978),
    "exon1": (978, 1070),   # 92 bp, starts ATG
    "intron1": (1070, 1200),
    "exon2": (1200, 1423),  # 223 bp
    "intron2": (1423, 2023),
    "exon3": (2023, 2152),  # 129 bp, ends TAA
    "utr3": (2152, 2750),
    "cds_length": 444,
}

CANONICAL_BETA_GLOBIN_PROTEIN = (
    "MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVLGAFSDGLAHLDNLKGTFAT"
    "LSELHCDKLHVDPENFRLLGNVLVCVLAHHFGKEFTPPVQAAYQKVVAGVANALAHKYH"
)


class ReferenceValidationError(Exception):
    """Raised when the loaded FASTA does not match expected HBB structure."""


class HBBReference:
    """Parsed, validated case-annotated HBB reference sequence.

    Example:
        ref = HBBReference(Path("reference/HBB_reference.fasta"))
        print(ref.length)         # 2750
        print(ref.region_of(978)) # 'exon1'
    """

    def __init__(self, fasta_path: Path) -> None:
        self._load_and_validate(fasta_path)

    def _load_and_validate(self, fasta_path: Path) -> None:
        records = list(SeqIO.parse(fasta_path, "fasta"))
        if len(records) != 1:
            raise ReferenceValidationError(
                f"Expected exactly 1 FASTA record, got {len(records)}"
            )

        raw = str(records[0].seq)
        if len(raw) != _EXPECTED["total_length"]:
            raise ReferenceValidationError(
                f"Reference length {len(raw)} != expected {_EXPECTED['total_length']}"
            )

        self.seq: str = raw  # preserves case
        self.upper_seq: str = raw.upper()
        self.length: int = len(raw)

        # Validate segment case patterns
        self._validate_segment("5' flank", *_EXPECTED["utr5"], expect_upper=False)
        self._validate_segment("exon1", *_EXPECTED["exon1"], expect_upper=True)
        self._validate_segment("intron1", *_EXPECTED["intron1"], expect_upper=False)
        self._validate_segment("exon2", *_EXPECTED["exon2"], expect_upper=True)
        self._validate_segment("intron2", *_EXPECTED["intron2"], expect_upper=False)
        self._validate_segment("exon3", *_EXPECTED["exon3"], expect_upper=True)
        self._validate_segment("3' flank", *_EXPECTED["utr3"], expect_upper=False)

        # Validate ATG start and TAA stop
        if self.upper_seq[978:981] != "ATG":
            raise ReferenceValidationError(
                f"Exon 1 does not start with ATG: got {self.upper_seq[978:981]!r}"
            )
        if self.upper_seq[2149:2152] != "TAA":
            raise ReferenceValidationError(
                f"Exon 3 does not end with TAA: got {self.upper_seq[2149:2152]!r}"
            )

        # Build structural attributes
        self.exons: list[tuple[int, int]] = [
            _EXPECTED["exon1"], _EXPECTED["exon2"], _EXPECTED["exon3"]
        ]
        self.introns: list[tuple[int, int]] = [
            _EXPECTED["intron1"], _EXPECTED["intron2"]
        ]
        self.utr5: tuple[int, int] = _EXPECTED["utr5"]
        self.utr3: tuple[int, int] = _EXPECTED["utr3"]

        # Concatenated CDS (exons only, uppercase)
        self.cds: str = "".join(
            self.upper_seq[s:e] for s, e in self.exons
        )
        if len(self.cds) != _EXPECTED["cds_length"]:
            raise ReferenceValidationError(
                f"CDS length {len(self.cds)} != expected {_EXPECTED['cds_length']}"
            )

        # Validate translation
        self._validate_protein()

        self.version: str = records[0].id
        logger.info("HBB reference loaded and validated: %s (%d bp)", self.version, self.length)

    def _validate_segment(
        self, name: str, start: int, end: int, expect_upper: bool
    ) -> None:
        seg = self.seq[start:end]
        if expect_upper:
            if not seg.isupper():
                bad = [(i + start, c) for i, c in enumerate(seg) if not c.isupper()]
                raise ReferenceValidationError(
                    f"{name} [{start},{end}) should be all UPPERCASE; "
                    f"first offenders: {bad[:5]}"
                )
        else:
            if not seg.islower():
                bad = [(i + start, c) for i, c in enumerate(seg) if not c.islower()]
                raise ReferenceValidationError(
                    f"{name} [{start},{end}) should be all lowercase; "
                    f"first offenders: {bad[:5]}"
                )

    def _validate_protein(self) -> None:
        from Bio.Seq import Seq

        protein = str(Seq(self.cds).translate())
        # protein ends with '*' for stop codon
        if not protein.endswith("*"):
            raise ReferenceValidationError(
                f"CDS translation does not end with stop codon: {protein[-5:]!r}"
            )
        translated = protein.rstrip("*")
        if translated != CANONICAL_BETA_GLOBIN_PROTEIN:
            # Show first difference
            for i, (a, b) in enumerate(zip(translated, CANONICAL_BETA_GLOBIN_PROTEIN)):
                if a != b:
                    raise ReferenceValidationError(
                        f"CDS translates to non-canonical β-globin at aa {i+1}: "
                        f"got {a!r}, expected {b!r}"
                    )
            raise ReferenceValidationError(
                f"CDS protein length mismatch: got {len(translated)}, "
                f"expected {len(CANONICAL_BETA_GLOBIN_PROTEIN)}"
            )

    def region_of(self, genomic_pos: int) -> str:
        """Return the region name for a 0-based genomic position.

        Returns one of: '5UTR', 'exon1', 'intron1', 'exon2', 'intron2', 'exon3', '3UTR'.

        Example:
            ref.region_of(978)   # 'exon1'
            ref.region_of(1070)  # 'intron1'
        """
        if not (0 <= genomic_pos < self.length):
            raise ValueError(f"Position {genomic_pos} out of range [0, {self.length})")

        s5, e5 = self.utr5
        if s5 <= genomic_pos < e5:
            return "5UTR"

        s1, e1 = self.exons[0]
        if s1 <= genomic_pos < e1:
            return "exon1"

        si1, ei1 = self.introns[0]
        if si1 <= genomic_pos < ei1:
            return "intron1"

        s2, e2 = self.exons[1]
        if s2 <= genomic_pos < e2:
            return "exon2"

        si2, ei2 = self.introns[1]
        if si2 <= genomic_pos < ei2:
            return "intron2"

        s3, e3 = self.exons[2]
        if s3 <= genomic_pos < e3:
            return "exon3"

        s3u, e3u = self.utr3
        if s3u <= genomic_pos < e3u:
            return "3UTR"

        raise ValueError(f"Position {genomic_pos} falls in no region (bug in reference structure)")
