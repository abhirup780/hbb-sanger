"""Genomic ↔ HGVS c. coordinate translation and HGVS string builder.

Conventions used throughout:
- Genomic positions are 0-based (Python slice convention).
- c. numbering follows HGVS v21:
    * Exonic: c.1 = A of ATG (genomic 978).  No c.0.
    * 5' of c.1: c.-1, c.-2, … (genomic 977, 976, …).
    * 3' of last CDS base: c.*1, c.*2, … (genomic 2152, 2153, …).
    * Intronic: c.N+M (M bases after exon N end) or c.N-M (M bases before
      exon N+1 start), whichever is closer.  At equal distance, use the
      upstream (donor-side) notation.
- Protein numbering:
    * "hgvs"   – Met = position 1  (HbS = p.Glu7Val)
    * "legacy" – mature β-globin, Met removed = position 1  (HbS = p.Glu6Val)
"""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING, Literal

if TYPE_CHECKING:
    from hbb_pipeline.models import Variant
    from hbb_pipeline.reference import HBBReference

logger = logging.getLogger(__name__)

# Standard genetic code (codon → single-letter AA)
_CODON_TABLE: dict[str, str] = {
    "TTT": "Phe", "TTC": "Phe", "TTA": "Leu", "TTG": "Leu",
    "CTT": "Leu", "CTC": "Leu", "CTA": "Leu", "CTG": "Leu",
    "ATT": "Ile", "ATC": "Ile", "ATA": "Ile", "ATG": "Met",
    "GTT": "Val", "GTC": "Val", "GTA": "Val", "GTG": "Val",
    "TCT": "Ser", "TCC": "Ser", "TCA": "Ser", "TCG": "Ser",
    "CCT": "Pro", "CCC": "Pro", "CCA": "Pro", "CCG": "Pro",
    "ACT": "Thr", "ACC": "Thr", "ACA": "Thr", "ACG": "Thr",
    "GCT": "Ala", "GCC": "Ala", "GCA": "Ala", "GCG": "Ala",
    "TAT": "Tyr", "TAC": "Tyr", "TAA": "Ter", "TAG": "Ter",
    "CAT": "His", "CAC": "His", "CAA": "Gln", "CAG": "Gln",
    "AAT": "Asn", "AAC": "Asn", "AAA": "Lys", "AAG": "Lys",
    "GAT": "Asp", "GAC": "Asp", "GAA": "Glu", "GAG": "Glu",
    "TGT": "Cys", "TGC": "Cys", "TGA": "Ter", "TGG": "Trp",
    "CGT": "Arg", "CGC": "Arg", "CGA": "Arg", "CGG": "Arg",
    "AGT": "Ser", "AGC": "Ser", "AGA": "Arg", "AGG": "Arg",
    "GGT": "Gly", "GGC": "Gly", "GGA": "Gly", "GGG": "Gly",
}


class CoordinateConversionError(Exception):
    """Raised when a coordinate cannot be converted."""


class CoordinateTranslator:
    """Translates between 0-based genomic positions and HGVS c. coordinates.

    Example:
        ref = HBBReference(Path("reference/HBB_reference.fasta"))
        t = CoordinateTranslator(ref)
        t.genomic_to_c(978)   # "1"
        t.genomic_to_c(977)   # "-1"
        t.genomic_to_c(1070)  # "92+1"
        t.genomic_to_c(2152)  # "*1"
    """

    def __init__(self, ref: "HBBReference") -> None:
        self._ref = ref
        # Exon boundaries (genomic, 0-based half-open)
        # exon1: [978, 1070)  → c.1..92
        # exon2: [1200, 1423) → c.93..315
        # exon3: [2023, 2152) → c.316..444
        self._exons = ref.exons      # list of (start, end)
        self._introns = ref.introns  # list of (start, end)

        # First and last c. position of each exon
        c_pos = 1
        self._exon_c_ranges: list[tuple[int, int]] = []
        for s, e in self._exons:
            length = e - s
            self._exon_c_ranges.append((c_pos, c_pos + length - 1))
            c_pos += length

        # genomic start of CDS (A of ATG)
        self._cds_start = self._exons[0][0]   # 978
        # genomic end of CDS (exclusive: last base of stop codon is at 2151, end=2152)
        self._cds_end = self._exons[-1][1]     # 2152

    # ------------------------------------------------------------------
    # Public interface
    # ------------------------------------------------------------------

    def genomic_to_c(self, genomic_pos: int) -> str:
        """Convert a 0-based genomic position to an HGVS c. coordinate string.

        Does NOT include the "c." prefix.

        Example:
            genomic_to_c(978)  == "1"
            genomic_to_c(977)  == "-1"
            genomic_to_c(1070) == "92+1"
            genomic_to_c(2152) == "*1"
        """
        # 5' UTR / promoter → negative numbering, -1 is immediately upstream of c.1
        if genomic_pos < self._cds_start:
            offset = self._cds_start - genomic_pos  # 1-based distance
            return f"-{offset}"

        # 3' UTR → *1, *2, …
        if genomic_pos >= self._cds_end:
            offset = genomic_pos - self._cds_end + 1  # 1-based
            return f"*{offset}"

        # Exonic?
        for (gs, ge), (cs, ce) in zip(self._exons, self._exon_c_ranges):
            if gs <= genomic_pos < ge:
                return str(cs + (genomic_pos - gs))

        # Intronic — find surrounding exons
        for i, (is_, ie) in enumerate(self._introns):
            if is_ <= genomic_pos < ie:
                # Previous exon: exons[i], last c. position = exon_c_ranges[i][1]
                prev_exon_end_genomic = self._exons[i][1]    # exclusive end
                next_exon_start_genomic = self._exons[i + 1][0]

                last_c_prev = self._exon_c_ranges[i][1]
                first_c_next = self._exon_c_ranges[i + 1][0]

                dist_from_donor = genomic_pos - prev_exon_end_genomic + 1  # ≥1
                dist_to_acceptor = next_exon_start_genomic - genomic_pos   # ≥1

                if dist_from_donor <= dist_to_acceptor:
                    return f"{last_c_prev}+{dist_from_donor}"
                else:
                    return f"{first_c_next}-{dist_to_acceptor}"

        raise CoordinateConversionError(
            f"Genomic position {genomic_pos} is within declared CDS range "
            f"[{self._cds_start}, {self._cds_end}) but maps to no exon or intron. "
            "Reference structure may be corrupt."
        )

    def c_to_genomic(self, c_coord: str) -> int:
        """Convert an HGVS c. coordinate string to a 0-based genomic position.

        Handles plain integers ("20"), intronic ("92+1", "93-1"), upstream ("-78"),
        and 3' UTR ("*1").

        Example:
            c_to_genomic("1")    == 978
            c_to_genomic("-1")   == 977
            c_to_genomic("92+1") == 1070
            c_to_genomic("93-1") == 1199
            c_to_genomic("*1")   == 2152
        """
        c = c_coord.strip()

        # 3' UTR: *N
        if c.startswith("*"):
            offset = int(c[1:])
            return self._cds_end + offset - 1

        # 5' UTR: -N
        if c.startswith("-"):
            offset = int(c[1:])
            return self._cds_start - offset

        # Intronic: N+M or N-M
        if "+" in c:
            base, delta = c.split("+")
            base_c = int(base)
            delta = int(delta)
            genomic_base = self._c_exon_to_genomic(base_c)
            return genomic_base + delta

        if "-" in c:
            base, delta = c.split("-")
            base_c = int(base)
            delta = int(delta)
            genomic_base = self._c_exon_to_genomic(base_c)
            return genomic_base - delta

        # Plain exon position
        return self._c_exon_to_genomic(int(c))

    def _c_exon_to_genomic(self, c_pos: int) -> int:
        """Convert a plain exon c. position (integer) to genomic 0-based."""
        for (gs, ge), (cs, ce) in zip(self._exons, self._exon_c_ranges):
            if cs <= c_pos <= ce:
                return gs + (c_pos - cs)
        raise CoordinateConversionError(
            f"c.{c_pos} does not map to any exon (valid ranges: "
            + ", ".join(f"c.{cs}–{ce}" for cs, ce in self._exon_c_ranges)
            + ")"
        )

    # ------------------------------------------------------------------
    # HGVS string builders
    # ------------------------------------------------------------------

    def build_hgvs_c(self, variant: "Variant") -> str:
        """Build the full HGVS c. string for a variant.

        SNV:    c.20A>T
        DEL:    c.27delG  /  c.126_129delCTTT
        INS:    c.27_28insG
        DUP:    c.27dupG  (preferred over ins for single-base duplications)

        Indels are right-aligned (3' rule) per HGVS standard.

        Example:
            build_hgvs_c(snv_at_978_plus_19)  # "c.20A>T"  (HbS)
        """
        pos_c = self.genomic_to_c(variant.ref_pos_genomic)
        vtype = variant.variant_type

        if vtype == "SNV":
            return f"c.{pos_c}{variant.ref_allele}>{variant.alt_allele}"

        if vtype == "DEL":
            # Apply HGVS 3' right-alignment rule before building the notation.
            right_start = self._right_align_del(variant.ref_pos_genomic, variant.ref_allele)
            ref = self._ref.upper_seq[right_start: right_start + len(variant.ref_allele)]
            pos_c = self.genomic_to_c(right_start)
            if len(ref) == 1:
                return f"c.{pos_c}del{ref}"
            end_c = self.genomic_to_c(right_start + len(ref) - 1)
            return f"c.{pos_c}_{end_c}del{ref}"

        if vtype == "INS":
            alt = variant.alt_allele
            # Check for duplication (inserted sequence equals preceding bases)
            end_c = self.genomic_to_c(variant.ref_pos_genomic + 1)
            preceding = self._ref.upper_seq[
                variant.ref_pos_genomic - len(alt) + 1: variant.ref_pos_genomic + 1
            ]
            if preceding.upper() == alt.upper():
                if len(alt) == 1:
                    return f"c.{pos_c}dup{alt}"
                return f"c.{pos_c}_{end_c}dup{alt}"
            return f"c.{pos_c}_{end_c}ins{alt}"

        raise CoordinateConversionError(f"Unknown variant type: {vtype!r}")

    def build_hgvs_p(
        self,
        variant: "Variant",
        numbering: Literal["hgvs", "legacy"],
    ) -> str | None:
        """Build the HGVS p. string for a coding variant.

        numbering='hgvs':   Met = aa 1  → HbS = p.Glu7Val
        numbering='legacy': mature β-globin (Met removed) → HbS = p.Glu6Val

        Returns None for non-coding variants.
        Returns "p.(=)" for synonymous variants.

        Example:
            build_hgvs_p(hbs_variant, "hgvs")   # "p.Glu7Val"
            build_hgvs_p(hbs_variant, "legacy")  # "p.Glu6Val"
        """
        region = variant.region
        if region not in ("exon1", "exon2", "exon3"):
            return None  # non-coding — no protein consequence derivable

        # Offset within CDS (0-based)
        cds_offset = self._genomic_to_cds_offset(variant.ref_pos_genomic)
        if cds_offset is None:
            return None

        codon_index = cds_offset // 3      # 0-based codon
        codon_offset = cds_offset % 3      # position within codon (0, 1, 2)

        # Reference codon + amino acid
        cds = self._ref.cds
        ref_codon = cds[codon_index * 3: codon_index * 3 + 3]
        if len(ref_codon) < 3:
            return None
        ref_aa = _CODON_TABLE.get(ref_codon)
        if ref_aa is None:
            return None

        # Alt codon — apply SNV only (indels handled separately)
        if variant.variant_type == "SNV":
            alt_codon = list(ref_codon)
            alt_codon[codon_offset] = variant.alt_allele
            alt_codon_str = "".join(alt_codon)
            alt_aa = _CODON_TABLE.get(alt_codon_str)
            if alt_aa is None:
                return None
        elif variant.variant_type in ("INS", "DEL"):
            aa_num = self._aa_number(codon_index, numbering)
            return self._frameshift_hgvs_p(variant, cds_offset, codon_index, ref_aa, aa_num)
        else:
            return None

        aa_num = self._aa_number(codon_index, numbering)

        if ref_aa == alt_aa:
            return "p.(=)"

        return f"p.{ref_aa}{aa_num}{alt_aa}"

    # ------------------------------------------------------------------
    # Private helpers
    # ------------------------------------------------------------------

    def _genomic_to_cds_offset(self, genomic_pos: int) -> int | None:
        """Return 0-based offset within the concatenated CDS, or None if non-coding."""
        offset = 0
        for gs, ge in self._exons:
            if gs <= genomic_pos < ge:
                return offset + (genomic_pos - gs)
            offset += ge - gs
        return None

    def _aa_number(self, codon_index: int, numbering: Literal["hgvs", "legacy"]) -> int:
        """Convert 0-based codon index to amino-acid number under the chosen convention."""
        hgvs_num = codon_index + 1  # Met = 1
        if numbering == "hgvs":
            return hgvs_num
        # legacy: mature protein, Met (aa 1) is cleaved → subtract 1
        return hgvs_num - 1

    def _frameshift_hgvs_p(
        self,
        variant: "Variant",
        cds_offset: int,
        codon_index: int,
        ref_aa: str,
        aa_num: int,
    ) -> str:
        """Build HGVS p. for a frameshift by translating the downstream mutant CDS.

        Applies the indel to the CDS, translates from the affected codon, and
        finds the new stop codon to produce p.RefNNewfsTerM notation
        (e.g. p.Glu7ValfsTer10 for HbS-adjacent frameshifts).
        Falls back to p.RefNfsTerM / p.RefNfs if no stop is found in frame.

        Translates downstream from the frameshift until a stop codon is found.
        """
        cds = self._ref.cds
        if variant.variant_type == "DEL":
            del_len = len(variant.ref_allele)
            mut_cds = cds[:cds_offset] + cds[cds_offset + del_len:]
        else:  # INS
            mut_cds = cds[:cds_offset] + variant.alt_allele + cds[cds_offset:]

        # Translate from the start of the first affected codon
        codon_start = codon_index * 3
        mut_from_codon = mut_cds[codon_start:]

        # First amino acid in the new reading frame
        new_first_codon = mut_from_codon[:3] if len(mut_from_codon) >= 3 else ""
        new_first_aa = _CODON_TABLE.get(new_first_codon, "?") if new_first_codon else "?"

        # Walk the new frame to find the new stop codon
        ter_pos: int | None = None
        for i in range(0, len(mut_from_codon) - 2, 3):
            codon = mut_from_codon[i:i + 3]
            if _CODON_TABLE.get(codon) == "Ter":
                ter_pos = i // 3 + 1  # 1-based count from frameshift position
                break

        if ter_pos is not None:
            return f"p.{ref_aa}{aa_num}{new_first_aa}fsTer{ter_pos}"
        return f"p.{ref_aa}{aa_num}{new_first_aa}fs"

    def _right_align_del(self, start: int, deleted: str) -> int:
        """Return the rightmost valid start position for a deletion per the HGVS 3' rule.

        Shifts the deletion window one base to the right as long as the base
        immediately before the window equals the base immediately after it
        (i.e. the deleted sequence can be cyclically rotated without changing
        the resulting sequence).  This normalises variants in repetitive regions
        to the 3'-most representation required by HGVS.

        Example — Cd41/42 -CTTT in a CTTTCTTT context:
            start=1233, deleted="CTTT" → shifts to the rightmost CTTT window.
        """
        ref_seq = self._ref.upper_seq
        n = len(deleted)
        pos = start
        while pos + n < len(ref_seq):
            if ref_seq[pos].upper() == ref_seq[pos + n].upper():
                pos += 1
            else:
                break
        return pos
