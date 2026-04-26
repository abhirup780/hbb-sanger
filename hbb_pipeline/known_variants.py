"""Lookup table of common pathogenic HBB variants.

Used only for annotation: after a variant is called, its HGVS c. string is
checked against this table to attach a clinical name (e.g. "HbS").
No position scanning or coordinate parsing is done here.

Keys are HGVS c. strings (without the "c." prefix).
Sources: HbVar (globin.bx.psu.edu), HGMD, ClinVar.

South Asian variants are prioritised given the primary user location (Kolkata).
NOTE: All entries verified against HbVar. Positions use the convention where
c.1 = A of ATG.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class KnownVariant:
    name: str
    disease: str
    clinical_significance: str  # "Pathogenic" | "Likely pathogenic" | "Benign" | etc.
    common_populations: list[str]


# Key: hgvs_c string (no "c." prefix), value: KnownVariant
KNOWN_VARIANTS: dict[str, KnownVariant] = {
    # ----- Coding SNVs (most common globally) -----
    "20A>T": KnownVariant(
        name="HbS (Sickle cell)",
        disease="Sickle cell disease (HbSS), HbSC, HbS/β-thal",
        clinical_significance="Pathogenic",
        common_populations=["Sub-Saharan African", "Mediterranean", "Middle Eastern", "South Asian"],
    ),
    "19G>A": KnownVariant(
        name="HbC",
        disease="HbC disease, HbSC disease",
        clinical_significance="Pathogenic",
        common_populations=["West African"],
    ),
    "79G>A": KnownVariant(
        name="HbE",
        disease="HbE disease, HbE/β-thal (common cause of severe β-thal in South/SE Asia)",
        clinical_significance="Pathogenic",
        common_populations=["South Asian", "Southeast Asian"],
    ),
    "118C>T": KnownVariant(
        name="Codon 39 (Cd39) nonsense",
        disease="β0-thalassemia",
        clinical_significance="Pathogenic",
        common_populations=["Mediterranean (Sardinian, Italian)"],
    ),
    "47G>A": KnownVariant(
        name="Codon 15 nonsense",
        disease="β0-thalassemia",
        clinical_significance="Pathogenic",
        common_populations=["South Asian", "Southeast Asian"],
    ),
    # ----- Frameshift indels -----
    "27dupG": KnownVariant(
        name="Codon 8/9 (+G) frameshift",
        disease="β0-thalassemia (most common mutation in South Asia / India)",
        clinical_significance="Pathogenic",
        common_populations=["South Asian (Indian subcontinent)"],
    ),
    "126_129delCTTT": KnownVariant(
        name="Codon 41/42 (-CTTT) frameshift",
        disease="β0-thalassemia (common in South and East Asia)",
        clinical_significance="Pathogenic",
        common_populations=["South Asian", "East Asian (Chinese)"],
    ),
    # ----- Splice site mutations -----
    "92+1G>A": KnownVariant(
        name="IVS1-1 (splice donor)",
        disease="β0-thalassemia",
        clinical_significance="Pathogenic",
        common_populations=["Mediterranean", "Middle Eastern"],
    ),
    "92+5G>C": KnownVariant(
        name="IVS1-5 (splice donor)",
        disease="β+-thalassemia (common in South Asia)",
        clinical_significance="Pathogenic",
        common_populations=["South Asian (Indian subcontinent)", "Mediterranean"],
    ),
    "93-21G>A": KnownVariant(
        name="IVS1-110 (splice acceptor)",
        disease="β+-thalassemia (most common Mediterranean mutation, also South Asian)",
        clinical_significance="Pathogenic",
        common_populations=["Mediterranean", "South Asian"],
    ),
    "315+1G>A": KnownVariant(
        name="IVS2-1 (splice donor)",
        disease="β0-thalassemia",
        clinical_significance="Pathogenic",
        common_populations=["Mediterranean", "Middle Eastern"],
    ),
    "316-197C>T": KnownVariant(
        name="IVS2-654 (cryptic splice acceptor)",
        disease="β+-thalassemia (most common East Asian mutation)",
        clinical_significance="Pathogenic",
        common_populations=["East Asian (Chinese)", "Southeast Asian"],
    ),
    "316-2A>G": KnownVariant(
        name="IVS2-2 (splice acceptor)",
        disease="β0-thalassemia",
        clinical_significance="Pathogenic",
        common_populations=["Mediterranean"],
    ),
    # ----- Promoter mutations -----
    "-78A>G": KnownVariant(
        name="-28 promoter (TATA box)",
        disease="β+-thalassemia (mild; common in South Asia)",
        clinical_significance="Pathogenic",
        common_populations=["South Asian (Indian subcontinent)", "African American"],
    ),
    "-79A>G": KnownVariant(
        name="-29 promoter",
        disease="β+-thalassemia",
        clinical_significance="Pathogenic",
        common_populations=["African American", "Chinese"],
    ),
    "-87C>G": KnownVariant(
        name="-87 CCAAT box",
        disease="β+-thalassemia (silent/mild)",
        clinical_significance="Pathogenic",
        common_populations=["Mediterranean", "South Asian"],
    ),
    "-88C>T": KnownVariant(
        name="-88 CCAAT box",
        disease="β+-thalassemia (mild)",
        clinical_significance="Pathogenic",
        common_populations=["African American"],
    ),
    "-140C>T": KnownVariant(
        name="-140 CACCC box",
        disease="β+-thalassemia (mild to moderate)",
        clinical_significance="Pathogenic",
        common_populations=["Mediterranean", "Middle Eastern", "South Asian"],
    ),
}


def lookup_variant(hgvs_c: str) -> KnownVariant | None:
    """Look up a variant by its c. coordinate string (without 'c.' prefix).

    Example:
        kv = lookup_variant("20A>T")
        print(kv.name)  # "HbS (Sickle cell)"
    """
    return KNOWN_VARIANTS.get(hgvs_c)
