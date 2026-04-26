# HBB Sanger Analysis

> Clinical-grade Sanger sequencing analysis for **β-thalassaemia** and **sickle cell disease** — variant calling, HGVS annotation, zygosity, and trace visualisation in a single web app.

[![CI](https://github.com/abhirup780/hbb-sanger/actions/workflows/ci.yml/badge.svg)](https://github.com/abhirup780/hbb-sanger/actions/workflows/ci.yml)
[![Python 3.11+](https://img.shields.io/badge/python-3.11%2B-blue)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Streamlit](https://img.shields.io/badge/Streamlit-Community%20Cloud-FF4B4B?logo=streamlit)](https://streamlit.io/cloud)

---

## Features

| Capability | Detail |
|---|---|
| **SNV detection** | Pre-alignment IUPAC coding — catches HET calls the basecaller silences |
| **Indel detection** | TIDE-based decomposition via [tracy](https://github.com/gear-genomics/tracy) (optional) |
| **HGVS annotation** | Programmatic c. and p. notation — no hardcoded variant lists |
| **Zygosity** | Dual-strand secondary peak ratios (threshold 0.25) |
| **Known variants** | Curated registry of clinically significant HBB variants with automatic lookup |
| **Trace viewer** | Interactive Plotly chromatograms with peak annotations |
| **Coverage map** | Genomic coverage visualised across the HBB locus |
| **QC report** | Phred scores, usable length, tracy availability |
| **Export** | Markdown report + JSON for downstream processing |

---

## Quickstart

### Web app (Streamlit)

```bash
pip install -r requirements.txt
streamlit run app.py
```

Upload a forward + reverse `.ab1` pair → click **Run Analysis**.

### Command-line

```bash
pip install -r requirements.txt
python cli.py run forward.ab1 reverse.ab1 --out report.md
python cli.py run forward.ab1 reverse.ab1 --json result.json
python cli.py validate-reference reference/HBB_reference.fasta
```

---

## Optional: tracy (heterozygous indel detection)

tracy implements TIDE-based decomposition for HET DEL/DUP detection.  
The pipeline falls back to alignment-only calling if tracy is not installed.

**Linux / macOS (conda):**
```bash
conda install -c bioconda tracy
```

**Windows (WSL):**
```bash
# inside WSL Ubuntu
wget https://github.com/gear-genomics/tracy/releases/latest/download/tracy_linux_x86_64bit.tar.gz
tar xf tracy_linux_x86_64bit.tar.gz
sudo install -m755 tracy /usr/local/bin/tracy
```

---

## Project structure

```
hbb-sanger/
├── app.py                  # Streamlit web app entry point
├── cli.py                  # Command-line interface (typer)
├── plots.py                # Plotly chromatogram visualisations
├── requirements.txt
│
├── hbb_pipeline/           # Core analysis engine
│   ├── alignment.py        # Pairwise alignment + IUPAC pre-coding
│   ├── coordinates.py      # Genomic ↔ HGVS coordinate translation
│   ├── heterozygosity.py   # Secondary peak detection and zygosity
│   ├── known_variants.py   # Known pathogenic variant registry
│   ├── models.py           # Pydantic data models
│   ├── parsing.py          # ABI trace parser + Mott quality trimming
│   ├── reference.py        # HBB reference FASTA loader + validator
│   ├── reporting.py        # Markdown/JSON report generation
│   ├── tracy_wrapper.py    # tracy subprocess wrapper + VCF parser
│   └── variants.py         # Variant calling + merge logic
│
├── reference/
│   └── HBB_reference.fasta # Case-annotated 2750 bp HBB reference
│
└── tests/                  # pytest test suite (55 tests)
```

---

## How it works

```
ABI file
   │
   ├─ apply_iupac_symbols()      ← Secondary peak scan, IUPAC coding (pre-trim)
   │    Basecaller Q1 at HET positions → boosted to Q20, Mott trimmer preserves them
   │
   ├─ trim_by_quality()          ← Mott's modified algorithm (Q20)
   │
   ├─ align_to_reference()       ← BioPython PairwiseAligner, local mode
   │
   ├─ build_consensus()          ← Dual-strand merge with quality reconciliation
   │
   ├─ call_variants_from_alignment()  ← SNV + indel calling from consensus
   │    IUPAC codes in consensus → HET SNV calls
   │
   ├─ run_tracy_on_pair()        ← Optional: TIDE decomposition for HET indels
   │    Concordant on both strands → high confidence, no review flag
   │
   └─ generate_report()          ← HGVS annotation, zygosity, known variant lookup
```

---

## Reference

The bundled `reference/HBB_reference.fasta` is a **case-annotated 2750 bp** sequence spanning the complete HBB locus:

| Region | Coordinates |
|---|---|
| 5′ flank | 0–977 |
| Exon 1 | 978–1071 |
| Intron 1 | 1072–1541 |
| Exon 2 | 1542–1694 |
| Intron 2 | 1695–2148 |
| Exon 3 | 2149–2279 |
| 3′ flank | 2280–2749 |

HGVS coordinates are relative to NM_000518.5 (CDS start = genomic 978).

---

## Disclaimer

> **For research use only.** This tool is not validated for clinical diagnostic use. All variant calls, particularly those flagged for manual review, must be confirmed by a certified clinical laboratory before any clinical decision is made.

---

## License

MIT © 2024 Abhirup Sarkar
