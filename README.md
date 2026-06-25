# StructPhenotypes

BMEIS Hackathon 2026 project

Team: Francesca Peccei (lead), Anna Sadilova, George Webber, Konrad Wagstyl

Overall aim is to produce a 3D map of given protein (e.g. SCN2A) and annotate amino acid locations according to phenotypes-of-interest and AlphaMissense prediction scores

## Quick Start: AlphaMissense Predictions

Get pathogenicity scores for all amino acid positions in your protein:

```bash
# One-time setup (download ~1.1GB data file of all alpha missense predictions)
bash setup_alphamissense.sh

# Get mean pathogenicity per position (default)
python structpheno/get_alphamissense.py SCN2A
# Output: data/alphamissense/SCN2A_mean_pathogenicity.csv (2005 positions)

# Or get all possible variants
python structpheno/get_alphamissense.py SCN2A --variants
# Output: data/alphamissense/SCN2A_alphamissense.csv (38,114 variants)
```

**Output columns** (default mode):
- `position` — Amino acid position
- `mean_pathogenicity` — Score 0-1 (0=benign, 1=pathogenic)
- `mean_class` — Classification: benign/ambiguous/pathogenic
- `min_pathogenicity`, `max_pathogenicity` — Range of scores
- `std_pathogenicity` — Variation at this position

Use any human gene name: `BRCA1`, `MECP2`, `TP53`, `HTT`, etc.

## Quick Start: gnomAD Variants

Get observed population variants (allele frequencies) for your protein:

```bash
# Fetch variants via the gnomAD API (no setup required)
python structpheno/get_gnomad.py SCN2A
# Output: data/gnomad/SCN2A/variants.csv (6560 variants)

# Also write a JSON report compatible with visualize.py
python structpheno/get_gnomad.py SCN2A --report
# Output: data/gnomad/SCN2A/report.json
```

**Output columns:**
- `variant_id`, `pos`, `ref`, `alt` — Genomic variant identifiers
- `consequence` — Predicted effect (e.g. missense_variant)
- `hgvsp` — Protein-level change
- `residue`, `ref_aa`, `alt_aa` — Amino acid position and substitution
- `allele_count`, `allele_number`, `allele_freq` — Population frequency

Use any human gene name: `BRCA1`, `MECP2`, `TP53`, `HTT`, etc.
