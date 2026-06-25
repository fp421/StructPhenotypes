"""Utilities to extract and summarise phenotype labels from ClinVar JSON files.

This focuses on broad phenotype names (the `phenotype_names` / `phenotypes.name`
fields in the ClinVar export) and reports unique phenotype counts and variant
lists for each phenotype.
"""

from __future__ import annotations

import json
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, Iterable, List


def load_clinvar_json(path: str | Path) -> Dict:
    path = Path(path)
    with path.open("r", encoding="utf-8") as fh:
        return json.load(fh)


def _phenotype_names_from_record(rec: Dict) -> List[str]:
    """Return list of phenotype name strings for one ClinVar record.

    Prefers the `phenotype_names` list if present; otherwise falls back to
    extracting `phenotypes[*].name`.
    """
    names = []
    if rec.get("phenotype_names"):
        for n in rec.get("phenotype_names", []):
            if n:
                names.append(str(n).strip())
    elif rec.get("phenotypes"):
        for p in rec.get("phenotypes", []):
            n = p.get("name") if isinstance(p, dict) else None
            if n:
                names.append(str(n).strip())
    return names


def phenotype_counts(clinvar: Dict) -> Dict[str, Dict]:
    """Return a mapping of phenotype -> {count, variants}.

    - `count` is the number of ClinVar records that list that phenotype.
    - `variants` is a list of example accessions (up to 20) that have the
      phenotype.
    """
    counts = Counter()
    variants_by_pheno = defaultdict(list)

    for rec in clinvar.get("records", []):
        accession = rec.get("accession") or rec.get("variation_id")
        pheno_names = list(dict.fromkeys(_phenotype_names_from_record(rec)))
        for name in pheno_names:
            counts[name] += 1
            if accession and len(variants_by_pheno[name]) < 20:
                variants_by_pheno[name].append(accession)

    result = {}
    for name, cnt in counts.most_common():
        result[name] = {"count": cnt, "variants": variants_by_pheno.get(name, [])}
    return result


def unique_phenotype_list(clinvar: Dict) -> List[str]:
    """Return sorted list of unique phenotype names in the ClinVar data."""
    return sorted(phenotype_counts(clinvar).keys())


def print_summary(clinvar_path: str | Path, top: int = 30) -> None:
    clinvar = load_clinvar_json(clinvar_path)
    pc = phenotype_counts(clinvar)
    total_variants = len(clinvar.get("records", []))
    total_unique = len(pc)

    print(f"ClinVar file: {clinvar_path}")
    print(f"Total variants: {total_variants}")
    print(f"Unique phenotype labels: {total_unique}\n")

    print(f"Top {top} phenotype labels by variant count:")
    for i, (name, info) in enumerate(list(pc.items())[:top], start=1):
        sample = ", ".join(info["variants"]) if info["variants"] else "-"
        print(f"{i:2d}. {name} — {info['count']} variants (examples: {sample})")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Summarise ClinVar phenotype labels")
    parser.add_argument("clinvar", nargs="?", default="data/clinvar/scn2a_clinvar.json")
    parser.add_argument("--top", type=int, default=30)
    args = parser.parse_args()

    print_summary(args.clinvar, top=args.top)
