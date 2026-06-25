"""Convert ClinVar JSON into structure-informed-rvas input TSV.

Output columns include:
- Variant ID: chr:pos:ref:alt
- ac_case: placeholder for GoF variants
- ac_control: placeholder for LoF variants
- ClinVar metadata columns for later filtering/annotation
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


NC_TO_CHR = {
    "NC_000001.11": "chr1",
    "NC_000002.12": "chr2",
    "NC_000003.12": "chr3",
    "NC_000004.12": "chr4",
    "NC_000005.10": "chr5",
    "NC_000006.12": "chr6",
    "NC_000007.14": "chr7",
    "NC_000008.11": "chr8",
    "NC_000009.12": "chr9",
    "NC_000010.11": "chr10",
    "NC_000011.10": "chr11",
    "NC_000012.12": "chr12",
    "NC_000013.11": "chr13",
    "NC_000014.9": "chr14",
    "NC_000015.10": "chr15",
    "NC_000016.10": "chr16",
    "NC_000017.11": "chr17",
    "NC_000018.10": "chr18",
    "NC_000019.10": "chr19",
    "NC_000020.11": "chr20",
    "NC_000021.9": "chr21",
    "NC_000022.11": "chr22",
    "NC_000023.11": "chrX",
    "NC_000024.10": "chrY",
    "NC_012920.1": "chrM",
}


def canonical_spdi_to_chrpos(spdi: str) -> str | None:
    if not spdi:
        return None

    parts = spdi.split(":")
    if len(parts) != 4:
        return None

    nc_accession, pos, ref, alt = parts
    chromosome = NC_TO_CHR.get(nc_accession)

    if chromosome is None:
        return None

    return f"{chromosome}:{pos}:{ref}:{alt}"


def load_clinvar_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def list_to_string(value) -> str:
    if not value:
        return ""
    if isinstance(value, list):
        return "; ".join(str(x) for x in value)
    return str(value)

def clinvar_significance_to_counts(significance: str) -> tuple[int, int]:
    sig = (significance or "").lower()

    # avoid accidentally classifying "likely benign" as pathogenic
    if "benign" in sig:
        return 0, 1

    if "pathogenic" in sig:
        return 1, 0

    return 0, 0

def clinvar_to_rows(clinvar_json: dict) -> list[dict]:
    rows = []

    records = clinvar_json.get("records", [])

    for record in records:
        spdi = record.get("canonical_spdi")
        variant_id = canonical_spdi_to_chrpos(spdi)

        if not variant_id:
            continue

        significance = record.get("significance", "")
        ac_case, ac_control = clinvar_significance_to_counts(significance)      

        rows.append(
            {
                "Variant ID": variant_id,

                # placeholders for now
                # later: GoF = ac_case 1, LoF = ac_control 1
                "ac_case": ac_case,
                "ac_control": ac_control,

                # useful metadata to keep
                "variation_id": record.get("variation_id", ""),
                "accession": record.get("accession", ""),
                "gene": record.get("gene", ""),
                "protein_change": record.get("protein_change", ""),
                "residue": record.get("residue", ""),
                "amino_acid": record.get("amino_acid", ""),
                "alternate_amino_acid": record.get("alternate_amino_acid", ""),
                "significance": record.get("significance", ""),
                "review_status": record.get("review_status", ""),
                "last_evaluated": record.get("last_evaluated", ""),
                "variant_type": record.get("variant_type", ""),
                "molecular_consequences": list_to_string(
                    record.get("molecular_consequences", [])
                ),
                "phenotype_names": list_to_string(record.get("phenotype_names", [])),
                "title": record.get("title", ""),
                "canonical_spdi": spdi,
            }
        )

    return rows


def write_tsv(rows: list[dict], output_file: Path, unique: bool = False) -> None:
    output_file.parent.mkdir(parents=True, exist_ok=True)

    if unique:
        seen = set()
        unique_rows = []

        for row in rows:
            variant_id = row["Variant ID"]
            if variant_id not in seen:
                unique_rows.append(row)
                seen.add(variant_id)

        rows = unique_rows

    if not rows:
        raise ValueError("No valid variants found.")

    fieldnames = list(rows[0].keys())

    with output_file.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Convert ClinVar JSON into structure-informed-rvas input TSV."
    )
    parser.add_argument("input_json", help="Path to ClinVar JSON file")
    parser.add_argument("output_tsv", help="Path to output TSV file")
    parser.add_argument(
        "--unique",
        action="store_true",
        help="Keep only one row per Variant ID",
    )

    args = parser.parse_args()

    clinvar_data = load_clinvar_json(Path(args.input_json))
    rows = clinvar_to_rows(clinvar_data)

    write_tsv(rows, Path(args.output_tsv), unique=args.unique)

    print(f"Saved {len(rows)} rows to {args.output_tsv}")


if __name__ == "__main__":
    main()