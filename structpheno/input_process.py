"""Create RVAS-compatible SCN1A input from functional variant data + all missense controls.

This script:
1. Matches SCN1A cDNA variants in:
   data/exp_variance/SCN1A/SCN1A_variants_functional.csv

   to ClinVar JSON records in:
   data/clinvar/SCN1A/clinvar.json

2. Converts matched ClinVar canonical_spdi values to:
   chr:pos:ref:alt

3. Adds those matched functional variants as cases.

4. Adds all possible SCN1A missense variants from:
   sir-reference-data/all_missense_variants_gr38.h5

   as controls.

Output:
   input/SCN1A_cases_vs_all_missense.tsv
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from pathlib import Path
from typing import Dict, Iterable, Optional
import hdf5plugin

import h5py
import pandas as pd


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

CDNA_MATCH_RE = re.compile(
    r"(c\.[0-9]+(?:_[0-9]+)?(?:[ACGT]+>[ACGT]+|delins[ACGT]+|del|ins[ACGT]+|dup)?)"
)


def canonical_spdi_to_chrpos(spdi: str) -> Optional[str]:
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


def normalize_function_call(value: str) -> str:
    if not value:
        return "unknown"

    value = value.strip()
    if not value:
        return "unknown"

    lc = value.lower()

    if "gof" in lc or "gain" in lc:
        if "mixed" in lc:
            return "GOF/Mixed"
        return "GOF"

    if "lof" in lc or "loss" in lc:
        if "mixed" in lc:
            return "LOF/Mixed"
        return "LOF"

    if "mixed" in lc:
        return "Mixed"

    if "neutral" in lc:
        return "Neutral"

    return value


def extract_cdna_from_text(text: str) -> Optional[str]:
    if not text:
        return None

    match = CDNA_MATCH_RE.search(text)
    if match:
        return match.group(1)

    return None


def load_clinvar_records(path: Path) -> Dict[str, dict]:
    with path.open("r", encoding="utf-8") as handle:
        data = json.load(handle)

    records = data.get("records", [])
    mapping: Dict[str, dict] = {}

    for record in records:
        title = record.get("title", "")
        cdna = record.get("cdna") or extract_cdna_from_text(title)

        if cdna:
            mapping[cdna] = record

    return mapping


def read_variant_rows(path: Path) -> Iterable[dict]:
    with path.open("r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            yield row


def build_case_row(
    row: dict,
    clinvar_record: dict,
    case_count: int,
    control_count: int,
) -> dict:
    spdi = clinvar_record.get("canonical_spdi", "")
    variant_id = canonical_spdi_to_chrpos(str(spdi))

    if variant_id is None:
        raise ValueError(f"Invalid canonical_spdi for ClinVar record: {spdi}")

    chrom, pos, ref, alt = variant_id.split(":")

    function_call = normalize_function_call(
        str(row.get("function_call", "") or row.get("function_basis", ""))
    )

    return {
        "Variant ID": variant_id,
        "chr": chrom,
        "pos": int(pos),
        "ref": ref,
        "alt": alt,
        "ac_case": case_count,
        "ac_control": control_count,
        "variant_source": "paper_functional_variant",
        "function_call": function_call,
        "cdna": row.get("cdna", ""),
        "protein_change": row.get("protein_change", ""),
        "clinvar_accession": clinvar_record.get("accession", ""),
        "clinvar_variation_id": clinvar_record.get("variation_id", ""),
        "clinvar_title": clinvar_record.get("title", ""),
    }


def decode_ascii_column(series: pd.Series) -> pd.Series:
    return series.apply(lambda x: x.decode("ascii") if isinstance(x, bytes) else str(x))


def load_all_missense_controls_for_uniprot(
    h5_path: Path,
    chrom: str,
    uniprot_id: str,
    control_count: int = 1,
) -> list[dict]:
    """Load all possible missense variants for one UniProt ID as controls."""

    with h5py.File(h5_path, "r") as f:
        ref_alt = f[f"{chrom}_ref_alt"][:]
        pdb_filename = f[f"{chrom}_filename"][:]
        uniprot_ids = f[f"{chrom}_uniprot_id"][:]
        positions = f[f"{chrom}_pos"][:]

    df = pd.DataFrame(
        {
            "ref": ref_alt[:, 0].flatten(),
            "alt": ref_alt[:, 1].flatten(),
            "aa_ref": ref_alt[:, 2].flatten(),
            "aa_alt": ref_alt[:, 3].flatten(),
            "pdb_filename": pdb_filename.flatten(),
            "uniprot_id": uniprot_ids.flatten(),
            "pos": positions[:, 0].flatten(),
            "aa_pos": positions[:, 1].flatten(),
            "aa_pos_file": positions[:, 2].flatten(),
        }
    )

    for col in ["ref", "alt", "aa_ref", "aa_alt", "pdb_filename", "uniprot_id"]:
        df[col] = decode_ascii_column(df[col])

    df = df[df["uniprot_id"] == uniprot_id].copy()

    rows = []

    for _, r in df.iterrows():
        rows.append(
            {
                "Variant ID": f"{chrom}:{int(r['pos'])}:{r['ref']}:{r['alt']}",
                "chr": chrom,
                "pos": int(r["pos"]),
                "ref": r["ref"],
                "alt": r["alt"],
                "ac_case": 0,
                "ac_control": control_count,
                "variant_source": "all_possible_missense_control",
                "function_call": "all_missense_control",
                "cdna": "",
                "protein_change": f"{r['aa_ref']}{int(r['aa_pos'])}{r['aa_alt']}",
                "clinvar_accession": "",
                "clinvar_variation_id": "",
                "clinvar_title": "",
            }
        )

    return rows


def deduplicate_rows(rows: list[dict]) -> list[dict]:
    """Merge duplicate Variant IDs by summing ac_case/ac_control and preserving metadata."""

    merged: dict[str, dict] = {}

    for row in rows:
        variant_id = row["Variant ID"]

        if variant_id not in merged:
            merged[variant_id] = row.copy()
            continue

        merged[variant_id]["ac_case"] += int(row.get("ac_case", 0))
        merged[variant_id]["ac_control"] += int(row.get("ac_control", 0))

        existing_source = merged[variant_id].get("variant_source", "")
        new_source = row.get("variant_source", "")

        if new_source and new_source not in existing_source:
            merged[variant_id]["variant_source"] = existing_source + ";" + new_source

        existing_function = merged[variant_id].get("function_call", "")
        new_function = row.get("function_call", "")

        if new_function and new_function not in existing_function:
            merged[variant_id]["function_call"] = existing_function + ";" + new_function

    return list(merged.values())


def write_tsv(path: Path, rows: Iterable[dict]) -> None:
    rows = list(rows)

    if not rows:
        print("No rows to write.", file=sys.stderr)
        return

    path.parent.mkdir(parents=True, exist_ok=True)

    fieldnames = [
        "Variant ID",
        "chr",
        "pos",
        "ref",
        "alt",
        "ac_case",
        "ac_control",
        "variant_source",
        "function_call",
        "cdna",
        "protein_change",
        "clinvar_accession",
        "clinvar_variation_id",
        "clinvar_title",
    ]

    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Create SCN1A RVAS input using functional variants as cases and all missense variants as controls."
    )

    parser.add_argument(
        "--exp-csv",
        default="data/exp_variance/SCN1A/SCN1A_variants_functional.csv",
        help="Path to SCN1A functional variant CSV.",
    )

    parser.add_argument(
        "--clinvar-json",
        default="data/clinvar/SCN1A/clinvar.json",
        help="Path to SCN1A ClinVar JSON.",
    )

    parser.add_argument(
        "--missense-h5",
        default="sir-reference-data/all_missense_variants_gr38.h5",
        help="Path to all_missense_variants_gr38.h5.",
    )

    parser.add_argument(
        "--output-tsv",
        default="input/SCN1A_cases_vs_all_missense.tsv",
        help="Path to output TSV file.",
    )

    parser.add_argument(
        "--chrom",
        default="chr2",
        help="Chromosome for SCN1A.",
    )

    parser.add_argument(
        "--uniprot-id",
        default="P35498",
        help="UniProt ID for SCN1A.",
    )

    parser.add_argument(
        "--case-count",
        type=int,
        default=1,
        help="Allele count assigned to each matched functional variant.",
    )

    parser.add_argument(
        "--control-count",
        type=int,
        default=1,
        help="Allele count assigned to each all-missense control variant.",
    )

    parser.add_argument(
        "--case-function",
        default="",
        help="Optional function_call filter for cases, e.g. GOF or LOF. Leave blank to include all matched functional variants.",
    )

    args = parser.parse_args()

    clinvar_map = load_clinvar_records(Path(args.clinvar_json))

    rows: list[dict] = []
    unmatched: list[str] = []

    for row in read_variant_rows(Path(args.exp_csv)):
        function_call = normalize_function_call(
            str(row.get("function_call", "") or row.get("function_basis", ""))
        )

        if args.case_function:
            wanted = args.case_function.lower()
            if wanted not in function_call.lower():
                continue

        cdna = str(row.get("cdna", "")).strip()

        if not cdna:
            label = row.get("title") or row.get("protein_change") or "<missing>"
            print(f"Skipping row with missing cdna: {label}", file=sys.stderr)
            continue

        record = clinvar_map.get(cdna)

        if record is None:
            unmatched.append(cdna)
            continue

        try:
            rows.append(
                build_case_row(
                    row=row,
                    clinvar_record=record,
                    case_count=args.case_count,
                    control_count=0,
                )
            )
        except ValueError as exc:
            print(f"Skipping cdna {cdna}: {exc}", file=sys.stderr)

    if unmatched:
        unique_unmatched = sorted(set(unmatched))
        print(
            f"Warning: {len(unmatched)} rows did not match ClinVar records. "
            f"Unique unmatched cDNAs: {unique_unmatched[:10]}"
            + ("..." if len(unique_unmatched) > 10 else ""),
            file=sys.stderr,
        )

    control_rows = load_all_missense_controls_for_uniprot(
        h5_path=Path(args.missense_h5),
        chrom=args.chrom,
        uniprot_id=args.uniprot_id,
        control_count=args.control_count,
    )

    rows.extend(control_rows)
    rows = deduplicate_rows(rows)

    write_tsv(Path(args.output_tsv), rows)

    n_case = sum(int(r["ac_case"]) for r in rows)
    n_control = sum(int(r["ac_control"]) for r in rows)

    print(f"Wrote {len(rows)} unique variants to {args.output_tsv}")
    print(f"Total ac_case: {n_case}")
    print(f"Total ac_control: {n_control}")


if __name__ == "__main__":
    main()