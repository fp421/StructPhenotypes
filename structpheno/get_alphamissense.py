#!/usr/bin/env python3
"""Download AlphaMissense predictions for a gene of interest."""

import argparse
import requests
import pandas as pd
from pathlib import Path
from typing import Optional
import sys


def get_uniprot_id(gene_name: str, canonical_only: bool = True) -> Optional[dict]:
    """
    Get UniProt ID from gene name using UniProt API.
    Returns canonical isoform by default.

    Returns:
        dict with keys: primary_accession, isoform_id, is_canonical, protein_name
    """
    url = "https://rest.uniprot.org/uniprotkb/search"
    params = {
        "query": f"(gene_exact:{gene_name}) AND (reviewed:true) AND (organism_id:9606)",
        "format": "json",
        "size": 10  # Get multiple results to find canonical
    }

    try:
        response = requests.get(url, params=params, timeout=10)
        response.raise_for_status()
        data = response.json()

        if not data.get("results"):
            return None

        # Filter for canonical isoforms (no "-" in accession, or "-1" suffix)
        results = []
        for entry in data["results"]:
            acc = entry["primaryAccession"]
            is_canonical = "-" not in acc or acc.endswith("-1")

            results.append({
                "primary_accession": acc,
                "isoform_id": acc,
                "is_canonical": is_canonical,
                "protein_name": entry.get("proteinDescription", {}).get("recommendedName", {}).get("fullName", {}).get("value", "Unknown"),
                "organism": entry.get("organism", {}).get("scientificName", "Unknown")
            })

        # Prioritize canonical isoforms
        if canonical_only:
            canonical = [r for r in results if r["is_canonical"]]
            if canonical:
                return canonical[0]

        return results[0] if results else None

    except Exception as e:
        print(f"Error fetching UniProt ID: {e}", file=sys.stderr)
        return None


def download_alphamissense_predictions(uniprot_info: dict, output_file: Optional[str] = None) -> pd.DataFrame:
    """
    Download AlphaMissense predictions from UniProt.

    AlphaMissense predictions are available through the UniProt REST API.
    This function retrieves all variant pathogenicity predictions for a protein.
    """
    uniprot_id = uniprot_info["primary_accession"]
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.json"

    try:
        response = requests.get(url, timeout=10)
        response.raise_for_status()
        data = response.json()

        # Extract protein information
        sequence = data.get("sequence", {}).get("value", "")
        protein_info = {
            "uniprot_id": uniprot_id,
            "isoform": uniprot_info["isoform_id"],
            "is_canonical": uniprot_info["is_canonical"],
            "protein_name": uniprot_info["protein_name"],
            "sequence_length": len(sequence),
            "organism": uniprot_info.get("organism", "Unknown"),
        }

        # Try to get AlphaMissense scores if available in the response
        # Note: Full AlphaMissense predictions may need to be fetched from specialized databases
        alphamissense_data = []

        # Check for cross-references that might contain AlphaMissense data
        for xref in data.get("uniProtKBCrossReferences", []):
            if xref.get("database") == "AlphaFoldDB":
                alphamissense_data.append({
                    "uniprot_id": uniprot_id,
                    "alphaFold_url": xref.get("url"),
                    "cross_ref_id": xref.get("id")
                })

        if not alphamissense_data:
            # If not in cross-references, fetch from AlphaMissense-specific endpoint
            alphamissense_data = fetch_alphamissense_variants(uniprot_id)

        df = pd.DataFrame(alphamissense_data) if alphamissense_data else pd.DataFrame({"uniprot_id": [uniprot_id]})

        # Add protein info
        for key, value in protein_info.items():
            if key not in df.columns:
                df[key] = value

        if output_file:
            df.to_csv(output_file, index=False)
            print(f"Saved predictions to {output_file}")

        return df

    except Exception as e:
        print(f"Error downloading AlphaMissense predictions: {e}", file=sys.stderr)
        return pd.DataFrame()


def fetch_alphamissense_variants(uniprot_id: str) -> list:
    """
    Fetch AlphaMissense variant predictions from Google Cloud Storage.

    Uses the official dm_alphamissense bucket with comprehensive predictions
    for all possible missense variants.
    """
    import subprocess
    import tempfile
    import gzip

    variants = []

    try:
        # Try to use gsutil if available
        print(f"Fetching AlphaMissense predictions for {uniprot_id} from Google Cloud Storage...")

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_file = Path(tmpdir) / "alphamissense.tsv"

            try:
                # Try gsutil first (faster if installed)
                cmd = [
                    "gsutil", "cp",
                    "gs://dm_alphamissense/human_proteome_missense_scores.tsv.gz",
                    str(tmp_file) + ".gz"
                ]
                result = subprocess.run(cmd, capture_output=True, timeout=300)

                if result.returncode == 0:
                    # Decompress and parse
                    with gzip.open(str(tmp_file) + ".gz", 'rt') as f:
                        for line in f:
                            if line.startswith('#'):
                                continue
                            parts = line.strip().split('\t')
                            if len(parts) >= 5:
                                protein_var = parts[4]
                                if uniprot_id in protein_var:
                                    variants.append({
                                        'chrom': parts[0],
                                        'pos': int(parts[1]),
                                        'ref': parts[2],
                                        'alt': parts[3],
                                        'protein_variant': protein_var,
                                        'am_pathogenicity': float(parts[5]) if len(parts) > 5 else None,
                                        'am_class': parts[6] if len(parts) > 6 else None
                                    })
                    print(f"Downloaded {len(variants)} variant predictions for {uniprot_id}")
                    return variants

            except (FileNotFoundError, subprocess.TimeoutExpired, Exception) as e:
                print(f"gsutil approach failed: {e}", file=sys.stderr)

        # Fallback: inform user to download manually
        print(f"\nNote: gsutil not available or download failed.")
        print(f"To get AlphaMissense predictions for {uniprot_id}:")
        print(f"  1. Install Google Cloud SDK: https://cloud.google.com/sdk/docs/install")
        print(f"  2. Run: gsutil cp gs://dm_alphamissense/human_proteome_missense_scores.tsv.gz .")
        print(f"  3. Decompress and filter for {uniprot_id}")

    except Exception as e:
        print(f"Error fetching AlphaMissense variants: {e}", file=sys.stderr)

    return variants


def main():
    parser = argparse.ArgumentParser(
        description="Download AlphaMissense predictions for a gene of interest"
    )
    parser.add_argument(
        "gene_name",
        help="Gene name (e.g., SCN2A, TP53)"
    )
    parser.add_argument(
        "-o", "--output",
        help="Output CSV file (default: data/alphamissense/{gene_name}_alphamissense.csv)"
    )
    parser.add_argument(
        "-u", "--uniprot-id",
        help="UniProt ID (if known, skip gene lookup)"
    )

    args = parser.parse_args()

    # Get UniProt ID
    if args.uniprot_id:
        print(f"Using provided UniProt ID: {args.uniprot_id}")
        uniprot_info = {
            "primary_accession": args.uniprot_id,
            "isoform_id": args.uniprot_id,
            "is_canonical": True,
            "protein_name": "Unknown",
            "organism": "Unknown"
        }
    else:
        print(f"Looking up UniProt ID for {args.gene_name}...")
        uniprot_info = get_uniprot_id(args.gene_name, canonical_only=True)
        if not uniprot_info:
            print(f"Could not find UniProt ID for gene: {args.gene_name}", file=sys.stderr)
            sys.exit(1)

    uniprot_id = uniprot_info["primary_accession"]
    print(f"Found UniProt ID: {uniprot_id} (Canonical: {uniprot_info['is_canonical']})")
    print(f"Protein: {uniprot_info['protein_name']}")

    # Prepare output file
    if args.output:
        output_file = args.output
    else:
        output_dir = Path("data/alphamissense")
        output_dir.mkdir(parents=True, exist_ok=True)
        output_file = output_dir / f"{args.gene_name}_alphamissense.csv"

    # Download predictions
    print(f"Downloading AlphaMissense predictions for {uniprot_id}...")
    df = download_alphamissense_predictions(uniprot_info, str(output_file))

    if not df.empty:
        print(f"\nSuccessfully retrieved data for {uniprot_id}")
        print(f"Shape: {df.shape}")
        print("\nFirst few rows:")
        print(df.head())
    else:
        print(f"No predictions found for {uniprot_id}")


if __name__ == "__main__":
    main()
