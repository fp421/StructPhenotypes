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
    Download AlphaMissense predictions for a protein.

    Fetches variant-level predictions from Google Cloud Storage and protein
    metadata from UniProt to create comprehensive annotation dataset.
    """
    uniprot_id = uniprot_info["primary_accession"]

    try:
        # Fetch protein metadata from UniProt
        url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.json"
        response = requests.get(url, timeout=10)
        response.raise_for_status()
        data = response.json()

        # Extract protein sequence and info
        sequence = data.get("sequence", {}).get("value", "")
        protein_info = {
            "uniprot_id": uniprot_id,
            "isoform": uniprot_info["isoform_id"],
            "is_canonical": uniprot_info["is_canonical"],
            "protein_name": uniprot_info["protein_name"],
            "sequence_length": len(sequence),
            "organism": uniprot_info.get("organism", "Unknown"),
        }

        # Fetch AlphaMissense variant predictions
        print("Retrieving variant predictions...")
        variants = fetch_alphamissense_variants(uniprot_id)

        if variants:
            df = pd.DataFrame(variants)
            # Add protein info
            for key, value in protein_info.items():
                df[key] = value
        else:
            # If no variants found, return metadata only
            df = pd.DataFrame([protein_info])

        if output_file:
            df.to_csv(output_file, index=False)
            print(f"Saved {len(df)} records to {output_file}")

        return df

    except Exception as e:
        print(f"Error downloading AlphaMissense predictions: {e}", file=sys.stderr)
        return pd.DataFrame()


def load_alphamissense_from_file(filepath: str, uniprot_id: str) -> list:
    """
    Load and filter AlphaMissense predictions from a local TSV file.

    Expected format (can be gzipped):
        #CHROM  POS     REF     ALT     protein_variant am_pathogenicity am_class
    """
    import gzip

    variants = []
    open_fn = gzip.open if str(filepath).endswith('.gz') else open

    try:
        with open_fn(filepath, 'rt') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                parts = line.strip().split('\t')
                if len(parts) >= 5:
                    protein_var = parts[4]
                    if uniprot_id in protein_var:
                        try:
                            variants.append({
                                'chrom': parts[0],
                                'pos': int(parts[1]),
                                'ref': parts[2],
                                'alt': parts[3],
                                'protein_variant': protein_var,
                                'am_pathogenicity': float(parts[5]) if len(parts) > 5 else None,
                                'am_class': parts[6] if len(parts) > 6 else None
                            })
                        except (ValueError, IndexError):
                            continue
        print(f"Loaded {len(variants)} variant predictions for {uniprot_id}")
    except Exception as e:
        print(f"Error loading AlphaMissense file: {e}", file=sys.stderr)

    return variants


def fetch_alphamissense_variants(uniprot_id: str) -> list:
    """
    Fetch AlphaMissense variant predictions from Google Cloud Storage.

    Uses the official dm_alphamissense bucket with comprehensive predictions
    for all possible missense variants.
    """
    import subprocess
    import tempfile

    variants = []

    print(f"Fetching AlphaMissense predictions for {uniprot_id}...")

    # Check if cached file exists locally
    cache_path = Path("data/alphamissense/human_proteome_missense_scores.tsv.gz")
    if cache_path.exists():
        print(f"Using cached AlphaMissense file: {cache_path}")
        return load_alphamissense_from_file(str(cache_path), uniprot_id)

    # Try to use gsutil if available
    try:
        cmd = ["gsutil", "--version"]
        result = subprocess.run(cmd, capture_output=True, timeout=5)
        has_gsutil = result.returncode == 0
    except FileNotFoundError:
        has_gsutil = False

    if has_gsutil:
        print("Downloading from Google Cloud Storage...")
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_file = Path(tmpdir) / "alphamissense.tsv.gz"
            try:
                cmd = [
                    "gsutil", "-m", "cp",
                    "gs://dm_alphamissense/human_proteome_missense_scores.tsv.gz",
                    str(tmp_file)
                ]
                result = subprocess.run(cmd, capture_output=True, timeout=600)
                if result.returncode == 0:
                    return load_alphamissense_from_file(str(tmp_file), uniprot_id)
                else:
                    print(f"gsutil failed: {result.stderr.decode()}", file=sys.stderr)
            except Exception as e:
                print(f"gsutil download failed: {e}", file=sys.stderr)

    # If no gsutil, provide instructions
    print(f"\n⚠️  AlphaMissense data not found. To fetch predictions for {uniprot_id}:")
    print("\nOption 1: Install Google Cloud SDK (recommended)")
    print("  1. brew install google-cloud-sdk  # or follow: https://cloud.google.com/sdk/docs/install")
    print("  2. gcloud auth application-default login")
    print("  3. gsutil cp gs://dm_alphamissense/human_proteome_missense_scores.tsv.gz data/alphamissense/")
    print("  4. Re-run this script\n")
    print("Option 2: Manual download")
    print("  1. Visit: https://console.cloud.google.com/storage/browser/dm_alphamissense")
    print("  2. Download: human_proteome_missense_scores.tsv.gz")
    print("  3. Place in: data/alphamissense/")
    print("  4. Re-run this script\n")

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
