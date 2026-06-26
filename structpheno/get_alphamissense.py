#!/usr/bin/env python3
"""Download AlphaMissense predictions for a gene of interest."""

import argparse
import requests
import pandas as pd
from pathlib import Path
from typing import Optional
import sys

try:
    from .paths import (
        alphamissense_mean_path,
        alphamissense_variants_path,
        normalize_gene,
    )
except ImportError:
    from paths import (
        alphamissense_mean_path,
        alphamissense_variants_path,
        normalize_gene,
    )


class AlphaMissenseRetriever:
    """Retrieve local-first AlphaMissense predictions for one gene."""

    def __init__(self, gene: str, uniprot_id: str | None = None):
        self.gene = normalize_gene(gene)
        self.uniprot_id = uniprot_id

    def get_alpha_missense_data(self) -> list[dict]:
        """Return mean AlphaMissense pathogenicity by residue position."""
        cached_mean = alphamissense_mean_path(self.gene)
        if cached_mean.exists():
            return _load_mean_pathogenicity_csv(cached_mean)

        uniprot_info = self._get_uniprot_info()
        if not uniprot_info:
            return []

        variants = download_alphamissense_predictions(uniprot_info, None)
        if variants.empty:
            return []

        mean_path = alphamissense_mean_path(self.gene)
        mean_path.parent.mkdir(parents=True, exist_ok=True)
        aggregated = aggregate_by_position(variants, str(mean_path))
        return aggregated.to_dict(orient="records")

    def _get_uniprot_info(self) -> Optional[dict]:
        if self.uniprot_id:
            return {
                "primary_accession": self.uniprot_id,
                "isoform_id": self.uniprot_id,
                "is_canonical": True,
                "protein_name": "Unknown",
                "organism": "Unknown",
            }
        return get_uniprot_id(self.gene, canonical_only=True)


def _load_mean_pathogenicity_csv(path: str | Path) -> list[dict]:
    """Load cached per-position mean AlphaMissense scores."""
    df = pd.read_csv(path)
    if "position" in df.columns and "residue" not in df.columns:
        df["residue"] = df["position"]
    if "mean_pathogenicity" in df.columns and "score" not in df.columns:
        df["score"] = df["mean_pathogenicity"]
    if "mean_class" in df.columns and "class" not in df.columns:
        df["class"] = df["mean_class"]
    df = df.astype(object).where(pd.notnull(df), None)
    return df.to_dict(orient="records")


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

        return df

    except Exception as e:
        print(f"Error downloading AlphaMissense predictions: {e}", file=sys.stderr)
        return pd.DataFrame()


def load_alphamissense_from_file(filepath: str, uniprot_id: str) -> list:
    """
    Load and filter AlphaMissense predictions from a local TSV file.

    Expected format (can be gzipped):
        uniprot_id  protein_variant  am_pathogenicity  am_class
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
                if len(parts) >= 4:
                    uid = parts[0]
                    if uid == uniprot_id:
                        try:
                            protein_var = parts[1]
                            # Extract position from variant (e.g., "M1A" -> position 1)
                            import re
                            match = re.search(r'(\d+)', protein_var)
                            pos = int(match.group(1)) if match else None

                            variants.append({
                                'uniprot_id': uid,
                                'protein_variant': protein_var,
                                'position': pos,
                                'am_pathogenicity': float(parts[2]) if len(parts) > 2 else None,
                                'am_class': parts[3] if len(parts) > 3 else None
                            })
                        except (ValueError, IndexError):
                            continue
        if variants:
            print(f"Loaded {len(variants)} variant predictions for {uniprot_id}")
        else:
            print(f"No variants found for {uniprot_id} in AlphaMissense data")
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
    cache_path = Path("data/alphamissense/AlphaMissense_aa_substitutions.tsv.gz")
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
                    "gs://dm_alphamissense/AlphaMissense_aa_substitutions.tsv.gz",
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
    print(f"\nWarning: AlphaMissense data not found. To fetch predictions for {uniprot_id}:")
    print("\nOption 1: Download using gsutil (recommended, faster)")
    print("  gsutil -m cp gs://dm_alphamissense/AlphaMissense_aa_substitutions.tsv.gz data/alphamissense/")
    print("  Then re-run this script\n")
    print("Option 2: Manual download")
    print("  1. Visit: https://console.cloud.google.com/storage/browser/dm_alphamissense")
    print("  2. Download: AlphaMissense_aa_substitutions.tsv.gz")
    print("  3. Place in: data/alphamissense/")
    print("  4. Re-run this script\n")

    return variants


def build_report(gene_name: str, df: pd.DataFrame) -> dict:
    """
    Build a report dictionary compatible with visualize.py.

    Converts aggregated AlphaMissense data into the format expected by
    the 3Dmol.js visualization script.
    """
    # Extract residue and score columns
    if 'position' in df.columns and 'mean_pathogenicity' in df.columns:
        missense_data = df[['position', 'mean_pathogenicity']].dropna()
        predictions = [
            {
                "residue": int(row['position']),
                "score": float(row['mean_pathogenicity'])
            }
            for _, row in missense_data.iterrows()
        ]
    else:
        predictions = []

    report = {
        "gene": gene_name,
        "alpha_missense": {
            "data": predictions
        }
    }
    return report


def aggregate_by_position(df: pd.DataFrame, output_file: Optional[str] = None) -> pd.DataFrame:
    """
    Calculate mean pathogenicity score per amino acid position.

    Groups all variants at each position and computes average pathogenicity.
    Useful for 3D structure visualization.
    """
    if df.empty or 'position' not in df.columns:
        return df

    # Group by position and calculate mean pathogenicity
    aggregated = df.groupby('position').agg({
        'am_pathogenicity': ['mean', 'min', 'max', 'std', 'count'],
        'uniprot_id': 'first',
    }).reset_index()

    # Flatten column names
    aggregated.columns = ['position', 'mean_pathogenicity', 'min_pathogenicity',
                          'max_pathogenicity', 'std_pathogenicity', 'num_variants',
                          'uniprot_id']

    # Determine class based on mean score
    aggregated['mean_class'] = aggregated['mean_pathogenicity'].apply(
        lambda x: 'pathogenic' if x > 0.75 else ('benign' if x < 0.25 else 'ambiguous')
    )

    if output_file:
        aggregated.to_csv(output_file, index=False)
        print(f"Saved aggregated predictions to {output_file}")

    return aggregated


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
    parser.add_argument(
        "--variants",
        action="store_true",
        help="Output all variants (default is mean pathogenicity per position)"
    )
    parser.add_argument(
        "--report",
        action="store_true",
        help="Output as JSON report format compatible with visualize.py"
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
        output_file = alphamissense_variants_path(args.gene_name)
        output_file.parent.mkdir(parents=True, exist_ok=True)

    # Download predictions
    print(f"Downloading AlphaMissense predictions for {uniprot_id}...")
    df = download_alphamissense_predictions(uniprot_info, None)

    if not df.empty:
        # Default: aggregate by position
        if not args.variants:
            agg_output = str(output_file) if args.output else (
                alphamissense_mean_path(args.gene_name)
            )
            df = aggregate_by_position(df, agg_output)
            print(f"Aggregated to {len(df)} positions")
            output_dir = Path(agg_output).parent
        else:
            # Save the full variant data
            output_dir = Path(output_file).parent if args.output else (
                alphamissense_variants_path(args.gene_name).parent
            )
            output_dir.mkdir(parents=True, exist_ok=True)
            variants_output = output_dir / "variants.csv"
            df.to_csv(variants_output, index=False)
            print(f"Saved {len(df)} variant records to {variants_output}")

        print(f"Output directory: {output_dir}")
        print(f"Shape: {df.shape}")
        if not args.report:
            print("\nFirst few rows:")
            print(df.head())
    else:
        print(f"No predictions found for {uniprot_id}")


if __name__ == "__main__":
    main()
