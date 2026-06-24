#!/usr/bin/env python3
"""Fetch variant data from gnomAD."""

import argparse
import requests
import pandas as pd
from pathlib import Path
from typing import Optional, List
import json
import sys


def query_gnomad_api(gene_symbol: str) -> dict:
    """Query gnomAD API for variants in a gene.

    Args:
        gene_symbol: Gene symbol (e.g., "SCN2A")

    Returns:
        Dictionary with variants and their frequencies
    """
    url = "https://gnomad.broadinstitute.org/api"

    query = """
    {
        gene(gene_symbol: "%s") {
            symbol
            variants(first: 10000) {
                edges {
                    node {
                        variant_id
                        chrom
                        pos
                        ref
                        alt
                        allele_type
                        exome {
                            ac
                            an
                            af
                        }
                        genome {
                            ac
                            an
                            af
                        }
                    }
                }
            }
        }
    }
    """ % gene_symbol

    try:
        response = requests.post(
            url,
            json={"query": query},
            headers={"Content-Type": "application/json"},
            timeout=30
        )
        response.raise_for_status()
        return response.json()
    except Exception as e:
        print(f"Error querying gnomAD API: {e}", file=sys.stderr)
        return {}


def parse_gnomad_variants(data: dict) -> List[dict]:
    """Parse gnomAD API response into variant records.

    Args:
        data: Response from query_gnomad_api

    Returns:
        List of variant dictionaries
    """
    variants = []

    if "errors" in data:
        print(f"API Error: {data['errors']}", file=sys.stderr)
        return variants

    try:
        gene = data.get("data", {}).get("gene", {})
        if not gene:
            return variants

        for edge in gene.get("variants", {}).get("edges", []):
            variant = edge.get("node", {})

            exome_af = variant.get("exome", {}).get("af")
            genome_af = variant.get("genome", {}).get("af")
            af = exome_af if exome_af is not None else genome_af

            if af is not None:
                variants.append({
                    "variant_id": variant.get("variant_id"),
                    "chrom": variant.get("chrom"),
                    "pos": variant.get("pos"),
                    "ref": variant.get("ref"),
                    "alt": variant.get("alt"),
                    "allele_type": variant.get("allele_type"),
                    "allele_freq": af,
                    "exome_af": exome_af,
                    "genome_af": genome_af,
                })

    except Exception as e:
        print(f"Error parsing variants: {e}", file=sys.stderr)

    return variants


def query_gnomad_bq(gene_symbol: str, limit: int = 10000) -> Optional[pd.DataFrame]:
    """Query gnomAD via BigQuery for large-scale analysis.

    Requires: pip install google-cloud-bigquery
    Requires: gcloud auth application-default login

    Args:
        gene_symbol: Gene symbol (e.g., "SCN2A")
        limit: Maximum number of variants to return

    Returns:
        DataFrame with variants, or None if BigQuery unavailable
    """
    try:
        from google.cloud import bigquery
    except ImportError:
        print("Error: google-cloud-bigquery not installed", file=sys.stderr)
        print("Install with: pip install google-cloud-bigquery", file=sys.stderr)
        return None

    try:
        client = bigquery.Client()

        query = f"""
        SELECT
            variant_id,
            chrom,
            pos,
            ref,
            alt,
            allele_type,
            exome.ac as exome_ac,
            exome.an as exome_an,
            exome.af as exome_af,
            genome.ac as genome_ac,
            genome.an as genome_an,
            genome.af as genome_af,
            most_severe_consequence,
            transcript.gene_symbol as gene_symbol
        FROM
            `bigquery-public-data.gnomad.v4_variants` as t,
            UNNEST(transcript) as transcript
        WHERE
            transcript.gene_symbol = '{gene_symbol}'
            AND variant_type = 'snv'
        ORDER BY
            pos
        LIMIT {limit}
        """

        print("Querying BigQuery (this may take a minute)...")
        results = client.query(query).to_dataframe()
        return results

    except Exception as e:
        print(f"Error querying BigQuery: {e}", file=sys.stderr)
        return None


def download_gnomad_variants(gene_symbol: str, use_bq: bool = False) -> pd.DataFrame:
    """Download gnomAD variants for a gene.

    Args:
        gene_symbol: Gene symbol (e.g., "SCN2A")
        use_bq: Use BigQuery instead of API (requires authentication)

    Returns:
        DataFrame with variants
    """
    if use_bq:
        print(f"Fetching {gene_symbol} variants from BigQuery...")
        df = query_gnomad_bq(gene_symbol)
        if df is not None:
            return df
        else:
            print("Falling back to API...")

    print(f"Fetching {gene_symbol} variants from gnomAD API...")
    data = query_gnomad_api(gene_symbol)
    variants = parse_gnomad_variants(data)

    if variants:
        df = pd.DataFrame(variants)
        return df
    else:
        return pd.DataFrame()


def build_gnomad_report(gene_name: str, df: pd.DataFrame) -> dict:
    """Build gnomAD report compatible with visualize.py format.

    Args:
        gene_name: Gene name
        df: DataFrame with variant data

    Returns:
        Report dictionary
    """
    if df.empty:
        return {"gene": gene_name, "gnomad": {"data": []}}

    if "pos" in df.columns:
        grouped = df.groupby("pos").agg({
            "allele_freq": ["mean", "min", "max"],
            "variant_id": "count"
        }).reset_index()

        grouped.columns = ["position", "mean_af", "min_af", "max_af", "num_variants"]

        predictions = [
            {
                "residue": int(row["position"]),
                "score": float(row["mean_af"]),
                "num_variants": int(row["num_variants"])
            }
            for _, row in grouped.iterrows()
        ]
    else:
        predictions = []

    return {
        "gene": gene_name,
        "gnomad": {
            "data": predictions,
            "source": "gnomAD v4",
            "description": "Allele frequency in healthy populations"
        }
    }


def main():
    parser = argparse.ArgumentParser(
        description="Download gnomAD variants for a gene of interest"
    )
    parser.add_argument(
        "gene_name",
        help="Gene name (e.g., SCN2A, TP53)"
    )
    parser.add_argument(
        "-o", "--output",
        help="Output CSV file (default: data/gnomad/{gene_name}_variants.csv)"
    )
    parser.add_argument(
        "--bq",
        action="store_true",
        help="Use BigQuery instead of API (requires gcloud auth)"
    )
    parser.add_argument(
        "--report",
        action="store_true",
        help="Output as JSON report format compatible with visualize.py"
    )

    args = parser.parse_args()

    output_dir = Path("data/gnomad") / args.gene_name.upper()
    output_dir.mkdir(parents=True, exist_ok=True)

    df = download_gnomad_variants(args.gene_name, use_bq=args.bq)

    if not df.empty:
        print(f"Retrieved {len(df)} variants for {args.gene_name}")

        if args.output:
            output_file = args.output
        else:
            output_file = output_dir / "variants.csv"

        df.to_csv(output_file, index=False)
        print(f"Saved to {output_file}")

        if args.report:
            report = build_gnomad_report(args.gene_name, df)
            report_file = output_dir / "report.json"
            report_file.write_text(json.dumps(report, indent=2))
            print(f"Saved report to {report_file}")

        print(f"\nFirst few rows:")
        print(df.head())
    else:
        print(f"No variants found for {args.gene_name}")


if __name__ == "__main__":
    main()
