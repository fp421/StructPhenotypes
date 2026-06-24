"""

This module contains functions to retrieve AlphaFold PDB data for a gene.

"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any, Optional
import requests

ALPHAFOLD_PDB_TEMPLATE = "https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_v{version}.pdb"


class AlphaFoldRetriever:
    """Retrieve AlphaFold PDB data for a gene."""

    def __init__(self, gene: str, organism_id: int = 9606):
        self.gene = gene.strip()
        self.organism_id = organism_id

    def get_uniprot_id(self, canonical_only: bool = True) -> Optional[str]:
        return get_uniprot_id(self.gene, self.organism_id, canonical_only=canonical_only)

    def download_pdb(self, output_path: Optional[str | Path] = None, version: int = 6
                     ) -> Path:
        """Download the AlphaFold PDB file for this gene."""
        uniprot_id = self.get_uniprot_id()
        if uniprot_id is None:
            raise ValueError(f"Could not find a UniProt ID for gene '{self.gene}'")
        return download_alphafold_pdb(uniprot_id, output_path=output_path, version=version)

    def get_alpha_fold_data(self, output_path: Optional[str | Path] = None, version: int = 6) -> Path:
        """Download the AlphaFold PDB file for this gene and return its path."""
        return self.download_pdb(output_path=output_path, version=version)


def get_uniprot_id(gene_name: str, organism_id: int = 9606, canonical_only: bool = True) -> Optional[str]:
    """Return the UniProt accession for a gene symbol."""
    url = "https://rest.uniprot.org/uniprotkb/search"
    params = {
        "query": f"(gene_exact:{gene_name}) AND (reviewed:true) AND (organism_id:{organism_id})",
        "format": "json",
        "size": 10,
    }

    response = requests.get(url, params=params, timeout=10)
    response.raise_for_status()
    data = response.json()

    results = data.get("results", [])
    if not results:
        return None

    accessions = []
    for entry in results:
        accession = entry.get("primaryAccession")
        if not accession:
            continue
        is_canonical = "-" not in accession or accession.endswith("-1")
        accessions.append({"accession": accession, "canonical": is_canonical})

    if canonical_only:
        canonical = [a for a in accessions if a["canonical"]]
        if canonical:
            return canonical[0]["accession"]

    return accessions[0]["accession"] if accessions else None


def download_alphafold_pdb(
    uniprot_id: str,
    output_path: Optional[str | Path] = None,
    version: int = 6,
) -> Path:
    """Download an AlphaFold PDB file for a UniProt accession."""
    if output_path is None:
        output_path = Path(f"data/alphafold/SCN2A/AF-{uniprot_id}-F1-model_v{version}.pdb")
    else:
        output_path = Path(output_path)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    url = ALPHAFOLD_PDB_TEMPLATE.format(uniprot_id=uniprot_id, version=version)

    response = requests.get(url, timeout=20)

    if response.status_code != 200:
        raise FileNotFoundError(
            f"AlphaFold PDB not found from {url} (status {response.status_code})."
        )

    output_path.write_bytes(response.content)
    return output_path





    
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--gene", required=True)
    parser.add_argument("--out")
    args = parser.parse_args()

    retriever = AlphaFoldRetriever(args.gene)
    pdb_path = retriever.download_pdb(output_path=args.out)

    print(f"Downloaded AlphaFold PDB for {args.gene} to {pdb_path}")
