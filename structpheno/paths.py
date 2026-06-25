"""Shared filesystem layout helpers for StructPhenotypes."""

from __future__ import annotations

from pathlib import Path


DATA_ROOT = Path("data")


def normalize_gene(gene: str) -> str:
    """Return a filesystem-friendly uppercase gene symbol."""
    cleaned = "".join(char if char.isalnum() or char in "-_" else "_" for char in gene.strip())
    return (cleaned or "GENE").upper()


def source_dir(source: str, gene: str) -> Path:
    """Return ``data/<source>/<GENE>``."""
    return DATA_ROOT / source / normalize_gene(gene)


def clinvar_json_path(gene: str) -> Path:
    """Return the canonical cached ClinVar JSON path."""
    return source_dir("clinvar", gene) / "clinvar.json"


def annotations_json_path(gene: str) -> Path:
    """Return the canonical preprocessed visualization annotations path."""
    return source_dir("annotations", gene) / "annotations.json"


def alphamissense_mean_path(gene: str) -> Path:
    """Return the canonical cached AlphaMissense mean pathogenicity CSV path."""
    return source_dir("alphamissense", gene) / "mean_pathogenicity.csv"


def alphamissense_variants_path(gene: str) -> Path:
    """Return the canonical cached AlphaMissense variant-level CSV path."""
    return source_dir("alphamissense", gene) / "variants.csv"


def alphafold_dir(gene: str) -> Path:
    """Return the canonical AlphaFold directory for a gene."""
    return source_dir("alphafold", gene)


def find_alphafold_pdb(gene: str) -> Path | None:
    """Return the first cached AlphaFold PDB for a gene, if present."""
    directory = alphafold_dir(gene)
    if not directory.exists():
        return None
    matches = sorted(directory.glob("*.pdb"))
    return matches[0] if matches else None


def alphafold_pdb_path(gene: str, uniprot_id: str, version: int = 6) -> Path:
    """Return the canonical AlphaFold PDB path for a UniProt accession."""
    return alphafold_dir(gene) / f"AF-{uniprot_id}-F1-model_v{version}.pdb"


def report_json_path(gene: str) -> Path:
    """Return the canonical combined report JSON path."""
    return source_dir("reports", gene) / "report.json"


def gnomad_variants_path(gene: str) -> Path:
    """Return the canonical cached gnomAD variant-level CSV path."""
    return source_dir("gnomad", gene) / "variants.csv"


def rvas_input_path(gene: str) -> Path:
    """Return the canonical structure-informed-rvas case/control input TSV path."""
    return source_dir("rvas_input", gene) / "clinvar_gnomad.tsv"
