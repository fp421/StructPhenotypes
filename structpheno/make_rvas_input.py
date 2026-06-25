"""Build a structure-informed-rvas case/control input TSV from ClinVar + gnomAD.

Cases are ClinVar pathogenic/likely-pathogenic missense SNVs; controls are gnomAD
variants. The output has the ``chr/pos/ref/alt/ac_case/ac_control`` columns that
``structure-informed-rvas/run.py --rvas-data-to-map`` expects. Both sources are
GRCh38. Only missense survives because run.py inner-joins against its reference
missense table, so gnomAD intronic/synonymous variants drop out during mapping.

Example:
    python -m structpheno.make_rvas_input SCN2A
    # -> data/rvas_input/SCN2A/clinvar_gnomad.tsv
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

try:
    from .paths import clinvar_json_path, gnomad_variants_path, rvas_input_path
except ImportError:
    from paths import clinvar_json_path, gnomad_variants_path, rvas_input_path


# Significance labels that count as "case", keyed by mode.
SIGNIFICANCE_SETS = {
    "plp": {"Pathogenic", "Likely pathogenic", "Pathogenic/Likely pathogenic"},
    "pathogenic": {"Pathogenic", "Pathogenic/Likely pathogenic"},
}


def refseq_to_chrom(accession: str) -> str | None:
    """Map a RefSeq chromosome accession (e.g. ``NC_000002.12``) to ``chr2``.

    NC_000001-22 -> chr1-22, NC_000023 -> chrX, NC_000024 -> chrY,
    NC_012920 -> chrM. Returns None for unrecognised accessions.
    """
    base = accession.split(".")[0]
    if not base.startswith("NC_"):
        return None
    if base == "NC_012920":
        return "chrM"
    try:
        n = int(base[3:])
    except ValueError:
        return None
    if 1 <= n <= 22:
        return f"chr{n}"
    if n == 23:
        return "chrX"
    if n == 24:
        return "chrY"
    return None


def build_cases(clinvar_path: Path, significance: set[str]) -> pd.DataFrame:
    """Return ClinVar P/LP missense SNVs as case rows (ac_case=1, ac_control=0).

    The ClinVar ``canonical_spdi`` is ``accession:position:deletion:insertion``
    with a 0-based position, so genomic pos = spdi_pos + 1. Only single-base
    ref/alt (SNVs) are kept, excluding indels/frameshifts.
    """
    records = json.loads(clinvar_path.read_text())["records"]
    rows = []
    for r in records:
        if "missense variant" not in r.get("molecular_consequences", []):
            continue
        if r.get("significance") not in significance:
            continue
        spdi = r.get("canonical_spdi")
        if not spdi or spdi.count(":") != 3:
            continue
        acc, pos0, ref, alt = spdi.split(":")
        if len(ref) != 1 or len(alt) != 1:  # SNVs only
            continue
        chrom = refseq_to_chrom(acc)
        if chrom is None:
            continue
        rows.append(
            {
                "chr": chrom,
                "pos": int(pos0) + 1,  # SPDI 0-based -> 1-based genomic
                "ref": ref,
                "alt": alt,
                "ac_case": 1,
                "ac_control": 0,
                "clinvar_residue": r.get("residue"),  # for downstream QC
            }
        )
    df = pd.DataFrame(rows)
    if not df.empty:
        df = df.drop_duplicates(subset=["chr", "pos", "ref", "alt"])
    return df


def build_controls(gnomad_path: Path) -> pd.DataFrame:
    """Return gnomAD variants as control rows (ac_case=0, ac_control=1).

    The gnomAD ``variant_id`` is already ``chrom-pos-ref-alt`` (chrom without the
    ``chr`` prefix). All consequences are kept; run.py keeps only missense during
    reference mapping. Binary weighting: each unique variant counts once.
    """
    g = pd.read_csv(gnomad_path)
    parts = g["variant_id"].str.split("-", expand=True)
    df = pd.DataFrame(
        {
            "chr": "chr" + parts[0],
            "pos": parts[1].astype(int),
            "ref": parts[2],
            "alt": parts[3],
            "ac_case": 0,
            "ac_control": 1,
        }
    )
    return df.drop_duplicates(subset=["chr", "pos", "ref", "alt"])


def make_rvas_input(
    gene: str,
    significance_mode: str = "plp",
    output: Path | None = None,
) -> Path:
    """Build and write the case/control TSV for ``gene``; return the output path."""
    significance = SIGNIFICANCE_SETS[significance_mode]
    df_case = build_cases(clinvar_json_path(gene), significance)
    df_ctrl = build_controls(gnomad_variants_path(gene))
    print(f"ClinVar {significance_mode} missense SNV cases: {len(df_case)}")
    print(f"gnomAD control variants (all consequences): {len(df_ctrl)}")

    out_cols = ["chr", "pos", "ref", "alt", "ac_case", "ac_control"]
    df = pd.concat(
        [df_case[out_cols] if not df_case.empty else df_case, df_ctrl[out_cols]],
        ignore_index=True,
    )

    out_path = Path(output) if output is not None else rvas_input_path(gene)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_path, sep="\t", index=False)
    print(f"Wrote {len(df)} rows -> {out_path}")

    # Sidecar mapping ClinVar's stated residue per case variant, used to QC the
    # SPDI->genomic offset after run.py mapping (mapped aa_pos should match).
    if not df_case.empty:
        qc_path = out_path.with_name(out_path.stem + "_clinvar_resid.tsv")
        df_case[["chr", "pos", "ref", "alt", "clinvar_residue"]].to_csv(
            qc_path, sep="\t", index=False
        )
        print(f"Wrote QC residue map -> {qc_path}")
    return out_path


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("gene", help="Gene symbol, e.g. SCN2A")
    parser.add_argument(
        "--significance",
        choices=sorted(SIGNIFICANCE_SETS),
        default="plp",
        help="Which ClinVar significance labels count as cases (default: plp).",
    )
    parser.add_argument(
        "-o",
        "--output",
        help="Output TSV path (default: data/rvas_input/<GENE>/clinvar_gnomad.tsv).",
    )
    args = parser.parse_args()
    make_rvas_input(args.gene, args.significance, args.output)


if __name__ == "__main__":
    main()
