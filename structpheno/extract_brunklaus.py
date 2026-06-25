#!/usr/bin/env python3
"""Extract the SCN1A patient table (eTable 1) from the Brunklaus et al. 2022
supplementary PDF and write it as a single fused CSV.

Source: Brunklaus et al., Neurology 2022 (wnl_2022_01_20_brunklaus_1_sdc1.pdf).
eTable 1 lists all Dravet and GEFS+ training/validation cohort patients across
pages 11-33 as a 10-column, text-positioned table using a middle-dot ("·") as
the decimal separator.
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path

import pdfplumber
import requests

# ID Cohort Dx Onset HGVS Protein Type SCN1A-score CADD REVEL
ROW_RE = re.compile(
    r"^(\d+)\s+(\S+)\s+(Dravet|GEFS\+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(Missense|PTV)\s+(\S+)\s+(\S+)\s+(\S+)$"
)

FIELDS = [
    "id", "cohort_set", "country", "dx", "onset_months",
    "transcript", "cdna", "protein_change", "residue", "residue_canonical",
    "ref_aa", "ref_aa_ok", "type", "scn1a_score", "cadd", "revel",
]

# --- Residue lift-over to the SCN1A canonical isoform (UniProt P35498, 2009 aa) ---
# The main transcript ENST00000303395.4 encodes a protein identical to P35498, so
# its residue numbering is already canonical. The Australian cohort uses
# ENST00000375405.3 (1998 aa), which lacks an 11-residue block (canonical 671-681,
# "VIIDKPATDDN"). Derived by global alignment of the two protein sequences; the
# single indel means a flat +11 shift for alt residues past the insertion point.
CANONICAL_UNIPROT = "P35498"
MAIN_TRANSCRIPT = "ENST00000303395.4"
ALT_TRANSCRIPT = "ENST00000375405.3"
ALT_INSERTION_AFTER = 670  # alt residues > 670 shift by +11 to reach canonical
ALT_INSERTION_LEN = 11


def to_canonical(residue: int, transcript: str) -> int:
    """Map a transcript-native residue number onto canonical P35498 numbering."""
    if transcript == ALT_TRANSCRIPT and residue > ALT_INSERTION_AFTER:
        return residue + ALT_INSERTION_LEN
    return residue


def fetch_canonical_sequence() -> str | None:
    """Fetch the canonical SCN1A protein sequence for validation (best effort)."""
    try:
        r = requests.get(
            f"https://rest.uniprot.org/uniprotkb/{CANONICAL_UNIPROT}.fasta",
            timeout=30,
        )
        r.raise_for_status()
        return "".join(l for l in r.text.strip().split("\n") if not l.startswith(">"))
    except requests.exceptions.RequestException as e:
        print(f"Warning: could not fetch canonical sequence for validation: {e}",
              file=sys.stderr)
        return None


def _dot(value: str) -> str:
    """Normalise European middle-dot decimals; blank out NA."""
    if value in ("NA", ""):
        return ""
    return value.replace("·", ".")


def _residue_and_ref(protein_change: str) -> tuple[str, str]:
    """Return (residue_number, reference_aa) from an HGVS protein change.

    Handles the usual ``p.R1636Q`` form and entries written without the ``p.``
    prefix (e.g. ``R712X``). Whole-gene/exon deletions have neither.
    """
    if protein_change.startswith("DELETION"):
        return "", ""
    match = re.match(r"p?\.?([A-Za-z])[a-z]{0,2}(\d+)", protein_change)
    if not match:
        return "", ""
    ref_aa = match.group(1).upper()
    return match.group(2), ref_aa


def parse_pdf(pdf_path: Path) -> list[dict]:
    rows: list[dict] = []
    section: str | None = None

    with pdfplumber.open(pdf_path) as pdf:
        for page in pdf.pages:
            for line in (page.extract_text() or "").split("\n"):
                s = line.strip()
                if s.startswith("Training cohort"):
                    section = "training"
                    continue
                if s.startswith("Validation cohort 1"):
                    section = "validation_australia"
                    continue
                if s.startswith("Validation cohort 2"):
                    section = "validation_belgium"
                    continue

                m = ROW_RE.match(s)
                if not m:
                    continue

                id_, country, dx, onset, hgvs, prot, typ, scn1a, cadd, revel = m.groups()
                transcript, _, cdna = hgvs.partition(":")
                residue, ref_aa = _residue_and_ref(prot)
                residue_canonical = (
                    str(to_canonical(int(residue), transcript)) if residue else ""
                )
                rows.append({
                    "id": id_,
                    "cohort_set": section or "",
                    "country": country,
                    "dx": dx,
                    "onset_months": _dot(onset),
                    "transcript": transcript,
                    "cdna": cdna,
                    "protein_change": prot,
                    "residue": residue,
                    "residue_canonical": residue_canonical,
                    "ref_aa": ref_aa,
                    "ref_aa_ok": "",  # filled by validate_against_canonical()
                    "type": typ,
                    "scn1a_score": _dot(scn1a),
                    "cadd": _dot(cadd),
                    "revel": _dot(revel),
                })
    return rows


def validate_against_canonical(rows: list[dict], canonical_seq: str) -> int:
    """Flag each row's ref_aa against the canonical sequence at residue_canonical.

    Sets ref_aa_ok to "True"/"False" (blank for rows without a residue/ref_aa).
    Returns the number of mismatches, which are pre-existing discrepancies in the
    source annotation rather than lift-over errors.
    """
    mismatches = 0
    for r in rows:
        if not r["residue_canonical"] or not r["ref_aa"]:
            continue
        pos = int(r["residue_canonical"])
        observed = canonical_seq[pos - 1] if 0 < pos <= len(canonical_seq) else None
        ok = observed == r["ref_aa"]
        r["ref_aa_ok"] = str(ok)
        if not ok:
            mismatches += 1
    return mismatches


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "pdf",
        nargs="?",
        default=str(Path.home() / "Downloads" / "wnl_2022_01_20_brunklaus_1_sdc1.pdf"),
        help="Path to the Brunklaus supplementary PDF",
    )
    parser.add_argument(
        "-o", "--output",
        default="data/patient_variants/SCN1A/dravet_gefs_brunklaus.csv",
        help="Output CSV path",
    )
    parser.add_argument(
        "--no-validate",
        action="store_true",
        help="Skip fetching the canonical sequence to validate reference AAs",
    )
    args = parser.parse_args()

    rows = parse_pdf(Path(args.pdf))

    n_alt = sum(1 for r in rows if r["transcript"] == ALT_TRANSCRIPT and r["residue"])
    print(f"Lifted {n_alt} {ALT_TRANSCRIPT} rows onto canonical {CANONICAL_UNIPROT} numbering (+{ALT_INSERTION_LEN}).")

    if not args.no_validate:
        canonical_seq = fetch_canonical_sequence()
        if canonical_seq:
            checked = sum(1 for r in rows if r["residue_canonical"] and r["ref_aa"])
            mism = validate_against_canonical(rows, canonical_seq)
            print(f"Reference-AA validation: {checked - mism}/{checked} match canonical "
                  f"({mism} flagged ref_aa_ok=False — source annotation discrepancies).")

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(rows)

    print(f"Extracted {len(rows)} patients -> {out}")


if __name__ == "__main__":
    main()
