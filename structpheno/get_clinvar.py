"""Retrieve ClinVar variant annotations for a gene."""

from __future__ import annotations

import csv
import json
import os
import re
import sys
import time
from typing import Any

import requests
import tqdm

try:
    from .paths import clinvar_json_path
except ImportError:
    from paths import clinvar_json_path

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

THREE_TO_ONE_AA = {
    "Ala": "A",
    "Arg": "R",
    "Asn": "N",
    "Asp": "D",
    "Cys": "C",
    "Gln": "Q",
    "Glu": "E",
    "Gly": "G",
    "His": "H",
    "Ile": "I",
    "Leu": "L",
    "Lys": "K",
    "Met": "M",
    "Phe": "F",
    "Pro": "P",
    "Ser": "S",
    "Thr": "T",
    "Trp": "W",
    "Tyr": "Y",
    "Val": "V",
    "Ter": "*",
}


class ClinVarRetriever:
    """Retrieve ClinVar variant summaries for one gene symbol."""

    def __init__(
        self,
        gene: str,
        *,
        email: str | None = None,
        api_key: str | None = None,
        tool: str = "StructPhenotypes",
        retmax: int = 500,
        single_gene: bool = True,
        request_delay: float | None = None,
        timeout: int = 60,
    ) -> None:
        """Initialize a ClinVar retriever.

        Args:
            gene: Gene symbol to query, for example ``SCN2A``.
            email: Contact email for NCBI E-utilities. Defaults to ``NCBI_EMAIL``.
            api_key: Optional NCBI API key. Defaults to ``NCBI_API_KEY``.
            tool: Tool name sent to NCBI.
            retmax: Number of IDs to fetch per ESearch page.
            single_gene: Restrict ClinVar search to records associated with one gene.
            request_delay: Delay between API calls. Defaults to NCBI-safe timing.
            timeout: Request timeout in seconds.
        """
        self.gene = gene.strip().upper()
        self.email = email or os.getenv("NCBI_EMAIL")
        self.api_key = api_key or os.getenv("NCBI_API_KEY")
        self.tool = tool
        self.retmax = retmax
        self.single_gene = single_gene
        self.timeout = timeout
        self.session = requests.Session()
        self.request_delay = request_delay if request_delay is not None else self._default_request_delay()

    def get_clinvar_data(self, max_records: int | None = None) -> list[dict[str, Any]]:
        """Return variant annotations with residue, significance, and phenotypes.

        Args:
            max_records: Optional cap for quick smoke tests. ``None`` retrieves all
                matching records returned by ClinVar.

        Returns:
            A list of dictionaries. Each dictionary is JSON-serializable and
            includes residue-level protein change data when ClinVar provides it.
        """
        variation_ids = self.get_variation_ids(max_records=max_records)
        summaries = self.get_summaries(variation_ids)
        return [self._summary_to_variant(summary) for summary in summaries]

    def save_clinvar_data(
        self,
        output_path: str | os.PathLike[str],
        *,
        data: list[dict[str, Any]] | None = None,
        max_records: int | None = None,
    ) -> str:
        """Save ClinVar annotations for future processing.

        The format is inferred from the file extension:
        ``.json`` writes metadata plus all records, ``.jsonl`` writes one record
        per line, and ``.csv`` writes a flattened table.

        Args:
            output_path: Destination path ending in ``.json``, ``.jsonl``, or ``.csv``.
            data: Optional existing ClinVar data to save. If omitted, data is fetched.
            max_records: Optional fetch cap used only when ``data`` is omitted.

        Returns:
            The saved output path as a string.
        """
        records = data if data is not None else self.get_clinvar_data(max_records=max_records)
        path = os.fspath(output_path)
        os.makedirs(os.path.dirname(path) or ".", exist_ok=True)

        extension = os.path.splitext(path)[1].lower()
        if extension == ".json":
            self._write_json(path, records)
        elif extension == ".jsonl":
            self._write_jsonl(path, records)
        elif extension == ".csv":
            self._write_csv(path, records)
        else:
            raise ValueError("Unsupported output format. Use .json, .jsonl, or .csv.")

        return path

    def get_variation_ids(self, max_records: int | None = None) -> list[str]:
        """Search ClinVar and return matching Variation IDs."""
        ids: list[str] = []
        retstart = 0
        target_count: int | None = None

        while target_count is None or len(ids) < target_count:
            page_size = self.retmax
            if max_records is not None:
                page_size = min(page_size, max_records - len(ids))
            if page_size <= 0:
                break

            payload = self._request_json(
                "esearch.fcgi",
                {
                    "db": "clinvar",
                    "term": self._search_term(),
                    "retmode": "json",
                    "retstart": retstart,
                    "retmax": page_size,
                },
            )

            result = payload.get("esearchresult", {})
            total_count = int(result.get("count", 0))
            if target_count is None:
                target_count = total_count if max_records is None else min(total_count, max_records)

            page_ids = result.get("idlist", [])
            ids.extend(page_ids)
            if not page_ids:
                break

            retstart += len(page_ids)

        return ids[:target_count]

    def get_summaries(self, variation_ids: list[str], batch_size: int = 200) -> list[dict[str, Any]]:
        """Fetch ClinVar ESummary records for Variation IDs."""
        summaries: list[dict[str, Any]] = []

        for start in tqdm.tqdm(range(0, len(variation_ids), batch_size), desc="Fetching summaries"):
            batch = variation_ids[start : start + batch_size]
            payload = self._request_json(
                "esummary.fcgi",
                {
                    "db": "clinvar",
                    "id": ",".join(batch),
                    "retmode": "json",
                },
            )

            result = payload.get("result", {})
            for uid in result.get("uids", []):
                summary = result.get(uid)
                if isinstance(summary, dict):
                    summaries.append(summary)

        return summaries

    def _search_term(self) -> str:
        term = f"{self.gene}[gene]"
        if self.single_gene:
            term += " AND single_gene[prop]"
        return term

    def _request_json(self, endpoint: str, params: dict[str, Any]) -> dict[str, Any]:
        request_params = {
            **params,
            "tool": self.tool,
        }
        if self.email:
            request_params["email"] = self.email
        if self.api_key:
            request_params["api_key"] = self.api_key

        response = self.session.get(
            f"{EUTILS}/{endpoint}",
            params=request_params,
            timeout=self.timeout,
        )
        response.raise_for_status()
        self._sleep_for_rate_limit()
        return response.json()

    def _sleep_for_rate_limit(self) -> None:
        if self.request_delay > 0:
            time.sleep(self.request_delay)

    def _default_request_delay(self) -> float:
        if self.api_key:
            return 0.11
        return 0.34

    def _summary_to_variant(self, summary: dict[str, Any]) -> dict[str, Any]:
        germline = summary.get("germline_classification", {})
        protein_change = str(summary.get("protein_change") or "")
        protein_parts = self._parse_protein_change(protein_change, str(summary.get("title") or ""))
        variation = self._primary_variation(summary)
        phenotypes = self._extract_phenotypes(germline)

        return {
            "variation_id": summary.get("uid"),
            "accession": summary.get("accession"),
            "accession_version": summary.get("accession_version"),
            "gene": summary.get("gene_sort") or self.gene,
            "title": summary.get("title"),
            "variant_type": variation.get("variant_type") or summary.get("obj_type"),
            "canonical_spdi": variation.get("canonical_spdi"),
            "protein_change": protein_change or None,
            "amino_acid": protein_parts["reference_amino_acid"],
            "alternate_amino_acid": protein_parts["alternate_amino_acid"],
            "residue": protein_parts["residue"],
            "significance": germline.get("description") or None,
            "review_status": germline.get("review_status") or None,
            "last_evaluated": germline.get("last_evaluated") or None,
            "phenotypes": phenotypes,
            "phenotype_names": [
                phenotype["name"]
                for phenotype in phenotypes
                if phenotype.get("name")
            ],
            "molecular_consequences": summary.get("molecular_consequence_list", []),
        }

    def _write_json(self, path: str, records: list[dict[str, Any]]) -> None:
        payload = {
            "gene": self.gene,
            "source": "ClinVar",
            "record_count": len(records),
            "records": records,
        }
        with open(path, "w", encoding="utf-8") as handle:
            json.dump(payload, handle, indent=2)
            handle.write("\n")

    def _write_jsonl(self, path: str, records: list[dict[str, Any]]) -> None:
        with open(path, "w", encoding="utf-8") as handle:
            for record in records:
                handle.write(json.dumps(record) + "\n")

    def _write_csv(self, path: str, records: list[dict[str, Any]]) -> None:
        fieldnames = [
            "variation_id",
            "accession",
            "accession_version",
            "gene",
            "residue",
            "amino_acid",
            "alternate_amino_acid",
            "protein_change",
            "significance",
            "review_status",
            "last_evaluated",
            "phenotype_names",
            "phenotypes_json",
            "molecular_consequences",
            "variant_type",
            "canonical_spdi",
            "title",
        ]

        with open(path, "w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames)
            writer.writeheader()
            for record in records:
                writer.writerow(self._flatten_for_csv(record))

    @staticmethod
    def _flatten_for_csv(record: dict[str, Any]) -> dict[str, Any]:
        flattened = dict(record)
        flattened["phenotype_names"] = "|".join(record.get("phenotype_names", []))
        flattened["phenotypes_json"] = json.dumps(record.get("phenotypes", []))
        flattened["molecular_consequences"] = "|".join(record.get("molecular_consequences", []))
        flattened.pop("phenotypes", None)
        return flattened

    @staticmethod
    def _primary_variation(summary: dict[str, Any]) -> dict[str, Any]:
        variation_set = summary.get("variation_set", [])
        if variation_set and isinstance(variation_set[0], dict):
            return variation_set[0]
        return {}

    @staticmethod
    def _extract_phenotypes(classification: dict[str, Any]) -> list[dict[str, Any]]:
        phenotypes = []
        for trait in classification.get("trait_set", []):
            phenotypes.append(
                {
                    "name": trait.get("trait_name"),
                    "xrefs": trait.get("trait_xrefs", []),
                }
            )
        return phenotypes

    @staticmethod
    def _parse_protein_change(protein_change: str, title: str) -> dict[str, Any]:
        parsed = {
            "reference_amino_acid": None,
            "alternate_amino_acid": None,
            "residue": None,
        }

        # Prefer the title's HGVS.p, which uses ClinVar's representative
        # transcript and matches the canonical protein numbering. The
        # protein_change field can list several isoform-specific changes
        # (e.g. "W118R, W903R, ..., W932R"); a plain search there would grab
        # the first (118) rather than the canonical one (932).
        three_letter_match = re.search(
            r"p\.\(?(?P<ref>[A-Z][a-z]{2})(?P<pos>\d+)(?P<alt>[A-Z][a-z]{2}|Ter|fs|del|dup|ins|=)?\)?",
            title,
        )
        if three_letter_match:
            parsed["reference_amino_acid"] = THREE_TO_ONE_AA.get(three_letter_match.group("ref"))
            alt = three_letter_match.group("alt")
            parsed["alternate_amino_acid"] = THREE_TO_ONE_AA.get(alt, alt) if alt else None
            parsed["residue"] = int(three_letter_match.group("pos"))
            return parsed

        # Fall back to protein_change only when the title lacks an HGVS.p.
        one_letter_match = re.search(
            r"(?P<ref>[A-Z*])(?P<pos>\d+)(?P<alt>[A-Z*]|fs|del|dup|ins|=)?",
            protein_change,
        )
        if one_letter_match:
            parsed["reference_amino_acid"] = one_letter_match.group("ref")
            parsed["alternate_amino_acid"] = one_letter_match.group("alt")
            parsed["residue"] = int(one_letter_match.group("pos"))

        return parsed


if __name__ == "__main__":
    retriever = ClinVarRetriever("SCN1A")
    if not retriever.email:
        print(
            "Tip: set NCBI_EMAIL=you@example.org to identify this tool to NCBI E-utilities.",
            file=sys.stderr,
        )

    data = retriever.get_clinvar_data(max_records=10_000)
    output_path = retriever.save_clinvar_data(
        "data/clinvar/scn1a_clinvar.json",
        data=data,
    )
    print(f"Retrieved {len(data)} ClinVar records for SCN1A.")
    print(f"Saved ClinVar records to {output_path}.")
