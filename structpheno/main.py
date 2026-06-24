"""Command-line entry point for StructPhenotypes."""

from __future__ import annotations

import argparse
import json
import sys
import webbrowser
from pathlib import Path
from typing import Any, Callable

try:
    from .get_alpha_fold import AlphaFoldRetriever
except ImportError:
    from get_alpha_fold import AlphaFoldRetriever

try:
    from .get_clinvar import ClinVarRetriever
except ImportError:
    from get_clinvar import ClinVarRetriever

try:
    from .visualize import visualize
except ImportError:
    from visualize import visualize

try:
    from .get_alphamissense import AlphaMissenseRetriever
except ImportError:
    try:
        from get_alphamissense import AlphaMissenseRetriever
    except ImportError:
        AlphaMissenseRetriever = None


def _retrieve_source(source_name: str, retrieve: Callable[[], Any]) -> dict[str, Any]:
    """Run one data-source retriever without letting it crash the whole report."""
    try:
        return {
            "status": "ok",
            "data": retrieve(),
            "error": None,
        }
    except Exception as exc:
        return {
            "status": "error",
            "data": [],
            "error": f"{source_name} retrieval failed: {exc}",
        }


def build_report(gene: str, viewer_path: Path | None = None) -> dict[str, Any]:
    """Build one combined report for a gene symbol."""
    gene = gene.strip()

    clinvar = _retrieve_source("ClinVar", lambda: ClinVarRetriever(gene).get_clinvar_data())

    alpha_fold = _retrieve_source("AlphaFold", lambda: AlphaFoldRetriever(gene).get_alpha_fold_data())

    alpha_missense = _retrieve_source("AlphaMissense", lambda: _get_alpha_missense_data(gene))

    report = {
        "gene": gene,
        "clinvar": clinvar,
        "alpha_fold": alpha_fold,
        "alpha_missense": alpha_missense,
    }
    report["visualization"] = _retrieve_source(
        "Visualization",
        lambda: visualize(report, viewer_path),
    )
    return report


def _get_alpha_missense_data(gene: str) -> Any:
    """Retrieve AlphaMissense data when the retriever class is available."""
    if AlphaMissenseRetriever is None:
        raise NotImplementedError(
            "AlphaMissenseRetriever is not available in get_alphamissense.py yet."
        )
    return AlphaMissenseRetriever(gene).get_alpha_missense_data()


def default_viewer_path(gene: str) -> Path:
    """Return the default local viewer path for a gene."""
    safe_gene = "".join(char if char.isalnum() or char in "-_" else "_" for char in gene.strip())
    safe_gene = safe_gene or "protein"
    return Path("outputs") / f"{safe_gene.lower()}_viewer.html"


def open_viewer(report: dict[str, Any]) -> None:
    """Open a generated viewer HTML file in the default browser."""
    visualization = report.get("visualization", {})
    visualization_data = visualization.get("data", {}) if isinstance(visualization, dict) else {}
    if not isinstance(visualization_data, dict):
        return

    html_path = visualization_data.get("html_path")
    if not html_path:
        visualization_data["html_opened"] = False
        visualization_data["html_open_error"] = "No viewer HTML path was generated."
        return

    try:
        opened = webbrowser.open(Path(html_path).resolve().as_uri())
    except Exception as exc:
        visualization_data["html_opened"] = False
        visualization_data["html_open_error"] = str(exc)
        return

    visualization_data["html_opened"] = opened
    visualization_data["html_open_error"] = None if opened else "Browser did not accept the viewer URL."


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Retrieve StructPhenotypes data for a gene.",
    )
    parser.add_argument(
        "gene",
        help="Gene symbol to query, for example SCN2A.",
    )
    parser.add_argument(
        "--out",
        type=Path,
        help="Optional path to write the JSON report.",
    )
    parser.add_argument(
        "--pretty",
        action="store_true",
        help="Print indented JSON to stdout.",
    )
    parser.add_argument(
        "--viewer",
        type=Path,
        help="Optional path to write a local 3Dmol.js HTML viewer.",
    )
    parser.add_argument(
        "--no-open",
        action="store_true",
        help="Generate the viewer HTML but do not open it in a browser.",
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    """Run the StructPhenotypes CLI."""
    args = parse_args(argv)
    viewer_path = args.viewer or default_viewer_path(args.gene)
    report = build_report(args.gene, viewer_path)

    if not args.no_open:
        open_viewer(report)

    json_indent = 2 if args.pretty else None
    json_text = json.dumps(report, indent=json_indent)

    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(json_text + "\n", encoding="utf-8")
    else:
        print(json_text)

    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
