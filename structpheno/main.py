"""Command-line entry point for StructPhenotypes."""

from __future__ import annotations

import argparse
import json
import logging
import sys
import webbrowser
from pathlib import Path
from typing import Any, Callable
try:
    from tqdm import tqdm
except ImportError:
    tqdm = None

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
    from .paths import (
        alphafold_dir,
        alphamissense_mean_path,
        annotations_json_path,
        clinvar_json_path,
        normalize_gene,
        report_json_path,
    )
except ImportError:
    from paths import (
        alphafold_dir,
        alphamissense_mean_path,
        annotations_json_path,
        clinvar_json_path,
        normalize_gene,
        report_json_path,
    )

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


def build_report(
    gene: str,
    viewer_path: Path | None = None,
    *,
    logger: logging.Logger | None = None,
) -> dict[str, Any]:
    """Build one combined report for a gene symbol."""
    gene = normalize_gene(gene)
    logger = logger or logging.getLogger(__name__)
    steps = [
        ("ClinVar", lambda: _get_clinvar_data(gene, logger)),
        ("AlphaFold", lambda: AlphaFoldRetriever(gene).get_alpha_fold_data()),
        ("AlphaMissense", lambda: _get_alpha_missense_data(gene, logger)),
    ]
    progress = tqdm(total=len(steps) + 1, desc=f"{gene} pipeline", unit="step") if tqdm else None

    logger.info("Starting StructPhenotypes pipeline for %s", gene)
    if progress:
        progress.set_postfix_str("ClinVar")
    clinvar = _retrieve_source("ClinVar", steps[0][1])
    _log_source_result(logger, "ClinVar", clinvar)
    if progress:
        progress.update(1)
        progress.set_postfix_str("AlphaFold")

    alpha_fold = _retrieve_source("AlphaFold", steps[1][1])
    _log_source_result(logger, "AlphaFold", alpha_fold)
    if progress:
        progress.update(1)
        progress.set_postfix_str("AlphaMissense")

    alpha_missense = _retrieve_source("AlphaMissense", steps[2][1])
    _log_source_result(logger, "AlphaMissense", alpha_missense)
    if progress:
        progress.update(1)
        progress.set_postfix_str("Visualization")

    report = {
        "gene": gene,
        "clinvar": clinvar,
        "alpha_fold": alpha_fold,
        "alpha_missense": alpha_missense,
        "cache_paths": {
            "clinvar": str(clinvar_json_path(gene)),
            "annotations": str(annotations_json_path(gene)),
            "alphafold": str(alphafold_dir(gene)),
            "alphamissense": str(alphamissense_mean_path(gene)),
            "report": str(report_json_path(gene)),
        },
    }
    report["visualization"] = _retrieve_source(
        "Visualization",
        lambda: visualize(report, viewer_path, annotations_json_path(gene)),
    )
    _log_source_result(logger, "Visualization", report["visualization"])
    if progress:
        progress.update(1)
        progress.close()
    logger.info("Finished pipeline for %s", gene)
    return report


def _get_clinvar_data(gene: str, logger: logging.Logger | None = None) -> dict[str, Any]:
    """Load cached ClinVar data or fetch and cache it."""
    logger = logger or logging.getLogger(__name__)
    path = clinvar_json_path(gene)
    retriever = ClinVarRetriever(gene)
    if path.exists():
        cached = json.loads(path.read_text(encoding="utf-8"))
        if _clinvar_cache_matches_retriever(cached, retriever):
            logger.info("Using cached ClinVar data from %s", path)
            return cached
        logger.info(
            "ClinVar cache at %s was made with an older or different query; fetching fresh data",
            path,
        )

    logger.info("ClinVar cache miss for %s; fetching from ClinVar", gene)
    records = retriever.get_clinvar_data()
    retriever.save_clinvar_data(path, data=records)
    logger.info("Saved ClinVar data to %s", path)
    return json.loads(path.read_text(encoding="utf-8"))


def _clinvar_cache_matches_retriever(
    payload: dict[str, Any],
    retriever: ClinVarRetriever,
) -> bool:
    """Return whether cached ClinVar data was fetched with the current query."""
    query = payload.get("query")
    if isinstance(query, dict):
        return query.get("term") == retriever.query_metadata()["term"]

    # Legacy SCN1A caches were created with ``single_gene[prop]``, which excludes
    # valid C-terminal variants overlapping LOC102724058.
    return retriever.gene != "SCN1A"


def _get_alpha_missense_data(gene: str, logger: logging.Logger | None = None) -> Any:
    """Retrieve AlphaMissense data when the retriever class is available."""
    logger = logger or logging.getLogger(__name__)
    if AlphaMissenseRetriever is None:
        raise NotImplementedError(
            "AlphaMissenseRetriever is not available in get_alphamissense.py yet."
        )
    logger.info("Loading AlphaMissense data for %s", gene)
    return AlphaMissenseRetriever(gene).get_alpha_missense_data()


def _source_record_count(payload: dict[str, Any]) -> int | None:
    data = payload.get("data")
    if isinstance(data, list):
        return len(data)
    if isinstance(data, dict):
        if isinstance(data.get("records"), list):
            return len(data["records"])
        record_count = data.get("record_count")
        if isinstance(record_count, int):
            return record_count
    return None


def _log_source_result(logger: logging.Logger, source_name: str, payload: dict[str, Any]) -> None:
    status = payload.get("status", "unknown")
    if status == "ok":
        record_count = _source_record_count(payload)
        if record_count is None:
            logger.info("%s complete", source_name)
        else:
            logger.info("%s complete (%s records)", source_name, record_count)
    else:
        logger.warning("%s failed: %s", source_name, payload.get("error"))


def _configure_logging() -> logging.Logger:
    logging.basicConfig(
        level=logging.INFO,
        format="[%(levelname)s] %(message)s",
    )
    return logging.getLogger("structpheno")


def _print_report_summary(
    report: dict[str, Any],
    report_path: Path,
    viewer_path: Path | None,
) -> None:
    print(f"Gene: {report['gene']}")
    print(f"Report JSON: {report_path}")
    if viewer_path:
        print(f"Viewer HTML: {viewer_path}")
    for source_name in ("clinvar", "alpha_fold", "alpha_missense", "visualization"):
        payload = report.get(source_name, {})
        if not isinstance(payload, dict):
            continue
        label = source_name.replace("_", " ").title()
        status = payload.get("status", "unknown")
        count = _source_record_count(payload)
        suffix = f" ({count} records)" if count is not None else ""
        print(f"- {label}: {status}{suffix}")


def default_viewer_path(gene: str) -> Path:
    """Return the default local viewer path for a gene."""
    return Path("outputs") / normalize_gene(gene) / "viewer.html"


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
    logger = _configure_logging()
    viewer_path = args.viewer or default_viewer_path(args.gene)
    report = build_report(args.gene, viewer_path, logger=logger)

    if not args.no_open:
        open_viewer(report)

    json_indent = 2 if args.pretty else None
    json_text = json.dumps(report, indent=json_indent)
    report_path = args.out or report_json_path(normalize_gene(args.gene))
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(json_text + "\n", encoding="utf-8")
    logger.info("Wrote report JSON to %s", report_path)
    _print_report_summary(report, report_path, viewer_path)

    alpha_fold_status = report.get("alpha_fold", {}).get("status")
    visualization_status = report.get("visualization", {}).get("status")
    if alpha_fold_status != "ok" or visualization_status != "ok":
        logger.error("Pipeline failed: AlphaFold structure/viewer generation did not complete successfully.")
        return 1

    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
