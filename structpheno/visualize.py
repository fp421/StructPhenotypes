"""Local 3Dmol.js visualization helpers for StructPhenotypes."""

from __future__ import annotations

import csv
from collections import Counter, defaultdict
import json
from html import escape
from pathlib import Path
import re
from typing import Any

try:
    from .paths import normalize_gene
except ImportError:
    from paths import normalize_gene


EXAMPLE_PDB = """\
ATOM      1  N   MET A   1       0.000   1.200   0.000  1.00 20.00           N
ATOM      2  CA  MET A   1       1.100   0.300   0.400  1.00 20.00           C
ATOM      3  C   MET A   1       2.300   1.100   0.900  1.00 20.00           C
ATOM      4  O   MET A   1       2.300   2.300   0.900  1.00 20.00           O
ATOM      5  N   GLU A   2       3.300   0.400   1.300  1.00 20.00           N
ATOM      6  CA  GLU A   2       4.600   1.000   1.800  1.00 20.00           C
ATOM      7  C   GLU A   2       5.500   0.000   2.500  1.00 20.00           C
ATOM      8  O   GLU A   2       5.200  -1.200   2.500  1.00 20.00           O
ATOM      9  N   LEU A   3       6.600   0.500   3.100  1.00 20.00           N
ATOM     10  CA  LEU A   3       7.600  -0.300   3.800  1.00 20.00           C
ATOM     11  C   LEU A   3       8.900   0.400   3.600  1.00 20.00           C
ATOM     12  O   LEU A   3       9.000   1.600   3.700  1.00 20.00           O
ATOM     13  N   LYS A   4       9.900  -0.400   3.300  1.00 20.00           N
ATOM     14  CA  LYS A   4      11.200   0.100   3.000  1.00 20.00           C
ATOM     15  C   LYS A   4      12.000  -0.900   2.200  1.00 20.00           C
ATOM     16  O   LYS A   4      11.600  -2.000   1.900  1.00 20.00           O
ATOM     17  N   ALA A   5      13.200  -0.400   1.900  1.00 20.00           N
ATOM     18  CA  ALA A   5      14.100  -1.200   1.100  1.00 20.00           C
ATOM     19  C   ALA A   5      15.300  -0.400   0.700  1.00 20.00           C
ATOM     20  O   ALA A   5      15.300   0.800   0.900  1.00 20.00           O
ATOM     21  N   GLY A   6      16.300  -1.100   0.100  1.00 20.00           N
ATOM     22  CA  GLY A   6      17.500  -0.500  -0.400  1.00 20.00           C
ATOM     23  C   GLY A   6      17.700  -0.700  -1.900  1.00 20.00           C
ATOM     24  O   GLY A   6      17.400  -1.700  -2.500  1.00 20.00           O
ATOM     25  N   ARG A   7      18.300   0.300  -2.500  1.00 20.00           N
ATOM     26  CA  ARG A   7      18.600   0.200  -3.900  1.00 20.00           C
ATOM     27  C   ARG A   7      20.100   0.100  -4.100  1.00 20.00           C
ATOM     28  O   ARG A   7      20.900   0.800  -3.500  1.00 20.00           O
ATOM     29  N   VAL A   8      20.400  -0.800  -4.900  1.00 20.00           N
ATOM     30  CA  VAL A   8      21.800  -1.000  -5.300  1.00 20.00           C
ATOM     31  C   VAL A   8      22.000  -0.800  -6.800  1.00 20.00           C
ATOM     32  O   VAL A   8      21.200  -1.200  -7.600  1.00 20.00           O
ATOM     33  N   SER A   9      23.100  -0.100  -7.100  1.00 20.00           N
ATOM     34  CA  SER A   9      23.500   0.100  -8.500  1.00 20.00           C
ATOM     35  C   SER A   9      24.700  -0.700  -8.900  1.00 20.00           C
ATOM     36  O   SER A   9      25.700  -0.700  -8.200  1.00 20.00           O
ATOM     37  N   THR A  10      24.600  -1.400 -10.000  1.00 20.00           N
ATOM     38  CA  THR A  10      25.700  -2.200 -10.500  1.00 20.00           C
ATOM     39  C   THR A  10      26.300  -1.500 -11.700  1.00 20.00           C
ATOM     40  O   THR A  10      25.600  -1.100 -12.600  1.00 20.00           O
ATOM     41  N   ILE A  11      27.600  -1.300 -11.700  1.00 20.00           N
ATOM     42  CA  ILE A  11      28.300  -0.700 -12.800  1.00 20.00           C
ATOM     43  C   ILE A  11      29.400  -1.600 -13.300  1.00 20.00           C
ATOM     44  O   ILE A  11      30.000  -2.300 -12.500  1.00 20.00           O
ATOM     45  N   GLN A  12      29.600  -1.500 -14.600  1.00 20.00           N
ATOM     46  CA  GLN A  12      30.600  -2.200 -15.400  1.00 20.00           C
ATOM     47  C   GLN A  12      31.400  -1.200 -16.200  1.00 20.00           C
ATOM     48  O   GLN A  12      30.900  -0.200 -16.600  1.00 20.00           O
TER
END
"""


EXAMPLE_ANNOTATIONS = {
    "clinvar_variants": [
        {
            "residue": 3,
            "amino_acid": "L",
            "significance": "Pathogenic",
            "phenotype": "Example seizure phenotype",
            "color": "#d73027",
        },
        {
            "residue": 7,
            "amino_acid": "R",
            "significance": "Likely pathogenic",
            "phenotype": "Example developmental phenotype",
            "color": "#fc8d59",
        },
        {
            "residue": 11,
            "amino_acid": "I",
            "significance": "Uncertain significance",
            "phenotype": "Example mixed phenotype",
            "color": "#fee08b",
        },
    ],
    "missense_predictions": [
        {"residue": 1, "score": 0.08},
        {"residue": 2, "score": 0.16},
        {"residue": 3, "score": 0.91},
        {"residue": 4, "score": 0.34},
        {"residue": 5, "score": 0.27},
        {"residue": 6, "score": 0.42},
        {"residue": 7, "score": 0.82},
        {"residue": 8, "score": 0.48},
        {"residue": 9, "score": 0.55},
        {"residue": 10, "score": 0.63},
        {"residue": 11, "score": 0.72},
        {"residue": 12, "score": 0.22},
    ],
}


PATHOGENIC_LABELS = {
    "Pathogenic",
    "Likely pathogenic",
    "Pathogenic/Likely pathogenic",
}

GENERIC_PHENOTYPES = {
    "not provided",
    "not specified",
    "See cases",
}

PHENOTYPE_PALETTE = [
    "#1f77b4",
    "#ff7f0e",
    "#2ca02c",
    "#d62728",
    "#9467bd",
    "#8c564b",
    "#e377c2",
    "#7f7f7f",
    "#bcbd22",
    "#17becf",
    "#393b79",
    "#637939",
    "#8c6d31",
    "#843c39",
    "#7b4173",
]

FUNCTION_CLASS_LABELS = {
    "gain-of-function": "Gain of function",
    "loss-of-function": "Loss of function",
}


def visualize(
    report: dict[str, Any],
    output_path: str | Path | None = None,
    annotation_path: str | Path | None = None,
    use_cached_annotations: bool = True,
) -> dict[str, Any]:
    """Prepare or write a local 3Dmol.js viewer for a StructPhenotypes report.

    Args:
        report: Combined report from ``main.build_report``.
        output_path: Optional HTML path. If omitted, no file is written.

    Returns:
        JSON-serializable metadata about the local viewer.
    """
    gene = str(report.get("gene", "unknown"))
    resolved_path = Path(output_path) if output_path else None
    resolved_annotation_path = (
        Path(annotation_path) if annotation_path else _default_annotation_path(gene, resolved_path)
    )
    annotations, annotation_source = _load_or_make_annotations(
        report,
        resolved_annotation_path,
        use_cached=use_cached_annotations,
    )
    pdb_text, structure_source = _pdb_text_for_report(report)

    metadata: dict[str, Any] = {
        "status": "stub ready",
        "viewer": "3Dmol.js",
        "gene": gene,
        "html_path": str(resolved_path) if resolved_path else None,
        "html_written": False,
        "annotation_path": str(resolved_annotation_path.resolve()),
        "annotation_source": annotation_source,
        "annotation_residue_count": annotations["residue_count"],
        "annotation_record_count": annotations["record_count"],
        "structure_source": structure_source,
        "example_structure": structure_source == "example",
        "data_sources": {
            "clinvar": _source_has_data(report, "clinvar"),
            "alpha_fold": _source_has_data(report, "alpha_fold"),
            "alpha_missense": _source_has_data(report, "alpha_missense"),
        },
        "controls": [
            "Base color mode: gray protein, default, AlphaMissense gradient/class, AlphaFold confidence, ClinVar pathogenicity, ClinVar/AlphaMissense surprise, disease contact network",
            "Rendering style: cartoon, trace, tube, cartoon plus sticks, sticks, lines, sphere/spacefill, cross, VDW/MS/SAS/SES surfaces, plus a smoother SAS-based surface option",
            "Phenotype checkbox overlay with per-phenotype color pickers and hidden pathogenicity override menus",
            "Global pathogenicity display/filter mode inherited by phenotype overlays",
            "Contact cutoff slider for distance-based disease neighborhoods",
            "Disease residue summary statistics",
            "Dynamic color legend",
            "Residue focus with 30 percent opacity for the rest of the protein",
        ],
    }

    if resolved_path:
        resolved_path.parent.mkdir(parents=True, exist_ok=True)
        resolved_path.write_text(_render_html(report, annotations, pdb_text), encoding="utf-8")
        metadata["html_path"] = str(resolved_path.resolve())
        metadata["html_written"] = True

    return metadata


def make_visualization_annotations(report: dict[str, Any]) -> dict[str, Any]:
    """Build residue-level annotations for the 3D viewer."""
    records = _extract_clinvar_records(report)
    gene = str(report.get("gene") or _infer_gene(records) or "unknown")
    experimental_function = _load_experimental_function_data(gene)
    records_by_residue = _group_records_by_residue(records)
    phenotype_counts = _count_phenotypes(records)
    phenotype_colors = _assign_phenotype_colors(phenotype_counts)
    missense_by_residue, missense_source = _build_missense_gradient(records_by_residue, report)
    confidence_by_residue, confidence_source = _extract_alphafold_confidence(report)

    residues: dict[str, dict[str, Any]] = {}
    residue_list: list[dict[str, Any]] = []
    all_residues = sorted(
        set(records_by_residue)
        | {int(residue) for residue in missense_by_residue}
        | {int(residue) for residue in confidence_by_residue}
    )
    for residue in all_residues:
        residue_key = str(residue)
        annotation = _build_residue_annotation(
            residue,
            records_by_residue.get(residue, []),
            phenotype_colors,
            missense_by_residue.get(residue_key),
            confidence_by_residue.get(residue_key),
            experimental_function,
        )
        residues[residue_key] = annotation
        residue_list.append(annotation)

    return {
        "gene": gene,
        "source": "StructPhenotypes visualization preprocessing",
        "annotation_version": 6,
        "record_count": len(records),
        "residue_count": len(residues),
        "phenotype_counts": dict(phenotype_counts),
        "phenotype_colors": phenotype_colors,
        "pathogenic_labels": sorted(PATHOGENIC_LABELS),
        "function_class_labels": FUNCTION_CLASS_LABELS,
        "function_class_counts": experimental_function["class_counts"],
        "residues": residues,
        "residue_list": residue_list,
        "experimental_function": {
            "source": experimental_function["source"],
            "path": experimental_function["path"],
            "paths": experimental_function.get("paths", []),
            "patient_row_count": experimental_function["patient_row_count"],
            "variant_count": experimental_function["variant_count"],
            "residue_count": experimental_function["residue_count"],
            "class_counts": experimental_function["class_counts"],
            "source_counts": experimental_function.get("source_counts", {}),
        },
        "missense_predictions": {
            "source": missense_source,
            "score_range": [0.0, 1.0],
            "by_residue": missense_by_residue,
        },
        "alphafold_confidence": {
            "source": confidence_source,
            "score_range": [0.0, 100.0],
            "by_residue": confidence_by_residue,
        },
    }


def save_visualization_annotations(
    annotations: dict[str, Any],
    output_path: str | Path,
) -> Path:
    """Save preprocessed visualization annotations as JSON."""
    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(annotations, indent=2) + "\n", encoding="utf-8")
    return path


def preprocess_clinvar_json_file(
    clinvar_json_path: str | Path,
    output_path: str | Path | None = None,
) -> dict[str, Any]:
    """Preprocess a saved ClinVar JSON file and save residue annotations."""
    input_path = Path(clinvar_json_path)
    report = json.loads(input_path.read_text(encoding="utf-8"))
    resolved_output = Path(output_path) if output_path else input_path.with_name(
        f"{input_path.stem}_annotations.json"
    )
    annotations = make_visualization_annotations(report)
    save_visualization_annotations(annotations, resolved_output)
    return annotations


def _load_or_make_annotations(
    report: dict[str, Any],
    annotation_path: Path,
    *,
    use_cached: bool,
) -> tuple[dict[str, Any], str]:
    """Load cached annotations when present, otherwise preprocess and save."""
    if use_cached and annotation_path.exists():
        annotations = json.loads(annotation_path.read_text(encoding="utf-8"))
        if _annotation_cache_is_usable(annotations, report):
            return annotations, "local"

    annotations = make_visualization_annotations(report)
    save_visualization_annotations(annotations, annotation_path)
    return annotations, "generated"


def _annotation_cache_is_usable(annotations: dict[str, Any], report: dict[str, Any]) -> bool:
    """Return whether cached annotations match the current report requirements."""
    expected_gene = str(report.get("gene", "unknown")).upper()
    cached_gene = str(annotations.get("gene", "unknown")).upper()
    if cached_gene != expected_gene:
        return False
    if int(annotations.get("annotation_version", 0)) < 6:
        return False

    if _source_has_data(report, "alpha_missense"):
        missense = annotations.get("missense_predictions", {})
        if missense.get("source") != "alphamissense":
            return False

    if _source_has_data(report, "alpha_fold"):
        confidence = annotations.get("alphafold_confidence", {})
        if confidence.get("source") != "alphafold_pdb_b_factor":
            return False

    function_source_paths = _experimental_function_source_paths(expected_gene)
    experimental_function = annotations.get("experimental_function", {})
    if function_source_paths:
        if experimental_function.get("source") in {None, "none"}:
            return False
    elif experimental_function.get("source") not in {None, "none"}:
        return False

    return True


def _pdb_text_for_report(report: dict[str, Any]) -> tuple[str, str]:
    """Return PDB text from AlphaFold metadata or an embedded demo structure."""
    alpha_fold = report.get("alpha_fold", {})
    if isinstance(alpha_fold, dict):
        data = alpha_fold.get("data", {})
        if isinstance(data, dict):
            pdb_path = data.get("pdb_path")
            if pdb_path:
                path = Path(pdb_path)
                if path.exists():
                    return path.read_text(encoding="utf-8", errors="replace"), "alphafold"
                raise FileNotFoundError(f"AlphaFold PDB path does not exist: {path}")

    return EXAMPLE_PDB, "example"


def _extract_alphafold_confidence(report: dict[str, Any]) -> tuple[dict[str, dict[str, Any]], str]:
    """Extract per-residue AlphaFold confidence from PDB B-factor values."""
    alpha_fold = report.get("alpha_fold", {})
    if not isinstance(alpha_fold, dict):
        return {}, "none"

    data = alpha_fold.get("data", {})
    if not isinstance(data, dict):
        return {}, "none"

    pdb_path = data.get("pdb_path")
    if not pdb_path:
        return {}, "none"

    path = Path(pdb_path)
    if not path.exists():
        return {}, "none"

    residue_scores: dict[int, list[float]] = defaultdict(list)
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if not line.startswith("ATOM"):
            continue
        residue = _coerce_residue(line[22:26].strip())
        score = _coerce_score_100(line[60:66].strip())
        if residue is not None and score is not None:
            residue_scores[residue].append(score)

    confidence: dict[str, dict[str, Any]] = {}
    for residue, scores in sorted(residue_scores.items()):
        mean_score = sum(scores) / len(scores)
        confidence[str(residue)] = {
            "residue": residue,
            "score": round(mean_score, 2),
            "class": _confidence_class(mean_score),
            "color": _confidence_color(mean_score),
        }

    return confidence, "alphafold_pdb_b_factor"


def _render_html_legacy(report: dict[str, Any], annotations: dict[str, Any], pdb_text: str) -> str:
    """Legacy renderer retained for compatibility; use _render_html for current viewers."""
    return _render_html(report, annotations, pdb_text)


def _render_html(report: dict[str, Any], annotations: dict[str, Any], pdb_text: str) -> str:
    """Render a standalone HTML viewer with data-driven 3Dmol.js controls."""
    gene = str(report.get("gene", "unknown"))
    structure_label = "AlphaFold PDB structure" if pdb_text != EXAMPLE_PDB else "embedded example protein structure"
    template = """<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>StructPhenotypes Viewer - __GENE__</title>
  <script src="https://3Dmol.org/build/3Dmol-min.js"></script>
  <style>
    :root {
      color-scheme: light;
      font-family: Arial, sans-serif;
      --border: #d8dee8;
      --ink: #1f2937;
      --muted: #607085;
      --panel: #f7f9fc;
      --blue: #2459d6;
    }
    body {
      margin: 0;
      color: var(--ink);
      background: #ffffff;
    }
    header {
      padding: 18px 22px 12px;
      border-bottom: 1px solid var(--border);
    }
    h1 {
      margin: 0 0 6px;
      font-size: 22px;
    }
    .subtitle {
      margin: 0;
      color: var(--muted);
      font-size: 14px;
    }
    .toolbar {
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(220px, 1fr));
      gap: 12px;
      align-items: end;
      padding: 14px 22px;
      background: var(--panel);
      border-bottom: 1px solid var(--border);
    }
    label {
      display: grid;
      gap: 5px;
      font-size: 12px;
      font-weight: 700;
      text-transform: uppercase;
      letter-spacing: 0.04em;
      color: var(--muted);
    }
      select,
      button,
      input[type="color"],
      input[type="range"],
      input[type="text"] {
        height: 36px;
        border: 1px solid var(--border);
        border-radius: 6px;
        background: #ffffff;
        color: var(--ink);
        font-size: 14px;
        padding: 0 10px;
      }
      input[type="range"] {
        padding: 0;
      }
    button {
      cursor: pointer;
      color: #ffffff;
      background: var(--blue);
      border-color: var(--blue);
      font-weight: 700;
    }
    #viewer {
      width: 100%;
      height: calc(100vh - 176px);
      min-height: 430px;
      position: relative;
    }
    #status {
      padding: 10px 12px;
      color: var(--muted);
      border-top: 1px solid var(--border);
      border-bottom: 1px solid var(--border);
      font-size: 14px;
    }
    .workspace {
      display: grid;
      grid-template-columns: minmax(0, 2fr) minmax(260px, 0.55fr) minmax(320px, 1fr);
      gap: 14px;
      padding: 16px 22px 24px;
      align-items: start;
    }
    .viewer-pane,
    .legend-pane,
    .control-pane {
      min-width: 0;
    }
    .legend-pane {
      position: sticky;
      top: 12px;
      max-height: calc(100vh - 24px);
      overflow: auto;
      border-left: 1px solid var(--border);
      border-right: 1px solid var(--border);
      padding: 0 12px;
    }
    .control-pane {
      display: grid;
      gap: 14px;
    }
    .panel {
      border: 1px solid var(--border);
      border-radius: 8px;
      padding: 12px;
      min-width: 0;
    }
    .panel h2 {
      margin: 0 0 8px;
      font-size: 15px;
    }
    .panel p,
    .panel li {
      color: var(--muted);
      font-size: 14px;
      line-height: 1.4;
    }
    .phenotype-panel {
      max-height: calc(100vh - 310px);
      overflow: auto;
      display: grid;
      gap: 7px;
      padding-right: 6px;
    }
    .query-panel {
      max-height: calc(100vh - 290px);
      overflow: auto;
      display: grid;
      gap: 8px;
      padding-right: 6px;
    }
    .query-row {
      border: 1px solid #edf1f7;
      border-radius: 8px;
      padding: 8px;
      display: grid;
      gap: 8px;
      background: #ffffff;
    }
    .query-row.is-placeholder {
      opacity: 0.58;
      background: #f8fafc;
    }
    .query-row.is-invalid {
      border-color: #dc2626;
      background: #fef2f2;
    }
    .query-row-top {
      display: grid;
      grid-template-columns: 20px 40px 84px minmax(0, 1fr) auto;
      gap: 8px;
      align-items: center;
    }
    .query-row-top input[type="checkbox"] {
      width: 15px;
      height: 15px;
    }
    .query-row-top input[type="color"] {
      width: 38px;
      min-width: 38px;
      height: 30px;
      padding: 2px;
    }
    .query-size {
      width: 100%;
      min-width: 0;
    }
    .query-input-wrap {
      position: relative;
      min-width: 0;
    }
    .query-input {
      width: 100%;
      font-family: Consolas, "Courier New", monospace;
      font-size: 13px;
    }
    .query-row-meta {
      display: flex;
      justify-content: space-between;
      gap: 8px;
      color: var(--muted);
      font-size: 12px;
      flex-wrap: wrap;
    }
    .query-error {
      color: #b91c1c;
      font-size: 12px;
    }
    .query-suggestions {
      position: absolute;
      left: 0;
      right: 0;
      top: calc(100% + 4px);
      z-index: 15;
      background: #ffffff;
      border: 1px solid var(--border);
      border-radius: 8px;
      box-shadow: 0 12px 24px rgba(15, 23, 42, 0.12);
      max-height: 180px;
      overflow: auto;
      padding: 4px 0;
    }
    .query-suggestion {
      padding: 6px 10px;
      font-size: 12px;
      color: var(--ink);
      cursor: pointer;
      font-family: Consolas, "Courier New", monospace;
    }
    .query-suggestion.is-active,
    .query-suggestion:hover {
      background: #eff6ff;
      color: #1d4ed8;
    }
    .query-examples {
      margin-top: 12px;
      display: grid;
      gap: 4px;
    }
    .query-examples code {
      display: block;
      padding: 6px 8px;
      border-radius: 6px;
      background: #f8fafc;
      border: 1px solid #edf1f7;
      font-size: 12px;
      overflow-wrap: anywhere;
    }
    .phenotype-row {
      display: grid;
      grid-template-columns: 22px 42px minmax(0, 1fr) auto;
      gap: 8px;
      align-items: center;
      border-bottom: 1px solid #edf1f7;
      padding: 4px 0;
    }
    .phenotype-row.has-options {
      grid-template-columns: 22px 42px minmax(0, 1fr) auto 28px;
    }
    .phenotype-row input[type="checkbox"] {
      width: 16px;
      height: 16px;
    }
    .phenotype-row input[type="color"] {
      width: 38px;
      min-width: 38px;
      padding: 2px;
    }
    .phenotype-name {
      color: var(--ink);
      font-size: 13px;
      overflow-wrap: anywhere;
    }
    .phenotype-count {
      color: var(--muted);
      font-size: 12px;
    }
      .phenotype-options {
        grid-column: 3 / -1;
        display: none;
        gap: 6px;
        padding: 2px 0 5px;
      }
      .phenotype-options.is-visible {
        display: grid;
      }
      .pathogenicity-row {
        display: grid;
        grid-template-columns: 20px minmax(0, 1fr) 40px auto;
        gap: 7px;
        align-items: center;
        color: var(--muted);
        font-size: 12px;
      }
      .pathogenicity-toggle {
        grid-template-columns: 20px 1fr;
        padding-bottom: 4px;
        border-bottom: 1px solid #edf1f7;
        margin-bottom: 4px;
      }
      .pathogenicity-row input[type="checkbox"] {
        width: 15px;
        height: 15px;
      }
      .pathogenicity-row input[type="color"] {
        width: 38px;
        min-width: 38px;
        height: 30px;
        padding: 2px;
      }
      .pathogenicity-count {
        font-size: 11px;
        color: var(--muted);
        white-space: nowrap;
      }
      .small-button {
        height: 28px;
        min-width: 28px;
        padding: 0;
        border-radius: 6px;
        background: #ffffff;
        border-color: var(--border);
        color: var(--muted);
      }
      .stats-grid {
        display: grid;
        grid-template-columns: repeat(auto-fit, minmax(150px, 1fr));
        gap: 8px;
        margin: 8px 0 14px;
      }
      .stat {
        border: 1px solid #edf1f7;
        border-radius: 6px;
        padding: 8px;
      }
      .stat strong {
        display: block;
        color: var(--ink);
        font-size: 18px;
      }
      .stat span {
        color: var(--muted);
        font-size: 12px;
      }
      .mini-list {
        margin: 8px 0 0;
        padding-left: 18px;
      }
    .legend {
      display: flex;
      flex-wrap: wrap;
      gap: 8px;
      margin-top: 8px;
    }
    .legend-item {
      display: inline-flex;
      gap: 6px;
      align-items: center;
      color: var(--muted);
      font-size: 13px;
    }
    .swatch {
      width: 14px;
      height: 14px;
      border-radius: 3px;
      border: 1px solid var(--border);
    }
    .colorbar {
      display: grid;
      gap: 6px;
      margin: 8px 0 14px;
    }
    .colorbar-track {
      height: 16px;
      border: 1px solid var(--border);
      border-radius: 999px;
      background: linear-gradient(90deg, #2b83ba, #ffffbf, #d7191c);
    }
    .colorbar-labels {
      display: flex;
      justify-content: space-between;
      color: var(--muted);
      font-size: 12px;
    }
      .legend-note {
        color: var(--muted);
        font-size: 13px;
        line-height: 1.4;
        margin: 0 0 8px;
      }
      .slider-value {
        color: var(--ink);
        font-size: 13px;
        font-weight: 400;
        text-transform: none;
        letter-spacing: 0;
      }
      .inline-field {
        display: grid;
        grid-template-columns: 1fr 84px;
        gap: 8px;
        align-items: center;
      }
      @media (max-width: 800px) {
        .workspace {
          grid-template-columns: 1fr;
        }
        .legend-pane {
          position: static;
          max-height: none;
          border-left: 0;
          border-right: 0;
          padding: 0;
        }
        #viewer {
          height: 58vh;
        }
      }
  </style>
</head>
<body>
  <header>
    <h1>StructPhenotypes Viewer: __GENE__</h1>
    <p class="subtitle">Local HTML viewer using 3Dmol.js and __STRUCTURE_LABEL__.</p>
  </header>

  <div class="toolbar">
    <label>
      Base color mode
      <select id="baseColorMode">
        <option value="gray">Gray protein</option>
        <option value="default">Default spectrum</option>
        <option value="missense-gradient">AlphaMissense score gradient</option>
        <option value="missense-class">AlphaMissense class</option>
        <option value="alphafold-confidence">AlphaFold confidence</option>
        <option value="clinvar-pathogenicity">ClinVar pathogenicity</option>
        <option value="surprise">Surprise: ClinVar vs AlphaMissense</option>
        <option value="contact-network">Disease contact network</option>
      </select>
    </label>
    <label>
      Rendering style
      <select id="renderStyle">
        <option value="cartoon">Cartoon</option>
        <option value="trace">Trace cartoon</option>
        <option value="tube">Tube cartoon</option>
        <option value="cartoon-stick">Cartoon + sticks</option>
        <option value="stick">Sticks</option>
        <option value="line">Lines</option>
        <option value="sphere">Sphere / spacefill</option>
        <option value="cross">Cross</option>
        <option value="surface-smooth">Smooth surface</option>
        <option value="surface-vdw">VDW surface</option>
        <option value="surface-ms">Molecular surface</option>
        <option value="surface-sas">Solvent-accessible surface</option>
        <option value="surface-ses">Solvent-excluded surface</option>
      </select>
    </label>
      <label>
        Focus residue
        <select id="focusResidue"></select>
      </label>
      <label>
        Background opacity
        <input id="backgroundOpacity" type="range" min="0.05" max="1" step="0.05" value="0.8">
        <span id="backgroundOpacityValue" class="slider-value">0.80</span>
      </label>
      <label>
        Contact cutoff
        <div class="inline-field">
          <input id="contactCutoff" type="range" min="4" max="12" step="0.5" value="8.0">
          <input id="contactCutoffText" type="text" value="8.0 A" inputmode="decimal" aria-label="Contact cutoff value">
        </div>
        <span id="contactCutoffValue" class="slider-value">8.0 A</span>
      </label>
      <button id="resetView" type="button">Reset view</button>
    </div>

  <section class="workspace">
    <div class="viewer-pane">
      <div id="viewer"></div>
      <div id="status">Loading viewer...</div>
    </div>
    <aside class="legend-pane">
      <h2>Color Scale</h2>
      <div id="colorLegend"></div>
    </aside>
    <div class="control-pane">
      <div class="panel">
        <h2>Query Highlights</h2>
        <p>Write residue queries to add colored markers without changing the protein coloring underneath.</p>
        <div id="queryControls" class="query-panel">
          <div class="query-row">
            <div class="query-row-top">
              <input class="query-enabled" type="checkbox" checked disabled>
              <input class="query-color" type="color" value="#1f77b4" disabled aria-label="Query color">
              <input class="query-size" type="number" min="0.6" max="50.0" step="0.1" value="5.0" disabled aria-label="Query marker size">
              <div class="query-input-wrap">
                <input class="query-input" type="text" spellcheck="false" placeholder="e.g. phenotype CONTAINS 'epileptic'" aria-label="Residue query">
              </div>
              <span class="pathogenicity-count">loading</span>
            </div>
          </div>
        </div>
        <div id="queryExamples" class="query-examples">
          <p class="legend-note">Queries are residue-level WHERE-style expressions. Strings should be quoted.</p>
          <p class="legend-note">Use <strong>AND</strong> for intersection, <strong>OR</strong> for union, and <strong>NOT</strong> for negation.</p>
          <p class="legend-note">Common fields: residue, phenotype, pathogenicity_class, function_class, missense_score, missense_class, alphafold_score, alphafold_class.</p>
          <code>all</code>
          <code>residue IN (42, 43, 44)</code>
          <code>residue BETWEEN 100 AND 120</code>
          <code>function_class = 'gain-of-function'</code>
          <code>NOT function_class IN ('gain-of-function', 'loss-of-function')</code>
          <code>phenotype CONTAINS 'epileptic' AND pathogenicity_class = 'pathogenic'</code>
        </div>
      </div>
      <div class="panel">
        <h2>Query Match Summary</h2>
        <div id="diseaseStats"></div>
        <h2>Residue Summary</h2>
        <ul id="variantList"></ul>
      </div>
    </div>
  </section>

  <script type="application/json" id="structpheno-annotations-data">__PREPROCESSED_JSON__</script>
  <script type="application/json" id="structpheno-pdb-data">__PDB_JSON__</script>

  <script type="text/plain" id="unused-generated-viewer-script">
    const REPORT = __REPORT_JSON__;
    const PREPROCESSED_ANNOTATIONS = __PREPROCESSED_JSON__;
    const PDB_TEXT = __PDB_JSON__;

    const MISSENSE_CLASS_COLORS = {
      benign: "#2ca25f",
      ambiguous: "#fee08b",
      pathogenic: "#d73027",
      stub: "#9ca3af",
    };
    const SURFACE_TYPES = {
      "surface-smooth": "SAS",
      "surface-vdw": "VDW",
      "surface-ms": "MS",
      "surface-sas": "SAS",
      "surface-ses": "SES",
    };
      const PATHOGENICITY_CLASSES = [
        { key: "pathogenic", label: "Pathogenic", color: "#b2182b" },
        { key: "likely-pathogenic", label: "Likely pathogenic", color: "#ef8a62" },
        { key: "uncertain", label: "Uncertain/conflicting", color: "#fddbc7" },
        { key: "benign", label: "Benign/likely benign", color: "#d1e5f0" },
      ];
      const CONTACT_NEIGHBORHOOD_LIMIT = 8;
      const CONTACT_CLUSTER_MIN_SIZE = 3;
      const QUERY_ROW_MIN_SIZE = 0.6;
      const QUERY_ROW_MAX_SIZE = 50.0;
      const QUERY_EXAMPLES = [
        "all",
        "residue = 42",
        "residue IN (42, 43, 44)",
        "residue BETWEEN 100 AND 120",
        "function_class = 'gain-of-function'",
        "function_class = 'gain-of-function' OR function_class = 'loss-of-function'",
        "NOT function_class IN ('gain-of-function', 'loss-of-function')",
        "phenotype CONTAINS 'epileptic encephalopathy' AND pathogenicity_class = 'pathogenic'",
        "missense_score >= 0.8 AND alphafold_class IN ('confident', 'very high')",
      ];
      let viewer = null;
      let model = null;
      let residueCenterCache = new Map();
      let contactNetworkLabels = [];
      let currentFocusResidue = null;
      let currentSceneKey = null;
      let queryRows = [];
      let queryRowIdCounter = 0;
      let autocompleteState = {
        rowId: null,
        items: [],
        activeIndex: 0,
        tokenStart: 0,
        tokenEnd: 0,
      };

      function init() {
        populateControls();
        populateVariantList();

        if (!window.$3Dmol) {
          setStatus("3Dmol.js did not load from the CDN. Check internet access or open this file in a browser that permits external scripts.");
          return;
        }

        viewer = $3Dmol.createViewer("viewer", { backgroundColor: "white" });
        model = viewer.addModel(PDB_TEXT, "pdb");
        if (!model && typeof viewer.getModel === "function") {
          model = viewer.getModel();
        }
        renderBasicProtein();
        viewer.zoomTo();
        viewer.render();

        try {
          applyStyle(true);
          viewer.zoomTo();
          viewer.render();
          setStatus(proteinLoadStatus("Protein loaded. Choose a base color mode, add query highlights, or focus a residue."));
        } catch (error) {
          const message = error && error.message ? error.message : String(error);
          setStatus(proteinLoadStatus(`Protein loaded with the basic style, but annotation styling failed: ${message}`));
          logError(error);
        }
      }

      function renderBasicProtein() {
        if (!viewer) {
          return;
        }
        viewer.setStyle({}, { cartoon: { color: "#8da0cb", opacity: selectedBackgroundOpacity() } });
      }

      function proteinLoadStatus(prefix) {
        const atomCount = model && typeof model.selectedAtoms === "function"
          ? (model.selectedAtoms({}) || []).length
          : null;
        if (atomCount === null) {
          return `${prefix} Structure text length: ${PDB_TEXT.length.toLocaleString()} characters.`;
        }
        return `${prefix} Parsed ${atomCount.toLocaleString()} atoms.`;
      }

      function populateControls() {
        populateQueryControls();
        populateResidueFocus();
        updateBackgroundOpacityLabel();
        updateContactCutoffLabel();

        document.getElementById("baseColorMode").addEventListener("change", () => applyStyle(false));
        document.getElementById("renderStyle").addEventListener("change", () => applyStyle(false));
        document.getElementById("focusResidue").addEventListener("change", () => handleFocusChange(true));
        document.getElementById("backgroundOpacity").addEventListener("input", () => {
          updateBackgroundOpacityLabel();
          applyStyle(false);
        });
        document.getElementById("contactCutoff").addEventListener("input", () => {
          updateContactCutoffLabel();
          applyStyle(false);
        });
        document.getElementById("contactCutoffText").addEventListener("change", () => {
          syncContactCutoffFromText();
          applyStyle(false);
        });
        document.getElementById("resetView").addEventListener("click", resetView);
      }

      function populateQueryControls() {
        initQueryRows();
        renderQueryControls();
      }

      function initQueryRows() {
        queryRows = [createQueryRow(false, ""), createQueryRow(true, "")];
        autocompleteState = {
          rowId: null,
          items: [],
          activeIndex: 0,
          tokenStart: 0,
          tokenEnd: 0,
        };
      }

      function createQueryRow(isPlaceholder, queryText) {
        const activeIndex = queryRows.filter((row) => !row.isPlaceholder).length;
        return {
          id: `query-row-${queryRowIdCounter++}`,
          enabled: !isPlaceholder,
          color: fallbackColor(activeIndex),
          size: 5.0,
          draftSize: "5.0",
          query: queryText,
          draftQuery: queryText,
          isPlaceholder,
          error: null,
          matchResidues: [],
          matchCount: 0,
        };
      }

      function renderQueryControls(focusRowId = null, caret = null) {
        ensureTrailingPlaceholderRow();
        evaluateQueryRows();

        const panel = document.getElementById("queryControls");
        const examples = document.getElementById("queryExamples");
        if (!panel || !examples) {
          return;
        }

        panel.innerHTML = queryRows.map((row) => {
          const countLabel = row.isPlaceholder ? "placeholder" : String(row.matchCount) + " hits";
          return `
          <div class="query-row${row.isPlaceholder ? " is-placeholder" : ""}${row.error ? " is-invalid" : ""}" data-row-id="${row.id}">
            <div class="query-row-top">
              <input class="query-enabled" type="checkbox" data-row-id="${row.id}" ${row.enabled ? "checked" : ""} ${row.isPlaceholder ? "disabled" : ""}>
              <input class="query-color" type="color" data-row-id="${row.id}" value="${row.color}" ${row.isPlaceholder ? "disabled" : ""} aria-label="Query color">
              <input class="query-size" type="number" min="${QUERY_ROW_MIN_SIZE}" max="${QUERY_ROW_MAX_SIZE}" step="0.1" value="${escapeAttr(row.draftSize ?? Number(row.size).toFixed(1))}" data-row-id="${row.id}" ${row.isPlaceholder ? "disabled" : ""} aria-label="Query marker size">
              <div class="query-input-wrap">
                <input class="query-input" type="text" spellcheck="false" placeholder="${row.isPlaceholder ? "Type a query to activate this row" : "e.g. phenotype CONTAINS 'epileptic'"}" value="${escapeAttr(row.draftQuery ?? row.query)}" data-row-id="${row.id}" aria-label="Residue query">
                ${renderAutocomplete(row.id)}
              </div>
              <span class="pathogenicity-count">${escapeHtml(countLabel)}</span>
            </div>
            <div class="query-row-meta">
              <span>${row.isPlaceholder ? "Placeholder row" : (row.enabled ? "Enabled" : "Disabled")}</span>
              ${row.error ? `<span class="query-error">${escapeHtml(row.error)}</span>` : ""}
            </div>
          </div>
        `;
        }).join("");

        examples.innerHTML = queryDocumentationHtml();

        panel.querySelectorAll(".query-enabled").forEach((input) => {
          input.addEventListener("change", handleQueryRowToggle);
        });
        panel.querySelectorAll(".query-color").forEach((input) => {
          input.addEventListener("input", handleQueryRowColorChange);
        });
        panel.querySelectorAll(".query-size").forEach((input) => {
          input.addEventListener("input", handleQueryRowSizeChange);
          input.addEventListener("keydown", handleQueryRowSizeKeydown);
        });
        panel.querySelectorAll(".query-input").forEach((input) => {
          input.addEventListener("input", handleQueryRowInput);
          input.addEventListener("keydown", handleQueryInputKeydown);
          input.addEventListener("blur", handleQueryInputBlur);
        });
        panel.querySelectorAll(".query-suggestion").forEach((item) => {
          item.addEventListener("mousedown", handleAutocompleteSelection);
        });

        if (focusRowId) {
          const input = panel.querySelector(`.query-input[data-row-id="${focusRowId}"]`);
          if (input) {
            input.focus();
            if (caret !== null) {
              const safeCaret = Math.max(0, Math.min(caret, input.value.length));
              input.setSelectionRange(safeCaret, safeCaret);
            }
          }
        }
      }

      function renderAutocomplete(rowId) {
        if (autocompleteState.rowId !== rowId || autocompleteState.items.length === 0) {
          return "";
        }
        return `
          <div class="query-suggestions">
            ${autocompleteState.items.map((item, index) => `
              <div class="query-suggestion${index === autocompleteState.activeIndex ? " is-active" : ""}" data-row-id="${rowId}" data-index="${index}">
                ${escapeHtml(item)}
              </div>
            `).join("")}
          </div>
        `;
      }

      function queryDocumentationHtml() {
        return `
          <p class="legend-note">Queries are residue-level WHERE-style expressions. They match residues, not individual variant records.</p>
          <p class="legend-note"><strong>Set logic:</strong> use <code>OR</code> for union, <code>AND</code> for intersection, <code>NOT</code> for negation, and parentheses for grouping.</p>
          <p class="legend-note"><strong>Fields:</strong> all, residue, variant_count, primary_phenotype, phenotype, primary_significance, pathogenicity_class, has_pathogenic, function_class, missense_score, missense_class, alphafold_score, alphafold_class.</p>
          <p class="legend-note"><strong>Operators:</strong> =, !=, &lt;, &lt;=, &gt;, &gt;=, IN (...), BETWEEN ... AND ..., CONTAINS.</p>
          <p class="legend-note"><strong>Common values:</strong> pathogenicity_class is 'pathogenic', 'likely-pathogenic', 'uncertain', or 'benign'. function_class is 'gain-of-function' or 'loss-of-function'.</p>
          <p class="legend-note">Examples:</p>
          ${QUERY_EXAMPLES.map((example) => `<code>${escapeHtml(example)}</code>`).join("")}
        `;
      }

      function ensureTrailingPlaceholderRow() {
        const nonPlaceholder = queryRows.filter((row) => !row.isPlaceholder);
        const placeholder = queryRows.filter((row) => row.isPlaceholder);
        if (placeholder.length === 0) {
          queryRows.push(createQueryRow(true, ""));
          return;
        }
        const trailingPlaceholder = placeholder[placeholder.length - 1];
        queryRows = nonPlaceholder.concat([trailingPlaceholder]);
      }

      function findQueryRow(rowId) {
        return queryRows.find((row) => row.id === rowId) || null;
      }

      function handleQueryRowToggle(event) {
        const row = findQueryRow(event.target.dataset.rowId);
        if (!row || row.isPlaceholder) {
          return;
        }
        row.enabled = Boolean(event.target.checked);
        applyStyle(false);
      }

      function handleQueryRowColorChange(event) {
        const row = findQueryRow(event.target.dataset.rowId);
        if (!row || row.isPlaceholder) {
          return;
        }
        row.color = event.target.value || row.color;
        applyStyle(false);
      }

      function handleQueryRowSizeChange(event) {
        const row = findQueryRow(event.target.dataset.rowId);
        if (!row || row.isPlaceholder) {
          return;
        }
        row.draftSize = event.target.value;
      }

      function handleQueryRowSizeKeydown(event) {
        if (event.key !== "Enter") {
          return;
        }
        event.preventDefault();
        const row = findQueryRow(event.target.dataset.rowId);
        if (!row || row.isPlaceholder) {
          return;
        }
        row.draftSize = event.target.value;
        row.size = normalizeQueryRowSize(event.target.value);
        row.draftSize = row.size.toFixed(1);
        event.target.value = row.draftSize;
        applyStyle(false);
      }

      function handleQueryRowInput(event) {
        const rowId = event.target.dataset.rowId;
        const row = findQueryRow(rowId);
        if (!row) {
          return;
        }
        row.draftQuery = event.target.value;
        updateAutocompleteForInput(event.target);
      }

      function commitQueryRowInput(input) {
        const rowId = input.dataset.rowId;
        const row = findQueryRow(rowId);
        if (!row) {
          return;
        }
        row.draftQuery = input.value;
        row.query = input.value;
        if (row.isPlaceholder && row.query.trim()) {
          row.isPlaceholder = false;
          row.enabled = true;
          queryRows.push(createQueryRow(true, ""));
        }
        renderQueryControls(rowId, inputCaret(input, row.query.length));
        applyStyle(false);
      }

      function handleQueryInputFocus(event) {
        updateAutocompleteForInput(event.target);
        renderQueryControls(event.target.dataset.rowId, inputCaret(event.target, event.target.value.length));
      }

      function handleQueryInputBlur() {
        window.setTimeout(() => {
          autocompleteState = {
            rowId: null,
            items: [],
            activeIndex: 0,
            tokenStart: 0,
            tokenEnd: 0,
          };
          renderQueryControls();
        }, 120);
      }

      function handleQueryInputKeydown(event) {
        const suggestionsActive = autocompleteState.rowId === event.target.dataset.rowId && autocompleteState.items.length > 0;
        if (event.key === "Enter") {
          event.preventDefault();
          commitQueryRowInput(event.target);
        } else if (suggestionsActive && event.key === "ArrowDown") {
          event.preventDefault();
          autocompleteState.activeIndex = (autocompleteState.activeIndex + 1) % autocompleteState.items.length;
          renderQueryControls(event.target.dataset.rowId, inputCaret(event.target, event.target.value.length));
        } else if (suggestionsActive && event.key === "ArrowUp") {
          event.preventDefault();
          autocompleteState.activeIndex = (autocompleteState.activeIndex - 1 + autocompleteState.items.length) % autocompleteState.items.length;
          renderQueryControls(event.target.dataset.rowId, inputCaret(event.target, event.target.value.length));
        } else if (suggestionsActive && event.key === "Tab") {
          event.preventDefault();
          applyAutocompleteChoice(event.target.dataset.rowId, autocompleteState.items[autocompleteState.activeIndex]);
        } else if (suggestionsActive && event.key === "Escape") {
          autocompleteState.items = [];
          autocompleteState.rowId = null;
          renderQueryControls(event.target.dataset.rowId, inputCaret(event.target, event.target.value.length));
        }
      }

      function inputCaret(input, fallback) {
        return typeof input.selectionStart === "number" ? input.selectionStart : fallback;
      }

      function handleAutocompleteSelection(event) {
        const rowId = event.currentTarget.dataset.rowId;
        const index = Number(event.currentTarget.dataset.index);
        if (!Number.isFinite(index) || !autocompleteState.items[index]) {
          return;
        }
        applyAutocompleteChoice(rowId, autocompleteState.items[index]);
      }

      function applyAutocompleteChoice(rowId, suggestion) {
        const row = findQueryRow(rowId);
        if (!row) {
          return;
        }
        const query = row.query;
        const tokenStart = autocompleteState.tokenStart;
        row.query = `${query.slice(0, autocompleteState.tokenStart)}${suggestion}${query.slice(autocompleteState.tokenEnd)}`;
        if (row.isPlaceholder && row.query.trim()) {
          row.isPlaceholder = false;
          row.enabled = true;
          queryRows.push(createQueryRow(true, ""));
        }
        autocompleteState = {
          rowId: rowId,
          items: [],
          activeIndex: 0,
          tokenStart: 0,
          tokenEnd: 0,
        };
        const nextCaret = tokenStart + suggestion.length;
        renderQueryControls(rowId, nextCaret);
        applyStyle(false);
      }

      function normalizeQueryRowSize(value) {
        const parsed = Number(value);
        if (!Number.isFinite(parsed)) {
          return 5.0;
        }
        return Math.max(QUERY_ROW_MIN_SIZE, Math.min(QUERY_ROW_MAX_SIZE, parsed));
      }

      function populateGlobalPathogenicityControls() {
        const panel = document.getElementById("globalPathogenicityControls");
        if (!panel) {
          return;
        }
        panel.innerHTML = PATHOGENICITY_CLASSES.map((pathClass) => `
          <label class="pathogenicity-row">
            <input class="global-pathogenicity-check" type="checkbox" data-class="${pathClass.key}" checked>
            <span>${escapeHtml(pathClass.label)}</span>
            <input class="global-pathogenicity-color" type="color" data-class="${pathClass.key}" value="${pathClass.color}" aria-label="${escapeHtml(pathClass.label)} color">
            <span class="pathogenicity-count">${pathogenicityClassResidueCount(pathClass.key)} residues</span>
          </label>
        `).join("");

        panel.querySelectorAll("input").forEach((input) => {
          input.addEventListener("change", handleGlobalPathogenicityChange);
        });
      }

      function populateFunctionControls() {
        const panel = document.getElementById("functionControls");
        if (!panel) {
          return;
        }
        const counts = PREPROCESSED_ANNOTATIONS.function_class_counts || {};
        const available = FUNCTION_CLASSES.filter((functionClass) => Number(counts[functionClass.key] || 0) > 0);
        const neitherCount = functionClassResidueCount("neither");
        if (available.length === 0) {
          panel.innerHTML = neitherCount > 0 ? `
            <label class="pathogenicity-row">
              <input class="function-check" type="checkbox" data-class="neither" checked>
              <span>Neither</span>
              <span class="pathogenicity-count">${neitherCount} residues</span>
            </label>
          ` : "<p>No experimental gain/loss-of-function calls were found for this gene.</p>";
        } else {
          panel.innerHTML = available.map((functionClass) => `
            <label class="pathogenicity-row">
              <input class="function-check" type="checkbox" data-class="${functionClass.key}" checked>
              <span>${escapeHtml(functionClass.label)}</span>
              <span class="pathogenicity-count">${functionClassResidueCount(functionClass.key)} residues</span>
            </label>
          `).join("") + `
            <label class="pathogenicity-row">
              <input class="function-check" type="checkbox" data-class="neither">
              <span>Neither</span>
              <span class="pathogenicity-count">${neitherCount} residues</span>
            </label>
          `;
        }
        panel.querySelectorAll("input").forEach((input) => {
          input.addEventListener("change", handleOverlayFilterChange);
        });
      }

      function populatePhenotypeControls() {
        const panel = document.getElementById("phenotypeControls");
      const counts = PREPROCESSED_ANNOTATIONS.phenotype_counts || {};
      const colors = PREPROCESSED_ANNOTATIONS.phenotype_colors || {};
      phenotypeEntries = Object.entries(counts)
        .sort((left, right) => right[1] - left[1])
        .slice(0, PHENOTYPE_LIMIT)
        .map(([name, count], index) => ({
          name,
          count,
          color: colors[name] || fallbackColor(index),
        }));

      if (phenotypeEntries.length === 0) {
        panel.innerHTML = "<p>No phenotype annotations found.</p>";
        return;
      }

        panel.innerHTML = `
          <label class="pathogenicity-row">
            <input id="selectAllPhenotypes" type="checkbox">
            <span>Select all phenotypes</span>
            <span class="pathogenicity-count">${phenotypeEntries.length} sets</span>
          </label>
        ` + phenotypeEntries.map((entry, index) => `
          <div class="phenotype-row has-options" data-index="${index}">
            <input class="phenotype-check" type="checkbox" data-index="${index}">
            <input class="phenotype-color" type="color" data-index="${index}" value="${entry.color}">
            <span class="phenotype-name">${escapeHtml(entry.name)}</span>
            <span class="phenotype-count">${entry.count}</span>
            <button class="small-button phenotype-options-toggle" type="button" data-index="${index}" aria-label="Show pathogenicity options for ${escapeHtml(entry.name)}">...</button>
            <div class="phenotype-options" data-index="${index}">
              ${phenotypePathogenicityOptions(index, entry.name)}
            </div>
        </div>
      `).join("");

        panel.querySelectorAll("input").forEach((input) => {
          if (input.classList.contains("phenotype-pathogenicity-color")) {
            input.addEventListener("input", () => {
              input.dataset.customized = "true";
            });
          }
          if (input.id === "selectAllPhenotypes") {
            return;
          }
          input.addEventListener("change", handleOverlayFilterChange);
        });
        panel.querySelectorAll(".phenotype-check").forEach((checkbox) => {
          checkbox.addEventListener("change", () => syncPhenotypeOptionsVisibility(Number(checkbox.dataset.index)));
        });
        panel.querySelectorAll(".phenotype-options-toggle").forEach((button) => {
          button.addEventListener("click", () => togglePhenotypeOptions(Number(button.dataset.index)));
        });
        const selectAll = document.getElementById("selectAllPhenotypes");
        if (selectAll) {
          selectAll.addEventListener("change", handleSelectAllPhenotypesChange);
        }
        syncSelectAllPhenotypesState();
      }

      function phenotypePathogenicityOptions(index, phenotype) {
        return `
          <label class="pathogenicity-row pathogenicity-toggle">
            <input class="phenotype-use-pathogenicity-colors" type="checkbox" data-index="${index}">
            <span>Use pathogenicity-specific colors for this phenotype</span>
          </label>
        ` + PATHOGENICITY_CLASSES.map((pathClass) => `
          <label class="pathogenicity-row">
            <input class="phenotype-pathogenicity-check" type="checkbox" data-index="${index}" data-class="${pathClass.key}" checked>
            <span>${escapeHtml(pathClass.label)}</span>
            <input class="phenotype-pathogenicity-color" type="color" data-index="${index}" data-class="${pathClass.key}" value="${pathClass.color}" aria-label="${escapeHtml(pathClass.label)} color for ${escapeHtml(phenotype)}">
            <span class="pathogenicity-count">${phenotypePathogenicityResidueCount(phenotype, pathClass.key)} residues</span>
          </label>
        `)
          .join("");
      }

      function syncPhenotypeOptionsVisibility(index) {
        const checkbox = document.querySelector(`.phenotype-check[data-index="${index}"]`);
        const options = document.querySelector(`.phenotype-options[data-index="${index}"]`);
        if (options && checkbox && checkbox.checked) {
          options.classList.add("is-visible");
        } else if (options) {
          options.classList.remove("is-visible");
        }
      }

      function togglePhenotypeOptions(index) {
        const options = document.querySelector(`.phenotype-options[data-index="${index}"]`);
        if (options) {
          options.classList.toggle("is-visible");
        }
      }

      function handleSelectAllPhenotypesChange(event) {
        const checked = Boolean(event && event.target && event.target.checked);
        document.querySelectorAll(".phenotype-check").forEach((checkbox) => {
          checkbox.checked = checked;
          syncPhenotypeOptionsVisibility(Number(checkbox.dataset.index));
          if (!checked) {
            const options = document.querySelector(`.phenotype-options[data-index="${checkbox.dataset.index}"]`);
            if (options) {
              options.classList.remove("is-visible");
            }
          }
        });
        syncSelectAllPhenotypesState();
        handleOverlayFilterChange();
      }

      function syncSelectAllPhenotypesState() {
        const selectAll = document.getElementById("selectAllPhenotypes");
        const checkboxes = Array.from(document.querySelectorAll(".phenotype-check"));
        if (!selectAll || checkboxes.length === 0) {
          return;
        }
        const checkedCount = checkboxes.filter((checkbox) => checkbox.checked).length;
        selectAll.checked = checkedCount > 0 && checkedCount === checkboxes.length;
        selectAll.indeterminate = checkedCount > 0 && checkedCount < checkboxes.length;
      }

    function populateResidueFocus() {
      const residueSelect = document.getElementById("focusResidue");
      const residues = getResidueAnnotations();
      residueSelect.innerHTML = '<option value="">No residue focus</option>' +
        residues
          .map((annotation) => `<option value="${annotation.residue}">Residue ${annotation.residue}</option>`)
          .join("");
    }

    function getResidueAnnotations() {
      return PREPROCESSED_ANNOTATIONS.residue_list || [];
    }

      function evaluateQueryRows() {
        queryRows.forEach((row) => {
          row.error = null;
          row.matchResidues = [];
          row.matchCount = 0;
          if (row.isPlaceholder || !row.query.trim()) {
            return;
          }
          try {
            const predicate = compileResidueQuery(row.query);
            const matches = getResidueAnnotations()
              .filter((annotation) => predicate(annotation))
              .map((annotation) => annotation.residue);
            row.matchResidues = matches;
            row.matchCount = matches.length;
          } catch (error) {
            row.error = error instanceof Error ? error.message : String(error);
          }
        });
      }

      function activeQueryRows() {
        evaluateQueryRows();
        return queryRows.filter((row) => !row.isPlaceholder && row.enabled && row.query.trim() && !row.error);
      }

      function filteredResidueAnnotations() {
        const uniqueResidues = new Set();
        activeQueryRows().forEach((row) => {
          row.matchResidues.forEach((residue) => uniqueResidues.add(residue));
        });
        return Array.from(uniqueResidues)
          .map((residue) => annotationForResidue(residue))
          .filter(Boolean)
          .sort((left, right) => left.residue - right.residue);
      }

      function queryMatchEntriesForResidue(annotation) {
        return activeQueryRows().filter((row) => row.matchResidues.includes(annotation.residue));
      }

      function queryMarkerRadius(annotation, rowSize) {
        const size = normalizeQueryRowSize(rowSize);
        const variantScale = annotation ? Math.min(1.0 + Math.log2((annotation.variant_count || 0) + 1) * 0.16, 1.8) : 1.0;
        return Math.max(0.18, size * 0.16 * variantScale);
      }

      function addQueryMarkers() {
        const rows = activeQueryRows();
        if (rows.length === 0) {
          return;
        }
        const uniqueResidues = filteredResidueAnnotations();
        uniqueResidues.forEach((annotation) => {
          const matchingRows = rows.filter((row) => row.matchResidues.includes(annotation.residue));
          if (matchingRows.length === 0) {
            return;
          }
          matchingRows.forEach((row, index) => {
            const radius = queryMarkerRadius(annotation, row.size);
            styleResidueMarker(
              annotation.residue,
              row.color,
              1.0,
              radius,
              phenotypeMarkerCenter(annotation.residue, `query:${row.id}`, index, matchingRows.length, radius)
            );
          });
        });
      }

      function querySummaryStats() {
        const rows = activeQueryRows();
        const unique = new Set();
        rows.forEach((row) => row.matchResidues.forEach((residue) => unique.add(residue)));
        return {
          activeRows: rows,
          uniqueResidueCount: unique.size,
        };
      }

      function queryFieldValue(annotation, field) {
        const normalizedField = String(field || "").toLowerCase();
        if (normalizedField === "residue") {
          return Number(annotation.residue);
        }
        if (normalizedField === "variant_count") {
          return Number(annotation.variant_count || 0);
        }
        if (normalizedField === "primary_phenotype") {
          return annotation.primary_phenotype || null;
        }
        if (normalizedField === "phenotype") {
          return annotation.phenotypes || [];
        }
        if (normalizedField === "primary_significance") {
          return annotation.primary_significance || null;
        }
        if (normalizedField === "pathogenicity_class") {
          return annotationPathogenicityClass(annotation);
        }
        if (normalizedField === "has_pathogenic") {
          return Boolean(annotation.has_pathogenic);
        }
        if (normalizedField === "function_class") {
          return annotation.function_classes || [];
        }
        if (normalizedField === "missense_score") {
          return annotation.missense ? Number(annotation.missense.score) : null;
        }
        if (normalizedField === "missense_class") {
          return annotation.missense ? annotation.missense.class : null;
        }
        if (normalizedField === "alphafold_score") {
          return annotation.alphafold_confidence ? Number(annotation.alphafold_confidence.score) : null;
        }
        if (normalizedField === "alphafold_class") {
          return annotation.alphafold_confidence ? annotation.alphafold_confidence.class : null;
        }
        throw new Error(`Unknown field '${field}'`);
      }

      function annotationPathogenicityClass(annotation) {
        if (!Number(annotation.variant_count || 0)) {
          return null;
        }
        return annotation.pathogenicity_class || significanceClassKey(annotation.primary_significance || "Benign");
      }

      function compileResidueQuery(queryText) {
        const tokens = tokenizeQuery(queryText);
        if (tokens.length === 0) {
          throw new Error("Enter a query.");
        }
        const parser = createQueryParser(tokens);
        const ast = parser.parseExpression();
        parser.expectEnd();
        return (annotation) => Boolean(evaluateQueryAst(ast, annotation));
      }

      function tokenizeQuery(queryText) {
        const tokens = [];
        let index = 0;
        while (index < queryText.length) {
          const char = queryText[index];
          if (/\s/.test(char)) {
            index += 1;
            continue;
          }
          if (char === "(" || char === ")" || char === ",") {
            tokens.push({ type: char, value: char });
            index += 1;
            continue;
          }
          if ((char === "!" || char === "<" || char === ">") && queryText[index + 1] === "=") {
            tokens.push({ type: "operator", value: `${char}=` });
            index += 2;
            continue;
          }
          if (char === "=" || char === "<" || char === ">") {
            tokens.push({ type: "operator", value: char });
            index += 1;
            continue;
          }
          if (char === "'" || char === "\"") {
            const quote = char;
            let cursor = index + 1;
            let value = "";
            while (cursor < queryText.length && queryText[cursor] !== quote) {
              value += queryText[cursor];
              cursor += 1;
            }
            if (cursor >= queryText.length) {
              throw new Error("Unterminated string literal.");
            }
            tokens.push({ type: "string", value });
            index = cursor + 1;
            continue;
          }
          const match = /^[A-Za-z_][A-Za-z0-9_-]*/.exec(queryText.slice(index));
          if (match) {
            const value = match[0];
            if (/^(true|false)$/i.test(value)) {
              tokens.push({ type: "boolean", value: value.toLowerCase() === "true" });
            } else {
              tokens.push({ type: "word", value });
            }
            index += value.length;
            continue;
          }
          const numberMatch = /^\d+(?:\.\d+)?/.exec(queryText.slice(index));
          if (numberMatch) {
            tokens.push({ type: "number", value: Number(numberMatch[0]) });
            index += numberMatch[0].length;
            continue;
          }
          throw new Error(`Unexpected token near '${queryText.slice(index, index + 12)}'.`);
        }
        return tokens;
      }

      function createQueryParser(tokens) {
        let cursor = 0;

        function peek(offset = 0) {
          return tokens[cursor + offset] || null;
        }

        function consume() {
          const token = tokens[cursor];
          cursor += 1;
          return token;
        }

        function matchesWord(word) {
          const token = peek();
          return token && token.type === "word" && token.value.toUpperCase() === word;
        }

        function parseExpression() {
          return parseOr();
        }

        function parseOr() {
          let node = parseAnd();
          while (matchesWord("OR")) {
            consume();
            node = { type: "or", left: node, right: parseAnd() };
          }
          return node;
        }

        function parseAnd() {
          let node = parseNot();
          while (matchesWord("AND")) {
            consume();
            node = { type: "and", left: node, right: parseNot() };
          }
          return node;
        }

        function parseNot() {
          if (matchesWord("NOT")) {
            consume();
            return { type: "not", value: parseNot() };
          }
          return parsePrimary();
        }

        function parsePrimary() {
          const token = peek();
          if (!token) {
            throw new Error("Unexpected end of query.");
          }
          if (token.type === "(") {
            consume();
            const expr = parseExpression();
            expect(")");
            return expr;
          }
          if (token.type === "word" && token.value.toLowerCase() === "all") {
            consume();
            return { type: "literal", value: true };
          }
          return parsePredicate();
        }

        function parsePredicate() {
          const fieldToken = consume();
          if (!fieldToken || fieldToken.type !== "word") {
            throw new Error("Expected a field name.");
          }
          const field = fieldToken.value;
          if (matchesWord("BETWEEN")) {
            consume();
            const start = parseValue();
            if (!matchesWord("AND")) {
              throw new Error("BETWEEN requires AND.");
            }
            consume();
            const end = parseValue();
            return { type: "between", field, start, end };
          }
          if (matchesWord("IN")) {
            consume();
            expect("(");
            const values = [];
            while (peek() && peek().type !== ")") {
              values.push(parseValue());
              if (peek() && peek().type === ",") {
                consume();
              } else {
                break;
              }
            }
            expect(")");
            return { type: "in", field, values };
          }
          if (matchesWord("CONTAINS")) {
            consume();
            return { type: "contains", field, value: parseValue() };
          }
          const operatorToken = consume();
          if (!operatorToken || operatorToken.type !== "operator") {
            throw new Error(`Expected an operator after '${field}'.`);
          }
          return { type: "compare", field, operator: operatorToken.value, value: parseValue() };
        }

        function parseValue() {
          const token = consume();
          if (!token) {
            throw new Error("Expected a value.");
          }
          if (token.type === "number" || token.type === "string" || token.type === "boolean") {
            return token.value;
          }
          if (token.type === "word") {
            return token.value;
          }
          throw new Error("Expected a literal value.");
        }

        function expect(type) {
          const token = consume();
          if (!token || token.type !== type) {
            throw new Error(`Expected '${type}'.`);
          }
        }

        function expectEnd() {
          if (cursor < tokens.length) {
            throw new Error(`Unexpected token '${tokens[cursor].value}'.`);
          }
        }

        return {
          parseExpression,
          expectEnd,
        };
      }

      function evaluateQueryAst(ast, annotation) {
        if (ast.type === "literal") {
          return ast.value;
        }
        if (ast.type === "or") {
          return evaluateQueryAst(ast.left, annotation) || evaluateQueryAst(ast.right, annotation);
        }
        if (ast.type === "and") {
          return evaluateQueryAst(ast.left, annotation) && evaluateQueryAst(ast.right, annotation);
        }
        if (ast.type === "not") {
          return !evaluateQueryAst(ast.value, annotation);
        }
        if (ast.type === "compare") {
          const fieldValue = queryFieldValue(annotation, ast.field);
          return compareFieldValue(fieldValue, ast.operator, ast.value);
        }
        if (ast.type === "contains") {
          return containsFieldValue(queryFieldValue(annotation, ast.field), ast.value);
        }
        if (ast.type === "in") {
          return inFieldValue(queryFieldValue(annotation, ast.field), ast.values);
        }
        if (ast.type === "between") {
          return betweenFieldValue(queryFieldValue(annotation, ast.field), ast.start, ast.end);
        }
        return false;
      }

      function compareFieldValue(fieldValue, operator, rawValue) {
        if (fieldValue === null || fieldValue === undefined) {
          return false;
        }
        if (Array.isArray(fieldValue)) {
          if (operator === "=") {
            return fieldValue.some((value) => valuesEqual(value, rawValue));
          }
          if (operator === "!=") {
            return fieldValue.every((value) => !valuesEqual(value, rawValue));
          }
          return false;
        }
        const left = normalizeComparable(fieldValue);
        const right = normalizeComparable(rawValue);
        if (operator === "=") {
          return valuesEqual(left, right);
        }
        if (operator === "!=") {
          return !valuesEqual(left, right);
        }
        if (typeof left === "number" && typeof right === "number") {
          if (operator === "<") {
            return left < right;
          }
          if (operator === "<=") {
            return left <= right;
          }
          if (operator === ">") {
            return left > right;
          }
          if (operator === ">=") {
            return left >= right;
          }
        }
        return false;
      }

      function containsFieldValue(fieldValue, rawValue) {
        const needle = String(rawValue || "").toLowerCase();
        if (!needle) {
          return false;
        }
        if (Array.isArray(fieldValue)) {
          return fieldValue.some((value) => String(value || "").toLowerCase().includes(needle));
        }
        return String(fieldValue || "").toLowerCase().includes(needle);
      }

      function inFieldValue(fieldValue, values) {
        if (Array.isArray(fieldValue)) {
          return fieldValue.some((value) => values.some((candidate) => valuesEqual(value, candidate)));
        }
        return values.some((candidate) => valuesEqual(fieldValue, candidate));
      }

      function betweenFieldValue(fieldValue, start, end) {
        if (fieldValue === null || fieldValue === undefined || Array.isArray(fieldValue)) {
          return false;
        }
        const value = normalizeComparable(fieldValue);
        const min = normalizeComparable(start);
        const max = normalizeComparable(end);
        if (typeof value === "number" && typeof min === "number" && typeof max === "number") {
          return value >= min && value <= max;
        }
        if (typeof value === "string" && typeof min === "string" && typeof max === "string") {
          return value >= min && value <= max;
        }
        return false;
      }

      function valuesEqual(left, right) {
        if (typeof left === "number" || typeof right === "number") {
          return Number(left) === Number(right);
        }
        if (typeof left === "boolean" || typeof right === "boolean") {
          return Boolean(left) === Boolean(right);
        }
        return String(left || "").toLowerCase() === String(right || "").toLowerCase();
      }

      function normalizeComparable(value) {
        if (typeof value === "number" || typeof value === "boolean") {
          return value;
        }
        if (value === null || value === undefined) {
          return value;
        }
        if (/^-?\d+(?:\.\d+)?$/.test(String(value))) {
          return Number(value);
        }
        if (/^(true|false)$/i.test(String(value))) {
          return String(value).toLowerCase() === "true";
        }
        return String(value).toLowerCase();
      }

      function queryAutocompleteValues() {
        const phenotypeValues = Object.keys(PREPROCESSED_ANNOTATIONS.phenotype_counts || {}).map((value) => `'${value}'`);
        const significanceValues = Array.from(new Set(
          getResidueAnnotations()
            .map((annotation) => annotation.primary_significance)
            .filter(Boolean)
        )).map((value) => `'${value}'`);
        const functionValues = ["'gain-of-function'", "'loss-of-function'"];
        const missenseClassValues = ["'benign'", "'ambiguous'", "'pathogenic'", "'stub'"];
        const alphafoldClassValues = ["'very low'", "'low'", "'confident'", "'very high'"];
        return [
          "all",
          "residue",
          "variant_count",
          "primary_phenotype",
          "phenotype",
          "primary_significance",
          "pathogenicity_class",
          "has_pathogenic",
          "function_class",
          "missense_score",
          "missense_class",
          "alphafold_score",
          "alphafold_class",
          "=",
          "!=",
          "<",
          "<=",
          ">",
          ">=",
          "AND",
          "OR",
          "NOT",
          "IN",
          "BETWEEN",
          "CONTAINS",
          "true",
          "false",
          "'pathogenic'",
          "'likely-pathogenic'",
          "'uncertain'",
          "'benign'",
          ...functionValues,
          ...missenseClassValues,
          ...alphafoldClassValues,
          ...significanceValues,
          ...phenotypeValues,
        ];
      }

      function updateAutocompleteForInput(input) {
        const rowId = input.dataset.rowId;
        const cursor = inputCaret(input, input.value.length);
        const token = currentQueryToken(input.value, cursor);
        const prefix = token.text.toLowerCase();
        const suggestions = queryAutocompleteValues()
          .filter((candidate, index, values) => values.indexOf(candidate) === index)
          .filter((candidate) => !prefix || candidate.toLowerCase().includes(prefix))
          .slice(0, 12);
        autocompleteState = {
          rowId,
          items: suggestions,
          activeIndex: 0,
          tokenStart: token.start,
          tokenEnd: token.end,
        };
      }

      function currentQueryToken(query, cursor) {
        let start = cursor;
        let end = cursor;
        while (start > 0 && /[A-Za-z0-9_'"-]/.test(query[start - 1])) {
          start -= 1;
        }
        while (end < query.length && /[A-Za-z0-9_'"-]/.test(query[end])) {
          end += 1;
        }
        return {
          text: query.slice(start, end),
          start,
          end,
        };
      }

      function selectedPhenotypes() {
        return new Map();
      }

      function selectedPhenotypeClasses(index) {
        return new Map();
      }

      function availableFunctionClassKeys() {
        return [];
      }

      function selectedFunctionClasses() {
        return new Set();
      }

      function functionFilterIsActive() {
        return false;
      }

      function variantMatchesFunctionFilter(variant) {
        return true;
      }

      function filteredVariants(annotation) {
        return (annotation.variants || []).filter((variant) => variantMatchesFunctionFilter(variant));
      }

      function annotationMatchesFunctionFilter(annotation) {
        return true;
      }

      function annotationHasKnownFunction(annotation) {
        if ((annotation.function_classes || []).length > 0) {
          return true;
        }
        return (annotation.variants || []).some((variant) => (variant.function_classes || []).length > 0);
      }

      function annotationMatchesSelectedPhenotypes(annotation, phenotypes = selectedPhenotypes()) {
        return true;
      }

      function annotationMatchesActiveFilters(annotation, phenotypes = selectedPhenotypes()) {
        return true;
      }

    function populateVariantList() {
      const list = document.getElementById("variantList");
      const residues = filteredResidueAnnotations();
      const displayed = residues.slice(0, 120);
      if (displayed.length === 0) {
        list.innerHTML = "<li>No residues match the active queries.</li>";
        return;
      }
      list.innerHTML = displayed
        .map((annotation) => `<li>Residue ${annotation.residue}: ${annotation.primary_significance || "annotated"}, ${annotation.primary_phenotype || "no phenotype"} (${queryMatchEntriesForResidue(annotation).length} query hits)</li>`)
        .join("");
    }

      function populateDiseaseStats() {
        const panel = document.getElementById("diseaseStats");
        if (!panel) {
          return;
        }
        const summary = querySummaryStats();
        const rows = summary.activeRows;
        const mode = document.getElementById("baseColorMode").value;
        const network = buildContactNetwork(selectedContactCutoff());
        const mainNeighborhoods = summarizeContactNeighborhoods(network);
        panel.innerHTML = `
          <div class="stats-grid">
            ${statItem(rows.length, "active query rows")}
            ${statItem(summary.uniqueResidueCount, "unique matched residues")}
            ${statItem(getResidueAnnotations().length, "total residues")}
            ${statItem(queryRows.filter((row) => !row.isPlaceholder && row.query.trim() && row.error).length, "invalid queries")}
          </div>
          <p class="legend-note">${mode === "contact-network"
            ? `Contact network at ${selectedContactCutoff().toFixed(1)} A: showing ${mainNeighborhoods.length} main clusters of ${network.components.length} total, largest cluster ${network.largestClusterSize} residues.`
            : "Colored query markers are overlaid on top of the protein and do not change the protein color mode."}</p>
          <ol class="mini-list">
            ${rows.map((row) => `<li><span class="swatch" style="background:${escapeAttr(row.color)}"></span> ${escapeHtml(row.query)}: ${row.matchCount} residues</li>`).join("") || "<li>No active query rows yet.</li>"}
          </ol>
        `;
      }

      function statItem(value, label) {
        return `<div class="stat"><strong>${value}</strong><span>${escapeHtml(label)}</span></div>`;
      }

      function applyStyle(zoomFocus) {
        if (!viewer) {
          return;
        }

        const focusResidue = selectedFocusResidue();
        const restOpacity = selectedBackgroundOpacity();
        clearSurfaces();
        clearShapes();
        viewer.setStyle({}, {});

        renderScene(restOpacity);
        currentFocusResidue = focusResidue;
        currentSceneKey = sceneKey(restOpacity);
        redrawShapeOverlays();
        applyFocus(focusResidue, zoomFocus);
        populateDiseaseStats();
        populateVariantList();
        updateLegend();
        viewer.render();
      }

      function handleOverlayFilterChange() {
        if (!viewer) {
          return;
        }
        syncSelectAllPhenotypesState();
        applyStyle(false);
      }

      function handleFocusChange(zoomFocus) {
        if (!viewer) {
          return;
        }
        const nextFocusResidue = selectedFocusResidue();
        const hadFocus = Boolean(currentFocusResidue);
        const willHaveFocus = Boolean(nextFocusResidue);

        if (hadFocus !== willHaveFocus) {
          applyStyle(zoomFocus);
          return;
        }

        currentFocusResidue = nextFocusResidue;
        clearShapes();
        redrawShapeOverlays();
        applyFocus(nextFocusResidue, zoomFocus);
        viewer.render();
      }

      function renderScene(opacity) {
        applyBaseColorMode(opacity);
      }

      function redrawShapeOverlays() {
        const opacity = selectedBackgroundOpacity();
        if (isSurfaceStyle()) {
          addSurfaceModeMarkers(opacity);
        }
        addQueryMarkers();
      }

      function sceneKey(opacity) {
        return JSON.stringify({
          renderStyle: document.getElementById("renderStyle").value,
          baseColorMode: document.getElementById("baseColorMode").value,
          opacity: Number(opacity.toFixed(2)),
          contactCutoff: Number(selectedContactCutoff().toFixed(1)),
        });
      }

    function applyBaseColorMode(opacity) {
      const mode = document.getElementById("baseColorMode").value;
      const residues = getResidueAnnotations();

      if (mode === "gray") {
        styleWholeProtein("#d9d9d9", opacity);
        setStatus("Neutral gray protein.");
        return;
      }

      if (mode === "default") {
        if (isSurfaceStyle()) {
          styleWholeProtein("#d9d9d9", opacity);
          setStatus("Surface styles use a neutral protein surface. Query markers are drawn on top.");
          return;
        }
        styleWholeProtein("spectrum", opacity);
        setStatus("Default spectrum coloring.");
        return;
      }

      styleWholeProtein("#d9d9d9", opacity);

      if (mode === "missense-gradient") {
        residues.forEach((annotation) => {
          const color = annotation.missense
            ? (annotation.missense.color || scoreToColor(annotation.missense.score))
            : "#4b5563";
          styleResidue(annotation.residue, color, opacity, 0.1);
        });
        setStatus("Coloring by AlphaMissense mean pathogenicity score from 0 to 1.");
      } else if (mode === "missense-class") {
        residues.forEach((annotation) => {
          const color = annotation.missense
            ? missenseClassColor(annotation.missense.class)
            : "#4b5563";
          styleResidue(annotation.residue, color, opacity, 0.1);
        });
        setStatus("Coloring by AlphaMissense class.");
      } else if (mode === "alphafold-confidence") {
        residues.forEach((annotation) => {
          const color = annotation.alphafold_confidence
            ? annotation.alphafold_confidence.color
            : "#4b5563";
          styleResidue(annotation.residue, color, opacity, 0.08);
        });
        setStatus("Coloring by AlphaFold model confidence from the PDB B-factor column.");
      } else if (mode === "clinvar-pathogenicity") {
        residues.forEach((annotation) => {
          const pathClass = bestIntrinsicPathogenicityClass(annotation);
          const color = pathClass ? pathogenicityClassColor(pathClass) : "#4b5563";
          styleResidue(annotation.residue, color, opacity, 0.14);
        });
        setStatus("Coloring ClinVar residues by strongest pathogenicity label.");
      } else if (mode === "surprise") {
        residues.forEach((annotation) => {
          const surprise = surpriseClass(annotation);
          const color = surprise ? surprise.color : "#4b5563";
          styleResidue(annotation.residue, color, opacity, 0.16);
        });
        setStatus("Coloring residues where ClinVar labels and AlphaMissense scores disagree.");
      } else if (mode === "contact-network") {
        styleWholeProtein("#d9d9d9", opacity);
        addContactNetwork();
        setStatus(`Showing disease neighborhoods with an ${selectedContactCutoff().toFixed(1)} A cutoff. Query markers remain overlaid.`);
      }
    }

      function addPhenotypeMarkers(phenotypes) {
        const activeGlobalClasses = selectedGlobalPathogenicityClasses();
        getResidueAnnotations().forEach((annotation) => {
          if (!annotationMatchesFunctionFilter(annotation)) {
            return;
          }
          const markerItems = [];
          (annotation.phenotypes || []).forEach((phenotype) => {
            if (!phenotypes.has(phenotype)) {
              return;
            }
            const phenotypeConfig = phenotypes.get(phenotype);
            PATHOGENICITY_CLASSES.forEach((pathClass) => {
              const globalConfig = activeGlobalClasses.get(pathClass.key);
              const phenotypeClassConfig = phenotypeConfig.classes.get(pathClass.key);
              if (!globalConfig || !globalConfig.enabled || !phenotypeClassConfig || !phenotypeClassConfig.enabled) {
                return;
              }
              if (!phenotypeHasPathogenicityClass(annotation, phenotype, pathClass.key)) {
                return;
              }
              markerItems.push({
                key: `${phenotype}:${pathClass.key}`,
                color: phenotypeConfig.usePathogenicityColors
                  ? (phenotypeClassConfig.color || globalConfig.color)
                  : phenotypeConfig.color,
              });
            });
          });
          if (markerItems.length === 0) {
            return;
          }
          const radius = markerRadius(annotation, 0.38);
          markerItems.forEach((item, index) => {
            styleResidueMarker(
              annotation.residue,
              item.color,
              1.0,
              radius,
              phenotypeMarkerCenter(annotation.residue, item.key, index, markerItems.length, radius)
            );
          });
        });
      }

    function styleWholeProtein(color, opacity) {
      if (isSurfaceStyle()) {
        const surfaceColor = color === "spectrum" ? "#d9d9d9" : color;
        const showOutline = document.getElementById("renderStyle").value !== "surface-smooth";
        viewer.setStyle(
          {},
          showOutline ? { line: { color: "#b8bec8", opacity: Math.min(opacity, 0.08) } } : {}
        );
        addProteinSurface(surfaceColor, surfaceOpacity(opacity));
        return;
      }
      viewer.setStyle({}, renderStyle(color, opacity));
    }

      function styleResidue(residue, color, opacity, radius) {
        const annotation = annotationForResidue(residue);
        if (isSurfaceStyle()) {
          styleResidueMarker(
            residue,
            color,
            1.0,
            markerRadius(annotation, smoothSurfaceStyle() ? Math.max(radius * 4.2, 0.52) : Math.max(radius * 3.1, 0.4))
          );
          return;
        }
        viewer.setStyle(
          { resi: residue },
          renderStyle(color, opacity, radius)
        );
      }

      function styleResidueMarker(residue, color, opacity = 1.0, radius = 0.36, centerOverride = null) {
        const center = centerOverride || residueCenter(residue);
        if (!center) {
          return;
        }
        viewer.addSphere({
          center,
        radius,
        color,
        opacity,
      });
      }

      function applyFocus(residue, zoomFocus) {
        if (!residue) {
          return;
        }
        const annotation = annotationForResidue(residue);
        styleResidueMarker(residue, "#2459d6", 1.0, markerRadius(annotation, 0.62));
        if (zoomFocus) {
          focusCameraOnResidue(residue);
        }
      }

    function renderStyle(color, opacity, radius = 0.12) {
      const style = document.getElementById("renderStyle").value;
      if (style === "trace") {
        return { cartoon: { color, opacity, style: "trace" } };
      }
      if (style === "tube") {
        return { cartoon: { color, opacity, tubes: true, thickness: 0.6 } };
      }
      if (style === "stick") {
        return { stick: { color, radius, opacity } };
      }
      if (style === "line") {
        return { line: { color, opacity } };
      }
      if (style === "sphere") {
        return { sphere: { color, radius: Math.max(radius * 2.8, 0.35), opacity } };
      }
      if (style === "cross") {
        return { cross: { color, radius: Math.max(radius * 1.8, 0.22), opacity } };
      }
      if (style === "cartoon-stick") {
        return {
          cartoon: { color, opacity },
          stick: { color, radius, opacity },
        };
      }
      return { cartoon: { color, opacity } };
    }

      function isSurfaceStyle() {
        return document.getElementById("renderStyle").value.startsWith("surface-");
      }

      function smoothSurfaceStyle() {
        return document.getElementById("renderStyle").value === "surface-smooth";
      }

    function addProteinSurface(color, opacity) {
      const style = document.getElementById("renderStyle").value;
      const surfaceName = SURFACE_TYPES[style] || "VDW";
      const surfaceType = $3Dmol.SurfaceType[surfaceName] || $3Dmol.SurfaceType.VDW;
      const surface = viewer.addSurface(
        surfaceType,
        {
          color,
          opacity,
          transparent: opacity < 1.0,
        },
        {}
      );
      if (surface && typeof surface.then === "function") {
        surface.then(() => viewer.render());
      }
    }

      function surfaceOpacity(opacity) {
        const style = document.getElementById("renderStyle").value;
        if (style === "surface-smooth") {
          return Math.min(Math.max(opacity, 0.38), 0.72);
        }
        return Math.min(opacity, 0.88);
      }

    function clearSurfaces() {
      if (viewer && typeof viewer.removeAllSurfaces === "function") {
        viewer.removeAllSurfaces();
      }
    }

      function clearShapes() {
        if (viewer && typeof viewer.removeAllShapes === "function") {
          viewer.removeAllShapes();
        }
        clearContactNetworkLabels();
      }

      function clearContactNetworkLabels() {
        if (!viewer) {
          contactNetworkLabels = [];
          return;
        }
        if (typeof viewer.removeLabel === "function") {
          contactNetworkLabels.forEach((label) => {
            viewer.removeLabel(label);
          });
        } else if (typeof viewer.removeAllLabels === "function") {
          viewer.removeAllLabels();
        }
        contactNetworkLabels = [];
      }

      function addContactNetwork() {
        const network = buildContactNetwork(selectedContactCutoff());
        const neighborhoods = summarizeContactNeighborhoods(network);
        if (neighborhoods.length === 0) {
          return;
        }
        neighborhoods.forEach((cluster, index) => {
          viewer.addSphere({
            center: cluster.center,
            color: cluster.color,
            opacity: 0.2,
            radius: contactNeighborhoodRadius(cluster.size),
          });
          viewer.addSphere({
            center: cluster.center,
            color: cluster.color,
            opacity: 0.92,
            radius: contactNodeRadius(cluster.size),
          });
          if (typeof viewer.addLabel === "function") {
            const label = viewer.addLabel(`${index + 1}: ${cluster.size}`, {
              position: cluster.center,
              backgroundColor: "rgba(255,255,255,0.72)",
              fontColor: "#1f2937",
              borderColor: "#d8dee8",
              borderThickness: 1,
              inFront: true,
              showBackground: true,
            });
            if (label) {
              contactNetworkLabels.push(label);
            }
          }
        });
      }

      function addSurfaceModeMarkers(opacity) {
        const mode = document.getElementById("baseColorMode").value;
        if (mode === "gray" || mode === "default" || mode === "contact-network") {
          return;
        }
        getResidueAnnotations().forEach((annotation) => {
          if (mode === "missense-gradient" && annotation.missense) {
            styleResidue(annotation.residue, annotation.missense.color || scoreToColor(annotation.missense.score), opacity, 0.1);
          } else if (mode === "missense-class" && annotation.missense) {
            styleResidue(annotation.residue, missenseClassColor(annotation.missense.class), opacity, 0.1);
          } else if (mode === "alphafold-confidence" && annotation.alphafold_confidence) {
            styleResidue(annotation.residue, annotation.alphafold_confidence.color, opacity, 0.08);
          } else if (mode === "clinvar-pathogenicity") {
            const pathClass = bestIntrinsicPathogenicityClass(annotation);
            if (pathClass) {
              styleResidue(annotation.residue, pathogenicityClassColor(pathClass), opacity, 0.14);
            }
          } else if (mode === "surprise") {
            const surprise = surpriseClass(annotation);
            if (surprise) {
              styleResidue(annotation.residue, surprise.color, opacity, 0.16);
            }
          }
        });
      }

      function updateLegend() {
        const legend = document.getElementById("colorLegend");
        const mode = document.getElementById("baseColorMode").value;
        const sections = [];

        sections.push(baseLegend(mode));
        const activeRows = activeQueryRows();
        if (activeRows.length > 0) {
          sections.push(`
            <p class="legend-note">Active query overlays:</p>
            <div class="legend">
              ${activeRows.map((row) => legendItem(row.color, `${row.query} (${row.matchCount})`)).join("")}
            </div>
          `);
        } else {
          sections.push(`<p class="legend-note">No active query overlays. Type a query in the right-hand panel to add colored residue markers.</p>`);
        }
        legend.innerHTML = sections.join("");
      }

      function baseLegend(mode) {
        if (mode === "gray") {
          return `
            <p class="legend-note">Neutral gray protein with no annotation coloring applied.</p>
            <div class="legend">${legendItem("#d9d9d9", "Protein")}</div>
          `;
        }
        if (mode === "missense-gradient") {
          return `
            <p class="legend-note">AlphaMissense mean pathogenicity score.</p>
            <div class="colorbar">
              <div class="colorbar-track" style="background: linear-gradient(90deg, #2b83ba, #ffffbf, #d7191c);"></div>
              <div class="colorbar-labels"><span>0 benign-like</span><span>0.5</span><span>1 pathogenic-like</span></div>
            </div>
          `;
        }
        if (mode === "missense-class") {
          return `
            <p class="legend-note">AlphaMissense classification by residue mean score.</p>
            <div class="legend">
              ${legendItem(MISSENSE_CLASS_COLORS.benign, "Benign")}
              ${legendItem(MISSENSE_CLASS_COLORS.ambiguous, "Ambiguous")}
              ${legendItem(MISSENSE_CLASS_COLORS.pathogenic, "Pathogenic")}
            </div>
          `;
        }
        if (mode === "alphafold-confidence") {
          return `
            <p class="legend-note">AlphaFold confidence parsed from PDB B-factor/pLDDT values.</p>
            <div class="colorbar">
              <div class="colorbar-track" style="background: linear-gradient(90deg, #d62728, #c65d5b, #8e6e83, #1f77b4);"></div>
              <div class="colorbar-labels"><span>0 very low</span><span>50 low</span><span>70 confident</span><span>100 very high</span></div>
            </div>
          `;
        }
        if (mode === "clinvar-pathogenicity") {
          return `
            <p class="legend-note">Strongest ClinVar pathogenicity label at each residue.</p>
            ${pathogenicityClassLegend()}
          `;
        }
        if (mode === "contact-network") {
          return `
            <p class="legend-note">Disease contact network across all ClinVar-annotated residues. Inner sphere size scales with cluster size, and query markers remain visible on top.</p>
            ${pathogenicityClassLegend()}
          `;
        }
        if (mode === "surprise") {
          return `
            <p class="legend-note">Surprise layer: residues where ClinVar labels and AlphaMissense mean scores point in different directions.</p>
            <div class="legend">
              ${legendItem("#7b3294", "ClinVar pathogenic, low AlphaMissense")}
              ${legendItem("#008837", "ClinVar benign, high AlphaMissense")}
              ${legendItem("#fdae61", "ClinVar uncertain, high AlphaMissense")}
            </div>
          `;
        }
        return `
            <p class="legend-note">Default 3Dmol spectrum coloring. Query residue markers can be overlaid on top.</p>
            <div class="colorbar">
              <div class="colorbar-track" style="background: linear-gradient(90deg, #304ffe, #00bcd4, #4caf50, #ffeb3b, #f44336);"></div>
              <div class="colorbar-labels"><span>N terminus</span><span>C terminus</span></div>
            </div>
          `;
      }

      function legendItem(color, label) {
      return `<span class="legend-item"><span class="swatch" style="background:${color}"></span>${escapeHtml(label)}</span>`;
      }

      function pathogenicityColorLegend(note) {
        return `
          <p class="legend-note">${escapeHtml(note)}.</p>
          ${pathogenicityClassLegend()}
        `;
      }

      function pathogenicityClassLegend() {
        const items = PATHOGENICITY_CLASSES
          .map((pathClass) => legendItem(pathogenicityClassColor(pathClass.key), pathClass.label))
          .join("");
        return `<div class="legend">${items || '<span class="legend-note">No pathogenicity classes selected.</span>'}</div>`;
      }

      function bestIntrinsicPathogenicityClass(annotation) {
        const labels = Object.keys(annotation.significance_counts || {});
        return labels
          .map((label) => significanceClassKey(label))
          .sort((left, right) => pathogenicityClassRank(right) - pathogenicityClassRank(left))[0] || null;
      }

      function pathogenicityClassColor(pathClassKey) {
        const pathClass = PATHOGENICITY_CLASSES.find((item) => item.key === pathClassKey);
        return pathClass ? pathClass.color : "#9ca3af";
      }

      function selectedGlobalPathogenicityClasses() {
        const classes = new Map();
        PATHOGENICITY_CLASSES.forEach((pathClass) => {
          const checkbox = document.querySelector(`.global-pathogenicity-check[data-class="${pathClass.key}"]`);
          const colorInput = document.querySelector(`.global-pathogenicity-color[data-class="${pathClass.key}"]`);
          classes.set(pathClass.key, {
            enabled: Boolean(checkbox && checkbox.checked),
            color: colorInput ? colorInput.value || pathClass.color : pathClass.color,
          });
        });
        return classes;
      }

      function selectedGlobalPathogenicityColor(pathClassKey) {
        const selectedClass = selectedGlobalPathogenicityClasses().get(pathClassKey);
        const defaultClass = pathogenicityClassDefaults().get(pathClassKey);
        return (selectedClass && selectedClass.color) || (defaultClass && defaultClass.color) || "#9ca3af";
      }

      function pathogenicityClassDefaults() {
        return new Map(PATHOGENICITY_CLASSES.map((pathClass) => [pathClass.key, pathClass]));
      }

      function handleGlobalPathogenicityChange(event) {
        if (event && event.target && event.target.classList && event.target.classList.contains("global-pathogenicity-color")) {
          const pathClassKey = event.target.dataset.class;
          document.querySelectorAll(`.phenotype-pathogenicity-color[data-class="${pathClassKey}"]`).forEach((input) => {
            if (!input.dataset.customized) {
              input.value = event.target.value;
            }
          });
        }
        handleOverlayFilterChange();
      }

      function phenotypeVariantSignificances(annotation, phenotype) {
        const variants = filteredVariants(annotation);
        const labels = variants
          .filter((variant) => (variant.phenotype_names || []).includes(phenotype))
          .map((variant) => variant.significance)
          .filter(Boolean);
        if (labels.length > 0) {
          return labels;
        }
        return Object.keys(annotation.significance_counts || {});
      }

      function phenotypeHasPathogenicityClass(annotation, phenotype, pathClassKey) {
        return phenotypeVariantSignificances(annotation, phenotype).some((label) => significanceClassKey(label) === pathClassKey);
      }

      function bestActivePathogenicityClass(annotation) {
        const activeClasses = selectedGlobalPathogenicityClasses();
        const labels = filteredVariants(annotation)
          .map((variant) => variant.significance)
          .filter(Boolean);
        return labels
          .map((label) => significanceClassKey(label))
          .filter((pathClassKey) => {
            const activeClass = activeClasses.get(pathClassKey);
            return activeClass && activeClass.enabled;
          })
          .sort((left, right) => pathogenicityClassRank(right) - pathogenicityClassRank(left))[0] || null;
      }

      function pathogenicityClassResidueCount(pathClassKey) {
        return getResidueAnnotations()
          .filter((annotation) => filteredVariants(annotation).some((variant) => significanceClassKey(variant.significance) === pathClassKey))
          .length;
      }

      function phenotypePathogenicityResidueCount(phenotype, pathClassKey) {
        return getResidueAnnotations()
          .filter((annotation) => (annotation.phenotypes || []).includes(phenotype) && phenotypeHasPathogenicityClass(annotation, phenotype, pathClassKey))
          .length;
      }

      function functionClassResidueCount(functionClassKey) {
        if (functionClassKey === "neither") {
          return getResidueAnnotations()
            .filter((annotation) => !annotationHasKnownFunction(annotation))
            .length;
        }
        return getResidueAnnotations()
          .filter((annotation) => (annotation.function_classes || []).includes(functionClassKey))
          .length;
      }

      function significanceClassKey(significance) {
        const severity = significanceSeverity(significance);
        if (severity >= 4) {
          return "pathogenic";
        }
        if (severity === 3) {
          return "likely-pathogenic";
        }
        if (severity === 2) {
          return "uncertain";
        }
        return "benign";
      }

      function pathogenicityClassRank(pathClassKey) {
        if (pathClassKey === "pathogenic") {
          return 4;
        }
        if (pathClassKey === "likely-pathogenic") {
          return 3;
        }
        if (pathClassKey === "uncertain") {
          return 2;
        }
        if (pathClassKey === "benign") {
          return 1;
        }
        return 0;
      }

      function buildContactNetwork(cutoff) {
        const nodes = getResidueAnnotations()
          .map((annotation) => {
            const dominantClass = bestIntrinsicPathogenicityClass(annotation);
            const center = residueCenter(annotation.residue);
            if (!dominantClass || !center) {
              return null;
            }
            return {
              residue: annotation.residue,
              annotation,
              center,
              classKeys: [dominantClass],
              dominantClass,
              color: pathogenicityClassColor(dominantClass),
              clusterSize: 1,
              degree: 0,
            };
          })
          .filter(Boolean);

        const adjacency = new Map(nodes.map((node) => [node.residue, new Set()]));
        const edges = [];
        for (let leftIndex = 0; leftIndex < nodes.length; leftIndex += 1) {
          for (let rightIndex = leftIndex + 1; rightIndex < nodes.length; rightIndex += 1) {
            const leftNode = nodes[leftIndex];
            const rightNode = nodes[rightIndex];
            const distance = centerDistance(leftNode.center, rightNode.center);
            if (distance > cutoff) {
              continue;
            }
            adjacency.get(leftNode.residue).add(rightNode.residue);
            adjacency.get(rightNode.residue).add(leftNode.residue);
            edges.push({
              start: leftNode.center,
              end: rightNode.center,
              residues: [leftNode.residue, rightNode.residue],
              distance,
              clusterSize: 1,
            });
          }
        }

        const components = [];
        const seen = new Set();
        nodes.forEach((node) => {
          if (seen.has(node.residue)) {
            return;
          }
          const stack = [node.residue];
          const componentResidues = [];
          seen.add(node.residue);
          while (stack.length > 0) {
            const residue = stack.pop();
            componentResidues.push(residue);
            adjacency.get(residue).forEach((neighbor) => {
              if (!seen.has(neighbor)) {
                seen.add(neighbor);
                stack.push(neighbor);
              }
            });
          }
          components.push(componentResidues);
        });

        const componentByResidue = new Map();
        components.forEach((componentResidues, componentIndex) => {
          componentResidues.forEach((residue) => {
            componentByResidue.set(residue, { index: componentIndex, size: componentResidues.length });
          });
        });

        nodes.forEach((node) => {
          const component = componentByResidue.get(node.residue);
          const neighbors = adjacency.get(node.residue);
          node.clusterSize = component ? component.size : 1;
          node.degree = neighbors ? neighbors.size : 0;
        });
        edges.forEach((edge) => {
          const component = componentByResidue.get(edge.residues[0]);
          edge.clusterSize = component ? component.size : 1;
        });

        return {
          nodes,
          edges,
          components,
          largestClusterSize: components.reduce((maximum, component) => Math.max(maximum, component.length), 0),
        };
      }

      function summarizeContactNeighborhoods(network) {
        const nodesByResidue = new Map(network.nodes.map((node) => [node.residue, node]));
        return network.components
          .map((componentResidues) => {
            const clusterNodes = componentResidues
              .map((residue) => nodesByResidue.get(residue))
              .filter(Boolean);
            if (clusterNodes.length === 0) {
              return null;
            }
            const center = averageCenters(clusterNodes.map((node) => node.center));
            const classCounts = new Map();
            clusterNodes.forEach((node) => {
              classCounts.set(node.dominantClass, (classCounts.get(node.dominantClass) || 0) + 1);
            });
            const sortedClassCounts = Array.from(classCounts.entries())
              .sort((left, right) => right[1] - left[1] || pathogenicityClassRank(right[0]) - pathogenicityClassRank(left[0]));
            const dominantClass = sortedClassCounts.length > 0 ? sortedClassCounts[0][0] : "uncertain";
            return {
              size: clusterNodes.length,
              center,
              color: pathogenicityClassColor(dominantClass),
              dominantClass,
              residues: componentResidues,
            };
          })
          .filter((cluster) => cluster && cluster.size >= CONTACT_CLUSTER_MIN_SIZE && cluster.center)
          .sort((left, right) => right.size - left.size)
          .slice(0, CONTACT_NEIGHBORHOOD_LIMIT);
      }

      function averageCenters(centers) {
        if (!centers || centers.length === 0) {
          return null;
        }
        let x = 0;
        let y = 0;
        let z = 0;
        centers.forEach((center) => {
          x += center.x;
          y += center.y;
          z += center.z;
        });
        return {
          x: x / centers.length,
          y: y / centers.length,
          z: z / centers.length,
        };
      }

      function centerDistance(left, right) {
        const dx = left.x - right.x;
        const dy = left.y - right.y;
        const dz = left.z - right.z;
        return Math.sqrt((dx * dx) + (dy * dy) + (dz * dz));
      }

      function significanceSeverity(significance) {
        const normalized = String(significance || "").toLowerCase();
        if (normalized.includes("pathogenic/likely pathogenic")) {
          return 4;
        }
        if (normalized.includes("pathogenic") && !normalized.includes("likely")) {
          return 4;
        }
        if (normalized.includes("likely pathogenic")) {
          return 3;
        }
        if (normalized.includes("uncertain") || normalized.includes("conflicting")) {
          return 2;
        }
        if (normalized.includes("likely benign")) {
          return 1;
        }
        if (normalized.includes("benign")) {
          return 0;
        }
        return 0;
      }

      function surpriseClass(annotation) {
        if (!annotation.missense) {
          return null;
        }
        const score = Number(annotation.missense.score);
        if (!Number.isFinite(score)) {
          return null;
        }
        const severity = Number(annotation.max_pathogenicity) || 0;
        if (severity >= 3 && score <= 0.35) {
          return { label: "ClinVar pathogenic, low AlphaMissense", color: "#7b3294" };
        }
        if (severity <= 1 && annotation.variant_count > 0 && score >= 0.75) {
          return { label: "ClinVar benign, high AlphaMissense", color: "#008837" };
        }
        if (severity === 2 && score >= 0.75) {
          return { label: "ClinVar uncertain, high AlphaMissense", color: "#fdae61" };
        }
        return null;
      }

    function selectedFocusResidue() {
      const value = document.getElementById("focusResidue").value;
      return value ? Number(value) : null;
    }

    function selectedBackgroundOpacity() {
      const value = Number(document.getElementById("backgroundOpacity").value);
      return Math.max(0.05, Math.min(1.0, Number.isFinite(value) ? value : 0.3));
    }

    function selectedMarkerScale() {
      const input = document.getElementById("markerScale");
      const value = Number(input ? input.value : null);
      return Math.max(0.5, Math.min(50.0, Number.isFinite(value) ? value : 5.0));
    }

    function selectedContactCutoff() {
      const value = Number(document.getElementById("contactCutoff").value);
      return Math.max(4.0, Math.min(12.0, Number.isFinite(value) ? value : 8.0));
    }

    function updateBackgroundOpacityLabel() {
      document.getElementById("backgroundOpacityValue").textContent = selectedBackgroundOpacity().toFixed(2);
    }

      function updateMarkerScaleLabel() {
        const valueLabel = document.getElementById("markerScaleValue");
        const valueInput = document.getElementById("markerScaleText");
        if (!valueLabel || !valueInput) {
          return;
        }
        const label = `${selectedMarkerScale().toFixed(1)}x`;
        valueLabel.textContent = label;
        valueInput.value = label;
      }

      function updateContactCutoffLabel() {
        const label = `${selectedContactCutoff().toFixed(1)} A`;
        document.getElementById("contactCutoffValue").textContent = label;
        document.getElementById("contactCutoffText").value = label;
      }

      function syncMarkerScaleFromText() {
        const input = document.getElementById("markerScaleText");
        const slider = document.getElementById("markerScale");
        if (!input || !slider) {
          return;
        }
        const parsed = Number.parseFloat(String(input.value).replace(/x/gi, "").trim());
        const normalized = Math.max(0.5, Math.min(50.0, Number.isFinite(parsed) ? parsed : selectedMarkerScale()));
        slider.value = normalized.toFixed(1);
        updateMarkerScaleLabel();
      }

      function syncContactCutoffFromText() {
        const input = document.getElementById("contactCutoffText");
        const parsed = Number.parseFloat(String(input.value).replace(/a/gi, "").trim());
        const normalized = Math.max(4.0, Math.min(12.0, Number.isFinite(parsed) ? parsed : selectedContactCutoff()));
        document.getElementById("contactCutoff").value = normalized.toFixed(1);
        updateContactCutoffLabel();
      }

    function resetView() {
      if (!viewer) {
        return;
      }
        document.getElementById("focusResidue").value = "";
        document.getElementById("backgroundOpacity").value = "0.3";
        updateBackgroundOpacityLabel();
      initQueryRows();
      renderQueryControls();
      document.getElementById("contactCutoff").value = "8.0";
      updateContactCutoffLabel();
      applyStyle(false);
      viewer.zoomTo();
      viewer.render();
    }

      function residueCenter(residue) {
        if (residueCenterCache.has(residue)) {
          return residueCenterCache.get(residue);
        }
        if (!model || typeof model.selectedAtoms !== "function") {
          return null;
        }
        const atoms = model.selectedAtoms({ resi: residue }) || [];
        if (atoms.length === 0) {
          residueCenterCache.set(residue, null);
          return null;
        }
        const alphaCarbon = atoms.find((atom) => String(atom.atom || "").toUpperCase() === "CA");
        if (alphaCarbon) {
          const center = {
            x: Number(alphaCarbon.x) || 0,
            y: Number(alphaCarbon.y) || 0,
            z: Number(alphaCarbon.z) || 0,
          };
          residueCenterCache.set(residue, center);
          return center;
        }
        let x = 0;
        let y = 0;
        let z = 0;
      atoms.forEach((atom) => {
        x += Number(atom.x) || 0;
        y += Number(atom.y) || 0;
        z += Number(atom.z) || 0;
      });
        const count = atoms.length;
        const center = { x: x / count, y: y / count, z: z / count };
        residueCenterCache.set(residue, center);
        return center;
      }

      function focusCameraOnResidue(residue) {
        if (typeof viewer.center === "function") {
          viewer.center({ resi: residue }, 350);
          if (typeof viewer.zoom === "function") {
            viewer.zoom(1.08, 250);
          }
          return;
        }
        viewer.zoomTo({ resi: residue });
        if (typeof viewer.zoom === "function") {
          viewer.zoom(1.55, 250);
        }
      }

      function phenotypeMarkerCenter(residue, phenotype, index, total, markerRadiusValue) {
        const center = residueCenter(residue);
        if (!center || total <= 1) {
          return center;
        }
        const ringRadius = Math.max(markerRadiusValue * 0.9, 0.22);
        const angle = phenotypeAngle(residue, phenotype, index, total);
        return {
          x: center.x + ringRadius * Math.cos(angle),
          y: center.y + ringRadius * Math.sin(angle),
          z: center.z,
        };
      }

      function phenotypeAngle(residue, phenotype, index, total) {
        const seed = stableHash(`${residue}:${phenotype}`);
        const baseAngle = (seed % 360) * (Math.PI / 180);
        return baseAngle + (2 * Math.PI * index) / total;
      }

      function stableHash(value) {
        let hash = 2166136261;
        const text = String(value);
        for (let index = 0; index < text.length; index += 1) {
          hash ^= text.charCodeAt(index);
          hash = Math.imul(hash, 16777619);
        }
        return hash >>> 0;
      }

    function annotationForResidue(residue) {
      return getResidueAnnotations().find((annotation) => Number(annotation.residue) === Number(residue)) || null;
    }

      function markerRadius(annotation, baseRadius = 0.36) {
        const variantCount = Math.max(1, Number(annotation ? annotation.variant_count : 1) || 1);
        const scale = selectedMarkerScale();
        const variantBoost = 1 + Math.log2(variantCount) * 0.25;
        return Math.max(0.16, baseRadius * variantBoost * scale);
      }

      function contactNodeRadius(clusterSize) {
        const scale = selectedMarkerScale();
        return Math.max(0.2, Math.min(1.5, (0.16 + (0.09 * Math.sqrt(Math.max(1, clusterSize)))) * scale));
      }

      function contactNeighborhoodRadius(clusterSize) {
        const scale = selectedMarkerScale();
        return Math.max(0.45, Math.min(3.8, (0.42 + (0.22 * Math.sqrt(Math.max(1, clusterSize)))) * scale));
      }

      function significanceColor(significance) {
        const normalized = String(significance || "").toLowerCase();
        if (normalized.includes("pathogenic") && !normalized.includes("likely")) {
          return "#b2182b";
        }
        if (normalized.includes("likely pathogenic")) {
          return "#ef8a62";
        }
        if (normalized.includes("uncertain") || normalized.includes("conflicting")) {
          return "#fddbc7";
        }
        if (normalized.includes("likely benign")) {
          return "#d1e5f0";
        }
        if (normalized.includes("benign")) {
          return "#d9d9d9";
        }
        return "#9ca3af";
      }

      function missenseClassColor(value) {
        const normalized = String(value || "").toLowerCase();
        return MISSENSE_CLASS_COLORS[normalized] || "#9ca3af";
      }

    function scoreToColor(score) {
      const clamped = Math.max(0, Math.min(1, Number(score) || 0));
      const red = Math.round(43 + (215 - 43) * clamped);
      const green = Math.round(131 + (25 - 131) * clamped);
      const blue = Math.round(186 + (28 - 186) * clamped);
      return `rgb(${red}, ${green}, ${blue})`;
    }

    function fallbackColor(index) {
      const palette = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b"];
      return palette[index % palette.length];
    }

    function escapeHtml(value) {
      return String(value == null ? "" : value)
        .replace(/&/g, "&amp;")
        .replace(/</g, "&lt;")
        .replace(/>/g, "&gt;")
        .replace(/"/g, "&quot;");
    }

    function escapeAttr(value) {
      return escapeHtml(value).replace(/'/g, "&#39;");
    }

    function setStatus(message) {
      document.getElementById("status").textContent = message;
    }

    function logError(error) {
      if (window.console && typeof window.console.error === "function") {
        window.console.error(error);
      }
    }

    window.addEventListener("load", () => {
      try {
        init();
      } catch (error) {
        const message = error && error.message ? error.message : String(error);
        setStatus(`Viewer startup error: ${message}`);
        logError(error);
      }
    });
  </script>
  <script>
    (function(){
      "use strict";
      var ann=JSON.parse(document.getElementById("structpheno-annotations-data").textContent||"{}");
      var pdb=JSON.parse(document.getElementById("structpheno-pdb-data").textContent||'""');
      var viewer=null, model=null, rows=[], rowId=0, centers={}, palette=["#1f77b4","#ff7f0e","#2ca02c","#d62728","#9467bd","#8c564b"];
      function id(x){return document.getElementById(x);} function rs(){return ann.residue_list||[];} function stat(m){id("status").textContent=m;} function esc(v){return String(v==null?"":v).replace(/&/g,"&amp;").replace(/</g,"&lt;").replace(/>/g,"&gt;").replace(/"/g,"&quot;");} function attr(v){return esc(v).replace(/'/g,"&#39;");}
      function init(){focusList(); rows=[row(false),row(true)]; docs(); renderRows(); bindMain(); if(!window.$3Dmol){stat("3Dmol.js did not load."); return;} viewer=$3Dmol.createViewer("viewer",{backgroundColor:"white"}); model=viewer.addModel(pdb,"pdb"); draw(); viewer.zoomTo(); viewer.render(); stat("Protein loaded. Query highlights add colored residue markers.");}
      function bindMain(){var a=["baseColorMode","renderStyle","focusResidue","backgroundOpacity"],i,e; for(i=0;i<a.length;i++){e=id(a[i]); if(e){e.addEventListener("change",draw); e.addEventListener("input",draw);}} if(id("resetView")){id("resetView").addEventListener("click",function(){id("focusResidue").value=""; if(viewer){viewer.zoomTo();} draw();});}}
      function focusList(){var h='<option value="">No residue focus</option>',r=rs(),i; for(i=0;i<r.length;i++){h+='<option value="'+r[i].residue+'">Residue '+r[i].residue+'</option>';} id("focusResidue").innerHTML=h;}
      function op(){var v=Number(id("backgroundOpacity").value); if(!isFinite(v)){v=.8;} v=Math.max(.05,Math.min(1,v)); if(id("backgroundOpacityValue")){id("backgroundOpacityValue").textContent=v.toFixed(2);} return v;}
      function row(ph){var n=0,i; for(i=0;i<rows.length;i++){if(!rows[i].ph){n++;}} return {id:"q"+(++rowId),ph:ph,en:!ph,color:palette[n%palette.length],size:5.0,sd:"5.0",q:"",draft:"",err:"",m:[]};}
      function placeholder(){var a=[],p=null,i; for(i=0;i<rows.length;i++){if(rows[i].ph){p=rows[i];}else{a.push(rows[i]);}} if(!p){p=row(true);} rows=a.concat([p]);}
      function find(rid){var i; for(i=0;i<rows.length;i++){if(rows[i].id===rid){return rows[i];}} return null;}
      function renderRows(focusRid,caret){placeholder(); evalRows(); var h="",i,r,inp,pos; for(i=0;i<rows.length;i++){r=rows[i]; h+='<div class="query-row'+(r.ph?' is-placeholder':'')+(r.err?' is-invalid':'')+'"><div class="query-row-top">'; h+='<input class="query-enabled" data-r="'+r.id+'" type="checkbox"'+(r.en?' checked':'')+(r.ph?' disabled':'')+'>'; h+='<input class="query-color" data-r="'+r.id+'" type="color" value="'+attr(r.color)+'"'+(r.ph?' disabled':'')+'>'; h+='<input class="query-size" data-r="'+r.id+'" type="number" min="0.6" max="50.0" step="0.1" value="'+attr(r.sd)+'"'+(r.ph?' disabled':'')+'>'; h+='<div class="query-input-wrap"><input class="query-input" data-r="'+r.id+'" type="text" value="'+attr(r.draft)+'" placeholder="'+(r.ph?'Type a query to activate this row':"e.g. phenotype CONTAINS 'epileptic'")+'"></div>'; h+='<span class="pathogenicity-count">'+(r.ph?'placeholder':r.m.length+' hits')+'</span></div><div class="query-row-meta"><span>'+(r.ph?'Placeholder row':(r.en?'Enabled':'Disabled'))+'</span>'+(r.err?'<span class="query-error">'+esc(r.err)+'</span>':'')+'</div></div>'; } id("queryControls").innerHTML=h; bindRows(); summary(); list(); if(focusRid){inp=document.querySelector('.query-input[data-r="'+focusRid+'"]'); if(inp){inp.focus(); if(typeof caret==="number"&&inp.setSelectionRange){pos=Math.max(0,Math.min(caret,inp.value.length)); inp.setSelectionRange(pos,pos);}}}}
      function qsa(s,e,fn){var x=document.querySelectorAll(s),i; for(i=0;i<x.length;i++){x[i].addEventListener(e,function(ev){fn(ev.target,ev);});}}
      function bindRows(){qsa(".query-enabled","change",function(e){var r=find(e.getAttribute("data-r")); if(r){r.en=e.checked; draw();}}); qsa(".query-color","input",function(e){var r=find(e.getAttribute("data-r")); if(r){r.color=e.value; draw();}}); qsa(".query-size","input",function(e){var r=find(e.getAttribute("data-r")); if(r){r.sd=e.value;}}); qsa(".query-size","keydown",function(e,ev){var r=find(e.getAttribute("data-r")),v; if(!r){return;} if(ev&&ev.key==="Enter"){if(ev.preventDefault){ev.preventDefault();} v=Math.max(.6,Math.min(50.0,Number(e.value)||5.0)); r.size=v; r.sd=v.toFixed(1); e.value=r.sd; draw();}}); qsa(".query-input","input",function(e){var r=find(e.getAttribute("data-r")); if(r){r.draft=e.value;}}); qsa(".query-input","keydown",function(e,ev){var r=find(e.getAttribute("data-r")); if(!r){return;} if(ev&&ev.key==="Enter"){if(ev.preventDefault){ev.preventDefault();} r.draft=e.value; r.q=e.value; if(r.ph&&r.q.replace(/\s/g,"")){r.ph=false; r.en=true; rows.push(row(true));} renderRows(r.id,typeof e.selectionStart==="number"?e.selectionStart:e.value.length); draw();}});}
      function style(c,o,rad){var s=id("renderStyle").value; rad=rad||.12; if(s==="stick"){return {stick:{color:c,radius:rad,opacity:o}};} if(s==="line"){return {line:{color:c,opacity:o}};} if(s==="sphere"){return {sphere:{color:c,radius:Math.max(rad*2.8,.35),opacity:o}};} if(s==="cartoon-stick"){return {cartoon:{color:c,opacity:o},stick:{color:c,radius:rad,opacity:o}};} return {cartoon:{color:c,opacity:o}};}
      function isSurface(){return String(id("renderStyle").value).indexOf("surface-")===0;}
      function surfaceType(){var s=id("renderStyle").value,t="VDW"; if(s==="surface-ms"){t="MS";} if(s==="surface-sas"||s==="surface-smooth"){t="SAS";} if(s==="surface-ses"){t="SES";} if(window.$3Dmol&&$3Dmol.SurfaceType&&$3Dmol.SurfaceType[t]){return $3Dmol.SurfaceType[t];} return window.$3Dmol&&$3Dmol.SurfaceType?$3Dmol.SurfaceType.VDW:null;}
      function addSurface(c,o){var st=surfaceType(),smooth=id("renderStyle").value==="surface-smooth",alpha;if(!st){return;} alpha=smooth?Math.max(.25,Math.min(.85,o)):Math.max(.08,Math.min(.9,o)); viewer.addSurface(st,{color:c,opacity:alpha,transparent:alpha<1},{}); if(!smooth){viewer.setStyle({},{line:{color:"#aeb8c2",opacity:Math.min(.14,o)}});}}
      function paintResidue(res,c,o){if(isSurface()){sphere(res,c,.34,0,1);}else{viewer.setStyle({resi:res},style(c,o,.12));}}
      function draw(){evalRows(); summary(); list(); if(!viewer){return;} var mode=id("baseColorMode").value,o=op(),r=rs(),i,c; if(viewer.removeAllShapes){viewer.removeAllShapes();} if(viewer.removeAllSurfaces){viewer.removeAllSurfaces();} viewer.setStyle({},{}); if(isSurface()){addSurface("#d9d9d9",o);}else{viewer.setStyle({},style(mode==="default"?"spectrum":"#d9d9d9",o,.12));} for(i=0;i<r.length;i++){c=null; if(mode==="missense-gradient"&&r[i].missense){c=r[i].missense.color||score(r[i].missense.score);} if(mode==="missense-class"&&r[i].missense){c=miss(r[i].missense["class"]);} if(mode==="alphafold-confidence"&&r[i].alphafold_confidence){c=r[i].alphafold_confidence.color;} if(mode==="clinvar-pathogenicity"){c=pathColor(path(r[i]));} if(mode==="surprise"){c=surprise(r[i]);} if(c){paintResidue(r[i].residue,c,o);}} markers(); focus(); viewer.render();}
      function center(res){var k=String(res),a,i,c; if(centers.hasOwnProperty(k)){return centers[k];} if(!model||!model.selectedAtoms){return null;} a=model.selectedAtoms({resi:res})||[]; if(!a.length){centers[k]=null; return null;} for(i=0;i<a.length;i++){if(String(a[i].atom||"").toUpperCase()==="CA"){c={x:+a[i].x||0,y:+a[i].y||0,z:+a[i].z||0}; centers[k]=c; return c;}} c={x:0,y:0,z:0}; for(i=0;i<a.length;i++){c.x+=+a[i].x||0; c.y+=+a[i].y||0; c.z+=+a[i].z||0;} c.x/=a.length; c.y/=a.length; c.z/=a.length; centers[k]=c; return c;}
      function sphere(res,c,rad,idx,total){var p=center(res),ang,ring; if(!p){return;} if(total>1){ang=2*Math.PI*idx/total; ring=Math.max(rad*.9,.22); p={x:p.x+ring*Math.cos(ang),y:p.y+ring*Math.sin(ang),z:p.z};} viewer.addSphere({center:p,color:c,radius:rad,opacity:1});}
      function markers(){var g={},i,j,r,res; for(i=0;i<rows.length;i++){r=rows[i]; if(!r.en||r.ph||r.err){continue;} for(j=0;j<r.m.length;j++){res=r.m[j]; g[res]=g[res]||[]; g[res].push(r);}} for(res in g){if(g.hasOwnProperty(res)){for(i=0;i<g[res].length;i++){sphere(+res,g[res][i].color,Math.max(.18,g[res][i].size*.16),i,g[res].length);}}}}
      function focus(){var res=Number(id("focusResidue").value); if(res){sphere(res,"#2459d6",.62,0,1); viewer.zoomTo({resi:res});}}
      function evalRows(){var r=rs(),i,j; for(i=0;i<rows.length;i++){rows[i].m=[]; rows[i].err=""; if(rows[i].ph||!rows[i].q.replace(/\s/g,"")){continue;} try{for(j=0;j<r.length;j++){if(match(rows[i].q,r[j])){rows[i].m.push(r[j].residue);}}}catch(e){rows[i].err=e.message||String(e);}}}
      function match(q,a){q=q.replace(/^\s+|\s+$/g,""); if(!q||/^all$/i.test(q)){return true;} return expr(q,a);}
      function split(q,w){var out=[],start=0,i=0,j,ch,quote="",depth=0,word,upper=String(w).toUpperCase(),between=false; while(i<q.length){ch=q.charAt(i); if(quote){if(ch===quote){quote="";} i++; continue;} if(ch==="'"||ch==='"'){quote=ch; i++; continue;} if(ch==="("){depth++; i++; continue;} if(ch===")"){depth=Math.max(0,depth-1); i++; continue;} if(depth===0&&/[A-Za-z_]/.test(ch)){j=i+1; while(j<q.length&&/[A-Za-z0-9_-]/.test(q.charAt(j))){j++;} word=q.slice(i,j).toUpperCase(); if(word==="BETWEEN"){between=true; i=j; continue;} if(word===upper){if(upper==="AND"&&between){between=false; i=j; continue;} out.push(q.slice(start,i).replace(/^\s+|\s+$/g,"")); start=j;} i=j; continue;} i++;} out.push(q.slice(start).replace(/^\s+|\s+$/g,"")); return out;}
      function outerParens(q){var i,depth=0,quote="",ch; q=q.replace(/^\s+|\s+$/g,""); if(q.charAt(0)!=="("||q.charAt(q.length-1)!==")"){return q;} for(i=0;i<q.length;i++){ch=q.charAt(i); if(quote){if(ch===quote){quote="";} continue;} if(ch==="'"||ch==='"'){quote=ch; continue;} if(ch==="("){depth++;} else if(ch===")"){depth--; if(depth===0&&i<q.length-1){return q;}}} return depth===0?q.slice(1,-1).replace(/^\s+|\s+$/g,""):q;}
      function expr(q,a){var p,i; q=outerParens(q); if(/^NOT\s+/i.test(q)){return !expr(q.replace(/^NOT\s+/i,""),a);} p=split(q,"OR"); if(p.length>1){for(i=0;i<p.length;i++){if(expr(p[i],a)){return true;}} return false;} p=split(q,"AND"); if(p.length>1){for(i=0;i<p.length;i++){if(!expr(p[i],a)){return false;}} return true;} return pred(q,a);}
      function strip(v){v=String(v).replace(/^\s+|\s+$/g,""); if((v.charAt(0)==="'"&&v.charAt(v.length-1)==="'")||(v.charAt(0)==='"'&&v.charAt(v.length-1)==='"')){return v.slice(1,-1);} return v;}
      function pred(q,a){var m,f,x,arr,i; m=/^([A-Za-z_][A-Za-z0-9_-]*)\s+CONTAINS\s+(.+)$/i.exec(q); if(m){return contains(val(a,m[1]),strip(m[2]));} m=/^([A-Za-z_][A-Za-z0-9_-]*)\s+BETWEEN\s+(.+)\s+AND\s+(.+)$/i.exec(q); if(m){x=+val(a,m[1]); return x>=+strip(m[2])&&x<=+strip(m[3]);} m=/^([A-Za-z_][A-Za-z0-9_-]*)\s+IN\s*\((.*)\)$/i.exec(q); if(m){f=val(a,m[1]); arr=m[2].split(/\s*,\s*/); for(i=0;i<arr.length;i++){if(eq(f,strip(arr[i]))){return true;}} return false;} m=/^([A-Za-z_][A-Za-z0-9_-]*)\s*(=|!=|<=|>=|<|>)\s*(.+)$/.exec(q); if(m){return cmp(val(a,m[1]),m[2],strip(m[3]));} throw new Error("Could not parse query");}
      function val(a,f){f=String(f).toLowerCase(); if(f==="residue"){return +a.residue;} if(f==="variant_count"){return +(a.variant_count||0);} if(f==="phenotype"){return a.phenotypes||[];} if(f==="primary_phenotype"){return a.primary_phenotype||null;} if(f==="primary_significance"){return a.primary_significance||null;} if(f==="pathogenicity_class"){return path(a);} if(f==="has_pathogenic"){return !!a.has_pathogenic;} if(f==="function_class"){return a.function_classes||[];} if(f==="missense_score"){return a.missense?+a.missense.score:null;} if(f==="missense_class"){return a.missense?a.missense["class"]:null;} if(f==="alphafold_score"){return a.alphafold_confidence?+a.alphafold_confidence.score:null;} if(f==="alphafold_class"){return a.alphafold_confidence?a.alphafold_confidence["class"]:null;} throw new Error("Unknown field '"+f+"'");}
      function contains(f,n){var i; n=String(n).toLowerCase(); if(Object.prototype.toString.call(f)==="[object Array]"){for(i=0;i<f.length;i++){if(String(f[i]).toLowerCase().indexOf(n)!==-1){return true;}} return false;} return String(f||"").toLowerCase().indexOf(n)!==-1;}
      function eq(f,x){var i; if(Object.prototype.toString.call(f)==="[object Array]"){for(i=0;i<f.length;i++){if(eq(f[i],x)){return true;}} return false;} if(!isNaN(Number(f))&&!isNaN(Number(x))){return Number(f)===Number(x);} return String(f).toLowerCase()===String(x).toLowerCase();} function cmp(f,op,x){if(op==="="){return eq(f,x);} if(op==="!="){return !eq(f,x);} f=Number(f); x=Number(x); if(!isFinite(f)||!isFinite(x)){return false;} if(op==="<"){return f<x;} if(op==="<="){return f<=x;} if(op===">"){return f>x;} if(op===">="){return f>=x;} return false;}
      function path(a){var s=+(a.max_pathogenicity||0); if(!(+(a.variant_count||0))){return null;} if(a.pathogenicity_class){return a.pathogenicity_class;} if(s>=4){return "pathogenic";} if(s===3){return "likely-pathogenic";} if(s===2){return "uncertain";} return "benign";} function pathColor(c){if(c==="pathogenic"){return "#b2182b";} if(c==="likely-pathogenic"){return "#ef8a62";} if(c==="uncertain"){return "#fddbc7";} if(c==="benign"){return "#d1e5f0";} return "#d9d9d9";} function miss(c){c=String(c||"").toLowerCase(); if(c==="benign"){return "#2ca25f";} if(c==="pathogenic"){return "#d73027";} if(c==="ambiguous"){return "#fee08b";} return "#9ca3af";} function score(s){var v=Math.max(0,Math.min(1,+s||0)); return "rgb("+Math.round(43+172*v)+", "+Math.round(131-106*v)+", "+Math.round(186-158*v)+")";} function surprise(a){var s=a.missense?+a.missense.score:NaN,p=+(a.max_pathogenicity||0); if(!isFinite(s)){return null;} if(p>=3&&s<=.35){return "#7b3294";} if(p<=1&&+(a.variant_count||0)>0&&s>=.75){return "#008837";} if(p===2&&s>=.75){return "#fdae61";} return null;}
      function docs(){id("queryExamples").innerHTML='<p class="legend-note">Queries are residue-level WHERE-style expressions. Strings may be quoted, but simple hyphenated values can be typed directly.</p><p class="legend-note"><strong>Set logic:</strong> OR is union, AND is intersection, NOT is negation; parentheses can group boolean clauses.</p><p class="legend-note"><strong>Fields:</strong> all, residue, variant_count, primary_phenotype, phenotype, primary_significance, pathogenicity_class, has_pathogenic, function_class, missense_score, missense_class, alphafold_score, alphafold_class.</p><p class="legend-note"><strong>Operators:</strong> =, !=, &lt;, &lt;=, &gt;, &gt;=, IN (...), BETWEEN ... AND ..., CONTAINS.</p><code>all</code><code>residue IN (42, 43, 44)</code><code>residue BETWEEN 100 AND 120</code><code>function_class = gain-of-function</code><code>NOT function_class IN (gain-of-function, loss-of-function)</code><code>phenotype CONTAINS epileptic AND pathogenicity_class = pathogenic</code>';}
      function summary(){var u={},a=0,b=0,i,j,r; for(i=0;i<rows.length;i++){r=rows[i]; if(r.err){b++;} if(!r.en||r.ph||r.err){continue;} a++; for(j=0;j<r.m.length;j++){u[r.m[j]]=true;}} id("diseaseStats").innerHTML='<div class="stats-grid"><div class="stat"><strong>'+a+'</strong><span>active query rows</span></div><div class="stat"><strong>'+Object.keys(u).length+'</strong><span>unique matched residues</span></div><div class="stat"><strong>'+rs().length+'</strong><span>total residues</span></div><div class="stat"><strong>'+b+'</strong><span>invalid queries</span></div></div>';}
      function list(){var r=rs(),h="",i,n=Math.min(80,r.length); for(i=0;i<n;i++){h+="<li>Residue "+r[i].residue+": "+esc(r[i].primary_significance||"annotated")+", "+esc(r[i].primary_phenotype||"no phenotype")+"</li>";} id("variantList").innerHTML=h;}
      window.addEventListener("load",function(){try{init();}catch(e){stat("Viewer startup error: "+(e.message||String(e))); if(window.console&&console.error){console.error(e);}}});
    }());
  </script></body>
</html>
"""
    report_summary = {
        "gene": gene,
        "structure_label": structure_label,
        "annotation_residue_count": annotations.get("residue_count"),
        "annotation_record_count": annotations.get("record_count"),
    }
    replacements = {
        "__GENE__": escape(gene),
        "__STRUCTURE_LABEL__": escape(structure_label),
        "__REPORT_JSON__": _safe_json(report_summary),
        "__PREPROCESSED_JSON__": _safe_json(_annotations_for_viewer_script(annotations)),
        "__PDB_JSON__": _safe_json(pdb_text),
    }
    for placeholder, value in replacements.items():
        template = template.replace(placeholder, value)
    return template


def _safe_json(value: Any) -> str:
    """Serialize data for embedding in a script block."""
    return json.dumps(value).replace("</", "<\\/")


def _annotations_for_viewer_script(annotations: dict[str, Any]) -> dict[str, Any]:
    """Return a compact residue-level annotation payload for the HTML viewer."""
    compact = {
        key: value
        for key, value in annotations.items()
        if key not in {"residues", "residue_list"}
    }
    compact["residue_list"] = [
        {key: value for key, value in residue.items() if key != "variants"}
        for residue in annotations.get("residue_list", [])
        if isinstance(residue, dict)
    ]
    return compact


def _default_annotation_path(gene: str, viewer_path: Path | None = None) -> Path:
    """Return the default saved annotation path for a viewer/report."""
    if viewer_path:
        return viewer_path.with_name(f"{viewer_path.stem}_annotations.json")

    safe_gene = "".join(char if char.isalnum() or char in "-_" else "_" for char in gene.strip())
    safe_gene = safe_gene or "protein"
    return Path("outputs") / f"{safe_gene.lower()}_annotations.json"


def _extract_clinvar_records(report: dict[str, Any]) -> list[dict[str, Any]]:
    """Extract ClinVar records from either saved ClinVar JSON or CLI report JSON."""
    if isinstance(report.get("records"), list):
        return [record for record in report["records"] if isinstance(record, dict)]

    clinvar = report.get("clinvar", {})
    if not isinstance(clinvar, dict):
        return []

    data = clinvar.get("data", [])
    if isinstance(data, list):
        return [record for record in data if isinstance(record, dict)]
    if isinstance(data, dict) and isinstance(data.get("records"), list):
        return [record for record in data["records"] if isinstance(record, dict)]
    return []


def _infer_gene(records: list[dict[str, Any]]) -> str | None:
    for record in records:
        gene = record.get("gene")
        if gene:
            return str(gene)
    return None


def _group_records_by_residue(records: list[dict[str, Any]]) -> dict[int, list[dict[str, Any]]]:
    records_by_residue: dict[int, list[dict[str, Any]]] = defaultdict(list)
    for record in records:
        residue = _coerce_residue(record.get("residue"))
        if residue is not None:
            records_by_residue[residue].append(record)
    return dict(records_by_residue)


def _count_phenotypes(records: list[dict[str, Any]]) -> Counter[str]:
    counts: Counter[str] = Counter()
    for record in records:
        for phenotype in record.get("phenotype_names", []):
            if phenotype:
                counts[str(phenotype)] += 1
    return counts


def _assign_phenotype_colors(phenotype_counts: Counter[str]) -> dict[str, str]:
    phenotype_colors: dict[str, str] = {}
    for index, (phenotype, _count) in enumerate(phenotype_counts.most_common()):
        phenotype_colors[phenotype] = PHENOTYPE_PALETTE[index % len(PHENOTYPE_PALETTE)]
    return phenotype_colors


def _build_residue_annotation(
    residue: int,
    records: list[dict[str, Any]],
    phenotype_colors: dict[str, str],
    missense_prediction: dict[str, Any] | None,
    confidence_prediction: dict[str, Any] | None,
    experimental_function: dict[str, Any],
) -> dict[str, Any]:
    residue_key = str(residue)
    variant_summaries = [
        _variant_summary(record, experimental_function)
        for record in records
    ]
    phenotypes = sorted(
        {
            str(phenotype)
            for record in records
            for phenotype in record.get("phenotype_names", [])
            if phenotype
        }
    )
    significance_counts = Counter(str(record.get("significance") or "Unclassified") for record in records)
    primary_phenotype = _choose_primary_phenotype(records)
    primary_significance = _choose_primary_significance(significance_counts)
    max_pathogenicity = max(
        (_significance_severity(significance) for significance in significance_counts),
        default=0,
    )
    pathogenicity_class = _pathogenicity_class(max_pathogenicity, len(records))
    function_class_counts = Counter()
    residue_function_entry = experimental_function["by_residue"].get(residue_key, {})
    if residue_function_entry.get("function_class_counts"):
        function_class_counts.update(residue_function_entry["function_class_counts"])
    else:
        for variant in variant_summaries:
            for function_class in variant.get("function_classes", []):
                function_class_counts[function_class] += 1
    for function_class in residue_function_entry.get("function_classes", []):
        function_class_counts.setdefault(function_class, 0)

    return {
        "residue": residue,
        "variant_count": len(records),
        "phenotypes": phenotypes,
        "primary_phenotype": primary_phenotype,
        "phenotype_color": phenotype_colors.get(primary_phenotype, "#9ca3af"),
        "significance_counts": dict(significance_counts),
        "primary_significance": primary_significance,
        "has_pathogenic": max_pathogenicity >= 3,
        "max_pathogenicity": max_pathogenicity,
        "pathogenicity_class": pathogenicity_class,
        "pathogenic_color": _pathogenicity_color(max_pathogenicity, len(records)),
        "function_classes": sorted(function_class_counts),
        "function_class_counts": dict(function_class_counts),
        "missense": missense_prediction,
        "alphafold_confidence": confidence_prediction,
        "variants": variant_summaries,
    }


def _variant_summary(record: dict[str, Any], experimental_function: dict[str, Any]) -> dict[str, Any]:
    function_entry = _experimental_function_for_record(record, experimental_function)
    return {
        "variation_id": record.get("variation_id"),
        "accession": record.get("accession"),
        "protein_change": record.get("protein_change"),
        "amino_acid": record.get("amino_acid"),
        "alternate_amino_acid": record.get("alternate_amino_acid"),
        "significance": record.get("significance"),
        "phenotype_names": record.get("phenotype_names", []),
        "molecular_consequences": record.get("molecular_consequences", []),
        "function_classes": function_entry.get("function_classes", []),
        "function_calls": function_entry.get("function_calls", []),
    }


def _experimental_function_for_record(record: dict[str, Any], experimental_function: dict[str, Any]) -> dict[str, Any]:
    by_variant = experimental_function["by_variant"]
    for key in _record_variant_lookup_keys(record):
        if key in by_variant:
            return by_variant[key]
    return {}


def _record_variant_lookup_keys(record: dict[str, Any]) -> list[str]:
    keys: list[str] = []
    protein_change = _normalize_variant_name(record.get("protein_change"))
    if protein_change:
        keys.append(protein_change)
    residue = _coerce_residue(record.get("residue"))
    amino_acid = str(record.get("amino_acid") or "").strip()
    alternate_amino_acid = str(record.get("alternate_amino_acid") or "").strip()
    if residue is not None and amino_acid and alternate_amino_acid:
        keys.append(_normalize_variant_name(f"{amino_acid}{residue}{alternate_amino_acid}"))
    return keys


def _exp_variance_patient_data_path(gene: str) -> Path:
    normalized_gene = normalize_gene(gene)
    return Path("data") / "exp_variance" / normalized_gene / f"{normalized_gene}_patient_data.csv"


def _scn1a_bosselman_path() -> Path:
    return Path("data") / "functional_variants" / "SCN1A" / "bosselman.csv"


def _scn1a_gof_ndeema_path() -> Path:
    return Path("data") / "patient_variants" / "SCN1A" / "gof_ndeema.csv"


def _scn1a_brunklaus_path() -> Path:
    return Path("data") / "patient_variants" / "SCN1A" / "dravet_gefs_brunklaus.csv"


def _experimental_function_source_paths(gene: str) -> list[Path]:
    normalized_gene = normalize_gene(gene)
    paths = []
    exp_variance_path = _exp_variance_patient_data_path(normalized_gene)
    if exp_variance_path.exists():
        paths.append(exp_variance_path)
    if normalized_gene == "SCN1A":
        for path in [_scn1a_bosselman_path(), _scn1a_gof_ndeema_path(), _scn1a_brunklaus_path()]:
            if path.exists():
                paths.append(path)
    return paths


def _load_experimental_function_data(gene: str) -> dict[str, Any]:
    normalized_gene = normalize_gene(gene)
    source_paths = _experimental_function_source_paths(normalized_gene)
    if not source_paths:
        return {
            "source": "none",
            "path": str(_exp_variance_patient_data_path(normalized_gene)),
            "paths": [],
            "patient_row_count": 0,
            "variant_count": 0,
            "residue_count": 0,
            "class_counts": {},
            "source_counts": {},
            "by_variant": {},
            "by_residue": {},
        }

    by_variant: dict[str, dict[str, Any]] = {}
    by_residue: dict[str, dict[str, Any]] = {}
    class_counts: Counter[str] = Counter()
    source_counts: Counter[str] = Counter()
    patient_row_count = 0

    def add_call(
        *,
        variant: Any,
        function_call: Any,
        source_label: str,
        residue: Any = None,
    ) -> bool:
        variant_key = _normalize_variant_name(variant)
        function_classes = _function_classes_from_call(function_call)
        if not variant_key or not function_classes:
            return False
        resolved_residue = _coerce_residue(residue)
        if resolved_residue is None:
            resolved_residue = _residue_from_variant_name(variant_key)
        call_label = str(function_call or "").strip()
        if source_label:
            call_label = f"{call_label} [{source_label}]" if call_label else source_label
        entry = by_variant.setdefault(
            variant_key,
            {"function_classes": set(), "function_calls": set(), "patient_count": 0},
        )
        entry["function_classes"].update(function_classes)
        if call_label:
            entry["function_calls"].add(call_label)
        entry["patient_count"] += 1
        source_counts[source_label or "unknown"] += 1
        for function_class in function_classes:
            class_counts[function_class] += 1
        if resolved_residue is not None:
            residue_entry = by_residue.setdefault(
                str(resolved_residue),
                {
                    "function_classes": set(),
                    "function_calls": set(),
                    "variant_keys": set(),
                    "function_class_counts": Counter(),
                },
            )
            residue_entry["function_classes"].update(function_classes)
            residue_entry["variant_keys"].add(variant_key)
            if call_label:
                residue_entry["function_calls"].add(call_label)
            for function_class in function_classes:
                residue_entry["function_class_counts"][function_class] += 1
        return True

    exp_variance_path = _exp_variance_patient_data_path(normalized_gene)
    if exp_variance_path.exists():
        with exp_variance_path.open("r", encoding="utf-8-sig", newline="") as handle:
            reader = csv.DictReader(handle)
            for row in reader:
                patient_row_count += 1
                add_call(
                    variant=row.get("variant"),
                    function_call=row.get("function_call"),
                    source_label="exp_variance_patient_data",
                )

    if normalized_gene == "SCN1A":
        bosselman_path = _scn1a_bosselman_path()
        if bosselman_path.exists():
            with bosselman_path.open("r", encoding="utf-8-sig", newline="") as handle:
                reader = csv.DictReader(handle)
                for row in reader:
                    patient_row_count += 1
                    add_call(
                        variant=row.get("protein_change"),
                        function_call=row.get("function_call"),
                        source_label="bosselman_functional",
                    )

        gof_path = _scn1a_gof_ndeema_path()
        if gof_path.exists():
            with gof_path.open("r", encoding="utf-8-sig", newline="") as handle:
                reader = csv.DictReader(handle)
                for row in reader:
                    patient_row_count += 1
                    add_call(
                        variant=row.get("variant"),
                        function_call=row.get("function_call"),
                        source_label="gof_ndeema",
                    )

        brunklaus_path = _scn1a_brunklaus_path()
        if brunklaus_path.exists():
            with brunklaus_path.open("r", encoding="utf-8-sig", newline="") as handle:
                reader = csv.DictReader(handle)
                for row in reader:
                    patient_row_count += 1
                    dx = str(row.get("dx") or "").strip()
                    if dx not in {"Dravet", "GEFS+"}:
                        continue
                    add_call(
                        variant=row.get("protein_change"),
                        residue=row.get("residue_canonical") or row.get("residue"),
                        function_call=f"{dx} inferred LOF",
                        source_label="brunklaus_dravet_gefs",
                    )

    return {
        "source": "combined_experimental_function_data",
        "path": str(source_paths[0]) if len(source_paths) == 1 else [str(path) for path in source_paths],
        "paths": [str(path) for path in source_paths],
        "patient_row_count": patient_row_count,
        "variant_count": len(by_variant),
        "residue_count": len(by_residue),
        "class_counts": dict(class_counts),
        "source_counts": dict(source_counts),
        "by_variant": {
            key: {
                "function_classes": sorted(value["function_classes"]),
                "function_calls": sorted(value["function_calls"]),
                "patient_count": value["patient_count"],
            }
            for key, value in by_variant.items()
        },
        "by_residue": {
            key: {
                "function_classes": sorted(value["function_classes"]),
                "function_calls": sorted(value["function_calls"]),
                "function_class_counts": dict(value["function_class_counts"]),
                "variant_keys": sorted(value["variant_keys"]),
            }
            for key, value in by_residue.items()
        },
    }


def _function_classes_from_call(function_call: Any) -> list[str]:
    normalized = str(function_call or "").strip().lower()
    classes: list[str] = []
    if "gof" in normalized or "gain" in normalized:
        classes.append("gain-of-function")
    if "lof" in normalized or "loss" in normalized:
        classes.append("loss-of-function")
    return classes


def _normalize_variant_name(value: Any) -> str | None:
    if value is None:
        return None
    text = str(value).strip()
    if not text:
        return None
    text = re.sub(r"^p\.", "", text, flags=re.IGNORECASE)
    text = re.sub(r"[^A-Za-z0-9=*]+", "", text)
    return text.upper() or None


def _residue_from_variant_name(variant_name: str | None) -> int | None:
    if not variant_name:
        return None
    match = re.search(r"(\d+)", variant_name)
    if not match:
        return None
    return int(match.group(1))


def _choose_primary_phenotype(records: list[dict[str, Any]]) -> str | None:
    counts: Counter[str] = Counter()
    fallback_counts: Counter[str] = Counter()
    for record in records:
        for phenotype in record.get("phenotype_names", []):
            if not phenotype:
                continue
            phenotype_name = str(phenotype)
            fallback_counts[phenotype_name] += 1
            if phenotype_name not in GENERIC_PHENOTYPES:
                counts[phenotype_name] += 1

    source = counts or fallback_counts
    if not source:
        return None
    return sorted(source.items(), key=lambda item: (-item[1], item[0]))[0][0]


def _choose_primary_significance(significance_counts: Counter[str]) -> str | None:
    if not significance_counts:
        return None
    return sorted(
        significance_counts.items(),
        key=lambda item: (-_significance_severity(item[0]), -item[1], item[0]),
    )[0][0]


def _significance_severity(significance: str | None) -> int:
    normalized = str(significance or "").lower()
    if "pathogenic/likely pathogenic" in normalized:
        return 4
    if "pathogenic" in normalized and "likely" not in normalized:
        return 4
    if "likely pathogenic" in normalized:
        return 3
    if "conflicting" in normalized:
        return 2
    if "uncertain" in normalized:
        return 2
    if "likely benign" in normalized:
        return 1
    if "benign" in normalized:
        return 0
    return 0


def _pathogenicity_class(severity: int, variant_count: int) -> str | None:
    if variant_count <= 0:
        return None
    if severity >= 4:
        return "pathogenic"
    if severity == 3:
        return "likely-pathogenic"
    if severity == 2:
        return "uncertain"
    return "benign"


def _pathogenicity_color(severity: int, variant_count: int) -> str:
    if variant_count <= 0:
        return "#d9d9d9"
    if severity >= 4:
        return "#b2182b"
    if severity == 3:
        return "#ef8a62"
    if severity == 2:
        return "#fddbc7"
    if severity == 1:
        return "#d1e5f0"
    return "#d9d9d9"


def _build_missense_gradient(
    records_by_residue: dict[int, list[dict[str, Any]]],
    report: dict[str, Any],
) -> tuple[dict[str, dict[str, Any]], str]:
    real_predictions = _extract_missense_predictions(report)
    if real_predictions:
        return real_predictions, "alphamissense"

    residues = sorted(records_by_residue)
    if not residues:
        return {}, "stub_linear_by_residue_position"

    minimum = residues[0]
    maximum = residues[-1]
    span = max(maximum - minimum, 1)
    return {
        str(residue): {
            "residue": residue,
            "score": round((residue - minimum) / span, 4),
            "class": "stub",
            "color": _score_to_color((residue - minimum) / span),
        }
        for residue in residues
    }, "stub_linear_by_residue_position"


def _extract_missense_predictions(report: dict[str, Any]) -> dict[str, dict[str, Any]]:
    alpha_missense = report.get("alpha_missense", {})
    if not isinstance(alpha_missense, dict):
        return {}

    data = alpha_missense.get("data", [])
    if not isinstance(data, list):
        return {}

    predictions: dict[str, dict[str, Any]] = {}
    for record in data:
        if not isinstance(record, dict):
            continue
        residue = _coerce_residue(record.get("residue") or record.get("position"))
        score = _coerce_score(
            _first_present(
                record,
                "score",
                "mean_pathogenicity",
                "am_pathogenicity",
                "pathogenicity_score",
            )
        )
        if residue is None or score is None:
            continue
        predictions[str(residue)] = {
            "residue": residue,
            "score": score,
            "class": record.get("am_class") or record.get("class") or record.get("mean_class"),
            "color": _score_to_color(score),
        }

    return predictions


def _first_present(record: dict[str, Any], *keys: str) -> Any:
    for key in keys:
        if key in record and record[key] is not None:
            return record[key]
    return None


def _score_to_color(score: float) -> str:
    clamped = max(0.0, min(1.0, score))
    low = (43, 131, 186)
    high = (215, 25, 28)
    red = round(low[0] + (high[0] - low[0]) * clamped)
    green = round(low[1] + (high[1] - low[1]) * clamped)
    blue = round(low[2] + (high[2] - low[2]) * clamped)
    return f"#{red:02x}{green:02x}{blue:02x}"


def _confidence_class(score: float) -> str:
    if score >= 90:
        return "very high"
    if score >= 70:
        return "confident"
    if score >= 50:
        return "low"
    return "very low"


def _confidence_color(score: float) -> str:
    normalized = max(0.0, min(100.0, score)) / 100.0
    low = (214, 39, 40)
    high = (31, 119, 180)
    red = round(low[0] + (high[0] - low[0]) * normalized)
    green = round(low[1] + (high[1] - low[1]) * normalized)
    blue = round(low[2] + (high[2] - low[2]) * normalized)
    return f"#{red:02x}{green:02x}{blue:02x}"


def _coerce_residue(value: Any) -> int | None:
    if value is None:
        return None
    try:
        residue = int(value)
    except (TypeError, ValueError):
        return None
    return residue if residue > 0 else None


def _coerce_score(value: Any) -> float | None:
    if value is None:
        return None
    try:
        score = float(value)
    except (TypeError, ValueError):
        return None
    return max(0.0, min(1.0, score))


def _coerce_score_100(value: Any) -> float | None:
    if value is None:
        return None
    try:
        score = float(value)
    except (TypeError, ValueError):
        return None
    return max(0.0, min(100.0, score))


def _source_has_data(report: dict[str, Any], source_name: str) -> bool:
    """Return whether a report source currently contains non-empty data."""
    source = report.get(source_name, {})
    if not isinstance(source, dict):
        return False

    data = source.get("data", [])
    return bool(data)
