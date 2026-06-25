"""Local 3Dmol.js visualization helpers for StructPhenotypes."""

from __future__ import annotations

from collections import Counter, defaultdict
import json
from html import escape
from pathlib import Path
from typing import Any


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
            "Base color mode: gray protein, default, AlphaMissense gradient/class, AlphaFold confidence, ClinVar pathogenicity",
            "Rendering style: cartoon, trace, tube, cartoon plus sticks, sticks, lines, sphere/spacefill, cross, VDW/MS/SAS/SES surfaces, plus a smoother SAS-based surface option",
            "Phenotype checkbox overlay with per-phenotype color pickers",
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
        )
        residues[residue_key] = annotation
        residue_list.append(annotation)

    return {
        "gene": gene,
        "source": "StructPhenotypes visualization preprocessing",
        "annotation_version": 1,
        "record_count": len(records),
        "residue_count": len(residues),
        "phenotype_counts": dict(phenotype_counts),
        "phenotype_colors": phenotype_colors,
        "pathogenic_labels": sorted(PATHOGENIC_LABELS),
        "residues": residues,
        "residue_list": residue_list,
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

    if _source_has_data(report, "alpha_missense"):
        missense = annotations.get("missense_predictions", {})
        if missense.get("source") != "alphamissense":
            return False

    if _source_has_data(report, "alpha_fold"):
        confidence = annotations.get("alphafold_confidence", {})
        if confidence.get("source") != "alphafold_pdb_b_factor":
            return False

    return True


def _pdb_text_for_report(report: dict[str, Any]) -> tuple[str, str]:
    """Return PDB text from AlphaFold metadata or raise if unavailable."""
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

    raise FileNotFoundError("AlphaFold structure is required for visualization but no PDB file was found.")


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
    """Render a standalone HTML viewer around example protein data."""
    gene = str(report.get("gene", "unknown"))
    structure_label = "AlphaFold PDB structure" if pdb_text != EXAMPLE_PDB else "embedded example protein structure"
    report_json = _safe_json(report)
    preprocessed_json = _safe_json(annotations)
    pdb_json = _safe_json(pdb_text)
    annotations_json = _safe_json(EXAMPLE_ANNOTATIONS)

    return f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>StructPhenotypes Viewer - {escape(gene)}</title>
  <script src="https://3Dmol.org/build/3Dmol-min.js"></script>
  <style>
    :root {{
      color-scheme: light;
      font-family: Arial, sans-serif;
      --border: #d8dee8;
      --ink: #1f2937;
      --muted: #607085;
      --panel: #f7f9fc;
      --blue: #2459d6;
    }}
    body {{
      margin: 0;
      color: var(--ink);
      background: #ffffff;
    }}
    header {{
      padding: 18px 22px 12px;
      border-bottom: 1px solid var(--border);
    }}
    h1 {{
      margin: 0 0 6px;
      font-size: 22px;
      font-weight: 700;
    }}
    .subtitle {{
      margin: 0;
      color: var(--muted);
      font-size: 14px;
    }}
    .toolbar {{
      display: flex;
      flex-wrap: wrap;
      gap: 12px;
      align-items: end;
      padding: 14px 22px;
      background: var(--panel);
      border-bottom: 1px solid var(--border);
    }}
    label {{
      display: grid;
      gap: 5px;
      font-size: 12px;
      font-weight: 700;
      text-transform: uppercase;
      letter-spacing: 0.04em;
      color: var(--muted);
    }}
    select,
    button {{
      min-width: 180px;
      height: 36px;
      border: 1px solid var(--border);
      border-radius: 6px;
      background: #ffffff;
      color: var(--ink);
      font-size: 14px;
      padding: 0 10px;
    }}
    button {{
      min-width: 110px;
      cursor: pointer;
      color: #ffffff;
      background: var(--blue);
      border-color: var(--blue);
      font-weight: 700;
    }}
    #viewer {{
      width: 100%;
      height: 68vh;
      min-height: 430px;
      position: relative;
    }}
    #status {{
      padding: 10px 22px;
      color: var(--muted);
      border-top: 1px solid var(--border);
      border-bottom: 1px solid var(--border);
      font-size: 14px;
    }}
    .details {{
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(230px, 1fr));
      gap: 14px;
      padding: 16px 22px 24px;
    }}
    .panel {{
      border: 1px solid var(--border);
      border-radius: 8px;
      padding: 12px;
    }}
    .panel h2 {{
      margin: 0 0 8px;
      font-size: 15px;
    }}
    .panel p,
    .panel li {{
      color: var(--muted);
      font-size: 14px;
      line-height: 1.4;
    }}
    .panel ul {{
      margin: 0;
      padding-left: 18px;
    }}
    @media (max-width: 700px) {{
      select,
      button {{
        min-width: 100%;
      }}
      .toolbar {{
        align-items: stretch;
      }}
      label {{
        width: 100%;
      }}
    }}
  </style>
</head>
<body>
  <header>
    <h1>StructPhenotypes Viewer: {escape(gene)}</h1>
    <p class="subtitle">Local HTML prototype using 3Dmol.js and {escape(structure_label)}.</p>
  </header>

  <div class="toolbar">
    <label>
      Color mode
      <select id="colorMode">
        <option value="default">Default protein view</option>
        <option value="phenotype">ClinVar phenotype</option>
        <option value="significance">ClinVar significance</option>
        <option value="missense">AlphaMissense gradient</option>
      </select>
    </label>
    <label>
      Phenotype filter
      <select id="phenotypeFilter"></select>
    </label>
    <label>
      Focus residue
      <select id="focusResidue"></select>
    </label>
    <button id="resetView" type="button">Reset view</button>
  </div>

  <div id="viewer"></div>
  <div id="status">Loading viewer...</div>

  <section class="details">
    <div class="panel">
      <h2>Current Prototype</h2>
      <p>This page uses preprocessed residue annotations and the best available local structure file.</p>
    </div>
    <div class="panel">
      <h2>Variant Layers</h2>
      <ul id="variantList"></ul>
    </div>
    <div class="panel">
      <h2>Next Data Hooks</h2>
      <p>ClinVar residues are preprocessed into residue-level annotations. AlphaMissense scores are currently represented by a stub residue gradient.</p>
    </div>
  </section>

  <script>
    const REPORT = {report_json};
    const PREPROCESSED_ANNOTATIONS = {preprocessed_json};
    const PDB_TEXT = {pdb_json};
    const ANNOTATIONS = {annotations_json};

    let viewer = null;

    function init() {{
      populateControls();
      populateVariantList();

      if (!window.$3Dmol) {{
        setStatus("3Dmol.js did not load. Check your internet connection for the viewer library.");
        return;
      }}

      viewer = $3Dmol.createViewer("viewer", {{ backgroundColor: "white" }});
      viewer.addModel(PDB_TEXT, "pdb");
      applyStyle();
      viewer.zoomTo();
      viewer.render();
      setStatus("Example protein loaded. Use the dropdowns to switch annotation layers.");
    }}

    function populateControls() {{
      const phenotypeSelect = document.getElementById("phenotypeFilter");
      const phenotypes = Object.keys(PREPROCESSED_ANNOTATIONS.phenotype_counts || {{}}).sort();
      phenotypeSelect.innerHTML = '<option value="all">All phenotypes</option>' +
        phenotypes.map((phenotype) => `<option value="${{escapeAttr(phenotype)}}">${{phenotype}}</option>`).join("");

      const residueSelect = document.getElementById("focusResidue");
      const residues = getResidueAnnotations();
      residueSelect.innerHTML = '<option value="">No residue focus</option>' +
        residues
          .map((annotation) => `<option value="${{annotation.residue}}">Residue ${{annotation.residue}}: ${{annotation.primary_significance || "annotated"}}</option>`)
          .join("");

      document.getElementById("colorMode").addEventListener("change", applyStyle);
      phenotypeSelect.addEventListener("change", applyStyle);
      residueSelect.addEventListener("change", focusResidue);
      document.getElementById("resetView").addEventListener("click", resetView);
    }}

    function getResidueAnnotations() {{
      const residues = PREPROCESSED_ANNOTATIONS.residue_list || [];
      if (residues.length > 0) {{
        return residues;
      }}

      return ANNOTATIONS.clinvar_variants.map((variant) => ({{
        residue: variant.residue,
        variant_count: 1,
        phenotypes: [variant.phenotype],
        primary_phenotype: variant.phenotype,
        phenotype_color: variant.color,
        primary_significance: variant.significance,
        has_pathogenic: true,
        max_pathogenicity: 4,
        pathogenic_color: significanceColor(variant.significance),
        missense: null,
      }}));
    }}

    function populateVariantList() {{
      const list = document.getElementById("variantList");
      const residues = getResidueAnnotations();
      const highlighted = residues.filter((annotation) => annotation.has_pathogenic).slice(0, 100);
      const displayed = highlighted.length > 0 ? highlighted : residues.slice(0, 100);
      list.innerHTML = displayed
        .map((annotation) => `<li>Residue ${{annotation.residue}}: ${{annotation.primary_significance || "annotated"}}, ${{annotation.primary_phenotype || "no phenotype"}} (${{annotation.variant_count}} variants)</li>`)
        .join("");
    }}

    function applyStyle() {{
      if (!viewer) {{
        return;
      }}

      const mode = document.getElementById("colorMode").value;
      const phenotype = document.getElementById("phenotypeFilter").value;

      viewer.setStyle({{}}, {{ cartoon: {{ color: "spectrum" }}, stick: {{ radius: 0.08 }} }});

      if (mode === "phenotype") {{
        colorClinVarPhenotypes(phenotype);
        setStatus("Coloring example ClinVar residues by phenotype.");
      }} else if (mode === "significance") {{
        colorClinVarSignificance();
        setStatus("Coloring example ClinVar residues by significance.");
      }} else if (mode === "missense") {{
        colorMissenseGradient();
        setStatus("Coloring example residues by AlphaMissense-like prediction score.");
      }} else {{
        setStatus("Default spectrum protein view.");
      }}

      focusResidue(false);
      viewer.render();
    }}

    function colorClinVarPhenotypes(phenotype) {{
      getResidueAnnotations().forEach((annotation) => {{
        if (phenotype !== "all" && !annotation.phenotypes.includes(phenotype)) {{
          return;
        }}
        setResidueStyle(annotation.residue, annotation.phenotype_color || "#9ca3af", 0.2);
      }});
    }}

    function colorClinVarSignificance() {{
      getResidueAnnotations().forEach((annotation) => {{
        if (annotation.max_pathogenicity <= 0) {{
          return;
        }}
        setResidueStyle(annotation.residue, annotation.pathogenic_color || "#999999", 0.2);
      }});
    }}

    function colorMissenseGradient() {{
      getResidueAnnotations().forEach((annotation) => {{
        if (!annotation.missense) {{
          return;
        }}
        setResidueStyle(annotation.residue, annotation.missense.color || scoreToColor(annotation.missense.score), 0.12);
      }});
    }}

    function setResidueStyle(residue, color, radius) {{
      viewer.setStyle(
        {{ chain: "A", resi: residue }},
        {{ cartoon: {{ color }}, stick: {{ color, radius }}, sphere: {{ color, radius: radius * 1.8 }} }}
      );
    }}

    function focusResidue(zoom = true) {{
      if (!viewer) {{
        return;
      }}

      const residueValue = document.getElementById("focusResidue").value;
      if (!residueValue) {{
        return;
      }}

      const residue = Number(residueValue);
      viewer.addStyle({{ chain: "A", resi: residue }}, {{ sphere: {{ color: "#2459d6", radius: 0.55 }} }});
      if (zoom) {{
        viewer.zoomTo({{ chain: "A", resi: residue }});
      }}
      viewer.render();
    }}

    function resetView() {{
      if (!viewer) {{
        return;
      }}
      document.getElementById("focusResidue").value = "";
      applyStyle();
      viewer.zoomTo();
      viewer.render();
    }}

    function significanceColor(significance) {{
      const normalized = significance.toLowerCase();
      if (normalized.includes("pathogenic") && !normalized.includes("likely")) {{
        return "#b2182b";
      }}
      if (normalized.includes("likely pathogenic")) {{
        return "#ef8a62";
      }}
      if (normalized.includes("uncertain")) {{
        return "#fddbc7";
      }}
      return "#999999";
    }}

    function scoreToColor(score) {{
      const clamped = Math.max(0, Math.min(1, score));
      const red = Math.round(255 * clamped);
      const green = Math.round(180 * (1 - clamped));
      const blue = Math.round(255 * (1 - clamped));
      return `rgb(${{red}}, ${{green}}, ${{blue}})`;
    }}

    function escapeAttr(value) {{
      return String(value).replaceAll("&", "&amp;").replaceAll('"', "&quot;").replaceAll("<", "&lt;");
    }}

    function setStatus(message) {{
      document.getElementById("status").textContent = message;
    }}

    window.addEventListener("load", init);
  </script>
</body>
</html>
"""


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
      height: 66vh;
      min-height: 430px;
      position: relative;
    }
    #status {
      padding: 10px 22px;
      color: var(--muted);
      border-top: 1px solid var(--border);
      border-bottom: 1px solid var(--border);
      font-size: 14px;
    }
    .details {
      display: grid;
      grid-template-columns: minmax(260px, 360px) 1fr;
      gap: 14px;
      padding: 16px 22px 24px;
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
      max-height: 260px;
      overflow: auto;
      display: grid;
      gap: 7px;
      padding-right: 6px;
    }
    .phenotype-row {
      display: grid;
      grid-template-columns: 22px 42px 1fr auto;
      gap: 8px;
      align-items: center;
      border-bottom: 1px solid #edf1f7;
      padding: 4px 0;
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
        .details {
          grid-template-columns: 1fr;
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
        <input id="backgroundOpacity" type="range" min="0.05" max="1" step="0.05" value="0.3">
        <span id="backgroundOpacityValue" class="slider-value">0.30</span>
      </label>
      <label>
        Residue marker size
        <div class="inline-field">
          <input id="markerScale" type="range" min="0.5" max="5" step="0.1" value="1.6">
          <input id="markerScaleText" type="text" value="1.6x" inputmode="decimal" aria-label="Residue marker size value">
        </div>
        <span id="markerScaleValue" class="slider-value">1.6x</span>
      </label>
      <button id="resetView" type="button">Reset view</button>
    </div>

  <div id="viewer"></div>
  <div id="status">Loading viewer...</div>

  <section class="details">
    <div class="panel">
      <h2>Phenotype Sets</h2>
        <p>Tick phenotypes to overlay colored residue markers on top of the current scaffold coloring.</p>
      <div id="phenotypeControls" class="phenotype-panel"></div>
    </div>
    <div class="panel">
      <h2>Legend</h2>
      <div id="colorLegend"></div>
      <h2>Residue Summary</h2>
      <ul id="variantList"></ul>
    </div>
  </section>

  <script>
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
    const PHENOTYPE_LIMIT = 80;

      let viewer = null;
      let model = null;
      let phenotypeEntries = [];
      let currentFocusResidue = null;
      let currentSceneKey = null;

      function init() {
        populateControls();
        populateVariantList();

      if (!window.$3Dmol) {
        setStatus("3Dmol.js did not load. Check your internet connection for the viewer library.");
        return;
      }

        viewer = $3Dmol.createViewer("viewer", { backgroundColor: "white" });
        viewer.addModel(PDB_TEXT, "pdb");
        model = viewer.getModel();
        applyStyle(true);
        viewer.zoomTo();
        viewer.render();
        setStatus("Protein loaded. Choose a base color mode, phenotype sets, or a focused residue.");
      }

      function populateControls() {
        populatePhenotypeControls();
        populateResidueFocus();
        updateBackgroundOpacityLabel();
        updateMarkerScaleLabel();

        document.getElementById("baseColorMode").addEventListener("change", () => applyStyle(false));
        document.getElementById("renderStyle").addEventListener("change", () => applyStyle(false));
        document.getElementById("focusResidue").addEventListener("change", () => handleFocusChange(true));
        document.getElementById("backgroundOpacity").addEventListener("input", () => {
          updateBackgroundOpacityLabel();
          applyStyle(false);
        });
        document.getElementById("markerScale").addEventListener("input", () => {
          updateMarkerScaleLabel();
          applyStyle(false);
        });
        document.getElementById("markerScaleText").addEventListener("change", () => {
          syncMarkerScaleFromText();
          applyStyle(false);
        });
        document.getElementById("resetView").addEventListener("click", resetView);
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

        panel.innerHTML = phenotypeEntries.map((entry, index) => `
          <div class="phenotype-row">
            <input class="phenotype-check" type="checkbox" data-index="${index}">
            <input class="phenotype-color" type="color" data-index="${index}" value="${entry.color}">
            <span class="phenotype-name">${escapeHtml(entry.name)}</span>
          <span class="phenotype-count">${entry.count}</span>
        </div>
      `).join("");

        panel.querySelectorAll("input").forEach((input) => {
          input.addEventListener("change", handlePhenotypeOverlayChange);
        });
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

    function selectedPhenotypes() {
      const selected = new Map();
      document.querySelectorAll(".phenotype-check:checked").forEach((checkbox) => {
        const index = Number(checkbox.dataset.index);
        const entry = phenotypeEntries[index];
        const colorInput = document.querySelector(`.phenotype-color[data-index="${index}"]`);
        if (entry && colorInput) {
          selected.set(entry.name, colorInput.value);
        }
      });
      return selected;
    }

    function populateVariantList() {
      const list = document.getElementById("variantList");
      const residues = getResidueAnnotations();
      const highlighted = residues.filter((annotation) => annotation.has_pathogenic).slice(0, 100);
      const displayed = highlighted.length > 0 ? highlighted : residues.slice(0, 100);
      list.innerHTML = displayed
        .map((annotation) => `<li>Residue ${annotation.residue}: ${annotation.primary_significance || "annotated"}, ${annotation.primary_phenotype || "no phenotype"} (${annotation.variant_count || 0} variants)</li>`)
        .join("");
    }

      function applyStyle(zoomFocus) {
        if (!viewer) {
          return;
        }

        const focusResidue = selectedFocusResidue();
        const restOpacity = focusResidue ? selectedBackgroundOpacity() : 1.0;
        clearSurfaces();
        clearShapes();
        viewer.setStyle({}, {});

        renderScene(restOpacity);
        currentFocusResidue = focusResidue;
        currentSceneKey = sceneKey(restOpacity);
        redrawShapeOverlays();
        applyFocus(focusResidue, zoomFocus);
        updateLegend();
        viewer.render();
      }

      function handlePhenotypeOverlayChange() {
        if (!viewer) {
          return;
        }
        const focusResidue = selectedFocusResidue();
        clearShapes();
        redrawShapeOverlays();
        applyFocus(focusResidue, false);
        updateLegend();
        viewer.render();
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
        const focusResidue = selectedFocusResidue();
        const opacity = focusResidue ? selectedBackgroundOpacity() : 1.0;
        const phenotypes = selectedPhenotypes();
        if (isSurfaceStyle()) {
          addSurfaceModeMarkers(opacity);
        }
        if (phenotypes.size > 0) {
          addPhenotypeMarkers(phenotypes);
        }
      }

      function sceneKey(opacity) {
        return JSON.stringify({
          renderStyle: document.getElementById("renderStyle").value,
          baseColorMode: document.getElementById("baseColorMode").value,
          opacity: Number(opacity.toFixed(2)),
        });
      }

    function applyBaseColorMode(opacity) {
      const mode = document.getElementById("baseColorMode").value;

      if (mode === "gray") {
        styleWholeProtein("#d9d9d9", opacity);
        setStatus("Neutral gray protein.");
        return;
      }

      if (mode === "default") {
        if (isSurfaceStyle()) {
          styleWholeProtein("#d9d9d9", opacity);
          setStatus("Surface styles use a gray protein base. Choose a data mode to overlay colored residues.");
          return;
        }
        styleWholeProtein("spectrum", opacity);
        setStatus("Default spectrum coloring.");
        return;
      }

      styleWholeProtein("#d9d9d9", opacity);

      if (mode === "missense-gradient") {
        getResidueAnnotations().forEach((annotation) => {
          if (annotation.missense) {
            styleResidue(annotation.residue, annotation.missense.color || scoreToColor(annotation.missense.score), opacity, 0.1);
          }
        });
        setStatus("Coloring by AlphaMissense mean pathogenicity score from 0 to 1.");
      } else if (mode === "missense-class") {
        getResidueAnnotations().forEach((annotation) => {
          if (annotation.missense) {
            styleResidue(annotation.residue, missenseClassColor(annotation.missense.class), opacity, 0.1);
          }
        });
        setStatus("Coloring by AlphaMissense class.");
      } else if (mode === "alphafold-confidence") {
        getResidueAnnotations().forEach((annotation) => {
          if (annotation.alphafold_confidence) {
            styleResidue(annotation.residue, annotation.alphafold_confidence.color, opacity, 0.08);
          }
        });
        setStatus("Coloring by AlphaFold model confidence from the PDB B-factor column.");
      } else if (mode === "clinvar-pathogenicity") {
        getResidueAnnotations().forEach((annotation) => {
          if (annotation.max_pathogenicity > 0) {
            styleResidue(annotation.residue, annotation.pathogenic_color || "#999999", opacity, 0.14);
          }
        });
        setStatus("Coloring ClinVar residues by strongest pathogenicity label.");
      }
    }

      function addPhenotypeMarkers(phenotypes) {
        getResidueAnnotations().forEach((annotation) => {
          const matchingPhenotypes = annotation.phenotypes
            .filter((phenotype) => phenotypes.has(phenotype))
            .sort();
          if (matchingPhenotypes.length === 0) {
            return;
          }
          const radius = markerRadius(annotation, 0.38);
          matchingPhenotypes.forEach((phenotype, index) => {
            styleResidueMarker(
              annotation.residue,
              phenotypes.get(phenotype),
              1.0,
              radius,
              phenotypeMarkerCenter(annotation.residue, phenotype, index, matchingPhenotypes.length, radius)
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
      }

      function addSurfaceModeMarkers(opacity) {
        const mode = document.getElementById("baseColorMode").value;
        if (mode === "gray" || mode === "default") {
          return;
        }
        getResidueAnnotations().forEach((annotation) => {
          if (mode === "missense-gradient" && annotation.missense) {
            styleResidue(annotation.residue, annotation.missense.color || scoreToColor(annotation.missense.score), opacity, 0.1);
          } else if (mode === "missense-class" && annotation.missense) {
            styleResidue(annotation.residue, missenseClassColor(annotation.missense.class), opacity, 0.1);
          } else if (mode === "alphafold-confidence" && annotation.alphafold_confidence) {
            styleResidue(annotation.residue, annotation.alphafold_confidence.color, opacity, 0.08);
          } else if (mode === "clinvar-pathogenicity" && annotation.max_pathogenicity > 0) {
            styleResidue(annotation.residue, annotation.pathogenic_color || "#999999", opacity, 0.14);
          }
        });
      }

      function updateLegend() {
        const legend = document.getElementById("colorLegend");
        const phenotypes = selectedPhenotypes();
        const mode = document.getElementById("baseColorMode").value;
        const sections = [];

        sections.push(baseLegend(mode));

        if (phenotypes.size > 0) {
          const selectedItems = Array.from(phenotypes.entries())
            .map(([phenotype, color]) => legendItem(color, phenotype))
            .join("");
          sections.push(`
            <p class="legend-note">Phenotype overlay is active. The scaffold keeps the current base coloring mode, while selected phenotypes appear as colored residue markers.</p>
            <div class="legend">${selectedItems}</div>
          `);
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
            <div class="legend">
              ${legendItem("#b2182b", "Pathogenic")}
              ${legendItem("#ef8a62", "Likely pathogenic")}
              ${legendItem("#fddbc7", "Uncertain/conflicting")}
              ${legendItem("#d1e5f0", "Likely benign")}
              ${legendItem("#d9d9d9", "Benign or unannotated")}
            </div>
          `;
        }
        return `
            <p class="legend-note">Default 3Dmol spectrum coloring. Phenotype residue markers can be overlaid on top.</p>
            <div class="colorbar">
              <div class="colorbar-track" style="background: linear-gradient(90deg, #304ffe, #00bcd4, #4caf50, #ffeb3b, #f44336);"></div>
              <div class="colorbar-labels"><span>N terminus</span><span>C terminus</span></div>
            </div>
          `;
      }

    function legendItem(color, label) {
      return `<span class="legend-item"><span class="swatch" style="background:${color}"></span>${escapeHtml(label)}</span>`;
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
      const value = Number(document.getElementById("markerScale").value);
      return Math.max(0.5, Math.min(5.0, Number.isFinite(value) ? value : 1.6));
    }

    function updateBackgroundOpacityLabel() {
      document.getElementById("backgroundOpacityValue").textContent = selectedBackgroundOpacity().toFixed(2);
    }

      function updateMarkerScaleLabel() {
        const label = `${selectedMarkerScale().toFixed(1)}x`;
        document.getElementById("markerScaleValue").textContent = label;
        document.getElementById("markerScaleText").value = label;
      }

      function syncMarkerScaleFromText() {
        const input = document.getElementById("markerScaleText");
        const parsed = Number.parseFloat(String(input.value).replace(/x/gi, "").trim());
        const normalized = Math.max(0.5, Math.min(5.0, Number.isFinite(parsed) ? parsed : selectedMarkerScale()));
        document.getElementById("markerScale").value = normalized.toFixed(1);
        updateMarkerScaleLabel();
      }

    function resetView() {
      if (!viewer) {
        return;
      }
        document.getElementById("focusResidue").value = "";
        document.getElementById("backgroundOpacity").value = "0.3";
      document.getElementById("markerScale").value = "1.6";
        updateBackgroundOpacityLabel();
      updateMarkerScaleLabel();
      document.querySelectorAll(".phenotype-check").forEach((checkbox) => {
        checkbox.checked = false;
      });
      applyStyle(false);
      viewer.zoomTo();
      viewer.render();
    }

      function residueCenter(residue) {
        if (!model || typeof model.selectedAtoms !== "function") {
          return null;
        }
        const atoms = model.selectedAtoms({ resi: residue }) || [];
        if (atoms.length === 0) {
          return null;
        }
        const alphaCarbon = atoms.find((atom) => String(atom.atom || "").toUpperCase() === "CA");
        if (alphaCarbon) {
          return {
            x: Number(alphaCarbon.x) || 0,
            y: Number(alphaCarbon.y) || 0,
            z: Number(alphaCarbon.z) || 0,
          };
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
        return { x: x / count, y: y / count, z: z / count };
      }

      function focusCameraOnResidue(residue) {
        viewer.zoomTo({ resi: residue });
        if (typeof viewer.zoom === "function") {
          viewer.zoom(0.99, 250);
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
      const variantCount = Math.max(1, Number(annotation?.variant_count) || 1);
      const scale = selectedMarkerScale();
      const variantBoost = 1 + Math.log2(variantCount) * 0.25;
      return Math.max(0.16, baseRadius * variantBoost * scale);
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
      return String(value ?? "")
        .replaceAll("&", "&amp;")
        .replaceAll("<", "&lt;")
        .replaceAll(">", "&gt;")
        .replaceAll('"', "&quot;");
    }

    function setStatus(message) {
      document.getElementById("status").textContent = message;
    }

    window.addEventListener("load", init);
  </script>
</body>
</html>
"""
    replacements = {
        "__GENE__": escape(gene),
        "__STRUCTURE_LABEL__": escape(structure_label),
        "__REPORT_JSON__": _safe_json(report),
        "__PREPROCESSED_JSON__": _safe_json(annotations),
        "__PDB_JSON__": _safe_json(pdb_text),
    }
    for placeholder, value in replacements.items():
        template = template.replace(placeholder, value)
    return template


def _safe_json(value: Any) -> str:
    """Serialize data for embedding in a script block."""
    return json.dumps(value).replace("</", "<\\/")


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
) -> dict[str, Any]:
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
        "pathogenic_color": _pathogenicity_color(max_pathogenicity),
        "missense": missense_prediction,
        "alphafold_confidence": confidence_prediction,
        "variants": [_variant_summary(record) for record in records],
    }


def _variant_summary(record: dict[str, Any]) -> dict[str, Any]:
    return {
        "variation_id": record.get("variation_id"),
        "accession": record.get("accession"),
        "protein_change": record.get("protein_change"),
        "amino_acid": record.get("amino_acid"),
        "alternate_amino_acid": record.get("alternate_amino_acid"),
        "significance": record.get("significance"),
        "phenotype_names": record.get("phenotype_names", []),
        "molecular_consequences": record.get("molecular_consequences", []),
    }


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


def _pathogenicity_color(severity: int) -> str:
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
