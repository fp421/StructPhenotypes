#!/usr/bin/env python3
"""Render SCN1A patient variants (Dravet, GEFS+, gain-of-function) coloured by
category on the AlphaFold structure (AF-P35498) as a standalone 3Dmol.js page.

Sources (all mapped to canonical UniProt P35498 numbering):
  - Dravet / GEFS+ : data/patient_variants/SCN1A/dravet_gefs_brunklaus.csv
                     (uses residue_canonical; see structpheno/extract_brunklaus.py)
  - Gain-of-function: data/patient_variants/SCN1A/gof_ndeema.csv
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from collections import defaultdict
from pathlib import Path

BRUNKLAUS_CSV = Path("data/patient_variants/SCN1A/dravet_gefs_brunklaus.csv")
GOF_CSV = Path("data/patient_variants/SCN1A/gof_ndeema.csv")
ALPHAFOLD_DIR = Path("data/alphafold/SCN1A")

# Category display order + colours (red = most severe LOF, blue = GOF).
CATEGORIES = [
    ("Dravet", "#d62728"),
    ("GEFS+", "#ff7f0e"),
    ("GOF", "#1f77b4"),
]


def _residue_from_change(protein_change: str) -> int | None:
    """Parse a canonical residue number from a protein change like p.T162I."""
    m = re.match(r"p?\.?[A-Za-z]{1,3}(\d+)", protein_change)
    return int(m.group(1)) if m else None


def _is_missense(protein_change: str) -> bool:
    """True for a single-residue substitution (e.g. p.T162I); excludes
    nonsense (X), frameshifts (fs), and in-frame indels (del/dup/ins)."""
    return bool(re.fullmatch(r"p?\.?[A-Za-z]{1,3}\d+[A-Za-z]{1,3}", protein_change)) \
        and "fs" not in protein_change and "X" not in protein_change


def load_variant_categories() -> dict[str, dict[int, list[str]]]:
    """Return {category: {residue: [variant labels]}} for the three categories.

    Missense variants only — protein-truncating variants (nonsense, frameshift,
    whole-gene/exon deletions) are excluded, as they do not map to a single
    structural residue.
    """
    cats: dict[str, dict[int, list[str]]] = {name: defaultdict(list) for name, _ in CATEGORIES}

    if BRUNKLAUS_CSV.exists():
        for r in csv.DictReader(BRUNKLAUS_CSV.open()):
            dx = r["dx"]
            res = r.get("residue_canonical")
            if dx not in cats or not res or r.get("type") != "Missense":
                continue
            cats[dx][int(res)].append(r["protein_change"])
    else:
        print(f"Warning: {BRUNKLAUS_CSV} not found", file=sys.stderr)

    if GOF_CSV.exists():
        for r in csv.DictReader(GOF_CSV.open()):
            if not _is_missense(r["variant"]):
                continue
            res = _residue_from_change(r["variant"])
            if res is not None:
                cats["GOF"][res].append(r["variant"])
    else:
        print(f"Warning: {GOF_CSV} not found", file=sys.stderr)

    return {name: dict(sorted(residues.items())) for name, residues in cats.items()}


def find_pdb() -> Path:
    """Locate the SCN1A AlphaFold PDB."""
    candidates = sorted(ALPHAFOLD_DIR.glob("AF-*.pdb")) + sorted(ALPHAFOLD_DIR.glob("*.pdb"))
    if not candidates:
        sys.exit(
            f"No AlphaFold PDB in {ALPHAFOLD_DIR}. Run:\n"
            f"  python structpheno/get_alpha_fold.py --gene SCN1A"
        )
    return candidates[0]


def build_payload(cats: dict[str, dict[int, list[str]]]) -> dict:
    """Shape category data for the page, including per-residue hover info."""
    categories = []
    residue_info: dict[int, dict] = {}
    for name, color in CATEGORIES:
        residues = cats.get(name, {})
        categories.append({
            "name": name,
            "color": color,
            "residues": sorted(residues),
            "residue_count": len(residues),
            "variant_count": sum(len(v) for v in residues.values()),
        })
        for res, variants in residues.items():
            info = residue_info.setdefault(res, {"residue": res, "categories": []})
            info["categories"].append({"name": name, "color": color, "variants": variants})
    return {"categories": categories, "residue_info": residue_info}


def render_html(pdb_text: str, payload: dict) -> str:
    pdb_json = json.dumps(pdb_text).replace("</", "<\\/")
    data_json = json.dumps(payload).replace("</", "<\\/")

    legend_rows = "".join(
        f'<label class="cat"><input type="checkbox" data-cat="{c["name"]}" checked>'
        f'<span class="swatch" style="background:{c["color"]}"></span>'
        f'{c["name"]} <span class="muted">({c["residue_count"]} residues, {c["variant_count"]} variants)</span>'
        f'</label>'
        for c in payload["categories"]
    )

    return f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>SCN1A variants on AlphaFold (P35498)</title>
<script src="https://3Dmol.org/build/3Dmol-min.js"></script>
<style>
  body {{ margin:0; font-family:Arial, sans-serif; color:#1f2937; }}
  header {{ padding:16px 22px 10px; border-bottom:1px solid #d8dee8; }}
  h1 {{ margin:0 0 4px; font-size:20px; }}
  .subtitle {{ margin:0; color:#607085; font-size:13px; }}
  .toolbar {{ display:flex; flex-wrap:wrap; gap:18px; align-items:center;
             padding:12px 22px; background:#f7f9fc; border-bottom:1px solid #d8dee8; }}
  .cat {{ display:flex; align-items:center; gap:7px; font-size:14px; font-weight:600; cursor:pointer; }}
  .swatch {{ width:14px; height:14px; border-radius:3px; display:inline-block; }}
  .muted {{ color:#607085; font-weight:400; }}
  button {{ height:32px; padding:0 12px; border:1px solid #2459d6; background:#2459d6;
           color:#fff; border-radius:6px; font-weight:600; cursor:pointer; }}
  #viewer {{ width:100%; height:74vh; min-height:460px; position:relative; }}
  #status {{ padding:9px 22px; color:#607085; font-size:13px; border-top:1px solid #d8dee8; }}
</style>
</head>
<body>
<header>
  <h1>SCN1A variants on AlphaFold structure</h1>
  <p class="subtitle">AF-P35498 &middot; missense only &middot; colour = phenotype (overlaps take the most-severe: Dravet &gt; GEFS+ &gt; GOF) &middot; sphere size &prop; observations at that residue &middot; hover for details</p>
</header>
<div class="toolbar">
  {legend_rows}
  <label class="cat"><input type="checkbox" id="sizeToggle" checked> Size by count</label>
  <label class="cat"><input type="checkbox" id="surfaceToggle"> Surface</label>
  <button id="reset">Reset view</button>
</div>
<div id="viewer"></div>
<div id="status">Loading…</div>
<script>
  const PDB = {pdb_json};
  const DATA = {data_json};
  let viewer = null;

  function setStatus(msg) {{ document.getElementById("status").textContent = msg; }}

  function enabledCategories() {{
    return [...document.querySelectorAll('input[data-cat]')]
      .filter(cb => cb.checked).map(cb => cb.dataset.cat);
  }}

  // Severity order for resolving overlaps (index 0 = most severe).
  function priorityIndex(name) {{ return DATA.categories.findIndex(c => c.name === name); }}

  function applyStyles() {{
    if (!viewer) return;
    viewer.setStyle({{}}, {{ cartoon: {{ color: "#d9dee6" }} }});
    const enabled = enabledCategories();
    const sizeByCount = document.getElementById("sizeToggle").checked;
    let shown = 0, maxCount = 0;
    // One sphere per residue: most-severe enabled category sets the colour,
    // radius grows with the number of patient observations at that residue.
    for (const key in DATA.residue_info) {{
      const info = DATA.residue_info[key];
      const present = info.categories.filter(c => enabled.includes(c.name));
      if (!present.length) continue;
      present.sort((a, b) => priorityIndex(a.name) - priorityIndex(b.name));
      const color = present[0].color;
      const count = present.reduce((n, c) => n + c.variants.length, 0);
      maxCount = Math.max(maxCount, count);
      const radius = sizeByCount ? (1.0 + 0.55 * Math.sqrt(count)) : 1.6;
      viewer.addStyle({{ resi: info.residue }}, {{ cartoon: {{ color }} }});
      viewer.addStyle({{ resi: info.residue, atom: "CA" }}, {{ sphere: {{ color, radius }} }});
      shown++;
    }}
    viewer.render();
    setStatus(`Showing ${{enabled.join(", ") || "no categories"}} — ${{shown}} residues`
      + (sizeByCount ? ` · sphere size ∝ √(observations), busiest = ${{maxCount}}.` : "."));
  }}

  function setupHover() {{
    viewer.setHoverable({{}}, true,
      function(atom) {{
        if (!atom) return;
        const info = DATA.residue_info[atom.resi];
        if (!info) return;
        const lines = info.categories.map(c =>
          `${{c.name}}: ${{c.variants.join(", ")}}`).join("\\n");
        viewer.addLabel(`Residue ${{atom.resi}}\\n${{lines}}`,
          {{ position: {{ x: atom.x, y: atom.y, z: atom.z }},
            backgroundColor: "black", backgroundOpacity: 0.8, fontSize: 12 }});
        viewer.render();
      }},
      function() {{ viewer.removeAllLabels(); viewer.render(); }});
  }}

  function init() {{
    if (!window.$3Dmol) {{ setStatus("3Dmol.js failed to load (check connection)."); return; }}
    viewer = $3Dmol.createViewer("viewer", {{ backgroundColor: "white" }});
    viewer.addModel(PDB, "pdb");
    applyStyles();
    setupHover();
    viewer.zoomTo();
    viewer.render();

    document.querySelectorAll('input[data-cat]').forEach(cb =>
      cb.addEventListener("change", applyStyles));
    document.getElementById("sizeToggle").addEventListener("change", applyStyles);
    document.getElementById("reset").addEventListener("click", () => {{
      viewer.zoomTo(); viewer.render();
    }});
    document.getElementById("surfaceToggle").addEventListener("change", e => {{
      viewer.removeAllSurfaces();
      if (e.target.checked) {{
        viewer.addSurface($3Dmol.SurfaceType.VDW, {{ opacity: 0.55, color: "white" }});
      }}
      viewer.render();
    }});
  }}
  window.addEventListener("load", init);
</script>
</body>
</html>
"""


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-o", "--output", default="outputs/scn1a_variants.html",
                        help="Output HTML path")
    parser.add_argument("--pdb", help="AlphaFold PDB path (default: auto-detect in data/alphafold/SCN1A)")
    args = parser.parse_args()

    cats = load_variant_categories()
    pdb_path = Path(args.pdb) if args.pdb else find_pdb()
    payload = build_payload(cats)

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(render_html(pdb_path.read_text(), payload), encoding="utf-8")

    print(f"Structure: {pdb_path}")
    for c in payload["categories"]:
        print(f"  {c['name']:7}: {c['residue_count']:4} residues, {c['variant_count']} variants")
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
