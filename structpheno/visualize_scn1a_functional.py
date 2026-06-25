#!/usr/bin/env python3
"""Render SCN1A funNCion functional predictions on the AlphaFold structure.

Source: data/functional_variants/SCN1A/bosselman.csv — the funNCion tool's
gain-of-function (GOF) confidence per variant (column ``funNCion_prediction``,
e.g. "GOF (83%)"). All entries are GOF, so residues are coloured on a single
sequential confidence ramp; training-cohort variants without a funNCion score
are shown in grey. Residue numbering is canonical UniProt P35498.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from collections import defaultdict
from pathlib import Path

try:
    from .visualize_scn1a_variants import find_pdb, _residue_from_change
except ImportError:
    from visualize_scn1a_variants import find_pdb, _residue_from_change

FUNCTIONAL_CSV = Path("data/functional_variants/SCN1A/bosselman.csv")

# Confidence ramp endpoints (GOF): light at low confidence -> dark red at high.
RAMP_LOW_PCT, RAMP_HIGH_PCT = 65, 100
RAMP_LOW_RGB, RAMP_HIGH_RGB = (255, 237, 160), (127, 0, 0)  # YlOrRd-style
NO_SCORE_COLOR = "#9aa3ad"  # training-cohort / unscored GOF variants


def parse_funncion(value: str) -> tuple[str | None, int | None]:
    """Parse 'GOF (83%)' -> ('GOF', 83); training/NA rows -> (None, None)."""
    m = re.match(r"(GOF|LOF|Mixed|Neutral)\s*\((\d+)%\)", value or "")
    return (m.group(1), int(m.group(2))) if m else (None, None)


def _score_color(pct: int) -> str:
    """Interpolate the confidence ramp for a percentage."""
    t = (pct - RAMP_LOW_PCT) / (RAMP_HIGH_PCT - RAMP_LOW_PCT)
    t = max(0.0, min(1.0, t))
    rgb = tuple(round(lo + (hi - lo) * t) for lo, hi in zip(RAMP_LOW_RGB, RAMP_HIGH_RGB))
    return "#{:02x}{:02x}{:02x}".format(*rgb)


def load_functional() -> dict[int, dict]:
    """Aggregate funNCion predictions by canonical residue."""
    if not FUNCTIONAL_CSV.exists():
        sys.exit(f"{FUNCTIONAL_CSV} not found")

    by_res: dict[int, dict] = {}
    for r in csv.DictReader(FUNCTIONAL_CSV.open()):
        res = _residue_from_change(r["protein_change"])
        if res is None:
            continue
        call, pct = parse_funncion(r["funNCion_prediction"])
        entry = by_res.setdefault(res, {"residue": res, "variants": [], "score": None, "call": call})
        entry["variants"].append({
            "variant": r["protein_change"],
            "funncion": r["funNCion_prediction"],
            "function_call": r["function_call"],
        })
        if pct is not None and (entry["score"] is None or pct > entry["score"]):
            entry["score"] = pct
            entry["call"] = call
    return by_res


def build_payload(by_res: dict[int, dict]) -> dict:
    residues = []
    for res, info in sorted(by_res.items()):
        residues.append({
            "residue": res,
            "score": info["score"],
            "call": info["call"],
            "color": _score_color(info["score"]) if info["score"] is not None else NO_SCORE_COLOR,
            "variants": info["variants"],
        })
    scored = [r for r in residues if r["score"] is not None]
    return {
        "residues": residues,
        "scored_count": len(scored),
        "unscored_count": len(residues) - len(scored),
        "ramp": {"low": RAMP_LOW_PCT, "high": RAMP_HIGH_PCT,
                 "low_color": _score_color(RAMP_LOW_PCT), "high_color": _score_color(RAMP_HIGH_PCT),
                 "no_score": NO_SCORE_COLOR},
    }


def render_html(pdb_text: str, payload: dict) -> str:
    pdb_json = json.dumps(pdb_text).replace("</", "<\\/")
    data_json = json.dumps(payload).replace("</", "<\\/")
    ramp = payload["ramp"]
    gradient = (f"linear-gradient(to right, {ramp['low_color']}, {ramp['high_color']})")

    return f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>SCN1A funNCion GOF scores on AlphaFold (P35498)</title>
<script src="https://3Dmol.org/build/3Dmol-min.js"></script>
<style>
  body {{ margin:0; font-family:Arial, sans-serif; color:#1f2937; }}
  header {{ padding:16px 22px 10px; border-bottom:1px solid #d8dee8; }}
  h1 {{ margin:0 0 4px; font-size:20px; }}
  .subtitle {{ margin:0; color:#607085; font-size:13px; }}
  .toolbar {{ display:flex; flex-wrap:wrap; gap:18px; align-items:center;
             padding:12px 22px; background:#f7f9fc; border-bottom:1px solid #d8dee8; font-size:14px; }}
  .ramp {{ display:flex; align-items:center; gap:8px; }}
  .bar {{ width:160px; height:14px; border-radius:3px; background:{gradient}; border:1px solid #c3c9d2; }}
  .swatch {{ width:14px; height:14px; border-radius:3px; display:inline-block; border:1px solid #c3c9d2; }}
  .ctl {{ display:flex; align-items:center; gap:7px; font-weight:600; cursor:pointer; }}
  .muted {{ color:#607085; font-weight:400; }}
  button {{ height:32px; padding:0 12px; border:1px solid #2459d6; background:#2459d6;
           color:#fff; border-radius:6px; font-weight:600; cursor:pointer; }}
  #viewer {{ width:100%; height:74vh; min-height:460px; position:relative; }}
  #status {{ padding:9px 22px; color:#607085; font-size:13px; border-top:1px solid #d8dee8; }}
</style>
</head>
<body>
<header>
  <h1>SCN1A funNCion functional scores</h1>
  <p class="subtitle">AF-P35498 &middot; gain-of-function confidence from the funNCion predictor &middot; sphere size &prop; score &middot; hover for details</p>
</header>
<div class="toolbar">
  <div class="ramp"><span class="muted">GOF confidence</span>
    <span class="muted">{ramp['low']}%</span><span class="bar"></span><span class="muted">{ramp['high']}%</span></div>
  <span><span class="swatch" style="background:{ramp['no_score']}"></span> no funNCion score (training)</span>
  <label class="ctl"><input type="checkbox" id="unscoredToggle" checked> Show unscored</label>
  <label class="ctl"><input type="checkbox" id="surfaceToggle"> Surface</label>
  <button id="reset">Reset view</button>
</div>
<div id="viewer"></div>
<div id="status">Loading…</div>
<script>
  const PDB = {pdb_json};
  const DATA = {data_json};
  let viewer = null;
  function setStatus(m) {{ document.getElementById("status").textContent = m; }}

  function applyStyles() {{
    if (!viewer) return;
    viewer.setStyle({{}}, {{ cartoon: {{ color: "#d9dee6" }} }});
    const showUnscored = document.getElementById("unscoredToggle").checked;
    let shown = 0;
    DATA.residues.forEach(r => {{
      if (r.score === null && !showUnscored) return;
      // size by score (scored), fixed small radius for unscored
      const radius = r.score === null ? 1.3 : 1.2 + 2.2 * ((r.score - {RAMP_LOW_PCT}) / ({RAMP_HIGH_PCT} - {RAMP_LOW_PCT}));
      viewer.addStyle({{ resi: r.residue }}, {{ cartoon: {{ color: r.color }} }});
      viewer.addStyle({{ resi: r.residue, atom: "CA" }}, {{ sphere: {{ color: r.color, radius: Math.max(1.1, radius) }} }});
      shown++;
    }});
    viewer.render();
    setStatus(`${{shown}} residues shown — ${{DATA.scored_count}} with funNCion GOF score, ${{DATA.unscored_count}} unscored.`);
  }}

  function setupHover() {{
    const byRes = {{}};
    DATA.residues.forEach(r => byRes[r.residue] = r);
    viewer.setHoverable({{}}, true,
      function(atom) {{
        const r = byRes[atom.resi];
        if (!r) return;
        const lines = r.variants.map(v => `${{v.variant}}: ${{v.funncion}}`).join("\\n");
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
    document.getElementById("unscoredToggle").addEventListener("change", applyStyles);
    document.getElementById("reset").addEventListener("click", () => {{ viewer.zoomTo(); viewer.render(); }});
    document.getElementById("surfaceToggle").addEventListener("change", e => {{
      viewer.removeAllSurfaces();
      if (e.target.checked) viewer.addSurface($3Dmol.SurfaceType.VDW, {{ opacity: 0.55, color: "white" }});
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
    parser.add_argument("-o", "--output", default="outputs/scn1a_functional.html",
                        help="Output HTML path")
    parser.add_argument("--pdb", help="AlphaFold PDB path (default: auto-detect)")
    args = parser.parse_args()

    by_res = load_functional()
    payload = build_payload(by_res)
    pdb_path = Path(args.pdb) if args.pdb else find_pdb()

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(render_html(pdb_path.read_text(), payload), encoding="utf-8")

    print(f"Structure: {pdb_path}")
    print(f"  funNCion-scored residues : {payload['scored_count']}")
    print(f"  unscored (training) GOF  : {payload['unscored_count']}")
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
