"""Build the analysis tables from the pixel digitisation (scripts/digitise_epc.py).

1) data/published/epc_measured_boards.csv: every measured marker (slides 15-17).
2) data/published/epc_fig7b_digitised.csv: EPC simulated L(board thickness) per family and
   inner-layer distance, extrapolated to exactly 31 and 62 mil along each straight line.
   Style assignment: at a given x the lower of two same-colour lines is the thinner inner
   distance (4 mil), the next 12 mil (monotonic in d1, WP010 text); the 26-mil dotted line is
   the third cluster at 57.5 mil.  The two vertical-loop lines nearly coincide; their two
   values at 31.5 mil are averaged (difference reported).
"""
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
d = pd.read_csv(ROOT / "data" / "published" / "digitised_apec13.csv")

m = d[d.slide.isin([15, 16, 17])].copy()
m["source"] = "Reusch APEC2013 slide " + m.slide.astype(str) + " (pixel-digitised, scripts/digitise_epc.py)"
m[["family", "L_epc_nH", "quantity", "value", "occluded", "slide", "source"]].to_csv(
    ROOT / "data" / "published" / "epc_measured_boards.csv", index=False)

s = d[d.slide == 14]
rows = []
for fam in ("optimal", "lateral", "vertical"):
    f = s[s.family == fam]
    a = sorted(f[f.board_thickness_mil == 31.5].value)
    b = sorted(f[f.board_thickness_mil == 61.0].value)
    c = sorted(f[f.board_thickness_mil == 57.5].value)
    if fam == "lateral" and len(a) == 3:          # solid green split in two by the crossing red line
        a = [np.mean(a[:2]), a[2]]
    if fam == "vertical":
        series = {4: (np.mean(a), b[0], abs(a[1] - a[0]))}
    else:
        series = {4: (a[0], b[0], 0.0), 12: (a[1], b[1], 0.0)}
    for d1, (y1, y2, spread) in series.items():
        slope = (y2 - y1) / (61.0 - 31.5)
        for T, y in ((31, y1 + slope * (31 - 31.5)), (62, y2 + slope * (62 - 61.0))):
            rows.append(dict(family=fam, inner_layer_mil=d1, board_thickness_mil=T, L_sim_nH=round(y, 3),
                             note=f"two near-coincident lines, spread {spread:.3f} nH" if spread else ""))
    if fam != "vertical":
        top = c[-1]                                # dotted 26-mil curve: highest cluster at 57.5 mil
        rows.append(dict(family=fam, inner_layer_mil=26, board_thickness_mil=62, L_sim_nH=round(top, 3),
                         note="dotted curve, nearly flat; value at 57.5 mil"))
out = pd.DataFrame(rows)
out["source"] = "EPC WP010 Fig. 7b = Reusch APEC2013 slide 14 (pixel-digitised, +-0.003 nH/px)"
out.to_csv(ROOT / "data" / "published" / "epc_fig7b_digitised.csv", index=False)
print(out.to_string())
