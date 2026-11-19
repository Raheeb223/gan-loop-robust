"""Tier-1b: openEMS (FDTD) vs FastHenry (PEEC) on the canonical Letellier loops.

Like-for-like rule: openEMS uses zero-thickness PEC sheets at the copper
centre-planes.  The FastHenry twin uses the *same centre-planes* with a 2 um
thin conductor evaluated at 100 MHz and 1 GHz (approaching the PEC / skin-limit).
The separate effect of real 35 um copper thickness is reported from FastHenry.

Output: results/tables/canonical_oems.csv
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
import fhrun  # noqa: E402
import geom  # noqa: E402
import oems  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "results" / "raw" / "canonical_oems"
TAB = ROOT / "results" / "tables"
H_CU = 0.035
H_THIN = 0.002

CASES = [  # id, w, l, e(gap between 35um copper faces), MagNet FEA nH
    ("L18_EPC2014_inner", 1.7, 9.5, 1.0, 4.31),
    ("L18_GS61008_inner", 7.0, 15.5, 1.0, 2.28),
    ("L18_EPC2014_topbot", 1.7, 9.5, 1.6, 5.72),
    ("L18_Fig6_inner", 7.5, 19.0, 0.53, 1.35),
]
LADDER = [0.4, 0.2, 0.1]
F_EVAL = np.linspace(20e6, 200e6, 19)


def main(only_fh: bool = False):
    rows = []
    out = TAB / "canonical_oems.csv"
    done = set()
    if out.exists():
        prev = pd.read_csv(out)
        rows = prev.to_dict("records")
        done = {(r["case"], r["solver"], r["mesh"]) for r in rows}
    for cid, w, l, e, fea in CASES:
        e_c = e + H_CU                       # centre-line spacing
        # FastHenry twin: thin conductor at identical centre-lines
        if (cid, "fasthenry_thin", 0.0625) not in done:
            g_thin = geom.two_plate_loop(w, l, e_c - H_THIN, h=H_THIN)
            r = fhrun.run(geom.to_fasthenry(g_thin, (1e6, 1e9), max_fil_w=0.0625, nh=1),
                          RAW / f"{cid}_fh_thin")
            for f, L in zip(r["freqs"], r["L"]):
                if f in (1e6, 1e8, 1e9):
                    rows.append(dict(case=cid, solver="fasthenry_thin", mesh=0.0625, f_Hz=f, L_nH=L * 1e9,
                                     fea_nH=fea, runtime_s=r["runtime_s"]))
        if only_fh:
            continue
        for fine in LADDER:
            if (cid, "openems_sheet", fine) in done:
                continue
            g = geom.two_plate_loop(w, l, e)
            res = oems.run_geometry(g, RAW / f"{cid}_oems_{fine}", f_max=2e9, f_eval=F_EVAL,
                                    fine=fine, margin=20, sheet=True)
            L0, fr, resid = oems.fit_L0(F_EVAL, res["L"], f_hi=200e6)
            rows.append(dict(case=cid, solver="openems_sheet", mesh=fine, f_Hz=0.0, L_nH=L0 * 1e9,
                             L_app_20MHz_nH=res["L"][0] * 1e9, L_app_100MHz_nH=float(np.interp(1e8, F_EVAL, res["L"])) * 1e9,
                             fr_GHz=fr / 1e9, fit_resid=resid, fea_nH=fea, runtime_s=res["runtime_s"],
                             ncells=res["ncells"]))
            pd.DataFrame(rows).to_csv(out, index=False)   # resume-safe
            print(cid, fine, f"L0={L0*1e9:.4f} nH  t={res['runtime_s']:.0f}s", flush=True)
    pd.DataFrame(rows).to_csv(out, index=False)


if __name__ == "__main__":
    main(only_fh="--fh-only" in sys.argv)
