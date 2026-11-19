"""Tier-1 verification of FastHenry on canonical geometries.

V1  straight bar: partial self-L vs Grover/Ruehli closed form.
V2  per-unit-length loop inductance of two anti-parallel plates vs the exact 2D
    GMD solution (end effects removed by differencing two lengths).
V3  full 3D two-plate loops published by Letellier et al. (IEEE TPEL 2018) with
    their commercial-FEA (Infolytica MagNet) values.
Each case is run on a filament-refinement ladder (mesh convergence).

Outputs: results/tables/canonical_fh.csv, results/raw/canonical_fh/<job>/...
"""
from __future__ import annotations

import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
import analytic as an  # noqa: E402
import fhrun  # noqa: E402
import geom  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "results" / "raw" / "canonical_fh"
TAB = ROOT / "results" / "tables"

# filament width ladder (mm) and thickness filaments
LADDER = [(0.5, 1), (0.25, 1), (0.125, 2), (0.0625, 3)]
FREQS = (1e3, 1e8)  # decade sweep: 1 kHz ~ DC check, 1 MHz = Letellier FEA freq, 100 MHz ~ ringing band

# Letellier 2018 published cases: (id, w, l, e, FEA nH, Eq13 nH, source location)
LETELLIER = [
    ("L18_EPC2014_inner", 1.7, 9.5, 1.0, 4.31, 4.67, "Table 1, col 1"),
    ("L18_GS61008_inner", 7.0, 15.5, 1.0, 2.28, 2.52, "Table 1, col 2"),
    ("L18_EPC2014_gate", 0.254, 7.2, 1.0, 8.81, 7.82, "Table 1, col 3"),
    ("L18_EPC2014_topbot", 1.7, 9.5, 1.6, 5.72, 6.20, "Fig. 14 table"),
    ("L18_Fig6_inner", 7.5, 19.0, 0.53, 1.35, None, "Sec. II.C / Fig. 6 (e=0.53 text, 0.59 caption)"),
    ("L18_Fig6_inner_e059", 7.5, 19.0, 0.59, 1.35, None, "Fig. 6 caption variant e=0.59"),
]


def _job(args):
    jid, inp, wd = args
    try:
        r = fhrun.run(inp, wd)
        return jid, {f: (L, R) for f, L, R in zip(r["freqs"], r["L"], r["R"])}, r["runtime_s"], ""
    except Exception as ex:  # keep going, report
        return jid, {}, 0.0, str(ex)


def main():
    jobs, meta = [], {}
    # V1 straight bar
    for fw, nh in LADDER:
        g = geom.straight_bar(10.0, 1.0, 0.035)
        jid = f"V1_bar_fw{fw}_nh{nh}"
        jobs.append((jid, geom.to_fasthenry(g, FREQS, max_fil_w=fw, nh=nh), RAW / jid))
        meta[jid] = dict(tier="V1", case="bar_10x1x0.035", fw=fw, nh=nh,
                         ref_nH=an.bar_partial_self_L(10e-3, 1e-3, 35e-6) * 1e9, ref="Grover/Ruehli")
    # V2 per-unit-length, difference of l=40 and l=20 mm
    for (cid, w, _l, e, *_rest) in LETELLIER[:4]:
        for fw, nh in LADDER:
            for l in (20.0, 40.0):
                g = geom.two_plate_loop(w, l, e)
                jid = f"V2_{cid}_l{int(l)}_fw{fw}_nh{nh}"
                jobs.append((jid, geom.to_fasthenry(g, FREQS, max_fil_w=fw, nh=nh), RAW / jid))
                meta[jid] = dict(tier="V2", case=cid, fw=fw, nh=nh, l=l, w=w, e=e)
    # V3 published full loops
    for (cid, w, l, e, fea, eq13, src) in LETELLIER:
        for fw, nh in LADDER:
            g = geom.two_plate_loop(w, l, e)
            jid = f"V3_{cid}_fw{fw}_nh{nh}"
            jobs.append((jid, geom.to_fasthenry(g, FREQS, max_fil_w=fw, nh=nh), RAW / jid))
            meta[jid] = dict(tier="V3", case=cid, fw=fw, nh=nh, w=w, l=l, e=e,
                             ref_nH=fea, ref="Letellier2018 MagNet FEA", eq13_nH=eq13, src=src,
                             exact2d_nH=an.two_bar_loop_L(l * 1e-3, w * 1e-3, 35e-6, e * 1e-3) * 1e9)

    rows = []
    with ProcessPoolExecutor(max_workers=8) as ex:
        for jid, res, rt, err in ex.map(_job, jobs):
            m = dict(meta[jid], job=jid, runtime_s=rt, error=err)
            for f, (L, R) in res.items():
                rows.append(dict(m, f_Hz=f, L_nH=L * 1e9, R_mohm=R * 1e3))
            if not res:
                rows.append(m)
                print("FAILED", jid, err)
    df = pd.DataFrame(rows)
    TAB.mkdir(parents=True, exist_ok=True)
    df.to_csv(TAB / "canonical_fh_raw.csv", index=False)

    # ---- summaries
    out = []
    v1 = df[df.tier == "V1"]
    for _, r in v1.iterrows():
        out.append(dict(tier="V1", case=r.case, fw=r.fw, nh=r.nh, f_Hz=r.f_Hz, fh_nH=r.L_nH,
                        ref_nH=r.ref_nH, ref=r.ref))
    v2 = df[df.tier == "V2"]
    for (case, fw, nh, f), grp in v2.groupby(["case", "fw", "nh", "f_Hz"]):
        L20 = grp[grp.l == 20.0].L_nH.iloc[0]
        L40 = grp[grp.l == 40.0].L_nH.iloc[0]
        w, e = grp.w.iloc[0], grp.e.iloc[0]
        per_mm = (L40 - L20) / 20.0
        ref = an.two_bar_loop_L_per_m(w * 1e-3, 35e-6, e * 1e-3) * 1e9 * 1e-3
        out.append(dict(tier="V2", case=case, fw=fw, nh=nh, f_Hz=f, fh_nH=per_mm, ref_nH=ref,
                        ref="exact 2D GMD (nH/mm)"))
    v3 = df[df.tier == "V3"]
    for _, r in v3.iterrows():
        out.append(dict(tier="V3", case=r.case, fw=r.fw, nh=r.nh, f_Hz=r.f_Hz, fh_nH=r.L_nH,
                        ref_nH=r.ref_nH, ref=r.ref, eq13_nH=r.eq13_nH, exact2d_nH=r.exact2d_nH,
                        src=r.src))
    s = pd.DataFrame(out)
    s["err_pct"] = 100 * (s.fh_nH - s.ref_nH) / s.ref_nH
    s.to_csv(TAB / "canonical_fh.csv", index=False)
    pd.set_option("display.width", 200)
    print(s.round(4).to_string())


if __name__ == "__main__":
    main()
