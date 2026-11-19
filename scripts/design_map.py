"""Design map for low-cost fab stackups (JLCPCB 4-layer 1.6 mm, data/published/jlcpcb_stackups.csv)
with a manufacturing-tolerance Monte Carlo.

Physics shortcut, verified by the sweep (results/tables/epc_sweep*.csv): optimal and
lateral loops depend on the L1-L2 dielectric d1 only (board thickness T has < 0.5 %
effect); the vertical loop depends on T only.  So per (family, stackup) we run
FastHenry at the converged HF settings at 0.9x, 1.0x, 1.1x of the sensitive variable,
fit a quadratic surrogate and sample it.

Tolerance model (ASSUMPTION, stated in the paper): JLCPCB states +-10 % on finished
board thickness and nothing per ply; we take d1 and T independently uniform in +-10 %.
Copper per JLC: outer 0.035 mm, inner 0.0152 mm.
Output: results/tables/design_map.csv, results/tables/design_map_points.csv
"""
from __future__ import annotations

import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
import fhrun  # noqa: E402
import geom  # noqa: E402
import layouts  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "results" / "raw" / "design_map"
TAB = ROOT / "results" / "tables"
PITCH, FW, NH = 0.15, 0.2, 5
FAMS = ["optimal", "lateral_stitched", "lateral", "vertical"]
SCALES = (0.9, 1.0, 1.1)


def job(a):
    fam, name, T, d1, s = a
    if fam == "vertical":
        T = T * s
    else:
        d1 = d1 * s
    p = layouts.Params(family=fam, T=T, d1=d1, cu=0.035, cu_in=0.0152, pitch=PITCH)
    r = fhrun.run(geom.to_fasthenry(layouts.build(p), (1e6, 1e8), max_fil_w=FW, nh=NH), RAW / p.tag())
    return dict(family=fam, stackup=name, scale=s, T=T, d1=d1, L_1MHz_nH=r["L"][0] * 1e9,
                L_100MHz_nH=r["L"][-1] * 1e9, R_100MHz_mohm=r["R"][-1] * 1e3, runtime_s=r["runtime_s"])


def main():
    st = pd.read_csv(ROOT / "data" / "published" / "jlcpcb_stackups.csv", comment="#")
    jobs = []
    for _, s in st.iterrows():
        T = s.total_derived_mm - s.cu_outer_mm          # L1-to-L4 copper centre distance
        for fam in FAMS:
            for k in SCALES:
                jobs.append((fam, s["name"], T, s.d1_mm, k))
    with ProcessPoolExecutor(max_workers=3) as ex:
        pts = pd.DataFrame(list(ex.map(job, jobs)))
    pts.to_csv(TAB / "design_map_points.csv", index=False)
    rng = np.random.default_rng(7)
    rows = []
    for (fam, name), g in pts.groupby(["family", "stackup"], sort=False):
        g = g.sort_values("scale")
        c = np.polyfit(g.scale, g.L_100MHz_nH, 2)
        u = rng.uniform(0.9, 1.1, 100_000)
        Ls = np.polyval(c, u)
        rows.append(dict(family=fam, stackup=name.replace("JLC04161H-", ""), d1_mm=g.d1[g.scale == 1.0].iloc[0],
                         L_nom_nH=g.L_100MHz_nH[g.scale == 1.0].iloc[0], L_mean_nH=Ls.mean(),
                         L_p2_5_nH=np.percentile(Ls, 2.5), L_p97_5_nH=np.percentile(Ls, 97.5),
                         sens_pct_per_10pct=100 * (g.L_100MHz_nH.iloc[2] - g.L_100MHz_nH.iloc[0]) / 2
                         / g.L_100MHz_nH.iloc[1]))
    d = pd.DataFrame(rows)
    d.to_csv(TAB / "design_map.csv", index=False)
    print(d.round(3).to_string())


if __name__ == "__main__":
    main()
