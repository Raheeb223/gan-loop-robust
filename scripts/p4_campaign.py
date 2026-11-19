"""Phase 4/5 FastHenry campaign: one space-filling set of runs over design AND tolerance
variables of the optimal power loop on a JLCPCB 1.6 mm four-layer board.

Plain language: instead of running the slow solver inside the uncertainty analysis and
the optimiser, we run it a few hundred times at well-spread input combinations once,
then train a fast statistical model (surrogate) on those runs (p4_uq.py).

Variables (ranges and sources)
  design      d1        L1-L2 prepreg 0.07-0.22 mm  (covers JLCPCB 1080/3313/2116/7628 = 76/99/116/210 um)
              n_via     vias per row 3-7            (7 = densest row with >= 0.28 mm hole-to-hole gap at 0.3 mm drill)
              l_vin, l_sw, l_gnd   copper lengths 0.8-1.8 / 0.6-1.6 / 0.8-1.8 mm
  tolerance   t_d1      x[0.9, 1.1]  finished-thickness tolerance +-10 % (JLCPCB capabilities page) applied to d1
              t_cu      x[0.9, 1.1]  outer copper thickness (ASSUMPTION; JLCPCB states no copper tolerance)
              t_cuin    x[0.9, 1.1]  inner copper thickness (ASSUMPTION)
              hole      0.22-0.43 mm finished hole for a 0.3 mm drill (JLCPCB PTH tolerance +0.13/-0.08 mm);
                        via bar width = 0.25 mm x hole / 0.3 mm
              bump      0.090-0.110 mm  +-10 % seating height (EPC2218 datasheet states 120 +- 12 um; same
                        relative tolerance assumed for EPC2015C)
              dx_die    +-0.075 mm  die placement offset (ASSUMPTION, set equal to the JLCPCB hole-position
                        tolerance +-0.075 mm as an order of magnitude)
Sets: 'train' (256, scrambled Sobol), 'cal' (64) and 'test' (64) independent LHS sets for
split-conformal calibration and coverage testing, 'mc_nominal' (40) pure tolerance Monte
Carlo at the nominal design (direct check of the surrogate-based uncertainty).
Fixed seeds; resume-safe.  Output: results/phase4/campaign.csv
Runtime: ~424 runs x 20-40 s / 4 workers ~ 40-60 min.
"""
from __future__ import annotations

import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import qmc

sys.path.insert(0, str(Path(__file__).parent))
import fhrun  # noqa: E402
import geom  # noqa: E402
import layouts  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "phase4"
RAW = ROOT / "results" / "raw" / "phase4"
PITCH, FW, NH = 0.15, 0.2, 5
CU, CU_IN, T_JLC = 0.035, 0.0152, 1.6 - 0.035
VARS = ["d1", "n_via", "l_vin", "l_sw", "l_gnd", "t_d1", "t_cu", "t_cuin", "hole", "bump", "dx_die"]
LO = np.array([0.07, 2.5, 0.8, 0.6, 0.8, 0.9, 0.9, 0.9, 0.22, 0.090, -0.075])
HI = np.array([0.22, 7.5, 1.8, 1.6, 1.8, 1.1, 1.1, 1.1, 0.43, 0.110, 0.075])
NOMINAL = dict(d1=0.099, n_via=5, l_vin=1.2, l_sw=1.0, l_gnd=1.0)
TOL = VARS[5:]


def params(x: dict) -> layouts.Params:
    return layouts.Params(family="optimal", T=T_JLC, d1=x["d1"] * x["t_d1"], cu=CU * x["t_cu"],
                          cu_in=CU_IN * x["t_cuin"], n_via=int(round(x["n_via"])), l_vin=x["l_vin"],
                          l_sw=x["l_sw"], l_gnd=x["l_gnd"], via_w=0.25 * x["hole"] / 0.3, bump=x["bump"],
                          dx_die=x["dx_die"], pitch=PITCH)


def run(a):
    set_, i, x = a
    p = params(x)
    try:
        r = fhrun.run(geom.to_fasthenry(layouts.build(p), (1e8, 1e8), max_fil_w=FW, nh=NH),
                      RAW / f"{set_}_{i:04d}", timeout=1800)
        return dict(set=set_, idx=i, **x, L_nH=r["L"][0] * 1e9, R_mohm=r["R"][0] * 1e3, runtime_s=r["runtime_s"],
                    err="")
    except Exception as ex:
        return dict(set=set_, idx=i, **x, err=str(ex)[:200])


def design():
    rows = []
    U = qmc.Sobol(len(VARS), scramble=True, seed=11).random(256)
    rows += [("train", i, u) for i, u in enumerate(U)]
    rows += [("cal", i, u) for i, u in enumerate(qmc.LatinHypercube(len(VARS), seed=12).random(64))]
    rows += [("test", i, u) for i, u in enumerate(qmc.LatinHypercube(len(VARS), seed=13).random(64))]
    out = []
    for s, i, u in rows:
        x = LO + u * (HI - LO)
        d = dict(zip(VARS, x))
        d["n_via"] = float(int(round(d["n_via"])))
        out.append((s, i, d))
    U = qmc.LatinHypercube(len(TOL), seed=14).random(40)
    for i, u in enumerate(U):
        d = dict(NOMINAL)
        d.update(zip(TOL, LO[5:] + u * (HI[5:] - LO[5:])))
        out.append(("mc_nominal", i, d))
    return out


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    f = OUT / "campaign.csv"
    prev = pd.read_csv(f) if f.exists() else pd.DataFrame()
    done = set() if prev.empty else {(r.set, r.idx) for r in prev[prev.err.fillna("") == ""].itertuples()}
    jobs = [j for j in design() if (j[0], j[1]) not in done]
    rows = [] if prev.empty else prev[prev.err.fillna("") == ""].to_dict("records")
    print("jobs to run:", len(jobs), flush=True)
    with ProcessPoolExecutor(max_workers=4) as ex:
        for k, r in enumerate(ex.map(run, jobs)):
            rows.append(r)
            if k % 20 == 0 or k == len(jobs) - 1:
                pd.DataFrame(rows).to_csv(f, index=False)
                print(k + 1, "/", len(jobs), r.get("L_nH"), r.get("runtime_s"), r.get("err"), flush=True)
    pd.DataFrame(rows).to_csv(f, index=False)


if __name__ == "__main__":
    main()
