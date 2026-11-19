"""End-to-end validation: open-tool loop inductance -> ngspice -> measured EPC switching data.

Data: EPC/Reusch boards (3 families x designs D1 31/4, D2 31/12, D3 62/4, D4 62/26 mil),
measured switch-node overshoot and dv/dt, pixel-digitised from APEC 2013 slides 16-17
(scripts/digitise_epc.py -> data/published/epc_measured_boards.csv).  In our model D1 and D3 have identical L for the
optimal and lateral families (L independent of T) and D1/D2, D3/D4 identical for
vertical (independent of d1), so boards with identical L are compared against the mean
of the corresponding measured points (groups below).

L_spice = L_copper(100 MHz, FastHenry, converged) + L_ESL(Cin bank)
  L_ESL = 0.35 nH / n_caps  (0402 MLCC typical ESL; ASSUMPTION), n = 4 (optimal,
  vertical: 4.1 mm bank) or 2 (lateral: 1.2 mm bank).

Unknown circuit parameters (not published by EPC): effective gate pull-up resistance
(drive strength incl. common-source effects) and HF damping resistance (MLCC ESR at
f_ring, skin/proximity, Coss loss, probe).  Calibrated on the OPTIMAL family only;
LATERAL and VERTICAL are hold-out predictions.
Outputs: results/tables/circuit_calibration_grid.csv, results/tables/circuit_validation.csv
"""
from __future__ import annotations

import itertools
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
import spice_sim as ss  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
TAB = ROOT / "results" / "tables"
MIL = 0.0254

# measured groups: (family, our design key, measured L_epc values to average)
GROUPS = [  # (family, our design (T, d1) mil, centre of EPC's L for the matching measured markers)
    ("optimal", (31, 4), 0.375),
    ("optimal", (31, 12), 0.64),
    ("optimal", (62, 26), 0.875),
    ("lateral", (31, 4), 0.66),
    ("lateral", (31, 12), 1.047),
    ("lateral", (62, 26), 1.374),
    ("vertical", (31, 4), 1.02),
    ("vertical", (62, 4), 1.713),
]
L_WIN = 0.05   # nH: markers within +-L_WIN of the centre belong to the group
N_CAPS = {"optimal": 4, "vertical": 4, "lateral": 2}
ESL_0402 = 0.35


def our_L(sweep: pd.DataFrame, fam, T_mil, d1_mil):
    """L at 100 MHz for (family, T, d1).  WP010 Design 2 (31 mil board, 12 mil inner
    distance) is geometrically impossible with 2 oz copper on all four layers
    (2 x 0.305 + 4 x 0.070 = 0.89 mm > 0.787 mm), so it is absent from our sweep; for the
    optimal and lateral loops L is independent of T (< 0.5 %, epc_sweep_hf.csv), so the
    nearest feasible T is used and flagged."""
    s = sweep[(sweep.family == fam) & ((sweep.d1 / MIL).round() == d1_mil)].copy()
    s["dT"] = ((s["T"] / MIL).round() - T_mil).abs()
    s = s.sort_values("dT")
    if s.dT.iloc[0] > 0 and fam == "vertical":
        raise ValueError("vertical depends on T; no substitution allowed")
    return float(s["L_1e+08_nH"].iloc[0]), float(s["R_1e+08_mohm"].iloc[0])


def measured(meas, fam, Lc):
    """Mean of the pixel-digitised measured markers of one group (data/published/
    epc_measured_boards.csv, long format: family, L_epc_nH, quantity, value)."""
    g = meas[(meas.family == fam) & ((meas.L_epc_nH - Lc).abs() <= L_WIN)]
    os_ = g[g.quantity == "overshoot_pct"].value
    dv = g[g.quantity == "dvdt_V_per_ns"].value
    return (os_.mean() if len(os_) else np.nan, dv.mean() if len(dv) else np.nan,
            g.L_epc_nH.mean(), len(os_), len(dv))


def sim(args):
    L_nH, rup, rdamp = args
    try:
        o = ss.run_bench(ss.Bench(L_loop=L_nH * 1e-9, r_pullup=rup, R_loop=rdamp, t_stop=150e-9))
        return o["overshoot_pct"], o["dvdt_V_per_ns"], o["f_ring"], o["q_ring"]
    except Exception:
        return np.nan, np.nan, np.nan, np.nan


def main():
    sweep = pd.read_csv(TAB / "epc_sweep_hf.csv")
    meas = pd.read_csv(ROOT / "data" / "published" / "epc_measured_boards.csv")
    rows = []
    for fam, (T, d1), Lcen in GROUPS:
        Lc, Rc = our_L(sweep, fam, T, d1)
        os_m, dv_m, Lepc, n_os, n_dv = measured(meas, fam, Lcen)
        rows.append(dict(family=fam, design=f"{T}/{d1}", L_cu_nH=Lc, R_cu_mohm=Rc,
                         L_spice_nH=Lc + ESL_0402 / N_CAPS[fam], L_epc_nH=Lepc, os_meas=os_m, dv_meas=dv_m,
                         n_os=n_os, n_dv=n_dv))
    B = pd.DataFrame(rows)
    # --- calibration grid on the optimal family
    RUP = [2.0, 3.0, 4.0, 5.0, 6.0, 8.0]
    RD = [0.01, 0.03, 0.06, 0.1, 0.15, 0.2, 0.3]
    opt = B[B.family == "optimal"]
    grid = list(itertools.product(RUP, RD))
    jobs = [(L, a, b) for (a, b) in grid for L in opt.L_spice_nH]
    with ProcessPoolExecutor(max_workers=3) as ex:
        res = list(ex.map(sim, jobs))
    G = []
    k = 0
    for (a, b) in grid:
        e_os, e_dv = [], []
        for _, r in opt.iterrows():
            os_p, dv_p, *_ = res[k]; k += 1
            if np.isfinite(r.os_meas):
                e_os.append((os_p - r.os_meas) / 10.0)           # 10 %-points scale
            if np.isfinite(r.dv_meas):
                e_dv.append((dv_p - r.dv_meas) / 1.0)             # 1 V/ns scale
        cost = np.sqrt(np.nanmean(np.square(e_os + e_dv)))
        G.append(dict(r_pullup=a, R_damp=b, cost=cost))
    G = pd.DataFrame(G).sort_values("cost")
    G.to_csv(TAB / "circuit_calibration_grid.csv", index=False)
    best = G.iloc[0]
    print("best calibration:", best.to_dict(), flush=True)
    # --- predictions for all boards: calibrated and nominal (uncalibrated) model
    with ProcessPoolExecutor(max_workers=3) as ex:
        cal = list(ex.map(sim, [(L, best.r_pullup, best.R_damp) for L in B.L_spice_nH]))
        nom = list(ex.map(sim, [(L, 2.7, 2e-3 + 1.5e-3) for L in B.L_spice_nH]))
    B["os_pred"] = [c[0] for c in cal]
    B["dv_pred"] = [c[1] for c in cal]
    B["fring_pred_MHz"] = [c[2] / 1e6 for c in cal]
    B["os_nominal"] = [c[0] for c in nom]
    B["dv_nominal"] = [c[1] for c in nom]
    B["role"] = np.where(B.family == "optimal", "calibration", "hold-out")
    B["r_pullup"] = best.r_pullup
    B["R_damp"] = best.R_damp
    B.to_csv(TAB / "circuit_validation.csv", index=False)
    pd.set_option("display.width", 220)
    print(B.round(2).to_string())
    h = B[B.role == "hold-out"]
    print("hold-out overshoot MAE %.1f %%-pts, dv/dt MAE %.2f V/ns" %
          (np.mean(np.abs(h.os_pred - h.os_meas)), np.nanmean(np.abs(h.dv_pred - h.dv_meas))))


if __name__ == "__main__":
    main()
