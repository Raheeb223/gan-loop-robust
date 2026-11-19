"""Phase 2: robustness of the circuit-level validation.

(a) Calibration error surface and leave-one-family-out (LOFO) cross-validation.
    Plain language: the original study tuned two unknown circuit numbers (gate pull-up
    resistance and a damping resistance) on the "optimal" boards and predicted the
    others.  Here every family takes a turn as the tuning set, so we see whether the
    good hold-out result depends on which family was used for tuning.
    Every board is simulated on the full 6 x 7 grid once (336 ngspice runs); the
    calibrations are then pure post-processing.

(b) dv/dt definition test.  EPC did not publish how dv/dt was measured.  For the
    calibrated model we compute four common definitions from the same simulated
    waveform (10-90 %, 20-80 %, 0-100 % of VIN, maximum slope) and compare their
    spread with the measured trend.

Outputs (results/phase2/): crossval_grid_all.csv, crossval_lofo.csv,
    dvdt_definitions.csv, fig_p2_error_surface.(pdf|png), fig_p2_dvdt_defs.(pdf|png)
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
OUT = ROOT / "results" / "phase2"
RUP = [2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 12.0, 16.0]
RD = [0.01, 0.03, 0.06, 0.1, 0.15, 0.2, 0.3]


def sim(args):
    L_nH, rup, rdamp = args
    for attempt in range(2):
        try:
            o = ss.run_bench(ss.Bench(L_loop=L_nH * 1e-9, r_pullup=rup, R_loop=rdamp, t_stop=150e-9))
            return o["overshoot_pct"], o["dvdt_V_per_ns"]
        except Exception:
            continue
    return np.nan, np.nan


def cost(df, os_col="os_p", dv_col="dv_p"):
    e = list((df[os_col] - df.os_meas).dropna() / 10.0) + list((df[dv_col] - df.dv_meas).dropna() / 1.0)
    return float(np.sqrt(np.mean(np.square(e))))


def dvdt_defs(t, v, vin, t_on):
    i0 = int(np.searchsorted(t, t_on))
    out = {}
    for name, lo, hi in (("10-90", 0.1, 0.9), ("20-80", 0.2, 0.8), ("0-100", 0.02, 1.0)):
        a, _ = ss._cross(t, v, lo * vin, i0)
        b, _ = ss._cross(t, v, hi * vin, i0)
        out[name] = (hi - lo) * vin / (b - a) * 1e-9 if np.isfinite(a) and np.isfinite(b) else np.nan
    tt, idx = np.unique(t[i0:], return_index=True)
    vv = v[i0:][idx]
    out["max_slope"] = float(np.max(np.diff(vv) / np.diff(tt))) * 1e-9
    return out


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    B = pd.read_csv(ROOT / "results" / "tables" / "circuit_validation.csv")
    grid = list(itertools.product(RUP, RD))
    gf = OUT / "crossval_grid_all.csv"
    old = pd.read_csv(gf) if gf.exists() else pd.DataFrame(columns=["r_pullup", "R_damp"])
    have = {(round(a, 6), round(b, 6)) for a, b in zip(old.r_pullup, old.R_damp)}
    todo = [(a, b) for (a, b) in grid if (round(a, 6), round(b, 6)) not in have]
    jobs = [(L, a, b) for (a, b) in todo for L in B.L_spice_nH]
    with ProcessPoolExecutor(max_workers=3) as ex:
        res = list(ex.map(sim, jobs, chunksize=4))
    rows, k = old.to_dict("records"), 0
    for (a, b) in todo:
        for _, r in B.iterrows():
            rows.append(dict(r_pullup=a, R_damp=b, family=r.family, design=r.design, L_spice_nH=r.L_spice_nH,
                             os_meas=r.os_meas, dv_meas=r.dv_meas, os_p=res[k][0], dv_p=res[k][1]))
            k += 1
    G = pd.DataFrame(rows)
    G.to_csv(OUT / "crossval_grid_all.csv", index=False)

    lofo = []
    for fam in ("optimal", "lateral", "vertical"):
        for target in ("overshoot+dvdt", "overshoot"):
            best = None
            for (a, b), g in G.groupby(["r_pullup", "R_damp"]):
                cal = g[g.family == fam]
                c = cost(cal) if target == "overshoot+dvdt" else cost(cal.assign(dv_meas=np.nan))
                if best is None or c < best[0]:
                    best = (c, a, b)
            c, a, b = best
            g = G[(G.r_pullup == a) & (G.R_damp == b)]
            h = g[g.family != fam]
            lofo.append(dict(calibration_family=fam, objective=target, r_pullup=a, R_damp=b, cal_cost=c,
                             interior=(a not in (RUP[0], RUP[-1])) and (b not in (RD[0], RD[-1])),
                             holdout_os_MAE=float(np.mean(np.abs(h.os_p - h.os_meas))),
                             holdout_dv_MAE=float(np.nanmean(np.abs(h.dv_p - h.dv_meas))),
                             n_holdout=len(h)))
    L = pd.DataFrame(lofo)
    L.to_csv(OUT / "crossval_lofo.csv", index=False)
    print(L.round(3).to_string())

    # ---- (b) dv/dt definitions with the original calibration (optimal family)
    cv = pd.read_csv(ROOT / "results" / "tables" / "circuit_validation.csv")
    rup, rd = float(cv.r_pullup.iloc[0]), float(cv.R_damp.iloc[0])
    drows = []
    for _, r in cv.iterrows():
        o = ss.run_bench(ss.Bench(L_loop=r.L_spice_nH * 1e-9, r_pullup=rup, R_loop=rd, t_stop=150e-9))
        d = dvdt_defs(o["t"], o["vsw"], 12.0, 5e-9 + 10e-9)
        drows.append(dict(family=r.family, design=r.design, L_spice_nH=r.L_spice_nH, dv_meas=r.dv_meas, **d))
    D = pd.DataFrame(drows)
    D.to_csv(OUT / "dvdt_definitions.csv", index=False)
    print(D.round(2).to_string())
    figures(G, L, D)


def figures(G, L, D):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axs = plt.subplots(1, 3, figsize=(10, 3.1), constrained_layout=True)
    for ax, fam in zip(axs, ("optimal", "lateral", "vertical")):
        Z = np.zeros((len(RUP), len(RD)))
        for i, a in enumerate(RUP):
            for j, b in enumerate(RD):
                g = G[(G.r_pullup == a) & (G.R_damp == b) & (G.family == fam)]
                Z[i, j] = cost(g)
        im = ax.imshow(Z, origin="lower", aspect="auto", cmap="viridis",
                       extent=[-0.5, len(RD) - 0.5, -0.5, len(RUP) - 0.5])
        r = L[(L.calibration_family == fam) & (L.objective == "overshoot+dvdt")].iloc[0]
        ax.plot(RD.index(r.R_damp), RUP.index(r.r_pullup), "w*", ms=12)
        ax.set_xticks(range(len(RD)), [f"{x:g}" for x in RD], fontsize=7)
        ax.set_yticks(range(len(RUP)), [f"{x:g}" for x in RUP], fontsize=7)
        ax.set_xlabel(r"$R_\mathrm{damp}$ ($\Omega$)")
        ax.set_title(f"calibrated on {fam}", fontsize=9)
        fig.colorbar(im, ax=ax, shrink=0.85, label="cost")
    axs[0].set_ylabel(r"$R_\mathrm{pull\text{-}up}$ ($\Omega$)")
    for ext in ("pdf", "png"):
        fig.savefig(OUT / f"fig_p2_error_surface.{ext}", dpi=200)
    fig, ax = plt.subplots(figsize=(4.2, 3.1), constrained_layout=True)
    Ds = D.sort_values("L_spice_nH")
    for c, m in (("10-90", "o"), ("20-80", "s"), ("0-100", "^"), ("max_slope", "d")):
        ax.plot(Ds.L_spice_nH, Ds[c], m + "-", label=c, ms=4)
    ax.plot(Ds.L_spice_nH, Ds.dv_meas, "k*", ms=9, label="measured [EPC]")
    ax.set_xlabel("loop inductance in SPICE (nH)")
    ax.set_ylabel("dv/dt (V/ns)")
    ax.legend(fontsize=7)
    for ext in ("pdf", "png"):
        fig.savefig(OUT / f"fig_p2_dvdt_defs.{ext}", dpi=200)


if __name__ == "__main__":
    main()
