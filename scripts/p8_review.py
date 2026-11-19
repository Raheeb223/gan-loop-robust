"""Robustness checks requested in review (extends Phases 3 and 4).

(a) Mesh convergence of the mutual inductances: M12, M13 (and L11) of the optimal (Kelvin, common)
    and vertical (Kelvin) layouts at four FastHenry discretizations, 100 MHz.
(b) Sensitivity of M13 to the Kelvin pick-off position s_gate (distance of the Kelvin contact
    from the gate pad along the die edge), which is a layout choice the paper fixes at 0.5 mm.
(c) False turn-on severity: the peak off-state gate voltage alone does not say whether the
    low-side channel conducts.  The recorded internal gate waveform is passed through the
    EPC2015C channel law of the device model (I = GM*s*ln(1+exp((V-VT0)/s)), saturation
    branch = upper bound) to obtain the shoot-through current and charge, both for the model
    as fitted and for a worst-case part whose threshold is shifted to the datasheet minimum
    (9 mA at 0.8 V).
(d) Prepreg tolerance assumption: Sobol indices recomputed on the surrogate with the prepreg
    tolerance narrowed to +-5 % and +-3 %, and with Gaussian inputs (limit = 3 sigma) instead of
    uniform ones.
Outputs: results/review/{mutual_convergence,kelvin_position,shoot_through,sobol_tolerance}.csv
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
OUT = ROOT / "results" / "review"
RAW = ROOT / "results" / "raw" / "review"
MIL = 0.0254
MESH = [(0.30, 0.40, 3), (0.15, 0.20, 5), (0.10, 0.13, 7), (0.075, 0.10, 7)]   # pitch, max_fil_w, nh
CASES = [("optimal", "kelvin"), ("optimal", "common"), ("vertical", "kelvin")]


def fh_job(a):
    fam, gate, pitch, fw, nh, s_gate, tag = a
    p = layouts.Params(family=fam, gate=gate, T=31 * MIL, d1=4 * MIL, pitch=pitch, s_gate=s_gate)
    try:
        r = fhrun.run(geom.to_fasthenry(layouts.build(p), (1e8, 1e8), max_fil_w=fw, nh=nh), RAW / tag, timeout=7200)
        Z = r["Z"][0]
        w = 2 * np.pi * 1e8
        L = Z.imag / w * 1e9
        return dict(family=fam, gate=gate, pitch=pitch, max_fil_w=fw, nh=nh, s_gate=s_gate, L11_nH=L[0, 0],
                    L33_nH=L[2, 2], M12_nH=0.5 * (L[0, 1] + L[1, 0]), M13_nH=0.5 * (L[0, 2] + L[2, 0]),
                    runtime_s=r["runtime_s"], err="")
    except Exception as ex:
        return dict(family=fam, gate=gate, pitch=pitch, max_fil_w=fw, nh=nh, s_gate=s_gate, err=str(ex)[:200])


# ---------------------------------------------------------------- (c) shoot-through
GM, S_, VT0 = 150.0, 0.15, 2.25            # EPC2015C channel law of epc_model.py


def i_ch(v, shift=0.0):
    x = (v - (VT0 + shift)) / S_
    return GM * S_ * np.logaddexp(0, x)


def vth_at_9mA(shift=0.0):
    from scipy.optimize import brentq
    return brentq(lambda v: i_ch(v, shift) - 9e-3, -2, 5)


def spice_job(kw):
    import spice_sim as ss
    o = ss.run_bench(ss.Bench(t_stop=150e-9, **kw))
    t, v, vsw = o["t"], o["vgs_ls_int"], o["vsw"]
    m = t >= 5e-9 + o["dead"]
    shift_min = 0.8 - vth_at_9mA()
    res = dict(vgs_ls_peak=float(v[m].max()))
    for name, sh in (("model", 0.0), ("vth_min", shift_min)):
        i = i_ch(v[m], sh)
        res[f"ipk_{name}_A"] = float(i.max())
        res[f"q_{name}_nC"] = float(np.trapezoid(i, t[m]) * 1e9)
        res[f"e_{name}_uJ"] = float(np.trapezoid(i * np.clip(vsw[m], 0, None), t[m]) * 1e6)
    return res


def shoot_through():
    import p3_spice as P3
    Lt = pd.read_csv(ROOT / "results" / "phase3" / "coupled_L100.csv")
    cal = pd.read_csv(ROOT / "results" / "tables" / "circuit_validation.csv").iloc[0]
    rup, rd = float(cal.r_pullup), float(cal.R_damp)
    jobs, meta = [], []
    for _, r in Lt[Lt.case.str.endswith("_sweep")].iterrows():
        L1 = r.L11_nH * 1e-9 + P3.ESL
        k12 = r.L12_nH * 1e-9 / np.sqrt(L1 * r.L22_nH * 1e-9)
        k13 = r.L13_nH * 1e-9 / np.sqrt(L1 * r.L33_nH * 1e-9)
        k23 = r.L23_nH / np.sqrt(r.L22_nH * r.L33_nH)
        base = dict(L_loop=L1, L_gate_hs=r.L22_nH * 1e-9, L_gate_ls=r.L33_nH * 1e-9, k_p_hs=k12, k_p_ls=k13,
                    k_hs_ls=k23, R_loop=rd)
        for drive, extra in (("calibrated_20A", dict(r_pullup=rup, iload=20.0)),
                             ("strong_1ohm_40A", dict(r_pullup=1.0, iload=40.0))):
            jobs.append(dict(base, **extra))
            meta.append(dict(case=r.case.replace("_sweep", ""), drive=drive))
    with ProcessPoolExecutor(max_workers=2) as ex:
        res = list(ex.map(spice_job, jobs))
    S = pd.DataFrame([dict(m, **x) for m, x in zip(meta, res)])
    S.attrs["vth_model"] = vth_at_9mA()
    return S


# ---------------------------------------------------------------- (d) Sobol vs tolerance assumption
def sobol_tolerance():
    from SALib.analyze import sobol
    from SALib.sample import sobol as sobol_sample
    import p5_opt
    import p4_uq
    from p4_campaign import HI, LO, NOMINAL, TOL
    sur = p5_opt.load_sur()
    rows = []
    for dist in ("unif", "norm"):
        for tol_d1 in (0.10, 0.05, 0.03):
            lo, hi = LO[5:].copy(), HI[5:].copy()
            lo[0], hi[0] = 1 - tol_d1, 1 + tol_d1
            if dist == "unif":
                bounds = [[a, b] for a, b in zip(lo, hi)]
                prob = dict(num_vars=len(TOL), names=TOL, bounds=bounds)
            else:                                   # Gaussian, +-limit = 3 sigma, truncated to the limits
                prob = dict(num_vars=len(TOL), names=TOL,
                            bounds=[[(a + b) / 2, (b - a) / 6] for a, b in zip(lo, hi)], dists=["norm"] * len(TOL))
            for d1n, name in ((0.099, "3313"), (0.210, "7628")):
                Xs = sobol_sample.sample(prob, 4096, calc_second_order=False, seed=3)
                Xs = np.clip(Xs, lo, hi)
                m, _ = sur.predict(p4_uq.X_at(dict(NOMINAL, d1=d1n), Xs))
                L = np.exp(m)
                Si = sobol.analyze(prob, L, calc_second_order=False, seed=3)
                lo95, hi95 = np.percentile(L, [2.5, 97.5])
                for k, v in enumerate(TOL):
                    rows.append(dict(dist=dist, tol_d1=tol_d1, stackup=name, variable=v, ST=Si["ST"][k],
                                     S1=Si["S1"][k], half_range95_pct=50 * (hi95 - lo95) / np.median(L)))
    return pd.DataFrame(rows)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    pd.set_option("display.width", 250)
    what = sys.argv[1:] or ["sobol", "shoot", "mesh", "kelvin"]
    if "sobol" in what:
        S = sobol_tolerance()
        S.to_csv(OUT / "sobol_tolerance.csv", index=False)
        print(S.pivot_table(index=["dist", "tol_d1", "stackup"], columns="variable", values="ST").round(3), flush=True)
        print(S.groupby(["dist", "tol_d1", "stackup"]).half_range95_pct.first().round(2), flush=True)
    if "shoot" in what:
        S = shoot_through()
        S.to_csv(OUT / "shoot_through.csv", index=False)
        print("model threshold at 9 mA:", round(S.attrs["vth_model"], 3), "V", flush=True)
        print(S.round(4).to_string(), flush=True)
    jobs = []
    if "mesh" in what:
        for fam, gate in CASES:
            for pitch, fw, nh in MESH:
                jobs.append((fam, gate, pitch, fw, nh, 0.5, f"mesh_{fam}_{gate}_{pitch}_{nh}"))
    if "kelvin" in what:
        for gate in ("kelvin", "common"):
            for sg in (0.25, 1.0, 1.5):
                jobs.append(("optimal", gate, 0.15, 0.2, 5, sg, f"sgate_optimal_{gate}_{sg}"))
    if jobs:
        with ProcessPoolExecutor(max_workers=2) as ex:
            R = pd.DataFrame(list(ex.map(fh_job, jobs)))
        R.to_csv(OUT / "mutual_convergence.csv", index=False)
        print(R.round(4).to_string(), flush=True)


if __name__ == "__main__":
    main()
