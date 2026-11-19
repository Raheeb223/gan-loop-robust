"""Phase 6 (C4): newer GaN devices - EPC2218 (100 V, chip-scale, 2026 datasheet) and
EPC2302 (100 V, 3 x 5 mm QFN, 2026 datasheet) - next to the validated EPC2015C.

(a) Device models from datasheet tables only (devices.py).  Check: the Coss(V) fit
    reproduces Coss(V_test), Coss,er and Coss,tr; the large-signal Coss ramp of the
    finished ngspice model is compared with the fitted curve.
(b) Coupled 3-port extraction (power + both gate loops) for each device footprint, on
    the JLCPCB 3313 1.6 mm stackup (35/15.2 um copper), optimal and vertical loops,
    Kelvin gate return (and common return for the optimal loop), at 100 MHz.
(c) Switching at the device's intended bus: 48 V / 20 A for the 100 V parts, 12 V / 20 A
    for the 40 V EPC2015C; calibrated drive (6 Ohm pull-up, 0.1 Ohm damping - taken over
    from the EPC2015C calibration, an ASSUMPTION) and a strong 1 Ohm / 40 A stress case.
    Outputs: overshoot, dv/dt, peak off-state gate voltage and false turn-on margin
    (datasheet VGS(TH),min = 0.8 V for all three).
Outputs: results/phase6/coss_fit.csv, coss_check.csv, coupled_L100.csv, switching.csv,
fig_p6_*.pdf/png
"""
from __future__ import annotations

import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
import devices  # noqa: E402
import fhrun  # noqa: E402
import geom  # noqa: E402
import layouts  # noqa: E402
import spice_sim as ss  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "phase6"
RAW = ROOT / "results" / "raw" / "phase6"
PITCH, FW, NH = 0.15, 0.2, 5
CU, CU_IN, T_JLC, D1 = 0.035, 0.0152, 1.6 - 0.035, 0.099
VTH_MIN = 0.8
BUS = {"EPC2015C": 12.0, "EPC2218": 48.0, "EPC2302": 48.0}
ESL = 0.35e-9 / 4


def coss_ramp(dev):
    _, sub = devices.subckt(dev)
    sn = sub.split("\n")[1].split()[1]
    vmax = 0.9 * devices.DEV[dev]["vds"]
    deck = f"""* coss ramp
{sub}
VD d 0 PWL(0 0 400n {vmax})
VG g 0 0
X1 d g 0 {sn}
.control
tran 0.5n 400n
let idd = -i(VD)
wrdata out.txt v(d) idd
.endc
.end
"""
    a = ss.run_deck(deck)["raw"]
    v, c = a[:, 1], a[:, 3] / (vmax / 400e-9) * 1e12
    return v, c


def fh_job(a):
    dev, fam, gate = a
    g0 = devices.DEV[dev]["geo"]
    p = layouts.Params(family=fam, gate=gate, T=T_JLC, d1=D1, cu=CU, cu_in=CU_IN, pitch=PITCH, **g0)
    r = fhrun.run(geom.to_fasthenry(layouts.build(p), (1e8, 1e8), max_fil_w=FW, nh=NH),
                  RAW / f"{dev}_{fam}_{gate}", timeout=3600)
    Z = r["Z"][0]
    w = 2 * np.pi * 1e8
    Lm = 0.5 * (Z.imag + Z.imag.T) / w
    return dict(device=dev, family=fam, gate=gate, L11_nH=Lm[0, 0] * 1e9, L22_nH=Lm[1, 1] * 1e9,
                L33_nH=Lm[2, 2] * 1e9, M12_nH=Lm[0, 1] * 1e9, M13_nH=Lm[0, 2] * 1e9, M23_nH=Lm[1, 2] * 1e9,
                R11_mohm=Z[0, 0].real * 1e3, runtime_s=r["runtime_s"])


def sp_job(kw):
    try:
        o = ss.run_bench(ss.Bench(t_stop=150e-9, **kw))
        return {k: o.get(k, np.nan) for k in ("overshoot_pct", "dvdt_V_per_ns", "f_ring", "vgs_ls_peak")}
    except Exception as ex:
        return dict(err=str(ex)[:200])


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    rows, chk = [], []
    for dev, d in devices.DEV.items():
        C0, Vj, m, res = devices.fit_coss(d["coss"], d["coss_er"], d["coss_tr"], d["v_test"])
        rows.append(dict(device=dev, C0_pF=C0, Vj_V=Vj, m=m, resid_coss_pct=100 * res[0], resid_er_pct=100 * res[1],
                         resid_tr_pct=100 * res[2]))
        if dev != "EPC2015C":
            v, c = coss_ramp(dev)
            for vv in np.array([0.03, 0.05, 0.1, 0.25, 0.4, 0.55, 0.7, 0.85]) * d["vds"]:
                cm = float(np.interp(vv, v, c))
                cf = float(devices.coss_eq(vv, C0, Vj, m))
                chk.append(dict(device=dev, V=vv, coss_model_pF=cm, coss_fit_pF=cf, dev_pct=100 * (cm / cf - 1)))
    pd.DataFrame(rows).to_csv(OUT / "coss_fit.csv", index=False)
    pd.DataFrame(chk).to_csv(OUT / "coss_check.csv", index=False)
    print(pd.DataFrame(rows).round(3).to_string(), "\n", pd.DataFrame(chk).round(2).to_string(), flush=True)

    jobs = [(dev, fam, gate) for dev in devices.DEV for fam, gate in
            (("optimal", "kelvin"), ("optimal", "common"), ("vertical", "kelvin"))]
    if (OUT / "coupled_L100.csv").exists():                 # solver results do not depend on the SPICE model
        L = pd.read_csv(OUT / "coupled_L100.csv")
    else:
        with ProcessPoolExecutor(max_workers=4) as ex:
            L = pd.DataFrame(list(ex.map(fh_job, jobs)))
        L.to_csv(OUT / "coupled_L100.csv", index=False)
    print(L.round(4).to_string(), flush=True)

    cal = pd.read_csv(ROOT / "results" / "tables" / "circuit_validation.csv").iloc[0]
    sj, meta = [], []
    for _, r in L.iterrows():
        L1 = r.L11_nH * 1e-9 + ESL
        k12 = r.M12_nH * 1e-9 / np.sqrt(L1 * r.L22_nH * 1e-9)
        k13 = r.M13_nH * 1e-9 / np.sqrt(L1 * r.L33_nH * 1e-9)
        k23 = r.M23_nH / np.sqrt(r.L22_nH * r.L33_nH)
        for drive, kw in (("calibrated_20A", dict(r_pullup=float(cal.r_pullup), iload=20.0)),
                          ("strong_1ohm_40A", dict(r_pullup=1.0, iload=40.0))):
            sj.append(dict(device=r.device, vin=BUS[r.device], L_loop=L1, R_loop=float(cal.R_damp),
                           L_gate_hs=r.L22_nH * 1e-9, L_gate_ls=r.L33_nH * 1e-9, k_p_hs=k12, k_p_ls=k13, k_hs_ls=k23,
                           **kw))
            meta.append(dict(device=r.device, family=r.family, gate=r.gate, drive=drive, vin=BUS[r.device],
                             L_spice_nH=L1 * 1e9, k13=k13))
    with ProcessPoolExecutor(max_workers=4) as ex:
        res = list(ex.map(sp_job, sj))
    S = pd.DataFrame([dict(m, **x) for m, x in zip(meta, res)])
    S["false_turn_on_margin_V"] = VTH_MIN - S.vgs_ls_peak
    S.to_csv(OUT / "switching.csv", index=False)
    pd.set_option("display.width", 250)
    print(S.round(3).to_string(), flush=True)
    figures(L, S)


def figures(L, S):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axs = plt.subplots(1, 2, figsize=(8, 3.0), constrained_layout=True)
    devs = list(devices.DEV)
    cases = [("optimal", "kelvin"), ("optimal", "common"), ("vertical", "kelvin")]
    for k, (fam, gate) in enumerate(cases):
        l = L[(L.family == fam) & (L.gate == gate)].set_index("device").loc[devs]
        axs[0].bar(np.arange(3) + (k - 1) * 0.27, l.L11_nH, 0.25, label=f"{fam}, {gate}")
        s = S[(S.family == fam) & (S.gate == gate) & (S.drive == "strong_1ohm_40A")].set_index("device").loc[devs]
        axs[1].bar(np.arange(3) + (k - 1) * 0.27, s.vgs_ls_peak, 0.25, label=f"{fam}, {gate}")
    axs[0].set_ylabel("power-loop $L$ at 100 MHz (nH)")
    axs[1].set_ylabel("peak off-state $V_{GS}$, 1 $\\Omega$/40 A (V)")
    axs[1].axhline(VTH_MIN, color="k", ls="--", lw=1)
    for ax in axs:
        ax.set_xticks(range(3), devs)
        ax.legend(fontsize=6)
    for ext in ("pdf", "png"):
        fig.savefig(OUT / f"fig_p6_devices.{ext}", dpi=200)


if __name__ == "__main__":
    main()
