"""Phase 3: from the extracted impedance matrix to circuit behaviour.

(a) Passive RL ladder ("Foster") fit of the power-loop impedance Z11(f), 1 MHz - 1 GHz:
        Z(s) = R0 + s Linf + sum_k  s r_k / (s + p_k)
    Each term s r/(s+p) is a resistor r in parallel with an inductor r/p, so the fitted
    model is a chain of ordinary R and L parts.  Poles p_k are fixed on a log grid and
    R0, Linf, r_k are found by non-negative least squares: every element is >= 0, so the
    model is passive by construction (a pole-fixed variant of vector fitting).
    Plain language: the copper's resistance and inductance change with frequency (skin
    and proximity effect); the ladder reproduces that in a normal SPICE netlist.
(b) K-coupled ngspice bench: power loop, high-side and low-side gate loops with the
    extracted 100 MHz self and mutual inductances (K = M / sqrt(L1 L2)).  Reports
    overshoot, dv/dt and the peak low-side gate voltage during the high-side turn-on;
    false turn-on margin = VGS(TH),min (datasheet, 0.8 V) - peak.
    Two drive cases: the calibrated EPC drive (6 Ohm pull-up, Section VII) and a strong
    1 Ohm driver at 40 A (stress case, fastest edges).
(c) dv/dt hypothesis test: dv/dt with and without the extracted gate coupling.
Outputs: results/phase3/foster_fit.csv, spice_coupled.csv, fig_p3_*.pdf/png
"""
from __future__ import annotations

import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import nnls

sys.path.insert(0, str(Path(__file__).parent))
import spice_sim as ss  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "phase3"
ESL = 0.35e-9 / 4          # 0402 bank ESL, as in the circuit validation
VTH_MIN = 0.8              # EPC2015C datasheet VGS(TH) min
POLES_HZ = np.logspace(6, 9.5, 8)


def foster_fit(f, Z):
    s = 2j * np.pi * f
    p = 2 * np.pi * POLES_HZ
    cols = [np.ones_like(s), s] + [s / (s + pk) for pk in p]
    A = np.array(cols).T
    w = 1 / np.abs(Z)
    Ar = np.vstack([(A.real.T * w).T, (A.imag.T * w).T])
    br = np.concatenate([Z.real * w, Z.imag * w])
    x, _ = nnls(Ar, br)
    Zf = A @ x
    return x, Zf


def ladder_netlist(x, R_extra=0.0):
    """Lloop (Linf) -> R0 -> parallel RL sections, from node vc to dh2."""
    R0, Linf, rk = x[0], x[1], x[2:]
    lines = [f"Lloop vc fl0 {Linf:.6e}", f"Rfl0 fl0 fl1 {R0 + R_extra:.6e}"]
    n = 1
    for r, pk in zip(rk, 2 * np.pi * POLES_HZ):
        if r <= 0:
            continue
        lines += [f"Rfs{n} fl{n} fl{n+1} {r:.6e}", f"Lfs{n} fl{n} fl{n+1} {r / pk:.6e}"]
        n += 1
    lines.append(f"Rfend fl{n} dh2 1e-9")
    return "\n".join(lines)


def sim(kw):
    try:
        o = ss.run_bench(ss.Bench(t_stop=150e-9, **kw))
        return {k: o.get(k, np.nan) for k in ("overshoot_pct", "dvdt_V_per_ns", "f_ring", "vgs_ls_peak",
                                                "vgs_ls_min")}
    except Exception as ex:
        return dict(err=str(ex)[:200])


def main():
    Z = pd.read_csv(OUT / "coupled_Z.csv")
    cal = pd.read_csv(ROOT / "results" / "tables" / "circuit_validation.csv").iloc[0]
    rup, rd = float(cal.r_pullup), float(cal.R_damp)
    sweeps = sorted(c for c in Z.case.unique() if c.endswith("_sweep"))
    # ---------------- (a) Foster fit of Z11
    frows, ladders = [], {}
    for c in sweeps:
        z = Z[(Z.case == c) & (Z.i == 0) & (Z.j == 0)].sort_values("f_Hz")
        f = z.f_Hz.values
        Zc = z.R_mohm.values * 1e-3 + 1j * 2 * np.pi * f * z.L_nH.values * 1e-9
        x, Zf = foster_fit(f, Zc)
        ladders[c] = x
        err = np.abs(Zf - Zc) / np.abs(Zc)
        frows.append(dict(case=c, max_rel_err_pct=100 * err.max(), R0_mohm=x[0] * 1e3, Linf_nH=x[1] * 1e9,
                          n_sections=int((x[2:] > 0).sum()), all_nonneg=bool((x >= 0).all()),
                          L_1MHz_nH=z.L_nH.values[0], L_1GHz_nH=z.L_nH.values[-1]))
    F = pd.DataFrame(frows)
    F.to_csv(OUT / "foster_fit.csv", index=False)
    print(F.round(4).to_string(), flush=True)

    # ---------------- (b) K-coupled benches
    z = Z[np.isclose(Z.f_Hz, 1e8)]
    jobs, meta = [], []
    for c in sorted(z.case.unique()):
        Lm = np.zeros((3, 3))
        for _, r in z[z.case == c].iterrows():
            Lm[int(r.i), int(r.j)] = r.L_nH * 1e-9
        Lm = 0.5 * (Lm + Lm.T)
        L1 = Lm[0, 0] + ESL
        k12 = Lm[0, 1] / np.sqrt(L1 * Lm[1, 1])
        k13 = Lm[0, 2] / np.sqrt(L1 * Lm[2, 2])
        k23 = Lm[1, 2] / np.sqrt(Lm[1, 1] * Lm[2, 2])
        R11 = float(Z[(Z.case == c) & np.isclose(Z.f_Hz, 1e8) & (Z.i == 0) & (Z.j == 0)].R_mohm.iloc[0]) * 1e-3
        base = dict(L_loop=L1, L_gate_hs=Lm[1, 1], L_gate_ls=Lm[2, 2])
        for drive, extra in (("calibrated_20A", dict(r_pullup=rup, R_loop=rd, iload=20.0)),
                             ("strong_1ohm_40A", dict(r_pullup=1.0, R_loop=rd, iload=40.0))):
            for coup, kk in (("coupled", dict(k_p_hs=k12, k_p_ls=k13, k_hs_ls=k23)), ("uncoupled", {})):
                jobs.append(dict(base, **extra, **kk))
                meta.append(dict(case=c, drive=drive, coupling=coup, L11_nH=Lm[0, 0] * 1e9, L22_nH=Lm[1, 1] * 1e9,
                                 L33_nH=Lm[2, 2] * 1e9, M12_nH=Lm[0, 1] * 1e9, M13_nH=Lm[0, 2] * 1e9,
                                 k12=k12, k13=k13, R11_mohm=R11 * 1e3))
    # ladder vs fixed-L comparison (calibrated drive, no gate coupling)
    for c, x in ladders.items():
        net = f"Lesl vc fl00 {ESL:.4e}\n" + ladder_netlist(x, R_extra=rd).replace("Lloop vc fl0 ", "Lloop fl00 fl0 ")
        jobs.append(dict(r_pullup=rup, R_loop=rd, iload=20.0, ploop_netlist=net, L_gate_hs=1e-9, L_gate_ls=1e-9))
        meta.append(dict(case=c, drive="calibrated_20A", coupling="foster_ladder"))
        L100 = float(z[(z.case == c) & (z.i == 0) & (z.j == 0)].L_nH.iloc[0]) * 1e-9
        jobs.append(dict(r_pullup=rup, R_loop=rd, iload=20.0, L_loop=L100 + ESL))
        meta.append(dict(case=c, drive="calibrated_20A", coupling="fixed_L100MHz"))
    with ProcessPoolExecutor(max_workers=4) as ex:
        res = list(ex.map(sim, jobs))
    S = pd.DataFrame([dict(m, **r) for m, r in zip(meta, res)])
    S["false_turn_on_margin_V"] = VTH_MIN - S.vgs_ls_peak
    S.to_csv(OUT / "spice_coupled.csv", index=False)
    pd.set_option("display.width", 250)
    print(S.drop(columns=[c for c in S.columns if c.startswith("err")], errors="ignore").round(3).to_string())
    figures(S)


def figures(S):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    s = S[(S.coupling == "coupled") & S.case.str.endswith("_sweep")].copy()
    s["label"] = s.case.str.replace("_sweep", "").str.replace("_", "\n")
    fig, axs = plt.subplots(1, 2, figsize=(8, 3.0), constrained_layout=True)
    for ax, drive in zip(axs, ("calibrated_20A", "strong_1ohm_40A")):
        d = s[s.drive == drive]
        ax.bar(range(len(d)), d.vgs_ls_peak, color=["C3" if "common" in c else "C0" for c in d.case])
        ax.axhline(VTH_MIN, color="k", ls="--", lw=1)
        ax.text(len(d) - 0.5, VTH_MIN, r" $V_\mathrm{GS(th),min}$", va="bottom", ha="right", fontsize=7)
        ax.set_xticks(range(len(d)), d.label, fontsize=6)
        ax.set_ylabel("peak low-side $V_{GS}$ (V)")
        ax.set_title(drive.replace("_", ", "), fontsize=9)
    for ext in ("pdf", "png"):
        fig.savefig(OUT / f"fig_p3_false_turn_on.{ext}", dpi=200)


if __name__ == "__main__":
    main()
