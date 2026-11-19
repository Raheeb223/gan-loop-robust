"""How reliable is the ringing method for 'measuring' PCB power-loop inductance?

Part A (method bias, simulation with KNOWN L):
  Simulate the hard-commutation bench (EPC2015C datasheet model, 12 V, 20 A) for
  L_true in 0.2..3 nH; read f_ring from the switch-node waveform exactly as one
  would from a scope; estimate L = 1 / ((2 pi f)^2 C) with the C an engineer
  would pick from the datasheet: Coss(VIN) typ, Coss(tr), Coss(er).
  -> bias of each choice vs L_true.

Part B (input uncertainty, Monte Carlo):
  (i)  datasheet Coss spread: EPC2015C Coss(20 V) = 710 pF typ, 1070 pF max
       -> part-to-part C uniform in [typ, max] (conservative; no min is given).
  (ii) Letellier et al. (TPEL 2018) decomposition of their GS61008P prototype:
       L_stray = 1.45 nH from ringing with Coss(50 V) = 280 pF, minus 2 x 0.1 nH
       package and 0.87 nH capacitor inductance -> L_pcb = 0.38 nH.
       We propagate +-30 % (uniform) on each subtracted quantity and on Coss and
       report the resulting distribution of L_pcb.
Outputs: results/tables/ring_bias.csv, results/tables/ring_mc.csv
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
import epc_model as em  # noqa: E402
import spice_sim as ss  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
TAB = ROOT / "results" / "tables"


def part_a(Ls=(0.2, 0.3, 0.4, 0.6, 0.8, 1.0, 1.3, 1.6, 2.0, 2.5, 3.0), **bench_kw):
    vin = 12.0
    c_vin = float(np.interp(vin, em.V_PTS, em.COSS_PF)) * 1e-12     # Coss(12 V) typ from Fig. 5a
    c_tr, c_er = 940e-12, 870e-12                                      # datasheet table (0-20 V)
    rows = []
    for L in Ls:
        o = ss.run_bench(ss.Bench(L_loop=L * 1e-9, t_stop=150e-9, **bench_kw))
        f = o["f_ring"]
        row = dict(L_true_nH=L, f_ring_MHz=f / 1e6, overshoot_pct=o["overshoot_pct"], q_ring=o["q_ring"])
        for name, c in (("Coss_VIN", c_vin), ("Coss_tr", c_tr), ("Coss_er", c_er)):
            Le = 1 / ((2 * np.pi * f) ** 2 * c) * 1e9
            row[f"L_est_{name}_nH"] = Le
            row[f"err_{name}_pct"] = 100 * (Le - L) / L
        # effective resonant capacitance actually seen by the ring
        row["C_eff_pF"] = 1 / ((2 * np.pi * f) ** 2 * L * 1e-9) * 1e12
        rows.append(row)
        print(row, flush=True)
    return pd.DataFrame(rows)


def part_b(n=200_000, seed=1):
    rng = np.random.default_rng(seed)
    out = {}
    # (i) Coss part-to-part spread: true C uniform in [typ, max]; estimate uses typ
    c_typ, c_max = 710e-12, 1070e-12
    c_true = rng.uniform(c_typ, c_max, n)
    out["coss_spread_Lest_over_Ltrue"] = c_true / c_typ     # L_est/L_true = C_true/C_typ
    # (ii) Letellier decomposition with +-30 % uniform on each input
    f_ring = 1 / (2 * np.pi * np.sqrt(1.45e-9 * 280e-12))  # back-computed from their numbers
    u = lambda: rng.uniform(0.7, 1.3, n)
    L_stray = 1 / ((2 * np.pi * f_ring) ** 2 * (280e-12 * u()))
    L_pcb = L_stray - 2 * 0.1e-9 * u() - 0.87e-9 * u()
    out["letellier_L_pcb_nH"] = L_pcb * 1e9
    q = lambda a: np.percentile(a, [2.5, 50, 97.5])
    rows = [dict(quantity="L_est/L_true from Coss typ when part is typ..max",
                 p2_5=q(out["coss_spread_Lest_over_Ltrue"])[0], p50=q(out["coss_spread_Lest_over_Ltrue"])[1],
                 p97_5=q(out["coss_spread_Lest_over_Ltrue"])[2]),
            dict(quantity="Letellier2018 L_pcb (nH), +-30% on Coss, L_pkg, L_cap",
                 p2_5=q(out["letellier_L_pcb_nH"])[0], p50=q(out["letellier_L_pcb_nH"])[1],
                 p97_5=q(out["letellier_L_pcb_nH"])[2],
                 frac_negative=float(np.mean(out["letellier_L_pcb_nH"] < 0)),
                 f_ring_MHz=f_ring / 1e6)]
    return pd.DataFrame(rows), out


if __name__ == "__main__":
    a = part_a()
    a.to_csv(TAB / "ring_bias.csv", index=False)
    b, raw = part_b()
    b.to_csv(TAB / "ring_mc.csv", index=False)
    np.save(TAB / "ring_mc_letellier_samples.npy", raw["letellier_L_pcb_nH"][:20000])
    print(b.to_string())
