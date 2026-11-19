"""Phase 3: independent full-wave (openEMS FDTD) check of the coupled extraction.

Three geometries (EPC Design 1, 31/4 mil): optimal + Kelvin gate return, optimal + common
(ground-plane) gate return, vertical + Kelvin.  Two runs each:
  run A: excite the power port, low-side gate port passive (1 kOhm) -> Z11, Z31
  run B: excite the low-side gate port, power port passive          -> Z33, Z13
The high-side gate port is left open in both, matching the open-circuit definition
of FastHenry's Z matrix.  Copper is PEC sheets at the layer centres (as in Tier 1c),
so the FastHenry reference is its 1 GHz value (current at the conductor surfaces).
Self inductances use the self-resonance fit of oems.fit_L0 over 20-200 MHz; mutual
inductances are Im(Z_transfer)/(2 pi f) averaged over 20-60 MHz.
Plain language: a second, completely different solver recomputes the same numbers.
Mesh: 0.2 mm in-plane and 0.04 mm in z in the copper region (four cells across the
0.17 mm L1-L2 centre gap; one cell underestimated L by ~15 %); mesh lines closer than 50 um are merged onto the sheet
planes (oems min_cell), which keeps the FDTD time step near 0.2 ps instead of the 2.6 fs
forced by micrometre gaps between unrelated copper edges; excitation 0-1 GHz, end at
-30 dB.  Each run takes roughly 15-30 min on the laptop.  Budget: 4 FDTD runs (run B only for the common-return case).  Output: results/phase3/openems_coupled.csv
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
import layouts  # noqa: E402
import oems  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "phase3"
RAW = ROOT / "results" / "raw" / "phase3_oems"
MIL = 0.0254
F_EVAL = np.linspace(20e6, 200e6, 19)
FINE = float(sys.argv[sys.argv.index("--fine") + 1]) if "--fine" in sys.argv else 0.2
MIN_CELL = 0.05


def gate_box(g, p, tag="S"):
    gx, ydrv, zt = g.nodes[f"g{tag}d"]
    xr = g.nodes[f"r{tag}d"][0]
    return ([gx + p.w_gate / 2, ydrv, zt], [xr - p.w_gate / 2, ydrv + p.w_gate, zt], "x")


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    out = OUT / f"openems_coupled_fine{FINE:g}_z0.04.csv"
    rows = pd.read_csv(out).to_dict("records") if out.exists() else []
    done = {(r["case"], r["run"]) for r in rows}
    for fam, gate in (("optimal", "kelvin"), ("optimal", "common"), ("vertical", "kelvin")):
        p = layouts.Params(family=fam, gate=gate, T=31 * MIL, d1=4 * MIL, pitch=0.15)
        g = layouts.build(p)
        gb = gate_box(g, p)
        case = f"{fam}_{gate}"
        runs = [("A_power", g.port_box, gb)] + ([("B_gateLS", gb, g.port_box)] if case == "optimal_common" else [])
        for run, exc, pas in runs:
            if (case, run) in done:
                continue
            r = oems.run_geometry(g, RAW / f"{case}_{run}_f{FINE:g}_z0.04", f_max=1e9, f_eval=F_EVAL, fine=FINE,
                                  margin=20, sheet=True, excite_box=exc, passive_boxes=[pas], passive_R=1e3,
                                  min_cell=MIN_CELL, fine_z=0.04, end_db=-30, max_steps=400_000)
            L0, fr, resid = oems.fit_L0(F_EVAL, r["L"], f_hi=200e6)
            zt = r["Z_transfer"][0]
            m = F_EVAL <= 60e6
            M = float(np.mean(zt.imag[m] / (2 * np.pi * F_EVAL[m])))
            rows.append(dict(case=case, run=run, fine_mm=FINE, L_self_nH=L0 * 1e9, fr_GHz=fr / 1e9,
                             fit_resid=resid, M_nH=M * 1e9, M_100MHz_nH=float(np.interp(1e8, F_EVAL, zt.imag /
                                                                               (2 * np.pi * F_EVAL))) * 1e9,
                             runtime_s=r["runtime_s"], ncells=r["ncells"]))
            pd.DataFrame(rows).to_csv(out, index=False)
            print(case, run, f"L={L0*1e9:.3f} nH M={M*1e9:.4f} nH t={r['runtime_s']:.0f}s cells={r['ncells']}",
                  flush=True)


if __name__ == "__main__":
    main()
