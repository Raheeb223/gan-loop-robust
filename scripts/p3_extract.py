"""Phase 3 (C1): coupled multi-port extraction of power loop, gate loops and their coupling.

Plain language: besides the power loop, each transistor has a small gate-drive loop.
If the gate loop shares copper with the power loop (common-source inductance) or is
magnetically coupled to it, the large power-loop di/dt induces a voltage in the gate
loop that can slow switching or falsely turn a transistor on.  FastHenry returns the
full port impedance matrix Z(f); with ports [power, gate_HS, gate_LS] its imaginary
part / (2 pi f) is the 3x3 loop inductance matrix: diagonal = self inductances,
off-diagonal = mutual inductances (M), which include the common-source part.

(a) Verification of mutual inductance: two parallel bars, partial mutual inductance
    from FastHenry vs Grover's closed form for parallel filaments.
(b) Coupled extraction on the EPC Design-1 stackup (31/4 mil, 2 oz) for the layout
    families with Kelvin and with common (ground-plane) gate return, swept 1 MHz - 1 GHz
    (2 points/decade) for the RL fit of p3_spice.py.
(c) Gate-trace length sweep (100 MHz) for the common-return low-side gate.
Outputs: results/phase3/grover_mutual.csv, coupled_Z.csv (long format), coupled_L100.csv
"""
from __future__ import annotations

import json
import sys
from concurrent.futures import ProcessPoolExecutor
from dataclasses import asdict
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
import fhrun  # noqa: E402
import geom  # noqa: E402
import layouts  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "phase3"
RAW = ROOT / "results" / "raw" / "phase3"
PITCH, FW, NH = 0.15, 0.2, 5
MIL = 0.0254
MU0_2PI = 2e-7  # H/m


def grover_M(l, d):
    """Partial mutual inductance of two parallel, aligned filaments (Grover 1946), SI."""
    return MU0_2PI * l * (np.log(l / d + np.sqrt(1 + (l / d) ** 2)) - np.sqrt(1 + (d / l) ** 2) + d / l)


def grover_job(d):
    l, w = 10.0, 0.035
    g = geom.Geometry(f"pair_d{d}")
    g.node("a1", 0, 0, 0); g.node("a2", l, 0, 0)
    g.node("b1", 0, d, 0); g.node("b2", l, d, 0)
    g.bar("a1", "a2", w, w); g.bar("b1", "b2", w, w)
    g.ports += [("a1", "a2"), ("b1", "b2")]
    r = fhrun.run(geom.to_fasthenry(g, (1e3, 1e3), max_fil_w=w, nh=1), RAW / f"grover_d{d}")
    Z = r["Z"][0]
    w_ = 2 * np.pi * 1e3
    M_fh = Z[0, 1].imag / w_
    M_gr = grover_M(l * 1e-3, d * 1e-3)
    return dict(l_mm=l, d_mm=d, M_fasthenry_nH=M_fh * 1e9, M_grover_nH=M_gr * 1e9,
                err_pct=100 * (M_fh - M_gr) / M_gr, L11_nH=Z[0, 0].imag / w_ * 1e9)


def coupled_job(a):
    p, freqs, ndec, tag = a
    g = layouts.build(p)
    r = fhrun.run(geom.to_fasthenry(g, freqs, ndec=ndec, max_fil_w=FW, nh=NH), RAW / tag)
    rows = []
    for f, Z in zip(r["freqs"], r["Z"]):
        w = 2 * np.pi * f
        for i in range(Z.shape[0]):
            for j in range(Z.shape[1]):
                rows.append(dict(case=tag, family=p.family, gate=p.gate, l_gate=p.l_gate, f_Hz=f, i=i, j=j,
                                 R_mohm=Z[i, j].real * 1e3, L_nH=Z[i, j].imag / w * 1e9))
    return rows, dict(case=tag, runtime_s=r["runtime_s"], params=json.dumps(asdict(p)))


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    with ProcessPoolExecutor(max_workers=4) as ex:
        G = pd.DataFrame(list(ex.map(grover_job, [0.25, 0.5, 1.0, 2.0, 4.0])))
        G.to_csv(OUT / "grover_mutual.csv", index=False)
        print(G.round(4).to_string(), flush=True)
        base = dict(T=31 * MIL, d1=4 * MIL, pitch=PITCH)
        jobs = []
        for fam, gate in (("optimal", "kelvin"), ("optimal", "common"), ("vertical", "kelvin"),
                          ("lateral", "kelvin"), ("lateral", "common"), ("lateral_stitched", "kelvin"),
                          ("lateral_stitched", "common")):
            p = layouts.Params(family=fam, gate=gate, **base)
            jobs.append((p, (1e6, 1e9), 2, f"{fam}_{gate}_sweep"))
        for lg in (1.0, 6.0):
            for gate in ("kelvin", "common"):
                p = layouts.Params(family="optimal", gate=gate, l_gate=lg, **base)
                jobs.append((p, (1e8, 1e8), 1, f"optimal_{gate}_lg{lg:g}"))
        res = list(ex.map(coupled_job, jobs))
    Z = pd.DataFrame([r for rows, _ in res for r in rows])
    Z.to_csv(OUT / "coupled_Z.csv", index=False)
    pd.DataFrame([m for _, m in res]).to_csv(OUT / "coupled_runs.csv", index=False)
    z = Z[np.isclose(Z.f_Hz, 1e8)]
    piv = z.pivot_table(index=["case", "family", "gate", "l_gate"], columns=["i", "j"], values="L_nH")
    piv.columns = [f"L{i+1}{j+1}_nH" for i, j in piv.columns]
    piv = piv.reset_index()
    # coupling coefficient power -> LS gate
    piv["k13"] = piv.L13_nH / np.sqrt(piv.L11_nH * piv.L33_nH)
    piv["k12"] = piv.L12_nH / np.sqrt(piv.L11_nH * piv.L22_nH)
    piv.to_csv(OUT / "coupled_L100.csv", index=False)
    pd.set_option("display.width", 250)
    print(piv.round(4).to_string())


if __name__ == "__main__":
    main()
