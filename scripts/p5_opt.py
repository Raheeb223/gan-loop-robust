"""Phase 5 (C3): tolerance-robust layout optimisation with NSGA-II, verified and exported.

Plain language: we let a genetic algorithm (NSGA-II, pymoo) search the layout choices of
the optimal power loop on a JLCPCB 1.6 mm board.  Each candidate is scored on
  f1  robust loop inductance: the 95th percentile of L(100 MHz) over the manufacturing
      tolerances (256 common random tolerance samples through the Phase-4 surrogate),
      raised by the 95 % split-conformal margin so it is a guaranteed-conservative bound;
  f2  minus the device-side copper area W (l_vin + l_sw + l_gnd) - a heat-spreading
      proxy: GaN chip-scale devices are cooled mainly through the board copper, so more
      copper around them helps, but it also lengthens the loop.  (A proxy, not a thermal
      simulation.)
Decision variables: JLCPCB prepreg (1080/3313/2116/7628 -> d1 = 76/99/116/210 um),
vias per row n_via (3-7), l_vin, l_sw, l_gnd.
Then: (i) three Pareto designs (lowest L, knee, largest copper) are re-checked with 24
direct FastHenry tolerance Monte Carlo runs each; (ii) their KiCad boards are generated,
DRC-checked and exported (Gerbers + drill).
Outputs: results/phase5/pareto.csv, verify.csv, picks.csv, kicad_drc.csv, fig_p5_pareto.(pdf|png);
kicad/boards_p5/<tag>/...
"""
from __future__ import annotations

import json
import pickle
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor
from dataclasses import replace
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
import p4_campaign as C  # noqa: E402
import p4_uq  # noqa: E402,F401  (needed to unpickle the surrogate)
from p4_uq import X_at, sample_tol  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "phase5"
PREPREG = [0.076, 0.099, 0.116, 0.210]
PREPREG_NAME = ["1080", "3313", "2116", "7628"]
W = 4.1


def load_sur():
    import __main__
    __main__.Surrogate = p4_uq.Surrogate      # the model was pickled by running p4_uq.py as __main__
    with open(ROOT / "results" / "phase4" / "gp_model.pkl", "rb") as fh:
        return pickle.load(fh)


def decode(x):
    return dict(d1=PREPREG[int(np.clip(round(x[0]), 0, 3))], n_via=float(int(np.clip(round(x[1]), 3, 7))),
                l_vin=float(x[2]), l_sw=float(x[3]), l_gnd=float(x[4]))


def make_problem(sur, tol):
    from pymoo.core.problem import ElementwiseProblem
    q95 = sur.q[0.05]

    class P(ElementwiseProblem):
        def __init__(self):
            super().__init__(n_var=5, n_obj=2, xl=np.array([0, 3, 0.8, 0.6, 0.8]),
                             xu=np.array([3, 7, 1.8, 1.6, 1.8]))

        def _evaluate(self, x, out, *a, **k):
            d = decode(x)
            mu, sd = sur.predict(X_at(d, tol))
            L95 = np.percentile(np.exp(mu + q95 * sd), 95)
            out["F"] = [L95, -W * (d["l_vin"] + d["l_sw"] + d["l_gnd"])]
    return P()


def verify_job(a):
    tag, d, i, tolx = a
    x = dict(d)
    x.update(zip(C.TOL, tolx))
    r = C.run((f"p5_{tag}", i, x))
    return dict(pick=tag, **r)


def kicad(picks):
    import kicad_all as K
    import kicad_prims
    rows = []
    for _, r in picks.iterrows():
        p = C.params(dict(d1=r.d1, n_via=r.n_via, l_vin=r.l_vin, l_sw=r.l_sw, l_gnd=r.l_gnd, t_d1=1, t_cu=1, t_cuin=1,
                          hole=0.3, bump=0.1, dx_die=0.0))
        tag = f"p5_{r.pick}_{PREPREG_NAME[PREPREG.index(r.d1)]}_v{int(r.n_via)}"
        bd = ROOT / "kicad" / "boards_p5" / tag
        (bd / "fab").mkdir(parents=True, exist_ok=True)
        js = bd / f"{tag}.json"
        js.write_text(json.dumps(kicad_prims.build(p), indent=1))
        pcb = bd / f"{tag}.kicad_pcb"
        K.sh(K.KPY, ROOT / "scripts" / "kicad_build.py", js, pcb)
        K.sh(K.CLI, "pcb", "drc", "--format", "json", "--severity-error", "-o", bd / "drc.json", pcb)
        dj = json.load(open(bd / "drc.json"))
        K.sh(K.CLI, "pcb", "export", "gerbers", "-o", str(bd / "fab") + "\\", pcb)
        K.sh(K.CLI, "pcb", "export", "drill", "-o", str(bd / "fab") + "\\", pcb)
        K.sh(K.CLI, "pcb", "export", "svg", "--layers", "F.Cu,Edge.Cuts", "--mode-single", "--fit-page-to-board",
             "-o", bd / f"{tag}_F.svg", pcb)
        rows.append(dict(pick=r.pick, tag=tag, drc_violations=len(dj.get("violations", [])),
                         unconnected=len(dj.get("unconnected_items", [])), fab_files=len(list((bd / "fab").iterdir()))))
    return pd.DataFrame(rows)


def main():
    from pymoo.algorithms.moo.nsga2 import NSGA2
    from pymoo.optimize import minimize
    OUT.mkdir(parents=True, exist_ok=True)
    sur = load_sur()
    tol = sample_tol(256, np.random.default_rng(31))
    res = minimize(make_problem(sur, tol), NSGA2(pop_size=80), ("n_gen", 120), seed=5, verbose=False)
    rows = []
    for x, F in zip(res.X, res.F):
        d = decode(x)
        rows.append(dict(**d, prepreg=PREPREG_NAME[PREPREG.index(d["d1"])], L95_robust_nH=F[0], copper_area_mm2=-F[1]))
    P = pd.DataFrame(rows).drop_duplicates().sort_values("L95_robust_nH").reset_index(drop=True)
    P.to_csv(OUT / "pareto.csv", index=False)
    print(P.round(3).to_string(), flush=True)
    # reference: the nominal design and the JLCPCB-default (7628) version of it
    refs = []
    for name, d in (("nominal_3313", dict(C.NOMINAL)), ("nominal_7628", dict(C.NOMINAL, d1=0.210))):
        mu, sd = sur.predict(X_at(d, tol))
        refs.append(dict(name=name, **d, L95_robust_nH=np.percentile(np.exp(mu + sur.q[0.05] * sd), 95),
                         copper_area_mm2=W * (d["l_vin"] + d["l_sw"] + d["l_gnd"])))
    pd.DataFrame(refs).to_csv(OUT / "references.csv", index=False)
    # picks: lowest L, knee (max distance to the line between the extremes, normalised), largest copper
    f1 = (P.L95_robust_nH - P.L95_robust_nH.min()) / np.ptp(P.L95_robust_nH)
    f2 = (P.copper_area_mm2.max() - P.copper_area_mm2) / np.ptp(P.copper_area_mm2)
    knee = int(np.argmin(np.hypot(f1, f2)))
    picks = pd.concat([P.iloc[[0]].assign(pick="minL"), P.iloc[[knee]].assign(pick="knee"),
                       P.iloc[[len(P) - 1]].assign(pick="maxCu")]).reset_index(drop=True)
    picks.to_csv(OUT / "picks.csv", index=False)
    # direct FastHenry verification (24 tolerance samples per pick)
    tolv = sample_tol(24, np.random.default_rng(41))
    jobs = [(r.pick, dict(d1=r.d1, n_via=r.n_via, l_vin=r.l_vin, l_sw=r.l_sw, l_gnd=r.l_gnd), i, tolv[i])
            for _, r in picks.iterrows() for i in range(len(tolv))]
    with ProcessPoolExecutor(max_workers=4) as ex:
        V = pd.DataFrame(list(ex.map(verify_job, jobs)))
    V.to_csv(OUT / "verify_runs.csv", index=False)
    summ = []
    for _, r in picks.iterrows():
        v = V[(V.pick == r.pick) & (V.err.fillna("") == "")]
        d = dict(d1=r.d1, n_via=r.n_via, l_vin=r.l_vin, l_sw=r.l_sw, l_gnd=r.l_gnd)
        lo, mid, hi = sur.interval(np.array([[v_[k] for k in C.VARS] for _, v_ in v.iterrows()]), 0.05)
        summ.append(dict(pick=r.pick, prepreg=r.prepreg, n_via=r.n_via, l_vin=r.l_vin, l_sw=r.l_sw, l_gnd=r.l_gnd,
                         L95_robust_pred_nH=r.L95_robust_nH, fh_mean_nH=v.L_nH.mean(), fh_p95_nH=np.percentile(v.L_nH, 95),
                         fh_max_nH=v.L_nH.max(), n_runs=len(v),
                         conformal_coverage=float(np.mean((v.L_nH.values >= lo) & (v.L_nH.values <= hi))),
                         gp_rmse_pct=float(100 * np.sqrt(np.mean((mid / v.L_nH.values - 1) ** 2)))))
    S = pd.DataFrame(summ)
    S.to_csv(OUT / "verify.csv", index=False)
    print(S.round(4).to_string(), flush=True)
    Kd = kicad(picks)
    Kd.to_csv(OUT / "kicad_drc.csv", index=False)
    print(Kd.to_string(), flush=True)
    figure(P, picks, pd.DataFrame(refs), S)


def figure(P, picks, refs, S):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(4.2, 3.1), constrained_layout=True)
    for pp, g in P.groupby("prepreg"):
        ax.plot(g.copper_area_mm2, g.L95_robust_nH, "o", ms=3.5, label=f"Pareto, {pp}")
    ax.plot(picks.copper_area_mm2, picks.L95_robust_nH, "k*", ms=11, label="selected designs")
    ax.plot(picks.copper_area_mm2, S.fh_p95_nH, "rx", ms=8, label="FastHenry 95th pct (24 runs)")
    for _, r in refs.iterrows():
        ax.plot(r.copper_area_mm2, r.L95_robust_nH, "s", mfc="none", ms=7,
                label=r["name"].replace("_", " "))
    ax.set_xlabel("copper area next to the devices (mm$^2$)")
    ax.set_ylabel("robust $L_{95}$ at 100 MHz (nH)")
    ax.legend(fontsize=6)
    for ext in ("pdf", "png"):
        fig.savefig(OUT / f"fig_p5_pareto.{ext}", dpi=200)


if __name__ == "__main__":
    main()
