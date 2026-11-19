"""Sensitivity of the ringing-derived PCB inductance of Letellier et al. (TPEL 2019)
to the assumed input tolerance.

Same decomposition as ring_uncertainty.py part (ii): L_stray from the ringing
frequency and Coss(50 V) = 280 pF, minus 2 x 0.1 nH package and 0.87 nH capacitor
inductance (published L_pcb = 0.38 nH). Each of the four inputs is perturbed by an
independent uniform factor in [1-u, 1+u] for u = 10, 20, 30 %. Plain language: we
ask how wide the answer becomes if each number the authors had to assume is only
known to within u percent. Output: results/revision/ring_sensitivity.csv
"""
import pathlib
import numpy as np
import pandas as pd

OUT = pathlib.Path(__file__).resolve().parents[1] / "results" / "revision"
OUT.mkdir(parents=True, exist_ok=True)
N, SEED = 200_000, 1

f_ring = 1 / (2 * np.pi * np.sqrt(1.45e-9 * 280e-12))
rows = []
for u_pct in (10, 20, 30):
    rng = np.random.default_rng(SEED)
    u = lambda: rng.uniform(1 - u_pct / 100, 1 + u_pct / 100, N)
    L_stray = 1 / ((2 * np.pi * f_ring) ** 2 * (280e-12 * u()))
    L_pcb = (L_stray - 2 * 0.1e-9 * u() - 0.87e-9 * u()) * 1e9
    lo, med, hi = np.percentile(L_pcb, [2.5, 50, 97.5])
    rows.append(dict(tolerance_pct=u_pct, L_pcb_median_nH=med, L_pcb_p2_5_nH=lo,
                     L_pcb_p97_5_nH=hi, rel_halfwidth_pct=100 * (hi - lo) / 2 / 0.38,
                     p_negative_pct=100 * np.mean(L_pcb < 0), n=N, seed=SEED))
df = pd.DataFrame(rows)
p = OUT / "ring_sensitivity.csv"
df.to_csv(p, index=False)
print(df.to_string(index=False))
