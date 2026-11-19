"""Attach calibrated-model overshoot predictions (12 V, 20 A, EPC2015C) to the design map."""
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))
import pandas as pd
import spice_sim as ss

ROOT = Path(__file__).resolve().parents[1]
TAB = ROOT / "results" / "tables"
N_CAPS = {"optimal": 4, "vertical": 4, "lateral": 2, "lateral_stitched": 2}


def run(args):
    L, rup, rd = args
    o = ss.run_bench(ss.Bench(L_loop=L * 1e-9, r_pullup=rup, R_loop=rd, t_stop=150e-9))
    return o["overshoot_pct"], o["f_ring"] / 1e6


if __name__ == "__main__":
    d = pd.read_csv(TAB / "design_map.csv")
    cal = pd.read_csv(TAB / "circuit_validation.csv").iloc[0]
    d["L_spice_nH"] = d.L_nom_nH + 0.35 / d.family.map(N_CAPS)
    jobs = [(L, cal.r_pullup, cal.R_damp) for L in d.L_spice_nH]
    with ProcessPoolExecutor(max_workers=4) as ex:
        res = list(ex.map(run, jobs))
    d["overshoot_pred_pct"] = [r[0] for r in res]
    d["f_ring_pred_MHz"] = [r[1] for r in res]
    d.to_csv(TAB / "design_map.csv", index=False)
    print(d[["family", "stackup", "L_nom_nH", "L_p2_5_nH", "L_p97_5_nH", "overshoot_pred_pct", "f_ring_pred_MHz"]].round(2).to_string())
