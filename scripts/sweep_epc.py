"""Tier-2: reproduce the WP010 / Reusch Fig. 7b study with FastHenry.

Sweep: family x board thickness (31, 39, 47, 62 mil) x inner-layer distance
(4, 12, 26 mil), 2 oz copper on all layers (WP010 builds).  Plus a mesh
convergence ladder on WP010 'Design 1' (31 mil / 4 mil).

Outputs results/tables/epc_sweep.csv, results/tables/epc_mesh.csv
"""
from __future__ import annotations

import sys
from concurrent.futures import ProcessPoolExecutor
from dataclasses import asdict, replace
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
import fhrun  # noqa: E402
import geom  # noqa: E402
import layouts  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "results" / "raw" / "epc_sweep"
TAB = ROOT / "results" / "tables"
MIL = 0.0254
FAMS = ["optimal", "vertical", "lateral", "lateral_stitched"]
T_MIL = [31, 39, 47, 62]
D1_MIL = [4, 12, 26]
# converged settings (results/tables/hf_ladder.csv): nh=5 changes L(100 MHz) < 0.3 % vs nh=7
PITCH, FW, NH = 0.15, 0.2, 5
DO_MESH = "--mesh" in sys.argv


def job(args):
    p, fw, nh, tag = args
    g = layouts.build(p)
    try:
        r = fhrun.run(geom.to_fasthenry(g, (1e6, 1e9), max_fil_w=fw, nh=nh), RAW / tag)
        out = dict(asdict(p), fw=fw, nh=nh, tag=tag, runtime_s=r["runtime_s"], err="")
        for f, L, R in zip(r["freqs"], r["L"], r["R"]):
            out[f"L_{f:.0e}_nH"] = L * 1e9
            out[f"R_{f:.0e}_mohm"] = R * 1e3
        return out
    except Exception as ex:
        return dict(asdict(p), fw=fw, nh=nh, tag=tag, err=str(ex))


def main():
    jobs = []
    for fam in FAMS:
        for t in T_MIL:
            for d in D1_MIL:
                T = t * MIL
                d1 = d * MIL
                # need room for two inner dielectrics + 4 copper layers + a core >= 0.1 mm
                if 2 * (d1 + 0.07) + 0.1 > T:
                    continue
                p = layouts.Params(family=fam, T=T, d1=d1, pitch=PITCH)
                jobs.append((p, FW, NH, f"{p.tag()}_nh{NH}"))
    mesh = []
    for fam in FAMS:
        for pitch, fw, nh in [(0.4, 0.6, 1), (0.2, 0.3, 1), (0.1, 0.15, 2), (0.05, 0.1, 3)]:
            p = layouts.Params(family=fam, pitch=pitch)
            mesh.append((p, fw, nh, f"mesh_{p.tag()}_fw{fw}_nh{nh}"))
    with ProcessPoolExecutor(max_workers=3) as ex:
        res = list(ex.map(job, jobs))
        pd.DataFrame(res).to_csv(TAB / "epc_sweep_hf.csv", index=False)
        print("sweep done", len(res), "failed", sum(bool(r.get("err")) for r in res), flush=True)
        if not DO_MESH:
            return
        res_m = list(ex.map(job, mesh))
        pd.DataFrame(res_m).to_csv(TAB / "epc_mesh.csv", index=False)
    print("mesh done", flush=True)


if __name__ == "__main__":
    main()
