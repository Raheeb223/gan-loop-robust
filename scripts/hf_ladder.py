"""HF (ringing-band) convergence of the layout inductance in the thickness discretisation.

At 100 MHz the skin depth in copper is 6.6 um << 70 um (2 oz), so the number of
filaments across the conductor thickness (nh, with geometric ratio 2 toward the
surfaces) controls L(100 MHz).  Ladder nh = 1, 3, 5, 7 at in-plane pitch 0.15 mm.
Output: results/tables/hf_ladder.csv
"""
from __future__ import annotations

import sys
from concurrent.futures import ProcessPoolExecutor
from dataclasses import asdict
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
import fhrun  # noqa: E402
import geom  # noqa: E402
import layouts  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "results" / "raw" / "hf_ladder"
TAB = ROOT / "results" / "tables"


def job(a):
    fam, nh = a
    p = layouts.Params(family=fam, pitch=0.15)
    tag = f"{p.tag()}_nh{nh}"
    r = fhrun.run(geom.to_fasthenry(layouts.build(p), (1e6, 1e9), max_fil_w=0.2, nh=nh), RAW / tag)
    out = dict(asdict(p), nh=nh, runtime_s=r["runtime_s"])
    for f, L, R in zip(r["freqs"], r["L"], r["R"]):
        out[f"L_{f:.0e}_nH"] = L * 1e9
        out[f"R_{f:.0e}_mohm"] = R * 1e3
    print(fam, nh, f"L100M={out['L_1e+08_nH']:.4f}", f"t={r['runtime_s']:.0f}s", flush=True)
    return out


if __name__ == "__main__":
    jobs = [(f, nh) for nh in (1, 3, 5, 7) for f in ("optimal", "vertical", "lateral", "lateral_stitched")]
    with ProcessPoolExecutor(max_workers=2) as ex:
        rows = list(ex.map(job, jobs))
    pd.DataFrame(rows).to_csv(TAB / "hf_ladder.csv", index=False)
