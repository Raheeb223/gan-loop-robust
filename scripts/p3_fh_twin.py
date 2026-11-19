"""FastHenry 'twin' of the openEMS coupled models (Phase 3 cross-check).

openEMS represents copper as zero-thickness PEC sheets at the layer centre planes, so the
dielectric gap it sees is the centre-to-centre distance (d1 + cu), not the facing-surface
gap d1.  Exactly as in the Tier-1c canonical comparison, FastHenry is therefore re-run
with 2 um thin copper placed on the same centre planes (d1' = d1 + cu - 2 um, bump height
raised by cu/2 - 1 um), which makes the two solvers see identical conductor positions.
Output: results/phase3/fh_twin.csv
"""
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
import fhrun  # noqa: E402
import geom  # noqa: E402
import layouts  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
MIL, CU, THIN = 0.0254, 0.070, 0.002


def job(a):
    fam, gate = a
    p = layouts.Params(family=fam, gate=gate, T=31 * MIL, d1=4 * MIL + CU - THIN, cu=THIN,
                       bump=layouts.BUMP + CU / 2 - THIN / 2, pitch=0.15)
    r = fhrun.run(geom.to_fasthenry(layouts.build(p), (1e8, 1e9), max_fil_w=0.2, nh=1),
                  ROOT / "results" / "raw" / "phase3" / f"twin_{fam}_{gate}")
    rows = []
    for f, Z in zip(r["freqs"], r["Z"]):
        L = 0.5 * (Z.imag + Z.imag.T) / (2 * np.pi * f) * 1e9
        rows.append(dict(case=f"{fam}_{gate}", f_Hz=f, L11_nH=L[0, 0], L33_nH=L[2, 2], M13_nH=L[0, 2],
                         M12_nH=L[0, 1], runtime_s=r["runtime_s"]))
    return rows


if __name__ == "__main__":
    with ProcessPoolExecutor(max_workers=3) as ex:
        res = list(ex.map(job, [("optimal", "kelvin"), ("optimal", "common"), ("vertical", "kelvin")]))
    D = pd.DataFrame([r for rr in res for r in rr])
    D.to_csv(ROOT / "results" / "phase3" / "fh_twin.csv", index=False)
    print(D.round(4).to_string())
