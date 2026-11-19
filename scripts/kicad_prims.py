"""Manufacturable copper primitives (for KiCad) from the same layout Params.

The solver Geometry (layouts.build) and these KiCad primitives come from one
parametric description.  Differences, by design:
  * KiCad copper under each FET is split into interdigitated fingers matching
    the EPC2015C land pattern (D bars 3,5,7,9,11; S bars 4,6,8,10; 180 um wide,
    400 um pitch, 1362 um long).  The solver represents the device region by
    plane copper up to the die mid-line plus die bridges (see layouts.py).
  * KiCad vias are round (0.3 mm drill / 0.6 mm pad); solver uses square bars.
Output: JSON-serialisable dict with rects, vias, pads (all mm, solver coordinates).
"""
from __future__ import annotations

import json
import sys
from dataclasses import asdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
import layouts  # noqa: E402

NET_OF_PLANE = {"vin": "VIN", "bvin": "VIN", "sw": "SW", "gnd": "GND", "gpad": "GND", "ret": "GND",
                "l2": "GND", "bgnd": "GND"}
FINGER_W = 0.18
FINGER_PITCH = 0.40
FINGER_L = 1.362


def build(p: layouts.Params) -> dict:
    g = layouts.build(p)
    zt, zb, z2 = layouts._z(p)
    layer_of = {round(zt, 6): "F.Cu", round(z2, 6): "In1.Cu", round(zb, 6): "B.Cu"}
    rects, vias, fps = [], [], []
    # device x-extents (from FET bridge contact nodes "Ta0", "Sa0")
    devs = []
    for tag in ("T", "S"):
        xa = g.nodes[f"{tag}a0"][0] - 0.3
        devs.append((tag, xa, xa + layouts.FET_L))
    for pl in g.planes:
        net = NET_OF_PLANE[pl.name]
        layer = layer_of[round(pl.z, 6)]
        x0, x1 = pl.x0, pl.x1
        if layer == "F.Cu" and pl.y0 < p.W / 2 - 1e-6 and pl.y1 > -p.W / 2 + 1e-6:
            # clip top-layer planes (that overlap the die in y) to the die edges;
            # interdigitated fingers carry the current under the die.
            # Finger-end to opposite-plane clearance = (FET_L - FINGER_L)/2 = 0.137 mm.
            for tag, a, b in devs:
                mid = 0.5 * (a + b)
                if x0 < mid < x1 + 0.2 and x1 <= mid:      # plane ends at the die mid-line
                    x1 = min(x1, a)
                if x1 > mid > x0 - 0.2 and x0 >= mid:      # plane starts at the die mid-line
                    x0 = max(x0, b)
        rects.append(dict(layer=layer, net=net, x0=x0, y0=pl.y0, x1=x1, y1=pl.y1))
    # lateral: copper seam joining the return strip to the GND elbow (solver: seam bars)
    seam = [n for n in g.nodes if n.startswith("seamA")]
    if seam:
        ret = next(pl for pl in g.planes if pl.name == "ret")
        b_sr = devs[1][2]
        rects.append(dict(layer="F.Cu", net="GND", x0=ret.x1 - 0.05, y0=ret.y0, x1=b_sr + 0.05, y1=ret.y1))
    # interdigitated fingers + FET pads (die centred at y=0, bars along x)
    nets = {"T": ("VIN", "SW"), "S": ("SW", "GND")}
    for tag, a, b in devs:
        dnet, snet = nets[tag]
        pads = []
        ys = [-(4 * FINGER_PITCH) + i * FINGER_PITCH for i in range(9)]  # 9 bars, centred
        for i, y in enumerate(ys):
            is_drain = (i % 2 == 0)       # bars 3,5,7,9,11 -> drain
            net = dnet if is_drain else snet
            px0, px1 = a + (layouts.FET_L - FINGER_L) / 2, b - (layouts.FET_L - FINGER_L) / 2
            pads.append(dict(num=str(3 + i), net=net, x0=px0, y0=y - FINGER_W / 2, x1=px1, y1=y + FINGER_W / 2))
            # finger copper from the plane edge to the far end of the pad
            if is_drain:
                rects.append(dict(layer="F.Cu", net=net, x0=a - 0.05, y0=y - FINGER_W / 2, x1=px1, y1=y + FINGER_W / 2))
            else:
                rects.append(dict(layer="F.Cu", net=net, x0=px0, y0=y - FINGER_W / 2, x1=b + 0.05, y1=y + FINGER_W / 2))
        # gate pad (datasheet pad 1, 180 x 560 um) on its own net, placed just outside the power
        # copper (|y| > W/2) so the power-loop geometry is unchanged; needed by tools that
        # identify FET topology from D/S/G nets (e.g. dcdc-parasitics)
        yg = -p.W / 2 - 0.30
        gx = 0.5 * (a + b)
        pads.append(dict(num="1", net="GH" if tag == "T" else "GL",
                         x0=gx - 0.28, y0=yg - 0.09, x1=gx + 0.28, y1=yg + 0.09))
        fps.append(dict(ref=f"Q_{tag}", value="EPC2015C", pads=pads))
    # vias from via-row bars
    for bar in g.bars:
        n1, n2 = g.nodes[bar.n1], g.nodes[bar.n2]
        if bar.n1.startswith(("vg", "vc", "vr", "vv")):   # via rows (see layouts._via_row)
            vias.append(dict(x=n1[0], y=n1[1], net="GND" if bar.n1.startswith(("vg", "vc", "vr")) else "VIN",
                             drill=0.3, dia=0.6))
    # Cin bank: 0402 pads across the port gap (bundle positions)
    cap = []
    if g.port_box is not None:
        (sx, sy, sz), (ex, ey, ez), ax = g.port_box
        layer = "F.Cu" if abs(sz - zt) < 1e-6 else "B.Cu"
        if ax == "x":
            n = int(round((ey - sy) / 1.0))
            for i in range(n):
                yc = sy + (i + 0.5) * (ey - sy) / n
                cap.append(dict(ref=f"C{i+1}", value="1uF 0402", layer=layer, x=0.5 * (sx + ex), y=yc, rot=0))
        else:
            n = int(round((ex - sx) / 0.6))
            for i in range(n):
                xc = sx + (i + 0.5) * (ex - sx) / n
                cap.append(dict(ref=f"C{i+1}", value="1uF 0402", layer=layer, x=xc, y=0.5 * (sy + ey), rot=90))
    xs = [r["x0"] for r in rects] + [r["x1"] for r in rects]
    ys = [r["y0"] for r in rects] + [r["y1"] for r in rects]
    outline = dict(x0=min(xs) - 1.5, y0=min(ys) - 1.5, x1=max(xs) + 1.5, y1=max(ys) + 1.5)
    return dict(params=asdict(p), rects=rects, vias=vias, fets=fps, caps=cap, outline=outline,
                stackup=dict(T=p.T, d1=p.d1, cu=p.cu))


if __name__ == "__main__":
    out = Path(sys.argv[1])
    out.mkdir(parents=True, exist_ok=True)
    for fam in ("optimal", "vertical", "lateral", "lateral_stitched"):
        for T, d1 in ((31, 4), (31, 12), (62, 4), (62, 26)):
            p = layouts.Params(family=fam, T=T * 0.0254, d1=d1 * 0.0254)
            (out / f"{p.tag()}.json").write_text(json.dumps(build(p), indent=1))
    print("wrote", len(list(out.glob('*.json'))), "primitive files")
