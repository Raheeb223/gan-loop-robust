"""Reproducible pixel-level digitisation of Reusch & Strydom, APEC 2013 slides 14-17.

The charts are vector graphics (Excel) rendered at 4x with pdfium.  Calibration:
the plot area is filled with RGB (255,255,153); its bounding box spans exactly the
axis limits printed on the slide.  Series colours are exact: blue (0,0,255),
green (0,128,0), red (192,0,0).

* Scatter slides 15-17 (MEASURED power loss, overshoot, dv/dt vs EPC's L): data
  markers are filled discs, thicker than the fitted trend lines and the text strokes.
  We take the Euclidean distance transform of each colour mask; local maxima with an
  inscribed radius > R_MIN px are disc centres.  Overlapping discs are separated by
  non-maximum suppression with a minimum spacing.
* Line slide 14 (EPC SIMULATED L vs board thickness): for each colour, at a target
  board thickness we gather the vertical pixel runs in +-W columns and keep y-clusters
  that persist across >= 35 % of the columns (dashed/dotted lines persist, text does
  not).  Line style is not identified; clusters are reported by y order.
Outputs: data/published/digitised_apec13.csv (+ annotated overlay PNGs for checking).
Pixel resolution: x 1 px = (x-range)/(plot width px), typically 0.001-0.03 units.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import ndimage as ndi

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "data" / "published" / "sources" / "Reusch_APEC2013_slides.pdf"
OUT = ROOT / "data" / "published" / "digitise"
COL = {"optimal": (0, 0, 255), "lateral": (0, 128, 0), "vertical": (192, 0, 0)}
YELLOW = (255, 255, 153)
AXES = {  # slide: (xmin, xmax, ymin, ymax, quantity)  -- read from the slide's tick labels
    14: (20, 70, 0.2, 2.2, "L_sim_nH"),
    15: (0.2, 2.0, 3.0, 4.0, "power_loss_W"),
    16: (0.2, 2.0, 20, 100, "overshoot_pct"),
    17: (0.2, 2.0, 0, 10, "dvdt_V_per_ns"),
}


def render(page, scale=4):
    import pypdfium2 as pdfium
    img = pdfium.PdfDocument(str(SRC))[page - 1].render(scale=scale).to_pil().convert("RGB")
    return np.asarray(img).astype(int)


def mask(a, rgb, tol=10):
    return np.all(np.abs(a - np.array(rgb)) <= tol, axis=2)


def plot_box(a):
    # gridlines split the yellow plot area into cells: bridge them by dilation,
    # take the largest merged region, then the bbox of the yellow pixels inside it
    m = mask(a, YELLOW, 2)
    md = ndi.binary_dilation(m, structure=np.ones((15, 15), bool))
    lab, n = ndi.label(md)
    sizes = ndi.sum(md, lab, range(1, n + 1))
    k = int(np.argmax(sizes)) + 1
    ys, xs = np.nonzero(m & (lab == k))
    return xs.min(), xs.max(), ys.min(), ys.max()


def to_data(px, py, box, ax):
    x0, x1, y0, y1 = box
    xmin, xmax, ymin, ymax, _ = ax
    return (xmin + (px - x0) / (x1 - x0) * (xmax - xmin),
            ymax - (py - y0) / (y1 - y0) * (ymax - ymin))


def discs(a, rgb, box, r_min=None):
    x0, x1, y0, y1 = box
    m = mask(a, rgb)
    pad = 40
    sub = np.zeros_like(m)
    sub[y0 - pad:y1 + pad, x0 - pad:x1 + pad] = m[y0 - pad:y1 + pad, x0 - pad:x1 + pad]
    d = ndi.distance_transform_edt(sub)
    if r_min is None:
        # line half-width = mode of the ridge values; discs are clearly larger
        ridge = d[(d > 1) & (d == ndi.maximum_filter(d, 3))]
        r_line = np.percentile(ridge, 60)
        r_min = max(1.6 * r_line, 6)
    peaks = (d == ndi.maximum_filter(d, 9)) & (d > r_min)
    ys, xs = np.nonzero(peaks)
    order = np.argsort(-d[ys, xs])
    keep = []
    for i in order:
        if all((xs[i] - xs[j]) ** 2 + (ys[i] - ys[j]) ** 2 > (1.2 * r_min) ** 2 for j in keep):
            keep.append(i)
    return [(xs[i], ys[i], d[ys[i], xs[i]]) for i in keep], r_min


def strip_glyphs(m, box, line_h_max=35, glyph_w_max=120):
    """Remove text glyphs: connected components taller than a line stroke but narrow.
    Chart strokes are either long (solid lines, dashes) or small (dots)."""
    x0, x1, y0, y1 = box
    sub = m.copy()
    sub[:y0 - 5] = False; sub[y1 + 5:] = False; sub[:, :x0 - 5] = False; sub[:, x1 + 5:] = False
    lab, n = ndi.label(sub)
    for k, sl in enumerate(ndi.find_objects(lab), start=1):
        h = sl[0].stop - sl[0].start
        w = sl[1].stop - sl[1].start
        if h > line_h_max and w < glyph_w_max:
            sub[sl][lab[sl] == k] = False
    return sub


def line_clusters(a, rgb, box, ax, x_target, half_w=12, min_frac=0.2, ytol=6):
    """Return [(y_data, duty)] for each line crossing column x_target.  duty = fraction of
    columns in the window where the stroke is present (solid ~1, dashed < 1, dotted lower)."""
    x0, x1, y0, y1 = box
    xmin, xmax = ax[0], ax[1]
    pc = int(round(x0 + (x_target - xmin) / (xmax - xmin) * (x1 - x0)))
    m = strip_glyphs(mask(a, rgb), box)
    pts = []
    cols = list(range(pc - half_w, pc + half_w + 1))
    for c in cols:
        lab, n = ndi.label(m[y0:y1 + 1, c])
        for k in range(1, n + 1):
            ys = np.nonzero(lab == k)[0]
            if len(ys) < 60:
                pts.append((c, y0 + ys.mean()))
    if not pts:
        return []
    pts.sort(key=lambda t: t[1])
    groups, cur = [], [pts[0]]
    for p in pts[1:]:
        if p[1] - cur[-1][1] <= ytol:
            cur.append(p)
        else:
            groups.append(cur); cur = [p]
    groups.append(cur)
    res = []
    for g in groups:
        duty = len({c for c, _ in g}) / len(cols)
        if duty >= min_frac:
            # value at the target column: linear fit through the stroke pixels
            cs = np.array([c for c, _ in g]); ys = np.array([y for _, y in g])
            yc = np.polyval(np.polyfit(cs, ys, 1), pc) if len(set(cs)) > 2 else ys.mean()
            res.append((round(to_data(pc, yc, box, ax)[1], 3), round(duty, 2)))
    return res


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    rows = []
    for slide in (15, 16, 17):
        a = render(slide)
        box = plot_box(a)
        ax = AXES[slide]
        for fam, rgb in COL.items():
            pts, rmin = discs(a, rgb, box)
            for (px, py, r) in pts:
                if not (box[0] - 5 <= px <= box[1] + 5 and box[2] - 5 <= py <= box[3] + 5):
                    continue
                x, y = to_data(px, py, box, ax)
                rows.append(dict(slide=slide, kind="measured", family=fam, L_epc_nH=round(x, 3),
                                 quantity=ax[4], value=round(y, 3), px=int(px), py=int(py), r_px=round(r, 1),
                                 occluded=bool(r < 16.0)))   # full markers: r = 17.5-18.4 px
        print(f"slide {slide}: plot box {box}, {sum(1 for r in rows if r['slide'] == slide)} markers", flush=True)
    a = render(14)
    box = plot_box(a)
    for fam, rgb in COL.items():
        for xt, hw in ((31.5, 12), (61.0, 12), (57.5, 60)):   # 57.5 mil: where the 26-mil dotted curves exist
            for y, duty in line_clusters(a, rgb, box, AXES[14], xt, half_w=hw):
                rows.append(dict(slide=14, kind="EPC simulation", family=fam, board_thickness_mil=xt,
                                 quantity="L_sim_nH", value=y, duty=duty))
    print(f"slide 14: plot box {box}", flush=True)
    d = pd.DataFrame(rows)
    d.to_csv(ROOT / "data" / "published" / "digitised_apec13.csv", index=False)
    pd.set_option("display.width", 200)
    print(d.drop(columns=[c for c in ("px", "py") if c in d]).to_string())


if __name__ == "__main__":
    main()
