"""openEMS (FDTD) loop-impedance extraction for a solver-neutral Geometry.

Method
------
* Every Bar of the Geometry becomes a PEC (or finite-conductivity) box of its true
  thickness; the mesh has lines on every copper face, so the conductor thickness
  is resolved exactly like FastHenry's.
* The port is an openEMS LumpedPort: a sheet that spans the gap between the
  two port nodes over the full conductor width (direction = the port vector).
* Excitation: Gaussian pulse 0..f_max.  Z(f) = V/I at the port.
* Loop inductance: L(f) = Im Z / (2 pi f), reported on a low-frequency plateau
  where the structure is electrically small (we report f and check flatness).

Optional dielectric: a slab between given z limits with eps_r (vacuum by default,
so that the comparison with magneto-quasi-static FastHenry is like-for-like).
"""
from __future__ import annotations

import os
import shutil
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
_OEMS = Path(os.environ.get("OPENEMS_INSTALL_PATH", ROOT / "tools" / "openEMS_dl" / "openEMS"))
os.add_dll_directory(str(_OEMS))
os.environ.setdefault("OPENEMS_INSTALL_PATH", str(_OEMS))

from CSXCAD import ContinuousStructure  # noqa: E402
from CSXCAD.SmoothMeshLines import SmoothMeshLines  # noqa: E402
from openEMS import openEMS  # noqa: E402
from openEMS.physical_constants import C0  # noqa: E402

UNIT = 1e-3  # mm
_THIRDS: list = []   # (axis, edge coordinate, outward sign) for thirds-rule meshing
_PASSIVE: list = []  # passive (non-excited) lumped ports of the current run
_MINCELL: list = [None]
_FINEZ: list = [None]
_PROTECT: dict = {0: set(), 1: set(), 2: set()}   # coordinates of zero-thickness sheets


def _merge_close(c, dmin, protect=()):
    """Replace clusters of sorted lines closer than dmin by one line: a protected
    coordinate (the plane of a zero-thickness sheet, which must stay on a mesh line or
    openEMS drops the sheet) if the cluster has one, else the cluster mean."""
    c = np.sort(np.asarray(c, float))
    prot = {round(x, 6) for x in protect}

    def pick(grp):
        p = [x for x in grp if round(x, 6) in prot]
        return p[0] if p else float(np.mean(grp))
    out, grp = [], [c[0]]
    for x in c[1:]:
        if x - grp[-1] < dmin:
            grp.append(x)
        else:
            out.append(pick(grp)); grp = [x]
    out.append(pick(grp))
    return np.array(out)


def run_geometry(g, workdir: Path, *, f_max: float = 2e9, f_eval=None, margin: float = 20.0,
                 max_cell: float | None = None, fine: float | None = None, kappa: float | None = None,
                 dielectric: tuple[float, float, float] | None = None, end_db: float = -40,
                 max_steps: int = 2_000_000, port_R: float = 50.0, boundary: str = "PML_8",
                 keep_files: bool = False, verbose: int = 0, sheet: bool = True,
                 threads: int = 0, thirds_y: bool = False, excite_box=None, passive_boxes=(),
                 passive_R: float = 1e4, min_cell: float | None = None,
                 fine_z: float | None = None) -> dict:
    """Simulate geometry `g` and return Z(f) at its first port.

    margin   : air margin around the copper bounding box (mm)
    max_cell : maximum cell size (mm); default lambda_min/20 capped at margin/4
    fine     : maximum cell size inside the copper bounding box (mm)
    kappa    : conductivity (S/m) for copper; None => PEC
    dielectric: (z_lo, z_hi, eps_r) slab across the whole domain, optional
    excite_box: (start, stop, axis) lumped port to excite instead of g.port_box (multi-port use)
    fine_z   : maximum z cell inside the copper bounding box (mm), for thin plane-pair gaps
    min_cell : merge mesh lines closer than this (mm); avoids micrometre cells (and femtosecond
               time steps) where unrelated copper edges nearly coincide.  None = off.
    passive_boxes: further (start, stop, axis) lumped ports, not excited, terminated in
               passive_R (>> |Z| of the loop, i.e. practically open); the result then also
               holds Z_transfer = V_passive / I_excited, the open-circuit transfer impedance.
    """
    workdir = Path(workdir)
    if workdir.exists():
        shutil.rmtree(workdir)
    workdir.mkdir(parents=True)
    f_eval = np.asarray(f_eval if f_eval is not None else np.linspace(f_max / 200, f_max / 4, 60))

    fdtd = openEMS(NrTS=max_steps, EndCriteria=10 ** (end_db / 10))
    if threads:
        fdtd.SetNumberOfThreads(threads)
    fdtd.SetGaussExcite(f_max / 2, f_max / 2)
    fdtd.SetBoundaryCond([boundary] * 6)
    csx = ContinuousStructure()
    fdtd.SetCSX(csx)
    grid = csx.GetGrid()
    grid.SetDeltaUnit(UNIT)

    if kappa is None:
        cu = csx.AddMetal("cu")
    else:
        cu = csx.AddMaterial("cu", kappa=kappa)

    _THIRDS.clear()
    for a in range(3):
        _PROTECT[a].clear()
    boxes = g.boxes()
    if sheet:
        # Collapse each bar's thickness to its centre plane (zero-thickness PEC sheet).
        # FDTD timestep is then set by the in-plane/gap mesh, not by 35 um copper.
        coll = []
        for (lo, hi), b in zip(boxes, g.bars):
            size = hi - lo
            ax = int(np.argmin(np.abs(size - b.h)))
            c = 0.5 * (lo[ax] + hi[ax])
            lo = lo.copy(); hi = hi.copy(); lo[ax] = hi[ax] = c
            coll.append((lo, hi))
            _PROTECT[ax].add(round(float(c), 6))
        boxes = coll
    if thirds_y:
        # free edges along y of horizontal sheets (in-plane width edges)
        _THIRDS.clear()
        for lo, hi in boxes:
            if abs(hi[2] - lo[2]) < 1e-12 and hi[1] - lo[1] > 1e-9:
                _THIRDS.extend([(1, lo[1], -1), (1, hi[1], +1)])
    for p in g.planes:
        if sheet:
            boxes.append((np.array([p.x0, p.y0, p.z]), np.array([p.x1, p.y1, p.z])))
            _PROTECT[2].add(round(float(p.z), 6))
        else:
            boxes.append((np.array([p.x0, p.y0, p.z - p.thick / 2]), np.array([p.x1, p.y1, p.z + p.thick / 2])))
    # snap every coordinate to a 1e-6 mm grid: a zero-thickness PEC sheet whose end
    # misses a mesh line by floating-point round-off is silently dropped by openEMS
    # ("Unused primitive"), which can open-circuit the loop.
    boxes = [(np.round(lo, 6), np.round(hi, 6)) for lo, hi in boxes]
    lines = {0: [], 1: [], 2: []}
    lo_all = np.min([b[0] for b in boxes], axis=0)
    hi_all = np.max([b[1] for b in boxes], axis=0)
    for lo, hi in boxes:
        cu.AddBox(lo.tolist(), hi.tolist(), priority=10)
        for a in range(3):
            lines[a] += [lo[a], hi[a]]

    if dielectric is not None:
        zlo, zhi, er = dielectric
        diel = csx.AddMaterial("diel", epsilon=er)
        diel.AddBox([lo_all[0] - margin, lo_all[1] - margin, zlo], [hi_all[0] + margin, hi_all[1] + margin, zhi],
                    priority=1)
        lines[2] += [zlo, zhi]

    passive = []
    for k, (start, stop, ax) in enumerate(passive_boxes):
        start = np.round(np.array(start, float), 6); stop = np.round(np.array(stop, float), 6)
        passive.append(fdtd.AddLumpedPort(2 + k, passive_R, start.tolist(), stop.tolist(), ax, 0.0, priority=5))
        for a in range(3):
            lines[a] += [start[a], stop[a]]
    _PASSIVE[:] = passive
    _MINCELL[0] = min_cell
    _FINEZ[0] = fine_z
    pbox = excite_box if excite_box is not None else g.port_box
    if pbox is not None:
        start, stop, ax = pbox
        start = np.round(np.array(start, float), 6); stop = np.round(np.array(stop, float), 6)
        port = fdtd.AddLumpedPort(1, port_R, start.tolist(), stop.tolist(), ax, 1.0, priority=5)
        for a in range(3):
            lines[a] += [start[a], stop[a]]
        return _finish(fdtd, csx, grid, port, lines, lo_all, hi_all, margin, max_cell, fine, f_max, f_eval,
                       workdir, verbose, keep_files)
    # port: sheet between the two port nodes, spanning the width of the conductor at n_plus
    p_n, m_n = g.ports[0]
    pp = np.array(g.nodes[p_n]); pm = np.array(g.nodes[m_n])
    bar_p = next(b for b in g.bars if p_n in (b.n1, b.n2))
    bar_m = next(b for b in g.bars if m_n in (b.n1, b.n2))
    d = pm - pp
    pdir = int(np.argmax(np.abs(d)))
    # width direction of the bar touching the + node (in-plane, perpendicular to bar axis)
    seg = np.array(g.nodes[bar_p.n2]) - np.array(g.nodes[bar_p.n1])
    wv = np.abs(np.array(bar_p.wdir if bar_p.wdir else (-seg[1], seg[0], 0.0), float))
    wv /= np.linalg.norm(wv)
    wax = int(np.argmax(wv))
    start = pp.copy(); stop = pm.copy()
    # move the sheet ends from conductor centre-lines to conductor surfaces
    sgn = np.sign(d[pdir])
    if not sheet:
        start[pdir] += sgn * bar_p.h / 2
        stop[pdir] -= sgn * bar_m.h / 2
    start[wax] -= bar_p.w / 2
    stop[wax] += bar_p.w / 2
    start = np.round(start, 6); stop = np.round(stop, 6)
    port = fdtd.AddLumpedPort(1, port_R, start.tolist(), stop.tolist(), "xyz"[pdir], 1.0, priority=5)
    for a in range(3):
        lines[a] += [start[a], stop[a]]
    return _finish(fdtd, csx, grid, port, lines, lo_all, hi_all, margin, max_cell, fine, f_max, f_eval,
                   workdir, verbose, keep_files)


def _finish(fdtd, csx, grid, port, lines, lo_all, hi_all, margin, max_cell, fine, f_max, f_eval,
            workdir, verbose, keep_files):
    lam_min = C0 / f_max / UNIT
    if max_cell is None:
        max_cell = min(lam_min / 20, margin / 4)
    if fine is None:
        fine = max_cell
    thirds = list(_THIRDS)
    for a in range(3):
        core = np.unique(np.round(lines[a], 9))
        if thirds:
            # openEMS 'thirds rule': for a free metal edge at e (outward sign s) put
            # lines at e - s*h/3 (inside) and e + s*2h/3 (outside), none on the edge.
            eds = [(e, s) for (ax, e, s) in thirds if ax == a]
            if eds:
                drop = {round(e, 9) for e, _ in eds}
                core = np.array([c for c in core if round(c, 9) not in drop])
                add = [e - s * fine / 3 for e, s in eds] + [e + s * 2 * fine / 3 for e, s in eds]
                core = np.unique(np.concatenate([core, add]))
        if _MINCELL[0]:
            core = _merge_close(core, _MINCELL[0], _PROTECT[a])
        # enforce `fine` resolution inside the copper bounding box
        fa = min(fine, _FINEZ[0]) if (a == 2 and _FINEZ[0]) else fine   # finer z mesh across thin dielectrics
        core = SmoothMeshLines(np.unique(np.concatenate([core, [lo_all[a], hi_all[a]]])), fa, 1.3)
        full = np.concatenate([core, [lo_all[a] - margin, hi_all[a] + margin]])
        full = SmoothMeshLines(np.unique(full), max_cell, 1.3)
        grid.SetLines("xyz"[a], full)

    ncell = int(np.prod([len(grid.GetLines(a)) for a in "xyz"]))
    t0 = time.perf_counter()
    fdtd.Run(str(workdir), verbose=verbose, cleanup=True)
    rt = time.perf_counter() - t0
    port.CalcPort(str(workdir), f_eval)
    Z = port.uf_tot / port.if_tot
    L = Z.imag / (2 * np.pi * f_eval)
    zt = []
    for pp in _PASSIVE:
        pp.CalcPort(str(workdir), f_eval)
        zt.append(pp.uf_tot / port.if_tot)
    if not keep_files:
        for p in workdir.glob("*.h5"):
            p.unlink()
    return dict(f=f_eval, Z=Z, L=L, R=Z.real, runtime_s=rt, ncells=ncell, Z_transfer=zt)


def fit_L0(f, L_app, f_hi: float | None = None):
    """Remove the self-resonance bias from apparent inductance.

    Below the first resonance a shorted loop behaves as jwL0 / (1 - (f/fr)^2),
    so 1/L_app = 1/L0 - f^2/(L0 fr^2) is linear in f^2.  Least-squares fit
    over f <= f_hi.  Returns (L0, fr, max relative residual).
    """
    fg = np.asarray(f) / 1e9; Ln = np.asarray(L_app) * 1e9          # GHz, nH (conditioning)
    m = fg <= ((f_hi or np.max(f)) / 1e9)
    A = np.vstack([np.ones(m.sum()), fg[m] ** 2]).T
    coef, *_ = np.linalg.lstsq(A, 1 / Ln[m], rcond=None)
    L0 = 1 / coef[0] * 1e-9
    fr = np.sqrt(-coef[0] / coef[1]) * 1e9 if coef[1] < 0 else np.inf
    resid = np.max(np.abs(1 / (A @ coef) - Ln[m]) / Ln[m])
    return float(L0), float(fr), float(resid)
