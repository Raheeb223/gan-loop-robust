"""Parametric EPC-style GaN half-bridge power-loop geometries.

Three families after Reusch & Strydom (APEC 2013 / IEEE TPEL 29(4) 2014) and
EPC WP010: LATERAL, VERTICAL and OPTIMAL power loops.  All three share the same
top-layer "chain"

     [ Cin ]  VIN copper | T-FET | SW copper | SR-FET | GND copper
     x -->

so that only the *return path* differs, exactly the comparison of WP010:

  optimal  : return on inner layer 1 (L2), directly under the chain, via rows
             at the Cin GND pad and after the SR source           (sep. ~ d1)
  vertical : Cin on the bottom layer under the devices, return on L4 (sep. ~ T)
             inner layers absent in the loop region (no shield)
  lateral  : Cin on top, return strip on top layer beside the chain (U-loop);
             L2 is a full "shield" plane tied to GND at ONE location (eddy
             currents only).  Variant `lateral_stitched` ties it at both ends.

ASSUMPTIONS (documented, not taken from EPC drawings, which are not published
with dimensions):
  * EPC2015C die 4.105 x 1.635 mm (datasheet); chain width W = 4.1 mm.
  * Flip-chip die: current crosses the device at bump height (100 um, datasheet
    seating height) through k parallel die "bridges" (one per D/S finger pair).
    Device-internal (package) inductance is NOT included (it belongs to the SPICE
    device model) -- identical treatment for all families.
  * Cin bank = distributed ideal port across a 0.5 mm pad gap (bank of 0402s).
    MLCC body ESL excluded (added in SPICE).
  * Vias: 0.3 mm finished hole -> modelled as 0.25 x 0.25 mm square bar, rows of
    `n_via` across the copper width.
All lengths mm.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass

import numpy as np

from geom import Geometry

FET_L = 1.635      # EPC2015C die length along x (mm)
BUMP = 0.100       # solder bump / seating height (mm)
VIA_W = 0.25       # square-bar equivalent of 0.3 mm via
DIE_H = 0.02       # die metallisation thickness used for the bridge (mm)


@dataclass
class Params:
    family: str = "optimal"         # optimal | vertical | lateral | lateral_stitched
    T: float = 0.787                # board thickness, copper centre L1 to L4 (mm) (31 mil)
    d1: float = 0.102               # dielectric L1-L2 (mm) (4 mil)
    cu: float = 0.070               # outer copper thickness (2 oz, WP010)
    cu_in: float | None = None      # inner copper thickness (None = same as outer)
    W: float = 4.1                  # chain width (mm)
    l_vin: float = 1.2              # VIN copper length before T (mm)
    l_sw: float = 1.0               # SW copper between devices (mm)
    l_gnd: float = 1.0              # GND copper after SR (mm)
    cap_w: float = 1.2              # Cin bank width (lateral) (mm) = VIN copper before the die (2x 0402)
    cap_gap: float = 0.5            # Cin pad gap (port) (mm)
    gap_ret: float = 0.5            # lateral: chain-to-return-strip gap (mm)
    w_ret: float = 2.0              # lateral: return strip width (mm)
    n_via: int = 5                  # vias per row
    n_bridge: int = 5               # die bridges per device
    pitch: float = 0.2              # plane discretisation (mm)
    margin_plane: float = 1.0       # inner plane overhang beyond loop (mm)
    # --- tolerance variables (Phase 4); defaults reproduce the original geometry
    fet_l: float = FET_L            # die / package length along the current (mm) (Phase 6 devices)
    bump: float = BUMP              # solder bump / seating height (mm)
    via_w: float = VIA_W            # via barrel, square-bar equivalent (mm)
    dx_die: float = 0.0             # die placement offset along x (mm)
    # --- gate loops (Phase 3); 'none' reproduces the original power-loop-only geometry
    gate: str = "none"              # none | kelvin | common  (common: LS driver return via to L2)
    l_gate: float = 3.0             # driver-to-gate-pad trace length (mm)
    w_gate: float = 0.3             # gate / kelvin trace width (mm)
    s_gate: float = 0.5             # gate-to-return trace centre spacing (mm)

    def tag(self) -> str:
        ci = "" if self.cu_in is None else f"i{self.cu_in:.4f}"
        return (f"{self.family}_T{self.T:.3f}_d{self.d1:.3f}_cu{self.cu:.3f}{ci}"
                f"_v{self.n_via}_p{self.pitch:.2f}")


def _z(p: Params):
    zt = p.T
    zb = 0.0
    ci = p.cu if p.cu_in is None else p.cu_in
    z2 = p.T - (p.d1 + p.cu / 2 + ci / 2)   # centre of L2 (d1 = dielectric between copper faces)
    return zt, zb, z2


def _ys(W, k, inset=0.25):
    return np.linspace(-W / 2 + inset, W / 2 - inset, k)


def _device(g: Geometry, tag, pl_in, pl_out, x0, p: Params, zt):
    """FET between plane pl_in (left) and pl_out (right), die spanning x0..x0+p.fet_l."""
    x0 = x0 + p.dx_die
    for i, y in enumerate(_ys(p.W, p.n_bridge, inset=p.W / (2 * p.n_bridge))):
        a = g.contact(pl_in, f"{tag}a{i}", x0 + 0.3, y)
        b = g.contact(pl_out, f"{tag}b{i}", x0 + p.fet_l - 0.3, y)
        zd = zt + p.cu / 2 + p.bump
        ua = g.node(f"{tag}ua{i}", x0 + 0.3, y, zd)
        ub = g.node(f"{tag}ub{i}", x0 + p.fet_l - 0.3, y, zd)
        wb = p.W / p.n_bridge * 0.8
        g.bar(a, ua, 0.18 * 1.5, 0.18 * 1.5, wdir=(0, 1, 0))    # solder bumps (bar ~ 0.27 mm)
        g.bar(ub, b, 0.18 * 1.5, 0.18 * 1.5, wdir=(0, 1, 0))
        g.bar(ua, ub, wb, DIE_H)                                # die bridge


def _via_row(g: Geometry, tag, pl_top, pl_bot, x, ycentre, width, n, via_w=VIA_W):
    for i, y in enumerate(_ys(width, n, inset=0.3) + ycentre):
        a = g.contact(pl_top, f"{tag}t{i}", x, y)
        b = g.contact(pl_bot, f"{tag}b{i}", x, y)
        g.bar(a, b, via_w, via_w, wdir=(0, 1, 0))


def _chain(g: Geometry, p: Params, zt, x_start):
    """Top-layer chain from x_start: VIN | T | SW | SR | GND. Returns planes and x marks."""
    W = p.W
    xT = x_start + p.l_vin
    xS = xT + p.fet_l + p.l_sw
    x_end = xS + p.fet_l + p.l_gnd
    mT, mS = xT + p.fet_l / 2, xS + p.fet_l / 2
    vin = g.plane("vin", x_start, -W / 2, mT - 0.1, W / 2, zt, p.cu, p.pitch)
    sw = g.plane("sw", mT + 0.1, -W / 2, mS - 0.1, W / 2, zt, p.cu, p.pitch)
    gnd = g.plane("gnd", mS + 0.1, -W / 2, x_end, W / 2, zt, p.cu, p.pitch)
    _device(g, "T", vin, sw, xT, p, zt)
    _device(g, "S", sw, gnd, xS, p, zt)
    return vin, sw, gnd, dict(xT=xT, xS=xS, x_end=x_end, mT=mT, mS=mS)


def _bundle_port(g: Geometry, plane_p, plane_m, pts_p, pts_m, port_box):
    pp = [g.contact(plane_p, f"pp{i}", x, y) for i, (x, y) in enumerate(pts_p)]
    pm = [g.contact(plane_m, f"pm{i}", x, y) for i, (x, y) in enumerate(pts_m)]
    g.equivs += [pp, pm]
    g.ports.append((pp[0], pm[0]))
    g.port_box = port_box


def build(p: Params) -> Geometry:
    zt, zb, z2 = _z(p)
    W = p.W
    g = Geometry(p.tag(), notes=str(asdict(p)))
    k = p.n_via
    if p.family == "optimal":
        # Cin GND pad | gap(port) | chain ...
        x_gp1 = 0.6
        x_chain = x_gp1 + p.cap_gap
        gpad = g.plane("gpad", 0.0, -W / 2, x_gp1, W / 2, zt, p.cu, p.pitch)
        vin, sw, gnd, m = _chain(g, p, zt, x_chain)
        l2 = g.plane("l2", -p.margin_plane, -W / 2 - p.margin_plane, m["x_end"] + p.margin_plane,
                     W / 2 + p.margin_plane, z2, p.cu if p.cu_in is None else p.cu_in, max(p.pitch, 0.25))
        _via_row(g, "vg", gnd, l2, m["x_end"] - 0.45, 0.0, W, k, via_w=p.via_w)
        _via_row(g, "vc", gpad, l2, 0.3, 0.0, W, k, via_w=p.via_w)
        ys = _ys(W, 6, 0.2)
        _bundle_port(g, vin, gpad, [(x_chain + 0.1, y) for y in ys], [(x_gp1 - 0.1, y) for y in ys],
                     ([x_gp1, -W / 2, zt], [x_chain, W / 2, zt], "x"))
    elif p.family == "vertical":
        vin, sw, gnd, m = _chain(g, p, zt, 0.0)
        # bottom: VIN pad under VIN copper / T, cap gap, GND plane to the end
        xv = 0.45
        xb1 = m["xT"] + 0.3            # cap gap starts under the T device
        bv = g.plane("bvin", 0.0, -W / 2, xb1, W / 2, zb, p.cu, p.pitch)
        bg = g.plane("bgnd", xb1 + p.cap_gap, -W / 2, m["x_end"], W / 2, zb, p.cu, p.pitch)
        _via_row(g, "vv", vin, bv, xv, 0.0, W, k, via_w=p.via_w)
        _via_row(g, "vg", gnd, bg, m["x_end"] - 0.45, 0.0, W, k, via_w=p.via_w)
        ys = _ys(W, 6, 0.2)
        _bundle_port(g, bv, bg, [(xb1 - 0.1, y) for y in ys], [(xb1 + p.cap_gap + 0.1, y) for y in ys],
                     ([xb1, -W / 2, zb], [xb1 + p.cap_gap, W / 2, zb], "x"))
    elif p.family in ("lateral", "lateral_stitched"):
        vin, sw, gnd_c, m = _chain(g, p, zt, 0.0)
        y_r0 = W / 2 + p.gap_ret
        y_r1 = y_r0 + p.w_ret
        # GND "elbow": plane after SR spanning chain + gap + return strip height
        g.planes.remove(gnd_c)
        gnd = g.plane("gnd", m["mS"] + 0.1, -W / 2, m["x_end"], y_r1, zt, p.cu, p.pitch)
        # re-attach SR bridge contacts that pointed at the removed plane
        gnd.contacts += gnd_c.contacts
        ret = g.plane("ret", 0.0, y_r0, m["mS"] - 0.15, y_r1, zt, p.cu, p.pitch)
        # seam between return strip and elbow (planes are separate rectangles)
        for i, y in enumerate(np.linspace(y_r0 + 0.2, y_r1 - 0.2, 5)):
            a = g.contact(ret, f"seamA{i}", m["mS"] - 0.25, y)
            b = g.contact(gnd, f"seamB{i}", m["mS"] + 0.2, y)
            g.bar(a, b, (p.w_ret - 0.4) / 5 * 0.9, p.cu)
        l2 = g.plane("l2", -p.margin_plane, -W / 2 - p.margin_plane, m["x_end"] + p.margin_plane,
                     y_r1 + p.margin_plane, z2, p.cu if p.cu_in is None else p.cu_in, max(p.pitch, 0.25))
        _via_row(g, "vg", gnd, l2, m["x_end"] - 0.45, 0.0, W, k, via_w=p.via_w)
        if p.family == "lateral_stitched":
            _via_row(g, "vr", ret, l2, 0.4, (y_r0 + y_r1) / 2, p.w_ret, max(2, k // 2), via_w=p.via_w)
        xs = np.linspace(0.2, p.cap_w - 0.2, 4)
        _bundle_port(g, vin, ret, [(x, W / 2 - 0.1) for x in xs], [(x, y_r0 + 0.1) for x in xs],
                     ([0.0, W / 2, zt], [p.cap_w, y_r0, zt], "y"))
    else:
        raise ValueError(p.family)
    if p.gate != "none":
        _gates(g, p, zt)
    return g


def _gates(g: Geometry, p: Params, zt):
    """Gate-drive loops of both FETs as extra ports (Phase 3).

    Gate pad (datasheet pad 1) sits just outside the power copper at y = -W/2 - 0.3 mm,
    x = pad-row centre, exactly where kicad_prims.py places it.  A w_gate trace runs
    l_gate in -y to the driver.  Inside the device the gate runner is modelled as thin
    die metal from the gate bump to the source-side bump row (node ub0), so the source
    bumps and source copper are shared with the power loop: that shared path is the
    common-source inductance and appears as the off-diagonal Z entries.
    Return path:  kelvin -> separate trace from the source copper next to the die back to
    the driver (both FETs);  common -> (low-side only) driver return through a via into
    the inner ground plane L2, so gate current shares the power return (no Kelvin).
    Ports appended after the power port: [power, gate_HS, gate_LS].
    """
    W = p.W
    l2 = next((pl for pl in g.planes if pl.name == "l2"), None)
    zd = zt + p.cu / 2 + p.bump
    yg = -W / 2 - 0.30
    ydrv = yg - p.l_gate
    if l2 is not None:          # inner ground plane continues under the driver area
        l2.y0 = min(l2.y0, ydrv - 0.6)
    for tag, src in (("T", "sw"), ("S", "gnd")):
        pl_src = next(pl for pl in g.planes if pl.name == src)
        x0 = g.nodes[f"{tag}a0"][0] - 0.3 - p.dx_die       # nominal pad row
        gx = x0 + p.fet_l / 2
        xs = g.nodes[f"{tag}ub0"][0]
        ys0 = g.nodes[f"{tag}ub0"][1]
        gd = g.node(f"g{tag}d", gx, ydrv, zt)
        gp = g.node(f"g{tag}p", gx, yg, zt)
        gu = g.node(f"g{tag}u", gx, yg, zd)
        gk = g.node(f"g{tag}k", xs, yg, zd)
        g.bar(gd, gp, p.w_gate, p.cu)                                   # gate trace
        g.bar(gp, gu, 0.18 * 1.5, 0.18 * 1.5, wdir=(1, 0, 0))           # gate bump
        g.bar(gu, gk, 0.1, DIE_H)                                       # on-die gate runner
        g.bar(gk, f"{tag}ub0", 0.1, DIE_H)
        xr = gx + p.s_gate
        rd = g.node(f"r{tag}d", xr, ydrv, zt)
        if p.gate == "common" and tag == "S":
            if l2 is None:
                raise ValueError("common gate return needs an inner ground plane")
            rv = g.node(f"r{tag}v", xr, yg - 0.3, zt)
            g.bar(rd, rv, p.w_gate, p.cu)
            c2 = g.contact(l2, f"r{tag}l2", xr, yg - 0.3)
            g.bar(rv, c2, p.via_w, p.via_w, wdir=(0, 1, 0))                # driver-ground via
        else:
            rc = g.contact(pl_src, f"r{tag}c", xr, -W / 2 + 0.05)          # Kelvin pick-off
            g.bar(rd, rc, p.w_gate, p.cu)
        g.ports.append((gd, rd))
