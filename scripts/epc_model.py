"""Datasheet-derived behavioural model of the EPC2015C eGaN FET for ngspice.

Every number below is read from the EPC2015C datasheet (rev. April 2021,
data/published/sources/EPC2015C_datasheet.pdf).  Curve values are visual reads
of the printed figures (estimated reading accuracy +-3 %).

Channel (Fig. 1, Fig. 2, Fig. 3):
    I_sat(Vgs) = GM * s * ln(1 + exp((Vgs - VT0)/s))         (smooth turn-on)
    I_ch       = I_sat * tanh(Vds / (I_sat * R_on(Vgs)))      (linear -> saturation)
  GM = 150 A/V, VT0 = 2.25 V, s = 0.15 V fit to Fig. 2 at Vds = 3 V:
    Id(2.55 V) ~ 50 A, Id(3.6 V) ~ 200 A.
  R_on(Vgs) = 3.2 mOhm * sqrt(2.75 V / Vov) -> 3.2 mOhm at Vgs = 5 V (typ. RDS(on)).
  Third quadrant: the lateral HEMT is symmetric; for Vds < 0 the roles of
  drain and source swap (Vgd controls), which yields the ~1.7-2 V reverse drop
  of Fig. 7 at moderate current (no body diode, Qrr = 0).

Capacitances (Fig. 5a, VGS = 0, linear scale), pF:
  Vds   0     5     10    15    20    25    30    35    40
  Coss  1140  1050  930   820   710   590   510   470   440   (Coss(20V)=710 typ: table)
  Crss  160   110   40    25    18    16    14    12    10    (Crss(20V)=18: table)
  Ciss ~ 1100 -> 980 (nearly flat) => Cgs = Ciss - Crss ~ 960 pF
Implemented as charge tables Q(V) = int C dV (charge-conserving) with
behavioural current sources I = dQ/dt.
Gate resistance RG = 0.3 Ohm (table).
"""
from __future__ import annotations

import numpy as np

V_PTS = np.array([0, 5, 10, 15, 20, 25, 30, 35, 40, 48], float)
COSS_PF = np.array([1140, 1050, 930, 820, 710, 590, 510, 470, 440, 420], float)
CRSS_PF = np.array([160, 110, 40, 25, 18, 16, 14, 12, 10, 9], float)
CGS_PF = 960.0
RG = 0.3


def charge_table(c_pf, v=V_PTS, scale=1.0):
    """Cumulative charge Q(V) (coulombs) for a piecewise-linear C(V), extended
    linearly for V < 0 with C(0)."""
    # resample C(V) on a 0.5 V grid so that the piecewise-linear Q(V) gives a
    # piecewise-constant C that tracks the curve (5 V knots would average it)
    vf = np.arange(v[0], v[-1] + 1e-9, 0.5)
    c = np.interp(vf, v, c_pf) * 1e-12 * scale
    v = vf
    q = np.concatenate([[0.0], np.cumsum(0.5 * (c[1:] + c[:-1]) * np.diff(v))])
    vv = np.concatenate([[-20.0], v])
    qq = np.concatenate([[-20.0 * c[0]], q])
    return vv, qq


def _pwl(vexpr, vv, qq):
    pts = ", ".join(f"{a:g}, {b:.6e}" for a, b in zip(vv, qq))
    return f"pwl({vexpr}, {pts})"


def subckt(name: str = "EPC2015C_DS", coss_scale: float = 1.0) -> str:
    """Return the .subckt text. Pins: D G S.  `coss_scale` scales Cds (for the
    datasheet typ->max sensitivity study)."""
    cds_pf = COSS_PF - CRSS_PF
    vv, q_ds = charge_table(cds_pf, scale=coss_scale)
    _, q_gd = charge_table(CRSS_PF)
    ich = ("( {I}*tanh(V(d,s)/sqrt(({I}*{R})*({I}*{R})+4e-4)) )")  # knee floor 20 mV
    isat_gs = "(150*0.15*ln(1+exp((V(gi,s)-2.25)/0.15)))"
    isat_gd = "(150*0.15*ln(1+exp((V(gi,d)-2.25)/0.15)))"
    ron_gs = "(3.2e-3*sqrt(2.75/max(V(gi,s)-2.25,0.05)))"
    ron_gd = "(3.2e-3*sqrt(2.75/max(V(gi,d)-2.25,0.05)))"
    fwd = ich.format(I=isat_gs, R=ron_gs)
    rev = ich.format(I=isat_gd, R=ron_gd)
    vf = np.arange(V_PTS[0], V_PTS[-1] + 1e-9, 0.5)
    c_ds = np.interp(vf, V_PTS, cds_pf) * coss_scale / 1000.0      # ratio to Cref = 1 nF
    c_gd = np.interp(vf, V_PTS, CRSS_PF) / 1000.0
    return f"""* {name}: datasheet-derived behavioural EPC2015C model (see scripts/epc_model.py)
.subckt {name} d g s
RG g gi {RG}
* channel: forward (Vds>=0) controlled by Vgs, reverse controlled by Vgd
Bch d s I = V(d,s) >= 0 ? {fwd} : {rev}
* nonlinear capacitances, capacitance-multiplier form:
*   I = (C(V)/Cref) * Cref dV/dt, with Cref dV/dt sensed in a linear 1 nF cap.
*   (ngspice B-source ddt() of a nonlinear charge was numerically unstable here.)
Eds xds 0 d s 1
Rxds xds xds2 1e-4
Cxds xds2 yds 1n
Vxds yds 0 0
Bcds d s I = {_pwl('V(d,s)', np.concatenate([[-20.0], vf]), np.concatenate([[c_ds[0]], c_ds]))} * i(Vxds)
Egd xgd 0 d gi 1
Rxgd xgd xgd2 1e-4
Cxgd xgd2 ygd 1n
Vxgd ygd 0 0
Bcgd d gi I = {_pwl('V(d,gi)', np.concatenate([[-20.0], vf]), np.concatenate([[c_gd[0]], c_gd]))} * i(Vxgd)
Cgs gi s {CGS_PF}p
* tiny leakage for DC convergence
Rleak d s 1e6
.ends {name}
"""


if __name__ == "__main__":
    print(subckt())
