"""Device registry for Phase 6: footprint geometry + datasheet-only ngspice models.

EPC2015C keeps the curve-digitised model of epc_model.py (validated earlier).  The
newer devices are modelled from datasheet TABLE values only (no curve reading), with
a method that is first validated on EPC2015C (see p6_devices.py):

  Coss(V) = C0 / (1 + V/Vj)^m,  three unknowns fixed by three datasheet numbers:
    Coss(V_test) (table), Coss,er = 2/V^2 int_0^V C v dv (energy related),
    Coss,tr = Qoss / V = 1/V int_0^V C dv (time related), all at V = V_test.
  Crss(V): the EPC2015C Crss curve shape, voltage axis scaled by the rating ratio and
    amplitude scaled to the table Crss at V_test (ASSUMPTION: same shape family).
  Cgs = Ciss - Crss at V_test.
  Channel: same smooth law as epc_model.py with transconductance scaled by
    RDS(on),typ ratio and threshold shifted by the VGS(TH),typ difference (ASSUMPTION).

Datasheet values (EPC, transcribed from the PDFs, revision dates below):
  EPC2015C  rev. Apr 2021 : 40 V, RDS(on) 3.2/4 mOhm, VTH 0.8/1.4/2.5 V @ 9 mA,
            Ciss/Crss/Coss @20 V; Coss 710 pF, Coss,tr 940 pF, Qoss 19 nC, RG 0.3 Ohm
  EPC2218   rev. 20 Feb 2026: 100 V, RDS(on) 2.4/3.2 mOhm, VTH 0.8/1.1/2.5 V @ 7 mA;
            @50 V: Ciss 1189, Crss 4.3, Coss 562 pF; Coss,er 740, Coss,tr 925 pF;
            Qoss 46 nC; RG 0.4 Ohm; QG 10.5 nC; die 3.5 x 1.95 mm, 7 solder bars
            (G + 4 S + 3 D) on 0.5 mm pitch, bar 250 um x 1625 um, bar height 120 +- 12 um
  EPC2302   rev. 29 Apr 2026: 100 V, RDS(on) 1.4/1.8 mOhm, VTH 0.8/1.3/2.5 V @ 14 mA;
            @50 V: Ciss 3200, Crss 7, Coss 1000 pF; Coss,er 1300, Coss,tr 1700 pF;
            Qoss 85 nC; RG 0.5 Ohm; QG 23 nC; QFN 3 x 5 x 0.65 mm, 7 pins
            (G + 3 S + 3 D interleaved), lead-frame thickness A3 = 0.25 mm
Footprint model (solver): the die/package spans fet_l along the current and the chain
width W across it; n_bridge = number of drain bars/pins.  For the QFN the vertical
connection height is taken as the lead-frame thickness (ASSUMPTION).
"""
from __future__ import annotations

import numpy as np
from scipy.optimize import least_squares

import epc_model

DEV = {
    "EPC2015C": dict(vds=40, v_test=20, ron=3.2e-3, vth=(0.8, 1.4, 2.5), ciss=None, crss=18, coss=710,
                     coss_er=870, coss_tr=940, qoss=19e-9, rg=0.3,
                     geo=dict(fet_l=1.635, W=4.1, n_bridge=5, bump=0.100)),
    "EPC2218": dict(vds=100, v_test=50, ron=2.4e-3, vth=(0.8, 1.1, 2.5), ciss=1189, crss=4.3, coss=562,
                    coss_er=740, coss_tr=925, qoss=46e-9, rg=0.4,
                    geo=dict(fet_l=1.95, W=3.5, n_bridge=3, bump=0.120)),
    "EPC2302": dict(vds=100, v_test=50, ron=1.4e-3, vth=(0.8, 1.3, 2.5), ciss=3200, crss=7, coss=1000,
                    coss_er=1300, coss_tr=1700, qoss=85e-9, rg=0.5,
                    geo=dict(fet_l=3.0, W=5.0, n_bridge=3, bump=0.250)),
}


def coss_eq(V, C0, Vj, m):
    return C0 / (1 + V / Vj) ** m


def fit_coss(coss, coss_er, coss_tr, v_test):
    """Fit C0, Vj, m to the three datasheet numbers (pF, V).  Returns (C0, Vj, m, residuals)."""
    v = np.linspace(0, v_test, 2001)

    def res(x):
        C0, Vj, m = np.exp(x)
        c = coss_eq(v, C0, Vj, m)
        er = 2 / v_test ** 2 * np.trapezoid(c * v, v)
        tr = np.trapezoid(c, v) / v_test
        return [(coss_eq(v_test, C0, Vj, m) - coss) / coss, (er - coss_er) / coss_er, (tr - coss_tr) / coss_tr]

    best = None
    for vj0 in (0.5, 2, 5, 20):
        for m0 in (0.3, 0.6, 1.0):
            r = least_squares(res, np.log([coss * 3, vj0, m0]), xtol=1e-12, ftol=1e-12)
            if best is None or r.cost < best.cost:
                best = r
    C0, Vj, m = np.exp(best.x)
    return C0, Vj, m, np.array(res(best.x))


def curves(name):
    """C(V) tables (pF) on the epc_model voltage grid scaled to the device rating."""
    d = DEV[name]
    k = d["vds"] / 40.0
    V = epc_model.V_PTS * k
    C0, Vj, m, _ = fit_coss(d["coss"], d["coss_er"], d["coss_tr"], d["v_test"])
    coss = coss_eq(V, C0, Vj, m)
    crss_shape = np.interp(V / k, epc_model.V_PTS, epc_model.CRSS_PF)
    crss = crss_shape * d["crss"] / np.interp(d["v_test"] / k, epc_model.V_PTS, epc_model.CRSS_PF)
    cgs = d["ciss"] - d["crss"]
    return V, coss, crss, cgs


def subckt(name: str, coss_scale: float = 1.0) -> tuple[str, str]:
    """(subckt name, .subckt text).  Pins D G S."""
    if name == "EPC2015C":
        return "EPC2015C_DS", epc_model.subckt("EPC2015C_DS", coss_scale=coss_scale)
    d = DEV[name]
    V, coss, crss, cgs = curves(name)
    cds = coss - crss
    # fine knots at low voltage, where the junction-like Coss(V) is steepest
    vf = np.unique(np.concatenate([np.arange(0, 10, 0.2), np.arange(10, V[-1] + 1e-9, 1.25)]))
    C0, Vj, m, _ = fit_coss(d["coss"], d["coss_er"], d["coss_tr"], d["v_test"])
    crss_f = np.interp(vf, V, crss)
    c_ds = (coss_eq(vf, C0, Vj, m) - crss_f) * coss_scale / 1000.0
    c_gd = crss_f / 1000.0
    gm = 150.0 * 3.2e-3 / d["ron"]
    vt0 = 2.25 + (d["vth"][1] - 1.4)
    ich = "( {I}*tanh(V(d,s)/sqrt(({I}*{R})*({I}*{R})+4e-4)) )"
    isat = lambda v: f"({gm:.4g}*0.15*ln(1+exp(({v}-{vt0:.4g})/0.15)))"
    ron = lambda v: f"({d['ron']:.4g}*sqrt(2.75/max({v}-{vt0:.4g},0.05)))"
    fwd = ich.format(I=isat("V(gi,s)"), R=ron("V(gi,s)"))
    rev = ich.format(I=isat("V(gi,d)"), R=ron("V(gi,d)"))

    def pwl(expr, c):
        xs = np.concatenate([[-20.0], vf]); ys = np.concatenate([[c[0]], c])
        return f"pwl({expr}, " + ", ".join(f"{a:g}, {b:.6e}" for a, b in zip(xs, ys)) + ")"

    sn = f"{name}_DS"
    return sn, f"""* {sn}: datasheet-table model (scripts/devices.py)
.subckt {sn} d g s
RG g gi {d['rg']}
Bch d s I = V(d,s) >= 0 ? {fwd} : {rev}
Eds xds 0 d s 1
Rxds xds xds2 1e-4
Cxds xds2 yds 1n
Vxds yds 0 0
Bcds d s I = {pwl('V(d,s)', c_ds)} * i(Vxds)
Egd xgd 0 d gi 1
Rxgd xgd xgd2 1e-4
Cxgd xgd2 ygd 1n
Vxgd ygd 0 0
Bcgd d gi I = {pwl('V(d,gi)', c_gd)} * i(Vxgd)
Cgs gi s {cgs}p
Rleak d s 1e6
.ends {sn}
"""
