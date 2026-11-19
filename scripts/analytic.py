"""Closed-form / semi-analytical inductance references.

All lengths in metres unless stated. These are *independent* of FastHenry and
openEMS and serve as the first verification tier.

References
----------
[Grover]  F. W. Grover, Inductance Calculations, Dover, 1946/2004.
[Ruehli]  A. E. Ruehli, "Inductance calculations in a complex integrated circuit
          environment," IBM J. Res. Dev., 16(5), 1972.  Partial self-inductance
          of a rectangular bar.
[Letellier2018] A. Letellier, M. Dubois, J. P. F. Trovao, H. Maher, "Calculation
          of PCB power loop stray inductance in GaN or high di/dt applications,"
          IEEE TPEL, 2018, DOI 10.1109/TPEL.2018.2826920.  Their Eq. (13) is the
          exact 2D (infinite-length) inductance of two anti-parallel rectangular
          conductors carrying uniform current.  We evaluate the same physics
          through geometric mean distances (GMD), which is mathematically
          equivalent for uniform current density, and check it against the
          numbers printed in their Table 1.
"""
from __future__ import annotations

import numpy as np

MU0 = 4e-7 * np.pi


def bar_partial_self_L(l: float, w: float, t: float) -> float:
    """Partial self-inductance of a straight rectangular bar (Grover/Ruehli approx.).

    L = (mu0 l / 2pi) [ ln(2l/(w+t)) + 1/2 + 0.2235 (w+t)/l ],  valid for l >> w,t.
    """
    return MU0 * l / (2 * np.pi) * (np.log(2 * l / (w + t)) + 0.5 + 0.2235 * (w + t) / l)


def parallel_plate_L(l: float, w: float, e: float) -> float:
    """Classical wide-plate approximation, Letellier Eq. (1): L = mu0 e l / w."""
    return MU0 * e * l / w


def ln_gmd_aligned(w: float, h: float, c: float) -> float:
    """ln(GMD) between two identical w x h rectangles whose centres are offset by c
    along the h direction (c = 0 gives the self-GMD).

    Exact reduction of the 4D mean-log-distance integral: the offset (u, v)
    between two uniformly distributed points has triangular density
    (w-|u|)(h-|v|)/(w h)^2, hence

        ln GMD = (1/(w h)^2) * int_{-w}^{w} int_{-h}^{h} (w-|u|)(h-|v|) ln sqrt(u^2+(v+c)^2) dv du

    The u-integral is done in closed form,
        F(A) = int_0^w (w-u) ln(u^2+A^2) du
             = (w^2/2 - A^2/2) ln(w^2+A^2) + A^2 ln(A^2)/2*... (see code)
    and the remaining 1D v-integral (log singularity only at v = -c when c < h)
    is done with adaptive quadrature, split at the singular point.
    """
    from scipy.integrate import quad

    def F(A: float) -> float:
        # int_0^w (w-u) ln(u^2 + A^2) du, closed form (A >= 0)
        A2 = A * A
        if A2 == 0.0:
            return w * w * (np.log(w) - 1.5)  # = int_0^w (w-u) 2 ln u du
        t1 = w * (w * np.log(w * w + A2) - 2 * w + 2 * A * np.arctan(w / A))  # w*int ln
        t2 = 0.5 * ((w * w + A2) * np.log(w * w + A2) - A2 * np.log(A2) - w * w)  # int u ln
        return t1 - t2

    def g(v: float) -> float:
        return (h - abs(v)) * F(abs(v + c))  # symmetric in u -> factor 2 applied below

    pts = sorted({-h, h, *([-c] if -h < -c < h else [])})
    total = 0.0
    for a, b in zip(pts[:-1], pts[1:]):
        val, _ = quad(g, a, b, limit=400, epsabs=0, epsrel=1e-11)
        total += val
    # factor 2 for u in [-w,0]; factor 1/2 converts ln(r^2) -> ln r
    return float(2 * total * 0.5 / (w * h) ** 2)


def two_bar_loop_L_per_m(w: float, h: float, e: float) -> float:
    """Exact 2D loop inductance per metre of two anti-parallel w x h bars.

    Gap between facing copper surfaces is `e`, so centre spacing is e + h.
    L' = (mu0/pi) [ ln GMD_12 - ln GMD_11 ]   (uniform current, DC limit).
    """
    return MU0 / np.pi * (ln_gmd_aligned(w, h, e + h) - ln_gmd_aligned(w, h, 0.0))


def two_bar_loop_L(l: float, w: float, h: float, e: float) -> float:
    """2D loop inductance times length (end effects neglected), as Letellier Eq. (13)."""
    return l * two_bar_loop_L_per_m(w, h, e)


def ring_L(f_ring: float, c_oss: float) -> float:
    """Loop inductance from ringing frequency: L = 1 / ((2 pi f)^2 C)  (Letellier Eq. 15)."""
    return 1.0 / ((2 * np.pi * f_ring) ** 2 * c_oss)


if __name__ == "__main__":
    # Self-check 1: known GMDs. Square: 0.44705 a; thin strip (h->0): 0.22313 w;
    # two points-like far rectangles: GMD -> centre distance.
    print("square  GMD/a =", np.exp(ln_gmd_aligned(1.0, 1.0, 0.0)), "(expect 0.44705)")
    print("strip   GMD/w =", np.exp(ln_gmd_aligned(1.0, 1e-6, 0.0)), "(expect 0.22313)")
    print("far     GMD/c =", np.exp(ln_gmd_aligned(1e-3, 1e-3, 1.0)), "(expect ~1.0)")
    # Self-check 2: reproduce Letellier Table 1, Eq. (13) column (h = 35 um).
    cases = [("EPC2014 power", 1.7e-3, 9.5e-3, 1e-3, 4.67),
             ("GS61008 power", 7e-3, 15.5e-3, 1e-3, 2.52),
             ("EPC2014 gate", 0.254e-3, 7.2e-3, 1e-3, 7.82),
             ("EPC2014 top-bottom", 1.7e-3, 9.5e-3, 1.6e-3, 6.20)]
    for name, w, l, e, pub in cases:
        L = two_bar_loop_L(l, w, 35e-6, e) * 1e9
        print(f"{name:20s} ours {L:6.3f} nH   Letellier(13) {pub:5.2f} nH   diff {100*(L-pub)/pub:+5.1f}%"
              f"   classical(1) {parallel_plate_L(l, w, e)*1e9:6.2f} nH")
