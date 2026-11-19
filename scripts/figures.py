"""Publication figures (IEEE single column 3.5 in; PDF + 300 dpi PNG).

Colour: validated categorical palette (dataviz skill reference palette, slots
1-4 in fixed order; CVD check PASS).  Aqua/yellow are < 3:1 against the
surface, so every series also carries a distinct marker and a direct label.
Each figure is generated only if its input table exists.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
TAB = ROOT / "results" / "tables"
FIG = ROOT / "results" / "figures"
FIG.mkdir(parents=True, exist_ok=True)

COL = {"optimal": "#2a78d6", "lateral": "#eb6834", "vertical": "#1baf7a", "lateral_stitched": "#eda100"}
MRK = {"optimal": "o", "lateral": "s", "vertical": "^", "lateral_stitched": "D"}
LAB = {"optimal": "Optimal", "lateral": "Lateral", "vertical": "Vertical", "lateral_stitched": "Lateral, stitched"}
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"

plt.rcParams.update({
    "font.family": "serif", "font.size": 8, "axes.labelsize": 8, "legend.fontsize": 7,
    "xtick.labelsize": 7, "ytick.labelsize": 7, "axes.edgecolor": INK2, "axes.labelcolor": INK,
    "xtick.color": INK2, "ytick.color": INK2, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.5,
    "axes.spines.top": False, "axes.spines.right": False, "lines.linewidth": 1.5, "lines.markersize": 4,
    "savefig.dpi": 300, "savefig.bbox": "tight", "legend.frameon": False,
})
W1 = 3.5


def save(fig, name):
    fig.savefig(FIG / f"{name}.pdf")
    fig.savefig(FIG / f"{name}.png")
    plt.close(fig)
    print("figure", name)


def fig_canonical():
    f = TAB / "canonical_fh.csv"
    if not f.exists():
        return
    d = pd.read_csv(f)
    v3 = d[(d.tier == "V3") & (d.fw == 0.0625) & (d.f_Hz == 1e6) & (~d.case.str.contains("Fig6"))]
    fig, ax = plt.subplots(figsize=(W1, 2.4))
    names = {"L18_EPC2014_inner": "EPC2014 inner", "L18_GS61008_inner": "GS61008 inner",
             "L18_EPC2014_gate": "EPC2014 gate", "L18_EPC2014_topbot": "EPC2014 top-bot."}
    x = np.arange(len(v3))
    ax.axhspan(-5, 5, color=GRID, alpha=0.6, lw=0)
    ax.plot(x, v3.err_pct, "o", color=COL["optimal"], label="FastHenry (this work)")
    ax.plot(x, 100 * (v3.eq13_nH - v3.ref_nH) / v3.ref_nH, "s", color=COL["lateral"], label="Letellier Eq. (13)")
    ax.plot(x, 100 * (v3.exact2d_nH - v3.ref_nH) / v3.ref_nH, "^", color=COL["vertical"], label="Exact 2-D (GMD)")
    ax.axhline(0, color=INK2, lw=0.8)
    ax.set_xticks(x, [names[c] for c in v3.case], rotation=15)
    ax.set_ylabel("Deviation from MagNet FEA (%)")
    ax.legend(loc="lower left", ncol=1)
    save(fig, "fig_canonical_vs_fea")


def fig_convergence():
    f = TAB / "hf_ladder.csv"
    if not f.exists():
        return
    d = pd.read_csv(f)
    fig, ax = plt.subplots(figsize=(W1, 2.3))
    for fam in ["optimal", "lateral_stitched", "lateral", "vertical"]:
        s = d[d.family == fam].sort_values("nh")
        ref = s.L_1e08_nH.iloc[-1] if "L_1e08_nH" in s else s["L_1e+08_nH"].iloc[-1]
        y = 100 * (s["L_1e+08_nH"] / ref - 1)
        ax.plot(s.nh, y, marker=MRK[fam], color=COL[fam], label=LAB[fam])
    ax.set_xlabel("Filaments across copper thickness, $n_h$")
    ax.set_ylabel("$L$(100 MHz) rel. to $n_h$=7 (%)")
    ax.set_xticks([1, 3, 5, 7])
    ax.legend(title="Layout (31 mil / 4 mil, 2 oz)", title_fontsize=7)
    save(fig, "fig_hf_convergence")


def fig_epc_trends():
    f = TAB / "epc_sweep_hf.csv"
    b = ROOT / "data" / "published" / "epc_fig7b_digitised.csv"
    if not f.exists():
        return
    d = pd.read_csv(f)
    d = d[d.err.isna() | (d.err == "")]
    d["T_mil"] = (d["T"] / 0.0254).round().astype(int)
    d["d1_mil"] = (d["d1"] / 0.0254).round().astype(int)
    bm = pd.read_csv(b).rename(columns={"family": "layout_type", "L_sim_nH": "value"})
    fig, axs = plt.subplots(1, 3, figsize=(7.16, 2.3), sharey=True)
    for ax, d1 in zip(axs, [4, 12, 26]):
        for fam in ["optimal", "lateral", "vertical", "lateral_stitched"]:
            s = d[(d.family == fam) & (d.d1_mil == d1)].sort_values("T_mil")
            if s.empty:
                continue
            ax.plot(s.T_mil, s["L_1e+08_nH"], marker=MRK[fam], color=COL[fam], label=LAB[fam])
            e = bm[(bm.layout_type == fam) & (bm.inner_layer_mil.astype(float) == d1)]
            if len(e):
                ax.plot(e.board_thickness_mil.astype(float), e.value.astype(float), ls="none", marker=MRK[fam],
                        mfc="white", mec=COL[fam], ms=6, mew=1.2)
        ax.set_title(f"Inner-layer distance {d1} mil", fontsize=8)
        ax.set_xlabel("Board thickness (mil)")
        ax.set_xlim(28, 65)
    axs[0].set_ylabel("Loop inductance at 100 MHz (nH)")
    h, l = axs[0].get_legend_handles_labels()
    from matplotlib.lines import Line2D
    h.append(Line2D([], [], ls="none", marker="o", mfc="white", mec=INK2, ms=6))
    l.append("EPC simulation [2] (digitised, open)")
    fig.legend(h, l, loc="upper center", ncol=5, bbox_to_anchor=(0.5, 1.08))
    save(fig, "fig_epc_trends")


def fig_ring():
    f = TAB / "ring_bias.csv"
    if not f.exists():
        return
    d = pd.read_csv(f)
    fig, ax = plt.subplots(figsize=(W1, 2.3))
    for c, name, k in (("Coss_VIN", "$C_{oss}(V_{IN})$ typ.", "optimal"), ("Coss_tr", "$C_{oss,tr}$", "lateral"),
                       ("Coss_er", "$C_{oss,er}$", "vertical")):
        ax.plot(d.L_true_nH, d[f"err_{c}_pct"], marker=MRK[k], color=COL[k], label=name)
    ax.axhline(0, color=INK2, lw=0.8)
    ax.set_xlabel("True loop inductance (nH)")
    ax.set_ylabel("Ringing-method error (%)")
    ax.legend(title="Capacitance assumed", title_fontsize=7, loc="center right", bbox_to_anchor=(1.0, 0.42))
    save(fig, "fig_ring_bias")


def fig_ring_mc():
    f = TAB / "ring_mc_letellier_samples.npy"
    if not f.exists():
        return
    s = np.load(f)
    fig, ax = plt.subplots(figsize=(W1, 2.0))
    ax.hist(s, bins=60, color=COL["optimal"], edgecolor="#fcfcfb", linewidth=0.3)
    top = ax.get_ylim()[1]
    ax.axvline(0, color=INK2, lw=0.8)
    for v, lab, ha, dx, yf in ((0.38, "reported 0.38 nH", "right", -3, 0.95), (0.53, "Eq. (14): 0.53 nH", "left", 3, 0.80)):
        ax.axvline(v, color=INK, lw=0.8, ls="--")
        ax.annotate(lab, (v, top * yf), xytext=(dx, 0), textcoords="offset points", fontsize=6.5, ha=ha,
                    color=INK, bbox=dict(fc="#fcfcfb", ec="none", pad=0.5))
    neg = 100 * np.mean(s < 0)
    ax.annotate(f"{neg:.1f} % of samples < 0 nH", (0.0, top * 0.30), xytext=(-4, 0), textcoords="offset points",
                fontsize=6.5, color=INK, ha="right", bbox=dict(fc="#fcfcfb", ec="none", pad=0.8))
    ax.set_xlabel("Inferred PCB loop inductance (nH)")
    ax.set_ylabel("MC samples")
    save(fig, "fig_ring_mc")


def fig_overshoot():
    f = TAB / "circuit_validation.csv"
    if not f.exists():
        return
    d = pd.read_csv(f)
    fig, axs = plt.subplots(1, 2, figsize=(7.16, 2.4))
    for fam in ["optimal", "lateral", "vertical"]:
        s = d[d.family == fam]
        axs[0].plot(s.os_meas, s.os_nominal, ls="none", marker=MRK[fam], mfc="white", mec=COL[fam], mew=1.0)
        axs[0].plot(s.os_meas, s.os_pred, ls="none", marker=MRK[fam], color=COL[fam], label=LAB[fam])
        axs[1].plot(s.dv_meas, s.dv_pred, ls="none", marker=MRK[fam], color=COL[fam], label=LAB[fam])
    axs[0].annotate("hollow: uncalibrated\ndatasheet-only model", (100, 128), fontsize=6.5, color=INK2, bbox=dict(fc="#fcfcfb", ec="none", pad=0.8))
    axs[0].annotate("filled: calibrated on optimal;\nlateral, vertical held out", (5, 112), fontsize=6.5, color=INK2)
    for ax, lim, lab in ((axs[0], (0, 160), "overshoot (%)"), (axs[1], (0, 14), "d$v$/d$t$ (V/ns)")):
        ax.plot(lim, lim, color=INK2, lw=0.8)
        ax.set_xlim(lim); ax.set_ylim(lim)
        ax.set_xlabel(f"Measured {lab} [EPC]")
        ax.set_ylabel(f"Predicted {lab}")
    axs[0].legend(loc="lower right")
    save(fig, "fig_circuit_validation")


def fig_designmap():
    f = TAB / "design_map.csv"
    if not f.exists():
        return
    d = pd.read_csv(f, dtype={"stackup": str})
    order = d.drop_duplicates("stackup").sort_values("d1_mm")
    xpos = {s: i for i, s in enumerate(order.stackup)}
    fig, ax = plt.subplots(figsize=(W1, 2.6))
    fams = ["optimal", "lateral_stitched", "lateral", "vertical"]
    for k, fam in enumerate(fams):
        s = d[d.family == fam]
        x = s.stackup.map(xpos) + (k - 1.5) * 0.12
        ax.errorbar(x, s.L_nom_nH, yerr=[s.L_nom_nH - s.L_p2_5_nH, s.L_p97_5_nH - s.L_nom_nH],
                    fmt=MRK[fam], color=COL[fam], label=LAB[fam], capsize=2, lw=1)
    ax.set_xticks(range(len(order)),
                  [f"{s}\n{1000 * t:.0f} µm" for s, t in zip(order.stackup, order.d1_mm)])
    ax.set_xlabel("JLCPCB 1.6 mm 4-layer stackup (L1–L2 prepreg)")
    ax.set_ylabel("Loop inductance, 100 MHz (nH)")
    ax.set_ylim(0, 2.05)
    ax.legend(ncol=2, loc="center left", bbox_to_anchor=(0.0, 0.57))
    ax.annotate("bars: 95 % range, ±10 % thickness tolerance", (0.0, 1.0), xycoords="axes fraction",
                fontsize=6.5, color=INK2, va="top")
    save(fig, "fig_design_map")


def fig_workflow():
    from matplotlib.patches import FancyBboxPatch
    fig, ax = plt.subplots(figsize=(7.16, 2.2))
    ax.set_xlim(0, 100); ax.set_ylim(0, 30); ax.axis("off")
    def box(x, y, w, h, t, c=GRID):
        ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.3,rounding_size=1.2", fc=c, ec=INK2, lw=0.8))
        ax.text(x + w / 2, y + h / 2, t, ha="center", va="center", fontsize=6.8, color=INK)
    def arr(x0, y0, x1, y1):
        ax.annotate("", (x1, y1), (x0, y0), arrowprops=dict(arrowstyle="-|>", color=INK2, lw=0.8))
    box(0.5, 10.5, 16.5, 9, "Parametric\nlayout (layouts.py)\nfamily, stackup,\ndevices, vias")
    box(24, 21, 17, 7, "KiCad 10 board\nDRC, Gerbers, drill")
    box(24, 11.5, 17, 7, "FastHenry (PEEC)\nL(f), R(f)")
    box(24, 2, 17, 7, "openEMS (FDTD)\nZ(f), cross-check")
    box(49, 11.5, 17, 7, "Loop L at f_ring\n+ ESL, package")
    box(49, 21, 17, 7, "EPC2015C model\nfrom datasheet")
    box(74, 11.5, 17, 7, "ngspice bench\novershoot, dv/dt,\nf_ring")
    box(74, 1.5, 24, 7.5, "Verified against: exact 2-D,\npublished FEA, EPC sim. + meas.")
    arr(16, 17, 24, 24.5); arr(16, 15, 24, 15); arr(16, 13, 24, 5.5)
    arr(41, 15, 49, 15); arr(41, 5.5, 49, 13)
    arr(66, 15, 74, 15); arr(66, 24.5, 76, 18.5)
    arr(82.5, 11.5, 82.5, 9)
    save(fig, "fig_workflow")


if __name__ == "__main__":
    for fn in (fig_workflow,fig_canonical, fig_convergence, fig_epc_trends, fig_ring, fig_ring_mc, fig_overshoot, fig_designmap):
        try:
            fn()
        except Exception as ex:  # report and continue
            print("FAILED", fn.__name__, ex)
