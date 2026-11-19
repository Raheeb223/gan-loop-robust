"""Manuscript versions of the Phase 3-5 figures, redrawn from the saved results in a uniform
IEEE style (Times, embedded TrueType fonts, colour-blind-safe colours that also differ in grey level).

Inputs : results/phase3/spice_coupled.csv, results/phase4/{campaign,sobol}.csv + gp_model.pkl,
         results/phase5/{pareto,picks,references,verify}.csv
Outputs: results/figures/fig_p3_false_turn_on, fig_p4_parity, fig_p4_sobol, fig_p4_uq, fig_p5_pareto (.pdf/.png)
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
import figures as F  # noqa: E402  (shared save())

ROOT = Path(__file__).resolve().parents[1]
plt = F.plt
W = 3.5  # IEEE column width (in)
STYLE = {
    "font.family": "serif", "font.serif": ["Times New Roman", "Times", "STIXGeneral"], "mathtext.fontset": "stix",
    "font.size": 8, "axes.labelsize": 8, "xtick.labelsize": 7.5, "ytick.labelsize": 7.5, "legend.fontsize": 7,
    "axes.edgecolor": "black", "axes.labelcolor": "black", "xtick.color": "black", "ytick.color": "black",
    "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6, "xtick.direction": "in",
    "ytick.direction": "in", "axes.grid": True, "grid.color": "0.85", "grid.linewidth": 0.4,
    "axes.spines.top": True, "axes.spines.right": True, "legend.frameon": False, "hatch.linewidth": 0.5,
    "pdf.fonttype": 42, "ps.fonttype": 42, "savefig.dpi": 300, "savefig.bbox": "tight", "savefig.pad_inches": 0.02,
}
BLUE, RED, ORANGE, GREEN = "#1f77b4", "#d62728", "#ff7f0e", "#2ca02c"
KEL = dict(fc=BLUE, ec="black", lw=0.4)      # Kelvin return
COM = dict(fc=RED, ec="black", lw=0.4)       # common return
P3313 = dict(fc=BLUE, ec="black", lw=0.4)
P7628 = dict(fc=ORANGE, ec="black", lw=0.4)


def false_turn_on():
    S = pd.read_csv(ROOT / "results" / "phase3" / "spice_coupled.csv")
    S = S[(S.coupling == "coupled") & S.case.str.endswith("_sweep")]
    fams = [("optimal", "Optimal"), ("lateral", "Lateral"), ("lateral_stitched", "Stitched"), ("vertical", "Vertical")]
    fig, axs = plt.subplots(1, 2, figsize=(W, 1.75), sharey=True)
    for ax, (drive, title) in zip(axs, (("calibrated_20A", "Calibrated, 20 A"), ("strong_1ohm_40A", "1 Ω, 40 A"))):
        d = S[S.drive == drive].set_index("case").vgs_ls_peak
        for i, (f, _) in enumerate(fams):
            if f"{f}_kelvin_sweep" in d:
                ax.bar(i - 0.19, d[f"{f}_kelvin_sweep"], 0.36, **KEL, label="Kelvin" if i == 0 else None)
            if f"{f}_common_sweep" in d:
                ax.bar(i + 0.19, d[f"{f}_common_sweep"], 0.36, **COM, label="common" if i == 0 else None)
        ax.axhline(0.8, color="black", ls="--", lw=0.8)
        ax.set_xticks(range(len(fams)), [n for _, n in fams], fontsize=7)
        ax.set_title(title, fontsize=8, pad=3)
        ax.grid(axis="x", visible=False)
        ax.set_ylim(0, 1.45)
    axs[1].text(3.5, 0.83, r"$V_\mathrm{GS(th),min}$", ha="right", va="bottom", fontsize=7)
    axs[0].set_ylabel("Peak off-state $V_\\mathrm{GS}$ (V)")
    axs[1].legend(loc="upper right", ncol=2, handlelength=1.2, columnspacing=0.8, borderaxespad=0.2)
    fig.tight_layout(pad=0.2, w_pad=0.6)
    F.save(fig, "fig_p3_false_turn_on")


def parity():
    import p5_opt
    from p4_campaign import VARS
    sur = p5_opt.load_sur()
    C = pd.read_csv(ROOT / "results" / "phase4" / "campaign.csv")
    te = C[C.err.isna() & (C.set == "test")]
    lo, mid, hi = sur.interval(te[VARS].values, 0.05)
    fig, ax = plt.subplots(figsize=(W, 2.15))
    ax.errorbar(te.L_nH, mid, yerr=[mid - lo, hi - mid], fmt="o", ms=2.8, mfc=BLUE, mec=BLUE, mew=0.6,
                ecolor="#7fb3d9", elinewidth=0.6, capsize=0, label="Test layout, 95 % conformal interval")
    a = [te.L_nH.min() * 0.95, te.L_nH.max() * 1.05]
    ax.plot(a, a, "k-", lw=0.6, label="Ideal")
    ax.set_xlim(a); ax.set_ylim(a)
    ax.set_xlabel("FastHenry $L$ (nH)")
    ax.set_ylabel("Surrogate $L$ (nH)")
    ax.legend(loc="upper left")
    fig.tight_layout(pad=0.2)
    F.save(fig, "fig_p4_parity")


def sobol():
    S = pd.read_csv(ROOT / "results" / "phase4" / "sobol.csv")
    lab = {"t_d1": "Prepreg\n$d_1$", "t_cu": "Outer\nCu", "t_cuin": "Inner\nCu", "hole": "Via\nhole",
           "bump": "Bump\nheight", "dx_die": "Die\noffset"}
    fig, ax = plt.subplots(figsize=(W, 1.85))
    for k, (st, sty, name) in enumerate((("3313 (nominal)", P3313, "3313 prepreg (nominal)"),
                                          ("7628 (default)", P7628, "7628 prepreg (default)"))):
        g = S[S.stackup == st]
        x = np.arange(len(g)) + (k - 0.5) * 0.38
        ax.bar(x, g.ST, 0.36, yerr=g.ST_conf, **sty, label=name,
               error_kw=dict(ecolor="black", elinewidth=0.6, capsize=1.5, capthick=0.6))
    g = S[S.stackup == "3313 (nominal)"]
    ax.set_xticks(range(len(g)), [lab[v] for v in g.variable])
    ax.set_ylabel("Total Sobol index $S_T$")
    ax.set_ylim(0, 1.0)
    ax.grid(axis="x", visible=False)
    ax.legend(loc="upper right")
    fig.tight_layout(pad=0.2)
    F.save(fig, "fig_p4_sobol")


def pareto():
    P = pd.read_csv(ROOT / "results" / "phase5" / "pareto.csv")
    picks = pd.read_csv(ROOT / "results" / "phase5" / "picks.csv")
    V = pd.read_csv(ROOT / "results" / "phase5" / "verify.csv").set_index("pick").loc[picks.pick]
    refs = pd.read_csv(ROOT / "results" / "phase5" / "references.csv")
    fig, ax = plt.subplots(figsize=(W, 2.2))
    ax.plot(P.copper_area_mm2, P.L95_robust_nH, "o", ms=2.2, color=BLUE, mec="none", label="Pareto front (1080 prepreg)")
    ax.plot(picks.copper_area_mm2, picks.L95_robust_nH, "*", ms=8, color="black", label="Selected designs")
    ax.plot(picks.copper_area_mm2, V.fh_p95_nH.values, "x", ms=6, mew=1.0, color=RED,
            label="FastHenry 95th percentile")
    mk = {"nominal_3313": ("s", "Nominal, 3313 prepreg", ORANGE), "nominal_7628": ("D", "Nominal, 7628 prepreg", GREEN)}
    for _, r in refs.iterrows():
        m, lab, c = mk.get(r["name"], ("^", r["name"], "black"))
        ax.plot(r.copper_area_mm2, r.L95_robust_nH, m, ms=5, mfc="white", mec=c, mew=1.0, label=lab)
    ax.set_xlabel("Copper area next to the devices (mm$^2$)")
    ax.set_ylabel("Robust $L_{95}$ at 100 MHz (nH)")
    ax.legend(loc="upper right", fontsize=6.5, handletextpad=0.4, labelspacing=0.3)
    fig.tight_layout(pad=0.2)
    F.save(fig, "fig_p5_pareto")


def uq_combined():
    """(a) surrogate parity on the test set, (b) total Sobol indices; one column."""
    import p5_opt
    from p4_campaign import VARS
    sur = p5_opt.load_sur()
    C = pd.read_csv(ROOT / "results" / "phase4" / "campaign.csv")
    te = C[C.err.isna() & (C.set == "test")]
    lo, mid, hi = sur.interval(te[VARS].values, 0.05)
    S = pd.read_csv(ROOT / "results" / "phase4" / "sobol.csv")
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(W, 1.9), gridspec_kw=dict(width_ratios=[1, 1.15]))
    a1.errorbar(te.L_nH, mid, yerr=[mid - lo, hi - mid], fmt="o", ms=2.2, mfc=BLUE, mec=BLUE, mew=0.5,
                ecolor="#7fb3d9", elinewidth=0.5, capsize=0)
    lim = [0.27, 0.62]
    a1.plot(lim, lim, "k-", lw=0.5)
    a1.set_xlim(lim); a1.set_ylim(lim)
    a1.set_xticks([0.3, 0.4, 0.5, 0.6]); a1.set_yticks([0.3, 0.4, 0.5, 0.6])
    a1.set_xlabel("FastHenry $L$ (nH)"); a1.set_ylabel("Surrogate $L$ (nH)", labelpad=1)
    a1.set_title("(a)", fontsize=8, pad=3)
    lab = {"t_d1": "$d_1$", "t_cu": r"Cu$_\mathrm{o}$", "t_cuin": r"Cu$_\mathrm{i}$", "hole": "Hole",
           "bump": "Bump", "dx_die": "Die"}
    for k, (st, sty, name) in enumerate((("3313 (nominal)", P3313, "3313"), ("7628 (default)", P7628, "7628"))):
        g = S[S.stackup == st]
        x = np.arange(len(g)) + (k - 0.5) * 0.4
        a2.bar(x, g.ST, 0.38, yerr=g.ST_conf, **sty, label=name,
               error_kw=dict(ecolor="black", elinewidth=0.5, capsize=1.2, capthick=0.5))
    g = S[S.stackup == "3313 (nominal)"]
    a2.set_xticks(range(len(g)), [lab[v] for v in g.variable], fontsize=7)
    a2.set_ylabel("Total Sobol index $S_T$", labelpad=1)
    a2.set_ylim(0, 1.0)
    a2.grid(axis="x", visible=False)
    a2.legend(loc="upper right", handlelength=1.2)
    a2.set_title("(b)", fontsize=8, pad=3)
    fig.tight_layout(pad=0.2, w_pad=0.8)
    F.save(fig, "fig_p4_uq")


if __name__ == "__main__":
    with plt.rc_context(STYLE):
        false_turn_on(); parity(); sobol(); uq_combined(); pareto()
