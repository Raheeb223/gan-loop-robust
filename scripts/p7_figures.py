"""Manuscript figures added for the extended study (workflow overview).

Output: results/figures/fig_workflow_v2.(pdf|png)
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
import figures as F  # noqa: E402  (shared style and save())

ARROW = "→"


def fig_workflow_v2():
    """Black-and-white block diagram of the workflow (IEEE two-column width)."""
    plt = F.plt
    rc = {"font.family": "serif", "font.serif": ["Times New Roman", "Times", "STIXGeneral", "DejaVu Serif"],
          "mathtext.fontset": "stix", "pdf.fonttype": 42, "ps.fonttype": 42}
    with plt.rc_context(rc):
        fig, ax = plt.subplots(figsize=(7.16, 2.55))
        ax.set_xlim(0, 100); ax.set_ylim(0, 35); ax.axis("off")
        B = {}

        def box(key, x, y, w, h, title, body, lw=0.7):
            ax.add_patch(plt.Rectangle((x, y), w, h, fc="white", ec="black", lw=lw))
            B[key] = (x, y, w, h)
            txt = title if not body else f"{title}\n{body}"
            t = ax.text(x + w / 2, y + h / 2, txt, ha="center", va="center", fontsize=7, color="black",
                        linespacing=1.15)
            if body:  # bold first line only: draw title separately
                t.remove()
                n = body.count("\n") + 2
                lh = 2.05  # line height in data units at 7 pt
                y0 = y + h / 2 + (n - 1) * lh / 2
                ax.text(x + w / 2, y0, title, ha="center", va="center", fontsize=7, fontweight="bold", color="black")
                ax.text(x + w / 2, y0 - lh, body, ha="center", va="top", fontsize=7, color="black",
                        linespacing=1.15)

        def R(k, f=0.5):  # point on right edge
            x, y, w, h = B[k]; return (x + w, y + f * h)

        def L(k, f=0.5):  # point on left edge
            x, y, w, h = B[k]; return (x, y + f * h)

        def arr(p0, p1, style="-|>"):
            ax.annotate("", p1, p0, arrowprops=dict(arrowstyle=style, color="black", lw=0.7,
                                                    shrinkA=0, shrinkB=0, mutation_scale=7))

        box("par", 0, 8.5, 15.8, 18, "Parametric layout", "family, stackup,\ndevice footprint,\ngate return,\nfabrication\ntolerances")
        box("kic", 19.5, 26.5, 18, 8, "KiCad 10 board", "DRC, Gerbers, drill")
        box("fh", 19.5, 13.5, 18, 9.5, "FastHenry (PEEC)", "3-port Z(f): power,\ngate HS, gate LS")
        box("oem", 19.5, 0.5, 18, 8, "openEMS (FDTD)", "independent check")
        box("lm", 42.5, 21.5, 19, 11, "Circuit model", "L matrix, passive RL\nladder, K-coupled\nngspice")
        box("cmp", 42.5, 0.5, 19, 14.5, "Solver campaign", f"design + tolerance\nvariables {ARROW} GP\nsurrogate, conformal\nerror bars")
        box("ovs", 66.5, 24, 16, 10.5, "Switching", "overshoot, dv/dt,\nfalse turn-on margin")
        box("sob", 66.5, 12.5, 16, 8.5, "Tolerance analysis", "Sobol indices,\nspread")
        box("opt", 66.5, 0.5, 16, 9.5, "Robust layout", f"NSGA-II {ARROW} verify\n{ARROW} KiCad boards")
        box("ver", 86.3, 5.5, 13.7, 24, "Verified against", "exact 2-D,\nGrover,\npublished FEA,\nEPC simulation,\nmeasured\nboards")

        arr(R("par", 0.85), L("kic")); arr(R("par", 0.5), L("fh")); arr(R("par", 0.15), L("oem"))
        x_m = 19.5 + 9
        arr((x_m, 8.5), (x_m, 13.5), "<|-|>")  # openEMS cross-checks FastHenry
        arr(R("fh", 0.7), L("lm", 0.4)); arr(R("fh", 0.3), L("cmp", 0.7))
        arr(R("lm"), L("ovs")); arr(R("cmp", 0.75), L("sob")); arr(R("cmp", 0.3), L("opt"))
        arr(R("ovs"), L("ver", 0.8)); arr(R("sob"), L("ver", 0.5)); arr(R("opt"), L("ver", 0.2))
        F.save(fig, "fig_workflow_v2")


if __name__ == "__main__":
    fig_workflow_v2()
