"""Manuscript figures added for the extended study (workflow overview).

Output: results/figures/fig_workflow_v2.(pdf|png)
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
import figures as F  # noqa: E402  (shared style, colours and save())

ARROW = "→"


def fig_workflow_v2():
    from matplotlib.patches import FancyBboxPatch
    plt = F.plt
    fig, ax = plt.subplots(figsize=(7.16, 2.6))
    ax.set_xlim(0, 100); ax.set_ylim(0, 35); ax.axis("off")

    def box(x, y, w, h, t):
        ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.3,rounding_size=1.2", fc=F.GRID, ec=F.INK2,
                                    lw=0.8))
        ax.text(x + w / 2, y + h / 2, t, ha="center", va="center", fontsize=5.8, color=F.INK)

    def arr(x0, y0, x1, y1):
        ax.annotate("", (x1, y1), (x0, y0), arrowprops=dict(arrowstyle="-|>", color=F.INK2, lw=0.8))

    box(0.5, 10.5, 15, 14, "Parametric layout:\nfamily, stackup,\ndevice footprint,\ngate return,\nfabrication\ntolerances")
    box(21, 26, 17, 7, "KiCad 10 board\nDRC, Gerbers, drill")
    box(21, 15, 17, 8, "FastHenry (PEEC)\n3-port Z(f): power,\ngate HS, gate LS")
    box(21, 4, 17, 8, "openEMS (FDTD)\nindependent check")
    box(44, 20, 18, 9, "L matrix, passive\nRL ladder,\nK-coupled ngspice")
    box(44, 4, 18, 12, f"FastHenry campaign:\ndesign + tolerance\nvariables {ARROW} GP\nsurrogate, conformal\nerror bars")
    box(67.5, 24, 15, 9, "Overshoot,\ndv/dt, false\nturn-on margin")
    box(67.5, 13, 15, 7, "Sobol indices,\ntolerance spread")
    box(67.5, 1.5, 15, 9, f"NSGA-II robust\nlayout {ARROW} verify\n{ARROW} KiCad boards")
    box(87, 7.5, 12.5, 19, "Verified against:\nexact 2-D,\nGrover,\npublished FEA,\nEPC simulation,\nmeasured\nboards")
    arr(15.8, 21, 21, 29); arr(15.8, 18, 21, 19); arr(15.8, 14, 21, 8)
    arr(38, 19, 44, 24); arr(38, 17, 44, 12); arr(38, 8, 44, 21)
    arr(62, 26, 67.5, 28); arr(62, 11, 67.5, 16.5); arr(62, 7, 67.5, 6)
    arr(82.5, 6, 87, 11); arr(82.5, 28, 87, 23)
    F.save(fig, "fig_workflow_v2")


if __name__ == "__main__":
    fig_workflow_v2()
