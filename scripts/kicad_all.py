"""Build, DRC-check and export fabrication data for every generated layout.

For each primitives JSON in kicad/prims/:
  kicad/boards/<tag>/<tag>.kicad_pcb (+ .kicad_pro)
  kicad/boards/<tag>/drc.json
  kicad/boards/<tag>/fab/  (Gerbers + Excellon drill)
  kicad/boards/<tag>/<tag>_F.svg, _In1.svg, _B.svg
Summary: results/tables/kicad_drc.csv
"""
from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
KICAD_BIN = Path(os.environ.get("KICAD_BIN", r"C:\Program Files\KiCad\10.0\bin"))
KPY = KICAD_BIN / "python.exe"
CLI = KICAD_BIN / "kicad-cli.exe"


def sh(*a):
    p = subprocess.run([str(x) for x in a], capture_output=True, text=True)
    if p.returncode != 0:
        raise RuntimeError(f"{a[0]} failed: {p.stderr[-800:]}")
    return p.stdout


def main():
    rows = []
    for js in sorted((ROOT / "kicad" / "prims").glob("*.json")):
        tag = js.stem
        bd = ROOT / "kicad" / "boards" / tag
        (bd / "fab").mkdir(parents=True, exist_ok=True)
        pcb = bd / f"{tag}.kicad_pcb"
        sh(KPY, ROOT / "scripts" / "kicad_build.py", js, pcb)
        sh(CLI, "pcb", "drc", "--format", "json", "--severity-error", "-o", bd / "drc.json", pcb)
        d = json.load(open(bd / "drc.json"))
        nviol = len(d.get("violations", []))
        nunc = len(d.get("unconnected_items", []))
        sh(CLI, "pcb", "export", "gerbers", "-o", str(bd / "fab") + os.sep, pcb)
        sh(CLI, "pcb", "export", "drill", "-o", str(bd / "fab") + os.sep, pcb)
        for lay, suf in (("F.Cu,Edge.Cuts", "F"), ("In1.Cu,Edge.Cuts", "In1"), ("B.Cu,Edge.Cuts", "B")):
            sh(CLI, "pcb", "export", "svg", "--layers", lay, "--mode-single", "--fit-page-to-board",
               "-o", bd / f"{tag}_{suf}.svg", pcb)
        nfab = len(list((bd / "fab").iterdir()))
        rows.append(dict(tag=tag, drc_violations=nviol, unconnected=nunc, fab_files=nfab))
        print(tag, "violations", nviol, "unconnected", nunc, "fab files", nfab, flush=True)
    pd.DataFrame(rows).to_csv(ROOT / "results" / "tables" / "kicad_drc.csv", index=False)


if __name__ == "__main__":
    main()
