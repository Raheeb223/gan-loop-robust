"""Baseline: the existing open KiCad->FastHenry extractor fl4p/dcdc-parasitics on our boards.

Configurations (all with lead_mm=0, i.e. copper-only, as its README recommends for
lead-inclusive device models):
  A  published tool, default input-capacitor port (one capacitor)
  B  published tool, 4-capacitor bank declared (--cin-loop-refs C1 C2 C3 C4)
  C  diagnostic copy: every drain/source pad of each FET connected (ideal on-die busbar), bank
Only two compatibility edits were made to the published tool for KiCad 10
(PCB_VIA::GetWidth() now needs a layer): lib/kicad_geom.py line 942, copper_dump.py line 68.
Output: results/tables/dcdc_benchmark.csv
"""
from __future__ import annotations

import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
ENV = dict(os.environ, PYTHONUTF8="1", KICAD_PY=r"C:\Program Files\KiCad\10.0\bin\python.exe",
           FASTHENRY=str(ROOT / "tools" / "bin" / "fasthenry.exe"))
BOARDS = {
    "optimal_31_4": "optimal_T0.787_d0.102_cu0.070_v5_p0.20",
    "vertical_31_4": "vertical_T0.787_d0.102_cu0.070_v5_p0.20",
}
CONFIGS = {
    "A_published_default": ("dcdc-parasitics", []),
    "B_published_bank": ("dcdc-parasitics", ["--cin-loop-refs", "C1", "C2", "C3", "C4"]),
    "C_all_pads_bank": ("dcdc_diag", ["--cin-loop-refs", "C1", "C2", "C3", "C4"]),
}
PITCHES = ["0.5", "0.25"]


def main():
    rows = []
    for bname, tag in BOARDS.items():
        pcb = ROOT / "kicad" / "boards" / tag / f"{tag}.kicad_pcb"
        for cname, (tool, extra) in CONFIGS.items():
            out = ROOT / "results" / "raw" / "dcdc_bench" / f"{bname}_{cname}"
            shutil.rmtree(out, ignore_errors=True)
            cmd = [sys.executable, str(ROOT / "tools" / tool / "extract_parasitics.py"), str(pcb),
                   "--sw", "SW", "--gnd", "GND", "--vin", "VIN", "--hs-ref", "Q_T", "--ls-ref", "Q_S",
                   "--hs-gate", "GH", "--ls-gate", "GL", "--allow-missing-gate-ports", "--lead-mm", "0",
                   "--cu-thickness", "0.07", "--pitch", *PITCHES, *extra, "-o", str(out)]
            try:
                p = subprocess.run(cmd, env=ENV, capture_output=True, text=True, timeout=2400)
                log = p.stdout + p.stderr
            except subprocess.TimeoutExpired as ex:
                log = (ex.stdout or "") if isinstance(ex.stdout, str) else ""
                log += "\nTIMEOUT"
            out.mkdir(parents=True, exist_ok=True)
            (out / "run.log").write_text(log, encoding="utf-8")
            found = re.findall(r"pitch\s+([\d.]+) mm : L_loop=\s*([\d.]+) nH", log)
            for pitch, L in found:
                rows.append(dict(board=bname, config=cname, pitch_mm=float(pitch), L_loop_nH=float(L)))
            if not found:
                err = [ln for ln in log.splitlines() if re.search(r"rror|invalid|TIMEOUT|failed", ln)]
                rows.append(dict(board=bname, config=cname, pitch_mm=None, L_loop_nH=None,
                                 error=(err[-1] if err else "no result")[:200]))
            print(bname, cname, found or rows[-1].get("error"), flush=True)
            pd.DataFrame(rows).to_csv(ROOT / "results" / "tables" / "dcdc_benchmark.csv", index=False)


if __name__ == "__main__":
    main()
