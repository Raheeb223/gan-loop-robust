"""Print the version of every tool in the chain and fail loudly if one is missing.
Also writes results/provenance.json (tool versions, platform, timestamp)."""
from __future__ import annotations

import json
import os
import platform
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(Path(__file__).parent))
KICAD_BIN = Path(os.environ.get("KICAD_BIN", r"C:\Program Files\KiCad\10.0\bin"))


def out(cmd, **kw):
    try:
        p = subprocess.run(cmd, capture_output=True, text=True, timeout=120, **kw)
        return (p.stdout + p.stderr).strip()
    except Exception as ex:
        return f"ERROR {ex}"


def main():
    prov = {"timestamp_utc": datetime.now(timezone.utc).isoformat(), "platform": platform.platform(),
            "python": sys.version}
    ok = True
    import fhrun
    s = ROOT / "results" / "raw" / "_verify"
    s.mkdir(parents=True, exist_ok=True)
    try:
        import analytic
        import geom
        r = fhrun.run(geom.to_fasthenry(geom.straight_bar(10, 1, 0.035), (1e3, 1e3)), s / "fh")
        ref = analytic.bar_partial_self_L(10e-3, 1e-3, 35e-6)
        prov["fasthenry"] = {"path": str(fhrun.FASTHENRY), "bar_L_nH": r["L"][0] * 1e9, "ref_nH": ref * 1e9}
        ok &= abs(r["L"][0] / ref - 1) < 0.01
    except Exception as ex:
        prov["fasthenry"] = f"MISSING/FAILED: {ex}"; ok = False
    import spice_sim
    v = out([str(spice_sim.NGSPICE), "-v"])
    prov["ngspice"] = next((ln.strip(" *") for ln in v.splitlines() if "ngspice-" in ln.lower()), "MISSING")
    ok &= "ngspice" in v.lower()
    try:
        import oems  # noqa: F401
        import openEMS
        prov["openEMS_python"] = getattr(openEMS, "__version__", "import OK")
    except Exception as ex:
        prov["openEMS_python"] = f"MISSING: {ex}"; ok = False
    prov["kicad_cli"] = out([str(KICAD_BIN / "kicad-cli.exe"), "version"])
    prov["kicad_python_pcbnew"] = out([str(KICAD_BIN / "python.exe"), "-c", "import pcbnew;print(pcbnew.Version())"])
    ok &= prov["kicad_cli"][:1].isdigit()
    import numpy, scipy, pandas, matplotlib  # noqa: E401
    prov["numpy"], prov["scipy"], prov["pandas"], prov["matplotlib"] = (numpy.__version__, scipy.__version__,
                                                                        pandas.__version__, matplotlib.__version__)
    gc = out(["git", "-C", str(ROOT), "rev-parse", "HEAD"])
    prov["git_commit"] = gc if len(gc) == 40 else "no commits yet"
    (ROOT / "results" / "provenance.json").write_text(json.dumps(prov, indent=2))
    print(json.dumps(prov, indent=2))
    if not ok:
        sys.exit("toolchain check FAILED")


if __name__ == "__main__":
    main()
