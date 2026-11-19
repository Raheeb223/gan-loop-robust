"""Thin wrapper around the FastHenry executable: run an .inp file and parse Zc.mat.

FastHenry writes the impedance matrix to Zc.mat in the working directory, one
matrix per frequency, in the format:

    Impedance matrix for frequency = 1000 1 x 1
      0.00123  +6.97e-06j

We run each job in its own directory so parallel runs never clash.
"""
from __future__ import annotations

import os
import re
import subprocess
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
FASTHENRY = Path(os.environ.get("FASTHENRY", ROOT / "tools" / "bin" / "fasthenry.exe"))

_num = r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?"


def parse_zc(path: Path) -> dict[float, np.ndarray]:
    """Return {frequency_Hz: complex impedance matrix} from a FastHenry Zc.mat file."""
    out: dict[float, np.ndarray] = {}
    lines = Path(path).read_text().splitlines()
    i = 0
    while i < len(lines):
        m = re.match(rf"Impedance matrix for frequency = ({_num})\s+(\d+) x (\d+)", lines[i].strip())
        if not m:
            i += 1
            continue
        f, n, k = float(m.group(1)), int(m.group(2)), int(m.group(3))
        rows = []
        for r in range(n):
            vals = re.findall(rf"({_num})\s*([-+]\s*{_num})j", lines[i + 1 + r])
            rows.append([complex(float(a), float(b.replace(" ", ""))) for a, b in vals])
        z = np.array(rows, dtype=complex)
        assert z.shape == (n, k), f"bad Zc.mat block at line {i}: {z.shape} vs {(n, k)}"
        out[f] = z
        i += 1 + n
    if not out:
        raise RuntimeError(f"No impedance matrices found in {path}")
    return out


def run(inp_text: str, workdir: Path, timeout: float = 7200.0, extra_args: list[str] | None = None) -> dict:
    """Write `inp_text` to workdir/model.inp, run FastHenry, return parsed results.

    Returns dict with keys: freqs (array), Z (array of matrices), L (self-L of port 1, H),
    R (ohm), runtime_s, stdout.
    """
    workdir = Path(workdir)
    workdir.mkdir(parents=True, exist_ok=True)
    inp = workdir / "model.inp"
    inp.write_text(inp_text)
    zc = workdir / "Zc.mat"
    if zc.exists():
        zc.unlink()
    t0 = time.perf_counter()
    proc = subprocess.run(
        [str(FASTHENRY), *(extra_args or []), "model.inp"],
        cwd=workdir, capture_output=True, text=True, timeout=timeout,
    )
    dt = time.perf_counter() - t0
    (workdir / "fasthenry.log").write_text(proc.stdout + "\n--- stderr ---\n" + proc.stderr)
    if proc.returncode != 0 or not zc.exists():
        raise RuntimeError(f"FastHenry failed (code {proc.returncode}) in {workdir}; see fasthenry.log")
    zmap = parse_zc(zc)
    freqs = np.array(sorted(zmap))
    Z = np.array([zmap[f] for f in freqs])
    z11 = Z[:, 0, 0]
    return {
        "freqs": freqs,
        "Z": Z,
        "L": z11.imag / (2 * np.pi * freqs),
        "R": z11.real,
        "runtime_s": dt,
        "stdout": proc.stdout,
    }
