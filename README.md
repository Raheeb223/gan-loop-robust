# gan-loop-robust

**Open-source, verified extraction of GaN half-bridge PCB power-loop inductance, from KiCad layout to switching behaviour.**

This repository is the reproducibility package for the accompanying manuscript (under review).
It contains all code, input data, generated KiCad boards, result tables and figures. One
command regenerates every number, table and figure using only free tools:
**KiCad 10, FastHenry 2, openEMS 0.37, ngspice 47 and Python 3.14**.

One parametric layout description (`scripts/layouts.py`) generates both
* manufacturable 4-layer **KiCad boards** (DRC-clean, with Gerbers and drill files), and
* **solver models** for FastHenry (PEEC) and openEMS (FDTD) with identical copper.

An **ngspice** model of the EPC2015C, built only from its datasheet, turns the extracted
inductances into switching waveforms.

## Verification summary

| Tier | Reference | Result |
|---|---|---|
| 1a | Exact 2-D geometric-mean-distance solution; Grover bar formula | FastHenry within 0.02 % and 0.13 % (DC) |
| 1b | Published commercial FEA (Letellier *et al.*, IEEE TPEL 2018), three power loops | FastHenry within 0.9–2.6 % |
| 1c | openEMS full-wave FDTD (independent method) | within 0.4–5.2 % of FastHenry |
| 2  | EPC commercial simulation (WP010 / Reusch & Strydom), three layout families | optimal −1.2 %, vertical +2.5 % (31/4 mil); all trends reproduced |
| 3  | EPC measured overshoot (8 boards) | calibrated on one layout family, other two predicted: hold-out MAE 7.7 percentage points |

Known differences larger than 10 % between this work and published data, with their
assessment, are listed in [`results/discrepancies.md`](results/discrepancies.md).

## Requirements

* Windows 10/11 (the scripts are PowerShell; the Python code itself is portable)
* Python 3.14, Git, KiCad 10 (default install location)
* About 2 GB of free disk space. No administrator rights are needed.

## Quick start

```powershell
powershell -ExecutionPolicy Bypass -File setup.ps1     # venv, openEMS, ngspice, FastHenry build, baseline tool
powershell -ExecutionPolicy Bypass -File run_all.ps1   # every table and figure
powershell -ExecutionPolicy Bypass -File run_all.ps1 -SkipOpenEMS   # same, without the ~2 h FDTD cross-check
```

`setup.ps1` installs the pinned Python packages from `requirements.txt` and compiles FastHenry
from source using the `ziglang` pip package (see `scripts/build_fasthenry.ps1` for why it is a
32-bit build). Tool versions are recorded in `results/provenance.json` on every run.

## Where each result comes from

| Result | Script(s) | Output | Runtime* |
|---|---|---|---|
| Closed-form references | `analytic.py` | console | < 1 min |
| FastHenry vs published FEA | `canonical_fh.py` | `results/tables/canonical_fh.csv` | 5 min |
| openEMS vs FastHenry | `canonical_oems.py` | `results/tables/canonical_oems.csv` | ~2 h |
| Comparison with EPC simulation, high-frequency convergence | `sweep_epc.py`, `hf_ladder.py` | `epc_sweep_hf.csv`, `epc_mesh.csv`, `hf_ladder.csv` | 26 min |
| KiCad boards and DRC | `kicad_prims.py`, `kicad_all.py` | `kicad/boards/`, `kicad_drc.csv` | few min |
| Existing open extractor (baseline) | `bench_dcdc.py` | `dcdc_benchmark.csv` | 10 min |
| Bias and uncertainty of ringing-based inductance | `ring_uncertainty.py` | `ring_bias.csv`, `ring_mc.csv` | 1 min |
| Circuit-level validation against measured boards | `circuit_validation.py` | `circuit_validation.csv`, `circuit_calibration_grid.csv` | 6 min |
| Design map for standard fabrication stack-ups | `design_map.py`, `design_map_overshoot.py` | `design_map.csv`, `design_map_points.csv` | 13 min |
| Figures | `figures.py` | `results/figures/*.pdf, *.png` | < 1 min |

\*Measured on a 4-core Intel i5-1035G1 laptop with 12 GB RAM. Fixed seeds are used throughout;
re-running the pipeline reproduces the stored tables bit-for-bit (solver runtimes excepted).

## Repository layout

```
scripts/          all code; each file starts with a docstring explaining the method and assumptions
  analytic.py       closed-form references (Grover bar, exact 2-D loop, ringing formula)
  geom.py           solver-neutral copper geometry and FastHenry exporter
  oems.py           openEMS model builder and impedance extraction
  layouts.py        EPC-style lateral / vertical / optimal power loops (parametric)
  kicad_prims.py    manufacturable copper primitives; kicad_build.py -> .kicad_pcb
  epc_model.py      datasheet-derived EPC2015C ngspice model; spice_sim.py -> commutation bench
data/published/   every external number used, with source, location and confidence (see SOURCES.md)
results/tables/   every computed number (CSV)
results/figures/  figures (PDF and PNG)
kicad/boards/     16 generated boards with DRC reports and fabrication files
patches/          changes applied to the third-party baseline tool
```

## Third-party software

* **fl4p/dcdc-parasitics** (baseline extractor) is cloned by `setup.ps1` at commit
  `12962bb` and is not redistributed here. Two patches are applied:
  * `patches/dcdc-parasitics_kicad10.patch`: two one-line KiCad 10 API compatibility fixes
    (`PCB_VIA::GetWidth` needs a layer argument);
  * `patches/dcdc_diag_all_pad_terminals.patch`: a diagnostic variant in which every pad of
    a device net is connected (used only for the "all pads connected" column of the baseline
    comparison).
* FastHenry 2 (MIT-style licence), openEMS (GPL-3), ngspice (BSD) and KiCad (GPL-3) are
  downloaded or built by `setup.ps1` and keep their own licences.

## Source documents

Published values were transcribed or pixel-digitised from EPC application material and IEEE
papers. The source documents are copyrighted and not included. The digitised values are, so
every downstream result is reproducible without them. See
[`data/published/SOURCES.md`](data/published/SOURCES.md) to re-run the digitisation step.

## Licence

Code: [MIT](LICENSE). Generated data, tables, figures and board files:
[CC BY 4.0](LICENSE-DATA).

## Citation

See [`CITATION.cff`](CITATION.cff). A DOI will be added on publication.
