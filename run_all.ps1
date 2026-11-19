# Regenerate every result, table and figure from scratch (Windows / PowerShell).
# Usage:  powershell -ExecutionPolicy Bypass -File run_all.ps1 [-SkipOpenEMS]
# Approximate run time on a 4-core laptop (i5-1035G1): FastHenry stages ~1.5 h,
# openEMS cross-check ~2 h (skip with -SkipOpenEMS), SPICE ~20 min.
param([switch]$SkipOpenEMS)
$ErrorActionPreference = 'Stop'
Set-Location $PSScriptRoot
$py = '.\env\Scripts\python.exe'
function Step($name, [scriptblock]$b) { Write-Host "=== $name" -ForegroundColor Cyan; & $b; if ($LASTEXITCODE) { throw "$name failed" } }

Step 'verify toolchain'            { & $py scripts\verify_install.py }
# Digitisation needs the copyrighted source PDFs (see data\published\SOURCES.md); the digitised CSVs are in the repo.
if (Test-Path data\published\sources\Reusch_APEC2013_slides.pdf) {
  Step 'digitise published EPC data' { & $py scripts\digitise_epc.py; & $py scripts\build_measured_table.py }
} else {
  Step 'build tables from digitised EPC data' { & $py scripts\build_measured_table.py }
}
Step 'tier-1 analytic self-check'  { & $py scripts\analytic.py }
Step 'tier-1 FastHenry canonical'  { & $py scripts\canonical_fh.py }
if (-not $SkipOpenEMS) { Step 'tier-1 openEMS canonical' { & $py scripts\canonical_oems.py } }
Step 'tier-2 WP010 sweep (coarse + mesh ladder)' { & $py scripts\sweep_epc.py --mesh }
Step 'tier-2 HF convergence ladder' { & $py scripts\hf_ladder.py }
Step 'KiCad boards, DRC, Gerbers'  { & $py scripts\kicad_prims.py kicad\prims; & $py scripts\kicad_all.py }
Step 'baseline: dcdc-parasitics'   { & $py scripts\bench_dcdc.py }
Step 'ringing-method uncertainty'  { & $py scripts\ring_uncertainty.py; & $py scripts\ring_sensitivity.py }
Step 'circuit validation'          { & $py scripts\circuit_validation.py }
Step 'fab design map + Monte Carlo' { & $py scripts\design_map.py; & $py scripts\design_map_overshoot.py }
Step 'figures'                     { & $py scripts\figures.py }
# ---- extensions (results/phase2 ... phase6); see each script's docstring for method and runtime
Step 'P2 calibration cross-validation, dv/dt definitions' { & $py scripts\p2_crossval.py }
Step 'P3 coupled power/gate extraction (~25 min)' { & $py scripts\p3_extract.py }
Step 'P3 RL ladder, K-coupled SPICE, false turn-on' { & $py scripts\p3_spice.py }
if (-not $SkipOpenEMS) {
  Step 'P3 FastHenry thin-copper twin of the openEMS models' { & $py scripts\p3_fh_twin.py }
  Step 'P3 openEMS coupled cross-check (4 runs, ~1 h each)' { & $py scripts\p3_openems.py }
}
Step 'P6 newer devices (EPC2218, EPC2302)' { & $py scripts\p6_devices.py }
Step 'P4 FastHenry campaign (~424 runs, ~1.5 h)' { & $py scripts\p4_campaign.py }
Step 'P4 surrogate, conformal intervals, Sobol' { & $py scripts\p4_uq.py }
Step 'P5 NSGA-II robust optimisation + KiCad boards' { & $py scripts\p5_opt.py }
if (Test-Path paper) { Step 'LaTeX source + Overleaf zip' { & $py scripts\make_latex.py; & $py scripts\make_overleaf_zip.py; & $py scripts\make_ojpel_zip.py } }
Write-Host 'done' -ForegroundColor Green
