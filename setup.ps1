# One-time toolchain setup on Windows 10/11, no admin rights needed.
# Prerequisites: Python 3.14 (python.org), Git, KiCad 10 (kicad.org; default install path).
# Installs into this folder: env\ (Python venv), tools\openEMS_dl, tools\ngspice_dl,
# tools\bin\fasthenry.exe, tools\dcdc-parasitics and tools\dcdc_diag (baseline extractor).
$ErrorActionPreference = 'Stop'
Set-Location $PSScriptRoot
New-Item -ItemType Directory -Force tools | Out-Null

# 1) Python venv + pinned packages (ziglang provides the C compiler for FastHenry)
if (-not (Test-Path env)) { python -m venv env }
$py = '.\env\Scripts\python.exe'
& $py -m pip install --upgrade pip
& $py -m pip install -r requirements.txt

# 2) openEMS 0.37.0 (official Windows build, includes CSXCAD/openEMS Python wheels)
if (-not (Test-Path tools\openEMS_dl)) {
  Invoke-WebRequest 'https://github.com/thliebig/openEMS-Project/releases/download/v0.37.0/openEMS_x64_v0.37.0_msvc.zip' -OutFile tools\openEMS.zip
  Expand-Archive tools\openEMS.zip -DestinationPath tools\openEMS_dl -Force
}
$tag = (& $py -c "import sys;print(f'cp{sys.version_info[0]}{sys.version_info[1]}')")
& $py -m pip install (Get-ChildItem "tools\openEMS_dl\openEMS\python\csxcad-*-$tag-*.whl").FullName (Get-ChildItem "tools\openEMS_dl\openEMS\python\openems-*-$tag-*.whl").FullName

# 3) ngspice 47 (official Windows build)
if (-not (Test-Path tools\ngspice_dl)) {
  curl.exe -sL -o tools\ngspice-47_64.7z 'https://downloads.sourceforge.net/project/ngspice/ng-spice-rework/47/ngspice-47_64.7z'
  & $py -c "import py7zr;py7zr.SevenZipFile(r'tools\ngspice-47_64.7z').extractall(r'tools\ngspice_dl')"
}

# 4) FastHenry2 from source (32-bit build, see scripts\build_fasthenry.ps1 for why)
& .\scripts\build_fasthenry.ps1

# 5) Baseline extractor fl4p/dcdc-parasitics at the commit used in the study.
#    patches\dcdc-parasitics_kicad10.patch : two one-line KiCad 10 API fixes (PCB_VIA::GetWidth needs a layer)
#    patches\dcdc_diag_all_pad_terminals.patch : diagnostic variant that connects every pad of a device net
$dcdcCommit = '12962bba11556a02ea4a3fd1c45ca3a476d5cb26'
if (-not (Test-Path tools\dcdc-parasitics)) {
  git clone https://github.com/fl4p/dcdc-parasitics.git tools\dcdc-parasitics
  git -C tools\dcdc-parasitics checkout $dcdcCommit
  git -C tools\dcdc-parasitics apply --ignore-whitespace ..\..\patches\dcdc-parasitics_kicad10.patch
}
if (-not (Test-Path tools\dcdc_diag)) {
  git clone tools\dcdc-parasitics tools\dcdc_diag
  git -C tools\dcdc_diag checkout $dcdcCommit
  git -C tools\dcdc_diag apply --ignore-whitespace ..\..\patches\dcdc-parasitics_kicad10.patch
  git -C tools\dcdc_diag apply --ignore-whitespace ..\..\patches\dcdc_diag_all_pad_terminals.patch
}

# 6) check everything
& $py scripts\verify_install.py
