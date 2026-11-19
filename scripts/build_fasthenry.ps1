# Build FastHenry2 (ediloren/FastHenry2) on Windows with the zig C compiler (pip install ziglang).
#
# Why 32-bit: the 1990s FastHenry C code assumes int/long can hold pointers and relies on
# implicit function declarations. On 64-bit Windows (LLP64, long = 32 bit) this truncates
# pointers and the binary crashes with an access violation. A 32-bit build is ABI-faithful
# to what the code expects. The 2 GB 32-bit address space is
# far more than the PCB power-loop problems in this project need.
#
# Patches (applied to a copy in tools/fh_build, upstream clone untouched):
#   1. induct.c: the old SysV `matherr` hook is renamed (it clashes with MinGW's _matherr macro).
#   2. win_shim.c: stubs for sbrk() and gethostname(), used only for diagnostics printouts.
$ErrorActionPreference = 'Stop'
$root = Split-Path -Parent $PSScriptRoot
$tools = Join-Path $root 'tools'
$py = Join-Path $root 'env\Scripts\python.exe'
if (-not (Test-Path "$tools\FastHenry2")) {
    git clone --depth 1 https://github.com/ediloren/FastHenry2.git "$tools\FastHenry2"
}
$build = "$tools\fh_build"
if (Test-Path $build) { Remove-Item -Recurse -Force $build }
Copy-Item -Recurse "$tools\FastHenry2\src\fasthenry" $build

$f = "$build\induct.c"
$c = Get-Content $f -Raw
$c = $c -replace 'int matherr\(exc\)\r?\nstruct exception \*exc;', "int fh_matherr_unused(exc)`nstruct _exception *exc;"
Set-Content $f $c -NoNewline -Encoding ascii

@'
/* Windows shim for FastHenry2 build (diagnostic-only POSIX calls). */
#include <string.h>
static char fh_dummy_heap;
char *sbrk(long inc) { (void)inc; return &fh_dummy_heap; }
int gethostname(char *name, int len) { strncpy(name, "windows", len); return 0; }
'@ | Set-Content "$build\win_shim.c" -Encoding ascii

Push-Location $build
try {
    $sp = 'spAllocate','spBuild','spFactor','spOutput','spSolve','spUtils','spCompat','spFortran' | ForEach-Object { "sparse\$_.c" }
    $skip = 'makeMakefile.c','bigmeshPre.c','bigmeshPre_direct.c','newPrecond_cube.c'
    $main = (Get-ChildItem -Filter *.c | Where-Object Name -notin $skip).FullName
    New-Item -ItemType Directory -Force "$tools\bin" | Out-Null
    $flags = @('-target','x86-windows-gnu','-O2','-std=gnu89','-fcommon','-w',
               '-Wno-error','-Wno-implicit-function-declaration','-Wno-int-conversion',
               '-Wno-incompatible-pointer-types','-Wno-implicit-int','-Wno-return-type',
               '-Isparse','-o',"$tools\bin\fasthenry.exe")
    & $py -m ziglang cc @flags @main @sp -lm
    if ($LASTEXITCODE -ne 0) { throw "FastHenry build failed" }
    Write-Host "Built $tools\bin\fasthenry.exe"
} finally { Pop-Location }
