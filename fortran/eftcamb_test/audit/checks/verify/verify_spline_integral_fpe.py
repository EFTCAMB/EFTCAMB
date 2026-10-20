"""Verifier (spline-integral-convention): sub-claim "with the default nodes (x1 = 0) the cubic init divides by zero
(FPE trap at 02_interpolation_linear_1D.f90:600)".

Compiles a driver against the AD (debug, -ffpe-trap=invalid,zero,overflow) build with the same trap flags and only
initialises a cubic spline (flag 9) from an ini that gives no node positions (code default x_i = (i-1)/(N-1)), i.e.
init_func_from_file + init_from_file, no integral() call. Two data sets: w = -0.9 (1+b > 0 -> 0**negative) and
Omega-like data y = x (b = 0 -> 0**-3) to show it hits any use of flag 9 with default nodes, not just w(a).
usage: /Users/marco/Software/miniforge3/bin/python verify_spline_integral_fpe.py
"""

import subprocess
from pathlib import Path

AUDIT = Path("/private/tmp/claude-501/-Users-marco-Desktop-IW-eftcamb-update-EFTCAMB/"
             "2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad/audit")
BUILD = AUDIT / "build" / "AD"
WORK = AUDIT / "work" / "verify" / "spline_integral" / "fpe"
LIB = "/Users/marco/Software/miniforge3/lib"

SRC = r"""
program vfpe
use precision
use IniObjects
use EFTCAMB_parametrizations_1D
implicit none
class(parametrized_function_1D), allocatable :: f
type(TIniFile) :: Ini
integer :: err
err = 0
call Ini%Open('v.ini')
call allocate_parametrized_1D_function(f, 9, 'EFTw', 'w', err, 0)
call f%init_func_from_file(Ini, err)
call f%init_from_file(Ini, err)
write(*, *) 'INIT DONE, value(0.5) =', f%value(0.5_dl)
end program
"""


def main():
    WORK.mkdir(parents=True, exist_ok=True)
    (WORK / "vfpe.f90").write_text(SRC)
    subprocess.run(["gfortran", "-g", "-O0", "-fbacktrace", "-ffpe-trap=invalid,zero,overflow", "-cpp",
                    "-ffree-line-length-none", "-fopenmp", f"-I{BUILD}/camb/Releaselib",
                    f"-I{BUILD}/forutils/Releaselib", "vfpe.f90", "-o", "vfpe",
                    str(AUDIT / "trees" / "AD" / "camb" / "camblib.so"), f"-Wl,-rpath,{LIB}"],
                   cwd=WORK, check=True, capture_output=True)
    for name, vals in (("w=-0.9", [-0.9] * 6), ("y=x (Omega-like)", [0.0, 0.2, 0.4, 0.6, 0.8, 1.0])):
        ini = "EFTw_Spline_Pixels = 6\nEFTw_null_value = 0.0\n" + "".join(
            f"EFTwv{i + 1} = {v}\n" for i, v in enumerate(vals))
        (WORK / "v.ini").write_text(ini)
        p = subprocess.run([str(WORK / "vfpe")], cwd=WORK, capture_output=True, text=True)
        frames = [ln.strip() for ln in p.stderr.splitlines() if ".f90:" in ln or "Floating" in ln
                  or "Program received" in ln]
        print(f"{name}: rc={p.returncode} stdout={p.stdout.strip()[:80]!r}")
        for fr in frames[:6]:
            print("    ", fr[:200])


if __name__ == "__main__":
    main()
