"""Repro (t4a-param): Pade parametrization (allocator flag 6) fourth_derivative drops the u' factor of the 4 u' g'''
Leibniz term (04p010_pade_series_1D.f90:419: `4._dl*(-ddddown/down**2 ...)` should be `4._dl*dup*(...)`).
Expected values are the exact fourth derivative of the same rational function computed with sympy.

usage: /Users/marco/Software/miniforge3/bin/python pade-fourth-derivative.py
"""

import subprocess
from pathlib import Path

import sympy as sp

AUDIT = Path("/private/tmp/claude-501/-Users-marco-Desktop-IW-eftcamb-update-EFTCAMB/"
             "2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad/audit")
BUILD = AUDIT / "build" / "A0"
WORK = AUDIT / "work" / "t4a-param" / "repro" / Path(__file__).stem
LIB = "/Users/marco/Software/miniforge3/lib"

SRC = r"""
program repro
use precision
use IniObjects
use EFTCAMB_parametrizations_1D
implicit none
class(parametrized_function_1D), allocatable :: f
type(TIniFile) :: Ini
integer :: err, i
real(dl) :: xs(5) = [0.1_dl, 0.3_dl, 0.5_dl, 0.8_dl, 1.0_dl]
character(len=512) :: path
call get_command_argument(1, path)
err = 0
call Ini%Open(trim(path))
call allocate_parametrized_1D_function(f, 6, 'EFTOmega', 'Omega', err, 0)
call f%init_func_from_file(Ini, err)
call f%init_from_file(Ini, err)
do i = 1, size(xs)
    write(*, '(5es26.17)') xs(i), f%value(xs(i)), f%second_derivative(xs(i)), f%third_derivative(xs(i)), &
        f%fourth_derivative(xs(i))
end do
end program
"""


def main():
    a0, up, down = 1.0, [0.1, 0.2, -0.3], [0.4]  # N = 2, D = 1 (denominator 1 + 0.4 (a - a0))
    ini = (f"EFTOmega_Pade_order_N = {len(up) - 1}\nEFTOmega_Pade_order_D = {len(down)}\nEFTOmegaa0 = {a0}\n"
           + "".join(f"EFTOmegaN{i} = {v}\n" for i, v in enumerate(up))
           + "".join(f"EFTOmegaD{i + 1} = {v}\n" for i, v in enumerate(down)))
    WORK.mkdir(parents=True, exist_ok=True)
    (WORK / "r.f90").write_text(SRC)
    (WORK / "r.ini").write_text(ini)
    subprocess.run(["gfortran", "-O0", "-cpp", "-ffree-line-length-none", "-fopenmp",
                    f"-I{BUILD}/camb/Releaselib", f"-I{BUILD}/forutils/Releaselib", "r.f90", "-o", "r",
                    str(BUILD / "camb/Releaselib/camblib.so"), f"-Wl,-rpath,{LIB}"], cwd=WORK, check=True,
                   capture_output=True)
    out = subprocess.run([str(WORK / "r"), str(WORK / "r.ini")], capture_output=True, text=True, check=True).stdout
    a = sp.symbols("a")
    expr = sum(sp.Rational(str(c)) * (a - a0) ** i for i, c in enumerate(up)) / (
        1 + sum(sp.Rational(str(c)) * (a - a0) ** (i + 1) for i, c in enumerate(down)))
    print("Omega(a) = (0.1 + 0.2 (a-1) - 0.3 (a-1)^2) / (1 + 0.4 (a-1))")
    print(f"{'a':>5} | {'d3 obs':>14} {'d3 exact':>14} | {'d4 obs':>14} {'d4 exact':>14} {'rel err d4':>11}")
    for line in out.splitlines():
        x, v, d2, d3, d4 = map(float, line.split())
        e3 = float(sp.diff(expr, a, 3).subs(a, x))
        e4 = float(sp.diff(expr, a, 4).subs(a, x))
        print(f"{x:5.2f} | {d3:14.8f} {e3:14.8f} | {d4:14.8f} {e4:14.8f} {abs(d4 - e4) / abs(e4):11.3e}")


if __name__ == "__main__":
    main()
