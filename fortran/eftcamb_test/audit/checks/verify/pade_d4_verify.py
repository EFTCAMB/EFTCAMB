"""Adversarial verification of finding pade-fourth-derivative (04p010_pade_series_1D.f90:419).

Independent of the t4a-param repro: different Pade orders (N=5, D=4, a0=0.6, so every Taylor derivative of
numerator and denominator is non-zero), different evaluation points, and three independent references:
  (1) exact sympy derivative of the rational function (exact rationals);
  (2) Richardson-extrapolated central finite difference of the code's own third_derivative (no sympy);
  (3) Python transcription of the Fortran expression as written (no u') and with u' inserted.
Controls: D=0 (pure polynomial, the buggy term vanishes) and N=0 (u'=0, error should be the full 4*T term).

usage: /Users/marco/Software/miniforge3/bin/python pade_d4_verify.py
"""

import json
import subprocess
from pathlib import Path

import sympy as sp

AUDIT = Path("/private/tmp/claude-501/-Users-marco-Desktop-IW-eftcamb-update-EFTCAMB/"
             "2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad/audit")
BUILD = AUDIT / "build" / "A0"
WORK = AUDIT / "work" / "verify" / "pade_d4"
LIB = "/Users/marco/Software/miniforge3/lib"

SRC = r"""
program pade_verify
use precision
use IniObjects
use EFTCAMB_parametrizations_1D
implicit none
class(parametrized_function_1D), allocatable :: f
type(TIniFile) :: Ini
integer :: err, i, nx
real(dl) :: x, h
character(len=512) :: path, buf
call get_command_argument(1, path)
call get_command_argument(2, buf)
read(buf, *) nx
err = 0
call Ini%Open(trim(path))
call allocate_parametrized_1D_function(f, 6, 'PTest', 'P', err, 0)
call f%init_func_from_file(Ini, err)
call f%init_from_file(Ini, err)
if (err /= 0) stop 'init error'
do i = 1, nx
    call get_command_argument(2 + i, buf)
    read(buf, *) x
    h = 2.e-3_dl
    write(*, '(12es26.17)') x, f%value(x), f%first_derivative(x), f%second_derivative(x), &
        f%third_derivative(x), f%fourth_derivative(x), &
        f%third_derivative(x + h), f%third_derivative(x - h), &
        f%third_derivative(x + h/2), f%third_derivative(x - h/2), h, 0._dl
end do
end program
"""


def build():
    WORK.mkdir(parents=True, exist_ok=True)
    (WORK / "d.f90").write_text(SRC)
    subprocess.run(["gfortran", "-O0", "-cpp", "-ffree-line-length-none", "-fopenmp",
                    f"-I{BUILD}/camb/Releaselib", f"-I{BUILD}/forutils/Releaselib", "d.f90", "-o", "d",
                    str(BUILD / "camb/Releaselib/camblib.so"), f"-Wl,-rpath,{LIB}"], cwd=WORK, check=True,
                   capture_output=True)


def run_case(name, a0, up, down, xs):
    ini = (f"PTest_Pade_order_N = {len(up) - 1}\nPTest_Pade_order_D = {len(down)}\nPTesta0 = {a0}\n"
           + "".join(f"PTestN{i} = {v}\n" for i, v in enumerate(up))
           + "".join(f"PTestD{i + 1} = {v}\n" for i, v in enumerate(down)))
    p = WORK / f"{name}.ini"
    p.write_text(ini)
    out = subprocess.run([str(WORK / "d"), str(p), str(len(xs))] + [repr(x) for x in xs], capture_output=True,
                         text=True, check=True).stdout
    a = sp.symbols("a")
    A0 = sp.Rational(str(a0))
    U = sum(sp.Rational(str(c)) * (a - A0) ** i for i, c in enumerate(up))
    G = 1 + sum(sp.Rational(str(c)) * (a - A0) ** (i + 1) for i, c in enumerate(down))
    F = U / G
    d = [sp.diff(F, a, k) for k in range(5)]
    dU = [sp.diff(U, a, k) for k in range(5)]
    dG = [sp.diff(G, a, k) for k in range(5)]
    rows = []
    for line in out.splitlines():
        v = list(map(float, line.split()))
        x, val, d1, d2, d3, d4, t3p, t3m, t3p2, t3m2, h = v[:11]
        xr = sp.Rational(repr(x))
        ex = [float(dk.subs(a, xr)) for dk in d]
        # Richardson-extrapolated central FD of the code's own d3 (error O(h^4))
        fd_h = (t3p - t3m) / (2 * h)
        fd_h2 = (t3p2 - t3m2) / h
        fd = (4 * fd_h2 - fd_h) / 3
        # transcription of the Fortran expression (line 417-421), with and without u'
        u, up1, up2, up3, up4 = (float(q.subs(a, xr)) for q in dU)
        g, g1, g2, g3, g4 = (float(q.subs(a, xr)) for q in dG)
        T = -g3 / g**2 - 6 * g1**3 / g**4 + 6 * g1 * g2 / g**3
        rest = (up4 / g - 4 * (up3 * g1 / g**2) + 6 * up2 * (2 * g1**2 / g**3 - g2 / g**2)
                + u * (-g4 / g**2 + 6 * g2**2 / g**3 + 24 * g1**4 / g**5 + 8 * g3 * g1 / g**3
                       - 36 * g2 * g1**2 / g**4))
        as_written = rest + 4 * T
        fixed = rest + 4 * up1 * T
        rows.append(dict(x=x, d0_rel=abs(val - ex[0]) / abs(ex[0]), d1_rel=abs(d1 - ex[1]) / abs(ex[1]),
                         d2_rel=abs(d2 - ex[2]) / abs(ex[2]), d3_rel=abs(d3 - ex[3]) / abs(ex[3]),
                         d4_code=d4, d4_exact=ex[4], d4_fd=fd, d4_rel_code=abs(d4 - ex[4]) / abs(ex[4]),
                         fd_vs_exact_rel=abs(fd - ex[4]) / abs(ex[4]),
                         transcribed_as_written_vs_code=abs(as_written - d4) / abs(d4),
                         transcribed_fixed_vs_exact=abs(fixed - ex[4]) / abs(ex[4]),
                         predicted_error=4 * (1 - up1) * T, observed_error=d4 - ex[4], uprime=up1))
    return rows


def main():
    build()
    xs = [0.05, 0.2, 0.45, 0.6, 0.85, 1.0, 1.15]
    cases = {
        "N5D4": (0.6, [0.05, -0.3, 0.7, -0.4, 0.25, -0.1], [0.3, -0.2, 0.15, -0.05]),
        "N0D2_uprime0": (1.0, [0.2], [0.5, 0.1]),
        "N5D0_control": (0.6, [0.05, -0.3, 0.7, -0.4, 0.25, -0.1], []),
    }
    res = {}
    for name, (a0, up, down) in cases.items():
        rows = run_case(name, a0, up, down, xs)
        res[name] = rows
        print(f"== {name}: a0={a0} N={up} D={down}")
        print(f"{'x':>5} {'d3 rel':>9} {'d4 code':>14} {'d4 exact':>14} {'d4 FD(d3)':>14} {'rel err':>9} "
              f"{'FD-ex':>9} {'tr(as-is)':>9} {'tr(fix)':>9} {'pred err':>12} {'obs err':>12}")
        for r in rows:
            print(f"{r['x']:5.2f} {r['d3_rel']:9.1e} {r['d4_code']:14.6f} {r['d4_exact']:14.6f} {r['d4_fd']:14.6f} "
                  f"{r['d4_rel_code']:9.2e} {r['fd_vs_exact_rel']:9.1e} {r['transcribed_as_written_vs_code']:9.1e} "
                  f"{r['transcribed_fixed_vs_exact']:9.1e} {r['predicted_error']:12.6f} {r['observed_error']:12.6f}")
    (WORK / "results.json").write_text(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
