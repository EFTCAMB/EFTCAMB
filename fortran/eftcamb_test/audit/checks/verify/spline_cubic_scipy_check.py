"""Verifier: compare EFTCAMB spline parametrization (flag 9) with scipy natural CubicSpline on curved data."""

import subprocess
from pathlib import Path

import numpy as np
from scipy.interpolate import CubicSpline

AUDIT = Path("/private/tmp/claude-501/-Users-marco-Desktop-IW-eftcamb-update-EFTCAMB/"
             "2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad/audit")
BUILD = AUDIT / "build" / "A0"
WORK = AUDIT / "work" / "verify" / "spline-cubic" / "scipy"
LIB = "/Users/marco/Software/miniforge3/lib"
SRC = r"""
program chk
use precision
use IniObjects
use EFTCAMB_parametrizations_1D
implicit none
class(parametrized_function_1D), allocatable :: f
type(TIniFile) :: Ini
integer :: err, i
real(dl) :: x
err = 0
call Ini%Open('c.ini')
call allocate_parametrized_1D_function(f, 9, 'EFTw', 'w', err, 0)
call f%init_func_from_file(Ini, err)
call f%init_from_file(Ini, err)
do i = 0, 400
    x = 0.1_dl + 0.9_dl*real(i, dl)/400._dl
    write(*, '(5es25.16)') x, f%value(x), f%first_derivative(x), f%second_derivative(x), f%third_derivative(x)
end do
end program
"""
xs = np.array([0.1, 0.25, 0.4, 0.55, 0.7, 0.85, 1.0])
vs = np.sin(3 * xs) + xs**2
WORK.mkdir(parents=True, exist_ok=True)
(WORK / "c.f90").write_text(SRC)
(WORK / "c.ini").write_text(f"EFTw_Spline_Pixels = {len(xs)}\nEFTw_null_value = 0\n" + "".join(
    f"EFTwx{i + 1} = {float(x)!r}\nEFTwv{i + 1} = {float(v)!r}\n" for i, (x, v) in enumerate(zip(xs, vs))))
subprocess.run(["gfortran", "-O0", "-cpp", "-ffree-line-length-none", "-fopenmp", f"-I{BUILD}/camb/Releaselib",
                f"-I{BUILD}/forutils/Releaselib", "c.f90", "-o", "c", str(BUILD / "camb/Releaselib/camblib.so"),
                f"-Wl,-rpath,{LIB}"], cwd=WORK, check=True)
out = np.loadtxt(subprocess.run(["./c"], cwd=WORK, capture_output=True, text=True, check=True).stdout.splitlines())
x, v, d1, d2, d3 = out.T
cs = CubicSpline(xs, vs, bc_type="natural")
for name, obs, k in (("value", v, 0), ("d1", d1, 1), ("d2", d2, 2), ("d3", d3, 3)):
    ref = cs(x, k)
    err = np.abs(obs - ref)
    i = np.argmax(err)
    print(f"{name}: max|EFTCAMB - scipy natural spline| = {err.max():.3e} at x={x[i]:.4f} "
          f"(obs {obs[i]:.5f}, ref {ref[i]:.5f}); interior (0.1<x<1) max = {err[1:-1].max():.3e}")
print("end nodes: value(0.1)=%.5f (exp %.5f), value(1.0)=%.5f (exp %.5f)" % (v[0], vs[0], v[-1], vs[-1]))
