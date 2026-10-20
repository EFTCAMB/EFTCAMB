"""Evaluate the spline5_parametrization_1D (allocator flag 10, order 5) directly for given nodes, A0 vs poly5-fix build.

Used to check what the pure-EFT Omega spline (log15_a2 case of verify_spline5_cl_impact.py) actually returns and to
separate the poly5 typo from the spline's own interpolation error.
usage: python verify_spline5_param_eval.py
"""

import math
import subprocess
from pathlib import Path

AUDIT = Path("/private/tmp/claude-501/-Users-marco-Desktop-IW-eftcamb-update-EFTCAMB/"
             "2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad/audit")
BUILDS = {"A0": AUDIT / "build" / "A0", "FIX": AUDIT / "work" / "verify" / "spline5_poly5" / "fix_poly5_build"}
WORK = AUDIT / "work" / "verify" / "spline5_poly5" / "param_driver"
LIB = "/Users/marco/Software/miniforge3/lib"
SRC = r"""
program p
use precision
use IniObjects
use EFTCAMB_parametrizations_1D
implicit none
class(parametrized_function_1D), allocatable :: f
type(TIniFile) :: Ini
integer :: err, i, n
real(dl) :: x
character(len=512) :: path
call get_command_argument(1, path)
err = 0
call Ini%Open(trim(path))
call allocate_parametrized_1D_function(f, 10, 'EFTOmega', 'Omega', err, 0)
call f%init_func_from_file(Ini, err)
call f%init_from_file(Ini, err)
read(*, *) n
do i = 1, n
    read(*, *) x
    write(*, '(3es26.17)') x, f%value(x), f%first_derivative(x)
end do
end program
"""


def evaluate(tag, nodes, vals, logx, xs):
    d = WORK / tag
    d.mkdir(parents=True, exist_ok=True)
    b = BUILDS[tag]
    (d / "p.f90").write_text(SRC)
    ini = (f"EFTOmega_Spline_Order = 5\nEFTOmega_Logx_Flag = {'T' if logx else 'F'}\n"
           f"EFTOmega_Spline_Pixels = {len(nodes)}\n"
           + "".join(f"EFTOmegax{i + 1} = {float(x)!r}\nEFTOmegav{i + 1} = {float(v)!r}\n" for i, (x, v) in enumerate(zip(nodes, vals))))
    (d / "p.ini").write_text(ini)
    subprocess.run(["gfortran", "-O0", "-cpp", "-ffree-line-length-none", "-fopenmp", f"-I{b}/camb/Releaselib",
                    f"-I{b}/forutils/Releaselib", "p.f90", "-o", "p", str(b / "camb/Releaselib/camblib.so"),
                    f"-Wl,-rpath,{LIB}"], cwd=d, check=True, capture_output=True)
    inp = f"{len(xs)}\n" + "".join(f"{x:.17e}\n" for x in xs)
    out = subprocess.run([str(d / "p"), str(d / "p.ini")], input=inp, capture_output=True, text=True, check=True).stdout
    return [tuple(map(float, l.split())) for l in out.splitlines() if l.strip()]


def main():
    nodes = [10.0 ** (-3 + 3 * i / 14) for i in range(15)]
    vals = [a**2 for a in nodes]
    xs = [10 ** (-3 + 3 * k / 280) for k in range(281)]
    res = {t: evaluate(t, nodes, vals, True, xs) for t in BUILDS}
    print("log15_a2: Omega = a^2 on 15 nodes 1e-3..1 (Logx=T)")
    for t, r in res.items():
        ev = [abs(v / x**2 - 1) for x, v, _ in r]
        ed = [abs(dv / (2 * x) - 1) for x, _, dv in r]
        iv, idd = max(range(len(ev)), key=ev.__getitem__), max(range(len(ed)), key=ed.__getitem__)
        print(f"  {t}: max |Omega/a^2-1| = {ev[iv]:.3e} at a={r[iv][0]:.4g};  max |Omega'/(2a)-1| = {ed[idd]:.3e} "
              f"at a={r[idd][0]:.4g};  over a>=0.01: {max(e for e, (x, _, _) in zip(ev, r) if x >= 0.01):.3e}, "
              f"{max(e for e, (x, _, _) in zip(ed, r) if x >= 0.01):.3e}")
    diff = max(abs(a[1] - b[1]) / a[0] ** 2 for a, b in zip(res["A0"], res["FIX"]))
    print(f"  A0 vs FIX value max |diff|/a^2 = {diff:.3e}")


if __name__ == "__main__":
    main()
