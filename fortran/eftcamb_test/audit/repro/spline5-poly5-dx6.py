"""Repro (t4a-param): quintic spline parametrization (allocator flag 10, EFT*_Spline_Order = 5) evaluates its value with
poly5 = c1 + c2 dx + c3 dx^2 + c4 dx^3 + c5 dx^4 + c6 dx^6   (02_bundled_algorithms/02p12_splines5.f90:391)
while dpoly5..d4poly5 differentiate c6 dx^5.  value() therefore does not interpolate the nodes, jumps at every
interior node and is inconsistent with first_derivative().

Checks (expected values come from the spline's own definition, no external spline needed):
  * value at a node == node value (interpolation property)
  * value continuous across interior nodes (left vs right limit)
  * first_derivative == centred finite difference of value

usage: /Users/marco/Software/miniforge3/bin/python spline5-poly5-dx6.py
"""

import math
import subprocess
from pathlib import Path

AUDIT = Path("/private/tmp/claude-501/-Users-marco-IW-eftcamb-update-EFTCAMB/"
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
integer :: err, i, n
real(dl) :: x
character(len=512) :: path
call get_command_argument(1, path)
err = 0
call Ini%Open(trim(path))
call allocate_parametrized_1D_function(f, 10, 'EFTw', 'w', err, 0)
call f%init_func_from_file(Ini, err)
call f%init_from_file(Ini, err)
read(*, *) n
do i = 1, n
    read(*, *) x
    write(*, '(3es26.17)') x, f%value(x), f%first_derivative(x)
end do
end program
"""


def run(ini: str, xs: list[float]) -> list[tuple[float, float, float]]:
    WORK.mkdir(parents=True, exist_ok=True)
    (WORK / "r.f90").write_text(SRC)
    (WORK / "r.ini").write_text(ini)
    subprocess.run(["gfortran", "-O0", "-cpp", "-ffree-line-length-none", "-fopenmp",
                    f"-I{BUILD}/camb/Releaselib", f"-I{BUILD}/forutils/Releaselib", "r.f90", "-o", "r",
                    str(BUILD / "camb/Releaselib/camblib.so"), f"-Wl,-rpath,{LIB}"], cwd=WORK, check=True,
                   capture_output=True)
    inp = f"{len(xs)}\n" + "".join(f"{x:.17e}\n" for x in xs)
    out = subprocess.run([str(WORK / "r"), str(WORK / "r.ini")], input=inp, capture_output=True, text=True,
                         check=True).stdout
    return [tuple(map(float, line.split())) for line in out.splitlines() if line.strip()]


def main():
    nodes = [1e-3, 0.01, 0.05, 0.1, 0.25, 0.4, 0.6, 0.8, 1.0]
    for logx in (False, True):
        vals = [0.2 * math.sin(4 * x) + 0.1 for x in nodes]
        ini = (f"EFTw_Spline_Order = 5\nEFTw_Logx_Flag = {'T' if logx else 'F'}\nEFTw_Spline_Pixels = {len(nodes)}\n"
               + "".join(f"EFTwx{i + 1} = {x}\nEFTwv{i + 1} = {v}\n" for i, (x, v) in enumerate(zip(nodes, vals))))
        eps = 1e-9
        xs = []
        for x in nodes[1:-1]:
            xs += [x * (1 - eps), x, x * (1 + eps)]
        mids = [math.sqrt(a * b) if logx else 0.5 * (a + b) for a, b in zip(nodes[:-1], nodes[1:])]
        for m in mids:
            h = 1e-6 * m
            xs += [m - h, m, m + h]
        r = run(ini, xs)
        print(f"Logx_Flag = {logx}: f = 0.2 sin(4a) + 0.1 sampled at {len(nodes)} nodes")
        print(f"{'node':>8} {'value(node)':>14} {'node value':>12} {'jump L/R':>10}")
        worst_i = worst_j = 0.0
        for k, x in enumerate(nodes[1:-1]):
            lo, mid, hi = r[3 * k], r[3 * k + 1], r[3 * k + 2]
            v = vals[k + 1]
            worst_i = max(worst_i, abs(mid[1] - v) / abs(v))
            worst_j = max(worst_j, abs(hi[1] - lo[1]) / abs(v))
            print(f"{x:8.3f} {mid[1]:14.9f} {v:12.9f} {abs(hi[1] - lo[1]) / abs(v):10.2e}")
        off = 3 * len(nodes[1:-1])
        worst_d = 0.0
        for k, m in enumerate(mids):
            a, b, c = r[off + 3 * k], r[off + 3 * k + 1], r[off + 3 * k + 2]
            fd = (c[1] - a[1]) / (c[0] - a[0])
            worst_d = max(worst_d, abs(b[2] - fd) / max(abs(fd), 1e-3))
        print(f"   observed: worst |value(node)-v|/|v| = {worst_i:.2e}, worst value jump at nodes = {worst_j:.2e}, "
              f"worst |d1 - FD(value)|/|FD| = {worst_d:.2e}")
        print("   expected: all three ~1e-9 or smaller (interpolating, C4-continuous quintic spline)\n")


if __name__ == "__main__":
    main()
