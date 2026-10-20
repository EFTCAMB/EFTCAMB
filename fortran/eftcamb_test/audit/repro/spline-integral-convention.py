"""Repro (t4a-param): integral() of the cubic spline (flag 9) and quintic spline (flag 10) parametrizations does not
follow the EFTCAMB convention used by every other w(a) parametrization and by its consumers
(designer f(R): grhov_t = grhov*integral(a), EFT_E_gfun = -(log(integral(a)) - 2 ln a)/3):

    I(a) = a^2 exp(-3 int_1^a (1 + w(a'))/a' da')        (constant w: a^(-1-3w); I(1) = 1)

  * cubic spline: returns a spline of exp(-3 int_{x_N}^{x} (1+w_lin)/x dx) of the piecewise-LINEAR w, without the a^2
    factor, normalised at the last node instead of a = 1, 0 outside the nodes and NaN when x1 = 0 (default nodes)
  * quintic spline (Logx_Flag = T required): sums bcoef(1:7) (bcoef(1) is the node position), uses (ln a1/ln a2)^3
    instead of (a2/a1)^3, never resets the running sum, misses the factor 3; -inf/Inf for nodes ending at a = 1
Data: constant w = -0.9 on all nodes (null_value = -0.9), so the exact answer is a^(-1-3w) = a^1.7.

usage: /Users/marco/Software/miniforge3/bin/python spline-integral-convention.py
"""

import subprocess
from pathlib import Path

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
integer :: err, i, flag
real(dl) :: xs(6) = [1.e-6_dl, 0.01_dl, 0.1_dl, 0.5_dl, 0.9_dl, 1._dl]
character(len=512) :: path, sflag
call get_command_argument(1, path)
call get_command_argument(2, sflag)
read(sflag, *) flag
err = 0
call Ini%Open(trim(path))
call allocate_parametrized_1D_function(f, flag, 'EFTw', 'w', err, 0)
call f%init_func_from_file(Ini, err)
call f%init_from_file(Ini, err)
do i = 1, size(xs)
    write(*, '(3es26.17)') xs(i), f%value(xs(i)), f%integral(xs(i))
end do
end program
"""


def main():
    WORK.mkdir(parents=True, exist_ok=True)
    (WORK / "r.f90").write_text(SRC)
    subprocess.run(["gfortran", "-O0", "-cpp", "-ffree-line-length-none", "-fopenmp",
                    f"-I{BUILD}/camb/Releaselib", f"-I{BUILD}/forutils/Releaselib", "r.f90", "-o", "r",
                    str(BUILD / "camb/Releaselib/camblib.so"), f"-Wl,-rpath,{LIB}"], cwd=WORK, check=True,
                   capture_output=True)
    cases = {
        "cubic spline, default-like nodes x=0..1": (9, [0.0, 0.2, 0.4, 0.6, 0.8, 1.0], ""),
        "cubic spline, nodes x=1e-8,0.2..1": (9, [1e-8, 0.2, 0.4, 0.6, 0.8, 1.0], ""),
        "quintic spline Logx, nodes 1e-8..1": (10, [10.0 ** e for e in range(-8, 1)], "EFTw_Logx_Flag = T\n"),
    }
    w = -0.9
    for name, (flag, xs, extra) in cases.items():
        ini = (f"EFTw_Spline_Pixels = {len(xs)}\nEFTw_null_value = {w}\n" + extra
               + "".join(f"EFTwx{i + 1} = {x}\nEFTwv{i + 1} = {w}\n" for i, x in enumerate(xs)))
        (WORK / "r.ini").write_text(ini)
        out = subprocess.run([str(WORK / "r"), str(WORK / "r.ini"), str(flag)], capture_output=True, text=True)
        print(f"{name} (w = {w}):")
        print(f"{'a':>8} {'value':>8} {'integral obs':>16} {'expected a^1.7':>16}")
        for line in out.stdout.splitlines():
            a, v, i_obs = map(float, line.split())
            print(f"{a:8.1e} {v:8.3f} {i_obs:16.8g} {a ** (-1 - 3 * w):16.8g}")
        if out.returncode:
            print("   program stopped:", out.stdout[-200:], out.stderr[-200:])
        print()


if __name__ == "__main__":
    main()
