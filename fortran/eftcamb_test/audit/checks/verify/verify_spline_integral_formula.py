"""Verifier (spline-integral-convention): formula-level attribution of the cubic (9) / quintic (10) spline integral().

Independent of the reporter's driver/data:
  1. own Fortran driver compiled against build/A0 (release lib) that prints value() and integral() for any ini;
  2. own numpy re-implementation of what the Fortran code computes (the "bug model"):
       cubic : yint_i = prod (x_i/x_{i+1})^(-3(1+b)) exp(-3(y_i-y_{i+1})), yint_N = 1, then a cubic spline (end 2nd
               derivative from the cubic through the first/last 4 points = spline3pars bctype 2) of yint, 0 outside
               (x_1, x_N)
       quintic (Logx): bcoef(:,i) = [N_i, c1..c6]; init  yint_i = yint_{i+1} (N_i/N_{i+1})^3 exp(-tmp), tmp cumulative
                       sum_j bcoef(j)*dN^j/j; integral yint_ip (N/N_ip)^3 exp(sum_j bcoef(j)*dN^j/j)
     and the correct EFTCAMB convention I(a) = a^2 exp(-3 int_1^a (1+w)/a' da') (scipy quad on the spline w).
Cases use non-degenerate data (CPL-like linear w, and a quintic node set ending at a=0.5 so that ln a_N != 0) so that
finite numbers can be matched digit by digit.
usage: /Users/marco/Software/miniforge3/bin/python verify_spline_integral_formula.py
"""

import re
import subprocess
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from scipy.interpolate import CubicSpline

AUDIT = Path("/private/tmp/claude-501/-Users-marco-Desktop-IW-eftcamb-update-EFTCAMB/"
             "2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad/audit")
BUILD = AUDIT / "build" / "A0"
WORK = AUDIT / "work" / "verify" / "spline_integral" / "formula"
LIB = "/Users/marco/Software/miniforge3/lib"

SRC = r"""
program vdrv
use precision
use IniObjects
use EFTCAMB_parametrizations_1D
implicit none
class(parametrized_function_1D), allocatable :: f
type(TIniFile) :: Ini
integer :: err, i, flag, n
real(dl) :: a
character(len=512) :: path, sflag
call get_command_argument(1, path)
call get_command_argument(2, sflag)
read(sflag, *) flag
err = 0
call Ini%Open(trim(path))
call allocate_parametrized_1D_function(f, flag, 'EFTw', 'w', err, 0)
call f%init_func_from_file(Ini, err)
call f%init_from_file(Ini, err)
open(10, file='a.txt')
read(10, *) n
do i = 1, n
    read(10, *) a
    write(*, '(3es26.17)') a, f%value(a), f%integral(a)
end do
end program
"""

A_PROBE = np.array([3e-6, 1e-4, 0.003, 0.05, 0.15, 0.3, 0.45, 0.55, 0.7, 0.85, 0.95, 1.0])


def build():
    WORK.mkdir(parents=True, exist_ok=True)
    (WORK / "vdrv.f90").write_text(SRC)
    subprocess.run(["gfortran", "-O0", "-cpp", "-ffree-line-length-none", "-fopenmp",
                    f"-I{BUILD}/camb/Releaselib", f"-I{BUILD}/forutils/Releaselib", "vdrv.f90", "-o", "vdrv",
                    str(BUILD / "camb/Releaselib/camblib.so"), f"-Wl,-rpath,{LIB}"], cwd=WORK, check=True)


def run(flag, xs, ws, extra="", a=A_PROBE):
    ini = (f"EFTw_Spline_Pixels = {len(xs)}\nEFTw_null_value = {float(ws[0])!r}\n" + extra
           + "".join(f"EFTwx{i + 1} = {float(x)!r}\nEFTwv{i + 1} = {float(w)!r}\n" for i, (x, w) in enumerate(zip(xs, ws))))
    (WORK / "v.ini").write_text(ini)
    (WORK / "a.txt").write_text(f"{len(a)}\n" + "".join(f"{float(x)!r}\n" for x in a))
    out = subprocess.run([str(WORK / "vdrv"), "v.ini", str(flag)], capture_output=True, text=True, cwd=WORK)
    # Fortran es format drops the 'E' for 3-digit exponents (e.g. 7.4+137)
    arr = np.array([[float(re.sub(r"(\d)([+-]\d{3})$", r"\1E\2", t)) for t in ln.split()]
                    for ln in out.stdout.splitlines()])
    return arr[:, 1], arr[:, 2]


def spline3_bc2(x, y):
    cl = np.polyfit(x[:4] - x[0], y[:4], 3)
    cr = np.polyfit(x[-4:] - x[-1], y[-4:], 3)
    return CubicSpline(x, y, bc_type=((2, 2 * cl[1]), (2, 2 * cr[1])))


def cubic_bug_model(xs, ys, a):
    n = len(xs)
    yint = np.ones(n)
    for i in range(n - 2, -1, -1):
        b = (xs[i] * ys[i + 1] - xs[i + 1] * ys[i]) / (xs[i] - xs[i + 1])
        yint[i] = yint[i + 1] * (xs[i] / xs[i + 1]) ** (-3 * (1 + b)) * np.exp(-3 * (ys[i] - ys[i + 1]))
    s = spline3_bc2(xs, yint)
    out = s(a)
    out[(a <= xs[0]) | (a >= xs[-1])] = 0.0
    return out


def quintic_bug_model_constw(lnx, w, a):
    """bcoef(:,i) = [N_i, w, 0, 0, 0, 0, 0] for constant data."""
    n = len(lnx)
    yint = np.ones(n)
    tmp = 0.0
    for i in range(n - 2, -1, -1):
        dx = lnx[i + 1] - lnx[i]
        tmp += lnx[i] * dx + w * dx ** 2 / 2
        yint[i] = yint[i + 1] * (lnx[i] / lnx[i + 1]) ** 3 * np.exp(-tmp)
    out = []
    for aa in a:
        x_ = np.log(aa)
        if x_ < lnx[0]:
            out.append(yint[0])
        elif x_ > lnx[-1]:
            out.append(1.0)
        else:
            ip = min(np.searchsorted(lnx, x_, side="right") - 1, n - 2)  # probes avoid the nodes (iix tie)
            dx = x_ - lnx[ip]
            t = lnx[ip] * dx + w * dx ** 2 / 2
            out.append(yint[ip] * (x_ / lnx[ip]) ** 3 * np.exp(t))
    return np.array(out)


def exact_I(wfun, a):
    return np.array([aa ** 2 * np.exp(-3 * quad(lambda t: (1 + wfun(t)) / t, 1.0, aa, limit=200,
                                                 points=None)[0]) for aa in a])


def show(title, a, cols):
    print(title)
    print(f"{'a':>9}" + "".join(f"{k:>16}" for k in cols))
    for i in range(len(a)):
        print(f"{a[i]:9.2e}" + "".join(f"{cols[k][i]:16.8g}" for k in cols))
    print()


def relmax(x, y, m):
    return float(np.max(np.abs(x[m] - y[m]) / np.abs(y[m])))


def main():
    build()
    res = {}
    # ---- case 1: cubic, CPL-like linear w(a) = -0.9 + 0.2 (1 - a) on 11 nodes 0.05 .. 1 ---------------------------
    w0, wa = -0.9, 0.2
    xs = np.linspace(0.05, 1.0, 11)
    ws = w0 + wa * (1 - xs)
    v, i_f = run(9, xs, ws)
    wcpl = lambda t: w0 + wa * (1 - t)  # noqa: E731
    i_ex = exact_I(wcpl, A_PROBE)
    i_bug = cubic_bug_model(xs, ws, A_PROBE)
    show("cubic (9), linear w = -0.9 + 0.2 (1-a), nodes 0.05..1", A_PROBE,
         {"w(a) fortran": v, "I fortran": i_f, "bug model": i_bug, "exact I": i_ex, "exact I/a^2": i_ex / A_PROBE ** 2})
    inside = (A_PROBE > xs[0]) & (A_PROBE < xs[-1])
    res["cubic_lin: fortran vs bug model (inside)"] = relmax(i_f, i_bug, inside)
    res["cubic_lin: fortran vs exact I (inside)"] = relmax(i_f, i_ex, inside)
    res["cubic_lin: fortran vs exact I/a^2 (inside)"] = relmax(i_f, i_ex / A_PROBE ** 2, inside)
    res["cubic_lin: fortran I(a=1)"] = float(i_f[-1])
    res["cubic_lin: fortran I(a<x1) max"] = float(np.max(np.abs(i_f[A_PROBE <= xs[0]])))
    res["cubic_lin: w(a) fortran vs exact (inside)"] = relmax(v, wcpl(A_PROBE), inside)

    # ---- case 2: quintic Logx, constant w = -0.9, nodes a = 1e-6 .. 0.5 (ln a_N = -0.69 != 0) ----------------------
    w = -0.9
    xq = np.array([1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 0.1, 0.3, 0.5])
    aq = np.array([3e-6, 2e-4, 0.003, 0.05, 0.15, 0.35, 0.45, 0.55, 0.7, 0.85, 0.95, 1.0])  # off-node probes
    v, i_f = run(10, xq, np.full(len(xq), w), "EFTw_Logx_Flag = T\n", a=aq)
    i_bug = quintic_bug_model_constw(np.log(xq), w, aq)
    i_ex = aq ** (-1 - 3 * w)
    show("quintic (10) Logx, constant w = -0.9, nodes 1e-6..0.5", aq,
         {"w(a) fortran": v, "I fortran": i_f, "bug model": i_bug, "exact a^1.7": i_ex})
    m = np.isfinite(i_bug) & (i_bug != 0)
    res["quintic_const: fortran vs bug model (finite pts)"] = relmax(i_f, i_bug, m)
    res["quintic_const: fortran vs exact (a<=0.5)"] = relmax(i_f, i_ex, aq <= 0.5)
    # discontinuity of the coded formula at the node a = 1e-4 (left vs right interval)
    lnq = np.log(xq)
    left = quintic_bug_model_constw(lnq, w, np.array([1e-4 * (1 - 1e-12)]))[0]
    yint_r = quintic_bug_model_constw(lnq, w, np.array([1e-4 * (1 + 1e-12)]))[0]
    res["quintic_const: coded I jump at node a=1e-4 (ratio)"] = float(left / yint_r)

    # ---- case 3: quintic Logx, the reporter's node layout but our own w(a): linear-in-ln a data -------------------
    xq2 = np.array([10.0 ** e for e in range(-8, 1)])
    v, i_f = run(10, xq2, -0.95 + 0.02 * np.log(xq2) / np.log(10), "EFTw_Logx_Flag = T\n")
    show("quintic (10) Logx, nodes 1e-8..1 (last node a = 1)", A_PROBE, {"w(a) fortran": v, "I fortran": i_f})
    res["quintic_lastnode1: finite I count"] = int(np.sum(np.isfinite(i_f)))

    for k, val in res.items():
        print(f"{k:55s} {val:.4g}")


if __name__ == "__main__":
    main()
