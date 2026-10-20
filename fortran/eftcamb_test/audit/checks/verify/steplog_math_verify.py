"""Adversarial verification (steplog-d4-and-integral): independent math check of the steplog parametrization.

Compiles a small Fortran driver against the frozen A0 build, evaluates value, d1..d4 and integral of the
steplog function (allocator flag 8) for several parameter sets, and compares against mpmath (50 digits):
  - derivatives: mpmath numerical differentiation of the code's own value() closed form
    f(a) = 0.5 (v2-v1) L / sqrt(1+L^2) + (v1+v2)/2, L = ln(a/aT)/delta
  - integral: EFTCAMB convention I(a) = a^2 rho_DE(a)/rho_DE(1) = a^2 exp(3 int_{ln a}^0 (1+w) dN), checked to
    agree with the constant (flag 1) parametrization (I = a^(-1-3w)), evaluated with mpmath.quad.
Also evaluates the proposed fix formulas (d4, integral) in mpmath and the control parametrizations' integral(1).
"""

import json
import subprocess
import sys
from pathlib import Path

import mpmath as mp

mp.mp.dps = 50

AUDIT = Path("/private/tmp/claude-501/-Users-marco-Desktop-IW-eftcamb-update-EFTCAMB/"
             "2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad/audit")
_ARG = sys.argv[1] if len(sys.argv) > 1 else "A0"
BUILD = Path(_ARG) if "/" in _ARG else AUDIT / "build" / _ARG
WORK = AUDIT / "work" / "verify" / "steplog_phys" / ("math_" + BUILD.name)
LIB = "/Users/marco/Software/miniforge3/lib"

SRC = r"""
program drv
use precision
use IniObjects
use EFTCAMB_parametrizations_1D
implicit none
class(parametrized_function_1D), allocatable :: f
type(TIniFile) :: Ini
integer :: err, i, flag, n
real(dl) :: x
character(len=512) :: path, sflag
call get_command_argument(1, path)
call get_command_argument(2, sflag)
read(sflag, *) flag
err = 0
call Ini%Open(trim(path))
call allocate_parametrized_1D_function(f, flag, 'EFTw', 'w', err, 0)
call f%init_func_from_file(Ini, err)
call f%init_from_file(Ini, err)
n = Ini%Read_Int('npts', 0)
do i = 1, n
    x = Ini%Read_Double('x' // trim(adjustl(itoa(i))), 1._dl)
    write(*, '(7es26.17)') x, f%value(x), f%first_derivative(x), f%second_derivative(x), &
        f%third_derivative(x), f%fourth_derivative(x), f%integral(x)
end do
contains
function itoa(k) result(s)
integer, intent(in) :: k
character(len=16) :: s
write(s, '(I0)') k
end function
end program
"""

XS = [1e-4, 0.01, 0.05, 0.2, 0.3, 0.385, 0.5, 0.642, 0.7, 0.9, 0.999999, 1.0]
PSETS = {
    "finding": (-1.0, -0.8, 0.3, 0.5),
    "peft_omega": (0.0, 1.5, 0.5, 0.5),
    "late_noNaN": (0.2, -0.1, 0.8, 0.5),   # aT exp(delta^2) = 1.027 > 1: no NaN on [0, 1]
    "aT_gt_1": (-1.2, -0.9, 1.5, 1.0),
    "narrow": (-1.0, -0.5, 0.6, 0.1),
}


def ref_value(a, v1, v2, at, dl):
    L = (mp.log(a) - mp.log(at)) / dl
    return (v2 - v1) / 2 * L / mp.sqrt(1 + L**2) + (v1 + v2) / 2


def ref_integral(a, v1, v2, at, dl):
    s = mp.quad(lambda n: 1 + ref_value(mp.e**n, v1, v2, at, dl), [mp.log(a), 0])
    return a**2 * mp.e**(3 * s)


def fix_d4(a, v1, v2, at, dl):
    t = (mp.log(a) - mp.log(at)) / dl
    return (v2 - v1) / 2 / a**4 * (15 * t * (3 - 4 * t**2) / (1 + t**2)**mp.mpf(4.5) / dl**4
                                  - 6 * (12 * t**2 - 3) / (1 + t**2)**mp.mpf(3.5) / dl**3
                                  - 33 * t / (1 + t**2)**mp.mpf(2.5) / dl**2
                                  - 6 / (1 + t**2)**mp.mpf(1.5) / dl)


def fix_integral(a, v1, v2, at, dl):
    return a**(-1 - mp.mpf(1.5) * (v1 + v2)) * mp.e**(-mp.mpf(1.5) * (v2 - v1) * (
        mp.sqrt(dl**2 + (mp.log(a) - mp.log(at))**2) - mp.sqrt(dl**2 + mp.log(at)**2)))


def rel(o, e):
    if o != o:
        return float("nan")
    return abs(o - float(e)) / max(abs(float(e)), 1e-300)


def run_driver(flag, keys, xs):
    tag = f"f{flag}_" + "_".join(f"{k}{v}" for k, v in keys.items())
    ini = WORK / f"{abs(hash(tag))}.ini"
    lines = [f"{k} = {v}" for k, v in keys.items()] + [f"npts = {len(xs)}"]
    lines += [f"x{i + 1} = {x!r}" for i, x in enumerate(xs)]
    ini.write_text("\n".join(lines) + "\n")
    out = subprocess.run([str(WORK / "drv"), str(ini), str(flag)], capture_output=True, text=True, check=True).stdout
    return [list(map(float, ln.split())) for ln in out.splitlines() if ln.strip()]


def main():
    WORK.mkdir(parents=True, exist_ok=True)
    (WORK / "drv.f90").write_text(SRC)
    subprocess.run(["gfortran", "-O0", "-cpp", "-ffree-line-length-none", "-fopenmp",
                    f"-I{BUILD}/camb/Releaselib", f"-I{BUILD}/forutils/Releaselib", "drv.f90", "-o", "drv",
                    str(BUILD / "camb/Releaselib/camblib.so"), f"-Wl,-rpath,{LIB}"], cwd=WORK, check=True)
    summary = {}
    for name, (v1, v2, at, dl) in PSETS.items():
        rows = run_driver(8, {"EFTw_v1": v1, "EFTw_v2": v2, "EFTw_at": at, "EFTw_delta": dl}, XS)
        mv = [mp.mpf(str(t)) for t in (v1, v2, at, dl)]
        print(f"\n== {name}: v1={v1} v2={v2} aT={at} delta={dl}; NaN onset predicted a > aT exp(delta^2) = "
              f"{at * float(mp.e**(dl**2)):.4f}")
        print(f"{'a':>9} | {'max rel d0..d3':>14} | {'d4 code':>13} {'d4 ref':>13} {'d4 fix':>13} | "
              f"{'I code':>13} {'I ref':>13} {'I fix':>13}")
        worst = {"d0-3": 0.0, "d4": 0.0, "d4_fix": 0.0, "I": 0.0, "I_fix": 0.0, "d4_nan": 0}
        for r in rows:
            x = mp.mpf(repr(r[0]))
            f = lambda a: ref_value(a, *mv)  # noqa: E731
            refs = [f(x)] + [mp.diff(f, x, n) for n in (1, 2, 3, 4)]
            e_low = max(rel(r[1 + n], refs[n]) if abs(refs[n]) > 1e-12 else abs(r[1 + n] - float(refs[n]))
                        for n in range(4))
            ri = ref_integral(x, *mv)
            f4, fi = fix_d4(x, *mv), fix_integral(x, *mv)
            worst["d0-3"] = max(worst["d0-3"], e_low)
            e4 = rel(r[5], refs[4])
            if e4 != e4:
                worst["d4_nan"] += 1
            else:
                worst["d4"] = max(worst["d4"], e4)
            worst["d4_fix"] = max(worst["d4_fix"], float(abs(f4 - refs[4]) / abs(refs[4])))
            worst["I"] = max(worst["I"], rel(r[6], ri))
            worst["I_fix"] = max(worst["I_fix"], float(abs(fi - ri) / ri))
            print(f"{r[0]:9.6g} | {e_low:14.3e} | {r[5]:13.6g} {float(refs[4]):13.6g} {float(f4):13.6g} | "
                  f"{r[6]:13.9f} {float(ri):13.9f} {float(fi):13.9f}")
        summary[name] = worst
        print("worst:", json.dumps(worst))
    # control: integral convention of other parametrizations at a = 1 and a = 0.5
    ctrl = run_driver(1, {"EFTw0": -0.9}, [0.5, 1.0])
    print(f"\ncontrol constant w=-0.9 (flag 1): I(0.5) = {ctrl[0][6]:.12f} (a^(-1-3w) = {0.5**(-1 + 2.7):.12f}),"
          f" I(1) = {ctrl[1][6]:.15f}")
    summary["control_constant_I1"] = ctrl[1][6]
    (WORK / "summary.json").write_text(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main()
