"""Repro (t4a-param): outside the node range [x_1, x_N] the cubic spline (flag 9) and quintic spline (flag 10)
parametrizations return value = null_value (a constant) but ALSO first..fourth derivative = null_value, instead of 0
(04p016_spline_parametrizations_1D.f90:380/404/426/448, 04p019_spline5_parametrizations_1D.f90:430/481/536/592).
For a w(a) spline with null_value = -1 this gives w' = w'' = w''' = -1 wherever the nodes do not reach.

1) function level: w = -1 on nodes 0.01..0.9 with null_value = -1 (the function is exactly -1 everywhere)
2) model level: Pure EFT 2_PEFT_Omega_lin_1.ini (EFTwDE = 0, w = -1) vs the same w = -1 given as
   EFTwDE = 9 / 10 splines on nodes 0.01..0.9 with null_value = -1: expected identical spectra.

usage: /Users/marco/Software/miniforge3/bin/python spline-derivatives-outside-nodes.py [--no-model]
"""

import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "checks"))
import harness  # noqa: E402

AUDIT = harness.AUDIT_ROOT
BUILD = AUDIT / "build" / "A0"
WORK = AUDIT / "work" / "t4a-param" / "repro" / Path(__file__).stem
LIB = "/Users/marco/Software/miniforge3/lib"
NODES = [0.01, 0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 0.9]

SRC = r"""
program repro
use precision
use IniObjects
use EFTCAMB_parametrizations_1D
implicit none
class(parametrized_function_1D), allocatable :: f
type(TIniFile) :: Ini
integer :: err, i, n, flag
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
read(*, *) n
do i = 1, n
    read(*, *) x
    write(*, '(6es26.17)') x, f%value(x), f%first_derivative(x), f%second_derivative(x), f%third_derivative(x), &
        f%fourth_derivative(x)
end do
end program
"""


def spline_keys(flag: int) -> dict:
    d = {"EFTw_Spline_Pixels": len(NODES), "EFTw_null_value": -1.0}
    if flag == 10:
        d.update({"EFTw_Spline_Order": 5, "EFTw_Logx_Flag": "F", "EFTw_Logy_Flag": "F"})
    d.update({f"EFTwx{i + 1}": x for i, x in enumerate(NODES)})
    d.update({f"EFTwv{i + 1}": -1.0 for i in range(len(NODES))})
    return d


def function_level():
    WORK.mkdir(parents=True, exist_ok=True)
    (WORK / "r.f90").write_text(SRC)
    subprocess.run(["gfortran", "-O0", "-cpp", "-ffree-line-length-none", "-fopenmp",
                    f"-I{BUILD}/camb/Releaselib", f"-I{BUILD}/forutils/Releaselib", "r.f90", "-o", "r",
                    str(BUILD / "camb/Releaselib/camblib.so"), f"-Wl,-rpath,{LIB}"], cwd=WORK, check=True,
                   capture_output=True)
    xs = [1e-4, 0.005, 0.4, 0.95, 1.0]
    bad = 0
    for flag in (9, 10):
        (WORK / f"f{flag}.ini").write_text("\n".join(f"{k} = {v}" for k, v in spline_keys(flag).items()) + "\n")
        inp = f"{len(xs)}\n" + "".join(f"{x:.17e}\n" for x in xs)
        out = subprocess.run([str(WORK / "r"), str(WORK / f"f{flag}.ini"), str(flag)], input=inp, capture_output=True,
                             text=True, check=True).stdout
        print(f"flag {flag}: w = -1 on nodes {NODES[0]}..{NODES[-1]}, null_value = -1 (expected: value -1, all "
              "derivatives 0)")
        print(f"{'a':>8} {'value':>9} {'d1':>9} {'d2':>9} {'d3':>9} {'d4':>9}")
        for line in out.splitlines():
            r = list(map(float, line.split()))
            wrong = any(abs(d) > 1e-8 for d in r[2:])
            bad += wrong
            print(f"{r[0]:8.4g} " + " ".join(f"{v:9.4f}" for v in r[1:]) + ("   <-- WRONG" if wrong else ""))
    print(f"   -> {bad} of {2 * len(xs)} points with non-zero derivatives of a constant function\n")


def model_level():
    ini = harness.PARAMS_DIR / "2_PEFT_Omega_lin_1.ini"
    ref = WORK / "ref.npz"
    st = harness.run_model("A0", ini, ref, threads=2)
    print("reference EFTwDE = 0:", st["status"])
    for flag in (9, 10):
        out = WORK / f"spline{flag}.npz"
        ov = {"EFTwDE": str(flag), **{k: str(v) for k, v in spline_keys(flag).items()}}
        s = harness.run_model("A0", ini, out, overrides=ov, threads=2)
        if s["status"] != "ok":
            print(f"EFTwDE = {flag}: {s['status']}  {(s.get('log_tail') or [''])[-1]}")
            continue
        c = harness.compare(harness.load(out), harness.load(ref))
        eft = sorted(c["eft"].items(), key=lambda kv: -kv[1])[:4]
        print(f"EFTwDE = {flag} (w = -1 as spline): worst Cl rel diff vs EFTwDE=0 = {c['worst_cl']:.3e}, "
              f"worst P(k) = {c['worst_pk']:.3e}, sigma8 = {c['sigma8']:.3e}  (expected ~0)")
        print("    largest EFT-function differences:", ", ".join(f"{k} {v:.2e}" for k, v in eft))


if __name__ == "__main__":
    function_level()
    if "--no-model" not in sys.argv:
        model_level()
