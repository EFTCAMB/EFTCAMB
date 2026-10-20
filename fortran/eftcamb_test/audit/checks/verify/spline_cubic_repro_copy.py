"""Repro (t4a-param): cubic spline parametrization (allocator flag 9) is wrong at the end nodes and its derivatives
are inconsistent with its value.

Nodes x = 0, 0.2, ..., 1 with values v = 0.1 + 0.5 x (exactly linear data, so the exact spline is the line itself):
  * value(x1) and value(xN) return null_value instead of v1, vN    (test `x <= x_initial .or. x >= x_final`)
  * first_derivative is a cubic spline of nodal derivatives whose first/last entries are forced to 0
    (02_interpolation_linear_1D.f90 InitDerivatives), so d/dx is wrong near both ends
  * second_derivative is ~0 where it should be... 0 (ok), but for curved data it is unrelated to value''
  * outside the node range every derivative returns null_value instead of 0
With --model also runs Pure EFT (2_PEFT_Omega_lin_1, Omega = a) with Omega given as this spline of exact nodes.

usage: /Users/marco/Software/miniforge3/bin/python spline-cubic-endpoints-derivatives.py [--model]
"""

import subprocess
import sys
from pathlib import Path

AUDIT = Path("/private/tmp/claude-501/-Users-marco-IW-eftcamb-update-EFTCAMB/"
             "2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad/audit")
BUILD = AUDIT / "build" / "A0"
WORK = AUDIT / "work" / "verify" / "spline-cubic" / "repro" / Path(__file__).stem
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
real(dl) :: xs(9) = [0._dl, 1._dl, 0.06_dl, 0.5_dl, 0.94_dl, 0.3_dl, 0.9_dl, -0.1_dl, 1.1_dl]
character(len=512) :: path
call get_command_argument(1, path)
err = 0
call Ini%Open(trim(path))
call allocate_parametrized_1D_function(f, 9, 'EFTw', 'w', err, 0)
call f%init_func_from_file(Ini, err)
call f%init_from_file(Ini, err)
do i = 1, size(xs)
    write(*, '(4es24.15)') xs(i), f%value(xs(i)), f%first_derivative(xs(i)), f%second_derivative(xs(i))
end do
end program
"""


def build_run(src: str, ini: str) -> list[list[float]]:
    WORK.mkdir(parents=True, exist_ok=True)
    (WORK / "r.f90").write_text(src)
    (WORK / "r.ini").write_text(ini)
    subprocess.run(["gfortran", "-O0", "-cpp", "-ffree-line-length-none", "-fopenmp",
                    f"-I{BUILD}/camb/Releaselib", f"-I{BUILD}/forutils/Releaselib", "r.f90", "-o", "r",
                    str(BUILD / "camb/Releaselib/camblib.so"), f"-Wl,-rpath,{LIB}"], cwd=WORK, check=True,
                   capture_output=True)
    out = subprocess.run([str(WORK / "r"), str(WORK / "r.ini")], capture_output=True, text=True, check=True).stdout
    return [list(map(float, line.split())) for line in out.splitlines() if line.strip()]


def main():
    xs = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
    for null in (0.0, -1.0):
        ini = f"EFTw_Spline_Pixels = 6\nEFTw_null_value = {null}\n" + "".join(
            f"EFTwx{i + 1} = {x}\nEFTwv{i + 1} = {0.1 + 0.5 * x}\n" for i, x in enumerate(xs))
        rows = build_run(SRC, ini)
        print(f"null_value = {null}:   f(x) = 0.1 + 0.5 x on [0,1] (exact cubic spline of linear data)")
        print(f"{'x':>6} | {'value obs':>12} {'exp':>8} | {'d1 obs':>12} {'exp':>6} | {'d2 obs':>10} {'exp':>4}")
        bad = 0
        for x, v, d1, d2 in rows:
            inside = 0.0 <= x <= 1.0
            ev = 0.1 + 0.5 * x if inside else null
            ed1 = 0.5 if inside else 0.0
            flag = "  <-- WRONG" if abs(v - ev) > 1e-9 or abs(d1 - ed1) > 1e-6 else ""
            bad += bool(flag)
            print(f"{x:6.2f} | {v:12.6f} {ev:8.4f} | {d1:12.6f} {ed1:6.2f} | {d2:10.3e} {0.0:4.1f}{flag}")
        print(f"   -> {bad} of {len(rows)} points wrong\n")
    if "--model" in sys.argv:
        sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "checks"))
        import harness
        out = WORK / "model"
        out.mkdir(parents=True, exist_ok=True)
        ini = harness.PARAMS_DIR / "2_PEFT_Omega_lin_1.ini"
        ov = {"PureEFTmodelOmega": "9", "EFTOmega_Spline_Pixels": "6"}
        ov.update({f"EFTOmegax{i + 1}": str(x) for i, x in enumerate(xs)})
        ov.update({f"EFTOmegav{i + 1}": str(x) for i, x in enumerate(xs)})
        s1 = harness.run_model("A0", ini, out / "spline.npz", overrides=ov)
        s2 = harness.run_model("A0", ini, out / "ref.npz")
        print("model status:", s1["status"], s2["status"])
        if s1["status"] == s2["status"] == "ok":
            c = harness.compare(out / "spline.npz", out / "ref.npz")
            print(f"Pure EFT Omega=a as spline vs Omega=a linear: worst Cl rel diff {c['worst_cl']:.4f} "
                  f"(expected ~0), worst P(k) {c['worst_pk']:.4f}, sigma8 {c['sigma8']:.4f}")


if __name__ == "__main__":
    main()
