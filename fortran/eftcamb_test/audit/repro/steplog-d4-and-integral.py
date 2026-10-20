"""Repro (t4a-param): steplog parametrization (allocator flag 8)
  1) fourth_derivative is a wrong expression (uses (1 - log(x/aT)/Delta^2)^(k/2) instead of (1 + (log(x/aT)/Delta)^2)):
     wrong everywhere and NaN for a > aT exp(Delta^2)  (04p018_steplog_parametrizations_1D.f90:340-352)
  2) integral(a) divides by log(x): at a = 1 exactly one factor becomes 1**Inf = 1 and the result is not 1
     (04p018_steplog_parametrizations_1D.f90:372-376).  integral is grhov_t/grhov in designer f(R)/ShiftSym/3+1.
Expected values: sympy derivative of the documented function; the EFTCAMB integral convention
I(a) = a^2 exp(-3 int_1^a (1+w)/a' da') (designer f(R): grhov_t = grhov*integral(a)), which is 1 at a = 1.

With --model also runs Pure EFT 2_PEFT_Omega_lin_1.ini with Omega given as the steplog function
(PureEFTmodelOmega = 8, v1 = 0, v2 = 1.5, aT = 0.5, delta = 0.5, EFT_mass_stability = F): EFTOmegaPPPP is NaN for
a > aT exp(delta^2) = 0.642, the NaN check of the timestep cache flags it and the model is rejected as unstable.
Expected: ok (a tree with the corrected fourth derivative runs it without error).

usage: /Users/marco/Software/miniforge3/bin/python steplog-d4-and-integral.py [--model]
"""

import math
import subprocess
from pathlib import Path

import sympy as sp
from scipy.integrate import quad

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
real(dl) :: xs(7) = [0.05_dl, 0.2_dl, 0.3_dl, 0.5_dl, 0.9_dl, 0.999999_dl, 1.0_dl]
character(len=512) :: path
call get_command_argument(1, path)
err = 0
call Ini%Open(trim(path))
call allocate_parametrized_1D_function(f, 8, 'EFTw', 'w', err, 0)
call f%init_func_from_file(Ini, err)
call f%init_from_file(Ini, err)
do i = 1, size(xs)
    write(*, '(4es26.17)') xs(i), f%value(xs(i)), f%fourth_derivative(xs(i)), f%integral(xs(i))
end do
end program
"""


def main():
    v1, v2, at, dl = -1.0, -0.8, 0.3, 0.5
    WORK.mkdir(parents=True, exist_ok=True)
    (WORK / "r.f90").write_text(SRC)
    (WORK / "r.ini").write_text(f"EFTw_v1 = {v1}\nEFTw_v2 = {v2}\nEFTw_at = {at}\nEFTw_delta = {dl}\n")
    subprocess.run(["gfortran", "-O0", "-cpp", "-ffree-line-length-none", "-fopenmp",
                    f"-I{BUILD}/camb/Releaselib", f"-I{BUILD}/forutils/Releaselib", "r.f90", "-o", "r",
                    str(BUILD / "camb/Releaselib/camblib.so"), f"-Wl,-rpath,{LIB}"], cwd=WORK, check=True,
                   capture_output=True)
    out = subprocess.run([str(WORK / "r"), str(WORK / "r.ini")], capture_output=True, text=True, check=True).stdout
    a = sp.symbols("a", positive=True)
    u = (sp.log(a) - sp.log(sp.Rational(str(at)))) / sp.Rational(str(dl))
    w = sp.Rational(str(v2 - v1)) / 2 * u / sp.sqrt(1 + u**2) + sp.Rational(str(v1 + v2)) / 2
    d4 = sp.lambdify(a, sp.diff(w, a, 4))
    wf = sp.lambdify(a, w)
    print(f"steplog w(a): v1={v1}, v2={v2}, aT={at}, delta={dl}  (NaN expected in code for a > {at * math.exp(dl**2):.3f})")
    print(f"{'a':>9} | {'d4 obs':>14} {'d4 exact':>14} | {'integral obs':>14} {'integral exact':>14}")
    for line in out.splitlines():
        x, v, c4, ci = map(float, line.split())
        e4 = d4(x)
        ei = x * x * math.exp(3 * quad(lambda n: 1 + wf(math.exp(n)), math.log(x), 0, epsabs=1e-13, epsrel=1e-13)[0])
        print(f"{x:9.6f} | {c4:14.6g} {e4:14.6g} | {ci:14.10f} {ei:14.10f}")


def model():
    import sys
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "checks"))
    import harness
    ini = harness.PARAMS_DIR / "2_PEFT_Omega_lin_1.ini"
    ov = {"EFT_mass_stability": "F", "PureEFTmodelOmega": "8", "EFTOmega_v1": "0.0", "EFTOmega_v2": "1.5",
          "EFTOmega_at": "0.5", "EFTOmega_delta": "0.5"}
    st = harness.run_model("A0", ini, WORK / "peft_steplog.npz", overrides=ov, threads=1)
    log = Path(st["log"]).read_text(errors="replace").splitlines() if st.get("log") else []
    nan = [t.strip() for t in log if "is NaN" in t][:3]
    print(f"Pure EFT, Omega = steplog(v1=0, v2=1.5, aT=0.5, delta=0.5): expected ok | observed {st['status']} "
          f"{' / '.join(nan)}")


if __name__ == "__main__":
    import sys
    main()
    if "--model" in sys.argv:
        model()
