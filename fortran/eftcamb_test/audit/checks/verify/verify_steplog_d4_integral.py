"""Adversarial verification of finding steplog-d4-and-integral (t4a-param).

Independent of the t4a-param repro:
  * own Fortran driver: instantiates steplog_parametrization_1D directly (no allocator, no ini) and calls
    init_parameters(array); also the constant and CPL parametrizations to establish the integral convention.
  * reference derivatives with mpmath.diff (50 digits, numerical differentiation of the closed form w(a)),
    not sympy; and a finite-difference derivative of the code's own third_derivative (no symbolic input).
  * reference integral by mpmath.quad of I(a) = a^2 exp(3 int_{ln a}^0 (1 + w(e^N)) dN)
    (convention cross-checked against the code's constant-w integral = a^(-1-3w)).

usage: python verify_steplog_d4_integral.py [--tree A0|<path to build dir>]
"""

import argparse
import json
import subprocess
from pathlib import Path

import mpmath as mp

AUDIT = Path("/private/tmp/claude-501/-Users-marco-IW-eftcamb-update-EFTCAMB/"
             "2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad/audit")
OUT = AUDIT / "work" / "verify" / "steplog"
LIB = "/Users/marco/Software/miniforge3/lib"

SRC = r"""
program vsteplog
use precision
use EFTCAMB_steplog_parametrizations_1D
use EFTCAMB_constant_parametrization_1D
use EFTCAMB_CPL_parametrizations_1D
implicit none
type(steplog_parametrization_1D) :: f
type(constant_parametrization_1D) :: c
type(CPL_parametrization_1D) :: cpl
real(dl) :: p(4, 5), xs(13), h
integer :: i, j
p(:, 1) = [0._dl, 1.5_dl, 0.5_dl, 0.5_dl]
p(:, 2) = [-1._dl, -0.8_dl, 0.3_dl, 0.5_dl]
p(:, 3) = [-0.9_dl, -1.1_dl, 0.7_dl, 1.0_dl]
p(:, 4) = [0.1_dl, 0.0_dl, 0.05_dl, 2.0_dl]
p(:, 5) = [-1.0_dl, -0.9_dl, 0.5_dl, 0.3_dl]
xs = [0.01_dl, 0.05_dl, 0.1_dl, 0.2_dl, 0.3_dl, 0.5_dl, 0.7_dl, 0.9_dl, 0.99_dl, 0.999999_dl, 0.999999999999_dl, &
      1.0_dl, 1.0_dl - epsilon(1._dl)]
call f%set_param_number()
do j = 1, 5
    call f%init_parameters(p(:, j))
    do i = 1, size(xs)
        h = 1.e-4_dl*xs(i)
        write(*, '(a,i2,8es26.17)') 'S', j, xs(i), f%value(xs(i)), f%first_derivative(xs(i)), &
            f%second_derivative(xs(i)), f%third_derivative(xs(i)), f%fourth_derivative(xs(i)), f%integral(xs(i)), &
            (-f%third_derivative(xs(i) + 2*h) + 8*f%third_derivative(xs(i) + h) - 8*f%third_derivative(xs(i) - h) &
             + f%third_derivative(xs(i) - 2*h))/(12*h)
    end do
end do
call c%set_param_number()
call c%init_parameters([-0.8_dl])
call cpl%set_param_number()
call cpl%init_parameters([-0.9_dl, 0.2_dl])
do i = 1, size(xs)
    write(*, '(a,3es26.17)') 'C', xs(i), c%integral(xs(i)), cpl%integral(xs(i))
end do
end program
"""

PARAMS = [(0.0, 1.5, 0.5, 0.5), (-1.0, -0.8, 0.3, 0.5), (-0.9, -1.1, 0.7, 1.0), (0.1, 0.0, 0.05, 2.0),
          (-1.0, -0.9, 0.5, 0.3)]


def build(builddir: Path) -> Path:
    OUT.mkdir(parents=True, exist_ok=True)
    tag = builddir.parent.name if builddir.name == "camb" else builddir.name
    (OUT / "vsteplog.f90").write_text(SRC)
    exe = OUT / f"vsteplog_{tag}"
    mods = builddir / "Releaselib"
    fmods = builddir.parent / "forutils" / "Releaselib"
    subprocess.run(["gfortran", "-O0", "-cpp", "-ffree-line-length-none", "-fopenmp", f"-I{mods}", f"-I{fmods}",
                    "-J", str(OUT), "vsteplog.f90", "-o", str(exe), str(mods / "camblib.so"), f"-Wl,-rpath,{LIB}"],
                   cwd=OUT, check=True)
    return exe


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--build", default=str(AUDIT / "build" / "A0" / "camb"))
    args = ap.parse_args()
    exe = build(Path(args.build))
    out = subprocess.run([str(exe)], capture_output=True, text=True, check=True).stdout.splitlines()
    mp.mp.dps = 50
    res = {"sets": [], "conv": []}
    for line in out:
        tok = line.split()
        if tok[0] == "S":
            j = int(tok[1]) - 1
            x, v, d1, d2, d3, d4, integ, d4fd = map(float, tok[2:])
            v1, v2, at, de = (mp.mpf(str(t)) for t in PARAMS[j])

            def w(a):
                u = (mp.log(a) - mp.log(at)) / de
                return (v2 - v1) / 2 * u / mp.sqrt(1 + u * u) + (v1 + v2) / 2

            xm = mp.mpf(x)
            ex = [w(xm)] + [mp.diff(w, xm, n) for n in range(1, 5)]
            ei = xm**2 * mp.exp(3 * mp.quad(lambda n: 1 + w(mp.exp(n)), [mp.log(xm), 0]))
            obs = [v, d1, d2, d3, d4]

            def rel(o, e):
                return float(abs(mp.mpf(o) - e) / max(abs(e), mp.mpf("1e-300"))) if o == o else float("nan")

            res["sets"].append({"set": j + 1, "params": PARAMS[j], "a": x, "obs": obs + [integ], "d4_fd_of_d3": d4fd,
                                "exact": [float(e) for e in ex] + [float(ei)],
                                "rel": [rel(o, e) for o, e in zip(obs, ex)] + [rel(integ, ei)],
                                "rel_fd_vs_exact": rel(d4fd, ex[4])})
        else:
            x, ci, cpli = map(float, tok[1:])
            res["conv"].append({"a": x, "const_obs": ci, "const_exact": x ** (-1 - 3 * (-0.8)), "cpl_obs": cpli})
    tag = Path(args.build).parent.name
    (OUT / f"verify_steplog_{tag}.json").write_text(json.dumps(res, indent=1))
    print(f"build: {args.build}")
    print("set  params                 a            | rel err d1      d2       d3       d4       I   "
          "| d4 obs          d4 exact        FD(d3)          | I obs        I exact")
    for r in res["sets"]:
        e = r["rel"]
        print(f"{r['set']:3d} {str(r['params']):24s} {r['a']:.12f} | " + " ".join(f"{t:8.1e}" for t in e[1:]) +
              f" | {r['obs'][4]:14.6g} {r['exact'][4]:14.6g} {r['d4_fd_of_d3']:14.6g} | {r['obs'][5]:.10f} "
              f"{r['exact'][5]:.10f}")
    print("convention: constant w=-0.8 integral(a) vs a^(-1-3w); CPL(-0.9,0.2) integral")
    for r in res["conv"]:
        print(f"  a={r['a']:.16f}  const {r['const_obs']:.12f} vs {r['const_exact']:.12f}   cpl {r['cpl_obs']:.12f}")


if __name__ == "__main__":
    main()
