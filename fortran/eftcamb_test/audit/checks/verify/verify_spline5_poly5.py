"""Independent verification of finding spline5-poly5-dx6 (02p12_splines5.f90:391).

Different code path from the t4a-param repro (which went through the spline5_parametrization_1D class):
  1. sympy: expand the Hermite basis used by spline5pars and confirm that c(6) is the coefficient of (x-x0)^5.
  2. Fortran driver linked against the frozen A0 camblib.so calling the low-level splines5 API directly
     (spline5pars + poly5/dpoly5..d4poly5 + the high-level spline5 / spline5ders):
       a. data sampled from an exact quintic polynomial: an interpolating quintic spline with "consistent"
          (bctype=2) end conditions must reproduce it to round-off everywhere;
       b. the coefficients returned by spline5pars, evaluated in Python with dx^5, must interpolate and be
          C4 continuous; evaluated with dx^6 they must reproduce the Fortran poly5 output (attribution).

usage: /Users/marco/Software/miniforge3/bin/python verify_spline5_poly5.py [build_dir]
"""

import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import sympy as sp

AUDIT = Path("/private/tmp/claude-501/-Users-marco-IW-eftcamb-update-EFTCAMB/"
             "2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad/audit")
# optional argv[1]: build dir to link against (default frozen A0); outputs go to a per-build subdir
BUILD = Path(sys.argv[1]) if len(sys.argv) > 1 else AUDIT / "build" / "A0"
WORK = AUDIT / "work" / "verify" / "spline5_poly5" / ("driver_" + BUILD.name)
LIB = "/Users/marco/Software/miniforge3/lib"

SRC = r"""
program drv
use precision
use splines5
implicit none
integer :: n, m, i, j
real(dl), allocatable :: xi(:), yi(:), c(:,:), xs(:), y(:), d1(:), d2(:), d3(:), d4(:), yh(:)
read(*,*) n
allocate(xi(n), yi(n), c(0:6, n-1))
do i = 1, n
  read(*,*) xi(i), yi(i)
end do
read(*,*) m
allocate(xs(m), y(m), d1(m), d2(m), d3(m), d4(m), yh(m))
do i = 1, m
  read(*,*) xs(i)
end do
call spline5pars(xi, yi, [2, 2], [0._dl, 0._dl, 0._dl, 0._dl], c)
do j = 1, n-1
  write(*,'(a,7es26.17)') 'C ', c(:, j)
end do
call spline5ders(xi, yi, xs, y, d1, d2, d3, d4)
yh = spline5(xi, yi, xs)
do i = 1, m
  write(*,'(a,7es26.17)') 'E ', xs(i), y(i), d1(i), d2(i), d3(i), d4(i), yh(i)
end do
end program
"""


def build() -> Path:
    WORK.mkdir(parents=True, exist_ok=True)
    (WORK / "drv.f90").write_text(SRC)
    subprocess.run(["gfortran", "-O0", "-cpp", "-ffree-line-length-none", "-fopenmp",
                    f"-I{BUILD}/camb/Releaselib", f"-I{BUILD}/forutils/Releaselib", "drv.f90", "-o", "drv",
                    str(BUILD / "camb/Releaselib/camblib.so"), f"-Wl,-rpath,{LIB}"], cwd=WORK, check=True)
    return WORK / "drv"


def run(exe: Path, xi, yi, xs):
    inp = f"{len(xi)}\n" + "".join(f"{a:.17e} {b:.17e}\n" for a, b in zip(xi, yi))
    inp += f"{len(xs)}\n" + "".join(f"{a:.17e}\n" for a in xs)
    out = subprocess.run([str(exe)], input=inp, capture_output=True, text=True, check=True).stdout
    C = np.array([list(map(float, l.split()[1:])) for l in out.splitlines() if l.startswith("C ")])
    E = np.array([list(map(float, l.split()[1:])) for l in out.splitlines() if l.startswith("E ")])
    return C, E


def sympy_check() -> dict:
    t, h = sp.symbols("t h", positive=True)
    c1, c2, c3, c4, c5, c6 = sp.symbols("c1:7")
    H1 = 1 - 10 * t**3 + 15 * t**4 - 6 * t**5
    H2 = t - 6 * t**3 + 8 * t**4 - 3 * t**5
    H3 = t**2 / 2 - 3 * t**3 / 2 + 3 * t**4 / 2 - t**5 / 2
    H4 = H3.subs(t, 1 - t)
    H5 = -H2.subs(t, 1 - t)
    H6 = H1.subs(t, 1 - t)
    p = sp.expand(c1 * H1 + c2 * H2 + c3 * H3 + c4 * H4 + c5 * H5 + c6 * H6)
    poly = sp.Poly(p, t)
    deg = poly.degree()
    # code's conversion: coefficient of dx^k = code_k ; dx = h t => coefficient of t^k = code_k * h^k
    code = {
        0: c1,
        1: c2 / h,
        2: c3 / h**2 / 2,
        3: (-10 * c1 - 6 * c2 - sp.Rational(3, 2) * c3 + c4 / 2 - 4 * c5 + 10 * c6) / h**3,
        4: (15 * c1 + 8 * c2 + sp.Rational(3, 2) * c3 - c4 + 7 * c5 - 15 * c6) / h**4,
        5: (-6 * c1 - 3 * c2 - c3 / 2 + c4 / 2 - 3 * c5 + 6 * c6) / h**5,
    }
    match = {k: sp.simplify(poly.coeff_monomial(t**k) - code[k] * h**k) == 0 for k in range(6)}
    return {"degree_in_t": deg, "code_coeff_k_matches_t^k": match,
            "coeff_t6": str(poly.coeff_monomial(t**6))}


def pieces(C, x, power6=False):
    """Evaluate spline from Fortran coefficients with the same interval rule (left interval at a node)."""
    xi0 = C[:, 0]
    out = []
    for xx in x:
        ip = max(0, min(len(xi0) - 1, np.searchsorted(xi0, xx, side="right") - 1))
        c = C[ip]
        dx = xx - c[0]
        p6 = dx**6 if power6 else dx**5
        out.append(c[1] + c[2] * dx + c[3] * dx**2 + c[4] * dx**3 + c[5] * dx**4 + c[6] * p6)
    return np.array(out)


def derivs_at(c, dx, power6=False):
    k = 6 if power6 else 5
    coeffs = [c[1], c[2], c[3], c[4], c[5]] + [0.0] * (k - 5) + [c[6]]
    P = np.polynomial.Polynomial(coeffs)
    return [P.deriv(m)(dx) if m else P(dx) for m in range(5)]


def main():
    res = {"sympy": sympy_check()}
    print("sympy:", res["sympy"])
    exe = build()

    # (a) exact quintic polynomial
    q = np.polynomial.Polynomial([0.3, -0.7, 1.1, 2.0, -1.5, 0.9])
    for label, xi in (("uniform", np.linspace(0.0, 1.0, 9)),
                      ("lognodes", np.log(np.array([1e-3, 1e-2, 0.05, 0.1, 0.25, 0.4, 0.6, 0.8, 1.0])))):
        yi = q(xi)
        xs = np.linspace(xi[0], xi[-1], 401)
        C, E = run(exe, xi, yi, xs)
        scale = np.max(np.abs(q(xs)))
        err_f = np.max(np.abs(E[:, 1] - q(xs))) / scale
        err_hl = np.max(np.abs(E[:, 6] - q(xs))) / scale
        err_d1 = np.max(np.abs(E[:, 2] - q.deriv(1)(xs))) / np.max(np.abs(q.deriv(1)(xs)))
        err_d4 = np.max(np.abs(E[:, 5] - q.deriv(4)(xs))) / np.max(np.abs(q.deriv(4)(xs)))
        err_py5 = np.max(np.abs(pieces(C, xs) - q(xs))) / scale
        err_py6_vs_f = np.max(np.abs(pieces(C, xs, True) - E[:, 1])) / scale
        r = {"value_poly5_relerr": err_f, "value_spline5_relerr": err_hl, "d1_relerr": err_d1,
             "d4_relerr": err_d4, "python_dx5_relerr": err_py5, "python_dx6_vs_fortran_poly5": err_py6_vs_f}
        res[f"quintic_{label}"] = r
        print(f"quintic data, {label} nodes:", {k: f"{v:.2e}" for k, v in r.items()})

    # (b) non-polynomial data: continuity of value and derivatives across interior nodes
    nodes = np.array([1e-3, 0.01, 0.05, 0.1, 0.25, 0.4, 0.6, 0.8, 1.0])
    for logx in (False, True):
        xi = np.log(nodes) if logx else nodes.copy()
        yi = 0.2 * np.sin(4 * nodes) + 0.1
        C, E = run(exe, xi, yi, xi)
        jumps = {5: np.zeros((len(xi) - 2, 5)), 6: np.zeros((len(xi) - 2, 5))}
        interp = {5: 0.0, 6: 0.0}
        for j in range(1, len(xi) - 1):
            for pw in (5, 6):
                L = derivs_at(C[j - 1], xi[j] - C[j - 1, 0], pw == 6)
                R = derivs_at(C[j], 0.0, pw == 6)
                jumps[pw][j - 1] = [abs(L[m] - R[m]) / max(abs(R[m]), 1e-12) for m in range(5)]
                interp[pw] = max(interp[pw], abs(L[0] - yi[j]) / abs(yi[j]))
        # Fortran poly5 at nodes (left interval for nodes 2,3 by iix) vs python dx^6 / dx^5
        f_node = E[:, 1]
        py6 = pieces(C, xi, True)
        py5 = pieces(C, xi, False)
        r = {
            "dx5_max_rel_jump_derivs0to4": jumps[5].max(axis=0).tolist(),
            "dx6_max_rel_jump_derivs0to4": jumps[6].max(axis=0).tolist(),
            "dx5_left_limit_interp_relerr": interp[5],
            "dx6_left_limit_interp_relerr": interp[6],
            "fortran_poly5_at_nodes_vs_node_values_max_rel": float(np.max(np.abs(f_node - yi) / np.abs(yi))),
            "python_dx6_vs_fortran_poly5_at_nodes": float(np.max(np.abs(py6 - f_node))),
            "python_dx5_at_nodes_vs_node_values_max_rel": float(np.max(np.abs(py5 - yi) / np.abs(yi))),
        }
        res[f"sin_logx_{logx}"] = r
        print(f"sin data, Logx={logx}:")
        for k, v in r.items():
            print(f"   {k}: {np.array2string(np.asarray(v), precision=2)}")
    (WORK / "result.json").write_text(json.dumps(res, indent=1, default=str))


if __name__ == "__main__":
    main()
