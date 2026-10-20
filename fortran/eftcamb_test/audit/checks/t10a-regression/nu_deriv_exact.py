"""Accuracy of the massive-neutrino pressure derivatives used by EFTCAMB (A0 and B1) against exact quadrature.

Reads the outputs of fortran/nu_deriv_test.f90 (columns: am P dP1_ana dP1_num dP2_ana dP2_num dP3_ana dP3_num
drho_ana drho_num rho), i.e. P(am) and its first three derivatives with respect to log(am) as computed by the tree's
ThermalNuBackground (pidot, pidotdot, pidotdotdot with adotoa=1, Hdot=Hdotdot=0), and compares them with
  P(am) = int q^4/(3 eps) f(q) dq / F,  eps = sqrt(q^2 + am^2),  f = 1/(e^q+1),  F = 7 pi^4/120
differentiated exactly (mpmath, 40 digits).  Errors are quoted relative to |P| (the natural scale of
d^nP/dlog(am)^n) and, for the derivatives, also relative to max(|P|, |d^nP|).

Usage: python nu_deriv_exact.py NU_A0.txt NU_B1.txt [--every 10]
"""

from __future__ import annotations

import argparse

import mpmath as mp
import numpy as np

mp.mp.dps = 40
F = 7 * mp.pi**4 / 120


def p_exact(am):
    return mp.quad(lambda q: q**4 / (3 * mp.sqrt(q * q + am * am)) / (mp.exp(q) + 1), [0, 1, 5, 20, 80, mp.inf]) / F


def derivs(am: float):
    x0 = mp.log(am)
    f = lambda x: p_exact(mp.exp(x))  # noqa: E731
    return [float(mp.diff(f, x0, n)) for n in range(4)]


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("nu_a0")
    p.add_argument("nu_b1")
    p.add_argument("--every", type=int, default=10)
    a = p.parse_args()
    A = np.loadtxt(a.nu_a0)
    B = np.loadtxt(a.nu_b1)
    idx = list(range(0, len(A), a.every))
    # always include the points nearest the A0 fit boundaries (0.3, 0.42, 70)
    for b in (0.3, 0.42, 70.0):
        idx += [int(np.argmin(np.abs(A[:, 0] - b)))]
    idx = sorted(set(idx))
    print(f"{'am':>10s} | {'A0 dP1':>9s} {'dP2':>9s} {'dP3':>9s} {'P':>9s} | {'B1 dP1':>9s} {'dP2':>9s} {'dP3':>9s} {'P':>9s}"
          "   (|tree-exact|/|P|)")
    worst = {"A0": np.zeros(4), "B1": np.zeros(4)}
    for i in idx:
        am = A[i, 0]
        ex = derivs(am)
        row = []
        for lab, T in (("A0", A), ("B1", B)):
            vals = [T[i, 2], T[i, 4], T[i, 6], T[i, 1]]
            exv = [ex[1], ex[2], ex[3], ex[0]]
            errs = np.array([abs(v - e) / abs(ex[0]) for v, e in zip(vals, exv)])
            worst[lab] = np.maximum(worst[lab], errs)
            row.append(errs)
        print(f"{am:10.4g} | " + " ".join(f"{e:9.2e}" for e in row[0]) + " | " + " ".join(f"{e:9.2e}" for e in row[1]))
    for lab in ("A0", "B1"):
        print(f"worst {lab}: dP1 {worst[lab][0]:.2e}  dP2 {worst[lab][1]:.2e}  dP3 {worst[lab][2]:.2e}  P {worst[lab][3]:.2e}")


if __name__ == "__main__":
    main()
