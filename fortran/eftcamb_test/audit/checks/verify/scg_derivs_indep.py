"""Independent verifier check for finding scg-higher-derivatives-wrong.

Builds the exact Scaling Cubic background flow in N = ln a with sympy (state x1, y1, y2, OmegaRad, ln adotoa, ln a),
takes exact total derivatives of h, f, EFTLambda (as adotoa^2 L) and EFTGamma2V, and compares them, at random state
points, with the transcribed Fortran formulas of 008p6_Scaling_Cubic.f90 (code's own lower derivatives fed in).
"""

import numpy as np
import sympy as sp

A, lam, b1, b2 = sp.symbols("A lam b1 b2")
x1, y1, y2, Om, lH, lna = sp.symbols("x1 y1 y2 Om lH lna")
s6 = sp.sqrt(6)
D = 1 + 6 * A**2 - 2 * A * lam
f = 3 * (-3 * A + s6 * (-1 + 2 * A * lam) * x1 + 3 * A * (1 - 2 * A * lam) * x1**2 + A * Om - 3 * A * y2**2
         + b2 * y2**2 - 3 * A * y1**2 + b1 * y1**2) / D
h = ((-1 + 2 * A * lam) * Om - 3 * (1 + 12 * A**2 - 2 * A * lam + 2 * s6 * A * (1 - 2 * A * lam) * x1
     + (1 - 2 * A * lam)**2 * x1**2 + (-1 + 2 * A * lam - 2 * A * b2) * y2**2 - y1**2 + 2 * A * lam * y1**2
     - 2 * A * b1 * y1**2)) / (2 * D)
flow = {x1: f / s6 - x1 * h, y1: -sp.sqrt(sp.Rational(3, 2)) * b1 * x1 * y1 - y1 * h,
        y2: -sp.sqrt(sp.Rational(3, 2)) * b2 * x1 * y2 - y2 * h, Om: -2 * Om * (2 + h), lH: 1 + h, lna: 1}


def dN(e):
    return sum(sp.diff(e, v) * r for v, r in flow.items())


H = sp.exp(lH)  # adotoa (conformal), d ln adotoa/dN = 1 + h
a = sp.exp(lna)
hp, hpp = dN(h), None
hpp = dN(hp)
fp = dN(f)
fpp = dN(fp)
X1p = flow[x1]
X1pp = dN(X1p)
X1ppp = dN(X1pp)
Lam = 3 * H**2 * (x1**2 - y1**2 - y2**2) - A * H**2 * (6 * lam * x1**2 + 2 * s6 * (x1 * h + X1p))
Lamdot = H * (dN(Lam) - 2 * Lam)  # a^2 d(Lambda)/dtau convention
Lamdotdot = H * (dN(Lamdot) - 2 * Lamdot)
G2 = -2 * s6 * A * H * x1 / a  # / h0 dropped (constant)
G2P = dN(G2) / a
G2PP = dN(G2P) / a
G2PPP = dN(G2PP) / a
Hdd = H * dN(H * dN(H))
Hddd = H * dN(Hdd)

rng = np.random.default_rng(1)
syms = [A, lam, b1, b2, x1, y1, y2, Om, lH, lna]
exact = [h, hp, hpp, f, fp, fpp, X1p, X1pp, X1ppp, Lamdot, Lamdotdot, G2PP, G2PPP, Hdd, Hddd]
fn = sp.lambdify(syms, exact, "numpy")
names = ["h", "hp", "hpp", "f", "fp", "fpp", "x1p", "x1pp", "x1ppp", "Ld", "Ldd", "G2PP", "G2PPP", "Hdd", "Hddd"]


def code(A, lam, b1, b2, x1, y1, y2, Om, lH, lna):
    s6 = np.sqrt(6.0)
    adotoa = np.exp(lH)
    a = np.exp(lna)
    Dd = 1 + 6 * A**2 - 2 * A * lam
    f = 3 * (-3 * A + s6 * (-1 + 2 * A * lam) * x1 + 3 * A * (1 - 2 * A * lam) * x1**2 + A * Om - 3 * A * y2**2
             + b2 * y2**2 - 3 * A * y1**2 + b1 * y1**2) / Dd
    h = ((-1 + 2 * A * lam) * Om - 3 * (1 + 12 * A**2 - 2 * A * lam + 2 * s6 * A * (1 - 2 * A * lam) * x1
         + (1 - 2 * A * lam)**2 * x1**2 + (-1 + 2 * A * lam - 2 * A * b2) * y2**2 - y1**2 + 2 * A * lam * y1**2
         - 2 * A * b1 * y1**2)) / (2 + 12 * A**2 - 4 * A * lam)
    x1p = f / s6 - x1 * h
    y1p = -np.sqrt(1.5) * b1 * x1 * y1 - y1 * h
    y2p = -np.sqrt(1.5) * b2 * x1 * y2 - y2 * h
    Omp = -2 * Om * (2 + h)
    fprime = (-3 * (-A * Omp + (-1 + 2 * A * lam) * (-s6 + 6 * A * x1) * x1p + 6 * A * y2 * y2p - 2 * b2 * y2 * y2p
                    + 6 * A * y1 * y1p - 2 * b1 * y1 * y1p)) / Dd
    hprime = (0.5 * Omp * (-1 + 2 * A * lam) + 3 * ((-1 + 2 * A * lam) * (s6 * A + x1 - 2 * A * lam * x1) * x1p
              + (1 - 2 * A * lam + 2 * A * b2) * y2 * y2p + (1 - 2 * A * lam + 2 * A * b1) * y1 * y1p)) / Dd
    x1pp = fprime / s6 - hprime * x1 - h * x1p
    y1pp = -np.sqrt(1.5) * b1 * (x1p * y1 + x1 * y1p) - hprime * y1 - h * y1p
    y2pp = -np.sqrt(1.5) * b2 * (x1p * y2 + x1 * y2p) - hprime * y2 - h * y2p
    Ompp = -2 * Omp * (2 + h) - 2 * Om * hprime
    D2 = 1 + 2 * A * (3 * A - lam)
    fpprime = -3 / D2 * (-A * Ompp + 6 * A * (2 * A * lam - 1) * x1p**2 + 2 * (3 * A - b1) * (y1p**2 + y1 * y1pp)
                         + 2 * (3 * A - b2) * (y2p**2 + y2 * y2pp))
    hpprime = -1 / D2 * (0.5 * Ompp * (2 * A * lam - 1) + 3 * (-(1 - 2 * A * lam)**2 * x1p**2 + (2 * A * lam - 1)
                         * (s6 * A + x1 - 2 * A * lam * x1) * x1pp + (1 - 2 * A * lam + 2 * A * b1) * (y1p**2 + y1 * y1pp)
                         + (1 + 2 * A * lam + 2 * A * b2) * (y2p**2 + y2 * y2pp)))
    x1ppp = np.sqrt(1 / 6) * fpprime - hpprime * x1 - h * x1pp - 2 * hprime * x1p
    Ld = -2 * adotoa**3 * (2 * s6 * A * x1 * h**2 + 3 * h * ((-1 + 2 * A * lam) * x1**2 + y1**2 + y2**2 + s6 * A * x1p)
                           + x1 * (s6 * A * hprime + 3 * (-1 + 2 * A * lam) * x1p) + 3 * y1 * y1p + 3 * y2 * y2p
                           + s6 * A * x1pp)
    q = 2 * A * lam - 1
    Ldd = -2 * adotoa**4 * (6 * s6 * A * x1 * h**3 + 3 * hprime * (y1**2 + y2**2 + x1**2 * q) + 4 * s6 * A * x1p * hprime
                            + 3 * q + A * s6 * (x1pp + x1ppp)
                            + h**2 * (2 * s6 * A * x1 + 9 * (y1**2 + y2**2) + 9 * q * x1**2 + 11 * s6 * x1p * A)
                            + x1p * (s6 * A * hprime + 3 * (x1p + x1pp) * q + s6 * A * hpprime)
                            + 3 * (y1 * (y1p + y1pp) + y1p**2 + y2 * (y2p + y2pp) + y2p**2)
                            + h * (3 * x1**2 * q + x1 * (15 * q + 7 * s6 * A * hprime) + 3 * s6 * a * (x1p + 2 * x1pp)
                                   + 3 * (y1**2 + y2**2 + 5 * y1 * y1p + 5 * y2 * y2p)))
    G2PP = -2 * s6 * A * adotoa / a**3 * (-x1p - h * x1 + 2 * x1p * h + hprime * x1 + h**2 * x1 + x1pp)
    G2PPP = -2 * s6 * A * adotoa / a**4 * (x1 * (2 * h - 3 * hprime - 3 * h**2 + 3 * h * hprime + h**3 + hpprime)
                                           + x1p * (2 - 6 * h + 3 * h * 2 + 3 * hprime) + x1pp * (3 * h - 3) + x1ppp)
    Hdot = adotoa**2 * (1 + h)
    # Hdotdot/Hdotdotdot from (h, hprime, hpprime): Hdot = H^2(1+h); d/dtau = H d/dN
    Hdd = adotoa**3 * (2 * (1 + h)**2 + hprime)
    Hddd = adotoa**4 * (3 * (1 + h) * (2 * (1 + h)**2 + hprime) + 4 * (1 + h) * hprime + hpprime)
    return [h, hprime, hpprime, f, fprime, fpprime, x1p, x1pp, x1ppp, Ld, Ldd, G2PP, G2PPP, Hdd, Hddd]


worst = {n: 0.0 for n in names}
for trial in range(200):
    pars = [rng.uniform(-0.5, 0.5), rng.uniform(-2, 2), rng.uniform(-2, 2), rng.uniform(-2, 2),
            rng.uniform(-0.6, 0.6), rng.uniform(-0.6, 0.6), rng.uniform(-0.6, 0.6), rng.uniform(0, 0.5),
            rng.uniform(-3, 1), rng.uniform(-4, 0)]
    ex = fn(*pars)
    co = code(*pars)
    for n, e, c in zip(names, ex, co):
        worst[n] = max(worst[n], abs(c - e) / max(abs(e), 1e-12))
for n in names:
    print(f"{n:7s} max rel diff code vs exact = {worst[n]:.3e}")

# y2-only check of hpprime sign: x1=y1=0, Om=0
pars = [0.1, 0.5, 0.3, 0.7, 0.0, 0.0, 0.4, 0.0, 0.0, -1.0]
ex, co = fn(*pars), code(*pars)
print("x1=y1=Om=0: hpprime code/exact =", co[2] / ex[2])
pars = [0.1, 0.5, 0.3, 0.7, 0.0, 0.4, 0.0, 0.0, 0.0, -1.0]
ex, co = fn(*pars), code(*pars)
print("x1=y2=Om=0: hpprime code/exact =", co[2] / ex[2], "(y1 coefficient right -> ratio -1 means pure sign error)")
