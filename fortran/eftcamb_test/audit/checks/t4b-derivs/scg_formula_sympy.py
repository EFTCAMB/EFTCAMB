"""Track t4b-derivs: offline symbolic check of the Scaling Cubic Galileon (008p6_Scaling_Cubic.f90) derivative formulas.

The model evolves x1, y1, y2 in N = ln a with
    x1' = f/sqrt6 - x1 h,  y_i' = -sqrt(3/2) beta_i x1 y_i - y_i h,  OmegaR' = -2 OmegaR (2 + h)
and f(x1, y1, y2, OmegaR), h(x1, y1, y2, OmegaR) as coded in EFTCAMBScalingCubicBackgroundEFTFunctions.  The code
hand-writes fprime, hprime, fpprime, hpprime, x1ppp and from them EFTLambdadot(dot), EFTGamma2P(PP/PPP) and Hdotdotdot.
Here every hand-written expression is transcribed verbatim from the Fortran and compared, at random points of
(x1, y1, y2, OmegaR) and for the default model parameters, with the exact derivatives obtained by sympy along the
ODE flow.  Conventions (03_abstract_EFTCAMB_cache.f90): adotoa = conformal H, Hdot = H^2 (1 + h);
EFTLambdadot = a^2 d(Lambda/a^2)/dtau = H (dLambda/dN - 2 Lambda), EFTLambdadotdot = H (d EFTLambdadot/dN - 2 EFTLambdadot);
EFTGamma2P(PP, PPP) = d^n Gamma2 / da^n.

Usage: python scg_formula_sympy.py
"""

from __future__ import annotations

import numpy as np
import sympy as sp

A, lam, b1, b2 = -0.3, 154.0, 100.0, 0.7  # defaults of EFTCAMBScalingCubicInitModelParametersFromFile
s6 = np.sqrt(6.0)

N = sp.symbols("N")
X1, Y1, Y2, OR, HH = (sp.Function(n)(N) for n in ("x1", "y1", "y2", "OR", "H"))
sA, sl, sb1, sb2 = sp.nsimplify(A), sp.nsimplify(lam), sp.nsimplify(b1), sp.nsimplify(b2)
S6 = sp.sqrt(6)


def f_expr(x1, y1, y2, Or):
    return 3 * (-3 * sA + S6 * (-1 + 2 * sA * sl) * x1 + 3 * sA * (1 - 2 * sA * sl) * x1**2 + sA * Or
                - 3 * sA * y2**2 + sb2 * y2**2 - 3 * sA * y1**2 + sb1 * y1**2) / (1 + 6 * sA**2 - 2 * sA * sl)


def h_expr(x1, y1, y2, Or):
    return ((-1 + 2 * sA * sl) * Or - 3 * (1 + 12 * sA**2 - 2 * sA * sl + 2 * S6 * sA * (1 - 2 * sA * sl) * x1
            + (1 - 2 * sA * sl) ** 2 * x1**2 + (-1 + 2 * sA * sl - 2 * sA * sb2) * y2**2 - y1**2
            + 2 * sA * sl * y1**2 - 2 * sA * sb1 * y1**2)) / (2 + 12 * sA**2 - 4 * sA * sl)


fN = f_expr(X1, Y1, Y2, OR)
hN = h_expr(X1, Y1, Y2, OR)
flow = {
    sp.Derivative(X1, N): fN / S6 - X1 * hN,
    sp.Derivative(Y1, N): -sp.sqrt(sp.Rational(3, 2)) * sb1 * X1 * Y1 - Y1 * hN,
    sp.Derivative(Y2, N): -sp.sqrt(sp.Rational(3, 2)) * sb2 * X1 * Y2 - Y2 * hN,
    sp.Derivative(OR, N): -2 * OR * (2 + hN),
    sp.Derivative(HH, N): HH * (1 + hN),
}


def dN(e):
    return sp.diff(e, N).subs(flow)


x1s, y1s, y2s, ors, Hs = sp.symbols("x1 y1 y2 Or H")
back = {X1: x1s, Y1: y1s, Y2: y2s, OR: ors, HH: Hs}


def lam_fn(e):
    return sp.lambdify((x1s, y1s, y2s, ors, Hs), e.subs(back), "numpy")


# exact derivatives along the flow
x1p_e = flow[sp.Derivative(X1, N)]
exact = {"f": fN, "h": hN, "x1p": x1p_e}
exact["hprime"] = dN(hN)
exact["hpprime"] = dN(exact["hprime"])
exact["fprime"] = dN(fN)
exact["fpprime"] = dN(exact["fprime"])
exact["x1pp"] = dN(x1p_e)
exact["x1ppp"] = dN(exact["x1pp"])
Lam = 3 * HH**2 * (X1**2 - Y1**2 - Y2**2) - sA * HH**2 * (6 * sl * X1**2 + 2 * S6 * (X1 * hN + x1p_e))
Ldot = HH * (dN(Lam) - 2 * Lam)
exact["EFTLambdadot"] = Ldot
exact["EFTLambdadotdot"] = HH * (dN(Ldot) - 2 * Ldot)
Hdot = HH**2 * (1 + hN)
Hdd = HH * dN(Hdot)
exact["Hdotdot"] = Hdd
exact["Hdotdotdot"] = HH * dN(Hdd)
# Gamma2 (a and h0 factored out: Gamma2 = -2 sqrt6 A H x1 / (a h0)); d/da = (1/a) d/dN
aS = sp.exp(N)
G2 = -2 * S6 * sA * HH * X1 / aS
exact["G2P"] = dN(G2) / aS
exact["G2PP"] = dN(dN(G2) / aS) / aS
exact["G2PPP"] = dN(dN(dN(G2) / aS) / aS) / aS
ex_fn = {k: lam_fn(v) for k, v in exact.items() if not k.startswith("G2")}
g2_fn = {k: sp.lambdify((x1s, y1s, y2s, ors, Hs, N), v.subs(back), "numpy") for k, v in exact.items()
         if k.startswith("G2")}


def code(x1, y1, y2, Or, H, a):
    """Verbatim transcription of the Fortran expressions (BackgroundEFTFunctions, SecondOrder, ComputeHubbleDer)."""
    D = 1 + 6 * A**2 - 2 * A * lam
    f = 3 * (-3 * A + s6 * (-1 + 2 * A * lam) * x1 + 3 * A * (1 - 2 * A * lam) * x1**2 + A * Or - 3 * A * y2**2
             + b2 * y2**2 - 3 * A * y1**2 + b1 * y1**2) / D
    h = ((-1 + 2 * A * lam) * Or - 3 * (1 + 12 * A**2 - 2 * A * lam + 2 * s6 * A * (1 - 2 * A * lam) * x1
         + (1 - 2 * A * lam) ** 2 * x1**2 + (-1 + 2 * A * lam - 2 * A * b2) * y2**2 - y1**2 + 2 * A * lam * y1**2
         - 2 * A * b1 * y1**2)) / (2 + 12 * A**2 - 4 * A * lam)
    x1p = f / s6 - x1 * h
    y1p = -np.sqrt(1.5) * b1 * x1 * y1 - y1 * h
    y2p = -np.sqrt(1.5) * b2 * x1 * y2 - y2 * h
    Orp = -2 * Or * (2 + h)
    fprime = (-3 * (-A * Orp + (-1 + 2 * A * lam) * (-s6 + 6 * A * x1) * x1p + 6 * A * y2 * y2p - 2 * b2 * y2 * y2p
                    + 6 * A * y1 * y1p - 2 * b1 * y1 * y1p)) / D
    hprime = (0.5 * Orp * (-1 + 2 * A * lam) + 3 * ((-1 + 2 * A * lam) * (s6 * A + x1 - 2 * A * lam * x1) * x1p
              + (1 - 2 * A * lam + 2 * A * b2) * y2 * y2p + (1 - 2 * A * lam + 2 * A * b1) * y1 * y1p)) / D
    x1pp = fprime / s6 - hprime * x1 - h * x1p
    y1pp = -np.sqrt(1.5) * b1 * (x1p * y1 + x1 * y1p) - hprime * y1 - h * y1p
    y2pp = -np.sqrt(1.5) * b2 * (x1p * y2 + x1 * y2p) - hprime * y2 - h * y2p
    Orpp = -2 * Orp * (2 + h) - 2 * Or * hprime
    D2 = 1 + 2 * A * (3 * A - lam)
    fpprime = -3 / D2 * (-A * Orpp + 6 * A * (2 * A * lam - 1) * x1p**2 + 2 * (3 * A - b1) * (y1p**2 + y1 * y1pp)
                         + 2 * (3 * A - b2) * (y2p**2 + y2 * y2pp))
    hpprime = -1 / D2 * (0.5 * Orpp * (2 * A * lam - 1) + 3 * (-(1 - 2 * A * lam) ** 2 * x1p**2
                         + (2 * A * lam - 1) * (s6 * A + x1 - 2 * A * lam * x1) * x1pp
                         + (1 - 2 * A * lam + 2 * A * b1) * (y1p**2 + y1 * y1pp)
                         + (1 + 2 * A * lam + 2 * A * b2) * (y2p**2 + y2 * y2pp)))
    x1ppp = np.sqrt(1 / 6) * fpprime - hpprime * x1 - h * x1pp - 2 * hprime * x1p
    adotoa = H
    Ldot = -2 * adotoa**3 * (2 * s6 * A * x1 * h**2 + 3 * h * ((-1 + 2 * A * lam) * x1**2 + y1**2 + y2**2
                             + s6 * A * x1p) + x1 * (s6 * A * hprime + 3 * (-1 + 2 * A * lam) * x1p)
                             + 3 * y1 * y1p + 3 * y2 * y2p + s6 * A * x1pp)
    Lddot = -2 * adotoa**4 * (6 * s6 * A * x1 * h**3 + 3 * hprime * (y1**2 + y2**2 + x1**2 * (2 * A * lam - 1))
                              + 4 * s6 * A * x1p * hprime + 3 * (2 * A * lam - 1) + A * s6 * (x1pp + x1ppp)
                              + h**2 * (2 * s6 * A * x1 + 9 * (y1**2 + y2**2) + 9 * (2 * A * lam - 1) * x1**2
                                        + 11 * s6 * x1p * A)
                              + x1p * (s6 * A * hprime + 3 * (x1p + x1pp) * (2 * A * lam - 1) + s6 * A * hpprime)
                              + 3 * (y1 * (y1p + y1pp) + y1p**2 + y2 * (y2p + y2pp) + y2p**2)
                              + h * (3 * x1**2 * (2 * A * lam - 1) + x1 * (15 * (2 * A * lam - 1) + 7 * s6 * A * hprime)
                                     + 3 * s6 * a * (x1p + 2 * x1pp) + 3 * (y1**2 + y2**2 + 5 * y1 * y1p + 5 * y2 * y2p)))
    Hdot = adotoa**2 * (1 + h)
    ht = adotoa * hprime
    htt = adotoa**2 * hpprime + Hdot * hprime
    Hdd = 2 * adotoa * Hdot * (1 + h) + adotoa**2 * ht
    Hddd = adotoa**2 * htt + 4 * adotoa * Hdot * ht + 2 * adotoa * Hdd * (1 + h) + 2 * Hdot**2 * (1 + h)
    G2P = -2 * s6 * A * adotoa * (h * x1 + x1p) / a**2
    G2PP = -2 * s6 * A * adotoa / a**3 * (-x1p - h * x1 + 2 * x1p * h + hprime * x1 + h**2 * x1 + x1pp)
    G2PPP = -2 * s6 * A * adotoa / a**4 * (x1 * (2 * h - 3 * hprime - 3 * h**2 + 3 * h * hprime + h**3 + hpprime)
                                          + x1p * (2 - 6 * h + 3 * h * 2 + 3 * hprime) + x1pp * (3 * h - 3) + x1ppp)
    return {"f": f, "h": h, "x1p": x1p, "hprime": hprime, "hpprime": hpprime, "fprime": fprime, "fpprime": fpprime,
            "x1pp": x1pp, "x1ppp": x1ppp, "EFTLambdadot": Ldot, "EFTLambdadotdot": Lddot, "Hdotdot": Hdd,
            "Hdotdotdot": Hddd, "G2P": G2P, "G2PP": G2PP, "G2PPP": G2PPP}


def fixed(x1, y1, y2, Or, H, a):
    """Proposed corrected expressions (same inputs as the Fortran): fpprime with the x1pp terms, hpprime with the
    overall sign and the y2 coefficient (1 - 2 A lambda + 2 A beta2) fixed, EFTLambdadotdot = -2 H^4 ((1+3h) K + K'),
    EFTGamma2PPP with 3 h**2 (not 3 h*2) in the x1p coefficient."""
    c = code(x1, y1, y2, Or, H, a)
    h, hprime, x1p, x1pp = c["h"], c["hprime"], c["x1p"], c["x1pp"]
    y1p = -np.sqrt(1.5) * b1 * x1 * y1 - y1 * h
    y2p = -np.sqrt(1.5) * b2 * x1 * y2 - y2 * h
    y1pp = -np.sqrt(1.5) * b1 * (x1p * y1 + x1 * y1p) - hprime * y1 - h * y1p
    y2pp = -np.sqrt(1.5) * b2 * (x1p * y2 + x1 * y2p) - hprime * y2 - h * y2p
    Orp = -2 * Or * (2 + h)
    Orpp = -2 * Orp * (2 + h) - 2 * Or * hprime
    D2 = 1 + 2 * A * (3 * A - lam)
    fpprime = 3 / D2 * (A * Orpp + (2 * A * lam - 1) * (s6 - 6 * A * x1) * x1pp + 6 * A * (1 - 2 * A * lam) * x1p**2
                        + 2 * (b1 - 3 * A) * (y1p**2 + y1 * y1pp) + 2 * (b2 - 3 * A) * (y2p**2 + y2 * y2pp))
    hpprime = 1 / D2 * (0.5 * Orpp * (2 * A * lam - 1) + 3 * (-(1 - 2 * A * lam) ** 2 * x1p**2
                        + (2 * A * lam - 1) * (s6 * A + x1 - 2 * A * lam * x1) * x1pp
                        + (1 - 2 * A * lam + 2 * A * b1) * (y1p**2 + y1 * y1pp)
                        + (1 - 2 * A * lam + 2 * A * b2) * (y2p**2 + y2 * y2pp)))
    x1ppp = np.sqrt(1 / 6) * fpprime - hpprime * x1 - h * x1pp - 2 * hprime * x1p
    K = (2 * s6 * A * x1 * h**2 + 3 * h * ((-1 + 2 * A * lam) * x1**2 + y1**2 + y2**2 + s6 * A * x1p)
         + x1 * (s6 * A * hprime + 3 * (-1 + 2 * A * lam) * x1p) + 3 * y1 * y1p + 3 * y2 * y2p + s6 * A * x1pp)
    Kp = (2 * s6 * A * (x1p * h**2 + 2 * x1 * h * hprime)
          + 3 * hprime * ((-1 + 2 * A * lam) * x1**2 + y1**2 + y2**2 + s6 * A * x1p)
          + 3 * h * (2 * (-1 + 2 * A * lam) * x1 * x1p + 2 * y1 * y1p + 2 * y2 * y2p + s6 * A * x1pp)
          + x1p * (s6 * A * hprime + 3 * (-1 + 2 * A * lam) * x1p) + x1 * (s6 * A * hpprime + 3 * (-1 + 2 * A * lam) * x1pp)
          + 3 * (y1p**2 + y1 * y1pp) + 3 * (y2p**2 + y2 * y2pp) + s6 * A * x1ppp)
    Lddot = -2 * H**4 * ((1 + 3 * h) * K + Kp)
    Hdot = H**2 * (1 + h)
    ht = H * hprime
    htt = H**2 * hpprime + Hdot * hprime
    Hdd = 2 * H * Hdot * (1 + h) + H**2 * ht
    Hddd = H**2 * htt + 4 * H * Hdot * ht + 2 * H * Hdd * (1 + h) + 2 * Hdot**2 * (1 + h)
    G2PPP = -2 * s6 * A * H / a**4 * (x1 * (2 * h - 3 * hprime - 3 * h**2 + 3 * h * hprime + h**3 + hpprime)
                                     + x1p * (2 - 6 * h + 3 * h**2 + 3 * hprime) + x1pp * (3 * h - 3) + x1ppp)
    c.update({"fpprime": fpprime, "hpprime": hpprime, "x1ppp": x1ppp, "EFTLambdadotdot": Lddot, "Hdotdotdot": Hddd,
              "G2PPP": G2PPP})
    return c


def main() -> None:
    rng = np.random.default_rng(1)
    pts = [(rng.uniform(-0.05, 0.05), rng.uniform(0, 0.05), rng.uniform(0, 0.8), rng.uniform(0, 0.9),
            rng.uniform(1e-4, 1e-2), rng.uniform(-6, 0)) for _ in range(40)]
    for label, impl in (("Fortran expression", code), ("proposed fix", fixed)):
        report(label, impl, pts)


def report(label, impl, pts) -> None:
    worst = {}
    for x1, y1, y2, Or, H, n in pts:
        c = impl(x1, y1, y2, Or, H, np.exp(n))
        for k, fn in ex_fn.items():
            e = float(fn(x1, y1, y2, Or, H))
            r = abs(c[k] - e) / max(abs(e), 1e-300)
            worst[k] = max(worst.get(k, 0.0), r)
        for k, fn in g2_fn.items():
            e = float(fn(x1, y1, y2, Or, H, n))
            r = abs(c[k] - e) / max(abs(e), 1e-300)
            worst[k] = max(worst.get(k, 0.0), r)
    print(f"Scaling Cubic Galileon: max relative difference, {label} vs exact derivative along the ODE flow")
    print(f"(40 random points, default parameters A={A}, lambda={lam}, beta1={b1}, beta2={b2}); expected ~1e-12")
    for k in ["h", "x1p", "fprime", "hprime", "x1pp", "fpprime", "hpprime", "x1ppp", "EFTLambdadot", "EFTLambdadotdot",
              "Hdotdot", "Hdotdotdot", "G2P", "G2PP", "G2PPP"]:
        flag = "OK " if worst[k] < 1e-9 else "BAD"
        print(f"  {flag} {k:16s} {worst[k]:.3e}")
    # targeted: y2 = 0 (early times) isolates the overall sign of hpprime
    x1, y1, y2, Or, H = 0.01, 0.02, 0.0, 0.5, 1e-3
    c = impl(x1, y1, y2, Or, H, 1e-3)
    e = float(ex_fn["hpprime"](x1, y1, y2, Or, H))
    print(f"  y2 = 0 point: hpprime {label} = {c['hpprime']:+.6e}, exact = {e:+.6e} (ratio {c['hpprime'] / e:+.6f})")


if __name__ == "__main__":
    main()
