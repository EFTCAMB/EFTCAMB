"""Verifier (t10a-fullmap-hdotdotdot-formula), independent method: exact analytic toy trajectory fed to the COMPILED
A0 routine EFTCAMBFullModelComputeHubbleDer (06p1_abstract_EFTCAMB_full_map.f90) through a small Fortran driver
(t10a_hddd_toy_driver.f90, linked against the frozen A0 camblib.so).

Construction (no CAMB run, no splines, no finite differences): pick analytic a(tau), Omega(a), EFTLambda L(tau).
Then H = a'/a and its exact derivatives Hd, Hdd, Hddd follow from sympy.  The pressure P = gpresm_t is chosen so that
the code's own Hdot equation holds exactly:  P = -2 D Hd - (1+O+2aO'+a^2 O'') H^2 - L,  D = 1+O+a O'/2.
Derivative conventions (checked against the code's radiation expressions for gpresdotm_t, gpresdotdotm_t and against
compute_additional_derivs for EFTLambdadotdot):  Pd = P' - 2HP, Pdd = Pd' - 2HPd, Ld = L' - 2HL, Ldd = Ld' - 2HLd.
Expected: Fortran Hdot == Hd, Hdotdot == Hdd, Hdotdotdot == Hddd to rounding (~1e-12 relative).

Also: a Python transcription of the shipped Hdotdotdot block (validated against the Fortran output) with each of the
four root-cause corrections (i)-(iv) switched on/off, to test the stated attribution.
Usage: python t10a_hddd_toy.py <driver executable> [out.json]
"""

import itertools
import json
import subprocess
import sys

import numpy as np
import sympy as sp

DRIVER = sys.argv[1]
OUT = sys.argv[2] if len(sys.argv) > 2 else None

tau, x = sp.symbols("tau x", positive=True)


def make_case(name, omega_expr, lam=1.0):
    """lam rescales conformal time (tau -> tau/lam): H ~ 1/lam, as in Mpc units when lam ~ 1e3-1e4."""
    s = tau / lam
    a = sp.Rational(3, 10) * s + sp.Rational(1, 5) * s**2 + sp.Rational(1, 20) * s**3
    L = (sp.Rational(3, 10) + sp.Rational(1, 10) * s**2 + sp.Rational(1, 50) * s**3) / lam**2
    Om = [sp.diff(omega_expr, x, k) for k in range(5)]
    Om_t = [o.subs(x, a) for o in Om]
    H = sp.diff(a, tau) / a
    Hd = sp.diff(H, tau)
    Hdd = sp.diff(Hd, tau)
    Hddd = sp.diff(Hdd, tau)
    D = 1 + Om_t[0] + a * Om_t[1] / 2
    P = -2 * D * Hd - (1 + Om_t[0] + 2 * a * Om_t[1] + a**2 * Om_t[2]) * H**2 - L
    Pd = sp.diff(P, tau) - 2 * H * P
    Pdd = sp.diff(Pd, tau) - 2 * H * Pd
    Ld = sp.diff(L, tau) - 2 * H * L
    Ldd = sp.diff(Ld, tau) - 2 * H * Ld
    rows, truth = [], []
    for t in (0.5, 1.0, 1.5, 2.0, 3.0):
        tv = sp.Float(t * lam, 40)
        ev = [sp.N(e.subs(tau, tv), 30) for e in (a, H, *Om_t, P, Pd, Pdd, L, Ld, Ldd)]
        rows.append([float(v) for v in ev])
        truth.append([float(sp.N(e.subs(tau, tv), 30)) for e in (Hd, Hdd, Hddd)])
    return name, np.array(rows), np.array(truth)


def run_driver(rows):
    inp = "\n".join(" ".join(f"{v:.17e}" for v in r) for r in rows) + "\n"
    out = subprocess.run([DRIVER], input=inp, capture_output=True, text=True, check=True).stdout
    return np.array([[float(v) for v in ln.split()] for ln in out.strip().splitlines()])


def shipped_variant(r, Hd, Hdd, fix_i=False, fix_ii=False, fix_iii=False, fix_iv=False):
    """Python transcription of 06p1 lines 141-163 (shipped), with optional corrections:
    (i) N1 = -a H (3 O1 + a O2)        [shipped: 0.5 a O2]
    (ii) '+' instead of '*' between 3.5 a O1 and 1.5 a^2 O2 (both occurrences)
    (iii) bracket multiplying N1 = Hdotdot numerator (Pd + 2HP + Ld + 2HL)
    (iv) -a H^2 (0.5 H^2 + 1.5 Hd) S3   [shipped: -a H^2 (H^2 + 3 H) S3]"""
    a, H, O0, O1, O2, O3, O4, P, Pd, Pdd, L, Ld, Ldd = r
    a2, H2, opo = a * a, H * H, 1 + O0
    denom = opo + 0.5 * a * O1
    N1 = -a * H * (3 * O1 + (a * O2 if fix_i else 0.5 * a * O2))
    T3 = (opo + 3.5 * a * O1 + 1.5 * a2 * O2) if fix_ii else (opo + 3.5 * a * O1 * 1.5 * a2 * O2)
    S3 = 3 * O1 + 4 * a * O2 + a2 * O3
    br_pl = (Pd + 2 * H * P + Ld + 2 * H * L) if fix_iii else (Pdd + 2 * H * Pd - 0.5 * (Ld + 2 * H * L))
    first = (-a * H2 * (0.5 * H2 + 1.5 * Hd) * S3) if fix_iv else (-a * H2 * (H2 + 3 * H) * S3)
    return (N1 / 2 / denom**2 * (-0.5 * a * H**3 * S3 - H * Hd * T3 - 0.5 * br_pl)
            + 1 / denom * (first - 0.5 * a2 * H**4 * (7 * O2 + 6 * a * O3 + a2 * O4) - (Hd**2 + H * Hdd) * T3
                           - 0.5 * a * H2 * Hd * (9 * O1 + 13 * a * O2 + 3 * a2 * O3)
                           - 0.5 * (Pdd + 4 * H * Pd + 2 * Hd * P + 4 * H2 * P)
                           - 0.5 * (Ldd + 4 * H * Ld + 2 * Hd * L + 4 * H2 * L)))


cases = [
    make_case("Omega varying (poly+sin), lam=1", sp.Rational(1, 10) * x**2 - sp.Rational(1, 20) * x**3
              + sp.Rational(1, 20) * sp.sin(2 * x)),
    make_case("Omega varying, lam=3e3 (Mpc-like H~1e-4)", sp.Rational(1, 10) * x**2 - sp.Rational(1, 20) * x**3
              + sp.Rational(1, 20) * sp.sin(2 * x), lam=3000),
    make_case("Omega small power law 0.01 a^2, lam=3e3", sp.Rational(1, 100) * x**2, lam=3000),
    make_case("Omega constant 0.1, lam=1", sp.Rational(1, 10) + 0 * x),
    make_case("Omega = 0, lam=3e3", 0 * x, lam=3000),
]

res = {}
print(f"driver: {DRIVER}")
for name, rows, truth in cases:
    f = run_driver(rows)
    rel = lambda u, v: float(np.max(np.abs(u - v) / np.abs(v)))  # noqa: E731
    r_hd, r_hdd, r_hddd, r_prop = rel(f[:, 0], truth[:, 0]), rel(f[:, 1], truth[:, 1]), rel(f[:, 2], truth[:, 2]), \
        rel(f[:, 3], truth[:, 2])
    ship_py = np.array([shipped_variant(r, f[i, 0], f[i, 1]) for i, r in enumerate(rows)])
    r_transcr = float(np.max(np.abs(ship_py - f[:, 2]) / np.abs(f[:, 2])))
    print(f"\n== {name}")
    print(f"  max rel |Hdot_F - Hd_exact|            = {r_hd:.2e}")
    print(f"  max rel |Hdotdot_F - Hdd_exact|        = {r_hdd:.2e}")
    print(f"  max rel |Hdotdotdot_F - Hddd_exact|    = {r_hddd:.2e}")
    print(f"  max rel |proposed_fix - Hddd_exact|    = {r_prop:.2e}")
    print(f"  python transcription of shipped vs Fortran Hdotdotdot: {r_transcr:.2e}")
    for i in (1, 3):
        print(f"   tau={rows[i][0]:.3g}(a)  Hddd exact {truth[i, 2]: .6e}  Fortran {f[i, 2]: .6e}  proposed {f[i, 3]: .6e}")
    attrib = {}
    for k in range(5):
        for combo in itertools.combinations(("i", "ii", "iii", "iv"), k):
            kw = {f"fix_{c}": True for c in combo}
            v = np.array([shipped_variant(r, f[i, 0], f[i, 1], **kw) for i, r in enumerate(rows)])
            attrib["+".join(combo) or "none"] = float(np.max(np.abs(v - truth[:, 2]) / np.abs(truth[:, 2])))
    print("  attribution (max rel error vs exact with listed corrections applied):")
    print("   " + ", ".join(f"{k}:{v:.1e}" for k, v in attrib.items()))
    res[name] = dict(Hdot=r_hd, Hdotdot=r_hdd, Hdotdotdot=r_hddd, proposed=r_prop, transcription=r_transcr,
                     attribution=attrib)
if OUT:
    json.dump(res, open(OUT, "w"), indent=1)
