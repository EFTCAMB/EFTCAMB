"""Exact (rational) symbolic evaluation of EFTCAMB's EFT_mu1/EFT_mu2 for constant-Gamma1 pure EFT on LCDM.

Uses the automatic transcription of EFTCAMBModelComputeMassFactors from gamma1_mu2_formula.py.
Units H = 1 (cosmic-time Hubble rate at the evaluation time): g = 4 H0^2 Gamma1 / H^2, h1 = Hdot/H^2,
h2 = Hddot/H^3, r = rho_dust/H^2 (8 pi G rho), q = rho_rad/H^2. Output: mu/H^2.
"""

import sys
from pathlib import Path

import sympy as sp

sys.path.insert(0, str(Path(__file__).resolve().parent))
import gamma1_mu2_formula as gf  # noqa: E402

a, g, h0 = sp.symbols("a g h0", positive=True)
h1, h2 = sp.symbols("h1 h2", real=True)
r, q = sp.symbols("r q", nonnegative=True)


def sqrt_custom(x):
    v = sp.sqrt(sp.factor(x))
    # the only sqrt in the formula is sqrt(F2^4/(3F2^2+2F1F3)^2); here 3F2^2+2F1F3 = 4g-12 < 0 for g < 3
    return v.replace(sp.Abs, lambda arg: -arg if sp.simplify(arg.subs(g, 0)) < 0 else arg)


def build():
    env = {"a": a, "h0": h0}
    for f in gf.CACHE_FIELDS:
        env["C_" + f] = sp.Integer(0)
    env["C_adotoa"] = a
    env["C_Hdot"] = a**2 * (h1 + 1)
    env["C_Hdotdot"] = a**3 * (h2 + 4 * h1 + 2)
    env["C_grhob_t"] = r * a**2
    env["C_EFTGamma1V"] = g / (4 * h0**2)
    code = gf.to_python(gf.fortran_body(), rational=True)
    ns = dict(env)
    ns["Rat"] = sp.Rational
    ns["sqrt"] = sqrt_custom
    exec(code, ns)
    return ns


def main() -> None:
    ns = build()
    for k in ("F1", "F2", "F3", "F1dot", "F2dot", "F3dot", "F2dotdot", "F3dotdot", "rho_d", "rho_d_dot"):
        print(f"  {k:9s} = {sp.simplify(ns[k])}")
    print("  3F2^2+F1F3 =", sp.simplify(3 * ns["F2"] ** 2 + ns["F1"] * ns["F3"]),
          "  3F2^2+2F1F3 =", sp.simplify(3 * ns["F2"] ** 2 + 2 * ns["F1"] * ns["F3"]))
    lcdm = {h2: sp.Rational(3, 2) * r + sp.Rational(8, 3) * q, h1: -r / 2 - sp.Rational(2, 3) * q}
    out = {}
    for mu in ("EFT_mu1", "EFT_mu2"):
        e = sp.factor(sp.together(ns["C_" + mu]))
        print(f"\n{mu} (generic background) = {e}")
        num, den = sp.fraction(e)
        print(f"{mu}: denominator = {sp.factor(den)}")
        print(f"{mu}: numerator at g=0 = {sp.factor(num.subs(g, 0))}")
        e_l = sp.factor(sp.together(e.subs(lcdm)))
        print(f"{mu} (LCDM: h1=-r/2-2q/3, h2=3r/2+8q/3) = {e_l}")
        ser = sp.series(e_l, g, 0, 1).removeO()
        print(f"{mu} (LCDM) Laurent series in g up to g^0: {sp.expand(ser)}")
        out[mu] = e_l
        e_m = sp.factor(e_l.subs(q, 0))
        print(f"{mu} (LCDM, no radiation) = {e_m}")
    return out


if __name__ == "__main__":
    main()
