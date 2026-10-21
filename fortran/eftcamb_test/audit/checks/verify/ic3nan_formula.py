"""Offline check of the EFT_IC_type=3 formulas, extracted verbatim from the Fortran source (no transcription):
 - float64 evaluation at Omega0 = Gamma30 = 0 (IEEE, like the release build)  -> NaN?
 - sympy limits along rays Gamma30 = alpha*Omega0 -> 0: removable (direction independent) or not?
 - compare with case 2 (Omega0 -> 0) and GR (IC1) leading terms.
"""

import re

import numpy as np
import sympy as sp

SRC = "/Users/marco/IW/eftcamb_update/EFTCAMB/fortran"


def statements(block):
    block = re.sub(r"!.*", "", block)
    block = re.sub(r"&\s*\n\s*", "", block)
    out = {}
    for line in block.splitlines():
        line = line.strip()
        m = re.match(r"(initv\(1, (i_\w+)\)|y\(EV%w_ix( \+ 1)?\))\s*=\s*(.+)$", line)
        if m:
            key = m.group(2) or ("pidot" if m.group(3) else "pi")
            out[key] = m.group(4)
    return out


def to_py(expr):
    e = expr.replace("_dl", "")
    for a, b in [("CP%eft_par_cache%h0_Mpc", "h0"), ("EV%Kf(1)", "Kf"), ("EV%eft_cache%tau", "tau"),
                 ("EV%eft_cache%EFTOmegaV", "Omega0")]:
        e = e.replace(a, b)
    e = re.sub(r"initv\(1, (i_\w+)\)", r"\1", e)
    return e


eq = open(SRC + "/equations.f90").read()
blk = eq[eq.index("CASE 3: Modified Gravity (Constant Gamma + Constant Omega)"):]
blk = blk[:blk.index("case default")]
s3 = statements(blk)
ic = open(SRC + "/eftcamb/09_EFTCAMB_IC.f90").read()
b9 = ic[ic.index("--- CASE 3: Constant Gamma + Constant Omega ---"):]
b9 = b9[:b9.index("case default")]
s9 = statements(b9)
b2 = eq[eq.index("CASE 2: Modified Gravity (Constant Omega)"):eq.index("CASE 3: Modified Gravity (Constant Gamma")]
s2 = statements(b2)
print("extracted equations.f90 case 3:", list(s3), "| 09 case 3:", list(s9), "| equations.f90 case 2:", list(s2))

O, G, G2, x, tau, h0, Kf, Rv, chi, al = sp.symbols("Omega0 Gamma30 Gamma20 x tau h0 Kf Rv chi alpha")
k2 = sp.Symbol("k2")
ns = {"Omega0": O, "Gamma30": G, "Gamma20": G2, "x": x, "x2": x**2, "x3": x**3, "tau": tau, "h0": h0, "Kf": Kf,
      "Rv": Rv, "chi": sp.Integer(1), "k2": k2}


def sym(e, extra=None):
    d = dict(ns)
    d.update(extra or {})
    return sp.sympify(to_py(e), locals=d)


# --- equations.f90 case 3 adiabatic ICs, in dependency order
vals3 = {}
for key in ["i_clxg", "i_clxr", "i_clxb", "i_clxc", "i_qg", "i_qr", "i_vb", "i_pir", "i_aj3r", "i_eta"]:
    vals3[key] = sym(s3[key], {k: v for k, v in vals3.items()})
vals2 = {}
for key in ["i_clxg", "i_clxr", "i_clxb", "i_clxc", "i_qg", "i_qr", "i_vb", "i_pir", "i_aj3r", "i_eta"]:
    vals2[key] = sym(s2[key], {k: v for k, v in vals2.items()})
pi3 = sym(s9["pi"])
pid3 = sym(s9["pidot"])

num = {x: 1e-3, tau: 0.5, h0: 2.2e-4, Kf: 1.0, Rv: 0.4089, G2: 0.0, k2: 4e-6}
print("\n[float64 IEEE evaluation at Omega0 = Gamma30 = 0]  (Gamma20 = 0, x = 1e-3, tau = 0.5)")
with np.errstate(all="ignore"):
    for name, ex in list(vals3.items()) + [("pi(09)", pi3), ("pidot(09)", pid3)]:
        f = sp.lambdify((O, G), ex.subs(num), "numpy")
        v0 = f(np.float64(0.0), np.float64(0.0))
        v1 = f(np.float64(1e-12), np.float64(0.0))
        print(f"  case3 {name:10s} at (0,0) = {v0!r:>24}   at (Om=1e-12, G3=0) = {v1: .6e}")

print("\n[limits along rays Gamma30 = alpha*Omega0, Omega0 -> 0, leading order in x]  vs case 2 (Om->0) and GR")
gr = {"i_clxg": -Kf / 3 * x**2, "i_eta": -2 * Kf}
for name in ["i_clxg", "i_qr", "i_pir", "i_eta"]:
    ex = vals3[name].subs(G2, 0)
    ray = sp.simplify(sp.limit(ex.subs(G, al * O), O, 0))
    c2 = sp.simplify(sp.limit(vals2[name], O, 0))
    print(f"  {name:7s} case3 ray limit: {sp.factor(ray)}   | case2 Om->0: {sp.factor(c2)}   | depends on alpha: "
          f"{ray.has(al)}")
for name, ex in [("pi(09)", pi3), ("pidot(09)", pid3)]:
    ray = sp.simplify(sp.limit(ex.subs(G2, 0).subs(G, al * O), O, 0))
    print(f"  {name:9s} case3 ray limit: {sp.factor(ray)}   | depends on alpha: {ray.has(al)}  "
          f"(alpha=0: {sp.factor(ray.subs(al, 0))}, alpha->inf: {sp.factor(sp.limit(ray, al, sp.oo))})")
# case 2 pi IC (09) at Omega -> 0 for reference
c2pi = sym(s9_2) if (s9_2 := None) else None
b92 = ic[ic.index("--- CASE 2: Constant Omega ---"):ic.index("--- CASE 3: Constant Gamma")]
s92 = statements(b92)
print("  09 case 2 pi at Omega0 = 0:", sp.factor(sym(s92["pi"]).subs(O, 0)), " pidot:", sp.factor(sym(s92["pidot"]).subs(O, 0)))
