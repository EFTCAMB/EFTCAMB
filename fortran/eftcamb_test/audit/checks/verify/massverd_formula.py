"""Verify (mass-stability-verdict-unreliable): independent symbolic test of the EFT_mu1/EFT_mu2 expressions.

The Fortran block of EFTCAMBModelComputeMassFactors (06_abstract_EFTCAMB_model.f90) is transcribed automatically to
sympy. For a canonical (minimally coupled) quintessence field (Omega = gamma_i = 0) the inputs are expressed through
cosmic-time physical quantities (8 pi G = m0^2 = 1):
    3 H^2 = rho_m + rho_r + rho_phi,  dH/dt = -(rho_m + 4/3 rho_r)/2 - c,  c = (rho_phi + p_phi)/2,
and the code's conformal-time cache entries (adotoa, Hdot, Hdotdot, EFTc, EFTcdot, EFTcdotdot) are built from them
exactly, so all derivatives are mutually consistent.  The near-LCDM limit is c = eps*c0, dc/dt = eps*c1,
d2c/dt2 = eps*c2, eps -> 0.  A physical mass eigenvalue (in units of H^2) must stay finite as eps -> 0.

Prints the leading power of eps in mu1, mu2 and its coefficient with and without radiation.

Run: python massverd_formula.py [path to 06_abstract_EFTCAMB_model.f90]
"""

import re
import sys

import sympy as sp

SRC = sys.argv[1] if len(sys.argv) > 1 else (
    "/private/tmp/claude-501/-Users-marco-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/"
    "scratchpad/audit/trees/A0/fortran/eftcamb/06_abstract_EFTCAMB_model.f90")


def transcribe(path):
    text = open(path).read()
    i = text.index("subroutine EFTCAMBModelComputeMassFactors")
    j = text.index("end subroutine EFTCAMBModelComputeMassFactors")
    body = text[i:j]
    body = body[body.index("rho_d = "):]
    lines = [ln.split("!")[0].rstrip() for ln in body.splitlines()]
    joined = ""
    for ln in lines:
        s = ln.strip()
        if not s:
            continue
        if joined.endswith("&"):
            joined = joined[:-1] + " " + s
        else:
            joined += "\n" + s
    stmts = [s for s in joined.split("\n") if s.strip()]
    out = []
    for s in stmts:
        s = s.replace("eft_cache%", "").replace("eft_par_cache%", "")
        s = re.sub(r"(\d+)\._dl", r"\1", s)
        s = re.sub(r"(\d+)\.(?!\d)", r"\1", s)
        s = s.replace("sqrt", "sp.sqrt")
        out.append(s)
    return out


stmts = transcribe(SRC)
names = ["a", "adotoa", "Hdot", "Hdotdot", "EFTc", "EFTcdot", "EFTcdotdot", "grhob_t", "grhoc_t", "h0_Mpc"]
zero = ["EFTOmegaV", "EFTOmegaP", "EFTOmegaPP", "EFTOmegaPPP", "EFTGamma1V", "EFTGamma1P", "EFTGamma1PP",
        "EFTGamma2V", "EFTGamma2P", "EFTGamma2PP", "EFTGamma3V", "EFTGamma3P", "EFTGamma3PP", "EFTGamma4V",
        "EFTGamma4P", "EFTGamma4PP"]

a, H, rm, rr, eps, c0, c1, c2 = sp.symbols("a H rho_m rho_r eps c0 c1 c2", positive=True)
c, ct, ctt = eps * c0, eps * c1, eps * c2
Hd = -(rm + sp.Rational(4, 3) * rr) / 2 - c                       # dH/dt
Hdd = (3 * H * rm + sp.Rational(16, 3) * H * rr) / 2 - ct         # d2H/dt2
env = {"sp": sp}
env.update({z: 0 for z in zero})
env.update({
    "a": a,
    "adotoa": a * H,
    "Hdot": a**2 * (H**2 + Hd),
    "Hdotdot": a**3 * (2 * H**3 + 4 * H * Hd + Hdd),
    "EFTc": c * a**2,
    "EFTcdot": a * ct * a**2,
    "EFTcdotdot": a**2 * (H * ct + ctt) * a**2,
    "grhob_t": rm * a**2,
    "grhoc_t": 0,
    "h0_Mpc": sp.Symbol("h0"),
})
for s in stmts:
    lhs, rhs = s.split("=", 1)
    lhs = lhs.strip()
    env[lhs] = eval(rhs, env)  # noqa: S307

mu1, mu2 = env["EFT_mu1"], env["EFT_mu2"]
# consistency: F3dot computed by the code equals d/dt of F3 (cosmic time) for these inputs
print("F1, F2, F3 =", sp.simplify(env["F1"]), sp.simplify(env["F2"]), sp.simplify(env["F3"]))


def leading(expr, subs):
    e = expr.subs(subs)
    e = e.replace(sp.Abs, lambda x: x if x.subs(eps, 0) > 0 else -x)
    e = sp.together(sp.nsimplify(e))
    num, den = sp.fraction(e)
    pn, pd = sp.Poly(sp.expand(num), eps), sp.Poly(sp.expand(den), eps)
    ln = min(m[0] for m in pn.monoms())
    ld = min(m[0] for m in pd.monoms())
    cn = pn.coeff_monomial(eps**ln)
    cd = pd.coeff_monomial(eps**ld)
    return ln - ld, sp.simplify(cn / cd)


# numbers: units H = 1 at a = 1 (physical), Omega_m = 0.3, c0..c2 generic O(1) numbers
for label, rrv in (("dust only (rho_r = 0)", 0), ("with radiation rho_r = 3e-3", sp.Rational(3, 1000))):
    subs = {a: 1, H: 1, rm: sp.Rational(9, 10), rr: rrv, c0: sp.Rational(1, 7), c1: sp.Rational(-2, 9),
            c2: sp.Rational(3, 11)}
    for nm, mu in (("mu1", mu1), ("mu2", mu2)):
        p, coef = leading(mu, subs)
        print(f"{label:30s} {nm}/H^2 ~ {float(coef):+.6g} * eps^{p}")
    # numeric values at finite eps to show the size
    for ev in (1e-1, 1e-2, 1e-3, 1e-4):
        s2 = dict(subs)
        s2[eps] = sp.Rational(ev).limit_denominator(10**8)
        print(f"    eps={ev:g}: mu1/H^2 = {float(mu1.subs(s2)):+.6g}, mu2/H^2 = {float(mu2.subs(s2)):+.6g}")

# general structure of the leading coefficient (symbolic in rho_r, other numbers fixed)
subs = {a: 1, H: 1, rm: sp.Rational(9, 10), c0: sp.Rational(1, 7), c1: sp.Rational(-2, 9), c2: sp.Rational(3, 11)}
for nm, mu in (("mu1", mu1), ("mu2", mu2)):
    p, coef = leading(mu, subs)
    print(f"symbolic in rho_r: {nm} leading power eps^{p}, coefficient = {sp.factor(coef)}")

# a-independence (physical inputs fixed, a = 1/2 instead of 1): the eigenvalues must not change
subs = {a: sp.Rational(1, 2), H: 1, rm: sp.Rational(9, 10), rr: 0, c0: sp.Rational(1, 7), c1: sp.Rational(-2, 9),
        c2: sp.Rational(3, 11), eps: sp.Rational(1, 100)}
subs1 = dict(subs)
subs1[a] = 1
print("a-independence: mu1(a=1/2) =", float(mu1.subs(subs)), " mu1(a=1) =", float(mu1.subs(subs1)))

# generic inconsistency of the background inputs with the dust-only relations: dH/dt -> dH/dt + d1, d2H/dt2 -> + d2
d1, d2 = sp.symbols("d1 d2")
env2 = dict(env)
Hd2 = -rm / 2 - c + d1
Hdd2 = 3 * H * rm / 2 - ct + d2
env2.update({"Hdot": a**2 * (H**2 + Hd2), "Hdotdot": a**3 * (2 * H**3 + 4 * H * Hd2 + Hdd2)})
for s in transcribe(SRC):
    lhs, rhs = s.split("=", 1)
    env2[lhs.strip()] = eval(rhs, env2)  # noqa: S307
subs = {a: 1, H: 1, rm: sp.Rational(9, 10), c0: sp.Rational(1, 7), c1: sp.Rational(-2, 9), c2: sp.Rational(3, 11)}
for nm in ("EFT_mu1", "EFT_mu2"):
    p, coef = leading(env2[nm], subs)
    print(f"background residuals d1 (dH/dt), d2 (d2H/dt2): {nm} ~ eps^{p} * [{sp.factor(coef)}]")

# EFTcPP must be d^2(EFTc)/da^2 (F3dot uses a*EFTcP with EFTcP = dEFTc/da; F3dotdot uses 2 a^2 H^2 EFTcPP / a^4).
# The code has EFTcPP = (EFTcdotdot/adotoa + (3 - Hdot/adotoa)*EFTcdot + 2*adotoa*EFTc)/(a^2*adotoa^2).
# Exact: d^2 EFTc/da^2 = (EFTcdotdot + (3*adotoa - Hdot/adotoa)*EFTcdot + 2*adotoa^2*EFTc)/(a^2*adotoa^2).
tt = sp.Symbol("t")
cfun = sp.Function("cphys")(tt)
afun = sp.Function("afun")(tt)
# conformal-time quantities written via cosmic time t: d/dtau = a d/dt
dtau = lambda f: afun * sp.diff(f, tt)  # noqa: E731
EFTc_t = cfun * afun**2
calH = dtau(afun) / afun
EFTcdot_t = dtau(cfun) * afun**2
EFTcdotdot_t = dtau(dtau(cfun)) * afun**2
Hdot_t = dtau(calH)
dda = lambda f: dtau(f) / dtau(afun)  # noqa: E731  d/da
exact = dda(dda(EFTc_t))
code = (EFTcdotdot_t / calH + (3 - Hdot_t / calH) * EFTcdot_t + 2 * calH * EFTc_t) / (afun**2 * calH**2)
fixed = (EFTcdotdot_t + (3 * calH - Hdot_t / calH) * EFTcdot_t + 2 * calH**2 * EFTc_t) / (afun**2 * calH**2)
print("EFTcPP check: exact - fixed =", sp.simplify(exact - fixed), "; exact - code == 0 ?",
      sp.simplify(exact - code) == 0)

stmts_fix = [("EFTcPP = (EFTcdotdot + (3*adotoa - Hdot/adotoa)*EFTcdot + 2*adotoa**2*EFTc)/(a**2*adotoa**2)"
              if s.startswith("EFTcPP") else s) for s in transcribe(SRC)]
envf = dict(env)
for s in stmts_fix:
    lhs, rhs = s.split("=", 1)
    envf[lhs.strip()] = eval(rhs, envf)  # noqa: S307
subs = {a: sp.Rational(1, 2), H: 1, rm: sp.Rational(9, 10), rr: 0, c0: sp.Rational(1, 7), c1: sp.Rational(-2, 9),
        c2: sp.Rational(3, 11), eps: sp.Rational(1, 100)}
subs1 = dict(subs)
subs1[a] = 1
print("with EFTcPP fixed, a-independence: mu1(a=1/2) =", float(envf["EFT_mu1"].subs(subs)), " mu1(a=1) =",
      float(envf["EFT_mu1"].subs(subs1)), "; mu2:", float(envf["EFT_mu2"].subs(subs)),
      float(envf["EFT_mu2"].subs(subs1)))
for label, rrv in (("dust only", 0), ("rho_r=3e-3", sp.Rational(3, 1000))):
    s = {a: 1, H: 1, rm: sp.Rational(9, 10), rr: rrv, c0: sp.Rational(1, 7), c1: sp.Rational(-2, 9),
         c2: sp.Rational(3, 11)}
    for nm in ("EFT_mu1", "EFT_mu2"):
        p, coef = leading(envf[nm], s)
        print(f"EFTcPP fixed, {label}: {nm}/H^2 ~ {float(coef):+.6g} eps^{p}")
