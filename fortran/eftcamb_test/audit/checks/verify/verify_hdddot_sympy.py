"""Verifier (independent, symbolic): is the full-mapping Hdotdotdot in 06p1_abstract_EFTCAMB_full_map.f90 the
conformal-time derivative of the shipped Hdotdot?

Conventions taken from the code (checked in 09_EFTCAMB_stability.f90 / 06p1 compute_additional_derivs):
  a' = a H,  H' = Hd, Hd' = Hdd, Hdd' = Hddd
  gpresm_t = P = a^2 p, gpresdotm_t = Pd = a^2 pdot, gpresdotdotm_t = Pdd = a^2 pddot  ->  P' = 2HP + Pd, Pd' = 2HPd + Pdd
  EFTLambda = L = a^2 Lambda, EFTLambdadot = Ld, EFTLambdadotdot = Ldd                ->  L' = 2HL + Ld, Ld' = 2HLd + Ldd
  EFTOmega(a): dOmega/dtau = a H Omega'
The shipped expressions are transcribed verbatim (with the Fortran operators) below.
"""

import sympy as sp

tau = sp.Symbol("tau")
a = sp.Function("a")(tau)
H = sp.Function("H")(tau)
Hd = sp.Function("Hd")(tau)
Hdd = sp.Function("Hdd")(tau)
Pm = sp.Function("P")(tau)
Pd = sp.Function("Pd")(tau)
Pdd = sp.Function("Pdd")(tau)
L = sp.Function("L")(tau)
Ld = sp.Function("Ld")(tau)
Ldd = sp.Function("Ldd")(tau)
x = sp.Symbol("x")
Om = sp.Function("Om")
O0, O1, O2, O3, O4, O5 = [sp.diff(Om(x), x, k).subs(x, a) if k else Om(a) for k in range(6)]

# symbols for "on-shell" evaluation
sa, sH, sHd, sHdd, sHddd, sP, sPd, sPdd, sPddd, sL, sLd, sLdd, sLddd = sp.symbols(
    "a H Hd Hdd Hddd P Pd Pdd Pddd L Ld Ldd Lddd")
so = sp.symbols("o0:6")
rules = {}


def to_symbols(e):
    """replace time derivatives using the chain rules above, then functions by symbols"""
    e = e.doit()
    reps = {
        sp.Derivative(a, tau): a * H,
        sp.Derivative(H, tau): Hd,
        sp.Derivative(Hd, tau): Hdd,
        sp.Derivative(Hdd, tau): sp.Symbol("Hddd"),
        sp.Derivative(Pm, tau): 2 * H * Pm + Pd,
        sp.Derivative(Pd, tau): 2 * H * Pd + Pdd,
        sp.Derivative(Pdd, tau): sp.Symbol("Pddd"),
        sp.Derivative(L, tau): 2 * H * L + Ld,
        sp.Derivative(Ld, tau): 2 * H * Ld + Ldd,
        sp.Derivative(Ldd, tau): sp.Symbol("Lddd"),
    }
    # Omega derivatives: d/dtau Om^(k)(a) = Om^(k+1)(a) a H, handled by sympy chain rule automatically through a(tau)
    for _ in range(3):
        e = e.subs(reps).doit()
    fun = {a: sa, H: sH, Hd: sHd, Hdd: sHdd, Pm: sP, Pd: sPd, Pdd: sPdd, L: sL, Ld: sLd, Ldd: sLdd}
    # Omega derivatives in terms of symbols
    for k in range(5, -1, -1):
        e = e.subs(sp.diff(Om(x), x, k).subs(x, a) if k else Om(a), so[k])
    e = e.subs(fun)
    # leftover Subs/Derivative objects of Om w.r.t. a(tau)
    return sp.expand(e)


a2 = a * a
opo = 1 + O0
adotoa2 = H**2
# ---- shipped Hdot, Hdotdot (06p1 lines 129-139) ----
Hdot_s = 1 / (opo + sp.Rational(1, 2) * a * O1) * (
    -sp.Rational(1, 2) * (opo + 2 * a * O1 + a2 * O2) * H**2 - sp.Rational(1, 2) * Pm - sp.Rational(1, 2) * L)
Hdotdot_s = 1 / (opo + sp.Rational(1, 2) * a * O1) * (
    -sp.Rational(1, 2) * a * H**3 * (3 * O1 + 4 * a * O2 + a2 * O3)
    - H * Hd * (opo + sp.Rational(7, 2) * a * O1 + sp.Rational(3, 2) * a2 * O2)
    - sp.Rational(1, 2) * (Pd + 2 * H * Pm + Ld + 2 * H * L))
# ---- shipped Hdotdotdot (06p1 lines 141-163), verbatim ----
denom = opo + sp.Rational(1, 2) * a * O1
N1 = -a * H * (3 * O1 + sp.Rational(1, 2) * a * O2)
Hdddot_s = (N1 / 2 / denom**2 * (
    -sp.Rational(1, 2) * a * H**3 * (3 * O1 + 4 * a * O2 + a2 * O3)
    - H * Hd * (opo + sp.Rational(7, 2) * a * O1 * sp.Rational(3, 2) * a2 * O2)
    - sp.Rational(1, 2) * (Pdd + 2 * H * Pd - sp.Rational(1, 2) * (Ld + 2 * H * L)))
    + 1 / denom * (-a * adotoa2 * (adotoa2 + 3 * H) * (3 * O1 + 4 * a * O2 + a2 * O3)
                   - sp.Rational(1, 2) * a2 * H**4 * (7 * O2 + 6 * a * O3 + a2 * O4)
                   - (Hd**2 + H * Hdd) * (opo + sp.Rational(7, 2) * a * O1 * sp.Rational(3, 2) * a2 * O2)
                   - sp.Rational(1, 2) * a * adotoa2 * Hd * (9 * O1 + 13 * a * O2 + 3 * a2 * O3)
                   - sp.Rational(1, 2) * (Pdd + 4 * H * Pd + 2 * Hd * Pm + 4 * adotoa2 * Pm)
                   - sp.Rational(1, 2) * (Ldd + 4 * H * Ld + 2 * Hd * L + 4 * adotoa2 * L)))

# 1) control: Hdotdot_s == d(Hdot_s)/dtau on shell (Hd -> Hdot_s)
dHdot = to_symbols(sp.diff(Hdot_s, tau))
c1 = sp.simplify((dHdot - to_symbols(Hdotdot_s)).subs(sHd, to_symbols(Hdot_s)))
print("control  d(Hdot)/dtau - Hdotdot (on shell)          :", c1)

# 2) truth: d(Hdotdot_s)/dtau, on shell (Hdd -> Hdotdot_s, Hd -> Hdot_s)
true3 = to_symbols(sp.diff(Hdotdot_s, tau))
ship3 = to_symbols(Hdddot_s)
onshell = {sHdd: to_symbols(Hdotdot_s)}
diff_ship = sp.simplify(sp.together((true3 - ship3).subs(onshell)))
print("shipped  d(Hdotdot)/dtau - Hdotdotdot (on shell) == 0 ?", diff_ship == 0)

# 3) proposed fix from the finding
D = opo + sp.Rational(1, 2) * a * O1
Ddot = a * H * (sp.Rational(3, 2) * O1 + sp.Rational(1, 2) * a * O2)
S3 = 3 * O1 + 4 * a * O2 + a2 * O3
T3 = opo + sp.Rational(7, 2) * a * O1 + sp.Rational(3, 2) * a2 * O2
N1p = (-sp.Rational(1, 2) * a * H**2 * (H**2 + 3 * Hd) * S3 - sp.Rational(1, 2) * a2 * H**4 * (7 * O2 + 6 * a * O3 + a2 * O4)
       - (Hd**2 + H * Hdd) * T3 - sp.Rational(1, 2) * a * H**2 * Hd * (9 * O1 + 13 * a * O2 + 3 * a2 * O3)
       - sp.Rational(1, 2) * (Pdd + 4 * H * Pd + 2 * Hd * Pm + 4 * H**2 * Pm)
       - sp.Rational(1, 2) * (Ldd + 4 * H * Ld + 2 * Hd * L + 4 * H**2 * L))
prop3 = to_symbols((N1p - Ddot * Hdd) / D)
diff_prop = sp.simplify(sp.together((true3 - prop3).subs(onshell)))
print("proposed d(Hdotdot)/dtau - Hdotdotdot (on shell)     :", diff_prop)

# 4) attribution: in the Omega=0 limit (quintessence, ADE), the only remaining error is the missing Ldd
zero_om = {s: 0 for s in so}
z = sp.simplify((true3 - ship3).subs(onshell).subs(zero_om))
print("shipped error at Omega=0 (should be 0 if Ldd supplied):", z)

# 5) attribution of each shipped defect: fix them one at a time and see what is left (Omega != 0)
def shipped_variant(fix_typo=False, fix_pref=False, fix_bracket=False, fix_ddot=False):
    t3 = T3 if fix_typo else (opo + sp.Rational(7, 2) * a * O1 * sp.Rational(3, 2) * a2 * O2)
    n1 = -a * H * (3 * O1 + (a * O2 if fix_pref else sp.Rational(1, 2) * a * O2))
    br = ((-sp.Rational(1, 2) * a * H**3 * S3 - H * Hd * t3 - sp.Rational(1, 2) * (Pd + 2 * H * Pm + Ld + 2 * H * L))
          if fix_bracket else
          (-sp.Rational(1, 2) * a * H**3 * S3 - H * Hd * t3
           - sp.Rational(1, 2) * (Pdd + 2 * H * Pd - sp.Rational(1, 2) * (Ld + 2 * H * L))))
    first = (sp.Rational(-1, 2) * a * H**2 * (H**2 + 3 * Hd) * S3) if fix_ddot else (-a * adotoa2 * (adotoa2 + 3 * H) * S3)
    return to_symbols(n1 / 2 / D**2 * br + 1 / D * (
        first - sp.Rational(1, 2) * a2 * H**4 * (7 * O2 + 6 * a * O3 + a2 * O4) - (Hd**2 + H * Hdd) * t3
        - sp.Rational(1, 2) * a * adotoa2 * Hd * (9 * O1 + 13 * a * O2 + 3 * a2 * O3)
        - sp.Rational(1, 2) * (Pdd + 4 * H * Pd + 2 * Hd * Pm + 4 * adotoa2 * Pm)
        - sp.Rational(1, 2) * (Ldd + 4 * H * Ld + 2 * Hd * L + 4 * adotoa2 * L)))


v_all = shipped_variant(True, True, True, True)
print("all four corrections applied, error on shell == 0 ?  :",
      sp.simplify(sp.together((true3 - v_all).subs(onshell))) == 0)
for name, kw in (("typo '*'->'+' in T3", dict(fix_pref=True, fix_bracket=True, fix_ddot=True)),
                 ("prefactor 0.5 a O'' -> a O''", dict(fix_typo=True, fix_bracket=True, fix_ddot=True)),
                 ("bracket (Pdd,-0.5L) -> N", dict(fix_typo=True, fix_pref=True, fix_ddot=True)),
                 ("-a H^2 (H^2+3H) -> -0.5 a H^2 (H^2+3Hd)", dict(fix_typo=True, fix_pref=True, fix_bracket=True))):
    v = shipped_variant(**kw)
    r = sp.simplify(sp.together((true3 - v).subs(onshell)))
    print(f"  leave only [{name}] unfixed -> residual is zero? {r == 0}")
