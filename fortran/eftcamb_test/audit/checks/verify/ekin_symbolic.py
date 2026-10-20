"""Independent symbolic check of EFT_kinetic (06_abstract_EFTCAMB_model.f90:690).

1. Dimensional homogeneity: assign length dimension L^-1 to adotoa (conformal H) and h0_Mpc, L^-2 to EFTc,
   dimensionless to a, Omega, Omega' (d/da), gamma_i.  Every term of a physical no-ghost coefficient must have the same
   dimension.
2. RPH mapping (007p2_RPH.f90): substitute the code's own EFT functions in terms of alpha_M/B/K/T and compare to
   18 M^3 H^2 (alpha_K + 6 alpha_B^2).
3. GLPV (Gleyzes-Langlois-Piazza-Vernizzi 2013) no-ghost condition (c + 2 M2^4)(M^2 f + 2 m4^2) + 3/4 (m3^3 -+ M^2 fdot)^2
   in EFTCAMB variables, with the relative sign of m3^3 vs M^2 fdot fixed by EFTCAMB's own pi equation
   (EFTpiA1 contains (adotoa Omega' + h0 gamma2)^2).
4. Compare with 36 Mt^2 EFTpiA1 (06_abstract_EFTCAMB_model.f90:515), with gamma3 = -gamma4 (the only case where the
   ghost check is applied, 09_EFTCAMB_stability.f90:301).
"""

import sympy as sp

a, H, h0, c, g1, g2, g3, g4, Om, OP = sp.symbols("a H h0 c gamma1 gamma2 gamma3 gamma4 Omega OmegaP", real=True)
lam = sp.symbols("lambda", positive=True)
Mt = 1 + Om - g4

# code as written (line 690-694)
code = 9 * Mt * (4 * c * Mt + 3 * H**2 * OP**2 * a**2 + a**2 * h0 * (h0 * (3 * g2**2 + 8 * g1 * Mt + 6 * H * g2 * OP)))
# proposed fix (single h0 on the cross term)
fixed = 9 * Mt * (4 * c * Mt + 3 * H**2 * OP**2 * a**2 + a**2 * h0 * (h0 * (3 * g2**2 + 8 * g1 * Mt) + 6 * H * g2 * OP))

print("1) dimensional homogeneity (scale H,h0 -> lam*H, lam*h0; c -> lam^2 c): each term must scale as lam^2")
for name, expr in (("code", code), ("fixed", fixed)):
    inner = sp.expand(expr / (9 * Mt))
    powers = set()
    for term in sp.Add.make_args(inner):
        scaled = term.subs({H: lam * H, h0: lam * h0, c: lam**2 * c}, simultaneous=True)
        powers.add(sp.degree(sp.simplify(scaled / term), lam))
    print(f"   {name}: powers of lambda over terms = {sorted(powers)}")

print("2) RPH mapping (alpha_T = 0 for brevity, then general alpha_T)")
M, aM, aB, aK, aT = sp.symbols("M alpha_M alpha_B alpha_K alpha_T", real=True)
# 007p2_RPH.f90: Omega = M-1 + aT*M ; a*Omega' = aM*M(1+aT) + a aT' M (take aT' = 0 -> a Omega' = aM M (1+aT))
# gamma2 = (2 aB M - a Omega') H / (h0 a);  gamma1 = (aK M H^2 - 2 c)/(4 h0^2 a^2);  gamma3 = -aT M; gamma4 = -gamma3
subs = {
    Om: M - 1 + aT * M,
    OP: aM * M * (1 + aT) / a,
    g2: (2 * aB * M - aM * M * (1 + aT)) * H / (h0 * a),
    g1: (aK * M * H**2 - 2 * c) / (4 * h0**2 * a**2),
    g4: aT * M,
}
target = 18 * M**3 * H**2 * (aK + 6 * aB**2)
fx = sp.simplify(fixed.subs(subs, simultaneous=True) - target)
cd = sp.factor(sp.simplify(code.subs(subs, simultaneous=True) - target))
print(f"   fixed - 18M^3H^2(aK+6aB^2) = {fx}")
print(f"   code  - 18M^3H^2(aK+6aB^2) = {cd}")
cd0 = sp.factor(sp.simplify(code.subs(subs, simultaneous=True).subs({aB: 0, aT: 0}) / (18 * M**3 * H**2)))
print(f"   code/(18 M^3 H^2) at alpha_B = alpha_T = 0: {cd0}   (limit h0 -> 0: {sp.limit(cd0, h0, 0)})")
print("   model B of the finding: aK=-0.2 ODE^2, aM=0.5 ODE  -> code ratio = 1 + 3*0.25/(-0.2) =", 1 + 3 * 0.25 / -0.2)

print("3) GLPV no-ghost: A = (c+2M2^4)(M^2 f + 2 m4^2) + 3/4 (m3^3 + M^2 fdot)^2 (EFTCAMB sign, see pi eq.)")
# conformal units: a^2 c/m0^2 = EFTc; a^2 M2^4/m0^2 = gamma1 a^2 h0^2; a m3^3/m0^2 = a gamma2 h0;
# a M^2 fdot/m0^2 = a Omega' adotoa (dOmega/dt = Omega' da/dt = Omega' * conformal H); M^2 f + 2 m4^2 = m0^2 Mt
glpv = (c + 2 * g1 * a**2 * h0**2) * Mt + sp.Rational(3, 4) * (a * g2 * h0 + a * OP * H) ** 2
print(f"   36 Mt * GLPV - fixed = {sp.simplify(36 * Mt * glpv - fixed)}")
print(f"   36 Mt * GLPV - code  = {sp.factor(sp.simplify(36 * Mt * glpv - code))}")

print("4) internal consistency with the pi-equation kinetic coefficient EFTpiA1 (line 515), gamma3 = -gamma4")
A1 = c + 2 * a**2 * h0**2 * g1 + sp.Rational(3, 2) * a**2 * (H * OP + h0 * g2) ** 2 / (2 * (1 + Om) + 3 * g3 + g4)
A1h = A1.subs(g3, -g4)
print(f"   36 Mt^2 A1 - fixed = {sp.simplify(36 * Mt**2 * A1h - fixed)}")
print(f"   36 Mt^2 A1 - code  = {sp.factor(sp.simplify(36 * Mt**2 * A1h - code))}")
