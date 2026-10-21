"""Verifier (t8a polterddot-nonflat-gr): analytic check with CAMB's own symbolic hierarchy (camb/symbolic.py).

Derives d^2 polter / dt^2 from camb.symbolic.J_eq(2) and E_eq(2) (CAMB's reference Boltzmann equations, which carry
Kf[L] on the l+1 coupling) and compares it with
  (a) EFTCAMB's replacement in fortran/equations.f90 (pigddot / E2ddot, no Kf(2)),
  (b) CAMB's original polterddot expression (commented in A0, live in B0), after eliminating the GR constraints
      is not needed: we compare the Kf(2)-dependent parts directly, keeping qgdot, octgdot, Edot(3) as symbols,
  (c) the Fortran evolved hierarchies at l=2 (pigdot uses denlk2(2) = 3/5 k Kf(2); Edot(2) uses k/3 without Kf(2)).
Pure sympy, no CAMB run.
"""

import sys

import sympy as sp

sys.path.insert(
    0,
    "/private/tmp/claude-501/-Users-marco-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/"
    "scratchpad/audit/trees/B0",
)
import camb.symbolic as cs  # noqa: E402

t, k = cs.t, cs.k
Kf2 = cs.Kf[2]

# symbols standing for the code variables
qg, qgdot, octg, octgdot, E2, E2dot, E3, E3dot, pig, pigdot = sp.symbols(
    "qg qgdot octg octgdot E2 E2dot E3 E3dot pig pigdot"
)
op, dop, sig, sigdot, polter, polterdot = sp.symbols("opacity dopacity sigma sigmadot polter polterdot")

Jeq = cs.J_eq(2)
Eeq = cs.E_eq(2)
print("camb.symbolic J_eq(2):", Jeq)
print("camb.symbolic E_eq(2):", Eeq)

pi_g, q_g, J3 = cs.pi_g, cs.q_g, cs._make_index_func("J", 3)
E_2, E_3 = cs._make_index_func("E", 2), cs._make_index_func("E", 3)
opacity, sigma, polterf = cs.opacity, cs.sigma, cs.polter

# second derivatives of pi_g and E_2 by differentiating the reference first-order equations
pigdd = sp.diff(Jeq.rhs, t)
E2dd = sp.diff(Eeq.rhs, t)
rep = {
    sp.Derivative(q_g, t): qgdot,
    sp.Derivative(J3, t): octgdot,
    sp.Derivative(pi_g, t): pigdot,
    sp.Derivative(E_2, t): E2dot,
    sp.Derivative(E_3, t): E3dot,
    sp.Derivative(opacity, t): dop,
    sp.Derivative(sigma, t): sigdot,
    sp.Derivative(polterf, t): polterdot,
}
rep2 = {q_g: qg, J3: octg, pi_g: pig, E_2: E2, E_3: E3, opacity: op, sigma: sig, polterf: polter}


def to_code(e):
    return sp.expand(e.subs(rep).subs(rep2))


ref_pigdd = to_code(pigdd)
ref_E2dd = to_code(E2dd)
ref_polterdd = sp.expand(ref_pigdd / 10 + sp.Rational(3, 5) * ref_E2dd)

# (a) EFTCAMB replacement (equations.f90:3858-3861)
eft_pigdd = (sp.Rational(2, 5) * k * qgdot - sp.Rational(3, 5) * k * octgdot - dop * (pig - polter)
             - op * (pigdot - polterdot) + sp.Rational(8, 15) * k * sigdot)
eft_E2dd = -dop * (E2 - polter) - op * (E2dot - polterdot) - k / 3 * E3dot
eft_polterdd = sp.expand(eft_pigdd / 10 + sp.Rational(3, 5) * eft_E2dd)

# (b) Kf(2)-carrying terms of CAMB's original polterddot (B0 equations.f90:2883-2888):
#     - 3/50 k octgdot Kf(2) - 1/5 k Kf(2) Edot(3); the remaining terms are Kf(2)-independent GR-constraint forms
camb_kf_part = -sp.Rational(3, 50) * k * octgdot * Kf2 - sp.Rational(1, 5) * k * Kf2 * E3dot

print("\nreference (camb.symbolic) pigddot  =", ref_pigdd)
print("reference (camb.symbolic) E2ddot   =", ref_E2dd)
d_pig = sp.simplify(eft_pigdd - ref_pigdd)
d_E2 = sp.simplify(eft_E2dd - ref_E2dd)
d_pol = sp.factor(sp.simplify(eft_polterdd - ref_polterdd))
print("\nEFTCAMB pigddot - reference  =", d_pig)
print("EFTCAMB E2ddot  - reference  =", d_E2)
print("EFTCAMB polterddot - ref     =", d_pol)
print("  ... at Kf(2)=1 (flat)      =", sp.simplify(d_pol.subs(Kf2, 1)))
ref_kf_part = sp.expand(ref_polterdd).coeff(Kf2) * Kf2
print("Kf(2) terms: reference =", ref_kf_part, "  CAMB original =", camb_kf_part,
      "  equal:", sp.simplify(ref_kf_part - camb_kf_part) == 0)

# (c) Fortran evolved hierarchy at l=2 (equations.f90:3399 and 3417)
denl2 = sp.Rational(1, 5)
fort_pigdot = 2 * k * denl2 * qg - 3 * k * Kf2 * denl2 * octg - op * (pig - polter) + sp.Rational(8, 15) * k * sig
fort_E2dot = -op * (E2 - polter) - k / 3 * E3
ref_pigdot = to_code(Jeq.rhs)
ref_E2dot = to_code(Eeq.rhs)
print("\nFortran pigdot (3399) - camb.symbolic =", sp.simplify(fort_pigdot - ref_pigdot))
print("Fortran E2dot  (3417) - camb.symbolic =", sp.simplify(fort_E2dot - ref_E2dot),
      "   <- upstream CAMB hierarchy omits Kf(2) at l=2 for E")
# time derivative of the Fortran evolved pigdot / E2dot (what polterddot should be for self-consistency)
fd_pig = sp.expand(fort_pigdot.subs({qg: qg, octg: octg}))
dpig_dt = (sp.Rational(2, 5) * k * qgdot - sp.Rational(3, 5) * k * Kf2 * octgdot - dop * (pig - polter)
           - op * (pigdot - polterdot) + sp.Rational(8, 15) * k * sigdot)
dE2_dt = -dop * (E2 - polter) - op * (E2dot - polterdot) - k / 3 * E3dot
print("EFTCAMB pigddot - d/dt(Fortran pigdot) =", sp.simplify(eft_pigdd - dpig_dt))
print("EFTCAMB E2ddot  - d/dt(Fortran E2dot)  =", sp.simplify(eft_E2dd - dE2_dt))
