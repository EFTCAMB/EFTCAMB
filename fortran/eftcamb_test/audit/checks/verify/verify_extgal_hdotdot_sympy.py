"""Adversarial verifier (t3-extgal-gr-point-background): symbolic check of the Extended Galileon Hdot / Hdotdot formulas.

Tracker: E^(2+S) = Omega_phi0 + Omega(a) E^S, E = adotoa/(a h0), adotoa = conformal Hubble.
Hdot = d adotoa/d tau = a adotoa d adotoa/da; Hdotdot = a adotoa d Hdot/da (implicit differentiation).
Compare with the expressions transcribed from 008p7_Extended_Galileon.f90 ComputeHubbleDer (lines ~525-541), and with
the same Hdotdot expression with its denominator ((S+2) H^2 A S Om)^3 replaced by (A S Om - (S+2) H^2)^3.
"""

import sympy as sp

a, h0, S, H, Om, Om1, Om2, ophi = sp.symbols("a h0 S H Om Om1 Om2 Omphi0", positive=True)
Of = sp.Function("O")(a)
Hf = sp.Function("Hc")(a)
E = Hf / (a * h0)
F = E ** (2 + S) - ophi - Of * E**S
dH = sp.solve(sp.Eq(sp.diff(F, a), 0), sp.diff(Hf, a))[0]
hdot = a * Hf * dH
hdotdot = a * Hf * sp.diff(hdot, a)
hdotdot = hdotdot.subs(sp.Derivative(Hf, a), dH)
sub = {sp.Derivative(Of, (a, 2)): Om2, sp.Derivative(Of, a): Om1}
hdot_s = sp.simplify(hdot.subs(sub).subs({Of: Om, Hf: H}))
hdotdot_s = hdotdot.subs(sub).subs({Of: Om, Hf: H})
# eliminate ophi using the tracker equation itself (code formulas do not contain Omega_phi0)
ophi_val = sp.solve(sp.Eq(F.subs({Of: Om, Hf: H}), 0), ophi)[0]
hdot_s = sp.simplify(hdot_s.subs(ophi, ophi_val))
hdotdot_s = sp.simplify(hdotdot_s.subs(ophi, ophi_val))

A = a**2 * h0**2
code_hdot = -(H**2 * (A * (a * Om1 - S * Om) + (S + 2) * H**2)) / (A * S * Om - (S + 2) * H**2)
num = (A * (S + 2) ** 2 * H**4 * (a * (5 * Om1 + a * Om2) - 6 * S * Om)
       + a**6 * h0**6 * S * Om * (-a**2 * (S + 2) * Om1**2 - 2 * S**2 * Om**2 + a * S * Om * (5 * Om1 + a * Om2))
       + a**4 * h0**4 * S * (S + 2) * H**2 * (a**2 * Om1**2 + 6 * S * Om**2 - 2 * a * Om * (5 * Om1 + a * Om2))
       + 2 * (S + 2) ** 3 * H**6)
code_hdotdot = H**3 / ((S + 2) * H**2 * a**2 * h0**2 * S * Om) ** 3 * num
fixed_hdotdot = H**3 / (A * S * Om - (S + 2) * H**2) ** 3 * num

print("Hdot code - exact          :", sp.simplify(code_hdot - hdot_s))
print("Hdotdot (fixed denom) - exact:", sp.simplify(fixed_hdotdot - hdotdot_s))
vals = {a: sp.Rational(1, 2), h0: sp.Rational(23349, 10**8), S: 2, Om: sp.Rational(23, 10), Om1: -sp.Rational(138, 10),
        Om2: sp.Rational(1104, 10)}
# a consistent H for these values: solve E^4 = 0.725 + Om E^2 at a=0.5
Ev = sp.nsolve(sp.Symbol("e") ** 4 - sp.Rational(725, 1000) - vals[Om] * sp.Symbol("e") ** 2, sp.Symbol("e"), 2)
vals[H] = Ev * vals[a] * vals[h0]
print("numbers at a=0.5, S=2, Om=2.3: exact Hdotdot =", sp.N(hdotdot_s.subs(vals)), " code =",
      sp.N(code_hdotdot.subs(vals)), " fixed-denominator =", sp.N(fixed_hdotdot.subs(vals)))
print("S->0 limit of exact Hdotdot finite:", sp.simplify(sp.limit(hdotdot_s, S, 0)) != sp.zoo)
print("code Hdotdot at S=0 (denominator):", sp.simplify(((S + 2) * H**2 * a**2 * h0**2 * S * Om).subs(S, 0)))
