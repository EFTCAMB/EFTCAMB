"""Independent derivation of the Extended (generalized cubic covariant) Galileon tracker background.

Source: Frusciante, Peirone, Atayde, De Felice, arXiv:1912.07586, eq. (16):

    (H/(a H0))^(2+s) = Omega_phi0 + [Omega_m0/a^3 + Omega_r0/a^4 + Omega_nu0 rho_nu/rho_nu0] (H/(a H0))^s

where H is the conformal Hubble rate (their eq. 3 uses conformal time; H = da/dtau / a).  Write E = H/(a H0) and
Omega_tot(a) for the square bracket.  Then F(E, a) = E^(2+s) - Omega_phi0 - Omega_tot(a) E^s = 0 defines E(a).

We derive, by implicit differentiation (no code transcription):
    H     = a H0 E
    Hdot  = dH/dtau   = a H dH/da
    Hddot = d2H/dtau2 = a H d(Hdot)/da
as functions of (a, H, Omega_tot, Omega_tot', Omega_tot''), and check them against a brute-force numerical solution
(brentq for E(a) + finite differences).  The resulting lambdified functions are saved for comparison with the
code's cached adotoa/Hdot/Hdotdot (see ext_galileon_compare.py).

Run: /Users/marco/Software/miniforge3/bin/python ext_galileon_background.py
"""

import numpy as np
import sympy as sp
from scipy.optimize import brentq

a, H0, s, Ophi = sp.symbols("a H0 s Omega_phi0", positive=True)
Ot, Otp, Otpp = sp.symbols("Omega_tot Omega_totp Omega_totpp", positive=True)
E = sp.symbols("E", positive=True)

# implicit function F(E, a) = 0 with Omega_tot(a) treated through its derivatives
F = E ** (2 + s) - Ophi - Ot * E**s
FE = sp.diff(F, E)
Fa = -E**s * Otp  # dF/da at fixed E (only Omega_tot depends on a)
dEda = sp.simplify(-Fa / FE)  # dE/da

# second derivative: d/da of dEda with E = E(a), Omega_tot = Omega_tot(a)
def dda(expr):
    return sp.diff(expr, E) * dEda + sp.diff(expr, Ot) * Otp + sp.diff(expr, Otp) * Otpp


d2Eda2 = sp.simplify(dda(dEda))

# conformal Hubble and its conformal-time derivatives
Hc = a * H0 * E
dHc_da = H0 * E + a * H0 * dEda
Hdot = sp.simplify(a * Hc * dHc_da)
# d/da of Hdot along the solution
dHdot_da = sp.diff(Hdot, a) + sp.diff(Hdot, E) * dEda + sp.diff(Hdot, Ot) * Otp + sp.diff(Hdot, Otp) * Otpp
Hddot = sp.simplify(a * Hc * dHdot_da)

print("dE/da      =", dEda)
print("Hdot       =", Hdot)

f_Hdot = sp.lambdify((a, H0, s, Ophi, E, Ot, Otp), Hdot, "numpy")
f_Hddot = sp.lambdify((a, H0, s, Ophi, E, Ot, Otp, Otpp), Hddot, "numpy")

# ---------------------------------------------------------------------------------------------
# brute-force check with an explicit Omega_tot(a) = Om/a^3 + Or/a^4 (no massive neutrinos)
Om, Or, OphiV, sV, H0V = 0.2747, 8.5e-5, 0.7252, 2.0, 2.335e-4


def omtot(x):
    return Om / x**3 + Or / x**4


def omtot_p(x):
    return -3 * Om / x**4 - 4 * Or / x**5


def omtot_pp(x):
    return 12 * Om / x**5 + 20 * Or / x**6


def E_of_a(x):
    return brentq(lambda e: e ** (2 + sV) - OphiV - omtot(x) * e**sV, 1e-3, 1e8, xtol=1e-15, rtol=1e-15)


def Hc_of_tau_derivs(x, h=1e-4):
    """Numerical Hdot, Hddot via finite differences in ln a: d/dtau = a H d/da = H d/dln a."""
    lx = np.log(x)
    Hs = np.array([x0 * H0V * E_of_a(x0) for x0 in np.exp(lx + h * np.arange(-2, 3))])
    Hc0 = Hs[2]
    dH = (Hs[0] - 8 * Hs[1] + 8 * Hs[3] - Hs[4]) / (12 * h)  # dH/dlna
    hdot = Hc0 * dH
    # Hdot on the stencil, then derivative again
    hd = []
    for j in range(-2, 3):
        xj = np.exp(lx + h * j)
        Hj = xj * H0V * E_of_a(xj)
        Hsj = np.array([x0 * H0V * E_of_a(x0) for x0 in np.exp(np.log(xj) + h * np.arange(-2, 3))])
        hd.append(Hj * (Hsj[0] - 8 * Hsj[1] + 8 * Hsj[3] - Hsj[4]) / (12 * h))
    hd = np.array(hd)
    hddot = Hc0 * (hd[0] - 8 * hd[1] + 8 * hd[3] - hd[4]) / (12 * h)
    return Hc0, hdot, hddot


worst1 = worst2 = 0.0
for x in [0.25, 0.4, 0.6, 0.8, 1.0]:
    e = E_of_a(x)
    Hc0, hd_num, hdd_num = Hc_of_tau_derivs(x)
    hd_sym = f_Hdot(x, H0V, sV, OphiV, e, omtot(x), omtot_p(x))
    hdd_sym = f_Hddot(x, H0V, sV, OphiV, e, omtot(x), omtot_p(x), omtot_pp(x))
    r1 = abs(hd_sym / hd_num - 1)
    r2 = abs(hdd_sym / hdd_num - 1)
    worst1, worst2 = max(worst1, r1), max(worst2, r2)
    print(f"a={x:4.2f} Hdot sym/num-1={r1:.2e}  Hddot sym/num-1={r2:.2e}")
print(f"CHECK [arXiv:1912.07586 eq.16, implicit differentiation vs brute force]: "
      f"worst Hdot {worst1:.2e}, worst Hddot {worst2:.2e} -> {'PASS' if max(worst1, worst2) < 1e-6 else 'FAIL'}")
