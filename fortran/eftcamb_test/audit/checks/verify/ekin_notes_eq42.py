"""Compare EFT_kinetic (06_abstract_EFTCAMB_model.f90:690) with the documented kinetic (no-ghost) condition of the
EFTCAMB numerical notes, arXiv:1405.3590 Sec. F, eqs. (42), (44)-(47), valid for gamma3 = -gamma4, gamma6 = 0:
    W2 [4 W1 W2 - W3^2] > 0
    W1 = c a^2/m0^2 + 2 H0^2 a^2 gamma1 - 3 H^2 (1+Omega) - 3 a H^2 Omega' + 3 H^2 gamma4 - 3 a H0 H gamma2
    W2 = -3 [(1+Omega) - gamma4]
    W3 = 6 H (1+Omega) + 3 a H Omega' - 6 H gamma4 + 3 a H0 gamma2
(H conformal Hubble = adotoa, H0 = h0_Mpc, c a^2/m0^2 = EFTc, Omega' = dOmega/da.)
"""

import sympy as sp

a, H, h0, c, g1, g2, g4, Om, OP = sp.symbols("a H h0 c gamma1 gamma2 gamma4 Omega OmegaP", real=True)
Mt = 1 + Om - g4
W1 = c + 2 * h0**2 * a**2 * g1 - 3 * H**2 * (1 + Om) - 3 * a * H**2 * OP + 3 * H**2 * g4 - 3 * a * h0 * H * g2
W2 = -3 * ((1 + Om) - g4)
W3 = 6 * H * (1 + Om) + 3 * a * H * OP - 6 * H * g4 + 3 * a * h0 * g2
notes = sp.expand(W2 * (4 * W1 * W2 - W3**2))

code = 9 * Mt * (4 * c * Mt + 3 * H**2 * OP**2 * a**2 + a**2 * h0 * (h0 * (3 * g2**2 + 8 * g1 * Mt + 6 * H * g2 * OP)))
fixed = 9 * Mt * (4 * c * Mt + 3 * H**2 * OP**2 * a**2 + a**2 * h0 * (h0 * (3 * g2**2 + 8 * g1 * Mt) + 6 * H * g2 * OP))

print("notes eq.(42) W2[4W1W2-W3^2] - fixed =", sp.simplify(notes - fixed))
print("notes eq.(42) W2[4W1W2-W3^2] - code  =", sp.factor(sp.simplify(notes - code)))
print("notes eq.(42) factorised            =", sp.factor(notes))
