"""Independent check for t3-fr-rgr-switch-gr-approach: quasi-static linear growth for designer f(R).

Uses only the background EFT functions written by EFTCAMB (Omega = f_R, its a-derivative, conformal H, H', H'',
rho_b, rho_c), which do not depend on the return-to-GR time, and solves
    D'' + (2 + dlnH/dlna) D' = (1/2) (grhob_t + grhoc_t)/Hconf^2 * mu(k, a) * D      (' = d/dlna)
with the standard f(R) quasi-static mu = (1+4Q)/((1+3Q)(1+f_R)), Q = (k/a)^2 f_RR/(1+f_R), f_RR = (dOmega/da)/(dR/da),
R = 6 (H' + H^2)/a^2 (conformal).  Starts at a_start in matter domination with D = a, D' = a.
Prints D_fR(k)^2/D_GR^2 - 1 at a = 1 for k = 0.3, 1, 3, 7.3 /Mpc, next to the CAMB P(k) ratios given on the command line.

Usage: python frrgr_qs_growth.py npz [a_start]
"""

import sys
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicSpline

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import harness  # noqa: E402

KS = (0.31, 1.0, 3.1, 7.3)


def main() -> None:
    d = harness.load(sys.argv[1])
    a_start = float(sys.argv[2]) if len(sys.argv) > 2 else 0.01
    a = d["eft_a_grid"]
    sel = (a >= a_start / 2) & (a <= 1.0)
    la = np.log(a[sel])
    H = d["eft_adotoa"][sel]
    Hd = d["eft_Hdot"][sel]
    Hdd = d["eft_Hdotdot"][sel]
    om = d["eft_EFTOmegaV"][sel]
    omp = d["eft_EFTOmegaP"][sel]
    gcb = d["eft_grhob_t"][sel] + d["eft_grhoc_t"][sel]
    aa = a[sel]
    R = 6 * (Hd + H**2) / aa**2
    dRdtau = 6 * (Hdd + 2 * H * Hd) / aa**2 - 2 * H * R
    dRda = dRdtau / (aa * H)
    # numerical cross-check of dR/da
    dRda_num = np.gradient(R, aa)
    fRR = omp / dRda
    print(f"dR/da analytic vs numeric, max rel diff over a in [{a_start},1]: "
          f"{np.max(np.abs(dRda_num[2:-2] / dRda[2:-2] - 1)):.2e}")
    print(f"f_R(a=1) = {om[-1]:.3e}, f_RR(a=1) = {fRR[-1]:.3e} Mpc^2, "
          f"Compton k_C(a=1) = sqrt((1+f_R)/(3 f_RR)) = {np.sqrt((1 + om[-1]) / (3 * fRR[-1])):.3f} /Mpc")
    sp = {n: CubicSpline(la, v) for n, v in
          (("H", H), ("Hd", Hd), ("om", om), ("fRR", fRR), ("gcb", gcb))}

    def rhs(lna, y, k, mg):
        h, hd = sp["H"](lna), sp["Hd"](lna)
        dlnH = hd / h**2 - 1.0
        src = 0.5 * sp["gcb"](lna) / h**2
        if mg:
            ex = np.exp(lna)
            fr = sp["om"](lna)
            q = (k / ex) ** 2 * sp["fRR"](lna) / (1 + fr)
            mu = (1 + 4 * q) / ((1 + 3 * q) * (1 + fr))
        else:
            mu = 1.0
        return [y[1], -(2 + dlnH) * y[1] + src * mu * y[0]]

    def grow(k, mg):
        s = solve_ivp(rhs, (np.log(a_start), 0.0), [a_start, a_start], args=(k, mg), rtol=1e-10, atol=1e-14,
                      method="DOP853")
        return s.y[0, -1]

    dgr = grow(1.0, False)
    out = []
    for k in KS:
        out.append((grow(k, True) / dgr) ** 2 - 1)
    print("QS growth D_fR^2/D_GR^2 - 1 at a=1: " + " ".join(f"k={k}: {v:+.3e}" for k, v in zip(KS, out)))


if __name__ == "__main__":
    main()
