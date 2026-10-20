"""Independent check for fr-rgr-threshold-high-k-pk: quasi-static (QS) linear growth of designer f(R) vs turn-on time.

Uses only the background EFT functions of one EFTCAMB run (Omega = f_R, dOmega/da, conformal H, H', H'', rho_b, rho_c;
these do not depend on the perturbation turn-on) and solves, for each k,
    D'' + (2 + dlnH/dlna) D' = (1/2) (grhob_t + grhoc_t)/Hconf^2 * mu(k, a) * D     (' = d/dlna)
with mu = 1 for a < a_on and the f(R) QS mu = (1+4Q)/((1+3Q)(1+f_R)), Q = (k/a)^2 f_RR/(1+f_R), f_RR = Omega_a/R_a,
R = 6 (H' + H^2)/a^2, for a >= a_on.  Starts at a_ini = 0.002 with D = a, D' = a (same for GR, so the ratio is
insensitive to the start).  Prints D_fR^2/D_GR^2 - 1 at a = 1 for several k and a_on, plus mu - 1 at a few epochs.

Usage: python frthr_qs.py RUN.npz a_on[,a_on...]
"""

import sys
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicSpline

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import harness  # noqa: E402

KS = (0.1, 0.3, 1.0, 3.0, 10.0)
A_INI = 0.002


def main() -> None:
    d = harness.load(sys.argv[1])
    a_ons = [float(x) for x in sys.argv[2].split(",")]
    a = d["eft_a_grid"]
    sel = (a >= A_INI / 2) & (a <= 1.0)
    aa, la = a[sel], np.log(a[sel])
    H, Hd, Hdd = d["eft_adotoa"][sel], d["eft_Hdot"][sel], d["eft_Hdotdot"][sel]
    om, omp = d["eft_EFTOmegaV"][sel], d["eft_EFTOmegaP"][sel]
    gcb = d["eft_grhob_t"][sel] + d["eft_grhoc_t"][sel]
    R = 6 * (Hd + H**2) / aa**2
    dRda = (6 * (Hdd + 2 * H * Hd) / aa**2 - 2 * H * R) / (aa * H)
    fRR = omp / dRda
    chk = np.gradient(R, aa)
    print(f"dR/da analytic vs numeric max rel diff (interior): {np.max(np.abs(chk[3:-3] / dRda[3:-3] - 1)):.1e}")
    print(f"f_R(a=1)={om[-1]:+.3e}  f_RR(a=1)={fRR[-1]:+.3e} Mpc^2  k_C(a=1)={np.sqrt((1 + om[-1]) / (3 * fRR[-1])):.3f}/Mpc")
    sp = {n: CubicSpline(la, v) for n, v in (("H", H), ("Hd", Hd), ("om", om), ("fRR", fRR), ("gcb", gcb))}

    def mu(lna, k):
        fr = sp["om"](lna)
        q = (k / np.exp(lna)) ** 2 * sp["fRR"](lna) / (1 + fr)
        return (1 + 4 * q) / ((1 + 3 * q) * (1 + fr))

    for ae in (0.01, 0.0239, 0.0471, 0.0931, 0.1833, 0.5, 1.0):
        print(f"  a={ae:<7g} f_R={sp['om'](np.log(ae)):+.2e}  mu-1: "
              + " ".join(f"k={k}:{mu(np.log(ae), k) - 1:+.2e}" for k in KS))

    def rhs(lna, y, k, lon):
        h, hd = sp["H"](lna), sp["Hd"](lna)
        m = mu(lna, k) if lna >= lon else 1.0
        return [y[1], -(hd / h**2 + 1.0) * y[1] + 0.5 * sp["gcb"](lna) / h**2 * m * y[0]]

    def grow(k, a_on):
        lon = np.log(a_on) if a_on is not None else 1.0
        s = solve_ivp(rhs, (np.log(A_INI), 0.0), [A_INI, A_INI], args=(k, lon), rtol=1e-10, atol=1e-14,
                      method="DOP853", max_step=0.01)
        return s.y[0, -1]

    dgr = grow(1.0, None)
    print("QS D_fR^2/D_GR^2 - 1 at a=1 (mu=1 before a_on):")
    for a_on in a_ons:
        print(f"  a_on={a_on:<8g} " + " ".join(f"k={k}:{(grow(k, a_on) / dgr) ** 2 - 1:+.4e}" for k in KS))


if __name__ == "__main__":
    main()
