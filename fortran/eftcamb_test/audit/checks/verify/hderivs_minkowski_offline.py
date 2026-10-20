"""Offline re-evaluation of the Minkowski-limit condition (06_abstract_EFTCAMB_model.f90:1104-1173,
09_EFTCAMB_stability.f90:410-417) from the reference cache fields in runs/A0/<model>.npz.

For each model we rebuild ghostpos = D1 and tachpos = -(Q1+...+Q9)/D1 with the code's formulas, using the stored
fields, once with EFTLambdadotdot / Hdotdotdot as stored by the model (what the stability routine sees: the cache is
re-initialised to 0 by EFTTestStability and nothing assigns them before compute_positivity_bounds), and once with
EFTLambdadotdot = spline derivative of EFTLambdadot (EFTCAMB convention) and Hdotdotdot = spline derivative of Hdotdot.
Gamma2PP/PPP and Gamma3PP/PPP/PPPP are taken as stored (they are 0 for the models used here, where gamma2 = gamma3 = 0;
this is checked and printed).  Counts are over 1e-3 <= a <= 0.95 (401-point log grid).

Usage: python hderivs_minkowski_offline.py
"""

from pathlib import Path

import numpy as np
from scipy.interpolate import CubicSpline

RUNS = Path("/private/tmp/claude-501/-Users-marco-Desktop-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/"
            "scratchpad/audit/runs/A0")
H0_MPC = None


def sd(lna, y):
    return CubicSpline(lna, y)(lna, 1)


def bounds(g, a, h0, Ldd, Hddd):
    H, Hd, Hdd = g["adotoa"], g["Hdot"], g["Hdotdot"]
    a2 = a * a
    G2, G2P, G2PP, G2PPP = g["EFTGamma2V"], g["EFTGamma2P"], g["EFTGamma2PP"], g["EFTGamma2PPP"]
    G3, G3P, G3PP, G3PPP, G3PPPP = (g["EFTGamma3V"], g["EFTGamma3P"], g["EFTGamma3PP"], g["EFTGamma3PPP"],
                                    g["EFTGamma3PPPP"])
    D1 = -2 * g["EFTc"] + 3 * a * H * h0 * G2 + (H**2 + 2 * Hd) * G3 + 2 * a * H**2 * G3P
    Q1 = (-H * g["EFTLambdadot"] + Ldd) * a2
    Q2 = h0 * (Hdd - H * Hd) * a2 * G2P / 2
    Q3 = 3 * h0 * H * Hd * a**3 * G2PP / 2
    Q4 = h0 * H**3 * a2**2 * G2PPP / 2
    Q5 = (-Hddd - 12 * H**4 + 18 * H**2 * Hd - 2 * Hd**2 + H * Hdd) * G3
    Q6 = (a / 2) * (-Hddd + 12 * H**4 - 5 * Hd**2 - 5 * H * Hdd) * G3P
    Q7 = -(a2 / 2) * (3 * Hd**2 + 8 * H**2 * Hd + 4 * Hdd * H) * G3PP
    Q8 = -(a**3) * H**2 * (3 * Hd + H**2) * G3PPP
    Q9 = -a2**2 * H**4 * G3PPPP / 2
    tach = -(Q1 + Q2 + Q3 + Q4 + Q5 + Q6 + Q7 + Q8 + Q9) / D1
    fail = (-D1 <= 0) | (-tach < 0)
    return D1, tach, fail


models = ["5_quint_exponential", "5_quint_inverse_power", "5_quint_cubic", "5_quint_axion", "5_quint_quadratic",
          "5_ADE_1", "5_Horava_eta_1", "5_Horava_eta_3", "5_Horava_eta_4", "5_Horava_SolSyst_2", "5_Horava_SolSyst_3"]
print(f"{'model':24s} {'max|g2|+|g3|':>12s} {'fail(code fields)':>18s} {'fail(fd Ldd,Hddd)':>18s}"
      f"  a=0.5: Ldd stored / fd ; tachpos code / fd")
for m in models:
    d = np.load(RUNS / f"{m}.npz", allow_pickle=True)
    a = d["eft_a_grid"]
    lna = np.log(a)
    g = {k[4:]: d[k] for k in d.files if k.startswith("eft_")}
    # h0_Mpc from H at a=1: adotoa(a=1) = H0 in 1/Mpc
    h0 = float(g["adotoa"][-1])
    H = g["adotoa"]
    Ldd_fd = H * sd(lna, g["EFTLambdadot"]) - 2 * H * g["EFTLambdadot"]
    Hddd_fd = H * sd(lna, g["Hdotdot"])
    sel = (a >= 1e-3) & (a <= 0.95)
    gam = float(np.max(np.abs(g["EFTGamma2V"])) + np.max(np.abs(g["EFTGamma3V"])))
    _, t0, f0 = bounds(g, a, h0, g["EFTLambdadotdot"], g["Hdotdotdot"])
    _, t1, f1 = bounds(g, a, h0, Ldd_fd, Hddd_fd)
    i5 = int(np.argmin(abs(a - 0.5)))
    print(f"{m:24s} {gam:12.1e} {int(f0[sel].sum()):>8d}/{int(sel.sum()):<9d} {int(f1[sel].sum()):>8d}/{int(sel.sum()):<9d}"
          f"  {g['EFTLambdadotdot'][i5]:.3e} / {Ldd_fd[i5]:.3e} ; {t0[i5]:.4e} / {t1[i5]:.4e}")
