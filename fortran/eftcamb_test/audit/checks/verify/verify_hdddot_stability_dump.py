"""Verifier: read the DebugEFTCAMB stability cache dump (ADD tree, EFT_mass_stability=T) and check, on the stability
path that feeds the positivity bounds, (1) Hdotdotdot vs a*H*d(Hdotdot)/da (spline in a), (2) EFTLambdadotdot as left by
compute_additional_derivs vs the true a*H*dLd/da - 2 H Ld and vs -2 H Ld. Usage: python ... <dump prefix>"""

import sys

import numpy as np
from scipy.interpolate import CubicSpline

pre = sys.argv[1]
frw = np.loadtxt(pre + "cache_FRW.dat")
bg = np.loadtxt(pre + "cache_BackgroundEFT.dat")
a, H, Hd, Hdd, Hddd = frw[:, 0], frw[:, 3], frw[:, 4], frw[:, 5], frw[:, 6]
assert np.allclose(a, bg[:, 0])
Ld, Ldd_dump = bg[:, 11], bg[:, 12]
m = np.concatenate([[True], np.diff(a) > 0])
a, H, Hd, Hdd, Hddd, Ld, Ldd_dump = (v[m] for v in (a, H, Hd, Hdd, Hddd, Ld, Ldd_dump))
num = a * H * CubicSpline(a, Hdd)(a, 1)
numc = a * H * CubicSpline(a, Hd)(a, 1)
Ldd_true = a * H * CubicSpline(a, Ld)(a, 1) - 2 * H * Ld
print(f"rows {len(a)}, a range {a.min():.3e}..{a.max():.3e}")
print(f"{'a':>8s} {'Hdddot dump':>12s} {'aH dHdd/da':>12s} {'Hdot ctrl rel':>13s} {'Ldd dump':>12s} {'Ldd true':>12s} {'-2HLd':>12s}")
for x in (0.05, 0.1, 0.2, 0.5, 0.9):
    i = np.argmin(abs(a - x))
    print(f"{a[i]:8.4f} {Hddd[i]:12.4e} {num[i]:12.4e} {abs(Hdd[i] - numc[i]) / abs(numc[i]):13.2e} "
          f"{Ldd_dump[i]:12.4e} {Ldd_true[i]:12.4e} {-2 * H[i] * Ld[i]:12.4e}")
sel = (a > 0.03) & (a < 0.95)
print("late-time (0.03<a<0.95) max|Hdddot-num|/max|num| =", f"{np.max(abs(Hddd - num)[sel]) / np.max(abs(num[sel])):.3e}")
print("late-time control max|Hdotdot-aH dHdot/da|/max  =", f"{np.max(abs(Hdd - numc)[sel]) / np.max(abs(numc[sel])):.3e}")
print("Ldd dump vs -2 H Ld  max rel =", f"{np.max(abs(Ldd_dump + 2 * H * Ld)[sel] / np.maximum(abs(Ldd_dump[sel]), 1e-300)):.3e}")
print("Ldd dump vs true     max rel =", f"{np.max(abs(Ldd_dump - Ldd_true)[sel]) / np.max(abs(Ldd_true[sel])):.3e}")
