"""Repro (track 10a): the generic full-mapping third Hubble derivative EFTCAMBFullModelComputeHubbleDer
(fortran/eftcamb/06p1_abstract_EFTCAMB_full_map.f90, eft_cache%Hdotdotdot) is wrong whenever Omega(a) is not constant.

Used by every full-mapping model that does not override compute_H_derivs (K-mouflage 008p3, quintessence 008p4,
ADE 008p2, Beyond Horndeski 008p5).  Hdotdotdot enters the tachyon/mass positivity check (06_abstract_EFTCAMB_model.f90
Q5, Q6 -> tach_pos_num) and the cache output.  The same code is in the pre-rebase EFTCAMB (not a regression).

Check: run 5_Kmouflage_1 (Omega ~ -0.09) in the frozen A0 build, take the EFT functions on a dense a-grid and compare
  Hdotdot     with d(Hdot)/dtau    = a*H * dHdot/da       (control: the second derivative is right)
  Hdotdotdot  with d(Hdotdot)/dtau = a*H * dHdotdot/da
  Hdotdotdot  with the analytic derivative of the Hdotdot formula (corrected expression, evaluated here in Python)
Expected: all three agree to ~1e-4.  Observed: the Fortran Hdotdotdot is off by O(1) or more; the corrected
expression matches the numerical derivative.

Run: /Users/marco/Software/miniforge3/bin/python t10a-fullmap-hdotdotdot-formula.py [TREE]
"""

import os
import sys
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")
TREE = sys.argv[1] if len(sys.argv) > 1 else (
    "/private/tmp/claude-501/-Users-marco-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/"
    "scratchpad/audit/trees/A0")
sys.path.insert(0, TREE)
import camb  # noqa: E402
import numpy as np  # noqa: E402
from scipy.interpolate import CubicSpline  # noqa: E402

INI = Path(__file__).resolve().parents[2] / "parameters" / "5_Kmouflage_1.ini"

cp = camb.read_ini(str(INI), no_validate=True)
cp.WantTransfer = False  # (WantCls = False segfaults in this setup, see t10a-upstream-get-background-windows-segfault)
res = camb.get_results(cp)
a = np.logspace(-3, 0, 3001)
_, f = cp.EFTCAMB.get_eft_functions(res, a)

H, Hd, Hdd, Hddd = f["adotoa"], f["Hdot"], f["Hdotdot"], f["Hdotdotdot"]
O, O1, O2, O3, O4 = (f[k] for k in ("EFTOmegaV", "EFTOmegaP", "EFTOmegaPP", "EFTOmegaPPP", "EFTOmegaPPPP"))
L, Ld, Ldd = f["EFTLambda"], f["EFTLambdadot"], f["EFTLambdadotdot"]
P, Pd, Pdd = f["gpresm_t"], f["gpresdotm_t"], f["gpresdotdotm_t"]

la = np.log(a)
Hdd_num = H * CubicSpline(la, Hd)(la, 1)
Hddd_num = H * CubicSpline(la, Hdd)(la, 1)

# corrected analytic derivative of Hdotdot = N2/D (same conventions as the Fortran code)
X = 3 * O1 + 4 * a * O2 + a**2 * O3
Xp = 7 * O2 + 6 * a * O3 + a**2 * O4
Y = 1 + O + 3.5 * a * O1 + 1.5 * a**2 * O2
Yp = 4.5 * O1 + 6.5 * a * O2 + 1.5 * a**2 * O3
D = 1 + O + 0.5 * a * O1
Ddot = a * H * (1.5 * O1 + 0.5 * a * O2)
N2 = -0.5 * a * H**3 * X - H * Hd * Y - 0.5 * (Pd + 2 * H * P + Ld + 2 * H * L)
N2dot = (-0.5 * a * H**4 * X - 1.5 * a * H**2 * Hd * X - 0.5 * a**2 * H**4 * Xp
         - (Hd**2 + H * Hdd) * Y - a * H**2 * Hd * Yp
         - 0.5 * (Pdd + 4 * H * Pd + 2 * Hd * P + 4 * H**2 * P)
         - 0.5 * (Ldd + 4 * H * Ld + 2 * Hd * L + 4 * H**2 * L))
Hddd_fix = N2dot / D - N2 * Ddot / D**2

sel = (a > 2e-3) & (a < 0.98)
scale = np.max(np.abs(Hddd_num[sel]))


def rel(x, y, s=None):
    s = np.max(np.abs(y[sel])) if s is None else s
    return np.max(np.abs(x[sel] - y[sel])) / s


print(f"5_Kmouflage_1 (A0), a in [2e-3, 0.98], Omega(a) in [{O[sel].min():.3f}, {O[sel].max():.3f}]")
print(f"  control  max|Hdotdot - dHdot/dtau|        / max|.| = {rel(Hdd, Hdd_num):.2e}")
print(f"  Fortran  max|Hdotdotdot - dHdotdot/dtau|  / max|.| = {rel(Hddd, Hddd_num, scale):.2e}")
print(f"  fixed    max|Hdddot_fixed - dHdotdot/dtau| / max|.| = {rel(Hddd_fix, Hddd_num, scale):.2e}")
for ai in (0.01, 0.1, 0.5, 0.9):
    i = int(np.argmin(np.abs(a - ai)))
    print(f"  a={a[i]:.3f}: Hdotdotdot Fortran {Hddd[i]: .4e}  numerical {Hddd_num[i]: .4e}  fixed {Hddd_fix[i]: .4e}")
