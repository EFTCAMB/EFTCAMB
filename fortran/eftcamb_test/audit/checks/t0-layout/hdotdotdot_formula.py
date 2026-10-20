"""Track 0 side check: validate the full-mapping Hdotdotdot formula (06p1_abstract_EFTCAMB_full_map.f90) in python.

Using the fields returned by get_eft_functions (plus finite-difference values for EFTOmegaPPPP and EFTLambdadotdot,
which are not computed on this path for K-mouflage), evaluate
  (a) the shipped formula, transcribed term by term,
  (b) the derivative of the shipped Hdotdot formula, Hdotdotdot = (dN/dtau - Hdotdot ddenom/dtau)/denom,
and compare with the numerical d(Hdotdot)/dtau.

Usage: python hdotdotdot_formula.py [--tree PATH] [--model 5_Kmouflage_1]
"""

import argparse
import os
import sys
import tempfile
from pathlib import Path

import numpy as np

ROOT = "/private/tmp/claude-501/-Users-marco-Desktop-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad/audit"
PARAMS = Path(__file__).resolve().parents[3] / "parameters"
ap = argparse.ArgumentParser()
ap.add_argument("--tree", default=f"{ROOT}/trees/A0")
ap.add_argument("--model", default="5_Kmouflage_1")
args = ap.parse_args()
os.environ.setdefault("OMP_NUM_THREADS", "2")
sys.path.insert(0, args.tree)
import camb  # noqa: E402

tmp = tempfile.mkdtemp(prefix="h3f_")
ini = os.path.join(tmp, "m.ini")
Path(ini).write_text(f"DEFAULT({PARAMS / (args.model + '.ini')})\noutput_root = {tmp}/m\nEFTCAMB_write_background = F\n"
                     "feedback_level = 0\n")
cp = camb.read_ini(ini, no_validate=True)
res = camb.get_background(cp)
a = np.array([0.05, 0.1, 0.2, 0.5, 0.9])
h = 1e-3


def F(x):
    return cp.EFTCAMB.get_eft_functions(res, x)[1]


C, P, M, P2, M2 = F(a), F(a * (1 + h)), F(a * (1 - h)), F(a * (1 + h / 2)), F(a * (1 - h / 2))


def dda(f):  # Richardson-extrapolated central difference in a
    d1 = (P[f] - M[f]) / (2 * h * a)
    d2 = (P2[f] - M2[f]) / (h * a)
    return (4 * d2 - d1) / 3


H, Hd, Hdd = C["adotoa"], C["Hdot"], C["Hdotdot"]
O, O1, O2, O3 = C["EFTOmegaV"], C["EFTOmegaP"], C["EFTOmegaPP"], C["EFTOmegaPPP"]
O4 = C["EFTOmegaPPPP"] if np.any(C["EFTOmegaPPPP"]) else dda("EFTOmegaPPP")
Pm, Pd, Pdd = C["gpresm_t"], C["gpresdotm_t"], C["gpresdotdotm_t"]
L, Ld = C["EFTLambda"], C["EFTLambdadot"]
Ldd = a * H * dda("EFTLambdadot") - 2 * H * Ld  # same definition as compute_additional_derivs
a2, H2 = a * a, H * H
opo = 1 + O
denom = opo + 0.5 * a * O1
S3 = 3 * O1 + 4 * a * O2 + a2 * O3
# (a) shipped formula, term by term (fortran/eftcamb/06p1_abstract_EFTCAMB_full_map.f90:141-163)
N1s = -a * H * (3 * O1 + 0.5 * a * O2)
shipped = (N1s / 2 / denom**2 * (-0.5 * a * H**3 * S3 - H * Hd * (opo + 3.5 * a * O1 * 1.5 * a2 * O2)
                                  - 0.5 * (Pdd + 2 * H * Pd - 0.5 * (Ld + 2 * H * L)))
           + 1 / denom * (-a * H2 * (H2 + 3 * H) * S3 - 0.5 * a2 * H**4 * (7 * O2 + 6 * a * O3 + a2 * O4)
                          - (Hd**2 + H * Hdd) * (opo + 3.5 * a * O1 * 1.5 * a2 * O2)
                          - 0.5 * a * H2 * Hd * (9 * O1 + 13 * a * O2 + 3 * a2 * O3)
                          - 0.5 * (Pdd + 4 * H * Pd + 2 * Hd * Pm + 4 * H2 * Pm)
                          - 0.5 * (Ldd + 4 * H * Ld + 2 * Hd * L + 4 * H2 * L)))
# (b) derivative of Hdotdot = N/denom
denomdot = a * H * (1.5 * O1 + 0.5 * a * O2)
T3 = opo + 3.5 * a * O1 + 1.5 * a2 * O2
Ndot = (-0.5 * a * H2 * (H2 + 3 * Hd) * S3 - 0.5 * a2 * H**4 * (7 * O2 + 6 * a * O3 + a2 * O4)
        - (Hd**2 + H * Hdd) * T3 - 0.5 * a * H2 * Hd * (9 * O1 + 13 * a * O2 + 3 * a2 * O3)
        - 0.5 * (Pdd + 4 * H * Pd + 2 * Hd * Pm + 4 * H2 * Pm) - 0.5 * (Ldd + 4 * H * Ld + 2 * Hd * L + 4 * H2 * L))
fixed = (Ndot - denomdot * Hdd) / denom
num = a * H * dda("Hdotdot")
print(f"model {args.model}, tree {args.tree}")
print(f"{'a':>6s} {'d(Hdotdot)/dtau':>16s} {'field':>12s} {'shipped(py)':>12s} {'corrected':>12s}")
for i in range(len(a)):
    print(f"{a[i]:6.2f} {num[i]:16.4e} {C['Hdotdotdot'][i]:12.4e} {shipped[i]:12.4e} {fixed[i]:12.4e}")
scale = np.max(np.abs(num))
print(f"max |field - numerical|/max|numerical|     = {np.max(np.abs(C['Hdotdotdot'] - num)) / scale:.3e}")
print(f"max |shipped(py) - numerical|/max|num|     = {np.max(np.abs(shipped - num)) / scale:.3e}")
print(f"max |corrected - numerical|/max|numerical| = {np.max(np.abs(fixed - num)) / scale:.3e}")
print(f"EFTLambdadotdot field (zero = not computed on this path): {C['EFTLambdadotdot'].tolist()}")
