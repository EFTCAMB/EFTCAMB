"""Repro: the third conformal-time derivative of the Hubble rate (timestep cache field Hdotdotdot) computed by
EFTCAMBFullModelComputeHubbleDer (fortran/eftcamb/06p1_abstract_EFTCAMB_full_map.f90:145) is not the derivative of
Hdotdot for full-mapping models that do not override compute_H_derivs (K-mouflage, FM quintessence, ADE).
The designer f(R) model never sets it (always 0). Hdot and Hdotdot are shown as controls.

The numerical derivative uses a Richardson-extrapolated central difference in a of the Hdotdot field returned by
get_eft_functions, converted to conformal time with d/dtau = a*adotoa d/da.

Usage: python fullmap-hdotdotdot-wrong.py [--tree PATH]
"""

import argparse
import os
import sys
import tempfile
from pathlib import Path

import numpy as np

ROOT = "/private/tmp/claude-501/-Users-marco-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad/audit"
PARAMS = Path(__file__).resolve().parents[2] / "parameters"

ap = argparse.ArgumentParser()
ap.add_argument("--tree", default=f"{ROOT}/trees/A0")
ap.add_argument("--models", default="5_Kmouflage_1,5_quint_exponential,5_ADE_1,4_Des_F_R_LCDM_1")
args = ap.parse_args()
os.environ.setdefault("OMP_NUM_THREADS", "2")
sys.path.insert(0, args.tree)
import camb  # noqa: E402

print("camb from", camb.__file__)
tmp = tempfile.mkdtemp(prefix="h3_")
a0 = np.array([0.05, 0.2, 0.5, 0.9])


def deriv_tau(cp, res, field, h):
    _, P = cp.EFTCAMB.get_eft_functions(res, a0 * (1 + h))
    _, M = cp.EFTCAMB.get_eft_functions(res, a0 * (1 - h))
    return (P[field] - M[field]) / (2 * h * a0)


for m in args.models.split(","):
    ini = os.path.join(tmp, f"{m}.ini")
    Path(ini).write_text(f"DEFAULT({PARAMS / (m + '.ini')})\noutput_root = {tmp}/{m}\n"
                         "EFTCAMB_write_background = F\nfeedback_level = 0\n")
    cp = camb.read_ini(ini, no_validate=True)
    res = camb.get_background(cp)
    _, A = cp.EFTCAMB.get_eft_functions(res, a0)
    conv = a0 * A["adotoa"]  # d/dtau = a adotoa d/da
    print(f"\n== {m} (EFTflag={cp.EFTCAMB.EFTflag}, model={cp.EFTCAMB.model_name()})")
    for f, df in (("adotoa", "Hdot"), ("Hdot", "Hdotdot"), ("Hdotdot", "Hdotdotdot")):
        h = 1e-3
        num = (4 * deriv_tau(cp, res, f, h / 2) - deriv_tau(cp, res, f, h)) / 3 * conv
        rel = np.abs(A[df] - num) / np.maximum(np.abs(num), 1e-300)
        print(f"  {df:11s} field   : " + " ".join(f"{v: .4e}" for v in A[df]))
        print(f"  {'':11s} d{f}/dtau: " + " ".join(f"{v: .4e}" for v in num) + f"   max rel diff = {rel.max():.2e}")
print("\nexpected: field == numerical derivative to <~1e-3 for all three rows; a = ", a0.tolist())
