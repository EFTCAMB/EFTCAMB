"""Verifier (ic1-pidot-transient-noisy-coefficients): compare early-time EFT background coefficients (get_eft_functions)
between A0 and B1 reference runs, and measure their smoothness (number of slope sign changes in a in [1e-8, 1e-5])."""

import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[5]
sys.path.insert(0, str(REPO / "fortran" / "eftcamb_test" / "audit" / "checks"))
import harness  # noqa: E402

R = harness.AUDIT_ROOT / "runs"
FIELDS = ["EFTOmegaV", "EFTc", "EFTLambda", "EFTGamma1V", "EFTGamma2V", "EFTGamma3V", "EFTGamma4V", "adotoa", "Hdot",
          "alphaB", "alphaK", "Meff2"]
for model in sys.argv[1:] or ["5_hdsk_TG_2", "5_Kmimic_1"]:
    print(f"=== {model}")
    da, db = harness.load(R / "A0" / f"{model}.npz"), harness.load(R / "B1" / f"{model}.npz")
    a = da["eft_a_grid"]
    m = (a >= 1e-8) & (a <= 1e-5)
    for f in FIELDS:
        x, y = da[f"eft_{f}"][m], db[f"eft_{f}"][m]
        nsa = int(np.sum(np.diff(np.sign(np.diff(x))) != 0))
        nsb = int(np.sum(np.diff(np.sign(np.diff(y))) != 0))
        rel = np.max(np.abs(x - y)) / max(np.max(np.abs(y)), 1e-300)
        print(f"  {f:11s} A0 range [{x.min():+.3e},{x.max():+.3e}] slope-sign-changes A0 {nsa:3d} B1 {nsb:3d} "
              f"max|A0-B1|/max|B1| = {rel:.2e}")
