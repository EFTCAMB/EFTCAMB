"""Verifier (ic1-pidot-eftpid-typo): attribute the A0/FIX pi'_ini ratio to the k^4 D1' term alone.

Reads the integrated pidot just after turn-on (a = a_on (1+1e-7)) from ic1typo_attractor.py output. pi_ini is identical
in both builds and E' = 0 in the IC (Gamma1 != 0), so pidot_ini is proportional to Den' with
    FIX: C' + k^2 D1' + k^4 D2'  (D2 = 0 for K-mouflage)      A0: C' + k^2 D1' + k^4 D1'.
Hence R(k) = pidot_A0/pidot_FIX = 1 + k^4/(r + k^2) with r = C'/D1'. Fit r at k=1, predict k=3 and k=10.
"""

import json
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[5]
sys.path.insert(0, str(REPO / "fortran" / "eftcamb_test" / "audit" / "checks"))
import harness  # noqa: E402

s = json.loads((harness.AUDIT_ROOT / "work/verify/ic1typo/attractor/summary_on1e-08.json").read_text())
R = {float(k): s["A0"][k]["pidot0"] / s["FIX"][k]["pidot0"] for k in s["A0"]}
r = 1 / (R[1.0] - 1) - 1  # from R(1) = 1 + 1/(r+1)
print(f"fitted r = C'/D1' = {r:.5f} from R(k=1) = {R[1.0]:.6f}")
for k in sorted(R):
    pred = 1 + k**4 / (r + k**2)
    print(f"k={k:5g}: observed R = {R[k]:.6f}   predicted (typo term only) = {pred:.6f}   rel diff {R[k] / pred - 1:+.2e}")
