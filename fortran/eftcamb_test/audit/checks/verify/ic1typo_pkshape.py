"""Verifier (ic1-pidot-eftpid-typo): print P(k)/P_GR - 1 on the native k grid for cached turn-on runs (A0 vs FIX),
to see whether the large ratios are smooth physics or spiky numerical features."""

import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[5]
sys.path.insert(0, str(REPO / "fortran" / "eftcamb_test" / "audit" / "checks"))
import harness  # noqa: E402

OUT = harness.AUDIT_ROOT / "work" / "verify" / "ic1typo" / "turnon"
model = sys.argv[1] if len(sys.argv) > 1 else "5_Kmimic_1"
ton = sys.argv[2] if len(sys.argv) > 2 else "1.e-8"
gr = harness.load(OUT / "GR_A0.npz")
a = harness.load(OUT / f"A0_{model}_on{ton}.npz")
f = harness.load(OUT / f"FIX_{model}_on{ton}.npz")
k = a["pk_k"]
print("A0/FIX native k grids identical:", np.array_equal(a["pk_k"], f["pk_k"]), " GR grid len", len(gr["pk_k"]))


def onk(d):
    return np.exp(np.interp(np.log(k), np.log(d["pk_k"]), np.log(d["pk"][0])))


ra = a["pk"][0] / onk(gr) - 1
rf = onk(f) / onk(gr) - 1
print(f"{model} turn_on={ton}: nk={len(k)} kmax={k.max():.3g}")
for i in range(0, len(k), max(1, len(k) // 80)):
    print(f"  k={k[i]:9.4g}  A0 {ra[i]:+.4e}  FIX {rf[i]:+.4e}")
