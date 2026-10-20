"""Verifier (ic1-pidot-transient-noisy-coefficients): P(k) of each K-mimic run relative to the reference
(a_ini = 1e-11, turn-on 1e-9, IC1), over the whole k grid: is the deviation broad-band or isolated spikes?"""

import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[5]
sys.path.insert(0, str(REPO / "fortran" / "eftcamb_test" / "audit" / "checks"))
import harness  # noqa: E402

OUT = harness.AUDIT_ROOT / "work" / "verify" / "ic1noisy" / "turnon"
ref = harness.load(OUT / "km_ic1_on1e-9_ai1e-11.npz")
kr = ref["pk_k"]


def onk(d, k):
    return np.exp(np.interp(np.log(k), np.log(d["pk_k"]), np.log(d["pk"][0])))


for f in sorted(OUT.glob("km_*.npz")):
    d = harness.load(f)
    k = d["pk_k"]
    r = d["pk"][0] / onk(ref, k) - 1
    big = np.abs(r) > 0.01
    idx = np.where(big)[0]
    desc = ""
    if len(idx):
        desc = f"first k={k[idx[0]]:.3g}, last k={k[idx[-1]]:.3g}, argmax k={k[np.argmax(np.abs(r))]:.3g}"
    print(f"{f.stem:28s} n_k={len(k)} n(|dP/P|>1%)={big.sum():3d} max|dP/P|={np.max(np.abs(r)):.3e}  {desc}")
    if len(sys.argv) > 1 and sys.argv[1] in f.stem:
        for i in range(len(k)):
            if k[i] > 1e-3:
                print(f"   k={k[i]:.4g} dP/P={r[i]:+.3e}")
