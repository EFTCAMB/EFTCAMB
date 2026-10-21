"""Verifier (ic1-pidot-transient-noisy-coefficients): attribution of the TG_2 (pi, IC1) vs TG_1 (delta_phi) mismatch.
Pairwise lensed-TT/EE max relative differences between existing runs (no CAMB runs)."""

import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[5]
sys.path.insert(0, str(REPO / "fortran" / "eftcamb_test" / "audit" / "checks"))
import harness  # noqa: E402

V = harness.AUDIT_ROOT / "work" / "verify" / "ic1noisy" / "turnon"
T = harness.AUDIT_ROOT / "work" / "t5-robust" / "repro_ic1_pidot_transient"
R = {"tg2_on1e-8": V / "tg2_A0.npz", "tg2_on1.01e-8": V / "tg2_A0_on1.01e-8.npz", "tg2_on2e-8": V / "tg2_A0_on2e-8.npz",
     "tg1_on1e-8": V / "tg1_A0.npz", "tg1_on2e-8": V / "tg1_A0_on2e-8.npz", "tg2_pidot0": T / "pidot0_TG2.npz",
     "tg1_pidot0": T / "pidot0_TG1.npz", "tg2_B1": V / "tg2_B1.npz", "tg1_B1": V / "tg1_B1.npz"}
D = {k: harness.load(v) for k, v in R.items() if v.exists()}


def diff(a, b, j):
    x, y = D[a]["cl_total"][2:2501, j], D[b]["cl_total"][2:2501, j]
    r = x / y - 1
    i = int(np.argmax(np.abs(r)))
    return r[i], i + 2


for a, b in [("tg2_on1e-8", "tg2_on2e-8"), ("tg2_on1e-8", "tg2_on1.01e-8"), ("tg2_on1e-8", "tg1_on1e-8"),
             ("tg2_on2e-8", "tg1_on2e-8"), ("tg2_pidot0", "tg1_pidot0"), ("tg2_on1e-8", "tg2_pidot0"),
             ("tg2_B1", "tg1_B1"), ("tg1_on1e-8", "tg1_B1"), ("tg2_on1e-8", "tg2_B1")]:
    if a in D and b in D:
        (tt, ltt), (ee, lee) = diff(a, b, 0), diff(a, b, 1)
        print(f"{a:14s} vs {b:14s}: dTT/TT {tt:+.3e} (ell {ltt:4d})  dEE/EE {ee:+.3e} (ell {lee:4d})")
# shape of the TG2-TG1 residual
r = D["tg2_on1e-8"]["cl_total"][2:2501, 0] / D["tg1_on1e-8"]["cl_total"][2:2501, 0] - 1
print("TG2/TG1-1 (TT) at ell 2,10,30,60,125,200,400,1000,2000:",
      " ".join(f"{r[l - 2]:+.2e}" for l in (2, 10, 30, 60, 125, 200, 400, 1000, 2000)))
for key in ("pk",):
    if "pk" in D["tg2_on1e-8"]:
        p2, p1 = D["tg2_on1e-8"]["pk"][0], D["tg1_on1e-8"]["pk"][0]
        k = D["tg2_on1e-8"]["pk_k"]
        rr = p2 / p1 - 1
        print("P(k) TG2/TG1-1 max", f"{np.max(np.abs(rr)):.2e} at k={k[np.argmax(np.abs(rr))]:.3g}")
        p0 = D["tg2_pidot0"]["pk"][0]
        print("P(k) TG2(pidot0)/TG1-1 max", f"{np.max(np.abs(p0 / p1 - 1)):.2e}")
