"""Verifier (ic1-pidot-transient-noisy-coefficients): re-analyse already produced npz (no CAMB runs).
Reads the reporter's repro outputs (work/t5-robust/repro_ic1_pidot_transient) and the verifier turn-on runs
(work/verify/ic1noisy/turnon)."""

import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[5]
sys.path.insert(0, str(REPO / "fortran" / "eftcamb_test" / "audit" / "checks"))
import harness  # noqa: E402

KS = (0.1, 0.5, 1.0, 2.0, 5.0, 10.0)


def pk_at(d, k):
    return float(np.exp(np.interp(np.log(k), np.log(d["pk_k"]), np.log(d["pk"][0]))))


def cl_diff(a, b):
    out = []
    for j in (0, 1):
        x, y = a["cl_total"][2:2501, j], b["cl_total"][2:2501, j]
        r = np.abs(x / y - 1)
        out.append((float(r.max()), int(np.argmax(r)) + 2))
    return out


for d in (harness.AUDIT_ROOT / "work" / "t5-robust" / "repro_ic1_pidot_transient",
          harness.AUDIT_ROOT / "work" / "verify" / "ic1noisy" / "turnon"):
    print(f"### {d}")
    gr = harness.load(d / "GR.npz")
    print("  P/P_GR-1 at k = " + ", ".join(f"{k:g}" for k in KS) + "   (pk_k max = %.1f)" % gr["pk_k"].max())
    for f in sorted(d.glob("*.npz")):
        if "Kmimic" not in f.name and not f.name.startswith("km"):
            continue
        x = harness.load(f)
        print(f"  {f.stem:28s} " + " ".join(f"{pk_at(x, k) / pk_at(gr, k) - 1:+.3e}" for k in KS)
              + f"  s8={float(x['sigma8'][0]):.4f}")
    for p, q in (("A0_TG2", "A0_TG1"), ("pidot0_TG2", "pidot0_TG1"), ("pidot0_TG1", "A0_TG1")):
        if (d / f"{p}.npz").exists() and (d / f"{q}.npz").exists():
            (tt, ltt), (ee, lee) = cl_diff(harness.load(d / f"{p}.npz"), harness.load(d / f"{q}.npz"))
            print(f"  {p} vs {q}: TT {tt:.2e} (ell {ltt}) EE {ee:.2e} (ell {lee})")
