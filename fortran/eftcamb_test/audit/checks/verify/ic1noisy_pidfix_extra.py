"""Verifier (ic1-pidot-transient-noisy-coefficients): separate this finding from the EFTpiD2 typo finding:
5_Kmimic_1 on the typo-fixed tree (work/t5-robust/A0_pidfix, pi' still from dfridr) at turn-on 1e-9 (= table start)
and 1e-7 (table node), plus A0 at 1e-9 with model_background_a_ini = 1e-9*0.9 (stencil inside the table)."""

import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[5]
sys.path.insert(0, str(REPO / "fortran" / "eftcamb_test" / "audit" / "checks"))
import harness  # noqa: E402

OUT = harness.AUDIT_ROOT / "work" / "verify" / "ic1noisy" / "turnon"
PIDFIX = str(harness.AUDIT_ROOT / "work" / "t5-robust" / "A0_pidfix")
KM = "5_Kmimic_1"
KS = (0.1, 0.5, 1.0, 2.0, 5.0, 10.0)
JOBS = [("km_pidfix_on1e-9", PIDFIX, {"EFTCAMB_turn_on_time": "1e-9"}),
        ("km_pidfix_on1e-7", PIDFIX, {"EFTCAMB_turn_on_time": "1e-7"}),
        ("km_ic1_on1e-9_ai9e-10", "A0", {"EFTCAMB_turn_on_time": "1e-9", "model_background_a_ini": "9e-10"})]


def job(j):
    name, tree, ov = j
    out = OUT / f"{name}.npz"
    if not out.exists():
        st = harness.run_model(tree, harness.PARAMS_DIR / f"{KM}.ini", out, overrides=ov, threads=1, timeout=600)
        if st["status"] != "ok":
            return name, st["status"]
    return name, harness.load(out)


def pk_at(d, k):
    return float(np.exp(np.interp(np.log(k), np.log(d["pk_k"]), np.log(d["pk"][0]))))


with ThreadPoolExecutor(3) as ex:
    R = dict(ex.map(job, JOBS))
gr = harness.load(OUT / "GR.npz")
print("P/P_GR-1 at k = " + ", ".join(map(str, KS)))
for name, d in R.items():
    if isinstance(d, str):
        print(name, d)
        continue
    print(f"  {name:24s} " + " ".join(f"{pk_at(d, k) / pk_at(gr, k) - 1:+.3e}" for k in KS) +
          f"  s8={float(d['sigma8'][0]):.4f}")
