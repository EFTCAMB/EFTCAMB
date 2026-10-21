"""Verifier (ic1-pidot-transient-noisy-coefficients): for the K-mimic configurations that misbehave on A0 (turn-on near a
background grid node, or 5000-point grid), does pi'_ini = 0 (t5's A0_icsw tree, EFT_T5_PIDOT0=1, which also has the
EFTpiD2 typo fixed) remove the problem?  Compares sigma8 and P/P_GR at k = 0.1, 0.5, 1, 2 /Mpc."""

import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[5]
sys.path.insert(0, str(REPO / "fortran" / "eftcamb_test" / "audit" / "checks"))
import harness  # noqa: E402

ICSW = harness.AUDIT_ROOT / "work" / "t5-robust" / "A0_icsw"
OUT = harness.AUDIT_ROOT / "work" / "verify" / "ic1noisy" / "turnon"
KS = (0.1, 0.5, 1.0, 2.0)
JOBS = [
    ("km_pidot0_on1e-8_np5000", {"EFTCAMB_turn_on_time": "1e-8", "model_background_num_points": "5000"}),
    ("km_pidot0_on3e-8", {"EFTCAMB_turn_on_time": "3e-8"}),
    ("km_pidot0_on1.003e-8", {"EFTCAMB_turn_on_time": "1.003e-8"}),
    ("km_pidot0_on1e-8", {"EFTCAMB_turn_on_time": "1e-8"}),
]


def job(a):
    name, ov = a
    out = OUT / f"{name}.npz"
    if not out.exists():
        st = harness.run_model(ICSW, harness.PARAMS_DIR / "5_Kmimic_1.ini", out, overrides=ov, threads=1,
                               timeout=900, extra_env={"EFT_T5_PIDOT0": "1"})
        if st["status"] != "ok":
            return name, st["status"]
    return name, harness.load(out)


def pk_at(d, k):
    return float(np.exp(np.interp(np.log(k), np.log(d["pk_k"]), np.log(d["pk"][0]))))


with ThreadPoolExecutor(4) as ex:
    res = dict(ex.map(job, JOBS))
gr = harness.load(OUT / "GR.npz")
ref = harness.load(OUT / "km_ic1_on1e-9_ai1e-11.npz")
for name, d in res.items():
    if isinstance(d, str):
        print(name, d)
        continue
    k = d["pk_k"]
    r = d["pk"][0] / np.exp(np.interp(np.log(k), np.log(ref["pk_k"]), np.log(ref["pk"][0]))) - 1
    print(f"{name:28s} " + "  ".join(f"{pk_at(d, kk) / pk_at(gr, kk) - 1:+.4e}" for kk in KS)
          + f"   s8 {float(d['sigma8'][0]):.5f}   max|P/P_ref-1| (k<=10) {np.max(np.abs(r[k <= 10])):.3e}")
