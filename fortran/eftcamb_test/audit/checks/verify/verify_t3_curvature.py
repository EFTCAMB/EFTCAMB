"""Verifier for t3-eft-curvature-ignored: GR-equivalent Pure EFT vs EFTflag=0 at omk!=0, on A0 and B1, same tree."""

import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import harness  # noqa: E402

OUT = harness.AUDIT_ROOT / "work" / "verify" / "curv"
PEFT = harness.PARAMS_DIR / "2_PEFT_Omega_const_1.ini"
GR = harness.PARAMS_DIR / "1_EFT_GR.ini"
COMMON = {"num_redshiftwindows": "0", "get_tensor_cls": "F", "l_max_scalar": "2200"}


def job(a):
    tree, ini, ov, name = a
    return harness.run_model(tree, ini, OUT / f"{name}.npz", overrides=ov, threads=1)


jobs = []
for tree in ("A0", "B1"):
    for omk in ("0.01", "0.05"):
        ov = {**COMMON, "omk": omk}
        jobs.append((tree, PEFT, {**ov, "PureEFTmodelOmega": "0", "EFTCAMB_skip_RGR": "F"}, f"{tree}_PEFT_omk{omk}"))
        jobs.append((tree, GR, ov, f"{tree}_GR_omk{omk}"))
with ThreadPoolExecutor(4) as ex:
    st = list(ex.map(job, jobs))
for s in st:
    print(s["status"], s.get("out"))
for tree in ("A0", "B1"):
    for omk in ("0.01", "0.05"):
        a = harness.load(OUT / f"{tree}_PEFT_omk{omk}.npz")
        b = harness.load(OUT / f"{tree}_GR_omk{omk}.npz")
        r = harness.compare(a, b)
        t = r["cl"]["total"]
        print(tree, omk, {k: f"{v:.2e}" for k, v in t.items()}, "pk", f"{r['worst_pk']:.2e}", "s8",
              f"{r['sigma8']:.2e}", "H", f"{r['bg_H']:.2e}")
        ca, cb = a["cl_total"][:, 0], b["cl_total"][:, 0]
        for ell in (2, 10, 30, 100, 220, 1000, 2000):
            print(f"   ell={ell} TT ratio {ca[ell] / cb[ell]:.4f}")
