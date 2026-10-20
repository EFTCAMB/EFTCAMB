"""Pre-rebase EFTCAMB (B1, CAMB 1.6.5): GR-limit Pure EFT vs B1 EFTflag=0 at omk=0 and 0.05 (is it a rebase regression?)."""

import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE))
import harness  # noqa: E402
import t3curv_lens_physics as L  # noqa: E402

JOBS = []
for omk in ("0", "0.05"):
    JOBS.append((L.PEFT, {**L.COMMON, "omk": omk, "PureEFTmodelOmega": "0"}, f"B1_PEFT0_omk{omk}"))
    JOBS.append((L.GR, {**L.COMMON, "omk": omk}, f"B1_GR_omk{omk}"))


def one(j):
    ini, ov, name = j
    out = L.OUT / f"{name}.npz"
    if out.exists():
        return {"status": "ok"}
    return harness.run_model("B1", ini, out, overrides=ov, threads=1, timeout=900)


with ThreadPoolExecutor(4) as ex:
    for s in ex.map(one, JOBS):
        if s["status"] != "ok":
            print("FAILED", s)
keys = ("TT", "TT_l>=30", "EE", "TE", "PP", "pk", "sigma8", "TT_l2_ratio")
print(f"{'case':40s} " + " ".join(f"{k[:9]:>9s}" for k in keys))
for omk in ("0", "0.05"):
    r = L.row(f"B1_PEFT0_omk{omk}", f"B1_GR_omk{omk}")
    print(f"{'B1 PEFT0 vs B1 GR omk=' + omk:40s} " + " ".join(f"{r[k]:9.2e}" for k in keys))
