"""Planck-sized curvature (|omk|=1e-3) with the EFT perturbation equations active (Pure EFT, constant Omega=1e-6,
i.e. GR to 1e-6) on unpatched A0, vs A0 EFTflag=0 at the same omk."""

import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE))
import harness  # noqa: E402
import t3curv_lens_physics as L  # noqa: E402

OMKS = ("0.001", "-0.001")


def one(omk):
    out = L.OUT / f"A0_PEFTom6_omk{omk}.npz"
    if out.exists():
        return {"status": "ok"}
    ov = {**L.COMMON, "omk": omk, "PureEFTmodelOmega": "1", "EFTOmega0": "1e-6"}
    return harness.run_model("A0", L.PEFT, out, overrides=ov, threads=1, timeout=900)


with ThreadPoolExecutor(2) as ex:
    for s in ex.map(one, OMKS):
        if s["status"] != "ok":
            print("FAILED", s)
keys = ("TT", "TT_l>=30", "EE", "TE", "PP", "pk", "sigma8", "TT_l2_ratio")
print(f"{'case':40s} " + " ".join(f"{k[:9]:>9s}" for k in keys))
for omk in OMKS:
    r = L.row(f"A0_PEFTom6_omk{omk}", f"A0_GR_omk{omk}")
    print(f"{'A0 PEFT Omega=1e-6 omk=' + omk:40s} " + " ".join(f"{r[k]:9.2e}" for k in keys))
