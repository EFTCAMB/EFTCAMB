"""Baseline for t3curv_lens_physics.py: A0 EFTflag=0 at the same omk, so the diagnostic-patch residual is measured
against the same tree's GR branch (removes the known A0-vs-B0 EFTflag=0 offset)."""

import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE))
import harness  # noqa: E402
import t3curv_lens_physics as L  # noqa: E402

OMKS = ["0", "0.001", "-0.001", "0.01", "0.05", "-0.05"]


def one(omk):
    out = L.OUT / f"A0_GR_omk{omk}.npz"
    if out.exists():
        return {"status": "ok"}
    return harness.run_model("A0", L.GR, out, overrides={**L.COMMON, "omk": omk}, threads=1, timeout=900)


with ThreadPoolExecutor(4) as ex:
    for s in ex.map(one, OMKS):
        if s["status"] != "ok":
            print("FAILED", s)
keys = ("TT", "TT_l>=30", "EE", "TE", "PP", "pk", "sigma8", "bg_H", "TT_l2_ratio")
print(f"{'case':52s} " + " ".join(f"{k[:9]:>9s}" for k in keys))
for omk in OMKS:
    for lab, a in (("unpatched A0 PEFT0", f"A0_PEFT0_omk{omk}"), ("diag-patched PEFT0", f"D_PEFT0_omk{omk}"),
                   ("A0 GR vs B0 GR", None)):
        if a is None:
            r = L.row(f"A0_GR_omk{omk}", f"B0_GR_omk{omk}")
            name = f"{lab} omk={omk}"
        else:
            if not (L.OUT / f"{a}.npz").exists():
                continue
            r = L.row(a, f"A0_GR_omk{omk}")
            name = f"{lab} vs A0 GR omk={omk}"
        print(f"{name:52s} " + " ".join(f"{r[k]:9.2e}" if r.get(k) is not None else f"{'-':>9s}" for k in keys))
