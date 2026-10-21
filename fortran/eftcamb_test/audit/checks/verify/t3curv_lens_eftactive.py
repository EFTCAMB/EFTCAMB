"""EFT-perturbations-active GR limit: Pure EFT with all EFT functions zero and EFTCAMB_skip_RGR=T, so derivs runs the EFT
Einstein equations from a=0.01. Unpatched A0 and the diagnostic tree (background curvature restored) vs A0 EFTflag=0.
A residual in the diagnostic tree at omk/=0 measures curvature missing from the EFT perturbation equations."""

import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE))
import harness  # noqa: E402
import t3curv_lens_physics as L  # noqa: E402

# variant "act": exact GR limit with RGR skipped (mass-stability check off: pi has no kinetic term in exact GR)
# variant "om6": constant Omega=1e-6 (> RGR threshold 1e-8, so the EFT equations are active from a=0.01)
JOBS = []
VAR = {
    "act": {"PureEFTmodelOmega": "0", "EFTCAMB_skip_RGR": "T", "EFT_mass_stability": "F"},
    "om6": {"PureEFTmodelOmega": "1", "EFTOmega0": "1e-6"},
}
OMKS = ("0", "0.01", "0.05", "-0.05")
for v, extra in VAR.items():
    for omk in OMKS:
        ov = {**L.COMMON, "omk": omk, **extra}
        JOBS.append(("A0", ov, f"A0_PEFT0{v}_omk{omk}"))
        JOBS.append((str(L.DIAG), ov, f"D_PEFT0{v}_omk{omk}"))


def one(j):
    tree, ov, name = j
    out = L.OUT / f"{name}.npz"
    if out.exists():
        return {"status": "ok"}
    return harness.run_model(tree, L.PEFT, out, overrides=ov, threads=1, timeout=900)


with ThreadPoolExecutor(4) as ex:
    for s in ex.map(one, JOBS):
        if s["status"] != "ok":
            print("FAILED", s.get("out"), s["status"], s.get("log_tail", [])[-4:])
keys = ("TT", "TT_l>=30", "EE", "TE", "PP", "pk", "sigma8", "TT_l2_ratio")
print(f"{'case':52s} " + " ".join(f"{k[:9]:>9s}" for k in keys))
for v in VAR:
    for omk in OMKS:
        for lab, a in ((f"A0 PEFT {v}", f"A0_PEFT0{v}_omk{omk}"), (f"diag PEFT {v}", f"D_PEFT0{v}_omk{omk}")):
            if not (L.OUT / f"{a}.npz").exists():
                print(f"{lab} omk={omk}: no output")
                continue
            r = L.row(a, f"A0_GR_omk{omk}")
            print(f"{lab + ' vs A0 GR omk=' + omk:52s} " + " ".join(f"{r[k]:9.2e}" for k in keys))
