"""Check 2b: are the EFT-modified potentials used in transfer outputs?  Weyl transfer consistency.

C_L^phiphi from CAMB's line-of-sight lensing source (EFT: -EFTLensing/k) vs a Limber integral of
P_Weyl(k,z) = get_matter_power_interpolator(var1=var2=Transfer_Weyl) (EFT: Transfer_Weyl = k EFTLensing/2), CAMBdemo
recipe, linear.  For a model with Sigma != 1 (PEFT Omega const: lensing ~ 1/(1+Omega)) a GR expression in either path
would show up as an O(Omega) mismatch; GR sets the baseline Limber accuracy.  Compared over 100 <= L <= 1500.
"""

from __future__ import annotations

import json
import sys

import numpy as np
from t8c_lib import WORK, load, model_ini, run_many

OUT = WORK / "c2b"
OV = {"get_tensor_cls": "F", "feedback_level": "0", "num_redshiftwindows": "0", "get_transfer": "F"}
MODELS = ["1_EFT_GR", "2_PEFT_Omega_const_1", "4_Des_F_R_LCDM_1", "5_hdsk_TG_1", "2_PEFT_gamma3_const_1"]
ARGS = {"set_for_lmax": {"lmax": 2000, "lens_potential_accuracy": 2, "nonlinear": False},
        "outputs": ["weyl_limber"]}


def jobs():
    return [dict(tree="A0", ini=model_ini(m), task="generic", out=OUT / f"{m}.npz", overrides=OV, args=ARGS)
            for m in MODELS]


def report():
    res = {}
    for m in MODELS:
        try:
            d = load(OUT / f"{m}.npz")
        except FileNotFoundError:
            st = json.loads((OUT / f"{m}.json").read_text())
            res[m] = st["status"] + " " + str(st.get("log_tail", [])[-2:])
            continue
        ls = d["limber_ls"]
        r = d["limber_pp"] / d["los_pp"] - 1
        sel = (ls >= 100) & (ls <= 1500)
        res[m] = {"max_abs_rel_100_1500": float(np.max(np.abs(r[sel]))), "mean_rel_100_1500": float(np.mean(r[sel])),
                  "rel_at": {int(L): float(r[int(L) - 2]) for L in (100, 300, 1000, 1500)}}
    print(json.dumps(res, indent=1))
    return res


if __name__ == "__main__":
    cmd = sys.argv[1] if len(sys.argv) > 1 else "run"
    if cmd == "run":
        for r in run_many(jobs(), nproc=int(sys.argv[2]) if len(sys.argv) > 2 else 4):
            print(r["out"].split("/")[-1], r["status"], r["elapsed"])
    report()
