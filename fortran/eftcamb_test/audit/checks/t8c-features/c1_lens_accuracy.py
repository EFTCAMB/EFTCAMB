"""Check 1: v2 automatic lens_potential_accuracy (None -> max(4, (lmax-1500)/500)) with EFT models.

For each model: set_for_lmax(lmax, lens_potential_accuracy=None, nonlinear=False) vs lens_potential_accuracy=8
(linear lensing in both, so only the k range differs).  Also runs the v2 default (nonlinear=None), which switches on
non-linear lensing.  GR (EFTflag=0) is the calibration baseline: an MG model is "under-resolved" if its auto-vs-8
difference is much larger than GR's.

usage: python c1_lens_accuracy.py [run|report] [lmax]
"""

from __future__ import annotations

import json
import sys

import numpy as np
from t8c_lib import WORK, cl_diffs, load, maxrel, model_ini, run_many

MODELS = ["1_EFT_GR", "4_Des_F_R_LCDM_1", "4_Des_F_R_wCDM_1", "5_hdsk_TG_1", "5_Kmouflage_1", "5_quint_exponential",
          "2_PEFT_gamma3_const_1", "2_PEFT_Omega_const_1"]
OV = {"get_tensor_cls": "F", "get_transfer": "F", "feedback_level": "0", "num_redshiftwindows": "0"}
OUT = WORK / "c1"


def jobs(lmax):
    js = []
    for m in MODELS:
        for tag, lpa, nl in (("auto", None, False), ("lpa8", 8, False), ("autoNL", None, None)):
            js.append(dict(tree="A0", ini=model_ini(m), task="generic", out=OUT / f"{m}_{tag}_{lmax}.npz",
                           overrides=OV, args={"set_for_lmax": {"lmax": lmax, "lens_potential_accuracy": lpa,
                                                                 "nonlinear": nl},
                                               "outputs": ["cmb"], "cl_lmax": lmax}))
    return js


def report(lmax):
    rows = {}
    for m in MODELS:
        try:
            a = load(OUT / f"{m}_auto_{lmax}.npz")
            b = load(OUT / f"{m}_lpa8_{lmax}.npz")
        except FileNotFoundError:
            rows[m] = "missing"
            continue
        d = cl_diffs(a["cl_total"], b["cl_total"], lmax=lmax)
        d["PP"] = maxrel(a["cl_lens_potential"][2:lmax + 1, 0], b["cl_lens_potential"][2:lmax + 1, 0])
        d["kmax_eta_auto"] = a["meta"]["max_eta_k"]
        d["kmax_eta_8"] = b["meta"]["max_eta_k"]
        try:
            c = load(OUT / f"{m}_autoNL_{lmax}.npz")
            dn = cl_diffs(c["cl_total"], a["cl_total"], lmax=lmax)
            d["NL_vs_lin"] = dn
            d["NonLinear_default"] = c["meta"]["NonLinear"]
        except FileNotFoundError:
            pass
        rows[m] = d
    print(json.dumps(rows, indent=1))
    return rows


if __name__ == "__main__":
    cmd = sys.argv[1] if len(sys.argv) > 1 else "run"
    lmax = int(sys.argv[2]) if len(sys.argv) > 2 else 3000
    if cmd == "run":
        sel = sys.argv[3].split(",") if len(sys.argv) > 3 else None
        js = [j for j in jobs(lmax) if sel is None or any(s in str(j["out"]) for s in sel)]
        for r in run_many(js, nproc=4):
            print(r["out"].split("/")[-1], r["status"], r["elapsed"])
    report(lmax)
