"""Check 3: non-linear models (Halofit/HMcode, SPk) with EFT on.

T1: every model x halofit_version runs with NonLinear_both and gives a finite non-linear ratio in a sane range.
T2: consistency in the no-MG limit with non-LCDM background: pure-EFT wCDM (2_PEFT_wDE_wCDM_1, w=-0.8, no
    modification of gravity) vs CAMB fluid wCDM (EFTflag=0, w=-0.8, cs2=1). Linear P(k) agree; the non-linear ratio
    P_nl/P_lin must agree too if Halofit/HMcode see the right w (EFTCAMB passes Effective_w_wa, which is (-1, 0)
    unless EFTCAMB_effective_w0wa=T and the model implements it -- only quintessence does).
T3: quintessence with EFTCAMB_effective_w0wa=T vs F.
"""

from __future__ import annotations

import json
import sys

import numpy as np
from t8c_lib import WORK, load, model_ini, run_many

OUT = WORK / "c3"
MODELS = ["1_EFT_GR", "2_PEFT_wDE_wCDM_1", "4_Des_F_R_LCDM_1", "5_hdsk_TG_1", "5_Kmouflage_1", "5_quint_exponential",
          "2_PEFT_gamma3_const_1"]
VERSIONS = ["mead2020", "takahashi", "mead"]
OV = {"get_tensor_cls": "F", "feedback_level": "0", "num_redshiftwindows": "0"}
MP = {"redshifts": [2.0, 1.0, 0.5, 0.0], "kmax": 10.0, "nonlinear": True}


def args(version):
    return {"matter_power": MP, "NonLinear": "NonLinear_both", "nonlinear_model": {"halofit_version": version},
            "set_for_lmax": {"lmax": 2500, "lens_potential_accuracy": 1},
            "outputs": ["cmb", "pk_lin", "pk_nl", "eftw0wa"], "cl_lmax": 2500}


def jobs():
    js = []
    for v in VERSIONS:
        for m in MODELS:
            js.append(dict(tree="A0", ini=model_ini(m), task="generic", out=OUT / f"{m}_{v}.npz", overrides=OV,
                           args=args(v)))
        js.append(dict(tree="A0", ini=model_ini("1_EFT_GR"), task="generic", out=OUT / f"GRw08_{v}.npz",
                       overrides=OV | {"w": "-0.8", "cs2_lam": "1"}, args=args(v)))
        js.append(dict(tree="A0", ini=model_ini("5_quint_exponential"), task="generic", out=OUT / f"quintW_{v}.npz",
                       overrides=OV | {"EFTCAMB_effective_w0wa": "T"}, args=args(v)))
    return js


def ratio(d):
    return d["pknl"] / d["pk"]


def report():
    res = {"T1": {}, "T2": {}, "T3": {}}
    for v in VERSIONS:
        for m in MODELS + ["GRw08", "quintW"]:
            try:
                d = load(OUT / f"{m}_{v}.npz")
            except FileNotFoundError:
                res["T1"][f"{m}:{v}"] = "missing"
                continue
            r = ratio(d)
            k = d["pk_k"]
            sel = k < 10
            res["T1"][f"{m}:{v}"] = {"finite": bool(np.all(np.isfinite(r))), "min": float(r[:, sel].min()),
                                     "max": float(r[:, sel].max()), "w0wa": d["meta"].get("eft_w0wa")}
        try:
            a, b = load(OUT / f"2_PEFT_wDE_wCDM_1_{v}.npz"), load(OUT / f"GRw08_{v}.npz")
            k = a["pk_k"]
            sel = (k > 0.05) & (k < 5)
            lin = np.abs(a["pk"][:, sel] / b["pk"][:, sel] - 1).max(axis=1)
            rat = np.abs(ratio(a)[:, sel] / ratio(b)[:, sel] - 1).max(axis=1)
            res["T2"][v] = {"z": list(a["pk_z"]), "lin_maxrel": lin.tolist(), "nlratio_maxrel": rat.tolist(),
                            "w0wa_eft": a["meta"].get("eft_w0wa")}
        except FileNotFoundError:
            res["T2"][v] = "missing"
        try:
            a, b = load(OUT / f"quintW_{v}.npz"), load(OUT / f"5_quint_exponential_{v}.npz")
            k = a["pk_k"]
            sel = (k > 0.05) & (k < 5)
            res["T3"][v] = {"nlratio_maxrel": np.abs(ratio(a)[:, sel] / ratio(b)[:, sel] - 1).max(axis=1).tolist(),
                            "w0wa_on": a["meta"].get("eft_w0wa"), "w0wa_off": b["meta"].get("eft_w0wa")}
        except FileNotFoundError:
            res["T3"][v] = "missing"
    print(json.dumps(res, indent=1))
    return res


if __name__ == "__main__":
    cmd = sys.argv[1] if len(sys.argv) > 1 else "run"
    if cmd == "run":
        for r in run_many(jobs(), nproc=int(sys.argv[2]) if len(sys.argv) > 2 else 4):
            print(r["out"].split("/")[-1], r["status"], r["elapsed"])
    report()
