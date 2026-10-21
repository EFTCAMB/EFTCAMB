"""Check 10: changed v2 accuracy defaults (AccuracyTarget=1 default, AccuracyTarget=0 ~ CAMB 1.x) with EFT on.

For each model, AccuracyTarget=0 vs default (1): max fractional CMB/P(k) change.  The v2 targeted accuracy changes
are model independent; an EFT model whose change is much larger than GR's flags EFT sensitivity to the new defaults.
Also a convergence reference (AccuracyBoost=2, lAccuracyBoost=2) to see which setting is closer.
"""

from __future__ import annotations

import json
import sys

from t8c_lib import WORK, cl_diffs, load, maxrel, model_ini, pk_diff, run_many

OUT = WORK / "c10"
MODELS = ["1_EFT_GR", "4_Des_F_R_LCDM_1", "5_hdsk_TG_1", "5_quint_exponential", "2_PEFT_gamma3_const_1",
          "5_Kmouflage_1"]
OV = {"get_tensor_cls": "F", "feedback_level": "0", "num_redshiftwindows": "0", "transfer_kmax": "2"}
# AccuracyTarget is a global (camb.config), read from an ini only by the command-line driver: set it from Python
VARIANTS = {"AT1": ({}, {}), "AT0": ({}, {"config": {"AccuracyTarget": 0}}),
            "boost2": ({"accuracy_boost": "2", "l_accuracy_boost": "2"}, {})}
ARGS = {"set_for_lmax": {"lmax": 2500, "lens_potential_accuracy": 1, "nonlinear": False},
        "outputs": ["cmb", "pk_lin"], "cl_lmax": 2500}


def jobs():
    return [dict(tree="A0", ini=model_ini(m), task="generic", out=OUT / f"{m}_{v}.npz", overrides=OV | ov,
                 args=ARGS | py) for m in MODELS for v, (ov, py) in VARIANTS.items()]


def report():
    res = {}
    for m in MODELS:
        try:
            a1, a0, b = (load(OUT / f"{m}_{v}.npz") for v in ("AT1", "AT0", "boost2"))
        except FileNotFoundError:
            res[m] = "missing"
            continue
        res[m] = {
            "AT0_vs_AT1": cl_diffs(a0["cl_total"], a1["cl_total"], lmax=2500) | {"pk": pk_diff(a0, a1)},
            "AT1_vs_boost2": cl_diffs(a1["cl_total"], b["cl_total"], lmax=2500) | {"pk": pk_diff(a1, b)},
            "AT0_vs_boost2": cl_diffs(a0["cl_total"], b["cl_total"], lmax=2500) | {"pk": pk_diff(a0, b)},
            "time": [a0["meta"]["time_results"], a1["meta"]["time_results"]],
        }
    print(json.dumps(res, indent=1))
    return res


if __name__ == "__main__":
    cmd = sys.argv[1] if len(sys.argv) > 1 else "run"
    if cmd == "run":
        for r in run_many(jobs(), nproc=int(sys.argv[2]) if len(sys.argv) > 2 else 4):
            print(r["out"].split("/")[-1], r["status"], r["elapsed"])
    report()
