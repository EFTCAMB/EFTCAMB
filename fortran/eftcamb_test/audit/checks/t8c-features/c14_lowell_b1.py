"""Check 14: is the low-ell TT sensitivity to accuracy settings of 2_PEFT_gamma3_const_1 (c10: default vs
accuracy_boost=2,l_accuracy_boost=2 differs by 9e-3 at l~13 while GR differs by 8e-4) new in the rebased code (A0)
or already present in the pre-rebase EFTCAMB (B1)?  Ini path (harness.run_model), lmax 2500.

usage: python c14_lowell_b1.py
"""

from __future__ import annotations

import json
from concurrent.futures import ThreadPoolExecutor

import numpy as np
from t8c_lib import WORK, harness, load, model_ini

OUT = WORK / "c14"
OV = {"get_tensor_cls": "F", "feedback_level": "0", "num_redshiftwindows": "0", "get_transfer": "F",
      "l_max_scalar": "2500"}
BOOST = {"accuracy_boost": "2", "l_accuracy_boost": "2"}
MODELS = ["2_PEFT_gamma3_const_1", "1_EFT_GR"]


def one(job):
    tree, m, tag, ov = job
    return harness.run_model(tree, model_ini(m), OUT / f"{tree}_{m}_{tag}.npz", overrides=OV | ov, threads=1)


jobs = [(t, m, tag, ov) for t in ("A0", "B1") for m in MODELS for tag, ov in (("def", {}), ("boost2", BOOST))]
with ThreadPoolExecutor(4) as ex:
    for j, r in zip(jobs, ex.map(one, jobs)):
        print(j[:3], r["status"], r.get("elapsed"))
res = {}
for t in ("A0", "B1"):
    for m in MODELS:
        a, b = load(OUT / f"{t}_{m}_def.npz"), load(OUT / f"{t}_{m}_boost2.npz")
        r = a["cl_unlensed_scalar"][2:41, 0] / b["cl_unlensed_scalar"][2:41, 0] - 1
        i = int(np.argmax(np.abs(r)))
        res[f"{t}:{m}"] = {"max_abs_TT_l2_40": float(abs(r[i])), "at_l": i + 2,
                           "l=2,5,10,13,20": [float(r[x - 2]) for x in (2, 5, 10, 13, 20)]}
print(json.dumps(res, indent=1))
