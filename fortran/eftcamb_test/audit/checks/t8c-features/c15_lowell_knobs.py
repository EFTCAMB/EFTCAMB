"""Check 15: which accuracy knob moves the low-ell TT of 2_PEFT_gamma3_const_1 (c10/c14: default vs
accuracy_boost=2 + l_accuracy_boost=2 differs by ~0.9% at l~13, GR by <1e-3).  Ini path, l_max_scalar=300,
unlensed scalar TT l=2..40 vs the full boost2 reference.  Trees A0 and B1.

usage: python c15_lowell_knobs.py [model]
"""

from __future__ import annotations

import json
import sys
from concurrent.futures import ThreadPoolExecutor

import numpy as np
from t8c_lib import WORK, harness, load, model_ini

MODEL = sys.argv[1] if len(sys.argv) > 1 else "2_PEFT_gamma3_const_1"
OUT = WORK / "c15"
OV = {"get_tensor_cls": "F", "feedback_level": "0", "num_redshiftwindows": "0", "get_transfer": "F",
      "l_max_scalar": "300", "do_lensing": "F"}
KNOBS = {
    "def": {},
    "ref": {"accuracy_boost": "2", "l_accuracy_boost": "2"},
    "acc2": {"accuracy_boost": "2"},
    "lacc2": {"l_accuracy_boost": "2"},
    "tstep2": {"TimeStepBoost": "2"},
    "bgstep2": {"BackgroundTimeStepBoost": "2"},
    "srck2": {"SourcekAccuracyBoost": "2"},
    "intk2": {"IntkAccuracyBoost": "2"},
    "inttol4": {"IntTolBoost": "4"},
    "tswitch2": {"TimeSwitchBoost": "2"},
}


def one(job):
    tree, tag, ov = job
    return harness.run_model(tree, model_ini(MODEL), OUT / f"{tree}_{MODEL}_{tag}.npz", overrides=OV | ov,
                             threads=1)


jobs = [(t, tag, ov) for t in ("A0", "B1") for tag, ov in KNOBS.items()]
with ThreadPoolExecutor(4) as ex:
    st = list(ex.map(one, jobs))
for j, r in zip(jobs, st):
    if r["status"] != "ok":
        print(j[:2], r["status"], r.get("log_tail", [])[-3:])
res = {}
for t in ("A0", "B1"):
    ref = load(OUT / f"{t}_{MODEL}_ref.npz")["cl_unlensed_scalar"][2:41, 0]
    for tag in KNOBS:
        try:
            x = load(OUT / f"{t}_{MODEL}_{tag}.npz")["cl_unlensed_scalar"][2:41, 0]
        except FileNotFoundError:
            continue
        r = x / ref - 1
        i = int(np.argmax(np.abs(r)))
        res[f"{t}:{tag}"] = f"max|TT/ref-1| l<=40 = {abs(r[i]):.2e} at l={i + 2}; l=2: {r[0]:+.2e} l=13: {r[11]:+.2e}"
print(json.dumps(res, indent=1))
