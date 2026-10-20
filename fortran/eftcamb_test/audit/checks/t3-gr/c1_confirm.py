"""t3-gr (resume): cheap reruns confirming the key part-(a) numbers, and the 3+1 designer GR point with the
`success` fix (R7 = A0 + p31_success).

(a) A0 vs B0 (release), and R1nofma (A0 + maxz50 + polter_orig + weyl_orig, -ffp-contract=off) vs B0nofma, base GR ini.
(b) 3p1 GR point on A0 and R7, compared with A0 EFTflag=0 (IntTolBoost=2).
"""

from __future__ import annotations

import os
from concurrent.futures import ThreadPoolExecutor

from common import GR_INI, INIS, WORK, fmt, harness, run

SUB = "c1"
P31 = INIS / "007p8-3p1.ini"
# GR: sigma1 = rho1 = 1 (constant), no ODE part, rho1N = rho2 = 0, w_DE = -1
P31_GR = {"EFTwDE": "0", "Designer3p1modelsigma1": "1", "Designer3p1sigma10": "1.0",
          "Designer3p1modelsigma1_ODE": "0", "Designer3p1modelrho1N": "0", "Designer3p1modelrho2": "0"}
JOBS = {
    "A0_base": ("A0", GR_INI, {}, {}),
    "B0_base": ("B0", GR_INI, {}, {}),
    "R1nofma_base": (WORK / "R1nofma", GR_INI, {}, {}),
    "B0nofma_base": (WORK / "B0nofma", GR_INI, {}, {}),
    "A0g_base": ("A0", GR_INI, {}, {"Accuracy.IntTolBoost": 2.0}),
    "A0_3p1GR": ("A0", P31, P31_GR, {}),
    "R7_3p1GR": (WORK / "R7", P31, P31_GR, {}),
    "R7_3p1GRforcedNS": (WORK / "R7", P31, {**P31_GR, "EFTCAMB_GR_threshold": "-1", "EFTCAMB_skip_stability": "T"},
                         {}),
    "R7_3p1": (WORK / "R7", P31, {}, {}),
    "A0_3p1GRb": ("A0", P31, P31_GR, {}),
    "R7_3p1GRb": (WORK / "R7", P31, P31_GR, {}),
    "R7_3p1GRbforcedNS": (WORK / "R7", P31, {**P31_GR, "EFTCAMB_GR_threshold": "-1",
                                             "EFTCAMB_skip_stability": "T"}, {}),
    "R7_3p1GRbforced": (WORK / "R7", P31, {**P31_GR, "EFTCAMB_GR_threshold": "-1"}, {}),
}


def one(item):
    name, (tree, ini, ov, ps) = item
    out = WORK / "runs" / SUB / f"{name}.npz"
    if out.exists():
        return name, out
    return name, run(tree, name, ini=ini, overrides=ov, python_set=ps, sub=SUB, threads=1, timeout=1200)


with ThreadPoolExecutor(int(os.environ.get("T3_JOBS", "3"))) as ex:
    paths = dict(ex.map(one, JOBS.items()))

def cmp(a, b):
    pa, pb = paths[a], paths[b]
    if not (pa.exists() and pb.exists()):
        print(f"{a} vs {b}: missing")
        return
    print(f"{a:18s} vs {b:14s}", fmt(harness.compare(pa, pb)))


cmp("A0_base", "B0_base")
cmp("R1nofma_base", "B0nofma_base")
for n in ("A0_3p1GR", "R7_3p1GR", "R7_3p1GRforcedNS", "A0_3p1GRb", "R7_3p1GRb", "R7_3p1GRbforcedNS",
          "R7_3p1GRbforced"):
    cmp(n, "A0g_base")
