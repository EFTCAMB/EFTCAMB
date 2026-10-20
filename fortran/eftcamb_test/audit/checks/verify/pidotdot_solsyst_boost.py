"""Verifier: A0 vs FX for 5_Horava_SolSyst_3 and 5_Horava_eta_4 at accuracy_boost=2."""
import sys
from concurrent.futures import ThreadPoolExecutor

import numpy as np

sys.path.insert(0, "/Users/marco/IW/eftcamb_update/EFTCAMB/fortran/eftcamb_test/audit/checks")
import harness  # noqa: E402

AUD = "/private/tmp/claude-501/-Users-marco-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad/audit/"
OUT = AUD + "work/verify/pidotdot/"
P = "/Users/marco/IW/eftcamb_update/EFTCAMB/fortran/eftcamb_test/parameters/"
trees = {"A0": AUD + "trees/A0", "FX": AUD + "work/t8a-equations/FX"}
models = ["5_Horava_SolSyst_3", "5_Horava_eta_4"]
jobs = [(t, m) for m in models for t in trees]
ov = {"accuracy_boost": 2, "l_accuracy_boost": 2}


def go(tm):
    t, m = tm
    return tm, harness.run_model(trees[t], P + m + ".ini", f"{OUT}{m}_{t}_b2.npz", overrides=ov, threads=1)


with ThreadPoolExecutor(4) as ex:
    for (t, m), st in ex.map(go, jobs):
        print(t, m, st["status"])
for m in models:
    a = harness.load(f"{OUT}{m}_A0_b2.npz")["cl_unlensed_scalar"]
    f = harness.load(f"{OUT}{m}_FX_b2.npz")["cl_unlensed_scalar"]
    d = (f[:, 0] - a[:, 0]) / f[:, 0]
    print(m, "(FX-A0)/FX TT:", {ell: f"{d[ell]:+.3e}" for ell in (2, 3, 4, 10, 30, 100, 300, 1000, 1500, 2000)})
    print(m, "EE max rel", float(np.max(np.abs((f[2:2001, 1] - a[2:2001, 1]) / f[2:2001, 1]))))
