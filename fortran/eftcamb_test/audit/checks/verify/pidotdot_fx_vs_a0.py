"""Verifier: A0 vs FX (pi'' fix) on Gamma4/Gamma5/Omega/Horava models via harness."""
import json
import sys
from concurrent.futures import ThreadPoolExecutor

sys.path.insert(0, "/Users/marco/IW/eftcamb_update/EFTCAMB/fortran/eftcamb_test/audit/checks")
import harness  # noqa: E402

AUD = "/private/tmp/claude-501/-Users-marco-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad/audit/"
OUT = AUD + "work/verify/pidotdot/"
P = "/Users/marco/IW/eftcamb_update/EFTCAMB/fortran/eftcamb_test/parameters/"
models = sys.argv[1:]
trees = {"A0": AUD + "trees/A0", "FX": AUD + "work/t8a-equations/FX"}
jobs = [(t, m) for m in models for t in trees]


def go(tm):
    t, m = tm
    return tm, harness.run_model(trees[t], P + m + ".ini", f"{OUT}{m}_{t}.npz", threads=1)


with ThreadPoolExecutor(4) as ex:
    for (t, m), st in ex.map(go, jobs):
        print(t, m, st["status"])
for m in models:
    c = harness.compare(harness.load(f"{OUT}{m}_A0.npz"), harness.load(f"{OUT}{m}_FX.npz"), lmax=2000)
    print(m, json.dumps({"cl_total": c["cl"].get("cl_total", c["cl"]), "worst_cl": c["worst_cl"],
                         "worst_pk": c["worst_pk"], "sigma8": c["sigma8"], "bit": c["bit_identical"]}, default=str)[:900])
