"""Nondeterminism of EFT_IC_type=3 runs on the frozen A0 release build caused by the uninitialised EV%eft_cache read in
`initial`.  Repeats the same run N times (same tree, ini, threads) and records status / sigma8 / bit-identity.
usage: python staleg_nondet.py model IC N threads [tree]"""

import json
import sys
from concurrent.futures import ThreadPoolExecutor

import numpy as np

sys.path.insert(0, "/Users/marco/IW/eftcamb_update/EFTCAMB/fortran/eftcamb_test/audit/checks")
import harness  # noqa: E402

model, ic, n, thr = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
tree = sys.argv[5] if len(sys.argv) > 5 else "A0"
env = json.loads(sys.argv[6]) if len(sys.argv) > 6 else {}
tag = sys.argv[7] if len(sys.argv) > 7 else tree
OUT = harness.AUDIT_ROOT / "work" / "verify" / "staleg" / "nondet"


def one(i):
    out = OUT / f"{model}_ic{ic}_t{thr}_{tag}_{i}.npz"
    return i, harness.run_model(tree, harness.PARAMS_DIR / f"{model}.ini", out, overrides={"EFT_IC_type": ic},
                                threads=thr, extra_env=env, timeout=1500)


with ThreadPoolExecutor(4 if thr == 1 else 2) as ex:
    sts = dict(ex.map(one, range(n)))
ref = None
res = []
for i in range(n):
    s = sts[i]["status"]
    row = {"i": i, "status": s}
    if s == "ok":
        d = harness.load(OUT / f"{model}_ic{ic}_t{thr}_{tag}_{i}.npz")
        row["sigma8"] = float(np.atleast_1d(d["sigma8"])[0])
        row["finite"] = bool(np.all(np.isfinite(d["cl_total"][2:2501])))
        if ref is None:
            ref = d
        else:
            row["bitident_vs_first_ok"] = harness.compare(d, ref)["bit_identical"]
    else:
        row["err"] = [x for x in sts[i].get("log_tail", []) if "RungeKutta" in x or "Error" in x][-1:]
        row["err"] = [e[:90] for e in row["err"]]
    res.append(row)
    print(json.dumps(row))
nok = sum(r["status"] == "ok" for r in res)
print(f"{model} IC{ic} tree={tag} threads={thr}: {nok}/{n} ok")
(OUT / f"summary_{model}_ic{ic}_t{thr}_{tag}.json").write_text(json.dumps(res, indent=1))
