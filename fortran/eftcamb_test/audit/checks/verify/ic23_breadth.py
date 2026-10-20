"""Breadth: run the frozen A0 build with EFT_IC_type = 2 or 3 on a set of test models and compare with the A0 reference
run (EFT_IC_type = 1, default settings).  usage: python ic23_breadth.py IC model [model ...]"""

import json
import sys

import numpy as np

sys.path.insert(0, "/Users/marco/IW/eftcamb_update/EFTCAMB/fortran/eftcamb_test/audit/checks")
import harness  # noqa: E402

ic = sys.argv[1]
models = sys.argv[2:]
OUT = harness.AUDIT_ROOT / "work" / "verify" / "ic23" / "breadth"
REF = harness.AUDIT_ROOT / "runs" / "A0"
inis = [harness.PARAMS_DIR / f"{m}.ini" for m in models]
st = harness.run_batch("A0", inis, OUT, jobs=4, suffix=f"_ic{ic}", overrides={"EFT_IC_type": ic}, threads=1,
                       timeout=1200)
rows = {}
for m, s in zip(models, st):
    row = {"status": s["status"]}
    if s["status"] == "ok":
        d = harness.load(OUT / f"{m}_ic{ic}.npz")
        r = harness.load(REF / f"{m}.npz")
        row["errors"] = d["meta"].get("errors")
        if "cl_total" in d and "pk" in d:
            fin = bool(np.all(np.isfinite(d["cl_total"])) and np.all(np.isfinite(d["pk"])))
            row["finite"] = fin
            c = harness.compare(d, r)
            row.update(TT=c["cl"]["total"]["TT"], EE=c["cl"]["total"]["EE"], pk0=c["pk"][0],
                       s8=float(d["sigma8"][0]), s8_ref=float(r["sigma8"][0]))
    else:
        row["tail"] = [x[:160] for x in s.get("log_tail", [])[-2:]]
    rows[m] = row
    print(f"IC{ic} {m:28s} {json.dumps(row, default=float)}")
(OUT / f"summary_ic{ic}_{models[0]}.json").write_text(json.dumps(rows, indent=1, default=float))
