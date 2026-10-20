"""CMB impact of output-pidotdot-zero: FX (A0 + proposed fix) vs frozen A0 reference runs, via the audit harness."""

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import harness  # noqa: E402
import numpy as np  # noqa: E402

AUDIT = Path("/private/tmp/claude-501/-Users-marco-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad/audit")
FX = AUDIT / "work/t8a-equations/FX"
OUT = AUDIT / "work/verify/pidotdot/cl_FX"
PARAMS = Path(__file__).resolve().parents[3] / "parameters"

models = sys.argv[1:]
inis = [PARAMS / f"{m}.ini" for m in models]
st = harness.run_batch(FX, inis, OUT, jobs=4, threads=1)
summary = {}
for m, s in zip(models, st):
    if s["status"] != "ok":
        print(m, s["status"], s.get("log_tail", [])[-5:])
        continue
    a = harness.load(AUDIT / f"runs/A0/{m}.npz")
    b = harness.load(OUT / f"{m}.npz")
    c = harness.compare(b, a)
    ca, cb = a["cl_unlensed_scalar"], b["cl_unlensed_scalar"]
    ell = np.arange(len(ca))
    d = (cb[:, 0] - ca[:, 0]) / np.where(ca[:, 0] != 0, ca[:, 0], 1)
    sel = {L: float(d[L]) for L in (2, 3, 4, 10, 30, 100, 1000) if L < len(d)}
    imax = 2 + int(np.argmax(np.abs(d[2:2501])))
    summary[m] = {"bit_identical": c["bit_identical"], "cl_total": c["cl"].get("cl_total"),
                  "unl_TT_rel_at": sel, "unl_TT_max": float(abs(d[imax])), "unl_TT_argmax": imax,
                  "worst_pk": c.get("worst_pk"), "sigma8": c.get("sigma8")}
    print(f"{m:22s} bit_identical={c['bit_identical']} total={c['cl'].get('cl_total')} "
          f"unlTT max={abs(d[imax]):.2e}@{imax} at(2,4,30,1000)={[f'{sel.get(L, 0):+.2e}' for L in (2, 4, 30, 1000)]} "
          f"pk={c.get('worst_pk')} s8={c.get('sigma8')}")
p = OUT / "summary.json"
old = json.loads(p.read_text()) if p.exists() else {}
old.update(summary)
p.write_text(json.dumps(old, indent=1, default=str))
