"""t3-gr (b): fine B0 scan of designer f(R) across the RGR switch (default EFTCAMB_GR_threshold = 1e-8), matter power
only. The physical P(k) is a continuous function of B0; a jump to exactly LCDM at the B0 where the RGR time crosses
a = 1 is an artefact of the return-to-GR approximation.
"""

import re
from concurrent.futures import ThreadPoolExecutor

import numpy as np

from common import PARAMS, RUNS, harness

INI = PARAMS / "4_Des_F_R_LCDM_1.ini"
SUB = RUNS / "b6"
PK = {"get_scalar_cls": "F", "get_tensor_cls": "F", "do_lensing": "F", "transfer_kmax": "10",
      "num_redshiftwindows": "0"}
VALS = np.geomspace(1.2e-7, 2e-8, 14)


def one(v):
    out = SUB / f"pkfine_{v:.3e}.npz"
    if not out.with_suffix(".json").exists():
        harness.run_model("A0", INI, out, overrides={**PK, "EFTB0": f"{v:.6e}"}, threads=1, timeout=300)
    return v, out


with ThreadPoolExecutor(4) as ex:
    res = list(ex.map(one, VALS))
ref = harness.load(SUB / "pk_GR.npz")
k = ref["pk_k"]
ik = [int(np.argmin(abs(k - x))) for x in (0.3, 1.0, 3.0, 7.0)]
print(f"{'B0':>10s} {'RGR a':>8s} " + " ".join(f"{'P(k=%.2g)-1' % k[i]:>12s}" for i in ik))
for v, out in res:
    st = harness.json.loads(out.with_suffix(".json").read_text())["status"]
    if st != "ok":
        print(f"{v:10.3e} {st}")
        continue
    d = harness.load(out)
    rgr = re.findall(r"Return to GR time: *([-+0-9.E]+)", out.with_suffix(".log").read_text())
    r = d["pk"][0] / ref["pk"][0] - 1
    print(f"{v:10.3e} {float(rgr[-1]):8.4f} " + " ".join(f"{r[i]:+12.3e}" for i in ik))
