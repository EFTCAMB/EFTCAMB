"""Verifier: designer f(R) small B0, P(k) vs EFTCAMB_GR_threshold (independent of the repro)."""
import re
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[5]
sys.path.insert(0, str(REPO / "fortran" / "eftcamb_test" / "audit" / "checks"))
import harness  # noqa: E402

INI = harness.PARAMS_DIR / "4_Des_F_R_LCDM_1.ini"
GR = harness.PARAMS_DIR / "1_EFT_GR.ini"
OUT = harness.AUDIT_ROOT / "work" / "verify" / "t3_fr_rgr"
PK = {"get_scalar_cls": "F", "get_tensor_cls": "F", "do_lensing": "F", "transfer_kmax": "10",
      "num_redshiftwindows": "0"}
B0S = [float(x) for x in sys.argv[1].split(",")]
THR = [float(x) for x in sys.argv[2].split(",")]
TIMEOUT = int(sys.argv[3]) if len(sys.argv) > 3 else 400


def job(a):
    name, ini, ov = a
    out = OUT / f"{name}.npz"
    t = time.time()
    st = harness.run_model("A0", ini, out, overrides={**PK, **ov}, python_set={"Accuracy.IntTolBoost": 2.0} if name == "GR" else {},
                           threads=1, timeout=TIMEOUT)
    return name, st, out, time.time() - t


jobs = [("GR", GR, {})]
for b in B0S:
    for th in THR:
        jobs.append((f"fR_{b:.2e}_thr{th:.0e}", INI, {"EFTB0": repr(b), "EFTCAMB_GR_threshold": repr(th)}))
with ThreadPoolExecutor(3) as ex:
    res = {n: (st, out, dt) for n, st, out, dt in ex.map(job, jobs)}
ref = harness.load(res["GR"][1])
k = ref["pk_k"]
ik = [int(np.argmin(abs(k - x))) for x in (0.3, 1.0, 3.0, 7.0)]
for b in B0S:
    for th in THR:
        st, out, dt = res[f"fR_{b:.2e}_thr{th:.0e}"]
        if st["status"] != "ok":
            print(f"B0={b:.2e} thr={th:.0e} {st['status']} t={dt:.0f}s")
            continue
        rgr = re.findall(r"Return to GR time: *([-+0-9.E]+)", out.with_suffix(".log").read_text())
        r = harness.load(out)["pk"][0] / ref["pk"][0] - 1
        print(f"B0={b:.2e} thr={th:.0e} RGR={float(rgr[-1]) if rgr else float('nan'):.4g} t={dt:.0f}s "
              + " ".join(f"k={k[i]:.2g}:{r[i]:+.3e}" for i in ik))
