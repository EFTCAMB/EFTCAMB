import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import numpy as np
sys.path.insert(0, "/Users/marco/IW/eftcamb_update/EFTCAMB/fortran/eftcamb_test/audit/checks")
import harness
OUT = harness.AUDIT_ROOT / "work" / "verify" / "pidtypo"
FIX = str(harness.AUDIT_ROOT / "work" / "t5-robust" / "A0_pidfix")
def job(a):
    name, tree, model, ov = a
    out = OUT / f"{name}.npz"
    st = harness.run_model(tree, harness.PARAMS_DIR / f"{model}.ini", out, overrides=ov, threads=1, timeout=900)
    assert st["status"] == "ok", st
    return name, harness.load(out)
ov = {"accuracy_boost": 2, "l_accuracy_boost": 2, "transfer_high_precision": "T"}
jobs = [("GR", "A0", "1_EFT_GR", ov), ("A0", "A0", "5_Kmimic_1", ov), ("FIX", FIX, "5_Kmimic_1", ov)]
with ThreadPoolExecutor(3) as ex:
    r = dict(ex.map(job, jobs))
pk = lambda d, k: float(np.exp(np.interp(np.log(k), np.log(d["pk_k"]), np.log(d["pk"][0]))))
for n in ("A0", "FIX"):
    print(n, "boost2 P/P_GR-1:", " ".join(f"{pk(r[n],k)/pk(r['GR'],k)-1:+.4f}" for k in (1, 2, 5, 10, 20)))
