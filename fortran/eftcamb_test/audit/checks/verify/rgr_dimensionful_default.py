"""Default threshold 1e-8, CPL w0=-0.97, wa=0 (designer mc5e): RGR vs skip_RGR."""
import re, sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import numpy as np
sys.path.insert(0, "/Users/marco/IW/eftcamb_update/EFTCAMB/fortran/eftcamb_test/audit/checks")
import harness
OUT = harness.AUDIT_ROOT / "work" / "verify" / "rgr_dimensionful"
M = "4_Des_mc5e_CPL_1"
def job(kv):
    n, ov = kv
    out = OUT / f"{n}.npz"
    st = harness.run_model("A0", harness.PARAMS_DIR / f"{M}.ini", out, overrides=ov, threads=1, timeout=600)
    rg = re.findall(r"Return to GR time:\s*(\S+)", out.with_suffix(".log").read_text(errors="replace"))
    return n, (st["status"], rg, harness.load(out))
base = {"EFTw0": "-0.97", "EFTwa": "0.0"}
with ThreadPoolExecutor(2) as ex:
    R = dict(ex.map(job, {"w097_def": base, "w097_skip": {**base, "EFTCAMB_skip_RGR": "T"}}.items()))
for n, (s, rg, d) in R.items():
    a = d["eft_a_grid"]; sel = a > 1e-2
    print(n, s, "RGR", rg[:1], "sigma8", d["sigma8"][0], "max|c/a^2|", np.abs(d["eft_EFTc"]/a**2)[sel].max())
c = harness.compare(R["w097_def"][2], R["w097_skip"][2])
print("default vs skip: worst_cl", c["worst_cl"], "sigma8", c["sigma8"])
