"""Independent check for rgr-threshold-dimensionful-drops-de-perturbations.
1) from the A0 reference npz compute RGR functions 5..8 raw and normalised by H^2 (H=adotoa/a) to show the
   model is far from GR in dimensionless terms while raw values sit below 1e-7.
2) attribution: thr=1e-7 + skip_RGR=T must reproduce thr=1e-8 result; thr=1e-8 + turn_on_time=1.1-like (0.999)
   must reproduce thr=1e-7 result (=> effect is entirely 'perturbations off').
3) unit dependence: same model, H0 lowered -> raw functions scale with H0^2.
"""
import re, sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import numpy as np
REPO = Path("/Users/marco/IW/eftcamb_update/EFTCAMB")
sys.path.insert(0, str(REPO / "fortran/eftcamb_test/audit/checks"))
import harness
OUT = harness.AUDIT_ROOT / "work" / "verify" / "rgr_dimensionful"
M = "4_Des_mc5e_CPL_1"
ref = harness.load(harness.AUDIT_ROOT / "runs" / "A0" / f"{M}.npz")
a = ref["eft_a_grid"]; H = ref["eft_adotoa"] / a
f5 = np.abs(ref["eft_EFTc"] / a**2); f6 = np.abs(ref["eft_EFTLambda"] / a**2 + ref["eft_grhov_t"][-1]/1.0)
sel = a >= 1e-2
print("raw max f5 |c/a^2| =", f5[sel].max(), " normalised max |c/a^2|/(3H^2) =", (f5 / (3 * H**2))[sel].max())
print("w-like: max |Lambda/a^2|/(3H^2) =", (np.abs(ref["eft_EFTLambda"]/a**2) / (3*H**2))[sel].max())
print("H0^2 =", H[-1]**2)
cases = {
    "thr8": {"EFTCAMB_GR_threshold": "1e-8"},
    "thr7": {"EFTCAMB_GR_threshold": "1e-7"},
    "thr7_skip": {"EFTCAMB_GR_threshold": "1e-7", "EFTCAMB_skip_RGR": "T"},
    "thr8_late": {"EFTCAMB_GR_threshold": "1e-8", "EFTCAMB_turn_on_time": "0.9999"},
}
def job(kv):
    n, ov = kv
    out = OUT / f"{n}.npz"
    st = harness.run_model("A0", harness.PARAMS_DIR / f"{M}.ini", out, overrides=ov, threads=1, timeout=600)
    rg = re.findall(r"Return to GR time:\s*(\S+)", out.with_suffix(".log").read_text(errors="replace"))
    return n, (st["status"], rg, harness.load(out) if st["status"] == "ok" else None)
with ThreadPoolExecutor(4) as ex:
    R = dict(ex.map(job, cases.items()))
for n, (s, rg, d) in R.items():
    print(n, s, "RGR", rg, "sigma8", None if d is None else d["sigma8"][0])
def cmp(x, y):
    c = harness.compare(R[x][2], R[y][2])
    return c["worst_cl"], c["sigma8"]
print("thr7_skip vs thr8:", cmp("thr7_skip", "thr8"))
print("thr8_late vs thr7:", cmp("thr8_late", "thr7"))
print("thr7 vs thr8:", cmp("thr7", "thr8"))
