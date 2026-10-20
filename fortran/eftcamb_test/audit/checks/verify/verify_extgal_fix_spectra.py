"""Adversarial verifier (t3-extgal-gr-point-background): spectra with candidate fixes (default ini settings).

Trees: A0, hddfix (Hdotdot denominator only), allfix (hddfix + nu in ComputeAdotoa + tracker at all a).
Checks: S=0 (stability ON) -> finite and equal to EFTflag=0 (A0 1_EFT_GR) only when all three causes are fixed;
S=2 (GCCG1) and S=0.1: size of each defect on C_l, P(k), sigma8, thetastar.
"""

import json
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[5]
sys.path.insert(0, str(REPO / "fortran/eftcamb_test/audit/checks"))
import harness  # noqa: E402

W = harness.AUDIT_ROOT / "work/verify/extgal_gr"
OUT = W / "spectra"
INI = REPO / "fortran/eftcamb_test/audit/inis/008p7-ext-galileon.ini"
GR = harness.PARAMS_DIR / "1_EFT_GR.ini"
TREES = {"A0": "A0", "hddfix": W / "hddfix", "allfix": W / "allfix"}
CASES = {"S0": {"ExtendedGalileon_S": "0", "ExtendedGalileon_q": "100"},
         "S0.1": {"ExtendedGalileon_S": "0.1"},
         "S2": {}}
jobs = [(t, c) for c in CASES for t in TREES]


def one(tc):
    t, c = tc
    out = OUT / f"{t}_{c}.npz"
    st = harness.run_model(TREES[t], INI, out, overrides=CASES[c], threads=1, timeout=900)
    return (t, c), st, out


with ThreadPoolExecutor(4) as ex:
    gr_st = ex.submit(harness.run_model, "A0", GR, OUT / "A0_GR.npz", threads=1, timeout=900)
    res = {k: (st, out) for k, st, out in ex.map(one, jobs)}
    gr_st = gr_st.result()
gr = harness.load(OUT / "A0_GR.npz")
summary = {}


def desc(d):
    m = d["meta"] if isinstance(d["meta"], dict) else json.loads(str(d["meta"]))
    return m.get("derived", {})


for (t, c), (st, out) in res.items():
    line = f"{t:7s} {c:5s} status={st['status']:12s}"
    if st["status"] != "ok":
        tail = [ln.strip() for ln in st.get("log_tail", []) if "unstable" in ln or "Nan" in ln]
        print(line, tail[-2:])
        summary[f"{t}_{c}"] = {"status": st["status"], "tail": tail[-2:]}
        continue
    d = harness.load(out)
    nan = bool(np.isnan(d["cl_total"]).any())
    rec = {"status": "ok", "cl_nan": nan}
    if not nan and c == "S0":
        cmp = harness.compare(d, gr)
        rec["vs_GR_worst_cl"] = cmp["worst_cl"]
        rec["vs_GR_worst_pk"] = cmp["worst_pk"]
        rec["vs_GR_sigma8"] = cmp["sigma8"]
        rec["thetastar_vs_GR"] = desc(d)["thetastar"] / desc(gr)["thetastar"] - 1
        line += (f" vs EFTflag=0: worst_cl={cmp['worst_cl']:.2e} worst_pk={cmp['worst_pk']:.2e} "
                 f"thetastar-1={rec['thetastar_vs_GR']:+.2e}")
    print(line + f" cl_nan={nan}")
    summary[f"{t}_{c}"] = rec
for c in CASES:
    for a, b in (("A0", "hddfix"), ("hddfix", "allfix"), ("A0", "allfix")):
        sa, sb = res[(a, c)][0], res[(b, c)][0]
        if sa["status"] == "ok" and sb["status"] == "ok":
            da, db = harness.load(res[(a, c)][1]), harness.load(res[(b, c)][1])
            if np.isnan(da["cl_total"]).any() or np.isnan(db["cl_total"]).any():
                continue
            cmp = harness.compare(da, db)
            tt = cmp["cl"]["total"]
            r = {"worst_cl": cmp["worst_cl"], "worst_pk": cmp["worst_pk"], "sigma8": cmp["sigma8"],
                 "TT": tt["TT"], "EE": tt["EE"], "thetastar": desc(da)["thetastar"] / desc(db)["thetastar"] - 1}
            summary[f"{c}:{a}_vs_{b}"] = r
            print(f"{c:5s} {a} vs {b}: worst_cl={r['worst_cl']:.2e} TT={r['TT']:.2e} EE={r['EE']:.2e} "
                  f"worst_pk={r['worst_pk']:.2e} sigma8={r['sigma8']} thetastar-1={r['thetastar']:+.2e}")
(OUT / "summary.json").write_text(json.dumps(summary, indent=1, default=str))
