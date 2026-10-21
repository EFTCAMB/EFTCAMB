"""Adversarial verifier (t3-extgal-gr-point-background): spectra impact of the massive-nu omission (item ii).

hdd_branch (Hdotdot fixed, tracker at all a, nu still read from the stale cache) vs allfix (nu computed inside
ComputeAdotoa), default ini settings (omnuh2 = 0.01), S = 0 and S = 2; also hdd_branch S=0 vs EFTflag=0.
"""

import json
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[5]
sys.path.insert(0, str(REPO / "fortran/eftcamb_test/audit/checks"))
import harness  # noqa: E402

W = harness.AUDIT_ROOT / "work/verify/extgal_gr"
OUT = W / "spectra"
INI = REPO / "fortran/eftcamb_test/audit/inis/008p7-ext-galileon.ini"
CASES = {"S0": {"ExtendedGalileon_S": "0", "ExtendedGalileon_q": "100"}, "S2": {}}
res = {}
for c, ov in CASES.items():
    out = OUT / f"hdd_branch_{c}.npz"
    st = harness.run_model(W / "hdd_branch", INI, out, overrides=ov, threads=2, timeout=900)
    print(c, st["status"])
    res[c] = out
gr = harness.load(OUT / "A0_GR.npz")
summ = {}
for c in CASES:
    a, b = harness.load(res[c]), harness.load(OUT / f"allfix_{c}.npz")
    cmp = harness.compare(a, b)
    da, db = a["meta"]["derived"], b["meta"]["derived"]
    r = {"worst_cl": cmp["worst_cl"], "cl_total": cmp["cl"]["total"], "worst_pk": cmp["worst_pk"],
         "sigma8": cmp["sigma8"], "bg_H": cmp["bg_H"], "thetastar": da["thetastar"] / db["thetastar"] - 1}
    summ[f"{c}:hdd_branch_vs_allfix"] = r
    print(f"{c} hdd_branch vs allfix: TT={r['cl_total']['TT']:.2e} EE={r['cl_total']['EE']:.2e} "
          f"worst_pk={r['worst_pk']:.2e} sigma8={r['sigma8']:.2e} bg_H={r['bg_H']:.2e} thetastar-1={r['thetastar']:+.2e}")
    if c == "S0":
        cg = harness.compare(a, gr)
        print(f"S0 hdd_branch vs EFTflag=0: TT={cg['cl']['total']['TT']:.2e} EE={cg['cl']['total']['EE']:.2e} "
              f"worst_pk={cg['worst_pk']:.2e}")
        summ["S0:hdd_branch_vs_GR"] = {"cl_total": cg["cl"]["total"], "worst_pk": cg["worst_pk"]}
(OUT / "nu_impact.json").write_text(json.dumps(summ, indent=1, default=str))
