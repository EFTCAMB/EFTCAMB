"""Adversarial verifier (t3-extgal-gr-point-background): compare S=0 runs of candidate-fix trees (built by a parallel
verifier in work/verify/extgal_gr/{hddfix,allfix}; hddfix = Hdotdot denominator only, allfix = hddfix + massive nu
computed in ComputeAdotoa + tracker at all a + real bisection tolerance) against EFTflag=0 (A0, 1_EFT_GR.ini).
Read-only on existing npz files."""

import json
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[5]
sys.path.insert(0, str(REPO / "fortran/eftcamb_test/audit/checks"))
import harness  # noqa: E402

SP = harness.AUDIT_ROOT / "work/verify/extgal_gr/spectra"
gr = harness.load(SP / "A0_GR.npz")


def der(d):
    m = d["meta"] if isinstance(d["meta"], dict) else json.loads(str(d["meta"]))
    return m.get("derived", {})


for n in ("hddfix_S0", "allfix_S0"):
    p = SP / f"{n}.npz"
    if not p.exists():
        print(n, "missing")
        continue
    d = harness.load(p)
    nan = bool(np.isnan(d["cl_total"]).any())
    print(f"== {n}: cl_nan={nan}  keys has pk: {'pk' in d}")
    if nan:
        continue
    c = harness.compare(d, gr)
    tot = c["cl"].get("total", c["cl"])
    print("   vs EFTflag=0: worst_cl", f"{c['worst_cl']:.3e}", "worst_pk", c["worst_pk"], "sigma8", c["sigma8"])
    print("   total TT/EE/TE:", {k: f"{v:.3e}" for k, v in tot.items()} if isinstance(tot, dict) else tot)
    for k in ("thetastar", "age", "zstar", "rstar", "DAstar"):
        if k in der(d) and k in der(gr):
            print(f"   {k}: ratio-1 = {der(d)[k] / der(gr)[k] - 1:+.3e}")
    z = np.array([1.0, 3.0, 4.1, 5.0, 9.0, 30.0])
    hr = np.interp(z, d["bg_z"], d["bg_H"]) / np.interp(z, gr["bg_z"], gr["bg_H"]) - 1
    print("   H(z)/H_GR-1 at z=", z.tolist(), ":", [f"{x:+.2e}" for x in hr])
