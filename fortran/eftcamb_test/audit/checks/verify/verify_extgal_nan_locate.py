"""Adversarial verifier (t3-extgal-gr-point-background): where do the NaNs come from once Hdotdot is fixed?

Runs S=0 and S=2 on hddfix / allfix with get_scale_evolution at k = 0.01, 0.1 /Mpc and reports the first a at which
each evolved field becomes non-finite, plus the eft_cache adotoa from the evolution vs get_eft_functions.
"""

import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[5]
sys.path.insert(0, str(REPO / "fortran/eftcamb_test/audit/checks"))
import harness  # noqa: E402

W = harness.AUDIT_ROOT / "work/verify/extgal_gr"
OUT = W / "nanloc"
INI = REPO / "fortran/eftcamb_test/audit/inis/008p7-ext-galileon.ini"
FAST = {"get_tensor_cls": "F", "do_lensing": "F", "l_max_scalar": "200", "get_transfer": "F",
        "num_redshiftwindows": "0"}
trees = sys.argv[1].split(",") if len(sys.argv) > 1 else ["hddfix"]
cases = {"S0": {"ExtendedGalileon_S": "0", "ExtendedGalileon_q": "100"}, "S2": {}}
for t in trees:
    for c, ov in cases.items():
        tp = "A0" if t == "A0" else W / t
        out = OUT / f"{t}_{c}.npz"
        st = harness.run_model(tp, INI, out, overrides={**FAST, **ov, "EFTCAMB_skip_stability": "T"},
                               evol_k=[0.01, 0.1], threads=1, timeout=600)
        if st["status"] != "ok":
            print(t, c, st["status"])
            continue
        d = harness.load(out)
        a = d["evol_a"][0] if d["evol_a"].ndim == 2 else d["evol_a"]  # field "a" overrides the grid
        print(f"== {t} {c}: cl NaN={bool(np.isnan(d['cl_total']).any())}")
        bad = []
        for key in sorted(k for k in d if k.startswith("evol_") and k not in ("evol_a", "evol_k")):
            v = d[key]
            if v.dtype.kind != "f":
                continue
            nf = ~np.isfinite(v)
            if nf.any():
                rows = np.argwhere(nf)
                r = rows[np.argmin(rows[:, 1])]
                bad.append((float(a[r[1]]), key, int(r[0])))
        bad.sort()
        for row in bad[:12]:
            print(f"   first non-finite at a={row[0]:.3e} k_index={row[2]}: {row[1]}")
        if not bad:
            print("   no non-finite evolved fields")
        if "evol_adotoa" in d:
            ad = d["evol_adotoa"][0]
            ae = np.interp(a, d["eft_a_grid"], d["eft_adotoa"])
            for av in (1e-4, 0.1, 0.19, 0.21, 0.5, 1.0):
                i = int(np.argmin(abs(a - av)))
                print(f"   a={a[i]:.3e}: evol adotoa / eft_functions adotoa - 1 = {ad[i] / ae[i] - 1:+.3e}")
