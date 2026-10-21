"""(iii) stability check interval [EFTCAMB_stability_time=1e-10, 1] vs designer background table range
[model_background_a_ini (default 1e-8), 1]: below a_ini the tables return their first value (Hdot frozen)."""
import sys, json, numpy as np
sys.path.insert(0, ".")
import t6lib
from boundary_scans import INIS
ini = INIS / "t6_RPH_KB_const.ini"
g = np.array([1e-10, 1e-9, 5e-9, 2e-8, 1e-7])
base = {"EFTCAMB_skip_stability": "T", "feedback_level": "0", "RPHkineticity0": "0.1", "RPHbraiding0": "-0.15"}
r1 = t6lib.probe("A0", ini, overrides=base, a_grid=g, feedback=0)
r2 = t6lib.probe("A0", ini, overrides={**base, "model_background_a_ini": "1e-12"}, a_grid=g, feedback=0)
for i, a in enumerate(g):
    h1 = (r1["Hdot"][i] - r1["adotoa"][i]**2) / r1["adotoa"][i]**2
    h2 = (r2["Hdot"][i] - r2["adotoa"][i]**2) / r2["adotoa"][i]**2
    print(f"a={a:.0e} dH/dt/H^2: default={h1: .4f} a_ini=1e-12: {h2: .4f} | EFT_gradient default={r1['EFT_gradient'][i]: .4e} a_ini=1e-12: {r2['EFT_gradient'][i]: .4e} | EFTc {r1['EFTc'][i]: .3e} {r2['EFTc'][i]: .3e}")
# verdict comparison over alpha_B
res = []
for B in np.round(np.arange(-0.20, -0.139, 0.0025), 5):
    v = []
    for ainit in (None, "1e-12"):
        ov = {"feedback_level": "0", "RPHkineticity0": "0.1", "RPHbraiding0": repr(float(B))}
        if ainit:
            ov["model_background_a_ini"] = ainit
        r = t6lib.probe("A0", ini, overrides=ov)
        v.append(r["bg_error"] is None)
    res.append((float(B), v[0], v[1]))
print("B, stable(default a_ini=1e-8), stable(a_ini=1e-12):")
for x in res:
    print("  ", x, "" if x[1] == x[2] else "<-- differs")
