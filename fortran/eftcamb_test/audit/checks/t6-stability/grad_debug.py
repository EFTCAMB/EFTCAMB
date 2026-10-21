import sys, numpy as np
sys.path.insert(0, ".")
import t6lib
from boundary_scans import INIS, gradient_numerator
ini = INIS / "t6_RPH_KB_const.ini"
g = np.logspace(-10, 0, 41)
for B in (-0.15, -0.2):
    r = t6lib.probe("A0", ini, overrides={"EFTCAMB_skip_stability": "T", "feedback_level": "0", "RPHkineticity0": "0.1",
                                         "RPHbraiding0": repr(B)}, a_grid=g, feedback=0)
    a, H = r["a"], r["adotoa"]
    h = (r["Hdot"] - H**2) / H**2
    q = (r["grhom_t"] + r["gpresm_t"]) / (2 * H**2 * r["Meff2"])
    nplus = gradient_numerator(r, B)
    nminus = gradient_numerator(r, -B)
    print("B", B, "alphaB cache", r["alphaB"][0], "Meff2", r["Meff2"][0], "alphaM", r["alphaM"][0])
    for i in range(0, 41, 4):
        print(f"  a={a[i]:.1e} h={h[i]: .4f} q={q[i]: .4f} N(+B)={nplus[i]: .4e} N(-B)={nminus[i]: .4e} code_grad={r['EFT_gradient'][i]: .4e}")
