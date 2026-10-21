"""(vi) Minkowski-limit tachyonic condition for full-mapping models with gamma2 = gamma3 = 0 (quintessence,
K-mouflage, ADE): tachpos = -Q1/D1 with Q1 = (EFTLambdadotdot - adotoa EFTLambdadot) a^2, D1 = ghostpos = -2 EFTc.
compute_positivity_bounds runs before compute_additional_derivs, and these models never set EFTLambdadotdot
themselves, so the code uses EFTLambdadotdot = 0.  Rebuild EFTLambdadotdot from a numerical derivative of the cached
EFTLambdadot (EFTCAMB convention: X stores Lambda a^2, so EFTLambdadotdot = d(EFTLambdadot)/dtau - 2 H EFTLambdadot)
and compare the Minkowski verdict.  Optional 2nd tree (candidate fix) is compared too."""
import json
import sys
import numpy as np
sys.path.insert(0, ".")
import t6lib, harness
trees = sys.argv[1:] or ["A0"]
models = ["5_quint_exponential", "5_quint_inverse_power", "5_quint_quadratic", "5_quint_cubic", "5_quint_axion",
          "5_Kmouflage_1", "5_ADE_1"]
g = np.logspace(-3, 0, 6001)
sl = slice(10, -10)
for m in models:
    ini = harness.PARAMS_DIR / f"{m}.ini"
    rec = {}
    for tree in trees:
        ov = {"EFTCAMB_skip_stability": "T", "feedback_level": "0", "EFT_minkowski_limit": "T"}
        r = t6lib.probe(tree, ini, overrides=ov, a_grid=g, feedback=0, cleanup_each=False)
        if r["status"] != "ok" or "a" not in r:
            rec[str(tree)[-6:]] = {"status": r["status"], "bg_error": r.get("bg_error")}
            continue
        a, H = r["a"], r["adotoa"]
        Ld = r["EFTLambdadot"]
        Ldd_num = H * np.gradient(Ld, np.log(a), edge_order=2) - 2 * H * Ld
        D1 = r["ghostpos"]
        g23 = float(np.abs(r["EFTGamma2V"]).max() + np.abs(r["EFTGamma3V"]).max())
        tach_code = r["tachpos"]
        tach_num = -(a**2) * (Ldd_num - H * Ld) / D1
        fail_code = (-D1 <= 0) | (-tach_code < 0)
        fail_num = (-D1 <= 0) | (-tach_num < 0)
        rel = np.abs(tach_code - tach_num)[sl] / np.maximum(np.abs(tach_num[sl]), 1e-300)
        rec[str(tree)[-6:]] = {
            "max|g2|+max|g3|": g23,
            "max|EFTLambdadotdot| in cache": float(np.abs(r["EFTLambdadotdot"]).max()),
            "median |Ldd_num|": float(np.median(np.abs(Ldd_num[sl]))),
            "tachpos rel err median": float(np.median(rel)), "tachpos rel err max": float(rel.max()),
            "minkowski fail pts code": int(fail_code[sl].sum()), "minkowski fail pts numeric": int(fail_num[sl].sum()),
            "n": int(len(a[sl])),
        }
    print(m, json.dumps(rec))
