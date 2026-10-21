"""RPH with alpha_i proportional to Omega_DE: EFT_kinetic = 18 M^3 H^2 (alpha_K + 6 alpha_B^2) is computed through
the cancellation 2 M c - 2 M c (gamma1 is built from c), so at early times its sign is roundoff.  Run the test models
with the default EFTCAMB_stability_time = 1e-10 instead of the 0.1 their inis set."""
import sys, json, numpy as np
sys.path.insert(0, ".")
import t6lib, harness
TREE = sys.argv[1] if len(sys.argv) > 1 else "A0"
for m in ("3_RPH_MB_OmegaDE_1", "3_RPH_MK_OmegaDE_1", "3_RPH_KT_OmegaDE_1"):
    ini = harness.PARAMS_DIR / f"{m}.ini"
    out = {}
    for st in ("0.1", "1e-10", "1e-4", "1e-3"):
        r = t6lib.probe(TREE, ini, overrides={"feedback_level": "0", "EFTCAMB_stability_time": st})
        out[st] = "stable" if r["bg_error"] is None else "unstable"
    g = np.logspace(-10, 0, 2001)
    r = t6lib.probe(TREE, ini, overrides={"EFTCAMB_skip_stability": "T", "feedback_level": "0"}, a_grid=g, feedback=0)
    H, M = r["adotoa"], r["Meff2"]
    D = r["alphaK"] + 6 * r["alphaB"] ** 2
    exact = 18 * M**3 * H**2 * D
    neg = r["EFT_kinetic"] < 0
    rel = np.abs(r["EFT_kinetic"] - exact) / np.maximum(np.abs(exact), 1e-300)
    out["n_kinetic_neg_of_2001"] = int(neg.sum())
    out["n_analytic_D_neg"] = int((D < 0).sum())
    out["last_a_kinetic_neg"] = float(g[neg].max()) if neg.any() else None
    out["max_rel_err_kinetic_a>1e-2"] = float(rel[g > 1e-2].max())
    out["max_rel_err_kinetic_a<1e-6"] = float(rel[g < 1e-6].max())
    print(m, json.dumps(out))
