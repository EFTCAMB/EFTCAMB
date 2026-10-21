"""(vi) Minkowski-limit / positivity coefficients for full-mapping models: EFTLambdadotdot and the third/fourth
derivatives of gamma2/gamma3 used by compute_positivity_bounds are only filled by compute_additional_derivs, which
EFTStabilityComputation calls after compute_positivity_bounds (and only with EFT_mass_stability).  Compare the code's
tachpos with the one rebuilt from a numerical derivative of the cached EFTLambdadot."""
import sys, json, numpy as np
sys.path.insert(0, ".")
import t6lib, harness
models = sys.argv[1:] or ["5_quint_exponential", "5_quint_quadratic", "5_quint_inverse_power", "5_Kmouflage_1",
                          "5_Kmouflage_2", "5_ADE_1", "5_Horava_eta_1"]
g = np.logspace(-4, 0, 4001)
res = {}
for m in models:
    ini = harness.PARAMS_DIR / f"{m}.ini"
    out = {}
    for mass in ("F", "T"):
        ov = {"EFTCAMB_skip_stability": "T", "feedback_level": "0", "EFT_minkowski_limit": "T",
              "EFT_positivity_bounds": "T", "EFT_mass_stability": mass}
        r = t6lib.probe("A0", ini, overrides=ov, a_grid=g, feedback=0, cleanup_each=False, check=True)
        if r["status"] != "ok" or "a" not in r:
            out[mass] = {"status": r["status"], "err": r.get("bg_error"), "stderr": r["stderr"][-3:]}
            continue
        a, H = r["a"], r["adotoa"]
        Ldd_num = H * np.gradient(r["EFTLambdadot"], np.log(a))  # d(Lambdadot)/dtau
        D1 = r["ghostpos"]
        gam23 = np.abs(r["EFTGamma2V"]).max() + np.abs(r["EFTGamma3V"]).max()
        tach_num = -(a**2) * (Ldd_num - H * r["EFTLambdadot"]) / D1
        sl = slice(5, -5)
        out[mass] = {
            "status": r["status"], "full_check_with_minkowski": r.get("full_check"),
            "max|gamma2|+max|gamma3|": float(gam23),
            "max|EFTLambdadotdot|_code": float(np.abs(r["EFTLambdadotdot"]).max()),
            "max|Ldd_numeric|": float(np.abs(Ldd_num[sl]).max()),
            "code_minkowski_fail_pts": int(((-r["ghostpos"] <= 0) | (-r["tachpos"] < 0)).sum()),
            "numeric_minkowski_fail_pts": int(((-D1[sl] <= 0) | (-tach_num[sl] < 0)).sum()),
            "tachpos_code_range": [float(r["tachpos"].min()), float(r["tachpos"].max())],
            "tachpos_numeric_range": [float(tach_num[sl].min()), float(tach_num[sl].max())],
        }
    res[m] = out
    print(m, json.dumps(out))
