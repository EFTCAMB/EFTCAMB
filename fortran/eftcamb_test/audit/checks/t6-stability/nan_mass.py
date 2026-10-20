import sys, json, numpy as np
sys.path.insert(0, ".")
import t6lib, harness
for m in ("2_PEFT_gamma3_lin_1", "2_PEFT_gamma3_const_1", "2_PEFT_gamma3_power_1", "2_PEFT_gamma5_const_1", "5_quint_exponential"):
    ini = harness.PARAMS_DIR / f"{m}.ini"
    g = np.logspace(-2, 0, 201)
    r = t6lib.probe("A0", ini, overrides={"EFTCAMB_skip_stability": "T", "feedback_level": "0", "EFT_mass_stability": "T"},
                    a_grid=g, feedback=0)
    nan = ~np.isfinite(r["EFT_mu1"]) | ~np.isfinite(r["EFT_mu2"])
    r2 = t6lib.probe("A0", ini, overrides={"feedback_level": "0", "EFT_mass_stability": "T"})
    others = {f: int((~np.isfinite(r[f])).sum()) for f in t6lib.DEFAULT_FIELDS if f not in ("EFT_mu1", "EFT_mu2") and (~np.isfinite(r[f])).any()}
    print(m, "pert_turn_on", r["pert_turn_on"], "NaN mu points:", int(nan.sum()), "/", len(g), "other NaN fields:", others,
          "| per-point verdict on NaN points:", set(r["verdict"][nan].tolist()), "| code verdict:", "stable" if r2["bg_error"] is None else r2["bg_error"])
