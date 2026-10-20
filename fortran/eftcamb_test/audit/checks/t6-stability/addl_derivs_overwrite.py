"""compute_additional_derivs (designer generic) overwrites the analytic EFTcdotdot / EFTGamma{1,2,3}PP that designer
models (Pure EFT, RPH) compute themselves, before compute_mass_factors uses them.  Compare the cache with
EFT_mass_stability = F (model's own values kept) vs T (overwritten), and also against numeric derivatives."""
import sys, json, numpy as np
sys.path.insert(0, ".")
import t6lib, harness
import os
TREE_T = os.environ.get("T6_TREE", "A0")
models = sys.argv[1:] or ["2_PEFT_gamma1_const_1", "2_PEFT_Omega_lin_1", "2_PEFT_gamma3_lin_1", "3_RPH_KT_power_1",
                          "3_RPH_MK_power_1", "3_RPH_MB_const_1"]
g = np.logspace(-1.5, 0, 1501)
fields = ["EFTcdotdot", "EFTGamma1PP", "EFTGamma2PP", "EFTGamma3PP"]
for m in models:
    ini = harness.PARAMS_DIR / f"{m}.ini"
    rF = t6lib.probe("A0", ini, overrides={"EFTCAMB_skip_stability": "T", "feedback_level": "0", "EFT_mass_stability": "F"},
                     a_grid=g, feedback=0, fields=t6lib.DEFAULT_FIELDS + ["EFTGamma1P", "EFTGamma1PP", "EFTcdotdot"])
    rT = t6lib.probe(TREE_T, ini, overrides={"EFTCAMB_skip_stability": "T", "feedback_level": "0", "EFT_mass_stability": "T"},
                     a_grid=g, feedback=0, fields=t6lib.DEFAULT_FIELDS + ["EFTGamma1P", "EFTGamma1PP", "EFTcdotdot"])
    a, H = rF["a"], rF["adotoa"]
    num = {"EFTcdotdot": H * np.gradient(rF["EFTcdot"], np.log(a)) - 2 * H * rF["EFTcdot"],
           "EFTGamma1PP": np.gradient(rF["EFTGamma1P"], a), "EFTGamma2PP": np.gradient(rF["EFTGamma2P"], a),
           "EFTGamma3PP": np.gradient(rF["EFTGamma3P"], a)}
    sl = slice(5, -5)
    out = {}
    for f in fields:
        sc = max(np.abs(num[f][sl]).max(), 1e-300)
        out[f] = {"own_vs_numeric": float(np.abs(rF[f][sl] - num[f][sl]).max() / sc),
                  "after_addl_derivs_vs_numeric": float(np.abs(rT[f][sl] - num[f][sl]).max() / sc),
                  "max|numeric|": float(sc)}
    x = lambda r: np.minimum(r["EFT_mu1"], r["EFT_mu2"]) * (a / H) ** 2
    out["min_mu/H^2 (code path, mass on)"] = float(x(rT)[a > rT["pert_turn_on"]].min()) if (a > rT["pert_turn_on"]).any() else None
    print(m, json.dumps(out))
