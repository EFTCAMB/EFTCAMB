import sys, numpy as np
sys.path.insert(0, ".")
import t6lib
ini = "/Users/marco/Desktop/IW/eftcamb_update/EFTCAMB/fortran/eftcamb_test/audit/inis/t6_RPH_KB_const.ini"
g = np.logspace(-10, 0, 41)
for K in (0.1, -0.05, -0.07):
    r = t6lib.probe("A0", ini, overrides={"EFTCAMB_skip_stability": "T", "feedback_level": "0", "RPHkineticity0": K},
                    a_grid=g, check=True, feedback=0)
    print(K, r["status"], r.get("bg_error"), r.get("full_check"))
    if "verdict" in r:
        D = r["alphaK"] + 6 * r["alphaB"] ** 2
        print("  D range", D.min(), D.max(), " kin<0:", (r["EFT_kinetic"] < 0).sum(), " grad<0:", (r["EFT_gradient"] < 0).sum(),
              " unstable:", (~r["verdict"]).sum(), " Meff2", r["Meff2"][[0, -1]], " AT", r["EFTAT"][[0, -1]])
        M = r["Meff2"]; H = r["adotoa"]
        print("  kin/(18 M^3 H^2 D) =", (r["EFT_kinetic"] / (18 * M**3 * H**2 * D))[[0, 20, 40]])
