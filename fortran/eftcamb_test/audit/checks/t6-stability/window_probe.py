import sys, numpy as np
sys.path.insert(0, ".")
import t6lib
ini = "/Users/marco/Desktop/IW/eftcamb_update/EFTCAMB/fortran/eftcamb_test/audit/inis/t6_PEFT_gamma1_taylor_window.ini"
grid = np.linspace(0.999, 1.0, 101)
r = t6lib.probe("A0", ini, overrides={"EFTCAMB_skip_stability": "T", "feedback_level": "0"}, a_grid=grid,
                check=True, feedback=0)
print(r["status"], r.get("bg_error"), "full_check:", r.get("full_check"))
bad = grid[~r["verdict"]]
print("unstable points on fine grid:", len(bad), bad.min() if len(bad) else None, bad.max() if len(bad) else None)
print("gamma1 at a=0.9995:", r["EFTGamma1V"][50], "kinetic:", r["EFT_kinetic"][50])
r2 = t6lib.probe("A0", ini, overrides={"feedback_level": "0"})
print("get_background without skip:", r2["status"], r2.get("bg_error"))
