"""Cross-check of the EFT_kinetic fix against the pi-equation coefficient EFTpiA1 that the perturbation code actually
evolves: for Horndeski-class models A1 = EFT_kinetic / (36 (1+Omega)^2) (ratio seen in all models with Omega' gamma2 = 0).
Prints r = 36 (1+Omega)^2 A1 / EFT_kinetic on A0 and on the fixed tree."""
import sys
import numpy as np
sys.path.insert(0, ".")
import t6lib, harness
FIX = str(t6lib.WORK / "fixkin")
g = np.logspace(-3, 0, 31)
for m in ("3_RPH_MK_OmegaDE_1", "3_RPH_MB_OmegaDE_1", "3_RPH_KT_OmegaDE_1", "3_RPH_MK_power_1", "4_Des_F_R_wCDM_2",
          "2_PEFT_Omega_exp_1", "5_Kmouflage_2"):
    ini = harness.PARAMS_DIR / f"{m}.ini"
    row = []
    for tree in ("A0", FIX):
        r = t6lib.probe(tree, ini, overrides={"EFTCAMB_skip_stability": "T", "feedback_level": "0"}, a_grid=g,
                        feedback=0)
        q = 36 * (1 + r["EFTOmegaV"]) ** 2 * r["EFTpiA1"] / r["EFT_kinetic"]
        row.append(f"[{q.min():.6f}, {q.max():.6f}]")
    print(f"{m:22s} A0 r in {row[0]:28s} fixed r in {row[1]}")
