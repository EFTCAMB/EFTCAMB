import sys, numpy as np
sys.path.insert(0, ".")
import t6lib, harness
ini = harness.PARAMS_DIR / "2_PEFT_gamma1_const_1.ini"
r = t6lib.probe("A0", ini, overrides={"EFTCAMB_skip_stability": "T", "feedback_level": "0"},
                a_grid=np.array([1e-10, 1e-3, 0.02, 0.5, 1.0]), check=True, feedback=0)
for k in ("status", "bg_error", "full_check", "pert_turn_on", "EFTflag"):
    print(k, r.get(k))
for f in ("a", "adotoa", "EFTGamma1V", "EFTAT", "EFTDT", "EFT_kinetic", "EFT_gradient", "EFT_mu1", "verdict"):
    print(f, r.get(f))
print(r["stderr"][-5:], r["log"][-5:])
