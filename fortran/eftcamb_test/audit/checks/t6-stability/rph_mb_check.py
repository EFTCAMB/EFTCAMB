import sys, numpy as np
sys.path.insert(0, ".")
import t6lib, harness
from scan_all import sampler_grid
ini = harness.PARAMS_DIR / "3_RPH_MB_OmegaDE_1.ini"
ov = {"EFTCAMB_skip_stability": "T", "feedback_level": "0"}
lin, log1, log2 = sampler_grid()
# exact sampler order: lin, then log1, cleanup in between as the code does
r = t6lib.probe("A0", ini, overrides=ov, a_grid=np.concatenate([lin]), feedback=0, cleanup_each=False, check=True)
print("lin: unstable", (~r["verdict"]).sum(), "full_check", r["full_check"])
r = t6lib.probe("A0", ini, overrides=ov, a_grid=log1, feedback=0, cleanup_each=False)
print("log1: unstable", (~r["verdict"]).sum(), "first unstable a", log1[~r["verdict"]][:3])
r = t6lib.probe("A0", ini, overrides=ov, a_grid=log1, feedback=3, cleanup_each=False, check=True)
print("log1 with feedback 3: unstable", (~r["verdict"]).sum(), "full", r["full_check"])
print([l for l in r["log"] if "nstab" in l or "stable" in l][-6:])
