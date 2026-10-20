import sys, numpy as np
sys.path.insert(0, ".")
import t6lib, harness
from scan_all import sampler_grid
ini = harness.PARAMS_DIR / "3_RPH_MB_OmegaDE_1.ini"
ov = {"EFTCAMB_skip_stability": "T", "feedback_level": "0"}
lin, log1, log2 = sampler_grid()
for tree in sys.argv[1:]:
    r = t6lib.probe(tree, ini, overrides=ov, a_grid=log1, feedback=0, cleanup_each=False, check=True)
    print(tree, r["status"], "log1 per-point unstable:", (~r["verdict"]).sum() if "verdict" in r else None,
          "full_check:", r.get("full_check"), r["stderr"][-3:])
    # restricted full check on [2e-10, 1e-5]: the linear loop then samples 1000 points there
    r = t6lib.probe(tree, ini, overrides=ov, a_grid=np.linspace(2e-10, 1e-5, 1000), feedback=0, cleanup_each=False,
                    check=True, astart=2e-10, aend=1e-5)
    print(tree, " restricted: per-point unstable", (~r["verdict"]).sum() if "verdict" in r else None,
          "full_check:", r.get("full_check"))
