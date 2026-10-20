import sys, json, numpy as np
sys.path.insert(0, ".")
import t6lib
from boundary_scans import INIS
from scan_all import reasons
ini = INIS / "t6_RPH_MK_hidden_ghost.ini"
g = np.logspace(-10, 0, 3001)
for k0 in ("-0.2", "-0.1", "-0.05"):
    ov = {"EFTCAMB_skip_stability": "T", "feedback_level": "0", "EFTCAMB_stability_time": "1e-10", "RPHkineticity_ODE0": k0}
    r = t6lib.probe("A0", ini, overrides=ov, a_grid=g, feedback=0, check=True, cleanup_each=False)
    rs = reasons(r)
    bad = ~r["verdict"]
    D = r["alphaK"]
    print(k0, "full_check", r["full_check"], {k: (int(v.sum()), float(g[v].min()), float(g[v].max())) for k, v in rs.items() if v.any()},
          "D<0 pts:", int((D < 0).sum()), "of", len(g))
    for st in ("0.1", "1e-3", "1e-6", "1e-10"):
        r2 = t6lib.probe("A0", ini, overrides={"feedback_level": "0", "EFTCAMB_stability_time": st, "RPHkineticity_ODE0": k0})
        print("    stability_time", st, "->", "stable" if r2["bg_error"] is None else "unstable")
