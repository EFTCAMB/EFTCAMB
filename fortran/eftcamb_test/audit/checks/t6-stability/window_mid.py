"""(ii) sampling resolution: ghost window |a - a0| < 2e-4 placed between the sampler's grid points at intermediate a
(a0 = 0.2998: no linear, log-from-astart or (fixed) log-from-aend point falls inside).  Checks A0 and the stale-y fix;
also with EFT_ghost_math_stability = T (sign change of A1 between samples)."""
import sys
import numpy as np
sys.path.insert(0, ".")
import t6lib
from scan_all import sampler_grid
ini = t6lib.HERE.parents[1] / "inis" / "t6_PEFT_gamma1_taylor_window.ini"
FIXY = str(t6lib.WORK / "fixy")
a0, hw = 0.2998, 2e-4
lin, log1, log3 = sampler_grid(1e-10, 1.0, fixed=True)
allg = np.concatenate([lin, log1, log3])
print("sampler points inside window:", int(((allg > a0 - hw) & (allg < a0 + hw)).sum()))
ov = {"EFTGamma1a0": repr(a0), "EFTGamma10": "-1.e-3", "EFTGamma12": "2.5e4", "feedback_level": "0"}
g = np.linspace(a0 - 3e-4, a0 + 3e-4, 61)
r = t6lib.probe("A0", ini, overrides={**ov, "EFTCAMB_skip_stability": "T"}, a_grid=g, feedback=0)
bad = g[~r["verdict"]]
print(f"per-point unstable: {len(bad)} in [{bad.min():.6f}, {bad.max():.6f}]  kinetic(a0) = {r['EFT_kinetic'][30]:.3e}")
for tree in ("A0", FIXY):
    for gm in ("F", "T"):
        rr = t6lib.probe(tree, ini, overrides={**ov, "EFT_ghost_math_stability": gm})
        print(f"{tree[-4:]:>4s} ghost_math={gm}: verdict = {'stable' if rr['bg_error'] is None else 'unstable'}")
