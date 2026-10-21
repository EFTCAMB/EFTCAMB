import sys, json, numpy as np
sys.path.insert(0, ".")
import t6lib
from boundary_scans import INIS
ini = INIS / "t6_RPH_MK_hidden_ghost.ini"
g = np.linspace(0.1, 1.0, 10)
r = t6lib.probe("A0", ini, overrides={"EFTCAMB_skip_stability": "T", "feedback_level": "0"}, a_grid=g, feedback=0)
H, M, a = r["adotoa"], r["Meff2"], r["a"]
D = r["alphaK"] + 6 * r["alphaB"] ** 2
Mt = 1 + r["EFTOmegaV"] - r["EFTGamma4V"]
H0 = r["adotoa"][-1]
c, g1, g2, OmP = r["EFTc"], r["EFTGamma1V"], r["EFTGamma2V"], r["EFTOmegaP"]
kin_fixed = 9 * Mt * (4 * c * Mt + 3 * H**2 * OmP**2 * a**2 + a**2 * H0 * (H0 * (3 * g2**2 + 8 * g1 * Mt) + 6 * H * g2 * OmP))
for i in range(len(g)):
    print(f"a={a[i]:.2f} alphaK={r['alphaK'][i]: .4f} alphaM={r['alphaM'][i]: .4f} D={D[i]: .4f} kin_code={r['EFT_kinetic'][i]: .3e} kin_fixed={kin_fixed[i]: .3e} 18M^3H^2D={18*M[i]**3*H[i]**2*D[i]: .3e} grad_code={r['EFT_gradient'][i]: .3e}")
for st in ("0.1", "1e-10"):
    r2 = t6lib.probe("A0", ini, overrides={"feedback_level": "0", "EFTCAMB_stability_time": st})
    print("stability_time", st, "code verdict:", "stable" if r2["bg_error"] is None else r2["bg_error"])
