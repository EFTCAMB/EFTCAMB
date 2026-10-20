import sys, numpy as np
sys.path.insert(0, ".")
import t6lib, harness
m = sys.argv[1]
ini = harness.PARAMS_DIR / f"{m}.ini"
g = np.logspace(-2, 0, 9)
r = t6lib.probe("A0", ini, overrides={"EFTCAMB_skip_stability": "T", "feedback_level": "0"}, a_grid=g, feedback=0)
H, M = r["adotoa"], r["Meff2"]
D = r["alphaK"] + 6 * r["alphaB"] ** 2
exact = 18 * M**3 * H**2 * D
Mt = 1 + r["EFTOmegaV"] - r["EFTGamma4V"]
for i, a in enumerate(g):
    print(f"a={a:.3f} kin={r['EFT_kinetic'][i]: .4e} 18M^3H^2D={exact[i]: .4e} ratio={r['EFT_kinetic'][i]/exact[i]: .4f} "
          f"Meff2={M[i]:.6f} 1+Om-g4={Mt[i]:.6f} aK={r['alphaK'][i]: .3e} aB={r['alphaB'][i]: .3e}")
