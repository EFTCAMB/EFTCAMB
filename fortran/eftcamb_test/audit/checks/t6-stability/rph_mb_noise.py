import sys, numpy as np
sys.path.insert(0, ".")
import t6lib, harness
ini = harness.PARAMS_DIR / "3_RPH_MB_OmegaDE_1.ini"
g = np.logspace(-10, -4, 25)
ov = {"EFTCAMB_skip_stability": "T", "feedback_level": "0"}
r1 = t6lib.probe("A0", ini, overrides=ov, a_grid=g, feedback=0, check=True)
r2 = t6lib.probe("A0", ini, overrides=ov, a_grid=g[::-1], feedback=0)
k1, k2 = r1["EFT_kinetic"], r2["EFT_kinetic"][::-1]
print("full_check", r1["full_check"])
print("same kinetic in both orders:", np.array_equal(k1, k2))
H = r1["adotoa"]; M = r1["Meff2"]
D = r1["alphaK"] + 6 * r1["alphaB"] ** 2
for i in range(0, 25, 3):
    print(f"a={g[i]:.2e} kin={k1[i]: .3e} analytic 18M^3H^2D={18*M[i]**3*H[i]**2*D[i]: .3e} c={r1['EFTc'][i]: .3e} alphaB={r1['alphaB'][i]: .3e} verdict={r1['verdict'][i]}")
