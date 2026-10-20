import sys, numpy as np
sys.path.insert(0, ".")
import t6lib, harness
m = sys.argv[1]
ini = harness.PARAMS_DIR / f"{m}.ini"
g = np.logspace(-2, 0, 5)
r = t6lib.probe("A0", ini, overrides={"EFTCAMB_skip_stability": "T", "feedback_level": "0", **dict(s.split("=") for s in sys.argv[2:])}, a_grid=g, feedback=0)
H, M, a = r["adotoa"], r["Meff2"], r["a"]
H0 = 70 / 299792.458
Om, OmP, c, g1, g2, g4 = r["EFTOmegaV"], r["EFTOmegaP"], r["EFTc"], r["EFTGamma1V"], r["EFTGamma2V"], r["EFTGamma4V"]
Mt = 1 + Om - g4
kin_py = 9 * Mt * (4 * c * Mt + 3 * H**2 * OmP**2 * a**2 + a**2 * H0 * (H0 * (3 * g2**2 + 8 * g1 * Mt) + 6 * H * g2 * OmP))
print("code formula re-evaluated / code:", kin_py / r["EFT_kinetic"])
print("4 g1 H0^2 a^2 vs aK M H^2 - 2c:", (4 * g1 * H0**2 * a**2) / (r["alphaK"] * M * H**2 - 2 * c))
print("H Om' + H0 g2 vs 2 aB M H / a:", (H * OmP + H0 * g2), 2 * r["alphaB"] * M * H / a)
print("terms: 4cM:", 4 * c * Mt, " 3a^2(HOm'+H0g2)^2:", 3 * a**2 * (H * OmP + H0 * g2) ** 2, " 8g1 M a^2 H0^2:", 8 * g1 * Mt * a**2 * H0**2)
