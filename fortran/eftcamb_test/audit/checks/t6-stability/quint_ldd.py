"""Consistency of the cached time derivatives: the cache stores X a^2 for X = c, Lambda and their conformal-time
derivatives, so dEFTX/dtau = EFTXdot + 2 H EFTX."""
import sys, numpy as np
sys.path.insert(0, ".")
import t6lib, harness
name = sys.argv[1] if len(sys.argv) > 1 else "5_quint_exponential"
ini = harness.PARAMS_DIR / f"{name}.ini"
g = np.logspace(-3, 0, 3001)
ov = {"EFTCAMB_skip_stability": "T", "feedback_level": "0", "EFT_mass_stability": "T"}
r = t6lib.probe("A0", ini, overrides=ov, a_grid=g, feedback=0, cleanup_each=False)
a, H = r["a"], r["adotoa"]
d = lambda f: H * np.gradient(r[f], np.log(a))
Ld_num = d("EFTLambda") - 2 * H * r["EFTLambda"]
Ldd_num = d("EFTLambdadot") - 2 * H * r["EFTLambdadot"]
cd_num = d("EFTc") - 2 * H * r["EFTc"]
print(name)
for x in (0.01, 0.1, 0.3, 0.5, 0.8, 0.99):
    i = np.argmin(np.abs(a - x))
    print(f"a={a[i]:.3f} Ldot={r['EFTLambdadot'][i]: .4e} num={Ld_num[i]: .4e} | Lddot={r['EFTLambdadotdot'][i]: .4e} num={Ldd_num[i]: .4e} | cdot={r['EFTcdot'][i]: .4e} num={cd_num[i]: .4e}")
