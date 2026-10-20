import sys, numpy as np
sys.path.insert(0, ".")
import t6lib, harness
from scan_all import reasons
ini = harness.PARAMS_DIR / "5_quint_exponential.ini"
g = np.unique(np.concatenate([np.logspace(-10, 0, 2001)]))
ov = {"EFTCAMB_skip_stability": "T", "feedback_level": "0", "EFT_minkowski_limit": "T", "EFT_mass_stability": "T"}
r = t6lib.probe("A0", ini, overrides=ov, a_grid=g, feedback=0, cleanup_each=False, check=True)
a, H = r["a"], r["adotoa"]
print("full_check", r["full_check"], "pert_turn_on", r["pert_turn_on"])
bad = ~r["verdict"]
print("per-point unstable:", bad.sum(), "range", a[bad].min() if bad.any() else None, a[bad].max() if bad.any() else None)
rs = reasons(r)
print({k: int(v.sum()) for k, v in rs.items()})
print("ghostpos==0 points:", int((r["ghostpos"] == 0).sum()), " c==0 points:", int((r["EFTc"] == 0).sum()))
m = a > r["pert_turn_on"]
tach_fixed = -(a**2) * (r["EFTLambdadotdot"] - H * r["EFTLambdadot"]) / r["ghostpos"]
print("after turn-on: code tachpos>0 (fail) pts:", int((r["tachpos"][m] > 0).sum()), "of", int(m.sum()),
      "; tachpos rebuilt with code's own EFTLambdadotdot >0 pts:", int((tach_fixed[m] > 0).sum()))
i = np.argmin(np.abs(a - 0.5))
print(f"a=0.5: Lambdadot={r['EFTLambdadot'][i]:.4e} Lambdadotdot(addl derivs)={r['EFTLambdadotdot'][i]:.4e} "
      f"H*Lambdadot={H[i]*r['EFTLambdadot'][i]:.4e} ghostpos={r['ghostpos'][i]:.4e} tachpos code={r['tachpos'][i]:.4e} "
      f"rebuilt={tach_fixed[i]:.4e}")
