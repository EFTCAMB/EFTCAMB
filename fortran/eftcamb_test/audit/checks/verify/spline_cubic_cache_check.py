"""Independent check of spline-cubic-endpoints-derivatives: inspect the EFT timestep cache of the Pure EFT run with
Omega given as a cubic spline of the exact line Omega=a (nodes 0..1) vs the analytic linear Omega=a run."""
import sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import harness

d = Path(sys.argv[1])
s, r = harness.load(d / "spline.npz"), harness.load(d / "ref.npz")
a = s["eft_a_grid"]
keys = sorted(k for k in s if k.startswith("eft_") and k != "eft_a_grid")
for k in ["eft_EFTOmegaV", "eft_EFTOmegaP", "eft_EFTOmegaPP", "eft_EFTOmegaPPP"]:
    if k not in s:
        print("missing", k, [x for x in keys if "Omega" in x][:20]); continue
    for aa in [1e-6, 1e-3, 0.01, 0.05, 0.1, 0.3, 0.5, 0.7, 0.9, 0.95, 0.99, 1.0]:
        i = np.argmin(abs(a - aa))
        print(f"{k:18s} a={a[i]:.3e} spline={s[k][i]: .6e} ref={r[k][i]: .6e}")
