import sys, json, numpy as np
from pathlib import Path
sys.path.insert(0, ".")
import t6lib
D = t6lib.WORK / "scan_A0"
for j in sorted(D.glob("*.json")):
    s = json.loads(j.read_text())
    if not s.get("flags", {}).get("EFT_mass_stability"):
        continue
    d = np.load(D / (s["model"] + ".npz"))
    a = d["a"]; m = a > s["pert_turn_on"]
    x1 = d["EFT_mu1"][m] * (a[m] / d["adotoa"][m]) ** 2
    x2 = d["EFT_mu2"][m] * (a[m] / d["adotoa"][m]) ** 2
    if not m.any():
        print(f"{s['model']:28s} pert_turn_on={s['pert_turn_on']} (no points after turn-on)"); continue
    i1, i2 = np.nanargmin(x1), np.nanargmin(x2)
    print(f"{s['model']:28s} min mu1/H^2={x1[i1]: .4e} at a={a[m][i1]:.3e}  min mu2/H^2={x2[i2]: .4e} at a={a[m][i2]:.3e} nan={np.isnan(x1).sum()+np.isnan(x2).sum()}")
