"""(v) ghost_math vs physical ghost: compare the sign of the pi-equation coefficient EFTpiA1 (whose sign change the
ghost_math check detects) with the sign of EFT_kinetic (the physical no-ghost coefficient), and A1 + k^2 A2 at k_max,
on the fine-grid caches of every test model (scan_A0*/<model>.npz)."""
import json
import sys
import numpy as np
sys.path.insert(0, ".")
import t6lib
out = {}
for tag in ("scan_A0_v2", "scan_A0"):
    for f in sorted((t6lib.WORK / tag).glob("*.npz")):
        if f.stem in out:
            continue
        d = np.load(f)
        A1, K, A2 = d["EFTpiA1"], d["EFT_kinetic"], d["EFTpiA2"]
        m = (A1 != 0) & (K != 0) & np.isfinite(A1) & np.isfinite(K)
        if not m.any():
            continue
        ratio = A1[m] / K[m]
        sc = A1[m] * K[m] < 0
        out[f.stem] = {"n": int(m.sum()), "sign_mismatch": int(sc.sum()),
                       "A1/kin range": [float(ratio.min()), float(ratio.max())],
                       "A1 sign changes": int((np.diff(np.sign(A1[m])) != 0).sum()),
                       "kin sign changes": int((np.diff(np.sign(K[m])) != 0).sum())}
bad = {k: v for k, v in out.items() if v["sign_mismatch"]}
print("models:", len(out), " with sign(A1) != sign(kinetic) somewhere:", len(bad))
for k, v in out.items():
    print(f"{k:30s} {json.dumps(v)}")
