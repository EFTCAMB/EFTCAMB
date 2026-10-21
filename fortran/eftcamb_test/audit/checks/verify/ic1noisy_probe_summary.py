"""Summarise ic1noisy_probe outputs (work/verify/ic1noisy/probe/*.json)."""

import json
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[5]
sys.path.insert(0, str(REPO / "fortran" / "eftcamb_test" / "audit" / "checks"))
import harness  # noqa: E402

D = harness.AUDIT_ROOT / "work" / "verify" / "ic1noisy" / "probe"
for f in sorted(D.glob("*.json")):
    r = json.loads(f.read_text())
    a = r["a_on"]
    print(f"=== {f.stem}: a_on={a:g}  A1={r['EFTpiA1']:.4g} B1={r['EFTpiB1']:.4g} C={r['EFTpiC']:.4g} "
          f"D1={r['EFTpiD1']:.4g} D2={r['EFTpiD2']:.3g} aH={r['adotoa']:.4g}")
    d = r["dense"]
    C, D1 = np.array(d["EFTpiC"]), np.array(d["EFTpiD1"])
    i0 = len(C) // 2
    print(f"   dense C (dx=5e-4) at -2..+2: " + " ".join(f"{c:+.3e}" for c in C[i0 - 2:i0 + 3]) +
          f" | median|C| over +-0.1 in ln a = {np.median(np.abs(C)):.3e}; slope sign changes C {int(np.sum(np.diff(np.sign(np.diff(C))) != 0))}"
          f" D1 {int(np.sum(np.diff(np.sign(np.diff(D1))) != 0))} (of {len(C)})")
    print("   k     M/(A aH^2)   B/(A aH)   -n_IC(D2)  dfridr-relerr  pi'/(aH pi) code   pi_code      pi_QS(code E)")
    for k, row in r["k"].items():
        ic = r["ic_actual"][k]
        print(f"   {float(k):5g} {row['M_over_A_H2']:11.3e} {row['B_over_AH']:10.3e} {-row['n_IC_D2']:+10.3e} "
              f"{row['dfridr_relerr_D2']:10.2e}     {ic['pidot_over_aHpi'][0]:+11.3e}   {ic['pi'][0]:+.3e}   "
              f"{ic['pi_QS'][0]:+.3e}")
