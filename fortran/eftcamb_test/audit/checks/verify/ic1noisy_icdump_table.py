"""Verifier (ic1-pidot-transient-noisy-coefficients): compact per-k table of the IC dump (see ic1noisy_icdump.py)."""

import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[5]
sys.path.insert(0, str(REPO / "fortran" / "eftcamb_test" / "audit" / "checks"))
import harness  # noqa: E402
from ic1noisy_icdump import COLS, OUT  # noqa: E402

KS = [float(x) for x in sys.argv[2:]] or [0.05, 0.1, 0.3, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0]
for name in sys.argv[1].split(","):
    rows = np.array([[float(x) for x in ln.split()[1:]] for ln in (OUT / f"{name}.log").read_text().splitlines()
                     if ln.startswith("ICDUMP ")])
    d = {c: rows[:, i] for i, c in enumerate(COLS)}
    o = np.argsort(d["k"])
    d = {c: v[o] for c, v in d.items()}
    k, a, h0, aH = d["k"], d["a"], d["h0"], d["adotoa"]
    den = d["C"] + k**2 * d["D1"] + k**4 * d["D2"]
    pi = d["pi"]
    pd_fix = h0 * (d["E"] / den**2 * a * aH * (d["cd"] + k**2 * d["dd1"]) - d["Edot"] / den)
    pd_sm = h0 * (d["E"] / den**2 * a * aH * (d["cw"] + k**2 * d["d1w"]) - d["Edot"] / den)
    print(f"=== {name}: k range {k.min():.3g}..{k.max():.3g}; IC a range {a.min():.4g}..{a.max():.4g}")
    print("      k        a_IC        C          k2*D1      dlnC/dlna(Ridders) (wide)   n_code    n_fixtypo  n_smooth")
    for kk in KS:
        i = int(np.argmin(np.abs(k - kk)))
        print(f"  {k[i]:8.4g} {a[i]:10.4g} {d['C'][i]:+10.3e} {k[i]**2*d['D1'][i]:+10.3e}  "
              f"{d['cd'][i]*a[i]/d['C'][i]:+10.3e} {d['cw'][i]*a[i]/d['C'][i]:+10.3e}  "
              f"{d['pidot'][i]/(aH[i]*pi[i]):+10.3e} {pd_fix[i]/(aH[i]*pi[i]):+10.3e} {pd_sm[i]/(aH[i]*pi[i]):+10.3e}")
