"""Verifier (ic1-pidot-transient-noisy-coefficients): is the K-mimic C(a) oscillation tied to the background
interpolation grid (model_background_num_points = 1000 on ln a in [ln a_ini, 0])?  Locates the local minima of |C(a)|
in the ICPROF profile (a (1 + j 1e-4), |j| <= 400) and their distance to the nearest grid node, in units of the grid
width. Also reports the relative peak-to-peak oscillation of C across the profile."""

import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[5]
sys.path.insert(0, str(REPO / "fortran" / "eftcamb_test" / "audit" / "checks"))
from ic1noisy_icdump import OUT  # noqa: E402

CASES = {"km_on1e-8": 1e-9, "km_on1.01e-8": 1e-9, "km_on1e-7": 1e-9, "km_on1e-9_ai1e-11": 1e-11, "tg2": 1e-10}
NPTS = {"tg2": 60000}
for name, a_ini in CASES.items():
    prof = np.array([[float(x) for x in ln.split()[1:]] for ln in (OUT / f"{name}.log").read_text().splitlines()
                     if ln.startswith("ICPROF ")])
    x, c, d1 = np.log(prof[:, 0]), prof[:, 1], prof[:, 2]
    n = NPTS.get(name, 1000)
    w = -np.log(a_ini) / (n - 1)
    pos = (x - np.log(a_ini)) / w  # grid coordinate
    ac = np.abs(c)
    mins = [i for i in range(1, len(c) - 1) if ac[i] < ac[i - 1] and ac[i] < ac[i + 1]]
    frac = [float(pos[i] - np.round(pos[i])) for i in mins]
    print(f"=== {name}: a in [{prof[0,0]:.4g},{prof[-1,0]:.4g}], grid width {w:.5f} (in ln a), "
          f"profile spans {pos[-1]-pos[0]:.2f} cells")
    print(f"  C: min {c.min():+.3e} max {c.max():+.3e} median {np.median(c):+.3e}; "
          f"D1 relative spread {np.ptp(d1)/np.median(np.abs(d1)):.3e} (smooth power law slope "
          f"{np.polyfit(x, np.log(np.abs(d1)), 1)[0]:+.3f})")
    print(f"  {len(mins)} local minima of |C|; offset from nearest node (cells): "
          + " ".join(f"{f:+.3f}" for f in frac[:30]))
