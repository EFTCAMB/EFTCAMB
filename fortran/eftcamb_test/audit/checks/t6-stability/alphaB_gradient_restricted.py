"""RPH alpha_B scalar-gradient boundary: analytic c_s^2 numerator evaluated on the cached background for
a >= model_background_a_ini (1e-8; below it the designer tables are clamped and Hdot is inconsistent), compared with
the code's bisected boundary (boundaries/scan_rph_alphaB_gradient_A0.json)."""
import json
import sys

import numpy as np

sys.path.insert(0, ".")
import t6lib  # noqa: E402
from boundary_scans import GRID, INIS, cached, gradient_numerator  # noqa: E402

r = cached("A0", INIS / "t6_RPH_KB_const.ini", {"RPHkineticity0": "0.1", "RPHbraiding0": "-0.1"})
m = r["a"] >= 1.0001e-8
rr = {k: (v[m] if isinstance(v, np.ndarray) and v.shape == r["a"].shape else v) for k, v in r.items()}
Bs = np.linspace(-0.35, -0.1, 250001)
lo, hi = -0.35, -0.1
for _ in range(60):
    mid = 0.5 * (lo + hi)
    if gradient_numerator(rr, mid).min() >= 0:
        hi = mid
    else:
        lo = mid
analytic = 0.5 * (lo + hi)
code = json.loads((t6lib.WORK / "boundaries" / "scan_rph_alphaB_gradient_A0.json").read_text())["code"]
amin = rr["a"][np.argmin(gradient_numerator(rr, analytic))]
print(json.dumps({"analytic_boundary_a>=1e-8": analytic, "code_boundary": code, "rel": abs(code - analytic) / abs(analytic),
                  "a_of_min": float(amin), "n_points": int(m.sum())}))
