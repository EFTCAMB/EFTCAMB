"""Run the spline (flag 9) Pure EFT Omega=a cases on the minimal-root-cause-patched tree and compare with the
analytic (PureEFTmodelOmega=2) run on A0. If the stated root cause is complete, differences drop to ~1e-12."""
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import harness  # noqa: E402
from verify_spline_cubic_model import CASES, INI, OUT  # noqa: E402

TREE = harness.AUDIT_ROOT / "work" / "verify" / "spline_cubic_fixmin"
names = ["A0_N6", "A0_N21", "A0_N10ext"]


def run(name):
    return name, harness.run_model(TREE, INI, OUT / f"fix_{name}.npz", overrides=CASES[name][1], threads=1)


with ThreadPoolExecutor(3) as ex:
    for name, st in ex.map(run, names):
        print(name, st["status"])
for name in names:
    d = harness.load(OUT / f"fix_{name}.npz")
    a = d["eft_a_grid"]
    m = a > 1e-3
    c = harness.compare(OUT / f"fix_{name}.npz", OUT / "A0_ref.npz")
    print(f"fix_{name}: max|OmegaV-a|={np.max(np.abs(d['eft_EFTOmegaV'][m] - a[m])):.2e} "
          f"max|OmegaP-1|={np.max(np.abs(d['eft_EFTOmegaP'][m] - 1)):.2e} "
          f"max|OmegaPP|={np.max(np.abs(d['eft_EFTOmegaPP'][m])):.2e} | worst Cl {c['worst_cl']:.3e} "
          f"P(k) {c['worst_pk']:.3e} sigma8 {c['sigma8']:.3e}")
