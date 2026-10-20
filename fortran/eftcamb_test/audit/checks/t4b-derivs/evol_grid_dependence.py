"""Track t4b-derivs: does get_scale_evolution output at fixed a depend on the other requested output times?

    python evol_grid_dependence.py TREE MODEL_INI K [fields...]
Prints pi (and other fields) at a few target a values, obtained with (1) only the targets requested,
(2) the targets embedded in a dense log grid, (3) the targets embedded in a 9-point stencil each.
"""

from __future__ import annotations

import os
import sys
import tempfile
from pathlib import Path

import numpy as np

tree, ini_in, k = sys.argv[1], sys.argv[2], float(sys.argv[3])
fields = sys.argv[4:] or ["pi", "pidot", "z", "etak", "clxc"]
os.environ["OMP_NUM_THREADS"] = "1"
sys.path.insert(0, tree)
sys.path.insert(1, str(Path(__file__).resolve().parents[1]))
import camb  # noqa: E402
import harness  # noqa: E402

tmp = Path(tempfile.mkdtemp(prefix="t4b_evgrid_"))
ini = tmp / "w.ini"
ini.write_text(f"DEFAULT({Path(ini_in).resolve()})\noutput_root = {tmp}/out\nEFTCAMB_write_background = F\n"
               "EFTCAMB_feedback_level = 0\nfeedback_level = 0\n")
cp = camb.read_ini(str(ini), no_validate=True)
if cp.EFTCAMB.EFTflag == 0 and int(harness.parse_ini(ini_in).get("EFTflag", "0")) != 0:
    cp.EFTCAMB.initialize_parameters(cp, harness.parse_ini(str(ini)), print_header=False)
if os.environ.get("T4B_TOLBOOST"):
    cp.Accuracy.IntTolBoost = float(os.environ["T4B_TOLBOOST"])
with harness.slot(1):
    res = camb.get_results(cp) if os.environ.get("T4B_FULL") else camb.get_background(cp)
targets = np.array([0.05, 0.15, 0.3, 0.6, 0.9])


def ev(a):
    a = np.sort(a)
    _, arr = cp.EFTCAMB.get_scale_evolution(res, [k], a)
    out = {f: np.interp(targets, arr["a"][0], arr[f][0]) for f in fields}
    # exact match rows
    idx = [int(np.argmin(np.abs(arr["a"][0] - t))) for t in targets]
    return {f: arr[f][0][idx] for f in fields}, arr["a"][0][idx]


with harness.slot(1):
    r1, a1 = ev(targets)
    r2, a2 = ev(np.unique(np.concatenate([targets, np.logspace(-2, 0, 300)])))
    st = np.concatenate([targets * np.exp(0.01 * j) for j in range(-4, 5)])
    r3, a3 = ev(st)
print("tree", tree, "model", Path(ini_in).stem, "k", k)
for f in fields:
    for i, t in enumerate(targets):
        print(f"{f:8s} a={t:.3f} only={r1[f][i]:+.6e} dense={r2[f][i]:+.6e} stencil={r3[f][i]:+.6e}  "
              f"(a: {a1[i]:.6f} {a2[i]:.6f} {a3[i]:.6f})")
