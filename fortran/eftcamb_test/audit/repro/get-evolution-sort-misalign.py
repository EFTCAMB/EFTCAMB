"""Repro: EFTCAMB.get_evolution sorts the requested conformal times (np.sort) and returns the outputs in sorted order,
so for non-increasing input the rows no longer correspond to the caller's times. get_scale_evolution with a
decreasing scale factor array (e.g. a = 1/(1+z) for increasing z) returns the evolution reversed, silently.

Usage: python get-evolution-sort-misalign.py [--tree PATH]
Prints requested vs returned a / tau and a physical field (adotoa) vs get_eft_functions at the requested a.
"""

import argparse
import os
import sys
import tempfile
from pathlib import Path

import numpy as np

ROOT = "/private/tmp/claude-501/-Users-marco-Desktop-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad/audit"
PARAMS = Path(__file__).resolve().parents[2] / "parameters"

ap = argparse.ArgumentParser()
ap.add_argument("--tree", default=f"{ROOT}/trees/A0")
args = ap.parse_args()
os.environ.setdefault("OMP_NUM_THREADS", "2")
sys.path.insert(0, args.tree)
import camb  # noqa: E402

tmp = tempfile.mkdtemp(prefix="evsort_")
ini = os.path.join(tmp, "fR.ini")
Path(ini).write_text(f"DEFAULT({PARAMS / '4_Des_F_R_LCDM_1.ini'})\noutput_root = {tmp}/fR\n"
                     "EFTCAMB_write_background = F\nfeedback_level = 0\n")
cp = camb.read_ini(ini, no_validate=True)
res = camb.get_background(cp)
print("camb from", camb.__file__)

z = np.array([0.0, 0.5, 1.0, 2.0, 5.0])  # increasing z -> decreasing a
a = 1 / (1 + z)
_, ev = cp.EFTCAMB.get_scale_evolution(res, [0.1], a)
_, ef = cp.EFTCAMB.get_eft_functions(res, a)
print("\nget_scale_evolution(results, k=[0.1], a=1/(1+z)), z =", z.tolist())
print(f"{'requested a':>12s} {'returned a':>12s} {'adotoa returned':>16s} {'adotoa at requested a':>22s}")
for i in range(len(a)):
    print(f"{a[i]:12.5f} {ev['a'][0, i]:12.5f} {ev['adotoa'][0, i]:16.6e} {ef['adotoa'][i]:22.6e}")
err = np.max(np.abs(ev["a"][0] - a) / a)
print(f"observed max |a_returned - a_requested|/a = {err:.3e}; expected <= 1e-6")

eta = res.conformal_time(np.array([3.0, 0.0, 1.0]))
_, ev2 = cp.EFTCAMB.get_evolution(res, [0.1], eta)
print("\nget_evolution(results, k=[0.1], eta) with eta =", np.round(eta, 3).tolist())
print("returned tau =", np.round(ev2["tau"][0], 3).tolist())
print(f"observed max |tau_returned - eta_requested| = {np.max(np.abs(ev2['tau'][0] - eta)):.3e} Mpc; expected 0")
