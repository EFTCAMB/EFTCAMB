"""Verify finding extgal-adotoa-stale-grhonu: low-ell Cl impact (the model's perturbation integration fails at the
default l_max on every tree, so only l <= 50 without transfer/lensing is computable). A0 vs FIX3, default stability.
"""

from __future__ import annotations

import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[5]
sys.path.insert(0, str(REPO / "fortran/eftcamb_test/audit/checks"))
import harness  # noqa: E402

SP = Path("/private/tmp/claude-501/-Users-marco-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45"
          "/scratchpad/audit")
TREES = {"A0": SP / "trees/A0", "FIX3": SP / "work/t2d-static/FIX3"}
OUT = SP / "work/verify/extgal/lowell"
INI = REPO / "fortran/eftcamb_test/audit/inis/008p7-ext-galileon.ini"
BASE = {"get_tensor_cls": "F", "get_transfer": "F", "do_lensing": "F", "l_max_scalar": "50",
        "num_redshiftwindows": "0"}
CASES = {"nu0.01": {}, "mnu0.06": {"omnuh2": "0.000645"},
         "massless": {"omnuh2": "0", "massive_neutrinos": "0", "massless_neutrinos": "3.046"}}


def main() -> None:
    jobs = [(c, t) for c in CASES for t in TREES]

    def one(j):
        c, t = j
        return c, t, harness.run_model(TREES[t], INI, OUT / f"{c}_{t}.npz", overrides={**BASE, **CASES[c]}, threads=1)

    with ThreadPoolExecutor(4) as ex:
        for c, t, st in ex.map(one, jobs):
            print(c, t, st["status"], st["elapsed"], st.get("log_tail", [""])[-3:] if st["status"] != "ok" else "")
    for c in CASES:
        a, b = OUT / f"{c}_A0.npz", OUT / f"{c}_FIX3.npz"
        if not (a.exists() and b.exists()):
            print(c, "missing")
            continue
        A, B = harness.load(a), harness.load(b)
        ca, cb = A["cl_unlensed_scalar"], B["cl_unlensed_scalar"]
        n = min(len(ca), len(cb), 51)
        r = (ca[2:n] - cb[2:n]) / cb[2:n]
        print(f"== {c}: A0/FIX3-1 unlensed scalar, l=2..{n - 1}: max|TT| = {np.max(np.abs(r[:, 0])):.3e} "
              f"(l={2 + np.argmax(np.abs(r[:, 0]))}), max|EE| = {np.max(np.abs(r[:, 1])):.3e}; "
              f"TT at l=2,10,30,50: {[f'{r[i - 2, 0]:+.2e}' for i in (2, 10, 30, n - 1)]}")


if __name__ == "__main__":
    main()
