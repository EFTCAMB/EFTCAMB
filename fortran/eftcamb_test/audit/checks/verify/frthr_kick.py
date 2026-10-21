"""Verify fr-rgr-threshold-high-k-pk: growth rate f = dln(delta_c)/dlna = -k z/(adotoa delta_c) (synchronous gauge)
around the late turn-on (a=0.0931) for the default-threshold run vs the early-turn-on (threshold 1e-12) run,
from the get_scale_evolution output of frthr_evol.py.  Also a finite-difference f from delta_c itself as a cross-check.

Usage: python frthr_kick.py
"""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import harness  # noqa: E402

OUT = harness.AUDIT_ROOT / "work" / "verify" / "frthr" / "evol"
KS = [1.0, 3.0, 10.0]

d = {n: harness.load(OUT / f"{n}.npz") for n in ("def", "thr12")}
a = d["def"]["evol_a"]
a = a[0] if a.ndim > 1 else a
sel = np.where((a > 0.07) & (a < 0.16))[0]
for ik, k in enumerate(KS):
    print(f"k = {k} /Mpc   f = -k z/(aH delta_c)     [finite-diff dln(dc)/dlna between this and next point]")
    print(f"{'a':>8s} {'f_def':>9s} {'f_thr12':>9s} {'fd_def':>9s} {'fd_thr12':>9s}")
    for i in sel:
        row = []
        for n in ("def", "thr12"):
            e = d[n]
            row.append(-k * e["evol_z"][ik, i] / (e["evol_adotoa"][ik, i] * e["evol_clxc"][ik, i]))
        fd = [np.log(d[n]["evol_clxc"][ik, i + 1] / d[n]["evol_clxc"][ik, i]) / np.log(a[i + 1] / a[i])
              for n in ("def", "thr12")]
        print(f"{a[i]:8.4f} {row[0]:9.4f} {row[1]:9.4f} {fd[0]:9.4f} {fd[1]:9.4f}")
