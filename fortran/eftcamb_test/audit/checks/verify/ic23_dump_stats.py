"""Summarise IC23DUMP lines (series vs attractor pi at the turn-on) from an instrumented run log."""

import sys

import numpy as np

rows = []
for line in open(sys.argv[1], errors="replace"):
    if line.startswith("IC23DUMP"):
        p = line.split()
        try:
            rows.append([float(x) for x in p[2:11]])
        except ValueError:
            pass
r = np.array(rows)
if r.size == 0:
    print("no dump lines")
    sys.exit()
k, tau, a, ktau, ps, pa, dps, dpa, om = r.T
o = np.argsort(k)
print(f"n={len(k)}  a_turnon={np.unique(np.round(a, 12))[:3]}  Omega={om[0]:.4g}")
print(f"{'k':>10} {'tau':>10} {'k*tau':>10} {'pi_series':>12} {'pi_attr':>12} {'ratio':>10} {'dpi_ser':>12} {'dpi_attr':>12}")
sel = o[np.unique(np.linspace(0, len(o) - 1, 14).astype(int))]
for i in sel:
    ratio = ps[i] / pa[i] if pa[i] != 0 else np.inf
    print(f"{k[i]:10.4g} {tau[i]:10.4g} {ktau[i]:10.4g} {ps[i]:12.4e} {pa[i]:12.4e} {ratio:10.3g} {dps[i]:12.4e} {dpa[i]:12.4e}")
with np.errstate(divide="ignore", invalid="ignore"):
    rat = np.abs(ps / pa)
fin = np.isfinite(rat)
print(f"|series/attractor|: min {rat[fin].min():.3g} median {np.median(rat[fin]):.3g} max {rat[fin].max():.3g}; "
      f"k*tau range {ktau.min():.3g}..{ktau.max():.3g}")
