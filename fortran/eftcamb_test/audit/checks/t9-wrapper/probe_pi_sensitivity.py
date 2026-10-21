"""Is the pi-field mismatch between get_evolution and the source dump just sensitivity to the output times?

Evaluates get_evolution twice at the dumped (k, tau): with the dumped taus and with taus perturbed by a relative
1e-11 (the rounding of the dump format), and reports the change of pi, pidot, Psi."""

import sys

import numpy as np

from t9lib import WORK, run_snippet

stem = sys.argv[1]
CODE = f"""
d = WORK / 'dbg' / {stem!r}
dump = np.load(d / 'dump_selected.npz')
pars = camb.read_ini(str(d / 'wrapper.ini'))
results = camb.get_results(pars)
k_all = dump['FRW__k']
res = {{}}
for kk in np.unique(k_all):
    s = k_all == kk
    taus = dump['FRW__tau'][s]
    fl, e1 = pars.EFTCAMB.get_evolution(results, float(kk), taus)
    fl, e2 = pars.EFTCAMB.get_evolution(results, float(kk), taus * (1 + 1e-11))
    fl, e3 = pars.EFTCAMB.get_evolution(results, float(kk), taus[::4])
    r = {{}}
    for f in ('pi', 'pidot', 'Psi', 'clxc'):
        x, y = e1[0][f], e2[0][f]
        r[f + '_perturbed_tau'] = float(np.max(np.abs(x - y)) / np.max(np.abs(x)))
        r[f + '_coarser_stops'] = float(np.max(np.abs(e1[0][f][::4] - e3[0][f])) / np.max(np.abs(x)))
        r[f + '_vs_dump'] = float(np.max(np.abs(x - dump[[n for n in dump.files if n.endswith('__' + f)][0]][s]))
                                  / np.max(np.abs(x)))
    res[float(kk)] = r
emit(res=res)
"""
r = run_snippet(CODE, tree_name="A0", threads=2, timeout=1500)
print(r["status"])
for kk, v in r["result"].get("res", {}).items():
    print(kk, {k: f"{x:.2e}" for k, x in v.items()})
