"""Track 9 (i) follow-up: at k ~ 1/Mpc the pi field returned by get_evolution changes O(1) with the set of output
times. Is this an integration-tolerance effect (converges with IntTolBoost) or a wrapper defect?

For one model and one k, evaluates pi, Psi on a late-time tau grid with (default tol, fine stops), (default tol,
every 4th stop), (IntTolBoost=8, fine), (IntTolBoost=8, coarse); prints where the differences are and the size of pi
relative to its maximum.

Usage: python probe_pi_highk.py MODEL [k]
"""

import sys

from t9lib import WORK, run_snippet

stem = sys.argv[1]
kk = float(sys.argv[2]) if len(sys.argv) > 2 else 1.3172836143
CODE = f"""
d = WORK / 'dbg' / {stem!r}
dump = np.load(d / 'dump_selected.npz')
s = dump['FRW__k'] == {kk!r}
taus = dump['FRW__tau'][s]
a = dump['FRW__a'][s]
res = {{}}
for boost in (1, 8):
    pars = camb.read_ini(str(d / 'wrapper.ini'))
    pars.Accuracy.IntTolBoost = boost
    results = camb.get_results(pars)
    fl, fine = pars.EFTCAMB.get_evolution(results, {kk!r}, taus)
    fl, coarse = pars.EFTCAMB.get_evolution(results, {kk!r}, taus[::4])
    res[boost] = (fine[0], coarse[0])
def rel(x, y, m):
    return np.abs(x - y) / m
out = {{}}
for f in ('pi', 'Psi', 'clxc'):
    f1, c1 = res[1][0][f], res[1][1][f]
    f8, c8 = res[8][0][f], res[8][1][f]
    m = np.max(np.abs(f8))
    late = a > 0.1
    out[f] = dict(
        max_abs=float(m),
        late_max_abs=float(np.max(np.abs(f8[late]))),
        tol1_fine_vs_coarse=float(np.max(rel(f1[::4], c1, m))),
        tol8_fine_vs_coarse=float(np.max(rel(f8[::4], c8, m))),
        tol1_fine_vs_tol8_fine=float(np.max(rel(f1, f8, m))),
        tol1_dump_vs_tol8_fine=float(np.max(rel(dump[[n for n in dump.files if n.endswith('__' + f)][0]][s], f8, m))),
        a_of_worst=float(a[::4][np.argmax(rel(f1[::4], c1, m))]),
        late_rel_tol1_fine_vs_tol8=float(np.max(np.abs(f1[late] - f8[late])) / np.max(np.abs(f8[late]))),
    )
emit(k={kk!r}, out=out, a_samples=[float(x) for x in a[::len(a)//8]],
     pi_tol8_samples=[float(x) for x in res[8][0]['pi'][::len(a)//8]],
     pi_tol1_samples=[float(x) for x in res[1][0]['pi'][::len(a)//8]])
"""
r = run_snippet(CODE, tree_name="A0", threads=2, timeout=2400)
print(r["status"])
import json  # noqa: E402

print(json.dumps(r["result"], indent=1))
if r["status"] != "ok":
    print("\n".join(r["tail"]))
