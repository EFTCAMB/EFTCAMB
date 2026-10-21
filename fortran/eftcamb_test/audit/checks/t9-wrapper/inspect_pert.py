"""Inspect the perturbation mismatch between get_evolution and the DEBUG cache dump for one model.

Re-runs the model on a tree (default A0) with the dbg_compare wrapper ini, evaluates get_evolution at the dumped
(k, tau) of dump_selected.npz, and prints the rows with the largest differences."""

import json
import sys

from t9lib import WORK, run_snippet

stem = sys.argv[1]
tree = sys.argv[2] if len(sys.argv) > 2 else "A0"
tag = tree.rstrip("/").split("/")[-1]
fields = sys.argv[3].split(",") if len(sys.argv) > 3 else ["pi", "Psi", "clxg", "z"]

CODE = f"""
d = WORK / 'dbg' / {stem!r}
dump = np.load(d / 'dump_selected.npz')
pars = camb.read_ini(str(d / 'wrapper.ini'))
results = camb.get_results(pars)
k_all = dump['FRW__k']
out = {{}}
for kk in np.unique(k_all):
    s = k_all == kk
    taus = dump['FRW__tau'][s]
    a = dump['FRW__a'][s]
    fl, ev = pars.EFTCAMB.get_evolution(results, float(kk), taus)
    rec = {{'a': a.tolist(), 'tau': taus.tolist()}}
    for f in {fields!r}:
        src = [n for n in dump.files if n.endswith('__' + f)][0]
        rec[f + '_dump'] = dump[src][s].tolist()
        rec[f + '_evol'] = ev[0][f].tolist()
    out[float(kk)] = rec
np.save(d / 'inspect_{tag}.npy', out, allow_pickle=True)
emit(done=True)
"""

r = run_snippet(CODE, tree_name=tree, threads=2, timeout=1500)
print(r["status"], r["result"], "\n".join(r["tail"][-5:]) if r["status"] != "ok" else "")
import numpy as np  # noqa: E402

out = np.load(WORK / "dbg" / stem / f"inspect_{tag}.npy", allow_pickle=True).item()
for kk, rec in out.items():
    a = np.array(rec["a"])
    print(f"k = {kk:.4g}: {len(a)} rows, a in [{a.min():.3g}, {a.max():.3g}]")
    for f in fields:
        x, y = np.array(rec[f + "_evol"]), np.array(rec[f + "_dump"])
        scale = np.max(np.abs(y)) or 1.0
        i = int(np.argmax(np.abs(x - y)))
        print(f"   {f:8s} maxdiff/max|dump| = {np.max(np.abs(x - y)) / scale:.3e} at a={a[i]:.4g} "
              f"(evol={x[i]:.6g}, dump={y[i]:.6g});  "
              + " ".join(f"a={a[j]:.3g}:{x[j]:.4g}/{y[j]:.4g}" for j in np.linspace(0, len(a) - 1, 6).astype(int)))
print(json.dumps({"tree": tree}))
