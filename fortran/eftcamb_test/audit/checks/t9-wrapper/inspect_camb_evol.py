"""Compare EFTCAMB get_evolution, upstream CAMB get_time_evolution and the DEBUG source dump at one k.

Separates wrapper defects from the intrinsic difference between a stand-alone evolution and the source evolution.
Usage: python inspect_camb_evol.py MODEL [tree]"""

import sys

import numpy as np

from t9lib import WORK, run_snippet

stem = sys.argv[1]
tree = sys.argv[2] if len(sys.argv) > 2 else "A0"
tag = tree.rstrip("/").split("/")[-1]
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
    fl, ev = pars.EFTCAMB.get_evolution(results, float(kk), taus)
    ct = results.get_time_evolution(float(kk), taus, ['delta_cdm', 'delta_photon', 'delta_baryon', 'a'])
    out[float(kk)] = dict(a=dump['FRW__a'][s], tau=taus,
                         clxc_dump=dump['EinsteinSol__clxc'][s], clxg_dump=dump['EinsteinSol__clxg'][s],
                         clxc_eft=ev[0]['clxc'], clxg_eft=ev[0]['clxg'],
                         clxc_camb=ct[:, 0], clxg_camb=ct[:, 1], a_camb=ct[:, 3])
np.save(d / 'inspect_camb_{tag}.npy', out, allow_pickle=True)
emit(done=True)
"""
r = run_snippet(CODE, tree_name=tree, threads=2, timeout=1500)
print(r["status"], r["result"], "\n".join(r["tail"][-6:]) if r["status"] != "ok" else "")
out = np.load(WORK / "dbg" / stem / f"inspect_camb_{tag}.npy", allow_pickle=True).item()


def rel(x, y):
    return float(np.max(np.abs(x - y)) / np.max(np.abs(y)))


for kk, o in out.items():
    print(f"k={kk:.4g}: clxc eft-vs-dump {rel(o['clxc_eft'], o['clxc_dump']):.2e}  camb-vs-dump "
          f"{rel(o['clxc_camb'], o['clxc_dump']):.2e}  eft-vs-camb {rel(o['clxc_eft'], o['clxc_camb']):.2e} | "
          f"clxg eft-vs-dump {rel(o['clxg_eft'], o['clxg_dump']):.2e}  camb-vs-dump "
          f"{rel(o['clxg_camb'], o['clxg_dump']):.2e}  eft-vs-camb {rel(o['clxg_eft'], o['clxg_camb']):.2e}")
