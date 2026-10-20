"""Dump get_scale_evolution arrays (A0 or FX) for one model to npz, for inspecting the pidotdot verification."""
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))
import harness  # noqa: E402

AUDIT = str(harness.AUDIT_ROOT) + "/"
TREES = {"A0": AUDIT + "trees/A0", "FX": AUDIT + "work/t8a-equations/FX", "P": AUDIT + "work/t8a-equations/P"}
OUT = AUDIT + "work/verify/pidotdot/dbg"
WORKER = r"""
import sys
import numpy as np
sys.path.insert(0, sys.argv[1])
import camb
pars = camb.read_ini(sys.argv[2], no_validate=True)
pars.WantTensors = False
pars.InitPower.r = 0
pars.set_for_lmax(600)
res = camb.get_results(pars)
ks = [float(x) for x in sys.argv[4].split(",")]
a = np.logspace(float(sys.argv[5]), 0, int(sys.argv[6]))
fields, arr = pars.EFTCAMB.get_scale_evolution(res, ks, a)
out = {f: np.array([arr[i][f] for i in range(len(ks))]) for f in fields}
out["a_req"] = a
out["ks"] = np.array(ks)
out["h0_mpc"] = pars.H0 / 299792.458
np.savez(sys.argv[3], **out)
print("ok")
"""


def dump(tree, model, ks="0.005,0.02,0.1", loga0="-1", n="2000", env_extra=None):
    os.makedirs(OUT, exist_ok=True)
    if model.endswith(".ini"):
        ini, model = model, os.path.basename(model)[:-4]
    else:
        ini = os.path.join(str(harness.PARAMS_DIR), model + ".ini")
    out = f"{OUT}/{model}_{tree}.npz"
    env = dict(os.environ, OMP_NUM_THREADS="1")
    env.pop("PYTHONPATH", None)
    env.update(env_extra or {})
    with harness.slot(1):
        p = subprocess.run([sys.executable, "-c", WORKER, TREES[tree], ini, out, ks, loga0, n], env=env,
                           capture_output=True, text=True, cwd=str(harness.PARAMS_DIR) + "/..", timeout=1200)
    if p.returncode:
        print(p.stderr[-2000:])
    return out


if __name__ == "__main__":
    from concurrent.futures import ThreadPoolExecutor

    tree = sys.argv[1]
    models = sys.argv[2:]
    with ThreadPoolExecutor(4) as ex:
        print(list(ex.map(lambda m: dump(tree, m), models)))
