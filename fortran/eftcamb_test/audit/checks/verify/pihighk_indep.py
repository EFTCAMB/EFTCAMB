"""Verifier (get-evolution-pi-highk-tolerance): independent convergence ladder for the pi sector at k ~ 1/Mpc.

Differences from the reporter's repro:
  * tolerance ladder IntTolBoost in {1, 2, 3, 4, 5.6, 8, 12} with reference 12 (not a single 1-vs-8 pair), to
    check that the deviation scales like the integrator tolerance (integration error signature) and that 8 is
    converged;
  * a second, upstream code path: CAMBdata.get_time_evolution (camb_python.f90 GetOutputEvolutionFork) for the
    Weyl potential (EFT: 0.5*k*EFTLensing), delta_cdm and etak, on the same conformal times;
  * GR control: same model ini with EFTflag = 0 through both paths;
  * a different k (0.5/Mpc) and a sparser set of output stops (every 7th time) to check sampling dependence.

Usage: python pihighk_indep.py TREE MODEL_INI OUT_JSON [--eftflag0]
Runs in-process under harness.slot(1) with OMP_NUM_THREADS=1.
"""

import json
import os
import sys
import tempfile
from pathlib import Path

os.environ["OMP_NUM_THREADS"] = "1"
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import harness  # noqa: E402

tree, ini, out_json = sys.argv[1], sys.argv[2], sys.argv[3]
gr = "--eftflag0" in sys.argv
sys.path.insert(0, tree)
import numpy as np  # noqa: E402

import camb  # noqa: E402

assert Path(camb.__file__).resolve().is_relative_to(Path(tree).resolve()), camb.__file__

BOOSTS = [1.0, 2.0, 3.0, 4.0, 5.6, 8.0, 12.0]
KS = [1.3172836143, 0.5]
A_FINE = np.logspace(np.log10(0.02), 0, 400)
A_COARSE = A_FINE[::7]  # exact subset: same times requested with fewer output stops
EFT_FIELDS = ["pi", "pidot", "Psi", "Phi", "mu"]
CAMB_VARS = ["Weyl", "delta_cdm", "etak"]

tmp = tempfile.mkdtemp(dir=os.path.dirname(os.path.abspath(out_json)))
fast = os.path.join(tmp, "model.ini")
extra = "EFTflag = 0\n" if gr else ""
with open(fast, "w") as f:
    f.write(
        f"DEFAULT({ini})\nl_max_scalar = 1500\nget_tensor_cls = F\nget_transfer = F\n"
        f"EFTCAMB_write_background = F\nfeedback_level = 0\n{extra}"
    )

res = {}
with harness.slot(1):
    for boost in BOOSTS:
        p = camb.read_ini(fast)
        p.Accuracy.IntTolBoost = boost
        r = camb.get_results(p)
        entry = {}
        for k in KS:
            for gname, agrid in (("fine", A_FINE), ("coarse", A_COARSE)):
                eta = r.conformal_time(1 / agrid - 1)
                fields, arr = p.EFTCAMB.get_scale_evolution(r, k, agrid)
                d = {f: np.asarray(arr[f][0], dtype=float) for f in EFT_FIELDS}
                ev = r.get_time_evolution(k, eta, vars=CAMB_VARS, lAccuracyBoost=0)
                for i, v in enumerate(CAMB_VARS):
                    d["camb_" + v] = np.asarray(ev[:, i], dtype=float)
                entry[(k, gname)] = d
        res[boost] = entry

ref = res[BOOSTS[-1]]
report = {"tree": tree, "ini": ini, "gr_control": gr, "boosts": BOOSTS, "ks": KS, "dev": {}, "samples": {}}
for k in KS:
    for gname in ("fine", "coarse"):
        key = f"k={k:.4g},{gname}"
        report["dev"][key] = {}
        for name in ref[(k, gname)]:
            refv = ref[(k, gname)][name]
            m = float(np.max(np.abs(refv)))
            report["dev"][key][name] = [
                float(np.max(np.abs(res[b][(k, gname)][name] - refv)) / m) if m > 0 else None for b in BOOSTS
            ]
# sampling dependence at default and at boost 8: same times, every 7th stop only
for b in (1.0, 8.0):
    for k in KS:
        key = f"sampling,boost={b},k={k:.4g}"
        report["dev"][key] = {}
        for name in ("pi", "Psi", "camb_Weyl"):
            fine = res[b][(k, "fine")][name]
            coarse = res[b][(k, "coarse")][name]
            m = float(np.max(np.abs(ref[(k, "fine")][name])))
            report["dev"][key][name] = float(np.max(np.abs(fine[::7] - coarse)) / m)
# sample values
k = KS[0]
for i in (0, 100, 200, 300, 350, 399):
    report["samples"][f"a={A_FINE[i]:.4f}"] = {
        name: [float(res[b][(k, "fine")][name][i]) for b in BOOSTS] for name in ("pi", "Psi", "camb_Weyl")
    }
with open(out_json, "w") as f:
    json.dump(report, f, indent=1)
for key, d in report["dev"].items():
    print(key)
    for name, v in d.items():
        if isinstance(v, list):
            print(f"   {name:14s}", " ".join(f"{x:9.2e}" if x is not None else "     None" for x in v))
        else:
            print(f"   {name:14s} {v}")
print("samples k=1.317 (boosts", BOOSTS, ")")
for a, d in report["samples"].items():
    print(a, {n: [f"{x:.3e}" for x in v] for n, v in d.items()})
