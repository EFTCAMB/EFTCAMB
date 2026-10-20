"""Independent verification of finding output-pidotdot-zero (EFTeomNdot/EFTeomM evaluated with pi''=0 at output).

Method (different from the repro, which compares EFTeomNdot with d(EFTeomN)/dtau):
  get_scale_evolution (python fork, yprime zeroed before derivs, same code path as output()) on
    A0 (code under audit) and FX (A0 + proposed fix, built by t8a), and compare time derivatives that must be
  consistent by construction:
    ISW   : EFTISW    vs d(EFTLensing)/dtau     (EFTISW is the time derivative of the Weyl source sigma'+eta k)
    PsiDot: PsiDot    vs dPsi/dtau              (PsiDot uses sigma'' -> EFTeomNdot)
    sigdd : sigmadotdot vs d(sigmadot)/dtau     (sigma'' uses EFTeomNdot)
    PhiDot: PhiDot    vs dPhi/dtau              (CONTROL: no sigma''/EFTeomNdot, must be consistent in both trees)
  Affected models (Gamma4 + 2 Gamma5 != 0): Horava, pure-EFT gamma4/gamma5 without Horndeski constraint.
  Controls (Gamma4 + 2 Gamma5 == 0): GR, Omega const, gamma3, RPH alpha_T (Horndeski).

Run: /Users/marco/Software/miniforge3/bin/python pidotdot_indep.py [model ...]
"""

import json
import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))
import harness  # noqa: E402

AUDIT = str(harness.AUDIT_ROOT) + "/"
TREES = {"A0": AUDIT + "trees/A0", "FX": AUDIT + "work/t8a-equations/FX"}
PARAMS = str(harness.PARAMS_DIR)
OUT = AUDIT + "work/verify/pidotdot"

AFFECTED = ["5_Horava_eta_3", "5_Horava_eta_4", "2_PEFT_gamma4_const_1", "2_PEFT_gamma5_const_1"]
CONTROLS = ["1_EFT_GR", "2_PEFT_Omega_const_1", "2_PEFT_gamma3_const_1", "3_RPH_KT_power_1"]

WORKER = r"""
import sys, json
import numpy as np
sys.path.insert(0, sys.argv[1])
import camb
pars = camb.read_ini(sys.argv[2], no_validate=True)
pars.WantTensors = False
pars.InitPower.r = 0
pars.set_for_lmax(600)
res = camb.get_results(pars)
ks = [0.005, 0.02, 0.1]
a = np.logspace(-1.0, 0, 2000)
fields, arr = pars.EFTCAMB.get_scale_evolution(res, ks, a)
out = {}
sl = slice(10, -10)
def resid(val, ref):
    return float(np.max(np.abs(val[sl] - ref[sl])) / np.max(np.abs(ref[sl])))
for ik, k in enumerate(ks):
    c = arr[ik]
    tau = c["tau"]
    g = lambda f: np.gradient(c[f], tau, edge_order=2)
    out[str(k)] = {
        "ISW": resid(c["EFTISW"], g("EFTLensing")),
        "PsiDot": resid(c["PsiDot"], g("Psi")),
        "sigdd": resid(c["sigmadotdot"], g("sigmadot")),
        "PhiDot": resid(c["PhiDot"], g("Phi")),
        "pidd_fd": resid(c["pidotdot"], g("pidot")),
        "G4p2G5_max": float(np.max(np.abs(c["EFTGamma4V"] + 2 * c["EFTGamma5V"]))),
        "ISW_today": float(c["EFTISW"][-1]),
    }
print("JSON" + json.dumps(out))
"""


def run_one(tree_name, model):
    env = dict(os.environ, OMP_NUM_THREADS="1")
    env.pop("PYTHONPATH", None)
    ini = os.path.join(PARAMS, model + ".ini")
    with harness.slot(1):
        p = subprocess.run([sys.executable, "-c", WORKER, TREES[tree_name], ini], env=env, capture_output=True,
                           text=True, cwd=PARAMS + "/..", preexec_fn=lambda: os.nice(10), timeout=1200)
    line = [ln for ln in p.stdout.splitlines() if ln.startswith("JSON")]
    if p.returncode or not line:
        return tree_name, model, {"error": p.stderr[-1500:]}
    return tree_name, model, json.loads(line[-1][4:])


def main():
    os.makedirs(OUT, exist_ok=True)
    models = sys.argv[1:] or AFFECTED + CONTROLS
    jobs = [(t, m) for m in models for t in TREES]
    results = {}
    with ThreadPoolExecutor(4) as ex:
        for t, m, r in ex.map(lambda tm: run_one(*tm), jobs):
            results.setdefault(m, {})[t] = r
    with open(os.path.join(OUT, "indep_" + "_".join(models)[:80] + ".json"), "w") as f:
        json.dump(results, f, indent=1)
    print("relative residual max|cached - FD| / max|FD|, a in [0.1,1]; columns per k = 0.005, 0.02, 0.1 /Mpc")
    for m in models:
        for t in TREES:
            r = results[m][t]
            if "error" in r:
                print(f"{m:24s} {t}: ERROR {r['error'][-300:]}")
                continue
            ks = list(r)
            row = "  ".join(
                f"{q}=" + ",".join(f"{r[k][q]:.1e}" for k in ks) for q in ("ISW", "PsiDot", "sigdd", "PhiDot", "pidd_fd")
            )
            print(f"{m:24s} {t}: |G4+2G5|max={r[ks[0]]['G4p2G5_max']:.1e}  {row}")


if __name__ == "__main__":
    main()
