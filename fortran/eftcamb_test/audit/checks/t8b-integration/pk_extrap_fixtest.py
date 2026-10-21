"""Repro (track t8b-integration): EFTCAMB MOD in camb/results.py get_matter_power_interpolator(extrap_kmax=...)
biases the power-law extrapolation beyond kmax, also for plain GR / LCDM (no EFT involved).

Compares, for LCDM at z=0, the extrapolated P(k) obtained with transfer kmax=2/Mpc (kh_max ~ 3.1 h/Mpc) against
the 'truth' computed directly with kmax=100/Mpc, for the code under audit (A0) and upstream CAMB v2 (B0).
Run: /Users/marco/Software/miniforge3/bin/python pk-extrap-slope-bias.py
"""

import json
import os
import subprocess
import sys

AUDIT = ("/private/tmp/claude-501/-Users-marco-IW-eftcamb-update-EFTCAMB/"
         "2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad/audit/trees/")

WORKER = r"""
import sys, json
import numpy as np
sys.path.insert(0, sys.argv[1])
import camb
def pars(kmax, nl):
    return camb.set_params(H0=67.5, ombh2=0.022, omch2=0.122, As=2.1e-9, ns=0.965, WantCls=False,
                           WantTransfer=True, kmax=kmax, redshifts=[0.0], NonLinear=nl)
out = {}
for nl in ("NonLinear_none", "NonLinear_both"):
    lo = camb.get_results(pars(2.0, nl))
    hi = camb.get_results(pars(100.0, nl))
    kh_lo = lo.get_linear_matter_power_spectrum(nonlinear=nl != "NonLinear_none")[0][-1]
    k = kh_lo * np.array([1.5, 2.0, 5.0, 10.0])
    P_ex = lo.get_matter_power_interpolator(nonlinear=nl != "NonLinear_none", extrap_kmax=30 * kh_lo, silent=True).P(0.0, k)
    P_true = hi.get_matter_power_interpolator(nonlinear=nl != "NonLinear_none", silent=True).P(0.0, k)
    out[nl] = {"kh_max": float(kh_lo), "k_over_kmax": [1.5, 2.0, 5.0, 10.0],
               "rel_err_extrapolation": (np.asarray(P_ex) / np.asarray(P_true) - 1).tolist()}
print("@@" + json.dumps(out))
"""

res = {}
for tree in ("A0", "B0", "FIX"):
    env = dict(os.environ, OMP_NUM_THREADS="2")
    env.pop("PYTHONPATH", None)
    p = subprocess.run([sys.executable, "-c", WORKER, (AUDIT + tree if tree != 'FIX' else '/private/tmp/claude-501/-Users-marco-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad/audit/work/t8b-integration/A0_fixpy')], capture_output=True, text=True, env=env)
    line = [ln for ln in p.stdout.splitlines() if ln.startswith("@@")]
    if not line:
        print(tree, "failed", p.stderr[-2000:])
        sys.exit(1)
    res[tree] = json.loads(line[0][2:])

for nl in ("NonLinear_none", "NonLinear_both"):
    print(f"{nl}: kh_max = {res['A0'][nl]['kh_max']:.3f} h/Mpc; fractional error of extrapolated P(k,z=0) vs direct kmax=100 run")
    print("   k/kh_max             ", res["A0"][nl]["k_over_kmax"])
    print("   observed  (A0, EFTCAMB)", ["%+.3f" % x for x in res["A0"][nl]["rel_err_extrapolation"]])
    print("   expected  (B0, upstream)", ["%+.3f" % x for x in res["B0"][nl]["rel_err_extrapolation"]])
for nl in ("NonLinear_none", "NonLinear_both"):
    print("   candidate fix (A0+py)  ", ["%+.3f" % x for x in res["FIX"][nl]["rel_err_extrapolation"]])
