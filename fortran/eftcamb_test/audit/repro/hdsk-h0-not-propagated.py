"""Repro (track t8b-integration): for Horndeski (EFTflag=5) models without shooting, the EFT background fixes
its own H0 (008p0_Horndeski.f90 sets params_cache%h0 = H(a=1) of the integrated background) but CAMB keeps the
input CP%H0, which it still uses for h-units, sigma8 (R = 8 Mpc/h), halofit/HMcode Omega_m = omega_m/h^2,
curvature, derived parameters and the theta -> H0 solver.

Observed: results.Params.H0 vs results.hubble_parameter(0) for the test ini 5_hdsk_TG_1; sigma8 from CAMB vs
sigma(R = 8/h_true Mpc). theta solve with the model's own theta.
Expected: H0 identical; sigma8 = sigma(8 Mpc/h) with the model's h.
Run: /Users/marco/Software/miniforge3/bin/python hdsk-h0-not-propagated.py
"""

import json
import os
import subprocess
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "checks"))
try:  # hold a machine-wide CPU slot of the audit harness when available
    from harness import slot as _slot
except Exception:  # noqa: BLE001
    from contextlib import nullcontext as _slot


def _run(*args, **kwargs):
    with _slot(1):
        return subprocess.run(*args, preexec_fn=lambda: os.nice(10), **kwargs)


TREE = ("/private/tmp/claude-501/-Users-marco-IW-eftcamb-update-EFTCAMB/"
        "2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad/audit/trees/A0")
PARAMS = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "parameters"))

WORKER = r"""
import sys, json
import numpy as np
sys.path.insert(0, sys.argv[1])
import camb

cp = camb.read_ini(sys.argv[2] + "/5_hdsk_TG_1.ini", no_validate=True)
cp.WantCls = True
cp.set_matter_power(redshifts=[0.0], kmax=2.0)
res = camb.get_results(cp)
H0_true = float(res.hubble_parameter(0.0))
s8 = float(res.get_sigma8_0())
sR = float(np.atleast_1d(res.get_sigmaR(8.0 / (H0_true / 100), hubble_units=False, z_indices=[-1]))[-1])
th = res.cosmomc_theta()
# theta as a function of the input H0 (the theta -> H0 solvers vary CP%H0)
th_of_H0 = {}
for H0 in (60.0, 80.0):
    cp.H0 = H0
    th_of_H0[str(H0)] = camb.get_background(cp, no_thermo=True).cosmomc_theta()
out = dict(CP_H0=res.Params.H0, Params_H0=res.Params.H0, H0_true=H0_true, sigma8_camb=s8, sigma8_true_h=sR,
           derived_H0=res.get_derived_params().get("H0", None), theta=th, th_of_H0=th_of_H0)
print("@@" + json.dumps(out))
"""

env = dict(os.environ, OMP_NUM_THREADS="1")
env.pop("PYTHONPATH", None)
p = _run([sys.executable, "-c", WORKER, TREE, PARAMS], capture_output=True, text=True, env=env)
line = [ln for ln in p.stdout.splitlines() if ln.startswith("@@")]
if not line:
    print(p.stdout[-3000:], p.stderr[-3000:])
    sys.exit(1)
r = json.loads(line[0][2:])
print("5_hdsk_TG_1 (EFTflag=5, Horndeski_shooting=F)")
print(f"   CAMB H0 (results.Params.H0)          : {r['Params_H0']:.4f} km/s/Mpc   (observed, used for h units/sigma8)")
print(f"   H(z=0) of the EFT background (dtauda): {r['H0_true']:.4f} km/s/Mpc   (expected CAMB H0)")
print(f"   sigma8 reported (R = 8 Mpc/h_CAMB)   : {r['sigma8_camb']:.4f}  (observed)")
print(f"   sigma(R = 8 Mpc/h_true)              : {r['sigma8_true_h']:.4f}  (expected sigma8)")
print(f"   relative error sigma8                : {r['sigma8_camb'] / r['sigma8_true_h'] - 1:+.3f}")
print(f"   cosmomc_theta at input H0 = 60 / 80   : {r['th_of_H0']['60.0']:.9f} / {r['th_of_H0']['80.0']:.9f}"
      "  (expected: theta depends on H0 as in GR; identical values mean theta -> H0 solvers cannot work)")
