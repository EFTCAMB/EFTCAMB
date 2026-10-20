"""Verifier: reproduce t8b's TG_3 theta round trip (set_cosmology(cosmomc_theta=...)) and print theta(H0) across the
bracket, to see whether H0_rec = 43.16 is a genuine root or a flat-function/noise artefact."""

import json
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))
import harness  # noqa: E402

TREE = str(harness.TREES["A0"])
INI = str(harness.PARAMS_DIR / "5_hdsk_TG_3.ini")
OUT = ("/private/tmp/claude-501/-Users-marco-Desktop-IW-eftcamb-update-EFTCAMB/"
       "2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad/audit/work/verify/hdsk_h0")
WORKER = r"""
import sys, json
import numpy as np
sys.path.insert(0, sys.argv[1])
import camb
cp = camb.read_ini(sys.argv[2], no_validate=True)
kw = dict(ombh2=cp.ombh2, omch2=cp.omch2, omk=cp.omk, omnuh2_active=cp.omnuh2, mnu=None,
          num_massive_neutrinos=max(1, int(cp.num_nu_massive)), nnu=cp.N_eff, TCMB=cp.TCMB, YHe=cp.YHe)
H0 = cp.H0
cp.set_cosmology(H0=H0, **kw)
th = camb.get_background(cp, no_thermo=True).cosmomc_theta()
out = {"theta_target": th, "scan": {}}
for h in (43.05, 43.1, 43.16, 43.3, 45.0, 60.0, 70.87, 100.0):
    cp.H0 = h
    b = camb.get_background(cp, no_thermo=True)
    out["scan"][str(h)] = [b.cosmomc_theta() - th, float(b.hubble_parameter(0.0))]
try:
    cp.set_cosmology(cosmomc_theta=th, **kw)
    out["H0_rec"] = cp.H0
except Exception as e:
    out["H0_rec"] = repr(e)
print("@@" + json.dumps(out))
"""
env = dict(os.environ, OMP_NUM_THREADS="1")
env.pop("PYTHONPATH", None)
with harness.slot(1):
    p = subprocess.run([sys.executable, "-c", WORKER, TREE, INI], capture_output=True, text=True, env=env, cwd=OUT,
                       preexec_fn=lambda: os.nice(10), timeout=900)
line = [ln for ln in p.stdout.splitlines() if ln.startswith("@@")]
print(json.dumps(json.loads(line[0][2:]), indent=1) if line else p.stdout[-2000:] + p.stderr[-2000:])
