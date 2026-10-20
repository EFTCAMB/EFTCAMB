"""Verifier: does Horndeski_shooting=T (suggested in the finding's fix text) make the background H0 follow CP%H0?
Background only, TG_1 with Horndeski_shooting=T, input H0 = ini value and 85.5373."""

import json
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))
import harness  # noqa: E402

TREE = str(harness.TREES["A0"])
INI = os.path.join(HERE, "hdsk_TG_1_shooting.ini")
OUT = ("/private/tmp/claude-501/-Users-marco-Desktop-IW-eftcamb-update-EFTCAMB/"
       "2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad/audit/work/verify/hdsk_h0")
WORKER = r"""
import sys, json
sys.path.insert(0, sys.argv[1])
import camb
out = {}
for H0 in (None, 60.0):
    cp = camb.read_ini(sys.argv[2], no_validate=True)
    if H0: cp.H0 = H0
    try:
        b = camb.get_background(cp, no_thermo=True)
        out[str(cp.H0)] = dict(Params_H0=b.Params.H0, cache_h0=float(b.Params.EFTCAMB_parameter_cache.h0),
                               H_z0=float(b.hubble_parameter(0.0)))
    except Exception as e:
        out[str(cp.H0)] = repr(e)
print("@@" + json.dumps(out))
"""
env = dict(os.environ, OMP_NUM_THREADS="1")
env.pop("PYTHONPATH", None)
with harness.slot(1):
    p = subprocess.run([sys.executable, "-c", WORKER, TREE, INI], capture_output=True, text=True, env=env, cwd=OUT,
                       preexec_fn=lambda: os.nice(10), timeout=900)
line = [ln for ln in p.stdout.splitlines() if ln.startswith("@@")]
print(json.dumps(json.loads(line[0][2:]), indent=1) if line else p.stdout[-2000:] + p.stderr[-2000:])
