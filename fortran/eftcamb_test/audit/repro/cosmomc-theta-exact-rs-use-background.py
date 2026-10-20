"""Repro (track t8b-integration): with EFTCAMB_use_background = T (Horndeski, EFTflag=5), CAMBdata_CosmomcTheta
(results.f90 EFTCAMB MOD) takes r_s from the EFT table compute_SoundHorizon, which integrates the EXACT baryon
loading R = 3 rho_b / (4 rho_g), instead of CAMB's cosmomc_theta definition dsound_da_approx (R = 3e4 a ombh2).
The same physical model therefore returns a theta_MC 0.1% lower when the table is used.

Observed: results.cosmomc_theta() with use_background=T vs F (identical background, H(z) equal to 1e-15).
Expected: equal (both should follow the CAMB theta_MC definition), i.e. the T value should match an
independent quadrature with the approximate R.
Run: /Users/marco/Software/miniforge3/bin/python cosmomc-theta-exact-rs-use-background.py
"""

import json
import os
import subprocess
import sys
import tempfile

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "checks"))
try:  # hold a machine-wide CPU slot of the audit harness when available
    from harness import slot as _slot
except Exception:  # noqa: BLE001
    from contextlib import nullcontext as _slot


def _run(*args, **kwargs):
    with _slot(1):
        return subprocess.run(*args, preexec_fn=lambda: os.nice(10), **kwargs)


TREE = ("/private/tmp/claude-501/-Users-marco-Desktop-IW-eftcamb-update-EFTCAMB/"
        "2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad/audit/trees/A0")
PARAMS = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "parameters"))

WORKER = r"""
import sys, json, os
import numpy as np
from scipy.integrate import quad
sys.path.insert(0, sys.argv[1])
import camb

out = {}
for ub in ("F", "T"):
    ini = os.path.join(sys.argv[3], "theta_ub%s.ini" % ub)
    open(ini, "w").write("DEFAULT(%s/5_hdsk_TG_1.ini)\nEFTCAMB_use_background = %s\nfeedback_level = 0\n" % (sys.argv[2], ub))
    cp = camb.read_ini(ini, no_validate=True)
    res = camb.get_background(cp)
    ombh2, omdmh2 = cp.ombh2, cp.omch2 + cp.omnuh2
    zs = 1048 * (1 + 0.00124 * ombh2**-0.738) * (1 + (0.0783 * ombh2**-0.238 / (1 + 39.5 * ombh2**0.763))
                                                 * (omdmh2 + ombh2) ** (0.560 / (1 + 21.1 * ombh2**1.81)))
    ast = 1 / (1 + zs)
    dtauda = lambda a: 1.0 / (a * a * float(res.h_of_z(1 / a - 1)))
    pts = np.logspace(-8, np.log10(ast), 40)
    rs = sum(quad(lambda a: dtauda(a) / np.sqrt(3 * (1 + 3.0e4 * a * ombh2)), pts[i], pts[i + 1], epsrel=1e-10)[0]
             for i in range(39))
    DA = float(res.angular_diameter_distance(zs)) / ast
    out[ub] = dict(theta_code=res.cosmomc_theta(), theta_def=rs / DA, H1=float(res.hubble_parameter(1.0)))
print("@@" + json.dumps(out))
"""

env = dict(os.environ, OMP_NUM_THREADS="1")
env.pop("PYTHONPATH", None)
with tempfile.TemporaryDirectory() as tmp:
    p = _run([sys.executable, "-c", WORKER, TREE, PARAMS, tmp], capture_output=True, text=True, env=env,
                       cwd=tmp)
line = [ln for ln in p.stdout.splitlines() if ln.startswith("@@")]
if not line:
    print(p.stdout[-2000:], p.stderr[-2000:])
    sys.exit(1)
r = json.loads(line[0][2:])
print("5_hdsk_TG_1 (EFTflag=5), cosmomc_theta():")
print(f"   use_background=F : {r['F']['theta_code']:.9f}   (CAMB-definition quadrature {r['F']['theta_def']:.9f})")
print(f"   use_background=T : {r['T']['theta_code']:.9f}   (CAMB-definition quadrature {r['T']['theta_def']:.9f})")
print(f"   H(z=1) F/T       : {r['F']['H1']:.10f} / {r['T']['H1']:.10f}")
print(f"   observed  theta(T)/theta(F) - 1 = {r['T']['theta_code'] / r['F']['theta_code'] - 1:+.3e}")
print("   expected  theta(T)/theta(F) - 1 = ~0 (|.| < 1e-5)")
print(f"   100*Delta theta = {100 * (r['T']['theta_code'] - r['F']['theta_code']):+.5f} (Planck sigma(100 theta_MC) ~ 0.0003)")
