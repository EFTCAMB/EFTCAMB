"""Repro: with EFTflag=0 (plain GR) the code under audit (A0) differs from plain CAMB v2 (B0) in non-flat models.

equations.f90 (derivs, output block, "EFTCAMB MOD START: polterddot from the photon equations ...") replaces CAMB's
polterddot for every model, including EFTflag=0. The replacement drops the curvature factors EV%Kf(2) that CAMB
uses on the octgdot and Edot(3) terms, so only the temperature quadrupole source changes: EE is identical, TT/TE
differ at low ell in curved models. In flat models both expressions agree to round-off (~1e-13).

Run: /Users/marco/Software/miniforge3/bin/python polterddot-nonflat-gr.py
"""

import json
import os
import subprocess
import sys

TREES = (
    "/private/tmp/claude-501/-Users-marco-Desktop-IW-eftcamb-update-EFTCAMB/"
    "2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad/audit/trees/"
)

try:  # hold a machine-wide CPU slot (audit harness) when available; the script also runs without it
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "checks"))
    import harness  # noqa: E402

    _slot = harness.slot
except ImportError:  # pragma: no cover
    import contextlib

    _slot = lambda n=1: contextlib.nullcontext()  # noqa: E731

WORKER = r"""
import sys, json
sys.path.insert(0, sys.argv[1])
import camb, numpy as np
omk = float(sys.argv[2])
pars = camb.CAMBparams()
pars.set_cosmology(H0=67.5, ombh2=0.0224, omch2=0.12, omk=omk, mnu=0.06, tau=0.06)
pars.InitPower.set_params(As=2.1e-9, ns=0.965)
pars.set_for_lmax(2000, lens_potential_accuracy=0)
pars.Reion.max_redshift = 50.0   # A0 changed the default to 1000; align with CAMB
if hasattr(pars, "EFTCAMB"):
    assert pars.EFTCAMB.EFTflag == 0
res = camb.get_results(pars)
cl = res.get_unlensed_scalar_cls(CMB_unit="muK", raw_cl=True)
print(json.dumps({"TT": cl[:31, 0].tolist(), "EE": cl[:31, 1].tolist(), "TE": cl[:31, 3].tolist()}))
"""


def run(tree, omk):
    env = dict(os.environ, OMP_NUM_THREADS="2")
    env.pop("PYTHONPATH", None)
    with _slot(2):
        out = subprocess.run([sys.executable, "-c", WORKER, TREES + tree, str(omk)], env=env, capture_output=True,
                             text=True, check=True)
    return json.loads(out.stdout.strip().splitlines()[-1])


for omk in (0.0, -0.05, 0.05):
    a, b = run("A0", omk), run("B0", omk)
    print(f"omk = {omk:+.2f}   (A0 = EFTCAMB with EFTflag=0, B0 = plain CAMB v2; expected A0/B0 - 1 ~ 1e-13)")
    for ell in (2, 3, 4, 5, 10, 20, 30):
        tt = a["TT"][ell] / b["TT"][ell] - 1
        ee = a["EE"][ell] / b["EE"][ell] - 1
        te = (a["TE"][ell] - b["TE"][ell]) / (b["TT"][ell] * b["EE"][ell]) ** 0.5
        print(f"  ell={ell:3d}  dTT/TT={tt:+.3e}  dEE/EE={ee:+.3e}  dTE/sqrt(TTEE)={te:+.3e}")
