"""Repro: EFT_IC_type = 3 returns all-NaN spectra, without any error, when Omega(a_eval) = Gamma3(a_eval) = 0.

equations.f90 (subroutine initial, case 3) and 09_EFTCAMB_IC.f90 (case 3) divide by
    3 Gamma30 + 5 Omega0 + 3 Gamma30 Omega0 + 3 Omega0^2      and      8 (3 Gamma30 (1 + Omega0) + Omega0 (5 + 3 Omega0))
which vanish for any model with Omega = Gamma3 = 0 at the turn-on (e.g. the pure-EFT Gamma1 test model).  The
adiabatic initial conditions of every k mode become 0/0 = NaN and camb.get_results returns NaN C_l, P(k), sigma8.

Run: /Users/marco/Software/miniforge3/bin/python ic3-nan-when-omega-gamma3-zero.py   (one CAMB run, ~1-5 min)
"""

import json
import os
import subprocess
import sys
import tempfile

ROOT = "/private/tmp/claude-501/-Users-marco-Desktop-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad/audit"
TREE = os.environ.get("REPRO_TREE", ROOT + "/trees/A0")
PARAMS = "/Users/marco/Desktop/IW/eftcamb_update/EFTCAMB/fortran/eftcamb_test/parameters"
CHILD = r"""
import sys, json, numpy as np
sys.path.insert(0, sys.argv[1])
import camb
cp = camb.read_ini(sys.argv[2], no_validate=True)
try:
    res = camb.get_results(cp)
except Exception as e:
    print('JSON' + json.dumps({'error': str(e)[:200]})); sys.exit(0)
cl = res.get_cmb_power_spectra(CMB_unit='muK', raw_cl=True)['total'][2:2501]
print('JSON' + json.dumps({'nan_frac': float(np.isnan(cl).mean()), 'sigma8': float(res.get_sigma8()[0]),
                           'TT_l100': float(cl[98, 0])}))
"""


def run(ic):
    tmp = tempfile.mkdtemp()
    ini = os.path.join(tmp, "w.ini")
    with open(ini, "w") as f:
        f.write(f"DEFAULT({PARAMS}/2_PEFT_gamma1_const_2.ini)\noutput_root = {tmp}/r\nfeedback_level = 0\n"
                f"EFTCAMB_write_background = F\nEFT_IC_type = {ic}\n")
    env = dict(os.environ, OMP_NUM_THREADS="2")
    env.pop("PYTHONPATH", None)
    p = subprocess.run([sys.executable, "-c", CHILD, TREE, ini], capture_output=True, text=True, cwd=tmp, env=env)
    for line in p.stdout.splitlines():
        if line.startswith("JSON"):
            return json.loads(line[4:])
    return {"error": f"rc={p.returncode} " + p.stderr[-300:]}


for ic in (1, 3):
    print(f"2_PEFT_gamma1_const_2 (Omega = 0, Gamma3 = 0) EFT_IC_type = {ic}: {run(ic)}")
print("expected: finite spectra for both (or an explicit error for IC 3); observed: IC 3 gives NaN with no error")
