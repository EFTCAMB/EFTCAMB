"""Repro (t9-wrapper): camb.set_params(cp=pars, ...) silently resets the EFTCAMB model of `pars` to GR.

camb/camb.py set_params (EFTCAMB MOD) always calls cp.EFTCAMB.initialize_parameters(cp, params) with only the
keywords of this call. initialize_EFTCAMB_from_py (11_EFTCAMB_python.f90:115-146) first resets cp.EFTCAMB to the
defaults and then reads the new keywords, so every EFT setting not repeated in the call is lost (EFTflag -> 0).
  case A: set_params(cp=pars, H0=pars.H0)            (change only a CAMB parameter, here even to the same value)
  case B: set_params(cp=pars, EFTOmega0=0.02)        (change one EFT parameter of the model)
Upstream CAMB documents cp as "use this CAMBparams instead of creating a new one": parameters not passed are kept.
Expected: the EFT model is kept (case A: spectra identical to read_ini + set_cosmology(H0=...), which is what
set_params does to the CAMB part; case B: same model with the new EFTOmega0).

Run: /Users/marco/Software/miniforge3/bin/python set-params-cp-drops-eft-model.py [--tree PATH]
"""

import argparse
import os
import subprocess
import sys
import tempfile
from pathlib import Path

REPO = Path(__file__).resolve().parents[4]
INI = REPO / "fortran" / "eftcamb_test" / "parameters" / "2_PEFT_Omega_const_1.ini"
TREE = ("/private/tmp/claude-501/-Users-marco-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/"
        "scratchpad/audit/trees/A0")

CHILD = r"""
import sys, os
sys.path.insert(0, sys.argv[1])
import numpy as np
import camb
ini, tmp = sys.argv[2], sys.argv[3]
fast = os.path.join(tmp, 'model.ini')
open(fast, 'w').write(f'DEFAULT({ini})\nl_max_scalar = 1200\nget_tensor_cls = F\nget_transfer = F\n'
                      'EFTCAMB_write_background = F\n')
def tt(p):
    return camb.get_results(p).get_cmb_power_spectra(CMB_unit='muK', raw_cl=True)['total'][2:1200, 0]
def state(p):
    return f'EFTflag={p.EFTCAMB.EFTflag} model={p.EFTCAMB.model_name()!r}'
# control: what set_params(cp, H0=...) does to the CAMB part (set_cosmology resets the other cosmological
# parameters to their defaults, upstream behaviour), applied directly; the EFT model must stay as it is.
q = camb.read_ini(fast)
q.set_cosmology(H0=q.H0)
t_ctrl = tt(q)
p = camb.read_ini(fast)
print('read_ini              :', state(p), 'EFTOmega0 =', p.EFTCAMB.param_values())
camb.set_params(cp=p, H0=p.H0)
print('A set_params(cp, H0)  :', state(p), f'| max|TT/TT_control-1| = {np.max(np.abs(tt(p) / t_ctrl - 1)):.3e} '
      '(expected 0; control = read_ini + set_cosmology(H0) with the EFT model untouched)', flush=True)
p = camb.read_ini(fast)
try:
    camb.set_params(cp=p, EFTOmega0=0.02)
    exc = 'no exception'
except Exception as e:
    exc = f'{type(e).__name__}: {e}'
print('B set_params(cp, EFTOmega0=0.02):', exc, '| afterwards', state(p), 'values', (p.EFTCAMB.param_values() if p.EFTCAMB.EFTflag else 'n/a (EFTflag=0)'),
      '(expected no exception, EFTflag=1 Standard Pure EFT with EFTOmega0=0.02)', flush=True)
"""

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--tree", default=TREE)
    a = ap.parse_args()
    env = dict(os.environ, OMP_NUM_THREADS="2")
    env.pop("PYTHONPATH", None)
    with tempfile.TemporaryDirectory() as tmp:
        r = subprocess.run([sys.executable, "-c", CHILD, a.tree, str(INI), tmp], cwd=tmp, env=env,
                           capture_output=True, text=True)
        print("\n".join(ln for ln in r.stdout.splitlines() if ln[:1] in "rAB"))
        if r.returncode:
            print("child failed rc =", r.returncode, r.stderr[-800:])
