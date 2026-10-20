"""Repro (t9-wrapper): CAMBparams.write_ini / camb.write_ini do not write any EFTCAMB key.

camb/_ini.py:_update_ini_state_from_params writes the CAMB parameters and the DarkEnergy/Reion/InitPower/Recomb/
NonLinear classes, but nothing for CAMBparams.EFTCAMB. For any EFTCAMB model:
  * write_ini(validate=True) (the default) writes the file and then raises CAMBValueError
    "Saved ini did not round-trip ... params.EFTCAMB.EFTflag differs: 0 != 1";
  * write_ini(validate=False) silently writes a GR (EFTflag = 0) ini: read_ini of it gives LCDM spectra.
Expected: the written ini reproduces the model (same EFT flags/parameters, bit-identical spectra).

Run: /Users/marco/Software/miniforge3/bin/python write-ini-drops-eftcamb.py [--tree PATH]
"""

import argparse
import os
import subprocess
import sys
import tempfile
from pathlib import Path

REPO = Path(__file__).resolve().parents[4]
TREE = ("/private/tmp/claude-501/-Users-marco-Desktop-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/"
        "scratchpad/audit/trees/A0")
INI = REPO / "fortran" / "eftcamb_test" / "parameters" / "2_PEFT_Omega_power_1.ini"

CHILD = r"""
import sys, os
sys.path.insert(0, sys.argv[1])
import numpy as np
import camb
ini, tmp = sys.argv[2], sys.argv[3]
fast = os.path.join(tmp, 'model.ini')
open(fast, 'w').write(f'DEFAULT({ini})\nl_max_scalar = 1200\nget_tensor_cls = F\nget_transfer = F\n'
                      'EFTCAMB_write_background = F\n')
p = camb.read_ini(fast)
print('original: EFTflag =', p.EFTCAMB.EFTflag, '| model =', p.EFTCAMB.model_name())
out = os.path.join(tmp, 'written.ini')
try:
    p.write_ini(out)
    print('write_ini(validate=True): no exception')
except Exception as e:
    print('write_ini(validate=True):', type(e).__name__, str(e)[:160])
p.write_ini(out, validate=False)
txt = open(out).read()
eft_lines = [l for l in txt.splitlines() if l.split('=')[0].strip() in p.EFTCAMB.read_parameters()]
print('EFTCAMB keys in written ini:', len(eft_lines), eft_lines[:5])
q = camb.read_ini(out)
print('re-read:  EFTflag =', q.EFTCAMB.EFTflag, '| model =', repr(q.EFTCAMB.model_name()))
c1 = camb.get_results(p).get_cmb_power_spectra(CMB_unit='muK', raw_cl=True)['total'][2:1000, 0]
c2 = camb.get_results(q).get_cmb_power_spectra(CMB_unit='muK', raw_cl=True)['total'][2:1000, 0]
print(f'max |TT_reread/TT_original - 1| (l<1000) = {np.max(np.abs(c2 / c1 - 1)):.3e}   (expected 0)')
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
        print("\n".join(line for line in r.stdout.splitlines() if not line.startswith((" ", "*"))
                        or "EFTCAMB keys" in line))
        if r.returncode:
            print("child failed rc =", r.returncode, r.stderr[-800:])
