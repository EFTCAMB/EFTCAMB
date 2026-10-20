"""Repro (t9-wrapper): a malformed value passed to EFTCAMB.initialize_parameters (or camb.set_params with EFT keys)
terminates the Python interpreter (Fortran ERROR STOP in forutils IniObjects) instead of raising an exception.

initialize_parameters converts every value with str() (only python bool becomes T/F) and feeds it to the module
TIniFile of 11_EFTCAMB_python.f90; a value the Fortran reader cannot parse calls TNameValueList_Error -> FailStop ->
'error stop'. Values that do this and are natural in Python:
  numpy bool (np.True_ -> 'True'), float for an integer flag (1.0 -> '1.0'), None, a non-numeric string,
  bytes (b'1' -> "b'1'").
The ini-file path (camb.read_ini) is protected by validate_ini_file (a subprocess); the dictionary path is not.
Expected: a Python exception (ValueError) and a usable interpreter.

Run: /Users/marco/Software/miniforge3/bin/python initialize-parameters-bad-value-kills-python.py [--tree PATH]
"""

import argparse
import os
import subprocess
import sys

TREE = ("/private/tmp/claude-501/-Users-marco-Desktop-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/"
        "scratchpad/audit/trees/A0")

CHILD = r"""
import sys
sys.path.insert(0, sys.argv[1])
import numpy as np
import camb
key, expr, via = sys.argv[2], sys.argv[3], sys.argv[4]
d = {'EFTflag': 1, 'PureEFTmodel': 1, 'PureEFTmodelOmega': 1, 'EFTOmega0': 0.01, 'EFTwDE': 0}
d[key] = eval(expr)
try:
    if via == 'set_params':
        camb.set_params(H0=67.5, ombh2=0.022, omch2=0.12, **d)
    else:
        p = camb.CAMBparams()
        p.set_cosmology(H0=67.5)
        p.EFTCAMB.initialize_parameters(p, d, print_header=False)
    print('RESULT no exception', flush=True)
except Exception as e:
    print('RESULT python exception', type(e).__name__, flush=True)
print('RESULT interpreter alive', flush=True)
"""

CASES = [
    ("EFT_ghost_stability", "np.True_", "initialize_parameters"),
    ("EFTflag", "1.0", "initialize_parameters"),
    ("EFTOmega0", "None", "initialize_parameters"),
    ("EFTOmega0", "'0.0l'", "initialize_parameters"),
    ("EFTOmega0", "b'0.01'", "initialize_parameters"),
    ("EFTwDE", "np.float64(0)", "set_params"),
]

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--tree", default=TREE)
    a = ap.parse_args()
    env = dict(os.environ, OMP_NUM_THREADS="1")
    env.pop("PYTHONPATH", None)
    for key, expr, via in CASES:
        r = subprocess.run([sys.executable, "-c", CHILD, a.tree, key, expr, via], env=env, capture_output=True,
                           text=True)
        res = [ln[7:] for ln in r.stdout.splitlines() if ln.startswith("RESULT")]
        err = [ln.strip() for ln in (r.stdout + r.stderr).splitlines() if "Error for key" in ln or "STOP" in ln]
        print(f"{via:22s} {key}={expr:14s}: rc={r.returncode} {res} {err}  (expected: python exception, alive)")
