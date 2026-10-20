"""Repro (t9-wrapper): EFTCAMB.read_parameters() -> initialize_parameters() round trip is broken for string values
and for numbers written with a Fortran D exponent.

camb/eftcamb.py read_parameters: autoconvert tries bool/int/float on the raw *bytes* value and otherwise returns the
bytes object. initialize_parameters then does str(value), which for bytes gives "b'...'":
  * output_root (present for every EFT model) comes back as b'...' and is re-fed as the string "b'...'";
  * '1.d-4' (EFTB0 in 4_Des_F_R_LCDM_2.ini; valid Fortran, rejected by python float()) comes back as b'1.d-4'
    and re-fed as "b'1.d-4'": the Fortran ini reader hits 'error reading double' and ERROR STOP kills the
    interpreter (no python exception).
Expected: read_parameters returns str for strings and float for D-exponent numbers, and feeding its output back to
initialize_parameters reproduces the model.

Run: /Users/marco/Software/miniforge3/bin/python read-parameters-roundtrip-bytes.py [--tree PATH]
"""

import argparse
import os
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[4]
PARAMS = REPO / "fortran" / "eftcamb_test" / "parameters"
TREE = ("/private/tmp/claude-501/-Users-marco-Desktop-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/"
        "scratchpad/audit/trees/A0")

CHILD = r"""
import sys
sys.path.insert(0, sys.argv[1])
import camb
ini = sys.argv[2]
p = camb.read_ini(ini)
d = p.EFTCAMB.read_parameters()
for k in ('output_root', 'EFTB0'):
    if k in d:
        print(f'read_parameters()[{k!r}] = {d[k]!r}  (type {type(d[k]).__name__})')
q = camb.CAMBparams()
q.set_cosmology(H0=67)
try:
    q.EFTCAMB.initialize_parameters(q, d, print_header=False)
    print('initialize_parameters(read_parameters()) OK; output_root now', repr(q.EFTCAMB.read_parameters()['output_root']))
except Exception as e:
    print('python exception:', type(e).__name__, e)
print('REACHED_END')
"""

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--tree", default=TREE)
    a = ap.parse_args()
    env = dict(os.environ, OMP_NUM_THREADS="1")
    env.pop("PYTHONPATH", None)
    for model in ("2_PEFT_Omega_power_1.ini", "4_Des_F_R_LCDM_2.ini"):
        print(f"--- {model}")
        r = subprocess.run([sys.executable, "-c", CHILD, a.tree, str(PARAMS / model)], env=env, capture_output=True,
                           text=True, cwd=str(REPO / "fortran"))
        out = [ln for ln in (r.stdout + r.stderr).splitlines() if ln.strip() and not ln.startswith((" *", "*"))]
        print("\n".join(out[-6:]))
        print(f"child return code = {r.returncode}; reached end = {'REACHED_END' in r.stdout}  "
              "(expected: 0 / True, output_root unchanged)")
