"""Repro: EFTCAMB library code calls MpiStop / stop / error stop, which terminates the Python interpreter instead of
raising camb.CAMBError (global_error_flag / error handler).

Cases (tree A0, valid CAMB inputs):
  1. 5_Kmouflage_1.ini with massive neutrinos (omnuh2 = 0.0006, the CAMB default 1 massive eigenstate):
     008p3_Kmouflage.f90:438  MpiStop('EFTCAMB error') in EFTCAMBKmouflageInitBackground.
  2. 2_PEFT_Omega_const_1.ini with EFT_IC_type = 4 (outside 1..3):
     equations.f90:2361  MpiStop inside initial(), i.e. inside the OpenMP k loop.
  3. 2_PEFT_Omega_const_1.ini with EFTCAMB_use_background = T:
     09_EFTCAMB_main.f90:250  MpiStop while reading the flags (python read_ini).
For each case the child process runs camb.read_ini + camb.get_results inside try/except camb.CAMBError and prints
CAUGHT if an exception reaches Python.  Expected: CAUGHT and the interpreter survives.
Usage: /Users/marco/Software/miniforge3/bin/python stop-kills-python.py [TREE]
"""

from __future__ import annotations

import os
import subprocess
import sys
import tempfile
from pathlib import Path

REPO = Path(__file__).resolve().parents[4]
TREE = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(
    "/private/tmp/claude-501/-Users-marco-Desktop-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/"
    "scratchpad/audit/trees/A0"
)
PAR = REPO / "fortran/eftcamb_test/parameters"
CASES = [
    ("Kmouflage + massive neutrinos", "5_Kmouflage_1.ini", {"omnuh2": "0.0006", "massive_neutrinos": "1"}),
    ("EFT_IC_type = 4", "2_PEFT_Omega_const_1.ini", {"EFT_IC_type": "4"}),
    ("EFTCAMB_use_background = T (non-Horndeski)", "2_PEFT_Omega_const_1.ini", {"EFTCAMB_use_background": "T"}),
]
CHILD = r"""
import sys
sys.path.insert(0, sys.argv[1])
import camb
try:
    pars = camb.read_ini(sys.argv[2], no_validate=True)
    pars.WantTensors = False
    res = camb.get_results(pars)
    print("NO-ERROR")
except camb.CAMBError as e:
    print("CAUGHT", repr(e)[:200])
print("INTERPRETER-ALIVE")
"""


def main():
    tmp = Path(tempfile.mkdtemp(prefix="stopkill_"))
    print(f"tree: {TREE}")
    for label, ini, extra in CASES:
        w = tmp / f"w_{len(label)}.ini"
        lines = [f"DEFAULT({PAR / ini})", "feedback_level = 1", "EFTCAMB_write_background = F", "get_tensor_cls = F",
                 "do_lensing = F", f"output_root = {tmp / 'out'}"] + [f"{k} = {v}" for k, v in extra.items()]
        w.write_text("\n".join(lines) + "\n")
        r = subprocess.run([sys.executable, "-c", CHILD, str(TREE), str(w)], capture_output=True, text=True, cwd=tmp,
                           env=dict(os.environ, OMP_NUM_THREADS="1"), timeout=900)
        alive = "INTERPRETER-ALIVE" in r.stdout
        caught = "CAUGHT" in r.stdout
        tail = [ln.strip() for ln in (r.stdout + r.stderr).splitlines() if ln.strip()][-3:]
        print(f"- {label}: exit code {r.returncode}, python exception caught: {caught}, interpreter alive after "
              f"the call: {alive}\n    last output: {tail}")
    print("expected for every case: python exception caught: True, interpreter alive: True")


if __name__ == "__main__":
    main()
