"""Repro: pickle / copy.deepcopy / repr+eval of a CAMBparams with EFTCAMB active drops the EFTCAMB model.

The python EFTCAMB class maps only the flag components of TEFTCAMB (the Fortran model, outroot and stored read
parameters are not mapped), so CAMB_Structure.__getstate__ saves the flags only. The restored parameters have
EFTflag /= 0 but TEFTCAMB%model unallocated, and get_background / get_results dereference it: SIGSEGV.
CAMBparams documents pickle and repr/eval as supported; multiprocessing pickles parameters.
Setting EFTflag from python without initialize_parameters hits the same unguarded dereference.

Usage: python eft-pickle-repr-loses-model.py [--tree PATH]
"""

import argparse
import os
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = "/private/tmp/claude-501/-Users-marco-Desktop-IW-eftcamb-update-EFTCAMB/2b9dddd0-ab26-46fb-8810-571fb7648a45/scratchpad/audit"
PARAMS = Path(__file__).resolve().parents[2] / "parameters"

CHILD = r'''
import sys, pickle, copy
sys.path.insert(0, sys.argv[1])
import camb
mode, ini = sys.argv[2], sys.argv[3]
if mode == "python_flag_only":
    p2 = camb.set_params(H0=70, ombh2=0.0226, omch2=0.112, mnu=0.06)
    p2.EFTCAMB.EFTflag = 3
else:
    cp = camb.read_ini(ini, no_validate=True)
    if mode == "reference":
        p2 = cp
    elif mode == "pickle":
        p2 = pickle.loads(pickle.dumps(cp))
    elif mode == "deepcopy":
        p2 = copy.deepcopy(cp)
    elif mode == "repr_eval":
        p2 = eval(repr(cp), {"camb": camb})
    elif mode == "copy":
        p2 = cp.copy()
print("STATE EFTflag=%d model_name=%r" % (p2.EFTCAMB.EFTflag, p2.EFTCAMB.model_name()), flush=True)
try:
    print("RESULT", camb.get_background(p2).hubble_parameter([0.5, 1.0]).tolist(), flush=True)
except Exception as e:
    print("RESULT exception", type(e).__name__, str(e)[:120], flush=True)
'''


def run(tree, mode, ini):
    env = dict(os.environ, OMP_NUM_THREADS="2")
    env.pop("PYTHONPATH", None)
    p = subprocess.run([sys.executable, "-c", CHILD, tree, mode, ini], capture_output=True, text=True, env=env,
                       cwd=os.path.dirname(ini))
    lines = [ln for ln in p.stdout.splitlines() if ln.startswith(("STATE", "RESULT"))]
    return p.returncode, lines


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tree", default=f"{ROOT}/trees/A0")
    a = ap.parse_args()
    tmp = tempfile.mkdtemp(prefix="pickle_")
    ini = os.path.join(tmp, "fR.ini")
    Path(ini).write_text(f"DEFAULT({PARAMS / '4_Des_F_R_LCDM_1.ini'})\noutput_root = {tmp}/fR\n"
                         "EFTCAMB_write_background = F\nfeedback_level = 0\n")
    print(f"tree: {a.tree}")
    for mode in ("reference", "copy", "pickle", "deepcopy", "repr_eval", "python_flag_only"):
        rc, lines = run(a.tree, mode, ini)
        status = "OK" if rc == 0 else (f"CRASH signal {-rc}" if rc < 0 else f"CRASH rc={rc}")
        print(f"{mode:17s} observed: {status:16s} {' | '.join(lines)}")
    print("expected: every round trip (copy/pickle/deepcopy/repr_eval) rc=0 with the reference H(z); "
          "python_flag_only: a CAMBError instead of a crash")


if __name__ == "__main__":
    main()
